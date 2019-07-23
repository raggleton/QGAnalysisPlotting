#!/usr/bin/env python


"""Folding testing
"""


from __future__ import print_function

import os
import sys
import argparse
from array import array
import numpy as np
import math

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()

import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat(".3f")


# Control plot output format
OUTPUT_FMT = "pdf"



def th2_to_tmatrixsparse(hist_A, oflow_x=False, oflow_y=False):
    """Convert TH2 to TMatrixSparse

    Mainly taken from TUnfold ctor: https://root.cern.ch/doc/master/TUnfold_8cxx_source.html#l01715

    Assumes histmap == kHistMapOutputHoriz
    """
    ncol = hist_A.GetNbinsX()
    if oflow_x:
        ncol += 2
    nrow = hist_A.GetNbinsY()
    if oflow_y:
        nrow += 2

    r = ROOT.TMatrixD(nrow, ncol)  #nrow, ncol

    for iy in range(0, nrow):
        for ix in range(0, ncol):
            z = hist_A.GetBinContent(ix, iy+1)

    return matrix


def th1_to_ndarray(hist_A, oflow_x=False):
    """Convert TH1 to numpy ndarray"""
    ncol = hist_A.GetNbinsX()
    if oflow_x:
        ncol += 2
    result = np.zeros(shape=(1, ncol), dtype=float)
    errors = np.zeros(shape=(1, ncol), dtype=float)
    # access via result[irow][icol]

    # Get ROOT indices to loop over
    x_start = 0 if oflow_x else 1
    x_end = hist_A.GetNbinsX()
    if oflow_x:
        x_end += 1

    # x_ind for numpy as always starts at 0
    # ix for ROOT
    for x_ind, ix in enumerate(range(x_start, x_end+1)):
        result[0][x_ind] = hist_A.GetBinContent(ix)
        errors[0][x_ind] = hist_A.GetBinError(ix)

    # check sparsity
    return result, errors


def th2_to_ndarray(hist_A, oflow_x=False, oflow_y=False):
    """Convert TH2 to numpy ndarray

    Don't use verison in common_utils - wrong axes?
    """
    ncol = hist_A.GetNbinsX()
    if oflow_x:
        ncol += 2
    nrow = hist_A.GetNbinsY()
    if oflow_y:
        nrow += 2

    result = np.zeros(shape=(nrow, ncol), dtype=float)
    errors = np.zeros(shape=(nrow, ncol), dtype=float)
    # access via result[irow][icol]

    # Get ROOT indices to loop over
    y_start = 0 if oflow_y else 1
    y_end = hist_A.GetNbinsY()
    if oflow_y:
        y_end += 1

    x_start = 0 if oflow_x else 1
    x_end = hist_A.GetNbinsX()
    if oflow_x:
        x_end += 1

    # y_ind, x_ind for numpy as always starts at 0
    # iy, ix for ROOT
    for y_ind, iy in enumerate(range(y_start, y_end+1)):
        for x_ind, ix in enumerate(range(x_start, x_end+1)):
            result[y_ind][x_ind] = hist_A.GetBinContent(ix, iy)
            errors[y_ind][x_ind] = hist_A.GetBinError(ix, iy)

    # check sparsity
    num_empty = np.count_nonzero(result == 0)
    num_entries = result.size
    sparsity = num_empty / float(num_entries)
    print("Converting TH2 to ndarray...")
    print("num_empty:", num_empty)
    print("num_entries:", num_entries)
    print("sparsity:", sparsity)
    if (sparsity > 0.5):
        print("Matrix has %d/%d empty entries - consider using sparse matrix (which I don't know how to do yet)" % (num_empty, num_entries))

    return result, errors


def ndarray_to_th1(nd_array, has_oflow_x=False):
    """Convert numpy ndarray row vector to TH1, with shape (1, nbins)

    Use has_oflow_x to include the under/overflow bins
    """
    nbinsx = nd_array.shape[1]
    nbins_hist = nbinsx
    if has_oflow_x:
        nbins_hist -= 2

    # need the 0.5 offset to match TUnfold
    h = ROOT.TH1F(cu.get_unique_str(), "", nbins_hist, 0.5, nbins_hist+0.5)

    x_start = 1
    x_end = nbins_hist

    if has_oflow_x:
        x_start = 0
        x_end = nbins_hist+1

    for x_ind, ix in enumerate(range(x_start, x_end+1)):
        h.SetBinContent(ix, nd_array[0][x_ind])
        #FIXME how to do errors
    return h


def ndarray_to_th2(nd_array, has_oflow_x=False, has_oflow_y=False):
    pass


def normalise_ndarray_by_row(matrix):
    matrix = matrix.T # makes life a bit easier
    for i in range(matrix.shape[0]):
        row_sum = matrix[i].sum()
        if row_sum != 0:
            matrix[i] = matrix[i] / row_sum
    return matrix.T


def normalise_ndarray_by_col(matrix):
    for i in range(matrix.shape[0]):
        row_sum = matrix[i].sum()
        if row_sum != 0:
            matrix[i] = matrix[i] / row_sum
    return matrix


def get_folded_hist(hist_mc_gen_reco_map, hist_mc_gen):
    oflow = True
    # Convert map to matrix
    response_matrix, response_matrix_err = th2_to_ndarray(hist_mc_gen_reco_map, oflow_x=oflow, oflow_y=oflow)

    # Normalise response_matrix so that bins represent prob to go from
    # given gen bin to a reco bin
    # ASSUMES GEN ON X AXIS!
    response_matrix = normalise_ndarray_by_row(response_matrix)

    # Convert hist to vector
    gen_vec, gen_vec_err = th1_to_ndarray(hist_mc_gen, oflow_x=oflow)

    # Multiply
    # Note that we need to transpose from row vecc to column vec
    folded_vec = response_matrix.dot(gen_vec.T)
    print("response_matrix.shape:", response_matrix.shape)
    print("gen_vec.shape", gen_vec.shape)
    print("folded_vec.shape", folded_vec.shape)

    # Convert vector to TH1
    folded_hist = ndarray_to_th1(folded_vec.T, has_oflow_x=oflow)

    return folded_hist


def generate_2d_canvas(size=(800, 600)):
    canv = ROOT.TCanvas(cu.get_unique_str(), "", *size)
    canv.SetTicks(1, 1)
    canv.SetRightMargin(1.5)
    canv.SetLeftMargin(0.9)
    return canv


def draw_response_matrix(rsp_map, region_name, variable_name, output_filename, draw_values=False):
    canv = generate_2d_canvas()
    canv.SetLogz()

    rsp_map.SetTitle("Response matrix, %s region, %s;Generator Bin;Detector Bin" % (region_name, variable_name))
    rsp_map.GetYaxis().SetTitleOffset(1.5)
    rsp_map.GetXaxis().SetTitleOffset(1.5)
    if draw_values:
        rsp_map.Draw("COLZ TEXT")
    else:
        rsp_map.Draw("COLZ")
    canv.SaveAs(output_filename)


def draw_folded_hists(hist_mc_folded, hist_mc_reco, hist_data_reco, output_filename, title=""):
    entries = []

    if hist_mc_folded:
        entries.append(
            Contribution(hist_mc_folded, label="Folded MC [detector-level]",
                         line_color=ROOT.kGreen+2, line_width=1,
                         marker_color=ROOT.kGreen+2, marker_size=0,
                         normalise_hist=False, subplot=hist_data_reco),
        )

    if hist_mc_reco:
        entries.append(
            Contribution(hist_mc_reco, label="Reco MC [detector-level]",
                         line_color=ROOT.kAzure+2, line_width=1,
                         marker_color=ROOT.kAzure+2, marker_size=0,
                         normalise_hist=False, subplot=hist_data_reco),
        )

    if hist_data_reco:
        entries.append(
            Contribution(hist_data_reco, label="Reco Data [detector-level]",
                         line_color=ROOT.kRed, line_width=0,
                         marker_color=ROOT.kRed, marker_size=0.6, marker_style=20,
                         normalise_hist=False),
        )

    plot = Plot(entries,
                what='hist',
                title=title,
                xtitle="Bin number",
                ytitle="N",
                subplot_type='ratio',
                subplot_title='MC/Data',
                subplot_limits=(0.25, 1.75))
    plot.default_canvas_size = (800, 600)
    plot.plot("NOSTACK HISTE")
    plot.main_pad.SetLogy(1)
    ymax = max(h.GetMaximum() for h in [hist_mc_folded, hist_mc_reco, hist_data_reco] if h)
    plot.container.SetMaximum(ymax * 100)
    plot.container.SetMinimum(1)
    plot.legend.SetY1NDC(0.77)
    plot.legend.SetX2NDC(0.85)
    plot.save(output_filename)


def draw_folded_hists_physical(hist_mc_folded, hist_mc_reco, hist_data_reco, output_filename, title="", xtitle="", logx=False, logy=False):
    entries = []

    if hist_mc_folded:
        entries.append(
            Contribution(hist_mc_folded, label="Folded MC [detector-level]",
                         line_color=ROOT.kGreen+2, line_width=1,
                         marker_color=ROOT.kGreen+2, marker_size=0,
                         normalise_hist=False, subplot=hist_data_reco),
        )

    if hist_mc_reco:
        entries.append(
            Contribution(hist_mc_reco, label="Reco MC [detector-level]",
                         line_color=ROOT.kAzure+2, line_width=1,
                         marker_color=ROOT.kAzure+2, marker_size=0,
                         normalise_hist=False, subplot=hist_data_reco),
        )

    if hist_data_reco:
        entries.append(
            Contribution(hist_data_reco, label="Reco Data [detector-level]",
                         line_color=ROOT.kRed, line_width=0,
                         marker_color=ROOT.kRed, marker_size=0.6, marker_style=20,
                         normalise_hist=False),
        )

    plot = Plot(entries,
                what='hist',
                title=title,
                xtitle=xtitle,
                ytitle="N",
                subplot_type='ratio',
                subplot_title='MC/Data',
                subplot_limits=(0.25, 1.75))
    plot.default_canvas_size = (800, 600)
    plot.plot("NOSTACK HISTE")
    if logx:
        plot.set_logx(True)
    if logy:
        plot.main_pad.SetLogy(1)
    ymax = max(h.GetMaximum() for h in [hist_mc_folded, hist_mc_reco, hist_data_reco] if h)
    plot.container.SetMaximum(ymax * 100)
    plot.container.SetMinimum(1)
    plot.legend.SetY1NDC(0.77)
    plot.legend.SetX2NDC(0.85)
    plot.save(output_filename)


def draw_projection_comparison(h_orig, h_projection, title, xtitle, output_filename, print_bin_comparison=True):
    """Draw 2 hists, h_orig the original, and h_projection the projection of a 2D hist"""

    # Check integrals
    int_orig = h_orig.Integral()
    int_proj = h_projection.Integral()
    if abs(int_orig - int_proj)/int_orig > 0.01:
        print("draw_projection_comparison: different integrals: %f vs %f" % (int_orig, int_proj))

    # Check bin-by-bin
    if print_bin_comparison:
        for i in range(1, h_orig.GetNbinsX()+1):
            value_orig = h_orig.GetBinContent(i)
            value_proj = h_projection.GetBinContent(i)
            if abs(value_orig - value_proj) > 1E-2:
                print("draw_projection_comparison: bin %s has different contents: %f vs %f" % (i, value_orig, value_proj))

    entries = [
        Contribution(h_orig, label="1D hist",
                     line_color=ROOT.kBlue, line_width=1,
                     marker_color=ROOT.kBlue, marker_size=0,
                     normalise_hist=False),
        Contribution(h_projection, label="Response map projection",
                     line_color=ROOT.kRed, line_width=1,
                     marker_color=ROOT.kRed, marker_size=0,
                     normalise_hist=False,
                     subplot=h_orig),
    ]
    plot = Plot(entries,
                what='hist',
                title=title,
                xtitle=xtitle,
                ytitle="N",
                subplot_type='ratio',
                subplot_title='Projection / 1D',
                subplot_limits=(0.999, 1.001))
    plot.default_canvas_size = (800, 600)
    plot.plot("NOSTACK HIST")
    plot.main_pad.SetLogy(1)
    ymax = max(h.GetMaximum() for h in [h_orig, h_projection])
    plot.container.SetMaximum(ymax * 10)
    # plot.container.SetMinimum(1E-8)
    plot.legend.SetY1NDC(0.77)
    plot.legend.SetX2NDC(0.85)
    plot.save(output_filename)


def convert_th1_physical_bins(hist, bin_values):
    if hist.GetNbinsX() != len(bin_values)-1:
        print("hist has", hist.GetNbinsX(), "bins")
        print("bin_values =", bin_values, ", has", len(bin_values), "entries")
        raise RuntimeError("Incorrect number of bin_values for convert_th1_physical_bins()")

    h_new = ROOT.TH1D(hist.GetName()+"Physical"+cu.get_unique_str(),
                      ';'.join([hist.GetTitle(), hist.GetXaxis().GetTitle(), hist.GetYaxis().GetTitle()]),
                      hist.GetNbinsX(), array('f', bin_values))
    for i in range(0, h_new.GetNbinsX()+2):
        h_new.SetBinContent(i, hist.GetBinContent(i))
        h_new.SetBinError(i, hist.GetBinError(i))
    return h_new


if __name__ == "__main__":
    # FOR Z+JETS:
    input_mc_dy_mgpythia_tfile = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_noZReweight/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root")

    input_singlemu_tfile = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu.root")

    # FOR DIJET:
    input_mc_qcd_mgpythia_tfile = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_noZReweight/uhh2.AnalysisModuleRunner.MC.MC_QCD.root")

    input_mc_qcd_pythia_tfile = cu.open_root_file("workdir_ak4puppi_pythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_noZReweight/uhh2.AnalysisModuleRunner.MC.MC_PYTHIA-QCD.root")
    # input_mc_qcd_herwig_tfile = cu.open_root_file("workdir_ak4puppi_herwig_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1RecoConstituents_V11JEC_JER_tUnfold/uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD.root")

    input_jetht_tfile = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5/uhh2.AnalysisModuleRunner.DATA.Data_JetHTZeroBias.root")


    regions = [
        {
            "name": "Dijet",
            "dirname": "Dijet_QG_Unfold_tighter",
            "label": "Dijet",
            "data_tfile": input_jetht_tfile,
            # "mc_tfile": input_mc_qcd_mgpythia_tfile,
            "mc_tfile": input_mc_qcd_mgpythia_tfile,
            "mc_tfile_other": input_mc_qcd_pythia_tfile,
            "tau_limits": {
                'jet_puppiMultiplicity': (1E-13, 1E-10),
                'jet_pTD': (1E-13, 1E-10),
                'jet_LHA': (1E-13, 1E-10),
                'jet_width': (1E-13, 1E-10),
                'jet_thrust': (1E-13, 1E-10),
                'jet_puppiMultiplicity_charged': (1E-13, 1E-10),
                'jet_pTD_charged': (1E-13, 1E-10),
                'jet_LHA_charged': (1E-13, 1E-10),
                'jet_width_charged': (1E-13, 1E-10),
                'jet_thrust_charged': (1E-13, 1E-10),
            },
        },
        {
            "name": "ZPlusJets",
            "dirname": "ZPlusJets_QG_Unfold",
            "label": "Z+jets",
            "data_tfile": input_singlemu_tfile,
            "mc_tfile": input_mc_dy_mgpythia_tfile,
            "tau_limits": {
                'jet_puppiMultiplicity': (1E-10, 1E-4),
                'jet_pTD': (1E-10, 1E-4),
                'jet_LHA': (1E-10, 1E-4),
                'jet_width': (1E-10, 1E-4),
                'jet_thrust': (1E-10, 1E-4),
                'jet_puppiMultiplicity_charged': (1E-10, 1E-4),
                'jet_pTD_charged': (1E-10, 1E-4),
                'jet_LHA_charged': (1E-10, 1E-4),
                'jet_width_charged': (1E-10, 1E-4),
                'jet_thrust_charged': (1E-10, 1E-4),
            },
        },
    ]

    # If True, use part of MC for response matrix, and separate part for unfolding
    # as independent test
    MC_split = True
    mc_append = "_split" if MC_split else "_all"

    # Subtract fakes from reconstructed hists, using MC fakes as template
    subtract_fakes = True
    sub_append = "_subFakes" if subtract_fakes else ""

    output_dir = "folding_better_target0p5%s%s" % (mc_append, sub_append)
    cu.check_dir_exists_create(output_dir)

    # TODO automate this
    jet_algo = "AK4 PF PUPPI"

    # Save hists etc to ROOT file for access later
    output_tfile = ROOT.TFile("%s/folding_result.root" % (output_dir), "RECREATE")
    print("Saving hists to ROOT file:", output_tfile.GetName())

    for region in regions[:]:

        # Setup pt bins
        # -------------
        # need different ones for Z+Jets region
        is_zpj = "ZPlusJets" in region['name']

        zpj_append = "_zpj" if is_zpj else ""

        pt_bin_edges_gen = qgc.PT_UNFOLD_DICT['signal%s_gen' % (zpj_append)]
        pt_bin_edges_reco = qgc.PT_UNFOLD_DICT['signal%s_reco' % (zpj_append)]

        pt_bin_edges_underflow_gen = qgc.PT_UNFOLD_DICT['underflow%s_gen' % (zpj_append)]
        pt_bin_edges_underflow_reco = qgc.PT_UNFOLD_DICT['underflow%s_reco' % (zpj_append)]

        all_pt_bin_edges_gen = list(pt_bin_edges_underflow_gen[:-1]) + list(pt_bin_edges_gen)
        all_pt_bin_edges_reco = list(pt_bin_edges_underflow_reco[:-1]) + list(pt_bin_edges_reco)

        # Do 1D pt distribution
        # ---------------------
        angle_shortname = "pt"
        # Get response matrix, gen hist, reco hist (MC & Data)
        hist_data_reco = cu.get_from_tfile(region['data_tfile'], "%s/hist_%s_reco_all" % (region['dirname'], angle_shortname))
        mc_hname_append = "split" if MC_split else "all"
        hist_mc_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_%s" % (region['dirname'], angle_shortname, mc_hname_append))
        hist_mc_gen = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_truth_%s" % (region['dirname'], angle_shortname, mc_hname_append))
        hist_mc_gen_reco_map = cu.get_from_tfile(region['mc_tfile'], "%s/tu_%s_GenReco_%s" % (region['dirname'], angle_shortname, mc_hname_append))

        # To account for fraction that went into MC if split
        split_scale_factor = 1/.2 if MC_split else 1.
        hist_mc_reco.Scale(split_scale_factor)
        hist_mc_gen.Scale(split_scale_factor)

        if subtract_fakes:
            # to construct our "fakes" template, we use the ratio as predicted by MC, and apply it to data
            # this way we ensure we don't have -ve values, and avoid any issue with cross sections
            hist_mc_fakes_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_fake_%s" % (region['dirname'], angle_shortname, mc_hname_append))
            hist_mc_fakes_reco.Scale(split_scale_factor)

            # Create fraction template
            hist_fakes_reco_frac = hist_mc_fakes_reco.Clone("hist_%s_fakes" % angle_shortname)
            hist_fakes_reco_frac.Divide(hist_mc_reco)

            hist_data_reco_bg_subtracted = hist_data_reco.Clone(hist_data_reco.GetName() + "_bgrSubtracted")
            hist_data_reco_fakes = hist_data_reco.Clone()
            hist_data_reco_fakes.Multiply(hist_fakes_reco_frac)
            hist_data_reco_bg_subtracted.Add(hist_data_reco_fakes, -1)

            hist_mc_reco_bg_subtracted = hist_mc_reco.Clone(hist_mc_reco.GetName() + "_bgrSubtracted")
            hist_mc_reco_bg_subtracted.Add(hist_mc_fakes_reco, -1)

        append = "%s_%s" % (region['name'], angle_shortname)  # common str to put on filenames, etc
        angle_str = "p_{T} [GeV]"

        # put plots in subdir, to avoid overcrowding
        this_output_dir = "%s/%s/%s" % (output_dir, region['name'], angle_shortname)
        cu.check_dir_exists_create(this_output_dir)

        draw_response_matrix(hist_mc_gen_reco_map,
                             region['label'],
                             angle_str,
                             draw_values=False,
                             output_filename="%s/response_map_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

        draw_response_matrix(cu.make_normalised_TH2(hist_mc_gen_reco_map, 'X', recolour=False),
                             region['label'],
                             angle_str + " (normalised by gen bin)",
                             draw_values=True,
                             output_filename="%s/response_map_%s_normX.%s" % (this_output_dir, append, OUTPUT_FMT))

        # Do folding, plot
        # ----------------
        hist_mc_folded = get_folded_hist(hist_mc_gen_reco_map, hist_mc_gen)

        hist_mc_reco_1d = hist_mc_reco_bg_subtracted if subtract_fakes else hist_mc_reco
        hist_data_reco_1d = hist_data_reco_bg_subtracted if subtract_fakes else hist_data_reco

        # plot in raw bin numbers
        draw_folded_hists(hist_mc_folded,
                          hist_mc_reco_1d,
                          hist_data_reco_1d,
                          output_filename="%s/folded_hist_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                          title="%s region, %s jets, %s" % (region['label'], jet_algo, angle_str))

        # Convert to physical bin limits
        bin_values = all_pt_bin_edges_reco
        hist_mc_gen_physical = convert_th1_physical_bins(hist_mc_gen, all_pt_bin_edges_gen)
        hist_mc_folded_physical = convert_th1_physical_bins(hist_mc_folded, bin_values)
        hist_mc_reco_physical = convert_th1_physical_bins(hist_mc_reco_bg_subtracted, bin_values)
        hist_data_reco_physical = convert_th1_physical_bins(hist_data_reco_bg_subtracted, bin_values)

        draw_folded_hists_physical(hist_mc_folded_physical,
                                   hist_mc_reco_physical,
                                   hist_data_reco_physical,
                                   output_filename="%s/folded_hist_physical_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                                   title="%s region, %s jets" % (region['label'], jet_algo),
                                   xtitle=angle_str,
                                   logx=True,
                                   logy=True)

        # Save hists to ROOT file for later use
        new_tdir = "%s/%s" % (region['name'], "pt")
        output_tfile.mkdir(new_tdir)
        this_tdir = output_tfile.Get(new_tdir)
        this_tdir.cd()
        this_tdir.WriteTObject(hist_mc_gen_reco_map, "response_map")
        this_tdir.WriteTObject(hist_mc_gen, "gen_mc")
        this_tdir.WriteTObject(hist_mc_reco_1d, "reco_mc")
        this_tdir.WriteTObject(hist_data_reco_1d, "reco_data")
        this_tdir.WriteTObject(hist_mc_gen_physical, "gen_mc_physical")
        this_tdir.WriteTObject(hist_mc_reco_physical, "reco_mc_physical")
        this_tdir.WriteTObject(hist_data_reco_physical, "reco_data_physical")

    output_tfile.Close()

