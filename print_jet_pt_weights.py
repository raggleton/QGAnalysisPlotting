#!/usr/bin/env python

"""
Plots pt vs variable, and variable vs weights

To spot large-weight events, but also how the variables correlate with pT
"""

from __future__ import print_function, division

import os
os.nice(10)
import sys
from array import array
import numpy as np
import math
from itertools import product, chain
from copy import copy, deepcopy

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp

# Use rootpy to throw exceptions on ROOT errors, but need DANGER enabled
import rootpy
import rootpy.logger.magic as M; M.DANGER.enabled = True

ROOT.gErrorIgnoreLevel = ROOT.kWarning
# ROOT.gErrorIgnoreLevel = ROOT.kInfo
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)  # VERY IMPORTANT - somewhere, closing a TFile for exp systs deletes a map...dunno why

palette_2D = ROOT.kBird
palette_2D = ROOT.kRainBow
palette_1D = ROOT.kRainBow
ROOT.gStyle.SetPalette(palette_1D)

pt_str = "p_{T}^{RecoJet}"
pt_genjet_str = "p_{T}^{GenJet}"
genht_str = "H_{T}^{Gen}"
qscale_str = "q scale"
ptHat_str = "#hat{p}_{T}"
PU_ptHat_str = "max(#hat{p}_{T}^{PU})"
jetkT_str = "k_{T}^{parton}"

metric_str_dict = {
    "pt_jet_genHT_ratio" : "%s / %s" % (pt_str, genht_str),
    "pt_genjet_genHT_ratio" : "%s / %s" % (pt_genjet_str, genht_str),

    "pt_jet_qScale_ratio" : "%s / %s" % (pt_str, qscale_str),
    "pt_genjet_qScale_ratio" : "%s / %s" % (pt_genjet_str, qscale_str),

    "pt_jet_ptHat_ratio" : "%s / %s" % (pt_str, ptHat_str),
    "pt_genjet_ptHat_ratio" : "%s / %s" % (pt_genjet_str, ptHat_str),

    "PU_ptHat_genHT_ratio" : "%s / %s" % (PU_ptHat_str, genht_str),

    "PU_ptHat_ptHat_ratio" : "%s / %s" % (PU_ptHat_str, ptHat_str),
    
    "PU_ptHat_jetkT_ratio" : "%s / %s" % (PU_ptHat_str, jetkT_str),

    "pt_jet_jetkT_ratio" : "%s / %s" % (pt_str, jetkT_str),
}


def get_var_str(histname):
    for k, v in metric_str_dict.items():
        if k in histname:
            return v
    return histname.replace("weight_vs_pt_vs_", "").replace("weight_vs_pt_genjet_vs_", "")


def do_var_vs_pt_plot(histname, input_filename, output_filename):
    ROOT.gStyle.SetPalette(palette_2D)
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    if h3d.GetEntries() == 0:
        return
    h2d = h3d.Project3D("zy")

    xlabel = h2d.GetXaxis().GetTitle()
    ylabel = h2d.GetYaxis().GetTitle()
    ylabel = get_var_str(histname)

    # find largest var value (ie row) that has a filled bin
    h2d_ndarray = cu.th2_to_ndarray(h2d)[0]

    xbins = np.array(cu.get_bin_edges(h2d, 'x'))
    ybins = np.array(cu.get_bin_edges(h2d, 'y'))

    # remove dodgy bins with 0 width cos I was an idiot and duplicated some bins
    n_deleted = 0
    # weight bin
    # xax = h2d.GetXaxis()
    # for ix in range(1, h2d.GetNbinsX()+1):
    #     if xax.GetBinWidth(ix) == 0:
    #         h2d_ndarray = np.delete(h2d_ndarray, ix-1-n_deleted, axis=1)
    #         xbins = np.delete(xbins, ix-1-n_deleted, axis=0)
    #         n_deleted += 1
    #         print("Deleting bin", ix)

    # pt bin
    # n_deleted = 0
    # yax = h2d.GetYaxis()
    # for iy in range(1, h2d.GetNbinsY()+1):
    #     if yax.GetBinWidth(iy) == 0:
    #         h2d_ndarray = np.delete(h2d_ndarray, iy-1-n_deleted, axis=0)
    #         ybins = np.delete(ybins, iy-1-n_deleted, axis=0)
    #         n_deleted += 1
    #         print("Deleting bin", iy)

    # nonzero returns (row #s)(col #s) of non-zero elements
    # we only want the largest row #
    max_filled_row_ind = int(np.nonzero(h2d_ndarray)[0].max())
    h2d = cu.ndarray_to_th2(h2d_ndarray, binsx=xbins, binsy=ybins)

    if "unweighted" in histname:
        h2d.SetTitle("Unweighted;%s;%s" % (xlabel, ylabel))
    else:
        h2d.SetTitle("Weighted;%s;%s" % (xlabel, ylabel))

    h2d.GetYaxis().SetRange(1, max_filled_row_ind+2)  # +1 as ROOT 1-indexed, +1 for padding
    h2d.GetYaxis().SetTitle(get_var_str(histname))

    xmin = 15 if "pt_genjet_vs" in histname else 30
    xmax = 300

    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetLeftMargin(0.12)
    canv.SetRightMargin(0.15)
    # canv.SetLogz()
    # canv.SetLogy()
    h2d_copy = h2d.Clone()
    # h2d_copy.Scale(1, "width")
    h2d_copy.Draw("COLZ")
    canv.SetLogx()
    h2d_copy.GetXaxis().SetMoreLogLabels()
    canv.SaveAs(output_filename)

    zoom_ymin, zoom_ymax = 0.1, 5

    h2d_copy.SetAxisRange(zoom_ymin, zoom_ymax,"Y")
    h2d_copy.SetAxisRange(xmin, xmax, "X")
    canv.SaveAs(output_filename.replace(".pdf", "_zoomY.pdf"))
    
    canv.SetLogz()
    canv.SaveAs(output_filename.replace(".pdf", "_zoomY_logZ.pdf"))

    canv.SetLogz(False)
    # h2d.Scale(1, "width")
    h2d_normed = cu.make_normalised_TH2(h2d, norm_axis='x', recolour=True)
    h2d_normed.Draw("COLZ")
    h2d_normed.GetXaxis().SetMoreLogLabels()
    # h2d_normed.SetMinimum(1E-5)
    h2d_normed.SetAxisRange(xmin, xmax, "X")
    canv.SaveAs(output_filename.replace(".pdf", "_normX.pdf"))
    
    h2d_normed.SetAxisRange(zoom_ymin, zoom_ymax,"Y")
    canv.SaveAs(output_filename.replace(".pdf", "_normX_zoomY.pdf"))

    # Do cumulative plot per column (ie fraction of events passing cut < y)
    h2d_ndarray_cumsum = h2d_ndarray.cumsum(axis=0)
    nonzero_mask = h2d_ndarray_cumsum[-1] > 0
    h2d_ndarray_cumsum[:, nonzero_mask] /= h2d_ndarray_cumsum[-1][nonzero_mask] # scale so total is 1
    
    h2d_cumsum = cu.ndarray_to_th2(h2d_ndarray_cumsum, binsx=xbins, binsy=ybins)
    # Get max row ind
    max_filled_row_ind = int(h2d_ndarray_cumsum.argmax(axis=0).max())
    h2d_cumsum.GetYaxis().SetRange(1, max_filled_row_ind+1)  # +1 as ROOT 1-indexed

    # ROOT.gStyle.SetPalette(ROOT.kBird)
    ylabel = "Fraction of events with " + ylabel + " < y"
    if "unweighted" in histname:
        h2d_cumsum.SetTitle("Unweighted;%s;%s" % (xlabel, ylabel))
    else:
        h2d_cumsum.SetTitle("Weighted;%s;%s" % (xlabel, ylabel))
    canv.Clear()
    canv.SetLogz(False)

    h2d_cumsum.SetContour(20)
    h2d_cumsum.Draw("CONT1Z")
    h2d_cumsum.SetAxisRange(xmin, xmax, "X")
    canv.SetLogx()
    h2d_cumsum.GetXaxis().SetMoreLogLabels()
    canv.SaveAs(output_filename.replace(".pdf", "_cumulY.pdf"))

    h2d_cumsum.SetAxisRange(zoom_ymin, zoom_ymax,"Y")
    canv.SaveAs(output_filename.replace(".pdf", "_cumulY_zoomY.pdf"))
    canv.Clear()

    h2d_normed.Draw("COL")
    h2d_cumsum.Draw("CONT1Z SAME")
    h2d_cumsum.SetAxisRange(xmin, xmax, "X")
    canv.SetLogx()
    h2d_cumsum.GetXaxis().SetMoreLogLabels()
    canv.SaveAs(output_filename.replace(".pdf", "_cumulY_normX.pdf"))

    h2d_cumsum.SetAxisRange(zoom_ymin, zoom_ymax,"Y")
    canv.SaveAs(output_filename.replace(".pdf", "_cumulY_normX_zoomY.pdf"))

    tf.Close()


def do_weight_vs_var_plot(histname, input_filename, output_filename):
    ROOT.gStyle.SetPalette(palette_2D)
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    if h3d.GetEntries() == 0:
        return
    h2d = h3d.Project3D("xz")
    if "unweighted" in histname:
        h2d.SetTitle("Unweighted")
    else:
        h2d.SetTitle("Weighted")
    h2d.GetXaxis().SetTitle(get_var_str(histname))

    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetLeftMargin(0.1)
    canv.SetRightMargin(0.15)
    canv.SetLogz()
    h2d.Draw("COLZ")
    canv.SaveAs(output_filename)
    tf.Close()


def do_weight_vs_var_plot_per_pt(histname, input_filename, output_filename):
    ROOT.gStyle.SetPalette(palette_2D)
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    if h3d.GetEntries() == 0:
        return

    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetLeftMargin(0.1)
    canv.SetRightMargin(0.15)
    canv.SetLogz()

    weight_str = "(unweighted)" if "unweighted" in histname else "(weighted)"
    for ibin in range(1, h3d.GetNbinsY()+1):
        h3d.GetYaxis().SetRange(ibin, ibin+1)
        h2d = h3d.Project3D("xz")
        if h2d.GetEntries() == 0:
            continue
        pt_low = h3d.GetYaxis().GetBinLowEdge(ibin)
        pt_high = h3d.GetYaxis().GetBinLowEdge(ibin+1)
        jet_str = pt_genjet_str if "_vs_pt_genjet_vs_" in histname else pt_str
        h2d.SetTitle("%g < %s < %g GeV %s" % (pt_low, jet_str, pt_high, weight_str))
        h2d.GetXaxis().SetTitle(get_var_str(histname))
        h2d.Draw("COLZ")
        jet_app = "_genjet" if "_vs_pt_genjet_vs" in histname else ""
        this_output_filename = output_filename.replace(".pdf", "_pt%s%gto%g.pdf" % (jet_app, pt_low, pt_high))
        canv.SaveAs(this_output_filename)
        canv.Clear()

    tf.Close()


def do_weight_vs_pt_plot(input_filename, output_filename):
    ROOT.gStyle.SetPalette(palette_2D)
    histname = "Weight_Presel/weight_vs_pt_vs_pt_jet_qScale_ratio"
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    if h3d.GetEntries() == 0:
        return
    h2d = h3d.Project3D("xy")

    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetLeftMargin(0.1)
    canv.SetRightMargin(0.15)
    canv.SetLogz()
    canv.SetLogx()
    h2d.Draw("COLZ")
    canv.SaveAs(output_filename)
    tf.Close()


def do_weight_vs_genjet_pt_plot(input_filename, output_filename):
    ROOT.gStyle.SetPalette(palette_2D)
    histname = "Weight_Presel/weight_vs_pt_vs_pt_genjet_qScale_ratio"
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    if h3d.GetEntries() == 0:
        return

    h2d = h3d.Project3D("xy")

    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetLeftMargin(0.1)
    canv.SetRightMargin(0.15)
    canv.SetLogz()
    canv.SetLogx()
    h2d.Draw("COLZ")
    canv.SaveAs(output_filename)
    tf.Close()


def do_jet_pt_with_var_cuts(histname, cuts, input_filename, output_filename):
    ROOT.gStyle.SetPalette(palette_1D)
    total = len(cuts) -1 + .1
    if len(cuts) <= 3:
        # ROOT.gStyle.SetPalette(ROOT.kCool)
        num_colours = ROOT.TColor.GetPalette().fN - 1
        print('num_colours:', num_colours)
        for index in range(len(cuts)):
            print(num_colours, index, len(cuts), index / len(cuts), num_colours * index / total)
            print(index, ROOT.TColor.GetColorPalette(int(num_colours * 1. * index / total)))
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    if h3d.GetEntries() == 0:
        return
    pt_hists = []
    for cut in cuts:
        max_bin = h3d.GetZaxis().FindFixBin(cut)
        # print("cut:", cut, "bin:", max_bin)
        h = h3d.ProjectionY("pt_var_lt_%g" % cut, 0, -1, 0, max_bin, "e")
        h2 = h.Clone()
        h2.Rebin(2)
        if h.GetEntries() > 0:
            h3 = qgp.hist_divide_bin_width(h2)
        pt_hists.append(h3)

    line_styles = [1, 2, 3]
    if len(cuts) <= 3:
        line_styles = [1]
    n_line_styles = len(line_styles)
    ref_ind = 0
    conts = [Contribution(h, label=" < %g" % cut,
                          line_color=cu.get_colour_seq(ind, total),
                          line_style=line_styles[ind % n_line_styles],
                          line_width=2,
                          marker_color=cu.get_colour_seq(ind, total),
                          subplot=pt_hists[ref_ind] if ind != ref_ind else None)
             for ind, (h, cut) in enumerate(zip(pt_hists, cuts))]

    jet_str = pt_genjet_str if "_vs_pt_genjet_vs_" in histname else pt_str
    weight_str = "(unweighted)" if "unweighted" in histname else "(weighted)"
    ratio_lims = (0.5, 2.5)
    ratio_lims = None
    plot = Plot(conts, what='hist',
                title='%s for cuts on %s %s' % (jet_str, get_var_str(histname), weight_str),
                xtitle=None,
                ytitle='N',
                # xlim=None, ylim=None,
                legend=True,
                subplot_type='ratio',
                subplot_title='* / var < %g' % cuts[ref_ind],
                subplot_limits=ratio_lims,
                has_data=False)
    plot.y_padding_max_log = 200
    plot.subplot_maximum_ceil = 4
    plot.subplot_maximum_floor = 1.02
    plot.subplot_minimum_ceil = 0.98
    plot.legend.SetY1(0.7)
    plot.legend.SetY2(0.89)
    plot.legend.SetX1(0.78)
    plot.legend.SetX2(0.88)
    plot.plot("NOSTACK HISTE", "NOSTACK HIST")
    plot.set_logx(True, do_more_labels=True)
    plot.set_logy(True, do_more_labels=False)
    plot.save(output_filename)


def do_jet_pt_rel_error_with_var_cuts(histname, cuts, input_filename, output_filename):
    ROOT.gStyle.SetPalette(palette_1D)
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    if h3d.GetEntries() == 0:
        return
    pt_hists = []
    for cut in cuts:
        max_bin = h3d.GetZaxis().FindFixBin(cut)
        # print("cut:", cut, "bin:", max_bin)
        h = h3d.ProjectionY("pt_var_lt_%g" % cut, 0, -1, 0, max_bin, "e")
        h2 = h.Clone()
        h2.Rebin(2)
        if h.GetEntries() > 0:
            h3 = qgp.hist_divide_bin_width(h2)
        # convert bin contents to bin error/bin contents
        for ibin in range(1, h2.GetNbinsX()+1):
            if h3.GetBinContent(ibin) == 0:
                continue
            h3.SetBinContent(ibin, h3.GetBinError(ibin) / h3.GetBinContent(ibin))
            h3.SetBinError(ibin, 0)
        pt_hists.append(h3)

    line_styles = [1, 2, 3]
    n_line_styles = len(line_styles)
    conts = [Contribution(h, label=" < %g" % cut,
                          line_color=cu.get_colour_seq(ind, len(cuts)),
                          line_style=line_styles[ind % n_line_styles],
                          line_width=2,
                          marker_color=cu.get_colour_seq(ind, len(cuts)),
                          subplot=pt_hists[-1])
             for ind, (h, cut) in enumerate(zip(pt_hists, cuts))]

    jet_str = pt_genjet_str if "_vs_pt_genjet_vs_" in histname else pt_str
    weight_str = "(unweighted)" if "unweighted" in histname else "(weighted)"
    ratio_lims = (0.98, 1.02) if "unweighted" in histname else None
    plot = Plot(conts, what='hist',
                title='%s for cuts on %s %s' % (jet_str, get_var_str(histname), weight_str),
                xtitle=None,
                ytitle='Relative error',
                # xlim=None, ylim=None,
                legend=True,
                subplot_type='ratio',
                subplot_title='* / var < %g' % cuts[-1],
                subplot_limits=ratio_lims,
                has_data=False)
    plot.y_padding_max_log = 200
    plot.subplot_maximum_ceil = 2
    plot.subplot_maximum_floor = 1.02
    plot.subplot_minimum_ceil = 0.98
    plot.legend.SetY1(0.7)
    plot.legend.SetY2(0.89)
    plot.legend.SetX1(0.78)
    plot.legend.SetX2(0.88)
    plot.plot("NOSTACK HISTE", "NOSTACK HIST")
    plot.set_logx(True, do_more_labels=True)
    plot.set_logy(True, do_more_labels=False)
    plot.save(output_filename)


def do_cut_scan_per_pt(histname, input_filename, output_filename):
    ROOT.gStyle.SetPalette(palette_1D)
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    h3d_unweighted = cu.get_from_tfile(tf, histname+"_unweighted")
    if h3d.GetEntries() == 0:
        return

    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetLeftMargin(0.12)
    canv.SetRightMargin(0.12)
    # canv.SetLogz()
    h2d = h3d.Project3D("zy") # var vs pt
    h2d_unweighted = h3d_unweighted.Project3D("zy") # var vs pt

    var_name = os.path.basename(histname).replace("weight_vs_pt_vs_", "").replace("weight_vs_pt_genjet_vs_", "")

    for ibin in range(3, h2d.GetNbinsX()+1):  # iterate over pt bins
        if h2d.Integral(ibin, ibin+1, 0, -1) == 0:
            # data.append(None)
            continue

        pt_low = h3d.GetYaxis().GetBinLowEdge(ibin)
        pt_high = h3d.GetYaxis().GetBinLowEdge(ibin+1)

        data = []
        data_unweighted = []
        # Do integral, error for increasingly looser cuts
        # find maximum in this pt bin
        # yes I probably should collapse to a 1D hist and use GetMaximumBin
        max_val, max_bin = 0, 0
        for icut in range(1, h2d.GetNbinsY()+1):
            val =  h2d.GetBinContent(ibin, icut)
            if val > max_val:
                max_val = val
                max_bin = icut

        for icut in range(max_bin + 5, h2d.GetNbinsY()+2):
            err = array('d', [0])
            count = h2d.IntegralAndError(ibin, ibin+1, 1, icut-1, err)
            if count == 0:
                continue
            data.append([count, err[0], h2d.GetYaxis().GetBinLowEdge(icut)])

            err = array('d', [0])
            count = h2d_unweighted.IntegralAndError(ibin, ibin+1, 1, icut, err)
            data_unweighted.append([count, err[0], h2d.GetYaxis().GetBinLowEdge(icut)])

        # Plot count, rel error vs cut value
        cuts = [d[2] for d in data]
        gr_count = ROOT.TGraph(len(data), array('d', cuts), array('d', [d[0] / data[-1][0] for d in data]))

        gr_count.SetMarkerColor(ROOT.kRed)
        gr_count.SetMarkerStyle(22)
        gr_count.SetLineColor(ROOT.kRed)
        gr_count.SetTitle("%g < p_{T} < %g GeV;%s cut (<);Count (relative to loosest cut)" % (pt_low, pt_high, get_var_str(histname)))

        gr_count_unweighted = ROOT.TGraph(len(data), array('d', cuts), array('d', [d[0] / data_unweighted[-1][0] for d in data_unweighted]))
        gr_count_unweighted.SetMarkerColor(ROOT.kBlack)
        gr_count_unweighted.SetMarkerStyle(23)
        gr_count_unweighted.SetLineColor(ROOT.kBlack)
        gr_count_unweighted.SetTitle("%g < p_{T} < %g GeV;%s cut (<);Count (relative to loosest cut)" % (pt_low, pt_high, get_var_str(histname)))

        leg = ROOT.TLegend(0.7, 0.5, 0.85, 0.65)
        leg.AddEntry(gr_count, "Weighted", "LP")
        leg.AddEntry(gr_count_unweighted, "Unweighted", "LP")

        # gr_rel_err = ROOT.TGraph(len(data), array('d', cuts), array('d', [(d[1] / d[0]) if d[0] != 0 else 0 for d in data ]))
        # gr_rel_err.SetMarkerColor(ROOT.kRed)
        # gr_rel_err.SetMarkerStyle(22)
        # gr_rel_err.SetLineColor(ROOT.kRed)
        # gr_rel_err.SetTitle("%g < p_{T} < %g GeV;%s cut (<);Rel. error" % (pt_low, pt_high, var_name))

        canv.SetLogy(False)

        gr_count.Draw("ALP")
        gr_count_unweighted.Draw("LP")
        gr_count.Draw("LP")
        leg.Draw()
        canv.SaveAs(output_filename.replace(".pdf", "_count_pt%gto%g.pdf" % (pt_low, pt_high)))

        # canv.Clear()
        # gr_rel_err.Draw("ALP")
        # canv.SetLogy()
        # gr_rel_err.GetYaxis().SetMoreLogLabels()
        # canv.SaveAs(output_filename.replace(".pdf", "_rel_err_pt%gto%g.pdf" % (pt_low, pt_high)))


    tf.Close()


def do_cut_roc_per_pt(histname, input_filename, output_filename):
    """Plot fractional # unweighted vs fraction # weighted, for different cuts
    Not a true ROC, but kinda like one
    """
    ROOT.gStyle.SetPalette(palette_1D)
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    h3d_unweighted = cu.get_from_tfile(tf, histname+"_unweighted")
    if h3d.GetEntries() == 0:
        return

    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetLeftMargin(0.12)
    canv.SetRightMargin(0.12)
    # canv.SetLogz()
    h2d = h3d.Project3D("zy") # var vs pt
    h2d_unweighted = h3d_unweighted.Project3D("zy") # var vs pt

    var_name = os.path.basename(histname).replace("weight_vs_pt_vs_", "").replace("weight_vs_pt_genjet_vs_", "")

    for ibin in range(3, h2d.GetNbinsX()+1):  # iterate over pt bins
        if h2d.Integral(ibin, ibin+1, 0, -1) == 0:
            # data.append(None)
            continue

        pt_low = h3d.GetYaxis().GetBinLowEdge(ibin)
        pt_high = h3d.GetYaxis().GetBinLowEdge(ibin+1)

        data = []
        data_unweighted = []
        # Do integral, error for increasingly looser cuts
        # find maximum in this pt bin
        # yes I probably should collapse to a 1D hist and use GetMaximumBin
        max_val, max_bin = 0, 0
        for icut in range(1, h2d.GetNbinsY()+1):
            val =  h2d.GetBinContent(ibin, icut)
            if val > max_val:
                max_val = val
                max_bin = icut

        for icut in range(max_bin + 5, h2d.GetNbinsY()+2, 2):
        # for icut in range(2, h2d.GetNbinsY()+2):
            err = array('d', [0])
            count = h2d.IntegralAndError(ibin, ibin+1, 1, icut-1, err)
            if count == 0:
                continue
            data.append([count, err[0], h2d.GetYaxis().GetBinLowEdge(icut)])

            err = array('d', [0])
            count = h2d_unweighted.IntegralAndError(ibin, ibin+1, 1, icut-1, err)
            data_unweighted.append([count, err[0], h2d.GetYaxis().GetBinLowEdge(icut)])

        cuts = np.array([d[2] for d in data][1:])
        # cuts = np.array([d[2] for d in data])

        weighted_fractions = np.array([abs(d[0]-dd[0]) / dd[0] for d, dd in zip(data[:-1], data[1:])])
        unweighted_fractions = np.array([abs(d[0]-dd[0]) / dd[0] for d, dd in zip(data_unweighted[:-1], data_unweighted[1:])])

        non_zero_mask = (unweighted_fractions>0) & (weighted_fractions>0)
        non_zero_weighted = weighted_fractions[non_zero_mask]
        weight_min_pow = math.floor(math.log10(min(non_zero_weighted))) if len(non_zero_weighted) > 0 else -10
        weight_max_pow = math.floor(math.log10(max(non_zero_weighted))) if len(non_zero_weighted) > 0 else 0
        assert(weight_max_pow>=weight_min_pow)

        non_zero_unweighted = unweighted_fractions[non_zero_mask]
        unweight_min_pow = math.floor(math.log10(min(non_zero_unweighted))) if len(non_zero_unweighted) > 0 else -10
        unweight_max_pow = math.floor(math.log10(max(non_zero_unweighted))) if len(non_zero_unweighted) > 0 else 0
        assert(unweight_max_pow>=unweight_min_pow)

        mask = unweighted_fractions < 10**(unweight_min_pow+1)  # last decade of unweighted drops
        mask &= weighted_fractions > 10**(weight_max_pow-1)  # largest decades of weighted drops

        if np.sum(mask) == 0:
            continue

        # weighted_fractions = np.array([d[0] / data[-1][0] for d in data])
        # unweighted_fractions = np.array([d[0] / data_unweighted[-1][0] for d in data_unweighted])

        unweighted_useful = unweighted_fractions[mask & non_zero_mask]
        weighted_useful = weighted_fractions[mask & non_zero_mask]
        if "pt_jet_genHT_ratio" in histname and pt_low == 800:
            print("weight_min_pow:", weight_min_pow)
            print("weight_max_pow:", weight_max_pow)
            print("unweight_min_pow:", unweight_min_pow)
            print("unweight_max_pow:", unweight_max_pow)
            print("unweight_max_pow:", unweight_max_pow)
            print("weighted_useful:", weighted_useful)
            print("unweighted_useful:", unweighted_useful)

        gr_count = ROOT.TGraph(len(unweighted_useful), unweighted_useful, weighted_useful)
        gr_count.SetMarkerColor(ROOT.kRed)
        gr_count.SetMarkerSize(0)
        gr_count.SetMarkerStyle(21)
        gr_count.SetLineColor(ROOT.kRed)
        gr_count.SetTitle("%s, %g < p_{T} < %g GeV;Relative unweighted count;Relative weighted count" % (get_var_str(histname), pt_low, pt_high))
        gr_count.SetTitle("%s, %g < p_{T} < %g GeV;Unweighted fractional drop;Weighted fractional drop" % (get_var_str(histname), pt_low, pt_high))


        # add annotations of cuts
        latexs = []
        for i, cut in enumerate(cuts[mask * non_zero_mask]):
            latex = ROOT.TLatex(gr_count.GetX()[i], gr_count.GetY()[i], " < %.2f" % cut)
            latex.SetTextSize(0.02)
            latex.SetTextColor(ROOT.kBlue)
            gr_count.GetListOfFunctions().Add(latex)
            latexs.append(latex)

        # canv.SetLogx(False)
        # canv.SetLogy(False)

        # ROOT.TGaxis.SetMaxDigits(2)

        # gr_count.Draw("ALP")

        # ROOT.TGaxis.SetMaxDigits(2)
        # unweighted_min = 0.9999

        # # Calculate differences between points
        # unweighted_diffs = unweighted_fractions[1:] - unweighted_fractions[:-1]
        # weighted_diffs = weighted_fractions[1:] - weighted_fractions[:-1]
        # big_diff_inds = []
        # for ind, (u, w) in enumerate(zip(unweighted_diffs, weighted_diffs)):
        #     # look for big diff in weighted frac, small diff in unweighted,
        #     # with a limit on the minimum size of unweighted frac
        #     # (only trying to remove a few events)
        #     if u > 0 and w / u > 100 and u < 0.005 and unweighted_fractions[ind] > unweighted_min:
        #         big_diff_inds.append(ind)

        # if "pt_jet_genHT_ratio" in histname and pt_low == 186:
        #     for u, w in zip(unweighted_diffs, weighted_diffs):
        #         print(u, w)
        #     print(big_diff_inds)

        # make graph of big diff points, add annotations of cuts
        # if len(big_diff_inds) > 0:
        #     gr_big_diffs = ROOT.TGraph(len(big_diff_inds), array('d', [unweighted_fractions[i+1] for i in big_diff_inds]), array('d', [weighted_fractions[i+1] for i in big_diff_inds]))
        #     gr_big_diffs.SetLineWidth(0)
        #     gr_big_diffs.SetMarkerColor(ROOT.kBlue)
        #     gr_big_diffs.SetMarkerStyle(25)
        #     latexs = []
        #     for i, ind in enumerate(big_diff_inds[:]):
        #         latex = ROOT.TLatex(gr_big_diffs.GetX()[i], gr_big_diffs.GetY()[i], " < %.2f" % cuts[ind+1])
        #         latex.SetTextSize(0.02)
        #         latex.SetTextColor(ROOT.kBlue)
        #         gr_big_diffs.GetListOfFunctions().Add(latex)
        #         latexs.append(latex)
        #     gr_big_diffs.Draw("*")

        # gr_count.GetXaxis().SetLimits(unweighted_min, 1)

        # find corresponding value for weighted to set axis range
        # weighted_min = 0
        # for ind, u in enumerate(unweighted_fractions):
        #     if u >= unweighted_min:
        #         weighted_min = weighted_fractions[ind-1]
        #         if ind == len(unweighted_fractions) - 1:
        #             weighted_min = 0
        #         break
        # gr_count.GetHistogram().SetMinimum(weighted_min*1.1 - 0.1)
        # gr_count.GetHistogram().SetMaximum(1)
        # canv.SaveAs(output_filename.replace(".pdf", "_count_pt%gto%g.pdf" % (pt_low, pt_high)))

        # do a version zoomed out
        canv.Clear()
        gr_count.SetMarkerSize(0.5)
        gr_count.Draw("AP")

        # unweighted_min = 0.
        # gr_count.GetXaxis().SetLimits(unweighted_min, 1)
        # weighted_min = 0
        # for ind, u in enumerate(unweighted_fractions):
        #     if u >= unweighted_min:
        #         weighted_min = weighted_fractions[ind-1]
        #         if ind == len(unweighted_fractions) - 1:
        #             weighted_min = 0
        #         break
        # gr_count.GetHistogram().SetMinimum(weighted_min*1.1 - 0.1)

        gr_count.GetXaxis().SetMoreLogLabels()
        gr_count.GetYaxis().SetMoreLogLabels()

        weight_min_pow = math.floor(math.log10(min(weighted_useful))) if len(weighted_useful) > 0 else -10
        weight_max_pow = math.floor(math.log10(max(weighted_useful))) if len(weighted_useful) > 0 else 0

        unweight_min_pow = math.floor(math.log10(min(unweighted_useful))) if len(unweighted_useful) > 0 else -10
        unweight_max_pow = math.floor(math.log10(max(unweighted_useful))) if len(unweighted_useful) > 0 else 0
        gr_count.GetHistogram().SetMinimum(10**weight_min_pow)
        gr_count.GetHistogram().SetMaximum(10**(weight_max_pow+1))
        gr_count.GetXaxis().SetLimits(10**unweight_min_pow, 10**(unweight_max_pow+1))
        canv.SetLogy()
        canv.SetLogx()
        canv.SaveAs(output_filename.replace(".pdf", "_count_pt%gto%g_all.pdf" % (pt_low, pt_high)))


    tf.Close()


if __name__ == "__main__":
    MAINDIR="/Volumes/Extreme SSD/Projects/QGAnalysis"
    histnames = [
        # "weight_vs_pt_vs_pt_jet_genHT_ratio",
        # "weight_vs_pt_genjet_vs_pt_jet_genHT_ratio",
        
        # "weight_vs_pt_vs_pt_genjet_genHT_ratio",
        # "weight_vs_pt_genjet_vs_pt_genjet_genHT_ratio",

        # "weight_vs_pt_vs_pt_jet_qScale_ratio",
        # "weight_vs_pt_genjet_vs_pt_jet_qScale_ratio",
        
        # "weight_vs_pt_vs_pt_genjet_qScale_ratio",
        # "weight_vs_pt_genjet_vs_pt_genjet_qScale_ratio",

        "weight_vs_pt_vs_pt_jet_ptHat_ratio",
        "weight_vs_pt_genjet_vs_pt_jet_ptHat_ratio",
        
        "weight_vs_pt_vs_pt_genjet_ptHat_ratio",
        "weight_vs_pt_genjet_vs_pt_genjet_ptHat_ratio",

        # "weight_vs_pt_vs_PU_ptHat_genHT_ratio",
        # "weight_vs_pt_genjet_vs_PU_ptHat_genHT_ratio",
        
        "weight_vs_pt_vs_PU_ptHat_ptHat_ratio",
        "weight_vs_pt_genjet_vs_PU_ptHat_ptHat_ratio",
        
        # "weight_vs_pt_vs_PU_ptHat_jetkT_ratio",
        # "weight_vs_pt_genjet_vs_PU_ptHat_jetkT_ratio",
        
        # "weight_vs_pt_vs_pt_jet_jetkT_ratio",
        # "weight_vs_pt_genjet_vs_pt_jet_jetkT_ratio",
    ]

    histnames_unweighted = [h + "_unweighted" for h in histnames]
    histnames.extend(histnames_unweighted)

    hist_dirs = ['Weight_Presel', 'Weight_Reco_sel'][1:]

    workdirs = [
        # "workdir_102X_v3data_v2mc_ak4puppi_data_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_weightHists",
        # "workdir_102X_v3data_v2mc_ak8puppi_data_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_weightHists",
        # "workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_weightHistsBig_noWeightCuts",
        # "workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_weightHists_PUWeightCuts",
        # "workdir_102X_v2_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_PUWeightCuts",
        # "workdir_102X_v2_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts",
        "workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning6_noPUcuts",
        "workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning6_onlyPUPtHatCut",
    ]


    qcd_filenames = [
        # qgc.QCD_FILENAME,
        # qgc.QCD_PYTHIA_ONLY_FILENAME,
        qgc.QCD_HERWIG_FILENAME,
    ]

    dy_filenames = [
        qgc.DY_FILENAME,
        qgc.DY_HERWIG_FILENAME,
        # 'uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL_JetKtMin170_PartonKtMin300.root'
        # 'uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_DYJetsToLL_M-50_HT-0to70.root'
    ]

    # ROOT.gStyle.SetPalette(ROOT.kPastel)

    for workdir in workdirs:
        print(".Doing", workdir)
        # for filename in dy_filenames:
        for filename in qcd_filenames:
        # for filename in chain(qcd_filenames, dy_filenames):

            ifile = os.path.join(MAINDIR, workdir, filename)
            if not os.path.isfile(ifile):
                continue

            for hist_dir in hist_dirs:
                print("...Doing", hist_dir)

                dir_append = hist_dir.lower().replace("weight_", "")
                odir = os.path.join(MAINDIR, workdir, "var_jet_pt_weights_%s_%s" % (dir_append, filename.replace("uhh2.AnalysisModuleRunner.MC.MC_", "").replace(".root", "")))
                cu.check_dir_exists_create(odir)

                # do_weight_vs_pt_plot(input_filename=ifile,
                #                      output_filename=os.path.join(odir, "weight_vs_pt.pdf"))

                # do_weight_vs_genjet_pt_plot(input_filename=ifile,
                #                             output_filename=os.path.join(odir, "weight_vs_genjet_pt.pdf"))

                for ind, base_hname in enumerate(histnames):
                    hname = os.path.join(hist_dir, base_hname)
                    
                    tf = cu.open_root_file(ifile)
                    if not cu.exists_in_file(tf, hname):
                        continue

                    do_var_vs_pt_plot(histname=hname,
                                      input_filename=ifile,
                                      output_filename=os.path.join(odir, base_hname.replace("weight_vs_", "") + ".pdf"))

                    # if "_vs_pt_" in hname:
                    #     # skip the _vs_pt_genjet_ since identical, and we collapse that axis anyway
                    #     do_weight_vs_var_plot(histname=hname,
                    #                           input_filename=ifile,
                    #                           output_filename=os.path.join(odir, base_hname.replace("_vs_pt_vs_", "_").replace("_vs_pt_genjet_vs_", "") + ".pdf"))

                    cuts = [10, 9, 8, 7, 6, 5, 4, 3, 2]
                    cuts = [10, 8, 6, 4, 3, 2]
                    cuts = [100, 75, 50, 25, 10, 8, 5, 2]
                    cuts = [25, 10, 8, 5, 4, 3, 2]
                    if "PU_ptHat" in hname:
                        cuts = [25, 1, 0.7]
                    if 'pt_genjet_genHT_ratio' in hname:
                        cuts.extend([1.5, 1])
                    if 'pt_jet_genHT_ratio' in hname:
                        cuts = [5, 2, 1.5, 1, 0.9, 0.8, 0.75]
                    do_jet_pt_with_var_cuts(histname=hname,
                                            cuts=cuts,
                                            input_filename=ifile,
                                            output_filename=os.path.join(odir, base_hname.replace("weight_vs_pt_vs_", "pt_cuts_").replace("weight_vs_pt_genjet_vs_", "pt_genjet_cuts_") + ".pdf"))

                    do_jet_pt_rel_error_with_var_cuts(histname=hname,
                                                      cuts=cuts,
                                                      input_filename=ifile,
                                                      output_filename=os.path.join(odir, base_hname.replace("weight_vs_pt_vs_", "pt_cuts_rel_error_").replace("weight_vs_pt_genjet_vs_", "pt_genjet_rel_error_") + ".pdf"))

                    # do_weight_vs_var_plot_per_pt(histname=hname,
                    #                              input_filename=ifile,
                    #                              output_filename=os.path.join(odir, base_hname.replace("_vs_pt_vs_", "_").replace("_vs_pt_genjet_vs_", "_") + ".pdf"))

                    # if "unweighted" not in hname and "_vs_pt_genjet_vs_" not in hname:
                    #     do_cut_scan_per_pt(histname=hname,
                    #                        input_filename=ifile,
                    #                        output_filename=os.path.join(odir, base_hname.replace("weight_vs_pt_vs_", "cuts_scan_") + ".pdf"))
                        # do_cut_roc_per_pt(histname=hname,
                        #                   input_filename=ifile,
                        #                   output_filename=os.path.join(odir, base_hname.replace("weight_vs_pt_vs_", "cuts_roc_") + ".pdf"))

                    tf.Close()
