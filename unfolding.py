#!/usr/bin/env python


"""TUnfold it all

Thanks to Ashley, Dennis
"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse
from array import array
import numpy as np
import math
import distutils
from distutils import util

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from unfolding_classes import TauScanner, LCurveScanner, MyUnfolder, MyUnfolderPlotter


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat(".3f")


# Control plot output format
OUTPUT_FMT = "pdf"


def calculate_chi2(hist_test, hist_truth, hist_covariance=None):
    """Calculate chi2 = (test - truth)T inv(covariance) (test - truth)

    Parameters
    ----------
    hist_test : TYPE
        Description
    hist_truth : TYPE
        Description
    hist_covariance : TH2, optional
        If None, assumes covariance = diagonals of hist_truth
    """
    pass
    # test_vector, test_err = th1_to_ndarray(hist_test)
    # truth_vector, truth_err = th1_to_ndarray(hist_truth)
    # covariance_matrix = th2_to_ndarray(hist_covariance)


def plot_simple_unfolded(unfolded, tau, reco, gen, fake, output_filename, title=""):
    """Simple plot of unfolded, reco, gen, by bin number (ie non physical axes)"""
    entries = []

    if reco:
        entries.append(
            Contribution(reco, label="Reco",
                         line_color=ROOT.kGreen+2, line_width=1,
                         marker_color=ROOT.kGreen+2, marker_size=0,
                         normalise_hist=False),
        )

    if gen:
        entries.append(
            Contribution(gen, label="Gen",
                         line_color=ROOT.kBlue, line_width=1,
                         marker_color=ROOT.kBlue, marker_size=0,
                         normalise_hist=False),
        )

    if fake:
        entries.append(
            Contribution(fake, label="Fakes",
                         line_color=ROOT.kOrange+4, line_width=1,
                         marker_color=ROOT.kOrange+4, marker_size=0,
                         normalise_hist=False),
        )

    if unfolded:
        entries.append(
            Contribution(unfolded, label="Unfolded (#tau = %.3g)" % (tau),
                         line_color=ROOT.kRed, line_width=0,
                         marker_color=ROOT.kRed, marker_size=0.6, marker_style=20,
                         normalise_hist=False, subplot=gen),
        )

    plot = Plot(entries,
                what='hist',
                title=title,
                xtitle="Bin number",
                ytitle="N",
                subplot_type='ratio',
                subplot_title='Data / MC',
                subplot_limits=(0.8, 1.2))
    plot.default_canvas_size = (800, 600)
    plot.plot("NOSTACK HISTE")
    plot.main_pad.SetLogy(1)
    ymax = max(h.GetMaximum() for h in [reco, gen, fake, unfolded] if h)
    plot.container.SetMaximum(ymax * 100)
    ymin = min(h.GetMinimum(1E-8) for h in [reco, gen, fake, unfolded] if h)
    plot.container.SetMinimum(ymin*0.1)
    plot.legend.SetY1NDC(0.77)
    plot.legend.SetX1NDC(0.65)
    plot.legend.SetX2NDC(0.88)
    plot.save(output_filename)


def plot_simple_detector(reco_data, reco_data_fake, reco_mc, reco_mc_fake, output_filename, title):
    """Plot detector-level quantities for data & MC, by bin number (ie non physical axes)"""
    entries = []

    if reco_mc:
        entries.append(
            Contribution(reco_mc, label="MC [detector-level]",
                         line_color=ROOT.kGreen+2, line_width=1,
                         marker_color=ROOT.kGreen+2, marker_size=0,
                         normalise_hist=False),
        )

    if reco_mc_fake:
        entries.append(
            Contribution(reco_mc_fake, label="MC fakes [detector-level]",
                         line_color=ROOT.kOrange+4, line_width=1,
                         marker_color=ROOT.kOrange+4, marker_size=0,
                         normalise_hist=False),
        )

    if reco_data_fake:
        entries.append(
            Contribution(reco_data_fake, label="Data fakes template [detector-level]",
                         line_color=ROOT.kMagenta+2, line_width=1,
                         marker_color=ROOT.kMagenta+2, marker_size=0,
                         normalise_hist=False),
        )

    if reco_data:
        entries.append(
            Contribution(reco_data, label="Data [detector-level]",
                         line_color=ROOT.kRed, line_width=0,
                         marker_color=ROOT.kRed, marker_size=0.6, marker_style=20,
                         normalise_hist=False, subplot=reco_mc),
        )

    plot = Plot(entries,
                what='hist',
                title=title,
                xtitle="Bin number",
                ytitle="N",
                subplot_type='ratio',
                subplot_title='Data / MC',
                subplot_limits=(0.8, 1.2))
    plot.default_canvas_size = (800, 600)
    plot.plot("NOSTACK HISTE")
    plot.main_pad.SetLogy(1)
    ymax = max(h.GetMaximum() for h in [reco_mc, reco_mc_fake, reco_data] if h)
    plot.container.SetMaximum(ymax * 100)
    ymin = min(h.GetMinimum(1E-10) for h in [reco_mc, reco_mc_fake, reco_data] if h)
    plot.container.SetMinimum(ymin*0.1)

    # plot.container.SetMinimum(0.001)
    plot.legend.SetY1NDC(0.77)
    plot.legend.SetX1NDC(0.65)
    plot.legend.SetX2NDC(0.88)
    plot.save(output_filename)


def plot_fake_fraction(hist, output_filename, title):
    entries = [Contribution(hist,
                            label='From MC',
                            line_color=ROOT.kBlue, line_width=1,
                            marker_color=ROOT.kBlue, marker_size=0)]
    plot = Plot(entries,
                what='hist',
                title=title,
                xtitle='Bin number',
                ytitle='Fake fraction')
    plot.default_canvas_size = (800, 600)
    plot.plot('NOSTACK HIST')
    plot.set_logy()
    plot.container.SetMaximum(hist.GetMaximum() * 100)
    plot.container.SetMinimum(hist.GetMinimum(1E-12) * 0.1)
    plot.save(output_filename)


def create_hist_with_errors(hist, err_matrix):
    hnew = hist.Clone(cu.get_unique_str())
    nbins = hist.GetNbinsX()
    for i in range(1, nbins+1):
        err = math.sqrt(err_matrix.GetBinContent(i, i))
        hnew.SetBinError(i, err)
    return hnew


def make_hist_from_diagonal_errors(h2d, do_sqrt=True):
    nbins = h2d.GetNbinsX()
    hnew = ROOT.TH1D("h_diag" + cu.get_unique_str(), "", nbins, 0, nbins)
    for i in range(1, nbins+1):
        err = h2d.GetBinContent(i, i)
        if do_sqrt and err > 0:
            err = math.sqrt(err)
        hnew.SetBinError(i, err)
    return hnew


def update_hist_bin_content(h_orig, h_to_be_updated):
    if h_orig.GetNbinsX() != h_to_be_updated.GetNbinsX():
        raise RuntimeError("Need same # x bins, %d vs %s" % (h_orig.GetNbinsX(), h_to_be_updated.GetNbinsX()))
    for i in range(0, h_orig.GetNbinsX()+2):
        h_to_be_updated.SetBinContent(i, h_orig.GetBinContent(i))
        # h_to_be_updated.SetBinError(i, h_orig.GetBinError(i))


def update_hist_bin_error(h_orig, h_to_be_updated):
    if h_orig.GetNbinsX() != h_to_be_updated.GetNbinsX():
        raise RuntimeError("Need same # x bins, %d vs %s" % (h_orig.GetNbinsX(), h_to_be_updated.GetNbinsX()))
    for i in range(0, h_orig.GetNbinsX()+2):
        h_to_be_updated.SetBinError(i, h_orig.GetBinError(i))


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
            if value_orig == 0 and value_proj == 0:
                continue
            rel_diff = abs((value_orig - value_proj)/max(abs(value_orig), abs(value_proj)))
            if rel_diff > 1E-3:
                print("draw_projection_comparison: bin %s has different contents: %f vs %f (rel diff %f)" % (i, value_orig, value_proj, rel_diff))

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


def draw_reco_folded(hist_folded, tau, hist_reco_data, hist_reco_mc, title, xtitle, output_filename, folded_subplot=None, folded_subplot_name=""):
    entries = []
    if hist_reco_mc:
        entries.append(
            Contribution(hist_reco_mc, label="MC (detector, bg-subtracted)",
                         line_color=ROOT.kBlack, line_width=1,
                         marker_color=ROOT.kBlack, marker_size=0,
                         normalise_hist=False)
            )
    if hist_reco_data:
        entries.append(
            Contribution(hist_reco_data, label="Data (detector, bg-subtracted)",
                         line_color=ROOT.kBlue, line_width=1,
                         marker_color=ROOT.kBlue, marker_size=0,
                         normalise_hist=False)
            )
    if hist_folded:
        entries.append(
            Contribution(hist_folded, label="Folded data (#tau = %.3g)" % (tau),
                         line_color=ROOT.kRed, line_width=1,
                         marker_color=ROOT.kRed, marker_size=0,
                         normalise_hist=False,
                         subplot=folded_subplot)
            )
    plot = Plot(entries,
                what='hist',
                title=title,
                xtitle=xtitle,
                ytitle="N",
                subplot_type='ratio' if folded_subplot else None,
                subplot_title='Folded / %s' % (folded_subplot_name),
                subplot_limits=(0.8, 1.2))
    plot.default_canvas_size = (800, 600)
    plot.plot("NOSTACK HIST")
    plot.main_pad.SetLogy(1)
    ymax = max(h.GetMaximum() for h in [hist_reco_data, hist_reco_mc, hist_folded] if h)
    plot.container.SetMaximum(ymax * 200)
    # plot.container.SetMinimum(1E-8)
    plot.legend.SetY1NDC(0.77)
    plot.legend.SetX2NDC(0.85)
    plot.save(output_filename)


def check_entries(entries, message=""):
    """Check that at least 1 Contribution has something in it"""
    has_entries = [c.obj.GetEntries() > 0 for c in entries]
    if not any(has_entries):
        if message:
            print("Skipping 0 entries (%s)" % (message))
        else:
            print("Skipping 0 entries")
        return False
    return True


def plot_uncertainty_shifts(total_hist, stat_hist, syst_shifts, systs, output_filename, title, angle_str):
    """
    Parameters
    ----------
    total_hist : TH1
        Distribution with total uncertainty
    stat_hist : TH1
        Distribution with stat-only uncertainties
    syst_shifts : list[TH1]
        Distributions with 1-sigma shift from systematics
    systs : list[dict]
        Dicts describing each syst
    output_filename : str
        Output plot filename
    title : str
        Title to put on plot
    angle_str : str
        Angle name for x axis

    Returns
    -------
    None
        If all hists are empty
    """
    entries = []
    hists = []
    for h, syst_dict in zip(syst_shifts, systs):
        h_fraction = h.Clone()
        h_fraction.Divide(total_hist)
        hists.append(h_fraction)
        for i in range(1, h_fraction.GetNbinsX()+1):
            h_fraction.SetBinContent(i, abs(h_fraction.GetBinContent(i)))
        c = Contribution(h_fraction,
                         label=syst_dict['label'],
                         line_color=syst_dict['colour'],
                         line_style=syst_dict.get('linestyle', 1),
                         line_width=2,
                         marker_size=0,
                         marker_color=syst_dict['colour'],
                         )
        entries.append(c)

    # Add systematic
    h_syst = stat_hist.Clone()
    h_total = total_hist.Clone()
    for i in range(1, h_syst.GetNbinsX()+1):
        if total_hist.GetBinContent(i) > 0:
            h_syst.SetBinContent(i, stat_hist.GetBinError(i) / total_hist.GetBinContent(i))
            h_total.SetBinContent(i, total_hist.GetBinError(i) / total_hist.GetBinContent(i))
        else:
            h_syst.SetBinContent(i, 0)
            h_total.SetBinContent(i, 0)
        h_syst.SetBinError(i, 0)
        h_total.SetBinError(i, 0)
    c_stat = Contribution(h_syst,
                         label="Stat.",
                         line_color=ROOT.kRed,
                         line_style=3,
                         line_width=3,
                         marker_size=0,
                         marker_color=ROOT.kRed,
                         )
    entries.append(c_stat)
    c_tot = Contribution(h_total,
                         label="Total",
                         line_color=ROOT.kBlack,
                         line_style=1,
                         line_width=3,
                         marker_size=0,
                         marker_color=ROOT.kBlack,
                         )
    entries.append(c_tot)

    if not check_entries(entries, "systematic shifts"):
        return
    plot = Plot(entries,
                what="hist",
                title=title,
                xtitle="Particle-level "+angle_str,
                ytitle="| Fractional shift |")
    plot.legend.SetX1(0.55)
    plot.legend.SetY1(0.68)
    plot.legend.SetX2(0.98)
    plot.legend.SetY2(0.88)
    plot.legend.SetNColumns(2)
    plot.y_padding_max_linear = 1.4
    plot.plot("NOSTACK HIST")
    plot.save(output_filename)

    plot.set_logy()
    plot.y_padding_max_log = 50
    plot._set_automatic_y_limits()
    log_filename, ext = os.path.splitext(output_filename)
    plot.save(log_filename+"_log"+ext)


def plot_bias_hist(bias_hist, title, output_filename):
    entries = [
        Contribution(bias_hist, label="Bias histogram",
                     line_color=ROOT.kBlack, line_width=1,
                     marker_color=ROOT.kBlack, marker_size=0,
                     normalise_hist=False),
    ]
    plot = Plot(entries,
                what='hist',
                title=title,
                xtitle="Generator bin",
                ytitle="Bias")
    plot.default_canvas_size = (800, 600)
    plot.plot("NOSTACK HIST")
    plot.set_logy()
    plot.legend.SetY1NDC(0.77)
    plot.legend.SetX2NDC(0.85)
    plot.save(output_filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",
                        help="Source directory with ROOT files")

    parser.add_argument("--angles",
                        choices=list(qgc.VAR_UNFOLD_DICT.keys()),
                        nargs='+',
                        help="Lambda angles to unfold")

    parser.add_argument("--doSummaryPlot",
                        type=lambda x:bool(distutils.util.strtobool(x)),
                        default=False,
                        help='Do summary plot')

    parser.add_argument("--outputDir",
                        default='',
                        help='Output directory')

    region_group = parser.add_argument_group('Region selection')
    region_group.add_argument("--doDijetCentral",
                              action='store_true',
                              help='Do unfolding for dijet (central) jets')
    region_group.add_argument("--doDijetForward",
                              action='store_true',
                              help='Do unfolding for dijet (forward) jets')
    region_group.add_argument("--doDijetCentralGroomed",
                              action='store_true',
                              help='Do unfolding for groomed dijet (central) jets')
    region_group.add_argument("--doDijetForwardGroomed",
                              action='store_true',
                              help='Do unfolding for groomed dijet (forward) jets')
    region_group.add_argument("--doZPJ",
                              action='store_true',
                              help='Do unfolding for Z+jet jets')
    region_group.add_argument("--doZPJGroomed",
                              action='store_true',
                              help='Do unfolding for groomed Z+jet jets')

    regularization_group = parser.add_argument_group('Regularization options')
    regularization_group.add_argument("--regularize",
                                      choices=['None', 'tau', 'L'],
                                      default='None',
                                      help='Regularization scheme')
    regularization_group.add_argument("--nScan",
                                      type=int,
                                      default=100,
                                      help='Number of scan points for regularization')
    regularization_group.add_argument("--biasFactor",
                                      type=float,
                                      default=0,
                                      help='Bias factor for regularization')

    mc_group = parser.add_argument_group('MC input options')
    mc_group.add_argument("--MCinput",
                          type=lambda x:bool(distutils.util.strtobool(x)),
                          default=False,
                          help='Unfold MC instead of data')

    mc_group.add_argument("--MCsplit",
                          type=lambda x:bool(distutils.util.strtobool(x)),
                          default=False,
                          help='Split MC between response & 1D reco, good for testing procedure')

    syst_group = parser.add_argument_group('Systematics options')
    syst_group.add_argument("--doExperimentalSysts",
                            type=lambda x:bool(distutils.util.strtobool(x)),
                            default=False,
                            help='Do experimental systematics (i.e. those that modify response matrix)')

    syst_group.add_argument("--doModelSysts",
                            type=lambda x:bool(distutils.util.strtobool(x)),
                            default=False,
                            help='Do model systematics (i.e. those that modify input to be unfolded)')

    syst_group.add_argument("--doPDFSysts",
                            type=lambda x:bool(distutils.util.strtobool(x)),
                            default=False,
                            help='Do pdf systematics (may be slow!)')

    syst_group.add_argument("--useAltResponse",
                            type=lambda x:bool(distutils.util.strtobool(x)),
                            default=False,
                            help='Use alternate response matrix to unfold')

    args = parser.parse_args()
    print(args)

    if not any([args.doDijetCentral, args.doDijetForward, args.doDijetCentralGroomed, args.doDijetForwardGroomed, args.doZPJ, args.doZPJGroomed]):
        raise RuntimeError("You need to specify at least one signal region e.g. --doDijetCentral")

    if not args.MCinput and args.doModelSysts:
        raise RuntimeError("You cannot do both model systs and run on data")

    # if args.useAltResponse and args.doExperimentalSysts:
    #     args.doExperimentalSysts = False
    #     print("You cannot use both --useAltResponse and --doExperimentalSysts: disabling doExperimentalSysts")

    # # TODO handle both?
    # if args.useAltResponse and args.doModelSysts:
    #     args.doExperimentalSysts = False
    #     print("You cannot use both --useAltResponse and --doModelSysts: disabling doModelSysts")

    # Setup files and regions to unfold
    # --------------------------------------------------------------------------
    src_dir = args.source
    src_dir_systs = os.path.join(src_dir, "systematics_files")

    regions = []

    if any([args.doDijetCentral, args.doDijetForward, args.doDijetCentralGroomed, args.doDijetForwardGroomed]):
        # FOR DIJET:
        input_mc_qcd_mgpythia_tfile = os.path.join(src_dir, qgc.QCD_FILENAME)
        # input_mc_qcd_pythia_tfile = cu.open_root_file(os.path.join(src_dir, qgc.QCD_PYTHIA_ONLY_FILENAME))
        input_mc_qcd_herwig_tfile = os.path.join(src_dir, qgc.QCD_HERWIG_FILENAME)
        # input_mc_qcd_herwig_tfile = os.path.join(src_dir, "uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD_PtReweight.root")

        input_jetht_tfile = os.path.join(src_dir, qgc.JETHT_ZB_FILENAME)

        dijet_region_dict_template = {
            "name": "Dijet",
            "dirname": "Dijet_QG_Unfold_central_tighter",
            "label": "Dijet",
            "data_tfile": input_jetht_tfile,
            "mc_tfile": input_mc_qcd_mgpythia_tfile,
            "mc_label": "MG+Pythia8",
            # "mc_tfile": input_mc_qcd_herwig_tfile,
            # "mc_label": "Herwig++",
            "alt_mc_tfile": input_mc_qcd_herwig_tfile,
            # "alt_mc_label": "Herwig++ (p_{T} reweight)",
            "alt_mc_label": "Herwig++",
            "tau_limits": {
                'jet_puppiMultiplicity': (1E-11, 1E-9),
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
            "experimental_systematics": [
                {
                    "label": "Neutral hadron up",
                    "tfile": os.path.join(src_dir_systs, 'neutralHadronShiftUp', qgc.QCD_FILENAME),
                    "colour": ROOT.kOrange-3,
                },
                {
                    "label": "Neutral hadron down",
                    "tfile": os.path.join(src_dir_systs, 'neutralHadronShiftDown', qgc.QCD_FILENAME),
                    "colour": ROOT.kOrange-3,
                    "linestyle": 2,
                },
                {
                    "label": "Photon up",
                    "tfile": os.path.join(src_dir_systs, 'photonShiftUp', qgc.QCD_FILENAME),
                    "colour": ROOT.kMagenta-3,
                },
                {
                    "label": "Photon down",
                    "tfile": os.path.join(src_dir_systs, 'photonShiftDown', qgc.QCD_FILENAME),
                    "colour": ROOT.kMagenta-3,
                    "linestyle": 2,
                },
                {
                    "label": "JEC up",
                    "tfile": os.path.join(src_dir_systs, 'jecsmear_directionUp', qgc.QCD_FILENAME),
                    "colour": ROOT.kGreen+2,
                },
                {
                    "label": "JEC down",
                    "tfile": os.path.join(src_dir_systs, 'jecsmear_directionDown', qgc.QCD_FILENAME),
                    "colour": ROOT.kGreen+2,
                    "linestyle": 2,
                },
                {
                    "label": "JER up",
                    "tfile": os.path.join(src_dir_systs, 'jersmear_directionUp', qgc.QCD_FILENAME),
                    "colour": ROOT.kOrange+3,
                },
                {
                    "label": "JER down",
                    "tfile": os.path.join(src_dir_systs, 'jersmear_directionDown', qgc.QCD_FILENAME),
                    "colour": ROOT.kOrange+3,
                    "linestyle": 2,
                },
                {
                    "label": "Pileup up",
                    "tfile": os.path.join(src_dir_systs, 'pileup_directionUp', qgc.QCD_FILENAME),
                    "colour": ROOT.kBlue-4,
                },
                {
                    "label": "Pileup down",
                    "tfile": os.path.join(src_dir_systs, 'pileup_directionDown', qgc.QCD_FILENAME),
                    "colour": ROOT.kBlue-4,
                    "linestyle": 2,
                },
            ],
            "model_systematics": [
                {
                    "label": "muR up, muF nominal",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNominal', qgc.QCD_FILENAME),
                    "colour": ROOT.kAzure,
                },
                {
                    "label": "muR down, muF nominal",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFNominal', qgc.QCD_FILENAME),
                    "colour": ROOT.kAzure+1,
                },
                {
                    "label": "muR nominal, muF up",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFUp', qgc.QCD_FILENAME),
                    "colour": ROOT.kAzure+2,
                },
                {
                    "label": "muR nominal, muF down",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFDown', qgc.QCD_FILENAME),
                    "colour": ROOT.kAzure+3,
                },
                {
                    "label": "muR down, muF down",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFDown', qgc.QCD_FILENAME),
                    "colour": ROOT.kAzure+4,
                },
                {
                    "label": "muR up, muF up",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFUp', qgc.QCD_FILENAME),
                    "colour": ROOT.kAzure+5,
                },
                {
                    "label": "Herwig++",
                    "tfile": input_mc_qcd_herwig_tfile,
                    "colour": ROOT.kOrange-3,
                },
            ],
            "pdf_systematics": [
                {
                    "label": "PDF",  # this is a tempalte entry, used for future
                    "tfile": os.path.join(src_dir_systs, 'PDFvariationsTrue', qgc.QCD_FILENAME),
                    "colour": ROOT.kCyan+2,
                    "variations": range(100),  # list of all the variation #s to be used
                },
            ]
        }

        if args.doDijetCentral:
            dijet_region_central_dict = dijet_region_dict_template.copy()
            dijet_region_central_dict['dirname'] = 'Dijet_QG_Unfold_central_tighter'
            dijet_region_central_dict['label'] = 'Dijet central'
            dijet_region_central_dict['name'] = 'Dijet_central'
            regions.append(dijet_region_central_dict)

        if args.doDijetForward:
            dijet_region_forward_dict = dijet_region_dict_template.copy()
            dijet_region_forward_dict['dirname'] = 'Dijet_QG_Unfold_forward_tighter'
            dijet_region_forward_dict['label'] = 'Dijet forward'
            dijet_region_forward_dict['name'] = 'Dijet_forward'
            regions.append(dijet_region_forward_dict)

        if args.doDijetCentralGroomed:
            dijet_region_central_groomed_dict = dijet_region_dict_template.copy()
            dijet_region_central_groomed_dict['dirname'] = 'Dijet_QG_Unfold_central_tighter_groomed'
            dijet_region_central_groomed_dict['label'] = 'Dijet central'
            dijet_region_central_groomed_dict['name'] = 'Dijet_central_groomed'
            regions.append(dijet_region_central_groomed_dict)

        if args.doDijetForwardGroomed:
            dijet_region_forward_groomed_dict = dijet_region_dict_template.copy()
            dijet_region_forward_groomed_dict['dirname'] = 'Dijet_QG_Unfold_forward_tighter_groomed'
            dijet_region_forward_groomed_dict['label'] = 'Dijet forward'
            dijet_region_forward_groomed_dict['name'] = 'Dijet_forward_groomed'
            regions.append(dijet_region_forward_groomed_dict)

    if any([args.doZPJ, args.doZPJGroomed]):
        # FOR Z+JETS:
        input_mc_dy_mgpythia_tfile = os.path.join(src_dir, qgc.DY_FILENAME)
        input_mc_dy_mgherwig_tfile = os.path.join(src_dir, qgc.DY_MG_HERWIG_FILENAME)
        input_mc_dy_herwig_tfile = os.path.join(src_dir, qgc.DY_HERWIG_FILENAME)

        input_singlemu_tfile = os.path.join(src_dir, qgc.SINGLE_MU_FILENAME)

        zpj_region_dict = {
            "name": "ZPlusJets",
            "dirname": "ZPlusJets_QG_Unfold",
            "label": "Z+jets",
            "data_tfile": input_singlemu_tfile,
            "mc_tfile": input_mc_dy_mgpythia_tfile,
            "mc_label": "MG+Pythia8",
            "alt_mc_tfile": input_mc_dy_mgherwig_tfile,
            "alt_mc_label": "MG+Herwig++",
            "tau_limits": {
                'jet_puppiMultiplicity': (1E-10, 1E-4),
                'jet_pTD': (1E-10, 1E-4),
                'jet_LHA': (1E-7, 1E-4),
                'jet_width': (1E-10, 1E-4),
                'jet_thrust': (1E-10, 1E-4),
                'jet_puppiMultiplicity_charged': (1E-10, 1E-4),
                'jet_pTD_charged': (1E-10, 1E-4),
                'jet_LHA_charged': (1E-10, 1E-4),
                'jet_width_charged': (1E-10, 1E-4),
                'jet_thrust_charged': (1E-10, 1E-4),
            },
            "experimental_systematics": [
                {
                    "label": "Neutral hadron up",
                    "tfile": os.path.join(src_dir_systs, 'neutralHadronShiftUp', qgc.DY_FILENAME),
                    "colour": ROOT.kOrange-3,
                },
                {
                    "label": "Neutral hadron down",
                    "tfile": os.path.join(src_dir_systs, 'neutralHadronShiftDown', qgc.DY_FILENAME),
                    "colour": ROOT.kOrange-3,
                    "linestyle": 2,
                },
                {
                    "label": "Photon up",
                    "tfile": os.path.join(src_dir_systs, 'photonShiftUp', qgc.DY_FILENAME),
                    "colour": ROOT.kMagenta-3,
                },
                {
                    "label": "Photon down",
                    "tfile": os.path.join(src_dir_systs, 'photonShiftDown', qgc.DY_FILENAME),
                    "colour": ROOT.kMagenta-3,
                    "linestyle": 2,
                },
                {
                    "label": "JEC up",
                    "tfile": os.path.join(src_dir_systs, 'jecsmear_directionUp', qgc.DY_FILENAME),
                    "colour": ROOT.kGreen+2,
                },
                {
                    "label": "JEC down",
                    "tfile": os.path.join(src_dir_systs, 'jecsmear_directionDown', qgc.DY_FILENAME),
                    "colour": ROOT.kGreen+2,
                    "linestyle": 2,
                },
                {
                    "label": "JER up",
                    "tfile": os.path.join(src_dir_systs, 'jersmear_directionUp', qgc.DY_FILENAME),
                    "colour": ROOT.kOrange+3,
                },
                {
                    "label": "JER down",
                    "tfile": os.path.join(src_dir_systs, 'jersmear_directionDown', qgc.DY_FILENAME),
                    "colour": ROOT.kOrange+3,
                    "linestyle": 2,
                },
                {
                    "label": "Pileup up",
                    "tfile": os.path.join(src_dir_systs, 'pileup_directionUp', qgc.DY_FILENAME),
                    "colour": ROOT.kBlue-4,
                },
                {
                    "label": "Pileup down",
                    "tfile": os.path.join(src_dir_systs, 'pileup_directionDown', qgc.DY_FILENAME),
                    "colour": ROOT.kBlue-4,
                    "linestyle": 2,
                },
            ],
            "model_systematics": [
                {
                    "label": "muR up, muF nominal",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNominal', qgc.DY_FILENAME),
                    "colour": ROOT.kAzure,
                },
                {
                    "label": "muR down, muF nominal",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFNominal', qgc.DY_FILENAME),
                    "colour": ROOT.kAzure+1,
                },
                {
                    "label": "muR nominal, muF up",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFUp', qgc.DY_FILENAME),
                    "colour": ROOT.kAzure+2,
                },
                {
                    "label": "muR nominal, muF down",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFDown', qgc.DY_FILENAME),
                    "colour": ROOT.kAzure+3,
                },
                {
                    "label": "muR down, muF down",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFDown', qgc.DY_FILENAME),
                    "colour": ROOT.kAzure+4,
                },
                {
                    "label": "muR up, muF up",
                    "tfile": os.path.join(src_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFUp', qgc.DY_FILENAME),
                    "colour": ROOT.kAzure+5,
                },
            ],
            "pdf_systematics": [
                {
                    "label": "PDF",
                    "tfile": os.path.join(src_dir_systs, 'PDFvariationsTrue', qgc.DY_FILENAME),
                    "colour": ROOT.kCyan+2,
                    "variations": range(100),
                },
            ]
        }

        if args.doZPJ:
            regions.append(zpj_region_dict)

        if args.doZPJGroomed:
            zpj_region_groomed_dict = zpj_region_dict.copy()
            zpj_region_groomed_dict['dirname'] = 'ZPlusJets_QG_Unfold_groomed'
            zpj_region_groomed_dict['name'] = 'ZPlusJets_groomed'
            zpj_region_groomed_dict['label'] = 'Z+jets'
            regions.append(zpj_region_groomed_dict)

    # Setup various options
    # --------------------------------------------------------------------------

    REGULARIZE = None
    REGULARIZE = "tau"
    # REGULARIZE = "L"
    REGULARIZE = args.regularize

    # Run with MC input instead of data
    MC_INPUT = args.MCinput
    mc_append = "_MC" if MC_INPUT else ""

    # If True, use part of MC for response matrix, and separate part for unfolding
    # as independent test
    # Should always be false for data
    MC_SPLIT = args.MCsplit if args.MCinput else False
    if MC_INPUT:
        mc_append += "_split" if MC_SPLIT else "_all"

    SUBTRACT_FAKES = True  # this should alwys be True
    sub_append = "_subFakes" if SUBTRACT_FAKES else ""

    append = ""

    # Add in systematics
    if args.doExperimentalSysts:
        append += "_experimentalSyst"

    if args.doModelSysts:
        append += "_modelSyst"

    if args.doPDFSysts:
        append += "_pdfSyst"

    if args.useAltResponse:
        append += "_altResponse"

    bias_str = "%g" % args.biasFactor
    bias_str = bias_str.replace(".", "p")
    output_dir = os.path.join(src_dir, "unfolding_better_regularise%s%s%s_densityModeBinWidth_constraintNone%s_signalRegionOnly_biasFactor%s" % (str(REGULARIZE).capitalize(), mc_append, sub_append, append, bias_str))
    # output_dir = os.path.join(src_dir, "unfolding_better_regularise%s%s%s_densityModeBinWidth_constraintNone%s_signalRegionOnly_biasFactor%s_HerwigNominal" % (str(REGULARIZE).capitalize(), mc_append, sub_append, append, bias_str))
    if args.outputDir:
        output_dir = args.outputDir
    cu.check_dir_exists_create(output_dir)

    jet_algo = "AK4 PF PUPPI"
    if "ak8puppi" in src_dir:
        jet_algo = "AK8 PF PUPPI"

    # Save hists etc to ROOT file for access later
    output_tfile = ROOT.TFile("%s/unfolding_result.root" % (output_dir), "RECREATE")

    angles = [a for a in qgc.COMMON_VARS if a.var in args.angles]
    print(angles)

    print("Running TUnfold version", ROOT.TUnfold.GetTUnfoldVersion())

    LAMBDA_VAR_DICTS = qgc.VAR_UNFOLD_DICT
    if 'target0p5' in src_dir:
        LAMBDA_VAR_DICTS = qgc.VAR_UNFOLD_DICT_TARGET0p5
    elif 'target0p6' in src_dir:
        LAMBDA_VAR_DICTS = qgc.VAR_UNFOLD_DICT_TARGET0p6


    # Do unfolding per signal region
    # --------------------------------------------------------------------------
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

        new_tdir = region['name']
        output_tfile.mkdir(new_tdir)
        region_tdir = output_tfile.Get(new_tdir)
        region_tdir.cd()

        # Do 1D unfolding of pt
        # ----------------------------------------------------------------------
        append = "%s_pt" % (region['name'])  # common str to put on filenames, etc
        # print("*"*80)
        # print("Region/var: %s" % (append))
        # print("*"*80)

        # hist_data_reco = cu.get_from_tfile(region['data_tfile'], "%s/hist_pt_reco_all" % (region['dirname']))
        mc_hname_append = "split" if MC_SPLIT else "all"
        if isinstance(region['mc_tfile'], str):
            region['mc_tfile'] = cu.open_root_file(region['mc_tfile'])
        # hist_mc_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_pt_reco_%s" % (region['dirname'], mc_hname_append))
        # hist_mc_gen = cu.get_from_tfile(region['mc_tfile'], "%s/hist_pt_truth_%s" % (region['dirname'], mc_hname_append))
        hist_mc_gen_pt = cu.get_from_tfile(region['mc_tfile'], "%s/hist_pt_truth_%s" % (region['dirname'], mc_hname_append))
        # hist_mc_gen_reco_map = cu.get_from_tfile(region['mc_tfile'], "%s/tu_pt_GenReco_%s" % (region['dirname'], mc_hname_append))
        # TODO!

        # Remake gen hist with physical bins & save to file
        all_pt_bins_gen = np.concatenate((pt_bin_edges_underflow_gen[:-1], pt_bin_edges_gen))
        hist_mc_gen_pt_physical = ROOT.TH1F("mc_gen_pt", ";p_{T}^{jet} [GeV];N", len(all_pt_bins_gen)-1, array('d', all_pt_bins_gen))
        update_hist_bin_content(h_orig=hist_mc_gen_pt, h_to_be_updated=hist_mc_gen_pt_physical)
        update_hist_bin_error(h_orig=hist_mc_gen_pt, h_to_be_updated=hist_mc_gen_pt_physical)
        region_tdir.WriteTObject(hist_mc_gen_pt_physical, "mc_gen_pt")

        # Do unfolding for each angle
        # ----------------------------------------------------------------------
        for angle in angles:
            angle_prepend = "groomed " if "groomed" in region['name'] else ""
            append = "%s_%s%s" % (region['name'], angle_prepend.lower(), angle.var)  # common str to put on filenames, etc
            angle_str = "%s%s (%s)" % (angle_prepend, angle.name, angle.lambda_str)

            print("*"*80)
            print("Region/var: %s" % (append))
            print("*"*80)

            new_tdir = "%s/%s" % (region['name'], angle.var)
            output_tfile.mkdir(new_tdir)
            this_tdir = output_tfile.Get(new_tdir)
            this_tdir.cd()

            # put plots in subdir, to avoid overcrowding
            this_output_dir = "%s/%s/%s" % (output_dir, region['name'], angle.var)
            cu.check_dir_exists_create(this_output_dir)

            # Setup MyUnfolder object to do unfolding etc
            # -------------------------------------------
            angle_bin_edges_reco = LAMBDA_VAR_DICTS[angle.var]['reco']
            angle_bin_edges_gen = LAMBDA_VAR_DICTS[angle.var]['gen']
            angle_shortname = angle.var.replace("jet_", "")

            if not isinstance(region['data_tfile'], ROOT.TFile):
                region['data_tfile'] = cu.open_root_file(region['data_tfile'])
            hist_data_reco = cu.get_from_tfile(region['data_tfile'], "%s/hist_%s_reco_all" % (region['dirname'], angle_shortname))
            mc_hname_append = "split" if MC_SPLIT else "all"
            hist_mc_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_%s" % (region['dirname'], angle_shortname, mc_hname_append))
            hist_mc_gen = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_truth_%s" % (region['dirname'], angle_shortname, mc_hname_append))
            hist_mc_gen_reco_map = cu.get_from_tfile(region['mc_tfile'], "%s/tu_%s_GenReco_%s" % (region['dirname'], angle_shortname, mc_hname_append))

            # Need to scale if using H++ as input
            # hist_mc_gen_reco_map.Scale(1E8)
            # hist_mc_gen.Scale(1E8)
            # hist_mc_reco.Scale(1E8)

            # Actual distribution to be unfolded
            reco_1d = hist_mc_reco.Clone() if MC_INPUT else hist_data_reco

            hist_fakes_reco_fraction = None

            if SUBTRACT_FAKES:
                # to construct our "fakes" template, we use the ratio as predicted by MC, and apply it to data
                # this way we ensure we don't have -ve values, and avoid any issue with cross sections
                hist_mc_fakes_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_fake_%s" % (region['dirname'], angle_shortname, mc_hname_append))
                # hist_mc_fakes_reco.Scale(1E8)
                hist_fakes_reco = hist_mc_fakes_reco.Clone("hist_%s_fakes" % angle_shortname)
                hist_fakes_reco.Divide(hist_mc_reco)

                # plot fake fraction before multiplting by 'data'
                hist_fakes_reco_fraction = hist_fakes_reco.Clone("hist_fakes_reco_fraction")
                plot_fake_fraction(hist_fakes_reco_fraction,
                                   output_filename="%s/fake_fraction_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                                   title="%s region, %s" % (region['label'], angle_str))

                hist_fakes_reco.Multiply(reco_1d)

                # background-subtracted reco hists, only for plotting purposes, not for TUnfold (that does background subtraction internally)
                reco_1d_bg_subtracted = reco_1d.Clone()
                reco_1d_bg_subtracted.Add(hist_fakes_reco, -1)

                chosen_bin = 15
                print("1D reco input without background subtraction:", reco_1d.GetBinContent(chosen_bin))
                print("1D reco input with background subtraction:", reco_1d_bg_subtracted.GetBinContent(chosen_bin))
                print("1D reco input fakes:", hist_fakes_reco.GetBinContent(chosen_bin))

                hist_data_reco_bg_subtracted = hist_data_reco.Clone(hist_data_reco.GetName() + "_bgrSubtracted")
                hist_data_reco_bg_subtracted.Add(hist_fakes_reco, -1)

                hist_mc_reco_bg_subtracted = hist_mc_reco.Clone(hist_mc_reco.GetName() + "_bgrSubtracted")
                hist_mc_reco_bg_subtracted.Add(hist_mc_fakes_reco, -1)  # should this be hist_fakes_reco? depends on what we want to see...

            # hist_mc_gen_reco_neutralUp_map = cu.get_from_tfile(region['mc_neutralUp_tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
            # hist_mc_gen_reco_neutralDown_map = cu.get_from_tfile(region['mc_neutralDown_tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))

            mc_hname_append = "_split" if MC_SPLIT else ""  # FIXME consistency in unfold hist module!
            hist_data_reco_gen_binning = cu.get_from_tfile(region['data_tfile'], "%s/hist_%s_reco_gen_binning" % (region['dirname'], angle_shortname))
            hist_mc_reco_gen_binning = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_gen_binning%s" % (region['dirname'], angle_shortname, mc_hname_append))

            # hist_data_reco_gen_binning.Scale(1e8)
            # hist_mc_reco_gen_binning.Scale(1e8)

            # Actual distribution to be unfolded, but with gen binning
            reco_1d_gen_binning = hist_mc_reco_gen_binning.Clone() if MC_INPUT else hist_data_reco_gen_binning
            # reco_1d_gen_binning = hist_mc_reco_gen_binning # only for testing that everything is setup OK

            if SUBTRACT_FAKES:
                mc_hname_append = "_split" if MC_SPLIT else ""  # FIXME consistency in unfold hist module!
                hist_mc_fakes_reco_gen_binning = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_fake_gen_binning%s" % (region['dirname'], angle_shortname, mc_hname_append))
                # create template as above, but with gen binning
                hist_fakes_reco_gen_binning = hist_mc_fakes_reco_gen_binning.Clone("hist_%s_fakes_gen_binning" % angle_shortname)
                hist_fakes_reco_gen_binning.Divide(hist_mc_reco_gen_binning)
                hist_fakes_reco_gen_binning.Multiply(reco_1d_gen_binning)

                # background-subtracted reco hists, only for plotting purposes, not for TUnfold (that does background subtraction internally)
                reco_1d_gen_binning_bg_subtracted = reco_1d_gen_binning.Clone()
                reco_1d_gen_binning_bg_subtracted.Add(hist_fakes_reco_gen_binning, -1)

                hist_data_reco_gen_binning_bg_subtracted = hist_data_reco_gen_binning.Clone(hist_data_reco_gen_binning.GetName() + "_bgrSubtracted")
                hist_data_reco_gen_binning_bg_subtracted.Add(hist_fakes_reco_gen_binning, -1)

                hist_mc_reco_gen_binning_bg_subtracted = hist_mc_reco_gen_binning.Clone(hist_mc_reco_gen_binning.GetName() + "_bgrSubtracted")
                hist_mc_reco_gen_binning_bg_subtracted.Add(hist_mc_fakes_reco_gen_binning, -1)  # should this be hist_fakes_reco_gen_binning? depends on what we want to see...

            # Setup unfolder object
            # ---------------------
            variable_name = "%s%s" % (angle_prepend, angle.name)
            unfolder = MyUnfolder(response_map=hist_mc_gen_reco_map,
                                  variable_bin_edges_reco=angle_bin_edges_reco,
                                  variable_bin_edges_gen=angle_bin_edges_gen,
                                  variable_name=variable_name,
                                  pt_bin_edges_reco=pt_bin_edges_reco,
                                  pt_bin_edges_gen=pt_bin_edges_gen,
                                  pt_bin_edges_underflow_reco=pt_bin_edges_underflow_reco,
                                  pt_bin_edges_underflow_gen=pt_bin_edges_underflow_gen,
                                  orientation=ROOT.TUnfold.kHistMapOutputHoriz,
                                  # constraintMode=ROOT.TUnfold.kEConstraintArea,
                                  constraintMode=ROOT.TUnfold.kEConstraintNone,
                                  regMode=ROOT.TUnfold.kRegModeCurvature,
                                  densityFlags=ROOT.TUnfoldDensity.kDensityModeBinWidth, # important as we have varying bin sizes!
                                  axisSteering='*[B]')

            unfolder.save_binning(txt_filename="%s/binning_scheme.txt" % (this_output_dir), print_xml=False)

            unfolder_plotter = MyUnfolderPlotter(unfolder)

            # Set what is to be unfolded
            # ------------------------------------------------------------------
            unfolder.set_input(reco_1d, args.biasFactor)

            # Add systematic errors as different response matrices
            # ------------------------------------------------------------------
            if args.doExperimentalSysts:
                chosen_rsp_bin = (18, 18)
                print("nominal response bin content for", chosen_rsp_bin, hist_mc_gen_reco_map.GetBinContent(*chosen_rsp_bin))
                for syst_ind, syst_dict in enumerate(region['experimental_systematics']):
                    if not isinstance(syst_dict['tfile'], ROOT.TFile):
                        syst_dict['tfile'] = cu.open_root_file(syst_dict['tfile'])
                    map_syst = cu.get_from_tfile(syst_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                    print("Adding systematic:", syst_dict['label'])
                    print("    syst bin", chosen_rsp_bin, map_syst.GetBinContent(*chosen_rsp_bin))
                    unfolder.tunfolder.AddSysError(map_syst, syst_dict['label'], unfolder.orientation, ROOT.TUnfoldDensity.kSysErrModeMatrix)

            # Subtract fakes (treat as background)
            # ------------------------------------------------------------------
            if SUBTRACT_FAKES:
                unfolder.subtract_background(hist_fakes_reco, "fakes")

            # Save important stuff to TFile
            # ------------------------------------------------------------------
            this_tdir.WriteTObject(unfolder.detector_binning)
            this_tdir.WriteTObject(unfolder.generator_binning)
            this_tdir.WriteTObject(reco_1d, "unfold_input")
            this_tdir.WriteTObject(reco_1d_gen_binning, "unfold_input_gen_binning")
            this_tdir.WriteTObject(hist_mc_reco, "mc_reco")
            this_tdir.WriteTObject(hist_mc_reco_gen_binning, "mc_reco_gen_binning")
            this_tdir.WriteTObject(hist_mc_gen, "mc_gen")
            this_tdir.WriteTObject(hist_mc_gen_reco_map, "response_map")  # response map
            # this_tdir.WriteTObject(hist_mc_gen_reco_map.RebinY(2), "response_map_rebin")  # response map

            if SUBTRACT_FAKES:
                this_tdir.WriteTObject(hist_fakes_reco, "fakes")
                this_tdir.WriteTObject(hist_fakes_reco_gen_binning, "fakes_gen_binning")
                this_tdir.WriteTObject(reco_1d_bg_subtracted, "unfold_input_bg_subtracted")
                this_tdir.WriteTObject(reco_1d_gen_binning_bg_subtracted, "unfold_input_gen_binning_bg_subtracted")
                this_tdir.WriteTObject(hist_mc_reco_bg_subtracted, "mc_reco_bg_subtracted")
                this_tdir.WriteTObject(hist_mc_reco_gen_binning_bg_subtracted, "mc_reco_gen_binning_bg_subtracted")

            # Plot bias vector
            title = "%s region, %s" % (region['label'], angle_str)
            unfolder_plotter.plot_bias_hist(output_dir=this_output_dir, append=append, title=title)

            # Do any regularization
            # ---------------------
            unfolder.print_condition_number()
            # tau = 1E-10
            tau = 0
            scan_mode = ROOT.TUnfoldDensity.kEScanTauRhoAvgSys
            scan_distribution = "generatordistribution"
            if REGULARIZE == "L":
                print("Regularizing with ScanLcurve, please be patient...")
                l_scanner = LCurveScanner()
                tau = l_scanner.scan_L(tunfolder=unfolder.tunfolder,
                                       n_scan=args.nScan,
                                       tau_min=region['tau_limits'][angle.var][0],
                                       tau_max=region['tau_limits'][angle.var][1])
                print("Found tau:", tau)
                l_scanner.plot_scan_L_curve(output_filename="%s/scanL_%s.%s" % (this_output_dir, unfolder.variable_name, OUTPUT_FMT))
                l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_%s.%s" % (this_output_dir, unfolder.variable_name, OUTPUT_FMT))

            elif REGULARIZE == "tau":
                print("Regularizing with ScanTau, please be patient...")
                tau_scanner = TauScanner()
                tau = tau_scanner.scan_tau(tunfolder=unfolder.tunfolder,
                                           n_scan=args.nScan,
                                           tau_min=region['tau_limits'][angle.var][0],
                                           tau_max=region['tau_limits'][angle.var][1],
                                           scan_mode=scan_mode,
                                           distribution=scan_distribution,
                                           axis_steering=unfolder.axisSteering)
                print("Found tau:", tau)
                tau_scanner.plot_scan_tau(output_filename="%s/scantau_%s.%s" % (this_output_dir, unfolder.variable_name, OUTPUT_FMT))
                this_tdir.WriteTObject(tau_scanner.graph_all_scan_points, "regularize_all_scan_points")
                this_tdir.WriteTObject(tau_scanner.graph_best_scan_point, "regularize_best_scan_point")

            # Do unfolding!
            # ---------------------
            unfolder.do_unfolding(tau)
            unfolded_1d = unfolder.get_output()
            unfolded_1d.SetName("unfolded_1d")
            chosen_bin = 18
            print("Bin %d:" % chosen_bin, unfolded_1d.GetBinContent(chosen_bin))
            print("original uncert:", unfolded_1d.GetBinError(chosen_bin))

            # Get various error matrices
            # ------------------------------------------------------------------
            # stat errors only - do before or after systematics?
            ematrix_input = unfolder.get_ematrix_input() # stat errors from input to be unfolded
            print("ematrix_input uncert:", ematrix_input.GetBinContent(chosen_bin, chosen_bin))
            this_tdir.WriteTObject(ematrix_input, "ematrix_input")

            ematrix_sys_uncorr = unfolder.get_ematrix_sys_uncorr() # stat errors in response matrix
            print("ematrix_sys_uncorr uncert:", ematrix_sys_uncorr.GetBinContent(chosen_bin, chosen_bin))
            this_tdir.WriteTObject(ematrix_sys_uncorr, "ematrix_sys_uncorr")

            ematrix_stat_sum = ematrix_input.Clone("ematrix_stat_sum")
            ematrix_stat_sum.Add(ematrix_sys_uncorr) # total 'stat' errors

            error_stat_1d = make_hist_from_diagonal_errors(ematrix_stat_sum, do_sqrt=True) # note that bin contents = 0, only bin errors are non-0
            this_tdir.WriteTObject(error_stat_1d, "error_stat_1d")
            print("stat uncert:", error_stat_1d.GetBinError(chosen_bin))

            ematrix_total = unfolder.get_ematrix_total()
            error_total_1d = make_hist_from_diagonal_errors(ematrix_total, do_sqrt=True) # note that bin contents = 0, only bin errors are non-0
            this_tdir.WriteTObject(ematrix_total, "ematrix_total_1d")
            this_tdir.WriteTObject(error_total_1d, "error_total_1d")
            print("total uncert:", error_total_1d.GetBinError(chosen_bin))

            # Update errors to big unfolded 1D
            update_hist_bin_error(h_orig=error_total_1d, h_to_be_updated=unfolded_1d)
            print("new uncert:", unfolded_1d.GetBinError(chosen_bin))
            this_tdir.WriteTObject(unfolded_1d)

            # Get shifts due to systematics
            # ------------------------------------------------------------------
            systematic_shift_hists = []
            if args.doExperimentalSysts:
                for syst_dict in region['experimental_systematics']:
                    syst_label = syst_dict['label']
                    h_syst = unfolder.tunfolder.GetDeltaSysSource(syst_label,
                                                                  '%s_%s_shift_%s' % (region['name'], angle.var, syst_label),
                                                                  "",
                                                                  "generator",  # must be the same as what's used in get_output
                                                                  unfolder.axisSteering,
                                                                  unfolder.use_axis_binning)
                    systematic_shift_hists.append(h_syst)
                    this_tdir.WriteTObject(h_syst)

            # Draw unified unfolded distributions
            # ------------------------------------------------------------------
            # unfolded, gen, and reco for comparison
            plot_simple_unfolded(unfolded=unfolded_1d,
                                 tau=tau,
                                 reco=None,
                                 gen=hist_mc_gen,
                                 fake=None,
                                 output_filename="%s/unfolded_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                                 title="%s region, %s" % (region['label'], angle_str))

            # reco using detector binning
            plot_simple_detector(reco_data=reco_1d,
                                 reco_mc=hist_mc_reco,
                                 reco_mc_fake=hist_mc_fakes_reco if SUBTRACT_FAKES else None,
                                 reco_data_fake=hist_fakes_reco if SUBTRACT_FAKES else None,
                                 output_filename="%s/detector_reco_binning_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                                 title="%s region, %s" % (region['label'], angle_str))

            if SUBTRACT_FAKES:
                # same plot but with background-subtracted reco
                plot_simple_detector(reco_data=reco_1d_bg_subtracted,
                                     reco_mc=hist_mc_reco_bg_subtracted,
                                     reco_mc_fake=hist_mc_fakes_reco,
                                     reco_data_fake=hist_fakes_reco,
                                     output_filename="%s/detector_reco_binning_bg_subtracted_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                                     title="%s region, %s" % (region['label'], angle_str))

            # reco using gen binning
            plot_simple_detector(reco_data=reco_1d_gen_binning,
                                 reco_mc=hist_mc_reco_gen_binning,
                                 reco_mc_fake=hist_mc_fakes_reco_gen_binning if SUBTRACT_FAKES else None,
                                 reco_data_fake=hist_fakes_reco_gen_binning if SUBTRACT_FAKES else None,
                                 output_filename="%s/detector_gen_binning_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                                 title="%s region, %s" % (region['label'], angle_str))

            if SUBTRACT_FAKES:
                # same but with background-subtracted
                plot_simple_detector(reco_data=reco_1d_gen_binning_bg_subtracted,
                                     reco_mc=hist_mc_reco_gen_binning_bg_subtracted,
                                     reco_mc_fake=hist_mc_fakes_reco_gen_binning,
                                     reco_data_fake=hist_fakes_reco_gen_binning,
                                     output_filename="%s/detector_gen_binning_bg_subtracted_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                                     title="%s region, %s" % (region['label'], angle_str))

            # Draw collapsed distributions
            # ---------------------
            # TODO!

            # Draw projections of response matrix vs 1D hist to check normalisation OK
            # Only makes sense if the same MC events go into matrix & 1D plot
            # ------------------------------------------------------------------
            if not MC_SPLIT:
                proj_reco = hist_mc_gen_reco_map.ProjectionY("proj_reco_%s" % (append))

                proj_gen = hist_mc_gen_reco_map.ProjectionX("proj_gen_%s" % (append))
                draw_projection_comparison(hist_mc_gen, proj_gen,
                                           title="%s\n%s region" % (jet_algo, region['label']),
                                           xtitle="%s, Generator binning" % (angle_str),
                                           output_filename="%s/projection_gen_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

                print("projection reco #bins:", proj_reco.GetNbinsX())
                print("response map # bins x:", hist_mc_gen_reco_map.GetNbinsX())
                print("response map # bins y:", hist_mc_gen_reco_map.GetNbinsY())
                if SUBTRACT_FAKES:
                    print("reco bg subtracted #bins:", hist_mc_reco_bg_subtracted.GetNbinsX())
                    print(proj_reco.GetNbinsX())
                    # Do the same but with backgrounds subtracted from the 1D
                    draw_projection_comparison(hist_mc_reco_bg_subtracted, proj_reco,
                                               title="%s\n%s region" % (jet_algo, region['label']),
                                               xtitle="%s, Detector binning" % (angle_str),
                                               output_filename="%s/projection_reco_bg_subtracted_%s.%s" % (this_output_dir, append, OUTPUT_FMT))
                else:
                    draw_projection_comparison(hist_mc_reco, proj_reco,
                           title="%s\n%s region" % (jet_algo, region['label']),
                           xtitle="%s, Detector binning" % (angle_str),
                           output_filename="%s/projection_reco_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                           print_bin_comparison=False)

            # Draw matrices
            # ------------------------------------------------------------------
            plot_args = dict(output_dir=this_output_dir, append=append)
            
            title = "Response matrix, %s region, %s" % (region['label'], angle_str)
            unfolder_plotter.draw_response_matrix(title=title, **plot_args)

            title = "Probability matrix, %s region, %s" % (region['label'], angle_str)
            unfolder_plotter.draw_probability_matrix(title=title, **plot_args)

            title = "Correlation matrix, %s region, %s" % (region['label'], angle_str)
            unfolder_plotter.draw_correlation_matrix(title=title, **plot_args)

            title = "Error matrix (input), %s region, %s" % (region['label'], angle_str)
            unfolder_plotter.draw_error_matrix_input(title=title, **plot_args)

            title = "Error matrix (sys uncorr), %s region, %s" % (region['label'], angle_str)
            unfolder_plotter.draw_error_matrix_sys_uncorr(title=title, **plot_args)

            title = "Error matrix (total), %s region, %s" % (region['label'], angle_str)
            unfolder_plotter.draw_error_matrix_total(title=title, **plot_args)

            # Do forward-folding to check unfolding
            # ------------------------------------------------------------------
            # Do it on the unfolded result
            hist_data_folded_unfolded = unfolder.get_folded_output(hist_name="folded_unfolded_%s" % (append))
            this_tdir.WriteTObject(hist_data_folded_unfolded, "folded_unfolded_1d")
            draw_reco_folded(hist_folded=hist_data_folded_unfolded,
                             tau=tau,
                             hist_reco_data=reco_1d_bg_subtracted if SUBTRACT_FAKES else reco_1d,
                             hist_reco_mc=hist_mc_reco_bg_subtracted if SUBTRACT_FAKES else hist_mc_reco,
                             folded_subplot=reco_1d_bg_subtracted if SUBTRACT_FAKES else reco_1d,
                             folded_subplot_name='Data',
                             title="%s\n%s region" % (jet_algo, region['label']),
                             xtitle="%s, Detector binning" % (angle_str),
                             output_filename="%s/folded_unfolded_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

            hist_gen_folded = None
            if MC_INPUT:
                hist_gen_folded = unfolder.get_folded_hist(hist_mc_gen)
                this_tdir.WriteTObject(hist_gen_folded, "folded_gen_1d")
                draw_reco_folded(hist_folded=hist_gen_folded,
                                 tau=0,
                                 hist_reco_data=None,
                                 hist_reco_mc=hist_mc_reco_bg_subtracted if SUBTRACT_FAKES else hist_mc_reco,
                                 folded_subplot=hist_mc_reco_bg_subtracted if SUBTRACT_FAKES else hist_mc_reco,
                                 folded_subplot_name='MC',
                                 title="%s\n%s region" % (jet_algo, region['label']),
                                 xtitle="%s, Detector binning" % (angle_str),
                                 output_filename="%s/folded_gen_%s.%s" % (this_output_dir, append, OUTPUT_FMT))


            # ------------------------------------------------------------------
            # Do unfolding using alternate response matrix
            # ------------------------------------------------------------------
            alt_unfolded_1d = None
            alt_hist_mc_gen = None
            alt_tau = 0
            if args.useAltResponse:
                print("*** Unfolding with alternate response matrix ***")
                if not isinstance(region['alt_mc_tfile'], ROOT.TFile):
                    region['alt_mc_tfile'] = cu.open_root_file(region['alt_mc_tfile'])
                alt_hist_mc_gen = cu.get_from_tfile(region['alt_mc_tfile'], "%s/hist_%s_truth_all" % (region['dirname'], angle_shortname))
                hist_mc_gen_reco_map_alt = cu.get_from_tfile(region['alt_mc_tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                alt_unfolder = MyUnfolder(response_map=hist_mc_gen_reco_map_alt,
                                          variable_bin_edges_reco=unfolder.variable_bin_edges_reco,
                                          variable_bin_edges_gen=unfolder.variable_bin_edges_gen,
                                          variable_name=unfolder.variable_name,
                                          pt_bin_edges_reco=unfolder.pt_bin_edges_reco,
                                          pt_bin_edges_gen=unfolder.pt_bin_edges_gen,
                                          pt_bin_edges_underflow_reco=unfolder.pt_bin_edges_underflow_reco,
                                          pt_bin_edges_underflow_gen=unfolder.pt_bin_edges_underflow_gen,
                                          orientation=unfolder.orientation,
                                          constraintMode=unfolder.constraintMode,
                                          regMode=unfolder.regMode,
                                          densityFlags=unfolder.densityFlags,
                                          axisSteering=unfolder.axisSteering)

                # Set what is to be unfolded
                # --------------------------------------------------------------
                alt_unfolder.set_input(reco_1d, args.biasFactor)

                # Subtract fakes (treat as background)
                # --------------------------------------------------------------
                if SUBTRACT_FAKES:
                    alt_unfolder.subtract_background(hist_fakes_reco, "fakes")

                # Save important stuff to TFile
                # --------------------------------------------------------------
                this_tdir.WriteTObject(hist_mc_gen_reco_map_alt, "alt_response_map")  # response map

                # Do any regularization
                # --------------------------------------------------------------
                alt_tau = 0
                if REGULARIZE == "L":
                    print("Regularizing alternative with ScanL, please be patient...")
                    alt_L_scanner = LCurveScanner()
                    alt_tau = l_scanner.scan_L(tunfolder=unfolder.tunfolder,
                                               n_scan=args.nScan,
                                               tau_min=region['tau_limits'][angle.var][0],
                                               tau_max=region['tau_limits'][angle.var][1])
                    print("Found tau:", alt_tau)
                    l_scanner.plot_scan_L_curve(output_filename="%s/scanL_alt_%s.%s" % (this_output_dir, unfolder.variable_name, OUTPUT_FMT))
                    l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_alt_%s.%s" % (this_output_dir, unfolder.variable_name, OUTPUT_FMT))

                elif REGULARIZE == "tau":
                    print("Regularizing alternative with ScanTau, please be patient...")
                    alt_tau_scanner = TauScanner()
                    alt_tau = alt_tau_scanner.scan_tau(tunfolder=alt_unfolder.tunfolder,
                                                       n_scan=args.nScan,
                                                       tau_min=region['tau_limits'][angle.var][0],
                                                       tau_max=region['tau_limits'][angle.var][1],
                                                       scan_mode=scan_mode,
                                                       distribution=scan_distribution,
                                                       axis_steering=alt_unfolder.axisSteering)
                    print("Found tau for alt matrix:", alt_tau)
                    alt_tau_scanner.plot_scan_tau(output_filename="%s/scantau_alt_%s.%s" % (this_output_dir, alt_unfolder.variable_name, OUTPUT_FMT))

                # Do unfolding!
                # --------------------------------------------------------------
                alt_unfolder.do_unfolding(alt_tau)
                alt_unfolded_1d = alt_unfolder.get_output()
                alt_unfolded_1d.SetName("alt_unfolded_1d")
                print("Bin %d:" % chosen_bin, alt_unfolded_1d.GetBinContent(chosen_bin))
                print("original uncert:", alt_unfolded_1d.GetBinError(chosen_bin))

                # Get error matrix to update errors
                # --------------------------------------------------------------
                alt_ematrix_total = alt_unfolder.get_ematrix_total()
                alt_error_total_1d = make_hist_from_diagonal_errors(ematrix_total, do_sqrt=True) # note that bin contents = 0, only bin errors are non-0
                this_tdir.WriteTObject(alt_ematrix_total, "alt_ematrix_total_1d")
                this_tdir.WriteTObject(alt_error_total_1d, "alt_error_total_1d")
                print("total uncert:", alt_error_total_1d.GetBinError(chosen_bin))

                # Update errors to big unfolded 1D
                update_hist_bin_error(h_orig=alt_error_total_1d, h_to_be_updated=alt_unfolded_1d)
                print("new uncert:", alt_unfolded_1d.GetBinError(chosen_bin))
                this_tdir.WriteTObject(alt_unfolded_1d)

            # ------------------------------------------------------------------
            # MODEL INPUT VARIATIONS
            # ------------------------------------------------------------------
            # For each model variation, we unfold using the same settings as
            # the nominal one, just changing the input 1D hist
            if args.doModelSysts:
                for ind, syst_dict in enumerate(region['model_systematics']):
                    syst_label = syst_dict['label']
                    syst_label_no_spaces = syst_dict['label'].replace(", ", "_").replace(" ", "_")

                    print("*** Unfolding with alternate input:", syst_label, "(%d/%d) ***" % (ind+1, len(region['model_systematics'])))
                    is_herwig = "Herwig" in syst_label

                    mc_hname_append = "split" if MC_SPLIT else "all"
                    if is_herwig:
                        # use all the stats!
                        mc_hname_append = "all"
                    if not isinstance(syst_dict['tfile'], ROOT.TFile):
                        syst_dict['tfile'] = cu.open_root_file(syst_dict['tfile'])
                    hist_syst_reco = cu.get_from_tfile(syst_dict['tfile'], "%s/hist_%s_reco_%s" % (region['dirname'], angle_shortname, mc_hname_append))
                    hist_syst_gen = cu.get_from_tfile(syst_dict['tfile'], "%s/hist_%s_truth_%s" % (region['dirname'], angle_shortname, mc_hname_append))


                    syst_unfolder = MyUnfolder(response_map=unfolder.response_map,
                                               variable_bin_edges_reco=unfolder.variable_bin_edges_reco,
                                               variable_bin_edges_gen=unfolder.variable_bin_edges_gen,
                                               variable_name=unfolder.variable_name,
                                               pt_bin_edges_reco=unfolder.pt_bin_edges_reco,
                                               pt_bin_edges_gen=unfolder.pt_bin_edges_gen,
                                               pt_bin_edges_underflow_reco=unfolder.pt_bin_edges_underflow_reco,
                                               pt_bin_edges_underflow_gen=unfolder.pt_bin_edges_underflow_gen,
                                               orientation=unfolder.orientation,
                                               constraintMode=unfolder.constraintMode,
                                               regMode=unfolder.regMode,
                                               densityFlags=unfolder.densityFlags,
                                               axisSteering=unfolder.axisSteering)

                    if is_herwig:
                        herwig_sf = 1E8  # should've done this earlier
                        hist_syst_reco.Scale(herwig_sf)
                        hist_syst_gen.Scale(herwig_sf)
                        # SetEpsMatrix ensures rank properly calculated when inverting
                        # "rank of matrix E 55 expect 170"
                        syst_unfolder.tunfolder.SetEpsMatrix(1E-18)

                    # Set what is to be unfolded
                    # --------------------------------------------------------------
                    syst_unfolder.set_input(hist_syst_reco, args.biasFactor)

                    # Subtract fakes (treat as background)
                    # --------------------------------------------------------------
                    if SUBTRACT_FAKES:
                        # Use the background template from the nominal MC
                        # (since we're only testing different input shapes,
                        # and our bkg estimate is always from MC)
                        hist_fakes_syst = hist_fakes_reco_fraction.Clone("hist_fakes_syst_%s" % syst_label_no_spaces)
                        hist_fakes_syst.Multiply(hist_syst_reco)
                        syst_unfolder.subtract_background(hist_fakes_syst, "fakes")

                    plot_simple_detector(reco_data=hist_syst_reco,
                                         reco_mc=hist_mc_reco,
                                         reco_mc_fake=None,
                                         reco_data_fake=None,
                                         output_filename="%s/detector_reco_binning_bg_subtracted_model_%s_%s.%s" % (this_output_dir, syst_label_no_spaces, append, OUTPUT_FMT),
                                         title="%s region, %s" % (region['label'], angle_str))

                    # Add systematic errors as different response matrices
                    # ----------------------------------------------------
                    if args.doExperimentalSysts:
                        chosen_rsp_bin = (18, 18)
                        print("nominal response bin content for", chosen_rsp_bin, syst_unfolder.response_map.GetBinContent(*chosen_rsp_bin))
                        for exp_dict in region['experimental_systematics']:
                            if not isinstance(exp_dict['tfile'], ROOT.TFile):
                                exp_dict['tfile'] = cu.open_root_file(exp_dict['tfile'])
                            map_syst = cu.get_from_tfile(exp_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                            print("Adding systematic:", exp_dict['label'])
                            print("    syst bin", chosen_rsp_bin, map_syst.GetBinContent(*chosen_rsp_bin))
                            syst_unfolder.tunfolder.AddSysError(map_syst, exp_dict['label'], syst_unfolder.orientation, ROOT.TUnfoldDensity.kSysErrModeMatrix)

                    # Save important stuff to TFile
                    # --------------------------------------------------------------
                    this_tdir.WriteTObject(hist_syst_reco, "syst_%s_unfold_input" % (syst_label_no_spaces))

                    # Do any regularization
                    # --------------------------------------------------------------
                    syst_output_dir = "%s/%s/%s" % (output_dir, region['name'], angle.var)
                    syst_tau = 0
                    if REGULARIZE == "L":
                        print("Regularizing systematic model with ScanL, please be patient...")
                        syst_l_scanner = LCurveScanner()
                        syst_tau = syst_l_scanner.scan_L(tunfolder=unfolder.tunfolder,
                                                         n_scan=args.nScan,
                                                         tau_min=region['tau_limits'][angle.var][0],
                                                         tau_max=region['tau_limits'][angle.var][1])
                        print("Found tau:", syst_tau)
                        syst_l_scanner.plot_scan_L_curve(output_filename="%s/scanL_syst_%s_%s.%s" % (this_output_dir, syst_label_no_spaces, unfolder.variable_name, OUTPUT_FMT))
                        syst_l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_syst_%s_%s.%s" % (this_output_dir, syst_label_no_spaces, unfolder.variable_name, OUTPUT_FMT))
                    elif REGULARIZE == "tau":
                        print("Regularizing systematic model with ScanTau, please be patient...")
                        syst_tau_scanner = TauScanner()
                        syst_tau = syst_tau_scanner.scan_tau(tunfolder=syst_unfolder.tunfolder,
                                                             n_scan=args.nScan,
                                                             tau_min=region['tau_limits'][angle.var][0],
                                                             tau_max=region['tau_limits'][angle.var][1],
                                                             scan_mode=scan_mode,
                                                             distribution=scan_distribution,
                                                             axis_steering=syst_unfolder.axisSteering)
                        print("Found tau for syst matrix:", syst_tau)
                        syst_tau_scanner.plot_scan_tau(output_filename="%s/scantau_syst_%s_%s.%s" % (this_output_dir, syst_label_no_spaces, syst_unfolder.variable_name, OUTPUT_FMT))

                    region['model_systematics'][ind]['tau'] = syst_tau

                    # Do unfolding!
                    # --------------------------------------------------------------
                    syst_unfolder.do_unfolding(syst_tau)
                    syst_unfolded_1d = syst_unfolder.get_output()
                    syst_unfolded_1d.SetName("syst_%s_unfolded_1d" % (syst_label_no_spaces))
                    print("Bin %d:" % (chosen_bin), syst_unfolded_1d.GetBinContent(chosen_bin))
                    print("original uncert:", syst_unfolded_1d.GetBinError(chosen_bin))

                    # Get error matrix to update errors
                    # --------------------------------------------------------------
                    syst_ematrix_total = syst_unfolder.get_ematrix_total()
                    syst_error_total_1d = make_hist_from_diagonal_errors(ematrix_total, do_sqrt=True) # note that bin contents = 0, only bin errors are non-0
                    this_tdir.WriteTObject(syst_ematrix_total, "syst_%s_ematrix_total_1d" % (syst_label_no_spaces))
                    this_tdir.WriteTObject(syst_error_total_1d, "syst_%s_error_total_1d" % (syst_label_no_spaces))
                    print("total uncert:", syst_error_total_1d.GetBinError(chosen_bin))

                    # Update errors to big unfolded 1D
                    update_hist_bin_error(h_orig=syst_error_total_1d, h_to_be_updated=syst_unfolded_1d)
                    print("new uncert:", syst_unfolded_1d.GetBinError(chosen_bin))
                    this_tdir.WriteTObject(syst_unfolded_1d)

                    region['model_systematics'][ind]['unfolded_1d'] = syst_unfolded_1d
                    region['model_systematics'][ind]['gen_1d'] = hist_syst_gen

            # ------------------------------------------------------------------
            # DO PDF VARIATIONS
            # ------------------------------------------------------------------
            if args.doPDFSysts:
                # first construct all new systematic variations
                original_pdf_dict = region['pdf_systematics'][0]
                original_pdf_dict['label'] = '_PDF_template'  # initial _ to ignore it later on
                tfile = original_pdf_dict['tfile']
                if not isinstance(tfile, ROOT.TFile):
                    tfile = cu.open_root_file(tfile)

                region['pdf_systematics']  = []
                for pdf_ind in original_pdf_dict['variations']:
                    mc_hname_append = "split" if MC_SPLIT else "all"
                    mc_hname_append = "all"
                    region['pdf_systematics'].append(
                        {
                            "label": "PDF_%d" % (pdf_ind),
                            "hist_reco": cu.get_from_tfile(tfile, "%s/hist_%s_reco_%s_PDF_%d" % (region['dirname'], angle_shortname, mc_hname_append, pdf_ind)),
                            "hist_gen": cu.get_from_tfile(tfile, "%s/hist_%s_gen_%s_PDF_%d" % (region['dirname'], angle_shortname, mc_hname_append, pdf_ind)),
                            "colour": ROOT.kCyan+2,
                        })


                # Now run over all variations like for model systs
                for ind, pdf_dict in enumerate(region['pdf_systematics']):
                    pdf_label = pdf_dict['label']
                    syst_label_no_spaces = pdf_label.replace(" ", "_")

                    if pdf_label.startswith("_"):
                        continue

                    print("*** Unfolding with alternate input:", pdf_label, "(%d/%d) ***" % (ind+1, len(region['pdf_systematics'])))

                    hist_syst_reco = pdf_dict['hist_reco']
                    hist_syst_gen = pdf_dict['hist_gen']

                    syst_unfolder = MyUnfolder(response_map=unfolder.response_map,
                                               variable_bin_edges_reco=unfolder.variable_bin_edges_reco,
                                               variable_bin_edges_gen=unfolder.variable_bin_edges_gen,
                                               variable_name=unfolder.variable_name,
                                               pt_bin_edges_reco=unfolder.pt_bin_edges_reco,
                                               pt_bin_edges_gen=unfolder.pt_bin_edges_gen,
                                               pt_bin_edges_underflow_reco=unfolder.pt_bin_edges_underflow_reco,
                                               pt_bin_edges_underflow_gen=unfolder.pt_bin_edges_underflow_gen,
                                               orientation=unfolder.orientation,
                                               constraintMode=unfolder.constraintMode,
                                               regMode=unfolder.regMode,
                                               densityFlags=unfolder.densityFlags,
                                               axisSteering=unfolder.axisSteering)

                    # Set what is to be unfolded
                    # --------------------------------------------------------------
                    syst_unfolder.set_input(hist_syst_reco, args.biasFactor)

                    # Subtract fakes (treat as background)
                    # --------------------------------------------------------------
                    if SUBTRACT_FAKES:
                        # Use the background template from the nominal MC
                        # (since we're only testing different input shapes,
                        # and our bkg estimate is always from MC)
                        hist_fakes_syst = hist_fakes_reco_fraction.Clone("hist_fakes_syst_%s" % syst_label_no_spaces)
                        hist_fakes_syst.Multiply(hist_syst_reco)
                        syst_unfolder.subtract_background(hist_fakes_syst, "fakes")

                    plot_simple_detector(reco_data=hist_syst_reco,
                                         reco_mc=hist_mc_reco,
                                         reco_mc_fake=None,
                                         reco_data_fake=None,
                                         output_filename="%s/detector_reco_binning_bg_subtracted_model_%s_%s.%s" % (this_output_dir, syst_label_no_spaces, append, OUTPUT_FMT),
                                         title="%s region, %s" % (region['label'], angle_str))

                    # Add systematic errors as different response matrices
                    # ----------------------------------------------------
                    if args.doExperimentalSysts:
                        chosen_rsp_bin = (18, 18)
                        print("nominal response bin content for", chosen_rsp_bin, syst_unfolder.response_map.GetBinContent(*chosen_rsp_bin))
                        for exp_dict in region['experimental_systematics']:
                            if not isinstance(exp_dict['tfile'], ROOT.TFile):
                                exp_dict['tfile'] = cu.open_root_file(exp_dict['tfile'])
                            map_syst = cu.get_from_tfile(exp_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                            print("Adding systematic:", exp_dict['label'])
                            print("    syst bin", chosen_rsp_bin, map_syst.GetBinContent(*chosen_rsp_bin))
                            syst_unfolder.tunfolder.AddSysError(map_syst, exp_dict['label'], syst_unfolder.orientation, ROOT.TUnfoldDensity.kSysErrModeMatrix)

                    # Save important stuff to TFile
                    # --------------------------------------------------------------
                    this_tdir.WriteTObject(hist_syst_reco, "syst_%s_unfold_input" % (syst_label_no_spaces))

                    # Do any regularization
                    # --------------------------------------------------------------
                    syst_output_dir = "%s/%s/%s" % (output_dir, region['name'], angle.var)
                    syst_tau = 0
                    if REGULARIZE == "L":
                        print("Regularizing systematic model with ScanL, please be patient...")
                        syst_l_scanner = LCurveScanner()
                        syst_tau = syst_l_scanner.scan_L(tunfolder=unfolder.tunfolder,
                                                         n_scan=args.nScan,
                                                         tau_min=region['tau_limits'][angle.var][0],
                                                         tau_max=region['tau_limits'][angle.var][1])
                        print("Found tau:", syst_tau)
                        syst_l_scanner.plot_scan_L_curve(output_filename="%s/scanL_syst_%s_%s.%s" % (this_output_dir, syst_label_no_spaces, unfolder.variable_name, OUTPUT_FMT))
                        syst_l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_syst_%s_%s.%s" % (this_output_dir, syst_label_no_spaces, unfolder.variable_name, OUTPUT_FMT))
                    elif REGULARIZE == "tau":
                        print("Regularizing systematic model with ScanTau, please be patient...")
                        syst_tau_scanner = TauScanner()
                        syst_tau = syst_tau_scanner.scan_tau(tunfolder=syst_unfolder.tunfolder,
                                                             n_scan=args.nScan,
                                                             tau_min=region['tau_limits'][angle.var][0],
                                                             tau_max=region['tau_limits'][angle.var][1],
                                                             scan_mode=scan_mode,
                                                             distribution=scan_distribution,
                                                             axis_steering=syst_unfolder.axisSteering)
                        print("Found tau for syst matrix:", syst_tau)
                        syst_tau_scanner.plot_scan_tau(output_filename="%s/scantau_syst_%s_%s.%s" % (this_output_dir, syst_label_no_spaces, syst_unfolder.variable_name, OUTPUT_FMT))

                    region['pdf_systematics'][ind]['tau'] = syst_tau

                    # Do unfolding!
                    # --------------------------------------------------------------
                    syst_unfolder.do_unfolding(syst_tau)
                    syst_unfolded_1d = syst_unfolder.get_output()
                    syst_unfolded_1d.SetName("syst_%s_unfolded_1d" % (syst_label_no_spaces))
                    print("Bin %d:" % (chosen_bin), syst_unfolded_1d.GetBinContent(chosen_bin))
                    print("original uncert:", syst_unfolded_1d.GetBinError(chosen_bin))

                    # Get error matrix to update errors
                    # --------------------------------------------------------------
                    syst_ematrix_total = syst_unfolder.get_ematrix_total()
                    syst_error_total_1d = make_hist_from_diagonal_errors(ematrix_total, do_sqrt=True) # note that bin contents = 0, only bin errors are non-0
                    this_tdir.WriteTObject(syst_ematrix_total, "syst_%s_ematrix_total_1d" % (syst_label_no_spaces))
                    this_tdir.WriteTObject(syst_error_total_1d, "syst_%s_error_total_1d" % (syst_label_no_spaces))
                    print("total uncert:", syst_error_total_1d.GetBinError(chosen_bin))

                    # Update errors to big unfolded 1D
                    update_hist_bin_error(h_orig=syst_error_total_1d, h_to_be_updated=syst_unfolded_1d)
                    print("new uncert:", syst_unfolded_1d.GetBinError(chosen_bin))
                    this_tdir.WriteTObject(syst_unfolded_1d)

                    region['pdf_systematics'][ind]['unfolded_1d'] = syst_unfolded_1d
                    region['pdf_systematics'][ind]['gen_1d'] = hist_syst_gen

            # ------------------------------------------------------------------
            # PLOTTING LOTS OF THINGS
            # ------------------------------------------------------------------

            # Some common plotting vars
            # ------------------------------------------------------------------
            detector_title = "Detector-level " + angle_str
            particle_title = "Particle-level " + angle_str
            normalised_differential_label = "#frac{1}{#sigma} #frac{d^{2}#sigma}{dp_{T} d%s}" % angle.lambda_str
            summary_1d_entries = []  # for final summary plot

            # Draw individual pt bin plots - GEN binning
            # ------------------------------------------------------------------
            for ibin_pt in range(0, len(pt_bin_edges_gen)-1):

                this_pt_bin_tdir = this_tdir.mkdir("gen_bin_%d" % (ibin_pt))
                print("Individual gen bin", ibin_pt, "=", pt_bin_edges_gen[ibin_pt], pt_bin_edges_gen[ibin_pt+1])

                # Produce 1D hists for this pt bin
                # --------------------------------------------------------------
                # Unfolded hists
                mc_gen_hist_bin = unfolder.get_var_hist_pt_binned(hist_mc_gen, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(mc_gen_hist_bin, "mc_gen_hist_bin")

                unfolded_hist_bin = unfolder.get_var_hist_pt_binned(unfolded_1d, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(unfolded_hist_bin, "unfolded_hist_bin")

                unfolded_hist_bin_stat_errors = unfolder.get_var_hist_pt_binned(error_stat_1d, ibin_pt, binning_scheme="generator")
                update_hist_bin_content(h_orig=unfolded_hist_bin, h_to_be_updated=unfolded_hist_bin_stat_errors)  # use stat errors, update central values
                this_pt_bin_tdir.WriteTObject(unfolded_hist_bin_stat_errors, "unfolded_hist_bin_stat_errors")

                unfolded_hist_bin_total_errors = unfolder.get_var_hist_pt_binned(error_total_1d, ibin_pt, binning_scheme="generator")
                update_hist_bin_content(h_orig=unfolded_hist_bin, h_to_be_updated=unfolded_hist_bin_total_errors)  # use total errors, update central values
                this_pt_bin_tdir.WriteTObject(unfolded_hist_bin_total_errors, "unfolded_hist_bin_total_errors") # same as unfolded_1d?

                # Reco level but with gen binning
                # For 'data'
                reco_hist_bin_gen_binning = unfolder.get_var_hist_pt_binned(reco_1d_gen_binning, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(reco_hist_bin_gen_binning, "reco_hist_bin_gen_binning")

                if SUBTRACT_FAKES:
                    reco_hist_bg_subtracted_bin_gen_binning = unfolder.get_var_hist_pt_binned(reco_1d_gen_binning_bg_subtracted, ibin_pt, binning_scheme="generator")
                    this_pt_bin_tdir.WriteTObject(reco_hist_bg_subtracted_bin_gen_binning, "reco_hist_bg_subtracted_bin_gen_binning")

                # For MC
                mc_reco_hist_bin_gen_binning = unfolder.get_var_hist_pt_binned(hist_mc_reco_gen_binning, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(mc_reco_hist_bin_gen_binning, "mc_reco_hist_bin_gen_binning")

                if SUBTRACT_FAKES:
                    mc_reco_hist_bg_subtracted_bin_gen_binning = unfolder.get_var_hist_pt_binned(hist_mc_reco_gen_binning_bg_subtracted, ibin_pt, binning_scheme="generator")
                    this_pt_bin_tdir.WriteTObject(mc_reco_hist_bg_subtracted_bin_gen_binning, "mc_reco_hist_bg_subtracted_bin_gen_binning")


                # Make lots of plots
                # ------------------------------------------------------------
                lw = 2
                # common hist settings
                title = "%s\n%s region\n%g < p_{T}^{jet} < %g GeV" % (jet_algo, region['label'], pt_bin_edges_gen[ibin_pt], pt_bin_edges_gen[ibin_pt+1])
                if "ptavebinning" in src_dir.lower():
                    title = "%s\n%s region\n%g < #LT p_{T}^{jet} #GT < %g GeV" % (jet_algo, region['label'], pt_bin_edges_gen[ibin_pt], pt_bin_edges_gen[ibin_pt+1])
                common_hist_args = dict(
                    what="hist",
                    title=title,
                    subplot_type='ratio',
                    # subplot_limits=(0.5, 1.5),
                    subplot_limits=(0.75, 1.25),
                )
                subplot_title = "Unfolded / Gen"

                # PLOT UNFOLDED DATA
                # --------------------------------------------------------------
                gen_colour = ROOT.kRed
                unfolded_basic_colour = ROOT.kAzure+7
                unfolded_stat_colour = ROOT.kGreen+2
                unfolded_total_colour = ROOT.kBlack

                def _modify_plot(this_plot):
                    this_plot.legend.SetX1(0.6)
                    this_plot.legend.SetY1(0.68)
                    this_plot.legend.SetX2(0.98)
                    this_plot.legend.SetY2(0.9)
                    this_plot.left_margin = 0.13

                # unnormalised version
                entries = [
                    Contribution(mc_gen_hist_bin,
                                 label="Generator",
                                 line_color=gen_colour, line_width=lw,
                                 marker_color=gen_colour, marker_size=0,
                                 normalise_hist=False),
                    Contribution(unfolded_hist_bin_total_errors,
                                 label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                 line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_gen_hist_bin,
                                 normalise_hist=False),
                    Contribution(unfolded_hist_bin_stat_errors,
                                 label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                 line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_stat_colour, marker_size=0,
                                 subplot=mc_gen_hist_bin,
                                 normalise_hist=False),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=particle_title,
                            ytitle="N",
                            subplot_title='Unfolded / gen',
                            **common_hist_args)
                _modify_plot(plot)
                plot.plot("NOSTACK E1")
                plot.save("%s/unfolded_unnormalised_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # now normalise each plot to unity
                # Note that this modifies e.g. mc_gen_hist_bin, so from this point
                # onwards it will be normalised to unity
                entries = [
                    Contribution(mc_gen_hist_bin,
                                 label="Generator (MG+Pythia8)",
                                 line_color=gen_colour, line_width=lw,
                                 marker_color=gen_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(unfolded_hist_bin_total_errors,
                                 label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                 line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_gen_hist_bin,
                                 normalise_hist=True),
                    Contribution(unfolded_hist_bin_stat_errors,
                                 label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                 line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_stat_colour, marker_size=0,
                                 normalise_hist=True),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=particle_title,
                            ytitle=normalised_differential_label,
                            subplot_title=subplot_title,
                            **common_hist_args)
                _modify_plot(plot)
                plot.plot("NOSTACK E1")
                # plot.save("%s/unfolded_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Do a version where divided by bin width
                # Note that these hists are already normalised to 1!
                # Do not use normalise_hist!
                mc_gen_hist_bin_div_bin_width = qgp.hist_divide_bin_width(mc_gen_hist_bin)
                unfolded_hist_bin_stat_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_stat_errors)
                unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_total_errors)
                entries = [
                    Contribution(mc_gen_hist_bin_div_bin_width,
                                 label="Generator (MG+Pythia8)",
                                 line_color=gen_colour, line_width=lw,
                                 marker_color=gen_colour, marker_size=0,
                                 normalise_hist=False),
                    Contribution(unfolded_hist_bin_total_errors_div_bin_width,
                                 label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                 line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_gen_hist_bin_div_bin_width,
                                 normalise_hist=False),
                    Contribution(unfolded_hist_bin_stat_errors_div_bin_width,
                                 label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                 line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_stat_colour, marker_size=0,
                                 normalise_hist=False),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=particle_title,
                            ytitle=normalised_differential_label,
                            subplot_title=subplot_title,
                            **common_hist_args)
                _modify_plot(plot)
                plot.plot("NOSTACK E1")
                plot.save("%s/unfolded_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                summary_1d_entries.append([
                    (mc_gen_hist_bin_div_bin_width,
                        dict(label=region['mc_label'],
                             line_color=gen_colour, line_width=lw,
                             marker_color=gen_colour, marker_size=0,
                             normalise_hist=True)),  # generator
                    (unfolded_hist_bin_total_errors_div_bin_width,
                        dict(label="Unfolded\n($\\tau = %.3g$)" % (tau),
                             line_color=unfolded_total_colour, line_width=lw, line_style=3,
                             marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=True)),  # unfolded data
                ])

                # Unfolded plots with alternate response matrix results as well
                # --------------------------------------------------------------
                # And alternate MC gen level to compare
                # (i.e. is the difference between unfolded results < diff at gen level?)
                if args.useAltResponse:
                    alt_mc_gen_hist_bin = alt_unfolder.get_var_hist_pt_binned(alt_hist_mc_gen, ibin_pt, binning_scheme="generator")
                    alt_unfolded_hist_bin_total_errors = alt_unfolder.get_var_hist_pt_binned(alt_unfolded_1d, ibin_pt, binning_scheme="generator")
                    this_pt_bin_tdir.WriteTObject(alt_unfolded_hist_bin_total_errors, "alt_unfolded_hist_bin_total_errors")

                    alt_gen_colour = ROOT.kBlack
                    alt_colour = ROOT.kBlue

                    entries = [
                        Contribution(mc_gen_hist_bin,
                                     label="Generator (%s)" % (region['mc_label']),
                                     line_color=gen_colour, line_width=lw,
                                     marker_color=gen_colour, marker_size=0,
                                     normalise_hist=True),
                        Contribution(alt_mc_gen_hist_bin,
                                     label="Generator (%s)" % (region['alt_mc_label']),
                                     line_color=alt_gen_colour, line_width=lw,
                                     marker_color=alt_gen_colour, marker_size=0,
                                     subplot=mc_gen_hist_bin,
                                     normalise_hist=True),
                        Contribution(unfolded_hist_bin_stat_errors,
                                     label="Unfolded (#tau = %.3g) (stat err)\n(%s response matrix)" % (tau, region['mc_label']),
                                     line_color=unfolded_stat_colour, line_width=lw, line_style=2,
                                     marker_color=unfolded_stat_colour, marker_size=0,
                                     subplot=mc_gen_hist_bin,
                                     normalise_hist=True),
                        Contribution(alt_unfolded_hist_bin_total_errors,
                                     label="Unfolded (#tau = %.3g) (stat err)\n(%s response matrix)" % (alt_tau, region['alt_mc_label']),
                                     line_color=alt_colour, line_width=lw, line_style=3,
                                     marker_color=alt_colour, marker_size=0,
                                     subplot=mc_gen_hist_bin,
                                     normalise_hist=True),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=particle_title,
                                ytitle=normalised_differential_label,
                                subplot_title='#splitline{Unfolded / Gen}{(%s)}' % (region['mc_label']),
                                **common_hist_args)
                    plot.legend.SetX1(0.55)
                    plot.legend.SetY1(0.72)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.88)
                    plot.plot("NOSTACK E1")
                    # plot.save("%s/unfolded_%s_alt_response_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                    # Do a version where divided by bin width
                    # Note that inputs are already normalised to 1
                    # Do not use normalise_hist!
                    alt_mc_gen_hist_bin_div_bin_width = qgp.hist_divide_bin_width(alt_mc_gen_hist_bin)
                    alt_unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(alt_unfolded_hist_bin_total_errors)
                    unfolded_hist_bin_stat_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_stat_errors)
                    # unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_total_errors)

                    entries = [
                        Contribution(mc_gen_hist_bin_div_bin_width,
                                     label="Generator (%s)" % (region['mc_label']),
                                     line_color=gen_colour, line_width=lw,
                                     marker_color=gen_colour, marker_size=0,
                                     normalise_hist=False),
                        Contribution(alt_mc_gen_hist_bin_div_bin_width,
                                     label="Generator (%s)" % (region['alt_mc_label']),
                                     line_color=alt_gen_colour, line_width=lw,
                                     marker_color=alt_gen_colour, marker_size=0,
                                     subplot=mc_gen_hist_bin_div_bin_width,
                                     normalise_hist=False),
                        Contribution(unfolded_hist_bin_stat_errors_div_bin_width,
                                     label="Unfolded (#tau = %.3g) (stat err)\n(%s response matrix)" % (tau, region['mc_label']),
                                     line_color=unfolded_stat_colour, line_width=lw, line_style=2,
                                     marker_color=unfolded_stat_colour, marker_size=0,
                                     subplot=mc_gen_hist_bin_div_bin_width,
                                     normalise_hist=False),
                        Contribution(alt_unfolded_hist_bin_total_errors_div_bin_width,
                                     label="Unfolded (#tau = %.3g) (stat err)\n(%s response matrix)" % (alt_tau, region['alt_mc_label']),
                                     line_color=alt_colour, line_width=lw, line_style=3,
                                     marker_color=alt_colour, marker_size=0,
                                     subplot=mc_gen_hist_bin_div_bin_width,
                                     normalise_hist=False),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=particle_title,
                                ytitle=normalised_differential_label,
                                subplot_title='#splitline{Unfolded / Gen}{(%s)}' % (region['mc_label']),
                                **common_hist_args)
                    plot.legend.SetX1(0.55)
                    plot.legend.SetY1(0.72)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.88)
                    plot.plot("NOSTACK E1")
                    plot.save("%s/unfolded_%s_alt_response_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Unfolded plots with variations in input model systematics plotted
                # --------------------------------------------------------------
                if args.doModelSysts:
                    syst_entries = []
                    syst_entries_div_bin_width = []
                    for syst_dict in region['model_systematics']:
                        syst_label = syst_dict['label']
                        syst_label_no_spaces = syst_dict['label'].replace(" ", "_")
                        syst_tau = syst_dict['tau']
                        syst_unfolded_1d = syst_dict.get('unfolded_1d', None)
                        if not syst_unfolded_1d:
                            continue
                        syst_unfolded_hist_bin_total_errors = unfolder.get_var_hist_pt_binned(syst_unfolded_1d, ibin_pt, binning_scheme="generator")
                        this_pt_bin_tdir.WriteTObject(syst_unfolded_hist_bin_total_errors, "syst_%s_unfolded_hist_bin_total_errors" % (syst_label_no_spaces))

                        syst_gen_1d = syst_dict.get('gen_1d', None)
                        if not syst_gen_1d:
                            continue
                        syst_gen_1d_bin = unfolder.get_var_hist_pt_binned(syst_gen_1d, ibin_pt, binning_scheme="generator")
                        this_pt_bin_tdir.WriteTObject(syst_gen_1d_bin, "syst_%s_gen_hist_bin" % (syst_label_no_spaces))

                        syst_entries.extend([
                            Contribution(syst_unfolded_hist_bin_total_errors,
                                         label="Unfolded (#tau = %.3g) (total err) (%s)" % (syst_tau, syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=1,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         subplot=syst_gen_1d_bin,
                                         normalise_hist=True),
                            Contribution(syst_gen_1d_bin,
                                         label="Generator (%s)" % (syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=2,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         normalise_hist=True),
                        ])
                        # already normalised to 1
                        # do not use normalise_hist!
                        syst_unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(syst_unfolded_hist_bin_total_errors)
                        syst_gen_1d_bin_div_bin_width = qgp.hist_divide_bin_width(syst_gen_1d_bin)
                        syst_entries_div_bin_width.extend([
                            Contribution(syst_unfolded_hist_bin_total_errors_div_bin_width,
                                         label="Unfolded (#tau = %.3g) (total err) (%s)" % (syst_tau, syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=1,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         subplot=syst_gen_1d_bin_div_bin_width,
                                         normalise_hist=False),
                            Contribution(syst_gen_1d_bin_div_bin_width,
                                         label="Generator (%s)" % (syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=2,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         normalise_hist=False),
                        ])

                    if len(syst_entries):
                        entries = [
                            Contribution(mc_gen_hist_bin,
                                         label="Generator (%s)" % (region['mc_label']),
                                         line_color=gen_colour, line_width=lw,
                                         marker_color=gen_colour, marker_size=0,
                                         normalise_hist=True),
                            Contribution(unfolded_hist_bin_total_errors,
                                         label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                         line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                         marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                         subplot=mc_gen_hist_bin,
                                         normalise_hist=True),
                            Contribution(unfolded_hist_bin_stat_errors,
                                         label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                         line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                         marker_color=unfolded_stat_colour, marker_size=0,
                                         normalise_hist=True),
                        ]
                        entries.extend(syst_entries)
                        if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                            continue
                        plot = Plot(entries,
                                    xtitle=particle_title,
                                    ytitle=normalised_differential_label,
                                    subplot_title=subplot_title,
                                    **common_hist_args)
                        plot.legend.SetX1(0.55)
                        plot.legend.SetY1(0.72)
                        plot.legend.SetX2(0.98)
                        plot.legend.SetY2(0.88)
                        plot.legend.SetNColumns(3)
                        plot.plot("NOSTACK E1")
                        # plot.save("%s/unfolded_%s_syst_model_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                        # Do a version where divided by bin width
                        mc_gen_hist_bin_div_bin_width = qgp.hist_divide_bin_width(mc_gen_hist_bin)
                        unfolded_hist_bin_stat_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_stat_errors)
                        unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_total_errors)
                        entries_div_bin_width = [
                            Contribution(mc_gen_hist_bin_div_bin_width,
                                         label="Generator (%s)" % (region['mc_label']),
                                         line_color=gen_colour, line_width=lw,
                                         marker_color=gen_colour, marker_size=0,
                                         normalise_hist=False),
                            Contribution(unfolded_hist_bin_total_errors_div_bin_width,
                                         label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                         line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                         marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                         subplot=mc_gen_hist_bin_div_bin_width,
                                         normalise_hist=False),
                            Contribution(unfolded_hist_bin_stat_errors_div_bin_width,
                                         label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                         line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                         marker_color=unfolded_stat_colour, marker_size=0,
                                         normalise_hist=False),
                        ]
                        entries_div_bin_width.extend(syst_entries_div_bin_width)
                        if not check_entries(entries_div_bin_width, "%s %d" % (append, ibin_pt)):
                            continue
                        plot = Plot(entries_div_bin_width,
                                    xtitle=particle_title,
                                    ytitle=normalised_differential_label,
                                    subplot_title=subplot_title,
                                    **common_hist_args)
                        plot.legend.SetX1(0.55)
                        plot.legend.SetY1(0.72)
                        plot.legend.SetX2(0.98)
                        plot.legend.SetY2(0.88)
                        plot.legend.SetNColumns(3)
                        plot.plot("NOSTACK E1")
                        plot.save("%s/unfolded_%s_syst_model_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Unfolded plots with PDF variations
                # --------------------------------------------------------------
                if args.doPDFSysts:
                    syst_entries = []
                    syst_entries_div_bin_width = []
                    for syst_dict in region['pdf_systematics']:
                        syst_label = syst_dict['label']
                        if syst_label.startswith("_"):
                            continue
                        syst_label_no_spaces = syst_dict['label'].replace(" ", "_")
                        syst_tau = syst_dict['tau']
                        syst_unfolded_1d = syst_dict.get('unfolded_1d', None)
                        if not syst_unfolded_1d:
                            continue
                        syst_unfolded_hist_bin_total_errors = unfolder.get_var_hist_pt_binned(syst_unfolded_1d, ibin_pt, binning_scheme="generator")
                        this_pt_bin_tdir.WriteTObject(syst_unfolded_hist_bin_total_errors, "syst_%s_unfolded_hist_bin_total_errors" % (syst_label_no_spaces))

                        syst_gen_1d = syst_dict.get('gen_1d', None)
                        if not syst_gen_1d:
                            continue
                        syst_gen_1d_bin = unfolder.get_var_hist_pt_binned(syst_gen_1d, ibin_pt, binning_scheme="generator")
                        this_pt_bin_tdir.WriteTObject(syst_gen_1d_bin, "syst_%s_gen_hist_bin" % (syst_label_no_spaces))

                        syst_entries.extend([
                            Contribution(syst_unfolded_hist_bin_total_errors,
                                         label="Unfolded (#tau = %.3g) (total err) (%s)" % (syst_tau, syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=1,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         subplot=syst_gen_1d_bin,
                                         normalise_hist=True),
                            Contribution(syst_gen_1d_bin,
                                         label="Generator (%s)" % (syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=2,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         normalise_hist=True),
                        ])
                        # already normalised to 1
                        # do not use normalise_hist!
                        syst_unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(syst_unfolded_hist_bin_total_errors)
                        syst_gen_1d_bin_div_bin_width = qgp.hist_divide_bin_width(syst_gen_1d_bin)
                        syst_entries_div_bin_width.extend([
                            Contribution(syst_unfolded_hist_bin_total_errors_div_bin_width,
                                         label="Unfolded (#tau = %.3g) (total err) (%s)" % (syst_tau, syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=1,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         subplot=syst_gen_1d_bin_div_bin_width,
                                         normalise_hist=False),
                            Contribution(syst_gen_1d_bin_div_bin_width,
                                         label="Generator (%s)" % (syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=2,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         normalise_hist=False),
                        ])

                    if len(syst_entries):
                        entries = [
                            Contribution(mc_gen_hist_bin,
                                         label="Generator (%s)" % (region['mc_label']),
                                         line_color=gen_colour, line_width=lw,
                                         marker_color=gen_colour, marker_size=0,
                                         normalise_hist=True),
                            Contribution(unfolded_hist_bin_total_errors,
                                         label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                         line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                         marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                         subplot=mc_gen_hist_bin,
                                         normalise_hist=True),
                            Contribution(unfolded_hist_bin_stat_errors,
                                         label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                         line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                         marker_color=unfolded_stat_colour, marker_size=0,
                                         normalise_hist=True),
                        ]
                        entries.extend(syst_entries)
                        if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                            continue
                        plot = Plot(entries,
                                    xtitle=particle_title,
                                    ytitle=normalised_differential_label,
                                    subplot_title=subplot_title,
                                    **common_hist_args)
                        plot.legend.SetX1(0.55)
                        plot.legend.SetY1(0.72)
                        plot.legend.SetX2(0.98)
                        plot.legend.SetY2(0.88)
                        plot.legend.SetNColumns(3)
                        plot.plot("NOSTACK E1")
                        # plot.save("%s/unfolded_%s_syst_model_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                        # Do a version where divided by bin width
                        mc_gen_hist_bin_div_bin_width = qgp.hist_divide_bin_width(mc_gen_hist_bin)
                        unfolded_hist_bin_stat_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_stat_errors)
                        unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_total_errors)
                        entries_div_bin_width = [
                            Contribution(mc_gen_hist_bin_div_bin_width,
                                         label="Generator (%s)" % (region['mc_label']),
                                         line_color=gen_colour, line_width=lw,
                                         marker_color=gen_colour, marker_size=0,
                                         normalise_hist=False),
                            Contribution(unfolded_hist_bin_total_errors_div_bin_width,
                                         label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                         line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                         marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                         subplot=mc_gen_hist_bin_div_bin_width,
                                         normalise_hist=False),
                            Contribution(unfolded_hist_bin_stat_errors_div_bin_width,
                                         label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                         line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                         marker_color=unfolded_stat_colour, marker_size=0,
                                         normalise_hist=False),
                        ]
                        entries_div_bin_width.extend(syst_entries_div_bin_width)
                        if not check_entries(entries_div_bin_width, "%s %d" % (append, ibin_pt)):
                            continue
                        plot = Plot(entries_div_bin_width,
                                    xtitle=particle_title,
                                    ytitle=normalised_differential_label,
                                    subplot_title=subplot_title,
                                    **common_hist_args)
                        plot.legend.SetX1(0.55)
                        plot.legend.SetY1(0.72)
                        plot.legend.SetX2(0.98)
                        plot.legend.SetY2(0.88)
                        plot.legend.SetNColumns(3)
                        plot.plot("NOSTACK E1")
                        plot.save("%s/unfolded_%s_pdf_model_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))


                # --------------------------------------------------------------
                # PLOT UNCERTAINTY SHIFTS
                # --------------------------------------------------------------
                systematic_shift_hists_bin = [unfolder.get_var_hist_pt_binned(h, ibin_pt, binning_scheme='generator')
                                              for h in systematic_shift_hists]
                unfolded_stat_error_bin = unfolder.get_var_hist_pt_binned(error_stat_1d, ibin_pt, binning_scheme="generator")
                unfolded_total_error_bin =  unfolder.get_var_hist_pt_binned(unfolded_1d, ibin_pt, binning_scheme="generator")
                plot_uncertainty_shifts(total_hist=unfolded_total_error_bin,
                                        stat_hist=unfolded_stat_error_bin,
                                        syst_shifts=systematic_shift_hists_bin,
                                        systs=region['experimental_systematics'],
                                        output_filename='%s/unfolded_systs_%s_bin_%d.%s' % (this_output_dir, append, ibin_pt, OUTPUT_FMT),
                                        title=title,
                                        angle_str=angle_str)

                # --------------------------------------------------------------
                # PLOT RECO-LEVEL DISTRIBUTIONS (with gen binning)
                # --------------------------------------------------------------
                # Reco level, generator-binning
                reco_mc_colour = ROOT.kGreen+2
                reco_data_colour = ROOT.kRed
                entries = [
                    Contribution(mc_reco_hist_bin_gen_binning,
                                 label="MC",
                                 line_color=reco_mc_colour, line_width=lw,
                                 marker_color=reco_mc_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(reco_hist_bin_gen_binning,
                                 label="Data",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_reco_hist_bin_gen_binning,
                                 normalise_hist=True),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Data / MC',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                # plot.save("%s/detector_gen_binning_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Again but divided by bin width
                mc_reco_hist_bin_gen_binning_div_bin_width = qgp.hist_divide_bin_width(mc_reco_hist_bin_gen_binning)
                reco_hist_bin_gen_binning_div_bin_width = qgp.hist_divide_bin_width(reco_hist_bin_gen_binning)
                entries = [
                    Contribution(mc_reco_hist_bin_gen_binning_div_bin_width,
                                 label="MC",
                                 line_color=reco_mc_colour, line_width=lw,
                                 marker_color=reco_mc_colour, marker_size=0,
                                 normalise_hist=False),
                    Contribution(reco_hist_bin_gen_binning_div_bin_width,
                                 label="Data",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_reco_hist_bin_gen_binning_div_bin_width,
                                 normalise_hist=False),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Data / MC',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/detector_gen_binning_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                if SUBTRACT_FAKES:
                    # Same but background-subtracted
                    entries = [
                        Contribution(mc_reco_hist_bg_subtracted_bin_gen_binning,
                                     label="MC (bg-subtracted)",
                                     line_color=reco_mc_colour, line_width=lw,
                                     marker_color=reco_mc_colour, marker_size=0,
                                     normalise_hist=True),
                        Contribution(reco_hist_bg_subtracted_bin_gen_binning,
                                     label="Data (bg-subtracted)",
                                     line_color=reco_data_colour, line_width=lw,
                                     marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                     subplot=mc_reco_hist_bg_subtracted_bin_gen_binning,
                                     normalise_hist=True),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=detector_title,
                                ytitle=normalised_differential_label,
                                subplot_title='Data / MC',
                                **common_hist_args)
                    plot.legend.SetX1(0.6)
                    plot.legend.SetY1(0.75)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.9)
                    plot.plot("NOSTACK E1")
                    # plot.save("%s/detector_gen_binning_bg_subtracted_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                    mc_reco_hist_bg_subtracted_bin_gen_binning_div_bin_width = qgp.hist_divide_bin_width(mc_reco_hist_bg_subtracted_bin_gen_binning)
                    reco_hist_bg_subtracted_bin_gen_binning_div_bin_width = qgp.hist_divide_bin_width(reco_hist_bg_subtracted_bin_gen_binning)
                    entries = [
                        Contribution(mc_reco_hist_bg_subtracted_bin_gen_binning_div_bin_width,
                                     label="MC (bg-subtracted)",
                                     line_color=reco_mc_colour, line_width=lw,
                                     marker_color=reco_mc_colour, marker_size=0,
                                     normalise_hist=False),
                        Contribution(reco_hist_bg_subtracted_bin_gen_binning_div_bin_width,
                                     label="Data (bg-subtracted)",
                                     line_color=reco_data_colour, line_width=lw,
                                     marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                     subplot=mc_reco_hist_bg_subtracted_bin_gen_binning_div_bin_width,
                                     normalise_hist=False),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=detector_title,
                                ytitle=normalised_differential_label,
                                subplot_title='Data / MC',
                                **common_hist_args)
                    plot.legend.SetX1(0.6)
                    plot.legend.SetY1(0.75)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.9)
                    plot.plot("NOSTACK E1")
                    plot.save("%s/detector_gen_binning_bg_subtracted_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))


            # Draw individual pt bin plots - RECO binning
            # ------------------------------------------------------------------
            for ibin_pt in range(0, len(pt_bin_edges_reco)-1):
                this_pt_bin_tdir = this_tdir.mkdir("reco_bin_%d" % (ibin_pt))
                print("Individual detector bin", ibin_pt, "=", pt_bin_edges_reco[ibin_pt], pt_bin_edges_reco[ibin_pt+1])

                # Get 1D hists
                # The folded unfolded result
                folded_unfolded_hist_bin_reco_binning = unfolder.get_var_hist_pt_binned(hist_data_folded_unfolded, ibin_pt, binning_scheme="detector")
                this_pt_bin_tdir.WriteTObject(folded_unfolded_hist_bin_reco_binning, "folded_unfolded_hist_bin_reco_binning")

                # here this is the thing to be unfolded, should be data or MC depending on MC_INPUT flag
                reco_hist_bin_reco_binning = unfolder.get_var_hist_pt_binned(reco_1d, ibin_pt, binning_scheme="detector")
                this_pt_bin_tdir.WriteTObject(reco_hist_bin_reco_binning, "reco_hist_bin_reco_binning")

                if SUBTRACT_FAKES:
                    reco_hist_bg_subtracted_bin_reco_binning = unfolder.get_var_hist_pt_binned(reco_1d_bg_subtracted, ibin_pt, binning_scheme="detector")
                    this_pt_bin_tdir.WriteTObject(reco_hist_bg_subtracted_bin_reco_binning, "reco_hist_bg_subtracted_bin_reco_binning")

                # Get MC hists
                mc_reco_hist_bin_reco_binning = unfolder.get_var_hist_pt_binned(hist_mc_reco, ibin_pt, binning_scheme="detector")
                this_pt_bin_tdir.WriteTObject(mc_reco_hist_bin_reco_binning, "mc_reco_hist_bin_reco_binning")

                if SUBTRACT_FAKES:
                    mc_reco_hist_bg_subtracted_bin_reco_binning = unfolder.get_var_hist_pt_binned(hist_mc_reco_bg_subtracted, ibin_pt, binning_scheme="detector")
                    this_pt_bin_tdir.WriteTObject(mc_reco_hist_bg_subtracted_bin_reco_binning, "mc_reco_hist_bg_subtracted_bin_reco_binning")

                # Do the plots
                # --------------------------------------------------------------
                # common hist settings
                lw = 2
                title = "%s\n%s region\n%g < p_{T}^{Jet} < %g GeV" % (jet_algo, region['label'], pt_bin_edges_reco[ibin_pt], pt_bin_edges_reco[ibin_pt+1])
                common_hist_args = dict(
                    what="hist",
                    title=title,
                    subplot_type='ratio',
                    # subplot_limits=(0.5, 1.5),
                    subplot_limits=(0.75, 1.25),
                )

                # Reco only, detector-binning
                entries = [
                    Contribution(mc_reco_hist_bin_reco_binning,
                                 label="MC",
                                 line_color=reco_mc_colour, line_width=lw,
                                 marker_color=reco_mc_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(reco_hist_bin_reco_binning,
                                 label="Data",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_reco_hist_bin_reco_binning,
                                 normalise_hist=True),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Data / MC',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                # plot.save("%s/detector_reco_binning_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # same but divided by bin width
                mc_reco_hist_bin_reco_binning_div_bin_width = qgp.hist_divide_bin_width(mc_reco_hist_bin_reco_binning)
                reco_hist_bin_reco_binning_div_bin_width = qgp.hist_divide_bin_width(reco_hist_bin_reco_binning)
                entries = [
                    Contribution(mc_reco_hist_bin_reco_binning_div_bin_width,
                                 label="MC",
                                 line_color=reco_mc_colour, line_width=lw,
                                 marker_color=reco_mc_colour, marker_size=0,
                                 normalise_hist=False),
                    Contribution(reco_hist_bin_reco_binning_div_bin_width,
                                 label="Data",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_reco_hist_bin_reco_binning_div_bin_width,
                                 normalise_hist=False),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Data / MC',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/detector_reco_binning_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                if SUBTRACT_FAKES:
                    # Same but background-subtracted
                    entries = [
                        Contribution(mc_reco_hist_bg_subtracted_bin_reco_binning,
                                     label="MC (bg-subtracted)",
                                     line_color=reco_mc_colour, line_width=lw,
                                     marker_color=reco_mc_colour, marker_size=0,
                                     normalise_hist=True),
                        Contribution(reco_hist_bg_subtracted_bin_reco_binning,
                                     label="Data (bg-subtracted)",
                                     line_color=reco_data_colour, line_width=lw,
                                     marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                     subplot=mc_reco_hist_bg_subtracted_bin_reco_binning,
                                     normalise_hist=True),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=detector_title,
                                ytitle=normalised_differential_label,
                                subplot_title='Data / MC',
                                **common_hist_args)
                    plot.legend.SetX1(0.6)
                    plot.legend.SetY1(0.75)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.9)
                    plot.plot("NOSTACK E1")
                    # plot.save("%s/detector_reco_binning_bg_subtracted_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                    # again but divided by bin width
                    mc_reco_hist_bg_subtracted_bin_reco_binning_div_bin_width = qgp.hist_divide_bin_width(mc_reco_hist_bg_subtracted_bin_reco_binning)
                    reco_hist_bg_subtracted_bin_reco_binning_div_bin_width = qgp.hist_divide_bin_width(reco_hist_bg_subtracted_bin_reco_binning)
                    entries = [
                        Contribution(mc_reco_hist_bg_subtracted_bin_reco_binning_div_bin_width,
                                     label="MC (bg-subtracted)",
                                     line_color=reco_mc_colour, line_width=lw,
                                     marker_color=reco_mc_colour, marker_size=0,
                                     normalise_hist=False),
                        Contribution(reco_hist_bg_subtracted_bin_reco_binning_div_bin_width,
                                     label="Data (bg-subtracted)",
                                     line_color=reco_data_colour, line_width=lw,
                                     marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                     subplot=mc_reco_hist_bg_subtracted_bin_reco_binning_div_bin_width,
                                     normalise_hist=False),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=detector_title,
                                ytitle=normalised_differential_label,
                                subplot_title='Data / MC',
                                **common_hist_args)
                    plot.legend.SetX1(0.6)
                    plot.legend.SetY1(0.75)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.9)
                    plot.plot("NOSTACK E1")
                    plot.save("%s/detector_reco_binning_bg_subtracted_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # PLOT FOLDED UNFOLDED DATA
                # --------------------------------------------------------------
                # Reco + folded, detector binning
                # FIXME: shouldn't this use RECO gen bin
                reco_mc_colour = ROOT.kGreen+2
                reco_data_colour = ROOT.kRed
                reco_folded_colour = ROOT.kAzure+1
                mc_reco_bin_hist = mc_reco_hist_bg_subtracted_bin_reco_binning if SUBTRACT_FAKES else mc_reco_hist_bin_reco_binning
                reco_bin_hist = reco_hist_bg_subtracted_bin_reco_binning if SUBTRACT_FAKES else reco_hist_bin_reco_binning

                entries = [
                    Contribution(mc_reco_bin_hist,
                                 label="MC (reco, bg-subtracted)" if SUBTRACT_FAKES else "MC (reco)",
                                 line_color=reco_mc_colour, line_width=lw,
                                 marker_color=reco_mc_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(reco_bin_hist,
                                 label="Data (reco, bg-subtracted)" if SUBTRACT_FAKES else "Data (reco)",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_size=0,
                                 subplot=mc_reco_bin_hist,
                                 normalise_hist=True),
                    Contribution(folded_unfolded_hist_bin_reco_binning,
                                 label="Folded unfolded data (#tau = %.3g)" % (tau),
                                 line_color=reco_folded_colour, line_width=lw,
                                 marker_color=reco_folded_colour, marker_size=0,
                                 subplot=mc_reco_bin_hist,
                                 normalise_hist=True),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Data / MC',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.72)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                # plot.save("%s/detector_folded_unfolded_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Folded, but only comparing data with data to check it is sane
                entries = [
                    Contribution(reco_bin_hist,
                                 label="Data (reco, bg-subtracted)" if SUBTRACT_FAKES else "Data (reco)",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(folded_unfolded_hist_bin_reco_binning,
                                 label="Folded unfolded data (#tau = %.3g)" % (tau),
                                 line_color=reco_folded_colour, line_width=lw,
                                 marker_color=reco_folded_colour, marker_size=0,
                                 subplot=reco_bin_hist,
                                 normalise_hist=True),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Folded unfolded / reco',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                # plot.save("%s/detector_folded_unfolded_only_data_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Same but divided by bin width
                # Do not normalise again!
                # data + folded + MC
                mc_reco_bin_hist_div_bin_width = qgp.hist_divide_bin_width(mc_reco_bin_hist)
                reco_bin_hist_div_bin_width = qgp.hist_divide_bin_width(reco_bin_hist)
                folded_unfolded_hist_bin_reco_binning_div_bin_width = qgp.hist_divide_bin_width(folded_unfolded_hist_bin_reco_binning)
                entries = [
                    Contribution(mc_reco_bin_hist_div_bin_width,
                                 label="MC (reco, bg-subtracted)" if SUBTRACT_FAKES else "MC (reco)",
                                 line_color=reco_mc_colour, line_width=lw,
                                 marker_color=reco_mc_colour, marker_size=0,
                                 normalise_hist=False),
                    Contribution(reco_bin_hist_div_bin_width,
                                 label="Data (reco, bg-subtracted)" if SUBTRACT_FAKES else "Data (reco)",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_size=0,
                                 subplot=mc_reco_bin_hist_div_bin_width,
                                 normalise_hist=False),
                    Contribution(folded_unfolded_hist_bin_reco_binning_div_bin_width,
                                 label="Folded unfolded data (#tau = %.3g)" % (tau),
                                 line_color=reco_folded_colour, line_width=lw,
                                 marker_color=reco_folded_colour, marker_size=0,
                                 subplot=mc_reco_bin_hist_div_bin_width,
                                 normalise_hist=False),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Data / MC',
                            **common_hist_args)
                plot.legend.SetX1(0.56)
                plot.legend.SetY1(0.72)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/detector_folded_unfolded_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # data + folded
                entries = [
                    Contribution(reco_bin_hist_div_bin_width,
                                 label="Data (reco, bg-subtracted)" if SUBTRACT_FAKES else "Data (reco)",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_size=0,
                                 normalise_hist=False),
                    Contribution(folded_unfolded_hist_bin_reco_binning_div_bin_width,
                                 label="Folded unfolded data (#tau = %.3g)" % (tau),
                                 line_color=reco_folded_colour, line_width=lw,
                                 marker_color=reco_folded_colour, marker_size=0,
                                 subplot=reco_bin_hist_div_bin_width,
                                 normalise_hist=False),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Folded unfolded / reco',
                            **common_hist_args)
                plot.legend.SetX1(0.56)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/detector_folded_unfolded_only_data_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # PLOT FOLDED GEN MC
                # --------------------------------------------------------------
                if MC_INPUT and hist_gen_folded:
                    # Get 1D hists
                    # The folded result
                    folded_gen_hist_bin_reco_binning = unfolder.get_var_hist_pt_binned(hist_gen_folded, ibin_pt, binning_scheme="detector")
                    this_pt_bin_tdir.WriteTObject(folded_gen_hist_bin_reco_binning, "folded_gen_hist_bin_reco_binning")

                    # Folded gen, comparing to original reco
                    entries = [
                        Contribution(mc_reco_bin_hist,
                                     label="MC (reco, bg-subtracted)" if SUBTRACT_FAKES else "MC (reco)",
                                     line_color=reco_data_colour, line_width=lw,
                                     marker_color=reco_data_colour, marker_size=0,
                                     normalise_hist=True),
                        Contribution(folded_gen_hist_bin_reco_binning,
                                     label="Folded particle-level MC",
                                     line_color=reco_folded_colour, line_width=lw,
                                     marker_color=reco_folded_colour, marker_size=0,
                                     subplot=reco_bin_hist,
                                     normalise_hist=True),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=detector_title,
                                ytitle=normalised_differential_label,
                                subplot_title='Folded gen / reco',
                                **common_hist_args)
                    plot.legend.SetX1(0.6)
                    plot.legend.SetY1(0.75)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.9)
                    plot.plot("NOSTACK E1")
                    # plot.save("%s/detector_folded_gen_only_data_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                    # Same but divided by bin width
                    # Do not normalise again!
                    mc_reco_bin_hist_div_bin_width = qgp.hist_divide_bin_width(mc_reco_bin_hist)
                    folded_gen_hist_bin_reco_binning_div_bin_width = qgp.hist_divide_bin_width(folded_gen_hist_bin_reco_binning)
                    entries = [
                        Contribution(mc_reco_bin_hist_div_bin_width,
                                     label="MC (reco, bg-subtracted)" if SUBTRACT_FAKES else "MC (reco)",
                                     line_color=reco_mc_colour, line_width=lw,
                                     marker_color=reco_mc_colour, marker_size=0,
                                     normalise_hist=False),
                        Contribution(folded_gen_hist_bin_reco_binning_div_bin_width,
                                     label="Folded particle-level MC",
                                     line_color=reco_folded_colour, line_width=lw,
                                     marker_color=reco_folded_colour, marker_size=0,
                                     subplot=mc_reco_bin_hist_div_bin_width,
                                     normalise_hist=False),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=detector_title,
                                ytitle=normalised_differential_label,
                                subplot_title='Folded gen / reco',
                                **common_hist_args)
                    plot.legend.SetX1(0.56)
                    plot.legend.SetY1(0.72)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.9)
                    plot.plot("NOSTACK E1")
                    plot.save("%s/detector_folded_gen_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

            # DO SUMMARY PLOT
            # ------------------------------------------------------------------
            if args.doSummaryPlot:
                marker = ""
                if "_" in angle.name or "^" in angle.name:
                    marker = "$"
                var_label = "Particle-level " + marker + angle.name + marker + " ($%s$)" % angle.lambda_str
                v = "%s_vs_pt" % (angle.var)
                bins = [(pt_bin_edges_gen[i], pt_bin_edges_gen[i+1]) for i in range(len(pt_bin_edges_gen)-1)]
                print("Making summary plot from pt bins:", bins)
                xlim = (50, 2000)
                if "ZPlusJets" in region['name']:
                    xlim = (50, 614)
                region_label = region['label'].replace("Dijet", "dijet")  # to ensure correct capitalisation
                qgp.do_mean_rms_summary_plot(summary_1d_entries, bins,
                                             "%s/%s_box_dijet_mpl.%s" % (this_output_dir, v, OUTPUT_FMT),
                                             var_label=var_label,
                                             xlim=xlim,
                                             region_title=region_label)

    print("Saved hists to", output_tfile.GetName())
