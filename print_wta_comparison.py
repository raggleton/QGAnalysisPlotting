#!/usr/bin/env python

"""Print main lambda variable plots, comparing computation with/without WTA axis definition"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
os.nice(10)

# My stuff
from comparator import Contribution, Plot
import qg_common as qgc
import qg_general_plots as qgp
import common_utils as cu

# For debugging
import sys
import tracers
# sys.settrace(tracers.trace_calls)
# sys.settrace(tracers.trace_calls_detail)
# sys.settrace(tracers.trace_calls_and_returns)

# import hunter
# hunter.trace(module='comparator')

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)


# Control plot output format
OUTPUT_FMT = "pdf"

TOTAL_LUMI = 35918


if __name__ == "__main__":
    pythia_no_wta_dir = "workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_noZReweight"
    pythia_wta_dir = "workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_noZReweight_wta_groomed"
    sources = [
        {"root_dir": pythia_no_wta_dir, 'label': " [MG+PY8]\n(normal axis)", "style": {'line_color': ROOT.kBlue, 'marker_color': ROOT.kBlue}},
        {"root_dir": pythia_wta_dir , 'label': " [MG+PY8]\n(WTA axis)", "style": {'line_color': ROOT.kRed, 'marker_color': ROOT.kRed}},
    ]

    # Do Z+jets region only
    qgp.do_all_exclusive_plots_comparison(sources=sources,
                                          var_list=qgc.COMMON_VARS[:5],  # can't do charged-only as doesn't exist for non-wta
                                          plot_dir=os.path.join(pythia_wta_dir, "wta_vs_normal_jet_axis_zpj_ungroomed"),
                                          dy_filename=qgc.DY_FILENAME,
                                          zpj_dirname="ZPlusJets_QG",
                                          qcd_filename=None,
                                          dj_cen_dirname=None,
                                          dj_fwd_dirname=None,
                                          subplot_type="ratio",
                                          subplot_title="#splitline{WTA/normal}{jet axis}",
                                          do_flav_tagged=False,
                                          pt_bins=qgc.PT_BINS,
                                          ofmt=OUTPUT_FMT)

    # Do Dijet region only
    qgp.do_all_exclusive_plots_comparison(sources=sources,
                                          var_list=qgc.COMMON_VARS[:5],
                                          plot_dir=os.path.join(pythia_wta_dir, "wta_vs_normal_jet_axis_dijet_ungroomed"),
                                          dy_filename=None,
                                          zpj_dirname=None,
                                          qcd_filename=qgc.QCD_FILENAME,
                                          dj_dirname="Dijet_QG_tighter",
                                          dj_cen_dirname=None,
                                          dj_fwd_dirname=None,
                                          subplot_type="ratio",
                                          subplot_title="#splitline{WTA/normal}{jet axis}",
                                          do_flav_tagged=False,
                                          pt_bins=qgc.PT_BINS,
                                          ofmt=OUTPUT_FMT)

    # Do Pileup comparison
    pu_bins = [(5, 15), (20, 25), (30, 40)]
    sources = []
    for ind, (pu_min, pu_max) in enumerate(pu_bins):
        sources.append({
            "root_dir": pythia_wta_dir,
            'label': ", PU %d-%d (ungroomed)" % (pu_min, pu_max),
            'zpj_dirname': "ZPlusJets_QG_PU_%d_to_%d" % (pu_min, pu_max),
            'dj_dirname': "Dijet_QG_PU_%d_to_%d" % (pu_min, pu_max),
            "style": {'line_style': 1, "line_width": 2},
            "zpj_style": {'line_color': qgc.DY_COLOURS[ind], 'fill_color': qgc.DY_COLOURS[ind], 'marker_color': qgc.DY_COLOURS[ind]},
            "qcd_style": {'line_color': qgc.QCD_COLOURS[ind], 'fill_color': qgc.QCD_COLOURS[ind], 'marker_color': qgc.QCD_COLOURS[ind]}
        })
    subplot_title = "#splitline{Ratio wrt}{PU %d-%d}" % (pu_bins[0][0], pu_bins[0][1])
    qgp.do_all_exclusive_plots_comparison(sources=sources, 
                                          var_list=qgc.COMMON_VARS, 
                                          dj_cen_dirname=None,
                                          dj_fwd_dirname=None,
                                          plot_dir=os.path.join(pythia_wta_dir, "plots_dy_vs_qcd_compare_pu_zpj"),
                                          pt_bins=qgc.PT_BINS, 
                                          subplot_type="ratio", 
                                          subplot_title=subplot_title, 
                                          do_flav_tagged=False)
    
    qgp.do_all_exclusive_plots_comparison(sources=sources, var_list=qgc.COMMON_VARS, 
                                          zpj_dirname=None,
                                          plot_dir=os.path.join(pythia_wta_dir, "plots_dy_vs_qcd_compare_pu_dijet"),
                                          pt_bins=qgc.PT_BINS, 
                                          subplot_type="ratio", 
                                          subplot_title=subplot_title, 
                                          dj_dirname="blah",
                                          dj_cen_dirname=None,
                                          dj_fwd_dirname=None,
                                          do_flav_tagged=False)