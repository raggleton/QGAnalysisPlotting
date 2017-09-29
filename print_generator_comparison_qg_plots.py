#!/usr/bin/env python

"""Print plots comparing different generators"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
from itertools import product

# My stuff
import qg_common as qgc
import qg_general_plots as qgg
import qg_delta_plots as qgd
import qg_roc_plots as qgr

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


# Control legend labelling
PYTHIA_LABEL = ", MG+Pythia"
HERWIG_LABEL = ", Herwig"

# Control output format
OUTPUT_FMT = "pdf"

def do_reco_generator_comparison_plots(pythia_dir, herwig_dir, plot_dir):
    """Compare reco jets from different generators"""
    
    # Distributions
    sources = [
        {"root_dir": pythia_dir, 'label': PYTHIA_LABEL, "style": {'line_style': 1}},
        {"root_dir": herwig_dir, 'label': HERWIG_LABEL, "style": {'line_style': 2}}

    ]
    qgg.do_all_exclusive_plots_comparison(sources=sources, var_list=qgc.COMMON_VARS,
                                          plot_dir=os.path.join(plot_dir, "plots_dy_vs_qcd_compare_generators"),
                                          subplot_type=None, do_flav_tagged=False, pt_bins=qgc.THEORY_PT_BINS,
                                          ofmt=OUTPUT_FMT)

    # Delta plots
    sources = [
        {"root_dir": pythia_dir, 'label': PYTHIA_LABEL, "style": {'line_style': 1, 'line_color': ROOT.kBlack}},
        {"root_dir": herwig_dir, 'label': HERWIG_LABEL, "style": {'line_style': 2, 'line_color': ROOT.kRed}}

    ]
    # qgd.do_pt_min_delta_plots(sources, var_list=qgc.COMMON_VARS)
    qgd.do_angularity_delta_plots(sources, plot_dir=os.path.join(plot_dir, "delta_angularities_compare_generators"),
                                  var_list=qgc.COMMON_VARS, pt_bins=qgc.THEORY_PT_BINS,
                                  ofmt=OUTPUT_FMT)



def do_gen_generator_comparison_plots(pythia_dir, herwig_dir, plot_dir):
    """Compare genjets from different generators"""

    # Distributions
    sources = [
        {"root_dir": pythia_dir, 'label': PYTHIA_LABEL, "style": {'line_style': 1}},
        {"root_dir": herwig_dir, 'label': HERWIG_LABEL, "style": {'line_style': 2}}
    ]
    qgg.do_all_exclusive_plots_comparison(sources=sources, var_list=qgc.COMMON_VARS, var_prepend="gen",
                                          plot_dir=os.path.join(plot_dir, "plots_dy_vs_qcd_gen_compare_generators"),
                                          zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR,
                                          pt_bins=qgc.THEORY_PT_BINS, subplot_type=None, do_flav_tagged=False,
                                          ofmt=OUTPUT_FMT)

    # Delta plots
    sources = [
        {"root_dir": pythia_dir, 'label': PYTHIA_LABEL, "style": {'line_color': ROOT.kBlack}},
        {"root_dir": herwig_dir, 'label': HERWIG_LABEL, "style": {'line_color': ROOT.kRed, 'line_style': 2}}
    ]
    qgd.do_pt_min_delta_plots(sources, var_list=qgc.COMMON_VARS, var_prepend="gen",
                          plot_dir=os.path.join(plot_dir, "deltas_ptMin_gen_compare_generators"),
                          zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR,
                          ofmt=OUTPUT_FMT)
    qgd.do_angularity_delta_plots(sources, var_list=qgc.COMMON_VARS, var_prepend="gen",
                                  plot_dir=os.path.join(plot_dir, "deltas_angularities_gen_compare_generators"),
                                  zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR, pt_bins=qgc.THEORY_PT_BINS,
                                  ofmt=OUTPUT_FMT)


def do_jet_algo_generator_comparison_plots(pythia_ak4_dir, herwig_ak4_dir, pythia_ak8_dir, herwig_ak8_dir, plot_dir):
    """Do ak4 vs ak8, for both pythia and herwig"""
    
    # Distributions
    sources = [
        {"root_dir": pythia_ak4_dir, 'label': PYTHIA_LABEL + " AK4", "style": {'line_style': 1}},
        {"root_dir": herwig_ak4_dir, 'label': HERWIG_LABEL + " AK4", "style": {'line_style': 2}},
        {"root_dir": pythia_ak8_dir, 'label': PYTHIA_LABEL + " AK8", "style": {'line_style': 1}, 
          "dy_style": {'line_color': qgc.DY_COLOURS[-2], 'fill_color': qgc.DY_COLOURS[-2]},
          "qcd_style": {'line_color': qgc.QCD_COLOURS[-2], 'fill_color': qgc.QCD_COLOURS[-2]}},
        {"root_dir": herwig_ak8_dir, 'label': HERWIG_LABEL + " AK8", "style": {'line_style': 2},
          "dy_style": {'line_color': qgc.DY_COLOURS[-2], 'fill_color': qgc.DY_COLOURS[-2]},
          "qcd_style": {'line_color': qgc.QCD_COLOURS[-2], 'fill_color': qgc.QCD_COLOURS[-2]}},
    ]

    # qgg.do_all_exclusive_plots_comparison(sources=sources, var_list=qgc.COMMON_VARS,
    #                                       plot_dir=os.path.join(plot_dir, "plots_dy_vs_qcd_compare_ak4_ak8_generators"),
    #                                       subplot_type=None, do_flav_tagged=False, pt_bins=qgc.THEORY_PT_BINS,
    #                                       ofmt=OUTPUT_FMT)
    sources = [
        {"root_dir": pythia_ak4_dir, 'label': PYTHIA_LABEL + " AK4", "style": {'line_width': 1, 'line_style': 1, 'line_color': ROOT.kBlue}},
        {"root_dir": herwig_ak4_dir, 'label': HERWIG_LABEL + " AK4", "style": {'line_width': 1, 'line_style': 1, 'line_color': ROOT.kRed}},
        {"root_dir": pythia_ak8_dir, 'label': PYTHIA_LABEL + " AK8", "style": {'line_width': 1, 'line_style': 2, 'line_color': ROOT.kBlue}}, 
        {"root_dir": herwig_ak8_dir, 'label': HERWIG_LABEL + " AK8", "style": {'line_width': 1, 'line_style': 2, 'line_color': ROOT.kRed}},
    ]
    qgd.do_angularity_delta_plots(sources, plot_dir=os.path.join(plot_dir, "delta_angularities_compare_ak4_ak8_generators"),
                                  var_list=qgc.COMMON_VARS, pt_bins=qgc.THEORY_PT_BINS, save_component_hists=True,
                                  ofmt=OUTPUT_FMT)
    # qgr.do_angularity_roc_plots(sources, plot_dir=os.path.join(plot_dir, "roc_angularities_compare_ak4_ak8_generators"),
    #                             pt_bins=qgc.THEORY_PT_BINS, var_list=qgc.COMMON_VARS)


if __name__ == '__main__':
    """
    ALGOS = ["ak4", "ak8"]
    PUS = ["chs", "puppi"]

    for algo, pu, in product(ALGOS, PUS):
        SETUP = algo + pu

        PYTHIA_DIR = "workdir_%s_mgpythia" % SETUP
        HERWIG_DIR = "workdir_%s_herwig" % SETUP
        PLOT_DIR = HERWIG_DIR

        do_reco_generator_comparison_plots(PYTHIA_DIR, HERWIG_DIR, PLOT_DIR)
        do_gen_generator_comparison_plots(PYTHIA_DIR, HERWIG_DIR, PLOT_DIR)
    do_jet_algo_generator_comparison_plots("workdir_ak4puppi_mgpythia", 
                                           "workdir_ak4puppi_herwig_reweight", 
                                           "workdir_ak8puppi_mgpythia", 
                                           "workdir_ak8puppi_herwig_reweight", 
                                           "workdir_ak8puppi_herwig_reweight")
