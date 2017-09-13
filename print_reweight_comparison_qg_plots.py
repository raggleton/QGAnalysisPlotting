#!/usr/bin/env python

"""Print plots comparing pythia, herwig, and herwig reweighted to match pythia pt spectrum."""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
from itertools import product

# My stuff
import qg_common as qgc
import qg_general_plots as qgg
import qg_delta_plots as qgd

# For debugging
import sys
import tracers
# sys.settrace(tracers.trace_calls)
# sys.settrace(tracers.trace_calls_detail)
# sys.settrace(tracers.trace_calls_and_returns)


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)

# Control legend labelling
REWEIGHT_LABEL = ", Reweighted to MG+Pythia"
NOT_REWEIGHT_LABEL = ", Not reweighted"

# Control output format
OUTPUT_FMT = "pdf"


def do_reco_reweight_comparison_plots(pythia_dir, herwig_dir, herwig_reweighted_dir, plot_dir):
    """Compare effect of reweighting to Pythia spectrum - reco jets"""
    
    # Just do Herwig vs Herwig reweighted
    sources = [
        {"root_dir": herwig_reweighted_dir, 'label': REWEIGHT_LABEL, "style": {'line_style': 1}},
        {"root_dir": herwig_dir, 'label': NOT_REWEIGHT_LABEL, "style": {'line_style': 2}}
    ]
    qgg.do_all_exclusive_plots_comparison(sources=sources, var_list=qgc.COMMON_VARS,
                                          plot_dir=os.path.join(plot_dir, "plots_dy_vs_qcd_compare_reweight"),
                                          pt_bins=qgc.THEORY_PT_BINS, subplot_type=None, do_flav_tagged=False,
                                          ofmt=OUTPUT_FMT)

    sources = [
        {"root_dir": herwig_reweighted_dir, 'label': REWEIGHT_LABEL, "style": {'line_color': ROOT.kGreen+2}},
        {"root_dir": herwig_dir, 'label': NOT_REWEIGHT_LABEL, "style": {'line_color': ROOT.kAzure}}
    ]
    # qgd.do_pt_min_delta_plots(sources, var_list=qgc.COMMON_VARS, plot_dir=os.path.join(plot_dir, "deltas_ptMin_compare_reweight"))
    qgd.do_angularity_delta_plots(sources, var_list=qgc.COMMON_VARS,
                                  plot_dir=os.path.join(plot_dir, "deltas_angularities_compare_reweight"),
                                  pt_bins=qgc.THEORY_PT_BINS,
                                  ofmt=OUTPUT_FMT)
    
    # Now just do Pythia vs Herwig reweighted
    sources = [
        {"root_dir": pythia_dir, 'label': ", MG + Pythia", "style": {'line_style': 1}},
        {"root_dir": herwig_reweighted_dir, 'label': ", Herwig, reweighted to MG+Pythia", "style": {'line_style': 1},
          "dy_style": {'line_color': ROOT.kRed, 'fill_color': ROOT.kRed},
          "qcd_style": {'line_color': ROOT.kBlue, 'fill_color': ROOT.kBlue}
        }
    ]
    qgg.do_all_exclusive_plots_comparison(sources=sources, var_list=qgc.COMMON_VARS,
                                          plot_dir=os.path.join(plot_dir, "plots_dy_vs_qcd_compare_reweight_compare_pythia_herwig_reweight"),
                                          pt_bins=qgc.THEORY_PT_BINS, 
                                          subplot_type=None, 
                                          do_flav_tagged=False,
                                          ofmt=OUTPUT_FMT)

    # Now put Pythia, Herwig, & Herwig reweighted on same plots
    sources = [
        {"root_dir": pythia_dir, 'label': ", MG + Pythia", "style": {'line_style': 1}},
        {"root_dir": herwig_reweighted_dir, 'label': ", Herwig, reweighted to MG+Pythia",
            "style": {'line_style': 2},
            "dy_style": {'line_color': ROOT.kRed, 'fill_color': ROOT.kRed},
            "qcd_style": {'line_color': ROOT.kBlue, 'fill_color': ROOT.kBlue}
        },
        {"root_dir": herwig_dir, 'label': ", Herwig, not reweighted", "style": {'line_style': 2}}
    ]
    qgg.do_all_exclusive_plots_comparison(sources=sources, var_list=qgc.COMMON_VARS,
                                          plot_dir=os.path.join(plot_dir, "plots_dy_vs_qcd_compare_reweight_compare_generators"),
                                          pt_bins=qgc.THEORY_PT_BINS, 
                                          # subplot_type="diff", 
                                          subplot_type=None, 
                                          # subplot_title="#splitline{#Delta wrt}{%s}" % sources[0]['label'],
                                          subplot_title="#Delta",
                                          do_flav_tagged=False,
                                          ofmt=OUTPUT_FMT)

    sources = [
        {"root_dir": pythia_dir, 'label': ", MG + Pythia", "style": {'line_color': ROOT.kBlack}},
        {"root_dir": herwig_reweighted_dir, 'label': REWEIGHT_LABEL, "style": {'line_color': ROOT.kGreen+2}},
        {"root_dir": herwig_dir, 'label': NOT_REWEIGHT_LABEL, "style": {'line_color': ROOT.kAzure}}
    ]
    # qgd.do_pt_min_delta_plots(sources, var_list=qgc.COMMON_VARS, plot_dir=os.path.join(plot_dir, "deltas_ptMin_compare_normal_and_reweight"))
    qgd.do_angularity_delta_plots(sources, var_list=qgc.COMMON_VARS,
                                  plot_dir=os.path.join(plot_dir, "deltas_angularities_compare_normal_and_reweight"),
                                  pt_bins=qgc.THEORY_PT_BINS,
                                  ofmt=OUTPUT_FMT)


def do_gen_reweight_comparison_plots(pythia_dir, herwig_dir, herwig_reweighted_dir, plot_dir):
    """Compare effect of reweighting to Pythia spectrum - genjets"""
    
    # Just do Herwig vs Herwig reweighted
    sources = [
        {"root_dir": herwig_reweighted_dir, 'label': REWEIGHT_LABEL, "style": {'line_style': 1}},
        {"root_dir": herwig_dir, 'label': NOT_REWEIGHT_LABEL, "style": {'line_style': 2}}
    ]
    qgg.do_all_exclusive_plots_comparison(sources=sources, var_list=qgc.COMMON_VARS, var_prepend="gen",
                                          plot_dir=os.path.join(plot_dir, "plots_dy_vs_qcd_gen_compare_reweight"),
                                          zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR,
                                          pt_bins=qgc.THEORY_PT_BINS, subplot_type=None, do_flav_tagged=False,
                                          ofmt=OUTPUT_FMT)

    sources = [
        {"root_dir": herwig_reweighted_dir, 'label': REWEIGHT_LABEL, "style": {'line_color': ROOT.kGreen+2}},
        {"root_dir": herwig_dir, 'label': NOT_REWEIGHT_LABEL, "style": {'line_color': ROOT.kAzure}}
    ]
    qgd.do_pt_min_delta_plots(sources, var_list=qgc.COMMON_VARS, var_prepend="gen",
                              plot_dir=os.path.join(plot_dir, "deltas_ptMin_gen_compare_reweight"),
                              zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR,
                              ofmt=OUTPUT_FMT)
    qgd.do_angularity_delta_plots(sources, var_list=qgc.COMMON_VARS, var_prepend="gen",
                                  plot_dir=os.path.join(plot_dir, "deltas_angularities_gen_compare_reweight"),
                                  zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR, pt_bins=qgc.THEORY_PT_BINS,
                                  ofmt=OUTPUT_FMT)

    # Now jsut do Pythia vs Herwig reweighted
    sources = [
        {"root_dir": pythia_dir, 'label': ", MG + Pythia", "style": {'line_style': 1}},
        {"root_dir": herwig_reweighted_dir, 'label': ", Herwig, reweighted to MG+Pythia", "style": {'line_style': 1},
          "dy_style": {'line_color': ROOT.kRed, 'fill_color': ROOT.kRed},
          "qcd_style": {'line_color': ROOT.kBlue, 'fill_color': ROOT.kBlue}
        }
    ]
    qgg.do_all_exclusive_plots_comparison(sources=sources, var_list=qgc.COMMON_VARS, var_prepend="gen",
                                          plot_dir=os.path.join(plot_dir, "plots_dy_vs_qcd_gen_compare_reweight_compare_pythia_herwig_reweight"),
                                          zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR,
                                          pt_bins=qgc.THEORY_PT_BINS, 
                                          subplot_type=None, 
                                          do_flav_tagged=False,
                                          ofmt=OUTPUT_FMT)


    # Now put Pythia, Herwig, & Herwig reweighted on same plots
    sources = [
        {"root_dir": pythia_dir, 'label': ", MG + Pythia", "style": {'line_style': 1}},
        {"root_dir": herwig_reweighted_dir, 'label': ", Herwig, reweighted to MG+Pythia",
            "style": {'line_style': 2},
            "dy_style": {'line_color': ROOT.kRed, 'fill_color': ROOT.kRed},
            "qcd_style": {'line_color': ROOT.kBlue, 'fill_color': ROOT.kBlue}
        },
        {"root_dir": herwig_dir, 'label': ", Herwig, not reweighted", "style": {'line_style': 2}}
    ]
    qgg.do_all_exclusive_plots_comparison(sources=sources, var_list=qgc.COMMON_VARS, var_prepend="gen",
                                          plot_dir=os.path.join(plot_dir, "plots_dy_vs_qcd_gen_compare_reweight_compare_generators"),
                                          zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR,
                                          pt_bins=qgc.THEORY_PT_BINS, 
                                          subplot_type=None, 
                                          # subplot_title="#splitline{#Delta wrt}{%s}" % sources[0]['label'],
                                          subplot_title="#Delta",
                                          do_flav_tagged=False,
                                          ofmt=OUTPUT_FMT)
    sources = [
        {"root_dir": pythia_dir, 'label': "Pythia", "style": {'line_color': ROOT.kBlack}},
        {"root_dir": herwig_reweighted_dir, 'label': REWEIGHT_LABEL, "style": {'line_color': ROOT.kGreen+2}},
        {"root_dir": herwig_dir, 'label': NOT_REWEIGHT_LABEL, "style": {'line_color': ROOT.kAzure}}
    ]
    qgd.do_pt_min_delta_plots(sources, var_list=qgc.COMMON_VARS, var_prepend="gen",
                              plot_dir=os.path.join(plot_dir, "deltas_ptMin_gen_compare_normal_and_reweight"),
                              zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR,
                              ofmt=OUTPUT_FMT)
    qgd.do_angularity_delta_plots(sources, var_list=qgc.COMMON_VARS, var_prepend="gen",
                                  plot_dir=os.path.join(plot_dir, "deltas_angularities_gen_compare_normal_and_reweight"),
                                  zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR, pt_bins=qgc.THEORY_PT_BINS,
                                  ofmt=OUTPUT_FMT)


if __name__ == '__main__':
    ALGOS = ["ak4", "ak8"]
    PUS = ["chs", "puppi"]

    for algo, pu, in product(ALGOS, PUS):
        SETUP = algo + pu

        PYTHIA_DIR = "workdir_%s_mgpythia" % SETUP
        HERWIG_DIR = "workdir_%s_herwig" % SETUP
        HERWIG_REWEIGHTED_DIR = "workdir_%s_herwig_reweight" % SETUP
        PLOT_DIR = HERWIG_REWEIGHTED_DIR

        do_reco_reweight_comparison_plots(PYTHIA_DIR, HERWIG_DIR, HERWIG_REWEIGHTED_DIR, PLOT_DIR)
        do_gen_reweight_comparison_plots(PYTHIA_DIR, HERWIG_DIR, HERWIG_REWEIGHTED_DIR, PLOT_DIR)
