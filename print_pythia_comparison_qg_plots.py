#!/usr/bin/env python

"""Print plots comparing MG+Pythia with plain Pythia."""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
from itertools import product

# My stuff
import qg_common as qgc
import qg_general_plots as qgg
import qg_flavour_plots as qgf

# For debugging
# If your code segfaults, try enabling one of these, sometimes it helps it to run...
import sys
import tracers
# sys.settrace(tracers.trace_calls)
# sys.settrace(tracers.trace_calls_detail)
# sys.settrace(tracers.trace_calls_and_returns)


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)


OUTPUT_FMT = "pdf"


def do_pythia_comparison_distribution_plots(mgpythia_dir, pythia_only_dir, plot_dir):
    """To compare mg+pythia vs pythia distributions"""
    # reco
    sources = [
        {"root_dir": mgpythia_dir, 'label': "Madgraph+Pythia", "style": {'line_style': 1}},
        {"root_dir": pythia_only_dir, 'label': "Pythia only", "style": {'line_style': 1, 'line_color': ROOT.kRed, 'fill_color': ROOT.kRed, }}
    ]
    qgg.do_all_exclusive_plots_comparison(sources, var_list=qgc.COMMON_VARS_WITH_FLAV,
                                          plot_dir=os.path.join(plot_dir, "mg_pythia_vs_pythia_only"),
                                          zpj_dirname="",
                                          subplot_type=None, do_flav_tagged=False,
                                          pt_bins=qgc.THEORY_PT_BINS)

    # gen
    sources = [
        {"root_dir": mgpythia_dir, 'label': "Madgraph+Pythia", "style": {'line_style': 1}},
        {"root_dir": pythia_only_dir, 'label': "Pythia only", "style": {'line_style': 1, 'line_color': ROOT.kRed, 'fill_color': ROOT.kRed}}
    ]
    qgg.do_all_exclusive_plots_comparison(sources, var_list=qgc.COMMON_VARS_WITH_FLAV, var_prepend="gen",
                                          plot_dir=os.path.join(plot_dir, "mg_pythia_vs_pythia_only_gen"),
                                          dj_dirname=qgc.DJ_GENJET_RDIR, zpj_dirname="",
                                          subplot_type=None, do_flav_tagged=False,
                                          pt_bins=qgc.THEORY_PT_BINS)


def do_pythia_comparison_flav_fractions_plots(mgpythia_dir, pythia_only_dir, plot_dir, algo, pus):
    """Make plots of gluon jet fraction as a function of pT, both gen & reco."""
    
    # flavour fractions
    input_files = [
        os.path.join(mgpythia_dir, 'uhh2.AnalysisModuleRunner.MC.MC_QCD_.root'),
        os.path.join(pythia_only_dir, 'uhh2.AnalysisModuleRunner.MC.MC_QCD_.root')
    ]
    # reco
    qgf.compare_flavour_fractions_vs_pt(input_files,
                                        dirnames=[qgc.DJ_RECOJET_RDIR, qgc.DJ_RECOJET_RDIR],
                                        labels=[qgc.QCD_Dijet_LABEL + " Madgraph+Pythia", qgc.QCD_Dijet_LABEL+" Pythia only"],
                                        flav="g", 
                                        output_filename="%s/flav_fractions/compare_g_frac.%s" % (plot_dir, OUTPUT_FMT),
                                        title="%s PF %s jets" % (algo, pus.upper()), var_prepend="")
    # gen
    qgf.compare_flavour_fractions_vs_pt(input_files,
                                        dirnames=[qgc.DJ_GENJET_RDIR, qgc.DJ_GENJET_RDIR],
                                        labels=[qgc.QCD_Dijet_LABEL + " Madgraph+Pythia", qgc.QCD_Dijet_LABEL+" Pythia only"],
                                        flav="g", 
                                        output_filename="%s/flav_fractions_gen/compare_g_frac.%s" % (plot_dir, OUTPUT_FMT),
                                        title="%s GenJets" % (algo), var_prepend="gen")


if __name__ == '__main__':

    ALGOS = ["ak4", "ak8"]
    PUS = ["chs", "puppi"]

    for algo, pu, in product(ALGOS, PUS):
        SETUP = algo + pu

        MGPYTHIA_DIR = "workdir_%s_mgpythia" % SETUP
        PYTHIA_ONLY_DIR = "workdir_%s_pythiaOnlyFlat" % SETUP
        PLOT_DIR = PYTHIA_ONLY_DIR

        do_pythia_comparison_distribution_plots(MGPYTHIA_DIR, PYTHIA_ONLY_DIR, PLOT_DIR)
        do_pythia_comparison_flav_fractions_plots(MGPYTHIA_DIR, PYTHIA_ONLY_DIR, PLOT_DIR, algo, pu)
