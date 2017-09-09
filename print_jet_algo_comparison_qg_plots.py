#!/usr/bin/env python

"""Print plots comparing different jet algorithms & cone sizes"""

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

# import hunter
# hunter.trace(module='comparator')

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)


# Control output format
OUTPUT_FMT = "pdf"


def do_jet_algo_comparison_plots(ak4_dir, ak8_dir, plot_dir="compare_jet_algo", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG",
                                 var_list=None, var_prepend="", pt_bins=None, subplot_type="diff"):
    """Do 1D plots comparing different jet algos"""
    sources = [
        {"root_dir": ak4_dir, 'label': ", AK4", "style": {'line_style': 1}},
        {"root_dir": ak8_dir, 'label': ", AK8", "style": {'line_style': 2}}
    ]
    qgg.do_all_exclusive_plots_comparison(sources, var_list=var_list,
                                          plot_dir=plot_dir,
                                          zpj_dirname=zpj_dirname, dj_dirname=dj_dirname,
                                          var_prepend=var_prepend, pt_bins=pt_bins,
                                          subplot_type=subplot_type, do_flav_tagged=True,
                                          ofmt=OUTPUT_FMT)

if __name__ == '__main__':
    for pus in ['chs', 'puppi']:
        # only need to edit this:
        TEMPLATE = "workdir_{algo}{pus}_mgpythia"
        
        AK4_DIR = TEMPLATE.format(algo="ak4", pus=pus)
        AK8_DIR = TEMPLATE.format(algo="ak8", pus=pus)
        PLOT_DIR = AK8_DIR
        
        # reco jets
        do_jet_algo_comparison_plots(AK4_DIR, AK8_DIR, var_list=qgc.COMMON_VARS,
                                     plot_dir=os.path.join(PLOT_DIR, "compare_jet_algo"))
        
        # gen jets
        do_jet_algo_comparison_plots(AK4_DIR, AK8_DIR, 
                                     var_list=qgc.COMMON_VARS, var_prepend="gen", 
                                     plot_dir=os.path.join(PLOT_DIR, "compare_jet_algo_gen"),
                                     zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR, 
                                     pt_bins=qgc.THEORY_PT_BINS, subplot_type=None)
