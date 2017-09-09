#!/usr/bin/env python

"""Print plots comparing different pileup subtraction methods"""

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


def do_chs_vs_puppi_plots(chs_dir, puppi_dir, plot_dir):
    sources = [
        {"root_dir": chs_dir, 'label': ", CHS", "style": {'line_style': 1}},
        {"root_dir": puppi_dir, 'label': ", PUPPI", "style": {'line_style': 3}}
    ]
    qgg.do_all_exclusive_plots_comparison(sources, var_list=qgc.COMMON_VARS,
                                          plot_dir=os.path.join(plot_dir, "chs_vs_puppi"),
                                          zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG",
                                          var_prepend="", pt_bins=None,
                                          subplot_type=None, do_flav_tagged=True, 
                                          ofmt=OUTPUT_FMT)


if __name__ == '__main__':
    for algo in ['ak4', 'ak8']:
        # only need to edit this:
        TEMPLATE = "workdir_{algo}{pus}_mgpythia"
        
        CHS_DIR = TEMPLATE.format(algo=algo, pus="chs")
        PUPPI_DIR = TEMPLATE.format(algo=algo, pus="puppi")
        PLOT_DIR = PUPPI_DIR
        do_chs_vs_puppi_plots(CHS_DIR, PUPPI_DIR, PLOT_DIR)
