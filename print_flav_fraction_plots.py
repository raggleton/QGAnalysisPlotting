#!/usr/bin/env python

"""Print flavour fraction plots/graphs"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os

# My stuff
from comparator import Contribution, Plot, grab_obj
import qg_common as qgc
import qg_general_plots as qgg
import qg_flavour_plots as qgf

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


def do_all_flavour_fraction_plots(root_dir, plot_dir="flav_fractions", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG", var_prepend="", flav_source=""):
    """Do plots of jet flavour fractions vs pT, for both Z+jets and dijets regions"""
    
    # Plots of all flavour fractions vs pT for a given sample/selection
    # Z+jets
    qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.DY_FILENAME), 
                                  title="Z+jets selection",
                                  dirname=zpj_dirname, flav_source=flav_source, var_prepend=var_prepend,
                                  output_filename="%s/zpj_flavour_fractions.%s" % (plot_dir, OUTPUT_FMT))

    # Dijets
    qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.QCD_FILENAME), 
                                  title="Dijet selection",
                                  dirname=dj_dirname, flav_source=flav_source, var_prepend=var_prepend,
                                  output_filename="%s/dj_flavour_fractions.%s" % (plot_dir, OUTPUT_FMT)) 

    dirnames = [dj_dirname, zpj_dirname]
    labels = ["Dijet", "Z+jets"]
    this_flav = "g"
    # Compare gluon fractions across samples/selections
    qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(root_dir, qgc.QCD_FILENAME), os.path.join(root_dir, qgc.DY_FILENAME)],
                                        dirnames=dirnames,
                                        labels=labels,
                                        flav=this_flav,
                                        output_filename="%s/g_flav_fraction_compare_bothjets.%s" % (plot_dir, OUTPUT_FMT),
                                        flav_source=flav_source,
                                        var_prepend=var_prepend)
    qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(root_dir, qgc.QCD_FILENAME), os.path.join(root_dir, qgc.DY_FILENAME)],
                                        dirnames=dirnames,
                                        labels=labels,
                                        flav=this_flav,
                                        output_filename="%s/g_flav_fraction_compare_jet1.%s" % (plot_dir, OUTPUT_FMT),
                                        flav_source=flav_source,
                                        var_prepend=var_prepend,
                                        which_jet="1",
                                        xtitle="p_{T}^{jet 1} [GeV]")
    qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(root_dir, qgc.QCD_FILENAME), os.path.join(root_dir, qgc.DY_FILENAME)],
                                        dirnames=dirnames,
                                        labels=labels,
                                        flav=this_flav,
                                        output_filename="%s/g_flav_fraction_compare_jet2.%s" % (plot_dir, OUTPUT_FMT),
                                        flav_source=flav_source,
                                        var_prepend=var_prepend,
                                        which_jet="2",
                                        xtitle="p_{T}^{jet 2} [GeV]")


if __name__ == '__main__':
    for root_dir in sys.argv[1:]:
        do_all_flavour_fraction_plots(root_dir, plot_dir=os.path.join(root_dir, "flav_fractions"), flav_source="genParton_")
