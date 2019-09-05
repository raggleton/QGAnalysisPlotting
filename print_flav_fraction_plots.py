#!/usr/bin/env python

"""Print flavour fraction plots/graphs"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
import argparse

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

    pt_bins = qgc.PT_BINS_INC_UFLOW
    # Plots of all flavour fractions vs pT for a given sample/selection
    if zpj_dirname:
        # Z+jets
        qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.DY_FILENAME),
                                      title="Z+jets selection",
                                      dirname=zpj_dirname,
                                      pt_bins=pt_bins,
                                      flav_source=flav_source,
                                      var_prepend=var_prepend,
                                      output_filename="%s/zpj_flavour_fractions.%s" % (plot_dir, OUTPUT_FMT))
    if dj_dirname:
        # Dijets
        qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.QCD_FILENAME),
                                      title="Dijet selection (both jets)",
                                      dirname=dj_dirname,
                                      pt_bins=pt_bins,
                                      flav_source=flav_source,
                                      var_prepend=var_prepend,
                                      output_filename="%s/dj_flavour_fractions.%s" % (plot_dir, OUTPUT_FMT))
        qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.QCD_FILENAME),
                                      title="Dijet selection (jet 1)",
                                      dirname=dj_dirname,
                                      pt_bins=pt_bins,
                                      flav_source=flav_source,
                                      var_prepend=var_prepend,
                                      which_jet="1",
                                      output_filename="%s/dj_flavour_fractions_jet1.%s" % (plot_dir, OUTPUT_FMT))
        qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.QCD_FILENAME),
                                      title="Dijet selection (jet 2)",
                                      dirname=dj_dirname,
                                      pt_bins=pt_bins,
                                      flav_source=flav_source,
                                      var_prepend=var_prepend,
                                      which_jet="2",
                                      output_filename="%s/dj_flavour_fractions_jet2.%s" % (plot_dir, OUTPUT_FMT))

    if dj_dirname and zpj_dirname:
        dirnames = [dj_dirname, zpj_dirname]
        labels = ["Dijet", "Z+jets"]
        this_flav = "g"
        # Compare gluon fractions across samples/selections
        qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(root_dir, qgc.QCD_FILENAME), os.path.join(root_dir, qgc.DY_FILENAME)],
                                            dirnames=dirnames,
                                            pt_bins=pt_bins,
                                            labels=labels,
                                            flav=this_flav,
                                            output_filename="%s/g_flav_fraction_compare_bothjets.%s" % (plot_dir, OUTPUT_FMT),
                                            flav_source=flav_source,
                                            var_prepend=var_prepend)
        qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(root_dir, qgc.QCD_FILENAME), os.path.join(root_dir, qgc.DY_FILENAME)],
                                            dirnames=dirnames,
                                            pt_bins=pt_bins,
                                            labels=labels,
                                            flav=this_flav,
                                            output_filename="%s/g_flav_fraction_compare_jet1.%s" % (plot_dir, OUTPUT_FMT),
                                            flav_source=flav_source,
                                            var_prepend=var_prepend,
                                            which_jet="1",
                                            xtitle="p_{T}^{jet 1} [GeV]")
        qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(root_dir, qgc.QCD_FILENAME), os.path.join(root_dir, qgc.DY_FILENAME)],
                                            dirnames=dirnames,
                                            pt_bins=pt_bins,
                                            labels=labels,
                                            flav=this_flav,
                                            output_filename="%s/g_flav_fraction_compare_jet2.%s" % (plot_dir, OUTPUT_FMT),
                                            flav_source=flav_source,
                                            var_prepend=var_prepend,
                                            which_jet="2",
                                            xtitle="p_{T}^{jet 2} [GeV]")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--dj", help="do dijet plots", action="store_true")
    parser.add_argument("--zpj", help="do z + jets plots", action="store_true")
    parser.add_argument("--genparton", help="Use genparton flavour instead of partonflavour", action="store_true")
    parser.add_argument("dirs", help="Workdir(s) with ROOT files", nargs="+")
    args = parser.parse_args()
    for root_dir in args.dirs:
        if args.genparton:
            do_all_flavour_fraction_plots(root_dir, 
                plot_dir=os.path.join(root_dir, "flav_fractions_genParton"), 
                flav_source="genParton_",
                zpj_dirname="ZPlusJets_QG" if args.zpj else None,
                dj_dirname="Dijet_QG" if args.dj else None)
        else:
            do_all_flavour_fraction_plots(root_dir, 
                plot_dir=os.path.join(root_dir, "flav_fractions"), 
                flav_source="",
                zpj_dirname="ZPlusJets_QG" if args.zpj else None,
                dj_dirname="Dijet_QG" if args.dj else None)
