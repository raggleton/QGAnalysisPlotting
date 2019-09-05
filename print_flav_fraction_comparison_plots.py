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


def do_all_flavour_fraction_plots(root_dir, plot_dir="flav_fractions", zpj_dirname="ZPlusJets_QG", dj_cen_dirname="Dijet_QG_central_tighter", dj_fwd_dirname="Dijet_QG_forward_tighter", var_prepend=""):
    """Do plots of jet flavour fractions vs pT, for both Z+jets and dijets regions"""
    pt_bins = qgc.PT_BINS_INC_UFLOW
    # Plots of all flavour fractions vs pT for a given sample/selection
    if zpj_dirname:
        # Z+jets
        qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.DY_FILENAME),
                                      title="Z+jets selection",
                                      dirname=zpj_dirname,
                                      pt_bins=pt_bins,
                                      var_prepend=var_prepend,
                                      which_jet="1",
                                      output_filename="%s/zpj_flavour_fractions.%s" % (plot_dir, OUTPUT_FMT))
    if dj_cen_dirname:
        # Dijets central
        qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.QCD_FILENAME),
                                      title="Dijet selection (central jet)",
                                      dirname=dj_cen_dirname,
                                      pt_bins=pt_bins,
                                      var_prepend=var_prepend,
                                      which_jet="1",
                                      output_filename="%s/dj_flavour_fractions_central_jet.%s" % (plot_dir, OUTPUT_FMT))

    if dj_fwd_dirname:
        # Dijets central
        qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.QCD_FILENAME),
                                      title="Dijet selection (forward jet)",
                                      dirname=dj_fwd_dirname,
                                      pt_bins=pt_bins,
                                      var_prepend=var_prepend,
                                      which_jet="1",
                                      output_filename="%s/dj_flavour_fractions_forward_jet.%s" % (plot_dir, OUTPUT_FMT))


    dirnames = [dj_cen_dirname, dj_fwd_dirname, zpj_dirname]
    labels = ["Dijet (central jet)", "Dijet (forward jet)", "Z+jets"]
    this_dirnames = []
    this_labels = []
    for d, l in zip(dirnames, labels):
        if d:
            this_dirnames.append(d)
            this_labels.append(l)
    for this_flav in ['g', 'u', 'd', '1-g']:
        # this_flav = "g"
        # Compare gluon fractions across samples/selections
        input_files = [os.path.join(root_dir, qgc.QCD_FILENAME) if "Dijet" in d else os.path.join(root_dir, qgc.DY_FILENAME) for d in this_dirnames]
        qgf.compare_flavour_fractions_vs_pt(input_files=input_files,
                                            dirnames=dirnames,
                                            pt_bins=pt_bins,
                                            labels=labels,
                                            flav=this_flav,
                                            output_filename="%s/%s_flav_fraction_compare_jet1.%s" % (plot_dir, this_flav, OUTPUT_FMT),
                                            var_prepend=var_prepend,
                                            which_jet="1",
                                            xtitle="p_{T}^{jet} [GeV]")


def do_flavour_fraction_input_comparison_plots(root_dirs, labels, plot_dir="flav_fractions_comparison", zpj_dirname="ZPlusJets_QG", dj_cen_dirname="Dijet_QG_central_tighter", dj_fwd_dirname="Dijet_QG_forward_tighter", var_prepend=""):
    """Do plots comparing several input dirs """
    pt_bins = qgc.PT_BINS_INC_UFLOW
    if dj_cen_dirname:
        this_flav = "g"
        dirnames = [dj_cen_dirname]*len(root_dirs)
        # Compare gluon fractions across samples/selections
        qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.QCD_FILENAME) for rd in root_dirs],
                                            dirnames=dirnames,
                                            pt_bins=pt_bins,
                                            labels=labels,
                                            flav=this_flav,
                                            output_filename="%s/dj_g_flav_fraction_compare_central_jet.%s" % (plot_dir, OUTPUT_FMT),
                                            var_prepend=var_prepend,
                                            which_jet="1",
                                            xtitle="p_{T}^{jet} [GeV]")
        this_flav = "1-g"
        dirnames = [dj_cen_dirname]*len(root_dirs)
        # Compare non-gluon fractions across samples/selections
        qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.QCD_FILENAME) for rd in root_dirs],
                                            dirnames=dirnames,
                                            pt_bins=pt_bins,
                                            labels=labels,
                                            flav=this_flav,
                                            output_filename="%s/dj_q_flav_fraction_compare_central_jet.%s" % (plot_dir, OUTPUT_FMT),
                                            var_prepend=var_prepend,
                                            which_jet="1",
                                            xtitle="p_{T}^{jet 1} [GeV]")
    if dj_fwd_dirname:
        this_flav = "g"
        dirnames = [dj_fwd_dirname]*len(root_dirs)
        # Compare gluon fractions across samples/selections
        qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.QCD_FILENAME) for rd in root_dirs],
                                            dirnames=dirnames,
                                            pt_bins=pt_bins,
                                            labels=labels,
                                            flav=this_flav,
                                            output_filename="%s/dj_g_flav_fraction_compare_forward_jet.%s" % (plot_dir, OUTPUT_FMT),
                                            var_prepend=var_prepend,
                                            which_jet="1",
                                            xtitle="p_{T}^{jet} [GeV]")
        this_flav = "1-g"
        dirnames = [dj_fwd_dirname]*len(root_dirs)
        # Compare non-gluon fractions across samples/selections
        qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.QCD_FILENAME) for rd in root_dirs],
                                            dirnames=dirnames,
                                            pt_bins=pt_bins,
                                            labels=labels,
                                            flav=this_flav,
                                            output_filename="%s/dj_q_flav_fraction_compare_forward_jet.%s" % (plot_dir, OUTPUT_FMT),
                                            var_prepend=var_prepend,
                                            which_jet="1",
                                            xtitle="p_{T}^{jet 1} [GeV]")

    if zpj_dirname:
        this_flav = "g"
        dirnames = [zpj_dirname] * len(root_dirs)
        # Compare gluon fractions across samples/selections
        qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.DY_FILENAME) for rd in root_dirs],
                                            dirnames=dirnames,
                                            pt_bins=pt_bins,
                                            labels=labels,
                                            flav=this_flav,
                                            output_filename="%s/zpj_g_flav_fraction_compare.%s" % (plot_dir, OUTPUT_FMT),
                                            var_prepend=var_prepend)
        this_flav = "1-g"
        dirnames = [zpj_dirname] * len(root_dirs)
        # Compare quark fractions across samples/selections
        qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.DY_FILENAME) for rd in root_dirs],
                                            dirnames=dirnames,
                                            pt_bins=pt_bins,
                                            labels=labels,
                                            flav=this_flav,
                                            output_filename="%s/zpj_q_flav_fraction_compare.%s" % (plot_dir, OUTPUT_FMT),
                                            var_prepend=var_prepend)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--dj", help="do dijet plots", action="store_true")
    parser.add_argument("--zpj", help="do z + jets plots", action="store_true")
    parser.add_argument("--dir", help="Directory to get plot from. Can be used multiple times", action="append")
    parser.add_argument("--dirlabel", help="Label to be given for dir. Must be used in conjunction with --dir, once per entry.", action="append")
    args = parser.parse_args()
    print(args)

    # One set of plots per input
    for adir in args.dir:
        do_all_flavour_fraction_plots(adir,
            plot_dir=os.path.join(adir, "flav_fractions_genParton"),
            # zpj_dirname=None,
            zpj_dirname="ZPlusJets_QG",
            dj_cen_dirname="Dijet_QG_central_tighter",
            dj_fwd_dirname="Dijet_QG_forward_tighter")

    # Now comparison across all inputs
    do_flavour_fraction_input_comparison_plots(args.dir,
                                               plot_dir=os.path.join(args.dir[0], "flav_fractions_comparison"),
                                               labels=args.dirlabel,
                                               dj_cen_dirname="Dijet_QG_central_tighter" if args.dj else None,
                                               dj_fwd_dirname="Dijet_QG_forward_tighter" if args.dj else None,
                                               zpj_dirname="ZPlusJets_QG" if args.zpj else None)
