#!/usr/bin/env python

"""Print flavour fraction plots/graphs"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
os.nice(10)
import argparse
import numpy as np

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

MG_SAMPLE = "MG5+Pythia8"
HPP_SAMPLE = "Herwig++"

def do_all_flavour_fraction_plots(root_dir,
                                  plot_dir="flav_fractions",
                                  use_gen=True,
                                  zpj_dirname="ZPlusJets_QG",
                                  dj_cen_dirname="Dijet_QG_central_tighter",
                                  dj_fwd_dirname="Dijet_QG_forward_tighter",
                                  var_prepend=""):
    """Do plots of jet flavour fractions vs pT, eta for both Z+jets and dijets regions"""
    # pt_bins = qgc.PT_BINS_INC_UFLOW
    pt_bins = qgc.PT_BINS_ZPJ

    # Plots of all flavour fractions vs pT for a given sample/selection
    # --------------------------------------------------------------------------
    if zpj_dirname:
        # Z+jets
        if os.path.isfile(os.path.join(root_dir, qgc.DY_FILENAME)):
            qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.DY_FILENAME),
                                          title=qgc.ZpJ_LABEL + "\n" + MG_SAMPLE,
                                          dirname=zpj_dirname,
                                          pt_bins=pt_bins,
                                          var_prepend=var_prepend,
                                          which_jet="1",
                                          output_filename="%s/zpj_flavour_fractions.%s" % (plot_dir, OUTPUT_FMT))

        if os.path.isfile(os.path.join(root_dir, qgc.DY_HERWIG_FILENAME)):
            qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.DY_HERWIG_FILENAME),
                                          title=qgc.ZpJ_LABEL + "\n" + HPP_SAMPLE,
                                          dirname=zpj_dirname,
                                          pt_bins=pt_bins,
                                          var_prepend=var_prepend,
                                          which_jet="1",
                                          output_filename="%s/zpj_flavour_fractions_herwigpp.%s" % (plot_dir, OUTPUT_FMT))
    if dj_cen_dirname:
        # Dijets central
        if os.path.isfile(os.path.join(root_dir, qgc.QCD_FILENAME)):
            qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.QCD_FILENAME),
                                          title=qgc.Dijet_CEN_LABEL + "\n" + MG_SAMPLE,
                                          dirname=dj_cen_dirname,
                                          pt_bins=pt_bins,
                                          var_prepend=var_prepend,
                                          which_jet="1",
                                          output_filename="%s/dj_flavour_fractions_central_jet.%s" % (plot_dir, OUTPUT_FMT))

        if os.path.isfile(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME)):
            qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME),
                                          title=qgc.Dijet_CEN_LABEL + "\n" + HPP_SAMPLE,
                                          dirname=dj_cen_dirname,
                                          pt_bins=pt_bins,
                                          var_prepend=var_prepend,
                                          which_jet="1",
                                          output_filename="%s/dj_flavour_fractions_central_jet_herwigpp.%s" % (plot_dir, OUTPUT_FMT))

    if dj_fwd_dirname:
        # Dijets central
        if os.path.isfile(os.path.join(root_dir, qgc.QCD_FILENAME)):
            qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.QCD_FILENAME),
                                          title=qgc.Dijet_FWD_LABEL + "\n" + MG_SAMPLE,
                                          dirname=dj_fwd_dirname,
                                          pt_bins=pt_bins,
                                          var_prepend=var_prepend,
                                          which_jet="1",
                                          output_filename="%s/dj_flavour_fractions_forward_jet.%s" % (plot_dir, OUTPUT_FMT))

        if os.path.isfile(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME)):
          qgf.do_flavour_fraction_vs_pt(input_file=os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME),
                                          title=qgc.Dijet_FWD_LABEL + "\n" + HPP_SAMPLE,
                                          dirname=dj_fwd_dirname,
                                          pt_bins=pt_bins,
                                          var_prepend=var_prepend,
                                          which_jet="1",
                                          output_filename="%s/dj_flavour_fractions_forward_jet_herwigpp.%s" % (plot_dir, OUTPUT_FMT))

    # Plots of all flavour fractions vs eta for a given sample/selection
    # --------------------------------------------------------------------------
    # Do for separate pt regions
    pt_regions = [
        {
            'append': "_lowPt",
            'title': '30 < p_{T}^{jet} < 100 GeV',
        },
        {
            'append': "_midPt",
            'title': '100 < p_{T}^{jet} < 250 GeV',
        },
        {
            'append': "_highPt",
            'title': '250 < p_{T}^{jet} < 500 GeV',
        },
        {
            'append': "_highPt2",
            'title': 'p_{T}^{jet} > 500 GeV',
        },
    ]
    end = 1.7
    interval = 0.25
    eta_bins = np.arange(-end, end+interval, interval)
    eta_bins = list(zip(eta_bins[:-1], eta_bins[1:])) # make pairwise bins
    print(eta_bins)
    for pt_region in pt_regions:
        eta_title = pt_region['title']
        append = pt_region['append']

        if zpj_dirname:
            # Z+jets
            if os.path.isfile(os.path.join(root_dir, qgc.DY_FILENAME)):
                qgf.do_flavour_fraction_vs_eta(input_file=os.path.join(root_dir, qgc.DY_FILENAME),
                                               title=qgc.ZpJ_LABEL +  "\n" + eta_title +  "\n" + MG_SAMPLE,
                                               dirname=zpj_dirname,
                                               eta_bins=eta_bins,
                                               var_prepend=var_prepend,
                                               which_jet="both",
                                               append=append,
                                               output_filename="%s/zpj_flavour_fractions_vs_eta%s.%s" % (plot_dir, append, OUTPUT_FMT))

            if os.path.isfile(os.path.join(root_dir, qgc.DY_HERWIG_FILENAME)):
                qgf.do_flavour_fraction_vs_eta(input_file=os.path.join(root_dir, qgc.DY_HERWIG_FILENAME),
                                               title=qgc.ZpJ_LABEL +  "\n" + eta_title + "\n" + HPP_SAMPLE,
                                               dirname=zpj_dirname,
                                               eta_bins=eta_bins,
                                               var_prepend=var_prepend,
                                               which_jet="both",
                                               append=append,
                                               output_filename="%s/zpj_flavour_fractions_vs_eta%s_herwigpp.%s" % (plot_dir, append, OUTPUT_FMT))

        if dj_cen_dirname:
            # Dijets central
            if os.path.isfile(os.path.join(root_dir, qgc.QCD_FILENAME)):
                qgf.do_flavour_fraction_vs_eta(input_file=os.path.join(root_dir, qgc.QCD_FILENAME),
                                               title=qgc.Dijet_CEN_LABEL +  "\n" + eta_title +  "\n" + MG_SAMPLE,
                                               dirname=dj_cen_dirname,
                                               eta_bins=eta_bins,
                                               var_prepend=var_prepend,
                                               which_jet="both",
                                               append=append,
                                               output_filename="%s/dj_flavour_fractions_central_jet_vs_eta%s.%s" % (plot_dir, append, OUTPUT_FMT))

            if os.path.isfile(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME)):
                qgf.do_flavour_fraction_vs_eta(input_file=os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME),
                                               title=qgc.Dijet_CEN_LABEL +  "\n" + eta_title + "\n" + HPP_SAMPLE,
                                               dirname=dj_cen_dirname,
                                               eta_bins=eta_bins,
                                               var_prepend=var_prepend,
                                               which_jet="both",
                                               append=append,
                                               output_filename="%s/dj_flavour_fractions_central_jet_vs_eta%s_herwigpp.%s" % (plot_dir, append, OUTPUT_FMT))

        if dj_fwd_dirname:
            # Dijets central
            if os.path.isfile(os.path.join(root_dir, qgc.QCD_FILENAME)):
                qgf.do_flavour_fraction_vs_eta(input_file=os.path.join(root_dir, qgc.QCD_FILENAME),
                                               title=qgc.Dijet_FWD_LABEL +  "\n" + eta_title +  "\n" + MG_SAMPLE,
                                               dirname=dj_fwd_dirname,
                                               eta_bins=eta_bins,
                                               var_prepend=var_prepend,
                                               which_jet="both",
                                               append=append,
                                               output_filename="%s/dj_flavour_fractions_forward_jet_vs_eta%s.%s" % (plot_dir, append, OUTPUT_FMT))

            if os.path.isfile(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME)):
                qgf.do_flavour_fraction_vs_eta(input_file=os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME),
                                               title=qgc.Dijet_FWD_LABEL +  "\n" + eta_title + "\n" + HPP_SAMPLE,
                                               dirname=dj_fwd_dirname,
                                               eta_bins=eta_bins,
                                               var_prepend=var_prepend,
                                               which_jet="both",
                                               append=append,
                                               output_filename="%s/dj_flavour_fractions_forward_jet_vs_eta%s_herwigpp.%s" % (plot_dir, append, OUTPUT_FMT))

    # Compare one flav fraction across samples/selections
    # --------------------------------------------------------------------------
    dirnames = [dj_cen_dirname, dj_fwd_dirname, zpj_dirname]
    for this_flav in ['g', 'u', 'd', '1-g'][:-1]:
        input_files = [os.path.join(root_dir, qgc.QCD_FILENAME) if "dijet" in d.lower()
                       else os.path.join(root_dir, qgc.DY_FILENAME)
                       for d in dirnames if d is not None]

        labels = [qgc.Dijet_CEN_LABEL, qgc.Dijet_FWD_LABEL, qgc.ZpJ_LABEL]

        qgf.compare_flavour_fractions_vs_pt(input_files=input_files,
                                            dirnames=dirnames,
                                            pt_bins=pt_bins,
                                            labels=labels,
                                            flav=this_flav,
                                            output_filename="%s/%s_flav_fraction_compare_jet1.%s" % (plot_dir, this_flav, OUTPUT_FMT),
                                            var_prepend=var_prepend,
                                            which_jet="1",
                                            xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(var_prepend))

        # Do for various n partons exclusively
        for n_partons in ["1", "2", "3", "4"]:
            qgf.compare_flavour_fractions_vs_pt(input_files=input_files,
                                                dirnames=dirnames,
                                                pt_bins=pt_bins,
                                                labels=labels,
                                                flav=this_flav,
                                                output_filename="%s/%s_flav_fraction_compare_jet1_npartons_%s.%s" % (plot_dir, this_flav, n_partons, OUTPUT_FMT),
                                                n_partons=n_partons,
                                                var_prepend=var_prepend,
                                                which_jet="1",
                                                title=n_partons+" outgoing parton"+("s" if n_partons != "1" else ""),
                                                xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(var_prepend))

        # Compare MG+Pythia vs Herwig++ for each region, flav
        for dname, label in zip(dirnames, labels):
            input_files = [os.path.join(root_dir, qgc.DY_FILENAME), os.path.join(root_dir, qgc.DY_HERWIG_FILENAME)]
            if 'dijet' in dname.lower():
                input_files = [os.path.join(root_dir, qgc.QCD_FILENAME), os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME)]

            labels = [MG_SAMPLE, HPP_SAMPLE]

            short_name = "dj_central"
            if "forward" in dname.lower():
                short_name = "dj_forward"
            if "zplusjet" in dname.lower():
                short_name = "zpj"

            qgf.compare_flavour_fractions_vs_pt(input_files=input_files,
                                                dirnames=[dname for i in input_files],
                                                pt_bins=pt_bins,
                                                labels=labels,
                                                flav=this_flav,
                                                output_filename="%s/%s_%s_flav_fraction_compare_mgpythia_vs_herwigpp_jet1.%s" % (plot_dir, this_flav, short_name, OUTPUT_FMT),
                                                var_prepend=var_prepend,
                                                which_jet="1",
                                                title=label,
                                                xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(var_prepend))

            # Do for various n partons exclusively
            for n_partons in ["1", "2", "3", "4"]:
                qgf.compare_flavour_fractions_vs_pt(input_files=input_files,
                                                    dirnames=[dname for i in input_files],
                                                    pt_bins=pt_bins,
                                                    labels=labels,
                                                    flav=this_flav,
                                                    output_filename="%s/%s_%s_flav_fraction_compare_mgpythia_vs_herwigpp_jet1_npartons_%s.%s" % (plot_dir, this_flav, short_name, n_partons, OUTPUT_FMT),
                                                    n_partons=n_partons,
                                                    var_prepend=var_prepend,
                                                    which_jet="1",
                                                    title=label + "\n" + n_partons+" outgoing parton"+("s" if n_partons != "1" else ""),
                                                    xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(var_prepend))


def do_flavour_fraction_input_comparison_plots(root_dirs,
                                               labels,
                                               plot_dir="flav_fractions_comparison",
                                               zpj_dirname="ZPlusJets_QG",
                                               dj_cen_dirname="Dijet_QG_central_tighter",
                                               dj_fwd_dirname="Dijet_QG_forward_tighter",
                                               var_prepend=""):
    """Do plots comparing several input dirs """
    # pt_bins = qgc.PT_BINS_INC_UFLOW
    pt_bins = qgc.PT_BINS_ZPJ
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
                                            title=qgc.Dijet_CEN_LABEL,
                                            xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(var_prepend))
        # this_flav = "1-g"
        # dirnames = [dj_cen_dirname]*len(root_dirs)
        # # Compare non-gluon fractions across samples/selections
        # qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.QCD_FILENAME) for rd in root_dirs],
        #                                     dirnames=dirnames,
        #                                     pt_bins=pt_bins,
        #                                     labels=labels,
        #                                     flav=this_flav,
        #                                     output_filename="%s/dj_q_flav_fraction_compare_central_jet.%s" % (plot_dir, OUTPUT_FMT),
        #                                     var_prepend=var_prepend,
        #                                     which_jet="1",
        #                                     xtitle="p_{T}^{jet 1} [GeV]")
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
                                            title=qgc.Dijet_FWD_LABEL,
                                            xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(var_prepend))
        # this_flav = "1-g"
        # dirnames = [dj_fwd_dirname]*len(root_dirs)
        # # Compare non-gluon fractions across samples/selections
        # qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.QCD_FILENAME) for rd in root_dirs],
        #                                     dirnames=dirnames,
        #                                     pt_bins=pt_bins,
        #                                     labels=labels,
        #                                     flav=this_flav,
        #                                     output_filename="%s/dj_q_flav_fraction_compare_forward_jet.%s" % (plot_dir, OUTPUT_FMT),
        #                                     var_prepend=var_prepend,
        #                                     which_jet="1",
        #                                     xtitle="p_{T}^{jet 1} [GeV]")

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
                                            title=qgc.ZpJ_LABEL,
                                            xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(var_prepend),
                                            var_prepend=var_prepend)
        # this_flav = "1-g"
        # dirnames = [zpj_dirname] * len(root_dirs)
        # # Compare quark fractions across samples/selections
        # qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.DY_FILENAME) for rd in root_dirs],
        #                                     dirnames=dirnames,
        #                                     pt_bins=pt_bins,
        #                                     labels=labels,
        #                                     flav=this_flav,
        #                                     output_filename="%s/zpj_q_flav_fraction_compare.%s" % (plot_dir, OUTPUT_FMT),
        #                                     var_prepend=var_prepend)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--dj",
                        help="Do dijet plots",
                        action="store_true")
    parser.add_argument("--zpj",
                        help="Do z + jets plots",
                        action="store_true")
    parser.add_argument("--gen",
                        help="Use genjet flavour",
                        action="store_true")
    parser.add_argument("--dir",
                        help="Directory to get plot from. Can be used multiple times",
                        action="append")
    parser.add_argument("--dirLabel",
                        help="Label to be given for dir. Must be used in conjunction with --dir, once per entry.",
                        action="append")
    args = parser.parse_args()
    print(args)

    # One set of plots per input
    for adir in args.dir:
        do_all_flavour_fraction_plots(adir,
                                      plot_dir=os.path.join(adir, "flav_fractions%s" % ("_gen" if args.gen else "")),
                                      var_prepend="gen" if args.gen else "",
                                      # zpj_dirname=None,
                                      zpj_dirname="ZPlusJets_QG" if args.zpj else None,
                                      dj_cen_dirname="Dijet_QG_central_tighter" if args.dj else None,
                                      dj_fwd_dirname="Dijet_QG_forward_tighter" if args.dj else None)

    if len(args.dir) > 1:
        # Now comparison across all inputs
        do_flavour_fraction_input_comparison_plots(args.dir,
                                                   plot_dir=os.path.join(args.dir[0], "flav_fractions_comparison%s" % ("_gen" if args.gen else "")),
                                                   labels=args.dirlabel,
                                                   var_prepend="gen" if args.gen else "",
                                                   dj_cen_dirname="Dijet_QG_central_tighter" if args.dj else None,
                                                   dj_fwd_dirname="Dijet_QG_forward_tighter" if args.dj else None,
                                                   zpj_dirname="ZPlusJets_QG" if args.zpj else None)
