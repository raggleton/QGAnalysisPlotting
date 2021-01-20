#!/usr/bin/env python

"""Print flavour fraction plots/graphs"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
os.nice(10)
import argparse
from array import array
import numpy as np

# My stuff
from comparator import Contribution, Plot
import qg_common as qgc
import qg_general_plots as qgg
import qg_flavour_plots as qgf
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
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# Control plot output format
OUTPUT_FMT = "pdf"

MG_SAMPLE = "MG5+Pythia8"
HPP_SAMPLE = "Herwig++"



def do_flav_nparton_fraction_vs_pt(input_file,
                                   dirname,
                                   flav,
                                   n_partons,
                                   pt_bins,
                                   output_filename,
                                   var_prepend='',
                                   which_jet='1',
                                   title=''):

    tfile = cu.open_root_file(input_file)
    hist_2d_total = tfile.Get("%s/%sjet%s_flavour_vs_pt" % (dirname, var_prepend, which_jet))
    hists_2d = [tfile.Get("%s/%sjet%s_flavour_vs_pt_npartons_%s" % (dirname, var_prepend, which_jet, n_parton))
                for n_parton in n_partons]
    # pdigid is x axis, pt is y axis
    # so we want the 1D projection on the y axis,
    # for the x bin corresponding to this flav
    flav_dict = {
        'd': 1,
        'u': 2,
        's': 3,
        'c': 4,
        'b': 5,
        't': 6,
        'g': 21,
        'unknown': 0,
    }
    flav_bin = hists_2d[0].GetXaxis().FindFixBin(flav_dict[flav])
    projs = [h.ProjectionY(cu.get_unique_str(), flav_bin, flav_bin, 'e')
             for h in hists_2d]
    hist_total = hist_2d_total.ProjectionY(cu.get_unique_str(), flav_bin, flav_bin, 'e')

    # rebin hists
    # add (0, x) bin for lowest bin, since hist has that
    bins = np.array([(0., qgc.PT_BINS_ZPJ_INC_UFLOW[0][0])] + qgc.PT_BINS_ZPJ_INC_UFLOW)
    # convert to only have each edge once, not pairs of edges
    bins = np.concatenate((bins[:,0], bins[-1:,1]))
    hist_total = hist_total.Rebin(len(bins)-1, cu.get_unique_str(), bins)
    rebinned_hists = [h.Rebin(len(bins)-1, cu.get_unique_str(), bins) for h in projs]

    # divide each by total to get fraction
    for h in rebinned_hists:
        h.Divide(hist_total)

    # create Contributions for Plot
    colours = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2]
    marker = cu.Marker()
    entries = [Contribution(h, label="%s parton%s" % (n_parton, 's' if n_parton != '1' else ''),
                            line_width=1, line_color=colours[int(n_parton)-1], line_style=1,
                            fill_color=colours[int(n_parton)-1], fill_style=0,
                            marker_size=1, marker_color=colours[int(n_parton)-1], marker_style=mark,
                            leg_draw_opt=None)
               for n_parton, h, mark in zip(n_partons, rebinned_hists, marker.cycle())
               if h.GetEntries() > 0]

    plot = Plot(entries,
                what='hist',
                xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(var_prepend),
                ytitle="Fraction",
                title=title,
                xlim=(pt_bins[0][0], pt_bins[-1][1]),
                ylim=[0, 1.4],
                has_data=False)
    plot.default_canvas_size = (600, 600)
    plot.plot("HIST E NOSTACK")
    plot.set_logx(do_more_labels=True, do_exponent=False)
    plot.legend.SetX1(0.55)
    plot.legend.SetY1(0.72)
    plot.legend.SetY2(0.85)
    plot.save(output_filename)


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
    if dj_fwd_dirname and zpj_dirname:
        dirnames = [dj_cen_dirname, dj_fwd_dirname, zpj_dirname]
        for this_flav in ['g', 'u', 'd', 's', 'c', 'b'][:]:
            input_files = [os.path.join(root_dir, qgc.QCD_FILENAME) if "dijet" in d.lower()
                           else os.path.join(root_dir, qgc.DY_FILENAME)
                           for d in dirnames]

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

    dirnames, labels = [], []
    if dj_cen_dirname:
        dirnames += [dj_cen_dirname, dj_fwd_dirname]
        labels += [qgc.Dijet_CEN_LABEL, qgc.Dijet_FWD_LABEL]
    if zpj_dirname:
        dirnames += [zpj_dirname]
        labels += [qgc.ZpJ_LABEL]

    # Compare MG+Pythia vs Herwig++ for each region, flav
    for dname, label in zip(dirnames, labels):
        input_files = [os.path.join(root_dir, qgc.DY_FILENAME), os.path.join(root_dir, qgc.DY_HERWIG_FILENAME)]
        if 'dijet' in dname.lower():
            input_files = [os.path.join(root_dir, qgc.QCD_FILENAME), os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME)]

        sample_labels = [MG_SAMPLE, HPP_SAMPLE]
        sample_short_names = ['mgpythia', 'herwigpp']

        short_name = "dj_central"
        if "forward" in dname.lower():
            short_name = "dj_forward"
        if "zplusjet" in dname.lower():
            short_name = "zpj"

        for this_flav in ['g', 'u', 'd', 's', 'c', 'b'][:]:

            qgf.compare_flavour_fractions_vs_pt(input_files=input_files,
                                                dirnames=[dname for i in input_files],
                                                pt_bins=pt_bins,
                                                labels=sample_labels,
                                                flav=this_flav,
                                                output_filename="%s/%s_%s_flav_fraction_compare_mgpythia_vs_herwigpp_jet1.%s" % (plot_dir, this_flav, short_name, OUTPUT_FMT),
                                                var_prepend=var_prepend,
                                                which_jet="1",
                                                title=label,
                                                xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(var_prepend))

            # Do for various n partons exclusively
            all_n_partons = ["2", "3", "4"] if 'dijet' in dname.lower() else ["1", "2", "3", "4"]
            for n_partons in all_n_partons:
                qgf.compare_flavour_fractions_vs_pt(input_files=input_files,
                                                    dirnames=[dname for i in input_files],
                                                    pt_bins=pt_bins,
                                                    labels=sample_labels,
                                                    flav=this_flav,
                                                    output_filename="%s/%s_%s_flav_fraction_compare_mgpythia_vs_herwigpp_jet1_npartons_%s.%s" % (plot_dir, this_flav, short_name, n_partons, OUTPUT_FMT),
                                                    n_partons=n_partons,
                                                    var_prepend=var_prepend,
                                                    which_jet="1",
                                                    title=label + "\n" + n_partons+" outgoing parton"+("s" if n_partons != "1" else ""),
                                                    xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(var_prepend))

            # Do all n partons on one plot, for brevity
            qgf.compare_flavour_fractions_vs_pt(input_files=input_files,
                                                dirnames=[dname for i in input_files],
                                                pt_bins=pt_bins,
                                                labels=sample_labels,
                                                flav=this_flav,
                                                output_filename="%s/%s_%s_flav_fraction_compare_mgpythia_vs_herwigpp_jet1_npartons_all.%s" % (plot_dir, this_flav, short_name, OUTPUT_FMT),
                                                n_partons=all_n_partons,
                                                var_prepend=var_prepend,
                                                which_jet="1",
                                                title=label + "\n%s outgoing partons" % (", ".join(["%s-" % n for n in all_n_partons])),
                                                xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(var_prepend))


            # Do a plot of fractions of n partons for this flavour / input etc
            for input_file, sample_label, sample_short_name in zip(input_files, sample_labels, sample_short_names):
                do_flav_nparton_fraction_vs_pt(input_file=input_file,
                                               dirname=dname,
                                               flav=this_flav,
                                               n_partons=['2', '3', '4'] if 'dijet' in dname.lower() else ['1', '2', '3', '4'],
                                               pt_bins=pt_bins,
                                               var_prepend=var_prepend,
                                               which_jet='1',
                                               title=label + "\n" + sample_label + '\n%s flavour jets' % this_flav,
                                               output_filename="%s/%s_%s_%s_nparton_fraction.%s" % (plot_dir, this_flav, short_name, sample_short_name, OUTPUT_FMT))



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
                                            title=qgc.Dijet_CEN_LABEL + "\nMG+PYTHIA8",
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
                                            title=qgc.Dijet_FWD_LABEL + "\nMG+PYTHIA8",
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
                                            title=qgc.ZpJ_LABEL + "\nMG+PYTHIA8",
                                            xtitle="p_{T}^{%s} [GeV]" % qgf.get_jet_str(var_prepend),
                                            var_prepend=var_prepend)

        qgf.compare_flavour_fractions_vs_pt(input_files=[os.path.join(rd, qgc.DY_HERWIG_FILENAME) for rd in root_dirs],
                                            dirnames=dirnames,
                                            pt_bins=pt_bins,
                                            labels=labels,
                                            flav=this_flav,
                                            output_filename="%s/zpj_hpp_g_flav_fraction_compare.%s" % (plot_dir, OUTPUT_FMT),
                                            title=qgc.ZpJ_LABEL + "\nHERWIG++",
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
                                                   labels=args.dirLabel,
                                                   var_prepend="gen" if args.gen else "",
                                                   dj_cen_dirname="Dijet_QG_central_tighter" if args.dj else None,
                                                   dj_fwd_dirname="Dijet_QG_forward_tighter" if args.dj else None,
                                                   zpj_dirname="ZPlusJets_QG" if args.zpj else None)
