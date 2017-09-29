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
from comparator import grab_obj

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
        {"root_dir": mgpythia_dir, 'label': ", Madgraph+Pythia", "style": {'line_style': 1}},
        {"root_dir": pythia_only_dir, 'label': ", Pythia only", "style": {'line_style': 1, 'line_color': ROOT.kRed, 'fill_color': ROOT.kRed, }}
    ]
    qgg.do_all_exclusive_plots_comparison(sources, var_list=qgc.COMMON_VARS_WITH_FLAV,
                                          plot_dir=os.path.join(plot_dir, "mg_pythia_vs_pythia_only"),
                                          zpj_dirname="",
                                          subplot_type="ratio", subplot_title="#splitline{Ratio wrt}{MG+Pythia}",
                                          do_flav_tagged=False,
                                          pt_bins=qgc.THEORY_PT_BINS)
    qgg.do_all_exclusive_plots_comparison(sources, var_list=qgc.COMMON_VARS_WITH_FLAV,
                                          plot_dir=os.path.join(plot_dir, "mg_pythia_vs_pythia_only"),
                                          zpj_dirname="",
                                          subplot_type="ratio", subplot_title="#splitline{Ratio wrt}{MG+Pythia}",
                                          do_flav_tagged=True,
                                          pt_bins=qgc.THEORY_PT_BINS)

    # gen
    sources = [
        {"root_dir": mgpythia_dir, 'label': ", Madgraph+Pythia", "style": {'line_style': 1}},
        {"root_dir": pythia_only_dir, 'label': ", Pythia only", "style": {'line_style': 1, 'line_color': ROOT.kRed, 'fill_color': ROOT.kRed}}
    ]
    qgg.do_all_exclusive_plots_comparison(sources, var_list=qgc.COMMON_VARS_WITH_FLAV[:-1], var_prepend="gen",
                                          plot_dir=os.path.join(plot_dir, "mg_pythia_vs_pythia_only_gen"),
                                          dj_dirname=qgc.DJ_GENJET_RDIR, zpj_dirname="",
                                          subplot_type="ratio", subplot_title="#splitline{Ratio wrt}{MG+Pythia}",
                                          do_flav_tagged=False,
                                          pt_bins=qgc.THEORY_PT_BINS)
    qgg.do_all_exclusive_plots_comparison(sources, var_list=qgc.COMMON_VARS_WITH_FLAV[:-1], var_prepend="gen",
                                          plot_dir=os.path.join(plot_dir, "mg_pythia_vs_pythia_only_gen"),
                                          dj_dirname=qgc.DJ_GENJET_RDIR, zpj_dirname="",
                                          subplot_type="ratio", subplot_title="#splitline{Ratio wrt}{MG+Pythia}",
                                          do_flav_tagged=True,
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


def do_component_plots(mgpythia_dir, pythia_only_dir, plot_dir):

    ptmin_vals = [50, 100, 200, 400, 800]
    plotnames = ["genjet_pt_vs_constituent_pt", "genjet_pt_vs_constituent_deta", "genjet_pt_vs_constituent_dphi", "genjet_pt_vs_constituent_dr"]
    for plotname, ptmin, renorm, logz in product(plotnames, ptmin_vals, ['Y', None], [True, False]):
        rebin = [1, 1]
        ylim = [ptmin, 2000]
        dirname = "Dijet_genjet_ptMin_%d" % ptmin
        this_hist = os.path.join(dirname, plotname)
        log_append = "_logZ" if logz else ""
        qgg.do_2D_plot(grab_obj(os.path.join(mgpythia_dir, qgc.QCD_FILENAME), this_hist),
                               output_filename="%s/constituent_plots/ptMin_%d/%s_mgpythia_norm%s%s.%s" % (plot_dir, ptmin, plotname, renorm, log_append, OUTPUT_FMT),
                               renorm_axis=renorm, title="MG+Pythia", rebin=rebin, recolour=True, ylim=ylim, logz=logz)
        qgg.do_2D_plot(grab_obj(os.path.join(pythia_only_dir, qgc.QCD_FILENAME), this_hist),
                               output_filename="%s/constituent_plots/ptMin_%d/%s_pythiaOnly_norm%s%s.%s" % (plot_dir, ptmin, plotname, renorm, log_append, OUTPUT_FMT),
                               renorm_axis=renorm, title="Pythia only", rebin=rebin, recolour=True, ylim=ylim, logz=logz)

    jet_pt_vals = [(100, 200), (800, 1000)]
    for (jet_pt_min, jet_pt_max), renorm, logz in product(jet_pt_vals, ['X', 'Y', None], [True, False]):
        rebin = [1, 1]
        ylim = [0, jet_pt_max]
        dirname = "Dijet_genjet"
        plotname = "genjet_constituent_pt_vs_constituent_dr_pt%dto%d" % (jet_pt_min, jet_pt_max)
        this_hist = os.path.join(dirname, plotname)
        log_append = "_logZ" if logz else ""
        qgg.do_2D_plot(grab_obj(os.path.join(mgpythia_dir, qgc.QCD_FILENAME), this_hist),
                               output_filename="%s/constituent_plots/%s_mgpythia_norm%s%s.%s" % (plot_dir, plotname, renorm, log_append, OUTPUT_FMT),
                               renorm_axis=renorm, title="MG+Pythia, %d < p_{T}^{j} < %d GeV" % (jet_pt_min, jet_pt_max), rebin=rebin, recolour=True, ylim=ylim, logz=logz)
        qgg.do_2D_plot(grab_obj(os.path.join(pythia_only_dir, qgc.QCD_FILENAME), this_hist),
                               output_filename="%s/constituent_plots/%s_pythiaOnly_norm%s%s.%s" % (plot_dir, plotname, renorm, log_append, OUTPUT_FMT),
                               renorm_axis=renorm, title="Pythia only, %d < p_{T}^{j} < %d GeV" % (jet_pt_min, jet_pt_max), rebin=rebin, recolour=True, ylim=ylim, logz=logz)

    # Do projection plots for bins & compare
    for (jet_pt_min, jet_pt_max) in jet_pt_vals:
        lw = 2
        dirname = "Dijet_genjet"
        plotname = "genjet_constituent_pt_vs_constituent_dr_pt%dto%d" % (jet_pt_min, jet_pt_max)
        this_hist = os.path.join(dirname, plotname)
        cut_pt_bins = [(0, 2), (5, 10), (20, 25), (50, 60), (100, 120), (500, 700)]

        for i, (start_val, end_val) in enumerate(cut_pt_bins):
            entries_normal = []    

            h2d_mgpythia = grab_obj(os.path.join(mgpythia_dir, qgc.QCD_FILENAME), this_hist)
            mgpythia_kwargs = dict(line_color=qgc.QCD_COLOUR, fill_color=qgc.QCD_COLOUR, label="MG+Pythia", line_width=lw)
            entries_normal.append((qgg.get_projection_plot(h2d_mgpythia, start_val, end_val, 'y'), mgpythia_kwargs))

            h2d_pythiaonly = grab_obj(os.path.join(pythia_only_dir, qgc.QCD_FILENAME), this_hist)
            pythiaonly_kwargs = dict(line_color=ROOT.kRed, fill_color=ROOT.kRed, label="Pythia only", line_width=lw)
            entries_normal.append((qgg.get_projection_plot(h2d_pythiaonly, start_val, end_val, 'y'), pythiaonly_kwargs))

            rebin = 4
            ofilename = "genjet_constituent_dr_binned_pt%dto%d_jetPt%dto%d" % (start_val, end_val, jet_pt_min, jet_pt_max)
            qgg.do_comparison_plot(entries_normal, "%s/constituent_plots/%s.%s" % (plot_dir, ofilename, OUTPUT_FMT),
                                   rebin=rebin, title="%d < p_{T}^{const.} < %d GeV, %d < p_{T}^{jet} < %d GeV" % (start_val, end_val, jet_pt_min, jet_pt_max),
                                   xtitle="#DeltaR(const., jet axis)",
                                   xlim=[0, 0.6], ylim=None, 
                                   subplot_type="ratio", subplot_title="#splitline{Ratio wrt}{MG+Pythia}")

        # do box n whisker plots
        canv = ROOT.TCanvas(ROOT.TUUID().AsString(), "", 800, 600)
        canv.SetTicks(1, 1)
        leg = ROOT.TLegend(0.6, 0.6, 0.88, 0.88)
        x_range = [0, 0.4]
        y_range = [0, 120]
        rebin = 4
        bar_offset = 0.2
        bar_width = 0.2
        title = "%d < p_{T}^{jet} < %d GeV" % (jet_pt_min, jet_pt_max)

        h2d_mgpythia = grab_obj(os.path.join(mgpythia_dir, qgc.QCD_FILENAME), this_hist)
        h2d_mgpythia.RebinX(rebin)
        h2d_mgpythia.SetLineColor(qgc.QCD_COLOUR)
        h2d_mgpythia.SetBarOffset(-bar_offset)
        h2d_mgpythia.SetBarWidth(bar_width)
        h2d_mgpythia.GetXaxis().SetRangeUser(*x_range)
        h2d_mgpythia.SetTitle(title)
        leg.AddEntry(h2d_mgpythia, "MG+Pythia",  "L")

        h2d_pythiaonly = grab_obj(os.path.join(pythia_only_dir, qgc.QCD_FILENAME), this_hist)
        h2d_pythiaonly.RebinX(rebin)
        h2d_pythiaonly.SetLineColor(ROOT.kRed)
        h2d_pythiaonly.SetBarOffset(bar_offset)
        h2d_pythiaonly.SetBarWidth(bar_width)
        h2d_pythiaonly.GetXaxis().SetRangeUser(*x_range)
        h2d_pythiaonly.SetTitle(title)
        leg.AddEntry(h2d_pythiaonly, "Pythia only",  "L")

        draw_opt = "candlex(00000311)"
        # draw_opt = "candlex"
        h2d_mgpythia.Draw(draw_opt)
        # h2d_pythiaonly.Draw(draw_opt + " same")
        h2d_mgpythia.GetYaxis().SetRangeUser(*y_range)
        # ROOT.gPad.GetYaxis().SetRangeUser(*y_range)
        # h2d_pythiaonly.GetYaxis().SetRangeUser(*y_range)
        leg.Draw()
        # ROOT.gPad.BuildLegend(0.6, 0.6, 0.88, 0.88, "", "L")

        ofilename = "genjet_constituent_pt_dr_box_pt%dto%d" % (jet_pt_min, jet_pt_max)
        canv.SaveAs("%s/constituent_plots/%s.%s" % (plot_dir, ofilename, OUTPUT_FMT))

        continue

        cut_dr_bins = [(0, 0.01), (0.01, 0.02), (0.02, 0.05), (0.05, 0.10), (0.15, 0.2), (0.3, 0.4), (0.4, 0.6)]

        for i, (start_val, end_val) in enumerate(cut_dr_bins):
            entries_normal = []    

            h2d_mgpythia = grab_obj(os.path.join(mgpythia_dir, qgc.QCD_FILENAME), this_hist)
            mgpythia_kwargs = dict(line_color=qgc.QCD_COLOUR, fill_color=qgc.QCD_COLOUR, label="MG+Pythia", line_width=lw)
            entries_normal.append((qgg.get_projection_plot(h2d_mgpythia, start_val, end_val, 'x'), mgpythia_kwargs))

            h2d_pythiaonly = grab_obj(os.path.join(pythia_only_dir, qgc.QCD_FILENAME), this_hist)
            pythiaonly_kwargs = dict(line_color=ROOT.kRed, fill_color=ROOT.kRed, label="Pythia only", line_width=lw)
            entries_normal.append((qgg.get_projection_plot(h2d_pythiaonly, start_val, end_val, 'x'), pythiaonly_kwargs))

            rebin = 4 if end_val <= 0.1 else 1
            ofilename = "genjet_constituent_pt_binned_dr%gto%g_jetPt%dto%d" % (start_val, end_val, jet_pt_min, jet_pt_max)
            xlim = [0, 100] if end_val <= 0.1 else [0, 20]
            qgg.do_comparison_plot(entries_normal, "%s/constituent_plots/%s.%s" % (plot_dir, ofilename, OUTPUT_FMT),
                                   rebin=rebin, title="%g < #DeltaR(const., jet) < %g, %d < p_{T}^{jet} < %d GeV" % (start_val, end_val, jet_pt_min, jet_pt_max),
                                   xtitle="p_{T}^{const.} [GeV]",
                                   xlim=xlim, 
                                   subplot_type="ratio", subplot_title="#splitline{Ratio wrt}{MG+Pythia}")



if __name__ == '__main__':

    ALGOS = ["ak4", "ak8"]
    PUS = ["chs", "puppi"]

    for algo, pu, in product(ALGOS, PUS):
        SETUP = algo + pu

        MGPYTHIA_DIR = "workdir_%s_mgpythia" % SETUP
        PYTHIA_ONLY_DIR = "workdir_%s_pythiaOnlyFlat_reweight" % SETUP
        PLOT_DIR = PYTHIA_ONLY_DIR

        do_pythia_comparison_distribution_plots(MGPYTHIA_DIR, PYTHIA_ONLY_DIR, PLOT_DIR)
        do_pythia_comparison_flav_fractions_plots(MGPYTHIA_DIR, PYTHIA_ONLY_DIR, PLOT_DIR, algo, pu)

    SETUP = "ak4"+"puppi"
    MGPYTHIA_DIR = "workdir_%s_mgpythia" % SETUP
    PYTHIA_ONLY_DIR = "workdir_%s_pythiaOnlyFlat_reweight" % SETUP
    PLOT_DIR = PYTHIA_ONLY_DIR
    do_component_plots(MGPYTHIA_DIR, PYTHIA_ONLY_DIR, PLOT_DIR)
