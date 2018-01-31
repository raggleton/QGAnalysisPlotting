#!/usr/bin/env python

"""Print plots comparing MG+Pythia with plain Pythia."""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
from itertools import product
import math
from array import array

# My stuff
import qg_common as qgc
import qg_general_plots as qgg
import qg_flavour_plots as qgf
import qg_delta_plots as qgd
from comparator import grab_obj, Contribution, Plot
import common_utils as cu

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
    return
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

    # Do 2D plots for diff ptMin folders
    jet_flavs = ["", "g", "q"]
    ptmin_vals = [50, 100, 200, 400, 800]
    ptmin_vals = [100, 800]
    # first the pt vs ones:
    plotnames = ["genjet_pt_vs_constituent_pt", "genjet_pt_vs_constituent_zi", "genjet_pt_vs_constituent_deta", "genjet_pt_vs_constituent_dphi", "genjet_pt_vs_constituent_dr"]
    for plotname, ptmin, renorm, logz, flav in product(plotnames, ptmin_vals, ['Y', None], [True, False], jet_flavs):
        rebin = [1, 1]
        ylim = [ptmin, 2000]
        dirname = "Dijet_genjet_ptMin_%d" % ptmin
        this_hist = os.path.join(dirname, flav+plotname)
        log_append = "_logZ" if logz else ""
        qgg.do_2D_plot(grab_obj(os.path.join(mgpythia_dir, qgc.QCD_FILENAME), this_hist),
                               output_filename="%s/constituent_plots/ptMin_%d/%s_mgpythia_norm%s%s.%s" % (plot_dir, ptmin, flav+plotname, renorm, log_append, OUTPUT_FMT),
                               renorm_axis=renorm, title="MG+Pythia", rebin=rebin, recolour=True, ylim=ylim, logz=logz)
        qgg.do_2D_plot(grab_obj(os.path.join(pythia_only_dir, qgc.QCD_FILENAME), this_hist),
                               output_filename="%s/constituent_plots/ptMin_%d/%s_pythiaOnly_norm%s%s.%s" % (plot_dir, ptmin, flav+plotname, renorm, log_append, OUTPUT_FMT),
                               renorm_axis=renorm, title="Pythia only", rebin=rebin, recolour=True, ylim=ylim, logz=logz)

    # Now LHA vs X ones
    plotnames = ["genjet_LHA_vs_zi", "genjet_LHA_vs_thetai"]
    for plotname, ptmin, renorm, logz, flav in product(plotnames, ptmin_vals, ['X', 'Y', None], [True, False], jet_flavs):
        rebin = [4, 4]
        xlim = [0, 0.5] if ("zi" in plotname and renorm != "X") else None
        ylim = None
        dirname = "Dijet_genjet_ptMin_%d" % ptmin
        this_hist = os.path.join(dirname, flav+plotname)
        log_append = "_logZ" if logz else ""
        qgg.do_2D_plot(grab_obj(os.path.join(mgpythia_dir, qgc.QCD_FILENAME), this_hist),
                               output_filename="%s/constituent_plots/ptMin_%d/%s_mgpythia_norm%s%s.%s" % (plot_dir, ptmin, flav+plotname, renorm, log_append, OUTPUT_FMT),
                               renorm_axis=renorm, title="MG+Pythia, %s-jet, ptMin = %d GeV" % (flav, ptmin), rebin=rebin, recolour=True, xlim=xlim, ylim=ylim, logz=logz)
        qgg.do_2D_plot(grab_obj(os.path.join(pythia_only_dir, qgc.QCD_FILENAME), this_hist),
                               output_filename="%s/constituent_plots/ptMin_%d/%s_pythiaOnly_norm%s%s.%s" % (plot_dir, ptmin, flav+plotname, renorm, log_append, OUTPUT_FMT),
                               renorm_axis=renorm, title="Pythia only, %s-jet, ptMin = %d GeV" % (flav, ptmin), rebin=rebin, recolour=True, xlim=xlim, ylim=ylim, logz=logz)


    def do_2D_plot_with_contours(obj, contour_obj, output_filename, renorm_axis=None, title=None, rebin=None, recolour=True, xlim=None, ylim=None, logx=False, logy=False, logz=False):
        """Like normal 2D plotter but contour_obj will be plotted using contours"""
        if rebin:
            obj.Rebin2D(*rebin)
        if renorm_axis:
            obj_renorm = cu.make_normalised_TH2(obj, renorm_axis, recolour)
        else:
            obj_renorm = obj
        if title:
            obj_renorm.SetTitle(title)
        canvas = ROOT.TCanvas(ROOT.TUUID().AsString(), "", 800, 800)
        canvas.SetTicks(1, 1)
        canvas.SetLeftMargin(0.13)
        canvas.SetBottomMargin(0.11)
        if logx:
            canvas.SetLogx(1)
        if logy:
            canvas.SetLogy(1)
        if logz:
            canvas.SetLogz(1)
        obj_renorm.Draw("COLZ")
        if xlim is not None:
            obj_renorm.GetXaxis().SetRangeUser(*xlim)
        if ylim is not None:
            obj_renorm.GetYaxis().SetRangeUser(*ylim)
        obj_renorm.GetYaxis().SetTitleOffset(1.7)
        obj_renorm.GetXaxis().SetTitleOffset(1.2)
        if contour_obj:
            contour_obj.Draw("CONT3 SAME")
        output_filename = os.path.abspath(output_filename)
        odir = os.path.dirname(output_filename)
        if not os.path.isdir(odir):
            os.makedirs(odir)
        canvas.SaveAs(output_filename)

    #  Do 2D plots for premade jet pt cut objs
    jet_pt_vals = [(100, 200), (800, 1000)]
    for (jet_pt_min, jet_pt_max), renorm, logz in product(jet_pt_vals, ['X', 'Y', None][2:], [True, False]):
        rebin = [1, 1]
        xlim = [0, 1]
        ylim = [0, 0.05]
        dirname = "Dijet_genjet"
        jet_type = "g"
        plotname = jet_type+"genjet_constituent_zi_vs_constituent_thetai_pt%dto%d" % (jet_pt_min, jet_pt_max)
        this_hist = os.path.join(dirname, plotname)
        log_append = "_logZ" if logz else ""
        contour_hist = grab_obj(os.path.join(mgpythia_dir, qgc.QCD_FILENAME), this_hist).Clone("contours")
        for xbin in range(1, contour_hist.GetNbinsX()+1):
            for ybin in range(1, contour_hist.GetNbinsY()+1):
                xval = 0.5*(contour_hist.GetXaxis().GetBinCenter(xbin)+contour_hist.GetXaxis().GetBinCenter(xbin+1))
                yval = 0.5*(contour_hist.GetYaxis().GetBinCenter(ybin)+contour_hist.GetYaxis().GetBinCenter(ybin+1))
                contour_hist.SetBinContent(xbin, ybin, math.sqrt(xval/0.4) * yval)
        levels = array('d', [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5])
        contour_hist.SetContour(len(levels), levels)
        contour_hist.SetLineColor(ROOT.kRed)
        do_2D_plot_with_contours(grab_obj(os.path.join(mgpythia_dir, qgc.QCD_FILENAME), this_hist),
                                 contour_hist,
                                 output_filename="%s/constituent_plots/%s_mgpythia_norm%s%s.%s" % (plot_dir, plotname, renorm, log_append, OUTPUT_FMT),
                                 renorm_axis=renorm, title="MG+Pythia, %d < p_{T}^{j} < %d GeV" % (jet_pt_min, jet_pt_max),
                                 rebin=rebin, recolour=True, xlim=xlim, ylim=ylim, logz=logz)
        do_2D_plot_with_contours(grab_obj(os.path.join(pythia_only_dir, qgc.QCD_FILENAME), this_hist),
                                 contour_hist,
                                 output_filename="%s/constituent_plots/%s_pythiaOnly_norm%s%s.%s" % (plot_dir, plotname, renorm, log_append, OUTPUT_FMT),
                                 renorm_axis=renorm, title="Pythia only, %d < p_{T}^{j} < %d GeV" % (jet_pt_min, jet_pt_max),
                                 rebin=rebin, recolour=True, xlim=xlim, ylim=ylim, logz=logz)


    # Do 1D projection plots
    for (jet_pt_min, jet_pt_max) in jet_pt_vals:
        lw = 2
        dirname = "Dijet_genjet"
        jet_type = "g"
        plotname = jet_type + "genjet_constituent_zi_vs_constituent_thetai_pt%dto%d" % (jet_pt_min, jet_pt_max)
        this_hist = os.path.join(dirname, plotname)
        
        # Do projection plots for various zi bins
        cut_zi_bins = [(0, 0.001), (0, 0.005), (0, 0.01), (0, 0.05), (0.05, 0.1), (0, 0.1), (0.1, 0.2), (0.3, 0.4), (0.4, 0.5), (0.8, 1.0), (0, 1)]

        for i, (start_val, end_val) in enumerate(cut_zi_bins):
            entries_normal = []    

            h2d_mgpythia = grab_obj(os.path.join(mgpythia_dir, qgc.QCD_FILENAME), this_hist)
            mgpythia_kwargs = dict(line_color=qgc.QCD_COLOUR, fill_color=qgc.QCD_COLOUR, label="MG+Pythia", line_width=lw)
            entries_normal.append((qgg.get_projection_plot(h2d_mgpythia, start_val, end_val, 'y'), mgpythia_kwargs))

            h2d_pythiaonly = grab_obj(os.path.join(pythia_only_dir, qgc.QCD_FILENAME), this_hist)
            pythiaonly_kwargs = dict(line_color=ROOT.kRed, fill_color=ROOT.kRed, label="Pythia only", line_width=lw)
            entries_normal.append((qgg.get_projection_plot(h2d_pythiaonly, start_val, end_val, 'y'), pythiaonly_kwargs))

            ofilename = jet_type+"genjet_constituent_thetai_binned_zi%gto%g_jetPt%dto%d" % (start_val, end_val, jet_pt_min, jet_pt_max)
            rebin = 5 if end_val <= 0.1 else 2
            xlim = None if end_val <= 0.1 else [0, 0.6]
            if end_val >= 1:
                rebin = 4
                xlim = [0, 1.2]
            qgg.do_comparison_plot(entries_normal, "%s/constituent_plots/%s.%s" % (plot_dir, ofilename, OUTPUT_FMT),
                                   rebin=rebin, title="%g < z_{i} < %g, %d < p_{T}^{jet} < %d GeV" % (start_val, end_val, jet_pt_min, jet_pt_max),
                                   xtitle="#theta_{i}",
                                   xlim=xlim, ylim=None, 
                                   subplot_type="ratio", subplot_title="#splitline{Ratio wrt}{MG+Pythia}")

        # Do projections hists for various theta_i bins
        cut_thetai_bins = [(0, 0.01), (0.01, 0.02), (0.02, 0.05), (0.05, 0.10), (0.1, 0.15), (0.15, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1), (1, 1.5), (0, 2)][-1:]
        deltas = []
        delta_components = []
        colours = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2, ROOT.kOrange-3, ROOT.kMagenta, ROOT.kAzure+1]

        for i, (start_val, end_val) in enumerate(cut_thetai_bins[::]):
            entries_normal = []    

            h2d_mgpythia = grab_obj(os.path.join(mgpythia_dir, qgc.QCD_FILENAME), this_hist)
            mgpythia_kwargs = dict(line_color=qgc.QCD_COLOUR, fill_color=qgc.QCD_COLOUR, label="MG+Pythia", line_width=lw)
            h_mgpythia = qgg.get_projection_plot(h2d_mgpythia, start_val, end_val, 'x')
            h_mgpythia.Rebin(4)
            entries_normal.append((h_mgpythia.Clone(ROOT.TUUID().AsString()), mgpythia_kwargs))

            h2d_pythiaonly = grab_obj(os.path.join(pythia_only_dir, qgc.QCD_FILENAME), this_hist)
            pythiaonly_kwargs = dict(line_color=ROOT.kRed, fill_color=ROOT.kRed, label="Pythia only", line_width=lw)
            h_pythiaonly = qgg.get_projection_plot(h2d_pythiaonly, start_val, end_val, 'x')
            h_pythiaonly.Rebin(4)
            entries_normal.append((h_pythiaonly.Clone(ROOT.TUUID().AsString()), pythiaonly_kwargs))

            # rebin = 5 if end_val <= 0.1 else 2
            rebin = None
            xlim = [0, 0.5] if end_val <= 0.1 else [0, 0.2]
            if end_val == 2:
                xlim = [0, 0.2]
            ofilename = jet_type+"genjet_constituent_zi_binned_thetai%gto%g_jetPt%dto%d" % (start_val, end_val, jet_pt_min, jet_pt_max)
            qgg.do_comparison_plot(entries_normal, "%s/constituent_plots/%s.%s" % (plot_dir, ofilename, OUTPUT_FMT),
                                   rebin=rebin, title="%g < #theta_{i} < %g, %d < p_{T}^{jet} < %d GeV" % (start_val, end_val, jet_pt_min, jet_pt_max),
                                   xtitle="z_{i}", xlim=xlim,
                                   subplot_type="ratio", subplot_title="#splitline{Ratio wrt}{MG+Pythia}")
            
            ddelta_hist = qgd.get_ddelta_plot(h_mgpythia, h_pythiaonly)

            this_colour = colours[i%len(colours)]
            line_style = 1 + i // len(colours)
            c = Contribution(ddelta_hist, label="%g-%g" % (start_val, end_val), line_width=1, line_style=line_style,
                             marker_color=this_colour, line_color=this_colour, fill_color=this_colour)
            delta_components.append(c)

            deltas.append(qgd.calculate_delta(ddelta_hist))

        p = Plot(delta_components, what="hist", ytitle="p.d.f", xlim=[0, 0.2], title="%d < p_{T}^{jet} < %d GeV" % (jet_pt_min, jet_pt_max))
        p.plot("NOSTACK HISTE")
        p.save("%s/constituent_plots/ddelta_thetai_jetPt%dto%d.%s" % (plot_dir, jet_pt_min, jet_pt_max, OUTPUT_FMT))

        gr = qgd.construct_deltas_graph(deltas[::])
        c = Contribution(gr, label="BLAH", marker_style=0)
        graph_contribs = []
        graph_contribs.append(c)
        bin_labels = ["{:g} - {:g}".format(*b) for b in cut_thetai_bins]
        qgd.do_deltas_plot(graph_contribs, 
                           "%s/constituent_plots/delta_thetai_jetPt%dto%d.%s" % (plot_dir, jet_pt_min, jet_pt_max, OUTPUT_FMT),
                           bin_labels=bin_labels, 
                           title="%d < p_{T}^{jet} < %d GeV" % (jet_pt_min, jet_pt_max), 
                           xtitle="#theta_{i}")


        # do box n whisker plots
        # canv = ROOT.TCanvas(ROOT.TUUID().AsString(), "", 800, 600)
        # canv.SetTicks(1, 1)
        # leg = ROOT.TLegend(0.6, 0.6, 0.88, 0.88)
        # x_range = [0, 0.4]
        # y_range = [0, 120]
        # rebin = 4
        # bar_offset = 0.2
        # bar_width = 0.2
        # title = "%d < p_{T}^{jet} < %d GeV" % (jet_pt_min, jet_pt_max)

        # h2d_mgpythia = grab_obj(os.path.join(mgpythia_dir, qgc.QCD_FILENAME), this_hist)
        # h2d_mgpythia.RebinX(rebin)
        # h2d_mgpythia.SetLineColor(qgc.QCD_COLOUR)
        # h2d_mgpythia.SetBarOffset(-bar_offset)
        # h2d_mgpythia.SetBarWidth(bar_width)
        # h2d_mgpythia.GetXaxis().SetRangeUser(*x_range)
        # h2d_mgpythia.SetTitle(title)
        # leg.AddEntry(h2d_mgpythia, "MG+Pythia",  "L")

        # h2d_pythiaonly = grab_obj(os.path.join(pythia_only_dir, qgc.QCD_FILENAME), this_hist)
        # h2d_pythiaonly.RebinX(rebin)
        # h2d_pythiaonly.SetLineColor(ROOT.kRed)
        # h2d_pythiaonly.SetBarOffset(bar_offset)
        # h2d_pythiaonly.SetBarWidth(bar_width)
        # h2d_pythiaonly.GetXaxis().SetRangeUser(*x_range)
        # h2d_pythiaonly.SetTitle(title)
        # leg.AddEntry(h2d_pythiaonly, "Pythia only",  "L")

        # draw_opt = "candlex(00000311)"
        # # draw_opt = "candlex"
        # h2d_mgpythia.Draw(draw_opt)
        # h2d_pythiaonly.Draw(draw_opt + " same")
        # # h2d_mgpythia.GetYaxis().SetRangeUser(*y_range)
        # # h2d_pythiaonly.GetYaxis().SetRangeUser(*y_range)
        # h2d_mgpythia.GetYaxis().SetLimits(y_range[0], y_range[1])
        # h2d_pythiaonly.GetYaxis().SetLimits(y_range[0], y_range[1])
        # # ROOT.gPad.GetYaxis().SetRangeUser(*y_range)
        # leg.Draw()
        # # ROOT.gPad.BuildLegend(0.6, 0.6, 0.88, 0.88, "", "L")

        # ofilename = "genjet_constituent_pt_dr_box_pt%dto%d" % (jet_pt_min, jet_pt_max)
        # canv.SaveAs("%s/constituent_plots/%s.%s" % (plot_dir, ofilename, OUTPUT_FMT))

        # continue


if __name__ == '__main__':

    ALGOS = ["ak4", "ak8"][0:1]
    PUS = ["chs", "puppi"][1:]

    for algo, pu, in product(ALGOS, PUS):
        SETUP = algo + pu

        MGPYTHIA_DIR = "workdir_%s_mgpythia" % SETUP
        PYTHIA_ONLY_DIR = "workdir_%s_pythiaOnlyFlat_reweight" % SETUP
        
        MGPYTHIA_DIR = "workdir_%s_mgpythia_newFlav_withLeptonOverlapVeto" % SETUP
        PYTHIA_ONLY_DIR = "workdir_%s_pythiaOnlyFlat_newFlav_withLeptonOverlapVeto" % SETUP
        
        PLOT_DIR = PYTHIA_ONLY_DIR

        do_pythia_comparison_distribution_plots(MGPYTHIA_DIR, PYTHIA_ONLY_DIR, PLOT_DIR)
        do_pythia_comparison_flav_fractions_plots(MGPYTHIA_DIR, PYTHIA_ONLY_DIR, PLOT_DIR, algo, pu)

    # SETUP = "ak4"+"puppi"
    # MGPYTHIA_DIR = "workdir_%s_mgpythia" % SETUP
    # PYTHIA_ONLY_DIR = "workdir_%s_pythiaOnlyFlat_reweight" % SETUP
    # PLOT_DIR = PYTHIA_ONLY_DIR
    # do_component_plots(MGPYTHIA_DIR, PYTHIA_ONLY_DIR, PLOT_DIR)
