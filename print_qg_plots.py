#!/usr/bin/env python

"""Print QG plots"""

import ROOT
from comparator import Contribution, Plot, grab_obj
from TDRStyle import TDR_Style
from MyStyle import My_Style
import common_utils as cu
from uuid import uuid4
import bisect
import numpy as np
from itertools import chain
import os


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
# TDR_Style.cd()
# My_Style.cd()


DY_COLOUR = 880
QCD_COLOUR = 867

DY_ZpJ_LABEL = "DY+j, Z+jets selection"
DY_ZpJ_GEN_LABEL = "DY+j, Z+jets selection (GenJets)"
DY_ZpJ_QFLAV_LABEL = "DY+j, Z+jets selection (uds-matched)"
DY_ZpJ_GFLAV_LABEL = "DY+j, Z+jets selection (g-matched)"

DY_Dijet_LABEL = "DY+j, Dijet selection"
DY_Dijet_GEN_LABEL = "DY+j, Dijet selection (GenJets)"
DY_Dijet_QFLAV_LABEL = "DY+j, Dijet selection (uds-matched)"
DY_Dijet_GFLAV_LABEL = "DY+j, Dijet selection (g-matched)"

QCD_ZpJ_LABEL = "QCD, Z+jets selection"
QCD_ZpJ_GEN_LABEL = "QCD, Z+jets selection (GenJets)"
QCD_ZpJ_QFLAV_LABEL = "QCD, Z+jets selection (uds-matched)"
QCD_ZpJ_GFLAV_LABEL = "QCD, Z+jets selection (g-matched)"

QCD_Dijet_LABEL = "QCD, Dijet selection"
QCD_Dijet_GEN_LABEL = "QCD, Dijet selection (GenJets)"
QCD_Dijet_QFLAV_LABEL = "QCD, Dijet selection (uds-matched)"
QCD_Dijet_GFLAV_LABEL = "QCD, Dijet selection (g-matched)"

COMMON_VARS = ['jet_LHA', 'jet_pTD', 'jet_width', 'jet_thrust', 'jet_multiplicity', 'jet_flavour', "jet_genParton_flavour"]

CHS_DIR = "workdir_ak4chs"
PUPPI_DIR = "workdir_ak4puppi"

ROOT_DIR = CHS_DIR
# ROOT_DIR = PUPPI_DIR

TITLE_STR = "[%s]" % ROOT_DIR.replace("workdir_", "")


PT_BINS = [(80, 100), (100, 200), (400, 500), (1000, 2000)]


def do_comparison_plot(entries, output_filename, rebin=1, **plot_kwargs):
    conts = [Contribution(ent[0], normalise_hist=True, fill_style=0, rebin_hist=rebin, **ent[1]) for ent in entries]
    p = Plot(conts, what="hist", subplot=conts[0], ytitle="p.d.f", subplot_type="diff", **plot_kwargs)
    draw_opt = "NOSTACK HISTE"
    p.legend.SetX1(0.55)
    p.legend.SetY1(0.6)
    p.plot(draw_opt)
    p.save(output_filename)


def get_projection_plot(h2d, start_val, end_val):
    y_axis = h2d.GetYaxis()
    bin_edges = [y_axis.GetBinLowEdge(i) for i in range(1, y_axis.GetNbins()+1)]
    bin_start = bisect.bisect_right(bin_edges, start_val)
    bin_end = bisect.bisect_right(bin_edges, end_val)
    hproj = h2d.ProjectionX(h2d.GetName()+"_projX"+str(uuid4()), bin_start, bin_end, "eo")
    return hproj


def do_2D_plot(obj, output_filename, renorm_axis=None, title=None, rebin=None, recolour=True):
    if rebin:
        obj.Rebin2D(*rebin)
    if renorm_axis:
        obj_renorm = cu.make_normalised_TH2(obj, renorm_axis, recolour)
    else:
        obj_renorm = obj
    if title:
        obj_renorm.SetTitle(title)
    canvas = ROOT.TCanvas("canv%s" % uuid4(), "", 800, 600)
    canvas.SetTicks(1, 1)
    obj_renorm.Draw("COLZ")
    obj_renorm.GetYaxis().SetTitleOffset(1.3)
    obj_renorm.GetXaxis().SetTitleOffset(1.1)
    odir = os.path.dirname(os.path.abspath(output_filename))
    if not os.path.isdir(odir):
        os.makedirs(odir)
    canvas.SaveAs(output_filename)


def do_all_2D_plots(plot_dir="plots_2d", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG", var_list=None, var_prepend=""):
    """Do 2D distributions"""
    var_list = var_list or COMMON_VARS
    for v in var_list:
        v = "%s%s_vs_pt" % (var_prepend, v)

        rebin = [2, 4]
        if v == "jet_multiplicity_vs_pt":
            rebin = [2, 4]
        elif v == "jet_thrust_vs_pt" or "flavour" in v:
            rebin = [1, 4]

        recolour = False if "flavour" in v else True

        for rn in ['Y']:  # different renormalisation axes
            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "%s/%s" % (zpj_dirname, v)),
                       output_filename="%s/%s/dy_zpj_%s_norm%s.pdf" % (ROOT_DIR, plot_dir, v, rn),
                       renorm_axis=rn, title=DY_ZpJ_LABEL, rebin=rebin, recolour=recolour)
            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "%s/%s" % (dj_dirname, v)),
                       output_filename="%s/%s/dy_dijet_%s_norm%s.pdf" % (ROOT_DIR, plot_dir, v, rn),
                       renorm_axis=rn, title=DY_Dijet_LABEL, rebin=rebin, recolour=recolour)

            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "%s/%s" % (dj_dirname, v)),
                       output_filename="%s/%s/qcd_dijet_%s_norm%s.pdf" % (ROOT_DIR, plot_dir, v, rn),
                       renorm_axis=rn, title=QCD_Dijet_LABEL, rebin=rebin, recolour=recolour)

            if "flavour" in v:
                continue

            # Flavour matched reco
            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "%s/q%s" % (zpj_dirname, v)),
                       output_filename="%s/%s/dy_zpj_%s_norm%s_qflavMatched.pdf" % (ROOT_DIR, plot_dir, v, rn),
                       renorm_axis=rn, title=DY_ZpJ_QFLAV_LABEL, rebin=rebin)

            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "%s/g%s" % (zpj_dirname, v)),
                       output_filename="%s/%s/dy_zpj_%s_norm%s_gflavMatched.pdf" % (ROOT_DIR, plot_dir, v, rn),
                       renorm_axis=rn, title=DY_ZpJ_GFLAV_LABEL, rebin=rebin)

            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "%s/g%s" % (dj_dirname, v)),
                       output_filename="%s/%s/dy_dijet_%s_norm%s_gflavMatched.pdf" % (ROOT_DIR, plot_dir, v, rn),
                       renorm_axis=rn, title=DY_Dijet_GFLAV_LABEL, rebin=rebin)

            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "%s/q%s" % (dj_dirname, v)),
                       output_filename="%s/%s/dy_dijet_%s_norm%s_qflavMatched.pdf" % (ROOT_DIR, plot_dir, v, rn),
                       renorm_axis=rn, title=DY_Dijet_QFLAV_LABEL, rebin=rebin)

            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "%s/g%s" % (dj_dirname, v)),
                       output_filename="%s/%s/qcd_dijet_%s_norm%s_gflavMatched.pdf" % (ROOT_DIR, plot_dir, v, rn),
                       renorm_axis=rn, title=QCD_Dijet_GFLAV_LABEL, rebin=rebin)

            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "%s/q%s" % (dj_dirname, v)),
                       output_filename="%s/%s/qcd_dijet_%s_norm%s_qflavMatched.pdf" % (ROOT_DIR, plot_dir, v, rn),
                       renorm_axis=rn, title=QCD_Dijet_QFLAV_LABEL, rebin=rebin)


def do_all_exclusive_plots(plot_dir="plots_dy_vs_qcd", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG",
                           var_list=None, var_prepend="", pt_bins=None):
    """Do pt/eta/nvtx binned 1D plots"""
    var_list = var_list or COMMON_VARS
    pt_bins = pt_bins or PT_BINS

    for v in var_list:
        v = "%s%s_vs_pt" % (var_prepend, v)

        h2d_dyj = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "%s/%s" % (zpj_dirname, v))
        h2d_qcd = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "%s/%s" % (dj_dirname, v))

        if "flavour" not in v:
            h2d_dyj_q = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "%s/q%s" % (zpj_dirname, v))
            h2d_qcd_g = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "%s/g%s" % (dj_dirname, v))

        lw = 2
        dy_kwargs = dict(line_color=DY_COLOUR, fill_color=DY_COLOUR, label=DY_ZpJ_LABEL, line_width=lw)
        qcd_kwargs = dict(line_color=QCD_COLOUR, fill_color=QCD_COLOUR, label=QCD_Dijet_LABEL, line_width=lw)

        dy_kwargs_q = dict(line_color=DY_COLOUR, fill_color=DY_COLOUR, label=DY_ZpJ_QFLAV_LABEL, line_width=lw)
        qcd_kwargs_g = dict(line_color=QCD_COLOUR, fill_color=QCD_COLOUR, label=QCD_Dijet_GFLAV_LABEL, line_width=lw)

        rebin = 2
        if v == "jet_multiplicity_vs_pt":
            rebin = 4
        elif "flavour" in v or "thrust" in v:
            rebin = 1

        xlim = None
        if "thrust" in v:
            xlim = (0, 0.5)

        for (start_val, end_val) in pt_bins:
            entries = [
                (get_projection_plot(h2d_dyj, start_val, end_val), dy_kwargs),
                (get_projection_plot(h2d_qcd, start_val, end_val), qcd_kwargs),
            ]

            do_comparison_plot(entries, "%s/%s/ptBinned/%s_pt%dto%d.pdf" % (ROOT_DIR, plot_dir, v, start_val, end_val),
                               rebin=rebin, title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val), xlim=xlim)

            if "flavour" in v:
                continue

            entries = [
                (get_projection_plot(h2d_dyj_q, start_val, end_val), dy_kwargs_q),
                (get_projection_plot(h2d_qcd_g, start_val, end_val), qcd_kwargs_g),
            ]

            do_comparison_plot(entries, "%s/%s/ptBinned/%s_pt%dto%d_flavMatched.pdf" % (ROOT_DIR, plot_dir, v, start_val, end_val),
                               rebin=rebin, title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val), xlim=xlim)


def do_all_flavour_fraction_plots(var_prepend="", plot_dir="flav_fractions", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG"):
    """Do plots of flavour fraction of jets as a function of PT, for various samples"""

    def get_flavour_fractions(input_file, dirname, which="", var_prepend=""):
        h2d_flav = grab_obj(input_file, "%s/%sjet_%sflavour_vs_pt" % (dirname, var_prepend, which))

        h2d_flav.Rebin2D(1, 5)
        y_axis = h2d_flav.GetYaxis()
        pt_bins_lower, pt_bins_upper = [], []
        flav_dict = {'d': [], 'u': [], 's': [], 'c': [], 'b': [] ,'t': [], 'g': []}
        for i in range(1, y_axis.GetNbins()-1):

            d_frac = h2d_flav.GetBinContent(2, i)
            u_frac = h2d_flav.GetBinContent(3, i)
            s_frac = h2d_flav.GetBinContent(4, i)
            c_frac = h2d_flav.GetBinContent(5, i)
            b_frac = h2d_flav.GetBinContent(6, i)
            t_frac = h2d_flav.GetBinContent(7, i)
            g_frac = h2d_flav.GetBinContent(22, i)

            total = d_frac+u_frac+s_frac+c_frac+b_frac+t_frac+g_frac

            if total == 0:
                continue

            pt_bins_lower.append(y_axis.GetBinLowEdge(i))
            pt_bins_upper.append(y_axis.GetBinLowEdge(i+1))

            flav_dict['d'].append(d_frac / total)
            flav_dict['u'].append(u_frac / total)
            flav_dict['s'].append(s_frac / total)
            flav_dict['c'].append(c_frac / total)
            flav_dict['b'].append(b_frac / total)
            flav_dict['t'].append(t_frac / total)
            flav_dict['g'].append(g_frac / total)

        x_bins = [0.5 * (x1+x2) for x1,x2 in zip(pt_bins_lower[:], pt_bins_upper[:])]

        return x_bins, flav_dict

    def do_flavour_fraction_vs_pt(input_file, dirname, output_filename, title="", which="", var_prepend=""):
        """Plot flavour fractions vs PT for one input file & dirname"""
        x_bins, flav_dict = get_flavour_fractions(input_file, dirname, which, var_prepend)
        # TODO: check if empy arrays
        gr_flav_u = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['u']))
        gr_flav_d = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['d']))
        gr_flav_s = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['s']))
        gr_flav_c = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['c']))
        gr_flav_b = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['b']))
        gr_flav_g = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['g']))
        # gr_flav_uds = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array([u+d+s for u, d, s in zip(flav_dict['u'], flav_dict['d'], flav_dict['s'])]))

        plot_u = Contribution(gr_flav_u, label="u frac", line_color=ROOT.kRed, marker_color=ROOT.kRed)
        plot_d = Contribution(gr_flav_d, label="d frac", line_color=ROOT.kBlue, marker_color=ROOT.kBlue)
        plot_s = Contribution(gr_flav_s, label="s frac", line_color=ROOT.kBlack, marker_color=ROOT.kBlack)
        plot_c = Contribution(gr_flav_c, label="c frac", line_color=ROOT.kGreen-3, marker_color=ROOT.kGreen-3)
        plot_b = Contribution(gr_flav_b, label="b frac", line_color=ROOT.kOrange, marker_color=ROOT.kOrange)
        plot_g = Contribution(gr_flav_g, label="g frac", line_color=ROOT.kViolet, marker_color=ROOT.kViolet)
        # plot_uds = Contribution(gr_flav_uds, label="uds frac", line_color=ROOT.kOrange+2, marker_color=ROOT.kOrange+2)

        p_flav = Plot([plot_u, plot_d, plot_s, plot_g, plot_c, plot_b], what='graph', xtitle="p_{T}^{jet} [GeV]", title=title)
        p_flav.plot("ALP")
        p_flav.save(output_filename)

    do_flavour_fraction_vs_pt("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, dj_dirname, "%s/%s/qcd_dijet_flav_frac.pdf" % (ROOT_DIR, plot_dir), title=QCD_Dijet_LABEL, var_prepend=var_prepend)
    do_flavour_fraction_vs_pt("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, dj_dirname, "%s/%s/dy_dijet_flav_frac.pdf" % (ROOT_DIR, plot_dir), title=DY_Dijet_LABEL, var_prepend=var_prepend)
    do_flavour_fraction_vs_pt("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, zpj_dirname, "%s/%s/dy_zpj_flav_frac.pdf" % (ROOT_DIR, plot_dir), title=DY_ZpJ_LABEL, var_prepend=var_prepend)

    if var_prepend != "gen":
        do_flavour_fraction_vs_pt("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, dj_dirname, "%s/%s/qcd_dijet_genParton_flav_frac.pdf" % (ROOT_DIR, plot_dir), title=QCD_Dijet_LABEL, which="genParton_")
        do_flavour_fraction_vs_pt("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, dj_dirname, "%s/%s/dy_dijet_genParton_flav_frac.pdf" % (ROOT_DIR, plot_dir), title=DY_Dijet_LABEL, which="genParton_")
        do_flavour_fraction_vs_pt("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, zpj_dirname, "%s/%s/dy_zpj_genParton_flav_frac.pdf" % (ROOT_DIR, plot_dir), title=DY_ZpJ_LABEL, which="genParton_")

    def compare_flavour_fractions_vs_pt(input_files, dirnames, labels, flav, output_filename, title="", which="", var_prepend=""):
        """Plot a specified flavour fraction vs pT for several sources.
        Each entry in input_files, dirnames, and labels corresponds to one line"""
        info = [get_flavour_fractions(ifile, sel, which=which, var_prepend=var_prepend) for ifile, sel in zip(input_files, dirnames)]
        contribs = []
        for i, (x_bins, fdict) in enumerate(info):
            if flav in ['u', 'd', 's', 'c', 'b', 't', 'g']:
                gr = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(fdict[flav]))
            else:
                gr = ROOT.TGraph(len(x_bins), np.array(x_bins), 1.-np.array(fdict[flav.replace("1-", '')]))
            c = Contribution(gr, label="%s" % (labels[i]), line_style=i+1)
            contribs.append(c)
        ytitle = "%s flavour fraction" % flav
        p = Plot(contribs, what='graph', xtitle="p_{T}^{jet} [GeV]", ytitle=ytitle, title=title)
        p.plot("ALP")
        p.canvas.SetLogx()
        p.save(output_filename)

    input_files = [
        "%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR,
        "%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR,
        "%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR
    ]
    compare_flavour_fractions_vs_pt(input_files, [dj_dirname, dj_dirname, zpj_dirname], [QCD_Dijet_LABEL, DY_Dijet_LABEL, DY_ZpJ_LABEL], "g", "%s/%s/compare_g_frac.pdf" % (ROOT_DIR, plot_dir), title=TITLE_STR, var_prepend=var_prepend)
    if var_prepend != "gen":
        compare_flavour_fractions_vs_pt(input_files, [dj_dirname, dj_dirname, zpj_dirname], [QCD_Dijet_LABEL, DY_Dijet_LABEL, DY_ZpJ_LABEL], "g", "%s/%s/compare_g_frac_genParton.pdf" % (ROOT_DIR, plot_dir), title=TITLE_STR, which="genParton_")


def do_chs_vs_puppi_plots():
    for v in ['jet_LHA', 'jet_pTD', 'jet_width', 'jet_thrust', 'jet_multiplicity']:
        v += "_vs_pt"

        h2d_dyj_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % CHS_DIR, "ZPlusJets_QG/q%s" % v)
        h2d_qcd_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % CHS_DIR, "Dijet_QG/g%s" % v)

        h2d_dyj_puppi = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % PUPPI_DIR, "ZPlusJets_QG/q%s" % v)
        h2d_qcd_puppi = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % PUPPI_DIR, "Dijet_QG/g%s" % v)

        lw = 2
        dy_kwargs_chs = dict(line_color=DY_COLOUR, fill_color=DY_COLOUR, label=DY_ZpJ_QFLAV_LABEL + " [CHS]", line_width=lw)
        qcd_kwargs_chs = dict(line_color=QCD_COLOUR, fill_color=QCD_COLOUR, label=QCD_Dijet_GFLAV_LABEL + " [CHS]", line_width=lw)

        puppi_ls = 3
        dy_kwargs_puppi = dict(line_color=DY_COLOUR+2, fill_color=DY_COLOUR+2, label=DY_ZpJ_QFLAV_LABEL + " [PUPPI]", line_style=puppi_ls, line_width=lw)
        qcd_kwargs_puppi = dict(line_color=QCD_COLOUR-3, fill_color=QCD_COLOUR-3, label=QCD_Dijet_GFLAV_LABEL + " [PUPPI]", line_style=puppi_ls, line_width=lw)

        rebin = 2
        xlim = None
        if "thrust" in v:
            rebin = 1
            xlim = (0, 0.5)

        for (start_val, end_val) in PT_BINS:
            entries = [
                (get_projection_plot(h2d_dyj_chs, start_val, end_val), dy_kwargs_chs),
                (get_projection_plot(h2d_qcd_chs, start_val, end_val), qcd_kwargs_chs),
                (get_projection_plot(h2d_dyj_puppi, start_val, end_val), dy_kwargs_puppi),
                (get_projection_plot(h2d_qcd_puppi, start_val, end_val), qcd_kwargs_puppi),
            ]

            do_comparison_plot(entries, "ak4_chs_vs_puppi/%s_pt%dto%d_flavMatched.pdf" % (v, start_val, end_val),
                               rebin=rebin, title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val), xlim=xlim)


def do_wrong_plots(var_prepend="", plot_dir="wrong_flavs", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG", pt_bins=None):
    # var_list = var_list or COMMON_VARS
    pt_bins = pt_bins or PT_BINS
    for v in ['jet_LHA', 'jet_pTD', 'jet_width', 'jet_thrust', 'jet_multiplicity']:
        v = "%s%s_vs_pt" % (var_prepend, v)

        h2d_dyj_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % CHS_DIR, "%s/q%s" % (zpj_dirname, v))
        h2d_dyj_wrong_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % CHS_DIR, "%s/g%s" % (zpj_dirname, v))
        h2d_dyj_qcd_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % CHS_DIR, "%s/q%s" % (dj_dirname, v))
        h2d_dyj_qcd_wrong_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % CHS_DIR, "%s/g%s" % (dj_dirname, v))
        h2d_qcd_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % CHS_DIR, "%s/g%s" % (dj_dirname, v))
        h2d_qcd_wrong_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % CHS_DIR, "%s/q%s" % (dj_dirname, v))

        lw = 1
        dy_kwargs_chs = dict(line_color=DY_COLOUR, fill_color=DY_COLOUR, label=DY_ZpJ_QFLAV_LABEL, line_width=lw)
        dy_kwargs_wrong_chs = dict(line_color=DY_COLOUR+4, fill_color=DY_COLOUR+4, label=DY_ZpJ_GFLAV_LABEL, line_width=lw, line_style=1)
        dy_kwargs_qcd_chs = dict(line_color=ROOT.kGreen+2, fill_color=ROOT.kGreen+2, label=DY_Dijet_QFLAV_LABEL, line_width=lw, line_style=1)
        dy_kwargs_qcd_wrong_chs = dict(line_color=ROOT.kOrange-1, fill_color=ROOT.kOrange-1, label=DY_Dijet_GFLAV_LABEL, line_width=lw, line_style=1)
        qcd_kwargs_chs = dict(line_color=QCD_COLOUR, fill_color=QCD_COLOUR, label=QCD_Dijet_GFLAV_LABEL, line_width=lw)
        qcd_kwargs_wrong_chs = dict(line_color=ROOT.kRed, fill_color=ROOT.kRed, label=QCD_Dijet_QFLAV_LABEL, line_width=lw, line_style=1)

        rebin = 2
        xlim = None
        if "thrust" in v:
            rebin = 1
            xlim = (0, 0.5)

        for (start_val, end_val) in pt_bins:
            entries = [
                (get_projection_plot(h2d_dyj_chs, start_val, end_val), dy_kwargs_chs),
                (get_projection_plot(h2d_dyj_wrong_chs, start_val, end_val), dy_kwargs_wrong_chs),
                (get_projection_plot(h2d_dyj_qcd_chs, start_val, end_val), dy_kwargs_qcd_chs),
                (get_projection_plot(h2d_dyj_qcd_wrong_chs, start_val, end_val), dy_kwargs_qcd_wrong_chs),
                (get_projection_plot(h2d_qcd_wrong_chs, start_val, end_val), qcd_kwargs_wrong_chs),
                (get_projection_plot(h2d_qcd_chs, start_val, end_val), qcd_kwargs_chs),
            ]

            do_comparison_plot(entries, "%s/%s/%s_pt%dto%d_flavMatched.pdf" % (ROOT_DIR, plot_dir, v, start_val, end_val),
                               rebin=rebin, title="%d < p_{T}^{jet} < %d GeV %s (\"wrong\" flavours)" % (start_val, end_val, TITLE_STR), xlim=xlim)


def do_reco_plots():
    global TITLE_STR
    TITLE_STR = "[%s]" % ROOT_DIR.replace("workdir_", "")
    # do_all_2D_plots()
    # do_all_exclusive_plots()
    # do_all_flavour_fraction_plots()
    # do_chs_vs_puppi_plots()
    do_wrong_plots()


def do_gen_plots():
    global TITLE_STR
    TITLE_STR = TITLE_STR.replace("chs", " GenJet").replace("puppi", " GenJet")
    do_all_2D_plots(var_list=COMMON_VARS[:-1], var_prepend="gen", plot_dir="plots_2d_gen",
                    zpj_dirname="ZPlusJets_genjet", dj_dirname="Dijet_genjet")
    do_all_exclusive_plots(var_list=COMMON_VARS[:-1], var_prepend="gen", plot_dir="plots_dy_vs_qcd_gen",
                           zpj_dirname="ZPlusJets_genjet", dj_dirname="Dijet_genjet", pt_bins=PT_BINS+[(100, 2000)])
    do_all_flavour_fraction_plots(var_prepend="gen", plot_dir="flav_fractions_gen",
                                  zpj_dirname="ZPlusJets_genjet", dj_dirname="Dijet_genjet")
    do_wrong_plots(var_prepend="gen", plot_dir="wrong_flavs_gen",
                   zpj_dirname="ZPlusJets_genjet", dj_dirname="Dijet_genjet", pt_bins=PT_BINS+[(100, 2000)])


if __name__ == '__main__':
    # do_reco_plots()
    do_gen_plots()
