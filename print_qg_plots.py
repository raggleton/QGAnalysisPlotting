#!/usr/bin/env python

"""Print QG plots"""

import ROOT
# from TDRStyle import TDR_Style
from MyStyle import My_Style
My_Style.cd()
from comparator import Contribution, Plot, grab_obj
import common_utils as cu
import bisect
import numpy as np
from itertools import chain
import os
from array import array
from collections import namedtuple

# For debugging
import sys
import tracers
# sys.settrace(tracers.trace_calls)
# sys.settrace(tracers.trace_calls_detail)
# sys.settrace(tracers.trace_calls_and_returns)


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)


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

Angle = namedtuple("Angle", ['var', 'kappa', 'beta', 'name', "lambda_str"])
COMMON_VARS = [
    Angle('jet_multiplicity', 0, 0, "Multiplicity", "#lambda_{0}^{0}"),
    Angle('jet_pTD', 2, 0, "(p_{T}^{D})^{2}", "#lambda_{0}^{2}"),
    Angle('jet_LHA', 1, 0.5, "LHA", "#lambda_{0.5}^{1}"),
    Angle('jet_width', 1, 1, "Width", "#lambda_{1}^{1}"),
    Angle('jet_thrust', 1, 2, "Thrust", "#lambda_{2}^{1}"),
    Angle('jet_flavour', 0, 0, "", ""),
    Angle("jet_genParton_flavour", 0, 0, "", "")
]


ZPJ_RECOJET_RDIR = "ZPlusJets_QG"
DJ_RECOJET_RDIR = "Dijet_QG"
ZPJ_GENJET_RDIR = "ZPlusJets_genjet"
DJ_GENJET_RDIR = "Dijet_genjet"

PYTHIA_AK4_DIR = "workdir_ak4chs"
HERWIG_AK4_DIR = "workdir_ak4chs_herwig"
HERWIG_AK4_REWEIGHTED_DIR = "workdir_ak4chs_herwig_reweight"

AK4_GENJET_DIR = PYTHIA_AK4_DIR
AK4_GENJET_DIR = HERWIG_AK4_DIR
AK4_GENJET_DIR = HERWIG_AK4_REWEIGHTED_DIR
# AK4_GENJET_DIR = "workdir_ak4chs_thirdjetcut"

AK8_GENJET_DIR = "workdir_ak4chs_ak8genjet"
AK8_GENJET_DIR = "workdir_ak4chs_ak8genjet_herwig"

CHS_DIR = AK4_GENJET_DIR
TITLE_STR = "ak4 GenJet"
# TITLE_STR = "ak4 GenJet (w/extra jet cut)"

# CHS_DIR = AK8_GENJET_DIR
# TITLE_STR = "ak8 GenJet"

PUPPI_DIR = "workdir_ak4puppi"

# Controls all the things!!
ROOT_DIR = CHS_DIR
# ROOT_DIR = PUPPI_DIR
# Use this for data plots
# TITLE_STR = "[%s]" % ROOT_DIR.replace("workdir_", "")


PT_BINS = [(80, 100), (100, 200), (400, 500), (1000, 2000), (80, 2000)]
# PT_BINS = [(80, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 700), (700, 1000), (1000, 2000)]
THEORY_PT_BINS = [(100, 200), (400, 500), (1000, 2000), (80, 2000)]
THEORY_PT_BINS = [(100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 800), (800, 1000), (1000, 2000), (80, 2000)]


def do_comparison_plot(entries, output_filename, rebin=1, **plot_kwargs):
    """Plot several different objects on a single plot

    entries : list of 2-tuples, with (object, dict), where the dict is a set of kwargs passed to the Contribution object
    plot_kwargs : any other kwargs to be passed to the Plot object ctor
    """
    conts = [Contribution(ent[0], normalise_hist=True, fill_style=0, rebin_hist=rebin, **ent[1]) for ent in entries]
    p = Plot(conts, what="hist", subplot=conts[0], ytitle="p.d.f", **plot_kwargs)
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
    hproj = h2d.ProjectionX(ROOT.TUUID().AsString(), bin_start, bin_end, "eo")
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
    canvas = ROOT.TCanvas(ROOT.TUUID().AsString(), "", 800, 800)
    canvas.SetTicks(1, 1)
    canvas.SetLeftMargin(0.13)
    canvas.SetBottomMargin(0.11)
    obj_renorm.Draw("COLZ")
    obj_renorm.GetYaxis().SetTitleOffset(1.7)
    obj_renorm.GetXaxis().SetTitleOffset(1.2)
    odir = os.path.dirname(os.path.abspath(output_filename))
    if not os.path.isdir(odir):
        os.makedirs(odir)
    canvas.SaveAs(output_filename)


def get_ddelta_plot(one_hist, other_hist):
    """Make ddelta/dX plot

    Calculated as:
    0.5 * (x-y)^2 / (x+y)
    """
    new_hist = one_hist.Clone(ROOT.TUUID().AsString())
    new_hist.Add(other_hist, -1.)
    new_hist.Multiply(new_hist)
    sum_hist = one_hist.Clone(ROOT.TUUID().AsString())
    sum_hist.Add(other_hist)
    new_hist.Divide(sum_hist)
    new_hist.Scale(0.5)
    return new_hist


def calculate_delta(one_hist, other_hist):
    dist = get_ddelta_plot(one_hist, other_hist)
    return dist.Integral()*dist.GetBinWidth(1)


def calculate_delta(diff_hist):
    return diff_hist.Integral()*diff_hist.GetBinWidth(1)


def do_all_2D_plots(plot_dir="plots_2d", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG", var_list=None, var_prepend=""):
    """Do 2D distributions"""
    var_list = var_list or COMMON_VARS
    for ang in var_list:
        v = "%s%s_vs_pt" % (var_prepend, ang.var)

        rebin = [2, 4]
        if v == "jet_multiplicity_vs_pt":
            rebin = [2, 4]
        elif v == "jet_thrust_vs_pt" or "flavour" in v:
            rebin = [1, 4]

        recolour = False if "flavour" in v else True

        for rn in ['Y', None]:  # different renormalisation axes
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
                           var_list=None, var_prepend="", pt_bins=None, subplot_type="diff"):
    """Do pt/eta/nvtx binned 1D plots"""
    var_list = var_list or COMMON_VARS[2:]
    pt_bins = pt_bins or PT_BINS

    for ang in var_list:
        v = "%s%s_vs_pt" % (var_prepend, ang.var)

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
            rebin = 2
        elif "flavour" in v or "thrust" in v:
            rebin = 1

        xlim = None
        if "thrust" in v or "pTD" in v:
            xlim = (0, 0.5)

        ylim = None
        if "flavour" in v:
            ylim = (0, 1)
        elif "LHA" in v:
            ylim = (0, 5)

        for (start_val, end_val) in pt_bins:
            entries = [
                (get_projection_plot(h2d_dyj, start_val, end_val), dy_kwargs),
                (get_projection_plot(h2d_qcd, start_val, end_val), qcd_kwargs),
            ]
            do_comparison_plot(entries, "%s/%s/ptBinned/%s_pt%dto%d.pdf" % (ROOT_DIR, plot_dir, v, start_val, end_val),
                               rebin=rebin, title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val),
                               xtitle=ang.name + " (" + ang.lambda_str + ")",
                               xlim=xlim, ylim=ylim, subplot_type=subplot_type)

            if "flavour" in v:
                continue
            continue
            entries = [
                (get_projection_plot(h2d_dyj_q, start_val, end_val), dy_kwargs_q),
                (get_projection_plot(h2d_qcd_g, start_val, end_val), qcd_kwargs_g),
            ]
            do_comparison_plot(entries, "%s/%s/ptBinned/%s_pt%dto%d_flavMatched.pdf" % (ROOT_DIR, plot_dir, v, start_val, end_val),
                               rebin=rebin, title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val),
                               xtitle=ang.name + " (" + ang.lambda_str + ")",
                               xlim=xlim, subplot_type=subplot_type)


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
        p = Plot(contribs, what='graph', xtitle="p_{T}^{jet} [GeV]", ytitle=ytitle, title=title, ylim=(0, 1))
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
    """Plot all the sample/selection.flavour combinations to check distributions indep of sample"""
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


def do_jet_algo_comparison_plots(plot_dir="compare_jet_algo", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG",
                                 var_list=None, var_prepend="", pt_bins=None, subplot_type="diff"):
    """Do pt/eta/nvtx binned 1D plots, comparing different jet algos"""
    var_list = var_list or COMMON_VARS
    pt_bins = pt_bins or PT_BINS

    for ang in var_list:
        v = "%s%s_vs_pt" % (var_prepend, ang.var)

        h2d_dyj_ak4 = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % AK4_GENJET_DIR, "%s/%s" % (zpj_dirname, v))
        h2d_dyj_ak8 = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % AK8_GENJET_DIR, "%s/%s" % (zpj_dirname, v))
        h2d_qcd_ak4 = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % AK4_GENJET_DIR, "%s/%s" % (dj_dirname, v))
        h2d_qcd_ak8 = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % AK8_GENJET_DIR, "%s/%s" % (dj_dirname, v))

        if "flavour" not in v:
            h2d_dyj_q_ak4 = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % AK4_GENJET_DIR, "%s/q%s" % (zpj_dirname, v))
            h2d_dyj_q_ak8 = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % AK8_GENJET_DIR, "%s/q%s" % (zpj_dirname, v))
            h2d_qcd_g_ak4 = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % AK4_GENJET_DIR, "%s/g%s" % (dj_dirname, v))
            h2d_qcd_g_ak8 = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % AK8_GENJET_DIR, "%s/g%s" % (dj_dirname, v))

        lw = 2
        dy_kwargs_ak4 = dict(line_color=DY_COLOUR, fill_color=DY_COLOUR, label=DY_ZpJ_LABEL + " [ak4]", line_width=lw)
        dy_kwargs_ak8 = dict(line_color=DY_COLOUR, fill_color=DY_COLOUR, label=DY_ZpJ_LABEL + " [ak8]", line_width=lw, line_style=2)
        qcd_kwargs_ak4 = dict(line_color=QCD_COLOUR, fill_color=QCD_COLOUR, label=QCD_Dijet_LABEL + " [ak4]", line_width=lw)
        qcd_kwargs_ak8 = dict(line_color=QCD_COLOUR, fill_color=QCD_COLOUR, label=QCD_Dijet_LABEL + " [ak8]", line_width=lw, line_style=2)

        dy_kwargs_q_ak4 = dict(line_color=DY_COLOUR, fill_color=DY_COLOUR, label=DY_ZpJ_QFLAV_LABEL + " [ak4]", line_width=lw)
        dy_kwargs_q_ak8 = dict(line_color=DY_COLOUR, fill_color=DY_COLOUR, label=DY_ZpJ_QFLAV_LABEL + " [ak8]", line_width=lw, line_style=2)
        qcd_kwargs_g_ak4 = dict(line_color=QCD_COLOUR, fill_color=QCD_COLOUR, label=QCD_Dijet_GFLAV_LABEL + " [ak4]", line_width=lw)
        qcd_kwargs_g_ak8 = dict(line_color=QCD_COLOUR, fill_color=QCD_COLOUR, label=QCD_Dijet_GFLAV_LABEL + " [ak8]", line_width=lw, line_style=2)

        rebin = 2
        if v == "jet_multiplicity_vs_pt":
            rebin = 4
        elif "flavour" in v or "thrust" in v or "pTD" in v:
            rebin = 1

        xlim = None
        if "thrust" in v or "pTD" in v:
            xlim = (0, 0.5)

        ylim = None
        if "flavour" in v:
            ylim = (0, 1)

        for (start_val, end_val) in pt_bins:
            entries = [
                (get_projection_plot(h2d_dyj_ak4, start_val, end_val), dy_kwargs_ak4),
                (get_projection_plot(h2d_qcd_ak4, start_val, end_val), qcd_kwargs_ak4),
                (get_projection_plot(h2d_dyj_ak8, start_val, end_val), dy_kwargs_ak8),
                (get_projection_plot(h2d_qcd_ak8, start_val, end_val), qcd_kwargs_ak8),
            ]

            do_comparison_plot(entries, "%s/%s/ptBinned/%s_pt%dto%d.pdf" % (ROOT_DIR, plot_dir, v, start_val, end_val),
                               rebin=rebin, title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val),
                               xlim=xlim, ylim=ylim, subplot_type=subplot_type)

            if "flavour" in v:
                continue

            entries = [
                (get_projection_plot(h2d_dyj_q_ak4, start_val, end_val), dy_kwargs_q_ak4),
                (get_projection_plot(h2d_dyj_q_ak8, start_val, end_val), dy_kwargs_q_ak8),
                (get_projection_plot(h2d_qcd_g_ak4, start_val, end_val), qcd_kwargs_g_ak4),
                (get_projection_plot(h2d_qcd_g_ak8, start_val, end_val), qcd_kwargs_g_ak8),
            ]

            do_comparison_plot(entries, "%s/%s/ptBinned/%s_pt%dto%d_flavMatched.pdf" % (ROOT_DIR, plot_dir, v, start_val, end_val),
                               rebin=rebin, title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val), xlim=xlim, subplot_type=subplot_type)


def plot_ddelta(ddelta_hist, output_filename, xtitle, ytitle, title=""):
    cont = Contribution(ddelta_hist)
    p = Plot([cont], what="hist", legend=None, xtitle=xtitle, ytitle=ytitle, title=title)
    p.plot("HISTE")
    p.save(output_filename)


def construct_deltas_graph(deltas):
    N = len(deltas)
    gr = ROOT.TGraphErrors(N, array('d', range(N)), array('d', deltas), array('d', [0.5]*N), array('d', [0.00001]*N))
    return gr


def do_deltas_plot(graph_contribs, output_filename, bin_labels, title="", xtitle=""):
    do_legend = len(graph_contribs) > 1
    p = Plot(graph_contribs, what="graph", title=title, xtitle=xtitle, ytitle="Separation #Delta", legend=do_legend)
    if do_legend:
        p.legend.SetX1(0.7)
        p.legend.SetY1(0.15)
        p.legend.SetY2(0.35)
    p.plot("AP")
    # p.container.GetXaxis().LabelsOption("h")
    xax = p.container.GetXaxis()
    for i, lab in enumerate(bin_labels):
        bin_ind = xax.FindBin(i)  # need this as they don't correspond at all!
        p.container.GetXaxis().SetBinLabel(bin_ind, lab)
    p.container.SetMinimum(0)
    p.save(output_filename)


def do_pt_min_delta_plots(sources, plot_dir="deltas_ptmin", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG", var_list=None, var_prepend="", flavour_tag=False, save_component_hists=False):
    """Do plots comparing power of different ptMin cuts"""
    var_list = var_list or COMMON_VARS
    ptmin_bins = [50, 100, 200, 400, 800][:-1]
    dyj_flav = "q" if flavour_tag else ""
    dj_flav = "g" if flavour_tag else ""
    output_append = "_flavMatched" if flavour_tag else ""
    for ang in var_list:
        v = "%s%s_vs_pt" % (var_prepend, ang.var)
        graph_contribs, bin_labels = [], []

        for source_ind, source in enumerate(sources):
            deltas, conts = [], []

            for ind, pt_min in enumerate(ptmin_bins, 1):
                h2d_dyj = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % source['root_dir'], "%s_ptMin_%d/%s%s" % (zpj_dirname, pt_min, dyj_flav, v))
                h2d_qcd = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % source['root_dir'], "%s_ptMin_%d/%s%s" % (dj_dirname, pt_min, dj_flav, v))
                start_val, end_val = 80, 2000
                h_dy = get_projection_plot(h2d_dyj, start_val, end_val)
                if (h_dy.Integral()>0):
                    h_dy.Scale(1./(h_dy.GetBinWidth(1)*h_dy.Integral()))

                h_qcd = get_projection_plot(h2d_qcd, start_val, end_val)
                if (h_qcd.Integral()>0):
                    h_qcd.Scale(1./(h_qcd.GetBinWidth(1)*h_qcd.Integral()))

                ddelta_hist = get_ddelta_plot(h_dy, h_qcd)
                conts.append(Contribution(ddelta_hist, line_width=1, line_style=ind, line_color=(ind*10)+44, fill_color=(ind*10)+44, label="p_{T}^{Min} = %d GeV" % pt_min, rebin_hist=2))
                deltas.append(calculate_delta(ddelta_hist))

                if source_ind == 0:
                    bin_labels.append("%d" % pt_min)

                if save_component_hists:
                    plot_ddelta(ddelta_hist, "%s/%s/%s_ddelta_ptMin_%d%s.pdf" % (source['root_dir'], plot_dir, ang.var, pt_min, output_append),
                                xtitle=ang.name + " (" + ang.lambda_str + ")", ytitle="d#Delta/d" + ang.lambda_str)

            if save_component_hists:
                p = Plot(conts, what="hist", xtitle=ang.name, ytitle="p.d.f")
                p.plot("NOSTACK HISTE")
                p.save("%s/%s/%s_ddelta_ptMin_comparison%s.pdf" % (source['root_dir'], plot_dir, ang.var, output_append))

            gr = construct_deltas_graph(deltas)
            gr.SetName(source.get("label", ""))
            c = Contribution(gr, label=source.get("label", ""), marker_style=0, line_width=2, **source.get("style", {}))
            graph_contribs.append(c)

        do_deltas_plot(graph_contribs, "%s/%s/ptMins_%s%s.pdf" % (ROOT_DIR, plot_dir, ang.var, output_append),
                       bin_labels=bin_labels, title=TITLE_STR + ", %s [%s]" % (ang.name, ang.lambda_str), xtitle="p_{T}^{min} [GeV]")


def do_angularity_delta_plots(sources, plot_dir="delta_angularities", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG", var_list=None, var_prepend="", pt_bins=None, flavour_tag=False, save_component_hists=False):
    """Do plots comparing power of different angularities"""
    var_list = var_list or COMMON_VARS
    pt_bins = pt_bins or PT_BINS
    h2d_dyj = None
    h2d_qcd = None
    dyj_flav = "q" if flavour_tag else ""
    dj_flav = "g" if flavour_tag else ""
    output_append = "_flavMatched" if flavour_tag else ""
    for (start_val, end_val) in pt_bins:
        graph_contribs, bin_labels = [], []

        for source_ind, source in enumerate(sources):
            deltas = []

            # construct a graph of angularities for this source
            for ind, ang in enumerate(var_list, 1):
                v = "%s%s_vs_pt" % (var_prepend, ang.var)

                h2d_dyj = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % source['root_dir'], "%s/%s%s" % (zpj_dirname, dyj_flav, v))
                h2d_qcd = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % source['root_dir'], "%s/%s%s" % (dj_dirname, dj_flav, v))

                h_dy = get_projection_plot(h2d_dyj, start_val, end_val)
                if (h_dy.Integral() > 0):
                    h_dy.Scale(1./(h_dy.GetBinWidth(1)*h_dy.Integral()))

                h_qcd = get_projection_plot(h2d_qcd, start_val, end_val)
                if (h_qcd.Integral() > 0):
                    h_qcd.Scale(1./(h_qcd.GetBinWidth(1)*h_qcd.Integral()))

                ddelta_hist = get_ddelta_plot(h_dy, h_qcd)
                deltas.append(calculate_delta(ddelta_hist))

                if source_ind == 0:
                    bin_labels.append("#splitline{%s}{%s}" % (ang.name, ang.lambda_str))

                if save_component_hists:
                    plot_ddelta(ddelta_hist, "%s/%s/angularities_pt%dto%d_ddelta_%s%s.pdf" % (source['root_dir'], plot_dir, start_val, end_val, ang.var, output_append),
                                xtitle=ang.name + " (" + ang.lambda_str + ")", ytitle="d#Delta/d" + ang.lambda_str)

            gr = construct_deltas_graph(deltas)
            gr.SetName(source.get("label", ""))
            c = Contribution(gr, label=source.get("label", ""), marker_style=0, line_width=2, **source.get("style", {}))
            graph_contribs.append(c)

        do_deltas_plot(graph_contribs, "%s/%s/angularities_pt%dto%d%s.pdf" % (ROOT_DIR, plot_dir, start_val, end_val, output_append),
                       bin_labels=bin_labels, title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val), xtitle="Angularity: (#kappa, #beta)")


def do_gen_reco_comparison_plots(var_list=None, gen_var_prepend="gen", reco_var_prepend="",
                                 plot_dir="plot_reco_gen", zpj_reco_dirname=ZPJ_RECOJET_RDIR, dj_reco_dirname=DJ_RECOJET_RDIR,
                                 zpj_gen_dirname=ZPJ_GENJET_RDIR, dj_gen_dirname=DJ_GENJET_RDIR, pt_bins=THEORY_PT_BINS, subplot_type=None):
    var_list = var_list or COMMON_VARS[:-1]
    for ang in var_list:
        v_reco = "%s%s_vs_pt" % (reco_var_prepend, ang.var)
        v_gen = "%s%s_vs_pt" % (gen_var_prepend, ang.var)

        h2d_dyj_reco = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "%s/%s" % (zpj_reco_dirname, v_reco))
        h2d_qcd_reco = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "%s/%s" % (dj_reco_dirname, v_reco))
        h2d_dyj_gen = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "%s/%s" % (zpj_gen_dirname, v_gen))
        h2d_qcd_gen = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "%s/%s" % (dj_gen_dirname, v_gen))

        if "flavour" not in v_reco:
            h2d_dyj_reco_q = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "%s/q%s" % (zpj_reco_dirname, v_reco))
            h2d_qcd_reco_g = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "%s/g%s" % (dj_reco_dirname, v_reco))
            h2d_dyj_gen_q = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "%s/q%s" % (zpj_gen_dirname, v_gen))
            h2d_qcd_gen_g = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "%s/g%s" % (dj_gen_dirname, v_gen))


        for (start_val, end_val) in pt_bins:
            lw = 2

            dy_reco_kwargs = dict(line_color=DY_COLOUR, fill_color=DY_COLOUR, label=DY_ZpJ_LABEL + " [RecoJet]", line_width=lw)
            qcd_reco_kwargs = dict(line_color=QCD_COLOUR, fill_color=QCD_COLOUR, label=QCD_Dijet_LABEL + " [RecoJet]", line_width=lw)

            dy_gen_kwargs = dict(line_color=DY_COLOUR, fill_color=DY_COLOUR, label=DY_ZpJ_LABEL + " [GenJet]", line_width=lw, line_style=2)
            qcd_gen_kwargs = dict(line_color=QCD_COLOUR, fill_color=QCD_COLOUR, label=QCD_Dijet_LABEL + " [GenJet]", line_width=lw, line_style=2)


            entries = [
                (get_projection_plot(h2d_dyj_reco, start_val, end_val), dy_reco_kwargs),
                (get_projection_plot(h2d_qcd_reco, start_val, end_val), qcd_reco_kwargs),
                (get_projection_plot(h2d_dyj_gen, start_val, end_val), dy_gen_kwargs),
                (get_projection_plot(h2d_qcd_gen, start_val, end_val), qcd_gen_kwargs)
            ]

            rebin = 2
            xlim = None
            ylim = None
            if "flavour" in v_reco:
                rebin = 1
                ylim = (0, 1)
            if "thrust" in v_reco:
                xlim = (0, 0.5)

            do_comparison_plot(entries, "%s/%s/%s_pt%dto%d.pdf" % (ROOT_DIR, plot_dir, ang.var, start_val, end_val),
                               rebin=rebin, title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val),
                               xtitle=ang.name + " (" + ang.lambda_str + ")",
                               xlim=xlim, ylim=ylim, subplot_type=subplot_type)

            # Do flavour-tagged comparison
            if "flavour" in v_reco:
                continue

            dy_reco_kwargs['label'] = DY_ZpJ_QFLAV_LABEL + " [RecoJet]"
            qcd_reco_kwargs['label'] = QCD_Dijet_GFLAV_LABEL + " [RecoJet]"
            dy_gen_kwargs['label'] = DY_ZpJ_QFLAV_LABEL + " [GenJet]"
            qcd_gen_kwargs['label'] = QCD_Dijet_GFLAV_LABEL + " [GenJet]"

            entries = [
                (get_projection_plot(h2d_dyj_reco_q, start_val, end_val), dy_reco_kwargs),
                (get_projection_plot(h2d_qcd_reco_g, start_val, end_val), qcd_reco_kwargs),
                (get_projection_plot(h2d_dyj_gen_q, start_val, end_val), dy_gen_kwargs),
                (get_projection_plot(h2d_qcd_gen_g, start_val, end_val), qcd_gen_kwargs)
            ]

            do_comparison_plot(entries, "%s/%s/%s_flavMatched_pt%dto%d.pdf" % (ROOT_DIR, plot_dir, ang.var, start_val, end_val),
                               rebin=rebin, title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val),
                               xtitle=ang.name + " (" + ang.lambda_str + ")",
                               xlim=xlim, ylim=ylim, subplot_type=subplot_type)


def do_all_exclusive_plots_comparison(sources, plot_dir="plots_dy_vs_qcd", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG",
                                      var_list=None, var_prepend="", pt_bins=None, subplot_type="diff", do_flav_tagged=True):
    """Do 1D plots, comparing various sources. FOr each source plots DY & QCD samples."""
    var_list = var_list or COMMON_VARS[2:]
    pt_bins = pt_bins or PT_BINS

    for ang in var_list:

        v = "%s%s_vs_pt" % (var_prepend, ang.var)
        for (start_val, end_val) in pt_bins:
            entries_normal, entries_flav = [], []

            # Get all plots
            for source in sources:

                h2d_dyj = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % source['root_dir'], "%s/%s" % (zpj_dirname, v))
                h2d_qcd = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % source['root_dir'], "%s/%s" % (dj_dirname, v))
                lw = 2
                dy_kwargs = dict(line_color=DY_COLOUR, fill_color=DY_COLOUR, label=DY_ZpJ_LABEL + " " + source.get('label', ''), line_width=lw)
                dy_kwargs.update(source.get('style', {}))
                dy_kwargs.update(source.get('dy_style', {}))

                qcd_kwargs = dict(line_color=QCD_COLOUR, fill_color=QCD_COLOUR, label=QCD_Dijet_LABEL + " " + source.get('label', ''), line_width=lw)
                qcd_kwargs.update(source.get('style', {}))
                qcd_kwargs.update(source.get('qcd_style', {}))

                entries_normal.append((get_projection_plot(h2d_dyj, start_val, end_val), dy_kwargs))
                entries_normal.append((get_projection_plot(h2d_qcd, start_val, end_val), qcd_kwargs))

                if not do_flav_tagged or "flavour" in v:
                    continue

                # Flav tagged plots
                h2d_dyj_q = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % source['root_dir'], "%s/q%s" % (zpj_dirname, v))
                h2d_qcd_g = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % source['root_dir'], "%s/g%s" % (dj_dirname, v))

                dy_kwargs_q = dict(line_color=DY_COLOUR, fill_color=DY_COLOUR, label=DY_ZpJ_QFLAV_LABEL + " " + source.get('label', ''), line_width=lw)
                dy_kwargs_q.update(source.get('style', {}))
                dy_kwargs_q.update(source.get('dy_style', {}))

                qcd_kwargs_g = dict(line_color=QCD_COLOUR, fill_color=QCD_COLOUR, label=QCD_Dijet_GFLAV_LABEL + " " + source.get('label', ''), line_width=lw)
                qcd_kwargs_g.update(source.get('style', {}))
                qcd_kwargs_g.update(source.get('qcd_style', {}))

                entries_flav.append((get_projection_plot(h2d_dyj_q, start_val, end_val), dy_kwargs_q))
                entries_flav.append((get_projection_plot(h2d_qcd_g, start_val, end_val), qcd_kwargs_g))

            rebin = 2
            if v == "jet_multiplicity_vs_pt":
                rebin = 2
            elif "flavour" in v or "thrust" in v or 'pTD' in v:
                rebin = 1

            xlim = None
            if "thrust" in v or "pTD" in v:
                xlim = (0, 0.5)

            ylim = None
            if "flavour" in v:
                ylim = (0, 1)
            elif "LHA" in v:
                ylim = (0, 5)

            do_comparison_plot(entries_normal, "%s/%s/ptBinned/%s_pt%dto%d.pdf" % (ROOT_DIR, plot_dir, v, start_val, end_val),
                               rebin=rebin, title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val),
                               xtitle=ang.name + " (" + ang.lambda_str + ")",
                               xlim=xlim, ylim=ylim, subplot_type=subplot_type)
            if do_flav_tagged:
                do_comparison_plot(entries_flav, "%s/%s/ptBinned/%s_pt%dto%d_flavMatched.pdf" % (ROOT_DIR, plot_dir, v, start_val, end_val),
                                   rebin=rebin, title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val),
                                   xtitle=ang.name + " (" + ang.lambda_str + ")",
                                   xlim=xlim, subplot_type=subplot_type)


def do_reco_plots():
    global TITLE_STR
    TITLE_STR = "[%s]" % ROOT_DIR.replace("workdir_", "")
    # do_all_2D_plots()
    do_all_exclusive_plots(subplot_type=None)
    # do_all_flavour_fraction_plots()
    # do_chs_vs_puppi_plots()
    # do_wrong_plots()
    sources = [
        {"root_dir": ROOT_DIR, 'label': "-", "style": {'line_style': 1}}
        # {"root_dir": ROOT_DIR, 'label': "Herwig", "style": {'line_style': 1}}
    ]
    # do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2])
    do_angularity_delta_plots(sources, var_list=COMMON_VARS[0:-2], pt_bins=THEORY_PT_BINS, save_component_hists=True)


def do_reco_comparison_plots():
    """Compare reco jets from different generators"""
    sources = [
        {"root_dir": PYTHIA_AK4_DIR, 'label': "Pythia", "style": {'line_style': 1}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Herwig", "style": {'line_style': 2}}

    ]
    do_all_exclusive_plots_comparison(sources=sources, var_list=COMMON_VARS[:-2],
                                      plot_dir="plots_dy_vs_qcd_compare_generators",
                                      subplot_type=None, do_flav_tagged=False, pt_bins=THEORY_PT_BINS)
    sources = [
        {"root_dir": PYTHIA_AK4_DIR, 'label': "Pythia", "style": {'line_style': 1, 'line_color': ROOT.kBlack}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Herwig", "style": {'line_style': 2, 'line_color': ROOT.kRed}}

    ]
    # do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2])
    do_angularity_delta_plots(sources, plot_dir="delta_angularities_compare_generators",
                              var_list=COMMON_VARS[0:-2], pt_bins=THEORY_PT_BINS)


def do_reco_reweight_comparison_plots():
    """Compare effect of reweighting to Pythia spectrum"""
    sources = [
        {"root_dir": HERWIG_AK4_REWEIGHTED_DIR, 'label': "Reweighted to Pythia", "style": {'line_style': 1}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Not reweighted", "style": {'line_style': 2}}
    ]
    do_all_exclusive_plots_comparison(sources=sources, var_list=COMMON_VARS[:-2],
                                      plot_dir="plots_dy_vs_qcd_compare_reweight",
                                      pt_bins=THEORY_PT_BINS, subplot_type=None, do_flav_tagged=False)

    sources = [
        {"root_dir": HERWIG_AK4_REWEIGHTED_DIR, 'label': "Reweighted to Pythia", "style": {'line_color': ROOT.kGreen}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Not reweighted", "style": {'line_color': ROOT.kAzure}}
    ]
    # do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2], plot_dir="deltas_ptMin_compare_reweight")
    do_angularity_delta_plots(sources, var_list=COMMON_VARS[:-2], plot_dir="deltas_angularities_compare_reweight", pt_bins=THEORY_PT_BINS)

    sources = [
        {"root_dir": PYTHIA_AK4_DIR, 'label': "Pythia", "style": {'line_color': ROOT.kBlack}},
        {"root_dir": HERWIG_AK4_REWEIGHTED_DIR, 'label': "Reweighted to Pythia", "style": {'line_color': ROOT.kGreen}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Not reweighted", "style": {'line_color': ROOT.kAzure}}
    ]
    # do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2], plot_dir="deltas_ptMin_compare_normal_and_reweight")
    do_angularity_delta_plots(sources, var_list=COMMON_VARS[:-2], plot_dir="deltas_angularities_compare_normal_and_reweight", pt_bins=THEORY_PT_BINS)


def do_gen_plots():
    global TITLE_STR
    TITLE_STR = "ak4 GenJet"
    # sources = [
    #     {"root_dir": ROOT_DIR, 'label': "", "style": {'line_style': 1}},
    # ]
    # do_all_exclusive_plots_comparison(sources=sources, var_list=COMMON_VARS[:-1], var_prepend="gen",
    #                                   plot_dir="plots_dy_vs_qcd_gen", zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR,
    #                                   pt_bins=THEORY_PT_BINS, subplot_type=None, do_flav_tagged=True)

    do_all_2D_plots(var_list=COMMON_VARS[:-1], var_prepend="gen", plot_dir="plots_2d_gen",
                    zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR)
    do_all_exclusive_plots(var_list=COMMON_VARS[:-1], var_prepend="gen", plot_dir="plots_dy_vs_qcd_gen",
                           zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR, pt_bins=THEORY_PT_BINS, subplot_type=None)
    # do_all_flavour_fraction_plots(var_prepend="gen", plot_dir="flav_fractions_gen",
    #                               zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR)
    # do_wrong_plots(var_prepend="gen", plot_dir="wrong_flavs_gen",
    #                zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR, pt_bins=THEORY_PT_BINS)
    # do_jet_algo_comparison_plots(var_list=COMMON_VARS[:-1], var_prepend="gen", plot_dir="compare_jet_algo",
    #                              zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR, pt_bins=THEORY_PT_BINS, subplot_type=None)
    # do_gen_reco_comparison_plots(var_list=COMMON_VARS[:-1], gen_var_prepend="gen", reco_var_prepend="",
    #                              plot_dir="plot_reco_gen", zpj_reco_dirname=ZPJ_RECOJET_RDIR, dj_reco_dirname=DJ_RECOJET_RDIR,
    #                              zpj_gen_dirname=ZPJ_GENJET_RDIR, dj_gen_dirname=DJ_GENJET_RDIR, pt_bins=THEORY_PT_BINS)

    # Separation plots
    sources = [
        {"root_dir": AK4_GENJET_DIR, 'label': "", "style": {'line_style': 1}},
    ]
    do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2], var_prepend="gen", plot_dir="deltas_ptMin_gen",
                          zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR, save_component_hists=True)
    do_angularity_delta_plots(sources, var_list=COMMON_VARS[:-2], var_prepend="gen", plot_dir="deltas_angularities_gen",
                          zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR, pt_bins=THEORY_PT_BINS, save_component_hists=True)
    # flav-tagged versions
    do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2], var_prepend="gen", plot_dir="deltas_ptMin_gen",
                          zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR, flavour_tag=True, save_component_hists=True)
    do_angularity_delta_plots(sources, var_list=COMMON_VARS[:-2], var_prepend="gen", plot_dir="deltas_angularities_gen",
                          zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR, pt_bins=THEORY_PT_BINS, flavour_tag=True, save_component_hists=True)



def do_gen_comparison_plots():
    """Compare genjets from different generators"""
    sources = [
        {"root_dir": PYTHIA_AK4_DIR, 'label': "Pythia", "style": {'line_style': 1}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Herwig", "style": {'line_style': 2}}
    ]
    do_all_exclusive_plots_comparison(sources=sources, var_list=COMMON_VARS[:-2], var_prepend="gen",
                                      plot_dir="plots_dy_vs_qcd_gen_compare_generators",
                                      zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR,
                                      pt_bins=THEORY_PT_BINS, subplot_type=None, do_flav_tagged=False)

    sources = [
        {"root_dir": PYTHIA_AK4_DIR, 'label': "Pythia", "style": {'line_color': ROOT.kBlack}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Herwig", "style": {'line_color': ROOT.kRed, 'line_style': 2}}
    ]
    do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2], var_prepend="gen", plot_dir="deltas_ptMin_gen_compare_generators",
                          zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR)
    do_angularity_delta_plots(sources, var_list=COMMON_VARS[:-2], var_prepend="gen", plot_dir="deltas_angularities_gen_compare_generators",
                          zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR, pt_bins=THEORY_PT_BINS)


def do_gen_reweight_comparison_plots():
    """Compare effect of reweighting to Pythia spectrum"""
    sources = [
        {"root_dir": HERWIG_AK4_REWEIGHTED_DIR, 'label': "Reweighted to Pythia", "style": {'line_style': 1}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Not reweighted", "style": {'line_style': 2}}
    ]
    do_all_exclusive_plots_comparison(sources=sources, var_list=COMMON_VARS[:-2], var_prepend="gen",
                                      plot_dir="plots_dy_vs_qcd_gen_compare_reweight",
                                      zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR,
                                      pt_bins=THEORY_PT_BINS, subplot_type=None, do_flav_tagged=False)

    sources = [
        {"root_dir": HERWIG_AK4_REWEIGHTED_DIR, 'label': "Reweighted to Pythia", "style": {'line_color': ROOT.kGreen}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Not reweighted", "style": {'line_color': ROOT.kAzure}}
    ]
    do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2], var_prepend="gen", plot_dir="deltas_ptMin_gen_compare_reweight",
                          zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR)
    do_angularity_delta_plots(sources, var_list=COMMON_VARS[:-2], var_prepend="gen", plot_dir="deltas_angularities_gen_compare_reweight",
                          zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR, pt_bins=THEORY_PT_BINS)

    sources = [
        {"root_dir": PYTHIA_AK4_DIR, 'label': "Pythia", "style": {'line_color': ROOT.kBlack}},
        {"root_dir": HERWIG_AK4_REWEIGHTED_DIR, 'label': "Reweighted to Pythia", "style": {'line_color': ROOT.kGreen}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Not reweighted", "style": {'line_color': ROOT.kAzure}}
    ]
    do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2], var_prepend="gen", plot_dir="deltas_ptMin_gen_compare_normal_and_reweight",
                          zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR)
    do_angularity_delta_plots(sources, var_list=COMMON_VARS[:-2], var_prepend="gen", plot_dir="deltas_angularities_gen_compare_normal_and_reweight",
                          zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR, pt_bins=THEORY_PT_BINS)


if __name__ == '__main__':
    do_reco_plots()
    do_reco_comparison_plots()
    do_reco_reweight_comparison_plots()
    do_gen_plots()
    do_gen_comparison_plots()
    do_gen_reweight_comparison_plots()
