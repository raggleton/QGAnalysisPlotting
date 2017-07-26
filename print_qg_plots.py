#!/usr/bin/env python

"""Print QG plots"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import bisect
import numpy as np
import os
from array import array

# My stuff
from comparator import Contribution, Plot, grab_obj
import common_utils as cu
from qg_common import *
from qg_general_plots import *
from qg_flavour_plots import *
from qg_delta_plots import *

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


PYTHIA_AK4_DIR = "workdir_ak4chs"
PYTHIA_ONLY_AK4_DIR = "workdir_ak4chs_pythiaOnlyFlat"
HERWIG_AK4_DIR = "workdir_ak4chs_herwig"
HERWIG_AK4_REWEIGHTED_DIR = "workdir_ak4chs_herwig_reweight"

AK4_GENJET_DIR = PYTHIA_AK4_DIR
AK4_GENJET_DIR = PYTHIA_ONLY_AK4_DIR
# AK4_GENJET_DIR = HERWIG_AK4_DIR
# AK4_GENJET_DIR = HERWIG_AK4_REWEIGHTED_DIR

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


def do_chs_vs_puppi_plots():
    sources = [
        {"root_dir": CHS_DIR, 'label': "CHS", "style": {'line_style': 1}},
        {"root_dir": PUPPI_DIR, 'label': "PUPPI", "style": {'line_style': 3}}
    ]
    do_all_exclusive_plots_comparison(sources, var_list=COMMON_VARS[2:],
                                      plot_dir=os.path.join(ROOT_DIR, "ak4_chs_vs_puppi"),
                                      zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG",
                                      var_prepend="", pt_bins=None,
                                      subplot_type=None, do_flav_tagged=True)


def do_wrong_plots(var_prepend="", plot_dir="wrong_flavs", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG", pt_bins=None):
    """Plot all the sample/selection.flavour combinations to check distributions indep of sample"""
    pt_bins = pt_bins or PT_BINS
    for v in ['jet_LHA', 'jet_pTD', 'jet_width', 'jet_thrust', 'jet_multiplicity']:
        v = "%s%s_vs_pt" % (var_prepend, v)

        h2d_dyj_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "%s/q%s" % (zpj_dirname, v))
        h2d_dyj_wrong_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "%s/g%s" % (zpj_dirname, v))
        h2d_dyj_qcd_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "%s/q%s" % (dj_dirname, v))
        h2d_dyj_qcd_wrong_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "%s/g%s" % (dj_dirname, v))
        h2d_qcd_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "%s/g%s" % (dj_dirname, v))
        h2d_qcd_wrong_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "%s/q%s" % (dj_dirname, v))

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
    """Do 1D plots comparing different jet algos"""
    sources = [
        {"root_dir": AK4_GENJET_DIR, 'label': "AK4", "style": {'line_style': 1}},
        {"root_dir": AK8_GENJET_DIR, 'label': "AK8", "style": {'line_style': 2}}
    ]
    do_all_exclusive_plots_comparison(sources, var_list=var_list,
                                      plot_dir=plot_dir,
                                      zpj_dirname=zpj_dirname, dj_dirname=dj_dirname,
                                      var_prepend=var_prepend, pt_bins=pt_bins,
                                      subplot_type=subplot_type, do_flav_tagged=True)


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


def do_reco_plots():
    global TITLE_STR
    TITLE_STR = "[%s]" % ROOT_DIR.replace("workdir_", "")
    # do_all_2D_plots()
    sources = [{"root_dir": ROOT_DIR, 'label': "", "style": {'line_style': 1}}]
    do_all_exclusive_plots_comparison(sources=sources, var_list=COMMON_VARS[:-1],
                                      plot_dir=os.path.join(ROOT_DIR, "plots_dy_vs_qcd"),
                                      pt_bins=THEORY_PT_BINS, subplot_type=None, do_flav_tagged=True)
    # do_all_flavour_fraction_plots()
    # do_chs_vs_puppi_plots()
    # do_wrong_plots()
    sources = [
        {"root_dir": ROOT_DIR, 'label': "-", "style": {'line_style': 1}}
        # {"root_dir": ROOT_DIR, 'label': "Herwig", "style": {'line_style': 1}}
    ]
    # do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2])
    do_angularity_delta_plots(sources, var_list=COMMON_VARS[0:-2], pt_bins=THEORY_PT_BINS, save_component_hists=True)


def do_reco_flav_split_plots():
    """Do plots comparing individual q flavours"""
    global TITLE_STR
    TITLE_STR = "[%s]" % ROOT_DIR.replace("workdir_", "")
    sources = [
        # {"root_dir": ROOT_DIR, 'label': "All", "style": {'line_style': 1}},
        {"root_dir": ROOT_DIR, 'label': "jet1 = u", "zpj_dirname": "ZPlusJets_QG_u", "dj_dirname": "Dijet_QG_u", "style": {'line_style': 1}},
        {"root_dir": ROOT_DIR, 'label': "jet1 = d", "zpj_dirname": "ZPlusJets_QG_d", "dj_dirname": "Dijet_QG_d", "style": {'line_style': 2}},
        {"root_dir": ROOT_DIR, 'label': "jet1 = s", "zpj_dirname": "ZPlusJets_QG_s", "dj_dirname": "Dijet_QG_s", "style": {'line_style': 3}}
    ]
    do_all_exclusive_plots_comparison(sources=sources, var_list=COMMON_VARS[:-1],
                                      plot_dir=os.path.join(ROOT_DIR, "plots_dy_vs_qcd_flav_split"),
                                      pt_bins=THEORY_PT_BINS, subplot_type=None, do_flav_tagged=False)


def do_reco_generator_comparison_plots():
    """Compare reco jets from different generators"""
    sources = [
        {"root_dir": PYTHIA_AK4_DIR, 'label': "Pythia", "style": {'line_style': 1}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Herwig", "style": {'line_style': 2}}

    ]
    do_all_exclusive_plots_comparison(sources=sources, var_list=COMMON_VARS[:-2],
                                      plot_dir=os.path.join(ROOT_DIR, "plots_dy_vs_qcd_compare_generators"),
                                      subplot_type=None, do_flav_tagged=False, pt_bins=THEORY_PT_BINS)
    sources = [
        {"root_dir": PYTHIA_AK4_DIR, 'label': "Pythia", "style": {'line_style': 1, 'line_color': ROOT.kBlack}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Herwig", "style": {'line_style': 2, 'line_color': ROOT.kRed}}

    ]
    # do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2])
    do_angularity_delta_plots(sources, plot_dir=os.path.join(ROOT_DIR, "delta_angularities_compare_generators"),
                              var_list=COMMON_VARS[0:-2], pt_bins=THEORY_PT_BINS)


def do_reco_reweight_comparison_plots():
    """Compare effect of reweighting to Pythia spectrum"""
    sources = [
        {"root_dir": HERWIG_AK4_REWEIGHTED_DIR, 'label': "Reweighted to Pythia", "style": {'line_style': 1}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Not reweighted", "style": {'line_style': 2}}
    ]
    do_all_exclusive_plots_comparison(sources=sources, var_list=COMMON_VARS[:-2],
                                      plot_dir=os.path.join(ROOT_DIR, "plots_dy_vs_qcd_compare_reweight"),
                                      pt_bins=THEORY_PT_BINS, subplot_type=None, do_flav_tagged=False)

    sources = [
        {"root_dir": HERWIG_AK4_REWEIGHTED_DIR, 'label': "Reweighted to Pythia", "style": {'line_color': ROOT.kGreen}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Not reweighted", "style": {'line_color': ROOT.kAzure}}
    ]
    # do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2], plot_dir=os.path.join(ROOT_DIR, "deltas_ptMin_compare_reweight"))
    do_angularity_delta_plots(sources, var_list=COMMON_VARS[:-2],
                              plot_dir=os.path.join(ROOT_DIR, "deltas_angularities_compare_reweight"),
                              pt_bins=THEORY_PT_BINS)

    sources = [
        {"root_dir": PYTHIA_AK4_DIR, 'label': "Pythia", "style": {'line_color': ROOT.kBlack}},
        {"root_dir": HERWIG_AK4_REWEIGHTED_DIR, 'label': "Reweighted to Pythia", "style": {'line_color': ROOT.kGreen}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Not reweighted", "style": {'line_color': ROOT.kAzure}}
    ]
    # do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2], plot_dir=os.path.join(ROOT_DIR, "deltas_ptMin_compare_normal_and_reweight"))
    do_angularity_delta_plots(sources, var_list=COMMON_VARS[:-2],
                              plot_dir=os.path.join(ROOT_DIR, "deltas_angularities_compare_normal_and_reweight"),
                              pt_bins=THEORY_PT_BINS)


def do_reco_pu_comparison_plots():
    """Compare by PU bins"""
    pu_bins = [(5, 15), (20, 25), (30, 40)]
    sources = []
    for ind, (pu_min, pu_max) in enumerate(pu_bins):
        sources.append({
            "root_dir": ROOT_DIR,
            'label': "PU %d-%d" % (pu_min, pu_max),
            'zpj_dirname': ZPJ_RECOJET_RDIR + "_PU_%d_to_%d" % (pu_min, pu_max),
            'dj_dirname': DJ_RECOJET_RDIR + "_PU_%d_to_%d" % (pu_min, pu_max),
            "style": {'line_style': 1, "line_width": 1},
            "dy_style": {'line_color': DY_COLOURS[ind], 'fill_color': DY_COLOURS[ind]},
            "qcd_style": {'line_color': QCD_COLOURS[ind], 'fill_color': QCD_COLOURS[ind]}
        })
    # do_all_exclusive_plots_comparison(sources=sources, var_list=COMMON_VARS[:-2], zpj_dirname=None,
    #                                   plot_dir=os.path.join(ROOT_DIR, "plots_dy_vs_qcd_compare_pu_dijet"),
    #                                   pt_bins=THEORY_PT_BINS, subplot_type="ratio", do_flav_tagged=False)
    # do_all_exclusive_plots_comparison(sources=sources, var_list=COMMON_VARS[:-2], dj_dirname=None,
    #                                   plot_dir=os.path.join(ROOT_DIR, "plots_dy_vs_qcd_compare_pu_zpj"),
    #                                   pt_bins=THEORY_PT_BINS, subplot_type="ratio", do_flav_tagged=False)


    for ind, s in enumerate(sources):
        sources[ind]['style']['line_width'] = 2
        sources[ind]['style']['line_color'] = DY_COLOURS[ind]
        if ind == 2:
            do_angularity_delta_plots(sources[ind:ind+1], var_list=COMMON_VARS[2:-2],
                              plot_dir=os.path.join(ROOT_DIR, "deltas_angularities_compare_pu_PU_%d_to_%d" % (pu_bins[ind][0], pu_bins[ind][1])),
                              pt_bins=THEORY_PT_BINS, save_component_hists=True)

    do_angularity_delta_plots(sources, var_list=COMMON_VARS[:-2],
                              plot_dir=os.path.join(ROOT_DIR, "deltas_angularities_compare_pu"),
                              pt_bins=THEORY_PT_BINS)


def do_gen_plots():
    global TITLE_STR
    TITLE_STR = "ak4 GenJet"
    sources = [{"root_dir": ROOT_DIR, 'label': "", "style": {'line_style': 1}}]
    do_all_exclusive_plots_comparison(sources=sources, var_list=COMMON_VARS[:-1], var_prepend="gen",
                                      plot_dir=os.path.join(ROOT_DIR, "plots_dy_vs_qcd_gen"),
                                      zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR,
                                      pt_bins=THEORY_PT_BINS, subplot_type=None, do_flav_tagged=True)
    do_all_2D_plots(var_list=COMMON_VARS[:-1], var_prepend="gen", plot_dir="plots_2d_gen",
                    zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR)
    # do_all_flavour_fraction_plots(var_prepend="gen", plot_dir="flav_fractions_gen",
    #                               zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR)
    # do_wrong_plots(var_prepend="gen", plot_dir="wrong_flavs_gen",
    #                zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR, pt_bins=THEORY_PT_BINS)
    # do_jet_algo_comparison_plots(var_list=COMMON_VARS[:-1], var_prepend="gen", plot_dir=os.path.join(ROOT_DIR, "compare_jet_algo"),
    #                              zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR, pt_bins=THEORY_PT_BINS, subplot_type=None)
    # do_gen_reco_comparison_plots(var_list=COMMON_VARS[:-1], gen_var_prepend="gen", reco_var_prepend="",
    #                              plot_dir="plot_reco_gen", zpj_reco_dirname=ZPJ_RECOJET_RDIR, dj_reco_dirname=DJ_RECOJET_RDIR,
    #                              zpj_gen_dirname=ZPJ_GENJET_RDIR, dj_gen_dirname=DJ_GENJET_RDIR, pt_bins=THEORY_PT_BINS)

    # Separation plots
    sources = [
        {"root_dir": AK4_GENJET_DIR, 'label': "", "style": {'line_style': 1}},
    ]
    do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2], var_prepend="gen",
                          plot_dir=os.path.join(ROOT_DIR, "deltas_ptMin_gen"),
                          zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR, save_component_hists=True)
    do_angularity_delta_plots(sources, var_list=COMMON_VARS[:-2], var_prepend="gen",
                              plot_dir=os.path.join(ROOT_DIR, "deltas_angularities_gen"),
                              zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR,
                              pt_bins=THEORY_PT_BINS, save_component_hists=True)
    # flav-tagged versions
    do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2], var_prepend="gen",
                          plot_dir=os.path.join(ROOT_DIR, "deltas_ptMin_gen"),
                          zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR,
                          flavour_tag=True, save_component_hists=True)
    do_angularity_delta_plots(sources, var_list=COMMON_VARS[:-2], var_prepend="gen",
                              plot_dir=os.path.join(ROOT_DIR, "deltas_angularities_gen"),
                              zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR,
                              pt_bins=THEORY_PT_BINS, flavour_tag=True, save_component_hists=True)



def do_gen_generator_comparison_plots():
    """Compare genjets from different generators"""
    sources = [
        {"root_dir": PYTHIA_AK4_DIR, 'label': "Pythia", "style": {'line_style': 1}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Herwig", "style": {'line_style': 2}}
    ]
    do_all_exclusive_plots_comparison(sources=sources, var_list=COMMON_VARS[:-2], var_prepend="gen",
                                      plot_dir=os.path.join(ROOT_DIR, "plots_dy_vs_qcd_gen_compare_generators"),
                                      zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR,
                                      pt_bins=THEORY_PT_BINS, subplot_type=None, do_flav_tagged=False)

    sources = [
        {"root_dir": PYTHIA_AK4_DIR, 'label': "Pythia", "style": {'line_color': ROOT.kBlack}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Herwig", "style": {'line_color': ROOT.kRed, 'line_style': 2}}
    ]
    do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2], var_prepend="gen",
                          plot_dir=os.path.join(ROOT_DIR, "deltas_ptMin_gen_compare_generators"),
                          zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR)
    do_angularity_delta_plots(sources, var_list=COMMON_VARS[:-2], var_prepend="gen",
                              plot_dir=os.path.join(ROOT_DIR, "deltas_angularities_gen_compare_generators"),
                              zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR, pt_bins=THEORY_PT_BINS)


def do_gen_reweight_comparison_plots():
    """Compare effect of reweighting to Pythia spectrum"""
    sources = [
        {"root_dir": HERWIG_AK4_REWEIGHTED_DIR, 'label': "Reweighted to Pythia", "style": {'line_style': 1}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Not reweighted", "style": {'line_style': 2}}
    ]
    do_all_exclusive_plots_comparison(sources=sources, var_list=COMMON_VARS[:-2], var_prepend="gen",
                                      plot_dir=os.path.join(ROOT_DIR, "plots_dy_vs_qcd_gen_compare_reweight"),
                                      zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR,
                                      pt_bins=THEORY_PT_BINS, subplot_type=None, do_flav_tagged=False)

    sources = [
        {"root_dir": HERWIG_AK4_REWEIGHTED_DIR, 'label': "Reweighted to Pythia", "style": {'line_color': ROOT.kGreen}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Not reweighted", "style": {'line_color': ROOT.kAzure}}
    ]
    do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2], var_prepend="gen",
                          plot_dir=os.path.join(ROOT_DIR, "deltas_ptMin_gen_compare_reweight"),
                          zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR)
    do_angularity_delta_plots(sources, var_list=COMMON_VARS[:-2], var_prepend="gen",
                              plot_dir=os.path.join(ROOT_DIR, "deltas_angularities_gen_compare_reweight"),
                              zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR, pt_bins=THEORY_PT_BINS)

    # Put Pythia, Herwig, & Herwig reweighted on same plots
    sources = [
        {"root_dir": PYTHIA_AK4_DIR, 'label': "Pythia", "style": {'line_style': 1}},
        {"root_dir": HERWIG_AK4_REWEIGHTED_DIR, 'label': "Herwig, reweighted to Pythia",
            "style": {'line_style': 2},
            "dy_style": {'line_color': ROOT.kRed, 'fill_color': ROOT.kRed},
            "qcd_style": {'line_color': ROOT.kBlue, 'fill_color': ROOT.kBlue}
        },
        {"root_dir": HERWIG_AK4_DIR, 'label': "Herwig, not reweighted", "style": {'line_style': 2}}
    ]
    do_all_exclusive_plots_comparison(sources=sources, var_list=COMMON_VARS[:-2], var_prepend="gen",
                                      plot_dir=os.path.join(ROOT_DIR, "plots_dy_vs_qcd_gen_compare_reweight_compare_generators"),
                                      zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR,
                                      pt_bins=THEORY_PT_BINS, subplot_type="diff", do_flav_tagged=False)
    sources = [
        {"root_dir": PYTHIA_AK4_DIR, 'label': "Pythia", "style": {'line_color': ROOT.kBlack}},
        {"root_dir": HERWIG_AK4_REWEIGHTED_DIR, 'label': "Reweighted to Pythia", "style": {'line_color': ROOT.kGreen}},
        {"root_dir": HERWIG_AK4_DIR, 'label': "Not reweighted", "style": {'line_color': ROOT.kAzure}}
    ]
    do_pt_min_delta_plots(sources, var_list=COMMON_VARS[0:-2], var_prepend="gen",
                          plot_dir=os.path.join(ROOT_DIR, "deltas_ptMin_gen_compare_normal_and_reweight"),
                          zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR)
    do_angularity_delta_plots(sources, var_list=COMMON_VARS[:-2], var_prepend="gen",
                              plot_dir=os.path.join(ROOT_DIR, "deltas_angularities_gen_compare_normal_and_reweight"),
                              zpj_dirname=ZPJ_GENJET_RDIR, dj_dirname=DJ_GENJET_RDIR, pt_bins=THEORY_PT_BINS)


def do_pythia_comparison_plots():
    """To compare mg+pythia vs pythia only"""
    # reco
    sources = [
        {"root_dir": PYTHIA_AK4_DIR, 'label': "Madgraph+Pythia", "style": {'line_style': 1}},
        {"root_dir": PYTHIA_ONLY_AK4_DIR, 'label': "Pythia only", "style": {'line_style': 2, 'line_color': ROOT.kRed, 'fill_color': ROOT.kRed, }}
    ]
    do_all_exclusive_plots_comparison(sources, var_list=COMMON_VARS[:-2],
                                      plot_dir=os.path.join(ROOT_DIR, "mg_pythia_vs_pythia_only"),
                                      zpj_dirname="",
                                      subplot_type=None, do_flav_tagged=False,
                                      pt_bins=THEORY_PT_BINS)

    # gen
    sources = [
        {"root_dir": PYTHIA_AK4_DIR, 'label': "Madgraph+Pythia", "style": {'line_style': 1}},
        {"root_dir": PYTHIA_ONLY_AK4_DIR, 'label': "Pythia only", "style": {'line_style': 2, 'line_color': ROOT.kRed, 'fill_color': ROOT.kRed}}
    ]
    do_all_exclusive_plots_comparison(sources, var_list=COMMON_VARS[:-2], var_prepend="gen",
                                      plot_dir=os.path.join(ROOT_DIR, "mg_pythia_vs_pythia_only_gen"),
                                      dj_dirname=DJ_GENJET_RDIR, zpj_dirname="",
                                      subplot_type=None, do_flav_tagged=False,
                                      pt_bins=THEORY_PT_BINS)

    # flavour fractions
    input_files = [
        os.path.join(PYTHIA_AK4_DIR, 'uhh2.AnalysisModuleRunner.MC.MC_QCD_.root'),
        os.path.join(PYTHIA_ONLY_AK4_DIR, 'uhh2.AnalysisModuleRunner.MC.MC_QCD_.root')
    ]
    # reco
    compare_flavour_fractions_vs_pt(input_files,
                                    [DJ_RECOJET_RDIR, DJ_RECOJET_RDIR],
                                    [QCD_Dijet_LABEL + " Madgraph+Pythia", QCD_Dijet_LABEL+" Pythia only"],
                                    "g", "%s/flav_fractions/compare_g_frac.pdf" % ROOT_DIR,
                                    title="ak4 PFCHS jets", var_prepend="")
    # gen
    compare_flavour_fractions_vs_pt(input_files,
                                    [DJ_GENJET_RDIR, DJ_GENJET_RDIR],
                                    [QCD_Dijet_LABEL + " Madgraph+Pythia", QCD_Dijet_LABEL+" Pythia only"],
                                    "g", "%s/flav_fractions_gen/compare_g_frac.pdf" % ROOT_DIR,
                                    title="ak4 GenJets", var_prepend="gen")


if __name__ == '__main__':
    do_reco_plots()
    do_reco_flav_split_plots()
    do_reco_generator_comparison_plots()
    do_reco_reweight_comparison_plots()
    do_reco_pu_comparison_plots()

    do_gen_plots()
    do_gen_generator_comparison_plots()
    do_gen_reweight_comparison_plots()
    do_pythia_comparison_plots()
