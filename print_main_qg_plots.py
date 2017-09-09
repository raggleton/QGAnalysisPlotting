#!/usr/bin/env python

"""Print main QG plots for a given sample.

Any other plots comparing more than one sample should go into its own script!
"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os

# My stuff
from comparator import Contribution, Plot, grab_obj
import qg_common as qgc
import qg_general_plots as qgg
import qg_flavour_plots as qgf
import qg_delta_plots as qgd
import qg_roc_plots as qgr

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


SETUP = "ak8puppi"

PYTHIA_AK4_DIR = "workdir_%s_mgpythia" % SETUP
PYTHIA_ONLY_AK4_DIR = "workdir_%s_pythiaOnlyFlat" % SETUP
HERWIG_AK4_DIR = "workdir_%s_herwig" % SETUP
HERWIG_AK4_REWEIGHTED_DIR = "workdir_%s_herwig_reweight" % SETUP

# Controls all the things!!
ROOT_DIR = PYTHIA_AK4_DIR

# Control plot output format
OUTPUT_FMT = "pdf"


def do_all_2D_plots(root_dir, plot_dir="plots_2d", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG", var_list=None, var_prepend=""):
    """Do 2D distributions"""
    var_list = var_list or qgc.COMMON_VARS_WITH_FLAV
    for ang in var_list:
        v = "%s%s_vs_pt" % (var_prepend, ang.var)

        rebin = [2, 4]
        if v == "jet_multiplicity_vs_pt":
            rebin = [2, 4]
        elif v == "jet_thrust_vs_pt" or "flavour" in v:
            rebin = [1, 4]

        recolour = False if "flavour" in v else True

        for rn in ['Y', None]:  # different renormalisation axes
            qgg.do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % root_dir, "%s/%s" % (zpj_dirname, v)),
                           output_filename="%s/%s/dy_zpj_%s_norm%s.%s" % (root_dir, plot_dir, v, rn, OUTPUT_FMT),
                           renorm_axis=rn, title=qgc.DY_ZpJ_LABEL, rebin=rebin, recolour=recolour)
            qgg.do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % root_dir, "%s/%s" % (dj_dirname, v)),
                           output_filename="%s/%s/dy_dijet_%s_norm%s.%s" % (root_dir, plot_dir, v, rn, OUTPUT_FMT),
                           renorm_axis=rn, title=qgc.DY_Dijet_LABEL, rebin=rebin, recolour=recolour)

            qgg.do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % root_dir, "%s/%s" % (dj_dirname, v)),
                           output_filename="%s/%s/qcd_dijet_%s_norm%s.%s" % (root_dir, plot_dir, v, rn, OUTPUT_FMT),
                           renorm_axis=rn, title=qgc.QCD_Dijet_LABEL, rebin=rebin, recolour=recolour)

            if "flavour" in v:
                continue

            # Flavour matched reco
            qgg.do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % root_dir, "%s/q%s" % (zpj_dirname, v)),
                           output_filename="%s/%s/dy_zpj_%s_norm%s_qflavMatched.%s" % (root_dir, plot_dir, v, rn, OUTPUT_FMT),
                           renorm_axis=rn, title=qgc.DY_ZpJ_QFLAV_LABEL, rebin=rebin)

            qgg.do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % root_dir, "%s/g%s" % (zpj_dirname, v)),
                           output_filename="%s/%s/dy_zpj_%s_norm%s_gflavMatched.%s" % (root_dir, plot_dir, v, rn, OUTPUT_FMT),
                           renorm_axis=rn, title=qgc.DY_ZpJ_GFLAV_LABEL, rebin=rebin)

            qgg.do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % root_dir, "%s/g%s" % (dj_dirname, v)),
                           output_filename="%s/%s/dy_dijet_%s_norm%s_gflavMatched.%s" % (root_dir, plot_dir, v, rn, OUTPUT_FMT),
                           renorm_axis=rn, title=qgc.DY_Dijet_GFLAV_LABEL, rebin=rebin)

            qgg.do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % root_dir, "%s/q%s" % (dj_dirname, v)),
                           output_filename="%s/%s/dy_dijet_%s_norm%s_qflavMatched.%s" % (root_dir, plot_dir, v, rn, OUTPUT_FMT),
                           renorm_axis=rn, title=qgc.DY_Dijet_QFLAV_LABEL, rebin=rebin)

            qgg.do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % root_dir, "%s/g%s" % (dj_dirname, v)),
                           output_filename="%s/%s/qcd_dijet_%s_norm%s_gflavMatched.%s" % (root_dir, plot_dir, v, rn, OUTPUT_FMT),
                           renorm_axis=rn, title=qgc.QCD_Dijet_GFLAV_LABEL, rebin=rebin)

            qgg.do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % root_dir, "%s/q%s" % (dj_dirname, v)),
                           output_filename="%s/%s/qcd_dijet_%s_norm%s_qflavMatched.%s" % (root_dir, plot_dir, v, rn, OUTPUT_FMT),
                           renorm_axis=rn, title=qgc.QCD_Dijet_QFLAV_LABEL, rebin=rebin)


def do_all_flavour_fraction_plots(root_dir, plot_dir="flav_fractions", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG", var_prepend="", flav_source=""):
    """Do plots of jet flavour fractions vs pT, for both Z+jets and dijets regions"""
    
    # Z+jets
    qgf.do_flavour_fraction_vs_pt(input_file="%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % root_dir, 
                                  dirname=zpj_dirname, which=flav_source, var_prepend=var_prepend,
                                  output_filename="%s/%s/zpj.%s" % (root_dir, plot_dir, OUTPUT_FMT))

    # Dijets
    qgf.do_flavour_fraction_vs_pt(input_file="%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % root_dir, 
                                  dirname=dj_dirname, which=flav_source, var_prepend=var_prepend,
                                  output_filename="%s/%s/dj.%s" % (root_dir, plot_dir, OUTPUT_FMT)) 



def do_wrong_plots(root_dir, var_prepend="", plot_dir="wrong_flavs", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG", pt_bins=None):
    """Plot all the sample/selection.flavour combinations to check distributions indep of sample"""
    pt_bins = pt_bins or qgc.PT_BINS
    for v in ['jet_LHA', 'jet_pTD', 'jet_width', 'jet_thrust', 'jet_multiplicity']:
        v = "%s%s_vs_pt" % (var_prepend, v)

        h2d_dyj_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % root_dir, "%s/q%s" % (zpj_dirname, v))
        h2d_dyj_wrong_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % root_dir, "%s/g%s" % (zpj_dirname, v))
        h2d_dyj_qcd_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % root_dir, "%s/q%s" % (dj_dirname, v))
        h2d_dyj_qcd_wrong_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % root_dir, "%s/g%s" % (dj_dirname, v))
        h2d_qcd_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % root_dir, "%s/g%s" % (dj_dirname, v))
        h2d_qcd_wrong_chs = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % root_dir, "%s/q%s" % (dj_dirname, v))

        lw = 1
        dy_kwargs_chs = dict(line_color=qgc.DY_COLOUR, fill_color=qgc.DY_COLOUR, label=qgc.DY_ZpJ_QFLAV_LABEL, line_width=lw)
        dy_kwargs_wrong_chs = dict(line_color=qgc.DY_COLOUR+4, fill_color=qgc.DY_COLOUR+4, label=qgc.DY_ZpJ_GFLAV_LABEL, line_width=lw, line_style=1)
        dy_kwargs_qcd_chs = dict(line_color=ROOT.kGreen+2, fill_color=ROOT.kGreen+2, label=qgc.DY_Dijet_QFLAV_LABEL, line_width=lw, line_style=1)
        dy_kwargs_qcd_wrong_chs = dict(line_color=ROOT.kOrange-1, fill_color=ROOT.kOrange-1, label=qgc.DY_Dijet_GFLAV_LABEL, line_width=lw, line_style=1)
        qcd_kwargs_chs = dict(line_color=qgc.QCD_COLOUR, fill_color=qgc.QCD_COLOUR, label=qgc.QCD_Dijet_GFLAV_LABEL, line_width=lw)
        qcd_kwargs_wrong_chs = dict(line_color=ROOT.kRed, fill_color=ROOT.kRed, label=qgc.QCD_Dijet_QFLAV_LABEL, line_width=lw, line_style=1)

        rebin = 2
        xlim = None
        if "thrust" in v:
            rebin = 1
            xlim = (0, 0.5)

        for (start_val, end_val) in pt_bins:
            entries = [
                (qgg.get_projection_plot(h2d_dyj_chs, start_val, end_val), dy_kwargs_chs),
                (qgg.get_projection_plot(h2d_dyj_wrong_chs, start_val, end_val), dy_kwargs_wrong_chs),
                (qgg.get_projection_plot(h2d_dyj_qcd_chs, start_val, end_val), dy_kwargs_qcd_chs),
                (qgg.get_projection_plot(h2d_dyj_qcd_wrong_chs, start_val, end_val), dy_kwargs_qcd_wrong_chs),
                (qgg.get_projection_plot(h2d_qcd_wrong_chs, start_val, end_val), qcd_kwargs_wrong_chs),
                (qgg.get_projection_plot(h2d_qcd_chs, start_val, end_val), qcd_kwargs_chs),
            ]

            qgg.do_comparison_plot(entries, "%s/%s/%s_pt%dto%d_flavMatched.%s" % (root_dir, plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                   rebin=rebin, title="%d < p_{T}^{jet} < %d GeV (\"wrong\" flavours)" % (start_val, end_val), xlim=xlim)


def do_gen_reco_comparison_plots(root_dir, var_list=None, gen_var_prepend="gen", reco_var_prepend="",
                                 plot_dir="plot_reco_gen", zpj_reco_dirname=qgc.ZPJ_RECOJET_RDIR, dj_reco_dirname=qgc.DJ_RECOJET_RDIR,
                                 zpj_gen_dirname=qgc.ZPJ_GENJET_RDIR, dj_gen_dirname=qgc.DJ_GENJET_RDIR, pt_bins=qgc.THEORY_PT_BINS, subplot_type=None):
    var_list = var_list or qgc.COMMON_VARS[:-1]
    for ang in var_list:
        v_reco = "%s%s_vs_pt" % (reco_var_prepend, ang.var)
        v_gen = "%s%s_vs_pt" % (gen_var_prepend, ang.var)

        h2d_dyj_reco = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % root_dir, "%s/%s" % (zpj_reco_dirname, v_reco))
        h2d_qcd_reco = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % root_dir, "%s/%s" % (dj_reco_dirname, v_reco))
        h2d_dyj_gen = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % root_dir, "%s/%s" % (zpj_gen_dirname, v_gen))
        h2d_qcd_gen = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % root_dir, "%s/%s" % (dj_gen_dirname, v_gen))

        if "flavour" not in v_reco:
            h2d_dyj_reco_q = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % root_dir, "%s/q%s" % (zpj_reco_dirname, v_reco))
            h2d_qcd_reco_g = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % root_dir, "%s/g%s" % (dj_reco_dirname, v_reco))
            h2d_dyj_gen_q = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % root_dir, "%s/q%s" % (zpj_gen_dirname, v_gen))
            h2d_qcd_gen_g = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % root_dir, "%s/g%s" % (dj_gen_dirname, v_gen))


        for (start_val, end_val) in pt_bins:
            lw = 2

            dy_reco_kwargs = dict(line_color=qgc.DY_COLOUR, fill_color=qgc.DY_COLOUR, label=qgc.DY_ZpJ_LABEL + " [RecoJet]", line_width=lw)
            qcd_reco_kwargs = dict(line_color=qgc.QCD_COLOUR, fill_color=qgc.QCD_COLOUR, label=qgc.QCD_Dijet_LABEL + " [RecoJet]", line_width=lw)

            dy_gen_kwargs = dict(line_color=qgc.DY_COLOUR, fill_color=qgc.DY_COLOUR, label=qgc.DY_ZpJ_LABEL + " [GenJet]", line_width=lw, line_style=2)
            qcd_gen_kwargs = dict(line_color=qgc.QCD_COLOUR, fill_color=qgc.QCD_COLOUR, label=qgc.QCD_Dijet_LABEL + " [GenJet]", line_width=lw, line_style=2)


            entries = [
                (qgg.get_projection_plot(h2d_dyj_reco, start_val, end_val), dy_reco_kwargs),
                (qgg.get_projection_plot(h2d_qcd_reco, start_val, end_val), qcd_reco_kwargs),
                (qgg.get_projection_plot(h2d_dyj_gen, start_val, end_val), dy_gen_kwargs),
                (qgg.get_projection_plot(h2d_qcd_gen, start_val, end_val), qcd_gen_kwargs)
            ]

            rebin = 2
            xlim = None
            ylim = None
            if "flavour" in v_reco:
                rebin = 1
                ylim = (0, 1)
            if "thrust" in v_reco:
                xlim = (0, 0.5)

            qgg.do_comparison_plot(entries, "%s/%s/%s_pt%dto%d.%s" % (root_dir, plot_dir, ang.var, start_val, end_val, OUTPUT_FMT),
                                   rebin=rebin, title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val),
                                   xtitle=ang.name + " (" + ang.lambda_str + ")",
                                   xlim=xlim, ylim=ylim, subplot_type=subplot_type)

            # Do flavour-tagged comparison
            if "flavour" in v_reco:
                continue

            dy_reco_kwargs['label'] = qgc.DY_ZpJ_QFLAV_LABEL + " [RecoJet]"
            qcd_reco_kwargs['label'] = qgc.QCD_Dijet_GFLAV_LABEL + " [RecoJet]"
            dy_gen_kwargs['label'] = qgc.DY_ZpJ_QFLAV_LABEL + " [GenJet]"
            qcd_gen_kwargs['label'] = qgc.QCD_Dijet_GFLAV_LABEL + " [GenJet]"

            entries = [
                (qgg.get_projection_plot(h2d_dyj_reco_q, start_val, end_val), dy_reco_kwargs),
                (qgg.get_projection_plot(h2d_qcd_reco_g, start_val, end_val), qcd_reco_kwargs),
                (qgg.get_projection_plot(h2d_dyj_gen_q, start_val, end_val), dy_gen_kwargs),
                (qgg.get_projection_plot(h2d_qcd_gen_g, start_val, end_val), qcd_gen_kwargs)
            ]

            qgg.do_comparison_plot(entries, "%s/%s/%s_flavMatched_pt%dto%d.%s" % (root_dir, plot_dir, ang.var, start_val, end_val, OUTPUT_FMT),
                                   rebin=rebin, title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val),
                                   xtitle=ang.name + " (" + ang.lambda_str + ")",
                                   xlim=xlim, ylim=ylim, subplot_type=subplot_type)


def do_reco_plots(root_dir):
    do_all_2D_plots(root_dir)
    sources = [{"root_dir": root_dir, 'label': "", "style": {'line_style': 1}}]
    qgg.do_all_exclusive_plots_comparison(sources=sources, var_list=qgc.COMMON_VARS_WITH_FLAV,
                                          plot_dir=os.path.join(root_dir, "plots_dy_vs_qcd"),
                                          pt_bins=qgc.THEORY_PT_BINS, subplot_type=None, do_flav_tagged=True)
    do_all_flavour_fraction_plots(root_dir, plot_dir="flav_fractions", flav_source="genParton_")
    do_wrong_plots(root_dir)
    do_reco_pu_comparison_plots(root_dir)

    # Separation plots
    qgd.do_angularity_delta_plots(sources, var_list=qgc.COMMON_VARS, pt_bins=qgc.THEORY_PT_BINS, 
                                  plot_dir=os.path.join(root_dir, 'delta_angularities'),
                                  save_component_hists=True)
    qgr.do_angularity_roc_plots(sources, var_list=qgc.COMMON_VARS, pt_bins=qgc.THEORY_PT_BINS,
                                plot_dir=os.path.join(root_dir, 'roc_angularities'))
    


def do_reco_pu_comparison_plots(root_dir):
    """Compare by PU bins"""
    pu_bins = [(5, 15), (20, 25), (30, 40)]
    sources = []
    for ind, (pu_min, pu_max) in enumerate(pu_bins):
        sources.append({
            "root_dir": root_dir,
            'label': "PU %d-%d" % (pu_min, pu_max),
            'zpj_dirname': qgc.ZPJ_RECOJET_RDIR + "_PU_%d_to_%d" % (pu_min, pu_max),
            'dj_dirname': qgc.DJ_RECOJET_RDIR + "_PU_%d_to_%d" % (pu_min, pu_max),
            "style": {'line_style': 1, "line_width": 1},
            "dy_style": {'line_color': qgc.DY_COLOURS[ind], 'fill_color': qgc.DY_COLOURS[ind]},
            "qcd_style": {'line_color': qgc.QCD_COLOURS[ind], 'fill_color': qgc.QCD_COLOURS[ind]}
        })
    qgg.do_all_exclusive_plots_comparison(sources=sources, var_list=qgc.COMMON_VARS, zpj_dirname=None,
                                          plot_dir=os.path.join(root_dir, "plots_dy_vs_qcd_compare_pu_dijet"),
                                          pt_bins=qgc.THEORY_PT_BINS, subplot_type="ratio", do_flav_tagged=False)
    qgg.do_all_exclusive_plots_comparison(sources=sources, var_list=qgc.COMMON_VARS, dj_dirname=None,
                                          plot_dir=os.path.join(root_dir, "plots_dy_vs_qcd_compare_pu_zpj"),
                                          pt_bins=qgc.THEORY_PT_BINS, subplot_type="ratio", do_flav_tagged=False)

    # Separation plots
    for ind, s in enumerate(sources):
        sources[ind]['style']['line_width'] = 2
        sources[ind]['style']['line_color'] = qgc.DY_COLOURS[ind]
        if ind == 2:
            qgd.do_angularity_delta_plots(sources[ind:ind+1], 
                                          var_list=qgc.COMMON_VARS,
                                          plot_dir=os.path.join(root_dir, "deltas_angularities_compare_pu_PU_%d_to_%d" % (pu_bins[ind][0], pu_bins[ind][1])),
                                          pt_bins=qgc.THEORY_PT_BINS, save_component_hists=True)

    qgd.do_angularity_delta_plots(sources, var_list=qgc.COMMON_VARS,
                                  plot_dir=os.path.join(root_dir, "deltas_angularities_compare_pu"),
                                  pt_bins=qgc.THEORY_PT_BINS)


def do_gen_plots(root_dir):
    # need to avoid genPArtonFlavour
    do_all_2D_plots(root_dir, var_list=qgc.COMMON_VARS_WITH_FLAV[:-1], var_prepend="gen", plot_dir="plots_2d_gen",
                    zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR)
    sources = [{"root_dir": root_dir, 'label': "", "style": {'line_style': 1}}]
    qgg.do_all_exclusive_plots_comparison(sources=sources, var_list=qgc.COMMON_VARS_WITH_FLAV[:-1], var_prepend="gen",
                                          plot_dir=os.path.join(root_dir, "plots_dy_vs_qcd_gen"),
                                          zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR,
                                          pt_bins=qgc.THEORY_PT_BINS, subplot_type=None, do_flav_tagged=True)
    do_all_flavour_fraction_plots(root_dir, var_prepend="gen", plot_dir="flav_fractions_gen", flav_source="",
                                  zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR)
    do_wrong_plots(root_dir, var_prepend="gen", plot_dir="wrong_flavs_gen",
                   zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR, pt_bins=qgc.THEORY_PT_BINS)
    # do_gen_reco_comparison_plots(var_list=qgc.COMMON_VARS[:-1], gen_var_prepend="gen", reco_var_prepend="",
    #                              plot_dir="plot_reco_gen", zpj_reco_dirname=qgc.ZPJ_RECOJET_RDIR, dj_reco_dirname=qgc.DJ_RECOJET_RDIR,
    #                              zpj_gen_dirname=qgc.ZPJ_GENJET_RDIR, dj_gen_dirname=qgc.DJ_GENJET_RDIR, pt_bins=qgc.THEORY_PT_BINS)

    qgr.do_angularity_roc_plots(sources, var_list=qgc.COMMON_VARS, pt_bins=qgc.THEORY_PT_BINS,
                                zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR, var_prepend="gen",
                                plot_dir=os.path.join(root_dir, 'roc_angularities_gen'))

    # Separation plots
    qgd.do_pt_min_delta_plots(sources, var_list=qgc.COMMON_VARS, var_prepend="gen",
                              plot_dir=os.path.join(root_dir, "deltas_ptMin_gen"),
                              zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR, save_component_hists=True)
    qgd.do_angularity_delta_plots(sources, var_list=qgc.COMMON_VARS, var_prepend="gen",
                                  plot_dir=os.path.join(root_dir, "deltas_angularities_gen"),
                                  zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR,
                                  pt_bins=qgc.THEORY_PT_BINS, save_component_hists=True)
    # flav-tagged versions
    if root_dir == HERWIG_AK4_DIR:
      return
    qgd.do_pt_min_delta_plots(sources, var_list=qgc.COMMON_VARS, var_prepend="gen",
                              plot_dir=os.path.join(root_dir, "deltas_ptMin_gen"),
                              zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR,
                              flavour_tag=True, save_component_hists=True)
    qgd.do_angularity_delta_plots(sources, var_list=qgc.COMMON_VARS, var_prepend="gen",
                                  plot_dir=os.path.join(root_dir, "deltas_angularities_gen"),
                                  zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR,
                                  pt_bins=qgc.THEORY_PT_BINS, flavour_tag=True, save_component_hists=True)


if __name__ == '__main__':
    do_reco_plots(ROOT_DIR)

    do_gen_plots(ROOT_DIR)
