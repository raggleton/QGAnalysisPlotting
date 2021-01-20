#!/usr/bin/env python

"""Print main lambda variable plots, comparing computation with/without grooming definition"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
os.nice(10)

# My stuff
from comparator import Contribution, Plot
import qg_common as qgc
import qg_general_plots as qgp
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


if __name__ == "__main__":
    # pythia_dir = "workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_noZReweight_wta_groomed"
    # pythia_dir = "workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_noZReweight_wta_groomed_fwdcenDijet"
    # pythia_dir = "workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_noZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto"

    pythia_dir = "workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet"
    data_dir = "workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet"

    # pythia_dir = "workdir_ak8chs_pythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_noZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto"
    # pythia_dir = "workdir_ak8puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet"
    # data_dir = "workdir_ak8puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet"

    mc_sources = [
        {"root_dir": pythia_dir ,
         "label": "MG+PYTHIA8 MC (ungroomed)",
         "style": {
             "line_color": ROOT.kRed,
             "marker_color": ROOT.kRed
         },
         "dy_filename": qgc.DY_FILENAME,
         "zpj_dirname": "ZPlusJets_QG",
         "qcd_filename": qgc.QCD_FILENAME,
         "dj_cen_dirname": "Dijet_QG_central_tighter",
         "dj_fwd_dirname": "Dijet_QG_forward_tighter"
        },

        {"root_dir": pythia_dir,
         "label": "MG+PYTHIA8 MC (groomed)",
         "style": {
             "line_color": ROOT.kBlue,
             "marker_color": ROOT.kBlue
         },
         "dy_filename": qgc.DY_FILENAME,
         "zpj_dirname": "ZPlusJets_QG_groomed",
         "qcd_filename": qgc.QCD_FILENAME,
         "dj_cen_dirname": "Dijet_QG_central_tighter_groomed",
         "dj_fwd_dirname": "Dijet_QG_forward_tighter_groomed"
        },
    ]
    data_sources = [
        {"root_dir": data_dir ,
         "label": "Data (ungroomed)",
         "style": {
             "line_color": ROOT.kBlack,
             "marker_color": ROOT.kBlack
         },
         "dy_filename": qgc.SINGLE_MU_FILENAME,
         "zpj_dirname": "ZPlusJets_QG",
         "qcd_filename": qgc.JETHT_ZB_FILENAME,
         "dj_cen_dirname": "Dijet_QG_central_tighter",
         "dj_fwd_dirname": "Dijet_QG_forward_tighter"
        },

        {"root_dir": data_dir,
         "label": "Data (groomed)",
         "style": {
             "line_color": ROOT.kGreen+2,
             "marker_color": ROOT.kGreen+2
         },
         "dy_filename": qgc.SINGLE_MU_FILENAME,
         "zpj_dirname": "ZPlusJets_QG_groomed",
         "qcd_filename": qgc.JETHT_ZB_FILENAME,
         "dj_cen_dirname": "Dijet_QG_central_tighter_groomed",
         "dj_fwd_dirname": "Dijet_QG_forward_tighter_groomed"
        },
    ]

    all_sources = [
        data_sources[0],
        mc_sources[0],
        data_sources[1],
        mc_sources[1],
        # {"root_dir": data_dir ,
        #  "label": "Data (ungroomed)",
        #  "style": {
        #      "line_color": ROOT.kBlack,
        #      "marker_color": ROOT.kBlack
        #  },
        #  "dy_filename": qgc.SINGLE_MU_FILENAME,
        #  "zpj_dirname": "ZPlusJets_QG",
        #  "qcd_filename": qgc.JETHT_ZB_FILENAME,
        #  "dj_cen_dirname": "Dijet_QG_central_tighter",
        #  "dj_fwd_dirname": "Dijet_QG_forward_tighter"
        # },
        # {"root_dir": pythia_dir ,
        #  "label": "MG+PYTHIA8 MC (ungroomed)",
        #  "style": {
        #      "line_color": ROOT.kRed,
        #      "marker_color": ROOT.kRed
        #  },
        #  "dy_filename": qgc.DY_FILENAME,
        #  "zpj_dirname": "ZPlusJets_QG",
        #  "qcd_filename": qgc.QCD_FILENAME,
        #  "dj_cen_dirname": "Dijet_QG_central_tighter",
        #  "dj_fwd_dirname": "Dijet_QG_forward_tighter"
        # },

        # {"root_dir": data_dir,
        #  "label": "Data (groomed)",
        #  "style": {
        #      "line_color": ROOT.kGreen+2,
        #      "marker_color": ROOT.kGreen+2
        #  },
        #  "dy_filename": qgc.SINGLE_MU_FILENAME,
        #  "zpj_dirname": "ZPlusJets_QG_groomed",
        #  "qcd_filename": qgc.JETHT_ZB_FILENAME,
        #  "dj_cen_dirname": "Dijet_QG_central_tighter_groomed",
        #  "dj_fwd_dirname": "Dijet_QG_forward_tighter_groomed"
        # },
        # {"root_dir": pythia_dir,
        #  "label": "MG+PYTHIA8 MC (groomed)",
        #  "style": {
        #      "line_color": ROOT.kBlue,
        #      "marker_color": ROOT.kBlue
        #  },
        #  "dy_filename": qgc.DY_FILENAME,
        #  "zpj_dirname": "ZPlusJets_QG_groomed",
        #  "qcd_filename": qgc.QCD_FILENAME,
        #  "dj_cen_dirname": "Dijet_QG_central_tighter_groomed",
        #  "dj_fwd_dirname": "Dijet_QG_forward_tighter_groomed"
        # },
    ]

    title = qgc.extract_jet_config(pythia_dir)
    # Do Z+jets region only
    # print("Doing Z+J ungroomed")
    # qgp.do_all_exclusive_plots_comparison(sources=all_sources[:2],
    #                                       var_list=qgc.COMMON_VARS[:],
    #                                       plot_dir=os.path.join(pythia_dir, "ungroomed_jet_zpj"),
    #                                       dy_filename=qgc.DY_FILENAME,
    #                                       qcd_filename=None,
    #                                       dj_cen_dirname=None,
    #                                       dj_fwd_dirname=None,
    #                                       show_region_labels=False,
    #                                       subplot_type="ratio",
    #                                       # subplot_title="#splitline{Groomed/}{ungroomed}",
    #                                       do_flav_tagged=False,
    #                                       pt_bins=qgc.PT_BINS,
    #                                       title=title+"\n"+qgc.ZpJ_LABEL,
    #                                       has_data=True,
    #                                       ofmt=OUTPUT_FMT)
    print("Doing Z+J groomed vs ungroomed")
    all_sources[2]['subplot'] = 0
    all_sources[3]['subplot'] = 1
    qgp.do_all_exclusive_plots_comparison(sources=all_sources,
                                          var_list=qgc.COMMON_VARS[:],
                                          plot_dir=os.path.join(pythia_dir, "groomed_vs_ungroomed_jet_zpj"),
                                          dy_filename=qgc.DY_FILENAME,
                                          qcd_filename=None,
                                          dj_cen_dirname=None,
                                          dj_fwd_dirname=None,
                                          show_region_labels=False,
                                          subplot_type="ratio",
                                          subplot_title="#splitline{Groomed/}{ungroomed}",
                                          do_flav_tagged=False,
                                          pt_bins=qgc.PT_BINS,
                                          title=title+"\n"+qgc.ZpJ_LABEL,
                                          has_data=True,
                                          ofmt=OUTPUT_FMT)

    # Do Z+jet and dijet regions altogether
    # qgp.do_all_exclusive_plots_comparison(sources=sources,
    #                                       var_list=qgc.COMMON_VARS[:],
    #                                       plot_dir=os.path.join(pythia_dir, "groomed_vs_ungroomed_jet_dijet_zpj"),
    #                                       subplot_type="ratio",
    #                                       subplot_title="#splitline{Groomed/}{ungroomed}",
    #                                       do_flav_tagged=False,
    #                                       pt_bins=qgc.PT_BINS,
    #                                       title=title,
    #                                       has_data=False,
    #                                       ofmt=OUTPUT_FMT)

    # qgp.do_all_exclusive_plots_comparison(sources=sources[0:1]+data_sources[0:1],
    #                                       var_list=qgc.COMMON_VARS[:],
    #                                       plot_dir=os.path.join(pythia_dir, "ungroomed_jet_dijet_zpj"),
    #                                       # subplot_type="ratio",
    #                                       # subplot_title="#splitline{Groomed/}{ungroomed}",
    #                                       do_flav_tagged=False,
    #                                       pt_bins=qgc.PT_BINS,
    #                                       title=title,
    #                                       has_data=False,
    #                                       ofmt=OUTPUT_FMT)

    # qgp.do_all_exclusive_plots_comparison(sources=sources[1:2]+data_sources[1:2],
    #                                       var_list=qgc.COMMON_VARS[:],
    #                                       plot_dir=os.path.join(pythia_dir, "groomed_jet_dijet_zpj"),
    #                                       # subplot_type="ratio",
    #                                       # subplot_title="#splitline{Groomed/}{ungroomed}",
    #                                       do_flav_tagged=False,
    #                                       pt_bins=qgc.PT_BINS,
    #                                       title=title,
    #                                       has_data=False,
    #                                       ofmt=OUTPUT_FMT)

    # Do Dijet region (central + forward)
    print("Doing dijet central + forward ungroomed")
    qgp.do_all_exclusive_plots_comparison(sources=data_sources[0:1]+mc_sources[0:1],
                                          var_list=qgc.COMMON_VARS[:],
                                          plot_dir=os.path.join(pythia_dir, "ungroomed_jet_dijet"),
                                          dy_filename=None,
                                          zpj_dirname=None,
                                          do_flav_tagged=False,
                                          pt_bins=qgc.PT_BINS,
                                          title=title,
                                          has_data=True,
                                          ofmt=OUTPUT_FMT)

    print("Doing dijet central + forward groomed")
    qgp.do_all_exclusive_plots_comparison(sources=data_sources[1:]+mc_sources[1:],
                                          var_list=qgc.COMMON_VARS[:],
                                          plot_dir=os.path.join(pythia_dir, "groomed_jet_dijet"),
                                          dy_filename=None,
                                          zpj_dirname=None,
                                          do_flav_tagged=False,
                                          pt_bins=qgc.PT_BINS,
                                          title=title,
                                          has_data=True,
                                          ofmt=OUTPUT_FMT)

    # qgp.do_all_exclusive_plots_comparison(sources=all_sources,
    #                                       var_list=qgc.COMMON_VARS[:],
    #                                       plot_dir=os.path.join(pythia_dir, "groomed_vs_ungroomed_jet_dijet"),
    #                                       dy_filename=None,
    #                                       zpj_dirname=None,
    #                                       qcd_filename=qgc.QCD_FILENAME,
    #                                       subplot_type="ratio",
    #                                       subplot_title="#splitline{Groomed/}{ungroomed}",
    #                                       do_flav_tagged=False,
    #                                       pt_bins=qgc.PT_BINS,
    #                                       title=title,
    #                                       has_data=True,
    #                                       ofmt=OUTPUT_FMT)

    # Do Dijet central region only
    print("Doing dijet central groomed vs ungroomed")
    qgp.do_all_exclusive_plots_comparison(sources=all_sources,
                                          var_list=qgc.COMMON_VARS[:],
                                          plot_dir=os.path.join(pythia_dir, "groomed_vs_ungroomed_jet_dijet_central"),
                                          zpj_dirname=None,
                                          qcd_filename=qgc.QCD_FILENAME,
                                          dj_fwd_dirname=None,
                                          subplot_type="ratio",
                                          subplot_title="#splitline{Groomed/}{ungroomed}",
                                          do_flav_tagged=False,
                                          pt_bins=qgc.PT_BINS,
                                          title=title+"\n"+qgc.Dijet_CEN_LABEL,
                                          has_data=True,
                                          show_region_labels=False,
                                          ofmt=OUTPUT_FMT)

    # Do Dijet forward region only
    print("Doing dijet forward groomed vs ungroomed")
    qgp.do_all_exclusive_plots_comparison(sources=all_sources,
                                          var_list=qgc.COMMON_VARS[:],
                                          plot_dir=os.path.join(pythia_dir, "groomed_vs_ungroomed_jet_dijet_forward"),
                                          zpj_dirname=None,
                                          qcd_filename=qgc.QCD_FILENAME,
                                          dj_cen_dirname=None,
                                          subplot_type="ratio",
                                          subplot_title="#splitline{Groomed/}{ungroomed}",
                                          do_flav_tagged=False,
                                          pt_bins=qgc.PT_BINS,
                                          title=title+"\n"+qgc.Dijet_FWD_LABEL,
                                          has_data=True,
                                          show_region_labels=False,
                                          ofmt=OUTPUT_FMT)

    # Do Pileup comparison
    pu_bins = [(5, 15), (20, 25), (30, 40)]
    ungroomed_sources = []
    groomed_sources = []
    for ind, (pu_min, pu_max) in enumerate(pu_bins):
        ungroomed_sources.append({
            "root_dir": pythia_dir,
            "label": "PU %d-%d" % (pu_min, pu_max),
            "zpj_dirname": "ZPlusJets_QG_PU_%d_to_%d" % (pu_min, pu_max),
            "dj_cen_dirname": "Dijet_QG_central_tighter_PU_%d_to_%d" % (pu_min, pu_max),
            "dj_fwd_dirname": "Dijet_QG_forward_tighter_PU_%d_to_%d" % (pu_min, pu_max),
            "style": {"line_style": 1, "line_width": 2},
            "zpj_style": {"line_color": qgc.DY_COLOURS[ind], "fill_color": qgc.DY_COLOURS[ind], "marker_color": qgc.DY_COLOURS[ind]},
            "qcd_cen_style": {"line_color": qgc.QCD_COLOURS[ind], "fill_color": qgc.QCD_COLOURS[ind], "marker_color": qgc.QCD_COLOURS[ind]},
            "qcd_fwd_style": {"line_color": qgc.QCD_COLOURS[ind], "fill_color": qgc.QCD_COLOURS[ind], "marker_color": qgc.QCD_COLOURS[ind], "line_style": 2},
        })
        groomed_sources.append({
            "root_dir": pythia_dir,
            "label": "PU %d-%d (groomed)" % (pu_min, pu_max),
            "zpj_dirname": "ZPlusJets_QG_PU_%d_to_%d_groomed" % (pu_min, pu_max),
            "dj_cen_dirname": "Dijet_QG_central_tighter_PU_%d_to_%d_groomed" % (pu_min, pu_max),
            "dj_fwd_dirname": "Dijet_QG_forward_tighter_PU_%d_to_%d_groomed" % (pu_min, pu_max),
            "style": {"line_style": 1, "line_width": 2},
            "zpj_style": {"line_color": qgc.DY_COLOURS[ind], "fill_color": qgc.DY_COLOURS[ind], "marker_color": qgc.DY_COLOURS[ind]},
            "qcd_cen_style": {"line_color": qgc.QCD_COLOURS[ind], "fill_color": qgc.QCD_COLOURS[ind], "marker_color": qgc.QCD_COLOURS[ind]},
            "qcd_fwd_style": {"line_color": qgc.QCD_COLOURS[ind], "fill_color": qgc.QCD_COLOURS[ind], "marker_color": qgc.QCD_COLOURS[ind], "line_style": 2},
        })
    subplot_title = "#splitline{Ratio wrt}{PU %d-%d}" % (pu_bins[0][0], pu_bins[0][1])

    # UNGROOMED
    # Z+jets plots
    print("Doing Z+J PU binned ungroomed")
    qgp.do_all_exclusive_plots_comparison(sources=ungroomed_sources,
                                          var_list=qgc.COMMON_VARS,
                                          dj_cen_dirname=None,
                                          dj_fwd_dirname=None,
                                          plot_dir=os.path.join(pythia_dir, "plots_dy_vs_qcd_compare_pu_zpj"),
                                          pt_bins=qgc.PT_BINS,
                                          subplot_type="ratio",
                                          subplot_title=subplot_title,
                                          title=title+"\n"+qgc.ZpJ_LABEL,
                                          has_data=False,
                                          show_region_labels=False,
                                          do_flav_tagged=False)
    # Dijet central
    print("Doing dijet central PU binned ungroomed")
    qgp.do_all_exclusive_plots_comparison(sources=ungroomed_sources,
                                          var_list=qgc.COMMON_VARS,
                                          zpj_dirname=None,
                                          dj_fwd_dirname=None,
                                          plot_dir=os.path.join(pythia_dir, "plots_dy_vs_qcd_compare_pu_dijet_central"),
                                          # qcd_filename=qgc.QCD_PYTHIA_ONLY_FILENAME,
                                          pt_bins=qgc.PT_BINS,
                                          subplot_type="ratio",
                                          subplot_title=subplot_title,
                                          title=title+"\n"+qgc.Dijet_CEN_LABEL,
                                          has_data=False,
                                          show_region_labels=False,
                                          do_flav_tagged=False)

    # Dijet forward
    print("Doing dijet forward PU binned ungroomed")
    qgp.do_all_exclusive_plots_comparison(sources=ungroomed_sources,
                                          var_list=qgc.COMMON_VARS,
                                          zpj_dirname=None,
                                          dj_cen_dirname=None,
                                          plot_dir=os.path.join(pythia_dir, "plots_dy_vs_qcd_compare_pu_dijet_forward"),
                                          # qcd_filename=qgc.QCD_PYTHIA_ONLY_FILENAME,
                                          pt_bins=qgc.PT_BINS,
                                          subplot_type="ratio",
                                          subplot_title=subplot_title,
                                          title=title+"\n"+qgc.Dijet_FWD_LABEL,
                                          has_data=False,
                                          show_region_labels=False,
                                          do_flav_tagged=False)

    # GROOMED
    # Z+jets plots
    print("Doing Z+J PU binned groomed")
    qgp.do_all_exclusive_plots_comparison(sources=groomed_sources,
                                          var_list=qgc.COMMON_VARS,
                                          dj_cen_dirname=None,
                                          dj_fwd_dirname=None,
                                          plot_dir=os.path.join(pythia_dir, "plots_dy_vs_qcd_compare_pu_zpj_groomed"),
                                          pt_bins=qgc.PT_BINS,
                                          subplot_type="ratio",
                                          subplot_title=subplot_title,
                                          title=title+"\n"+qgc.ZpJ_LABEL,
                                          has_data=False,
                                          show_region_labels=False,
                                          do_flav_tagged=False)
    # Dijet central
    print("Doing dijet central PU binned groomed")
    qgp.do_all_exclusive_plots_comparison(sources=groomed_sources,
                                          var_list=qgc.COMMON_VARS,
                                          zpj_dirname=None,
                                          dj_fwd_dirname=None,
                                          plot_dir=os.path.join(pythia_dir, "plots_dy_vs_qcd_compare_pu_dijet_central_groomed"),
                                          # qcd_filename=qgc.QCD_PYTHIA_ONLY_FILENAME,
                                          pt_bins=qgc.PT_BINS,
                                          subplot_type="ratio",
                                          subplot_title=subplot_title,
                                          title=title+"\n"+qgc.Dijet_CEN_LABEL,
                                          has_data=False,
                                          show_region_labels=False,
                                          do_flav_tagged=False)

    # Dijet forward
    print("Doing dijet forward PU binned groomed")
    qgp.do_all_exclusive_plots_comparison(sources=groomed_sources,
                                          var_list=qgc.COMMON_VARS,
                                          zpj_dirname=None,
                                          dj_cen_dirname=None,
                                          plot_dir=os.path.join(pythia_dir, "plots_dy_vs_qcd_compare_pu_dijet_forward_groomed"),
                                          # qcd_filename=qgc.QCD_PYTHIA_ONLY_FILENAME,
                                          pt_bins=qgc.PT_BINS,
                                          subplot_type="ratio",
                                          subplot_title=subplot_title,
                                          title=title+"\n"+qgc.Dijet_FWD_LABEL,
                                          has_data=False,
                                          show_region_labels=False,
                                          do_flav_tagged=False)

