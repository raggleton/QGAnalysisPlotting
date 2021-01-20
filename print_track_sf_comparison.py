#!/usr/bin/env python

"""Print main QG plots comparing with/without track SF, and impact of variations.

NB all at reco level
"""

from __future__ import print_function

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

    NO_SF_DIR = "workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet_Zreweight_noZjet2Cut_newBinningFixed"
    NOMINAL_SF_DIR = "workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet_Zreweight_noZjet2Cut_newBinningFixed_trkSFMyFormula"
    SF_UP_DIR = "workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet_Zreweight_noZjet2Cut_newBinningFixed_trkSFMyFormula/systematics_files/track_directionUp"
    SF_DOWN_DIR = "workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet_Zreweight_noZjet2Cut_newBinningFixed_trkSFMyFormula/systematics_files/track_directionDown"

    QCD_NO_SF_TFILE = cu.TFileCacher(os.path.join(NO_SF_DIR, qgc.QCD_FILENAME))
    QCD_NOMINAL_SF_TFILE = cu.TFileCacher(os.path.join(NOMINAL_SF_DIR, qgc.QCD_FILENAME))
    QCD_SF_UP_TFILE = cu.TFileCacher(os.path.join(SF_UP_DIR, qgc.QCD_FILENAME))
    QCD_SF_DOWN_TFILE = cu.TFileCacher(os.path.join(SF_DOWN_DIR, qgc.QCD_FILENAME))

    DY_NO_SF_TFILE = cu.TFileCacher(os.path.join(NO_SF_DIR, qgc.DY_FILENAME))
    DY_NOMINAL_SF_TFILE = cu.TFileCacher(os.path.join(NOMINAL_SF_DIR, qgc.DY_FILENAME))
    DY_SF_UP_TFILE = cu.TFileCacher(os.path.join(SF_UP_DIR, qgc.DY_FILENAME))
    DY_SF_DOWN_TFILE = cu.TFileCacher(os.path.join(SF_DOWN_DIR, qgc.DY_FILENAME))

    # QG variable plots
    pt_bins = qgc.PT_BINS[:]
    var_list = qgc.COMMON_VARS
    var_prepend = ""

    radius, pus = cu.get_jet_config_from_dirname(NO_SF_DIR)
    jet_str = "AK%s PF %s" % (radius.upper(), pus.upper())

    for gr_append in ["", "_groomed"]:
        if gr_append == "_groomed":
            print("Doing groomed plots...")
        else:
            print("Doing ungroomed plots...")

        zpj_dirname = "ZPlusJets_QG%s" % (gr_append)
        dj_cen_dirname = "Dijet_QG_central_tighter%s" % (gr_append)
        dj_fwd_dirname = "Dijet_QG_forward_tighter%s" % (gr_append)

        for ang in var_list[:]:
            print("...Doing", ang.name)
            v = "%s%s_vs_pt" % (var_prepend, ang.var)

            zpj_1d_entries = []
            dijet_cen_1d_entries = []
            dijet_fwd_1d_entries = []

            for pt_ind, (start_val, end_val) in enumerate(pt_bins):
                zpj_no_vs_nominal_sf_entries = []
                zpj_nominal_vs_syst_sf_entries = []
                dijet_cen_no_vs_nominal_sf_entries = []
                dijet_cen_nominal_vs_syst_sf_entries = []
                dijet_fwd_no_vs_nominal_sf_entries = []
                dijet_fwd_nominal_vs_syst_sf_entries = []

                # Get all plots
                lw = 2
                msize = 1.1
                data_line_width = lw

                ####################
                # Z+JETS REGION
                ####################
                zpj_obj_name = "%s/%s" % (zpj_dirname, v)

                # NO SF
                zpj_no_sf_hist = qgp.get_projection_plot(DY_NO_SF_TFILE[zpj_obj_name], start_val, end_val)
                dy_kwargs_no_sf = dict(line_color=qgc.SINGLE_MU_COLOUR, line_width=data_line_width, fill_color=qgc.SINGLE_MU_COLOUR,
                                       marker_color=qgc.SINGLE_MU_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=0,
                                       label="No TRK SF")
                zpj_no_vs_nominal_sf_entries.append((zpj_no_sf_hist, dy_kwargs_no_sf))

                # NOMINAL TK SF
                zpj_nominal_sf_hist = qgp.get_projection_plot(DY_NOMINAL_SF_TFILE[zpj_obj_name], start_val, end_val)
                dy_kwargs_nominal_sf = dict(line_color=qgc.DY_COLOUR, line_width=lw, fill_color=qgc.DY_COLOUR,
                                            marker_color=qgc.DY_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=0,
                                            label="Nominal TRK SF",
                                            subplot=zpj_no_sf_hist)
                zpj_no_vs_nominal_sf_entries.append((zpj_nominal_sf_hist, dy_kwargs_nominal_sf))

                dy_kwargs_nominal_sf2 = dict(line_color=qgc.DY_COLOUR, line_width=lw, fill_color=qgc.DY_COLOUR,
                                            marker_color=qgc.DY_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=0,
                                            label="Nominal TRK SF")
                zpj_nominal_vs_syst_sf_entries.append((zpj_nominal_sf_hist.Clone(), dy_kwargs_nominal_sf2))  # clone since we don't want to rebin twice

                # TK SF UP
                dy_col_up = qgc.DY_COLOURS[2]
                dy_kwargs_sf_up = dict(line_color=dy_col_up, line_width=lw, fill_color=dy_col_up,
                                       marker_color=dy_col_up, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=0,
                                       label="TRK SF Up",
                                       subplot=zpj_nominal_sf_hist)
                zpj_sf_up_hist = qgp.get_projection_plot(DY_SF_UP_TFILE[zpj_obj_name], start_val, end_val)
                zpj_nominal_vs_syst_sf_entries.append((zpj_sf_up_hist, dy_kwargs_sf_up))

                # TK SF DOWN
                dy_col_down = qgc.DY_COLOURS[3]
                dy_kwargs_sf_down = dict(line_color=dy_col_down, line_width=lw, fill_color=dy_col_down,
                                         marker_color=dy_col_down, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=0,
                                         label="TRK SF Down",
                                         subplot=zpj_nominal_sf_hist)
                zpj_sf_down_hist = qgp.get_projection_plot(DY_SF_DOWN_TFILE[zpj_obj_name], start_val, end_val)
                zpj_nominal_vs_syst_sf_entries.append((zpj_sf_down_hist, dy_kwargs_sf_down))

                zpj_1d_entries.append([
                    (zpj_no_sf_hist, dy_kwargs_no_sf),
                    (zpj_nominal_sf_hist, dy_kwargs_nominal_sf),
                    (zpj_sf_up_hist, dy_kwargs_sf_up),
                    (zpj_sf_down_hist, dy_kwargs_sf_down),
                ])

                ####################
                # DIJET CENTRAL REGION
                ####################
                cen_obj_name = "%s/%s" % (dj_cen_dirname, v)

                # NO TK SF
                dijet_cen_no_sf_hist = qgp.get_projection_plot(QCD_NO_SF_TFILE[cen_obj_name], start_val, end_val)
                qcd_cen_kwargs_no_sf = dict(line_color=qgc.JETHT_COLOUR, line_width=data_line_width, fill_color=qgc.JETHT_COLOUR,
                                            marker_color=qgc.JETHT_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                            label="No TRK SF")
                dijet_cen_no_vs_nominal_sf_entries.append((dijet_cen_no_sf_hist, qcd_cen_kwargs_no_sf))

                # NOMINAL TK SF
                dijet_cen_nominal_sf_hist = qgp.get_projection_plot(QCD_NOMINAL_SF_TFILE[cen_obj_name], start_val, end_val)

                qcd_cen_kwargs_mc = dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                                         marker_color=qgc.QCD_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                         label="Nominal TRK SF",
                                         subplot=dijet_cen_no_sf_hist)
                dijet_cen_no_vs_nominal_sf_entries.append((dijet_cen_nominal_sf_hist, qcd_cen_kwargs_mc))

                qcd_cen_kwargs_mc2 = dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                                          marker_color=qgc.QCD_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                          label="Nominal TRK SF")
                dijet_cen_nominal_vs_syst_sf_entries.append((dijet_cen_nominal_sf_hist.Clone(), qcd_cen_kwargs_mc2))

                # TK SF UP
                col_up = qgc.QCD_COLOURS[2]
                qcd_cen_kwargs_sf_up = dict(line_color=col_up, line_width=lw, fill_color=col_up,
                                            marker_color=col_up, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                            label="TRK SF Up",
                                            subplot=dijet_cen_nominal_sf_hist)
                dijet_cen_sf_up_hist = qgp.get_projection_plot(QCD_SF_UP_TFILE[cen_obj_name], start_val, end_val)
                dijet_cen_nominal_vs_syst_sf_entries.append((dijet_cen_sf_up_hist, qcd_cen_kwargs_sf_up))

                # TK SF DOWN
                col_down = qgc.QCD_COLOURS[3]
                qcd_cen_kwargs_sf_down = dict(line_color=col_down, line_width=lw, fill_color=col_down,
                                          marker_color=col_down, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                          label="TRK SF Down",
                                          subplot=dijet_cen_nominal_sf_hist)
                dijet_cen_sf_down_hist = qgp.get_projection_plot(QCD_SF_DOWN_TFILE[cen_obj_name], start_val, end_val)
                dijet_cen_nominal_vs_syst_sf_entries.append((dijet_cen_sf_down_hist, qcd_cen_kwargs_sf_down))

                dijet_cen_1d_entries.append([
                    (dijet_cen_no_sf_hist, qcd_cen_kwargs_no_sf),
                    (dijet_cen_nominal_sf_hist, qcd_cen_kwargs_mc),
                    (dijet_cen_sf_up_hist, qcd_cen_kwargs_sf_up),
                    (dijet_cen_sf_down_hist, qcd_cen_kwargs_sf_down),
                ])

                ####################
                # DIJET FORWARD REGION
                ####################
                fwd_obj_name = "%s/%s" % (dj_fwd_dirname, v)

                # NO TK SF
                dijet_fwd_no_sf_hist = qgp.get_projection_plot(QCD_NO_SF_TFILE[fwd_obj_name], start_val, end_val)
                qcd_fwd_kwargs_no_sf = dict(line_color=qgc.JETHT_COLOUR, line_width=data_line_width, fill_color=qgc.JETHT_COLOUR,
                                            marker_color=qgc.JETHT_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                            label="No TRK SF")
                dijet_fwd_no_vs_nominal_sf_entries.append((dijet_fwd_no_sf_hist, qcd_fwd_kwargs_no_sf))

                # NOMINAL TK SF
                dijet_fwd_nominal_sf_hist = qgp.get_projection_plot(QCD_NOMINAL_SF_TFILE[fwd_obj_name], start_val, end_val)

                qcd_fwd_kwargs_mc = dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                                         marker_color=qgc.QCD_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                         label="Nominal TRK SF",
                                         subplot=dijet_fwd_no_sf_hist)
                dijet_fwd_no_vs_nominal_sf_entries.append((dijet_fwd_nominal_sf_hist, qcd_fwd_kwargs_mc))

                qcd_fwd_kwargs_mc2 = dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                                          marker_color=qgc.QCD_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                          label="Nominal TRK SF")
                dijet_fwd_nominal_vs_syst_sf_entries.append((dijet_fwd_nominal_sf_hist.Clone(), qcd_fwd_kwargs_mc2))

                # TK SF UP
                qcd_fwd_kwargs_sf_up = dict(line_color=col_up, line_width=lw, fill_color=col_up,
                                            marker_color=col_up, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                            label="TRK SF Up",
                                            subplot=dijet_fwd_nominal_sf_hist)
                dijet_fwd_sf_up_hist = qgp.get_projection_plot(QCD_SF_UP_TFILE[fwd_obj_name], start_val, end_val)
                dijet_fwd_nominal_vs_syst_sf_entries.append((dijet_fwd_sf_up_hist, qcd_fwd_kwargs_sf_up))

                # TK SF DOWN
                qcd_fwd_kwargs_sf_down = dict(line_color=col_down, line_width=lw, fill_color=col_down,
                                          marker_color=col_down, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                          label="TRK SF Down",
                                          subplot=dijet_fwd_nominal_sf_hist)
                dijet_fwd_sf_down_hist = qgp.get_projection_plot(QCD_SF_DOWN_TFILE[fwd_obj_name], start_val, end_val)
                dijet_fwd_nominal_vs_syst_sf_entries.append((dijet_fwd_sf_down_hist, qcd_fwd_kwargs_sf_down))

                dijet_fwd_1d_entries.append([
                    (dijet_fwd_no_sf_hist, qcd_fwd_kwargs_no_sf),
                    (dijet_fwd_nominal_sf_hist, qcd_fwd_kwargs_mc),
                    (dijet_fwd_sf_up_hist, qcd_fwd_kwargs_sf_up),
                    (dijet_fwd_sf_down_hist, qcd_fwd_kwargs_sf_down),
                ])

                #################
                # SETUP PLOTTING
                #################
                rebin = 2
                v_lower = v.lower()
                if "multiplicity" in v_lower:
                    rebin = 2
                elif "flavour" in v_lower or "thrust" in v_lower:
                    rebin = 1
                elif 'ptd' in v_lower:
                    rebin = 2

                xlim = None
                if "width" in v_lower or "ptd" in v_lower:
                    xlim = (0, 1)
                elif"thrust" in v_lower:
                    xlim = (0, 0.5)
                elif "multiplicity" in v_lower and "ak4" in NOMINAL_SF_DIR.lower():
                    if end_val <= 150:
                        xlim = (0, 50)
                    else:
                        xlim = (0, 80)

                ylim = None
                if "flavour" in v_lower:
                    ylim = (0, 1)
                # elif "lha" in v_lower:
                    # ylim = (0, 5)

                # Do nominal vs no SF
                # --------------------------------------------------------------
                plot_dir = os.path.join(NOMINAL_SF_DIR, "plots_trk_sf_no_vs_nominal%s" % (gr_append))

                subplot_title = "Nominal / No SF"
                subplot_limits = (0.5, 1.5)

                xlabel = ang.name + " (" + ang.lambda_str + ")"
                if gr_append is not "":
                    xlabel = "Groomed " + ang.name + " (" + ang.lambda_str + ")"

                # dj central only
                qgp.do_comparison_plot(dijet_cen_no_vs_nominal_sf_entries,
                                       "%s/ptBinned/%s_pt%dto%d_dijet_central.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                       rebin=rebin,
                                       title="%d < p_{T}^{jet} < %d GeV\n%s\n%s" % (start_val, end_val, jet_str, qgc.Dijet_CEN_LABEL),
                                       xtitle=xlabel,
                                       xlim=xlim,
                                       ylim=ylim,
                                       has_data=False,
                                       subplot_type='ratio',
                                       subplot_title=subplot_title,
                                       subplot_limits=subplot_limits)

                # dj forward only
                qgp.do_comparison_plot(dijet_fwd_no_vs_nominal_sf_entries,
                                       "%s/ptBinned/%s_pt%dto%d_dijet_forward.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                       rebin=rebin,
                                       title="%d < p_{T}^{jet} < %d GeV\n%s\n%s" % (start_val, end_val, jet_str, qgc.Dijet_FWD_LABEL),
                                       xtitle=xlabel,
                                       xlim=xlim,
                                       ylim=ylim,
                                       has_data=False,
                                       subplot_type='ratio',
                                       subplot_title=subplot_title,
                                       subplot_limits=subplot_limits)

                # zpj only
                qgp.do_comparison_plot(zpj_no_vs_nominal_sf_entries,
                                       "%s/ptBinned/%s_pt%dto%d_zpj.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                       rebin=rebin,
                                       title="%d < p_{T}^{jet} < %d GeV\n%s\n%s" % (start_val, end_val, jet_str, qgc.ZpJ_LABEL),
                                       xtitle=xlabel,
                                       xlim=xlim,
                                       ylim=ylim,
                                       has_data=False,
                                       subplot_type='ratio',
                                       subplot_title=subplot_title,
                                       subplot_limits=subplot_limits)

                # Do nominal vs no SF
                # --------------------------------------------------------------
                plot_dir = os.path.join(NOMINAL_SF_DIR, "plots_trk_sf_nominal_vs_systs%s" % (gr_append))

                subplot_title = "Syst / Nominal"
                subplot_limits = (0.5, 1.5)

                # dj central only
                qgp.do_comparison_plot(dijet_cen_nominal_vs_syst_sf_entries,
                                       "%s/ptBinned/%s_pt%dto%d_dijet_central.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                       rebin=rebin,
                                       title="%d < p_{T}^{jet} < %d GeV\n%s\n%s" % (start_val, end_val, jet_str, qgc.Dijet_CEN_LABEL),
                                       xtitle=xlabel,
                                       xlim=xlim,
                                       ylim=ylim,
                                       has_data=False,
                                       subplot_type='ratio',
                                       subplot_title=subplot_title,
                                       subplot_limits=subplot_limits)

                # dj forward only
                qgp.do_comparison_plot(dijet_fwd_nominal_vs_syst_sf_entries,
                                       "%s/ptBinned/%s_pt%dto%d_dijet_forward.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                       rebin=rebin,
                                       title="%d < p_{T}^{jet} < %d GeV\n%s\n%s" % (start_val, end_val, jet_str, qgc.Dijet_FWD_LABEL),
                                       xtitle=xlabel,
                                       xlim=xlim,
                                       ylim=ylim,
                                       has_data=False,
                                       subplot_type='ratio',
                                       subplot_title=subplot_title,
                                       subplot_limits=subplot_limits)

                # zpj only
                qgp.do_comparison_plot(zpj_nominal_vs_syst_sf_entries,
                                       "%s/ptBinned/%s_pt%dto%d_zpj.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                       rebin=rebin,
                                       title="%d < p_{T}^{jet} < %d GeV\n%s\n%s" % (start_val, end_val, jet_str, qgc.ZpJ_LABEL),
                                       xtitle=xlabel,
                                       xlim=xlim,
                                       ylim=ylim,
                                       has_data=False,
                                       subplot_type='ratio',
                                       subplot_title=subplot_title,
                                       subplot_limits=subplot_limits)


            # Do overall summary plots across all pt bins
            # ------------------------------------------------------------------
            ylim_mean = None
            if "width" in v_lower or "ptd" in v_lower:
                ylim_mean = (0, 0.4)
            elif"thrust" in v_lower:
                ylim_mean = (0, 0.5)
            elif "multiplicity" in v_lower and "ak4" in NOMINAL_SF_DIR.lower():
                ylim_mean = (0, 100)
                ylim_mean = (0, 80)
                if end_val < 150:
                    ylim_mean = (0, 50)
            ylim_mean=None
            ylim_rms=None

            # Setup variable names for MPL
            marker = ""
            if "_" in ang.name or "^" in ang.name:
                marker = "$"
            var_label = marker + ang.name + marker + " ($%s$)" % ang.lambda_str
            if gr_append != "":
                var_label = "Groomed " + marker + ang.name + marker + " ($%s$)" % ang.lambda_str

            qgp.do_mean_rms_summary_plot(dijet_cen_1d_entries,
                                         pt_bins,
                                         "%s/ptBinned/%s_box_dijet_cen_mpl.%s" % (plot_dir, v, OUTPUT_FMT),
                                         var_label=var_label,
                                         xlim=(50, 2000),
                                         ylim_mean=ylim_mean,
                                         ylim_rms=ylim_rms,
                                         region_title="%s jets in %s" % (jet_str, qgc.Dijet_CEN_LABEL),
                                         has_data=False,
                                         numerator_label="Var.",
                                         denominator_label="No SF")

            qgp.do_mean_rms_summary_plot(dijet_fwd_1d_entries,
                                         pt_bins,
                                         "%s/ptBinned/%s_box_dijet_fwd_mpl.%s" % (plot_dir, v, OUTPUT_FMT),
                                         var_label=var_label,
                                         xlim=(50, 2000),
                                         ylim_mean=ylim_mean,
                                         ylim_rms=ylim_rms,
                                         region_title="%s jets in %s" % (jet_str, qgc.Dijet_FWD_LABEL),
                                         has_data=False,
                                         numerator_label="Var.",
                                         denominator_label="No SF")

            # # zpj_1d_entries[i][j] is the jth sample in the ith pt bin
            qgp.do_mean_rms_summary_plot(zpj_1d_entries,
                                         pt_bins,
                                         "%s/ptBinned/%s_box_zpj_mpl.%s" % (plot_dir, v, OUTPUT_FMT),
                                         var_label=var_label,
                                         xlim=(50, 614),
                                         ylim_mean=ylim_mean,
                                         ylim_rms=ylim_rms,
                                         region_title="%s jets in %s" % (jet_str, qgc.ZpJ_LABEL),
                                         has_data=False,
                                         numerator_label="Var.",
                                         denominator_label="No SF")


