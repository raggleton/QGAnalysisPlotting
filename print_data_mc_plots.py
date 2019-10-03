#!/usr/bin/env python

"""Print main QG plots for a given sample.

Any other plots comparing more than one sample should go into its own script!

USe addZeroBiasJetHT.py first to sum jetht and zerobias
"""

from __future__ import print_function

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
os.nice(10)

# My stuff
from comparator import Contribution, Plot, grab_obj
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

TOTAL_LUMI = 35918


def do_plots(root_dir):
    # QG variable plots
    pt_bins = qgc.PT_BINS[:]
    var_list = qgc.COMMON_VARS
    var_prepend = ""

    radius, pus = cu.get_jet_config_from_dirname(root_dir)
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
            zpj_2d_entries = []
            dijet_cen_2d_entries = []
            dijet_fwd_2d_entries = []

            zpj_1d_entries = []
            dijet_cen_1d_entries = []
            dijet_fwd_1d_entries = []

            for pt_ind, (start_val, end_val) in enumerate(pt_bins):
                entries = []
                zpj_entries = []
                dijet_cen_entries = []
                dijet_fwd_entries = []

                # Get all plots
                lw = 2
                msize = 1.1
                data_line_width = lw

                ####################
                # Z+JETS REGION
                ####################

                # SINGLE MU DATA
                h2d_dyj_data = grab_obj(os.path.join(root_dir, qgc.SINGLE_MU_FILENAME), "%s/%s" % (zpj_dirname, v))
                dy_kwargs_data = dict(line_color=qgc.SINGLE_MU_COLOUR, line_width=data_line_width, fill_color=qgc.SINGLE_MU_COLOUR,
                                      marker_color=qgc.SINGLE_MU_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=msize,
                                      label="Data")
                zpj_data_hist = qgp.get_projection_plot(h2d_dyj_data, start_val, end_val)
                entries.append((zpj_data_hist, dy_kwargs_data))
                zpj_entries.append((zpj_data_hist, dy_kwargs_data))
                if pt_ind == 0:
                    zpj_2d_entries.append((h2d_dyj_data, dy_kwargs_data))

                # PYTHIA DY MC
                h2d_dyj_mc = grab_obj(os.path.join(root_dir, qgc.DY_FILENAME), "%s/%s" % (zpj_dirname, v))
                dy_kwargs_mc = dict(line_color=qgc.DY_COLOUR, line_width=lw, fill_color=qgc.DY_COLOUR,
                                    marker_color=qgc.DY_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=0,
                                    label="MG+PY8",
                                    subplot=zpj_data_hist)
                entries.append((qgp.get_projection_plot(h2d_dyj_mc, start_val, end_val), dy_kwargs_mc))
                zpj_entries.append((qgp.get_projection_plot(h2d_dyj_mc, start_val, end_val), dy_kwargs_mc))
                if pt_ind == 0:
                    zpj_2d_entries.append((h2d_dyj_mc, dy_kwargs_mc))

                # HERWIG++ DY
                if end_val < 255:
                    h2d_dyj_mc2 = grab_obj(os.path.join(root_dir, qgc.DY_HERWIG_FILENAME), "%s/%s" % (zpj_dirname, v))
                    col3 = qgc.DY_COLOURS[2]
                    dy_kwargs_mc2 = dict(line_color=col3, line_width=lw, fill_color=col3,
                                         marker_color=col3, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=0,
                                         label="H++",
                                         subplot=zpj_data_hist)
                    entries.append((qgp.get_projection_plot(h2d_dyj_mc2, start_val, end_val), dy_kwargs_mc2))
                    zpj_entries.append((qgp.get_projection_plot(h2d_dyj_mc2, start_val, end_val), dy_kwargs_mc2))
                    if pt_ind == 0:
                        zpj_2d_entries.append((h2d_dyj_mc2, dy_kwargs_mc2))
                else:
                    zpj_entries.append(None)

                # MG+HERWIG++ DY
                if end_val < 151:
                    h2d_dyj_mc3 = grab_obj(os.path.join(root_dir, qgc.DY_MG_HERWIG_FILENAME), "%s/%s" % (zpj_dirname, v))
                    col4 = qgc.DY_COLOURS[3]
                    dy_kwargs_mc3 = dict(line_color=col4, line_width=lw, fill_color=col4,
                                         marker_color=col4, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=0,
                                         label="MG+H++",
                                         subplot=zpj_data_hist)
                    entries.append((qgp.get_projection_plot(h2d_dyj_mc3, start_val, end_val), dy_kwargs_mc3))
                    zpj_entries.append((qgp.get_projection_plot(h2d_dyj_mc3, start_val, end_val), dy_kwargs_mc3))
                    if pt_ind == 0:
                        zpj_2d_entries.append((h2d_dyj_mc3, dy_kwargs_mc3))
                else:
                    zpj_entries.append(None)

                ####################
                # DIJET CENTRAL REGION
                ####################

                # JETHT/ZEROBIAS DATA
                h2d_qcd_cen_data = grab_obj(os.path.join(root_dir, qgc.JETHT_ZB_FILENAME), "%s/%s" % (dj_cen_dirname, v))  # use already merged jetht+zb
                # h2d_zb_data = grab_obj(os.path.join(root_dir, qgc.ZB_FILENAME), "%s/%s" % (dj_cen_dirname, v))
                # h2d_zb_data.Scale(1235009.27580634)
                # h2d_qcd_cen_data.Add(h2d_zb_data)
                qcd_cen_kwargs_data = dict(line_color=qgc.JETHT_COLOUR, line_width=data_line_width, fill_color=qgc.JETHT_COLOUR,
                                           marker_color=qgc.JETHT_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=msize,
                                           label="Data")
                dijet_cen_data_hist = qgp.get_projection_plot(h2d_qcd_cen_data, start_val, end_val)
                entries.append((dijet_cen_data_hist, qcd_cen_kwargs_data))
                dijet_cen_entries.append((dijet_cen_data_hist, qcd_cen_kwargs_data))
                if pt_ind == 0:
                    dijet_cen_2d_entries.append((h2d_qcd_cen_data, qcd_cen_kwargs_data))

                # MG+PYTHIA QCD MC
                h2d_qcd_cen_mc = grab_obj(os.path.join(root_dir, qgc.QCD_FILENAME), "%s/%s" % (dj_cen_dirname, v))
                qcd_cen_kwargs_mc = dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                                         marker_color=qgc.QCD_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                         label="MG+PY8",
                                         subplot=dijet_cen_data_hist)
                # h2d_qcd_cen_mc.Scale(35860)
                dijet_mgpy_hist = qgp.get_projection_plot(h2d_qcd_cen_mc, start_val, end_val)
                entries.append((dijet_mgpy_hist, qcd_cen_kwargs_mc))
                dijet_cen_entries.append((dijet_mgpy_hist, qcd_cen_kwargs_mc))
                if pt_ind == 0:
                    dijet_cen_2d_entries.append((h2d_qcd_cen_mc, qcd_cen_kwargs_mc))

                # PYTHIA ONLY
                col = qgc.QCD_COLOURS[2]
                h2d_qcd_cen_mc2 = grab_obj(os.path.join(root_dir, qgc.QCD_PYTHIA_ONLY_FILENAME), "%s/%s" % (dj_cen_dirname, v))
                qcd_cen_kwargs_mc2 = dict(line_color=col, line_width=lw, fill_color=col,
                                          marker_color=col, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                          label="PY8",
                                          subplot=dijet_cen_data_hist)
                # h2d_qcd_cen_mc2.Scale(35860)
                entries.append((qgp.get_projection_plot(h2d_qcd_cen_mc2, start_val, end_val), qcd_cen_kwargs_mc2))
                dijet_cen_entries.append((qgp.get_projection_plot(h2d_qcd_cen_mc2, start_val, end_val), qcd_cen_kwargs_mc2))
                if pt_ind == 0:
                    dijet_cen_2d_entries.append((h2d_qcd_cen_mc2, qcd_cen_kwargs_mc2))

                # HERWIG++ QCD
                col2 = qgc.QCD_COLOURS[3]
                h2d_qcd_cen_mc3 = grab_obj(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME), "%s/%s" % (dj_cen_dirname, v))
                qcd_cen_kwargs_mc3 = dict(line_color=col2, line_width=lw, fill_color=col2,
                                          marker_color=col2, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                          label="H++",
                                          subplot=dijet_cen_data_hist)
                h2d_qcd_cen_mc3.Scale(TOTAL_LUMI)
                dijet_hpp_hist = qgp.get_projection_plot(h2d_qcd_cen_mc3, start_val, end_val)
                entries.append((dijet_hpp_hist, qcd_cen_kwargs_mc3))
                dijet_cen_entries.append((dijet_hpp_hist, qcd_cen_kwargs_mc3))
                if pt_ind == 0:
                    dijet_cen_2d_entries.append((h2d_qcd_cen_mc3, qcd_cen_kwargs_mc3))

                ####################
                # DIJET FORWARD REGION
                ####################

                # JETHT/ZEROBIAS DATA
                h2d_qcd_fwd_data = grab_obj(os.path.join(root_dir, qgc.JETHT_ZB_FILENAME), "%s/%s" % (dj_cen_dirname, v))  # use already merged jetht+zb
                # h2d_zb_data = grab_obj(os.path.join(root_dir, qgc.ZB_FILENAME), "%s/%s" % (dj_cen_dirname, v))
                # h2d_zb_data.Scale(1235009.27580634)
                # h2d_qcd_cen_data.Add(h2d_zb_data)
                qcd_fwd_kwargs_data = dict(line_color=qgc.JETHT_COLOUR, line_width=data_line_width, fill_color=qgc.JETHT_COLOUR,
                                           marker_color=qgc.JETHT_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=msize,
                                           label="Data")
                dijet_fwd_data_hist = qgp.get_projection_plot(h2d_qcd_fwd_data, start_val, end_val)
                entries.append((dijet_fwd_data_hist, qcd_fwd_kwargs_data))
                dijet_fwd_entries.append((dijet_fwd_data_hist, qcd_fwd_kwargs_data))
                if pt_ind == 0:
                    dijet_fwd_2d_entries.append((h2d_qcd_fwd_data, qcd_fwd_kwargs_data))

                # MG+PYTHIA QCD MC
                h2d_qcd_fwd_mc = grab_obj(os.path.join(root_dir, qgc.QCD_FILENAME), "%s/%s" % (dj_fwd_dirname, v))
                qcd_fwd_kwargs_mc = dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                                         marker_color=qgc.QCD_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                         label="MG+PY8",
                                         subplot=dijet_fwd_data_hist)
                # h2d_qcd_fwd_mc.Scale(35860)
                dijet_mgpy_hist = qgp.get_projection_plot(h2d_qcd_fwd_mc, start_val, end_val)
                entries.append((dijet_mgpy_hist, qcd_fwd_kwargs_mc))
                dijet_fwd_entries.append((dijet_mgpy_hist, qcd_fwd_kwargs_mc))
                if pt_ind == 0:
                    dijet_fwd_2d_entries.append((h2d_qcd_fwd_mc, qcd_fwd_kwargs_mc))

                # PYTHIA ONLY
                col = qgc.QCD_COLOURS[2]
                h2d_qcd_fwd_mc2 = grab_obj(os.path.join(root_dir, qgc.QCD_PYTHIA_ONLY_FILENAME), "%s/%s" % (dj_fwd_dirname, v))
                qcd_fwd_kwargs_mc2 = dict(line_color=col, line_width=lw, fill_color=col,
                                          marker_color=col, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                          label="PY8",
                                          subplot=dijet_fwd_data_hist)
                # h2d_qcd_fwd_mc2.Scale(35860)
                entries.append((qgp.get_projection_plot(h2d_qcd_fwd_mc2, start_val, end_val), qcd_fwd_kwargs_mc2))
                dijet_fwd_entries.append((qgp.get_projection_plot(h2d_qcd_fwd_mc2, start_val, end_val), qcd_fwd_kwargs_mc2))
                if pt_ind == 0:
                    dijet_fwd_2d_entries.append((h2d_qcd_fwd_mc2, qcd_fwd_kwargs_mc2))

                # HERWIG++ QCD
                col2 = qgc.QCD_COLOURS[3]
                h2d_qcd_fwd_mc3 = grab_obj(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME), "%s/%s" % (dj_fwd_dirname, v))
                qcd_fwd_kwargs_mc3 = dict(line_color=col2, line_width=lw, fill_color=col2,
                                          marker_color=col2, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                          label="H++",
                                          subplot=dijet_fwd_data_hist)
                h2d_qcd_fwd_mc3.Scale(TOTAL_LUMI)
                dijet_hpp_hist = qgp.get_projection_plot(h2d_qcd_fwd_mc3, start_val, end_val)
                entries.append((dijet_hpp_hist, qcd_fwd_kwargs_mc3))
                dijet_fwd_entries.append((dijet_hpp_hist, qcd_fwd_kwargs_mc3))
                if pt_ind == 0:
                    dijet_fwd_2d_entries.append((h2d_qcd_fwd_mc3, qcd_fwd_kwargs_mc3))

                zpj_1d_entries.append(zpj_entries)
                dijet_cen_1d_entries.append(dijet_cen_entries)
                dijet_fwd_1d_entries.append(dijet_fwd_entries)

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
                elif "multiplicity" in v_lower and "ak4" in root_dir.lower():
                    if end_val <= 150:
                        xlim = (0, 50)
                    else:
                        xlim = (0, 80)

                ylim = None
                if "flavour" in v_lower:
                    ylim = (0, 1)
                # elif "lha" in v_lower:
                    # ylim = (0, 5)

                plot_dir = os.path.join(root_dir, "plots_dy_vs_qcd_mc_vs_data%s" % (gr_append))

                subplot_title = "MC / Data"
                subplot_limits = (0.5, 1.5)

                xlabel = ang.name + " (" + ang.lambda_str + ")"
                if gr_append is not "":
                    xlabel = "Groomed " + ang.name + " (" + ang.lambda_str + ")"

                # dj+zpj
                # qgp.do_comparison_plot(entries,
                #                        "%s/ptBinned/%s_pt%dto%d.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                #                        rebin=rebin,
                #                        title="%d < p_{T}^{jet} < %d GeV\n%s" % (start_val, end_val, jet_str),
                #                        xtitle=ang.name + " (" + ang.lambda_str + ")",
                #                        xlim=xlim, ylim=ylim,
                #                        subplot_type='ratio',
                #                        subplot_title=subplot_title,
                #                        subplot_limits=subplot_limits)

                # dj central only
                qgp.do_comparison_plot(dijet_cen_entries,
                                       "%s/ptBinned/%s_pt%dto%d_dijet_central.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                       rebin=rebin,
                                       title="%d < p_{T}^{jet} < %d GeV\n%s\n%s" % (start_val, end_val, jet_str, qgc.Dijet_CEN_LABEL),
                                       xtitle=xlabel,
                                       xlim=xlim,
                                       ylim=ylim,
                                       subplot_type='ratio',
                                       subplot_title=subplot_title,
                                       subplot_limits=subplot_limits)

                # dj forward only
                qgp.do_comparison_plot(dijet_fwd_entries,
                                       "%s/ptBinned/%s_pt%dto%d_dijet_forward.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                       rebin=rebin,
                                       title="%d < p_{T}^{jet} < %d GeV\n%s\n%s" % (start_val, end_val, jet_str, qgc.Dijet_FWD_LABEL),
                                       xtitle=xlabel,
                                       xlim=xlim,
                                       ylim=ylim,
                                       subplot_type='ratio',
                                       subplot_title=subplot_title,
                                       subplot_limits=subplot_limits)

                # zpj only
                qgp.do_comparison_plot(zpj_entries,
                                       "%s/ptBinned/%s_pt%dto%d_zpj.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                       rebin=rebin,
                                       title="%d < p_{T}^{jet} < %d GeV\n%s\n%s" % (start_val, end_val, jet_str, qgc.ZpJ_LABEL),
                                       xtitle=xlabel,
                                       xlim=xlim,
                                       ylim=ylim,
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
            elif "multiplicity" in v_lower and "ak4" in root_dir.lower():
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

            qgp.do_mean_rms_summary_plot(dijet_cen_1d_entries[:],
                                         pt_bins[:],
                                         "%s/ptBinned/%s_box_dijet_cen_mpl.%s" % (plot_dir, v, OUTPUT_FMT),
                                         var_label=var_label,
                                         xlim=(50, 2000),
                                         ylim_mean=ylim_mean,
                                         ylim_rms=ylim_rms,
                                         region_title="%s jets in dijet (central)" % (jet_str))

            qgp.do_mean_rms_summary_plot(dijet_fwd_1d_entries[:],
                                         pt_bins[:],
                                         "%s/ptBinned/%s_box_dijet_fwd_mpl.%s" % (plot_dir, v, OUTPUT_FMT),
                                         var_label=var_label,
                                         xlim=(50, 2000),
                                         ylim_mean=ylim_mean,
                                         ylim_rms=ylim_rms,
                                         region_title="%s jets in dijet (forward)" % (jet_str))

            # zpj_1d_entries[i][j] is the jth sample in the ith pt bin
            qgp.do_mean_rms_summary_plot(zpj_1d_entries[:],
                                         pt_bins[:],
                                         "%s/ptBinned/%s_box_zpj_mpl.%s" % (plot_dir, v, OUTPUT_FMT),
                                         var_label=var_label,
                                         xlim=(50, 614),
                                         ylim_mean=ylim_mean,
                                         ylim_rms=ylim_rms,
                                         region_title="%s jets in Z+jets" % (jet_str))

            # qgp.do_box_plot(dijet_2d_entries,
            #                 "%s/ptBinned/%s_box_dijet.%s" % (plot_dir, v, OUTPUT_FMT),
            #                 ylim=ylim, xlim=(0, 300), transpose=True)

            # qgp.do_box_plot(zpj_2d_entries,
            #                 "%s/ptBinned/%s_box_zpj.%s" % (plot_dir, v, OUTPUT_FMT),
            #                 ylim=ylim, xlim=(0, 300), transpose=True)


if __name__ == "__main__":
    parser = qgc.get_parser()
    args = parser.parse_args()

    for workdir in args.workdirs:
        do_plots(workdir)
