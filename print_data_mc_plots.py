#!/usr/bin/env python

"""Print main QG lambda plots for a given sample.

Any other plots comparing more than one sample should go into its own script!

Make sure jetht and zerobias are hadded into same file
"""

from __future__ import print_function

import os
os.nice(10)
import ROOT
from MyStyle import My_Style
My_Style.cd()
import argparse

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

TOTAL_LUMI = 35918


def do_plots(root_dir, title):
    # QG variable plots
    pt_bins = qgc.PT_BINS[:]
    print(pt_bins)
    var_list = qgc.COMMON_VARS
    var_prepend = ""

    radius, pus = cu.get_jet_config_from_dirname(root_dir)
    jet_str = "AK%s PF %s" % (radius.upper(), pus.upper())

    if radius == "8":
        pt_bins = qgc.PT_BINS[2:]
        print(pt_bins)

    single_mu_tfile = cu.TFileCacher(os.path.join(root_dir, qgc.SINGLE_MU_FILENAME))
    dy_tfile = cu.TFileCacher(os.path.join(root_dir, qgc.DY_FILENAME))
    dy_hpp_tfile = cu.TFileCacher(os.path.join(root_dir, qgc.DY_HERWIG_FILENAME))
    jetht_zb_tfile = cu.TFileCacher(os.path.join(root_dir, qgc.JETHT_ZB_FILENAME))
    qcd_tfile = cu.TFileCacher(os.path.join(root_dir, qgc.QCD_FILENAME))
    qcd_hpp_tfile = cu.TFileCacher(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME))

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
                mc_msize = 1E-3  # small enough to not be seen but not 0 - otherwise the horizontal lines on error bars don't get drawn

                mgpy_label = "MG5+Pythia8"
                hpp_label = "Herwig++"

                ####################
                # Z+JETS REGION
                ####################

                # SINGLE MU DATA
                if single_mu_tfile:
                    h2d_dyj_data = single_mu_tfile.get("%s/%s" % (zpj_dirname, v))
                    dy_kwargs_data = dict(line_color=qgc.SINGLE_MU_COLOUR, line_width=data_line_width, fill_color=qgc.SINGLE_MU_COLOUR,
                                          marker_color=qgc.SINGLE_MU_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=msize*0.7,
                                          label="Data")
                    zpj_data_hist = qgp.get_projection_plot(h2d_dyj_data, start_val, end_val)
                    entries.append((zpj_data_hist, dy_kwargs_data))
                    zpj_entries.append((zpj_data_hist, dy_kwargs_data))
                    if pt_ind == 0:
                        zpj_2d_entries.append((h2d_dyj_data, dy_kwargs_data))

                # PYTHIA DY MC
                if dy_tfile:
                    h2d_dyj_mc = dy_tfile.get("%s/%s" % (zpj_dirname, v))
                    dy_kwargs_mc = dict(line_color=qgc.DY_COLOUR, line_width=lw, fill_color=qgc.DY_COLOUR,
                                        marker_color=qgc.DY_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=mc_msize,
                                        label=mgpy_label,
                                        subplot=zpj_data_hist)
                    entries.append((qgp.get_projection_plot(h2d_dyj_mc, start_val, end_val), dy_kwargs_mc))
                    zpj_entries.append((qgp.get_projection_plot(h2d_dyj_mc, start_val, end_val), dy_kwargs_mc))
                    if pt_ind == 0:
                        zpj_2d_entries.append((h2d_dyj_mc, dy_kwargs_mc))

                # HERWIG++ DY
                if dy_hpp_tfile:
                    h2d_dyj_mc_hpp = dy_hpp_tfile.get("%s/%s" % (zpj_dirname, v))
                    col_hpp = qgc.DY_COLOURS[2]
                    dy_kwargs_mc_hpp = dict(line_color=col_hpp, line_width=lw, fill_color=col_hpp,
                                            marker_color=col_hpp, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=mc_msize,
                                            label=hpp_label,
                                            subplot=zpj_data_hist)
                    entries.append((qgp.get_projection_plot(h2d_dyj_mc_hpp, start_val, end_val), dy_kwargs_mc_hpp))
                    zpj_entries.append((qgp.get_projection_plot(h2d_dyj_mc_hpp, start_val, end_val), dy_kwargs_mc_hpp))
                    if pt_ind == 0:
                        zpj_2d_entries.append((h2d_dyj_mc_hpp, dy_kwargs_mc_hpp))

                # MG+HERWIG++ DY
                # if end_val < 151:
                #     h2d_dyj_mc3 = get("%s/%s" % (zpj_dirname, v))
                #     col4 = qgc.DY_COLOURS[3]
                #     dy_kwargs_mc3 = dict(line_color=col4, line_width=lw, fill_color=col4,
                #                          marker_color=col4, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=0,
                #                          label="MG+Herwig++",
                #                          subplot=zpj_data_hist)
                #     entries.append((qgp.get_projection_plot(h2d_dyj_mc3, start_val, end_val), dy_kwargs_mc3))
                #     zpj_entries.append((qgp.get_projection_plot(h2d_dyj_mc3, start_val, end_val), dy_kwargs_mc3))
                #     if pt_ind == 0:
                #         zpj_2d_entries.append((h2d_dyj_mc3, dy_kwargs_mc3))
                # else:
                #     zpj_entries.append(None)

                ####################
                # DIJET CENTRAL REGION
                ####################

                # JETHT/ZEROBIAS DATA
                if jetht_zb_tfile:
                    h2d_qcd_cen_data = jetht_zb_tfile.get("%s/%s" % (dj_cen_dirname, v))  # use already merged jetht+zb
                    qcd_cen_kwargs_data = dict(line_color=qgc.JETHT_COLOUR, line_width=data_line_width, fill_color=qgc.JETHT_COLOUR,
                                               marker_color=qgc.JETHT_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=msize,
                                               label="Data")
                    dijet_cen_data_hist = qgp.get_projection_plot(h2d_qcd_cen_data, start_val, end_val)
                    entries.append((dijet_cen_data_hist, qcd_cen_kwargs_data))
                    dijet_cen_entries.append((dijet_cen_data_hist, qcd_cen_kwargs_data))
                    if pt_ind == 0:
                        dijet_cen_2d_entries.append((h2d_qcd_cen_data, qcd_cen_kwargs_data))

                # MG+PYTHIA QCD MC
                if qcd_tfile:
                    h2d_qcd_cen_mc = qcd_tfile.get("%s/%s" % (dj_cen_dirname, v))
                    qcd_cen_kwargs_mc = dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                                             marker_color=qgc.QCD_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=mc_msize,
                                             label=mgpy_label,
                                             subplot=dijet_cen_data_hist)
                    dijet_mgpy_hist = qgp.get_projection_plot(h2d_qcd_cen_mc, start_val, end_val)
                    entries.append((dijet_mgpy_hist, qcd_cen_kwargs_mc))
                    dijet_cen_entries.append((dijet_mgpy_hist, qcd_cen_kwargs_mc))
                    if pt_ind == 0:
                        dijet_cen_2d_entries.append((h2d_qcd_cen_mc, qcd_cen_kwargs_mc))

                # PYTHIA ONLY
                # if qcd_py_tfile:
                #     col = qgc.QCD_COLOURS[2]
                #     h2d_qcd_cen_mc2 = qcd_py_tfile.get("%s/%s" % (dj_cen_dirname, v))
                #     qcd_cen_kwargs_mc2 = dict(line_color=col, line_width=lw, fill_color=col,
                #                               marker_color=col, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=mc_msize,
                #                               label="Pythia8",
                #                               subplot=dijet_cen_data_hist)
                #     entries.append((qgp.get_projection_plot(h2d_qcd_cen_mc2, start_val, end_val), qcd_cen_kwargs_mc2))
                #     dijet_cen_entries.append((qgp.get_projection_plot(h2d_qcd_cen_mc2, start_val, end_val), qcd_cen_kwargs_mc2))
                #     if pt_ind == 0:
                #         dijet_cen_2d_entries.append((h2d_qcd_cen_mc2, qcd_cen_kwargs_mc2))

                # HERWIG++ QCD
                if qcd_hpp_tfile:
                    h2d_qcd_cen_mc_hpp = qcd_hpp_tfile.get("%s/%s" % (dj_cen_dirname, v))
                    qcd_cen_kwargs_mc_hpp = dict(line_color=qgc.HERWIGPP_QCD_COLOUR, line_width=lw, fill_color=qgc.HERWIGPP_QCD_COLOUR,
                                                 marker_color=qgc.HERWIGPP_QCD_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=mc_msize,
                                                 label=hpp_label,
                                                 subplot=dijet_cen_data_hist)
                    dijet_hpp_hist = qgp.get_projection_plot(h2d_qcd_cen_mc_hpp, start_val, end_val)
                    entries.append((dijet_hpp_hist, qcd_cen_kwargs_mc_hpp))
                    dijet_cen_entries.append((dijet_hpp_hist, qcd_cen_kwargs_mc_hpp))
                    if pt_ind == 0:
                        dijet_cen_2d_entries.append((h2d_qcd_cen_mc_hpp, qcd_cen_kwargs_mc_hpp))

                ####################
                # DIJET FORWARD REGION
                ####################

                # JETHT/ZEROBIAS DATA
                if jetht_zb_tfile:
                    h2d_qcd_fwd_data = jetht_zb_tfile.get("%s/%s" % (dj_fwd_dirname, v))  # use already merged jetht+zb
                    qcd_fwd_kwargs_data = dict(line_color=qgc.JETHT_COLOUR, line_width=data_line_width, fill_color=qgc.JETHT_COLOUR,
                                               marker_color=qgc.JETHT_COLOUR, marker_style=cu.Marker.get('triangleDown'), marker_size=msize,
                                               label="Data")
                    dijet_fwd_data_hist = qgp.get_projection_plot(h2d_qcd_fwd_data, start_val, end_val)
                    entries.append((dijet_fwd_data_hist, qcd_fwd_kwargs_data))
                    dijet_fwd_entries.append((dijet_fwd_data_hist, qcd_fwd_kwargs_data))
                    if pt_ind == 0:
                        dijet_fwd_2d_entries.append((h2d_qcd_fwd_data, qcd_fwd_kwargs_data))

                # MG+PYTHIA QCD MC
                if qcd_tfile:
                    h2d_qcd_fwd_mc = qcd_tfile.get("%s/%s" % (dj_fwd_dirname, v))
                    qcd_fwd_kwargs_mc = dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                                             marker_color=qgc.QCD_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=mc_msize,
                                             label=mgpy_label,
                                             subplot=dijet_fwd_data_hist)
                    dijet_mgpy_hist = qgp.get_projection_plot(h2d_qcd_fwd_mc, start_val, end_val)
                    entries.append((dijet_mgpy_hist, qcd_fwd_kwargs_mc))
                    dijet_fwd_entries.append((dijet_mgpy_hist, qcd_fwd_kwargs_mc))
                    if pt_ind == 0:
                        dijet_fwd_2d_entries.append((h2d_qcd_fwd_mc, qcd_fwd_kwargs_mc))

                # PYTHIA ONLY
                # if qcd_py_tfile:
                    # col = qgc.QCD_COLOURS[2]
                    # h2d_qcd_fwd_mc2 = qcd_py_tfile.get("%s/%s" % (dj_fwd_dirname, v))
                    # qcd_fwd_kwargs_mc2 = dict(line_color=col, line_width=lw, fill_color=col,
                    #                           marker_color=col, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=mc_msize,
                    #                           label="Pythia8",
                    #                           subplot=dijet_fwd_data_hist)
                    # entries.append((qgp.get_projection_plot(h2d_qcd_fwd_mc2, start_val, end_val), qcd_fwd_kwargs_mc2))
                    # dijet_fwd_entries.append((qgp.get_projection_plot(h2d_qcd_fwd_mc2, start_val, end_val), qcd_fwd_kwargs_mc2))
                    # if pt_ind == 0:
                    #     dijet_fwd_2d_entries.append((h2d_qcd_fwd_mc2, qcd_fwd_kwargs_mc2))

                # HERWIG++ QCD
                if qcd_hpp_tfile:
                    h2d_qcd_fwd_mc_hpp = qcd_hpp_tfile.get("%s/%s" % (dj_fwd_dirname, v))
                    qcd_fwd_kwargs_mc_hpp = dict(line_color=qgc.HERWIGPP_QCD_COLOUR, line_width=lw, fill_color=qgc.HERWIGPP_QCD_COLOUR,
                                                 marker_color=qgc.HERWIGPP_QCD_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=mc_msize,
                                                 label=hpp_label,
                                                 subplot=dijet_fwd_data_hist)
                    dijet_hpp_hist = qgp.get_projection_plot(h2d_qcd_fwd_mc_hpp, start_val, end_val)
                    entries.append((dijet_hpp_hist, qcd_fwd_kwargs_mc_hpp))
                    dijet_fwd_entries.append((dijet_hpp_hist, qcd_fwd_kwargs_mc_hpp))
                    if pt_ind == 0:
                        dijet_fwd_2d_entries.append((h2d_qcd_fwd_mc_hpp, qcd_fwd_kwargs_mc_hpp))

                zpj_1d_entries.append(zpj_entries)
                dijet_cen_1d_entries.append(dijet_cen_entries)
                dijet_fwd_1d_entries.append(dijet_fwd_entries)

                #################
                # SETUP PLOTTING
                #################
                rebin = 2
                v_lower = v.lower()
                if "multiplicity" in v_lower:
                    rebin = 1
                # elif "flavour" in v_lower or "thrust" in v_lower:
                #     rebin = 1
                # elif 'ptd' in v_lower:
                #     rebin = 4
                # elif 'lha_charged' in v_lower:
                #     rebin = 4
                # rebin = 1

                xlim = None
                if "width" in v_lower or "ptd" in v_lower:
                    xlim = (0, 1)
                elif "thrust" in v_lower:
                    xlim = (0, 0.5)
                # elif "multiplicity" in v_lower and "ak4" in root_dir.lower():
                #     if end_val <= 150:
                #         xlim = (0, 50)
                #     else:
                #         xlim = (0, 80)

                auto_xlim = False
                if "multiplicity" in v_lower or "thrust" in v_lower:
                    auto_xlim = True

                ylim = [0, None]
                if "flavour" in v_lower:
                    ylim = (0, 1)
                # elif "lha" in v_lower:
                    # ylim = (0, 5)

                plot_dir = os.path.join(root_dir, "plots_lambda_mc_vs_data%s" % (gr_append))

                subplot_title = "Simulation / Data"
                subplot_limits = (0, 2)

                xlabel = ang.name + " (" + ang.lambda_str + ")"
                if gr_append is not "":
                    xlabel = "Groomed " + ang.name + " (" + ang.lambda_str + ")"

                def _title(region_str, start_val, end_val):
                    pt_var_str = "p_{T}^{jet}"
                    s = (("{jet_algo}\n"
                          "{region_label}\n"
                          "{bin_edge_low:g} < {pt_str} < {bin_edge_high:g} GeV")
                          .format(
                            jet_algo=jet_str,
                            region_label=region_str,
                            pt_str=pt_var_str,
                            bin_edge_low=start_val,
                            bin_edge_high=end_val))
                    if title is not None:
                        s = "%s\n%s" % (s, title)
                    return s

                has_data = True
                draw_opt = "NOSTACK HIST E1"
                # dj central only
                # dont' try to do multiple signal regions per plot, it looks rubbish
                if len(dijet_cen_entries) > 0:
                    qgp.do_comparison_plot(dijet_cen_entries,
                                           "%s/ptBinned/%s_pt%dto%d_dijet_central.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                           rebin=rebin,
                                           draw_opt=draw_opt,
                                           title=_title(qgc.Dijet_CEN_LABEL, start_val, end_val),
                                           xtitle=xlabel,
                                           xlim='auto' if auto_xlim else xlim, # don't use calc_auto_xlim, since do_comparison_plot will rebin it anyway
                                           ylim=ylim,
                                           data_first=has_data,
                                           mean_rel_error=0.4,
                                           subplot_type='ratio',
                                           subplot_title=subplot_title,
                                           subplot_limits=subplot_limits,
                                           lumi=cu.get_lumi_str(do_dijet=True, do_zpj=False))

                # dj forward only
                if len(dijet_fwd_entries) > 0:
                    qgp.do_comparison_plot(dijet_fwd_entries,
                                           "%s/ptBinned/%s_pt%dto%d_dijet_forward.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                           rebin=rebin,
                                           draw_opt=draw_opt,
                                           title=_title(qgc.Dijet_FWD_LABEL, start_val, end_val),
                                           xtitle=xlabel,
                                           xlim='auto' if auto_xlim else xlim,
                                           ylim=ylim,
                                           data_first=has_data,
                                           mean_rel_error=0.4,
                                           subplot_type='ratio',
                                           subplot_title=subplot_title,
                                           subplot_limits=subplot_limits,
                                           lumi=cu.get_lumi_str(do_dijet=True, do_zpj=False))

                # zpj only
                if len(zpj_entries) > 0:
                    if start_val > 149:
                        # rebin *= 2
                        rebin += 1
                        # find nearest divisor
                        while (zpj_entries[0][0].GetNbinsX() % rebin != 0):
                            rebin += 1
                    if "multiplicity" in v_lower:
                        if start_val > 300:
                            rebin += 1
                            # find nearest divisor
                            while (zpj_entries[0][0].GetNbinsX() % rebin != 0):
                                rebin += 1
                    qgp.do_comparison_plot(zpj_entries,
                                           "%s/ptBinned/%s_pt%dto%d_zpj.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                           rebin=rebin,
                                           draw_opt=draw_opt,
                                           title=_title(qgc.ZpJ_LABEL, start_val, end_val),
                                           xtitle=xlabel,
                                           xlim=qgp.calc_auto_xlim([d[0] for d in zpj_entries]) if auto_xlim else xlim,
                                           ylim=ylim,
                                           data_first=has_data,
                                           mean_rel_error=0.4,
                                           subplot_type='ratio',
                                           subplot_title=subplot_title,
                                           subplot_limits=subplot_limits,
                                           lumi=cu.get_lumi_str(do_dijet=False, do_zpj=True))


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

            # construct a dataframe from the hists
            # df = construct_dataframe_from_hists(dijet_cen_entries=dijet_cen_1d_entries,
            #                                     dijet_fwd_entries=dijet_fwd_entries,
            #                                     zpj_entries=zpj_entries)

            # pt_low = 88 if "8" in radius else 50 
            pt_low = pt_bins[0][0]
            qgp.do_mean_rms_summary_plot(dijet_cen_1d_entries,
                                         pt_bins,
                                         "%s/ptBinned/%s_box_dijet_cen_mpl.%s" % (plot_dir, v, OUTPUT_FMT),
                                         var_label=var_label,
                                         xlim=(pt_low, 4000),
                                         ylim_mean=ylim_mean,
                                         ylim_rms=ylim_rms,
                                         region_title="%s jets in %s" % (jet_str, qgc.Dijet_CEN_LABEL.lower()))

            qgp.do_mean_rms_summary_plot(dijet_fwd_1d_entries,
                                         pt_bins,
                                         "%s/ptBinned/%s_box_dijet_fwd_mpl.%s" % (plot_dir, v, OUTPUT_FMT),
                                         var_label=var_label,
                                         xlim=(pt_low, 4000),
                                         ylim_mean=ylim_mean,
                                         ylim_rms=ylim_rms,
                                         region_title="%s jets in %s" % (jet_str, qgc.Dijet_FWD_LABEL.lower()))

            # zpj_1d_entries[i][j] is the jth sample in the ith pt bin
            qgp.do_mean_rms_summary_plot(zpj_1d_entries,
                                         pt_bins,
                                         "%s/ptBinned/%s_box_zpj_mpl.%s" % (plot_dir, v, OUTPUT_FMT),
                                         var_label=var_label,
                                         xlim=(pt_low, 800),
                                         ylim_mean=ylim_mean,
                                         ylim_rms=ylim_rms,
                                         region_title="%s jets in %s" % (jet_str, qgc.ZpJ_LABEL))

            # # qgp.do_box_plot(dijet_2d_entries,
            # #                 "%s/ptBinned/%s_box_dijet.%s" % (plot_dir, v, OUTPUT_FMT),
            # #                 ylim=ylim, xlim=(0, 300), transpose=True)

            # # qgp.do_box_plot(zpj_2d_entries,
            # #                 "%s/ptBinned/%s_box_zpj.%s" % (plot_dir, v, OUTPUT_FMT),
            # #                 ylim=ylim, xlim=(0, 300), transpose=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workdir',
                        help='Workdir with ROOT files to process. '
                        'Each directory must have all the necessary ROOT files')
    # parser.add_argument("-o", "--output",
    #                     help="Directory to put output plot dirs into. Default is workdir.",
    #                     default=None)
    parser.add_argument("--doExperimentalSysts",
                        help="Also draw experimental systematic variations",
                        action='store_true')
    parser.add_argument("--doModelSysts",
                        help="Also draw model systematic variations",
                        action='store_true')

    parser.add_argument("--title",
                        help="Add optional label to plot")

    args = parser.parse_args()

    do_plots(args.workdir, title=args.title)
#             do_experimental_systs=args.doExperimentalSysts,
#             do_model_systs=args.doModelSysts)
