#!/usr/bin/env python

"""Print main QG plots for a given sample, comparing with/without Z k factors

The standard DYJetsToLL file is assuemd to be weighted

Requires that the unweighted file be called uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_noZReweight.root
in the smae dir
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


def make_1d_plot(root_dir, histname, jet_str, gen_only=False):
    entries = []
    lw = 2
    data_line_width = lw
    msize = 1.1
    zpj_data_hist = None
    if not gen_only:
        zpj_data_hist = grab_obj(os.path.join(root_dir, qgc.SINGLE_MU_FILENAME), histname)
        SINGLE_MU_COLOUR = ROOT.kBlack
        dy_kwargs_data = dict(line_color=SINGLE_MU_COLOUR, line_width=2, fill_color=SINGLE_MU_COLOUR,
                              marker_color=SINGLE_MU_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=msize,
                              label="Data")
        entries.append((zpj_data_hist, dy_kwargs_data))

    # PYTHIA DY MC WITH k FACTORS
    h_dyj_mc = grab_obj(os.path.join(root_dir, qgc.DY_FILENAME), histname)
    dy_kwargs_mc = dict(line_color=qgc.DY_COLOUR, line_width=lw, fill_color=qgc.DY_COLOUR,
                        marker_color=qgc.DY_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=0,
                        label="MG+PY8 (with k-factor)",
                        subplot=zpj_data_hist)
    entries.append((h_dyj_mc, dy_kwargs_mc))

    # WITHOUT K FACTORS
    h_dyj_mc_noWeight = grab_obj(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_noZReweight.root"), histname)
    dy_kwargs_mc = dict(line_color=qgc.QCD_COLOURS[0], line_width=lw, fill_color=qgc.QCD_COLOURS[0],
                        marker_color=qgc.QCD_COLOURS[0], marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=0,
                        label="MG+PY8 (without k-factor)",
                        subplot=zpj_data_hist)
    entries.append((h_dyj_mc_noWeight, dy_kwargs_mc))

    # Must do un normed first, since it manipulates original objects
    rebin = 10
    p = qgp.make_comparison_plot_ingredients(entries, rebin=rebin, normalise_hist=False, mean_rel_error=0.4,
                                            title="%s\n%s" % (jet_str, qgc.ZpJ_LABEL),
                                            xlim=[30, 4000],
                                            ylim=[1E-1, 1E8],
                                            subplot_type='ratio',
                                            subplot_title="MC / Data",
                                            subplot_limits=(0.5, 1.5))
    output_filename = os.path.join(root_dir, "plots_dy_vs_qcd_mc_vs_data_z_reweight", "zpj_%s.%s" % (histname.replace("/", "_"), OUTPUT_FMT))
    draw_opt = "NOSTACK HISTE"
    p.plot(draw_opt)
    p.set_logx()
    p.set_logy(do_more_labels=False)
    dirname = os.path.dirname(output_filename)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    p.save(output_filename)

    p = qgp.make_comparison_plot_ingredients(entries, rebin=1, normalise_hist=True, mean_rel_error=0.4,
                                            title="%s\n%s" % (jet_str, qgc.ZpJ_LABEL),
                                            xlim=[30, 4000],
                                            ylim=[1E-6, 50],
                                            subplot_type='ratio',
                                            subplot_title="MC / Data",
                                            subplot_limits=(0.5, 1.5))
    output_filename = os.path.join(root_dir, "plots_dy_vs_qcd_mc_vs_data_z_reweight", "zpj_%s_normed.%s" % (histname.replace("/", "_"), OUTPUT_FMT))
    p.plot(draw_opt)
    p.set_logx()
    p.set_logy(do_more_labels=False)
    p.save(output_filename)


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

        # do a pt plot just for sanity
        if gr_append == "":
            make_1d_plot(root_dir, "ZPlusJets_gen/pt_jet1", jet_str, gen_only=True)
            make_1d_plot(root_dir, "ZPlusJets_gen/pt_mumu", jet_str, gen_only=True)
            make_1d_plot(root_dir, "ZPlusJets/pt_jet1", jet_str)
            make_1d_plot(root_dir, "ZPlusJets/pt_mumu", jet_str)

        continue

        zpj_dirname = "ZPlusJets_QG%s" % (gr_append)

        for ang in var_list[:]:
            print("...Doing", ang.name)
            v = "%s%s_vs_pt" % (var_prepend, ang.var)
            zpj_2d_entries = []

            zpj_1d_entries = []

            for pt_ind, (start_val, end_val) in enumerate(pt_bins):
                entries = []
                zpj_entries = []

                # Get all plots
                lw = 2
                msize = 1.1
                data_line_width = 0

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

                # PYTHIA DY MC WITH k FACTORS
                h2d_dyj_mc = grab_obj(os.path.join(root_dir, qgc.DY_FILENAME), "%s/%s" % (zpj_dirname, v))
                dy_kwargs_mc = dict(line_color=qgc.DY_COLOUR, line_width=lw, fill_color=qgc.DY_COLOUR,
                                    marker_color=qgc.DY_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=0,
                                    label="MG+PY8 (w/ k-factor)",
                                    subplot=zpj_data_hist)
                entries.append((qgp.get_projection_plot(h2d_dyj_mc, start_val, end_val), dy_kwargs_mc))
                zpj_entries.append((qgp.get_projection_plot(h2d_dyj_mc, start_val, end_val), dy_kwargs_mc))
                if pt_ind == 0:
                    zpj_2d_entries.append((h2d_dyj_mc, dy_kwargs_mc))

                # WIHTOUT K FACTORS
                h2d_dyj_mc_noWeight = grab_obj(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_noZReweight.root"), "%s/%s" % (zpj_dirname, v))
                dy_kwargs_mc_noWeight = dict(line_color=qgc.DY_COLOURS[0], line_width=lw, fill_color=qgc.DY_COLOURS[0], line_style=2,
                                    marker_color=qgc.DY_COLOURS[0], marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=0,
                                    label="MG+PY8 (w/out k-factor)",
                                    subplot=zpj_data_hist)
                entries.append((qgp.get_projection_plot(h2d_dyj_mc_noWeight, start_val, end_val), dy_kwargs_mc_noWeight))
                zpj_entries.append((qgp.get_projection_plot(h2d_dyj_mc_noWeight, start_val, end_val), dy_kwargs_mc_noWeight))
                if pt_ind == 0:
                    zpj_2d_entries.append((h2d_dyj_mc_noWeight, dy_kwargs_mc_noWeight))

                # HERWIG++ DY
                # if end_val < 255:
                #     h2d_dyj_mc2 = grab_obj(os.path.join(root_dir, qgc.DY_HERWIG_FILENAME), "%s/%s" % (zpj_dirname, v))
                #     col3 = qgc.DY_COLOURS[2]
                #     dy_kwargs_mc2 = dict(line_color=col3, line_width=lw, fill_color=col3,
                #                          marker_color=col3, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=0,
                #                          label="H++",
                #                          subplot=zpj_data_hist)
                #     entries.append((qgp.get_projection_plot(h2d_dyj_mc2, start_val, end_val), dy_kwargs_mc2))
                #     zpj_entries.append((qgp.get_projection_plot(h2d_dyj_mc2, start_val, end_val), dy_kwargs_mc2))
                #     if pt_ind == 0:
                #         zpj_2d_entries.append((h2d_dyj_mc2, dy_kwargs_mc2))
                # else:
                #     zpj_entries.append(None)

                # MG+HERWIG++ DY
                # if end_val < 151:
                #     h2d_dyj_mc3 = grab_obj(os.path.join(root_dir, qgc.DY_MG_HERWIG_FILENAME), "%s/%s" % (zpj_dirname, v))
                #     col4 = qgc.DY_COLOURS[3]
                #     dy_kwargs_mc3 = dict(line_color=col4, line_width=lw, fill_color=col4,
                #                          marker_color=col4, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=0,
                #                          label="MG+H++",
                #                          subplot=zpj_data_hist)
                #     entries.append((qgp.get_projection_plot(h2d_dyj_mc3, start_val, end_val), dy_kwargs_mc3))
                #     zpj_entries.append((qgp.get_projection_plot(h2d_dyj_mc3, start_val, end_val), dy_kwargs_mc3))
                #     if pt_ind == 0:
                #         zpj_2d_entries.append((h2d_dyj_mc3, dy_kwargs_mc3))
                # else:
                #     zpj_entries.append(None)

                zpj_1d_entries.append(zpj_entries)

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

                plot_dir = os.path.join(root_dir, "plots_dy_vs_qcd_mc_vs_data%s_z_reweight" % (gr_append))

                subplot_title = "MC / Data"
                subplot_limits = (0.5, 1.5)

                xlabel = ang.name + " (" + ang.lambda_str + ")"
                if gr_append is not "":
                    xlabel = "Groomed " + ang.name + " (" + ang.lambda_str + ")"

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


            # zpj_1d_entries[i][j] is the jth sample in the ith pt bin
            qgp.do_mean_rms_summary_plot(zpj_1d_entries[:],
                                         pt_bins[:],
                                         "%s/ptBinned/%s_box_zpj_mpl.%s" % (plot_dir, v, OUTPUT_FMT),
                                         var_label=var_label,
                                         xlim=(50, 614),
                                         ylim_mean=ylim_mean,
                                         ylim_rms=ylim_rms,
                                         region_title="%s jets in Z+jets" % (jet_str))


if __name__ == "__main__":
    parser = qgc.get_parser()
    args = parser.parse_args()

    for workdir in args.workdirs:
        do_plots(workdir)
