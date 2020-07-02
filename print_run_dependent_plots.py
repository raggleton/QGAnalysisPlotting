#!/usr/bin/env python

"""
Plot stuff by Run period e.g. prefiring checks
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

TOTAL_LUMI = 35918


def do_plots(root_dir):
    # QG variable plots
    pt_bins = qgc.PT_BINS[:]
    var_list = qgc.COMMON_VARS
    var_prepend = ""

    radius, pus = cu.get_jet_config_from_dirname(root_dir)
    jet_str = "AK%s PF %s" % (radius.upper(), pus.upper())

    for gr_append in ["", "_groomed"][:1]:
        if gr_append == "_groomed":
            print("Doing groomed plots...")
        else:
            print("Doing ungroomed plots...")

        zpj_dirname = "ZPlusJets_QG%s" % (gr_append)
        dj_cen_dirname = "Dijet_QG_central_tighter%s" % (gr_append)
        dj_fwd_dirname = "Dijet_QG_forward_tighter%s" % (gr_append)

        # tfiles = {
        #     "data_run_B": cu.open_root_file(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunB.root")),
        #     "data_run_C": cu.open_root_file(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunC.root")),
        #     "data_run_D": cu.open_root_file(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunD.root")),
        #     "data_run_E": cu.open_root_file(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunE.root")),
        #     "data_run_F": cu.open_root_file(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunF.root")),
        #     "data_run_G": cu.open_root_file(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunG.root")),
        #     "data_run_H": cu.open_root_file(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunH.root")),
        # }

        tfiles = {
            "data_run_B": cu.open_root_file(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_ZeroBias_RunB.root")),
            "data_run_C": cu.open_root_file(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_ZeroBias_RunC.root")),
            "data_run_D": cu.open_root_file(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_ZeroBias_RunD.root")),
            "data_run_E": cu.open_root_file(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_ZeroBias_RunE.root")),
            "data_run_F": cu.open_root_file(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_ZeroBias_RunF.root")),
            "data_run_G": cu.open_root_file(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_ZeroBias_RunG.root")),
            "data_run_H": cu.open_root_file(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_ZeroBias_RunH.root")),
        }

        for ang in var_list[:]:
            print("...Doing", ang.name)
            v = "%s%s_vs_pt" % (var_prepend, ang.var)

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

                colours = [
                    ROOT.kBlack,
                    ROOT.kRed,
                    ROOT.kBlue,
                    ROOT.kGreen+2,
                    ROOT.kOrange+1,
                    ROOT.kAzure+1,
                    ROOT.kGreen+3
                ]

                # --------------------------------------------------------------
                # Create reference hists, for Run G+H
                # --------------------------------------------------------------
                ####################
                # DIJET CENTRAL REGION
                ####################
                h2d_qcd_cen_data_g = cu.get_from_tfile(tfiles['data_run_G'], "%s/%s" % (dj_cen_dirname, v))
                h2d_qcd_cen_data_h = cu.get_from_tfile(tfiles['data_run_H'], "%s/%s" % (dj_cen_dirname, v))
                colour = colours[0]
                qcd_cen_kwargs_data = dict(line_color=colour, line_width=data_line_width, fill_color=colour,
                                           marker_color=colour, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                           label="Runs G+H")
                dijet_cen_data_hist_g = qgp.get_projection_plot(h2d_qcd_cen_data_g, start_val, end_val)
                dijet_cen_data_hist_h = qgp.get_projection_plot(h2d_qcd_cen_data_h, start_val, end_val)
                dijet_cen_data_hist_g.Add(dijet_cen_data_hist_h)
                dijet_cen_entries.append((dijet_cen_data_hist_g, qcd_cen_kwargs_data))

                ####################
                # DIJET FORWARD REGION
                ####################
                h2d_qcd_fwd_data_g = cu.get_from_tfile(tfiles['data_run_G'], "%s/%s" % (dj_fwd_dirname, v))
                h2d_qcd_fwd_data_h = cu.get_from_tfile(tfiles['data_run_H'], "%s/%s" % (dj_fwd_dirname, v))
                qcd_fwd_kwargs_data = dict(line_color=colour, line_width=data_line_width, fill_color=colour,
                                           marker_color=colour, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                           label="Runs G+H")
                dijet_fwd_data_hist_g = qgp.get_projection_plot(h2d_qcd_fwd_data_g, start_val, end_val)
                dijet_fwd_data_hist_h = qgp.get_projection_plot(h2d_qcd_fwd_data_h, start_val, end_val)
                dijet_fwd_data_hist_g.Add(dijet_fwd_data_hist_h)
                dijet_fwd_entries.append((dijet_fwd_data_hist_g, qcd_fwd_kwargs_data))


                # --------------------------------------------------------------
                # FOR SUMMED B-F:
                # --------------------------------------------------------------
                dijet_cen_data_hist_btof, dijet_fwd_data_hist_btof = None, None
                runs = ["B", "C", "D", "E", "F"]
                for run in runs:

                    # DIJET CENTRAL REGION
                    h2d_qcd_cen_data = cu.get_from_tfile(tfiles['data_run_%s' % run], "%s/%s" % (dj_cen_dirname, v))
                    this_dijet_cen_data_hist = qgp.get_projection_plot(h2d_qcd_cen_data, start_val, end_val)
                    if dijet_cen_data_hist_btof is None:
                        dijet_cen_data_hist_btof = this_dijet_cen_data_hist.Clone()
                    else:
                        dijet_cen_data_hist_btof.Add(this_dijet_cen_data_hist)


                    # DIJET FORWARD REGION
                    h2d_qcd_fwd_data = cu.get_from_tfile(tfiles['data_run_%s' % run], "%s/%s" % (dj_fwd_dirname, v))  # use already merged jetht+zb
                    this_dijet_fwd_data_hist = qgp.get_projection_plot(h2d_qcd_fwd_data, start_val, end_val)
                    if dijet_fwd_data_hist_btof is None:
                        dijet_fwd_data_hist_btof = this_dijet_fwd_data_hist.Clone()
                    else:
                        dijet_fwd_data_hist_btof.Add(this_dijet_fwd_data_hist)

                colour = colours[1]
                mark = cu.Marker.get('circle')
                btof_cen_kwargs_data = dict(line_color=colour, line_width=data_line_width, fill_color=colour,
                                            marker_color=colour, marker_style=mark, marker_size=msize,
                                            label="Runs %s-%s" % (runs[0], runs[-1]),
                                            subplot=dijet_cen_data_hist_g)
                dijet_cen_entries.append((dijet_cen_data_hist_btof, btof_cen_kwargs_data))

                btof_fwd_kwargs_data = dict(line_color=colour, line_width=data_line_width, fill_color=colour,
                                            marker_color=colour, marker_style=mark, marker_size=msize,
                                            label="Runs %s-%s" % (runs[0], runs[-1]),
                                            subplot=dijet_fwd_data_hist_g)
                dijet_fwd_entries.append((dijet_fwd_data_hist_btof, btof_fwd_kwargs_data))

                # FOR INDIVIDUAL RUN PERIODS:
                # runs = ["B", "C", "D", "E", "F"]
                # for run, colour, mark in zip(runs, colours[1:], cu.Marker().cycle()):

                #     ####################
                #     # DIJET CENTRAL REGION
                #     ####################
                #     # JETHT DATA
                #     h2d_qcd_cen_data = cu.get_from_tfile(tfiles['data_run_%s' % run], "%s/%s" % (dj_cen_dirname, v))
                #     qcd_cen_kwargs_data = dict(line_color=colour, line_width=data_line_width, fill_color=colour,
                #                                marker_color=colour, marker_style=mark, marker_size=msize,
                #                                label="Run %s" % run,
                #                                subplot=dijet_cen_data_hist_g)
                #     dijet_cen_data_hist = qgp.get_projection_plot(h2d_qcd_cen_data, start_val, end_val)
                #     dijet_cen_entries.append((dijet_cen_data_hist, qcd_cen_kwargs_data))

                #     ####################
                #     # DIJET FORWARD REGION
                #     ####################
                #     # JETHT DATA
                #     h2d_qcd_fwd_data = cu.get_from_tfile(tfiles['data_run_%s' % run], "%s/%s" % (dj_fwd_dirname, v))  # use already merged jetht+zb
                #     qcd_fwd_kwargs_data = dict(line_color=colour, line_width=data_line_width, fill_color=colour,
                #                                marker_color=colour, marker_style=mark, marker_size=msize,
                #                                label="Run %s" % run,
                #                                subplot=dijet_fwd_data_hist_g)
                #     dijet_fwd_data_hist = qgp.get_projection_plot(h2d_qcd_fwd_data, start_val, end_val)
                #     dijet_fwd_entries.append((dijet_fwd_data_hist, qcd_fwd_kwargs_data))

                #     ####################
                #     # Z+JETS REGION
                #     ####################

                #     # SINGLE MU DATA
                #     # h2d_dyj_data = cu.get_from_tfile(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.DATA.Data_SingleMu_Run%s.root" % run), "%s/%s" % (zpj_dirname, v))
                #     # dy_kwargs_data = dict(line_color=qgc.colour, line_width=data_line_width, fill_color=qgc.colour,
                #     #                       marker_color=qgc.colour, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=msize,
                #     #                       label="SingleMu Run %s" % run)
                #     # zpj_data_hist = qgp.get_projection_plot(h2d_dyj_data, start_val, end_val)
                #     # zpj_entries.append((zpj_data_hist, dy_kwargs_data))

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

                plot_dir = os.path.join(root_dir, "plots_run_dependent%s" % (gr_append))
                plot_dir = os.path.join(root_dir, "plots_run_dependent_btof%s" % (gr_append))

                subplot_title = "Run * / G+H"
                subplot_title = "B-F / G+H"
                subplot_limits = (0.5, 1.5)

                xlabel = ang.name + " (" + ang.lambda_str + ")"
                if gr_append is not "":
                    xlabel = "Groomed " + ang.name + " (" + ang.lambda_str + ")"

                # dj central only
                dijet_cen_entries = dijet_cen_entries[1:] + dijet_cen_entries[0:1]  # move G+H to last
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
                dijet_fwd_entries = dijet_fwd_entries[1:] + dijet_fwd_entries[0:1]  # move G+H to last
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

                # # zpj only
                # qgp.do_comparison_plot(zpj_entries,
                #                        "%s/ptBinned/%s_pt%dto%d_zpj.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                #                        rebin=rebin,
                #                        title="%d < p_{T}^{jet} < %d GeV\n%s\n%s" % (start_val, end_val, jet_str, qgc.ZpJ_LABEL),
                #                        xtitle=xlabel,
                #                        xlim=xlim,
                #                        ylim=ylim,
                #                        subplot_type='ratio',
                #                        subplot_title=subplot_title,
                #                        subplot_limits=subplot_limits)


if __name__ == "__main__":
    parser = qgc.get_parser()
    args = parser.parse_args()

    for workdir in args.workdirs:
        do_plots(workdir)
