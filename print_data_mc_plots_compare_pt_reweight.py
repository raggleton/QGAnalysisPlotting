#!/usr/bin/env python

"""Print main QG dijet plots for a given sample, comparing pythia8 to herwig++,
the latter with and without pt reweighting to match the spectrum

Assumes the standard Herwig QCD file is not reweighted, and there is another
file called uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD_PtReweight.root,
that does have pt reweighting

Use addZeroBiasJetHT.py first to sum jetht and zerobias
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
    dijet_cen_data_hist = None

    if not gen_only:
        # JETHT/ZEROBIAS DATA
        dijet_cen_data_hist = grab_obj(os.path.join(root_dir, qgc.JETHT_ZB_FILENAME), histname)  # use already merged jetht+zb
        qcd_cen_kwargs_data = dict(line_color=qgc.JETHT_COLOUR, line_width=data_line_width, fill_color=qgc.JETHT_COLOUR,
                               marker_color=qgc.JETHT_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=msize,
                               label="Data")
        entries.append((dijet_cen_data_hist, qcd_cen_kwargs_data))

    # MG+PYTHIA QCD MC
    dijet_mgpy_hist = grab_obj(os.path.join(root_dir, qgc.QCD_FILENAME), histname)
    qcd_cen_kwargs_mc = dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                             marker_color=qgc.QCD_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                             label="MG+PY8",
                             subplot=dijet_cen_data_hist)
    entries.append((dijet_mgpy_hist, qcd_cen_kwargs_mc))

    # HERWIG++ QCD
    col2 = qgc.QCD_COLOURS[3]
    dijet_hpp_hist = grab_obj(os.path.join(root_dir, qgc.QCD_HERWIG_FILENAME), histname)
    qcd_cen_kwargs_mc3 = dict(line_color=col2, line_width=lw, fill_color=col2,
                              marker_color=col2, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                              label="H++",
                              subplot=dijet_cen_data_hist)
    dijet_hpp_hist.Scale(TOTAL_LUMI)
    entries.append((dijet_hpp_hist, qcd_cen_kwargs_mc3))

    # HERWIG++ QCD WITH PT REWEIGHT
    col3 = ROOT.kRed
    dijet_hpp_hist_reweight = grab_obj(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD_PtReweight.root"), histname)
    qcd_cen_kwargs_mc4 = dict(line_color=col3, line_width=lw, fill_color=col3, line_style=2,
                              marker_color=col3, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                              label="H++ (p_{T} reweighted)",
                              subplot=dijet_cen_data_hist)
    dijet_hpp_hist_reweight.Scale(TOTAL_LUMI)
    entries.append((dijet_hpp_hist_reweight, qcd_cen_kwargs_mc4))

    # Must do un normed first, since it manipulates original objects
    rebin = 1
    p = qgp.make_comparison_plot_ingredients(entries, rebin=rebin, normalise_hist=False, mean_rel_error=0.4,
                                            title="%s\n%s" % (jet_str, qgc.Dijet_LABEL),
                                            xlim=[0, 4000],
                                            ylim=[1E-1, 1E8],
                                            subplot_type='ratio',
                                            subplot_title="MC / Data",
                                            subplot_limits=(0.5, 1.5))
    output_filename = os.path.join(root_dir, "plots_dijet_vs_qcd_mc_vs_data_pt_reweight", "dj_%s.%s" % (histname.replace("/", "_"), OUTPUT_FMT))
    draw_opt = "NOSTACK HISTE"
    p.plot(draw_opt)
    p.set_logx()
    p.set_logy()
    dirname = os.path.dirname(output_filename)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    p.save(output_filename)

    p = qgp.make_comparison_plot_ingredients(entries, rebin=rebin, normalise_hist=True, mean_rel_error=0.4,
                                            title="%s\n%s" % (jet_str, qgc.Dijet_LABEL),
                                            xlim=[10, 4000],
                                            ylim=[1E-10, 50],
                                            subplot_type='ratio',
                                            subplot_title="MC / Data",
                                            subplot_limits=(0.5, 1.5))
    output_filename = os.path.join(root_dir, "plots_dijet_vs_qcd_mc_vs_data_pt_reweight", "dj_%s_normed.%s" % (histname.replace("/", "_"), OUTPUT_FMT))
    p.plot(draw_opt)
    p.set_logx()
    p.set_logy()
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

        if gr_append == "":
            make_1d_plot(root_dir, "Dijet_gen/pt_dijet_ave", jet_str, gen_only=True)
            make_1d_plot(root_dir, "Dijet_tighter/pt_jet", jet_str)

        dj_cen_dirname = "Dijet_QG_central_tighter%s" % (gr_append)
        dj_fwd_dirname = "Dijet_QG_forward_tighter%s" % (gr_append)

        for ang in var_list[:]:
            print("...Doing", ang.name)
            v = "%s%s_vs_pt" % (var_prepend, ang.var)
            dijet_cen_2d_entries = []
            dijet_fwd_2d_entries = []

            dijet_cen_1d_entries = []
            dijet_fwd_1d_entries = []

            for pt_ind, (start_val, end_val) in enumerate(pt_bins):
                entries = []
                dijet_cen_entries = []
                dijet_fwd_entries = []

                # Get all plots
                lw = 2
                msize = 1.1
                data_line_width = lw

                ####################
                # DIJET CENTRAL REGION
                ####################

                # JETHT/ZEROBIAS DATA
                h2d_qcd_cen_data = grab_obj(os.path.join(root_dir, qgc.JETHT_ZB_FILENAME), "%s/%s" % (dj_cen_dirname, v))  # use already merged jetht+zb
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
                # col = qgc.QCD_COLOURS[2]
                # h2d_qcd_cen_mc2 = grab_obj(os.path.join(root_dir, qgc.QCD_PYTHIA_ONLY_FILENAME), "%s/%s" % (dj_cen_dirname, v))
                # qcd_cen_kwargs_mc2 = dict(line_color=col, line_width=lw, fill_color=col,
                #                           marker_color=col, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                #                           label="PY8",
                #                           subplot=dijet_cen_data_hist)
                # # h2d_qcd_cen_mc2.Scale(35860)
                # entries.append((qgp.get_projection_plot(h2d_qcd_cen_mc2, start_val, end_val), qcd_cen_kwargs_mc2))
                # dijet_cen_entries.append((qgp.get_projection_plot(h2d_qcd_cen_mc2, start_val, end_val), qcd_cen_kwargs_mc2))
                # if pt_ind == 0:
                #     dijet_cen_2d_entries.append((h2d_qcd_cen_mc2, qcd_cen_kwargs_mc2))

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

                # HERWIG++ QCD WITH PT REWEIGHT
                col3 = ROOT.kRed
                h2d_qcd_cen_mc4 = grab_obj(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD_PtReweight.root"), "%s/%s" % (dj_cen_dirname, v))
                qcd_cen_kwargs_mc4 = dict(line_color=col3, line_width=lw, fill_color=col3,
                                          marker_color=col3, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                          label="H++ (p_{T} reweighted)",
                                          subplot=dijet_cen_data_hist)
                h2d_qcd_cen_mc4.Scale(TOTAL_LUMI)
                dijet_hpp_hist_reweight = qgp.get_projection_plot(h2d_qcd_cen_mc4, start_val, end_val)
                entries.append((dijet_hpp_hist_reweight, qcd_cen_kwargs_mc4))
                dijet_cen_entries.append((dijet_hpp_hist_reweight, qcd_cen_kwargs_mc4))
                if pt_ind == 0:
                    dijet_cen_2d_entries.append((h2d_qcd_cen_mc4, qcd_cen_kwargs_mc4))

                ####################
                # DIJET FORWARD REGION
                ####################

                # JETHT/ZEROBIAS DATA
                h2d_qcd_fwd_data = grab_obj(os.path.join(root_dir, qgc.JETHT_ZB_FILENAME), "%s/%s" % (dj_fwd_dirname, v))  # use already merged jetht+zb
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
                # col = qgc.QCD_COLOURS[2]
                # h2d_qcd_fwd_mc2 = grab_obj(os.path.join(root_dir, qgc.QCD_PYTHIA_ONLY_FILENAME), "%s/%s" % (dj_fwd_dirname, v))
                # qcd_fwd_kwargs_mc2 = dict(line_color=col, line_width=lw, fill_color=col,
                #                           marker_color=col, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                #                           label="PY8",
                #                           subplot=dijet_fwd_data_hist)
                # # h2d_qcd_fwd_mc2.Scale(35860)
                # entries.append((qgp.get_projection_plot(h2d_qcd_fwd_mc2, start_val, end_val), qcd_fwd_kwargs_mc2))
                # dijet_fwd_entries.append((qgp.get_projection_plot(h2d_qcd_fwd_mc2, start_val, end_val), qcd_fwd_kwargs_mc2))
                # if pt_ind == 0:
                #     dijet_fwd_2d_entries.append((h2d_qcd_fwd_mc2, qcd_fwd_kwargs_mc2))

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

                col3 = ROOT.kRed
                h2d_qcd_fwd_mc4 = grab_obj(os.path.join(root_dir, "uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD_PtReweight.root"), "%s/%s" % (dj_fwd_dirname, v))
                qcd_fwd_kwargs_mc4 = dict(line_color=col3, line_width=lw, fill_color=col3,
                                          marker_color=col3, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                          label="H++ (p_{T} reweighted)",
                                          subplot=dijet_fwd_data_hist)
                h2d_qcd_fwd_mc4.Scale(TOTAL_LUMI)
                dijet_hpp_hist_reweight = qgp.get_projection_plot(h2d_qcd_fwd_mc4, start_val, end_val)
                entries.append((dijet_hpp_hist_reweight, qcd_fwd_kwargs_mc4))
                dijet_fwd_entries.append((dijet_hpp_hist_reweight, qcd_fwd_kwargs_mc4))
                if pt_ind == 0:
                    dijet_fwd_2d_entries.append((h2d_qcd_fwd_mc4, qcd_fwd_kwargs_mc4))

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

                plot_dir = os.path.join(root_dir, "plots_dijet_vs_qcd_mc_vs_data%s_pt_reweight" % (gr_append))

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


if __name__ == "__main__":
    parser = qgc.get_parser()
    args = parser.parse_args()

    for workdir in args.workdirs:
        do_plots(workdir)
