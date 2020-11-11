#!/usr/bin/env python

"""Print plots comparing gen-level lambdas across different flavours"""

from __future__ import print_function

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
os.nice(10)
import argparse

# My stuff
from comparator import Contribution, Plot, grab_obj
import qg_common as qgc
import qg_general_plots as qgp
import common_utils as cu


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


# Control output format
OUTPUT_FMT = "pdf"



def do_plots(root_dir):
    # QG variable plots
    pt_bins = qgc.PT_BINS[:]
    var_list = qgc.COMMON_VARS

    radius, pus = cu.get_jet_config_from_dirname(root_dir)
    jet_str = "AK%s" % (radius.upper())

    dy_tfile, qcd_tfile = None, None
    
    if os.path.isfile(os.path.join(root_dir, qgc.DY_FILENAME)):
        dy_tfile = cu.TFileCacher(os.path.join(root_dir, qgc.DY_FILENAME))
    
    if os.path.isfile(os.path.join(root_dir, qgc.QCD_FILENAME)):
        qcd_tfile = cu.TFileCacher(os.path.join(root_dir, qgc.QCD_FILENAME))


    for gr_append in ["", "_groomed"]:
        if gr_append == "_groomed":
            print("Doing groomed plots...")
        else:
            print("Doing ungroomed plots...")

        # GEnJets matched to recojet, reco selection
        zpj_dirname = "ZPlusJets_QG%s" % (gr_append)
        dj_cen_dirname = "Dijet_QG_central_tighter%s" % (gr_append)
        dj_fwd_dirname = "Dijet_QG_forward_tighter%s" % (gr_append)

        # Gen Jets from gen selection
        zpj_dirname = "ZPlusJets_QG_gen%s" % (gr_append)
        dj_cen_dirname = "Dijet_QG_gen_central%s" % (gr_append)
        dj_fwd_dirname = "Dijet_QG_gen_forward%s" % (gr_append)

        for ang in var_list[:]:
            print("...Doing", ang.name)
            var_template = "{prefix}%s_vs_pt" % (ang.var)  # {prefix} is for flavour

            for pt_ind, (start_val, end_val) in enumerate(pt_bins):
                zpj_entries = []
                dijet_cen_entries = []
                dijet_fwd_entries = []

                # Get all plots
                lw = 2
                msize = 1.1
                mc_msize = 1E-3  # small enough to not be seen but not 0 - otherwise the horizontal lines on error bars don't get drawn
                mc_msize = msize  # small enough to not be seen but not 0 - otherwise the horizontal lines on error bars don't get drawn

                flavours_dicts = [
                    {
                        "prefix": "q",
                        "label": "u/d/s"
                    },
                    {
                        "prefix": "g",
                        "label": "g"
                    },
                    {
                        "prefix": "bc",
                        "label": "b/c"
                    },

                ]

                if dy_tfile:
                    # Z+JETS REGION
                    marker = cu.Marker(qgc.DY_MARKER)
                    for flav_dict, color, mark in zip(flavours_dicts, qgc.DY_COLOURS, marker.cycle()):
                        h2d_dyj_mc = dy_tfile.get("%s/%s" % (zpj_dirname, var_template.format(**flav_dict)))
                        dy_kwargs_mc = dict(line_color=color, line_width=lw, fill_color=color,
                                            marker_color=color, marker_style=mark, marker_size=mc_msize,
                                            label="%s flavour" % (flav_dict['label']),
                                            subplot=None)
                        zpj_entries.append((qgp.get_projection_plot(h2d_dyj_mc, start_val, end_val), dy_kwargs_mc))


                if qcd_tfile:
                    # DIJET CENTRAL REGION
                    marker = cu.Marker(qgc.QCD_MARKER)
                    for flav_dict, color, mark in zip(flavours_dicts, qgc.QCD_COLOURS, marker.cycle()):
                        h2d_qcd_cen_mc = qcd_tfile.get("%s/%s" % (dj_cen_dirname, var_template.format(**flav_dict)))
                        qcd_cen_kwargs_mc = dict(line_color=color, line_width=lw, fill_color=color,
                                                 marker_color=color, marker_style=mark, marker_size=mc_msize,
                                                 label="%s flavour" % (flav_dict['label']),
                                                 subplot=None)
                        dijet_mgpy_hist = qgp.get_projection_plot(h2d_qcd_cen_mc, start_val, end_val)
                        dijet_cen_entries.append((dijet_mgpy_hist, qcd_cen_kwargs_mc))

                    # DIJET FORWARD REGION
                    marker = cu.Marker(qgc.QCD_MARKER)
                    for flav_dict, color, mark in zip(flavours_dicts, qgc.QCD_COLOURS, marker.cycle()):
                        h2d_qcd_fwd_mc = qcd_tfile.get("%s/%s" % (dj_fwd_dirname, var_template.format(**flav_dict)))
                        qcd_fwd_kwargs_mc = dict(line_color=color, line_width=lw, fill_color=color,
                                                 marker_color=color, marker_style=mark, marker_size=mc_msize,
                                                 label="%s flavour" % (flav_dict['label']),
                                                 subplot=None)
                        dijet_mgpy_hist = qgp.get_projection_plot(h2d_qcd_fwd_mc, start_val, end_val)
                        dijet_fwd_entries.append((dijet_mgpy_hist, qcd_fwd_kwargs_mc))

                #################
                # SETUP PLOTTING
                #################
                rebin = 2
                v_lower = var_template.lower()
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

                plot_dir = os.path.join(root_dir, "plots_gen_lambda_flav_comparison%s" % (gr_append))

                subplot_title = "Simulation / Data"
                subplot_limits = (0, 2)

                xlabel = ang.name + " (" + ang.lambda_str + ")"
                if gr_append is not "":
                    xlabel = "Groomed " + ang.name + " (" + ang.lambda_str + ")"

                def _title(region_str, start_val, end_val):
                    pt_var_str = "p_{T}^{jet}"
                    s = (("{jet_algo}\n"
                          "{region_label}\n"
                          "{mc_sample}\n"
                          "{bin_edge_low:g} < {pt_str} < {bin_edge_high:g} GeV")
                          .format(
                            jet_algo=jet_str,
                            region_label=region_str,
                            mc_sample="MG5+Pythia8",
                            pt_str=pt_var_str,
                            bin_edge_low=start_val,
                            bin_edge_high=end_val))
                    return s

                has_data = False
                draw_opt = "NOSTACK HIST E1"
                # dj central only
                # dont' try to do multiple signal regions per plot, it looks rubbish
                if len(dijet_cen_entries) > 0:
                    qgp.do_comparison_plot(dijet_cen_entries,
                                           "%s/ptBinned/%s_pt%dto%d_dijet_central.%s" % (plot_dir, ang.var, start_val, end_val, OUTPUT_FMT),
                                           rebin=rebin,
                                           draw_opt=draw_opt,
                                           title=_title(qgc.Dijet_CEN_LABEL, start_val, end_val),
                                           xtitle=xlabel,
                                           xlim='auto' if auto_xlim else xlim, # don't use calc_auto_xlim, since do_comparison_plot will rebin it anyway
                                           ylim=ylim,
                                           data_first=has_data,
                                           mean_rel_error=0.4,
                                           subplot_type=None,
                                           lumi=cu.get_lumi_str(do_dijet=False, do_zpj=True)) # full lumi as just MC

                # dj forward only
                if len(dijet_fwd_entries) > 0:
                    qgp.do_comparison_plot(dijet_fwd_entries,
                                           "%s/ptBinned/%s_pt%dto%d_dijet_forward.%s" % (plot_dir, ang.var, start_val, end_val, OUTPUT_FMT),
                                           rebin=rebin,
                                           draw_opt=draw_opt,
                                           title=_title(qgc.Dijet_FWD_LABEL, start_val, end_val),
                                           xtitle=xlabel,
                                           xlim='auto' if auto_xlim else xlim,
                                           ylim=ylim,
                                           data_first=has_data,
                                           mean_rel_error=0.4,
                                           subplot_type=None,
                                           lumi=cu.get_lumi_str(do_dijet=False, do_zpj=True)) # full lumi as just MC

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
                                           "%s/ptBinned/%s_pt%dto%d_zpj.%s" % (plot_dir, ang.var, start_val, end_val, OUTPUT_FMT),
                                           rebin=rebin,
                                           draw_opt=draw_opt,
                                           title=_title(qgc.ZpJ_LABEL, start_val, end_val),
                                           xtitle=xlabel,
                                           xlim=qgp.calc_auto_xlim([d[0] for d in zpj_entries]) if auto_xlim else xlim,
                                           ylim=ylim,
                                           data_first=has_data,
                                           mean_rel_error=0.4,
                                           subplot_type=None,
                                           lumi=cu.get_lumi_str(do_dijet=False, do_zpj=True))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workdir',
                        help='Workdir with ROOT files to process. '
                        'Each directory must have all the necessary ROOT files')
    args = parser.parse_args()

    do_plots(args.workdir)
