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


# Control plot output format
OUTPUT_FMT = "pdf"

TOTAL_LUMI = 35918


def do_plots(root_dir):
    # QG variable plots
    pt_bins = qgc.PT_BINS
    # pt_bins = qgc.THEORY_PT_BINS
    sources = [{"root_dir": root_dir, 'label': "", "style": {'line_style': 1}}]
    var_list = qgc.COMMON_VARS
    var_prepend = ""

    zpj_dirname = "ZPlusJets_QG"
    dj_dirname = "Dijet_QG_tighter"

    for ang in var_list:

        v = "%s%s_vs_pt" % (var_prepend, ang.var)
        for (start_val, end_val) in pt_bins:
            entries = []
            dijet_entries = []
            zpj_entries = []
            # Get all plots
            for ind, source in enumerate(sources):
                lw = 2
                msize = 1.1
                data_line_width = lw

                ####################
                # Z+JETS REGION
                ####################

                # SINGLE MU DATA
                h2d_dyj_data = grab_obj(os.path.join(source['root_dir'], qgc.SINGLE_MU_FILENAME), "%s/%s" % (zpj_dirname, v))
                dy_kwargs_data = dict(line_color=qgc.SINGLE_MU_COLOUR, line_width=data_line_width, fill_color=qgc.SINGLE_MU_COLOUR,
                                      marker_color=qgc.SINGLE_MU_COLOUR, marker_style=qgc.DY_MARKER+ind, marker_size=msize,
                                      label=qgc.SINGLE_MU_LABEL)
                entries.append((qgp.get_projection_plot(h2d_dyj_data, start_val, end_val), dy_kwargs_data))
                zpj_entries.append((qgp.get_projection_plot(h2d_dyj_data, start_val, end_val), dy_kwargs_data))

                # PYTHIA DY MC
                h2d_dyj_mc = grab_obj(os.path.join(source['root_dir'], qgc.DY_FILENAME), "%s/%s" % (zpj_dirname, v))
                dy_kwargs_mc = dict(line_color=qgc.DY_COLOUR, line_width=lw, fill_color=qgc.DY_COLOUR,
                                    marker_color=qgc.DY_COLOUR, marker_style=qgc.DY_MARKER+ind, marker_size=0,
                                    label=qgc.DY_ZpJ_LABEL + " [MG+PY8]")
                entries.append((qgp.get_projection_plot(h2d_dyj_mc, start_val, end_val), dy_kwargs_mc))
                zpj_entries.append((qgp.get_projection_plot(h2d_dyj_mc, start_val, end_val), dy_kwargs_mc))

                # HERWIG++ DY
                h2d_dyj_mc2 = grab_obj(os.path.join(source['root_dir'], qgc.DY_HERWIG_FILENAME), "%s/%s" % (zpj_dirname, v))
                col3 = qgc.DY_COLOURS[2]
                dy_kwargs_mc = dict(line_color=col3, line_width=lw, fill_color=col3,
                                    marker_color=col3, marker_style=qgc.DY_MARKER+ind, marker_size=0,
                                    label=qgc.DY_ZpJ_LABEL + " [H++]")
                entries.append((qgp.get_projection_plot(h2d_dyj_mc2, start_val, end_val), dy_kwargs_mc))
                zpj_entries.append((qgp.get_projection_plot(h2d_dyj_mc2, start_val, end_val), dy_kwargs_mc))

                # MG+HERWIG++ DY
                h2d_dyj_mc3 = grab_obj(os.path.join(source['root_dir'], qgc.DY_MG_HERWIG_FILENAME), "%s/%s" % (zpj_dirname, v))
                col4 = qgc.DY_COLOURS[3]
                dy_kwargs_mc = dict(line_color=col4, line_width=lw, fill_color=col4,
                                    marker_color=col4, marker_style=qgc.DY_MARKER+ind, marker_size=0,
                                    label=qgc.DY_ZpJ_LABEL + " [MG+H++]")
                entries.append((qgp.get_projection_plot(h2d_dyj_mc3, start_val, end_val), dy_kwargs_mc))
                zpj_entries.append((qgp.get_projection_plot(h2d_dyj_mc3, start_val, end_val), dy_kwargs_mc))

                ####################
                # DIJET REGION
                ####################

                # JETHT/ZEROBIAS DATA
                h2d_qcd_data = grab_obj(os.path.join(source['root_dir'], qgc.JETHT_ZB_FILENAME), "%s/%s" % (dj_dirname, v))  # use already merged jetht+zb
                # h2d_zb_data = grab_obj(os.path.join(source['root_dir'], qgc.ZB_FILENAME), "%s/%s" % (dj_dirname, v))
                # h2d_zb_data.Scale(1235009.27580634)
                # h2d_qcd_data.Add(h2d_zb_data)
                qcd_kwargs_data = dict(line_color=qgc.JETHT_COLOUR, line_width=data_line_width, fill_color=qgc.JETHT_COLOUR,
                                       marker_color=qgc.JETHT_COLOUR, marker_style=qgc.QCD_MARKER+ind, marker_size=msize,
                                       label=qgc.JETHT_ZB_LABEL)
                entries.append((qgp.get_projection_plot(h2d_qcd_data, start_val, end_val), qcd_kwargs_data))
                dijet_entries.append((qgp.get_projection_plot(h2d_qcd_data, start_val, end_val), qcd_kwargs_data))

                # MG+PYTHIA QCD MC
                h2d_qcd_mc = grab_obj(os.path.join(source['root_dir'], qgc.QCD_FILENAME), "%s/%s" % (dj_dirname, v))
                qcd_kwargs_mc = dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                                     marker_color=qgc.QCD_COLOUR, marker_style=qgc.QCD_MARKER+ind, marker_size=0,
                                     label=qgc.QCD_Dijet_LABEL + " [MG+PY8]")
                # h2d_qcd_mc.Scale(35860)
                entries.append((qgp.get_projection_plot(h2d_qcd_mc, start_val, end_val), qcd_kwargs_mc))
                dijet_entries.append((qgp.get_projection_plot(h2d_qcd_mc, start_val, end_val), qcd_kwargs_mc))

                # PYTHIA ONLY
                col = qgc.QCD_COLOURS[2]
                h2d_qcd_mc2 = grab_obj(os.path.join(source['root_dir'], qgc.QCD_PYTHIA_ONLY_FILENAME), "%s/%s" % (dj_dirname, v))
                qcd_kwargs_mc2 = dict(line_color=col, line_width=lw, fill_color=col,
                                     marker_color=col, marker_style=qgc.QCD_MARKER+ind, marker_size=0,
                                     label=qgc.QCD_Dijet_LABEL + " [PY8]")
                # h2d_qcd_mc2.Scale(35860)
                entries.append((qgp.get_projection_plot(h2d_qcd_mc2, start_val, end_val), qcd_kwargs_mc2))
                dijet_entries.append((qgp.get_projection_plot(h2d_qcd_mc2, start_val, end_val), qcd_kwargs_mc2))

                # HERWIG++ QCD
                col2 = qgc.QCD_COLOURS[3]
                h2d_qcd_mc3 = grab_obj(os.path.join(source['root_dir'], qgc.QCD_HERWIG_FILENAME), "%s/%s" % (dj_dirname, v))
                qcd_kwargs_mc3 = dict(line_color=col2, line_width=lw, fill_color=col2,
                                     marker_color=col2, marker_style=qgc.QCD_MARKER+ind, marker_size=0,
                                     label=qgc.QCD_Dijet_LABEL + " [H++]")
                h2d_qcd_mc3.Scale(TOTAL_LUMI)
                entries.append((qgp.get_projection_plot(h2d_qcd_mc3, start_val, end_val), qcd_kwargs_mc3))
                dijet_entries.append((qgp.get_projection_plot(h2d_qcd_mc3, start_val, end_val), qcd_kwargs_mc3))


            rebin = 2
            v_lower = v.lower()
            if "multiplicity" in v_lower:
                rebin = 1
            elif "flavour" in v_lower or "thrust" in v_lower or 'ptd' in v_lower:
                rebin = 1

            xlim = None
            if "width" in v_lower or "ptd" in v_lower:
                xlim = (0, 0.75)
            elif"thrust" in v_lower:
                xlim = (0, 0.5)
            elif "multiplicity" in v_lower and "ak4" in sources[0]['root_dir'].lower():
                xlim = (0, 100)
                xlim = (0, 80)

            ylim = None
            if "flavour" in v_lower:
                ylim = (0, 1)
            elif "lha" in v_lower:
                ylim = (0, 5)
                ylim = None

            plot_dir = os.path.join(root_dir, "plots_dy_vs_qcd_mc_vs_data")
            radius, pus = cu.get_jet_config_from_dirname(root_dir)
            jet_str = "AK%s PF %s" % (radius.upper(), pus.upper())
            qgp.do_comparison_plot(entries, "%s/ptBinned/%s_pt%dto%d.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                   rebin=rebin,
                                   title="%d < p_{T}^{jet} < %d GeV\n%s" % (start_val, end_val, jet_str),
                                   xtitle=ang.name + " (" + ang.lambda_str + ")",
                                   xlim=xlim, ylim=ylim,
                                   subplot_type=None,
                                   subplot_title=None,
                                   subplot_limits=(0, 2))

            qgp.do_comparison_plot(dijet_entries, 
                                   "%s/ptBinned/%s_pt%dto%d_dijet.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                   rebin=rebin,
                                   title="%d < p_{T}^{jet} < %d GeV\n%s" % (start_val, end_val, jet_str),
                                   xtitle=ang.name + " (" + ang.lambda_str + ")",
                                   xlim=xlim, ylim=ylim,
                                   subplot_type=None,
                                   subplot_title=None,
                                   subplot_limits=(0, 2))

            qgp.do_comparison_plot(zpj_entries,
                                   "%s/ptBinned/%s_pt%dto%d_zpj.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                   rebin=rebin,
                                   title="%d < p_{T}^{jet} < %d GeV\n%s" % (start_val, end_val, jet_str),
                                   xtitle=ang.name + " (" + ang.lambda_str + ")",
                                   xlim=xlim, ylim=ylim,
                                   subplot_type=None,
                                   subplot_title=None,
                                   subplot_limits=(0, 2))


if __name__ == "__main__":
    parser = qgc.get_parser()
    args = parser.parse_args()

    for workdir in args.workdirs:
        do_plots(workdir)
