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


# Control plot output format
OUTPUT_FMT = "pdf"


def do_plots(root_dir):
    # Kinematic plots

    # QG variable plots
    pt_bins = qgc.PT_BINS
    pt_bins = qgc.THEORY_PT_BINS
    sources = sources = [{"root_dir": root_dir, 'label': "", "style": {'line_style': 1}}]
    var_list = qgc.COMMON_VARS
    var_prepend = ""

    zpj_dirname = "ZPlusJets_QG"
    dj_dirname = "Dijet_QG"

    for ang in var_list:

        v = "%s%s_vs_pt" % (var_prepend, ang.var)
        for (start_val, end_val) in pt_bins:
            entries = []

            # Get all plots
            for ind, source in enumerate(sources):
                lw = 2
                msize = 1.1

                # MC sources
                h2d_dyj_mc = grab_obj(os.path.join(source['root_dir'], qgc.DY_FILENAME), "%s/%s" % (zpj_dirname, v))
                dy_kwargs_mc = dict(line_color=qgc.DY_COLOUR, line_width=lw, fill_color=qgc.DY_COLOUR,
                                    marker_color=qgc.DY_COLOUR, marker_style=qgc.DY_MARKER+ind, marker_size=0,
                                    label=qgc.DY_ZpJ_LABEL)
                entries.append((qgp.get_projection_plot(h2d_dyj_mc, start_val, end_val), dy_kwargs_mc))

                h2d_qcd_mc = grab_obj(os.path.join(source['root_dir'], qgc.QCD_FILENAME), "%s/%s" % (dj_dirname, v))
                qcd_kwargs_mc = dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                                     marker_color=qgc.QCD_COLOUR, marker_style=qgc.QCD_MARKER+ind, marker_size=0,
                                     label=qgc.QCD_Dijet_LABEL)
                entries.append((qgp.get_projection_plot(h2d_qcd_mc, start_val, end_val), qcd_kwargs_mc))

                # Data sources
                h2d_dyj_data = grab_obj(os.path.join(source['root_dir'], qgc.SINGLE_MU_FILENAME), "%s/%s" % (zpj_dirname, v))
                dy_kwargs_data = dict(line_color=qgc.SINGLE_MU_COLOUR, line_width=0, fill_color=qgc.SINGLE_MU_COLOUR,
                                      marker_color=qgc.SINGLE_MU_COLOUR, marker_style=qgc.DY_MARKER+ind, marker_size=msize,
                                      label=qgc.SINGLE_MU_LABEL)
                entries.append((qgp.get_projection_plot(h2d_dyj_data, start_val, end_val), dy_kwargs_data))

                h2d_qcd_data = grab_obj(os.path.join(source['root_dir'], qgc.JETHT_FILENAME), "%s/%s" % (dj_dirname, v))
                qcd_kwargs_data = dict(line_color=qgc.JETHT_COLOUR, line_width=0, fill_color=qgc.JETHT_COLOUR,
                                       marker_color=qgc.JETHT_COLOUR, marker_style=qgc.QCD_MARKER+ind, marker_size=msize,
                                       label=qgc.JETHT_LABEL)
                entries.append((qgp.get_projection_plot(h2d_qcd_data, start_val, end_val), qcd_kwargs_data))

            rebin = 2
            if "multiplicity" in v:
                rebin = 2
            elif "flavour" in v or "thrust" in v or 'pTD' in v:
                rebin = 1

            xlim = None
            if "width" in v or "thrust" in v or "pTD" in v:
                xlim = (0, 0.5)
            elif "multiplicity" in v.lower() and "ak4" in sources[0]['root_dir'].lower():
                xlim = (0, 100)

            ylim = None
            if "flavour" in v:
                ylim = (0, 1)
            elif "LHA" in v:
                ylim = (0, 5)

            plot_dir = os.path.join(root_dir, "plots_dy_vs_qcd_mc_vs_data")
            qgp.do_comparison_plot(entries, "%s/ptBinned/%s_pt%dto%d.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                   rebin=rebin,
                                   title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val),
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
