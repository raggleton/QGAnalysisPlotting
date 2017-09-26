"""Common sets of functions to make plots of separation (delta)"""


import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
from array import array
import numpy as np
from collections import OrderedDict

# My stuff
from comparator import Contribution, Plot, grab_obj
import common_utils as cu
from qg_common import *
from qg_general_plots import *

# For debugging
import sys
import tracers
# sys.settrace(tracers.trace_calls)
# sys.settrace(tracers.trace_calls_detail)
# sys.settrace(tracers.trace_calls_and_returns)


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)


def make_roc_graph(hist_signal, hist_background):
    x_ax = hist_signal.GetXaxis()
    x_min = x_ax.GetXmin()
    x_max = x_ax.GetXmax()
    n_bins = x_ax.GetNbins()
    signal_total = hist_signal.Integral()
    background_total = hist_background.Integral()
    signal_effs, background_effs = [], []

    # forward is integral from 0 to x
    signal_cum_forward = hist_signal.GetCumulative(True, "_cumFwd")
    if (signal_cum_forward.GetMaximum() > 0):
        signal_cum_forward.Scale(1./signal_cum_forward.GetMaximum())
    
    background_cum_forward = hist_background.GetCumulative(True, "_cumFwd")
    if (background_cum_forward.GetMaximum() > 0):
        background_cum_forward.Scale(1./background_cum_forward.GetMaximum())

    # backward is from x to 1
    signal_cum_backward = hist_signal.GetCumulative(False, "_cumBwd")
    if (signal_cum_backward.GetMaximum() > 0):
        signal_cum_backward.Scale(1./signal_cum_backward.GetMaximum())
    background_cum_backward = hist_background.GetCumulative(False, "_cumBwd")
    if (background_cum_backward.GetMaximum() > 0):
        background_cum_backward.Scale(1./background_cum_backward.GetMaximum())

    # p = Plot([Contribution(signal_cum_forward, label="signalFwd", line_color=ROOT.kRed),
    #           Contribution(background_cum_forward, label="bkgFwd", line_color=ROOT.kBlue)],
    #           what="hist")
    # p.plot("HIST NOSTACK")
    # p.set_logy()
    # p.save("roc_cum_fwd.pdf")
    # p = Plot([Contribution(signal_cum_backward, label="signalBwd", line_color=ROOT.kRed),
    #           Contribution(background_cum_backward, label="bkgBwd", line_color=ROOT.kBlue)],
    #           what="hist")
    # p.plot("HIST NOSTACK")
    # p.save("roc_cum_bwd.pdf")

    for i in range(2, n_bins+1):
        s_fwd = signal_cum_forward.GetBinContent(i)
        b_fwd = background_cum_forward.GetBinContent(i)

        s_bwd = signal_cum_backward.GetBinContent(i)
        b_bwd = background_cum_backward.GetBinContent(i)

        if s_fwd == 0 and b_fwd == 0:
            signal_effs.append(s_bwd)
            background_effs.append(b_bwd)

        elif s_bwd == 0 and b_bwd == 0:
            signal_effs.append(s_fwd)
            background_effs.append(b_fwd)

        elif b_bwd == 0:
            signal_effs.append(s_bwd)
            background_effs.append(b_bwd)

        elif b_fwd == 0:
            signal_effs.append(s_fwd)
            background_effs.append(b_fwd)

        #  Determine if cutting below or above the value is better
        elif (s_fwd/b_fwd) > (s_bwd/b_bwd):
            signal_effs.append(s_fwd)
            background_effs.append(b_fwd)
        else:
            signal_effs.append(s_bwd)
            background_effs.append(b_bwd)

    # Sort by x axis value so pretty line
    rd = {b:s for b, s in zip(background_effs, signal_effs)}
    sorted_be = sorted(rd.keys())
    sorted_se = [rd[b] for b in sorted_be]

    gr = ROOT.TGraph(len(sorted_be), array('d', sorted_be), array('d', sorted_se))
    return gr


def do_roc_plot(hist_signal, hist_background, output_filename):
    """"Make a single ROC plot"""
    gr = make_roc_graph(hist_signal, hist_background)
    cont = Contribution(gr, marker_style=21)
    p = Plot([cont], "graph", xtitle="#epsilon_{ g}", ytitle="#epsilon_{ q}", xlim=[0, 1], ylim=[0, 1], legend=False)
    p.plot("AL")
    p.save(output_filename)


def do_angularity_roc_plots(sources, plot_dir="roc_angularities_roc", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG",
                             pt_bins=None, var_list=None, var_prepend="", flavour_tag=False):
    """Do roc plots, with different angularities on the same plot"""
    var_list = var_list or COMMON_VARS
    pt_bins = pt_bins or PT_BINS
    zpj_flav = "q" if flavour_tag else ""
    dj_flav = "g" if flavour_tag else ""
    output_append = "_flavMatched" if flavour_tag else ""

    for (start_val, end_val) in pt_bins:
        graph_contribs = []

        for source_ind, source in enumerate(sources):

            # construct a graph of angularities for this source
            for ind, ang in enumerate(var_list, 1):
                v = "%s%s_vs_pt" % (var_prepend, ang.var)

                h2d_dyj = grab_obj(os.path.join(source['root_dir'], qgc.DY_FILENAME), "%s/%s%s" % (source.get('zpj_dirname', zpj_dirname), zpj_flav, v))
                h2d_qcd = grab_obj(os.path.join(source['root_dir'], qgc.QCD_FILENAME), "%s/%s%s" % (source.get('dj_dirname', dj_dirname), dj_flav, v))

                h_dy = get_projection_plot(h2d_dyj, start_val, end_val)
                if (h_dy.Integral() > 0):
                    h_dy.Scale(1./(h_dy.GetBinWidth(1)*h_dy.Integral()))

                h_qcd = get_projection_plot(h2d_qcd, start_val, end_val)
                if (h_qcd.Integral() > 0):
                    h_qcd.Scale(1./(h_qcd.GetBinWidth(1)*h_qcd.Integral()))

                roc_gr = make_roc_graph(h_dy, h_qcd)

                ang_label = "%s [%s]" % (ang.name, ang.lambda_str)

                roc_gr.SetName(source.get("label", ""))
                # if 'line_color' not in source.get('style', {}):
                source['style']['line_color'] = ang.colour

                c = Contribution(roc_gr, label=source.get("label", "")+" "+ang_label,
                                 marker_style=0, marker_color=ang.colour, **source.get("style", {}))
                graph_contribs.append(c)

        basic_gr = ROOT.TGraph(2, array('d', [0, 1]), array('d', [0, 1]))
        graph_contribs.insert(0, Contribution(basic_gr, line_color=ROOT.kBlack, marker_color=ROOT.kBlack, line_style=2, label="y=x"))
        p = Plot(graph_contribs, what='graph', xtitle="#epsilon_{ g}", ytitle="#epsilon_{ q}", title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val),
                 xlim=[1E-1, 1], ylim=[0, 1], legend=True)
        p.legend.SetY1(0.15)
        p.legend.SetY2(0.4)
        p.plot("AL")
        p.save("%s/angularities_roc_pt%dto%d%s.pdf" % (plot_dir, start_val, end_val, output_append))
