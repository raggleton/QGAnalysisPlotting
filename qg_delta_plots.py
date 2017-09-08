"""Common sets of functions to make plots of separation (delta)"""


import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
from array import array

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


def get_ddelta_plot(one_hist, other_hist):
    """Make ddelta/dX plot

    Calculated as:
    0.5 * (x-y)^2 / (x+y)
    """
    new_hist = one_hist.Clone(ROOT.TUUID().AsString())
    new_hist.Add(other_hist, -1.)
    new_hist.Multiply(new_hist)
    sum_hist = one_hist.Clone(ROOT.TUUID().AsString())
    sum_hist.Add(other_hist)
    new_hist.Divide(sum_hist)
    new_hist.Scale(0.5)
    return new_hist


def calculate_delta(one_hist, other_hist):
    dist = get_ddelta_plot(one_hist, other_hist)
    return dist.Integral()*dist.GetBinWidth(1)


def calculate_delta(diff_hist):
    return diff_hist.Integral()*diff_hist.GetBinWidth(1)


def plot_ddelta(ddelta_hist, output_filename, xtitle, ytitle, title=""):
    cont = Contribution(ddelta_hist)
    p = Plot([cont], what="hist", legend=None, xtitle=xtitle, ytitle=ytitle, title=title)
    p.plot("HISTE")
    p.save(output_filename)


def construct_deltas_graph(deltas):
    N = len(deltas)
    gr = ROOT.TGraphErrors(N, array('d', range(N)), array('d', deltas), array('d', [0.5]*N), array('d', [0.00001]*N))
    return gr


def do_deltas_plot(graph_contribs, output_filename, bin_labels, title="", xtitle=""):
    do_legend = len(graph_contribs) > 1
    p = Plot(graph_contribs, what="graph", title=title, xtitle=xtitle, ytitle="Separation #Delta", legend=do_legend)
    if do_legend:
        p.legend.SetX1(0.7)
        p.legend.SetY1(0.15)
        p.legend.SetY2(0.35)
    p.plot("AP")
    # p.container.GetXaxis().LabelsOption("h")
    xax = p.container.GetXaxis()
    for i, lab in enumerate(bin_labels):
        bin_ind = xax.FindBin(i)  # need this as they don't correspond at all!
        p.container.GetXaxis().SetBinLabel(bin_ind, lab)
    p.container.SetMinimum(0)
    p.save(output_filename)


def do_pt_min_delta_plots(sources, plot_dir="deltas_ptmin", 
                          zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG",
                          var_list=None, var_prepend="", flavour_tag=False, 
                          save_component_hists=False, ofmt="pdf"):
    """Do plots comparing power of different ptMin cuts"""
    var_list = var_list or COMMON_VARS
    ptmin_bins = [50, 100, 200, 400, 800][:-1]
    zpj_flav = "q" if flavour_tag else ""
    dj_flav = "g" if flavour_tag else ""
    output_append = "_flavMatched" if flavour_tag else ""
    for ang in var_list:
        v = "%s%s_vs_pt" % (var_prepend, ang.var)
        graph_contribs, bin_labels = [], []

        for source_ind, source in enumerate(sources):
            deltas, conts = [], []

            for ind, pt_min in enumerate(ptmin_bins, 1):
                h2d_dyj = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % source['root_dir'], "%s_ptMin_%d/%s%s" % (source.get('zpj_dirname', zpj_dirname), pt_min, zpj_flav, v))
                h2d_qcd = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % source['root_dir'], "%s_ptMin_%d/%s%s" % (source.get('dj_dirname', dj_dirname), pt_min, dj_flav, v))
                start_val, end_val = 80, 2000
                h_dy = get_projection_plot(h2d_dyj, start_val, end_val)
                if (h_dy.Integral()>0):
                    h_dy.Scale(1./(h_dy.GetBinWidth(1)*h_dy.Integral()))

                h_qcd = get_projection_plot(h2d_qcd, start_val, end_val)
                if (h_qcd.Integral()>0):
                    h_qcd.Scale(1./(h_qcd.GetBinWidth(1)*h_qcd.Integral()))

                ddelta_hist = get_ddelta_plot(h_dy, h_qcd)
                conts.append(Contribution(ddelta_hist, line_width=1, line_style=ind, line_color=(ind*10)+44, fill_color=(ind*10)+44, label="p_{T}^{Min} = %d GeV" % pt_min, rebin_hist=2))
                deltas.append(calculate_delta(ddelta_hist))

                if source_ind == 0:
                    bin_labels.append("%d" % pt_min)

                if save_component_hists:
                    plot_ddelta(ddelta_hist, "%s/%s/%s_ddelta_ptMin_%d%s.%s" % (source['root_dir'], plot_dir, ang.var, pt_min, output_append, ofmt),
                                xtitle=ang.name + " (" + ang.lambda_str + ")", ytitle="d#Delta/d" + ang.lambda_str)

            if save_component_hists:
                p = Plot(conts, what="hist", xtitle=ang.name, ytitle="p.d.f")
                p.plot("NOSTACK HISTE")
                p.save("%s/%s/%s_ddelta_ptMin_comparison%s.%s" % (source['root_dir'], plot_dir, ang.var, output_append, ofmt))

            gr = construct_deltas_graph(deltas)
            gr.SetName(source.get("label", ""))
            if 'style' in source and 'line_width' not in source['style']:
                source['style']['line_width'] = 2
            c = Contribution(gr, label=source.get("label", ""), marker_style=0, **source.get("style", {}))
            graph_contribs.append(c)

        do_deltas_plot(graph_contribs, 
                       "%s/ptMins_%s%s.%s" % (plot_dir, ang.var, output_append, ofmt),
                       bin_labels=bin_labels, 
                       title="%s [%s]" % (ang.name, ang.lambda_str), 
                       xtitle="p_{T}^{min} [GeV]")


def do_angularity_delta_plots(sources, plot_dir="delta_angularities", 
                              zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG", 
                              var_list=None, var_prepend="", pt_bins=None, 
                              flavour_tag=False, save_component_hists=False, ofmt="pdf"):
    """Do plots comparing power of different angularities"""
    var_list = var_list or COMMON_VARS
    pt_bins = pt_bins or PT_BINS
    h2d_dyj = None
    h2d_qcd = None
    zpj_flav = "q" if flavour_tag else ""
    dj_flav = "g" if flavour_tag else ""
    output_append = "_flavMatched" if flavour_tag else ""
    for (start_val, end_val) in pt_bins:
        graph_contribs, bin_labels = [], []

        for source_ind, source in enumerate(sources):
            deltas = []

            # construct a graph of angularities for this source
            for ind, ang in enumerate(var_list, 1):
                v = "%s%s_vs_pt" % (var_prepend, ang.var)

                h2d_dyj = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % source['root_dir'], "%s/%s%s" % (source.get('zpj_dirname', zpj_dirname), zpj_flav, v))
                h2d_qcd = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % source['root_dir'], "%s/%s%s" % (source.get('dj_dirname', dj_dirname), dj_flav, v))

                h_dy = get_projection_plot(h2d_dyj, start_val, end_val)
                if (h_dy.Integral() > 0):
                    h_dy.Scale(1./(h_dy.GetBinWidth(1)*h_dy.Integral()))

                h_qcd = get_projection_plot(h2d_qcd, start_val, end_val)
                if (h_qcd.Integral() > 0):
                    h_qcd.Scale(1./(h_qcd.GetBinWidth(1)*h_qcd.Integral()))

                ddelta_hist = get_ddelta_plot(h_dy, h_qcd)
                deltas.append(calculate_delta(ddelta_hist))

                if source_ind == 0:
                    bin_labels.append("#splitline{%s}{%s}" % (ang.name, ang.lambda_str))

                if save_component_hists:
                    plot_ddelta(ddelta_hist, "%s/%s/angularities_pt%dto%d_ddelta_%s%s.%s" % (source['root_dir'], plot_dir, start_val, end_val, ang.var, output_append, ofmt),
                                xtitle=ang.name + " (" + ang.lambda_str + ")", ytitle="d#Delta/d" + ang.lambda_str)

            gr = construct_deltas_graph(deltas)
            gr.SetName(source.get("label", ""))
            if 'style' in source and 'line_width' not in source['style']:
                source['style']['line_width'] = 2
            c = Contribution(gr, label=source.get("label", ""), marker_style=0, **source.get("style", {}))
            graph_contribs.append(c)

        do_deltas_plot(graph_contribs, 
                       "%s/angularities_pt%dto%d%s.%s" % (plot_dir, start_val, end_val, output_append, ofmt),
                       bin_labels=bin_labels, 
                       title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val), 
                       xtitle="Angularity: (#kappa, #beta)")

