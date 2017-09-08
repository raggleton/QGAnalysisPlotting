"""Common functions for doing QG flavour info plots"""


import ROOT
from MyStyle import My_Style
My_Style.cd()
import bisect
import numpy as np
import os
from array import array

# My stuff
from comparator import Contribution, Plot, grab_obj
import common_utils as cu
from qg_common import *
from qg_general_plots import get_projection_plot

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


def get_flavour_fractions(input_file, dirname, which="", var_prepend=""):
    h2d_flav = grab_obj(input_file, "%s/%sjet_%sflavour_vs_pt" % (dirname, var_prepend, which))

    h2d_flav.Rebin2D(1, 5)
    y_axis = h2d_flav.GetYaxis()
    pt_bins_lower, pt_bins_upper = [], []
    flav_dict = {'d': [], 'u': [], 's': [], 'c': [], 'b': [] ,'t': [], 'g': []}
    for i in range(1, y_axis.GetNbins()-1):

        d_frac = h2d_flav.GetBinContent(2, i)
        u_frac = h2d_flav.GetBinContent(3, i)
        s_frac = h2d_flav.GetBinContent(4, i)
        c_frac = h2d_flav.GetBinContent(5, i)
        b_frac = h2d_flav.GetBinContent(6, i)
        t_frac = h2d_flav.GetBinContent(7, i)
        g_frac = h2d_flav.GetBinContent(22, i)

        total = d_frac+u_frac+s_frac+c_frac+b_frac+t_frac+g_frac

        if total == 0:
            continue

        pt_bins_lower.append(y_axis.GetBinLowEdge(i))
        pt_bins_upper.append(y_axis.GetBinLowEdge(i+1))

        flav_dict['d'].append(d_frac / total)
        flav_dict['u'].append(u_frac / total)
        flav_dict['s'].append(s_frac / total)
        flav_dict['c'].append(c_frac / total)
        flav_dict['b'].append(b_frac / total)
        flav_dict['t'].append(t_frac / total)
        flav_dict['g'].append(g_frac / total)

    x_bins = [0.5 * (x1+x2) for x1,x2 in zip(pt_bins_lower[:], pt_bins_upper[:])]

    return x_bins, flav_dict


def compare_flavour_fractions_vs_pt(input_files, dirnames, labels, flav, output_filename, title="", which="", var_prepend=""):
    """Plot a specified flavour fraction vs pT for several sources.
    Each entry in input_files, dirnames, and labels corresponds to one line"""
    info = [get_flavour_fractions(ifile, sel, which=which, var_prepend=var_prepend) for ifile, sel in zip(input_files, dirnames)]
    contribs = []
    for i, (x_bins, fdict) in enumerate(info):
        if flav in ['u', 'd', 's', 'c', 'b', 't', 'g']:
            gr = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(fdict[flav]))
        else:
            gr = ROOT.TGraph(len(x_bins), np.array(x_bins), 1.-np.array(fdict[flav.replace("1-", '')]))
        c = Contribution(gr, label="%s" % (labels[i]), line_style=i+1)
        contribs.append(c)
    ytitle = "%s flavour fraction" % flav
    p = Plot(contribs, what='graph', xtitle="p_{T}^{jet} [GeV]", ytitle=ytitle, title=title, ylim=(0, 1))
    p.legend.SetX1(0.4)
    p.plot("ALP")
    p.canvas.SetLogx()
    p.save(output_filename)


def do_flavour_fraction_vs_pt(input_file, dirname, output_filename, title="", which="", var_prepend=""):
    """Plot flavour fractions vs PT for one input file & dirname in the ROOT file"""
    x_bins, flav_dict = get_flavour_fractions(input_file, dirname, which, var_prepend)
    # TODO: check if empy arrays
    gr_flav_u = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['u']))
    gr_flav_d = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['d']))
    gr_flav_s = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['s']))
    gr_flav_c = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['c']))
    gr_flav_b = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['b']))
    gr_flav_g = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['g']))
    # gr_flav_uds = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array([u+d+s for u, d, s in zip(flav_dict['u'], flav_dict['d'], flav_dict['s'])]))

    plot_u = Contribution(gr_flav_u, label="u frac", line_color=ROOT.kRed, marker_color=ROOT.kRed)
    plot_d = Contribution(gr_flav_d, label="d frac", line_color=ROOT.kBlue, marker_color=ROOT.kBlue)
    plot_s = Contribution(gr_flav_s, label="s frac", line_color=ROOT.kBlack, marker_color=ROOT.kBlack)
    plot_c = Contribution(gr_flav_c, label="c frac", line_color=ROOT.kGreen-3, marker_color=ROOT.kGreen-3)
    plot_b = Contribution(gr_flav_b, label="b frac", line_color=ROOT.kOrange, marker_color=ROOT.kOrange)
    plot_g = Contribution(gr_flav_g, label="g frac", line_color=ROOT.kViolet, marker_color=ROOT.kViolet)
    # plot_uds = Contribution(gr_flav_uds, label="uds frac", line_color=ROOT.kOrange+2, marker_color=ROOT.kOrange+2)

    p_flav = Plot([plot_u, plot_d, plot_s, plot_g, plot_c, plot_b], what='graph', xtitle="p_{T}^{jet} [GeV]", title=title)
    p_flav.plot("ALP")
    p_flav.save(output_filename)
