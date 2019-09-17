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


def get_flavour_fractions(input_file, dirname, pt_bins, flav_source="", var_prepend="", which_jet="both"):
    """Get dict of flav : [fraction for specified pt bins] for a given input file & directory.

    flav_source : str
        Which flavour to use, either "" or "genParton_"
    which_jet : "both", 1, 2
        Which jet to look at - both, or one of the two in particular
    """
    jet_str = "" if which_jet == "both" else str(which_jet)
    h2d_flav = grab_obj(input_file, "%s/%sjet%s_%sflavour_vs_pt" % (dirname, var_prepend, jet_str, flav_source))

    flav_dict = {'d': [], 'u': [], 's': [], 'c': [], 'b': [] ,'t': [], 'g': [], 'unknown': []}
    for (pt_min, pt_max) in pt_bins:
        h_flav = get_projection_plot(h2d_flav, pt_min, pt_max)

        unknown_frac = h_flav.GetBinContent(1)
        d_frac = h_flav.GetBinContent(2)
        u_frac = h_flav.GetBinContent(3)
        s_frac = h_flav.GetBinContent(4)
        c_frac = h_flav.GetBinContent(5)
        b_frac = h_flav.GetBinContent(6)
        t_frac = h_flav.GetBinContent(7)
        g_frac = h_flav.GetBinContent(22)

        total = d_frac+u_frac+s_frac+c_frac+b_frac+t_frac+g_frac+unknown_frac

        if total == 0:
            total = 1  # so no Inf

        flav_dict['unknown'].append(unknown_frac / total)
        flav_dict['d'].append(d_frac / total)
        flav_dict['u'].append(u_frac / total)
        flav_dict['s'].append(s_frac / total)
        flav_dict['c'].append(c_frac / total)
        flav_dict['b'].append(b_frac / total)
        flav_dict['t'].append(t_frac / total)
        flav_dict['g'].append(g_frac / total)

    return flav_dict


def compare_flavour_fractions_vs_pt(input_files, dirnames, pt_bins, labels, flav, output_filename, title="", flav_source="", var_prepend="", which_jet="both", xtitle="p_{T}^{jet} [GeV]"):
    """Plot a specified flavour fraction vs pT for several sources.
    Each entry in input_files, dirnames, and labels corresponds to one line

    TODO: fix this - bit stupid input format
    """
    bin_centers = [0.5*(x[0]+x[1]) for x in pt_bins]
    bin_widths = [0.5*(x[1]-x[0]) for x in pt_bins]
    info = [get_flavour_fractions(ifile, sel, pt_bins, flav_source=flav_source, var_prepend=var_prepend, which_jet=(which_jet if "Dijet" in sel else "both")) 
            for ifile, sel in zip(input_files, dirnames)]
    contribs = []
    N = len(bin_centers)
    colours = [ROOT.kRed, ROOT.kBlack, ROOT.kBlue, ROOT.kGreen-3]
    for i, fdict in enumerate(info):
        if flav in ['u', 'd', 's', 'c', 'b', 't', 'g']:
            gr = ROOT.TGraphErrors(N, np.array(bin_centers), np.array(fdict[flav]), np.array(bin_widths), np.zeros(N))
        else:
            gr = ROOT.TGraphErrors(N, np.array(bin_centers), 1.-np.array(fdict[flav.replace("1-", '')]), np.array(bin_widths), np.zeros(N))
        c = Contribution(gr, label="%s" % (labels[i]), line_color=colours[i], marker_style=20+i, marker_color=colours[i])
        contribs.append(c)
    ytitle = "%s flavour fraction" % flav
    p = Plot(contribs, what='graph', xtitle=xtitle, ytitle=ytitle, title=title, ylim=(0, 1))
    # p.legend.SetX1(0.4)
    p.plot("ALP")
    p.set_logx()
    p.save(output_filename)


def do_flavour_fraction_vs_pt(input_file, dirname, pt_bins, output_filename, title="", flav_source="", var_prepend="", which_jet="both"):
    """Plot all flavour fractions vs PT for one input file & dirname in the ROOT file"""
    bin_centers = [0.5*(x[0]+x[1]) for x in pt_bins]
    bin_widths = [0.5*(x[1]-x[0]) for x in pt_bins]
    flav_dict = get_flavour_fractions(input_file, dirname, pt_bins, flav_source, var_prepend, which_jet)
    
    # TODO: check if empty arrays
    N = len(bin_centers)
    gr_flav_u = ROOT.TGraphErrors(N, np.array(bin_centers), np.array(flav_dict['u']), np.array(bin_widths), np.zeros(N))
    gr_flav_d = ROOT.TGraphErrors(N, np.array(bin_centers), np.array(flav_dict['d']), np.array(bin_widths), np.zeros(N))
    gr_flav_s = ROOT.TGraphErrors(N, np.array(bin_centers), np.array(flav_dict['s']), np.array(bin_widths), np.zeros(N))
    gr_flav_c = ROOT.TGraphErrors(N, np.array(bin_centers), np.array(flav_dict['c']), np.array(bin_widths), np.zeros(N))
    gr_flav_b = ROOT.TGraphErrors(N, np.array(bin_centers), np.array(flav_dict['b']), np.array(bin_widths), np.zeros(N))
    gr_flav_g = ROOT.TGraphErrors(N, np.array(bin_centers), np.array(flav_dict['g']), np.array(bin_widths), np.zeros(N))
    gr_flav_unknown = ROOT.TGraphErrors(N, np.array(bin_centers), np.array(flav_dict['unknown']), np.array(bin_widths), np.zeros(N))

    plot_u = Contribution(gr_flav_u, label="u", line_color=ROOT.kRed, marker_color=ROOT.kRed, marker_style=20)
    plot_d = Contribution(gr_flav_d, label="d", line_color=ROOT.kBlue, marker_color=ROOT.kBlue, marker_style=21)
    plot_s = Contribution(gr_flav_s, label="s", line_color=ROOT.kBlack, marker_color=ROOT.kBlack, marker_style=22)
    plot_c = Contribution(gr_flav_c, label="c", line_color=ROOT.kGreen-3, marker_color=ROOT.kGreen-3, marker_style=23)
    plot_b = Contribution(gr_flav_b, label="b", line_color=ROOT.kOrange, marker_color=ROOT.kOrange, marker_style=24)
    plot_g = Contribution(gr_flav_g, label="g", line_color=ROOT.kViolet, marker_color=ROOT.kViolet, marker_style=25)
    plot_unknown = Contribution(gr_flav_unknown, label="unknown", line_color=ROOT.kGray+1, marker_color=ROOT.kGray+1, marker_style=26)

    p_flav = Plot([plot_u, plot_d, plot_s, plot_c, plot_b, plot_g, plot_unknown], what='graph', 
                  xtitle="p_{T}^{jet} [GeV]", ytitle="Fraction", title=title, ylim=[0, 1])
    p_flav.plot("ALP")
    p_flav.save(output_filename)
