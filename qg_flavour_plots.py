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


def get_jet_str(thing):
    if "gen" in thing.lower():
        return "genjet"
    else:
        return "jet"


def get_flavour_fractions(input_file, dirname, pt_bins, var_prepend="", which_jet="both"):
    """Get dict of flav : [fraction for specified pt bins] for a given input file & directory.

    which_jet : "both", 1, 2
        Which jet to look at - both, or one of the two in particular
    """
    jet_str = "" if which_jet == "both" else str(which_jet)
    h2d_flav = grab_obj(input_file, "%s/%sjet%s_flavour_vs_pt" % (dirname, var_prepend, jet_str))

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

def get_flavour_efficiencies(input_file, dirname, pt_bins, var_prepend="", which_jet="both"):
    """Get dict of flav : [TEfficiency for specified pt bins] for a given input file & directory.

    which_jet : "both", 1, 2
        Which jet to look at - both, or one of the two in particular
    """
    jet_str = "" if which_jet == "both" else str(which_jet)
    h2d_flav = grab_obj(input_file, "%s/%sjet%s_flavour_vs_pt" % (dirname, var_prepend, jet_str))

    flav_dict = {'d': [], 'u': [], 's': [], 'c': [], 'b': [] ,'t': [], 'g': [], 'unknown': [], 'total': []}

    for (pt_min, pt_max) in pt_bins:
        h_flav = get_projection_plot(h2d_flav, pt_min, pt_max)

        total_err = array('d', [-1.])
        total2 = h_flav.IntegralAndError(1, h_flav.GetNbinsX()+1, total_err)
        total_err = total_err[0]

        unknown_num = h_flav.GetBinContent(1)
        d_num = h_flav.GetBinContent(2)
        u_num = h_flav.GetBinContent(3)
        s_num = h_flav.GetBinContent(4)
        c_num = h_flav.GetBinContent(5)
        b_num = h_flav.GetBinContent(6)
        t_num = h_flav.GetBinContent(7)
        g_num = h_flav.GetBinContent(22)

        total = d_num+u_num+s_num+c_num+b_num+t_num+g_num+unknown_num

        if not cu.same_floats(total, total2):
            raise RuntimeError("totals dont match: %.9g vs %.9g" % (total, total2))

        unknown_err = h_flav.GetBinError(1)
        d_err = h_flav.GetBinError(2)
        u_err = h_flav.GetBinError(3)
        s_err = h_flav.GetBinError(4)
        c_err = h_flav.GetBinError(5)
        b_err = h_flav.GetBinError(6)
        t_err = h_flav.GetBinError(7)
        g_err = h_flav.GetBinError(22)

        # print(pt_min, pt_max, b_num, b_err, total2, total_err)
        # print(pt_min, pt_max, g_num, g_err, total2, total_err)

        flav_dict['unknown'].append([unknown_num, unknown_err])
        flav_dict['d'].append([d_num, d_err])
        flav_dict['u'].append([u_num, u_err])
        flav_dict['s'].append([s_num, s_err])
        flav_dict['c'].append([c_num, c_err])
        flav_dict['b'].append([b_num, b_err])
        flav_dict['t'].append([t_num, t_err])
        flav_dict['g'].append([g_num, g_err])
        flav_dict['total'].append([total, total_err])


    flav_str_dict = {
        'u': 'Up quark',
        'd': 'Down quark',
        'c': 'Charm quark',
        's': 'Strange quark',
        'b': 'Bottom quark',
        't': 'Top quark',
        'g': 'Gluon',
        '1-g': 'Non-gluon',
    }
    flav_eff_dict = {
        'd': None,
        'u': None,
        's': None,
        'c': None,
        'b': None,
        't': None,
        'g': None,
        'unknown': None,
    }

    bin_edges = [p[0] for p in pt_bins]
    bin_edges += pt_bins[-1][1:]
    bin_edges = array('d', bin_edges)

    h_total = ROOT.TH1D("h_total", ";p_{T}^{%s} [GeV];N" % (get_jet_str(var_prepend)), len(pt_bins), bin_edges)
    for ind, val in enumerate(flav_dict['total'], 1):
        h_total.SetBinContent(ind, val[0])
        h_total.SetBinError(ind, val[1])

    for flav_key in [k for k in flav_dict.keys() if k != 'total']:
        h_flav = ROOT.TH1D("h_%s" % flav_key, ";p_{T}^{%s} [GeV];N" % (get_jet_str(var_prepend)).title(), len(pt_bins), bin_edges)
        for ind, val in enumerate(flav_dict[flav_key], 1):
            h_flav.SetBinContent(ind, val[0])
            h_flav.SetBinError(ind, val[1])

        eff_flav = ROOT.TEfficiency(h_flav, h_total)
        flav_eff_dict[flav_key] = eff_flav

        # gr_flav = ROOT.TGraphAsymmErrors(h_flav, h_total)
        # flav_eff_dict[flav_key] = gr_flav

        # print('bin12 val:', flav_key, val[0], val[1])
        # print('bin12 total:', flav_dict['total'][-1])
        # print('bin12 teff:', flav_key, eff_flav.GetEfficiency(12), eff_flav.GetEfficiencyErrorLow(12), eff_flav.GetEfficiencyErrorUp(12))
        # print('bin12 tgraph:', flav_key, gr_flav.GetErrorYlow(12), gr_flav.GetErrorYhigh(12))

    return flav_eff_dict


def compare_flavour_fractions_vs_pt(input_files, dirnames, pt_bins, labels, flav, output_filename, title="", var_prepend="", which_jet="both", xtitle="p_{T}^{jet} [GeV]"):
    """Plot a specified flavour fraction vs pT for several sources.
    Each entry in input_files, dirnames, and labels corresponds to one line

    TODO: fix this - bit stupid input format
    """
    bin_centers = [0.5*(x[0]+x[1]) for x in pt_bins]
    bin_widths = [0.5*(x[1]-x[0]) for x in pt_bins]
    info = [get_flavour_efficiencies(ifile, sel, pt_bins, var_prepend=var_prepend, which_jet=(which_jet if "Dijet" in sel else "both"))
            for ifile, sel in zip(input_files, dirnames)]
    contribs = []
    N = len(bin_centers)
    colours = [ROOT.kRed, ROOT.kBlack, ROOT.kBlue, ROOT.kGreen-3]

    for i, fdict in enumerate(info):
        if flav in ['u', 'd', 's', 'c', 'b', 't', 'g']:
            obj = fdict[flav].CreateGraph()
        else:
            raise RuntimeError("Robin broke 1-X functionality")
            obj = ROOT.TGraphErrors(N, np.array(bin_centers), 1.-np.array(fdict[flav.replace("1-", '')]), np.array(bin_widths), np.zeros(N))

        c = Contribution(obj,
                         label="%s" % (labels[i]),
                         line_color=colours[i], line_width=1,
                         marker_style=20+i, marker_color=colours[i], marker_size=1,
                         leg_draw_opt="LEP")
        contribs.append(c)

    flav_str_dict = {
        'u': 'Up quark',
        'd': 'Down quark',
        'c': 'Charm quark',
        's': 'Strange quark',
        'b': 'Bottom quark',
        't': 'Top quark',
        'g': 'Gluon',
        '1-g': 'Non-gluon',
    }
    flav_str = flav_str_dict[flav]
    ytitle = "Fraction of %s flavour %ss" % (flav_str.lower(), get_jet_str(var_prepend))
    p = Plot(contribs,
             what='graph',
             xtitle=xtitle,
             ytitle=ytitle,
             title=title,
             xlim=(pt_bins[0][0], pt_bins[-1][1]),
             ylim=(0, 1),
             has_data=False)
    p.default_canvas_size = (600, 600)
    p.plot("ALP")
    p.main_pad.SetBottomMargin(0.16)
    p.get_modifier().GetXaxis().SetTitleOffset(1.4)
    p.get_modifier().GetXaxis().SetTitleSize(.045)
    p.legend.SetX1(0.56)
    p.legend.SetY1(0.6)
    p.legend.SetY2(0.87)
    p.set_logx(do_more_labels=False)
    p.save(output_filename)


def do_flavour_fraction_vs_pt(input_file, dirname, pt_bins, output_filename, title="", var_prepend="", which_jet="both"):
    """Plot all flavour fractions vs PT for one input file & dirname in the ROOT file"""
    info = get_flavour_efficiencies(input_file, dirname, pt_bins, var_prepend=var_prepend, which_jet=(which_jet if "Dijet" in dirname else "both"))

    plot_u = Contribution(info['u'].CreateGraph(), label="Up", line_color=ROOT.kRed, marker_color=ROOT.kRed, marker_style=20, leg_draw_opt="LEP")
    plot_d = Contribution(info['d'].CreateGraph(), label="Down", line_color=ROOT.kBlue, marker_color=ROOT.kBlue, marker_style=21, leg_draw_opt="LEP")
    plot_s = Contribution(info['s'].CreateGraph(), label="Strange", line_color=ROOT.kBlack, marker_color=ROOT.kBlack, marker_style=22, leg_draw_opt="LEP")
    plot_c = Contribution(info['c'].CreateGraph(), label="Charm", line_color=ROOT.kGreen-3, marker_color=ROOT.kGreen-3, marker_style=23, leg_draw_opt="LEP")
    plot_b = Contribution(info['b'].CreateGraph(), label="Bottom", line_color=ROOT.kOrange-3, marker_color=ROOT.kOrange-3, marker_style=33, leg_draw_opt="LEP")
    plot_g = Contribution(info['g'].CreateGraph(), label="Gluon", line_color=ROOT.kViolet, marker_color=ROOT.kViolet, marker_style=29, leg_draw_opt="LEP")
    plot_unknown = Contribution(info['unknown'].CreateGraph(), label="Unknown", line_color=ROOT.kGray+1, marker_color=ROOT.kGray+1, marker_style=26, leg_draw_opt="LEP")

    p_flav = Plot([plot_d, plot_u, plot_s, plot_c, plot_b, plot_g, plot_unknown],
                  what='graph',
                  xtitle="p_{T}^{%s} [GeV]" % get_jet_str(var_prepend).title(),
                  ytitle="Fraction",
                  title=title,
                  xlim=(pt_bins[0][0], pt_bins[-1][1]),
                  ylim=[0, 1],
                  has_data=False)
    p_flav.default_canvas_size = (600, 600)
    p_flav.plot("ALP")
    p_flav.main_pad.SetBottomMargin(0.16)
    p_flav.get_modifier().GetXaxis().SetTitleOffset(1.4)
    p_flav.get_modifier().GetXaxis().SetTitleSize(.045)
    p_flav.set_logx(do_more_labels=False)
    p_flav.legend.SetX1(0.56)
    p_flav.legend.SetY1(0.6)
    p_flav.legend.SetY2(0.87)
    p_flav.save(output_filename)
