#!/usr/bin/env python

"""Print main QG plots, comparing systematics samples

"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
import argparse
from array import array
from copy import deepcopy

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


def do_plots(nominal_dir, plot_dir, neutral_hadron_shift_up_dir=None, neutral_hadron_shift_down_dir=None):
    # QG variable plots
    pt_bins = qgc.PT_BINS
    # pt_bins = qgc.THEORY_PT_BINS
    var_list = qgc.COMMON_VARS
    var_prepend = ""

    dj_dirname = "Dijet_QG_tighter"

    rebin_dict = {
        "jet_puppiMultiplicity": [0.0, 6.0, 9.0, 12.0, 18.0, 150.0],
        'jet_multiplicity': [0.0, 6.0, 9.0, 12.0, 18.0, 150.0],
        'jet_pTD': [0.0, 0.1, 0.13, 0.17, 0.23, 0.33, 0.51, 0.85, 1.0],
        'jet_LHA': [0.0, 0.29, 0.37, 0.44, 0.5, 0.56, 0.62, 0.68, 0.75, 1.0],
        'jet_width': [0.0, 0.12, 0.18, 0.24, 0.3, 0.36, 0.43, 0.51, 1.0],
        'jet_thrust': [0.0, 0.04, 0.08, 0.12, 0.17, 0.24, 0.33, 1.0],
    }

    for ang in var_list[:]:

        v = "%s%s_vs_pt" % (var_prepend, ang.var)
        
        for pt_ind, (start_val, end_val) in enumerate(pt_bins[:]):
            dijet_entries = []
            # Get all plots
            lw = 2

            ####################
            # DIJET REGION
            ####################

            # NOMINAL MC
            h2d_qcd_mc = grab_obj(os.path.join(nominal_dir, qgc.QCD_PYTHIA_ONLY_FILENAME), "%s/%s" % (dj_dirname, v))
            qcd_kwargs_mc = dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                                 marker_color=qgc.QCD_COLOUR, marker_style=qgc.QCD_MARKER, marker_size=0,
                                 label=qgc.QCD_Dijet_LABEL + " [nominal]")
            nominal_hist = qgp.get_projection_plot(h2d_qcd_mc, start_val, end_val)
            dijet_entries.append((nominal_hist, qcd_kwargs_mc))

            # NEUTRAL SHIFT UP
            col = qgc.QCD_COLOURS[2]
            h2d_qcd_mc2 = grab_obj(os.path.join(neutral_hadron_shift_up_dir, qgc.QCD_PYTHIA_ONLY_FILENAME), "%s/%s" % (dj_dirname, v))
            qcd_kwargs_mc2 = dict(line_color=col, line_width=lw, fill_color=col,
                                 marker_color=col, marker_style=qgc.QCD_MARKER, marker_size=0,
                                 label=qgc.QCD_Dijet_LABEL + " [Neutral Hadron Shift Up]", subplot=nominal_hist)
            dijet_entries.append((qgp.get_projection_plot(h2d_qcd_mc2, start_val, end_val), qcd_kwargs_mc2))
            
            # NEUTRAL SHIFT DOWN
            col2 = qgc.QCD_COLOURS[3]
            h2d_qcd_mc3 = grab_obj(os.path.join(neutral_hadron_shift_down_dir, qgc.QCD_PYTHIA_ONLY_FILENAME), "%s/%s" % (dj_dirname, v))
            qcd_kwargs_mc3 = dict(line_color=col2, line_width=lw, fill_color=col2,
                                 marker_color=col2, marker_style=qgc.QCD_MARKER, marker_size=0,
                                 label=qgc.QCD_Dijet_LABEL + " [Neutral Hadron Shift Down]", subplot=nominal_hist)
            # h2d_qcd_mc3.Scale(TOTAL_LUMI)
            dijet_entries.append((qgp.get_projection_plot(h2d_qcd_mc3, start_val, end_val), qcd_kwargs_mc3))

            rebin = 2
            v_lower = v.lower()
            if "multiplicity" in v_lower:
                rebin = 2
            elif "flavour" in v_lower or "thrust" in v_lower or 'ptd' in v_lower:
                rebin = 1
            elif "ptd" in v_lower:
                rebin = 5

            xlim = None
            if "width" in v_lower or "ptd" in v_lower:
                xlim = [0, 1]
            elif"thrust" in v_lower:
                xlim = [0, 0.5]
            elif "multiplicity" in v_lower and "ak4" in nominal_dir.lower():
                xlim = [0, 100]
                xlim = [0, 80]
            if xlim:
                xlim[1] = max(xlim[1], nominal_hist.GetBinLowEdge(nominal_hist.GetNbinsX()+1))

            ylim = None
            if "flavour" in v_lower:
                ylim = (0, 1)
            elif "lha" in v_lower:
                ylim = (0, 5)
                ylim = None

            # plot_dir = os.path.join(root_dir, "plots_dy_vs_qcd_mc_vs_data")
            radius, pus = cu.get_jet_config_from_dirname(nominal_dir)
            jet_str = "AK%s PF %s" % (radius.upper(), pus.upper())
            subplot_title = "Variation / nominal"
            subplot_limits = (0.9, 1.1)

            dijet_entries_rebin = []
            this_rebins = rebin_dict[ang.var]
            rebin_hist_norminal = None
            for hist, kwargs in dijet_entries:
                rebin_hist = hist.Rebin(len(this_rebins)-1, hist.GetName()+"Rebin", array('d', this_rebins))
                if not rebin_hist_norminal:
                    rebin_hist_norminal = rebin_hist
                new_kwargs = deepcopy(kwargs)
                if 'subplot' in new_kwargs:
                    new_kwargs['subplot'] = rebin_hist_norminal
                dijet_entries_rebin.append((rebin_hist, new_kwargs))

            qgp.do_comparison_plot(dijet_entries, 
                                   "%s/ptBinned/%s_pt%dto%d_dijet.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                   rebin=rebin,
                                   title="%d < p_{T}^{jet} < %d GeV\n%s" % (start_val, end_val, jet_str),
                                   xtitle=ang.name + " (" + ang.lambda_str + ")",
                                   xlim=xlim, ylim=ylim,
                                   subplot_type='ratio',
                                   subplot_title=subplot_title,
                                   subplot_limits=subplot_limits,
                                   has_data=False)
            
            qgp.do_comparison_plot(dijet_entries_rebin, 
                                   "%s/ptBinned/%s_pt%dto%d_dijet_rebin.%s" % (plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                   title="%d < p_{T}^{jet} < %d GeV\n%s" % (start_val, end_val, jet_str),
                                   xtitle=ang.name + " (" + ang.lambda_str + ")",
                                   xlim=xlim, ylim=ylim,
                                   subplot_type='ratio',
                                   subplot_title=subplot_title,
                                   subplot_limits=subplot_limits,
                                   has_data=False)

       


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--nominal", 
                        help="Directory name for nominal files")
    parser.add_argument("--neutralHadronShiftUp",
                        help="Directory name for neutral hadron scale shift up files",
                        default=None)
    parser.add_argument("--neutralHadronShiftDown",
                        help="Directory name for neutral hadron scale shift down files",
                        default=None)
    parser.add_argument("--outputDir", help="Directory for output file")
    args = parser.parse_args()
    do_plots(args.nominal, args.outputDir, args.neutralHadronShiftUp, args.neutralHadronShiftDown)

