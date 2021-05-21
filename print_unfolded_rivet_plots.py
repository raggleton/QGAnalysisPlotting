#!/usr/bin/env python


"""
Do all the unfolding plots: per pT bin, per lambda bin, summary plot
"""


from __future__ import print_function, division

import os
import sys
import argparse
import math
from array import array
import pandas as pd
from copy import copy, deepcopy
import numpy as np
import scipy
from itertools import chain
from collections import defaultdict
import inspect
import warnings

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot

import yoda

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from my_unfolder import MyUnfolder, HistBinChopper, unpickle_region, PtVarPerPtBinning
from my_unfolder_plotter import MyUnfolderPlotter
from unfolding_config import setup_regions_from_argparse
import metric_calculators as metrics
from do_summary_unfolding_plots import SAMPLE_STYLE_DICTS
from extract_rivet_summary_stats import dataframe_yoda_key
import rivet_naming as rn
from do_unfolding_plots import (Setup, PLOT_STYLES, calc_chi2_stats, scale_ematrix_by_bin_widths,
                                get_correlated_mean_err, get_correlated_rms_err,
                                get_uncorrelated_mean_err, get_uncorrelated_rms_err)


My_Style.cd()
pd.set_option('display.max_columns', None)
os.nice(10)

# monkey-patch warning formatter
warnings.formatwarning = cu._formatwarning

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()


class YodaGenPtBinnedPlotter():
    def __init__(self, setup, bins, hist_bin_chopper, unfolder):
        self.plots_dict = defaultdict(list)  # use to cache plot names per function name

        self.line_width = 2
        self.plot_styles = PLOT_STYLES

        self.setup = setup
        self.region = setup.region  # just to make life easier
        self.bins = bins
        self.hist_bin_chopper = hist_bin_chopper

        self.pt_bin_plot_args = dict(
            what="hist",
            xtitle=self.setup.particle_title,
            has_data=self.setup.has_data,
            # ylim=[0, None],
            subplot_type='ratio',
            subplot_title="* / %s" % (self.region['mc_label']),
            # subplot_limits=(0, 2) if self.setup.has_data else (0.75, 1.25),
            subplot_limits=(0.5, 2.5) if self.setup.has_data else (0.75, 1.25),
        )
        self.unfolder = unfolder

        self.rivet_entries = []

    def save_plot(self, plot, filename):
        plot.save(filename)
        func_name = inspect.currentframe().f_back.f_code.co_name
        self.plots_dict[func_name].append(filename)

    @staticmethod
    def check_entries(entries, message=""):
        """Check that at least 1 Contribution has something in it, that is non-zero"""
        has_entries = [c.obj.GetEntries() > 0 for c in entries]
        if not any(has_entries):
            if message:
                print("Skipping 0 entries (%s)" % message)
            else:
                print("Skipping 0 entries")
            return False

        max_bin = max([c.obj.GetMaximum() for c in entries])
        min_bin = min([c.obj.GetMinimum() for c in entries])
        if max_bin == min_bin:
            if message:
                print("Skipping min=max hists (%s)" % message)
            else:
                print("Skipping min=max hists")
            return False

        return True

    def _modify_plot(self, this_plot):
        if self.setup.output_fmt == "gif":
            _make_thumbnail_canvas(this_plot)
        this_plot.legend.SetX1(0.6)
        this_plot.legend.SetY1(0.7)
        this_plot.legend.SetX2(0.98)
        this_plot.legend.SetY2(0.88)
        this_plot.left_margin = 0.16
        this_plot.y_padding_max_linear = 1.8
        this_plot.lumi = self.setup.lumi

    def _modify_plot_paper(self, this_plot):
        self._modify_plot(this_plot)
        this_plot.is_preliminary = self.setup.is_preliminary

    def get_pt_bin_title(self, bin_edge_low, bin_edge_high):
        title = (("{jet_algo}\n"
                  "{region_label} region\n"
                  "{bin_edge_low:g} < {pt_str} < {bin_edge_high:g} GeV")
                 .format(
                    jet_algo=self.setup.jet_algo,
                    region_label=self.region['label'],
                    pt_str=self.setup.pt_var_str,
                    bin_edge_low=bin_edge_low,
                    bin_edge_high=bin_edge_high
                ))
        return title

    def add_rivet_entry(self, yoda_dict, style_dict):
        self.rivet_entries.append({'yoda_dict': yoda_dict,
                                   'style_dict': style_dict})

    def plot_unfolded_with_yoda_normalised(self, do_chi2=False, do_zoomed=True):
        data_total_errors_style = dict(label="Data (total unc.)",
                                       line_color=self.plot_styles['unfolded_total_colour'],
                                       line_width=self.line_width,
                                       line_style=1,
                                       marker_color=self.plot_styles['unfolded_total_colour'],
                                       marker_style=cu.Marker.get('circle'),
                                       marker_size=self.plot_styles['unfolded_marker_size'],
                                       leg_draw_opt="LEP")
        data_stat_errors_style = dict(label="Data (stat. unc.)",
                                      line_color=self.plot_styles['unfolded_stat_colour'],
                                      line_width=self.line_width,
                                      line_style=1,
                                      marker_color=self.plot_styles['unfolded_stat_colour'],
                                      marker_style=cu.Marker.get('circle'),
                                      marker_size=0.0001,
                                      leg_draw_opt="LEP")  # you need a non-0 marker to get the horizontal bars at the end of errors

        mc_style = dict(label=self.region['mc_label'],
                        line_color=self.plot_styles['gen_colour'],
                        line_width=self.line_width,
                        marker_color=self.plot_styles['gen_colour'],
                        marker_size=self.plot_styles['gen_marker_size'],
                        marker_style=self.plot_styles['gen_marker'],
                        leg_draw_opt="LEP" if self.plot_styles['gen_marker_size'] > 0 else "LE")

        rivet_path, rivet_region, rivet_radius, rivet_lambda, rivet_pt_bins = get_matching_rivet_setup(self.setup)

        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            mc_gen_hist_bin = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', **hbc_args)
            unfolded_hist_bin_stat_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', **hbc_args)
            unfolded_hist_bin_total_errors = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', **hbc_args)

            # Get RIVET hists, which are absolute counts, so need normalising
            rivet_hist_name = '/%s/%s' % (rivet_path, rn.get_plot_name(rivet_radius, rivet_region, rivet_lambda, rivet_pt_bins[ibin]))
            rivet_hists = [qgp.normalise_hist_divide_bin_width(yoda.root.to_root(ent['yoda_dict'][rivet_hist_name])) for ent in self.rivet_entries]

            # Create copy of data to go on top of stat unc,
            # but remove vertical error bar so we can see the stat unc
            # Note that you CAN'T set it to 0, otherwise vertical lines connecting
            # bins start being drawn. Instead set it to some super small value.
            unfolded_hist_bin_total_errors_marker_noerror = unfolded_hist_bin_total_errors.Clone()  # clone to avoid restyling the original as well
            for i in range(1, unfolded_hist_bin_total_errors_marker_noerror.GetNbinsX()+1):
                unfolded_hist_bin_total_errors_marker_noerror.SetBinError(i, 1E-100)

            data_entries = [
                Contribution(unfolded_hist_bin_total_errors, **data_total_errors_style),
                Contribution(unfolded_hist_bin_stat_errors, **data_stat_errors_style),
                # do data with black marker to get it on top
                Contribution(unfolded_hist_bin_total_errors_marker_noerror, **data_total_errors_style),
            ]

            # For subplot to ensure only MC errors drawn, not MC+data
            data_no_errors = unfolded_hist_bin_total_errors_marker_noerror.Clone()
            cu.remove_th1_errors(data_no_errors)

            this_mc_style = deepcopy(mc_style)

            rivet_styles = []
            for ind, _ in enumerate(rivet_hists):
                s_dict = self.rivet_entries[ind]['style_dict']
                rivet_styles.append(dict(label=s_dict['label'],
                                         line_color=s_dict['color'],
                                         line_width=self.line_width,
                                         marker_color=s_dict['color'],
                                         marker_size=s_dict.get('marker_size', self.plot_styles['gen_marker_size']),
                                         marker_style=s_dict['marker_style'],
                                         leg_draw_opt="LEP" if self.plot_styles['gen_marker_size'] > 0 else "LE"))

            # Calculate chi2 between data and MCs if desired
            if do_chi2:
                # print("unfolded_alt_truth bin", ibin)
                ematrix = self.hist_bin_chopper.get_pt_bin_normed_div_bin_width(self.unfolder.total_ematrix_name, **hbc_args)
                # stats are chi2, ndof, p
                mc_stats = calc_chi2_stats(unfolded_hist_bin_total_errors, mc_gen_hist_bin, ematrix)
                # print(mc_stats)
                # print(alt_mc_stats)
                nbins = sum([1 for i in range(1, unfolded_hist_bin_total_errors.GetNbinsX()+1)
                             if unfolded_hist_bin_total_errors.GetBinContent(i) != 0])
                # reduced_chi2 = mc_stats[0] / nbins
                # alt_reduced_chi2 = alt_mc_stats[0] / nbins

                n_sig_fig = 2
                chi2_template = "\n#lower[-0.1]{{(#chi^{{2}} / N_{{bins}} = {chi2:g} / {nbins:d})}}"
                this_mc_style['label'] += chi2_template.format(chi2=cu.nsf(mc_stats[0], n_sig_fig), nbins=nbins)

                for ind, h in enumerate(rivet_hists):
                    this_stats = calc_chi2_stats(unfolded_hist_bin_total_errors, h, ematrix)
                    rivet_styles[ind]['label'] += chi2_template.format(chi2=cu.nsf(this_stats[0], n_sig_fig), nbins=nbins)

            mc_entries = [
                Contribution(mc_gen_hist_bin, subplot=data_no_errors, **this_mc_style),
            ]

            for h, s_dict in zip(rivet_hists, rivet_styles):
                mc_entries.append(Contribution(h, subplot=data_no_errors, **s_dict))

            entries = [
                # Draw MC
                *mc_entries,
                # Draw data after to put on top of MC
                *data_entries
            ]

            func_name = cu.get_current_func_name()
            if not self.check_entries(entries, "%s bin %d" % (func_name, ibin)):
                return

            ymin = 0
            if np.any(cu.th1_to_ndarray(unfolded_hist_bin_total_errors)[0] < 0):
                ymin = None  # let it do its thing and auto calc ymin
            max_rel_err = 0.5 if "multiplicity" in self.setup.angle.var.lower() else -1
            plot = Plot(entries,
                        ytitle=self.setup.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        legend=True,
                        xlim=qgp.calc_auto_xlim(entries[2:3], max_rel_err=0.5),  # set x lim to where data is non-0
                        ylim=[ymin, None],
                        **self.pt_bin_plot_args)

            plot.subplot_title = qgc.SIM_DATA_STR
            self._modify_plot_paper(plot)

            # disable adding objects to legend & drawing - we'll do it manually
            plot.do_legend = False
            plot.legend.SetTextSize(0.03)
            plot.legend.SetY1(0.6)
            plot.legend.SetX1(0.57)
            plot.legend.SetX2(0.93)
            if len(entries) > 4:
                # if lots of entries, try auto-expand
                plot.legend.SetY1(0.6 - (0.02*(len(entries)-4)))
            # plot.legend.SetEntrySeparation(0.005)
            subplot_draw_opts = "NOSTACK E1"
            plot.plot("NOSTACK E1", subplot_draw_opts)

            dummy_graphs = qgp.do_fancy_legend(chain(data_entries[:2], mc_entries), plot, use_splitline=False)

            plot.canvas.cd()
            plot.legend.Draw()

            # Create hists for data with error region for ratio
            # Easiest way to get errors right is to do data (with 0 errors)
            # and divide by data (with errors), as if you had MC = data with 0 error
            data_stat_ratio = data_no_errors.Clone()
            data_stat_ratio.Divide(unfolded_hist_bin_stat_errors)
            data_stat_ratio.SetFillStyle(3245)
            data_stat_ratio.SetFillColor(self.plot_styles['unfolded_stat_colour'])
            data_stat_ratio.SetLineWidth(0)
            data_stat_ratio.SetMarkerSize(0)

            data_total_ratio = data_no_errors.Clone()
            data_total_ratio.Divide(unfolded_hist_bin_total_errors)
            data_total_ratio.SetFillStyle(3254)
            data_total_ratio.SetFillColor(self.plot_styles['unfolded_total_colour'])
            data_total_ratio.SetLineWidth(0)
            data_total_ratio.SetMarkerSize(0)

            # now draw the data error shaded area
            # this is a bit hacky - basically draw them on the ratio pad,
            # then redraw the existing hists & line to get them ontop
            # note that we use "same" for all - this is to keep the original axes
            # (we may want to rethink this later?)
            plot.subplot_pad.cd()
            draw_opt = "E2 SAME"
            data_stat_ratio.Draw(draw_opt)
            data_total_ratio.Draw(draw_opt)
            plot.subplot_line.Draw()
            plot.subplot_container.Draw("SAME" + subplot_draw_opts)

            # Add subplot legend
            x_left = 0.25
            y_bottom = 0.75
            width = 0.67
            height = 0.15
            plot.subplot_legend = ROOT.TLegend(x_left, y_bottom, x_left+width, y_bottom+height)
            plot.subplot_legend.AddEntry(data_total_ratio, qgc.DATA_TOTAL_UNC_STR, "F")
            plot.subplot_legend.AddEntry(data_stat_ratio, qgc.DATA_STAT_UNC_STR, "F")
            plot.subplot_legend.SetTextSize(0.085)
            plot.subplot_legend.SetFillStyle(0)
            plot.subplot_legend.SetNColumns(2)
            plot.subplot_legend.Draw()
            plot.canvas.cd()

            stp = self.setup
            fname = f"unfolded_{stp.append}_rivet_bin_{ibin:d}_divBinWidth{stp.paper_str}.{stp.output_fmt}"
            self.save_plot(plot, os.path.join(stp.output_dir, fname))

            # Do version with small x values only
            if do_zoomed:
                if self.setup.angle.var in ["jet_thrust_charged", "jet_width_charged", "jet_thrust", "jet_width"]:
                    # plot.ylim = (1E-5)
                    plot.y_padding_max_log = 50
                    plot.y_padding_min_log = 0.5
                    plot.ylim = None
                    plot.set_logy(do_exponent=False, do_more_labels=False)
                    fname = f"unfolded_{stp.append}_alt_truth_bin_{ibin:d}_divBinWidth_logY.{stp.output_fmt}"
                    self.save_plot(plot, os.path.join(stp.output_dir, fname))

                if self.setup.angle.var in ["jet_LHA_charged", "jet_thrust_charged", "jet_width_charged", "jet_thrust", "jet_width"]:
                    bin_edges = cu.get_bin_edges(mc_gen_hist_bin, 'x')
                    # get the bin edge thats smallest between 0.2, and 5th bin
                    bin_lt_lim = [x for x in bin_edges if x < 0.2][-1]
                    upper_bin = min(bin_edges[5], bin_lt_lim)
                    plot2 = Plot(entries,
                                 ytitle=self.setup.pt_bin_normalised_differential_label,
                                 title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                                 xlim=(0, upper_bin),
                                 **self.pt_bin_plot_args)
                    self._modify_plot(plot2)
                    plot2.subplot_title = "* / Generator"
                    plot2.plot("NOSTACK E1")
                    # plot2.set_logx(do_exponent=False)
                    fname = f"unfolded_{stp.append}_rivet_bin_{ibin:d}_divBinWidth_lowX.{stp.output_fmt}"
                    self.save_plot(plot2, os.path.join(stp.output_dir, fname))


def do_binned_plots_per_region_angle(setup, rivet_files, rivet_labels):
    """Do individual binned plots"""
    region = setup.region

    unfolder = region['unfolder']
    hbc = unfolder.hist_bin_chopper

    # Iterate through pt bins - gen binning
    print("Doing YodaGenPtBinnedPlotter...")
    pt_bins = unfolder.binning_handler.get_pt_bins(binning_scheme='generator', is_signal_region=True)
    gen_pt_binned_plotter = YodaGenPtBinnedPlotter(setup=setup,
                                                   bins=pt_bins,
                                                   hist_bin_chopper=hbc,
                                                   unfolder=unfolder)

    # Add RIVET & style dicts
    for rf, rl in zip(rivet_files, rivet_labels):
        yoda_dict = yoda.read(rf)
        gen_pt_binned_plotter.add_rivet_entry(yoda_dict, style_dict=SAMPLE_STYLE_DICTS[dataframe_yoda_key(rl)])

    gen_pt_binned_plotter.plot_unfolded_with_yoda_normalised(do_chi2=True, do_zoomed=True)

    return hbc


def get_matching_rivet_setup(setup):
    is_dijet = 'Dijet' in setup.region['name']
    path = rn.DIJET_PATH if is_dijet else rn.ZPJ_PATH
    regions = rn.DIJET_REGIONS if is_dijet else rn.ZPJ_REGIONS
    rivet_region = next(r for r in regions if r.name == setup.region['name'])
    rivet_radius = next(r for r in rn.JET_RADII if r.name == setup.jet_algo)
    rivet_lambda = next(r for r in rn.LAMBDA_VARS if r.hist_name == setup.angle.var)
    pt_bins = rn.PT_BINS_DIJET if is_dijet else rn.PT_BINS_ZPJ
    return path, rivet_region, rivet_radius, rivet_lambda, pt_bins


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",
                        help="Source directory (should be the one made by unfolding.py)")
    parser.add_argument("--angles",
                        choices=list(qgc.VAR_UNFOLD_DICT['ungroomed'].keys()) + ["all"],
                        nargs='+',
                        help="Lambda angles to unfold, or 'all' for all of them (default).",
                        default=["all"])

    parser.add_argument("--yodaInputDijet",
                        action='append',
                        default=[],
                        help='Yoda input file (from dijet plugin)')
    parser.add_argument("--yodaInputZPJ",
                        action='append',
                        default=[],
                        help='Yoda input file (from Z+Jet plugin)')
    parser.add_argument("--yodaLabel",
                        action='append',
                        default=[],
                        help='Yoda input file label')

    parser.add_argument("--paper",
                        action='store_true',
                        help='Don\'t add "Preliminary" to plots')

    region_group = parser.add_argument_group('Region selection')
    region_group.add_argument("--doAllRegions",
                              action='store_true',
                              help='Do unfolding for all regions (dijet, Z+J, groomed, ungroomed)')
    region_group.add_argument("--doDijetCentral",
                              action='store_true',
                              help='Do unfolding for dijet (central) jets')
    region_group.add_argument("--doDijetForward",
                              action='store_true',
                              help='Do unfolding for dijet (forward) jets')
    region_group.add_argument("--doDijetCentralGroomed",
                              action='store_true',
                              help='Do unfolding for groomed dijet (central) jets')
    region_group.add_argument("--doDijetForwardGroomed",
                              action='store_true',
                              help='Do unfolding for groomed dijet (forward) jets')
    region_group.add_argument("--doZPJ",
                              action='store_true',
                              help='Do unfolding for Z+jet jets')
    region_group.add_argument("--doZPJGroomed",
                              action='store_true',
                              help='Do unfolding for groomed Z+jet jets')

    args = parser.parse_args()

    if (len(args.yodaInputDijet) != len(args.yodaLabel)
        and len(args.yodaInputZPJ) != len(args.yodaLabel)):
        raise RuntimeError("Number of --yodaInputDijet/yodaInputZPJ must match number of --yodaLabel")

    if len(args.yodaInputDijet) == 0:
        raise RuntimeError("Need at least one YODA input")

    if args.doAllRegions:
        for x in ['doDijetCentral', 'doDijetForward', 'doDijetCentralGroomed', 'doDijetForwardGroomed', 'doZPJ', 'doZPJGroomed']:
            setattr(args, x, True)

    regions = setup_regions_from_argparse(args)
    if len(regions) == 0:
        raise RuntimeError("Need at least 1 region")

    if args.angles[0] == "all":
        angles = qgc.COMMON_VARS
    else:
        angles = [a for a in qgc.COMMON_VARS if a.var in args.angles]

    if len(angles) == 0:
        raise RuntimeError("Need at least 1 angle")

    jet_algo = "AK4"
    if "ak8puppi" in args.source:
        jet_algo = "AK8"

    has_data = not ('_MC_all' in args.source or '_MC_split' in args.source)

    num_all_iterations = len(regions) * len(angles)
    counter = 1

    # Iterate through regions & variables
    for region in regions:

        region_dir = os.path.join(args.source, region['name'])
        if not os.path.isdir(region_dir):
            print("! Warning ! cannot find region dir", region_dir, '- skipping region')
            continue

        for angle in angles:
            append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name
            print("*"*120)
            print("Algo/Region/Var: %s %s %s (%d/%d)" % (jet_algo, region['name'], angle.var, counter, num_all_iterations))
            print("*"*120)
            counter += 1

            angle_output_dir = "%s/%s/%s" % (args.source, region['name'], angle.var)
            if not os.path.isdir(angle_output_dir):
                print("! Warning ! cannot find angle dir", angle_output_dir, '- skipping angle', angle.var)
                continue

            # Get region dict from pickle file
            pickle_filename = os.path.join(angle_output_dir, "unfolding_result.pkl")
            unpickled_region = unpickle_region(pickle_filename)
            if unpickled_region is None:
                continue

            # check
            if region['name'] != unpickled_region['name']:
                raise RuntimeError("Mismatch region name")

            # use region dict from unpickling
            # don't use update(), mega slow
            this_region = unpickled_region
            this_region['label'] = region['label']

            # MAKE ALL THE PLOTS
            # ------------------------------------------------------------------
            setup = Setup(jet_algo=jet_algo,
                          region=this_region,
                          angle=angle,
                          output_dir=angle_output_dir,
                          has_data=has_data,
                          is_preliminary=not args.paper)

            hist_bin_chopper = this_region['unfolder'].hist_bin_chopper

            rivet_files = args.yodaInputZPJ if 'ZPlusJets' in region['name'] else args.yodaInputDijet

            print("...........................................................")
            print(" Doing lots of plots per pt/lambda bin...")
            print("...........................................................")
            hist_bin_chopper = do_binned_plots_per_region_angle(setup,
                                                                rivet_files=rivet_files,
                                                                rivet_labels=args.yodaLabel)

            # cleanup object, as it uses loads of memory
            if num_all_iterations > 1:
                print("...tidying up...")
                del hist_bin_chopper
                del unpickled_region
                del setup


if __name__ == "__main__":
    main()
