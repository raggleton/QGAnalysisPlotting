#!/usr/bin/env python


"""Compare result of unfolding"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse
from array import array
import numpy as np
import math
import distutils
from distutils import util
from itertools import product
from pprint import pprint

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from my_unfolder import MyUnfolder, HistBinChopper, unpickle_region
from my_unfolder_plotter import MyUnfolderPlotter
from unfolding_config import setup_regions
from do_unfolding_plots import PLOT_STYLES, Setup

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat(".3f")


# Control plot output format
OUTPUT_FMT = "pdf"

class ComparisonPlotter():
    def __init__(self, bins, setup1, hist_bin_chopper1, unfolder1, setup2, hist_bin_chopper2, unfolder2):
        self.bins = bins
        self.setup1 = setup1
        self.region1 = setup1.region  # just to make life easier
        self.hist_bin_chopper1 = hist_bin_chopper1
        self.unfolder1 = unfolder1

        self.setup2 = setup2
        self.region2 = setup2.region
        self.hist_bin_chopper2 = hist_bin_chopper2
        self.unfolder2 = unfolder2

        self.line_width = 2
        self.plot_styles = PLOT_STYLES
        self.pt_bin_plot_args = dict(
            what="hist",
            xtitle=self.setup1.particle_title,
            has_data=self.setup1.has_data or self.setup2.has_data,
            ylim=[0, None],
            subplot_type='ratio',
            # subplot_limits=(0, 2) if self.setup.has_data else (0.75, 1.25),
            # subplot_limits=(0.5, 2.5),
        )

    @staticmethod
    def _modify_plot(this_plot):
        this_plot.legend.SetX1(0.6)
        this_plot.legend.SetY1(0.7)
        this_plot.legend.SetX2(0.98)
        this_plot.legend.SetY2(0.88)
        this_plot.left_margin = 0.16
        this_plot.y_padding_max_linear = 1.8

    @staticmethod
    def check_entries(entries, message=""):
        """Check that at least 1 Contribution has something in it, that is non-zero"""
        has_entries = [c.obj.GetEntries() > 0 for c in entries]
        if not any(has_entries):
            if message:
                print("Skipping 0 entries (%s)" % (message))
            else:
                print("Skipping 0 entries")
            return False

        max_bin = max([c.obj.GetMaximum() for c in entries])
        min_bin = min([c.obj.GetMinimum() for c in entries])
        if max_bin == min_bin:
            if message:
                print("Skipping min=max hists (%s)" % (message))
            else:
                print("Skipping min=max hists")
            return False

        return True

    def get_pt_bin_title(self, bin_edge_low, bin_edge_high):
        region_str = "{jet_algo1}, {region1_label} region\n"
        if self.region1['label'] != self.region2['label']:
            region_str += "{jet_algo2}, {region2_label} region\n"

        title = ((region_str + 
                  "{bin_edge_low:g} < {pt_str} < {bin_edge_high:g} GeV")
                 .format(
                    jet_algo1=self.setup1.jet_algo,
                    region1_label=self.region1['label'],
                    jet_algo2=self.setup2.jet_algo,
                    region2_label=self.region2['label'],
                    pt_str=self.setup1.pt_var_str,
                    bin_edge_low=bin_edge_low,
                    bin_edge_high=bin_edge_high
                ))
        return title

    def plot_unfolded_normalised_pt_bin_offset(self, bin_offset_2=0):

        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            # print(bin_edge_low, bin_edge_high)
            ind2 = ibin+bin_offset_2
            if ind2 < 0:
                # print("...skipping")
                continue

            hbc1_args = dict(ind=ibin, binning_scheme='generator')
            unfolded1_hist_bin_stat_errors = self.hist_bin_chopper1.get_pt_bin_normed_div_bin_width('unfolded_stat_err', **hbc1_args)

            hbc2_args = dict(ind=ind2, binning_scheme='generator')
            unfolded2_hist_bin_stat_errors = self.hist_bin_chopper2.get_pt_bin_normed_div_bin_width('unfolded_stat_err', **hbc2_args)

            entries = [
                Contribution(unfolded1_hist_bin_stat_errors,
                             label="Data (stat. unc.)\n%s" % self.setup1.label,
                             line_color=self.plot_styles['unfolded_stat_colour'],
                             line_width=self.line_width,
                             line_style=1,
                             marker_color=self.plot_styles['unfolded_stat_colour'],
                             marker_style=cu.Marker.get('circle'),
                             marker_size=0.75),
                Contribution(unfolded2_hist_bin_stat_errors,
                             label="Data (stat. unc.)\n%s" % self.setup2.label,
                             line_color=self.plot_styles['unfolded_unreg_colour'],
                             line_width=self.line_width,
                             line_style=1,
                             marker_color=self.plot_styles['unfolded_unreg_colour'],
                             marker_style=cu.Marker.get('square', filled=False),
                             marker_size=0.75,
                             subplot=unfolded1_hist_bin_stat_errors),
            ]
            if not self.check_entries(entries, "plot_unfolded_normalised_pt_bin_offset %d" % (ibin)):
                return

            plot = Plot(entries,
                        ytitle=self.setup1.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        xlim=qgp.calc_auto_xlim(entries),
                        subplot_limits=(0.8, 1.2),
                        **self.pt_bin_plot_args)
            self._modify_plot(plot)
            plot.subplot_title = "* / %s" % (self.region1['label'])
            plot.plot("NOSTACK E1")
            plot.save("%s/compare_unfolded_%s_bin_%d_divBinWidth.%s" % (self.setup1.output_dir, self.setup1.append, ibin, self.setup1.output_fmt))


def do_comparison_binned_plots_per_region_angle(setup1, setup2):
    """"""
    # Do plots of unfolded hists
    region1 = setup1.region
    unfolder1 = region1['unfolder']
    hist_bin_chopper1 = unfolder1.hist_bin_chopper

    region2 = setup2.region
    unfolder2 = region2['unfolder']
    hist_bin_chopper2 = unfolder2.hist_bin_chopper

    plotter = ComparisonPlotter(bins=unfolder1.pt_bin_edges_gen,
                                setup1=setup1,
                                hist_bin_chopper1=hist_bin_chopper1,
                                unfolder1=unfolder1,
                                setup2=setup2,
                                hist_bin_chopper2=hist_bin_chopper2,
                                unfolder2=unfolder2)

    plotter.plot_unfolded_normalised_pt_bin_offset(bin_offset_2=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sourceDir1",
                        help="Source directory with unfolded results")
    parser.add_argument("label1",
                        help="Label for source 1")
    parser.add_argument("sourceDir2",
                        help="Another source directory with unfolded results")
    parser.add_argument("label2",
                        help="Label for source 2")
    parser.add_argument("--outputDir",
                        required=True,
                        help="Output directory for plots")
    parser.add_argument("--angles",
                        choices=list(qgc.VAR_UNFOLD_DICT['ungroomed'].keys()) + ["all"],
                        nargs='+',
                        help="Lambda angles to unfold, or 'all' for all of them (default).",
                        default=["all"])
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


    if args.doAllRegions:
        for x in ['doDijetCentral', 'doDijetForward', 'doDijetCentralGroomed', 'doDijetForwardGroomed', 'doZPJ', 'doZPJGroomed']:
            setattr(args, x, True)

    regions = setup_regions(do_dijet_central=args.doDijetCentral,
                            do_dijet_forward=args.doDijetForward,
                            do_dijet_central_groomed=args.doDijetCentralGroomed,
                            do_dijet_forward_groomed=args.doDijetForwardGroomed,
                            do_zpj=args.doZPJ,
                            do_zpj_groomed=args.doZPJGroomed,
                            source=args.sourceDir1)
    if len(regions) == 0:
        raise RuntimeError("Need at least 1 region")

    if args.angles[0] == "all":
        angles = qgc.COMMON_VARS
    else:
        angles = [a for a in qgc.COMMON_VARS if a.var in args.angles]

    if len(angles) == 0:
        raise RuntimeError("Need at least 1 angle")

    cu.check_dir_exists_create(args.outputDir)

    num_all_iterations = len(regions) * len(angles)
    counter = 1

    has_data1 = not ('_MC_all' in args.sourceDir1 or '_MC_split' in args.sourceDir1)
    has_data2 = not ('_MC_all' in args.sourceDir2 or '_MC_split' in args.sourceDir2)

    jet_algo1 = "AK4"
    if "ak8puppi" in args.sourceDir1:
        jet_algo1 = "AK8"

    jet_algo2 = "AK4"
    if "ak8puppi" in args.sourceDir2:
        jet_algo2 = "AK8"

    # Iterate through regions & variables
    for region in regions:

        region_dir1 = os.path.join(args.sourceDir1, region['name'])
        if not os.path.isdir(region_dir1):
            print("! Warning ! cannot find region dir", region_dir1, '- skipping region')
            continue

        region_dir2 = os.path.join(args.sourceDir2, region['name'])
        if not os.path.isdir(region_dir2):
            print("! Warning ! cannot find region dir", region_dir2, '- skipping region')
            continue

        for angle in angles:
            append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name
            print("*"*120)
            print("Algo/Region/Var: %s %s %s (%d/%d)" % (jet_algo1, region['name'], angle.var, counter, num_all_iterations))
            print("*"*120)
            counter += 1

            angle_dir1 = "%s/%s" % (region_dir1, angle.var)
            if not os.path.isdir(angle_dir1):
                print("! Warning ! cannot find angle dir", angle_dir1, '- skipping angle', angle.var)
                continue
            angle_dir2 = "%s/%s" % (region_dir2, angle.var)
            if not os.path.isdir(angle_dir2):
                print("! Warning ! cannot find angle dir", angle_dir2, '- skipping angle', angle.var)
                continue

            # Get region dicts from pickle file
            pickle_filename1 = os.path.join(angle_dir1, "unfolding_result.pkl")
            unpickled_region1 = unpickle_region(pickle_filename1)
            if unpickled_region1 is None:
                continue

            pickle_filename2 = os.path.join(angle_dir2, "unfolding_result.pkl")
            unpickled_region2 = unpickle_region(pickle_filename2)
            if unpickled_region2 is None:
                continue

            # check
            if region['name'] != unpickled_region1['name']:
                raise RuntimeError("Mismatch region1 name")
            if region['name'] != unpickled_region2['name']:
                raise RuntimeError("Mismatch region2 name")


            setup1 = Setup(jet_algo=jet_algo1,
                           region=unpickled_region1,
                           angle=angle,
                           output_dir=os.path.join(args.outputDir, region['name'], angle.var),
                           has_data=has_data1)
            setup1.label = args.label1  # add custom label
            setup2 = Setup(jet_algo=jet_algo2,
                           region=unpickled_region2,
                           angle=angle,
                           output_dir=os.path.join(args.outputDir, region['name'], angle.var),
                           has_data=has_data2)
            setup2.label = args.label2

            do_comparison_binned_plots_per_region_angle(setup1, setup2)
