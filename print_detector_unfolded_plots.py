#!/usr/bin/env python


"""
Plot detector & unfolding distributions together for data only

Does dijet/zpj separately, and together

"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse
from copy import copy, deepcopy
from itertools import chain
import numpy as np
from pprint import pprint
from array import array

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from my_unfolder import MyUnfolder, HistBinChopper, unpickle_region
from my_unfolder_plotter import MyUnfolderPlotter
from unfolding_config import get_dijet_config, get_zpj_config
from do_unfolding_plots import PLOT_COLOURS
import metric_calculators as metrics

import rootpy

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()


# Define own linestyle for smaller plots
# linestyle 2 (dashed) has too big dashes
ROOT.gStyle.SetLineStyleString(22, "8 6")
ROOT.gStyle.SetLineStyleString(23, "5 10")


def scale_ematrix_by_bin_widths(ematrix, widths):
    this_widths = widths.reshape(len(widths), 1)
    return ematrix * this_widths * this_widths.T


class DijetZPJGenPtBinnedPlotter(object):
    """Class to do plots of dijet and/or ZPJ with generator pt binning"""

    def __init__(self, jet_algo, angle, bins, dijet_region, zpj_region, output_dir, is_groomed):
        self.jet_algo = jet_algo
        self.has_data = has_data
        self.angle = angle
        self.dijet_region = dijet_region
        if self.dijet_region is not None:
            self.dijet_hbc = self.dijet_region['unfolder'].hist_bin_chopper
            self.dijet_hbc.add_obj("input_hist_gen_binning_bg_subtracted", dijet_region['unfolder'].input_hist_gen_binning_bg_subtracted)
        self.zpj_region = zpj_region
        if self.zpj_region is not None:
            self.zpj_hbc = self.zpj_region['unfolder'].hist_bin_chopper
            self.zpj_hbc.add_obj("input_hist_gen_binning_bg_subtracted", zpj_region['unfolder'].input_hist_gen_binning_bg_subtracted)
        self.bins = bins
        self.is_groomed = is_groomed

        # for plot axis titles
        self.pt_var_str = "p_{T}^{jet}"
        self.pt_str = self.pt_var_str + " [GeV]"

        angle_prepend = "Groomed " if self.is_groomed else ""
        this_angle_name = angle.name
        if (angle_prepend != ""
            and 'LHA' not in this_angle_name
            and "_{T}" not in this_angle_name
            and "PUPPI" not in this_angle_name):
            # lower case if Groomed..., but be careful of e.g. pTD, LHA
            this_angle_name = this_angle_name[0].lower() + this_angle_name[1:]

        self.angle_str = "{prepend}{name} ({lambda_str})".format(prepend=angle_prepend,
                                                                 name=this_angle_name,
                                                                 lambda_str=angle.lambda_str)
        # self.particle_title = "Particle-level " + self.angle_str
        self.particle_title = self.angle_str

        self.detector_title = "Detector-level {prepend}{name} ({lambda_str})".format(prepend=angle_prepend.lower(),  # lower as no longer first
                                                                                     name=this_angle_name,
                                                                                     lambda_str=angle.lambda_str)
        self.pt_bin_normalised_differential_label = "#frac{1}{d#sigma/dp_{T}} #frac{d^{2}#sigma}{dp_{T} d%s}" % (angle.lambda_str)
        self.pt_bin_detector_normalised_differential_label = "#frac{1}{dN/dp_{T}} #frac{d^{2}N}{dp_{T} d%s}" % (angle.lambda_str)
        self.pt_bin_normalised_differential_times_width_label = "#frac{1}{d#sigma/dp_{T}} #frac{d^{2}#sigma}{dp_{T} d%s}" % (angle.lambda_str)
        self.pt_bin_unnormalised_differential_label = "#frac{1}{dN/dp_{T}} #frac{d^{2}N}{dp_{T} d%s}" % (angle.lambda_str)  # FIXME
        self.pt_bin_unnormalised_differential_label = "N"
        self.lambda_bin_normalised_differential_label = "#frac{1}{d#sigma/d%s} #frac{d^{2}#sigma}{dp_{T} d%s}" % (angle.lambda_str, angle.lambda_str)
        self.lambda_bin_unnormalised_differential_label = "#frac{1}{dN/d%s} #frac{d^{2}N}{dp_{T} d%s}" % (angle.lambda_str, angle.lambda_str)  # FIXME
        self.lambda_bin_unnormalised_differential_label = "N"
        self.output_dir = output_dir
        self.output_fmt = 'pdf'
        angle_prepend = "groomed_" if self.is_groomed else ""
        self.append = "%s%s" % (angle_prepend, angle.var)  # common str to put on filenames, etc

        self.plot_colours = dict(
            dijet_colour=ROOT.kBlack,
            zpj_colour=ROOT.kBlue
        )
        self.line_width = 2
        self.line_style_detector = 22


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
        title = (("{jet_algo}\n"
                  "{bin_edge_low:g} < {pt_str} < {bin_edge_high:g} GeV")
                 .format(
                    jet_algo=self.jet_algo,
                    pt_str=self.pt_var_str,
                    bin_edge_low=bin_edge_low,
                    bin_edge_high=bin_edge_high
                ))
        return title

    def plot_detector_unfolded(self, do_dijet=True, do_zpj=True):
        for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(self.bins[:-1], self.bins[1:])):
            hbc_args = dict(ind=ibin, binning_scheme='generator')
            dijet_entries = []
            # kerning necessary as it puts massive space around #pm
            # but the first #kern doesn't affect the space after #pm as much (?!),
            # so I have to add another one with harder kerning
            mean_template = 'Mean = {:.2f}#kern[-0.2dx]{{ #pm}}#kern[-0.5dx]{{ }}{:.1g}'
            if do_dijet:
                # get detector-level data
                detector_hist = self.dijet_hbc.get_pt_bin_normed_div_bin_width('input_hist_gen_binning_bg_subtracted', **hbc_args)
                # Get mean to add to plot
                contents, errors = cu.th1_to_ndarray(detector_hist)
                widths = cu.get_th1_bin_widths(detector_hist)
                centers = cu.get_th1_bin_centers(detector_hist)
                # need to multiply by widths, since the original hist has bin contents divided by width
                bin_areas = contents*widths
                # print("Detector bin areas:", bin_areas)
                # print("Detector bin areas * centers:", bin_areas * centers)
                areas, centers = metrics.hist_values_to_uarray(bin_areas=bin_areas, bin_centers=centers, bin_errors=errors*widths)
                mean_u = metrics.calc_mean_ucert(areas, centers)
                detector_mean, detector_mean_err = mean_u.nominal_value, mean_u.std_dev

                # get unfolded data
                unfolded_hist = self.dijet_hbc.get_pt_bin_normed_div_bin_width("unfolded", **hbc_args)
                contents, _ = cu.th1_to_ndarray(unfolded_hist)
                contents = contents.reshape(-1)  # necessary for jax to avoid shape discrepancy
                centers = cu.get_th1_bin_centers(unfolded_hist)
                bin_areas = contents*widths
                # print("Unfolded bin areas:", bin_areas)
                # print("Unfolded bin areas * centers:", bin_areas * centers)

                cov_matrix = self.dijet_hbc.get_pt_bin_normed_div_bin_width(self.dijet_region['unfolder'].total_ematrix_name, **hbc_args)
                unfolded_ematrix, _ = cu.th2_to_ndarray(cov_matrix)
                # need to scale ematrix by bin areas
                unfolded_ematrix = scale_ematrix_by_bin_widths(unfolded_ematrix, widths)
                unfolded_mean = float(metrics.calc_mean_jax(bin_areas, centers))
                unfolded_mean_err = float(metrics.calc_mean_cov_matrix_jax(bin_areas, centers, unfolded_ematrix))

                dijet_entries.append(Contribution(detector_hist,
                                            label='Detector-level\n%s' % (mean_template.format(detector_mean, detector_mean_err)),
                                            line_color=self.plot_colours['dijet_colour'],
                                            line_width=self.line_width,
                                            line_style=self.line_style_detector,
                                            marker_color=self.plot_colours['dijet_colour'],
                                            marker_style=cu.Marker.get('circle', filled=False),
                                            marker_size=0.75))
                dijet_entries.append(Contribution(unfolded_hist,
                                            label='Particle-level\n%s' % (mean_template.format(unfolded_mean, unfolded_mean_err)),
                                            line_color=self.plot_colours['dijet_colour'],
                                            line_width=self.line_width,
                                            line_style=1,
                                            marker_color=self.plot_colours['dijet_colour'],
                                            marker_style=cu.Marker.get('circle', filled=True),
                                            marker_size=0.75))
            zpj_entries = []
            if do_zpj:
                # get detector-level data
                detector_hist = self.zpj_hbc.get_pt_bin_normed_div_bin_width('input_hist_gen_binning_bg_subtracted', **hbc_args)
                # Get mean to add to plot
                contents, errors = cu.th1_to_ndarray(detector_hist)
                widths = cu.get_th1_bin_widths(detector_hist)
                centers = cu.get_th1_bin_centers(detector_hist)
                # need to multiply by widths, since the original hist has bin contents divided by width
                bin_areas = contents*widths
                areas, centers = metrics.hist_values_to_uarray(bin_areas=contents*widths, bin_centers=centers, bin_errors=errors*widths)
                mean_u = metrics.calc_mean_ucert(areas, centers)
                detector_mean, detector_mean_err = mean_u.nominal_value, mean_u.std_dev

                # get unfolded data
                unfolded_hist = self.zpj_hbc.get_pt_bin_normed_div_bin_width("unfolded", **hbc_args)
                contents, _ = cu.th1_to_ndarray(unfolded_hist)
                contents = contents.reshape(-1)  # necessary for jax to avoid shape discrepancy
                centers = cu.get_th1_bin_centers(unfolded_hist)
                bin_areas = contents*widths

                cov_matrix = self.zpj_hbc.get_pt_bin_normed_div_bin_width(self.zpj_region['unfolder'].total_ematrix_name, **hbc_args)
                unfolded_ematrix, _ = cu.th2_to_ndarray(cov_matrix)
                # need to scale ematrix by bin areas
                unfolded_ematrix = scale_ematrix_by_bin_widths(unfolded_ematrix, widths)
                unfolded_mean = float(metrics.calc_mean_jax(bin_areas, centers))
                unfolded_mean_err = float(metrics.calc_mean_cov_matrix_jax(bin_areas, centers, unfolded_ematrix))


                zpj_entries.append(Contribution(detector_hist,
                                            label='Detector-level\n%s' % (mean_template.format(detector_mean, detector_mean_err)),
                                            line_color=self.plot_colours['zpj_colour'],
                                            line_width=self.line_width,
                                            line_style=self.line_style_detector,
                                            marker_color=self.plot_colours['zpj_colour'],
                                            marker_style=cu.Marker.get('square', filled=False),
                                            marker_size=0.75))
                zpj_entries.append(Contribution(unfolded_hist,
                                            label='Particle-level\n%s' % (mean_template.format(unfolded_mean, unfolded_mean_err)),
                                            line_color=self.plot_colours['zpj_colour'],
                                            line_width=self.line_width,
                                            line_style=1,
                                            marker_color=self.plot_colours['zpj_colour'],
                                            marker_style=cu.Marker.get('square', filled=True),
                                            marker_size=0.75))

            all_entries = list(chain(dijet_entries, zpj_entries))
            plot = Plot(all_entries,
                        ytitle=self.pt_bin_normalised_differential_label,
                        title=self.get_pt_bin_title(bin_edge_low, bin_edge_high),
                        legend=True,
                        xlim=qgp.calc_auto_xlim(all_entries),  # set x lim to where data is non-0
                        what="hist",
                        xtitle=self.particle_title,
                        has_data=self.has_data,
                        ylim=[0, None],
                        )
            plot.default_canvas_size = (600, 600)
            plot.left_margin = 0.2
            plot.left_title_offset_fudge_factor = 8
            plot.y_padding_max_linear = 1.8

            # disable adding objects to legend & drawing - we'll do it manually
            # since we want proper error bar
            plot.do_legend = False
            plot.splitline_legend = False
            # plot.legend.SetFillColor(ROOT.kRed)
            # plot.legend.SetFillStyle(1001)
            plot.plot("NOSTACK  E1")
            plot.get_modifier().GetYaxis().SetTitleOffset(plot.get_modifier().GetYaxis().GetTitleOffset()*1.5)

            # Create dummy graphs with the same styling to put into the legend
            dummy_gr = ROOT.TGraphErrors(1, array('d', [1]), array('d', [1]), array('d', [1]), array('d', [1]))
            dummies = []  # to stop garbage collection
            label_height = 0.03
            legend_height = 0.17
            legend_x1 = 0.58
            legend_x2 = 0.85
            label_left_offset = 0.01
            label_text_size = 0.037
            label_top = 0.85
            legend_text_size = 0.032
            inter_region_offset = 0.025
            if do_dijet:
                dijet_legend = plot.legend.Clone()
                dijet_legend.SetX1(legend_x1)
                dijet_legend.SetX2(legend_x2)
                dijet_legend.SetY1(label_top-label_height-legend_height)
                dijet_legend.SetY2(label_top-label_height+0.01)
                dijet_legend.SetTextSize(legend_text_size)
                # Add text with region label
                dijet_pt = ROOT.TPaveText(legend_x1-label_left_offset, label_top-label_height, legend_x2-label_left_offset, label_top, "NDC NB")
                dijet_pt.SetFillStyle(0)
                dijet_pt.SetBorderSize(0)
                text = dijet_pt.AddText(self.dijet_region['label'] + " region")
                text.SetTextAlign(11)
                text.SetTextFont(62)
                text.SetTextSize(label_text_size)
                dummies.append(dijet_pt)
                dummies.append(text)
                dijet_pt.Draw()
                for cont in dijet_entries:
                    this_dummy_entry = dummy_gr.Clone()
                    cont.update_obj_styling(this_dummy_entry)
                    if "\n" in cont.label:
                        parts = cont.label.split("\n")
                        dijet_legend.AddEntry(this_dummy_entry, "#splitline{%s}{%s}" % (parts[0], parts[1]), "LEP")
                    else:
                        dijet_legend.AddEntry(this_dummy_entry, cont.label, "LEP")
                    dummies.append(this_dummy_entry)  # to avoid garbage collection
                dijet_legend.Draw()
                dummies.append(dijet_legend)
                # setup for Z+J
                label_top -= label_height+legend_height+inter_region_offset

            if do_zpj:
                zpj_legend = plot.legend.Clone()
                zpj_legend.SetX1(legend_x1)
                zpj_legend.SetX2(legend_x2)
                zpj_legend.SetY1(label_top-label_height-legend_height)
                zpj_legend.SetY2(label_top-label_height+0.01)
                zpj_legend.SetTextSize(legend_text_size)
                # Add text with region label
                zpj_pt = ROOT.TPaveText(legend_x1-label_left_offset, label_top-label_height, legend_x2-label_left_offset, label_top, "NDC NB")
                zpj_pt.SetFillStyle(0)
                zpj_pt.SetBorderSize(0)
                text = zpj_pt.AddText(self.zpj_region['label'] + " region")
                text.SetTextAlign(11)
                text.SetTextFont(62)
                text.SetTextSize(label_text_size)
                dummies.append(zpj_pt)
                dummies.append(text)
                zpj_pt.Draw()
                for cont in zpj_entries:
                    this_dummy_entry = dummy_gr.Clone()
                    cont.update_obj_styling(this_dummy_entry)
                    if "\n" in cont.label:
                        parts = cont.label.split("\n")
                        zpj_legend.AddEntry(this_dummy_entry, "#splitline{%s}{%s}" % (parts[0], parts[1]), "LEP")
                    else:
                        zpj_legend.AddEntry(this_dummy_entry, cont.label, "LEP")
                    dummies.append(this_dummy_entry)
                zpj_legend.Draw()
                dummies.append(zpj_legend)

            parts = ['detector_unfolded',
                     'dijet' if do_dijet else None,
                     'zpj' if do_zpj else None,
                     self.append,
                     'bin_%d' % ibin,
                     'divBinWidth.%s' % self.output_fmt]
            filename = '_'.join([x for x in parts if x])
            plot.save("%s/%s" % (self.output_dir, filename))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",
                        help="Source directory (should be the one made by unfolding.py)")
    parser.add_argument("--angles",
                        choices=list(qgc.VAR_UNFOLD_DICT.keys()) + ["all"],
                        nargs='+',
                        help="Lambda angles to unfold, or 'all' for all of them (default).",
                        default=["all"])
    args = parser.parse_args()

    if args.angles[0] == "all":
        angles = qgc.COMMON_VARS
    else:
        angles = [a for a in qgc.COMMON_VARS if a.var in args.angles]

    if len(angles) == 0:
        raise RuntimeError("Need at least 1 angle")

    jet_algo = "AK4 PF PUPPI"
    if "ak8puppi" in args.source:
        jet_algo = "AK8 PF PUPPI"

    has_data = not ('_MC_all' in args.source or '_MC_split' in args.source)

    grooming_status = [False, True][:1]
    num_all_iterations = len(grooming_status) * len(angles) * 2 # 2 for dijet & Z+Jet in each angle/grooming variation
    counter = 1

    # Iterate through (un)groomed, variables
    for do_grooming in grooming_status:

        for angle in angles:

            unpickled_regions = []

            # get all regions for this plot
            dijet_region = get_dijet_config(args.source, central=True, groomed=do_grooming)
            zpj_region = get_zpj_config(args.source, groomed=do_grooming)

            def _get_unpicked_region_dict(region):
                region_dir = os.path.join(args.source, region['name'])
                if not os.path.isdir(region_dir):
                    raise RuntimeError("! Warning ! cannot find region dir", region_dir, '- skipping region')

                global counter
                append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name
                print("*"*120)
                print("Algo/Region/Var: %s %s %s (%d/%d)" % (jet_algo, region['name'], angle.var, counter, num_all_iterations))
                print("*"*120)
                counter += 1

                angle_output_dir = "%s/%s/%s" % (args.source, region['name'], angle.var)
                if not os.path.isdir(angle_output_dir):
                    raise RuntimeError("! Warning ! cannot find angle dir", angle_output_dir, '- skipping angle', angle.var)

                # Get region dict from pickle file
                pickle_filename = os.path.join(angle_output_dir, "unfolding_result.pkl")
                unpickled_region = unpickle_region(pickle_filename)
                if unpickled_region is None:
                    raise RuntimeError("Unpickling for %s gave None" % append)

                # check
                if region['name'] != unpickled_region['name']:
                    raise RuntimeError("Mismatch region name")

                return unpickled_region

            dijet_region = _get_unpicked_region_dict(dijet_region)
            zpj_region = _get_unpicked_region_dict(zpj_region)

            # Now make plots
            groom_str = "_groomed" if do_grooming else ""
            output_dir = "%s/detector_unfolded_dijet_zpj%s/%s" % (args.source, groom_str, angle.var)
            binned_plotter = DijetZPJGenPtBinnedPlotter(jet_algo=jet_algo,
                                                        angle=angle,
                                                        bins=dijet_region['unfolder'].pt_bin_edges_gen,
                                                        dijet_region=dijet_region,
                                                        zpj_region=zpj_region,
                                                        output_dir=output_dir,
                                                        is_groomed=do_grooming)

            # dijet-only
            # binned_plotter.plot_detector_unfolded(do_dijet=True, do_zpj=False)

            # ZPJ-only
            # binned_plotter.plot_detector_unfolded(do_dijet=False, do_zpj=True)

            # dijet+ZPJ
            binned_plotter.plot_detector_unfolded(do_dijet=True, do_zpj=True)


