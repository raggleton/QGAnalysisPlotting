#!/usr/bin/env python


"""
Do all the unfolding plots: per pT bin, per lambda bin, summary plot
"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse
import pandas as pd
import numpy as np
from itertools import product, chain
from array import array
from math import sqrt
from copy import copy

import jax.numpy as np
from jax import grad, jit

# make things blazingly fast
import uproot

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from my_unfolder import MyUnfolder, HistBinChopper, unpickle_region, unpack_slim_unfolding_root_file
from my_unfolder_plotter import MyUnfolderPlotter
from unfolding_config import get_dijet_config, get_zpj_config


ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()

# import hunter
# hunter.trace(module='my_unfolder', action=hunter.CallPrinter)


COMMON_STYLE_DICT = {
    "line_width": 2,
    "dijet_cen_color": ROOT.kBlack,
    "dijet_fwd_color": ROOT.kRed,
    "zpj_color": ROOT.kBlue,
    "jet_pt_units_str": "Jet p_{T} [GeV]",

    "data_line_style": 1,
    "data_color": ROOT.kBlack,

    "mc_line_style": 2,
    # format: 3ijk,
    # i=distance between lines,
    # j=angle between 0 and 90 degrees (5 = not drawn),
    # k=angle between 90 and 180 degrees (5 = not drawn)
    "mc_fill_style": 3445,
    "mc_color": ROOT.kBlue+1,

    "mc_alt_line_style": 3,
    "mc_alt_fill_style": 3454,
    "mc_alt_color": ROOT.kAzure+1,
}


def create_angle_label(angle, do_groomed=False):
    angle_prepend = "Groomed " if do_groomed else ""
    if "charged" in angle.var:
        if angle_prepend != "":
            angle_prepend += "c"
        else:
            angle_prepend += "C"
        angle_prepend += "harged-only "

    # this_angle_name = angle.name
    # if (angle_prepend != ""
    #     and 'LHA' not in this_angle_name
    #     and "_{T}" not in this_angle_name
    #     and "PUPPI" not in this_angle_name):
    #     # lower case if Groomed..., but be careful of e.g. pTD, LHA
    #     this_angle_name = this_angle_name[0].lower() + this_angle_name[1:]

    # for plot axis titles
    angle_str = "{prepend}{lambda_str}".format(prepend=angle_prepend,
                                                 lambda_str=angle.lambda_str)
    return angle_str


class SummaryPlotter(object):
    """Do lots of summary plots"""

    def __init__(self, jet_algos, regions, angles, pt_bins_dijet, pt_bins_zpj, df, output_dir, has_data):
        if len(jet_algos) == 0:
            raise RuntimeError("jet_algos is empty")
        self.jet_algos = jet_algos
        if len(regions) == 0:
            raise RuntimeError("regions is empty")
        self.regions = regions
        if len(angles) == 0:
            raise RuntimeError("angles is empty")
        self.angles = angles
        self.pt_bins_dijet = pt_bins_dijet
        self.pt_bins_zpj = pt_bins_zpj
        self.df = df
        self.output_fmt = 'pdf'
        self.output_dir = output_dir
        self.has_data = has_data
        self.is_preliminary = True
        self.mc_label = 'MG5+Pythia8'
        self.alt_mc_label = 'Herwig++'

    @staticmethod
    def data_to_hist(data, data_err, bins):
        if len(data) != len(bins)-1:
            print("data:", data)
            print("bins:", bins)
            raise RuntimeError("len(data) != len(bins)-1")
        h = ROOT.TH1D("h_"+cu.get_unique_str(), "", len(bins)-1, array('d', bins))
        for ind, (y, err) in enumerate(zip(data, data_err), 1):
            h.SetBinContent(ind, y)
            h.SetBinError(ind, err)
        return h

    @staticmethod
    def _generate_filename_prefix(do_dijet, do_zpj):
        this_str = ""
        if do_dijet:
            this_str += "dijet_"
        if do_zpj:
            this_str += "zpj_"
        return this_str

    def plot_dijet_zpj_means_vs_pt_all(self):
        """Plot mean vs pt for dijet (cen+fwd) & Z+jet on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_zpj_means_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_means_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_zpj_means_vs_pt_all' % self.output_dir)

    def plot_dijet_means_vs_pt_all(self):
        """Plot mean vs pt for dijet (cen+fwd) on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_means_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_means_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_means_vs_pt_all' % self.output_dir, do_zpj=False)

    def plot_zpj_means_vs_pt_all(self):
        """Plot mean vs pt for zpj on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_zpj_means_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_means_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_zpj_means_vs_pt_all' % self.output_dir, do_dijet=False)

    def plot_dijet_zpj_means_vs_pt_one_angle_one_jet(self, angle, jet_algo, do_groomed, output_dir, do_zpj=True, do_dijet=True):
        """Do plot of mean lambda vs pt, for the dijet cen+fwd and zpj regions, for a given angle/jet algo/grooming"""
        df = self.df
        mask = ((df['angle'] == angle.var) & (df['jet_algo'] == jet_algo['name']) & (df['isgroomed'] == do_groomed))
        if not mask.any():
            return

        region_name = 'Dijet_central'
        if do_groomed:
            region_name += "_groomed"
        dijet_central_data = df[mask & (df['region'] == region_name)]
        dijet_central_hist = self.data_to_hist(dijet_central_data['mean'], dijet_central_data['mean_err'], self.pt_bins_dijet)
        dijet_central_hist_truth = self.data_to_hist(dijet_central_data['mean_truth'], dijet_central_data['mean_err_truth'], self.pt_bins_dijet)
        dijet_central_hist_alt_truth = self.data_to_hist(dijet_central_data['mean_alt_truth'], dijet_central_data['mean_err_alt_truth'], self.pt_bins_dijet)

        # if "LHA" in angle.var:
        #     print(dijet_central_data)
        #     for i in range(1, dijet_central_hist.GetNbinsX()+1):
        #         print(i, dijet_central_hist.GetBinContent(i), dijet_central_hist.GetBinError(i))


        region_name = 'Dijet_forward'
        if do_groomed:
            region_name += "_groomed"
        dijet_forward_data = df[mask & (df['region'] == region_name)]
        dijet_forward_hist = self.data_to_hist(dijet_forward_data['mean'], dijet_forward_data['mean_err'], self.pt_bins_dijet)
        dijet_forward_hist_truth = self.data_to_hist(dijet_forward_data['mean_truth'], dijet_forward_data['mean_err_truth'], self.pt_bins_dijet)
        dijet_forward_hist_alt_truth = self.data_to_hist(dijet_forward_data['mean_alt_truth'], dijet_forward_data['mean_err_alt_truth'], self.pt_bins_dijet)

        dijet_cen_col = COMMON_STYLE_DICT['dijet_cen_color']
        dijet_fwd_col = COMMON_STYLE_DICT['dijet_fwd_color']
        zpj_col = COMMON_STYLE_DICT['zpj_color']

        # Create copy with 0 error bars, for the ratio subplot
        # The data will have its own error region
        dijet_central_hist_no_errors = dijet_central_hist.Clone()
        dijet_forward_hist_no_errors = dijet_forward_hist.Clone()
        for h in [dijet_central_hist_no_errors, dijet_forward_hist_no_errors]:
            for i in range(1, h.GetNbinsX()+1):
                h.SetBinError(i, 0)

        # Create hists for data with error reigon for ratio
        # Easiest way to get errors right is to do data (with 0 errors)
        # and divide by data (with errors), as if you had MC = data with 0 error
        dijet_central_hist_ratio_error = dijet_central_hist_no_errors.Clone()
        dijet_central_hist_ratio_error.Divide(dijet_central_hist)
        dijet_central_hist_ratio_error.SetFillStyle(3245)
        dijet_central_hist_ratio_error.SetFillColor(dijet_cen_col)
        dijet_central_hist_ratio_error.SetLineWidth(0)
        dijet_central_hist_ratio_error.SetMarkerSize(0)

        dijet_forward_hist_ratio_error = dijet_forward_hist_no_errors.Clone()
        dijet_forward_hist_ratio_error.Divide(dijet_forward_hist)
        dijet_forward_hist_ratio_error.SetFillStyle(3254)
        dijet_forward_hist_ratio_error.SetFillColor(dijet_fwd_col)
        dijet_forward_hist_ratio_error.SetLineWidth(0)
        dijet_forward_hist_ratio_error.SetMarkerSize(0)

        if do_zpj:
            region_name = 'ZPlusJets'
            if do_groomed:
                region_name += "_groomed"
            # drop last pt bin as massive error
            zpj_data = df[mask & (df['region'] == region_name) & (df['pt_bin'] < (len(self.pt_bins_zpj)-3))]
            zpj_hist = self.data_to_hist(zpj_data['mean'], zpj_data['mean_err'], self.pt_bins_zpj[:-2])
            zpj_hist_truth = self.data_to_hist(zpj_data['mean_truth'], zpj_data['mean_err_truth'], self.pt_bins_zpj[:-2])
            zpj_hist_alt_truth = self.data_to_hist(zpj_data['mean_alt_truth'], zpj_data['mean_err_alt_truth'], self.pt_bins_zpj[:-2])
            # 0 error bar hists & ratio = 1 hists for subplot
            zpj_hist_no_errors = zpj_hist.Clone()
            for i in range(1, zpj_hist_no_errors.GetNbinsX()+1):
                zpj_hist_no_errors.SetBinError(i, 0)
            zpj_hist_ratio_error = zpj_hist_no_errors.Clone()
            zpj_hist_ratio_error.Divide(zpj_hist)
            zpj_hist_ratio_error.SetFillStyle(3003)
            zpj_hist_ratio_error.SetFillColor(zpj_col)
            zpj_hist_ratio_error.SetLineWidth(0)
            zpj_hist_ratio_error.SetMarkerSize(0)

        m_size = 1
        lw = COMMON_STYLE_DICT['line_width']
        entries = []
        # Spaces in legend labels are important for padding
        # Add data
        if do_dijet:
            entries.extend([
                Contribution(dijet_central_hist, label=' Dijet (central)',
                             line_color=dijet_cen_col, line_width=lw,
                             marker_color=dijet_cen_col, marker_style=cu.Marker.get('circle', True), marker_size=m_size),
                Contribution(dijet_forward_hist, label=' Dijet (forward)',
                             line_color=dijet_fwd_col, line_width=lw,
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', True), marker_size=m_size),
            ])
        if do_zpj:
            entries.extend([
                Contribution(zpj_hist, label=' Z+jet',
                             line_color=zpj_col, line_width=lw,
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', True), marker_size=m_size),
            ])
        # Add nominal MC
        if do_dijet:
            entries.extend([
                # TODO: make truth plotting optional, also plot alternate generators
                Contribution(dijet_central_hist_truth, label='#splitline{ Dijet (central)  }{ [%s]}' % (self.mc_label),
                             line_color=dijet_cen_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_line_style'],
                             marker_color=dijet_cen_col, marker_style=cu.Marker.get('circle', False), marker_size=0,
                             subplot=dijet_central_hist_no_errors),
                Contribution(dijet_forward_hist_truth, label='#splitline{ Dijet (forward)  }{ [%s]}' % (self.mc_label),
                             line_color=dijet_fwd_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_line_style'],
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', False), marker_size=0,
                             subplot=dijet_forward_hist_no_errors),
            ])
        if do_zpj:
            entries.extend([
                Contribution(zpj_hist_truth, label='#splitline{ Z+jet  }{ [%s]}' % (self.mc_label),
                             line_color=zpj_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_line_style'],
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', False), marker_size=0,
                             subplot=zpj_hist_no_errors),
            ])
        # add alt MC
        if do_dijet:
            entries.extend([
                Contribution(dijet_central_hist_alt_truth, label='#splitline{ Dijet (central)  }{ [%s]}' % (self.alt_mc_label),
                             line_color=dijet_cen_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                             marker_color=dijet_cen_col, marker_style=cu.Marker.get('circle', False), marker_size=0,
                             subplot=dijet_central_hist_no_errors),
                Contribution(dijet_forward_hist_alt_truth, label='#splitline{ Dijet (forward)  }{ [%s]}' % (self.alt_mc_label),
                             line_color=dijet_fwd_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', False), marker_size=0,
                             subplot=dijet_forward_hist_no_errors),
            ])
        if do_zpj:
            entries.extend([
                Contribution(zpj_hist_alt_truth, label='#splitline{ Z+jet}{ [%s]}' % (self.alt_mc_label),
                             line_color=zpj_col, line_width=2, line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', False), marker_size=0,
                             subplot=zpj_hist_no_errors),
            ])

        # for plot axis titles
        angle_str = "#LT %s #GT" % create_angle_label(angle, do_groomed)

        h_max = max([c.obj.GetMaximum() for c in entries])
        h_min = min([c.obj.GetMinimum(1E-10) for c in entries])
        h_range = h_max - h_min
        ylim = (max(0, h_min-(h_range*0.2)), h_max + (h_range*0.8))
        plot = Plot(entries,
                    what='hist',
                    xtitle=COMMON_STYLE_DICT['jet_pt_units_str'],
                    ytitle=angle_str,
                    title="%s jets" % (jet_algo['label']),
                    # ylim=(0, h_max*1.75),
                    # ylim=(h_min*0.75, h_max*1.5),
                    ylim=ylim,
                    has_data=self.has_data,
                    is_preliminary=self.is_preliminary,
                    subplot_type='ratio',
                    subplot_title='MC / Data',
                    subplot_limits=(0.5, 1.5) if self.has_data else (0.9, 1.1)
                    )
        # plot.default_canvas_size = (700, 600)
        plot.title_start_y = 0.85
        plot.title_left_offset = 0.05
        plot.title_font_size = 0.035
        plot.legend.SetX1(0.55)
        plot.legend.SetX2(0.78)
        plot.legend.SetY1(0.68)
        plot.legend.SetY2(0.92)
        if len(entries) > 3:
            plot.legend.SetNColumns(2)
            plot.legend.SetX1(0.50)
            plot.legend.SetY1(0.68)
            plot.legend.SetX2(0.92)
            plot.legend.SetY2(0.92)
            # plot.legend.SetBorderSize(1)
            # plot.legend.SetLineColor(ROOT.kBlack)
            plot.title_left_offset = 0.02
        if len(entries) > 6:
            plot.legend.SetNColumns(3)
        plot.legend.SetY2(0.87)
        plot.left_margin = 0.16
        plot.subplot_line_style = 1
        plot.y_padding_max_linear = 1.9
        subplot_draw_opts = "NOSTACK E1"
        plot.plot("NOSTACK E1", subplot_draw_opts)
        # plot.get_modifier().GetYaxis().SetTitleOffset(plot.get_modifier().GetYaxis().GetTitleOffset()*1.1)
        plot.set_logx(do_more_labels=False)

        # now draw the data error shaded area
        # this is a bit hacky - basically draw them on the ratio pad,
        # then redraw the existing hists & line to get them ontop
        # note that we use "same" for all - this is to keep the original axes
        # (we may want to rethink this later?)
        plot.subplot_pad.cd()
        draw_opt = "E2 SAME"
        if do_dijet:
            dijet_central_hist_ratio_error.Draw(draw_opt)
            dijet_forward_hist_ratio_error.Draw(draw_opt)
        if do_zpj:
            zpj_hist_ratio_error.Draw(draw_opt)
        plot.subplot_container.Draw("SAME" + subplot_draw_opts)
        plot.subplot_line.Draw()
        plot.canvas.cd()

        prefix = self._generate_filename_prefix(do_dijet, do_zpj)
        groomed_str = '_groomed' if do_groomed else ''
        plot.save("%s/%smeans_vs_pt_%s%s_%s.%s" % (output_dir, prefix, angle.var, groomed_str, jet_algo['name'], self.output_fmt))

    def plot_dijet_zpj_rms_vs_pt_all(self):
        """Plot RMS vs pt for dijet (cen+fwd) & Z+jet on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_zpj_rms_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_rms_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_zpj_rms_vs_pt_all' % self.output_dir)

    def plot_dijet_rms_vs_pt_all(self):
        """Plot RMS vs pt for dijet (cen+fwd) on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_rms_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_rms_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_rms_vs_pt_all' % self.output_dir, do_zpj=False)

    def plot_zpj_rms_vs_pt_all(self):
        """Plot RMS vs pt for zpj on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_zpj_rms_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_rms_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_zpj_rms_vs_pt_all' % self.output_dir, do_dijet=False)

    def plot_dijet_zpj_rms_vs_pt_one_angle_one_jet(self, angle, jet_algo, do_groomed, output_dir, do_zpj=True, do_dijet=True):
        """Do plot of RMS lambda vs pt, for the dijet cen+fwd and zpj regions, for a given angle/jet algo/grooming"""
        df = self.df
        mask = ((df['angle'] == angle.var) & (df['jet_algo'] == jet_algo['name']) & (df['isgroomed'] == do_groomed))
        if not mask.any():
            return

        region_name = 'Dijet_central'
        if do_groomed:
            region_name += "_groomed"
        dijet_central_data = df[mask & (df['region'] == region_name)]
        dijet_central_hist = self.data_to_hist(dijet_central_data['rms'], dijet_central_data['rms_err'], self.pt_bins_dijet)
        dijet_central_hist_truth = self.data_to_hist(dijet_central_data['rms_truth'], dijet_central_data['rms_err_truth'], self.pt_bins_dijet)
        dijet_central_hist_alt_truth = self.data_to_hist(dijet_central_data['rms_alt_truth'], dijet_central_data['rms_err_alt_truth'], self.pt_bins_dijet)

        region_name = 'Dijet_forward'
        if do_groomed:
            region_name += "_groomed"
        dijet_forward_data = df[mask & (df['region'] == region_name)]
        dijet_forward_hist = self.data_to_hist(dijet_forward_data['rms'], dijet_forward_data['rms_err'], self.pt_bins_dijet)
        dijet_forward_hist_truth = self.data_to_hist(dijet_forward_data['rms_truth'], dijet_forward_data['rms_err_truth'], self.pt_bins_dijet)
        dijet_forward_hist_alt_truth = self.data_to_hist(dijet_forward_data['rms_alt_truth'], dijet_forward_data['rms_err_alt_truth'], self.pt_bins_dijet)

        dijet_cen_col = COMMON_STYLE_DICT['dijet_cen_color']
        dijet_fwd_col = COMMON_STYLE_DICT['dijet_fwd_color']
        zpj_col = COMMON_STYLE_DICT['zpj_color']

        # Create copy with 0 error bars, for the ratio subplot
        # The data will have its own error region
        dijet_central_hist_no_errors = dijet_central_hist.Clone()
        dijet_forward_hist_no_errors = dijet_forward_hist.Clone()
        for h in [dijet_central_hist_no_errors, dijet_forward_hist_no_errors]:
            for i in range(1, h.GetNbinsX()+1):
                h.SetBinError(i, 0)

        # Create hists for data with error reigon for ratio
        # Easiest way to get errors right is to do data (with 0 errors)
        # and divide by data (with errors), as if you had MC = data with 0 error
        dijet_central_hist_ratio_error = dijet_central_hist_no_errors.Clone()
        dijet_central_hist_ratio_error.Divide(dijet_central_hist)
        dijet_central_hist_ratio_error.SetFillStyle(3245)
        dijet_central_hist_ratio_error.SetFillColor(dijet_cen_col)
        dijet_central_hist_ratio_error.SetLineWidth(0)
        dijet_central_hist_ratio_error.SetMarkerSize(0)

        dijet_forward_hist_ratio_error = dijet_forward_hist_no_errors.Clone()
        dijet_forward_hist_ratio_error.Divide(dijet_forward_hist)
        dijet_forward_hist_ratio_error.SetFillStyle(3254)
        dijet_forward_hist_ratio_error.SetFillColor(dijet_fwd_col)
        dijet_forward_hist_ratio_error.SetLineWidth(0)
        dijet_forward_hist_ratio_error.SetMarkerSize(0)

        if do_zpj:
            region_name = 'ZPlusJets'
            if do_groomed:
                region_name += "_groomed"
            zpj_data = df[mask & (df['region'] == region_name) & (df['pt_bin'] < (len(self.pt_bins_zpj)-3))]
            zpj_hist = self.data_to_hist(zpj_data['rms'], zpj_data['rms_err'], self.pt_bins_zpj[:-2])
            zpj_hist_truth = self.data_to_hist(zpj_data['rms_truth'], zpj_data['rms_err_truth'], self.pt_bins_zpj[:-2])
            zpj_hist_alt_truth = self.data_to_hist(zpj_data['rms_alt_truth'], zpj_data['rms_err_alt_truth'], self.pt_bins_zpj[:-2])
            # 0 error bar hists & ratio = 1 hists for subplot
            zpj_hist_no_errors = zpj_hist.Clone()
            for i in range(1, zpj_hist_no_errors.GetNbinsX()+1):
                zpj_hist_no_errors.SetBinError(i, 0)
            zpj_hist_ratio_error = zpj_hist_no_errors.Clone()
            zpj_hist_ratio_error.Divide(zpj_hist)
            zpj_hist_ratio_error.SetFillStyle(3003)
            zpj_hist_ratio_error.SetFillColor(zpj_col)
            zpj_hist_ratio_error.SetLineWidth(0)
            zpj_hist_ratio_error.SetMarkerSize(0)

        m_size = 1
        lw = COMMON_STYLE_DICT['line_width']
        entries = []
        # Spaces in legend labels are important for padding
        # Add data
        if do_dijet:
            entries.extend([
                Contribution(dijet_central_hist, label=' Dijet (central)',
                             line_color=dijet_cen_col, line_width=lw,
                             marker_color=dijet_cen_col, marker_style=cu.Marker.get('circle', True), marker_size=m_size),
                Contribution(dijet_forward_hist, label=' Dijet (forward)',
                             line_color=dijet_fwd_col, line_width=lw,
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', True), marker_size=m_size),
            ])
        if do_zpj:
            entries.extend([
                Contribution(zpj_hist, label=' Z+jet',
                             line_color=zpj_col, line_width=lw,
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', True), marker_size=m_size),
            ])
        # Add nominal MC
        if do_dijet:
            entries.extend([
                # TODO: make truth plotting optional, also plot alternate generators
                Contribution(dijet_central_hist_truth, label='#splitline{ Dijet (central)  }{ [%s]}' % (self.mc_label),
                             line_color=dijet_cen_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_line_style'],
                             marker_color=dijet_cen_col, marker_style=cu.Marker.get('circle', False), marker_size=0,
                             subplot=dijet_central_hist_no_errors),
                Contribution(dijet_forward_hist_truth, label='#splitline{ Dijet (forward)  }{ [%s]}' % (self.mc_label),
                             line_color=dijet_fwd_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_line_style'],
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', False), marker_size=0,
                             subplot=dijet_forward_hist_no_errors),
            ])
        if do_zpj:
            entries.extend([
                Contribution(zpj_hist_truth, label='#splitline{ Z+jet  }{ [%s]}' % (self.mc_label),
                             line_color=zpj_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_line_style'],
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', False), marker_size=0,
                             subplot=zpj_hist_no_errors),
            ])
        # add alt MC
        if do_dijet:
            entries.extend([
                Contribution(dijet_central_hist_alt_truth, label='#splitline{ Dijet (central)  }{ [%s]}' % (self.alt_mc_label),
                             line_color=dijet_cen_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                             marker_color=dijet_cen_col, marker_style=cu.Marker.get('circle', False), marker_size=0,
                             subplot=dijet_central_hist_no_errors),
                Contribution(dijet_forward_hist_alt_truth, label='#splitline{ Dijet (forward)  }{ [%s]}' % (self.alt_mc_label),
                             line_color=dijet_fwd_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', False), marker_size=0,
                             subplot=dijet_forward_hist_no_errors),
            ])
        if do_zpj:
            entries.extend([
                Contribution(zpj_hist_alt_truth, label='#splitline{ Z+jet}{ [%s]}' % (self.alt_mc_label),
                             line_color=zpj_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', False), marker_size=0,
                             subplot=zpj_hist_no_errors),
            ])

        angle_str = "RMS %s" % create_angle_label(angle, do_groomed)

        h_max = max([c.obj.GetMaximum() for c in entries])
        plot = Plot(entries,
                    what='hist',
                    xtitle=COMMON_STYLE_DICT['jet_pt_units_str'],
                    ytitle=angle_str,
                    title="%s jets" % (jet_algo['label']),
                    # ylim=(0, h_max*1.75),
                    has_data=self.has_data,
                    is_preliminary=self.is_preliminary,
                    subplot_type='ratio',
                    subplot_title='MC / Data',
                    subplot_limits=(0.5, 1.5) if self.has_data else (0.9, 1.1))
        # plot.default_canvas_size = (700, 600)
        plot.title_start_y = 0.85
        plot.title_left_offset = 0.05
        plot.title_font_size = 0.035
        plot.legend.SetX1(0.55)
        plot.legend.SetX2(0.78)
        plot.legend.SetY1(0.68)
        plot.legend.SetY2(0.92)
        if len(entries) > 3:
            plot.legend.SetNColumns(2)
            plot.legend.SetX1(0.50)
            plot.legend.SetY1(0.68)
            plot.legend.SetX2(0.92)
            plot.legend.SetY2(0.92)
            # plot.legend.SetBorderSize(1)
            # plot.legend.SetLineColor(ROOT.kBlack)
            plot.title_left_offset = 0.02
        if len(entries) > 6:
            plot.legend.SetNColumns(3)
        plot.legend.SetY2(0.87)
        plot.left_margin = 0.16
        plot.subplot_line_style = 1
        plot.y_padding_max_linear = 1.9
        # plot.default_canvas_size = (600, 800)
        subplot_draw_opts = "NOSTACK E1"
        plot.plot("NOSTACK E1", subplot_draw_opts)
        # plot.get_modifier().GetYaxis().SetTitleOffset(plot.get_modifier().GetYaxis().GetTitleOffset()*1.1)
        plot.set_logx(do_more_labels=False)

         # now draw the data error shaded area
        # this is a bit hacky - basically draw them on the ratio pad,
        # then redraw the existing hists & line to get them ontop
        # note that we use "same" for all - this is to keep the original axes
        # (we may want to rethink this later?)
        plot.subplot_pad.cd()
        draw_opt = "E2 SAME"
        if do_dijet:
            dijet_central_hist_ratio_error.Draw(draw_opt)
            dijet_forward_hist_ratio_error.Draw(draw_opt)
        if do_zpj:
            zpj_hist_ratio_error.Draw(draw_opt)
        plot.subplot_container.Draw("SAME" + subplot_draw_opts)
        plot.subplot_line.Draw()
        plot.canvas.cd()

        prefix = self._generate_filename_prefix(do_dijet, do_zpj)
        groomed_str = '_groomed' if do_groomed else ''
        plot.save("%s/%srms_vs_pt_%s%s_%s.%s" % (output_dir, prefix, angle.var, groomed_str, jet_algo['name'], self.output_fmt))

    def plot_dijet_zpj_delta_vs_pt_all(self):
        """Plot delta vs pt for dijet (cen+fwd) & Z+jet on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_zpj_delta_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_delta_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_zpj_delta_vs_pt_all' % self.output_dir)

    def plot_dijet_delta_vs_pt_all(self):
        """Plot mean vs pt for dijet (cen+fwd) on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_delta_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_delta_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_delta_vs_pt_all' % self.output_dir, do_zpj=False)

    def plot_zpj_delta_vs_pt_all(self):
        """Plot mean vs pt for zpj on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_zpj_delta_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_delta_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_zpj_delta_vs_pt_all' % self.output_dir, do_dijet=False)

    def plot_dijet_zpj_delta_vs_pt_one_angle_one_jet(self, angle, jet_algo, do_groomed, output_dir, do_zpj=True, do_dijet=True):
        """Do plot of delta lambda vs pt, for the dijet cen+fwd and zpj regions, for a given angle/jet algo/grooming"""
        df = self.df
        mask = ((df['angle'] == angle.var) & (df['jet_algo'] == jet_algo['name']) & (df['isgroomed'] == do_groomed))
        if not mask.any():
            return

        region_name = 'Dijet_central'
        if do_groomed:
            region_name += "_groomed"
        dijet_central_data = df[mask & (df['region'] == region_name)]
        dijet_central_hist_truth = self.data_to_hist(dijet_central_data['delta_truth'], dijet_central_data['delta_err_truth'], self.pt_bins_dijet)
        dijet_central_hist_alt_truth = self.data_to_hist(dijet_central_data['delta_alt_truth'], dijet_central_data['delta_err_alt_truth'], self.pt_bins_dijet)

        region_name = 'Dijet_forward'
        if do_groomed:
            region_name += "_groomed"
        dijet_forward_data = df[mask & (df['region'] == region_name)]
        dijet_forward_hist_truth = self.data_to_hist(dijet_forward_data['delta_truth'], dijet_forward_data['delta_err_truth'], self.pt_bins_dijet)
        dijet_forward_hist_alt_truth = self.data_to_hist(dijet_forward_data['delta_alt_truth'], dijet_forward_data['delta_err_alt_truth'], self.pt_bins_dijet)

        dijet_cen_col = COMMON_STYLE_DICT['dijet_cen_color']
        dijet_fwd_col = COMMON_STYLE_DICT['dijet_fwd_color']
        zpj_col = COMMON_STYLE_DICT['zpj_color']

        if do_zpj:
            region_name = 'ZPlusJets'
            if do_groomed:
                region_name += "_groomed"
            zpj_data = df[mask & (df['region'] == region_name) & (df['pt_bin'] < (len(self.pt_bins_zpj)-3))]
            zpj_hist_truth = self.data_to_hist(zpj_data['delta_truth'], zpj_data['delta_err_truth'], self.pt_bins_zpj[:-2])
            zpj_hist_alt_truth = self.data_to_hist(zpj_data['delta_alt_truth'], zpj_data['delta_err_alt_truth'], self.pt_bins_zpj[:-2])

        m_size = 1
        lw = COMMON_STYLE_DICT['line_width']
        entries = []
        # Spaces in legend labels are important for padding
        # Add nominal MC
        if do_dijet:
            entries.extend([
                # TODO: make truth plotting optional, also plot alternate generators
                Contribution(dijet_central_hist_truth, label='#splitline{ Dijet (central)  }{ [%s]}' % (self.mc_label),
                             line_color=dijet_cen_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_line_style'],
                             marker_color=dijet_cen_col, marker_style=cu.Marker.get('circle', False), marker_size=0,
                             ),
                Contribution(dijet_forward_hist_truth, label='#splitline{ Dijet (forward)  }{ [%s]}' % (self.mc_label),
                             line_color=dijet_fwd_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_line_style'],
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', False), marker_size=0,
                             ),
            ])
        if do_zpj:
            entries.extend([
                Contribution(zpj_hist_truth, label='#splitline{ Z+jet  }{ [%s]}' % (self.mc_label),
                             line_color=zpj_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_line_style'],
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', False), marker_size=0,
                             ),
            ])
        # add alt MC
        if do_dijet:
            entries.extend([
                Contribution(dijet_central_hist_alt_truth, label='#splitline{ Dijet (central)  }{ [%s]}' % (self.alt_mc_label),
                             line_color=dijet_cen_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                             marker_color=dijet_cen_col, marker_style=cu.Marker.get('circle', False), marker_size=0,
                             ),
                Contribution(dijet_forward_hist_alt_truth, label='#splitline{ Dijet (forward)  }{ [%s]}' % (self.alt_mc_label),
                             line_color=dijet_fwd_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', False), marker_size=0,
                             ),
            ])
        if do_zpj:
            entries.extend([
                Contribution(zpj_hist_alt_truth, label='#splitline{ Z+jet}{ [%s]}' % (self.alt_mc_label),
                             line_color=zpj_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', False), marker_size=0,
                             ),
            ])

        angle_str = "#Delta, %s" % create_angle_label(angle, do_groomed)

        h_max = max([c.obj.GetMaximum() for c in entries])
        plot = Plot(entries,
                    what='hist',
                    xtitle=COMMON_STYLE_DICT['jet_pt_units_str'],
                    ytitle=angle_str,
                    title="%s jets" % (jet_algo['label']),
                    ylim=(0, h_max*1.75),
                    has_data=self.has_data,
                    is_preliminary=self.is_preliminary)
        # plot.default_canvas_size = (700, 600)
        plot.title_start_y = 0.85
        plot.title_left_offset = 0.05
        plot.title_font_size = 0.035
        plot.legend.SetX1(0.55)
        plot.legend.SetX2(0.78)
        plot.legend.SetY1(0.68)
        plot.legend.SetY2(0.92)
        if len(entries) > 3:
            plot.legend.SetNColumns(2)
            plot.legend.SetX1(0.50)
            plot.legend.SetY1(0.68)
            plot.legend.SetX2(0.92)
            plot.legend.SetY2(0.92)
            # plot.legend.SetBorderSize(1)
            # plot.legend.SetLineColor(ROOT.kBlack)
            plot.title_left_offset = 0.02
        if len(entries) > 6:
            plot.legend.SetNColumns(3)
        plot.legend.SetY2(0.87)
        plot.left_margin = 0.16
        plot.subplot_line_style = 1
        # plot.y_padding_max_linear = 1.9
        # plot.y_padding_min_linear = 20
        # plot.default_canvas_size = (600, 800)
        plot.plot("NOSTACK HIST E1")
        # plot.get_modifier().GetYaxis().SetTitleOffset(plot.get_modifier().GetYaxis().GetTitleOffset()*1.1)
        plot.set_logx(do_more_labels=False)

        prefix = self._generate_filename_prefix(do_dijet, do_zpj)
        groomed_str = '_groomed' if do_groomed else ''
        plot.save("%s/%sdelta_vs_pt_%s%s_%s.%s" % (output_dir, prefix, angle.var, groomed_str, jet_algo['name'], self.output_fmt))

    @staticmethod
    def _make_hist_from_values(value_error_pairs, bins=None, title="", name="", bin_names=None):
        name = name or cu.get_unique_str()
        if bins is not None:
            nbins = len(bins)-1
            if len(value_error_pairs) != nbins:
                return ValueError("len(value) != len(bins)-1")
            h = ROOT.TH1D(name, title, bins)
        else:
            nbins = len(value_error_pairs)
            h = ROOT.TH1D(name, title, nbins, 0, nbins)
        if bin_names is not None:
            if len(bin_names) != nbins:
                raise ValueError("len(bin_names) != len(value_error_pairs)")
        for i, (v, e) in enumerate(value_error_pairs, 1):
            if bin_names is not None:
                h.GetXaxis().SetBinLabel(i, bin_names[i-1])
            h.SetBinContent(i, v)
            h.SetBinError(i, e)
        # h.SetDirectory(0)
        return h

    @staticmethod
    def _style_data_hist(hist):
        hist.SetLineStyle(COMMON_STYLE_DICT['data_line_style'])

    @staticmethod
    def _style_mc_hist(hist):
        hist.SetLineStyle(COMMON_STYLE_DICT['mc_line_style'])
        hist.SetLineColor(COMMON_STYLE_DICT['mc_color'])
        hist.SetFillStyle(COMMON_STYLE_DICT['mc_fill_style'])
        hist.SetFillColor(COMMON_STYLE_DICT['mc_color'])
        hist.SetMarkerColor(COMMON_STYLE_DICT['mc_color'])
        hist.SetMarkerSize(0)

    @staticmethod
    def _style_alt_mc_hist(hist):
        hist.SetLineStyle(COMMON_STYLE_DICT['mc_alt_line_style'])
        hist.SetLineColor(COMMON_STYLE_DICT['mc_alt_color'])
        hist.SetFillStyle(COMMON_STYLE_DICT['mc_alt_fill_style'])
        hist.SetFillColor(COMMON_STYLE_DICT['mc_alt_color'])
        hist.SetMarkerColor(COMMON_STYLE_DICT['mc_alt_color'])
        hist.SetMarkerSize(0)

    @staticmethod
    def calc_hists_max_min(hists):
        y_up, y_down = -999999, 999999
        for h in hists:
            y_up = max(max([h.GetBinContent(i) + h.GetBinError(i) for i in range(1, h.GetNbinsX()+1)]), y_up)
            y_down = min(min([h.GetBinContent(i) - h.GetBinError(i) for i in range(1, h.GetNbinsX()+1)]), y_down)
        return y_up, y_down

    def construct_mean_rms_hist_groups(self, selections):
        """Construct mean & RMS hists for plots

        See plot_mean_rms_bins_summary() docstring about `selections` arg.

        Each of the mean/RMS returns are lists. Each item corresponds to the
        respective item in the outer level of the `selections` arg.
        (e.g. each angle).
        Each of the hist list items is itself a list, corresponding to the hists
        for [data, nominal MC, alt MC].

        Also does styling of hists (as easier here)
        """
        mean_hists, rms_hists = [], []
        for selection_group in selections:

            mean_entries_data = []
            mean_entries_mc = []
            mean_entries_alt_mc = []

            rms_entries_data = []
            rms_entries_mc = []
            rms_entries_alt_mc = []

            bin_names = []

            for query, label in selection_group['selections']:
                # print(query)
                results = self.df.query(query)
                if len(results.index) != 1:
                    print(query)
                    print(results)
                    raise ValueError("Got != 1 results, check query")

                mean_entries_data.append([results['mean'], results['mean_err']])
                mean_entries_mc.append([results['mean_truth'], results['mean_err_truth']])
                mean_entries_alt_mc.append([results['mean_alt_truth'], results['mean_err_alt_truth']])

                rms_entries_data.append([results['rms'], results['rms_err']])
                rms_entries_mc.append([results['rms_truth'], results['rms_err_truth']])
                rms_entries_alt_mc.append([results['rms_alt_truth'], results['rms_err_alt_truth']])

                bin_names.append(label)

            hist_mean_data = self._make_hist_from_values(mean_entries_data, bin_names=bin_names)
            self._style_data_hist(hist_mean_data)
            hist_mean_mc = self._make_hist_from_values(mean_entries_mc, bin_names=bin_names)
            self._style_mc_hist(hist_mean_mc)
            hist_mean_alt_mc = self._make_hist_from_values(mean_entries_alt_mc, bin_names=bin_names)
            self._style_alt_mc_hist(hist_mean_alt_mc)

            mean_hists.append([hist_mean_data, hist_mean_mc, hist_mean_alt_mc])

            hist_rms_data = self._make_hist_from_values(rms_entries_data, bin_names=bin_names)
            self._style_data_hist(hist_rms_data)
            hist_rms_mc = self._make_hist_from_values(rms_entries_mc, bin_names=bin_names)
            self._style_mc_hist(hist_rms_mc)
            hist_rms_alt_mc = self._make_hist_from_values(rms_entries_alt_mc, bin_names=bin_names)
            self._style_alt_mc_hist(hist_rms_alt_mc)

            rms_hists.append([hist_rms_data, hist_rms_mc, hist_rms_alt_mc])

        return mean_hists, rms_hists

    def plot_mean_rms_bins_summary(self, selections, output_file):
        """Make plot of mean & RMS for various selections, showing data, mc, alt mc

        `selections` is a multi-level list
        The first level is of dicts,
            {'label': str, 'selections': list}
        where each dict represents an angle (i.e. its own column).

        The 'selections' key of that dict then provides a list of tuples,
        where each is a bin in the histograms that make up each plot column.

        Each tuple is of the form (query_str, bin_label_str),
        where the query_str gets passed to the dataframe to retrieve a single entry.
        If 0 or >1 entries are found, this raises a ValueError.
        From that entry, we get mean ± error, and RMS ± error.
        That bin in the histogram is then labelled with bin_label_str.

        In this way, we build up the contents of the histograms for each angle.

        The plotting part is complicated, since we need 2 pads (upper & lower)
        for each column. These need to overlap, since otherwise it will chop
        off the y axes labels. However, we want the top x axis of the lower pad,
        and the lower x axis of the upper pad, to align vertically. So we have
        to be careful with the exact pad & margin sizes.

        Note that the margins basically determine where the axes go in the pad,
        and are a fraction of the pad width/height. Axis labels must go in the
        pad margin - nothing can be drawn outside the pad itself.

        We also manually put in the y axis titles, and the angle names, as
        individual text elements, just so we can control their position exactly.
        """
        # Get data for plots
        mean_hists, rms_hists = self.construct_mean_rms_hist_groups(selections)

        gc_stash = [] # to stop stuff being deleted

        # Setup canvas and pads
        canvas = ROOT.TCanvas("c_"+cu.get_unique_str(), "", 1200, 600)
        canvas.SetBottomMargin(0.0)
        canvas.SetTopMargin(0.0)
        canvas.SetLeftMargin(0.0)
        canvas.SetRightMargin(0.0)
        canvas.cd()

        # For each selection group, create 2 pads, one for mean, one for RMS
        mean_pads, rms_pads = [], []
        n_pads = len(selections)

        # gap between right end of plots and edge of canvas, used for legend
        right_margin = 0.12
        # pad_left_titles_gap = 0.01 # gap between pad_left_titles and all plots
        pad_to_pad_gap = 0.005  # gap between plot pad columns
        # how far in from the left the first plotting pad starts. used for y axis title
        left_margin = 0.04
        # figure out width per pad - get total width available, then divide by number of pads
        pad_width = (1 - left_margin - right_margin - (pad_to_pad_gap*(n_pads - 1))) / n_pads

        pad_offset_bottom = 0.01 # spacing between bottom of RMS pad and canvas edge
        pad_offset_top = 0.08  # spacing between top of mean pad and canvas edge - for CMS and lumi text

        # per-pad margins: these determine where the hist axes lie,
        # and are fractions of the **pad** width/height, not the global canvas
        pad_right_margin = 0.02
        pad_left_margin = 0.2

        # bottom margin includes space for x axis labels
        # note that we apply it BOTH to the upper and lower pads,
        # even though the labels are on the lower pads, since we need the pads
        # to be exactly the same size, in order to get man things the same size,
        # e.g. ticks, labels, hatching, all of which depend on pad size,
        # and not histogram axes size!
        pad_bottom_margin = 0.49

        # extra bit to add to the top margins of lower and upper pads
        # to ensure y axis numbers don't get cut off
        pad_top_margin = 0.012

        # pad height is constrained by the available space
        # (i.e. after pad_offset_top and pad_offset_bottom),
        # and the fact that we want the pads to overlap exactly by both the
        # top and bottom margins, to ensure that the x axes align vertically
        pad_height = (1 - pad_offset_top - pad_offset_bottom) / (2 - pad_bottom_margin - pad_top_margin)

        for isel, selection_group in enumerate(selections):
            canvas.cd()
            pad_start_x = left_margin + (isel*pad_to_pad_gap) + (isel*pad_width)
            pad_end_x = pad_start_x + pad_width
            # Create pad for mean hist - upper half of this column
            pad_mean = ROOT.TPad(cu.get_unique_str(), "", pad_start_x, 1-pad_offset_top-pad_height, pad_end_x, 1-pad_offset_top)
            ROOT.SetOwnership(pad_mean, False)
            pad_mean.SetFillColor(isel+2)
            pad_mean.SetFillStyle(3004)
            pad_mean.SetFillStyle(4000)
            pad_mean.SetTopMargin(pad_top_margin)
            pad_mean.SetBottomMargin(pad_bottom_margin)
            pad_mean.SetRightMargin(pad_right_margin)
            pad_mean.SetLeftMargin(pad_left_margin)
            pad_mean.SetTicks(1, 1)
            pad_mean.Draw()
            mean_pads.append(pad_mean)

            canvas.cd()
            # Create pad for mean hist - lower half of this column
            pad_rms = ROOT.TPad(cu.get_unique_str(), "", pad_start_x, pad_offset_bottom, pad_end_x, pad_offset_bottom+pad_height)
            ROOT.SetOwnership(pad_rms, False)
            pad_rms.SetFillColor(isel+2)
            pad_rms.SetFillStyle(3003)
            pad_rms.SetFillStyle(4000)
            pad_rms.SetTopMargin(pad_top_margin)
            pad_rms.SetBottomMargin(pad_bottom_margin)
            pad_rms.SetRightMargin(pad_right_margin)
            pad_rms.SetLeftMargin(pad_left_margin)
            pad_rms.SetTicks(1, 1)
            pad_rms.Draw()
            rms_pads.append(pad_rms)


        # Now draw the histograms
        for isel, (selection_group, mean_pad, mean_hist_group, rms_pad, rms_hist_group) \
            in enumerate(zip(selections, mean_pads, mean_hists,  rms_pads, rms_hists)):
            mean_pad.cd()

            mean_hist_group[2].Draw("E1")
            mean_hist_group[1].Draw("E1 SAME")
            mean_hist_group[0].Draw("E1 SAME")

            mean_draw_hist = mean_hist_group[-1]
            # remove x axis label
            xax = mean_draw_hist.GetXaxis()
            xax.SetLabelSize(0)
            xax.SetTitleSize(0)

            factor = 1.5
            xax.SetTickLength(xax.GetTickLength()*factor)

            yax = mean_draw_hist.GetYaxis()
            label_size_fudge = 1.6
            yax.SetLabelSize(yax.GetLabelSize()*factor*label_size_fudge)
            label_offset_fudge = 4
            yax.SetLabelOffset(yax.GetLabelOffset()*factor*label_offset_fudge)
            tick_fudge = 4
            yax.SetTickLength(yax.GetTickLength()*factor*tick_fudge)
            n_divisions = 1005  # fewer big ticks so less chance of overlapping with lower plot
            # n_divisions = 510  # default
            yax.SetNdivisions(n_divisions)

            # Set range using bin contents + error bars
            y_up, y_down = self.calc_hists_max_min(mean_hist_group)
            y_range = y_up - y_down
            down_padding = 0.2 * y_range
            up_padding = 0.2 * y_range
            mean_draw_hist.SetMinimum(y_down - down_padding)
            mean_draw_hist.SetMaximum(y_up + up_padding)

            rms_pad.cd()
            rms_hist_group[2].Draw("E1")
            rms_hist_group[1].Draw("E1 SAME")
            rms_hist_group[0].Draw("E1 SAME")

            rms_draw_hist = rms_hist_group[-1]
            xax = rms_draw_hist.GetXaxis()
            xax.CenterLabels()
            xax.LabelsOption("v")

            xax.SetTickLength(xax.GetTickLength()*factor)
            xax.SetLabelSize(xax.GetLabelSize()*factor*2)
            xax.SetLabelOffset(xax.GetLabelOffset()*factor)

            yax = rms_draw_hist.GetYaxis()
            yax.SetLabelSize(yax.GetLabelSize()*factor*label_size_fudge)
            yax.SetLabelOffset(yax.GetLabelOffset()*factor*label_offset_fudge)
            yax.SetTickLength(yax.GetTickLength()*factor*tick_fudge)  # extra bit of fudging, probably becasue pads are sligthly different sizes
            yax.SetNdivisions(n_divisions)

            # Set range using bin contents + error bars
            y_up, y_down = self.calc_hists_max_min(rms_hist_group)
            y_range = y_up - y_down
            down_padding = 0.2 * y_range
            up_padding = 0.22 * y_range
            rms_draw_hist.SetMinimum(y_down - down_padding)
            rms_draw_hist.SetMaximum(y_up + up_padding)

            # Draw variable name
            var_latex = ROOT.TLatex()
            var_latex.SetTextAlign(ROOT.kHAlignCenter + ROOT.kVAlignBottom)
            var_latex.SetTextSize(0.1)
            # var_latex.SetTextFont(42)
            # these are relative to the RMS pad! not the canvas
            # put in middle of the plot (if 0.5, woud look off-centre)
            var_x = 0.5*(1-pad_right_margin-pad_left_margin) + pad_left_margin
            var_y = 0.03  # do by eye
            var_latex.DrawLatexNDC(var_x, var_y, selection_group['label'])
            gc_stash.append(var_latex)

        canvas.cd()

        # Add legend
        # Bit of a hack here - to get the fill hashing style the same as in the plots,
        # we have to put it in a pad of the same size as a plotting pad
        # This is because the hashing separation is scaled by the pad size
        # So we can't just put the legend in the global canvas.
        # We also can't modify the fill style of the legend entries (I tried, it does nothing)
        leg_y_top = 0.93
        leg_x_right = 1-0.0
        leg_pad = ROOT.TPad("leg_pad_"+cu.get_unique_str(), "", leg_x_right-mean_pads[0].GetAbsWNDC(), leg_y_top-mean_pads[0].GetAbsHNDC(), leg_x_right, leg_y_top)
        ROOT.SetOwnership(leg_pad, False)  # important! otherwise seg fault
        # leg_pad.SetFillColor(ROOT.kYellow)
        # leg_pad.SetFillStyle(3004)
        leg_pad.SetFillStyle(4000)
        leg_pad.SetLeftMargin(0)
        leg_pad.SetRightMargin(0)
        leg_pad.SetTopMargin(0)
        leg_pad.SetBottomMargin(0)
        leg_pad.Draw()
        leg_pad.cd()
        gc_stash.append(leg_pad)
        leg = ROOT.TLegend(0.28, mean_pads[0].GetBottomMargin(), 1, 1)
        leg.AddEntry(mean_hists[0][0], "Data" ,"EL")
        le_mc = leg.AddEntry(mean_hists[0][1], self.mc_label ,"L")
        le_mc_alt = leg.AddEntry(mean_hists[0][2], self.alt_mc_label ,"L")
        leg.SetTextSize(0.08)
        leg.Draw()

        canvas.cd()

        # Add mean, RMS text
        text_width = left_margin
        text_x = (0.5 * left_margin) - (0.5*text_width)
        # let ROOT center it, just make the box the height of the axis
        text_y_end = 1 - pad_offset_top - (pad_top_margin*mean_pads[0].GetAbsHNDC())
        text_y = 1 - pad_offset_top - mean_pads[0].GetAbsHNDC() + (mean_pads[0].GetAbsHNDC()*pad_bottom_margin)
        mean_text = ROOT.TPaveText(text_x, text_y, text_x+text_width, text_y_end, "NDC NB")
        mean_text.SetFillStyle(0)
        mean_text.SetBorderSize(0)
        t_mean = mean_text.AddText("Mean")
        t_mean.SetTextAngle(90)
        text_size = 0.04
        t_mean.SetTextSize(text_size)
        mean_text.Draw()

        # let ROOT center it, just make the box the height of the axis
        text_y = (rms_pads[0].GetAbsHNDC() * pad_bottom_margin) + pad_offset_bottom
        text_y_end = (rms_pads[0].GetAbsHNDC() * (1-pad_top_margin)) + pad_offset_bottom
        rms_text = ROOT.TPaveText(text_x, text_y, text_x+text_width, text_y_end, "NDC NB")
        rms_text.SetFillStyle(0)
        rms_text.SetBorderSize(0)
        t_mean = rms_text.AddText("RMS")
        t_mean.SetTextAngle(90)
        text_size = 0.04
        t_mean.SetTextSize(text_size)
        rms_text.Draw()

        # Add CMS text
        cms_latex = ROOT.TLatex()
        cms_latex.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        cms_latex.SetTextFont(42)
        cms_latex.SetTextSize(0.04)
        # Get the text sitting just above the axes of the mean plot
        # Axes end inside the mean pad at (1-top_margin), but this has
        # to be scaled to canvas NDC
        # Then add a little extra spacing ontop to separate text from axes line
        latex_height = 1 - pad_offset_top - (mean_pads[0].GetAbsHNDC() * mean_pads[0].GetTopMargin()) + 0.02

        # Want it to start at the left edge of the first plot
        start_x = left_margin + (pad_width*pad_left_margin)
        # # start_x = 100
        if self.is_preliminary:
            if self.has_data:
                cms_latex.DrawLatexNDC(start_x, latex_height, "#font[62]{CMS}#font[52]{ Preliminary}")
            else:
                cms_latex.DrawLatexNDC(start_x, latex_height, "#font[62]{CMS}#font[52]{ Preliminary Simulation}")
        else:
            if self.has_data:
                cms_latex.DrawLatexNDC(start_x, latex_height, "#font[62]{CMS}")
            else:
                cms_latex.DrawLatexNDC(start_x, latex_height, "#font[62]{CMS}#font[52]{ Simulation}")
        cms_latex.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignBottom)
        # # Get the lumi text aligned to right edge of axes
        # # i.e. 1-pad_right_margin, but remember to scale by pad width
        end_x = 1 - right_margin - (mean_pads[0].GetAbsWNDC() * pad_right_margin)
        end_x = 0.985  # to match legend
        cms_latex.DrawLatexNDC(end_x, latex_height, " 35.9 fb^{-1} (13 TeV)")
        gc_stash.append(cms_latex)

        canvas.Update()
        canvas.SaveAs(output_file)


    def construct_delta_hist_groups(self, selections):
        """Construct delta hists for plots (nominal MC, & alt MC)

        See plot_delta_bins_summary() docstring about `selections` arg.

        Returns a list. Each item corresponds to the respective item in the
        outer level of the `selections` arg. (e.g. each angle).

        Each of the hist list items is itself a list, corresponding to the hists
        for [nominal MC, alt MC].

        Also does styling of hists (as easier here)
        """
        delta_hists = []
        for selection_group in selections:

            entries_nominal = []
            entries_alt = []

            bin_names = []

            for query, label in selection_group['selections']:
                # print(query)
                results = self.df.query(query)
                if len(results.index) != 1:
                    print(query)
                    print(results)
                    raise ValueError("Got != 1 results, check query")

                # if any([np.isnan(v) for v in [results['delta_truth'], results['delta_err_truth']]]):
                #     raise ValueError("Got a nan in delta_truth, delta_err_truth")
                # if any([np.isnan(v) for v in [results['delta_alt_truth'], results['delta_err_alt_truth']]]):
                #     raise ValueError("Got a nan in delta_truth, delta_err_truth")
                entries_nominal.append([results['delta_truth'], results['delta_err_truth']])
                entries_alt.append([results['delta_alt_truth'], results['delta_err_alt_truth']])

                bin_names.append(label)

            hist_delta_nominal = self._make_hist_from_values(entries_nominal, bin_names=bin_names)
            self._style_mc_hist(hist_delta_nominal)

            hist_delta_alt = self._make_hist_from_values(entries_alt, bin_names=bin_names)
            self._style_alt_mc_hist(hist_delta_alt)

            delta_hists.append([hist_delta_nominal, hist_delta_alt])

        return delta_hists


    def plot_delta_bins_summary(self, selections, output_file):
        """"""
        # Get data for plots
        delta_hists = self.construct_delta_hist_groups(selections)

        gc_stash = [] # to stop stuff being deleted

        # Setup canvas and pads
        # Setup same dimension etc as mean & RMS plot?
        canvas = ROOT.TCanvas("c_"+cu.get_unique_str(), "", 1200, 400)
        canvas.SetBottomMargin(0.0)
        canvas.SetTopMargin(0.0)
        canvas.SetLeftMargin(0.0)
        canvas.SetRightMargin(0.0)
        canvas.cd()

        # For each selection group, create a pad
        pads = []
        n_pads = len(selections)

        # gap between right end of plots and edge of canvas, used for legend
        right_margin = 0.12
        # pad_left_titles_gap = 0.01 # gap between pad_left_titles and all plots
        pad_to_pad_gap = 0.005  # gap between plot pad columns
        # how far in from the left the first plotting pad starts. used for y axis title
        left_margin = 0.04
        # figure out width per pad - get total width available, then divide by number of pads
        pad_width = (1 - left_margin - right_margin - (pad_to_pad_gap*(n_pads - 1))) / n_pads

        pad_offset_bottom = 0.01 # spacing between bottom of RMS pad and canvas edge
        pad_offset_top = 0.08  # spacing between top of mean pad and canvas edge - for CMS and lumi text

        # per-pad margins: these determine where the hist axes lie,
        # and are fractions of the **pad** width/height, not the global canvas
        pad_right_margin = 0.02
        pad_left_margin = 0.25

        # bottom margin includes space for x axis labels
        pad_bottom_margin = 0.49

        # extra bit to add to the top margins of lower and upper pads
        # to ensure y axis numbers don't get cut off
        pad_top_margin = 0.012

        for isel, selection_group in enumerate(selections):
            canvas.cd()
            pad_start_x = left_margin + (isel*pad_to_pad_gap) + (isel*pad_width)
            pad_end_x = pad_start_x + pad_width
            # Create pad
            pad = ROOT.TPad(cu.get_unique_str(), "", pad_start_x, pad_offset_bottom, pad_end_x, 1-pad_offset_top)
            ROOT.SetOwnership(pad, False)
            pad.SetFillColor(isel+2)
            pad.SetFillStyle(3004)
            pad.SetFillStyle(4000)
            pad.SetTopMargin(pad_top_margin)
            pad.SetBottomMargin(pad_bottom_margin)
            pad.SetRightMargin(pad_right_margin)
            pad.SetLeftMargin(pad_left_margin)
            pad.SetTicks(1, 1)
            pad.Draw()
            pads.append(pad)

        # Now draw the histograms
        for isel, (selection_group, pad, hist_group) in enumerate(zip(selections, pads, delta_hists)):
            pad.cd()
            hist_group[0].Draw("E1")
            hist_group[1].Draw("E1 SAME")

            draw_hist = hist_group[0]

            xax = draw_hist.GetXaxis()
            xax.CenterLabels()
            xax.LabelsOption("v")

            factor = 1.5
            label_size_fudge = 1.6
            label_offset_fudge = 4
            tick_fudge = 4

            xax.SetTickLength(xax.GetTickLength()*factor)
            xax.SetLabelSize(xax.GetLabelSize()*factor*2)
            xax.SetLabelOffset(xax.GetLabelOffset()*factor)

            yax = draw_hist.GetYaxis()
            yax.SetLabelSize(yax.GetLabelSize()*factor*label_size_fudge)
            yax.SetLabelOffset(yax.GetLabelOffset()*factor*label_offset_fudge)
            yax.SetTickLength(yax.GetTickLength()*factor*tick_fudge)  # extra bit of fudging, probably becasue pads are sligthly different sizes
            n_divisions = 1005  # fewer big ticks so less chance of overlapping with lower plot
            n_divisions = 510
            yax.SetNdivisions(n_divisions)

            # Set range using bin contents + error bars
            y_up, y_down = self.calc_hists_max_min(hist_group)
            y_range = y_up - y_down
            down_padding = 0.2 * y_range
            up_padding = 0.22 * y_range
            draw_hist.SetMinimum(y_down - down_padding)
            draw_hist.SetMaximum(y_up + up_padding)

            # Draw variable name
            var_latex = ROOT.TLatex()
            var_latex.SetTextAlign(ROOT.kHAlignCenter + ROOT.kVAlignBottom)
            var_latex.SetTextSize(0.1)
            # var_latex.SetTextFont(42)
            # these are relative to the RMS pad! not the canvas
            # put in middle of the plot (if 0.5, woud look off-centre)
            var_x = 0.5*(1-pad_right_margin-pad_left_margin) + pad_left_margin
            var_y = 0.03  # do by eye
            var_latex.DrawLatexNDC(var_x, var_y, selection_group['label'])
            gc_stash.append(var_latex)

        canvas.cd()

        # Add legend
        # Bit of a hack here - to get the fill hashing style the same as in the plots,
        # we have to put it in a pad of the same size as a plotting pad
        # This is because the hashing separation is scaled by the pad size
        # So we can't just put the legend in the global canvas.
        # We also can't modify the fill style of the legend entries (I tried, it does nothing)
        leg_y_top = 0.93
        leg_x_right = 1-0.0
        leg_pad = ROOT.TPad("leg_pad_"+cu.get_unique_str(), "", leg_x_right-pads[0].GetAbsWNDC(), leg_y_top-pads[0].GetAbsHNDC(), leg_x_right, leg_y_top)
        ROOT.SetOwnership(leg_pad, False)  # important! otherwise seg fault
        # leg_pad.SetFillColor(ROOT.kYellow)
        # leg_pad.SetFillStyle(3004)
        leg_pad.SetFillStyle(4000)
        leg_pad.SetLeftMargin(0)
        leg_pad.SetRightMargin(0)
        leg_pad.SetTopMargin(0)
        leg_pad.SetBottomMargin(0)
        leg_pad.Draw()
        leg_pad.cd()
        gc_stash.append(leg_pad)
        leg = ROOT.TLegend(0.28, pads[0].GetBottomMargin(), 1, 1)
        le_mc = leg.AddEntry(delta_hists[0][0], self.mc_label ,"L")
        le_mc_alt = leg.AddEntry(delta_hists[0][1], self.alt_mc_label ,"L")
        leg.SetTextSize(0.08)
        leg.Draw()
        canvas.cd()

        # Add delta text
        text_width = left_margin
        text_x = (0.5 * left_margin) - (0.5*text_width)
        # let ROOT center it, just make the box the height of the axis
        text_y_end = 1 - pad_offset_top - (pad_top_margin*pads[0].GetAbsHNDC())
        text_y = 1 - pad_offset_top - pads[0].GetAbsHNDC() + (pads[0].GetAbsHNDC()*pad_bottom_margin)
        delta_text = ROOT.TPaveText(text_x, text_y, text_x+text_width, text_y_end, "NDC NB")
        delta_text.SetFillStyle(0)
        delta_text.SetBorderSize(0)
        t_delta = delta_text.AddText("#Delta")
        t_delta.SetTextAngle(90)
        text_size = 0.06
        t_delta.SetTextSize(text_size)
        delta_text.Draw()

        # Add CMS text
        cms_latex = ROOT.TLatex()
        cms_latex.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        cms_latex.SetTextFont(42)
        cms_latex.SetTextSize(0.05)
        # Get the text sitting just above the axes of the mean plot
        # Axes end inside the mean pad at (1-top_margin), but this has
        # to be scaled to canvas NDC
        # Then add a little extra spacing ontop to separate text from axes line
        latex_height = 1 - pad_offset_top - (pads[0].GetAbsHNDC() * pads[0].GetTopMargin()) + 0.02

        # Want it to start at the left edge of the first plot
        start_x = left_margin + (pad_width*pad_left_margin)
        # # start_x = 100
        if self.is_preliminary:
            if self.has_data:
                cms_latex.DrawLatexNDC(start_x, latex_height, "#font[62]{CMS}#font[52]{ Preliminary}")
            else:
                cms_latex.DrawLatexNDC(start_x, latex_height, "#font[62]{CMS}#font[52]{ Preliminary Simulation}")
        else:
            if self.has_data:
                cms_latex.DrawLatexNDC(start_x, latex_height, "#font[62]{CMS}")
            else:
                cms_latex.DrawLatexNDC(start_x, latex_height, "#font[62]{CMS}#font[52]{ Simulation}")
        cms_latex.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignBottom)
        # Get the lumi text aligned to right edge of axes
        # i.e. 1-pad_right_margin, but remember to scale by pad width
        end_x = 1 - right_margin - (pads[0].GetAbsWNDC() * pad_right_margin)
        end_x = 0.985  # to match legend
        cms_latex.DrawLatexNDC(end_x, latex_height, " 35.9 fb^{-1} (13 TeV)")
        gc_stash.append(cms_latex)

        canvas.Update()
        canvas.SaveAs(output_file)

# ------------------------------------------------------------------------------
# VARIOUS CALCULATION METHODS
# ------------------------------------------------------------------------------

def get_bin_widths(hist):
    return hist.edges[1:] - hist.edges[:-1]


def hist_to_arrays(hist):
    """Convert histogram bins to arrays

    Assumes `hist` is an uproot.TH1 object, not a ROOT.TH1

    Note that this returns bin *areas* and not heights: assume bin contents
    already divided by bin width.

    Note that errors here are multiplied by the bin width (since they are
    assumed to have been created by originally dividing by the bin width)
    """
    bin_widths = get_bin_widths(hist)
    bin_areas = hist.values * bin_widths
    bin_centers = hist.edges[:-1] + (0.5*bin_widths)
    bin_errors = np.sqrt(hist.variances) * bin_widths
    return bin_areas, bin_widths, bin_centers, bin_errors


def scale_ematrix_by_bin_widths(ematrix, widths):
    this_widths = widths.reshape(len(widths), 1)
    return ematrix * this_widths * this_widths.T

# --------------------------------
# Functions to calculate mean
# --------------------------------


def check_hist_for_negatives(hist):
    areas, widths, centers, errors = hist_to_arrays(hist)
    for i, x in enumerate(areas, 1):
        if x < 0:
            raise ValueError("Area of bin %d = %f" % (i, x))


def calc_hist_mean(bin_areas, bin_centers):
    """Calculate mean of hist from value arrays.

    Must use np.X functions for e.g. sum(), square(), to ensure jax can differentiate it
    """
    return np.sum(bin_areas * bin_centers) / np.sum(bin_areas)


mean_differential = jit(grad(calc_hist_mean, argnums=0))


def calc_hist_mean_uncorrelated_error(bin_areas, bin_centers, bin_errors):
    """Calculate error on mean, assuming uncorrelated errors.

    Uses propagation of uncertainty bia partial differentials,
    calculated automatically using jax.
    """
    # differential wrt bin_areas
    diffs = mean_differential(bin_areas, bin_centers)
    err_sq = np.sum(np.square((diffs * bin_errors)))
    return np.sqrt(err_sq)


def calc_hist_mean_correlated_error(bin_areas, bin_centers, error_matrix):
    """Get error on mean, assuming covariance matrix error_matrix"""
    diffs = mean_differential(bin_areas, bin_centers)
    sum_sq = diffs @ error_matrix @ diffs
    return np.sqrt(sum_sq)


def calc_hist_mean_and_uncorrelated_error(hist):
    """Calculate from hist both mean and its error,
    assuming uncorrelated uncertainties

    Parameters
    ----------
    hist : TH1

    Returns
    -------
    float, float
    """
    areas, widths, centers, errors = hist_to_arrays(hist)
    mean = calc_hist_mean(areas, centers)
    err = calc_hist_mean_uncorrelated_error(areas, centers, errors)
    return float(mean), float(err)


def calc_hist_mean_and_correlated_error(hist, ematrix):
    areas, widths, centers, errors = hist_to_arrays(hist)
    mean = calc_hist_mean(areas, centers)
    err = calc_hist_mean_correlated_error(areas, centers, ematrix)
    return float(mean), float(err)

# --------------------------------
# Functions to calculate RMS
# --------------------------------

def calc_hist_rms(bin_areas, bin_centers):
    """Calculate RMS of hist from value arrays.

    Must use np.X functions for e.g. sum(), square(), to ensure jax can differentiate it
    """
    mean = calc_hist_mean(bin_areas, bin_centers)
    sum_sq = np.sum(np.square((bin_areas * bin_centers) - mean))
    return np.sqrt(sum_sq / np.sum(bin_areas))


rms_differential = jit(grad(calc_hist_rms, argnums=0))


def calc_hist_rms_uncorrelated_error(bin_areas, bin_centers, bin_errors):
    """Calculate error on RMS, assuming uncorrelated errors.

    Uses propagation of uncertainty bia partial differentials,
    calculated automatically using jax.
    """
    # differential wrt bin_areas
    diffs = rms_differential(bin_areas, bin_centers)
    err_sq = np.sum(np.square((diffs * bin_errors)))
    return np.sqrt(err_sq)


def calc_hist_rms_correlated_error(bin_areas, bin_centers, error_matrix):
    """Get error on rms, assuming covariance matrix error_matrix"""
    diffs = rms_differential(bin_areas, bin_centers)
    sum_sq = diffs @ error_matrix @ diffs
    return np.sqrt(sum_sq)


def calc_hist_rms_and_uncorrelated_error(hist):
    """Calculate from hist both RMS and its error,
    assuming uncorrelated uncertainties

    Parameters
    ----------
    hist : TH1

    Returns
    -------
    float, float
    """
    areas, widths, centers, errors = hist_to_arrays(hist)
    rms = calc_hist_rms(areas, centers)
    err = calc_hist_rms_uncorrelated_error(areas, centers, errors)
    return float(rms), float(err)


def calc_hist_rms_and_correlated_error(hist, ematrix):
    areas, widths, centers, errors = hist_to_arrays(hist)
    rms = calc_hist_rms(areas, centers)
    err = calc_hist_rms_correlated_error(areas, centers, ematrix)
    return float(rms), float(err)

# --------------------------------
# Functions to calculate delta
# --------------------------------

def calc_hist_delta_and_error(hist_a, ematrix_a, hist_b):
    """Calculate delta between hists, along with its error

    Defined as 0.5 * integral[ (a - b)^2 / (a+b) ]
    """
    areas_a, widths_a, centers_a, errors_a = hist_to_arrays(hist_a)
    areas_b, widths_b, centers_b, errors_b = hist_to_arrays(hist_b)
    delta = calc_hist_delta(areas_a, areas_b)
    err = calc_hist_delta_correlated_error(areas_a, ematrix_a, areas_b, errors_b)
    return float(delta), float(err)


def calc_hist_delta(areas_a, areas_b):
    # do I need bin areas or densities?
    # I guess since by definition sum(area_a) = 1, areas are needed?!
    integrand = np.true_divide(np.square(areas_a - areas_b), areas_a + areas_b)
    # nan_to_num important as divide gives nans if both 0
    delta = 0.5 * np.sum(np.nan_to_num(integrand))
    return delta


delta_diff = jit(grad(calc_hist_delta, argnums=[0, 1]))

def calc_hist_delta_uncorrelated_error(areas_a, errors_a, areas_b, errors_b):
    pass


def calc_hist_delta_correlated_error(areas_a, ematrix_a, areas_b, errors_b):
    diffs_a, diffs_b = delta_diff(areas_a, areas_b)
    # need to do nan_to_num since the differential can return nan...
    # not sure how to fix "properly" though
    diffs_a = np.nan_to_num(diffs_a)
    diffs_b = np.nan_to_num(diffs_b)
    # for the total, we need to do
    # diffs_a * ematrix_a * diffs_a + diffs_b*errors_b*diffs_b,
    # since the errors on a and b have no connections, we can get away with this.
    err_a_sq = diffs_a.T @ ematrix_a @ diffs_a
    err_b_sq = np.sum(np.square((diffs_b * errors_b)))
    return np.sqrt(err_a_sq + err_b_sq)


def unpack_slim_unfolding_root_file_uproot(input_tfile, region_name, angle_name, pt_bins):
    tdir = "%s/%s" % (region_name, angle_name)
    indices = range(len(pt_bins)-1)

    unfolding_stat_err_hists = [
        input_tfile["%s/unfolded_stat_err_norm_divBinWidth_%d" % (tdir, ibin)]
        for ibin in indices
    ]
    unfolding_total_err_hists = [
        input_tfile["%s/unfolded_norm_divBinWidth_%d" % (tdir, ibin)]
        for ibin in indices
    ]
    unfolding_total_err_ematrices = [
        input_tfile["%s/unfolded_total_ematrix_norm_divBinWidth_%d" % (tdir, ibin)]
        for ibin in indices
    ]
    truth_hists = [
        input_tfile["%s/hist_truth_norm_divBinWidth_%d" % (tdir, ibin)]
        for ibin in indices
    ]
    alt_truth_hists = [
        input_tfile["%s/alt_hist_truth_norm_divBinWidth_%d" % (tdir, ibin)]
        for ibin in indices
    ]

    return dict(
        unfolding_stat_err_hists=unfolding_stat_err_hists,
        unfolding_total_err_hists=unfolding_total_err_hists,
        unfolding_total_ematrices=unfolding_total_err_ematrices,
        truth_hists=truth_hists,
        alt_truth_hists=alt_truth_hists,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ak4source",
                        help="Source directory for AK4 jets (should be the one made by unfolding.py")
    parser.add_argument("--ak8source",
                        help="Source directory for AK8 jets (should be the one made by unfolding.py")
    parser.add_argument("--h5input",
                        help="Read data from H5 input file (from previous running of this script)")
    parser.add_argument("--h5output",
                        default=None,
                        help=("Store data as H5 output file (ignored if --h5input used). "
                              "Default if <outputDir>/store.h5"))
    parser.add_argument("--outputDir",
                        default=None,
                        help='Output directory (default is the source dir')
    args = parser.parse_args()

    # Get data
    if not any([args.h5input, args.ak4source, args.ak8source]):
        raise RuntimeError("Need one of --h5input or --ak4input/--ak8input")

    if args.h5input is None:
        # ----------------------------------------------------------------------
        # READ IN DATA FROM UNFOLDING ROOT FILES
        # ----------------------------------------------------------------------
        if not args.outputDir and args.ak4source:
            args.outputDir = os.path.join(args.ak4source, 'SummaryPlots')
        elif not args.outputDir and args.ak8source:
            args.outputDir = os.path.join(args.ak8source, 'SummaryPlots')

        cu.check_dir_exists_create(args.outputDir)

        results_dicts = []

        jet_algos = []
        if args.ak4source:
            jet_algos.append({'src': args.ak4source, 'label': 'AK4 PUPPI', 'name': 'ak4puppi'})
        if args.ak8source:
            jet_algos.append({'src': args.ak8source, 'label': 'AK8 PUPPI', 'name': 'ak8puppi'})

        for jet_algo in jet_algos:
            print("--- Jet algo ---:", jet_algo['label'])
            source_dir = jet_algo['src']
            regions = [
                get_dijet_config(source_dir, central=True, groomed=False),
                get_dijet_config(source_dir, central=False, groomed=False),
                get_dijet_config(source_dir, central=True, groomed=True),
                get_dijet_config(source_dir, central=False, groomed=True),
                get_zpj_config(source_dir, groomed=False),
                get_zpj_config(source_dir, groomed=True),
            ]

            for region in regions:
                region_dir = os.path.join(source_dir, region['name'])
                if not os.path.isdir(region_dir):
                    print("! Warning ! cannot find region dir", region_dir, '- skipping region')
                    continue

                angles = qgc.COMMON_VARS[:]
                for angle in angles:
                    angle_output_dir = "%s/%s" % (region_dir, angle.var)
                    if not os.path.isdir(angle_output_dir):
                        print("! Warning ! cannot find angle dir", angle_output_dir, '- skipping angle', angle.var)
                        continue

                    this_region = copy(region)
                    # Get region dict from pickle file
                    # pickle_filename = os.path.join(angle_output_dir, "unfolding_result.pkl")
                    # unpickled_region = unpickle_region(pickle_filename)

                    # # # check
                    # if this_region['name'] != unpickled_region['name']:
                    #     raise RuntimeError("Mismatch region name")

                    # this_region.update(unpickled_region)

                    # Get bare necessary hists from slim ROOT file
                    # Using pickle one is much slower
                    root_filename = os.path.join(angle_output_dir, "unfolding_result_slim.root")
                    pt_bins = qgc.PT_UNFOLD_DICT['signal_zpj_gen'] if 'ZPlusJets' in this_region['name'] else qgc.PT_UNFOLD_DICT['signal_gen']
                    uproot_file = uproot.open(root_filename)
                    unfolding_dict = unpack_slim_unfolding_root_file_uproot(uproot_file, this_region['name'], angle.var, pt_bins)

                    # common str to put on filenames, etc.
                    # don't need angle_prepend as 'groomed' in region name
                    append = "%s_%s" % (this_region['name'], angle.var)
                    print("*"*120)
                    print("Region/var: %s" % (append))
                    print("*"*120)

                    # ----------------------------------------------------------
                    # CALCULATE STATS FOR EACH PT BIN
                    # ----------------------------------------------------------
                    # Iterate through pt bins, get lambda histogram for that bin,
                    # derive metrics from it, save
                    for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(pt_bins[:-1], pt_bins[1:])):
                        # print("   done pt bin", ibin)

                        # Handle nominal MC hist -> metrics
                        mc_gen_hist_bin = unfolding_dict['truth_hists'][ibin]
                        try:
                            check_hist_for_negatives(mc_gen_hist_bin)
                        except ValueError as e:
                            print("-ve value for MC hist in pt bin", ibin, ":", bin_edge_low, "-", bin_edge_high)
                            raise e
                        mc_gen_hist_bin_mean, mc_gen_hist_bin_mean_err = calc_hist_mean_and_uncorrelated_error(mc_gen_hist_bin)
                        mc_gen_hist_bin_rms, mc_gen_hist_bin_rms_err = calc_hist_rms_and_uncorrelated_error(mc_gen_hist_bin)

                        # Handle alt MC hist -> metrics
                        alt_mc_gen_hist_bin = unfolding_dict['alt_truth_hists'][ibin]
                        try:
                            check_hist_for_negatives(alt_mc_gen_hist_bin)
                        except ValueError as e:
                            print("-ve value for alt MC hist in pt bin", ibin, ":", bin_edge_low, "-", bin_edge_high)
                            raise e
                        alt_mc_gen_hist_bin_mean, alt_mc_gen_hist_bin_mean_err = calc_hist_mean_and_uncorrelated_error(alt_mc_gen_hist_bin)
                        alt_mc_gen_hist_bin_rms, alt_mc_gen_hist_bin_rms_err = calc_hist_rms_and_uncorrelated_error(alt_mc_gen_hist_bin)

                        # Handle unfolded data hist -> metrics
                        unfolded_hist_bin_total_errors = unfolding_dict['unfolding_total_err_hists'][ibin]
                        try:
                            check_hist_for_negatives(unfolded_hist_bin_total_errors)
                        except ValueError as e:
                            print("-ve value for data hist in pt bin", ibin, ":", bin_edge_low, "-", bin_edge_high)
                            raise e
                        ematrix = scale_ematrix_by_bin_widths(unfolding_dict['unfolding_total_ematrices'][ibin].values, get_bin_widths(unfolded_hist_bin_total_errors))
                        unfolded_hist_bin_total_errors_rms, unfolded_hist_bin_total_errors_rms_err = calc_hist_rms_and_correlated_error(unfolded_hist_bin_total_errors, ematrix)
                        unfolded_hist_bin_total_errors_mean, unfolded_hist_bin_total_errors_mean_err = calc_hist_mean_and_correlated_error(unfolded_hist_bin_total_errors, ematrix)

                        delta_nominal, delta_nominal_err = calc_hist_delta_and_error(unfolded_hist_bin_total_errors, ematrix, mc_gen_hist_bin)
                        delta_alt, delta_alt_err = calc_hist_delta_and_error(unfolded_hist_bin_total_errors, ematrix, alt_mc_gen_hist_bin)

                        result_dict = {
                            'jet_algo': jet_algo['name'],
                            'region': this_region['name'], # TODO remove "_groomed"?
                            'isgroomed': 'groomed' in this_region['name'].lower(),
                            'pt_bin': ibin,
                            'angle': angle.var,

                            'mean': unfolded_hist_bin_total_errors_mean,
                            'mean_err': unfolded_hist_bin_total_errors_mean_err, # FIXME

                            'mean_truth': mc_gen_hist_bin_mean,
                            'mean_err_truth': mc_gen_hist_bin_mean_err,

                            'mean_alt_truth': alt_mc_gen_hist_bin_mean,
                            'mean_err_alt_truth': alt_mc_gen_hist_bin_mean_err,

                            'rms': unfolded_hist_bin_total_errors_rms,
                            'rms_err': unfolded_hist_bin_total_errors_rms_err, #FIXME

                            'rms_truth': mc_gen_hist_bin_rms,
                            'rms_err_truth': mc_gen_hist_bin_rms_err,

                            'rms_alt_truth': alt_mc_gen_hist_bin_rms,
                            'rms_err_alt_truth': alt_mc_gen_hist_bin_rms_err,

                            'delta_truth': delta_nominal,
                            'delta_err_truth': delta_nominal_err,

                            'delta_alt_truth': delta_alt,
                            'delta_err_alt_truth': delta_alt_err,
                        }
                        results_dicts.append(result_dict)

                    # important to keep memory footprint small
                    # del unpickled_region
                    del this_region

        df = pd.DataFrame(results_dicts)
        df['jet_algo'] = df['jet_algo'].astype('category')
        df['region'] = df['region'].astype('category')
        df['angle'] = df['angle'].astype('category')
        print(df.head())
        print(df.tail())
        print(len(df.index), 'entries in dataframe')
        print(df.dtypes)

        if args.h5output is None:
            args.h5output = os.path.join(args.outputDir, "store.h5")
        print("Saving dataframe to", args.h5output)
        # need format='table' to store category dtype
        df.to_hdf(args.h5output, key='df', format='table')

    else:
        # ----------------------------------------------------------------------
        # READ IN DATA FROM H5 FILE
        # -----------------------------------------------------------------------
        print("Reading in from data existing HDF5 file...")
        if not args.outputDir:
            args.outputDir = os.path.dirname(os.path.abspath(args.h5input))

        with pd.HDFStore(args.h5input) as store:
            df = store['df']
        print(df.head())
        print("# entries:", len(df.index))

    # Filter only regions/algos/angles in the dataframe, since it could have
    # been modified earlier
    all_jet_algos = [
        {'src': args.ak4source, 'label': 'AK4 PUPPI', 'name': 'ak4puppi'},
        {'src': args.ak8source, 'label': 'AK8 PUPPI', 'name': 'ak8puppi'}
    ]
    jet_algos = [j for j in all_jet_algos if j['name'] in df['jet_algo'].unique()]
    # print("Plotting jet_algos:", jet_algos)

    all_regions = [
        get_dijet_config('', central=True, groomed=False),
        get_dijet_config('', central=False, groomed=False),
        get_dijet_config('', central=True, groomed=True),
        get_dijet_config('', central=False, groomed=True),
        get_zpj_config('', groomed=False),
        get_zpj_config('', groomed=True),
    ]
    regions = [r for r in all_regions if r['name'] in df['region'].unique()]
    # print("Plotting regions:", regions)

    angles = [a for a in qgc.COMMON_VARS if a.var in df['angle'].unique()]
    # print("Plotting angles:", angles)
    charged_only_angles = [a for a in angles if "charged" in a.var]
    charged_and_neutral_angles = [a for a in angles if "charged" not in a.var]
    # charged_and_neutral_angles[0] = charged_and_neutral_angles[1]

    # --------------------------------------------------------------------------
    # Do all the plotting
    # --------------------------------------------------------------------------
    print("*"*120)
    print("* Plotting time!")
    print("*"*120)
    plotter = SummaryPlotter(jet_algos,
                             regions,
                             angles,
                             qgc.PT_UNFOLD_DICT['signal_gen'],
                             qgc.PT_UNFOLD_DICT['signal_zpj_gen'],
                             df,
                             args.outputDir,
                             has_data=True)

    has_dijet = any(["Dijet" in r['name'] for r in regions])
    has_zpj = any(["ZPlusJet" in r['name'] for r in regions])

    # For summary plots
    low_pt = 120
    high_pt = 614
    # high_pt = 800

    ak4_str = "AK4"
    ak8_str = "AK8"

    if has_dijet:
        plotter.plot_dijet_means_vs_pt_all()
        plotter.plot_dijet_rms_vs_pt_all()
        plotter.plot_dijet_delta_vs_pt_all()

        pt_bins = qgc.PT_UNFOLD_DICT['signal_gen']
        low_pt_bin = np.where(pt_bins == low_pt)[0][0]
        low_pt_str = "[%g, %g] GeV" % (low_pt, pt_bins[low_pt_bin+1])

        high_pt_bin = np.where(pt_bins == high_pt)[0][0]
        high_pt_str = "[%g, %g] GeV" % (high_pt, pt_bins[high_pt_bin+1])

        # Create selection queries for mean/RMS/delta summary plots
        selections = []
        for angle in charged_and_neutral_angles:
            this_angle_str = "%s (%s)" % (angle.name, angle.lambda_str)
            # this_angle_str = "%s" % (angle.lambda_str)
            this_selection = [
                ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_central" & angle=="%s"' % (low_pt_bin, angle.var),
                    "{jet_str}, {pt_str}".format(jet_str=ak4_str, pt_str=low_pt_str)),
                    # "#splitline{{{jet_str}}}{{{pt_str}}}".format(jet_str=ak4_str, pt_str=low_pt_str)),

                ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_central" & angle=="%s"' % (high_pt_bin, angle.var),
                    "{jet_str}, {pt_str}".format(jet_str=ak4_str, pt_str=high_pt_str)),
                    # "#splitline{{{jet_str}}}{{{pt_str}}}".format(jet_str=ak4_str, pt_str=high_pt_str)),

                ('jet_algo=="ak8puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_central" & angle=="%s"' % (low_pt_bin, angle.var),
                    "{jet_str}, {pt_str}".format(jet_str=ak8_str, pt_str=low_pt_str)),
                    # "#splitline{{{jet_str}}}{{{pt_str}}}".format(jet_str=ak8_str, pt_str=low_pt_str)),

                ('jet_algo=="ak8puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_central" & angle=="%s_charged"' % (low_pt_bin, angle.var),
                    "#splitline{{{jet_str}, {pt_str}}}{{Charged-only}}".format(jet_str=ak8_str, pt_str=low_pt_str)),
                    # "#splitline{{#splitline{{{jet_str}}}{{{pt_str}}}}}{{Charged-only}}".format(jet_str=ak8_str, pt_str=low_pt_str)),

                ('jet_algo=="ak4puppi" & pt_bin==%d & isgroomed  & region=="Dijet_central_groomed" & angle=="%s"' % (low_pt_bin, angle.var),
                    "#splitline{{{jet_str}, {pt_str}}}{{Groomed}}".format(jet_str=ak4_str, pt_str=low_pt_str)),
            ]
            selections.append({'label': this_angle_str, 'selections': this_selection})

        plotter.plot_mean_rms_bins_summary(
            selections=selections,
            output_file=os.path.join(args.outputDir, "dijet_mean_rms_summary.pdf")
        )

        plotter.plot_delta_bins_summary(
            selections=selections,
            output_file=os.path.join(args.outputDir, "dijet_delta_summary.pdf")
        )

    if has_zpj:
        plotter.plot_zpj_means_vs_pt_all()
        plotter.plot_zpj_rms_vs_pt_all()
        plotter.plot_zpj_delta_vs_pt_all()

        pt_bins = qgc.PT_UNFOLD_DICT['signal_zpj_gen']
        low_pt_bin = np.where(pt_bins == low_pt)[0][0]
        low_pt_str = "[%g, %g] GeV" % (low_pt, pt_bins[low_pt_bin+1])

        high_pt_bin = np.where(pt_bins == high_pt)[0][0]
        high_pt_str = "[%g, %g] GeV" % (high_pt, pt_bins[high_pt_bin+1])

        selections = []
        for angle in charged_and_neutral_angles:
            this_angle_str = "%s (%s)" % (angle.name, angle.lambda_str)
            # this_angle_str = "%s" % (angle.lambda_str)
            this_selection = [
                ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s"' % (low_pt_bin, angle.var),
                    "{jet_str}, {pt_str}".format(jet_str=ak4_str, pt_str=low_pt_str)),
                    # "#splitline{{{jet_str}}}{{{pt_str}}}".format(jet_str=ak4_str, pt_str=low_pt_str)),

                ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s"' % (high_pt_bin, angle.var),
                    "{jet_str}, {pt_str}".format(jet_str=ak4_str, pt_str=high_pt_str)),
                    # "#splitline{{{jet_str}}}{{{pt_str}}}".format(jet_str=ak4_str, pt_str=high_pt_str)),

                ('jet_algo=="ak8puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s"' % (low_pt_bin, angle.var),
                    "{jet_str}, {pt_str}".format(jet_str=ak8_str, pt_str=low_pt_str)),
                    # "#splitline{{{jet_str}}}{{{pt_str}}}".format(jet_str=ak8_str, pt_str=low_pt_str)),

                ('jet_algo=="ak8puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s_charged"' % (low_pt_bin, angle.var),
                    "#splitline{{{jet_str}, {pt_str}}}{{Charged-only}}".format(jet_str=ak8_str, pt_str=low_pt_str)),
                    # "#splitline{{#splitline{{{jet_str}}}{{{pt_str}}}}}{{Charged-only}}".format(jet_str=ak8_str, pt_str=low_pt_str)),

                ('jet_algo=="ak4puppi" & pt_bin==%d & isgroomed  & region=="ZPlusJets_groomed" & angle=="%s"' % (low_pt_bin, angle.var),
                    "#splitline{{{jet_str}, {pt_str}}}{{Groomed}}".format(jet_str=ak4_str, pt_str=low_pt_str)),
            ]
            selections.append({'label': this_angle_str, 'selections': this_selection})

        plotter.plot_mean_rms_bins_summary(
            selections=selections,
            output_file=os.path.join(args.outputDir, "zpj_mean_rms_summary.pdf")
        )

        plotter.plot_delta_bins_summary(
            selections=selections,
            output_file=os.path.join(args.outputDir, "zpj_delta_summary.pdf")
        )

    # if has_dijet and has_zpj:
    #     plotter.plot_dijet_zpj_means_vs_pt_all()
        # plotter.plot_dijet_zpj_rms_vs_pt_all()

