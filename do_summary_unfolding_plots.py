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

        if "LHA" in angle.var:
            print(dijet_central_data)
            for i in range(1, dijet_central_hist.GetNbinsX()+1):
                print(i, dijet_central_hist.GetBinContent(i), dijet_central_hist.GetBinError(i))

        region_name = 'Dijet_forward'
        if do_groomed:
            region_name += "_groomed"
        dijet_forward_data = df[mask & (df['region'] == region_name)]
        dijet_forward_hist = self.data_to_hist(dijet_forward_data['mean'], dijet_forward_data['mean_err'], self.pt_bins_dijet)
        dijet_forward_hist_truth = self.data_to_hist(dijet_forward_data['mean_truth'], dijet_forward_data['mean_err_truth'], self.pt_bins_dijet)
        dijet_forward_hist_alt_truth = self.data_to_hist(dijet_forward_data['mean_alt_truth'], dijet_forward_data['mean_err_alt_truth'], self.pt_bins_dijet)

        if do_zpj:
            region_name = 'ZPlusJets'
            if do_groomed:
                region_name += "_groomed"
            # drop last pt bin as massive error
            zpj_data = df[mask & (df['region'] == region_name) & (df['pt_bin'] < (len(self.pt_bins_zpj)-3))]
            zpj_hist = self.data_to_hist(zpj_data['mean'], zpj_data['mean_err'], self.pt_bins_zpj[:-2])
            zpj_hist_truth = self.data_to_hist(zpj_data['mean_truth'], zpj_data['mean_err_truth'], self.pt_bins_zpj[:-2])
            zpj_hist_alt_truth = self.data_to_hist(zpj_data['mean_alt_truth'], zpj_data['mean_err_alt_truth'], self.pt_bins_zpj[:-2])

        dijet_cen_col = COMMON_STYLE_DICT['dijet_cen_color']
        dijet_fwd_col = COMMON_STYLE_DICT['dijet_fwd_color']
        zpj_col = COMMON_STYLE_DICT['zpj_color']
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
                             subplot=dijet_central_hist),
                Contribution(dijet_forward_hist_truth, label='#splitline{ Dijet (forward)  }{ [%s]}' % (self.mc_label),
                             line_color=dijet_fwd_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_line_style'],
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', False), marker_size=0,
                             subplot=dijet_forward_hist),
            ])
        if do_zpj:
            entries.extend([
                Contribution(zpj_hist_truth, label='#splitline{ Z+jet  }{ [%s]}' % (self.mc_label),
                             line_color=zpj_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_line_style'],
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', False), marker_size=0,
                             subplot=zpj_hist),
            ])
        # add alt MC
        if do_dijet:
            entries.extend([
                Contribution(dijet_central_hist_alt_truth, label='#splitline{ Dijet (central)  }{ [%s]}' % (self.alt_mc_label),
                             line_color=dijet_cen_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                             marker_color=dijet_cen_col, marker_style=cu.Marker.get('circle', False), marker_size=0,
                             subplot=dijet_central_hist),
                Contribution(dijet_forward_hist_alt_truth, label='#splitline{ Dijet (forward)  }{ [%s]}' % (self.alt_mc_label),
                             line_color=dijet_fwd_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', False), marker_size=0,
                             subplot=dijet_forward_hist),
            ])
        if do_zpj:
            entries.extend([
                Contribution(zpj_hist_alt_truth, label='#splitline{ Z+jet}{ [%s]}' % (self.alt_mc_label),
                             line_color=zpj_col, line_width=2, line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', False), marker_size=0,
                             subplot=zpj_hist),
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
        plot.plot("NOSTACK E1 X0", "NOSTACK E1")
        # plot.get_modifier().GetYaxis().SetTitleOffset(plot.get_modifier().GetYaxis().GetTitleOffset()*1.1)
        plot.set_logx(do_more_labels=False)
        groomed_str = '_groomed' if do_groomed else ''
        plot.save("%s/dijet_zpj_means_vs_pt_%s%s_%s.%s" % (output_dir, angle.var, groomed_str, jet_algo['name'], self.output_fmt))

    def plot_dijet_zpj_rms_vs_pt_all(self):
        """Plot mean vs pt for dijet (cen+fwd) & Z+jet on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_zpj_rms_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_rms_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_zpj_rms_vs_pt_all' % self.output_dir)

    def plot_dijet_rms_vs_pt_all(self):
        """Plot mean vs pt for dijet (cen+fwd) on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_rms_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_rms_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_rms_vs_pt_all' % self.output_dir, do_zpj=False)

    def plot_zpj_rms_vs_pt_all(self):
        """Plot mean vs pt for zpj on one plot,
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

        if do_zpj:
            region_name = 'ZPlusJets'
            if do_groomed:
                region_name += "_groomed"
            zpj_data = df[mask & (df['region'] == region_name) & (df['pt_bin'] < (len(self.pt_bins_zpj)-3))]
            zpj_hist = self.data_to_hist(zpj_data['rms'], zpj_data['rms_err'], self.pt_bins_zpj[:-2])
            zpj_hist_truth = self.data_to_hist(zpj_data['rms_truth'], zpj_data['rms_err_truth'], self.pt_bins_zpj[:-2])
            zpj_hist_alt_truth = self.data_to_hist(zpj_data['rms_alt_truth'], zpj_data['rms_err_alt_truth'], self.pt_bins_zpj[:-2])

        # unify this?
        dijet_cen_col = COMMON_STYLE_DICT['dijet_cen_color']
        dijet_fwd_col = COMMON_STYLE_DICT['dijet_fwd_color']
        zpj_col = COMMON_STYLE_DICT['zpj_color']
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
                             subplot=dijet_central_hist),
                Contribution(dijet_forward_hist_truth, label='#splitline{ Dijet (forward)  }{ [%s]}' % (self.mc_label),
                             line_color=dijet_fwd_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_line_style'],
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', False), marker_size=0,
                             subplot=dijet_forward_hist),
            ])
        if do_zpj:
            entries.extend([
                Contribution(zpj_hist_truth, label='#splitline{ Z+jet  }{ [%s]}' % (self.mc_label),
                             line_color=zpj_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_line_style'],
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', False), marker_size=0,
                             subplot=zpj_hist),
            ])
        # add alt MC
        if do_dijet:
            entries.extend([
                Contribution(dijet_central_hist_alt_truth, label='#splitline{ Dijet (central)  }{ [%s]}' % (self.alt_mc_label),
                             line_color=dijet_cen_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                             marker_color=dijet_cen_col, marker_style=cu.Marker.get('circle', False), marker_size=0,
                             subplot=dijet_central_hist),
                Contribution(dijet_forward_hist_alt_truth, label='#splitline{ Dijet (forward)  }{ [%s]}' % (self.alt_mc_label),
                             line_color=dijet_fwd_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', False), marker_size=0,
                             subplot=dijet_forward_hist),
            ])
        if do_zpj:
            entries.extend([
                Contribution(zpj_hist_alt_truth, label='#splitline{ Z+jet}{ [%s]}' % (self.alt_mc_label),
                             line_color=zpj_col, line_width=lw, line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', False), marker_size=0,
                             subplot=zpj_hist),
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
        plot.plot("NOSTACK E1 X0", "NOSTACK E1")
        # plot.get_modifier().GetYaxis().SetTitleOffset(plot.get_modifier().GetYaxis().GetTitleOffset()*1.1)
        plot.set_logx(do_more_labels=False)
        groomed_str = '_groomed' if do_groomed else ''
        plot.save("%s/dijet_zpj_rms_vs_pt_%s%s_%s.%s" % (output_dir, angle.var, groomed_str, jet_algo['name'], self.output_fmt))

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
            mean_pad.cd()  # this causes a segfault for no reason?

            mean_hist_group[2].Draw("E1")
            mean_hist_group[1].Draw("E1 SAME")
            mean_hist_group[0].Draw("E1 SAME")

            mean_draw_hist = mean_hist_group[2]
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

            rms_draw_hist = rms_hist_group[2]
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

# ------------------------------------------------------------------------------
# VARIOUS CALCULATION METHODS
# ------------------------------------------------------------------------------

def hist_to_arrays(hist):
    """Convert histogram bins to arrays

    Note that this returns bin *areas* and not heights: assume bin contents
    already divided by bin width.

    Note that errors here are multiplied by the bin width (since they are
    assumed to have been created by originally dividing by the bin width)
    """
    bin_areas = np.array([hist.GetBinContent(ibin)*hist.GetBinWidth(ibin) for ibin in range(1, hist.GetNbinsX()+1)])
    bin_centers = np.array([hist.GetBinCenter(ibin) for ibin in range(1, hist.GetNbinsX()+1)])
    bin_widths = np.array([hist.GetBinWidth(ibin) for ibin in range(1, hist.GetNbinsX()+1)])
    bin_errors = np.array([hist.GetBinError(ibin)*hist.GetBinWidth(ibin) for ibin in range(1, hist.GetNbinsX()+1)])
    return bin_areas, bin_widths, bin_centers, bin_errors


def calc_hist_mean(hist):
    """Manual calculation of mean of normalised histogram, taking into account bin width

    Note that this is exactly the formula for a weighted mean, where the bin areas are the weights
    (since already divided by bin width), and the bin centres are the "data"
    https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
    """
    bin_areas = np.array([hist.GetBinContent(ibin)*hist.GetBinWidth(ibin) for ibin in range(1, hist.GetNbinsX()+1)])
    bin_centers = np.array([hist.GetBinCenter(ibin) for ibin in range(1, hist.GetNbinsX()+1)])
    return float(calc_hist_mean_from_values(bin_areas, bin_centers))


def calc_hist_mean_uncorrelated_error(hist):
    """Get error on mean using error bars, assuming 0 correlation between errors
    i.e. statistical only.

    Assumes aready divided by bin width!

    if mean = sum_{bins} (bin_center * bin_area) / N
    where N = sum_{bins} (bin_area)
    then err^2 = sum_{bins} [ ( (bin_center * bin_error / N) - (mean/N) ) * bin_error * bin_width ]^2

    (if bin_error has already been divided by the bin width)

    since our error is on both the numerator and denominator

    **RECOMMEND USING calc_hist_mean_uncorrelated_error_from_values() INSTEAD**

    Parameters
    ----------
    hist : ROOT.TH1
        Histogram
    """
    mean = calc_hist_mean(hist)
    N = sum([hist.GetBinContent(ibin)*hist.GetBinWidth(ibin) for ibin in range(1, hist.GetNbinsX()+1)])
    err_sq = sum([pow( ((hist.GetBinCenter(ibin)/N) - (mean/N) )*hist.GetBinWidth(ibin)*hist.GetBinError(ibin), 2) for ibin in range(1, hist.GetNbinsX()+1) ])
    return sqrt(err_sq)


def calc_hist_mean_error(hist, covariance_matrix, scale_factor):
    """Get error on the mean using covariance matrix

    hist is a TH1
    covariance_matrix should be a 2D numpy array for this bin

    Uses the fact that if f = A+B,
    then sigma^2_f = sigma^2_A + sigma^2_B + 2*sigma_A*sigma_B

    And also if f = aA, where a const, then sigma^2_f = a^2 sigma^2_A
    Used for bin width scaling, and the fact that mean = sum (bin_i / Nbins)

    And if f = A / B, then
    sigma^2_f = f^2 ( (sigma_A / A)^2 + (sigma_B / B)^2 - 2 sigma_AB / (A*B) )


    """
    # f = hist.GetMean()
    # f2 = f**2

    sum_sq = 0
    sum_reduced = 0
    # print("Doing mean error")
    for ibin in range(1, hist.GetNbinsX()+1):
        # bin_h_i = hist.GetBinContent(ibin)
        bin_x_i = hist.GetBinCenter(ibin)

        # sum_reduced += (bin_h_i) # / hist.GetBinWidth(ibin))
        for jbin in range(ibin, hist.GetNbinsX()+1):
            # bin_h_j = hist.GetBinContent(jbin)
            bin_x_j = hist.GetBinCenter(jbin)
            # Account for bin width scaling, and normalisation factor
            # total_bin_width = hist.GetBinWidth(ibin) * hist.GetBinWidth(jbin) * scale_factor * scale_factor
            this_err2 = 0
            if ibin == jbin:
                this_err2 = bin_x_i*bin_x_j*covariance_matrix[ibin-1, ibin-1]
                # this_err2 = covariance_matrix[ibin-1, ibin-1]
            else:
                # pass
                this_err2 = 2*bin_x_i*bin_x_j*covariance_matrix[ibin-1, jbin-1]  # should this be sqrt?
                # this_err2 = 2*covariance_matrix[ibin-1, jbin-1]  # should this be sqrt?

            sum_sq += this_err2


    # scale for total entries
    # print('final', sum_sq, hist.Integral())
    return sqrt(sum_sq) / hist.Integral("width")
    # return sqrt(sum_sq) / sum_reduced


def calc_hist_mean_from_values(bin_areas, bin_centers):
    """Calculate mean of hist from value arrays.

    Must use np.X functions for e.g. sum(), square(), to ensure jax can differentiate it

    Parameters
    ----------
    bin_areas : TYPE
        Description
    bin_centers : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    return np.sum(bin_areas * bin_centers) / np.sum(bin_areas)


mean_differential = jit(grad(calc_hist_mean_from_values, argnums=0))


def calc_hist_mean_uncorrelated_error_from_values(bin_areas, bin_centers, bin_errors):
    """Calculate error on mean, assuming uncorrelated errors.

    Uses propagation of uncertainty bia partial differentials,
    calculated automatically using jax.

    Parameters
    ----------
    bin_areas : TYPE
        Description
    bin_centers : TYPE
        Description
    bin_errors : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    # differential wrt bin_areas
    diffs = mean_differential(bin_areas, bin_centers)
    err_sq = np.sum(np.square((diffs * bin_errors)))
    return np.sqrt(err_sq)


def calc_hist_mean_correlated_error_from_values(bin_areas, bin_widths, bin_centers, error_matrix):
    """Get error on mean, assuming covariance matrix error_matrix

    Note that we scale the error matrix by the bin width (since we also convert
    from bin height to bin area)

    Parameters
    ----------
    bin_areas : TYPE
        Description
    bin_widths : TYPE
        Description
    bin_centers : TYPE
        Description
    error_matrix : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    diffs = mean_differential(bin_areas, bin_centers)
    sum_sq = (diffs*bin_widths) @ error_matrix @ (diffs*bin_widths)
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
    mean = calc_hist_mean_from_values(areas, centers)
    err = calc_hist_mean_uncorrelated_error_from_values(areas, centers, errors)
    return float(mean), float(err)


def calc_hist_mean_and_correlated_error(hist, ematrix):
    areas, widths, centers, errors = hist_to_arrays(hist)
    mean = calc_hist_mean_from_values(areas, centers)
    err = calc_hist_mean_correlated_error_from_values(areas, widths, centers, ematrix)
    return float(mean), float(err)


def calc_hist_rms(hist):
    """Calculate std deviation for hist

    RMS^2 = (1/N) sum[ (bin_area - mean)^2 ]

    """
    areas, widths, centers, errors = hist_to_arrays(hist)
    return float(calc_hist_rms_from_values(areas, centers))
    # mean = calc_hist_mean(hist)
    # sum_sq = sum([pow((hist.GetBinWidth(ibin)*hist.GetBinContent(ibin)*hist.GetBinCenter(ibin) - mean), 2)
    #               for ibin in range(1, hist.GetNbinsX()+1)])
    # N = sum([hist.GetBinWidth(ibin)*hist.GetBinContent(ibin)
    #          for ibin in range(1, hist.GetNbinsX()+1)])
    # if N == 0:
    #     raise ValueError("N == 0")
    # return sqrt(sum_sq / N)


def calc_hist_rms_from_values(bin_areas, bin_centers):
    """Calculate RMS of hist from value arrays.

    Must use np.X functions for e.g. sum(), square(), to ensure jax can differentiate it

    Parameters
    ----------
    bin_areas : TYPE
        Description
    bin_centers : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    mean = calc_hist_mean_from_values(bin_areas, bin_centers)
    sum_sq = np.sum(np.square((bin_areas * bin_centers) - mean))
    return np.sqrt(sum_sq / np.sum(bin_areas))


rms_differential = jit(grad(calc_hist_rms_from_values, argnums=0))


def calc_hist_rms_uncorrelated_error_from_values(bin_areas, bin_centers, bin_errors):
    """Calculate error on RMS, assuming uncorrelated errors.

    Uses propagation of uncertainty bia partial differentials,
    calculated automatically using jax.

    Parameters
    ----------
    bin_areas : TYPE
        Description
    bin_centers : TYPE
        Description
    bin_errors : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    # differential wrt bin_areas
    diffs = rms_differential(bin_areas, bin_centers)
    err_sq = np.sum(np.square((diffs * bin_errors)))
    return np.sqrt(err_sq)


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
    rms = calc_hist_rms_from_values(areas, centers)
    err = calc_hist_rms_uncorrelated_error_from_values(areas, centers, errors)
    return float(rms), float(err)


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
                    input_tfile = cu.TFileCacher(root_filename)
                    pt_bins = qgc.PT_UNFOLD_DICT['signal_zpj_gen'] if 'ZPlusJets' in this_region['name'] else qgc.PT_UNFOLD_DICT['signal_gen']
                    unfolding_dict = unpack_slim_unfolding_root_file(input_tfile, this_region['name'], angle.var, pt_bins)

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
                        # mc_gen_hist_bin = hbc.get_pt_bin_normed_div_bin_width('hist_truth', ibin)
                        mc_gen_hist_bin = unfolding_dict['truth_hists'][ibin]
                        mc_gen_hist_bin_mean, mc_gen_hist_bin_mean_err = calc_hist_mean_and_uncorrelated_error(mc_gen_hist_bin)
                        mc_gen_hist_bin_rms, mc_gen_hist_bin_rms_err = calc_hist_rms_and_uncorrelated_error(mc_gen_hist_bin)

                        # Handle alt MC hist -> metrics
                        # alt_mc_gen_hist_bin = hbc.get_pt_bin_normed_div_bin_width('alt_hist_truth', ibin)
                        alt_mc_gen_hist_bin = unfolding_dict['alt_truth_hists'][ibin]
                        alt_mc_gen_hist_bin_mean, alt_mc_gen_hist_bin_mean_err = calc_hist_mean_and_uncorrelated_error(alt_mc_gen_hist_bin)
                        alt_mc_gen_hist_bin_rms, alt_mc_gen_hist_bin_rms_err = calc_hist_rms_and_uncorrelated_error(alt_mc_gen_hist_bin)

                        # Handle unfolded data hist -> metrics
                        # unfolded_hist_bin_total_errors = hbc.get_pt_bin_normed_div_bin_width('unfolded', ibin)
                        unfolded_hist_bin_total_errors = unfolding_dict['unfolding_total_err_hists'][ibin]
                        # unfolded_hist_bin_total_errors_mean, unfolded_hist_bin_total_errors_mean_err = calc_hist_mean_and_uncorrelated_error(unfolded_hist_bin_total_errors)
                        unfolded_hist_bin_total_errors_rms, unfolded_hist_bin_total_errors_rms_err = calc_hist_rms_and_uncorrelated_error(unfolded_hist_bin_total_errors)
                        unfolded_total_ematrix, _ = cu.th2_to_ndarray(unfolding_dict['unfolding_total_ematrices'][ibin])
                        unfolded_hist_bin_total_errors_mean, unfolded_hist_bin_total_errors_mean_err = calc_hist_mean_and_correlated_error(unfolded_hist_bin_total_errors, unfolded_total_ematrix)

                        print("uncorrelated value: ", calc_hist_mean_and_uncorrelated_error(unfolded_hist_bin_total_errors))
                        print("correlated value: ", calc_hist_mean_and_correlated_error(unfolded_hist_bin_total_errors, unfolded_total_ematrix))

                        # print("my mean error:", calc_hist_mean_error(mc_gen_hist_bin, this_cov_matrix, mc_gen_hist_bin_unnorm.Integral()))
                        # print("ROOT get mean error:", mc_gen_hist_bin.GetMeanError())

                        # print("ROOT mean:", mc_gen_hist_bin.GetMean())
                        # print("my mean:", calc_hist_mean(mc_gen_hist_bin))

                        # areas, centers, errors = hist_to_arrays(mc_gen_hist_bin)
                        # print("my mean:", calc_hist_mean_from_values(areas, centers))
                        # print("my mean error:", calc_hist_mean_uncorrelated_error(mc_gen_hist_bin))
                        # print("my mean error jax:", calc_hist_mean_uncorrelated_error_from_values(areas, centers, errors))

                        # print("ROOT RMS:", mc_gen_hist_bin.GetRMS())
                        # print("my RMS:", calc_hist_rms(mc_gen_hist_bin))
                        # print("my RMS from value:", calc_hist_rms_from_values(areas, centers))
                        # print("my RMS error jax:", calc_hist_rms_uncorrelated_error_from_values(areas, centers, errors))

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

    # --------------------------------------------------------------------------
    # Do all the plotting
    # --------------------------------------------------------------------------
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
        # plotter.plot_dijet_means_vs_pt_all()
        # plotter.plot_dijet_rms_vs_pt_all()

        pt_bins = qgc.PT_UNFOLD_DICT['signal_gen']
        low_pt_bin = np.where(pt_bins == low_pt)[0][0]
        low_pt_str = "[%g, %g] GeV" % (low_pt, pt_bins[low_pt_bin+1])

        high_pt_bin = np.where(pt_bins == high_pt)[0][0]
        high_pt_str = "[%g, %g] GeV" % (high_pt, pt_bins[high_pt_bin+1])

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

                # FIXME need groomed ak4 instead
                # ('jet_algo=="ak8puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_central" & angle=="%s_charged"' % (low_pt_bin, angle.var),
                    # "#splitline{{{jet_str}, {pt_str}}}{{Charged-only}}".format(jet_str=ak8_str, pt_str=low_pt_str)),
                    # "#splitline{{#splitline{{{jet_str}}}{{{pt_str}}}}}{{Charged-only}}".format(jet_str=ak8_str, pt_str=low_pt_str)),

                ('jet_algo=="ak4puppi" & pt_bin==%d & isgroomed  & region=="Dijet_central_groomed" & angle=="%s"' % (low_pt_bin, angle.var),
                    "#splitline{{{jet_str}, {pt_str}}}{{Groomed}}".format(jet_str=ak4_str, pt_str=low_pt_str)),
            ]
            selections.append({'label': this_angle_str, 'selections': this_selection})

        plotter.plot_mean_rms_bins_summary(
            selections=selections,
            output_file=os.path.join(args.outputDir, "dijet_mean_rms_summary.pdf")
        )

    if has_zpj:
        # plotter.plot_zpj_means_vs_pt_all()
        # plotter.plot_zpj_rms_vs_pt_all()

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

                # FIXME need groomed ak4 instead
                # ('jet_algo=="ak8puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s_charged"' % (low_pt_bin, angle.var),
                    # "#splitline{{{jet_str}, {pt_str}}}{{Charged-only}}".format(jet_str=ak8_str, pt_str=low_pt_str)),
                    # "#splitline{{#splitline{{{jet_str}}}{{{pt_str}}}}}{{Charged-only}}".format(jet_str=ak8_str, pt_str=low_pt_str)),

                ('jet_algo=="ak4puppi" & pt_bin==%d & isgroomed  & region=="ZPlusJets_groomed" & angle=="%s"' % (low_pt_bin, angle.var),
                    "#splitline{{{jet_str}, {pt_str}}}{{Groomed}}".format(jet_str=ak4_str, pt_str=low_pt_str)),
            ]
            selections.append({'label': this_angle_str, 'selections': this_selection})

        plotter.plot_mean_rms_bins_summary(
            selections=selections,
            output_file=os.path.join(args.outputDir, "zpj_mean_rms_summary.pdf")
        )


    # if has_dijet and has_zpj:
    #     plotter.plot_dijet_zpj_means_vs_pt_all()
        # plotter.plot_dijet_zpj_rms_vs_pt_all()

