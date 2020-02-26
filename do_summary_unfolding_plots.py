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
from itertools import product
from array import array
from math import sqrt

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from my_unfolder import MyUnfolder, unfolder_from_tdir, HistBinChopper, unpack_unfolding_root_file
from my_unfolder_plotter import MyUnfolderPlotter
from unfolding_config import get_dijet_config, get_zpj_config


ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()


class SummaryPlotter(object):
    """Do lots of summary plots"""

    def __init__(self, jet_algos, regions, angles, pt_bins_dijet, pt_bins_zpj, df, output_dir, has_data):
        self.jet_algos = jet_algos
        self.regions = regions
        self.angles = angles
        self.pt_bins_dijet = pt_bins_dijet
        self.pt_bins_zpj = pt_bins_zpj
        self.df = df
        self.output_fmt = 'pdf'
        self.output_dir = output_dir
        self.has_data = has_data
        self.mc_label = 'MG+Pythia8'
        self.alt_mc_label = 'Herwig++'

    @staticmethod
    def data_to_hist(data, data_err, bins):
        if len(data) != len(bins)-1:
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
        for jet_algo, angle, groomed in product(self.jet_algos[:1], self.angles[:], [False, True][:]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_means_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_zpj_means_vs_pt_all' % self.output_dir)

    def plot_dijet_means_vs_pt_all(self):
        """Plot mean vs pt for dijet (cen+fwd) on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_zpj_means_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos[:1], self.angles[:], [False, True][:]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_means_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_means_vs_pt_all' % self.output_dir, do_zpj=False)

    def plot_zpj_means_vs_pt_all(self):
        """Plot mean vs pt for zpj on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_zpj_means_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos[:1], self.angles[:], [False, True][:]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_means_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_zpj_means_vs_pt_all' % self.output_dir, do_dijet=False)

    def plot_dijet_zpj_means_vs_pt_one_angle_one_jet(self, angle, jet_algo, do_groomed, output_dir, do_zpj=True, do_dijet=True):
        """Do plot of mean lambda vs pt, for the dijet cen+fwd and zpj regions, for a given angle/jet algo/grooming"""
        df = self.df
        mask = ((df['angle'] == angle.var) & (df['jet_algo'] == jet_algo['name']) & (df['isgroomed'] == do_groomed))

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

        # unify this?
        dijet_cen_col = ROOT.kBlack
        dijet_fwd_col = ROOT.kRed
        zpj_col = ROOT.kBlue
        m_size = 1
        entries = []
        if do_dijet:
            entries.extend([
                Contribution(dijet_central_hist, label='Dijet (central)',
                             line_color=dijet_cen_col, line_width=2,
                             marker_color=dijet_cen_col, marker_style=cu.Marker.get('circle', True), marker_size=m_size),
                Contribution(dijet_forward_hist, label='Dijet (forward)',
                             line_color=dijet_fwd_col, line_width=2,
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', True), marker_size=m_size),
            ])
        if do_zpj:
            entries.extend([
                Contribution(zpj_hist, label='Z+jet',
                             line_color=zpj_col, line_width=2,
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', True), marker_size=m_size),
            ])
        if do_dijet:
            entries.extend([
                # TODO: make truth plotting optional, also plot alternate generators
                Contribution(dijet_central_hist_truth, label='Dijet (central) [%s]' % (self.mc_label),
                             line_color=dijet_cen_col, line_width=2, line_style=2,
                             marker_color=dijet_cen_col, marker_style=cu.Marker.get('circle', False), marker_size=0,
                             subplot=dijet_central_hist),
                Contribution(dijet_forward_hist_truth, label='Dijet (forward) [%s]' % (self.mc_label),
                             line_color=dijet_fwd_col, line_width=2, line_style=2,
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', False), marker_size=0,
                             subplot=dijet_forward_hist),
            ])
        if do_zpj:
            entries.extend([
                Contribution(zpj_hist_truth, label='Z+jet [%s]' % (self.mc_label),
                             line_color=zpj_col, line_width=2, line_style=2,
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', False), marker_size=0,
                             subplot=zpj_hist),
            ])
        if do_dijet:
            entries.extend([
                Contribution(dijet_central_hist_alt_truth, label='Dijet (central) [%s]' % (self.alt_mc_label),
                             line_color=dijet_cen_col, line_width=2, line_style=3,
                             marker_color=dijet_cen_col, marker_style=cu.Marker.get('circle', False), marker_size=0,
                             subplot=dijet_central_hist),
                Contribution(dijet_forward_hist_alt_truth, label='Dijet (forward) [%s]' % (self.alt_mc_label),
                             line_color=dijet_fwd_col, line_width=2, line_style=3,
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', False), marker_size=0,
                             subplot=dijet_forward_hist),
            ])
        if do_zpj:
            entries.extend([
                Contribution(zpj_hist_alt_truth, label='Z+jet [%s]' % (self.alt_mc_label),
                             line_color=zpj_col, line_width=2, line_style=3,
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', False), marker_size=0,
                             subplot=zpj_hist),
            ])

        angle_prepend = "Groomed " if do_groomed else ""
        this_angle_name = angle.name
        if (angle_prepend != ""
            and 'LHA' not in this_angle_name
            and "_{T}" not in this_angle_name
            and "PUPPI" not in this_angle_name):
            # lower case if Groomed..., but be careful of e.g. pTD, LHA
            this_angle_name = this_angle_name[0].lower() + this_angle_name[1:]
        # for plot axis titles
        angle_str = "#LT {prepend}{name}  {lambda_str} #GT".format(prepend=angle_prepend,
                                                             name=this_angle_name,
                                                             lambda_str=angle.lambda_str)
        h_max = max([c.obj.GetMaximum() for c in entries])
        plot = Plot(entries,
                    what='hist',
                    xtitle='Jet p_{T} [GeV]',
                    ytitle=angle_str,
                    title="%s jets" % (jet_algo['label']),
                    ylim=(0, h_max*1.75),
                    has_data=self.has_data,
                    subplot_type='ratio',
                    subplot_title='MC / Data',
                    subplot_limits=(0.5, 1.5) if self.has_data else (0.9, 1.1)
                    )
        plot.title_start_y = 0.85
        plot.title_left_offset = 0.05
        plot.title_font_size = 0.035
        plot.legend.SetX1(0.55)
        plot.legend.SetX2(0.88)
        plot.legend.SetY1(0.7)
        plot.legend.SetY2(0.85+0.02)
        plot.left_margin = 0.16
        # plot.default_canvas_size = (600, 800)
        plot.plot("NOSTACK E1 X0", "NOSTACK E1")
        # plot.get_modifier().GetYaxis().SetTitleOffset(plot.get_modifier().GetYaxis().GetTitleOffset()*1.1)
        plot.set_logx(do_more_labels=False)
        groomed_str = '_groomed' if do_groomed else ''
        plot.save("%s/dijet_zpj_means_vs_pt_%s%s_%s.%s" % (output_dir, angle.var, groomed_str, jet_algo['name'], self.output_fmt))

    def plot_dijet_zpj_rms_vs_pt_all(self):
        """Plot mean vs pt for dijet (cen+fwd) & Z+jet on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_zpj_rms_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos[:1], self.angles[:], [False, True][:]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_rms_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_zpj_rms_vs_pt_all' % self.output_dir)

    def plot_dijet_rms_vs_pt_all(self):
        """Plot mean vs pt for dijet (cen+fwd) on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_zpj_rms_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos[:1], self.angles[:], [False, True][:]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_rms_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_rms_vs_pt_all' % self.output_dir, do_zpj=False)

    def plot_zpj_rms_vs_pt_all(self):
        """Plot mean vs pt for zpj on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_zpj_rms_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos[:1], self.angles[:], [False, True][:]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_rms_vs_pt_one_angle_one_jet(angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_zpj_rms_vs_pt_all' % self.output_dir, do_dijet=False)

    def plot_dijet_zpj_rms_vs_pt_one_angle_one_jet(self, angle, jet_algo, do_groomed, output_dir, do_zpj=True, do_dijet=True):
        """Do plot of RMS lambda vs pt, for the dijet cen+fwd and zpj regions, for a given angle/jet algo/grooming"""
        df = self.df
        mask = ((df['angle'] == angle.var) & (df['jet_algo'] == jet_algo['name']) & (df['isgroomed'] == do_groomed))

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
        dijet_cen_col = ROOT.kBlack
        dijet_fwd_col = ROOT.kRed
        zpj_col = ROOT.kBlue
        m_size = 1
        mc_label='MG+Pythia8'
        entries = []
        if do_dijet:
            entries.extend([
                Contribution(dijet_central_hist, label='Dijet (central)',
                             line_color=dijet_cen_col, line_width=2,
                             marker_color=dijet_cen_col, marker_style=cu.Marker.get('circle', True), marker_size=m_size),
                Contribution(dijet_forward_hist, label='Dijet (forward)',
                             line_color=dijet_fwd_col, line_width=2,
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', True), marker_size=m_size),
            ])
        if do_zpj:
            entries.extend([

                Contribution(zpj_hist, label='Z+jet',
                             line_color=zpj_col, line_width=2,
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', True), marker_size=m_size),
            ])
        if do_dijet:
            entries.extend([
                # TODO: make truth plotting optional, also plot alternate generators
                Contribution(dijet_central_hist_truth, label='Dijet (central) [%s]' % (self.mc_label),
                             line_color=dijet_cen_col, line_width=2, line_style=2,
                             marker_color=dijet_cen_col, marker_style=cu.Marker.get('circle', False), marker_size=0,
                             subplot=dijet_central_hist),
                Contribution(dijet_forward_hist_truth, label='Dijet (forward) [%s]' % (self.mc_label),
                             line_color=dijet_fwd_col, line_width=2, line_style=2,
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', False), marker_size=0,
                             subplot=dijet_forward_hist),
            ])

        if do_zpj:
            entries.extend([
                Contribution(zpj_hist_truth, label='Z+jet [%s]' % (self.mc_label),
                             line_color=zpj_col, line_width=2, line_style=2,
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', False), marker_size=0,
                             subplot=zpj_hist),
                ])
        if do_dijet:
            entries.extend([
                Contribution(dijet_central_hist_alt_truth, label='Dijet (central) [%s]' % (self.alt_mc_label),
                             line_color=dijet_cen_col, line_width=2, line_style=3,
                             marker_color=dijet_cen_col, marker_style=cu.Marker.get('circle', False), marker_size=0,
                             subplot=dijet_central_hist),
                Contribution(dijet_forward_hist_alt_truth, label='Dijet (forward) [%s]' % (self.alt_mc_label),
                             line_color=dijet_fwd_col, line_width=2, line_style=3,
                             marker_color=dijet_fwd_col, marker_style=cu.Marker.get('square', False), marker_size=0,
                             subplot=dijet_forward_hist),
            ])
        if do_zpj:
            entries.extend([
                Contribution(zpj_hist_alt_truth, label='Z+jet [%s]' % (self.alt_mc_label),
                             line_color=zpj_col, line_width=2, line_style=3,
                             marker_color=zpj_col, marker_style=cu.Marker.get('triangleUp', False), marker_size=0,
                             subplot=zpj_hist),
            ])

        angle_prepend = "Groomed " if do_groomed else ""
        this_angle_name = angle.name
        if (angle_prepend != ""
            and 'LHA' not in this_angle_name
            and "_{T}" not in this_angle_name
            and "PUPPI" not in this_angle_name):
            # lower case if Groomed..., but be careful of e.g. pTD, LHA
            this_angle_name = this_angle_name[0].lower() + this_angle_name[1:]
        # for plot axis titles
        angle_str = "RMS {prepend}{name}  {lambda_str}".format(prepend=angle_prepend,
                                                             name=this_angle_name,
                                                             lambda_str=angle.lambda_str)
        h_max = max([c.obj.GetMaximum() for c in entries])
        plot = Plot(entries,
                    what='hist',
                    xtitle='Jet p_{T} [GeV]',
                    ytitle=angle_str,
                    title="%s jets" % (jet_algo['label']),
                    # ylim=(0, h_max*1.75),
                    has_data=self.has_data,
                    subplot_type='ratio',
                    subplot_title='MC / Data',
                    subplot_limits=(0.5, 1.5) if self.has_data else (0.9, 1.1))
        plot.title_start_y = 0.85
        plot.title_left_offset = 0.05
        plot.title_font_size = 0.035
        plot.legend.SetX1(0.55)
        plot.legend.SetY1(0.7)
        plot.legend.SetX2(0.88)
        plot.legend.SetY2(0.85+0.02)
        plot.left_margin = 0.16
        plot.y_padding_min_linear = 0.5
        # plot.default_canvas_size = (600, 800)
        plot.plot("NOSTACK E1 X0", "NOSTACK E1")
        # plot.get_modifier().GetYaxis().SetTitleOffset(plot.get_modifier().GetYaxis().GetTitleOffset()*1.1)
        plot.set_logx(do_more_labels=False)
        groomed_str = '_groomed' if do_groomed else ''
        plot.save("%s/dijet_zpj_rms_vs_pt_%s%s_%s.%s" % (output_dir, angle.var, groomed_str, jet_algo['name'], self.output_fmt))


def calc_hist_mean(hist):
    return sum([hist.GetBinContent(ibin)*hist.GetBinWidth(ibin)*hist.GetBinCenter(ibin) for ibin in range(1, hist.GetNbinsX()+1)]) / sum([hist.GetBinContent(ibin)*hist.GetBinWidth(ibin) for ibin in range(1, hist.GetNbinsX()+1)])


def calc_hist_mean_error(hist, covariance_matrix, scale_factor):
    """Get error on the mean using covariance matrix

    covariance_matrix shoudl be a 2D numpy array for this bin

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
    return sqrt(sum_sq) / hist.Integral()
    # return sqrt(sum_sq) / sum_reduced


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ak4source",
                        help="Source directory for AK4 jets (should be the one made by unfolding.py")
    parser.add_argument("--ak8source",
                        help="Source directory for AK8 jets (should be the one made by unfolding.py")
    parser.add_argument("--outputDir",
                        default=None,
                        help='Output directory (default is the source dir')
    args = parser.parse_args()
    if not args.outputDir:
        args.outputDir = os.path.join(args.ak4source, 'SummaryPlots')

    results_dicts = []

    jet_algos = [
        {'src': args.ak4source, 'label': 'AK4 PUPPI', 'name': 'ak4puppi'},
    ]
    if args.ak8source:
        jet_algos.append({'src': args.ak8source, 'label': 'AK8 PUPPI', 'name': 'ak8puppi'})

    for jet_algo in jet_algos:
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

            angles = qgc.COMMON_VARS
            for angle in angles:
                angle_output_dir = "%s/%s" % (region_dir, angle.var)
                if not os.path.isdir(angle_output_dir):
                    print("! Warning ! cannot find angle dir", angle_output_dir, '- skipping angle', angle.var)
                    continue

                # TODO: put this in a method / class
                # Check if ROOT file exists
                root_filename = os.path.join(angle_output_dir, "unfolding_result.root")
                if not os.path.isfile(root_filename):
                    print("! Warning ! cannot fine unfolding ROOT file", root_filename, ' - skipping angle')

                print('Unpacking', region['name'], angle.name)
                append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name
                input_tfile = cu.TFileCacher(root_filename)  # keep this here otherwise crashes
                unpack_dict = unpack_unfolding_root_file(input_tfile, region, angle, do_alt_response=True, do_model_systs=False, do_pdf_systs=False)
                unfolder = unpack_dict['unfolder']
                unreg_unfolder = unpack_dict['unreg_unfolder']
                alt_unfolder = unpack_dict['alt_unfolder']
                alt_hist_truth = unpack_dict['alt_hist_truth']
                print(alt_hist_truth)
                hbc = HistBinChopper(unfolder)
                hbc.add_obj("unfolded", unfolder.unfolded)
                hbc.add_obj("unfolded_stat_err", unfolder.unfolded_stat_err)
                hbc.add_obj("hist_truth", unfolder.hist_truth)
                hbc.add_obj("alt_hist_truth", alt_hist_truth)

                # Iterate through pt bins, get lambda histogram for that bin,
                # derive metrics from it, save
                key = 'signal_gen'
                if 'ZPlusJets' in region['name']:
                    key = 'signal_zpj_gen'
                pt_bins = qgc.PT_UNFOLD_DICT[key]

                for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(pt_bins[:-1], pt_bins[1:])):
                    # mc_gen_hist_bin_unnorm_div_bin_width = hbc.get_pt_bin_div_bin_width('hist_truth', ibin)
                    # mc_gen_hist_bin = hbc.get_pt_bin_normed_div_bin_width('hist_truth', ibin)
                    mc_gen_hist_bin_unnorm = hbc.get_pt_bin('hist_truth', ibin)
                    # print('unnorm', mc_gen_hist_bin_unnorm.Integral())
                    # print('unnorm_div_bin_width', mc_gen_hist_bin_unnorm_div_bin_width.Integral())
                    # print('norm', mc_gen_hist_bin.Integral())
                    # unfolded_hist_bin_stat_errors = hbc.get_pt_bin_normed_div_bin_width('unfolded_stat_err', ibin)
                    # unfolded_hist_bin_total_errors = hbc.get_pt_bin_normed_div_bin_width('unfolded', ibin)
                    unfolded_hist_bin_unnorm = hbc.get_pt_bin('unfolded', ibin)

                    # alt_mc_gen_hist_bin = hbc.get_pt_bin_normed_div_bin_width('alt_hist_truth', ibin)
                    alt_mc_gen_hist_bin_unnorm = hbc.get_pt_bin('alt_hist_truth', ibin)

                    # GetMean() equivalent to
                    # sum([hist.GetBinContent(ibin)*hist.GetBinCenter(ibin) for ibin in range(1, hist.GetNbinsX()+1)]) /
                    #  sum([hist.GetBinContent(ibin) for ibin in range(1, hist.GetNbinsX()+1)])

                    # print("ROOT mean:", mc_gen_hist_bin.GetMean())
                    # print("my mean:", calc_hist_mean(mc_gen_hist_bin))
                    start_lambda = unfolder.variable_bin_edges_gen[0]*1.0000001
                    end_lambda = unfolder.variable_bin_edges_gen[-2]*1.0000001
                    this_pt = bin_edge_low*1.0000001
                    start_bin = unfolder.generator_distribution.GetGlobalBinNumber(start_lambda, this_pt)
                    end_bin = unfolder.generator_distribution.GetGlobalBinNumber(end_lambda, this_pt)
                    # -1 as ROOT bins start at 1, numpy starts at 0
                    this_cov_matrix = unfolder.ematrix_total_ndarray[start_bin-1:end_bin, start_bin-1:end_bin]
                    # print("my mean error:", calc_hist_mean_error(mc_gen_hist_bin, this_cov_matrix, mc_gen_hist_bin_unnorm.Integral()))
                    # print("ROOT get mean error:", mc_gen_hist_bin.GetMeanError())


                    result_dict = {
                        'jet_algo': jet_algo['name'],
                        'region': region['name'],
                        'isgroomed': 'groomed' in region['name'].lower(),
                        'pt_bin': ibin,
                        'angle': angle.var,
                        # 'mean': calc_hist_mean(unfolded_hist_bin_total_errors),
                        'mean': unfolded_hist_bin_unnorm.GetMean(),
                        # 'mean_err': unfolded_hist_bin_total_errors.GetMeanError(),
                        'mean_err': calc_hist_mean_error(unfolded_hist_bin_unnorm, this_cov_matrix, 1),

                        'mean_truth': mc_gen_hist_bin_unnorm.GetMean(),
                        # 'mean_err_truth': mc_gen_hist_bin.GetMeanError(),
                        'mean_err_truth': 0,

                        'mean_alt_truth': alt_mc_gen_hist_bin_unnorm.GetMean(),
                        'mean_err_alt_truth': 0,

                        'rms': unfolded_hist_bin_unnorm.GetRMS(),
                        'rms_err': unfolded_hist_bin_unnorm.GetRMSError(),  # FIXME

                        'rms_truth': mc_gen_hist_bin_unnorm.GetRMS(),
                        'rms_err_truth': 0,

                        'rms_alt_truth': alt_mc_gen_hist_bin_unnorm.GetRMS(),
                        'rms_err_alt_truth': 0,
                    }
                    results_dicts.append(result_dict)

    df = pd.DataFrame(results_dicts)
    df['jet_algo'] = df['jet_algo'].astype('category')
    df['region'] = df['region'].astype('category')
    df['angle'] = df['angle'].astype('category')
    print(df.head())
    print(len(df.index), 'entries in dataframe')
    print(df.dtypes)

    plotter = SummaryPlotter(jet_algos,
                             regions,
                             angles,
                             qgc.PT_UNFOLD_DICT['signal_gen'],
                             qgc.PT_UNFOLD_DICT['signal_zpj_gen'],
                             df,
                             args.outputDir,
                             has_data='_DATA_' in args.ak4source)
    plotter.plot_dijet_zpj_means_vs_pt_all()
    plotter.plot_dijet_means_vs_pt_all()
    plotter.plot_zpj_means_vs_pt_all()

    plotter.plot_dijet_zpj_rms_vs_pt_all()
    plotter.plot_dijet_rms_vs_pt_all()
    plotter.plot_zpj_rms_vs_pt_all()
