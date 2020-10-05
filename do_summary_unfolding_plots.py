#!/usr/bin/env python


"""
Do all the unfolding plots: per pT bin, per lambda bin, summary plot

Input is HDF5 file made by extract_unfolding_summary_stats.py
"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse
import pandas as pd
import numpy as np
from itertools import product
from array import array
from copy import copy

import yoda

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from unfolding_config import get_dijet_config, get_zpj_config
import rivet_naming as rn
import metric_calculators as metrics
from extract_rivet_summary_stats import get_dataframe_from_yoda_inputs, dataframe_yoda_key, convert_df_types


ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()

# import hunter
# hunter.trace(module='my_unfolder', action=hunter.CallPrinter)

# Define own linestyle for smaller plots
# linestyle 2 (dashed) has too big dashes
ROOT.gStyle.SetLineStyleString(22, "8 6")
ROOT.gStyle.SetLineStyleString(23, "5 10")

COMMON_STYLE_DICT = {
    "line_width": 2,
    "dijet_cen_color": ROOT.kBlack,
    "dijet_fwd_color": ROOT.kRed,
    "zpj_color": ROOT.kBlue,
    "jet_pt_units_str": "Jet p_{T} [GeV]",

    "data_line_style": 1,
    "data_color": ROOT.kBlack,
    "data_marker_style": cu.Marker.get("circle", filled=True),

    "marker_size": 1,

    # "mc_line_style": 22,
    "mc_line_style": 1,
    # format: 3ijk,
    # i=distance between lines,
    # j=angle between 0 and 90 degrees (5 = not drawn),
    # k=angle between 90 and 180 degrees (5 = not drawn)
    "mc_fill_style": 3445,
    "mc_color": ROOT.kBlue+1,
    "mc_marker_style": cu.Marker.get('square', filled=False),

    # "mc_alt_line_style": 22,
    "mc_alt_line_style": 1,
    "mc_alt_fill_style": 3454,
    "mc_alt_color": ROOT.kAzure+1,
    "mc_alt_marker_style": cu.Marker.get('triangleUp', filled=False),
}


# some pre-defined sample dicts for YODA files
SAMPLE_STYLE_DICTS = {
    "Sherpa": {
        "color": ROOT.kGreen+2,
        "label": "Sherpa",
        "marker_style": cu.Marker.get('triangleDown', filled=False),
        "marker_size": COMMON_STYLE_DICT['marker_size'] * 1,
    },

    "Herwig7 CH3": {
        "color": ROOT.kOrange-3,
        "label": "Herwig7 CH3",
        "marker_style": cu.Marker.get('star', filled=False),
        "marker_size": COMMON_STYLE_DICT['marker_size'] * 1.2, # these shapes always come out small
    },

    "Herwig7 CH3 alphaS=0.136": {
        "color": ROOT.kRed-2,
        "label": "Herwig7 CH3\n#alpha_{S} = 0.136",
        "marker_style": cu.Marker.get('diamond', filled=False),
        "marker_size": COMMON_STYLE_DICT['marker_size'] * 1.2, # these shapes always come out small
    },

    "Pythia8 CP5": {
        "color": ROOT.kMagenta-7,
        "label": "Pythia8 CP5",
        "marker_style": cu.Marker.get('crossX', filled=False),
        "marker_size": COMMON_STYLE_DICT['marker_size'] * 1.2, # these shapes always come out small
    }
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

    def __init__(self, jet_algos, regions, angles, pt_bins_dijet, pt_bins_zpj, df, output_dir, has_data, only_yoda_data):
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
        self.mc_label = 'MG5+Pythia8\nCUETP8M1'
        self.alt_mc_label = 'Herwig++'
        self.other_samples = []
        self.only_yoda_data = only_yoda_data
        self.filename_append = "_onlyYodaData" if self.only_yoda_data else ""

    def add_sample(self, key, label, style_dict):
        self.other_samples.append({
            "key": key,
            "label": label,
            "style_dict": style_dict
        })

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
    def _generate_filename_prefix(do_dijet_cen, do_dijet_fwd, do_zpj):
        this_str = ""
        if do_dijet_cen:
            this_str += "dijet_cen_"
        if do_dijet_fwd:
            this_str += "dijet_fwd_"
        if do_zpj:
            this_str += "zpj_"
        return this_str

    def plot_dijet_zpj_metric_vs_pt_one_angle_one_jet(self, metric, angle, jet_algo, do_groomed, output_dir, do_zpj=True, do_dijet_cen=True, do_dijet_fwd=True):
        """Do plot of lambda metric vs pt, for any combo of the dijet cen/fwd and zpj regions, for a given angle/jet algo/grooming"""
        if metric not in ["mean", "rms", "delta"]:
            raise ValueError("metric must be one of 'mean', 'rms', 'delta'")

        only_one_region = sum([do_dijet_cen, do_dijet_fwd, do_zpj]) == 1

        df = self.df
        mask = ((df['angle'] == angle.var) & (df['jet_algo'] == jet_algo['name']) & (df['isgroomed'] == do_groomed))
        if not mask.any():
            return

        # create hists for dijet central
        dijet_cen_col = COMMON_STYLE_DICT['dijet_cen_color']
        other_samples_dijet_central_hists = []
        dijet_central_hist_no_errors = None
        if do_dijet_cen:
            region_name = 'Dijet_central'
            if do_groomed:
                region_name += "_groomed"
            dijet_central_data = df[mask & (df['region'] == region_name)]

            if not self.only_yoda_data:
                dijet_central_hist_truth = self.data_to_hist(dijet_central_data['%s_truth' % metric], dijet_central_data['%s_err_truth' % metric], self.pt_bins_dijet)
                dijet_central_hist_alt_truth = self.data_to_hist(dijet_central_data['%s_alt_truth' % metric], dijet_central_data['%s_err_alt_truth' % metric], self.pt_bins_dijet)

            for sample in self.other_samples:
                hist = self.data_to_hist(dijet_central_data['%s_%s' % (metric, sample['key'])], dijet_central_data['%s_err_%s' % (metric, sample['key'])], self.pt_bins_dijet)
                other_samples_dijet_central_hists.append(hist)

            if metric != 'delta':
                dijet_central_hist = self.data_to_hist(dijet_central_data[metric], dijet_central_data['%s_err' % metric], self.pt_bins_dijet)
                # Create copy with 0 error bars, for the ratio subplot
                # The data will have its own error region
                dijet_central_hist_no_errors = dijet_central_hist.Clone()
                cu.remove_th1_errors(dijet_central_hist_no_errors)
                # Create hists for data with error reigon for ratio
                # Easiest way to get errors right is to do data (with 0 errors)
                # and divide by data (with errors), as if you had MC = data with 0 error
                dijet_central_hist_ratio_error = dijet_central_hist_no_errors.Clone()
                dijet_central_hist_ratio_error.Divide(dijet_central_hist)
                dijet_central_hist_ratio_error.SetFillStyle(3245)
                dijet_central_hist_ratio_error.SetFillColor(ROOT.kBlack if only_one_region else dijet_cen_col)
                dijet_central_hist_ratio_error.SetLineWidth(0)
                dijet_central_hist_ratio_error.SetMarkerSize(0)

        # create hists for dijet forward
        dijet_fwd_col = COMMON_STYLE_DICT['dijet_fwd_color']
        other_samples_dijet_forward_hists = []
        dijet_forward_hist_no_errors = None
        if do_dijet_fwd:
            region_name = 'Dijet_forward'
            if do_groomed:
                region_name += "_groomed"
            dijet_forward_data = df[mask & (df['region'] == region_name)]

            if not self.only_yoda_data:
                dijet_forward_hist_truth = self.data_to_hist(dijet_forward_data['%s_truth' % metric], dijet_forward_data['%s_err_truth' % metric], self.pt_bins_dijet)
                dijet_forward_hist_alt_truth = self.data_to_hist(dijet_forward_data['%s_alt_truth' % metric], dijet_forward_data['%s_err_alt_truth' % metric], self.pt_bins_dijet)

            for sample in self.other_samples:
                hist = self.data_to_hist(dijet_forward_data['%s_%s' % (metric, sample['key'])], dijet_forward_data['%s_err_%s' % (metric, sample['key'])], self.pt_bins_dijet)
                other_samples_dijet_forward_hists.append(hist)

            if metric != 'delta':
                dijet_forward_hist = self.data_to_hist(dijet_forward_data[metric], dijet_forward_data['%s_err' % metric], self.pt_bins_dijet)
                dijet_forward_hist_no_errors = dijet_forward_hist.Clone()
                cu.remove_th1_errors(dijet_forward_hist.Clone())

                dijet_forward_hist_ratio_error = dijet_forward_hist_no_errors.Clone()
                dijet_forward_hist_ratio_error.Divide(dijet_forward_hist)
                dijet_forward_hist_ratio_error.SetFillStyle(3254)
                dijet_forward_hist_ratio_error.SetFillColor(ROOT.kBlack if only_one_region else dijet_fwd_col)
                dijet_forward_hist_ratio_error.SetLineWidth(0)
                dijet_forward_hist_ratio_error.SetMarkerSize(0)

        # create hists for Z+jets
        zpj_col = COMMON_STYLE_DICT['zpj_color']
        other_samples_zpj_hists = []
        zpj_hist_no_errors = None
        if do_zpj:
            region_name = 'ZPlusJets'
            if do_groomed:
                region_name += "_groomed"
            # drop last pt bin as massive error
            pt_bins = self.pt_bins_zpj[:-2]
            zpj_data = df[mask & (df['region'] == region_name) & (df['pt_bin'] < len(pt_bins)-1)]
            # the -1 is because the last entry in pt_bins is the upper edge of the last bin
            if not self.only_yoda_data:
                zpj_hist_truth = self.data_to_hist(zpj_data['%s_truth' % metric], zpj_data['%s_err_truth' % metric], pt_bins)
                zpj_hist_alt_truth = self.data_to_hist(zpj_data['%s_alt_truth' % metric], zpj_data['%s_err_alt_truth' % metric], pt_bins)

            for sample in self.other_samples:
                key = sample['key']
                hist = self.data_to_hist(zpj_data['%s_%s' % (metric, key)], zpj_data['%s_err_%s' % (metric, key)], pt_bins)
                other_samples_zpj_hists.append(hist)

            if metric != 'delta':
                zpj_hist = self.data_to_hist(zpj_data[metric], zpj_data['%s_err' % metric], pt_bins)
                # 0 error bar hists & ratio = 1 hists for subplot
                zpj_hist_no_errors = zpj_hist.Clone()
                cu.remove_th1_errors(zpj_hist_no_errors)
                zpj_hist_ratio_error = zpj_hist_no_errors.Clone()
                zpj_hist_ratio_error.Divide(zpj_hist)
                zpj_hist_ratio_error.SetFillStyle(3003)
                zpj_hist_ratio_error.SetFillStyle(3254)
                zpj_hist_ratio_error.SetFillColor(ROOT.kBlack if only_one_region else zpj_col)
                zpj_hist_ratio_error.SetLineWidth(0)
                zpj_hist_ratio_error.SetMarkerSize(0)

        m_size = COMMON_STYLE_DICT['marker_size']
        lw = COMMON_STYLE_DICT['line_width']
        entries = []

        # Create dummy graphs with the same styling to put into the legend,
        # since graphs get error bar ends, but hists don't
        dummy_gr = ROOT.TGraphErrors(1, array('d', [1]), array('d', [1]), array('d', [1]), array('d', [1]))
        dummy_entries = []
        # NB Spaces in legend labels are important for padding

        # Add data
        if metric != 'delta':
            if do_dijet_cen:
                cont_args = dict(label='Data' if only_one_region else ' Dijet (central)',
                                 leg_draw_opt="LEP",
                                 line_color=ROOT.kBlack if only_one_region else dijet_cen_col,
                                 line_width=lw,
                                 marker_color=ROOT.kBlack if only_one_region else dijet_cen_col,
                                 marker_style=COMMON_STYLE_DICT['data_marker_style'],
                                 marker_size=m_size)
                entries.append(Contribution(dijet_central_hist, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

            if do_dijet_fwd:
                cont_args = dict(label='Data' if only_one_region else ' Dijet (forward)',
                                 leg_draw_opt="LEP",
                                 line_color=ROOT.kBlack if only_one_region else dijet_fwd_col,
                                 line_width=lw,
                                 marker_color=ROOT.kBlack if only_one_region else dijet_fwd_col,
                                 marker_style=COMMON_STYLE_DICT['data_marker_style'],
                                 marker_size=m_size)
                entries.append(Contribution(dijet_forward_hist, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

            if do_zpj:
                cont_args = dict(label='Data' if only_one_region else ' Z+jets',
                                 leg_draw_opt="LEP",
                                 line_color=ROOT.kBlack if only_one_region else zpj_col,
                                 line_width=lw,
                                 marker_color=ROOT.kBlack if only_one_region else zpj_col,
                                 marker_style=COMMON_STYLE_DICT['data_marker_style'],
                                 marker_size=m_size)
                entries.append(Contribution(zpj_hist, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

        if not self.only_yoda_data:
            # Add nominal MC
            if do_dijet_cen:
                cont_args = dict(label=self.mc_label if only_one_region else '#splitline{ Dijet (central)  }{ [%s]}' % (self.mc_label),
                                 line_color=COMMON_STYLE_DICT['mc_color'],
                                 line_width=lw,
                                 line_style=COMMON_STYLE_DICT['mc_line_style'],
                                 marker_color=COMMON_STYLE_DICT['mc_color'],
                                 marker_style=COMMON_STYLE_DICT['mc_marker_style'],
                                 marker_size=m_size,
                                 leg_draw_opt="LE" if m_size == 0 else "EP",
                                 subplot=dijet_central_hist_no_errors)
                entries.append(Contribution(dijet_central_hist_truth, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

            if do_dijet_fwd:
                cont_args = dict(label=self.mc_label if only_one_region else '#splitline{ Dijet (forward)  }{ [%s]}' % (self.mc_label),
                                 line_color=COMMON_STYLE_DICT['mc_color'],
                                 line_width=lw,
                                 line_style=COMMON_STYLE_DICT['mc_line_style'],
                                 marker_color=COMMON_STYLE_DICT['mc_color'],
                                 marker_style=COMMON_STYLE_DICT['mc_marker_style'],
                                 marker_size=m_size,
                                 leg_draw_opt="LE" if m_size == 0 else "EP",
                                 subplot=dijet_forward_hist_no_errors)
                entries.append(Contribution(dijet_forward_hist_truth, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

            if do_zpj:
                cont_args = dict(label=self.mc_label if only_one_region else '#splitline{ Z+jets  }{ [%s]}' % (self.mc_label),
                                 line_color=COMMON_STYLE_DICT['mc_color'],
                                 line_width=lw,
                                 line_style=COMMON_STYLE_DICT['mc_line_style'],
                                 marker_color=COMMON_STYLE_DICT['mc_color'],
                                 marker_style=COMMON_STYLE_DICT['mc_marker_style'],
                                 marker_size=m_size,
                                 leg_draw_opt="LE" if m_size == 0 else "EP",
                                 subplot=zpj_hist_no_errors)
                entries.append(Contribution(zpj_hist_truth, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

            # add alt MC
            if do_dijet_cen:
                cont_args = dict(label=self.alt_mc_label if only_one_region else '#splitline{ Dijet (central)  }{ [%s]}' % (self.alt_mc_label),
                                 line_color=COMMON_STYLE_DICT['mc_alt_color'],
                                 line_width=lw,
                                 line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                                 marker_color=COMMON_STYLE_DICT['mc_alt_color'],
                                 marker_style=COMMON_STYLE_DICT['mc_alt_marker_style'],
                                 marker_size=m_size,
                                 leg_draw_opt="LE" if m_size == 0 else "EP",
                                 subplot=dijet_central_hist_no_errors)
                entries.append(Contribution(dijet_central_hist_alt_truth, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

            if do_dijet_fwd:
                cont_args = dict(label=self.alt_mc_label if only_one_region else '#splitline{ Dijet (forward)  }{ [%s]}' % (self.alt_mc_label),
                                 line_color=COMMON_STYLE_DICT['mc_alt_color'],
                                 line_width=lw,
                                 line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                                 marker_color=COMMON_STYLE_DICT['mc_alt_color'],
                                 marker_style=COMMON_STYLE_DICT['mc_alt_marker_style'],
                                 marker_size=m_size,
                                 leg_draw_opt="LE" if m_size == 0 else "EP",
                                 subplot=dijet_forward_hist_no_errors)
                entries.append(Contribution(dijet_forward_hist_alt_truth, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

            if do_zpj:
                cont_args = dict(label=self.alt_mc_label if only_one_region else '#splitline{ Z+jets}{ [%s]}' % (self.alt_mc_label),
                                 line_color=COMMON_STYLE_DICT['mc_alt_color'],
                                 line_width=lw,
                                 line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                                 marker_color=COMMON_STYLE_DICT['mc_alt_color'],
                                 marker_style=COMMON_STYLE_DICT['mc_alt_marker_style'],
                                 marker_size=m_size,
                                 leg_draw_opt="LE" if m_size == 0 else "EP",
                                 subplot=zpj_hist_no_errors)
                entries.append(Contribution(zpj_hist_alt_truth, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

        # Add other samples: dijet cen
        if do_dijet_cen:
            marker = cu.Marker(shape='triangleDown')
            for sample, hist, mark in zip(self.other_samples, other_samples_dijet_central_hists, marker.cycle()):
                style_dict = sample['style_dict']
                color = style_dict.get('color', COMMON_STYLE_DICT['mc_alt_color'])
                fancy_label = style_dict.get("label", sample['label'])
                marker_size = style_dict.get('marker_size', m_size)
                cont_args = dict(label=fancy_label if only_one_region else '#splitline{ Dijet (central) }{ [%s]}' % (fancy_label),
                                 line_color=color,
                                 line_width=lw,
                                 line_style=style_dict.get('line_style', COMMON_STYLE_DICT['mc_alt_line_style']),
                                 marker_color=color,
                                 marker_style=style_dict.get('marker_style', mark),
                                 marker_size=marker_size,
                                 leg_draw_opt="LE" if marker_size == 0 else "EP",
                                 subplot=dijet_central_hist_no_errors)
                entries.append(Contribution(hist, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

        # Add other samples: dijet fwd
        if do_dijet_fwd:
            marker = cu.Marker(shape='triangleDown')
            for sample, hist, mark in zip(self.other_samples, other_samples_dijet_forward_hists, marker.cycle()):
                style_dict = sample['style_dict']
                color = style_dict.get('color', COMMON_STYLE_DICT['mc_alt_color'])
                fancy_label = style_dict.get("label", sample['label'])
                marker_size = style_dict.get('marker_size', m_size)
                cont_args = dict(label=fancy_label if only_one_region else '#splitline{ Dijet (forward) }{ [%s]}' % (fancy_label),
                                 line_color=color,
                                 line_width=lw,
                                 line_style=style_dict.get('line_style', COMMON_STYLE_DICT['mc_alt_line_style']),
                                 marker_color=color,
                                 marker_style=style_dict.get('marker_style', mark),
                                 marker_size=marker_size,
                                 leg_draw_opt="LE" if marker_size == 0 else "EP",
                                 subplot=dijet_forward_hist_no_errors)
                entries.append(Contribution(hist, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

        # Add other samples: Z+jet
        if do_zpj:
            marker = cu.Marker(shape='triangleDown')
            for sample, hist, mark in zip(self.other_samples, other_samples_zpj_hists, marker.cycle()):
                style_dict = sample['style_dict']
                color = style_dict.get('color', COMMON_STYLE_DICT['mc_alt_color'])
                fancy_label = style_dict.get("label", sample['label'])
                marker_size = style_dict.get('marker_size', m_size)
                cont_args = dict(label=fancy_label if only_one_region else '#splitline{ Z+jets}{ [%s]}' % (fancy_label),
                                 line_color=color,
                                 line_width=lw,
                                 line_style=style_dict.get('line_style', COMMON_STYLE_DICT['mc_alt_line_style']),
                                 marker_color=color,
                                 marker_style=style_dict.get('marker_style', mark),
                                 marker_size=marker_size,
                                 leg_draw_opt="LE" if marker_size == 0 else "EP",
                                 subplot=zpj_hist_no_errors)
                entries.append(Contribution(hist, **cont_args))
                dummy_entries.append(Contribution(dummy_gr.Clone(), **cont_args))

        # for plot axis titles
        if metric == "mean":
            angle_str = "#LT %s #GT" % create_angle_label(angle, do_groomed)
        elif metric == "rms":
            angle_str = "RMS %s" % create_angle_label(angle, do_groomed)
        elif metric == "delta":
            angle_str = "#Delta, %s" % create_angle_label(angle, do_groomed)

        region_str = ""
        if only_one_region:
            # only 1 signal region, put it on the title
            region_str = "\n"
            if do_dijet_cen:
                region_str += qgc.Dijet_CEN_LABEL
            elif do_dijet_fwd:
                region_str += qgc.Dijet_FWD_LABEL
            elif do_zpj:
                region_str += qgc.ZpJ_LABEL

        h_max = max([c.obj.GetMaximum() for c in entries])
        h_min = min([c.obj.GetMinimum(1E-10) for c in entries])
        h_range = h_max - h_min
        ylim = (max(0, h_min-(h_range*0.2)), h_max + (h_range*0.8))
        plot = Plot(entries,
                    what='hist',
                    xtitle=COMMON_STYLE_DICT['jet_pt_units_str'],
                    ytitle=angle_str,
                    title="%s jets%s" % (jet_algo['label'], region_str),
                    # ylim=(0, h_max*1.75),
                    # ylim=(h_min*0.75, h_max*1.5),
                    ylim=ylim,
                    has_data=self.has_data,
                    is_preliminary=self.is_preliminary,
                    subplot_type='ratio' if metric != 'delta' else None,
                    subplot_title='Simulation / Data',
                    subplot_limits=(0.5, 1.5) if self.has_data else (0.9, 1.1)
                    )
        # plot.default_canvas_size = (700, 600)
        plot.title_start_y = 0.85
        plot.title_left_offset = 0.05
        plot.title_font_size = 0.035
        if only_one_region:
            plot.legend.SetX1(0.6)
            plot.legend.SetX2(0.88)
            plot.legend.SetY1(0.72)
        else:
            plot.legend.SetX1(0.55)
            plot.legend.SetX2(0.78)
            plot.legend.SetY1(0.68)
        plot.legend.SetY2(0.92)
        if len(entries) > 3:
            # TODO: scale with numberof entries
            plot.legend.SetNColumns(2)
            plot.legend.SetX1(0.58)
            plot.legend.SetY1(0.6)
            plot.legend.SetX2(0.92)
            plot.legend.SetY2(0.92)
            # plot.legend.SetBorderSize(1)
            # plot.legend.SetLineColor(ROOT.kBlack)
            plot.title_left_offset = 0.03
        if len(entries) > 8:
            plot.legend.SetNColumns(3)
        plot.legend.SetY2(0.87)
        plot.left_margin = 0.16
        plot.subplot_line_style = 1
        plot.y_padding_max_linear = 1.9

        plot.do_legend = False  # do it ourselves manually
        subplot_draw_opts = "NOSTACK E1"
        plot.plot("NOSTACK E1", subplot_draw_opts)
        # plot.get_modifier().GetYaxis().SetTitleOffset(plot.get_modifier().GetYaxis().GetTitleOffset()*1.1)
        plot.set_logx(do_more_labels=False)

        # Do legend manually with graphs to get the right bars on the ends (sigh)
        for cont in dummy_entries:
            this_label = cont.label
            if '\n' not in this_label:
                plot.legend.AddEntry(cont.obj, this_label, cont.leg_draw_opt)
            else:
                # use #splitline instead of \n,
                # since the latter with > 1 columns will go in the wrong place
                parts = this_label.split("\n")
                if len(parts) > 2:
                    raise RuntimeError("Cannot handle > newlines in legend yet - requires nested #splitline")
                plot.legend.AddEntry(cont.obj, "#splitline{%s}{%s}" % (parts[0], parts[1]), cont.leg_draw_opt)

        plot.canvas.cd()
        plot.legend.Draw()

        if metric != 'delta':
            # now draw the data error shaded area
            # this is a bit hacky - basically draw them on the ratio pad,
            # then redraw the existing hists & line to get them ontop
            # note that we use "same" for all - this is to keep the original axes
            # (we may want to rethink this later?)
            plot.subplot_pad.cd()
            plot.subplot_legend = ROOT.TLegend(0.25, 0.75, 0.47, 0.9)
            plot.subplot_legend.SetFillStyle(0)
            draw_opt = "E2 SAME"
            if do_dijet_cen:
                dijet_central_hist_ratio_error.Draw(draw_opt)
                plot.subplot_legend.AddEntry(dijet_central_hist_ratio_error, "Data uncert.%s" % (" (central)" if do_dijet_fwd else ""), "F")
                if do_dijet_fwd:
                    plot.subplot_legend.SetNColumns(2)
                    plot.subplot_legend.SetX2(0.8)
            if do_dijet_fwd:
                dijet_forward_hist_ratio_error.Draw(draw_opt)
                plot.subplot_legend.AddEntry(dijet_forward_hist_ratio_error, "Data uncert.%s" % (" (forward)" if do_dijet_cen else ""), "F")
            if do_zpj:
                plot.subplot_legend.AddEntry(zpj_hist_ratio_error, "Data uncert.", "F")
                zpj_hist_ratio_error.Draw(draw_opt)
            plot.subplot_container.Draw("SAME" + subplot_draw_opts)
            plot.subplot_line.Draw()
            plot.subplot_legend.Draw()
            plot.canvas.cd()

        prefix = self._generate_filename_prefix(do_dijet_cen, do_dijet_fwd, do_zpj)
        groomed_str = '_groomed' if do_groomed else ''
        plot.save("{output_dir}/{prefix}{metric}_vs_pt_{angle_var}{groomed_str}_{algo}{append}.{fmt}"
                    .format(
                        output_dir=output_dir,
                        prefix=prefix,
                        metric=metric,
                        angle_var=angle.var,
                        groomed_str=groomed_str,
                        algo=jet_algo['name'],
                        fmt=self.output_fmt,
                        append=self.filename_append))

    def plot_dijet_zpj_means_vs_pt_all(self):
        """Plot mean vs pt for dijet (cen+fwd) & Z+jet on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_zpj_means_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('mean', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_zpj_means_vs_pt_all' % self.output_dir)

    def plot_dijet_means_vs_pt_all(self):
        """Plot mean vs pt for dijet (cen+fwd, cen, fwd) on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_means_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            # self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('mean', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_means_vs_pt_all' % self.output_dir, do_zpj=False)
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('mean', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_cen_means_vs_pt_all' % self.output_dir, do_zpj=False, do_dijet_fwd=False)
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('mean', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_fwd_means_vs_pt_all' % self.output_dir, do_zpj=False, do_dijet_cen=False)

    def plot_zpj_means_vs_pt_all(self):
        """Plot mean vs pt for zpj on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_zpj_means_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('mean', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_zpj_means_vs_pt_all' % self.output_dir, do_dijet_cen=False, do_dijet_fwd=False)

    def plot_dijet_zpj_rms_vs_pt_all(self):
        """Plot RMS vs pt for dijet (cen+fwd) & Z+jet on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_zpj_rms_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('rms', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_zpj_rms_vs_pt_all' % self.output_dir)

    def plot_dijet_rms_vs_pt_all(self):
        """Plot RMS vs pt for dijet (cen+fwd, cen, fwd) on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_rms_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            # self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('rms', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_rms_vs_pt_all' % self.output_dir, do_zpj=False)
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('rms', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_cen_rms_vs_pt_all' % self.output_dir, do_zpj=False, do_dijet_fwd=False)
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('rms', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_fwd_rms_vs_pt_all' % self.output_dir, do_zpj=False, do_dijet_cen=False)

    def plot_zpj_rms_vs_pt_all(self):
        """Plot RMS vs pt for zpj on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_zpj_rms_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('rms', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_zpj_rms_vs_pt_all' % self.output_dir, do_dijet_cen=False, do_dijet_fwd=False)

    def plot_dijet_zpj_delta_vs_pt_all(self):
        """Plot delta vs pt for dijet (cen+fwd) & Z+jet on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_zpj_delta_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('delta', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_zpj_delta_vs_pt_all' % self.output_dir)

    def plot_dijet_delta_vs_pt_all(self):
        """Plot mean vs pt for dijet (cen+fwd, cen, fwd) on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_dijet_delta_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            # self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('delta', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_delta_vs_pt_all' % self.output_dir, do_zpj=False)
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('delta', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_cen_delta_vs_pt_all' % self.output_dir, do_zpj=False, do_dijet_fwd=False)
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('delta', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_dijet_fwd_delta_vs_pt_all' % self.output_dir, do_zpj=False, do_dijet_cen=False)

    def plot_zpj_delta_vs_pt_all(self):
        """Plot mean vs pt for zpj on one plot,
        per angle, per jet algo, per groomed/ungroomed"""
        print('plot_zpj_delta_vs_pt_all...')
        for jet_algo, angle, groomed in product(self.jet_algos, self.angles, [False, True]):
            print("  ...doing", jet_algo['label'], angle.name, 'groomed' if groomed else 'ungroomed')
            self.plot_dijet_zpj_metric_vs_pt_one_angle_one_jet('delta', angle, jet_algo, do_groomed=groomed, output_dir='%s/plot_zpj_delta_vs_pt_all' % self.output_dir, do_dijet_cen=False, do_dijet_fwd=False)

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
    def _style_hist(hist, line_style=None, color=None, fill_style=None,
                    marker_size=None, marker_style=None, **kwargs):
        if line_style is not None:
            hist.SetLineStyle(line_style)
        if fill_style is not None:
            hist.SetFillStyle(fill_style)
        if color is not None:
            hist.SetLineColor(color)
            hist.SetFillColor(color)
            hist.SetMarkerColor(color)
        if marker_size is not None:
            hist.SetMarkerSize(marker_size)
        if marker_style is not None:
            hist.SetMarkerStyle(marker_style)

    def _style_data_hist(self, hist):
        self._style_hist(hist,
                         color=COMMON_STYLE_DICT['data_color'],
                         line_style=COMMON_STYLE_DICT['data_line_style'],
                         marker_size=COMMON_STYLE_DICT['marker_size'],
                         marker_style=COMMON_STYLE_DICT['data_marker_style'])

    def _style_mc_hist(self, hist):
        self._style_hist(hist,
                         line_style=COMMON_STYLE_DICT['mc_line_style'],
                         color=COMMON_STYLE_DICT['mc_color'],
                         fill_style=COMMON_STYLE_DICT['mc_fill_style'],
                         marker_size=COMMON_STYLE_DICT['marker_size'],
                         marker_style=COMMON_STYLE_DICT['mc_marker_style'])

    def _style_alt_mc_hist(self, hist):
        self._style_hist(hist,
                         line_style=COMMON_STYLE_DICT['mc_alt_line_style'],
                         color=COMMON_STYLE_DICT['mc_alt_color'],
                         fill_style=COMMON_STYLE_DICT['mc_alt_fill_style'],
                         marker_size=COMMON_STYLE_DICT['marker_size'],
                         marker_style=COMMON_STYLE_DICT['mc_alt_marker_style'])

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
        for [data, nominal MC, alt MC, other samples, ...].

        Also does styling of hists (as easier here)
        """
        mean_hists, rms_hists = [], []
        for selection_group in selections:

            mean_entries_data = []
            mean_entries_mc = []
            mean_entries_alt_mc = []
            # each top-level index is a sample
            mean_entries_other_samples = [[] for _ in self.other_samples]  # for loop needed to create indepedent lists

            rms_entries_data = []
            rms_entries_mc = []
            rms_entries_alt_mc = []
            rms_entries_other_samples = [[] for _ in self.other_samples]

            bin_names = []

            for query, label in selection_group['selections']:
                # print(query)
                results = self.df.query(query)
                if len(results.index) != 1:
                    print(query)
                    print(results)
                    raise ValueError("Got != 1 results, check query")

                mean_entries_data.append([results['mean'].item(), results['mean_err'].item()])

                rms_entries_data.append([results['rms'].item(), results['rms_err'].item()])

                if not self.only_yoda_data:
                    mean_entries_mc.append([results['mean_truth'].item(), results['mean_err_truth'].item()])
                    mean_entries_alt_mc.append([results['mean_alt_truth'].item(), results['mean_err_alt_truth'].item()])

                    rms_entries_mc.append([results['rms_truth'].item(), results['rms_err_truth'].item()])
                    rms_entries_alt_mc.append([results['rms_alt_truth'].item(), results['rms_err_alt_truth'].item()])

                for ind, sample in enumerate(self.other_samples):
                    key = sample['key']
                    mean_entries_other_samples[ind].append([results['mean_%s' % key].item(), results['mean_err_%s' % key].item()])
                    rms_entries_other_samples[ind].append([results['rms_%s' % key].item(), results['rms_err_%s' % key].item()])

                bin_names.append(label)

            # check for dodgy values
            def _check_values(values, label):
                arr = np.array(values)
                if np.any(arr == 0):
                    print("0 in", label)
                if np.all(arr == 0):
                    print("0 for all", label)
                if np.any(arr < 0):
                    print("<0 in", label)
                if np.all(arr < 0):
                    print("<0 for all", label)
                if np.isnan(arr).any():
                    print("nan in", label)
                    print(arr)
                if np.isinf(arr).any():
                    print("inf in", label)
                    print(arr)

            _check_values(mean_entries_data, "mean_entries_data")
            hist_mean_data = self._make_hist_from_values(mean_entries_data, bin_names=bin_names)
            self._style_data_hist(hist_mean_data)

            hists = [hist_mean_data]

            if not self.only_yoda_data:
                _check_values(mean_entries_mc, "mean_entries_mc")
                hist_mean_mc = self._make_hist_from_values(mean_entries_mc, bin_names=bin_names)
                self._style_mc_hist(hist_mean_mc)
                hists.append(hist_mean_mc)

                _check_values(mean_entries_alt_mc, "mean_entries_alt_mc")
                hist_mean_alt_mc = self._make_hist_from_values(mean_entries_alt_mc, bin_names=bin_names)
                self._style_alt_mc_hist(hist_mean_alt_mc)
                hists.append(hist_mean_alt_mc)

            for sample, entries in zip(self.other_samples, mean_entries_other_samples):
                _check_values(entries, "mean_entries_%s" % sample['label'])
                hist = self._make_hist_from_values(entries, bin_names=bin_names)
                self._style_hist(hist, **sample['style_dict'])
                hists.append(hist)

            mean_hists.append(hists)

            _check_values(rms_entries_data, "rms_entries_data")
            hist_rms_data = self._make_hist_from_values(rms_entries_data, bin_names=bin_names)
            self._style_data_hist(hist_rms_data)

            hists = [hist_rms_data]

            if not self.only_yoda_data:
                _check_values(rms_entries_mc, "rms_entries_mc")
                hist_rms_mc = self._make_hist_from_values(rms_entries_mc, bin_names=bin_names)
                self._style_mc_hist(hist_rms_mc)
                hists.append(hist_rms_mc)

                _check_values(rms_entries_alt_mc, "rms_entries_alt_mc")
                hist_rms_alt_mc = self._make_hist_from_values(rms_entries_alt_mc, bin_names=bin_names)
                self._style_alt_mc_hist(hist_rms_alt_mc)
                hists.append(hist_rms_alt_mc)


            for sample, entries in zip(self.other_samples, rms_entries_other_samples):
                _check_values(entries, "rms_entries_%s" % sample['label'])
                hist = self._make_hist_from_values(entries, bin_names=bin_names)
                self._style_hist(hist, **sample['style_dict'])
                hists.append(hist)

            rms_hists.append(hists)

        return mean_hists, rms_hists

    def plot_mean_rms_bins_summary(self, selections, output_file, legend_header=None):
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
        print("Plotting mean_rms_bins_summary for", legend_header)
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
        right_margin = 0.15
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

            # Draw hists - reverse order, so data last
            for ind, mhg in enumerate(mean_hist_group[::-1]):
                draw_opt = "E1"
                if ind != 0:
                    draw_opt += " SAME"
                mhg.Draw(draw_opt)

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
            # Draw hists - reverse order, so data last
            for ind, rhg in enumerate(rms_hist_group[::-1]):
                draw_opt = "E1"
                if ind != 0:
                    draw_opt += " SAME"
                rhg.Draw(draw_opt)

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
            up_padding = 0.2 * y_range
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
        # Note that this only matters if you have fill styles - if you have lines
        # you can get away with other sizing
        leg_y_top = 0.93
        # leg_pad = ROOT.TPad("leg_pad_"+cu.get_unique_str(), "", leg_x2-mean_pads[0].GetAbsWNDC(), leg_y_top-mean_pads[0].GetAbsHNDC(), leg_x2, leg_y_top)
        leg_left = mean_pads[-1].GetAbsXlowNDC() + mean_pads[-1].GetAbsWNDC() + pad_to_pad_gap
        leg_right = 1-0.02
        leg_y_bottom = leg_y_top-(1.*mean_pads[0].GetAbsHNDC())
        leg_pad = ROOT.TPad("leg_pad_"+cu.get_unique_str(), "", leg_left, leg_y_bottom, leg_right, leg_y_top)
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
        leg = ROOT.TLegend(0., mean_pads[0].GetBottomMargin(), 1, 1)

        pt = None
        if legend_header:
            # Add title to legend
            # Add ability to do multiple lines by splitting on \n
            # Assumes first line most important, so bolded
            # Dont account for blank lines
            num_header_lines = len([x for x in legend_header.split("\n") if len(x) > 0])
            line_height = 0.1
            offset = num_header_lines * line_height
            # move legend down by the height of the new TPaveText
            leg.SetY1(leg.GetY1()-offset)
            leg.SetY2(leg.GetY2()-offset)
            pt = ROOT.TPaveText(leg.GetX1(), leg.GetY2(), leg.GetX2(), leg_y_top, "NDC NB")
            pt.SetFillStyle(0)
            pt.SetBorderSize(0)
            for line_ind, line in enumerate(legend_header.split("\n")):
                text = pt.AddText(line)
                text.SetTextAlign(11)
                if line_ind == 0:
                    text.SetTextFont(62)
                    text.SetTextSize(0.1)
                else:
                    text.SetTextFont(42)
                    text.SetTextSize(0.09)
            pt.Draw()

        # Replace legend markers with graph to get correct error bar endings
        # Yes this is ridiculous
        dummy_gr = ROOT.TGraphErrors(1, array('d', [1]), array('d', [1]), array('d', [1]), array('d', [1]))
        dummy_data = dummy_gr.Clone()
        self._style_data_hist(dummy_data)
        dummy_mc = dummy_gr.Clone()
        self._style_mc_hist(dummy_mc)
        dummy_alt_mc = dummy_gr.Clone()
        self._style_alt_mc_hist(dummy_alt_mc)

        def _add_entry(obj, label, option):
            parts = label.split("\n")
            if len(parts) == 1:
                leg.AddEntry(obj, label, option)
            elif len(parts) == 2:
                leg.AddEntry(obj, "#splitline{%s}{%s}" % (parts[0], parts[1]), option)
            else:
                raise RuntimeError("Cannot handle > 1 \n in legend label")

            # for line_ind, line in enumerate(label.split("\n")):
            #     if line_ind == 0:
            #         leg.AddEntry(obj, line, option)
            #     else:
            #         leg.AddEntry(0, line, "")

        _add_entry(dummy_data, "Data", "P")
        if not self.only_yoda_data:
            _add_entry(dummy_mc, self.mc_label, "P")
            _add_entry(dummy_alt_mc, self.alt_mc_label, "P")

        dummy_other = []  # keep reference
        for sample in self.other_samples:
            dummy_sample = dummy_gr.Clone()
            self._style_hist(dummy_sample, **sample['style_dict'])
            dummy_other.append(dummy_sample)
            _add_entry(dummy_sample, sample['style_dict'].get('label', sample['label']), "P")

        # Add a dummy entry, otherwise it won't print the label of the last entry
        # No idea why - seems correlated with having > 2 lines in the legend header?
        # Absolute mess
        leg.AddEntry(0, "", "")
        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        # leg.SetFillColor(ROOT.kGreen)
        # leg.SetFillStyle(3004)
        leg.SetTextSize(0.1)
        leg.SetTextAlign(12)
        leg.SetEntrySeparation(0.08)
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

                # FIXME: add other_samples

                bin_names.append(label)

            hist_delta_nominal = self._make_hist_from_values(entries_nominal, bin_names=bin_names)
            self._style_mc_hist(hist_delta_nominal)

            hist_delta_alt = self._make_hist_from_values(entries_alt, bin_names=bin_names)
            self._style_alt_mc_hist(hist_delta_alt)

            delta_hists.append([hist_delta_nominal, hist_delta_alt])

        return delta_hists


    def plot_delta_bins_summary(self, selections, output_file, legend_header=None):
        """"""
        print("Plotting delta_bins_summary for", legend_header)
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
        right_margin = 0.15
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
            # FIXME add other_samples
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
        leg_left = pads[-1].GetAbsXlowNDC() + pads[-1].GetAbsWNDC() + pad_to_pad_gap
        leg_right = 1-0.02
        leg_y_bottom = leg_y_top-(1.*pads[0].GetAbsHNDC())
        leg_pad = ROOT.TPad("leg_pad_"+cu.get_unique_str(), "", leg_left, leg_y_bottom, leg_right, leg_y_top)
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
        leg = ROOT.TLegend(0., pads[0].GetBottomMargin(), 1, 1)

        pt = None
        if legend_header:
            # Add title to legend
            # Add ability to do multiple lines by splitting on \n
            # Assumes first line most important, so bolded
            # Dont account for blank lines
            num_header_lines = len([x for x in legend_header.split("\n") if len(x) > 0])
            line_height = 0.1
            offset = num_header_lines * line_height
            # move legend down by the height of the new TPaveText
            leg.SetY1(leg.GetY1()-offset)
            leg.SetY2(leg.GetY2()-offset)
            pt = ROOT.TPaveText(leg.GetX1(), leg.GetY2(), leg.GetX2(), leg_y_top, "NDC NB")
            pt.SetFillStyle(0)
            pt.SetBorderSize(0)
            for line_ind, line in enumerate(legend_header.split("\n")):
                text = pt.AddText(line)
                text.SetTextAlign(11)
                if line_ind == 0:
                    text.SetTextFont(62)
                    text.SetTextSize(0.1)
                else:
                    text.SetTextFont(42)
                    text.SetTextSize(0.09)
            pt.Draw()

        # Replace legend markers with graph to get correct error bar endings
        # Yes this is ridiculous
        dummy_gr = ROOT.TGraphErrors(1, array('d', [1]), array('d', [1]), array('d', [1]), array('d', [1]))
        dummy_mc = dummy_gr.Clone()
        self._style_mc_hist(dummy_mc)
        dummy_alt_mc = dummy_gr.Clone()
        self._style_alt_mc_hist(dummy_alt_mc)
        leg.AddEntry(dummy_mc, self.mc_label, "EL")
        leg.AddEntry(dummy_alt_mc, self.alt_mc_label, "EL")
        # Add a dummy entry, otherwise it won't print the label of the last entry
        # No idea why - seems correlated with having > 2 lines in the legend header?
        # Absolute mess
        leg.AddEntry(0, "", "")
        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        # leg.SetFillColor(ROOT.kGreen)
        # leg.SetFillStyle(3004)
        leg.SetTextSize(0.1)
        leg.SetTextAlign(12)
        leg.SetEntrySeparation(0.08)
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
    parser.add_argument("--h5input",
                        help="Read analysis data from H5 input file (from running extract_unfolding_summary_stats.py)")
    parser.add_argument("--h5inputRivet",
                        help="Read RIVET data from H5 input file (from running extract_rivet_summary_stats.py)")
    parser.add_argument("--outputDir",
                        default=None,
                        help='Output directory (default is the source dir)')
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
    parser.add_argument("--onlyYodaData",
                        action='store_true',
                        help='Only plot data & Yoda inputs (ignore MG+Py, H++)')
    args = parser.parse_args()

    # Get input data
    if not any([args.h5input, args.h5inputRivet, args.yodaInputDijet, args.yodaInputZPJ]):
        raise RuntimeError("Need one of --h5input, --h5inputRivet, --yodaInputDijet, or --yodaInputZPJ")

    if (len(args.yodaInputDijet) != len(args.yodaLabel)
        and len(args.yodaInputZPJ) != len(args.yodaLabel)):
        raise RuntimeError("Number of --yodaInputDijet/yodaInputZPJ must match number of --yodaLabel")

    if not args.outputDir:
        if args.h5input:
            args.outputDir = os.path.dirname(os.path.abspath(args.h5input))
        else:
            args.outputDir = os.getcwd()
        print("Setting output dir to", args.outputDir)

    # ----------------------------------------------------------------------
    # Read in data from h5 file
    # -----------------------------------------------------------------------
    print("Reading in unfolding data from existing HDF5 file...")
    with pd.HDFStore(args.h5input) as store:
        df = store['df']
    print(df.head())
    print("# entries:", len(df.index))

    if args.h5inputRivet:
        print("Reading in RIVET data from existing HDF5 file...")
        with pd.HDFStore(args.h5inputRivet) as store:
            df_rivet = store['df']
        print(df_rivet.head())
        print("# Rivet entries:", len(df_rivet.index))
        # Figure out YODA entries from column names
        mean_columns = [c.replace("mean_err_", '') for c in df_rivet.columns if 'mean_err_' in c]
        print(mean_columns)
        args.yodaLabel = mean_columns
        df = pd.merge(df, df_rivet, how='outer')

    # -----------------------------------------------------------------------
    # Get stats from YODA files, add to dataframe
    # -----------------------------------------------------------------------
    if len(args.yodaInputDijet) > 0:
        df_rivet = get_dataframe_from_yoda_inputs(zip(args.yodaInputDijet, args.yodaInputZPJ, args.yodaLabel))
        df = pd.merge(df, df_rivet, how='outer')

    convert_df_types(df)
    print(df.columns)
    print(df.head())
    print(df.dtypes)
    print(len(df.index), 'entries in dataframe')

    # Filter only regions/algos/angles in the dataframe, since it could have
    # been modified earlier
    all_jet_algos = [
        {'label': 'AK4 PUPPI', 'name': 'ak4puppi'},
        {'label': 'AK8 PUPPI', 'name': 'ak8puppi'}
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
                             has_data=True,
                             only_yoda_data=args.onlyYodaData)

    has_dijet = any(["Dijet" in r['name'] for r in regions])
    has_zpj = any(["ZPlusJet" in r['name'] for r in regions])

    # Add extra samples
    for ylabel in args.yodaLabel:
        if ylabel not in SAMPLE_STYLE_DICTS:
            print("No entry found in SAMPLE_STYLE_DICTS for", ylabel, ", using defaults")
        plotter.add_sample(key=dataframe_yoda_key(ylabel),
                           label=ylabel,
                           style_dict=SAMPLE_STYLE_DICTS.get(ylabel, dict()))

    # For summary plots
    low_pt = 120
    high_pt = 614
    high_pt = 800
    high_pt = 1000

    ak4_str = "AK4"
    ak8_str = "AK8"

    DO_X_VS_PT_PLOTS = True
    DO_OVERALL_SUMMARY_PLOTS = True

    filename_append = plotter.filename_append

    if has_dijet:
        if DO_X_VS_PT_PLOTS:
            plotter.plot_dijet_means_vs_pt_all()
            plotter.plot_dijet_rms_vs_pt_all()
            plotter.plot_dijet_delta_vs_pt_all()

        pt_bins = qgc.PT_UNFOLD_DICT['signal_gen']
        low_pt_bin = np.where(pt_bins == low_pt)[0][0]
        low_pt_str = "[%g, %g] GeV" % (low_pt, pt_bins[low_pt_bin+1])

        high_pt_bin = np.where(pt_bins == high_pt)[0][0]
        high_pt_str = "[%g, %g] GeV" % (high_pt, pt_bins[high_pt_bin+1])

        if DO_OVERALL_SUMMARY_PLOTS:
            # DIJET CENTRAL
            # ---------------------------------------------

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

                    ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_central" & angle=="%s_charged"' % (low_pt_bin, angle.var),
                        "#splitline{{{jet_str}, {pt_str}}}{{Charged-only}}".format(jet_str=ak4_str, pt_str=low_pt_str)),
                        # "#splitline{{#splitline{{{jet_str}}}{{{pt_str}}}}}{{Charged-only}}".format(jet_str=ak8_str, pt_str=low_pt_str)),

                    ('jet_algo=="ak4puppi" & pt_bin==%d & isgroomed  & region=="Dijet_central_groomed" & angle=="%s"' % (low_pt_bin, angle.var),
                        "#splitline{{{jet_str}, {pt_str}}}{{Groomed}}".format(jet_str=ak4_str, pt_str=low_pt_str)),
                ]
                selections.append({'label': this_angle_str, 'selections': this_selection})

            legend_header = "Dijet (central) region"
            plotter.plot_mean_rms_bins_summary(
                selections=selections,
                legend_header=legend_header,
                output_file=os.path.join(args.outputDir, "dijet_central_mean_rms_summary%s.pdf" % (filename_append))
            )

            legend_header = "Dijet (central) region"
            plotter.plot_delta_bins_summary(
                selections=selections,
                legend_header=legend_header,
                output_file=os.path.join(args.outputDir, "dijet_central_delta_summary%s.pdf" % (filename_append))
            )

            legend_header = "Gluon-enriched jets\nDijet (central) region"
            plotter.plot_mean_rms_bins_summary(
                selections=selections,
                legend_header=legend_header,
                output_file=os.path.join(args.outputDir, "gluon_mean_rms_summary%s.pdf" % (filename_append))
            )

            legend_header = "Gluon-enriched jets\nDijet (central) region"
            plotter.plot_delta_bins_summary(
                selections=selections,
                legend_header=legend_header,
                output_file=os.path.join(args.outputDir, "gluon_delta_summary%s.pdf" % (filename_append))
            )

            # DIJET FORWARD
            # ---------------------------------------------

            # Create selection queries for mean/RMS/delta summary plots
            selections = []
            for angle in charged_and_neutral_angles:
                this_angle_str = "%s (%s)" % (angle.name, angle.lambda_str)
                # this_angle_str = "%s" % (angle.lambda_str)
                this_selection = [
                    ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_forward" & angle=="%s"' % (low_pt_bin, angle.var),
                        "{jet_str}, {pt_str}".format(jet_str=ak4_str, pt_str=low_pt_str)),
                        # "#splitline{{{jet_str}}}{{{pt_str}}}".format(jet_str=ak4_str, pt_str=low_pt_str)),

                    ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_forward" & angle=="%s"' % (high_pt_bin, angle.var),
                        "{jet_str}, {pt_str}".format(jet_str=ak4_str, pt_str=high_pt_str)),
                        # "#splitline{{{jet_str}}}{{{pt_str}}}".format(jet_str=ak4_str, pt_str=high_pt_str)),

                    ('jet_algo=="ak8puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_forward" & angle=="%s"' % (low_pt_bin, angle.var),
                        "{jet_str}, {pt_str}".format(jet_str=ak8_str, pt_str=low_pt_str)),
                        # "#splitline{{{jet_str}}}{{{pt_str}}}".format(jet_str=ak8_str, pt_str=low_pt_str)),

                    ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_forward" & angle=="%s_charged"' % (low_pt_bin, angle.var),
                        "#splitline{{{jet_str}, {pt_str}}}{{Charged-only}}".format(jet_str=ak4_str, pt_str=low_pt_str)),
                        # "#splitline{{#splitline{{{jet_str}}}{{{pt_str}}}}}{{Charged-only}}".format(jet_str=ak8_str, pt_str=low_pt_str)),

                    ('jet_algo=="ak4puppi" & pt_bin==%d & isgroomed  & region=="Dijet_forward_groomed" & angle=="%s"' % (low_pt_bin, angle.var),
                        "#splitline{{{jet_str}, {pt_str}}}{{Groomed}}".format(jet_str=ak4_str, pt_str=low_pt_str)),
                ]
                selections.append({'label': this_angle_str, 'selections': this_selection})

            legend_header = "Dijet (forward) region"
            plotter.plot_mean_rms_bins_summary(
                selections=selections,
                legend_header=legend_header,
                output_file=os.path.join(args.outputDir, "dijet_forward_mean_rms_summary%s.pdf" % (filename_append))
            )

            plotter.plot_delta_bins_summary(
                selections=selections,
                legend_header=legend_header,
                output_file=os.path.join(args.outputDir, "dijet_forward_delta_summary%s.pdf" % (filename_append))
            )

    if has_zpj:
        if DO_X_VS_PT_PLOTS:
            plotter.plot_zpj_means_vs_pt_all()
            plotter.plot_zpj_rms_vs_pt_all()
            plotter.plot_zpj_delta_vs_pt_all()

        pt_bins = qgc.PT_UNFOLD_DICT['signal_zpj_gen']
        low_pt_bin = np.where(pt_bins == low_pt)[0][0]
        low_pt_str = "[%g, %g] GeV" % (low_pt, pt_bins[low_pt_bin+1])

        high_pt = 326
        high_pt_bin = np.where(pt_bins == high_pt)[0][0]
        high_pt_str = "[%g, %g] GeV" % (high_pt, pt_bins[high_pt_bin+1])

        if DO_OVERALL_SUMMARY_PLOTS:
            selections = []
            for angle in charged_and_neutral_angles:
                this_angle_str = "%s (%s)" % (angle.name, angle.lambda_str)
                # this_angle_str = "%s" % (angle.lambda_str)
                this_selection = [
                    ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s"' % (low_pt_bin, angle.var),
                        "{jet_str}, {pt_str}".format(jet_str=ak4_str, pt_str=low_pt_str)),
                        # "#splitline{{{jet_str}}}{{{pt_str}}}".format(jet_str=ak4_str, pt_str=low_pt_str)),

                    # ignore high pt bin as not useful - same composition as dijet but fewer stats
                    # ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s"' % (high_pt_bin, angle.var),
                    #     "{jet_str}, {pt_str}".format(jet_str=ak4_str, pt_str=high_pt_str)),
                    #     # "#splitline{{{jet_str}}}{{{pt_str}}}".format(jet_str=ak4_str, pt_str=high_pt_str)),

                    ('jet_algo=="ak8puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s"' % (low_pt_bin, angle.var),
                        "{jet_str}, {pt_str}".format(jet_str=ak8_str, pt_str=low_pt_str)),
                        # "#splitline{{{jet_str}}}{{{pt_str}}}".format(jet_str=ak8_str, pt_str=low_pt_str)),

                    ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s_charged"' % (low_pt_bin, angle.var),
                        "#splitline{{{jet_str}, {pt_str}}}{{Charged-only}}".format(jet_str=ak4_str, pt_str=low_pt_str)),
                        # "#splitline{{#splitline{{{jet_str}}}{{{pt_str}}}}}{{Charged-only}}".format(jet_str=ak8_str, pt_str=low_pt_str)),

                    ('jet_algo=="ak4puppi" & pt_bin==%d & isgroomed  & region=="ZPlusJets_groomed" & angle=="%s"' % (low_pt_bin, angle.var),
                        "#splitline{{{jet_str}, {pt_str}}}{{Groomed}}".format(jet_str=ak4_str, pt_str=low_pt_str)),
                ]
                selections.append({'label': this_angle_str, 'selections': this_selection})

            legend_header = "Z+jets region"

            plotter.plot_mean_rms_bins_summary(
                selections=selections,
                legend_header=legend_header,
                output_file=os.path.join(args.outputDir, "zpj_mean_rms_summary%s.pdf" % (filename_append))
            )

            plotter.plot_delta_bins_summary(
                selections=selections,
                legend_header=legend_header,
                output_file=os.path.join(args.outputDir, "zpj_delta_summary%s.pdf" % (filename_append))
            )

    if has_dijet and has_zpj and DO_OVERALL_SUMMARY_PLOTS:
        low_pt = 120
        pt_bins = qgc.PT_UNFOLD_DICT['signal_gen']
        low_pt_bin = np.where(pt_bins == low_pt)[0][0]
        low_pt_str = "[%g, %g] GeV" % (low_pt, pt_bins[low_pt_bin+1])

        high_pt = 1000
        high_pt_bin = np.where(pt_bins == high_pt)[0][0]
        high_pt_str = "[%g, %g] GeV" % (high_pt, pt_bins[high_pt_bin+1])

        # DIJET CENTRAL
        # ---------------------------------------------

        # Create selection queries for mean/RMS/delta summary plots
        selections = []
        for angle in charged_and_neutral_angles:
            this_angle_str = "%s (%s)" % (angle.name, angle.lambda_str)
            # this_angle_str = "%s" % (angle.lambda_str)
            this_selection = [
                ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s"' % (low_pt_bin, angle.var),
                    "{jet_str}, {pt_str}".format(jet_str=ak4_str, pt_str=low_pt_str)),
                    # "#splitline{{{jet_str}}}{{{pt_str}}}".format(jet_str=ak4_str, pt_str=low_pt_str)),

                ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="Dijet_forward" & angle=="%s"' % (high_pt_bin, angle.var),
                    "{jet_str}, {pt_str}".format(jet_str=ak4_str, pt_str=high_pt_str)),
                    # "#splitline{{{jet_str}}}{{{pt_str}}}".format(jet_str=ak4_str, pt_str=high_pt_str)),

                ('jet_algo=="ak8puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s"' % (low_pt_bin, angle.var),
                    "{jet_str}, {pt_str}".format(jet_str=ak8_str, pt_str=low_pt_str)),
                    # "#splitline{{{jet_str}}}{{{pt_str}}}".format(jet_str=ak8_str, pt_str=low_pt_str)),

                ('jet_algo=="ak4puppi" & pt_bin==%d & ~isgroomed & region=="ZPlusJets" & angle=="%s_charged"' % (low_pt_bin, angle.var),
                    "#splitline{{{jet_str}, {pt_str}}}{{Charged-only}}".format(jet_str=ak4_str, pt_str=low_pt_str)),
                    # "#splitline{{#splitline{{{jet_str}}}{{{pt_str}}}}}{{Charged-only}}".format(jet_str=ak8_str, pt_str=low_pt_str)),

                ('jet_algo=="ak4puppi" & pt_bin==%d & isgroomed  & region=="ZPlusJets_groomed" & angle=="%s"' % (low_pt_bin, angle.var),
                    "#splitline{{{jet_str}, {pt_str}}}{{Groomed}}".format(jet_str=ak4_str, pt_str=low_pt_str)),
            ]
            selections.append({'label': this_angle_str, 'selections': this_selection})

        legend_header = "Quark-enriched jets\n%s:\nZ+jets region\n%s:\nDijet (forward) region" % (low_pt_str, high_pt_str)
        plotter.plot_mean_rms_bins_summary(
            selections=selections,
            legend_header=legend_header,
            output_file=os.path.join(args.outputDir, "quark_mean_rms_summary%s.pdf" % (filename_append))
        )
        plotter.plot_delta_bins_summary(
            selections=selections,
            legend_header=legend_header,
            output_file=os.path.join(args.outputDir, "quark_delta_summary%s.pdf" % (filename_append))
        )
