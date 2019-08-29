#!/usr/bin/env python

"""Do lambda corrections"""


import argparse
from MyStyle import My_Style
My_Style.cd()
import os
from itertools import product
from array import array
from math import sqrt
import bisect
import numpy as np
from collections import OrderedDict


import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
# ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)


# My stuff
import qg_general_plots as qgg
import qg_common as qgc
from comparator import Contribution, Plot
import common_utils as cu

ROOT.gStyle.SetPaintTextFormat(".3f")

# Control output format
OUTPUT_FMT = "pdf"


def get_hist_metric(h_reco, metric="rawmean"):
    metric = metric.lower()
    options = ["rawmean", "median", "gausfit"]
    if metric not in options:
        raise RuntimeError("metric must be one of: %s" % (", ".join(options)))
    if metric == "rawmean":
        return h_reco.GetMean(), h_reco.GetMeanError()
    elif metric == "median":
        quantiles = array('d', [0.25, 0.5, 0.75])
        values = array('d', [0, 0, 0])
        h_reco.GetQuantiles(len(quantiles), values, quantiles)
        quantile_ideal_lower_edge, median, quantile_ideal_upper_edge = values
        return median, 1.253 * h_reco.GetMeanError() / sqrt(h_reco.GetEffectiveEntries())
    elif metric == "gausfit":
        return 1, 1


def get_corrections(h2d_response, h2d_rel_response, metric):
    """"""
    # assumes x axis is GEN
    gen_bins = cu.get_bin_edges(h2d_response, 'x')
    x_values, corrections = [], []
    reco_hists, response_hists = [], []
    for bin_ind, (bin_low, bin_high) in enumerate(zip(gen_bins[:-1], gen_bins[1:])):
        print(bin_ind, ":", bin_low, bin_high)

        h_reco = qgg.get_projection_plot(h2d_response, bin_low, bin_high, cut_axis='x')
        h_reco.SetName("variable: %g - %g" % (bin_low, bin_high))
        h_reco.SetTitle("variable: %g - %g" % (bin_low, bin_high))

        h_response = qgg.get_projection_plot(h2d_rel_response, bin_low, bin_high, cut_axis='x')
        h_response.SetName("variable: %g - %g" % (bin_low, bin_high))
        h_response.SetTitle("variable: %g - %g" % (bin_low, bin_high))

        if h_reco.GetEntries() < 10 or h_response.GetEntries() < 10:
            continue
        reco_value, reco_err = get_hist_metric(h_reco, "rawmean")
        response_value, response_err = get_hist_metric(h_response, metric)
        if response_value > 0:
            corr = 1./response_value
            corr_err = response_err / (response_value ** 2)
        else:
            raise RuntimeError("response_value for %g-%g is <= 0" % (bin_low, bin_high))
        x_values.append([reco_value, reco_err])
        corrections.append([corr, corr_err])
        reco_hists.append(h_reco)
        response_hists.append(h_response)

    return np.array(x_values), np.array(corrections), reco_hists, response_hists


def plot_corrections(corrections_graph, xtitle, title, output_filename):
    # print(x_values)
    # print(corrections)
    conts = [Contribution(corrections_graph, label="Correction", 
                          line_color=ROOT.kRed, marker_color=ROOT.kRed, 
                          line_width=2, marker_size=1, marker_style=20, 
                          fit_match_style=False)
    ]
    # plot = Plot(conts, what='graph', xtitle=xtitle, ytitle="Correction", title=title, has_data=False, ylim=[0, 2.5])
    plot = Plot(conts, what='both', xtitle=xtitle, ytitle="Correction", title=title, has_data=False, ylim=[0, 2.5])
    plot.plot("ALP")
    plot.save(output_filename)


def plot_response_hists(response_hists, output_dir):
    cu.check_dir_exists_create(output_dir)
    for ind, h in enumerate(response_hists):
        c = ROOT.TCanvas(cu.get_unique_str(), "", 600, 800)
        h.Draw("HISTE")
        c.SaveAs(os.path.join(output_dir, "response_hist_%d.%s" % (ind, OUTPUT_FMT)))


def fit_corrections(corrections_graph, start_fit, end_fit, fitfunc):
    # fit to corrections graph
    
    fit_func = None
    if fitfunc == "linear":
        fit_func = ROOT.TF1("fit_func"+cu.get_unique_str(), "[0]+[1]*x", start_fit, end_fit)
        fit_func.SetParameter(0, 1.2)
        fit_func.SetParameter(1, 0.005)

    elif fitfunc == "log":
        fit_func = ROOT.TF1("fit_func"+cu.get_unique_str(), "([1]*log((x-[0])/[2]))", start_fit, end_fit)
        fit_func.SetParameter(0, 1)
        fit_func.SetParameter(1, 1)
        fit_func.SetParameter(2, 1)
    else:
        raise RuntimeError("Don't know how to handle fitfunc = %s" % fitfunc)
    
    fit_result = corrections_graph.Fit(fit_func, "VRS", "", start_fit, end_fit)
    return fit_result

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        help='Input ROOT file to process.')
    parser.add_argument("--outputDir",
                        help="Directory to put output plot dirs into",
                        default=None)
    parser.add_argument("--outputFile",
                        help="Output ROOT file for things",
                        default=None)
    parser.add_argument("--fitFunc",
                        help="Fit function to use",
                        choices=["linear", "log"])

    args = parser.parse_args()

    input_dir, input_basename = os.path.split(args.input)
    default_plot_dir = os.path.join(input_dir, "rebinning_"+os.path.splitext(input_basename)[0])
    plot_dir = args.outputDir if args.outputDir else default_plot_dir

    # default_output_root_filename = os.path.join(input_dir, "rebinned_" + input_basename)
    # output_root_filename = args.outputFile if args.outputFile else default_output_root_filename
    # cu.check_dir_exists_create(os.path.dirname(os.path.abspath(output_root_filename)))
    # output_tfile = cu.open_root_file(output_root_filename, 'RECREATE')

    source_plot_dir_name = None
    region_label = None
    if "qcd" in args.input.lower():
        source_plot_dir_name = "Dijet_QG_tighter"
        region_label = qgc.QCD_Dijet_LABEL
    elif "dyjetstoll" in args.input.lower():
        source_plot_dir_name = "ZPlusJets_QG"
        region_label = qgc.DY_ZpJ_LABEL
    else:
        raise RuntimeError("No idea which region we're using")

    pt_regions = [
        # {
        #     "append": "",
        #     "title": "p_{T}^{Reco} > 30 GeV",
        # },
        # {
        #     "append": "_lowPt",
        #     "title": "30 < p_{T}^{Reco} < 100 GeV",
        # },
        {
            "append": "_midPt",
            "title": "100 < p_{T}^{Reco} < 250 GeV",
        },
        {
            "append": "_highPt",
            "title": "p_{T}^{Reco} > 250 GeV",
        },
    ]

    input_tfile = cu.open_root_file(args.input)

    for angle in qgc.COMMON_VARS[:]: # only care about multiplicity

        for pt_region_dict in pt_regions[:]:

            var_dict = {
                "name": "%s/%s%s" % (source_plot_dir_name, angle.var, pt_region_dict['append']),
                "var_label": "%s (%s)" % (angle.name, angle.lambda_str),
                "title": "%s\n%s" % (region_label, pt_region_dict['title']),
            }

            # Need both to do corrections
            # this is reco vs gen
            full_var_name = var_dict['name'] + "_response"
            # relative response is reco/gen vs gen
            full_var_name_rel = var_dict['name'] + "_rel_response"

            h2d_orig = cu.get_from_tfile(input_tfile, full_var_name)
            h2d_rel_orig = cu.get_from_tfile(input_tfile, full_var_name_rel)

            h2d_orig.Rebin2D(4, 4)
            h2d_rel_orig.Rebin2D(4, 4)

            metric = "rawmean"
            metric = "median"
            x_values, corrections, reco_hists, response_hists = get_corrections(h2d_orig, h2d_rel_orig, metric)

            corrections_graph = ROOT.TGraphErrors(len(x_values),
                                                  array('d', x_values[:,0]), 
                                                  array('d', corrections[:,0]), 
                                                  array('d', x_values[:,1]), 
                                                  array('d', corrections[:,1]))

            start_fit, end_fit = 0, 1000
            # fit to corrections
            if args.fitFunc == "linear":
                start_fit, end_fit = 20, 90
                start_fit = max(start_fit, x_values[:,0].min())
                end_fit = min(end_fit, x_values[:,0].max())
            
            elif args.fitFunc == "log":
                start_fit = x_values[3,0]
                end_fit = x_values[-3,0]
            # print(start_fit, end_fit)
            
            fit_result = fit_corrections(corrections_graph, start_fit, end_fit, fitfunc=args.fitFunc)

            output_filename = "%s/%s_corrections_%s.%s" % (args.outputDir, var_dict['name'], metric, OUTPUT_FMT)
            plot_corrections(corrections_graph,
                             xtitle='RECO '+var_dict['var_label'], 
                             title=var_dict['title'],
                             output_filename=output_filename)
            
            plot_response_hists(response_hists, 
                                output_dir=os.path.join(args.outputDir, "%s_correction_hists" % var_dict['name']))

    # output_tfile.Close()
    input_tfile.Close()

