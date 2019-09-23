#!/usr/bin/env python

"""Figure out optimial binning for lambda variables.

Currently uses minimum required purity & stability in response matrix.

Can also rebin other inputs according to the binning scheme from the main input.
"""


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
from copy import deepcopy
os.nice(10)

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


def calc_variable_binning(h2d, plot_dir, metric):
    metric = metric.lower()
    this_h2d = h2d.Clone(h2d.GetName()+"Clone")
    # this_h2d.Rebin2D(2, 2)
    reco_bin_edges = cu.get_bin_edges(this_h2d, 'Y')
    gen_bin_edges = cu.get_bin_edges(this_h2d, 'X')
    print(reco_bin_edges)

    # For first reco bin, get 1D hist of gen level variable
    bin_start, bin_end = 0, 1
    lower_edge = reco_bin_edges[bin_start]
    bins = []
    counter = 0
    while bin_end < len(reco_bin_edges) and counter < 100:
        print("***ITERATION", counter)
        print("bin_start, bin_end:", bin_start, bin_end)

        if bin_end == len(reco_bin_edges)-1:
            # got to the end, add one last bin
            bins.append((reco_bin_edges[bin_start], reco_bin_edges[-1]))
            break

        hproj = this_h2d.ProjectionX(h2d.GetName()+"_%d_%d" % (bin_start, bin_end), bin_start, bin_end, "eo")

        if hproj.GetEntries() < 100:
            # assume we're at the high end of the scale
            if bin_start > 5:
                print("WARNING: Too few entries, adding one final wide bin & quitting")
                bins.append((reco_bin_edges[bin_start], reco_bin_edges[-1]))
                break
            else:
                # assume we're at the low end of the scale
                print("WARNING: Too few entries, increasing bin width")
                bin_end += 1
                counter += 1
                continue

        # Rough way to find width of distribution
        mean = hproj.GetMean()
        width = sqrt(hproj.GetRMS())

        # if metric == "gausfit":
        # Fit the distribution to actually find the width
        peak = hproj.GetBinCenter(hproj.GetMaximumBin())
        fit_mod = 1
        fgaus = ROOT.TF1("fgaus", "gaus", peak-fit_mod*hproj.GetRMS(), peak+fit_mod*hproj.GetRMS());
        fgaus.SetNpx(1000)
        fgaus.SetParameter('mean', peak)
        fgaus.SetParameter('sigma', width)
        fgaus.SetLineColor(ROOT.kRed)
        fgaus.SetLineWidth(1)
        fit_res = hproj.Fit("fgaus", "QSR")
        if fit_res.Status() == 0:
            mean = fit_res.Parameter(1)
            width = fit_res.Parameter(2)

        print("Mean: %.3f" % mean, "width: %.3f" % width, "[fit]" if fit_res.Status() == 0 else "[raw]")

        # if bin_end == len(reco_bin_edges)-1:
        #     # deal with last bin
        #     bins.append((reco_bin_edges[bin_start], reco_bin_edges[-1]))

        if mean < reco_bin_edges[bin_start]:
            if bin_end == len(reco_bin_edges)-1:
                # deal with last bin
                bins.append((reco_bin_edges[bin_start], reco_bin_edges[-1]))
            else:
                print("Mean < lower edge, so let's make bin wider", bin_start, bin_end)
                bin_end += 1
                counter += 1
                continue

        # Ideal bin edges
        factor = 0.85
        ideal_lower_edge = max(0, mean - (factor * width))
        ideal_upper_edge = mean + (factor * width)

        if metric == "quantile":
            desired_width = 0.7
            quantiles = array('d', [(1-desired_width)/2., 1-((1-desired_width)), 0.5])
            values = array('d', [0, 0, 0])
            hproj.GetQuantiles(len(quantiles), values, quantiles)
            quantile_ideal_lower_edge, quantile_ideal_upper_edge, median = values
            # ideal_upper_edge = max(quantile_ideal_upper_edge, ideal_upper_edge)
            bin_mid = 0.5*(reco_bin_edges[bin_start]+reco_bin_edges[bin_start])
            quantile_width = values[1] - values[0]
            # ideal_upper_edge =  bin_mid + (quantile_width * bin_mid / median)
            ideal_upper_edge =  reco_bin_edges[bin_start] + (quantile_width * bin_mid / median)

        print("Ideal bin:", ideal_lower_edge, ideal_upper_edge)

        if bin_start == 0 and ideal_lower_edge > 0:
            # We need an extra bin to accommodate the region between 0 and the low edge
            # print("bin_start", bin_start, "bin_end", bin_end, "mean", mean, "width", width, "lower_edge", ideal_lower_edge)
            # print("WARNING, lower_edge > 0, adding an extra bin")

            # new_bin_end = bisect.bisect_right(reco_bin_edges, ideal_lower_edge)

            # Absorb the region 0 - ideal_lower_edge into this bin
            new_bin_end = bisect.bisect_right(reco_bin_edges, ideal_upper_edge)
            upper_edge = reco_bin_edges[new_bin_end]
            bins.append((0, upper_edge))
            print("Adding bin", bins[-1])

            c = ROOT.TCanvas("c"+cu.get_unique_str(), "", 800, 600)
            hproj.Draw()
            output_filename = os.path.join(plot_dir, h2d.GetName() + "_bin%dto%d.%s" % (bin_start, new_bin_end, OUTPUT_FMT))
            c.SaveAs(output_filename)

            bin_start = new_bin_end
            bin_end = bin_start + 1
        else:
            # Just use the upper edge from previous bin as lower edge
            # TODO check this?
            lower_edge = reco_bin_edges[bin_start]
            new_bin_end = bisect.bisect_right(reco_bin_edges, ideal_upper_edge)

            if new_bin_end == bin_start:
                # raise RuntimeError("new_bin_end == bin_start, got stuck")
                print("new_bin_end == bin_start, got stuck, increasing by 1")
                new_bin_end += 1

            if new_bin_end >= len(reco_bin_edges):
                new_bin_end = len(reco_bin_edges)-1

            upper_edge = reco_bin_edges[new_bin_end]
            bins.append((lower_edge, upper_edge))
            print("Adding bin", bins[-1])

            c = ROOT.TCanvas("c"+cu.get_unique_str(), "", 800, 600)
            hproj.Draw()
            output_filename = os.path.join(plot_dir, h2d.GetName() + "_bin%dto%d.%s" % (bin_start, new_bin_end, OUTPUT_FMT))
            c.SaveAs(output_filename)

            if (upper_edge == reco_bin_edges[-1]):
                break

            bin_start = new_bin_end
            bin_end = bin_start + 1

            print("new bin_start:", bin_start)

            counter += 1

    print("BINS:", bins)
    return bins


def renorm(arr2d, axis):
    # create version where each axis summed to 1
    # use where and out args to ensure nans are made into 0s
    summed = arr2d.sum(axis=axis, keepdims=True)
    return np.divide(arr2d, summed, where=summed!=0, out=np.zeros_like(arr2d))


def concat_row(arr2d, row_ind):
    # concat row row_ind + row_ind+1
    nrows, ncols = arr2d.shape
    if row_ind > nrows - 2:
        raise RuntimeError("Cannot concat row [%d] as only %d rows in matrix" % (row_ind, nrows))
    arr2d_new = np.zeros(shape=(nrows-1, ncols), dtype=float)
    new_row = arr2d[row_ind] + arr2d[row_ind+1]
    arr2d_new[row_ind] = new_row
    # fill in new matrix
    if row_ind > 0:
        # do that bit before the new row
        arr2d_new[:row_ind, ] = arr2d[:row_ind, ]
    if row_ind < nrows - 2:
        arr2d_new[row_ind+1:, :] = arr2d[row_ind+2:, :]
    return arr2d_new


def calc_variable_binning_purity_stability(h2d, purity_goal=0.4, stability_goal=0.4):
    """Determine binning by combining neighbouring bins until we reach desired purity & stability"""
    arr2d, _ = cu.th2_to_arr(h2d)
    reco_bin_edges = cu.get_bin_edges(h2d, 'Y')
    gen_bin_edges = cu.get_bin_edges(h2d, 'X')

    new_bin_edges = np.array(reco_bin_edges).reshape(1, len(reco_bin_edges))

    # assumes both axes have same dimension!
    bin_ind = 0
    counter = 0  # safety measure
    while bin_ind < len(arr2d)-1 and counter < 10000:
        counter += 1
        arr2d_renormx = renorm(arr2d, axis=0)
        arr2d_renormy = renorm(arr2d, axis=1)
        purity = arr2d_renormy[bin_ind][bin_ind]
        stability = arr2d_renormx[bin_ind][bin_ind]
        if purity > purity_goal and stability > stability_goal:
            # print("found bin")
            print("bin_ind:", bin_ind, "purity: %.3f" % purity, "stability: %.3f" % stability)
            bin_ind += 1
            continue
        else:
            # print("combining bin", bin_ind, "/", len(arr2d))
            # combine rows & columns (use transpose for latter)
            arr2d = concat_row(arr2d, bin_ind)
            arr2d = concat_row(arr2d.T, bin_ind).T
            new_bin_edges = np.delete(new_bin_edges, bin_ind+1)  # keep track of new binning
            continue

    # keep [1] to be same as next [0], otherwise you lose a bin later when
    # making the rebinned TH2
    these_bins = [list(x) for x in zip(new_bin_edges[:-1], new_bin_edges[1:])]
    # manual hack for last bin
    these_bins[-1][1] = reco_bin_edges[-1]
    print(these_bins)
    return these_bins


def rebin_2d_hist(h2d, new_binning_x, new_binning_y):
    """Rebin a 2D histogram according to specific bin edges for x & y axes

    new_binning_x, new_binning_y are lists of tuple pairs of bin edges
    e.g. [(0, 1), (1, 4), (4, 10)]
    """
    bin_edges_x = [b[0] for b in new_binning_x]
    bin_edges_x.append(new_binning_x[-1][1])

    bin_edges_y = [b[0] for b in new_binning_y]
    bin_edges_y.append(new_binning_y[-1][1])

    print("rebin_2d_hist, new axes:", bin_edges_x, bin_edges_y)

    new_hist = ROOT.TH2D(
        h2d.GetName()+"Rebin",
        ';'.join([h2d.GetTitle(), h2d.GetXaxis().GetTitle(), h2d.GetYaxis().GetTitle()]),
        len(new_binning_x),
        array('d', bin_edges_x),
        len(new_binning_y),
        array('d', bin_edges_y)
    )

    nbins_x_orig = h2d.GetNbinsX()
    bins_x_orig = cu.get_bin_edges(h2d_orig, 'X')

    nbins_y_orig = h2d.GetNbinsY()
    bins_y_orig = cu.get_bin_edges(h2d_orig, 'Y')

    # TODO: handle under/overflow
    for xind, xbin in enumerate(new_binning_x, 1):  # these are tuples
        for yind, ybin in enumerate(new_binning_y, 1):
            # Find all the old bins that correspond to this new bin, get their contents
            new_bin_content = 0
            new_bin_error = 0
            # TODO handle error - reset to sqrt(N)?
            for xind_orig, xbin_orig in enumerate(bins_x_orig[:-1], 1):
                for yind_orig, ybin_orig in enumerate(bins_y_orig[:-1], 1):
                    if (xbin_orig >= xbin[0] and xbin_orig < xbin[1] and
                        ybin_orig >= ybin[0] and ybin_orig < ybin[1]):
                        new_bin_content += h2d.GetBinContent(xind_orig, yind_orig)
                        # print("For new bin", xbin, ybin, "using contents from", xbin_orig, bins_x_orig[xind_orig], ybin_orig, bins_y_orig[yind_orig])
                        # print("orig bin", xind_orig, yind_orig)

            new_hist.SetBinContent(xind, yind, new_bin_content)
            # print("Setting bin", xind, yind)
            new_hist.SetBinError(xind, yind, new_bin_error)

    return new_hist


def make_rebinned_2d_hist(h2d, new_binning, use_half_width_y=False):
    """Rebin 2D histogram using new binning.

    new_binning is list of tuple pairs of bin edges
    e.g. [(0, 1), (1, 4), (4, 10)]

    If use_half_width_y=False, uses new_binning for both x (gen) & y (reco) axes
    If True, creates bins that are half the width of new_binning for
    the y axes (reco)
    """
    if use_half_width_y:
        # create half width bins from new_binning
        reco_binning = []
        reco_bin_edges = cu.get_bin_edges(h2d, 'Y')
        for s, e in new_binning:
            ideal_mid = (s+e)/2.
            mid = reco_bin_edges[bisect.bisect_left(reco_bin_edges, ideal_mid)]
            reco_binning.append((s, mid))
            reco_binning.append((mid, e))
        return rebin_2d_hist(h2d, new_binning, reco_binning)
    else:
        return rebin_2d_hist(h2d, new_binning, new_binning)


def make_plots(h2d, var_dict, plot_dir, append="",
               plot_migrations=True, plot_reco=True, plot_gen=True):
    """Plot a 2D hist, with copies renormalised by row and by column.

    Also optionally plot migrations as 1D plot,
    as well as 1D projections for reco (y) & gen (x)
    """
    # Plot original 2D map, no renormalizing by axis
    canv = ROOT.TCanvas("c"+cu.get_unique_str(), "", 700, 600)
    canv.SetTicks(1, 1)
    if var_dict.get("log", False):
        canv.SetLogx()
        canv.SetLogy()
    pad = ROOT.gPad
    pad.SetBottomMargin(0.12)
    pad.SetLeftMargin(0.13)
    pad.SetRightMargin(0.12)
    xtitle_offset = 1.4
    ytitle_offset = xtitle_offset * 1.1
    h2d.SetTitleOffset(xtitle_offset, 'X')
    h2d.SetTitleOffset(ytitle_offset, 'Y')
    h2d.SetMinimum(1E-3)
    if var_dict.get('log', False):
        h2d.GetXaxis().SetLimits(1, 150)
        h2d.GetYaxis().SetLimits(1, 150)
    if "title" in var_dict:
        h2d.SetTitle(var_dict['title'].replace("\n", ", "))
    h2d.Draw("COLZ")
    output_filename = os.path.join(plot_dir, var_dict['name']+"_%s.%s" % (append, OUTPUT_FMT))
    output_dir = os.path.dirname(output_filename)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    canv.SaveAs(output_filename)

    canv.SetLogz()
    output_filename = os.path.join(plot_dir, var_dict['name']+"_%s_logZ.%s" % (append, OUTPUT_FMT))
    canv.SaveAs(output_filename)

    # Plot 2D map, renormalized by row
    canv.SetLogz(0)
    h2d_renorm_y = cu.make_normalised_TH2(h2d, 'Y', recolour=False, do_errors=False)
    marker_size = 0.8
    h2d_renorm_y.SetMarkerSize(marker_size)
    h2d_renorm_y.SetMaximum(1)
    draw_opt = "COLZ"
    if plot_migrations:
        draw_opt += " TEXT45"
    h2d_renorm_y.SetMinimum(1E-3)
    h2d_renorm_y.Draw(draw_opt)
    xtitle_offset = 1.5
    h2d_renorm_y.SetTitleOffset(xtitle_offset, 'X')

    yax = h2d_renorm_y.GetYaxis()
    h2d_renorm_y.SetTitleOffset(ytitle_offset, 'Y')
    canv.Update()

    canv.SaveAs(os.path.join(plot_dir, "%s_%s_renormY_linZ.%s" % (var_dict['name'], append, OUTPUT_FMT)))
    canv.SetLogz()
    h2d_renorm_y.SetMaximum(1)
    h2d_renorm_y.SetMinimum(1E-3)
    canv.SaveAs(os.path.join(plot_dir, "%s_%s_renormY_logZ.%s" % (var_dict['name'], append, OUTPUT_FMT)))

    # Plot 2D map, renormalized by column
    canv.Clear()
    canv.SetLogz(0)
    h2d_renorm_x = cu.make_normalised_TH2(h2d, 'X', recolour=False, do_errors=False)
    h2d_renorm_x.SetMarkerSize(marker_size)
    h2d_renorm_x.SetMaximum(1)
    h2d_renorm_x.SetMinimum(1E-3)
    h2d_renorm_x.Draw(draw_opt)

    h2d_renorm_x.SetTitleOffset(xtitle_offset, 'X')

    yax = h2d_renorm_x.GetYaxis()
    h2d_renorm_x.SetTitleOffset(ytitle_offset, 'Y')
    yax.SetMoreLogLabels()
    canv.Update()
    canv.SaveAs(os.path.join(plot_dir, "%s_%s_renormX_linZ.%s" % (var_dict['name'], append, OUTPUT_FMT)))
    canv.SetLogz()
    h2d_renorm_x.SetMaximum(1)
    h2d_renorm_x.SetMinimum(1E-3)
    canv.SaveAs(os.path.join(plot_dir, "%s_%s_renormX_logZ.%s" % (var_dict['name'], append, OUTPUT_FMT)))

    # Plot migrations as 1D hists
    if plot_migrations:
        output_filename = os.path.join(plot_dir, "%s_%s_migration_summary.%s" % (var_dict['name'], append, OUTPUT_FMT))
        qgg.make_migration_summary_plot(h2d_renorm_x=h2d_renorm_x,
                                        h2d_renorm_y=h2d_renorm_y,
                                        xlabel=var_dict['var_label'],
                                        output_filename=output_filename,
                                        title=var_dict.get('title', ''))

    if plot_reco or plot_gen:
        conts = []
        if plot_reco:
            h_reco = h2d.ProjectionY(cu.get_unique_str(), 0, -1, "e")
            conts.append(Contribution(h_reco, label="Reco", normalise_hist=True,
                                      line_color=ROOT.kRed, line_width=2))
        if plot_gen:
            h_gen = h2d.ProjectionX(cu.get_unique_str(), 0, -1, "e")
            conts.append(Contribution(h_gen, label="Gen", normalise_hist=True,
                                      line_color=ROOT.kBlue, line_width=2))

        plot = Plot(conts, what='hist', has_data=False,
                    title=var_dict.get('title', ''),
                    xtitle=var_dict['var_label'], ytitle='p.d.f.')
        plot.plot("NOSTACK HISTE")
        bits = []
        if plot_reco:
            bits.append("reco")
        if plot_gen:
            bits.append("gen")
        content = "_".join(bits)
        output_filename = os.path.join(plot_dir, "%s_%s_1d_%s.%s" % (var_dict['name'], append, content, OUTPUT_FMT))
        plot.save(output_filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        help='Input ROOT file to process.')
    parser.add_argument("--outputDir",
                        help="Directory to put output plot dirs into",
                        default=None)
    parser.add_argument("--outputFile",
                        help="Output ROOT file for rebinned 2D hists",
                        default=None)
    # following useful for comparisons of migrations for systematics
    parser.add_argument("--rebinThisInput",
                        help="Apply new binning to these input file(s)",
                        default=None,
                        action="append")
    parser.add_argument("--rebinThisLabel",
                        help="Labels for files to be rebinned. Must match entry in --rebinThisInput.",
                        default=None,
                        action="append")
    parser.add_argument("--rebinThisOutputFile",
                        help="Output ROOT file for rebinned 2D hists. Must match entry in --rebinThisInput.",
                        default=None,
                        action="append")
    parser.add_argument("--target",
                        help="Target purity & stability as a fraction",
                        default=0.4,
                        type=float)
    # acceptable_metrics = ['gausfit', 'quantile']
    # parser.add_argument("--metric",
    #                     help="Metric for deciding bin width.",
    #                     default=acceptable_metrics[0],
    #                     choices=acceptable_metrics)
    args = parser.parse_args()

    input_dir, input_basename = os.path.split(args.input)
    default_plot_dir = os.path.join(input_dir, "rebinning_"+os.path.splitext(input_basename)[0])
    plot_dir = args.outputDir if args.outputDir else default_plot_dir

    default_output_root_filename = os.path.join(plot_dir, "rebinned_" + input_basename)
    output_root_filename = args.outputFile if args.outputFile else default_output_root_filename
    cu.check_dir_exists_create(os.path.dirname(os.path.abspath(output_root_filename)))
    output_tfile = cu.open_root_file(output_root_filename, 'RECREATE')

    source_plot_dir_names = None
    region_labels = None
    if "qcd" in args.input.lower():
        source_plot_dir_names = ["Dijet_QG_central_tighter", "Dijet_QG_forward_tighter", "Dijet_QG_central_tighter_groomed", "Dijet_QG_forward_tighter_groomed"]
        region_labels = [qgc.QCD_Dijet_CEN_LABEL, qgc.QCD_Dijet_FWD_LABEL, qgc.QCD_Dijet_CEN_GROOMED_LABEL, qgc.QCD_Dijet_FWD_GROOMED_LABEL]
    elif "dyjetstoll" in args.input.lower():
        source_plot_dir_names = ["ZPlusJets_QG"]
        region_labels = [qgc.DY_ZpJ_LABEL]
    else:
        raise RuntimeError("No idea which region we're using")

    rebin_results_dict = OrderedDict()

    pt_regions = [
        {
            "append": "",
            "title": "p_{T}^{Reco} > 30 GeV",
        },
        {
            "append": "_lowPt",
            "title": "30 < p_{T}^{Reco} < 100 GeV",
        },
        {
            "append": "_midPt",
            "title": "100 < p_{T}^{Reco} < 250 GeV",
        },
        {
            "append": "_highPt",
            "title": "p_{T}^{Reco} > 250 GeV",
        },
        # {
        #     "append": "_charged",
        #     "title": "p_{T}^{Reco} > 30 GeV (only charged constits.)",
        # },
        # {
        #     "append": "_charged_lowPt",
        #     "title": "30 < p_{T}^{Reco} < 100 GeV (only charged constits.)",
        # },
        # {
        #     "append": "_charged_midPt",
        #     "title": "100 < p_{T}^{Reco} < 250 GeV (only charged constits.)",
        # },
        # {
        #     "append": "_charged_highPt",
        #     "title": "p_{T}^{Reco} > 250 GeV (only charged constits.)",
        # },
    ]

    input_tfile = cu.open_root_file(args.input)

    for source_plot_dir_name, region_label in zip(source_plot_dir_names, region_labels):

        for angle in qgc.COMMON_VARS[:]:

            response_maps_dict = {} # store 2D maps for later

            for pt_region_dict in pt_regions[:]:
                var_dict = {
                    "name": "%s/%s%s" % (source_plot_dir_name, angle.var, pt_region_dict['append']),
                    "var_label": "%s (%s)" % (angle.name, angle.lambda_str),
                    "title": "%s\n%s" % (region_label, pt_region_dict['title']),
                }

                do_rel_response = False  # Use relative response instead - only makes sense with quantiles/gaussian fit?
                full_var_name = var_dict['name']
                if do_rel_response:
                    full_var_name += "_rel_response"
                else:
                    full_var_name += "_response"

                print(full_var_name)

                h2d_orig = cu.get_from_tfile(input_tfile, full_var_name)

                # Make plots with original fine equidistant binning
                # -------------------------------------------------
                make_plots(h2d_orig, var_dict, plot_dir=plot_dir, append="orig", plot_migrations=False)

                response_maps_dict[var_dict['name']] = h2d_orig

                # metric = "gausfit"
                # metric = "quantile"
                # new_binning = calc_variable_binning(h2d_orig, plot_dir, args.metric)

                # Calculate new binning
                # ---------------------
                # goal = 0.5
                goal = args.target
                new_binning = calc_variable_binning_purity_stability(h2d_orig, purity_goal=goal, stability_goal=goal)
                rebin_results_dict[var_dict['name']] = new_binning

                # Rebin 2D heatmap
                # ----------------
                h2d_rebin = make_rebinned_2d_hist(h2d_orig, new_binning, use_half_width_y=False)

                # Plot with new binning
                # ---------------------
                make_plots(h2d_rebin, var_dict, plot_dir=plot_dir, append="rebinned",
                           plot_migrations=True, plot_reco=True, plot_gen=True)

                # Cache renormed plots here for migration plots
                h2d_renorm_x = cu.make_normalised_TH2(h2d_rebin, 'X', recolour=False, do_errors=False)
                output_tfile.WriteTObject(h2d_renorm_x)  # we want renormalised by col for doing folding
                response_maps_dict[var_dict['name']+"_renormX"] = h2d_renorm_x

                h2d_renorm_y = cu.make_normalised_TH2(h2d_rebin, 'Y', recolour=False, do_errors=False)
                output_tfile.WriteTObject(h2d_renorm_y)  # save anyway for completeness
                response_maps_dict[var_dict['name']+"_renormY"] = h2d_renorm_y

                # Now rebin any other input files with the same hist using the new binning
                # ------------------------------------------------------------------------
                if args.rebinThisInput and len(args.rebinThisInput) > 0:
                    contributions = qgg.migration_plot_components(h2d_renorm_x, h2d_renorm_y, var_dict['var_label'])

                    lens = [len(x) for x in [args.rebinThisInput, args.rebinThisLabel, args.rebinThisOutputFile]]
                    if any(l != lens[0] for l in lens[1:]):
                        raise RuntimeError("Number of --rebinThisInput, --rebinThisLabel, --rebinThisOutputFile should be the same - need one of each per extra input file")

                    for ind, (other_input, other_label, other_output_filename) in enumerate(zip(args.rebinThisInput, args.rebinThisLabel, args.rebinThisOutputFile)):
                        print(other_label)
                        tfile_other = cu.open_root_file(other_input)
                        tfile_other_out = cu.open_root_file(other_output_filename, 'RECREATE')

                        h2d_other = cu.get_from_tfile(tfile_other, full_var_name)
                        h2d_rebin_other = make_rebinned_2d_hist(h2d_other, new_binning, use_half_width_y=False)


                        make_plots(h2d_rebin_other, var_dict, plot_dir=plot_dir+"_"+other_label.replace(" ", "_"), append="rebinned", plot_migrations=True)

                        h2d_renorm_x_other = cu.make_normalised_TH2(h2d_rebin_other, 'X', recolour=False, do_errors=False)
                        tfile_other_out.WriteTObject(h2d_renorm_x_other)
                        h2d_renorm_y_other = cu.make_normalised_TH2(h2d_rebin_other, 'Y', recolour=False, do_errors=False)
                        contributions_other = qgg.migration_plot_components(h2d_renorm_x_other, h2d_renorm_y_other, var_dict['var_label'])
                        for ind, c in enumerate(contributions_other):
                            c.obj.SetLineStyle(ind+2)
                            c.label += " [%s]" % other_label
                            c.subplot = contributions[ind].obj
                        contributions.extend(contributions_other)
                        tfile_other_out.Close()

                    # Plot all these rebinned extra inputs alongside the main map
                    # -----------------------------------------------------------
                    for c in contributions:
                        c.obj.SetLineWidth(2)

                    binning = [h2d_renorm_x.GetXaxis().GetBinLowEdge(bin_ind) for bin_ind in range(1, h2d_renorm_x.GetNbinsX()+2)]
                    xlim = [binning[0], binning[-1]]

                    # plot = Plot(contributions, what='hist', xlim=xlim, ylim=[1e-3, 2], xtitle=var_dict['var_label'], has_data=False)
                    plot = Plot(contributions, what='hist', xlim=xlim, ylim=[0, 1.25],
                                xtitle=var_dict['var_label'], has_data=False,
                                subplot_type='ratio', subplot_title='Syst / nominal',
                                subplot_limits=[0.9, 1.1])
                    y1 = 0.15
                    y1 = 0.65
                    plot.legend.SetNColumns(len(args.rebinThisInput)+1)
                    plot.legend.SetX1(0.15)
                    plot.legend.SetX2(0.95)
                    plot.legend.SetY1(y1)
                    plot.legend.SetY2(y1+0.25)
                    plot.legend.SetTextSize(0.015)
                    plot.plot("NOSTACK HISTE")
                    plot.legend.SetFillStyle(1001)
                    plot.legend.SetFillColorAlpha(ROOT.kWhite, 0.75)
                    if var_dict.get('log', None):
                        plot.set_logx()
                    # plot.set_logy()
                    plot.main_pad.cd()
                    lines = []
                    # Boundary lines
                    # for val in [1, 0.5, 1e-1, 1e-2, 1e-3]:
                    #     line = ROOT.TLine(xlim[0], val, xlim[1], val)
                    #     line.SetLineStyle(2)
                    #     line.SetLineColor(ROOT.kGray+2)
                    #     lines.append(line)
                    #     line.Draw("same")
                    output_filename = os.path.join(plot_dir, "%s_combined_migration_summary.%s" % (var_dict['name'], OUTPUT_FMT))
                    plot.save(output_filename)

                    # Do another version where you plot the relative difference wrt nominal for each bin
                        # Contribution(hist_purity, label="Purity (gen right)", line_color=col_purity, marker_color=col_purity),
                        # Contribution(hist_stability, label="Stability (reco right)", line_color=col_stability, marker_color=col_stability),
                        # Contribution(hist_xfer_down, label="-1 reco bin", line_color=col_xfer_down, marker_color=col_xfer_down),
                        # Contribution(hist_xfer_down2, label="-2 reco bin", line_color=col_xfer_down2, marker_color=col_xfer_down2),
                        # # Contribution(hist_xfer_down3, label="3 lower reco bin", line_color=col_xfer_down3, marker_color=col_xfer_down3),
                        # Contribution(hist_xfer_up, label="+1 reco bin", line_color=col_xfer_up, marker_color=col_xfer_up),
                        # Contribution(hist_xfer_up2,
                    # new_contributions = []
                    # ref_cont = contributions[0]
                    # for cont in contributions[1:]:
                    #     new_cont = deepcopy(cont)
                    #     for ref, var in zip(ref_cont, cont):
                    #         new_obj = var.obj.Clone()
                    #         new_obj.Divide(ref)

            # Apply binning scheme dervied from one pT region to others
            # ---------------------------------------------------------
            rebin_other_pt_regions = True
            reference_pt_region = "_midPt"  # correspond to 'append' key in a dict
            # reference_pt_region = "_lowPt"  # correspond to 'append' key in a dict
            if rebin_other_pt_regions:

                this_pt_dict = [x for x in pt_regions if x['append'] == reference_pt_region][0]
                print("Rebinning other pt regions using:")
                print(this_pt_dict)

                # God this is awful, such disconnect between rebinning, names, hists

                # First find rebinning scheme
                reference_binning = None
                for hname, h2d in response_maps_dict.items():
                    if "_renorm" in hname:
                        continue
                    if reference_pt_region in hname:
                        reference_binning = rebin_results_dict[hname]
                        break

                if not reference_binning:
                    raise RuntimeError("Something went wrong looking up reference_binning")

                # Now rebin all other 2D hists
                for hname, h2d in response_maps_dict.items():
                    if reference_pt_region in hname or "_renorm" in hname:
                        continue
                    print("Rebinning", hname)
                    h2d_rebin = make_rebinned_2d_hist(h2d, reference_binning, use_half_width_y=False)
                    var_dict = {
                        "name": hname,
                        "var_label": "%s (%s)" % (angle.name, angle.lambda_str),
                        "title": h2d.GetTitle() + "\n(rebinned for %s)" % this_pt_dict['title'],
                    }
                    make_plots(h2d_rebin, var_dict, plot_dir=plot_dir,
                               append="rebinned_for%s" % (reference_pt_region),
                               plot_migrations=True)

                # Now rebin any other input files with the same hist using the new binning
                # ------------------------------------------------------------------------
                if args.rebinThisInput and len(args.rebinThisInput) > 0:
                    # contributions = qgg.migration_plot_components(h2d_renorm_x, h2d_renorm_y, var_dict['var_label'])

                    lens = [len(x) for x in [args.rebinThisInput, args.rebinThisLabel, args.rebinThisOutputFile]]
                    if any(l != lens[0] for l in lens[1:]):
                        raise RuntimeError("Number of --rebinThisInput, --rebinThisLabel, --rebinThisOutputFile should be the same - need one of each per extra input file")


                    for pt_region_dict in pt_regions[:]:

                        var_dict = {
                            "name": "%s/%s%s" % (source_plot_dir_name, angle.var, pt_region_dict['append']),
                            "var_label": "%s (%s)" % (angle.name, angle.lambda_str),
                            "title": "%s\n%s\n(rebinned for %s)" % (region_label, pt_region_dict['title'], this_pt_dict['title']),
                        }
                        contributions = []

                        full_var_name = var_dict['name']
                        if do_rel_response:
                            full_var_name += "_rel_response"
                        else:
                            full_var_name += "_response"

                        print(var_dict)

                        for ind, (other_input, other_label, other_output_filename) in enumerate(zip(args.rebinThisInput, args.rebinThisLabel, args.rebinThisOutputFile)):
                            print(other_label)
                            tfile_other = cu.open_root_file(other_input)
                            tfile_other_out = cu.open_root_file(other_output_filename, 'RECREATE')

                            h2d_other = cu.get_from_tfile(tfile_other, full_var_name)
                            h2d_rebin_other = make_rebinned_2d_hist(h2d_other, reference_binning, use_half_width_y=False)
                            print(h2d_other)
                            print(h2d_rebin_other)

                            this_var_dict = deepcopy(var_dict)
                            this_var_dict['title'] += "\n[%s]" % other_label
                            make_plots(h2d_rebin_other, this_var_dict, plot_dir=plot_dir+"_"+other_label.replace(" ", "_"), append="rebinned_for%s" % (reference_pt_region), plot_migrations=True)

                            h2d_renorm_x_other = cu.make_normalised_TH2(h2d_rebin_other, 'X', recolour=False, do_errors=False)
                            tfile_other_out.WriteTObject(h2d_renorm_x_other)
                            h2d_renorm_y_other = cu.make_normalised_TH2(h2d_rebin_other, 'Y', recolour=False, do_errors=False)
                            contributions_other = qgg.migration_plot_components(h2d_renorm_x_other, h2d_renorm_y_other, var_dict['var_label'])
                            for ind, c in enumerate(contributions_other):
                                # c.obj.SetLineStyle(ind+2)
                                c.label += " [%s]" % other_label
                                # c.subplot = contributions[ind].obj
                            contributions.extend(contributions_other)
                            tfile_other_out.Close()

                        # Plot all these rebinned extra inputs alongside the main map
                        # -----------------------------------------------------------
                        for c in contributions:
                            c.obj.SetLineWidth(2)

                        binning = [contributions[0].obj.GetXaxis().GetBinLowEdge(bin_ind) for bin_ind in range(1, contributions[0].obj.GetNbinsX()+2)]
                        xlim = [binning[0], binning[-1]]
                        # xlim = None

                        # plot = Plot(contributions, what='hist', xlim=xlim, ylim=[1e-3, 2], xtitle=var_dict['var_label'], has_data=False)
                        # plot = Plot(contributions, what='hist', xlim=xlim, ylim=[0, 1.25],
                        plot = Plot(contributions, what='hist', xlim=xlim, ylim=[5e-3, 5],
                                    xtitle=var_dict['var_label'], has_data=False,
                                    title=var_dict['title'],
                                    # subplot_type='ratio', subplot_title='Syst / nominal',
                                    # subplot_limits=[0.9, 1.1])
                                    )
                        y1 = 0.15
                        y1 = 0.65
                        plot.legend.SetNColumns(len(args.rebinThisInput)+1)
                        y1 = 0.15
                        plot.legend.SetX1(0.5)
                        plot.legend.SetY1(y1)
                        plot.legend.SetY2(y1+0.25)
                        plot.legend.SetTextSize(0.015)
                        plot.plot("NOSTACK HISTE")
                        plot.legend.SetFillStyle(1001)
                        plot.legend.SetFillColorAlpha(ROOT.kWhite, 0.75)
                        if var_dict.get('log', None):
                            plot.set_logx()
                        plot.set_logy()
                        plot.main_pad.cd()
                        lines = []
                        # Boundary lines
                        for val in [1, 0.5, 1e-1, 1e-2]:
                            line = ROOT.TLine(xlim[0], val, xlim[1], val)
                            line.SetLineStyle(2)
                            line.SetLineColor(ROOT.kGray+2)
                            lines.append(line)
                            line.Draw("same")
                        output_filename = os.path.join(plot_dir+"_"+other_label.replace(" ", "_"), "%s_combined_migration_summary_rebinned_for%s.%s" % (var_dict['name'], reference_pt_region, OUTPUT_FMT))
                        plot.save(output_filename)


    output_tfile.Close()
    input_tfile.Close()

    # Save new binning to txt file
    # ----------------------------
    output_txt = os.path.splitext(args.input)[0] + ".txt"
    parts = os.path.split(output_txt)
    output_txt = os.path.join(args.outputDir, 'binning_'+parts[1])
    with open(output_txt, 'w') as fout:
        for k, v in rebin_results_dict.items():
            # turn pairs of bin edges into one list
            bins = [vv[0] for vv in v]
            bins.append(v[-1][1])
            fout.write("%s: %s\n" % (k, bins))

    print("saved new binning to", output_txt)
    print("saved rebinned 2D hists to", output_root_filename)

