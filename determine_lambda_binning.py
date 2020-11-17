#!/usr/bin/env python

"""Figure out optimial binning for lambda variables.

Currently uses minimum required purity & stability in response matrix.

Can also rebin other inputs according to the binning scheme from the main input.
"""

from __future__ import print_function, division

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
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# My stuff
import qg_general_plots as qgp
import qg_common as qgc
from comparator import Contribution, Plot
import common_utils as cu
import metric_calculators as metrics

ROOT.gStyle.SetPaintTextFormat(".3f")

# Control output format
OUTPUT_FMT = "pdf"


def calc_variable_binning_metric(h2d, plot_dir, metric):
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
    # If you get "TypeError: No loop matching the specified signature and casting was found for ufunc true_divide"
    # then you should create your arr2d with 'type=float'.
    # I don't understand why: https://github.com/numpy/numpy/issues/10565
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


def get_1D_bins_from_arr2d(arr2d, axis='gen'):
    return arr2d.sum(axis=0 if axis.lower() == 'gen' else 1)


def get_unsmooth_bins(arr2d, axis='gen'):
    """Find bins that are not smooth, ie they are spikey"""
    # Get 1D hist, by summing over bins of other axis in 2D hist
    bins = get_1D_bins_from_arr2d(arr2d, axis)
    diffs = np.diff(bins)
    diff_signs = np.sign(np.diff(bins))
    bad_diff_bins = []
    for i in range(len(diff_signs)-2):
        if i == 0:
            # edge case: consider i, i+1
            # but avoid cases where binning is coarse - could be genuine
            # 1st bin and not a spike i.e. +ve gradient up to peak
            if diff_signs[i] != diff_signs[i+1] and diffs[i] < 0:
                print('>>>>! Bin', i, 'diff wrong sign:', diffs[i], diffs[i+1])
                bad_diff_bins.append(i)
        else:
            # consider i-1, i, i+1
            if (diff_signs[i] != diff_signs[i-1] and
                diff_signs[i] != diff_signs[i+1]):
                print('>>>>! Bin', i, 'diff wrong sign:', diffs[i-1], diffs[i], diffs[i+1])
                bad_diff_bins.append(i)
    return bad_diff_bins


def printable_bins(bins):
    """Stop printing out floating-point errors

    Trick by converting to str, then back to float

    Parameters
    ----------
    bins : TYPE
        Description

    Returns
    -------
    list[[bin edge, bin edge]]
    """
    s = []
    fmt = '%.3g'
    for (a, b) in bins:
        s.append([float(fmt % a), float(fmt % b)])
    return s


def calc_variable_binning_purity_stability(h2d,
                                           purity_goal=0.4,
                                           stability_goal=0.4,
                                           integer_binning=False,
                                           spike_smoothing=True):
    """Determine binning by combining neighbouring bins until we reach desired purity & stability

    If integer_binning, then ensure only even gaps between bins, such that when
    they are halved, you get an equal proportion on each side. This ensures
    the binned variable doesn't have spikes.
    It also subtracts 0.5 from each bin edge to ensure integers fall in middle of bin.

    e.g. 10 - 14: becomes [10, 12, 14] when halved, so covers [[10, 12), [12, 14)],
    which mean for integers: [[10, 11], [12, 13]], and bin widths [2, 2].
    Thus for a uniform true distribution, get uniform binned distribution
    (equal number of integers fall into each bin, each bin has same bin width)

    e.g. 10 - 13: becomes [10, 11.5, 13] when halved, so covers [[10, 11.5), [11.5, 13)]
    which means for integers: [[10, 11], [12]], and bin width [1.5, 1.5].
    So here the first bin would get double the number of events (for a uniform distribution),
    despite both bin widths being the same.

    Note that it doesn't matter if the bin edges are integers or half-integers;
    but they must be consistent.
    """
    arr2d_orig, _ = cu.th2_to_arr(h2d)
    reco_bin_edges = cu.get_bin_edges(h2d, 'Y')
    # gen_bin_edges = cu.get_bin_edges(h2d, 'X')
    # print(reco_bin_edges)
    new_bin_edges = np.array(reco_bin_edges[:])

    # smooth_bins = False
    # smooth_bin_counter = 0 # safety
    # bad_diff_bins, bad_diff_bin_edges = [], []
    # while not smooth_bins and smooth_bin_counter < 100:
    #     smooth_bin_counter += 1

    # assumes both axes have same dimension!
    arr2d = np.copy(arr2d_orig)
    bin_ind = 0
    counter = 0  # safety measure
    while bin_ind < len(arr2d)-1 and counter < 10000:
        counter += 1
        arr2d_renormx = renorm(arr2d, axis=0)  # renormed per x (gen) bin
        arr2d_renormy = renorm(arr2d, axis=1)  # renormed per y (reco) bin
        purity = arr2d_renormy[bin_ind][bin_ind]  # fraction in a reco bin that are actually from that gen bin
        stability = arr2d_renormx[bin_ind][bin_ind] # fraction in a gen bin that make it to that reco bin
        has_even_interval = (new_bin_edges[bin_ind+1]-new_bin_edges[bin_ind]) % 2 == 0
        # print(purity, stability)

        if (purity > purity_goal and
            stability > stability_goal and
            ((integer_binning and has_even_interval) or not integer_binning)):
            # FOUND OUR DESIRED BIN
            # print("found bin")
            # print(new_bin_edges)
            print("bin_ind:", bin_ind, " = ", new_bin_edges[bin_ind], "-", new_bin_edges[bin_ind+1], "purity: %.3f" % purity, "stability: %.3f" % stability)
            bin_ind += 1
            continue
        else:
            # print("combining bin", bin_ind, "/", len(arr2d))
            # combine rows & columns (use transpose for latter)
            arr2d = concat_row(arr2d, bin_ind)
            arr2d = concat_row(arr2d.T, bin_ind).T
            # keep track of new binning by removing the bin edge we just tried
            # and found failed the purity/stability requirements
            new_bin_edges = np.delete(new_bin_edges, bin_ind+1)
            # print(new_bin_edges)
            continue

    # handle last bin: if it fails target, just merge last 2 bins together
    if purity < purity_goal or stability < stability_goal:
        # print("Purity &/ stability not meeting target in last bin - merging ")
        # remove the penultimate entry, which would have been the lower edge of the last bin
        # with the insufficient purity/stability
        new_bin_edges = np.delete(new_bin_edges, -2)
    print(new_bin_edges)

    if spike_smoothing:
        # Check for smooth/continuous distribution: if not, then need to reconsider
        # TODO check either of reco/gen, or both?
        axis = 'reco'
        bad_diff_bins = get_unsmooth_bins(arr2d, axis=axis)
        # handle 0 especially, since we can just combine the first few bins

        bad_counter = 0
        while 0 in bad_diff_bins and bad_counter < 1:
            bad_counter += 1

            new_arr2d = np.copy(arr2d)
            zero_counter = 0
            hist_bins = get_1D_bins_from_arr2d(arr2d, axis)
            # find end bin - which has a height similar to 1st bin, add a bit of padding
            start_bin = 0
            # +1 since offset by 1 (for factor < 0, np.where will always pick 0th entry, so we check [1:])
            # +1 to well-cover the range
            end_bin = int(np.where(hist_bins[1:] > 1.*hist_bins[0])[0][0]) + 1
            print('end_bin:', end_bin)

            # Now do clustering, to try and get similar heights
            bins_to_cluster = hist_bins[start_bin:end_bin]
            print("len(bins_to_cluster)", len(bins_to_cluster))
            if len(bins_to_cluster) == 2:
                # just merge the first two
                row_ind = 0
                print("Joining", row_ind)
                new_arr2d = concat_row(new_arr2d, row_ind)
                new_arr2d = concat_row(new_arr2d.T, row_ind).T
                new_bin_edges = np.delete(new_bin_edges, row_ind+1)

            else:
                print('new_bin_edges[start_bin:end_bin]', new_bin_edges[start_bin:end_bin])

                # the smallest number of bins is 1
                # the largest number is constrained by the 0th bin (or largest?)
                # that we are trying to merge, since we don't want it much larger
                # than the average
                sum_bins = bins_to_cluster.sum()
                num_target_bins = np.arange(1, len(bins_to_cluster), 1)
                # num_target_bins = np.array([2])
                averages = sum_bins / num_target_bins
                mask = bins_to_cluster.max() < 1.2 * averages
                averages = averages[mask]
                num_target_bins = num_target_bins[mask]

                results = []

                print('bins_to_cluster', bins_to_cluster)

                # iterate over different number of target bins
                for n_target, ave in zip(num_target_bins, averages):

                    this_arr2d = np.copy(new_arr2d)
                    this_new_bin_edges = np.copy(new_bin_edges)

                    print('n_target', n_target)
                    print('target', ave)

                    # figure out most even binning for given target number of bins
                    cumul_bins = bins_to_cluster.cumsum()
                    ideal_cumul = np.linspace(0, sum_bins, n_target+1)
                    indices = np.searchsorted(cumul_bins, ideal_cumul, side='right')  # right to ensure last bin done
                    print('cumul_bins', cumul_bins)
                    print('ideal_cumul', ideal_cumul)

                    # sometimes get 0, 0 as first indices, when 1st bin is too big for
                    # the final 1st bin
                    # this makes the 2nd bin huge, so override that
                    if indices[0] == 0 and indices[1] == 0:
                        indices[1] = 1

                    print('indices', indices)

                    # do rebinning of these bins
                    new_data = []
                    for g, gn in zip(indices[:-1], indices[1:]):
                        new_data.append(bins_to_cluster[g:gn].sum())
                    print("rebinned bins", new_data)

                    # now do rebinning of th2, bin edges
                    row_ind = start_bin-1
                    for g, gn in zip(indices[:-1], indices[1:]):
                        print("outer Joining", g)
                        row_ind += 1
                        rebin_these = range(g, gn-1) # -1 because we merge 2 bins into 1, so only need to concat one fewer
                        # if len(rebin_these) <= 1:
                            # print('...skipping')
                            # continue
                        for gg in rebin_these:
                            print("Joining", row_ind)
                            this_arr2d = concat_row(this_arr2d, row_ind)
                            this_arr2d = concat_row(this_arr2d.T, row_ind).T
                            this_new_bin_edges = np.delete(this_new_bin_edges, row_ind+1)

                    this_hist = get_1D_bins_from_arr2d(this_arr2d, axis)
                    # Calculate std deviation across bins of interest
                    this_std = np.std(np.diff(this_hist[start_bin:end_bin+1]))
                    # Calculate max absolute diff across bins of interest
                    this_max_diff = np.abs(np.diff(this_hist[start_bin:end_bin+1])).max()
                    results.append({
                        'arr2d': this_arr2d,
                        'bin_edges': this_new_bin_edges,
                        'std_dev': this_std,
                        'max_diff': this_max_diff
                    })

                # Choose best binning, based on some metric (min diff, min std dev)
                min_metric = 1E100
                min_metric_n = -1
                for n_bins, r in zip(num_target_bins, results):
                    metric = r['max_diff']
                    if metric < min_metric:  # favour smaller # bins
                        min_metric = metric
                        min_metric_n = n_bins
                    print(n_bins, r['bin_edges'], r['std_dev'], r['max_diff'])
                # min_metric_n = 3
                ind = np.where(num_target_bins == min_metric_n)[0][0]
                new_arr2d = results[ind]['arr2d']
                new_bin_edges = results[ind]['bin_edges']
                print("chosen nbins", min_metric_n)

            bad_diff_bins = get_unsmooth_bins(new_arr2d, axis=axis)

            if 0 in bad_diff_bins:
                print("!!!! I couldn't get rid of bad 1st bin")
            else:
                print("My new bins:", new_bin_edges)


    if integer_binning:
        # ensure integers fall in centre of bin
        # this assumes the original bin edges were integers!
        # if not, comment out this line
        new_bin_edges-= 0.5

    # convert to bin edge pairs, keeping pair[1] to be same as next_pair[0],
    # otherwise you lose a bin later when making the rebinned TH2
    these_bins = [list(x) for x in zip(new_bin_edges[:-1], new_bin_edges[1:])]
    # manual hack for last bin
    last_bin = reco_bin_edges[-1]
    if integer_binning:
        last_bin -= 0.5
    these_bins[-1][1] = last_bin
    print("Final binning:", printable_bins(these_bins))
    return these_bins


def rebin_2d_hist(h2d, new_binning_x, new_binning_y):
    """Rebin a 2D histogram according to specific bin edges for x & y axes

    new_binning_x, new_binning_y are lists of tuple pairs of bin edges
    e.g. [(0, 1), (1, 4), (4, 10)]
    """
    print("rebinning...")
    # convert pairs of bins to list of edges, including upper edge of last bin
    bin_edges_x = [b[0] for b in new_binning_x]
    bin_edges_x.append(new_binning_x[-1][1])

    bin_edges_y = [b[0] for b in new_binning_y]
    bin_edges_y.append(new_binning_y[-1][1])

    # print("rebin_2d_hist, new axes:", bin_edges_x, bin_edges_y)

    new_h2d = ROOT.TH2D(
        h2d.GetName()+"Rebin",
        ';'.join([h2d.GetTitle(), h2d.GetXaxis().GetTitle(), h2d.GetYaxis().GetTitle()]),
        len(new_binning_x),
        array('d', bin_edges_x),
        len(new_binning_y),
        array('d', bin_edges_y)
    )

    # Get original bin edges
    bins_x_orig = cu.get_bin_edges(h2d, 'X')
    bins_y_orig = cu.get_bin_edges(h2d, 'Y')

    # Get original bin contents
    arr, err = cu.th2_to_arr(h2d, do_errors=True)

    # Get map of old bin edges -> new bin edges
    # -1 to start at 0
    bin_groups_x = np.digitize(bins_x_orig, bin_edges_x) - 1
    bin_groups_y = np.digitize(bins_y_orig, bin_edges_y) - 1

    # Count cumulative number of entries in each group (remove last bin as upper edge)
    group_counts_x = np.bincount(bin_groups_x)[:-1].cumsum()
    group_counts_y = np.bincount(bin_groups_y)[:-1].cumsum()

    group_counts_x = np.insert(group_counts_x, 0, 0)
    group_counts_y = np.insert(group_counts_y, 0, 0)

    # Iterate over each group, sum, set new TH2 contents
    for xind, (xl, xh) in enumerate(zip(group_counts_x[:-1], group_counts_x[1:]), 1):
        for yind, (yl, yh) in enumerate(zip(group_counts_y[:-1], group_counts_y[1:]), 1):
            new_bin_content = arr[yl:yh,xl:xh].sum()
            new_bin_err = np.sqrt(np.power(err[yl:yh,xl:xh], 2).sum())
            new_h2d.SetBinContent(xind, yind, new_bin_content)
            new_h2d.SetBinError(xind, yind, new_bin_err)

    print("...done rebinning")
    return new_h2d


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


def make_plots(h2d,
               var_dict,
               plot_dir,
               append="",
               plot_maps=True,
               plot_migrations=True,
               plot_reco=True,
               plot_gen=True,
               true_mean=None,
               true_rms=None,
               div_bin_width=True):
    """Plot a 2D hist, with copies renormalised by row and by column.

    Also optionally plot migrations as 1D plot,
    as well as 1D projections for reco (y) & gen (x)
    """
    # Plot original 2D map, no renormalizing by axis
    h2d_renorm_x, h2d_renorm_y = None, None
    if plot_maps or plot_migrations:
        h2d_renorm_y = cu.make_normalised_TH2(h2d, 'Y', recolour=False, do_errors=True)
        h2d_renorm_x = cu.make_normalised_TH2(h2d, 'X', recolour=False, do_errors=True)

    if plot_maps:
        canv = ROOT.TCanvas("c"+cu.get_unique_str(), "", 700, 600)
        canv.SetTicks(1, 1)
        if var_dict.get("log", False):
            canv.SetLogx()
            canv.SetLogy()
        pad = ROOT.gPad
        pad.SetBottomMargin(0.12)
        pad.SetLeftMargin(0.13)
        pad.SetRightMargin(0.16)
        xtitle_offset = 1.4
        ytitle_offset = xtitle_offset * 1.1
        ztitle_offset = xtitle_offset * 0.9
        h2d.SetTitleOffset(xtitle_offset, 'X')
        h2d.SetTitleOffset(ytitle_offset, 'Y')
        h2d.SetTitleOffset(ztitle_offset, 'Z')
        h2d.SetMinimum(1E-3)
        if var_dict.get('log', False):
            h2d.GetXaxis().SetLimits(1, 150)
            h2d.GetYaxis().SetLimits(1, 150)
        title = "%s;%s (GEN);%s (RECO)" % (var_dict.get('title', '').replace("\n", ", "), var_dict['var_label'], var_dict['var_label'])
        h2d.SetTitle(title)
        h2d.Draw("COLZ")
        old_font_size = ROOT.gStyle.GetTitleFontSize()
        if "#splitline" in var_dict.get('title', ''):
            # need to shrink it down a bit
            pad.SetTopMargin(0.16)
            ROOT.gStyle.SetTitleFontSize(0.035)
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
        marker_size = 0.8
        h2d_renorm_y.SetMarkerSize(marker_size)
        h2d_renorm_y.SetMaximum(1)
        h2d_renorm_y.SetZTitle("P(GEN bin | RECO bin)")
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
        h2d_renorm_x.SetMarkerSize(marker_size)
        h2d_renorm_x.SetMaximum(1)
        h2d_renorm_x.SetMinimum(1E-3)
        h2d_renorm_x.SetZTitle("P(RECO bin | GEN bin)")
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

        ROOT.gStyle.SetTitleFontSize(old_font_size)

    # Plot migrations as 1D hists
    if plot_migrations:
        output_filename = os.path.join(plot_dir, "%s_%s_migration_summary.%s" % (var_dict['name'], append, OUTPUT_FMT))
        qgp.make_migration_summary_plot(h2d_renorm_x=h2d_renorm_x,
                                        h2d_renorm_y=h2d_renorm_y,
                                        xlabel=var_dict['var_label'],
                                        output_filename=output_filename,
                                        title=var_dict.get('title', ''))

    # Do 1D hist(s)
    if plot_reco or plot_gen:
        conts = []
        if plot_reco:
            h_reco = h2d.ProjectionY(cu.get_unique_str(), 0, -1, "e")
            h_reco.Scale(1./h_reco.Integral())
            # Get rebinned mean/RMS to add to plot
            contents, errors = cu.th1_to_ndarray(h_reco)
            centers = cu.get_th1_bin_centers(h_reco)
            areas, centers = metrics.hist_values_to_uarray(bin_areas=contents, bin_centers=centers, bin_errors=errors)
            # print(areas, areas.sum())
            mean_u = metrics.calc_mean_ucert(areas, centers)
            mean, mean_err = mean_u.nominal_value, mean_u.std_dev
            rms_u = metrics.calc_rms_ucert(areas, centers)
            rms, rms_err = rms_u.nominal_value, rms_u.std_dev

            h_reco_div_bin_width = h_reco
            if div_bin_width:
                h_reco_div_bin_width = qgp.hist_divide_bin_width(h_reco)
            conts.append(Contribution(h_reco_div_bin_width,
                                      label="Reco\n(mean = %.3f+-%.3f)\n(RMS = %.3f+-%.3f)" % (mean, mean_err, rms, rms_err),
                                      normalise_hist=False,
                                      marker_color=ROOT.kRed,
                                      line_color=ROOT.kRed, line_width=2))
        if plot_gen:
            h_gen = h2d.ProjectionX(cu.get_unique_str(), 0, -1, "e")
            h_gen.Scale(1./h_gen.Integral())
            # Get rebinned mean/RMS to add to plot
            contents, errors = cu.th1_to_ndarray(h_gen)
            centers = cu.get_th1_bin_centers(h_gen)
            areas, centers = metrics.hist_values_to_uarray(bin_areas=contents, bin_centers=centers, bin_errors=errors)
            mean_u = metrics.calc_mean_ucert(areas, centers)
            mean, mean_err = mean_u.nominal_value, mean_u.std_dev
            rms_u = metrics.calc_rms_ucert(areas, centers)
            rms, rms_err = rms_u.nominal_value, rms_u.std_dev

            h_gen_div_bin_width = h_gen
            if div_bin_width:
                h_gen_div_bin_width = qgp.hist_divide_bin_width(h_gen)
            conts.append(Contribution(h_gen_div_bin_width,
                                      label="Gen\n(mean = %.3f+-%.3f)\n(RMS = %.3f+-%.3f)" % (mean, mean_err, rms, rms_err),
                                      normalise_hist=False,
                                      marker_color=ROOT.kBlue,
                                      line_color=ROOT.kBlue, line_width=2, line_style=2 if plot_reco else 1))

        plot = Plot(conts,
                    what='hist',
                    has_data=False,
                    title=var_dict.get('title', ''),
                    xtitle=var_dict['var_label'],
                    ytitle='p.d.f.' if div_bin_width else "#Delta N / N")
        plot.default_canvas_size = (700, 600)
        plot.legend.SetX1(0.65)
        plot.legend.SetY1(0.65)
        plot.plot("NOSTACK HISTE")
        bits = []
        if plot_reco:
            bits.append("reco")
        if plot_gen:
            bits.append("gen")
        if div_bin_width:
            bits.append("div_bin_width")
        content = "_".join(bits)
        output_filename = os.path.join(plot_dir, "%s_%s_1d_%s.%s" % (var_dict['name'], append, content, OUTPUT_FMT))
        plot.save(output_filename)


def get_summed_hists(tfile, hist_names):
    first_hist = cu.get_from_tfile(tfile, hist_names[0])
    if len(hist_names) == 1:
        return first_hist
    else:
        first_hist = first_hist.Clone()
        for hname in hist_names:
            this_hist = cu.get_from_tfile(tfile, hname)
            first_hist.Add(this_hist)
        return first_hist


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
                        default=0.5,
                        type=float)
    parser.add_argument("--spikeSmoothing",
                        help="Smooth spike @ 0",
                        action='store_true')
    # acceptable_metrics = ['gausfit', 'quantile']
    # parser.add_argument("--metric",
    #                     help="Metric for deciding bin width.",
    #                     default=acceptable_metrics[0],
    #                     choices=acceptable_metrics)
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        raise IOError("Cannot find input file %s" % args.input)

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
        # Source directory(ies) for 2D hists - if more than one in the tuple,
        # it will auto sum over the releveant map across the dirs
        source_plot_dir_names = [("Dijet_QG_tighter",),
                                 ("Dijet_QG_central_tighter", "Dijet_QG_forward_tighter"),
                                 ("Dijet_QG_central_tighter_groomed", "Dijet_QG_forward_tighter_groomed")][1:]
        region_labels = [qgc.Dijet_LABEL,
                         qgc.Dijet_LABEL,
                         qgc.Dijet_GROOMED_LABEL][1:]

    elif "dyjetstoll" in args.input.lower():
        source_plot_dir_names = [("ZPlusJets_QG"), ("ZPlusJets_QG_groomed")]
        region_labels = [qgc.DY_ZpJ_LABEL]

    else:
        raise RuntimeError("No idea which region we're using")

    jet_str = "AK4 PUPPI"
    if "ak8puppi" in args.input:
        jet_str = "AK8 PUPPI"

    rebin_results_dict = OrderedDict()

    pt_regions = [
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
    ]

    do_rel_response = False  # Use relative response instead - only makes sense with quantiles/gaussian fit?

    input_tfile = cu.open_root_file(args.input)

    for angle in qgc.COMMON_VARS[:]:

        integer_binning = "multiplicity" in angle.var.lower()

        all_response_maps_dict = {}

        for source_plot_dir_name_list, region_label in zip(source_plot_dir_names, region_labels):
            var_prepend = "Groomed " if "groomed" in source_plot_dir_name_list[0] else ""
            total_plot_dir_name = "+".join(source_plot_dir_name_list)

            response_maps_dict = {} # store 2D maps for later

            for pt_region_dict in pt_regions[:]:
                var_dict = {
                    "name": "%s/%s%s" % (total_plot_dir_name, angle.var, pt_region_dict['append']),
                    "var_label": "%s%s (%s)" % (var_prepend, angle.name, angle.lambda_str),
                    "title": "%s jets: %s\n%s" % (jet_str, region_label, pt_region_dict['title']),
                }

                h_append = "_rel_response" if do_rel_response else "_response"
                names = ["%s/%s%s%s" % (sname, angle.var, pt_region_dict['append'], h_append) for sname in source_plot_dir_name_list]
                h2d_orig = get_summed_hists(input_tfile, names)

                # Make plots with original fine equidistant binning
                # -------------------------------------------------
                make_plots(h2d_orig, var_dict, plot_dir=plot_dir, append="orig",
                           plot_migrations=False) # don't do a div_bin_width=False version, as equidistant binning

                response_maps_dict[var_dict['name']] = h2d_orig

                # metric = "gausfit"
                # metric = "quantile"
                # new_binning = calc_variable_binning_metric(h2d_orig, plot_dir, args.metric)

                # Calculate new binning
                # ---------------------
                # goal = 0.5
                goal = args.target
                print("Calculating binning for", var_dict['name'])
                new_binning = calc_variable_binning_purity_stability(h2d_orig,
                                                                     purity_goal=goal,
                                                                     stability_goal=goal,
                                                                     integer_binning=integer_binning,
                                                                     spike_smoothing='charged' in angle.var and args.spikeSmoothing)
                rebin_results_dict[var_dict['name']] = new_binning

                # Rebin 2D heatmap
                # ----------------
                h2d_rebin = make_rebinned_2d_hist(h2d_orig, new_binning, use_half_width_y=False)

                # Plot with new binning
                # ---------------------
                make_plots(h2d_rebin, var_dict, plot_dir=plot_dir, append="rebinned")
                make_plots(h2d_rebin, var_dict, plot_dir=plot_dir, append="rebinned",
                           plot_maps=False,
                           plot_migrations=False,
                           div_bin_width=False)

                # Cache renormed plots here for migration plots
                h2d_renorm_x = cu.make_normalised_TH2(h2d_rebin, 'X', recolour=False, do_errors=True)
                output_tfile.WriteTObject(h2d_renorm_x)  # we want renormalised by col for doing folding
                response_maps_dict[var_dict['name']+"_renormX"] = h2d_renorm_x

                h2d_renorm_y = cu.make_normalised_TH2(h2d_rebin, 'Y', recolour=False, do_errors=True)
                output_tfile.WriteTObject(h2d_renorm_y)  # save anyway for completeness
                response_maps_dict[var_dict['name']+"_renormY"] = h2d_renorm_y

                # Now rebin any other input files with the same hist using the new binning
                # ------------------------------------------------------------------------
                if args.rebinThisInput and len(args.rebinThisInput) > 0:
                    contributions = qgp.migration_plot_components(h2d_renorm_x, h2d_renorm_y, var_dict['var_label'])

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

                        h2d_renorm_x_other = cu.make_normalised_TH2(h2d_rebin_other, 'X', recolour=False, do_errors=True)
                        tfile_other_out.WriteTObject(h2d_renorm_x_other)
                        h2d_renorm_y_other = cu.make_normalised_TH2(h2d_rebin_other, 'Y', recolour=False, do_errors=True)
                        contributions_other = qgp.migration_plot_components(h2d_renorm_x_other, h2d_renorm_y_other, var_dict['var_label'])
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

            all_response_maps_dict.update(response_maps_dict)

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
                        "var_label": "%s%s (%s)" % (var_prepend, angle.name, angle.lambda_str),
                        "title": "#splitline{%s}{(rebinned for %s)}" % (h2d.GetTitle(), this_pt_dict['title']),
                    }
                    make_plots(h2d_rebin, var_dict, plot_dir=plot_dir,
                               append="rebinned_for%s" % (reference_pt_region),
                               plot_migrations=True)
                    make_plots(h2d_rebin, var_dict, plot_dir=plot_dir,
                               append="rebinned_for%s" % (reference_pt_region),
                               plot_maps=False,
                               plot_migrations=False,
                               div_bin_width=False)

                # Apply binning scheme derived from one ungroomed region to groomed
                # And use the reference pt region
                # ---------------------------------------------------------
                rebin_other_groomed_regions = True
                if "groomed" in source_plot_dir_name_list[0] and rebin_other_groomed_regions:
                    # Find rebinning scheme
                    ungroomed_source_plot_dir_name = "+".join([s.replace("_groomed", "") for s in source_plot_dir_name_list])
                    key = "%s/%s%s" % (ungroomed_source_plot_dir_name, angle.var, reference_pt_region)
                    print("key", key)
                    if key not in rebin_results_dict:
                        raise RuntimeError("Missing key %s in dict" % key)

                    reference_binning = rebin_results_dict[key]

                    # Now rebin all other 2D hists
                    for hname, h2d in all_response_maps_dict.items():
                        # print("hname", hname)
                        # print("source_plot_dir_name", source_plot_dir_name)
                        # print("angle.var", angle.var)
                        if ("_renorm" in hname
                            or angle.var not in hname
                            or source_plot_dir_name_list[0] not in hname):
                            continue

                        print("Rebinning groomed", hname)
                        h2d_rebin = make_rebinned_2d_hist(h2d, reference_binning, use_half_width_y=False)
                        var_dict = {
                            "name": hname,
                            "var_label": "Groomed %s (%s)" % (angle.name, angle.lambda_str),
                            "title": "#splitline{%s}{(rebinned for %s, ungroomed)}" % (h2d.GetTitle(), this_pt_dict['title']),
                        }
                        make_plots(h2d_rebin, var_dict, plot_dir=plot_dir,
                                   append="rebinned_for%s_ungroomed" % (reference_pt_region),
                                   plot_migrations=True)
                        make_plots(h2d_rebin, var_dict, plot_dir=plot_dir,
                                   append="rebinned_for%s_ungroomed" % (reference_pt_region),
                                   plot_maps=False,
                                   plot_migrations=False,
                                   div_bin_width=False)

                # Now rebin any other input files with the same hist using the new binning
                # ------------------------------------------------------------------------
                if args.rebinThisInput and len(args.rebinThisInput) > 0:
                    # contributions = qgp.migration_plot_components(h2d_renorm_x, h2d_renorm_y, var_dict['var_label'])

                    lens = [len(x) for x in [args.rebinThisInput, args.rebinThisLabel, args.rebinThisOutputFile]]
                    if any(l != lens[0] for l in lens[1:]):
                        raise RuntimeError("Number of --rebinThisInput, --rebinThisLabel, --rebinThisOutputFile should be the same - need one of each per extra input file")


                    for pt_region_dict in pt_regions[:]:
                        var_dict = {
                            "name": "%s/%s%s" % (total_plot_dir_name, angle.var, pt_region_dict['append']),
                            "var_label": "%s%s (%s)" % (var_prepend, angle.name, angle.lambda_str),
                            "title": "%s\n%s\n(rebinned for %s)" % (region_label, pt_region_dict['title'], this_pt_dict['title']),
                        }
                        contributions = []

                        print(var_dict)
                        h_append = "_rel_response" if do_rel_response else "_response"
                        names = ["%s/%s%s%s" % (sname, angle.var, pt_region_dict['append'], h_append) for sname in source_plot_dir_name_list]

                        for ind, (other_input, other_label, other_output_filename) in enumerate(zip(args.rebinThisInput, args.rebinThisLabel, args.rebinThisOutputFile)):
                            print(other_label)
                            tfile_other = cu.open_root_file(other_input)
                            tfile_other_out = cu.open_root_file(other_output_filename, 'RECREATE')

                            h2d_orig = get_summed_hists(tfile_other, names)
                            h2d_rebin_other = make_rebinned_2d_hist(h2d_other, reference_binning, use_half_width_y=False)
                            print(h2d_other)
                            print(h2d_rebin_other)

                            this_var_dict = deepcopy(var_dict)
                            this_var_dict['title'] += "\n[%s]" % other_label
                            make_plots(h2d_rebin_other, this_var_dict,
                                       plot_dir=plot_dir+"_"+other_label.replace(" ", "_"),
                                       append="rebinned_for%s" % (reference_pt_region),
                                       plot_migrations=True)
                            make_plots(h2d_rebin_other, this_var_dict,
                                       plot_dir=plot_dir+"_"+other_label.replace(" ", "_"),
                                       append="rebinned_for%s" % (reference_pt_region),
                                       plot_migrations=False,
                                       plot_maps=False,
                                       div_bin_width=False)

                            h2d_renorm_x_other = cu.make_normalised_TH2(h2d_rebin_other, 'X', recolour=False, do_errors=True)
                            tfile_other_out.WriteTObject(h2d_renorm_x_other)
                            h2d_renorm_y_other = cu.make_normalised_TH2(h2d_rebin_other, 'Y', recolour=False, do_errors=True)
                            contributions_other = qgp.migration_plot_components(h2d_renorm_x_other, h2d_renorm_y_other, var_dict['var_label'])
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

