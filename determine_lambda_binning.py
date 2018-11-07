#!/usr/bin/env python

"""Figure out binnining for lambda variables based on width, etc"""


import argparse
# from MyStyle import My_Style
# My_Style.cd()
import os
from itertools import product
from array import array
from math import sqrt
import bisect
import numpy as np


import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
# ROOT.gStyle.SetOptStat(0)
# ROOT.gStyle.SetOptFit(1)


# My stuff
from comparator import Contribution, Plot
import qg_common as qgc
import qg_general_plots as qgg
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

        print("Mean:", mean, "width:", width, "fit" if fit_res.Status() == 0 else "Raw")


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
            desired_width = 0.85
            quantiles = array('d', [(1-desired_width)/2., 1-((1-desired_width))])
            values = array('d', [0, 0])
            hproj.GetQuantiles(len(quantiles), values, quantiles)
            quantile_ideal_lower_edge, quantile_ideal_upper_edge = values
            ideal_upper_edge = max(quantile_ideal_upper_edge, ideal_upper_edge)

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
                raise RuntimeError("new_bin_end == bin_start, got stuck")

            if new_bin_end >= len(reco_bin_edges):
                new_bin_end = len(reco_bin_edges)-1

            upper_edge = reco_bin_edges[new_bin_end]
            bins.append((lower_edge, upper_edge))
            print("Adding bin", bins[-1])
            
            c = ROOT.TCanvas("c"+cu.get_unique_str(), "", 800, 600)
            hproj.Draw()
            output_filename = os.path.join(plot_dir, h2d.GetName() + "_bin%dto%d.%s" % (bin_start, bin_end, OUTPUT_FMT))
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


def calc_variable_binning_other(h2d):
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
        if purity > 0.4 and stability > 0.4:
            print("found bin")
            print("bin_ind:", bin_ind, "purity: %.3f" % purity, "stability: %.3f" % stability)
            bin_ind += 1
            continue
        else:
            print("combining bin", bin_ind, "/", len(arr2d))
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
    bin_edges_x = [b[0] for b in new_binning_x]
    bin_edges_x.append(new_binning_x[-1][1])

    bin_edges_y = [b[0] for b in new_binning_y]
    bin_edges_y.append(new_binning_y[-1][1])
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


def make_rebinned_plot(h2d, new_binning, use_half_width_y=False):
    # define reco binning as being half of gen binning
    if use_half_width_y:
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        nargs='+',
                        help='Input ROOT files to process. '
                        'Several dirs can be specified here, separated by a space.')
    parser.add_argument("-o", "--output", help="Directory to put output plot dirs into", default=None)
    acceptable_metrics = ['gausfit', 'quantile']
    parser.add_argument("--metric", help="Metric for deciding bin width.",
                        default=acceptable_metrics[0],
                        choices=acceptable_metrics)
    args = parser.parse_args()

    for in_file in args.input:
        default_plot_dir = os.path.join(os.path.dirname(in_file), "rebinning_"+os.path.splitext(os.path.basename(in_file))[0])
        plot_dir = args.output if args.output else default_plot_dir
        # cu.check_dir_exists_create(plot_dir)

        do_these = [
            {
                "name": "Dijet_QG_tighter/jet_puppiMultiplicity",
                "var_label": "PUPPI Multiplicity (#lambda_{0}^{0} (PUPPI))",
                "log": True,
            },
            {
                "name": "Dijet_QG_tighter/jet_LHA",
                "var_label": "LHA (#lambda_{0.5}^{1})"
            },
            {
                "name": "Dijet_QG_tighter/jet_pTD",
                "var_label": "p_{T}^{D} (#lambda_{0}^{2})"
            },
            {
                "name": "Dijet_QG_tighter/jet_width",
                "var_label": "Width (#lambda_{1}^{1})"
            },
            {
                "name": "Dijet_QG_tighter/jet_thrust",
                "var_label": "Thrust (#lambda_{2}^{1})"
            },
            # charged vars
            {
                "name": "Dijet_QG_tighter/jet_puppiMultiplicity_charged",
                "var_label": "PUPPI Multiplicity (#lambda_{0}^{0} (PUPPI)) [charged]",
                "log": True,
            },
            {
                "name": "Dijet_QG_tighter/jet_LHA_charged",
                "var_label": "LHA (#lambda_{0.5}^{1}) [charged only]"
            },
            {
                "name": "Dijet_QG_tighter/jet_pTD_charged",
                "var_label": "p_{T}^{D} (#lambda_{0}^{2}) [charged only]"
            },
            {
                "name": "Dijet_QG_tighter/jet_width_charged",
                "var_label": "Width (#lambda_{1}^{1}) [charged only]"
            },
            {
                "name": "Dijet_QG_tighter/jet_thrust_charged",
                "var_label": "Thrust (#lambda_{2}^{1}) [charged only]"
            },
        ][:2]

        for var_dict in do_these:
            do_rel_response = False
            full_var_name = var_dict['name']
            if do_rel_response:
                full_var_name += "_rel_response"
            else:
                full_var_name += "_response"

                print(full_var_name)

            tfile = cu.open_root_file(in_file)
            h2d_orig = cu.get_from_tfile(tfile, full_var_name)

            # metric = "gausfit"
            # metric = "quantile"
            # new_binning = calc_variable_binning(h2d_orig, plot_dir, args.metric)

            new_binning = calc_variable_binning_other(h2d_orig)


            h2d_rebin = make_rebinned_plot(h2d_orig, new_binning, use_half_width_y=False)

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
            h2d_rebin.SetTitleOffset(xtitle_offset, 'X')
            h2d_rebin.SetTitleOffset(ytitle_offset, 'Y')
            h2d_rebin.Draw("COLZ")
            output_filename = os.path.join(plot_dir, var_dict['name']+"_rebinned.%s" % (OUTPUT_FMT))
            output_dir = os.path.dirname(output_filename)
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir)
            canv.SaveAs(output_filename)

            canv.SetLogz()
            output_filename = os.path.join(plot_dir, var_dict['name']+"_rebinned_logZ.%s" % (OUTPUT_FMT))
            canv.SaveAs(output_filename)

            # renorm by row
            canv.SetLogz(0)
            h2d_renorm_y = cu.make_normalised_TH2(h2d_rebin, 'Y', recolour=False, do_errors=False)
            marker_size = 0.8
            h2d_renorm_y.SetMarkerSize(0.5)
            h2d_renorm_y.SetMarkerSize(marker_size)
            h2d_renorm_y.SetMaximum(1)
            draw_opt = "COLZ TEXT45"
            h2d_renorm_y.Draw(draw_opt)
            xtitle_offset = 1.5
            h2d_renorm_y.SetTitleOffset(xtitle_offset, 'X')

            yax = h2d_renorm_y.GetYaxis()
            h2d_renorm_y.SetTitleOffset(ytitle_offset, 'Y')
            canv.Update()

            canv.SaveAs(os.path.join(plot_dir, "%s_rebinned_renormY_linZ.%s" % (var_dict['name'], OUTPUT_FMT)))
            canv.SetLogz()
            h2d_renorm_y.SetMaximum(1)
            h2d_renorm_y.SetMinimum(1E-3)
            canv.SaveAs(os.path.join(plot_dir, "%s_rebinned_renormY_logZ.%s" % (var_dict['name'], OUTPUT_FMT)))

            # renorm by column
            canv.Clear()
            canv.SetLogz(0)
            h2d_renorm_x = cu.make_normalised_TH2(h2d_rebin, 'X', recolour=False, do_errors=False)
            h2d_renorm_x.SetMarkerSize(0.5)
            h2d_renorm_x.SetMarkerSize(marker_size)
            h2d_renorm_x.SetMaximum(1)
            h2d_renorm_x.Draw(draw_opt)

            h2d_renorm_x.SetTitleOffset(xtitle_offset, 'X')

            yax = h2d_renorm_x.GetYaxis()
            h2d_renorm_x.SetTitleOffset(ytitle_offset, 'Y')
            yax.SetMoreLogLabels()
            canv.Update()
            canv.SaveAs(os.path.join(plot_dir, "%s_rebinned_renormX_linZ.%s" % (var_dict['name'], OUTPUT_FMT)))
            canv.SetLogz()
            h2d_renorm_x.SetMaximum(1)
            h2d_renorm_x.SetMinimum(1E-3)
            canv.SaveAs(os.path.join(plot_dir, "%s_rebinned_renormX_logZ.%s" % (var_dict['name'], OUTPUT_FMT)))

            # Plot migrations
            output_filename = os.path.join(plot_dir, "%s_migration_summary.%s" % (var_dict['name'], OUTPUT_FMT))
            qgg.make_migration_summary_plot(h2d_renorm_x, h2d_renorm_y, var_dict['var_label'], output_filename)
