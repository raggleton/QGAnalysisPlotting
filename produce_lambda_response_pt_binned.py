#!/usr/bin/env python


"""Do lambda reponse/binning as a function of pt"""


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
from comparator import Contribution, Plot
import common_utils as cu

ROOT.gStyle.SetPaintTextFormat(".3f")

# Control output format
OUTPUT_FMT = "pdf"



# def renorm(arr2d, axis):
#     # create version where each axis summed to 1
#     # use where and out args to ensure nans are made into 0s
#     summed = arr2d.sum(axis=axis, keepdims=True)
#     return np.divide(arr2d, summed, where=summed!=0, out=np.zeros_like(arr2d))


# def concat_row(arr2d, row_ind):
#     # concat row row_ind + row_ind+1
#     nrows, ncols = arr2d.shape
#     if row_ind > nrows - 2:
#         raise RuntimeError("Cannot concat row [%d] as only %d rows in matrix" % (row_ind, nrows))
#     arr2d_new = np.zeros(shape=(nrows-1, ncols), dtype=float)
#     new_row = arr2d[row_ind] + arr2d[row_ind+1]
#     arr2d_new[row_ind] = new_row
#     # fill in new matrix
#     if row_ind > 0:
#         # do that bit before the new row
#         arr2d_new[:row_ind, ] = arr2d[:row_ind, ]
#     if row_ind < nrows - 2:
#         arr2d_new[row_ind+1:, :] = arr2d[row_ind+2:, :]
#     return arr2d_new


# def calc_variable_binning_other(h2d):
#     """here we just keep combining neighbouring bins until we reach desired purity & stability"""
#     arr2d, _ = cu.th2_to_arr(h2d)
#     reco_bin_edges = cu.get_bin_edges(h2d, 'Y')
#     gen_bin_edges = cu.get_bin_edges(h2d, 'X')

#     new_bin_edges = np.array(reco_bin_edges).reshape(1, len(reco_bin_edges))

#     # assumes both axes have same dimension!
#     bin_ind = 0
#     counter = 0  # safety measure
#     while bin_ind < len(arr2d)-1 and counter < 10000:
#         counter += 1
#         arr2d_renormx = renorm(arr2d, axis=0)
#         arr2d_renormy = renorm(arr2d, axis=1)
#         purity = arr2d_renormy[bin_ind][bin_ind]
#         stability = arr2d_renormx[bin_ind][bin_ind]
#         if purity > 0.4 and stability > 0.4:
#             # print("found bin")
#             print("bin_ind:", bin_ind, "purity: %.3f" % purity, "stability: %.3f" % stability)
#             bin_ind += 1
#             continue
#         else:
#             # print("combining bin", bin_ind, "/", len(arr2d))
#             # combine rows & columns (use transpose for latter)
#             arr2d = concat_row(arr2d, bin_ind)
#             arr2d = concat_row(arr2d.T, bin_ind).T
#             new_bin_edges = np.delete(new_bin_edges, bin_ind+1)  # keep track of new binning
#             continue

#     # keep [1] to be same as next [0], otherwise you lose a bin later when
#     # making the rebinned TH2
#     these_bins = [list(x) for x in zip(new_bin_edges[:-1], new_bin_edges[1:])]
#     # manual hack for last bin
#     these_bins[-1][1] = reco_bin_edges[-1]
#     print(these_bins)
#     return these_bins


# def rebin_2d_hist(h2d, new_binning_x, new_binning_y):
#     bin_edges_x = [b[0] for b in new_binning_x]
#     bin_edges_x.append(new_binning_x[-1][1])

#     bin_edges_y = [b[0] for b in new_binning_y]
#     bin_edges_y.append(new_binning_y[-1][1])

#     print("rebin_2d_hist, new axes:", bin_edges_x, bin_edges_y)

#     new_hist = ROOT.TH2D(
#         h2d.GetName()+"Rebin",
#         ';'.join([h2d.GetTitle(), h2d.GetXaxis().GetTitle(), h2d.GetYaxis().GetTitle()]),
#         len(new_binning_x),
#         array('d', bin_edges_x),
#         len(new_binning_y),
#         array('d', bin_edges_y)
#     )

#     nbins_x_orig = h2d.GetNbinsX()
#     bins_x_orig = cu.get_bin_edges(h2d_orig, 'X')

#     nbins_y_orig = h2d.GetNbinsY()
#     bins_y_orig = cu.get_bin_edges(h2d_orig, 'Y')

#     # TODO: handle under/overflow
#     for xind, xbin in enumerate(new_binning_x, 1):  # these are tuples
#         for yind, ybin in enumerate(new_binning_y, 1):
#             # Find all the old bins that correspond to this new bin, get their contents
#             new_bin_content = 0
#             new_bin_error = 0
#             # TODO handle error - reset to sqrt(N)?
#             for xind_orig, xbin_orig in enumerate(bins_x_orig[:-1], 1):
#                 for yind_orig, ybin_orig in enumerate(bins_y_orig[:-1], 1):
#                     if (xbin_orig >= xbin[0] and xbin_orig < xbin[1] and
#                         ybin_orig >= ybin[0] and ybin_orig < ybin[1]):
#                         new_bin_content += h2d.GetBinContent(xind_orig, yind_orig)
#                         # print("For new bin", xbin, ybin, "using contents from", xbin_orig, bins_x_orig[xind_orig], ybin_orig, bins_y_orig[yind_orig])
#                         # print("orig bin", xind_orig, yind_orig)

#             new_hist.SetBinContent(xind, yind, new_bin_content)
#             # print("Setting bin", xind, yind)
#             new_hist.SetBinError(xind, yind, new_bin_error)

#     return new_hist


# def make_rebinned_plot(h2d, new_binning, use_half_width_y=False):
#     # define reco binning as being half of gen binning
#     if use_half_width_y:
#         reco_binning = []
#         reco_bin_edges = cu.get_bin_edges(h2d, 'Y')
#         for s, e in new_binning:
#             ideal_mid = (s+e)/2.
#             mid = reco_bin_edges[bisect.bisect_left(reco_bin_edges, ideal_mid)]
#             reco_binning.append((s, mid))
#             reco_binning.append((mid, e))
#         return rebin_2d_hist(h2d, new_binning, reco_binning)
#     else:
#         return rebin_2d_hist(h2d, new_binning, new_binning)


def make_plots(h2d, var_dict, plot_dir, append="", plot_migrations=True):
    """Plot a 2D hist, with copies renormalised by row and by column.
    Also optionally plot migrations as 1D plot.
    """
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
        h2d.SetTitle(var_dict['title'])
    h2d.Draw("COLZ")
    output_filename = os.path.join(plot_dir, var_dict['name']+"_%s.%s" % (append, OUTPUT_FMT))
    output_dir = os.path.dirname(output_filename)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    canv.SaveAs(output_filename)

    canv.SetLogz()
    output_filename = os.path.join(plot_dir, var_dict['name']+"_%s_logZ.%s" % (append, OUTPUT_FMT))
    canv.SaveAs(output_filename)

    # renorm by row
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

    # renorm by column
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

    # Plot migrations
    if plot_migrations:
        output_filename = os.path.join(plot_dir, "%s_migration_summary.%s" % (var_dict['name'], OUTPUT_FMT))
        qgg.make_migration_summary_plot(h2d_renorm_x, h2d_renorm_y, var_dict['var_label'], output_filename)


def extract_parts_from_name(name):
    parts = name.split("_")
    stem, ind0, ind1 = parts[:-2], parts[-2], parts[-1]
    return ['_'.join(stem), int(ind0), int(ind1)]


def sum_over_inds(tobj, template, bins):

    final = None
    for b in bins:
        if final is None:
            final = cu.get_from_tfile(tobj, template.format(b)).Clone(template.format("all"))
        else:
            final.Add(cu.get_from_tfile(tobj, template.format(b)))
    return final


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        help='Input ROOT file to process.'
                        'Several dirs can be specified here, separated by a space.')
    parser.add_argument("--outputDir",
                        help="Directory to put output plot dirs into",
                        default=None)
    parser.add_argument("--outputFile",
                        help="Output ROOT file for rebinned 2D hists",
                        default=None)
    args = parser.parse_args()

    input_dir, input_basename = os.path.split(args.input)
    default_plot_dir = os.path.join(input_dir, "rebinning_"+os.path.splitext(input_basename)[0])
    plot_dir = args.outputDir if args.outputDir else default_plot_dir

    default_output_root_filename = os.path.join(input_dir, "rebinned_" + input_basename)
    output_root_filename = args.outputFile if args.outputFile else default_output_root_filename
    output_tfile = cu.open_root_file(output_root_filename, 'RECREATE')

    source_plot_dir_name = None
    if "qcd" in args.input.lower():
        source_plot_dir_name = "Dijet_QG_tighter"

    if "dyjetstoll" in args.input.lower():
        source_plot_dir_name = "ZPlusJets_QG"

    do_these = [
        {
            "name": "jet_puppiMultiplicity",
            "var_label": "PUPPI Multiplicity (#lambda_{0}^{0} (PUPPI))",
            "log": True,
        },
        {
            "name": "jet_lha",
            "var_label": "LHA (#lambda_{0.5}^{1})"
        },
        {
            "name": "jet_pTD",
            "var_label": "p_{T}^{D} (#lambda_{0}^{2})"
        },
        {
            "name": "jet_width",
            "var_label": "Width (#lambda_{1}^{1})"
        },
        {
            "name": "jet_thrust",
            "var_label": "Thrust (#lambda_{2}^{1})"
        },
    ][:]

    rebin_results_dict = OrderedDict()

    for var_dict in do_these:
        do_rel_response = False
        full_var_name = var_dict['name']
        if do_rel_response:
            full_var_name += "_rel_response"
        else:
            full_var_name += "_response"

        print(full_var_name)

        tfile = cu.open_root_file(args.input)
        
        plot_names = cu.get_list_of_element_names(tfile.Get(source_plot_dir_name))
        plot_names = [p for p in plot_names if full_var_name+"_bin_" in p]
        # print(plot_names)
        
        var_stem = None
        first_ind = []
        second_ind = []
        for p in plot_names:
            stem, ind1, ind2 = extract_parts_from_name(p)
            if var_stem is None:
                var_stem = stem
            first_ind.append(ind1)
            second_ind.append(ind2)

        # First index is for gen pt bins, second is for reco pt bins
        first_ind = sorted(list(set(first_ind)))
        second_ind = sorted(list(set(second_ind)))

        # for binning by reco
        sum_over = first_ind
        bin_by = second_ind
        title_part = 1
        # for binning by gen
        sum_over = second_ind
        bin_by = first_ind
        title_part = 0

        for bin_ind in bin_by:
            template = "%s/%s_{}_%d" % (source_plot_dir_name, var_stem, bin_ind)
            template = "%s/%s_%d_{}" % (source_plot_dir_name, var_stem, bin_ind)
            # Sum over gen pt bins
            summed_hist = sum_over_inds(tfile, 
                                        template=template,
                                        bins=sum_over)
            this_var_dict = deepcopy(var_dict)
            this_var_dict['name'] = summed_hist.GetName()
            title = summed_hist.GetTitle()
            new_title = title.split(",")[title_part].strip()
            this_var_dict['title'] = new_title
            make_plots(summed_hist, this_var_dict, plot_dir=plot_dir, append="orig", plot_migrations=False)
            summed_hist.SetName(this_var_dict['name'].split("/")[-1])
            output_tfile.WriteTObject(summed_hist)
    
        # new_binning = calc_variable_binning_other(h2d_orig)
        # rebin_results_dict[var_dict['name']] = new_binning

        # h2d_rebin = make_rebinned_plot(h2d_orig, new_binning, use_half_width_y=False)

        # make_plots(h2d_rebin, var_dict, plot_dir=plot_dir, append="rebinned", plot_migrations=True)

        # h2d_renorm_x = cu.make_normalised_TH2(h2d_rebin, 'X', recolour=False, do_errors=False)
        # output_tfile.WriteTObject(h2d_renorm_x)  # we want renormalised by col for doing folding
        
        # h2d_renorm_y = cu.make_normalised_TH2(h2d_rebin, 'Y', recolour=False, do_errors=False)

        # contributions = qgg.migration_plot_components(h2d_renorm_x, h2d_renorm_y, var_dict['var_label'])


    output_tfile.Close()

    # # Save new binning to txt file
    # output_txt = os.path.splitext(args.input)[0] + ".txt"
    # parts = os.path.split(output_txt)
    # output_txt = os.path.join(input_dir, 'binning_'+parts[1])
    # with open(output_txt, 'w') as fout:
    #     for k, v in rebin_results_dict.items():
    #         fout.write("%s: %s\n" % (k, v))

    # print("saved new binning to", output_txt)
    # print("saved rebinned 2D hists to", output_root_filename)

