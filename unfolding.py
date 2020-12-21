#!/usr/bin/env python


"""TUnfold it all

Thanks to Ashley, Dennis

Can profile memory usage with memory_profiler package:

mprof run --interval 0.1 --python unfolding.py <args>
"""


from __future__ import print_function, division

import os
import sys
from array import array
import numpy as np
import math
import json
from itertools import product, chain
from copy import copy, deepcopy
from pprint import pprint
from bisect import bisect_right
import warnings
from functools import partial

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot

# my packages
import common_utils as cu
import qg_common as qgc
from my_unfolder import MyUnfolder, pickle_region, unpickle_region, ExpSystematic, TruthTemplateMaker, BinningHandler, PtVarBinning, PtVarPerPtBinning, InputHandler
from my_unfolder_plotter import MyUnfolderPlotter
from unfolding_regularisation_classes import TauScanner, LCurveScanner
from unfolding_config import get_dijet_config, get_zpj_config
from unfolding_logistics import get_unfolding_argparser, get_unfolding_output_dir, sanitise_args, AREA_OPT_DICT
from do_unfolding_plots import Setup, do_binned_plots_per_region_angle, do_all_big_normalised_1d_plots_per_region_angle

# Use rootpy to throw exceptions on ROOT errors, but need DANGER enabled
# Doesn't work with python 3.8 for now
if not (sys.version_info.major == 3 and sys.version_info.minor >= 8):
    import rootpy
    # import rootpy.logger.magic as M; M.DANGER.enabled = True


# monkey-patch warning formatter
warnings.formatwarning = cu._formatwarning

My_Style.cd()

ROOT.gErrorIgnoreLevel = ROOT.kWarning
# ROOT.gErrorIgnoreLevel = ROOT.kInfo
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)  # VERY IMPORTANT - somewhere, closing a TFile for exp systs deletes a map...dunno why

os.nice(10)

# Control plot output format
OUTPUT_FMT = "pdf"


# When using memory_profiler/mprof, handy to have @profile to mark functions
# But when running normally we want to pass through without manually commenting out,
# so define our own decorator instead that does nothing
# Note that if this module is imported elsewhere,
# profile will be in __builtins__, not directly in locals(),
# so we check multiple places
if not any(['profile' in locals(),
            'profile' in globals(),
            (isinstance(globals()['__builtins__'], dict) and 'profile' in globals()['__builtins__']),
            'profile' in dir(__builtins__)]):

    print("I have no memory_profiler @profile decorator in unfolding, creating my own instead")

    def profile(func):
        return func


def rm_large_rel_error_bins_th1(hist, relative_err_threshold=-1):
    """Reset bins in 1D hist to 0 if error/contents exceeds a certain value"""
    if relative_err_threshold < 0:
        return hist
    new_hist = hist.Clone()
    new_hist.SetDirectory(0)
    for ix in range(0, hist.GetNbinsX()+2):
        val = hist.GetBinContent(ix)
        err = hist.GetBinError(ix)
        if val == 0:
            continue
        if abs(1.*err/val) > relative_err_threshold:
            new_hist.SetBinContent(ix, 0)
            new_hist.SetBinError(ix, 0)
    return new_hist


def rm_large_rel_error_bins_th2(hist, relative_err_threshold=-1):
    """Reset bins in 2D hist to 0 if error/contents exceeds a certain value"""
    if relative_err_threshold < 0:
        return hist
    new_hist = hist.Clone()
    new_hist.SetDirectory(0)
    for ix in range(0, hist.GetNbinsX()+2):
        for iy in range(0, hist.GetNbinsY()+2):
            val = hist.GetBinContent(ix, iy)
            err = hist.GetBinError(ix, iy)
            if val == 0:
                continue
            if (abs(1.*err/val)) > relative_err_threshold:
                new_hist.SetBinContent(ix, iy, 0)
                new_hist.SetBinError(ix, iy, 0)
    return new_hist


def draw_projection_comparison(h_orig, h_projection, title, xtitle, output_filename, do_bin_comparison=True):
    """Draw 2 hists, h_orig the original, and h_projection the projection of a 2D hist

    Also compares total integrals, and bin-by-bin (if do_bin_comparison = True)
    """
    entries = [
        Contribution(h_orig, label="1D hist",
                     line_color=ROOT.kBlue, line_width=1,
                     marker_color=ROOT.kBlue, marker_size=0,
                     normalise_hist=False),
        Contribution(h_projection, label="Response map projection",
                     line_color=ROOT.kRed, line_width=1,
                     marker_color=ROOT.kRed, marker_size=0,
                     normalise_hist=False,
                     subplot=h_orig),
    ]
    plot = Plot(entries,
                what='hist',
                title=title,
                xtitle=xtitle,
                ytitle="N",
                subplot_type='ratio',
                subplot_title='Projection / 1D',
                subplot_limits=(0.999, 1.001)
                )
    plot.default_canvas_size = (800, 600)
    plot.plot("NOSTACK HIST")
    plot.main_pad.SetLogy(1)
    ymax = max(h.GetMaximum() for h in [h_orig, h_projection])
    plot.container.SetMaximum(ymax * 10)
    # plot.container.SetMinimum(1E-8)
    plot.legend.SetY1NDC(0.77)
    plot.legend.SetX2NDC(0.85)
    plot.save(output_filename)

    # Check integrals
    int_orig = h_orig.Integral()
    int_proj = h_projection.Integral()
    if abs(int_orig - int_proj)/int_orig > 0.01:
        warnings.warn(cu.pcolors.WARNING + "draw_projection_comparison: different integrals: %f vs %f" % (int_orig, int_proj) + cu.pcolors.ENDC)

    # Check bin-by-bin
    if do_bin_comparison:
        for i in range(1, h_orig.GetNbinsX()+1):
            value_orig = h_orig.GetBinContent(i)
            value_proj = h_projection.GetBinContent(i)
            if value_orig == 0 and value_proj == 0:
                continue
            rel_diff = abs((value_orig - value_proj)/max(abs(value_orig), abs(value_proj)))
            if rel_diff > 1E-5:
                # print("draw_projection_comparison: bin %s has different contents: %f vs %f (rel diff %f)" % (i, value_orig, value_proj, rel_diff))
                raise ValueError(
                    "draw_projection_comparison: bin %s has different contents: "
                    "hist: %f vs projection: %f (abs diff %f, rel diff %f)" % (
                     i, value_orig, value_proj, value_orig - value_proj, rel_diff))


# def fill_empty_bins(response_map,
#                     variable_bin_edges_reco,
#                     variable_bin_edges_gen,
#                     variable_name,
#                     pt_bin_edges_reco,
#                     pt_bin_edges_gen,
#                     pt_bin_edges_underflow_reco,
#                     pt_bin_edges_underflow_gen):
#     """Fill in empty bins in the response map with some small value

#     Testing to see if it helps stabilise unfolding
#     """
#     new_map = response_map.Clone()

#     generator_binning, detector_binning = MyUnfolder.construct_tunfold_binning(variable_bin_edges_reco,
#                                                                                variable_bin_edges_gen,
#                                                                                variable_name,
#                                                                                pt_bin_edges_reco,
#                                                                                pt_bin_edges_gen,
#                                                                                pt_bin_edges_underflow_reco,
#                                                                                pt_bin_edges_underflow_gen)

#     generator_binning_uflow = generator_binning.FindNode("generatordistribution_underflow")
#     generator_binning_main = generator_binning.FindNode("generatordistribution")

#     detector_binning_uflow = detector_binning.FindNode("detectordistribution_underflow")
#     detector_binning_main = detector_binning.FindNode("detectordistribution")

#     # in each gen pt, detector pt bin, find the largest and smallest bin counts
#     # then go through again and fill in the empty bins with some very small value,
#     # that is relatively smaller than the smallest bin is compared to the largest bin
#     # (we're just trying to ensure non-zero, not physically super sensible)
#     n_bins = (len(variable_bin_edges_gen)-1)*(len(variable_bin_edges_reco)-1)
#     for ibin_pt_gen, gen_pt in enumerate(chain(pt_bin_edges_underflow_gen[:-1], pt_bin_edges_gen)):
#         for ibin_pt_reco, reco_pt in enumerate(chain(pt_bin_edges_underflow_reco[:-1], pt_bin_edges_reco)):
#             largest_bin_count = -99999
#             smallest_bin_count = 9E99
#             empty_bin_indices = []

#             for gen_var in variable_bin_edges_gen[:-1]:
#                 # urghhhhhh have to manually choose the TUnfoldBinning object
#                 this_gen_binning = generator_binning_uflow if gen_pt < pt_bin_edges_gen[0] else generator_binning_main
#                 gen_bin_num = this_gen_binning.GetGlobalBinNumber(gen_var+0.001, gen_pt+0.001)
#                 for reco_var in variable_bin_edges_reco[:-1]:
#                     this_reco_binning = detector_binning_uflow if reco_pt < pt_bin_edges_reco[0] else detector_binning_main
#                     reco_bin_num = this_reco_binning.GetGlobalBinNumber(reco_var+0.001, reco_pt+0.001)
#                     val = response_map.GetBinContent(gen_bin_num, reco_bin_num)
#                     if val == 0:
#                         empty_bin_indices.append([gen_bin_num, reco_bin_num])
#                     elif val < smallest_bin_count:
#                         smallest_bin_count = val
#                     elif val > largest_bin_count:
#                         largest_bin_count = val

#             # if ibin_pt_gen < 4 and ibin_pt_reco < 4:
#             #     print(ibin_pt_gen, ibin_pt_reco, ":", len(empty_bin_indices), n_bins, largest_bin_count, smallest_bin_count, smallest_bin_count/largest_bin_count)

#             is_near_diagonal = abs(ibin_pt_gen - (ibin_pt_reco/2)) < 3
#             if (largest_bin_count > 0) and (smallest_bin_count > 0) and (((1. * len(empty_bin_indices) / n_bins) < 0.75) or is_near_diagonal):
#                 ratio = smallest_bin_count / largest_bin_count
#                 new_val = ratio * ratio * smallest_bin_count
#                 print("Filling", ibin_pt_gen, ibin_pt_reco, smallest_bin_count * ratio)
#                 for x, y in empty_bin_indices:
#                     new_map.SetBinContent(x, y, new_val)
#                     new_map.SetBinError(x, y, 0)

#     return new_map


def calc_background(hist, bg_fraction_hist):
    """Calculate background from hist,
    given the background fractions in bg_fraction_hist"""
    if hist.GetNbinsX() != bg_fraction_hist.GetNbinsX():
        raise ValueError("Mimsmatch in number of bins in calc_background()")
    bg_hist = hist.Clone(cu.get_unique_str())
    bg_hist.Multiply(bg_fraction_hist)
    return bg_hist


def subtract_background(hist, bg_fraction_hist):
    """Calculate and subtract background from hist,
    given the background fractions in bg_fraction_hist

    Returns boths the bg-subtracted hist, and the calculated bg-hist
    """
    bg_hist = calc_background(hist, bg_fraction_hist)
    new_hist = hist.Clone(cu.get_unique_str())
    new_hist.Add(bg_hist, -1)
    return new_hist, bg_hist


def find_disconnected_output_bins(response_map):
    """Find columns in response map with no entries over any Y bins"""
    project_x = response_map.ProjectionX()  # assumed gen on x axis
    disconnected_bins = []
    for ix in range(1, project_x.GetNbinsX()+1):
        val = project_x.GetBinContent(ix)
        if val == 0:
            disconnected_bins.append(ix)
    return disconnected_bins


def find_disconnected_input_bins(response_map):
    """Find rows in response map with no entries over any X bins"""
    # assumed gen on x axis
    # don't include under/overflow gen bins
    # project_y = response_map.ProjectionY(cu.get_unique_str(), 1, response_map.GetNbinsX())
    project_y = response_map.ProjectionY()
    disconnected_bins = []
    for ix in range(1, project_y.GetNbinsX()+1):
        val = project_y.GetBinContent(ix)
        if val == 0:
            disconnected_bins.append(ix)
    return disconnected_bins


def merge_th1_bins(h, bin_list, new_bin_edges=None):
    """Merges each bin in bin_list with its leftwards neighbour. Can specify new bin edges.

    Parameters
    ----------
    h : ROOT.TH1
        Histogram to rebin
    bin_list : list[int]
        List of bin indices (1-index, since it's ROOT)
    new_bin_edges : None, optional
        New bin edges for returned TH1D.
        If None, uses bin edges from original hist, just deleting those that
        are in bin_list.

    Returns
    -------
    ROOT.TH1D
    """
    if not bin_list or len(bin_list) == 0:
        return h
    bin_edges_orig = cu.get_bin_edges(h, 'x')
    # print("orig", bin_edges_orig)
    if new_bin_edges is None:
        new_bin_edges = [x for i, x in enumerate(bin_edges_orig, 1)
                         if i not in bin_list]
    # print("new", new_bin_edges)
    h_new = ROOT.TH1D(h.GetName() + cu.get_unique_str(),
                      ";".join([h.GetTitle(), h.GetXaxis().GetTitle(), h.GetYaxis().GetTitle()]),
                      len(new_bin_edges)-1,
                      array('d', new_bin_edges))
    for ix in range(0, h.GetNbinsX()+2):
        val = h.GetBinContent(ix)
        err = h.GetBinError(ix)
        pos = bisect_right(bin_list, ix)
        ix_new = ix - pos
        # print("ix", ix, "ix_new", ix_new)
        h_new.SetBinContent(ix_new, val + h_new.GetBinContent(ix_new))
        h_new.SetBinError(ix_new, np.hypot(err, h_new.GetBinError(ix_new)))
    return h_new


def merge_th2_bins(h, bin_list_x, bin_list_y, new_bin_edges_x=None, new_bin_edges_y=None):
    """Merges each bin in bin_list with its leftwards (index-1) neighbour. Can specify new bin edges.

    Parameters
    ----------
    h : ROOT.TH2D
        Histogram to rebin
    bin_list_x : list[int]
        List of bin indices (1-index, since it's ROOT).
        If empty or None, no rebinning of this axis.
    bin_list_y : list[int]
        List of bin indices (1-index, since it's ROOT).
        If empty or None, no rebinning of this axis.
    new_bin_edges_x : ndarray, list[float], optional
        New x bin edges for returned TH2D.
        If None, uses bin edges from original hist, just deleting those that
        are in bin_list.
    new_bin_edges_y : ndarray, list[float], optional
        New y bin edges for returned TH2D.
        If None, uses bin edges from original hist, just deleting those that
        are in bin_list.

    Returns
    -------
    ROOT.TH2D
    """
    if not bin_list_x or len(bin_list_x) == 0:
        h_new = h
    else:
        # do x bin merging first
        bin_edges_x_orig = cu.get_bin_edges(h, 'x')
        bin_edges_y_orig = cu.get_bin_edges(h, 'y')
        # print("origx", bin_edges_x_orig)
        # print("bin_list_x", bin_list_x)
        if new_bin_edges_x is None:
            new_bin_edges_x = [x for i, x in enumerate(bin_edges_x_orig, 1)
                               if i not in bin_list_x]
        # print("newx", new_bin_edges_x)
        h_new = ROOT.TH2D(h.GetName() + cu.get_unique_str(),
                          ";".join([h.GetTitle(), h.GetXaxis().GetTitle(), h.GetYaxis().GetTitle()]),
                          len(new_bin_edges_x)-1,
                          array('d', new_bin_edges_x),
                          len(bin_edges_y_orig)-1,
                          array('d', bin_edges_y_orig))

        for iy in range(0, h.GetNbinsY()+2):
            for ix in range(0, h.GetNbinsX()+2):
                val = h.GetBinContent(ix, iy)
                err = h.GetBinError(ix, iy)
                pos = bisect_right(bin_list_x, ix)
                ix_new = ix - pos
                h_new.SetBinContent(ix_new, iy, val + h_new.GetBinContent(ix_new, iy))
                h_new.SetBinError(ix_new, iy, np.hypot(err, h_new.GetBinError(ix_new, iy)))

    if not bin_list_y or len(bin_list_y) == 0:
        return h_new

    # now do y bin merging
    bin_edges_x_orig = cu.get_bin_edges(h_new, 'x')
    # print("origx", bin_edges_x_orig)
    bin_edges_y_orig = cu.get_bin_edges(h, 'y')
    # print("origy", bin_edges_y_orig)
    if new_bin_edges_y is None:
        new_bin_edges_y = [y for i, y in enumerate(bin_edges_y_orig, 1)
                           if i not in bin_list_y]
    # print("newy", new_bin_edges_y)
    h_new2 = ROOT.TH2D(h.GetName() + cu.get_unique_str(),
                       ";".join([h.GetTitle(), h.GetXaxis().GetTitle(), h.GetYaxis().GetTitle()]),
                       len(bin_edges_x_orig)-1,
                       array('d', bin_edges_x_orig),
                       len(new_bin_edges_y)-1,
                       array('d', new_bin_edges_y))

    for ix in range(0, h_new.GetNbinsX()+2):
        for iy in range(0, h_new.GetNbinsY()+2):
            val = h_new.GetBinContent(ix, iy)
            err = h_new.GetBinError(ix, iy)
            pos = bisect_right(bin_list_x, iy)
            iy_new = iy - pos
            h_new2.SetBinContent(ix, iy_new, val + h_new2.GetBinContent(ix, iy_new))
            h_new2.SetBinError(ix, iy_new, np.hypot(err, h_new2.GetBinError(ix, iy_new)))
    return h_new2


def merge_th1_bin_pairs(h, bin_pairs, new_bin_edges):
    if not bin_pairs or len(bin_pairs) == 0:
        return h
    h_new = ROOT.TH1D(h.GetName() + cu.get_unique_str(),
                      ";".join([h.GetTitle(), h.GetXaxis().GetTitle(), h.GetYaxis().GetTitle()]),
                      len(new_bin_edges)-1,
                      array('d', new_bin_edges))

    for ix in range(0, h_new.GetNbinsX()+2):
        # copy original data first, but only the subset of overlapping bins
        h_new.SetBinContent(ix, h.GetBinContent(ix))
        h_new.SetBinError(ix, h.GetBinError(ix))

    # now do the extra bins that don't appear in h_new,
    # updating the bins in the new hist
    for ix_src, ix_dest in bin_pairs:
        h_new.SetBinContent(ix_dest, h_new.GetBinContent(ix_dest) + h.GetBinContent(ix_src))
        h_new.SetBinError(ix_dest, np.hypot(h_new.GetBinError(ix_dest), h.GetBinError(ix_src)))
    return h_new


def merge_th2_bin_pairs(h, bin_pairs_x, bin_pairs_y, new_bin_edges_x, new_bin_edges_y):
    if not bin_pairs_x or len(bin_pairs_x) == 0:
        print("Skipping merging x pairs")
        h_new = h
    else:
        # do x bin merging first
        bin_edges_y_orig = cu.get_bin_edges(h, 'y')
        h_new = ROOT.TH2D(h.GetName() + cu.get_unique_str(),
                          ";".join([h.GetTitle(), h.GetXaxis().GetTitle(), h.GetYaxis().GetTitle()]),
                          len(new_bin_edges_x)-1,
                          array('d', new_bin_edges_x),
                          len(bin_edges_y_orig)-1,
                          array('d', bin_edges_y_orig))

        for iy in range(0, h_new.GetNbinsY()+2):
            for ix in range(0, h_new.GetNbinsX()+2):
                # copy original data first, but only the subset of overlapping bins
                h_new.SetBinContent(ix, iy, h.GetBinContent(ix, iy))
                h_new.SetBinError(ix, iy, h.GetBinError(ix, iy))

                # now do the extra bins that don't appear in h_new,
                # updating the bins in the new hist
                for ix_src, ix_dest in bin_pairs_x:
                    h_new.SetBinContent(ix_dest, iy, h_new.GetBinContent(ix_dest, iy) + h.GetBinContent(ix_src, iy))
                    h_new.SetBinError(ix_dest, iy, np.hypot(h_new.GetBinError(ix_dest, iy), h.GetBinError(ix_src, iy)))

    if not bin_pairs_y or len(bin_pairs_y) == 0:
        print("Skipping merging y pairs")
        return h_new

    # now do y bin merging
    bin_edges_x_orig = cu.get_bin_edges(h_new, 'x')
    h_new2 = ROOT.TH2D(h.GetName() + cu.get_unique_str(),
                       ";".join([h.GetTitle(), h.GetXaxis().GetTitle(), h.GetYaxis().GetTitle()]),
                       len(bin_edges_x_orig)-1,
                       array('d', bin_edges_x_orig),
                       len(new_bin_edges_y)-1,
                       array('d', new_bin_edges_y))

    for ix in range(0, h_new2.GetNbinsX()+2):
        for iy in range(0, h_new2.GetNbinsY()+2):
            # copy original data first, but only the subset of overlapping bins
            h_new2.SetBinContent(ix, iy, h_new.GetBinContent(ix, iy))
            h_new2.SetBinError(ix, iy, h_new.GetBinError(ix, iy))

        # now do the extra bins that don't appear in h_new,
        # updating the bins in the new hist
        for iy_src, iy_dest in bin_pairs_y:
            h_new2.SetBinContent(ix, iy_dest, h_new2.GetBinContent(ix, iy_dest) + h_new.GetBinContent(ix, iy_src))
            h_new2.SetBinError(ix, iy_dest, np.hypot(h_new2.GetBinError(ix, iy_dest), h_new.GetBinError(ix, iy_src)))

    return h_new2


def setup_merged_last_pt_binning(orig_binning_handler):
    """Setup BinningHandler for merged last pt bins, along with bin pairs to merge.

    These can be used in merge_th1_bin_pairs(), etc

    Parameters
    ----------
    orig_binning_handler : BinningHandler
        To get original binning scheme.
        Must have PtVarBinning objects for generator_/detector_binning,
        otherwise it won't work.

    Returns
    -------
    BinningHandler, list[[int, int]], list[[int, int]]
        New BinningHandler, list of gen bin pairs to merge, list of reco bin pairs to merge

    Raises
    ------
    TypeError
        If BinningHandler doesn't have PtVarBinning objects
    ValueError
        If mismatch in variable bins in different pt bins

    """
    # construct bin pairs to merge
    # generator
    orig_gen_binning = orig_binning_handler.get_binning_scheme('generator')
    src_pt = orig_gen_binning.get_pt_bins(is_signal_region=True)[-1]
    var_bins_src = orig_gen_binning.get_variable_bins(src_pt)
    dest_pt = orig_gen_binning.get_pt_bins(is_signal_region=True)[-2]
    var_bins_dest = orig_gen_binning.get_variable_bins(dest_pt)

    if set(var_bins_dest) != set(var_bins_src):
        raise ValueError("Cannot merge last pt bin - different variable bins")

    gen_bin_pairs = []

    these_var_bins = var_bins_src if orig_gen_binning.var_of else var_bins_src[:-1]
    # do all to include overflow bin
    for var in these_var_bins:
        src_bin = orig_gen_binning.physical_bin_to_global_bin(pt=src_pt, var=var + 1E-6)
        dest_bin = orig_gen_binning.physical_bin_to_global_bin(pt=dest_pt, var=var + 1E-6)
        gen_bin_pairs.append([src_bin, dest_bin])

    # create new binning objects
    if not isinstance(orig_gen_binning, PtVarBinning):
        raise TypeError("Expected orig_gen_binning to be type PtVarBinning - cannot setup otherwise")
    generator_binning_merged = PtVarBinning(variable_bin_edges=orig_gen_binning.variable_bin_edges,
                                            variable_name=orig_gen_binning.variable_name,
                                            pt_bin_edges_signal=orig_gen_binning.pt_bin_edges_signal,
                                            pt_bin_edges_underflow=orig_gen_binning.pt_bin_edges_underflow,
                                            binning_name=orig_gen_binning.binning_name,
                                            binning_underflow_name=orig_gen_binning.binning_underflow_name,
                                            binning_signal_name=orig_gen_binning.binning_signal_name,
                                            var_uf=orig_gen_binning.var_uf,
                                            var_of=orig_gen_binning.var_of,
                                            pt_uf=orig_gen_binning.pt_uf,
                                            pt_of=False)

    # detector
    orig_reco_binning = orig_binning_handler.get_binning_scheme('detector')
    reco_bin_pairs = []

    # do twice: 1500-inf -> 408-954, then 954-1500 -> 408-954
    src_pt = orig_reco_binning.get_pt_bins(is_signal_region=True)[-1]
    var_bins_src = orig_reco_binning.get_variable_bins(src_pt)
    src2_pt = orig_reco_binning.get_pt_bins(is_signal_region=True)[-2]
    dest_pt = orig_reco_binning.get_pt_bins(is_signal_region=True)[-3]
    var_bins_dest = orig_reco_binning.get_variable_bins(dest_pt)

    if set(var_bins_dest) != set(var_bins_src):
        raise ValueError("Cannot merge last pt bin - different variable bins")
    # do all to include overflow bin
    these_var_bins = var_bins_src if orig_reco_binning.var_of else var_bins_src[:-1]

    for var in these_var_bins:
        src_bin = orig_reco_binning.physical_bin_to_global_bin(pt=src_pt, var=var + 1E-6)
        src2_bin = orig_reco_binning.physical_bin_to_global_bin(pt=src2_pt, var=var + 1E-6)
        dest_bin = orig_reco_binning.physical_bin_to_global_bin(pt=dest_pt, var=var + 1E-6)
        reco_bin_pairs.append([src_bin, dest_bin])
        reco_bin_pairs.append([src2_bin, dest_bin])

    # create new binning objects
    if not isinstance(orig_reco_binning, PtVarBinning):
        raise TypeError("Expected orig_reco_binning to be type PtVarBinning - cannot setup otherwise")
    detector_binning_merged = PtVarBinning(variable_bin_edges=orig_reco_binning.variable_bin_edges,
                                           variable_name=orig_reco_binning.variable_name,
                                           pt_bin_edges_signal=np.concatenate((orig_reco_binning.pt_bin_edges_signal[:-2], orig_reco_binning.pt_bin_edges_signal[-1:])),
                                           pt_bin_edges_underflow=orig_reco_binning.pt_bin_edges_underflow,
                                           binning_name=orig_reco_binning.binning_name,
                                           binning_underflow_name=orig_reco_binning.binning_underflow_name,
                                           binning_signal_name=orig_reco_binning.binning_signal_name,
                                           var_uf=orig_reco_binning.var_uf,
                                           var_of=orig_reco_binning.var_of,
                                           pt_uf=orig_reco_binning.pt_uf,
                                           pt_of=False)

    binning_handler_merged = BinningHandler(generator_ptvar_binning=generator_binning_merged,
                                            detector_ptvar_binning=detector_binning_merged)

    return binning_handler_merged, gen_bin_pairs, reco_bin_pairs


def get_bins_to_merge_probability_stats(pmatrix, min_num=10, allow_zero_entries=True, use_neff_entries=False, max_p=0.7, max_rel_err=0.6):
    """Figure out which gen bins to merge based on poor number of events & error bars
    in given column of probability matrix

    Parameters
    ----------
    pmatrix : ROOT.TH2
        Description
    min_num : int, optional
        Minimum number of non-zero entries in each column
    allow_zero_entries : bool, optional
        If True allow matrix columns with 0 contents, otherwise include them to be merged
    use_neff_entries : bool, optional
        If True use number of effective entries with min_num,
        otherwise use number of entries
    max_p : float, optional
        Maximum probability in one bin in each column
    max_rel_err : float, optional
        Maximum relative uncertainty across all bins in each column

    Returns
    -------
    list[int]
        Bin numbers
    """
    print("Checking probability matrix stats...")
    bad_bins = []
    for ix in range(1, pmatrix.GetNbinsX()+1):
        value_err_pairs = []
        # Get all the non-zero entries in the prob matrix for this column
        for iy in range(1, pmatrix.GetNbinsY()+1):
            val = pmatrix.GetBinContent(ix, iy)
            if val == 0:
                continue
            err = pmatrix.GetBinError(ix, iy)
            value_err_pairs.append([val, abs(err/val)])

        n_eff = 0
        if len(value_err_pairs) > 0:
            # effective number of entries - useful when dealing with weighted elements
            # For a weighted histogram, this number corresponds to the hypothetical
            # number of unweighted entries a histogram would need to have the same statistical power.
            n_eff = sum([v for v, _ in value_err_pairs]) / sum([v**2 for v, _ in value_err_pairs])
            print("p matrix gen bin", ix, "len value_err_pairs:", len(value_err_pairs), "n_eff:", n_eff)

        def _print_val_err_pairs(pairs):
            return ",".join(["{} ± {}".format(v, e) for v, e in pairs])

        count_entries = n_eff if use_neff_entries else len(value_err_pairs)
        if not allow_zero_entries and count_entries == 0:
            # handle case of 0 entries separately
            print("p matrix gen bin", ix, "has 0 values")
            bad_bins.append(ix)
        elif 0 < count_entries < min_num:
            print("p matrix gen bin", ix, "has <", min_num, "values:", _print_val_err_pairs(value_err_pairs))
            bad_bins.append(ix)
        elif any([v > max_p for v, e in value_err_pairs]):
            print("p matrix gen bin", ix, "has bin with value >", max_p, "values:", __print_val_err_pairs(value_err_pairs))
            bad_bins.append(ix)
        elif any([e > max_rel_err for v, e in value_err_pairs]):
            print("p matrix gen bin", ix, "has bin with fractional error >", max_rel_err, "values:", _print_val_err_pairs(value_err_pairs))
            bad_bins.append(ix)
    return bad_bins


def get_bins_to_merge_nonnegative(h, binning_obj, max_bin=-1, ensure_nonnegative=True):
    """Figure out which bins to merge to ensure bin contents >= 0

    Parameters
    ----------
    h : TH1D
        Histogram
    binning_obj : PtVarBinning, PtVarPerPtBinning
        Object to get binning info from, to ensure we don't cross pt bin boundary
    max_bin : int, optional
        Only consider bins less than max_bin if > 0.
    ensure_nonnegative : bool, optional
        If True, figues out which set of adjacent bins must be merged to get
        a result > 0. Otherwise returns the individual -ve bins

    Returns
    -------
    list[int]
        List of bin numbers that should be merged with their lefthand neighbour
        e.g. with merge_th1_bins()
    """
    bins_to_merge = []
    if max_bin < 0:
        max_bin = h.GetNbinsX()
    contents = np.array([h.GetBinContent(i) for i in range(1, max_bin+1)])

    def _is_boundary_bin(bin_number):
        pt_curr = binning_obj.global_bin_to_physical_bin(bin_number).pt[0]
        pt_nxt = binning_obj.global_bin_to_physical_bin(bin_number-1).pt[0]
        return pt_nxt != pt_curr

    for bin_num, x in enumerate(contents, 1):
        this_merge_set = []
        if x < 0:
            negative_bin_value = x
            negative_bin_error = h.GetBinError(bin_num)
            print("Found -ve bin", bin_num, binning_obj.global_bin_to_physical_bin(bin_num), "=", negative_bin_value, "±", negative_bin_error)

            if _is_boundary_bin(bin_num):
                pt_current = binning_obj.global_bin_to_physical_bin(bin_num).pt[0]
                pt_next = binning_obj.global_bin_to_physical_bin(bin_num-1).pt[0]
                warnings.warn("Cannot merge bin {0} across pT boundary {1} -> {2}, skipping".format(bin_num, pt_current, pt_next))
                continue

            this_merge_set.append(bin_num)

            if ensure_nonnegative:
                # figure out how many adjacent bins you need to sum to make it +ve
                summed_bins = negative_bin_value

                # iterate over lefthand bins, but check you don't cross a pt boundary
                iter_ind = 1
                while summed_bins < 0:
                    this_bin = bin_num - iter_ind

                    if _is_boundary_bin(this_bin):
                        pt_current = binning_obj.global_bin_to_physical_bin(this_bin).pt[0]
                        pt_next = binning_obj.global_bin_to_physical_bin(this_bin-1).pt[0]
                        warnings.warn("Cannot merge bin {0} across pT boundary {1} -> {2}, not including".format(this_bin, pt_current, pt_next))
                        break

                    summed_bins += contents[this_bin-1]  # -1 as numpy array
                    print("..new sum", summed_bins)
                    if summed_bins < 0:
                        this_merge_set.append(this_bin)

                    iter_ind += 1

                this_merge_set = this_merge_set[::-1]  # reverse order so ascending
                print("Found set that sums > 0:", this_merge_set, "=", summed_bins)

            bins_to_merge.extend(this_merge_set)
    return sorted(list(set(bins_to_merge)))


def setup_merged_bin_binning(gen_bins_to_merge, reco_bins_to_merge, orig_binning_handler):
    """Setup binning scheme with arbitray gen & reco bins merged.

    Returns BinningHandler with PtVarPerPtBinning objects
    """
    gen_bins_to_merge = gen_bins_to_merge or []
    print("gen bins to merge:", gen_bins_to_merge)
    for b in gen_bins_to_merge:
        print(b, ":", orig_binning_handler.global_bin_to_physical_bin(b, "generator"))
    reco_bins_to_merge = reco_bins_to_merge or []
    print("reco bins to merge:", reco_bins_to_merge)
    for b in reco_bins_to_merge:
        print(b, ":", orig_binning_handler.global_bin_to_physical_bin(b, "detector"))

    # create configs for PerPtBinning
    # remove bins that need removing
    gen_underflow_config = []
    gen_pt_bins_uflow = orig_binning_handler.get_pt_bins(binning_scheme='generator', is_signal_region=False)
    for pt_low, pt_high in zip(gen_pt_bins_uflow[:-1], gen_pt_bins_uflow[1:]):
        var_bins = list(orig_binning_handler.get_variable_bins(pt_low, "generator"))
        for b in gen_bins_to_merge:
            phys_bin = orig_binning_handler.global_bin_to_physical_bin(b, "generator")
            if phys_bin.pt[0] != pt_low:
                continue
            var = phys_bin.var[0]
            if var in var_bins:
                var_bins.remove(var)
        if len(var_bins) > 0:
            gen_underflow_config.append([(pt_low, pt_high), var_bins])

    gen_signal_config = []
    gen_pt_bins_signal = orig_binning_handler.get_pt_bins(binning_scheme='generator', is_signal_region=True)
    for pt_low, pt_high in zip(gen_pt_bins_signal[:-1], gen_pt_bins_signal[1:]):
        var_bins = list(orig_binning_handler.get_variable_bins(pt_low, "generator"))
        for b in gen_bins_to_merge:
            phys_bin = orig_binning_handler.global_bin_to_physical_bin(b, "generator")
            if phys_bin.pt[0] != pt_low:
                continue
            var = phys_bin.var[0]
            if var in var_bins:
                var_bins.remove(var)
        if len(var_bins) > 0:
            gen_signal_config.append([(pt_low, pt_high), var_bins])

    # Add pt overflow bin explicitly
    has_pt_of = orig_binning_handler.get_binning_scheme('generator').pt_of
    if has_pt_of:
        pt_low, pt_high = gen_pt_bins_signal[-1], 13000
        var_bins = list(orig_binning_handler.get_variable_bins(pt_low, "generator"))
        for b in gen_bins_to_merge:
            phys_bin = orig_binning_handler.global_bin_to_physical_bin(b, "generator")
            if phys_bin.pt[0] != pt_low:
                continue
            var = phys_bin.var[0]
            if var in var_bins:
                var_bins.remove(var)
        if len(var_bins) > 0:
            gen_signal_config.append([(pt_low, pt_high), var_bins])

    reco_underflow_config = []
    reco_pt_bins_uflow = orig_binning_handler.get_pt_bins(binning_scheme='detector', is_signal_region=False)
    for pt_low, pt_high in zip(reco_pt_bins_uflow[:-1], reco_pt_bins_uflow[1:]):
        var_bins = list(orig_binning_handler.get_variable_bins(pt_low, "detector"))
        for b in reco_bins_to_merge:
            phys_bin = orig_binning_handler.global_bin_to_physical_bin(b, "detector")
            if phys_bin.pt[0] != pt_low:
                continue
            var = phys_bin.var[0]
            if var in var_bins:
                var_bins.remove(var)
        if len(var_bins) > 0:
            reco_underflow_config.append([(pt_low, pt_high), var_bins])

    reco_signal_config = []
    reco_pt_bins_signal = orig_binning_handler.get_pt_bins(binning_scheme='detector', is_signal_region=True)
    for pt_low, pt_high in zip(reco_pt_bins_signal[:-1], reco_pt_bins_signal[1:]):
        var_bins = list(orig_binning_handler.get_variable_bins(pt_low, "detector"))
        for b in reco_bins_to_merge:
            phys_bin = orig_binning_handler.global_bin_to_physical_bin(b, "detector")
            if phys_bin.pt[0] != pt_low:
                continue
            var = phys_bin.var[0]
            if var in var_bins:
                var_bins.remove(var)
        if len(var_bins) > 0:
            reco_signal_config.append([(pt_low, pt_high), var_bins])

    # Add pt overflow bin explicitly
    has_pt_of = orig_binning_handler.get_binning_scheme('detector').pt_of
    if has_pt_of:
        pt_low, pt_high = reco_pt_bins_signal[-1], 13000
        var_bins = list(orig_binning_handler.get_variable_bins(pt_low, "detector"))
        for b in reco_bins_to_merge:
            phys_bin = orig_binning_handler.global_bin_to_physical_bin(b, "detector")
            if phys_bin.pt[0] != pt_low:
                continue
            var = phys_bin.var[0]
            if var in var_bins:
                var_bins.remove(var)
        if len(var_bins) > 0:
            reco_signal_config.append([(pt_low, pt_high), var_bins])

    print('gen_underflow_config:'); pprint(gen_underflow_config, width=200)
    print('gen_signal_config:'); pprint(gen_signal_config, width=200)
    print('reco_underflow_config:'); pprint(reco_underflow_config, width=200)
    print('reco_signal_config:'); pprint(reco_signal_config, width=200)

    gen_perpt_binning = PtVarPerPtBinning(orig_binning_handler.variable_name,
                                          gen_underflow_config,
                                          gen_signal_config,
                                          orig_binning_handler.generator_ptvar_binning.binning_name,
                                          orig_binning_handler.generator_ptvar_binning.binning_underflow_name,
                                          orig_binning_handler.generator_ptvar_binning.binning_signal_name,
                                          var_uf=orig_binning_handler.generator_ptvar_binning.var_uf,
                                          var_of=orig_binning_handler.generator_ptvar_binning.var_of)

    reco_perpt_binning = PtVarPerPtBinning(orig_binning_handler.variable_name,
                                           reco_underflow_config,
                                           reco_signal_config,
                                           orig_binning_handler.detector_ptvar_binning.binning_name,
                                           orig_binning_handler.detector_ptvar_binning.binning_underflow_name,
                                           orig_binning_handler.detector_ptvar_binning.binning_signal_name,
                                           var_uf=orig_binning_handler.detector_ptvar_binning.var_uf,
                                           var_of=orig_binning_handler.detector_ptvar_binning.var_of)

    binning_handler_merged = BinningHandler(generator_ptvar_binning=gen_perpt_binning,
                                            detector_ptvar_binning=reco_perpt_binning)

    return binning_handler_merged


MERGE_JSON_FILENAME = "merged_bins_history.json"


def save_merge_bin_history_to_file(merge_bin_history, output_dir):
    cu.check_dir_exists_create(output_dir)
    with open(os.path.join(output_dir, MERGE_JSON_FILENAME), 'w') as f:
        json.dump(merge_bin_history, f, indent=4)


def get_merge_bin_history_from_file(merge_dir):
    with open(os.path.join(merge_dir, MERGE_JSON_FILENAME)) as f:
        merge_bin_history = json.load(f)
    return merge_bin_history


# To be able to export to XML, since getting std::ostream from python is impossible?
my_binning_xml_code = """
class BinningXMLExporter {
public:
    BinningXMLExporter() {cout << " Creating BinningXMLExporter " << endl;}

    static Int_t ExportXML(const TUnfoldBinning &binning,
                           std::string dirname,
                           std::string filename,
                           Bool_t writeHeader,
                           Bool_t writeFooter,
                           Int_t indent=0) {

        ofstream dtdFile(TString::Format("%s/tunfoldbinning.dtd", dirname.c_str()));
        TUnfoldBinningXML::WriteDTD(dtdFile);
        dtdFile.close();

        ofstream xmlfile;
        xmlfile.open(TString::Format("%s/%s", dirname.c_str(), filename.c_str()), ios::out);
        Int_t result = TUnfoldBinningXML::ExportXML(binning, xmlfile, writeHeader, writeFooter, indent);
        xmlfile.close();
        return result;
    }
};
"""
ROOT.gInterpreter.ProcessLine(my_binning_xml_code)


def setup_regions_from_args(args):
    """Create region dict(s) from args"""
    regions = []
    src_dir = args.source

    if any([args.doDijetCentral, args.doDijetForward, args.doDijetCentralGroomed, args.doDijetForwardGroomed]):
        tau_limits = {
            'jet_puppiMultiplicity': (1E-3, 1E2),
            'jet_pTD': (1E-1, 1E3),
            'jet_LHA': (1E-1, 1E3),
            'jet_width': (1E-1, 1E3),
            'jet_thrust': (1E-1, 1E3),
            'jet_puppiMultiplicity_charged': (1E-1, 1E3),
            'jet_pTD_charged': (1E-1, 1E3),
            'jet_LHA_charged': (1E-1, 1E3),
            'jet_width_charged': (1E-1, 1E3),
            'jet_thrust_charged': (1E-1, 1E3),
        }

        if args.doDijetCentral:
            dijet_region_central_dict = get_dijet_config(src_dir, central=True, groomed=False)
            dijet_region_central_dict['tau_limits'] = tau_limits
            regions.append(dijet_region_central_dict)

        if args.doDijetForward:
            dijet_region_forward_dict = get_dijet_config(src_dir, central=False, groomed=False)
            dijet_region_forward_dict['tau_limits'] = tau_limits
            regions.append(dijet_region_forward_dict)

        if args.doDijetCentralGroomed:
            dijet_region_central_groomed_dict = get_dijet_config(src_dir, central=True, groomed=True)
            dijet_region_central_groomed_dict['tau_limits'] = tau_limits
            regions.append(dijet_region_central_groomed_dict)

        if args.doDijetForwardGroomed:
            dijet_region_forward_groomed_dict = get_dijet_config(src_dir, central=False, groomed=True)
            dijet_region_forward_groomed_dict['tau_limits'] = tau_limits
            regions.append(dijet_region_forward_groomed_dict)

    if any([args.doZPJ, args.doZPJGroomed]):
        tau_limits = {
            'jet_puppiMultiplicity': (1E-1, 1E3),
            'jet_pTD': (1E-1, 1E3),
            'jet_LHA': (1E-1, 1E3),
            'jet_width': (1E-1, 1E3),
            'jet_thrust': (1E-1, 1E3),
            'jet_puppiMultiplicity_charged': (1E-1, 1E3),
            'jet_pTD_charged': (1E-1, 1E3),
            'jet_LHA_charged': (1E-1, 1E3),
            'jet_width_charged': (1E-1, 1E3),
            'jet_thrust_charged': (1E-1, 1E3),
        }

        if args.doZPJ:
            zpj_region_dict = get_zpj_config(src_dir, groomed=False)
            zpj_region_dict['tau_limits'] = tau_limits
            regions.append(zpj_region_dict)

        if args.doZPJGroomed:
            zpj_region_groomed_dict = get_zpj_config(src_dir, groomed=True)
            zpj_region_groomed_dict['tau_limits'] = tau_limits
            regions.append(zpj_region_groomed_dict)

    return regions


def modify_region_using_args(region, args):
    """Modify region according to input args"""
    if not args.doJackknifeInput:
        region['jackknife_input_variations'] = []

    if not args.doJackknifeResponse:
        region['jackknife_response_variations'] = []

    if args.doExperimentalSysts or args.doExperimentalSystsFromFile:
        # Remove the lumi one if we have no backgrounds, or the user has not said to remove backgrounds
        region['experimental_systematics'] = [syst_dict for syst_dict in region['experimental_systematics']
                                              if not ('lumi' in syst_dict['label'].lower()
                                                      and (len(region.get('backgrounds', [])) == 0
                                                           or not args.subtractBackgrounds))]
        if args.doExperimentalSystsOnlyHerwig:
            # only herwig related systs
            region['experimental_systematics'] = [s for s in region['experimental_systematics']
                                                  if 'herwig' in s['label'].lower() or 'shower' in s['label'].lower()]

    else:
        region['experimental_systematics'] = []

    if not (args.doScaleSysts or args.doScaleSystsFromFile):
        region['scale_systematics'] = []

    if not (args.doModelSysts or args.doModelSystsFromFile):
        region['model_systematics'] = []

    elif args.doModelSystsOnlyHerwig:
        # only herwig related systs
        region['model_systematics'] = [s for s in region['model_systematics']
                                       if 'herwig' in s['label'].lower()]
    elif args.doModelSystsOnlyScale:
        # only scale related systs
        region['model_systematics'] = [s for s in region['model_systematics']
                                       if 'mu' in s['label'].lower()]

    elif args.doModelSystsNotScale:
        # only non-scale related systs
        region['model_systematics'] = [s for s in region['model_systematics']
                                       if 'mu' not in s['label'].lower()]

    if not (args.doPDFSysts or args.doPDFSystsFromFile):
        region['pdf_systematics'] = []


@profile
def main():
    """Main function to run setup, unfolding, etc.

    In a function, because it's handy for profiling
    """
    parser = get_unfolding_argparser(description=__doc__)
    args = parser.parse_args()
    print("")
    print(args)
    print("")
    sanitise_args(args)
    output_dir = get_unfolding_output_dir(args)
    cu.check_dir_exists_create(output_dir)
    print("")
    print(args)

    print("")
    print("Outputting to", output_dir)
    print("")

    src_dir = args.source

    # Setup regions to unfold
    # --------------------------------------------------------------------------
    regions = setup_regions_from_args(args)

    # Setup various options
    # --------------------------------------------------------------------------
    REGULARIZE = args.regularize
    MC_INPUT = args.MCinput
    MC_SPLIT = args.MCsplit if args.MCinput else False

    jet_algo = "AK4"
    if "ak8puppi" in src_dir:
        jet_algo = "AK8"

    if args.angles[0] == "all":
        angles = qgc.COMMON_VARS
    else:
        angles = [a for a in qgc.COMMON_VARS if a.var in args.angles]

    print("Running TUnfold version", ROOT.TUnfold.GetTUnfoldVersion())

    LAMBDA_VAR_DICTS = qgc.VAR_UNFOLD_DICT

    # Do unfolding per signal region
    # --------------------------------------------------------------------------
    for orig_region in regions[:]:

        # Setup pt bins
        # -------------
        # need different ones for Z+Jets region
        is_zpj = "ZPlusJets" in orig_region['name']

        zpj_append = "_zpj" if is_zpj else ""

        pt_bin_edges_gen = qgc.PT_UNFOLD_DICT['signal%s_gen' % zpj_append]
        pt_bin_edges_reco = qgc.PT_UNFOLD_DICT['signal%s_reco' % zpj_append]

        pt_bin_edges_underflow_gen = qgc.PT_UNFOLD_DICT['underflow%s_gen' % zpj_append]
        pt_bin_edges_underflow_reco = qgc.PT_UNFOLD_DICT['underflow%s_reco' % zpj_append]

        # Modify systematics & other keys as necessary
        # ----------------------------------------------------------------------
        modify_region_using_args(orig_region, args)

        # Do unfolding for each angle
        # ----------------------------------------------------------------------
        for angle in angles:
            region = copy(orig_region)  # make copy since we might modify it later, e.g. PDF, and want same start for each angle
            is_groomed = "groomed" in region['name']

            if isinstance(region['mc_tfile'], str):
                region['mc_tfile'] = cu.open_root_file(region['mc_tfile'])

            angle_prepend = "groomed " if is_groomed else ""
            append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name
            angle_str = "%s%s (%s)" % (angle_prepend, qgc.lower_angle_name(angle), angle.lambda_str)

            print("*"*120)
            print("Region/var: %s %s" % (region['name'], angle.var))
            print("*"*120)

            # put plots in subdir, to avoid overcrowding
            this_output_dir = "%s/%s/%s" % (output_dir, region['name'], angle.var)
            cu.check_dir_exists_create(this_output_dir)

            # Save minimal set of hists etc to ROOT file for access later
            # For everything, we pickle it
            output_tfile_slim = ROOT.TFile("%s/unfolding_result_slim.root" % this_output_dir, "RECREATE")
            new_tdir = "%s/%s" % (region['name'], angle.var)
            output_tfile_slim.mkdir(new_tdir)
            this_slim_tdir = output_tfile_slim.Get(new_tdir)
            this_slim_tdir.cd()

            # Setup MyUnfolder object to do unfolding etc
            # -------------------------------------------
            groom_str = "groomed" if is_groomed else "ungroomed"
            angle_bin_edges_reco = LAMBDA_VAR_DICTS[groom_str][angle.var]['reco']
            angle_bin_edges_gen = LAMBDA_VAR_DICTS[groom_str][angle.var]['gen']
            angle_shortname = angle.var.replace("jet_", "")

            # Setup binning scheme
            # ------------------------------------------------------------------
            print('reco lambda binning', angle_bin_edges_reco)
            print('gen lambda binning', angle_bin_edges_gen)
            print('reco pt binning', pt_bin_edges_reco)
            print('gen pt binning', pt_bin_edges_gen)
            print('reco uflow pt binning', pt_bin_edges_underflow_reco)
            print('gen uflow pt binning', pt_bin_edges_underflow_gen)

            variable_name = "%s%s" % (angle_prepend, angle.name)
            var_uf = False
            var_of = True
            pt_uf = False
            pt_of = True
            generator_binning = PtVarBinning(variable_bin_edges=angle_bin_edges_gen,
                                             variable_name=variable_name,
                                             pt_bin_edges_signal=pt_bin_edges_gen,
                                             pt_bin_edges_underflow=pt_bin_edges_underflow_gen,
                                             binning_name="generator",
                                             binning_underflow_name="generatordistribution_underflow",
                                             binning_signal_name="generatordistribution",
                                             var_uf=var_uf,
                                             var_of=var_of,
                                             pt_uf=pt_uf,
                                             pt_of=pt_of)

            detector_binning = PtVarBinning(variable_bin_edges=angle_bin_edges_reco,
                                            variable_name=variable_name,
                                            pt_bin_edges_signal=pt_bin_edges_reco,
                                            pt_bin_edges_underflow=pt_bin_edges_underflow_reco,
                                            binning_name="detector",
                                            binning_underflow_name="detectordistribution_underflow",
                                            binning_signal_name="detectordistribution",
                                            var_uf=var_uf,
                                            var_of=var_of,
                                            pt_uf=pt_uf,
                                            pt_of=pt_of)

            binning_handler = BinningHandler(generator_ptvar_binning=generator_binning,
                                             detector_ptvar_binning=detector_binning)

            def zero_column_th2(h2d, column_num):
                for iy in range(0, h2d.GetNbinsY()+2):
                    if h2d.GetBinContent(column_num, iy) != 0:
                        print("Zeroing (%d, %d)" % (column_num, iy))
                    h2d.SetBinContent(column_num, iy, 0)
                    h2d.SetBinError(column_num, iy, 0)

            def zero_row_th2(h2d, row_num):
                for ix in range(0, h2d.GetNbinsX()+2):
                    h2d.SetBinContent(ix, row_num, 0)
                    h2d.SetBinError(ix, row_num, 0)

            # Remove last pt bin entries - detector level
            def zero_last_pt_bin_th1_reco(h):
                start_bin = binning_handler.physical_bin_to_global_bin(pt=pt_bin_edges_reco[-1]+1, var=0, binning_scheme='detector')
                end_bin = binning_handler.physical_bin_to_global_bin(pt=pt_bin_edges_reco[-1]+1, var=9999, binning_scheme='detector')
                for rb in range(start_bin, end_bin + 1):
                    h.SetBinContent(rb, 0)
                    h.SetBinError(rb, 0)

            def zero_last_pt_bin_th2_reco(h):
                start_bin = binning_handler.physical_bin_to_global_bin(pt=pt_bin_edges_reco[-1]+1, var=0, binning_scheme='detector')
                end_bin = binning_handler.physical_bin_to_global_bin(pt=pt_bin_edges_reco[-1]+1, var=9999, binning_scheme='detector')
                for ix in range(0, h.GetNbinsX()+2):
                    for rb in range(start_bin, end_bin+1):
                        h.SetBinContent(ix, rb, 0)
                        h.SetBinError(ix, rb, 0)

            def zero_last_pt_bin_th2_gen(h):
                start_bin = binning_handler.physical_bin_to_global_bin(pt=pt_bin_edges_gen[-1]+1, var=0, binning_scheme='generator')
                end_bin = binning_handler.physical_bin_to_global_bin(pt=pt_bin_edges_gen[-1]+1, var=9999, binning_scheme='generator')
                for gb in range(start_bin, end_bin + 1):
                    for iy in range(0, h.GetNbinsY()+2):
                        h.SetBinContent(gb, iy, 0)
                        h.SetBinError(gb, iy, 0)

            # Get lots of hists from input files
            # ------------------------------------------------------------------
            mc_hname_append = "split" if MC_SPLIT else "all"
            hist_mc_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_%s" % (region['dirname'], angle_shortname, mc_hname_append))
            hist_mc_gen = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_truth_%s" % (region['dirname'], angle_shortname, mc_hname_append))

            hist_mc_reco = rm_large_rel_error_bins_th1(hist_mc_reco, relative_err_threshold=args.relErr)
            hist_mc_gen = rm_large_rel_error_bins_th1(hist_mc_gen, relative_err_threshold=args.relErr)

            # _all versions used in TruthTemplateMaker and fakes, since we always want the full stats
            hist_mc_reco_all = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_all" % (region['dirname'], angle_shortname))
            hist_mc_gen_all = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_truth_all" % (region['dirname'], angle_shortname))
            hist_mc_gen_reco_map = cu.get_from_tfile(region['mc_tfile'], "%s/tu_%s_GenReco_%s" % (region['dirname'], angle_shortname, mc_hname_append))

            hist_mc_gen_all = rm_large_rel_error_bins_th1(hist_mc_gen_all, relative_err_threshold=args.relErr)
            hist_mc_reco_all = rm_large_rel_error_bins_th1(hist_mc_reco_all, relative_err_threshold=args.relErr)

            hist_data_reco = None
            if not MC_INPUT:
                if not isinstance(region['data_tfile'], ROOT.TFile):
                    region['data_tfile'] = cu.open_root_file(region['data_tfile'])
                hist_data_reco = cu.get_from_tfile(region['data_tfile'], "%s/hist_%s_reco_all" % (region['dirname'], angle_shortname))

            # Actual distribution to be unfolded
            reco_1d = hist_mc_reco.Clone() if MC_INPUT else hist_data_reco
            reco_1d = rm_large_rel_error_bins_th1(reco_1d, relative_err_threshold=args.relErr)

            # to construct our "fakes" template, we use the ratio as predicted by MC using all stats,
            # and apply it to data. this way we ensure we don't have -ve values,
            # and avoid any issue with cross sections
            hist_mc_fakes_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_fake_all" % (region['dirname'], angle_shortname))
            hist_mc_fakes_reco = rm_large_rel_error_bins_th1(hist_mc_fakes_reco, relative_err_threshold=args.relErr)

            hist_fake_fraction = hist_mc_fakes_reco.Clone("hist_%s_fakes_reco_fraction" % angle_shortname)
            hist_fake_fraction.Divide(rm_large_rel_error_bins_th1(hist_mc_reco_all, relative_err_threshold=args.relErr))

            if args.zeroLastPtBin:
                zero_last_pt_bin_th1_reco(hist_mc_reco)
                zero_last_pt_bin_th1_reco(hist_mc_reco_all)
                # zero_last_pt_bin_th1_gen(hist_mc_gen)
                # zero_last_pt_bin_th1_gen(hist_mc_gen_all)
                zero_last_pt_bin_th1_reco(reco_1d)
                zero_last_pt_bin_th2_reco(hist_mc_gen_reco_map)
                zero_last_pt_bin_th2_gen(hist_mc_gen_reco_map)
                zero_last_pt_bin_th1_reco(hist_mc_fakes_reco)
                zero_last_pt_bin_th1_reco(hist_fake_fraction)

            input_handler_args = dict(
                input_hist=reco_1d,
                hist_truth=hist_mc_gen,
                hist_mc_reco=hist_mc_reco,
                hist_mc_fakes=hist_mc_fakes_reco,
                hist_fake_fraction=hist_fake_fraction
            )

            # Gen binning versions
            # --------------------
            mc_hname_append = "_split" if MC_SPLIT else ""  # FIXME consistency in unfold hist module!
            hist_data_reco_gen_binning = None
            if not MC_INPUT:
                hist_data_reco_gen_binning = cu.get_from_tfile(region['data_tfile'], "%s/hist_%s_reco_gen_binning" % (region['dirname'], angle_shortname))
            hist_mc_reco_gen_binning = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_gen_binning%s" % (region['dirname'], angle_shortname, mc_hname_append))
            hist_mc_reco_gen_binning_all = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_gen_binning" % (region['dirname'], angle_shortname))

            # Actual distribution to be unfolded, but with gen binning
            reco_1d_gen_binning = hist_mc_reco_gen_binning.Clone() if MC_INPUT else hist_data_reco_gen_binning

            # hist_fake_fraction_gen_binning = None
            # create template as above, but with gen binning
            hist_mc_fakes_reco_gen_binning = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_fake_gen_binning" % (region['dirname'], angle_shortname))
            hist_fake_fraction_gen_binning = hist_mc_fakes_reco_gen_binning.Clone("hist_%s_fakes_fraction_gen_binning" % angle_shortname)
            hist_fake_fraction_gen_binning.Divide(hist_mc_reco_gen_binning_all)

            input_handler_gen_binning_args = dict(
                input_hist=reco_1d_gen_binning,
                hist_truth=None,
                hist_mc_reco=hist_mc_reco_gen_binning,
                hist_mc_fakes=hist_mc_fakes_reco_gen_binning,
                hist_fake_fraction=hist_fake_fraction_gen_binning,
            )

            # Setup unfolder object
            # ---------------------
            axis_steering = '*[B]'
            # axis_steering = '*[]'
            if args.regularizeAxis == 'pt':
                axis_steering = 'pt[B];%s[N]' % variable_name
            elif args.regularizeAxis == 'angle':
                axis_steering = 'pt[N];%s[B]' % variable_name

            # disconnected_output_bins = find_disconnected_output_bins(hist_mc_gen_reco_map)
            # print("disconnected output bins", disconnected_output_bins)

            # disconnected_input_bins = find_disconnected_input_bins(hist_mc_gen_reco_map)
            # print("disconnected input bins", disconnected_input_bins)

            unfolder_args = dict(
                orientation=ROOT.TUnfold.kHistMapOutputHoriz,
                constraintMode=AREA_OPT_DICT[args.areaConstraint],
                # regMode=ROOT.TUnfold.kRegModeCurvature,
                # densityFlags=ROOT.TUnfoldDensity.kDensityModeBinWidth, # important as we have varying bin sizes!
                regMode=ROOT.TUnfold.kRegModeNone,
                densityFlags=ROOT.TUnfoldDensity.kDensityModeBinWidthAndUser,  # doesn't actually matter as RegModNone
                distribution='generatordistribution',  # the one to use for actual final regularisation/unfolding
                axisSteering=axis_steering
            )

            unfolder = MyUnfolder(response_map=rm_large_rel_error_bins_th2(hist_mc_gen_reco_map, args.relErr),
                                  binning_handler=binning_handler,
                                  **unfolder_args)

            # Look for unconstrained output bins,
            # and set the corresponding entries in the response matrix to 0
            # unconstrained_bins = unfolder.check_input(reco_1d)
            # for gen_bin, reco_bins in unconstrained_bins:
            #     print("unconstrained_output:", gen_bin, "=", binning_handler.global_bin_to_physical_bin(gen_bin, 'generator'), "depends on reco bins:")
            #     for r in reco_bins:
            #         print('    reco bin', r, "=", binning_handler.global_bin_to_physical_bin(r, 'detector'))
            #         hist_mc_gen_reco_map.SetBinContent(gen_bin, r, 0)
            #         hist_mc_gen_reco_map.SetBinError(gen_bin, r, 0)
            #
            # unfolder = MyUnfolder(response_map=hist_mc_gen_reco_map,
            #                       binning_handler=binning_handler,
            #                       **unfolder_args)

            # Save binning to file
            # unfolder.save_binning(txt_filename="%s/binning_scheme.txt" % (this_output_dir), print_xml=False)
            # ROOT.BinningXMLExporter.ExportXML(unfolder.detector_binning, this_output_dir, "detector_binning.xml", True, True, 2)
            # ROOT.BinningXMLExporter.ExportXML(unfolder.generator_binning, this_output_dir, "generator_binning.xml", True, True, 2)

            # SetEpsMatrix ensures rank properly calculated when inverting
            # Needed if you get message "rank of matrix E 55 expect 170"
            # And unfolded looks wacko
            eps_matrix = 1E-60
            unfolder.SetEpsMatrix(eps_matrix)
            print("Running with eps =", unfolder.GetEpsMatrix())

            # Set what is to be unfolded
            # ------------------------------------------------------------------
            def input_args(**kwargs):
                # make into func so can generate new copy of hists each time,
                # instead of lots of references - needed for bg_subtracted ones
                d = dict(input_handler=InputHandler(**input_handler_args),
                         input_handler_gen_binning=InputHandler(**input_handler_gen_binning_args),
                         bias_factor=args.biasFactor,
                         error_unconstrained_truth_bins=False)
                d.update(kwargs)
                return d

            unfolder.set_input(**input_args())

            # Store any functions used to rebin/merge TH1 and TH2s, so we can apply them later
            th1_merge_reco_funcs = []
            th1_merge_gen_funcs = []
            th2_merge_funcs = []

            # Merge last pt bin with 2nd to last, if desired
            # ------------------------------------------------------------------
            def setup_merged_last_pt_bin_unfolder(orig_unfolder):
                """Setup unfolder with last pt bin merged with 2nd to last"""
                this_binning_handler = orig_unfolder.binning_handler

                # Create new BinningHandler, bin edges, to use in merging funcs
                binning_handler_merged, gen_bin_pairs, reco_bin_pairs = setup_merged_last_pt_binning(this_binning_handler)

                new_gen_bins = np.arange(0.5, orig_unfolder.hist_truth.GetNbinsX() - len(gen_bin_pairs) + 0.5 + 1, 1) # new bin edges to match TUnfold style
                new_reco_bins = np.arange(0.5, orig_unfolder.input_hist.GetNbinsX() - len(reco_bin_pairs) + 0.5 + 1, 1) # new bin edges to match TUnfold style

                # Setup functions to do merging, since we always use same args
                th1_merge_reco = partial(merge_th1_bin_pairs, bin_pairs=reco_bin_pairs, new_bin_edges=new_reco_bins)
                th1_merge_gen = partial(merge_th1_bin_pairs, bin_pairs=gen_bin_pairs, new_bin_edges=new_gen_bins)
                th2_merge = partial(merge_th2_bin_pairs, bin_pairs_x=gen_bin_pairs, new_bin_edges_x=new_gen_bins,
                                    bin_pairs_y=reco_bin_pairs, new_bin_edges_y=new_reco_bins)

                th1_merge_reco_funcs.append(th1_merge_reco)
                th1_merge_gen_funcs.append(th1_merge_gen)
                th2_merge_funcs.append(th2_merge)

                response_map_merged = th2_merge(orig_unfolder.response_map)

                print('orig response map dim:', orig_unfolder.response_map.GetNbinsX(), orig_unfolder.response_map.GetNbinsY())
                print('rebinned response map dim:', response_map_merged.GetNbinsX(), response_map_merged.GetNbinsY())

                # Setup new unfolder
                new_unfolder = MyUnfolder(response_map=response_map_merged,
                                          binning_handler=binning_handler_merged,
                                          **unfolder_args)

                new_unfolder.SetEpsMatrix(unfolder.GetEpsMatrix())

                orig_input_handler = orig_unfolder.input_handler
                input_handler_merged = InputHandler(input_hist=th1_merge_reco(orig_input_handler.input_hist),
                                                    hist_truth=th1_merge_gen(orig_input_handler.hist_truth),
                                                    hist_mc_reco=th1_merge_reco(orig_input_handler.hist_mc_reco),
                                                    hist_mc_fakes=th1_merge_reco(orig_input_handler.hist_mc_fakes))

                orig_input_handler_gen_binning = orig_unfolder.input_handler_gen_binning
                input_handler_gen_binning_merged = InputHandler(input_hist=th1_merge_gen(orig_input_handler_gen_binning.input_hist),
                                                                hist_truth=None,
                                                                hist_mc_reco=th1_merge_gen(orig_input_handler_gen_binning.hist_mc_reco),
                                                                hist_mc_fakes=th1_merge_gen(orig_input_handler_gen_binning.hist_mc_fakes))

                # Set what is to be unfolded
                # ------------------------------------------------------------------
                new_unfolder.set_input(input_handler=input_handler_merged,
                                       input_handler_gen_binning=input_handler_gen_binning_merged,
                                       bias_factor=args.biasFactor,
                                       error_unconstrained_truth_bins=False)

                return new_unfolder

            if args.mergeLastPtBin:
                print(cu.pcolors.OKBLUE, "Doing last pt bin merging...", cu.pcolors.ENDC)
                unfolder = setup_merged_last_pt_bin_unfolder(unfolder)
                unfolder.check_input()

                # update input args templates for future setups
                input_handler_args = dict(
                    input_hist=unfolder.input_hist,
                    hist_truth=unfolder.hist_truth,
                    hist_mc_reco=unfolder.hist_mc_reco,
                    hist_mc_fakes=unfolder.hist_mc_fakes
                )
                input_handler_gen_binning_args = dict(
                    input_hist=unfolder.input_hist_gen_binning,
                    hist_truth=None,
                    hist_mc_reco=unfolder.hist_mc_reco_gen_binning,
                    hist_mc_fakes=unfolder.hist_mc_fakes_gen_binning
                )

            # Merge any bad columns in the probability matrix
            # ------------------------------------------------------------------
            def setup_merged_bin_unfolder(gen_bins_to_merge, reco_bins_to_merge, orig_unfolder):
                """"Merge bins, and setup MyUnfolder with pt-dependent binning"""
                gen_bins_to_merge = gen_bins_to_merge or []
                reco_bins_to_merge = reco_bins_to_merge or []
                binning_handler_merged = setup_merged_bin_binning(gen_bins_to_merge, reco_bins_to_merge, orig_unfolder.binning_handler)

                # Create new inputs
                new_gen_bins = np.arange(0.5, orig_unfolder.hist_truth.GetNbinsX() - len(gen_bins_to_merge) + 0.5 + 1, 1) # new bin edges to match TUnfold style
                new_reco_bins = np.arange(0.5, orig_unfolder.input_hist.GetNbinsX() - len(reco_bins_to_merge) + 0.5 + 1, 1) # new bin edges to match TUnfold style

                # Setup functions to do merging, since we always use same args
                th1_merge_reco = partial(merge_th1_bins, bin_list=reco_bins_to_merge, new_bin_edges=new_reco_bins)
                th1_merge_gen = partial(merge_th1_bins, bin_list=gen_bins_to_merge, new_bin_edges=new_gen_bins)
                th2_merge = partial(merge_th2_bins, bin_list_x=gen_bins_to_merge, new_bin_edges_x=new_gen_bins,
                                    bin_list_y=reco_bins_to_merge, new_bin_edges_y=new_reco_bins)

                th1_merge_reco_funcs.append(th1_merge_reco)
                th1_merge_gen_funcs.append(th1_merge_gen)
                th2_merge_funcs.append(th2_merge)

                response_map_merged = th2_merge(orig_unfolder.response_map)
                print('orig response map dim:', orig_unfolder.response_map.GetNbinsX(), orig_unfolder.response_map.GetNbinsY())
                print('rebinned response map dim:', response_map_merged.GetNbinsX(), response_map_merged.GetNbinsY())

                new_unfolder = MyUnfolder(response_map=response_map_merged,
                                          binning_handler=binning_handler_merged,
                                          **unfolder_args)

                new_unfolder.SetEpsMatrix(eps_matrix)

                orig_input_handler = orig_unfolder.input_handler
                input_handler_merged = InputHandler(input_hist=th1_merge_reco(orig_input_handler.input_hist),
                                                    hist_truth=th1_merge_gen(orig_input_handler.hist_truth),
                                                    hist_mc_reco=th1_merge_reco(orig_input_handler.hist_mc_reco),
                                                    hist_mc_fakes=th1_merge_reco(orig_input_handler.hist_mc_fakes))

                orig_input_handler_gen_binning = orig_unfolder.input_handler_gen_binning
                input_handler_gen_binning_merged = InputHandler(input_hist=th1_merge_gen(orig_input_handler_gen_binning.input_hist),
                                                                hist_truth=None,
                                                                hist_mc_reco=th1_merge_gen(orig_input_handler_gen_binning.hist_mc_reco),
                                                                hist_mc_fakes=th1_merge_gen(orig_input_handler_gen_binning.hist_mc_fakes))

                # Set what is to be unfolded
                # ------------------------------------------------------------------
                new_unfolder.set_input(input_handler=input_handler_merged,
                                       input_handler_gen_binning=input_handler_gen_binning_merged,
                                       bias_factor=args.biasFactor,
                                       error_unconstrained_truth_bins=False)

                return new_unfolder

            if args.mergeBadResponseBins:
                print(cu.pcolors.OKBLUE, "Doing bad probability bin merging...", cu.pcolors.ENDC)
                gen_bins_to_merge = get_bins_to_merge_probability_stats(unfolder.get_probability_matrix())
                unfolder = setup_merged_bin_unfolder(gen_bins_to_merge, None, unfolder)
                unfolder.check_input()

                 # update input args templates for future setups
                input_handler_args = dict(
                    input_hist=unfolder.input_hist,
                    hist_truth=unfolder.hist_truth,
                    hist_mc_reco=unfolder.hist_mc_reco,
                    hist_mc_fakes=unfolder.hist_mc_fakes
                )
                input_handler_gen_binning_args = dict(
                    input_hist=unfolder.input_hist_gen_binning,
                    hist_truth=None,
                    hist_mc_reco=unfolder.hist_mc_reco_gen_binning,
                    hist_mc_fakes=unfolder.hist_mc_fakes_gen_binning
                )

            # update some other hists
            for merge_func in th1_merge_reco_funcs:
                hist_mc_reco = merge_func(hist_mc_reco)

            for merge_func in th1_merge_reco_funcs:
                hist_mc_gen = merge_func(hist_mc_gen)

            unfolder.hist_bin_chopper.add_obj('hist_truth', unfolder.hist_truth)

            # Setup plotter
            # ------------------------------------------------------------------
            unfolder_plotter = MyUnfolderPlotter(unfolder,
                                                 is_data=not MC_INPUT,
                                                 lumi=cu.get_lumi_str(do_dijet="Dijet" in region['name'],
                                                                      do_zpj="ZPlusJets" in region['name']))
            plot_args = dict(output_dir=this_output_dir, append=append)

            # Add systematic errors as different response matrices, use TUnfold machinery
            # ------------------------------------------------------------------
            if args.doExperimentalSysts and not args.doExperimentalSystsAsAltResponse:
                chosen_rsp_bin = (78, 150)
                print("nominal response bin content for", chosen_rsp_bin, unfolder.response_map.GetBinContent(*chosen_rsp_bin))

                for syst_ind, syst_dict in enumerate(region['experimental_systematics']):
                    print("Adding systematic:", syst_dict['label'])
                    if 'factor' in syst_dict:
                        # special case for e.g. lumi - we construct a reponse hist, and add it using relative mode
                        rel_map = unfolder.response_map.Clone(syst_dict['label']+"Map")
                        for xbin, ybin in product(range(1, rel_map.GetNbinsX()+1), range(1, rel_map.GetNbinsY()+1)):
                            rel_map.SetBinContent(xbin, ybin, syst_dict['factor'])
                            rel_map.SetBinError(xbin, ybin, 0)
                        unfolder.add_sys_error(rel_map, syst_dict['label'], ROOT.TUnfoldDensity.kSysErrModeRelative)

                    # elif 'herwig' in syst_dict['label'].lower():
                    #     if not isinstance(syst_dict['tfile'], ROOT.TFile):
                    #         syst_dict['tfile'] = cu.open_root_file(syst_dict['tfile'])
                    #     map_syst = cu.get_from_tfile(syst_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                    #     map_syst.Scale(unfolder.response_map.Integral() / map_syst.Integral())
                    #     print("    syst bin", chosen_rsp_bin, map_syst.GetBinContent(*chosen_rsp_bin))
                    #     unfolder.add_sys_error(rm_large_rel_error_bins_th2(map_syst), syst_dict['label'], ROOT.TUnfoldDensity.kSysErrModeMatrix)

                    else:
                        if not isinstance(syst_dict['tfile'], ROOT.TFile):
                            syst_dict['tfile'] = cu.open_root_file(syst_dict['tfile'])
                        map_syst = cu.get_from_tfile(syst_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                        # update binnings
                        for merge_func in th2_merge_funcs:
                            map_syst = merge_func(map_syst)
                        print("    syst bin", chosen_rsp_bin, map_syst.GetBinContent(*chosen_rsp_bin))
                        unfolder.add_sys_error(rm_large_rel_error_bins_th2(map_syst), syst_dict['label'], ROOT.TUnfoldDensity.kSysErrModeMatrix)

                        # Plot the response matrix for this systematic
                        syst_label_no_spaces = cu.no_space_str(syst_dict['label'])
                        output_filename = "%s/response_map_syst_%s_%s.%s" % (this_output_dir, syst_label_no_spaces, append, unfolder_plotter.output_fmt)
                        title = "%s, %s region, %s, %s" % (jet_algo, region['label'], angle_str, syst_dict['label'])
                        unfolder_plotter.draw_2d_hist(map_syst,
                                                      title=title,
                                                      output_filename=output_filename,
                                                      logz=True,
                                                      draw_bin_lines_x=True,
                                                      draw_bin_lines_y=True,
                                                      canvas_size=(800, 700),
                                                      xtitle='Generator bin',
                                                      ytitle='Detector bin')

            # Subtract fakes (treat as background)
            # ------------------------------------------------------------------
            hist_fakes_reco = unfolder.input_handler.calc_fake_hist(unfolder.input_handler.input_hist)
            unfolder.subtract_background(hist_fakes_reco, "Signal fakes")
            unfolder.input_handler.hist_mc_reco_bg_subtracted = unfolder.input_handler.subtract_fake_hist(unfolder.input_handler.hist_mc_reco)[0]

            hist_fakes_reco_gen_binning = unfolder.input_handler_gen_binning.calc_fake_hist(unfolder.input_handler_gen_binning.input_hist)
            unfolder.subtract_background_gen_binning(hist_fakes_reco_gen_binning, "Signal fakes")
            unfolder.input_handler_gen_binning.hist_mc_reco_bg_subtracted = unfolder.input_handler_gen_binning.subtract_fake_hist(unfolder.input_handler_gen_binning.hist_mc_reco)[0]

            # Subtract actual backgrounds if necessary
            # ------------------------------------------------------------------
            background_reco_1d = None
            background_gen_1d = None
            if "backgrounds" in region and args.subtractBackgrounds:
                for bg_ind, bg_dict in enumerate(region['backgrounds']):
                    print("Subtracting", bg_dict['name'], 'background')
                    if not isinstance(bg_dict['tfile'], ROOT.TFile):
                        bg_dict['tfile'] = cu.open_root_file(bg_dict['tfile'])

                    mc_hname_append = "split" if MC_SPLIT else "all"
                    bg_hist = cu.get_from_tfile(bg_dict['tfile'], "%s/hist_%s_reco_%s" % (region['dirname'], angle_shortname, mc_hname_append))
                    bg_hist_gen = cu.get_from_tfile(bg_dict['tfile'], "%s/hist_%s_truth_%s" % (region['dirname'], angle_shortname, mc_hname_append))
                    for merge_func in th1_merge_reco_funcs:
                        bg_hist = merge_func(bg_hist)
                    for merge_func in th1_merge_gen_funcs:
                        bg_hist_gen = merge_func(bg_hist_gen)
                    bg_dict['hist'] = bg_hist
                    bg_dict['hist_gen'] = bg_hist_gen

                    # keep one big hist
                    if not background_reco_1d:
                        background_reco_1d = bg_hist.Clone()
                        background_gen_1d = bg_hist_gen.Clone()
                    else:
                        background_reco_1d.Add(bg_hist, bg_dict.get('rate', 1.))
                        background_gen_1d.Add(bg_hist_gen, bg_dict.get('rate', 1.))

                    unfolder.subtract_background(hist=bg_hist,
                                                 name=bg_dict['name'],
                                                 scale=bg_dict.get('rate', 1.),
                                                 scale_err=bg_dict.get('rate_unc', 0.))

            title = "%s\n%s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_background_fractions(title=title, **plot_args)

            unfolder.print_condition_number()

            # Setup hists from alternative generator
            # ------------------------------------------------------------------
            alt_unfolder = None
            alt_hist_mc_gen = None  # mc truth of the alternate generator used to make reponse matrix
            alt_hist_mc_reco = None  # reco of the alternate generator used to make reponse matrix
            alt_hist_mc_reco_bg_subtracted = None
            alt_hist_mc_reco_bg_subtracted_gen_binning = None

            # Do this outside the if statement, since we might use it later in plotting e.g. for data
            if not isinstance(region['alt_mc_tfile'], ROOT.TFile) and os.path.isfile(region['alt_mc_tfile']):
                region['alt_mc_tfile'] = cu.open_root_file(region['alt_mc_tfile'])

                alt_hist_mc_gen = cu.get_from_tfile(region['alt_mc_tfile'], "%s/hist_%s_truth_all" % (region['dirname'], angle_shortname))
                alt_hist_mc_reco = cu.get_from_tfile(region['alt_mc_tfile'], "%s/hist_%s_reco_all" % (region['dirname'], angle_shortname))

                alt_scale = unfolder.hist_truth.Integral() / alt_hist_mc_gen.Integral()
                alt_hist_mc_gen.Scale(alt_scale)
                alt_hist_mc_reco.Scale(alt_scale)

                # Update binnings
                for merge_func in th1_merge_gen_funcs:
                    alt_hist_mc_gen = merge_func(alt_hist_mc_gen)
                for merge_func in th1_merge_reco_funcs:
                    alt_hist_mc_reco = merge_func(alt_hist_mc_reco)

                # fakes-subtracted version
                alt_hist_mc_reco_bg_subtracted, alt_hist_fakes = unfolder.input_handler.subtract_fake_hist(alt_hist_mc_reco)

                # gen-binned versions of detector-level plots
                alt_hist_mc_reco_gen_binning = cu.get_from_tfile(region['alt_mc_tfile'], "%s/hist_%s_reco_gen_binning" % (region['dirname'], angle_shortname))
                alt_hist_mc_reco_gen_binning.Scale(alt_scale)
                for merge_func in th1_merge_gen_funcs:
                    alt_hist_mc_reco_gen_binning = merge_func(alt_hist_mc_reco_gen_binning)
                alt_hist_mc_reco_bg_subtracted_gen_binning, alt_hist_fakes_gen_binning = unfolder.input_handler_gen_binning.subtract_fake_hist(alt_hist_mc_reco_gen_binning)

                unfolder.hist_bin_chopper.add_obj('alt_hist_truth', alt_hist_mc_gen)

            # Do any regularization
            # ---------------------
            unreg_unfolder = None
            L_matrix_entries = []
            tau = 0
            scan_mode = ROOT.TUnfoldDensity.kEScanTauRhoAvgSys
            scan_mode = ROOT.TUnfoldDensity.kEScanTauRhoAvg
            scan_distribution = unfolder.distribution

            if REGULARIZE != "None":
                print("Doing preliminary unregularised unfolding...")
                # Do an unregularised version first for comparison,
                # also potentially for factor/bias
                unreg_unfolder = MyUnfolder(response_map=unfolder.response_map,
                                            binning_handler=binning_handler,
                                            **unfolder_args)

                unreg_unfolder.SetEpsMatrix(eps_matrix)

                unreg_unfolder_plotter = MyUnfolderPlotter(unreg_unfolder, is_data=not MC_INPUT)

                # Do the unregularised unfolding to get an idea of bin contents
                # and uncertainties
                # Set what is to be unfolded
                # ------------------------------------------------------------------
                unreg_unfolder.set_input(**input_args(bias_factor=0))

                # For now, ignore experimental systematics

                # Subtract fakes (treat as background)
                # ------------------------------------------------------------------
                for key, h_bkg in unfolder.backgrounds.items():
                    unreg_unfolder.subtract_background(h_bkg, key)
                for key, h_bkg in unfolder.backgrounds_gen_binning.items():
                    unreg_unfolder.subtract_background_gen_binning(h_bkg, key)

                # Do the unfolding
                # ------------------------------------------------------------------
                unreg_unfolder.do_unfolding(0)
                unreg_unfolder.get_output(hist_name="unreg_unfolded_1d")
                unreg_unfolder._post_process()

                region['unreg_unfolder'] = unreg_unfolder

                if not (alt_hist_mc_gen and alt_hist_mc_reco_bg_subtracted):
                    raise RuntimeError("Cannot create truth template for regularisation as alt MC missing")

                truth_template = hist_mc_gen_all

                if args.biasVector == "template" or args.biasVector == "templateOtherProc":
                    # Create truth template by fitting MC to data @ detector level
                    # --------------------------------------------------------------
                    # Fit the two MC templates to data to get their fractions
                    # Then use the same at truth level
                    # This template will allow us to setup a more accurate L matrix,
                    # and a bias hist
                    template_maker = TruthTemplateMaker(binning_handler=binning_handler,
                                                        output_dir=this_output_dir)

                    # thing to be fitted
                    template_maker.set_input(unreg_unfolder.input_hist_gen_binning_bg_subtracted)

                    # templates to do the fitting, and to create the truth-level distribution
                    template_maker.add_mc_template(name=region['mc_label'],
                                                   hist_reco=unreg_unfolder.input_handler.hist_mc_reco_bg_subtracted,
                                                   hist_gen=hist_mc_gen_all,
                                                   colour=ROOT.kRed)
                    template_maker.add_mc_template(name=region['alt_mc_label'],
                                                   hist_reco=alt_hist_mc_reco_bg_subtracted_gen_binning,
                                                   hist_gen=alt_hist_mc_gen,
                                                   colour=ROOT.kViolet+1)

                    if args.biasVector == "templateOtherProc":
                        # add the "other" process (e.g. QCD for DY)
                        # this accounts for the assumption that not only are the shapes wrong
                        # (which we account for by using MG+Pythia & H++),
                        # but the fraction wrong (account for by using the other process)
                        if isinstance(region['mc_otherProc_tfile'], str):
                            print("Opening", region['mc_otherProc_tfile'])
                            region['mc_otherProc_tfile'] = cu.open_root_file(region['mc_otherProc_tfile'])

                        hist_mc_otherProc_gen_all = cu.get_from_tfile(region['mc_otherProc_tfile'], "%s/hist_%s_truth_all" % (region['dirname_otherProc'], angle_shortname))
                        hist_mc_otherProc_reco_gen_binning_all = cu.get_from_tfile(region['mc_otherProc_tfile'], "%s/hist_%s_reco_gen_binning" % (region['dirname_otherProc'], angle_shortname))
                        hist_mc_otherProc_fakes_reco_gen_binning = cu.get_from_tfile(region['mc_otherProc_tfile'], "%s/hist_%s_reco_fake_gen_binning" % (region['dirname_otherProc'], angle_shortname))
                        # FIXME need to apply th1_merge_funcs
                        hist_otherProc_fake_fraction_gen_binning = hist_mc_otherProc_fakes_reco_gen_binning.Clone("hist_%s_fakes_fraction_gen_binning" % angle_shortname)
                        hist_otherProc_fake_fraction_gen_binning.Divide(hist_mc_otherProc_reco_gen_binning_all)
                        hist_mc_otherProc_reco_gen_binning_all_bg_subtracted, hist_otherProc_fakes_reco_all_gen_binning = subtract_background(hist_mc_otherProc_reco_gen_binning_all, hist_otherProc_fake_fraction_gen_binning)

                        template_maker.add_mc_template(name=region['mc_otherProc_label'],
                                                       hist_reco=hist_mc_otherProc_reco_gen_binning_all_bg_subtracted,
                                                       hist_gen=hist_mc_otherProc_gen_all,
                                                       colour=ROOT.kRed+3)

                        if isinstance(region['alt_mc_otherProc_tfile'], str):
                            print("Opening", region['alt_mc_otherProc_tfile'])
                            region['alt_mc_otherProc_tfile'] = cu.open_root_file(region['alt_mc_otherProc_tfile'])

                        alt_hist_mc_otherProc_gen_all = cu.get_from_tfile(region['alt_mc_otherProc_tfile'], "%s/hist_%s_truth_all" % (region['dirname_otherProc'], angle_shortname))
                        alt_hist_mc_otherProc_reco_gen_binning_all = cu.get_from_tfile(region['alt_mc_otherProc_tfile'], "%s/hist_%s_reco_gen_binning" % (region['dirname_otherProc'], angle_shortname))
                        alt_hist_mc_otherProc_fakes_reco_gen_binning = cu.get_from_tfile(region['alt_mc_otherProc_tfile'], "%s/hist_%s_reco_fake_gen_binning" % (region['dirname_otherProc'], angle_shortname))
                        alt_hist_otherProc_fake_fraction_gen_binning = alt_hist_mc_otherProc_fakes_reco_gen_binning.Clone("hist_%s_fakes_fraction_gen_binning" % angle_shortname)
                        alt_hist_otherProc_fake_fraction_gen_binning.Divide(alt_hist_mc_otherProc_reco_gen_binning_all)
                        alt_hist_mc_otherProc_reco_gen_binning_all_bg_subtracted, hist_otherProc_fakes_reco_all_gen_binning = subtract_background(alt_hist_mc_otherProc_reco_gen_binning_all, alt_hist_otherProc_fake_fraction_gen_binning)

                        template_maker.add_mc_template(name=region['alt_mc_otherProc_label'],
                                                       hist_reco=alt_hist_mc_otherProc_reco_gen_binning_all_bg_subtracted,
                                                       hist_gen=alt_hist_mc_otherProc_gen_all,
                                                       colour=ROOT.kViolet+4)

                    truth_template = template_maker.create_template()

                    if MC_INPUT:
                        # do check by fitting at gen level and comparing
                        template_maker.set_input_gen(unreg_unfolder.hist_truth)
                        template_maker.check_template_at_gen()

                elif args.biasVector == "alttruth":
                    truth_template = alt_hist_mc_gen

                unfolder.truth_template = truth_template
                unfolder.hist_bin_chopper.add_obj("truth_template", truth_template)
                unreg_unfolder.truth_template = truth_template
                unreg_unfolder.hist_bin_chopper.add_obj("truth_template", truth_template)

                # Draw our new template
                ocs = [
                    Contribution(alt_hist_mc_gen, label=region['alt_mc_label'],
                                 line_color=ROOT.kViolet+1,
                                 marker_color=ROOT.kViolet+1,
                                 subplot=unreg_unfolder.hist_truth),
                    Contribution(truth_template, label="Template",
                                 line_color=ROOT.kAzure+1,
                                 marker_color=ROOT.kAzure+1,
                                 subplot=unreg_unfolder.hist_truth),
                ]
                # title = "%s\n%s region, %s" % (jet_algo, region['label'], angle_str)
                unreg_unfolder_plotter.draw_unfolded_1d(do_gen=True,
                                                        do_unfolded=True,
                                                        other_contributions=ocs,
                                                        output_dir=this_output_dir,
                                                        append='unreg_with_template',
                                                        title='',
                                                        subplot_title="* / Generator",
                                                        subplot_limits=(0, 2))

                # Setup our L matrix
                # ------------------------------------------------------------------
                # ref_hist = unreg_unfolder.hist_truth
                unfolder.SetBias(truth_template)
                unfolder.setup_L_matrix_curvature(ref_hist=truth_template, axis=args.regularizeAxis)

                # Get bin factors from an unregularised unfolding first,
                # to compensate for the fact that the shape differs between data & MC
                # bin_factors = unreg_unfolder.calculate_pt_bin_factors(which='gen') # calculate factors to get uniform pt spectrum
                # bin_factors = unreg_unfolder.calculate_pt_bin_factors(which='unfolded') # calculate factors to get uniform pt spectrum
                # bin_widths = unreg_unfolder.get_gen_bin_widths() # mapping {global bin number : (lambda bin width, pt bin width)}

                # loop over existing regularisation conditions, since we want to modify them
                # in our main unfolder
                # orig_Lmatrix = unreg_unfolder.GetL("orig_Lmatrix_%s" % (append), "", unreg_unfolder.use_axis_binning)
                # for iy in range(1, orig_Lmatrix.GetNbinsY()+1):
                #     # Look for gen bin number where values start for this regularisation row
                #     left_bin, mid_bin, right_bin = 0, 0, 0
                #     left_bin_val, mid_bin_val, right_bin_val = 0, 0, 0
                #     for ix in range(1, orig_Lmatrix.GetNbinsX()+1):
                #         bin_content = orig_Lmatrix.GetBinContent(ix, iy)
                #         if bin_content != 0:
                #             if left_bin == 0:
                #                 left_bin = ix
                #                 left_bin_val = bin_content
                #                 continue
                #             elif mid_bin == 0:
                #                 mid_bin = ix
                #                 mid_bin_val = bin_content
                #                 continue
                #             else:
                #                 right_bin = ix
                #                 right_bin_val = bin_content
                #                 break # got em all

                    # Things to try:
                    # - Ignore the original reg. condition: just use bin_factors
                    # - Same, but with pt bin_width compensated
                    # - Modify existing reg.condition with bin_factors
                    # - Same, but compensate for pt bin width
                    # - Weighting my median relative error? e.g. care more about
                    # regions/bins with larger error bars

                    # Rescale accord to pT, and according to pT bin width
                    # since the original was divided by both pt bin width and lambda bin width
                    # Doesn't matter left or right for bin widht - only care about pt bin width
                    # pt_factor = bin_factors[mid_bin] * bin_widths[left_bin][1]
                    # pt_factor = bin_factors[mid_bin]

                    # - signs since RegularizeCurvature also adds in a - sign,
                    # and we want to match the original sign (which is +ve)
                    # scale_left = -pt_factor
                    # scale_right = -pt_factor
                    # scale_left = -left_bin_val * pt_factor
                    # scale_right = -right_bin_val * pt_factor

                    # print("Adding regularisation rule: nR=%d, gen bins: [%d - %d], factors: [%f, %f, %f]" % (iy, left_bin, right_bin, -scale_left, 2*(scale_left+scale_right), -scale_right) )
                    # unfolder.RegularizeCurvature(left_bin, mid_bin, right_bin, scale_left, scale_right)

                if REGULARIZE == "L":
                    print("Regularizing with ScanLcurve, please be patient...")
                    l_scanner = LCurveScanner()
                    tau = l_scanner.scan_L(tunfolder=unfolder,
                                           n_scan=args.nScan,
                                           tau_min=region['tau_limits'][angle.var][0],
                                           tau_max=region['tau_limits'][angle.var][1])
                    print("Found tau:", tau)
                    l_scanner.plot_scan_L_curve(output_filename="%s/scanL_%s.%s" % (this_output_dir, append, OUTPUT_FMT))
                    l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

                elif REGULARIZE == "tau":
                    print("Regularizing with ScanTau, please be patient...")
                    tau_scanner = TauScanner()
                    tau = tau_scanner.scan_tau(tunfolder=unfolder,
                                               n_scan=args.nScan,
                                               tau_min=region['tau_limits'][angle.var][0],
                                               tau_max=region['tau_limits'][angle.var][1],
                                               scan_mode=scan_mode,
                                               distribution=scan_distribution,
                                               axis_steering=unfolder.axisSteering)
                    print("Found tau:", tau)
                    tau_scanner.plot_scan_tau(output_filename="%s/scantau_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

                title = "L matrix %s %s region %s" % (jet_algo, region['label'], angle_str)
                unfolder_plotter.draw_L_matrix(title=title, **plot_args)
                title = "L^{T}L matrix %s %s region %s" % (jet_algo, region['label'], angle_str)
                unfolder_plotter.draw_L_matrix_squared(title=title, **plot_args)
                title = "L * (x - bias vector)\n%s\n%s region\n%s" % (jet_algo, region['label'], angle_str)
                unfolder_plotter.draw_Lx_minus_bias(title=title, **plot_args)

            # Do unfolding!
            # ------------------------------------------------------------------
            unfolder.do_unfolding(tau)
            unfolder.get_output(hist_name="unfolded_1d")

            # Do lots of extra gubbins, like caching matrices,
            # creating unfolded hists with different levels of uncertainties,
            # ------------------------------------------------------------------
            unfolder._post_process()

            # Check result with numpy
            # unfolder.do_numpy_comparison(output_dir=this_output_dir)

            # Figure out if we should merge any bins (i.e. -ve ones)
            # ------------------------------------------------------------------
            # have to do it this way to ensure combined merged bins end up correct
            # TODO is it even possible to convert into a 1D list?
            merge_bin_history = []
            merge_dir = os.path.join(this_output_dir, 'merge_bins')

            if args.mergeBins or args.mergeBinsFromFile:
                unfolder_merged = unfolder
                reco_bins_to_merge = []
                tau_merged = 0

                if args.mergeBins:
                    print(cu.pcolors.OKBLUE, "Doing bin merging...", cu.pcolors.ENDC)
                    # pt_overflow_bin = binning_handler.get_first_pt_overflow_global_bin("generator")
                    gen_bins_to_merge = get_bins_to_merge_probability_stats(unfolder.get_probability_matrix())

                    # gen_bins_to_merge2 = get_bins_to_merge_nonnegative(h=unfolder.get_output(),
                    #                                                   binning_obj=binning_handler.get_binning_scheme('generator'),
                    #                                                   ensure_nonnegative=False)
                    # gen_bins_to_merge.extend(gen_bins_to_merge2)
                    counter = 0
                    max_iterations = 1
                    while len(gen_bins_to_merge) > 0 and counter < max_iterations:
                        print("********** MERGE LOOP:", counter, "***********")
                        counter += 1
                        if len(gen_bins_to_merge) > 0:
                            merge_bin_history.append(gen_bins_to_merge)

                        unfolder_merged = setup_merged_bin_unfolder(gen_bins_to_merge, reco_bins_to_merge, unfolder_merged)
                        unfolder_merged.check_input()
                        unfolder_merged.do_unfolding(tau_merged)
                        # gen_bins_to_merge = get_bins_to_merge_probability_stats(unfolder_merged.get_probability_matrix())
                        gen_bins_to_merge = get_bins_to_merge_nonnegative(h=unfolder_merged.get_output(),
                                                                           binning_obj=unfolder_merged.binning_handler.get_binning_scheme('generator'),
                                                                           ensure_nonnegative=False)
                        # gen_bins_to_merge.extend(gen_bins_to_merge2)
                        reco_bins_to_merge = []

                    print("merge_bin_history:", merge_bin_history)

                    save_merge_bin_history_to_file(merge_bin_history, merge_dir)

                elif args.mergeBinsFromFile:
                    merge_bin_history = get_merge_bin_history_from_file(os.path.join(args.mergeBinsFromFile, region['name'], angle.var, 'merge_bins'))
                    print("Loaded merge_bin_history:", merge_bin_history)
                    for gen_bins_to_merge in merge_bin_history:
                        unfolder_merged = setup_merged_bin_unfolder(gen_bins_to_merge, reco_bins_to_merge, unfolder_merged)
                        unfolder_merged.do_unfolding(tau_merged)

                unfolder_plotter2 = MyUnfolderPlotter(unfolder_merged,
                                                      is_data=not MC_INPUT,
                                                      lumi=cu.get_lumi_str(do_dijet="Dijet" in region['name'],
                                                                           do_zpj="ZPlusJets" in region['name']))
                plot_args2 = dict(output_dir=merge_dir, append=append)

                # Do lots of extra gubbins, like caching matrices,
                # creating unfolded hists with different levels of uncertainties,
                # ------------------------------------------------------------------
                unfolder_merged._post_process()
                unfolder_merged.print_condition_number()

                # Draw big 1D distributions
                # ------------------------------------------------------------------
                title = "%s\n%s region, %s" % (jet_algo, region['label'], angle_str)
                unfolder_plotter2.draw_unfolded_1d(title=title,
                                                   do_gen=True,
                                                   do_logy=True,
                                                   **plot_args2)

                # reco using detector binning
                unfolder_plotter2.draw_detector_1d(do_reco_mc=True,
                                                   do_reco_data=not MC_INPUT,
                                                   title=title,
                                                   **plot_args2)

                # reco using gen binning
                unfolder_plotter2.draw_generator_1d(do_reco_data=not MC_INPUT,
                                                    do_reco_data_bg_sub=False,
                                                    do_reco_bg=False,
                                                    do_reco_mc=True,
                                                    do_reco_mc_bg_sub=False,
                                                    do_truth_mc=True,
                                                    title=title,
                                                    **plot_args2)

                # same plot but with bg-subtracted reco
                unfolder_plotter2.draw_detector_1d(do_reco_data_bg_sub=not MC_INPUT,
                                                   do_reco_bg=True,
                                                   do_reco_mc_bg_sub=True,
                                                   output_dir=plot_args2['output_dir'],
                                                   append='bg_fakes_subtracted_%s' % append,
                                                   title=title)

                # same but with generator-binning
                unfolder_plotter2.draw_generator_1d(do_reco_data=False,
                                                    do_reco_data_bg_sub=not MC_INPUT,
                                                    do_reco_bg=True,
                                                    do_reco_mc=False,
                                                    do_reco_mc_bg_sub=True,
                                                    do_truth_mc=True,
                                                    output_dir=plot_args2['output_dir'],
                                                    append='bg_fakes_subtracted_%s' % append,
                                                    title=title)
                print("unfolder_plotter2 done")

                title = "Response matrix, %s, %s region, %s" % (jet_algo, region['label'], angle_str)
                unfolder_plotter2.draw_response_matrix(title=title, **plot_args2)

                # title = "Response matrix normed by detector p_{T} bin, %s, %s region, %s" % (jet_algo, region['label'], angle_str)
                # unfolder_plotter2.draw_response_matrix_normed_by_detector_pt(title=title, **plot_args2)

                title = ("#splitline{Probability matrix, %s, %s region, %s}{Condition number: #sigma_{max} / #sigma_{min} = %.3g / %.3g = %g}"
                         % (jet_algo, region['label'], angle_str, unfolder.sigma_max, unfolder.sigma_min, unfolder.condition_number))
                unfolder_plotter2.draw_probability_matrix(title=title, **plot_args2)

                unfolder_merged.setup_normalised_results_per_pt_bin()
                unfolder_merged.setup_absolute_results_per_pt_bin()
                region['unfolder'] = unfolder_merged
                setup2 = Setup(jet_algo=jet_algo,
                               region=region,
                               angle=angle,
                               output_dir=plot_args2['output_dir'],
                               has_data=not MC_INPUT)
                hbc2 = do_binned_plots_per_region_angle(setup2,
                                                        do_binned_gen_pt=True,
                                                        do_binned_gen_lambda=False,
                                                        do_binned_reco_pt=True)

                region['unfolder'] = unfolder

            prof_end_nominal()

            # ------------------------------------------------------------------
            # CALCULATE JACKKNIFED UNCERTAINTIES
            # ------------------------------------------------------------------
            if args.doJackknifeInput:
                # first construct all new systematic variations dicts
                original_jk_dict = region['jackknife_input_variations'][0]
                original_jk_dict['label'] = '_jackknife_template'  # initial _ to ignore it later on
                # tfile = original_jk_dict['tfile']
                # if not isinstance(tfile, ROOT.TFile):
                #     tfile = cu.open_root_file(tfile)

                region['jackknife_input_variations'] = []
                num_vars = len(original_jk_dict['variations'])
                input_tfile = region['mc_tfile'] if MC_INPUT else region['data_tfile']
                for jk_ind in original_jk_dict['variations']:
                    region['jackknife_input_variations'].append(
                        {
                            "label": "Jackknife_input_%d" % jk_ind,
                            "input_reco": cu.get_from_tfile(input_tfile, "%s/hist_%s_reco_all_jackknife_%d" % (region['dirname'], angle_shortname, jk_ind)),
                            "input_gen": cu.get_from_tfile(input_tfile, "%s/hist_%s_truth_all_jackknife_%d" % (region['dirname'], angle_shortname, jk_ind)) if MC_INPUT else hist_mc_gen,
                            "colour": cu.get_colour_seq(jk_ind, num_vars),
                        })

                # Now run over all variations, unfolding the various inputs
                # but with the nominal response matrices
                for jk_ind, jk_dict in enumerate(region['jackknife_input_variations']):
                    jk_label = jk_dict['label']
                    jk_label_no_spaces = cu.no_space_str(jk_label)

                    print("*" * 80)
                    print("*** Unfolding with jackknife input:", jk_label, "(%d/%d) ***" % (jk_ind+1, len(region['jackknife_input_variations'])))
                    print("*" * 80)

                    jk_unfolder = MyUnfolder(response_map=unfolder.response_map,
                                             binning_handler=binning_handler,
                                             **unfolder_args)

                    jk_unfolder.SetEpsMatrix(eps_matrix)

                    jk_unfolder_plotter = MyUnfolderPlotter(jk_unfolder, is_data=not MC_INPUT)
                    jk_output_dir = os.path.join(this_output_dir, "jackknife_input", jk_label_no_spaces)
                    jk_plot_args = dict(output_dir=jk_output_dir, append=append)

                    # Set what is to be unfolded
                    # --------------------------------------------------------------
                    jk_input_hist = jk_dict['input_reco']

                    jk_hist_mc_reco = rm_large_rel_error_bins_th1(jk_input_hist.Clone(), args.relErr)

                    # fakes-subtracted version
                    jk_hist_mc_reco_bg_subtracted, jk_hist_fakes = subtract_background(jk_hist_mc_reco, hist_fake_fraction)

                    # gen-binned versions of detector-level plots
                    # this is tricky - they don't exist in the ROOT file, so we'll have to construct it ourselves
                    # luckily only needed for regularisation template so not too crucial
                    jk_hist_mc_reco_gen_binning = jk_unfolder.convert_reco_binned_hist_to_gen_binned(jk_hist_mc_reco)

                    # fakes-subtracted version
                    jk_hist_mc_reco_bg_subtracted_gen_binning, jk_hist_fakes_gen_binning = subtract_background(jk_hist_mc_reco_gen_binning, hist_fake_fraction_gen_binning)

                    input_handler_jk = InputHandler(input_hist=jk_input_hist,
                                                    hist_truth=jk_dict['input_gen'].Clone(),
                                                    hist_mc_reco=jk_hist_mc_reco,
                                                    hist_mc_fakes=unfolder.input_handler.hist_mc_fakes,
                                                    hist_fake_fraction=hist_fake_fraction)

                    input_handler_gen_binning_jk = InputHandler(input_hist=jk_hist_mc_reco_gen_binning,
                                                                hist_truth=None,
                                                                hist_mc_reco=jk_hist_mc_reco_gen_binning,
                                                                hist_mc_fakes=unfolder.input_handler_gen_binning.hist_mc_fakes,
                                                                hist_fake_fraction=hist_fake_fraction_gen_binning)

                    jk_unfolder.set_input(input_handler=input_handler_jk,
                                           input_handler_gen_binning=input_handler_gen_binning_jk,
                                           bias_factor=args.biasFactor,
                                           error_unconstrained_truth_bins=False)

                    # Subtract fakes (treat as background)
                    # --------------------------------------------------------------
                    jk_unfolder.subtract_background(jk_hist_fakes, "Signal fakes", scale=1., scale_err=0.0)
                    jk_unfolder.subtract_background_gen_binning(jk_hist_fakes_gen_binning, "Signal fakes", scale=1.)

                    # Do regularisation
                    # --------------------------------------------------------------
                    # Since different input need to redo the template creation
                    jk_tau = 0
                    if REGULARIZE != "None":
                        jk_truth_template = jk_dict['input_gen']
                        if args.biasVector == "template":
                            # Create truth template by fitting MC to data @ detector level
                            # --------------------------------------------------------------
                            # Fit the two MC templates to data to get their fractions
                            # Then use the same at truth level
                            # This template will allow us to setup a more accurate L matrix,
                            # and a bias hist
                            jk_template_maker = TruthTemplateMaker(binning_handler=binning_handler,
                                                                   output_dir=jk_output_dir)

                            jk_template_maker.set_input(jk_unfolder.input_hist_gen_binning_bg_subtracted)

                            jk_template_maker.add_mc_template(name=region['mc_label'],
                                                              hist_reco=unfolder.input_handler_gen_binning.hist_mc_reco_bg_subtracted,
                                                              hist_gen=hist_mc_gen_all,
                                                              colour=ROOT.kRed)
                            jk_template_maker.add_mc_template(name=region['alt_mc_label'],
                                                              hist_reco=alt_hist_mc_reco_bg_subtracted_gen_binning,
                                                              hist_gen=alt_hist_mc_gen,
                                                              colour=ROOT.kViolet+1)

                            jk_truth_template = jk_template_maker.create_template()

                            if MC_INPUT:
                                # do check by fitting at gen level and comparing
                                jk_template_maker.set_input_gen(jk_unfolder.hist_truth)
                                jk_template_maker.check_template_at_gen()

                        elif args.biasVector == "alttruth":
                            jk_truth_template = alt_hist_mc_gen

                        jk_unfolder.truth_template = jk_truth_template
                        jk_unfolder.hist_bin_chopper.add_obj("truth_template", jk_unfolder.truth_template)
                        jk_unfolder.hist_bin_chopper.add_obj("alt_hist_truth", alt_hist_mc_gen)

                        jk_unfolder.SetBias(jk_unfolder.truth_template)
                        jk_unfolder.setup_L_matrix_curvature(ref_hist=jk_unfolder.truth_template, axis=args.regularizeAxis)

                        if REGULARIZE == "L":
                            print("Regularizing with ScanLcurve, please be patient...")
                            jk_l_scanner = LCurveScanner()
                            jk_tau = jk_l_scanner.scan_L(tunfolder=jk_unfolder,
                                                         n_scan=args.nScan,
                                                         tau_min=region['tau_limits'][angle.var][0],
                                                         tau_max=region['tau_limits'][angle.var][1])
                            print("Found tau:", jk_tau)
                            jk_l_scanner.plot_scan_L_curve(output_filename="%s/scanL_%s.%s" % (jk_output_dir, append, OUTPUT_FMT))
                            jk_l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_%s.%s" % (jk_output_dir, append, OUTPUT_FMT))

                        elif REGULARIZE == "tau":
                            print("Regularizing with ScanTau, please be patient...")
                            jk_tau_scanner = TauScanner()
                            jk_tau = jk_tau_scanner.scan_tau(tunfolder=jk_unfolder,
                                                             n_scan=args.nScan,
                                                             tau_min=region['tau_limits'][angle.var][0],
                                                             tau_max=region['tau_limits'][angle.var][1],
                                                             scan_mode=scan_mode,
                                                             distribution=scan_distribution,
                                                             axis_steering=unfolder.axisSteering)
                            print("Found tau:", jk_tau)
                            jk_tau_scanner.plot_scan_tau(output_filename="%s/scantau_%s.%s" % (jk_output_dir, append, OUTPUT_FMT))

                        title = "L matrix %s %s region %s" % (jet_algo, region['label'], angle_str)
                        jk_unfolder_plotter.draw_L_matrix(title=title, **jk_plot_args)
                        title = "L^{T}L matrix %s %s region %s" % (jet_algo, region['label'], angle_str)
                        jk_unfolder_plotter.draw_L_matrix_squared(title=title, **jk_plot_args)
                        title = "L * (x - bias vector)\n%s\n%s region\n%s" % (jet_algo, region['label'], angle_str)
                        jk_unfolder_plotter.draw_Lx_minus_bias(title=title, **jk_plot_args)

                    # Do unfolding!
                    # --------------------------------------------------------------
                    jk_unfolder.do_unfolding(jk_tau)
                    jk_unfolder.get_output(hist_name="%s_unfolded_1d" % jk_label_no_spaces)
                    jk_unfolder._post_process()
                    jk_unfolder.setup_normalised_results_per_pt_bin()

                    # Plot 1D results
                    # --------------------------------------------------------------
                    jk_title = "%s\n%s region, %s\n%s input" % (jet_algo, region['label'], angle_str, jk_label)
                    jk_unfolder_plotter.draw_unfolded_1d(title=jk_title, **jk_plot_args)

                    if REGULARIZE != "None":
                        # Draw our new template alongside unfolded
                        ocs = [
                            Contribution(alt_hist_mc_gen, label=region['alt_mc_label'],
                                         line_color=ROOT.kViolet+1,
                                         marker_color=ROOT.kViolet+1,
                                         subplot=jk_unfolder.hist_truth),
                            Contribution(jk_unfolder.truth_template, label="Template",
                                         line_color=ROOT.kAzure+1,
                                         marker_color=ROOT.kAzure+1,
                                         subplot=jk_unfolder.hist_truth),
                        ]
                        title = "%s\n%s region, %s" % (jet_algo, region['label'], angle_str)
                        jk_unfolder_plotter.draw_unfolded_1d(do_gen=True,
                                                             do_unfolded=True,
                                                             other_contributions=ocs,
                                                             output_dir=jk_output_dir,
                                                             append='jk_with_template',
                                                             title='',
                                                             subplot_title="* / Generator")

                        jk_unfolder_plotter.draw_bias_vector(title=jk_title, **jk_plot_args)
                        jk_unfolder_plotter.draw_x_minus_bias(title=jk_title, **jk_plot_args)
                        jk_unfolder_plotter.draw_Lx_minus_bias(title=jk_title, **jk_plot_args)

                    # save memory, in the Unfolder already
                    del region['jackknife_input_variations'][jk_ind]['input_reco']
                    del region['jackknife_input_variations'][jk_ind]['input_gen']

                    region['jackknife_input_variations'][jk_ind]['unfolder'] = jk_unfolder

                # ------------------------------------------------------------------
                # Update main input stat uncert with jackknife variations
                # ------------------------------------------------------------------
                unfolder.update_input_stat_uncert_from_jackknife(region['jackknife_input_variations'])
                unfolder.create_normalised_jackknife_input_uncertainty_per_pt_bin(region['jackknife_input_variations'])

                # ------------------------------------------------------------------
                # Big absolute plot with all jackknife variations
                # ------------------------------------------------------------------
                if len(region['jackknife_input_variations']) > 0:
                    # Do a big absolute 1D plots for sanity
                    # Compared to nominal
                    jk_contributions = [
                        Contribution(mdict['unfolder'].get_unfolded_with_ematrix_stat(),
                                     label=mdict['label'],
                                     line_color=mdict['colour'], line_style=1, line_width=1,
                                     marker_color=mdict['colour'], marker_size=0, marker_style=21,
                                     subplot=unfolder.get_unfolded_with_ematrix_stat())
                        for mdict in region['jackknife_input_variations']
                    ]
                    unfolder_plotter.draw_unfolded_1d(do_unfolded=True, do_gen=False,
                                                      output_dir=this_output_dir,
                                                      append='jackknife_input_%s' % append,
                                                      title=title,
                                                      other_contributions=jk_contributions,
                                                      subplot_title='#splitline{Variation /}{nominal}')

                    # Compared to respective truths
                    jk_contributions = [
                        Contribution(mdict['unfolder'].get_unfolded_with_ematrix_stat(),
                                     label=mdict['label'],
                                     line_color=mdict['colour'], line_style=1, line_width=1,
                                     marker_color=mdict['colour'], marker_size=0, marker_style=21,
                                     subplot=mdict['unfolder'].hist_truth)
                        for mdict in region['jackknife_input_variations']
                    ]
                    unfolder_plotter.draw_unfolded_1d(do_unfolded=True, do_gen=False,
                                                      output_dir=this_output_dir,
                                                      append='jackknife_input_vs_gen_%s' % append,
                                                      title=title,
                                                      other_contributions=jk_contributions,
                                                      subplot_title='#splitline{Variation /}{nominal}')

            if args.doJackknifeResponse:
                # first construct all new systematic variations dicts
                original_jk_dict = region['jackknife_response_variations'][0]
                original_jk_dict['label'] = '_jackknife_template'  # initial _ to ignore it later on
                tfile = original_jk_dict['tfile']
                if not isinstance(tfile, ROOT.TFile):
                    tfile = cu.open_root_file(tfile)

                region['jackknife_response_variations'] = []
                num_vars = len(original_jk_dict['variations'])
                for jk_ind in original_jk_dict['variations']:
                    region['jackknife_response_variations'].append(
                        {
                            "label": "Jackknife_response_%d" % jk_ind,
                            "response_map": cu.get_from_tfile(tfile, "%s/tu_%s_GenReco_all_jackknife_%d" % (region['dirname'], angle_shortname, jk_ind)),
                            "colour": cu.get_colour_seq(jk_ind, num_vars),
                        })

                # Now run over all variations, unfolding the nominal inputs
                # but with the various response matrices
                for jk_ind, jk_dict in enumerate(region['jackknife_response_variations']):
                    jk_label = jk_dict['label']
                    jk_label_no_spaces = cu.no_space_str(jk_label)

                    print("*" * 80)
                    print("*** Unfolding with jackknife response matrix:", jk_label, "(%d/%d) ***" % (jk_ind+1, len(region['jackknife_response_variations'])))
                    print("*" * 80)

                    jk_unfolder = MyUnfolder(response_map=rm_large_rel_error_bins_th2(jk_dict['response_map'], args.relErr),
                                             binning_handler=binning_handler,
                                             **unfolder_args)

                    jk_unfolder.SetEpsMatrix(eps_matrix)

                    jk_unfolder_plotter = MyUnfolderPlotter(jk_unfolder, is_data=not MC_INPUT)
                    jk_output_dir = os.path.join(this_output_dir, "jackknife_response", jk_label_no_spaces)
                    jk_plot_args = dict(output_dir=jk_output_dir, append=append)
                    cu.check_dir_exists_create(jk_output_dir)

                    # Set what is to be unfolded
                    # --------------------------------------------------------------
                    # Same input as nominal unfolder, since we only change responsematrix
                    jk_unfolder.set_input(**input_args())

                    # Subtract fakes (treat as background), same as nominal
                    # --------------------------------------------------------------
                    jk_unfolder.subtract_background(hist_fakes_reco, "Signal fakes", scale=1., scale_err=0.0)
                    jk_unfolder.subtract_background_gen_binning(hist_fakes_reco_gen_binning, "Signal fakes", scale=1.)

                    # Do regularisation
                    # --------------------------------------------------------------
                    # since the input is the same as the main unfolder's, we can
                    # use the same L matrix & bias vector
                    jk_tau = 0
                    if REGULARIZE != "None":
                        jk_unfolder.SetBias(unfolder.truth_template)
                        for L_args in unfolder.L_matrix_entries:
                            jk_unfolder.AddRegularisationCondition(*L_args)

                        if REGULARIZE == "L":
                            print("Regularizing with ScanLcurve, please be patient...")
                            jk_l_scanner = LCurveScanner()
                            jk_tau = jk_l_scanner.scan_L(tunfolder=jk_unfolder,
                                                         n_scan=args.nScan,
                                                         tau_min=region['tau_limits'][angle.var][0],
                                                         tau_max=region['tau_limits'][angle.var][1])
                            print("Found tau:", jk_tau)
                            jk_l_scanner.plot_scan_L_curve(output_filename="%s/scanL_%s.%s" % (jk_output_dir, append, OUTPUT_FMT))
                            jk_l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_%s.%s" % (jk_output_dir, append, OUTPUT_FMT))

                        elif REGULARIZE == "tau":
                            print("Regularizing with ScanTau, please be patient...")
                            jk_tau_scanner = TauScanner()
                            jk_tau = jk_tau_scanner.scan_tau(tunfolder=jk_unfolder,
                                                             n_scan=args.nScan,
                                                             tau_min=region['tau_limits'][angle.var][0],
                                                             tau_max=region['tau_limits'][angle.var][1],
                                                             scan_mode=scan_mode,
                                                             distribution=scan_distribution,
                                                             axis_steering=unfolder.axisSteering)
                            print("Found tau:", jk_tau)
                            jk_tau_scanner.plot_scan_tau(output_filename="%s/scantau_%s.%s" % (jk_output_dir, append, OUTPUT_FMT))

                    # Do unfolding!
                    # --------------------------------------------------------------
                    jk_unfolder.do_unfolding(jk_tau)
                    jk_unfolder.get_output(hist_name="%s_unfolded_1d" % jk_label_no_spaces)
                    jk_unfolder._post_process()
                    jk_unfolder.setup_normalised_results_per_pt_bin()

                    # Plot 1D results
                    # --------------------------------------------------------------
                    jk_title = "%s\n%s region, %s\n%s response matrix" % (jet_algo, region['label'], angle_str, jk_label)
                    jk_unfolder_plotter.draw_unfolded_1d(title=jk_title, **jk_plot_args)

                    del region['jackknife_response_variations'][jk_ind]['response_map']  # save memory

                    region['jackknife_response_variations'][jk_ind]['unfolder'] = jk_unfolder

                # ------------------------------------------------------------------
                # Update main response stat uncert with jackknife variations
                # ------------------------------------------------------------------
                unfolder.update_stat_response_from_jackknife(region['jackknife_response_variations'])
                unfolder.create_normalised_jackknife_response_uncertainty_per_pt_bin(region['jackknife_response_variations'])

                # ------------------------------------------------------------------
                # Big absolute plot with all jackknife variations
                # ------------------------------------------------------------------
                if len(region['jackknife_response_variations']) > 0:
                    # Do a big absolute 1D plots for sanity
                    jk_contributions = [
                        Contribution(mdict['unfolder'].get_unfolded_with_ematrix_stat(),
                                     label=mdict['label'],
                                     line_color=mdict['colour'], line_style=1, line_width=1,
                                     marker_color=mdict['colour'], marker_size=0, marker_style=21,
                                     subplot=unfolder.get_unfolded_with_ematrix_stat())
                        for mdict in region['jackknife_response_variations']
                    ]
                    unfolder_plotter.draw_unfolded_1d(do_unfolded=True, do_gen=False,
                                                      output_dir=this_output_dir,
                                                      append='jackknife_response_%s' % append,
                                                      title=title,
                                                      other_contributions=jk_contributions,
                                                      subplot_title='#splitline{Variation /}{nominal}')

            # Calculate exp systs manually, by unfolding with them separately
            # ------------------------------------------------------------------
            if args.doExperimentalSysts and args.doExperimentalSystsAsAltResponse:
                for syst_ind, syst_dict in enumerate(region['experimental_systematics']):
                    print("*******************************************************")
                    print("Adding systematic:", syst_dict['label'])
                    print("*******************************************************")
                    if not isinstance(syst_dict['tfile'], ROOT.TFile):
                        syst_dict['tfile'] = cu.open_root_file(syst_dict['tfile'])
                    map_syst = cu.get_from_tfile(syst_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                    for merge_func in th2_merge_funcs:
                        map_syst = merge_func(map_syst)

                    if args.zeroLastPtBin:
                        zero_last_pt_bin_th2_reco(map_syst)
                        zero_last_pt_bin_th2_gen(map_syst)

                    this_syst = ExpSystematic(label=syst_dict['label'], syst_map=None)

                    # construct unfolder like original but with this response matrix, do unfolding
                    exp_syst_unfolder = MyUnfolder(response_map=rm_large_rel_error_bins_th2(map_syst, args.relErr),
                                                   binning_handler=unfolder.binning_handler,
                                                   **unfolder_args)

                    exp_syst_unfolder.SetEpsMatrix(eps_matrix)

                    # Set what is to be unfolded - same as main unfolder
                    # --------------------------------------------------------------
                    exp_syst_unfolder.set_input(**input_args())

                    # Subtract fakes (treat as background), same as nominal
                    # --------------------------------------------------------------
                    exp_syst_unfolder.subtract_background(hist_fakes_reco, "fakes")
                    exp_syst_unfolder.subtract_background_gen_binning(hist_fakes_reco_gen_binning, "fakes")

                    # Setup regularisation
                    # --------------------------------------------------------------
                    exp_tau = 0
                    if REGULARIZE != "None":
                        exp_syst_unfolder.SetBias(unfolder.truth_template)
                        exp_syst_unfolder.setup_L_matrix_curvature(ref_hist=unfolder.truth_template, axis=args.regularizeAxis)

                        if REGULARIZE == "L":
                            print("Regularizing with ScanLcurve, please be patient...")
                            exp_l_scanner = LCurveScanner()
                            exp_tau = exp_l_scanner.scan_L(tunfolder=exp_syst_unfolder,
                                                           n_scan=args.nScan,
                                                           tau_min=region['tau_limits'][angle.var][0],
                                                           tau_max=region['tau_limits'][angle.var][1])
                            print("Found tau:", exp_tau)
                            exp_l_scanner.plot_scan_L_curve(output_filename="%s/scanL_%s_%s.%s" % (this_output_dir, append, this_syst.label_no_spaces, OUTPUT_FMT))
                            exp_l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_%s_%s.%s" % (this_output_dir, append, this_syst.label_no_spaces, OUTPUT_FMT))

                        elif REGULARIZE == "tau":
                            print("Regularizing with ScanTau, please be patient...")
                            exp_tau_scanner = TauScanner()
                            exp_tau = exp_tau_scanner.scan_tau(tunfolder=exp_syst_unfolder,
                                                               n_scan=args.nScan,
                                                               tau_min=region['tau_limits'][angle.var][0],
                                                               tau_max=region['tau_limits'][angle.var][1],
                                                               scan_mode=scan_mode,
                                                               distribution=scan_distribution,
                                                               axis_steering=exp_syst_unfolder.axisSteering)
                            print("Found tau:", exp_tau)
                            exp_tau_scanner.plot_scan_tau(output_filename="%s/scantau_%s_%s.%s" % (this_output_dir, append, this_syst.label_no_spaces, OUTPUT_FMT))

                    # Do unfolding!
                    # --------------------------------------------------------------
                    exp_syst_unfolder.do_unfolding(exp_tau)
                    exp_syst_unfolder.get_output(hist_name="exp_syst_%s_unfolded_1d" % this_syst.label_no_spaces)
                    exp_syst_unfolder._post_process()

                    # Do plots
                    # --------------------------------------------------------------
                    exp_syst_unfolder_plotter = MyUnfolderPlotter(exp_syst_unfolder, is_data=not MC_INPUT)
                    exp_syst_output_dir = this_output_dir+"/expSyst%s" % this_syst.label_no_spaces
                    exp_syst_plot_args = dict(output_dir=exp_syst_output_dir,
                                              append=append)
                    # if SUBTRACT_FAKES:
                    #     title = "%s\n%s region, %s, %s response map" % (jet_algo, region['label'], angle_str, this_syst.label)
                    #     exp_syst_unfolder_plotter.draw_detector_1d(do_reco_data_bg_sub=not MC_INPUT,
                    #                                           do_reco_bg=SUBTRACT_FAKES,
                    #                                           do_reco_mc_bg_sub=True,
                    #                                           output_dir=exp_syst_output_dir,
                    #                                           append='bg_fakes_subtracted_%s' % append,
                    #                                           title=title)

                    #     # same but with generator-binning
                    #     exp_syst_unfolder_plotter.draw_generator_1d(do_reco_data=False,
                    #                                            do_reco_data_bg_sub=not MC_INPUT,
                    #                                            do_reco_bg=True,
                    #                                            do_reco_mc=False,
                    #                                            do_reco_mc_bg_sub=True,
                    #                                            do_truth_mc=True,
                    #                                            output_dir=exp_syst_output_dir,
                    #                                            append='bg_fakes_subtracted_%s' % append,
                    #                                            title=title)

                    exp_syst_title = "%s\n%s region, %s, %s response map" % (jet_algo, region['label'], angle_str, this_syst.label)
                    exp_syst_unfolder_plotter.draw_unfolded_1d(title=exp_syst_title, **exp_syst_plot_args)

                    title = "Correlation matrix, %s, %s region, %s, %s response map" % (jet_algo, region['label'], angle_str, this_syst.label)
                    exp_syst_unfolder_plotter.draw_correlation_matrix(title=title, draw_values=False, **exp_syst_plot_args)

                    title = "Response matrix, %s, %s region, %s, %s" % (jet_algo, region['label'], angle_str, this_syst.label)
                    exp_syst_unfolder_plotter.draw_response_matrix(title=title, **exp_syst_plot_args)

                    title = "Probability matrix, %s, %s region, %s, %s" % (jet_algo, region['label'], angle_str, this_syst.label)
                    exp_syst_unfolder_plotter.draw_probability_matrix(title=title, **exp_syst_plot_args)

                    col_num = exp_syst_unfolder.binning_handler.physical_bin_to_global_bin(pt=151, var=0.6, binning_scheme='generator')
                    title = "Probability entries for gen bin %d, %s, %s region, %s, %s" % (col_num, jet_algo, region['label'], angle_str, this_syst.label)
                    exp_syst_unfolder_plotter.draw_probability_column(col_num, title=title, **exp_syst_plot_args)

                    exp_syst_unfolder.setup_normalised_results_per_pt_bin()

                    region['experimental_systematics'][syst_ind]['unfolder'] = exp_syst_unfolder

                    # Calculate shift in absolute due to systematic
                    # --------------------------------------------------------------
                    exp_syst_unfolded = exp_syst_unfolder.get_output()
                    nominal_unfolded = unfolder.get_output()
                    shift = exp_syst_unfolded.Clone()
                    shift.Add(nominal_unfolded, -1)
                    this_syst.syst_shifted = exp_syst_unfolded
                    this_syst.syst_shift = shift
                    this_syst.syst_ematrix = cu.shift_to_covariance(shift)
                    unfolder.exp_systs.append(this_syst)

            # Calculate experimental uncertainty shifts using results from another unfolding
            # ------------------------------------------------------------------
            if args.doExperimentalSystsFromFile is not None:
                print("Getting experimental systematics from another file...")
                angle_output_dir = "%s/%s/%s" % (args.doExperimentalSystsFromFile, region['name'], angle.var)
                this_pkl_filename = os.path.join(angle_output_dir, "unfolding_result.pkl")
                if not os.path.isfile(this_pkl_filename):
                    raise IOError("Cannot find systematics file, %s" % this_pkl_filename)

                exp_syst_region = unpickle_region(this_pkl_filename)

                reference_unfolder = exp_syst_region['unfolder']
                # ref_unfolded = reference_unfolder.unfolded

                # for exp_syst_dict in exp_syst_region.experimental_systematics:
                for exp_syst in reference_unfolder.exp_systs:
                    # For each systematic source, we figure out the
                    # relative shift compared to the original nominal result,
                    # for each normalised distribution (i.e. per pt bin).
                    # We then apply it to our new nominal result, and store it
                    # Note that this is different to taking the total shift,
                    # calculating its fractional diff, and then then applying it
                    # to the new nominal result?
                    # exp_syst = reference_unfolder.get_exp_syst(exp_syst_dict['label'])
                    this_exp_syst = ExpSystematic(label=exp_syst.label,
                                                  syst_map=None,
                                                  # copy the old shift/shifted, although no longer relevant here
                                                  syst_shift=exp_syst.syst_shift,
                                                  syst_shifted=exp_syst.syst_shifted)
                    unfolder.exp_systs.append(this_exp_syst)

                    # Add these to to HistBinChopper for later, but it isn't used
                    # Just to bypass internal checks that it exists in its cached objects
                    # when e.g. get_pt_bin_normed_div_bin_width() called
                    unfolder.hist_bin_chopper.add_obj(exp_syst.syst_shifted_label, unfolder.unfolded)
                    unfolder.hist_bin_chopper.add_obj(exp_syst.syst_shift_label, unfolder.unfolded)
                    unfolder.hist_bin_chopper.add_obj(exp_syst.syst_ematrix_label, unfolder.unfolded)

                    bins = unfolder.pt_bin_edges_gen
                    for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(bins[:-1], bins[1:])):
                        # get old nominal shape
                        hbc_args = dict(ind=ibin, binning_scheme='generator')
                        ref_nominal_hist = reference_unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width("unfolded", **hbc_args)

                        # calc rel shift
                        rel_syst_shift_hist = reference_unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width(exp_syst.syst_shift_label, **hbc_args).Clone()
                        rel_syst_shift_hist.Divide(ref_nominal_hist)

                        # get new nominal shape
                        new_syst_shift = unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width("unfolded", **hbc_args).Clone()

                        # apply rel shift to that to get shifted normalised hists
                        new_syst_shift.Multiply(rel_syst_shift_hist)
                        new_hist = unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width("unfolded", **hbc_args).Clone()
                        new_hist.Add(new_syst_shift)

                        # store shifted and shift hists
                        cu.remove_th1_errors(new_hist)
                        key = unfolder.hist_bin_chopper._generate_key(exp_syst.syst_shifted_label,
                                                                      ind=ibin,
                                                                      axis='pt',
                                                                      do_norm=True,
                                                                      do_div_bin_width=True,
                                                                      binning_scheme='generator')
                        unfolder.hist_bin_chopper._cache[key] = new_hist

                        cu.remove_th1_errors(new_syst_shift)
                        key = unfolder.hist_bin_chopper._generate_key(exp_syst.syst_shift_label,
                                                                      ind=ibin,
                                                                      axis='pt',
                                                                      do_norm=True,
                                                                      do_div_bin_width=True,
                                                                      binning_scheme='generator')
                        unfolder.hist_bin_chopper._cache[key] = new_syst_shift

                        # calculate error matrix
                        syst_ematrix = cu.shift_to_covariance(new_syst_shift)
                        key = unfolder.hist_bin_chopper._generate_key(exp_syst.syst_ematrix_label,
                                                                      ind=ibin,
                                                                      axis='pt',
                                                                      do_norm=True,
                                                                      do_div_bin_width=True,
                                                                      binning_scheme='generator')
                        unfolder.hist_bin_chopper._cache[key] = syst_ematrix

                # update region info
                region['experimental_systematics'] = [syst_dict for syst_dict in exp_syst_region['experimental_systematics']
                                                      if syst_dict['label'] in reference_unfolder.get_all_exp_syst_labels()]

                # cleanup, save memory
                del exp_syst_region

            else:
                if not args.jacobian:
                    unfolder.setup_normalised_experimental_systs_per_pt_bin()

            prof_end_exp_systs()

            # Draw big 1D distributions
            # ------------------------------------------------------------------
            title = "%s\n%s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_unfolded_1d(output_dir=this_output_dir,
                                              append=append,
                                              title=title,
                                              do_gen=True,
                                              do_logy=True)

            # reco using detector binning
            unfolder_plotter.draw_detector_1d(do_reco_mc=True,
                                              do_reco_data=not MC_INPUT,
                                              output_dir=this_output_dir,
                                              append=append,
                                              title=title)

            # reco using gen binning
            unfolder_plotter.draw_generator_1d(do_reco_data=not MC_INPUT,
                                               do_reco_data_bg_sub=False,
                                               do_reco_bg=False,
                                               do_reco_mc=True,
                                               do_reco_mc_bg_sub=False,
                                               do_truth_mc=True,
                                               output_dir=this_output_dir,
                                               append=append,
                                               title=title)

            # same plot but with bg-subtracted reco
            unfolder_plotter.draw_detector_1d(do_reco_data_bg_sub=not MC_INPUT,
                                              do_reco_bg=True,
                                              do_reco_mc_bg_sub=True,
                                              output_dir=this_output_dir,
                                              append='bg_fakes_subtracted_%s' % append,
                                              title=title)

            # same but with generator-binning
            unfolder_plotter.draw_generator_1d(do_reco_data=False,
                                               do_reco_data_bg_sub=not MC_INPUT,
                                               do_reco_bg=True,
                                               do_reco_mc=False,
                                               do_reco_mc_bg_sub=True,
                                               do_truth_mc=True,
                                               output_dir=this_output_dir,
                                               append='bg_fakes_subtracted_%s' % append,
                                               title=title)

            # Draw projections of response matrix vs 1D hist to check normalisation OK
            # Only makes sense if the same MC events go into matrix & 1D plot
            # ------------------------------------------------------------------
            # FIXME should work for mergeBins
            if not MC_SPLIT and args.relErr < 0 and not args.mergeBins and not args.mergeLastPtBin:
                # on gen axis
                proj_gen = unfolder.response_map.ProjectionX("proj_gen_%s" % (append))
                draw_projection_comparison(unfolder.hist_truth, proj_gen,
                                           title="%s\n%s region" % (jet_algo, region['label']),
                                           xtitle="%s, Generator binning" % (angle_str),
                                           output_filename="%s/projection_gen_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

                # on detector axis
                proj_reco = unfolder.response_map.ProjectionY("proj_reco_%s" % (append))
                draw_projection_comparison(unfolder.hist_mc_reco_bg_subtracted, proj_reco,
                                           title="%s\n%s region" % (jet_algo, region['label']),
                                           xtitle="%s, Detector binning" % (angle_str),
                                           output_filename="%s/projection_reco_bg_subtracted_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

            # Draw big 1D plot of unfolded shifted exp systs
            # ------------------------------------------------------------------
            if len(region['experimental_systematics']) > 0 and MC_INPUT:
                # Do a big absolute 1D plot for sanity
                syst_contributions = [
                    Contribution(unfolder.get_syst_shifted_hist(sdict['label'], unfolder.get_unfolded_with_ematrix_stat()),
                                 label=sdict['label'],
                                 line_color=sdict['colour'], line_style=2 if 'down' in sdict['label'].lower() else 1, line_width=1,
                                 marker_color=sdict['colour'], marker_size=0, marker_style=21,
                                 subplot=unfolder.get_unfolded_with_ematrix_stat())
                    for sdict in region['experimental_systematics']
                ]
                unfolder_plotter.draw_unfolded_1d(do_unfolded=True, do_gen=False,
                                                  output_dir=this_output_dir,
                                                  append='exp_systs_%s' % append,
                                                  title=title,
                                                  other_contributions=syst_contributions,
                                                  subplot_title='#splitline{Variation /}{nominal}')

            # Draw matrices
            # ------------------------------------------------------------------
            if REGULARIZE != "None":
                unfolder_plotter.draw_bias_vector(title=title, **plot_args)
                unfolder_plotter.draw_x_minus_bias(title=title, **plot_args)

            title = "Response matrix, %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_response_matrix(title=title, **plot_args)

            title = "Response matrix normed by detector p_{T} bin, %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_response_matrix_normed_by_detector_pt(title=title, **plot_args)

            title = ("#splitline{Probability matrix, %s, %s region, %s}{Condition number: #sigma_{max} / #sigma_{min} = %.3g / %.3g = %g}"
                     % (jet_algo, region['label'], angle_str, unfolder.sigma_max, unfolder.sigma_min, unfolder.condition_number))
            unfolder_plotter.draw_probability_matrix(title=title, **plot_args)

            col_num = unfolder.binning_handler.physical_bin_to_global_bin(pt=151, var=0.6, binning_scheme='generator')
            title = "Probability entries for gen bin %d, %s, %s region, %s" % (col_num, jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_probability_column(col_num, title=title, **plot_args)

            title = "%s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_failed_reco(title=title, **plot_args)

            title = "Correlation matrix, %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_correlation_matrix(title=title, draw_values=True, **plot_args)
            unfolder_plotter.draw_correlation_matrix(title=title, draw_values=False, **plot_args)

            title = "Error matrix (input), %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_error_matrix_input(title=title, **plot_args)

            title = "Error matrix (statistical input + backgrounds), %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_error_matrix_stat(title=title, **plot_args)

            title = "Error matrix (stat. response matrix), %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_error_matrix_stat_response(title=title, **plot_args)

            if args.doExperimentalSysts:
                for syst_dict in region['experimental_systematics']:
                    title = "Error matrix (%s systematic), %s, %s region, %s" % (syst_dict['label'], jet_algo, region['label'], angle_str)
                    unfolder_plotter.draw_error_matrix_syst(syst_dict['label'], title=title, **plot_args)

            if REGULARIZE != "None":
                title = "Error matrix (tau), %s, %s region, %s" % (jet_algo, region['label'], angle_str)
                unfolder_plotter.draw_error_matrix_tau(title=title, **plot_args)

            title = "Error matrix (total), %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_error_matrix_tunfold_total(title=title, **plot_args)

            # Do forward-folding to check unfolding
            # ------------------------------------------------------------------
            # Do it on the unfolded result
            title = "%s\n%s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_truth_unfolded_folded(draw_truth_folded=False, title=title, **plot_args)

            # Do some bottom-line tests
            # ------------------------------------------------------------------
            # smeared_chi2, smeared_ndf, smeared_p = unfolder.calculate_chi2(one_hist=unfolder.get_folded_mc_truth(),
            #                                                                other_hist=unfolder.hist_mc_reco_bg_subtracted,
            #                                                                cov_inv_matrix=unfolder.get_vyy_inv_ndarray(),
            #                                                                # cov_inv_matrix=unfolder.get_vyy_inv_no_bg_subtraction_ndarray(),
            #                                                                detector_space=True,
            #                                                                ignore_underflow_bins=True,
            #                                                                debugging_dir=os.path.join(this_output_dir, 'smeared_chi2_debug'))
            # print('smeared chi2:', smeared_chi2, smeared_ndf, smeared_chi2/smeared_ndf, smeared_p)

            # unfolded_chi2, unfolded_ndf, unfolded_p = unfolder.calculate_chi2(one_hist=unfolder.unfolded,
            #                                                                   other_hist=unfolder.hist_truth,
            #                                                                   cov_inv_matrix=unfolder.get_vxx_inv_ndarray(),
            #                                                                   detector_space=False,
            #                                                                   ignore_underflow_bins=True,
            #                                                                   debugging_dir=os.path.join(this_output_dir, 'unfolded_chi2_debug'))
            # print('unfolded chi2:', unfolded_chi2, unfolded_ndf, unfolded_chi2/unfolded_ndf, unfolded_p)

            # ------------------------------------------------------------------
            # UNFOLDING WITH ALTERNATIVE RESPONSE MATRIX
            # ------------------------------------------------------------------
            if args.useAltResponse:
                print("*" * 80)
                print("*** Unfolding with alternate response matrix ***")
                print("*" * 80)

                hist_mc_gen_reco_map_alt = cu.get_from_tfile(region['alt_mc_tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                hist_mc_gen_reco_map_alt.Scale(unfolder.response_map.Integral() / hist_mc_gen_reco_map_alt.Integral())  # just for display purposes, doesn't affect result
                for merge_func in th2_merge_funcs:
                        hist_mc_gen_reco_map_alt = merge_func(hist_mc_gen_reco_map_alt)

                alt_unfolder = MyUnfolder(response_map=rm_large_rel_error_bins_th2(hist_mc_gen_reco_map_alt, args.relErr),
                                          binning_handler=unfolder.binning_handler,
                                          **unfolder_args)

                # SetEpsMatrix ensures rank properly calculated when inverting
                # Needed if you get message "rank of matrix E 55 expect 170"
                # And unfolded looks wacko
                alt_unfolder.SetEpsMatrix(eps_matrix)

                alt_unfolder_plotter = MyUnfolderPlotter(alt_unfolder, is_data=not MC_INPUT)
                alt_output_dir = this_output_dir+"/altResponse"
                alt_plot_args = dict(output_dir=alt_output_dir,
                                     append=append)

                # Only plot response matrix
                # --------------------------------------------------------------
                title = ("#splitline{Probability matrix, %s region, %s, %s}{Condition number: #sigma_{max} / #sigma_{min} = %.3g / %.3g = %g}"
                            % (region['label'], angle_str, region['alt_mc_label'], unfolder.sigma_max, unfolder.sigma_min, unfolder.condition_number))
                alt_unfolder_plotter.draw_probability_matrix(title=title, **alt_plot_args)

                title = "Response matrix, %s, %s region, %s, %s" % (jet_algo, region['label'], angle_str, region['alt_mc_label'])
                alt_unfolder_plotter.draw_response_matrix(title=title, **alt_plot_args)

                alt_title = "%s\n%s region, %s, %s response map" % (jet_algo, region['label'], angle_str, region['alt_mc_label'])
                alt_unfolder_plotter.draw_failed_reco(title=alt_title, **alt_plot_args)

                # Compare with the nominal one
                ocs = [Contribution(unfolder.get_failed_reco(as_fraction=True),
                                    label=region['mc_label'],
                                    line_color=ROOT.kBlack,
                                    subplot=alt_unfolder.get_failed_reco(as_fraction=True)
                                    )]
                alt_unfolder_plotter.draw_failed_reco(title=alt_title,
                                                      output_dir=alt_output_dir,
                                                      other_contributions=ocs,
                                                      append="compare_nominal_"+append)

                # Set what is to be unfolded - same as main unfolder
                # --------------------------------------------------------------
                alt_unfolder.set_input(**input_args())

                # Subtract fakes (treat as background), same as nominal
                # --------------------------------------------------------------
                alt_unfolder.subtract_background(hist_fakes_reco, "fakes")
                alt_unfolder.subtract_background_gen_binning(hist_fakes_reco_gen_binning, "fakes")

                # Do regularisation
                # --------------------------------------------------------------
                # since the input is the same as the main unfolder's, we can
                # use the same L matrix & bias vector
                alt_tau = 0
                if REGULARIZE != "None":
                    alt_unfolder.SetBias(unfolder.truth_template)
                    for L_args in unfolder.L_matrix_entries:
                        alt_unfolder.AddRegularisationCondition(*L_args)

                    if REGULARIZE == "L":
                        print("Regularizing with ScanLcurve, please be patient...")
                        alt_l_scanner = LCurveScanner()
                        alt_tau = alt_l_scanner.scan_L(tunfolder=alt_unfolder,
                                                       n_scan=args.nScan,
                                                       tau_min=region['tau_limits'][angle.var][0],
                                                       tau_max=region['tau_limits'][angle.var][1])
                        print("Found tau:", alt_tau)
                        alt_l_scanner.plot_scan_L_curve(output_filename="%s/scanL_%s.%s" % (alt_output_dir, append, OUTPUT_FMT))
                        alt_l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_%s.%s" % (alt_output_dir, append, OUTPUT_FMT))

                    elif REGULARIZE == "tau":
                        print("Regularizing with ScanTau, please be patient...")
                        alt_tau_scanner = TauScanner()
                        alt_tau = alt_tau_scanner.scan_tau(tunfolder=alt_unfolder,
                                                           n_scan=args.nScan,
                                                           tau_min=region['tau_limits'][angle.var][0],
                                                           tau_max=region['tau_limits'][angle.var][1],
                                                           scan_mode=scan_mode,
                                                           distribution=scan_distribution,
                                                           axis_steering=unfolder.axisSteering)
                        print("Found tau:", alt_tau)
                        alt_tau_scanner.plot_scan_tau(output_filename="%s/scantau_%s.%s" % (alt_output_dir, append, OUTPUT_FMT))

                # Do unfolding!
                # --------------------------------------------------------------
                alt_unfolder.do_unfolding(alt_tau)
                alt_unfolder.get_output(hist_name="alt_unfolded_1d")
                alt_unfolder._post_process()

                # Draw 1D & 2D plots
                # --------------------------------------------------------------
                alt_unfolder_plotter.draw_detector_1d(do_reco_data_bg_sub=not MC_INPUT,
                                                      do_reco_bg=True,
                                                      do_reco_mc_bg_sub=True,
                                                      output_dir=alt_output_dir,
                                                      append='bg_fakes_subtracted_%s' % append,
                                                      title=title)

                # same but with generator-binning
                alt_unfolder_plotter.draw_generator_1d(do_reco_data=False,
                                                       do_reco_data_bg_sub=not MC_INPUT,
                                                       do_reco_bg=True,
                                                       do_reco_mc=False,
                                                       do_reco_mc_bg_sub=True,
                                                       do_truth_mc=True,
                                                       output_dir=alt_output_dir,
                                                       append='bg_fakes_subtracted_%s' % append,
                                                       title=title)

                alt_unfolder_plotter.draw_unfolded_1d(title=alt_title, **alt_plot_args)

                title = "Correlation matrix, %s, %s region, %s, %s response map" % (jet_algo, region['label'], angle_str, region['alt_mc_label'])
                alt_unfolder_plotter.draw_correlation_matrix(title=title, draw_values=False, **alt_plot_args)

                title = "Error matrix (statistical input + backgrounds), %s, %s region, %s, %s response map" % (jet_algo, region['label'], angle_str, region['alt_mc_label'])
                alt_unfolder_plotter.draw_error_matrix_stat(title=title, **alt_plot_args)

                title = "Error matrix (total), %s, %s region, %s, %s response map" % (jet_algo, region['label'], angle_str, region['alt_mc_label'])
                alt_unfolder_plotter.draw_error_matrix_tunfold_total(title=title, **alt_plot_args)

                # Do some more bottom-line tests:
                # --------------------------------------------------------------
                # Do before saving to file otherwise objects get deleted
                # if alt_unfolder:
                #     smeared_alt_chi2, smeared_alt_ndf, smeared_alt_p = unfolder.calculate_chi2(one_hist=alt_unfolder.get_folded_mc_truth(),
                #                                                                                other_hist=unfolder.hist_mc_reco_bg_subtracted,
                #                                                                                cov_inv_matrix=unfolder.get_vyy_inv_ndarray(),
                #                                                                                # cov_inv_matrix=unfolder.get_vyy_inv_no_bg_subtraction_ndarray(),
                #                                                                                detector_space=True,
                #                                                                                ignore_underflow_bins=True,
                #                                                                                debugging_dir=os.path.join(this_output_dir, 'smeared_alt_chi2_debug'))
                #     print('smeared chi2 (alt MC):', smeared_alt_chi2, smeared_alt_ndf, smeared_alt_chi2/smeared_alt_ndf, smeared_alt_p)

                # print(unfolder.unfolded.Integral())
                # print(alt_hist_mc_gen.Integral())
                # if alt_hist_mc_gen:
                #     unfolded_alt_chi2, unfolded_alt_ndf, unfolded_alt_p = unfolder.calculate_chi2(one_hist=unfolder.unfolded,
                #                                                                                   other_hist=alt_hist_mc_gen,
                #                                                                                   cov_inv_matrix=unfolder.get_vxx_inv_ndarray(),
                #                                                                                   detector_space=False,
                #                                                                                   ignore_underflow_bins=True,
                #                                                                                   debugging_dir=os.path.join(this_output_dir, 'unfolded_alt_chi2_debug'))
                #     print('unfolded chi2 (alt MC):', unfolded_alt_chi2, unfolded_alt_ndf, unfolded_alt_chi2/unfolded_alt_ndf, unfolded_alt_p)

                alt_unfolder.setup_normalised_results_per_pt_bin()

                region['alt_unfolder'] = alt_unfolder

            # Bit gnarly - have to save this stuff manually
            region["alt_hist_mc_gen"] = alt_hist_mc_gen
            region["alt_hist_mc_reco"] = alt_hist_mc_reco
            region["alt_hist_mc_reco_bg_subtracted"] = alt_hist_mc_reco_bg_subtracted
            region["alt_hist_mc_reco_bg_subtracted_gen_binning"] = alt_hist_mc_reco_bg_subtracted_gen_binning

            # ------------------------------------------------------------------
            # SCALE SYST VARIATIONS
            # ------------------------------------------------------------------
            # For each scale variation, we unfold the nominal input using the
            # variation's response matrix. Then we take the envelope as the final
            # scale uncertainty
            if args.doScaleSysts:
                for ind, scale_dict in enumerate(region['scale_systematics']):
                    scale_label = scale_dict['label']
                    scale_label_no_spaces = cu.no_space_str(scale_dict['label'])

                    print("*" * 80)
                    print("*** Unfolding scale variation:", scale_label, "(%d/%d) ***" % (ind+1, len(region['scale_systematics'])))
                    print("*" * 80)

                    if isinstance(scale_dict['tfile'], str):
                        scale_dict['tfile'] = cu.open_root_file(scale_dict['tfile'])
                    scale_map = cu.get_from_tfile(scale_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                    for merge_func in th2_merge_funcs:
                        scale_map = merge_func(scale_map)
                    scale_unfolder = MyUnfolder(response_map=rm_large_rel_error_bins_th2(scale_map, args.relErr),
                                                binning_handler=unfolder.binning_handler,
                                                **unfolder_args)

                    # Needed because some of the variations struggle to unfold
                    # Even 1E-18 wouldn't work - needs to be v.small
                    scale_unfolder.SetEpsMatrix(eps_matrix)

                    scale_unfolder_plotter = MyUnfolderPlotter(scale_unfolder, is_data=not MC_INPUT)
                    scale_output_dir = this_output_dir+"/scaleSyst/"+scale_label_no_spaces
                    scale_plot_args = dict(output_dir=scale_output_dir,
                                           append=append)

                    # Set what is to be unfolded
                    # --------------------------------------------------------------
                    # Same input as nominal unfolder, since we only change response matrix
                    scale_unfolder.set_input(**input_args())

                    # Subtract fakes (treat as background), same as nominal
                    # --------------------------------------------------------------
                    scale_unfolder.subtract_background(hist_fakes_reco, "Signal fakes", scale=1., scale_err=0.0)
                    scale_unfolder.subtract_background_gen_binning(hist_fakes_reco_gen_binning, "Signal fakes", scale=1.)

                    # Do regularisation
                    # --------------------------------------------------------------
                    # since the input is the same as the main unfolder's, we can
                    # use the same L matrix & bias vector
                    scale_tau = 0
                    if REGULARIZE != "None":
                        scale_unfolder.SetBias(unfolder.truth_template)
                        for L_args in unfolder.L_matrix_entries:
                            scale_unfolder.AddRegularisationCondition(*L_args)

                        if REGULARIZE == "L":
                            print("Regularizing with ScanLcurve, please be patient...")
                            scale_l_scanner = LCurveScanner()
                            scale_tau = scale_l_scanner.scan_L(tunfolder=scale_unfolder,
                                                               n_scan=args.nScan,
                                                               tau_min=region['tau_limits'][angle.var][0],
                                                               tau_max=region['tau_limits'][angle.var][1])
                            print("Found tau:", scale_tau)
                            scale_l_scanner.plot_scan_L_curve(output_filename="%s/scanL_%s.%s" % (scale_output_dir, append, OUTPUT_FMT))
                            scale_l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_%s.%s" % (scale_output_dir, append, OUTPUT_FMT))

                        elif REGULARIZE == "tau":
                            print("Regularizing with ScanTau, please be patient...")
                            scale_tau_scanner = TauScanner()
                            scale_tau = scale_tau_scanner.scan_tau(tunfolder=scale_unfolder,
                                                                   n_scan=args.nScan,
                                                                   tau_min=region['tau_limits'][angle.var][0],
                                                                   tau_max=region['tau_limits'][angle.var][1],
                                                                   scan_mode=scan_mode,
                                                                   distribution=scan_distribution,
                                                                   axis_steering=unfolder.axisSteering)
                            print("Found tau:", scale_tau)
                            scale_tau_scanner.plot_scan_tau(output_filename="%s/scantau_%s.%s" % (scale_output_dir, append, OUTPUT_FMT))

                    scale_dict['tau'] = scale_tau

                    # Do unfolding!
                    # --------------------------------------------------------------
                    scale_unfolder.do_unfolding(scale_tau)
                    scale_unfolder.get_output(hist_name="scale_%s_unfolded_1d" % scale_label_no_spaces)
                    scale_unfolder._post_process()
                    if not args.jacobian:
                        scale_unfolder.setup_normalised_results_per_pt_bin()

                    # Plot absolute 1D result
                    scale_title = "%s\n%s region, %s, %s input" % (jet_algo, region['label'], angle_str, scale_label)
                    scale_unfolder_plotter.draw_unfolded_1d(title=scale_title, **scale_plot_args)
                    scale_dict['unfolder'] = scale_unfolder

                    scale_unfolder.slim_down(keep_1d_hists=True)
                    cu.close_tfile(scale_dict['tfile'])

                    scale_dict['unfolder'] = scale_unfolder

                if not args.jacobian:
                    unfolder.create_normalised_scale_syst_uncertainty_per_pt_bin(region['scale_systematics'])
                    unfolder.create_normalised_scale_syst_ematrices_per_pt_bin()

                    unfolder.create_scale_syst_uncertainty_per_pt_bin(region['scale_systematics'])

                    unfolder.create_scale_syst_uncertainty_per_lambda_bin(region['scale_systematics'])

            # ------------------------------------------------------------------
            # LOAD SCALE VARIATIONS FROM ANOTHER FILE
            # ------------------------------------------------------------------
            # Load scale systs from another reference file, and calc fractional
            # uncertainty, apply to this unfolded result
            if args.doScaleSystsFromFile is not None:
                scale_dir = "%s/%s/%s" % (args.doScaleSystsFromFile, region['name'], angle.var)

                this_pkl_filename = os.path.join(scale_dir, "unfolding_result.pkl")
                if not os.path.isfile(this_pkl_filename):
                    raise IOError("Cannot find scale systematics file, %s" % this_pkl_filename)
                scale_syst_region = unpickle_region(this_pkl_filename)

                # update original region object with the scale syst info from the reference file
                region['scale_systematics'] = scale_syst_region['scale_systematics']

                # Add dummy object to hist_bin_chopper for later, so we can directly manipulate the cache
                # this replaces the procedure in create_normalised_scale_syst_uncertainty_per_pt_bin()
                unfolder.hist_bin_chopper.add_obj(unfolder.scale_uncert_name, unfolder.get_unfolded_with_ematrix_stat())
                unfolder.hist_bin_chopper.add_obj("unfolded_stat_err", unfolder.get_unfolded_with_ematrix_stat())

                for ibin_pt in range(len(unfolder.pt_bin_edges_gen[:-1])):
                    # Calculate scale uncertainty by taking relative uncertainty
                    # from reference file (which was already calculated in the 1st running),
                    # and applying to this result
                    key = unfolder.hist_bin_chopper._generate_key(unfolder.scale_uncert_name,
                                                                  ind=ibin_pt,
                                                                  axis='pt',
                                                                  do_norm=True,
                                                                  do_div_bin_width=True,
                                                                  binning_scheme='generator')
                    ref_scale_syst = scale_syst_region['unfolder'].hist_bin_chopper._cache[key]
                    scale_syst = unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width("unfolded_stat_err", ibin_pt, binning_scheme='generator').Clone()
                    for i in range(1, scale_syst.GetNbinsX()+1):
                        if ref_scale_syst.GetBinContent(i) != 0:
                            rel_err = ref_scale_syst.GetBinError(i) / ref_scale_syst.GetBinContent(i)
                        else:
                            rel_err = 0
                        scale_syst.SetBinError(i, rel_err * scale_syst.GetBinContent(i))
                    unfolder.hist_bin_chopper._cache[key] = scale_syst

                unfolder.create_normalised_scale_syst_ematrices_per_pt_bin()

                # cleanup, save memory
                del scale_syst_region

            # ------------------------------------------------------------------
            # BIG ABSOLUTE PLOT WITH ALL SCALE VARIATIONS
            # ------------------------------------------------------------------
            if len(region['scale_systematics']) > 0:
                # Do a big absolute 1D plots for sanity
                scale_contributions = [
                    Contribution(mdict['unfolder'].get_unfolded_with_ematrix_stat(),
                                 label=mdict['label'],
                                 line_color=mdict['colour'], line_style=1, line_width=1,
                                 marker_color=mdict['colour'], marker_size=0, marker_style=21,
                                 subplot=unfolder.get_unfolded_with_ematrix_stat())
                    for mdict in region['scale_systematics']
                ]
                unfolder_plotter.draw_unfolded_1d(do_unfolded=True, do_gen=False,
                                                  output_dir=this_output_dir,
                                                  append='scale_systs_%s' % append,
                                                  title=title,
                                                  other_contributions=scale_contributions,
                                                  subplot_title='#splitline{Variation /}{nominal}')

            prof_end_scale_systs()

            # ------------------------------------------------------------------
            # MODEL INPUT VARIATIONS
            # ------------------------------------------------------------------
            # For each model variation, we unfold using the same settings as
            # the nominal one, just changing the input 1D hist
            if args.doModelSysts:
                syst_entries = []
                for ind, syst_dict in enumerate(region['model_systematics']):
                    syst_label = syst_dict['label']
                    syst_label_no_spaces = cu.no_space_str(syst_dict['label'])

                    print("*" * 80)
                    print("*** Unfolding with model syst input:", syst_label, "(%d/%d) ***" % (ind+1, len(region['model_systematics'])))
                    print("*" * 80)

                    # is_herwig = "Herwig" in syst_label
                    mc_hname_append = "all"

                    if not isinstance(syst_dict['tfile'], ROOT.TFile):
                        syst_dict['tfile'] = cu.open_root_file(syst_dict['tfile'])

                    hist_syst_reco = cu.get_from_tfile(syst_dict['tfile'], "%s/hist_%s_reco_%s" % (region['dirname'], angle_shortname, mc_hname_append))
                    hist_syst_reco_gen_binning = cu.get_from_tfile(syst_dict['tfile'], "%s/hist_%s_reco_gen_binning" % (region['dirname'], angle_shortname))
                    hist_syst_gen = cu.get_from_tfile(syst_dict['tfile'], "%s/hist_%s_truth_%s" % (region['dirname'], angle_shortname, mc_hname_append))

                    syst_unfolder = MyUnfolder(response_map=unfolder.response_map,
                                               binning_handler=unfolder.binning_handler,
                                               **unfolder_args)

                    syst_output_dir = this_output_dir+"/modelSyst_"+syst_label_no_spaces
                    syst_unfolder_plotter = MyUnfolderPlotter(syst_unfolder, is_data=False)
                    syst_plot_args = dict(output_dir=syst_output_dir, append=append)

                    # SetEpsMatrix ensures rank properly calculated when inverting
                    # Needed if you get message "rank of matrix E 55 expect 170"
                    # And unfolded looks wacko
                    syst_unfolder.SetEpsMatrix(eps_matrix)

                    # because we only care about shape, not overall normalisation
                    # (which can artificially increase/decrease errors)
                    # we normalise to the nominal integral
                    # Note that we use the scaling from gen level, to take
                    # into account any reco-dependent efficiencies
                    # TODO: is this right?
                    sf = hist_mc_gen.Integral() / hist_syst_gen.Integral()
                    hist_syst_reco = rm_large_rel_error_bins_th1(hist_syst_reco, args.relErr)
                    hist_syst_reco.Scale(sf)
                    hist_syst_gen.Scale(sf)
                    hist_syst_reco_gen_binning.Scale(sf)

                    for merge_func in th1_merge_reco_funcs:
                        hist_syst_reco = merge_func(hist_syst_reco)
                    for merge_func in th1_merge_gen_funcs:
                        hist_syst_gen = merge_func(hist_syst_gen)
                        hist_syst_reco_gen_binning = merge_func(hist_syst_reco_gen_binning)

                    # Set what is to be unfolded
                    # --------------------------------------------------------------
                    # Use the background template from the nominal MC
                    # (since we're only testing different input shapes,
                    # and our bkg estimate is always from MC)
                    input_handler_syst = InputHandler(input_hist=hist_syst_reco,
                                                      hist_truth=hist_syst_gen,
                                                      hist_mc_reco=hist_syst_reco,
                                                      hist_mc_fakes=hist_mc_fakes_reco)

                    input_handler_gen_binning_syst = InputHandler(input_hist=hist_syst_reco_gen_binning,
                                                                  hist_truth=None,
                                                                  hist_mc_reco=hist_syst_reco_gen_binning,
                                                                  hist_mc_fakes=hist_mc_fakes_reco_gen_binning)

                    # Set what is to be unfolded
                    # ------------------------------------------------------------------
                    syst_unfolder.set_input(input_handler=input_handler_syst,
                                            input_handler_gen_binning=input_handler_gen_binning_syst,
                                            bias_factor=args.biasFactor,
                                            error_unconstrained_truth_bins=False)

                    # Subtract fakes (treat as background)
                    # --------------------------------------------------------------
                    syst_unfolder.subtract_background(input_handler_syst.calc_fake_hist(hist_syst_reco), "fakes")
                    syst_unfolder.subtract_background_gen_binning(input_handler_gen_binning_syst.calc_fake_hist(hist_syst_reco_gen_binning), "fakes")

                    # also show nominal bg-subtracted input for comparison
                    ocs = [
                        Contribution(unfolder.input_hist_bg_subtracted,
                                     label='Nominal unfolding input (bg-subtracted)',
                                     line_color=ROOT.kRed, line_width=1)
                    ]
                    syst_unfolder_plotter.draw_detector_1d(do_reco_data=False,
                                                           do_reco_data_bg_sub=False,
                                                           do_reco_bg=True,
                                                           do_reco_mc=False,
                                                           do_reco_mc_bg_sub=True,
                                                           other_contributions=ocs,
                                                           output_dir=syst_plot_args['output_dir'],
                                                           append='bg_fakes_subtracted_%s' % append,
                                                           title="%s region, %s, %s" % (region['label'], angle_str, syst_label))

                    # Do regularisation
                    # --------------------------------------------------------------
                    # since the input isn't the same as the main unfolder's,
                    # we can have to create new template, L matrix, etc
                    syst_tau = 0
                    if REGULARIZE != "None":

                        # Create truth template by fitting MC to data @ detector level
                        # --------------------------------------------------------------
                        # Fit the two MC templates to data to get their fractions
                        # Then use the same at truth level
                        # This template will allow us to setup a more accurate L matrix,
                        # and a bias hist
                        syst_template_maker = TruthTemplateMaker(binning_handler=unfolder.binning_handler,
                                                                 output_dir=syst_output_dir)

                        syst_template_maker.set_input(syst_unfolder.input_hist_gen_binning_bg_subtracted)

                        syst_template_maker.add_mc_template(name=region['mc_label'],
                                                            hist_reco=input_handler_gen_binning_syst.hist_mc_reco_bg_subtracted,
                                                            hist_gen=hist_mc_gen_all,
                                                            colour=ROOT.kRed)
                        syst_template_maker.add_mc_template(name=region['alt_mc_label'],
                                                            hist_reco=alt_hist_mc_reco_bg_subtracted_gen_binning,
                                                            hist_gen=alt_hist_mc_gen,
                                                            colour=ROOT.kViolet+1)

                        syst_truth_template = syst_template_maker.create_template()
                        syst_unfolder.truth_template = syst_truth_template
                        syst_unfolder.hist_bin_chopper.add_obj("truth_template", syst_unfolder.truth_template)
                        syst_unfolder.hist_bin_chopper.add_obj("alt_hist_truth", alt_hist_mc_gen)

                        if MC_INPUT:
                            # do check by fitting at gen level and comparing
                            syst_template_maker.set_input_gen(syst_unfolder.hist_truth)
                            syst_template_maker.check_template_at_gen()

                        syst_unfolder.SetBias(syst_unfolder.truth_template)
                        syst_unfolder.setup_L_matrix_curvature(ref_hist=syst_unfolder.truth_template, axis=args.regularizeAxis)

                        if REGULARIZE == "L":
                            print("Regularizing with ScanLcurve, please be patient...")
                            syst_l_scanner = LCurveScanner()
                            syst_tau = syst_l_scanner.scan_L(tunfolder=syst_unfolder,
                                                             n_scan=args.nScan,
                                                             tau_min=region['tau_limits'][angle.var][0],
                                                             tau_max=region['tau_limits'][angle.var][1])
                            print("Found tau:", syst_tau)
                            syst_l_scanner.plot_scan_L_curve(output_filename="%s/scanL_%s.%s" % (syst_output_dir, append, OUTPUT_FMT))
                            syst_l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_%s.%s" % (syst_output_dir, append, OUTPUT_FMT))

                        elif REGULARIZE == "tau":
                            print("Regularizing with ScanTau, please be patient...")
                            syst_tau_scanner = TauScanner()
                            syst_tau = syst_tau_scanner.scan_tau(tunfolder=syst_unfolder,
                                                                 n_scan=args.nScan,
                                                                 tau_min=region['tau_limits'][angle.var][0],
                                                                 tau_max=region['tau_limits'][angle.var][1],
                                                                 scan_mode=scan_mode,
                                                                 distribution=scan_distribution,
                                                                 axis_steering=unfolder.axisSteering)
                            print("Found tau:", syst_tau)
                            syst_tau_scanner.plot_scan_tau(output_filename="%s/scantau_%s.%s" % (syst_output_dir, append, OUTPUT_FMT))

                        title = "L matrix %s %s region %s" % (jet_algo, region['label'], angle_str)
                        syst_unfolder_plotter.draw_L_matrix(title=title, **syst_plot_args)
                        title = "L^{T}L matrix %s %s region %s" % (jet_algo, region['label'], angle_str)
                        syst_unfolder_plotter.draw_L_matrix_squared(title=title, **syst_plot_args)
                        title = "L * (x - bias vector)\n%s\n%s region\n%s" % (jet_algo, region['label'], angle_str)
                        syst_unfolder_plotter.draw_Lx_minus_bias(title=title, **syst_plot_args)

                    region['model_systematics'][ind]['tau'] = syst_tau

                    # Do unfolding!
                    # --------------------------------------------------------------
                    syst_unfolder.do_unfolding(syst_tau)
                    syst_unfolder.get_output(hist_name="syst_%s_unfolded_1d" % syst_label_no_spaces)
                    syst_unfolder._post_process()
                    if not args.jacobian:
                        syst_unfolder.setup_normalised_results_per_pt_bin()

                    syst_title = "%s\n%s region, %s, %s input" % (jet_algo, region['label'], angle_str, syst_label)
                    syst_unfolder_plotter.draw_unfolded_1d(title=syst_title, **syst_plot_args)

                    if REGULARIZE != "None":
                        # Draw our new template alongside unfolded
                        ocs = [
                            Contribution(alt_hist_mc_gen, label=region['alt_mc_label'],
                                         line_color=ROOT.kViolet+1,
                                         marker_color=ROOT.kViolet+1,
                                         subplot=syst_unfolder.hist_truth),
                            Contribution(syst_unfolder.truth_template, label="Template",
                                         line_color=ROOT.kAzure+1,
                                         marker_color=ROOT.kAzure+1,
                                         subplot=syst_unfolder.hist_truth),
                        ]
                        title = "%s\n%s region, %s" % (jet_algo, region['label'], angle_str)
                        syst_unfolder_plotter.draw_unfolded_1d(do_gen=True,
                                                               do_unfolded=True,
                                                               other_contributions=ocs,
                                                               output_dir=syst_output_dir,
                                                               append='syst_with_template',
                                                               title='',
                                                               subplot_title="* / Generator")

                        syst_unfolder_plotter.draw_bias_vector(title=syst_title, **syst_plot_args)
                        syst_unfolder_plotter.draw_x_minus_bias(title=syst_title, **syst_plot_args)

                    region['model_systematics'][ind]['unfolder'] = syst_unfolder

                    cu.close_tfile(syst_dict['tfile'])
                    del syst_dict['tfile']

                    # Do 1D plot of nominal vs syst unfolded
                    # --------------------------------------------------------------
                    # entries = []
                    # # add nominal
                    # label = 'MC' if MC_INPUT else "data"
                    # entries.append(
                    #     Contribution(unfolder.unfolded, label="Unfolded (#tau = %.3g)" % (unfolder.tau),
                    #                  line_color=ROOT.kRed, line_width=1,
                    #                  marker_color=ROOT.kRed, marker_size=0.6, marker_style=20,
                    #                  subplot_line_color=ROOT.kRed, subplot_line_width=1,
                    #                  subplot_marker_color=ROOT.kRed, subplot_marker_size=0, subplot_marker_style=20,
                    #                  normalise_hist=False, subplot=unfolder.hist_truth),
                    # )

                    # entries.append(
                    #     Contribution(unfolder.hist_truth, label="Generator",
                    #                  line_color=ROOT.kBlue, line_width=1,
                    #                  marker_color=ROOT.kBlue, marker_size=0,
                    #                  normalise_hist=False),
                    # )
                    # # add systematic
                    # entries.append(
                    #     Contribution(syst_unfolder.unfolded, label="Unfolded %s (#tau = %.3g)" % (syst_label, syst_unfolder.tau),
                    #                  line_color=syst_dict['colour'], line_width=1,
                    #                  marker_color=syst_dict['colour'], marker_size=0.6, marker_style=20+ind+1,
                    #                  subplot_line_color=syst_dict['colour'], subplot_line_width=1,
                    #                  subplot_marker_color=syst_dict['colour'], subplot_marker_size=0, subplot_marker_style=20,
                    #                  normalise_hist=False, subplot=syst_unfolder.hist_truth),
                    # )
                    # syst_entries.append(entries[-1])

                    # entries.append(
                    #     Contribution(syst_unfolder.hist_truth, label="Generator (%s)" % (syst_label),
                    #                  line_color=syst_dict['colour']+2, line_width=1, line_style=1,
                    #                  marker_color=syst_dict['colour']+2, marker_size=0,
                    #                  normalise_hist=False),
                    # )
                    # syst_entries.append(entries[-1])

                    # title = "%s\n%s region, %s, %s input" % (jet_algo, region['label'], angle_str, syst_label)
                    # plot = Plot(entries,
                    #             what='hist',
                    #             title=title,
                    #             xtitle="Generator bin number",
                    #             ytitle="N",
                    #             subplot_type='ratio',
                    #             subplot_title='Unfolded / Gen',
                    #             subplot_limits=(0, 2),
                    #             has_data=not MC_INPUT)
                    # plot.default_canvas_size = (800, 600)
                    # plot.plot("NOSTACK HISTE")
                    # plot.set_logy(do_more_labels=False, override_check=True)
                    # ymax = max([o.GetMaximum() for o in plot.contributions_objs])
                    # plot.container.SetMaximum(ymax * 200)
                    # ymin = max([o.GetMinimum(1E-10) for o in plot.contributions_objs])
                    # plot.container.SetMinimum(ymin*0.01)
                    # l, t = syst_unfolder_plotter.draw_pt_binning_lines(plot, which='gen', axis='x',
                    #                                                    do_underflow=True,
                    #                                                    do_labels_inside=True,
                    #                                                    do_labels_outside=False,
                    #                                                    labels_inside_align='lower'
                    #                                                    )
                    # plot.legend.SetY1NDC(0.77)
                    # plot.legend.SetY2NDC(0.88)
                    # plot.legend.SetX1NDC(0.65)
                    # plot.legend.SetX2NDC(0.88)
                    # output_filename = "%s/unfolded_1d_modelSyst_%s.%s" % (syst_output_dir, syst_label_no_spaces, syst_unfolder_plotter.output_fmt)
                    # plot.save(output_filename)

            # Load model (scale) systs from another reference file, and calc fractional
            # uncertainty, apply to this unfolded result
            if args.doModelSystsFromFile is not None:
                model_dir = "%s/%s/%s" % (args.doModelSystsFromFile, region['name'], angle.var)

                this_pkl_filename = os.path.join(model_dir, "unfolding_result.pkl")
                if not os.path.isfile(this_pkl_filename):
                    raise IOError("Cannot find model systematics file, %s" % this_pkl_filename)
                model_syst_region = unpickle_region(this_pkl_filename)

                # update original region object with the model syst info from the reference file
                region['model_systematics'] = model_syst_region['model_systematics']

                # Add dummy object to hist_bin_chopper for later, so we can directly manipulate the cache
                # this replaces the procedure in create_normalised_scale_syst_uncertainty_per_pt_bin()
                unfolder.hist_bin_chopper.add_obj(unfolder.scale_uncert_name, unfolder.get_unfolded_with_ematrix_stat())
                unfolder.hist_bin_chopper.add_obj("unfolded_stat_err", unfolder.get_unfolded_with_ematrix_stat())

                for ibin_pt in range(len(unfolder.pt_bin_edges_gen[:-1])):
                    # Calculate scale uncertainty by taking relative uncertainty
                    # from reference file (which was already calculated in the 1st running),
                    # and applying to this result
                    key = unfolder.hist_bin_chopper._generate_key(unfolder.scale_uncert_name,
                                                                  ind=ibin_pt,
                                                                  axis='pt',
                                                                  do_norm=True,
                                                                  do_div_bin_width=True,
                                                                  binning_scheme='generator')
                    ref_scale_syst = model_syst_region['unfolder'].hist_bin_chopper._cache[key]
                    scale_syst = unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width("unfolded_stat_err", ibin_pt, binning_scheme='generator').Clone()
                    for i in range(1, scale_syst.GetNbinsX()+1):
                        if ref_scale_syst.GetBinContent(i) != 0:
                            rel_err = ref_scale_syst.GetBinError(i) / ref_scale_syst.GetBinContent(i)
                        else:
                            rel_err = 0
                        scale_syst.SetBinError(i, rel_err * scale_syst.GetBinContent(i))
                    unfolder.hist_bin_chopper._cache[key] = scale_syst

                unfolder.create_normalised_scale_syst_ematrices_per_pt_bin()

            if len(region['model_systematics']) > 0 and MC_INPUT:
                # Do a big absolute 1D plots for sanity
                # Detector level
                model_contributions = [
                    Contribution(mdict['unfolder'].input_hist_bg_subtracted,
                                 label=mdict['label'],
                                 line_color=mdict['colour'], line_style=1, line_width=1,
                                 marker_color=mdict['colour'], marker_size=0, marker_style=21,
                                 subplot=unfolder.input_hist_bg_subtracted)
                    for mdict in region['model_systematics']
                ]
                unfolder_plotter.draw_detector_1d(do_reco_mc_bg_sub=True,
                                                  output_dir=this_output_dir,
                                                  append='bg_fakes_subtracted_model_systs_%s' % append,
                                                  title=title,
                                                  other_contributions=model_contributions)
                # Truth level vs nominal
                model_contributions = [
                    Contribution(mdict['unfolder'].get_unfolded_with_ematrix_stat(),
                                 label=mdict['label'],
                                 line_color=mdict['colour'], line_style=1, line_width=1,
                                 marker_color=mdict['colour'], marker_size=0, marker_style=21,
                                 subplot=unfolder.get_unfolded_with_ematrix_stat())
                    for mdict in region['model_systematics']
                ]
                unfolder_plotter.draw_unfolded_1d(do_unfolded=True, do_gen=False,
                                                  output_dir=this_output_dir,
                                                  append='model_systs_%s' % append,
                                                  title=title,
                                                  other_contributions=model_contributions,
                                                  subplot_title='#splitline{Variation /}{nominal}')
                # Truth vs own gen
                model_contributions = [
                    Contribution(mdict['unfolder'].get_unfolded_with_ematrix_stat(),
                                 label=mdict['label'],
                                 line_color=mdict['colour'], line_style=1, line_width=1,
                                 marker_color=mdict['colour'], marker_size=0, marker_style=21,
                                 subplot=mdict['unfolder'].hist_truth)
                    for mdict in region['model_systematics']
                ]
                unfolder_plotter.draw_unfolded_1d(do_unfolded=True, do_gen=False,
                                                  output_dir=this_output_dir,
                                                  append='model_systs_vs_gen_%s' % append,
                                                  title=title,
                                                  other_contributions=model_contributions,
                                                  subplot_title='#splitline{Unfolded /}{Generator}')

            prof_end_model_systs()

            # ------------------------------------------------------------------
            # DO PDF VARIATIONS
            # ------------------------------------------------------------------
            if args.doPDFSysts:
                # first construct all new systematic variations dicts
                original_pdf_dict = region['pdf_systematics'][0]
                original_pdf_dict['label'] = '_PDF_template'  # initial _ to ignore it later on
                pdf_tfile = original_pdf_dict['tfile']
                if not isinstance(pdf_tfile, ROOT.TFile):
                    pdf_tfile = cu.open_root_file(pdf_tfile)

                region['pdf_systematics'] = []
                num_vars = len(original_pdf_dict['variations'])
                for pdf_ind in original_pdf_dict['variations']:
                    region['pdf_systematics'].append(
                        {
                            "label": "PDF_%d" % pdf_ind,
                            "response_map": "%s/tu_%s_GenReco_all_PDF_%d" % (region['dirname'], angle_shortname, pdf_ind),  # only store obj name, not obj, as we don't need all simultaneously
                            "colour": cu.get_colour_seq(pdf_ind, num_vars)
                        })

                # Now run over all variations, unfolding the nominal inputs
                # but with the various response matrices
                for ind, pdf_dict in enumerate(region['pdf_systematics']):
                    pdf_label = pdf_dict['label']
                    pdf_label_no_spaces = cu.no_space_str(pdf_label)

                    print("*" * 80)
                    print("*** Unfolding with PDF variation:", pdf_label, "(%d/%d) ***" % (ind+1, len(region['pdf_systematics'])))
                    print("*" * 80)
                    pdf_response_map = cu.get_from_tfile(pdf_tfile, pdf_dict['response_map'])
                    for merge_func in th2_merge_funcs:
                        pdf_response_map = merge_func(pdf_response_map)
                    pdf_unfolder = MyUnfolder(response_map=pdf_response_map,
                                              binning_handler=unfolder.binning_handler,
                                              **unfolder_args)

                    # Needed because some fo the PDF variations struggle to unfold
                    # Even 1E-18 wouldn't work - needs to be v.small
                    pdf_unfolder.SetEpsMatrix(eps_matrix)

                    pdf_unfolder_plotter = MyUnfolderPlotter(pdf_unfolder, is_data=not MC_INPUT)
                    pdf_output_dir = this_output_dir+"/pdfSyst/"+pdf_label_no_spaces
                    pdf_plot_args = dict(output_dir=pdf_output_dir, append=append)

                    # Set what is to be unfolded
                    # --------------------------------------------------------------
                    # Same input as nominal unfolder, since we only change response matrix
                    pdf_unfolder.set_input(**input_args())

                    # Subtract fakes (treat as background)
                    # --------------------------------------------------------------
                    pdf_unfolder.subtract_background(hist_fakes_reco, "Signal fakes", scale=1., scale_err=0.0)
                    pdf_unfolder.subtract_background_gen_binning(hist_fakes_reco_gen_binning, "Signal fakes", scale=1.)

                    # Do any regularization
                    # --------------------------------------------------------------
                    # Setup L matrix
                    pdf_tau = 0
                    if REGULARIZE != "None":
                        pdf_unfolder.SetBias(unfolder.truth_template)
                        for L_args in unfolder.L_matrix_entries:
                            pdf_unfolder.AddRegularisationCondition(*L_args)

                        if REGULARIZE == "L":
                            print("Regularizing with ScanLcurve, please be patient...")
                            pdf_l_scanner = LCurveScanner()
                            pdf_tau = pdf_l_scanner.scan_L(tunfolder=pdf_unfolder,
                                                           n_scan=args.nScan,
                                                           tau_min=region['tau_limits'][angle.var][0],
                                                           tau_max=region['tau_limits'][angle.var][1])
                            print("Found tau:", pdf_tau)
                            pdf_l_scanner.plot_scan_L_curve(output_filename="%s/scanL_%s.%s" % (pdf_output_dir, append, OUTPUT_FMT))
                            pdf_l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_%s.%s" % (pdf_output_dir, append, OUTPUT_FMT))

                        elif REGULARIZE == "tau":
                            print("Regularizing with ScanTau, please be patient...")
                            pdf_tau_scanner = TauScanner()
                            pdf_tau = pdf_tau_scanner.scan_tau(tunfolder=pdf_unfolder,
                                                               n_scan=args.nScan,
                                                               tau_min=region['tau_limits'][angle.var][0],
                                                               tau_max=region['tau_limits'][angle.var][1],
                                                               scan_mode=scan_mode,
                                                               distribution=scan_distribution,
                                                               axis_steering=unfolder.axisSteering)
                            print("Found tau:", pdf_tau)
                            pdf_tau_scanner.plot_scan_tau(output_filename="%s/scantau_%s.%s" % (pdf_output_dir, append, OUTPUT_FMT))

                    region['pdf_systematics'][ind]['tau'] = pdf_tau

                    # Do unfolding!
                    # --------------------------------------------------------------
                    pdf_unfolder.do_unfolding(pdf_tau)
                    pdf_unfolder.get_output(hist_name="syst_%s_unfolded_1d" % pdf_label_no_spaces)
                    pdf_unfolder._post_process()

                    if not args.jacobian:
                        pdf_unfolder.setup_normalised_results_per_pt_bin()

                    pdf_title = "%s\n%s region, %s\n%s response matrix" % (jet_algo, region['label'], angle_str, pdf_label)
                    pdf_unfolder_plotter.draw_unfolded_1d(title=pdf_title, **pdf_plot_args)

                    prof_end_one_pdf()

                    pdf_unfolder.slim_down()
                    region['pdf_systematics'][ind]['unfolder'] = pdf_unfolder

                    prof_end_one_pdf_tidy()

                if not args.jacobian:
                    unfolder.create_normalised_pdf_syst_uncertainty_per_pt_bin(region['pdf_systematics'])
                    unfolder.create_normalised_pdf_syst_ematrices_per_pt_bin()

                    unfolder.create_pdf_syst_uncertainty_per_pt_bin(region['pdf_systematics'])

                    unfolder.create_pdf_syst_uncertainty_per_lambda_bin(region['pdf_systematics'])

                cu.close_tfile(pdf_tfile)

            prof_end_pdf_systs()

            # Load PDF syst from another reference file, and calc fractional
            # uncertainty, apply to this unfolded result
            if args.doPDFSystsFromFile is not None:
                pdf_dir = "%s/%s/%s" % (args.doPDFSystsFromFile, region['name'], angle.var)

                this_pkl_filename = os.path.join(pdf_dir, "unfolding_result.pkl")
                if not os.path.isfile(this_pkl_filename):
                    raise IOError("Cannot find PDF systematics file, %s" % this_pkl_filename)
                pdf_syst_region = unpickle_region(this_pkl_filename)

                # update original region object from reference file
                region['pdf_systematics'] = pdf_syst_region['pdf_systematics']

                # Add dummy object to hist_bin_chopper for later, so we can directly manipulate the cache
                unfolder.hist_bin_chopper.add_obj(unfolder.pdf_uncert_name, unfolder.get_unfolded_with_ematrix_stat())
                unfolder.hist_bin_chopper.add_obj("unfolded_stat_err", unfolder.get_unfolded_with_ematrix_stat())

                for ibin_pt in range(len(unfolder.pt_bin_edges_gen[:-1])):
                    # Calculate PDF uncertainty by taking relative uncertainty
                    # from reference file (which has already been calculated in the 1st running),
                    # and applying to this result
                    key = unfolder.hist_bin_chopper._generate_key(unfolder.pdf_uncert_name,
                                                                  ind=ibin_pt,
                                                                  axis='pt',
                                                                  do_norm=True,
                                                                  do_div_bin_width=True,
                                                                  binning_scheme='generator')
                    ref_pdf_syst = pdf_syst_region['unfolder'].hist_bin_chopper._cache[key]
                    pdf_syst = unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width("unfolded_stat_err", ibin_pt, binning_scheme='generator').Clone()
                    for i in range(1, pdf_syst.GetNbinsX()+1):
                        if ref_pdf_syst.GetBinContent(i) != 0:
                            rel_err = ref_pdf_syst.GetBinError(i) / ref_pdf_syst.GetBinContent(i)
                        else:
                            rel_err = 0
                        pdf_syst.SetBinError(i, rel_err * pdf_syst.GetBinContent(i))
                    unfolder.hist_bin_chopper._cache[key] = pdf_syst

                if not args.jacobian:
                    unfolder.create_normalised_pdf_syst_ematrices_per_pt_bin()

                # cleanup, save memory
                del pdf_syst_region

            if len(region['pdf_systematics']) > 0 and MC_INPUT:
                # Do a big absolute 1D plot for sanity
                # Don't care about inputs, all the same as nominal
                # pdf_contributions = [
                #     Contribution(pdict['unfolder'].input_hist_bg_subtracted,
                #                  label=pdict['label'],
                #                  line_color=pdict['colour'], line_style=1, line_width=1,
                #                  marker_color=pdict['colour'], marker_size=0, marker_style=21,
                #                  subplot=unfolder.input_hist_bg_subtracted)
                #     for pdict in region['pdf_systematics']
                # ]
                # unfolder_plotter.draw_detector_1d(do_reco_mc_bg_sub=True,
                #                                   output_dir=this_output_dir,
                #                                   append='bg_fakes_subtracted_pdf_systs_%s' % append,
                #                                   title=title,
                #                                   other_contributions=pdf_contributions)

                pdf_contributions = [
                    Contribution(pdict['unfolder'].get_unfolded_with_ematrix_stat(),
                                 label=pdict['label'],
                                 line_color=pdict['colour'], line_style=1, line_width=1,
                                 marker_color=pdict['colour'], marker_size=0, marker_style=21,
                                 subplot=unfolder.get_unfolded_with_ematrix_stat())
                    for pdict in region['pdf_systematics']
                ]
                unfolder_plotter.draw_unfolded_1d(do_unfolded=True, do_gen=False,
                                                  output_dir=this_output_dir,
                                                  append='pdf_systs_%s' % append,
                                                  title=title,
                                                  other_contributions=pdf_contributions,
                                                  subplot_title='#splitline{Variation /}{nominal}')
                pdf_contributions = [
                    Contribution(pdict['unfolder'].get_unfolded_with_ematrix_stat(),
                                 label=pdict['label'],
                                 line_color=pdict['colour'], line_style=1, line_width=1,
                                 marker_color=pdict['colour'], marker_size=0, marker_style=21,
                                 subplot=pdict['unfolder'].hist_truth)
                    for pdict in region['pdf_systematics']
                ]
                unfolder_plotter.draw_unfolded_1d(do_unfolded=True, do_gen=False,
                                                  output_dir=this_output_dir,
                                                  append='pdf_systs_vs_gen_%s' % append,
                                                  title=title,
                                                  other_contributions=pdf_contributions,
                                                  subplot_title='#splitline{Unfolded /}{generator}')

            # ------------------------------------------------------------------
            # Finally update absolute/normalised results
            # ------------------------------------------------------------------
            if not args.jacobian:
                unfolder.setup_normalised_results_per_pt_bin()

                unfolder.setup_absolute_results_per_pt_bin()
                if not any([args.mergeBins, args.mergeBinsFromFile, args.mergeBadResponseBins]):
                    unfolder.setup_absolute_results_per_lambda_bin()

            else:
                unfolder.create_normalisation_jacobian_np()
                unfolder_plotter.draw_jacobian(title="Jacobian", **plot_args)

                if len(region['scale_systematics']) > 0:
                    unfolder.create_absolute_scale_uncertainty(region['scale_systematics'])
                    unfolder_plotter.draw_error_matrix_scale(title="Scale error matrix", **plot_args)

                if len(region['pdf_systematics']) > 0:
                    unfolder.create_absolute_pdf_uncertianty(region['pdf_systematics'])
                    unfolder_plotter.draw_error_matrix_pdf(title="PDF error matrix", **plot_args)

                unfolder.create_absolute_total_uncertainty()
                unfolder_plotter.draw_error_matrix_total_abs(title="Total absolute error matrix", **plot_args)

                unfolder.setup_absolute_results()
                unfolder.normalise_all_systs()
                unfolder_plotter.draw_error_matrix_total_norm(title="Total normalised error matrix", **plot_args)

                unfolder.setup_normalised_results()

            region['unfolder'] = unfolder
            prof_end_all_processing()

            # ------------------------------------------------------------------
            # Save everything to pickle / TFile
            # ------------------------------------------------------------------
            print("")
            print("region sizes:")
            print("-"*80)
            cu.print_dict_item_sizes(region, recursive=True)
            print("-"*80)

            print("")
            print("unfolder attr sizes:")
            print("-"*80)
            cu.print_dict_item_sizes(unfolder.__dict__, recursive=True)
            print("-"*80)

            prof_begin_save_to_pickle()
            pickle_filename = os.path.join("%s/unfolding_result.pkl" % this_output_dir)
            pickle_region(region, pickle_filename, infos=False, convert_tfile_to_str=True)
            print(">> Saved to pickle file", pickle_filename)
            print("")

            prof_begin_save_to_root()
            print(">> Saving unfolder to ROOT file")
            # print(unfolder.hist_bin_chopper.objects)
            # print(unfolder.hist_bin_chopper._cache)
            unfolder.save_unfolded_binned_hists_to_tfile(this_slim_tdir)

            # test the pickle file by un-pickling it
            # print("Testing pickled file...")
            # data = unpickle_region(pickle_filename)
            # print("")
            # print("...unpickled data:")
            # print("    ", data)
            # print("")
            # # print("...data['unfolder'].hist_bin_chopper.objects:")
            # # print("    ", data['unfolder'].hist_bin_chopper.objects)
            # # print("")
            # print("unpickled region attr sizes:")
            # print("-"*80)
            # cu.print_dict_item_sizes(data)
            # print("-"*80)
            # print("unpickled unfolder attr sizes:")
            # print("-"*80)
            # cu.print_dict_item_sizes(data['unfolder'].__dict__)
            # print("-"*80)

            # ------------------------------------------------------------------
            # PLOT LOTS OF THINGS
            # ------------------------------------------------------------------
            if not args.noBinnedPlots:
                setup = Setup(jet_algo=jet_algo,
                              region=region,
                              angle=angle,
                              output_dir=this_output_dir,
                              has_data=not MC_INPUT)
                hbc = do_binned_plots_per_region_angle(setup,
                                                       do_binned_gen_pt=True,
                                                       do_binned_gen_lambda=not args.jacobian,
                                                       do_binned_reco_pt=not args.jacobian)

                do_all_big_normalised_1d_plots_per_region_angle(setup, hbc)

            print("Saved minimal hists to", output_tfile_slim.GetName())
            output_tfile_slim.Close()


# these are dummy functions for use with mprof, to mark various points in the program flow
# takes advantage of the fact that memory_profiler marks times when functions called
@profile
def prof_end_nominal():
    pass

@profile
def prof_end_exp_systs():
    pass

@profile
def prof_end_scale_systs():
    pass

@profile
def prof_end_model_systs():
    pass

@profile
def prof_end_pdf_systs():
    pass

@profile
def prof_end_one_pdf():
    pass

@profile
def prof_end_one_pdf_tidy():
    pass

@profile
def prof_end_all_processing():
    pass

@profile
def prof_begin_save_to_pickle():
    pass

@profile
def prof_begin_save_to_root():
    pass


if __name__ == "__main__":
    main()
