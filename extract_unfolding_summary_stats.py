#!/usr/bin/env python


"""
Calculate summary stats from unfolding analysis: per pT bin, per lambda bin, summary plot

Output is saved to HDF5 file, to be used in do_summary_unfolding_plots.py
"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse
import pandas as pd
from copy import copy

# make things blazingly fast
import uproot

import yoda

import ROOT

# my packages
import common_utils as cu
import qg_common as qgc
from my_unfolder import unpickle_region, unpack_slim_unfolding_root_file
from unfolding_config import get_dijet_config, get_zpj_config
import metric_calculators as metrics


ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()


def scale_ematrix_by_bin_widths(ematrix, widths):
    this_widths = widths.reshape(len(widths), 1)
    return ematrix * this_widths * this_widths.T


def check_hist_for_negatives(hist):
    """Check if any bins are < 0

    Parameters
    ----------
    hist : uproot.TH1
    """
    areas, widths, centers, errors = metrics.uproot_th1_to_arrays(hist)
    for i, x in enumerate(areas, 1):
        if x < 0:
            print("WARNING:", hist.name, " has area of bin %d = %f" % (i, x))
            # raise ValueError("Area of bin %d = %f" % (i, x))


def unpack_slim_unfolding_root_file_uproot(input_tfile, region_name, angle_name, pt_bins):
    tdir = "%s/%s" % (region_name, angle_name)
    indices = range(len(pt_bins)-1 if isinstance(pt_bins[0], float) else len(pt_bins))  # adjust if pairs of bins, or individual bin edges

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
    parser.add_argument("--h5output",
                        default='unfolding_summary.h5',
                        help=("Output HDF5 filename. Default is 'unfolding_summary.h5'"))
    args = parser.parse_args()

    # Get input data
    if not any([args.ak4source, args.ak8source]):
        raise RuntimeError("Need --ak4input/--ak8input")

    # ----------------------------------------------------------------------
    # READ IN DATA FROM UNFOLDING ROOT FILES
    # ----------------------------------------------------------------------
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
                    mc_gen_hist_bin_mean, mc_gen_hist_bin_mean_err = metrics.calc_hist_mean_and_uncorrelated_error_jax(mc_gen_hist_bin)
                    mc_gen_hist_bin_rms, mc_gen_hist_bin_rms_err = metrics.calc_hist_rms_and_uncorrelated_error_jax(mc_gen_hist_bin)

                    # Handle alt MC hist -> metrics
                    alt_mc_gen_hist_bin = unfolding_dict['alt_truth_hists'][ibin]
                    try:
                        check_hist_for_negatives(alt_mc_gen_hist_bin)
                    except ValueError as e:
                        print("-ve value for alt MC hist in pt bin", ibin, ":", bin_edge_low, "-", bin_edge_high)
                        raise e
                    alt_mc_gen_hist_bin_mean, alt_mc_gen_hist_bin_mean_err = metrics.calc_hist_mean_and_uncorrelated_error_jax(alt_mc_gen_hist_bin)
                    alt_mc_gen_hist_bin_rms, alt_mc_gen_hist_bin_rms_err = metrics.calc_hist_rms_and_uncorrelated_error_jax(alt_mc_gen_hist_bin)

                    # Handle unfolded data hist -> metrics
                    unfolded_hist_bin_total_errors = unfolding_dict['unfolding_total_err_hists'][ibin]
                    try:
                        check_hist_for_negatives(unfolded_hist_bin_total_errors)
                    except ValueError as e:
                        print("-ve value for data hist in pt bin", ibin, ":", bin_edge_low, "-", bin_edge_high)
                        raise e
                    ematrix = scale_ematrix_by_bin_widths(unfolding_dict['unfolding_total_ematrices'][ibin].values, metrics.get_uproot_th1_bin_widths(unfolded_hist_bin_total_errors))
                    unfolded_hist_bin_total_errors_mean, unfolded_hist_bin_total_errors_mean_err = metrics.calc_hist_mean_and_correlated_error_jax(unfolded_hist_bin_total_errors, ematrix)
                    unfolded_hist_bin_total_errors_rms, unfolded_hist_bin_total_errors_rms_err = metrics.calc_hist_rms_and_correlated_error_jax(unfolded_hist_bin_total_errors, ematrix)

                    delta_nominal, delta_nominal_err = metrics.calc_hist_delta_and_error_jax(unfolded_hist_bin_total_errors, ematrix, mc_gen_hist_bin)
                    delta_alt, delta_alt_err = metrics.calc_hist_delta_and_error_jax(unfolded_hist_bin_total_errors, ematrix, alt_mc_gen_hist_bin)

                    result_dict = {
                        'jet_algo': jet_algo['name'],
                        'region': this_region['name'], # TODO remove "_groomed"?
                        'isgroomed': 'groomed' in this_region['name'].lower(),
                        'pt_bin': ibin,
                        'angle': angle.var,

                        'mean': unfolded_hist_bin_total_errors_mean,
                        'mean_err': unfolded_hist_bin_total_errors_mean_err,

                        'mean_truth': mc_gen_hist_bin_mean,
                        'mean_err_truth': mc_gen_hist_bin_mean_err,

                        'mean_alt_truth': alt_mc_gen_hist_bin_mean,
                        'mean_err_alt_truth': alt_mc_gen_hist_bin_mean_err,

                        'rms': unfolded_hist_bin_total_errors_rms,
                        'rms_err': unfolded_hist_bin_total_errors_rms_err,

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

    if len(results_dicts) == 0:
        raise ValueError("No entries to go into dataframe!")

    df = pd.DataFrame(results_dicts)
    df['jet_algo'] = df['jet_algo'].astype('category')
    df['region'] = df['region'].astype('category')
    df['angle'] = df['angle'].astype('category')
    print(df.head())
    print(df.tail())
    print(len(df.index), 'entries in dataframe')
    print(df.dtypes)

    print("Saving dataframe to", args.h5output)
    cu.check_dir_exists_create(os.path.dirname(args.h5output))
    # need format='table' to store category dtype
    df.to_hdf(args.h5output, key='df', format='table')
