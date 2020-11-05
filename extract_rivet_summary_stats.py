#!/usr/bin/env python


"""
Calculate summary stats from RIVET YODA files: per pT bin, per lambda bin, summary plot

Output is saved to HDF5 file, to be used in do_summary_unfolding_plots.py
"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse
import pandas as pd
import numpy as np
import uproot
import yoda

# my packages
import common_utils as cu
import rivet_naming as rn
import metric_calculators as metrics
from extract_unfolding_summary_stats import unpack_slim_unfolding_root_file_uproot, scale_ematrix_by_bin_widths, check_hist_for_negatives


def normalize_areas(areas, errors):
    """Normalise the areas & errors such that they integrate to unity

    TODO: use YODA's Histo1D::normalize() ?

    Parameters
    ----------
    areas : np.array
        Bin areas
    errors : np.array
        Bin errors
    """
    integral = areas.sum()
    if integral > 0:
        scale_factor = 1./integral
        areas *= scale_factor
        errors *= scale_factor


def get_yoda_stats_dict(input_filename,
                        key_label,
                        data_ak4_dirname=None,
                        data_ak8_dirname=None,
                        ignore_missing=False):
    """Get summary statistics from YODA file full of histograms

    Parameters
    ----------
    input_filename : str
        YODA filename
    key_label : str
        Label to append to column names, e.g. mean_X, rms_err_X
    data_ak4_dirname : None, optional
        Name of directory with AK4 data results for delta calculation
    data_ak8_dirname : None, optional
        Name of directory with AK8 data results for delta calculation
    ignore_missing : bool, optional
        If True, fill in 0 values for missing hists; otherwise raise KeyError

    Returns
    -------
    list(dict)
        List of dicts. Each dict represents one radius/region/angle combination,
        with keys/values as metric names/values

    Raises
    ------
    KeyError
        Description
    """
    yoda_dict = yoda.read(input_filename)
    hname = list(yoda_dict.keys())[0]
    is_dijet = rn.DIJET_PATH in hname

    regions = rn.DIJET_REGIONS if is_dijet else rn.ZPJ_REGIONS
    pt_bins = rn.PT_BINS_DIJET if is_dijet else rn.PT_BINS_ZPJ
    path = rn.DIJET_PATH if is_dijet else rn.ZPJ_PATH
    results_dicts = []
    # print("parsing", input_filename)
    # print("regions:", regions)
    for jet_radius in rn.JET_RADII:
        algo_name = "%spuppi" % jet_radius.name.lower()
        for region in regions:
            for lambda_var in rn.LAMBDA_VARS:

                # Get data results for delta calculation
                unfolding_dict = None
                data_dirname = None
                if jet_radius.name.lower() == "ak4" and data_ak4_dirname:
                    data_dirname = data_ak4_dirname
                elif jet_radius.name.lower() == "ak8" and data_ak8_dirname:
                    data_dirname = data_ak8_dirname
                if data_dirname:
                    angle_output_dir = "%s/%s" % (region.name, lambda_var.hist_name)
                    root_filename = os.path.join(data_dirname, angle_output_dir, "unfolding_result_slim.root")
                    uproot_file = uproot.open(root_filename)
                    unfolding_dict = unpack_slim_unfolding_root_file_uproot(uproot_file, region.name, lambda_var.hist_name, pt_bins)

                # Iterate over each pt bin for this radius/region/var, get metrics
                for ibin, pt_bin in enumerate(pt_bins):
                    # print(jet_radius, lambda_var, region, ibin, pt_bin)

                    yoda_name = rn.get_plot_name(jet_radius, region, lambda_var, pt_bin)
                    yoda_name = "/%s/%s" % (path, yoda_name)
                    if yoda_name not in yoda_dict:
                        if ignore_missing:
                            print("Warning, missing", yoda_name, "filling with 0s")
                            result_dict = {
                                'jet_algo': algo_name,
                                'region': region.name,
                                'isgroomed': region.is_groomed,
                                'pt_bin': ibin,
                                'angle': lambda_var.hist_name,

                                'mean_%s' % key_label: 0,
                                'mean_err_%s' % key_label: 0,

                                'rms_%s' % key_label: 0,
                                'rms_err_%s' % key_label: 0,

                                'delta_%s' % key_label: 0,
                                'delta_err_%s' % key_label: 0,
                            }
                            results_dicts.append(result_dict)
                        else:
                            raise KeyError("Missing hist %s" % yoda_name)

                    hist = yoda_dict[yoda_name]
                    areas, widths, centers, errors = metrics.yoda_hist_to_values(hist)
                    normalize_areas(areas, errors)
                    areas, centers = metrics.hist_values_to_uarray(areas, centers, errors)
                    mean_u = metrics.calc_mean_ucert(areas, centers)
                    mean, mean_err = mean_u.nominal_value, mean_u.std_dev
                    rms_u = metrics.calc_rms_ucert(areas, centers)
                    rms, rms_err = rms_u.nominal_value, rms_u.std_dev

                    if "d11-x02-y13" in yoda_name:
                        print("d11-x02-y13:", mean, mean_err)

                    if not (0 < mean < 50):
                        print("WARNING bad mean:", mean, yoda_name)

                    delta, delta_err = 0, 0
                    if unfolding_dict:
                        unfolded_hist_bin_total_errors = unfolding_dict['unfolding_total_err_hists'][ibin]
                        ematrix = scale_ematrix_by_bin_widths(unfolding_dict['unfolding_total_ematrices'][ibin].values,
                                                              metrics.get_uproot_th1_bin_widths(unfolded_hist_bin_total_errors))
                        check_hist_for_negatives(unfolded_hist_bin_total_errors)
                        areas_a, widths_a, centers_a, errors_a = metrics.uproot_th1_to_arrays(unfolded_hist_bin_total_errors)
                        areas_b, widths_b, centers_b, errors_b = metrics.yoda_hist_to_values(hist)
                        print(centers_a)
                        print(widths_a)
                        print(centers_b)
                        print(widths_b)
                        delta = metrics.calc_delta_jax(areas_a, areas_b)
                        err = metrics.calc_delta_correlated_error_jax(areas_a, ematrix, areas_b, errors_b)

                    result_dict = {
                        'jet_algo': algo_name,
                        'region': region.name,
                        'isgroomed': region.is_groomed,
                        'pt_bin': ibin,
                        'angle': lambda_var.hist_name,

                        'mean_%s' % key_label: mean,
                        'mean_err_%s' % key_label: mean_err,

                        'rms_%s' % key_label: rms,
                        'rms_err_%s' % key_label: rms_err,

                        'delta_%s' % key_label: delta,
                        'delta_err_%s' % key_label: delta_err,
                    }
                    results_dicts.append(result_dict)
    return results_dicts


def dataframe_yoda_key(title):
    """Sanitise key for e.g. pandas column names (no spaces, periods, etc)"""
    return title.replace(" ", "_").replace("=", "_").replace(".", "p")


def convert_df_types(df):
    """Convert dataframe column types"""
    df['jet_algo'] = df['jet_algo'].astype('category')
    df['isgroomed'] = df['isgroomed'].astype('bool')
    df['region'] = df['region'].astype('category')
    df['angle'] = df['angle'].astype('category')


def create_yoda_dataframe(yoda_stats_dicts):
    """Create dataframe entries from list of dictionaries of stats"""
    df = pd.DataFrame(yoda_stats_dicts)
    convert_df_types(df)
    return df


def get_dataframe_from_yoda_inputs(yoda_input_groups,
                                   data_ak4_dirname=None,
                                   data_ak8_dirname=None,
                                   ignore_missing=False):
    """Create dataframe from yoda input files & labels

    Each entry in yoda_input_groups is a [dijet filename, Z+J filename, label]
    """
    df = None
    for yin_dijet, yin_zpj, ylabel in yoda_input_groups:
        yoda_stats_dicts = []
        # handle dijet first
        print("Processing", ylabel, yin_dijet)
        stats_dict_dijet = get_yoda_stats_dict(yin_dijet,
                                               key_label=dataframe_yoda_key(ylabel),
                                               data_ak4_dirname=data_ak4_dirname,
                                               data_ak8_dirname=data_ak8_dirname,
                                               ignore_missing=ignore_missing)
        yoda_stats_dicts.extend(stats_dict_dijet)

        # then Z+J
        print("Processing", ylabel, yin_zpj)
        stats_dict_zpj = get_yoda_stats_dict(yin_zpj,
                                             key_label=dataframe_yoda_key(ylabel),
                                             data_ak4_dirname=data_ak4_dirname,
                                             data_ak8_dirname=data_ak8_dirname,
                                             ignore_missing=ignore_missing)
        yoda_stats_dicts.extend(stats_dict_zpj)

        # Convert to dataframe, merge in with main dataframe
        df_yoda = create_yoda_dataframe(yoda_stats_dicts)
        if df is None:
            df = df_yoda
        else:
            df = pd.merge(df, df_yoda, how='outer')
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5output",
                        default='rivet_summary.h5',
                        help=("Output HDF5 filename. Default is 'rivet_summary.h5'"))
    parser.add_argument("--ak4source",
                        default=None,
                        help="Source directory for AK4 jets (should be the one made by unfolding.py")
    parser.add_argument("--ak8source",
                        default=None,
                        help="Source directory for AK8 jets (should be the one made by unfolding.py")
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
    parser.add_argument("--ignoreMissing",
                        action='store_true',
                        help='Ignore missing bins/values')
    args = parser.parse_args()

    if (len(args.yodaInputDijet) != len(args.yodaLabel)
        and len(args.yodaInputZPJ) != len(args.yodaLabel)):
        raise RuntimeError("Number of --yodaInputDijet/yodaInputZPJ must match number of --yodaLabel")

    if len(args.yodaInputDijet) == 0:
        raise RuntimeError("Need at least one YODA input")

    # -----------------------------------------------------------------------
    # Get stats from YODA files, add to dataframe
    # -----------------------------------------------------------------------
    df = get_dataframe_from_yoda_inputs(zip(args.yodaInputDijet, args.yodaInputZPJ, args.yodaLabel),
                                        data_ak4_dirname=args.ak4source,
                                        data_ak8_dirname=args.ak8source,
                                        ignore_missing=args.ignoreMissing)
    convert_df_types(df)
    print(df.head())
    print(df.tail())
    print(len(df.index), 'entries in dataframe')
    print(df.dtypes)

    # Save dataframe to file
    # -----------------------------------------------------------------------
    print("Saving dataframe to", args.h5output)
    cu.check_dir_exists_create(os.path.dirname(args.h5output))
    # need format='table' to store category dtype
    df.to_hdf(args.h5output, key='df', format='table')
