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

import yoda

# my packages
import common_utils as cu
import rivet_naming as rn
import metric_calculators as metrics


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


def get_yoda_stats_dict(input_filename, key_label, reference_hist=None):
    """Get summary statistics from YODA file full of histograms

    Parameters
    ----------
    input_filename : str
        YODA filename
    key_label : str
        Label to append to column names, e.g. mean_X, rms_err_X
    reference_hist : None, optional
        TODO

    Returns
    -------
    list(dict)
        List of dicts. Each dict represents one radius/region/angle combination,
        with keys/values as metric names/values
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
                for ibin, pt_bin in enumerate(pt_bins):
                    yoda_name = rn.get_plot_name(jet_radius, region, lambda_var, pt_bin)
                    yoda_name = "/%s/%s" % (path, yoda_name)
                    hist = yoda_dict[yoda_name]
                    areas, widths, centers, errors = metrics.yoda_hist_to_values(hist)
                    normalize_areas(areas, errors)
                    areas, centers = metrics.hist_values_to_uarray(areas, centers, errors)
                    mean_u = metrics.calc_mean_ucert(areas, centers)
                    mean, mean_err = mean_u.nominal_value, mean_u.std_dev
                    rms_u = metrics.calc_rms_ucert(areas, centers)
                    rms, rms_err = rms_u.nominal_value, rms_u.std_dev

                    # FIXME: need delta calculation, needs reference hist
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

                        'delta_%s' % key_label: 0,
                        'delta_err_%s' % key_label: 0,
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


def get_dataframe_from_yoda_inputs(yoda_input_groups):
    """Create dataframe from yoda input files & labels

    Each entry in yoda_input_groups is a [dijet filename, Z+J filename, label]
    """
    df = None
    for yin_dijet, yin_zpj, ylabel in yoda_input_groups:
        yoda_stats_dicts = []
        # handle dijet first
        print("Processing", ylabel, yin_dijet)
        stats_dict_dijet = get_yoda_stats_dict(yin_dijet, key_label=dataframe_yoda_key(ylabel))
        yoda_stats_dicts.extend(stats_dict_dijet)

        # then Z+J
        print("Processing", ylabel, yin_zpj)
        stats_dict_zpj = get_yoda_stats_dict(yin_zpj, key_label=dataframe_yoda_key(ylabel))
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
                        default='summary.h5',
                        help=("Output HDF5 filename. Default is 'summary.h5'"))
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
    args = parser.parse_args()

    if (len(args.yodaInputDijet) != len(args.yodaLabel)
        and len(args.yodaInputZPJ) != len(args.yodaLabel)):
        raise RuntimeError("Number of --yodaInputDijet/yodaInputZPJ must match number of --yodaLabel")

    if len(args.yodaInputDijet) == 0:
        raise RuntimeError("Need at least one YODA input")

    # -----------------------------------------------------------------------
    # Get stats from YODA files, add to dataframe
    # -----------------------------------------------------------------------
    df = get_dataframe_from_yoda_inputs(zip(args.yodaInputDijet, args.yodaInputZPJ, args.yodaLabel))
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
