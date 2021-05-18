#!/usr/bin/env python

"""
Create table of summary stats
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
from functools import partial

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
from extract_rivet_summary_stats import dataframe_yoda_key, convert_df_types


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5input",
                        help="Read analysis data from H5 input file "
                             "(from running extract_unfolding_summary_stats.py)",
                        required=True)
    args = parser.parse_args()

    # Read in data from h5 files
    # -----------------------------------------------------------------------
    print("Reading in unfolding data from existing HDF5 file...")
    if not os.path.isfile(args.h5input):
        raise IOError("Cannot find --h5input file")

    with pd.HDFStore(args.h5input) as store:
        df = store['df']
    print(df.head())
    print("# entries:", len(df.index))

    # yoda_labels = []
    # if args.h5inputRivet:
    #     print("Reading in RIVET data from existing HDF5 file...")
    #     with pd.HDFStore(args.h5inputRivet) as store:
    #         df_rivet = store['df']
    #     print(df_rivet.head())
    #     print("# Rivet entries:", len(df_rivet.index))
    #     # Figure out YODA entries from column names
    #     mean_columns = [c.replace("mean_err_", '') for c in df_rivet.columns if 'mean_err_' in c]
    #     print(mean_columns)
    #     df = pd.merge(df, df_rivet, how='outer')
    #     # sort manually, but check the ones we want are actually in the dataframe
    #     yoda_labels_ideal = ['Pythia8_CP2', 'Pythia8_CP5', 'Herwig7_CH3', 'Sherpa']
    #     yoda_labels = []
    #     for yl in yoda_labels_ideal:
    #         if yl not in mean_columns:
    #             warnings.warn("Missing yoda input %s from rivet file" % yl)
    #         else:
    #             yoda_labels.append(yl)
    #     print("Setting yoda_labels to", yoda_labels)

    convert_df_types(df)
    print(df.columns)
    print(df.head())
    print(df.dtypes)
    print(len(df.index), 'entries in dataframe')

    # Select necessary columns
    # --------------------------------------------------------------------------
    # Create table for each jet algo, region, grooming
    jet_algos = df['jet_algo'].unique().tolist()
    regions = df['region'].unique().tolist()
    angles = df['angle'].unique().tolist()

    cols = [
        'pt_bin',
        'mean_with_err',
        'rms_with_err',
    ]

    def rename_var(v):
        for angle in qgc.COMMON_VARS:
            if angle.var == v:
                return "$"+angle.mathmode+"$"
        return v
    
    def float_to_str(num, n_sig_fig=2):
        """Convert num to str with some number of significant figures"""
        numstr = float(("{0:.%ie}" % (n_sig_fig-1)).format(num))
        s = '{}'.format(numstr)
        # add on a missing 0 at the end to show n dp
        s_non_zero = s.strip("0.")
        if len(s_non_zero) != n_sig_fig:
            s += "0"
        return s

    def float_to_str_match_precision(num, err):
        """Convert num to str, such that it has the same precision as err,
        when err is converted to have some set number of sig fig"""
        s = float_to_str(err).replace("0.", "")
        dp = len(s)
        fmt = "{:.%df}" % dp
        return fmt.format(num)

    def get_pt_str(row):
        pt_bins = qgc.PT_UNFOLD_DICT['signal_zpj_gen']
        if "Dijet" in row['region']:
            pt_bins = qgc.PT_UNFOLD_DICT['signal_gen']
        pt_bin_lo = pt_bins[row['pt_bin']]
        pt_bin_hi = pt_bins[row['pt_bin']+1]
        return "%g--%g" % (pt_bin_lo, pt_bin_hi)

    charged_angles = [a for a in angles if '_charged' in a]
    chargedneutral_angles = [a for a in angles if '_charged' not in a]

    for algo, region in product(jet_algos, regions):
        for angles in [chargedneutral_angles, charged_angles]:
            table_data = []
            for angle in angles:
                algo_mask = df['jet_algo'] == algo
                region_mask = df['region'] == region
                angle_mask = df['angle'] == angle
                # regions include "_groomed" regions, for some reason

                this_df = df.loc[algo_mask & region_mask & angle_mask, :]

                if len(this_df.index) == 0:
                    raise RuntimeError("No entries for config %s %s" % (algo, region))

                this_df.loc[:, 'mean_with_err'] = this_df.apply(lambda row: r"$ %s \pm %s$" % (float_to_str_match_precision(row['mean'], row['mean_err']), float_to_str(row['mean_err'])), axis=1)
                this_df.loc[:, 'rms_with_err'] = this_df.apply(lambda row: r"$ %s \pm %s$" % (float_to_str_match_precision(row['rms'], row['rms_err']), float_to_str(row['rms_err'])), axis=1)
                this_df.loc[:, 'pt_bin'] = this_df.apply(get_pt_str, axis=1)

                # Rename some columns
                this_df = this_df.loc[:, cols].rename(
                    {
                        'pt_bin': '$p_{T}$ bin [GeV]',
                        'mean_with_err': r'Mean $\pm$ error',
                        'rms_with_err': r'RMS $\pm$ error',
                    },
                    axis='columns')

                algo_str_dict = {'ak4puppi': 'AK4', 'ak8puppi': 'AK8'}
                region_str_dict = {k: k.replace("_", " ").replace(" groomed", "").replace("ZPlusJets", "Z+Jet") for k in regions}
                grooming_str_dict = {True: "with", False: "without"}

                latex_table = this_df.to_latex(index=False, escape=False, column_format='lcc')
                latex_table = latex_table.replace(r"\toprule", "").replace(r"\midrule", r"\hline").replace(r"\bottomrule", "")
                latex_table = latex_table.replace(r"\begin{table}", "").replace(r"\centering", "").replace(r"\caption{}", "").replace(r"\end{table}", "")

                caption = rename_var(angle)
                table_data.append((latex_table, caption))

            output_filename = "summary_stats_{algo}_{region}_{angle}.tex".format(algo=algo, region=region, angle="charged" if 'charged' in angles[0] else "chargedneutral")
            caption = "Summary statistics for data, for {} jets in the {} region, {} grooming".format(algo_str_dict[algo], region_str_dict[region], grooming_str_dict["_groomed" in region])
            with open(output_filename, "w") as f:
                f.write(r"\begin{table}[htp]"+"\n")
                f.write(r"\topcaption{"+caption+"}"+"\n")
                f.write(r"\centering"+"\n")
                for ind, data in enumerate(table_data):
                    # wrap in resizebox to shrink it to fit 2 columns
                    f.write(r"\resizebox{0.4\columnwidth}{!}{\subfloat[" + data[1] + "]{" + "\n")
                    f.write(data[0])
                    f.write("}}\n")
                    if ind != len(table_data)-1:
                        f.write(r"\qquad"+"\n")

                f.write(r"\end{table}"+"\n")

            print("Written", output_filename)
