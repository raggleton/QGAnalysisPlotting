#!/usr/bin/env python


"""Do & print out bottom-line tests on unfolding output"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse
import pandas as pd
from copy import copy
import numpy as np
import scipy
from scipy import optimize
from functools import partial

# for webpage
from jinja2 import Environment, FileSystemLoader

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()

# my packages
import common_utils as cu
import qg_common as qgc
from my_unfolder import unpickle_region
from unfolding_config import setup_regions_from_argparse

pd.set_option('display.max_columns', 0)

# Use rootpy to throw exceptions on ROOT errors, but need DANGER enabled
# Doesn't work with python 3.8 for now
if not (sys.version_info.major == 3 and sys.version_info.minor >= 8):
    import rootpy
    # import rootpy.logger.magic as M; M.DANGER.enabled = True

ROOT.gErrorIgnoreLevel = ROOT.kError
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()


class Setup(object):
    """Loads of common consts, useful for plotting etc"""

    def __init__(self, jet_algo, region, angle, output_dir='.', has_data=False, is_ave_pt_binning=False):
        self.jet_algo = jet_algo
        self.region = region
        do_zpj = "ZPlusJets" in region['name']
        do_dijet = "Dijet" in region['name']
        self.lumi = cu.get_lumi_str(do_dijet=do_dijet, do_zpj=do_zpj)
        self.pt_var_str = "#LT p_{T}^{jet} #GT" if is_ave_pt_binning else "p_{T}^{jet}"
        self.pt_str = self.pt_var_str + " [GeV]"
        self.has_data = has_data
        self.angle = angle
        angle_prepend = "Groomed " if "groomed" in region['name'] else ""
        lower_angle_name = qgc.lower_angle_name(angle)
        # for plot axis titles
        self.angle_str = "{prepend}{name} ({lambda_str})".format(prepend=angle_prepend,
                                                                 name=lower_angle_name if angle_prepend != "" else angle.name,
                                                                 lambda_str=angle.lambda_str)
        # self.particle_title = "Particle-level " + self.angle_str
        self.particle_title = self.angle_str

        angle_prepend = "groomed " if "groomed" in region['name'] else ""
        self.detector_title = "Detector-level {prepend}{name} ({lambda_str})".format(prepend=angle_prepend,
                                                                                     name=lower_angle_name,
                                                                                     lambda_str=angle.lambda_str)
        self.pt_bin_normalised_differential_label = "#frac{1}{d#sigma/dp_{T}} #frac{d^{2}#sigma}{dp_{T} d%s}" % (angle.lambda_str)
        self.pt_bin_detector_normalised_differential_label = "#frac{1}{dN/dp_{T}} #frac{d^{2}N}{dp_{T} d%s}" % (angle.lambda_str)
        self.pt_bin_normalised_differential_times_width_label = "#frac{1}{d#sigma/dp_{T}} #frac{d^{2}#sigma}{dp_{T} d%s}" % (angle.lambda_str)
        self.pt_bin_unnormalised_differential_label = "#frac{1}{dN/dp_{T}} #frac{d^{2}N}{dp_{T} d%s}" % (angle.lambda_str)  # FIXME
        self.pt_bin_unnormalised_differential_label = "N"
        self.lambda_bin_normalised_differential_label = "#frac{1}{d#sigma/d%s} #frac{d^{2}#sigma}{dp_{T} d%s}" % (angle.lambda_str, angle.lambda_str)
        self.lambda_bin_unnormalised_differential_label = "#frac{1}{dN/d%s} #frac{d^{2}N}{dp_{T} d%s}" % (angle.lambda_str, angle.lambda_str)  # FIXME
        self.lambda_bin_unnormalised_differential_label = "N"
        self.output_dir = output_dir
        self.output_fmt = 'pdf'
        self.append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name


def calc_chi2_stats(one_hist, other_hist, cov_matrix):
    one_vec, one_err = cu.th1_to_ndarray(one_hist, False)
    # print(one_err)
    other_vec, _ = cu.th1_to_ndarray(other_hist, False)
    delta = one_vec - other_vec
    if isinstance(cov_matrix, ROOT.TH2):
        v, _ = cu.th2_to_ndarray(cov_matrix)
    else:
        v = cov_matrix
    # print("delta:", delta)
    # v = np.diag(np.diag(v))  # turn off correlations
    # print("v:", v)
    try:
        v_inv = np.linalg.inv(v)
    except np.linalg.LinAlgError:
        print("Trying pseudo-inverse instead")
        v_inv = np.linalg.pinv(v, rcond=1E-30)
    inter = v_inv.dot(delta.T)
    # print("parts:", delta * inter.T)
    chi2 = delta.dot(inter)[0][0]
    ndof = delta.shape[1]
    p = 1-scipy.stats.chi2.cdf(chi2, int(ndof))
    return chi2, ndof, p


def get_null_bins(h):
    null_bins = []
    if isinstance(h, (ROOT.TH1, ROOT.TH2)):
        if isinstance(h, ROOT.TH2):
            h_proj = h.ProjectionX()
        else:
            h_proj = h
        for ix in range(1, h_proj.GetNbinsX()+1):
            if h_proj.GetBinContent(ix) == 0:
                null_bins.append(ix)
        return null_bins
    else:
        proj = h.sum(axis=0)
        return np.where(proj == 0)[0]


def remove_null_bins(arr, null_bins):
    if len(arr.shape) > 1:
        for ax in range(len(arr.shape)):
            if arr.shape[ax] > 1:
                arr = np.delete(arr, null_bins, axis=ax)
    else:
        arr = np.delete(arr, null_bins)
    return arr


def get_bottom_line_stats(setup, no_null_bins=True):
    """Construct dict of bottom-line (i.e. chi2) stats for this region/angle combo"""
    unfolder = setup.region['unfolder']

    # ABSOLUTE VERSION
    # Not sure if this makes sense? Since subject to overall normalisation issues
    # But we can't do it on the Jacobian-transformed normalised version,
    # since the Jacobian is singular (i.e. non-invertible)

    # Do smeared chi2
    # Do folding for nominal MC
    # --------------------------------------------------------------------------
    folded_mc_truth_gen_binning = unfolder.get_ndarray_signal_region_no_overflow(
                                      unfolder.convert_reco_binned_hist_to_gen_binned(
                                          unfolder.get_folded_mc_truth()
                                      ),
                                      xbinning='generator',
                                      ybinning='generator'
                                  )

    input_hist_gen_binning_bg_subtracted_signal = unfolder.get_ndarray_signal_region_no_overflow(
                                                      unfolder.input_hist_gen_binning_bg_subtracted,
                                                      xbinning='generator',
                                                      ybinning='generator'
                                                  )



    vyy = unfolder.make_diag_cov_hist_from_errors(
              h1d=unfolder.input_hist_gen_binning_bg_subtracted,
              inverse=False
          )
    vyy_gen_binning = unfolder.get_ndarray_signal_region_no_overflow(
                          vyy,
                          xbinning='generator',
                          ybinning='generator'
                      )


    vyy_inv = unfolder.make_diag_cov_hist_from_errors(
                  h1d=unfolder.input_hist_gen_binning_bg_subtracted,
                  inverse=True
              )

    vyy_inv_gen_binning = unfolder.get_ndarray_signal_region_no_overflow(
                              vyy_inv,
                              xbinning='generator',
                              ybinning='generator'
                          )

    if no_null_bins:
        null_bins = get_null_bins(vyy_gen_binning)
        print('null bins smeared space:', null_bins)
        vyy_gen_binning = remove_null_bins(vyy_gen_binning, null_bins)
        vyy_inv_gen_binning = remove_null_bins(vyy_inv_gen_binning, null_bins)
        folded_mc_truth_gen_binning = remove_null_bins(folded_mc_truth_gen_binning, null_bins)
        input_hist_gen_binning_bg_subtracted_signal = remove_null_bins(input_hist_gen_binning_bg_subtracted_signal, null_bins)

    # folded_mc_truth_gen_binning = unfolder.convert_reco_binned_hist_to_gen_binned(unfolder.get_folded_mc_truth())
    # input_hist_gen_binning_bg_subtracted_signal = unfolder.input_hist_gen_binning_bg_subtracted
    # vyy_gen_binning = vyy
    # vyy_inv_gen_binning = vyy_inv

    # do fit to scale MC to minimise chi2 (since arbitrary normalisation, want the one that minimises it)
    chi2_kwargs = dict(
        detector_space=False,
        ignore_underflow_bins=False,
        has_underflow=False,
        has_overflow=False,
    )

    def scaled_chi2(f, ref_hist, scaled_hist, cov_inv_matrix):
        """Generic function to calculate chi2 between ref_hist and f*scaled_hist,
        with conv_in_matrix as inverse covariance matrix"""
        # f is scaling factor applied to scaled_hist only
        this_scaled_hist = scaled_hist * f
        chi2, _, _ = unfolder.calculate_chi2(one_hist=this_scaled_hist,
                                             other_hist=ref_hist,
                                             cov_inv_matrix=cov_inv_matrix,
                                             cov_matrix=None,
                                             debugging_dir=None,
                                             **chi2_kwargs)
        # print("Trying f:", f, "=", chi2)
        return chi2

    # hmm is this right? we scale the MC, but not Vyy on data?
    smeared_scaled_chi2 = partial(scaled_chi2,
                                  ref_hist=input_hist_gen_binning_bg_subtracted_signal,
                                  scaled_hist=folded_mc_truth_gen_binning,
                                  cov_inv_matrix=vyy_inv_gen_binning)
    smeared_fit = optimize.minimize_scalar(smeared_scaled_chi2, bounds=(0.1, 10.))
    folded_scale = smeared_fit.x
    print("folded scale:", folded_scale)
    folded_mc_truth_gen_binning *= folded_scale

    smeared_chi2, smeared_ndf, smeared_p = unfolder.calculate_chi2(one_hist=folded_mc_truth_gen_binning,
                                                                   # one_hist=unfolder.get_folded_mc_truth(),
                                                                   # other_hist=unfolder.input_hist_bg_subtracted,
                                                                   other_hist=input_hist_gen_binning_bg_subtracted_signal,
                                                                   # cov_inv_matrix=unfolder.get_vyy_inv_ndarray(),
                                                                   # cov_inv_matrix=unfolder.get_vyy_inv_no_bg_subtraction_ndarray(),
                                                                   cov_inv_matrix=vyy_inv_gen_binning,
                                                                   cov_matrix=vyy_gen_binning,
                                                                   debugging_dir=os.path.join(setup.output_dir, "smeared_chi2_gen_binning_signal"),
                                                                   **chi2_kwargs)
    print('smeared chi2, ndf, chi2/ndf, p:', smeared_chi2, smeared_ndf, smeared_chi2/smeared_ndf, smeared_p)


    # Do folding for Alt MC
    # --------------------------------------------------------------------------
    folded_alt_mc_truth_gen_binning = unfolder.get_ndarray_signal_region_no_overflow(
                                          unfolder.convert_reco_binned_hist_to_gen_binned(
                                            unfolder.fold_generator_level(setup.region['alt_hist_mc_gen'])
                                          ),
                                          xbinning='generator',
                                          ybinning=None
                                      )
    # folded_alt_mc_truth_gen_binning = unfolder.convert_reco_binned_hist_to_gen_binned(unfolder.fold_generator_level(setup.region['alt_hist_mc_gen']))

    if no_null_bins:
        folded_alt_mc_truth_gen_binning = remove_null_bins(folded_alt_mc_truth_gen_binning, null_bins)

    alt_smeared_scaled_chi2 = partial(scaled_chi2,
                                      ref_hist=input_hist_gen_binning_bg_subtracted_signal,
                                      scaled_hist=folded_alt_mc_truth_gen_binning,
                                      cov_inv_matrix=vyy_inv_gen_binning)
    alt_smeared_fit = optimize.minimize_scalar(alt_smeared_scaled_chi2, bounds=(0.1, 10.))
    alt_folded_scale = alt_smeared_fit.x
    print("alt folded scale:", alt_folded_scale)
    folded_alt_mc_truth_gen_binning *= alt_folded_scale

    # do for alt model
    smeared_alt_chi2, smeared_alt_ndf, smeared_alt_p = unfolder.calculate_chi2(one_hist=folded_alt_mc_truth_gen_binning,
                                                                               # one_hist=unfolder.fold_generator_level(setup.region['alt_hist_mc_gen']),
                                                                               # other_hist=unfolder.input_hist_bg_subtracted,
                                                                               other_hist=input_hist_gen_binning_bg_subtracted_signal,
                                                                               # cov_inv_matrix=unfolder.get_vyy_inv_ndarray(),
                                                                               # cov_inv_matrix=unfolder.get_vyy_inv_no_bg_subtraction_ndarray(),
                                                                               cov_inv_matrix=vyy_inv_gen_binning,
                                                                               debugging_dir=os.path.join(setup.output_dir, "smeared_alt_chi2_gen_binning_signal"),
                                                                               # debugging_dir=None,
                                                                               **chi2_kwargs)
    print('smeared alt chi2, ndf, chi2/ndf, p:', smeared_alt_chi2, smeared_alt_ndf, smeared_alt_chi2/smeared_alt_ndf, smeared_alt_p)


    # Do unfolding chi2 for nominal MC
    # --------------------------------------------------------------------------
    unfolded_signal = unfolder.get_ndarray_signal_region_no_overflow(
                          unfolder.unfolded,
                          xbinning='generator',
                          ybinning=None
                      )
    # print("unfolded_signal.shape:", unfolded_signal.shape)

    hist_truth_signal = unfolder.get_ndarray_signal_region_no_overflow(
                            unfolder.hist_truth,
                            xbinning='generator',
                            ybinning=None
                        )

    ematrix_input = unfolder.get_ematrix_input().Clone()
    ematrix_rsp = unfolder.get_ematrix_stat_response()
    ematrix_input.Add(ematrix_rsp)

    vxx_total_signal = unfolder.get_ndarray_signal_region_no_overflow(
                           ematrix_input,
                           # unfolder.get_vxx_ndarray(),  # tunfold's version
                           xbinning='generator',
                           ybinning='generator'
                       )
    if no_null_bins:
        null_bins = get_null_bins(vxx_total_signal)
        print('null bins unfolded space:', null_bins)
        unfolded_signal = remove_null_bins(unfolded_signal, null_bins)
        hist_truth_signal = remove_null_bins(hist_truth_signal, null_bins)
        vxx_total_signal = remove_null_bins(vxx_total_signal, null_bins)

    vxx_total_signal_inv = np.linalg.pinv(vxx_total_signal, 1E-200)  # no need to cut down to signal region

    # tunfold's version
    # vxx_total_signal_inv = unfolder.get_ndarray_signal_region_no_overflow(
    #                           unfolder.get_vxx_inv_ndarray(),
    #                           xbinning='generator',
    #                           ybinning='generator'
    #                       )

    # unfolded_scaled_chi2 = partial(scaled_chi2,
    #                                ref_hist=unfolded_signal,
    #                                scaled_hist=hist_truth_signal,
    #                                cov_inv_matrix=vxx_total_signal_inv)
    # unfolded_fit = optimize.minimize_scalar(unfolded_scaled_chi2, bounds=(0.1, 10.))
    # signal_scale = unfolded_fit.x
    signal_scale = folded_scale  #TODO: should this be the same as folded scale? Is v.close in practice, but not exactly the same?
    print("signal scale:", signal_scale)
    hist_truth_signal *= signal_scale

    unfolded_chi2, unfolded_ndf, unfolded_p = unfolder.calculate_chi2(one_hist=unfolded_signal,
                                                                      other_hist=hist_truth_signal,
                                                                      cov_inv_matrix=vxx_total_signal_inv,
                                                                      cov_matrix=vxx_total_signal,
                                                                      debugging_dir=os.path.join(setup.output_dir, "unfolded_chi2_signal_full_inv"),
                                                                      # debugging_dir=None,
                                                                      **chi2_kwargs)
    print('unfolded chi2, ndf, chi2/ndf, p:', unfolded_chi2, unfolded_ndf, unfolded_chi2/unfolded_ndf, unfolded_p)


    # Do unfolding chi2 for alt MC
    # --------------------------------------------------------------------------
    hist_alt_truth_signal = unfolder.get_ndarray_signal_region_no_overflow(
                            setup.region['alt_hist_mc_gen'],
                            xbinning='generator',
                            ybinning=None
                        )

    # unfolded_alt_scaled_chi2 = partial(scaled_chi2,
    #                                    ref_hist=unfolded_signal,
    #                                    scaled_hist=hist_alt_truth_signal,
    #                                    cov_inv_matrix=vxx_total_signal_inv)
    # unfolded_alt_fit = optimize.minimize_scalar(unfolded_alt_scaled_chi2, bounds=(0.1, 10.))
    # alt_signal_scale = unfolded_alt_fit.x

    if no_null_bins:
        hist_alt_truth_signal = remove_null_bins(hist_alt_truth_signal, null_bins)

    alt_signal_scale = alt_folded_scale
    print("alt signal scale:", alt_signal_scale)
    hist_alt_truth_signal *= alt_signal_scale

    unfolded_alt_chi2, unfolded_alt_ndf, unfolded_alt_p = unfolder.calculate_chi2(one_hist=unfolded_signal,
                                                                                  other_hist=hist_alt_truth_signal,
                                                                                  cov_inv_matrix=vxx_total_signal_inv,
                                                                                  debugging_dir=os.path.join(setup.output_dir, "unfolded_alt_chi2_signal_full_inv"),
                                                                                  # debugging_dir=None,
                                                                                  **chi2_kwargs)
    print('unfolded alt chi2, ndf, chi2/ndf, p:', unfolded_alt_chi2, unfolded_alt_ndf, unfolded_alt_chi2/unfolded_alt_ndf, unfolded_alt_p)

    return {
                "region": setup.region['label'],
                "is_groomed": "groomed" in setup.region['name'],
                "angle": setup.angle.var,

                "smeared_chi2": smeared_chi2,
                "smeared_ndf": smeared_ndf,
                "smeared_chi2_ov_ndf": smeared_chi2/smeared_ndf,
                "smeared_p": smeared_p,

                "smeared_alt_chi2": smeared_alt_chi2,
                "smeared_alt_ndf": smeared_alt_ndf,
                "smeared_alt_chi2_ov_ndf": smeared_alt_chi2/smeared_alt_ndf,
                "smeared_alt_p": smeared_alt_p,

                "unfolded_chi2": unfolded_chi2,
                "unfolded_ndf": unfolded_ndf,
                "unfolded_chi2_ov_ndf": unfolded_chi2/unfolded_ndf,
                "unfolded_p": unfolded_p,

                "unfolded_alt_chi2": unfolded_alt_chi2,
                "unfolded_alt_ndf": unfolded_alt_ndf,
                "unfolded_alt_chi2_ov_ndf": unfolded_alt_chi2/unfolded_alt_ndf,
                "unfolded_alt_p": unfolded_alt_p,
            }


def print_chi2_latex_table(df):
    df_sorted = df.sort_values(by=['region', 'angle'])
    # print(r'\begin{tabular}{c|c|c|c|c|c|c|c|c|c|c|c}')
    print(r'\begin{tabular}{c|c|c|c|c|c|c|c}')
    print(("Region & "
           "Angle & "
           r"$\chi_{\text{unfolded}}^2 / N_{DoF}$ & "
           # r"$p_{\text{unfolded}}$ & "
           r"$\chi_{\text{unfolded, alt. model}}^2 / N_{DoF}$ & "
           # r"$p_{\text{unfolded, alt. model}}$ & "
           r"$\chi_{\text{smeared}}^2 / N_{DoF}$ & "
           # r"$p_{\text{smeared}}$ & "
           r"$\chi_{\text{smeared, alt. model}}^2 / N_{DoF}$ & "
           # r"$p_{\text{smeared, alt. model}}$ & "
           r"$\chi_{\text{unfolded}}^2 < \chi_{\text{smeared}}^2$ & "
           r"$\chi_{\text{unfolded, alt. model}}^2 < \chi_{\text{smeared, alt. model}}^2$ \\"))
    print(r"\hline \\")

    angle_lookup = {a.var: a.mathmode for a in qgc.COMMON_VARS}
    for row in df_sorted.itertuples():
        angle_name = row.angle.replace("jet_", "").replace("_charged", " (charged)")
        angle_name = "$"+angle_lookup[row.angle]+"$"
        if row.is_groomed:
            angle_name = "Groomed " + angle_name

        # print(dir(row))
        print(("{region} & "
               "{angle_name} & "
               "{unfolded_chi2:.3e} / {unfolded_ndf} = {unfolded_chi2_ov_ndf:.3g} & "
               # "{unfolded_p:.3e} & "
               "{unfolded_alt_chi2:.3e} / {unfolded_alt_ndf} = {unfolded_alt_chi2_ov_ndf:.3g} & "
               # "{unfolded_alt_p:.3e} & "
               "{smeared_chi2:.3e} / {smeared_ndf} = {smeared_chi2_ov_ndf:.3g} & "
               # "{smeared_p:.3e} & "
               "{smeared_alt_chi2:.3e} / {smeared_alt_ndf} = {smeared_alt_chi2_ov_ndf:.3g} & "
               # "{smeared_alt_p:.3e} & "
               "{result} & "
               "{result_alt} \\\\").format(angle_name=angle_name,
                                           result=r"$\checkmark$" if row.unfolded_chi2_ov_ndf < row.smeared_chi2_ov_ndf else r"$\times$",
                                           result_alt=r"$\checkmark$" if row.unfolded_alt_chi2_ov_ndf < row.smeared_alt_chi2_ov_ndf else r"$\times$",
                                           **row._asdict()))
    print(r'\end{tabular}')


def print_chi2_user_table(df):
    """Print chi2 results to screen as quick pass/fail table"""
    df_sorted = df.sort_values(by=['region', 'angle'])
    df_sorted['pass nominal?'] = df_sorted['unfolded_chi2_ov_ndf'] < df_sorted['smeared_chi2_ov_ndf']
    df_sorted['pass alt.?'] = df_sorted['unfolded_alt_chi2_ov_ndf'] < df_sorted['smeared_alt_chi2_ov_ndf']
    # use emoji for visual speed & clarity
    check = u'\u2705'
    cross = u'\u274C'
    danger = u'\u2757'
    df_sorted['pass nominal?'].replace({True: check, False: cross}, inplace=True)
    df_sorted['pass alt.?'].replace({True: check, False: cross}, inplace=True)
    bad_nominal = df_sorted['unfolded_chi2'] < 0
    df_sorted['pass nominal?'][bad_nominal] =  danger
    bad_alt = df_sorted['unfolded_alt_chi2'] < 0
    df_sorted['pass alt.?'][bad_alt] =  danger
    cols = ['region', 'is_groomed', 'angle', 'pass nominal?', 'pass alt.?']
    print(df_sorted.to_string(index=False, columns=cols))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",
                        help="Source directory (should be the one made by unfolding.py)")
    parser.add_argument("--angles",
                        choices=list(qgc.VAR_UNFOLD_DICT['ungroomed'].keys()) + ["all"],
                        nargs='+',
                        help="Lambda angles to unfold, or 'all' for all of them (default).",
                        default=["all"])

    parser.add_argument("--noNullBins",
                        help='Remove null bins - this is important for inversion!',
                        action='store_true')

    region_group = parser.add_argument_group('Region selection')
    region_group.add_argument("--doAllRegions",
                              action='store_true',
                              help='Do unfolding for all regions (dijet, Z+J, groomed, ungroomed)')
    region_group.add_argument("--doDijetCentral",
                              action='store_true',
                              help='Do unfolding for dijet (central) jets')
    region_group.add_argument("--doDijetForward",
                              action='store_true',
                              help='Do unfolding for dijet (forward) jets')
    region_group.add_argument("--doDijetCentralGroomed",
                              action='store_true',
                              help='Do unfolding for groomed dijet (central) jets')
    region_group.add_argument("--doDijetForwardGroomed",
                              action='store_true',
                              help='Do unfolding for groomed dijet (forward) jets')
    region_group.add_argument("--doZPJ",
                              action='store_true',
                              help='Do unfolding for Z+jet jets')
    region_group.add_argument("--doZPJGroomed",
                              action='store_true',
                              help='Do unfolding for groomed Z+jet jets')

    args = parser.parse_args()

    if args.doAllRegions:
        for x in ['doDijetCentral', 'doDijetForward', 'doDijetCentralGroomed', 'doDijetForwardGroomed', 'doZPJ', 'doZPJGroomed']:
            setattr(args, x, True)

    regions = setup_regions_from_argparse(args)
    if len(regions) == 0:
        raise RuntimeError("Need at least 1 region")

    if args.angles[0] == "all":
        angles = qgc.COMMON_VARS
    else:
        angles = [a for a in qgc.COMMON_VARS if a.var in args.angles]

    if len(angles) == 0:
        raise RuntimeError("Need at least 1 angle")

    # jet_algo = "AK4 PF PUPPI"
    jet_algo = "AK4"
    if "ak8puppi" in args.source:
        # jet_algo = "AK8 PF PUPPI"
        jet_algo = "AK8"

    has_data = not ('_MC_all' in args.source or '_MC_split' in args.source)

    all_chi2_stats = []

    num_all_iterations = len(regions) * len(angles)
    counter = 1

    # Iterate through regions & variables
    for region in regions:

        region_dir = os.path.join(args.source, region['name'])
        if not os.path.isdir(region_dir):
            print("! Warning ! cannot find region dir", region_dir, '- skipping region')
            continue

        for angle in angles:
            append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name
            print("*"*120)
            print("Algo/Region/Var: %s %s %s (%d/%d)" % (jet_algo, region['name'], angle.var, counter, num_all_iterations))
            print("*"*120)
            counter += 1

            angle_output_dir = "%s/%s/%s" % (args.source, region['name'], angle.var)
            if not os.path.isdir(angle_output_dir):
                print("! Warning ! cannot find angle dir", angle_output_dir, '- skipping angle', angle.var)
                continue

            # Get region dict from pickle file
            pickle_filename = os.path.join(angle_output_dir, "unfolding_result.pkl")
            unpickled_region = unpickle_region(pickle_filename)
            if unpickled_region is None:
                continue

            # check
            if region['name'] != unpickled_region['name']:
                raise RuntimeError("Mismatch region name")

            # use region dict from unpickling
            # don't use update(), mega slow
            this_region = unpickled_region
            this_region['label'] = region['label']

            # Calc the stats for this configuration
            # ------------------------------------------------------------------
            setup = Setup(jet_algo=jet_algo,
                          region=this_region,
                          angle=angle,
                          output_dir=angle_output_dir,
                          has_data=has_data)

            all_chi2_stats.append(get_bottom_line_stats(setup, no_null_bins=args.noNullBins))

            # cleanup object, as it uses loads of memory
            if num_all_iterations > 1:
                print("...tidying up...")
                del unpickled_region
                del setup


    if len(all_chi2_stats) > 0:
        df_stats = pd.DataFrame(all_chi2_stats)
        df_stats['region'] = df_stats['region'].astype('category')
        df_stats['angle'] = df_stats['angle'].astype('category')
        print(df_stats.head())
        print("")
        print_chi2_user_table(df_stats)
        print("")
        print_chi2_latex_table(df_stats)

if __name__ == "__main__":
    main()
