#!/usr/bin/env python

"""Print main QG lambda plots, comparing samples.

Does comparison plots for data & MC, separately.

Make sure jetht and zerobias are hadded into same file
"""

from __future__ import print_function

import os
os.nice(10)
import ROOT
from MyStyle import My_Style
My_Style.cd()
import argparse
from copy import deepcopy

# My stuff
from comparator import Contribution, Plot
import qg_common as qgc
import qg_general_plots as qgp
import common_utils as cu

# For debugging
import sys
import tracers
# sys.settrace(tracers.trace_calls)
# sys.settrace(tracers.trace_calls_detail)
# sys.settrace(tracers.trace_calls_and_returns)

# import hunter
# hunter.trace(module='comparator')

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


# Control plot output format
OUTPUT_FMT = "pdf"

TOTAL_LUMI = 35918


def exists_in_all(filename, dirs):
    return all([os.path.isfile(os.path.join(d, filename)) for d in dirs])


def do_comparison_plots(workdir_label_pairs, output_dir):
    dirnames = [w[0] for w in workdir_label_pairs]

    # templates, we'll change the filename/dir as per instance
    total_len = len(workdir_label_pairs)
    mark = cu.Marker()
    sources = [
       {
           "root_dir": wd,
           'label': label,
           "style": {
               'line_style': 1,
               'line_color': cu.get_colour_seq(ind, total_len),
               'marker_color': cu.get_colour_seq(ind, total_len),
               'marker_style': m,
               'marker_size': 0.75,
           }

       }
       for ind, ((wd, label), m) in enumerate(zip(workdir_label_pairs, mark.cycle(cycle_filling=True)))
    ]
    # print(sources)

    # COMPARE NOMINAL QCD
    if exists_in_all(qgc.QCD_FILENAME, dirnames):
        # qgp.do_all_exclusive_plots_comparison(sources,
        #                                       var_list=qgc.COMMON_VARS,
        #                                       pt_bins=qgc.PT_BINS,
        #                                       qcd_filename=qgc.QCD_FILENAME,
        #                                       dj_cen_dirname="Dijet_QG_central_tighter",
        #                                       dj_fwd_dirname=None,
        #                                       zpj_dirname=None,
        #                                       plot_dir=os.path.join(output_dir, "plots_qcd_compare_dijet_central"),
        #                                       subplot_type="ratio", # will use the 1st entry by default
        #                                       subplot_title="* / %s" % (sources[0]['label']),
        #                                       do_flav_tagged=False,
        #                                       has_data=False)

        # qgp.do_all_exclusive_plots_comparison(sources,
        #                                       var_list=qgc.COMMON_VARS,
        #                                       pt_bins=qgc.PT_BINS,
        #                                       qcd_filename=qgc.QCD_FILENAME,
        #                                       dj_cen_dirname=None,
        #                                       dj_fwd_dirname="Dijet_QG_forward_tighter",
        #                                       zpj_dirname=None,
        #                                       plot_dir=os.path.join(output_dir, "plots_qcd_compare_dijet_forward"),
        #                                       subplot_type="ratio", # will use the 1st entry by default
        #                                       subplot_title="* / %s" % (sources[0]['label']),
        #                                       do_flav_tagged=False,
        #                                       has_data=False)

        # do both dijet jets together
        # add gen level
        gen_qcd_sources = [deepcopy(sources[0])] # only 1st one for now
        for gd in gen_qcd_sources:
            gd.update({
                'dj_fwd_dirname': 'Dijet_QG_gen',
                'label': gd['label'] + ' [Gen]'
            })
            style_dict = gd['style']
            style_dict.update({'line_style': 2})
            style_dict.update({'line_color': style_dict['line_color']+5})
            style_dict.update({'marker_color': style_dict['marker_color']+5})
            style_dict.update({'marker_style': style_dict['marker_style']+1})

        qgp.do_all_exclusive_plots_comparison(sources + gen_qcd_sources,
                                              var_list=qgc.COMMON_VARS,
                                              pt_bins=qgc.PT_BINS,
                                              qcd_filename=qgc.QCD_FILENAME,
                                              dj_cen_dirname=None,
                                              dj_fwd_dirname="Dijet_QG_tighter",
                                              zpj_dirname=None,
                                              plot_dir=os.path.join(output_dir, "plots_qcd_compare_dijet"),
                                              subplot_type="ratio", # will use the 1st entry by default
                                              subplot_title="* / %s" % (sources[0]['label']),
                                              do_flav_tagged=False,
                                              has_data=False,
                                              title=qgc.Dijet_LABEL,
                                              show_region_labels=False)

    # COMPARE NOMINAL DY
    if exists_in_all(qgc.DY_FILENAME, dirnames):
        qgp.do_all_exclusive_plots_comparison(sources,
                                              var_list=qgc.COMMON_VARS,
                                              pt_bins=qgc.PT_BINS,
                                              qcd_filename=qgc.DY_FILENAME,
                                              dj_cen_dirname=None,
                                              dj_fwd_dirname=None,
                                              zpj_dirname="ZPlusJets_QG",
                                              plot_dir=os.path.join(output_dir, "plots_dy_compare"),
                                              subplot_type="ratio", # will use the 1st entry by default
                                              subplot_title="* / %s" % (sources[0]['label']),
                                              do_flav_tagged=False,
                                              has_data=False,
                                              title=qgc.ZpJ_LABEL,
                                              show_region_labels=False)


    # COMPARE JETHT+ZEROBIAS
    if exists_in_all(qgc.JETHT_ZB_FILENAME, dirnames):
        qgp.do_all_exclusive_plots_comparison(sources,
                                              var_list=qgc.COMMON_VARS,
                                              pt_bins=qgc.PT_BINS,
                                              qcd_filename=qgc.JETHT_ZB_FILENAME,
                                              dj_cen_dirname="Dijet_QG_central_tighter",
                                              dj_fwd_dirname=None,
                                              zpj_dirname=None,
                                              plot_dir=os.path.join(output_dir, "plots_jetht_zb_compare_dijet_central"),
                                              subplot_type="ratio", # will use the 1st entry by default
                                              subplot_title="* / %s" % (sources[0]['label']),
                                              do_flav_tagged=False,
                                              has_data=True)

        qgp.do_all_exclusive_plots_comparison(sources,
                                              var_list=qgc.COMMON_VARS,
                                              pt_bins=qgc.PT_BINS,
                                              qcd_filename=qgc.JETHT_ZB_FILENAME,
                                              dj_cen_dirname=None,
                                              dj_fwd_dirname="Dijet_QG_forward_tighter",
                                              zpj_dirname=None,
                                              plot_dir=os.path.join(output_dir, "plots_jetht_zb_compare_dijet_forward"),
                                              subplot_type="ratio", # will use the 1st entry by default
                                              subplot_title="* / %s" % (sources[0]['label']),
                                              do_flav_tagged=False,
                                              has_data=True)

    # COMPARE SINGLEMU

    # COMPARE HERWIG++ QCD
    if exists_in_all(qgc.QCD_HERWIG_FILENAME, dirnames):
        qgp.do_all_exclusive_plots_comparison(sources,
                                              var_list=qgc.COMMON_VARS,
                                              pt_bins=qgc.PT_BINS,
                                              qcd_filename=qgc.QCD_HERWIG_FILENAME,
                                              dj_cen_dirname="Dijet_QG_central_tighter",
                                              dj_fwd_dirname=None,
                                              zpj_dirname=None,
                                              plot_dir=os.path.join(output_dir, "plots_qcd_herwig_compare_dijet_central"),
                                              subplot_type="ratio", # will use the 1st entry by default
                                              subplot_title="* / %s" % (sources[0]['label']),
                                              subplot_limits=[0.9, 1.1],
                                              do_flav_tagged=False,
                                              has_data=False)

        qgp.do_all_exclusive_plots_comparison(sources,
                                              var_list=qgc.COMMON_VARS,
                                              pt_bins=qgc.PT_BINS,
                                              qcd_filename=qgc.QCD_HERWIG_FILENAME,
                                              dj_cen_dirname=None,
                                              dj_fwd_dirname="Dijet_QG_forward_tighter",
                                              zpj_dirname=None,
                                              plot_dir=os.path.join(output_dir, "plots_qcd_herwig_compare_dijet_forward"),
                                              subplot_type="ratio", # will use the 1st entry by default
                                              subplot_title="* / %s" % (sources[0]['label']),
                                              subplot_limits=[0.9, 1.1],
                                              do_flav_tagged=False,
                                              has_data=False)

    # COMPARE HERWIG++ DY


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',
                        required=True,
                        action='append',
                        help='Workdir with ROOT files to process. '
                        'Each --workdir must have a corresponding --label')
    parser.add_argument('--label',
                        required=True,
                        action='append',
                        help='Label for workdir.')
    parser.add_argument("-o", "--output",
                        help="Directory to put output plot dirs into. Default is workdir.",
                        default=None)
    # parser.add_argument("--doExperimentalSysts",
    #                     help="Also draw experimental systematic variations",
    #                     action='store_true')
    # parser.add_argument("--doModelSysts",
    #                     help="Also draw model systematic variations",
    #                     action='store_true')

    args = parser.parse_args()
    print(args)

    if args.output is None:
        args.output = args.workdir[0]
        print("No output dir specified, using", args.output)

    do_comparison_plots(list(zip(args.workdir, args.label)), output_dir=args.output)
