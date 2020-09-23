#!/usr/bin/env python

"""Make 1D lambda hists & rebin according to custom rebinning.
Also make 1D migration summary plots.

Uses txt file from determine_lambda_binning.py, or python file (qg_common.py)
"""


import argparse
from MyStyle import My_Style
My_Style.cd()
import os
os.nice(10)
import sys
from itertools import product
from array import array
from math import sqrt
import bisect
import numpy as np
from collections import OrderedDict
from copy import deepcopy
from ast import literal_eval

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
# ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# My stuff
import qg_general_plots as qgg
import qg_common as qgc
from comparator import Contribution, Plot
import common_utils as cu
from determine_lambda_binning import rebin_2d_hist, make_rebinned_2d_hist, make_plots, get_summed_hists

ROOT.gStyle.SetPaintTextFormat(".3f")

# Control output format
OUTPUT_FMT = "pdf"


def load_binning_dict(binning_file):
    """Load binning schema from txt file into dict"""
    if not os.path.isfile(binning_file):
        raise IOError("Cannot find binning file %s" % binning_file)

    binning_dict = OrderedDict()

    with open(binning_file) as f:
        for line in f.readlines():
            this_line = line.strip()
            if not this_line:
                continue
            name, bins = line.split(":", 1)
            name = name.strip()
            bins = literal_eval(bins.strip()) # vector of individual bins
            # convert to pairs of bin edges
            bin_pairs = [[x,y] for x, y in zip(bins[:-1], bins[1:])]
            binning_dict[name] = bin_pairs

    return binning_dict


def get_reference_pt_region(var_name):
    """Get reference pt region

    In a function cos it may get complex
    """
    var_name_lower = var_name.lower()
    if "multiplicity" in var_name_lower:
        return "_highPt"

    return "_midPt"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',
                        help='Input ROOT file to rebin.')
    parser.add_argument('--binningFile',
                        help='Txt file with binning. Default is binning in qg_common.py')
    parser.add_argument("--outputDir",
                        help="Directory to put output plot dirs into",
                        default=None)
    parser.add_argument("--outputFile",
                        help="Output ROOT file for rebinned 2D hists",
                        default=None)
    ref_choices = ["groomed", "ungroomed"]
    parser.add_argument("--ungroomedRef",
                        help="Reference region for ungroomed vars",
                        choices=ref_choices,
                        default="ungroomed")
    parser.add_argument("--groomedRef",
                        help="Reference region for groomed vars",
                        choices=ref_choices,
                        default="groomed")
    parser.add_argument("--target",
                        help="Target purity/stability to put on plot",
                        default=None)
    args = parser.parse_args()
    print(args)

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
        source_plot_dir_names = ["Dijet_QG_central_tighter", "Dijet_QG_forward_tighter", 
                                 "Dijet_QG_central_tighter_groomed", "Dijet_QG_forward_tighter_groomed", 
                                 ("Dijet_QG_central_tighter", "Dijet_QG_forward_tighter"), # do sum over these regions
                                 ("Dijet_QG_central_tighter_groomed", "Dijet_QG_forward_tighter_groomed")][-2:]
        region_labels = [qgc.QCD_Dijet_CEN_LABEL, qgc.QCD_Dijet_FWD_LABEL, 
                         qgc.QCD_Dijet_CEN_GROOMED_LABEL, qgc.QCD_Dijet_FWD_GROOMED_LABEL, 
                         qgc.Dijet_LABEL,
                         qgc.Dijet_GROOMED_LABEL][-2:]
    elif "dyjetstoll" in args.input.lower():
        source_plot_dir_names = ["ZPlusJets_QG", "ZPlusJets_QG_groomed"]
        region_labels = [qgc.DY_ZpJ_LABEL, qgc.DY_ZpJ_GROOMED_LABEL]
    else:
        raise RuntimeError("No idea which region we're using")

    use_txt_binning = False
    if args.binningFile is not None:
        use_txt_binning = True
        binning_dict = load_binning_dict(args.binningFile)

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

    input_tfile = cu.open_root_file(args.input)

    # Iterate over all variables, regions, groomed/ungroomed, pt regions
    for angle, (source_plot_dir_name, region_label), pt_region_dict in product(qgc.COMMON_VARS[:],
                                                                               zip(source_plot_dir_names, region_labels),
                                                                               pt_regions):
        multi_region = not isinstance(source_plot_dir_name, str)
        if multi_region:
            var_prepend = "groomed " if "groomed" in source_plot_dir_name[0] else ""
            total_plot_dir_name = "+".join(source_plot_dir_name)
        else:
            var_prepend = "groomed " if "groomed" in source_plot_dir_name else ""
            total_plot_dir_name = source_plot_dir_name
        print(angle, source_plot_dir_name, region_label, pt_region_dict['title'])

        var_dict = {
            "name": "%s/%s%s" % (total_plot_dir_name, angle.var, pt_region_dict['append']),
            "var_label": "%s%s (%s)" % (var_prepend, angle.name, angle.lambda_str),
            "title": "", # setup later
        }
        if multi_region:
            names = ["%s/%s%s_response" % (sname, angle.var, pt_region_dict['append']) for sname in source_plot_dir_name]
            h2d_orig = get_summed_hists(input_tfile, names)
        else:
            h2d_orig = cu.get_from_tfile(input_tfile, var_dict['name'] + "_response")

        # Get desired binning
        # -------------------
        if use_txt_binning:
            # This has all the logic for deciding what is the "reference" region
            key = var_dict['name']
            print(key)

            # Use one pt region for all pt regions
            ref_pt_region = get_reference_pt_region(var_dict['name'])
            key = key.replace(pt_region_dict['append'], ref_pt_region)
            print(key)

            ref_pt_region_dict = [x for x in pt_regions if x['append'] == ref_pt_region][0]
            var_dict["title"] = "%s\n%s\nRebinned for %s (central+forward)" % (region_label, pt_region_dict['title'], ref_pt_region_dict['title'])

            # Use same for central/forward
            # Use Dijet ones for Z+Jets
            ref_region = "Dijet_QG_central_tighter+Dijet_QG_forward_tighter"

            # Here we want separate groomed/ungroomed binnings
            if "groomed" in total_plot_dir_name:
                if args.groomedRef == "groomed":
                    # Use groomed for groomed
                    ref_region = ref_region.replace("_tighter", "_tighter_groomed")
                    var_dict['title'] += " (groomed)"
                else:
                    # Use ungroomed for groomed
                    var_dict['title'] += " (ungroomed)"
            else:
                if args.ungroomedRef == "groomed":
                    # Use groomed for ungroomed as well
                    ref_region = ref_region.replace("_tighter", "_tighter_groomed")
                    var_dict['title'] += " (groomed)"
                else:
                    # Use ungroomed for ungroomed
                    var_dict['title'] += " (ungroomed)"

            if args.target:
                var_dict['title'] += " (target %s)" % (args.target)

            key = key.replace(total_plot_dir_name, ref_region)
            print(key)
            new_binning = binning_dict[key]
            print(new_binning)

        else:
            # Get common binning scheme from qg_common
            var_binning_key = angle.var
            new_binning = qgc.VAR_UNFOLD_DICT[var_binning_key]['gen']
            # convert to pairs of bin edges
            new_binning = list(zip(new_binning[:-1], new_binning[1:]))
            print(new_binning)

            var_dict["title"] = "%s\n%s\nRebinned using existing binning" % (region_label, pt_region_dict['title'])
            if args.target:
                var_dict['title'] += " (target %s)" % (args.target)

        # Rebin 2D heatmap
        # ----------------
        h2d_rebin = make_rebinned_2d_hist(h2d_orig, new_binning, use_half_width_y=False)

        # Plot with new binning
        # ---------------------
        make_plots(h2d_rebin,
                   var_dict,
                   plot_dir=plot_dir,
                   append="rebinned_for%s" % (ref_pt_region) if use_txt_binning else "rebinned",
                   plot_migrations=True,
                   plot_reco=True,
                   plot_gen=True)

    input_tfile.Close()
