#!/usr/bin/env python

"""Make 1D lambda hists & rebin according to custom rebinning. 
Also set names ready for RIVET/YODA.
"""


import ROOT
import os
import argparse
from array import array

# My stuff
import qg_common as qgc
import qg_general_plots as qgp
import common_utils as cu

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)


def create_yoda_hist_name(region, angle_name, pt_bin):
    """Produce YODA hist name"""
    region_num = 1 if "dijet" in region.lower() else 2
    angle_names = [a.var for a in qgc.COMMON_VARS]
    angle_num = angle_names.index(angle_name)
    pt_bin_num = qgc.PT_BINS.index(pt_bin)
    return "d{:02d}-x{:02d}-y{:02d}".format(region_num, angle_num, pt_bin_num)


def make_1D_rebin_hists(input_filename, plot_dirname, output_filename):
    
    # QG variable plots

    in_f = cu.open_root_file(input_filename)
    out_f = ROOT.TFile(output_filename, "RECREATE")

    for ang_ind, ang in enumerate(qgc.COMMON_VARS):
        if ang.var not in qgc.ANGLE_REBIN_DICT:
            continue

        in_f.cd()

        var_prepend = ""
        obj_name = "%s%s_vs_pt" % (var_prepend, ang.var)
        h2d = cu.get_from_tfile(in_f, "%s/%s" % (plot_dirname, obj_name))

        for pt_ind, (start_val, end_val) in enumerate(qgc.PT_BINS):
            hist = qgp.get_projection_plot(h2d, start_val, end_val)
            this_rebins = qgc.ANGLE_REBIN_DICT[ang.var]
            # new_name = "%s_Pt%sto%d" % (ang.var, start_val, end_val)
            new_name = create_yoda_hist_name(plot_dirname, ang.var, (start_val, end_val))
            rebin_hist = hist.Rebin(len(this_rebins)-1, new_name, array('d', this_rebins))
            out_f.WriteTObject(rebin_hist)  # saves faffing with cd()

    out_f.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input file", required=True)
    parser.add_argument("--plotDirname", help="Plot directory name in input file", required=True)
    parser.add_argument("--output", help="Output file", required=True)
    args = parser.parse_args()

    make_1D_rebin_hists(args.input, args.plotDirname, args.output)
