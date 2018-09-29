#!/usr/bin/env python

"""
Add 2 files together, applying scaling to each histogram.

Designed for adding ZeroBias + JetHT datasets 
(ZB requires scaling as prescaled, and JetHT prescaling handled in main analysis)
"""

from __future__ import print_function
import os
import sys
import ROOT
import argparse


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()


def get_list_of_objects(tobj):
    return sorted([x.GetName() for x in tobj.GetListOfKeys()])


def merge_zb_jetht(zb_filename, jetht_filename, output_filename, factor=1.):
    if not os.path.isfile(zb_filename):
        raise IOError("Cannot find %s" % zb_filename)
    if not os.path.isfile(jetht_filename):
        raise IOError("Cannot find %s" % jetht_filename)
    
    zb_file = ROOT.TFile(zb_filename)
    jetht_file = ROOT.TFile(jetht_filename)
    out_file = ROOT.TFile(output_filename, 'RECREATE')
    for f in [zb_file, jetht_file, out_file]:
        if f.IsZombie():
            raise IOError("Couldn't open ROOT file %s" % f.GetName)

    interesting_types = [ROOT.TH1D, ROOT.TH1F, ROOT.TH1I, 
                         ROOT.TH2D, ROOT.TH2F, ROOT.TH2I]

    # Only iterate over directories common to both files
    zb_dirs = set(get_list_of_objects(zb_file))
    jetht_dirs = set(get_list_of_objects(jetht_file))
    common_dirs = sorted(list(zb_dirs & jetht_dirs))
    # print(common_dirs)

    common_dirs_new = []
    for dirname in common_dirs:
        zb = zb_file.Get(dirname)
        if isinstance(zb, ROOT.TDirectory):
            common_dirs_new.append(dirname)
            # skip any top level hists for now
    
    common_dirs = common_dirs_new
    # Iterate over each dir, and in each iterate over each common hist
    for dirname in common_dirs:
        out_dir = out_file.mkdir(dirname)

        zb_tdir = zb_file.Get(dirname)
        jetht_tdir = jetht_file.Get(dirname)        

        zb_hists = set(get_list_of_objects(zb_tdir))
        jetht_hists = set(get_list_of_objects(jetht_tdir))
        common_hists = sorted(list(zb_hists & jetht_hists))
        # print(common_hists)

        for hname in common_hists:
            zb_hist = zb_tdir.Get(hname)
            jetht_hist = jetht_tdir.Get(hname)

            if not zb_hist:
                raise IOError("Couldn't get hist %s from ZB file" % os.path.join(dirname, hname))
            if not jetht_hist:
                raise IOError("Couldn't get hist %s from JetHT file" % os.path.join(dirname, hname))

            if not any(isinstance(zb_hist, t) for t in interesting_types):
                continue

            zb_clone = zb_hist.Clone()
            zb_clone.Sumw2()
            jetht_clone = jetht_hist.Clone()
            jetht_clone.Sumw2()
            jetht_clone.Add(zb_clone, factor)
            # print("Writing", jetht_clone.GetName())
            out_dir.cd()
            jetht_clone.Write()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--zb", required=True)
    parser.add_argument("--jetht", required=True)
    parser.add_argument("--output", required=True)
    ZB_SCALE_FACTOR = 35918219492.947 / 29048.362
    parser.add_argument("--zbFactor", type=float, default=ZB_SCALE_FACTOR)
    args = parser.parse_args()

    merge_zb_jetht(zb_filename=args.zb, jetht_filename=args.jetht, 
                   output_filename=args.output, factor=args.zbFactor)
