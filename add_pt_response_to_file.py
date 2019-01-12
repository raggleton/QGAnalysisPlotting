#!/usr/bin/env python

"""Extract pT response hist from one ROOT file, and add it to another ROOT file

Needed to construct one response ROOT file for use in RIVET.
"""

import argparse
import os


import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
# ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)

# My stuff
import common_utils as cu


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", help="Input ROOT file with response hist", required=True)
    parser.add_argument("--inputDir", help="TDirectory in ROOT file", required=True)
    parser.add_argument("--output", help="Output ROOT file to append hist to", required=True)
    args = parser.parse_args()

    hist = cu.grab_obj_from_file(args.input, "%s/pt_jet_response" % (args.inputDir))
    hist.SetName("pt_jet_response")
    outf = cu.open_root_file(args.output, "UPDATE")
    hist.Write()
    outf.Close()
