#!/usr/bin/env python

"""
Simple unpickling script to test different compression algorithms

Runs unpickling several times.
Profile by running with cProfile:

python -m cProfile -o unpickling.prof time_unpickling.py

Then you can use pstats to analyze it, e.g.:

import pstats
from pstats import SortKey
p = pstats.Stats('unpickling.prof')
p.sort_stats(SortKey.TIME).print_stats(10)

"""


from __future__ import print_function, division
import os
os.nice(10)
import gc
import lzma
import pickle
import gzip

import ROOT
from my_unfolder import MyUnfolder, pickle_region, unpickle_region
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)  # VERY IMPORTANT - somewhere, closing a TFile for exp systs deletes a map...dunno why
ROOT.gErrorIgnoreLevel = ROOT.kError


def unpickle_lzma(filename):
    with lzma.open(filename, 'r') as f:
        unpickled_region = pickle.load(f)
    return unpickled_region

def unpickle_lzma4(filename):
    with lzma.open(filename, 'r') as f:
        unpickled_region = pickle.load(f)
    return unpickled_region

def unpickle_gzip(filename):
    with gzip.open(filename, 'r') as f:
        unpickled_region = pickle.load(f)
    return unpickled_region



if __name__ == "__main__":
    lzma_filename = "dummyLZMA/unfolding_regularizeNone_MC_all_subFakes_densityModeBinWidth_constraintNone_pdfSystNoExperimentalSyst_signalRegionOnly/Dijet_central/jet_puppiMultiplicity/unfolding_result.pkl"
    for i in range(20):
        unpickle_lzma(lzma_filename)

    lzma_filename = "dummyLZMA4/unfolding_regularizeNone_MC_all_subFakes_densityModeBinWidth_constraintNone_pdfSystNoExperimentalSyst_signalRegionOnly/Dijet_central/jet_puppiMultiplicity/unfolding_result.pkl"
    for i in range(20):
        unpickle_lzma4(lzma_filename)

    gzip_filename = "dummyGZIP/unfolding_regularizeNone_MC_all_subFakes_densityModeBinWidth_constraintNone_pdfSystNoExperimentalSyst_signalRegionOnly/Dijet_central/jet_puppiMultiplicity/unfolding_result.pkl"
    for i in range(20):
        unpickle_gzip(gzip_filename)
