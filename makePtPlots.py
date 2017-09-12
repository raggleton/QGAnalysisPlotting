#!/usr/bin/env python

"""Print basic pT plots comparing 2 samples, here called pythia and herwig,
but could of course be anything.
"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
from comparator import Contribution, Plot, grab_obj
import os
from itertools import product

# For debugging
import sys
import tracers
# sys.settrace(tracers.trace_calls)
# sys.settrace(tracers.trace_calls_detail)
# sys.settrace(tracers.trace_calls_and_returns)


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)


def do_pt_plot(pythia_dir, herwig_dir, selection, hist_name, output_name, title=""):
    h_pythia = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_%s_.root" % (pythia_dir, selection), hist_name)
    h_herwig = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_%s_.root" % (herwig_dir, selection), hist_name)

    c_pythia = Contribution(h_pythia, label="MG+Pythia", line_color=ROOT.kBlue, marker_color=ROOT.kBlue, fill_color=ROOT.kBlue, normalise_hist=True)
    # c_herwig = Contribution(h_herwig, label="Herwig", line_color=ROOT.kRed, normalise_hist=True)
    c_herwig = Contribution(h_herwig, label="Pythia only", line_color=ROOT.kRed, marker_color=ROOT.kRed, fill_color=ROOT.kRed, normalise_hist=True)

    p = Plot([c_pythia, c_herwig], what="hist", legend=True, subplot_type='ratio', subplot=c_pythia, title=title)
    p.plot("NOSTACK HISTE")
    p.main_pad.SetLogy()
    p.container.SetMinimum(1E-12)
    p.subplot_container.SetMaximum(1.5)
    p.subplot_container.SetMinimum(0)
    p.canvas.Update()
    p.save(output_name)


if __name__ == "__main__":

    ALGOS = ["ak4", "ak8"]
    PUS = ["chs", "puppi"]

    for algo, pu, in product(ALGOS, PUS):
        SETUP = algo + pu

        PYTHIA_DIR = "workdir_%s_mgpythia" % (SETUP)
        HERWIG_DIR = "workdir_%s_herwig" % (SETUP)
        HERWIG_DIR = "workdir_%s_herwig_reweight" % (SETUP)
        # HERWIG_DIR = "workdir_%s_pythiaOnlyFlat" % (SETUP)

        odir = HERWIG_DIR

        # title = ""   # if not using reweighted samples
        title = ", Reweighted"  # if using reweighted samples

        do_pt_plot(PYTHIA_DIR, HERWIG_DIR, "QCD", "Dijet/pt_jet1", os.path.join(odir, "qcd_pt_jet1.pdf"), "Dijet"+title)
        do_pt_plot(PYTHIA_DIR, HERWIG_DIR, "QCD", "Dijet/pt_jet2", os.path.join(odir, "qcd_pt_jet2.pdf"), "Dijet"+title)
        do_pt_plot(PYTHIA_DIR, HERWIG_DIR, "QCD", "Dijet_QG/jet_pt", os.path.join(odir, "qcd_pt_jet_used.pdf"), "Dijet"+title)
        do_pt_plot(PYTHIA_DIR, HERWIG_DIR, "QCD", "Dijet_genjet/genjet_pt", os.path.join(odir, "qcd_pt_genjet_used.pdf"), "Dijet"+title)

        if "pythiaOnlyFlat" not in HERWIG_DIR:
            do_pt_plot(PYTHIA_DIR, HERWIG_DIR, "DyJetsToLL", "ZPlusJets/pt_jet1", os.path.join(odir, "zpj_pt_jet1.pdf"), "Z+jets"+title)
            do_pt_plot(PYTHIA_DIR, HERWIG_DIR, "DyJetsToLL", "ZPlusJets_QG/jet_pt", os.path.join(odir, "zpj_pt_jet_used.pdf"), "Z+jets"+title)
            do_pt_plot(PYTHIA_DIR, HERWIG_DIR, "DyJetsToLL", "ZPlusJets_genjet/genjet_pt", os.path.join(odir, "zpj_pt_genjet_used.pdf"), "Z+jets"+title)
