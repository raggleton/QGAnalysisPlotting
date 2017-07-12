#!/usr/bin/env python

"""Print basic plots comparing Pythia & herwig distributions"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
from comparator import Contribution, Plot, grab_obj
import common_utils as cu
import numpy as np
import os

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

PYTHIA_DIR = "workdir_ak4chs"
# HERWIG_DIR = "workdir_ak4chs_herwig"
HERWIG_DIR = "workdir_ak4chs_herwig_reweight"


def do_pt_plot(selection, hist_name, output_name, title=""):
    h_pythia = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_%s_.root" % (PYTHIA_DIR, selection), hist_name)
    h_herwig = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_%s_.root" % (HERWIG_DIR, selection), hist_name)

    c_pythia = Contribution(h_pythia, label="Pythia", line_color=ROOT.kBlue, marker_color=ROOT.kBlue, fill_color=ROOT.kBlue, normalise_hist=True)
    c_herwig = Contribution(h_herwig, label="Herwig", line_color=ROOT.kRed, normalise_hist=True)

    p = Plot([c_pythia, c_herwig], what="hist", legend=True, subplot_type='ratio', subplot=c_pythia, title=title)
    p.plot("NOSTACK HISTE")
    p.main_pad.SetLogy()
    p.container.SetMinimum(1E-10)
    p.subplot_container.SetMaximum(1.5)
    p.subplot_container.SetMinimum(0)
    p.canvas.Update()
    p.save(output_name)


if __name__ == "__main__":
    odir = "Pythia_Herwig_reweight"
    title = ", Reweighted"
    do_pt_plot("QCD", "Dijet/pt_jet1", os.path.join(odir, "pythia_herwig_qcd_pt_jet1.pdf"), "Dijet"+title)
    do_pt_plot("QCD", "Dijet/pt_jet2", os.path.join(odir, "pythia_herwig_qcd_pt_jet2.pdf"), "Dijet"+title)
    do_pt_plot("QCD", "Dijet_QG/jet_pt", os.path.join(odir, "pythia_herwig_qcd_pt_jet_used.pdf"), "Dijet"+title)
    do_pt_plot("QCD", "Dijet_genjet/genjet_pt", os.path.join(odir, "pythia_herwig_qcd_pt_genjet_used.pdf"), "Dijet"+title)

    do_pt_plot("DyJetsToLL", "ZPlusJets/pt_jet1", os.path.join(odir, "pythia_herwig_zpj_pt_jet1.pdf"), "Z+jets"+title)
    do_pt_plot("DyJetsToLL", "ZPlusJets_QG/jet_pt", os.path.join(odir, "pythia_herwig_zpj_pt_jet_used.pdf"), "Z+jets"+title)
    do_pt_plot("DyJetsToLL", "ZPlusJets_genjet/genjet_pt", os.path.join(odir, "pythia_herwig_zpj_pt_genjet_used.pdf"), "Z+jets"+title)
