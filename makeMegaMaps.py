#!/usr/bin/env python

"""Make mega response maps"""

import ROOT
# from MyStyle import My_Style
# My_Style.cd()
import os
from collections import OrderedDict
import sys
from array import array
from bisect import bisect_left
from copy import deepcopy
from math import sqrt
from itertools import product


# My stuff
# from comparator import Contribution, Plot, grab_obj
# import qg_common as qgc
# import qg_general_plots as qgg
import common_utils as cu

# For debugging
# import sys
# import tracers
# sys.settrace(tracers.trace_calls)
# sys.settrace(tracers.trace_calls_detail)
# sys.settrace(tracers.trace_calls_and_returns)

# import hunter
# hunter.trace(module='comparator')

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit()
# ROOT.gStyle.SetStatColor()

# Control plot output format
OUTPUT_FMT = "png"


def make_mega_maps(filename, var, out_filename):
    rootfile = cu.open_root_file(filename)
    outfile = cu.open_root_file(out_filename, "RECREATE")

    plot_dirname = "Dijet_QG_tighter"
    plot_dir = cu.get_from_tfile(rootfile, plot_dirname)
    var_stem = "jet_%s_response_bin_" % (var)
    all_plot_names = [x for x in cu.get_list_of_element_names(plot_dir) if var_stem in x]
    n_bins = int(sqrt(len(all_plot_names)))
    if len(all_plot_names) == 0:
        raise RuntimeError("No response plots in file")

    mult_hists = {histname : cu.get_from_tfile(rootfile, os.path.join(plot_dirname, histname)) for histname in all_plot_names}
    
    example_hist = mult_hists[all_plot_names[0]]
    n_bins_mult_x = example_hist.GetNbinsX()
    x_min = example_hist.GetXaxis().GetXmin()
    x_max = example_hist.GetXaxis().GetXmax()
    n_bins_mult_y = example_hist.GetNbinsY()
    y_min = example_hist.GetYaxis().GetXmin()
    y_max = example_hist.GetYaxis().GetXmax()
    h_mega = ROOT.TH2F("mega", ";"+example_hist.GetXaxis().GetTitle()+";"+example_hist.GetYaxis().GetTitle(), 
                       n_bins_mult_x*n_bins, x_min, x_max*n_bins,
                       n_bins_mult_y*n_bins, y_min, y_max*n_bins)
    print("mega hist has", h_mega.GetNbinsX(), "xbins and", h_mega.GetNbinsY(), "ybins")

    do_inds = list(range(0, n_bins))
    for xind in do_inds:
        for yind in do_inds:
            histname = "jet_%s_response_bin_%s_%s" % (var, xind, yind)
            print(histname)
            this_hist = mult_hists[histname]
            x_offset = (xind-do_inds[0])*n_bins_mult_x
            y_offset = (yind-do_inds[0])*n_bins_mult_y
            for x, y in product(range(1, n_bins_mult_x+1), range(1, n_bins_mult_y+1)):
                h_mega.SetBinContent(x+x_offset, y+y_offset, this_hist.GetBinContent(x, y))

    del mult_hists
    canv = ROOT.TCanvas("c", "", 5000, 5000)
    canv.SetLogz(1)
    # h_mega.Draw("COLZ")
    # canv.SaveAs(os.path.join(os.path.dirname(filename), "mega_map_%s.%s" % (var, OUTPUT_FMT)))
    # outfile.cd()
    # h_mega.Write()

    canv.Clear()
    h_renormx = cu.make_normalised_TH2(h_mega, "X", False)
    h_renormx.SetMinimum(1E-5)
    h_renormx.Draw("COLZ")
    # h_renormx.SetMinimum()
    canv.SaveAs(os.path.join(os.path.dirname(filename), "mega_map_%s_renormX.%s" % (var, OUTPUT_FMT)))
    outfile.cd()
    h_renormx.Write()
    del h_renormx

    # canv.Clear()
    # h_renormy = cu.make_normalised_TH2(h_mega, "Y", False)
    # h_renormy.Draw("COLZ")
    # # h_renormx.SetMinimum()
    # canv.SaveAs(os.path.join(os.path.dirname(filename), "mega_map_%s_renormY.%s" % (var, OUTPUT_FMT)))

    rootfile.Close()


if __name__ == "__main__":
    # make_mega_maps("workdir_ak4puppi_pythia_newFlav_binByAve_trigBinningBetter_jetAsymCut_multResponse/uhh2.AnalysisModuleRunner.MC.MC_PYTHIA-QCD.root", "multiplicity")
    # make_mega_maps("workdir_ak4puppi_pythia_newFlav_binByAve_trigBinningBetter_jetAsymCut_multResponse/uhh2.AnalysisModuleRunner.MC.MC_PYTHIA-QCD.root", "puppiMultiplicity")
    # make_mega_maps("workdir_ak4puppi_pythia_newFlav_binByAve_trigBinningBetter_jetAsymCut_multResponse/uhh2.AnalysisModuleRunner.MC.MC_PYTHIA-QCD.root", "lha")
    # make_mega_maps("workdir_ak4puppi_mgpythia_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1Constituents_V11JEC_JER/uhh2.AnalysisModuleRunner.MC.MC_QCD.root", "lha", "workdir_ak4puppi_mgpythia_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1Constituents_V11JEC_JER/mega_map.root")
    make_mega_maps("workdir_ak4puppi_mgpythia_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1Constituents_V11JEC_JER_murUp/uhh2.AnalysisModuleRunner.MC.MC_QCD.root", "lha", "workdir_ak4puppi_mgpythia_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1Constituents_V11JEC_JER_murUp/mega_map.root")
