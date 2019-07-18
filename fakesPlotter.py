#!/usr/bin/env python

"""Print plots relating to fakes in unfolding"""


import argparse
import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
import sys

# My stuff
from comparator import Contribution, Plot, grab_obj
import qg_common as qgc
# import qg_general_plots as qgg
import common_utils as cu

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit()
ROOT.gStyle.SetPaintTextFormat(".3f")
# ROOT.gStyle.SetStatColor()

# Control plot output format
OUTPUT_FMT = "pdf"


def print_cutflow(hist, output_filename, logy=False):
    """Print cutflwo hist to file, scaling bins such that
    they correspond to fractions of the 1st bin"""
    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetGrid()
    hist.Scale(1/hist.GetBinContent(1))
    hist.Draw("HIST TEXT")
    hist.GetYaxis().SetTitle("Fraction")
    if not logy:
        hist.SetMinimum(0)
        hist.SetMaximum(1.2)
    else:
        canv.SetLogy()
        hist.SetMinimum(1E-3)
        hist.SetMaximum(10)
        canv.Update()
    canv.SaveAs(output_filename)


def print_hist_comparison(entries, plots_kwargs, output_filename):
    """Print multiple hists on canvas, no rescaling, no stacking

    entries: list[(object, kwargs for Contribution)]
    """
    conts = [Contribution(e[0], **e[1]) for e in entries]
    logy = plots_kwargs.get('logy', False)
    if "logy" in plots_kwargs:
        del plots_kwargs['logy']
    plot = Plot(conts, what="hist", **plots_kwargs)
    plot.plot("HISTE NOSTACK")
    if logy:
        plot.set_logy()
    plot.save(output_filename)


def get_fraction_passing_cut(hist, cut_value_min=None, cut_value_max=None):
    """Get integral in region cut_value_min < x < cut_value_max"""
    min_bin = 0
    if cut_value_min:
        min_bin = hist.FindBin(cut_value_min)
    max_bin = hist.GetNbinsX()
    if cut_value_max:
        max_bin = hist.FindBin(cut_value_max) -1 # or no -1?
    number = hist.Integral(min_bin, max_bin)
    return number / hist.Integral()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Input file with plots")
    parser.add_argument("--outputDir", help="Directory to put output plot dirs into", default="fakes")
    args = parser.parse_args()

    cu.check_dir_exists_create(args.outputDir)

    in_tfile = cu.open_root_file(args.input)
    raw_filename = os.path.basename(args.input)
    if "QCD" in raw_filename:
        print_cutflow(cu.get_from_tfile(in_tfile, "cf_DijetGenSel"), os.path.join(args.outputDir, "cutflow_DijetGenSel.%s" % OUTPUT_FMT))
        print_cutflow(cu.get_from_tfile(in_tfile, "cf_DijetGenSelPassReco"), os.path.join(args.outputDir, "cutflow_DijetGenSelPassReco.%s" % OUTPUT_FMT))
        unfold_fakes_counter = cu.get_from_tfile(in_tfile, "Dijet_QG_Unfold_tighter/fakes_counter")
        print_cutflow(unfold_fakes_counter, os.path.join(args.outputDir, "cutflow_DijetUnfoldFakes.%s" % OUTPUT_FMT), logy=True)
        unfold_fakes_counter_scaled = unfold_fakes_counter.Clone()
        unfold_fakes_counter_scaled.Scale(1./unfold_fakes_counter_scaled.GetBinContent(1))
        # printout values per bin
        print("Fraction fakes in unfolding:")
        for i in range(2, unfold_fakes_counter_scaled.GetNbinsX()+1):
            print(unfold_fakes_counter_scaled.GetXaxis().GetBinLabel(i), ":", unfold_fakes_counter_scaled.GetBinContent(i))
        # can't use normal get_fraction_passing_cut as scales weirdly
        print("Fraction of fakes in unfolding:", print(unfold_fakes_counter_scaled.Integral(unfold_fakes_counter_scaled.FindBin(0), 4)))  # 0 as weird bin numbering, 1st bin is -1
        lw = 2
        print_hist_comparison(
            entries=[
                # (cu.get_from_tfile(in_tfile, "GenJetsPresel/eta_1"), {"line_color": ROOT.kBlack, "marker_color": ROOT.kBlack, "line_width": lw, "label": "passGen || passReco"}),
                (cu.get_from_tfile(in_tfile, "GenJetsPassDijetReco/eta_1"), {"line_color": ROOT.kRed, "marker_color": ROOT.kRed, "line_width": lw, "label": "pass dijet reco cuts"})
            ],
            plots_kwargs={"xtitle": "#eta_{GenJet 1}", "ytitle": "N", "has_data": False},
            output_filename=os.path.join(args.outputDir, "dijet_genjet_eta1.%s" % OUTPUT_FMT)
        )
        print_hist_comparison(
            entries=[
                # (cu.get_from_tfile(in_tfile, "GenJetsPresel/eta_2"), {"line_color": ROOT.kBlack, "marker_color": ROOT.kBlack, "line_width": lw, "label": "passGen || passReco"}),
                (cu.get_from_tfile(in_tfile, "GenJetsPassDijetReco/eta_2"), {"line_color": ROOT.kRed, "marker_color": ROOT.kRed, "line_width": lw, "label": "pass dijet reco cuts"})
            ],
            plots_kwargs={"xtitle": "#eta_{GenJet 2}", "ytitle": "N", "has_data": False},
            output_filename=os.path.join(args.outputDir, "dijet_genjet_eta2.%s" % OUTPUT_FMT)
        )
        eta_cut = 1.8
        print("Fraction failing eta1 cut:", 1-get_fraction_passing_cut(cu.get_from_tfile(in_tfile, "GenJetsPassDijetReco/eta_1"), -eta_cut, eta_cut))
        print("Fraction failing eta2 cut:", 1-get_fraction_passing_cut(cu.get_from_tfile(in_tfile, "GenJetsPassDijetReco/eta_2"), -eta_cut, eta_cut))
        print_hist_comparison(
            entries=[
                # (cu.get_from_tfile(in_tfile, "GenJetsPresel/eta_3"), {"line_color": ROOT.kBlack, "marker_color": ROOT.kBlack, "line_width": lw, "label": "passGen || passReco"}),
                (cu.get_from_tfile(in_tfile, "GenJetsPassDijetReco/eta_3"), {"line_color": ROOT.kRed, "marker_color": ROOT.kRed, "line_width": lw, "label": "pass dijet reco cuts"})
            ],
            plots_kwargs={"xtitle": "#eta_{GenJet 3}", "ytitle": "N", "has_data": False},
            output_filename=os.path.join(args.outputDir, "dijet_genjet_eta3.%s" % OUTPUT_FMT)
        )
        print_hist_comparison(
            entries=[
                (cu.get_from_tfile(in_tfile, "GenJetsPresel/pt_1"), {"line_color": ROOT.kBlack, "marker_color": ROOT.kBlack, "line_width": lw, "label": "passGen || passReco"}),
                (cu.get_from_tfile(in_tfile, "GenJetsPassDijetReco/pt_1"), {"line_color": ROOT.kRed, "marker_color": ROOT.kRed, "line_width": lw, "label": "passReco"})
            ],
            plots_kwargs={"xtitle": "p_{T}^{GenJet 1} [GeV]", "ytitle": "N", "has_data": False, "xlim": [0, 200], "ylim": [1E7, 5E12], "logy": True},
            output_filename=os.path.join(args.outputDir, "dijet_genjet_pt1.%s" % OUTPUT_FMT)
        )
        print_hist_comparison(
            entries=[
                (cu.get_from_tfile(in_tfile, "GenJetsPresel/pt_2"), {"line_color": ROOT.kBlack, "marker_color": ROOT.kBlack, "line_width": lw, "label": "passGen || passReco"}),
                (cu.get_from_tfile(in_tfile, "GenJetsPassDijetReco/pt_2"), {"line_color": ROOT.kRed, "marker_color": ROOT.kRed, "line_width": lw, "label": "passReco"})
            ],
            plots_kwargs={"xtitle": "p_{T}^{GenJet 2} [GeV]", "ytitle": "N", "has_data": False, "xlim": [0, 200], "ylim": [1E7, 5E12], "logy": True},
            output_filename=os.path.join(args.outputDir, "dijet_genjet_pt2.%s" % OUTPUT_FMT)
        )
        pt_cut = 10
        print("Fraction failing pt1 cut:", 1-get_fraction_passing_cut(cu.get_from_tfile(in_tfile, "GenJetsPassDijetReco/pt_1"), pt_cut, None))
        print("Fraction failing pt2 cut:", 1-get_fraction_passing_cut(cu.get_from_tfile(in_tfile, "GenJetsPassDijetReco/pt_2"), pt_cut, None))
        print_hist_comparison(
            entries=[
                (cu.get_from_tfile(in_tfile, "GenJetsPresel/pt_3"), {"line_color": ROOT.kBlack, "marker_color": ROOT.kBlack, "line_width": lw, "label": "passGen || passReco"}),
                (cu.get_from_tfile(in_tfile, "GenJetsPassDijetReco/pt_3"), {"line_color": ROOT.kRed, "marker_color": ROOT.kRed, "line_width": lw, "label": "passReco"})
            ],
            plots_kwargs={"xtitle": "p_{T}^{GenJet 3} [GeV]", "ytitle": "N", "has_data": False, "xlim": [0, 200], "ylim": [1E7, 5E12], "logy": True},
            output_filename=os.path.join(args.outputDir, "dijet_genjet_pt3.%s" % OUTPUT_FMT)
        )
    elif "DYJetsToLL" in raw_filename:
        print_cutflow(cu.get_from_tfile(in_tfile, "cf_ZPlusJetsGenSel"), os.path.join(args.outputDir, "cutflow_ZPlusJetsGenSel.%s" % OUTPUT_FMT))
        print_cutflow(cu.get_from_tfile(in_tfile, "cf_ZPlusJetsGenSelPassReco"), os.path.join(args.outputDir, "cutflow_ZPlusJetsGenSelPassReco.%s" % OUTPUT_FMT))
        print_cutflow(cu.get_from_tfile(in_tfile, "ZPlusJets_QG_Unfold/fakes_counter"), os.path.join(args.outputDir, "cutflow_ZPlusJetsUnfoldFakes.%s" % OUTPUT_FMT))

