#!/usr/bin/env python


"""
Print PF matching plots
"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse
from math import sqrt

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp


# ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat(".3f")
ROOT.gStyle.SetHistTopMargin(0.)


PF_DICT = {
    # 0: "Unknown",
    1: "Charged hadron",
    2: "Electron",
    3: "Muon",
    4: "Photon",
    5: "Neutral hadron",
}


def do_3d_plot(hist, pf_name, output_filename):
    canvas = ROOT.TCanvas(ROOT.TUUID().AsString(), "", 800, 600)
    canvas.SetTicks(1, 1)
    canvas.SetLeftMargin(0.13)
    canvas.SetBottomMargin(0.11)
    hist.SetTitle("Gen-PF matching for %s PF particles;Charged hadron GenParticle p_{T} [GeV];PF particle p_{T} [GeV];#DeltaR(Gen, PF)" % (pf_name))
    hist.GetXaxis().SetTitleOffset(1.75)
    hist.GetYaxis().SetTitleOffset(1.75)
    hist.GetZaxis().SetTitleOffset(1.5)
    hist.Draw("scat=0.1")# this does nothing?
    ROOT.gPad.SetTickx(1)
    ROOT.gPad.SetTicky(1)
    canvas.SaveAs(output_filename)


def do_pt_response_map(hist, pf_name, output_filename):
    """Plot heatmap of PF vs genparticle pT"""
    h2d = hist.Project3D('yx')
    h2d.SetName(cu.get_unique_str())
    qgp.do_2D_plot(h2d, output_filename, title="Gen-PF matching for %s PF particles;Charged hadron GenParticle p_{T} [GeV];PF particle p_{T} [GeV]" % (pf_name))
    # do a log Z version
    stem, ext = os.path.splitext(output_filename)
    log_output_filename = stem + "_logZ" + ext
    qgp.do_2D_plot(h2d, log_output_filename, title="Gen-PF matching for %s PF particles;Charged hadron GenParticle p_{T} [GeV];PF particle p_{T} [GeV]" % (pf_name), logz=True)


def do_deltaR_vs_gen_pT_map(hist, pf_name, output_filename):
    """Plot heatmap of deltaR vs genparticle pT"""
    h2d = hist.Project3D('zx')
    h2d.SetName(cu.get_unique_str())
    qgp.do_2D_plot(h2d, output_filename, title="Gen-PF matching for %s PF particles;Charged hadron GenParticle p_{T} [GeV];#DeltaR(Gen, PF)" % (pf_name))
    # do a log Z version
    stem, ext = os.path.splitext(output_filename)
    log_output_filename = stem + "_logZ" + ext
    qgp.do_2D_plot(h2d, log_output_filename, title="Gen-PF matching for %s PF particles;Charged hadron GenParticle p_{T} [GeV];#DeltaR(Gen, PF)" % (pf_name), logz=True)


def do_pf_fraction_plot(hist_map, pt_bins, output_filename):
    """Plot PF particle type fractioin for matches, binned by GenParticle pT"""
    entries = []
    for pt_low, pt_high, mark in zip(pt_bins[:-1], pt_bins[1:], cu.Marker().cycle()):
        values = {}
        for pf_ind, (pf_name, hist) in hist_map.items():
            ax = hist.GetXaxis()
            binx1 = ax.FindBin(pt_low)
            binx2 = ax.FindBin(pt_high)-1
            if pt_high == ax.GetBinUpEdge(ax.GetLast()):
                binx2 = ax.GetLast()
            biny1 = 1
            biny2 = hist.GetNbinsY()
            binz1 = 1
            binz2 = hist.GetNbinsZ()
            values[pf_ind] = hist.Integral(binx1, binx2, biny1, biny2, binz1, binz2)  # integral includes the last bin

        sum_values = sum(values.values())
        fracs = {k:(v/sum_values) for k, v in values.items()}

        h = ROOT.TH1D("h_pt_bin_%gto%g" % (pt_low, pt_high), "", len(values), 0, len(values))
        ax = h.GetXaxis()
        for ind, k in enumerate(sorted(fracs.keys()), 1):
            h.SetBinContent(ind, fracs[k])
            h.SetBinError(ind, sqrt(values[k]) / sum_values)
            ax.SetBinLabel(ind, hist_map[k][0]);

        c = Contribution(h,
                         label='%g < GenParticle p_{T} < %g GeV' % (pt_low, pt_high),
                         line_width=1,
                         marker_size=0.75, marker_style=mark,
                         normalise_hist=False)
        entries.append(c)

    ROOT.gStyle.SetPalette(55)
    plot = Plot(entries, 'hist',
                xtitle='PF particle type',
                ytitle='Fraction matched as type',
                ylim=(1E-3, 2),
                has_data=False)
    plot.default_canvas_size = (800, 600)
    plot.plot("NOSTACK PMC PLC HISTE")
    plot.set_logy(do_more_labels=False)
    plot.save(output_filename)
    ROOT.gStyle.SetPalette(ROOT.kViridis)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",
                        help="Source file")
    parser.add_argument("--outputDir",
                        help="Output dir (default is source dir)")
    args = parser.parse_args()

    if not os.path.isfile(args.source):
        raise IOError("source file doesn't exist")

    in_tfile = cu.TFileCacher(args.source)
    is_data = '.DATA.' in os.path.basename(args.source)
    if is_data:
        print("This script is only useful for MC files, exiting")
        exit()

    is_dijet = '_QCD.root' in os.path.basename(args.source)

    if not args.outputDir:
        in_dir = os.path.dirname(args.source)
        args.outputDir = os.path.join(in_dir, 'pf_matching_plots')
        if is_dijet:
            args.outputDir += "_dijet"
        else:
            args.outputDir += "_zpj"

    if not os.path.isdir(args.outputDir):
        os.makedirs(args.outputDir)

    print("Plots produced in", args.outputDir)

    jet_str = "AK4 PUPPI jets"
    if 'ak8puppi' in args.source:
        jet_str = "AK8 PUPPI jets"

    # DIJET
    if is_dijet:
        for pf_ind, pf_name in PF_DICT.items():
            print(pf_name)
            # 3D plot
            do_3d_plot(in_tfile.Get("MCTrackScaleFactor_PF_%d" % pf_ind),
                       pf_name=pf_name,
                       output_filename=os.path.join(args.outputDir, 'genpt_vs_recopt_vs_deltaR_%d_%s.pdf' % (pf_ind, pf_name.replace(" ", "_"))))

            # Print deltaR

            # Print deltaR vs gen pT
            do_deltaR_vs_gen_pT_map(in_tfile.Get("MCTrackScaleFactor_PF_%d" % pf_ind),
                                    pf_name=pf_name,
                                    output_filename=os.path.join(args.outputDir, 'deltaR_vs_genpT_map_%d_%s.pdf' % (pf_ind, pf_name.replace(" ", "_"))))

            # Print response 2D map
            do_pt_response_map(in_tfile.Get("MCTrackScaleFactor_PF_%d" % pf_ind),
                            pf_name=pf_name,
                            output_filename=os.path.join(args.outputDir, 'response_map_%d_%s.pdf' % (pf_ind, pf_name.replace(" ", "_"))))

        # print PF fractions as a function of gen particle pt
        hist_map = {pf_ind: [pf_name, in_tfile.Get("MCTrackScaleFactor_PF_%d" % pf_ind)] for pf_ind, pf_name in PF_DICT.items()}
        pt_bins = [0, 1, 2, 3, 4, 5, 10, 15, 20]
        do_pf_fraction_plot(hist_map,
                            pt_bins=pt_bins,
                            output_filename=os.path.join(args.outputDir, 'pf_fraction.pdf'))

    # Z+Jet
    else:
        pass