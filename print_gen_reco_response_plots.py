#!/usr/bin/env python

"""Print plots of gen vs reco response plots"""


import argparse
# from MyStyle import My_Style
# My_Style.cd()
import os
from itertools import product
from array import array

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)


# My stuff
from comparator import Contribution, Plot
# import qg_common as qgc
# import qg_general_plots as qgg
import common_utils as cu

ROOT.gStyle.SetPaintTextFormat(".3f")

# Control output format
OUTPUT_FMT = "pdf"


def do_pt_response_plot(tdir, plot_dir):
    """Do 2D plots of genjet pt vs recojet pt"""
    h2d = cu.get_from_tfile(tdir, "pt_jet_response")
    renorm_h2d = cu.make_normalised_TH2(h2d, 'X', recolour=False)
    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    pad = ROOT.gPad
    pad.SetBottomMargin(0.12)
    pad.SetLeftMargin(0.13)
    pad.SetRightMargin(0.12)
    canv.SetTicks(1, 1)
    canv.SetLogx()
    canv.SetLogy()
    renorm_h2d.SetMarkerSize(0.5)
    renorm_h2d.Draw("COLZ TEXT")
    
    xax = renorm_h2d.GetXaxis()
    upper_lim = renorm_h2d.GetXaxis().GetBinUpEdge(renorm_h2d.GetXaxis().GetLast())
    # upper_lim = 5000
    title_offset = 1.5
    renorm_h2d.SetTitleOffset(title_offset, 'X')
    xax.SetRangeUser(10, upper_lim)
    xax.SetMoreLogLabels()

    yax = renorm_h2d.GetYaxis()
    renorm_h2d.SetTitleOffset(title_offset*1.15, 'Y')
    yax.SetRangeUser(10, upper_lim)
    yax.SetMoreLogLabels()
    canv.Update()

    plot_dir = os.path.join(plot_dir, tdir.GetName())
    cu.check_dir_exists_create(plot_dir)
    canv.SaveAs(os.path.join(plot_dir, "pt_jet_response_renormX_linZ.%s" % (OUTPUT_FMT)))
    canv.SetLogz()
    renorm_h2d.SetMinimum(1E-3)
    canv.SaveAs(os.path.join(plot_dir, "pt_jet_response_renormX_logZ.%s" % (OUTPUT_FMT)))

    # Now do plot of purity, etc
    stability_values, stability_errs = [], []
    xfer_down_values, xfer_down_errs = [], []
    xfer_down2_values, xfer_down2_errs = [], []
    xfer_up_values, xfer_up_errs = [], []
    xfer_up2_values, xfer_up2_errs = [], []

    binning = [renorm_h2d.GetXaxis().GetBinLowEdge(bin_ind) for bin_ind in range(1, renorm_h2d.GetNbinsX()+2)]

    hist_stability = ROOT.TH1F("stability"+cu.get_unique_str(), ";p_{T}^{Gen} [GeV];Fraction", len(binning)-1, array('d', binning)) 
    hist_xfer_down = ROOT.TH1F("xfer_down"+cu.get_unique_str(), ";p_{T}^{Gen} [GeV];Fraction", len(binning)-1, array('d', binning)) 
    hist_xfer_down2 = ROOT.TH1F("xfer_down2"+cu.get_unique_str(), ";p_{T}^{Gen} [GeV];Fraction", len(binning)-1, array('d', binning)) 
    hist_xfer_down3 = ROOT.TH1F("xfer_down3"+cu.get_unique_str(), ";p_{T}^{Gen} [GeV];Fraction", len(binning)-1, array('d', binning)) 
    hist_xfer_up = ROOT.TH1F("xfer_up"+cu.get_unique_str(), ";p_{T}^{Gen} [GeV];Fraction", len(binning)-1, array('d', binning)) 
    hist_xfer_up2 = ROOT.TH1F("xfer_up2"+cu.get_unique_str(), ";p_{T}^{Gen} [GeV];Fraction", len(binning)-1, array('d', binning)) 
    hist_xfer_up3 = ROOT.TH1F("xfer_up3"+cu.get_unique_str(), ";p_{T}^{Gen} [GeV];Fraction", len(binning)-1, array('d', binning)) 

    for bin_ind in range(1, renorm_h2d.GetNbinsX()+1):
        hist_stability.SetBinContent(bin_ind, renorm_h2d.GetBinContent(bin_ind, bin_ind))
        hist_stability.SetBinError(bin_ind, renorm_h2d.GetBinError(bin_ind, bin_ind))
        
        hist_xfer_down.SetBinContent(bin_ind, renorm_h2d.GetBinContent(bin_ind, bin_ind-1))
        hist_xfer_down.SetBinError(bin_ind, renorm_h2d.GetBinError(bin_ind, bin_ind-1))
        hist_xfer_down2.SetBinContent(bin_ind, renorm_h2d.GetBinContent(bin_ind, bin_ind-2))
        hist_xfer_down2.SetBinError(bin_ind, renorm_h2d.GetBinError(bin_ind, bin_ind-2))
        hist_xfer_down3.SetBinContent(bin_ind, renorm_h2d.GetBinContent(bin_ind, bin_ind-3))
        hist_xfer_down3.SetBinError(bin_ind, renorm_h2d.GetBinError(bin_ind, bin_ind-3))
        
        hist_xfer_up.SetBinContent(bin_ind, renorm_h2d.GetBinContent(bin_ind, bin_ind+1))
        hist_xfer_up.SetBinError(bin_ind, renorm_h2d.GetBinError(bin_ind, bin_ind+1))
        hist_xfer_up2.SetBinContent(bin_ind, renorm_h2d.GetBinContent(bin_ind, bin_ind+2))
        hist_xfer_up2.SetBinError(bin_ind, renorm_h2d.GetBinError(bin_ind, bin_ind+2))
        hist_xfer_up3.SetBinContent(bin_ind, renorm_h2d.GetBinContent(bin_ind, bin_ind+3))
        hist_xfer_up3.SetBinError(bin_ind, renorm_h2d.GetBinError(bin_ind, bin_ind+3))
   
    col_stability = ROOT.kBlack
    col_xfer_down = ROOT.kRed
    col_xfer_down2 = ROOT.kMagenta
    col_xfer_down3 = ROOT.kGreen
    col_xfer_up = ROOT.kBlue-4
    col_xfer_up2 = ROOT.kAzure+8
    col_xfer_up3 = ROOT.kOrange
    contributions = [
        Contribution(hist_stability, label="Correct bin", line_color=col_stability, marker_color=col_stability),
        Contribution(hist_xfer_down, label="1 lower bin", line_color=col_xfer_down, marker_color=col_xfer_down),
        Contribution(hist_xfer_down2, label="2 lower bin", line_color=col_xfer_down2, marker_color=col_xfer_down2),
        # Contribution(hist_xfer_down3, label="3 lower bin", line_color=col_xfer_down3, marker_color=col_xfer_down3),
        Contribution(hist_xfer_up, label="1 higher bin", line_color=col_xfer_up, marker_color=col_xfer_up),
        Contribution(hist_xfer_up2, label="2 higher bin", line_color=col_xfer_up2, marker_color=col_xfer_up2),
        # Contribution(hist_xfer_up3, label="3 higher bin", line_color=col_xfer_up3, marker_color=col_xfer_up3),
    ]
    xlim = [10, binning[-1]]
    plot = Plot(contributions, what='hist', xlim=xlim, ylim=[1e-4, 2])
    y1 = 0.15
    plot.legend.SetY1(y1)
    plot.legend.SetY2(y1+0.2)
    plot.plot("NOSTACK HISTE")
    plot.set_logx()
    plot.set_logy()
    plot.main_pad.cd()
    lines = []
    for val in [1, 0.5, 1e-1, 1e-2, 1e-3]:
        line = ROOT.TLine(xlim[0], val, xlim[1], val)
        line.SetLineStyle(2)
        lines.append(line)
        line.Draw("same")
    plot.save(os.path.join(plot_dir, 'migration_summary.%s' % (OUTPUT_FMT)))


def do_jet_index_plots(tdir, plot_dir):
    """Do 2D plots of genjet index vs recojet index"""
    plot_dir = os.path.join(plot_dir, tdir.GetName())
    cu.check_dir_exists_create(plot_dir)
    stem = "genjet_ind_recojet_ind_pt_"
    for plot_name in cu.get_list_of_element_names(tdir):
        if not plot_name.startswith(stem):
            continue
        h2d = cu.get_from_tfile(tdir, plot_name)
        # renorm_h2d = cu.make_normalised_TH2(h2d, 'X')
        renorm_h2d = h2d
        canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
        pad = ROOT.gPad
        pad.SetBottomMargin(0.12)
        pad.SetLeftMargin(0.13)
        pad.SetRightMargin(0.12)
        canv.SetTicks(1, 1)
        renorm_h2d.Draw("COLZ")

        title_offset = 1.5
        renorm_h2d.SetTitleOffset(title_offset, 'X')
        renorm_h2d.SetTitleOffset(title_offset*1.15, 'Y')

        canv.SaveAs(os.path.join(plot_dir, "%s_renormX_linZ.%s" % (plot_name, OUTPUT_FMT)))
        # canv.SetLogz()
        # renorm_h2d.SetMinimum(1E-3)
        # canv.SaveAs(os.path.join(plot_dir, "%s_renormX_logZ.%s" % (plot_name, OUTPUT_FMT)))



def do_response_plots(in_file, plot_dir, do_these_dirs=None):
    tfile = cu.open_root_file(in_file)
    dirs = cu.get_list_of_element_names(tfile)

    for mydir in dirs:
        if do_these_dirs and mydir not in do_these_dirs:
            continue
        
        do_pt_response_plot(tfile.Get(mydir), plot_dir=plot_dir)

        do_jet_index_plots(tfile.Get(mydir), plot_dir=plot_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input', 
                        nargs='+',
                        help='Input ROOT files to process. '
                        'Several dirs can be specified here, separated by a space.')
    parser.add_argument("-o", "--output", help="Directory to put output plot dirs into", default=None)
    args = parser.parse_args()
    print(args)

    for in_file in args.input:
        default_plot_dir = os.path.join(os.path.dirname(in_file), "response_"+os.path.splitext(os.path.basename(in_file))[0])
        plot_dir = args.output if args.output else default_plot_dir
        cu.check_dir_exists_create(plot_dir)
        do_these_dirs = [
            'Dijet',
            'Dijet_tighter',
            'ZPlusJets',
        ]
        if "_QCD_" in in_file:
            do_these_dirs = do_these_dirs[:2]
        if "_DYJetsToLL_" in in_file:
            do_these_dirs = do_these_dirs[2:]
        do_response_plots(in_file, plot_dir=plot_dir, do_these_dirs=do_these_dirs)
