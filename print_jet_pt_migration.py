#!/usr/bin/env python

"""Print jet pt response matrix, and migration stats"""

from __future__ import print_function, division

import argparse
from MyStyle import My_Style
My_Style.cd()
import os
os.nice(10)


import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
# ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1111)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# My stuff
import qg_general_plots as qgp
import qg_common as qgc
import common_utils as cu

ROOT.gStyle.SetPaintTextFormat(".3f")

# Control output format
OUTPUT_FMT = "pdf"


def plot_jet_pt_response_matrix(h2d, h2d_renorm_x, h2d_renorm_y, title, plot_dir):
    """Plot 2D response matrix, absolute, and normed by x & y axes"""
    canv = ROOT.TCanvas("c"+cu.get_unique_str(), "", 700, 600)
    canv.SetTicks(1, 1)
    canv.SetLogx()
    canv.SetLogy()
    pad = ROOT.gPad
    pad.SetBottomMargin(0.12)
    pad.SetLeftMargin(0.13)
    pad.SetRightMargin(0.19)

    def _setup_map(hist_2d):
        xax = hist_2d.GetXaxis()
        xax.SetMoreLogLabels()
        xax.SetNoExponent()
        yax = hist_2d.GetYaxis()
        yax.SetMoreLogLabels()
        yax.SetNoExponent()
        xtitle_offset = 1.4
        ytitle_offset = xtitle_offset * 1.2
        ztitle_offset = xtitle_offset
        hist_2d.SetTitleOffset(xtitle_offset, 'X')
        hist_2d.SetTitleOffset(ytitle_offset, 'Y')
        hist_2d.SetTitleOffset(ztitle_offset, 'Z')

    _setup_map(h2d)
    h2d.SetMinimum(1E-3)
    var_label = 'p_{T}^{jet} [GeV]'
    title = "%s;%s (GEN);%s (RECO);N" % (title, var_label, var_label)
    h2d.SetTitle(title)
    h2d.Draw("COLZ")
    old_font_size = ROOT.gStyle.GetTitleFontSize()
    output_filename = os.path.join(plot_dir, "jet_pt.%s" % (OUTPUT_FMT))
    if not os.path.isdir(plot_dir):
        os.makedirs(plot_dir)
    canv.SaveAs(output_filename)

    canv.SetLogz()
    output_filename = os.path.join(plot_dir, "jet_pt_logZ.%s" % (OUTPUT_FMT))
    canv.SaveAs(output_filename)

    # Plot 2D map, renormalized by row
    canv.SetLogz(0)
    marker_size = 0.8
    _setup_map(h2d_renorm_y)
    h2d_renorm_y.SetTitle(title)
    h2d_renorm_y.SetMarkerSize(marker_size)
    h2d_renorm_y.SetMaximum(1)
    h2d_renorm_y.SetZTitle("P(GEN bin | RECO bin)")
    draw_opt = "COLZ"
    plot_migrations = True
    if plot_migrations:
        draw_opt += " TEXT45"
    h2d_renorm_y.Draw(draw_opt)
    canv.SaveAs(os.path.join(plot_dir, "jet_pt_renormY_linZ.%s" % (OUTPUT_FMT)))

    canv.SetLogz()
    h2d_renorm_y.SetMinimum(1E-3)
    canv.SaveAs(os.path.join(plot_dir, "jet_pt_renormY_logZ.%s" % (OUTPUT_FMT)))

    # Plot 2D map, renormalized by column
    canv.Clear()
    canv.SetLogz(0)
    _setup_map(h2d_renorm_x)
    h2d_renorm_x.SetTitle(title)
    h2d_renorm_x.SetMarkerSize(marker_size)
    h2d_renorm_x.SetMaximum(1)
    h2d_renorm_x.SetZTitle("P(RECO bin | GEN bin)")
    h2d_renorm_x.Draw(draw_opt)
    canv.SaveAs(os.path.join(plot_dir, "jet_pt_renormX_linZ.%s" % (OUTPUT_FMT)))

    canv.SetLogz()
    h2d_renorm_x.SetMinimum(1E-3)
    canv.SaveAs(os.path.join(plot_dir, "jet_pt_renormX_logZ.%s" % (OUTPUT_FMT)))

    ROOT.gStyle.SetTitleFontSize(old_font_size)


def do_jet_pt_migration_plot(input_filename, directory, title, output_dir):
    """Do migration stats plot"""
    tfile = cu.open_root_file(input_filename)
    h2d = tfile.Get("%s/jet_pt_vs_genjet_pt" % directory)

    h2d_new = h2d
    h2d_renorm_y = cu.make_normalised_TH2(h2d_new, 'Y', recolour=False, do_errors=True)
    h2d_renorm_x = cu.make_normalised_TH2(h2d_new, 'X', recolour=False, do_errors=True)

    # Plot 2D response matrix
    plot_jet_pt_response_matrix(h2d_new, h2d_renorm_x, h2d_renorm_y, title, output_dir)

    # Do migration metrics
    xlabel = 'p_{T}^{jet} [GeV]'
    qgp.make_migration_summary_plot(h2d_renorm_x,
                                    h2d_renorm_y,
                                    xlabel,
                                    title=title,
                                    log_var=True,
                                    output_filename=os.path.join(output_dir, "jet_pt_migration_summary.pdf"),
                                    do_reco_updown2=False,
                                    do_gen_updown=False)

    tfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('workdir',
                        help='Directory to process.')
    parser.add_argument("--outputDir",
                        help="Directory to put output plot dirs into",
                        default=None)
    args = parser.parse_args()

    if not args.outputDir:
        args.outputDir = os.path.join(args.workdir, "jet_pt_migration")
        print("Setting --outputDir to", args.outputDir)

    files_to_process = [
        (qgc.QCD_FILENAME, "QCD_MGPythia"),
        (qgc.DY_FILENAME, "DY_MGPythia"),
        (qgc.QCD_HERWIG_FILENAME, "QCD_Herwig"),
        (qgc.DY_HERWIG_FILENAME, "DY_Herwig"),
    ]

    for filename, file_label in files_to_process:
        input_file = os.path.join(args.workdir, filename)
        if not os.path.isfile(input_file):
            print("Cannot find input file %s - skipping" % args.input)
            continue

        sample_title = file_label.replace("_", " ")

        if "qcd" in file_label.lower():
            dirname = "Dijet_QG_central_tighter"
            do_jet_pt_migration_plot(input_file,
                                     dirname,
                                     title="%s (%s)" % (qgc.Dijet_CEN_LABEL, sample_title),
                                     output_dir=os.path.join(args.outputDir, '_'.join([dirname, file_label])))

            dirname = "Dijet_QG_forward_tighter"
            do_jet_pt_migration_plot(input_file,
                                     dirname,
                                     title="%s (%s)" % (qgc.Dijet_FWD_LABEL, sample_title),
                                     output_dir=os.path.join(args.outputDir, '_'.join([dirname, file_label])))
        else:
            dirname = "ZPlusJets_QG"
            do_jet_pt_migration_plot(input_file,
                                     dirname,
                                     title="%s (%s)" % (qgc.ZpJ_LABEL, sample_title),
                                     output_dir=os.path.join(args.outputDir, '_'.join([dirname, file_label])))
