#!/usr/bin/env python

"""Print plots comparing different matched flavours in different samples with different selections"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
os.nice(10)

# My stuff
from comparator import grab_obj
import qg_common as qgc
import qg_general_plots as qgg

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)

# Control plot output format
OUTPUT_FMT = "pdf"


def do_wrong_plots(root_dir, var_prepend="", plot_dir="wrong_flavs", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG", pt_bins=None, title=""):
    """Plot all the sample/selection/flavour combinations to check distributions indep of sample"""
    pt_bins = pt_bins or qgc.THEORY_PT_BINS
    plot_vars = ['jet_LHA', 'jet_pTD', 'jet_width', 'jet_thrust', 'jet_multiplicity']
    if "puppi" in root_dir.lower():
        plot_vars.append('jet_puppiMultiplicity')
    for v in plot_vars:
        v = "%s%s_vs_pt" % (var_prepend, v)

        h2d_dyj_chs = grab_obj(os.path.join(root_dir, qgc.DY_FILENAME), "%s/q%s" % (zpj_dirname, v))
        h2d_dyj_wrong_chs = grab_obj(os.path.join(root_dir, qgc.DY_FILENAME), "%s/g%s" % (zpj_dirname, v))
        h2d_dyj_qcd_chs = grab_obj(os.path.join(root_dir, qgc.DY_FILENAME), "%s/q%s" % (dj_dirname, v))
        h2d_dyj_qcd_wrong_chs = grab_obj(os.path.join(root_dir, qgc.DY_FILENAME), "%s/g%s" % (dj_dirname, v))
        h2d_qcd_chs = grab_obj(os.path.join(root_dir, qgc.QCD_FILENAME), "%s/g%s" % (dj_dirname, v))
        h2d_qcd_wrong_chs = grab_obj(os.path.join(root_dir, qgc.QCD_FILENAME), "%s/q%s" % (dj_dirname, v))

        lw = 1
        dy_kwargs_chs = dict(line_color=qgc.DY_COLOUR, fill_color=qgc.DY_COLOUR, label=qgc.DY_ZpJ_QFLAV_LABEL, line_width=lw, marker_style=20, marker_color=qgc.DY_COLOUR)
        dy_kwargs_wrong_chs = dict(line_color=qgc.DY_COLOUR+4, fill_color=qgc.DY_COLOUR+4, label=qgc.DY_ZpJ_GFLAV_LABEL, line_width=lw, line_style=1, marker_style=21, marker_color=qgc.DY_COLOUR+4)
        dy_kwargs_qcd_chs = dict(line_color=ROOT.kGreen+2, fill_color=ROOT.kGreen+2, label=qgc.DY_Dijet_QFLAV_LABEL, line_width=lw, line_style=1, marker_style=22, marker_color=ROOT.kGreen+2)
        dy_kwargs_qcd_wrong_chs = dict(line_color=ROOT.kOrange-1, fill_color=ROOT.kOrange-1, label=qgc.DY_Dijet_GFLAV_LABEL, line_width=lw, line_style=1, marker_style=23, marker_color=ROOT.kOrange-1)
        qcd_kwargs_chs = dict(line_color=qgc.QCD_COLOUR, fill_color=qgc.QCD_COLOUR, label=qgc.QCD_Dijet_GFLAV_LABEL, line_width=lw, marker_style=25, marker_color=qgc.QCD_COLOUR)
        qcd_kwargs_wrong_chs = dict(line_color=ROOT.kRed, fill_color=ROOT.kRed, label=qgc.QCD_Dijet_QFLAV_LABEL, line_width=lw, line_style=1, marker_style=24, marker_color=ROOT.kRed)

        rebin = 2
        xlim = None
        if "thrust" in v:
            rebin = 1
            xlim = (0, 0.5)
        if "multiplicity" in v.lower():
            rebin = 2
            xlim = (0, 100)

        for (start_val, end_val) in pt_bins:
            entries = [
                (qgg.get_projection_plot(h2d_dyj_chs, start_val, end_val), dy_kwargs_chs),
                (qgg.get_projection_plot(h2d_dyj_wrong_chs, start_val, end_val), dy_kwargs_wrong_chs),
                # (qgg.get_projection_plot(h2d_dyj_qcd_chs, start_val, end_val), dy_kwargs_qcd_chs),
                # (qgg.get_projection_plot(h2d_dyj_qcd_wrong_chs, start_val, end_val), dy_kwargs_qcd_wrong_chs),
                (qgg.get_projection_plot(h2d_qcd_wrong_chs, start_val, end_val), qcd_kwargs_wrong_chs),
                (qgg.get_projection_plot(h2d_qcd_chs, start_val, end_val), qcd_kwargs_chs),
            ]

            this_title = "%d < p_{T}^{jet} < %d GeV" % (start_val, end_val)
            if title != "":
                this_title += ", %s" % title
            qgg.do_comparison_plot(entries, "%s/%s/%s_pt%dto%d_flavMatched.%s" % (root_dir, plot_dir, v, start_val, end_val, OUTPUT_FMT),
                                   rebin=rebin, title=this_title, xlim=xlim)


if __name__ == '__main__':
    parser = qgc.get_parser()
    parser.add_argument('--title', help="Extra title to add to plots", default="")
    parser.add_argument('--reco', help="Do for reco jets", action='store_true')
    parser.add_argument('--gen', help="Do for gen jets", action='store_true')
    args = parser.parse_args()
    
    for workdir in args.workdirs:
        if args.reco:
            do_wrong_plots(workdir, title=args.title)
        if args.gen:
            do_wrong_plots(root_dir, var_prepend="gen", plot_dir="wrong_flavs_gen",
                           zpj_dirname=qgc.ZPJ_GENJET_RDIR, dj_dirname=qgc.DJ_GENJET_RDIR, pt_bins=qgc.THEORY_PT_BINS)
    