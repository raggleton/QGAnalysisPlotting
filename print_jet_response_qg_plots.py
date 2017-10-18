#!/usr/bin/env python

"""Print plots of jet response (reco pt / gen pt) split by pt, flavour, etc"""


import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
from itertools import product

# My stuff
from comparator import Contribution, Plot, grab_obj
import qg_common as qgc
import qg_general_plots as qgg
import common_utils as cu


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)

# Control output format
OUTPUT_FMT = "pdf"


# 2D plot, with various normalisations, logZ settings, box n whiskers plot
def do_all_2D_plots(root_dir, plot_dir="plots_2d", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG", var_list=None, var_prepend=""):
    for rn in ['Y', None]:
        pass
        # ZPJ
        # ALL
        # Q
        # G
        
        # DJ
        # ALL
        # Q
        # G
        


def do_projection_plots(root_dirs, plot_dir="response_plots", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG", flav_matched=False):
    """Plot both/either q and g jet repsonse, for various genjet pt bins. Contributions from all root_dirs are shown on same plot.
    
    flav_matched : If True, use hists with explicit gen-parton flavour matching
    """
    pt_bins = [(30, 40), (50, 70), (90, 100), (150, 200), (250, 300), (500, 600), (1000, 2000)]
    for (pt_min, pt_max) in pt_bins:
        
        lw = 2

        entries = []

        for root_dir in root_dirs:
            if zpj_dirname:
                flav_str = "q" if flav_matched else ""
                h2d_dyj = grab_obj(os.path.join(root_dir, qgc.DY_FILENAME), "%s/%s" % (zpj_dirname, "%sjet_response_vs_genjet_pt" % (flav_str)))
                obj = qgg.get_projection_plot(h2d_dyj, pt_min, pt_max)

                dy_reco_kwargs = dict(line_color=qgc.DY_COLOUR, fill_color=qgc.DY_COLOUR, line_width=lw,
                                       label=qgc.DY_ZpJ_QFLAV_LABEL if flav_matched else qgc.DY_ZpJ_LABEL)
                if len(root_dirs) > 1:
                    dy_reco_kwargs['label'] += " ["+root_dir+"]"
                entries.append((obj, dy_reco_kwargs))

            if dj_dirname:
                flav_str = "g" if flav_matched else ""
                h2d_qcd = grab_obj(os.path.join(root_dir, qgc.QCD_FILENAME), "%s/%s" % (dj_dirname, "%sjet_response_vs_genjet_pt" % (flav_str)))
                obj = qgg.get_projection_plot(h2d_qcd, pt_min, pt_max)

                qcd_reco_kwargs = dict(line_color=qgc.QCD_COLOUR, fill_color=qgc.QCD_COLOUR, line_width=lw,
                                       label=qgc.QCD_Dijet_GFLAV_LABEL if flav_matched else qgc.QCD_Dijet_LABEL)
                if len(root_dirs) > 1:
                    qcd_reco_kwargs['label'] += " ["+root_dir+"]"
                entries.append((obj, qcd_reco_kwargs))

        flav_str = "_flavMatched" if flav_matched else ""
        output_filename = os.path.join(plot_dir, "jet_response_ptGen%dto%d%s.%s" % (pt_min, pt_max, flav_str, OUTPUT_FMT))
        plot = qgg.make_comparison_plot_ingredients(entries, rebin=1, 
                                                    title="%d < p_{T}^{GenJet} < %d GeV" % (pt_min, pt_max),
                                                    xtitle="Response (p_{T}^{Reco} / p_{T}^{Gen})", xlim=(0, 3))
        plot.plot("HISTE NOSTACK")
        max_y = plot.container.GetHistogram().GetMaximum()
        line = ROOT.TLine(1, 0, 1, max_y)
        line.SetLineWidth(2)
        line.SetLineStyle(2)
        line.SetLineColor(13)
        line.Draw()
        plot.save(output_filename)


if __name__ == '__main__':
    parser = qgc.get_parser()
    args = parser.parse_args()

    if args.output is None:
        args.output = args.workdirs[0]

    # Do 2D plots
    for workdir in args.workdirs:
        do_all_2D_plots(workdir)

    # Do 1D comparison plots, without and with flavour matching
    app = "_comparison" if len(args.workdirs) > 1 else ""
    do_projection_plots(args.workdirs, plot_dir=os.path.join(args.output, "response_plots%s" % (app)))
    do_projection_plots(args.workdirs, plot_dir=os.path.join(args.output, "response_plots%s" % (app)), flav_matched=True)
