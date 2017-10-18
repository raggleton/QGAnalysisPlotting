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
def do_all_2D_plots(root_dir, plot_dir="plots_2d", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG"):
    """Do 2D plots of response vs pt, with various normalisations, etc
    
    plot_dir : output dir for plots
    """
    line = ROOT.TLine(1, 0, 1, 2000)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.SetLineColor(13)

    for flav_matched in [True, False]:
        for renorm, logz in product(['Y', None], [True, False]):

            log_append = "_logZ" if logz else ""


            if zpj_dirname:
                flav_str = "q" if flav_matched else ""
                h2d_dyj = grab_obj(os.path.join(root_dir, qgc.DY_FILENAME), "%s/%s" % (zpj_dirname, "%sjet_response_vs_genjet_pt" % (flav_str)))
                flav_append = "_%sflavMatched" % flav_str if flav_matched else ""
                output_filename = os.path.join(plot_dir, "response_vs_genjet_pt_zpj%s_norm%s%s.%s" % (flav_append, renorm, log_append, OUTPUT_FMT))
                title = qgc.DY_ZpJ_LABEL
                if flav_matched:
                    title += " (%s-matched)" % flav_str
                zlim = [1E-5, 1] if logz and renorm else None
                qgg.do_2D_plot(h2d_dyj, output_filename, renorm_axis=renorm, title=title, rebin=None, recolour=True, xlim=None, ylim=None, zlim=zlim, logz=logz, other_things_to_draw=[line])
                # output_filename = os.path.join(plot_dir, "response_vs_genjet_pt_zpj%s_norm%s%s_box.%s" % (flav_append, renorm, log_append, OUTPUT_FMT))
                # qgg.do_2D_plot(h2d_dyj, output_filename, draw_opt="CANDLEY1", renorm_axis=renorm, title=title, rebin=None, recolour=True, xlim=None, ylim=None, zlim=zlim, logz=logz, other_things_to_draw=[line])

            if dj_dirname:
                flav_str = "g" if flav_matched else ""
                h2d_qcd = grab_obj(os.path.join(root_dir, qgc.QCD_FILENAME), "%s/%s" % (dj_dirname, "%sjet_response_vs_genjet_pt" % (flav_str)))
                flav_append = "_%sflavMatched" % flav_str if flav_matched else ""
                output_filename = os.path.join(plot_dir, "response_vs_genjet_pt_dj%s_norm%s%s.%s" % (flav_append, renorm, log_append, OUTPUT_FMT))
                title = qgc.QCD_Dijet_LABEL
                if flav_matched:
                    title += " (%s-matched)" % flav_str
                zlim = [1E-5, 1] if logz and renorm else None
                qgg.do_2D_plot(h2d_qcd, output_filename, renorm_axis=renorm, title=title, rebin=None, recolour=True, xlim=None, ylim=None, zlim=zlim, logz=logz, other_things_to_draw=[line])
                # output_filename = os.path.join(plot_dir, "response_vs_genjet_pt_dj%s_norm%s%s_box.%s" % (flav_append, renorm, log_append, OUTPUT_FMT))
                # qgg.do_2D_plot(h2d_qcd, output_filename, draw_opt="CANDLEY1", renorm_axis=renorm, title=title, rebin=None, recolour=True, xlim=None, ylim=None, zlim=zlim, logz=logz, other_things_to_draw=[line])

        # Do box plot comparing them both
        canvas = ROOT.TCanvas(ROOT.TUUID().AsString(), "", 800, 800)
        canvas.SetTicks(1, 1)
        canvas.SetLeftMargin(0.13)
        canvas.SetBottomMargin(0.11)
        h2d_qcd.SetTitle("Dijet & Z+Jets%s" % (" (flav matched)" if flav_matched else ""))
        h2d_qcd.SetLineColor(ROOT.kBlue)
        h2d_qcd.SetFillStyle(0)
        draw_opt = "CANDLEY(00000311)"
        h2d_qcd.Draw(draw_opt)
        h2d_qcd.GetXaxis().SetRangeUser(0.7, 2.5)
        h2d_dyj.SetLineColor(ROOT.kRed)
        h2d_dyj.SetFillStyle(0)
        h2d_dyj.Draw(draw_opt + " SAME")
        h2d_qcd.GetYaxis().SetTitleOffset(1.7)
        h2d_qcd.GetXaxis().SetTitleOffset(1.2)
        leg = ROOT.TLegend(0.5, 0.6, 0.88, 0.88)
        leg.AddEntry(h2d_qcd, qgc.DY_ZpJ_QFLAV_LABEL if flav_matched else qgc.DY_ZpJ_LABEL, "L")
        leg.AddEntry(h2d_dyj, qgc.QCD_Dijet_GFLAV_LABEL if flav_matched else qgc.QCD_Dijet_LABEL, "L")
        leg.Draw()
        line.SetLineStyle(1)
        line.Draw()
        output_filename = os.path.join(plot_dir, "response_vs_genjet_pt_both%s_box.%s" % ("_flavMatched" if flav_matched else "", OUTPUT_FMT))
        canvas.SaveAs(output_filename)



def do_projection_plots(root_dirs, plot_dir="response_plots", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG", flav_matched=False):
    """Plot both/either q and g jet repsonse, for various genjet pt bins. Contributions from all root_dirs are shown on same plot.
    
    flav_matched : If True, use hists with explicit gen-parton flavour matching
    plot_dir : output dir for plots
    """
    pt_bins = [(20, 40), (60, 80), (100, 120), (160, 200), (260, 300), (500, 600), (1000, 2000)]
    for (pt_min, pt_max) in pt_bins:
        
        lw = 2 if len(root_dirs) == 1 else 1

        entries = []

        for ind, root_dir in enumerate(root_dirs):
            if zpj_dirname:
                flav_str = "q" if flav_matched else ""
                h2d_dyj = grab_obj(os.path.join(root_dir, qgc.DY_FILENAME), "%s/%s" % (zpj_dirname, "%sjet_response_vs_genjet_pt" % (flav_str)))
                obj = qgg.get_projection_plot(h2d_dyj, pt_min, pt_max)

                col = qgc.DY_COLOURS[ind]
                dy_reco_kwargs = dict(line_color=col, fill_color=col, line_width=lw,
                                      label=qgc.DY_ZpJ_QFLAV_LABEL if flav_matched else qgc.DY_ZpJ_LABEL)
                if len(root_dirs) > 1:
                    dy_reco_kwargs['label'] += " ["+root_dir+"]"
                entries.append((obj, dy_reco_kwargs))

            if dj_dirname:
                flav_str = "g" if flav_matched else ""
                h2d_qcd = grab_obj(os.path.join(root_dir, qgc.QCD_FILENAME), "%s/%s" % (dj_dirname, "%sjet_response_vs_genjet_pt" % (flav_str)))
                obj = qgg.get_projection_plot(h2d_qcd, pt_min, pt_max)

                col = qgc.QCD_COLOURS[ind]
                qcd_reco_kwargs = dict(line_color=col, fill_color=col, line_width=lw,
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
        app = "_comparison" if len(args.workdirs) > 1 else ""
        args.output = os.path.join(args.workdirs[0], "response_plots%s" % (app))

    # Do 2D plots
    for workdir in args.workdirs:
        do_all_2D_plots(workdir, plot_dir=os.path.join(workdir, "response_2d"))

    # Do 1D comparison plots, without and with flavour matching
    do_projection_plots(args.workdirs, plot_dir=args.output)
    do_projection_plots(args.workdirs, plot_dir=args.output, flav_matched=True)
