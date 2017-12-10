#!/usr/bin/env python

"""Print plots of jet response (reco pt / gen pt) split by pt, flavour, etc"""


import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
from itertools import product
from array import array

# My stuff
from comparator import Contribution, Plot, grab_obj
import qg_common as qgc
import qg_general_plots as qgg
import common_utils as cu


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)

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
        h2d_qcd.SetLineColor(qgc.QCD_COLOUR)
        h2d_qcd.SetFillStyle(0)
        draw_opt = "CANDLEY(00000311)"
        h2d_qcd.Draw(draw_opt)
        h2d_qcd.GetXaxis().SetRangeUser(0.7, 2.5)
        h2d_dyj.SetTitle("Dijet & Z+Jets%s" % (" (flav matched)" if flav_matched else ""))
        h2d_dyj.SetLineColor(qgc.DY_COLOUR)
        h2d_dyj.SetFillStyle(0)
        h2d_dyj.Draw(draw_opt + " SAME")
        h2d_qcd.GetYaxis().SetTitleOffset(1.7)
        h2d_qcd.GetXaxis().SetTitleOffset(1.2)
        leg = ROOT.TLegend(0.5, 0.6, 0.88, 0.88)
        leg.AddEntry(h2d_dyj, qgc.DY_ZpJ_QFLAV_LABEL if flav_matched else qgc.DY_ZpJ_LABEL, "L")
        leg.AddEntry(h2d_qcd, qgc.QCD_Dijet_GFLAV_LABEL if flav_matched else qgc.QCD_Dijet_LABEL, "L")
        leg.Draw()
        line.SetLineStyle(1)
        line.Draw()
        output_filename = os.path.join(plot_dir, "response_vs_genjet_pt_both%s_box.%s" % ("_flavMatched" if flav_matched else "", OUTPUT_FMT))
        canvas.SaveAs(output_filename)


def determine_fit_range(hist):
    """Determine lower & upper limit of fit range
    
    Parameters
    ----------
    hist : TH1
        Description
    
    Returns
    -------
    tuple
        (lower limit, upper limit)
    """
    # mean = hist.GetMean()
    mean = hist.GetBinCenter(hist.GetMaximumBin())
    rms = hist.GetRMS()
    return (mean - 1.*rms, mean + 1.*rms)


def do_gaus_fit(hist):
    """Do a Gaussian fit to histogram

    Parameters
    ----------
    hist : TH1
        Histogram to fit to
    """
    func_name = hist.GetName()+"_f1"
    func_name = "gausFit"
    fit_range = determine_fit_range(hist)
    func = ROOT.TF1(func_name, "gaus", fit_range[0], fit_range[1])
    # func.SetParameters(hist.GetMaximum(), hist.GetMean(), hist.GetRMS())
    fit_result = hist.Fit(func_name, "ERS", "L")
    print "fit result:", int(fit_result)


def fit_results_to_str(fit):
    """Turn fit results into str, lines split by \n
    
    Parameters
    ----------
    fit : TF1
        Description
    
    Returns
    -------
    str
        Description
    
    """
    parts = []
    chi2 = fit.GetChisquare()
    ndf = fit.GetNDF()
    if ndf > 0:
        parts.append("chi2/ndof: %.3e/%d = %.3e" % (chi2, ndf, chi2/ndf))
    else:
        parts.append("chi2/ndof: %.3e/0 = Inf" % (chi2))
    parts.append("prob: %.3e" % fit.GetProb())
    for i in range(fit.GetNpar()):
        parts.append("%s: %.3e #pm %.3e" % (fit.GetParName(i), fit.GetParameter(i), fit.GetParError(i)))
    return "\n".join(parts)


def fits_to_graph(pt_bins, fits):
    """Summary
    
    Parameters
    ----------
    pt_bins : [(float, float)]
        pt bin edges
    fits : TYPE
        Description
    
    Returns
    -------
    TGraphErrors
        Description
    """
    n_fits = len(fits)
    x, ex = [], []
    for pt_min, pt_max in pt_bins:
        center = 0.5*(pt_min+pt_max)
        x.append(center)
        half_width = pt_max - center
        ex.append(half_width)

    y, ey = [], []
    for fit in fits:
        for i in range(fit.GetNpar()):
            if fit.GetParName(i) == "Mean":
                y.append(fit.GetParameter(i))
                ey.append(fit.GetParError(i))
    gr = ROOT.TGraphErrors(n_fits, array('d', x), array('d', y), array('d', ex), array('d', ey))
    return gr


def do_response_graph(pt_bins, zpj_fits, dj_fits, title, output_filename):
    """Create and plot graph from fit results
    
    Parameters
    ----------
    zpj_fits : TF1
        Description
    dj_fits : TF1
        Description
    output_filename : str
        Description
    """
    gr_zpj = fits_to_graph(pt_bins, zpj_fits)
    gr_qcd = fits_to_graph(pt_bins, dj_fits)
    conts = [
        Contribution(gr_zpj, label=qgc.DY_ZpJ_LABEL, line_color=qgc.DY_COLOUR, marker_color=qgc.DY_COLOUR, marker_style=22),
        Contribution(gr_qcd, label=qgc.QCD_Dijet_LABEL, line_color=qgc.QCD_COLOUR, marker_color=qgc.QCD_COLOUR, marker_style=23)
    ]
    xmin = pt_bins[0][0]
    xmax = pt_bins[-1][1]
    plot = Plot(conts, what="graph", legend=True, xlim=(xmin, xmax),
                title=title,
                xtitle="p_{T}^{GenJet} [GeV]", 
                ytitle="Mean fitted response #pm fit error")
    plot.plot("ALP")
    plot.set_logx()
    line_center = ROOT.TLine(xmin, 1, xmax, 1)
    line_center.SetLineStyle(2)
    line_center.Draw("SAME")
    line_upper = ROOT.TLine(xmin, 1.1, xmax, 1.1)
    line_upper.SetLineStyle(2)
    line_upper.Draw("SAME")
    line_lower = ROOT.TLine(xmin, 0.9, xmax, 0.9)
    line_lower.SetLineStyle(2)
    line_lower.Draw("SAME")
    plot.save(output_filename)


def do_projection_plots(root_dirs, plot_dir="response_plots", zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG", flav_matched=False, do_fits=False):
    """Plot both/either q and g jet repsonse, for various genjet pt bins. Contributions from all root_dirs are shown on same plot.
    
    flav_matched : If True, use hists with explicit gen-parton flavour matching
    plot_dir : output dir for plots
    do_fits : If Ture, try and do a Gaussian fit to each 1D hist
    """
    pt_bins = [(20, 40), (40, 60), (60, 80), (80, 100), (100, 120), (120, 140), (160, 200), (200, 240), (240, 300), (300, 400), (400, 500), (500, 600), (600, 800), (800, 1000), (1000, 2000)]
    print_pt_bins = [(20, 40), (40, 60), (60, 80), (100, 120), (160, 200), (240, 300), (400, 500), (500, 600), (800, 1000), (1000, 2000)]
    print_pt_bins = pt_bins
    zpj_fits = []
    dj_fits = []
    for (pt_min, pt_max) in pt_bins:
        print pt_min, pt_max
        if pt_min > pt_max:
            raise RuntimeError("pt_min < pt_max!")

        lw = 2 if len(root_dirs) == 1 else 1

        plot_entries = []
        do_plot = (pt_min, pt_max) in print_pt_bins
        rebin = 1

        for ind, root_dir in enumerate(root_dirs):
            if zpj_dirname:
                flav_str = "q" if flav_matched else ""
                h2d_dyj = grab_obj(os.path.join(root_dir, qgc.DY_FILENAME), "%s/%s" % (zpj_dirname, "%sjet_response_vs_genjet_pt" % (flav_str)))
                obj = qgg.get_projection_plot(h2d_dyj, pt_min, pt_max)
                if obj.Integral() > 0:
                    obj.Rebin(rebin)
                    obj.Scale(1./obj.Integral())
                    col = qgc.DY_COLOURS[ind]
                    dy_reco_kwargs = dict(line_color=col, fill_color=col, line_width=lw, marker_color=col, marker_style=qgc.DY_MARKER,
                                          label=qgc.DY_ZpJ_QFLAV_LABEL if flav_matched else qgc.DY_ZpJ_LABEL)
                    if len(root_dirs) > 1:
                        dy_reco_kwargs['label'] += " ["+root_dir+"]"

                    if do_fits:
                        do_gaus_fit(obj)
                        fit = obj.GetFunction("gausFit")
                        dy_reco_kwargs['label'] += "\n"
                        dy_reco_kwargs['label'] += fit_results_to_str(fit)
                        zpj_fits.append(fit)

                    if do_plot:
                        plot_entries.append((obj, dy_reco_kwargs))

            if dj_dirname:
                flav_str = "g" if flav_matched else ""
                h2d_qcd = grab_obj(os.path.join(root_dir, qgc.QCD_FILENAME), "%s/%s" % (dj_dirname, "%sjet_response_vs_genjet_pt" % (flav_str)))
                obj = qgg.get_projection_plot(h2d_qcd, pt_min, pt_max)
                if obj.Integral() > 0:
                    obj.Rebin(rebin)
                    obj.Scale(1./obj.Integral())
                    col = qgc.QCD_COLOURS[ind]
                    qcd_reco_kwargs = dict(line_color=col, fill_color=col, line_width=lw, marker_color=col, marker_style=qgc.QCD_MARKER,
                                           label=qgc.QCD_Dijet_GFLAV_LABEL if flav_matched else qgc.QCD_Dijet_LABEL)
                    if len(root_dirs) > 1:
                        qcd_reco_kwargs['label'] += " ["+root_dir+"]"

                    if do_fits:
                        do_gaus_fit(obj)
                        fit = obj.GetFunction("gausFit")
                        qcd_reco_kwargs['label'] += "\n"
                        qcd_reco_kwargs['label'] += fit_results_to_str(fit)
                        dj_fits.append(fit)

                    if do_plot:
                        plot_entries.append((obj, qcd_reco_kwargs))

        flav_str = "_flavMatched" if flav_matched else ""

        if do_plot:
            output_filename = os.path.join(plot_dir, "jet_response_ptGen%dto%d%s.%s" % (pt_min, pt_max, flav_str, OUTPUT_FMT))
            plot = qgg.make_comparison_plot_ingredients(plot_entries, rebin=1, normalise_hist=False,
                                                    title="%d < p_{T}^{GenJet} < %d GeV" % (pt_min, pt_max),
                                                    xtitle="Response (p_{T}^{Reco} / p_{T}^{Gen})", xlim=(0.25, 1.75))
            plot.legend.SetX1(0.65)
            plot.legend.SetX2(0.95)
            plot.legend.SetY1(0.5)
            plot.legend.SetY2(0.89)
            plot.plot("E NOSTACK")
            max_y = plot.container.GetHistogram().GetMaximum()
            line = ROOT.TLine(1, 0, 1, max_y)
            line.SetLineWidth(2)
            line.SetLineStyle(2)
            line.SetLineColor(13)
            line.Draw()
            plot.save(output_filename)

    do_response_graph(pt_bins, zpj_fits, dj_fits, title="Flav Matched" if flav_matched else "Not flav matched",
                      output_filename=os.path.join(plot_dir, "gr_jet_response%s.%s" % (flav_str, OUTPUT_FMT)))


if __name__ == '__main__':
    parser = qgc.get_parser()
    parser.add_argument("--comparison", action='store_true', 
                        help="Only do comparison plots. Otherwise only do 2D plots for each dir.")
    args = parser.parse_args()
    print args

    if not args.comparison:
        # Do 2D plots
        for workdir in args.workdirs:
            do_all_2D_plots(workdir, plot_dir=os.path.join(workdir, "response_2d"))
            do_projection_plots([workdir], plot_dir=os.path.join(workdir, "response_plots"), do_fits=True)
            do_projection_plots([workdir], plot_dir=os.path.join(workdir, "response_plots"), do_fits=True, flav_matched=True)
    else:
        if args.output is None:
            app = "_comparison" if len(args.workdirs) > 1 else ""
            args.output = os.path.join(args.workdirs[0], "response_plots%s" % (app))
        # Do 1D comparison plots, without and with flavour matching
        do_projection_plots(args.workdirs, plot_dir=args.output)
        do_projection_plots(args.workdirs, plot_dir=args.output, flav_matched=True)
