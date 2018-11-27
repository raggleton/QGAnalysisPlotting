#!/usr/bin/env python

"""Print plots of gen vs reco response plots"""


import argparse
# from MyStyle import My_Style
# My_Style.cd()
import os
from itertools import product
from array import array
from math import sqrt

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)


# My stuff
from comparator import Contribution, Plot
import qg_common as qgc
import qg_general_plots as qgg
import common_utils as cu

ROOT.gStyle.SetPaintTextFormat(".3f")

# Control output format
OUTPUT_FMT = "pdf"


def do_response_plot(tdir, plot_dir, var_name, xlabel, log_var=False, rebinx=1, rebiny=1, do_migration_summary_plots=True, do_resolution_plots=True, save_response_hists=False):
    """Do 2D plots of genjet pt vs recojet pt"""
    h2d_orig = cu.get_from_tfile(tdir, var_name)
    h2d = h2d_orig.Clone(cu.get_unique_str())
    h2d.RebinX(rebinx)
    h2d.RebinY(rebiny)

    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    draw_opt = "COLZ TEXT"
    draw_opt = "COLZ"
    pad = ROOT.gPad
    pad.SetBottomMargin(0.12)
    pad.SetLeftMargin(0.13)
    pad.SetRightMargin(0.12)
    canv.SetTicks(1, 1)
    if log_var:
        canv.SetLogx()
        canv.SetLogy()

    # un normalised
    h2d.Draw(draw_opt)
    xax = h2d.GetXaxis()
    upper_lim = xax.GetBinUpEdge(xax.GetLast())
    # print(upper_lim)
    # upper_lim = 5000
    title_offset = 1.5
    h2d.SetTitleOffset(title_offset, 'X')
    xax.SetMoreLogLabels()

    yax = h2d.GetYaxis()
    h2d.SetTitleOffset(title_offset*1.15, 'Y')
    yax.SetMoreLogLabels()
    canv.Update()

    plot_dir = os.path.join(plot_dir, tdir.GetName())
    cu.check_dir_exists_create(plot_dir)
    canv.SaveAs(os.path.join(plot_dir, "%s_linZ.%s" % (var_name, OUTPUT_FMT)))
    canv.SetLogz()
    canv.SaveAs(os.path.join(plot_dir, "%s_logZ.%s" % (var_name, OUTPUT_FMT)))

    #  renorm by row (reco bins)
    canv.SetLogz(0)
    h2d_renorm_y = cu.make_normalised_TH2(h2d, 'Y', recolour=False, do_errors=True)
    h2d_renorm_y.SetMarkerSize(0.5)
    h2d_renorm_y.SetMaximum(1)
    h2d_renorm_y.Draw(draw_opt)
    xax = h2d_renorm_y.GetXaxis()
    upper_lim = xax.GetBinUpEdge(xax.GetLast())
    print(upper_lim)
    # upper_lim = 5000
    title_offset = 1.5
    h2d_renorm_y.SetTitleOffset(title_offset, 'X')
    # xax.SetRangeUser(0, upper_lim)
    xax.SetMoreLogLabels()

    yax = h2d_renorm_y.GetYaxis()
    h2d_renorm_y.SetTitleOffset(title_offset*1.15, 'Y')
    # yax.SetRangeUser(0, upper_lim)
    yax.SetMoreLogLabels()
    canv.Update()

    canv.SaveAs(os.path.join(plot_dir, "%s_renormY_linZ.%s" % (var_name, OUTPUT_FMT)))
    canv.SetLogz()
    h2d_renorm_y.SetMaximum(1)
    h2d_renorm_y.SetMinimum(1E-3)
    canv.SaveAs(os.path.join(plot_dir, "%s_renormY_logZ.%s" % (var_name, OUTPUT_FMT)))

    # renorm by column (gen bins)
    canv.Clear()
    canv.SetLogz(0)
    h2d_renorm_x = cu.make_normalised_TH2(h2d, 'X', recolour=False, do_errors=True)
    h2d_renorm_x.SetMarkerSize(0.5)
    h2d_renorm_x.SetMaximum(1)
    h2d_renorm_x.Draw(draw_opt)

    xax = h2d_renorm_x.GetXaxis()
    upper_lim = xax.GetBinUpEdge(xax.GetLast())
    # upper_lim = 5000
    title_offset = 1.5
    h2d_renorm_x.SetTitleOffset(title_offset, 'X')
    # xax.SetRangeUser(0, upper_lim)
    xax.SetMoreLogLabels()

    yax = h2d_renorm_x.GetYaxis()
    h2d_renorm_x.SetTitleOffset(title_offset*1.15, 'Y')
    # yax.SetRangeUser(0, upper_lim)
    yax.SetMoreLogLabels()
    canv.Update()
    canv.SaveAs(os.path.join(plot_dir, "%s_renormX_linZ.%s" % (var_name, OUTPUT_FMT)))
    canv.SetLogz()
    h2d_renorm_x.SetMaximum(1)
    h2d_renorm_x.SetMinimum(1E-3)
    canv.SaveAs(os.path.join(plot_dir, "%s_renormX_logZ.%s" % (var_name, OUTPUT_FMT)))

    # Now do plot of purity, etc
    if do_migration_summary_plots:
        qgg.make_migration_summary_plot(h2d_renorm_x,
                                        h2d_renorm_y,
                                        xlabel=xlabel,
                                        output_filename=os.path.join(plot_dir, '%s_migration_summary.%s' % (var_name, OUTPUT_FMT)),
                                        log_var=log_var)

    # Do resolution plots
    if do_resolution_plots:
        res_rms, rel_res_rms = make_resolution_plots(h2d_orig,
                                                     xlabel=xlabel,
                                                     output_filename=os.path.join(plot_dir, '%s_resolution_summary_rms.%s' % (var_name, OUTPUT_FMT)),
                                                     do_fit=False,
                                                     do_rms=True,
                                                     log_var=log_var,
                                                     save_response_hists=False)

        res_quantiles, rel_res_quantiles = make_resolution_plots(h2d_orig,
                                                                 xlabel=xlabel,
                                                                 output_filename=os.path.join(plot_dir, '%s_resolution_summary_quantiles.%s' % (var_name, OUTPUT_FMT)),
                                                                 do_fit=False,
                                                                 do_rms=False,
                                                                 quantiles=None,  # auto determine
                                                                 log_var=log_var,
                                                                 save_response_hists=save_response_hists)

        # compare RMS and quantile results
        conts = [
            Contribution(res_rms, label="#sqrt{RMS} #pm #frac{#sqrt{#delta RMS}}{#LT %s #GT}" % (xlabel), line_color=ROOT.kRed, marker_color=ROOT.kRed),
            Contribution(res_quantiles, label="68% quantile", line_color=ROOT.kBlue, marker_color=ROOT.kBlue)
        ]
        ylabel = "#frac{#sigma(RECO/GEN)}{GEN}" if "_rel_" in var_name else "#frac{#sigma(RECO)}{GEN}"
        res_plot = Plot(conts, what='graph',
                        xtitle=xlabel,
                        ytitle=ylabel,
                        legend=True,
                        ylim=[0, 3]
                        )
        res_plot.legend.SetX1(0.5)
        res_plot.legend.SetX2(0.9)
        res_plot.legend.SetY1(0.7)
        res_plot.legend.SetY2(0.85)
        res_plot.plot('ALP')
        res_plot.save(os.path.join(plot_dir, '%s_resolution_rms_vs_quantiles.%s' % (var_name, OUTPUT_FMT)))


def do_jet_index_plots(tdir, plot_dir):
    """Do 2D plots of genjet index vs recojet index"""
    plot_dir = os.path.join(plot_dir, tdir.GetName())
    cu.check_dir_exists_create(plot_dir)
    stem = "genjet_ind_recojet_ind_pt_"
    for plot_name in cu.get_list_of_element_names(tdir):
        if not plot_name.startswith(stem):
            continue
        h2d = cu.get_from_tfile(tdir, plot_name)
        h2d.SetTitle(h2d.GetTitle())
        renorm_h2d = cu.make_normalised_TH2(h2d, 'X', recolour=False)
        # renorm_h2d = h2d
        canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
        pad = ROOT.gPad
        pad.SetBottomMargin(0.12)
        pad.SetLeftMargin(0.13)
        pad.SetRightMargin(0.12)
        canv.SetTicks(1, 1)
        renorm_h2d.Draw("COLZ TEXT")
        renorm_h2d.SetMaximum(1)
        renorm_h2d.SetMinimum(0)
        title_offset = 1.5
        renorm_h2d.SetTitleOffset(title_offset, 'X')
        renorm_h2d.SetTitleOffset(title_offset*1.15, 'Y')

        canv.SaveAs(os.path.join(plot_dir, "%s_renormX_linZ.%s" % (plot_name, OUTPUT_FMT)))
        # canv.SetLogz()
        # renorm_h2d.SetMinimum(1E-3)
        # canv.SaveAs(os.path.join(plot_dir, "%s_renormX_logZ.%s" % (plot_name, OUTPUT_FMT)))

        canv.Clear()
        renorm_h2d = cu.make_normalised_TH2(h2d, 'Y', recolour=False)
        renorm_h2d.Draw("COLZ TEXT")
        renorm_h2d.SetMaximum(1)
        renorm_h2d.SetMinimum(0)
        title_offset = 1.5
        renorm_h2d.SetTitleOffset(title_offset, 'X')
        renorm_h2d.SetTitleOffset(title_offset*1.15, 'Y')
        canv.SaveAs(os.path.join(plot_dir, "%s_renormY_linZ.%s" % (plot_name, OUTPUT_FMT)))


def do_pt_transfer_plot(tdir, plot_dir):
    """Plot ratio between pt bins of the spectrum. Check to make sure xfer factor << drop in pt"""
    plot_dir = os.path.join(plot_dir, tdir.GetName())
    cu.check_dir_exists_create(plot_dir)
    hist_name = "pt_jet_response_binning"
    h = cu.get_from_tfile(tdir, hist_name)
    binning = [h.GetXaxis().GetBinLowEdge(bin_ind) for bin_ind in range(1, h.GetNbinsX()+1)]
    hist_factors = ROOT.TH1F("hist_factors"+cu.get_unique_str(), ";p_{T}^{Reco} [GeV];Fraction rel to previous bin", len(binning)-1, array('d', binning))
    for bin_ind in range(2, h.GetNbinsX()+1):
        cont = h.GetBinContent(bin_ind)
        cont_prev = h.GetBinContent(bin_ind-1)
        if cont == 0 or cont_prev == 0:
            continue
        factor = cont / cont_prev
        hist_factors.SetBinContent(bin_ind, factor)
        hist_factors.SetBinError(bin_ind, 0)

    col_purity = ROOT.kBlack
    conts = [
        Contribution(hist_factors, label="Factor relative to previous bin", line_color=col_purity, marker_color=col_purity),
        # Contribution(hist_purity, label="Purity (gen in right bin)", line_color=col_purity, marker_color=col_purity),
    ]
    xlim = [30, binning[-1]]
    plot = Plot(conts,  what='hist', xlim=xlim)
    plot.plot()
    plot.set_logx()
    plot.save(os.path.join(plot_dir, 'pt_migration_factors.%s' % (OUTPUT_FMT)))


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
    multiplier = 1.5
    return (mean - multiplier*rms, mean + multiplier*rms)


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
    fit_result = hist.Fit(func_name, "ERSQ", "L")
    # print("fit result:", int(fit_result))


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


def make_resolution_plots(h2d, xlabel, output_filename, do_fit=True, do_rms=True, quantiles=None, log_var=False, save_response_hists=False):
    """Make graph of resolution vs variables.

    Also optionally save all input histograms to file.
    """
    one_sigma = 0.682689
    quantiles = quantiles or [0.5*(1-one_sigma), 1 - 0.5*(1-one_sigma)]
    ax = h2d.GetXaxis()
    bin_edges = [ax.GetBinLowEdge(i) for i in range(1, ax.GetNbins()+1)]

    bin_centers, sigmas, sigmas_unc = [], [], []
    rel_sigmas, rel_sigmas_unc = [], []
    # bin_centers = [ax.GetBinCenter(i) for i in range(1, ax.GetNbins()+1)]

    for var_min, var_max in zip(bin_edges[:-1], bin_edges[1:]):
        h_projection = qgg.get_projection_plot(h2d, var_min, var_max, cut_axis='x')
        if h_projection.GetEffectiveEntries() < 20:
            continue
        # h_projection.Rebin(rebin)
        h_projection.Scale(1./h_projection.Integral())

        bin_centers.append(0.5*(var_max+var_min))
        if do_fit:
            do_gaus_fit(h_projection)
            fit = h_projection.GetFunction("gausFit")
            # label += "\n"
            # label += fit_results_to_str(fit)
            # bin_centers.append(fit.GetParameter(1))
            sigmas.append(fit.GetParameter(2))
            rel_sigmas.append(fit.GetParameter(2)/bin_centers[-1])
            sigmas_unc.append(fit.GetParError(2))
            rel_sigmas_unc.append(fit.GetParError(2)/bin_centers[-1])
        else:
            if do_rms:
                sigmas.append(sqrt(h_projection.GetRMS()))
                rel_sigmas.append(sqrt(h_projection.GetRMS())/bin_centers[-1])
                sigmas_unc.append(sqrt(h_projection.GetRMSError()))
                rel_sigmas_unc.append(sqrt(h_projection.GetRMSError())/bin_centers[-1])
            elif quantiles:
                if len(quantiles) != 2:
                    raise RuntimeError("Need 2 quantiles")
                q = array('d', quantiles)
                results = array('d', [0.]*len(quantiles))
                h_projection.GetQuantiles(len(quantiles), results, q)
                sigmas.append(results[1] - results[0])
                sigmas_unc.append(0)
                rel_sigmas.append((results[1] - results[0])/bin_centers[-1])
                rel_sigmas_unc.append(0)
            else:
                raise RuntimeError("Need either do_fit, do_rms, or 2-tuple in quantiles")


        if save_response_hists:
            xlabel =  h_projection.GetXaxis().GetTitle()
            cont = Contribution(h_projection, label="GEN: %g-%g" % (var_min, var_max))
            p = Plot([cont], what='hist')
            p.plot('HISTE')
            rsp_filename = os.path.abspath(output_filename.replace(".%s" % OUTPUT_FMT, "_hist%gto%g.%s" % (var_min, var_max, OUTPUT_FMT)))
            rsp_dir = os.path.dirname(rsp_filename)
            rsp_file = os.path.basename(rsp_filename)
            p.save(os.path.join(rsp_dir, "responseHists", rsp_file))

    gr = ROOT.TGraphErrors(len(bin_centers), array('d', bin_centers), array('d', sigmas), array('d', [0]*len(bin_centers)), array('d', sigmas_unc))
    gr_cont = Contribution(gr, label="")
    ylabel = ""
    if do_fit:
        ylabel = "Fit #sigma"
    elif do_rms:
        ylabel = "#sqrt{RMS}"
    elif quantiles:
        ylabel = "Central %g" % one_sigma
    plot = Plot([gr_cont], what='graph', xtitle=xlabel, ytitle=ylabel, xlim=[bin_edges[0], bin_edges[-1]], ylim=[0, max(sigmas)*1.2], legend=False)
    plot.plot()
    if log_var:
        plot.set_logx()
    plot.save(output_filename)

    gr_rel = ROOT.TGraphErrors(len(bin_centers), array('d', bin_centers), array('d', rel_sigmas), array('d', [0]*len(bin_centers)), array('d', rel_sigmas_unc))
    gr_rel_cont = Contribution(gr_rel, label="")
    ylabel = "Relative %s" % ylabel
    plot = Plot([gr_rel_cont], what='graph', xtitle=xlabel, ytitle=ylabel, xlim=[bin_edges[0], bin_edges[-1]], ylim=[min(rel_sigmas)/1.2, max(rel_sigmas)*1.2], legend=False)
    plot.plot()
    plot.set_logy()
    if log_var:
        plot.set_logx()
    plot.save(output_filename.replace(".%s" % OUTPUT_FMT, "_relative.%s" % OUTPUT_FMT))

    return gr, gr_rel


def do_response_plots(in_file, plot_dir, do_these=None):
    tfile = cu.open_root_file(in_file)
    for full_var_name, xlabel, log_var, rebin in do_these:
        mydir, myvar = full_var_name.split("/")
        # reco vs gen
        do_response_plot(tfile.Get(mydir),
            plot_dir=plot_dir,
            var_name=myvar+"_response",
            xlabel=xlabel,
            log_var=log_var,
            rebinx=rebin,
            rebiny=rebin,
            do_migration_summary_plots=True,
            do_resolution_plots=True,
            save_response_hists=False
        )

        rebiny = 10 if "multiplicity" in myvar.lower() else 5

        # relative respone (reco/gen) on y axis
        do_response_plot(tfile.Get(mydir),
            plot_dir=plot_dir,
            var_name=myvar+"_rel_response",
            xlabel=xlabel,
            log_var=log_var,
            rebinx=rebin,
            rebiny=rebiny,
            do_migration_summary_plots=False,
            do_resolution_plots=True,
            save_response_hists=True
        )

        # do_jet_index_plots(tfile.Get(mydir), plot_dir=plot_dir)

        # do_pt_transfer_plot(tfile.Get(mydir), plot_dir=plot_dir)


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
        do_these = None

        if "QCD" in in_file:
            do_these = [
                # ("Dijet_tighter/pt_jet", "p_{T}^{Gen} [GeV]", True, 1),
                # ("Dijet_QG_tighter/jet_puppiMultiplicity", "PUPPI Multiplicity (#lambda_{0}^{0} (PUPPI))", False, 5),
                # ("Dijet_QG_tighter/jet_multiplicity", "Multiplicity (#lambda_{0}^{0})", False, 5),
                ("Dijet_QG_tighter/jet_LHA", "LHA (#lambda_{0.5}^{1})", False, 1),
                ("Dijet_QG_tighter/jet_pTD", "p_{T}^{D} (#lambda_{0}^{2})", False, 2),
                ("Dijet_QG_tighter/jet_width", "Width (#lambda_{1}^{1})", False, 2),
                ("Dijet_QG_tighter/jet_thrust", "Thrust (#lambda_{2}^{1})", False, 2),

                ("Dijet_QG_tighter/jet_puppiMultiplicity_charged", "PUPPI Multiplicity (#lambda_{0}^{0} (PUPPI)) [charged]", False, 5),
                # ("Dijet_QG_tighter/jet_multiplicity_charged", "Multiplicity (#lambda_{0}^{0}) [charged only]", False, 5),
                ("Dijet_QG_tighter/jet_LHA_charged", "LHA (#lambda_{0.5}^{1}) [charged only]", False, 2),
                ("Dijet_QG_tighter/jet_pTD_charged", "p_{T}^{D} (#lambda_{0}^{2}) [charged only]", False, 2),
                ("Dijet_QG_tighter/jet_width_charged", "Width (#lambda_{1}^{1}) [charged only]", False, 2),
                ("Dijet_QG_tighter/jet_thrust_charged", "Thrust (#lambda_{2}^{1}) [charged only]", False, 2),
            ][:1]

        if "DYJetsToLL" in in_file:
            pass

        do_response_plots(in_file, plot_dir=plot_dir, do_these=do_these)
