"""Common functions for doing simple QG plots"""


import ROOT
from MyStyle import My_Style
My_Style.cd()
import bisect
import os
from array import array

# My stuff
from comparator import Contribution, Plot, grab_obj
import common_utils as cu
import qg_common as qgc

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.patches as patches


mpl.rcParams['font.size'] = 18
mpl.rcParams["font.family"] = "arial"


# For debugging
import sys
import tracers
# sys.settrace(tracers.trace_calls)
# sys.settrace(tracers.trace_calls_detail)
# sys.settrace(tracers.trace_calls_and_returns)


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)


def make_comparison_plot_ingredients(entries, rebin=1, normalise_hist=True, mean_rel_error=1.0, **plot_kwargs):
    """Make the Plot object for a comparison plot.

    User can then add other elements to the plot.

    Parameters
    ----------
    entries : list[(object, dict)]
        List of ROOT object & it's configuration dict, where the dict is a set of kwargs passed to the Contribution object
    rebin : int, optional
        Rebin factor
    normalise_hist : bool, optional
        Normalise each histogram's integral to unity
    mean_rel_error : float, optional
        Remove contributions that have a mean realtive error more than this value
        mean realtive error = mean of all the error/bin contents
    **plot_kwargs
        Any other kwargs to be passed to the Plot object contructor

    Returns
    -------
    Plot
        Plot object to be modified, plotted, etc

    Raises
    ------
    RuntimeError
        If there are 0 contributions
    """
    conts = [Contribution(ent[0], normalise_hist=normalise_hist, rebin_hist=rebin, **ent[1])
             for ent in entries
             if cu.get_hist_mean_rel_error(ent[0]) < mean_rel_error and ent[0].Integral() > 0]
    do_legend = len(conts) > 1
    if len(conts) == 0:
        raise RuntimeError("0 contributions for this plot")
    do_subplot = any(c.subplot for c in conts)
    if (len(conts) == 1 or not do_subplot) and "subplot_type" in plot_kwargs:
        plot_kwargs['subplot_type'] = None
    p = Plot(conts, what="hist", ytitle="p.d.f", legend=do_legend, **plot_kwargs)
    if do_legend:
        p.legend.SetX1(0.5)
        p.legend.SetX2(1.02)
        if len(entries) > 4:
            p.legend.SetY1(0.6)
        else:
            p.legend.SetY1(0.65)
        p.legend.SetY2(0.88)
    return p


def do_comparison_plot(entries, output_filename, rebin=1, **plot_kwargs):
    """Plot several different objects on a single plot

    entries : list of 2-tuples, with (object, dict), where the dict is a set of kwargs passed to the Contribution object
    plot_kwargs : any other kwargs to be passed to the Plot object ctor
    """
    try:
        p = make_comparison_plot_ingredients(entries, rebin=rebin, mean_rel_error=0.4, **plot_kwargs)
        draw_opt = "NOSTACK HISTE"
        p.plot(draw_opt)
        # p.container.SetMaximum(min(5, p.container.GetMaximum()))
        dirname = os.path.dirname(output_filename)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
        p.save(output_filename)
    except RuntimeError as e:
        print("Skipping")
        print(e)


def get_projection_plot(h2d, start_val, end_val, cut_axis='y'):
    cut_axis = cut_axis.lower()
    if cut_axis == "y":
        axis = h2d.GetYaxis()
    elif cut_axis == "x":
        axis = h2d.GetXaxis()
    else:
        raise RuntimeError("cut_axis must be x or y")
    bin_edges = [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins()+2)]
    bin_start = bisect.bisect_right(bin_edges, start_val)
    bin_end = bisect.bisect_right(bin_edges, end_val)
    if cut_axis == "y":
        hproj = h2d.ProjectionX(ROOT.TUUID().AsString(), bin_start, bin_end, "eo")
        return hproj
    elif cut_axis == "x":
        hproj = h2d.ProjectionY(ROOT.TUUID().AsString(), bin_start, bin_end, "eo")
        return hproj


def do_2D_plot(obj, output_filename, draw_opt="COLZ", renorm_axis=None, title=None, rebin=None, recolour=True, xlim=None, ylim=None, zlim=None, logz=False, other_things_to_draw=None):
    """Print a 2D hist to file.

    renorm_axis : [None, "X", "Y"]
        Allow renormalising along a given axis

    recolour : bool
        Only used if renorm_axis != None
        If True, resets colour so max bin in each row/col (depending on renorm_axis) has value 1.

    other_things_to_draw : list of other objects that will be drawn after the hist

    """
    if rebin:
        obj.Rebin2D(*rebin)
    if renorm_axis:
        obj_renorm = cu.make_normalised_TH2(obj, renorm_axis, recolour)
    else:
        obj_renorm = obj
    if title:
        obj_renorm.SetTitle(title)
    canvas = ROOT.TCanvas(ROOT.TUUID().AsString(), "", 800, 800)
    canvas.SetTicks(1, 1)
    canvas.SetLeftMargin(0.13)
    canvas.SetBottomMargin(0.11)
    if logz:
        canvas.SetLogz(1)
    obj_renorm.Draw(draw_opt)
    if xlim is not None:
        obj_renorm.GetXaxis().SetRangeUser(*xlim)
    if ylim is not None:
        obj_renorm.GetYaxis().SetRangeUser(*ylim)
    obj_renorm.GetYaxis().SetTitleOffset(1.7)
    obj_renorm.GetXaxis().SetTitleOffset(1.2)
    if zlim:
        obj_renorm.SetMinimum(zlim[0])
        obj_renorm.SetMaximum(zlim[1])
    if other_things_to_draw:
        for x in other_things_to_draw:
            x.Draw()
    output_filename = os.path.abspath(output_filename)
    odir = os.path.dirname(output_filename)
    if not os.path.isdir(odir):
        os.makedirs(odir)
    canvas.SaveAs(output_filename)


def do_all_exclusive_plots_comparison(sources, plot_dir="plots_dy_vs_qcd",
                                      dy_filename=qgc.DY_FILENAME, qcd_filename=qgc.QCD_FILENAME,
                                      zpj_dirname="ZPlusJets_QG", dj_dirname="Dijet_QG",
                                      var_list=None, var_prepend="", pt_bins=None,
                                      subplot_type=None, subplot_title=None,
                                      do_flav_tagged=True, ofmt="pdf"):
    """Do 1D plots, comparing various sources. For each source plots DY & QCD samples. If zpj_dirname or dj_dirname blank, not plotted.

    Relies on QCD sample file being called uhh2.AnalysisModuleRunner.MC.MC_QCD_.root,
    and the DYJetsToLL one being called uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root
    """
    var_list = var_list or qgc.COMMON_VARS_WITH_FLAV
    pt_bins = pt_bins or qgc.PT_BINS

    if isinstance(sources, dict):
        sources = [sources]

    for ang in var_list:

        v = "%s%s_vs_pt" % (var_prepend, ang.var)
        for (start_val, end_val) in pt_bins:
            entries_normal, entries_flav = [], []

            # Get all plots
            for ind, source in enumerate(sources):
                lw = 2
                msize = 1.1

                if zpj_dirname:
                    h2d_dyj = grab_obj(os.path.join(source['root_dir'], dy_filename),
                                       "%s/%s" % (source.get('zpj_dirname', zpj_dirname), v))
                    dy_kwargs = dict(line_color=qgc.DY_COLOUR, line_width=lw, fill_color=qgc.DY_COLOUR,
                                     label=qgc.DY_ZpJ_LABEL + source.get('label', ''),
                                     marker_color=qgc.DY_COLOUR, marker_style=qgc.DY_MARKER+ind, marker_size=msize)
                    dy_kwargs.update(source.get('style', {}))
                    dy_kwargs.update(source.get('dy_style', {}))
                    entries_normal.append((get_projection_plot(h2d_dyj, start_val, end_val), dy_kwargs))

                if dj_dirname:
                    h2d_qcd = grab_obj(os.path.join(source['root_dir'], qcd_filename),
                                       "%s/%s" % (source.get('dj_dirname', dj_dirname), v))
                    qcd_kwargs = dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                                      label=qgc.QCD_Dijet_LABEL + source.get('label', ''),
                                      marker_color=qgc.QCD_COLOUR, marker_style=qgc.QCD_MARKER+ind, marker_size=msize)
                    qcd_kwargs.update(source.get('style', {}))
                    qcd_kwargs.update(source.get('qcd_style', {}))
                    entries_normal.append((get_projection_plot(h2d_qcd, start_val, end_val), qcd_kwargs))

                if not do_flav_tagged or "flavour" in v:
                    continue

                # Flav tagged plots
                if zpj_dirname:
                    h2d_dyj_q = grab_obj(os.path.join(source['root_dir'], dy_filename),
                                         "%s/q%s" % (source.get('zpj_dirname', zpj_dirname), v))
                    dy_kwargs_q = dict(line_color=qgc.DY_COLOUR, line_width=lw, fill_color=qgc.DY_COLOUR,
                                       label=qgc.DY_ZpJ_QFLAV_LABEL + source.get('label', ''),
                                       # label=qgc.QCD_Dijet_QFLAV_LABEL + source.get('label', ''),
                                       marker_color=qgc.DY_COLOUR, marker_style=qgc.DY_MARKER+ind, marker_size=msize)
                    dy_kwargs_q.update(source.get('style', {}))
                    dy_kwargs_q.update(source.get('dy_style', {}))
                    entries_flav.append((get_projection_plot(h2d_dyj_q, start_val, end_val), dy_kwargs_q))

                if dj_dirname:
                    h2d_qcd_g = grab_obj(os.path.join(source['root_dir'], qcd_filename),
                                         "%s/g%s" % (source.get('dj_dirname', dj_dirname), v))
                    qcd_kwargs_g = dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                                        label=qgc.QCD_Dijet_GFLAV_LABEL + source.get('label', ''),
                                        marker_color=qgc.QCD_COLOUR, marker_style=qgc.QCD_MARKER+ind, marker_size=msize)
                    qcd_kwargs_g.update(source.get('style', {}))
                    qcd_kwargs_g.update(source.get('qcd_style', {}))
                    entries_flav.append((get_projection_plot(h2d_qcd_g, start_val, end_val), qcd_kwargs_g))

            rebin = 2
            if "multiplicity" in v:
                rebin = 2
            elif "flavour" in v or "thrust" in v or 'pTD' in v:
                rebin = 1

            xlim = None
            if "width" in v or "thrust" in v or "pTD" in v:
                xlim = (0, 0.5)
            elif "multiplicity" in v.lower() and "ak4" in sources[0]['root_dir'].lower():
                xlim = (0, 100)

            ylim = None
            if "flavour" in v:
                ylim = (0, 1)
            elif "LHA" in v:
                ylim = (0, 5)

            do_comparison_plot(entries_normal, "%s/ptBinned/%s_pt%dto%d.%s" % (plot_dir, v, start_val, end_val, ofmt),
                               rebin=rebin, title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val),
                               xtitle=ang.name + " (" + ang.lambda_str + ")",
                               xlim=xlim, ylim=ylim, subplot_type=subplot_type, subplot_title=subplot_title, subplot_limits=(0, 2))

            if do_flav_tagged and "flavour" not in v:
                do_comparison_plot(entries_flav, "%s/ptBinned/%s_pt%dto%d_flavMatched.%s" % (plot_dir, v, start_val, end_val, ofmt),
                                   rebin=rebin, title="%d < p_{T}^{jet} < %d GeV" % (start_val, end_val),
                                   xtitle=ang.name + " (" + ang.lambda_str + ")",
                                   xlim=xlim, subplot_type=subplot_type, subplot_title=subplot_title, subplot_limits=(0, 2))



def transpose_2d_hist(hist):
    hnew = ROOT.TH2D(hist.GetName()+"transpose",
                     ";".join([hist.GetTitle(), hist.GetYaxis().GetTitle(), hist.GetXaxis().GetTitle()]),
                     hist.GetNbinsY(), hist.GetYaxis().GetXmin(), hist.GetYaxis().GetXmax(),
                     hist.GetNbinsX(), hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax()
                     )
    for i in range(1, hist.GetNbinsX()+1):
        for j in range(1, hist.GetNbinsY()+1):
            hnew.SetBinContent(j, i, hist.GetBinContent(i, j))
            hnew.SetBinError(j, i, hist.GetBinError(i, j))
    return hnew


def rescale_plot_labels(container, factor):
    # What a pile of poop, why does ROOT scale all these sizes?
    container.GetXaxis().SetLabelSize(container.GetXaxis().GetLabelSize()/factor)
    container.GetXaxis().SetTitleSize(container.GetXaxis().GetTitleSize()/factor)
    container.GetXaxis().SetTitleOffset(container.GetXaxis().GetTitleOffset()*factor)  # doesn't seem to work?
    container.GetXaxis().SetTickLength(container.GetXaxis().GetTickLength()/factor)

    container.GetYaxis().SetLabelSize(container.GetYaxis().GetLabelSize()/factor)
    container.GetYaxis().SetTitleSize(container.GetYaxis().GetTitleSize()/factor)
    container.GetYaxis().SetTitleOffset(container.GetYaxis().GetTitleOffset()*factor)


def do_box_plot(entries, output_filename, xlim=None, ylim=None, transpose=False):
    """Create box n whickers plot from 2D hists

    Parameters
    ----------
    entries : list[(TH2, dict)]
        Description
    output_filename : str
        Output filename
    """
    hists2d = []  # keep references alive
    canvas = ROOT.TCanvas(cu.get_unique_str(), "", 800, 800)
    canvas.SetTicks(1, 1)
    right_margin = 0.03
    top_margin = 0.1
    subplot_pad_height = 0.32
    subplot_pad_fudge = 0.01  # to get non-overlapping subplot axis
    subplot_type = 'ratio'
    # subplot_type = None
    if subplot_type:
        main_pad = ROOT.TPad("main_pad", "", 0, subplot_pad_height+subplot_pad_fudge, 1, 1)
        ROOT.SetOwnership(main_pad, False)
        main_pad.SetTicks(1, 1)
        main_pad.SetBottomMargin(2*subplot_pad_fudge)
        main_pad.SetTopMargin(top_margin / (1-subplot_pad_height))
        main_pad.SetRightMargin(right_margin / (1-subplot_pad_height))
        canvas.cd()
        main_pad.Draw()
        subplot_pad = ROOT.TPad("subplot_pad", "", 0, 0, 1, subplot_pad_height-subplot_pad_fudge)
        ROOT.SetOwnership(subplot_pad, False)
        subplot_pad.SetTicks(1, 1)
        subplot_pad.SetFillColor(0)
        subplot_pad.SetFillStyle(0)
        subplot_pad.SetTopMargin(4*subplot_pad_fudge)
        # subplot_pad.SetRightMargin(right_margin)
        subplot_pad.SetBottomMargin(0.35)
        canvas.cd()
        subplot_pad.Draw()
    else:
        main_pad = ROOT.TPad("main_pad", "", 0, 0, 1, 1)
        ROOT.SetOwnership(main_pad, False)
        main_pad.SetRightMargin(right_margin)
        main_pad.SetTopMargin(top_margin)
        main_pad.SetTicks(1, 1)
        main_pad.Draw()

    leg = ROOT.TLegend(0.5, 0.7, 0.92, 0.88)

    main_pad.cd()

    median_hists = []
    lower_hists = []
    upper_hists = []

    quantiles = array('d', [0.25, 0.5, 0.75])

    for ind, ent in enumerate(entries[:]):
        quantile_values = array('d', [0., 0., 0.])

        # rebin_hist = ent[0].RebinX(20, cu.get_unique_str())
        # rebin_hist = ent[0].RebinX(20)
        rebin_hist = ent[0]
        if transpose:
            rebin_hist = transpose_2d_hist(rebin_hist)

        median_hist = ROOT.TH1D("mean_%d" % ind,
                               "",
                               rebin_hist.GetNbinsX(),
                               rebin_hist.GetXaxis().GetBinLowEdge(1),
                               rebin_hist.GetXaxis().GetBinLowEdge(rebin_hist.GetNbinsX()+1))
        lower_hist = ROOT.TH1D("lower_%d" % ind,
                               "",
                               rebin_hist.GetNbinsX(),
                               rebin_hist.GetXaxis().GetBinLowEdge(1),
                               rebin_hist.GetXaxis().GetBinLowEdge(rebin_hist.GetNbinsX()+1))
        upper_hist = ROOT.TH1D("upper_%d" % ind,
                               "",
                               rebin_hist.GetNbinsX(),
                               rebin_hist.GetXaxis().GetBinLowEdge(1),
                               rebin_hist.GetXaxis().GetBinLowEdge(rebin_hist.GetNbinsX()+1))

        for i in range(1, rebin_hist.GetNbinsX()):
            projection = rebin_hist.ProjectionY("_py%s"%cu.get_unique_str(), i, i+1)
            projection.GetQuantiles(len(quantiles), quantile_values, quantiles)
            median_hist.SetBinContent(i, quantile_values[1])
            lower_hist.SetBinContent(i, quantile_values[0])
            upper_hist.SetBinContent(i, quantile_values[2])

        median_hist.SetLineColor(ent[1]['line_color'])
        median_hist.SetFillColor(ent[1]['line_color'])
        median_hist.SetFillStyle(0)
        median_hists.append(median_hist)

        lower_hist.SetLineColor(ent[1]['line_color'])
        lower_hist.SetLineStyle(2)
        lower_hist.SetFillColor(ent[1]['line_color'])
        lower_hist.SetFillStyle(0)
        lower_hists.append(lower_hist)

        upper_hist.SetLineColor(ent[1]['line_color'])
        upper_hist.SetLineStyle(3)
        upper_hist.SetFillColor(ent[1]['line_color'])
        upper_hist.SetFillStyle(0)
        upper_hists.append(upper_hist)


        rebin_hist.SetBarWidth(0.1)
        offset = 0.011
        half = (int) (len(entries) / 2)
        factor = -1 if transpose else 1
        rebin_hist.SetBarOffset(factor * (half * offset - (ind * offset)))

        rebin_hist.SetLineColor(ent[1]['line_color'])

        rebin_hist.SetFillColor(ent[1]['line_color'])
        rebin_hist.SetFillStyle(0)

        if xlim and len(xlim) == 2:
            rebin_hist.SetAxisRange(xlim[0], xlim[1], 'X')
        if ylim and len(ylim) == 2:
            rebin_hist.SetAxisRange(ylim[0], ylim[1], 'Y')

        leg.AddEntry(rebin_hist, ent[1]['label'], "LP")
        ax = "X" if transpose else "Y"
        draw_opt = "CANDLE"+ax+"(00011311)"
        if ind > 0:
            draw_opt += "SAME"
        rebin_hist.Draw(draw_opt)
        hists2d.append(rebin_hist)

    canvas.cd()
    leg.Draw()
    cms_latex = ROOT.TLatex()
    cms_latex.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    cms_latex.SetTextFont(42)
    cms_latex.SetTextSize(0.035)
    latex_height = 0.92
    # cms_latex.DrawLatex(0.14, latex_height, "#font[62]{CMS}#font[52]{ Preliminary}")
    cms_latex.DrawLatexNDC(0.14, latex_height, "#font[62]{CMS}#font[52]{ Preliminary}")
    # cms_latex.DrawLatex(0.14, latex_height, "#font[62]{CMS}#font[52]{ Preliminary Simulation}")
    # cms_latex.DrawLatex(0.14, latex_height, "#font[62]{CMS}")
    cms_latex.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignBottom)
    cms_latex.DrawLatexNDC(0.97, latex_height, " 35.9 fb^{-1} (13 TeV)")

    subplot_title = None
    subplot_pad.cd()
    ratio_hst = ROOT.THStack("hst", ";"+hists2d[0].GetXaxis().GetTitle()+";#splitline{Ratio of median}{(MC / Data)}")
    ratio_hst = ROOT.THStack("hst", ";"+hists2d[0].GetXaxis().GetTitle()+";#splitline{Ratio of quantiles}{(MC / Data)}")
    if subplot_type:
        # # Get rid of main plot x axis labels
        # modifier.GetHistogram().GetXaxis().SetLabelSize(0)
        # modifier.GetXaxis().SetLabelSize(0)
        # modifier.GetHistogram().GetXaxis().SetTitleSize(0)
        # modifier.GetXaxis().SetTitleSize(0)


        if len(median_hists) > 1:
            for qh in median_hists[1:]:
                qh.Divide(median_hists[0])
                ratio_hst.Add(qh)

            for qh in lower_hists[1:]:
                qh.Divide(lower_hists[0])
                ratio_hst.Add(qh)

            for qh in upper_hists[1:]:
                qh.Divide(upper_hists[0])
                ratio_hst.Add(qh)
            ratio_hst.Draw("NOSTACK HIST")

        # subplot_container.Draw(draw_opts)

        # if subplot_title == None:
        #     if (subplot_type == "ratio"):
        #         subplot_title = "#splitline{Ratio vs}{%s}" % (subplot.label)
        #     elif (subplot_type == "diff"):
        #         subplot_title = "#splitline{Difference}{vs %s}" % (subplot.label)
        #     elif (subplot_type == "ddelta"):
        #         subplot_title = "d#Delta/d#lambda"

        # ratio_hst.SetTitle(";%s;%s" % (xtitle, subplot_title))


        # if xlim:
            ratio_hst.GetXaxis().SetRangeUser(0, 320)

            ratio_hst.SetMinimum(0.9)  # use this, not SetRangeUser()
            ratio_hst.SetMaximum(1.1)  # use this, not SetRangeUser()

            xax = ratio_hst.GetXaxis()
            subplot_line = ROOT.TLine(0, 1., 320, 1.)
            subplot_line.SetLineStyle(4)
            subplot_line.SetLineWidth(1)
            subplot_line.SetLineColor(ROOT.kBlack)
            subplot_line.Draw()

            rescale_plot_labels(ratio_hst, subplot_pad_height)
            ratio_hst.GetXaxis().SetTitleOffset(ratio_hst.GetXaxis().GetTitleOffset()*3)
            ratio_hst.GetYaxis().SetNdivisions(505)

    canvas.SaveAs(output_filename)


def extract_sample_name(label):
    parts = label.split("[")
    new_label = parts[-1]
    new_label = new_label.replace("]", "")
    return new_label


def do_box_plot_mpl(entries, bins, output_filename, var_label="", xlim=None, ylim=None, region_title=""):
    """Create box n whickers plot from 1D hists using matplotlib

    Parameters
    ----------
    entries : list[list[(TH1, dict)]]
        Description
    output_filename : str
        Output filename
    """
    # Get the statistics from our 1D hists
    # entries is a list, where each element corresponds to a pt bin
    # each element is itself a list, with different entries for the different datasets being compared

    quantiles = array('d', [0.25, 0.5, 0.75])

    lower_quants = []
    median_quants = []
    upper_quants = []

    # only do the bins we need, otherwise axis limits go awry
    if xlim:
        bins = [b for b in bins if b[1] <= xlim[1]]

    for ind, bin_entry in enumerate(entries):
        if ind == len(bins):
            break
        this_lower = []
        this_median = []
        this_upper = []
        for entry in bin_entry:
            quantile_values = array('d', [0., 0., 0.])
            hist = entry[0]
            if cu.get_hist_mean_rel_error(hist) < 0.4:
                hist.GetQuantiles(len(quantiles), quantile_values, quantiles)
                this_lower.append(quantile_values[0])
                this_median.append(quantile_values[1])
                this_upper.append(quantile_values[2])
            else:
                this_lower.append(np.nan)
                this_median.append(np.nan)
                this_upper.append(np.nan)
        lower_quants.append(this_lower)
        median_quants.append(this_median)
        upper_quants.append(this_upper)

    # convert to numpy array to maniuplate easier
    lower_quants = np.array(lower_quants)
    median_quants = np.array(median_quants)
    upper_quants = np.array(upper_quants)

    # Do all plot - start by setting up figure
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(16, 8), sharex=True)
    plt.subplots_adjust(left=0.1, right=0.78, top=0.9, hspace=0.1)

    n_datasets = len(entries[0])
    bin_widths = [(x[1]-x[0])/(n_datasets) for x in bins]

    # Convert ROOT colors to hex for matplotlib
    cols = []
    for ind, entry in enumerate(entries[0]):
        cols.append(ROOT.gROOT.GetColor(entry[1]['line_color']).AsHexString())

    # main plot
    plot_objs = []
    # loop over datasets
    for ind, entry in enumerate(entries[0]):  # assumes all are in 0th element!
        bxpstats = []
        for med, low, upper in zip(median_quants[:,ind], lower_quants[:,ind], upper_quants[:,ind]):
            bxpstats.append({
                "med": med,
                "q1": low,
                "q3": upper,
                "whislo": low,
                "whishi": upper,
            })
        boxprops = {
            "color": cols[ind],
            "linewidth": 2,
        }
        x_pos = [(x[0] + (0.5+ind)*binw) for x, binw in zip(bins, bin_widths)]
        # DO NOT USE axes[0,0], get indexerror
        # Make boxplot
        things = axes[0].bxp(bxpstats, positions=x_pos, widths=bin_widths,
                            vert=True, showfliers=False, medianprops=boxprops,
                            boxprops=boxprops, showcaps=False)

        # Explicitly set the label and store for the legend
        things['boxes'][0].set_label(extract_sample_name(entry[1]['label']).split(',')[0])
        plot_objs.append(things['boxes'][0])

        # loc = plticker.LogLocator(base=10.0, subs=(1.0, 2.0, 5.0))
        # axes[1].xaxis.set_major_locator(loc)
        # axes[0].xaxis.set_minor_locator(loc)
        # axes[0].set_ylim(bottom=0)

    leg_loc = (1.05, 0.1)
    labels = [e[1]['label'] for e in entries[0]]
    axes[0].legend(handles=plot_objs, loc=leg_loc, fancybox=False, edgecolor='white')

    plt.xscale('log')

    # convert from ROOTism to latex
    new_var_label = var_label.replace("#", "\\")
    axes[0].set_ylabel('%s' % new_var_label)
    axes[-1].set_xlabel('$p_{T}^{\mathrm{jet}}$ [GeV]')

    axes[0].set_xlim(*xlim)
    if ylim:
        axes[0].set_ylim(*ylim)
    if axes[0].get_ylim()[0] > 0:
        axes[0].set_ylim(bottom=0)

    axes[-1].set_xlim(*xlim)

    # subplot - ratio of medians
    ref_data = None
    x_pos = [0.5*sum(x) for x in bins]
    bin_widths = [0.5*(x[1]-x[0]) for x in bins]
    markers = ['x', 'o', 'd', '^']
    for ind, entry in enumerate(entries[0]):
        if ind == 0:
            ref_data = median_quants[:, ind]
        else:
            this_data = median_quants[:, ind] / ref_data
            axes[1].errorbar(x_pos, this_data, xerr=bin_widths,
                             marker=markers[ind],
                             label=extract_sample_name(entry[1]['label']),
                             color=cols[ind],
                             linewidth=0, elinewidth=2)
    axes[1].axhline(1.0, color='gray', linestyle='dashed')
    axes[1].set_ylabel("Ratio of medians\n(MC / Data)")
    ylim = axes[1].get_ylim()
    lim = 2
    if ylim[1] > lim:
        axes[1].set_ylim(top=lim)

    # subplot - ratio of lower/upper quantiles
    # ref_data_lower = None
    # ref_data_upper = None
    # for ind, entry in enumerate(entries[0]):
    #     if ind == 0:
    #         ref_data_lower = lower_quants[:, ind]
    #         ref_data_upper = upper_quants[:, ind]
    #     else:
    #         this_data_lower = lower_quants[:, ind] / ref_data_lower
    #         this_data_upper = upper_quants[:, ind] / ref_data_upper
    #         axes[2].errorbar(x_pos, this_data_lower, xerr=bin_widths,
    #                          marker=markers[ind],
    #                          label=extract_sample_name(entry[1]['label']) + " [25%]",
    #                          color=cols[ind],
    #                          linewidth=0, elinewidth=2)

    #         axes[2].errorbar(x_pos, this_data_upper, xerr=bin_widths,
    #                          marker=markers[ind],
    #                          label=extract_sample_name(entry[1]['label']) + " [75%]",
    #                          color=cols[ind],
    #                          linewidth=0, elinewidth=2, linestyle='dashed', markerfacecolor='white')
    # ylim = axes[2].get_ylim()
    # lim = 2
    # if ylim[1] > lim:
    #     axes[2].set_ylim(top=lim)
    # axes[2].axhline(1.0, color='gray', linestyle='dashed')
    # axes[2].set_ylabel("Ratio of upper\n& lower quantiles\n(MC / Data)")
    # axes[2].legend(loc=leg_loc, fancybox=False)

    # Subplot - ratio of width (75%-25%)
    ref_data = None
    for ind, entry in enumerate(entries[0]):
        if ind == 0:
            ref_data = upper_quants[:, ind]-lower_quants[:, ind]
        else:
            this_data = (upper_quants[:, ind]-lower_quants[:, ind]) / ref_data
            axes[2].errorbar(x_pos, this_data, xerr=bin_widths,
                             marker=markers[ind],
                             label=extract_sample_name(entry[1]['label']),
                             color=cols[ind],
                             linewidth=0, elinewidth=2)

    ylim = axes[2].get_ylim()
    lim = 2
    if ylim[1] > lim:
        axes[2].set_ylim(top=lim)
    axes[2].axhline(1.0, color='gray', linestyle='dashed')
    axes[2].set_ylabel("Ratio of central\n50% quantile width\n(MC / Data)")
    # axes[2].legend(loc=leg_loc, fancybox=False)


    axes[0].set_title("25%, 50%, 75% quantiles for {} region".format(region_title), pad=40)

    # Draw bin delineators
    for ax in axes.flat:
        ax.tick_params(axis='x', which='major', direction='in', top=True, bottom=True, length=10)
        ax.tick_params(axis='x', which='minor', direction='in', top=True, bottom=True, length=7)
        ax.tick_params(axis='y', which='major', direction='in', left=True, right=True, length=10)
        ax.tick_params(axis='y', which='minor', direction='in', left=True, right=True, length=7)

        this_ylim = ax.get_ylim()
        for ind, this_bin in enumerate(bins):
            # if xlim and this_bin[1] < xlim[0]:
            #     continue
            # if xlim and this_bin[0] > xlim[1]:
            #     continue
            # if ind % 2 == 0:
            #     continue
            # ax.add_patch(patches.Rectangle(
            #                         (this_bin[1], this_ylim[0]),
            #                         width=this_bin[1]-this_bin[0],
            #                         height=this_ylim[1]-this_ylim[0],
            #                         color='lightgray', zorder=-9999999
            #                         ))
            ax.axvline(this_bin[0], color='lightgray', linestyle='dashed')
    plt.savefig(output_filename)


def migration_plot_components(h2d_renorm_x, h2d_renorm_y, xlabel):
    stability_values, stability_errs = [], []
    xfer_down_values, xfer_down_errs = [], []
    xfer_down2_values, xfer_down2_errs = [], []
    xfer_up_values, xfer_up_errs = [], []
    xfer_up2_values, xfer_up2_errs = [], []

    binning = [h2d_renorm_x.GetXaxis().GetBinLowEdge(bin_ind) for bin_ind in range(1, h2d_renorm_x.GetNbinsX()+2)]

    n_bins = len(binning) - 1
    label = ";%s;Fraction" % xlabel
    hist_purity = ROOT.TH1F("purity"+cu.get_unique_str(), xlabel, n_bins, array('d', binning))
    hist_stability = ROOT.TH1F("stability"+cu.get_unique_str(), xlabel, n_bins, array('d', binning))
    hist_xfer_down = ROOT.TH1F("xfer_down"+cu.get_unique_str(), xlabel, n_bins, array('d', binning))
    hist_xfer_down2 = ROOT.TH1F("xfer_down2"+cu.get_unique_str(), xlabel, n_bins, array('d', binning))
    hist_xfer_down3 = ROOT.TH1F("xfer_down3"+cu.get_unique_str(), xlabel, n_bins, array('d', binning))
    hist_xfer_up = ROOT.TH1F("xfer_up"+cu.get_unique_str(), xlabel, n_bins, array('d', binning))
    hist_xfer_up2 = ROOT.TH1F("xfer_up2"+cu.get_unique_str(), xlabel, n_bins, array('d', binning))
    hist_xfer_up3 = ROOT.TH1F("xfer_up3"+cu.get_unique_str(), xlabel, n_bins, array('d', binning))

    for bin_ind in range(1, h2d_renorm_x.GetNbinsX()+1):
        hist_purity.SetBinContent(bin_ind, h2d_renorm_y.GetBinContent(bin_ind, bin_ind))
        hist_purity.SetBinError(bin_ind, h2d_renorm_y.GetBinError(bin_ind, bin_ind))

        hist_stability.SetBinContent(bin_ind, h2d_renorm_x.GetBinContent(bin_ind, bin_ind))
        hist_stability.SetBinError(bin_ind, h2d_renorm_x.GetBinError(bin_ind, bin_ind))

        hist_xfer_down.SetBinContent(bin_ind, h2d_renorm_x.GetBinContent(bin_ind, bin_ind-1))
        hist_xfer_down.SetBinError(bin_ind, h2d_renorm_x.GetBinError(bin_ind, bin_ind-1))
        hist_xfer_down2.SetBinContent(bin_ind, h2d_renorm_x.GetBinContent(bin_ind, bin_ind-2))
        hist_xfer_down2.SetBinError(bin_ind, h2d_renorm_x.GetBinError(bin_ind, bin_ind-2))
        hist_xfer_down3.SetBinContent(bin_ind, h2d_renorm_x.GetBinContent(bin_ind, bin_ind-3))
        hist_xfer_down3.SetBinError(bin_ind, h2d_renorm_x.GetBinError(bin_ind, bin_ind-3))

        hist_xfer_up.SetBinContent(bin_ind, h2d_renorm_x.GetBinContent(bin_ind, bin_ind+1))
        hist_xfer_up.SetBinError(bin_ind, h2d_renorm_x.GetBinError(bin_ind, bin_ind+1))
        hist_xfer_up2.SetBinContent(bin_ind, h2d_renorm_x.GetBinContent(bin_ind, bin_ind+2))
        hist_xfer_up2.SetBinError(bin_ind, h2d_renorm_x.GetBinError(bin_ind, bin_ind+2))
        hist_xfer_up3.SetBinContent(bin_ind, h2d_renorm_x.GetBinContent(bin_ind, bin_ind+3))
        hist_xfer_up3.SetBinError(bin_ind, h2d_renorm_x.GetBinError(bin_ind, bin_ind+3))

    # col_purity = ROOT.kViolet-3
    # col_purity = ROOT.kOrange-3
    col_purity = ROOT.kGreen+1
    col_stability = ROOT.kBlack
    col_xfer_down = ROOT.kRed
    col_xfer_down2 = ROOT.kMagenta
    col_xfer_down3 = ROOT.kGreen
    col_xfer_up = ROOT.kBlue-4
    col_xfer_up2 = ROOT.kAzure+8
    col_xfer_up3 = ROOT.kOrange
    contributions = [
        Contribution(hist_purity, label="Purity (gen right)", line_color=col_purity, marker_color=col_purity),
        Contribution(hist_stability, label="Stability (reco right)", line_color=col_stability, marker_color=col_stability),
        Contribution(hist_xfer_down, label="-1 reco bin", line_color=col_xfer_down, marker_color=col_xfer_down),
        Contribution(hist_xfer_down2, label="-2 reco bin", line_color=col_xfer_down2, marker_color=col_xfer_down2),
        # Contribution(hist_xfer_down3, label="3 lower reco bin", line_color=col_xfer_down3, marker_color=col_xfer_down3),
        Contribution(hist_xfer_up, label="+1 reco bin", line_color=col_xfer_up, marker_color=col_xfer_up),
        Contribution(hist_xfer_up2, label="+2 reco bin", line_color=col_xfer_up2, marker_color=col_xfer_up2),
        # Contribution(hist_xfer_up3, label="3 higher reco bin", line_color=col_xfer_up3, marker_color=col_xfer_up3),
    ]
    return contributions


def make_migration_summary_plot(h2d_renorm_x, h2d_renorm_y, xlabel, output_filename, log_var=False):
    contributions = migration_plot_components(h2d_renorm_x, h2d_renorm_y, xlabel)
    for c in contributions:
        c.obj.SetLineWidth(2)
    binning = [h2d_renorm_x.GetXaxis().GetBinLowEdge(bin_ind) for bin_ind in range(1, h2d_renorm_x.GetNbinsX()+2)]
    xlim = [binning[0], binning[-1]]
    plot = Plot(contributions, what='hist', xlim=xlim, ylim=[1e-3, 2], xtitle=xlabel, has_data=False)
    y1 = 0.15
    plot.legend.SetX1(0.5)
    plot.legend.SetY1(y1)
    plot.legend.SetY2(y1+0.25)
    plot.plot("NOSTACK HISTE")
    plot.legend.SetFillStyle(1001)
    plot.legend.SetFillColorAlpha(ROOT.kWhite, 0.75)
    if log_var:
        plot.set_logx()
    plot.set_logy()
    plot.main_pad.cd()
    lines = []
    for val in [1, 0.5, 1e-1, 1e-2, 1e-3]:
        line = ROOT.TLine(xlim[0], val, xlim[1], val)
        line.SetLineStyle(2)
        line.SetLineColor(ROOT.kGray+2)
        lines.append(line)
        line.Draw("same")
    plot.save(output_filename)
