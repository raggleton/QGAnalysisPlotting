"""Common functions for doing simple QG plots"""


import ROOT
from MyStyle import My_Style
My_Style.cd()
import bisect
import os

# My stuff
from comparator import Contribution, Plot, grab_obj
import common_utils as cu
import qg_common as qgc

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
        Remove contributions that have a 
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
        p.legend.SetX2(0.97)
        if len(entries) > 4:
            p.legend.SetY1(0.6)
        else:
            p.legend.SetY1(0.7)
        p.legend.SetY2(0.9)
    return p


def do_comparison_plot(entries, output_filename, rebin=1, **plot_kwargs):
    """Plot several different objects on a single plot

    entries : list of 2-tuples, with (object, dict), where the dict is a set of kwargs passed to the Contribution object
    plot_kwargs : any other kwargs to be passed to the Plot object ctor
    """
    p = make_comparison_plot_ingredients(entries, rebin=rebin, mean_rel_error=0.5, **plot_kwargs)
    draw_opt = "NOSTACK HISTE"
    p.plot(draw_opt)
    # p.container.SetMaximum(min(5, p.container.GetMaximum()))
    p.save(output_filename)


def get_projection_plot(h2d, start_val, end_val, cut_axis='y'):
    cut_axis = cut_axis.lower()
    if cut_axis == "y":
        axis = h2d.GetYaxis()
    elif cut_axis == "x":
        axis = h2d.GetXaxis()
    else:
        raise RuntimeError("cut_axis must be x or y")
    bin_edges = [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins()+1)]
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
