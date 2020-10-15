"""Common functions for doing simple QG plots"""


import ROOT
from MyStyle import My_Style
My_Style.cd()
import bisect
import os
from array import array
from math import sqrt, log10

# My stuff
from comparator import Contribution, Plot, grab_obj
import common_utils as cu
import qg_common as qgc
import metric_calculators as metrics

import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import matplotlib.patches as patches


mpl.rcParams['font.size'] = 16
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


def make_comparison_plot_ingredients(entries, rebin=1, normalise_hist=True, mean_rel_error=1.0, data_first=False, **plot_kwargs):
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
        Rebin contributions that have a mean relative error more than this value,
        where "mean relative error" = mean across all (error/bin contents) in a hist
    data_first : bool, optional
        If True, then the first entry is data, and we do special things like the ratio plot
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
    entries = [ent for ent in entries if ent]

    def _get_this_xlim():
        data_hist = entries[0][0]
        return (calc_auto_xlim([data_hist]) if (data_first and data_hist.GetEntries()>0)
                                            else calc_auto_xlim([e[0] for e in entries]))

    # first figure out if there is a plot where the mean relative error
    #  is greater than a certain amount
    orig_rebin = rebin
    if mean_rel_error > 0:
        big_mean_rel_err = any([cu.get_hist_mean_rel_error(ent[0]) > mean_rel_error
                                # and ent[0].Integral() > 0
                                for ent in entries])
        nbins = entries[0][0].GetNbinsX()

        # Get actual visible number of bins - that's what will look good/bad
        # only care about data if it's there
        xlim = _get_this_xlim()
        bin_ind_min = entries[0][0].FindFixBin(xlim[0])
        bin_ind_max = entries[0][0].FindFixBin(xlim[1])
        nbins_vis = bin_ind_max - bin_ind_min

        # keep looking for the next rebin factor until the mean_rel_error is small enough,
        # or we have < 50 visible bins, but make sure we have at least a few bins
        min_number_bins = 4
        max_number_bins = 50
        while (big_mean_rel_err or (nbins_vis / rebin) > max_number_bins) and ((nbins_vis / rebin) > min_number_bins):
            # If this is the case, rebin the plot. Start by making bins bigger,
            # but check if divisor
            rebin += 1

            if entries[0][0].GetNbinsX() % rebin != 0:
                continue

            # now rebin with this value and re-check mean rel error
            hist_copies = [ent[0].Clone(cu.get_unique_str()) for ent in entries]
            for h in hist_copies:
                h.Rebin(rebin)
            big_mean_rel_err = any([cu.get_hist_mean_rel_error(h) > mean_rel_error
                                    for h in hist_copies])
            # print("New rebin:", rebin)

        # the final rebin value may not be an exact divisor, in which case we reduce it until it is
        while nbins % rebin != 0:
            rebin -= 1

        # print("final rebin:", rebin)

    conts = [Contribution(ent[0], normalise_hist=normalise_hist, rebin_hist=rebin, **ent[1])
             for ent in entries]

    # update auto xlim if necessary
    if plot_kwargs.get("xlim", None) == "auto":
        plot_kwargs['xlim'] = _get_this_xlim()

    do_legend = len(conts) > 1
    if len(conts) == 0:
        raise RuntimeError("0 contributions for this plot")

    if data_first:
        plot_kwargs['has_data'] = True

        # Modify data Contribution to have no  lines
        # We'll redraw it at the end with the proper ones ontop
        conts[0].line_width = 0
        conts[0].marker_size = 0.01
        conts[0].update_obj_styling()

        # For subplot to ensure only MC errors drawn, not MC+data
        data_obj = entries[0][0]
        data_no_errors = data_obj.Clone()
        cu.remove_th1_errors(data_no_errors)
        for i, cont in enumerate(conts[1:], 1):
            if cont.subplot == data_obj:
                cont.subplot = data_no_errors

    # Automatically correct plot_kwargs if we only have 1 contribution (don't want subplot)
    do_subplot = any(c.subplot for c in conts) or "subplot" in plot_kwargs
    if (len(conts) == 1 or not do_subplot) and "subplot_type" in plot_kwargs:
        plot_kwargs['subplot_type'] = None
        print("Not doing subplot")

    # If 'subplot' in plot_kwargs, then Plot expects it to be a Contribution object
    # But we only make those here
    # So here we figure out which contribution matches the subplot, if it exists
    if "subplot" in plot_kwargs and plot_kwargs['subplot'] is not None:
        # FIXME check type, could match obj or Contribution
        subplot_cont = [c for c in conts if c.obj == plot_kwargs['subplot']]
        if len(subplot_cont) > 0:
            plot_kwargs['subplot'] = subplot_cont[0]
        else:
            print("Couldn't identify subplot object, no subplot")
            plot_kwargs['subplot_type'] = None
            del plot_kwargs['subplot']

    ytitle = "p.d.f." if normalise_hist else "N"
    ytitle = "#DeltaN/N" if normalise_hist else "N"
    p = Plot(conts, what="hist",
             ytitle=ytitle,
             legend=do_legend,
             **plot_kwargs)
    if do_legend:
        # ensure legend big enough, but not too big, depending on how long entries are
        max_leg_str = max([len(c.label) for c in conts])
        if max_leg_str < 12:
            p.legend.SetX1(0.65)
        else:
            p.legend.SetX1(0.55)

        p.legend.SetX2(0.99)
        if len(entries) > 4:
            p.legend.SetY1(0.6)
        elif len(entries) > 6:
            p.legend.SetY1(0.5)
        else:
            p.legend.SetY1(0.72)
        p.legend.SetY2(0.88)

    if mean_rel_error > 0:
        rebin = orig_rebin
    return p


def do_comparison_plot(entries, output_filename, rebin=1, draw_opt="NOSTACK HIST E", data_first=False, mean_rel_error=0.4, logx=False, logy=False, **plot_kwargs):
    """Plot several different objects on a single plot

    entries : list of 2-tuples, with (object, dict), where the dict is a set of kwargs passed to the Contribution object
    plot_kwargs : any other kwargs to be passed to the Plot object ctor
    """
    try:

        plot = make_comparison_plot_ingredients(entries,
                                                rebin=rebin,
                                                mean_rel_error=mean_rel_error,
                                                data_first=data_first,
                                                **plot_kwargs)
        plot.subplot_maximum_ceil = 5

        if data_first:
            # we'll do the filling of legend ourselves
            plot.do_legend = False

        if logy:
            # figure out some limits based on logarithmic range
            ymin = min([e[0].GetMinimum(1E-20) for e in entries if e[0].Integral() > 0])
            ymax = min([e[0].GetMaximum() for e in entries if e[0].Integral() > 0])
            pow_min = log10(ymin)
            pow_max = log10(ymax)
            pow_range = pow_max - pow_min
            # print(ymin, ymax)
            # print(pow_min, pow_max, pow_range)
            pow_extra = pow_range / 2
            plot.y_padding_max_log = 5 * 10**pow_extra
            # print(plot.y_padding_max_log)

        plot.plot(draw_opt)
        plot.get_modifier().GetYaxis().SetTitleOffset(plot.get_modifier().GetYaxis().GetTitleOffset()*1.1)
        
        if logx:
            plot.set_logx(do_more_labels=False)
        if logy:
            plot.set_logy(do_more_labels=False)

        # Special case if data first object
        if data_first:
            # Redraw data on top
            data_hist = plot.contributions[0].obj.Clone()
            data_hist.SetLineWidth(entries[0][1].get('line_width', 1))
            data_hist.SetMarkerSize(entries[0][1].get('marker_size', 1))
            plot.main_pad.cd()
            data_draw_opt = "E1 X0 SAME"
            data_hist.Draw(data_draw_opt)

            # Create dummy graphs with the same styling to put into the legend
            # Using graphs we can get the correct endings on the TLegend entries (!)
            # Yes a massive faff for something so trivial
            dummy_gr = ROOT.TGraphErrors(1, array('d', [1]), array('d', [1]), array('d', [1]), array('d', [1]))
            # dummy_hist = ROOT.TGraphErrors(1, array('d', [1]), array('d', [1]), array('d', [1]), array('d', [1]))
            # Add them to the legend and draw it
            dummy_conts = []  # stop premature deletion
            for i, entry in enumerate(entries):
                leg_draw_opt = "LE"
                # check if line_width > 0?
                if i == 0:
                    if "X0" in data_draw_opt:
                        leg_draw_opt = "E"
                    if data_hist.GetMarkerSize() > 0:
                        leg_draw_opt += "P"
                else:
                    if entry[1].get('marker_size', 0) > 0:
                        leg_draw_opt += "P"
                # print(entry[1])

                cont = Contribution(dummy_gr.Clone(), leg_draw_opt=leg_draw_opt, **entry[1])
                dummy_conts.append(cont)

                # Split text by newline \n
                # Add an entry for each line
                for i, substr in enumerate(cont.label.split("\n")):
                    if i == 0:
                        plot.legend.AddEntry(cont.obj, substr, cont.leg_draw_opt)
                    else:
                        plot.legend.AddEntry(0, substr, "")

            plot.canvas.cd()
            plot.legend.Draw()

            # Do the subplot uncertainty shading
            data_no_errors = entries[0][0].Clone()
            cu.remove_th1_errors(data_no_errors)

            # Create hists for data with error region for ratio
            # Easiest way to get errors right is to do data (with 0 errors)
            # and divide by data (with errors), as if you had MC = data with 0 error

            # data_stat_ratio = data_no_errors.Clone()
            # data_stat_ratio.Divide(unfolded_hist_bin_stat_errors)
            # data_stat_ratio.SetFillStyle(3245)
            # data_stat_ratio.SetFillColor(self.plot_colours['unfolded_stat_colour'])
            # data_stat_ratio.SetLineWidth(0)
            # data_stat_ratio.SetMarkerSize(0)

            data_total_ratio = data_no_errors.Clone()
            data_total_ratio.Divide(entries[0][0])
            data_total_ratio.SetFillStyle(3254)
            data_total_ratio.SetFillColor(entries[0][1]['fill_color'])
            data_total_ratio.SetLineWidth(0)
            data_total_ratio.SetMarkerSize(0)

            if plot_kwargs.get("subplot_type", None) is not None:
                # now draw the data error shaded area
                # this is a bit hacky - basically draw them on the ratio pad,
                # then redraw the existing hists & line to get them ontop
                # note that we use "same" for all - this is to keep the original axes
                # (we may want to rethink this later?)
                plot.subplot_pad.cd()
                data_draw_opt = "E2 SAME"  # need SAME otherwise axis get rescaled
                # data_stat_ratio.Draw(data_draw_opt)
                data_total_ratio.Draw(data_draw_opt)
                plot.subplot_container.Draw("SAME " + draw_opt)
                plot.subplot_line.Draw()
                plot.canvas.cd()

        dirname = os.path.dirname(output_filename)
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
        plot.save(output_filename)

    except RuntimeError as e:
        print("Skipping:", e)


def get_projection_plot(h2d, start_val, end_val, cut_axis='y'):
    """Get projection plot from h2d, only covering from
    start_val (incl) to end_val (not incl) on the other axis (specified by cut_axis)"""

    cut_axis = cut_axis.lower()
    if cut_axis == "y":
        axis = h2d.GetYaxis()
    elif cut_axis == "x":
        axis = h2d.GetXaxis()
    else:
        raise RuntimeError("cut_axis must be x or y")

    # bin_edges = [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins()+2)]

    # # insert under and overflow
    # bin_edges.insert(0, -999999999999999)
    # bin_edges.append(999999999999999)

    # print(bin_edges)
    # bin_start = bisect.bisect_left(bin_edges, start_val)
    # bin_end = bisect.bisect_left(bin_edges, end_val)-1

    # print("bin_edges[bin_start]", bin_edges[bin_start])
    # print("bin_edges[bin_end]", bin_edges[bin_end])
    # # if end_val == bin_edges[bin_end]:
    #     # bin_end -= 1
    # # elif bin_end == bin_start:
    # #     bin_end = bin_start+1
    # # else:
    #     # bin_end -= 1

    start_bin = axis.FindBin(start_val)
    end_bin = axis.FindBin(end_val)
    # TH2::Projection* does first to last bin INCLUSIVE so we need to account for that
    if axis.GetBinLowEdge(end_bin) == end_val and end_bin != 0:
        end_bin -=1

    if cut_axis == "y":
        hproj = h2d.ProjectionX(ROOT.TUUID().AsString(), start_bin, end_bin, "eo")
        return hproj
    elif cut_axis == "x":
        hproj = h2d.ProjectionY(ROOT.TUUID().AsString(), start_bin, end_bin, "eo")
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
    canvas.SetRightMargin(0.16)
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


def do_all_exclusive_plots_comparison(sources,
                                      plot_dir="plots_dy_vs_qcd",
                                      dy_filename=qgc.DY_FILENAME,
                                      zpj_dirname="ZPlusJets_QG",
                                      qcd_filename=qgc.QCD_FILENAME,
                                      dj_dirname=None,  # backwards-compatible
                                      dj_cen_dirname="Dijet_QG_central_tighter",
                                      dj_fwd_dirname="Dijet_QG_forward_tighter",
                                      var_list=None,
                                      var_prepend="",
                                      pt_bins=None,
                                      title="",
                                      subplot_type=None,
                                      subplot_title=None,
                                      do_flav_tagged=True,
                                      show_region_labels=True,
                                      has_data=False,
                                      data_first=False,
                                      ofmt="pdf"):
    """Do 1D plots, comparing various sources. Looping over pt_bins, for each source, plots DY & QCD samples.

    `sources` is a list of dicts. Each dict has structure:
    {
        'root_dir': directory with ROOT file
        'label': append to put in legend after main label (optional)
        'dy_filename': name of ROOT file with Z+jets hists (optional)
        'qcd_filename': name of ROOT file with dijet hists (optional)
        'zpj_dirname': name of TDirectory with Z+jets hists (optional)
        'dj_cen_dirname': name of TDirectory with dijet central hists (optional)
        'dj_fwd_dirname': name of TDirectory with dijet foward hists (optional)
        'style': dict of style options for this Contribution (optional)
        'zpj_style': dict of style options for this Contribution if Z+jets (optional, applied after 'style' dict)
        'qcd_cen_style': dict of style options for this Contribution if dijet central (optional, applied after 'style' dict)
        'qcd_fwd_style': dict of style options for this Contribution if dijet forward (optional, applied after 'style' dict)
        'subplot': index of dict that should be the subplot
    }

    If zpj_dirname, dj_cen_dirname, dj_fwd_dirname blank, those contributions will not be plotted.

    `title` is appended after pt bin title.
    `show_region_labels` adds region label to legend
    """
    var_list = var_list or qgc.COMMON_VARS
    pt_bins = pt_bins or qgc.PT_BINS

    if isinstance(sources, dict):
        sources = [sources]

    for ang in var_list:

        v = "%s%s_vs_pt" % (var_prepend, ang.var)
        for (start_val, end_val) in pt_bins:
            entries_normal, entries_flav = [], []

            # some default sizings
            lw = 2
            msize = 1.1

            marker_cycle_settings = dict(only_cycle_filling=True)
            if len(sources) == 3:
                marker_cycle_settings = dict(cycle_filling=False)
            elif len(sources) > 3:
                marker_cycle_settings = dict(cycle_filling=True)

            # Get all plots, grouped by signal region
            if zpj_dirname:
                marker_zpj = cu.Marker(qgc.DY_MARKER)
                marker_iter = marker_zpj.cycle(**marker_cycle_settings)
                for ind, source in enumerate(sources):
                    if not source:
                        continue
                    marker = next(marker_iter)
                    h2d_dyj = grab_obj(os.path.join(source['root_dir'], source.get('dy_filename', dy_filename)),
                                       "%s/%s" % (source.get('zpj_dirname', zpj_dirname), v))
                    label_parts = [qgc.ZpJ_LABEL if show_region_labels else "", source.get('label', '')]
                    dy_kwargs = dict(line_color=qgc.DY_COLOUR, line_width=lw, fill_color=qgc.DY_COLOUR,
                                     label='\n'.join([l for l in label_parts if l]),
                                     marker_color=qgc.DY_COLOUR, marker_style=marker, marker_size=msize)
                    dy_kwargs.update(source.get('style', {}))
                    dy_kwargs.update(source.get('zpj_style', {}))
                    entries_normal.append((get_projection_plot(h2d_dyj, start_val, end_val), dy_kwargs))
                    del h2d_dyj

            if dj_cen_dirname:
                marker_qcd_cen = cu.Marker(qgc.QCD_CEN_MARKER, filled=True)
                marker_iter = marker_qcd_cen.cycle(**marker_cycle_settings)
                for ind, source in enumerate(sources):
                    if not source:
                        continue
                    marker = next(marker_iter)
                    # print(marker_qcd_cen.fill_state)
                    h2d_qcd_cen = grab_obj(os.path.join(source['root_dir'], source.get('qcd_filename', qcd_filename)),
                                           "%s/%s" % (source.get('dj_cen_dirname', dj_cen_dirname), v))
                    label_parts = [qgc.Dijet_CEN_LABEL if show_region_labels else "", source.get('label', '')]
                    qcd_cen_kwargs = dict(line_color=qgc.QCD_CEN_COLOUR, line_width=lw,
                                          fill_color=qgc.QCD_CEN_COLOUR,
                                          label='\n'.join([l for l in label_parts if l]),
                                          marker_color=qgc.QCD_CEN_COLOUR, marker_style=marker, marker_size=msize)
                    qcd_cen_kwargs.update(source.get('style', {}))
                    qcd_cen_kwargs.update(source.get('qcd_cen_style', {}))
                    entries_normal.append((get_projection_plot(h2d_qcd_cen, start_val, end_val), qcd_cen_kwargs))
                    del h2d_qcd_cen

            if dj_fwd_dirname:
                marker_qcd_fwd = cu.Marker(qgc.QCD_FWD_MARKER, filled=True)
                marker_iter = marker_qcd_fwd.cycle(**marker_cycle_settings)
                for ind, source in enumerate(sources):
                    if not source:
                        continue
                    marker = next(marker_iter)
                    h2d_qcd_fwd = grab_obj(os.path.join(source['root_dir'], source.get('qcd_filename', qcd_filename)),
                                           "%s/%s" % (source.get('dj_fwd_dirname', dj_fwd_dirname), v))
                    label_parts = [qgc.Dijet_FWD_LABEL if show_region_labels else "", source.get('label', '')]
                    qcd_fwd_kwargs = dict(line_color=qgc.QCD_FWD_COLOUR, line_width=lw,
                                          fill_color=qgc.QCD_FWD_COLOUR,
                                          label='\n'.join([l for l in label_parts if l]),
                                          marker_color=qgc.QCD_FWD_COLOUR, marker_style=marker, marker_size=msize)
                    qcd_fwd_kwargs.update(source.get('style', {}))
                    qcd_fwd_kwargs.update(source.get('qcd_fwd_style', {}))
                    entries_normal.append((get_projection_plot(h2d_qcd_fwd, start_val, end_val), qcd_fwd_kwargs))
                    del h2d_qcd_fwd

            if dj_dirname:
                for ind, source in enumerate(sources):
                    if not source:
                        continue
                    h2d_qcd = grab_obj(os.path.join(source['root_dir'], source.get('qcd_filename', qcd_filename)),
                                           "%s/%s" % (source.get('dj_dirname', dj_dirname), v))
                    label_parts = [qgc.Dijet_LABEL if show_region_labels else "", source.get('label', '')]
                    qcd_kwargs = dict(line_color=qgc.QCD_COLOUR, line_width=lw,
                                      fill_color=qgc.QCD_COLOUR,
                                      label='\n'.join([l for l in label_parts if l]),
                                      marker_color=qgc.QCD_COLOUR, marker_style=qgc.QCD_MARKER+ind, marker_size=msize)
                    qcd_kwargs.update(source.get('style', {}))
                    qcd_kwargs.update(source.get('qcd_style', {}))
                    entries_normal.append((get_projection_plot(h2d_qcd, start_val, end_val), qcd_kwargs))
                    del h2d_qcd

            do_flav_plot = do_flav_tagged and "flavour" not in v

            # Flav tagged plots
            if zpj_dirname and do_flav_plot:
                for ind, source in enumerate(sources):
                    if not source:
                        continue
                    h2d_dyj_q = grab_obj(os.path.join(source['root_dir'], source.get('dy_filename', dy_filename)),
                                         "%s/q%s" % (source.get('zpj_dirname', zpj_dirname), v))
                    dy_kwargs_q = dict(line_color=qgc.DY_COLOUR, line_width=lw, fill_color=qgc.DY_COLOUR,
                                       label=qgc.DY_ZpJ_QFLAV_LABEL + source.get('label', ''),
                                       # label=qgc.QCD_Dijet_QFLAV_LABEL + source.get('label', ''),
                                       marker_color=qgc.DY_COLOUR, marker_style=qgc.DY_MARKER+ind, marker_size=msize)
                    dy_kwargs_q.update(source.get('style', {}))
                    dy_kwargs_q.update(source.get('zpj_style', {}))
                    entries_flav.append((get_projection_plot(h2d_dyj_q, start_val, end_val), dy_kwargs_q))

            if dj_cen_dirname and do_flav_plot:
                for ind, source in enumerate(sources):
                    if not source:
                        continue
                    h2d_qcd_cen_g = grab_obj(os.path.join(source['root_dir'], source.get('qcd_filename', qcd_filename)),
                                             "%s/g%s" % (source.get('dj_cen_dirname', dj_cen_dirname), v))
                    qcd_cen_kwargs_g = dict(line_color=qgc.QCD_CEN_COLOUR, line_width=lw,
                                            fill_color=qgc.QCD_CEN_COLOUR,
                                            label=qgc.QCD_Dijet_CEN_GFLAV_LABEL + source.get('label', ''),
                                            marker_color=qgc.QCD_CEN_COLOUR, marker_style=qgc.QCD_CEN_MARKER+ind, marker_size=msize)
                    qcd_cen_kwargs_g.update(source.get('style', {}))
                    qcd_cen_kwargs_g.update(source.get('qcd_cen_style', {}))
                    entries_flav.append((get_projection_plot(h2d_qcd_cen_g, start_val, end_val), qcd_cen_kwargs_g))

                    h2d_qcd_cen_q = grab_obj(os.path.join(source['root_dir'], source.get('qcd_filename', qcd_filename)),
                                             "%s/q%s" % (source.get('dj_dirname', dj_dirname), v))
                    qcd_cen_kwargs_q = dict(line_color=qgc.QCD_CEN_COLOUR, line_width=lw,
                                            fill_color=qgc.QCD_CEN_COLOUR,
                                            label=qgc.QCD_Dijet_CEN_QFLAV_LABEL + source.get('label', ''),
                                            marker_color=qgc.QCD_CEN_COLOUR, marker_style=qgc.QCD_CEN_MARKER+ind+1, marker_size=msize)
                    qcd_cen_kwargs_q.update(source.get('style', {}))
                    qcd_cen_kwargs_q.update(source.get('qcd_cen_style', {}))
                    entries_flav.append((get_projection_plot(h2d_qcd_cen_q, start_val, end_val), qcd_cen_kwargs_q))

            if dj_fwd_dirname and do_flav_plot:
                for ind, source in enumerate(sources):
                    if not source:
                        continue
                    h2d_qcd_fwd_g = grab_obj(os.path.join(source['root_dir'], source.get('qcd_filename', qcd_filename)),
                                             "%s/g%s" % (source.get('dj_fwd_dirname', dj_dirname), v))
                    qcd_fwd_kwargs_g = dict(line_color=qgc.QCD_FWD_COLOUR, line_width=lw,
                                            fill_color=qgc.QCD_FWD_COLOUR,
                                            label=qgc.QCD_Dijet_FWD_GFLAV_LABEL + source.get('label', ''),
                                            marker_color=qgc.QCD_FWD_COLOUR, marker_style=qgc.QCD_FWD_MARKER+ind, marker_size=msize)
                    qcd_fwd_kwargs_g.update(source.get('style', {}))
                    qcd_fwd_kwargs_g.update(source.get('qcd_fwd_style', {}))
                    entries_flav.append((get_projection_plot(h2d_qcd_fwd_g, start_val, end_val), qcd_fwd_kwargs_g))

                    h2d_qcd_fwd_q = grab_obj(os.path.join(source['root_dir'], source.get('qcd_filename', qcd_filename)),
                                             "%s/q%s" % (source.get('dj_fwd_dirname', dj_dirname), v))
                    qcd_fwd_kwargs_q = dict(line_color=qgc.QCD_FWD_COLOUR, line_width=lw,
                                            fill_color=qgc.QCD_FWD_COLOUR,
                                            label=qgc.QCD_Dijet_FWD_QFLAV_LABEL + source.get('label', ''),
                                            marker_color=qgc.QCD_FWD_COLOUR, marker_style=qgc.QCD_FWD_MARKER+ind+1, marker_size=msize)
                    qcd_fwd_kwargs_q.update(source.get('style', {}))
                    qcd_fwd_kwargs_q.update(source.get('qcd_fwd_style', {}))
                    entries_flav.append((get_projection_plot(h2d_qcd_fwd_q, start_val, end_val), qcd_fwd_kwargs_q))

            # Now we've gone through all entries, we can update with the subplots (if any),
            # where the subplot has been specified by an index
            start_ind = -1
            for dn in [zpj_dirname, dj_cen_dirname, dj_fwd_dirname, dj_dirname]:
                if dn:
                    start_ind += 1
                    for ind, source in enumerate(sources):
                        if 'subplot' in source and isinstance(source['subplot'], int):
                            subplot_ind = source['subplot']
                            entries_normal[start_ind + ind][1]['subplot'] = entries_normal[start_ind+subplot_ind][0]
                    start_ind += ind

            # some default plot options
            rebin = 2
            v_lower = v.lower()
            if "multiplicity" in v_lower:
                rebin = 2
            elif "flavour" in v_lower or "thrust" in v_lower:
                rebin = 1
            elif 'ptd' in v_lower:
                rebin = 2

            xlim = None
            if "width" in v_lower or "thrust" in v_lower: # or "pTD" in v_lower:
                xlim = (0, 0.5)
            elif "multiplicity" in v_lower and "ak4" in sources[0]['root_dir'].lower():
                xlim = (0, 100)

            ylim = None
            if "flavour" in v_lower:
                ylim = (0, 1)
            # elif "LHA" in v:
            #     ylim = (0, 5)

            subplot = None
            num_entries_with_subplot = len([e for e in entries_normal if e[1].get('subplot', None)])
            if subplot_type != None and len(entries_normal) > 0 and num_entries_with_subplot == 0:
                # print("creating common subplot obj")
                subplot = entries_normal[0][0]
                if subplot.GetEntries() < 1:
                    subplot = entries_normal[1][0]
                    if subplot.GetEntries() < 1:
                        subplot = entries_normal[2][0]
                        if subplot.GetEntries() < 1:
                            subplot = None

            plot_title = "%d < p_{T}^{jet} < %d GeV" % (start_val, end_val)
            if title:
                plot_title += "\n%s" % title
            x_title = ang.name + " (" + ang.lambda_str + ")"
            do_comparison_plot(entries_normal,
                               output_filename="%s/ptBinned/%s_pt%dto%d.%s" % (plot_dir, v, start_val, end_val, ofmt),
                               rebin=rebin,
                               title=plot_title,
                               xtitle=x_title,
                               xlim=xlim,
                               ylim=ylim,
                               has_data=has_data,
                               data_first=data_first,
                               subplot=subplot,
                               subplot_type=subplot_type,
                               subplot_title=subplot_title,
                               subplot_limits=None)


            if do_flav_plot:
                subplot = None
                if subplot_type != None and len(entries_flav) > 0:
                    subplot = entries_flav[0][0]

                do_comparison_plot(entries_flav,
                                   output_filename="%s/ptBinned/%s_pt%dto%d_flavMatched.%s" % (plot_dir, v, start_val, end_val, ofmt),
                                   rebin=rebin,
                                   title=plot_title,
                                   xtitle=x_title,
                                   xlim=xlim,
                                   data_first=data_first,
                                   has_data=has_data,
                                   subplot=subplot,
                                   subplot_type=subplot_type,
                                   subplot_title=subplot_title,
                                   subplot_limits=(0, 2))


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


def calc_mean_error(hist, covariance_matrix=None):
    """Get error on the mean, using covariance matrix if available

    covariance_matrix should be a 2D numpy array for this bin
    If None, then assumes no correlation, and that the diagonal elements are square of bin error

    Uses the fact that if f = A+B,
    then sigma^2_f = sigma^2_A + sigma^2_B + 2*sigma_A*sigma_B
    """
    sum_sq = 0
    for ibin in range(1, hist.GetNbinsX()+1):
        # bin_h_i = hist.GetBinContent(ibin)
        bin_x_i = hist.GetBinCenter(ibin)

        for jbin in range(ibin, hist.GetNbinsX()+1):
            # bin_h_j = hist.GetBinContent(jbin)
            # Account for bin width scaling, and normalisation factor
            # total_bin_width = hist.GetBinWidth(ibin) * hist.GetBinWidth(jbin) * scale_factor * scale_factor
            this_err2 = 0
            if ibin == jbin:
                if covariance_matrix:
                    this_err2 = bin_x_i*bin_x_i*covariance_matrix[ibin-1, ibin-1]
                else:
                    this_err2 = bin_x_i*bin_x_i*pow(hist.GetBinError(ibin), 2)
            elif covariance_matrix:
                # pass
                bin_x_j = hist.GetBinCenter(jbin)
                this_err2 = 2*bin_x_i*bin_x_j*covariance_matrix[ibin-1, jbin-1]  # should this be sqrt?

            sum_sq += this_err2

    # scale for total entries
    # print('final', sum_sq, hist.Integral())
    return sqrt(sum_sq) / hist.Integral()


def do_mean_rms_summary_plot(entries,
                             bins,
                             output_filename,
                             var_label="",
                             xlim=None,
                             ylim_mean=None,
                             ylim_rms=None,
                             region_title="",
                             has_data=True,
                             numerator_label="MC",
                             denominator_label="Data"):
    """Create summary plot from 1D hists using matplotlib

    Parameters
    ----------
    entries : list[list[(TH1, dict)]]
        Histograms to be used.
        entries[i][j] is the jth sample for the ith pt bin
    bins : list[(float, float)]
        Bin edges
    output_filename : str
        Output filename
    var_label : str, optional
        Variable name to put in title
    xlim : None, optional
        Limit range on x axis
    ylim_mean : None, optional
        Limit range on y axis for mean
    ylim_rms : None, optional
        Limit range on y axis for rms
    region_title : str, optional
        Name of region to put in title
    has_data : bool, optional
        True if data is in plot, False puts "Simulation" in top-left text
    denominator_label : str, optional
        Label for ratio titles for denominator
    numerator_label : str, optional
        Label for ratio titles for numerator
    """
    # Get the statistics from our 1D hists
    # only do the bins we need, otherwise axis limits go awry
    if xlim:
        bins = [b for b in bins if b[1] <= xlim[1]]

    # same structure as entries:
    # means[i][j] is the jth dataset for the ith pt bin
    data_means = []
    data_mean_errs = []
    data_rms = []
    data_rms_errs = []

    for pt_ind, pt_bin_entry in enumerate(entries):  # iterate over pt bins
        if pt_ind == len(bins):
            break
        this_mean = [] # hold data for all datasets, for this pt bin
        this_mean_err = []
        this_rms = []
        this_rms_err = []

        for dataset_entry in pt_bin_entry:  # iterate over datasets
            # here deal with missing, or bad hists
            # also skip low stats hists
            if not dataset_entry or (dataset_entry[0].GetEntries() < 100):
                this_mean.append(None)
                this_mean_err.append(None)
                this_rms.append(None)
                this_rms_err.append(None)
                continue

            hist = dataset_entry[0]

            contents, errors = cu.th1_to_ndarray(hist)
            # this hist is normalised to unity, so contents.sum() = 1
            centers = cu.get_th1_bin_centers(hist)
            areas, centers = metrics.hist_values_to_uarray(bin_areas=contents, bin_centers=centers, bin_errors=errors)
            mean_u = metrics.calc_mean_ucert(areas, centers)
            mean, mean_err = mean_u.nominal_value, mean_u.std_dev
            rms_u = metrics.calc_rms_ucert(areas, centers)
            rms, rms_err = rms_u.nominal_value, rms_u.std_dev

            this_mean.append(mean)
            this_mean_err.append(mean_err)
            this_rms.append(rms)
            this_rms_err.append(rms_err)

        data_means.append(np.array(this_mean, dtype=float))
        data_mean_errs.append(np.array(this_mean_err, dtype=float))
        data_rms.append(np.array(this_rms, dtype=float))
        data_rms_errs.append(np.array(this_rms_err, dtype=float))

    # convert to numpy array to maniuplate easier
    data_means = np.array(data_means)
    data_mean_errs = np.array(data_mean_errs)
    data_rms = np.array(data_rms)
    data_rms_errs = np.array(data_rms_errs)

    # Start by setting up figure
    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(12, 10), sharex=True)
    plt.subplots_adjust(left=0.1, right=0.8, top=0.85, hspace=0.2)
    mean_plot_ind, mean_ratio_plot_ind, rms_plot_ind, rms_ratio_plot_ind = 0, 1, 2, 3

    # Turn on minor ticks
    for ax in axes.flat:
        ax.yaxis.set_minor_locator(plticker.AutoMinorLocator())

    n_datasets = len(entries[0])
    bin_widths = [(x[1]-x[0])/(n_datasets) for x in bins]

    # Convert ROOT colors to hex for matplotlib
    cols = []
    for dataset_entry in entries[0]:
        cols.append(ROOT.gROOT.GetColor(dataset_entry[1]['line_color']).AsHexString())

    markers = ['x', 'o', 'd', '^']

    x_pos = [0.5*sum(x) for x in bins]
    bin_widths = [0.5*(x[1]-x[0]) for x in bins]
    plot_obj = []

    # Create & plot graphs
    ref_data_mean = None
    ref_data_rms = None
    for dataset_ind, dataset_entry in enumerate(entries[0]):
        common_plot_settings = dict(
            xerr=bin_widths,
            marker=markers[dataset_ind],
            label=extract_sample_name(dataset_entry[1]['label']),
            color=cols[dataset_ind],
            markersize=8,
            linewidth=0,
            elinewidth=2
        )

        # Plot means
        this_data = data_means[:, dataset_ind]
        this_data_err = data_mean_errs[:, dataset_ind]
        this_plot_obj = axes[mean_plot_ind].errorbar(x_pos, this_data,
                                                     yerr=this_data_err,
                                                     **common_plot_settings)
        plot_obj.append(this_plot_obj)

        # Plot RMS
        this_data = data_rms[:, dataset_ind]
        this_data_err = data_rms_errs[:, dataset_ind]
        axes[rms_plot_ind].errorbar(x_pos, this_data,
                                    yerr=this_data_err,
                                    **common_plot_settings)

        if dataset_ind == 0:
            ref_data_mean = data_means[:, dataset_ind]
            ref_data_rms = data_rms[:, dataset_ind]
        else:
            # Plot ratio of means
            this_data_mean = data_means[:, dataset_ind] / ref_data_mean
            axes[mean_ratio_plot_ind].errorbar(x_pos, this_data_mean,
                                               **common_plot_settings)
            # Plot ratio of RMS
            this_data_rms = data_rms[:, dataset_ind] / ref_data_rms
            axes[rms_ratio_plot_ind].errorbar(x_pos, this_data_rms,
                                              **common_plot_settings)

    # Add legend
    axes[mean_plot_ind].legend(handles=plot_obj,
                               loc=(1.01, 0.05),
                               fancybox=False,
                               edgecolor='white')

    # Set axis limits for mean
    axes[mean_plot_ind].set_xlim(*xlim)
    if ylim_mean:
        axes[mean_plot_ind].set_ylim(*ylim_mean)
    # if axes[mean_plot_ind].get_ylim()[0] > 0:
    #     axes[mean_plot_ind].set_ylim(bottom=0)
    axes[-1].set_xlim(*xlim)

    # Set axis limits for rms
    axes[rms_plot_ind].set_xlim(*xlim)
    if ylim_rms:
        axes[rms_plot_ind].set_ylim(*ylim_rms)
    # if axes[rms_plot_ind].get_ylim()[0] > 0:
    #     axes[rms_plot_ind].set_ylim(bottom=0)
    axes[-1].set_xlim(*xlim)

    plt.xscale('log')

    # Cap mean ratio maximum y to 2
    ratio_ylim = axes[mean_ratio_plot_ind].get_ylim()
    lim = 2
    if ratio_ylim[1] > lim:
        axes[mean_ratio_plot_ind].set_ylim(top=lim)
    # Add line at 1
    axes[mean_ratio_plot_ind].axhline(1.0, color='gray', linestyle='dashed')

    # Cap rms ratio maximum y to 2
    ylim = axes[rms_ratio_plot_ind].get_ylim()
    lim = 2
    if ylim[1] > lim:
        axes[rms_ratio_plot_ind].set_ylim(top=lim)
    # Add line at 1
    axes[rms_ratio_plot_ind].axhline(1.0, color='gray', linestyle='dashed')

    # Set overall plot title
    new_var_label = var_label.replace("#", "\\")
    axes[mean_plot_ind].set_title("Summary for {} region:\n{}".format(region_title, new_var_label), pad=40)

    # Add CMS text
    text_height = 1.1
    axes[mean_plot_ind].text(0, text_height, "CMS", fontweight='bold', horizontalalignment='left', transform=axes[mean_plot_ind].transAxes)
    # have to do bold bit separately
    if has_data:
        axes[mean_plot_ind].text(0, text_height, "         Preliminary", horizontalalignment='left', transform=axes[mean_plot_ind].transAxes)
    else:
        axes[mean_plot_ind].text(0, text_height, "         Preliminary Simulation", horizontalalignment='left', transform=axes[mean_plot_ind].transAxes)
    # Add lumi text
    axes[mean_plot_ind].text(1, text_height, "35.9 $\\mathrm{fb}^{\\mathrm{-1}}$ (13 TeV)", horizontalalignment='right', transform=axes[mean_plot_ind].transAxes)

    # Set axis titles
    axes[mean_plot_ind].set_ylabel('Mean $\\pm$ error')
    axes[mean_ratio_plot_ind].set_ylabel("Ratio of means\n(%s / %s)" % (numerator_label, denominator_label))
    axes[rms_plot_ind].set_ylabel('RMS $\\pm$ error')
    axes[rms_ratio_plot_ind].set_ylabel("Ratio of RMS\n(%s / %s)" % (numerator_label, denominator_label))
    axes[-1].set_xlabel('$p_{\mathrm{T}}^{\mathrm{jet}}$ [GeV]')

    # Draw bin delineators
    for ax in axes.flat:
        ax.tick_params(axis='x', which='major', direction='in', top=True, bottom=True, length=10)
        ax.tick_params(axis='x', which='minor', direction='in', top=True, bottom=True, length=7)
        ax.tick_params(axis='y', which='major', direction='in', left=True, right=True, length=10)
        ax.tick_params(axis='y', which='minor', direction='in', left=True, right=True, length=7)

        this_ylim = ax.get_ylim()
        for ind, this_bin in enumerate(bins):
            ax.axvline(this_bin[0], color='lightgray', linestyle='dashed')

    plt.tight_layout()

    # Save to file
    odir = os.path.dirname(os.path.abspath(output_filename))
    if not os.path.isdir(odir):
        os.makedirs(odir)
    plt.savefig(output_filename)
    plt.close(fig)
    # plt.clf()


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
        Contribution(hist_purity, label="Purity (gen correct)", line_color=col_purity, marker_color=col_purity),
        Contribution(hist_stability, label="Stability (reco correct)", line_color=col_stability, marker_color=col_stability),
        Contribution(hist_xfer_down, label="-1 reco bin", line_color=col_xfer_down, marker_color=col_xfer_down, line_style=2),
        Contribution(hist_xfer_down2, label="-2 reco bin", line_color=col_xfer_down2, marker_color=col_xfer_down2, line_style=3),
        # Contribution(hist_xfer_down3, label="3 lower reco bin", line_color=col_xfer_down3, marker_color=col_xfer_down3),
        Contribution(hist_xfer_up, label="+1 reco bin", line_color=col_xfer_up, marker_color=col_xfer_up, line_style=2),
        Contribution(hist_xfer_up2, label="+2 reco bin", line_color=col_xfer_up2, marker_color=col_xfer_up2, line_style=3),
        # Contribution(hist_xfer_up3, label="3 higher reco bin", line_color=col_xfer_up3, marker_color=col_xfer_up3),
    ]
    return contributions


def make_migration_summary_plot(h2d_renorm_x, h2d_renorm_y, xlabel, output_filename, title="", log_var=False):
    contributions = migration_plot_components(h2d_renorm_x, h2d_renorm_y, xlabel)
    for c in contributions:
        c.obj.SetLineWidth(2)
    binning = [h2d_renorm_x.GetXaxis().GetBinLowEdge(bin_ind) for bin_ind in range(1, h2d_renorm_x.GetNbinsX()+2)]
    xlim = [binning[0], binning[-1]]
    plot = Plot(contributions, what='hist', xlim=xlim, ylim=[1e-3, 10], xtitle=xlabel, has_data=False, title=title)
    plot.default_canvas_size = (700, 600)
    y1 = 0.15
    plot.legend.SetX1(0.55)
    plot.legend.SetY1(y1)
    plot.legend.SetY2(y1+0.2)
    plot.plot("NOSTACK HISTE")
    plot.legend.SetFillStyle(1001)
    plot.legend.SetFillColorAlpha(ROOT.kWhite, 0.75)
    if log_var:
        plot.set_logx()
    plot.set_logy()
    plot.main_pad.cd()
    lines = []
    for val in [1, 0.5, 1e-1, 1e-2]:
        line = ROOT.TLine(xlim[0], val, xlim[1], val)
        line.SetLineStyle(2)
        line.SetLineColor(ROOT.kGray+2)
        lines.append(line)
        line.Draw("same")
    plot.save(output_filename)


def normalise_hist(h):
    if h.Integral() > 0:
        h.Scale(1./h.Integral())


def hist_divide_bin_width(h):
    """Create copy of hist, but each bin's contents is divide by the bin width"""
    h_new = h.Clone(h.GetName()+"DivideBinWidth")
    h_new.Scale(1., 'width')
    # if any bin has 0 entries before, it will now have nan, so we need to manually fix that
    # arghhhhhhhhhhhhhhhhhh
    for i in range(1, h_new.GetNbinsX()+1):
        if np.isnan(h_new.GetBinContent(i)):
            h_new.SetBinContent(i, 0)
            h_new.SetBinError(i, 0)
    h_new.SetEntries(h.GetEntries())  # needed as the default is to replace with integral
    return h_new


def normalise_hist_divide_bin_width(h):
    """Create copy of hist, but each bin's contents is divide by the bin width"""
    h_new = h.Clone(h.GetName()+"DivideBinWidth")
    normalise_hist(h_new)
    h_new = hist_divide_bin_width(h_new)
    return h_new


def calc_auto_xlim(entries):
    """Figure out x axis range that includes all non-0 bins from all Contributions or TH1s"""
    if len(entries) == 0:
        return None
    if not isinstance(entries[0], Contribution) and not isinstance(entries[0], ROOT.TH1):
        raise TypeError("calc_auto_xlim: `entries` should be a list of Contributions or TH1s")

    x_min = 9999999999
    x_max = -9999999999
    for ent in entries:
        if isinstance(ent, ROOT.TH1):
            obj = ent
        else:
            obj = ent.obj
        if not isinstance(obj, ROOT.TH1):
            raise TypeError("Cannot handle obj type %s in calc_auto_xlim" % type(obj))
        xax = obj.GetXaxis()
        nbins = obj.GetNbinsX()
        found_min = False
        for i in range(1, nbins+1):
            val = obj.GetBinContent(i)
            if not found_min:
                # find the first non-empty bin
                if val != 0:
                    x_min = min(xax.GetBinLowEdge(i), x_min)
                    found_min = True
        for i in range(nbins, 0, -1):
            val = obj.GetBinContent(i)
            # find the last empty bin, work from rightmost bin
            if val != 0:
                x_max = max(xax.GetBinLowEdge(i+1), x_max)
                break

    if x_max < x_min:
        # didn't find high 0 bin, so set to last bin
        x_max = xax.GetBinLowEdge(nbins+1)
    if x_max > x_min:
        return (x_min, x_max)
    else:
        print("calc_auto_xlim: max < min wtf")
        # default x max
        return None


def do_fancy_legend(contributions, plot, use_splitline=True):
    """Fill legend with contributions, using TGraph instead to get marker done correctly

    Actually puts horizontal error bar in legend, whereas for a TH1 it doesn't
    """
    dummy_graphs = []
    dummy_gr = ROOT.TGraphErrors(1, array('d', [1]), array('d', [1]), array('d', [1]), array('d', [1]))
    for cont in contributions:
        new_gr = dummy_gr.Clone()
        cont.update_obj_styling(new_gr)
        dummy_graphs.append(new_gr)

        if '\n' not in cont.label:
            plot.legend.AddEntry(new_gr, cont.label, cont.leg_draw_opt)
        else:
            if not use_splitline:
                for label_ind, label_part in enumerate(cont.label.split("\n")):
                    obj = new_gr if label_ind == 0 else 0
                    draw_opt = cont.leg_draw_opt if label_ind == 0 else ""
                    plot.legend.AddEntry(obj, label_part, draw_opt)
            else:
                part1, part2 = cont.label.split("\n")
                new_label = "#splitline{%s}{%s}" % (part1, part2)
                plot.legend.AddEntry(new_gr, new_label, cont.leg_draw_opt)
    return dummy_graphs
