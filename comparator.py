"""
Classes for use when comparing functions/graphs

Example usage:

>>> A = Contribution("a.root", "li1corr_eta_0_0.348")
>>> B = Contribution("b.root", "li1corr_eta_0_0.348")
>>> p = Plot([A, B], "graphfunction", xlim=[0, 50])
>>> p.plot()
>>> p.save("AB.pdf")

"""


import os
import numpy as np
import ROOT
import common_utils as cu
from MyStyle import My_Style


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
# ROOT.gStyle.SetPalette(55)
My_Style.cd()


class MultiFunc(object):
    """Class to handle multiple TF1s, like a TMultiGraph, because ROOT doesn't.

    NOT SOPHISTICATED, JUST A CONTAINER.
    """
    def __init__(self):
        self.funcs = []

    def Add(self, func):
        self.funcs.append(func)

    def Draw(self, option=None):
        """Draw all the TF1s, adjust the axis sizes properly"""
        # if ROOT.gPad:
        #     if not ROOT.gPad.IsEditable():
        #         ROOT.gROOT.MakeDefCanvas()
        #     if "a" in option.lower():
        #         ROOT.gPad.Clear()
        # Draw, and figure out max/min of axes for some auto-ranging
        x_low, x_high = ROOT.Double(0), ROOT.Double(0)
        x_min, x_max = 999, -999
        y_min, y_max = 999, -999
        option = option or ""
        for i, f in enumerate(self.funcs):
            if i == 0:
                f.Draw(option)
            else:
                f.Draw(option + "SAME")

            f.GetRange(x_low, x_high)
            x_min = x_low if x_low < x_min else x_min
            x_max = x_high if x_high > x_max else x_max
            y_min = f.GetMinimum(x_low, x_high) if f.GetMinimum(x_low, x_high) < y_min else y_min
            y_max = f.GetMaximum(x_low, x_high) if f.GetMaximum(x_low, x_high) > y_max else y_max
        if x_max < x_min:
            raise Exception("MultiFunc: x_min > x_max")
        if y_max < y_min:
            raise Exception("MultiFunc: y_min > y_max")

        self.funcs[0].GetXaxis().SetRangeUser(x_min * 0.95, x_max * 1.05)
        self.funcs[0].GetYaxis().SetRangeUser(y_min * 0.95, y_max * 1.05)

    def Mod(self):
        """If we want to do any sort of styling, need to access the first drawn func"""
        return self.funcs[0]


def grab_obj(file_name, obj_name):
    """Get object names obj_name from ROOT file file_name"""
    # TODO: checks!
    input_file = cu.open_root_file(file_name)
    obj = cu.get_from_tfile(input_file, obj_name)
    # print("Getting", obj_name, "from", file_name)
    if isinstance(obj, (ROOT.TH1, ROOT.TGraph, ROOT.TH2)):
        # THIS ORDER IS VERY IMPORTANT TO AVOID MEMORY LEAKS
        new_obj = obj.Clone(ROOT.TUUID().AsString())
        new_obj.SetDirectory(0)
        input_file.Close()
        return new_obj
    else:
        return obj

class Contribution(object):
    """Basic class to handle information about one contribution to a canvas."""

    def __init__(self, obj, label=None,
                 line_width=1, line_color=ROOT.kRed, line_style=1,
                 fill_color=ROOT.kRed, fill_style=0,
                 marker_size=1, marker_color=ROOT.kRed, marker_style=1,
                 subplot_line_width=None, subplot_line_color=None, subplot_line_style=None,
                 subplot_fill_color=None, subplot_fill_style=None,
                 subplot_marker_size=None, subplot_marker_color=None, subplot_marker_style=None,
                 fit_match_style=True,
                 normalise_hist=False, divide_bin_width=False,
                 rebin_hist=None, subplot=None):
        """
        Parameters
        ----------
        obj : TH1, TGraph, ...
            Object to be plotted
        label : str
            Title of contribution, to be used in legend.
        line_width : int
            Width of line.
        line_color : int or ROOT color
            Color of line.
        line_style : int
            Line style.
        fill_color : int or ROOT color
            Color of fill.
        fill_style : int
            Fill style.
        marker_size : int
            Size of marker
        marker_color : int or ROOT color
            Color of markers.
        marker_style : int
            Marker style.
        fit_match_style : bool
            If True, apply line_color etc to fit function as well
        normalise_hist : bool
            If a histogram, normalise so integral = 1
        divide_bin_width : bool, optional
            If True, divide by bin width when normalising.
            Only makes sense when used with normalise_hist=True
        rebin_hist : int, None
            If a histogram, specify the number of bins to be grouped together.
        subplot : ROOT.TObject
            Object to use for subplot, if there is one
        """
        self.obj = obj
        self.label = label or ""
        self.line_width = line_width
        self.line_color = line_color
        self.line_style = line_style
        self.fill_color = fill_color
        self.fill_style = fill_style
        self.marker_size = marker_size
        self.marker_color = marker_color
        self.marker_style = marker_style
        self.subplot = subplot

        # subplot-specific stylings - inherit from main stylings by default
        self.subplot_line_width = subplot_line_width or self.line_width
        self.subplot_line_color = subplot_line_color or self.line_color
        self.subplot_line_style = subplot_line_style or self.line_style
        self.subplot_fill_color = subplot_fill_color or self.fill_color
        self.subplot_fill_style = subplot_fill_style or self.fill_style
        self.subplot_marker_size = subplot_marker_size or self.marker_size
        self.subplot_marker_color = subplot_marker_color or self.marker_color
        self.subplot_marker_style = subplot_marker_style or self.marker_style

        self.obj.SetLineWidth(self.line_width)
        self.obj.SetLineColor(self.line_color)
        self.obj.SetLineStyle(self.line_style)
        self.obj.SetFillColor(self.fill_color)
        self.obj.SetFillStyle(self.fill_style)
        self.obj.SetMarkerSize(self.marker_size)
        self.obj.SetMarkerColor(self.marker_color)
        self.obj.SetMarkerStyle(self.marker_style)

        # Match fit to hist styles
        if fit_match_style and self.obj.GetListOfFunctions().GetSize() == 1:
            func = self.obj.GetListOfFunctions().At(0)
            func.SetLineStyle(line_style)
            func.SetLineWidth(line_width)
            func.SetLineColor(line_color)

        if rebin_hist and rebin_hist != 1:
            self.obj.Rebin(rebin_hist) # Does this handle 2D hists?
        if normalise_hist and obj.Integral() != 0:
            if divide_bin_width:
                # check if all same bin width
                if sum([obj.GetBinWidth(i)/obj.GetBinWidth(1) for i in range(2, obj.GetNbinsX())]) < 1E-5:
                    self.obj.Scale(1./(obj.GetBinWidth(1) * obj.Integral()))
                else:
                    # deal with varaible bin width, much more involved
                    raise RuntimeError("Can't deal with variable bin width and divide_bin_width")
            else:
                self.obj.Scale(1./(obj.Integral()))
        if isinstance(self.obj, ROOT.TH1) or isinstance(self.obj, ROOT.TH2):
            self.obj.SetDirectory(0)

    def __eq__(self, other):
        return self.obj == other.obj


class ZeroContributions(Exception):
    """Own exception for when there are 0 Contributions (or they all have 0 entries)"""
    def __init__(self, message):
        super(ZeroContributions, self).__init__()
        self.message = message


class Plot(object):
    """
    Basic class to handle information about one plot,
    which can have several contributions.
    """

    def __init__(self, contributions=None, what="graph",
                 title=None, xtitle=None, ytitle=None, xlim=None, ylim=None,
                 legend=True, extend=False,
                 subplot=None, subplot_type=None, subplot_title=None, subplot_limits=None,
                 has_data=True):
        """
        contributions: list
            List of Contribution objects.
        what: str
            Options are "graph", "function", "both".
        title: str
            Title of plot.
        xtitle: str
            X axis title.
        ytitle: str
            Y axis title.
        xlim: list
            Limits of x axis. If None then determines suitable limits.
        ylim: list
            Limits of y axis. If None then determines suitable limits.
        legend: bool
            Include legend on plot.
        extend: bool
            Extend functions to cover the whole x axis range.
        subplot : object, None
            If not None, draws a plot below the main hist with
            all plots compared to the object provided as the argument here.
        subplot_type : str
            The method of comparison in the subplot: ratio, difference, or ddelta/dlambda
        subplot_limits : (float, float), optional
            Set hard limits on subplot y axis range. If None, will choose some
            vaguely sensible ones
        """
        self.contributions = contributions if contributions else []
        self.contributions_objs = []
        options = ['hist', 'graph', 'function', 'both']
        if what not in options:
            raise RuntimeError("`what` argument must be one of %s" % options)
        self.plot_what = what
        self.title = title or ""
        self.title_start_y = 0.87
        self.title_diff_y = 0.04
        self.title_left_offset = 0.03
        self.title_font_size = 0.03
        self.cms_text_font_size = 0.035
        self.xtitle = xtitle
        self.ytitle = ytitle
        self.xlim = xlim
        self.ylim = ylim
        self.y_padding_max_linear = 1.6  # factor to auto extend y upper limit for linear scale
        self.y_padding_max_log = 10  # factor to auto extend y upper limit for log scale
        self.y_padding_min_linear = 1.4 # factor to auto extend y lower limit for linear scale
        self.y_padding_min_log = 0.1  # factor to auto extend y lower limit for log scale
        self.do_legend = legend
        self.legend = ROOT.TLegend(0.65, 0.6, 0.94, 0.85) if legend else None
        self.reverse_legend = False
        self.do_extend = extend
        self.container = None
        self.canvas = None
        self.default_canvas_size = (600, 800)
        # self.default_canvas_size = (800, 600)
        self.right_margin = 0.04
        self.left_margin = 0.12 # use ROOT default
        self.text_left_offset = self.left_margin * (0.055/0.12)
        self.text_left_offset = 0
        self.top_margin = 0.1
        self.main_pad = None
        self.subplot = subplot
        if subplot_type and subplot_type not in ['ratio', 'diff', "ddelta"]:
            raise RuntimeError("subplot_type must be one of None, ratio, diff, or ddelta")
        self.subplot_type = subplot_type
        if not self.subplot_type:
            self.subplot = None
        self.subplot_container = None
        self.subplot_contributions = []
        self.subplot_pad = None
        self.subplot_limits = subplot_limits
        self.subplot_pad_height = 0.32
        self.subplot_pad_fudge = 0.01  # to get non-overlapping subplot axis
        self.subplot_line = None  # need this to remain visible...
        self.subplot_title = subplot_title
        self.subplot_line_style = 2
        self.subplot_line_width = 2
        self.subplot_line_color = ROOT.kBlack
        self.is_preliminary = True
        self.has_data = has_data

    def add_contribution(self, *contribution):
        """Add Contribution to Plot. Can be single item or list."""
        self.contributions.extend(*contribution)

    def _create_container(self):
        """Create a container object for lots of the same TObject"""
        if self.plot_what in ["graph", "both"]:
            self.container = ROOT.TMultiGraph(ROOT.TUUID().AsString(), "")
            self.subplot_container = ROOT.TMultiGraph(ROOT.TUUID().AsString(), "")
        elif self.plot_what == "function":
            self.container = MultiFunc()
            self.subplot_container = MultiFunc()
        elif self.plot_what == "hist":
            self.container = ROOT.THStack(ROOT.TUUID().AsString(), "")
            self.subplot_container = ROOT.THStack(ROOT.TUUID().AsString(), "")

        if self.container:
            ROOT.SetOwnership(self.container, False)
        if self.subplot_container:
            ROOT.SetOwnership(self.subplot_container, False)


    def _populate_container_and_legend(self):
        """Add objects to the container, and to the legend"""
        for contrib in self.contributions:
            if self.plot_what == "graph":
                # if drawing only graph, we need to remove the fit if there is one
                if (contrib.obj.GetListOfFunctions().GetSize() > 0):
                    contrib.obj.GetListOfFunctions().Remove(contrib.obj.GetListOfFunctions().At(0))
            else:
                if isinstance(contrib.obj, ROOT.TH1) and contrib.obj.GetEntries() == 0:
                    continue
                # if drawing function, extend range as per user's request
                if self.do_extend:
                    if self.plot_what == 'function':
                        contrib.obj.SetRange(self.xlim[0], self.xlim[1])
                    elif self.plot_what == "both":
                        contrib.obj.GetListOfFunctions().At(0).SetRange(self.xlim[0], self.xlim[1])
            self.container.Add(contrib.obj)
            self.contributions_objs.append(contrib.obj)

            # Add contributions for the subplot
            if self.subplot_type:
                if self.plot_what != 'hist':
                    print("WARNING: cannot do suplot for anything but hists, skipping")
                    self.subplot_type = None
                    continue
                new_hist = None
                if self.subplot:
                    # Use one reference object for all entries
                    subplot_obj = self.subplot.obj.Clone()
                    if contrib != self.subplot:
                        new_hist = contrib.obj.Clone()
                        if (self.subplot_type == "ratio"):
                            new_hist.Divide(subplot_obj)
                        elif (self.subplot_type == "diff"):
                            new_hist.Add(subplot_obj, -1.)
                        elif (self.subplot_type == "ddelta"):
                            # Do the differntial delta spectrum, see 1704.03878
                            new_hist.Add(subplot_obj, -1.)
                            new_hist.Multiply(new_hist)
                            sum_hist = contrib.obj.Clone()
                            sum_hist.Add(subplot_obj)
                            new_hist.Divide(sum_hist)
                            new_hist.Scale(0.5)
                elif contrib.subplot is not None:
                    new_hist = contrib.obj.Clone()
                    ref_obj = contrib.subplot
                    if (self.subplot_type == "ratio"):
                        new_hist.Divide(ref_obj)
                    elif (self.subplot_type == "diff"):
                        new_hist.Add(ref_obj, -1.)
                    elif (self.subplot_type == "ddelta"):
                        # Do the differntial delta spectrum, see 1704.03878
                        new_hist.Add(ref_obj, -1.)
                        new_hist.Multiply(new_hist)
                        sum_hist = contrib.obj.Clone()
                        sum_hist.Add(ref_obj)
                        new_hist.Divide(sum_hist)
                        new_hist.Scale(0.5)
                if new_hist:
                    # style it
                    new_hist.SetLineWidth(contrib.subplot_line_width)
                    new_hist.SetLineColor(contrib.subplot_line_color)
                    new_hist.SetLineStyle(contrib.subplot_line_style)
                    new_hist.SetFillColor(contrib.subplot_fill_color)
                    new_hist.SetFillStyle(contrib.subplot_fill_style)
                    new_hist.SetMarkerSize(contrib.subplot_marker_size)
                    new_hist.SetMarkerColor(contrib.subplot_marker_color)
                    new_hist.SetMarkerStyle(contrib.subplot_marker_style)
                    self.subplot_container.Add(new_hist)
                    self.subplot_contributions.append(new_hist)

        # do legend separately, since we might want to reverse the order
        # e.g stacking
        iterable = self.contributions
        if self.reverse_legend:
            iterable = self.contributions[::-1]
        for contrib in iterable:
            if self.do_legend and contrib.label != "" and contrib.label is not None:
                opt = "LP"
                if self.plot_what == "hist":
                    if contrib.marker_size != 0:
                        opt += "P"
                    if contrib.line_width != 0:
                        opt += "L"
                    if contrib.fill_style != 0:
                        opt += "F"
                # Split text by newline \n
                # Add an entry for each line
                for i, substr in enumerate(contrib.label.split("\n")):
                    if i == 0:
                        self.legend.AddEntry(contrib.obj, substr, opt)
                    else:
                        self.legend.AddEntry(0, substr, "")


    def _style_legend(self):
        # self.legend.SetBorderSize(0)
        self.legend.SetFillStyle(0)

    def _rescale_plot_labels(self, container, factor):
        # Eurgh, why does ROOT scale all these sizes?
        xax = container.GetXaxis()
        xax.SetLabelSize(xax.GetLabelSize()/factor)
        xax.SetTitleSize(xax.GetTitleSize()/factor)
        xax.SetTitleOffset(xax.GetTitleOffset()*factor)  # doesn't seem to work?
        xax.SetTickLength(xax.GetTickLength()/factor)

        yax = container.GetYaxis()
        yax.SetLabelSize(yax.GetLabelSize()/factor)
        yax.SetTitleSize(yax.GetTitleSize()/factor)
        # magic numbers: 0.1 is the default margin, but scaling against that gives too much, so we knock it down by a bit
        yax.SetTitleOffset(yax.GetTitleOffset()*factor*(0.7*self.left_margin/0.1))
        # container.GetYaxis().SetTickLength(0.03/factor)

    def set_logx(self, state=True, do_more_labels=True):
        # Call AFTER plot()
        self.main_pad.SetLogx(int(state))
        if self.subplot_pad:
            self.subplot_pad.SetLogx(int(state))
        if do_more_labels:
            if self.container:
                ax = self.container.GetXaxis()
                if ax:
                    ax.SetMoreLogLabels()
            if self.subplot_container:
                ax = self.subplot_container.GetXaxis()
                if ax:
                    ax.SetMoreLogLabels()

    def set_logy(self, state=True, do_more_labels=True, override_check=False):
        # Call AFTER plot()
        try:
            y_low = self.container.GetHistogram().GetMinimum()
            if y_low < 0:
                if override_check:
                    print("Warning in set_logy: minimum < 0, is %g" % y_low)
                else:
                    print("Cannot set_logy as minimum is %g" % y_low)
                    return
            self.main_pad.SetLogy(int(state))
        except AttributeError as e:
            print("")
            print(e)
            print("")
            print("You should run set_logy() after plot()")
            print("")
            raise

        # Don't make subplot log y...we never want that
        if self.container:
            ax = self.container.GetYaxis()
            if ax and state and do_more_labels:
                ax.SetMoreLogLabels()

        # update y limits since padding different for log/lin
        if not self.ylim:
            self._set_automatic_y_limits()

    def set_logz(self, state=True):
        self.main_pad.SetLogz(int(state))

    def get_modifier(self):
        if self.plot_what != 'function':
            modifier = self.container
        else:
            modifier = self.container.Mod()
        return modifier

    def _set_automatic_y_limits(self):
        # urgh why doesnt THStack.GetMaximum() return the actual maximum
        # GetYaxis().GetXmax() doesnt work either

        modifier = self.get_modifier()
        ymax = max([o.GetMaximum() for o in self.contributions_objs])
        if self.main_pad.GetLogy():
            modifier.SetMaximum(ymax * self.y_padding_max_log)
        else:
            modifier.SetMaximum(ymax * self.y_padding_max_linear)

        # this is tricky... how to handle various cases like -ve, ignoring 0s
        if isinstance(self.contributions_objs[0], ROOT.TH1):
            ymin = min([o.GetMinimum() for o in self.contributions_objs])
            if ymin >= 0:
                ymin = min([o.GetMinimum(1E-20) for o in self.contributions_objs])
        else:
            # TGraph doesn't have the argument
            ymin = min([o.GetMinimum() for o in self.contributions_objs])

        if self.main_pad.GetLogy():
            if ymin > 0:
                modifier.SetMinimum(ymin * self.y_padding_min_log)
        else:
            if ymin < 0:
                modifier.SetMinimum(ymin * self.y_padding_min_linear)
            elif ymin > 0:
                modifier.SetMinimum(ymin / self.y_padding_min_linear)

    def plot(self, draw_opts=None, subplot_draw_opts=None):
        """Make the plot.

        draw_opts: str
            Same options as you would pass to Draw() in ROOT, for the main plot
        subplot_draw_opts: str
            Same options as you would pass to Draw() in ROOT, for the subplot
        """
        # Now add all the contributions to the container, styling as we go
        if len(self.contributions) == 0:
            raise ZeroContributions("Contributions list is empty")

        if self.plot_what == 'hist':
            # check hists have data in them
            has_entries = [c.obj.GetEntries() > 0 for c in self.contributions]
            if not any(has_entries):
                raise ZeroContributions("All contributions have 0 entries")
        elif self.plot_what == 'graph':
            # check graphs have data in them
            has_entries = [c.obj.GetN() > 0 for c in self.contributions]
            if not any(has_entries):
                raise ZeroContributions("All contributions have 0 entries")

        if not self.container:
            # First make a container
            self._create_container()

            # Now add all the contributions to the container, styling as we go
            self._populate_container_and_legend()

        # Set default drawing opts
        # need different default drawing options for TF1s vs TGraphs
        # (ALP won't actually display axis for TF1)
        if draw_opts is None:
            if self.plot_what in ["graph", 'both']:
                draw_opts = "ALP"
            elif self.plot_what in ['function']:
                draw_opts = ""
            elif self.plot_what == "hist":
                draw_opts = "NOSTACK"

        if "STACK" in draw_opts and not "NOSTACK" in draw_opts:
            draw_opts = draw_opts.replace("STACK", "")
            print("WARNING: 'stack' is not a valid draw option - the default is "
                  "stacking, so will draw with options '%s'" % (draw_opts))

        if not subplot_draw_opts:
            subplot_draw_opts = draw_opts

        # Need a canvas before drawing
        # If we have "SAME" then we want to draw this Plot object
        # on the current canvas
        # Otherwise, we create a new one
        if not self.canvas:
            if "SAME" in draw_opts:
                self.canvas = ROOT.gPad
                self.canvas.cd()
                # print("Using existing canvas", self.canvas.GetName())
            else:
                self.canvas = ROOT.TCanvas(ROOT.TUUID().AsString(), "", *self.default_canvas_size)
                self.canvas.SetTicks(1, 1)

                if self.subplot_type:
                    self.main_pad = ROOT.TPad("main_pad", "", 0, self.subplot_pad_height+self.subplot_pad_fudge, 1, 1)
                    ROOT.SetOwnership(self.main_pad, False)
                    self.main_pad.SetTicks(1, 1)
                    self.main_pad.SetBottomMargin(3*self.subplot_pad_fudge)
                    self.top_margin /= (1-self.subplot_pad_height)
                    self.main_pad.SetTopMargin(self.top_margin)
                    self.right_margin /= (1-self.subplot_pad_height)
                    self.main_pad.SetRightMargin(self.right_margin)
                    self.left_margin /= (1-self.subplot_pad_height)
                    self.main_pad.SetLeftMargin(self.left_margin)
                    self.canvas.cd()
                    self.main_pad.Draw()
                    self.subplot_pad = ROOT.TPad("subplot_pad", "", 0, 0, 1, self.subplot_pad_height-self.subplot_pad_fudge)
                    ROOT.SetOwnership(self.subplot_pad, False)
                    self.subplot_pad.SetTicks(1, 1)
                    self.subplot_pad.SetFillColor(0)
                    self.subplot_pad.SetFillStyle(0)
                    self.subplot_pad.SetTopMargin(3*self.subplot_pad_fudge)
                    self.subplot_pad.SetRightMargin(self.right_margin)
                    self.subplot_pad.SetBottomMargin(0.32)
                    self.subplot_pad.SetLeftMargin(self.left_margin)
                    self.canvas.cd()
                    self.subplot_pad.Draw()
                else:
                    self.main_pad = ROOT.TPad("main_pad", "", 0, 0, 1, 1)
                    ROOT.SetOwnership(self.main_pad, False)
                    self.main_pad.SetTopMargin(self.top_margin)
                    self.main_pad.SetRightMargin(self.right_margin)
                    self.main_pad.SetLeftMargin(self.left_margin)
                    self.main_pad.SetTicks(1, 1)
                    self.main_pad.Draw()

        self.canvas.cd()
        self.main_pad.cd()

        # Plot the container
        self.container.Draw(draw_opts)

        # Customise
        modifier = self.get_modifier()

        # Use the x/y axis labels from the first contribution as defaults
        if not self.xtitle:
            self.xtitle = self.contributions_objs[0].GetXaxis().GetTitle()
        if not self.ytitle:
            self.ytitle = self.contributions_objs[0].GetYaxis().GetTitle()

        modifier.SetTitle(";%s;%s" % (self.xtitle, self.ytitle))

        for ob in self.contributions_objs:
            # who knows which actually sets the thing?
            ob.GetXaxis().SetTitle(self.xtitle)
            ob.GetYaxis().SetTitle(self.ytitle)

        if self.xlim and all([x is not None for x in self.xlim]):
            if self.plot_what == "graph":
                modifier.GetXaxis().SetLimits(*self.xlim)
            else:
                modifier.GetXaxis().SetRangeUser(*self.xlim)
        if self.ylim and all([x is not None for x in self.ylim]):
            # dont use the SetLimits for graphs, that doesnt work properly.
            # no idea, ROOT is fucking insane
            modifier.GetYaxis().SetRangeUser(*self.ylim)
            modifier.SetMinimum(self.ylim[0])  # use this, not SetRangeUser()
            modifier.SetMaximum(self.ylim[1])
        else:
            # add some padding to the y limits
            self.canvas.Update()
            self._set_automatic_y_limits()

        # Draw it again to update
        if self.plot_what == "graph":
            self.container.Draw(draw_opts)

        # Plot legend
        if self.do_legend:
            self._style_legend()
            self.canvas.cd()
            self.legend.Draw()
            self.main_pad.cd()


        # Add CMS text
        self.canvas.cd()
        cms_latex = ROOT.TLatex()
        cms_latex.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        cms_latex.SetTextFont(42)
        cms_latex.SetTextSize(self.cms_text_font_size)
        latex_height = 0.915
        # left_offset = (self.left_margin - 0.08)  # to account for left margin, magic numbers ahoy

        start_x = self.left_margin + self.text_left_offset
        if self.is_preliminary:
            if self.has_data:
                cms_latex.DrawLatex(start_x, latex_height, "#font[62]{CMS}#font[52]{ Preliminary}")
            else:
                cms_latex.DrawLatex(start_x, latex_height, "#font[62]{CMS}#font[52]{ Preliminary Simulation}")
        else:
            if self.has_data:
                cms_latex.DrawLatex(start_x, latex_height, "#font[62]{CMS}")
            else:
                cms_latex.DrawLatex(start_x, latex_height, "#font[62]{CMS}#font[52]{ Simulation}")
        # cms_latex.DrawLatex(0.14, latex_height, "#font[62]{CMS}#font[52]{ Preliminary}")
        # cms_latex.DrawLatex(0.14, latex_height, "#font[62]{CMS}")
        cms_latex.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignBottom)
        cms_latex.DrawLatex(0.94, latex_height, " 35.9 fb^{-1} (13 TeV)")

        # Add title to plot
        text_latex = ROOT.TLatex()
        text_latex.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignTop)
        text_latex.SetTextFont(42)
        text_latex.SetTextSize(self.title_font_size)
        for ind, line in enumerate(self.title.split('\n')):
            text_latex.DrawLatex(self.left_margin + self.title_left_offset + self.text_left_offset, self.title_start_y - (ind*self.title_diff_y), line)

        # Do subplot
        if self.subplot_type:
            self._rescale_plot_labels(modifier, 1-self.subplot_pad_height)

            # Get rid of main plot x axis labels
            modifier.GetHistogram().GetXaxis().SetLabelSize(0)
            modifier.GetXaxis().SetLabelSize(0)
            modifier.GetHistogram().GetXaxis().SetTitleSize(0)
            modifier.GetXaxis().SetTitleSize(0)

            self.subplot_pad.cd()
            self.subplot_container.Draw(subplot_draw_opts)

            if self.subplot_title == None:
                if (self.subplot_type == "ratio"):
                    self.subplot_title = "#splitline{Ratio vs}{%s}" % (self.subplot.label)
                elif (self.subplot_type == "diff"):
                    self.subplot_title = "#splitline{Difference}{vs %s}" % (self.subplot.label)
                elif (self.subplot_type == "ddelta"):
                    self.subplot_title = "d#Delta/d#lambda"

            self.subplot_container.SetTitle(";%s;%s" % (self.xtitle, self.subplot_title))


            if self.xlim:
                self.subplot_container.GetXaxis().SetRangeUser(*self.xlim)

            if self.subplot_type == "ratio":
                # Set ratio y limits
                self.subplot_container.SetMinimum(0)  # use this, not SetRangeUser()
                if self.subplot_limits:
                    # User-specified limits
                    self.subplot_container.SetMinimum(self.subplot_limits[0])  # use this, not SetRangeUser()
                    self.subplot_container.SetMaximum(self.subplot_limits[1])  # use this, not SetRangeUser()
                else:
                    # A more intelligent way?
                    # Make sure that the upper limit is the largest bin of the contributions,
                    # so long as it is within 1.5 and some upper limit
                    harrays = [cu.th1_to_arr(h)[0] for h in self.subplot_contributions]
                    try:
                        bin_maxs = [np.max(arr[np.nonzero(arr)]) for arr in harrays]
                    except ValueError:
                        bin_maxs = [0]
                    if len(bin_maxs) == 0:
                        bin_maxs = [0]
                    self.subplot_container.SetMaximum(min(2.5, max(1.5, 1.*max(bin_maxs))))

                    # Make sure the lower limit is the smallest bin of the contributions,
                    # so long as it is within 0 and 0.5
                    try:
                        bin_mins = [np.min(arr[np.nonzero(arr)]) for arr in harrays]
                    except ValueError:
                        bin_mins = [0]
                    if len(bin_mins) == 0:
                        bin_mins = [0]
                    self.subplot_container.SetMinimum(min(0.5, 0.75*min(bin_mins)))

                # Draw a line at 1
                xax = modifier.GetXaxis()
                ax_min, ax_max = xax.GetXmin(), xax.GetXmax()
                self.subplot_line = ROOT.TLine(xax.GetXmin(), 1., xax.GetXmax(), 1.)
                if self.xlim and all([x is not None for x in self.xlim]):
                    # use user-set limits but ensure within in actual limits
                    self.subplot_line = ROOT.TLine(max(ax_min, self.xlim[0]), 1., min(ax_max, self.xlim[1]), 1.)
                self.subplot_line.SetLineStyle(self.subplot_line_style)
                self.subplot_line.SetLineWidth(self.subplot_line_width)
                self.subplot_line.SetLineColor(self.subplot_line_color)
                self.subplot_line.Draw()

            # Some resizing of subplot things
            self._rescale_plot_labels(self.subplot_container, self.subplot_pad_height)
            self.subplot_container.GetXaxis().SetTitleOffset(self.subplot_container.GetXaxis().GetTitleOffset()*3)
            self.subplot_container.GetYaxis().SetNdivisions(505)

            # If subplot y label is 2 lines, let's scale it down a bit
            if "splitline" in self.subplot_title:
                self.subplot_container.GetYaxis().SetTitleSize(self.subplot_container.GetYaxis().GetTitleSize()*0.8)
                self.subplot_container.GetYaxis().SetTitleOffset(self.subplot_container.GetYaxis().GetTitleOffset()*1.3)

            self.subplot_pad.Update()
            self.canvas.Update()


    def save(self, filename):
        """Save the plot to file. Do some check to make sure dir exists."""
        filename = os.path.abspath(filename)
        if not os.path.isdir(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        self.canvas.SaveAs(filename)
        # del self.canvas
