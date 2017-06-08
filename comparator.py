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
import ROOT
import common_utils as cu
from uuid import uuid4
from MyStyle import My_Style


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
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

        print x_min, x_max
        print y_min, y_max

        self.funcs[0].GetXaxis().SetRangeUser(x_min * 0.95, x_max * 1.05)
        self.funcs[0].GetYaxis().SetRangeUser(y_min * 0.95, y_max * 1.05)

    def Mod(self):
        """If we want to do any sort of styling, need to access the first drawn func"""
        return self.funcs[0]


def grab_obj(file_name, obj_name):
    """Get object names obj_name from ROOT file file_name"""
    # TODO: checks!
    input_file = cu.open_root_file(file_name)
    obj = cu.get_from_file(input_file, obj_name)
    obj.SetDirectory(0)  # Ownership kludge
    input_file.Close()
    return obj.Clone()


class Contribution(object):
    """Basic class to handle information about one contribution to a canvas."""

    def __init__(self, obj, label=None,
                 line_width=1, line_color=ROOT.kRed, line_style=1,
                 fill_color=ROOT.kRed, fill_style=1,
                 marker_size=1, marker_color=ROOT.kRed, marker_style=1,
                 normalise_hist=False, rebin_hist=None):
        """
        obj : TH1, TGraph, ...
            Object to be plotted
        label: str
            Title of contribution, to be used in legend.
        line_width: int
            Width of line.
        line_color: int or ROOT color
            Color of line.
        line_style: int
            Line style.
        fill_color: int or ROOT color
            Color of fill.
        fill_style: int
            Fill style.
        marker_size: int
            Size of marker
        marker_color: int or ROOT color
            Color of markers.
        marker_style: int
            Marker style.
        normalise_hist : bool
            If a histogram, normalise so integral = 1
        rebin_hist : int, None
            If a histogram, specify the number of bins to be grouped together.
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

        self.obj.SetLineWidth(self.line_width)
        self.obj.SetLineColor(self.line_color)
        self.obj.SetLineStyle(self.line_style)
        self.obj.SetFillColor(self.fill_color)
        self.obj.SetFillStyle(self.fill_style)
        self.obj.SetMarkerSize(self.marker_size)
        self.obj.SetMarkerColor(self.marker_color)
        self.obj.SetMarkerStyle(self.marker_style)
        if rebin_hist:
            self.obj.Rebin(rebin_hist) # Does this handle 2D hists?
        if normalise_hist:
            self.obj.Scale(1./obj.Integral())
        if isinstance(self.obj, ROOT.TH1) or isinstance(self.obj, ROOT.TH2):
            self.obj.SetDirectory(0)

    def __eq__(self, other):
        return self.obj == other.obj


class Plot(object):
    """
    Basic class to handle information about one plot,
    which can have several contributions.
    """

    def __init__(self, contributions=None, what="graph",
                 title=None, xtitle=None, ytitle=None, xlim=None, ylim=None,
                 legend=True, extend=False,
                 subplot=None, subplot_type="ratio"):
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
        """
        self.contributions = contributions if contributions else []
        self.contributions_objs = []
        options = ['hist', 'graph', 'function', 'both']
        if what not in options:
            raise RuntimeError("`what` argument must be one of %s" % options)
        self.plot_what = what
        self.title = title or ""
        self.xtitle = xtitle  # don't do or "" trick here, see later
        self.ytitle = ytitle  # don't do or "" trick here, see later
        self.xlim = xlim
        self.ylim = ylim
        self.do_legend = legend
        self.legend = ROOT.TLegend(0.65, 0.6, 0.94, 0.85) if legend else None
        self.do_extend = extend
        self.container = None
        self.canvas = None
        self.default_canvas_size = (800, 800)
        self.main_pad = None
        self.subplot = subplot
        if subplot_type not in ['ratio', 'diff', "ddelta"]:
            raise RuntimeError("subplot_type must be ratio, diff, or ddelta")
        self.subplot_type = subplot_type
        self.subplot_container = None
        self.subplot_contributions = []
        self.subplot_pad = None
        self.subplot_ratio_lim = (0, 2)
        self.subplot_pad_height = 0.32
        self.subplot_pad_fudge = 0.01  # to get non-overlapping subplot axis
        self.subplot_line = None  # need this to remain visible...

    def add_contribution(self, *contribution):
        """Add Contribution to Plot. Can be single item or list."""
        self.contributions.extend(*contribution)

    def create_container(self):
        """Create a container object for lots of the same TObject"""
        if self.plot_what in ["graph", "both"]:
            self.container = ROOT.TMultiGraph("mg%s" % uuid4(), "")
            self.subplot_container = ROOT.TMultiGraph("mg_ratio%s" % uuid4(), "")
        elif self.plot_what == "function":
            self.container = MultiFunc()
            self.subplot_container = MultiFunc()
        elif self.plot_what == "hist":
            self.container = ROOT.THStack("hst%s" % uuid4(), "")
            self.subplot_container = ROOT.THStack("hst_ratio%s" % uuid4(), "")

    def populate_container_and_legend(self):
        """Add objects to the container, and to the legend"""
        for c in self.contributions:
            try:
                obj = c.obj.Clone()
            except IOError:
                print "Couldn't get", c.obj_name, 'from', c.file_name
                continue

            if self.plot_what == "graph":
                # if drawing only graph, we need to remove the fit if there is one
                obj.GetListOfFunctions().Remove(obj.GetListOfFunctions().At(0))
            else:
                # if drawing function, extend range as per user's request
                if self.do_extend:
                    if self.plot_what == 'function':
                        obj.SetRange(self.xlim[0], self.xlim[1])
                    elif self.plot_what == "both":
                        obj.GetListOfFunctions().At(0).SetRange(self.xlim[0], self.xlim[1])

            self.container.Add(obj)
            self.contributions_objs.append(obj)

            if self.do_legend:
                opt = "LF" if self.plot_what == "hist" else "LP"
                self.legend.AddEntry(obj, c.label, opt)

            # Add contributions for the subplot
            if self.subplot:
                subplot_obj = self.subplot.obj.Clone()
                if c != self.subplot:
                    new_hist = obj.Clone()
                    if (self.subplot_type == "ratio"):
                        new_hist.Divide(subplot_obj)
                    elif (self.subplot_type == "diff"):
                        new_hist.Add(subplot_obj, -1.)
                    elif (self.subplot_type == "ddelta"):
                        # Do the differntial delta spectrum, see 1704.03878
                        new_hist.Add(subplot_obj, -1.)
                        new_hist.Multiply(new_hist)
                        sum_hist = obj.Clone()
                        sum_hist.Add(subplot_obj)
                        new_hist.Divide(sum_hist)
                        new_hist.Scale(0.5)
                    self.subplot_container.Add(new_hist)
                    self.subplot_contributions.append(new_hist)

    def style_legend(self):
        # self.legend.SetBorderSize(0)
        self.legend.SetFillStyle(0)

    def rescale_plot_labels(self, container, factor):
        # What a pile of wank, why does ROOT scale all these sizes?
        container.GetXaxis().SetLabelSize(0.03/factor)
        container.GetXaxis().SetTitleSize(0.04/factor)
        # container.GetXaxis().SetTitleOffset(1.0)  # doesn't seem to work?
        container.GetXaxis().SetTickLength(0.02/factor)

        container.GetYaxis().SetLabelSize(0.03/factor)
        container.GetYaxis().SetTitleSize(0.04/factor)
        container.GetYaxis().SetTitleOffset(1.2*factor)
        # container.GetYaxis().SetTickLength(0.03/factor)

    def plot(self, draw_opts=None):
        """Make the plot.

        draw_opts: str
            Same options as you would pass to Draw() in ROOT.
        """
        if not self.container:
            # First make a container
            self.create_container()

            # Now add all the contributions to the container, styling as we go
            if len(self.contributions) == 0:
                raise UnboundLocalError("contributions list is empty")

            self.populate_container_and_legend()

        # Set default drawing opts
        # need different default drawing options for TF1s vs TGraphs
        # (ALP won't actually display axis for TF1)
        if draw_opts is None:
            if self.plot_what in ["graph", 'both']:
                draw_opts = "ALP"
            elif self.plot_what in ['function']:
                draw_opts = ""
            elif self.plot_what == "hist":
                draw_opts = ""

        # Need a canvas before drawing
        # If we have "SAME" then we want to draw this Plot object
        # on the current canvas
        # Otherwise, we create a new one
        if not self.canvas:
            if "SAME" in draw_opts:
                self.canvas = ROOT.gPad
                self.canvas.cd()
                print "Using existing canvas", self.canvas.GetName()
            else:
                self.canvas = ROOT.TCanvas("canv%s" % uuid4(), "", *self.default_canvas_size)
                self.canvas.SetTicks(1, 1)
                right_margin = 0.03
                if self.subplot:
                    self.main_pad = ROOT.TPad("main_pad", "", 0, self.subplot_pad_height+self.subplot_pad_fudge, 1, 1)
                    self.main_pad.SetTicks(1, 1)
                    self.main_pad.SetBottomMargin(2*self.subplot_pad_fudge)
                    self.main_pad.SetTopMargin(0.12)
                    self.main_pad.SetRightMargin(right_margin)
                    self.canvas.cd()
                    self.main_pad.Draw()
                    self.subplot_pad = ROOT.TPad("subplot_pad", "", 0, 0, 1, self.subplot_pad_height-self.subplot_pad_fudge)
                    self.subplot_pad.SetTicks(1, 1)
                    self.subplot_pad.SetFillColor(0)
                    self.subplot_pad.SetFillStyle(0)
                    self.subplot_pad.SetTopMargin(0.02)
                    self.subplot_pad.SetRightMargin(right_margin)
                    self.subplot_pad.SetBottomMargin(0.3)
                    self.canvas.cd()
                    self.subplot_pad.Draw()
                else:
                    self.main_pad = ROOT.TPad("main_pad", "", 0, 0, 1, 1)
                    self.main_pad.SetRightMargin(right_margin)
                    self.main_pad.SetTicks(1, 1)
                    self.main_pad.Draw()

        self.canvas.cd()
        self.main_pad.cd()

        # Plot the container
        self.container.Draw(draw_opts)

        # Customise
        if self.plot_what != 'function':
            modifier = self.container
        else:
            modifier = self.container.Mod()

        # Use the x/y axis labels from the first contribution as defaults
        if not self.xtitle:
            self.xtitle = self.contributions_objs[0].GetXaxis().GetTitle()
        if not self.ytitle:
            self.ytitle = self.contributions_objs[0].GetYaxis().GetTitle()

        modifier.SetTitle("%s;%s;%s" % (self.title, self.xtitle, self.ytitle))

        if self.xlim:
            modifier.GetXaxis().SetRangeUser(*self.xlim)
        if self.ylim:
            modifier.GetYaxis().SetRangeUser(*self.ylim)

        # Plot legend
        self.style_legend()
        if self.do_legend:
            self.legend.Draw()

        if self.subplot:
            self.rescale_plot_labels(modifier, 1-self.subplot_pad_height)

            # Get rid of main plot x axis labels
            modifier.GetHistogram().GetXaxis().SetLabelSize(0)
            modifier.GetXaxis().SetLabelSize(0)
            modifier.GetHistogram().GetXaxis().SetTitleSize(0)
            modifier.GetXaxis().SetTitleSize(0)

            self.subplot_pad.cd()
            self.subplot_container.Draw(draw_opts)

            if (self.subplot_type == "ratio"):
                self.subplot_container.SetTitle(";%s;Ratio" % (self.xtitle))
            elif (self.subplot_type == "diff"):
                self.subplot_container.SetTitle(";%s;Difference" % (self.xtitle))
            elif (self.subplot_type == "ddelta"):
                self.subplot_container.SetTitle(";%s;d#Delta/d#lambda" % (self.xtitle))

            self.rescale_plot_labels(self.subplot_container, self.subplot_pad_height)
            self.subplot_container.GetYaxis().SetNdivisions(505)

            if self.subplot_type == "ratio":
                # self.subplot_container.SetMinimum(self.subplot_ratio_lim[0])  # use this, not SetRangeUser()
                # self.subplot_container.SetMaximum(self.subplot_ratio_lim[1])
                self.subplot_line = ROOT.TLine(0., 1., 1., 1.,)
                self.subplot_line.SetLineStyle(2)
                self.subplot_line.SetLineWidth(2)
                self.subplot_line.SetLineColor(ROOT.kBlack)
                self.subplot_line.Draw()

            self.subplot_pad.Update()
            self.canvas.Update()


    def save(self, filename):
        """Save the plot to file. Do some check to make sure dir exists."""
        filename = os.path.abspath(filename)
        if not os.path.isdir(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))
        self.canvas.SaveAs(filename)
