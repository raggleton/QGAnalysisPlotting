"""
Main class to handle plotting of MyUnfolder
"""


from __future__ import print_function, division

from array import array
import numpy as np
import math
import os
from itertools import chain

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()

import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp

# This doesn't seem to work...sigh
np.set_printoptions(edgeitems=3,infstr='Infinity',
                    linewidth=75, nanstr='nan', precision=8,
                    suppress=False, threshold=1000, formatter=None)

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat(".3f")
ROOT.gStyle.SetHistTopMargin(0.)


class MyUnfolderPlotter(object):

    """Class to do all sort of plots about a given MyUnfolder obj"""

    def __init__(self, unfolder):
        self.unfolder = unfolder
        self.output_fmt = 'pdf'

    @staticmethod
    def generate_2d_canvas(size=(800, 600)):
        canv = ROOT.TCanvas(cu.get_unique_str(), "", *size)
        canv.SetTicks(1, 1)
        canv.SetLeftMargin(0.11)
        canv.SetRightMargin(0.12)
        return canv

    # @staticmethod
    def draw_2d_hist(self, h2d, title, output_filename,
                     logz=True, z_min=None, z_max=None,
                     xtitle="Generator bin",
                     ytitle="Generator bin",
                     canvas_size=(800, 600),
                     draw_bin_lines_x=False,
                     draw_bin_lines_y=False,
                     draw_values=False):
        canv = MyUnfolderPlotter.generate_2d_canvas(canvas_size)
        if draw_bin_lines_x:
            canv.SetBottomMargin(0.18)
        if draw_bin_lines_y:
            canv.SetLeftMargin(0.18)
        if logz:
            canv.SetLogz()
        this_title = "%s;%s;%s" % (title, xtitle, ytitle)
        h2d.SetTitle(this_title)
        orig_title_font_size = ROOT.gStyle.GetTitleFontSize()
        if 'splitline' in title:
            ROOT.gStyle.SetTitleFontSize(0.025)
            ROOT.gStyle.SetTitleAlign(23)

        h2d.GetXaxis().SetTitleOffset(1.5)
        h2d.GetYaxis().SetTitleOffset(1.5)
        if draw_bin_lines_x:
            h2d.GetXaxis().SetTitleOffset(2.2)
            h2d.GetXaxis().SetLabelOffset(999)
        if draw_bin_lines_y:
            h2d.GetYaxis().SetTitleOffset(2.1)
            h2d.GetYaxis().SetLabelOffset(999)
        if draw_values:
            h2d.SetMarkerSize(0.5)
            h2d.Draw("COLZ TEXT45")
        else:
            h2d.Draw("COLZ")
        if z_max:
            h2d.SetMaximum(z_max)
        if z_min:
            h2d.SetMinimum(z_min)
        elif logz:
            h2d.SetMinimum(h2d.GetMinimum(1E-40) / 10.)

        if draw_bin_lines_x:
            lx, tx = self.draw_pt_binning_lines(h2d,
                                                which='gen' if 'gen' in xtitle.lower() else 'reco',
                                                axis='x',
                                                do_underflow=True,
                                                do_labels_inside=False,
                                                do_labels_outside=True)

        if draw_bin_lines_y:
            ly, ty = self.draw_pt_binning_lines(h2d,
                                                which='gen' if 'gen' in ytitle.lower() else 'reco',
                                                axis='y',
                                                do_underflow=True,
                                                do_labels_inside=False,
                                                do_labels_outside=True)

        odir = os.path.dirname(os.path.abspath(output_filename))
        cu.check_dir_exists_create(odir)
        canv.SaveAs(output_filename)
        if 'splitline' in title:
            ROOT.gStyle.SetTitleFontSize(orig_title_font_size)


    def plot_bias_vector(self, output_dir='.', append="", title=""):
        """Plot bias vector used in regularisation (if it exists)"""
        if append != "":
            append = "_%s" % (append)

        bias_hist = self.unfolder.get_bias_vector()
        entries = [
            Contribution(bias_hist,
                         label="Bias histogram",
                         line_color=ROOT.kBlack, line_width=1,
                         marker_color=ROOT.kBlack, marker_size=0,
                         normalise_hist=False),
        ]
        plot = Plot(entries,
                    what='hist',
                    title=title,
                    xtitle="Generator bin",
                    ytitle="Bias")
        plot.default_canvas_size = (800, 600)
        plot.plot("NOSTACK HIST")
        plot.set_logy()
        plot.legend.SetY1NDC(0.77)
        plot.legend.SetX2NDC(0.85)
        output_filename = "%s/bias_hist%s.%s" % (output_dir, append, self.output_fmt)
        plot.save(output_filename)

    def draw_response_matrix(self, output_dir='.', append="", title=""):
        output_filename = "%s/response_map_%s.%s" % (output_dir, append, self.output_fmt)
        self.draw_2d_hist(self.unfolder.response_map,
                          title=title,
                          output_filename=output_filename,
                          xtitle='Generator bin', ytitle='Detector bin',
                          draw_bin_lines_x=True,
                          draw_bin_lines_y=True,
                          canvas_size=(800, 700))

    def draw_covariance_matrix(self, output_dir='.', append="", title=""):
        output_filename = "%s/covariance_map_%s.%s" % (output_dir, append, self.output_fmt)
        self.draw_2d_hist(self.unfolder.get_covariance_matrix(),
                          title=title,
                          output_filename=output_filename,
                          draw_bin_lines_x=True,
                          draw_bin_lines_y=True,
                          canvas_size=(800, 700),
                          xtitle='Generator bin', ytitle='Generator bin')

    def draw_probability_matrix(self, output_dir='.', append="", title=""):
        output_filename = "%s/probability_map_%s.%s" % (output_dir, append, self.output_fmt)
        self.draw_2d_hist(self.unfolder.get_probability_matrix(),
                          title=title,
                          output_filename=output_filename,
                          # z_min=1E-4, z_max=1,
                          logz=False,
                          draw_bin_lines_x=True,
                          draw_bin_lines_y=True,
                          canvas_size=(800, 700),
                          xtitle='Generator bin',
                          ytitle='Detector bin')

    # TODO: generalise to some "draw_2d_hist()"?
    def draw_error_matrix_input(self, output_dir='.', append="", title=""):
        output_filename = "%s/err_map_stat_input_%s.%s" % (output_dir, append, self.output_fmt)
        self.draw_2d_hist(self.unfolder.get_ematrix_input(),
                          title=title,
                          output_filename=output_filename,
                          draw_bin_lines_x=True,
                          draw_bin_lines_y=True,
                          canvas_size=(800, 700))

    def draw_error_matrix_stat(self, output_dir='.', append="", title=""):
        output_filename = "%s/err_map_stat_%s.%s" % (output_dir, append, self.output_fmt)
        self.draw_2d_hist(self.unfolder.get_ematrix_stat(),
                          title=title,
                          output_filename=output_filename,
                          draw_bin_lines_x=True,
                          draw_bin_lines_y=True,
                          canvas_size=(800, 700))

    def draw_error_matrix_stat_response(self, output_dir='.', append="", title=""):
        output_filename = "%s/err_map_stat_response_%s.%s" % (output_dir, append, self.output_fmt)
        self.draw_2d_hist(self.unfolder.get_ematrix_stat_response(),
                          title=title,
                          output_filename=output_filename,
                          draw_bin_lines_x=True,
                          draw_bin_lines_y=True,
                          canvas_size=(800, 700))

    def draw_error_matrix_tau(self, output_dir='.', append="", title=""):
        output_filename = "%s/err_map_tau_%s.%s" % (output_dir, append, self.output_fmt)
        self.draw_2d_hist(self.unfolder.get_ematrix_tau(),
                          title=title,
                          output_filename=output_filename,
                          draw_bin_lines_x=True,
                          draw_bin_lines_y=True,
                          canvas_size=(800, 700))

    def draw_error_matrix_syst(self, syst_label, output_dir='.', append="", title=""):
        syst_label_no_spaces = syst_label.replace(" ", "_")
        output_filename = "%s/err_map_syst_%s_%s.%s" % (output_dir, syst_label_no_spaces, append, self.output_fmt)
        self.draw_2d_hist(self.unfolder.get_ematrix_syst(syst_label),
                          title=title,
                          output_filename=output_filename,
                          draw_bin_lines_x=True,
                          draw_bin_lines_y=True,
                          canvas_size=(800, 700))

    def draw_error_matrix_total(self, output_dir='.', append="", title=""):
        output_filename = "%s/err_map_total_%s.%s" % (output_dir, append, self.output_fmt)
        self.draw_2d_hist(self.unfolder.get_ematrix_total(),
                          title=title,
                          output_filename=output_filename,
                          draw_bin_lines_x=True,
                          draw_bin_lines_y=True,
                          canvas_size=(800, 700))

    def draw_correlation_matrix(self, draw_values=False, draw_binning_lines=True, output_dir='.', append="", title=""):
        # Custom colour scheme - french flag colouring
        NRGBs = 3
        NCont = 99
        stops = [ 0.00, 0.53, 1.00 ]
        red = [ 0.00, 1.00, 1.00]
        green = [ 0.00, 1.00, 0.00 ]
        blue = [ 1.00, 1.00, 0.00 ]
        stopsArray = array('d', stops)
        redArray = array('d', red)
        greenArray = array('d', green)
        blueArray = array('d', blue)
        ROOT.TColor.CreateGradientColorTable(NRGBs, stopsArray, redArray, greenArray, blueArray, NCont)

        canv = self.generate_2d_canvas((800, 700))
        canv.SetLeftMargin(0.11)
        canv.SetRightMargin(0.12)
        if draw_binning_lines:
            canv.SetBottomMargin(0.18)
            canv.SetLeftMargin(0.18)
        corr_map = self.unfolder.get_rhoij_total()
        corr_map.SetTitle(title+";Generator bin;Generator bin")
        corr_map.GetYaxis().SetTitleOffset(1)
        # corr_map.GetXaxis().SetTitleOffset(1.5)
        corr_map.SetMinimum(-1)
        corr_map.SetMaximum(1)
        corr_map.SetMarkerSize(0.05)
        val_str = ""
        if draw_values:
            val_str = "_withValues"
            corr_map.Draw("COL1Z TEXT45")
        else:
            corr_map.Draw("COL1Z")
        if draw_binning_lines:
            corr_map.GetXaxis().SetTitleOffset(2.2)
            corr_map.GetXaxis().SetLabelOffset(999)
            corr_map.GetYaxis().SetTitleOffset(2.1)
            corr_map.GetYaxis().SetLabelOffset(999)
            lx, tx = self.draw_pt_binning_lines(corr_map, which='gen', axis='x', do_labels_inside=False, do_labels_outside=True)
            ly, ty = self.draw_pt_binning_lines(corr_map, which='gen', axis='y', do_labels_inside=False, do_labels_outside=True)
        output_filename = "%s/rho_map%s_%s.%s" % (output_dir, val_str, append, self.output_fmt)
        canv.SaveAs(output_filename)
        ROOT.gStyle.SetPalette(ROOT.kBird)

    def draw_background_fractions(self, output_dir='.', append="", title=""):
        """Do plot of individual background fractions, a cumulative one, and one with all sources non-cumulative"""
        all_contributions = []
        frac_min, frac_max = 3E-4, 3
        cu.check_dir_exists_create(output_dir)
        # Do individual plots for each bg source
        for bg_name, bg_hist in self.unfolder.backgrounds.items():
            this_fraction = bg_hist.Clone()
            this_fraction.Divide(self.unfolder.input_hist)
            entries = [Contribution(this_fraction,
                        label=bg_name,
                        line_color=ROOT.kBlue, line_width=1,
                        marker_color=ROOT.kBlue, marker_size=0)]
            all_contributions.extend(entries)
            plot = Plot(entries,
                        what='hist',
                        title=title,
                        xtitle='Detector bin number',
                        ytitle='Background fraction',
                        ylim=(frac_min, frac_max))
            plot.default_canvas_size = (800, 600)
            plot.plot('NOSTACK HIST')
            plot.set_logy()
            plot.legend.SetX1(0.75)
            plot.legend.SetY1(0.75)
            l, t = self.draw_pt_binning_lines(plot, which='reco', axis='x',
                                              do_underflow=True,
                                              do_labels_inside=True,
                                              do_labels_outside=False)
            output_filename = "%s/bg_fraction_%s_%s.%s" % (output_dir, bg_name.replace(" ", "_").lower(), append, self.output_fmt)
            plot.save(output_filename)

        # Now do one with all stacked
        # sort by ascending size
        all_contributions = sorted(all_contributions, key=lambda x: x.obj.Integral(), reverse=False)
        plot = Plot(all_contributions,
                    what='hist',
                    title=title,
                    xtitle='Detector bin number',
                    ytitle='Cumulative background fraction',
                    ylim=(frac_min, frac_max))
        plot.default_canvas_size = (800, 600)
        plot.reverse_legend = True
        ROOT.gStyle.SetPalette(ROOT.kCMYK)
        # ROOT.gStyle.SetPalette(ROOT.kPastel)
        # ROOT.gStyle.SetPalette(ROOT.kVisibleSpectrum)
        plot.plot("HIST PLC PMC")
        plot.set_logy()
        plot.legend.SetNColumns(2)
        plot.legend.SetY1(0.75)
        plot.legend.SetY2(0.88)
        l, t = self.draw_pt_binning_lines(plot, which='reco', axis='x',
                                          do_underflow=True,
                                          do_labels_inside=True,
                                          do_labels_outside=False)
        output_filename = "%s/bg_fraction_all_stack_%s.%s" % (output_dir, append, self.output_fmt)
        plot.save(output_filename)

        # And one not stacked
        plot.ytitle = "Individual background fraction"
        plot.plot("NOSTACK HIST PLC")
        l, t = self.draw_pt_binning_lines(plot, which='reco', axis='x',
                                          do_underflow=True,
                                          do_labels_inside=True,
                                          do_labels_outside=False)
        output_filename = "%s/bg_fraction_all_nostack_%s.%s" % (output_dir, append, self.output_fmt)
        plot.save(output_filename)
        ROOT.gStyle.SetPalette(ROOT.kBird)

    def draw_L_matrix(self, output_dir='.', append="", title=""):
        """Draw L matrix used for regularisation"""
        Lmatrix = self.unfolder.tunfolder.GetL("hist_Lmatrix_%s" % (append), title)
        # Custom colour scheme - french flag colouring
        NRGBs = 5
        NCont = 256
        # Set max & min such that 0 is in the middle
        stops = [ 0.00, 0.49, 0.5, 0.51, 1.00 ]
        red = [ 0.00, 1.00, 1.00, 1.00, 1.00]
        green = [ 0.00, 1.00, 1.00, 1.00, 0.00 ]
        blue = [ 1.00, 1.00, 1.00, 1.00, 0.00 ]
        stopsArray = array('d', stops)
        redArray = array('d', red)
        greenArray = array('d', green)
        blueArray = array('d', blue)
        ROOT.TColor.CreateGradientColorTable(NRGBs, stopsArray, redArray, greenArray, blueArray, NCont)

        canv = MyUnfolderPlotter.generate_2d_canvas()
        this_title = "%s;%s;%s" % (title, "Generator bin", "n_{R} bin")
        Lmatrix.SetTitle(this_title)
        h_lim = max(abs(Lmatrix.GetMinimum()), Lmatrix.GetMaximum())
        Lmatrix.SetMaximum(-h_lim)
        Lmatrix.SetMaximum(h_lim)
        Lmatrix.GetYaxis().SetTitleOffset(1.5)
        Lmatrix.GetXaxis().SetTitleOffset(1.5)
        Lmatrix.Draw("COL1 Z CJUST")
        l, t = self.draw_pt_binning_lines(Lmatrix, which='gen', axis='x',
                                          do_underflow=True,
                                          do_labels_inside=True,
                                          labels_inside_align='higher',
                                          do_labels_outside=False)
        output_filename = "%s/L_matrix_%s.%s" % (output_dir, append, self.output_fmt)
        cu.check_dir_exists_create(output_dir)
        canv.SaveAs(output_filename)
        ROOT.gStyle.SetPalette(ROOT.kBird)

    def draw_L_matrix_squared(self, output_dir='.', append="", title=""):
        """Draw L^T L, used for regularisation"""
        output_filename = "%s/L_matrix_squared_%s.%s" % (output_dir, append, self.output_fmt)
        Lmatrix = ROOT.TUnfoldBinning.CreateHistogramOfMigrations(
            self.unfolder.generator_binning,
            self.unfolder.generator_binning,
            "hist_Lmatrix_squared_%s" % (append),
            True,
            True,
            title)
        self.unfolder.tunfolder.GetLsquared(Lmatrix)
        print("Lmatrix squared:", Lmatrix.GetNbinsX())
        print("Lmatrix squared:", Lmatrix.GetNbinsY())
        self.draw_2d_hist(Lmatrix, title=title, output_filename=output_filename,
                          draw_values=False,
                          xtitle="Generator bin", ytitle="Generator bin")

    def draw_Lx_minus_bias(self, output_dir='.', append="", title=""):
        """Draw L * (x - bias_scale*bias_vec)

        NB can only do after doing ScanTau or ScanLCurve, otherwise crashes
        """
        output_filename = "%s/Lx_minus_bias_%s.%s" % (output_dir, append, self.output_fmt)
        hist_Lx_minus_bias = self.unfolder.tunfolder.GetLxMinusBias("hist_Lx_minus_bias_%s" % (append), title)
        conts = [Contribution(hist_Lx_minus_bias)]
        plot = Plot(conts,
                    what='hist',
                    title=title,
                    xtitle="n_{R} bin",
                    ytitle="L(x-f_{b}x_{0})"
                    )
        plot.y_padding_min_linear = 0.8
        plot.left_margin = 0.15
        plot.plot("HISTE")
        plot.save(output_filename)

    def draw_pt_binning_lines(self, plot, which, axis, do_underflow=True, do_labels_inside=True, do_labels_outside=False, labels_inside_align='lower'):
        """Draw lines marking pt bins

        You MUST store the return lists of lines and text objects,
        otherwise they will not appear on the saved plot - python GC cleans them up

        Parameters
        ----------
        plot : Plot or TH2
            Plot object to add lines/text to
        which : str
            "gen" or "reco"
        axis : str
            "x" or "y"
        do_underflow : bool, optional
            Include underflow bins
        do_labels_inside : bool, optional
            Add pT bin labels outside the main plot area
        do_labels_outside : bool, optional
            Add pT bin labels outside the main plot area
        labels_inside_align : str, optional
            'lower' or 'higher' ie against bottom/left, or top/right

        Raises
        ------
        ArgumentError
            Description

        Returns
        -------
        list, list
            List of TLines, TTexts
        """
        which = which.lower()
        if which not in ['gen', 'reco']:
            raise ArgumentError("'which' should be 'gen' or 'reco'")
        axis = axis.lower()
        if axis not in ['x', 'y']:
            raise ArgumentError("'axis' should be 'x' or 'y'")
        labels_inside_align = labels_inside_align.lower()
        if labels_inside_align not in ['lower', 'higher']:
            raise ArgumentError("'labels_inside_align' should be 'lower' or 'higher'")

        # setup bins, etc
        signal_pt_bins = []
        underflow_pt_bins = []
        variable_bins = []
        this_binning = None
        dist_name = None

        if which == 'gen':
            signal_pt_bins = self.unfolder.pt_bin_edges_gen
            underflow_pt_bins = self.unfolder.pt_bin_edges_underflow_gen
            variable_bins = self.unfolder.variable_bin_edges_gen
            this_binning = self.unfolder.generator_binning
            dist_name = "generatordistribution"
        else:
            signal_pt_bins = self.unfolder.pt_bin_edges_reco
            underflow_pt_bins = self.unfolder.pt_bin_edges_underflow_reco
            variable_bins = self.unfolder.variable_bin_edges_reco
            this_binning = self.unfolder.detector_binning
            dist_name = "detectordistribution"

        all_pt_bins = []
        if do_underflow:
            all_pt_bins.extend(underflow_pt_bins[:-1])
        all_pt_bins.extend(signal_pt_bins)

        axis_min, axis_max = 0, 0  # to determine start, end of line
        if isinstance(plot, Plot):
            obj = plot.container.GetHistogram()
            if axis == 'x':
                axis_low, axis_high = obj.GetMinimum(), obj.GetMaximum()
            else:
                axis_low, axis_high = obj.GetXaxis().GetXmin(), obj.GetXaxis().GetXmax()
        else:
            if axis == 'x':
                axis_low, axis_high = plot.GetYaxis().GetXmin(), plot.GetYaxis().GetXmax()
            else:
                axis_low, axis_high = plot.GetXaxis().GetXmin(), plot.GetXaxis().GetXmax()

        lines = []  # otherwise python kills them
        texts = []
        first_var = variable_bins[0]
        last_var = variable_bins[-1]
        if isinstance(plot, Plot):
            plot.main_pad.cd()
        # add line + text for each pt bin
        for pt_val, pt_val_upper in zip(all_pt_bins[:-1], all_pt_bins[1:]):
            # convert value to bin number
            # what a load of rubbish, why cant I just ask generator_binning for it?!
            binning = this_binning.FindNode("%s_underflow" % dist_name) if pt_val < signal_pt_bins[0] else this_binning.FindNode(dist_name)
            pt_bin = binning.GetGlobalBinNumber(first_var+0.000001, pt_val+0.01) - 0.5 # -0.5 for offset, since the bins are centered on the bin number (e.g. bin 81 goes from 80.5 to 81.5)

            if pt_val > all_pt_bins[0]:  # skip first value, as it probably aligns with lower axis limit
                if axis == 'x':
                    line = ROOT.TLine(pt_bin, axis_low, pt_bin, axis_high)
                else:
                    line = ROOT.TLine(axis_low, pt_bin, axis_high, pt_bin)
                line.SetLineStyle(2)
                line.SetLineColor(14)
                line.Draw()
                lines.append(line)

            # draw text for pt bins
            if do_labels_inside and axis == 'x':
                pt_bin_higher = binning.GetGlobalBinNumber(last_var-0.00001, pt_val+0.01) - 0.5
                pt_bin_interval = pt_bin_higher - pt_bin
                text_x = pt_bin + 0.3*(pt_bin_interval)
                axis_range = axis_high - axis_low
                log_axis_range = math.log10(axis_high) - math.log10(axis_low)
                # assumes logarithmic?
                if ROOT.gPad.GetLogy():
                    y_start = math.pow(10, 0.03*log_axis_range + math.log10(axis_low))
                    y_end = 10*axis_low
                else:
                    y_start = (0.02*axis_range) + axis_low
                    y_end = 10*axis_low

                if labels_inside_align == 'higher':
                    if ROOT.gPad.GetLogy():
                        y_start = axis_high/10
                        y_end = axis_high/1.05
                    else:
                        y_start = axis_high*0.8
                        y_end = axis_high*0.9

                text = ROOT.TPaveText(text_x, y_start, text_x + .5*pt_bin_interval, y_start)
                text.SetBorderSize(0)
                text.SetFillStyle(0)
                contents = '%d < p_{T} < %d GeV' % (pt_val, pt_val_upper)
                if pt_val_upper in underflow_pt_bins:
                    contents += " (underflow)"
                # You have to style the thing that is returned by AddText, not the TPaveText obj itself
                # WTAF
                t = text.AddText(contents)  # t is a TText
                t.SetTextColor(14)
                t.SetTextAngle(89)
                t.SetTextSize(0.0275)
                if (isinstance(plot, Plot) and not plot.subplot_pad) or not isinstance(plot, Plot):
                    t.SetTextSize(0.02)  # account for the fact that a subplot makes things smaller
                if labels_inside_align == 'lower':
                    t.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignTop)
                else:
                    t.SetTextAlign(ROOT.kHAlignCenter + ROOT.kVAlignTop)
                text.Draw()
                texts.append(text)

            if do_labels_outside and axis == 'x':
                pt_bin_higher = binning.GetGlobalBinNumber(last_var-0.00001, pt_val+0.01) - 0.5
                pt_bin_interval = pt_bin_higher - pt_bin
                text_x = pt_bin + 0.15*(pt_bin_interval)
                text = ROOT.TPaveText(text_x, 0.5*axis_low, text_x + .35*pt_bin_interval, 0.55*axis_low)  # urgh at some point it jsut ignores this and puts it against the axis
                text.SetBorderSize(0)
                text.SetFillStyle(0)
                contents = '[%d, %d] GeV   ' % (pt_val, pt_val_upper)  # spaces for alignment, since the text gets stuck against the axis
                if pt_val_upper in underflow_pt_bins:
                    contents = "#splitline{%s}{(underflow)   }" % (contents)
                # You have to style the thing that is returned by AddText, not the TPaveText obj itself
                # WTAF
                t = text.AddText(contents)  # t is a TText
                t.SetTextColor(14)
                t.SetTextAngle(89)
                t.SetTextSize(0.015)
                t.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignCenter)
                text.Draw()
                texts.append(text)

            if do_labels_outside and axis == 'y':
                pt_bin_higher = binning.GetGlobalBinNumber(last_var-0.00001, pt_val+0.01) - 0.5
                pt_bin_interval = pt_bin_higher - pt_bin
                text_y = pt_bin + 0.15*(pt_bin_interval)
                text = ROOT.TPaveText(0.5*axis_low, text_y, 0.5*axis_low, text_y + .35*pt_bin_interval)  # urgh at some point it jsut ignores this and puts it against the axis
                text.SetBorderSize(0)
                text.SetFillStyle(0)
                contents = '[%d, %d] GeV   ' % (pt_val, pt_val_upper)  # spaces for alignment, since the text gets stuck against the axis
                if pt_val_upper in underflow_pt_bins:
                    contents = contents.strip()
                    contents += " (underflow)   "
                # You have to style the thing that is returned by AddText, not the TPaveText obj itself
                # WTAF
                t = text.AddText(contents)  # t is a TText
                t.SetTextColor(14)
                # t.SetTextAngle(89)
                t.SetTextSize(0.015)
                t.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignCenter)
                text.Draw()
                texts.append(text)

        # lines on subplot if available
        if isinstance(plot, Plot) and plot.subplot_pad and axis == 'x':
                plot.subplot_pad.cd()
                y_low, y_high = plot.subplot_container.GetHistogram().GetMinimum(), plot.subplot_container.GetHistogram().GetMaximum()  # BINGO
                for pt_val in chain(self.unfolder.pt_bin_edges_underflow_gen[1:-1], self.unfolder.pt_bin_edges_gen[:-1]):
                    # convert value to bin number
                    # what a load of rubbish, why cant I just ask generator_binning for it?!
                    binning = this_binning.FindNode("%s_underflow" % dist_name) if pt_val < signal_pt_bins[0] else this_binning.FindNode(dist_name)
                    pt_bin = binning.GetGlobalBinNumber(first_var+0.000001, pt_val+0.01) - 0.5 # -0.5 for offset, since the bins are centered on the bin number (e.g. bin 81 goes from 80.5 to 81.5)
                    line = ROOT.TLine(pt_bin, y_low, pt_bin, y_high)
                    line.SetLineStyle(2)
                    line.SetLineColor(14)
                    line.Draw()
                    lines.append(line)

        return lines, texts

    def draw_unfolded_1d(self, do_gen=True, do_unfolded=True, is_data=True, output_dir='.', append='', title=''):
        """Simple plot of unfolded & gen, by bin number (ie non physical axes)"""
        entries = []

        if do_unfolded and self.unfolder.unfolded:
            subplot = self.unfolder.hist_truth if (do_gen and self.unfolder.hist_truth) else None
            label = 'data' if is_data else 'MC'
            entries.append(
                Contribution(self.unfolder.unfolded, label="Unfolded %s (#tau = %.3g)" % (label, self.unfolder.tau),
                             line_color=ROOT.kRed, line_width=1,
                             marker_color=ROOT.kRed, marker_size=0.6, marker_style=20,
                             subplot_line_color=ROOT.kRed, subplot_line_width=1,
                             subplot_marker_color=ROOT.kRed, subplot_marker_size=0, subplot_marker_style=20,
                             normalise_hist=False, subplot=subplot),
            )

        if do_gen and self.unfolder.hist_truth:
            entries.append(
                Contribution(self.unfolder.hist_truth, label="Generator",
                             line_color=ROOT.kBlue, line_width=1,
                             marker_color=ROOT.kBlue, marker_size=0,
                             normalise_hist=False),
            )

        plot = Plot(entries,
                    what='hist',
                    title=title,
                    xtitle="Generator bin number",
                    ytitle="N",
                    subplot_type='ratio',
                    subplot_title='#splitline{Unfolded %s /}{MC Gen}' % label,
                    subplot_limits=(0, 2),
                    has_data=is_data)
        plot.default_canvas_size = (800, 600)
        # plot.text_left_offset = 0.05  # have to bodge this
        plot.plot("NOSTACK HISTE", "NOSTACK HISTE")
        # plot.main_pad.SetLogy(1)
        plot.set_logy(do_more_labels=False)
        ymax = max([o.GetMaximum() for o in plot.contributions_objs])
        plot.container.SetMaximum(ymax * 500)
        ymin = min(c.obj.GetMinimum(1E-8) for c in entries)
        plot.container.SetMinimum(ymin*0.01)  # space for bin labels as well
        plot.legend.SetY1NDC(0.77)
        plot.legend.SetX1NDC(0.65)
        plot.legend.SetX2NDC(0.88)
        # draw pt bin lines
        plot.main_pad.cd()
        lines, text = self.draw_pt_binning_lines(plot, which='gen', axis='x', do_underflow=True)
        output_filename = "%s/unfolded_%s.%s" % (output_dir, append, self.output_fmt)
        plot.save(output_filename)

    def draw_detector_1d(self,
                        do_reco_data=False,
                        do_reco_data_bg_sub=False,
                        do_reco_bg=False,
                        do_reco_mc=False,
                        do_reco_mc_bg_sub=False,
                        output_dir='.', append="", title=""):
        """Plot detector-binned quantities for data & MC, by bin number (ie non physical axes)"""
        entries = []

        reco_bg = self.unfolder.get_total_background()
        if do_reco_bg and reco_bg:
            entries.append(
                Contribution(reco_bg, label="Background [detector-level]",
                             line_color=ROOT.kMagenta+2, line_width=1,
                             marker_color=ROOT.kMagenta+2, marker_size=0,
                             normalise_hist=False),
            )

        reco_data = self.unfolder.input_hist
        reco_mc = self.unfolder.hist_mc_reco
        if do_reco_data and reco_data:
            entries.append(
                Contribution(reco_data, label="Data [detector-level]",
                             line_color=ROOT.kGray+2, line_width=0,
                             marker_color=ROOT.kGray+2, marker_size=0.6, marker_style=20,
                             normalise_hist=False, subplot=reco_mc if do_reco_mc else None,
                             subplot_line_width=1),
            )

        reco_data_bg_sub = self.unfolder.input_hist_bg_subtracted
        reco_mc_bg_sub = self.unfolder.hist_mc_reco_bg_subtracted
        if do_reco_data_bg_sub and reco_data_bg_sub:
            entries.append(
                Contribution(reco_data_bg_sub, label="Data bg-subtracted [detector-level]",
                             line_color=ROOT.kBlack, line_width=0,
                             marker_color=ROOT.kBlack, marker_size=0.6, marker_style=20,
                             normalise_hist=False, subplot=reco_mc_bg_sub if do_reco_mc_bg_sub else None,
                             subplot_line_width=1),
            )

        if do_reco_mc and reco_mc:
            entries.append(
                Contribution(reco_mc, label="MC [detector-level]",
                             line_color=ROOT.kAzure+2, line_width=1,
                             marker_color=ROOT.kAzure+2, marker_size=0,
                             normalise_hist=False),
            )

        if do_reco_mc_bg_sub and reco_mc_bg_sub:
            entries.append(
                Contribution(reco_mc_bg_sub, label="MC bg-subtracted [detector-level]",
                             line_color=ROOT.kAzure+4, line_width=1,
                             marker_color=ROOT.kAzure+4, marker_size=0,
                             normalise_hist=False),
            )

        plot = Plot(entries,
                    what='hist',
                    title=title,
                    xtitle="Detector bin number",
                    ytitle="N",
                    subplot_type='ratio' if ((do_reco_mc and do_reco_data) or (do_reco_mc_bg_sub and do_reco_data_bg_sub)) else None,
                    subplot_title='Data / MC',
                    subplot_limits=(0.8, 1.2),
                    has_data=(do_reco_mc or do_reco_data_bg_sub))
        plot.default_canvas_size = (800, 600)
        plot.plot("NOSTACK HISTE")
        plot.set_logy(do_more_labels=False)
        ymax = max([o.GetMaximum() for o in plot.contributions_objs])
        plot.container.SetMaximum(ymax * 200)
        ymin = max([o.GetMinimum(1E-10) for o in plot.contributions_objs])
        plot.container.SetMinimum(ymin*0.01)
        l, t = self.draw_pt_binning_lines(plot, which='reco', axis='x',
                                          do_underflow=True,
                                          do_labels_inside=True,
                                          do_labels_outside=False,
                                          labels_inside_align='lower'
                                          )
        # # plot.container.SetMinimum(0.001)
        plot.legend.SetY1NDC(0.77)
        plot.legend.SetX1NDC(0.65)
        plot.legend.SetX2NDC(0.88)
        if append != "":
            append = "_" + append
        output_filename = "%s/detector_reco_binning%s.%s" % (output_dir, append, self.output_fmt)
        plot.save(output_filename)

    def draw_generator_1d(self,
                          do_reco_data=False,
                          do_reco_data_bg_sub=False,
                          do_reco_bg=False,
                          do_reco_mc=False,
                          do_reco_mc_bg_sub=False,
                          do_truth_mc=False,
                          output_dir='.', append="", title=""):
        """Plot generator-binned quantities for data & MC, by bin number (ie non physical axes)"""
        entries = []

        truth_mc = self.unfolder.hist_truth
        if do_truth_mc and truth_mc:
            entries.append(
                Contribution(truth_mc, label="MC truth [generator-level]",
                             line_color=ROOT.kRed, line_width=1,
                             marker_color=ROOT.kRed, marker_size=0,
                             normalise_hist=False),
            )

        reco_bg = self.unfolder.get_total_background_gen_binning()
        if do_reco_bg and reco_bg:
            entries.append(
                Contribution(reco_bg, label="Background [detector-level]",
                             line_color=ROOT.kMagenta+2, line_width=1,
                             marker_color=ROOT.kMagenta+2, marker_size=0,
                             normalise_hist=False),
            )

        reco_data = self.unfolder.input_hist_gen_binning
        reco_mc = self.unfolder.hist_mc_reco_gen_binning
        if do_reco_data and reco_data:
            entries.append(
                Contribution(reco_data, label="Data [detector-level]",
                             line_color=ROOT.kGray+2, line_width=0,
                             marker_color=ROOT.kGray+2, marker_size=0.6, marker_style=20,
                             normalise_hist=False, subplot=reco_mc if do_reco_mc else None,
                             subplot_line_width=1),
            )

        reco_data_bg_sub = self.unfolder.input_hist_gen_binning_bg_subtracted
        reco_mc_bg_sub = self.unfolder.hist_mc_reco_gen_binning_bg_subtracted
        if do_reco_data_bg_sub and reco_data_bg_sub:
            entries.append(
                Contribution(reco_data_bg_sub, label="Data bg-subtracted [detector-level]",
                             line_color=ROOT.kBlack, line_width=0,
                             marker_color=ROOT.kBlack, marker_size=0.6, marker_style=20,
                             normalise_hist=False, subplot=reco_mc_bg_sub if do_reco_mc_bg_sub else None,
                             subplot_line_width=1),
            )

        if do_reco_mc and reco_mc:
            entries.append(
                Contribution(reco_mc, label="MC [detector-level]",
                             line_color=ROOT.kAzure+2, line_width=1,
                             marker_color=ROOT.kAzure+2, marker_size=0,
                             normalise_hist=False),
            )

        if do_reco_mc_bg_sub and reco_mc_bg_sub:
            entries.append(
                Contribution(reco_mc_bg_sub, label="MC bg-subtracted [detector-level]",
                             line_color=ROOT.kAzure+4, line_width=1,
                             marker_color=ROOT.kAzure+4, marker_size=0,
                             normalise_hist=False),
            )

        plot = Plot(entries,
                    what='hist',
                    title=title,
                    xtitle="Generator bin number",
                    ytitle="N",
                    subplot_type='ratio' if ((do_reco_mc and do_reco_data) or (do_reco_mc_bg_sub and do_reco_data_bg_sub)) else None,
                    subplot_title='#splitline{Data / MC}{(detector)}',
                    subplot_limits=(0.8, 1.2),
                    has_data=(do_reco_data or do_reco_data_bg_sub))
        plot.default_canvas_size = (800, 600)
        plot.plot("NOSTACK HISTE")
        plot.set_logy(do_more_labels=False)
        ymax = max([o.GetMaximum() for o in plot.contributions_objs])
        plot.container.SetMaximum(ymax * 200)
        ymin = max([o.GetMinimum(1E-10) for o in plot.contributions_objs])
        plot.container.SetMinimum(ymin*0.01)
        l, t = self.draw_pt_binning_lines(plot, which='gen', axis='x',
                                          do_underflow=True,
                                          do_labels_inside=True,
                                          do_labels_outside=False,
                                          labels_inside_align='lower'
                                          )
        # # plot.container.SetMinimum(0.001)
        plot.legend.SetY1NDC(0.77)
        plot.legend.SetX1NDC(0.65)
        plot.legend.SetX2NDC(0.88)
        output_filename = "%s/detector_gen_binning_%s.%s" % (output_dir, append, self.output_fmt)
        plot.save(output_filename)

    def draw_unfolded_folded(self, output_dir='.', append="", title=""):
        """Draw big 1D plot of folded unfolded + original input (bg-subtracted if possible)"""
        entries = []
        reco_data_colour = ROOT.kRed
        reco_folded_colour = ROOT.kAzure+1
        reco_mc_colour = ROOT.kGreen+2

        hist_reco = self.unfolder.input_hist_bg_subtracted
        is_bg_subtracted = True
        if not hist_reco:
            hist_reco = self.input_hist
            is_bg_subtracted = False
        if hist_reco:
            entries.append(
                Contribution(hist_reco,
                             label="Unfolding input (background-subtracted)" if is_bg_subtracted else "Unfolding input",
                             line_color=reco_data_colour, line_width=1,
                             marker_color=reco_data_colour, marker_size=0.6, marker_style=20,
                             normalise_hist=False)
                )

        hist_folded = self.unfolder.get_folded_unfolded()
        if hist_folded:
            entries.append(
                Contribution(hist_folded, label="Folded unfolded result (#tau = %.3g)" % (self.unfolder.tau),
                             line_color=reco_folded_colour, line_width=1,
                             marker_color=reco_folded_colour, marker_size=0,
                             normalise_hist=False,
                             subplot=hist_reco)
                )

        plot = Plot(entries,
                    what='hist',
                    title=title,
                    xtitle="Detector binning",
                    ytitle="N",
                    subplot_type='ratio',
                    subplot_title='#splitline{Folded /}{Unfolding input}',
                    subplot_limits=(0.8, 1.2))
        plot.default_canvas_size = (800, 600)
        plot.plot("NOSTACK HIST")
        plot.main_pad.SetLogy(1)
        ymax = max([o.GetMaximum() for o in plot.contributions_objs])
        plot.container.SetMaximum(ymax * 200)
        # plot.container.SetMinimum(1E-8)
        plot.legend.SetY1NDC(0.77)
        plot.legend.SetX1NDC(0.6)
        plot.legend.SetX2NDC(0.88)
        plot.legend.SetY2NDC(0.87)

        l, t = self.draw_pt_binning_lines(plot, which='reco', axis='x',
                                          do_underflow=True,
                                          do_labels_inside=True,
                                          do_labels_outside=False)
        output_filename = "%s/folded_unfolded_%s.%s" % (output_dir, append, self.output_fmt)
        plot.save(output_filename)

    def draw_truth_folded(self, output_dir='.', append="", title=""):
        """Draw big 1D plot of folded truth + original input (bg-subtracted if possible)"""
        # TODO assimilate with draw_unfolded_folded()
        entries = []

        reco_data_colour = ROOT.kRed
        reco_folded_colour = ROOT.kAzure+1
        reco_mc_colour = ROOT.kGreen+2

        hist_reco = self.unfolder.input_hist_bg_subtracted
        is_bg_subtracted = True
        if not hist_reco:
            hist_reco = self.input_hist
            is_bg_subtracted = False
        if hist_reco:
            entries.append(
                Contribution(hist_reco,
                             label="Unfolding input (background-subtracted)" if is_bg_subtracted else "Unfolding input",
                             line_color=reco_mc_colour, line_width=1,
                             marker_color=reco_mc_colour,  marker_size=0.6, marker_style=20,
                             normalise_hist=False)
                )

        hist_folded = self.unfolder.get_folded_mc_truth()
        if hist_folded:
            entries.append(
                Contribution(hist_folded, label="Folded generator",
                             line_color=reco_folded_colour, line_width=1,
                             marker_color=reco_folded_colour, marker_size=0,
                             normalise_hist=False,
                             subplot=hist_reco)
                )

        plot = Plot(entries,
                    what='hist',
                    title=title,
                    xtitle="Detector binning",
                    ytitle="N",
                    subplot_type='ratio',
                    subplot_title='#splitline{Folded gen/}{Unfolding input}',
                    subplot_limits=(0.8, 1.2))
        plot.default_canvas_size = (800, 600)
        plot.plot("NOSTACK HIST")
        plot.main_pad.SetLogy(1)
        ymax = max([o.GetMaximum() for o in plot.contributions_objs])
        plot.container.SetMaximum(ymax * 200)
        # plot.container.SetMinimum(1E-8)
        plot.legend.SetY1NDC(0.77)
        plot.legend.SetX1NDC(0.6)
        plot.legend.SetX2NDC(0.88)
        plot.legend.SetY2NDC(0.87)


        l, t = self.draw_pt_binning_lines(plot, which='reco', axis='x',
                                          do_underflow=True,
                                          do_labels_inside=True,
                                          do_labels_outside=False)
        output_filename = "%s/folded_gen_%s.%s" % (output_dir, append, self.output_fmt)
        plot.save(output_filename)
