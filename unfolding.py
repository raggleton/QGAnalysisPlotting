#!/usr/bin/env python


"""TUnfold it all

Thanks to Ashley, Dennis
"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse
from array import array
import numpy as np
import math
import distutils
from distutils import util

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()
ROOT.gErrorIgnoreLevel = ROOT.kWarning


import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat(".3f")


# Control plot output format
OUTPUT_FMT = "pdf"


class TauScanner(object):
    """Class to handle doing ScanTau on a TUnfoldBinning object,
    since it produces many associated objects & values.

    Can also plot results from scanning.
    """
    def __init__(self):
        self.scan_results = ROOT.MakeNullPointer(ROOT.TSpline)
        self.l_curve = ROOT.MakeNullPointer(ROOT.TGraph)
        self.log_tau_x = ROOT.MakeNullPointer(ROOT.TSpline)
        self.log_tau_y = ROOT.MakeNullPointer(ROOT.TSpline)
        self.ind_best_point = 0
        self.tau = 0
        self.scan_mode = ""
        self.graph_all_scan_points = None
        self.graph_best_scan_point = None

    def scan_tau(self, tunfolder, n_scan, tau_min, tau_max, scan_mode, distribution, axis_steering):
        self.ind_best_point = tunfolder.ScanTau(n_scan,
                                                tau_min,
                                                tau_max,
                                                self.scan_results,
                                                scan_mode,
                                                distribution,
                                                axis_steering,
                                                self.l_curve,
                                                self.log_tau_x,
                                                self.log_tau_y)
        self.tau = tunfolder.GetTau()
        self._process_results(scan_mode)
        print("scan_tau value is {}".format(self.tau))
        print("chi**2 A {:3.1f} + chi**2 L {:3.1f} / NDOF {:3.1f} ".format(tunfolder.GetChi2A(),
                                                                           tunfolder.GetChi2L(),
                                                                           tunfolder.GetNdf()))
        return self.tau

    def _process_results(self, scan_mode):
        """Create graphs etc from ScanTau output

        User shouldn't call this, only internal
        """
        # these need to be arrays to go into TGraph later
        # t here is log_10(tau)
        x, y, t, rho = array('d'), array('d'), array('d'), array('d')

        # Get best scan point & make graph of it
        t0 = ROOT.Double(0.0)
        rho0 = ROOT.Double(0.0)
        self.scan_results.GetKnot(self.ind_best_point, t0, rho0)
        t.append(t0)
        rho.append(rho0)
        self.graph_best_scan_point = ROOT.TGraph(1, t, rho)

        # x0 = ROOT.Double(0.0)
        # y0 = ROOT.Double(0.0)
        # self.l_curve.GetPoint(self.ind_best_point, x0, y0)
        # x.append(x0)
        # y.append(y0)
        print("t[0] =", t[0])
        print("rho[0] =", rho[0])
        # print("x[0] =", x[0])
        # print("y[0] =", y[0])
        print("10^log_10(tau) = tau =", math.pow(10., float(t0)))

        # Make graph of all the points scanned
        t_all, rho_all = array('d'), array('d')
        n_scan = self.scan_results.GetNp()
        for i in range(n_scan):
            tt = ROOT.Double(0.0)
            rr = ROOT.Double(0.0)
            self.scan_results.GetKnot(i, tt, rr)
            t_all.append(tt)
            rho_all.append(rr)

        self.graph_all_scan_points = ROOT.TGraph(int(n_scan), t_all, rho_all)

        tau_mode_dict = {
            ROOT.TUnfoldDensity.kEScanTauRhoAvg: "average (stat+bgr) global correlation (#rho)",
            ROOT.TUnfoldDensity.kEScanTauRhoAvgSys: "average (stat+bgr+sys) global correlation (#rho)",
            ROOT.TUnfoldDensity.kEScanTauRhoMax: "maximum (stat+bgr) global correlation (#rho)",
            ROOT.TUnfoldDensity.kEScanTauRhoMaxSys: "maximum (stat+bgr+sys) global correlation (#rho)",
            ROOT.TUnfoldDensity.kEScanTauRhoSquareAvg: "average (stat+bgr) global correlation (#rho) squared",
            ROOT.TUnfoldDensity.kEScanTauRhoSquareAvgSys: "average (stat+bgr+sys) global correlation (#rho) squared",
        }
        self.graph_all_scan_points.SetTitle("Optimization of Regularization Parameter, #tau : Scan of {}".format(tau_mode_dict[scan_mode]))

    def plot_scan_tau(self, output_filename):
        """Plot graph of scan results, and optimum tau"""
        canv_tau_scan = ROOT.TCanvas("canv_tau_scan"+str(self.tau), "canv_tau_scan"+str(self.tau))

        self.graph_all_scan_points.SetLineColor(ROOT.kBlue+3)
        self.graph_all_scan_points.Draw()

        self.graph_best_scan_point.SetMarkerColor(ROOT.kRed)
        self.graph_best_scan_point.Draw("* same")

        self.graph_all_scan_points.GetXaxis().SetTitle("log_{10}(#tau)")
        self.graph_all_scan_points.GetYaxis().SetTitle(" #rho")

        leg = ROOT.TLegend(0.1, 0.6, 0.25, 0.89)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.026)
        leg.AddEntry(self.graph_all_scan_points, 'Scan over #tau', 'l')
        leg.AddEntry(self.graph_best_scan_point, 'Chosen point: #tau = {}'.format(self.tau), 'P')
        leg.Draw()

        canv_tau_scan.Print(output_filename)


class MyUnfolder(object):
    # Control plot output format
    OUTPUT_FMT = "pdf"

    def __init__(self,
                 response_map,  # 2D GEN-RECO heatmap
                 variable_bin_edges_reco, # 'variable' refers to e.g. ptD, LHA
                 variable_bin_edges_gen, # reco for detector binnig, gen for generator (final) binning
                 variable_name,
                 pt_bin_edges_reco,
                 pt_bin_edges_gen,
                 pt_bin_edges_underflow_reco,
                 pt_bin_edges_underflow_gen,
                 orientation=ROOT.TUnfold.kHistMapOutputHoriz,
                 constraintMode=ROOT.TUnfold.kEConstraintArea,
                 regMode=ROOT.TUnfold.kRegModeCurvature,
                 densityFlags=ROOT.TUnfoldDensity.kDensityModeBinWidth,
                 axisSteering='*[b]'):

        self.response_map = response_map
        self.variable_name = variable_name
        self.variable_name_safe = variable_name.replace(" ", "_")

        self.variable_bin_edges_reco = variable_bin_edges_reco
        self.nbins_variable_reco = len(variable_bin_edges_reco)-1 if variable_bin_edges_reco is not None else 0
        self.variable_bin_edges_gen = variable_bin_edges_gen
        self.nbins_variable_gen = len(variable_bin_edges_gen)-1 if variable_bin_edges_gen is not None else 0

        self.pt_bin_edges_reco = pt_bin_edges_reco
        self.nbins_pt_reco = len(pt_bin_edges_reco)-1
        self.pt_bin_edges_gen = pt_bin_edges_gen
        self.nbins_pt_gen = len(pt_bin_edges_gen)-1

        self.pt_bin_edges_underflow_reco = pt_bin_edges_underflow_reco
        self.nbins_pt_underflow_reco = len(pt_bin_edges_underflow_reco)-1
        self.pt_bin_edges_underflow_gen = pt_bin_edges_underflow_gen
        self.nbins_pt_underflow_gen = len(pt_bin_edges_underflow_gen)-1

        # Binning setup here MUST match how it was setup in making the input files, otherwise
        # you will have untold pain and suffering!
        # TODO read in from XML
        var_uf, var_of = False, False
        pt_uf, pt_of = False, False  # handle pt uder/over flow ourselves
        self.detector_binning = ROOT.TUnfoldBinning("detector")
        self.detector_distribution_underflow = self.detector_binning.AddBinning("detectordistribution_underflow")
        if self.variable_bin_edges_reco is not None:
            self.detector_distribution_underflow.AddAxis(self.variable_name, self.nbins_variable_reco, self.variable_bin_edges_reco, var_uf, var_of)
        self.detector_distribution_underflow.AddAxis("pt", self.nbins_pt_underflow_reco, self.pt_bin_edges_underflow_reco, pt_uf, pt_of)

        self.detector_distribution = self.detector_binning.AddBinning("detectordistribution")
        if self.variable_bin_edges_reco is not None:
            self.detector_distribution.AddAxis(self.variable_name, self.nbins_variable_reco, self.variable_bin_edges_reco, var_uf, var_of)
        self.detector_distribution.AddAxis("pt", self.nbins_pt_reco, self.pt_bin_edges_reco, pt_uf, pt_of)


        self.generator_binning = ROOT.TUnfoldBinning("generator")

        self.generator_distribution_underflow = self.generator_binning.AddBinning("generatordistribution_underflow")
        if self.variable_bin_edges_gen is not None:
            self.generator_distribution_underflow.AddAxis(self.variable_name, self.nbins_variable_gen, self.variable_bin_edges_gen, var_uf, var_of)
        self.generator_distribution_underflow.AddAxis("pt", self.nbins_pt_underflow_gen, self.pt_bin_edges_underflow_gen, pt_uf, pt_of)
        self.generator_distribution = self.generator_binning.AddBinning("generatordistribution")
        if self.variable_bin_edges_gen is not None:
            self.generator_distribution.AddAxis(self.variable_name, self.nbins_variable_gen, self.variable_bin_edges_gen, var_uf, var_of)
        self.generator_distribution.AddAxis("pt", self.nbins_pt_gen, self.pt_bin_edges_gen, pt_uf, pt_of)

        self.orientation = orientation
        self.constraintMode = constraintMode
        self.regMode = regMode
        self.densityFlags = densityFlags
        self.axisSteering = axisSteering

        self.tunfolder = ROOT.TUnfoldDensity(response_map,
                                            orientation,
                                            regMode,
                                            constraintMode,
                                            densityFlags,
                                            self.generator_binning,
                                            self.detector_binning)
                                            # hmm comment these out when scanL or scanTau?
                                            # "generator",
                                            # self.axisSteering)

        self.use_axis_binning = False  # for things like get_probability_matrix()

    def save_binning(self, print_xml=True, txt_filename=None):
        """Save binning scheme to txt and/or print XML to screen"""
        if txt_filename:
            with open(txt_filename, 'w') as of:
                of.write("GEN BINNING\n")
                of.write("--------------------\n")
                for name, region in [
                    ("generator_distribution", self.generator_distribution),
                    ("generator_distribution_underflow", self.generator_distribution_underflow)]:
                    of.write("%s: bins %d - %d\n" % (name, region.GetStartBin(), region.GetEndBin()))
                of.write("\nDETECTOR BINNING\n")
                of.write("--------------------\n")
                for name, region in [
                    ("detector_distribution", self.detector_distribution),
                    ("detector_distribution_underflow", self.detector_distribution_underflow)]:
                    of.write("%s: bins %d - %d\n" % (name, region.GetStartBin(), region.GetEndBin()))

        if print_xml:
            # don't know how to create a ofstream in python :( best we can do is ROOT.cout
            ROOT.TUnfoldBinningXML.ExportXML(self.detector_binning, ROOT.cout, True, False)
            ROOT.TUnfoldBinningXML.ExportXML(self.generator_binning, ROOT.cout, False, True)

    def setInput(self, hist):
        self.tunfolder.SetInput(hist)

    def doScanL(self, output_dir, n_scan=300, tau_min=1E-14, tau_max=0.1):
        """Figure out best tau by doing scan over L curve
        Taken from Ashley's code
        """
        scannedlcurve = ROOT.MakeNullPointer(ROOT.TGraph)
        logTauX = ROOT.MakeNullPointer(ROOT.TSpline3)
        logTauY = ROOT.MakeNullPointer(ROOT.TSpline3)
        logTauCurvature = ROOT.MakeNullPointer(ROOT.TSpline3)

        best = self.tunfolder.ScanLcurve(n_scan,
                                        tau_min,
                                        tau_max,
                                        scannedlcurve,
                                        logTauX,
                                        logTauY,
                                        logTauCurvature)

        xx, yy, tt = array('d'), array('d'), array('d')
        # save graphs with one point to visualize best choice of tau
        X0, Y0, T0 = ROOT.Double(0.0), ROOT.Double(0.0), ROOT.Double(0.0)

        for i in range(n_scan):
            logTauX.GetKnot(i, T0, X0)
            logTauY.GetKnot(i, T0, Y0)
            tt.append(T0)
            xx.append(X0)
            yy.append(Y0)
        print("printing X0 Y0 and T0 : ")
        print(T0)
        print(X0)
        print(Y0)

        bestLcurve = ROOT.TGraph(len(xx), xx, yy)  # int(nScan),xx,yy);
        bestLogTauLogChi2 = ROOT.TGraph(best, xx, yy)
        tau = self.tunfolder.GetTau()
        print("doScanL value is {}".format(tau))
        print("Returned the TUnfoldDensity object.")

        canv_lScan = ROOT.TCanvas("canv_lScan"+str(tau), "canv_lScan"+str(tau))

        bestLcurve.SetLineColor(ROOT.kBlue+2)

        bestLogTauLogChi2.SetMarkerColor(ROOT.kRed)

        bestLcurve.Draw("")
        bestLogTauLogChi2.Draw("* same")
        bestLcurve.GetXaxis().SetTitle("log_{10}(L_{1})")
        bestLcurve.GetYaxis().SetTitle("log_{10}(#frac{L_{2}}{#tau^{2}})")  # frac{Theory}{Unfolded MC}

        bestLcurve.SetTitle("Optimization of Regularization Parameter, #tau : Scan of L Curve")

        leg1 = ROOT.TLegend(0.5, 0.6, 0.7, 0.94)
        leg1.SetFillColor(0)
        leg1.SetFillStyle(0)
        leg1.SetBorderSize(0)
        leg1.SetTextSize(0.026)
        leg1.AddEntry(bestLcurve, 'L Curve Scan: #tau = {}'.format(tau), 'l')
        leg1.Draw()

        canv_lScan.Modified()
        canv_lScan.Update()

        canv_lScan.Print("%s/scanL_%s.%s" % (output_dir, self.variable_name, self.OUTPUT_FMT))
        return tau

    def do_unfolding(self, tau):
        print(">>> Unfolding with tau =", tau)
        self.tunfolder.DoUnfold(tau)

    def get_output(self):
        """Get 1D unfolded histogram covering all bins"""
        # print("( " + str(self.tunfolder.GetChi2A()) + " + " + str(self.tunfolder.GetChi2L()) + ") / " + str(self.tunfolder.GetNdf()))
        # use "generator" for signal + underflow region, "generatordistribution" for only signal region
        return self.tunfolder.GetOutput("unfolded_" + cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)


    def get_bias(self):
        return self.tunfolder.GetBias("bias_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)

    def get_ematrix_input(self):
        return self.tunfolder.GetEmatrixInput("ematrix_input_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)

    def get_ematrix_sys_uncorr(self):
        return self.tunfolder.GetEmatrixSysUncorr("ematrix_sys_uncorr_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)

    def get_ematrix_total(self):
        return self.tunfolder.GetEmatrixTotal("ematrix_total_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)

    def get_rhoij_total(self):
        return self.tunfolder.GetRhoIJtotal("rhoij_total_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)

    def get_probability_matrix(self):
        return self.tunfolder.GetProbabilityMatrix("prob_matrix_"+cu.get_unique_str(), "", self.use_axis_binning)

    def get_var_hist_pt_binned(self, hist1d, ibin_pt, binning_scheme='generator'):
        """Get hist of variable for given pt bin from massive 1D hist that TUnfold makes"""
        binning = self.generator_binning.FindNode("generatordistribution") if binning_scheme == "generator" else self.detector_binning.FindNode("detectordistribution")
        var_bins = np.array(binning.GetDistributionBinning(0))
        pt_bins = np.array(binning.GetDistributionBinning(1))

        # print("var_bins:", var_bins)
        # print("pt_bins:", pt_bins)
        bin_num = binning.GetGlobalBinNumber(0.001, 49)
        # print("Global bin num for (pt, lambda) = (49, 0.001) => %d" % (bin_num))

        # need the -1 on ibin_pt, as it references an array index, whereas ROOT bins start at 1
        h = ROOT.TH1D("h_%d_%s" % (ibin_pt, cu.get_unique_str()), "", len(var_bins)-1, var_bins)
        for var_ind, var_value in enumerate(var_bins[:-1], 1):
            this_val = var_value * 1.001  # ensure its inside
            bin_num = binning.GetGlobalBinNumber(this_val, pt_bins[ibin_pt]*1.001)
            # print("Global bin num for (pt, lambda) = (%.3f, %.3f) => %d" % (pt_bins[ibin_pt]*1.001, this_val, bin_num))
            h.SetBinContent(var_ind, hist1d.GetBinContent(bin_num))
            h.SetBinError(var_ind, hist1d.GetBinError(bin_num))
            # print("Bin:", bin_num, this_val, pt_bins[ibin_pt], "=", hist1d.GetBinContent(bin_num), "+-", hist1d.GetBinError(bin_num))
        return h

    def get_folded_output(self, hist_name):
        return self.tunfolder.GetFoldedOutput(hist_name)


def generate_2d_canvas(size=(800, 600)):
    canv = ROOT.TCanvas(cu.get_unique_str(), "", *size)
    canv.SetTicks(1, 1)
    canv.SetRightMargin(1.5)
    canv.SetLeftMargin(0.9)
    return canv


def draw_response_matrix(rsp_map, region_name, variable_name, output_filename):
    canv = generate_2d_canvas()
    canv.SetLogz()

    rsp_map.SetTitle("Response matrix, %s region, %s;Generator Bin;Detector Bin" % (region_name, variable_name))
    rsp_map.GetYaxis().SetTitleOffset(1.5)
    rsp_map.GetXaxis().SetTitleOffset(1.5)
    rsp_map.Draw("COLZ")
    canv.SaveAs(output_filename)


def draw_probability_matrix(prob_map, region_name, variable_name, output_filename):
    canv = generate_2d_canvas()
    canv.SetLogz()

    title = "Probability map, %s region, %s;Generator Bin;Detector Bin" % (region_name, variable_name)
    prob_map.SetTitle(title)
    prob_map.GetYaxis().SetTitleOffset(1.5)
    prob_map.GetXaxis().SetTitleOffset(1.5)
    prob_map.Draw("COLZ")
    canv.SaveAs(output_filename)


def draw_error_matrix_input(err_map, region_name, variable_name, output_filename):
    canv = generate_2d_canvas()
    canv.SetLogz()

    title = "Error matrix (input), %s region, %s;Generator bin;Generator bin" % (region_name, variable_name)
    err_map.SetTitle(title)
    err_map.GetYaxis().SetTitleOffset(1.5)
    err_map.GetXaxis().SetTitleOffset(1.5)
    err_map.Draw("COLZ")
    canv.SaveAs(output_filename)


def draw_error_matrix_sys_uncorr(err_map, region_name, variable_name, output_filename):
    canv = generate_2d_canvas()
    canv.SetLogz()

    title = "Error matrix (sys uncorr), %s region, %s;Generator bin;Generator bin" % (region_name, variable_name)
    err_map.SetTitle(title)
    err_map.GetYaxis().SetTitleOffset(1.5)
    err_map.GetXaxis().SetTitleOffset(1.5)
    err_map.Draw("COLZ")
    canv.SaveAs(output_filename)


def draw_error_matrix_total(err_map, region_name, variable_name, output_filename):
    canv = generate_2d_canvas()
    canv.SetLogz()

    title = "Error matrix (total), %s region, %s;Generator bin;Generator bin" % (region_name, variable_name)
    err_map.SetTitle(title)
    err_map.GetYaxis().SetTitleOffset(1.5)
    err_map.GetXaxis().SetTitleOffset(1.5)
    err_map.Draw("COLZ")
    canv.SaveAs(output_filename)


def draw_correlation_matrix(corr_map, region_name, variable_name, output_filename):
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

    canv = generate_2d_canvas()

    title = "Correlation map, %s region, %s;Generator Bin;Generator Bin" % (region_name, variable_name)
    corr_map.SetTitle(title)
    corr_map.GetYaxis().SetTitleOffset(1.5)
    corr_map.GetXaxis().SetTitleOffset(1.5)
    # ROOT.gStyle.SetPalette(ROOT.kLightTemperature)
    corr_map.SetMinimum(-1)
    corr_map.SetMaximum(1)
    corr_map.SetMarkerSize(0.05)
    corr_map.Draw("COLZ0 TEXT45")
    canv.SaveAs(output_filename)
    ROOT.gStyle.SetPalette(ROOT.kBird)


def plot_simple_unfolded(unfolded, tau, reco, gen, fake, output_filename, title=""):
    """Simple plot of unfolded, reco, gen, by bin number (ie non physical axes)"""
    entries = []

    if reco:
        entries.append(
            Contribution(reco, label="Reco",
                         line_color=ROOT.kGreen+2, line_width=1,
                         marker_color=ROOT.kGreen+2, marker_size=0,
                         normalise_hist=False),
        )

    if gen:
        entries.append(
            Contribution(gen, label="Gen",
                         line_color=ROOT.kBlue, line_width=1,
                         marker_color=ROOT.kBlue, marker_size=0,
                         normalise_hist=False),
        )

    if fake:
        entries.append(
            Contribution(fake, label="Fakes",
                         line_color=ROOT.kOrange+4, line_width=1,
                         marker_color=ROOT.kOrange+4, marker_size=0,
                         normalise_hist=False),
        )

    if unfolded:
        entries.append(
            Contribution(unfolded, label="Unfolded (#tau = %.3g)" % (tau),
                         line_color=ROOT.kRed, line_width=0,
                         marker_color=ROOT.kRed, marker_size=0.6, marker_style=20,
                         normalise_hist=False, subplot=gen),
        )

    plot = Plot(entries,
                what='hist',
                title=title,
                xtitle="Bin number",
                ytitle="N",
                subplot_type='ratio',
                subplot_title='Data / MC',
                subplot_limits=(0.8, 1.2))
    plot.default_canvas_size = (800, 600)
    plot.plot("NOSTACK HISTE")
    plot.main_pad.SetLogy(1)
    ymax = max(h.GetMaximum() for h in [reco, gen, fake, unfolded] if h)
    plot.container.SetMaximum(ymax * 100)
    plot.container.SetMinimum(0.1)
    plot.legend.SetY1NDC(0.77)
    plot.legend.SetX1NDC(0.65)
    plot.legend.SetX2NDC(0.88)
    plot.save(output_filename)


def plot_simple_detector(reco_data, reco_data_fake, reco_mc, reco_mc_fake, output_filename, title):
    """Plot detector-level quantities for data & MC, by bin number (ie non physical axes)"""
    entries = []

    if reco_mc:
        entries.append(
            Contribution(reco_mc, label="MC [detector-level]",
                         line_color=ROOT.kGreen+2, line_width=1,
                         marker_color=ROOT.kGreen+2, marker_size=0,
                         normalise_hist=False),
        )

    if reco_mc_fake:
        entries.append(
            Contribution(reco_mc_fake, label="MC fakes [detector-level]",
                         line_color=ROOT.kOrange+4, line_width=1,
                         marker_color=ROOT.kOrange+4, marker_size=0,
                         normalise_hist=False),
        )

    if reco_data_fake:
        entries.append(
            Contribution(reco_data_fake, label="Data fakes template [detector-level]",
                         line_color=ROOT.kMagenta+2, line_width=1,
                         marker_color=ROOT.kMagenta+2, marker_size=0,
                         normalise_hist=False),
        )

    if reco_data:
        entries.append(
            Contribution(reco_data, label="Data [detector-level]",
                         line_color=ROOT.kRed, line_width=0,
                         marker_color=ROOT.kRed, marker_size=0.6, marker_style=20,
                         normalise_hist=False, subplot=reco_mc),
        )

    plot = Plot(entries,
                what='hist',
                title=title,
                xtitle="Bin number",
                ytitle="N",
                subplot_type='ratio',
                subplot_title='Data / MC',
                subplot_limits=(0.8, 1.2))
    plot.default_canvas_size = (800, 600)
    plot.plot("NOSTACK HISTE")
    plot.main_pad.SetLogy(1)
    ymax = max(h.GetMaximum() for h in [reco_mc, reco_mc_fake, reco_data] if h)
    plot.container.SetMaximum(ymax * 100)
    plot.container.SetMinimum(0.1)
    plot.legend.SetY1NDC(0.77)
    plot.legend.SetX1NDC(0.65)
    plot.legend.SetX2NDC(0.88)
    plot.save(output_filename)


def plot_fake_fraction(hist, output_filename, title):
    entries = [Contribution(hist,
                            label='From MC',
                            line_color=ROOT.kBlue, line_width=1,
                            marker_color=ROOT.kBlue, marker_size=0)]
    plot = Plot(entries,
                what='hist',
                title=title,
                xtitle='Bin number',
                ytitle='Fake fraction')
    plot.default_canvas_size = (800, 600)
    plot.plot('NOSTACK HIST')
    plot.set_logy()
    plot.container.SetMaximum(hist.GetMaximum() * 100)
    plot.container.SetMinimum(hist.GetMinimum(1E-12) * 0.1)
    plot.save(output_filename)


def create_hist_with_errors(hist, err_matrix):
    hnew = hist.Clone(cu.get_unique_str())
    nbins = hist.GetNbinsX()
    for i in range(1, nbins+1):
        err = math.sqrt(err_matrix.GetBinContent(i, i))
        hnew.SetBinError(i, err)
    return hnew


def make_hist_from_diagonal_errors(h2d, do_sqrt=True):
    nbins = h2d.GetNbinsX()
    hnew = ROOT.TH1D("h_diag" + cu.get_unique_str(), "", nbins, 0, nbins)
    for i in range(1, nbins+1):
        err = h2d.GetBinContent(i, i)
        if do_sqrt and err > 0:
            err = math.sqrt(err)
        hnew.SetBinError(i, err)
    return hnew


def update_hist_bin_content(h_orig, h_to_be_updated):
    if h_orig.GetNbinsX() != h_to_be_updated.GetNbinsX():
        raise RuntimeError("Need same # x bins, %d vs %s" % (h_orig.GetNbinsX(), h_to_be_updated.GetNbinsX()))
    for i in range(0, h_orig.GetNbinsX()+2):
        h_to_be_updated.SetBinContent(i, h_orig.GetBinContent(i))
        # h_to_be_updated.SetBinError(i, h_orig.GetBinError(i))


def update_hist_bin_error(h_orig, h_to_be_updated):
    if h_orig.GetNbinsX() != h_to_be_updated.GetNbinsX():
        raise RuntimeError("Need same # x bins, %d vs %s" % (h_orig.GetNbinsX(), h_to_be_updated.GetNbinsX()))
    for i in range(0, h_orig.GetNbinsX()+2):
        h_to_be_updated.SetBinError(i, h_orig.GetBinError(i))


def draw_projection_comparison(h_orig, h_projection, title, xtitle, output_filename, print_bin_comparison=True):
    """Draw 2 hists, h_orig the original, and h_projection the projection of a 2D hist"""

    # Check integrals
    int_orig = h_orig.Integral()
    int_proj = h_projection.Integral()
    if abs(int_orig - int_proj)/int_orig > 0.01:
        print("draw_projection_comparison: different integrals: %f vs %f" % (int_orig, int_proj))

    # Check bin-by-bin
    if print_bin_comparison:
        for i in range(1, h_orig.GetNbinsX()+1):
            value_orig = h_orig.GetBinContent(i)
            value_proj = h_projection.GetBinContent(i)
            if value_orig == 0 and value_proj == 0:
                continue
            rel_diff = abs((value_orig - value_proj)/max(abs(value_orig), abs(value_proj)))
            if rel_diff > 1E-3:
                print("draw_projection_comparison: bin %s has different contents: %f vs %f (rel diff %f)" % (i, value_orig, value_proj, rel_diff))

    entries = [
        Contribution(h_orig, label="1D hist",
                     line_color=ROOT.kBlue, line_width=1,
                     marker_color=ROOT.kBlue, marker_size=0,
                     normalise_hist=False),
        Contribution(h_projection, label="Response map projection",
                     line_color=ROOT.kRed, line_width=1,
                     marker_color=ROOT.kRed, marker_size=0,
                     normalise_hist=False,
                     subplot=h_orig),
    ]
    plot = Plot(entries,
                what='hist',
                title=title,
                xtitle=xtitle,
                ytitle="N",
                subplot_type='ratio',
                subplot_title='Projection / 1D',
                subplot_limits=(0.999, 1.001))
    plot.default_canvas_size = (800, 600)
    plot.plot("NOSTACK HIST")
    plot.main_pad.SetLogy(1)
    ymax = max(h.GetMaximum() for h in [h_orig, h_projection])
    plot.container.SetMaximum(ymax * 10)
    # plot.container.SetMinimum(1E-8)
    plot.legend.SetY1NDC(0.77)
    plot.legend.SetX2NDC(0.85)
    plot.save(output_filename)


def draw_reco_folded(hist_folded, tau, hist_reco_data, hist_reco_mc, title, xtitle, output_filename):
    entries = [
        Contribution(hist_reco_mc, label="MC (detector, bg-subtracted)",
                     line_color=ROOT.kBlack, line_width=1,
                     marker_color=ROOT.kBlack, marker_size=0,
                     normalise_hist=False),
        Contribution(hist_reco_data, label="Data (detector, bg-subtracted)",
                     line_color=ROOT.kBlue, line_width=1,
                     marker_color=ROOT.kBlue, marker_size=0,
                     normalise_hist=False),
        Contribution(hist_folded, label="Unfolded data (#tau = %.3g, folded)" % (tau),
                     line_color=ROOT.kRed, line_width=1,
                     marker_color=ROOT.kRed, marker_size=0,
                     normalise_hist=False,
                     subplot=hist_reco_data),
    ]
    plot = Plot(entries,
                what='hist',
                title=title,
                xtitle=xtitle,
                ytitle="N",
                subplot_type='ratio',
                subplot_title='Folded / Data',
                subplot_limits=(0.8, 1.2))
    plot.default_canvas_size = (800, 600)
    plot.plot("NOSTACK HIST")
    plot.main_pad.SetLogy(1)
    ymax = max(h.GetMaximum() for h in [hist_reco_data, hist_reco_mc, hist_folded])
    plot.container.SetMaximum(ymax * 50)
    # plot.container.SetMinimum(1E-8)
    plot.legend.SetY1NDC(0.77)
    plot.legend.SetX2NDC(0.85)
    plot.save(output_filename)


def check_entries(entries, message=""):
    """Check that at least 1 Contribution has something in it"""
    has_entries = [c.obj.GetEntries() > 0 for c in entries]
    if not any(has_entries):
        if message:
            print("Skipping 0 entries (%s)" % (message))
        else:
            print("Skipping 0 entries")
        return False
    return True


def plot_uncertainty_shifts(total_hist, stat_hist, syst_shifts, systs, output_filename, title, angle_str):
    """
    Parameters
    ----------
    total_hist : TH1
        Distribution with total uncertainty
    stat_hist : TH1
        Distribution with stat-only uncertainties
    syst_shifts : list[TH1]
        Distributions with 1-sigma shift from systematics
    systs : list[dict]
        Dicts describing each syst
    output_filename : str
        Output plot filename
    title : str
        Title to put on plot
    angle_str : str
        Angle name for x axis

    Returns
    -------
    None
        If all hists are empty
    """
    entries = []
    hists = []
    for h, syst_dict in zip(syst_shifts, systs):
        h_fraction = h.Clone()
        h_fraction.Divide(total_hist)
        hists.append(h_fraction)
        for i in range(1, h_fraction.GetNbinsX()+1):
            h_fraction.SetBinContent(i, abs(h_fraction.GetBinContent(i)))
        c = Contribution(h_fraction,
                         label=syst_dict['label'],
                         line_color=syst_dict['colour'],
                         line_style=syst_dict.get('linestyle', 1),
                         line_width=2,
                         marker_size=0,
                         marker_color=syst_dict['colour'],
                         )
        entries.append(c)

    # Add systematic
    h_syst = stat_hist.Clone()
    h_total = total_hist.Clone()
    for i in range(1, h_syst.GetNbinsX()+1):
        if total_hist.GetBinContent(i) > 0:
            h_syst.SetBinContent(i, stat_hist.GetBinError(i) / total_hist.GetBinContent(i))
            h_total.SetBinContent(i, total_hist.GetBinError(i) / total_hist.GetBinContent(i))
        else:
            h_syst.SetBinContent(i, 0)
            h_total.SetBinContent(i, 0)
        h_syst.SetBinError(i, 0)
        h_total.SetBinError(i, 0)
    c_stat = Contribution(h_syst,
                         label="Stat.",
                         line_color=ROOT.kRed,
                         line_style=3,
                         line_width=3,
                         marker_size=0,
                         marker_color=ROOT.kRed,
                         )
    entries.append(c_stat)
    c_tot = Contribution(h_total,
                         label="Total",
                         line_color=ROOT.kBlack,
                         line_style=1,
                         line_width=3,
                         marker_size=0,
                         marker_color=ROOT.kBlack,
                         )
    entries.append(c_tot)

    if not check_entries(entries, "systematic shifts"):
        return
    plot = Plot(entries,
                what="hist",
                title=title,
                xtitle="Particle-level "+angle_str,
                ytitle="| Fractional shift |")
    plot.legend.SetX1(0.55)
    plot.legend.SetY1(0.68)
    plot.legend.SetX2(0.98)
    plot.legend.SetY2(0.88)
    plot.legend.SetNColumns(2)
    plot.y_padding_max_linear = 1.4
    plot.plot("NOSTACK HIST")
    # plot.set_logy()
    plot.save(output_filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",
                        help="Source directory with ROOT files")
    parser.add_argument("--doDijetCentral",
                        action='store_true',
                        help='Do unfolding for dijet (central) jets')
    parser.add_argument("--doDijetForward",
                        action='store_true',
                        help='Do unfolding for dijet (forward) jets')
    parser.add_argument("--doDijetCentralGroomed",
                        action='store_true',
                        help='Do unfolding for groomed dijet (central) jets')
    parser.add_argument("--doDijetForwardGroomed",
                        action='store_true',
                        help='Do unfolding for groomed dijet (forward) jets')
    parser.add_argument("--doZPJ",
                        action='store_true',
                        help='Do unfolding for Z+jet jets')
    parser.add_argument("--doZPJGroomed",
                        action='store_true',
                        help='Do unfolding for groomed Z+jet jets')

    parser.add_argument("--regularize",
                        choices=['None', 'tau', 'L'],
                        default='None',
                        help='Regularization scheme')

    parser.add_argument("--MCinput",
                        type=lambda x:bool(distutils.util.strtobool(x)),
                        default=False,
                        help='Unfold MC instead of data')

    parser.add_argument("--MCsplit",
                        type=lambda x:bool(distutils.util.strtobool(x)),
                        default=False,
                        help='Split MC between response & 1D reco, good for testing procedure')

    parser.add_argument("--doSysts",
                        type=lambda x:bool(distutils.util.strtobool(x)),
                        default=True,
                        help='Do systematics')

    parser.add_argument("--doSummaryPlot",
                        type=lambda x:bool(distutils.util.strtobool(x)),
                        default=False,
                        help='Do summary plot')

    parser.add_argument("--useAltResponse",
                        type=lambda x:bool(distutils.util.strtobool(x)),
                        default=False,
                        help='Use alternate response matrix to unfold')

    parser.add_argument("--outputDir",
                        default='',
                        help='Output directory')

    parser.add_argument("--angles",
                        choices=list(qgc.VAR_UNFOLD_DICT.keys()),
                        nargs='+',
                        help="Lambda angles to unfold")

    args = parser.parse_args()
    print(args)

    if not any([args.doDijetCentral, args.doDijetForward, args.doDijetCentralGroomed, args.doDijetForwardGroomed, args.doZPJ, args.doZPJGroomed]):
        raise RuntimeError("You need to specify at least one signal region e.g. --doDijetCentral")

    if args.useAltResponse and args.doSysts:
        args.doSysts = False
        print("You cannot use both --useAltResponse and --doSysts: disabling doSysts")

    # Setup files and regions to unfold
    # --------------------------------------------------------------------------
    src_dir = "workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet"
    src_dir = args.source
    src_dir_systs = os.path.join(src_dir, "systematics_files")

    regions = []

    if any([args.doDijetCentral, args.doDijetForward, args.doDijetCentralGroomed, args.doDijetForwardGroomed]):
        # FOR DIJET:
        input_mc_qcd_mgpythia_tfile = cu.open_root_file(os.path.join(src_dir, qgc.QCD_FILENAME))
        input_mc_qcd_pythia_tfile = cu.open_root_file(os.path.join(src_dir, qgc.QCD_PYTHIA_ONLY_FILENAME))
        input_mc_qcd_herwig_tfile = cu.open_root_file(os.path.join(src_dir, qgc.QCD_HERWIG_FILENAME))

        input_jetht_tfile = cu.open_root_file(os.path.join(src_dir, qgc.JETHT_ZB_FILENAME))

        dijet_region_dict_template = {
            "name": "Dijet",
            "dirname": "Dijet_QG_Unfold_central_tighter",
            "label": "Dijet",
            "data_tfile": input_jetht_tfile,
            "mc_tfile": input_mc_qcd_mgpythia_tfile,
            "mc_label": "MG+Pythia8",
            "alt_mc_tfile": input_mc_qcd_herwig_tfile,
            "alt_mc_label": "Herwig++",
            "tau_limits": {
                'jet_puppiMultiplicity': (1E-13, 1E-10),
                'jet_pTD': (1E-13, 1E-10),
                'jet_LHA': (1E-13, 1E-10),
                'jet_width': (1E-13, 1E-10),
                'jet_thrust': (1E-13, 1E-10),
                'jet_puppiMultiplicity_charged': (1E-13, 1E-10),
                'jet_pTD_charged': (1E-13, 1E-10),
                'jet_LHA_charged': (1E-13, 1E-10),
                'jet_width_charged': (1E-13, 1E-10),
                'jet_thrust_charged': (1E-13, 1E-10),
            },
            "systematics": [
                {
                    "label": "Neutral hadron up",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'neutralHadronShiftUp', qgc.QCD_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kOrange-3,
                },
                {
                    "label": "Neutral hadron down",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'neutralHadronShiftDown', qgc.QCD_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kOrange-3,
                    "linestyle": 2,
                },
                {
                    "label": "Photon up",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'photonShiftUp', qgc.QCD_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kMagenta-3,
                },
                {
                    "label": "Photon down",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'photonShiftDown', qgc.QCD_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kMagenta-3,
                    "linestyle": 2,
                },
                {
                    "label": "JEC up",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'jecsmear_directionUp', qgc.QCD_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kGreen+2,
                },
                {
                    "label": "JEC down",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'jecsmear_directionDown', qgc.QCD_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kGreen+2,
                    "linestyle": 2,
                },
                {
                    "label": "JER up",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'jersmear_directionUp', qgc.QCD_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kOrange+3,
                },
                {
                    "label": "JER down",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'jersmear_directionDown', qgc.QCD_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kOrange+3,
                    "linestyle": 2,
                },
                {
                    "label": "Pileup up",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'pileup_directionUp', qgc.QCD_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kBlue-4,
                },
                {
                    "label": "Pileup down",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'pileup_directionDown', qgc.QCD_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kBlue-4,
                    "linestyle": 2,
                },
            ]
        }

        if args.doDijetCentral:
            dijet_region_central_dict = dijet_region_dict_template.copy()
            dijet_region_central_dict['dirname'] = 'Dijet_QG_Unfold_central_tighter'
            dijet_region_central_dict['label'] = 'Dijet central'
            dijet_region_central_dict['name'] = 'Dijet_central'
            regions.append(dijet_region_central_dict)

        if args.doDijetForward:
            dijet_region_forward_dict = dijet_region_dict_template.copy()
            dijet_region_forward_dict['dirname'] = 'Dijet_QG_Unfold_forward_tighter'
            dijet_region_forward_dict['label'] = 'Dijet forward'
            dijet_region_forward_dict['name'] = 'Dijet_forward'
            regions.append(dijet_region_forward_dict)

        if args.doDijetCentralGroomed:
            dijet_region_central_groomed_dict = dijet_region_dict_template.copy()
            dijet_region_central_groomed_dict['dirname'] = 'Dijet_QG_Unfold_central_tighter_groomed'
            dijet_region_central_groomed_dict['label'] = 'Dijet central'
            dijet_region_central_groomed_dict['name'] = 'Dijet_central_groomed'
            regions.append(dijet_region_central_groomed_dict)

        if args.doDijetForwardGroomed:
            dijet_region_forward_groomed_dict = dijet_region_dict_template.copy()
            dijet_region_forward_groomed_dict['dirname'] = 'Dijet_QG_Unfold_forward_tighter_groomed'
            dijet_region_forward_groomed_dict['label'] = 'Dijet forward'
            dijet_region_forward_groomed_dict['name'] = 'Dijet_forward_groomed'
            regions.append(dijet_region_forward_groomed_dict)

    if any([args.doZPJ, args.doZPJGroomed]):
        # FOR Z+JETS:
        input_mc_dy_mgpythia_tfile = cu.open_root_file(os.path.join(src_dir, qgc.DY_FILENAME))
        input_mc_dy_mgherwig_tfile = cu.open_root_file(os.path.join(src_dir, qgc.DY_MG_HERWIG_FILENAME))
        input_mc_dy_herwig_tfile = cu.open_root_file(os.path.join(src_dir, qgc.DY_HERWIG_FILENAME))

        input_singlemu_tfile = cu.open_root_file(os.path.join(src_dir, qgc.SINGLE_MU_FILENAME))

        zpj_region_dict = {
            "name": "ZPlusJets",
            "dirname": "ZPlusJets_QG_Unfold",
            "label": "Z+jets",
            "data_tfile": input_singlemu_tfile,
            "mc_tfile": input_mc_dy_mgpythia_tfile,
            "alt_mc_tfile": input_mc_dy_mgherwig_tfile,
            "alt_mc_label": "MG+Herwig++",
            "tau_limits": {
                'jet_puppiMultiplicity': (1E-10, 1E-4),
                'jet_pTD': (1E-10, 1E-4),
                'jet_LHA': (1E-10, 1E-4),
                'jet_width': (1E-10, 1E-4),
                'jet_thrust': (1E-10, 1E-4),
                'jet_puppiMultiplicity_charged': (1E-10, 1E-4),
                'jet_pTD_charged': (1E-10, 1E-4),
                'jet_LHA_charged': (1E-10, 1E-4),
                'jet_width_charged': (1E-10, 1E-4),
                'jet_thrust_charged': (1E-10, 1E-4),
            },
            "systematics": [
                {
                    "label": "Neutral hadron up",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'neutralHadronShiftUp', qgc.DY_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kOrange-3,
                },
                {
                    "label": "Neutral hadron down",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'neutralHadronShiftDown', qgc.DY_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kOrange-3,
                    "linestyle": 2,
                },
                {
                    "label": "Photon up",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'photonShiftUp', qgc.DY_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kMagenta-3,
                },
                {
                    "label": "Photon down",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'photonShiftDown', qgc.DY_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kMagenta-3,
                    "linestyle": 2,
                },
                {
                    "label": "JEC up",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'jecsmear_directionUp', qgc.DY_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kGreen+2,
                },
                {
                    "label": "JEC down",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'jecsmear_directionDown', qgc.DY_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kGreen+2,
                    "linestyle": 2,
                },
                {
                    "label": "JER up",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'jersmear_directionUp', qgc.DY_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kOrange+3,
                },
                {
                    "label": "JER down",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'jersmear_directionDown', qgc.DY_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kOrange+3,
                    "linestyle": 2,
                },
                {
                    "label": "Pileup up",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'pileup_directionUp', qgc.DY_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kBlue-4,
                },
                {
                    "label": "Pileup down",
                    "tfile": cu.open_root_file(os.path.join(src_dir_systs, 'pileup_directionDown', qgc.DY_FILENAME)),
                    "systematic_type": "response",
                    "colour": ROOT.kBlue-4,
                    "linestyle": 2,
                },
            ]
        }

        if args.doZPJ:
            regions.append(zpj_region_dict)

        if args.doZPJGroomed:
            zpj_region_groomed_dict = zpj_region_dict.copy()
            zpj_region_groomed_dict['dirname'] = 'ZPlusJets_QG_Unfold_groomed'
            zpj_region_groomed_dict['name'] = 'ZPlusJets_groomed'
            zpj_region_groomed_dict['label'] = 'Z+jets'
            regions.append(zpj_region_groomed_dict)

    # Setup various options
    # --------------------------------------------------------------------------

    REGULARIZE = None
    REGULARIZE = "tau"
    # REGULARIZE = "L"
    REGULARIZE = args.regularize

    # Run with MC input instead of data
    MC_INPUT = args.MCinput
    mc_append = "_MC" if MC_INPUT else ""

    # If True, use part of MC for response matrix, and separate part for unfolding
    # as independent test
    # Should always be false for data
    MC_SPLIT = args.MCsplit if args.MCinput else False
    if MC_INPUT:
        mc_append += "_split" if MC_SPLIT else "_all"

    SUBTRACT_FAKES = True  # this should alwys be True
    sub_append = "_subFakes" if SUBTRACT_FAKES else ""

    append = ""

    # Add in systematics
    DO_SYSTS = args.doSysts
    if DO_SYSTS:
        append += "_syst"

    if args.useAltResponse:
        append += "_altResponse"

    output_dir = os.path.join(src_dir, "unfolding_better_regularise%s_target0p5%s%s_densityModeBinWidth_constraintArea%s_signalRegionOnly" % (str(REGULARIZE).capitalize(), mc_append, sub_append, append))
    if args.outputDir:
        output_dir = args.outputDir
    cu.check_dir_exists_create(output_dir)

    jet_algo = "AK4 PF PUPPI"
    if "ak8puppi" in src_dir:
        jet_algo = "AK8 PF PUPPI"

    # Save hists etc to ROOT file for access later
    output_tfile = ROOT.TFile("%s/unfolding_result.root" % (output_dir), "RECREATE")

    angles = [a for a in qgc.COMMON_VARS if a.var in args.angles]
    print(angles)

    print("Running TUnfold version", ROOT.TUnfold.GetTUnfoldVersion())

    # Do unfolding per signal region
    # --------------------------------------------------------------------------
    for region in regions[:]:

        # Setup pt bins
        # -------------
        # need different ones for Z+Jets region
        is_zpj = "ZPlusJets" in region['name']

        zpj_append = "_zpj" if is_zpj else ""

        pt_bin_edges_gen = qgc.PT_UNFOLD_DICT['signal%s_gen' % (zpj_append)]
        pt_bin_edges_reco = qgc.PT_UNFOLD_DICT['signal%s_reco' % (zpj_append)]

        pt_bin_edges_underflow_gen = qgc.PT_UNFOLD_DICT['underflow%s_gen' % (zpj_append)]
        pt_bin_edges_underflow_reco = qgc.PT_UNFOLD_DICT['underflow%s_reco' % (zpj_append)]

        new_tdir = region['name']
        output_tfile.mkdir(new_tdir)
        region_tdir = output_tfile.Get(new_tdir)
        region_tdir.cd()

        # Do 1D unfolding of pt
        # ----------------------------------------------------------------------
        append = "%s_pt" % (region['name'])  # common str to put on filenames, etc
        # print("*"*80)
        # print("Region/var: %s" % (append))
        # print("*"*80)

        # hist_data_reco = cu.get_from_tfile(region['data_tfile'], "%s/hist_pt_reco_all" % (region['dirname']))
        mc_hname_append = "split" if MC_SPLIT else "all"
        # hist_mc_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_pt_reco_%s" % (region['dirname'], mc_hname_append))
        # hist_mc_gen = cu.get_from_tfile(region['mc_tfile'], "%s/hist_pt_truth_%s" % (region['dirname'], mc_hname_append))
        hist_mc_gen_pt = cu.get_from_tfile(region['mc_tfile'], "%s/hist_pt_truth_%s" % (region['dirname'], mc_hname_append))
        # hist_mc_gen_reco_map = cu.get_from_tfile(region['mc_tfile'], "%s/tu_pt_GenReco_%s" % (region['dirname'], mc_hname_append))
        # TODO!

        # Remake gen hist with physical bins & save to file
        all_pt_bins_gen = np.concatenate((pt_bin_edges_underflow_gen[:-1], pt_bin_edges_gen))
        hist_mc_gen_pt_physical = ROOT.TH1F("mc_gen_pt", ";p_{T}^{jet} [GeV];N", len(all_pt_bins_gen)-1, array('d', all_pt_bins_gen))
        update_hist_bin_content(h_orig=hist_mc_gen_pt, h_to_be_updated=hist_mc_gen_pt_physical)
        update_hist_bin_error(h_orig=hist_mc_gen_pt, h_to_be_updated=hist_mc_gen_pt_physical)
        region_tdir.WriteTObject(hist_mc_gen_pt_physical, "mc_gen_pt")

        # Do unfolding for each angle
        # ----------------------------------------------------------------------
        for angle in angles:
            angle_prepend = "groomed " if "groomed" in region['name'] else ""
            append = "%s_%s%s" % (region['name'], angle_prepend.lower(), angle.var)  # common str to put on filenames, etc
            angle_str = "%s%s (%s)" % (angle_prepend, angle.name, angle.lambda_str)

            print("*"*80)
            print("Region/var: %s" % (append))
            print("*"*80)

            new_tdir = "%s/%s" % (region['name'], angle.var)
            output_tfile.mkdir(new_tdir)
            this_tdir = output_tfile.Get(new_tdir)
            this_tdir.cd()

            # put plots in subdir, to avoid overcrowding
            this_output_dir = "%s/%s/%s" % (output_dir, region['name'], angle.var)
            cu.check_dir_exists_create(this_output_dir)

            # Setup MyUnfolder object to do unfolding etc
            # -------------------------------------------
            angle_bin_edges_reco = qgc.VAR_UNFOLD_DICT[angle.var]['reco']
            angle_bin_edges_gen = qgc.VAR_UNFOLD_DICT[angle.var]['gen']
            angle_shortname = angle.var.replace("jet_", "")

            hist_data_reco = cu.get_from_tfile(region['data_tfile'], "%s/hist_%s_reco_all" % (region['dirname'], angle_shortname))
            mc_hname_append = "split" if MC_SPLIT else "all"
            hist_mc_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_%s" % (region['dirname'], angle_shortname, mc_hname_append))
            hist_mc_gen = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_truth_%s" % (region['dirname'], angle_shortname, mc_hname_append))
            hist_mc_gen_reco_map = cu.get_from_tfile(region['mc_tfile'], "%s/tu_%s_GenReco_%s" % (region['dirname'], angle_shortname, mc_hname_append))

            # Actual distribution to be unfolded
            reco_1d = hist_mc_reco.Clone() if MC_INPUT else hist_data_reco

            if SUBTRACT_FAKES:
                # to construct our "fakes" template, we use the ratio as predicted by MC, and apply it to data
                # this way we ensure we don't have -ve values, and avoid any issue with cross sections
                hist_mc_fakes_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_fake_%s" % (region['dirname'], angle_shortname, mc_hname_append))
                hist_fakes_reco = hist_mc_fakes_reco.Clone("hist_%s_fakes" % angle_shortname)
                hist_fakes_reco.Divide(hist_mc_reco)

                # plot fake fraction before multiplting by 'data'
                plot_fake_fraction(hist_fakes_reco,
                                   output_filename="%s/fake_fraction_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                                   title="%s region, %s" % (region['label'], angle_str))

                hist_fakes_reco.Multiply(reco_1d)

                # background-subtracted reco hists, only for plotting purposes, not for TUnfold (that does background subtraction internally)
                reco_1d_bg_subtracted = reco_1d.Clone()
                reco_1d_bg_subtracted.Add(hist_fakes_reco, -1)

                chosen_bin = 15
                print("1D reco input without background subtraction:", reco_1d.GetBinContent(chosen_bin))
                print("1D reco input with background subtraction:", reco_1d_bg_subtracted.GetBinContent(chosen_bin))
                print("1D reco input fakes:", hist_fakes_reco.GetBinContent(chosen_bin))

                hist_data_reco_bg_subtracted = hist_data_reco.Clone(hist_data_reco.GetName() + "_bgrSubtracted")
                hist_data_reco_bg_subtracted.Add(hist_fakes_reco, -1)

                hist_mc_reco_bg_subtracted = hist_mc_reco.Clone(hist_mc_reco.GetName() + "_bgrSubtracted")
                hist_mc_reco_bg_subtracted.Add(hist_mc_fakes_reco, -1)  # should this be hist_fakes_reco? depends on what we want to see...

            # hist_mc_gen_reco_neutralUp_map = cu.get_from_tfile(region['mc_neutralUp_tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
            # hist_mc_gen_reco_neutralDown_map = cu.get_from_tfile(region['mc_neutralDown_tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))

            mc_hname_append = "_split" if MC_SPLIT else ""  # FIXME consistency in unfold hist module!
            hist_data_reco_gen_binning = cu.get_from_tfile(region['data_tfile'], "%s/hist_%s_reco_gen_binning" % (region['dirname'], angle_shortname))
            hist_mc_reco_gen_binning = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_gen_binning%s" % (region['dirname'], angle_shortname, mc_hname_append))

            # Actual distribution to be unfolded, but with gen binning
            reco_1d_gen_binning = hist_mc_reco_gen_binning.Clone() if MC_INPUT else hist_data_reco_gen_binning
            # reco_1d_gen_binning = hist_mc_reco_gen_binning # only for testing that everything is setup OK

            if SUBTRACT_FAKES:
                mc_hname_append = "_split" if MC_SPLIT else ""  # FIXME consistency in unfold hist module!
                hist_mc_fakes_reco_gen_binning = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_fake_gen_binning%s" % (region['dirname'], angle_shortname, mc_hname_append))
                # create template as above, but with gen binning
                hist_fakes_reco_gen_binning = hist_mc_fakes_reco_gen_binning.Clone("hist_%s_fakes_gen_binning" % angle_shortname)
                hist_fakes_reco_gen_binning.Divide(hist_mc_reco_gen_binning)
                hist_fakes_reco_gen_binning.Multiply(reco_1d_gen_binning)

                # background-subtracted reco hists, only for plotting purposes, not for TUnfold (that does background subtraction internally)
                reco_1d_gen_binning_bg_subtracted = reco_1d_gen_binning.Clone()
                reco_1d_gen_binning_bg_subtracted.Add(hist_fakes_reco_gen_binning, -1)

                hist_data_reco_gen_binning_bg_subtracted = hist_data_reco_gen_binning.Clone(hist_data_reco_gen_binning.GetName() + "_bgrSubtracted")
                hist_data_reco_gen_binning_bg_subtracted.Add(hist_fakes_reco_gen_binning, -1)

                hist_mc_reco_gen_binning_bg_subtracted = hist_mc_reco_gen_binning.Clone(hist_mc_reco_gen_binning.GetName() + "_bgrSubtracted")
                hist_mc_reco_gen_binning_bg_subtracted.Add(hist_mc_fakes_reco_gen_binning, -1)  # should this be hist_fakes_reco_gen_binning? depends on what we want to see...

            # Setup unfolder object
            # ---------------------
            variable_name = "%s%s" % (angle_prepend, angle.name)
            unfolder = MyUnfolder(response_map=hist_mc_gen_reco_map,
                                  variable_bin_edges_reco=angle_bin_edges_reco,
                                  variable_bin_edges_gen=angle_bin_edges_gen,
                                  variable_name=variable_name,
                                  pt_bin_edges_reco=pt_bin_edges_reco,
                                  pt_bin_edges_gen=pt_bin_edges_gen,
                                  pt_bin_edges_underflow_reco=pt_bin_edges_underflow_reco,
                                  pt_bin_edges_underflow_gen=pt_bin_edges_underflow_gen,
                                  orientation=ROOT.TUnfold.kHistMapOutputHoriz,
                                  constraintMode=ROOT.TUnfold.kEConstraintArea,
                                  # constraintMode=ROOT.TUnfold.kEConstraintNone,
                                  regMode=ROOT.TUnfold.kRegModeCurvature,
                                  densityFlags=ROOT.TUnfoldDensity.kDensityModeBinWidth,
                                  # densityFlags=ROOT.TUnfoldDensity.kDensityModeNone,
                                  axisSteering='*[B]')

            unfolder.save_binning(txt_filename="%s/binning_scheme.txt" % (this_output_dir), print_xml=False)

            # Subtract fakes (treat as background)
            # ------------------------------------
            if SUBTRACT_FAKES:
                unfolder.tunfolder.SubtractBackground(hist_fakes_reco, "fakes")

            # Add systematic errors as different response matrices
            # ----------------------------------------------------
            if DO_SYSTS:
                chosen_rsp_bin = (18, 18)
                print("nominal response bin content for", chosen_rsp_bin, hist_mc_gen_reco_map.GetBinContent(*chosen_rsp_bin))
                for syst_dict in region['systematics']:
                    if 'tfile' in syst_dict and syst_dict['systematic_type'] == 'response':
                        map_syst = cu.get_from_tfile(syst_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                        print("Adding systematic:", syst_dict['label'])
                        print("    syst bin", chosen_rsp_bin, map_syst.GetBinContent(*chosen_rsp_bin))
                        unfolder.tunfolder.AddSysError(map_syst, syst_dict['label'], unfolder.orientation, ROOT.TUnfoldDensity.kSysErrModeMatrix)

            # Set what is to be unfolded
            # ---------------------
            unfolder.setInput(reco_1d)

            # Save important stuff to TFile
            # ------------------------------------------------------------------
            this_tdir.WriteTObject(unfolder.detector_binning)
            this_tdir.WriteTObject(unfolder.generator_binning)
            this_tdir.WriteTObject(reco_1d, "unfold_input")
            this_tdir.WriteTObject(reco_1d_gen_binning, "unfold_input_gen_binning")
            this_tdir.WriteTObject(hist_mc_reco, "mc_reco")
            this_tdir.WriteTObject(hist_mc_reco_gen_binning, "mc_reco_gen_binning")
            this_tdir.WriteTObject(hist_mc_gen, "mc_gen")
            this_tdir.WriteTObject(hist_mc_gen_reco_map, "response_map")  # response map
            # this_tdir.WriteTObject(hist_mc_gen_reco_map.RebinY(2), "response_map_rebin")  # response map

            if SUBTRACT_FAKES:
                this_tdir.WriteTObject(hist_fakes_reco, "fakes")
                this_tdir.WriteTObject(hist_fakes_reco_gen_binning, "fakes_gen_binning")
                this_tdir.WriteTObject(reco_1d_bg_subtracted, "unfold_input_bg_subtracted")
                this_tdir.WriteTObject(reco_1d_gen_binning_bg_subtracted, "unfold_input_gen_binning_bg_subtracted")
                this_tdir.WriteTObject(hist_mc_reco_bg_subtracted, "mc_reco_bg_subtracted")
                this_tdir.WriteTObject(hist_mc_reco_gen_binning_bg_subtracted, "mc_reco_gen_binning_bg_subtracted")

            # Do any regularization
            # ---------------------
            # tau = 1E-10
            tau = 0
            n_scan = 100
            scan_mode = ROOT.TUnfoldDensity.kEScanTauRhoAvgSys
            scan_distribution = "generatordistribution"
            if REGULARIZE == "L":
                print("Regularizing with ScanL, please be patient...")
                tau = unfolder.doScanL(output_dir=this_output_dir, n_scan=n_scan,
                                       tau_min=1E-14, tau_max=1E-4)
                print("Found tau:", tau)
            elif REGULARIZE == "tau":
                print("Regularizing with ScanTau, please be patient...")
                tau_scanner = TauScanner()
                tau = tau_scanner.scan_tau(tunfolder=unfolder.tunfolder,
                                           n_scan=n_scan,
                                           tau_min=region['tau_limits'][angle.var][0],
                                           tau_max=region['tau_limits'][angle.var][1],
                                           scan_mode=scan_mode,
                                           distribution=scan_distribution,
                                           axis_steering=unfolder.axisSteering)
                print("Found tau:", tau)
                tau_scanner.plot_scan_tau(output_filename="%s/scantau_%s.%s" % (this_output_dir, unfolder.variable_name, OUTPUT_FMT))
                this_tdir.WriteTObject(tau_scanner.graph_all_scan_points, "regularize_all_scan_points")
                this_tdir.WriteTObject(tau_scanner.graph_best_scan_point, "regularize_best_scan_point")

            # Do unfolding!
            # ---------------------
            unfolder.do_unfolding(tau)
            unfolded_1d = unfolder.get_output()
            unfolded_1d.SetName("unfolded_1d")
            chosen_bin = 18
            print("Bin %d:" % chosen_bin, unfolded_1d.GetBinContent(chosen_bin))
            print("original uncert:", unfolded_1d.GetBinError(chosen_bin))

            # Get various error matrices
            # ------------------------------------------------------------------
            # stat errors only - do before or after systematics?
            ematrix_input = unfolder.get_ematrix_input() # stat errors from input to be unfolded
            print("ematrix_input uncert:", ematrix_input.GetBinContent(chosen_bin, chosen_bin))
            this_tdir.WriteTObject(ematrix_input, "ematrix_input")

            ematrix_sys_uncorr = unfolder.get_ematrix_sys_uncorr() # stat errors in response matrix
            print("ematrix_sys_uncorr uncert:", ematrix_sys_uncorr.GetBinContent(chosen_bin, chosen_bin))
            this_tdir.WriteTObject(ematrix_sys_uncorr, "ematrix_sys_uncorr")

            ematrix_stat_sum = ematrix_input.Clone("ematrix_stat_sum")
            ematrix_stat_sum.Add(ematrix_sys_uncorr) # total 'stat' errors

            error_stat_1d = make_hist_from_diagonal_errors(ematrix_stat_sum, do_sqrt=True) # note that bin contents = 0, only bin errors are non-0
            this_tdir.WriteTObject(error_stat_1d, "error_stat_1d")
            print("stat uncert:", error_stat_1d.GetBinError(chosen_bin))

            ematrix_total = unfolder.get_ematrix_total()
            error_total_1d = make_hist_from_diagonal_errors(ematrix_total, do_sqrt=True) # note that bin contents = 0, only bin errors are non-0
            this_tdir.WriteTObject(ematrix_total, "ematrix_total_1d")
            this_tdir.WriteTObject(error_total_1d, "error_total_1d")
            print("total uncert:", error_total_1d.GetBinError(chosen_bin))

            # Update errors to big unfolded 1D
            update_hist_bin_error(h_orig=error_total_1d, h_to_be_updated=unfolded_1d)
            print("new uncert:", unfolded_1d.GetBinError(chosen_bin))
            this_tdir.WriteTObject(unfolded_1d)

            # Get shifts due to systematics
            # ------------------------------------------------------------------
            systematic_shift_hists = []
            if DO_SYSTS:
                for syst_dict in region['systematics']:
                    if syst_dict['systematic_type'] == 'response':
                        syst_label = syst_dict['label']
                        h_syst = unfolder.tunfolder.GetDeltaSysSource(syst_label,
                                                                      '%s_%s_shift_%s' % (region['name'], angle.var, syst_label),
                                                                      "",
                                                                      "generator",  # must be the same as what's used in get_output
                                                                      unfolder.axisSteering,
                                                                      unfolder.use_axis_binning)
                        systematic_shift_hists.append(h_syst)
                        this_tdir.WriteTObject(h_syst)

            # Draw unified unfolded distributions
            # ------------------------------------------------------------------
            # unfolded, gen, and reco for comparison
            plot_simple_unfolded(unfolded=unfolded_1d,
                                 tau=tau,
                                 reco=None,
                                 gen=hist_mc_gen,
                                 fake=None,
                                 output_filename="%s/unfolded_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                                 title="%s region, %s" % (region['label'], angle_str))

            # reco using detector binning
            plot_simple_detector(reco_data=reco_1d,
                                 reco_mc=hist_mc_reco,
                                 reco_mc_fake=hist_mc_fakes_reco if SUBTRACT_FAKES else None,
                                 reco_data_fake=hist_fakes_reco if SUBTRACT_FAKES else None,
                                 output_filename="%s/detector_reco_binning_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                                 title="%s region, %s" % (region['label'], angle_str))

            if SUBTRACT_FAKES:
                # same plot but with background-subtracted reco
                plot_simple_detector(reco_data=reco_1d_bg_subtracted,
                                     reco_mc=hist_mc_reco_bg_subtracted,
                                     reco_mc_fake=hist_mc_fakes_reco,
                                     reco_data_fake=hist_fakes_reco,
                                     output_filename="%s/detector_reco_binning_bg_subtracted_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                                     title="%s region, %s" % (region['label'], angle_str))

            # reco using gen binning
            plot_simple_detector(reco_data=reco_1d_gen_binning,
                                 reco_mc=hist_mc_reco_gen_binning,
                                 reco_mc_fake=hist_mc_fakes_reco_gen_binning if SUBTRACT_FAKES else None,
                                 reco_data_fake=hist_fakes_reco_gen_binning if SUBTRACT_FAKES else None,
                                 output_filename="%s/detector_gen_binning_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                                 title="%s region, %s" % (region['label'], angle_str))

            if SUBTRACT_FAKES:
                # same but with background-subtracted
                plot_simple_detector(reco_data=reco_1d_gen_binning_bg_subtracted,
                                     reco_mc=hist_mc_reco_gen_binning_bg_subtracted,
                                     reco_mc_fake=hist_mc_fakes_reco_gen_binning,
                                     reco_data_fake=hist_fakes_reco_gen_binning,
                                     output_filename="%s/detector_gen_binning_bg_subtracted_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                                     title="%s region, %s" % (region['label'], angle_str))

            # Draw collapsed distributions
            # ---------------------
            # TODO!

            # Draw projections of response matrix vs 1D hist to check normalisation OK
            # Only makes sense if the same MC events go into matrix & 1D plot
            # ------------------------------------------------------------------
            if not MC_SPLIT:
                proj_reco = hist_mc_gen_reco_map.ProjectionY("proj_reco_%s" % (append))

                proj_gen = hist_mc_gen_reco_map.ProjectionX("proj_gen_%s" % (append))
                draw_projection_comparison(hist_mc_gen, proj_gen,
                                           title="%s\n%s region" % (jet_algo, region['label']),
                                           xtitle="%s, Generator binning" % (angle_str),
                                           output_filename="%s/projection_gen_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

                print("projection reco #bins:", proj_reco.GetNbinsX())
                print("response map # bins x:", hist_mc_gen_reco_map.GetNbinsX())
                print("response map # bins y:", hist_mc_gen_reco_map.GetNbinsY())
                if SUBTRACT_FAKES:
                    print("reco bg subtracted #bins:", hist_mc_reco_bg_subtracted.GetNbinsX())
                    print(proj_reco.GetNbinsX())
                    # Do the same but with backgrounds subtracted from the 1D
                    draw_projection_comparison(hist_mc_reco_bg_subtracted, proj_reco,
                                               title="%s\n%s region" % (jet_algo, region['label']),
                                               xtitle="%s, Detector binning" % (angle_str),
                                               output_filename="%s/projection_reco_bg_subtracted_%s.%s" % (this_output_dir, append, OUTPUT_FMT))
                else:
                    draw_projection_comparison(hist_mc_reco, proj_reco,
                           title="%s\n%s region" % (jet_algo, region['label']),
                           xtitle="%s, Detector binning" % (angle_str),
                           output_filename="%s/projection_reco_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                           print_bin_comparison=False)

            # Draw matrices
            # ---------------------
            draw_response_matrix(unfolder.response_map,
                                 region['label'],
                                 angle_str,
                                 "%s/response_map_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

            prob_map = unfolder.get_probability_matrix()
            this_tdir.WriteTObject(prob_map, "prob_matrix")
            draw_probability_matrix(prob_map,
                                    region['label'],
                                    angle_str,
                                    "%s/probability_map_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

            rhoij = unfolder.get_rhoij_total()
            this_tdir.WriteTObject(rhoij, "rhoij")
            draw_correlation_matrix(rhoij,
                                    region['label'],
                                    angle_str,
                                    "%s/rho_map_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

            draw_error_matrix_input(ematrix_input,
                                    region['label'],
                                    angle_str,
                                    "%s/err_map_sys_input_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

            draw_error_matrix_sys_uncorr(ematrix_sys_uncorr,
                                         region['label'],
                                         angle_str,
                                         "%s/err_map_sys_uncorr_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

            draw_error_matrix_total(ematrix_total,
                                    region['label'],
                                    angle_str,
                                    "%s/err_map_total_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

            # Do forward-folding to check unfolding
            # ------------------------------------------------------------------
            hist_data_folded = unfolder.get_folded_output(hist_name="folded_%s" % (append))
            this_tdir.WriteTObject(hist_data_folded, "folded_1d")
            draw_reco_folded(hist_folded=hist_data_folded,
                             tau=tau,
                             hist_reco_data=reco_1d_bg_subtracted if SUBTRACT_FAKES else reco_1d,
                             hist_reco_mc=hist_mc_reco_bg_subtracted if SUBTRACT_FAKES else hist_mc_reco,
                             title="%s\n%s region" % (jet_algo, region['label']),
                             xtitle="%s, Detector binning" % (angle_str),
                             output_filename="%s/folded_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

            summary_1d_entries = []  # for final summary plot

            # Do unfolding again, this time using alternate response matrix
            # ------------------------------------------------------------------
            alt_unfolded_1d = None
            alt_hist_mc_gen = None
            alt_tau = 0
            if args.useAltResponse:
                print("*** Unfolding with alternate response matrix ***")
                alt_hist_mc_gen = cu.get_from_tfile(region['alt_mc_tfile'], "%s/hist_%s_truth_all" % (region['dirname'], angle_shortname))
                hist_mc_gen_reco_map_alt = cu.get_from_tfile(region['alt_mc_tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                alt_unfolder = MyUnfolder(response_map=hist_mc_gen_reco_map_alt,
                                          variable_bin_edges_reco=unfolder.variable_bin_edges_reco,
                                          variable_bin_edges_gen=unfolder.variable_bin_edges_gen,
                                          variable_name=unfolder.variable_name,
                                          pt_bin_edges_reco=unfolder.pt_bin_edges_reco,
                                          pt_bin_edges_gen=unfolder.pt_bin_edges_gen,
                                          pt_bin_edges_underflow_reco=unfolder.pt_bin_edges_underflow_reco,
                                          pt_bin_edges_underflow_gen=unfolder.pt_bin_edges_underflow_gen,
                                          orientation=unfolder.orientation,
                                          constraintMode=unfolder.constraintMode,
                                          regMode=unfolder.regMode,
                                          densityFlags=unfolder.densityFlags,
                                          axisSteering=unfolder.axisSteering)

                # Subtract fakes (treat as background)
                # --------------------------------------------------------------
                if SUBTRACT_FAKES:
                    alt_unfolder.tunfolder.SubtractBackground(hist_fakes_reco, "fakes")

                # Set what is to be unfolded
                # --------------------------------------------------------------
                alt_unfolder.setInput(reco_1d)

                # Save important stuff to TFile
                # --------------------------------------------------------------
                this_tdir.WriteTObject(hist_mc_gen_reco_map_alt, "alt_response_map")  # response map

                # Do any regularization
                # --------------------------------------------------------------
                alt_tau = 0
                if REGULARIZE == "L":
                    print("Regularizing alternative with ScanL, please be patient...")
                    alt_tau = alt_unfolder.doScanL(output_dir=this_output_dir, n_scan=100,
                                           tau_min=1E-14, tau_max=1E-4)
                    print("Found tau:", tau)
                elif REGULARIZE == "tau":
                    print("Regularizing alternative with ScanTau, please be patient...")
                    alt_tau_scanner = TauScanner()
                    alt_tau = alt_tau_scanner.scan_tau(tunfolder=alt_unfolder.tunfolder,
                                                       n_scan=n_scan,
                                                       tau_min=region['tau_limits'][angle.var][0],
                                                       tau_max=region['tau_limits'][angle.var][1],
                                                       scan_mode=scan_mode,
                                                       distribution=scan_distribution,
                                                       axis_steering=alt_unfolder.axisSteering)
                    print("Found tau for alt matrix:", alt_tau)
                    alt_tau_scanner.plot_scan_tau(output_filename="%s/scantau_alt_%s.%s" % (this_output_dir, alt_unfolder.variable_name, OUTPUT_FMT))

                # Do unfolding!
                # --------------------------------------------------------------
                alt_unfolder.do_unfolding(alt_tau)
                alt_unfolded_1d = alt_unfolder.get_output()
                alt_unfolded_1d.SetName("alt_unfolded_1d")
                print("Bin %d:" % chosen_bin, alt_unfolded_1d.GetBinContent(chosen_bin))
                print("original uncert:", alt_unfolded_1d.GetBinError(chosen_bin))

                # Get error matrix to update errors
                # --------------------------------------------------------------
                alt_ematrix_total = alt_unfolder.get_ematrix_total()
                alt_error_total_1d = make_hist_from_diagonal_errors(ematrix_total, do_sqrt=True) # note that bin contents = 0, only bin errors are non-0
                this_tdir.WriteTObject(alt_ematrix_total, "alt_ematrix_total_1d")
                this_tdir.WriteTObject(alt_error_total_1d, "alt_error_total_1d")
                print("total uncert:", alt_error_total_1d.GetBinError(chosen_bin))

                # Update errors to big unfolded 1D
                update_hist_bin_error(h_orig=alt_error_total_1d, h_to_be_updated=alt_unfolded_1d)
                print("new uncert:", alt_unfolded_1d.GetBinError(chosen_bin))
                this_tdir.WriteTObject(alt_unfolded_1d)

            # ------------------------------------------------------------------
            # PLOTTING LOTS OF THINGS
            # ------------------------------------------------------------------

            # Some common plotting vars
            # ------------------------------------------------------------------
            detector_title = "Detector-level " + angle_str
            particle_title = "Particle-level " + angle_str
            normalised_differential_label = "#frac{1}{#sigma} #frac{d#sigma}{d%s}" % angle.lambda_str

            # Draw individual pt bin plots - GEN binning
            # ------------------------------------------------------------------
            for ibin_pt in range(0, len(pt_bin_edges_gen)-1):

                this_pt_bin_tdir = this_tdir.mkdir("gen_bin_%d" % (ibin_pt))
                print("Individual gen bin", ibin_pt, "=", pt_bin_edges_gen[ibin_pt], pt_bin_edges_gen[ibin_pt+1])

                # Produce 1D hists for this pt bin
                # --------------------------------------------------------------
                # Unfolded hists
                mc_gen_hist_bin = unfolder.get_var_hist_pt_binned(hist_mc_gen, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(mc_gen_hist_bin, "mc_gen_hist_bin")

                unfolded_hist_bin = unfolder.get_var_hist_pt_binned(unfolded_1d, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(unfolded_hist_bin, "unfolded_hist_bin")

                unfolded_hist_bin_stat_errors = unfolder.get_var_hist_pt_binned(error_stat_1d, ibin_pt, binning_scheme="generator")
                update_hist_bin_content(h_orig=unfolded_hist_bin, h_to_be_updated=unfolded_hist_bin_stat_errors)  # use stat errors, update central values
                this_pt_bin_tdir.WriteTObject(unfolded_hist_bin_stat_errors, "unfolded_hist_bin_stat_errors")

                unfolded_hist_bin_total_errors = unfolder.get_var_hist_pt_binned(error_total_1d, ibin_pt, binning_scheme="generator")
                update_hist_bin_content(h_orig=unfolded_hist_bin, h_to_be_updated=unfolded_hist_bin_total_errors)  # use total errors, update central values
                this_pt_bin_tdir.WriteTObject(unfolded_hist_bin_total_errors, "unfolded_hist_bin_total_errors") # same as unfolded_1d?

                # Reco level but with gen binning
                # For 'data'
                reco_hist_bin_gen_binning = unfolder.get_var_hist_pt_binned(reco_1d_gen_binning, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(reco_hist_bin_gen_binning, "reco_hist_bin_gen_binning")

                if SUBTRACT_FAKES:
                    reco_hist_bg_subtracted_bin_gen_binning = unfolder.get_var_hist_pt_binned(reco_1d_gen_binning_bg_subtracted, ibin_pt, binning_scheme="generator")
                    this_pt_bin_tdir.WriteTObject(reco_hist_bg_subtracted_bin_gen_binning, "reco_hist_bg_subtracted_bin_gen_binning")

                # For MC
                mc_reco_hist_bin_gen_binning = unfolder.get_var_hist_pt_binned(hist_mc_reco_gen_binning, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(mc_reco_hist_bin_gen_binning, "mc_reco_hist_bin_gen_binning")

                if SUBTRACT_FAKES:
                    mc_reco_hist_bg_subtracted_bin_gen_binning = unfolder.get_var_hist_pt_binned(hist_mc_reco_gen_binning_bg_subtracted, ibin_pt, binning_scheme="generator")
                    this_pt_bin_tdir.WriteTObject(mc_reco_hist_bg_subtracted_bin_gen_binning, "mc_reco_hist_bg_subtracted_bin_gen_binning")


                # Make lots of plots
                # ------------------------------------------------------------
                lw = 2
                # common hist settings
                title = "%s\n%s region\n%g < p_{T}^{jet} < %g GeV" % (jet_algo, region['label'], pt_bin_edges_gen[ibin_pt], pt_bin_edges_gen[ibin_pt+1])
                common_hist_args = dict(
                    what="hist",
                    title=title,
                    subplot_type='ratio',
                    subplot_limits=(0.5, 1.5),
                )

                # PLOT UNFOLDED DATA
                # --------------------------------------------------------------
                gen_colour = ROOT.kRed
                unfolded_basic_colour = ROOT.kAzure+7
                unfolded_stat_colour = ROOT.kGreen+2
                unfolded_total_colour = ROOT.kBlack

                # unnormalised version
                entries = [
                    Contribution(mc_gen_hist_bin,
                                 label="Generator",
                                 line_color=gen_colour, line_width=lw,
                                 marker_color=gen_colour, marker_size=0,
                                 normalise_hist=False),
                    Contribution(unfolded_hist_bin_total_errors,
                                 label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                 line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_gen_hist_bin,
                                 normalise_hist=False),
                    Contribution(unfolded_hist_bin_stat_errors,
                                 label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                 line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_stat_colour, marker_size=0,
                                 subplot=mc_gen_hist_bin,
                                 normalise_hist=False),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=particle_title,
                            ytitle="N",
                            subplot_title='Unfolded / gen',
                            **common_hist_args
                            )
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.68)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/unfolded_unnormalised_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # now normalise each plot to unity
                # Note that this modifies e.g. mc_gen_hist_bin, so from this point
                # onwards it will be normalised to unity
                entries = [
                    Contribution(mc_gen_hist_bin,
                                 label="Generator (MG+Pythia8)",
                                 line_color=gen_colour, line_width=lw,
                                 marker_color=gen_colour, marker_size=0,
                                 subplot=unfolded_hist_bin_total_errors,
                                 normalise_hist=True),
                    Contribution(unfolded_hist_bin_total_errors,
                                 label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                 line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                 normalise_hist=True),
                    Contribution(unfolded_hist_bin_stat_errors,
                                 label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                 line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_stat_colour, marker_size=0,
                                 normalise_hist=True),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=particle_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Gen / Unfolded',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.68)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/unfolded_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Do a version where divided by bin width
                # Note that these hists are already normalised to 1!
                # Do not use normalise_hist!
                mc_gen_hist_bin_div_bin_width = qgp.hist_divide_bin_width(mc_gen_hist_bin)
                unfolded_hist_bin_stat_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_stat_errors)
                unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_total_errors)
                entries = [
                    Contribution(mc_gen_hist_bin_div_bin_width,
                                 label="Generator (MG+Pythia8)",
                                 line_color=gen_colour, line_width=lw,
                                 marker_color=gen_colour, marker_size=0,
                                 subplot=unfolded_hist_bin_total_errors_div_bin_width,
                                 normalise_hist=False),
                    Contribution(unfolded_hist_bin_total_errors_div_bin_width,
                                 label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                 line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                 normalise_hist=False),
                    Contribution(unfolded_hist_bin_stat_errors_div_bin_width,
                                 label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                 line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                 marker_color=unfolded_stat_colour, marker_size=0,
                                 normalise_hist=False),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=particle_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Gen / Unfolded',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.68)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/unfolded_%s_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                summary_1d_entries.append([
                    (unfolded_hist_bin_total_errors_div_bin_width,
                        dict(label="Unfolded\n($\\tau = %.3g$)" % (tau),
                             line_color=unfolded_total_colour, line_width=lw, line_style=3,
                             marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                             subplot=mc_gen_hist_bin,
                             normalise_hist=True)),  # unfolded data
                    (mc_gen_hist_bin_div_bin_width,
                        dict(label="MG+Pythia8",
                             line_color=gen_colour, line_width=lw,
                             marker_color=gen_colour, marker_size=0,
                             normalise_hist=True)),  # generator
                ])

                # Unfolded plots with alternate response matrix results as well
                # --------------------------------------------------------------
                # And alternate MC gen level to compare
                # (i.e. is the difference between unfolded results < diff at gen level?)
                if args.useAltResponse:
                    alt_mc_gen_hist_bin = alt_unfolder.get_var_hist_pt_binned(alt_hist_mc_gen, ibin_pt, binning_scheme="generator")
                    alt_unfolded_hist_bin_total_errors = alt_unfolder.get_var_hist_pt_binned(alt_unfolded_1d, ibin_pt, binning_scheme="generator")
                    this_pt_bin_tdir.WriteTObject(alt_unfolded_hist_bin_total_errors, "alt_unfolded_hist_bin_total_errors")

                    alt_gen_colour = ROOT.kBlack
                    alt_colour = ROOT.kBlue

                    entries = [
                        Contribution(mc_gen_hist_bin,
                                     label="Generator (%s)" % (region['mc_label']),
                                     line_color=gen_colour, line_width=lw,
                                     marker_color=gen_colour, marker_size=0,
                                     # subplot=unfolded_hist_bin_total_errors,
                                     normalise_hist=True),
                        Contribution(alt_mc_gen_hist_bin,
                                     label="Generator (%s)" % (region['alt_mc_label']),
                                     line_color=alt_gen_colour, line_width=lw,
                                     marker_color=alt_gen_colour, marker_size=0,
                                     # subplot=unfolded_hist_bin_total_errors,
                                     normalise_hist=True),
                        # Contribution(unfolded_hist_bin_total_errors,
                        #              label="Unfolded (#tau = %.3g) (total err)" % (tau),
                        #              line_color=unfolded_total_colour, line_width=lw, line_style=1,
                        #              marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                        #              subplot=mc_gen_hist_bin,
                        #              normalise_hist=True),
                        Contribution(unfolded_hist_bin_stat_errors,
                                     label="Unfolded (#tau = %.3g) (stat err)\n(%s response matrix)" % (tau, region['mc_label']),
                                     line_color=unfolded_stat_colour, line_width=lw, line_style=2,
                                     marker_color=unfolded_stat_colour, marker_size=0,
                                     subplot=mc_gen_hist_bin,
                                     normalise_hist=True),
                        Contribution(alt_unfolded_hist_bin_total_errors,
                                     label="Unfolded (#tau = %.3g) (stat err)\n(%s response matrix)" % (alt_tau, region['alt_mc_label']),
                                     line_color=alt_colour, line_width=lw, line_style=3,
                                     marker_color=alt_colour, marker_size=0,
                                     subplot=mc_gen_hist_bin,
                                     normalise_hist=True),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=particle_title,
                                ytitle=normalised_differential_label,
                                subplot_title='#splitline{Unfolded / Gen}{(%s)}' % (region['mc_label']),
                                **common_hist_args)
                    plot.legend.SetX1(0.55)
                    plot.legend.SetY1(0.72)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.88)
                    plot.plot("NOSTACK E1")
                    plot.save("%s/unfolded_%s_alt_response_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                    # Do a version where divided by bin width
                    # Note that inputs are already normalised to 1
                    # Do not use normalise_hist!
                    alt_mc_gen_hist_bin_div_bin_width = qgp.hist_divide_bin_width(alt_mc_gen_hist_bin)
                    alt_unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(alt_unfolded_hist_bin_total_errors)
                    unfolded_hist_bin_stat_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_stat_errors)
                    # unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_total_errors)

                    entries = [
                        Contribution(mc_gen_hist_bin_div_bin_width,
                                     label="Generator (%s)" % (region['mc_label']),
                                     line_color=gen_colour, line_width=lw,
                                     marker_color=gen_colour, marker_size=0,
                                     # subplot=unfolded_hist_bin_total_errors_div_bin_width,
                                     normalise_hist=False),
                        Contribution(alt_mc_gen_hist_bin_div_bin_width,
                                     label="Generator (%s)" % (region['alt_mc_label']),
                                     line_color=alt_gen_colour, line_width=lw,
                                     marker_color=alt_gen_colour, marker_size=0,
                                     # subplot=unfolded_hist_bin_total_errors_div_bin_width,
                                     normalise_hist=False),
                        # Contribution(unfolded_hist_bin_total_errors_div_bin_width,
                        #              label="Unfolded (#tau = %.3g) (total err)" % (tau),
                        #              line_color=unfolded_total_colour, line_width=lw, line_style=1,
                        #              marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                        #              subplot=mc_gen_hist_bin_div_bin_width,
                        #              normalise_hist=True),
                        Contribution(unfolded_hist_bin_stat_errors_div_bin_width,
                                     label="Unfolded (#tau = %.3g) (stat err)\n(%s response matrix)" % (tau, region['mc_label']),
                                     line_color=unfolded_stat_colour, line_width=lw, line_style=2,
                                     marker_color=unfolded_stat_colour, marker_size=0,
                                     subplot=mc_gen_hist_bin_div_bin_width,
                                     normalise_hist=False),
                        Contribution(alt_unfolded_hist_bin_total_errors_div_bin_width,
                                     label="Unfolded (#tau = %.3g) (stat err)\n(%s response matrix)" % (alt_tau, region['alt_mc_label']),
                                     line_color=alt_colour, line_width=lw, line_style=3,
                                     marker_color=alt_colour, marker_size=0,
                                     subplot=mc_gen_hist_bin_div_bin_width,
                                     normalise_hist=False),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=particle_title,
                                ytitle=normalised_differential_label,
                                subplot_title='#splitline{Unfolded / Gen}{(%s)}' % (region['mc_label']),
                                **common_hist_args)
                    plot.legend.SetX1(0.55)
                    plot.legend.SetY1(0.72)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.88)
                    plot.plot("NOSTACK E1")
                    plot.save("%s/unfolded_%s_alt_response_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Unfolded plots with variations in input model systematics plotted
                # --------------------------------------------------------------
                if DO_SYSTS:
                    syst_entries = []
                    syst_entries_div_bin_width = []
                    for syst_dict in region['systematics']:
                        if syst_dict['systematic_type'] != 'model':
                            continue
                        syst_label = syst_dict['label']
                        syst_label_no_spaces = syst_dict['label'].replace(" ", "_")
                        syst_tau = syst_dict['tau']
                        syst_unfolded_1d = syst_dict['unfolded_1d']
                        syst_unfolded_hist_bin_total_errors = unfolder.get_var_hist_pt_binned(syst_unfolded_1d, ibin_pt, binning_scheme="generator")
                        this_pt_bin_tdir.WriteTObject(syst_unfolded_hist_bin_total_errors, "syst_%s_unfolded_hist_bin_total_errors" % (syst_label_no_spaces))

                        syst_entries.append(
                            Contribution(syst_unfolded_hist_bin_total_errors,
                                         label="Unfolded (#tau = %.3g) (total err)\n(%s)" % (syst_tau, syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=1,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         subplot=mc_gen_hist_bin,
                                         normalise_hist=True),
                        )
                        # already normalised to 1
                        # do not use normalise_hist!
                        syst_unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(syst_unfolded_hist_bin_total_errors)
                        syst_entries_div_bin_width.append(
                            Contribution(syst_unfolded_hist_bin_total_errors_div_bin_width,
                                         label="Unfolded (#tau = %.3g) (total err)\n(%s)" % (syst_tau, syst_label),
                                         line_color=syst_dict['colour'], line_width=lw, line_style=1,
                                         marker_color=syst_dict['colour'], marker_size=0,
                                         subplot=mc_gen_hist_bin_div_bin_width,
                                         normalise_hist=False),
                        )

                    entries = [
                        Contribution(mc_gen_hist_bin,
                                     label="Generator (%s)" % (region['mc_label']),
                                     line_color=gen_colour, line_width=lw,
                                     marker_color=gen_colour, marker_size=0,
                                     normalise_hist=True),
                        Contribution(unfolded_hist_bin_total_errors,
                                     label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                     line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                     marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                     normalise_hist=True),
                        Contribution(unfolded_hist_bin_stat_errors,
                                     label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                     line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                     marker_color=unfolded_stat_colour, marker_size=0,
                                     normalise_hist=True),
                    ]
                    entries.extend(syst_entries)
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=particle_title,
                                ytitle=normalised_differential_label,
                                subplot_title='#splitline{Unfolded / Gen}{(%s)}' % (region['mc_label']),
                                **common_hist_args)
                    plot.legend.SetX1(0.55)
                    plot.legend.SetY1(0.72)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.88)
                    plot.legend.SetNColumns(3)
                    plot.plot("NOSTACK E1")
                    plot.save("%s/unfolded_%s_syst_model_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                    # Do a version where divided by bin width
                    mc_gen_hist_bin_div_bin_width = qgp.hist_divide_bin_width(mc_gen_hist_bin)
                    unfolded_hist_bin_stat_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_stat_errors)
                    unfolded_hist_bin_total_errors_div_bin_width = qgp.hist_divide_bin_width(unfolded_hist_bin_total_errors)
                    entries_div_bin_width = [
                        Contribution(mc_gen_hist_bin_div_bin_width,
                                     label="Generator (%s)" % (region['mc_label']),
                                     line_color=gen_colour, line_width=lw,
                                     marker_color=gen_colour, marker_size=0,
                                     normalise_hist=False),
                        Contribution(unfolded_hist_bin_total_errors_div_bin_width,
                                     label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                     line_color=unfolded_total_colour, line_width=lw, line_style=1,
                                     marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                     normalise_hist=False),
                        Contribution(unfolded_hist_bin_stat_errors_div_bin_width,
                                     label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                     line_color=unfolded_stat_colour, line_width=lw, line_style=1,
                                     marker_color=unfolded_stat_colour, marker_size=0,
                                     normalise_hist=False),
                    ]
                    entries_div_bin_width.extend(syst_entries_div_bin_width)
                    if not check_entries(entries_div_bin_width, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries_div_bin_width,
                                xtitle=particle_title,
                                ytitle=normalised_differential_label,
                                subplot_title='#splitline{Unfolded / Gen}{(%s)}' % (region['mc_label']),
                                **common_hist_args)
                    plot.legend.SetX1(0.55)
                    plot.legend.SetY1(0.72)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.88)
                    plot.legend.SetNColumns(3)
                    plot.plot("NOSTACK E1")
                    plot.save("%s/unfolded_%s_syst_model_bin_%d_divBinWidth.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))


                # --------------------------------------------------------------
                # PLOT UNCERTAINTY SHIFTS
                # --------------------------------------------------------------
                systematic_shift_hists_bin = [unfolder.get_var_hist_pt_binned(h, ibin_pt, binning_scheme='generator')
                                              for h in systematic_shift_hists]
                unfolded_stat_error_bin = unfolder.get_var_hist_pt_binned(error_stat_1d, ibin_pt, binning_scheme="generator")
                unfolded_total_error_bin =  unfolder.get_var_hist_pt_binned(unfolded_1d, ibin_pt, binning_scheme="generator")
                plot_uncertainty_shifts(total_hist=unfolded_total_error_bin,
                                        stat_hist=unfolded_stat_error_bin,
                                        syst_shifts=systematic_shift_hists_bin,
                                        systs=region['systematics'],
                                        output_filename='%s/unfolded_systs_%s_bin_%d.%s' % (this_output_dir, append, ibin_pt, OUTPUT_FMT),
                                        title=title,
                                        angle_str=angle_str)

                # --------------------------------------------------------------
                # PLOT RECO-LEVEL DISTRIBUTIONS (with gen binning)
                # --------------------------------------------------------------
                # Reco level, generator-binning
                reco_mc_colour = ROOT.kGreen+2
                reco_data_colour = ROOT.kRed
                entries = [
                    Contribution(mc_reco_hist_bin_gen_binning,
                                 label="MC",
                                 line_color=reco_mc_colour, line_width=lw,
                                 marker_color=reco_mc_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(reco_hist_bin_gen_binning,
                                 label="Data",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_reco_hist_bin_gen_binning,
                                 normalise_hist=True),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Data / MC',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/detector_gen_binning_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                if SUBTRACT_FAKES:
                    # Same but background-subtracted
                    entries = [
                        Contribution(mc_reco_hist_bg_subtracted_bin_gen_binning,
                                     label="MC (bg-subtracted)",
                                     line_color=reco_mc_colour, line_width=lw,
                                     marker_color=reco_mc_colour, marker_size=0,
                                     normalise_hist=True),
                        Contribution(reco_hist_bg_subtracted_bin_gen_binning,
                                     label="Data (bg-subtracted)",
                                     line_color=reco_data_colour, line_width=lw,
                                     marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                     subplot=mc_reco_hist_bg_subtracted_bin_gen_binning,
                                     normalise_hist=True),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=detector_title,
                                ytitle=normalised_differential_label,
                                subplot_title='Data / MC',
                                **common_hist_args)
                    plot.legend.SetX1(0.6)
                    plot.legend.SetY1(0.75)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.9)
                    plot.plot("NOSTACK E1")
                    plot.save("%s/detector_gen_binning_bg_subtracted_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))


            # Draw individual pt bin plots - RECO binning
            # ------------------------------------------------------------------
            for ibin_pt in range(0, len(pt_bin_edges_reco)-1):
                this_pt_bin_tdir = this_tdir.mkdir("reco_bin_%d" % (ibin_pt))
                print("Individual detector bin", ibin_pt, "=", pt_bin_edges_reco[ibin_pt], pt_bin_edges_reco[ibin_pt+1])

                # Get 1D hists
                # The folded unfolded result
                folded_unfolded_hist_bin_reco_binning = unfolder.get_var_hist_pt_binned(hist_data_folded, ibin_pt, binning_scheme="detector")
                this_pt_bin_tdir.WriteTObject(folded_unfolded_hist_bin_reco_binning, "folded_unfolded_hist_bin_reco_binning")

                # here this is the thing to be unfolded, should be data or MC depending on MC_INPUT flag
                reco_hist_bin_reco_binning = unfolder.get_var_hist_pt_binned(reco_1d, ibin_pt, binning_scheme="detector")
                this_pt_bin_tdir.WriteTObject(reco_hist_bin_reco_binning, "reco_hist_bin_reco_binning")

                if SUBTRACT_FAKES:
                    reco_hist_bg_subtracted_bin_reco_binning = unfolder.get_var_hist_pt_binned(reco_1d_bg_subtracted, ibin_pt, binning_scheme="detector")
                    this_pt_bin_tdir.WriteTObject(reco_hist_bg_subtracted_bin_reco_binning, "reco_hist_bg_subtracted_bin_reco_binning")

                # Get MC hists
                mc_reco_hist_bin_reco_binning = unfolder.get_var_hist_pt_binned(hist_mc_reco, ibin_pt, binning_scheme="detector")
                this_pt_bin_tdir.WriteTObject(mc_reco_hist_bin_reco_binning, "mc_reco_hist_bin_reco_binning")

                if SUBTRACT_FAKES:
                    mc_reco_hist_bg_subtracted_bin_reco_binning = unfolder.get_var_hist_pt_binned(hist_mc_reco_bg_subtracted, ibin_pt, binning_scheme="detector")
                    this_pt_bin_tdir.WriteTObject(mc_reco_hist_bg_subtracted_bin_reco_binning, "mc_reco_hist_bg_subtracted_bin_reco_binning")

                # Do the plots
                # --------------------------------------------------------------
                # common hist settings
                lw = 1
                title = "%s\n%s region\n%g < p_{T}^{Jet} < %g GeV" % (jet_algo, region['label'], pt_bin_edges_reco[ibin_pt], pt_bin_edges_reco[ibin_pt+1])
                common_hist_args = dict(
                    what="hist",
                    title=title,
                    subplot_type='ratio',
                    subplot_limits=(0.5, 1.5),
                )

                # Reco only, detector-binning
                entries = [
                    Contribution(mc_reco_hist_bin_reco_binning,
                                 label="MC",
                                 line_color=reco_mc_colour, line_width=lw,
                                 marker_color=reco_mc_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(reco_hist_bin_reco_binning,
                                 label="Data",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_reco_hist_bin_reco_binning,
                                 normalise_hist=True),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Data / MC',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/detector_reco_binning_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                if SUBTRACT_FAKES:
                    # Same but background-subtracted
                    entries = [
                        Contribution(mc_reco_hist_bg_subtracted_bin_reco_binning,
                                     label="MC (bg-subtracted)",
                                     line_color=reco_mc_colour, line_width=lw,
                                     marker_color=reco_mc_colour, marker_size=0,
                                     normalise_hist=True),
                        Contribution(reco_hist_bg_subtracted_bin_reco_binning,
                                     label="Data (bg-subtracted)",
                                     line_color=reco_data_colour, line_width=lw,
                                     marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                     subplot=mc_reco_hist_bg_subtracted_bin_reco_binning,
                                     normalise_hist=True),
                    ]
                    if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                        continue
                    plot = Plot(entries,
                                xtitle=detector_title,
                                ytitle="p.d.f.",
                                subplot_title='Data / MC',
                                **common_hist_args)
                    plot.legend.SetX1(0.6)
                    plot.legend.SetY1(0.75)
                    plot.legend.SetX2(0.98)
                    plot.legend.SetY2(0.9)
                    plot.plot("NOSTACK E1")
                    plot.save("%s/detector_reco_binning_bg_subtracted_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # PLOT FOLDED UNFOLDED DATA
                # --------------------------------------------------------------
                # Reco + folded, detector binning
                # FIXME: shouldn't this use RECO gen bin
                reco_mc_colour = ROOT.kGreen+2
                reco_data_colour = ROOT.kRed
                reco_folded_colour = ROOT.kAzure+1
                entries = [
                    Contribution(mc_reco_hist_bg_subtracted_bin_reco_binning if SUBTRACT_FAKES else mc_reco_hist_bin_reco_binning,
                                 label="MC (reco, bg-subtracted)" if SUBTRACT_FAKES else "MC (reco)",
                                 line_color=reco_mc_colour, line_width=lw,
                                 marker_color=reco_mc_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(reco_hist_bg_subtracted_bin_reco_binning if SUBTRACT_FAKES else reco_hist_bin_reco_binning,
                                 label="Data (reco, bg-subtracted)" if SUBTRACT_FAKES else "Data (reco)",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_size=0,
                                 subplot=mc_reco_hist_bin_reco_binning,
                                 normalise_hist=True),
                    Contribution(folded_unfolded_hist_bin_reco_binning,
                                 label="Folded unfolded data (#tau = %.3g)" % (tau),
                                 line_color=reco_folded_colour, line_width=lw,
                                 marker_color=reco_folded_colour, marker_size=0,
                                 subplot=mc_reco_hist_bin_reco_binning,
                                 normalise_hist=True),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Data / MC',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.72)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/detector_folded_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Folded, but only comparing data with data to check it is sane
                entries = [
                    Contribution(reco_hist_bg_subtracted_bin_reco_binning if SUBTRACT_FAKES else reco_hist_bin_reco_binning,
                                 label="Data (reco, bg-subtracted)" if SUBTRACT_FAKES else "Data (reco)",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(folded_unfolded_hist_bin_reco_binning,
                                 label="Folded unfolded data (#tau = %.3g)" % (tau),
                                 line_color=reco_folded_colour, line_width=lw,
                                 marker_color=reco_folded_colour, marker_size=0,
                                 subplot=reco_hist_bin_reco_binning,
                                 normalise_hist=True),
                ]
                if not check_entries(entries, "%s %d" % (append, ibin_pt)):
                    continue
                plot = Plot(entries,
                            xtitle=detector_title,
                            ytitle=normalised_differential_label,
                            subplot_title='Folded unfolded / reco',
                            **common_hist_args)
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/detector_folded_only_data_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))


            # DO SUMMARY PLOT
            # ------------------------------------------------------------------
            if args.doSummaryPlot:
                marker = ""
                if "_" in angle.name or "^" in angle.name:
                    marker = "$"
                var_label = "Particle-level " + marker + angle.name + marker + " ($%s$)" % angle.lambda_str
                v = "%s_vs_pt" % (angle.var)
                bins = [(pt_bin_edges_gen[i], pt_bin_edges_gen[i+1]) for i in range(len(pt_bin_edges_gen)-1)]
                print("Making summary plot from pt bins:", bins)
                xlim = (50, 2000)
                if "ZPlusJets" in region['name']:
                    xlim = (50, 614)
                qgp.do_mean_rms_summary_plot(summary_1d_entries, bins,
                                             "%s/%s_box_dijet_mpl.%s" % (this_output_dir, v, OUTPUT_FMT),
                                             var_label=var_label,
                                             xlim=xlim,
                                             region_title=region['label'].lower())

    print("Saved hists to", output_tfile.GetName())
