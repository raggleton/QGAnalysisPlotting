#!/usr/bin/env python


"""TUnfold it all

Thanks to Ashley, Dennis
"""


from __future__ import print_function

import os
import sys
import argparse
from array import array
import numpy as np
import math

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()

import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)


# Control plot output format
OUTPUT_FMT = "pdf"


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

        self.variable_bin_edges_reco = variable_bin_edges_reco
        self.nbins_variable_reco = len(variable_bin_edges_reco)-1
        self.variable_bin_edges_gen = variable_bin_edges_gen
        self.nbins_variable_gen = len(variable_bin_edges_gen)-1

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
        self.detector_distribution_underflow.AddAxis(self.variable_name, self.nbins_variable_reco, self.variable_bin_edges_reco, var_uf, var_of)
        self.detector_distribution_underflow.AddAxis("pt", self.nbins_pt_underflow_reco, self.pt_bin_edges_underflow_reco, pt_uf, pt_of)
        
        self.detector_distribution = self.detector_binning.AddBinning("detectordistribution")
        self.detector_distribution.AddAxis(self.variable_name, self.nbins_variable_reco, self.variable_bin_edges_reco, var_uf, var_of)
        self.detector_distribution.AddAxis("pt", self.nbins_pt_reco, self.pt_bin_edges_reco, pt_uf, pt_of)


        self.generator_binning = ROOT.TUnfoldBinning("generator")

        self.generator_distribution_underflow = self.generator_binning.AddBinning("generatordistribution_underflow")
        self.generator_distribution_underflow.AddAxis(self.variable_name, self.nbins_variable_gen, self.variable_bin_edges_gen, var_uf, var_of)
        self.generator_distribution_underflow.AddAxis("pt", self.nbins_pt_underflow_gen, self.pt_bin_edges_underflow_gen, pt_uf, pt_of)
        self.generator_distribution = self.generator_binning.AddBinning("generatordistribution")
        self.generator_distribution.AddAxis(self.variable_name, self.nbins_variable_gen, self.variable_bin_edges_gen, var_uf, var_of)
        self.generator_distribution.AddAxis("pt", self.nbins_pt_gen, self.pt_bin_edges_gen, pt_uf, pt_of)

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

    def doScanTau(self, output_dir, n_scan=300, tau_min=1E-14, tau_max=1E-4, scan_mode=ROOT.TUnfoldDensity.kEScanTauRhoMax):
        """Figure out best tau by scanning tau curve
        Taken from Ashley's code
        """

        # Graphs to save output from scan over tau
        scanresults = ROOT.MakeNullPointer(ROOT.TSpline)
        lCurve = ROOT.MakeNullPointer(ROOT.TGraph)
        LogTauX = ROOT.MakeNullPointer(ROOT.TSpline)
        LogTauY = ROOT.MakeNullPointer(ROOT.TSpline)

        iBestAvg = self.tunfolder.ScanTau(n_scan,
                                         tau_min,
                                         tau_max,
                                         scanresults,
                                         scan_mode,
                                         "generator",
                                         self.axisSteering,
                                         lCurve,
                                         LogTauX,
                                         LogTauY)
        tau = self.tunfolder.GetTau()

        x, y, t, rho = array('d'), array('d'), array('d'), array('d')
        t0 = ROOT.Double(0.0)
        rho0 = ROOT.Double(0.0)
        scanresults.GetKnot(iBestAvg, t0, rho0)
        t.append(t0)
        rho.append(rho0)

        x0 = ROOT.Double(0.0)
        y0 = ROOT.Double(0.0)
        lCurve.GetPoint(iBestAvg, x0, y0)
        x.append(x0)
        y.append(y0)
        print("t[0] = {}, type = {}".format(t[0], type(t)))
        print("rho[0] = {}".format(rho[0]))
        print("x[0] = {}".format(x[0]))
        print("y[0] = {}".format(y[0]))
        bestRhoLogTau = ROOT.TGraph(1, t, rho)
        print("10^log_10(tau) = tau = {}".format(math.pow(10., float(t0))))
        bestLCurve = ROOT.TGraph(1, x, y)

        tAll, rhoAll = array('d'), array('d')
        for i in range(n_scan):
            tt = ROOT.Double(0.0)
            rr = ROOT.Double(0.0)
            scanresults.GetKnot(i, tt, rr)
            tAll.append(tt)
            rhoAll.append(rr)

        knots = ROOT.TGraph(int(n_scan), tAll, rhoAll)
        print("chi**2 A {:3.1f} + chi**2 L {:3.1f} / NDOF {:3.1f} ".format(self.tunfolder.GetChi2A(),
                                                                           self.tunfolder.GetChi2L(),
                                                                           self.tunfolder.GetNdf()))

        print("doScanTau value is {}".format(tau))
        print("Returned the TUnfoldDensity object.")

        canv_tauScan = ROOT.TCanvas("canv_tauScan"+str(tau), "canv_tauScan"+str(tau))

        knots.SetLineColor(ROOT.kBlue+3)

        bestRhoLogTau.SetMarkerColor(ROOT.kRed)

        knots.Draw()
        bestRhoLogTau.Draw("* same")
        knots.GetXaxis().SetTitle("log_{10}(#tau)")
        knots.GetYaxis().SetTitle(" #rho")
        tauoptname = 'NONE'
        if scan_mode == ROOT.TUnfoldDensity.kEScanTauRhoMax:
            tauoptname = 'Max'
        else:
            tauoptname = 'Average'
        knots.SetTitle("Optimization of Regularization Parameter, #tau : Scan of {} #rho".format(tauoptname))

        leg1 = ROOT.TLegend(0.59, 0.6, 0.74, 0.89)
        leg1.SetFillColor(0)
        leg1.SetFillStyle(0)
        leg1.SetBorderSize(0)
        leg1.SetTextSize(0.026)
        leg1.AddEntry(knots, '{} #rho: #tau = {}'.format(tauoptname,  tau), 'l')
        leg1.Draw()

        canv_tauScan.Modified()
        canv_tauScan.Update()
        canv_tauScan.Print("%s/scantau_%s.%s" % (output_dir, self.variable_name, self.OUTPUT_FMT))
        return tau


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

        print("var_bins:", var_bins)
        print("pt_bins:", pt_bins)
        bin_num = binning.GetGlobalBinNumber(0.001, 49)
        print("Global bin num for (pt, lambda) = (49, 0.001) => %d" % (bin_num))

        # need the -1 on ibin_pt, as it references an array index, whereas ROOT bins start at 1
        h = ROOT.TH1D("h_%d_%s" % (ibin_pt, cu.get_unique_str()), "", len(var_bins)-1, var_bins)
        for var_ind, var_value in enumerate(var_bins[:-1], 1):
            this_val = var_value * 1.001  # ensure its inside
            bin_num = binning.GetGlobalBinNumber(this_val, pt_bins[ibin_pt]*1.001)
            print("Global bin num for (pt, lambda) = (%.3f, %.3f) => %d" % (pt_bins[ibin_pt]*1.001, this_val, bin_num))
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
    corr_map.Draw("COLZ0")
    canv.SaveAs(output_filename)
    ROOT.gStyle.SetPalette(ROOT.kBird)


def plot_simple_unfolded(unfolded, tau, reco, gen, output_filename, title=""):
    """Simple plot of unfolded, reco, gen, by bin number (ie non physical axes)"""
    canv_unfold = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv_unfold.SetLogy()
    canv_unfold.SetTicks(1, 1)
    leg = ROOT.TLegend(0.7, 0.7, 0.88, 0.88)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    hst = ROOT.THStack("hst", "%s;Bin Number;N" % (title))

    reco.SetLineColor(ROOT.kGreen+2)
    reco.SetMarkerColor(ROOT.kGreen+2)
    hst.Add(reco)
    leg.AddEntry(reco, "Reco", "L")

    gen.SetLineColor(ROOT.kBlue)
    hst.Add(gen)
    leg.AddEntry(gen, "Gen", "L")

    unfolded.SetLineColor(ROOT.kRed)
    unfolded.SetLineWidth(0)
    unfolded.SetMarkerColor(ROOT.kRed)
    unfolded.SetMarkerSize(0.6)
    unfolded.SetMarkerStyle(20)
    hst.Add(unfolded)
    leg.AddEntry(unfolded, "Unfolded (#tau = %.3g)" % (tau), "LP")

    hst.Draw("NOSTACK HISTE")
    leg.Draw()
    hst.SetMinimum(1E-2)
    hst.SetMaximum(5*max([h.GetMaximum() for h in [unfolded, reco, gen]]))
    canv_unfold.SaveAs(output_filename)


def plot_simple_detector(reco_data, reco_mc, output_filename, title):
    """Plot detector-level quantities for data & MC, by bin number (ie non physical axes)"""
    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetLogy()
    canv.SetTicks(1, 1)
    leg = ROOT.TLegend(0.7, 0.7, 0.88, 0.88)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    hst = ROOT.THStack("hst", "%s;Bin Number;N" % (title))

    reco_mc.SetLineColor(ROOT.kGreen+2)
    reco_mc.SetMarkerColor(ROOT.kGreen+2)
    hst.Add(reco_mc)
    leg.AddEntry(reco_mc, "MC [detector-level]", "L")

    reco_data.SetLineColor(ROOT.kRed)
    reco_data.SetLineWidth(0)
    reco_data.SetMarkerColor(ROOT.kRed)
    reco_data.SetMarkerSize(0.6)
    reco_data.SetMarkerStyle(20)
    hst.Add(reco_data)
    leg.AddEntry(reco_data, "Data [detector-level]", "LP")

    hst.Draw("NOSTACK HISTE")
    leg.Draw()
    hst.SetMinimum(1E-2)
    hst.SetMaximum(5*max([h.GetMaximum() for h in [reco_data, reco_mc]]))
    canv.SaveAs(output_filename)


def create_hist_with_errors(hist, err_matrix):
    hnew = hist.Clone(cu.get_unique_str())
    nbins = hist.GetNbinsX()
    for i in range(1, nbins+1):
        err = math.sqrt(err_matrix.GetBinContent(i, i))
        hnew.SetBinError(i, err)
    return hnew


def make_hist_from_diagonals(h2d, do_sqrt=True):
    nbins = h2d.GetNbinsX()
    hnew = ROOT.TH1D("h_diag" + cu.get_unique_str(), "", nbins, 0, nbins)
    for i in range(1, nbins+1):
        err = h2d.GetBinContent(i, i)
        if do_sqrt:
            err = math.sqrt(err)
        hnew.SetBinError(i, err)
    return hnew


def update_hist_bin_content(h_orig, h_to_be_updated):
    if h_orig.GetNbinsX() != h_to_be_updated.GetNbinsX():
        raise RuntimeError("Need same # x bins")
    for i in range(1, h_orig.GetNbinsX()+1):
        h_to_be_updated.SetBinContent(i, h_orig.GetBinContent(i))


def draw_projection_comparison(h_orig, h_projection, title, xtitle, output_filename):
    """Draw 2 hists, h_orig the original, and h_projection the projection of a 2D hist"""

    # Check integrals
    int_orig = h_orig.Integral()
    int_proj = h_projection.Integral()
    if abs(int_orig - int_proj)/int_orig > 0.01:
        print("draw_projection_comparison: different integrals: %f vs %f" % (int_orig, int_proj))

    # Check bin-by-bin
    for i in range(1, h_orig.GetNbinsX()+1):
        value_orig = h_orig.GetBinContent(i)
        value_proj = h_projection.GetBinContent(i)
        if abs(value_orig - value_proj) > 1E-2:
            print("draw_projection_comparison: bin %s has different contents: %f vs %f" % (i, value_orig, value_proj))

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
        Contribution(hist_reco_mc, label="MC (detector)",
                     line_color=ROOT.kBlack, line_width=1,
                     marker_color=ROOT.kBlack, marker_size=0,
                     normalise_hist=False),
        Contribution(hist_reco_data, label="Data (detector)",
                     line_color=ROOT.kBlue, line_width=1,
                     marker_color=ROOT.kBlue, marker_size=0,
                     normalise_hist=False),
        Contribution(hist_folded, label="Data (#tau = %.3g, folded)" % (tau),
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


if __name__ == "__main__":
    # FOR Z+JETS:
    input_mc_dy_mgpythia_tfile = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_noZReweight/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root")
    
    # input_mc_dy_mgherwig_tfile = cu.open_root_file("workdir_ak4puppi_herwig_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1RecoConstituents_V11JEC_JER_tUnfold/uhh2.AnalysisModuleRunner.MC.MC_MG_HERWIG_DYJetsToLL.root")
    # input_mc_dy_herwig_tfile = cu.open_root_file("workdir_ak4puppi_herwig_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1RecoConstituents_V11JEC_JER_tUnfold/uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL.root")

    # input_mc_dy_mgpythia_neutralUp_tfile = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1RecoConstituents_V11JEC_JER_tUnfold_neutralUp/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root")
    # input_mc_dy_mgpythia_neutralDown_tfile = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1RecoConstituents_V11JEC_JER_tUnfold_neutralDown/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root")

    input_singlemu_tfile = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu.root")

    # FOR DIJET:
    input_mc_qcd_mgpythia_tfile = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_noZReweight/uhh2.AnalysisModuleRunner.MC.MC_QCD.root")
    
    # input_mc_qcd_pythia_tfile = cu.open_root_file("workdir_ak4puppi_pythia_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1RecoConstituents_V11JEC_JER_tUnfold/uhh2.AnalysisModuleRunner.MC.MC_PYTHIA-QCD.root")
    # input_mc_qcd_herwig_tfile = cu.open_root_file("workdir_ak4puppi_herwig_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1RecoConstituents_V11JEC_JER_tUnfold/uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD.root")

    input_jetht_tfile = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5/uhh2.AnalysisModuleRunner.DATA.Data_JetHTZeroBias.root")

    # input_mc_qcd_mgpythia_neutralUp_tfile = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1RecoConstituents_V11JEC_JER_tUnfold_neutralUp/uhh2.AnalysisModuleRunner.MC.MC_QCD.root")
    # input_mc_qcd_mgpythia_neutralDown_tfile = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1RecoConstituents_V11JEC_JER_tUnfold_neutralDown/uhh2.AnalysisModuleRunner.MC.MC_QCD.root")

    regions = [
        {
            "name": "Dijet",
            "dirname": "Dijet_QG_unfold_tighter",
            "label": "Dijet",
            "data_tfile": input_jetht_tfile,
            "mc_tfile": input_mc_qcd_mgpythia_tfile,
            # "mc_neutralUp_tfile": input_mc_qcd_mgpythia_neutralUp_tfile,
            # "mc_neutralDown_tfile": input_mc_qcd_mgpythia_neutralDown_tfile,
            "tau_limits": {
                'jet_puppiMultiplicity': (1E-13, 1E-8),
                'jet_pTD': (1E-13, 1E-8),
                'jet_LHA': (1E-13, 1E-8),
                'jet_width': (1E-13, 1E-8),
                'jet_thrust': (1E-13, 1E-8),
                'jet_puppiMultiplicity_charged': (1E-16, 1E-8),
                'jet_pTD_charged': (1E-13, 1E-8),
                'jet_LHA_charged': (1E-13, 1E-8),
                'jet_width_charged': (1E-13, 1E-8),
                'jet_thrust_charged': (1E-13, 1E-8),
            },
        },
        {
            "name": "ZPlusJets",
            "dirname": "ZPlusJets_QG_unfold",
            "label": "Z+jets",
            "data_tfile": input_singlemu_tfile,
            "mc_tfile": input_mc_dy_mgpythia_tfile,
            # "mc_neutralUp_tfile": input_mc_dy_mgpythia_neutralUp_tfile,
            # "mc_neutralDown_tfile": input_mc_dy_mgpythia_neutralDown_tfile,
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
        },
    ]

    regularise = None
    # regularise = "tau"
    # regularise = "L"

    # Run with MC input instead of data
    MC_input = False
    mc_append = "_MC" if MC_input else ""
    
    output_dir = "unfolding_regularise%s_target0p5_uflowFirst%s" % (regularise, mc_append)
    cu.check_dir_exists_create(output_dir)

    # TODO automate this
    jet_algo = "AK4 PF PUPPI"

    # Save hists etc to ROOT file for access later
    output_tfile = ROOT.TFile("%s/unfolding_result.root" % (output_dir), "RECREATE")

    for region in regions[:]:

        # Setup pt bins
        # -------------
        # need different ones for Z+Jets region
        is_zpj = "ZPlusJets" in region['name']

        zpj_append = "_zpj" if is_zpj else ""

        pt_bin_edges_gen = qgc.PT_UNFOLD_DICT['signal%s_gen' % (zpj_append)]
        pt_bin_edges_reco = qgc.PT_UNFOLD_DICT['signal%s_reco' % (zpj_append)]
        # pt_bin_edges_reco = qgc.construct_fine_binning(pt_bin_edges_gen)

        pt_bin_edges_underflow_gen = qgc.PT_UNFOLD_DICT['underflow%s_gen' % (zpj_append)]
        pt_bin_edges_underflow_reco = qgc.PT_UNFOLD_DICT['underflow%s_reco' % (zpj_append)]
        # pt_bin_edges_underflow_reco = qgc.construct_fine_binning(pt_bin_edges_underflow_gen)

        for angle in qgc.COMMON_VARS[:]:
            append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc

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

            hist_data_reco = cu.get_from_tfile(region['data_tfile'], "%s/hist_%s_reco_new" % (region['dirname'], angle_shortname))
            hist_mc_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_new" % (region['dirname'], angle_shortname))
            hist_mc_gen = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_truth_new" % (region['dirname'], angle_shortname))
            hist_mc_gen_reco_map = cu.get_from_tfile(region['mc_tfile'], "%s/tu_%s_GenReco_new" % (region['dirname'], angle_shortname))

            # hist_mc_gen_reco_neutralUp_map = cu.get_from_tfile(region['mc_neutralUp_tfile'], "%s/tu_%s_GenReco_new" % (region['dirname'], angle_shortname))
            # hist_mc_gen_reco_neutralDown_map = cu.get_from_tfile(region['mc_neutralDown_tfile'], "%s/tu_%s_GenReco_new" % (region['dirname'], angle_shortname))

            hist_data_reco_gen_binning = cu.get_from_tfile(region['data_tfile'], "%s/hist_%s_reco_gen_binning_new" % (region['dirname'], angle_shortname))
            hist_mc_reco_gen_binning = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_gen_binning_new" % (region['dirname'], angle_shortname))

            # Setup unfolder object
            # ---------------------
            unfolder = MyUnfolder(response_map=hist_mc_gen_reco_map,
                                  variable_bin_edges_reco=angle_bin_edges_reco,
                                  variable_bin_edges_gen=angle_bin_edges_gen,
                                  variable_name=angle.name,
                                  pt_bin_edges_reco=pt_bin_edges_reco,
                                  pt_bin_edges_gen=pt_bin_edges_gen,
                                  pt_bin_edges_underflow_reco=pt_bin_edges_underflow_reco,
                                  pt_bin_edges_underflow_gen=pt_bin_edges_underflow_gen,
                                  orientation=ROOT.TUnfold.kHistMapOutputHoriz,
                                  # constraintMode=ROOT.TUnfold.kEConstraintArea,
                                  constraintMode=ROOT.TUnfold.kEConstraintNone,
                                  regMode=ROOT.TUnfold.kRegModeCurvature,
                                  # densityFlags=ROOT.TUnfoldDensity.kDensityModeBinWidth,
                                  densityFlags=ROOT.TUnfoldDensity.kDensityModeNone,
                                  axisSteering='*[b]')

            unfolder.save_binning(txt_filename="%s/binning_scheme.txt" % (this_output_dir), print_xml=False)

            # unfolder.tunfolder.AddSysError(hist_mc_gen_reco_neutralUp_map, "NeutralUp", ROOT.TUnfold.kHistMapOutputHoriz, ROOT.TUnfoldDensity.kSysErrModeMatrix)
            # unfolder.tunfolder.AddSysError(hist_mc_gen_reco_neutralDown_map, "NeutralDown", ROOT.TUnfold.kHistMapOutputHoriz, ROOT.TUnfoldDensity.kSysErrModeMatrix)

            # Set what is to be unfolded
            # ---------------------
            reco_1d = hist_mc_reco if MC_input else hist_data_reco
            unfolder.setInput(reco_1d)
            this_tdir.WriteTObject(reco_1d, "unfold_input")

            reco_1d_gen_binning = hist_data_reco_gen_binning
            # reco_1d_gen_binning = hist_mc_reco_gen_binning # only for testing that everything is setup OK
            this_tdir.WriteTObject(reco_1d_gen_binning, "reco_1d_gen_binning")

            # Do any regularisation
            # ---------------------
            # tau = 1E-10
            tau = 0
            if regularise == "L":
                tau = unfolder.doScanL(output_dir=this_output_dir, n_scan=100,
                                       tau_min=1E-14, tau_max=1E-4)
            elif regularise == "tau":
                tau = unfolder.doScanTau(output_dir=this_output_dir, n_scan=100,
                                         tau_min=region['tau_limits'][angle.var][0],
                                         tau_max=region['tau_limits'][angle.var][1],
                                         scan_mode=ROOT.TUnfoldDensity.kEScanTauRhoAvgSys)

            # Do unfolding!
            # ---------------------
            unfolder.do_unfolding(tau)
            unfolded_1d = unfolder.get_output()
            unfolded_1d.SetName("unfolded_1d")
            this_tdir.WriteTObject(unfolded_1d)

            # stat errors only - do before or after systematics?
            ematrix_input = unfolder.get_ematrix_input() # stat errors from input to be unfolded
            this_tdir.WriteTObject(ematrix_input, "ematrix_input")
            ematrix_sys_uncorr = unfolder.get_ematrix_sys_uncorr() # stat errors in response matrix
            this_tdir.WriteTObject(ematrix_sys_uncorr, "ematrix_sys_uncorr")

            ematrix_stat_sum = ematrix_input.Clone("ematrix_stat_sum")
            ematrix_stat_sum.Add(ematrix_sys_uncorr)

            error_stat_1d = make_hist_from_diagonals(ematrix_stat_sum, do_sqrt=True)

            ematrix_total = unfolder.get_ematrix_total()
            error_total_1d = make_hist_from_diagonals(ematrix_total, do_sqrt=True)
            this_tdir.WriteTObject(ematrix_total, "ematrix_total_1d")

            angle_str = "%s (%s)" % (angle.name, angle.lambda_str)

            # Draw unified unfolded distributions
            # ---------------------
            plot_simple_unfolded(unfolded=unfolded_1d,
                                 tau=tau,
                                 reco=reco_1d,
                                 gen=hist_mc_gen,
                                 output_filename="%s/unfolded_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                                 title="%s region, %s" % (region['label'], angle_str))

            plot_simple_detector(reco_data=reco_1d_gen_binning,
                                 reco_mc=hist_mc_reco_gen_binning,
                                 output_filename="%s/detector_gen_binning_%s.%s" % (this_output_dir, append, OUTPUT_FMT),
                                 title="%s region, %s" % (region['label'], angle_str))

            # Draw collapsed distributions
            # ---------------------
            # unfolded_1d_pt = unfolder.get_output_pt_1d()
            # hist_mc_gen

            # unfolded_1d_lambda = unfolder.get_output_lambda_1d()

            # Draw projections of response matrix vs 1D hist to check normalisation OK
            # ---------------------
            proj_reco = hist_mc_gen_reco_map.ProjectionY("proj_reco_%s" % (append))
            draw_projection_comparison(hist_mc_reco, proj_reco,
                                       title="%s\n%s region" % (jet_algo, region['label']),
                                       xtitle="%s, Detector binning" % (angle_str),
                                       output_filename="%s/projection_reco_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

            proj_gen = hist_mc_gen_reco_map.ProjectionX("proj_gen_%s" % (append))
            draw_projection_comparison(hist_mc_gen, proj_gen,
                                       title="%s\n%s region" % (jet_algo, region['label']),
                                       xtitle="%s, Generator binning" % (angle_str),
                                       output_filename="%s/projection_gen_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

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
            # ----------------------------
            hist_data_folded = unfolder.get_folded_output(hist_name="folded_%s" % (append))
            this_tdir.WriteTObject(hist_data_folded, "folded_1d")
            draw_reco_folded(hist_folded=hist_data_folded,
                             tau=tau,
                             hist_reco_data=reco_1d,
                             hist_reco_mc=hist_mc_reco,
                             title="%s\n%s region" % (jet_algo, region['label']),
                             xtitle="%s, Detector binning" % (angle_str),
                             output_filename="%s/folded_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

            # Draw individual pt bin plots
            # ----------------------------
            for ibin_pt in range(0, len(pt_bin_edges_gen)-1):

                this_pt_bin_tdir = this_tdir.mkdir("bin_%d" % (ibin_pt))

                # Produce 1D hist for this pt bin
                gen_hist_bin = unfolder.get_var_hist_pt_binned(hist_mc_gen, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(gen_hist_bin, "gen_hist_bin")

                unfolded_hist_bin = unfolder.get_var_hist_pt_binned(unfolded_1d, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(unfolded_hist_bin, "unfolded_hist_bin")

                unfolded_hist_bin_stat_errors = unfolder.get_var_hist_pt_binned(error_stat_1d, ibin_pt, binning_scheme="generator")
                update_hist_bin_content(unfolded_hist_bin, unfolded_hist_bin_stat_errors)
                this_pt_bin_tdir.WriteTObject(unfolded_hist_bin_stat_errors, "unfolded_hist_bin_stat_errors")

                unfolded_hist_bin_total_errors = unfolder.get_var_hist_pt_binned(error_total_1d, ibin_pt, binning_scheme="generator")
                update_hist_bin_content(unfolded_hist_bin, unfolded_hist_bin_total_errors)
                this_pt_bin_tdir.WriteTObject(unfolded_hist_bin_total_errors, "unfolded_hist_bin_total_errors")

                data_reco_hist_bin_gen_binning = unfolder.get_var_hist_pt_binned(hist_data_reco_gen_binning, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(data_reco_hist_bin_gen_binning, "data_reco_hist_bin_gen_binning")

                mc_reco_hist_bin_gen_binning = unfolder.get_var_hist_pt_binned(hist_mc_reco_gen_binning, ibin_pt, binning_scheme="generator")
                this_pt_bin_tdir.WriteTObject(mc_reco_hist_bin_gen_binning, "mc_reco_hist_bin_gen_binning")

                data_folded_hist_bin_reco_binning = unfolder.get_var_hist_pt_binned(hist_data_folded, ibin_pt, binning_scheme="detector")
                this_pt_bin_tdir.WriteTObject(data_folded_hist_bin_reco_binning, "data_folded_hist_bin_reco_binning")

                data_reco_hist_bin_reco_binning = unfolder.get_var_hist_pt_binned(reco_1d, ibin_pt, binning_scheme="detector")
                this_pt_bin_tdir.WriteTObject(data_reco_hist_bin_reco_binning, "data_reco_hist_bin_reco_binning")

                mc_reco_hist_bin_reco_binning = unfolder.get_var_hist_pt_binned(hist_mc_reco, ibin_pt, binning_scheme="detector")
                this_pt_bin_tdir.WriteTObject(mc_reco_hist_bin_reco_binning, "mc_reco_hist_bin_reco_binning")

                # print hist bins for check
                for n in range(1, gen_hist_bin.GetNbinsX()+1):
                    print("Bin", n)
                    print("gen_hist:", gen_hist_bin.GetBinContent(n), "+-", gen_hist_bin.GetBinError(n))
                    print("unfolded_hist:", unfolded_hist_bin.GetBinContent(n), "+-", unfolded_hist_bin.GetBinError(n))
                    print("unfolded_hist_bin_stat_errors:", unfolded_hist_bin_stat_errors.GetBinContent(n), "+-", unfolded_hist_bin_stat_errors.GetBinError(n))
                    print("unfolded_hist_bin_total_errors:", unfolded_hist_bin_total_errors.GetBinContent(n), "+-", unfolded_hist_bin_total_errors.GetBinError(n))

                # Make lots of plots
                # ------------------------------------------------------------
                # Unfolded only
                lw = 1
                gen_colour = ROOT.kBlue
                unfolded_basic_colour = ROOT.kAzure+10
                unfolded_stat_colour = ROOT.kGreen+1
                unfolded_total_colour = ROOT.kBlack

                # unnormalised version
                entries = [
                    Contribution(gen_hist_bin, label="Generator",
                                 line_color=gen_colour, line_width=lw,
                                 marker_color=gen_colour, marker_size=0,
                                 normalise_hist=False),
                    Contribution(unfolded_hist_bin, label="Unfolded (#tau = %.3g)" % (tau),
                                 line_color=unfolded_basic_colour, line_width=lw,
                                 marker_color=unfolded_basic_colour, marker_size=0,
                                 subplot=gen_hist_bin,
                                 normalise_hist=False),
                    Contribution(unfolded_hist_bin_stat_errors, label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                 line_color=unfolded_stat_colour, line_width=lw, line_style=2,
                                 marker_color=unfolded_stat_colour, marker_size=0,
                                 subplot=gen_hist_bin,
                                 normalise_hist=False),
                    Contribution(unfolded_hist_bin_total_errors, label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                 line_color=unfolded_total_colour, line_width=lw, line_style=3,
                                 marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                 subplot=gen_hist_bin,
                                 normalise_hist=False),
                ]
                has_entries = [c.obj.GetEntries() > 0 for c in entries]
                if not any(has_entries):
                    print("Skipping 0 entries in", append, ibin_pt)
                    continue
                title = "%s\n%s region\n%g < p_{T}^{Gen} < %g GeV" % (jet_algo, region['label'], pt_bin_edges_gen[ibin_pt], pt_bin_edges_gen[ibin_pt+1])
                plot = Plot(entries,
                            what="hist",
                            title=title,
                            xtitle="Generator-level "+angle_str,
                            ytitle="N",
                            subplot_type='ratio',
                            subplot_title='Unfolded / gen',
                            subplot_limits=(0.8, 1.2))
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.68)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/unfolded_unnormalised_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))
                
                # now normalise each plot to unity
                # Note that this modifies e.g. gen_hist_bin, so from this point
                # onwards it will be normalised to unity
                entries = [
                    Contribution(gen_hist_bin, label="Generator",
                                 line_color=gen_colour, line_width=lw,
                                 marker_color=gen_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(unfolded_hist_bin, label="Unfolded (#tau = %.3g)" % (tau),
                                 line_color=unfolded_basic_colour, line_width=lw,
                                 marker_color=unfolded_basic_colour, marker_size=0,
                                 subplot=gen_hist_bin,
                                 normalise_hist=True),
                    Contribution(unfolded_hist_bin_stat_errors, label="Unfolded (#tau = %.3g) (stat err)" % (tau),
                                 line_color=unfolded_stat_colour, line_width=lw, line_style=2,
                                 marker_color=unfolded_stat_colour, marker_size=0,
                                 subplot=gen_hist_bin,
                                 normalise_hist=True),
                    Contribution(unfolded_hist_bin_total_errors, label="Unfolded (#tau = %.3g) (total err)" % (tau),
                                 line_color=unfolded_total_colour, line_width=lw, line_style=3,
                                 marker_color=unfolded_total_colour, marker_style=20, marker_size=0.75,
                                 subplot=gen_hist_bin,
                                 normalise_hist=True),
                ]
                has_entries = [c.obj.GetEntries() > 0 for c in entries]
                if not any(has_entries):
                    print("Skipping 0 entries in", append, ibin_pt)
                    continue
                ytitle= "p.d.f."
                plot = Plot(entries,
                            what="hist",
                            title=title,
                            xtitle="Generator-level "+angle_str,
                            ytitle=ytitle,
                            subplot_type='ratio',
                            subplot_title='Unfolded / gen',
                            subplot_limits=(0.8, 1.2))
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.68)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/unfolded_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))



                # Reco only, generator-binning
                reco_mc_colour = ROOT.kGreen+2
                reco_data_colour = ROOT.kRed
                entries = [
                    Contribution(mc_reco_hist_bin_gen_binning, label="MC",
                                 line_color=reco_mc_colour, line_width=lw,
                                 marker_color=reco_mc_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(data_reco_hist_bin_gen_binning, label="Data",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                                 subplot=mc_reco_hist_bin_gen_binning,
                                 normalise_hist=True),
                ]
                has_entries = [c.obj.GetEntries() > 0 for c in entries]
                if not any(has_entries):
                    print("Skipping 0 entries in", append, ibin_pt)
                    continue
                title = "%s\n%s region\n%g < p_{T}^{Gen} < %g GeV" % (jet_algo, region['label'], pt_bin_edges_gen[ibin_pt], pt_bin_edges_gen[ibin_pt+1])
                plot = Plot(entries,
                            what="hist",
                            title=title,
                            xtitle="Detector-level " + angle_str,
                            ytitle=ytitle,
                            subplot_type='ratio',
                            subplot_title='Data / MC',
                            subplot_limits=(0.8, 1.2))
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/detector_gen_binning_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Reco + folded, detector binning
                reco_mc_colour = ROOT.kGreen+2
                reco_data_colour = ROOT.kRed
                reco_folded_colour = ROOT.kAzure+1
                entries = [
                    Contribution(mc_reco_hist_bin_reco_binning, label="MC",
                                 line_color=reco_mc_colour, line_width=lw,
                                 marker_color=reco_mc_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(data_reco_hist_bin_reco_binning, label="Data (reco)",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_size=0,
                                 subplot=mc_reco_hist_bin_reco_binning,
                                 normalise_hist=True),
                    Contribution(data_folded_hist_bin_reco_binning, label="Data (#tau = %.3g, folded)" % (tau),
                                 line_color=reco_folded_colour, line_width=lw,
                                 marker_color=reco_folded_colour, marker_size=0,
                                 subplot=mc_reco_hist_bin_reco_binning,
                                 normalise_hist=True),
                ]
                has_entries = [c.obj.GetEntries() > 0 for c in entries]
                if not any(has_entries):
                    print("Skipping 0 entries in", append, ibin_pt)
                    continue
                title = "%s\n%s region\n%g < p_{T}^{Gen} < %g GeV" % (jet_algo, region['label'], pt_bin_edges_gen[ibin_pt], pt_bin_edges_gen[ibin_pt+1])
                plot = Plot(entries,
                            what="hist",
                            title=title,
                            xtitle="Detector-level " + angle_str,
                            ytitle=ytitle,
                            subplot_type='ratio',
                            subplot_title='Data / MC',
                            subplot_limits=(0.8, 1.2))
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.72)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/detector_folded_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Folded, but only comparing data with data to check it is sane
                entries = [
                    Contribution(data_reco_hist_bin_reco_binning, label="Data (reco)",
                                 line_color=reco_data_colour, line_width=lw,
                                 marker_color=reco_data_colour, marker_size=0,
                                 normalise_hist=True),
                    Contribution(data_folded_hist_bin_reco_binning, label="Data (#tau = %.3g, folded)" % (tau),
                                 line_color=reco_folded_colour, line_width=lw,
                                 marker_color=reco_folded_colour, marker_size=0,
                                 subplot=data_reco_hist_bin_reco_binning,
                                 normalise_hist=True),
                ]
                has_entries = [c.obj.GetEntries() > 0 for c in entries]
                if not any(has_entries):
                    print("Skipping 0 entries in", append, ibin_pt)
                    continue
                title = "%s\n%s region\n%g < p_{T}^{Gen} < %g GeV" % (jet_algo, region['label'], pt_bin_edges_gen[ibin_pt], pt_bin_edges_gen[ibin_pt+1])
                plot = Plot(entries,
                            what="hist",
                            title=title,
                            xtitle="Detector-level " + angle_str,
                            ytitle=ytitle,
                            subplot_type='ratio',
                            subplot_title='Unfolded / reco',
                            subplot_limits=(0.8, 1.2))
                plot.legend.SetX1(0.6)
                plot.legend.SetY1(0.75)
                plot.legend.SetX2(0.98)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/detector_folded_only_data_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))

                # Both together, bit messy, not so useful
                # entries = [
                #     Contribution(gen_hist_bin, label="Generator",
                #                  line_color=ROOT.kBlue, line_width=lw,
                #                  marker_color=ROOT.kBlue, marker_size=0,
                #                  normalise_hist=True),
                #     Contribution(unfolded_hist_bin_total_errors, label="Unfolded (#tau = %.3g) (total err)" % (tau),
                #                  line_color=ROOT.kBlack, line_width=lw, line_style=1,
                #                  marker_color=ROOT.kBlack, marker_style=20, marker_size=0.75,
                #                  subplot=gen_hist_bin,
                #                  normalise_hist=True),
                #     Contribution(mc_reco_hist_bin_gen_binning, label="MC",
                #                  line_color=reco_mc_colour, line_width=lw,
                #                  marker_color=reco_mc_colour, marker_size=0,
                #                  normalise_hist=True),
                #     Contribution(data_reco_hist_bin_gen_binning, label="Data",
                #                  line_color=reco_data_colour, line_width=lw,
                #                  marker_color=reco_data_colour, marker_style=20, marker_size=0.75,
                #                  subplot=mc_reco_hist_bin_gen_binning,
                #                  normalise_hist=True),
                # ]
                # has_entries = [c.obj.GetEntries() > 0 for c in entries]
                # if not any(has_entries):
                #     print("Skipping 0 entries in", append, ibin_pt)
                #     continue
                # title = "%s\n%s region\n%g < p_{T}^{Gen} < %g GeV" % (jet_algo, region['label'], pt_bin_edges_gen[ibin_pt], pt_bin_edges_gen[ibin_pt+1])
                # plot = Plot(entries,
                #             what="hist",
                #             title=title,
                #             xtitle=angle_str,
                #             ytitle=ytitle,
                #             subplot_type='ratio',
                #             subplot_title='Data / MC',
                #             subplot_limits=(0.8, 1.2))
                # plot.legend.SetX1(0.6)
                # plot.legend.SetY1(0.68)
                # plot.legend.SetX2(0.98)
                # plot.legend.SetY2(0.9)
                # plot.plot("NOSTACK E1")
                # plot.save("%s/unfolded_detector_%s_bin_%d.%s" % (this_output_dir, append, ibin_pt, OUTPUT_FMT))


