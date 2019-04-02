#!/usr/bin/env python

"""Use TUnfold it all"""

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
                 variable_bin_edges, # for e.g. ptD, LHA
                 variable_bin_edges_coarse,
                 variable_name,
                 pt_bin_edges,
                 pt_bin_edges_coarse,
                 orientation=ROOT.TUnfold.kHistMapOutputHoriz,
                 constraintMode=ROOT.TUnfold.kEConstraintArea,
                 regMode=ROOT.TUnfold.kRegModeCurvature,
                 densityFlags=ROOT.TUnfoldDensity.kDensityModeBinWidth,
                 axisSteering='*[UuOob]'):

        self.response_map = response_map
        self.variable_name = variable_name

        self.variable_bin_edges = variable_bin_edges
        self.nbins_variable = len(variable_bin_edges)-1
        self.variable_bin_edges_coarse = variable_bin_edges_coarse
        self.nbins_variable_coarse = len(variable_bin_edges_coarse)-1
        self.pt_bin_edges = pt_bin_edges
        self.nbins_pt = len(pt_bin_edges)-1
        self.pt_bin_edges_coarse = pt_bin_edges_coarse
        self.nbins_pt_coarse = len(pt_bin_edges_coarse)-1

        self.detector_binning = ROOT.TUnfoldBinning("detector")
        detector_distribution_LHA = self.detector_binning.AddBinning("detectordistribution")
        detector_distribution_LHA.AddAxis(self.variable_name, self.nbins_variable, self.variable_bin_edges, False, False)
        detector_distribution_LHA.AddAxis("pt", self.nbins_pt, self.pt_bin_edges, False, True)

        self.generator_binning = ROOT.TUnfoldBinning("generator")
        generator_distribution_LHA = self.generator_binning.AddBinning("generatordistribution")
        generator_distribution_LHA.AddAxis(self.variable_name, self.nbins_variable_coarse, self.variable_bin_edges_coarse, False, False)
        generator_distribution_LHA.AddAxis("pt", self.nbins_pt_coarse, self.pt_bin_edges_coarse, True, True)

        self.axisSteering = axisSteering
        self.unfolder = ROOT.TUnfoldDensity(response_map, orientation, regMode, constraintMode, densityFlags,
                                            self.generator_binning, self.detector_binning,
                                            "generatordistribution", axisSteering)

        self.unfolded = None  # for later

    def setInput(self, hist):
        self.unfolder.SetInput(hist)

    def doScanTau(self, n_scan=300, scan_mode=ROOT.TUnfoldDensity.kEScanTauRhoMax):
        """Figure out best tau by scanning tau curve"""

        # Graphs to save output from scan over tau
        scanresults = ROOT.MakeNullPointer(ROOT.TSpline)
        lCurve = ROOT.MakeNullPointer(ROOT.TGraph)
        LogTauX = ROOT.MakeNullPointer(ROOT.TSpline)
        LogTauY = ROOT.MakeNullPointer(ROOT.TSpline)

        tau_min, tau_max = 0., 0.  # let it decide range for us

        iBestAvg = self.unfolder.ScanTau(n_scan,
                                         tau_min,
                                         tau_max,
                                         scanresults,
                                         scan_mode,
                                         "generatordistribution",
                                         axisSteering,
                                         lCurve,
                                         LogTauX,
                                         LogTauY)
        tau = self.unfolder.GetTau()

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
        int1 = 1
        bestRhoLogTau = ROOT.TGraph(int1, t, rho)
        print("10^log_10(tau) = tau = {}".format(math.pow(10., float(t0))))
        bestLCurve = ROOT.TGraph(int(1), x, y)

        tAll, rhoAll = array('d'), array('d')
        for i in range(n_scan):
            tt = ROOT.Double(0.0)
            rr = ROOT.Double(0.0)
            scanresults.GetKnot(i, tt, rr)
            tAll.append(tt)
            rhoAll.append(rr)

        knots = ROOT.TGraph(int(n_scan), tAll, rhoAll)
        print("chi**2 A {:3.1f} + chi**2 L {:3.1f} / NDOF {:3.1f} ".format(self.unfolder.GetChi2A(),
                                                                           self.unfolder.GetChi2L(),
                                                                           self.unfolder.GetNdf()))

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
        leg1.SetBorderSize(0)
        leg1.SetTextSize(0.026)
        leg1.AddEntry(knots, '{} #rho: #tau = {}'.format(tauoptname,  tau), 'l')
        leg1.Draw()

        canv_tauScan.Modified()
        canv_tauScan.Update()
        canv_tauScan.Print("scantau_%s.%s" % (self.variable_name, self.OUTPUT_FMT))
        return tau


    def doScanL(self, n_scan=300):
        """Figure out best tau by doing scan over L curve"""
        tau_min, tau_max = 0., 0.  # let it decide range for us
        scannedlcurve = ROOT.MakeNullPointer(ROOT.TGraph)
        logTauX = ROOT.MakeNullPointer(ROOT.TSpline3)
        logTauY = ROOT.MakeNullPointer(ROOT.TSpline3)
        logTauCurvature = ROOT.MakeNullPointer(ROOT.TSpline3)

        best = self.unfolder.ScanLcurve(n_scan,
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
        tau = self.unfolder.GetTau()
        print("doScanL value is {}".format(tau))
        print("Returned the TUnfoldDensity object.")

        canv_lScan = ROOT.TCanvas("canv_lScan"+str(tau), "canv_lScan"+str(tau))

        bestLcurve.SetLineColor(ROOT.kBlue+2)

        bestLogTauLogChi2.SetMarkerColor(ROOT.kRed)

        bestLcurve.Draw("")
        bestLogTauLogChi2.Draw("* same")
        bestLcurve.GetXaxis().SetTitle("log_{10}(L_{1})")
        bestLcurve.GetYaxis().SetTitle("log_{10}(#frac{L_{2}}{#tau^2})")  # frac{Theory}{Unfolded MC}

        bestLcurve.SetTitle("Optimization of Regularization Parameter, #tau : Scan of L Curve")

        leg1 = ROOT.TLegend(0.5, 0.6, 0.7, 0.94)
        leg1.SetFillColor(0)
        leg1.SetBorderSize(0)
        leg1.SetTextSize(0.026)
        leg1.AddEntry(bestLcurve, 'L Curve Scan: #tau = {}'.format(tau), 'l')
        leg1.Draw()

        canv_lScan.Modified()
        canv_lScan.Update()

        canv_lScan.Print("scanL_%s.%s" % (self.variable_name, self.OUTPUT_FMT))
        return tau

    def doUnfolding(self, tau):
        self.unfolder.DoUnfold(tau)

        print("( " + str(self.unfolder.GetChi2A()) + " + " + str(self.unfolder.GetChi2L()) + ") / " + str(self.unfolder.GetNdf()))
        print("Tau : ", tau)
        print("Tau : ", self.unfolder.GetTau())
        self.unfolded = self.unfolder.GetOutput("unfolded", "", "generator", self.axisSteering, False)
        # FIXME: do errors properly, not nullptr
        self.unfolded_2d = self.generator_binning.ExtractHistogram("unfolded2D", self.unfolded, ROOT.MakeNullPointer(ROOT.TH2), True, self.axisSteering)
        return self.unfolded

    def get_unfolded_var_hist_pt_binned(self, ibin_pt):
        """Get 1D histogram of our unfolded variable for a given pT bin

        Parameters
        ----------
        ibin_pt : int
            pT bin #
        """
        h_2d = self.unfolded_2d
        ax_y = h_2d.GetYaxis()
        pt_low = ax_y.GetBinLowEdge(ibin_pt)
        pt_high = ax_y.GetBinLowEdge(ibin_pt+1)
        hname = "h_unfolded_%d" % (ibin_pt)
        h_1d = h_2d.ProjectionX(hname, ibin_pt, ibin_pt+1)
        h_1d.SetName(hname)
        h_1d.SetTitle("%g < p_{T}^{Gen} < %g GeV;%s;N" % (pt_low, pt_high, self.variable_name))
        return h_1d


def plot_simple_unfolded(unfolded, reco, gen, output_filename):
    """Simple plot of unfolded, reco, gen, by bin number (ie non physical axes)"""
    canv_unfold = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv_unfold.SetLogy()
    canv_unfold.SetTicks(1, 1)
    leg = ROOT.TLegend(0.7, 0.7, 0.88, 0.88)
    hst = ROOT.THStack("hst", ";Bin Number;N")

    gen.SetLineColor(ROOT.kBlue)
    hst.Add(gen)
    leg.AddEntry(gen, "Gen", "L")

    unfolded.SetLineColor(ROOT.kRed)
    unfolded.SetLineWidth(0)
    unfolded.SetMarkerColor(ROOT.kRed)
    unfolded.SetMarkerSize(0.6)
    unfolded.SetMarkerStyle(20)
    hst.Add(unfolded)
    leg.AddEntry(unfolded, "Unfolded", "LP")

    reco.SetLineColor(ROOT.kGreen+2)
    hst.Add(reco)
    leg.AddEntry(reco, "Reco", "L")

    hst.Draw("NOSTACK HISTE")
    leg.Draw()
    hst.SetMinimum(1E-6)
    hst.SetMaximum(1E12)
    canv_unfold.Draw()
    canv_unfold.SaveAs(output_filename)



if __name__ == "__main__":
    input_mc_qcd_tfile = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1Constituents_V11JEC_JER_tUnfold/uhh2.AnalysisModuleRunner.MC.MC_QCD.root")
    input_mc_dy_tfile = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1Constituents_V11JEC_JER_tUnfold/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root")
    input_jetht_tfile = cu.open_root_file("workdir_ak4puppi_data_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1Constituents_V11JEC_JER_tUnfold/uhh2.AnalysisModuleRunner.DATA.Data_JetHTZeroBias.root")
    input_singlemu_tfile = cu.open_root_file("workdir_ak4puppi_data_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1Constituents_V11JEC_JER_tUnfold/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu.root")

    # for region in ["Dijet", "ZPlusJets"]:
        # for variable in qgc.COMMON_VARS:


    hist_mc_gen_reco_map_LHA = cu.get_from_tfile(input_mc_qcd_tfile, "Dijet_QG_tighter/histLHAGenRecnew")
    hist_mc_reco_LHA = cu.get_from_tfile(input_mc_qcd_tfile, "Dijet_QG_tighter/histLHAReconew")
    hist_mc_gen_LHA = cu.get_from_tfile(input_mc_qcd_tfile, "Dijet_QG_tighter/histLHATruthnew")

    lha_bin_edges = np.array([0.0, 0.29, 0.37, 0.44, 0.5, 0.56, 0.62, 0.68, 0.75, 1.0])
    lha_bin_edges_coarse = np.array([0.0, 0.37, 0.5, 0.62, 0.75, 1.0])

    pt_bin_edges = np.array([0, 29, 38, 50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 1000, 1300, 1700, 2200, 3000, 4000, 5000, 10000], dtype='d')
    pt_bin_edges_coarse = np.array([0, 29, 50, 88, 150, 254, 408, 614, 1000, 1700, 3000, 5000, 10000], dtype='d')

    unfolder_LHA = MyUnfolder(response_map=hist_mc_gen_reco_map_LHA,
                              variable_bin_edges=lha_bin_edges,
                              variable_bin_edges_coarse=lha_bin_edges_coarse,
                              variable_name="LHA",
                              pt_bin_edges=pt_bin_edges,
                              pt_bin_edges_coarse=pt_bin_edges_coarse,
                              orientation=ROOT.TUnfold.kHistMapOutputHoriz,
                              constraintMode=ROOT.TUnfold.kEConstraintArea,
                              regMode=ROOT.TUnfold.kRegModeCurvature,
                              densityFlags=ROOT.TUnfoldDensity.kDensityModeBinWidth,
                              axisSteering='*[b]')
    unfolder_LHA.setInput(hist_mc_reco_LHA)
    # tau = unfolder_LHA.doScanL()
    tau = 0
    unfolded = unfolder_LHA.doUnfolding(tau)
    plot_simple_unfolded(unfolded, hist_mc_reco_LHA, hist_mc_gen_LHA, "unfolded_LHA.%s" % (OUTPUT_FMT))
    # unfolder_LHA.plot_()

    for ibin_pt in range(1, len(pt_bin_edges_coarse)):
        unfolded_LHA_bin = unfolder_LHA.get_unfolded_var_hist_pt_binned(ibin_pt)
        gen_LHA_bin = unfolded_LHA_bin.Clone(cu.get_unique_str())
        entries = [
            Contribution(gen_LHA_bin, label="Generator",
                         line_color=ROOT.kBlue, line_width=2,
                         marker_color=ROOT.kBlue),
            Contribution(unfolded_LHA_bin, label="Unfolded",
                         line_color=ROOT.kRed, line_width=0,
                         marker_color=ROOT.kRed, marker_style=20,
                         subplot=gen_LHA_bin),
        ]
        plot = Plot(entries, "hist", title=unfolded_LHA_bin.GetTitle(), subplot_type='ratio', subplot_title='Unfolded / gen')
        plot.legend.SetY1(0.75)
        plot.legend.SetY2(0.9)
        plot.plot("NOSTACK HISTE")
        plot.save("unfolded_LHA_bin_%d.pdf" % (ibin_pt))


