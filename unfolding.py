#!/usr/bin/env python


"""TUnfold it all

Thanks to Ashley
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
        var_uf, var_of = True, False
        pt_uf, pt_of = True, False  # handle pt uder/over flow ourselves
        self.detector_binning = ROOT.TUnfoldBinning("detector")
        detector_distribution = self.detector_binning.AddBinning("detectordistribution")
        detector_distribution.AddAxis(self.variable_name, self.nbins_variable_reco, self.variable_bin_edges_reco, var_uf, var_of)
        detector_distribution.AddAxis("pt", self.nbins_pt_reco, self.pt_bin_edges_reco, pt_uf, pt_of)

        detector_distribution_underflow = self.detector_binning.AddBinning("detectordistribution_underflow")
        detector_distribution_underflow.AddAxis(self.variable_name, self.nbins_variable_reco, self.variable_bin_edges_reco, var_uf, var_of)
        detector_distribution_underflow.AddAxis("pt", self.nbins_pt_underflow_reco, self.pt_bin_edges_underflow_reco, pt_uf, pt_of)

        self.generator_binning = ROOT.TUnfoldBinning("generator")
        generator_distribution = self.generator_binning.AddBinning("generatordistribution")
        generator_distribution.AddAxis(self.variable_name, self.nbins_variable_gen, self.variable_bin_edges_gen, var_uf, var_of)
        generator_distribution.AddAxis("pt", self.nbins_pt_gen, self.pt_bin_edges_gen, pt_uf, pt_of)

        generator_distribution_underflow = self.generator_binning.AddBinning("generatordistribution_underflow")
        generator_distribution_underflow.AddAxis(self.variable_name, self.nbins_variable_gen, self.variable_bin_edges_gen, var_uf, var_of)
        generator_distribution_underflow.AddAxis("pt", self.nbins_pt_underflow_gen, self.pt_bin_edges_underflow_gen, pt_uf, pt_of)

        self.axisSteering = axisSteering
        self.unfolder = ROOT.TUnfoldDensity(response_map,
                                            orientation,
                                            regMode,
                                            constraintMode,
                                            densityFlags,
                                            self.generator_binning,
                                            self.detector_binning,
                                            # hmm comment these out when scanL or scanTau?
                                            "generator",
                                            self.axisSteering)

        self.use_axis_binning = True  # for things like get_probability_matrix()

    def setInput(self, hist):
        self.unfolder.SetInput(hist)

    def doScanTau(self, n_scan=300, scan_mode=ROOT.TUnfoldDensity.kEScanTauRhoMax):
        """Figure out best tau by scanning tau curve
        Taken from Ashley's code
        """

        # Graphs to save output from scan over tau
        scanresults = ROOT.MakeNullPointer(ROOT.TSpline)
        lCurve = ROOT.MakeNullPointer(ROOT.TGraph)
        LogTauX = ROOT.MakeNullPointer(ROOT.TSpline)
        LogTauY = ROOT.MakeNullPointer(ROOT.TSpline)

        tau_min, tau_max = 1E-15, 0.9

        iBestAvg = self.unfolder.ScanTau(n_scan,
                                         tau_min,
                                         tau_max,
                                         scanresults,
                                         scan_mode,
                                         "generator",
                                         self.axisSteering,
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
        leg1.SetFillStyle(0)
        leg1.SetBorderSize(0)
        leg1.SetTextSize(0.026)
        leg1.AddEntry(knots, '{} #rho: #tau = {}'.format(tauoptname,  tau), 'l')
        leg1.Draw()

        canv_tauScan.Modified()
        canv_tauScan.Update()
        canv_tauScan.Print("scantau_%s.%s" % (self.variable_name, self.OUTPUT_FMT))
        return tau


    def doScanL(self, n_scan=300):
        """Figure out best tau by doing scan over L curve
        Taken from Ashley's code
        """
        tau_min, tau_max = 0., 0.  # let it decide range for us
        tau_min, tau_max = 1E-14, 0.1
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

        canv_lScan.Print("scanL_%s.%s" % (self.variable_name, self.OUTPUT_FMT))
        return tau

    def do_unfolding(self, tau):
        print("Unfolding with tau =", tau)
        self.unfolder.DoUnfold(tau)
        # print("( " + str(self.unfolder.GetChi2A()) + " + " + str(self.unfolder.GetChi2L()) + ") / " + str(self.unfolder.GetNdf()))
        unfolded_1d = self.unfolder.GetOutput("unfolded_" + cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)  # use "generator" for signal + underflow region, "generatordistribution" for only signal region
        return unfolded_1d

    def get_bias(self):
        return self.unfolder.GetBias("bias_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)

    def get_ematrix_input(self):
        return self.unfolder.GetEmatrixInput("ematrix_input_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)

    def get_ematrix_sys_uncorr(self):
        return self.unfolder.GetEmatrixSysUncorr("ematrix_sys_uncorr_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)

    def get_ematrix_total(self):
        return self.unfolder.GetEmatrixTotal("ematrix_total_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)

    def get_rhoij_total(self):
        return self.unfolder.GetRhoIJtotal("rhoij_total_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)

    def get_probability_matrix(self):
        return self.unfolder.GetProbabilityMatrix("prob_matrix_"+cu.get_unique_str(), "", self.use_axis_binning)

    def get_var_hist_pt_binned(self, hist1d, ibin_pt, binning='generator'):
        """Get hist of variable for given pt bin from massive 1D hist that TUnfold makes"""
        binning = self.generator_binning.FindNode("generatordistribution") if binning == "generator" else self.detector_binning.FindNode("detectordistribution")
        var_bins = np.array(binning.GetDistributionBinning(0))
        pt_bins = np.array(binning.GetDistributionBinning(1))

        # need the -1 on ibin_pt, as it references an array index, whereas ROOT bins start at 1
        h = ROOT.TH1D("h_%d_%s" % (ibin_pt, cu.get_unique_str()), "", len(var_bins)-1, var_bins)
        for var_ind, var_value in enumerate(var_bins[:-1], 1):
            this_val = var_value * 1.001  # ensure its inside
            bin_num = binning.GetGlobalBinNumber(this_val, pt_bins[ibin_pt]*1.001)
            h.SetBinContent(var_ind, hist1d.GetBinContent(bin_num))
            h.SetBinError(var_ind, hist1d.GetBinError(bin_num))
            # print("Bin:", bin_num, this_val, pt_bins[ibin_pt], "=", hist1d.GetBinContent(bin_num), "+-", hist1d.GetBinError(bin_num))
        return h

    def draw_response_matrix(self, region_name, output_filename):
        canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
        canv.SetTicks(1, 1)
        canv.SetLogz()
        canv.SetRightMargin(1.5)
        canv.SetLeftMargin(0.9)
        rsp_map = self.response_map.Clone(cu.get_unique_str())
        rsp_map.SetTitle("Response matrix, %s region, %s;Generator Bin;Detector Bin" % (region_name, self.variable_name))
        rsp_map.GetYaxis().SetTitleOffset(1.5)
        rsp_map.GetXaxis().SetTitleOffset(1.5)
        rsp_map.Draw("COLZ")
        canv.SaveAs(output_filename)

    def draw_probability_matrix(self, region_name, output_filename):
        canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
        canv.SetTicks(1, 1)
        canv.SetLogz()
        canv.SetRightMargin(1.5)
        canv.SetLeftMargin(0.9)
        prob_map = self.get_probability_matrix()
        title = "Probability map, %s region, %s;Generator Bin;Detector Bin" % (region_name, self.variable_name)
        prob_map.SetTitle(title)
        prob_map.GetYaxis().SetTitleOffset(1.5)
        prob_map.GetXaxis().SetTitleOffset(1.5)
        prob_map.Draw("COLZ")
        canv.SaveAs(output_filename)

    def draw_error_matrix_input(self, region_name, output_filename):
        canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
        canv.SetTicks(1, 1)
        canv.SetLogz()
        canv.SetRightMargin(1.5)
        canv.SetLeftMargin(0.9)
        err_map = self.get_ematrix_input()
        title = "Error matrix (input), %s region, %s;Generator bin;Generator bin" % (region_name, self.variable_name)
        err_map.SetTitle(title)
        err_map.GetYaxis().SetTitleOffset(1.5)
        err_map.GetXaxis().SetTitleOffset(1.5)
        err_map.Draw("COLZ")
        canv.SaveAs(output_filename)

    def draw_error_matrix_sys_uncorr(self, region_name, output_filename):
        canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
        canv.SetTicks(1, 1)
        canv.SetLogz()
        canv.SetRightMargin(1.5)
        canv.SetLeftMargin(0.9)
        err_map = self.get_ematrix_sys_uncorr()
        title = "Error matrix (sys uncorr), %s region, %s;Generator bin;Generator bin" % (region_name, self.variable_name)
        err_map.SetTitle(title)
        err_map.GetYaxis().SetTitleOffset(1.5)
        err_map.GetXaxis().SetTitleOffset(1.5)
        err_map.Draw("COLZ")
        canv.SaveAs(output_filename)

    def draw_error_matrix_total(self, region_name, output_filename):
        canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
        canv.SetTicks(1, 1)
        canv.SetLogz()
        canv.SetRightMargin(1.5)
        canv.SetLeftMargin(0.9)
        err_map = self.get_ematrix_total()
        title = "Error matrix (total), %s region, %s;Generator bin;Generator bin" % (region_name, self.variable_name)
        err_map.SetTitle(title)
        err_map.GetYaxis().SetTitleOffset(1.5)
        err_map.GetXaxis().SetTitleOffset(1.5)
        err_map.Draw("COLZ")
        canv.SaveAs(output_filename)

    def draw_correlation_matrix(self, region_name, output_filename):
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

        canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
        canv.SetTicks(1, 1)
        # canv.SetLogz()
        canv.SetRightMargin(1.5)
        canv.SetLeftMargin(0.9)
        corr_map = self.get_rhoij_total()
        title = "Correlation map, %s region, %s;Generator Bin;Generator Bin" % (region_name, self.variable_name)
        corr_map.SetTitle(title)
        corr_map.GetYaxis().SetTitleOffset(1.5)
        corr_map.GetXaxis().SetTitleOffset(1.5)
        # ROOT.gStyle.SetPalette(ROOT.kLightTemperature)
        corr_map.SetMinimum(-1)
        corr_map.SetMaximum(1)
        corr_map.Draw("COLZ0")
        canv.SaveAs(output_filename)
        ROOT.gStyle.SetPalette(ROOT.kBird)



def plot_simple_unfolded(unfolded, reco, gen, output_filename, title=""):
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
    leg.AddEntry(unfolded, "Unfolded", "LP")

    hst.Draw("NOSTACK HISTE")
    leg.Draw()
    hst.SetMinimum(1E-2)
    hst.SetMaximum(1E12)
    canv_unfold.SaveAs(output_filename)


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


if __name__ == "__main__":
    input_mc_dy_tfile = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1RecoConstituents_V11JEC_JER_tUnfold/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root")
    input_singlemu_tfile = cu.open_root_file("workdir_ak4puppi_data_withAllResponses_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfold/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu.root")

    input_mc_qcd_tfile = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1RecoConstituents_V11JEC_JER_tUnfold/uhh2.AnalysisModuleRunner.MC.MC_QCD.root")
    # input_mc_qcd_tfile = cu.open_root_file("workdir_ak4puppi_pythia_newFlav_withAllResponses_jetAsymCut_chargedResp_pt1RecoConstituents_V11JEC_JER_tUnfold/uhh2.AnalysisModuleRunner.MC.MC_PYTHIA-QCD.root")
    input_jetht_tfile = cu.open_root_file("workdir_ak4puppi_data_withAllResponses_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfold/uhh2.AnalysisModuleRunner.DATA.Data_JetHTZeroBias.root")

    regions = [
        {
            "name": "Dijet",
            "dirname": "Dijet_QG_tighter",
            "label": "Dijet",
            "data_tfile": input_jetht_tfile,
            "mc_tfile": input_mc_qcd_tfile,
        },
        {
            "name": "ZPlusJets",
            "dirname": "ZPlusJets_QG",
            "label": "Z+jets",
            "data_tfile": input_singlemu_tfile,
            "mc_tfile": input_mc_dy_tfile,
        },
    ]

    output_dir = "unfolding"
    cu.check_dir_exists_create(output_dir)

    # TODO automate this
    jet_algo = "AK4 PF PUPPI"

    for region in regions[:]:

        # Setup pt bins
        # -------------
        # need different ones for Z+Jets region
        is_zpj = "ZPlusJets" in region['name']

        zpj_append = "_zpj" if is_zpj else ""

        pt_bin_edges_gen = qgc.PT_UNFOLD_DICT['signal%s_gen' % (zpj_append)]
        pt_bin_edges_reco = qgc.construct_fine_binning(pt_bin_edges_gen)

        pt_bin_edges_underflow_gen = qgc.PT_UNFOLD_DICT['underflow%s_gen' % (zpj_append)]
        pt_bin_edges_underflow_reco = qgc.construct_fine_binning(pt_bin_edges_underflow_gen)

        for angle in qgc.COMMON_VARS[2:3]:

            # Setup MyUnfolder object to do unfolding etc
            # -------------------------------------------
            angle_bin_edges_reco = qgc.VAR_UNFOLD_DICT[angle.var]['reco']
            angle_bin_edges_gen = qgc.VAR_UNFOLD_DICT[angle.var]['gen']
            angle_shortname = angle.var.replace("jet_", "")

            hist_data_reco = cu.get_from_tfile(region['data_tfile'], "%s/hist_%s_reco_new" % (region['dirname'], angle_shortname))
            hist_mc_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_new" % (region['dirname'], angle_shortname))
            hist_mc_gen = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_truth_new" % (region['dirname'], angle_shortname))
            hist_mc_gen_reco_map = cu.get_from_tfile(region['mc_tfile'], "%s/tu_%s_GenReco_new" % (region['dirname'], angle_shortname))

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

            # Set what is to be unfolded
            # ---------------------
            reco_1d = hist_data_reco
            reco_1d = hist_mc_reco
            unfolder.setInput(reco_1d)

            # Do any regularisation
            # ---------------------
            # tau = unfolder.doScanL()
            # tau = unfolder.doScanTau(300, ROOT.TUnfoldDensity.kEScanTauRhoAvgSys)
            # tau = 1E-10
            tau = 0

            # Do unfolding!
            # ---------------------
            unfolded_1d = unfolder.do_unfolding(tau)
            # stat errosrs only
            error_sys_uncorr_1d = make_hist_from_diagonals(unfolder.get_ematrix_sys_uncorr(), do_sqrt=True)


            # Draw unified unfolded distributions
            # ---------------------
            plot_simple_unfolded(unfolded=unfolded_1d,
                                 reco=reco_1d,
                                 gen=hist_mc_gen,
                                 output_filename="%s/unfolded_%s_%s.%s" % (output_dir, region['name'], angle.var, OUTPUT_FMT),
                                 title="%s region, %s" % (region['label'], angle.name))

            # Draw matrices
            # -------------
            unfolder.draw_response_matrix(region['name'], "%s/response_map_%s_%s.%s" % (output_dir, region['name'], angle.var, OUTPUT_FMT))
            unfolder.draw_probability_matrix(region['name'], "%s/probability_map_%s_%s.%s" % (output_dir, region['name'], angle.var, OUTPUT_FMT))
            unfolder.draw_correlation_matrix(region['name'], "%s/corr_map_%s_%s.%s" % (output_dir, region['name'], angle.var, OUTPUT_FMT))
            unfolder.draw_error_matrix_input(region['name'], "%s/err_map_sys_input_%s_%s.%s" % (output_dir, region['name'], angle.var, OUTPUT_FMT))
            unfolder.draw_error_matrix_sys_uncorr(region['name'], "%s/err_map_sys_uncorr_%s_%s.%s" % (output_dir, region['name'], angle.var, OUTPUT_FMT))
            unfolder.draw_error_matrix_total(region['name'], "%s/err_map_total_%s_%s.%s" % (output_dir, region['name'], angle.var, OUTPUT_FMT))

            # continue

            # Draw individual pt bin plots
            # ----------------------------
            for ibin_pt in range(0, len(pt_bin_edges_gen)-1):

                gen_hist_bin = unfolder.get_var_hist_pt_binned(hist_mc_gen, ibin_pt, binning="generator")
                unfolded_hist_bin = unfolder.get_var_hist_pt_binned(unfolded_1d, ibin_pt, binning="generator")
                unfolded_hist_bin_errors = unfolder.get_var_hist_pt_binned(error_sys_uncorr_1d, ibin_pt, binning="generator")
                update_hist_bin_content(unfolded_hist_bin, unfolded_hist_bin_errors)

                for n in range(1, gen_hist_bin.GetNbinsX()+1):
                    print("Bin", n)
                    print("gen_hist:", gen_hist_bin.GetBinContent(n), "+-", gen_hist_bin.GetBinError(n))
                    print("unfolded_hist:", unfolded_hist_bin.GetBinContent(n), "+-", unfolded_hist_bin.GetBinError(n))
                    print("unfolded_hist_bin_errors:", unfolded_hist_bin_errors.GetBinContent(n), "+-", unfolded_hist_bin_errors.GetBinError(n))

                entries = [
                    Contribution(gen_hist_bin, label="Generator",
                                 line_color=ROOT.kBlue, line_width=2,
                                 marker_color=ROOT.kBlue, marker_size=0,
                                 normalise_hist=True),
                    Contribution(unfolded_hist_bin, label="Unfolded",
                                 line_color=ROOT.kRed, line_width=1,
                                 marker_color=ROOT.kRed, marker_style=20, marker_size=0,
                                 subplot=gen_hist_bin,
                                 normalise_hist=True),
                    Contribution(unfolded_hist_bin_errors, label="Unfolded (stat err)",
                                 line_color=ROOT.kGreen+2, line_width=1,
                                 marker_color=ROOT.kGreen+2, marker_style=20, marker_size=0,
                                 # subplot=gen_hist_bin,
                                 normalise_hist=True),
                ]
                has_entries = [c.obj.GetEntries() > 0 for c in entries]
                if not any(has_entries):
                    print("Skipping 0 entries in", region['name'], angle.var, ibin_pt)
                    continue
                title = "%s\n%s region\n%g < p_{T}^{Gen} < %g GeV" % (jet_algo, region['label'], pt_bin_edges_gen[ibin_pt], pt_bin_edges_gen[ibin_pt+1])
                plot = Plot(entries,
                            what="hist",
                            title=title,
                            xtitle=angle.name,
                            ytitle='p.d.f.',
                            subplot_type='ratio',
                            subplot_title='Unfolded / gen',
                            subplot_limits=(0.8, 1.2))
                plot.legend.SetY1(0.75)
                plot.legend.SetY2(0.9)
                plot.plot("NOSTACK E1")
                plot.save("%s/unfolded_%s_%s_bin_%d.%s" % (output_dir, region['name'], angle.var, ibin_pt, OUTPUT_FMT))


