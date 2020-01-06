#!/usr/bin/env python

"""
Classes to be used in main unfolding script/routines
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
        self.tau = 0
        self.scan_mode = ""
        self.graph_all_scan_points = None
        self.graph_best_scan_point = None

    def scan_tau(self, tunfolder, n_scan, tau_min, tau_max, scan_mode, distribution, axis_steering):
        ind_best_point = tunfolder.ScanTau(n_scan,
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
        self._process_results(scan_mode, ind_best_point)
        print("scan_tau value is {}".format(self.tau))
        print("chi**2 A {:3.1f} + chi**2 L {:3.1f} / NDOF {:3.1f} ".format(tunfolder.GetChi2A(),
                                                                           tunfolder.GetChi2L(),
                                                                           tunfolder.GetNdf()))
        return self.tau

    def _process_results(self, scan_mode, ind_best_point):
        """Create graphs etc from ScanTau output

        User shouldn't call this, only internal
        """

        # Get best scan point & make graph of it
        # t here is log_10(tau)
        t, rho = array('d'), array('d')  # array obj needed to make TGraph
        t0 = ROOT.Double(0.0)
        rho0 = ROOT.Double(0.0)
        self.scan_results.GetKnot(ind_best_point, t0, rho0)
        t.append(t0)
        rho.append(rho0)
        self.graph_best_scan_point = ROOT.TGraph(1, t, rho)

        print("t[0] =", t[0])
        print("rho[0] =", rho[0])
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
        canv_tau_scan = ROOT.TCanvas("canv_tau_scan_"+str(self.tau), "canv_tau_scan_"+str(self.tau))

        self.graph_all_scan_points.SetLineColor(ROOT.kBlue+3)
        self.graph_all_scan_points.Draw()

        self.graph_best_scan_point.SetMarkerColor(ROOT.kRed)
        self.graph_best_scan_point.Draw("* same")

        self.graph_all_scan_points.GetXaxis().SetTitle("log_{10}(#tau)")
        self.graph_all_scan_points.GetYaxis().SetTitle(" #rho")

        leg = ROOT.TLegend(0.2, 0.6, 0.35, 0.89)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.026)
        leg.AddEntry(self.graph_all_scan_points, 'Scan over #tau', 'l')
        leg.AddEntry(self.graph_best_scan_point, 'Chosen point: #tau = {}'.format(self.tau), 'P')
        leg.Draw()

        canv_tau_scan.Print(output_filename)

    def save_to_tfile(self, tfile):
        tfile.WriteTObject(self.graph_all_scan_points, "regularize_all_scan_points")
        tfile.WriteTObject(self.graph_best_scan_point, "regularize_best_scan_point")


class LCurveScanner(object):
    """Class to handle doing ScanLcurve on a TUnfoldBinning object,
    since it produces many associated objects & values.

    Can also plot results from scanning.
    """
    def __init__(self):
        self.scanned_l_curve = ROOT.MakeNullPointer(ROOT.TGraph)
        self.log_tau_x = ROOT.MakeNullPointer(ROOT.TSpline3) # spline of L-curve x-coord as a func of log_10(tau)
        self.log_tau_y = ROOT.MakeNullPointer(ROOT.TSpline3) # spline of L-curve y-coord as a func of log_10(tau)
        self.log_tau_curvature = ROOT.MakeNullPointer(ROOT.TSpline3)
        self.graph_log_tau_curvature = None  # to hold graph of log_tau_curvature
        self.graph_log_tau_curvature_best = None
        self.tau = 0
        self.graph_best_scan_point = None # in terms of LcurveY vs LcurveX

    def scan_L(self, tunfolder, n_scan, tau_min, tau_max):
        ind_best_point = tunfolder.ScanLcurve(n_scan,
                                              tau_min,
                                              tau_max,
                                              self.scanned_l_curve,
                                              self.log_tau_x,
                                              self.log_tau_y,
                                              self.log_tau_curvature)
        self.tau = tunfolder.GetTau()
        self._process_results(ind_best_point)
        return self.tau

    def _process_results(self, ind_best_point):
        """Create graphs etc from ScanLcurve output

        User shouldn't call this, only internal
        """
        # Get best scan point & make graph of it
        t_0 = ROOT.Double(0.0)  # is log_10(tau)
        x_0 = ROOT.Double(0.0)
        y_0 = ROOT.Double(0.0)
        self.log_tau_x.GetKnot(ind_best_point, t_0, x_0)
        self.log_tau_y.GetKnot(ind_best_point, t_0, y_0)
        self.graph_best_scan_point = ROOT.TGraph(1, array('d', [x_0]), array('d', [y_0]))

        # Create graph of curvature
        t_all, c_all = array('d'), array('d')
        n_scan = self.log_tau_curvature.GetNp()
        for i in range(n_scan):
            t = ROOT.Double(0.0)  # is log_10(tau)
            c = ROOT.Double(0.0)
            self.log_tau_curvature.GetKnot(i, t, c)
            t_all.append(t)
            c_all.append(c)

        self.graph_log_tau_curvature = ROOT.TGraph(n_scan, t_all, c_all)

        # Get best scan point in terms of curvature vs log(tau)
        # you cannot use the index, it doesn't correspond to this graph
        c_0 = self.log_tau_curvature.Eval(t_0)
        self.graph_log_tau_curvature_best = ROOT.TGraph(1, array('d', [t_0]), array('d', [c_0]))

    def plot_scan_L_curve(self, output_filename):
        """Plot graph of scan results, and optimum tau"""
        canv_L_scan = ROOT.TCanvas("canv_L_scan_"+str(self.tau), "canv_L_scan_"+str(self.tau))

        self.scanned_l_curve.SetTitle("Optimization of Regularization Parameter, #tau : Scan of L curve")
        self.scanned_l_curve.SetLineColor(ROOT.kBlue+3)
        self.scanned_l_curve.Draw()

        self.graph_best_scan_point.SetMarkerColor(ROOT.kRed)
        self.graph_best_scan_point.Draw("* same")

        self.scanned_l_curve.GetXaxis().SetTitle("log_{10}(L_{1})")
        self.scanned_l_curve.GetYaxis().SetTitle("log_{10}(#frac{L_{2}}{#tau^{2}})")

        leg = ROOT.TLegend(0.5, 0.6, 0.85, 0.89)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.026)
        leg.AddEntry(self.scanned_l_curve, 'Scan over #tau', 'l')
        leg.AddEntry(self.graph_best_scan_point, 'Chosen point: #tau = {}'.format(self.tau), 'P')
        leg.Draw()

        canv_L_scan.Print(output_filename)

    def plot_scan_L_curvature(self, output_filename):
        """Plot graph of L curvature & optimum point"""
        canv_L_curvature = ROOT.TCanvas("canv_L_curvature_"+str(self.tau), "canv_L_curvature_"+str(self.tau))

        self.graph_log_tau_curvature.SetTitle("Optimization of Regularization Parameter, #tau : Scan of L curvature")
        self.graph_log_tau_curvature.SetLineColor(ROOT.kBlue+3)
        self.graph_log_tau_curvature.Draw()
        self.graph_log_tau_curvature.GetXaxis().SetTitle("log_{10}(#tau)")
        self.graph_log_tau_curvature.GetYaxis().SetTitle("L-curve curvature C")

        self.graph_log_tau_curvature_best.SetLineColor(ROOT.kRed)
        self.graph_log_tau_curvature_best.SetMarkerColor(ROOT.kRed)
        self.graph_log_tau_curvature_best.Draw("* same")

        leg = ROOT.TLegend(0.5, 0.6, 0.85, 0.89)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.026)
        leg.AddEntry(self.graph_log_tau_curvature, 'Curvature', 'l')
        leg.AddEntry(self.graph_log_tau_curvature_best, 'Chosen point: #tau = {}'.format(self.tau), 'P')
        leg.Draw()

        canv_L_curvature.Print(output_filename)

    def save_to_tfile(self, tfile):
        tfile.WriteTObject(self.scanned_l_curve, "scanned_l_curve")
        tfile.WriteTObject(self.graph_best_scan_point, "graph_best_scan_point")
        tfile.WriteTObject(self.graph_log_tau_curvature, "graph_log_tau_curvature")
        tfile.WriteTObject(self.graph_log_tau_curvature_best, "graph_log_tau_curvature_best")


class MyUnfolder(object):
    """Main class to handle unfolding input/outputs, all the associated objects"""


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
                 distribution='generatordistribution',
                 axisSteering='*[b]'):

        self.response_map = response_map
        self.response_map_matrix = self.th2_to_tmatrixd(response_map)
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
        self.distribution = distribution
        self.axisSteering = axisSteering

        self.tunfolder = ROOT.TUnfoldDensity(response_map,
                                             self.orientation,
                                             self.regMode,
                                             self.constraintMode,
                                             self.densityFlags,
                                             self.generator_binning,
                                             self.detector_binning,
                                             # hmm these take preference over whatever is use for scantau?
                                             self.distribution,
                                             self.axisSteering)

        self.use_axis_binning = False  # for things like get_probability_matrix()

        # hists that will be assigned later

        # reco, to be unfolded (data or MC)
        self.input_hist = None
        self.input_hist_bg_subtracted = None

        # reco, but using gen binning
        self.input_hist_gen_binning = None
        self.input_hist_gen_binning_bg_subtracted = None

        # fakes
        # for total bg, use get_total_background(),
        # or get_total_background_gen_binning()

        # reco MC
        self.hist_mc_reco = None
        self.hist_mc_reco_bg_subtracted = None

        # reco MC, but using gen binning
        self.hist_mc_reco_gen_binning = None
        self.hist_mc_reco_gen_binning_bg_subtracted = None

        # generator-level MC truth
        self.hist_truth = None  # gen truth

        self.tau = 0  # to be set by user later, via TauScanner or LCurveScanner
        self.backgrounds = {}  # gets filled with subtract_background()
        self.backgrounds_gen_binning = {}  # gets filled with subtract_background_gen_binning()

        self.unfolded = None  # set in get_output()
        self.unfolded_stat_err = None  # set in get_unfolded_with_ematrix_stat()

        self.syst_maps = {}  # gets filled with add_sys_error()
        self.syst_shifts = {}  # gets filled with get_sys_shift()

        # use "generator" for signal + underflow region, "generatordistribution" for only signal region
        self.output_distribution_name = "generator"

        self.folded_unfolded = None  # set in get_folded_unfolded()
        self.folded_mc_truth = None  # set in get_folded_mc_truth()

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

    @staticmethod
    def _check_save_to_tfile(tfile, obj, name=None):
        if obj:
            if name is None:
                name = obj.GetName()
            tfile.WriteTObject(obj, name)
        else:
            print("Not saving", name, "as does not exist")

    def save_to_tfile(self, tfile):
        """Save important stuff to TFile/TDirectory"""
        self._check_save_to_tfile(tfile, self.detector_binning)
        self._check_save_to_tfile(tfile, self.generator_binning)
        self._check_save_to_tfile(tfile, self.response_map, "response_map")

        # all 1D input gen/reco etc hists
        self._check_save_to_tfile(tfile, self.input_hist, "input_hist")
        self._check_save_to_tfile(tfile, self.input_hist_bg_subtracted, "input_hist_bg_subtracted")
        self._check_save_to_tfile(tfile, self.input_hist_gen_binning, "input_hist_gen_binning")
        self._check_save_to_tfile(tfile, self.input_hist_gen_binning_bg_subtracted, "input_hist_gen_binning_bg_subtracted")
        self._check_save_to_tfile(tfile, self.hist_truth, "hist_truth")
        self._check_save_to_tfile(tfile, self.hist_mc_reco, "hist_mc_reco")
        self._check_save_to_tfile(tfile, self.hist_mc_reco_bg_subtracted, "hist_mc_reco_bg_subtracted")
        self._check_save_to_tfile(tfile, self.hist_mc_reco_gen_binning, "hist_mc_reco_gen_binning")
        self._check_save_to_tfile(tfile, self.hist_mc_reco_gen_binning_bg_subtracted, "hist_mc_reco_gen_binning_bg_subtracted")

        # save all backgrounds (incl fakes)
        for name, hist in self.backgrounds.items():
            self._check_save_to_tfile(tfile, hist, "background_reco_binning_%s" % name.replace(" ", "_"))
        for name, hist in self.backgrounds_gen_binning.items():
            self._check_save_to_tfile(tfile, hist, "background_gen_binning_%s" % name.replace(" ", "_"))

        # save error matrices
        self._check_save_to_tfile(tfile, self.ematrix_input, "ematrix_input")
        self._check_save_to_tfile(tfile, self.ematrix_sys_uncorr, "ematrix_sys_uncorr")
        self._check_save_to_tfile(tfile, self.ematrix_total, "ematrix_total")
        self._check_save_to_tfile(tfile, self.ematrix_stat_sum, "ematrix_stat_sum")

        # save other matrices
        self._check_save_to_tfile(tfile, self.rhoij_total, "rhoij_total")
        # self._check_save_to_tfile(tfile, self.covariance_matrix, "covariance_matrix")
        self._check_save_to_tfile(tfile, self.probability_matrix, "probability_matrix")

        # save systematic response matrices
        for name, syst_map in self.syst_maps.items():
            self._check_save_to_tfile(tfile, syst_map, "syst_map_%s" % name.replace(" ", "_"))

        # save systematic shifts
        for name, syst_shift in self.syst_shifts.items():
            self._check_save_to_tfile(tfile, syst_shift, "syst_shift_%s" % name.replace(" ", "_"))

        # Folded things
        self._check_save_to_tfile(tfile, self.folded_unfolded, "folded_unfolded")
        self._check_save_to_tfile(tfile, self.folded_mc_truth, "folded_mc_truth")

        # Save unfolded
        self._check_save_to_tfile(tfile, self.unfolded, "unfolded")
        self._check_save_to_tfile(tfile, self.unfolded_stat_err, "unfolded_stat_err")

    def set_input(self,
                  input_hist,
                  input_hist_gen_binning=None,
                  hist_truth=None,
                  hist_mc_reco=None,
                  hist_mc_reco_bg_subtracted=None,
                  hist_mc_reco_gen_binning=None,
                  hist_mc_reco_gen_binning_bg_subtracted=None,
                  bias_factor=0):
        """Set hist to be unfolded, plus other basic hists

        Also allow other args to be passed to TUnfoldSys::SetInput
        """
        self.input_hist = input_hist.Clone()
        self.input_hist_bg_subtracted = input_hist.Clone()
        self.input_hist_gen_binning = input_hist_gen_binning.Clone()
        self.input_hist_gen_binning_bg_subtracted = input_hist_gen_binning.Clone()
        self.hist_truth = hist_truth
        self.hist_mc_reco = hist_mc_reco
        self.hist_mc_reco_bg_subtracted = hist_mc_reco_bg_subtracted
        self.hist_mc_reco_gen_binning = hist_mc_reco_gen_binning
        self.hist_mc_reco_gen_binning_bg_subtracted = hist_mc_reco_gen_binning_bg_subtracted
        self.tunfolder.SetInput(input_hist, bias_factor)

    def subtract_background(self, hist, name, scale=1.0, scale_err=0.0):
        """Subtract background source from input hist"""
        # Save into dict of components - needed? since TUnfoldDensity does this as well
        self.backgrounds[name] = hist.Clone()
        self.backgrounds[name].Scale(scale)
        # Also save total input subtracted
        self.input_hist_bg_subtracted.Add(hist, -1*scale)
        self.tunfolder.SubtractBackground(hist.Clone(), name, scale, scale_err)

    def subtract_background_gen_binning(self, hist, name, scale=1.0, scale_err=0.0):
        """Subtract background source with gen binning from input hist

        NB doesn't affect TUnfold, only for own book keeping
        """
        # Save into dict of components
        self.backgrounds_gen_binning[name] = hist.Clone()
        self.backgrounds_gen_binning[name].Scale(scale)
        # Also save total input subtracted
        self.input_hist_gen_binning_bg_subtracted.Add(hist, -1*scale)

    def get_total_background(self):
        """Get total cumulative background"""
        total_bg_hist = None
        for name, hist in self.backgrounds.items():
            if total_bg_hist is None:
                total_bg_hist = hist.Clone()
            else:
                total_bg_hist.Add(hist)
        return total_bg_hist

    def get_total_background_gen_binning(self):
        """Get total cumulative background with generator binning"""
        total_bg_hist = None
        for name, hist in self.backgrounds_gen_binning.items():
            if total_bg_hist is None:
                total_bg_hist = hist.Clone()
            else:
                total_bg_hist.Add(hist)
        return total_bg_hist

    def do_unfolding(self, tau):
        """Carry out unfolding with given regularisastion parameter"""
        print(">>> Unfolding with tau =", tau)
        self.tau = tau
        self.tunfolder.DoUnfold(tau)

    def add_sys_error(self, map_syst, name, mode):
        """Add systematic error via response map, arguments as per AddSysError()"""
        self.syst_maps[name] = map_syst
        self.tunfolder.AddSysError(map_syst, name, self.orientation, ROOT.TUnfoldDensity.kSysErrModeMatrix)

    def get_delta_sys_shift(self, syst_label, hist_name, hist_title=""):
        if syst_label not in self.syst_maps:
            raise KeyError("No systematic %s, only have: %s" % (syst_label, ", ".join(self.syst_maps.keys())))
        hist = self.tunfolder.GetDeltaSysSource(syst_label,
                                                hist_name,
                                                hist_title,
                                                self.output_distribution_name, # must be the same as what's used in get_output
                                                self.axisSteering,
                                                self.use_axis_binning)
        self.syst_shifts[syst_label] = hist
        return hist

    def get_output(self, hist_name='unfolded', update_with_ematrix_total=False):
        """Get 1D unfolded histogram covering all bins"""
        print("Ndf:", self.tunfolder.GetNdf())
        self.Ndf = self.tunfolder.GetNdf()
        print("Npar:", self.tunfolder.GetNpar())
        self.Npar = self.tunfolder.GetNpar()
        print("chi2sys:", self.tunfolder.GetChi2Sys())
        self.chi2sys = self.tunfolder.GetChi2Sys()
        print("chi2A:", self.tunfolder.GetChi2A())
        self.chi2A = self.tunfolder.GetChi2A()
        print("chi2L:", self.tunfolder.GetChi2L())
        self.chi2L = self.tunfolder.GetChi2L()

        # print("( " + str(self.tunfolder.GetChi2A()) + " + " + str(self.tunfolder.GetChi2L()) + ") / " + str(self.tunfolder.GetNdf()))

        self.unfolded = self.tunfolder.GetOutput(hist_name, "", self.output_distribution_name, "*[]", self.use_axis_binning)
        # if update_with_ematrix_total:
        #     self.update_unfolded_with_ematrix_total()
        return self.unfolded

    def _post_process(self):
        """Do some standard things & store various things that are done after unfolding"""
        self.get_ematrix_input()
        self.get_ematrix_sys_uncorr()
        self.get_ematrix_stat()
        self.get_ematrix_total()
        self.get_rhoij_total()
        self.get_probability_matrix()
        self.update_unfolded_with_ematrix_total()
        self.get_unfolded_with_ematrix_stat()
        self.get_folded_unfolded()
        self.get_folded_mc_truth()

    @staticmethod
    def make_hist_from_diagonal_errors(h2d, do_sqrt=True):
        nbins = h2d.GetNbinsX()
        hnew = ROOT.TH1D("h_diag" + cu.get_unique_str(), "", nbins, 0, nbins)
        for i in range(1, nbins+1):
            err = h2d.GetBinContent(i, i)
            if do_sqrt and err > 0:
                err = math.sqrt(err)
            hnew.SetBinError(i, err)
        return hnew

    @staticmethod
    def update_hist_bin_error(h_orig, h_to_be_updated):
        if h_orig.GetNbinsX() != h_to_be_updated.GetNbinsX():
            raise RuntimeError("Need same # x bins, %d vs %s" % (h_orig.GetNbinsX(), h_to_be_updated.GetNbinsX()))
        for i in range(0, h_orig.GetNbinsX()+2):
            h_to_be_updated.SetBinError(i, h_orig.GetBinError(i))

    def update_unfolded_with_ematrix_total(self):
        """Update unfolded hist with total errors from total error matrix"""
        ematrix_total = self.get_ematrix_total()
        error_total_1d = self.make_hist_from_diagonal_errors(ematrix_total, do_sqrt=True) # note that bin contents = 0, only bin errors are non-0
        self.update_hist_bin_error(h_orig=error_total_1d, h_to_be_updated=self.unfolded)

    def get_unfolded_with_ematrix_stat(self):
        """Make copy of unfolded, but only stat errors"""
        if getattr(self, 'unfolded_stat_err', None) is None:
            error_1d = self.make_hist_from_diagonal_errors(self.get_ematrix_stat(), do_sqrt=True) # note that bin contents = 0, only bin errors are non-0
            self.unfolded_stat_err = self.unfolded.Clone("unfolded_stat_err")
            self.update_hist_bin_error(h_orig=error_1d, h_to_be_updated=self.unfolded_stat_err)
        return self.unfolded_stat_err

    def get_bias_vector(self):
        if getattr(self, "bias_vector", None) is None:
            self.bias_vector = self.tunfolder.GetBias("bias_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)
        return self.bias_vector

    def get_ematrix_input(self):
        if getattr(self, "ematrix_input", None) is None:
            self.ematrix_input = self.tunfolder.GetEmatrixInput("ematrix_input_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)
        return self.ematrix_input

    def get_ematrix_sys_uncorr(self):
        if getattr(self, "ematrix_sys_uncorr", None) is None:
            self.ematrix_sys_uncorr = self.tunfolder.GetEmatrixSysUncorr("ematrix_sys_uncorr_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)
        return self.ematrix_sys_uncorr

    def get_ematrix_total(self):
        if getattr(self, "ematrix_total", None) is None:
            self.ematrix_total = self.tunfolder.GetEmatrixTotal("ematrix_total_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)
        return self.ematrix_total

    def get_ematrix_stat(self):
        if getattr(self, 'ematrix_stat', None) is None:
            self.ematrix_stat_sum = self.get_ematrix_input().Clone("ematrix_stat_sum")
            self.ematrix_stat_sum.Add(self.get_ematrix_sys_uncorr()) # total 'stat' errors
        return self.ematrix_stat_sum

    def get_rhoij_total(self):
        if getattr(self, "rhoij_total", None) is None:
            self.rhoij_total = self.tunfolder.GetRhoIJtotal("rhoij_total_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)
        return self.rhoij_total

    def get_probability_matrix(self):
        if getattr(self, "probability_matrix", None) is None:
            self.probability_matrix = self.tunfolder.GetProbabilityMatrix("prob_matrix_"+cu.get_unique_str(), "", self.use_axis_binning)
        return self.probability_matrix

    def get_covariance_matrix(self):
        if getattr(self, "covariance_matrix", None) is None:
            self.covariance_matrix = self.tunfolder.GetVxx()
        return self.covariance_matrix

    def get_var_hist_pt_binned(self, hist1d, ibin_pt, binning_scheme='generator'):
        """Get hist of variable for given pt bin from massive 1D hist that TUnfold makes"""
        # FIXME: assume no underflow?!
        binning = self.generator_binning.FindNode("generatordistribution") if binning_scheme == "generator" else self.detector_binning.FindNode("detectordistribution")
        var_bins = np.array(binning.GetDistributionBinning(0))
        pt_bins = np.array(binning.GetDistributionBinning(1))

        # print("var_bins:", var_bins)
        # print("pt_bins:", pt_bins)
        # bin_num = binning.GetGlobalBinNumber(0.001, 51)
        # print("Global bin num for (pt, lambda) = (51, 0.001) => %d" % (bin_num))
        # print("This bin goes from %g to %g" % (hist1d.GetXaxis().GetBinLowEdge(bin_num), hist1d.GetXaxis().GetBinLowEdge(bin_num+1)))

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

    def get_folded_unfolded(self):
        # don't use getfoldedoutput, because it doesn't have the updated errors from the total error matrix
        # so we'll have to do it ourselves
        # 1. Make unfolded hist into TVector/TMatrix

        # 2. Make response 2d hist into matrix

        # 3. Multiply the two, convert to TH1
        if getattr(self, 'folded_unfolded', None) is None:
            self.folded_unfolded = self.tunfolder.GetFoldedOutput("folded_unfolded")
        return self.folded_unfolded

    @staticmethod
    def th2_to_tmatrixd(hist, include_uflow=False, include_oflow=False):
        n_rows = hist.GetNbinsY()
        n_cols = hist.GetNbinsX()

        # ignore for now as too complicated
        # if include_uflow:
        #     n_rows += 1
        #     n_cols += 1
        # if include_oflow:
        #     n_rows += 1
        #     n_cols += 1

        # taken from https://root.cern.ch/doc/master/TH2_8cxx_source.html#l03739
        m = ROOT.TMatrixD(n_rows, n_cols)
        ilow = m.GetRowLwb()
        iup  = m.GetRowUpb()
        jlow = m.GetColLwb()
        jup  = m.GetColUpb()
        for i in range(ilow, iup+1):
            for j in range(jlow, jup+1):
                m[i,j] = hist.GetBinContent(j-jlow+1,i-ilow+1)
        return m

    @staticmethod
    def calculate_singular_max_min(matrix):
        """Calculate max & min singular values condition number as per StatsComm guidelines

        These are found using TDecompSVD.
        (we ignore the builtin condition() method as it calculates it differently)
        """
        if matrix.GetNcols() > matrix.GetNrows():
            raise RuntimeError("Condition number only for matrix where # rows >= # cols")

        svd = ROOT.TDecompSVD(matrix)
        sig = svd.GetSig()  # by construction, ordered descending
        sigma_max = sig[0]
        sigma_min = max(0, sig[sig.GetNrows()-1])
        print("sigma_max:", sigma_max, "sigma_min:", sigma_min)
        return sigma_max, sigma_min

    def print_condition_number(self):
        """Store & print response matrix condition number and some advice

        Defined as sigma_max / max(0, sigma_min), where sigma_{max/min} are the
        largest/smallest singular values.
        These are also stored for later usage if needed (since expensive to calc)
        """
        if getattr(self, 'condition_number', None) is None:
            sigma_max, sigma_min = self.calculate_singular_max_min(self.th2_to_tmatrixd(self.get_probability_matrix()))
            if sigma_min == 0:
                # avoid DivisionError
                print("Minmum singular value = 0, condition number = Infinity")
                num = np.inf
            else:
                num = sigma_max / sigma_min
            self.sigma_max = sigma_max
            self.sigma_min = sigma_min
            self.condition_number = num

        print("Condition number:", self.condition_number)
        if self.condition_number < 50:
            print(" - You probably shouldn't regularize this")
        elif self.condition_number > 1E5:
            print(" - You probably should regularize this")
        else:
            print(" - You probably should look into regularization")

    @staticmethod
    def th1_to_ndarray(hist_A, oflow_x=False):
        """Convert TH1 to numpy ndarray"""
        ncol = hist_A.GetNbinsX()
        if oflow_x:
            ncol += 2
        result = np.zeros(shape=(1, ncol), dtype=float)
        errors = np.zeros(shape=(1, ncol), dtype=float)

        # Get ROOT indices to loop over
        x_start = 0 if oflow_x else 1
        x_end = hist_A.GetNbinsX()
        if oflow_x:
            x_end += 1

        # x_ind for numpy as always starts at 0
        # ix for ROOT
        for x_ind, ix in enumerate(range(x_start, x_end+1)):
            result[0][x_ind] = hist_A.GetBinContent(ix)
            errors[0][x_ind] = hist_A.GetBinError(ix)

        # check sparsity
        return result, errors

    @staticmethod
    def ndarray_to_th1(nd_array, has_oflow_x=False):
        """Convert numpy ndarray row vector to TH1, with shape (1, nbins)

        Use has_oflow_x to include the under/overflow bins
        """
        nbinsx = nd_array.shape[1]
        nbins_hist = nbinsx
        if has_oflow_x:
            nbins_hist -= 2

        # need the 0.5 offset to match TUnfold
        h = ROOT.TH1F(cu.get_unique_str(), "", nbins_hist, 0.5, nbins_hist+0.5)

        x_start = 1
        x_end = nbins_hist

        if has_oflow_x:
            x_start = 0
            x_end = nbins_hist+1

        for x_ind, ix in enumerate(range(x_start, x_end+1)):
            h.SetBinContent(ix, nd_array[0][x_ind])
            h.SetBinError(ix, math.sqrt(nd_array[0][x_ind]))
            #FIXME how to do errors
        return h

    @staticmethod
    def th2_to_ndarray(hist_A, oflow_x=False, oflow_y=False):
        """Convert TH2 to numpy ndarray

        Don't use verison in common_utils - wrong axes?
        """
        ncol = hist_A.GetNbinsX()
        if oflow_x:
            ncol += 2
        nrow = hist_A.GetNbinsY()
        if oflow_y:
            nrow += 2

        result = np.zeros(shape=(nrow, ncol), dtype=float)
        errors = np.zeros(shape=(nrow, ncol), dtype=float)
        # access via result[irow][icol]

        # Get ROOT indices to loop over
        y_start = 0 if oflow_y else 1
        y_end = hist_A.GetNbinsY()
        if oflow_y:
            y_end += 1

        x_start = 0 if oflow_x else 1
        x_end = hist_A.GetNbinsX()
        if oflow_x:
            x_end += 1

        # y_ind, x_ind for numpy as always starts at 0
        # iy, ix for ROOT
        for y_ind, iy in enumerate(range(y_start, y_end+1)):
            for x_ind, ix in enumerate(range(x_start, x_end+1)):
                result[y_ind][x_ind] = hist_A.GetBinContent(ix, iy)
                errors[y_ind][x_ind] = hist_A.GetBinError(ix, iy)

        # check sparsity
        num_empty = np.count_nonzero(result == 0)
        num_entries = result.size
        sparsity = num_empty / float(num_entries)
        # print("Converting TH2 to ndarray...")
        # print("num_empty:", num_empty)
        # print("num_entries:", num_entries)
        # print("sparsity:", sparsity)
        if (sparsity > 0.5):
            print("Matrix has %d/%d empty entries - consider using sparse matrix (which I don't know how to do yet)" % (num_empty, num_entries))

        return result, errors

    @staticmethod
    def normalise_ndarray(matrix, by):
        if by == 'col':
            matrix = matrix.T # makes life a bit easier
        for i in range(matrix.shape[0]):
            row_sum = matrix[i].sum()
            if row_sum != 0:
                matrix[i] = matrix[i] / row_sum
        if by == 'col':
            return matrix.T
        else:
            return matrix

    def get_folded_hist(self, hist_gen):
        """Fold hist_gen using the stored respone matrix, ie do matrix * vector"""

        # TODO: proper error propagation

        oflow = True
        # Convert map to matrix
        response_matrix, response_matrix_err = self.th2_to_ndarray(self.response_map, oflow_x=oflow, oflow_y=oflow)

        # Normalise response_matrix so that bins represent prob to go from
        # given gen bin to a reco bin
        # ASSUMES GEN ON X AXIS!
        norm_by = 'col' if self.orientation == ROOT.TUnfold.kHistMapOutputHoriz else 'row'
        response_matrix_normed = self.normalise_ndarray(response_matrix, by='col')

        # Convert hist to vector
        gen_vec, gen_vec_err = self.th1_to_ndarray(hist_gen, oflow_x=oflow)

        # Multiply
        # Note that we need to transpose from row vec to column vec
        folded_vec = response_matrix_normed.dot(gen_vec.T)

        # Convert vector to TH1
        folded_hist = self.ndarray_to_th1(folded_vec.T, has_oflow_x=oflow)

        return folded_hist

    def get_folded_mc_truth(self):
        """Get response_matrix * MC truth"""
        if getattr(self, 'folded_mc_truth', None) is None:
            self.folded_mc_truth = self.get_folded_hist(self.hist_truth)
        return self.folded_mc_truth


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
                          xtitle='Generator bin', ytitle='Detector bin')

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
        output_filename = "%s/err_map_sys_input_%s.%s" % (output_dir, append, self.output_fmt)
        self.draw_2d_hist(self.unfolder.get_ematrix_input(),
                          title=title,
                          output_filename=output_filename,
                          draw_bin_lines_x=True,
                          draw_bin_lines_y=True,
                          canvas_size=(800, 700))

    def draw_error_matrix_sys_uncorr(self, output_dir='.', append="", title=""):
        output_filename = "%s/err_map_sys_uncorr_%s.%s" % (output_dir, append, self.output_fmt)
        self.draw_2d_hist(self.unfolder.get_ematrix_sys_uncorr(),
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

        hist_reco = self.unfolder.input_hist_bg_subtracted
        is_bg_subtracted = True
        if not hist_reco:
            hist_reco = self.input_hist
            is_bg_subtracted = False
        if hist_reco:
            entries.append(
                Contribution(hist_reco,
                             label="Unfolding input (background-subtracted)" if is_bg_subtracted else "Unfolding input",
                             line_color=ROOT.kBlue, line_width=1,
                             marker_color=ROOT.kBlue, marker_size=0,
                             normalise_hist=False)
                )

        hist_folded = self.unfolder.get_folded_unfolded()
        if hist_folded:
            entries.append(
                Contribution(hist_folded, label="Folded unfolded result (#tau = %.3g)" % (self.unfolder.tau),
                             line_color=ROOT.kRed, line_width=1,
                             marker_color=ROOT.kRed, marker_size=0,
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

        hist_reco = self.unfolder.input_hist_bg_subtracted
        is_bg_subtracted = True
        if not hist_reco:
            hist_reco = self.input_hist
            is_bg_subtracted = False
        if hist_reco:
            entries.append(
                Contribution(hist_reco,
                             label="Unfolding input (background-subtracted)" if is_bg_subtracted else "Unfolding input",
                             line_color=ROOT.kBlue, line_width=1,
                             marker_color=ROOT.kBlue, marker_size=0,
                             normalise_hist=False)
                )

        hist_folded = self.unfolder.get_folded_mc_truth()
        if hist_folded:
            entries.append(
                Contribution(hist_folded, label="Folded generator",
                             line_color=ROOT.kRed, line_width=1,
                             marker_color=ROOT.kRed, marker_size=0,
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
        output_filename = "%s/folded_gen_%s.%s" % (output_dir, append, self.output_fmt)
        plot.save(output_filename)
