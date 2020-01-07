"""
Classes to be used when determining regularisation in unfolding
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
