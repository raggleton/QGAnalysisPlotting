#!/usr/bin/env python

"""
Classes to be used in main unfolding script/routines
"""


from __future__ import print_function, division

from array import array
import numpy as np
import math

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

    def setInput(self, *args):
        self.tunfolder.SetInput(*args)

    def do_unfolding(self, tau):
        print(">>> Unfolding with tau =", tau)
        self.tunfolder.DoUnfold(tau)

    def get_output(self):
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

    def get_covariance_matrix(self):
        return self.tunfolder.GetVxx()

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
    def calculate_condition_num(matrix):
        """Calculate condition number as per StatsComm guidelines

        Defined as sigma_max / max(0, sigma_min), where sigma_{max/min} are the
        largest/smallest singular values.

        These are found using TDecompSVD.
        (we ignore the builtin condition() method as it calculates it differently)
        """
        if matrix.GetNcols() > matrix.GetNrows():
            raise RuntimeError("Condition number only for matrix where # rows >= # cols")

        svd = ROOT.TDecompSVD(matrix)
        sig = svd.GetSig()  # by construction, ordered descending
        sigma_max = sig[0]
        sigma_min = max(0, sig[sig.GetNrows()-1])
        if sigma_min == 0:
            print("Minmum singular value = 0, condition number = Infinity")
            return 999999999999999999999
        print("sigma_max:", sigma_max, "sigma_min:", sigma_min)
        return sigma_max / sigma_min

    def print_condition_number(self):
        """Print response matrix condition number and some advice"""
        num = self.calculate_condition_num(self.response_map_matrix)
        print("Condition number:", num)
        if num < 50:
            print(" - You probably shouldn't regularize this")
        elif num > 1E5:
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

