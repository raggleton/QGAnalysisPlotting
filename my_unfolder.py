"""
Main class to handle unfolding
"""


from __future__ import print_function, division

from array import array
import numpy as np
import math
import os
from itertools import chain
import scipy
from scipy import stats

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

# Load my derived class
with open("MyTUnfoldDensity.cpp") as f:
    code = f.read()
    ROOT.gInterpreter.ProcessLine(code)


class MyUnfolder(ROOT.MyTUnfoldDensity):
    """Main class to handle unfolding input/outputs, all the associated objects

    Derived from MyTUnfoldDensity to get access to protected methods/vars
    """
    # for saving/loading
    _simple_attr = [
        # input & truth histograms
        "input_hist",
        "input_hist_bg_subtracted",
        "input_hist_gen_binning",
        "input_hist_gen_binning_bg_subtracted",
        "hist_truth",
        "hist_mc_reco",
        "hist_mc_reco_bg_subtracted",
        "hist_mc_reco_gen_binning",
        "hist_mc_reco_gen_binning_bg_subtracted",
        # save other matrices
        "rhoij_total",
        "probability_matrix",
        # save error matrices
        "ematrix_input",
        "ematrix_stat_response",
        "ematrix_stat",
        "ematrix_tau",
        "ematrix_total",
        # save folded things
        "folded_unfolded",
        "folded_mc_truth",
        # unfolded
        "unfolded",
        "unfolded_stat_err",
    ]

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
        # check no uflow
        if orientation == ROOT.TUnfold.kHistMapOutputHoriz:
            if response_map.GetBinContent(0, 0) != 0 or response_map.GetBinContent(0, 1) != 0:
                raise RuntimeError("Your response_map has entries in 0th gen bin - this means you've got unintended underflow!")
        elif orientation == ROOT.TUnfold.kHistMapOutputVert:
            if response_map.GetBinContent(0, 1) != 0 or response_map.GetBinContent(1, 0) != 0:
                raise RuntimeError("Your response_map has entries in 0th gen bin - this means you've got unintended underflow!")
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
        var_uf, var_of = False, True
        pt_uf, pt_of = False, True  # handle pt under/over flow ourselves
        self.detector_binning = ROOT.TUnfoldBinning("detector")

        self.detector_distribution_underflow = self.detector_binning.AddBinning("detectordistribution_underflow")
        if self.variable_bin_edges_reco is not None:
            self.detector_distribution_underflow.AddAxis(self.variable_name, self.nbins_variable_reco, self.variable_bin_edges_reco, var_uf, var_of)
        self.detector_distribution_underflow.AddAxis("pt", self.nbins_pt_underflow_reco, self.pt_bin_edges_underflow_reco, False, False)

        self.detector_distribution = self.detector_binning.AddBinning("detectordistribution")
        if self.variable_bin_edges_reco is not None:
            self.detector_distribution.AddAxis(self.variable_name, self.nbins_variable_reco, self.variable_bin_edges_reco, var_uf, var_of)
        self.detector_distribution.AddAxis("pt", self.nbins_pt_reco, self.pt_bin_edges_reco, False, pt_of)


        self.generator_binning = ROOT.TUnfoldBinning("generator")

        self.generator_distribution_underflow = self.generator_binning.AddBinning("generatordistribution_underflow")
        if self.variable_bin_edges_gen is not None:
            self.generator_distribution_underflow.AddAxis(self.variable_name, self.nbins_variable_gen, self.variable_bin_edges_gen, var_uf, var_of)
        self.generator_distribution_underflow.AddAxis("pt", self.nbins_pt_underflow_gen, self.pt_bin_edges_underflow_gen, pt_uf, False)

        self.generator_distribution = self.generator_binning.AddBinning("generatordistribution")
        if self.variable_bin_edges_gen is not None:
            self.generator_distribution.AddAxis(self.variable_name, self.nbins_variable_gen, self.variable_bin_edges_gen, var_uf, var_of)
        self.generator_distribution.AddAxis("pt", self.nbins_pt_gen, self.pt_bin_edges_gen, False, pt_of)

        self.orientation = orientation
        self.constraintMode = constraintMode
        self.regMode = regMode
        self.densityFlags = densityFlags
        self.distribution = distribution
        self.axisSteering = axisSteering

        ROOT.TUnfoldDensity.__init__(self,
                                     self.response_map,
                                     self.orientation,
                                     self.regMode,
                                     self.constraintMode,
                                     self.densityFlags,
                                     self.generator_binning,
                                     self.detector_binning,
                                     # hmm these take preference over whatever is use for scantau?
                                     self.distribution,
                                     self.axisSteering)

        self.use_axis_binning = False  # for things like get_probability_matrix()...but doesn't seem to do anything?!

        # self.probability_ndarray = self.response_matrix_to_probability_array(self.response_map)
        # self.probability_ndarray, _ = self.th2_to_ndarray(self.get_probability_matrix(), oflow_x=False, oflow_y=False)
        # self.probability_ndarray, _ = None, None

        # hists that will be assigned later
        # TODO: change to properties? although still need to cache somewhere

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
        self.syst_shifts = {}  # gets filled with get_sys_shift(), just shift in unfolded value from specific syst
        self.systs_shifted = {}  # gets filled with get_syst_shifted_hist(), holds total unfolded with syst shift
        self.syst_ematrices = {}  # gets filled with get_ematrix_syst(), holds ematrix for each systeamtic

        # use "generator" for signal + underflow region, "generatordistribution" for only signal region
        self.output_distribution_name = "generator"

        self.folded_unfolded = None  # set in get_folded_unfolded()
        self.folded_unfolded_tunfold = None  # set in get_folded_unfolded()
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
        if obj is not None:  # be careful of int = 0!
            # for "simple" types, need to do some trickery
            if isinstance(obj, type(np.array([1.]))):
                if len(obj.shape) > 1:
                    raise RuntimeError("cannot save numpy unless 1D")
                n = ROOT.TVectorD(len(obj))
                for ind, val in enumerate(obj):
                    n[ind] = val
                tfile.WriteTObject(n, name)
            elif isinstance(obj, (bool, int, float)):
                n = ROOT.TVectorD(1)  # the only class that works for these types
                n[0] = obj
                tfile.WriteTObject(n, name)
            elif isinstance(obj, str):
                n = ROOT.TNamed(name, obj)
                tfile.WriteTObject(n, name)
            else:
                if name is None:
                    name = obj.GetName()
                tfile.WriteTObject(obj, name)
        else:
            print("Not saving", name, "as does not exist")

    def save_to_tfile(self, tfile):
        """Save important stuff to TFile/TDirectory"""
        self._check_save_to_tfile(tfile, self.detector_binning, "detector_binning")
        self._check_save_to_tfile(tfile, self.generator_binning, "generator_binning")
        self._check_save_to_tfile(tfile, self.response_map, "response_map")
        self._check_save_to_tfile(tfile, self.orientation, "orientation")
        self._check_save_to_tfile(tfile, self.constraintMode, "constraintMode")
        self._check_save_to_tfile(tfile, self.regMode, "regMode")
        self._check_save_to_tfile(tfile, self.densityFlags, "densityFlags")
        self._check_save_to_tfile(tfile, self.distribution, "distribution")
        self._check_save_to_tfile(tfile, self.axisSteering, "axisSteering")

        self._check_save_to_tfile(tfile, self.variable_bin_edges_reco, "variable_bin_edges_reco")
        self._check_save_to_tfile(tfile, self.variable_bin_edges_gen, "variable_bin_edges_gen")
        self._check_save_to_tfile(tfile, self.variable_name, "variable_name")

        self._check_save_to_tfile(tfile, self.pt_bin_edges_reco, "pt_bin_edges_reco")
        self._check_save_to_tfile(tfile, self.pt_bin_edges_gen, "pt_bin_edges_gen")

        self._check_save_to_tfile(tfile, self.pt_bin_edges_underflow_reco, "pt_bin_edges_underflow_reco")
        self._check_save_to_tfile(tfile, self.pt_bin_edges_underflow_gen, "pt_bin_edges_underflow_gen")

        # handle most of the simple hists
        for name in MyUnfolder._simple_attr:
            self._check_save_to_tfile(tfile, getattr(self, name), name)

        # save all backgrounds (incl fakes)
        for name, hist in self.backgrounds.items():
            self._check_save_to_tfile(tfile, hist, "background_reco_binning_%s" % name.replace(" ", "_"))
        for name, hist in self.backgrounds_gen_binning.items():
            self._check_save_to_tfile(tfile, hist, "background_gen_binning_%s" % name.replace(" ", "_"))

        # save systematic response matrices
        for name, syst_map in self.syst_maps.items():
            self._check_save_to_tfile(tfile, syst_map, "syst_map_%s" % name.replace(" ", "_"))

        # save systematic error matrices
        for name, syst_ematrix in self.syst_ematrices.items():
            self._check_save_to_tfile(tfile, syst_ematrix, "syst_ematrix_%s" % name.replace(" ", "_"))

        # save systematic shifts
        for name, syst_shift in self.syst_shifts.items():
            self._check_save_to_tfile(tfile, syst_shift, "syst_shift_%s" % name.replace(" ", "_"))

        # save systematic shifted hists (yes this is a bit wasteful)
        for name, syst_shift in self.systs_shifted.items():
            self._check_save_to_tfile(tfile, syst_shift, "syst_shifted_unfolded_%s" % name.replace(" ", "_"))

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
        self.input_hist_bg_subtracted = input_hist.Clone() if input_hist else None
        self.input_hist_gen_binning = input_hist_gen_binning.Clone() if input_hist_gen_binning else None
        self.input_hist_gen_binning_bg_subtracted = input_hist_gen_binning.Clone() if input_hist_gen_binning else None
        self.hist_truth = hist_truth
        self.hist_mc_reco = hist_mc_reco
        self.hist_mc_reco_bg_subtracted = hist_mc_reco_bg_subtracted
        self.hist_mc_reco_gen_binning = hist_mc_reco_gen_binning
        self.hist_mc_reco_gen_binning_bg_subtracted = hist_mc_reco_gen_binning_bg_subtracted
        self.SetInput(input_hist, bias_factor)

    def subtract_background(self, hist, name, scale=1.0, scale_err=0.0):
        """Subtract background source from input hist"""
        # Save into dict of components - needed? since TUnfoldDensity does this as well
        self.backgrounds[name] = hist.Clone()
        self.backgrounds[name].Scale(scale)
        # Also save total input subtracted
        self.input_hist_bg_subtracted.Add(hist, -1*scale)
        self.SubtractBackground(hist.Clone(), name, scale, scale_err)

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
        self.DoUnfold(tau)

    def calculate_pt_bin_factors(self, which):
        """Calculate bin factors to account for falling distributions when regularising

        NB done according to signal region - excludes underflow!

        which : str
            Choose which histogram to use for integrals, 'gen' or 'unfolded'
        """

        # Tricky - need to get counts in spectrum, ideally with gen binning
        # FIXME use input_hist instead with detector binning?
        which = which.lower()
        if which not in ['unfolded', 'gen']:
            raise ArgumentError("calculate_pt_bin_factors: 'which' arg should be 'unfolded' or 'gen'")

        hist = None
        if which == 'gen':
            if self.input_hist_gen_binning_bg_subtracted is None:
                raise RuntimeError("Need input_hist_gen_binning_bg_subtracted to be able to do calculate_pt_bin_factors")
            hist = self.input_hist_gen_binning_bg_subtracted

        if which == 'unfolded':
            if self.unfolded is None:
                raise RuntimeError("Need unfolded to be able to do calculate_pt_bin_factors")
            hist = self.unfolded

        # Get integral of 1st pt bin in signal region
        gen_node = self.generator_binning.FindNode('generatordistribution')
        first_var = self.variable_bin_edges_gen[0]
        last_var = self.variable_bin_edges_gen[-1]
        pt_val = self.pt_bin_edges_gen[0]
        start_bin = gen_node.GetGlobalBinNumber(first_var+0.00001, pt_val+0.001)
        end_bin = gen_node.GetGlobalBinNumber(last_var-0.00001, pt_val+0.001)
        first_bin_integral = hist.Integral(start_bin, end_bin)  # ROOTs integral is inclusive of last bin

        bin_factors = {}
        # Add 1s for the 1st pt bin
        for var in self.variable_bin_edges_gen[:-1]:
            this_bin = gen_node.GetGlobalBinNumber(var+0.00001, pt_val+0.001)
            bin_factors[this_bin] = 1

        # Iterate through pt bins, figure out integral, scale according to first bin
        for pt_val in self.pt_bin_edges_gen[1:-1]:
            start_bin = gen_node.GetGlobalBinNumber(first_var+0.00001, pt_val+0.001)
            end_bin = gen_node.GetGlobalBinNumber(last_var-0.00001, pt_val+0.001)
            integral = hist.Integral(start_bin, end_bin)
            sf = first_bin_integral / integral

            # Store bin factor for each lambda bin
            for var in self.variable_bin_edges_gen[:-1]:
                this_bin = gen_node.GetGlobalBinNumber(var+0.00001, pt_val+0.001)
                bin_factors[this_bin] = sf

        return bin_factors

    def get_gen_bin_widths(self):
        results = {}
        # loop through the real axes, convert to global bin number, store width
        # This is because TUnfoldBinning has no *public* method to go the other
        # way...ToAxisBins() is *protected* FFS
        for lambda_ind, lambda_var in enumerate(self.variable_bin_edges_gen[:-1]):
            # underflow region
            for pt_ind, pt_val in enumerate(self.pt_bin_edges_underflow_gen[:-1]):
                global_bin = self.generator_distribution_underflow.GetGlobalBinNumber(lambda_var+0.00000001, pt_val+0.0000001)
                lambda_width = self.variable_bin_edges_gen[lambda_ind+1] - self.variable_bin_edges_gen[lambda_ind]
                pt_width = self.pt_bin_edges_underflow_gen[pt_ind+1] - self.pt_bin_edges_underflow_gen[pt_ind]
                results[global_bin] = (lambda_width, pt_width)
            # signal region
            for pt_ind, pt_val in enumerate(self.pt_bin_edges_gen[:-1]):
                global_bin = self.generator_distribution.GetGlobalBinNumber(lambda_var+0.00000001, pt_val+0.0000001)
                lambda_width = self.variable_bin_edges_gen[lambda_ind+1] - self.variable_bin_edges_gen[lambda_ind]
                pt_width = self.pt_bin_edges_gen[pt_ind+1] - self.pt_bin_edges_gen[pt_ind]
                results[global_bin] = (lambda_width, pt_width)
        return results

    def add_sys_error(self, map_syst, name, mode):
        """Add systematic error via response map, arguments as per AddSysError()"""
        self.syst_maps[name] = map_syst
        self.AddSysError(map_syst, name, self.orientation, ROOT.TUnfoldDensity.kSysErrModeMatrix)
        self.syst_shifts[name] = None  # setup for get_delta_sys_shift
        self.systs_shifted[name] = None  # setup for get_syst_shifted_hist
        self.syst_ematrices[name] = None  # setup for get_ematrix_syst

    def get_delta_sys_shift(self, syst_label):
        """Get shift in result due to a particular systeamtic

        Label must be same as used to add it in add_sys_error()
        """
        if syst_label not in self.syst_shifts:
            raise KeyError("No systematic %s, only have: %s" % (syst_label, ", ".join(self.syst_shifts.keys())))
        if self.syst_shifts[syst_label] is None:
            hist = self.GetDeltaSysSource(syst_label,
                                                    "syst_shift_%s" % (syst_label.replace(" ", "_")),
                                                    "",
                                                    self.output_distribution_name, # must be the same as what's used in get_output
                                                    self.axisSteering,
                                                    self.use_axis_binning)
            self.syst_shifts[syst_label] = hist  # cache shifts
        return self.syst_shifts[syst_label]

    def get_syst_shifted_hist(self, syst_label, unfolded=None):
        """Get unfolded hist, shifted by a given systematic source

        Can specify unfolded hist, default is the one with all errors
        """
        if syst_label not in self.syst_shifts:
            raise KeyError("No systematic %s, only have: %s" % (syst_label, ", ".join(self.syst_shifts.keys())))
        if self.systs_shifted[syst_label] is None:
            hist_shift = self.get_delta_sys_shift(syst_label).Clone('syst_shifted_unfolded_%s' % syst_label.replace(" ", "_"))
            unfolded = unfolded or self.unfolded
            hist_shift.Add(unfolded)  # TODO what about errors?
            self.systs_shifted[syst_label] = hist_shift
        return self.systs_shifted[syst_label]

    def get_output(self, hist_name='unfolded'):
        """Get 1D unfolded histogram covering all bins"""
        print("Ndf:", self.GetNdf())
        self.Ndf = self.GetNdf()
        print("Npar:", self.GetNpar())
        self.Npar = self.GetNpar()
        print("chi2sys:", self.GetChi2Sys())
        self.chi2sys = self.GetChi2Sys()
        print("chi2A:", self.GetChi2A())
        self.chi2A = self.GetChi2A()
        print("chi2L:", self.GetChi2L())
        self.chi2L = self.GetChi2L()

        self.unfolded = self.GetOutput(hist_name, "", self.output_distribution_name, "*[]", self.use_axis_binning)
        return self.unfolded

    def _post_process(self):
        """Do some standard things & store various things that are done after unfolding"""
        self.get_ematrix_input()
        self.get_ematrix_stat_response()
        self.get_ematrix_stat()
        self.get_ematrix_tau()
        self.get_ematrix_total()
        self.get_rhoij_total()
        self.get_probability_matrix()
        self.update_unfolded_with_ematrix_total()
        self.get_unfolded_with_ematrix_stat()
        self.get_folded_unfolded()
        self.get_folded_mc_truth()
        for syst_label in self.syst_shifts.keys():
            # setup all the internal maps
            self.get_delta_sys_shift(syst_label)
            self.get_ematrix_syst(syst_label)
            self.get_syst_shifted_hist(syst_label)

    @staticmethod
    def make_hist_from_diagonal_errors(h2d, do_sqrt=True):
        """Make 1D hist, with errors set to diagonal elements from h2d

        Can be TH2 or numpy.ndarray, cos we have to use both
        Yes that is majorly wack
        """
        if isinstance(h2d, ROOT.TH2):
            nbins = h2d.GetNbinsX()
            hnew = ROOT.TH1D("h_diag" + cu.get_unique_str(), "", nbins, 0, nbins)
            for i in range(1, nbins+1):
                err = h2d.GetBinContent(i, i)
                if do_sqrt and err > 0:
                    err = math.sqrt(err)
                hnew.SetBinError(i, err)
            return hnew
        elif isinstance(h2d, np.ndarray):
            nbins = h2d.shape[0]
            hnew = ROOT.TH1D("h_diag" + cu.get_unique_str(), "", nbins, 0, nbins)
            for i in range(1, nbins+1):
                err = h2d[i-1, i-1]
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
        error_total_1d = self.make_hist_from_diagonal_errors(self.get_ematrix_total(), do_sqrt=True) # note that bin contents = 0, only bin errors are non-0
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
            self.bias_vector = self.GetBias("bias_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)
        return self.bias_vector

    def get_ematrix_input(self):
        """Get error matrix due to statistics from thing being unfolded"""
        if getattr(self, "ematrix_input", None) is None:
            self.ematrix_input = self.GetEmatrixInput("ematrix_input_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)
        return self.ematrix_input

    def get_ematrix_stat_response(self):
        """Statistical uncertainty error matrix from response matrix, should be considered a systematic uncert"""
        if getattr(self, "ematrix_stat_response", None) is None:
            self.ematrix_stat_response = self.GetEmatrixSysUncorr("ematrix_stat_response_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)
        return self.ematrix_stat_response

    def get_ematrix_total(self):
        """Total error matrix, from stat+systs"""
        if getattr(self, "ematrix_total", None) is None:
            self.ematrix_total = self.GetEmatrixTotal("ematrix_total_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)
        return self.ematrix_total

    def get_ematrix_stat(self):
        """Get total statitical error matrix (from input being unfolded + background sources, including fakes)"""
        if getattr(self, 'ematrix_stat', None) is None:
            # Have to manually create hist first, awkward
            this_binning = self.generator_binning.FindNode('generator')
            # I cannot figure out how to make the int** object for bin_map
            # So we are trusting that the default args for title and axisSteering are correct
            # Gnahhhhhhh
            self.ematrix_stat = this_binning.CreateErrorMatrixHistogram("ematrix_stat_"+cu.get_unique_str(), self.use_axis_binning) #, bin_map, "", "*[]")
            self.GetEmatrix(self.ematrix_stat)
        return self.ematrix_stat

    def get_ematrix_syst(self, syst_label):
        """Get error matrix from a systematic source"""
        if syst_label not in self.syst_shifts:
            raise KeyError("No systematic %s, only have: %s" % (syst_label, ", ".join(self.syst_shifts.keys())))
        if self.syst_ematrices[syst_label] is None:
            # Have to manually create hist first, awkward
            this_binning = self.generator_binning.FindNode('generator')
            # I cannot figure out how to make the int** object for bin_map
            # So we are trusting that the default args for title and axisSteering are correct
            # Gnahhhhhhh
            syst_label_no_spaces = syst_label.replace(" ", "_")
            hist = this_binning.CreateErrorMatrixHistogram("ematrix_syst_%s_%s" % (syst_label_no_spaces, cu.get_unique_str()), self.use_axis_binning) #, bin_map, "", "*[]")
            self.GetEmatrixSysSource(hist, syst_label)
            self.syst_ematrices[syst_label] = hist
        return self.syst_ematrices[syst_label]

    def get_ematrix_tau(self):
        """Get error matrix due to regularisation uncertainty"""
        if getattr(self, 'ematrix_tau', None) is None:
            # Have to manually create hist first, awkward
            this_binning = self.generator_binning.FindNode('generator')
            # I cannot figure out how to make the int** object for bin_map
            # So we are trusting that the default args for title and axisSteering are correct
            # Gnahhhhhhh
            self.ematrix_tau = this_binning.CreateErrorMatrixHistogram("ematrix_tau_"+cu.get_unique_str(), self.use_axis_binning) #, bin_map, "", "*[]")
            self.GetEmatrixSysTau(self.ematrix_tau)
        return self.ematrix_tau

    def get_ematrix_total_inv(self):
        """Total error matrix inverted, from stat+systs"""
        return self.InvertMSparseSymmPos(self.GetSummedErrorMatrixXX(), False)  # TMatrixDSparse

    def get_rhoij_total(self):
        if getattr(self, "rhoij_total", None) is None:
            self.rhoij_total = self.GetRhoIJtotal("rhoij_total_"+cu.get_unique_str(), "", "generator", "*[]", self.use_axis_binning)
        return self.rhoij_total

    def get_probability_matrix(self):
        if getattr(self, "probability_matrix", None) is None:
            self.probability_matrix = self.GetProbabilityMatrix("prob_matrix_"+cu.get_unique_str(), "", self.use_axis_binning)
        return self.probability_matrix

    @property
    def probability_ndarray(self):
        cached_attr_name = '_probability_ndarray'
        if not hasattr(self, cached_attr_name):
            arr, _ = self.th2_to_ndarray(self.get_probability_matrix())
            setattr(self, cached_attr_name, arr)
        return getattr(self, cached_attr_name)

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

    def get_pt_hist_var_binned(self, hist1d, ibin_var, binning_scheme='generator'):
        """Get hist of pt for given variable bin from massive 1D hist that TUnfold makes"""
        # FIXME: assume no underflow?!
        binning = self.generator_binning.FindNode("generatordistribution") if binning_scheme == "generator" else self.detector_binning.FindNode("detectordistribution")
        var_bins = np.array(binning.GetDistributionBinning(0))
        pt_bins = np.array(binning.GetDistributionBinning(1))

        # need the -1 on ibin_var, as it references an array index, whereas ROOT bins start at 1
        h = ROOT.TH1D("h_%d_%s" % (ibin_var, cu.get_unique_str()), "", len(pt_bins)-1, pt_bins)
        for pt_ind, pt_value in enumerate(pt_bins[:-1], 1):
            this_val = pt_value * 1.001  # ensure its inside
            bin_num = binning.GetGlobalBinNumber(var_bins[ibin_var]*1.001, this_val)
            h.SetBinContent(pt_ind, hist1d.GetBinContent(bin_num))
            h.SetBinError(pt_ind, hist1d.GetBinError(bin_num))
        return h

    # @staticmethod
    def tmatrixdsparse_to_ndarray(self, matrix):
        ndarr = np.zeros(shape=(matrix.GetNrows(), matrix.GetNcols()))

        rows_A = matrix.GetRowIndexArray()
        cols_A = matrix.GetColIndexArray()
        data_A = matrix.GetMatrixArray()
        for iy in range(matrix.GetNrows()):
            for indexA in range(rows_A[iy], rows_A[iy+1]):
                ix = cols_A[indexA]
                # print([x for x in self.GetXToHist()])
                # TODO: care about orientation?
                ndarr[iy, ix] = data_A[indexA]
        return ndarr

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
        result = np.zeros(shape=(1, ncol), dtype=np.float64)
        errors = np.zeros(shape=(1, ncol), dtype=np.float64)

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
            h.SetBinError(ix, math.sqrt(abs(nd_array[0][x_ind])))
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

        result = np.zeros(shape=(nrow, ncol), dtype=np.float64)
        errors = np.zeros(shape=(nrow, ncol), dtype=np.float64)
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
        # num_empty = np.count_nonzero(result == 0)
        # num_entries = result.size
        # sparsity = num_empty / float(num_entries)
        # print("Converting TH2 to ndarray...")
        # print("num_empty:", num_empty)
        # print("num_entries:", num_entries)
        # print("sparsity:", sparsity)
        # if (sparsity > 0.5):
        #     print("Matrix has %d/%d empty entries - consider using sparse matrix (which I don't know how to do yet)" % (num_empty, num_entries))

        return result, errors

    @staticmethod
    def ndarray_to_th2(data):
        nbinsy, nbinsx = data.shape
        binsx = array('d', list(range(1, nbinsx+2)))
        binsy = array('d', list(range(1, nbinsy+2)))
        h = ROOT.TH2D(cu.get_unique_str(), "", nbinsx, binsx, nbinsy, binsy)
        for ix in range(nbinsx):
            for iy in range(nbinsy):
                h.SetBinContent(ix+1, iy+1, data[iy,ix])
        return h

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

    def response_matrix_to_probability_array(self, response_map):
        """Convert response map to probability matrix

        Non-trivial, since to normalise one must also take into account underflow bins
        """
        prob_array, prob_err = self.th2_to_ndarray(response_map, oflow_x=False, oflow_y=False)
        # need one with oflow bins to normalise properly
        prob_array_oflow, prob_err = self.th2_to_ndarray(response_map, oflow_x=False, oflow_y=True)
        axis = 0 if self.orientation == ROOT.TUnfold.kHistMapOutputHoriz else 1
        sum_over_y = prob_array_oflow.sum(axis=axis)
        sum_over_y[sum_over_y==0] = 9999999  # avoid division errors
        result = prob_array / sum_over_y

        return result

    def get_folded_unfolded(self):
        # don't use getfoldedoutput, because it doesn't have the updated errors from the total error matrix
        # so we'll have to do it ourselves
        # 1. Make unfolded hist into TVector/TMatrix

        # 2. Make response 2d hist into matrix

        # 3. Multiply the two, convert to TH1

        # Get the TUnfold one for reference, although its errors will be wrong
        if getattr(self, 'folded_unfolded_tunfold', None) is None:
            self.folded_unfolded_tunfold = self.GetFoldedOutput("folded_unfolded_tunf")

        if getattr(self, 'folded_unfolded', None) is None:
            oflow = False

            # Get unfolded results as array
            unfolded_vector, _ = self.th1_to_ndarray(self.unfolded, oflow)

            # Multiply
            # Note that we need to transpose from row vec to column vec
            folded_vec = self.probability_ndarray.dot(unfolded_vector.T)

            # Convert vector to TH1
            self.folded_unfolded = self.ndarray_to_th1(folded_vec.T, has_oflow_x=oflow)

            # Error propagation: if y = Ax, with covariance matrices Vyy and Vxx,
            # respectively, then Vyy = (A*Vxx)*A^T
            unfolded_covariance_matrix, _ = self.th2_to_ndarray((self.get_ematrix_total()), oflow_x=oflow, oflow_y=oflow)
            result = self.probability_ndarray.dot(unfolded_covariance_matrix)
            folded_covariance = result.dot(self.probability_ndarray.T)
            folded_errors = self.make_hist_from_diagonal_errors(folded_covariance)
            self.update_hist_bin_error(h_orig=folded_errors, h_to_be_updated=self.folded_unfolded)

        return self.folded_unfolded

    def get_folded_hist(self, hist_gen):
        """Fold hist_gen using the stored response matrix, ie do matrix * vector"""
        oflow = False
        # Convert hist to vector
        gen_vec, gen_vec_err = self.th1_to_ndarray(hist_gen, oflow_x=oflow)

        # Multiply
        # Note that we need to transpose from row vec to column vec
        folded_vec = self.probability_ndarray.dot(gen_vec.T)

        # Convert vector to TH1
        folded_hist = self.ndarray_to_th1(folded_vec.T, has_oflow_x=oflow)

        # Error propagation: if y = Ax, with covariance matrices Vyy and Vxx,
        # respectively, then Vyy = (A*Vxx)*A^T
        result = self.probability_ndarray.dot(self.construct_covariance_matrix(hist_gen))
        folded_covariance = result.dot(self.probability_ndarray.T)
        folded_errors = self.make_hist_from_diagonal_errors(folded_covariance)
        self.update_hist_bin_error(h_orig=folded_errors, h_to_be_updated=folded_hist)

        return folded_hist

    def get_folded_mc_truth(self):
        """Get response_matrix * MC truth"""
        if getattr(self, 'folded_mc_truth', None) is None:
            self.folded_mc_truth = self.get_folded_hist(self.hist_truth)
        return self.folded_mc_truth

    def construct_covariance_matrix(self, hist):
        """Make covariance hist from 1D hist. Diagonals are square of error"""
        # Just to numpy array for now, can always use ndarray_to_th2()
        nbins = hist.GetNbinsX()
        cov_matrix_ndarray = np.zeros(shape=(nbins, nbins))  # important, otherwise can end up with wacko values
        for ind in range(1, nbins):
            err = hist.GetBinError(ind)
            err2 = pow(err, 2)
            cov_matrix_ndarray[ind-1, ind-1] = err2
        return cov_matrix_ndarray

    def calculate_smeared_chi2(self, ignore_underflow_bins=True):
        # TODO: should this all be divided by bin width?
        folded_vec, _ = self.th1_to_ndarray(self.get_folded_mc_truth(), False)
        reco_bg_subtracted_vec, _ = self.th1_to_ndarray(self.hist_mc_reco_bg_subtracted, False)
        delta = reco_bg_subtracted_vec - folded_vec
        # print(delta.shape)
        # print(delta)
        if ignore_underflow_bins:
            # print('start bin', self.detector_distribution.GetStartBin())
            first_signal_bin = self.detector_distribution.GetStartBin()
            # print('hist_mc_reco_bg_subtracted', self.hist_mc_reco_bg_subtracted.GetBinContent(first_signal_bin))
            # print('reco_bg_subtracted_vec', reco_bg_subtracted_vec[0, first_signal_bin-1])
            # print('delta', delta[0, first_signal_bin-1])
            delta[0][:first_signal_bin-1] = 0. # subtract 1 as numpy indices start at 0, hists start at 1
            # print(delta[0, 0:first_signal_bin-1])
            # print(delta[0, first_signal_bin-1:first_signal_bin+5])
            # print(reco_bg_subtracted_vec[0,first_signal_bin-5: first_signal_bin+5])
        vyy_inv = self.tmatrixdsparse_to_ndarray(self.GetVyyInv())
        print('vyy_inv', vyy_inv)
        inter = vyy_inv.dot(delta.T)
        # print(inter.shape)
        chi2 = delta.dot(inter)[0][0]
        ndof = len(delta[0][first_signal_bin-1:])
        p = 1-scipy.stats.chi2.cdf(chi2, int(ndof))
        return chi2, ndof, p

    def calculate_unfolded_chi2(self, ignore_underflow_bins=True):
        """(unfolded - truth)^T V_unfolded^{-1} (unfolded - truth)"""
        oflow = False
        unfolded_vec, _ = self.th1_to_ndarray(self.unfolded, oflow)
        gen_vec, _ = self.th1_to_ndarray(self.hist_truth, oflow_x=oflow)
        delta = unfolded_vec - gen_vec
        if ignore_underflow_bins:
            first_signal_bin = self.generator_distribution.GetStartBin()
            print('first_signal_bin', first_signal_bin)
            print(self.generator_distribution.GetEndBin())
            print('truth', self.hist_truth.GetBinContent(first_signal_bin))
            print('truth', gen_vec[0][first_signal_bin-1:])
            delta[0][:first_signal_bin-1] = 0

        # Have to manually create hist first, awkward
        this_binning = self.generator_binning.FindNode('generator')
        # I cannot figure out how to make the int** object for bin_map
        # So we are trusting that the default args for title and axisSteering are correct
        # Gnahhhhhhh
        ematrix_inv = this_binning.CreateErrorMatrixHistogram("ematrix_vxxinv_"+cu.get_unique_str(), self.use_axis_binning) #, bin_map, "", "*[]")
        self.ErrorMatrixToHist(ematrix_inv, self.GetVxxInv())
        # cov = self.tmatrixdsparse_to_ndarray(self.GetVxxInv())
        cov, _ = self.th2_to_ndarray(ematrix_inv)
        inter = cov.dot(delta.T)
        chi2 = delta.dot(inter)[0][0]
        ndof = len(delta[0][first_signal_bin-1:])
        p = 1-scipy.stats.chi2.cdf(chi2, int(ndof))
        return chi2, ndof, p



def unfolder_from_tdir(tdir):
    """Recover MyUnfolder from tdirectory

    Massive pain, but can't figure out how to serialize via pickle or ROOT
    """
    unfolder = MyUnfolder(response_map=cu.get_from_tfile(tdir, "response_map"),
                          variable_bin_edges_reco=np.array(cu.get_from_tfile(tdir, "variable_bin_edges_reco")),
                          variable_bin_edges_gen=np.array(cu.get_from_tfile(tdir, "variable_bin_edges_gen")),
                          variable_name=str(cu.get_from_tfile(tdir, "variable_name").GetTitle()),
                          pt_bin_edges_reco=np.array(cu.get_from_tfile(tdir, "pt_bin_edges_reco")),
                          pt_bin_edges_gen=np.array(cu.get_from_tfile(tdir, "pt_bin_edges_gen")),
                          pt_bin_edges_underflow_reco=np.array(cu.get_from_tfile(tdir, "pt_bin_edges_underflow_reco")),
                          pt_bin_edges_underflow_gen=np.array(cu.get_from_tfile(tdir, "pt_bin_edges_underflow_gen")),
                          orientation=int(cu.get_from_tfile(tdir, "orientation")[0]),
                          constraintMode=int(cu.get_from_tfile(tdir, "constraintMode")[0]),
                          regMode=int(cu.get_from_tfile(tdir, "regMode")[0]),
                          densityFlags=int(cu.get_from_tfile(tdir, "densityFlags")[0]),
                          distribution=str(cu.get_from_tfile(tdir, "distribution").GetTitle()),
                          axisSteering=str(cu.get_from_tfile(tdir, "axisSteering").GetTitle()))

    obj_names = cu.get_list_of_element_names(tdir)
    for name in obj_names:
        obj = tdir.Get(name)
        if name.startswith("background_reco_binning_"):
            bg_name = name.replace("background_reco_binning_", "").replace("_", " ")
            unfolder.backgrounds[bg_name] = obj

        elif name.startswith("background_gen_binning_"):
            bg_name = name.replace("background_gen_binning_", "").replace("_", " ")
            unfolder.backgrounds_gen_binning[bg_name] = obj

        elif name.startswith("syst_map_"):
            syst_name = name.replace("syst_map_", "").replace("_", " ")
            unfolder.syst_maps[syst_name] = obj

        elif name.startswith("syst_ematrix_"):
            syst_name = name.replace("syst_ematrix_", "").replace("_", " ")
            unfolder.syst_ematrices[syst_name] = obj

        elif name.startswith("syst_shift_"):
            syst_name = name.replace("syst_shift_", "").replace("_", " ")
            unfolder.syst_shifts[syst_name] = obj

        elif name.startswith("syst_shifted_unfolded_"):
            syst_name = name.replace("syst_shifted_unfolded_", "").replace("_", " ")
            unfolder.systs_shifted[syst_name] = obj

    for attr_name in MyUnfolder._simple_attr:
        obj = tdir.Get(attr_name)
        setattr(unfolder, attr_name, obj)

    return unfolder
