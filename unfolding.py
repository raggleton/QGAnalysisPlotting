#!/usr/bin/env python


"""TUnfold it all

Thanks to Ashley, Dennis
"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
from array import array
import numpy as np
import math
from itertools import product, chain
from copy import copy, deepcopy

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from my_unfolder import MyUnfolder, pickle_region, unpickle_region, ExpSystematic
from my_unfolder_plotter import MyUnfolderPlotter
from unfolding_regularisation_classes import TauScanner, LCurveScanner
from unfolding_config import get_dijet_config, get_zpj_config
from do_unfolding_plots import Setup, do_binned_plots_per_region_angle, do_all_big_normalised_1d_plots_per_region_angle
from unfolding_logistics import get_unfolding_argparser, get_unfolding_output_dir, sanitise_args

ROOT.gErrorIgnoreLevel = ROOT.kWarning
# ROOT.gErrorIgnoreLevel = ROOT.kInfo
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)  # VERY IMPORTANT - somewhere, closing a TFile for exp systs deletes a map...dunno why

# Control plot output format
OUTPUT_FMT = "pdf"


def rm_large_rel_error_bins(hist, relative_err_threshold=-1):
    """Reset bins in 2D hist to 0 if error/contents exceeds a certain value"""
    if relative_err_threshold < 0:
        return hist
    new_hist = hist.Clone()
    new_hist.SetDirectory(0)
    for ix in range(hist.GetNbinsX()+1):
        for iy in range(hist.GetNbinsY()+1):
            val = hist.GetBinContent(ix, iy)
            err = hist.GetBinError(ix, iy)
            if val == 0:
                continue
            if (1.*err/val) > relative_err_threshold:
                new_hist.SetBinContent(ix, iy, 0)
                new_hist.SetBinError(ix, iy, 0)
    return new_hist


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


def fill_empty_bins(response_map,
                    variable_bin_edges_reco,
                    variable_bin_edges_gen,
                    variable_name,
                    pt_bin_edges_reco,
                    pt_bin_edges_gen,
                    pt_bin_edges_underflow_reco,
                    pt_bin_edges_underflow_gen):
    """Fill in empty bins in the response map with some small value

    Testing to see if it helps stabilise unfolding
    """
    new_map = response_map.Clone()

    generator_binning, detector_binning = MyUnfolder.construct_tunfold_binning(variable_bin_edges_reco,
                                                                               variable_bin_edges_gen,
                                                                               variable_name,
                                                                               pt_bin_edges_reco,
                                                                               pt_bin_edges_gen,
                                                                               pt_bin_edges_underflow_reco,
                                                                               pt_bin_edges_underflow_gen)

    generator_binning_uflow = generator_binning.FindNode("generatordistribution_underflow")
    generator_binning_main = generator_binning.FindNode("generatordistribution")

    detector_binning_uflow = detector_binning.FindNode("detectordistribution_underflow")
    detector_binning_main = detector_binning.FindNode("detectordistribution")

    # in each gen pt, detector pt bin, find the largest and smallest bin counts
    # then go through again and fill in the empty bins with some very small value,
    # that is relatively smaller than the smallest bin is compared to the largest bin
    # (we're just trying to ensure non-zero, not physically super sensible)
    n_bins = (len(variable_bin_edges_gen)-1)*(len(variable_bin_edges_reco)-1)
    for ibin_pt_gen, gen_pt in enumerate(chain(pt_bin_edges_underflow_gen[:-1], pt_bin_edges_gen)):
        for ibin_pt_reco, reco_pt in enumerate(chain(pt_bin_edges_underflow_reco[:-1], pt_bin_edges_reco)):
            largest_bin_count = -99999
            smallest_bin_count = 9E99
            empty_bin_indices = []

            for gen_var in variable_bin_edges_gen[:-1]:
                # urghhhhhh have to manually choose the TUnfoldBinning object
                this_gen_binning = generator_binning_uflow if gen_pt < pt_bin_edges_gen[0] else generator_binning_main
                gen_bin_num = this_gen_binning.GetGlobalBinNumber(gen_var*1.0001, gen_pt*1.0001)
                for reco_var in variable_bin_edges_reco[:-1]:
                    this_reco_binning = detector_binning_uflow if reco_pt < pt_bin_edges_reco[0] else detector_binning_main
                    reco_bin_num = this_reco_binning.GetGlobalBinNumber(reco_var*1.0001, reco_pt*1.0001)
                    val = response_map.GetBinContent(gen_bin_num, reco_bin_num)
                    if val == 0:
                        empty_bin_indices.append([gen_bin_num, reco_bin_num])
                    elif val < smallest_bin_count:
                        smallest_bin_count = val
                    elif val > largest_bin_count:
                        largest_bin_count = val

            # if ibin_pt_gen < 4 and ibin_pt_reco < 4:
            #     print(ibin_pt_gen, ibin_pt_reco, ":", len(empty_bin_indices), n_bins, largest_bin_count, smallest_bin_count, smallest_bin_count/largest_bin_count)

            is_near_diagonal = abs(ibin_pt_gen - (ibin_pt_reco/2)) < 3
            if (largest_bin_count > 0) and (smallest_bin_count > 0) and (((1. * len(empty_bin_indices) / n_bins) < 0.75) or is_near_diagonal):
                ratio = smallest_bin_count / largest_bin_count
                new_val = ratio * ratio * smallest_bin_count
                print("Filling", ibin_pt_gen, ibin_pt_reco, smallest_bin_count * ratio)
                for x, y in empty_bin_indices:
                    new_map.SetBinContent(x, y, new_val)
                    new_map.SetBinError(x, y, 0)

    return new_map


# To be able to export to XML, since getting std::ostream from python is impossible?
my_binning_xml_code = """
class BinningXMLExporter {
public:
    BinningXMLExporter() {cout << " Creating BinningXMLExporter " << endl;}

    static Int_t ExportXML(const TUnfoldBinning &binning,
                           std::string dirname,
                           std::string filename,
                           Bool_t writeHeader,
                           Bool_t writeFooter,
                           Int_t indent=0) {

        ofstream dtdFile(TString::Format("%s/tunfoldbinning.dtd", dirname.c_str()));
        TUnfoldBinningXML::WriteDTD(dtdFile);
        dtdFile.close();

        ofstream xmlfile;
        xmlfile.open(TString::Format("%s/%s", dirname.c_str(), filename.c_str()), ios::out);
        Int_t result = TUnfoldBinningXML::ExportXML(binning, xmlfile, writeHeader, writeFooter, indent);
        xmlfile.close();
        return result;
    }
};
"""
ROOT.gInterpreter.ProcessLine(my_binning_xml_code)


if __name__ == "__main__":
    parser = get_unfolding_argparser(description=__doc__)
    args = parser.parse_args()
    print("")
    print(args)
    print("")
    sanitise_args(args)
    output_dir = get_unfolding_output_dir(args)
    cu.check_dir_exists_create(output_dir)
    print("")
    print(args)

    print("")
    print("Outputting to", output_dir)
    print("")

    src_dir = args.source

    # Setup files and regions to unfold
    # --------------------------------------------------------------------------
    regions = []

    # Configure Dijet regions
    if any([args.doDijetCentral, args.doDijetForward, args.doDijetCentralGroomed, args.doDijetForwardGroomed]):
        # actually these are all pretty similar...
        tau_limits_central = {
            'jet_puppiMultiplicity': (1E-9, 1E-6) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_pTD': (1E-12, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            # 'jet_LHA': (1E-11, 1E-7) if args.regularizeAxis == 'angle' else (1E-13, 1E-10),
            'jet_LHA': (1E-13, 1E-9) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_width': (1E-12, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_thrust': (1E-12, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_puppiMultiplicity_charged': (1E-12, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_pTD_charged': (1E-12, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_LHA_charged': (1E-10, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_width_charged': (1E-13, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_thrust_charged': (1E-12, 1E-9) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
        }

        tau_limits_central_groomed = {
            'jet_puppiMultiplicity': (1E-9, 1E-6) if args.regularizeAxis == "angle" else (1E-1, 1E3),
            'jet_pTD': (1E-12, 1E-8) if args.regularizeAxis == "angle" else (1E-1, 1E3),
            'jet_LHA': (1E-11, 1E-8) if args.regularizeAxis == 'angle' else (1E-13, 1E-10),
            'jet_width': (1E-12, 1E-8) if args.regularizeAxis == "angle" else (1E-1, 1E3),
            'jet_thrust': (1E-12, 1E-8) if args.regularizeAxis == "angle" else (1E-1, 1E3),
            'jet_puppiMultiplicity_charged': (1E-12, 1E-8) if args.regularizeAxis == "angle" else (1E-1, 1E3),
            'jet_pTD_charged': (1E-12, 1E-8) if args.regularizeAxis == "angle" else (1E-1, 1E3),
            'jet_LHA_charged': (1E-10, 1E-8) if args.regularizeAxis == "angle" else (1E-1, 1E3),
            'jet_width_charged': (1E-13, 1E-8) if args.regularizeAxis == "angle" else (1E-1, 1E3),
            'jet_thrust_charged': (1E-12, 1E-9) if args.regularizeAxis == "angle" else (1E-1, 1E3),
        }

        tau_limits_forward = {
            'jet_puppiMultiplicity': (1E-9, 1E-6) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_pTD': (1E-12, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_LHA': (1E-11, 1E-8) if args.regularizeAxis == 'angle' else (1E-13, 1E-10),
            'jet_width': (1E-12, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_thrust': (1E-12, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_puppiMultiplicity_charged': (1E-12, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_pTD_charged': (1E-12, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_LHA_charged': (1E-10, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_width_charged': (1E-13, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_thrust_charged': (1E-12, 1E-9) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
        }

        tau_limits_forward_groomed = {
            'jet_puppiMultiplicity': (1E-9, 1E-6) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_pTD': (1E-12, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_LHA': (1E-11, 1E-8) if args.regularizeAxis == 'angle' else (1E-13, 1E-10),
            'jet_width': (1E-12, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_thrust': (1E-12, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_puppiMultiplicity_charged': (1E-12, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_pTD_charged': (1E-12, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_LHA_charged': (1E-10, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_width_charged': (1E-13, 1E-8) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_thrust_charged': (1E-12, 1E-9) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
        }

        if args.doDijetCentral:
            dijet_region_central_dict = get_dijet_config(src_dir, central=True, groomed=False)
            dijet_region_central_dict['tau_limits'] = tau_limits_central
            regions.append(dijet_region_central_dict)

        if args.doDijetForward:
            dijet_region_forward_dict = get_dijet_config(src_dir, central=False, groomed=False)
            dijet_region_forward_dict['tau_limits'] = tau_limits_forward
            regions.append(dijet_region_forward_dict)

        if args.doDijetCentralGroomed:
            dijet_region_central_groomed_dict = get_dijet_config(src_dir, central=True, groomed=True)
            dijet_region_central_groomed_dict['tau_limits'] = tau_limits_central_groomed
            regions.append(dijet_region_central_groomed_dict)

        if args.doDijetForwardGroomed:
            dijet_region_forward_groomed_dict = get_dijet_config(src_dir, central=False, groomed=True)
            dijet_region_forward_groomed_dict['tau_limits'] = tau_limits_forward_groomed
            regions.append(dijet_region_forward_groomed_dict)

    if any([args.doZPJ, args.doZPJGroomed]):
        # FOR Z+JETS:
        tau_limits = {
            'jet_puppiMultiplicity': (1E-5, 1E-2) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_pTD': (1E-6, 1E-4) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            # 'jet_LHA': (1E-5, 1E-3) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            # 'jet_LHA': (1E-7, 1E-4) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            # 'jet_LHA': (1E-8, 1E-4) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            # 'jet_LHA': (1E-6, 1E-3) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_LHA': (1E-2, 1E2) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_width': (1E-5, 1E-2) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_thrust': (1E-6, 1E-2) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_puppiMultiplicity_charged': (1E-6, 1E-2) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_pTD_charged': (1E-6, 1E-2) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_LHA_charged': (1E-5, 1E-2) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_width_charged': (1E-6, 1E-3) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_thrust_charged': (1E-8, 1E-5) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
        }

        tau_limits_groomed = {
            'jet_puppiMultiplicity': (1E-5, 1E-2) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_pTD': (1E-7, 1E-3) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_LHA': (1E-5, 1E-3) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_width': (1E-5, 1E-2) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_thrust': (1E-6, 1E-2) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_puppiMultiplicity_charged': (1E-5, 1E-2) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_pTD_charged': (1E-6, 1E-2) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_LHA_charged': (1E-5, 1E-2) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_width_charged': (1E-6, 1E-3) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
            'jet_thrust_charged': (1E-7, 1E-5) if args.regularizeAxis == 'angle' else (1E-1, 1E3),
        }

        if args.doZPJ:
            zpj_region_dict = get_zpj_config(src_dir, groomed=False)
            zpj_region_dict['tau_limits'] = tau_limits
            regions.append(zpj_region_dict)

        if args.doZPJGroomed:
            zpj_region_groomed_dict = get_zpj_config(src_dir, groomed=True)
            zpj_region_groomed_dict['tau_limits'] = tau_limits_groomed
            regions.append(zpj_region_groomed_dict)

    # Setup various options
    # --------------------------------------------------------------------------
    REGULARIZE = args.regularize
    MC_INPUT = args.MCinput
    MC_SPLIT = args.MCsplit if args.MCinput else False
    SUBTRACT_FAKES = True  # this should alwys be True

    jet_algo = "AK4 PF PUPPI"
    if "ak8puppi" in src_dir:
        jet_algo = "AK8 PF PUPPI"

    if args.angles[0] == "all":
        angles = qgc.COMMON_VARS
    else:
        angles = [a for a in qgc.COMMON_VARS if a.var in args.angles]

    print("Running TUnfold version", ROOT.TUnfold.GetTUnfoldVersion())

    LAMBDA_VAR_DICTS = qgc.VAR_UNFOLD_DICT
    if 'target0p5' in src_dir:
        LAMBDA_VAR_DICTS = qgc.VAR_UNFOLD_DICT_TARGET0p5
    elif 'target0p6' in src_dir:
        LAMBDA_VAR_DICTS = qgc.VAR_UNFOLD_DICT_TARGET0p6


    # Do unfolding per signal region
    # --------------------------------------------------------------------------
    for orig_region in regions[:]:

        # Setup pt bins
        # -------------
        # need different ones for Z+Jets region
        is_zpj = "ZPlusJets" in orig_region['name']

        zpj_append = "_zpj" if is_zpj else ""

        pt_bin_edges_gen = qgc.PT_UNFOLD_DICT['signal%s_gen' % (zpj_append)]
        pt_bin_edges_reco = qgc.PT_UNFOLD_DICT['signal%s_reco' % (zpj_append)]

        pt_bin_edges_underflow_gen = qgc.PT_UNFOLD_DICT['underflow%s_gen' % (zpj_append)]
        pt_bin_edges_underflow_reco = qgc.PT_UNFOLD_DICT['underflow%s_reco' % (zpj_append)]

        # Modify systematics as necessary
        # ----------------------------------------------------------------------

        if args.doExperimentalSysts or args.doExperimentalSystsFromFile:
            # Remove the lumi one if we have no backgrounds, or the user has not said to remove backgrounds
            orig_region['experimental_systematics'] = [syst_dict for syst_dict in orig_region['experimental_systematics']
                                                       if not ('lumi' in syst_dict['label'].lower()
                                                               and (len(orig_region.get('backgrounds', [])) == 0
                                                                    or not args.subtractBackgrounds))]
            if args.doExperimentalSystsOnlyHerwig:
                # only herwig related systs
                orig_region['experimental_systematics'] = [s for s in orig_region['experimental_systematics']
                                                           if 'herwig' in s['label'].lower()]

        else:
            orig_region['experimental_systematics'] = []

        if not (args.doScaleSysts or args.doScaleSystsFromFile):
            orig_region['scale_systematics'] = []

        if not (args.doModelSysts or args.doModelSystsFromFile):
            orig_region['model_systematics'] = []

        elif args.doModelSystsOnlyHerwig:
            # only herwig related systs
            orig_region['model_systematics'] = [s for s in orig_region['model_systematics']
                                                if 'herwig' in s['label'].lower()]
        elif args.doModelSystsOnlyScale:
            # only scale related systs
            orig_region['model_systematics'] = [s for s in orig_region['model_systematics']
                                                if 'mu' in s['label'].lower()]

        if not (args.doPDFSysts or args.doPDFSystsFromFile):
            orig_region['pdf_systematics'] = []

        # Do 1D unfolding of pt
        # TODO!
        # ----------------------------------------------------------------------
        # hist_data_reco = cu.get_from_tfile(orig_region['data_tfile'], "%s/hist_pt_reco_all" % (orig_region['dirname']))
        # mc_hname_append = "split" if MC_SPLIT else "all"
        # if isinstance(orig_region['mc_tfile'], str):
        #     orig_region['mc_tfile'] = cu.open_root_file(orig_region['mc_tfile'])
        # hist_mc_reco = cu.get_from_tfile(orig_region['mc_tfile'], "%s/hist_pt_reco_%s" % (orig_region['dirname'], mc_hname_append))
        # hist_mc_gen = cu.get_from_tfile(orig_region['mc_tfile'], "%s/hist_pt_truth_%s" % (orig_region['dirname'], mc_hname_append))
        # hist_mc_gen_pt = cu.get_from_tfile(orig_region['mc_tfile'], "%s/hist_pt_truth_%s" % (orig_region['dirname'], mc_hname_append))
        # hist_mc_gen_reco_map = cu.get_from_tfile(orig_region['mc_tfile'], "%s/tu_pt_GenReco_%s" % (orig_region['dirname'], mc_hname_append))

        # Do unfolding for each angle
        # ----------------------------------------------------------------------
        for angle in angles:
            region = copy(orig_region)  # make copy since we might modify it later, e.g. PDF, and want same start for each angle

            if isinstance(region['mc_tfile'], str):
                region['mc_tfile'] = cu.open_root_file(region['mc_tfile'])

            angle_prepend = "groomed " if "groomed" in region['name'] else ""
            append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name
            this_angle_name = angle.name
            if (angle_prepend != ""
                and this_angle_name != 'LHA'
                and "_{T}" not in this_angle_name
                and "PUPPI" not in this_angle_name):
                # lower case if Groomed..., but be careful of e.g. pTD, LHA
                this_angle_name = this_angle_name[0].lower() + this_angle_name[1:]
            angle_str = "%s%s (%s)" % (angle_prepend, this_angle_name, angle.lambda_str)

            print("*"*120)
            print("Region/var: %s" % (append))
            print("*"*120)

            # put plots in subdir, to avoid overcrowding
            this_output_dir = "%s/%s/%s" % (output_dir, region['name'], angle.var)
            cu.check_dir_exists_create(this_output_dir)

            # Save minimal set of hists etc to ROOT file for access later
            # For everything, we pickle it
            output_tfile_slim = ROOT.TFile("%s/unfolding_result_slim.root" % (this_output_dir), "RECREATE")
            new_tdir = "%s/%s" % (region['name'], angle.var)
            output_tfile_slim.mkdir(new_tdir)
            this_slim_tdir = output_tfile_slim.Get(new_tdir)
            this_slim_tdir.cd()

            # Setup MyUnfolder object to do unfolding etc
            # -------------------------------------------
            angle_bin_edges_reco = LAMBDA_VAR_DICTS[angle.var]['reco']
            angle_bin_edges_gen = LAMBDA_VAR_DICTS[angle.var]['gen']
            angle_shortname = angle.var.replace("jet_", "")

            mc_hname_append = "split" if MC_SPLIT else "all"
            hist_mc_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_%s" % (region['dirname'], angle_shortname, mc_hname_append))
            hist_mc_gen = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_truth_%s" % (region['dirname'], angle_shortname, mc_hname_append))
            hist_mc_gen_reco_map = cu.get_from_tfile(region['mc_tfile'], "%s/tu_%s_GenReco_%s" % (region['dirname'], angle_shortname, mc_hname_append))

            hist_data_reco = None
            if not MC_INPUT:
                if not isinstance(region['data_tfile'], ROOT.TFile):
                    region['data_tfile'] = cu.open_root_file(region['data_tfile'])
                hist_data_reco = cu.get_from_tfile(region['data_tfile'], "%s/hist_%s_reco_all" % (region['dirname'], angle_shortname))

            # Need to scale if using H++ as input
            # hist_mc_gen_reco_map.Scale(1E8)
            # hist_mc_gen.Scale(1E8)
            # hist_mc_reco.Scale(1E8)

            # Actual distribution to be unfolded
            reco_1d = hist_mc_reco.Clone() if MC_INPUT else hist_data_reco

            hist_fakes_reco_fraction = None

            if SUBTRACT_FAKES:
                # to construct our "fakes" template, we use the ratio as predicted by MC, and apply it to data
                # this way we ensure we don't have -ve values, and avoid any issue with cross sections
                hist_mc_fakes_reco = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_fake_%s" % (region['dirname'], angle_shortname, mc_hname_append))
                # hist_mc_fakes_reco.Scale(1E8)
                hist_fakes_reco = hist_mc_fakes_reco.Clone("hist_%s_fakes" % angle_shortname)
                hist_fakes_reco.Divide(hist_mc_reco)

                # plot fake fraction before multiplting by 'data'
                hist_fakes_reco_fraction = hist_fakes_reco.Clone("hist_fakes_reco_fraction")
                hist_fakes_reco.Multiply(reco_1d)

                # background-subtracted reco hists, only for plotting purposes, not for TUnfold (that does background subtraction internally)
                reco_1d_bg_subtracted = reco_1d.Clone()
                reco_1d_bg_subtracted.Add(hist_fakes_reco, -1)

                if not MC_INPUT:
                    hist_data_reco_bg_subtracted = hist_data_reco.Clone(hist_data_reco.GetName() + "_bgrSubtracted")
                    hist_data_reco_bg_subtracted.Add(hist_fakes_reco, -1)

                hist_mc_reco_bg_subtracted = hist_mc_reco.Clone(hist_mc_reco.GetName() + "_bgrSubtracted")
                hist_mc_reco_bg_subtracted.Add(hist_mc_fakes_reco, -1)  # should this be hist_fakes_reco? depends on what we want to see...

            mc_hname_append = "_split" if MC_SPLIT else ""  # FIXME consistency in unfold hist module!
            hist_data_reco_gen_binning = None
            if not MC_INPUT:
                hist_data_reco_gen_binning = cu.get_from_tfile(region['data_tfile'], "%s/hist_%s_reco_gen_binning" % (region['dirname'], angle_shortname))
            hist_mc_reco_gen_binning = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_gen_binning%s" % (region['dirname'], angle_shortname, mc_hname_append))

            # Actual distribution to be unfolded, but with gen binning
            reco_1d_gen_binning = hist_mc_reco_gen_binning.Clone() if MC_INPUT else hist_data_reco_gen_binning

            hist_fakes_reco_gen_binning = None
            hist_fakes_reco_fraction_gen_binning = None
            if SUBTRACT_FAKES:
                mc_hname_append = "_split" if MC_SPLIT else ""  # FIXME consistency in unfold hist module!
                hist_mc_fakes_reco_gen_binning = cu.get_from_tfile(region['mc_tfile'], "%s/hist_%s_reco_fake_gen_binning%s" % (region['dirname'], angle_shortname, mc_hname_append))
                # create template as above, but with gen binning
                hist_fakes_reco_gen_binning = hist_mc_fakes_reco_gen_binning.Clone("hist_%s_fakes_gen_binning" % angle_shortname)
                hist_fakes_reco_gen_binning.Divide(hist_mc_reco_gen_binning)
                hist_fakes_reco_fraction_gen_binning = hist_fakes_reco_gen_binning.Clone("hist_fakes_reco_fraction_gen_binning")

                hist_fakes_reco_gen_binning.Multiply(reco_1d_gen_binning)

                # background-subtracted reco hists, only for plotting purposes, not for TUnfold (that does background subtraction internally)
                reco_1d_gen_binning_bg_subtracted = reco_1d_gen_binning.Clone()
                reco_1d_gen_binning_bg_subtracted.Add(hist_fakes_reco_gen_binning, -1)

                if not MC_INPUT:
                    hist_data_reco_gen_binning_bg_subtracted = hist_data_reco_gen_binning.Clone(hist_data_reco_gen_binning.GetName() + "_bgrSubtracted")
                    hist_data_reco_gen_binning_bg_subtracted.Add(hist_fakes_reco_gen_binning, -1)

                hist_mc_reco_gen_binning_bg_subtracted = hist_mc_reco_gen_binning.Clone(hist_mc_reco_gen_binning.GetName() + "_bgrSubtracted")
                hist_mc_reco_gen_binning_bg_subtracted.Add(hist_mc_fakes_reco_gen_binning, -1)  # should this be hist_fakes_reco_gen_binning? depends on what we want to see...

            # Setup unfolder object
            # ---------------------
            variable_name = "%s%s" % (angle_prepend, angle.name)
            axis_steering = '*[B]'
            if args.regularizeAxis == 'pt':
                axis_steering = 'pt[B];%s[N]' % variable_name
            elif args.regularizeAxis == 'angle':
                axis_steering = 'pt[N];%s[B]' % variable_name

            # setup one unfolder to get binning
            # TODO: move creation of TUnfoldbinning into own func
            dummy_unfolder = MyUnfolder(response_map=rm_large_rel_error_bins(hist_mc_gen_reco_map),
                                        variable_bin_edges_reco=angle_bin_edges_reco,
                                        variable_bin_edges_gen=angle_bin_edges_gen,
                                        variable_name=variable_name,
                                        pt_bin_edges_reco=pt_bin_edges_reco,
                                        pt_bin_edges_gen=pt_bin_edges_gen,
                                        pt_bin_edges_underflow_reco=pt_bin_edges_underflow_reco,
                                        pt_bin_edges_underflow_gen=pt_bin_edges_underflow_gen,
                                        orientation=ROOT.TUnfold.kHistMapOutputHoriz,
                                        constraintMode=args.areaConstraint,
                                        # regMode=ROOT.TUnfold.kRegModeCurvature,
                                        # densityFlags=ROOT.TUnfoldDensity.kDensityModeBinWidth, # important as we have varying bin sizes!
                                        regMode=ROOT.TUnfold.kRegModeNone,
                                        densityFlags=ROOT.TUnfoldDensity.kDensityModeBinWidthAndUser, # doesn't actually matter as RegModNone
                                        distribution='generatordistribution',  # the one to use for actual final regularisation/unfolding
                                        axisSteering=axis_steering)

            # fill in empty bins in response_map to ensure better unfolding
            new_response_map = fill_empty_bins(hist_mc_gen_reco_map,
                                               variable_bin_edges_reco=angle_bin_edges_reco,
                                               variable_bin_edges_gen=angle_bin_edges_gen,
                                               variable_name=variable_name,
                                               pt_bin_edges_reco=pt_bin_edges_reco,
                                               pt_bin_edges_gen=pt_bin_edges_gen,
                                               pt_bin_edges_underflow_reco=pt_bin_edges_underflow_reco,
                                               pt_bin_edges_underflow_gen=pt_bin_edges_underflow_gen,)


            unfolder = MyUnfolder(response_map=rm_large_rel_error_bins(hist_mc_gen_reco_map),
                                  variable_bin_edges_reco=angle_bin_edges_reco,
                                  variable_bin_edges_gen=angle_bin_edges_gen,
                                  variable_name=variable_name,
                                  pt_bin_edges_reco=pt_bin_edges_reco,
                                  pt_bin_edges_gen=pt_bin_edges_gen,
                                  pt_bin_edges_underflow_reco=pt_bin_edges_underflow_reco,
                                  pt_bin_edges_underflow_gen=pt_bin_edges_underflow_gen,
                                  orientation=ROOT.TUnfold.kHistMapOutputHoriz,
                                  constraintMode=args.areaConstraint,
                                  # regMode=ROOT.TUnfold.kRegModeCurvature,
                                  # densityFlags=ROOT.TUnfoldDensity.kDensityModeBinWidth, # important as we have varying bin sizes!
                                  regMode=ROOT.TUnfold.kRegModeNone,
                                  densityFlags=ROOT.TUnfoldDensity.kDensityModeBinWidthAndUser, # doesn't actually matter as RegModNone
                                  distribution='generatordistribution',  # the one to use for actual final regularisation/unfolding
                                  axisSteering=axis_steering)

            # Save binning to file
            unfolder.save_binning(txt_filename="%s/binning_scheme.txt" % (this_output_dir), print_xml=False)
            ROOT.BinningXMLExporter.ExportXML(unfolder.detector_binning, this_output_dir, "detector_binning.xml", True, True, 2)
            ROOT.BinningXMLExporter.ExportXML(unfolder.generator_binning, this_output_dir, "generator_binning.xml", True, True, 2)

            unfolder_plotter = MyUnfolderPlotter(unfolder, is_data=not MC_INPUT)
            plot_args = dict(output_dir=this_output_dir, append=append)

            is_herwig = "Herwig" in region['mc_label']
            is_pythia8 = region['mc_label'] == "Pythia8"  # not MG+Pythia9
            if is_herwig or is_pythia8:
                # SetEpsMatrix ensures rank properly calculated when inverting
                # Needed if you get message "rank of matrix E 55 expect 170"
                # And unfolded looks wacko
                unfolder.SetEpsMatrix(1E-18)

            eps_matrix = 1E-160
            unfolder.SetEpsMatrix(eps_matrix)
            print("Running with eps =", unfolder.GetEpsMatrix())

            # Set what is to be unfolded
            # ------------------------------------------------------------------
            unfolder.set_input(input_hist=reco_1d,
                               input_hist_gen_binning=reco_1d_gen_binning,
                               hist_truth=hist_mc_gen,
                               hist_mc_reco=hist_mc_reco,
                               hist_mc_reco_bg_subtracted=hist_mc_reco_bg_subtracted,
                               hist_mc_reco_gen_binning=hist_mc_reco_gen_binning,
                               hist_mc_reco_gen_binning_bg_subtracted=hist_mc_reco_gen_binning_bg_subtracted,
                               bias_factor=args.biasFactor)

            unfolder.hist_bin_chopper.add_obj('hist_truth', unfolder.hist_truth)

            # Add systematic errors as different response matrices
            # ------------------------------------------------------------------
            # if args.doExperimentalSysts:
            #     chosen_rsp_bin = (35, 35)
            #     print("nominal response bin content for", chosen_rsp_bin, unfolder.response_map.GetBinContent(*chosen_rsp_bin))

            #     for syst_ind, syst_dict in enumerate(region['experimental_systematics']):
            #         print("Adding systematic:", syst_dict['label'])
            #         if 'factor' in syst_dict:
            #             # special case for e.g. lumi - we construct a reponse hist, and add it using relative mode
            #             rel_map = unfolder.response_map.Clone(syst_dict['label']+"Map")
            #             for xbin, ybin in product(range(1, rel_map.GetNbinsX()+1), range(1, rel_map.GetNbinsY()+1)):
            #                 rel_map.SetBinContent(xbin, ybin, syst_dict['factor'])
            #                 rel_map.SetBinError(xbin, ybin, 0)
            #             unfolder.add_sys_error(rel_map, syst_dict['label'], ROOT.TUnfoldDensity.kSysErrModeRelative)

            #         # elif 'herwig' in syst_dict['label'].lower():
            #         #     if not isinstance(syst_dict['tfile'], ROOT.TFile):
            #         #         syst_dict['tfile'] = cu.open_root_file(syst_dict['tfile'])
            #         #     map_syst = cu.get_from_tfile(syst_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
            #         #     map_syst.Scale(unfolder.response_map.Integral() / map_syst.Integral())
            #         #     print("    syst bin", chosen_rsp_bin, map_syst.GetBinContent(*chosen_rsp_bin))
            #         #     unfolder.add_sys_error(rm_large_rel_error_bins(map_syst), syst_dict['label'], ROOT.TUnfoldDensity.kSysErrModeMatrix)

            #         else:
            #             if not isinstance(syst_dict['tfile'], ROOT.TFile):
            #                 syst_dict['tfile'] = cu.open_root_file(syst_dict['tfile'])
            #             map_syst = cu.get_from_tfile(syst_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
            #             print("    syst bin", chosen_rsp_bin, map_syst.GetBinContent(*chosen_rsp_bin))
            #             unfolder.add_sys_error(rm_large_rel_error_bins(map_syst), syst_dict['label'], ROOT.TUnfoldDensity.kSysErrModeMatrix)

            #         # Plot the reponse matrix for this systematic
            #         syst_label_no_spaces = cu.no_space_str(syst_dict['label'])
            #         output_filename = "%s/response_map_syst_%s_%s.%s" % (this_output_dir, syst_label_no_spaces, append, unfolder_plotter.output_fmt)
            #         title = "%s, %s region, %s, %s" % (jet_algo, region['label'], angle_str, syst_dict['label'])
            #         unfolder_plotter.draw_2d_hist(unfolder.get_exp_syst(syst_dict['label']).syst_map,
            #                                       title=title,
            #                                       output_filename=output_filename,
            #                                       logz=True,
            #                                       draw_bin_lines_x=True,
            #                                       draw_bin_lines_y=True,
            #                                       canvas_size=(800, 700),
            #                                       xtitle='Generator bin',
            #                                       ytitle='Detector bin')

            # Subtract fakes (treat as background)
            # ------------------------------------------------------------------
            if SUBTRACT_FAKES:
                unfolder.subtract_background(hist_fakes_reco, "Signal fakes", scale=1., scale_err=0.0)
                unfolder.subtract_background_gen_binning(hist_fakes_reco_gen_binning, "Signal fakes", scale=1., scale_err=0.0)

            # Subtract actual backgrounds if necessary
            # ------------------------------------------------------------------
            background_reco_1d = None
            background_gen_1d = None
            if "backgrounds" in region and args.subtractBackgrounds:
                for bg_ind, bg_dict in enumerate(region['backgrounds']):
                    print("Subtracting", bg_dict['name'], 'background')
                    if not isinstance(bg_dict['tfile'], ROOT.TFile):
                        bg_dict['tfile'] = cu.open_root_file(bg_dict['tfile'])

                    mc_hname_append = "split" if MC_SPLIT else "all"
                    bg_hist = cu.get_from_tfile(bg_dict['tfile'], "%s/hist_%s_reco_%s" % (region['dirname'], angle_shortname, mc_hname_append))
                    bg_hist_gen = cu.get_from_tfile(bg_dict['tfile'], "%s/hist_%s_truth_%s" % (region['dirname'], angle_shortname, mc_hname_append))
                    bg_dict['hist'] = bg_hist
                    bg_dict['hist_gen'] = bg_hist

                    # keep one big hist
                    if not background_reco_1d:
                        background_reco_1d = bg_hist.Clone()
                        background_gen_1d = bg_hist_gen.Clone()
                    else:
                        background_reco_1d.Add(bg_hist, bg_dict.get('rate', 1.))
                        background_gen_1d.Add(bg_hist_gen, bg_dict.get('rate', 1.))

                    unfolder.subtract_background(hist=bg_hist,
                                                 name=bg_dict['name'],
                                                 scale=bg_dict.get('rate', 1.),
                                                 scale_err=bg_dict.get('rate_unc', 0.))

            title = "%s\n%s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_background_fractions(title=title, **plot_args)

            # Do any regularization
            # ---------------------
            unfolder.print_condition_number()

            unreg_unfolder = None
            unreg_unfolded_1d = None
            L_matrix_entries = []
            if REGULARIZE != "None":
                print("Doing preliminary unregularised unfolding...")
                # To setup the L matrix correctly, we have to rescale
                # the default one by the inverse of the pT spectrum
                # This means we first need a copy of the L matrix, so make a
                # dummy unfolder, ensuring it's setup to make L
                # We also need to do an unregularised unfolding first to get
                # the correct pt factors, since data spectrum != MC
                unreg_unfolder = MyUnfolder(response_map=unfolder.response_map,
                                            variable_bin_edges_reco=unfolder.variable_bin_edges_reco,
                                            variable_bin_edges_gen=unfolder.variable_bin_edges_gen,
                                            variable_name=unfolder.variable_name,
                                            pt_bin_edges_reco=unfolder.pt_bin_edges_reco,
                                            pt_bin_edges_gen=unfolder.pt_bin_edges_gen,
                                            pt_bin_edges_underflow_reco=unfolder.pt_bin_edges_underflow_reco,
                                            pt_bin_edges_underflow_gen=unfolder.pt_bin_edges_underflow_gen,
                                            orientation=unfolder.orientation,
                                            constraintMode=unfolder.constraintMode,
                                            regMode=ROOT.TUnfold.kRegModeCurvature,
                                            densityFlags=ROOT.TUnfoldDensity.kDensityModeBinWidth, # important as we have varying bin sizes!
                                            distribution=unfolder.distribution,
                                            axisSteering=unfolder.axisSteering)

                unreg_unfolder.SetEpsMatrix(eps_matrix)

                # Do the unregularised unfolding to get an idea of bin contents
                # and uncertainties
                # Set what is to be unfolded
                # ------------------------------------------------------------------
                unreg_unfolder.set_input(input_hist=reco_1d,
                                         input_hist_gen_binning=reco_1d_gen_binning,
                                         hist_truth=hist_mc_gen,
                                         hist_mc_reco=hist_mc_reco,
                                         hist_mc_reco_bg_subtracted=hist_mc_reco_bg_subtracted,
                                         hist_mc_reco_gen_binning=hist_mc_reco_gen_binning,
                                         hist_mc_reco_gen_binning_bg_subtracted=hist_mc_reco_gen_binning_bg_subtracted,
                                         bias_factor=0)


                # For now, ignore experimental systematics

                # Subtract fakes (treat as background)
                # ------------------------------------------------------------------
                if SUBTRACT_FAKES:
                    unreg_unfolder.subtract_background(hist_fakes_reco, "Signal fakes", scale=1., scale_err=0.0)
                    unreg_unfolder.subtract_background_gen_binning(hist_fakes_reco_gen_binning, "Signal fakes", scale=1., scale_err=0.0)

                # Subtract actual backgrounds if necessary
                # ------------------------------------------------------------------
                if "backgrounds" in region and args.subtractBackgrounds:
                    for bg_ind, bg_dict in enumerate(region['backgrounds']):
                        print("Subtracting", bg_dict['name'], 'background')
                        if not isinstance(bg_dict['tfile'], ROOT.TFile):
                            bg_dict['tfile'] = cu.open_root_file(bg_dict['tfile'])

                        mc_hname_append = "split" if MC_SPLIT else "all"
                        bg_hist = cu.get_from_tfile(bg_dict['tfile'], "%s/hist_%s_reco_%s" % (region['dirname'], angle_shortname, mc_hname_append))
                        bg_hist_gen = cu.get_from_tfile(bg_dict['tfile'], "%s/hist_%s_truth_%s" % (region['dirname'], angle_shortname, mc_hname_append))
                        bg_dict['hist'] = bg_hist
                        bg_dict['hist_gen'] = bg_hist

                        unreg_unfolder.subtract_background(hist=bg_hist,
                                                     name=bg_dict['name'],
                                                     scale=bg_dict.get('rate', 1.),
                                                     scale_err=bg_dict.get('rate_unc', 0.))

                unreg_unfolder.do_unfolding(0)
                unreg_unfolder.get_output(hist_name="unreg_unfolded_1d")
                unreg_unfolder._post_process()
                unreg_unfolded_1d = unreg_unfolder.unfolded

                region['unreg_unfolder'] = unreg_unfolder

                orig_Lmatrix = unreg_unfolder.GetL("orig_Lmatrix_%s" % (append), "", unreg_unfolder.use_axis_binning)
                xax = orig_Lmatrix.GetXaxis()
                # Get bin factors from an unregularised unfolding first,
                # to compensate for the fact that the shape differs between data & MC
                bin_factors = unreg_unfolder.calculate_pt_bin_factors(which='gen') # calculate factors to get uniform pt spectrum
                # bin_factors = unreg_unfolder.calculate_pt_bin_factors(which='unfolded') # calculate factors to get uniform pt spectrum
                bin_widths = unreg_unfolder.get_gen_bin_widths() # mapping {global bin number : (lambda bin width, pt bin width)}

                print(unreg_unfolder.variable_bin_edges_gen)

                ref_hist = unreg_unfolder.hist_truth
                unfolded_max = ref_hist.GetMaximum()

                gen_node = unfolder.generator_binning.FindNode('generatordistribution')
                for ilambda in range(len(unfolder.variable_bin_edges_gen[:-1])):
                    for ipt in range(len(unfolder.pt_bin_edges_gen[:-3])):
                        pt_cen = unfolder.pt_bin_edges_gen[ipt+1] + 0.000001  # add a tiny bit to make sure we're in the bin properly (I can never remember if included or not)
                        # lambda_cen = unfolder.variable_bin_edges_gen[ilambda+1] + 0.000001  # add a tiny bit to make sure we're in the bin properly (I can never remember if included or not)
                        lambda_cen = unfolder.variable_bin_edges_gen[ilambda] + 0.000001  # add a tiny bit to make sure we're in the bin properly (I can never remember if included or not)

                        bin_ind_pt_down = gen_node.GetGlobalBinNumber(lambda_cen, unfolder.pt_bin_edges_gen[ipt] + 0.000001)
                        bin_ind_pt_up = gen_node.GetGlobalBinNumber(lambda_cen, unfolder.pt_bin_edges_gen[ipt+2] + 0.000001)

                        bin_ind_cen = gen_node.GetGlobalBinNumber(lambda_cen, pt_cen)

                        # bin_ind_var_down = gen_node.GetGlobalBinNumber(unfolder.variable_bin_edges_gen[ilambda], pt_cen)
                        # bin_ind_var_up = gen_node.GetGlobalBinNumber(unfolder.variable_bin_edges_gen[ilambda+2], pt_cen)

                        # print("Adding L matrix entry")
                        # print(lambda_cen, (unfolder.pt_bin_edges_gen[ipt], unfolder.pt_bin_edges_gen[ipt+1], unfolder.pt_bin_edges_gen[ipt+2]))
                        # print(lambda_cen, (unfolder.pt_bin_edges_gen[ipt], unfolder.pt_bin_edges_gen[ipt+1], unfolder.pt_bin_edges_gen[ipt+2]))
                        # pt_bin_width_down = pt_bin_edges_gen[ipt+1] - pt_bin_edges_gen[ipt]
                        # pt_bin_width_up = pt_bin_edges_gen[ipt+2] - pt_bin_edges_gen[ipt+1]
                        # factor = (pt_bin_width_down + pt_bin_width_up)
                        # value_pt_down = bin_factors[bin_ind_pt_down]
                        # value_pt_up = bin_factors[bin_ind_pt_up]
                        # ref_hist = unreg_unfolder.unfolded

                        ref_hist = unreg_unfolder.hist_truth

                        val_down = ref_hist.GetBinContent(bin_ind_pt_down)
                        value_pt_down = 1./val_down if val_down != 0 else 0

                        val_up = ref_hist.GetBinContent(bin_ind_pt_up)
                        value_pt_up = 1./val_up if val_up != 0 else 0

                        # value_pt_down = bin_factors[bin_ind_pt_down]
                        # value_pt_up = bin_factors[bin_ind_pt_up]
                        value_pt_cen = - (value_pt_down + value_pt_up)
                        # print(bin_ind_pt_down, value_pt_down, bin_ind_cen, value_pt_cen, bin_ind_pt_up, value_pt_up)
                        L_args = [bin_ind_pt_down, value_pt_down, bin_ind_cen, value_pt_cen, bin_ind_pt_up, value_pt_up]
                        L_matrix_entries.append(L_args)
                        unfolder.AddRegularisationCondition(*L_args)

                        # value_pt_down = unfolded_max/ref_hist.GetBinContent(bin_ind_pt_down)
                        # value_pt_up = unfolded_max/ref_hist.GetBinContent(bin_ind_pt_up)
                        # value_var_down = unfolded_max/ref_hist.GetBinContent(bin_ind_var_down)
                        # value_var_up = unfolded_max/ref_hist.GetBinContent(bin_ind_var_up)
                        # value_cen = - (value_pt_down + value_pt_up + value_var_down + value_var_up)
                        # print(bin_ind_pt_down, value_pt_down, bin_ind_cen, value_cen, bin_ind_pt_up, value_pt_up)
                        # indices = [bin_ind_pt_down, bin_ind_var_down, bin_ind_cen, bin_ind_pt_up, bin_ind_var_up]
                        # row_data = [value_pt_down, value_var_down, value_cen, value_pt_up, value_var_up]
                        # unfolder.AddRegularisationCondition(5, array('i', indices), array('d', row_data))
                        # print(indices, row_data)

                # loop over existing regularisation conditions, since we want to modify them
                # in our main unfolder
                # for iy in range(1, orig_Lmatrix.GetNbinsY()+1):
                #     # Look for gen bin number where values start for this regularisation row
                #     left_bin, mid_bin, right_bin = 0, 0, 0
                #     left_bin_val, mid_bin_val, right_bin_val = 0, 0, 0
                #     for ix in range(1, orig_Lmatrix.GetNbinsX()+1):
                #         bin_content = orig_Lmatrix.GetBinContent(ix, iy)
                #         if bin_content != 0:
                #             if left_bin == 0:
                #                 left_bin = ix
                #                 left_bin_val = bin_content
                #                 continue
                #             elif mid_bin == 0:
                #                 mid_bin = ix
                #                 mid_bin_val = bin_content
                #                 continue
                #             else:
                #                 right_bin = ix
                #                 right_bin_val = bin_content
                #                 break # got em all

                    # Things to try:
                    # - Ignore the original reg. condition: just use bin_factors
                    # - Same, but with pt bin_width compensated
                    # - Modify existing reg.condition with bin_factors
                    # - Same, but compensate for pt bin width
                    # - Weighting my median relative error? e.g. care more about
                    # regions/bins with larger error bars

                    # Rescale accord to pT, and according to pT bin width
                    # since the original was divided by both pt bin width and lambda bin width
                    # Doesn't matter left or right for bin widht - only care about pt bin width
                    # pt_factor = bin_factors[mid_bin] * bin_widths[left_bin][1]
                    # pt_factor = bin_factors[mid_bin]

                    # - signs since RegularizeCurvature also adds in a - sign,
                    # and we want to match the original sign (which is +ve)
                    # scale_left = -pt_factor
                    # scale_right = -pt_factor
                    # scale_left = -left_bin_val * pt_factor
                    # scale_right = -right_bin_val * pt_factor

                    # print("Adding regularisation rule: nR=%d, gen bins: [%d - %d], factors: [%f, %f, %f]" % (iy, left_bin, right_bin, -scale_left, 2*(scale_left+scale_right), -scale_right) )
                    # unfolder.RegularizeCurvature(left_bin, mid_bin, right_bin, scale_left, scale_right)

            # exit()
            tau = 3.322366546731066e-10
            tau = 0
            scan_mode = ROOT.TUnfoldDensity.kEScanTauRhoAvgSys
            scan_mode = ROOT.TUnfoldDensity.kEScanTauRhoAvg
            scan_distribution = unfolder.distribution
            if REGULARIZE == "L":
                print("Regularizing with ScanLcurve, please be patient...")
                l_scanner = LCurveScanner()
                tau = l_scanner.scan_L(tunfolder=unfolder,
                                       n_scan=args.nScan,
                                       tau_min=region['tau_limits'][angle.var][0]/10,
                                       tau_max=region['tau_limits'][angle.var][1])
                                       # tau_min=tau/10,
                                       # tau_max=tau*10)
                print("Found tau:", tau)
                l_scanner.plot_scan_L_curve(output_filename="%s/scanL_%s.%s" % (this_output_dir, append, OUTPUT_FMT))
                l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

            elif REGULARIZE == "tau":
                print("Regularizing with ScanTau, please be patient...")
                tau_scanner = TauScanner()
                tau = tau_scanner.scan_tau(tunfolder=unfolder,
                                           n_scan=args.nScan,
                                           tau_min=region['tau_limits'][angle.var][0],
                                           tau_max=region['tau_limits'][angle.var][1],
                                           scan_mode=scan_mode,
                                           distribution=scan_distribution,
                                           axis_steering=unfolder.axisSteering)
                print("Found tau:", tau)
                tau_scanner.plot_scan_tau(output_filename="%s/scantau_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

            if REGULARIZE != "None":
                title = "L matrix %s %s region %s" % (jet_algo, region['label'], angle_str)
                unfolder_plotter.draw_L_matrix(title=title, **plot_args)
                title = "L^{T}L matrix, %s, %s region, %s" % (jet_algo, region['label'], angle_str)
                unfolder_plotter.draw_L_matrix_squared(title=title, **plot_args)
                title = "L * (x - bias vector) %s %s region %s" % (jet_algo, region['label'], angle_str)
                unfolder_plotter.draw_Lx_minus_bias(title=title, **plot_args)

            # Do unfolding!
            # ------------------------------------------------------------------
            unfolder.do_unfolding(tau)
            unfolder.get_output(hist_name="unfolded_1d")

            # Do lots of extra gubbins, like caching matrices,
            # creating unfolded hists with different levels of uncertianties,
            # ------------------------------------------------------------------
            unfolder._post_process()

            # Check result with numpy
            unfolder.do_numpy_comparison(output_dir=this_output_dir)
            # Calculate experimental uncertainty shifts using results from another unfolding
            # ------------------------------------------------------------------
            ref_tfile_exp = None
            if args.doExperimentalSystsFromFile is not None:
                print("Getting experimental systematics from another file...")
                angle_output_dir = "%s/%s/%s" % (args.doExperimentalSystsFromFile, region['name'], angle.var)
                this_pkl_filename = os.path.join(angle_output_dir, "unfolding_result.pkl")
                if not os.path.isfile(this_pkl_filename):
                    raise IOError("Cannot find systematics file, %s" % this_pkl_filename)

                exp_syst_region = unpickle_region(this_pkl_filename)

                reference_unfolder = exp_syst_region['unfolder']
                ref_unfolded = reference_unfolder.unfolded

                for exp_syst in reference_unfolder.exp_systs:
                    # For each systematic source, we figure out the
                    # relative shift compared to the original nominal result,
                    # for each normalised distribution (i.e. per pt bin).
                    # We then apply it to our new nominal result, and store it
                    # Note that this is different to taking the total shift,
                    # calculating its fractional diff, and then then applying it
                    # to the new nominal result?

                    this_exp_syst = ExpSystematic(label=exp_syst.label,
                                                  syst_map=exp_syst.syst_map,
                                                  # copy the old shift/shifted, although no longer relevant here
                                                  syst_shift=exp_syst.syst_shift,
                                                  syst_shifted=exp_syst.syst_shifted)
                    unfolder.exp_systs.append(this_exp_syst)

                    # Add these to to HistBinChopper for later, but it isn't used
                    # Just to bypass internal checks that it exists in its cached objects
                    # when e.g. get_pt_bin_normed_div_bin_width() called
                    unfolder.hist_bin_chopper.add_obj(exp_syst.syst_shifted_label, unfolder.unfolded)
                    unfolder.hist_bin_chopper.add_obj(exp_syst.syst_shift_label, unfolder.unfolded)
                    unfolder.hist_bin_chopper.add_obj(exp_syst.syst_ematrix_label, unfolder.unfolded)

                    bins = unfolder.pt_bin_edges_gen
                    for ibin, (bin_edge_low, bin_edge_high) in enumerate(zip(bins[:-1], bins[1:])):
                        # get old nominal shape
                        hbc_args = dict(ind=ibin, binning_scheme='generator')
                        ref_nominal_hist = reference_unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width("unfolded", **hbc_args)

                        # calc rel shift
                        rel_syst_shift_hist = reference_unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width(exp_syst.syst_shift_label, **hbc_args).Clone()
                        rel_syst_shift_hist.Divide(ref_nominal_hist)

                        # get new nominal shape
                        new_syst_shift = unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width("unfolded", **hbc_args).Clone()

                        # apply rel shift to that to get shifted normalised hists
                        new_syst_shift.Multiply(rel_syst_shift_hist)
                        new_hist = unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width("unfolded", **hbc_args).Clone()
                        new_hist.Add(new_syst_shift)

                        # store shifted and shift hists
                        unfolder.remove_error_bars(new_hist)
                        key = unfolder.hist_bin_chopper._generate_key(exp_syst.syst_shifted_label,
                                                                      ind=ibin,
                                                                      axis='pt',
                                                                      do_norm=True,
                                                                      do_div_bin_width=True,
                                                                      binning_scheme='generator')
                        unfolder.hist_bin_chopper._cache[key] = new_hist

                        unfolder.remove_error_bars(new_syst_shift)
                        key = unfolder.hist_bin_chopper._generate_key(exp_syst.syst_shift_label,
                                                                      ind=ibin,
                                                                      axis='pt',
                                                                      do_norm=True,
                                                                      do_div_bin_width=True,
                                                                      binning_scheme='generator')
                        unfolder.hist_bin_chopper._cache[key] = new_syst_shift

                        # calculate error matrix
                        syst_ematrix = cu.shift_to_covariance(new_syst_shift)
                        key = unfolder.hist_bin_chopper._generate_key(exp_syst.syst_ematrix_label,
                                                                      ind=ibin,
                                                                      axis='pt',
                                                                      do_norm=True,
                                                                      do_div_bin_width=True,
                                                                      binning_scheme='generator')
                        unfolder.hist_bin_chopper._cache[key] = syst_ematrix


                # update region info
                region['experimental_systematics'] = [syst_dict for syst_dict in exp_syst_region['experimental_systematics']
                                                      if syst_dict['label'] in reference_unfolder.get_all_exp_syst_labels()]

            else:
                unfolder.setup_normalised_experimental_systs()

            # Draw big 1D distributions
            # ------------------------------------------------------------------
            title = "%s\n%s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_unfolded_1d(output_dir=this_output_dir, append=append, title=title)

            # reco using detector binning
            unfolder_plotter.draw_detector_1d(do_reco_mc=True,
                                              do_reco_data=not MC_INPUT,
                                              output_dir=this_output_dir,
                                              append=append,
                                              title=title)

            # reco using gen binning
            unfolder_plotter.draw_generator_1d(do_reco_data=not MC_INPUT,
                                               do_reco_data_bg_sub=False,
                                               do_reco_bg=False,
                                               do_reco_mc=True,
                                               do_reco_mc_bg_sub=False,
                                               do_truth_mc=True,
                                               output_dir=this_output_dir,
                                               append=append,
                                               title=title)

            if SUBTRACT_FAKES:
                # same plot but with bg-subtracted reco (incl fakes)
                unfolder_plotter.draw_detector_1d(do_reco_data_bg_sub=not MC_INPUT,
                                                  do_reco_bg=SUBTRACT_FAKES,
                                                  do_reco_mc_bg_sub=True,
                                                  output_dir=this_output_dir,
                                                  append='bg_fakes_subtracted_%s' % append,
                                                  title=title)

                # same but with generator-binning
                unfolder_plotter.draw_generator_1d(do_reco_data=False,
                                                   do_reco_data_bg_sub=not MC_INPUT,
                                                   do_reco_bg=True,
                                                   do_reco_mc=False,
                                                   do_reco_mc_bg_sub=True,
                                                   do_truth_mc=True,
                                                   output_dir=this_output_dir,
                                                   append='bg_fakes_subtracted_%s' % append,
                                                   title=title)

            # Draw projections of response matrix vs 1D hist to check normalisation OK
            # Only makes sense if the same MC events go into matrix & 1D plot
            # ------------------------------------------------------------------
            if not MC_SPLIT:
                proj_reco = unfolder.response_map.ProjectionY("proj_reco_%s" % (append))

                proj_gen = unfolder.response_map.ProjectionX("proj_gen_%s" % (append))
                draw_projection_comparison(unfolder.hist_truth, proj_gen,
                                           title="%s\n%s region" % (jet_algo, region['label']),
                                           xtitle="%s, Generator binning" % (angle_str),
                                           output_filename="%s/projection_gen_%s.%s" % (this_output_dir, append, OUTPUT_FMT))

                # print("projection reco #bins:", proj_reco.GetNbinsX())
                # print("response map # bins x:", unfolder.response_map.GetNbinsX())
                # print("response map # bins y:", unfolder.response_map.GetNbinsY())
                if SUBTRACT_FAKES:
                    # print("reco bg subtracted #bins:", hist_mc_reco_bg_subtracted.GetNbinsX())
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

            if len(region['experimental_systematics']) > 0 and MC_INPUT:
                # Do a big absolute 1D plot for sanity
                syst_contributions = [
                    Contribution(unfolder.get_syst_shifted_hist(sdict['label'], unfolder.get_unfolded_with_ematrix_stat()),
                                 label=sdict['label'],
                                 line_color=sdict['colour'], line_style=2 if 'down' in sdict['label'].lower() else 1, line_width=1,
                                 marker_color=sdict['colour'], marker_size=0, marker_style=21,
                                 subplot=unfolder.get_unfolded_with_ematrix_stat())
                    for sdict in region['experimental_systematics']
                ]
                unfolder_plotter.draw_unfolded_1d(do_unfolded=True, do_gen=False,
                                                  output_dir=this_output_dir,
                                                  append='exp_systs_%s' % append,
                                                  title=title,
                                                  other_contributions=syst_contributions,
                                                  subplot_title='#splitline{Variation /}{nominal}')

            # Draw matrices
            # ------------------------------------------------------------------
            unfolder_plotter.plot_bias_vector(title=title, **plot_args)

            title = "Response matrix, %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_response_matrix(title=title, **plot_args)

            title = "Response matrix normed by detector p_{T} bin, %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_response_matrix_normed_by_detector_pt(title=title, **plot_args)

            title = ("#splitline{Probability matrix, %s, %s region, %s}{Condition number: #sigma_{max} / #sigma_{min} = %.3g / %.3g = %g}"
                        % (jet_algo, region['label'], angle_str, unfolder.sigma_max, unfolder.sigma_min, unfolder.condition_number))
            unfolder_plotter.draw_probability_matrix(title=title, **plot_args)

            title = "%s, %s region, %s"% (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_failed_reco(title=title, **plot_args)

            title = "Correlation matrix, %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_correlation_matrix(title=title, draw_values=True, **plot_args)
            unfolder_plotter.draw_correlation_matrix(title=title, draw_values=False, **plot_args)

            title = "Error matrix (input), %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_error_matrix_input(title=title, **plot_args)

            title = "Error matrix (statistical input + backgrounds), %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_error_matrix_stat(title=title, **plot_args)

            title = "Error matrix (stat. response matrix), %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_error_matrix_stat_response(title=title, **plot_args)

            if args.doExperimentalSysts:
                for syst_dict in region['experimental_systematics']:
                    title = "Error matrix (%s systematic), %s, %s region, %s" % (syst_dict['label'], jet_algo, region['label'], angle_str)
                    unfolder_plotter.draw_error_matrix_syst(syst_dict['label'], title=title, **plot_args)

            if REGULARIZE != "None":
                title = "Error matrix (tau), %s, %s region, %s" % (jet_algo, region['label'], angle_str)
                unfolder_plotter.draw_error_matrix_tau(title=title, **plot_args)

            title = "Error matrix (total), %s, %s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_error_matrix_total(title=title, **plot_args)

            # Do forward-folding to check unfolding
            # ------------------------------------------------------------------
            # Do it on the unfolded result
            title = "%s\n%s region, %s" % (jet_algo, region['label'], angle_str)
            unfolder_plotter.draw_unfolded_folded(title=title, **plot_args)

            if MC_INPUT:
                # Folded MC truth
                title = "%s\n%s region, %s" % (jet_algo, region['label'], angle_str)
                unfolder_plotter.draw_truth_folded(title=title, **plot_args)


            # Do some bottom-line tests
            # ------------------------------------------------------------------
            # smeared_chi2, smeared_ndf, smeared_p = unfolder.calculate_chi2(one_hist=unfolder.get_folded_mc_truth(),
            #                                                                other_hist=unfolder.hist_mc_reco_bg_subtracted,
            #                                                                cov_inv_matrix=unfolder.get_vyy_inv_ndarray(),
            #                                                                # cov_inv_matrix=unfolder.get_vyy_inv_no_bg_ndarray(),
            #                                                                detector_space=True,
            #                                                                ignore_underflow_bins=True,
            #                                                                debugging_dir=os.path.join(this_output_dir, 'smeared_chi2_debug'))
            # print('smeared chi2:', smeared_chi2, smeared_ndf, smeared_chi2/smeared_ndf, smeared_p)

            # unfolded_chi2, unfolded_ndf, unfolded_p = unfolder.calculate_chi2(one_hist=unfolder.unfolded,
            #                                                                   other_hist=unfolder.hist_truth,
            #                                                                   cov_inv_matrix=unfolder.get_vxx_inv_ndarray(),
            #                                                                   detector_space=False,
            #                                                                   ignore_underflow_bins=True,
            #                                                                   debugging_dir=os.path.join(this_output_dir, 'unfolded_chi2_debug'))
            # print('unfolded chi2:', unfolded_chi2, unfolded_ndf, unfolded_chi2/unfolded_ndf, unfolded_p)

            # ------------------------------------------------------------------
            # UNFOLDING WITH ALTERNATIVE RESPONSE MATRIX
            # ------------------------------------------------------------------
            alt_unfolder = None
            alt_hist_mc_gen = None  # mc truth of the alternate generator used to make reponse matrix
            alt_hist_mc_reco = None  # reco of the alternate generator used to make reponse matrix
            alt_hist_mc_reco_bg_subtracted = None
            alt_hist_mc_reco_bg_subtracted_gen_binning = None

            # Do this outside the if statement, since we might use it later in plotting e.g. for data
            if not isinstance(region['alt_mc_tfile'], ROOT.TFile) and os.path.isfile(region['alt_mc_tfile']):
                region['alt_mc_tfile'] = cu.open_root_file(region['alt_mc_tfile'])

                alt_hist_mc_gen = cu.get_from_tfile(region['alt_mc_tfile'], "%s/hist_%s_truth_all" % (region['dirname'], angle_shortname))
                alt_hist_mc_reco = cu.get_from_tfile(region['alt_mc_tfile'], "%s/hist_%s_reco_all" % (region['dirname'], angle_shortname))

                alt_scale = unfolder.hist_truth.Integral() / alt_hist_mc_gen.Integral()
                alt_hist_mc_gen.Scale(alt_scale)
                alt_hist_mc_reco.Scale(alt_scale)

                # fakes-subtracted version
                alt_hist_fakes = hist_fakes_reco_fraction.Clone("hist_fakes_alt")
                alt_hist_fakes.Multiply(alt_hist_mc_reco)
                alt_hist_mc_reco_bg_subtracted = alt_hist_mc_reco.Clone()
                alt_hist_mc_reco_bg_subtracted.Add(alt_hist_fakes, -1)

                # gen-binned versions of detector-level plots
                alt_hist_mc_reco_gen_binning = cu.get_from_tfile(region['alt_mc_tfile'], "%s/hist_%s_reco_gen_binning" % (region['dirname'], angle_shortname))
                alt_hist_mc_reco_gen_binning.Scale(alt_scale)
                alt_hist_fakes_gen_binning = hist_fakes_reco_fraction_gen_binning.Clone("hist_fakes_alt_gen_binning")
                alt_hist_fakes_gen_binning.Multiply(alt_hist_mc_reco_gen_binning)
                alt_hist_mc_reco_bg_subtracted_gen_binning = alt_hist_mc_reco_gen_binning.Clone()
                alt_hist_mc_reco_bg_subtracted_gen_binning.Add(alt_hist_fakes_gen_binning, -1)

                unfolder.hist_bin_chopper.add_obj('alt_hist_truth', alt_hist_mc_gen)

            if args.useAltResponse:
                print("*" * 80)
                print("*** Unfolding with alternate response matrix ***")
                print("*" * 80)

                hist_mc_gen_reco_map_alt = cu.get_from_tfile(region['alt_mc_tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                hist_mc_gen_reco_map_alt.Scale(unfolder.response_map.Integral() / hist_mc_gen_reco_map_alt.Integral())  # just for display purposes, doesn't affect result

                alt_unfolder = MyUnfolder(response_map=rm_large_rel_error_bins(hist_mc_gen_reco_map_alt),
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
                                          distribution=unfolder.distribution,
                                          axisSteering=unfolder.axisSteering)

                is_herwig = "Herwig" in region['alt_mc_label']
                is_pythia8 = region['alt_mc_label'] == "Pythia8"  # not MG+Pythia9
                if is_herwig or is_pythia8:
                    # SetEpsMatrix ensures rank properly calculated when inverting
                    # Needed if you get message "rank of matrix E 55 expect 170"
                    # And unfolded looks wacko
                    alt_unfolder.SetEpsMatrix(1E-18)

                alt_unfolder.SetEpsMatrix(eps_matrix)

                alt_unfolder_plotter = MyUnfolderPlotter(alt_unfolder, is_data=not MC_INPUT)
                alt_output_dir = this_output_dir+"/altResponse"
                alt_plot_args = dict(output_dir=alt_output_dir,
                                     append=append)

                # Only plot response matrix
                # --------------------------------------------------------------
                title = ("#splitline{Probability matrix, %s region, %s, %s}{Condition number: #sigma_{max} / #sigma_{min} = %.3g / %.3g = %g}"
                            % (region['label'], angle_str, region['alt_mc_label'], unfolder.sigma_max, unfolder.sigma_min, unfolder.condition_number))
                alt_unfolder_plotter.draw_probability_matrix(title=title, **alt_plot_args)

                title = "Response matrix, %s, %s region, %s, %s" % (jet_algo, region['label'], angle_str, region['alt_mc_label'])
                alt_unfolder_plotter.draw_response_matrix(title=title, **alt_plot_args)

                # Set what is to be unfolded - same as main unfolder
                # --------------------------------------------------------------
                alt_unfolder.set_input(input_hist=reco_1d,
                                       input_hist_gen_binning=reco_1d_gen_binning,
                                       hist_truth=unfolder.hist_truth.Clone(),
                                       hist_mc_reco=unfolder.hist_mc_reco.Clone(),
                                       hist_mc_reco_bg_subtracted=unfolder.hist_mc_reco_bg_subtracted.Clone(),
                                       hist_mc_reco_gen_binning=unfolder.hist_mc_reco_gen_binning.Clone(),
                                       hist_mc_reco_gen_binning_bg_subtracted=unfolder.hist_mc_reco_gen_binning_bg_subtracted.Clone(),
                                       bias_factor=args.biasFactor)

                # Subtract fakes (treat as background)
                # --------------------------------------------------------------
                if SUBTRACT_FAKES:
                    alt_unfolder.subtract_background(hist_fakes_reco, "fakes")

                # Do any regularization
                # --------------------------------------------------------------
                # Setup L matrix
                if REGULARIZE != "None":
                    gen_node = unfolder.generator_binning.FindNode('generatordistribution')
                    for ilambda in range(len(unfolder.variable_bin_edges_gen[:-1])):
                        for ipt in range(len(unfolder.pt_bin_edges_gen[:-3])):
                            pt_cen = unfolder.pt_bin_edges_gen[ipt+1] + 0.000001  # add a tiny bit to make sure we're in the bin properly (I can never remember if included or not)
                            # lambda_cen = unfolder.variable_bin_edges_gen[ilambda+1] + 0.000001  # add a tiny bit to make sure we're in the bin properly (I can never remember if included or not)
                            lambda_cen = unfolder.variable_bin_edges_gen[ilambda] + 0.000001  # add a tiny bit to make sure we're in the bin properly (I can never remember if included or not)

                            bin_ind_pt_down = gen_node.GetGlobalBinNumber(lambda_cen, unfolder.pt_bin_edges_gen[ipt] + 0.000001)
                            bin_ind_pt_up = gen_node.GetGlobalBinNumber(lambda_cen, unfolder.pt_bin_edges_gen[ipt+2] + 0.000001)
                            bin_ind_cen = gen_node.GetGlobalBinNumber(lambda_cen, pt_cen)

                            val_down = unfolder.hist_truth.GetBinContent(bin_ind_pt_down)
                            value_pt_down = 1./val_down if val_down != 0 else 0

                            val_up = unfolder.hist_truth.GetBinContent(bin_ind_pt_down)
                            value_pt_up = 1./val_up if val_up != 0 else 0
                            value_pt_cen = - (value_pt_down + value_pt_up)

                            alt_unfolder.AddRegularisationCondition(bin_ind_pt_down, value_pt_down, bin_ind_cen, value_pt_cen, bin_ind_pt_up, value_pt_up)

                # Scan for best regularisation strength
                alt_tau = 0
                if REGULARIZE == "L":
                    print("Regularizing alternative with ScanL, please be patient...")
                    alt_L_scanner = LCurveScanner()
                    alt_tau = alt_l_scanner.scan_L(tunfolder=alt_unfolder,
                                               n_scan=args.nScan,
                                               tau_min=region['tau_limits'][angle.var][0],
                                               tau_max=region['tau_limits'][angle.var][1])
                    print("Found tau:", alt_tau)
                    alt_l_scanner.plot_scan_L_curve(output_filename="%s/scanL_alt_%s.%s" % (alt_output_dir, unfolder.variable_name, OUTPUT_FMT))
                    alt_l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_alt_%s.%s" % (alt_output_dir, unfolder.variable_name, OUTPUT_FMT))

                elif REGULARIZE == "tau":
                    print("Regularizing alternative with ScanTau, please be patient...")
                    alt_tau_scanner = TauScanner()
                    alt_tau = alt_tau_scanner.scan_tau(tunfolder=alt_unfolder,
                                                       n_scan=args.nScan,
                                                       tau_min=region['tau_limits'][angle.var][0],
                                                       tau_max=region['tau_limits'][angle.var][1],
                                                       scan_mode=scan_mode,
                                                       distribution=scan_distribution,
                                                       axis_steering=alt_unfolder.axisSteering)
                    print("Found tau for alt matrix:", alt_tau)
                    alt_tau_scanner.plot_scan_tau(output_filename="%s/scantau_alt_%s.%s" % (alt_output_dir, alt_unfolder.variable_name, OUTPUT_FMT))

                if REGULARIZE != "None":
                    title = "L matrix, %s region, %s, alt. response (%s)" % (region['label'], angle_str, region['alt_mc_label'])
                    alt_unfolder_plotter.draw_L_matrix(title=title, **alt_plot_args)
                    title = "L^{T}L matrix, %s region, %s, alt. response (%s)" % (region['label'], angle_str, region['alt_mc_label'])
                    alt_unfolder_plotter.draw_L_matrix_squared(title=title, **alt_plot_args)
                    title = "L * (x - bias vector), %s region, %s,  alt. response (%s)" % (region['label'], angle_str, region['alt_mc_label'])
                    alt_unfolder_plotter.draw_Lx_minus_bias(title=title, **alt_plot_args)

                # Do unfolding!
                # --------------------------------------------------------------
                alt_unfolder.do_unfolding(alt_tau)
                alt_unfolder.get_output(hist_name="alt_unfolded_1d")
                alt_unfolder._post_process()

                if SUBTRACT_FAKES:
                    title = "%s\n%s region, %s, %s response map" % (jet_algo, region['label'], angle_str, region['alt_mc_label'])
                    alt_unfolder_plotter.draw_detector_1d(do_reco_data_bg_sub=not MC_INPUT,
                                                          do_reco_bg=SUBTRACT_FAKES,
                                                          do_reco_mc_bg_sub=True,
                                                          output_dir=alt_output_dir,
                                                          append='bg_fakes_subtracted_%s' % append,
                                                          title=title)

                    # same but with generator-binning
                    alt_unfolder_plotter.draw_generator_1d(do_reco_data=False,
                                                           do_reco_data_bg_sub=not MC_INPUT,
                                                           do_reco_bg=True,
                                                           do_reco_mc=False,
                                                           do_reco_mc_bg_sub=True,
                                                           do_truth_mc=True,
                                                           output_dir=alt_output_dir,
                                                           append='bg_fakes_subtracted_%s' % append,
                                                           title=title)

                alt_title = "%s\n%s region, %s, %s response map" % (jet_algo, region['label'], angle_str, region['alt_mc_label'])
                alt_unfolder_plotter.draw_unfolded_1d(title=alt_title, **alt_plot_args)

                title = "Correlation matrix, %s, %s region, %s, %s response map" % (jet_algo, region['label'], angle_str, region['alt_mc_label'])
                alt_unfolder_plotter.draw_correlation_matrix(title=title, draw_values=False, **alt_plot_args)

                title = "Error matrix (statistical input + backgrounds), %s, %s region, %s, %s response map" % (jet_algo, region['label'], angle_str, region['alt_mc_label'])
                alt_unfolder_plotter.draw_error_matrix_stat(title=title, **alt_plot_args)

                title = "Error matrix (total), %s, %s region, %s, %s response map" % (jet_algo, region['label'], angle_str, region['alt_mc_label'])
                alt_unfolder_plotter.draw_error_matrix_total(title=title, **alt_plot_args)

                # Do some more bottom-line tests:
                # --------------------------------------------------------------
                # Do before saving to file otherwise objects get deleted
                # if alt_unfolder:
                #     smeared_alt_chi2, smeared_alt_ndf, smeared_alt_p = unfolder.calculate_chi2(one_hist=alt_unfolder.get_folded_mc_truth(),
                #                                                                                other_hist=unfolder.hist_mc_reco_bg_subtracted,
                #                                                                                cov_inv_matrix=unfolder.get_vyy_inv_ndarray(),
                #                                                                                # cov_inv_matrix=unfolder.get_vyy_inv_no_bg_ndarray(),
                #                                                                                detector_space=True,
                #                                                                                ignore_underflow_bins=True,
                #                                                                                debugging_dir=os.path.join(this_output_dir, 'smeared_alt_chi2_debug'))
                #     print('smeared chi2 (alt MC):', smeared_alt_chi2, smeared_alt_ndf, smeared_alt_chi2/smeared_alt_ndf, smeared_alt_p)

                # print(unfolder.unfolded.Integral())
                # print(alt_hist_mc_gen.Integral())
                # if alt_hist_mc_gen:
                #     unfolded_alt_chi2, unfolded_alt_ndf, unfolded_alt_p = unfolder.calculate_chi2(one_hist=unfolder.unfolded,
                #                                                                                   other_hist=alt_hist_mc_gen,
                #                                                                                   cov_inv_matrix=unfolder.get_vxx_inv_ndarray(),
                #                                                                                   detector_space=False,
                #                                                                                   ignore_underflow_bins=True,
                #                                                                                   debugging_dir=os.path.join(this_output_dir, 'unfolded_alt_chi2_debug'))
                #     print('unfolded chi2 (alt MC):', unfolded_alt_chi2, unfolded_alt_ndf, unfolded_alt_chi2/unfolded_alt_ndf, unfolded_alt_p)

                alt_unfolder.setup_normalised_results()

                region['alt_unfolder'] = alt_unfolder

            # Bit gnarly - have to save this stuff manually
            region["alt_hist_mc_gen"] = alt_hist_mc_gen
            region["alt_hist_mc_reco"] = alt_hist_mc_reco
            region["alt_hist_mc_reco_bg_subtracted"] = alt_hist_mc_reco_bg_subtracted
            region["alt_hist_mc_reco_bg_subtracted_gen_binning"] = alt_hist_mc_reco_bg_subtracted_gen_binning

            # ------------------------------------------------------------------
            # SCALE SYST VARIATIONS
            # ------------------------------------------------------------------
            # For each scale variation, we unfold the nominal input using the
            # variation's response matrix. Then we take the envelope as the final
            # scale uncertainty
            if args.doScaleSysts:
                for ind, scale_dict in enumerate(region['scale_systematics']):
                    scale_label = scale_dict['label']
                    scale_label_no_spaces = cu.no_space_str(scale_dict['label'])

                    print("*" * 80)
                    print("*** Unfolding scale variation:", scale_label, "(%d/%d) ***" % (ind+1, len(region['scale_systematics'])))
                    print("*" * 80)

                    if isinstance(scale_dict['tfile'], str):
                        scale_dict['tfile'] = cu.open_root_file(scale_dict['tfile'])
                    scale_dict['response_map'] = cu.get_from_tfile(scale_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))

                    scale_unfolder = MyUnfolder(response_map=scale_dict['response_map'],
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
                                                distribution=unfolder.distribution,
                                                axisSteering=unfolder.axisSteering)

                    # Needed beacuse some of the variations struggle to unfold
                    # Even 1E-18 wouldn't work - needs to be v.small
                    scale_unfolder.SetEpsMatrix(eps_matrix)

                    scale_unfolder_plotter = MyUnfolderPlotter(scale_unfolder, is_data=not MC_INPUT)
                    scale_output_dir = this_output_dir+"/scaleSyst/"+scale_label_no_spaces
                    scale_plot_args = dict(output_dir=scale_output_dir,
                                           append=append)

                    # Set what is to be unfolded
                    # --------------------------------------------------------------
                    # Same input as nominal unfolder, since we only change responsematrix
                    scale_unfolder.set_input(input_hist=reco_1d,
                                           input_hist_gen_binning=reco_1d_gen_binning,
                                           hist_truth=hist_mc_gen,
                                           hist_mc_reco=hist_mc_reco,
                                           hist_mc_reco_bg_subtracted=hist_mc_reco_bg_subtracted,
                                           hist_mc_reco_gen_binning=hist_mc_reco_gen_binning,
                                           hist_mc_reco_gen_binning_bg_subtracted=hist_mc_reco_gen_binning_bg_subtracted,
                                           bias_factor=args.biasFactor)

                    # Subtract fakes (treat as background)
                    # --------------------------------------------------------------
                    if SUBTRACT_FAKES:
                        scale_unfolder.subtract_background(hist_fakes_reco, "Signal fakes", scale=1., scale_err=0.0)
                        scale_unfolder.subtract_background_gen_binning(hist_fakes_reco_gen_binning, "Signal fakes", scale=1., scale_err=0.0)

                    # Do any regularization
                    # --------------------------------------------------------------
                    scale_tau = 0
                    if REGULARIZE != "None":
                        raise RuntimeError("Robin hasn't implemented regularization for scale systs yet")

                    scale_dict['tau'] = scale_tau

                    # Do unfolding!
                    # --------------------------------------------------------------
                    scale_unfolder.do_unfolding(scale_tau)
                    scale_unfolder.get_output(hist_name="scale_%s_unfolded_1d" % (scale_label_no_spaces))
                    scale_unfolder._post_process()
                    scale_unfolder.setup_normalised_results()

                    # Plot absolute 1D result
                    scale_title = "%s\n%s region, %s, %s input" % (jet_algo, region['label'], angle_str, scale_label)
                    scale_unfolder_plotter.draw_unfolded_1d(title=scale_title, **scale_plot_args)
                    scale_dict['unfolder'] = scale_unfolder

                    del scale_dict['response_map']  # save memory
                    scale_dict['unfolder'] = scale_unfolder

                unfolder.create_normalised_scale_syst_uncertainty(region['scale_systematics'])
                unfolder.create_normalised_scale_syst_ematrices()

            # ------------------------------------------------------------------
            # LOAD SCALE VARIATIONS FROM ANOTHER FILE
            # ------------------------------------------------------------------
            # Load scale systs from another reference file, and calc fractional
            # uncertainty, apply to this unfolded result
            if args.doScaleSystsFromFile is not None:
                scale_dir = "%s/%s/%s" % (args.doScaleSystsFromFile, region['name'], angle.var)

                this_pkl_filename = os.path.join(scale_dir, "unfolding_result.pkl")
                if not os.path.isfile(this_pkl_filename):
                    raise IOError("Cannot find scale systematics file, %s" % this_pkl_filename)
                scale_syst_region = unpickle_region(this_pkl_filename)

                # update original region object with the scale syst info from the reference file
                region['scale_systematics'] = scale_syst_region['scale_systematics']

                # Add dummy object to hist_bin_chopper for later, so we can directly manipulate the cache
                # this replaces the procedure in create_normalised_scale_syst_uncertainty()
                unfolder.hist_bin_chopper.add_obj(unfolder.scale_uncert_name, unfolder.get_unfolded_with_ematrix_stat())
                unfolder.hist_bin_chopper.add_obj("unfolded_stat_err", unfolder.get_unfolded_with_ematrix_stat())

                for ibin_pt in range(len(unfolder.pt_bin_edges_gen[:-1])):
                    # Calculate scale uncertainty by taking relative uncertainty
                    # from reference file (which was already calculated in the 1st running),
                    # and applying to this result
                    key = unfolder.hist_bin_chopper._generate_key(unfolder.scale_uncert_name,
                                                                  ind=ibin_pt,
                                                                  axis='pt',
                                                                  do_norm=True,
                                                                  do_div_bin_width=True,
                                                                  binning_scheme='generator')
                    ref_scale_syst = scale_syst_region['unfolder'].hist_bin_chopper._cache[key]
                    scale_syst = unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width("unfolded_stat_err", ibin_pt, binning_scheme='generator').Clone()
                    for i in range(1, scale_syst.GetNbinsX()+1):
                        if ref_scale_syst.GetBinContent(i) != 0:
                            rel_err = ref_scale_syst.GetBinError(i) / ref_scale_syst.GetBinContent(i)
                        else:
                            rel_err = 0
                        scale_syst.SetBinError(i, rel_err * scale_syst.GetBinContent(i))
                    unfolder.hist_bin_chopper._cache[key] = scale_syst

                unfolder.create_normalised_scale_syst_ematrices()

            # ------------------------------------------------------------------
            # BIG ABSOLUTE PLOT WITH ALL SCALE VARIATIONS
            # ------------------------------------------------------------------
            if len(region['scale_systematics']) > 0:
                # Do a big absolute 1D plots for sanity
                scale_contributions = [
                    Contribution(mdict['unfolder'].get_unfolded_with_ematrix_stat(),
                                 label=mdict['label'],
                                 line_color=mdict['colour'], line_style=1, line_width=1,
                                 marker_color=mdict['colour'], marker_size=0, marker_style=21,
                                 subplot=unfolder.get_unfolded_with_ematrix_stat())
                    for mdict in region['scale_systematics']
                ]
                unfolder_plotter.draw_unfolded_1d(do_unfolded=True, do_gen=False,
                                                  output_dir=this_output_dir,
                                                  append='scale_systs_%s' % append,
                                                  title=title,
                                                  other_contributions=scale_contributions,
                                                  subplot_title='#splitline{Variation /}{nominal}')

            # ------------------------------------------------------------------
            # MODEL INPUT VARIATIONS
            # ------------------------------------------------------------------
            # For each model variation, we unfold using the same settings as
            # the nominal one, just changing the input 1D hist
            if args.doModelSysts:
                syst_entries = []
                for ind, syst_dict in enumerate(region['model_systematics']):
                    syst_label = syst_dict['label']
                    syst_label_no_spaces = cu.no_space_str(syst_dict['label'])

                    print("*" * 80)
                    print("*** Unfolding with alternate input:", syst_label, "(%d/%d) ***" % (ind+1, len(region['model_systematics'])))
                    print("*" * 80)

                    is_herwig = "Herwig" in syst_label

                    mc_hname_append = "split" if MC_SPLIT else "all"
                    if is_herwig:
                        # use all the stats!
                        mc_hname_append = "all"
                    if not isinstance(syst_dict['tfile'], ROOT.TFile):
                        syst_dict['tfile'] = cu.open_root_file(syst_dict['tfile'])
                    hist_syst_reco = cu.get_from_tfile(syst_dict['tfile'], "%s/hist_%s_reco_%s" % (region['dirname'], angle_shortname, mc_hname_append))
                    hist_syst_gen = cu.get_from_tfile(syst_dict['tfile'], "%s/hist_%s_truth_%s" % (region['dirname'], angle_shortname, mc_hname_append))

                    syst_unfolder = MyUnfolder(response_map=unfolder.response_map,
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
                                               distribution=unfolder.distribution,
                                               axisSteering=unfolder.axisSteering)

                    syst_output_dir = this_output_dir+"/modelSyst_"+syst_label_no_spaces
                    syst_unfolder_plotter = MyUnfolderPlotter(syst_unfolder, is_data=False)
                    syst_plot_args = dict(output_dir=syst_output_dir,
                                          append=append)

                    if is_herwig:
                        # SetEpsMatrix ensures rank properly calculated when inverting
                        # Needed if you get message "rank of matrix E 55 expect 170"
                        # And unfolded looks wacko
                        syst_unfolder.SetEpsMatrix(1E-18)

                    syst_unfolder.SetEpsMatrix(eps_matrix)

                    # because we only care about shape, not overall normalisation
                    # (which can artificially increase/decrease errors)
                    # we normalise to the nominal integral
                    # Note that we use the scaling from gen level, to take
                    # into account any reco-dependent efficiencies
                    # TODO: is this right?
                    sf = hist_mc_gen.Integral() / hist_syst_gen.Integral()
                    hist_syst_reco.Scale(sf)
                    hist_syst_gen.Scale(sf)

                    if SUBTRACT_FAKES:
                        # Use the background template from the nominal MC
                        # (since we're only testing different input shapes,
                        # and our bkg estimate is always from MC)
                        hist_fakes_syst = hist_fakes_reco_fraction.Clone("hist_fakes_syst_%s" % syst_label_no_spaces)
                        hist_fakes_syst.Multiply(hist_syst_reco)

                    hist_mc_reco_bg_subtracted = hist_syst_reco.Clone()
                    hist_mc_reco_bg_subtracted.Add(hist_fakes_syst, -1)

                    # Set what is to be unfolded
                    # --------------------------------------------------------------
                    syst_unfolder.set_input(input_hist=hist_syst_reco,
                                            hist_truth=hist_syst_gen,
                                            hist_mc_reco=hist_syst_reco,
                                            hist_mc_reco_bg_subtracted=hist_mc_reco_bg_subtracted,
                                            bias_factor=args.biasFactor)

                    # Subtract fakes (treat as background)
                    # --------------------------------------------------------------
                    if SUBTRACT_FAKES:
                        syst_unfolder.subtract_background(hist_fakes_syst, "fakes")

                    # also show nominal bg-subtracted input for comparison
                    ocs = [
                        Contribution(unfolder.input_hist_bg_subtracted,
                                     label='Nominal unfolding input (bg-subtracted)',
                                     line_color=ROOT.kRed, line_width=1)
                    ]
                    syst_unfolder_plotter.draw_detector_1d(do_reco_data=False,
                                                           do_reco_data_bg_sub=False,
                                                           do_reco_bg=True,
                                                           do_reco_mc=False,
                                                           do_reco_mc_bg_sub=True,
                                                           other_contributions=ocs,
                                                           output_dir=syst_plot_args['output_dir'],
                                                           append='bg_fakes_subtracted_%s' % append,
                                                           title="%s region, %s, %s" % (region['label'], angle_str, syst_label))

                    # Do any regularization
                    # --------------------------------------------------------------
                    # Setup L matrix
                    if REGULARIZE != "None":
                        # gen_node = unfolder.generator_binning.FindNode('generatordistribution')
                        # for ilambda in range(len(unfolder.variable_bin_edges_gen[:-1])):
                        #     for ipt in range(len(unfolder.pt_bin_edges_gen[:-3])):
                        #         pt_cen = unfolder.pt_bin_edges_gen[ipt+1] + 0.000001  # add a tiny bit to make sure we're in the bin properly (I can never remember if included or not)
                        #         # lambda_cen = unfolder.variable_bin_edges_gen[ilambda+1] + 0.000001  # add a tiny bit to make sure we're in the bin properly (I can never remember if included or not)
                        #         lambda_cen = unfolder.variable_bin_edges_gen[ilambda] + 0.000001  # add a tiny bit to make sure we're in the bin properly (I can never remember if included or not)

                        #         bin_ind_pt_down = gen_node.GetGlobalBinNumber(lambda_cen, unfolder.pt_bin_edges_gen[ipt] + 0.000001)
                        #         bin_ind_pt_up = gen_node.GetGlobalBinNumber(lambda_cen, unfolder.pt_bin_edges_gen[ipt+2] + 0.000001)
                        #         bin_ind_cen = gen_node.GetGlobalBinNumber(lambda_cen, pt_cen)

                        #         val_down = unfolder.hist_truth.GetBinContent(bin_ind_pt_down)
                        #         value_pt_down = 1./val_down if val_down != 0 else 0

                        #         val_up = unfolder.hist_truth.GetBinContent(bin_ind_pt_down)
                        #         value_pt_up = 1./val_up if val_up != 0 else 0
                        #         value_pt_cen = - (value_pt_down + value_pt_up)

                                # syst_unfolder.AddRegularisationCondition(bin_ind_pt_down, value_pt_down, bin_ind_cen, value_pt_cen, bin_ind_pt_up, value_pt_up)
                        for L_args in L_matrix_entries:
                            syst_unfolder.AddRegularisationCondition(*L_args)

                    # Scan for best regularisation strength
                    syst_tau = 0
                    if REGULARIZE == "L":
                        print("Regularizing systematic model with ScanL, please be patient...")
                        syst_l_scanner = LCurveScanner()
                        syst_tau = syst_l_scanner.scan_L(tunfolder=syst_unfolder,
                                                         n_scan=args.nScan,
                                                         tau_min=region['tau_limits'][angle.var][0],
                                                         tau_max=region['tau_limits'][angle.var][1])
                        print("Found tau:", syst_tau)
                        syst_l_scanner.plot_scan_L_curve(output_filename="%s/scanL_syst_%s_%s.%s" % (syst_output_dir, syst_label_no_spaces, unfolder.variable_name, OUTPUT_FMT))
                        syst_l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_syst_%s_%s.%s" % (syst_output_dir, syst_label_no_spaces, unfolder.variable_name, OUTPUT_FMT))

                    elif REGULARIZE == "tau":
                        print("Regularizing systematic model with ScanTau, please be patient...")
                        syst_tau_scanner = TauScanner()
                        syst_tau = syst_tau_scanner.scan_tau(tunfolder=syst_unfolder,
                                                             n_scan=args.nScan,
                                                             tau_min=region['tau_limits'][angle.var][0],
                                                             tau_max=region['tau_limits'][angle.var][1],
                                                             scan_mode=scan_mode,
                                                             distribution=scan_distribution,
                                                             axis_steering=syst_unfolder.axisSteering)
                        print("Found tau for syst matrix:", syst_tau)
                        syst_tau_scanner.plot_scan_tau(output_filename="%s/scantau_syst_%s_%s.%s" % (syst_output_dir, syst_label_no_spaces, syst_unfolder.variable_name, OUTPUT_FMT))

                    region['model_systematics'][ind]['tau'] = syst_tau

                    if REGULARIZE != "None":
                        title = "L matrix, %s region, %s,\n%s" % (region['label'], angle_str, syst_label)
                        syst_unfolder_plotter.draw_L_matrix(title=title, **syst_plot_args)
                        title = "L^{T}L matrix, %s region, %s,\n%s" % (region['label'], angle_str, syst_label)
                        syst_unfolder_plotter.draw_L_matrix_squared(title=title, **syst_plot_args)
                        title = "L * (x - bias vector), %s region, %s,\n%s" % (region['label'], angle_str, syst_label)
                        syst_unfolder_plotter.draw_Lx_minus_bias(title=title, **syst_plot_args)

                    # Do unfolding!
                    # --------------------------------------------------------------
                    syst_unfolder.do_unfolding(syst_tau)
                    syst_unfolder.get_output(hist_name="syst_%s_unfolded_1d" % (syst_label_no_spaces))
                    syst_unfolder._post_process()
                    syst_unfolder.setup_normalised_results()

                    syst_title = "%s\n%s region, %s, %s input" % (jet_algo, region['label'], angle_str, syst_label)
                    syst_unfolder_plotter.draw_unfolded_1d(title=syst_title, **syst_plot_args)
                    region['model_systematics'][ind]['unfolder'] = syst_unfolder

                    # Do 1D plot of nominal vs syst unfolded
                    # --------------------------------------------------------------
                    entries = []
                    # add nominal
                    label = 'MC' if MC_INPUT else "data"
                    entries.append(
                        Contribution(unfolder.unfolded, label="Unfolded (#tau = %.3g)" % (unfolder.tau),
                                     line_color=ROOT.kRed, line_width=1,
                                     marker_color=ROOT.kRed, marker_size=0.6, marker_style=20,
                                     subplot_line_color=ROOT.kRed, subplot_line_width=1,
                                     subplot_marker_color=ROOT.kRed, subplot_marker_size=0, subplot_marker_style=20,
                                     normalise_hist=False, subplot=unfolder.hist_truth),
                    )

                    entries.append(
                        Contribution(unfolder.hist_truth, label="Generator",
                                     line_color=ROOT.kBlue, line_width=1,
                                     marker_color=ROOT.kBlue, marker_size=0,
                                     normalise_hist=False),
                    )
                    # add systematic
                    entries.append(
                        Contribution(syst_unfolder.unfolded, label="Unfolded %s (#tau = %.3g)" % (syst_label, syst_unfolder.tau),
                                     line_color=syst_dict['colour'], line_width=1,
                                     marker_color=syst_dict['colour'], marker_size=0.6, marker_style=20+ind+1,
                                     subplot_line_color=syst_dict['colour'], subplot_line_width=1,
                                     subplot_marker_color=syst_dict['colour'], subplot_marker_size=0, subplot_marker_style=20,
                                     normalise_hist=False, subplot=syst_unfolder.hist_truth),
                    )
                    syst_entries.append(entries[-1])

                    entries.append(
                        Contribution(syst_unfolder.hist_truth, label="Generator (%s)" % (syst_label),
                                     line_color=syst_dict['colour']+2, line_width=1, line_style=1,
                                     marker_color=syst_dict['colour']+2, marker_size=0,
                                     normalise_hist=False),
                    )
                    syst_entries.append(entries[-1])

                    title = "%s\n%s region, %s, %s input" % (jet_algo, region['label'], angle_str, syst_label)
                    plot = Plot(entries,
                                what='hist',
                                title=title,
                                xtitle="Generator bin number",
                                ytitle="N",
                                subplot_type='ratio',
                                subplot_title='Unfolded / Gen',
                                subplot_limits=(0, 2),
                                has_data=not MC_INPUT)
                    plot.default_canvas_size = (800, 600)
                    plot.plot("NOSTACK HISTE")
                    plot.set_logy(do_more_labels=False, override_check=True)
                    ymax = max([o.GetMaximum() for o in plot.contributions_objs])
                    plot.container.SetMaximum(ymax * 200)
                    ymin = max([o.GetMinimum(1E-10) for o in plot.contributions_objs])
                    plot.container.SetMinimum(ymin*0.01)
                    l, t = syst_unfolder_plotter.draw_pt_binning_lines(plot, which='gen', axis='x',
                                                                       do_underflow=True,
                                                                       do_labels_inside=True,
                                                                       do_labels_outside=False,
                                                                       labels_inside_align='lower'
                                                                       )
                    plot.legend.SetY1NDC(0.77)
                    plot.legend.SetY2NDC(0.88)
                    plot.legend.SetX1NDC(0.65)
                    plot.legend.SetX2NDC(0.88)
                    output_filename = "%s/unfolded_1d_modelSyst_%s.%s" % (syst_output_dir, syst_label_no_spaces, syst_unfolder_plotter.output_fmt)
                    plot.save(output_filename)

            # Load model (scale) systs from another reference file, and calc fractional
            # uncertainty, apply to this unfolded result
            if args.doModelSystsFromFile is not None:
                model_dir = "%s/%s/%s" % (args.doModelSystsFromFile, region['name'], angle.var)

                this_pkl_filename = os.path.join(model_dir, "unfolding_result.pkl")
                if not os.path.isfile(this_pkl_filename):
                    raise IOError("Cannot find model systematics file, %s" % this_pkl_filename)
                model_syst_region = unpickle_region(this_pkl_filename)

                # update original region object with the model syst info from the reference file
                region['model_systematics'] = model_syst_region['model_systematics']

                # Add dummy object to hist_bin_chopper for later, so we can directly manipulate the cache
                # this replaces the procedure in create_normalised_scale_syst_uncertainty()
                unfolder.hist_bin_chopper.add_obj(unfolder.scale_uncert_name, unfolder.get_unfolded_with_ematrix_stat())
                unfolder.hist_bin_chopper.add_obj("unfolded_stat_err", unfolder.get_unfolded_with_ematrix_stat())

                for ibin_pt in range(len(unfolder.pt_bin_edges_gen[:-1])):
                    # Calculate scale uncertainty by taking relative uncertainty
                    # from reference file (which was already calculated in the 1st running),
                    # and applying to this result
                    key = unfolder.hist_bin_chopper._generate_key(unfolder.scale_uncert_name,
                                                                  ind=ibin_pt,
                                                                  axis='pt',
                                                                  do_norm=True,
                                                                  do_div_bin_width=True,
                                                                  binning_scheme='generator')
                    ref_scale_syst = model_syst_region['unfolder'].hist_bin_chopper._cache[key]
                    scale_syst = unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width("unfolded_stat_err", ibin_pt, binning_scheme='generator').Clone()
                    for i in range(1, scale_syst.GetNbinsX()+1):
                        if ref_scale_syst.GetBinContent(i) != 0:
                            rel_err = ref_scale_syst.GetBinError(i) / ref_scale_syst.GetBinContent(i)
                        else:
                            rel_err = 0
                        scale_syst.SetBinError(i, rel_err * scale_syst.GetBinContent(i))
                    unfolder.hist_bin_chopper._cache[key] = scale_syst

                unfolder.create_normalised_scale_syst_ematrices()

            if len(region['model_systematics']) > 0 and MC_INPUT:
                # Do a big absolute 1D plots for sanity
                model_contributions = [
                    Contribution(mdict['unfolder'].input_hist_bg_subtracted,
                                 label=mdict['label'],
                                 line_color=mdict['colour'], line_style=1, line_width=1,
                                 marker_color=mdict['colour'], marker_size=0, marker_style=21,
                                 subplot=unfolder.input_hist_bg_subtracted)
                    for mdict in region['model_systematics']
                ]
                unfolder_plotter.draw_detector_1d(do_reco_mc_bg_sub=True,
                                                  output_dir=this_output_dir,
                                                  append='bg_fakes_subtracted_model_systs_%s' % append,
                                                  title=title,
                                                  other_contributions=model_contributions)
                model_contributions = [
                    Contribution(mdict['unfolder'].get_unfolded_with_ematrix_stat(),
                                 label=mdict['label'],
                                 line_color=mdict['colour'], line_style=1, line_width=1,
                                 marker_color=mdict['colour'], marker_size=0, marker_style=21,
                                 subplot=unfolder.get_unfolded_with_ematrix_stat())
                    for mdict in region['model_systematics']
                ]
                unfolder_plotter.draw_unfolded_1d(do_unfolded=True, do_gen=False,
                                                  output_dir=this_output_dir,
                                                  append='model_systs_%s' % append,
                                                  title=title,
                                                  other_contributions=model_contributions,
                                                  subplot_title='#splitline{Variation /}{nominal}')


            # ------------------------------------------------------------------
            # DO PDF VARIATIONS
            # ------------------------------------------------------------------
            if args.doPDFSysts:
                # first construct all new systematic variations dicts
                original_pdf_dict = region['pdf_systematics'][0]
                original_pdf_dict['label'] = '_PDF_template'  # initial _ to ignore it later on
                tfile = original_pdf_dict['tfile']
                if not isinstance(tfile, ROOT.TFile):
                    tfile = cu.open_root_file(tfile)

                region['pdf_systematics']  = []
                for pdf_ind in original_pdf_dict['variations']:
                    region['pdf_systematics'].append(
                        {
                            "label": "PDF_%d" % (pdf_ind),
                            "response_map": cu.get_from_tfile(tfile, "%s/tu_%s_GenReco_all_PDF_%d" % (region['dirname'], angle_shortname, pdf_ind)),
                            "colour": ROOT.kCyan+2,
                        })

                # Now run over all variations, unfolding the nominal inputs
                # but with the various response matrices
                for ind, pdf_dict in enumerate(region['pdf_systematics']):
                    pdf_label = pdf_dict['label']
                    pdf_label_no_spaces = cu.no_space_str(pdf_label)

                    print("*" * 80)
                    print("*** Unfolding with alternate input:", pdf_label, "(%d/%d) ***" % (ind+1, len(region['pdf_systematics'])))
                    print("*" * 80)

                    pdf_unfolder = MyUnfolder(response_map=pdf_dict['response_map'],
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
                                              distribution=unfolder.distribution,
                                              axisSteering=unfolder.axisSteering)

                    # Needed beacuse some fo the PDF variations struggle to unfold
                    # Even 1E-18 wouldn't work - needs to be v.small
                    pdf_unfolder.SetEpsMatrix(eps_matrix)

                    pdf_unfolder_plotter = MyUnfolderPlotter(pdf_unfolder, is_data=not MC_INPUT)
                    pdf_output_dir = this_output_dir+"/pdfSyst/"+pdf_label_no_spaces
                    pdf_plot_args = dict(output_dir=pdf_output_dir,
                                         append=append)

                    # Set what is to be unfolded
                    # --------------------------------------------------------------
                    # Same input as nominal unfolder, since we only change responsematrix
                    pdf_unfolder.set_input(input_hist=reco_1d,
                                           input_hist_gen_binning=reco_1d_gen_binning,
                                           hist_truth=hist_mc_gen,
                                           hist_mc_reco=hist_mc_reco,
                                           hist_mc_reco_bg_subtracted=hist_mc_reco_bg_subtracted,
                                           hist_mc_reco_gen_binning=hist_mc_reco_gen_binning,
                                           hist_mc_reco_gen_binning_bg_subtracted=hist_mc_reco_gen_binning_bg_subtracted,
                                           bias_factor=args.biasFactor)

                    # Subtract fakes (treat as background)
                    # --------------------------------------------------------------
                    if SUBTRACT_FAKES:
                        pdf_unfolder.subtract_background(hist_fakes_reco, "Signal fakes", scale=1., scale_err=0.0)
                        pdf_unfolder.subtract_background_gen_binning(hist_fakes_reco_gen_binning, "Signal fakes", scale=1., scale_err=0.0)

                    # Do any regularization
                    # --------------------------------------------------------------
                    # Setup L matrix
                    if REGULARIZE != "None":
                        # gen_node = unfolder.generator_binning.FindNode('generatordistribution')
                        # for ilambda in range(len(unfolder.variable_bin_edges_gen[:-1])):
                        #     for ipt in range(len(unfolder.pt_bin_edges_gen[:-3])):
                        #         pt_cen = unfolder.pt_bin_edges_gen[ipt+1] + 0.000001  # add a tiny bit to make sure we're in the bin properly (I can never remember if included or not)
                        #         # lambda_cen = unfolder.variable_bin_edges_gen[ilambda+1] + 0.000001  # add a tiny bit to make sure we're in the bin properly (I can never remember if included or not)
                        #         lambda_cen = unfolder.variable_bin_edges_gen[ilambda] + 0.000001  # add a tiny bit to make sure we're in the bin properly (I can never remember if included or not)

                        #         bin_ind_pt_down = gen_node.GetGlobalBinNumber(lambda_cen, unfolder.pt_bin_edges_gen[ipt] + 0.000001)
                        #         bin_ind_pt_up = gen_node.GetGlobalBinNumber(lambda_cen, unfolder.pt_bin_edges_gen[ipt+2] + 0.000001)
                        #         bin_ind_cen = gen_node.GetGlobalBinNumber(lambda_cen, pt_cen)

                        #         val_down = unfolder.hist_truth.GetBinContent(bin_ind_pt_down)
                        #         value_pt_down = 1./val_down if val_down != 0 else 0

                        #         val_up = unfolder.hist_truth.GetBinContent(bin_ind_pt_down)
                        #         value_pt_up = 1./val_up if val_up != 0 else 0
                        #         value_pt_cen = - (value_pt_down + value_pt_up)

                        #         pdf_unfolder.AddRegularisationCondition(bin_ind_pt_down, value_pt_down, bin_ind_cen, value_pt_cen, bin_ind_pt_up, value_pt_up)

                        for L_args in L_matrix_entries:
                            syst_unfolder.AddRegularisationCondition(*L_args)


                    pdf_tau = 0
                    if REGULARIZE == "L":
                        print("Regularizing PDF systematic with ScanL, please be patient...")
                        pdf_l_scanner = LCurveScanner()
                        pdf_tau = pdf_l_scanner.scan_L(tunfolder=pdf_unfolder,
                                                       n_scan=args.nScan,
                                                       tau_min=region['tau_limits'][angle.var][0],
                                                       tau_max=region['tau_limits'][angle.var][1])
                        print("Found tau:", pdf_tau)
                        pdf_l_scanner.plot_scan_L_curve(output_filename="%s/scanL_syst_%s_%s.%s" % (pdf_output_dir, pdf_label_no_spaces, unfolder.variable_name, OUTPUT_FMT))
                        pdf_l_scanner.plot_scan_L_curvature(output_filename="%s/scanLcurvature_syst_%s_%s.%s" % (pdf_output_dir, pdf_label_no_spaces, unfolder.variable_name, OUTPUT_FMT))

                    elif REGULARIZE == "tau":
                        print("Regularizing PDF systematic with ScanTau, please be patient...")
                        pdf_tau_scanner = TauScanner()
                        pdf_tau = pdf_tau_scanner.scan_tau(tunfolder=pdf_unfolder,
                                                             n_scan=args.nScan,
                                                             tau_min=region['tau_limits'][angle.var][0],
                                                             tau_max=region['tau_limits'][angle.var][1],
                                                             scan_mode=scan_mode,
                                                             distribution=scan_distribution,
                                                             axis_steering=pdf_unfolder.axisSteering)
                        print("Found tau for syst matrix:", pdf_tau)
                        pdf_tau_scanner.plot_scan_tau(output_filename="%s/scantau_syst_%s_%s.%s" % (pdf_output_dir, pdf_label_no_spaces, pdf_unfolder.variable_name, OUTPUT_FMT))

                    region['pdf_systematics'][ind]['tau'] = pdf_tau

                    # Do unfolding!
                    # --------------------------------------------------------------
                    pdf_unfolder.do_unfolding(pdf_tau)
                    pdf_unfolder.get_output(hist_name="syst_%s_unfolded_1d" % (pdf_label_no_spaces))
                    pdf_unfolder._post_process()

                    pdf_unfolder.setup_normalised_results()

                    pdf_title = "%s\n%s region, %s\n%s response matrix" % (jet_algo, region['label'], angle_str, pdf_label)
                    pdf_unfolder_plotter.draw_unfolded_1d(title=pdf_title, **pdf_plot_args)

                    del region['pdf_systematics'][ind]['response_map']  # save memory
                    region['pdf_systematics'][ind]['unfolder'] = pdf_unfolder

                unfolder.create_normalised_pdf_syst_uncertainty(region['pdf_systematics'])
                unfolder.create_normalised_pdf_syst_ematrices()

            # Load PDF syst from another reference file, and calc fractional
            # uncertainty, apply to this unfolded result
            if args.doPDFSystsFromFile is not None:
                pdf_dir = "%s/%s/%s" % (args.doPDFSystsFromFile, region['name'], angle.var)

                this_pkl_filename = os.path.join(pdf_dir, "unfolding_result.pkl")
                if not os.path.isfile(this_pkl_filename):
                    raise IOError("Cannot find PDF systematics file, %s" % this_pkl_filename)
                pdf_syst_region = unpickle_region(this_pkl_filename)

                # update original region object from reference file
                region['pdf_systematics'] = pdf_syst_region['pdf_systematics']

                # Add dummy object to hist_bin_chopper for later, so we can directly manipulate the cache
                unfolder.hist_bin_chopper.add_obj(unfolder.pdf_uncert_name, unfolder.get_unfolded_with_ematrix_stat())
                unfolder.hist_bin_chopper.add_obj("unfolded_stat_err", unfolder.get_unfolded_with_ematrix_stat())

                for ibin_pt in range(len(unfolder.pt_bin_edges_gen[:-1])):
                    # Calculate PDF uncertainty by taking relative uncertainty
                    # from reference file (which has already been calculated in the 1st running),
                    # and applying to this result
                    key = unfolder.hist_bin_chopper._generate_key(unfolder.pdf_uncert_name,
                                                                  ind=ibin_pt,
                                                                  axis='pt',
                                                                  do_norm=True,
                                                                  do_div_bin_width=True,
                                                                  binning_scheme='generator')
                    ref_pdf_syst = pdf_syst_region['unfolder'].hist_bin_chopper._cache[key]
                    pdf_syst = unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width("unfolded_stat_err", ibin_pt, binning_scheme='generator').Clone()
                    for i in range(1, pdf_syst.GetNbinsX()+1):
                        if ref_pdf_syst.GetBinContent(i) != 0:
                            rel_err = ref_pdf_syst.GetBinError(i) / ref_pdf_syst.GetBinContent(i)
                        else:
                            rel_err = 0
                        pdf_syst.SetBinError(i, rel_err * pdf_syst.GetBinContent(i))
                    unfolder.hist_bin_chopper._cache[key] = pdf_syst

                unfolder.create_normalised_pdf_syst_ematrices()


            if len(region['pdf_systematics']) > 0 and MC_INPUT:
                # Do a big absolute 1D plot for sanity
                pdf_contributions = [
                    Contribution(pdict['unfolder'].input_hist_bg_subtracted,
                                 label=pdict['label'],
                                 line_color=pdict['colour'], line_style=1, line_width=1,
                                 marker_color=pdict['colour'], marker_size=0, marker_style=21,
                                 subplot=unfolder.input_hist_bg_subtracted)
                    for pdict in region['pdf_systematics']
                ]
                unfolder_plotter.draw_detector_1d(do_reco_mc_bg_sub=True,
                                                  output_dir=this_output_dir,
                                                  append='bg_fakes_subtracted_pdf_systs_%s' % append,
                                                  title=title,
                                                  other_contributions=pdf_contributions)

                pdf_contributions = [
                    Contribution(pdict['unfolder'].get_unfolded_with_ematrix_stat(),
                                 label=pdict['label'],
                                 line_color=pdict['colour'], line_style=1, line_width=1,
                                 marker_color=pdict['colour'], marker_size=0, marker_style=21,
                                 subplot=unfolder.get_unfolded_with_ematrix_stat())
                    for pdict in region['pdf_systematics']
                ]
                unfolder_plotter.draw_unfolded_1d(do_unfolded=True, do_gen=False,
                                                  output_dir=this_output_dir,
                                                  append='pdf_systs_%s' % append,
                                                  title=title,
                                                  other_contributions=pdf_contributions,
                                                  subplot_title='#splitline{Variation /}{nominal}')

            # ------------------------------------------------------------------
            # Finally update normalised results
            # ------------------------------------------------------------------
            unfolder.setup_normalised_results()

            region['unfolder'] = unfolder

            # ------------------------------------------------------------------
            # Save everything to pickle / TFile
            # ------------------------------------------------------------------
            print("")
            print("unfolder attr sizes:")
            print("-"*80)
            cu.print_dict_item_sizes(unfolder.__dict__, recursive=True)
            print("-"*80)

            pickle_filename = os.path.join("%s/unfolding_result.pkl" % (this_output_dir))
            pickle_region(region, pickle_filename, infos=True, convert_tfile_to_str=True)
            print(">> Saved to pickle file", pickle_filename)
            print("")

            print(">> Saving unfolder to ROOT file")
            # print(unfolder.hist_bin_chopper.objects)
            # print(unfolder.hist_bin_chopper._cache)
            unfolder.save_unfolded_binned_hists_to_tfile(this_slim_tdir)

            # test the pickle file by un-pickling it
            print("Testing pickled file...")
            data = unpickle_region(pickle_filename)
            print("")
            print("...unpickled data:")
            print("    ", data)
            print("")
            # print("...data['unfolder'].hist_bin_chopper.objects:")
            # print("    ", data['unfolder'].hist_bin_chopper.objects)
            # print("")
            print("unpickled region attr sizes:")
            print("-"*80)
            cu.print_dict_item_sizes(data)
            print("-"*80)
            print("unpickled unfolder attr sizes:")
            print("-"*80)
            cu.print_dict_item_sizes(data['unfolder'].__dict__)
            print("-"*80)

            # ------------------------------------------------------------------
            # PLOT LOTS OF THINGS
            # ------------------------------------------------------------------
            if not args.noBinnedPlots:
                setup = Setup(jet_algo=jet_algo,
                              region=region,
                              angle=angle,
                              output_dir=this_output_dir,
                              has_data=not MC_INPUT)
                hbc = do_binned_plots_per_region_angle(setup,
                                                       do_binned_gen_pt=True,
                                                       do_binned_gen_lambda=True,
                                                       do_binned_reco_pt=True)

                do_all_big_normalised_1d_plots_per_region_angle(setup, hbc)

    print("Saved minimal hists to", output_tfile_slim.GetName())
    output_tfile_slim.Close()
