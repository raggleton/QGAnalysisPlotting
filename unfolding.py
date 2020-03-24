#!/usr/bin/env python


"""TUnfold it all

Thanks to Ashley, Dennis
"""


from __future__ import print_function, division

import os
os.nice(10)
from array import array
import numpy as np
import math
from itertools import product
from copy import copy, deepcopy

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from my_unfolder import MyUnfolder, unpack_unfolding_root_file
from my_unfolder_plotter import MyUnfolderPlotter
from unfolding_regularisation_classes import TauScanner, LCurveScanner
from unfolding_config import get_dijet_config, get_zpj_config
from do_unfolding_plots import Setup, do_all_plots_per_region_angle, do_all_big_1d_plots_per_region_angle
from unfolding_logistics import get_unfolding_argparser, get_unfolding_output_dir, sanitise_args

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()


# Control plot output format
OUTPUT_FMT = "pdf"


def calculate_chi2(hist_test, hist_truth, hist_covariance=None):
    """Calculate the chi2 between 2 histograms, given a covariance matrix

    = (hist1 - hist2)T * cov^-1 * (hist1 - hist2)

    Parameters
    ----------
    hist1 : list, np.array
        List of bin values for one histogram
    hist2 : list, np.array
        List of bin values for other histogram
    cov : np.array
        Covariance matrix

    Returns
    -------
    float
    """
    diff = hist1 - hist2
    diff = diff.reshape(len(diff), 1)
    # for now, hack the inversion, since we know it's diagonal
    inv_cov = np.zeros_like(cov)
    # if True:
    #     for i, e in enumerate(range(cov.shape[0])):
    #         cov_entry = cov[i][i]
    #         if cov_entry != 0:
    #             inv_cov[i][i] = 1./cov_entry
    #         else:
    #             inv_cov[i][i] = 0
    # else:
    inv_cov = np.linalg.inv(cov)
    # print('inv_cov', inv_cov)
    part = np.dot(inv_cov, diff)
    # print('part', part)
    result = np.dot(diff.T, part)
    # print('result', result)
    return result[0][0]


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


def rm_large_rel_error_bins(hist, relative_err_threshold=-1):
    """Reset bins in 2D hist to 0 if error/contents exceeds a certain value"""
    if relative_err_threshold < 0:
        return hist
    new_hist = hist.Clone()
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


def calc_auto_xlim(entries):
    """Figure out x axis range that includes all non-0 bins from all Contributions"""
    if len(entries) == 0:
        return None
    if not isinstance(entries[0], Contribution):
        raise TypeError("calc_auto_xlim: `entries` should be a list of Contributions")

    x_min = 9999999999
    x_max = -9999999999
    for ent in entries:
        obj = ent.obj
        if not isinstance(obj, ROOT.TH1):
            raise TypeError("Cannot handle obj in calc_auto_xlim")
        xax = obj.GetXaxis()
        nbins = obj.GetNbinsX()
        found_min = False
        for i in range(1, nbins+1):
            val = obj.GetBinContent(i)
            if not found_min:
                # find the first non-empty bin
                if val != 0:
                    x_min = min(xax.GetBinLowEdge(i), x_min)
                    found_min = True
            else:
                # find the first empty bin
                if val == 0:
                    x_max = max(xax.GetBinLowEdge(i+1), x_max)
                    break
    if x_max > x_min:
        return (x_min, x_max)
    else:
        # default x max
        return None
        # return (x_min, xax.GetBinLowEdge(nbins+1))


def plot_uncertainty_shifts(total_hist, stat_hist, syst_unfold_hist, systs_shifted, systs, output_filename, title, angle_str):
    """Plot fractional uncertainty shift for a given pt bin

    Systs shifted should be the shifted distributions, not the shifts themselves.
    They are normalised to unity, as are total_hist and stat_hist.
    After this, the fractional diff is calculated, i.e. for *normalised*
    cross sections.

    Parameters
    ----------
    total_hist : TH1
        Distribution with total uncertainty.
    stat_hist : TH1
        Distribution with stat-only uncertainties.
    syst_unfold_hist : TH1
        Distribution with uncertainties = systematics from unfolding (not other systeamtics sources)
    systs_shifted : list[TH1]
        Unfolded distributions with 1-sigma shift from systematics.
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
    total_hist_div_bin_width = qgp.normalise_hist_divide_bin_width(total_hist)
    stat_hist_div_bin_width = qgp.normalise_hist_divide_bin_width(stat_hist)
    syst_unfold_hist_div_bin_width = qgp.normalise_hist_divide_bin_width(syst_unfold_hist)
    for h, syst_dict in zip(systs_shifted, systs):
        # for each syst, normalise and divide by bin width,
        # then subtract total hist, and divide by it
        h_fraction = h.Clone()
        h_fraction = qgp.normalise_hist_divide_bin_width(h_fraction)
        h_fraction.Add(total_hist_div_bin_width, -1)
        h_fraction.Divide(total_hist_div_bin_width)
        hists.append(h_fraction)
        # set to abs values
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

    # Create total and stat error hists
    h_stat = stat_hist_div_bin_width.Clone()
    h_syst = syst_unfold_hist_div_bin_width.Clone()
    h_total = total_hist_div_bin_width.Clone()
    for i in range(1, h_stat.GetNbinsX()+1):
        if total_hist_div_bin_width.GetBinContent(i) > 0:
            h_stat.SetBinContent(i, stat_hist_div_bin_width.GetBinError(i) / total_hist_div_bin_width.GetBinContent(i))
            h_syst.SetBinContent(i, syst_unfold_hist_div_bin_width.GetBinError(i) / total_hist_div_bin_width.GetBinContent(i))
            h_total.SetBinContent(i, total_hist_div_bin_width.GetBinError(i) / total_hist_div_bin_width.GetBinContent(i))
        else:
            h_stat.SetBinContent(i, 0)
            h_syst.SetBinContent(i, 0)
            h_total.SetBinContent(i, 0)
        h_stat.SetBinError(i, 0)
        h_syst.SetBinError(i, 0)
        h_total.SetBinError(i, 0)
    c_stat = Contribution(h_stat,
                         label="Input stats",
                         line_color=ROOT.kRed,
                         line_style=3,
                         line_width=3,
                         marker_size=0,
                         marker_color=ROOT.kRed,
                         )
    entries.append(c_stat)
    c_syst = Contribution(h_syst,
                         label="Response matrix stats",
                         line_color=ROOT.kGray+2,
                         line_style=3,
                         line_width=3,
                         marker_size=0,
                         marker_color=ROOT.kGray+2,
                         )
    entries.append(c_syst)
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
                xlim=calc_auto_xlim(entries),
                xtitle="Particle-level "+angle_str,
                ytitle="| Fractional shift |")
    plot.legend.SetX1(0.55)
    plot.legend.SetY1(0.68)
    plot.legend.SetX2(0.98)
    plot.legend.SetY2(0.88)
    plot.legend.SetNColumns(2)
    plot.left_margin = 0.18
    plot.y_padding_max_linear = 1.4
    plot.plot("NOSTACK HIST")
    plot.save(output_filename)

    plot.y_padding_max_log = 50
    plot.set_logy()
    plot.get_modifier().SetMinimum(1E-4)
    log_filename, ext = os.path.splitext(output_filename)
    plot.save(log_filename+"_log"+ext)

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
    # print(angles)

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
        mc_hname_append = "split" if MC_SPLIT else "all"
        if isinstance(orig_region['mc_tfile'], str):
            orig_region['mc_tfile'] = cu.open_root_file(orig_region['mc_tfile'])
        # hist_mc_reco = cu.get_from_tfile(orig_region['mc_tfile'], "%s/hist_pt_reco_%s" % (orig_region['dirname'], mc_hname_append))
        # hist_mc_gen = cu.get_from_tfile(orig_region['mc_tfile'], "%s/hist_pt_truth_%s" % (orig_region['dirname'], mc_hname_append))
        hist_mc_gen_pt = cu.get_from_tfile(orig_region['mc_tfile'], "%s/hist_pt_truth_%s" % (orig_region['dirname'], mc_hname_append))
        # hist_mc_gen_reco_map = cu.get_from_tfile(orig_region['mc_tfile'], "%s/tu_pt_GenReco_%s" % (orig_region['dirname'], mc_hname_append))

        # Remake gen hist with physical bins & save to file
        # all_pt_bins_gen = np.concatenate((pt_bin_edges_underflow_gen[:-1], pt_bin_edges_gen))
        # hist_mc_gen_pt_physical = ROOT.TH1F("mc_gen_pt", ";p_{T}^{jet} [GeV];N", len(all_pt_bins_gen)-1, array('d', all_pt_bins_gen))
        # update_hist_bin_content(h_orig=hist_mc_gen_pt, h_to_be_updated=hist_mc_gen_pt_physical)
        # update_hist_bin_error(h_orig=hist_mc_gen_pt, h_to_be_updated=hist_mc_gen_pt_physical)

        # Do unfolding for each angle
        # ----------------------------------------------------------------------
        for angle in angles:
            region = copy(orig_region)  # make copy since we might modify it later, e.g. PDF, and want same start for each angle

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

            # Save hists etc to ROOT file for access later
            output_tfile = ROOT.TFile("%s/unfolding_result.root" % (this_output_dir), "RECREATE")

            new_tdir = "%s/%s" % (region['name'], angle.var)
            output_tfile.mkdir(new_tdir)
            this_tdir = output_tfile.Get(new_tdir)
            this_tdir.cd()
            # this_tdir.WriteTObject(hist_mc_gen_pt_physical, "mc_gen_pt")


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

                chosen_bin = 15
                print("1D reco input without background subtraction:", reco_1d.GetBinContent(chosen_bin))
                print("1D reco input with background subtraction:", reco_1d_bg_subtracted.GetBinContent(chosen_bin))
                print("1D reco input fakes:", hist_fakes_reco.GetBinContent(chosen_bin))

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
            # unfolder.SetEpsMatrix(1E-10)

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

            # Add systematic errors as different response matrices
            # ------------------------------------------------------------------
            if args.doExperimentalSysts:
                chosen_rsp_bin = (35, 35)
                print("nominal response bin content for", chosen_rsp_bin, unfolder.response_map.GetBinContent(*chosen_rsp_bin))

                for syst_ind, syst_dict in enumerate(region['experimental_systematics']):
                    print("Adding systematic:", syst_dict['label'])
                    if 'factor' in syst_dict:
                        # special case for e.g. lumi - we construct a reponse hist, and add it using relative mode
                        rel_map = unfolder.response_map.Clone(syst_dict['label']+"Map")
                        for xbin, ybin in product(range(1, rel_map.GetNbinsX()+1), range(1, rel_map.GetNbinsY()+1)):
                            rel_map.SetBinContent(xbin, ybin, syst_dict['factor'])
                            rel_map.SetBinError(xbin, ybin, 0)
                        unfolder.add_sys_error(rel_map, syst_dict['label'], ROOT.TUnfoldDensity.kSysErrModeRelative)

                    # elif 'herwig' in syst_dict['label'].lower():
                    #     if not isinstance(syst_dict['tfile'], ROOT.TFile):
                    #         syst_dict['tfile'] = cu.open_root_file(syst_dict['tfile'])
                    #     map_syst = cu.get_from_tfile(syst_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                    #     map_syst.Scale(unfolder.response_map.Integral() / map_syst.Integral())
                    #     print("    syst bin", chosen_rsp_bin, map_syst.GetBinContent(*chosen_rsp_bin))
                    #     unfolder.add_sys_error(rm_large_rel_error_bins(map_syst), syst_dict['label'], ROOT.TUnfoldDensity.kSysErrModeMatrix)

                    else:
                        if not isinstance(syst_dict['tfile'], ROOT.TFile):
                            syst_dict['tfile'] = cu.open_root_file(syst_dict['tfile'])
                        map_syst = cu.get_from_tfile(syst_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                        print("    syst bin", chosen_rsp_bin, map_syst.GetBinContent(*chosen_rsp_bin))
                        unfolder.add_sys_error(rm_large_rel_error_bins(map_syst), syst_dict['label'], ROOT.TUnfoldDensity.kSysErrModeMatrix)

                    # Plot the reponse matrix for this systematic
                    syst_label_no_spaces = cu.no_space_str(syst_dict['label'])
                    output_filename = "%s/response_map_syst_%s_%s.%s" % (this_output_dir, syst_label_no_spaces, append, unfolder_plotter.output_fmt)
                    title = "%s, %s region, %s, %s" % (jet_algo, region['label'], angle_str, syst_dict['label'])
                    unfolder_plotter.draw_2d_hist(unfolder.syst_maps[syst_dict['label']],
                                                  title=title,
                                                  output_filename=output_filename,
                                                  logz=True,
                                                  draw_bin_lines_x=True,
                                                  draw_bin_lines_y=True,
                                                  canvas_size=(800, 700),
                                                  xtitle='Generator bin',
                                                  ytitle='Detector bin')

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

                if is_herwig or is_pythia8:
                    # SetEpsMatrix ensures rank properly calculated when inverting
                    # Needed if you get message "rank of matrix E 55 expect 170"
                    # And unfolded looks wacko
                    unreg_unfolder.SetEpsMatrix(1E-18)

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
                l_scanner.save_to_tfile(this_tdir)

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
                tau_scanner.save_to_tfile(this_tdir)

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
            unfolded_1d = unfolder.get_output(hist_name="unfolded_1d")
            chosen_bin = 18
            print("Bin %d:" % chosen_bin, unfolded_1d.GetBinContent(chosen_bin))
            print("original uncert:", unfolded_1d.GetBinError(chosen_bin))

            # Calculate experimental uncertainty shifts using results from another unfolding
            # ------------------------------------------------------------------
            ref_tfile_exp = None
            if args.doExperimentalSystsFromFile is not None:
                print("Getting experimental systematics from another file...")
                angle_output_dir = "%s/%s/%s" % (args.doExperimentalSystsFromFile, region['name'], angle.var)
                this_root_filename = os.path.join(angle_output_dir, "unfolding_result.root")
                if not os.path.isfile(this_root_filename):
                    raise IOError("Cannot find systematics file, %s" % this_root_filename)
                ref_tfile_exp = cu.TFileCacher(this_root_filename)
                region_copy = deepcopy(region)
                unpack_dict = unpack_unfolding_root_file(ref_tfile_exp, region_copy, angle)
                reference_unfolder = unpack_dict['unfolder']
                ref_unfolded = reference_unfolder.unfolded
                for syst_label in reference_unfolder.syst_maps.keys():
                    # Note that this is potentially a bit dodgy - the
                    # underlying TUnfoldensity object has no knowledge of this syst
                    # Get shift on absolute result
                    new_delta_shift = reference_unfolder.get_delta_sys_shift(syst_label).Clone()
                    # Calc relative difference
                    new_delta_shift.Divide(ref_unfolded)
                    # Apply to our nominal unfolded result & store
                    # We only need the delta shift (like normal),
                    # since _post_process() calcualte everything from that
                    new_delta_shift.Multiply(unfolder.unfolded)
                    unfolder.syst_shifts[syst_label] = new_delta_shift
                    unfolder.syst_maps[syst_label] = reference_unfolder.syst_maps[syst_label]
                    unfolder.systs_shifted[syst_label] = None  # gets calculated in get_syst_shifted_hist()
                    unfolder.syst_ematrices[syst_label] = None

                # update region info
                # TODO what if the config has fewer than in the reference unfolder?
                region['experimental_systematics'] = [syst_dict for syst_dict in region_copy['experimental_systematics']
                                                      if syst_dict['label'] in reference_unfolder.syst_maps.keys()]

            # Do lots of extra gubbins, like caching matrices,
            # creating unfolded hists with different levels of uncertianties,
            # ------------------------------------------------------------------
            unfolder._post_process()

            # Get various error matrices
            # ------------------------------------------------------------------
            # stat errors only - do before or after systematics?
            print("stat uncert:", unfolder.get_unfolded_with_ematrix_stat().GetBinError(chosen_bin))
            print("new uncert:", unfolder.get_output().GetBinError(chosen_bin))

            # hist1, err1 = cu.th1_to_arr(unfolded_1d)
            # hist2, err2 = cu.th1_to_arr(hist_mc_gen)
            # cov, cov_err = cu.th2_to_arr(ematrix_total)
            # chi2 = calculate_chi2(hist1, hist2, cov)
            # print("my chi2 =", chi2)

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

                print("projection reco #bins:", proj_reco.GetNbinsX())
                print("response map # bins x:", unfolder.response_map.GetNbinsX())
                print("response map # bins y:", unfolder.response_map.GetNbinsY())
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
            smeared_chi2, smeared_ndf, smeared_p = unfolder.calculate_chi2(one_hist=unfolder.get_folded_mc_truth(),
                                                                           other_hist=unfolder.hist_mc_reco_bg_subtracted,
                                                                           cov_inv_matrix=unfolder.get_vyy_inv_ndarray(),
                                                                           # cov_inv_matrix=unfolder.get_vyy_inv_no_bg_ndarray(),
                                                                           detector_space=True,
                                                                           ignore_underflow_bins=True,
                                                                           debugging_dir=os.path.join(this_output_dir, 'smeared_chi2_debug'))
            print('smeared chi2:', smeared_chi2, smeared_ndf, smeared_chi2/smeared_ndf, smeared_p)

            unfolded_chi2, unfolded_ndf, unfolded_p = unfolder.calculate_chi2(one_hist=unfolder.unfolded,
                                                                              other_hist=unfolder.hist_truth,
                                                                              cov_inv_matrix=unfolder.get_vxx_inv_ndarray(),
                                                                              detector_space=False,
                                                                              ignore_underflow_bins=True,
                                                                              debugging_dir=os.path.join(this_output_dir, 'unfolded_chi2_debug'))
            print('unfolded chi2:', unfolded_chi2, unfolded_ndf, unfolded_chi2/unfolded_ndf, unfolded_p)

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

            this_tdir.cd()
            alt_tdir = this_tdir.mkdir("alt_response_%s" % cu.no_space_str(region['alt_mc_label']))
            alt_tdir.cd()

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

                # this_tdir.cd()
                # alt_tdir = this_tdir.mkdir("alt_response_%s" % cu.no_space_str(region['alt_mc_label']))
                # alt_tdir.cd()

                is_herwig = "Herwig" in region['alt_mc_label']
                is_pythia8 = region['alt_mc_label'] == "Pythia8"  # not MG+Pythia9
                if is_herwig or is_pythia8:
                    # SetEpsMatrix ensures rank properly calculated when inverting
                    # Needed if you get message "rank of matrix E 55 expect 170"
                    # And unfolded looks wacko
                    alt_unfolder.SetEpsMatrix(1E-18)

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
                    alt_l_scanner.save_to_tfile(alt_tdir)

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
                    alt_tau_scanner.save_to_tfile(alt_tdir)

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
                alt_unfolded_1d = alt_unfolder.get_output(hist_name="alt_unfolded_1d")
                print("Bin %d:" % chosen_bin, alt_unfolded_1d.GetBinContent(chosen_bin))
                print("original uncert:", alt_unfolded_1d.GetBinError(chosen_bin))
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
                if alt_unfolder:
                    smeared_alt_chi2, smeared_alt_ndf, smeared_alt_p = unfolder.calculate_chi2(one_hist=alt_unfolder.get_folded_mc_truth(),
                                                                                               other_hist=unfolder.hist_mc_reco_bg_subtracted,
                                                                                               cov_inv_matrix=unfolder.get_vyy_inv_ndarray(),
                                                                                               # cov_inv_matrix=unfolder.get_vyy_inv_no_bg_ndarray(),
                                                                                               detector_space=True,
                                                                                               ignore_underflow_bins=True,
                                                                                               debugging_dir=os.path.join(this_output_dir, 'smeared_alt_chi2_debug'))
                    print('smeared chi2 (alt MC):', smeared_alt_chi2, smeared_alt_ndf, smeared_alt_chi2/smeared_alt_ndf, smeared_alt_p)

                print(unfolder.unfolded.Integral())
                print(alt_hist_mc_gen.Integral())
                if alt_hist_mc_gen:
                    unfolded_alt_chi2, unfolded_alt_ndf, unfolded_alt_p = unfolder.calculate_chi2(one_hist=unfolder.unfolded,
                                                                                                  other_hist=alt_hist_mc_gen,
                                                                                                  cov_inv_matrix=unfolder.get_vxx_inv_ndarray(),
                                                                                                  detector_space=False,
                                                                                                  ignore_underflow_bins=True,
                                                                                                  debugging_dir=os.path.join(this_output_dir, 'unfolded_alt_chi2_debug'))
                    print('unfolded chi2 (alt MC):', unfolded_alt_chi2, unfolded_alt_ndf, unfolded_alt_chi2/unfolded_alt_ndf, unfolded_alt_p)

                region['alt_unfolder'] = alt_unfolder

                # Save important stuff to TFile
                # --------------------------------------------------------------
                alt_unfolder.save_to_tfile(alt_tdir)

            # Bit gnarly - have to save this stuff manually
            alt_tdir.WriteTObject(alt_hist_mc_gen, "alt_hist_mc_gen")
            alt_tdir.WriteTObject(alt_hist_mc_reco, "alt_hist_mc_reco")
            alt_tdir.WriteTObject(alt_hist_mc_reco_bg_subtracted, "alt_hist_mc_reco_bg_subtracted")
            alt_tdir.WriteTObject(alt_hist_mc_reco_bg_subtracted_gen_binning, "alt_hist_mc_reco_bg_subtracted_gen_binning")

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

                    this_tdir.cd()
                    syst_tdir_name = "modelSyst_"+syst_label_no_spaces
                    syst_tdir = this_tdir.mkdir(syst_tdir_name)
                    syst_tdir.cd()

                    syst_output_dir = this_output_dir+"/modelSyst_"+syst_label_no_spaces
                    syst_unfolder_plotter = MyUnfolderPlotter(syst_unfolder, is_data=False)
                    syst_plot_args = dict(output_dir=syst_output_dir,
                                          append=append)

                    if is_herwig:
                        # SetEpsMatrix ensures rank properly calculated when inverting
                        # Needed if you get message "rank of matrix E 55 expect 170"
                        # And unfolded looks wacko
                        syst_unfolder.SetEpsMatrix(1E-18)

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

                    # Add systematic errors as different response matrices
                    # ----------------------------------------------------
                    # if args.doExperimentalSysts:
                    #     chosen_rsp_bin = (18, 18)
                    #     print("nominal response bin content for", chosen_rsp_bin, syst_unfolder.response_map.GetBinContent(*chosen_rsp_bin))
                    #     for exp_dict in region['experimental_systematics']:
                    #         if not isinstance(exp_dict['tfile'], ROOT.TFile):
                    #             exp_dict['tfile'] = cu.open_root_file(exp_dict['tfile'])
                    #         map_syst = cu.get_from_tfile(exp_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                    #         print("Adding systematic:", exp_dict['label'])
                    #         print("    syst bin", chosen_rsp_bin, map_syst.GetBinContent(*chosen_rsp_bin))
                    #         syst_unfolder.add_sys_error(rm_large_rel_error_bins(map_syst), exp_dict['label'], ROOT.TUnfoldDensity.kSysErrModeMatrix)


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
                        syst_l_scanner.save_to_tfile(syst_tdir)

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
                        syst_tau_scanner.save_to_tfile(syst_tdir)

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
                    syst_unfolded_1d = syst_unfolder.get_output(hist_name="syst_%s_unfolded_1d" % (syst_label_no_spaces))
                    print("Bin %d:" % (chosen_bin), syst_unfolded_1d.GetBinContent(chosen_bin))
                    print("original uncert:", syst_unfolded_1d.GetBinError(chosen_bin))
                    syst_unfolder._post_process()

                    syst_title = "%s\n%s region, %s, %s input" % (jet_algo, region['label'], angle_str, syst_label)
                    syst_unfolder_plotter.draw_unfolded_1d(title=syst_title, **syst_plot_args)

                    region['model_systematics'][ind]['unfolder'] = syst_unfolder

                    # Save important stuff to TFile
                    # --------------------------------------------------------------
                    syst_unfolder.save_to_tfile(syst_tdir)

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

                # Update main unfolder with these unfolded variations
                # --------------------------------------------------------------
                unfolder.create_normalised_scale_syst_uncertainty(region['model_systematics'])

                # Do big 1D plot of nominal & all systs
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

                entries.extend(syst_entries)

                title = "%s\n%s region, %s" % (jet_algo, region['label'], angle_str)
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
                # # plot.container.SetMinimum(0.001)
                plot.legend.SetY1NDC(0.67)
                plot.legend.SetY2NDC(0.88)
                plot.legend.SetX1NDC(0.65)
                plot.legend.SetX2NDC(0.88)
                output_filename = "%s/unfolded_1d_modelSyst_%s.%s" % (this_output_dir, append, syst_unfolder_plotter.output_fmt)
                plot.save(output_filename)


            # Load model (scale) systs from another reference file, and calc fractional
            # uncertainty, apply to this unfolded result
            ref_tfile_model = None
            if args.doModelSystsFromFile is not None:
                model_dir = "%s/%s/%s" % (args.doModelSystsFromFile, region['name'], angle.var)
                this_root_filename = os.path.join(model_dir, "unfolding_result.root")
                # do NOT use save variable name as for exp or pdf systs...things will disappear
                ref_tfile_model = cu.TFileCacher(this_root_filename)
                region_copy = orig_region.copy()
                reference_dict = unpack_unfolding_root_file(ref_tfile_model,
                                                            region_copy,
                                                            angle,
                                                            do_alt_response=False,
                                                            do_model_systs=True,
                                                            do_pdf_systs=False)
                # update original region object
                region['model_systematics'] = region_copy['model_systematics']

                # Add dummy object to hist_bin_chopper for later, so we can directly manipulate the cache
                uncert_name = "scale_uncert"
                unfolder.hist_bin_chopper.add_obj(uncert_name, unfolder.get_unfolded_with_ematrix_stat())
                unfolder.hist_bin_chopper.add_obj("unfolded_stat_err", unfolder.get_unfolded_with_ematrix_stat())

                for ibin_pt in range(len(unfolder.pt_bin_edges_gen[:-1])):
                    # Calculate scale uncertainty by taking relative uncertainty
                    # from reference file, and applying to this result
                    key = unfolder.hist_bin_chopper._generate_key(uncert_name,
                                                                  ind=ibin_pt,
                                                                  axis='pt',
                                                                  do_norm=True,
                                                                  do_div_bin_width=True,
                                                                  binning_scheme='generator')
                    ref_scale_syst = reference_dict['unfolder'].hist_bin_chopper._cache[key]
                    scale_syst = unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width("unfolded_stat_err", ibin_pt, binning_scheme='generator').Clone()
                    for i in range(1, scale_syst.GetNbinsX()+1):
                        if ref_scale_syst.GetBinContent(i) != 0:
                            rel_err = ref_scale_syst.GetBinError(i) / ref_scale_syst.GetBinContent(i)
                        else:
                            rel_err = 0
                        scale_syst.SetBinError(i, rel_err * scale_syst.GetBinContent(i))
                    unfolder.hist_bin_chopper._cache[key] = scale_syst

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
                # first construct all new systematic variations
                original_pdf_dict = region['pdf_systematics'][0]
                original_pdf_dict['label'] = '_PDF_template'  # initial _ to ignore it later on
                tfile = original_pdf_dict['tfile']
                if not isinstance(tfile, ROOT.TFile):
                    tfile = cu.open_root_file(tfile)

                region['pdf_systematics']  = []
                for pdf_ind in original_pdf_dict['variations']:
                    mc_hname_append = "split" if MC_SPLIT else "all"
                    mc_hname_append = "all"  # want all stats since we are already independent of the response matrix
                    region['pdf_systematics'].append(
                        {
                            "label": "PDF_%d" % (pdf_ind),
                            "hist_reco": cu.get_from_tfile(tfile, "%s/hist_%s_reco_%s_PDF_%d" % (region['dirname'], angle_shortname, mc_hname_append, pdf_ind)),
                            "hist_gen": cu.get_from_tfile(tfile, "%s/hist_%s_gen_%s_PDF_%d" % (region['dirname'], angle_shortname, mc_hname_append, pdf_ind)),
                            "colour": ROOT.kCyan+2,
                        })

                    if mc_hname_append == 'all' and MC_SPLIT:
                        # Since the nominal MC only has 20% of the stats,
                        # need to scale this as well otherwise it will look weird
                        region['pdf_systematics'][-1]['hist_reco'].Scale(0.2)
                        region['pdf_systematics'][-1]['hist_gen'].Scale(0.2)

                # Now run over all variations like for model systs
                for ind, pdf_dict in enumerate(region['pdf_systematics']):
                    pdf_label = pdf_dict['label']
                    pdf_label_no_spaces = cu.no_space_str(pdf_label)

                    if pdf_label.startswith("_"):
                        continue

                    print("*" * 80)
                    print("*** Unfolding with alternate input:", pdf_label, "(%d/%d) ***" % (ind+1, len(region['pdf_systematics'])))
                    print("*" * 80)

                    hist_pdf_reco = pdf_dict['hist_reco']
                    hist_pdf_gen = pdf_dict['hist_gen']

                    pdf_unfolder = MyUnfolder(response_map=unfolder.response_map,
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

                    this_tdir.cd()
                    pdf_tdir_name = "pdfSyst_"+pdf_label_no_spaces
                    pdf_tdir = this_tdir.mkdir(pdf_tdir_name)
                    pdf_tdir.cd()

                    pdf_unfolder_plotter = MyUnfolderPlotter(pdf_unfolder, is_data=False)
                    pdf_output_dir = this_output_dir+"/pdfSyst/"+pdf_label_no_spaces
                    pdf_plot_args = dict(output_dir=pdf_output_dir,
                                         append=append)

                    if SUBTRACT_FAKES:
                        # Use the background template from the nominal MC
                        # (since we're only testing different input shapes,
                        # and our bkg estimate is always from MC)
                        hist_fakes_pdf = hist_fakes_reco_fraction.Clone("hist_fakes_pdf_%s" % pdf_label_no_spaces)
                        hist_fakes_pdf.Multiply(hist_pdf_reco)

                    hist_pdf_reco_bg_subtracted = hist_pdf_reco.Clone()
                    hist_pdf_reco_bg_subtracted.Add(hist_fakes_pdf, -1)

                    # Set what is to be unfolded
                    # --------------------------------------------------------------
                    pdf_unfolder.set_input(input_hist=hist_pdf_reco,
                                           hist_truth=hist_pdf_gen,
                                           hist_mc_reco=hist_pdf_reco,
                                           hist_mc_reco_bg_subtracted=hist_pdf_reco_bg_subtracted,
                                           bias_factor=args.biasFactor)

                    # Subtract fakes (treat as background)
                    # --------------------------------------------------------------
                    if SUBTRACT_FAKES:
                        pdf_unfolder.subtract_background(hist_fakes_pdf, "fakes")

                    # also show nominal bg-subtracted input for comparison
                    ocs = [
                        Contribution(unfolder.input_hist_bg_subtracted,
                                     label='Nominal unfolding input (bg-subtracted)',
                                     line_color=ROOT.kRed, line_width=1)
                    ]
                    pdf_unfolder_plotter.draw_detector_1d(do_reco_data=False,
                                                          do_reco_data_bg_sub=False,
                                                          do_reco_bg=True,
                                                          do_reco_mc=False,
                                                          do_reco_mc_bg_sub=True,
                                                          other_contributions=ocs,
                                                          output_dir=pdf_plot_args['output_dir'],
                                                          append='bg_fakes_subtracted_%s' % append,
                                                          title="%s region, %s, %s" % (region['label'], angle_str, pdf_label))

                    # Add systematic errors as different response matrices
                    # ----------------------------------------------------
                    # if args.doExperimentalSysts:
                    #     chosen_rsp_bin = (18, 18)
                    #     print("nominal response bin content for", chosen_rsp_bin, pdf_unfolder.response_map.GetBinContent(*chosen_rsp_bin))
                    #     for exp_dict in region['experimental_systematics']:
                    #         print("Adding systematic:", exp_dict['label'])
                    #         if 'factor' in exp_dict:
                    #             # special case for e.g. lumi - we construct a reponse hist, and add it using relative mode
                    #             rel_map = pdf_unfolder.response_map.Clone(exp_dict['label']+"MapPDF")
                    #             for xbin, ybin in product(range(1, rel_map.GetNbinsX()+1), range(1, rel_map.GetNbinsY()+1)):
                    #                 rel_map.SetBinContent(xbin, ybin, exp_dict['factor'])
                    #                 rel_map.SetBinError(xbin, ybin, 0)
                    #             pdf_unfolder.add_sys_error(rel_map, exp_dict['label'], ROOT.TUnfoldDensity.kSysErrModeRelative)
                    #         else:
                    #             if not isinstance(exp_dict['tfile'], ROOT.TFile):
                    #                 exp_dict['tfile'] = cu.open_root_file(exp_dict['tfile'])
                    #             map_syst = cu.get_from_tfile(exp_dict['tfile'], "%s/tu_%s_GenReco_all" % (region['dirname'], angle_shortname))
                    #             print("    syst bin", chosen_rsp_bin, map_syst.GetBinContent(*chosen_rsp_bin))
                    #             pdf_unfolder.add_sys_error(map_syst, exp_dict['label'], ROOT.TUnfoldDensity.kSysErrModeMatrix)

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
                        pdf_l_scanner.save_to_tfile(pdf_tdir)

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
                        pdf_tau_scanner.save_to_tfile(pdf_tdir)

                    region['pdf_systematics'][ind]['tau'] = pdf_tau

                    # Do unfolding!
                    # --------------------------------------------------------------
                    pdf_unfolder.do_unfolding(pdf_tau)
                    pdf_unfolded_1d = pdf_unfolder.get_output(hist_name="syst_%s_unfolded_1d" % (pdf_label_no_spaces))
                    print("Bin %d:" % (chosen_bin), pdf_unfolded_1d.GetBinContent(chosen_bin))
                    print("original uncert:", pdf_unfolded_1d.GetBinError(chosen_bin))
                    pdf_unfolder._post_process()

                    pdf_title = "%s\n%s region, %s, %s input" % (jet_algo, region['label'], angle_str, pdf_label)
                    pdf_unfolder_plotter.draw_unfolded_1d(title=pdf_title, **pdf_plot_args)

                    region['pdf_systematics'][ind]['unfolder'] = pdf_unfolder

                    # Save important stuff to TFile
                    # ----------------------------------------------------------
                    pdf_unfolder.save_to_tfile(pdf_tdir)

                unfolder.create_normalised_pdf_syst_uncertainty(region['pdf_systematics'])

            # Load PDF syst from another reference file, and calc fractional
            # uncertainty, apply to this unfolded result
            ref_tfile_pdf = None
            if args.doPDFSystsFromFile is not None:
                pdf_dir = "%s/%s/%s" % (args.doPDFSystsFromFile, region['name'], angle.var)
                this_root_filename = os.path.join(pdf_dir, "unfolding_result.root")
                # do NOT use save variable name as for exp or model systs...things will disappear
                ref_tfile_pdf = cu.TFileCacher(this_root_filename)
                region_copy_pdf = orig_region.copy()
                reference_dict = unpack_unfolding_root_file(ref_tfile_pdf,
                                                            region_copy_pdf,
                                                            angle,
                                                            do_alt_response=False,
                                                            do_pdf_systs=True,
                                                            do_model_systs=False)
                # update original region object
                region['pdf_systematics'] = region_copy_pdf['pdf_systematics']

                # Add dummy object to hist_bin_chopper for later, so we can directly manipulate the cache
                uncert_name = "pdf_uncert"
                unfolder.hist_bin_chopper.add_obj(uncert_name, unfolder.get_unfolded_with_ematrix_stat())
                unfolder.hist_bin_chopper.add_obj("unfolded_stat_err", unfolder.get_unfolded_with_ematrix_stat())

                for ibin_pt in range(len(unfolder.pt_bin_edges_gen[:-1])):
                    # Calculate PDF uncertainty by taking relative uncertainty
                    # from reference file, and applying to this result
                    key = unfolder.hist_bin_chopper._generate_key(uncert_name,
                                                                  ind=ibin_pt,
                                                                  axis='pt',
                                                                  do_norm=True,
                                                                  do_div_bin_width=True,
                                                                  binning_scheme='generator')
                    ref_pdf_syst = reference_dict['unfolder'].hist_bin_chopper._cache[key]
                    pdf_syst = unfolder.hist_bin_chopper.get_pt_bin_normed_div_bin_width("unfolded_stat_err", ibin_pt, binning_scheme='generator').Clone()
                    for i in range(1, pdf_syst.GetNbinsX()+1):
                        if ref_pdf_syst.GetBinContent(i) != 0:
                            rel_err = ref_pdf_syst.GetBinError(i) / ref_pdf_syst.GetBinContent(i)
                        else:
                            rel_err = 0
                        pdf_syst.SetBinError(i, rel_err * pdf_syst.GetBinContent(i))
                    unfolder.hist_bin_chopper._cache[key] = pdf_syst


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

            # Save everything to TFile
            print("DONE...saving unfolder to ROOT file")
            unfolder.save_to_tfile(this_tdir)

            # ------------------------------------------------------------------
            # PLOT LOTS OF THINGS
            # ------------------------------------------------------------------
            if not args.noBinnedPlots:
                setup = Setup(jet_algo=jet_algo,
                              region=region,
                              angle=angle,
                              output_dir=this_output_dir,
                              has_data=not MC_INPUT)
                unfold_dict = dict(
                    unfolder=unfolder,
                    unreg_unfolder=unreg_unfolder,
                    alt_unfolder=alt_unfolder,
                    alt_hist_truth=alt_hist_mc_gen,
                    alt_hist_reco=alt_hist_mc_reco,
                    alt_hist_reco_bg_subtracted=alt_hist_mc_reco_bg_subtracted,
                    alt_hist_reco_bg_subtracted_gen_binning=alt_hist_mc_reco_bg_subtracted_gen_binning,
                )
                hbc = do_all_plots_per_region_angle(setup, unfold_dict)

                do_all_big_1d_plots_per_region_angle(setup,
                                                     unfold_dict,
                                                     hbc)

            # DO SUMMARY PLOT
            # ------------------------------------------------------------------
            # if args.doSummaryPlot:
            #     marker = ""
            #     if "_" in angle.name or "^" in angle.name:
            #         marker = "$"
            #     var_label = "Particle-level " + marker + angle.name + marker + " ($%s$)" % angle.lambda_str
            #     v = "%s_vs_pt" % (angle.var)
            #     bins = [(pt_bin_edges_gen[i], pt_bin_edges_gen[i+1]) for i in range(len(pt_bin_edges_gen)-1)]
            #     print("Making summary plot from pt bins:", bins)
            #     xlim = (50, 614) if "ZPlusJets" in region['name'] else (50, 2000)
            #     region_label = region['label'].replace("Dijet", "dijet")  # to ensure correct capitalisation
            #     qgp.do_mean_rms_summary_plot(summary_1d_entries, bins,
            #                                  "%s/%s_box_dijet_mpl.%s" % (this_output_dir, v, OUTPUT_FMT),
            #                                  var_label=var_label,
            #                                  xlim=xlim,
            #                                  region_title=region_label)

    print("Saved hists to", output_tfile.GetName())
