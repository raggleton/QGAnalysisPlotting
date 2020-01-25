#!/usr/bin/env python


"""
Do all the unfolding plots: per pT bin, per lambda bin, summary plot
"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp
from my_unfolder import MyUnfolder, unfolder_from_tdir
from my_unfolder_plotter import MyUnfolderPlotter
from unfolding_config import get_dijet_config, get_zpj_config


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat(".3f")
ROOT.gStyle.SetHistTopMargin(0.)


def setup_regions(args):
    regions = []
    if args.doDijetCentral:
        dijet_region_central_dict = get_dijet_config(args.source, central=True, groomed=False)
        regions.append(dijet_region_central_dict)

    if args.doDijetForward:
        dijet_region_forward_dict = get_dijet_config(args.source, central=False, groomed=False)
        regions.append(dijet_region_forward_dict)

    if args.doDijetCentralGroomed:
        dijet_region_central_groomed_dict = get_dijet_config(args.source, central=True, groomed=True)
        regions.append(dijet_region_central_groomed_dict)

    if args.doDijetForwardGroomed:
        dijet_region_forward_groomed_dict = get_dijet_config(args.source, central=False, groomed=True)
        regions.append(dijet_region_forward_groomed_dict)

    if args.doZPJ:
        zpj_region_dict = get_zpj_config(args.source, groomed=False)
        regions.append(zpj_region_dict)

    if args.doZPJGroomed:
        zpj_region_groomed_dict = get_zpj_config(args.source, groomed=True)
        regions.append(zpj_region_groomed_dict)
    return regions


def unpack_unfolding_root_file(input_tfile, region, angle):
    input_tdir_name = "%s/%s" % (region['name'], angle.var)
    input_tdir = input_tfile.Get(input_tdir_name)
    cu.check_root_obj(input_tdir)
    unfolder = unfolder_from_tdir(input_tdir)
    print("Loaded main unfolder")

    list_of_obj = cu.get_list_of_element_names(input_tdir)

    # Get unregularised unfolder, if available
    unreg_tdir = [x for x in list_of_obj if x.startswith("unreg_unfolder")]
    unreg_unfolder = None
    if len(unreg_tdir) == 1:
        unreg_unfolder = unfolder_from_tdir(input_tfile.Get(os.path.join(input_tdir_name, unreg_tdir[0])))
        print("Loaded comparison unregularised unfolder")

    # Get alternate response object, if it exists
    alt_tdir = [x for x in list_of_obj if x.startswith("alt_response_")]
    alt_unfolder = None
    alt_unfolder_name = None
    alt_hist_truth = None
    if len(alt_tdir)  == 1:
        alt_unfolder = unfolder_from_tdir(input_tfile.Get(os.path.join(input_tdir_name, alt_tdir[0])))
        region['alt_unfolder'] = alt_unfolder
        alt_unfolder_name = alt_tdir[0].replace("alt_response_", "").replace("_", " ")
        if region['alt_mc_label'] != alt_unfolder_name:
            raise RuntimeError("Bad unpacking of alt response unfolder: expected %s, got %s" % (region['alt_mc_label'], alt_unfolder_name))
        print("Loaded alt unfolder")
        alt_hist_truth = input_tfile.Get(os.path.join(input_tdir_name, alt_tdir[0], "alt_hist_mc_gen"))
    if len(alt_tdir) > 1:
        raise RuntimeError(">1 alt_response?! %s" % (alt_tdir))

    # Get model systs
    print(list_of_obj)
    model_tdirs = [x for x in list_of_obj if x.startswith("modelSyst_")]
    if len(model_tdirs) > 0:
        for model_tdir_name in model_tdirs:
            syst_name = model_tdir_name.replace("modelSyst_", "").replace("_", " ")
            this_one = [x for x in region['model_systematics'] if x['label'] == syst_name]
            if len(this_one) == 0:
                print("No entry for model systematic", syst_name, "- skipping")
                continue
            # TODO: check it agrees with region dict?
            this_one[0]['unfolder'] = unfolder_from_tdir(input_tfile.Get(os.path.join(input_tdir_name, model_tdir_name)))
            print("Loaded", len(model_tdirs), "model systematic unfolders")
    # remove entries without an unfolder
    region['model_systematics'] = [k for k in region['model_systematics']
                                 if k.get('unfolder', None) is not None]

    # Get PDF systs
    # For some reason, this is done as a list instead of dict
    pdf_tdirs = [x for x in list_of_obj if x.startswith("pdfSyst_")]
    if len(pdf_tdirs) > 0:
        # Remove original, construct all other
        region['pdf_systematics'] = []

        for pdf_tdir_name in pdf_tdirs:
            pdf_name = pdf_tdir_name.replace("pdfSyst_", "")
            region['pdf_systematics'].append({
                'label': pdf_name,
                'unfolder': unfolder_from_tdir(input_tfile.Get(os.path.join(input_tdir_name, pdf_tdir_name))),
                'colour': ROOT.kCyan+2,
            })
        print("Loaded", len(pdf_tdirs), "PDF systematic unfolders")
    # remove entries without an unfolder
    region['pdf_systematics'] = [k for k in region['pdf_systematics']
                                 if k.get('unfolder', None) is not None]

    return dict(
        unfolder=unfolder,
        alt_unfolder=alt_unfolder,
        alt_hist_truth=alt_hist_truth
    )

class HistBinChopper(object):

    def __init__(self, unfolder):
        self.unfolder = unfolder
        self.objects = {}
        self._cache = {}

    def add_obj(self, name, obj):
        # TODO: allow overwrite?
        self.objects[name] = obj

    def get_pt_bin(self, name, ind):
        if name not in self.objects:
            raise KeyError("No %s in objects" % name)
        this_name = name+"_pt_bin_%d" % (ind)
        if this_name not in self._cache:
            self._cache[this_name] = self.unfolder.get_var_hist_pt_binned(self.objects[name], ind, binning_scheme="generator")
        return self._cache[this_name]

    def get_pt_bin_div_bin_width(self, name, ind):
        if name not in self.objects:
            raise KeyError("No %s in objects" % name)
        this_name = name+"_pt_bin_%d_divBinWidth" % (ind)
        if this_name not in self._cache:
            self._cache[this_name] = qgp.hist_divide_bin_width(self.get_pt_bin(name, ind))
        return self._cache[this_name]

    def get_pt_bin_normed_div_bin_width(self, name, ind):
        if name not in self.objects:
            raise KeyError("No %s in objects" % name)
        this_name = name+"_pt_bin_%d_norm_divBinWidth" % (ind)
        if this_name not in self._cache:
            self._cache[this_name] = qgp.normalise_hist_divide_bin_width(self.get_pt_bin(name, ind))
        return self._cache[this_name]

    def get_lambda_bin(self, name, ind):
        if name not in self.objects:
            raise KeyError("No %s in objects" % name)
        this_name = name+"_lambda_bin_%d" % (ind)
        if this_name not in self._cache:
            self._cache[this_name] = self.unfolder.get_pt_hist_var_binned(self.objects[name], ind, binning_scheme="generator")
        return self._cache[this_name]

    def get_lambda_bin_div_bin_width(self, name, ind):
        if name not in self.objects:
            raise KeyError("No %s in objects" % name)
        this_name = name+"_lambda_bin_%d_divBinWidth" % (ind)
        if this_name not in self._cache:
            self._cache[this_name] = qgp.hist_divide_bin_width(self.get_lambda_bin(name, ind))
        return self._cache[this_name]

    def get_lambda_bin_normed_div_bin_width(self, name, ind):
        if name not in self.objects:
            raise KeyError("No %s in objects" % name)
        this_name = name+"_lambda_bin_%d_norm_divBinWidth" % (ind)
        if this_name not in self._cache:
            self._cache[this_name] = qgp.normalise_hist_divide_bin_width(self.get_lambda_bin(name, ind))
        return self._cache[this_name]


class Setup(object):

    def __init__(self, jet_algo, region, angle, output_dir='.', has_data=False, is_ave_pt_binning=False, bin_edge_low=-1, bin_edge_high=-1):
        self.jet_algo = jet_algo
        self.region = region
        self.pt_str = "#LT p_{T}^{jet} #GT" if is_ave_pt_binning else "p_{T}^{jet}"
        self.has_data = has_data
        self.angle = angle
        angle_prepend = "groomed " if "groomed" in region['name'] else ""
        this_angle_name = angle.name
        if (angle_prepend != ""
            and this_angle_name != 'LHA'
            and "_{T}" not in this_angle_name
            and "PUPPI" not in this_angle_name):
            # lower case if Groomed..., but be careful of e.g. pTD, LHA
            this_angle_name = this_angle_name[0].lower() + this_angle_name[1:]
        self.angle_str = "{prepend}{name} ({lambda_str})".format(prepend=angle_prepend,
                                                                 name=this_angle_name,
                                                                 lambda_str=angle.lambda_str)
        self.particle_title = "Particle-level " + self.angle_str
        self.detector_title = "Detector-level " + self.angle_str
        self.pt_bin_normalised_differential_label = "#frac{1}{d#sigma/dp_{T}} #frac{d^{2}#sigma}{dp_{T} d%s}" % (angle.lambda_str)
        self.lambda_bin_normalised_differential_label = "#frac{1}{d#sigma/d%s}} #frac{d^{2}#sigma}{dp_{T} d%s}" % (angle.lambda_str, angle.lambda_str)
        self.bin_edge_low = bin_edge_low
        self.bin_edge_high = bin_edge_high
        self.output_dir = output_dir
        self.output_fmt = 'pdf'
        self.append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name


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

#Should all this stuff be put in a class?
PT_BIN_PLOT_ARGS = dict(
    what="hist",
    subplot_type='ratio',
    subplot_title="Unfolded / Gen",
    subplot_limits=(0.75, 1.25),
)


PLOT_COLOURS = dict(
    gen_colour=ROOT.kRed,
    unfolded_basic_colour=ROOT.kAzure+7,
    unfolded_stat_colour=ROOT.kAzure+7,
    unfolded_total_colour=ROOT.kBlack,
    unfolded_unreg_colour=ROOT.kViolet+2,
    alt_gen_colour=ROOT.kViolet+1,
    alt_unfolded_colour=ROOT.kBlue+1,
)


LW = 2

def get_pt_bin_title(setup, region, bin_edge_low, bin_edge_high):
    title = (("{jet_algo}\n"
              "{region_label} region\n"
              "{bin_edge_low:g} < {pt_str} < {bin_edge_high:g} GeV")
             .format(
                jet_algo=setup.jet_algo,
                region_label=region['label'],
                pt_str=setup.pt_str,
                bin_edge_low=bin_edge_low,
                bin_edge_high=bin_edge_high
            ))
    return title


def plot_unfolded_unnormalised_pt_bin(ibin_pt, unfolder, hist_bin_chopper, setup):
    mc_gen_hist_bin = hist_bin_chopper.get_pt_bin_div_bin_width('hist_truth', ibin_pt)
    unfolded_hist_bin_stat_errors = hist_bin_chopper.get_pt_bin_div_bin_width('unfolded_stat_err', ibin_pt)
    unfolded_hist_bin_total_errors = hist_bin_chopper.get_pt_bin_div_bin_width('unfolded', ibin_pt)

    # unnormalised version
    entries = [
        Contribution(mc_gen_hist_bin,
                     label="Generator (%s)" % (region['mc_label']),
                     line_color=PLOT_COLOURS['gen_colour'], line_width=LW,
                     marker_color=PLOT_COLOURS['gen_colour'], marker_size=0,
                     normalise_hist=False),
        Contribution(unfolded_hist_bin_total_errors,
                     label="Unfolded (#tau = %.3g) (total err)" % (unfolder.tau),
                     line_color=PLOT_COLOURS['unfolded_total_colour'], line_width=LW, line_style=1,
                     marker_color=PLOT_COLOURS['unfolded_total_colour'],# marker_style=20, marker_size=0.75,
                     subplot=mc_gen_hist_bin,
                     normalise_hist=False),
        Contribution(unfolded_hist_bin_stat_errors,
                     label="Unfolded (#tau = %.3g) (stat err)" % (unfolder.tau),
                     line_color=PLOT_COLOURS['unfolded_stat_colour'], line_width=LW, line_style=1,
                     marker_color=PLOT_COLOURS['unfolded_stat_colour'], marker_style=20, marker_size=0.75,
                     subplot=mc_gen_hist_bin,
                     normalise_hist=False),
    ]
    if not check_entries(entries, "plot_unfolded_unnormalised_pt_bin %d" % (ibin_pt)):
        return
    plot = Plot(entries,
                xtitle=setup.particle_title,
                ytitle="N",
                title=get_pt_bin_title(setup, region, unfolder.pt_bin_edges_gen[ibin_pt], unfolder.pt_bin_edges_gen[ibin_pt+1]),
                has_data=setup.has_data,
                **PT_BIN_PLOT_ARGS)
    _modify_plot(plot)
    plot.plot("NOSTACK E1")
    plot.save("%s/unfolded_unnormalised_%s_bin_%d_divBinWidth.%s" % (setup.output_dir, setup.append, ibin_pt, setup.output_fmt))


def plot_unfolded_normalised_pt_bin(ibin_pt, unfolder, hist_bin_chopper, setup):
    mc_gen_hist_bin = hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin_pt)
    unfolded_hist_bin_stat_errors = hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', ibin_pt)
    unfolded_hist_bin_total_errors = hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin_pt)

    entries = [
        Contribution(mc_gen_hist_bin,
                     label="Generator (%s)" % (region['mc_label']),
                     line_color=PLOT_COLOURS['gen_colour'], line_width=LW,
                     marker_color=PLOT_COLOURS['gen_colour'], marker_size=0,
                     normalise_hist=False),
        Contribution(unfolded_hist_bin_total_errors,
                     label="Unfolded (#tau = %.3g) (total err)" % (unfolder.tau),
                     line_color=PLOT_COLOURS['unfolded_total_colour'], line_width=LW, line_style=1,
                     marker_color=PLOT_COLOURS['unfolded_total_colour'],# marker_style=20, marker_size=0.75,
                     subplot=mc_gen_hist_bin,
                     normalise_hist=False),
        Contribution(unfolded_hist_bin_stat_errors,
                     label="Unfolded (#tau = %.3g) (stat err)" % (unfolder.tau),
                     line_color=PLOT_COLOURS['unfolded_stat_colour'], line_width=LW, line_style=1,
                     marker_color=PLOT_COLOURS['unfolded_stat_colour'], marker_style=20, marker_size=0.75,
                     subplot=mc_gen_hist_bin,
                     normalise_hist=False),
    ]
    if not check_entries(entries, "plot_unfolded_normalised_pt_bin %d" % (ibin_pt)):
        return

    plot = Plot(entries,
                xtitle=setup.particle_title,
                ytitle=setup.pt_bin_normalised_differential_label,
                title=get_pt_bin_title(setup, region, unfolder.pt_bin_edges_gen[ibin_pt], unfolder.pt_bin_edges_gen[ibin_pt+1]),
                has_data=setup.has_data,
                **PT_BIN_PLOT_ARGS)
    _modify_plot(plot)
    plot.plot("NOSTACK E1")
    plot.save("%s/unfolded_%s_bin_%d_divBinWidth.%s" % (setup.output_dir, setup.append, ibin_pt, setup.output_fmt))


def plot_unfolded_with_alt_truth_normalised_pt_bin(ibin_pt, unfolder, alt_truth, hist_bin_chopper, setup):
    mc_gen_hist_bin = hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin_pt)
    unfolded_hist_bin_stat_errors = hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded_stat_err', ibin_pt)
    unfolded_hist_bin_total_errors = hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin_pt)

    hist_bin_chopper.add_obj('alt_hist_truth', alt_truth)
    alt_mc_gen_hist_bin = hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_hist_truth', ibin_pt)

    entries = [
        Contribution(mc_gen_hist_bin,
                     label="Generator (%s)" % (region['mc_label']),
                     line_color=PLOT_COLOURS['gen_colour'], line_width=LW,
                     marker_color=PLOT_COLOURS['gen_colour'], marker_size=0,
                     normalise_hist=False),
        Contribution(alt_mc_gen_hist_bin,
                     label="Generator (%s)" % (region['alt_mc_label']),
                     line_color=PLOT_COLOURS['alt_gen_colour'], line_width=LW, line_style=2,
                     marker_color=PLOT_COLOURS['alt_gen_colour'], marker_size=0,
                     subplot=mc_gen_hist_bin,
                     normalise_hist=False),
        Contribution(unfolded_hist_bin_total_errors,
                     label="Unfolded (#tau = %.3g) (total err)" % (unfolder.tau),
                     line_color=PLOT_COLOURS['unfolded_total_colour'], line_width=LW, line_style=1,
                     marker_color=PLOT_COLOURS['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                     subplot=mc_gen_hist_bin,
                     normalise_hist=False),
        Contribution(unfolded_hist_bin_stat_errors,
                     label="Unfolded (#tau = %.3g) (stat err)" % (unfolder.tau),
                     line_color=PLOT_COLOURS['unfolded_stat_colour'], line_width=LW, line_style=1,
                     marker_color=PLOT_COLOURS['unfolded_stat_colour'], marker_style=20, marker_size=0.75,
                     subplot=mc_gen_hist_bin,
                     normalise_hist=False),
    ]
    if not check_entries(entries, "plot_unfolded_with_alt_truth_normalised_pt_bin %d" % (ibin_pt)):
        return
    plot = Plot(entries,
                xtitle=setup.particle_title,
                ytitle=setup.pt_bin_normalised_differential_label,
                title=get_pt_bin_title(setup, region, unfolder.pt_bin_edges_gen[ibin_pt], unfolder.pt_bin_edges_gen[ibin_pt+1]),
                has_data=setup.has_data,
                **PT_BIN_PLOT_ARGS)
    _modify_plot(plot)
    plot.plot("NOSTACK E1")
    plot.save("%s/unfolded_%s_alt_truth_bin_%d_divBinWidth.%s" % (setup.output_dir, setup.append, ibin_pt, setup.output_fmt))


def plot_unfolded_with_unreg_normalised_pt_bin(ibin_pt, unfolder, unreg_unfolder, hist_bin_chopper, setup):
    hist_bin_chopper.add_obj("unreg_unfolded", unreg_unfolded.unfolded)
    mc_gen_hist_bin = hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin_pt)
    unreg_unfolded_hist_bin_total_errors = hist_bin_chopper.get_pt_bin_normed_div_bin_width('unreg_unfolded', ibin_pt)
    unfolded_hist_bin_total_errors = hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin_pt)

    entries = [
        Contribution(mc_gen_hist_bin,
                     label="Generator (%s)" % (region['mc_label']),
                     line_color=PLOT_COLOURS['gen_colour'], line_width=LW,
                     marker_color=PLOT_COLOURS['gen_colour'], marker_size=0,
                     normalise_hist=False),
        Contribution(unfolded_hist_bin_total_errors,
                     label="Unfolded (#tau = %.3g) (total err)" % (unfolder.tau),
                     line_color=PLOT_COLOURS['unfolded_total_colour'], line_width=LW, line_style=1,
                     marker_color=PLOT_COLOURS['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                     subplot=mc_gen_hist_bin,
                     normalise_hist=False),
        Contribution(unreg_unfolded_hist_bin_total_errors,
                     label="Unfolded (#tau = 0) (total err)",
                     line_color=PLOT_COLOURS['unfolded_unreg_colour'], line_width=LW, line_style=1,
                     marker_color=PLOT_COLOURS['unfolded_unreg_colour'], #marker_style=20, marker_size=0.75,
                     subplot=mc_gen_hist_bin,
                     normalise_hist=False),
    ]
    if not check_entries(entries, "plot_unfolded_with_unreg_normalised_pt_bin %d" % (ibin_pt)):
        return
    plot = Plot(entries,
                xtitle=setup.particle_title,
                ytitle=setup.pt_bin_normalised_differential_label,
                title=get_pt_bin_title(setup, region, unfolder.pt_bin_edges_gen[ibin_pt], unfolder.pt_bin_edges_gen[ibin_pt+1]),
                has_data=setup.has_data,
                **PT_BIN_PLOT_ARGS)
    _modify_plot(plot)
    plot.plot("NOSTACK E1")
    plot.save("%s/unfolded_%s_with_unreg_bin_%d_divBinWidth.%s" % (setup.output_dir, setup.append, ibin_pt, setup.output_fmt))


def plot_unfolded_with_alt_response_normalised_pt_bin(ibin_pt, unfolder, alt_unfolder, hist_bin_chopper, setup):
    hist_bin_chopper.add_obj("alt_unfolded", alt_unfolder.unfolded)
    hist_bin_chopper.add_obj("alt_hist_truth", alt_unfolder.hist_truth)

    mc_gen_hist_bin = hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin_pt)
    unfolded_hist_bin_total_errors = hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin_pt)
    alt_unfolded_hist_bin_total_errors = hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_unfolded', ibin_pt)
    alt_mc_gen_hist_gin = hist_bin_chopper.get_pt_bin_normed_div_bin_width('alt_hist_truth', ibin_pt)

    entries = [
        Contribution(mc_gen_hist_bin,
                     label="Generator (%s)" % (region['mc_label']),
                     line_color=PLOT_COLOURS['gen_colour'], line_width=LW,
                     marker_color=PLOT_COLOURS['gen_colour'], marker_size=0,
                     normalise_hist=False),
        # Contribution(alt_mc_gen_hist_bin,
        #              label="Generator (%s)" % (region['alt_mc_label']),
        #              line_color=alt_gen_colour, line_width=LW, line_style=2,
        #              marker_color=alt_gen_colour, marker_size=0,
        #              normalise_hist=False),
        Contribution(unfolded_hist_bin_total_errors,
                     label="Unfolded (#tau = %.3g) (total err)\n(%s response matrix)" % (unfolder.tau, region['mc_label']),
                     line_color=PLOT_COLOURS['unfolded_total_colour'], line_width=LW, line_style=1,
                     marker_color=PLOT_COLOURS['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                     subplot=mc_gen_hist_bin,
                     normalise_hist=False),
        Contribution(alt_unfolded_hist_bin_total_errors,
                     label="Unfolded (#tau = %.3g) (total err)\n(%s response matrix)" % (alt_unfolder.tau, region['alt_mc_label']),
                     line_color=PLOT_COLOURS['alt_unfolded_colour'], line_width=LW, line_style=1,
                     marker_color=PLOT_COLOURS['alt_unfolded_colour'], #marker_style=20, marker_size=0.75,
                     subplot=mc_gen_hist_bin,
                     normalise_hist=False),
    ]
    if not check_entries(entries, "plot_unfolded_with_unreg_normalised_pt_bin %d" % (ibin_pt)):
        return
    plot = Plot(entries,
                xtitle=setup.particle_title,
                ytitle=setup.pt_bin_normalised_differential_label,
                title=get_pt_bin_title(setup, region, unfolder.pt_bin_edges_gen[ibin_pt], unfolder.pt_bin_edges_gen[ibin_pt+1]),
                has_data=setup.has_data,
                **PT_BIN_PLOT_ARGS)
    _modify_plot(plot)
    plot.plot("NOSTACK E1")
    plot.save("%s/unfolded_%s_with_alt_response_bin_%d_divBinWidth.%s" % (setup.output_dir, setup.append, ibin_pt, setup.output_fmt))


def plot_unfolded_with_model_systs_normalised_pt_bin(ibin_pt, unfolder, model_systs, hist_bin_chopper, setup):
    syst_entries = []
    for syst_dict in model_systs:
        syst_unfolder = syst_dict['unfolder']
        syst_label = syst_dict['label']
        syst_label_no_spaces = syst_dict['label'].replace(" ", "_")

        hbc.add_obj('model_syst_%s_unfolded' % (syst_label_no_spaces), syst_unfolder.unfolded)
        hbc.add_obj('model_syst_%s_hist_truth' % (syst_label_no_spaces), syst_unfolder.hist_truth)

        syst_unfolded_hist_bin = hist_bin_chopper.get_pt_bin_normed_div_bin_width('model_syst_%s_unfolded' % (syst_label_no_spaces), ibin_pt)
        syst_gen_hist_bin = hist_bin_chopper.get_pt_bin_normed_div_bin_width('model_syst_%s_hist_truth' % (syst_label_no_spaces), ibin_pt)

        syst_entries.extend([
            Contribution(syst_unfolded_hist_bin,
                         label="Unfolded (#tau = %.3g) (total err) (%s)" % (syst_unfolder.tau, syst_label),
                         line_color=syst_dict['colour'], line_width=LW, line_style=1,
                         marker_color=syst_dict['colour'], marker_size=0,
                         subplot=syst_gen_hist_bin,
                         normalise_hist=True),
            Contribution(syst_gen_hist_bin,
                         label="Generator (%s)" % (syst_label),
                         line_color=syst_dict['colour'], line_width=LW, line_style=2,
                         marker_color=syst_dict['colour'], marker_size=0,
                         normalise_hist=True),
        ])

    # add nominal ones last
    mc_gen_hist_bin = hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin_pt)
    unfolded_hist_bin_total_errors = hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin_pt)

    syst_entries.extend([
        Contribution(mc_gen_hist_bin,
                     label="Generator (%s)" % (region['mc_label']),
                     line_color=PLOT_COLOURS['gen_colour'], line_width=LW,
                     marker_color=PLOT_COLOURS['gen_colour'], marker_size=0,
                     normalise_hist=False),
        Contribution(unfolded_hist_bin_total_errors,
                     label="Unfolded (#tau = %.3g) (total err)" % (unfolder.tau),
                     line_color=PLOT_COLOURS['unfolded_total_colour'], line_width=LW, line_style=1,
                     marker_color=PLOT_COLOURS['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                     subplot=mc_gen_hist_bin,
                     normalise_hist=False),
    ])
    if not check_entries(syst_entries, "plot_unfolded_with_model_systs_normalised_pt_bin %d" % (ibin_pt)):
        return
    plot = Plot(syst_entries,
                xtitle=setup.particle_title,
                ytitle=setup.pt_bin_normalised_differential_label,
                title=get_pt_bin_title(setup, region, unfolder.pt_bin_edges_gen[ibin_pt], unfolder.pt_bin_edges_gen[ibin_pt+1]),
                has_data=setup.has_data,
                **PT_BIN_PLOT_ARGS)
    plot.legend.SetX1(0.55)
    plot.legend.SetY1(0.72)
    plot.legend.SetX2(0.98)
    plot.legend.SetY2(0.88)
    if len(syst_entries) > 5:
        plot.legend.SetNColumns(2)
    plot.plot("NOSTACK E1")
    plot.save("%s/unfolded_%s_syst_model_bin_%d_divBinWidth.%s" % (setup.output_dir, setup.append, ibin_pt, setup.output_fmt))


def plot_unfolded_with_pdf_systs_normalised_pt_bin(ibin_pt, unfolder, pdf_systs, hist_bin_chopper, setup):
    pdf_entries = []
    for pdf_dict in pdf_systs:
        pdf_unfolder = pdf_dict['unfolder']
        pdf_label = pdf_dict['label']
        pdf_label_no_spaces = pdf_dict['label'].replace(" ", "_")

        hbc.add_obj('pdf_syst_%s_unfolded' % (pdf_label_no_spaces), pdf_unfolder.unfolded)
        hbc.add_obj('pdf_syst_%s_hist_truth' % (pdf_label_no_spaces), pdf_unfolder.hist_truth)

        pdf_unfolded_hist_bin = hist_bin_chopper.get_pt_bin_normed_div_bin_width('pdf_syst_%s_unfolded' % (pdf_label_no_spaces), ibin_pt)
        pdf_gen_hist_bin = hist_bin_chopper.get_pt_bin_normed_div_bin_width('pdf_syst_%s_hist_truth' % (pdf_label_no_spaces), ibin_pt)

        pdf_entries.extend([
            Contribution(pdf_unfolded_hist_bin,
                         label="Unfolded (#tau = %.3g) (total err) (%s)" % (pdf_unfolder.tau, pdf_label),
                         line_color=pdf_dict['colour'], line_width=LW, line_style=1,
                         marker_color=pdf_dict['colour'], marker_size=0,
                         subplot=pdf_gen_hist_bin,
                         normalise_hist=True),
            Contribution(pdf_gen_hist_bin,
                         label="Generator (%s)" % (pdf_label),
                         line_color=pdf_dict['colour'], line_width=LW, line_style=2,
                         marker_color=pdf_dict['colour'], marker_size=0,
                         normalise_hist=True),
        ])

    # add nominal ones last
    mc_gen_hist_bin = hist_bin_chopper.get_pt_bin_normed_div_bin_width('hist_truth', ibin_pt)
    unfolded_hist_bin_total_errors = hist_bin_chopper.get_pt_bin_normed_div_bin_width('unfolded', ibin_pt)

    pdf_entries.extend([
        Contribution(mc_gen_hist_bin,
                     label="Generator (%s)" % (region['mc_label']),
                     line_color=PLOT_COLOURS['gen_colour'], line_width=LW,
                     marker_color=PLOT_COLOURS['gen_colour'], marker_size=0,
                     normalise_hist=False),
        Contribution(unfolded_hist_bin_total_errors,
                     label="Unfolded (#tau = %.3g) (total err)" % (unfolder.tau),
                     line_color=PLOT_COLOURS['unfolded_total_colour'], line_width=LW, line_style=1,
                     marker_color=PLOT_COLOURS['unfolded_total_colour'], #marker_style=20, marker_size=0.75,
                     subplot=mc_gen_hist_bin,
                     normalise_hist=False),
    ])
    if not check_entries(pdf_entries, "plot_unfolded_with_model_systs_normalised_pt_bin %d" % (ibin_pt)):
        return
    plot = Plot(pdf_entries,
                xtitle=setup.particle_title,
                ytitle=setup.pt_bin_normalised_differential_label,
                title=get_pt_bin_title(setup, region, unfolder.pt_bin_edges_gen[ibin_pt], unfolder.pt_bin_edges_gen[ibin_pt+1]),
                has_data=setup.has_data,
                **PT_BIN_PLOT_ARGS)
    plot.legend.SetX1(0.55)
    plot.legend.SetY1(0.72)
    plot.legend.SetX2(0.98)
    plot.legend.SetY2(0.88)
    if len(pdf_entries) > 5:
        plot.legend.SetNColumns(2)
    plot.plot("NOSTACK E1")
    plot.save("%s/unfolded_%s_pdf_model_bin_%d_divBinWidth.%s" % (setup.output_dir, setup.append, ibin_pt, setup.output_fmt))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",
                        help="Source directory (should be the one made by unfolding.py")
    parser.add_argument("--angles",
                        choices=list(qgc.VAR_UNFOLD_DICT.keys()) + ["all"],
                        nargs='+',
                        help="Lambda angles to unfold, or 'all' for all of them (default).",
                        default=["all"])
    parser.add_argument("--outputDir",
                        default=None,
                        help='Output directory (default is the source dir')

    region_group = parser.add_argument_group('Region selection')
    region_group.add_argument("--doDijetCentral",
                              action='store_true',
                              help='Do unfolding for dijet (central) jets')
    region_group.add_argument("--doDijetForward",
                              action='store_true',
                              help='Do unfolding for dijet (forward) jets')
    region_group.add_argument("--doDijetCentralGroomed",
                              action='store_true',
                              help='Do unfolding for groomed dijet (central) jets')
    region_group.add_argument("--doDijetForwardGroomed",
                              action='store_true',
                              help='Do unfolding for groomed dijet (forward) jets')
    region_group.add_argument("--doZPJ",
                              action='store_true',
                              help='Do unfolding for Z+jet jets')
    region_group.add_argument("--doZPJGroomed",
                              action='store_true',
                              help='Do unfolding for groomed Z+jet jets')


    args = parser.parse_args()

    regions = setup_regions(args)

    if args.angles[0] == "all":
        angles = qgc.COMMON_VARS
    else:
        angles = [a for a in qgc.COMMON_VARS if a.var in args.angles]

    jet_algo = "AK4 PF PUPPI"
    if "ak8puppi" in args.source:
        jet_algo = "AK8 PF PUPPI"

    has_data = not ('_MC_all' in args.source or '_MC_split' in args.source)

    # Store all things for final summary plots
    all_binned_hists = []

    # Iterate through regions & variables
    for region in regions:
        region_dir = os.path.join(args.source, region['name'])
        if not os.path.isdir(region_dir):
            print("! Warning ! cannot find region dir", region_dir, '- skipping region')
            continue

        for angle in angles:
            angle_output_dir = "%s/%s/%s" % (args.source, region['name'], angle.var)
            if not os.path.isdir(angle_output_dir):
                print("! Warning ! cannot find angle dir", angle_output_dir, '- skipping angle', angle.var)
                continue

            # TODO: put this in a method / class
            # Check if ROOT file exists
            root_filename = os.path.join(angle_output_dir, "unfolding_result.root")
            if not os.path.isfile(root_filename):
                print("! Warning ! cannot fine unfolding ROOT file", root_filename, ' - skipping angle')

            append = "%s_%s" % (region['name'], angle.var)  # common str to put on filenames, etc. don't need angle_prepend as 'groomed' in region name
            print("*"*120)
            print("Region/var: %s" % (append))
            print("*"*120)

            input_tfile = cu.TFileCacher(root_filename)  # keep this here otherwise crashes
            unpack_dict = unpack_unfolding_root_file(input_tfile, region, angle)
            unfolder = unpack_dict['unfolder']
            alt_unfolder = unpack_dict['alt_unfolder']
            alt_hist_truth = unpack_dict['alt_hist_truth']

            # Note that experimental systs are only different response matrices, and are stored in the main unfolder
            # we have
            # - unfolder
            # - unreg_unfolder
            # - alt_unfolder
            # - model_systematics[name] unfolder
            # - pdf[ind] unfolder

            # MAKE ALL THE PLOTS
            # ------------------------------------------------------------------

            # Big 1D plots to compare things
            hbc = HistBinChopper(unfolder)
            hbc.add_obj("unfolded", unfolder.unfolded)
            hbc.add_obj("unfolded_stat_err", unfolder.unfolded_stat_err)
            hbc.add_obj("hist_truth", unfolder.hist_truth)

            setup = Setup(jet_algo=jet_algo,
                          region=region,
                          angle=angle,
                          output_dir=angle_output_dir,
                          has_data=has_data)

            # Iterate through pt bins - gen binning
            for ibin_pt in range(0, len(unfolder.pt_bin_edges_gen)-1):
                print("Doing gen pt bin", ibin_pt, "=", unfolder.pt_bin_edges_gen[ibin_pt], "->", unfolder.pt_bin_edges_gen[ibin_pt+1])
                plot_unfolded_unnormalised_pt_bin(ibin_pt=ibin_pt,
                                                  unfolder=unfolder,
                                                  hist_bin_chopper=hbc,
                                                  setup=setup)
                plot_unfolded_normalised_pt_bin(ibin_pt=ibin_pt,
                                                unfolder=unfolder,
                                                hist_bin_chopper=hbc,
                                                setup=setup)

                if alt_hist_truth:
                    plot_unfolded_with_alt_truth_normalised_pt_bin(ibin_pt=ibin_pt,
                                                                   unfolder=unfolder,
                                                                   alt_truth=alt_hist_truth,
                                                                   hist_bin_chopper=hbc,
                                                                   setup=setup)

                if unfolder.tau > 0 and unreg_unfolder:
                    plot_unfolded_with_unreg_normalised_pt_bin(ibin_pt=ibin_pt,
                                                               unfolder=unfolder,
                                                               unreg_unfolder=unreg_unfolder,
                                                               hist_bin_chopper=hbc,
                                                               setup=setup)

                if alt_unfolder:
                    plot_unfolded_with_alt_response_normalised_pt_bin(ibin_pt=ibin_pt,
                                                                      unfolder=unfolder,
                                                                      alt_unfolder=alt_unfolder,
                                                                      hist_bin_chopper=hbc,
                                                                      setup=setup)
                if len(region['model_systematics']) > 0:
                    print(region['model_systematics'])
                    plot_unfolded_with_model_systs_normalised_pt_bin(ibin_pt=ibin_pt,
                                                                     unfolder=unfolder,
                                                                     model_systs=region['model_systematics'],
                                                                     hist_bin_chopper=hbc,
                                                                     setup=setup)

            # Iterate through lambda bins - gen binning
            for ibin_pt in range(0, len(unfolder.variable_bin_edges_gen)-1):
                pass

            # Iterate through pt bins - reco binning
            for ibin_pt in range(0, len(unfolder.pt_bin_edges_reco)-1):
                pass

            # Iterate through lambda bins - reco binning
            for ibin_pt in range(0, len(unfolder.variable_bin_edges_reco)-1):
                pass
