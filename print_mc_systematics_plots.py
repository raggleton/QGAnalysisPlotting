#!/usr/bin/env python


"""Print main QG plots, comparing systematics samples. Only for MC"""


import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
os.nice(10)
import argparse
from array import array
from copy import deepcopy

# My stuff
from comparator import Contribution, Plot
import qg_common as qgc
import qg_general_plots as qgp
import common_utils as cu

# For debugging
import sys
import tracers
# sys.settrace(tracers.trace_calls)
# sys.settrace(tracers.trace_calls_detail)
# sys.settrace(tracers.trace_calls_and_returns)

# import hunter
# hunter.trace(module='comparator')

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# Control plot output format
OUTPUT_FMT = "pdf"

TOTAL_LUMI = 35918


def do_plots(nominal_dir, herwig_dir, output_dir,
             syst_ups=None, syst_downs=None):
    """
    nominal_dir is without any systematic
    herwig_dir is the herwig directory for comparison
    syst_ups and syst_downs are lists of (dirname, label) to add to the plot
    """
    lw = 2

    # COMPARE QCD
    if os.path.isfile(os.path.join(nominal_dir, qgc.QCD_FILENAME)):
        sources = [
           {
               "root_dir": nominal_dir,
               'label': "Nominal MG5+Pythia8",
               "style": {
                   'line_style': 1,
                   'line_width': lw,
                   'line_color': qgc.QCD_COLOUR,
                   'fill_color': qgc.QCD_COLOUR,
                   'marker_color': qgc.QCD_COLOUR,
                   'marker_style': cu.Marker.get(qgc.QCD_MARKER),
                   'marker_size': 0,
               }

           }
        ]

        # Add systematics
        # SHIFT UP
        if syst_ups:
            for ind, (syst_up_dir, syst_up_label, syst_up_col) in enumerate(syst_ups):
                this_source = deepcopy(sources[0])
                this_source['root_dir'] = syst_up_dir
                this_source.update(
                    dict(style=dict(line_color=syst_up_col, line_width=lw, fill_color=syst_up_col, line_style=2,
                                     marker_color=syst_up_col, marker_style=cu.Marker.get('triangleUp'), marker_size=0),
                         label=syst_up_label,
                         subplot=0)
                )
                sources.append(this_source)

        # SHIFT DOWN
        if syst_downs:
            for ind, (syst_down_dir, syst_down_label, syst_down_col) in enumerate(syst_downs):
                this_source = deepcopy(sources[0])
                this_source['root_dir'] = syst_down_dir
                this_source.update(
                    dict(style=dict(line_color=syst_down_col, line_width=lw, fill_color=syst_down_col, line_style=3,
                                    marker_color=syst_down_col, marker_style=cu.Marker.get('triangleDown'), marker_size=0),
                         label=syst_down_label, subplot=0)
                )
                sources.append(this_source)

        subplot_lim = (0.75, 1.25)

        # qgp.do_all_exclusive_plots_comparison(sources,
        #                                       var_list=qgc.COMMON_VARS,
        #                                       pt_bins=qgc.PT_BINS,
        #                                       qcd_filename=qgc.QCD_FILENAME,
        #                                       dj_cen_dirname="Dijet_QG_central_tighter",
        #                                       dj_fwd_dirname=None,
        #                                       zpj_dirname=None,
        #                                       plot_dir=os.path.join(output_dir, "plots_qcd_compare_dijet_central"),
        #                                       subplot_type="ratio",
        #                                       subplot_title="Variation / nominal",
        #                                       subplot_limits=subplot_lim,
        #                                       do_flav_tagged=False,
        #                                       has_data=False)

        # qgp.do_all_exclusive_plots_comparison(sources,
        #                                       var_list=qgc.COMMON_VARS,
        #                                       pt_bins=qgc.PT_BINS,
        #                                       qcd_filename=qgc.QCD_FILENAME,
        #                                       dj_cen_dirname=None,
        #                                       dj_fwd_dirname="Dijet_QG_forward_tighter",
        #                                       zpj_dirname=None,
        #                                       plot_dir=os.path.join(output_dir, "plots_qcd_compare_dijet_forward"),
        #                                       subplot_type="ratio",
        #                                       subplot_title="Variation / nominal",
        #                                       subplot_limits=subplot_lim,
        #                                       do_flav_tagged=False,
        #                                       has_data=False)


        qgp.do_all_exclusive_plots_comparison(sources,
                                              var_list=qgc.COMMON_VARS,
                                              pt_bins=qgc.PT_BINS,
                                              qcd_filename=qgc.QCD_FILENAME,
                                              dj_cen_dirname=None,
                                              dj_fwd_dirname="Dijet_QG_tighter",
                                              zpj_dirname=None,
                                              plot_dir=os.path.join(output_dir, "plots_qcd_compare_dijet"),
                                              subplot_type="ratio", # will use the 1st entry by default
                                              subplot_title="Variation / nominal",
                                              subplot_limits=subplot_lim,
                                              do_flav_tagged=False,
                                              has_data=False,
                                              title=qgc.Dijet_LABEL,
                                              show_region_labels=False)


    # COMPARE DY
    if os.path.isfile(os.path.join(nominal_dir, qgc.DY_FILENAME)):
        sources = [
           {
               "root_dir": nominal_dir,
               'label': "Nominal MG5+Pythia8",
               "style": {
                   'line_style': 1,
                   'line_width': lw,
                   'line_color': qgc.DY_COLOUR,
                   'fill_color': qgc.DY_COLOUR,
                   'marker_color': qgc.DY_COLOUR,
                   'marker_style': cu.Marker.get(qgc.DY_MARKER),
                   'marker_size': 0,
               }

           }
        ]

        # Add systematics
        # SHIFT UP
        if syst_ups:
            for ind, (syst_up_dir, syst_up_label, syst_up_col) in enumerate(syst_ups):
                this_source = deepcopy(sources[0])
                this_source['root_dir'] = syst_up_dir
                this_source.update(
                    dict(style=dict(line_color=syst_up_col, line_width=lw, fill_color=syst_up_col, line_style=2,
                                     marker_color=syst_up_col, marker_style=cu.Marker.get('triangleUp'), marker_size=0),
                         label=syst_up_label,
                         subplot=0)
                )
                sources.append(this_source)

        # SHIFT DOWN
        if syst_downs:
            for ind, (syst_down_dir, syst_down_label, syst_down_col) in enumerate(syst_downs):
                this_source = deepcopy(sources[0])
                this_source['root_dir'] = syst_down_dir
                this_source.update(
                    dict(style=dict(line_color=syst_down_col, line_width=lw, fill_color=syst_down_col, line_style=3,
                                    marker_color=syst_down_col, marker_style=cu.Marker.get('triangleDown'), marker_size=0),
                         label=syst_down_label, subplot=0)
                )
                sources.append(this_source)

        subplot_lim = (0.75, 1.25)

        qgp.do_all_exclusive_plots_comparison(sources,
                                              var_list=qgc.COMMON_VARS,
                                              pt_bins=qgc.PT_BINS,
                                              qcd_filename=qgc.DY_FILENAME,
                                              dj_cen_dirname=None,
                                              dj_fwd_dirname=None,
                                              zpj_dirname="ZPlusJets_QG",
                                              plot_dir=os.path.join(output_dir, "plots_dy_compare"),
                                              subplot_type="ratio",
                                              subplot_title="Variation / nominal",
                                              subplot_limits=subplot_lim,
                                              do_flav_tagged=False,
                                              has_data=False,
                                              title=qgc.ZpJ_LABEL,
                                              show_region_labels=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--nominal",
                        help="Directory name for nominal files")

    parser.add_argument("--doHerwig",
                        help="Do nominal Herwig++ files for comparison",
                        action='store_true')

    parser.add_argument("--doChargedHadron",
                        help="Do charged hadron scale shifts",
                        action='store_true')
    parser.add_argument("--doNeutralHadron",
                        help="Do neutral hadron scale shifts",
                        action='store_true')
    parser.add_argument("--doPhoton",
                        help="Do photon scale shifts",
                        action='store_true')
    parser.add_argument("--doJes",
                        help="Do jes scale shifts",
                        action='store_true')
    parser.add_argument("--doJer",
                        help="Do jer shifts",
                        action='store_true')
    parser.add_argument("--doPileup",
                        help="Do pileup shifts",
                        action='store_true')
    parser.add_argument("--doTracking",
                        help="Do tracking shifts",
                        action='store_true')

    parser.add_argument("--doMuScales",
                        help="Do mu scale shifts",
                        action='store_true')

    parser.add_argument("--doAll",
                        help="Do all systematics",
                        action='store_true')

    parser.add_argument("--outputDir",
                        help="Directory for output file",
                        default=None)
    args = parser.parse_args()

    systs_dir = os.path.join(args.nominal, 'systematics_files')

    if args.doAll:
        all_opts = [
            "doHerwig",
            "doChargedHadron",
            "doNeutralHadron",
            "doPhoton",
            "doJes",
            "doJer",
            "doPileup",
            "doTracking",
            "doMuScales"
        ]
        for opt in all_opts:
            setattr(args, opt, True)

    herwig_dir = args.nominal if args.doHerwig else None

    if args.outputDir is None:
        args.outputDir = os.path.join(args.nominal, 'print_mc_systematics')

    print(args)

    all_up_shifts = []
    all_down_shifts = []

    if args.doHerwig:
        print('Doing Herwig')
        do_plots(args.nominal,
                 herwig_dir=herwig_dir,
                 output_dir=os.path.join(args.outputDir, 'herwig'),
                 syst_ups=[],
                 syst_downs=[]
                 )

    if args.doChargedHadron:
        print('Doing ChargedHadron')
        this_up_shifts = [(os.path.join(systs_dir, 'chargedHadronShiftUp'), "Charged hadron energy scale up", ROOT.kAzure+1)]
        this_down_shifts = [(os.path.join(systs_dir, 'chargedHadronShiftDown'), "Charged hadron energy scale down", ROOT.kAzure+1)]
        do_plots(args.nominal,
                 herwig_dir=herwig_dir,
                 output_dir=os.path.join(args.outputDir, 'chargedHadron'),
                 syst_ups=this_up_shifts,
                 syst_downs=this_down_shifts
                 )
        all_up_shifts.extend(this_up_shifts)
        all_down_shifts.extend(this_down_shifts)

    if args.doNeutralHadron:
        print('Doing NeutralHadron')
        this_up_shifts = [(os.path.join(systs_dir, 'neutralHadronShiftUp'), "Neutral hadron energy scale up", ROOT.kOrange-4)]
        this_down_shifts = [(os.path.join(systs_dir, 'neutralHadronShiftDown'), "Neutral hadron energy scale down", ROOT.kOrange-4)]
        do_plots(args.nominal,
                 herwig_dir=herwig_dir,
                 output_dir=os.path.join(args.outputDir, 'neutralHadron'),
                 syst_ups=this_up_shifts,
                 syst_downs=this_down_shifts
                 )
        all_up_shifts.extend(this_up_shifts)
        all_down_shifts.extend(this_down_shifts)

    if args.doPhoton:
        print('Doing Photon')
        this_up_shifts = [(os.path.join(systs_dir, 'photonShiftUp'), "Photon energy scale up", ROOT.kMagenta-3)]
        this_down_shifts = [(os.path.join(systs_dir, 'photonShiftDown'), "Photon energy scale down", ROOT.kMagenta-3)]
        do_plots(args.nominal,
                 herwig_dir=herwig_dir,
                 output_dir=os.path.join(args.outputDir, 'photon'),
                 syst_ups=this_up_shifts,
                 syst_downs=this_down_shifts
                 )
        all_up_shifts.extend(this_up_shifts)
        all_down_shifts.extend(this_down_shifts)

    if args.doJes:
        print('Doing Jes')
        this_up_shifts = [(os.path.join(systs_dir, 'jecsmear_directionUp'), "Jet energy scale up", ROOT.kGreen+2)]
        this_down_shifts = [(os.path.join(systs_dir, 'jecsmear_directionDown'), "Jet energy scale down", ROOT.kGreen+2)]
        do_plots(args.nominal,
                 herwig_dir=herwig_dir,
                 output_dir=os.path.join(args.outputDir, 'jes'),
                 syst_ups=this_up_shifts,
                 syst_downs=this_down_shifts
                 )
        all_up_shifts.extend(this_up_shifts)
        all_down_shifts.extend(this_down_shifts)

    if args.doJer:
        print('Doing Jer')
        this_up_shifts = [(os.path.join(systs_dir, 'jersmear_directionUp'), "Jet energy resolution up", ROOT.kOrange+3)]
        this_down_shifts = [(os.path.join(systs_dir, 'jersmear_directionDown'), "Jet energy resolution down", ROOT.kOrange+3)]
        do_plots(args.nominal,
                 herwig_dir=herwig_dir,
                 output_dir=os.path.join(args.outputDir, 'jer'),
                 syst_ups=this_up_shifts,
                 syst_downs=this_down_shifts
                 )
        all_up_shifts.extend(this_up_shifts)
        all_down_shifts.extend(this_down_shifts)

    if args.doPileup:
        print('Doing Pileup')
        this_up_shifts = [(os.path.join(systs_dir, 'pileup_directionUp'), "Pileup up", ROOT.kBlue-4)]
        this_down_shifts = [(os.path.join(systs_dir, 'pileup_directionDown'), "Pileup down", ROOT.kBlue-4)]
        do_plots(args.nominal,
                 herwig_dir=herwig_dir,
                 output_dir=os.path.join(args.outputDir, 'pileup'),
                 syst_ups=this_up_shifts,
                 syst_downs=this_down_shifts
                 )
        all_up_shifts.extend(this_up_shifts)
        all_down_shifts.extend(this_down_shifts)

    if args.doTracking:
        print('Doing Tracking')
        this_up_shifts = [(os.path.join(systs_dir, 'track_directionUp'), "Tracking up", ROOT.kMagenta+3)]
        this_down_shifts = [(os.path.join(systs_dir, 'track_directionDown'), "Tracking down", ROOT.kMagenta+3)]
        do_plots(args.nominal,
                 herwig_dir=herwig_dir,
                 output_dir=os.path.join(args.outputDir, 'tracking'),
                 syst_ups=this_up_shifts,
                 syst_downs=this_down_shifts
                 )
        all_up_shifts.extend(this_up_shifts)
        all_down_shifts.extend(this_down_shifts)

    if args.doMuScales:
        print('Doing MuScales')
        this_up_shifts = [
            (os.path.join(systs_dir, 'ScaleVariationMuRUp_ScaleVariationMuFNominal'), "#mu_{R} up", ROOT.kAzure+0),
            (os.path.join(systs_dir, 'ScaleVariationMuRNominal_ScaleVariationMuFUp'), "#mu_{F} up", ROOT.kAzure+2),
            (os.path.join(systs_dir, 'ScaleVariationMuRUp_ScaleVariationMuFUp'), "#mu_{R} & #mu_{F} up", ROOT.kAzure+5),
        ]
        this_down_shifts = [
            (os.path.join(systs_dir, 'ScaleVariationMuRDown_ScaleVariationMuFNominal'), "#mu_{R} down", ROOT.kAzure+1),
            (os.path.join(systs_dir, 'ScaleVariationMuRNominal_ScaleVariationMuFDown'), "#mu_{F} down", ROOT.kAzure+3),
            (os.path.join(systs_dir, 'ScaleVariationMuRDown_ScaleVariationMuFDown'), "#mu_{R} & #mu_{F} down", ROOT.kAzure+4),
        ]
        do_plots(args.nominal,
                 herwig_dir=herwig_dir,
                 output_dir=os.path.join(args.outputDir, 'scale'),
                 syst_ups=this_up_shifts,
                 syst_downs=this_down_shifts
                 )
        # don't add to all...too many contributions

    # do all shifts
    print("Doing all chosen")
    do_plots(args.nominal,
             herwig_dir=herwig_dir,
             output_dir=os.path.join(args.outputDir, 'all'),
             syst_ups=all_up_shifts,
             syst_downs=all_down_shifts)
