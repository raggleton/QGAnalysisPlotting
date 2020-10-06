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
from comparator import Contribution, Plot, grab_obj
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

    # QG variable plots
    pt_bins = qgc.PT_BINS
    # pt_bins = qgc.THEORY_PT_BINS
    var_list = qgc.COMMON_VARS
    var_prepend = ""

    plot_dirname = "Dijet_QG_tighter"
    # plot_dirname = "Dijet_QG"
    plot_dirname = "ZPlusJets_QG"

    # setup TFiles here to cache things
    main_filename = qgc.QCD_FILENAME
    main_filename = qgc.DY_FILENAME
    nominal_tfile = cu.TFileCacher(os.path.join(nominal_dir, main_filename))
    herwig_filename = qgc.QCD_HERWIG_FILENAME
    herwig_filename = qgc.DY_HERWIG_FILENAME
    herwig_tfile = cu.TFileCacher(os.path.join(herwig_dir, herwig_filename)) if herwig_dir else None
    mg_herwig_tfile = cu.TFileCacher(os.path.join(herwig_dir, qgc.DY_MG_HERWIG_FILENAME)) if herwig_dir else None
    pythia_tfile = cu.TFileCacher(os.path.join(nominal_dir, qgc.QCD_PYTHIA_ONLY_FILENAME))
    # data_tfile = cu.TFileCacher(os.path.join(nominal_dir, qgc.JETHT_ZB_FILENAME))
    data_tfile = cu.TFileCacher(os.path.join(nominal_dir, qgc.SINGLE_MU_FILENAME))

    syst_up_tfiles = [cu.TFileCacher(os.path.join(syst_dir, main_filename)) for (syst_dir, _, _) in syst_ups]
    syst_down_tfiles = [cu.TFileCacher(os.path.join(syst_dir, main_filename)) for (syst_dir, _, _) in syst_downs]

    for ang in var_list:
        if ang.var not in qgc.ANGLE_REBIN_DICT:
            continue

        v = "%s%s_vs_pt" % (var_prepend, ang.var)

        for pt_ind, (start_val, end_val) in enumerate(pt_bins[:]):
            dijet_entries = []
            # Get all plots
            lw = 2

            ####################
            # DIJET REGION
            ####################

            # NOMINAL MC
            h2d_qcd_mc = nominal_tfile.Get("%s/%s" % (plot_dirname, v))
            qcd_kwargs_mc = dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                                 marker_color=qgc.QCD_COLOUR, marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=0,
                                 label="MG5+Pythia8")
            nominal_hist = qgp.get_projection_plot(h2d_qcd_mc, start_val, end_val)
            dijet_entries.append((nominal_hist, qcd_kwargs_mc))

            # HERWIG MC
            if herwig_dir is not None:
                h2d_herwig_qcd_mc = herwig_tfile.Get("%s/%s" % (plot_dirname, v))
                colh = 797
                colh = ROOT.kBlue-3
                qcd_herwig_kwargs_mc = dict(line_color=colh, line_width=lw, fill_color=colh,
                                     marker_color=colh, marker_style=cu.Marker.get('square'), marker_size=0,
                                     label="Herwig++",
                                     subplot=nominal_hist)
                dijet_entries.append((qgp.get_projection_plot(h2d_herwig_qcd_mc, start_val, end_val), qcd_herwig_kwargs_mc))

                h2d_mg_herwig_qcd_mc = mg_herwig_tfile.Get("%s/%s" % (plot_dirname, v))
                colh = ROOT.kGreen+2
                qcd_mg_herwig_kwargs_mc = dict(line_color=colh, line_width=lw, fill_color=colh,
                                     marker_color=colh, marker_style=cu.Marker.get('square'), marker_size=0,
                                     label="MG+Herwig++",
                                     subplot=nominal_hist)
                dijet_entries.append((qgp.get_projection_plot(h2d_mg_herwig_qcd_mc, start_val, end_val), qcd_mg_herwig_kwargs_mc))

            # h2d_pythia_qcd_mc = pythia_tfile.Get("%s/%s" % (plot_dirname, v))
            # colp = ROOT.kRed
            # qcd_pythia_kwargs_mc = dict(line_color=colp, line_width=lw, fill_color=colp,
            #                      marker_color=colp, marker_style=cu.Marker.get('square'), marker_size=0,
            #                      label="Pythia",
            #                      subplot=nominal_hist)
            # dijet_entries.append((qgp.get_projection_plot(h2d_pythia_qcd_mc, start_val, end_val), qcd_pythia_kwargs_mc))

            # SHIFT UP
            if syst_ups:
                for ind, (syst_up_dir, syst_up_label, syst_up_col) in enumerate(syst_ups):
                    h2d_qcd_mc2 = syst_up_tfiles[ind].Get("%s/%s" % (plot_dirname, v))
                    qcd_kwargs_mc2 = dict(line_color=syst_up_col, line_width=lw, fill_color=syst_up_col, line_style=2,
                                          marker_color=syst_up_col, marker_style=cu.Marker.get('triangleUp'), marker_size=0,
                                          label=syst_up_label,
                                          subplot=nominal_hist)
                    dijet_entries.append((qgp.get_projection_plot(h2d_qcd_mc2, start_val, end_val), qcd_kwargs_mc2))

            # SHIFT DOWN
            if syst_downs:
                for ind, (syst_down_dir, syst_down_label, syst_down_col) in enumerate(syst_downs):
                    h2d_qcd_mc3 = syst_down_tfiles[ind].Get("%s/%s" % (plot_dirname, v))
                    qcd_kwargs_mc3 = dict(line_color=syst_down_col, line_width=lw, fill_color=syst_down_col, line_style=3,
                                          marker_color=syst_down_col, marker_style=cu.Marker.get('triangleDown'), marker_size=0,
                                          label=syst_down_label,
                                          subplot=nominal_hist)
                    dijet_entries.append((qgp.get_projection_plot(h2d_qcd_mc3, start_val, end_val), qcd_kwargs_mc3))

            # JETHT-ZB DATA
            h2d_qcd_data = data_tfile.Get("%s/%s" % (plot_dirname, v))
            qcd_kwargs_data = dict(line_color=ROOT.kBlack, line_width=lw, fill_color=ROOT.kBlack,
                                   marker_color=ROOT.kBlack, marker_style=cu.Marker.get('circle'), marker_size=1,
                                   label="Data", subplot=nominal_hist)
            data_hist = qgp.get_projection_plot(h2d_qcd_data, start_val, end_val)
            dijet_entries.append((data_hist, qcd_kwargs_data))

            # rebin here only affects constant bin width plots
            # other set of plots uses binning structure defined at the top
            rebin = 5
            v_lower = v.lower()
            if "multiplicity" in v_lower:
                rebin = 2
            elif "flavour" in v_lower or "thrust" in v_lower:
                rebin = 2
            elif "ptd" in v_lower:
                rebin = 5

            xlim = None
            if "width" in v_lower or "ptd" in v_lower:
                xlim = [0, 1]
            elif"thrust" in v_lower:
                xlim = [0, 0.5]
            elif "multiplicity" in v_lower and "ak4" in nominal_dir.lower():
                xlim = [0, 100]
                xlim = [0, 80]
                xlim = [
                    [0, 40],
                    [0, 40],
                    [0, 40],
                    [0, 50],
                    [0, 50],
                    [0, 60],
                    [0, 60],
                    [0, 75],
                    [0, 75],
                    [0, 80],
                    [0, 80],
                    [0, 90],
                    [0, 90],
                    [0, 100],
                    [0, 100]
                ][pt_ind]

            # if xlim:
            #     xlim[1] = max(xlim[1], nominal_hist.GetBinLowEdge(nominal_hist.GetNbinsX()+1))

            ylim = None
            if "flavour" in v_lower:
                ylim = (0, 1)
            elif "lha" in v_lower:
                ylim = (0, 5)
                ylim = None

            # plot_dir = os.path.join(root_dir, "plots_dy_vs_qcd_mc_vs_data")
            radius, pus = cu.get_jet_config_from_dirname(nominal_dir)
            jet_str = "AK%s PF %s" % (radius.upper(), pus.upper())
            subplot_title = "Variation / nominal"
            subplot_limits = (0.9, 1.1)

            qgp.do_comparison_plot(dijet_entries,
                                   "%s/%s/ptBinned/%s_pt%dto%d_dijet.%s" % (output_dir, plot_dirname, v, start_val, end_val, OUTPUT_FMT),
                                   rebin=rebin,
                                   draw_opt="NOSTACK HIST E",
                                   title="%d < p_{T}^{jet} < %d GeV\n%s" % (start_val, end_val, jet_str),
                                   xtitle=ang.name + " (" + ang.lambda_str + ")",
                                   xlim=xlim, ylim=ylim,
                                   subplot_type='ratio',
                                   subplot_title=subplot_title,
                                   subplot_limits=subplot_limits,
                                   has_data=False)

            # dijet_entries_rebin = []
            # this_rebins = qgc.ANGLE_REBIN_DICT[ang.var]
            # rebin_hist_norminal = None
            # for hist, kwargs in dijet_entries:
            #     rebin_hist = hist.Rebin(len(this_rebins)-1, hist.GetName()+"Rebin", array('d', this_rebins))
            #     if not rebin_hist_norminal:
            #         rebin_hist_norminal = rebin_hist
            #     new_kwargs = deepcopy(kwargs)
            #     if 'subplot' in new_kwargs:
            #         new_kwargs['subplot'] = rebin_hist_norminal
            #     dijet_entries_rebin.append((rebin_hist, new_kwargs))

            # qgp.do_comparison_plot(dijet_entries_rebin,
            #                        "%s/ptBinned/%s_pt%dto%d_dijet_rebin.%s" % (output_dir, v, start_val, end_val, OUTPUT_FMT),
            #                        title="%d < p_{T}^{jet} < %d GeV\n%s" % (start_val, end_val, jet_str),
            #                        xtitle=ang.name + " (" + ang.lambda_str + ")",
            #                        xlim=None, ylim=ylim,
            #                        subplot_type='ratio',
            #                        subplot_title=subplot_title,
            #                        subplot_limits=subplot_limits,
            #                        has_data=False)




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
