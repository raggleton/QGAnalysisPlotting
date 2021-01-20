#!/usr/bin/env python

"""Print kinematics plots, comparing samples.

Does comparison plots for data & MC, separately.

Make sure jetht and zerobias are hadded into same file
"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
os.nice(10)
from itertools import product
import numpy as np
np.seterr(all='raise')
from array import array
from uuid import uuid1
import argparse
from copy import copy, deepcopy

# My stuff
from comparator import Contribution, Plot
import qg_common as qgc
import qg_general_plots as qgg
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
ROOT.gStyle.SetOptStat(0)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

# Control output format
OUTPUT_FMT = "pdf"


def exists_in_all(filename, dirs):
    return all([os.path.isfile(os.path.join(d, filename)) for d in dirs])


def get_list_of_obj(directory):
    key_list = directory.GetListOfKeys()
    return sorted([x.GetName() for x in key_list])


def find_largest_filled_bin(hist):
    for i in range(hist.GetNbinsX(), 0, -1):
        n = hist.GetBinContent(i)
        if n != 0:
            return i, hist.GetBinLowEdge(i+1)
    return 1, hist.GetBinLowEdge(1)


def find_first_filled_bin(hist):
    for i in range(1, hist.GetNbinsX()+1):
        n = hist.GetBinContent(i)
        if n != 0:
            return i, hist.GetBinLowEdge(i)
    return 1, hist.GetBinLowEdge(1)



def do_1D_plot(hists, output_filename, components_styles_dicts=None,
               draw_opts="NOSTACK HISTE", do_ratio=True, logx=False, logy=False,
               normalise_hists=True, title=""):

    if (len(hists) != len(components_styles_dicts)):
        raise RuntimeError("# hists != # components_styles_dicts (%d vs %d)" % (len(hists), len(components_styles_dicts)))

    hists = [h.Clone(cu.get_unique_str()) for h in hists]
    contributions = [Contribution(hist, 
                                  normalise_hist=normalise_hists, 
                                  label=csd.get('label', 'x'), 
                                  **csd.get('style', {}))
                     for hist, csd in zip(hists, components_styles_dicts)]

    if len(contributions) == 0:
        return

    # Ignore if all empty objs
    total_entries = sum(c.obj.GetEntries() for c in contributions)
    if total_entries == 0:
        print("WARNING: all plots have 0 entries")
        return

    min_val = min([h.GetMinimum(0) for h in hists])
    max_val = max([h.GetMaximum() for h in hists])
    # print("Auto y limits:", min_val, max_val)
    if logy:
        ylim = [0.5*min_val, 50*max_val]
    else:
        # ylim = [0.5*min_val, 1.5*max_val]
        ylim = [0, 1.5*max_val]

    # Auto calc x limits to avoid lots of empty bins
    high_bin = max([find_largest_filled_bin(h)[0] for h in hists])
    low_bin = min([find_first_filled_bin(h)[0] for h in hists])
    xlim = [hists[0].GetBinLowEdge(low_bin-1), hists[0].GetBinLowEdge(high_bin+2)]

    p = Plot(contributions, what='hist',
             ytitle="#DeltaN/N" if normalise_hists else "N",
             title=title,
             xlim=xlim,
             ylim=ylim,
             subplot_type="ratio" if do_ratio else None,
             subplot_title="* / %s" % contributions[0].label,
             subplot=contributions[0],
             # subplot_limits=(0.5, 1.5),
             subplot_limits=(0.25, 1.75) if logy else (0.5, 1.5),
             )
    # p.legend.SetX1(0.55)
    # # p.legend.SetX2(0.95)
    # p.legend.SetY1(0.7)
    # p.legend.SetY2(0.85)
    p.legend.SetX1(0.5)
    p.legend.SetX2(0.97)
    if len(contributions) > 4:
        p.legend.SetY1(0.6)
    else:
        p.legend.SetY1(0.7)
    p.legend.SetY2(0.9)
    p.plot(draw_opts)

    if logy:
        p.set_logy(do_more_labels=False)
    if logx:
        p.set_logx(do_more_labels=False)

    # p.save(os.path.join(output_dir, obj_name+".%s" % (OUTPUT_FMT)))
    p.save(output_filename)


def do_all_1D_projection_plots_in_dir(directories,
                                      output_dir,
                                      components_styles_dicts=None,
                                      draw_opts="NOSTACK HISTE",
                                      do_ratio=True,
                                      normalise_hists=True,
                                      jet_config_str="",
                                      bin_by=None):
    """
    Given a set of TDirs, loop over all 2D hists, do projection hists for chosen bins, and plot all TDir contributions on a canvas for comparison.

    components_styles_dicts should be a list of style dicts, one for each directory/contribution to a plot.
    This should include label for this component.

    """
    if len(directories) != len(components_styles_dicts):
        raise IndexError("Need same number of style dicts and directories")

    list_of_obj = [get_list_of_obj(d) for d in directories]

    # check all have same list of plots
    if not all(x == list_of_obj[0] for x in list_of_obj):
        print("Different number of object in the TDirectorys")

    common_list_of_obj = set(list_of_obj[0])
    for l in list_of_obj[1:]:
        common_list_of_obj = common_list_of_obj & set(l)
    common_list_of_obj = sorted(list(common_list_of_obj))

    pt_bins = qgc.PT_BINS

    for obj_name in common_list_of_obj:

        if "flav" in obj_name:
            continue
        if obj_name in [
                        'eta_jet1_vs_eta_jet2',
                        'phi_jet1_vs_pt_jet1', 'phi_jet2_vs_pt_jet1',
                        # 'reliso_mu1_vs_pt_jet1', 'reliso_mu2_vs_pt_jet1',
                        # 'dphi_mumu_jet1_vs_pt_jet1',
                        # 'dphi_mumu_vs_pt_jet1',
                        'pt_jet1_z_pt_jet2_z_ratio',
                        # 'met_sig_vs_pt_jet1'
                        'weight_vs_puHat_genHT_ratio',
                        ]:
            continue

        if "genjet_ind_recojet_ind" in obj_name:
            continue

        print("Doing:", obj_name)

        objs = [d.Get(obj_name).Clone(obj_name + str(uuid1())) for d in directories]

        # Ignore TH1s
        if not isinstance(objs[0], (ROOT.TH2F, ROOT.TH2D, ROOT.TH2I)):
            logx = obj_name in ["pt_jet", "pt_jet1", "pt_jet2", "pt_mumu", 'gen_ht', 'pt_jet_response_binning', 'pt_genjet_response_binning', 'pt_jet1_unweighted', 'pt_jet_unweighted']
            rebin = 1
            xlim = None
            if obj_name in ['pt_jet', 'pt_jet1', 'pt_jet2', 'pt_mumu']:
                rebin = 10
                xlim = [30, None]
            for obj in objs:
                obj.Rebin(rebin)
            do_1D_plot(objs,
                       components_styles_dicts=components_styles_dicts,
                       draw_opts=draw_opts, do_ratio=do_ratio, normalise_hists=normalise_hists, logy=True,
                       title=jet_config_str, logx=logx,
                       output_filename=os.path.join(output_dir, obj_name+".%s" % (OUTPUT_FMT)))
        else:

            for pt_min, pt_max in pt_bins:
                # print(pt_min, pt_max)
                rebin = 1
                # exceptions...why didn't I pick the same number of bins...
                do_not_rebin = any([
                    "n_jets" in obj_name,
                    "n_mu" in obj_name,
                    "met_sig" in obj_name,
                    # obj_name.startswith('dphi_mumu'),
                    obj_name.startswith('pt_jet3_frac'),
                    obj_name.startswith('pt_jet1_jet2_ratio'),
                    obj_name.startswith('pt_jet1_z_ratio'),
                    obj_name.startswith('pt_jet2_z_ratio'),
                    obj_name.startswith('dphi_jet1_z_vs_pt_jet1'),
                    # obj_name.startswith('m_jj'),
                ])
                if not do_not_rebin:
                    # if objs[0].GetNbinsX() % 5 == 0 and objs[0].GetNbinsX() >= 100:
                    #     rebin = 5
                    if objs[0].GetNbinsX() % 4 == 0 and objs[0].GetNbinsX() >= 80:
                        rebin = 2
                    elif objs[0].GetNbinsX() % 3 == 0 and objs[0].GetNbinsX() >= 60:
                        rebin = 3
                if obj_name.startswith("m_jj"):
                    rebin = 2
                if "reliso" in obj_name:
                    rebin = 2
                hists = [qgg.get_projection_plot(ob, pt_min, pt_max).Rebin(rebin) for ob in objs]

                if bin_by == "ave":
                    title = "#splitline{%s}{%d < #LT p_{T}^{jet} #GT < %d GeV}" % (jet_config_str, pt_min, pt_max)
                elif bin_by == "Z":
                    title = "#splitline{%s}{%d < p_{T}^{Z} < %d GeV}" % (jet_config_str, pt_min, pt_max)

                logx = 'pt_jet_response' in obj_name

                do_1D_plot(hists, components_styles_dicts=components_styles_dicts,
                           draw_opts=draw_opts, do_ratio=do_ratio,
                           normalise_hists=normalise_hists,
                           logx=logx,
                           logy=False,
                           title=title,
                           output_filename=os.path.join(output_dir, obj_name+"_pt%dto%d.%s" % (pt_min, pt_max, OUTPUT_FMT)))


def do_comparison_plots(workdir_label_pairs, output_dir):
    dirnames = [w[0] for w in workdir_label_pairs]

    # templates, we'll change the filename/dir as per instance
    total_len = len(workdir_label_pairs)
    mark = cu.Marker()
    sources = [
       {
           # "root_dir": wd,
           'label': label,
           "style": {
               'line_style': 1,
               'line_color': cu.get_colour_seq(ind, total_len),
               'marker_color': cu.get_colour_seq(ind, total_len),
               'marker_style': m,
               'marker_size': 0.75,
           }

       }
       for ind, ((wd, label), m) in enumerate(zip(workdir_label_pairs, mark.cycle()))
    ]

    jet_config_str = qgc.extract_jet_config(dirnames[0])
    if len(dirnames) >1 and qgc.extract_jet_config(dirnames[1]) != jet_config_str:
        print("Conflicting jet config str, not adding")
        jet_config_str = None

    # COMPARE NOMINAL QCD
    if exists_in_all(qgc.QCD_FILENAME, dirnames):
        root_files = [cu.open_root_file(os.path.join(d, qgc.QCD_FILENAME)) for d in dirnames]
        directories = [cu.get_from_tfile(rf, "Dijet_tighter") for rf in root_files]

        this_sources = deepcopy(sources)
        for s in this_sources:
            s['label'] = "QCD [MG+PY8] [%s]" % s['label']

        do_all_1D_projection_plots_in_dir(directories=directories,
                                          components_styles_dicts=this_sources,
                                          output_dir=os.path.join(output_dir, "plots_qcd_compare_dijet_tighter_kinematics_normalised"),
                                          jet_config_str=jet_config_str,
                                          normalise_hists=True,
                                          bin_by='ave')

        directories = [cu.get_from_tfile(rf, "Dijet_eta_ordered") for rf in root_files]
        do_all_1D_projection_plots_in_dir(directories=directories,
                                          components_styles_dicts=this_sources,
                                          output_dir=os.path.join(output_dir, "plots_qcd_compare_dijet_eta_ordered_kinematics_normalised"),
                                          jet_config_str=jet_config_str,
                                          normalise_hists=True,
                                          bin_by='ave')


    # COMPARE NOMINAL DY
    if exists_in_all(qgc.DY_FILENAME, dirnames):
        root_files = [cu.open_root_file(os.path.join(d, qgc.DY_FILENAME)) for d in dirnames]
        directories = [cu.get_from_tfile(rf, "ZPlusJets") for rf in root_files]

        this_sources = deepcopy(sources)
        for s in this_sources:
            s['label'] = "Z+Jet [MG+PY8] [%s]" % s['label']

        do_all_1D_projection_plots_in_dir(directories=directories,
                                          components_styles_dicts=this_sources,
                                          output_dir=os.path.join(output_dir, "plots_dy_compare_kinematics_absolute"),
                                          jet_config_str=jet_config_str,
                                          normalise_hists=False,
                                          bin_by='Z')
        
        do_all_1D_projection_plots_in_dir(directories=directories,
                                          components_styles_dicts=this_sources,
                                          output_dir=os.path.join(output_dir, "plots_dy_compare_kinematics_normalised"),
                                          jet_config_str=jet_config_str,
                                          normalise_hists=True,
                                          bin_by='Z')

    # COMPARE JETHT+ZEROBIAS
    if exists_in_all(qgc.JETHT_ZB_FILENAME, dirnames):
        root_files = [cu.open_root_file(os.path.join(d, qgc.JETHT_ZB_FILENAME)) for d in dirnames]
        directories = [cu.get_from_tfile(rf, "Dijet_tighter") for rf in root_files]

        this_sources = deepcopy(sources)
        for s in this_sources:
            s['label'] = "Data [%s]" % s['label']

        directories = [cu.get_from_tfile(rf, "Dijet_Presel") for rf in root_files]
        do_all_1D_projection_plots_in_dir(directories=directories,
                                          components_styles_dicts=this_sources,
                                          output_dir=os.path.join(output_dir, "plots_jetht_zb_compare_dijet_presel_kinematics_normalised"),
                                          jet_config_str=jet_config_str,
                                          normalise_hists=False,
                                          bin_by='ave')

        directories = [cu.get_from_tfile(rf, "Dijet_tighter") for rf in root_files]
        do_all_1D_projection_plots_in_dir(directories=directories,
                                          components_styles_dicts=this_sources,
                                          output_dir=os.path.join(output_dir, "plots_jetht_zb_compare_dijet_tighter_kinematics_normalised"),
                                          jet_config_str=jet_config_str,
                                          normalise_hists=True,
                                          bin_by='ave')

        do_all_1D_projection_plots_in_dir(directories=directories,
                                          components_styles_dicts=this_sources,
                                          output_dir=os.path.join(output_dir, "plots_jetht_zb_compare_dijet_tighter_kinematics_abs"),
                                          jet_config_str=jet_config_str,
                                          normalise_hists=False,
                                          bin_by='ave')


    # COMPARE SINGLEMU

    # COMPARE HERWIG++ QCD
    if exists_in_all(qgc.QCD_HERWIG_FILENAME, dirnames):
        root_files = [cu.open_root_file(os.path.join(d, qgc.QCD_HERWIG_FILENAME)) for d in dirnames]
        directories = [cu.get_from_tfile(rf, "Dijet_tighter") for rf in root_files]

        this_sources = deepcopy(sources)
        for s in this_sources:
            s['label'] = "QCD [H++] [%s]" % s['label']

        # do_all_1D_projection_plots_in_dir(directories=directories,
        #                                   components_styles_dicts=this_sources,
        #                                   output_dir=os.path.join(output_dir, "plots_qcd_hpp_compare_dijet_tighter_kinematics_normalised"),
        #                                   jet_config_str=jet_config_str,
        #                                   normalise_hists=True,
        #                                   bin_by='ave')

        # directories = [cu.get_from_tfile(rf, "Dijet_eta_ordered") for rf in root_files]
        # do_all_1D_projection_plots_in_dir(directories=directories,
        #                                   components_styles_dicts=this_sources,
        #                                   output_dir=os.path.join(output_dir, "plots_qcd_hpp_compare_dijet_eta_ordered_kinematics_normalised"),
        #                                   jet_config_str=jet_config_str,
        #                                   normalise_hists=True,
        #                                   bin_by='ave')
        
        do_all_1D_projection_plots_in_dir(directories=directories,
                                          components_styles_dicts=this_sources,
                                          output_dir=os.path.join(output_dir, "plots_qcd_hpp_compare_dijet_tighter_kinematics_absolute"),
                                          jet_config_str=jet_config_str,
                                          normalise_hists=False,
                                          bin_by='ave')

        directories = [cu.get_from_tfile(rf, "Dijet_eta_ordered") for rf in root_files]
        do_all_1D_projection_plots_in_dir(directories=directories,
                                          components_styles_dicts=this_sources,
                                          output_dir=os.path.join(output_dir, "plots_qcd_hpp_compare_dijet_eta_ordered_kinematics_absolute"),
                                          jet_config_str=jet_config_str,
                                          normalise_hists=False,
                                          bin_by='ave')

    # COMPARE HERWIG++ DY
    if exists_in_all(qgc.DY_HERWIG_FILENAME, dirnames):
        root_files = [cu.open_root_file(os.path.join(d, qgc.DY_HERWIG_FILENAME)) for d in dirnames]
        directories = [cu.get_from_tfile(rf, "ZPlusJets") for rf in root_files]

        this_sources = deepcopy(sources)
        for s in this_sources:
            s['label'] = "Z+Jet [H++] [%s]" % s['label']

        do_all_1D_projection_plots_in_dir(directories=directories,
                                          components_styles_dicts=this_sources,
                                          output_dir=os.path.join(output_dir, "plots_dy_hpp_compare_kinematics_absolute"),
                                          jet_config_str=jet_config_str,
                                          normalise_hists=False,
                                          bin_by='Z')
        
        do_all_1D_projection_plots_in_dir(directories=directories,
                                          components_styles_dicts=this_sources,
                                          output_dir=os.path.join(output_dir, "plots_dy_hpp_compare_kinematics_normalised"),
                                          jet_config_str=jet_config_str,
                                          normalise_hists=True,
                                          bin_by='Z')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',
                        required=True,
                        action='append',
                        help='Workdir with ROOT files to process. '
                        'Each --workdir must have a corresponding --label')
    parser.add_argument('--label',
                        required=True,
                        action='append',
                        help='Label for workdir.')
    parser.add_argument("-o", "--output",
                        help="Directory to put output plot dirs into. Default is workdir.",
                        default=None)
    args = parser.parse_args()
    print(args)

    if args.output is None:
        args.output = args.workdir[0]
        print("No output dir specified, using", args.output)

    # for workdir in args.workdirs:
    #     do_dijet_distribution(workdir)
    #     do_zpj_distributions(workdir)

    do_comparison_plots(list(zip(args.workdir, args.label)), output_dir=args.output)
    sys.exit(0)
