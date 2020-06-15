#!/usr/bin/env python

"""produce plots comparing different kinematics between data & MC"""

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

# My stuff
from comparator import Contribution, Plot, grab_obj
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

    hists = [h.Clone(h.GetName() + str(uuid1())) for h in hists]
    contributions = [Contribution(hist, normalise_hist=normalise_hists, **csd)
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
             subplot_title="MC / Data",
             subplot=contributions[0],
             # subplot_limits=(0.5, 1.5),
             subplot_limits=(0, 2) if logy else (0.5, 1.5),
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


def do_all_1D_projection_plots_in_dir(directories, output_dir, components_styles_dicts=None,
                                      draw_opts="NOSTACK HISTE", do_ratio=True,
                                      normalise_hists=True,
                                      jet_config_str="",
                                      signal_mask=None,
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


def do_dijet_distributions(root_dir):
    """Do plots comparing different different inputs in dijet region"""
    # root_files = [qgc.ZB_FILENAME, qgc.JETHT_FILENAME, qgc.QCD_FILENAME][:]
    # root_files = [qgc.ZB_FILENAME, qgc.JETHT_FILENAME, qgc.QCD_PYTHIA_ONLY_FILENAME][:]
    # root_files = [qgc.ZB_FILENAME, qgc.JETHT_FILENAME, qgc.QCD_FILENAME][:]
    # root_files = [qgc.ZB_FILENAME, qgc.JETHT_FILENAME, qgc.QCD_FILENAME, qgc.QCD_PYTHIA_ONLY_FILENAME][:]
    root_files = [qgc.JETHT_ZB_FILENAME, qgc.QCD_FILENAME, qgc.QCD_PYTHIA_ONLY_FILENAME, qgc.QCD_HERWIG_FILENAME][:]
    # root_files = [qgc.JETHT_ZB_FILENAME, qgc.QCD_FILENAME, qgc.QCD_HERWIG_FILENAME][:]
    # root_files = [qgc.ZB_FILENAME, qgc.JETHT_FILENAME, qgc.QCD_FILENAME, qgc.ZPJ_ALL_FILENAME][:]
    root_files = [cu.open_root_file(os.path.join(root_dir, r)) for r in root_files]

    # herwig_dir = "workdir_ak4chs_herwig_newFlav_withPUreweight_withMuSF"
    # herwig_dir = "workdir_ak4chs_herwig_newFlav_withPUreweight_withMuSF_noExtraJetCuts"
    # root_files.append(cu.open_root_file(os.path.join(herwig_dir, qgc.QCD_FILENAME)))

    # directories = [cu.get_from_tfile(rf, "Dijet") for rf in root_files]
    directories = [cu.get_from_tfile(rf, "Dijet_tighter") for rf in root_files[:]]
    # directories.extend([cu.get_from_tfile(rf, "Dijet_tighter") for rf in root_files[1:]])
    mc_col = qgc.QCD_COLOUR
    mc_col2 = qgc.QCD_COLOURS[2]
    mc_col3 = qgc.QCD_COLOURS[3]
    data_col = qgc.JETHT_COLOUR
    zb_col = ROOT.kGreen+2
    msize = 0.75
    csd = [
        {"label": "Data", "line_color": data_col, "fill_color": data_col, "marker_color": data_col, "marker_style": 20, "fill_style": 0, "marker_size": msize},
        {"label": "QCD MC [MG+PY8]", "line_color": mc_col, "fill_color": mc_col, "marker_color": mc_col, "marker_style": 22, "fill_style": 0, "marker_size": msize},
        {"label": "QCD MC [PY8]", "line_color": mc_col2, "fill_color": mc_col2, "marker_color": mc_col2, "marker_style": 21, "fill_style": 0, "marker_size": msize},
        {"label": "QCD MC [H++]", "line_color": mc_col3, "fill_color": mc_col3, "marker_color": mc_col3, "marker_style": 23, "fill_style": 0, "marker_size": msize},
    ]
    jet_config_str = qgc.extract_jet_config(root_dir)

    # Compare yields
    # do_all_1D_projection_plots_in_dir(directories=directories,
    #                                   output_dir=os.path.join(root_dir, "Dijet_data_mc_kin_comparison_absolute"),
    # #                                   # output_dir=os.path.join(root_dir, "Dijet_data_mc_kin_comparison_absolute_pythiaOnly"),
    #                                   # output_dir=os.path.join(root_dir, "Dijet_data_mc_kin_comparison_absolute_both"),
    #                                   # output_dir=os.path.join(root_dir, "Dijet_data_mc_kin_comparison_absolute_everything"),
    # #                                   # output_dir=os.path.join(root_dir, "Dijet_data_mc_kin_comparison_absolute_all"),
    #                                   components_styles_dicts=csd,
    #                                   jet_config_str=jet_config_str,
    #                                   normalise_hists=False,
    #                                   bin_by='ave')

    # Compare shapes
    do_all_1D_projection_plots_in_dir(directories=directories,
                                      output_dir=os.path.join(root_dir, "Dijet_data_mc_kin_comparison_normalised"),
                                      # output_dir=os.path.join(root_dir, "Dijet_data_mc_kin_comparison_normalised_pythiaOnly"),
                                      # output_dir=os.path.join(root_dir, "Dijet_data_mc_kin_comparison_normalised_both"),
                                      # output_dir=os.path.join(root_dir, "Dijet_data_mc_kin_comparison_normalised_all"),
                                      components_styles_dicts=csd,
                                      jet_config_str=jet_config_str,
                                      bin_by='ave')


def do_zpj_distributions(root_dir):
    """Do plots comparing different different inputs in Z+jet region"""
    root_files = [qgc.SINGLE_MU_FILENAME, qgc.DY_FILENAME, qgc.DY_HERWIG_FILENAME, qgc.DY_MG_HERWIG_FILENAME]
    root_files = [cu.open_root_file(os.path.join(root_dir, r)) for r in root_files]

    directories = [cu.get_from_tfile(rf, "ZPlusJets") for rf in root_files]
    mc_col = qgc.DY_COLOUR
    mc_col2 = qgc.DY_COLOURS[2]
    mc_col3 = qgc.DY_COLOURS[3]
    data_col = qgc.SINGLE_MU_COLOUR
    msize = 0.75
    csd = [
        {"label": "Data", "line_color": data_col, "fill_color": data_col, "marker_color": data_col, "marker_style": 20, "fill_style": 0, "marker_size": msize},
        {"label": "DY+Jets MC [MG+PY8]", "line_color": mc_col, "fill_color": mc_col, "marker_color": mc_col, "marker_style": 21, "fill_style": 0, "marker_size": msize},
        {"label": "DY+Jets MC [H++]", "line_color": mc_col3, "fill_color": mc_col3, "marker_color": mc_col3, "marker_style": 23, "fill_style": 0, "marker_size": msize},
        {"label": "DY+Jets MC [MG+H++]", "line_color": mc_col2, "fill_color": mc_col2, "marker_color": mc_col2, "marker_style": 22, "fill_style": 0, "marker_size": msize},
    ]
    jet_config_str = qgc.extract_jet_config(root_dir)

    # Compare yields
    # do_all_1D_projection_plots_in_dir(directories=directories,
    #                                   output_dir=os.path.join(root_dir, "ZPlusJets_data_mc_kin_comparison_absolute"),
    #                                   # output_dir=os.path.join(root_dir, "ZPlusJets_data_mc_kin_comparison_absolute_compareKFactor"),
    #                                   # output_dir=os.path.join(root_dir, "ZPlusJets_data_mc_kin_comparison_absolute_all"),
    #                                   components_styles_dicts=csd,
    #                                   jet_config_str=jet_config_str,
    #                                   normalise_hists=False,
    #                                   bin_by='Z')

    # Compare shapes
    do_all_1D_projection_plots_in_dir(directories=directories,
                                      output_dir=os.path.join(root_dir, "ZPlusJets_data_mc_kin_comparison_normalised_compare"),
                                      # output_dir=os.path.join(root_dir, "ZPlusJets_data_mc_kin_comparison_normalised_compare_KFactor"),
                                      # output_dir=os.path.join(root_dir, "ZPlusJets_data_mc_kin_comparison_normalised_all"),
                                      jet_config_str=jet_config_str,
                                      components_styles_dicts=csd,
                                      normalise_hists=True,
                                      bin_by='Z')

    # Preselection hists
    # directories_presel = directories = [cu.get_from_tfile(rf, "ZPlusJets_Presel") for rf in root_files]
    # do_all_1D_projection_plots_in_dir(directories=directories,
    #                                   output_dir=os.path.join(root_dir, "ZPlusJets_Presel_data_mc_kin_comparison_normalised_compare"),
    #                                   # output_dir=os.path.join(root_dir, "ZPlusJets_data_mc_kin_comparison_normalised_compare_KFactor"),
    #                                   # output_dir=os.path.join(root_dir, "ZPlusJets_data_mc_kin_comparison_normalised_all"),
    #                                   jet_config_str=jet_config_str,
    #                                   components_styles_dicts=csd,
    #                                   normalise_hists=True,
    #                                   bin_by='Z')



if __name__ == "__main__":
    parser = qgc.get_parser()
    args = parser.parse_args()

    for workdir in args.workdirs:
            do_dijet_distributions(workdir)

            do_zpj_distributions(workdir)

    sys.exit(0)
