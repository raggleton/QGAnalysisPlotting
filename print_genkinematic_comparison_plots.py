#!/usr/bin/env python

"""produce plots comparing gen level kinematics"""

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


def do_1D_plot(hists,
               output_filename,
               components_styles_dicts=None,
               draw_opts="NOSTACK HISTE",
               do_ratio=True,
               logx=False,
               logy=False,
               normalise_hists=True,
               title=""):

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
             ytitle="p.d.f." if normalise_hists else "N",
             title=title,
             xlim=xlim,
             ylim=ylim,
             subplot_type="ratio" if do_ratio else None,
             subplot_title="Herwig / PY8",
             subplot=contributions[0],
             subplot_limits=(0.5, 1.5),
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
        p.set_logy()
    if logx:
        p.set_logx()

    p.save(output_filename)


def do_all_1D_plots_in_dir(directories,
                           output_dir,
                           components_styles_dicts=None,
                           draw_opts="NOSTACK HISTE",
                           do_ratio=True,
                           normalise_hists=True,
                           jet_config_str=""):
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

    for obj_name in common_list_of_obj:
        print("Doing:", obj_name)

        objs = [d.Get(obj_name).Clone(obj_name + str(uuid1())) for d in directories]

        # Ignore TH1s
        if isinstance(objs[0], (ROOT.TH1F, ROOT.TH1D, ROOT.TH1I)):
            logx = obj_name in ["pt_jet",
                                "pt_jet1",
                                "pt_jet2",
                                "pt_mumu",
                                'gen_ht',
                                'pt_mu1',
                                'pt_mu2',
                                'ptHat',
                               ]
            # logy = obj_name not in ['eta_jet1',
            #                         'eta_jet2',
            #                         'eta_mu1',
            #                         'eta_mu2',
            #                         'n_jets',
            #                         'n_mu',
            #                         'pt_jet1_z_ratio'
            #                    ]
            do_1D_plot(objs,
                       components_styles_dicts=components_styles_dicts,
                       draw_opts=draw_opts,
                       do_ratio=do_ratio,
                       normalise_hists=normalise_hists,
                       logy=False,
                       title=jet_config_str,
                       logx=logx,
                       output_filename=os.path.join(output_dir, obj_name+".%s" % (OUTPUT_FMT)))
            do_1D_plot(objs,
                       components_styles_dicts=components_styles_dicts,
                       draw_opts=draw_opts,
                       do_ratio=do_ratio,
                       normalise_hists=normalise_hists,
                       logy=True,
                       title=jet_config_str,
                       logx=logx,
                       output_filename=os.path.join(output_dir, obj_name+"_logY.%s" % (OUTPUT_FMT)))
        else:
            print("Not handing non-1D hist:", obj_name)


def do_dijet_gen_distributions(root_dir):
    """Do plots comparing different different inputs in dijet region"""
    root_files = [qgc.QCD_FILENAME, qgc.QCD_PYTHIA_ONLY_FILENAME, qgc.QCD_HERWIG_FILENAME][:]
    root_files = [cu.open_root_file(os.path.join(root_dir, r)) for r in root_files]

    directories = [cu.get_from_tfile(rf, "Dijet_tighter") for rf in root_files[:]]
    mc_col = qgc.QCD_COLOUR
    mc_col2 = qgc.QCD_COLOURS[2]
    mc_col3 = qgc.QCD_COLOURS[3]
    msize = 1
    lw = 2
    csd = [
        {"label": "QCD MC [MG+PY8]", "line_color": mc_col, "fill_color": mc_col, "marker_color": mc_col, "marker_style": 22, "fill_style": 0, "marker_size": msize, 'line_width': lw},
        {"label": "QCD MC [PY8]", "line_color": mc_col2, "fill_color": mc_col2, "marker_color": mc_col2, "marker_style": 21, "fill_style": 0, "marker_size": msize, 'line_width': lw},
        {"label": "QCD MC [H++]", "line_color": mc_col3, "fill_color": mc_col3, "marker_color": mc_col3, "marker_style": 23, "fill_style": 0, "marker_size": msize, 'line_width': lw},
    ]
    jet_config_str = qgc.extract_jet_config(root_dir)

    # Compare shapes
    do_all_1D_plots_in_dir(directories=directories,
                           output_dir=os.path.join(root_dir, "Dijet_gen_kin_comparison_normalised"),
                           components_styles_dicts=csd,
                           jet_config_str=jet_config_str,
                           normalise_hists=True)


def do_zpj_gen_distributions(root_dir):
    """Do plots comparing different different inputs in Z+jet region"""
    root_files = [qgc.DY_FILENAME, qgc.DY_HERWIG_FILENAME, qgc.DY_MG_HERWIG_FILENAME]
    root_files = [cu.open_root_file(os.path.join(root_dir, r)) for r in root_files]

    directories = [cu.get_from_tfile(rf, "ZPlusJets_gen") for rf in root_files]
    mc_col = qgc.DY_COLOUR
    mc_col2 = qgc.DY_COLOURS[2]
    mc_col3 = qgc.DY_COLOURS[3]
    msize = 1
    lw = 2
    csd = [
        {"label": "DY+Jets MC [MG+PY8]", "line_color": mc_col, "fill_color": mc_col, "marker_color": mc_col, "marker_style": 21, "fill_style": 0, "marker_size": msize, 'line_width': lw},
        {"label": "DY+Jets MC [H++]", "line_color": mc_col3, "fill_color": mc_col3, "marker_color": mc_col3, "marker_style": 23, "fill_style": 0, "marker_size": msize, 'line_width': lw},
        {"label": "DY+Jets MC [MG+H++]", "line_color": mc_col2, "fill_color": mc_col2, "marker_color": mc_col2, "marker_style": 22, "fill_style": 0, "marker_size": msize, 'line_width': lw},
    ]
    jet_config_str = qgc.extract_jet_config(root_dir)

    # Compare shapes
    do_all_1D_plots_in_dir(directories=directories,
                           output_dir=os.path.join(root_dir, "ZPlusJets_gen_kin_comparison_normalised_compare"),
                           jet_config_str=jet_config_str,
                           components_styles_dicts=csd,
                           normalise_hists=True)


if __name__ == "__main__":
    parser = qgc.get_parser()
    args = parser.parse_args()

    for workdir in args.workdirs:
            # do_dijet_gen_distributions(workdir)

            do_zpj_gen_distributions(workdir)

    sys.exit(0)
