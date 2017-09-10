#!/usr/bin/env python

"""Print plots comparing different kinematics between different flavours"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
from itertools import product

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

ROOT_DIR = "workdir_ak4chs_mgpythia"

# Control output format
OUTPUT_FMT = "pdf"


def get_list_of_obj(directory):
    key_list = directory.GetListOfKeys()
    return [x.GetName() for x in key_list]


def do_all_1D_plots_in_dir(directories, output_dir, components_styles_dicts=None):
    """
    Given a set of TDirs, loop over all 1D hists, and for each plot all TDirs on a canvas.

    components_styles_dicts should be a list of style dicts, one for each directory/contribution to a plot. 
    This should include label for this component.
    """
    if len(directories) != len(components_styles_dicts):
        raise IndexError("Need same number of style dicts and directories")

    list_of_obj = [get_list_of_obj(d) for d in directories]
    print list_of_obj[0]

    # check all have same list of plots
    if not all(x == list_of_obj[0] for x in list_of_obj):
        raise RuntimeError("Different number of object in the TDirectorys")

    for obj_name in list_of_obj[0]:
        objs = [d.Get(obj_name) for d in directories]

        # Ignore TH2s
        if not isinstance(objs[0], (ROOT.TH1F, ROOT.TH1D, ROOT.TH1I)):
            print obj_name, "is not a TH1"
            continue

        contributions = [Contribution(ob, normalise_hist=True, **csd) 
                         for ob, csd in zip(objs, components_styles_dicts)]
        p = Plot(contributions, what='hist', ytitle="p.d.f.", subplot_type="ratio", subplot=contributions[0])
        p.plot("NOSTACK HISTE")
        # EURGH FIXME
        if "pt_jet" in obj_name and "ratio" not in obj_name:
            p.set_logy()
        p.save(os.path.join(output_dir, obj_name + "." + OUTPUT_FMT))


def do_dijet_distributions():
    """Do plots comparing different jet flavs in dijet region"""
    dir_names = ["Dijet_Presel_gg", "Dijet_Presel_qg", "Dijet_Presel_qq"]
    root_file = cu.open_root_file(os.path.join(ROOT_DIR, qgc.QCD_FILENAME))
    directories = [cu.get_from_file(root_file, dn) for dn in dir_names]
    gg_col = ROOT.kRed
    qg_col = ROOT.kGreen+2
    qq_col = ROOT.kBlue
    csd = [
        {"label": "gg", "line_color": gg_col, "fill_color": gg_col, "marker_color": gg_col},
        {"label": "qg", "line_color": qg_col, "fill_color": qg_col, "marker_color": qg_col},
        {"label": "qq", "line_color": qq_col, "fill_color": qq_col, "marker_color": qq_col}
    ]
    do_all_1D_plots_in_dir(directories=directories, 
                           output_dir=os.path.join(ROOT_DIR, "Dijet_kin_comparison"),
                           components_styles_dicts=csd)


def do_zpj_distributions():
    """Do plots comparing different jet flavs in z+jets region"""
    dir_names = ["ZPlusJets_Presel_q", "ZPlusJets_Presel_g"]
    root_file = cu.open_root_file(os.path.join(ROOT_DIR, qgc.DY_FILENAME))
    directories = [cu.get_from_file(root_file, dn) for dn in dir_names]
    g_col = ROOT.kRed
    q_col = ROOT.kBlue
    csd = [
        {"label": "q", "line_color": q_col, "fill_color": q_col, "marker_color": q_col},
        {"label": "g", "line_color": g_col, "fill_color": g_col, "marker_color": g_col}
    ]
    do_all_1D_plots_in_dir(directories=directories, 
                           output_dir=os.path.join(ROOT_DIR, "ZpJ_kin_comparison"),
                           components_styles_dicts=csd)

if __name__ == "__main__":
    do_dijet_distributions()
    do_zpj_distributions()
