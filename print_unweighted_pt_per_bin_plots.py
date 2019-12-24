#!/usr/bin/env python


"""Make unweighted pt plots

Also make plots per HT bin, or per trigger, or per Run period
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

import hunter
# hunter.trace(module='comparator', action=hunter.CodePrinter)

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat(".3f")


# Control plot output format
OUTPUT_FMT = "pdf"


def do_mc_plot(mc_entries, hist_name, output_file, xlim=None, ylim=None):
    components = []
    hst = ROOT.THStack("hst", "")
    for ent in mc_entries:
        if 'tfile' not in ent:
            ent['tfile'] = cu.open_root_file(ent['filename'])
        ent['hist'] = cu.get_from_tfile(ent['tfile'], hist_name)
        components.append(
            Contribution(ent['hist'],
                         fill_color=ent['color'],
                         line_color=ent['color'],
                         marker_color=ent['color'],
                         marker_size=0,
                         line_width=2,
                         label=ent['label'],
                         rebin_hist=5
                        )
        )
    plot = Plot(components,
                what='hist',
                has_data=False,
                xlim=xlim,
                ylim=ylim,
                ytitle="Unweighted N" if "unweighted" in hist_name else 'N')
    # plot.y_padding_min_log = 10 if 'unweighted' in hist_name else 10
    plot.default_canvas_size = (700, 600)
    plot.legend.SetNColumns(2)
    plot.legend.SetX1(0.55)
    plot.legend.SetY1(0.7)
    plot.legend.SetY2(0.88)
    plot.plot("HISTE")
    plot.set_logx()
    plot.set_logy(do_more_labels=False)
    plot.save(output_file)

    stem, ext = os.path.splitext(output_file)
    plot.plot("HISTE NOSTACK")
    plot.set_logx()
    plot.set_logy(do_more_labels=False)
    plot.save(stem+"_nostack" + ext)


if __name__ == "__main__":
    parser = qgc.get_parser()
    args = parser.parse_args()
    input_dir = args.workdirs[0]

    # DATA plots
    # --------------------------------------------------------------------------
    # Unweighted pt, showing contributions from different triggers


    data_run_entries = []

    # QCD HT plots
    # --------------------------------------------------------------------------
    # Unweighted pt, showing contributions from each HT bin
    mc_ht_entries = [
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT50to100.root' % (input_dir),
            'label': 'HT50to100',
            'color': ROOT.kAzure,
        },
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT100to200.root' % (input_dir),
            'label': 'HT100to200',
            'color': ROOT.kOrange-2,
        },
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT200to300.root' % (input_dir),
            'label': 'HT200to300',
            'color': ROOT.kGreen+2,
        },
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT300to500.root' % (input_dir),
            'label': 'HT300to500',
            'color': ROOT.kMagenta+1,
        },
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT500to700.root' % (input_dir),
            'label': 'HT500to700',
            'color': ROOT.kCyan,
        },
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT700to1000.root' % (input_dir),
            'label': 'HT700to1000',
            'color': ROOT.kRed,
        },
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT1000to1500.root' % (input_dir),
            'label': 'HT1000to1500',
            'color': ROOT.kAzure+6,
        },
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT1500to2000.root' % (input_dir),
            'label': 'HT1500to2000',
            'color': ROOT.kOrange+3
        },
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT2000toInf.root' % (input_dir),
            'label': 'HT2000toInf',
            'color': ROOT.kGreen-7,
        },
    ]

    hist_names = [
        "Dijet_Presel/pt_jet_unweighted", "Dijet_Presel/pt_jet1_unweighted",
        "Dijet_Presel/pt_jet", "Dijet_Presel/pt_jet1",
        "Dijet_tighter/pt_jet_unweighted", "Dijet_tighter/pt_jet1_unweighted",
        "Dijet_tighter/pt_jet", "Dijet_tighter/pt_jet1",
    ]
    for hname in hist_names:
        output_filename = "%s/%s.%s" % (args.output, hname.replace("/", "-"), OUTPUT_FMT)
        do_mc_plot(mc_ht_entries,
                   hist_name=hname,
                   output_file=output_filename,
                   xlim=[30, 2000],
                   ylim=[10, 1E12])
