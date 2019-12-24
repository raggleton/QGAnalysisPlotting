#!/usr/bin/env python


"""Make unweighted pt plots for QCD, data, DYJets

Also make plots per HT bin, or per trigger, or per Run period
"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse
from copy import deepcopy

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


def do_plot(entries, output_file, hist_name=None, xlim=None, ylim=None, rebin=2, is_data=True):
    components = []
    do_unweighted = any(["unweighted" in e.get('hist_name', hist_name) for e in entries])
    for ent in entries:
        if 'tfile' not in ent:
            ent['tfile'] = cu.open_root_file(ent['filename'])
        ent['hist'] = cu.get_from_tfile(ent['tfile'], ent.get('hist_name', hist_name))
        if not do_unweighted and 'scale' in ent:
            ent['hist'].Scale(ent.get('scale', 1))
        components.append(
            Contribution(ent['hist'],
                         fill_color=ent['color'],
                         line_color=ent['color'],
                         marker_color=ent['color'],
                         marker_size=0,
                         line_width=2,
                         label=ent['label'],
                         rebin_hist=rebin
                        )
        )
    plot = Plot(components,
                what='hist',
                has_data=is_data,
                xlim=xlim,
                ylim=ylim,
                ytitle="Unweighted N" if do_unweighted else 'N')
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

    # do non-stacked version
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
    # Have to add in ZB manually
    zb_entry = {
        'filename': "%s/%s" % (input_dir, qgc.ZB_FILENAME),
        'label': 'HLT_ZeroBias',
        'color': ROOT.kMagenta-9,
        'scale': 35918219492.947 / 29048.362
    }
    jet_ht_filename = "%s/%s" % (input_dir, qgc.JETHT_FILENAME)
    jet_ht_entries = [
        {
            'filename': jet_ht_filename,
            'ind': '0',
            'label': "HLT_PFJet40",
            'color': ROOT.kAzure,
        },
        {
            'filename': jet_ht_filename,
            'ind': '1',
            'label': "HLT_PFJet60",
            'color': ROOT.kOrange-2,
        },
        {
            'filename': jet_ht_filename,
            'ind': '2',
            'label': "HLT_PFJet80",
            'color': ROOT.kGreen+2,
        },
        {
            'filename': jet_ht_filename,
            'ind': '3',
            'label': "HLT_PFJet140",
            'color': ROOT.kMagenta+1,
        },
        {
            'filename': jet_ht_filename,
            'ind': '4',
            'label': "HLT_PFJet200",
            'color': ROOT.kCyan,
        },
        {
            'filename': jet_ht_filename,
            'ind': '5',
            'label': "HLT_PFJet260",
            'color': ROOT.kRed,
        },
        {
            'filename': jet_ht_filename,
            'ind': '6',
            'label': "HLT_PFJet320",
            'color': ROOT.kAzure+6,
        },
        {
            'filename': jet_ht_filename,
            'ind': '7',
            'label': "HLT_PFJet400",
            'color': ROOT.kOrange+3
        },
        {
            'filename': jet_ht_filename,
            'ind': '8',
            'label': "HLT_PFJet450",
            'color': ROOT.kGreen-7,
        },
    ]
    zb_hist_names = ["Dijet_tighter/pt_jet1", "Dijet_tighter/pt_jet1_unweighted"]
    jet_ht_hist_names = ["Dijet_jet_hist_{ind}/pt_1", "Dijet_jet_hist_unweighted_{ind}/pt_1"]
    for zb_name, ht_name in  zip(zb_hist_names, jet_ht_hist_names):
        this_zb_entry = deepcopy(zb_entry)
        this_zb_entry['hist_name'] = zb_name
        this_data_entries = [this_zb_entry]
        for ent in jet_ht_entries:
            this_entry = deepcopy(ent)
            this_entry['hist_name'] = ht_name.format(ind=this_entry['ind'])
            this_data_entries.append(this_entry)

        output_filename = "%s/DataJetHTZB-%s.%s" % (args.output, zb_name.split("/", 1)[1], OUTPUT_FMT)
        do_plot(this_data_entries,
                output_file=output_filename,
                xlim=[30, 2E3],
                ylim=[10, 1E8] if 'unweighted' in ht_name else [1, 1E12])

    # QCD HT plots
    # --------------------------------------------------------------------------
    # Unweighted pt, showing contributions from each HT bin
    qcd_ht_entries = [
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
    # for hname in hist_names:
    #     output_filename = "%s/%s.%s" % (args.output, hname.replace("/", "-"), OUTPUT_FMT)
    #     do_plot(qcd_ht_entries,
    #             hist_name=hname,
    #             output_file=output_filename,
    #             xlim=[30, 2000],
    #             ylim=[10, 1E12],
    #             is_data=False)

    # DY HT plots
    # --------------------------------------------------------------------------
    # Unweighted pt, showing contributions from each HT bin
    dy_ht_entries = [
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_DYJetsToLL_M-50_HT-0to70.root' % (input_dir),
            'label': 'HT0to70',
            'color': ROOT.kAzure,
        },
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_DYJetsToLL_M-50_HT-70to100.root' % (input_dir),
            'label': 'HT70to100',
            'color': ROOT.kOrange-2,
        },
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_DYJetsToLL_M-50_HT-100to200.root' % (input_dir),
            'label': 'HT100to200',
            'color': ROOT.kGreen+2,
        },
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_DYJetsToLL_M-50_HT-200to400.root' % (input_dir),
            'label': 'HT200to400',
            'color': ROOT.kMagenta+1,
        },
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_DYJetsToLL_M-50_HT-400to600.root' % (input_dir),
            'label': 'HT400to600',
            'color': ROOT.kCyan,
        },
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_DYJetsToLL_M-50_HT-600to800.root' % (input_dir),
            'label': 'HT600to800',
            'color': ROOT.kRed,
        },
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_DYJetsToLL_M-50_HT-800to1200.root' % (input_dir),
            'label': 'HT800to1200',
            'color': ROOT.kAzure+6,
        },
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_DYJetsToLL_M-50_HT-1200to2500.root' % (input_dir),
            'label': 'HT1200to2500',
            'color': ROOT.kOrange+3
        },
        {
            'filename': '%s/uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_DYJetsToLL_M-50_HT-2500toInf.root' % (input_dir),
            'label': 'HT2500toInf',
            'color': ROOT.kGreen-7,
        },
    ]

    hist_names = [
        "ZPlusJets_Presel/pt_jet1_unweighted", "ZPlusJets_Presel/pt_jet1",
        "ZPlusJets/pt_jet1_unweighted", "ZPlusJets/pt_jet1",
    ]
    # for hname in hist_names:
    #     output_filename = "%s/%s.%s" % (args.output, hname.replace("/", "-"), OUTPUT_FMT)
    #     do_plot(dy_ht_entries,
    #             hist_name=hname,
    #             output_file=output_filename,
    #             xlim=[30, 2000],
    #             ylim=[0.1, 1E8],
    #             is_data=False)
