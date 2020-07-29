#!/usr/bin/env python


"""Some comaprison plots"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
from array import array
import numpy as np
import math
from itertools import product, chain
from copy import copy, deepcopy

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp

# Use rootpy to throw exceptions on ROOT errors, but need DANGER enabled
import rootpy
import rootpy.logger.magic as M; M.DANGER.enabled = True

ROOT.gErrorIgnoreLevel = ROOT.kWarning
# ROOT.gErrorIgnoreLevel = ROOT.kInfo
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)  # VERY IMPORTANT - somewhere, closing a TFile for exp systs deletes a map...dunno why


MAINDIR="/Volumes/Extreme SSD/Projects/QGAnalysis"

ak4_v2_v2_dir = os.path.join(MAINDIR, "workdir_80X_ak4puppi_data_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1")
ak8_v2_v2_dir = os.path.join(MAINDIR, "workdir_80X_ak8puppi_data_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1")

# ak4_v3_v2_dir = os.path.join(MAINDIR, "workdir_102X_data_80X_mc_ak4puppi_data_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1")
# ak8_v3_v2_dir = os.path.join(MAINDIR, "workdir_102X_data_80X_mc_ak8puppi_data_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1")

ak4_v3_v2_dir = os.path.join(MAINDIR, "workdir_102X_v3data_v2mc_ak4puppi_data_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1")
ak8_v3_v2_dir = os.path.join(MAINDIR, "workdir_102X_v3data_v2mc_ak8puppi_data_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1")

ak4_old_v2_noJER_dir = os.path.join(MAINDIR, "workdir_80X_ak4puppi_mgpythia_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_noJER")
ak8_old_v2_noJER_dir = os.path.join(MAINDIR, "workdir_80X_ak8puppi_mgpythia_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_noJER")

ak4_v3_v2_noJER_dir = os.path.join(MAINDIR, "workdir_102X_v2_ak4puppi_mgpythia_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_noJER")
ak8_v3_v2_noJER_dir = os.path.join(MAINDIR, "workdir_102X_v2_ak8puppi_mgpythia_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_noJER")

ak4_v3_v3_dir = os.path.join(MAINDIR, "workdir_102X_ak4puppi_data_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1")
ak8_v3_v3_dir = os.path.join(MAINDIR, "workdir_102X_ak8puppi_data_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1")

def do_data_mc_plot(dirname, histname, output_filename, **plot_kwargs):
    data_file = cu.open_root_file(os.path.join(dirname, qgc.JETHT_ZB_FILENAME))
    # qcd_file = cu.open_root_file(os.path.join(dirname, qgc.QCD_FILENAME))
    qcd_file = cu.open_root_file(os.path.join(dirname, qgc.QCD_PYTHIA_ONLY_FILENAME))

    data_hist = cu.get_from_tfile(data_file, histname)
    qcd_hist = cu.get_from_tfile(qcd_file, histname)
    conts = [
        Contribution(data_hist, label="Data", line_color=ROOT.kBlack),
        Contribution(qcd_hist, label="QCD MC", line_color=ROOT.kRed, subplot=data_hist)
    ]
    plot = Plot(conts,
                what='hist',
                ytitle="N",
                subplot_type="ratio",
                subplot_title="Simulation / data",
                **plot_kwargs)
    plot.y_padding_max_log = 500
    plot.legend.SetY1(0.7)
    plot.plot("NOSTACK HIST E")
    plot.set_logx(do_more_labels=False)
    plot.set_logy(do_more_labels=False)

    plot.save(output_filename)


def do_genht_plot(dirname, output_filename, **plot_kwargs):
    qcd_file = cu.open_root_file(os.path.join(dirname, qgc.QCD_FILENAME))
    histname = "Dijet_gen/gen_ht"
    qcd_hist = cu.get_from_tfile(qcd_file, histname)
    conts = [
        Contribution(qcd_hist, label="QCD MC", line_color=ROOT.kRed)
    ]
    plot = Plot(conts,
                what='hist',
                ytitle="N",
                **plot_kwargs)
    plot.y_padding_max_log = 500
    plot.legend.SetY1(0.7)
    plot.plot("NOSTACK HIST E")
    plot.set_logx(do_more_labels=False)
    plot.set_logy(do_more_labels=False)

    plot.save(output_filename)


def do_genht_comparison_plot(dirname_label_pairs, output_filename, **plot_kwargs):
    qcd_files = [cu.open_root_file(os.path.join(dl[0], qgc.QCD_FILENAME)) for dl in dirname_label_pairs]
    histname = "Dijet_gen/gen_ht"
    qcd_hists = [cu.get_from_tfile(qf, histname) for qf in qcd_files]
    N = len(dirname_label_pairs)
    conts = [
        Contribution(qcd_hists[i], label=lab, marker_color=cu.get_colour_seq(i, N), 
                     line_color=cu.get_colour_seq(i, N), line_style=i+1, line_width=2,
                     subplot=qcd_hists[0] if i != 0 else None)
        for i, (d, lab) in enumerate(dirname_label_pairs)
    ]
    plot = Plot(conts,
                what='hist',
                ytitle="N",
                subplot_limits=(0.75, 1.25),
                subplot_type="ratio",
                subplot_title="* / %s" % (dirname_label_pairs[0][1]),
                **plot_kwargs)
    plot.y_padding_max_log = 500
    plot.legend.SetY1(0.7)
    plot.plot("NOSTACK HIST E")
    plot.set_logx(do_more_labels=False)
    plot.set_logy(do_more_labels=False)

    plot.save(output_filename)


def do_mc_pt_comparison_plot(dirname_label_pairs, output_filename, **plot_kwargs):
    # qcd_files = [cu.open_root_file(os.path.join(dl[0], qgc.QCD_FILENAME)) for dl in dirname_label_pairs]
    qcd_files = [cu.open_root_file(os.path.join(dl[0], qgc.QCD_PYTHIA_ONLY_FILENAME)) for dl in dirname_label_pairs]
    histname =  "Dijet_tighter/pt_jet1"
    qcd_hists = [cu.get_from_tfile(qf, histname) for qf in qcd_files]
    N = len(dirname_label_pairs)
    conts = [
        Contribution(qcd_hists[i], label=lab, marker_color=cu.get_colour_seq(i, N), 
                     line_color=cu.get_colour_seq(i, N), line_style=(i%3)+1 , line_width=2,
                     rebin_hist=1,
                     subplot=qcd_hists[0] if i != 0 else None)
        for i, (d, lab) in enumerate(dirname_label_pairs)
    ]
    plot = Plot(conts,
                what='hist',
                ytitle="N",
                subplot_limits=(0.5, 1.5),
                subplot_type="ratio",
                subplot_title="* / %s" % (dirname_label_pairs[0][1]),
                **plot_kwargs)
    plot.y_padding_max_log = 500
    plot.legend.SetY1(0.7)
    plot.plot("NOSTACK HIST E")
    plot.set_logx(do_more_labels=False)
    plot.set_logy(do_more_labels=False)

    plot.save(output_filename)


def do_pthat_comparison_plot(dirname_label_pairs, output_filename, **plot_kwargs):
    qcd_files = [cu.open_root_file(os.path.join(dl[0], qgc.QCD_PYTHIA_ONLY_FILENAME)) for dl in dirname_label_pairs]
    histname = "Dijet_gen/ptHat"
    qcd_hists = [cu.get_from_tfile(qf, histname) for qf in qcd_files]
    N = len(dirname_label_pairs)
    pthat_rebin = array('d', [15, 30, 50, 80, 120, 170, 300, 470, 600, 800, 1000, 1400, 1800, 2400, 3200, 5000])
    nbins = len(pthat_rebin)-1
    qcd_hists = [h.Rebin(nbins, cu.get_unique_str(), pthat_rebin) for h in qcd_hists]
    conts = [
        Contribution(qcd_hists[i], label=lab, marker_color=cu.get_colour_seq(i, N), 
                     line_color=cu.get_colour_seq(i, N), line_style=i+1, line_width=2,
                     subplot=qcd_hists[0] if i != 0 else None)
        for i, (d, lab) in enumerate(dirname_label_pairs)
    ]
    plot = Plot(conts,
                what='hist',
                ytitle="N",
                subplot_limits=(0.75, 1.25),
                subplot_type="ratio",
                subplot_title="* / %s" % (dirname_label_pairs[0][1]),
                **plot_kwargs)
    plot.y_padding_max_log = 500
    plot.legend.SetY1(0.7)
    plot.plot("NOSTACK HIST E")
    plot.set_logx(do_more_labels=False)
    plot.set_logy(do_more_labels=False)

    plot.save(output_filename)


def do_jetht_trigger_comparison_plot(dirname_label_pairs, output_dir, append="", title="", **plot_kwargs):
    # Unweighted pt, showing contributions from different triggers
    # Have to add in ZB manually
    zb_entry = {
        'label': 'HLT_ZeroBias',
        'color': ROOT.kMagenta-9,
        # 'scale': 35918219492.947 / 29048.362
        'scale': 1
    }

    jet_ht_entries = [
        {
            'ind': '0',
            'label': "PFJet40",
            'color': ROOT.kRed,
        },
        {
            'ind': '1',
            'label': "PFJet60",
            'color': ROOT.kBlue,
        },
        {
            'ind': '2',
            'label': "PFJet80",
            'color': ROOT.kGreen+2,
        },
        {
            'ind': '3',
            'label': "PFJet140",
            'color':  ROOT.kViolet+5,
        },
        {
            'ind': '4',
            'label': "PFJet200",
            'color': ROOT.kOrange,
        },
        {
            'ind': '5',
            'label': "PFJet260",
            'color': ROOT.kTeal,
        },
        {
            'ind': '6',
            'label': "PFJet320",
            'color': ROOT.kViolet,
        },
        {
            'ind': '7',
            'label': "PFJet400",
            'color': ROOT.kOrange-6
        },
        {
            'ind': '8',
            'label': "PFJet450",
            'color': ROOT.kAzure+1,
        },
    ]
    
    # PT JET 1
    zb_hist_names = ["Dijet_jet_hist_0/pt_1", "Dijet_jet_hist_unweighted_0/pt_1"]
    jet_ht_hist_names = ["Dijet_jet_hist_{ind}/pt_1", "Dijet_jet_hist_unweighted_{ind}/pt_1"]
    
    zb_root_files = [cu.open_root_file(os.path.join(dl[0], qgc.ZB_FILENAME)) for dl in dirname_label_pairs]
    jetht_root_files = [cu.open_root_file(os.path.join(dl[0], qgc.JETHT_FILENAME)) for dl in dirname_label_pairs]
    N = len(dirname_label_pairs)
    rebin = 2
    for zb_name, ht_name in  zip(zb_hist_names, jet_ht_hist_names):
        # add zeero bias ones
        this_data_entries = [
            Contribution(cu.get_from_tfile(zb_root_files[i], zb_name), 
                         label=zb_entry['label'] + ": " + l, 
                         marker_color=zb_entry['color'],
                         line_color=zb_entry['color'],
                         line_style=1+i,
                         rebin_hist=rebin,
                         )
            for i, (d, l) in enumerate(dirname_label_pairs)
        ]
        for c in this_data_entries[1:]:
            c.subplot = this_data_entries[0].obj

        # # add jet ht ones
        for ent in jet_ht_entries:
            histname = ht_name.format(ind=ent['ind'])
            new_entries = [
                Contribution(cu.get_from_tfile(jetht_root_files[i], histname), 
                             label=ent['label'] + ": " + l, 
                             marker_color=ent['color'], 
                             line_color=ent['color'], 
                             line_style=1+i, 
                             rebin_hist=rebin,
                             )
                for i, (d, l) in enumerate(dirname_label_pairs)
            ]
            for c in new_entries[1:]:
                c.subplot = new_entries[0].obj

            this_data_entries.extend(new_entries)

        plot = Plot(this_data_entries, what='hist',
                    title=title,
                    ytitle="N",
                    xtitle="p_{T}^{jet 1} [GeV]",
                    xlim=[30, 1000],
                    ylim=[1E3, None],
                    # ylim=[10, 1E8] if 'unweighted' in ht_name else [1, 1E12],
                    subplot_type='ratio',
                    subplot_title='* / %s' % dirname_label_pairs[0][1],
                    **plot_kwargs
                    )
        plot.default_canvas_size = (800, 600)
        plot.subplot_maximum = 5
        plot.y_padding_max_log = 500
        plot.legend.SetY1(0.7)
        plot.legend.SetY2(0.88)
        plot.legend.SetX1(0.5)
        plot.legend.SetNColumns(2)
        plot.plot("NOSTACK HISTE")
        plot.set_logx(do_more_labels=False)
        plot.set_logy(do_more_labels=False)
        output_filename = "%s/DataJetHTZB-pt_jet1%s%s.pdf" % (output_dir, "_unweighted" if 'unweighted' in zb_name else "", append)
        plot.save(output_filename)

    # ETA JET 1
    zb_hist_names = ["Dijet_jet_hist_unweighted_0/eta_1"]
    jet_ht_hist_names = ["Dijet_jet_hist_unweighted_{ind}/eta_1"]
    
    zb_root_files = [cu.open_root_file(os.path.join(dl[0], qgc.ZB_FILENAME)) for dl in dirname_label_pairs]
    jetht_root_files = [cu.open_root_file(os.path.join(dl[0], qgc.JETHT_FILENAME)) for dl in dirname_label_pairs]
    N = len(dirname_label_pairs)
    rebin = 2
    for zb_name, ht_name in  zip(zb_hist_names, jet_ht_hist_names):
        # add zero bias ones
        this_data_entries = [
            Contribution(cu.get_from_tfile(zb_root_files[i], zb_name), 
                         label=zb_entry['label'] + ": " + l, 
                         marker_color=zb_entry['color'],
                         line_color=zb_entry['color'],
                         line_style=1+i,
                         rebin_hist=rebin,
                         )
            for i, (d, l) in enumerate(dirname_label_pairs)
        ]
        for c in this_data_entries[1:]:
            c.subplot = this_data_entries[0].obj

        # plot zb
        plot = Plot(this_data_entries, what='hist',
                    title=title,
                    xtitle="y^{jet 1}",
                    ytitle="N",
                    subplot_type='ratio',
                    subplot_title='* / %s' % dirname_label_pairs[0][1],
                    # subplot_limits=(0, 5),
                    **plot_kwargs
                    )
        plot.subplot_maximum = 5
        plot.default_canvas_size = (800, 600)
        plot.y_padding_max_log = 500
        plot.legend.SetY1(0.7)
        plot.legend.SetY2(0.88)
        plot.legend.SetX1(0.5)
        plot.legend.SetNColumns(2)
        plot.plot("NOSTACK HISTE")
        output_filename = "%s/DataZB-eta_jet1%s%s.pdf" % (output_dir, "_unweighted" if 'unweighted' in zb_name else "", append)
        plot.save(output_filename)

        # add jet ht ones
        for ent in jet_ht_entries:
            histname = ht_name.format(ind=ent['ind'])
            this_data_entries = [
                Contribution(cu.get_from_tfile(jetht_root_files[i], histname), 
                             label=ent['label'] + ": " + l, 
                             marker_color=ent['color'], 
                             line_color=ent['color'], 
                             line_style=1+i, 
                             rebin_hist=rebin,
                             )
                for i, (d, l) in enumerate(dirname_label_pairs)
            ]
            for c in this_data_entries[1:]:
                c.subplot = this_data_entries[0].obj

            plot = Plot(this_data_entries, what='hist',
                        title=title,
                        xtitle="y^{jet 1}",
                        ytitle="N",
                        # xlim=[30, 1000],
                        # ylim=[10, 1E8] if 'unweighted' in ht_name else [1, 1E12],
                        subplot_type='ratio',
                        subplot_title='* / %s' % dirname_label_pairs[0][1],
                        **plot_kwargs
                        )
            plot.default_canvas_size = (800, 600)
            plot.y_padding_max_log = 500
            plot.legend.SetY1(0.7)
            plot.legend.SetY2(0.88)
            plot.legend.SetX1(0.5)
            plot.legend.SetNColumns(2)
            plot.plot("NOSTACK HISTE")
            output_filename = "%s/DataJetHTZB-%s_eta_jet1%s%s.pdf" % (output_dir, ent['label'], "_unweighted" if 'unweighted' in zb_name else "", append)
            plot.save(output_filename)


def do_zerobias_per_run_comparison_plot(dirname_label_pairs, output_dir, append="", title="", **plot_kwargs):
    runs = [
        (qgc.ZEROBIAS_RUNB_FILENAME, 'B'),
        (qgc.ZEROBIAS_RUNC_FILENAME, 'C'),
        (qgc.ZEROBIAS_RUND_FILENAME, 'D'),
        (qgc.ZEROBIAS_RUNE_FILENAME, 'E'),
        (qgc.ZEROBIAS_RUNF_FILENAME, 'F'),
        (qgc.ZEROBIAS_RUNG_FILENAME, 'G'),
        (qgc.ZEROBIAS_RUNH_FILENAME, 'H'),
    ]
    zb_entry = {
        'label': 'HLT_ZeroBias',
        'color': ROOT.kMagenta-9,
        # 'scale': 35918219492.947 / 29048.362
        'scale': 1
    }

    for filename, run_period in runs:
        zb_root_files = [cu.open_root_file(os.path.join(dl[0], filename)) for dl in dirname_label_pairs]

        # PT JET 1
        zb_hist_names = ["Dijet_jet_hist_0/pt_1", "Dijet_jet_hist_unweighted_0/pt_1"][1:]
        N = len(dirname_label_pairs)
        rebin = 2
        for zb_name in zb_hist_names:
            # add zeero bias ones
            this_data_entries = [
                Contribution(cu.get_from_tfile(zb_root_files[i], zb_name), 
                             label=zb_entry['label'] + " Run %s: " % run_period + l, 
                             marker_color=zb_entry['color'],
                             line_color=zb_entry['color'],
                             line_style=1+i,
                             rebin_hist=rebin,
                             )
                for i, (d, l) in enumerate(dirname_label_pairs)
            ]
            for c in this_data_entries[1:]:
                c.subplot = this_data_entries[0].obj

            plot = Plot(this_data_entries, what='hist',
                        title=title,
                        xtitle="p_{T}^{jet 1} [GeV]",
                        ytitle="N",
                        xlim=[30, 1000],
                        ylim=[1E3, None],
                        # ylim=[10, 1E8] if 'unweighted' in ht_name else [1, 1E12],
                        subplot_type='ratio',
                        subplot_title='* / %s' % dirname_label_pairs[0][1],
                        **plot_kwargs
                        )
            plot.default_canvas_size = (800, 600)
            plot.y_padding_max_log = 500
            plot.legend.SetY1(0.7)
            plot.legend.SetY2(0.88)
            plot.legend.SetX1(0.5)
            plot.legend.SetNColumns(2)
            plot.plot("NOSTACK HISTE")
            plot.set_logx(do_more_labels=False)
            plot.set_logy(do_more_labels=False)
            output_filename = "%s/DataJetHTZB-pt_jet1%s_Run%s%s.pdf" % (output_dir, "_unweighted" if 'unweighted' in zb_name else "", run_period, append)
            plot.save(output_filename)

        # ETA JET 1
        zb_hist_names = ["Dijet_jet_hist_unweighted_0/eta_1"]
        
        N = len(dirname_label_pairs)
        rebin = 2
        for zb_name in zb_hist_names:
            # add zero bias ones
            this_data_entries = [
                Contribution(cu.get_from_tfile(zb_root_files[i], zb_name), 
                             label=zb_entry['label'] + " Run %s: " % run_period + l, 
                             marker_color=zb_entry['color'],
                             line_color=zb_entry['color'],
                             line_style=1+i,
                             rebin_hist=rebin,
                             )
                for i, (d, l) in enumerate(dirname_label_pairs)
            ]
            for c in this_data_entries[1:]:
                c.subplot = this_data_entries[0].obj

            # plot zb
            plot = Plot(this_data_entries, what='hist',
                        title=title,
                        xtitle="y^{jet 1}",
                        ytitle="N",
                        subplot_type='ratio',
                        subplot_title='* / %s' % dirname_label_pairs[0][1],
                        # subplot_limits=(0, 5),
                        **plot_kwargs
                        )
            plot.subplot_maximum = 5
            plot.default_canvas_size = (800, 600)
            plot.y_padding_max_log = 500
            plot.legend.SetY1(0.7)
            plot.legend.SetY2(0.88)
            plot.legend.SetX1(0.5)
            plot.legend.SetNColumns(2)
            plot.plot("NOSTACK HISTE")
            output_filename = "%s/DataZB-eta_jet1%s_Run%s%s.pdf" % (output_dir, "_unweighted" if 'unweighted' in zb_name else "", run_period, append)
            plot.save(output_filename)



def do_jet_pt_vs_genht_plot(dirname, output_filename, title=""):
    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    qcd_file = cu.open_root_file(os.path.join(dirname, qgc.QCD_FILENAME))
    histname = "Dijet_tighter/pt_jet_vs_genHT"
    h = cu.get_from_tfile(qcd_file, histname)
    h.SetTitle(title)
    h.Draw("COLZ")
    h.GetXaxis().SetRangeUser(0, 200)
    h.GetYaxis().SetRangeUser(0, 200)
    canv.SaveAs(output_filename)


if __name__ == "__main__":
    odir = "ak4_ak8_data_mc_comparisons_new"
    histname = "Dijet_tighter/pt_jet1"
    pt_xlim = (20, 1000)
    do_data_mc_plot(ak4_v2_v2_dir,
                    histname,
                    output_filename=os.path.join(odir, "pythia_pt_jet1_ak4_v2_v2.pdf"),
                    title="AK4, v2 data, v2 MC, old",
                    xlim=pt_xlim)
    do_data_mc_plot(ak4_v3_v2_dir,
                    histname,
                    output_filename=os.path.join(odir, "pythia_pt_jet1_ak4_v3_v2.pdf"),
                    title="AK4, v3 data, v2 MC, new",
                    xlim=pt_xlim)
    # do_data_mc_plot(ak4_v3_v3_dir,
    #                 histname,
    #                 output_filename=os.path.join(odir, "pt_jet1_ak4_v3_v3.pdf"),
    #                 title="AK4, v3 data, v3 MC",
    #                 xlim=pt_xlim)
    # do_data_mc_plot(ak4_v3_v2_noJER_dir,
    #                 histname,
    #                 output_filename=os.path.join(odir, "pt_jet1_ak4_v3_v2_noJER.pdf"),
    #                 title="AK4, v3 data, v2 MC, no JER",
    #                 xlim=pt_xlim)
    # do_data_mc_plot(ak4_old_v2_noJER_dir,
    #                 histname,
    #                 output_filename=os.path.join(odir, "pt_jet1_ak4_v2_v2_old_noJER.pdf"),
    #                 title="AK4, v3 data, v2 MC, old, no JER",
    #                 xlim=pt_xlim)

    do_data_mc_plot(ak8_v2_v2_dir,
                    histname,
                    output_filename=os.path.join(odir, "pythia_pt_jet1_ak8_v2_v2.pdf"),
                    title="AK8, v2 data, v2 MC, old",
                    xlim=pt_xlim)
    do_data_mc_plot(ak8_v3_v2_dir,
                    histname,
                    output_filename=os.path.join(odir, "pythia_pt_jet1_ak8_v3_v2.pdf"),
                    title="AK8, v3 data, v2 MC, new",
                    xlim=pt_xlim)
    # do_data_mc_plot(ak8_v3_v3_dir,
    #                 histname,
    #                 output_filename=os.path.join(odir, "pt_jet1_ak8_v3_v3.pdf"),
    #                 title="AK8, v3 data, v3 MC",
    #                 xlim=pt_xlim)
    # do_data_mc_plot(ak8_v3_v2_noJER_dir,
    #                 histname,
    #                 output_filename=os.path.join(odir, "pt_jet1_ak8_v3_v2_noJER.pdf"),
    #                 title="AK8, v3 data, v2 MC, no JER",
    #                 xlim=pt_xlim)
    # do_data_mc_plot(ak8_old_v2_noJER_dir,
    #                 histname,
    #                 output_filename=os.path.join(odir, "pt_jet1_ak8_v2_v2_old_noJER.pdf"),
    #                 title="AK8, v2 data, v2 MC, old, no JER",
    #                 xlim=pt_xlim)

    # do_genht_comparison_plot([
    #                             (ak4_v2_v2_dir, "v2 MC old"),
    #                             (ak4_v3_v2_dir, "v2 MC new"),
    #                             (ak4_v3_v3_dir, "v3 MC")
    #                          ],
    #                          output_filename=os.path.join(odir, "genht_compare_ak4.pdf"),
    #                          title="AK4",
    #                          xlim=pt_xlim)
    # do_genht_comparison_plot([
    #                               (ak8_v2_v2_dir, "v2 MC old"),
    #                               (ak8_v3_v2_dir, "v2 MC new"),
    #                               (ak8_v3_v3_dir, "v3 MC")
    #                          ],
    #                          output_filename=os.path.join(odir, "genht_compare_ak8.pdf"),
    #                          title="AK8",
    #                          xlim=pt_xlim)


    do_mc_pt_comparison_plot([
                                (ak8_v2_v2_dir, "v2 MC old"),
                                (ak8_v3_v2_dir, "v2 MC new"),
                                # (ak8_v3_v3_dir, "v3 MC"),
                                # (ak8_old_v2_noJER_dir, "v2 MC old (noJER)"),
                                # (ak8_v3_v2_noJER_dir, "v2 MC new (noJER)"),
                          ],
                          output_filename=os.path.join(odir, "pythia_mc_pt_compare_ak8.pdf"),
                          title="AK8",
                          xlim=pt_xlim)

    do_mc_pt_comparison_plot([
                                (ak4_v2_v2_dir, "v2 MC old"),
                                (ak4_v3_v2_dir, "v2 MC new"),
                                # (ak4_v3_v3_dir, "v3 MC"),
                                # (ak4_old_v2_noJER_dir, "v2 MC old (noJER)"),
                                # (ak4_v3_v2_noJER_dir, "v2 MC new (noJER)"),
                          ],
                          output_filename=os.path.join(odir, "pythia_mc_pt_compare_ak4.pdf"),
                          title="AK4",
                          xlim=pt_xlim)

    pthat_xlim = (10, 5000)
    do_pthat_comparison_plot([
                                (ak4_v2_v2_dir, "v2 MC old"),
                                (ak4_v3_v2_dir, "v2 MC new"),
                             ],
                             output_filename=os.path.join(odir, "pthat_compare_ak4.pdf"),
                             title="AK4",
                             xlim=pthat_xlim)
    do_pthat_comparison_plot([
                                  (ak8_v2_v2_dir, "v2 MC old"),
                                  (ak8_v3_v2_dir, "v2 MC new"),
                             ],
                             output_filename=os.path.join(odir, "pthat_compare_ak8.pdf"),
                             title="AK8",
                             xlim=pthat_xlim)

    # do_genht_plot(ak4_v2_v2_dir,
    #               output_filename=os.path.join(odir, "genht_ak4_v2.pdf"),
    #               title="AK4, v2 MC",
    #               xlim=pt_xlim)
    # do_genht_plot(ak4_v3_v3_dir,
    #               output_filename=os.path.join(odir, "genht_ak4_v3.pdf"),
    #               title="AK4, v3 MC",
    #               xlim=pt_xlim)
    # do_genht_plot(ak4_v3_v2_dir,
    #               output_filename=os.path.join(odir, "genht_ak4_v2_new.pdf"),
    #               title="AK4, v2 MC",
    #               xlim=pt_xlim)
    
    # do_genht_plot(ak8_v2_v2_dir,
    #               output_filename=os.path.join(odir, "genht_ak8_v2.pdf"),
    #               title="AK8, v2 MC",
    #               xlim=pt_xlim)
    # do_genht_plot(ak8_v3_v3_dir,
    #               output_filename=os.path.join(odir, "genht_ak8_v3.pdf"),
    #               title="AK8, v3 MC",
    #               xlim=pt_xlim)
    # do_genht_plot(ak8_v3_v2_dir,
    #               output_filename=os.path.join(odir, "genht_ak8_v2_new.pdf"),
    #               title="AK8, v3 MC",
    #               xlim=pt_xlim)

    # # compare per-trigger
    # do_jetht_trigger_comparison_plot([
    #                                     (ak4_v2_v2_dir, "v2 data"),
    #                                     (ak4_v3_v3_dir, "v3 data"),
    #                                  ],
    #                                  output_dir=odir,
    #                                  append="_AK4",
    #                                  title="AK4",
    #                                  )
    # do_jetht_trigger_comparison_plot([
    #                                     (ak8_v2_v2_dir, "v2 data"),
    #                                     (ak8_v3_v3_dir, "v3 data"),
    #                                  ],
    #                                  output_dir=odir,
    #                                  append="_AK8",
    #                                  title="AK8",
    #                                  )

    # do_zerobias_per_run_comparison_plot([
    #                                     (ak4_v2_v2_dir, "v2 data"),
    #                                     (ak4_v3_v3_dir, "v3 data"),
    #                                  ],
    #                                  output_dir=odir,
    #                                  append="_AK4",
    #                                  title="AK4",
    #                                  )
    # do_zerobias_per_run_comparison_plot([
    #                                     (ak8_v2_v2_dir, "v2 data"),
    #                                     (ak8_v3_v3_dir, "v3 data"),
    #                                  ],
    #                                  output_dir=odir,
    #                                  append="_AK8",
    #                                  title="AK8",
    #                                  )

    # do_jet_pt_vs_genht_plot(ak4_v2_v2_dir,
    #                         output_filename=os.path.join(odir, "jet_pt_vs_genht_ak4_v2_old.pdf"),
    #                         title="AK4 V2 old MC")
    # do_jet_pt_vs_genht_plot(ak4_v3_v2_dir,
    #                         output_filename=os.path.join(odir, "jet_pt_vs_genht_ak4_v2_new.pdf"),
    #                         title="AK4 V2 MC")
    # do_jet_pt_vs_genht_plot(ak4_v3_v3_dir,
    #                         output_filename=os.path.join(odir, "jet_pt_vs_genht_ak4_v3.pdf"),
    #                         title="AK4 V3 MC")

    # do_jet_pt_vs_genht_plot(ak8_v2_v2_dir,
    #                         output_filename=os.path.join(odir, "jet_pt_vs_genht_ak8_v2_old.pdf"),
    #                         title="AK8 V2 old MC")
    # do_jet_pt_vs_genht_plot(ak8_v3_v2_dir,
    #                         output_filename=os.path.join(odir, "jet_pt_vs_genht_ak8_v2_new.pdf"),
    #                         title="AK8 V2 MC")
    # do_jet_pt_vs_genht_plot(ak8_v3_v3_dir,
    #                         output_filename=os.path.join(odir, "jet_pt_vs_genht_ak8_v3.pdf"),
    #                         title="AK8 V3 MC")
