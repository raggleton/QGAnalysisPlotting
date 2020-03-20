#!/usr/bin/env python


"""Plots for stitching Herwig DY samples to check them."""


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

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetTitleXOffset(1.5)

# Control plot output format
OUTPUT_FMT = "pdf"


def style_incl_hist(hist, rebin=2):
    col = ROOT.kBlack
    hist.SetLineColor(col)
    hist.SetMarkerColor(col)
    hist.Rebin(rebin)


def style_low_pt_hist(hist, rebin=2):
    col = ROOT.kBlue
    hist.SetLineColor(col)
    hist.SetMarkerColor(col)
    hist.Rebin(rebin)


def style_high_pt_hist(hist, rebin=2):
    col = ROOT.kRed
    hist.SetLineColor(col)
    hist.SetMarkerColor(col)
    hist.Rebin(rebin)


def style_merged_hist(hist, rebin=2):
    col = ROOT.kGreen+2
    hist.SetLineColor(col)
    hist.SetMarkerColor(col)
    hist.Rebin(rebin)


def create_legend(entries):
    leg = ROOT.TLegend(0.7, 0.7, 0.88, 0.88)
    for ent in entries:
        leg.AddEntry(ent[0], ent[1])
    return leg


def do_plot(incl_hist, low_pt_hist, high_pt_hist, merged_hist, output_filename, logx=False, logy=False, xmin=None, xmax=None):
    c = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    c.SetTicks(1, 1)
    if logx:
        c.SetLogx()
    if logy:
        c.SetLogy()

    xax = incl_hist.GetXaxis()
    xmin = xmin or xax.GetXmin()
    xmax = xmax or xax.GetXmax()

    first_drawing = True
    leg_entries = []
    if incl_hist is not None:
        incl_hist.SetAxisRange(xmin, xmax, 'X')
        draw_opt = ""
        incl_hist.Draw(draw_opt)
        first_drawing = False
        leg_entries.append([incl_hist, 'Inclusive'])

    if low_pt_hist is not None and high_pt_hist is not None:
        low_pt_hist.SetAxisRange(xmin, xmax, 'X')
        high_pt_hist.SetAxisRange(xmin, xmax, 'X')
        hst = ROOT.THStack(cu.get_unique_str(), "")
        hst.Add(low_pt_hist)
        hst.Add(high_pt_hist)
        draw_opt = ""
        if not first_drawing:
            draw_opt += " SAME"
        hst.Draw(draw_opt)
        first_drawing = False
        leg_entries.append([low_pt_hist, 'Inclusive, qScale < 200 GeV'])
        leg_entries.append([high_pt_hist, 'My sample, qScale > 200 GeV'])

    if merged_hist is not None:
        merged_hist.SetAxisRange(xmin, xmax, 'X')
        draw_opt = ""
        if not first_drawing:
            draw_opt += " SAME"
        merged_hist.Draw(draw_opt)
        first_drawing = False
        leg_entries.append([merged_hist, 'Combined inclusive + high pT'])

    leg = create_legend(leg_entries)
    leg.Draw()

    c.SaveAs(output_filename)
    return c


if __name__ == "__main__":
    workdir = "workdir_ak4puppi_data_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noPtHatCut_noPtReweight_noZjet2Cut_zPt30_trkSF_wtaAK_fixPrescales_sameGenCuts"
    tfile_incl = cu.TFileCacher(os.path.join(workdir, qgc.DY_HERWIG_INCL_FILENAME))
    tfile_low_pt = cu.TFileCacher(os.path.join(workdir, qgc.DY_HERWIG_LOW_PT_FILENAME))
    tfile_high_pt = cu.TFileCacher(os.path.join(workdir, qgc.DY_HERWIG_HIGH_PT_FILENAME))
    tfile_merged = cu.TFileCacher(os.path.join(workdir, qgc.DY_HERWIG_LOW_HIGH_PT_FILENAME))

    # q scale plot
    qscale_name = "ZPlusJets_gen/q_scale"
    h_qscale_incl = tfile_incl.Get(qscale_name)
    style_incl_hist(h_qscale_incl)
    h_qscale_low_pt = tfile_low_pt.Get(qscale_name)
    style_low_pt_hist(h_qscale_low_pt)
    h_qscale_high_pt = tfile_high_pt.Get(qscale_name)
    style_high_pt_hist(h_qscale_high_pt)
    h_qscale_merged = tfile_merged.Get(qscale_name)
    style_merged_hist(h_qscale_merged)
    output_filename = os.path.join(workdir, "qscale.pdf")
    do_plot(h_qscale_incl, h_qscale_low_pt, h_qscale_high_pt, None, output_filename, logy=True)

    # genjet pt plot
    pt_genjet_name = "ZPlusJets_gen/pt_jet1"
    pt_genjet_incl = tfile_incl.Get(pt_genjet_name)
    style_incl_hist(pt_genjet_incl)
    pt_genjet_low_pt = tfile_low_pt.Get(pt_genjet_name)
    style_low_pt_hist(pt_genjet_low_pt)
    pt_genjet_high_pt = tfile_high_pt.Get(pt_genjet_name)
    style_high_pt_hist(pt_genjet_high_pt)
    pt_genjet_merged = tfile_merged.Get(pt_genjet_name)
    style_merged_hist(pt_genjet_merged)
    output_filename = os.path.join(workdir, "pt_genjet.pdf")
    do_plot(pt_genjet_incl, pt_genjet_low_pt, pt_genjet_high_pt, None, output_filename, logy=True, xmax=1000)

    # recojet pt plot
    pt_jet_name = "ZPlusJets/pt_jet1"
    rebin = 4
    h_pt_jet_incl = tfile_incl.Get(pt_jet_name)
    style_incl_hist(h_pt_jet_incl, rebin)
    h_pt_jet_low_pt = tfile_low_pt.Get(pt_jet_name)
    style_low_pt_hist(h_pt_jet_low_pt, rebin)
    h_pt_jet_high_pt = tfile_high_pt.Get(pt_jet_name)
    style_high_pt_hist(h_pt_jet_high_pt, rebin)
    h_pt_jet_merged = tfile_merged.Get(pt_jet_name)
    style_merged_hist(h_pt_jet_merged, rebin)
    output_filename = os.path.join(workdir, "pt_jet.pdf")
    do_plot(h_pt_jet_incl, h_pt_jet_low_pt, h_pt_jet_high_pt, None, output_filename, logy=True, logx=False, xmin=20, xmax=1000)

    # gen jet kT
    kt_name = "ZPlusJets/genjet_kt"
    rebin = 4
    h_kt_incl = tfile_incl.Get(kt_name)
    style_incl_hist(h_kt_incl, rebin)
    h_kt_low_pt = tfile_low_pt.Get(kt_name)
    style_low_pt_hist(h_kt_low_pt, rebin)
    h_kt_high_pt = tfile_high_pt.Get(kt_name)
    style_high_pt_hist(h_kt_high_pt, rebin)
    h_kt_merged = tfile_merged.Get(kt_name)
    style_merged_hist(h_kt_merged, rebin)
    output_filename = os.path.join(workdir, "genjet_kt.pdf")
    do_plot(h_kt_incl, h_kt_low_pt, h_kt_high_pt, None, output_filename, logy=True, logx=False, xmin=20, xmax=1000)

    # ptZ plot
    ptZ_name = "ZPlusJets/pt_mumu"
    rebin = 4
    h_ptZ_incl = tfile_incl.Get(ptZ_name)
    style_incl_hist(h_ptZ_incl, rebin)
    h_ptZ_low_pt = tfile_low_pt.Get(ptZ_name)
    style_low_pt_hist(h_ptZ_low_pt, rebin)
    h_ptZ_high_pt = tfile_high_pt.Get(ptZ_name)
    style_high_pt_hist(h_ptZ_high_pt, rebin)
    h_ptZ_merged = tfile_merged.Get(ptZ_name)
    style_merged_hist(h_ptZ_merged, rebin)
    output_filename = os.path.join(workdir, "ptZ.pdf")
    do_plot(h_ptZ_incl, h_ptZ_low_pt, h_ptZ_high_pt, None, output_filename, logy=True, logx=False, xmin=30, xmax=1000)
