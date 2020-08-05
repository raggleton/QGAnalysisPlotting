#!/usr/bin/env python

"""
Plots pt vs variable, and variable vs weights

To spot large-weight events, but also how the variables correlate with pT
"""

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


pt_str = "p_{T}^{RecoJet}"
pt_genjet_str = "p_{T}^{GenJet}"
genht_str = "H_{T}^{Gen}"
qscale_str = "q scale"
ptHat_str = "#hat{p}_{T}"
PU_ptHat_str = "max(#hat{p}_{T}^{PU})"

metric_str_dict = {
    "pt_jet_genHT_ratio" : "%s / %s" % (pt_str, genht_str),
    "pt_genjet_genHT_ratio" : "%s / %s" % (pt_genjet_str, genht_str),

    "pt_jet_qScale_ratio" : "%s / %s" % (pt_str, qscale_str),
    "pt_genjet_qScale_ratio" : "%s / %s" % (pt_genjet_str, qscale_str),

    "pt_jet_ptHat_ratio" : "%s / %s" % (pt_str, ptHat_str),
    "pt_genjet_ptHat_ratio" : "%s / %s" % (pt_genjet_str, ptHat_str),

    "PU_ptHat_genHT_ratio" : "%s / %s" % (PU_ptHat_str, genht_str),

    "PU_ptHat_ptHat_ratio" : "%s / %s" % (PU_ptHat_str, ptHat_str),
}


def get_var_str(histname):
    for k, v in metric_str_dict.items():
        if k in histname:
            return v
    return histname.replace("weight_vs_pt_vs_", "").replace("weight_vs_pt_genjet_vs_", "")


def do_var_vs_pt_plot(histname, input_filename, output_filename):
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    if h3d.GetEntries() == 0:
        return
    h2d = h3d.Project3D("zy")
    if "unweighted" in histname:
        h2d.SetTitle("Unweighted")
    else:
        h2d.SetTitle("Weighted")
    h2d.GetYaxis().SetTitle(get_var_str(histname))
    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetLeftMargin(0.1)
    canv.SetRightMargin(0.15)
    canv.SetLogz()
    h2d.Draw("COLZ")
    canv.SetLogx()
    h2d.GetXaxis().SetMoreLogLabels()
    canv.SaveAs(output_filename)

    canv.SetLogz(False)
    h2d_normed = cu.make_normalised_TH2(h2d, norm_axis='x')
    h2d_normed.Draw("COLZ")
    # h2d_normed.SetMinimum(1E-5)
    canv.SaveAs(output_filename.replace(".pdf", "_normX.pdf"))

    tf.Close()


def do_weight_vs_var_plot(histname, input_filename, output_filename):
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    if h3d.GetEntries() == 0:
        return
    h2d = h3d.Project3D("xz")
    if "unweighted" in histname:
        h2d.SetTitle("Unweighted")
    else:
        h2d.SetTitle("Weighted")
    h2d.GetXaxis().SetTitle(get_var_str(histname))

    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetLeftMargin(0.1)
    canv.SetRightMargin(0.15)
    canv.SetLogz()
    h2d.Draw("COLZ")
    canv.SaveAs(output_filename)
    tf.Close()


def do_weight_vs_var_plot_per_pt(histname, input_filename, output_filename):
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    if h3d.GetEntries() == 0:
        return

    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetLeftMargin(0.1)
    canv.SetRightMargin(0.15)
    canv.SetLogz()

    weight_str = "(unweighted)" if "unweighted" in histname else "(weighted)"
    for ibin in range(1, h3d.GetNbinsY()+1):
        h3d.GetYaxis().SetRange(ibin, ibin+1)
        h2d = h3d.Project3D("xz")
        if h2d.GetEntries() == 0:
            continue
        pt_low = h3d.GetYaxis().GetBinLowEdge(ibin)
        pt_high = h3d.GetYaxis().GetBinLowEdge(ibin+1)
        jet_str = pt_genjet_str if "_vs_pt_genjet_vs_" in histname else pt_str
        h2d.SetTitle("%g < %s < %g GeV %s" % (pt_low, jet_str, pt_high, weight_str))
        h2d.GetXaxis().SetTitle(get_var_str(histname))
        h2d.Draw("COLZ")
        jet_app = "_genjet" if "_vs_pt_genjet_vs" in histname else ""
        this_output_filename = output_filename.replace(".pdf", "_pt%s%gto%g.pdf" % (jet_app, pt_low, pt_high))
        canv.SaveAs(this_output_filename)
        canv.Clear()

    tf.Close()


def do_weight_vs_pt_plot(input_filename, output_filename):
    histname = "Weight_Presel/weight_vs_pt_vs_pt_jet_qScale_ratio"
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    if h3d.GetEntries() == 0:
        return
    h2d = h3d.Project3D("xy")

    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetLeftMargin(0.1)
    canv.SetRightMargin(0.15)
    canv.SetLogz()
    canv.SetLogx()
    h2d.Draw("COLZ")
    canv.SaveAs(output_filename)
    tf.Close()


def do_weight_vs_genjet_pt_plot(input_filename, output_filename):
    histname = "Weight_Presel/weight_vs_pt_vs_pt_genjet_qScale_ratio"
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    if h3d.GetEntries() == 0:
        return

    h2d = h3d.Project3D("xy")

    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetLeftMargin(0.1)
    canv.SetRightMargin(0.15)
    canv.SetLogz()
    canv.SetLogx()
    h2d.Draw("COLZ")
    canv.SaveAs(output_filename)
    tf.Close()


def do_jet_pt_with_var_cuts(histname, cuts, input_filename, output_filename):
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    if h3d.GetEntries() == 0:
        return
    pt_hists = []
    for cut in cuts:
        max_bin = h3d.GetZaxis().FindFixBin(cut)
        print("cut:", cut, "bin:", max_bin)
        h = h3d.ProjectionY("pt_var_lt_%g" % cut, 0, -1, 0, max_bin, "e")
        h2 = h.Clone()
        if h.GetEntries() > 0:
            h3 = qgp.hist_divide_bin_width(h2)
        pt_hists.append(h3)

    conts = [Contribution(h, label=" < %g" % cut,
                          line_color=cu.get_colour_seq(ind, len(cuts)),
                          marker_color=cu.get_colour_seq(ind, len(cuts)),
                          subplot=pt_hists[-1])
             for ind, (h, cut) in enumerate(zip(pt_hists, cuts))]

    jet_str = pt_genjet_str if "_vs_pt_genjet_vs_" in histname else pt_str
    weight_str = "(unweighted)" if "unweighted" in histname else "(weighted)"
    ratio_lims = None
    plot = Plot(conts, what='hist',
                title='%s for cuts on %s %s' % (jet_str, get_var_str(histname), weight_str),
                xtitle=None,
                ytitle='N',
                # xlim=None, ylim=None,
                legend=True,
                subplot_type='ratio',
                subplot_title='* / var < %g' % cuts[-1],
                subplot_limits=ratio_lims,
                has_data=False)
    plot.y_padding_max_log = 200
    plot.subplot_maximum_ceil = 5
    plot.subplot_maximum_floor = 1.02
    plot.subplot_minimum_ceil = 0.98
    plot.legend.SetY1(0.7)
    plot.legend.SetY2(0.89)
    plot.legend.SetX1(0.78)
    plot.legend.SetX2(0.88)
    plot.plot("NOSTACK HISTE")
    plot.set_logx(True, do_more_labels=True)
    plot.set_logy(True, do_more_labels=False)
    plot.save(output_filename)


def do_jet_pt_rel_error_with_var_cuts(histname, cuts, input_filename, output_filename):
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    if h3d.GetEntries() == 0:
        return
    pt_hists = []
    for cut in cuts:
        max_bin = h3d.GetZaxis().FindFixBin(cut)
        print("cut:", cut, "bin:", max_bin)
        h = h3d.ProjectionY("pt_var_lt_%g" % cut, 0, -1, 0, max_bin, "e")
        h2 = h.Clone()
        if h.GetEntries() > 0:
            h3 = qgp.hist_divide_bin_width(h2)
        # convert bin contents to bin error/bin contents
        for ibin in range(1, h2.GetNbinsX()+1):
            if h3.GetBinContent(ibin) == 0:
                continue
            h3.SetBinContent(ibin, h3.GetBinError(ibin) / h3.GetBinContent(ibin))
            h3.SetBinError(ibin, 0)
        pt_hists.append(h3)

    conts = [Contribution(h, label=" < %g" % cut,
                          line_color=cu.get_colour_seq(ind, len(cuts)),
                          marker_color=cu.get_colour_seq(ind, len(cuts)),
                          subplot=pt_hists[-1])
             for ind, (h, cut) in enumerate(zip(pt_hists, cuts))]

    jet_str = pt_genjet_str if "_vs_pt_genjet_vs_" in histname else pt_str
    weight_str = "(unweighted)" if "unweighted" in histname else "(weighted)"
    ratio_lims = None
    plot = Plot(conts, what='hist',
                title='%s for cuts on %s %s' % (jet_str, get_var_str(histname), weight_str),
                xtitle=None,
                ytitle='Relative error',
                # xlim=None, ylim=None,
                legend=True,
                subplot_type='ratio',
                subplot_title='* / var < %g' % cuts[-1],
                subplot_limits=ratio_lims,
                has_data=False)
    plot.y_padding_max_log = 200
    plot.subplot_maximum_ceil = 5
    plot.subplot_maximum_floor = 1.02
    plot.subplot_minimum_ceil = 0.98
    plot.legend.SetY1(0.7)
    plot.legend.SetY2(0.89)
    plot.legend.SetX1(0.78)
    plot.legend.SetX2(0.88)
    plot.plot("NOSTACK HISTE")
    plot.set_logx(True, do_more_labels=True)
    plot.set_logy(True, do_more_labels=False)
    plot.save(output_filename)


def do_cut_scan_per_pt(histname, input_filename, output_filename):
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    h3d_unweighted = cu.get_from_tfile(tf, histname+"_unweighted")
    if h3d.GetEntries() == 0:
        return

    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetLeftMargin(0.12)
    canv.SetRightMargin(0.12)
    # canv.SetLogz()
    h2d = h3d.Project3D("zy") # var vs pt
    h2d_unweighted = h3d_unweighted.Project3D("zy") # var vs pt

    var_name = os.path.basename(histname).replace("weight_vs_pt_vs_", "").replace("weight_vs_pt_genjet_vs_", "")

    for ibin in range(3, h2d.GetNbinsX()+1):  # iterate over pt bins
        if h2d.Integral(ibin, ibin+1, 0, -1) == 0:
            # data.append(None)
            continue

        pt_low = h3d.GetYaxis().GetBinLowEdge(ibin)
        pt_high = h3d.GetYaxis().GetBinLowEdge(ibin+1)

        data = []
        data_unweighted = []
        # Do integral, error for increasingly looser cuts
        # find maximum in this pt bin
        # yes I probably should collapse to a 1D hist and use GetMaximumBin
        max_val, max_bin = 0, 0
        for icut in range(1, h2d.GetNbinsY()+1):
            val =  h2d.GetBinContent(ibin, icut)
            if val > max_val:
                max_val = val
                max_bin = icut

        for icut in range(max_bin + 5, h2d.GetNbinsY()+2):
            err = array('d', [0])
            count = h2d.IntegralAndError(ibin, ibin+1, 1, icut-1, err)
            if count == 0:
                continue
            data.append([count, err[0], h2d.GetYaxis().GetBinLowEdge(icut)])

            err = array('d', [0])
            count = h2d_unweighted.IntegralAndError(ibin, ibin+1, 1, icut, err)
            data_unweighted.append([count, err[0], h2d.GetYaxis().GetBinLowEdge(icut)])

        # Plot count, rel error vs cut value
        cuts = [d[2] for d in data]
        gr_count = ROOT.TGraph(len(data), array('d', cuts), array('d', [d[0] / data[-1][0] for d in data]))

        gr_count.SetMarkerColor(ROOT.kRed)
        gr_count.SetMarkerStyle(22)
        gr_count.SetLineColor(ROOT.kRed)
        gr_count.SetTitle("%g < p_{T} < %g GeV;%s cut (<);Count (relative to loosest cut)" % (pt_low, pt_high, get_var_str(histname)))

        gr_count_unweighted = ROOT.TGraph(len(data), array('d', cuts), array('d', [d[0] / data_unweighted[-1][0] for d in data_unweighted]))
        gr_count_unweighted.SetMarkerColor(ROOT.kBlack)
        gr_count_unweighted.SetMarkerStyle(23)
        gr_count_unweighted.SetLineColor(ROOT.kBlack)
        gr_count_unweighted.SetTitle("%g < p_{T} < %g GeV;%s cut (<);Count (relative to loosest cut)" % (pt_low, pt_high, get_var_str(histname)))

        leg = ROOT.TLegend(0.7, 0.5, 0.85, 0.65)
        leg.AddEntry(gr_count, "Weighted", "LP")
        leg.AddEntry(gr_count_unweighted, "Unweighted", "LP")

        # gr_rel_err = ROOT.TGraph(len(data), array('d', cuts), array('d', [(d[1] / d[0]) if d[0] != 0 else 0 for d in data ]))
        # gr_rel_err.SetMarkerColor(ROOT.kRed)
        # gr_rel_err.SetMarkerStyle(22)
        # gr_rel_err.SetLineColor(ROOT.kRed)
        # gr_rel_err.SetTitle("%g < p_{T} < %g GeV;%s cut (<);Rel. error" % (pt_low, pt_high, var_name))

        canv.SetLogy(False)

        gr_count.Draw("ALP")
        gr_count_unweighted.Draw("LP")
        gr_count.Draw("LP")
        leg.Draw()
        canv.SaveAs(output_filename.replace(".pdf", "_count_pt%gto%g.pdf" % (pt_low, pt_high)))

        # canv.Clear()
        # gr_rel_err.Draw("ALP")
        # canv.SetLogy()
        # gr_rel_err.GetYaxis().SetMoreLogLabels()
        # canv.SaveAs(output_filename.replace(".pdf", "_rel_err_pt%gto%g.pdf" % (pt_low, pt_high)))


    tf.Close()


def do_cut_roc_per_pt(histname, input_filename, output_filename):
    """Plot fractional # unweighted vs fraction # weighted, for different cuts
    Not a true ROC, but kinda like one
    """
    tf = cu.open_root_file(input_filename)
    h3d = cu.get_from_tfile(tf, histname)
    h3d_unweighted = cu.get_from_tfile(tf, histname+"_unweighted")
    if h3d.GetEntries() == 0:
        return

    canv = ROOT.TCanvas(cu.get_unique_str(), "", 800, 600)
    canv.SetTicks(1, 1)
    canv.SetLeftMargin(0.12)
    canv.SetRightMargin(0.12)
    # canv.SetLogz()
    h2d = h3d.Project3D("zy") # var vs pt
    h2d_unweighted = h3d_unweighted.Project3D("zy") # var vs pt

    var_name = os.path.basename(histname).replace("weight_vs_pt_vs_", "").replace("weight_vs_pt_genjet_vs_", "")

    for ibin in range(3, h2d.GetNbinsX()+1):  # iterate over pt bins
        if h2d.Integral(ibin, ibin+1, 0, -1) == 0:
            # data.append(None)
            continue

        pt_low = h3d.GetYaxis().GetBinLowEdge(ibin)
        pt_high = h3d.GetYaxis().GetBinLowEdge(ibin+1)

        data = []
        data_unweighted = []
        # Do integral, error for increasingly looser cuts
        # find maximum in this pt bin
        # yes I probably should collapse to a 1D hist and use GetMaximumBin
        max_val, max_bin = 0, 0
        for icut in range(1, h2d.GetNbinsY()+1):
            val =  h2d.GetBinContent(ibin, icut)
            if val > max_val:
                max_val = val
                max_bin = icut

        # for icut in range(max_bin + 5, h2d.GetNbinsY()+2):
        for icut in range(2, h2d.GetNbinsY()+2):
            err = array('d', [0])
            count = h2d.IntegralAndError(ibin, ibin+1, 1, icut-1, err)
            if count == 0:
                continue
            data.append([count, err[0], h2d.GetYaxis().GetBinLowEdge(icut)])

            err = array('d', [0])
            count = h2d_unweighted.IntegralAndError(ibin, ibin+1, 1, icut-1, err)
            data_unweighted.append([count, err[0], h2d.GetYaxis().GetBinLowEdge(icut)])

        cuts = [d[2] for d in data]

        weighted_fractions = np.array([d[0] / data[-1][0] for d in data])
        unweighted_fractions = np.array([d[0] / data_unweighted[-1][0] for d in data_unweighted])
        gr_count = ROOT.TGraph(len(data), unweighted_fractions, weighted_fractions)
        gr_count.SetMarkerColor(ROOT.kRed)
        gr_count.SetMarkerSize(0)
        gr_count.SetMarkerStyle(21)
        gr_count.SetLineColor(ROOT.kRed)
        gr_count.SetTitle("%s, %g < p_{T} < %g GeV;Relative unweighted count;Relative weighted count" % (get_var_str(histname), pt_low, pt_high))

        canv.SetLogx(False)
        canv.SetLogy(False)

        ROOT.TGaxis.SetMaxDigits(2)

        gr_count.Draw("ALP")

        ROOT.TGaxis.SetMaxDigits(2)
        unweighted_min = 0.9999

        # Calculate differences between points
        unweighted_diffs = unweighted_fractions[1:] - unweighted_fractions[:-1]
        weighted_diffs = weighted_fractions[1:] - weighted_fractions[:-1]
        big_diff_inds = []
        for ind, (u, w) in enumerate(zip(unweighted_diffs, weighted_diffs)):
            # look for big diff in weighted frac, small diff in unweighted,
            # with a limit on the minimum size of unweighted frac
            # (only trying to remove a few events)
            if u > 0 and w / u > 100 and u < 0.005 and unweighted_fractions[ind] > unweighted_min:
                big_diff_inds.append(ind)

        # if "pt_jet_genHT_ratio" in histname and pt_low == 186:
        #     for u, w in zip(unweighted_diffs, weighted_diffs):
        #         print(u, w)
        #     print(big_diff_inds)

        # make graph of big diff points, add annotations of cuts
        if len(big_diff_inds) > 0:
            gr_big_diffs = ROOT.TGraph(len(big_diff_inds), array('d', [unweighted_fractions[i+1] for i in big_diff_inds]), array('d', [weighted_fractions[i+1] for i in big_diff_inds]))
            gr_big_diffs.SetLineWidth(0)
            gr_big_diffs.SetMarkerColor(ROOT.kBlue)
            gr_big_diffs.SetMarkerStyle(25)
            latexs = []
            for i, ind in enumerate(big_diff_inds[:]):
                latex = ROOT.TLatex(gr_big_diffs.GetX()[i], gr_big_diffs.GetY()[i], " < %.2f" % cuts[ind+1])
                latex.SetTextSize(0.02)
                latex.SetTextColor(ROOT.kBlue)
                gr_big_diffs.GetListOfFunctions().Add(latex)
                latexs.append(latex)
            gr_big_diffs.Draw("*")

        gr_count.GetXaxis().SetLimits(unweighted_min, 1)

        # find corresponding value for weighted to set axis range
        weighted_min = 0
        for ind, u in enumerate(unweighted_fractions):
            if u >= unweighted_min:
                weighted_min = weighted_fractions[ind-1]
                if ind == len(unweighted_fractions) - 1:
                    weighted_min = 0
                break
        gr_count.GetHistogram().SetMinimum(weighted_min*1.1 - 0.1)
        gr_count.GetHistogram().SetMaximum(1)
        canv.SaveAs(output_filename.replace(".pdf", "_count_pt%gto%g.pdf" % (pt_low, pt_high)))

        # do a version zoomed out
        canv.Clear()
        gr_count.Draw("ALP")
        unweighted_min = 0.95
        gr_count.GetXaxis().SetLimits(unweighted_min, 1)
        weighted_min = 0
        for ind, u in enumerate(unweighted_fractions):
            if u >= unweighted_min:
                weighted_min = weighted_fractions[ind-1]
                if ind == len(unweighted_fractions) - 1:
                    weighted_min = 0
                break
        gr_count.GetHistogram().SetMinimum(weighted_min*1.1 - 0.1)
        gr_count.GetHistogram().SetMaximum(1)
        canv.SaveAs(output_filename.replace(".pdf", "_count_pt%gto%g_all.pdf" % (pt_low, pt_high)))


    tf.Close()


if __name__ == "__main__":
    MAINDIR="/Volumes/Extreme SSD/Projects/QGAnalysis"
    histnames = [
        "weight_vs_pt_vs_pt_jet_genHT_ratio",
        "weight_vs_pt_vs_pt_genjet_genHT_ratio",
        "weight_vs_pt_genjet_vs_pt_genjet_genHT_ratio",

        "weight_vs_pt_vs_pt_jet_qScale_ratio",
        "weight_vs_pt_vs_pt_genjet_qScale_ratio",
        "weight_vs_pt_genjet_vs_pt_genjet_qScale_ratio",

        "weight_vs_pt_vs_pt_jet_ptHat_ratio",
        "weight_vs_pt_vs_pt_genjet_ptHat_ratio",
        "weight_vs_pt_genjet_vs_pt_genjet_ptHat_ratio",

        "weight_vs_pt_vs_pt_jet_genHT_ratio_unweighted",
        "weight_vs_pt_vs_pt_genjet_genHT_ratio_unweighted",
        "weight_vs_pt_genjet_vs_pt_genjet_genHT_ratio_unweighted",

        "weight_vs_pt_vs_pt_jet_qScale_ratio_unweighted",
        "weight_vs_pt_vs_pt_genjet_qScale_ratio_unweighted",
        "weight_vs_pt_genjet_vs_pt_genjet_qScale_ratio_unweighted",

        "weight_vs_pt_vs_pt_jet_ptHat_ratio_unweighted",
        "weight_vs_pt_vs_pt_genjet_ptHat_ratio_unweighted",
        "weight_vs_pt_genjet_vs_pt_genjet_ptHat_ratio_unweighted",

        "weight_vs_pt_vs_PU_ptHat_genHT_ratio",
        "weight_vs_pt_vs_PU_ptHat_genHT_ratio_unweighted",

        "weight_vs_pt_vs_PU_ptHat_ptHat_ratio",
        "weight_vs_pt_vs_PU_ptHat_ptHat_ratio_unweighted",
    ]

    hist_dirs = ['Weight_Presel', 'Weight_Reco_sel']

    workdirs = [
        "workdir_102X_v3data_v2mc_ak4puppi_data_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_weightHists",
        "workdir_102X_v3data_v2mc_ak8puppi_data_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_weightHists"
    ]


    qcd_filenames = [
        qgc.QCD_FILENAME,
        qgc.QCD_PYTHIA_ONLY_FILENAME,
        qgc.QCD_HERWIG_FILENAME,
    ]

    # ROOT.gStyle.SetPalette(ROOT.kPastel)

    for workdir in workdirs:
        for qcd_filename in qcd_filenames:

            ifile = os.path.join(MAINDIR, workdir, qcd_filename)

            for hist_dir in hist_dirs:

                dir_append = hist_dir.lower().replace("weight_", "")
                odir = os.path.join(MAINDIR, workdir, "var_jet_pt_weights_%s_%s" % (dir_append, qcd_filename.replace("uhh2.AnalysisModuleRunner.MC.MC_", "").replace(".root", "")))
                cu.check_dir_exists_create(odir)

                # do_weight_vs_pt_plot(input_filename=ifile,
                #                      output_filename=os.path.join(odir, "weight_vs_pt.pdf"))

                # do_weight_vs_genjet_pt_plot(input_filename=ifile,
                #                             output_filename=os.path.join(odir, "weight_vs_genjet_pt.pdf"))

                for ind, base_hname in enumerate(histnames):
                    hname = os.path.join(hist_dir, base_hname)
                    # do_var_vs_pt_plot(histname=hname,
                    #                   input_filename=ifile,
                    #                   output_filename=os.path.join(odir, base_hname.replace("weight_vs_", "") + ".pdf"))

                    # if "_vs_pt_" in hname:
                    #     # skip the _vs_pt_genjet_ since identical, and we collapse that axis anyway
                    #     do_weight_vs_var_plot(histname=hname,
                    #                           input_filename=ifile,
                    #                           output_filename=os.path.join(odir, base_hname.replace("_vs_pt_vs_", "_").replace("_vs_pt_genjet_vs_", "") + ".pdf"))

                    cuts = [10, 9, 8, 7, 6, 5, 4, 3, 2]
                    cuts = [10, 8, 6, 4, 3, 2]
                    if "PU_ptHat" in hname:
                        cuts.extend([1.5, 1, 0.5])
                    if 'pt_genjet_genHT_ratio' in hname:
                        cuts.extend([1.5, 1])
                    if 'pt_jet_genHT_ratio' in hname:
                        cuts.extend([1.5, 1])
                    # do_jet_pt_with_var_cuts(histname=hname,
                    #                         cuts=cuts,
                    #                         input_filename=ifile,
                    #                         output_filename=os.path.join(odir, base_hname.replace("weight_vs_pt_vs_", "pt_cuts_").replace("weight_vs_pt_genjet_vs_", "pt_genjet_cuts_") + ".pdf"))

                    # do_jet_pt_rel_error_with_var_cuts(histname=hname,
                    #                                   cuts=cuts,
                    #                                   input_filename=ifile,
                    #                                   output_filename=os.path.join(odir, base_hname.replace("weight_vs_pt_vs_", "pt_cuts_rel_error_").replace("weight_vs_pt_genjet_vs_", "pt_genjet_rel_error_") + ".pdf"))

                    # do_weight_vs_var_plot_per_pt(histname=hname,
                    #                              input_filename=ifile,
                    #                              output_filename=os.path.join(odir, base_hname.replace("_vs_pt_vs_", "_").replace("_vs_pt_genjet_vs_", "_") + ".pdf"))

                    if "unweighted" not in hname and "_vs_pt_genjet_vs_" not in hname:
                    #     do_cut_scan_per_pt(histname=hname,
                    #                        input_filename=ifile,
                    #                        output_filename=os.path.join(odir, base_hname.replace("weight_vs_pt_vs_", "cuts_scan_") + ".pdf"))
                        do_cut_roc_per_pt(histname=hname,
                                          input_filename=ifile,
                                          output_filename=os.path.join(odir, base_hname.replace("weight_vs_pt_vs_", "cuts_roc_") + ".pdf"))
