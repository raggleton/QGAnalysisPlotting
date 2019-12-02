#!/usr/bin/env python


"""
Make Z+jet data-background plots
"""


from __future__ import print_function

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
os.nice(10)

# My stuff
# from comparator import Contribution, Plot, grab_obj
# import qg_common as qgc
# import qg_general_plots as qgg
import common_utils as cu


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
# ROOT.gErrorIgnoreLevel = ROOT.kWarning


tfile_data = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu.root")
tfile_dy = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root")
# tfile_dy = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root")
tfile_ww = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet/uhh2.AnalysisModuleRunner.MC.MC_WW.root")
tfile_wz = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet/uhh2.AnalysisModuleRunner.MC.MC_WZ.root")
tfile_zz = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet/uhh2.AnalysisModuleRunner.MC.MC_ZZ.root")
tfile_tt = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet/uhh2.AnalysisModuleRunner.MC.MC_TTBAR.root")
# tfile_st_t = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet/uhh2.AnalysisModuleRunner.MC.MC_ST_T.root")
# tfile_st_tw = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet/uhh2.AnalysisModuleRunner.MC.MC_ST_TW.root")
# tfile_st_s = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet/uhh2.AnalysisModuleRunner.MC.MC_ST_S.root")


tfile_data = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_noZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu.root")
tfile_dy = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_noZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root")
tfile_ww = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_noZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut/uhh2.AnalysisModuleRunner.MC.MC_WW.root")
tfile_wz = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_noZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut/uhh2.AnalysisModuleRunner.MC.MC_WZ.root")
tfile_zz = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_noZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut/uhh2.AnalysisModuleRunner.MC.MC_ZZ.root")
tfile_tt = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_noZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut/uhh2.AnalysisModuleRunner.MC.MC_TTBAR.root")


class Entry(object):
    def __init__(self, hist, label, is_data=False, is_bkg=False):
        self.hist = hist
        self.label = label
        self.is_data = is_data
        self.is_bkg = is_bkg

    @property
    def integral(self):
        return self.hist.Integral()


def get_x_bin_low_edge(entries, index):
    return entries[0].hist.GetXaxis().GetBinLowEdge(index)


def print_data_bg_stats(data_entries, mc_entries):
    """Print various stats about data & MC/BG composition"""
    data_total = sum(h.integral for h in data_entries)
    mc_total = sum(h.integral for h in mc_entries)
    bg_total = sum(h.integral for h in mc_entries if h.is_bkg)
    print('--'*80)
    print('DATA TOTAL:', data_total)
    print('MC TOTAL:', mc_total)
    print('BKG TOTAL:', bg_total)
    print('MC / DATA:', mc_total / data_total)
    print('BKG / DATA:', bg_total / data_total)
    print('BKG / MC:', bg_total / mc_total)
    print('--'*80)
    data_bins = [sum(ent.hist.Integral(i, i+1) for ent in data_entries) for i in range(1, data_entries[0].hist.GetNbinsX()+1)]
    bg_bins = [sum(ent.hist.Integral(i, i+1) for ent in mc_entries if ent.is_bkg) for i in range(1, mc_entries[0].hist.GetNbinsX()+1)]
    # calculate bin-by-bin fractions of bg / data, but only for pt < 1000 to avoid big errors
    # TODO include errors in decision?
    fracs = [(bg / data) if (data != 0 and get_x_bin_low_edge(data_entries, ind+1) < 1000) else 0
             for ind, (bg, data) in enumerate(zip(bg_bins, data_bins))]
    max_frac = max(fracs)
    max_frac_index = fracs.index(max_frac)
    # print(fracs)
    print("LARGEST BG/DATA:", max_frac)
    print("LARGEST BG/DATA BIN INDEX:", max_frac_index+1)
    print("LARGEST BG/DATA BIN:", get_x_bin_low_edge(data_entries, max_frac_index+1), "->", get_x_bin_low_edge(data_entries, max_frac_index+2))
    print('--'*80)


def make_data_mc_plot(entries, hist_name, x_label, output_filename, rebin=1,
                      do_logx=True, x_min=None, x_max=None,
                      do_logy=True, y_min=1E-1, y_max=1E6):
    """Make data-MC plot

    Automatically sorts contributions by integral.
    Entries with 'is_data' = True are not stacked, and plotted as points.
    Entriies with 'is_data' = False are stacked together, and plotted as filled areas

    Parameters
    ----------
    entries : list[dict]
        Things to be added. Each entry is a dict, that must have keys
        'tfile', 'label', 'is_data', 'style'.
        The value for 'style' is another dict, with options for line/marker/fill
    hist_name : str
        Name of hist object in TFiles to be plotted
    x_label : str
        x axis label
    output_filename : st
        Output filename
    rebin : int, optional
        To rebin hists
    do_logx : bool, optional
        Logarithmic x-axis
    x_min : None, optional
        x-axis minimum
    x_max : None, optional
        x-axis maximum
    do_logy : bool, optional
        Logarithmic y-axis
    y_min : float, optional
        y-axis minimum
    y_max : float, optional
        y-axis maximum
    """

    data_entries = []
    mc_entries = []

    for ent in entries:
        hist = hist=cu.get_from_tfile(ent['tfile'], hist_name)
        hist.SetLineColor(ent['style']['line_color'])
        hist.SetLineStyle(ent['style'].get('line_style', 1))

        hist.SetMarkerColor(ent['style']['marker_color'])
        hist.SetMarkerStyle(ent['style'].get('marker_style', 19))
        hist.SetMarkerSize(ent['style'].get('marker_size', 1))

        hist.SetFillColor(ent['style']['fill_color'])
        hist.SetFillStyle(ent['style'].get('fill_style', 1001))

        hist.Rebin(rebin)
        e = Entry(hist=hist,
                  label=ent['label'],
                  is_data=ent['is_data'],
                  is_bkg=ent.get('is_bkg', False if ent['is_data'] else True))
        if ent['is_data']:
            data_entries.append(e)
        else:
            mc_entries.append(e)

    # Sort so largest contribution first
    data_entries.sort(key=lambda entry: entry.integral, reverse=True)
    mc_entries.sort(key=lambda entry: entry.integral, reverse=True)

    for ent in data_entries:
        print(ent.label, ent.integral)
    for ent in mc_entries:
        print(ent.label, ent.integral)

    print_data_bg_stats(data_entries, mc_entries)

    canv = ROOT.TCanvas("c"+cu.get_unique_str(), "", 800, 600)
    canv.SetTicks(1, 1)
    if do_logx:
        canv.SetLogx()
    if do_logy:
        canv.SetLogy()
    canv.SetLeftMargin(0.12)

    hst = ROOT.THStack("hst", ";%s;N" % (x_label))
    for ent in mc_entries[::-1]:
        # Add in reverse order to add smallest first (bottom of stack)
        hst.Add(ent.hist)

    leg = ROOT.TLegend(0.7, 0.65, 0.88, 0.88)
    for ent in data_entries:
        leg.AddEntry(ent.hist, "%s (%g)" % (ent.label, ent.integral), "P")

    for ent in mc_entries:
        # Add largest first
        leg.AddEntry(ent.hist, "%s (%g)" % (ent.label, ent.integral), "F")

    hst.Draw("hist")
    hst.SetMaximum(y_max)
    hst.SetMinimum(y_min)
    if x_min and x_max:
        hst.GetXaxis().SetLimits(x_min, x_max)
    hst.GetYaxis().SetTitleOffset(1.2)
    for ent in data_entries:
        ent.hist.Draw("SAME E")
    leg.Draw()

    canv.cd()
    cms_latex = ROOT.TLatex()
    cms_latex.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    cms_latex.SetTextFont(42)
    cms_latex.SetTextSize(0.035)
    latex_height = 0.91
    start_x = 0.12
    cms_latex.DrawLatexNDC(start_x, latex_height, "#font[62]{CMS}#font[52]{ Preliminary}")
    cms_latex.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignBottom)
    cms_latex.DrawLatexNDC(0.95, latex_height, " 35.9 fb^{-1} (13 TeV)")

    canv.Modified()
    canv.Update()
    canv.SaveAs(output_filename)


COMPONENTS = [

{
    'tfile': tfile_data,
    'label': "Data",
    'is_data': True,
    'style': {
        'fill_color': ROOT.kBlack,
        'fill_style': 0,
        'marker_color': ROOT.kBlack,
        'marker_size': 1,
        'marker_style': 20,
        'line_color': ROOT.kBlack,
        'line_width': 0,
    }
},
{
    'tfile': tfile_dy,
    'label': "DY#rightarrowLL",
    'is_data': False,
    'is_bkg': False,
    'style': {
        'fill_color': ROOT.kAzure+6,
        'marker_color': ROOT.kAzure+6,
        'marker_size': 0,
        'line_color': ROOT.kAzure+6,
        'line_width': 0,
    }
},
{
    'tfile': tfile_ww,
    'label': "WW",
    'is_data': False,
    'style': {
        'fill_color': ROOT.kGreen+1,
        'marker_color': ROOT.kGreen+1,
        'marker_size': 0,
        'line_color': ROOT.kGreen+1,
        'line_width': 0,
    }
},
{
    'tfile': tfile_wz,
    'label': "WZ",
    'is_data': False,
    'style': {
        'fill_color': ROOT.kRed,
        'marker_color': ROOT.kRed,
        'marker_size': 0,
        'line_color': ROOT.kRed,
        'line_width': 0,
    }
},
{
    'tfile': tfile_zz,
    'label': "ZZ",
    'is_data': False,
    'style': {
        'fill_color': ROOT.kBlue,
        'marker_color': ROOT.kBlue,
        'marker_size': 0,
        'line_color': ROOT.kBlue,
        'line_width': 0,
    }
},
{
    'tfile': tfile_tt,
    'label': "t#bar{t}",
    'is_data': False,
    'style': {
        'fill_color': ROOT.kOrange,
        'marker_color': ROOT.kOrange,
        'marker_size': 0,
        'line_color': ROOT.kOrange,
        'line_width': 0,
    }
},
# {
#     'tfile': tfile_st_t,
#     'label': "Single-top t",
#     'is_data': False,
#     'style': {
#         'fill_color': ROOT.kCyan,
#         'marker_color': ROOT.kCyan,
#         'marker_size': 0,
#         'line_color': ROOT.kCyan,
#         'line_width': 0,
#     }
# },
# {
#     'tfile': tfile_st_tW,
#     'label': "Single-top tW",
#     'is_data': False,
#     'style': {
#         'fill_color': ROOT.kGray,
#         'marker_color': ROOT.kGray,
#         'marker_size': 0,
#         'line_color': ROOT.kGray,
#         'line_width': 0,
#     }
# },
# {
#     'tfile': tfile_st_s,
#     'label': "Single-top s",
#     'is_data': False,
#     'style': {
#         'fill_color': ROOT.kMagenta,
#         'marker_color': ROOT.kMagenta,
#         'marker_size': 0,
#         'line_color': ROOT.kMagenta,
#         'line_width': 0,
#     }
# },

]


make_data_mc_plot(COMPONENTS,
                  hist_name="ZPlusJets/pt_jet1",
                  x_label="p_{T}^{jet} [GeV]",
                  output_filename="workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_noZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut/zpj_ptJ_noKfactor.pdf",
                  rebin=1,
                  do_logx=True, x_min=20, x_max=3E3,
                  do_logy=True, y_min=1E-1, y_max=1E6)

make_data_mc_plot(COMPONENTS,
                  hist_name="ZPlusJets/pt_mumu",
                  x_label="p_{T}^{#mu#mu} [GeV]",
                  output_filename="workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_noZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut/zpj_ptZ_noKfactor.pdf",
                  rebin=1,
                  do_logx=True, x_min=10, x_max=3E3,
                  do_logy=True, y_min=1E-1, y_max=1E6)
