#!/usr/bin/env python


"""
Make Z+jet data-background plots
"""


from __future__ import print_function

from copy import deepcopy
import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
os.nice(10)

# My stuff
# from comparator import Contribution, Plot, grab_obj
import qg_common as qgc
import qg_general_plots as qgp
import common_utils as cu


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
# ROOT.gErrorIgnoreLevel = ROOT.kWarning


# zpj_dir = "workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet"
# tfile_data = cu.open_root_file("%s/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu.root" % (zpj_dir))
# tfile_ww = cu.open_root_file("%s/uhh2.AnalysisModuleRunner.MC.MC_WW.root" % (zpj_dir))
# tfile_wz = cu.open_root_file("%s/uhh2.AnalysisModuleRunner.MC.MC_WZ.root" % (zpj_dir))
# tfile_zz = cu.open_root_file("%s/uhh2.AnalysisModuleRunner.MC.MC_ZZ.root" % (zpj_dir))
# tfile_tt = cu.open_root_file("%s/uhh2.AnalysisModuleRunner.MC.MC_TTBAR.root" % (zpj_dir))
# tfile_st_t = cu.open_root_file("%s/uhh2.AnalysisModuleRunner.MC.MC_ST_T.root" % (zpj_dir))
# tfile_st_tw = cu.open_root_file("%s/uhh2.AnalysisModuleRunner.MC.MC_ST_TW.root" % (zpj_dir))
# tfile_st_s = cu.open_root_file("%s/uhh2.AnalysisModuleRunner.MC.MC_ST_S.root" % (zpj_dir))
# zpj_dir = "workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto"
# tfile_dy = cu.open_root_file("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root" % (zpj_dir))
# tfile_dy = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root")

zpj_dir = "workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_noZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut"
tfile_data = cu.open_root_file("%s/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu.root" % (zpj_dir))
tfile_ww = cu.open_root_file("%s/uhh2.AnalysisModuleRunner.MC.MC_WW.root" % (zpj_dir))
tfile_wz = cu.open_root_file("%s/uhh2.AnalysisModuleRunner.MC.MC_WZ.root" % (zpj_dir))
tfile_zz = cu.open_root_file("%s/uhh2.AnalysisModuleRunner.MC.MC_ZZ.root" % (zpj_dir))
tfile_tt = cu.open_root_file("%s/uhh2.AnalysisModuleRunner.MC.MC_TTBAR.root" % (zpj_dir))
# zpj_dir = "workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut"
# tfile_dy = cu.open_root_file("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root" % (zpj_dir))
tfile_dy = cu.open_root_file("workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root")


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


def _rescale_plot_labels(container, factor, left_margin):
    # Eurgh, why does ROOT scale all these sizes?
    container.GetXaxis().SetLabelSize(container.GetXaxis().GetLabelSize()/factor)
    container.GetXaxis().SetTitleSize(container.GetXaxis().GetTitleSize()/factor)
    container.GetXaxis().SetTitleOffset(container.GetXaxis().GetTitleOffset()*factor)  # doesn't seem to work?
    container.GetXaxis().SetTickLength(container.GetXaxis().GetTickLength()/factor)

    container.GetYaxis().SetLabelSize(container.GetYaxis().GetLabelSize()/factor)
    container.GetYaxis().SetTitleSize(container.GetYaxis().GetTitleSize()/factor)
    # magic numbers: 0.1 is the default margin, but scaling against that gives too much, so we knock it down by a bit
    container.GetYaxis().SetTitleOffset(container.GetYaxis().GetTitleOffset()*factor*(0.85*left_margin/0.1))
    # container.GetYaxis().SetTickLength(0.03/factor)


def make_data_mc_plot(entries, hist_name, x_label, output_filename, rebin=1,
                      title="",
                      do_logx=True, x_min=None, x_max=None,
                      do_logy=True, y_min=None, y_max=None,
                      do_compare_shapes=False):
    """Make data-MC plot with ratio subplot

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
    title : str, optional
        Optional title to put on plot
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
    do_compare_shapes : bool, optional
        Do not stack anything, just compare all shapes (i.e normalise each to unity)
    """

    data_entries = []
    mc_entries = []

    for ent in entries:
        # use histogram if it comes with it, otherwise get from tfile
        hist = ent.get('hist', None)
        if not hist:
            hist = cu.get_from_tfile(ent['tfile'], hist_name)
        if do_compare_shapes:
            if hist.Integral() > 0:
                hist.Scale(1./hist.Integral())
        hist.SetName(hist.GetName() + "_" + ent['label'])
        hist.SetLineColor(ent['style']['line_color'])
        hist.SetLineStyle(ent['style'].get('line_style', 1))

        hist.SetMarkerColor(ent['style']['marker_color'])
        hist.SetMarkerStyle(ent['style'].get('marker_style', 19))
        hist.SetMarkerSize(ent['style'].get('marker_size', 1))

        hist.SetFillColor(ent['style']['fill_color'])

        hist.SetFillStyle(ent['style'].get('fill_style', 1001))
        if do_compare_shapes:
            hist.SetFillStyle(0)
            hist.SetLineWidth(1)

        hist.Rebin(rebin)
        e = Entry(hist=hist,
                  label=ent['label'],
                  is_data=ent['is_data'],
                  is_bkg=ent.get('is_bkg', False if ent['is_data'] else True))
        if ent['is_data']:
            data_entries.append(e)
        else:
            mc_entries.append(e)

    if not do_compare_shapes:
        # Sort so largest contribution first
        data_entries.sort(key=lambda entry: entry.integral, reverse=True)
        mc_entries.sort(key=lambda entry: entry.integral, reverse=True)

    for ent in data_entries:
        print(ent.label, ent.integral)
    for ent in mc_entries:
        print(ent.label, ent.integral)

    print_data_bg_stats(data_entries, mc_entries)

    # Set up canvases etc
    canv = ROOT.TCanvas("c"+cu.get_unique_str(), "", 800, 600)
    canv.cd()
    right_margin = 0.03
    left_margin = 0.08 # use ROOT default
    top_margin = 0.1
    subplot_pad_height = 0.32
    subplot_pad_fudge = 0.01  # to get non-overlapping subplot axis
    main_pad = ROOT.TPad("main_pad", "", 0, subplot_pad_height+subplot_pad_fudge, 1, 1)
    ROOT.SetOwnership(main_pad, False)
    main_pad.SetTicks(1, 1)
    main_pad.SetBottomMargin(2*subplot_pad_fudge)
    main_pad.SetTopMargin(top_margin / (1-subplot_pad_height))
    main_pad.SetRightMargin(right_margin / (1-subplot_pad_height))
    main_pad.SetLeftMargin(left_margin / (1-subplot_pad_height))
    canv.cd()
    main_pad.Draw()
    subplot_pad = ROOT.TPad("subplot_pad", "", 0, 0, 1, subplot_pad_height-subplot_pad_fudge)
    ROOT.SetOwnership(subplot_pad, False)
    subplot_pad.SetTicks(1, 1)
    subplot_pad.SetFillColor(0)
    subplot_pad.SetFillStyle(0)
    subplot_pad.SetTopMargin(4*subplot_pad_fudge)
    subplot_pad.SetRightMargin(right_margin / (1-subplot_pad_height))
    subplot_pad.SetBottomMargin(0.35)
    subplot_pad.SetLeftMargin(left_margin / (1-subplot_pad_height))
    canv.cd()
    subplot_pad.Draw()

    if do_logx:
        main_pad.SetLogx()
        subplot_pad.SetLogx()
    if do_logy:
        main_pad.SetLogy()

    canv.cd()
    main_pad.cd()

    # Do main plot: MC stack, data, legend
    hst = ROOT.THStack("hst", ";%s;N" % (x_label))
    for ent in mc_entries[::-1]:
        # Add in reverse order to add smallest first (bottom of stack)
        hst.Add(ent.hist)

    leg = ROOT.TLegend(0.7, 0.65, 0.88, 0.88)
    for ent in data_entries:
        extra = "" if do_compare_shapes else " (%g)" % ent.integral
        leg.AddEntry(ent.hist, "%s%s" % (ent.label, extra), "P")

    for ent in mc_entries:
        # Add largest first
        extra = "" if do_compare_shapes else " (%g)" % ent.integral
        leg.AddEntry(ent.hist, "%s%s" % (ent.label, extra), "L" if do_compare_shapes else "F")

    if do_compare_shapes:
        hst.Draw("histe nostack")
    else:
        hst.Draw("hist")

    hst_stack = hst.GetStack()
    # Get sum hist - GetStack makes a cumulative stack
    mc_total = hst_stack.Last().Clone(hst_stack.Last().GetName()+"Ratio")

    if y_min is not None and y_max is not None:
        hst.SetMaximum(y_max)
        hst.SetMinimum(y_min)
    else:
        max_y = max(mc_total.GetMaximum(), max([ent.hist.GetMaximum() for ent in data_entries]))
        min_y = min(mc_total.GetMinimum(1E-14), min([ent.hist.GetMinimum(1E-14) for ent in data_entries]))
        print(min_y)
        if do_logy:
            hst.SetMaximum(10 * max_y)
            hst.SetMinimum(0.5 * min_y)
        else:
            hst.SetMaximum(2 * max_y)
            hst.SetMinimum(0.5 * min_y)

    if x_min is not None and x_max is not None:
        hst.GetXaxis().SetRangeUser(x_min, x_max)
        # hst.GetXaxis().SetLimits(x_min, x_max)  # SetLimits obeys you, SetRangeUser does some rounding?! But SetLimits makes poitns that don't align in the ratio plot
    # hst.GetYaxis().SetTitleOffset(1.2)
    for ent in data_entries:
        ent.hist.Draw("SAME E")

    canv.cd()
    leg.Draw()

    # Some text stuff
    cms_latex = ROOT.TLatex()
    cms_latex.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    cms_latex.SetTextFont(42)
    cms_latex.SetTextSize(0.035)
    latex_height = 0.915
    start_x = 0.12
    cms_latex.DrawLatexNDC(start_x, latex_height, "#font[62]{CMS}#font[52]{ Preliminary}")
    cms_latex.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignBottom)
    cms_latex.DrawLatexNDC(0.95, latex_height, " 35.9 fb^{-1} (13 TeV)")

    # Add title to plot
    text_latex = ROOT.TLatex()
    text_latex.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignTop)
    text_latex.SetTextFont(42)
    text_latex.SetTextSize(0.03)
    start_y = 0.87
    diff_y = 0.05
    if title != "":
        for ind, line in enumerate(title.split('\n')):
            text_latex.DrawLatex(start_x + 0.04, start_y - (ind*diff_y), line)

    _rescale_plot_labels(hst, 1-subplot_pad_height, left_margin)
    # Get rid of main plot x axis labels
    hst.GetHistogram().GetXaxis().SetLabelSize(0)
    hst.GetXaxis().SetLabelSize(0)
    hst.GetHistogram().GetXaxis().SetTitleSize(0)
    hst.GetXaxis().SetTitleSize(0)

    # Construct & plot subplot
    subplot_pad.cd()
    ratio_hists = []
    if do_compare_shapes:
        # compare all mc shapes to data
        ratio_hists = [ent.hist.Clone(ent.hist.GetName()+"Ratio") for ent in mc_entries]
        for ind, hist_ratio in enumerate(ratio_hists):
            hist_ratio.Divide(data_entries[0].hist)

            draw_opt = "HISTE"
            if ind > 0:
                draw_opt += " SAME"
            hist_ratio.Draw(draw_opt)

            hist_ratio.SetTitle(";%s;%s" % (hist_ratio.GetXaxis().GetTitle(), "Data / MC"))
            if x_min is not None and x_max is not None:
                hist_ratio.GetXaxis().SetRangeUser(x_min, x_max)

            hist_ratio.SetMinimum(0)  # use this, not SetRangeUser()
            hist_ratio.SetMaximum(2)  # use this, not SetRangeUser()

    else:
        # compare all data shapes to mc sum
        ratio_hists = [ent.hist.Clone(ent.hist.GetName()+"Ratio") for ent in data_entries]

        for ind, hist_ratio in enumerate(ratio_hists):
            hist_ratio.Divide(mc_total)
            # NB already styled from originals

            draw_opt = "E"
            if ind > 0:
                draw_opt += " SAME"
            hist_ratio.Draw(draw_opt)

            hist_ratio.SetTitle(";%s;%s" % (hist_ratio.GetXaxis().GetTitle(), "%s / MC" % (data_entries[ind].label)))
            if x_min is not None and x_max is not None:
                hist_ratio.GetXaxis().SetRangeUser(x_min, x_max)

            hist_ratio.SetMinimum(0.5)  # use this, not SetRangeUser()
            hist_ratio.SetMaximum(1.5)  # use this, not SetRangeUser()

    ratio_modifier = ratio_hists[0]

    # Draw a line at 1
    xax = ratio_modifier.GetXaxis()
    subplot_line = ROOT.TLine(xax.GetXmin(), 1., xax.GetXmax(), 1.)  # make sure you use whatever ROOT has decided the x limits should be
    subplot_line.SetLineStyle(2)
    subplot_line.SetLineWidth(2)
    subplot_line.SetLineColor(ROOT.kBlack)
    subplot_line.Draw()

    subplot_line_up = ROOT.TLine(xax.GetXmin(), 1.2, xax.GetXmax(), 1.2)
    subplot_line_up.SetLineStyle(3)
    subplot_line_up.SetLineWidth(2)
    subplot_line_up.SetLineColor(ROOT.kBlack)
    subplot_line_up.Draw()

    subplot_line_down = ROOT.TLine(xax.GetXmin(), 0.8, xax.GetXmax(), 0.8)
    subplot_line_down.SetLineStyle(3)
    subplot_line_down.SetLineWidth(2)
    subplot_line_down.SetLineColor(ROOT.kBlack)
    subplot_line_down.Draw()

    # Some resizing of subplot things
    _rescale_plot_labels(ratio_modifier, subplot_pad_height, left_margin)  # last factor is a fudge. no idea why
    ratio_modifier.GetXaxis().SetTitleOffset(ratio_modifier.GetXaxis().GetTitleOffset()*3)
    ratio_modifier.GetYaxis().SetNdivisions(505)

    canv.Modified()
    canv.Update()
    canv.SaveAs(output_filename)


def make_binned_data_mc_plots(entries, hist_name, x_label,
                              bins, bin_variable,
                              output_filename, rebin=1,
                              do_logx=False, x_min=None, x_max=None,
                              do_logy=False, y_min=None, y_max=None,
                              do_compare_shapes=False):
    """Do data/MC plots, but binned by some variable"""

    # Make a copy of entries with the big 2D hsit
    original_entries = deepcopy(entries)
    for ind, ent in enumerate(entries):
        # original_entries[ind]['tfile'] = ent['tfile']
        original_entries[ind]['hist'] = cu.get_from_tfile(ent['tfile'], hist_name)

    # For each bin, make a copy, but replace the hist value with the 1D projection for this bin
    for bin_ind, (bin_low, bin_high) in enumerate(bins):
        print("----- bin", bin_ind)
        bin_entries = deepcopy(original_entries)
        for ind, ent in enumerate(original_entries):
            bin_entries[ind]['hist'] = qgp.get_projection_plot(ent['hist'], bin_low, bin_high, 'y')
        title = "%s #in [%g, %g]" % (bin_variable, bin_low, bin_high)
        output_filename_stem, ext = os.path.splitext(output_filename)
        this_output_filename = "%s_%d%s" % (output_filename_stem, bin_ind, ext)
        make_data_mc_plot(bin_entries,
                          hist_name=hist_name,
                          x_label=x_label,
                          title=title,
                          output_filename=this_output_filename,
                          rebin=rebin,
                          do_logx=do_logx, x_min=x_min, x_max=x_max,
                          do_logy=do_logy, y_min=y_min, y_max=y_max,
                          do_compare_shapes=do_compare_shapes)


if __name__ == "__main__":
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
                      output_filename="%s/zpj_ptJ_Kfactor.pdf" % (zpj_dir),
                      rebin=5,
                      do_logx=True, x_min=20, x_max=3E3,
                      do_logy=True, y_min=1E-1, y_max=1E6)

    make_data_mc_plot(COMPONENTS,
                      hist_name="ZPlusJets/pt_jet1",
                      x_label="p_{T}^{jet} [GeV]",
                      output_filename="%s/zpj_ptJ_Kfactor_shapes.pdf" % (zpj_dir),
                      rebin=10,
                      do_compare_shapes=True,
                      do_logx=True, x_min=20, x_max=3E3,
                      do_logy=True)

    make_data_mc_plot(COMPONENTS,
                      hist_name="ZPlusJets/pt_jet1",
                      x_label="p_{T}^{jet} [GeV]",
                      output_filename="%s/zpj_ptJ_Kfactor_shapes_linY.pdf" % (zpj_dir),
                      rebin=10,
                      do_compare_shapes=True,
                      do_logx=True, x_min=20, x_max=3E3,
                      do_logy=False, y_min=0, y_max=0.4)

    make_data_mc_plot(COMPONENTS,
                      hist_name="ZPlusJets/pt_mumu",
                      x_label="p_{T}^{#mu#mu} [GeV]",
                      output_filename="%s/zpj_ptZ_Kfactor.pdf" % (zpj_dir),
                      rebin=2,
                      do_logx=True, x_min=10, x_max=3E3,
                      do_logy=True, y_min=1E-1, y_max=1E6)
    
    make_data_mc_plot(COMPONENTS,
                      hist_name="ZPlusJets/pt_mumu",
                      x_label="p_{T}^{#mu#mu} [GeV]",
                      output_filename="%s/zpj_ptZ_Kfactor_shapes.pdf" % (zpj_dir),
                      rebin=10,
                      do_compare_shapes=True,
                      do_logx=True, x_min=10, x_max=3E3,
                      do_logy=True)
    
    make_data_mc_plot(COMPONENTS,
                      hist_name="ZPlusJets/pt_mumu",
                      x_label="p_{T}^{#mu#mu} [GeV]",
                      output_filename="%s/zpj_ptZ_Kfactor_shapes_linY.pdf" % (zpj_dir),
                      rebin=10,
                      do_compare_shapes=True,
                      do_logx=True, x_min=10, x_max=3E3,
                      do_logy=False, y_min=0, y_max=0.3)

    # exit()
    # Do pt-binned plots
    make_binned_data_mc_plots(COMPONENTS,
                             hist_name="ZPlusJets/pt_jet1_z_ratio_vs_pt_z",
                             bins=qgc.PT_BINS_ZPJ,
                             bin_variable="p_{T}^{#mu#mu} [GeV]",
                             x_label="p_{T}^{jet 1} / p_{T}^{#mu#mu}",
                             output_filename="%s/zpj_ptJ_ptZ_ratio_binned_by_ptZ_Kfactor.pdf" % (zpj_dir),
                             rebin=2,
                             do_logx=False,
                             do_logy=True,
                             do_compare_shapes=False)

    make_binned_data_mc_plots(COMPONENTS[:-1],
                             hist_name="ZPlusJets/pt_jet1_z_ratio_vs_pt_z",
                             bins=qgc.PT_BINS_ZPJ,
                             bin_variable="p_{T}^{#mu#mu} [GeV]",
                             x_label="p_{T}^{jet 1} / p_{T}^{#mu#mu}",
                             output_filename="%s/zpj_ptJ_ptZ_ratio_binned_by_ptZ_Kfactor_shapes.pdf" % (zpj_dir),
                             rebin=2,
                             do_logx=False,
                             do_logy=True,
                             do_compare_shapes=True)

    make_binned_data_mc_plots(COMPONENTS,
                             hist_name="ZPlusJets/jet1_z_asym_vs_pt_z",
                             bins=qgc.PT_BINS_ZPJ,
                             bin_variable="p_{T}^{#mu#mu} [GeV]",
                             x_label="(p_{T}^{jet 1} - p_{T}^{#mu#mu}) / (p_{T}^{jet 1} + p_{T}^{#mu#mu})",
                             output_filename="%s/zpj_ptJ_ptZ_asym_binned_by_ptZ_Kfactor.pdf" % (zpj_dir),
                             rebin=2,
                             do_logx=False,
                             do_logy=True,
                             do_compare_shapes=False)

    make_binned_data_mc_plots(COMPONENTS[:-1],
                             hist_name="ZPlusJets/jet1_z_asym_vs_pt_z",
                             bins=qgc.PT_BINS_ZPJ,
                             bin_variable="p_{T}^{#mu#mu} [GeV]",
                             x_label="(p_{T}^{jet 1} - p_{T}^{#mu#mu}) / (p_{T}^{jet 1} + p_{T}^{#mu#mu})",
                             output_filename="%s/zpj_ptJ_ptZ_asym_binned_by_ptZ_Kfactor_shapes.pdf" % (zpj_dir),
                             rebin=2,
                             do_logx=False,
                             do_logy=True,
                             do_compare_shapes=True)

    make_binned_data_mc_plots(COMPONENTS,
                             hist_name="ZPlusJets/pt_jet1_z_ratio_vs_pt_jet1",
                             bins=qgc.PT_BINS_ZPJ,
                             bin_variable="p_{T}^{jet 1} [GeV]",
                             x_label="p_{T}^{jet 1} / p_{T}^{#mu#mu}",
                             output_filename="%s/zpj_ptJ_ptZ_ratio_binned_by_ptJ_Kfactor.pdf" % (zpj_dir),
                             rebin=2,
                             do_logx=False,
                             do_logy=True)


    make_binned_data_mc_plots(COMPONENTS[:-1],
                             hist_name="ZPlusJets/pt_jet1_z_ratio_vs_pt_jet1",
                             bins=qgc.PT_BINS_ZPJ,
                             bin_variable="p_{T}^{jet 1} [GeV]",
                             x_label="p_{T}^{jet 1} / p_{T}^{#mu#mu}",
                             output_filename="%s/zpj_ptJ_ptZ_ratio_binned_by_ptJ_Kfactor_shapes.pdf" % (zpj_dir),
                             rebin=2,
                             do_logx=False,
                             do_logy=True,
                             do_compare_shapes=True)

    make_binned_data_mc_plots(COMPONENTS,
                             hist_name="ZPlusJets/jet1_z_asym_vs_pt_jet1",
                             bins=qgc.PT_BINS_ZPJ,
                             bin_variable="p_{T}^{jet 1} [GeV]",
                             x_label="(p_{T}^{jet 1} - p_{T}^{#mu#mu}) / (p_{T}^{jet 1} + p_{T}^{#mu#mu})",
                             output_filename="%s/zpj_ptJ_ptZ_asym_binned_by_ptJ_Kfactor.pdf" % (zpj_dir),
                             rebin=2,
                             do_logx=False,
                             do_logy=True)
    
    make_binned_data_mc_plots(COMPONENTS[:-1],
                             hist_name="ZPlusJets/jet1_z_asym_vs_pt_jet1",
                             bins=qgc.PT_BINS_ZPJ,
                             bin_variable="p_{T}^{jet 1} [GeV]",
                             x_label="(p_{T}^{jet 1} - p_{T}^{#mu#mu}) / (p_{T}^{jet 1} + p_{T}^{#mu#mu})",
                             output_filename="%s/zpj_ptJ_ptZ_asym_binned_by_ptJ_Kfactor_shapes.pdf" % (zpj_dir),
                             rebin=2,
                             do_logx=False,
                             do_logy=True,
                             do_compare_shapes=True)

