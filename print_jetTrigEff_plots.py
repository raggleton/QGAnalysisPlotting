#!/usr/bin/env python

"""Print main QG plots for a given sample.

Any other plots comparing more than one sample should go into its own script!
"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
from collections import OrderedDict
import sys
from array import array
from bisect import bisect_left
from copy import deepcopy

# My stuff
# from comparator import Contribution, Plot, grab_obj
# import qg_common as qgc
# import qg_general_plots as qgg

# For debugging
# import sys
# import tracers
# sys.settrace(tracers.trace_calls)
# sys.settrace(tracers.trace_calls_detail)
# sys.settrace(tracers.trace_calls_and_returns)

# import hunter
# hunter.trace(module='comparator')

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit()
# ROOT.gStyle.SetStatColor()

# Control plot output format
OUTPUT_FMT = "pdf"


trig_info = OrderedDict()
trig_info['HLT_PFJet40'] = {
    'threshold': 40.,
    'prescale': 135426.215343,
    'color': ROOT.kRed
}
trig_info['HLT_PFJet60'] = {
    'threshold': 60.,
    'prescale': 49891.9453547,
    'color': ROOT.kBlue,
    'fit_all': (lambda isFatJet: not isFatJet)
}
trig_info['HLT_PFJet80'] = {
    'threshold': 80.,
    'prescale': 13120.4895678,
    'color': ROOT.kGreen+2,
    'fit_all': (lambda isFatJet: not isFatJet)
}
trig_info['HLT_PFJet140'] = {
    'threshold': 140.,
    'prescale': 1496.44452961,
    'color': ROOT.kViolet+5,
    'fit_all': (lambda isFatJet: not isFatJet)
}
trig_info['HLT_PFJet200'] = {
    'threshold': 200.,
    'prescale': 348.686346954,
    'color': ROOT.kOrange,
    'fit_all': (lambda isFatJet: not isFatJet)
}
trig_info['HLT_PFJet260'] = {
    'threshold': 260.,
    'prescale': 61.0210313345,
    'color': ROOT.kTeal,
    'fit_all': (lambda isFatJet: not isFatJet)
}
trig_info['HLT_PFJet320'] = {
    'threshold': 320.,
    'prescale': 20.446914767,
    'color': ROOT.kViolet,
    'fit_all': (lambda isFatJet: not isFatJet)
}
trig_info['HLT_PFJet400'] = {
    'threshold': 400.,
    'prescale': 3*2.38456,
    'color': ROOT.kOrange-6,
    'fit_all': (lambda isFatJet: not isFatJet)
}
trig_info['HLT_PFJet450'] = {
    'threshold': 450.,
    'prescale': 1.00010464076,
    'color': ROOT.kAzure+1,
    'fit_all': (lambda isFatJet: False) 
}
trig_info['HLT_PFJet500'] = {
    'threshold': 500.,
    'prescale': 1.00010464076,
    'color': ROOT.kSpring-8,
    'fit_all': (lambda isFatJet: False)
}

all_results = OrderedDict()

def do_custom_rebin(hist, newname, lower_limit, factor):
    """Makes a rebinned histogram above lower_limit by grouping together bins with given factor"""
    print("custom rebin:", lower_limit, factor)
    if factor == 1:
        return hist.Clone(newname)

    nbins = hist.GetNbinsX()
    bins = [hist.GetBinLowEdge(i) for i in range(1, nbins+2)]

    # figure out sensible lower_limit if not in list of bin edges
    if lower_limit not in bins:
        print('WARNING: lower_limit not found in bin edges')
        # find the closest bin edge to the desired value
        ind = bisect_left(bins, lower_limit)
        if ind == 0:
            lower_limit = bins[0]
        if ind == len(bins):
            lower_limit = bins[-1]
        lower = bins[ind-1]
        higher = bins[ind]
        if (lower_limit-lower) < (higher - lower_limit):
            lower_limit = lower
        else:
            lower_limit = higher
        print("Adjusted lower_limit to nearest value =", lower_limit)

    # ensure integer multiple of factor bins to be regrouped
    rebin_remainder = (nbins-bins.index(lower_limit)) % factor
    if rebin_remainder != 0:
        print("WARNING: factor must be a divisor with no remainder. nbins: ", nbins-bins.index(lower_limit), "factor:", factor)
        lower_limit = bins[bins.index(lower_limit)+rebin_remainder]
        print("Will adjust lower_limit to higher value to make this so. New lower_limit = ", lower_limit)

    lower_limit_ind = bins.index(lower_limit)
    # original bins at low x
    new_bins = bins[:lower_limit_ind]
    # regrouped bins at higher x
    new_bins += [bins[i] for i in range(lower_limit_ind, nbins+2, factor)]
    hnew = hist.Rebin(len(new_bins)-1, newname, array('d', new_bins))
    return hnew


def do_trig_plots(input_filename, output_dir, title="", eta_min=-2.4, eta_max=2.4, append=""):
    """Do efficiencies and fits for all triggers"""
    if not os.path.isfile(input_filename):
        raise IOError("No input file", input_filename)

    f = ROOT.TFile(input_filename)
    dir_name = "PFJet"

    is_fat_jet = "AK8" in input_filename.upper()

    # Get singlemu hist
    h_all = f.Get(dir_name + "/pt_vs_eta_all")

    # Figure out bins for eta edges
    eta_min_bin, eta_max_bin = 0, h_all.GetNbinsY()+1
    yax = h_all.GetYaxis()

    eta_min_bin = yax.FindBin(eta_min)
    eta_max_bin = yax.FindBin(eta_max)
    # don't want to include the upper bin if the value is at the low edge
    if yax.GetBinLowEdge(eta_max_bin) == eta_max:
        eta_max_bin -= 1

    h_all_pt = h_all.ProjectionX("allPT", eta_min_bin, eta_max_bin)
    h_all_pt.Sumw2()
    h_all_pt.SetFillColor(17)
    h_all_pt.SetLineColor(17)

    this_trig_info = deepcopy(trig_info)
    for name, info in this_trig_info.iteritems():
        # for each trig, jet 2d pt vs eta hist, project into 1D pt distribution
        # then create efficiency hist using single mu hist
        rebin_factor = 2 if info['threshold'] > 100 else 2  # rough rebinning across all pt
        h2d = f.Get(dir_name + "/pt_vs_eta_%s_v*" % name)
        info['h2d'] = h2d
        info['hpt'] = h2d.ProjectionX(name+"PT", eta_min_bin, eta_max_bin)
        info['hpt'].Sumw2()
        info['hpt'].Scale(info['prescale']/list(this_trig_info.values())[-1]['prescale'])  # normalise it
        info['hpt'].SetLineColor(info['color'])
        info['hpt'].SetMarkerColor(info['color'])

        this_hpt = info['hpt'].Clone(info['hpt'].GetName()+"Clone").Rebin(rebin_factor)

        # rebin the jet pt hist at higher pt where it plateaus
        higher_pt_rebin_factor = 10
        higher_pt_rebin_limit = info['threshold'] * 1.4
        this_hpt_rebin = do_custom_rebin(this_hpt, this_hpt.GetName()+"CustomRebin", higher_pt_rebin_limit, higher_pt_rebin_factor)

        # rebin the muon pt hist in same way
        this_h_all_pt = h_all_pt.Clone(h_all_pt.GetName()+name)
        this_h_all_pt.Rebin(rebin_factor)
        h_all_pt_rebin = do_custom_rebin(this_h_all_pt, this_h_all_pt.GetName()+"CustomRebin", higher_pt_rebin_limit, higher_pt_rebin_factor)

        info['heff'] = this_hpt_rebin.Clone(this_hpt_rebin.GetName() + "Eff")
        info['heff'].Sumw2()
        info['heff'].Divide(info['heff'], h_all_pt_rebin, 1, 1, "B")
        new_title = info['heff'].GetTitle().replace("_v*", "")+" "+append
        info['heff'].SetTitle(new_title+";Leading jet p_{T} [GeV];#epsilon")
        # info['heff'] = ROOT.TEfficiency(info['hpt'], h_all_pt)  # cant use as > 1 due to prescaling

    # return
    # plot pt distributions
    hst = ROOT.THStack("hst", append+";Jet p_{T} [GeV];N")
    leg = ROOT.TLegend(0.5, 0.5, 0.88, 0.88)

    cms_text = ROOT.TPaveText(0.14, 0.9, 0.4, 0.92, "NDC")
    cms_text.AddText("CMS Preliminary 35.864 fb^{-1}")
    cms_text.SetFillStyle(0)
    cms_text.SetBorderSize(0)
    cms_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
    cms_text.SetTextFont(63)
    cms_text.SetTextSize(18)

    jet_text = ROOT.TPaveText(0.6, 0.9, 0.9, 0.92, "NDC")
    jet_text.AddText(title)
    jet_text.SetFillStyle(0)
    jet_text.SetBorderSize(0)
    jet_text.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignBottom)
    jet_text.SetTextFont(63)
    jet_text.SetTextSize(18)

    rebin_factor = 5
    hst.Add(h_all_pt.Rebin(rebin_factor))
    leg.AddEntry(h_all_pt, "HLT_IsoMu24 || HLT_IsoTkMu24" , "F")
    for name, info in this_trig_info.iteritems():
        hst.Add(info['hpt'].Rebin(rebin_factor))
        leg.AddEntry(info['hpt'], name, "L")

    c = ROOT.TCanvas("c1", "", 800, 600)
    c.SetTicks(1, 1)
    # c.SetLogx()
    hst.Draw("HISTE NOSTACK")
    # hst.GetXaxis().SetRangeUser(0, 600)
    hst.SetMinimum(10**-1)
    leg.Draw()
    cms_text.Draw()
    jet_text.Draw()
    c.SetLogy()
    c.SaveAs(output_dir + "/pt_trig_%s.%s" % (append, OUTPUT_FMT))

    # plot effs
    for name, info in this_trig_info.iteritems():
        c = ROOT.TCanvas("ceff"+name, "", 800, 600)
        c.SetTicks(1, 1)
        # c.SetLogy()

        info['heff'].SetMarkerStyle(22)
        # info['heff'].SetTitle(name)
        info['heff'].SetMaximum(1.5)
        info['heff'].SetMinimum(0)
        info['heff'].GetXaxis().SetRangeUser(0, min(6*info['threshold'], 2000))
        info['heff'].Draw()

        # Do fit
        lower_threshold = info['threshold']/3.
        higher_threshold = info['threshold']*3.
        eff_fit = ROOT.TF1("eff_%s" % name, '[3]*([0] + 0.5 * (1-[0]) * (1 + erf((x-[1])/[2])))', lower_threshold, higher_threshold)
        eff_fit.SetParName(0, 'a')
        eff_fit.SetParName(1, 'mu')
        eff_fit.SetParName(2, 'sigma')
        eff_fit.SetParName(3, 'N')
        eff_fit.SetLineColor(ROOT.kBlack)
        eff_fit.SetLineWidth(1)
        eff_fit.SetParameter('a', 0)
        eff_fit.SetParameter('mu', info['threshold'])
        eff_fit.SetParameter('sigma', info['threshold']/10)
        eff_fit.SetParameter('N', 1)
        eff_fit.SetNpx(5000)
        fit_result = info['heff'].Fit(eff_fit, 'RSEM')
        info['heff'].Draw("")

        ROOT.gPad.Modified()
        ROOT.gPad.Update()

        # Update fit by increasing lower limit to really capture the high efficiency region
        if not info.get("fit_all", lambda x: True)(is_fat_jet):
            fit_factor = 0.4 if is_fat_jet else 0.75
            eff_fit.SetRange(eff_fit.GetX(fit_factor*eff_fit.GetParameter("N")), higher_threshold*1.)
            fit_result = info['heff'].Fit(eff_fit, 'RSEM')
            ROOT.gPad.Modified()
            ROOT.gPad.Update()

        # Draw fit stats
        stats_box = info['heff'].FindObject("stats")
        stats_box.SetFillColor(ROOT.kWhite)
        stats_box.SetBorderSize(0)
        stats_box.SetFillStyle(0)
        stats_box.SetX1NDC(0.62)
        stats_box.SetX2NDC(0.88)
        stats_box.SetY1NDC(0.75)
        stats_box.SetY2NDC(0.88)

        # Add in info about 99% relative efficiency
        good_eff = 0.99 * eff_fit.GetParameter("N")
        good_eff_pt = eff_fit.GetX(good_eff)
        info['good_eff_pt'] = good_eff_pt
        eff_text = ROOT.TPaveText(0.63, 0.65, 0.88, 0.73, "NDC")
        eff_text.AddText("#epsilon = 0.99 #times %.3f" % eff_fit.GetParameter("N"))
        eff_text.AddText("@ p_{T} = %3.f GeV" % good_eff_pt)
        eff_text.SetFillStyle(0)
        eff_text.SetBorderSize(0)
        eff_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        eff_text.Draw()

        cms_text.Draw()
        jet_text.Draw()

        c.SaveAs(output_dir + "/eff_%s_%s.%s" % (name, append, OUTPUT_FMT))

    # make graph of fully efficiency pt vs threshold
    thresholds = [info['threshold'] for info in this_trig_info.itervalues()]
    fully_eff_pt = [info['good_eff_pt'] for info in this_trig_info.itervalues()]
    gr = ROOT.TGraph(len(thresholds), array('d', thresholds), array('d', fully_eff_pt))
    c = ROOT.TCanvas("cgr", "", 800, 600)
    c.SetTicks(1, 1)
    gr.SetTitle(";Trigger threshold [GeV];99% efficiency p_{T} [GeV]")
    gr.SetMarkerStyle(20)
    # do a pol1 fit
    thres_fit = ROOT.TF1("f1", "pol1", thresholds[0], thresholds[-1])
    thres_fit.SetLineColor(ROOT.kRed)
    thres_fit.SetLineWidth(1)
    thres_fit.SetLineStyle(2)
    status = gr.Fit(thres_fit, "RSE")
    gr.Draw("ALP")
    c.Modified()
    c.Update()
    stats_box = gr.FindObject("stats")
    stats_box.SetFillColor(ROOT.kWhite)
    stats_box.SetBorderSize(0)
    stats_box.SetFillStyle(0)
    stats_box.SetX1NDC(0.62)
    stats_box.SetX2NDC(0.88)
    stats_box.SetY1NDC(0.25)
    stats_box.SetY2NDC(0.38)
    cms_text.Draw()
    jet_text.Draw()
    c.SaveAs(output_dir + "/fully_eff_pt_vs_threshold_%s.%s" % (append, OUTPUT_FMT))

    return this_trig_info


def do_plots_and_comparisons(inputs):
    regions = (
        [-2.4, -1.6, "endcapMinus"],
        [1.6, 2.4, "endcapPlus"],
        [-1.6, 1.6, "barrel"]
    )
    for eta_min, eta_max, append in regions:
        for filename, title in do_these:
            results = do_trig_plots(filename, os.path.dirname(filename), title, eta_min, eta_max, append)
            all_results[title] = results

        mg = ROOT.TMultiGraph()
        c = ROOT.TCanvas("cmg"+append, "", 800, 600)
        c.SetTicks(1, 1)
        leg = ROOT.TLegend(0.7, 0.2, 0.88, 0.38)
        colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kOrange+1]
        for ind, (name, result) in enumerate(all_results.iteritems()):
            thresholds = [info['threshold'] for info in result.itervalues()]
            fully_eff_pt = [info['good_eff_pt'] for info in result.itervalues()]
            g = ROOT.TGraph(len(thresholds), array('d', thresholds), array('d', fully_eff_pt))
            g.SetMarkerColor(colors[ind])
            g.SetMarkerStyle(20+ind)
            g.SetLineColor(colors[ind])
            mg.Add(g)
            leg.AddEntry(g, name, "LP")
        mg.SetTitle(";Trigger threshold [GeV];99% efficiency p_{T} [GeV]")
        mg.Draw("ALP")
        leg.Draw()
        cms_text = ROOT.TPaveText(0.14, 0.9, 0.4, 0.92, "NDC")
        cms_text.AddText("CMS Preliminary 35.864 fb^{-1}")
        cms_text.SetFillStyle(0)
        cms_text.SetBorderSize(0)
        cms_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        cms_text.SetTextFont(63)
        cms_text.SetTextSize(18)
        cms_text.Draw()
        c.SaveAs("comparingTriggers_%s.%s" % (append, OUTPUT_FMT))


if __name__ == "__main__":
    do_these = [
        ('workdir_ak4chs_jettrig/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu_JetTrig.root', 'AK4 CHS'),
        ('workdir_ak4puppi_jettrig/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu_JetTrig.root', 'AK4 PUPPI'),
        ('workdir_ak8chs_jettrig/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu_JetTrig.root', 'AK8 CHS'),
        ('workdir_ak8puppi_jettrig/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu_JetTrig.root', 'AK8 PUPPI'),
    ]
    do_plots_and_comparisons(do_these)
