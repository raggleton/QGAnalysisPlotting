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
    'color': ROOT.kBlue
}
trig_info['HLT_PFJet80'] = {
    'threshold': 80.,
    'prescale': 13120.4895678,
    'color': ROOT.kGreen+2
}
trig_info['HLT_PFJet140'] = {
    'threshold': 140.,
    'prescale': 1496.44452961,
    'color': ROOT.kViolet+5
}
trig_info['HLT_PFJet200'] = {
    'threshold': 200.,
    'prescale': 348.686346954,
    'color': ROOT.kOrange
}
trig_info['HLT_PFJet260'] = {
    'threshold': 260.,
    'prescale': 61.0210313345,
    'color': ROOT.kTeal
}
trig_info['HLT_PFJet320'] = {
    'threshold': 320.,
    'prescale': 20.446914767,
    'color': ROOT.kViolet
}
trig_info['HLT_PFJet400'] = {
    'threshold': 400.,
    'prescale': 2.38456,
    'color': ROOT.kOrange-6
}
trig_info['HLT_PFJet450'] = {
    'threshold': 450.,
    'prescale': 1.00010464076,
    'color': ROOT.kAzure+1
}
trig_info['HLT_PFJet500'] = {
    'threshold': 500.,
    'prescale': 1.00010464076,
    'color': ROOT.kSpring-8
}


def do_custom_rebin(hist, newname, lower_limit, factor):
    """Makes a rebinned histogram above lower_limit by grouping together bins with given factor"""
    print "custom rebin:", lower_limit, factor
    nbins = hist.GetNbinsX()
    bins = [hist.GetBinLowEdge(i) for i in range(1, nbins+2)]

    # figure out sensible lower_limit if not in list of bin edges
    if lower_limit not in bins:
        print 'WARNING: lower_limit not found in bin edges'
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
        print "Adjusted lower_limit to nearest value =", lower_limit

    # ensure integer multiple of factor bins to be regrouped
    rebin_remainder = (nbins-bins.index(lower_limit)) % factor
    if rebin_remainder != 0:
        print "WARNING: factor must be a divisor with no remainder. nbins: ", nbins-bins.index(lower_limit), "factor:", factor
        lower_limit = bins[bins.index(lower_limit)+rebin_remainder]
        print "Will adjust lower_limit to higher value to make this so. New lower_limit = ", lower_limit

    lower_limit_ind = bins.index(lower_limit)
    # original bins at low x
    new_bins = bins[:lower_limit_ind]
    # regrouped bins at higher x
    new_bins += [bins[i] for i in range(lower_limit_ind, nbins+2, factor)]
    hnew = hist.Rebin(len(new_bins)-1, newname, array('d', new_bins))
    return hnew


def do_trig_plots(input_filename, output_dir, title=""):
    """Do efficiencies and fits for all triggers"""
    if not os.path.isfile(input_filename):
        raise IOError("No input file", input_filename)

    f = ROOT.TFile(input_filename)
    dir_name = "PFJet"

    # Get singlemu hist
    h_all = f.Get(dir_name + "/pt_vs_eta_all")

    # Figure out bins for eta edges
    eta_min_bin, eta_max_bin = 0, h_all.GetNbinsY()+1
    yax = h_all.GetYaxis()
    eta_min, eta_max = -2.4, 2.4
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
        rebin_factor = 4 if info['threshold'] > 100 else 1  # rough rebinning across all pt
        h2d = f.Get(dir_name + "/pt_vs_eta_%s_v*" % name)
        info['h2d'] = h2d
        info['hpt'] = h2d.ProjectionX(name+"PT", eta_min_bin, eta_max_bin)
        info['hpt'].Sumw2()
        info['hpt'].Scale(info['prescale']/this_trig_info.values()[-1]['prescale'])  # normalise it
        info['hpt'].SetLineColor(info['color'])
        info['hpt'].SetMarkerColor(info['color'])

        this_hpt = info['hpt'].Clone(info['hpt'].GetName()+"Clone").Rebin(rebin_factor)

        # rebin the jet pt hist at higher pt where it plateaus
        higher_pt_rebin_factor = 5
        higher_pt_rebin_limit = info['threshold'] * 1.4
        this_hpt_rebin = do_custom_rebin(this_hpt, this_hpt.GetName()+"CustomRebin", higher_pt_rebin_limit, higher_pt_rebin_factor)

        # rebin the muon pt hist in same way
        this_h_all_pt = h_all_pt.Clone(h_all_pt.GetName()+name)
        this_h_all_pt.Rebin(rebin_factor)
        h_all_pt_rebin = do_custom_rebin(this_h_all_pt, this_h_all_pt.GetName()+"CustomRebin", higher_pt_rebin_limit, higher_pt_rebin_factor)

        info['heff'] = this_hpt_rebin.Clone(this_hpt_rebin.GetName() + "Eff")
        info['heff'].Divide(h_all_pt_rebin)
        info['heff'].SetTitle(info['heff'].GetTitle()+";Leading jet p_{T} [GeV];#epsilon")
        # info['heff'] = ROOT.TEfficiency(info['hpt'], h_all_pt)  # cant use as > 1 due to prescaling

    # return
    # plot pt distributions
    hst = ROOT.THStack("hst", ";Jet p_{T} [GeV];N")
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

    rebin_factor = 4
    hst.Add(h_all_pt.Rebin(rebin_factor))
    leg.AddEntry(h_all_pt, "HLT_IsoMu24 || HLT_IsoTkMu24" , "F")
    for name, info in this_trig_info.iteritems():
        hst.Add(info['hpt'].Rebin(rebin_factor))
        leg.AddEntry(info['hpt'], name, "L")


    c = ROOT.TCanvas("c1", "", 800, 600)
    c.SetTicks(1, 1)
    hst.Draw("HISTE NOSTACK")
    # hst.GetXaxis().SetRangeUser(0, 600)
    hst.SetMinimum(10**-1)
    leg.Draw()
    cms_text.Draw()
    jet_text.Draw()
    c.SetLogy()
    c.SaveAs(output_dir + "/pt_trig.%s" % (OUTPUT_FMT))

    # plot effs
    for name, info in this_trig_info.iteritems():
        c = ROOT.TCanvas("ceff"+name, "", 800, 600)
        c.SetTicks(1, 1)
        # c.SetLogy()

        info['heff'].SetMarkerStyle(22)
        info['heff'].SetTitle(name)
        info['heff'].SetMaximum(1.5)
        info['heff'].SetMinimum(0)
        info['heff'].GetXaxis().SetRangeUser(0, min(10*info['threshold'], 2000))
        info['heff'].Draw()

        # Do fit
        eff_fit = ROOT.TF1("eff_%s" % name, '[3]*([0] + 0.5 * (1-[0]) * (1 + erf((x-[1])/[2])))', info['threshold']/3., info['threshold']*3.)
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
        fit_result = info['heff'].Fit(eff_fit, 'VRSEM')
        info['heff'].Draw()

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

        c.SaveAs(output_dir + "/eff_%s.%s" % (name, OUTPUT_FMT))

    # make graph of fully efficiency pt vs threshold
    thresholds = [info['threshold'] for info in this_trig_info.itervalues()]
    fully_eff_pt = [info['good_eff_pt'] for info in this_trig_info.itervalues()]
    gr = ROOT.TGraph(len(thresholds), array('d', thresholds), array('d', fully_eff_pt))
    c = ROOT.TCanvas("cgr", "", 800, 600)
    c.SetTicks(1, 1)
    gr.SetTitle(";Trigger threshold [GeV];99% efficiency p_{T} [GeV]")
    gr.SetMarkerStyle(20)
    gr.Draw("ALP")
    cms_text.Draw()
    jet_text.Draw()
    c.SaveAs(output_dir + "/fully_eff_pt_vs_threshold.%s" % OUTPUT_FMT)

    return this_trig_info


if __name__ == "__main__":
    do_these = [
        ('workdir_ak4chs_jettrig/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu_JetTrig.root', 'AK4 CHS'),
        ('workdir_ak8chs_jettrig/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu_JetTrig.root', 'AK8 CHS'),
        ('workdir_ak4puppi_jettrig/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu_JetTrig.root', 'AK4 PUPPI'),
        ('workdir_ak8puppi_jettrig/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu_JetTrig.root', 'AK8 PUPPI'),
    ]

    for filename, title in do_these:
        do_trig_plots(filename, os.path.dirname(filename), title)

