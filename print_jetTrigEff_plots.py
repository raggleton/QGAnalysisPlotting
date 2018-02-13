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

# My stuff
from comparator import Contribution, Plot, grab_obj
import qg_common as qgc
import qg_general_plots as qgg

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
ROOT.gStyle.SetOptFit()
# ROOT.gStyle.SetStatColor()

# Control plot output format
OUTPUT_FMT = "pdf"


trig_info = OrderedDict()
trig_info['HLT_PFJet40'] = {
    'prescale': 135426.215343,
    'color': ROOT.kRed
}
trig_info['HLT_PFJet60'] = {
    'prescale': 49891.9453547,
    'color': ROOT.kBlue
}
trig_info['HLT_PFJet80'] = {
    'prescale': 13120.4895678,
    'color': ROOT.kGreen+2
}
trig_info['HLT_PFJet140'] = {
    'prescale': 1496.44452961,
    'color': ROOT.kViolet+5
}
trig_info['HLT_PFJet200'] = {
    'prescale': 348.686346954,
    'color': ROOT.kOrange
}
trig_info['HLT_PFJet260'] = {
    'prescale': 61.0210313345,
    'color': ROOT.kTeal
}
trig_info['HLT_PFJet320'] = {
    'prescale': 20.446914767,
    'color': ROOT.kViolet
}
trig_info['HLT_PFJet400'] = {
    'prescale': 2.38456,
    'color': ROOT.kOrange-6
}
trig_info['HLT_PFJet450'] = {
    'prescale': 1.00010464076,
    'color': ROOT.kAzure+1
}
trig_info['HLT_PFJet500'] = {
    'prescale': 1.00010464076,
    'color': ROOT.kSpring-9
}


def do_trig_plots(input_filename, output_dir):
    f = ROOT.TFile(input_filename)
    dir_name = "PFJet"
    
    rebin_factor = 2

    # Get singlemu hist
    h_all = f.Get(dir_name + "/pt_vs_eta_all")
    h_all_pt = h_all.ProjectionX("allPT").Rebin(rebin_factor)
    h_all_pt.SetFillColor(17)
    h_all_pt.SetLineColor(17)

    for name in trig_info:
        # for each trig, jet 2d pt vs eta hist, project into 1D pt distribution
        # then create efficiency hist using single mu hist
        h2d = f.Get(dir_name + "/pt_vs_eta_%s_v*" % name)
        trig_info[name]['h2d'] = h2d
        trig_info[name]['hpt'] = h2d.ProjectionX(name+"PT")
        trig_info[name]['hpt'].Sumw2()
        trig_info[name]['hpt'].Scale(trig_info[name]['prescale']/trig_info.values()[-1]['prescale'])  # normalise it
        trig_info[name]['hpt'].SetLineColor(trig_info[name]['color'])
        trig_info[name]['hpt'].SetMarkerColor(trig_info[name]['color'])
        trig_info[name]['hpt'].Rebin(rebin_factor)
        trig_info[name]['heff'] = trig_info[name]['hpt'].Clone(trig_info[name]['hpt'].GetName() + "Eff")
        trig_info[name]['heff'].Divide(h_all_pt)
        trig_info[name]['heff'].SetTitle(trig_info[name]['heff'].GetTitle()+";Leading jet p_{T} [GeV];#epsilon")
        # trig_info[name]['heff'] = ROOT.TEfficiency(trig_info[name]['hpt'], h_all_pt)  # cant use as > 1 due to prescaling

    # plot pt distributions
    hst = ROOT.THStack("hst", ";Jet p_{T} [GeV];N")
    leg = ROOT.TLegend(0.5, 0.5, 0.88, 0.88)
    cms_text = ROOT.TPaveText(0.6, 0.9, 0.9, 0.92, "NDC")
    cms_text.AddText("CMS Preliminary 35.864 fb^{-1}")
    cms_text.SetFillStyle(0)
    cms_text.SetBorderSize(0)
    cms_text.SetTextAlign(ROOT.kHAlignRight + ROOT.kVAlignBottom)
    cms_text.SetTextFont(63)
    cms_text.SetTextSize(18)
    hst.Add(h_all_pt)
    leg.AddEntry(h_all_pt, "HLT_IsoMu24 || HLT_IsoTkMu24" , "F")
    for name, info in trig_info.iteritems():
        hst.Add(info['hpt'])
        leg.AddEntry(info['hpt'], name, "L")


    c = ROOT.TCanvas("c1", "", 800, 600)
    c.SetTicks(1, 1)
    hst.Draw("HISTE NOSTACK")
    # hst.GetXaxis().SetRangeUser(0, 600)
    hst.SetMinimum(10**-1)
    leg.Draw()
    cms_text.Draw()
    c.SetLogy()
    c.SaveAs(output_dir + "/pt_trig.%s" % (OUTPUT_FMT))

    # plot effs
    for name, info in trig_info.iteritems():
        c = ROOT.TCanvas("ceff"+name, "", 800, 600)
        c.SetTicks(1, 1)
        # c.SetLogy()
        trig_value = float(name.replace("HLT_PFJet", ''))
        
        info['heff'].SetMarkerStyle(22)
        info['heff'].SetTitle(name)
        info['heff'].SetMaximum(1.5)
        info['heff'].SetMinimum(0)
        info['heff'].GetXaxis().SetRangeUser(0, min(10*trig_value, 2000))
        info['heff'].Draw()

        # Do fit
        eff_fit = ROOT.TF1("eff_%s" % name, '[3]*([0] + 0.5 * (1-[0]) * (1 + erf((x-[1])/[2])))', trig_value/3., trig_value*3.)
        eff_fit.SetParName(0, 'a')
        eff_fit.SetParName(1, 'mu')
        eff_fit.SetParName(2, 'sigma')
        eff_fit.SetParName(3, 'N')
        eff_fit.SetLineColor(ROOT.kBlack)
        eff_fit.SetLineWidth(1)
        eff_fit.SetParameter('a', 0)
        eff_fit.SetParameter('mu', trig_value)
        eff_fit.SetParameter('sigma', trig_value/10)
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
        good_eff = eff_fit.GetX(0.99*eff_fit.GetParameter("N"))
        eff_text = ROOT.TPaveText(0.63, 0.65, 0.88, 0.73, "NDC")
        eff_text.AddText("#epsilon = 0.99 #times %.3f" % eff_fit.GetParameter("N"))
        eff_text.AddText("@ p_{T} = %3.f GeV" % good_eff)
        eff_text.SetFillStyle(0)
        eff_text.SetBorderSize(0)
        eff_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        eff_text.Draw()
        cms_text.Draw()

        c.SaveAs(output_dir + "/eff_%s.%s" % (name, OUTPUT_FMT))


if __name__ == "__main__":
    for filename in sys.argv[1:]:
        do_trig_plots(filename, os.path.dirname(filename))

