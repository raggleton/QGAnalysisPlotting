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
    'prescale': 0,
    'lumi': 264313.95,
    'color': ROOT.kRed
}
trig_info['HLT_PFJet60'] = {
    'threshold': 60.,
    'prescale': 0,
    'lumi': 718756.331,
    'color': ROOT.kBlue,
    'fit_all': (lambda isFatJet: not isFatJet)
}
trig_info['HLT_PFJet80'] = {
    'threshold': 80.,
    'prescale': 0,
    'lumi': 2727537.612,
    'color': ROOT.kGreen+2,
    'fit_all': (lambda isFatJet: not isFatJet)
}
trig_info['HLT_PFJet140'] = {
    'threshold': 140.,
    'prescale': 0,
    'lumi': 23949596.439,
    'color': ROOT.kViolet+5,
    'fit_all': (lambda isFatJet: not isFatJet)
}
trig_info['HLT_PFJet200'] = {
    'threshold': 200.,
    'prescale': 0,
    'lumi': 102671470.596,
    'color': ROOT.kOrange,
    'fit_all': (lambda isFatJet: not isFatJet)
}
trig_info['HLT_PFJet260'] = {
    'threshold': 260.,
    'prescale': 0,
    'lumi': 587533374.608,
    'color': ROOT.kTeal,
    'fit_all': (lambda isFatJet: not isFatJet)
}
trig_info['HLT_PFJet320'] = {
    'threshold': 320.,
    'prescale': 0,
    'lumi': 1753875063.756,
    'color': ROOT.kViolet,
    'fit_all': (lambda isFatJet: not isFatJet)
}
trig_info['HLT_PFJet400'] = {
    'threshold': 400.,
    'prescale': 0,
    'lumi': 5139581289.335,
    'color': ROOT.kOrange-6,
    'fit_all': (lambda isFatJet: not isFatJet)
}
trig_info['HLT_PFJet450'] = {
    'threshold': 450.,
    'prescale': 0,
    'lumi': 35918219492.947,
    'color': ROOT.kAzure+1,
    'fit_all': (lambda isFatJet: False) 
}
trig_info['HLT_PFJet500'] = {
    'threshold': 500.,
    'prescale': 0,
    'lumi': 35918219492.947,
    'color': ROOT.kSpring-8,
    'fit_all': (lambda isFatJet: False)
}

total_lumi = trig_info['HLT_PFJet500']['lumi']
for trig_name in trig_info:
    trig_info[trig_name]['prescale'] = total_lumi / trig_info[trig_name]['lumi']


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


def do_trig_plots_vs_singlemu(input_filename, output_dir, title="", eta_min=-2.4, eta_max=2.4, append=""):
    """Do efficiencies and fits for all triggers"""
    if not os.path.isfile(input_filename):
        raise IOError("No input file", input_filename)

    f = ROOT.TFile(input_filename)
    dir_name = "SingleMuRef"

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

    h_all_pt = h_all.ProjectionX("allPT", eta_min_bin, eta_max_bin, "e")
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
        info['hpt'] = h2d.ProjectionX(name+"PT", eta_min_bin, eta_max_bin, "e")
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
    c.SaveAs(output_dir + "/pt_trig_%s_new.%s" % (append, OUTPUT_FMT))

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

        c.SaveAs(output_dir + "/eff_%s_%s_new.%s" % (name, append, OUTPUT_FMT))

    # make graph of fully efficiency pt vs threshold
    # thresholds = [info['threshold'] for info in this_trig_info.itervalues()]
    # fully_eff_pt = [info['good_eff_pt'] for info in this_trig_info.itervalues()]
    # gr = ROOT.TGraph(len(thresholds), array('d', thresholds), array('d', fully_eff_pt))
    # c = ROOT.TCanvas("cgr", "", 800, 600)
    # c.SetTicks(1, 1)
    # gr.SetTitle(";Trigger threshold [GeV];99% efficiency p_{T} [GeV]")
    # gr.SetMarkerStyle(20)
    # # do a pol1 fit
    # thres_fit = ROOT.TF1("f1", "pol1", thresholds[0], thresholds[-1])
    # thres_fit.SetLineColor(ROOT.kRed)
    # thres_fit.SetLineWidth(1)
    # thres_fit.SetLineStyle(2)
    # status = gr.Fit(thres_fit, "RSE")
    # gr.Draw("ALP")
    # c.Modified()
    # c.Update()
    # stats_box = gr.FindObject("stats")
    # stats_box.SetFillColor(ROOT.kWhite)
    # stats_box.SetBorderSize(0)
    # stats_box.SetFillStyle(0)
    # stats_box.SetX1NDC(0.62)
    # stats_box.SetX2NDC(0.88)
    # stats_box.SetY1NDC(0.25)
    # stats_box.SetY2NDC(0.38)
    # cms_text.Draw()
    # jet_text.Draw()
    # c.SaveAs(output_dir + "/fully_eff_pt_vs_threshold_%s.%s" % (append, OUTPUT_FMT))

    return this_trig_info


def do_trig_plots_vs_prevjet(input_filename, output_dir, title="", eta_min=-2.4, eta_max=2.4, append=""):
    """Do efficiencies and fits for all triggers"""
    if not os.path.isfile(input_filename):
        raise IOError("No input file", input_filename)

    f = ROOT.TFile(input_filename)

    is_fat_jet = "AK8" in input_filename.upper()

    this_trig_info = deepcopy(trig_info)

    # first, get the denominator hists for each trigger, ie use plain fire HLT_PFX
    for denom_name, num_name in zip(list(this_trig_info.keys())[:-1], list(this_trig_info.keys())[1:]):
        # num_dict = this_trig_info[num_name]
        denom_dict = this_trig_info[denom_name]        

        dir_name = denom_name+"_v*Ref"
        h2d_denom = f.Get(dir_name + "/pt_vs_eta_all")
        
        # Figure out bins for eta edges
        eta_min_bin, eta_max_bin = 0, h2d_denom.GetNbinsY()+1
        yax = h2d_denom.GetYaxis()

        eta_min_bin = yax.FindBin(eta_min)
        eta_max_bin = yax.FindBin(eta_max)
        # don't want to include the upper bin if the value is at the low edge
        if yax.GetBinLowEdge(eta_max_bin) == eta_max:
            eta_max_bin -= 1
        
        hpt_denom = h2d_denom.ProjectionX(denom_name, eta_min_bin, eta_max_bin, "e")
        hpt_denom.Sumw2()
        print(denom_dict['threshold'], denom_dict['prescale'] )
        hpt_denom.Scale(denom_dict['prescale'])  
        hpt_denom.Rebin(4)  

        hpt_denom.SetLineColor(denom_dict['color'])
        hpt_denom.SetMarkerColor(denom_dict['color'])
        this_trig_info[denom_name]['hpt'] = hpt_denom

    # now we can do the conditional hists
    for denom_name, num_name in zip(list(this_trig_info.keys())[:-1], list(this_trig_info.keys())[1:]):
        num_dict = this_trig_info[num_name]
        denom_dict = this_trig_info[denom_name]        

        dir_name = denom_name+"_v*Ref"
        print(dir_name)
        h2d_num = f.Get(dir_name + "/pt_vs_eta_%s_v*" % num_name)
        hpt_num = h2d_num.ProjectionX(num_name+"And"+denom_name, eta_min_bin, eta_max_bin, "e")
        hpt_num.Sumw2()
        hpt_num.SetLineColor(num_dict['color'])
        hpt_num.SetMarkerColor(num_dict['color'])

        print('Scaling', hpt_num.GetName(), "by", denom_dict['prescale'], "*", num_dict['prescale'])
        hpt_num.Scale(denom_dict['prescale']*num_dict['prescale'])
        hpt_num.Rebin(4)

        num_dict['hpt_cond'] = hpt_num

        num_dict['heff'] = hpt_num.Clone(hpt_num.GetName() + "Eff")
        num_dict['heff'].Sumw2()
        # num_dict['grEff'] = ROOT.TGraphAsymmErrors(this_hpt_num_rebin, this_hpt_denom_rebin)
        
        this_hpt_denom_rebin = denom_dict['hpt']

        print("Creating efficiency from", num_dict['heff'].GetName(), "and", this_hpt_denom_rebin.GetName())
        num_dict['heff'].Divide(num_dict['heff'], this_hpt_denom_rebin, 1, 1, "B")
        new_title = num_dict['heff'].GetTitle().replace("_v*", "")+" relative to "+denom_name.replace("_v*", "")+" ["+append+"]"
        num_dict['heff'].SetTitle(new_title+";Leading jet p_{T} [GeV];#epsilon")


    # for denom_name, num_name in zip(list(this_trig_info.keys())[:-1], list(this_trig_info.keys())[1:]):
    #     denom_dict = this_trig_info[denom_name]
    #     num_dict = this_trig_info[num_name]
    #     print("Doing", num_name, "against", denom_name)

    #     # for each trig, jet 2d pt vs eta hist, project into 1D pt distribution
    #     # then create efficiency hist using previous jet hist hist
    #     rebin_factor = 2 if num_dict['threshold'] > 100 else 2  # rough rebinning across all pt
    #     rebin_factor = 1

    #     # denominator
    #     # ------------------------
    #     dir_name = denom_name+"_v*Ref"
    #     h2d_denom = f.Get(dir_name + "/pt_vs_eta_all")
        
    #     # Figure out bins for eta edges
    #     eta_min_bin, eta_max_bin = 0, h2d_denom.GetNbinsY()+1
    #     yax = h2d_denom.GetYaxis()

    #     eta_min_bin = yax.FindBin(eta_min)
    #     eta_max_bin = yax.FindBin(eta_max)
    #     # don't want to include the upper bin if the value is at the low edge
    #     if yax.GetBinLowEdge(eta_max_bin) == eta_max:
    #         eta_max_bin -= 1
        
    #     # last_prescale = float(list(this_trig_info.values())[-1]['prescale'])
        
    #     total_lumi = list(trig_info.values())[-1]['lumi']

    #     hpt_denom = h2d_denom.ProjectionX(denom_name, eta_min_bin, eta_max_bin, "e")
    #     hpt_denom.Sumw2()
    #     # scale up to account for prescale, assumes last trigger is unprescaled
    #     hpt_denom.Scale(total_lumi / denom_dict['lumi'])  
    #     hpt_denom.SetLineColor(num_dict['color'])
    #     hpt_denom.SetMarkerColor(num_dict['color'])
    #     denom_dict['hpt'] = hpt_denom

    #     this_hpt_denom = hpt_denom.Clone(hpt_denom.GetName()+"Clone").Rebin(rebin_factor)

    #     # rebin the jet pt hist at higher pt where it plateaus
    #     # higher_pt_rebin_factor = 10
    #     # higher_pt_rebin_limit = num_dict['threshold'] * 1.4
    #     # this_hpt_denom_rebin = do_custom_rebin(this_hpt_denom, this_hpt_denom.GetName()+"CustomRebin", higher_pt_rebin_limit, higher_pt_rebin_factor)
    #     this_hpt_denom_rebin = this_hpt_denom

    #     # numerator (next trigger)
    #     # ------------------------
    #     h2d_num = f.Get(dir_name + "/pt_vs_eta_%s_v*" % num_name)
    #     hpt_num = h2d_num.ProjectionX(num_name+"And"+denom_name, eta_min_bin, eta_max_bin, "e")
    #     hpt_num.Sumw2()
    #     # hpt_num.Scale((num_dict['prescale']/last_prescale)/denom_dict['prescale'])  # normalise it
    #     # hpt_num.Scale(num_dict['prescale']/denom_dict['prescale'])  # normalise it
    #     # hpt_num.Scale(total_lumi / num_dict['lumi'])  # normalise it
    #     hpt_num.Scale((total_lumi / denom_dict['lumi']))  # normalise it
    #     hpt_num.SetLineColor(num_dict['color'])
    #     hpt_num.SetMarkerColor(num_dict['color'])
    #     num_dict['hpt_cond'] = hpt_num

    #     this_hpt_num = hpt_num.Clone(hpt_num.GetName()+"Clone").Rebin(rebin_factor)

    #     # rebin the jet pt hist at higher pt where it plateaus
    #     # this_hpt_num_rebin = do_custom_rebin(this_hpt_num, this_hpt_num.GetName()+"CustomRebin", higher_pt_rebin_limit, higher_pt_rebin_factor)
    #     this_hpt_num_rebin = this_hpt_num

    #     num_dict['heff'] = this_hpt_num_rebin.Clone(this_hpt_num_rebin.GetName() + "Eff")
    #     num_dict['heff'].Sumw2()
    #     # num_dict['grEff'] = ROOT.TGraphAsymmErrors(this_hpt_num_rebin, this_hpt_denom_rebin)
        
    #     print("Creating efficiency from", num_dict['heff'].GetName(), "and", this_hpt_denom_rebin.GetName())
    #     num_dict['heff'].Divide(num_dict['heff'], this_hpt_denom_rebin, 1, 1, "B")
    #     new_title = num_dict['heff'].GetTitle().replace("_v*", "")+" relative to "+denom_name.replace("_v*", "")+" ["+append+"]"
    #     num_dict['heff'].SetTitle(new_title+";Leading jet p_{T} [GeV];#epsilon")
        
        # new_title = num_dict['grEff'].GetTitle().replace("_v*", "")+" "+append
        # num_dict['grEff'].SetTitle(new_title+";Leading jet p_{T} [GeV];#epsilon")

    # return
    
    # plot pt distribution of normal triggers
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

    rebin_factor = 1
    for name, info in this_trig_info.items():
        print(info)
        if 'hpt' in info:
            hst.Add(info['hpt'].Rebin(rebin_factor))
            leg.AddEntry(info['hpt'], name, "L")

    c = ROOT.TCanvas("c1", "", 800, 600)
    c.SetTicks(1, 1)
    hst.Draw("HISTE NOSTACK")
    # hst.GetXaxis().SetRangeUser(0, 600)
    hst.SetMinimum(10**2)
    leg.Draw()
    cms_text.Draw()
    jet_text.Draw()
    c.SetLogy()
    c.SaveAs(output_dir + "/pt_trig_%s_prevJet.%s" % (append, OUTPUT_FMT))

    # plot pt distribution of conditional triggers
    c_cond = ROOT.TCanvas("c1_cond", "", 800, 600)
    hst_cond = ROOT.THStack("hst_cond", append+";Jet p_{T} [GeV];N")
    leg = ROOT.TLegend(0.5, 0.5, 0.88, 0.88)
    for name, info in this_trig_info.items():
        # print(info)
        if 'hpt_cond' in info:
            # info['hpt_cond'].SetLineStyle(2)
            hst_cond.Add(info['hpt_cond'].Rebin(rebin_factor))
            leg.AddEntry(info['hpt_cond'], info['hpt_cond'].GetName(), "L")
        
        # if 'hpt' in info:
            # hst_cond.Add(info['hpt']) #.Rebin(rebin_factor))
            # leg.AddEntry(info['hpt'], info['hpt'].GetName(), "L")

    c_cond.SetTicks(1, 1)
    hst_cond.Draw("HISTE NOSTACK")
    # hst_cond.GetXaxis().SetRangeUser(0, 600)
    hst_cond.SetMinimum(10**2)
    leg.Draw()
    cms_text.Draw()
    jet_text.Draw()
    c_cond.SetLogy()
    c_cond.SaveAs(output_dir + "/pt_trig_%s_prevJet_cond.%s" % (append, OUTPUT_FMT))


    # return

    # plot effs
    for ind, (name, info) in enumerate(this_trig_info.items()):
        # skip first trigger as no previous one to compare to
        if ind == 0:
            continue

        c = ROOT.TCanvas("ceff"+name, "", 800, 600)
        c.SetTicks(1, 1)
        c.SetLogy()

        info['heff'].SetMarkerStyle(20)
        # info['heff'].SetTitle(name)
        # info['heff'].SetMaximum(1.5)
        # info['heff'].SetMinimum(0)
        info['heff'].GetXaxis().SetRangeUser(0, min(6*info['threshold'], 2000))
        info['heff'].Draw()
        # info['grEff'].Draw()

        # Do fit
        lower_threshold = info['threshold']/2.
        higher_threshold = info['threshold']*3.
        eff_fit = ROOT.TF1("eff_%s" % name, '[3]*([0] + 0.5 * (1-[0]) * (1 + erf((x-[1])/[2])))', lower_threshold, higher_threshold)
        eff_fit.SetParName(0, 'a')
        eff_fit.SetParName(1, 'mu')
        eff_fit.SetParLimits(1, lower_threshold, higher_threshold)  # enforce +ve parameter values
        eff_fit.SetParName(2, 'sigma')
        eff_fit.SetParLimits(2, 5, 500)
        eff_fit.SetParName(3, 'N')
        eff_fit.SetParLimits(3, 0.00001, 100)
        eff_fit.SetLineColor(ROOT.kBlack)
        eff_fit.SetLineWidth(1)
        eff_fit.SetParameter('a', 0)
        eff_fit.SetParameter('mu', info['threshold'])
        eff_fit.SetParameter('sigma', info['threshold']/10)
        eff_fit.SetParameter('N', 1)
        eff_fit.SetNpx(5000)
        fit_result = info['heff'].Fit(eff_fit, 'RSEQ')
        info['heff'].Draw("")

        ROOT.gPad.Modified()
        ROOT.gPad.Update()

        # # Update fit by increasing lower limit to really capture the high efficiency region
        # # if not info.get("fit_all", lambda x: True)(is_fat_jet):
        fit_factor = 0.4 if is_fat_jet else 0.4
        eff_fit.SetRange(eff_fit.GetX(fit_factor*eff_fit.GetParameter("N")), higher_threshold*1.)
        fit_result = info['heff'].Fit(eff_fit, 'RSEMQ')
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

        c.SaveAs(output_dir + "/eff_prevJet_%s_%s.%s" % (name, append, OUTPUT_FMT))

    # make graph of fully efficiency pt vs threshold
    # thresholds = [info['threshold'] for info in this_trig_info.itervalues()]
    # fully_eff_pt = [info['good_eff_pt'] for info in this_trig_info.itervalues()]
    # gr = ROOT.TGraph(len(thresholds), array('d', thresholds), array('d', fully_eff_pt))
    # c = ROOT.TCanvas("cgr", "", 800, 600)
    # c.SetTicks(1, 1)
    # gr.SetTitle(";Trigger threshold [GeV];99% efficiency p_{T} [GeV]")
    # gr.SetMarkerStyle(20)
    # # do a pol1 fit
    # thres_fit = ROOT.TF1("f1", "pol1", thresholds[0], thresholds[-1])
    # thres_fit.SetLineColor(ROOT.kRed)
    # thres_fit.SetLineWidth(1)
    # thres_fit.SetLineStyle(2)
    # status = gr.Fit(thres_fit, "RSE")
    # gr.Draw("ALP")
    # c.Modified()
    # c.Update()
    # stats_box = gr.FindObject("stats")
    # stats_box.SetFillColor(ROOT.kWhite)
    # stats_box.SetBorderSize(0)
    # stats_box.SetFillStyle(0)
    # stats_box.SetX1NDC(0.62)
    # stats_box.SetX2NDC(0.88)
    # stats_box.SetY1NDC(0.25)
    # stats_box.SetY2NDC(0.38)
    # cms_text.Draw()
    # jet_text.Draw()
    # c.SaveAs(output_dir + "/fully_eff_prevJet_pt_vs_threshold_%s.%s" % (append, OUTPUT_FMT))

    return this_trig_info


def do_plots_and_comparisons(inputs, vs="SingleMu"):
    regions = (
        [-4.7, 4.7, "all"],
        # [-2.4, 2.4, "center"],
        # [-2.4, -1.6, "endcapMinus"],
        # [1.6, 2.4, "endcapPlus"],
        # [-1.6, 1.6, "barrel"]
    )
    all_results = OrderedDict()
    for eta_min, eta_max, append in regions:
        for filename, title in do_these:
            if vs == "SingleMu": 
                results = do_trig_plots_vs_singlemu(filename, os.path.dirname(filename), title, eta_min, eta_max, append)
                all_results[title] = results
            elif vs == "PrevJet":
                results = do_trig_plots_vs_prevjet(filename, os.path.dirname(filename), title, eta_min, eta_max, append)
                all_results[title] = results

        # mg = ROOT.TMultiGraph()
        # c = ROOT.TCanvas("cmg"+append, "", 800, 600)
        # c.SetTicks(1, 1)
        # leg = ROOT.TLegend(0.7, 0.2, 0.88, 0.38)
        # colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, ROOT.kOrange+1]
        # for ind, (name, result) in enumerate(all_results.iteritems()):
        #     thresholds = [info['threshold'] for info in result.itervalues()]
        #     fully_eff_pt = [info['good_eff_pt'] for info in result.itervalues()]
        #     g = ROOT.TGraph(len(thresholds), array('d', thresholds), array('d', fully_eff_pt))
        #     g.SetMarkerColor(colors[ind])
        #     g.SetMarkerStyle(20+ind)
        #     g.SetLineColor(colors[ind])
        #     mg.Add(g)
        #     leg.AddEntry(g, name, "LP")
        # mg.SetTitle(";Trigger threshold [GeV];99% efficiency p_{T} [GeV]")
        # mg.Draw("ALP")
        # leg.Draw()
        # cms_text = ROOT.TPaveText(0.14, 0.9, 0.4, 0.92, "NDC")
        # cms_text.AddText("CMS Preliminary 35.864 fb^{-1}")
        # cms_text.SetFillStyle(0)
        # cms_text.SetBorderSize(0)
        # cms_text.SetTextAlign(ROOT.kHAlignLeft + ROOT.kVAlignBottom)
        # cms_text.SetTextFont(63)
        # cms_text.SetTextSize(18)
        # cms_text.Draw()
        # c.SaveAs("comparingTriggers_%s.%s" % (append, OUTPUT_FMT))


if __name__ == "__main__":
    do_these = [
        ('workdir_ak4chs_jettrig/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu_JetTrig.root', 'AK4 CHS'),
        ('workdir_ak4puppi_jettrig/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu_JetTrig.root', 'AK4 PUPPI'),
        ('workdir_ak8chs_jettrig/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu_JetTrig.root', 'AK8 CHS'),
        ('workdir_ak8puppi_jettrig/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu_JetTrig.root', 'AK8 PUPPI'),
    ]
    # do_plots_and_comparisons(do_these, vs="SingleMu")
    do_these = [
        ('workdir_ak4chs_jettrig/uhh2.AnalysisModuleRunner.DATA.Data_JetHT.root', 'AK4 CHS'),
        ('workdir_ak4puppi_jettrig/uhh2.AnalysisModuleRunner.DATA.Data_JetHT.root', 'AK4 PUPPI'),
        # ('workdir_ak8chs_jettrig/uhh2.AnalysisModuleRunner.DATA.Data_JetHT.root', 'AK8 CHS'),
        ('workdir_ak8puppi_jettrig/uhh2.AnalysisModuleRunner.DATA.Data_JetHT.root', 'AK8 PUPPI'),
    ]
    do_plots_and_comparisons(do_these, vs="PrevJet")
