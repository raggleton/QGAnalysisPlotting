#!/usr/bin/env python

"""Count # events passing trigger above certain threshold"""

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
import common_utils as cu

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit()
# ROOT.gStyle.SetStatColor()


trig_info = OrderedDict()

# AK8 limits with normal PFJet trigs
normal_trig_info = OrderedDict()

normal_trig_info['HLT_PFJet40'] = {
    'nominal_threshold': 40.,
    'actual_threshold': 85.,
    'lumi': 264313.95,
    'color': ROOT.kRed
}
normal_trig_info['HLT_PFJet60'] = {
    'nominal_threshold': 60.,
    'actual_threshold': 120.,
    'lumi': 718756.331,
    'color': ROOT.kBlue,
}
normal_trig_info['HLT_PFJet80'] = {
    'nominal_threshold': 80.,
    'actual_threshold': 148.,
    'lumi': 2727537.612,
    'color': ROOT.kGreen+2,
}
normal_trig_info['HLT_PFJet140'] = {
    'nominal_threshold': 140.,
    'actual_threshold': 233.,
    'lumi': 23949596.439,
    'color': ROOT.kViolet+5,
}
normal_trig_info['HLT_PFJet200'] = {
    'nominal_threshold': 200.,
    'actual_threshold': 326.,
    'lumi': 102671470.596,
    'color': ROOT.kOrange,
}
normal_trig_info['HLT_PFJet260'] = {
    'nominal_threshold': 260.,
    'actual_threshold': 405.,
    'lumi': 587533374.608,
    'color': ROOT.kTeal,
}
normal_trig_info['HLT_PFJet320'] = {
    'nominal_threshold': 320.,
    'actual_threshold': 477.,
    'lumi': 1753875063.756,
    'color': ROOT.kViolet,
}
normal_trig_info['HLT_PFJet400'] = {
    'nominal_threshold': 400.,
    'actual_threshold': 574.,
    'lumi': 5139581289.335,
    'color': ROOT.kOrange-6,
}
normal_trig_info['HLT_PFJet450'] = {
    'nominal_threshold': 450.,
    'actual_threshold': 614.,
    'lumi': 35918219492.947,
    'color': ROOT.kAzure+1,
}

# AK8 limits with AK8PFJet trigs
ak8_trig_info = OrderedDict()
ak8_trig_info['HLT_AK8PFJet40'] = {
    'nominal_threshold': 40.,
    'actual_threshold': 70.,
    'lumi': 49176.854,
    'color': ROOT.kRed
}
ak8_trig_info['HLT_AK8PFJet60'] = {
    'nominal_threshold': 60.,
    'actual_threshold': 88.,
    'lumi': 324768.681,
    'color': ROOT.kBlue,
}
ak8_trig_info['HLT_AK8PFJet80'] = {
    'nominal_threshold': 80.,
    'actual_threshold': 129.,
    'lumi': 994563.916,
    'color': ROOT.kGreen+2,
}
ak8_trig_info['HLT_AK8PFJet140'] = {
    'nominal_threshold': 140.,
    'actual_threshold': 185.,
    'lumi': 10005654.493,
    'color': ROOT.kViolet+5,
}
ak8_trig_info['HLT_AK8PFJet200'] = {
    'nominal_threshold': 200.,
    'actual_threshold': 252.,
    'lumi': 84893034.092,
    'color': ROOT.kOrange,
}
ak8_trig_info['HLT_AK8PFJet260'] = {
    'nominal_threshold': 260.,
    'actual_threshold': 320.,
    'lumi': 512841061.578,
    'color': ROOT.kTeal,
}
ak8_trig_info['HLT_AK8PFJet320'] = {
    'nominal_threshold': 320.,
    'actual_threshold': 390.,
    'lumi': 1510155111.062,
    'color': ROOT.kViolet,
}
ak8_trig_info['HLT_AK8PFJet400'] = {
    'nominal_threshold': 400.,
    'actual_threshold': 479.,
    'lumi': 4544785568.903,
    'color': ROOT.kOrange-6,
}
ak8_trig_info['HLT_AK8PFJet450'] = {
    'nominal_threshold': 450.,
    'actual_threshold': 536.,
    'lumi': 33182262109.421,
    'color': ROOT.kAzure+1,
}

# all lumis in ub-1
TOTAL_LUMI = 35918219492.947
ZB_LUMI = 29048.362


def do_counting(zb_input_file, jetht_input_file, trig_info, eta_min=-2.4, eta_max=2.4, title="", output_filename="pt.pdf"):
    # first get num event passing ZeroBias trigger, up to the first trigger threshold
    h_zb = cu.get_from_tfile(zb_input_file, "ZeroBiasRef/pt_vs_eta_all")

    # Figure out bins for eta edges
    eta_min_bin, eta_max_bin = 0, h_zb.GetNbinsY()+1
    yax = h_zb.GetYaxis()

    eta_min_bin = yax.FindBin(eta_min)
    eta_max_bin = yax.FindBin(eta_max)
    # don't want to include the upper bin if the value is at the low edge
    if yax.GetBinLowEdge(eta_max_bin) == eta_max:
        eta_max_bin -= 1

    # Cut for this region of eta and pt phase space
    pt_min = 50
    trig_names = list(trig_info.keys())
    pt_max = trig_info[trig_names[0]]['actual_threshold']
    cut_name = "zbcut"
    this_pt_cut = ROOT.TCutG(cut_name, 5)
    this_pt_cut.SetPoint(0, pt_min, eta_min)
    this_pt_cut.SetPoint(1, pt_max, eta_min)
    this_pt_cut.SetPoint(2, pt_max, eta_max)
    this_pt_cut.SetPoint(3, pt_min, eta_max)
    this_pt_cut.SetPoint(4, pt_min, eta_min)

    h_zb_pt = h_zb.ProjectionX("ZeroBias", eta_min_bin, eta_max_bin, "["+cut_name+"]")
    h_zb_pt.Sumw2()
    h_zb_pt.SetFillColor(17)
    h_zb_pt.SetLineColor(17)

    print("WARNING: rounding to nearest bin width", h_zb_pt.GetBinWidth(1))

    hst = ROOT.THStack("hst", title+";Leading jet p_{T} [GeV];N")
    hst.Add(h_zb_pt)

    total = []
    this_num_entries = h_zb_pt.Integral()
    total.append(("ZeroBias", this_num_entries))

    leg = ROOT.TLegend(0.6, 0.6, 0.88, 0.88)
    leg.AddEntry(h_zb_pt, "ZeroBias (p_{T}>%g GeV): %g" % (pt_min, this_num_entries), "LF")

    # Now go through each of the triggers and get their hists
    for this_trig_name, next_trig_name in zip(trig_names[:-1], trig_names[1:]):
        h2d = cu.get_from_tfile(jetht_input_file, this_trig_name+"_v*Ref/pt_vs_eta_all")
        pt_min = trig_info[this_trig_name]['actual_threshold']
        pt_max = trig_info[next_trig_name]['actual_threshold']
        cut_name = this_trig_name+"cut"
        this_pt_cut = ROOT.TCutG(cut_name, 5)
        this_pt_cut.SetPoint(0, pt_min, eta_min)
        this_pt_cut.SetPoint(1, pt_max, eta_min)
        this_pt_cut.SetPoint(2, pt_max, eta_max)
        this_pt_cut.SetPoint(3, pt_min, eta_max)
        this_pt_cut.SetPoint(4, pt_min, eta_min)

        h_pt = h2d.ProjectionX(this_trig_name, eta_min_bin, eta_max_bin, "["+cut_name+"]")
        h_pt.Sumw2()
        h_pt.SetFillColor(trig_info[this_trig_name]['color'])
        h_pt.SetLineColor(trig_info[this_trig_name]['color'])

        this_num_entries = h_pt.Integral()
        total.append((this_trig_name, this_num_entries))

        hst.Add(h_pt)
        leg.AddEntry(h_pt, this_trig_name + ": %g" % (this_num_entries), "LF")


    # Do the last trigger manually
    last_trig_name = trig_names[-1]
    h2d = cu.get_from_tfile(jetht_input_file, last_trig_name+"_v*Ref/pt_vs_eta_all")
    pt_min = trig_info[last_trig_name]['actual_threshold']
    pt_max = 10000
    cut_name = last_trig_name+"cut"
    this_pt_cut = ROOT.TCutG(cut_name, 5)
    this_pt_cut.SetPoint(0, pt_min, eta_min)
    this_pt_cut.SetPoint(1, pt_max, eta_min)
    this_pt_cut.SetPoint(2, pt_max, eta_max)
    this_pt_cut.SetPoint(3, pt_min, eta_max)
    this_pt_cut.SetPoint(4, pt_min, eta_min)

    h_pt = h2d.ProjectionX(last_trig_name, eta_min_bin, eta_max_bin, "["+cut_name+"]")
    h_pt.Sumw2()
    h_pt.SetFillColor(trig_info[last_trig_name]['color'])
    h_pt.SetLineColor(trig_info[last_trig_name]['color'])

    this_num_entries = h_pt.Integral()
    total.append((last_trig_name, this_num_entries))

    hst.Add(h_pt)
    leg.AddEntry(h_pt, last_trig_name + ": %g" % (this_num_entries), "LF")

    c = ROOT.TCanvas("c", "", 800, 600)
    c.SetTicks(1, 1)
    # c.SetLogx()
    c.SetLogy()
    hst.Draw("HISTE NOSTACK")
    hst.SetMinimum(10)
    hst.SetMaximum(1000000)
    total_count = sum(x[1] for x in total)
    leg.AddEntry(0, "Total: "+str(total_count), "")
    leg.Draw()
    c.SaveAs(output_filename)

    for count in total:
        print(count[0], "=", count[1])
    print("Total:", total_count)


if __name__ == "__main__":
    
    zb_file = cu.open_root_file('workdir_ak8puppi_jettrig_withAK8trig/uhh2.AnalysisModuleRunner.DATA.Data_ZeroBias.root')
    jetht_file = cu.open_root_file('workdir_ak8puppi_jettrig_withAK8trig/uhh2.AnalysisModuleRunner.DATA.Data_JetHT.root')

    # count for normal triggers
    do_counting(zb_file, jetht_file, normal_trig_info, -2.4, 2.4, "AK8PUPPI + HLT_PFJet*", "workdir_ak8puppi_jettrig_withAK8trig/ak8puppi_normaltrig.pdf")
    
    # count for AK8 specific triggers
    do_counting(zb_file, jetht_file, ak8_trig_info, -2.4, 2.4, "AK8PUPPI + HLT_AK8PFJet*", "workdir_ak8puppi_jettrig_withAK8trig/ak8puppi_ak8trig.pdf")
