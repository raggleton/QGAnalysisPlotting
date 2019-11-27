#!/usr/bin/env python

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
ROOT.gErrorIgnoreLevel = ROOT.kWarning


tfile_dy = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root")
tfile_ww = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet/uhh2.AnalysisModuleRunner.MC.MC_WW.root")
tfile_wz = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet/uhh2.AnalysisModuleRunner.MC.MC_WZ.root")
tfile_zz = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet/uhh2.AnalysisModuleRunner.MC.MC_ZZ.root")
tfile_data = cu.open_root_file("workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet/uhh2.AnalysisModuleRunner.DATA.Data_SingleMu.root")

h_dy = cu.get_from_tfile(tfile_dy, "ZPlusJets/pt_jet1")
h_ww = cu.get_from_tfile(tfile_ww, "ZPlusJets/pt_jet1")
h_wz = cu.get_from_tfile(tfile_wz, "ZPlusJets/pt_jet1")
h_zz = cu.get_from_tfile(tfile_zz, "ZPlusJets/pt_jet1")
h_data = cu.get_from_tfile(tfile_data, "ZPlusJets/pt_jet1")

print(h_dy.Integral())
print(h_ww.Integral())
print(h_wz.Integral())
print(h_zz.Integral())
print(h_data.Integral())
h_dy.SetFillColor(ROOT.kCyan+2)
h_dy.SetLineColor(ROOT.kCyan+2)
h_dy.SetFillStyle(1001)

h_ww.SetFillColor(ROOT.kGreen)
h_ww.SetLineColor(ROOT.kGreen)
h_ww.SetFillStyle(1001)

h_wz.SetFillColor(ROOT.kRed)
h_wz.SetLineColor(ROOT.kRed)
h_wz.SetFillStyle(1001)

h_zz.SetFillColor(ROOT.kBlue)
h_zz.SetLineColor(ROOT.kBlue)
h_zz.SetFillStyle(1001)

h_data.SetMarkerColor(ROOT.kBlack)
h_data.SetLineColor(ROOT.kBlack)
h_data.SetMarkerSize(1)
h_data.SetMarkerStyle(20)

canv = ROOT.TCanvas("c", "", 800, 600)
canv.SetTicks(1, 1)
canv.SetLogy()
canv.SetLogx()

hst = ROOT.THStack("hst", ";p_{T}^{jet} [GeV];N")
hst.Add(h_ww)
hst.Add(h_zz)
hst.Add(h_wz)
hst.Add(h_dy)

leg = ROOT.TLegend(0.6, 0.6, 0.88, 0.88)
leg.AddEntry(h_data, "Data (%g)" % (h_data.Integral()), "P")
leg.AddEntry(h_dy, "DY (%g)" % (h_dy.Integral()), "F")
leg.AddEntry(h_wz, "WZ (%g)" % (h_wz.Integral()), "F")
leg.AddEntry(h_zz, "ZZ (%g)" % (h_zz.Integral()), "F")
leg.AddEntry(h_ww, "WW (%g)" % (h_ww.Integral()), "F")

hst.Draw("hist")
hst.SetMaximum(1E7)
hst.SetMinimum(1E-1)
h_data.Draw("SAME E")
leg.Draw()

canv.SaveAs("zpj_pt.pdf")

