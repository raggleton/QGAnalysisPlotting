#!/usr/bin/env python

import ROOT

fileinfos = [
('uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT50to100.root', '50-100', ROOT.kRed),
('uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT100to200.root', '100-200', ROOT.kBlue),
('uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT200to300.root', '200-300', ROOT.kGreen),
('uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT300to500.root', '300-500', ROOT.kOrange),
('uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT500to700.root', '500-700', ROOT.kMagenta),
('uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT700to1000.root', '700-1000', ROOT.kAzure-2),
('uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT1000to1500.root', '1000-1500', ROOT.kCyan),
('uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT1500to2000.root', '1500-2000', ROOT.kGreen+3),
('uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_QCD_HT2000toInf.root', '2000-Inf', ROOT.kYellow-2),
]

c = ROOT.TCanvas ("c", "", 800, 600)
c.SetLogy(1)
c.SetTicks(1, 1)
hst = ROOT.THStack("hst", ";p_{T}^{jet 1} / gen H_{T}; p.d.f.")
leg = ROOT.TLegend(0.6, 0.6, 0.88, 0.88)

files = []

for filename, label, colour in fileinfos:
    f = ROOT.TFile(filename)
    files.append(f)
    hist = f.Get("Dijet_Presel/pt_jet_genHT_ratio")
    hist.SetLineColor(colour)
    hist.SetMarkerColor(colour)
    hist.Scale(1./hist.Integral())
    hst.Add(hist)
    leg.AddEntry(hist, label, "LP")

hst.Draw("HISTE NOSTACK")
hst.SetMaximum(2)
hst.SetMinimum(1E-10)
hst.GetXaxis().SetLimits(0, 15)
leg.Draw()
c.SaveAs("compareHT.pdf")


