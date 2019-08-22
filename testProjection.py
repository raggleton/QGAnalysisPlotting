#!/usr/bin/env python

import qg_general_plots as qgp

import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)

h = ROOT.TH2D("h", "", 10, 0, 10, 10, 0, 10)
for i in range(-1, 11):
    val = i+0.1
    print("Fill with value", val)
    h.Fill(val, val, 1)

print("original contents")
for i in range(0, h.GetNbinsX()+4):
    print(i, h.GetXaxis().GetBinLowEdge(i), h.GetBinContent(i, i))

pmin, pmax = 2, 2.3
hproj = qgp.get_projection_plot(h, pmin, pmax, 'y')
print("Projection for values", pmin, pmax)
for i in range(0, hproj.GetNbinsX()+2):
    print(i, hproj.GetBinLowEdge(i), hproj.GetBinContent(i))

pmin, pmax = 2, 3
hproj = qgp.get_projection_plot(h, pmin, pmax, 'y')
print("Projection for values", pmin, pmax)
for i in range(0, hproj.GetNbinsX()+2):
    print(i, hproj.GetBinLowEdge(i), hproj.GetBinContent(i))


pmin, pmax = 2, 4
hproj2 = qgp.get_projection_plot(h, pmin, pmax, 'y')
print("Projection for values", pmin, pmax)
for i in range(0, hproj2.GetNbinsX()+2):
    print(i, hproj2.GetBinLowEdge(i), hproj2.GetBinContent(i))


# test overflow
pmin, pmax = 10, 20
hprojof = qgp.get_projection_plot(h, pmin, pmax, 'y')
print("Projection for values", pmin, pmax)
for i in range(0, hprojof.GetNbinsX()+2):
    print(i, hprojof.GetBinLowEdge(i), hprojof.GetBinContent(i))

# test underflow
pmin, pmax = -2, -1
hprojof = qgp.get_projection_plot(h, pmin, pmax, 'y')
print("Projection for values", pmin, pmax)
for i in range(0, hprojof.GetNbinsX()+2):
    print(i, hprojof.GetBinLowEdge(i), hprojof.GetBinContent(i))


canv = ROOT.TCanvas("c", "", 800, 600)
h.Draw("COLZ")
canv.SaveAs("h2d.pdf")
canv.Clear()

hproj.Draw()
canv.SaveAs("hproj.pdf")

# hprojof.Draw()
# canv.SaveAs("hprojof.pdf")
