#!/usr/bin/env python


"""
How to fit a TH1 using multiple template TH1s

Taken from https://root-forum.cern.ch/t/hist-template-fit/15218/2
"""

from functools import partial
import ROOT

ROOT.gStyle.SetOptStat(1111)
ROOT.gStyle.SetOptFit(1111)

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()


# Create template TH1s
h1 = ROOT.TH1D("h1", "", 100, -10, 10)
h2 = ROOT.TH1D("h2", "", 100, -10, 10)

n = 100000;
for i in range(n): h1.Fill(ROOT.gRandom.Gaus(-3));
for i in range(n): h2.Fill(ROOT.gRandom.Gaus(+3));
h1.Scale(1./n)
h2.Scale(1./n)


# Define fit function to be used in TH1::Fit, that really uses the templates
# note that TH1::Fit expects only the x and pars arguments - we can use
# functools.partial to tell it which hist1 and hist2 to use
def double_gaussian(x, pars, hist1, hist2):
    xx = x[0]; # use x[1] to get 2nd dimension, x[2] for 3rd ...
    # the fit parameters, i.e. the histogram weights
    w1 = pars[0];
    w2 = pars[1];

    # get content of the histograms for this point
    y1 = hist1.GetBinContent(hist1.GetXaxis().FindFixBin(xx));
    y2 = hist2.GetBinContent(hist2.GetXaxis().FindFixBin(xx));

    return w1*y1 + w2*y2


# Create data hist to be fitted
h = ROOT.TH1D("h", "data", 100, -10, 10)
for i in range(100):
    h.Fill(ROOT.gRandom.Gaus(-3))
for i in range(100):
    h.Fill(ROOT.gRandom.Gaus(+3))

# create a fit function from a python function
# note the use of functools.partial
f = ROOT.TF1("double_gaussian", partial(double_gaussian, hist1=h1, hist2=h2), -10, 10, 2)
f.SetNpx(10000)
h.Fit(f, "L")

# draw
canv = ROOT.TCanvas("c", "", 800, 600)
h.Draw()
canv.SaveAs("test_fit.pdf")
