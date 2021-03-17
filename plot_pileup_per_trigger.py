#!/usr/bin/env python

"""Plot pileup profile per trigger path

Create them with doPileupProfilesPerTrigger.sh
"""

import os
os.nice(10)
import numpy as np
import ROOT

import common_utils as cu
from MyStyle import My_Style
My_Style.cd()
from comparator import Plot, Contribution


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
# ROOT.gStyle.SetPalette(55)


def do_pileup_plot(input_dir, trigger_names, output_filename):
    # get histograms
    hists = [cu.grab_obj_from_file(os.path.join(input_dir, 'MyDataPileupHistogram_%s.root' % t), 'pileup') 
             for t in trigger_names]
    h_up = cu.grab_obj_from_file(os.path.join(input_dir, 'MyDataPileupHistogram_PFJet500_72383.root'), 'pileup')
    h_down = cu.grab_obj_from_file(os.path.join(input_dir, 'MyDataPileupHistogram_PFJet500_66017.root'), 'pileup')
    h_down3 = cu.grab_obj_from_file(os.path.join(input_dir, 'MyDataPileupHistogram_PFJet500_59650.root'), 'pileup')
    ratio_up = h_up.Clone()
    ratio_up.Divide(hists[-1])
    ratio_down = h_down.Clone()
    ratio_down.Divide(hists[-1])
    # hists.append(h_up)
    # hists.append(h_down)
    # trigger_names.append('72.383')
    # trigger_names.append('66.017')
    # Create contributions
    mark = cu.Marker()
    n_hists = len(hists)
    conts = [Contribution(h, label=t, 
                          line_width=2,
                          line_color=cu.get_colour_seq(ind, n_hists),
                          marker_color=cu.get_colour_seq(ind, n_hists),
                          marker_style=m,
                          normalise_hist=True)
             for ind, (h, t, m) in enumerate(zip(hists, trigger_names, mark.cycle(cycle_filling=True)))]
    

    conts.insert(-1, 
        Contribution(h_up, label='72.383 (+4.6%)',
                     line_width=2,
                     line_color=ROOT.kRed,
                     marker_color=ROOT.kRed,
                     normalise_hist=True,
                     )
    )
    conts.insert(-1, 
        Contribution(h_down, label='66.017 (-4.6%)',
                     line_width=2,
                     line_color=ROOT.kMagenta,
                     marker_color=ROOT.kMagenta,
                     normalise_hist=True,
                     )
    )
    conts.insert(-1, 
        Contribution(h_down3, label='59.650 (-13.8%)',
                     line_width=2,
                     line_style=1,
                     line_color=ROOT.kMagenta+3,
                     marker_color=ROOT.kMagenta+3,
                     normalise_hist=True,
                     )
    )
    print([h.Integral() for h in hists])
    plot = Plot(conts, what='hist',
                xtitle='Pileup',
                lumi=cu.get_lumi_str(do_dijet=True, do_zpj=True),
                subplot_type='ratio',
                subplot=conts[-1],
                )
    plot.subplot_maximum_ceil = 4.5
    plot.plot("NOSTACK HISTE", "NOSTACK HIST")
    plot.subplot_pad.cd()
    
    # To shade region between up and down hists, need to create graph with
    # the error bars as the up/down variations
    
    x_low = np.array([ratio_up.GetBinLowEdge(i) for i in range(1, ratio_up.GetNbinsX()+1)], dtype=float)
    x_high = np.array([ratio_up.GetBinLowEdge(i) for i in range(2, ratio_up.GetNbinsX()+2)], dtype=float)
    x_mid = 0.5*(x_low + x_high)
    x_err = x_high - x_low
    
    y_high = np.array([max(ratio_up.GetBinContent(i), ratio_down.GetBinContent(i)) for i in range(1, ratio_up.GetNbinsX()+1)])
    y_high -= 1
    y_low = np.array([min(ratio_up.GetBinContent(i), ratio_down.GetBinContent(i)) for i in range(1, ratio_up.GetNbinsX()+1)])
    y_low = 1 - y_low
    y_mid = np.array([1. for i in range(1, ratio_up.GetNbinsX()+1)])
    
    gr = ROOT.TGraphAsymmErrors(len(x_mid), x_mid, y_mid, x_err, x_err, y_low, y_high)
    gr.SetFillColor(ROOT.kGray+2)
    gr.SetFillStyle(3254)
    gr.Draw('2 SAME')
    plot.subplot_container.Draw("SAME NOSTACK HIST")
    plot.subplot_line.Draw()
    plot.save(output_filename)


if  __name__ == "__main__":
    ak4_triggers = [
        "ZeroBias",
        "PFJet40",
        "PFJet60",
        "PFJet80",
        "PFJet140",
        "PFJet200",
        "PFJet260",
        "PFJet320",
        "PFJet400",
        "PFJet450",
        # "PFJet500",
    ]

    ak8_triggers = [
        "ZeroBias",
        "AK8PFJet40",
        "AK8PFJet60",
        "AK8PFJet80",
        "AK8PFJet140",
        "AK8PFJet200",
        "AK8PFJet260",
        "AK8PFJet320",
        "AK8PFJet400",
        "AK8PFJet450",
    ]

    input_dir = "/Volumes/Extreme SSD/Projects/QGAnalysis/jet_trig_pu/"
    # do ak4 plot
    do_pileup_plot(input_dir, ak4_triggers, 'pileup_per_trigger_ak4.pdf')

    # do ak8 plot
    do_pileup_plot(input_dir, ak8_triggers, 'pileup_per_trigger_ak8.pdf')
