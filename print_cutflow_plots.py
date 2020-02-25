#!/usr/bin/env python


"""
Print cutflow plots
"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse

import ROOT
from MyStyle import My_Style
from comparator import Contribution, Plot
My_Style.cd()

# my packages
import common_utils as cu
import qg_common as qgc
import qg_general_plots as qgp


# ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()


def plot_cutflow_hist(hist, output_filename, title='', has_data=False):
    """Plot one cutflow histogram. Normalises so bins are fractions of the first"""
    frac_hist = hist.Clone()
    first_bin = hist.GetBinContent(1)
    for i in range(1, frac_hist.GetNbinsX()+1):
        frac_hist.SetBinContent(i, hist.GetBinContent(i) / first_bin)
        frac_hist.SetBinError(i, hist.GetBinError(i) / first_bin)

    col = ROOT.kBlue
    entry = [Contribution(frac_hist, label='',
                          line_color=col, line_width=2, line_style=1,
                          marker_color=col, marker_size=0.75, marker_style=cu.Marker.get('circle'),
                          normalise_hist=False)]

    hmax = frac_hist.GetMaximum()
    hmin = frac_hist.GetMinimum(0)
    hdiff = hmax-hmin
    ymin = max(0, hmin - (hdiff*0.1))
    ymax = hmax + (hdiff*0.25)

    plot = Plot(entry, 'hist',
                xtitle='',
                ytitle='Fraction',
                title=title,
                ylim=(ymin, ymax),
                legend=False,
                has_data=has_data)
    plot.default_canvas_size = (800, 600)
    plot.plot("NOSTACK HISTE")
    plot.save(output_filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",
                        help="Source file")
    parser.add_argument("--outputDir",
                        help="Output dir (default is source dir)")
    args = parser.parse_args()

    if not os.path.isfile(args.source):
        raise IOError("source file doesn't exist")

    in_tfile = cu.TFileCacher(args.source)
    is_data = '.DATA.' in os.path.basename(args.source)
    is_dijet = '_QCD.root' in os.path.basename(args.source) or '_JetHT_ZeroBias.root' in os.path.basename(args.source)

    if not args.outputDir:
        in_dir = os.path.dirname(args.source)
        args.outputDir = os.path.join(in_dir, 'cutflow_plots')
        if is_data:
            args.outputDir += "_data"
        else:
            args.outputDir += "_MC"
        if is_dijet:
            args.outputDir += "_dijet"
        else:
            args.outputDir += "_zpj"

    print("Plots produced in", args.outputDir)

    jet_str = "AK4 PUPPI jets"
    if 'ak8puppi' in args.source:
        jet_str = "AK8 PUPPI jets"

    # DIJET
    if is_dijet:
        if not is_data:
            # Plot dijet reco cutflow, having passed gen selection (i.e. why inefficient)
            plot_cutflow_hist(in_tfile.Get("cf_DijetSelTighterPassGenCutFlow"),
                              output_filename=os.path.join(args.outputDir, 'cutflow_dijet_reco_sel_tighter_passGen.pdf'),
                              title='Dijet reco selection cutflow, having passed gen selection\n%s' % (jet_str),
                              has_data=False)

            # Plot reco cutflow (no gen requirements)
            plot_cutflow_hist(in_tfile.Get("cf_DijetSelTighterCutFlow"),
                              output_filename=os.path.join(args.outputDir, 'cutflow_dijet_reco_sel_tighter.pdf'),
                              title='Dijet reco selection cutflow (no requirements on gen selection)\n%s' % (jet_str),
                              has_data=False)

            # Plot dijet gen cutflow, having passed reco selection (i.e. why fakes)
            plot_cutflow_hist(in_tfile.Get("cf_DijetGenSelPassRecoCutFlow"),
                              output_filename=os.path.join(args.outputDir, 'cutflow_dijet_gen_sel_passReco.pdf'),
                              title='Dijet gen selection cutflow, having passed reco selection\n%s' % (jet_str),
                              has_data=False)
            # Plot dijet gen cutflow (no reco requirements)
            plot_cutflow_hist(in_tfile.Get("cf_DijetGenSelCutFlow"),
                              output_filename=os.path.join(args.outputDir, 'cutflow_dijet_gen_sel.pdf'),
                              title='Dijet gen selection cutflow (no requirements on reco selection)\n%s' % (jet_str),
                              has_data=False)

            # Plot unfolding fakes cause
            for region in ['central', 'forward']:
                plot_cutflow_hist(in_tfile.Get("Dijet_QG_Unfold_%s_tighter/fakes_counter" % (region)),
                                  output_filename=os.path.join(args.outputDir, 'cutflow_dijet_unfold_%s_fakes.pdf' % (region)),
                                  title='Dijet (%s) unfolding fakes (ungroomed, charged+neutral vars)\n%s' % (region, jet_str),
                                  has_data=False)
                plot_cutflow_hist(in_tfile.Get("Dijet_QG_Unfold_%s_tighter/fakes_countercharged_" % (region)),
                                  output_filename=os.path.join(args.outputDir, 'cutflow_dijet_unfold_%s_charged_fakes.pdf' % (region)),
                                  title='Dijet (%s) unfolding fakes (ungroomed, charged-only vars)\n%s' % (region, jet_str),
                                  has_data=False)
                plot_cutflow_hist(in_tfile.Get("Dijet_QG_Unfold_%s_tighter_groomed/fakes_counter" % (region)),
                                  output_filename=os.path.join(args.outputDir, 'cutflow_dijet_unfold_%s_groomed_fakes.pdf' % (region)),
                                  title='Dijet (%s) unfolding fakes (groomed, charged+neutral vars)\n%s' % (region, jet_str),
                                  has_data=False)
                plot_cutflow_hist(in_tfile.Get("Dijet_QG_Unfold_%s_tighter_groomed/fakes_countercharged_" % (region)),
                                  output_filename=os.path.join(args.outputDir, 'cutflow_dijet_unfold_%s_groomed_charged_fakes.pdf' % (region)),
                                  title='Dijet (%s) unfolding fakes (groomed, charged-only vars)\n%s' % (region, jet_str),
                                  has_data=False)

        else:
            # Plot reco cutflow (no gen requirements)
            plot_cutflow_hist(in_tfile.Get("cf_DijetSelTighterCutFlow"),
                              output_filename=os.path.join(args.outputDir, 'cutflow_dijet_reco_sel_tighter.pdf'),
                              title='Dijet reco selection cutflow\n%s' % (jet_str),
                              has_data=True)
    # Z+Jet
    else:
        if not is_data:
            # Plot dijet reco cutflow, having passed gen selection (i.e. why inefficient)
            plot_cutflow_hist(in_tfile.Get("cf_ZPlusJetsSelPassGenCutFlow"),
                              output_filename=os.path.join(args.outputDir, 'cutflow_ZPlusJets_reco_sel_passGen.pdf'),
                              title='Z+jet reco selection cutflow, having passed gen selection\n%s' % (jet_str),
                              has_data=False)

            # Plot reco cutflow (no gen requirements)
            plot_cutflow_hist(in_tfile.Get("cf_ZPlusJetsSelCutFlow"),
                              output_filename=os.path.join(args.outputDir, 'cutflow_ZPlusJets_reco_sel.pdf'),
                              title='Z+jet reco selection cutflow (no requirements on gen selection)\n%s' % (jet_str),
                              has_data=False)

            # Plot dijet gen cutflow, having passed reco selection (i.e. why fakes)
            plot_cutflow_hist(in_tfile.Get("cf_ZPlusJetsGenSelPassRecoCutFlow"),
                              output_filename=os.path.join(args.outputDir, 'cutflow_ZPlusJets_gen_sel_passReco.pdf'),
                              title='Z+jet gen selection cutflow, having passed reco selection\n%s' % (jet_str),
                              has_data=False)
            # Plot dijet gen cutflow (no reco requirements)
            plot_cutflow_hist(in_tfile.Get("cf_ZPlusJetsGenSelCutFlow"),
                              output_filename=os.path.join(args.outputDir, 'cutflow_ZPlusJets_gen_sel.pdf'),
                              title='Z+jet gen selection cutflow (no requirements on reco selection)\n%s' % (jet_str),
                              has_data=False)

            # Plot unfolding fakes cause
            plot_cutflow_hist(in_tfile.Get("ZPlusJets_QG_Unfold_tighter/fakes_counter"),
                              output_filename=os.path.join(args.outputDir, 'cutflow_ZPlusJets_unfold_fakes.pdf'),
                              title='Z+jet unfolding fakes (ungroomed, charged+neutral vars)\n%s' % (jet_str),
                              has_data=False)
            plot_cutflow_hist(in_tfile.Get("ZPlusJets_QG_Unfold_tighter/fakes_countercharged_"),
                              output_filename=os.path.join(args.outputDir, 'cutflow_ZPlusJets_unfold_charged_fakes.pdf'),
                              title='Z+jet unfolding fakes (ungroomed, charged-only vars)\n%s' % (jet_str),
                              has_data=False)
            plot_cutflow_hist(in_tfile.Get("ZPlusJets_QG_Unfold_tighter_groomed/fakes_counter"),
                              output_filename=os.path.join(args.outputDir, 'cutflow_ZPlusJets_unfold_groomed_fakes.pdf'),
                              title='Z+jet unfolding fakes (groomed, charged+neutral vars)\n%s' % (jet_str),
                              has_data=False)
            plot_cutflow_hist(in_tfile.Get("ZPlusJets_QG_Unfold_tighter_groomed/fakes_countercharged_"),
                              output_filename=os.path.join(args.outputDir, 'cutflow_ZPlusJets_unfold_groomed_charged_fakes.pdf'),
                              title='Z+jet unfolding fakes (groomed, charged-only vars)\n%s' % (jet_str),
                              has_data=False)

        else:
            # Plot reco cutflow (no gen requirements)
            plot_cutflow_hist(in_tfile.Get("cf_ZPlusJetsSelCutFlow"),
                              output_filename=os.path.join(args.outputDir, 'cutflow_ZPlusJets_reco_sel.pdf'),
                              title='Z+jet reco selection cutflow\n%s' % (jet_str),
                              has_data=True)
