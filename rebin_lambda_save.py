#!/usr/bin/env python

"""Make 1D lambda hists & rebin according to custom rebinning. 
Also set names ready for RIVET/YODA.
"""


import ROOT
import os
import argparse
from array import array

# My stuff
import qg_common as qgc
import qg_general_plots as qgp
import common_utils as cu

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)


def create_yoda_hist_name(region, angle_name, pt_bin):
    """Produce YODA hist name"""
    region_num = 1 if "dijet" in region.lower() else 2
    angle_names = [a.var for a in qgc.COMMON_VARS]
    angle_num = angle_names.index(angle_name)+1
    pt_bin_num = qgc.PT_BINS.index(pt_bin)+1
    return "d{:02d}-x{:02d}-y{:02d}".format(region_num, angle_num, pt_bin_num)

YODA_PLOT_TEMPLATE = """BEGIN PLOT /CMS_2018_PAS_SMP_18_QGX/{plotname}
Title=$\sqrt{{s}} = 13$ TeV, $\mathcal{{L}}_{{\mathrm{{int}}}} = 35.9~\mathrm{{fb}}^{{-1}}$, Dijet selection, {title}
XLabel={xlabel}
YLabel={ylabel}
XMin={xmin}
XMax={xmax}
LogX={logx}
LogY={logy}
LegendXPos=0.65
LegendYPos=0.9
Legend=1
ErrorBands=1
ErrorBars=1
ErrorBandOpacity=0.5
END PLOT

"""
"""
BEGIN HISTOGRAM /CMS_2018_PAS_SMP_18_QGX/{plotname}
ErrorBars=1
ErrorBandOpacity=1
ErrorType=stat
PolyMarker=o
ConnectBins=0
END HISTOGRAM
"""


def make_1D_rebin_hists(input_filename, plot_dirname, output_filename):
    
    # QG variable plots

    in_f = cu.open_root_file(input_filename)
    out_f = ROOT.TFile(output_filename, "RECREATE")

    plot_file = open(output_filename.replace(".root", ".plot"), "w")

    for ang_ind, ang in enumerate(qgc.COMMON_VARS[:5]):
        if ang.var not in qgc.ANGLE_REBIN_DICT:
            continue

        in_f.cd()

        var_prepend = ""
        obj_name = "%s%s_vs_pt" % (var_prepend, ang.var)
        if var_prepend == "gen" and "multiplicity" in ang.var.lower():
            obj_name = obj_name.replace("puppiMultiplicity", "multiplicity")
        h2d = cu.get_from_tfile(in_f, "%s/%s" % (plot_dirname, obj_name))

        for pt_ind, (start_val, end_val) in enumerate(qgc.PT_BINS):
            hist = qgp.get_projection_plot(h2d, start_val, end_val, 'y')  # y cos lambda is on x axis
            this_rebins = qgc.ANGLE_REBIN_DICT[ang.var]
            # new_name = "%s_Pt%sto%d" % (ang.var, start_val, end_val)
            yoda_name = create_yoda_hist_name(plot_dirname, ang.var, (start_val, end_val))
            rebin_hist = hist.Rebin(len(this_rebins)-1, yoda_name, array('d', this_rebins))
            # Normalise
            if rebin_hist.Integral() > 0:
                rebin_hist.Scale(1./rebin_hist.Integral())
            # Divide bin contents by bin width
            for bin_ind in range(1, rebin_hist.GetNbinsX()+1):
                val = rebin_hist.GetBinContent(bin_ind)
                width = rebin_hist.GetBinWidth(bin_ind)
                rebin_hist.SetBinContent(bin_ind, val/width)
            
            out_f.WriteTObject(rebin_hist)  # saves faffing with cd()

            # create yoda plot info
            xlabel = "%s $(%s)$" % (ang.name, ang.lambda_str.replace("#", "\\"))
            title_str = "$%d < p_{T}^{jet} < %d$ GeV" % (start_val, end_val)
            xmax = 1
            if "multiplicity" in ang.var.lower():
                xmax = 50
            this_plot_info = YODA_PLOT_TEMPLATE.format(plotname=yoda_name, xlabel=xlabel, title=title_str, ylabel="p.d.f", logx=0, logy=0, xmax=xmax, xmin=0)
            plot_file.write(this_plot_info)

    # Add in pt hist
    h_pt = cu.get_from_tfile(in_f, "Dijet_tighter/pt_jet_response_binning")
    pt_name = "d01-x20-y01"
    h_pt.SetName(pt_name)
    this_plot_info = YODA_PLOT_TEMPLATE.format(plotname=pt_name, xlabel="$\langle p_{T}^{j}\\rangle$ GeV", title="", ylabel="N", logy=1, logx=1, xmin=10, xmax=1000)
    plot_file.write(this_plot_info)
    out_f.WriteTObject(h_pt)

    out_f.Close()
    plot_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input file", required=True)
    parser.add_argument("--plotDirname", help="Plot directory name in input file", required=True)
    parser.add_argument("--output", help="Output file", required=True)
    args = parser.parse_args()

    make_1D_rebin_hists(args.input, args.plotDirname, args.output)
