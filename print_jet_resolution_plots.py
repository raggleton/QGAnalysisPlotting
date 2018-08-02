#!/usr/bin/env python

"""Print plots of jet resolution (reco pt - gen pt / gen pt) split by pt, etc"""


import argparse
from MyStyle import My_Style
My_Style.cd()
import os
from itertools import product
from array import array

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)


# My stuff
from comparator import Contribution, Plot, grab_obj
import qg_common as qgc
import qg_general_plots as qgg
import common_utils as cu



# Control output format
OUTPUT_FMT = "pdf"


def determine_fit_range(hist):
    """Determine lower & upper limit of fit range
    
    Parameters
    ----------
    hist : TH1
        Description
    
    Returns
    -------
    tuple
        (lower limit, upper limit)
    """
    # mean = hist.GetMean()
    mean = hist.GetBinCenter(hist.GetMaximumBin())
    rms = hist.GetRMS()
    multiplier = 1.5
    return (mean - multiplier*rms, mean + multiplier*rms)


def do_gaus_fit(hist):
    """Do a Gaussian fit to histogram

    Parameters
    ----------
    hist : TH1
        Histogram to fit to
    """
    func_name = hist.GetName()+"_f1"
    func_name = "gausFit"
    fit_range = determine_fit_range(hist)
    func = ROOT.TF1(func_name, "gaus", fit_range[0], fit_range[1])
    # func.SetParameters(hist.GetMaximum(), hist.GetMean(), hist.GetRMS())
    fit_result = hist.Fit(func_name, "ERSQ", "L")
    # print("fit result:", int(fit_result))


def fit_results_to_str(fit):
    """Turn fit results into str, lines split by \n
    
    Parameters
    ----------
    fit : TF1
        Description
    
    Returns
    -------
    str
        Description
    
    """
    parts = []
    chi2 = fit.GetChisquare()
    ndf = fit.GetNDF()
    if ndf > 0:
        parts.append("chi2/ndof: %.3e/%d = %.3e" % (chi2, ndf, chi2/ndf))
    else:
        parts.append("chi2/ndof: %.3e/0 = Inf" % (chi2))
    parts.append("prob: %.3e" % fit.GetProb())
    for i in range(fit.GetNpar()):
        parts.append("%s: %.3e #pm %.3e" % (fit.GetParName(i), fit.GetParameter(i), fit.GetParError(i)))
    return "\n".join(parts)


def do_projection_plots(in_file, plot_dir, do_fit=True, skip_dirs=None):
    hist_name = "pt_jet_response"
    tfile = cu.open_root_file(in_file)
    dirs = cu.get_list_of_element_names(tfile)


    for mydir in dirs:
        if skip_dirs and mydir in skip_dirs:
            continue

        if hist_name not in cu.get_list_of_element_names(tfile.Get(mydir)):
            continue

        print("Doing", mydir)

        h2d = cu.grab_obj_from_file(in_file, "%s/%s" % (mydir, hist_name))

        ax = h2d.GetXaxis()
        bin_edges = [ax.GetBinLowEdge(i) for i in range(1, ax.GetNbins()+2)]
        
        bin_centers, sigmas, sigmas_unc = [], [], []

        for pt_min, pt_max in zip(bin_edges[:-1], bin_edges[1:]):
            obj = qgg.get_projection_plot(h2d, pt_min, pt_max, cut_axis='x')
            if obj.GetEffectiveEntries() < 20:
                continue
            # obj.Rebin(rebin)
            obj.Scale(1./obj.Integral())

            label = "%s < p_{T}^{Gen} < %s GeV" % (str(pt_min), str(pt_max))
            if do_fit:
                do_gaus_fit(obj)
                fit = obj.GetFunction("gausFit")
                label += "\n"
                label += fit_results_to_str(fit)
                # bin_centers.append(fit.GetParameter(1))
                bin_centers.append(0.5*(pt_max+pt_min))
                sigmas.append(fit.GetParameter(2))
                sigmas_unc.append(fit.GetParError(2))

            # output_filename = os.path.join(plot_dir, "%s_%s_ptGen%sto%s.%s" % (mydir, hist_name, str(pt_min), str(pt_max), OUTPUT_FMT))

            # cont = Contribution(obj, label=label)
            # delta = pt_max - pt_min
            # # xlim = (pt_min - 10*delta, pt_max + 10*delta)
            # xlim = (obj.GetMean()-3*obj.GetRMS(), obj.GetMean()+3*obj.GetRMS())
            # ylim = (0, obj.GetMaximum()*1.1)
            # plot = Plot([cont], what='hist', 
            #             xtitle="p_{T}^{Reco} [GeV]", xlim=xlim, ylim=ylim)
            # plot.plot()  # don't use histe as it wont draw the fit
            # plot.save(output_filename)

        gr = ROOT.TGraphErrors(len(bin_centers), array('d', bin_centers), array('d', sigmas), array('d', [0]*len(bin_centers)), array('d', sigmas_unc))
        factor = 0.2
        gr_ideal = ROOT.TGraphErrors(len(bin_centers), array('d', bin_centers), array('d', [factor * pt for pt in bin_centers]), array('d', [0]*len(bin_centers)), array('d', [0]*len(bin_centers)))
        gr_cont = Contribution(gr, label='Measured')
        gr_ideal_cont = Contribution(gr_ideal, label=str(factor)+'*p_{T}', line_color=ROOT.kBlue, marker_color=ROOT.kBlue)
        plot = Plot([gr_cont, gr_ideal_cont], what='graph', xtitle="p_{T}^{Reco}", ytitle="#sigma [GeV]", ylim=[0, 100], xlim=[10, 4000])
        plot.plot()
        plot.set_logx()
        output_filename = os.path.join(plot_dir, "%s_%s_sigma_plot.%s" % (mydir, hist_name, OUTPUT_FMT))
        plot.save(output_filename)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input', 
                        nargs='+',
                        help='Input ROOT files to process. '
                        'Several dirs can be specified here, separated by a space.')
    parser.add_argument("-o", "--output", help="Directory to put output plot dirs into", default=None)
    args = parser.parse_args()
    print(args)

    for in_file in args.input:
        default_plot_dir = os.path.join(os.path.dirname(in_file), "resolution_"+os.path.splitext(os.path.basename(in_file))[0])
        plot_dir = args.output if args.output else default_plot_dir
        skip_dirs = [
            'Dijet_Presel',
            'Dijet_Presel_g_unknown',
            'Dijet_Presel_gg',
            'Dijet_Presel_gq',
            'Dijet_Presel_qq',
            'Dijet_Presel_qg',
            'Dijet_Presel_q_unknown',
            'Dijet_Presel_unknown',
            'Dijet_Presel_unknown_unknown',
            'Dijet_Presel_unknown_unknown',
            'Dijet_Presel_unknown_q',
            'Dijet_Presel_unknown_g',
            'Dijet_Presel_highPt',
            'Dijet_highPt',
            'Dijet_Presel_g_unknown_highPt',
            'Dijet_Presel_gg_highPt',
            'Dijet_Presel_gq_highPt',
            'Dijet_Presel_qq_highPt',
            'Dijet_Presel_qg_highPt',
            'Dijet_Presel_q_unknown_highPt',
            'Dijet_Presel_unknown_highPt',
            'Dijet_Presel_unknown_unknown_highPt',
            'Dijet_Presel_unknown_unknown_highPt',
            'Dijet_Presel_unknown_q_highPt',
            'Dijet_Presel_unknown_g_highPt',
            'ZPlusJets_Presel',
            'ZPlusJets_Presel_q',
            'ZPlusJets_Presel_g',
            'ZPlusJets_Presel_unknown',
            'SFrame',
            'cf_metfilters_raw',
            'cf_metfilters',
        ]
        if "_QCD_" in in_file:
            skip_dirs.append("ZPlusJets")
        do_projection_plots(in_file, plot_dir=plot_dir, skip_dirs=skip_dirs, do_fit=True)
