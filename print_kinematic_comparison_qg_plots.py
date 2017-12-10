#!/usr/bin/env python

"""Print plots comparing different kinematics between different flavours"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
from itertools import product
import numpy as np
np.seterr(all='raise')
from array import array

# My stuff
from comparator import Contribution, Plot, grab_obj
import qg_common as qgc
import qg_general_plots as qgg
import common_utils as cu


# For debugging
import sys
import tracers
# sys.settrace(tracers.trace_calls)
# sys.settrace(tracers.trace_calls_detail)
# sys.settrace(tracers.trace_calls_and_returns)

# import hunter
# hunter.trace(module='comparator')

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)

# Control output format
OUTPUT_FMT = "pdf"


def get_list_of_obj(directory):
    key_list = directory.GetListOfKeys()
    return [x.GetName() for x in key_list]


def do_all_1D_projection_plots_in_dir(directories, output_dir, components_styles_dicts=None, 
                                      draw_opts="NOSTACK HISTE", do_ratio=True, 
                                      normalise_hists=True, filter_noisy=True):
    """
    Given a set of TDirs, loop over all 2D hists, do projection hists for chosen bins, and plot all TDir contributions on a canvas for comparison.

    components_styles_dicts should be a list of style dicts, one for each directory/contribution to a plot. 
    This should include label for this component.

    filter_noisy aims to remove low stats/large error components that muddle things up
    """
    if len(directories) != len(components_styles_dicts):
        raise IndexError("Need same number of style dicts and directories")

    list_of_obj = [get_list_of_obj(d) for d in directories]
    print list_of_obj[0]

    # check all have same list of plots
    if not all(x == list_of_obj[0] for x in list_of_obj):
        print "Different number of object in the TDirectorys"

    common_list_of_obj = set(list_of_obj[0])
    for l in list_of_obj[1:]:
        common_list_of_obj = common_list_of_obj & set(l)


    pt_bins = [(20, 40), (40, 60), (60, 80), (100, 120), (160, 200), (260, 300), (500, 600), (1000, 2000)]

    for obj_name in common_list_of_obj:
        objs = [d.Get(obj_name) for d in directories]
        # Ignore TH1s
        if not isinstance(objs[0], (ROOT.TH2F, ROOT.TH2D, ROOT.TH2I)):
            print obj_name, "is not a TH2"
            continue
        
        if "flav" in obj_name:
            continue

        if obj_name in ['eta_jet1_vs_eta_jet2', 'phi_jet1_vs_pt_jet1', 'phi_jet2_vs_pt_jet1', 
                        'reliso_mu1_vs_pt_jet1', 'reliso_mu2_vs_pt_jet1', 'dphi_mumu_jet1_vs_pt_jet1', 
                        'dphi_mumu_vs_pt_jet1', 'pt_jet1_z_pt_jet2_z_ratio']:
            continue

        for pt_min, pt_max in pt_bins:
            # print pt_min, pt_max
            rebin = 1
            # exceptions...why didn't i pick the same number of bins...
            do_not_rebin = any([
                "n_jets" in obj_name, 
                "n_mu" in obj_name, 
                "met_sig" in obj_name, 
                obj_name.startswith('dphi_mumu'), 
                obj_name.startswith('pt_jet3_frac'), 
                obj_name.startswith('pt_jet1_jet2_ratio'),
                obj_name.startswith('pt_jet2_z_ratio')
            ])
            if not do_not_rebin:
                if objs[0].GetNbinsX() % 5 == 0:
                    rebin = 5
                elif objs[0].GetNbinsX() % 4 == 0:
                    rebin = 4
                elif objs[0].GetNbinsX() % 3 == 0:
                    rebin = 3
            contributions = [Contribution(qgg.get_projection_plot(ob, pt_min, pt_max), 
                                          normalise_hist=normalise_hists, rebin_hist=rebin, **csd) 
                             for ob, csd in zip(objs, components_styles_dicts)]

            if len(contributions) == 0:
                continue
            # Ignore empty objs
            total_entries = sum(c.obj.GetEntries() for c in contributions)
            if total_entries == 0:
                continue

            if filter_noisy:
                # filter contributions to ensure odd low stat ones don't dominate the scale
                # combination of factors to remove noisy samples
                mean_errs = []
                max_over_mean_errs = []
                rel_err_vars = []
                for cont in contributions:
                    obj = cont.obj
                    errs = [obj.GetBinError(i) for i in range(obj.GetNbinsX()+1)]
                    if obj.GetEntries() > 0:
                        mean_errs.append(np.mean(errs))
                        max_over_mean_errs.append(np.max(errs) / np.mean(errs))
                        rel_err_vars.append(np.std(errs) / np.mean(errs))
                        # print "obj", cont.label
                        # print 'mean errs', np.mean(errs)
                        # print 'rel max err', np.max(errs) / np.mean(errs)
                        # print 'stddev err', np.std(errs)
                        # print 'rel std dev err', np.std(errs) / np.mean(errs)
                    else:
                        # Dud values if 0 entries
                        mean_errs.append(9999999)
                        max_over_mean_errs.append(9999999)
                        rel_err_vars.append(9999999)
                        
                ref_mean_err = np.median(mean_errs)
                ref_rel_err_var = np.median(rel_err_vars)
                # print '-'*20
                # print 'mean mean err', ref_mean_err
                # print 'mean var', ref_rel_err_var
                contributions = [cont for cont, merr, rev, mom 
                                 in zip(contributions, mean_errs, rel_err_vars, max_over_mean_errs) 
                                 if (merr < 2.5*ref_mean_err) and (rev < 5*ref_rel_err_var or mom<5)]

            ylim = None
            # if "pt_jet" in obj_name and "ratio" not in obj_name and "frac" not in obj_name:
            #     ylim = [1E-9, 1E-1]
            title = "%d < p_{T}^{jet 1} < %d GeV" % (pt_min, pt_max)
            p = Plot(contributions, what='hist', ytitle="p.d.f." if normalise_hists else "N", 
                     title=title,
                     subplot_type="ratio" if do_ratio else None, 
                     subplot_title="#splitline{Ratio wrt}{%s}" % contributions[0].label,
                     subplot=contributions[0], ylim=ylim)
            p.legend.SetX1(0.7)
            p.legend.SetX2(0.9)
            p.plot(draw_opts)
            # EURGH FIXME
            # if "pt_jet" in obj_name and "ratio" not in obj_name and "frac" not in obj_name:
            #     p.set_logy()
            p.save(os.path.join(output_dir, obj_name+"_pt%dto%d.%s" % (pt_min, pt_max, OUTPUT_FMT)))


def do_dijet_distributions(root_dir, dir_append=""):
    """Do plots comparing different jet flavs in dijet region"""
    dir_names = ["Dijet_Presel_gg", "Dijet_Presel_qg", "Dijet_Presel_gq", "Dijet_Presel_qq", 
                 "Dijet_Presel_q_unknown", "Dijet_Presel_g_unknown",
                 "Dijet_Presel_unknown_q", "Dijet_Presel_unknown_g",
                 "Dijet_Presel_unknown_unknown"
                ]
    remove_unknowns = dir_append == "_highPt"
    dir_names = [d + dir_append for d in dir_names if not ("unknown" in d and remove_unknowns)]
    root_file = cu.open_root_file(os.path.join(root_dir, qgc.QCD_FILENAME))
    
    directories = [cu.get_from_file(root_file, dn) for dn in dir_names]
    gg_col = ROOT.kRed
    qg_col = ROOT.kGreen+2
    gq_col = ROOT.kBlack
    qq_col = ROOT.kBlue
    unknown_cols = [ROOT.kOrange+1, ROOT.kOrange+4, ROOT.kPink+6, ROOT.kViolet, ROOT.kAzure+1]
    csd = [
        {"label": "gg", "line_color": gg_col, "fill_color": gg_col, "marker_color": gg_col, "marker_style": 20},
        {"label": "qg", "line_color": qg_col, "fill_color": qg_col, "marker_color": qg_col, "marker_style": 21},
        {"label": "gq", "line_color": gq_col, "fill_color": gq_col, "marker_color": gq_col, "marker_style": 22},
        {"label": "qq", "line_color": qq_col, "fill_color": qq_col, "marker_color": qq_col, "marker_style": 23},
        {"label": "1:q  2:unknown", "line_color": unknown_cols[0], "fill_color": unknown_cols[0], "marker_color": unknown_cols[0], "marker_style": 29, "marker_size": 1.3},
        {"label": "1:g  2:unknown", "line_color": unknown_cols[1], "fill_color": unknown_cols[1], "marker_color": unknown_cols[1], "marker_style": 33, "marker_size": 1.3},
        {"label": "1:unknown  2:q", "line_color": unknown_cols[2], "fill_color": unknown_cols[2], "marker_color": unknown_cols[2], "marker_style": 41, "marker_size": 1.3},
        {"label": "1:unknown  2:g", "line_color": unknown_cols[3], "fill_color": unknown_cols[3], "marker_color": unknown_cols[3], "marker_style": 34, "marker_size": 1.2},
        {"label": "1:unknown  2:unknown", "line_color": unknown_cols[4], "fill_color": unknown_cols[4], "marker_color": unknown_cols[4], "marker_style": 47, "marker_size": 1.2},
    ]
    csd = [c for c in csd if not ("unknown" in c['label'] and remove_unknowns)]

    # Compare shapes
    do_all_1D_projection_plots_in_dir(directories=directories, 
                                      output_dir=os.path.join(root_dir, "Dijet_kin_comparison_normalised%s" % dir_append),
                                      components_styles_dicts=csd)
    # Compare yields
    do_all_1D_projection_plots_in_dir(directories=directories, 
                                      output_dir=os.path.join(root_dir, "Dijet_kin_comparison_absolute%s" % dir_append),
                                      components_styles_dicts=csd,
                                      normalise_hists=False,
                                      filter_noisy=False)

    # do jet1 vs jet2 flav
    output_filename = os.path.join(root_dir, "Dijet_kin_comparison_2d", "flav_jet1_jet2%s.%s" % (dir_append, OUTPUT_FMT))
    h2d = cu.get_from_file(root_file, "Dijet_Presel%s/flav_jet1_jet2" % dir_append)
    if h2d.Integral() > 0:
        h2d.Scale(1./h2d.Integral())
        qgg.do_2D_plot(h2d, output_filename, draw_opt="COLZ", logz=True, zlim=[1E-4, 1])
    
    # do jet1 vs jet2 eta
    ellipse_eff_dict = {}
    deta_eff_dict = {}
    for dname in dir_names:
        output_filename = os.path.join(root_dir, "Dijet_kin_comparison_2d%s" % dir_append, "eta_jet1_eta_jet2_%s.%s" % (dname.replace("Dijet_Presel_", ""), OUTPUT_FMT))
        h2d = cu.get_from_file(root_file, "%s/eta_jet1_vs_eta_jet2" % dname)
        h2d.Rebin2D(2, 2)
        # h2d.Scale(1./h2d.Integral())
        title = dname.replace("Dijet_Presel_", "")
        vline = ROOT.TLine(0, -5, 0, 5)
        vline.SetLineStyle(2)
        hline = ROOT.TLine(-5, 0, 5, 0)
        hline.SetLineStyle(2)
        lines = [hline, vline]
        print dname
        """
        # Do ellipse efficiencies for diff size ellipses
        effs = []
        for ind, factor in enumerate(np.arange(0.1, 2.2, 0.1), 1):
            ellipse = ROOT.TEllipse(0, 0, factor*1.7*np.sqrt(2), factor*1*np.sqrt(2), 0, 360, 45)
            integral = get_integral_under_ellipse(h2d, ellipse)
            print "Fraction under ellipse with factor", factor, " = ", integral / h2d.Integral()
            ellipse.SetFillStyle(0)
            ellipse.SetLineStyle(ind)
            lines.append(ellipse)
            effs.append(integral)
        ellipse_eff_dict[dname] = effs
        """
        qgg.do_2D_plot(h2d, output_filename, draw_opt="COLZ", logz=False, title=title, other_things_to_draw=lines)
        qgg.do_2D_plot(h2d, output_filename.replace(".%s" % OUTPUT_FMT, "_logZ.%s" % OUTPUT_FMT), 
                       draw_opt="COLZ", logz=True, title=title, other_things_to_draw=lines)
        """ 
        effs = []
        for ind, deta in enumerate(np.arange(0.1, 2.5, 0.1), 1):
            integral = get_integral_under_deta(h2d, deta)
            effs.append(integral)
        deta_eff_dict[dname] = effs
        """ 

    """
    output_filename = os.path.join(root_dir, "Dijet_kin_comparison_2d", "ellipse_roc_jet2.%s" % OUTPUT_FMT)
    do_roc_plot(ellipse_eff_dict, output_filename)

    output_filename = os.path.join(root_dir, "Dijet_kin_comparison_2d", "deta_roc_jet2.%s" % OUTPUT_FMT)
    do_roc_plot(deta_eff_dict, output_filename)
    """

def do_roc_plot(eff_dict, output_filename):
    """Turn dict of efficiencies into ROC plot
    
    Parameters
    ----------
    eff_dict : TYPE
        Description
    output_filename : TYPE
        Description
    """
    signal_effs = [0] * len(eff_dict.values()[0])
    bkg_effs = [0] * len(eff_dict.values()[0])
    for k, v in eff_dict.iteritems():
        if k.replace("Dijet_Presel_", "").lstrip("_unknown_").lstrip("_q").startswith("g"):
            print k
            for ind, eff in enumerate(v):
                signal_effs[ind] += eff
        else:
            for ind, eff in enumerate(v):
                bkg_effs[ind] += eff
    # divide by totals
    for ind, (s, b) in enumerate(zip(signal_effs, bkg_effs)):
        total = s + b
        signal_effs[ind] /= total
        bkg_effs[ind] /= total
    print signal_effs
    print bkg_effs

    gr = ROOT.TGraph(len(signal_effs), array('d', bkg_effs), array('d', signal_effs))
    cont = Contribution(gr, marker_style=21)
    p = Plot([cont], "graph", xtitle="fraction jet2!=g", ytitle="fraction jet2=g", xlim=[0, 1], ylim=[0, 1], legend=False)
    p.plot("AP")
    p.save(output_filename)


def get_integral_under_deta(hist, deta):
    integral = 0
    for xbin in range(1, hist.GetNbinsX()+1):
        for ybin in range(1, hist.GetNbinsY()+1):
            contents = hist.GetBinContent(xbin, ybin)
            if contents == 0:
                continue
            x_center = hist.GetXaxis().GetBinCenter(xbin)
            y_center = hist.GetYaxis().GetBinCenter(ybin)
            if (abs(y_center - x_center) < deta):
                integral += contents
    return integral


def get_integral_under_ellipse(hist, ellipse):
    """Get integral of 2d hist under ellipse
    
    Parameters
    ----------
    hist : TH2
        Description
    ellipse : Tellipse
        Description
    
    Returns
    -------
    float
    """
    old_fill_style = ellipse.GetFillStyle()
    ellipse.SetFillStyle(1001)
    integral = 0
    for xbin in range(1, hist.GetNbinsX()+1):
        for ybin in range(1, hist.GetNbinsY()+1):
            contents = hist.GetBinContent(xbin, ybin)
            if contents == 0:
                continue
            x_center = hist.GetXaxis().GetBinCenter(xbin)
            y_center = hist.GetYaxis().GetBinCenter(ybin)
            a = ellipse.GetR1()
            b = ellipse.GetR2()
            in_ellipse = ((pow((x_center+y_center), 2) / (a*a)) + (pow((x_center-y_center), 2) / (b*b))) <= 2;
            if in_ellipse:
                integral += contents
    ellipse.SetFillStyle(old_fill_style)
    return integral


def do_zpj_distributions(root_dir, dir_append=""):
    """Do plots comparing different jet flavs in z+jets region"""
    dir_names = ["ZPlusJets_Presel_q", "ZPlusJets_Presel_g", "ZPlusJets_Presel_unknown"]
    dir_names = [d+dir_append for d in dir_names]
    root_file = cu.open_root_file(os.path.join(root_dir, qgc.DY_FILENAME))
    directories = [cu.get_from_file(root_file, dn) for dn in dir_names]
    g_col = ROOT.kRed
    q_col = ROOT.kBlue
    unknown_col = ROOT.kViolet
    csd = [
        {"label": "q", "line_color": q_col, "fill_color": q_col, "marker_color": q_col, "marker_style": 20},
        {"label": "g", "line_color": g_col, "fill_color": g_col, "marker_color": g_col, "marker_style": 21},
        {"label": "unknown", "line_color": unknown_col, "fill_color": unknown_col, "marker_color": unknown_col, "marker_style": 22}
    ]
    # Compare shapes
    do_all_1D_projection_plots_in_dir(directories=directories, 
                                      output_dir=os.path.join(root_dir, "ZpJ_kin_comparison_normalised%s" % dir_append),
                                      components_styles_dicts=csd)
    # Compare yields
    do_all_1D_projection_plots_in_dir(directories=directories, 
                                      output_dir=os.path.join(root_dir, "ZpJ_kin_comparison_absolute%s" % dir_append),
                                      components_styles_dicts=csd,
                                      normalise_hists=False,
                                      filter_noisy=False)
    
    for dname in dir_names:
        output_filename = os.path.join(root_dir, "ZpJ_kin_comparison_2d%s" % dir_append, "pt_jet1_z_pt_jet2_z_ratio_%s.%s" % (dname.replace("ZPlusJets_Presel_", ""), OUTPUT_FMT))
        h2d = cu.get_from_file(root_file, "%s/pt_jet1_z_pt_jet2_z_ratio" % dname)
        # h2d.Scale(1./h2d.Integral())
        title = dname.replace("ZPlusJets_Presel_", "")
        qgg.do_2D_plot(h2d, output_filename, draw_opt="COLZ", logz=False, title=title)
        qgg.do_2D_plot(h2d, output_filename.replace(".%s" % OUTPUT_FMT, "_logZ.%s" % OUTPUT_FMT), 
                       draw_opt="COLZ", logz=True, title=title)

if __name__ == "__main__":
    parser = qgc.get_parser()
    args = parser.parse_args()

    for workdir in args.workdirs:
        if os.path.isfile(os.path.join(workdir, qgc.QCD_FILENAME)):
            do_dijet_distributions(workdir)
            do_dijet_distributions(workdir, "_highPt")
        if os.path.isfile(os.path.join(workdir, qgc.DY_FILENAME)):
            do_zpj_distributions(workdir)
            do_zpj_distributions(workdir, "_highPt")

