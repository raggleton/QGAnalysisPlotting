#!/usr/bin/env python

"""Print QG plots"""

import ROOT
from comparator import Contribution, Plot, grab_obj
from TDRStyle import TDR_Style
from MyStyle import My_Style
import common_utils as cu
from uuid import uuid4
import bisect
import numpy as np


ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
# TDR_Style.cd()
My_Style.cd()


DY_COLOUR = 880
QCD_COLOUR = 867

DY_ZpJ_LABEL = "DY+j, Z+jets selection"
DY_ZpJ_GEN_LABEL = "DY+j, Z+jets selection (GenJets)"
DY_ZpJ_QFLAV_LABEL = "DY+j, Z+jets selection (uds-matched)"
DY_ZpJ_GFLAV_LABEL = "DY+j, Z+jets selection (g-matched)"

DY_Dijet_LABEL = "DY+j, Dijet selection"
DY_Dijet_GEN_LABEL = "DY+j, Dijet selection (GenJets)"
DY_Dijet_QFLAV_LABEL = "DY+j, Dijet selection (uds-matched)"
DY_Dijet_GFLAV_LABEL = "DY+j, Dijet selection (g-matched)"

QCD_ZpJ_LABEL = "QCD, Z+jets selection"
QCD_ZpJ_GEN_LABEL = "QCD, Z+jets selection (GenJets)"
QCD_ZpJ_QFLAV_LABEL = "QCD, Z+jets selection (uds-matched)"
QCD_ZpJ_GFLAV_LABEL = "QCD, Z+jets selection (g-matched)"

QCD_Dijet_LABEL = "QCD, Dijet selection"
QCD_Dijet_GEN_LABEL = "QCD, Dijet selection (GenJets)"
QCD_Dijet_QFLAV_LABEL = "QCD, Dijet selection (uds-matched)"
QCD_Dijet_GFLAV_LABEL = "QCD, Dijet selection (g-matched)"


ROOT_DIR = "workdir_ak4chs"
# ROOT_DIR = "workdir_ak4puppi"

TITLE_STR = "[%s]" % ROOT_DIR.replace("workdir_", "")


def plot_dy_vs_qcd(dy_obj, qcd_obj, output_filename, xtitle=None, title=None, rebin=1, dy_kwargs=None, qcd_kwargs=None):
    """Plot DY vs QCD sample"""
    dy_kwargs = dy_kwargs or {}
    qcd_kwargs = qcd_kwargs or {}
    DYJ = Contribution(dy_obj, normalise_hist=True, fill_style=1, rebin_hist=rebin, **dy_kwargs)
    QCD = Contribution(qcd_obj, normalise_hist=True, fill_style=1, rebin_hist=rebin, **qcd_kwargs)
    p = Plot([DYJ, QCD], "hist", ratio_subplot=DYJ, xtitle=xtitle, ytitle="p.d.f", title=title)
    draw_opt = "NOSTACK HISTE"
    p.plot(draw_opt)
    hst = p.container
    hst.GetYaxis().SetTitleOffset(1.2)
    hst.GetXaxis().SetTitleOffset(1.1)
    p.plot(draw_opt)
    p.save(output_filename)


def do_all_inclusive_plots():
    """Do plots inclusive over all pt, eta, nvtx"""

    for v in ['LHA', 'pTD', 'width', 'thrust', 'multiplicity']:
        rebin = 1
        if v == "multiplicity":
            rebin = 2

        plot_dir = "plots_dy_vs_qcd"
        dy_kwargs = dict(line_color=DY_COLOUR, fill_color=DY_COLOUR)
        qcd_kwargs = dict(line_color=QCD_COLOUR, fill_color=QCD_COLOUR)

        # all reco jets
        dy_kwargs['label'] = DY_ZpJ_LABEL
        qcd_kwargs['label'] = QCD_Dijet_LABEL
        plot_dy_vs_qcd(dy_obj=grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "ZPlusJets_QG/jet_%s" % v),
                       qcd_obj=grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "Dijet_QG/jet_%s" % v),
                       output_filename="%s/%s/jet_%s.pdf" % (ROOT_DIR, plot_dir, v),
                       rebin=rebin, dy_kwargs=dy_kwargs, qcd_kwargs=qcd_kwargs)

        # "correct" flav matched reco jets
        dy_kwargs['label'] = DY_ZpJ_QFLAV_LABEL
        qcd_kwargs['label'] = QCD_Dijet_GFLAV_LABEL
        plot_dy_vs_qcd(dy_obj=grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "ZPlusJets_QG/qjet_%s" % v),
                       qcd_obj=grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "Dijet_QG/gjet_%s" % v),
                       output_filename="%s/%s/jet_%s_flav_matched.pdf" % (ROOT_DIR, plot_dir, v),
                       rebin=rebin, dy_kwargs=dy_kwargs, qcd_kwargs=qcd_kwargs)

        # gen jets (no flav matching)
        dy_kwargs['label'] = DY_ZpJ_GEN_LABEL
        qcd_kwargs['label'] = QCD_Dijet_GEN_LABEL
        plot_dy_vs_qcd(dy_obj=grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "ZPlusJets_QG/genjet_%s" % v),
                       qcd_obj=grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "Dijet_QG/genjet_%s" % v),
                       output_filename="%s/%s/genjet_%s.pdf" % (ROOT_DIR, plot_dir, v),
                       rebin=rebin, dy_kwargs=dy_kwargs, qcd_kwargs=qcd_kwargs)


def do_2D_plot(obj, output_filename, renorm_axis=None, title=None, rebin=None):
    if rebin:
        obj.Rebin2D(*rebin)
    if renorm_axis:
        obj_renorm = cu.make_normalised_TH2(obj, renorm_axis)
    else:
        obj_renorm = obj
    if title:
        obj_renorm.SetTitle(title)
    canvas = ROOT.TCanvas("canv%s" % uuid4(), "", 800, 600)
    canvas.SetTicks(1, 1)
    obj_renorm.Draw("COLZ")
    obj_renorm.GetYaxis().SetTitleOffset(1.3)
    obj_renorm.GetXaxis().SetTitleOffset(1.1)
    canvas.SaveAs(output_filename)


def do_all_2D_plots():
    """Do 2D distributions"""
    for v in ['jet_LHA', 'jet_pTD', 'jet_width', 'jet_thrust', 'jet_multiplicity', 'jet_flavour', "jet_genParton_flavour"]:
        v += "_vs_pt"
        rebin = [1, 4]
        if v == "jet_multiplicity_vs_pt":
            rebin = [2, 4]
        for rn in ['Y']:  # different renormalisation axes
            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "ZPlusJets_QG/%s" % v),
                       output_filename="%s/plots_2d/dy_zpj_%s_norm%s.pdf" % (ROOT_DIR, v, rn),
                       renorm_axis=rn, title=DY_ZpJ_LABEL, rebin=rebin)
            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "Dijet_QG/%s" % v),
                       output_filename="%s/plots_2d/dy_dijet_%s_norm%s.pdf" % (ROOT_DIR, v, rn),
                       renorm_axis=rn, title=DY_Dijet_LABEL, rebin=rebin)
            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "ZPlusJets_QG/%s" % v),
                       output_filename="%s/plots_2d/qcd_zpj_%s_norm%s.pdf" % (ROOT_DIR, v, rn),
                       renorm_axis=rn, title=QCD_ZpJ_LABEL, rebin=rebin)
            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "Dijet_QG/%s" % v),
                       output_filename="%s/plots_2d/qcd_dijet_%s_norm%s.pdf" % (ROOT_DIR, v, rn),
                       renorm_axis=rn, title=QCD_Dijet_LABEL, rebin=rebin)

            if v in ["jet_flavour_vs_pt", "jet_genParton_flavour_vs_pt"]:
                continue

            # Flavour matched reco
            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "ZPlusJets_QG/q%s" % v),
                       output_filename="%s/plots_2d/dy_zpj_%s_norm%s_flavMatched.pdf" % (ROOT_DIR, v, rn),
                       renorm_axis=rn, title=DY_ZpJ_QFLAV_LABEL, rebin=rebin)
            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "Dijet_QG/g%s" % v),
                       output_filename="%s/plots_2d/dy_dijet_%s_norm%s_flavMatched.pdf" % (ROOT_DIR, v, rn),
                       renorm_axis=rn, title=DY_Dijet_GFLAV_LABEL, rebin=rebin)
            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "ZPlusJets_QG/q%s" % v),
                       output_filename="%s/plots_2d/qcd_zpj_%s_norm%s_flavMatched.pdf" % (ROOT_DIR, v, rn),
                       renorm_axis=rn, title=QCD_ZpJ_QFLAV_LABEL, rebin=rebin)
            do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "Dijet_QG/g%s" % v),
                       output_filename="%s/plots_2d/qcd_dijet_%s_norm%s_flavMatched.pdf" % (ROOT_DIR, v, rn),
                       renorm_axis=rn, title=QCD_Dijet_GFLAV_LABEL, rebin=rebin)

            # GenJets
            # do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "ZPlusJets_QG/gen%s" % v),
            #            output_filename="%s/plots_2d/dy_zpj_%s_norm%s_genJets.pdf" % (ROOT_DIR, v, rn),
            #            renorm_axis=rn, title=DY_ZpJ_GEN_LABEL)
            # do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "Dijet_QG/gen%s" % v),
            #            output_filename="%s/plots_2d/dy_dijet_%s_norm%s_genJets.pdf" % (ROOT_DIR, v, rn),
            #            renorm_axis=rn, title=DY_Dijet_GEN_LABEL)
            # do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "ZPlusJets_QG/gen%s" % v),
            #            output_filename="%s/plots_2d/qcd_zpj_%s_norm%s_genJets.pdf" % (ROOT_DIR, v, rn),
            #            renorm_axis=rn, title=QCD_ZpJ_GEN_LABEL)
            # do_2D_plot(grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "Dijet_QG/gen%s" % v),
            #            output_filename="%s/plots_2d/qcd_dijet_%s_norm%s_genJets.pdf" % (ROOT_DIR, v, rn),
            #            renorm_axis=rn, title=QCD_Dijet_GEN_LABEL)



def do_all_exclusive_plots():
    """Do pt/eta/nvtx binned 1D plots"""
    plot_dir = "plots_dy_vs_qcd"
    dy_kwargs = dict(line_color=DY_COLOUR, fill_color=DY_COLOUR)
    qcd_kwargs = dict(line_color=QCD_COLOUR, fill_color=QCD_COLOUR)
    for v in ['jet_LHA', 'jet_pTD', 'jet_width', 'jet_thrust', 'jet_multiplicity', 'jet_flavour', "jet_genParton_flavour"]:
        v += "_vs_pt"

        # All reco jets
        dy_kwargs['label'] = DY_ZpJ_LABEL
        qcd_kwargs['label'] = QCD_Dijet_LABEL
        h2d_dyj = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "ZPlusJets_QG/%s" % v)
        h2d_qcd = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "Dijet_QG/%s" % v)

        y_axis = h2d_dyj.GetYaxis()
        bin_edges = [y_axis.GetBinLowEdge(i) for i in range(1, y_axis.GetNbins()+1)]

        def do_projection_plot(start_val, end_val, h2d_dyj, h2d_qcd, output_filename, rebin_y=1):
            """Find the best set of y bin edges fro given start_val & end_val, and plot the projection onto the x axis"""
            bin_start = bisect.bisect_right(bin_edges, start_val)
            bin_end = bisect.bisect_right(bin_edges, end_val)
            hproj_dyj = h2d_dyj.ProjectionX(v+"_dyj_projX", bin_start, bin_end, "eo")
            hproj_qcd = h2d_qcd.ProjectionX(v+"_qcd_projX", bin_start, bin_end, "eo")
            plot_dy_vs_qcd(hproj_dyj, hproj_qcd, output_filename.format(start_val=start_val, end_val=end_val),
                           xtitle=None, title="%d < p_{T}^{jet} < %d GeV %s" % (start_val, end_val, TITLE_STR),
                           rebin=rebin_y, dy_kwargs=dy_kwargs, qcd_kwargs=qcd_kwargs)

        rebin = 2 if v=="jet_multiplicity_vs_pt" else 1
        do_projection_plot(100, 200, h2d_dyj, h2d_qcd, "%s/%s/ptBinned/%s_pt{start_val}to{end_val}.pdf" % (ROOT_DIR, plot_dir, v), rebin)
        do_projection_plot(100, 2000, h2d_dyj, h2d_qcd, "%s/%s/ptBinned/%s_pt{start_val}to{end_val}.pdf" % (ROOT_DIR, plot_dir, v), rebin)
        do_projection_plot(400, 500, h2d_dyj, h2d_qcd, "%s/%s/ptBinned/%s_pt{start_val}to{end_val}.pdf" % (ROOT_DIR, plot_dir, v), rebin)
        do_projection_plot(1000, 2000, h2d_dyj, h2d_qcd, "%s/%s/ptBinned/%s_pt{start_val}to{end_val}.pdf" % (ROOT_DIR, plot_dir, v), rebin)

        if v in ["jet_flavour_vs_pt", "jet_genParton_flavour_vs_pt"]:
            continue

        # Flavour matched jets
        dy_kwargs['label'] = DY_ZpJ_QFLAV_LABEL
        qcd_kwargs['label'] = QCD_Dijet_GFLAV_LABEL
        h2d_dyj = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "ZPlusJets_QG/q%s" % v)
        h2d_qcd = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "Dijet_QG/g%s" % v)

        do_projection_plot(100, 200, h2d_dyj, h2d_qcd, "%s/%s/ptBinned/%s_pt{start_val}to{end_val}_flavMatched.pdf" % (ROOT_DIR, plot_dir, v), rebin)
        do_projection_plot(100, 2000, h2d_dyj, h2d_qcd, "%s/%s/ptBinned/%s_pt{start_val}to{end_val}_flavMatched.pdf" % (ROOT_DIR, plot_dir, v), rebin)
        do_projection_plot(400, 500, h2d_dyj, h2d_qcd, "%s/%s/ptBinned/%s_pt{start_val}to{end_val}_flavMatched.pdf" % (ROOT_DIR, plot_dir, v), rebin)
        do_projection_plot(1000, 2000, h2d_dyj, h2d_qcd, "%s/%s/ptBinned/%s_pt{start_val}to{end_val}_flavMatched.pdf" % (ROOT_DIR, plot_dir, v), rebin)

        # Gen Jets
        dy_kwargs['label'] = DY_ZpJ_GEN_LABEL
        qcd_kwargs['label'] = QCD_Dijet_GEN_LABEL
        h2d_dyj = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "ZPlusJets_QG/gen%s" % v)
        h2d_qcd = grab_obj("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "Dijet_QG/gen%s" % v)

        do_projection_plot(100, 200, h2d_dyj, h2d_qcd, "%s/%s/ptBinned/gen%s_pt{start_val}to{end_val}.pdf" % (ROOT_DIR, plot_dir, v), rebin)
        do_projection_plot(100, 2000, h2d_dyj, h2d_qcd, "%s/%s/ptBinned/gen%s_pt{start_val}to{end_val}.pdf" % (ROOT_DIR, plot_dir, v), rebin)
        do_projection_plot(400, 500, h2d_dyj, h2d_qcd, "%s/%s/ptBinned/gen%s_pt{start_val}to{end_val}.pdf" % (ROOT_DIR, plot_dir, v), rebin)
        do_projection_plot(1000, 2000, h2d_dyj, h2d_qcd, "%s/%s/ptBinned/gen%s_pt{start_val}to{end_val}.pdf" % (ROOT_DIR, plot_dir, v), rebin)


def do_all_flavour_fraction_plots():
    """Do plots of flavour fraction of jets as a function of PT, for various samples"""

    def get_flavour_fractions(input_file, selection, which=None):
        which = which or ""
        h2d_flav = grab_obj(input_file, "%s_QG/jet_%sflavour_vs_pt" % (selection, which))

        h2d_flav.Rebin2D(1, 5)
        y_axis = h2d_flav.GetYaxis()
        pt_bins_lower, pt_bins_upper = [], []
        flav_dict = {'d': [], 'u': [], 's': [], 'c': [], 'b': [] ,'t': [], 'g': []}
        for i in range(1, y_axis.GetNbins()-1):

            d_frac = h2d_flav.GetBinContent(2, i)
            u_frac = h2d_flav.GetBinContent(3, i)
            s_frac = h2d_flav.GetBinContent(4, i)
            c_frac = h2d_flav.GetBinContent(5, i)
            b_frac = h2d_flav.GetBinContent(6, i)
            t_frac = h2d_flav.GetBinContent(7, i)
            g_frac = h2d_flav.GetBinContent(22, i)

            total = d_frac+u_frac+s_frac+c_frac+b_frac+t_frac+g_frac

            if total == 0:
                continue

            pt_bins_lower.append(y_axis.GetBinLowEdge(i))
            pt_bins_upper.append(y_axis.GetBinLowEdge(i+1))

            flav_dict['d'].append(d_frac / total)
            flav_dict['u'].append(u_frac / total)
            flav_dict['s'].append(s_frac / total)
            flav_dict['c'].append(c_frac / total)
            flav_dict['b'].append(b_frac / total)
            flav_dict['t'].append(t_frac / total)
            flav_dict['g'].append(g_frac / total)

        x_bins = [0.5 * (x1+x2) for x1,x2 in zip(pt_bins_lower[:], pt_bins_upper[:])]

        return x_bins, flav_dict

    def do_flavour_fraction_vs_pt(input_file, selection, output_filename, title=None, which=None):
        """Plot flavour fractions vs PT for one input file & selection"""
        x_bins, flav_dict = get_flavour_fractions(input_file, selection, which)
        # TODO: check if empy arrays
        gr_flav_u = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['u']))
        gr_flav_d = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['d']))
        gr_flav_s = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['s']))
        gr_flav_c = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['c']))
        gr_flav_b = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['b']))
        gr_flav_g = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(flav_dict['g']))
        # gr_flav_uds = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array([u+d+s for u, d, s in zip(flav_dict['u'], flav_dict['d'], flav_dict['s'])]))

        plot_u = Contribution(gr_flav_u, label="u frac", line_color=ROOT.kRed, marker_color=ROOT.kRed)
        plot_d = Contribution(gr_flav_d, label="d frac", line_color=ROOT.kBlue, marker_color=ROOT.kBlue)
        plot_s = Contribution(gr_flav_s, label="s frac", line_color=ROOT.kBlack, marker_color=ROOT.kBlack)
        plot_c = Contribution(gr_flav_c, label="c frac", line_color=ROOT.kGreen-3, marker_color=ROOT.kGreen-3)
        plot_b = Contribution(gr_flav_b, label="b frac", line_color=ROOT.kOrange, marker_color=ROOT.kOrange)
        plot_g = Contribution(gr_flav_g, label="g frac", line_color=ROOT.kViolet, marker_color=ROOT.kViolet)
        # plot_uds = Contribution(gr_flav_uds, label="uds frac", line_color=ROOT.kOrange+2, marker_color=ROOT.kOrange+2)

        p_flav = Plot([plot_u, plot_d, plot_s, plot_g, plot_c, plot_b], what='graph', xtitle="p_{T}^{jet} [GeV]", title=title)
        p_flav.plot("ALP")
        p_flav.save(output_filename)

    do_flavour_fraction_vs_pt("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "Dijet", "%s/qcd_dijet_flav_frac.pdf" % (ROOT_DIR), title=QCD_Dijet_LABEL)
    do_flavour_fraction_vs_pt("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "Dijet", "%s/dy_dijet_flav_frac.pdf" % (ROOT_DIR), title=DY_Dijet_LABEL)
    do_flavour_fraction_vs_pt("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "ZPlusJets", "%s/dy_zpj_flav_frac.pdf" % (ROOT_DIR), title=DY_ZpJ_LABEL)

    do_flavour_fraction_vs_pt("%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR, "Dijet", "%s/qcd_dijet_genParton_flav_frac.pdf" % (ROOT_DIR), title=QCD_Dijet_LABEL, which="genParton_")
    do_flavour_fraction_vs_pt("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "Dijet", "%s/dy_dijet_genParton_flav_frac.pdf" % (ROOT_DIR), title=DY_Dijet_LABEL, which="genParton_")
    do_flavour_fraction_vs_pt("%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR, "ZPlusJets", "%s/dy_zpj_genParton_flav_frac.pdf" % (ROOT_DIR), title=DY_ZpJ_LABEL, which="genParton_")

    def compare_flavour_fractions_vs_pt(input_files, selections, labels, flav, output_filename, title=None, which=None):
        """Plot a specified flavour fraction vs pT for several sources.
        Each entry in input_files, selections, and labels corresponds to one line"""
        info = [get_flavour_fractions(ifile, sel, which) for ifile, sel in zip(input_files, selections)]
        contribs = []
        for i, (x_bins, fdict) in enumerate(info):
            if flav in ['u', 'd', 's', 'c', 'b', 't', 'g']:
                gr = ROOT.TGraph(len(x_bins), np.array(x_bins), np.array(fdict[flav]))
            else:
                gr = ROOT.TGraph(len(x_bins), np.array(x_bins), 1.-np.array(fdict[flav.replace("1-", '')]))
            c = Contribution(gr, label="%s" % (labels[i]), line_style=i+1)
            contribs.append(c)
        ytitle = "%s flavour fraction" % flav
        p = Plot(contribs, what='graph', xtitle="p_{T}^{jet} [GeV]", ytitle=ytitle, title=title)
        p.plot("ALP")
        p.save(output_filename)

    input_files = [
        "%s/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" % ROOT_DIR,
        "%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR,
        "%s/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" % ROOT_DIR
    ]
    compare_flavour_fractions_vs_pt(input_files, ["Dijet", "Dijet", "ZPlusJets"], [QCD_Dijet_LABEL, DY_Dijet_LABEL, DY_ZpJ_LABEL], "1-g", "%s/compare_q_frac.pdf" % ROOT_DIR, title=TITLE_STR)
    compare_flavour_fractions_vs_pt(input_files, ["Dijet", "Dijet", "ZPlusJets"], [QCD_Dijet_LABEL, DY_Dijet_LABEL, DY_ZpJ_LABEL], "1-g", "%s/compare_q_frac_genParton.pdf" % ROOT_DIR, title=TITLE_STR, which="genParton_")


if __name__ == '__main__':
    # do_all_inclusive_plots()
    # do_all_2D_plots()
    do_all_exclusive_plots()
    do_all_flavour_fraction_plots()
