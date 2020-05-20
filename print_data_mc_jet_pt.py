#!/usr/bin/env python

"""
Do jet pT distributions for paper plots
"""

import os
os.nice(10)
import ROOT
from MyStyle import My_Style
My_Style.cd()
import sys
from array import array
import numpy as np

# My stuff
from comparator import Contribution, Plot, grab_obj
import qg_common as qgc
import qg_general_plots as qgg
import common_utils as cu

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gErrorIgnoreLevel = ROOT.kWarning


# Control output format
OUTPUT_FMT = "pdf"


def compare_bins(histA, histB):
    h1Array = histA.GetXaxis().GetXbins()
    h2Array = histB.GetXaxis().GetXbins()
    if h1Array.fN != h2Array.fN:
        raise ValueError("h1Array.fN != h2Array.fN")

    for i in range(h1Array.fN):
        if not ROOT.TMath.AreEqualRel(h1Array.GetAt(i), h2Array.GetAt(i), 1E-10):
            raise ValueError("not equal: %g %g" % (h1Array.GetAt(i), h2Array.GetAt(i)))

def do_jet_pt_plot(entries,
                   output_filename,
                   xlim=None,
                   ylim=None,
                   title="",
                   rebin=1,
                   data_first=True,
                   normalise_hists=True,
                   subplot_limits=None
                   ):
    # entries = [ent for ent in entries if ent]

    conts = [Contribution(ent[0], rebin_hist=rebin, **ent[1])  # Don't use noramlise_hist - screws up with unequal binning
             for ent in entries]

    do_legend = len(conts) > 1
    if len(conts) == 0:
        raise RuntimeError("0 contributions for this plot")

    if data_first:
        # Modify data Contribution to have no horizontal lines
        # We'll redraw it at the end with the proper bars ontop
        conts[0].line_width = 0
        conts[0].marker_size = 0.01
        conts[0]._update_styles()

        # For subplot to ensure only MC errors drawn, not MC+data
        data_obj = entries[0][0]
        # data_obj = conts[0].obj
        data_no_errors = data_obj.Clone()
        cu.remove_th1_errors(data_no_errors)
        for i, cont in enumerate(conts[1:], 1):
            if cont.subplot == data_obj:
                cont.subplot = data_no_errors
                # compare_bins(cont.obj, cont.subplot)

    plot = Plot(conts,
                what='hist',
                legend=do_legend,
                has_data=data_first,
                xlim=xlim,
                ylim=ylim,
                xtitle="p_{T}^{jet} [GeV]",
                # ytitle="#DeltaN/N" if normalise_hists else "N",
                ytitle="#frac{dN}{dp_{T}} [events / GeV]",
                title=title,
                subplot_type='ratio',
                subplot_title="Simulation / data",
                subplot_limits=subplot_limits)
    plot.y_padding_max_log = 1E4
    plot.legend.SetX1(0.6)
    plot.legend.SetY1(0.68)
    plot.legend.SetX2(0.98)
    plot.legend.SetY2(0.88)
    plot.left_margin = 0.16

    if data_first:
        # we'll do the filling of legend ourselves
        plot.do_legend = False

    draw_opt = "NOSTACKE HIST E ]["
    plot.plot(draw_opt)

    # avoid x title hitting labels
    plot.subplot_container.GetXaxis().SetLabelSize(plot.subplot_container.GetXaxis().GetLabelSize()*0.9)

    plot.set_logx(do_more_labels=False)
    plot.set_logy(do_more_labels=False)

    # Special case if data first object
    if data_first:
        # Redraw data on top
        data_hist = plot.contributions[0].obj.Clone()
        data_hist.SetLineWidth(entries[0][1].get('line_width', 1))
        data_hist.SetMarkerSize(entries[0][1].get('marker_size', 1))
        plot.main_pad.cd()
        data_draw_opt = "E1 X0 SAME"
        data_hist.Draw(data_draw_opt)

        # Create dummy graphs with the same styling to put into the legend
        # Using graphs we can get the correct endings on the TLegend entries (!)
        # Yes a massive faff for something so trivial
        dummy_gr = ROOT.TGraphErrors(1, array('d', [1]), array('d', [1]), array('d', [1]), array('d', [1]))
        dummy_hist = ROOT.TGraphErrors(1, array('d', [1]), array('d', [1]), array('d', [1]), array('d', [1]))
        # Add them to the legend and draw it
        dummy_conts = []  # stop premature deletion
        for i, entry in enumerate(entries):
            leg_draw_opt = "LE"
            # check if line_width > 0?
            if i == 0:
                if "X0" in data_draw_opt:
                    leg_draw_opt = "E"
                if data_hist.GetMarkerSize() > 0:
                    leg_draw_opt += "P"

            cont = Contribution(dummy_gr.Clone(), leg_draw_opt=leg_draw_opt, **entry[1])
            dummy_conts.append(cont)
            plot.legend.AddEntry(cont.obj, cont.label, cont.leg_draw_opt)

        plot.canvas.cd()
        plot.legend.Draw()

        # Do the subplot uncertainty shading
        data_no_errors = entries[0][0].Clone()
        cu.remove_th1_errors(data_no_errors)

        # Create hists for data with error region for ratio
        # Easiest way to get errors right is to do data (with 0 errors)
        # and divide by data (with errors), as if you had MC = data with 0 error

        # data_stat_ratio = data_no_errors.Clone()
        # data_stat_ratio.Divide(unfolded_hist_bin_stat_errors)
        # data_stat_ratio.SetFillStyle(3245)
        # data_stat_ratio.SetFillColor(self.plot_colours['unfolded_stat_colour'])
        # data_stat_ratio.SetLineWidth(0)
        # data_stat_ratio.SetMarkerSize(0)

        data_total_ratio = data_no_errors.Clone()
        # compare_bins(data_total_ratio, entries[0][0])
        data_total_ratio.Divide(entries[0][0])
        data_total_ratio.SetFillStyle(3254)
        data_total_ratio.SetFillColor(entries[0][1]['fill_color'])
        data_total_ratio.SetLineWidth(0)
        data_total_ratio.SetMarkerSize(0)

        # now draw the data error shaded area
        # this is a bit hacky - basically draw them on the ratio pad,
        # then redraw the existing hists & line to get them ontop
        # note that we use "same" for all - this is to keep the original axes
        # (we may want to rethink this later?)
        plot.subplot_pad.cd()
        data_draw_opt = "E2 SAME ]["  # need SAME otherwise axis get rescaled
        # data_stat_ratio.Draw(data_draw_opt)
        data_total_ratio.Draw(data_draw_opt)
        plot.subplot_container.Draw("SAME " + draw_opt)
        plot.subplot_line.Draw()
        plot.canvas.cd()

    plot.save(output_filename)


def _rebin_scale(hist, binning, normalise=False):
    if binning is not None:
        new_hist = hist.Rebin(len(binning)-1, hist.GetName()+"_rebin" + cu.get_unique_str(), binning)
    else:
        new_hist = hist.Clone(hist.GetName()+"_rebin_scale" + cu.get_unique_str())
    factor = 1.
    if normalise and new_hist.Integral() != 0:
        factor = new_hist.Integral()
    # Urgh this breaks .Integral() somehow
    new_hist.Scale(1./factor, "width")
    return new_hist


def tunfold_to_phsyical_bins(hist, bins, divide_by_bin_width=True):
    if hist.GetNbinsX() != len(bins)-1:
        raise ValueError("tunfold_to_phsyical_bins: bins not correct size (%d vs %d)" % (hist.GetNbinsX(), len(bins)-1))

    new_hist = ROOT.TH1D(hist.GetName()+"_relabel" + cu.get_unique_str(), "", len(bins)-1, bins)
    for ix in range(1, hist.GetNbinsX()+1):
        bin_width = bins[ix] - bins[ix-1]
        if not divide_by_bin_width:
            bin_width = 1.
        # print(ix, hist.GetBinError(ix), hist.GetBinError(ix) / bin_width)
        new_hist.SetBinContent(ix, hist.GetBinContent(ix) / bin_width)
        new_hist.SetBinError(ix, hist.GetBinError(ix) / bin_width)
    new_hist.SetEntries(hist.GetEntries())
    return new_hist


def do_dijet_pt_plots(workdir):

    data_tfile = cu.open_root_file(os.path.join(workdir, qgc.JETHT_ZB_FILENAME))
    mg_tfile = cu.open_root_file(os.path.join(workdir, qgc.QCD_FILENAME))
    py_tfile = cu.open_root_file(os.path.join(workdir, qgc.QCD_PYTHIA_ONLY_FILENAME))
    hpp_tfile = cu.open_root_file(os.path.join(workdir, qgc.QCD_HERWIG_FILENAME))

    for region_shortname, region_label in [("central", qgc.Dijet_CEN_LABEL), ("forward", qgc.Dijet_FWD_LABEL)]:

        # histname = "Dijet_QG_%s_tighter/jet_pt" % (region_shortname)  # fine equidistant binning
        histname = "Dijet_QG_Unfold_%s_tighter/hist_pt_reco_all" % (region_shortname) # tunfold binning

        lw = 2
        msize = 1.1
        data_line_width = lw
        mc_msize = 1E-3  # small enough to not be seen but not 0 - otherwise the horizontal lines on error bars don't get drawn

        mgpy_label = "MG5+Pythia8"
        hpp_label = "Herwig++"
        col_hpp = qgc.QCD_COLOURS[3]

        data_hist = cu.get_from_tfile(data_tfile, histname)
        mg_hist = cu.get_from_tfile(mg_tfile, histname)
        py_hist = cu.get_from_tfile(py_tfile, histname)
        hpp_hist = cu.get_from_tfile(hpp_tfile, histname)

        # For use with tunfold binning, which just has bin indices as x values instead of physical values
        all_pt_bins = np.append(qgc.PT_UNFOLD_DICT['underflow_reco'][:-1], qgc.PT_UNFOLD_DICT['signal_reco'])
        all_pt_bins = np.append(all_pt_bins, 8000)  # the overflow bin
        print(all_pt_bins)
        data_hist = tunfold_to_phsyical_bins(data_hist, all_pt_bins, divide_by_bin_width=True)
        mg_hist = tunfold_to_phsyical_bins(mg_hist, all_pt_bins, divide_by_bin_width=True)
        py_hist = tunfold_to_phsyical_bins(py_hist, all_pt_bins, divide_by_bin_width=True)
        hpp_hist = tunfold_to_phsyical_bins(hpp_hist, all_pt_bins, divide_by_bin_width=True)

        # Scale to data
        mg_hist.Scale(data_hist.Integral()/mg_hist.Integral())
        py_hist.Scale(data_hist.Integral()/py_hist.Integral())
        hpp_hist.Scale(data_hist.Integral()/hpp_hist.Integral())


        # custom_bins = np.concatenate([
        #     np.arange(50, 200, 10),
        #     np.arange(200, 300, 20),
        #     np.arange(300, 1000., 50),
        #     np.arange(1000, 2100., 100),
        #     # np.array([500, 600, 700, 800, 900, 1000], dtype='d')
        #     ])
        # print('dijet binning:', custom_bins)

        # can't use this, screws up nentries for some reason
        # data_hist = _rebin_scale(data_hist, binning=custom_bins, normalise=False)
        # mg_hist = _rebin_scale(mg_hist, binning=custom_bins, normalise=False)
        # py_hist = _rebin_scale(py_hist, binning=custom_bins, normalise=False)
        # hpp_hist = _rebin_scale(hpp_hist, binning=custom_bins, normalise=False)

        # hpp_hist.Scale(23.76 / )
        # hpp_hist.Scale(1E6)
        # for ih in range(1, 4):
            # print(hpp_hist.GetBinContent(ih) / data_hist.GetBinContent(ih))

        entries = [
            # DATA
            [
                data_hist,
                dict(line_color=qgc.JETHT_COLOUR, line_width=data_line_width, fill_color=qgc.JETHT_COLOUR,
                     marker_color=qgc.JETHT_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=msize*0.7,
                     label="Data")
            ],

            # MG5+PYTHIA8 MC
            [
                mg_hist,
                dict(line_color=qgc.QCD_COLOUR, line_width=lw, fill_color=qgc.QCD_COLOUR,
                     marker_color=qgc.QCD_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=mc_msize,
                     label=mgpy_label,
                     subplot=data_hist)
            ],

            # PYTHIA-only MC
            [
                py_hist,
                dict(line_color=qgc.QCD_COLOURS[2], line_width=lw, fill_color=qgc.QCD_COLOURS[2],
                     marker_color=qgc.QCD_COLOURS[2], marker_style=cu.Marker.get(qgc.QCD_MARKER), marker_size=mc_msize,
                     label="Pythia8",
                     subplot=data_hist)
            ],

            # HERWIG++
            [
                hpp_hist,
                dict(line_color=col_hpp, line_width=lw, fill_color=col_hpp,
                     marker_color=col_hpp, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=mc_msize,
                     label=hpp_label,
                     subplot=data_hist)
            ]
        ]

        radius, pus = cu.get_jet_config_from_dirname(workdir)
        jet_str = "AK%s PF %s" % (radius.upper(), pus.upper())
        title = "{jet_algo}\n{region_label}".format(jet_algo=jet_str,
                                                      region_label=region_label)

        do_jet_pt_plot(entries,
                       output_filename=os.path.join(workdir, "data_mc_jet_pt/Dijet_%s/jet_pt.%s" % (region_shortname, OUTPUT_FMT)),
                       rebin=1,
                       xlim=(30, 4000),
                       ylim=None,
                       title=title,
                       subplot_limits=(0, 2.5),
                       data_first=True,
                       normalise_hists=False)



def do_zpj_pt_plots(workdir):

    single_mu_tfile = cu.open_root_file(os.path.join(workdir, qgc.SINGLE_MU_FILENAME))
    mg_dy_tfile = cu.open_root_file(os.path.join(workdir, qgc.DY_FILENAME))
    hpp_dy_tfile = cu.open_root_file(os.path.join(workdir, qgc.DY_HERWIG_LOW_HIGH_PT_FILENAME))

    # histname = "ZPlusJets_QG/jet_pt"  # fine equidistant binning
    histname = "ZPlusJets_QG_Unfold/hist_pt_reco_all"  # tunfold binning

    lw = 2
    msize = 1.1
    data_line_width = lw
    mc_msize = 1E-3  # small enough to not be seen but not 0 - otherwise the horizontal lines on error bars don't get drawn

    mgpy_label = "MG5+Pythia8"
    hpp_label = "Herwig++"
    col_hpp = qgc.DY_COLOURS[2]

    data_hist = cu.get_from_tfile(single_mu_tfile, histname)
    mg_hist = cu.get_from_tfile(mg_dy_tfile, histname)
    hpp_hist = cu.get_from_tfile(hpp_dy_tfile, histname)

    # For use with tunfold binning, which just has bin indices as x values instead of physical values
    all_pt_bins = np.append(qgc.PT_UNFOLD_DICT['underflow_zpj_reco'][:-1], qgc.PT_UNFOLD_DICT['signal_zpj_reco'])
    all_pt_bins = np.append(all_pt_bins, 8000)  # the overflow bin
    print(all_pt_bins)
    data_hist = tunfold_to_phsyical_bins(data_hist, all_pt_bins, divide_by_bin_width=True)
    mg_hist = tunfold_to_phsyical_bins(mg_hist, all_pt_bins, divide_by_bin_width=True)
    hpp_hist = tunfold_to_phsyical_bins(hpp_hist, all_pt_bins, divide_by_bin_width=True)

    # Scale to data
    mg_hist.Scale(data_hist.Integral()/mg_hist.Integral())
    hpp_hist.Scale(data_hist.Integral()/hpp_hist.Integral())

    # custom_bins = np.concatenate([
    #     np.arange(50, 100, 10),
    #     np.arange(100, 200, 20),
    #     np.arange(200, 300, 25),
    #     np.arange(300, 500., 50),
    #     np.arange(500, 1100., 100),
    #     # np.array([500, 600, 700, 800, 900, 1000], dtype='d')
    #     ])

    # custom_bins = np.concatenate([
    #     np.arange(50, 200, 10),
    #     np.arange(200, 300, 20),
    #     np.arange(300, 500., 50),
    #     np.arange(500, 1100., 100),
    #     ])
    # print('zpj binning:', custom_bins)

    # data_hist = _rebin_scale(data_hist, custom_bins)
    # mg_hist = _rebin_scale(mg_hist, custom_bins)
    # hpp_hist = _rebin_scale(hpp_hist, custom_bins)

    entries = [
        # SINGLE MU DATA
        [
            data_hist,
            dict(line_color=qgc.SINGLE_MU_COLOUR, line_width=data_line_width, fill_color=qgc.SINGLE_MU_COLOUR,
                 marker_color=qgc.SINGLE_MU_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=msize*0.7,
                 label="Data")
        ],

        # MG5+PYTHIA8 MC
        [
            mg_hist,
            dict(line_color=qgc.DY_COLOUR, line_width=lw, fill_color=qgc.DY_COLOUR,
                 marker_color=qgc.DY_COLOUR, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=mc_msize,
                 label=mgpy_label,
                 subplot=data_hist)
        ],

        # HERWIG++
        [
            hpp_hist,
            dict(line_color=col_hpp, line_width=lw, fill_color=col_hpp,
                 marker_color=col_hpp, marker_style=cu.Marker.get(qgc.DY_MARKER), marker_size=mc_msize,
                 label=hpp_label,
                 subplot=data_hist)
        ]
    ]

    print("data_hist.Integral('width'):", data_hist.Integral("width"))

    radius, pus = cu.get_jet_config_from_dirname(workdir)
    jet_str = "AK%s PF %s" % (radius.upper(), pus.upper())
    title = "{jet_algo}\n{region_label}".format(jet_algo=jet_str,
                                                  region_label=qgc.ZpJ_LABEL)

    do_jet_pt_plot(entries,
                   output_filename=os.path.join(workdir, "data_mc_jet_pt/ZPlusJets/jet_pt.%s" % (OUTPUT_FMT)),
                   rebin=1,
                   xlim=(30, 1000),
                   ylim=None,
                   title=title,
                   data_first=True,
                   subplot_limits=(0, 2),
                   normalise_hists=False)


if __name__ == "__main__":
    parser = qgc.get_parser()
    args = parser.parse_args()

    for workdir in args.workdirs:
        do_dijet_pt_plots(workdir)
        do_zpj_pt_plots(workdir)

    sys.exit(0)