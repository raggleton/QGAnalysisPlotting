#!/usr/bin/env python

"""Print 2D flavour correlation plots. Iterates over every hist in every dir"""

import ROOT
from MyStyle import My_Style
My_Style.cd()
import os
from uuid import uuid1

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
ROOT.gStyle.SetPaintTextFormat('.2e')

My_Style.SetPadLeftMargin(0.1)
My_Style.SetPadRightMargin(0.12)

# Control output format
OUTPUT_FMT = "pdf"


def get_list_of_obj(thing):
    key_list = thing.GetListOfKeys()
    return [x.GetName() for x in key_list]


def extract_pt_eta_from_name(name):
    parts = name.split("_")
    pt_edges, eta_edges = [], []
    region = parts[0]
    for item in parts:
        if item.startswith("Pt"):
            pt_edges = item.replace("Pt", "").split("to")
        elif item.startswith("Eta"):
            eta_edges = item.replace("Eta", "").split("to")
    return region, pt_edges, eta_edges


def print_plots(directory, output_dir, title=""):
    plot_names = get_list_of_obj(directory)
    cu.check_dir_exists_create(output_dir)
    for pname in plot_names:
        hist = directory.Get(pname)
        c = ROOT.TCanvas(str(uuid1()), "", 600, 600)
        c.SetLogz()
        hist.SetTitle(title)
        hist.Draw("COLZ TEXT89")
        c.SaveAs(os.path.join(output_dir, pname + "." + OUTPUT_FMT))

        c.Clear()
        hist_normed = cu.make_normalised_TH2(hist, "Y", False)
        hist_normed.SetMinimum(1E-4)
        hist_normed.Draw("COLZ TEXT89")
        c.SaveAs(os.path.join(output_dir, pname + "_normed." + OUTPUT_FMT))


def make_plots(input_filename, output_dir):
    input_file = cu.open_root_file(input_filename)
    forbidden = ['SFrame', 'cf_metfilters_raw', 'cf_metfilters']
    directories = [d for d in get_list_of_obj(input_file) if d not in forbidden]
    print(directories)
    for d in directories:
        region, pt_edges, eta_edges = extract_pt_eta_from_name(d)
        title = "%s, %s < p_{T} < %s GeV, %s < |#eta| < %s" % (region, *pt_edges, *eta_edges)
        print_plots(input_file.Get(d), os.path.join(output_dir, d), title)


if __name__ == "__main__":
    parser = qgc.get_parser()
    args = parser.parse_args()

    for workdir in args.workdirs:
        print(workdir)

        # if os.path.isfile(os.path.join(workdir, qgc.QCD_FILENAME)):
        #     make_plots(os.path.join(workdir, qgc.QCD_FILENAME), os.path.join(workdir, "QCD"))

        if os.path.isfile(os.path.join(workdir, qgc.DY_FILENAME)):
            make_plots(os.path.join(workdir, qgc.DY_FILENAME), os.path.join(workdir, "DY"))


    sys.exit(0)
