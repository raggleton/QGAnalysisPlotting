#!/usr/bin/env python

"""Create config for rebinning slides

Use with beamer-plot-slides
"""


import json
from itertools import product

import sys
sys.path.append("..")
import qg_common as qgc

if __name__ == "__main__":
    pt_regions = [
        {
            "append": "_lowPt",
            "title": "30 < p_{T}^{\\text{Reco}} < 100~\\text{GeV}",
        },
        {
            "append": "_midPt",
            "title": "100 < p_{T}^{\\text{Reco}} < 250~\\text{GeV}",
        },
        {
            "append": "_highPt",
            "title": "p_{T}^{\\text{Reco}} > 250~\\text{GeV}",
        },
    ]

    current_binning_dir = "/Volumes/Extreme SSD/Projects/QGAnalysis/workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_fixLambda/rebinning_uhh2.AnalysisModuleRunner.MC.MC_QCD/"
    optimal_binning_dir = "/Volumes/Extreme SSD/Projects/QGAnalysis/workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_fixLambda/lambda_binning/"
    # this is using new binning but ungroomed for all
    rebinned_binning_dir = "/Volumes/Extreme SSD/Projects/QGAnalysis/workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_fixLambda/rebin_ungroomed_optimal/"

    contents = []
    for angle, pt_region_dict in product(qgc.COMMON_VARS[:], pt_regions):
        pt_append = pt_region_dict['append']
        extra = "for_midPt_" if pt_append != "_midPt" else ""
        this_slide = {
            "title": f"${angle.mathmode}:~{pt_region_dict['title']}$",
            "plots": [
                # [f"{current_binning_dir}/Dijet_QG_central_tighter/{angle.var}{pt_append}_rebinned_1d_reco_gen", "Using existing binning (central)"],
                # [f"{current_binning_dir}/Dijet_QG_forward_tighter/{angle.var}{pt_append}_rebinned_1d_reco_gen", "Using existing binning (forward)"],
                [f"{current_binning_dir}/Dijet_QG_central_tighter+Dijet_QG_forward_tighter/{angle.var}{pt_append}_rebinned_1d_reco_gen", "Using existing binning (central+forward)"],
                [f"{optimal_binning_dir}/Dijet_QG_central_tighter+Dijet_QG_forward_tighter/{angle.var}{pt_append}_rebinned_{extra}1d_reco_gen", "Optimal binning using new sample (central+forward)"],
                [f"{optimal_binning_dir}/Dijet_QG_central_tighter+Dijet_QG_forward_tighter/{angle.var}{pt_append}_orig_1d_reco_gen", "Original"],
                # [f"{current_binning_dir}/Dijet_QG_central_tighter/{angle.var}{pt_append}_rebinned_migration_summary", ""],
                # [f"{current_binning_dir}/Dijet_QG_forward_tighter/{angle.var}{pt_append}_rebinned_migration_summary", ""],
                [f"{current_binning_dir}/Dijet_QG_central_tighter+Dijet_QG_forward_tighter/{angle.var}{pt_append}_rebinned_migration_summary", ""],
                [f"{optimal_binning_dir}/Dijet_QG_central_tighter+Dijet_QG_forward_tighter/{angle.var}{pt_append}_rebinned_{extra}migration_summary", ""],
            ],
        }
        contents.append(this_slide)
        
        # add groomed plots
        extra = "for_midPt_"
        this_slide = {
            "title": f"${angle.mathmode}:~{pt_region_dict['title']}$ (groomed)",
            "plots": [
                # [f"{current_binning_dir}/Dijet_QG_central_tighter/{angle.var}{pt_append}_rebinned_1d_reco_gen", "Using existing binning (central)"],
                # [f"{current_binning_dir}/Dijet_QG_forward_tighter/{angle.var}{pt_append}_rebinned_1d_reco_gen", "Using existing binning (forward)"],
                [f"{current_binning_dir}/Dijet_QG_central_tighter_groomed+Dijet_QG_forward_tighter_groomed/{angle.var}{pt_append}_rebinned_1d_reco_gen", "Using existing binning (central+forward)"],
                [f"{rebinned_binning_dir}/Dijet_QG_central_tighter_groomed+Dijet_QG_forward_tighter_groomed/{angle.var}{pt_append}_rebinned_{extra}1d_reco_gen", "Optimal binning using new sample (central+forward, ungroomed)"],
                [f"{optimal_binning_dir}/Dijet_QG_central_tighter_groomed+Dijet_QG_forward_tighter_groomed/{angle.var}{pt_append}_orig_1d_reco_gen", "Original"],
                # [f"{current_binning_dir}/Dijet_QG_central_tighter/{angle.var}{pt_append}_rebinned_migration_summary", ""],
                # [f"{current_binning_dir}/Dijet_QG_forward_tighter/{angle.var}{pt_append}_rebinned_migration_summary", ""],
                [f"{current_binning_dir}/Dijet_QG_central_tighter_groomed+Dijet_QG_forward_tighter_groomed/{angle.var}{pt_append}_rebinned_migration_summary", ""],
                [f"{rebinned_binning_dir}/Dijet_QG_central_tighter_groomed+Dijet_QG_forward_tighter_groomed/{angle.var}{pt_append}_rebinned_{extra}migration_summary", ""],
            ],
        }
        contents.append(this_slide)

    with open("configuration.json") as f:
        template = f.read()

    contents_json = json.dumps(contents, indent=4)
    this_contents = template.replace("@CONTENTS", contents_json)

    with open("rebinning.json", "w") as f:
        f.write(this_contents)
