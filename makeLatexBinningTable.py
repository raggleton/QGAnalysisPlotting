#!/usr/bin/env python


"""
Create latex snippets for binning tables etc to put in AN, etc
"""


from __future__ import print_function

import qg_common as qgc

# do this manually because ROOT uses a weird latex
aliases = {
    "jet_puppiMultiplicity": r"\multi",
    "jet_LHA": r"\lha",
    "jet_pTD": r"\ptd",
    "jet_width": r"\width",
    "jet_thrust": r"\thrust",
}

def make_values_strings(values):
    return ["%g" % x for x in values]


def print_variable_table(var_dict):
    for angle_name, angle_binning_dict in var_dict.items():
        if "_charged" in angle_name:
            continue
        angle = [a for a in qgc.COMMON_VARS if a.var == angle_name]
        if len(angle) == 0:
            raise RuntimeError("Cannot find angle %s" % angle_name)
        if len(angle) > 1:
            raise RuntimeError("Found too many angles %s" % angle_name)

        # do charged+neutral
        angle = angle[0]
        angle_name_str = angle.name
        if "ptd" in angle_name.lower():
            angle_name_str = "$%s$" % angle_name_str
        # gen binning
        gen_latex_line = "\\multirow{{2}}{{*}}{name} & Generator & {binning} \\\\"
        this_gen_entry = {
            "name": "{%s ($%s$)}" % (angle_name_str, aliases[angle_name]),
            "binning": ", ".join(make_values_strings(angle_binning_dict['gen'])),
        }
        print(gen_latex_line.format(**this_gen_entry))

        # reco binning
        reco_latex_line = "& Detector & {binning} \\\\ [\\cmsTabSkip]"
        this_reco_entry = {
            "binning": ", ".join(make_values_strings(angle_binning_dict['reco'])),
        }
        print(reco_latex_line.format(**this_reco_entry))

        # do charged-only
        charged_angle_name = angle_name + "_charged"
        v_charged = var_dict[charged_angle_name]
        # gen
        this_gen_entry = {
            "name": "{Charged %s ($%s$)}" % (angle_name_str, aliases[angle_name]),
            "binning": ", ".join(make_values_strings(v_charged['gen'])),
        }
        print(gen_latex_line.format(**this_gen_entry))
        # reco
        this_reco_entry = {
            "binning": ", ".join(make_values_strings(v_charged['reco'])),
        }
        print(reco_latex_line.format(**this_reco_entry))


def print_pt_lists():
    pt_bins_gen = list(qgc.PT_UNFOLD_DICT['underflow_gen'])+list(qgc.PT_UNFOLD_DICT['signal_gen'])[1:]
    print(", ".join(['%g' % x for x in pt_bins_gen]) + "\\GeV")
    pt_bins_reco = list(qgc.PT_UNFOLD_DICT['underflow_reco'])+list(qgc.PT_UNFOLD_DICT['signal_reco'])[1:]
    print(", ".join(['%g' % x for x in pt_bins_reco]) + "\\GeV")


if __name__ == "__main__":
    print_variable_table(qgc.VAR_UNFOLD_DICT)
    print_pt_lists()
