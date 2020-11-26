#!/usr/bin/env python

"""Create markdown snippet of binning"""

import qg_common as qgc


if __name__ == "__main__":
    contents = []

    # Add lambda bins
    grooming_headers = [
        ["ungroomed", "## Ungroomed binning"],
        ["groomed", "## Groomed binning"],
    ]

    charged_headers = [
        ["", "### Charged+neutral lambda variables"],
        ["_charged", "### Charged-only lambda variables"]
    ]

    variables = [
        ['jet_puppiMultiplicity', "Multiplicity (kappa=0, beta=0)"],
        ['jet_pTD', "pTD (kappa=2, beta=0)"],
        ['jet_LHA', "LHA (kappa=1, beta=0.5)"],
        ['jet_width', "Width (kappa=1, beta=1)"],
        ['jet_thrust', "Thrust (kappa=1, beta=2)"],
    ]

    for grooming_str, grooming_title in grooming_headers:
        contents.append(grooming_title + "\n")
        for charged_str, charged_title in charged_headers:
            contents.append(charged_title + "\n")
            for var_str, var_title in variables:
                bins = list(qgc.VAR_UNFOLD_DICT[grooming_str][var_str + charged_str]['gen'])
                contents.append("- {var_title}: {bins}\n".format(var_title=var_title, bins=bins))

    # Add pT bins
    contents.append("## pT bins\n")
    contents.append("- Dijet: {}\n".format(list(qgc.PT_UNFOLD_DICT['signal_gen'])))
    contents.append("- Z+jet: {}\n".format(list(qgc.PT_UNFOLD_DICT['signal_zpj_gen'])))

    print('\n'.join(contents))
