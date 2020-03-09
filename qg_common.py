"""Common vars for doing QG plots"""

from __future__ import print_function, division

import numpy as np
import argparse
from collections import namedtuple

# Use ints for ROOT colours as importing ROOT seems OTT

# Purple & blue look ~ same to red-green colourblind people
# (V.Dark blue & light purple do look a bit diff)
# Red and blue does look different
# DY_COLOUR =  880 # kViolet
DY_COLOUR =  628 # kViolet
QCD_COLOUR = 867 # kAzure+7, for backwards compatibility
QCD_CEN_COLOUR = 867 #
QCD_FWD_COLOUR = 868 #

# SINGLE_MU_COLOUR =  884 # kViolet+3
SINGLE_MU_COLOUR =  634
JETHT_COLOUR = 600 # kBlue
JETHT_CEN_COLOUR = 601 # kBlue
JETHT_FWD_COLOUR = 600 # kBlue
ZB_COLOUR = 416+1 # kGreen+1

# You should use markers as well for colourblindness
DY_MARKER = 20
DY_MARKER = 'circle'

QCD_MARKER = 22
QCD_MARKER = 'triangleUp'

QCD_CEN_MARKER = 22
QCD_CEN_MARKER = 'triangleUp'

QCD_FWD_MARKER = 23
QCD_FWD_MARKER = 'triangleDown'

# When comparing e.g. PU bins
# FIXME: update for colourblindness
# DY_COLOURS = [880, 881, 884, 873]
DY_COLOURS = [880, 881, 797, 902]
DY_COLOURS *= 3
QCD_COLOURS = [867, 600, 853, 882]
QCD_COLOURS *= 3  # for lots

PT_BINS = [(80, 100), (100, 200), (400, 500), (1000, 2000), (80, 2000)]
PT_BINS = [(80, 100), (100, 120), (200, 250), (400, 500), (1000, 2000), (80, 2000)]
THEORY_PT_BINS = [(100, 200), (400, 500), (1000, 2000), (80, 2000)]
THEORY_PT_BINS = [(100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 800), (800, 1000), (1000, 2000), (80, 2000)]
# THEORY_PT_BINS = [(20, 40), (40, 60), (60, 80), (80, 100), (100, 120),
THEORY_PT_BINS = [(20, 25), (40, 50), (50, 55), (60, 70), (70, 80), (80, 90), (90, 100), (100, 120),
                  (120, 160), (160, 200), (200, 220), (280, 300), (300, 350), (400, 410),
                  (400, 500), (500, 600), (500, 550), (600, 800), (800, 1000),
                  (1000, 1400), (1400, 2000)]
# PT_BINS = [
# (0, 29), (29,38), (38,50), (50,65), (65,88), (88,120), (120,150), (150,186), (186,254),
# (254,326), (326,408), (408,481), (481,614), (614,800), (800,1000), (1000,1300),
# (1300,1700), (1700,2200), (2200,3000), (3000,4000), (4000,5000), (5000,10000),
# ]

PT_BINS = [
(15, 23), (23, 30), (30, 38), (38, 50), (50, 65), (65, 88), (88, 120), (120, 150),
(150, 186), (186, 254), (254, 326), (326, 408), (408, 481), (481, 614), (614, 800),
(800, 1000), (1000, 2000), (2000, 6500)
]

PT_BINS_ZPJ = [
(15, 23), (23, 30), (30, 38), (38, 50), (50, 65), (65, 88), (88, 120), (120, 150),
(150, 186), (186, 254), (254, 326), (326, 408), (408, 481), (481, 614), (614, 800),
(800, 6500)
]

# PT_BINS = [
# (29,38), (38,50), (50,65), (65,88), (88,114), (114,145), (145,180), (180,254),
# (254,318), (318,400), (400,500), (500,625), (625,800), (800,1000), (1000,1300),
# (1300,1700), (1700,2200), (2200,3000), (3000,4000), (4000,5000), (5000,10000),
# ]
#

PT_BINS = [
(50, 65), (65, 88), (88, 120), (120, 150), (150, 186),
(186, 254), (254, 326), (326, 408), (408, 481), (481, 614),
(614, 800), (800, 1000), (1000, 1500), (1500, 2000), (2000, 6500)
]

PT_BINS_ZPJ = [
(50, 65), (65, 88), (88, 120), (120, 150), (150, 186),
(186, 254), (254, 326), (326, 408), (408, 481), (481, 614),
(614, 800), (800, 6500)
]

PT_BINS_INC_UFLOW = [(30, 38), (38, 50)]
PT_BINS_INC_UFLOW += PT_BINS


def construct_fine_binning(coarse_bin_edges):
    fine_bin_edges = []
    for x, y in zip(coarse_bin_edges[:-1], coarse_bin_edges[1:]):
        fine_bin_edges.append(x)
        fine_bin_edges.append(0.5*(x+y))
    fine_bin_edges.append(coarse_bin_edges[-1])
    return np.array(fine_bin_edges, dtype='d')


def construct_all_fine_binnings(var_dict):
    # Construct fine binning from splitting coarser bins
    # Coarser bins dervied from determine_lambda_binning.py
    for angle_name, angle_dict in var_dict.items():
        if angle_dict['reco'] is None:
            angle_dict['reco'] = construct_fine_binning(angle_dict['gen'])
        var_dict[angle_name] = angle_dict


Angle = namedtuple("Angle", ['var', 'kappa', 'beta', 'name', "lambda_str", "colour"])
COMMON_VARS_WITH_FLAV = [
    # Angle("jet_multiplicity", 0, 0, "Multiplicity", "#lambda_{0}^{0}", 2),
    Angle("jet_puppiMultiplicity", 0, 0, "Multiplicity", "#lambda_{0}^{0}", 2),
    Angle('jet_pTD', 2, 0, "(p_{T}^{D})^{2}", "#lambda_{0}^{2}", 418),
    Angle('jet_LHA', 1, 0.5, "LHA", "#lambda_{0.5}^{1}", 600),
    Angle('jet_width', 1, 1, "Width", "#lambda_{1}^{1}", 861),
    Angle('jet_thrust', 1, 2, "Thrust", "#lambda_{2}^{1}", 617),
    # charged-only constit
    # Angle("jet_multiplicity_charged", 0, 0, "Multiplicity (charged)", "#lambda_{0}^{0}", 2),
    Angle("jet_puppiMultiplicity_charged", 0, 0, "Multiplicity (charged-only)", "#lambda_{0}^{0}", 2),
    Angle('jet_pTD_charged', 2, 0, "(p_{T}^{D})^{2} (charged-only)", "#lambda_{0}^{2}", 418),
    Angle('jet_LHA_charged', 1, 0.5, "LHA (charged-only)", "#lambda_{0.5}^{1}", 600),
    Angle('jet_width_charged', 1, 1, "Width (charged-only)", "#lambda_{1}^{1}", 861),
    Angle('jet_thrust_charged', 1, 2, "Thrust (charged-only)", "#lambda_{2}^{1}", 617),
    # Angle('jet_multiplicity', 0, 0, "Multiplicity", "#lambda_{0}^{0}", 2),  # leave last as we don't want it in our RIVET numbering
    Angle('jet_flavour', 0, 0, "Flavour", "PDGID", 7),
    Angle("jet_genParton_flavour", 0, 0, "Flavour", "PDGID", 8)
]

COMMON_VARS = COMMON_VARS_WITH_FLAV[:-2]

ANGLE_REBIN_DICT = {
    "jet_puppiMultiplicity": [0.0, 6.0, 9.0, 12.0, 18.0, 150.0],
    'jet_multiplicity': [0.0, 6.0, 9.0, 12.0, 18.0, 150.0],
    'jet_pTD': [0.0, 0.1, 0.13, 0.17, 0.23, 0.33, 0.51, 0.85, 1.0],
    'jet_LHA': [0.0, 0.29, 0.37, 0.44, 0.5, 0.56, 0.62, 0.68, 0.75, 1.0],
    'jet_width': [0.0, 0.12, 0.18, 0.24, 0.3, 0.36, 0.43, 0.51, 1.0],
    'jet_thrust': [0.0, 0.04, 0.08, 0.12, 0.17, 0.24, 0.33, 1.0],
}

PT_UNFOLD_DICT = {
    # "signal_gen": np.array([50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 1000, 6500], dtype='d'),
    "signal_gen": np.array([50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 1000, 4000], dtype='d'),
    # "signal_gen": np.array([50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 1000, 2000], dtype='d'),
    "underflow_gen": np.array([15, 30, 38, 50], dtype='d'),
    "signal_zpj_gen": np.array([50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 2000], dtype='d'),
    "underflow_zpj_gen": np.array([15, 30, 38, 50], dtype='d'),
}


for pt_name in list(PT_UNFOLD_DICT.keys()):
    new_name = pt_name.replace("_gen", "_reco")
    PT_UNFOLD_DICT[new_name] = construct_fine_binning(PT_UNFOLD_DICT[pt_name])

# remove the first two bins, 15, 22.5
PT_UNFOLD_DICT['underflow_reco'] = PT_UNFOLD_DICT['underflow_reco'][2:]
PT_UNFOLD_DICT['underflow_zpj_reco'] = PT_UNFOLD_DICT['underflow_zpj_reco'][2:]


VAR_UNFOLD_DICT = {
    'jet_multiplicity': {
        'gen': np.array([1, 5, 10, 13, 19, 25, 35, 50, 75, 100, 150], dtype='d'),
        'reco': None  # reco gets set later
    },
     'jet_puppiMultiplicity': {
        'gen': np.array([1, 5, 10, 13, 19, 25, 35, 50, 75, 100, 150], dtype='d'),
         'reco': None
     },
     'jet_pTD': {
        'gen': np.array([0.0, 0.09, 0.12, 0.16, 0.21, 0.29, 0.43, 0.7, 1.0], dtype='d'),
         'reco': None
     },
     'jet_LHA': {
        'gen': np.array([0.0, 0.14, 0.29, 0.37, 0.44, 0.5, 0.56, 0.62, 0.68, 0.75, 1.0], dtype='d'),
         'reco': None
     },
     'jet_width': {
        'gen': np.array([0.0, 0.11, 0.17, 0.23, 0.29, 0.35, 0.42, 0.6, 1.0], dtype='d'),
         'reco': None
     },
     'jet_thrust': {
        'gen': np.array([0.0, 0.04, 0.08, 0.12, 0.17, 0.24, 0.33, 0.66, 1.0], dtype='d'),
        'reco': None
    }
}

# target 0.5 values
VAR_UNFOLD_DICT = {
    'jet_puppiMultiplicity': {
        'gen': np.array([0, 9, 15, 22, 35, 50, 75, 100, 150], dtype='d'),
        'reco': None
    },
    'jet_pTD': {
        'gen': np.array([0.0, 0.09, 0.14, 0.25, 1.0], dtype='d'),
        'reco': None
    },
    'jet_LHA': {
        'gen': np.array([0.0, 0.12, 0.25, 0.33, 0.39, 0.45, 0.51, 0.57, 0.62, 0.66, 0.7, 0.76, 1.0], dtype='d'),
        'reco': None
    },
    'jet_width': {
        'gen': np.array([0.0, 0.09, 0.14, 0.19, 0.25, 0.31, 0.37, 0.43, 0.48, 0.54, 1.0], dtype='d'),
        'reco': None
    },
    'jet_thrust': {
        'gen': np.array([0.0, 0.03, 0.06, 0.1, 0.15, 0.2, 0.25, 0.31, 0.47, 1.0], dtype='d'),
        'reco': None
    },

    'jet_puppiMultiplicity_charged': {
        'gen': np.array([0, 9, 15, 22, 35, 50, 75, 100, 150], dtype='d'),
        'reco': None
    },
    'jet_pTD_charged': {
        'gen': np.array([0.0, 0.09, 0.12, 0.15, 0.19, 0.24, 0.31, 0.4, 0.53, 0.73, 1.0], dtype='d'),
        'reco': None
    },
    'jet_LHA_charged': {
        'gen': np.array([0.0, 0.1, 0.2, 0.26, 0.32, 0.37, 0.42, 0.47, 0.52, 0.57, 0.62, 0.66, 0.7, 0.75, 1.0], dtype='d'),
        'reco': None
    },
    'jet_width_charged': {
        'gen': np.array([0.0, 0.06, 0.09, 0.13, 0.17, 0.21, 0.26, 0.31, 0.36, 0.41, 0.46, 0.52, 0.6, 1.0], dtype='d'),
        'reco': None
    },
    'jet_thrust_charged': {
        'gen': np.array([0.0, 0.02, 0.04, 0.07, 0.11, 0.16, 0.22, 0.28, 0.37, 1.0], dtype='d'),
        'reco': None
    },
}

# WTA axis, no dijet eta split
VAR_UNFOLD_DICT_TARGET0p5 = {
    'jet_puppiMultiplicity': {
        'gen': np.array([0, 9, 15, 22, 35, 50, 75, 100, 150], dtype='d'),
        'reco': None
    },
    'jet_pTD': {
        'gen': np.array([0.0, 0.09, 0.14, 0.25, 1.0], dtype='d'),
        'reco': None
    },
    'jet_LHA': {
        'gen': np.array([0.0, 0.16, 0.23, 0.3, 0.36, 0.42, 0.49, 0.56, 0.63, 0.72, 1.0], dtype='d'),
        'reco': None
    },
    'jet_width': {
        'gen': np.array([0.0, 0.08, 0.13, 0.18, 0.24, 0.31, 0.39, 0.47, 0.58, 1.0], dtype='d'),
        'reco': None
    },
    'jet_thrust': {
        'gen': np.array([0.0, 0.04, 0.07, 0.12, 0.19, 0.27, 0.37, 0.52, 1.0], dtype='d'),
        'reco': None
    },

    'jet_puppiMultiplicity_charged': {
        'gen': np.array([0, 9, 15, 22, 35, 50, 75, 100, 150], dtype='d'),
        'reco': None
    },
    'jet_pTD_charged': {
        'gen': np.array([0.0, 0.09, 0.12, 0.15, 0.19, 0.24, 0.31, 0.4, 0.53, 0.73, 1.0], dtype='d'),
        'reco': None
    },
    'jet_LHA_charged': {
        'gen': np.array([0.0, 0.09, 0.14, 0.19, 0.24, 0.3, 0.36, 0.43, 0.51, 0.6, 0.71, 1.0], dtype='d'),
        'reco': None
    },
    'jet_width_charged': {
        'gen': np.array([0.0, 0.04, 0.07, 0.1, 0.14, 0.18, 0.23, 0.29, 0.37, 0.46, 0.58, 0.99, 1.0], dtype='d'),
        'reco': None
    },
    'jet_thrust_charged': {
        'gen': np.array([0.0, 0.02, 0.04, 0.07, 0.11, 0.17, 0.25, 0.36, 0.52, 1.0], dtype='d'),
        'reco': None
    },
}

# target 0.75 values
# VAR_UNFOLD_DICT = {

#     "jet_LHA": {
#         'gen': np.array([0.0, 0.35, 0.5, 0.64, 1.0], dtype='d'),
#         'reco': None
#     },
#     "jet_LHA_charged": {
#         'gen': np.array([0.0, 0.29, 0.41, 0.53, 0.66, 1.0], dtype='d'),
#         'reco': None
#     },
#     "jet_puppiMultiplicity": {
#         'gen': np.array([0, 22, 44, 150], dtype='d'),
#         'reco': None
#     },
#     "jet_puppiMultiplicity_charged": {
#         'gen': np.array([0, 16, 32, 150], dtype='d'),
#         'reco': None
#     },
#     "jet_pTD": {
#         'gen': np.array([0.0, 0.15, 1.0], dtype='d'),
#         'reco': None
#     },
#     "jet_pTD_charged": {
#         'gen': np.array([0.0, 0.15, 0.26, 0.6, 1.0], dtype='d'),
#         'reco': None
#     },
#     "jet_thrust": {
#         'gen': np.array([0.0, 0.06, 0.17, 0.31, 1.0], dtype='d'),
#         'reco': None
#     },
#     "jet_thrust_charged": {
#         'gen': np.array([0.0, 0.05, 0.15, 0.45, 1.0], dtype='d'),
#         'reco': None
#     },
#     "jet_width": {
#         'gen': np.array([0.0, 0.16, 0.3, 0.45, 1.0], dtype='d'),
#         'reco': None
#     },
#     "jet_width_charged": {
#         'gen': np.array([0.0, 0.11, 0.2, 0.33, 0.49, 1.0], dtype='d'),
#         'reco': None
#     },
# }

# Target 0.6, using mainly forward ungroomed jets
VAR_UNFOLD_DICT_TARGET0p6 = {
    'jet_puppiMultiplicity': {
        'gen': np.array([0.0, 11.0, 18.0, 25.0, 150.0 ], dtype='d'),
        'reco': None
    },
    'jet_pTD': {
        'gen': np.array([0.0, 0.1, 0.17, 0.4, 1.0], dtype='d'),
        'reco': None
    },
    'jet_LHA': {
        'gen': np.array([0.0, 0.19, 0.29, 0.38, 0.49, 0.6, 0.75, 1.0], dtype='d'),
        'reco': None
    },
    'jet_width': {
        'gen': np.array([0.0, 0.12, 0.2, 0.32, 0.45, 0.64, 1.0], dtype='d'),
        'reco': None
    },
    'jet_thrust': {
        'gen': np.array([0.0, 0.06, 0.16, 0.3, 0.51, 1.0], dtype='d'),
        'reco': None
    },

    'jet_puppiMultiplicity_charged': {
        'gen': np.array([0.0, 4.0, 6.0, 9.0, 12.0, 16.0, 23.0, 150.0], dtype='d'),
        'reco': None
    },
    'jet_pTD_charged': {
        'gen': np.array([0.0, 0.09, 0.12, 0.16, 0.21, 0.29, 0.41, 0.6, 1.0], dtype='d'),
        'reco': None
    },
    'jet_LHA_charged': {
        'gen': np.array([0.0, 0.11, 0.17, 0.24, 0.32, 0.42, 0.54, 0.69, 1.0], dtype='d'),
        'reco': None
    },
    'jet_width_charged': {
        'gen': np.array([0.0, 0.04, 0.07, 0.11, 0.16, 0.22, 0.3, 0.41, 0.57, 1.0], dtype='d'),
        'reco': None
    },
    'jet_thrust_charged': {
        'gen': np.array([0.0, 0.01, 0.02, 0.04, 0.07, 0.11, 0.18, 0.29, 0.48, 1.0], dtype='d'),
        'reco': None
    },
}

# new version of WTA, cen+fwd dijet, target 0.5
VAR_UNFOLD_DICT_TARGET0p5 = {
    'jet_puppiMultiplicity': {
        # 'gen': np.array([0.0, 10, 15, 20, 27, 50, 75, 100, 150], dtype='d'),
        'gen': np.array([0.0, 10, 15, 20, 30, 50, 75, 100, 150], dtype='d'),
        'reco': None
    },
    'jet_pTD': {
        # 'gen': np.array([0.0, 0.07, 0.1, 0.15, 0.24, 0.45, 1.0], dtype='d'),
        'gen': np.array([0.0, 0.07, 0.1, 0.15, 0.24, 1.0], dtype='d'),
        'reco': None
    },
    'jet_LHA': {
        # 'gen': np.array([0.0, 0.14, 0.22, 0.29, 0.35, 0.42, 0.49, 0.56, 0.64, 0.75, 1.0], dtype='d'),
        # 'gen': np.array([0.0, 0.14, 0.22, 0.29, 0.35, 0.42, 0.49, 0.56, 0.64, 1.0], dtype='d'),
        'gen': np.array([0.0, 0.14, 0.21, 0.28, 0.34, 0.4, 0.47, 0.54, 0.61, 1.0], dtype='d'),
        'reco': None
    },
    'jet_width': {
        # 'gen': np.array([0.0, 0.09, 0.145, 0.205, 0.28, 0.36, 0.445, 0.545, 1.0 ], dtype='d'),
        'gen': np.array([ 0.0, 0.09, 0.145, 0.205, 0.28, 0.355, 0.43, 0.515, 1.0], dtype='d'),
        'reco': None
    },
    'jet_thrust': {
        # 'gen': np.array([0.0, 0.04, 0.08, 0.145, 0.225, 0.32, 0.445, 1.0], dtype='d'),
        'gen': np.array([0.0, 0.04, 0.08, 0.14, 0.195, 0.25, 1.0], dtype='d'),
        'reco': None
    },

    'jet_puppiMultiplicity_charged': {
        # 'gen': np.array([0.0, 2.0, 3.0, 4.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 18.0, 21.0, 25.0, 32.0, 91, 150.0], dtype='d'),
        'gen': np.array([0.0, 3.0, 5.0, 7.0, 10.0, 13.0, 16.0, 20.0, 30, 50, 75, 100, 150.0], dtype='d'),
        'reco': None
    },
    'jet_pTD_charged': {
        # 'gen': np.array([0.0, 0.08, 0.1, 0.12, 0.14, 0.17, 0.2, 0.24, 0.29, 0.34, 0.4, 0.47, 0.55, 0.65, 0.76, 0.89, 1.0], dtype='d'),
        'gen': np.array([0.0, 0.07, 0.09, 0.11, 0.14, 0.18, 0.23, 0.3, 0.39, 0.51, 0.64, 1.0], dtype='d'),
        'reco': None
    },
    'jet_LHA_charged': {
        # 'gen': np.array([0.0, 0.05, 0.09, 0.12, 0.15, 0.18, 0.22, 0.26, 0.31, 0.36, 0.42, 0.49, 0.56, 0.64, 0.75, 1.0], dtype='d'),
        # 'gen': np.array([0.0, 0.05, 0.09, 0.12, 0.15, 0.18, 0.22, 0.26, 0.31, 0.36, 0.42, 0.49, 0.56, 0.64, 1.0], dtype='d'),
        'gen': np.array([0.0, 0.05, 0.09, 0.14, 0.19, 0.25, 0.32, 0.39, 0.47, 0.55, 0.63, 0.73, 1.0], dtype='d'),
        'reco': None
    },
    'jet_width_charged': {
        # 'gen': np.array([0.0, 0.015, 0.025, 0.035, 0.05, 0.065, 0.085, 0.105, 0.13, 0.16, 0.195, 0.235, 0.28, 0.335, 0.4, 0.48, 0.585, 1.0], dtype='d'),
        'gen': np.array([0.0, 0.01, 0.025, 0.045, 0.075, 0.11, 0.155, 0.21, 0.275, 0.34, 0.415, 0.495, 0.595, 1.0], dtype='d'),
        'reco': None
    },
    'jet_thrust_charged': {
        # 'gen': np.array([0.0, 0.005, 0.01, 0.015, 0.025, 0.035, 0.05, 0.065, 0.085, 0.115, 0.15, 0.2, 0.265, 0.355, 0.48, 0.725, 1.], dtype='d'),
        'gen': np.array([0.0, 0.005, 0.015, 0.03, 0.06, 0.105, 0.16, 0.23, 0.315, 0.425, 1.0], dtype='d'),
        'reco': None
    },
}


construct_all_fine_binnings(VAR_UNFOLD_DICT_TARGET0p5)
construct_all_fine_binnings(VAR_UNFOLD_DICT_TARGET0p6)

VAR_UNFOLD_DICT = VAR_UNFOLD_DICT_TARGET0p5
# VAR_UNFOLD_DICT = VAR_UNFOLD_DICT_TARGET0p6

# Common labels for legends etc
ZpJ_LABEL = "Z+jets region"
DY_ZpJ_LABEL = "DY+jets MC, Z+jets region"
DY_ZpJ_GROOMED_LABEL = "DY+jets MC, Z+jets region\n(groomed)"
DY_ZpJ_GEN_LABEL = "DY+jets MC, Z+jets region\n(GenJets)"
DY_ZpJ_QFLAV_LABEL = "DY+jets MC, Z+jets region\n(uds-matched)"
DY_ZpJ_GFLAV_LABEL = "DY+jets MC, Z+jets region\n(g-matched)"

DY_Dijet_LABEL = "DY+jets MC, Dijet region"
DY_Dijet_GEN_LABEL = "DY+jets MC, Dijet region\n(GenJets)"
DY_Dijet_QFLAV_LABEL = "DY+jets MC, Dijet region\n(uds-matched)"
DY_Dijet_GFLAV_LABEL = "DY+jets MC, Dijet region\n(g-matched)"

QCD_ZpJ_LABEL = "QCD MC, Z+jet region"
QCD_ZpJ_GEN_LABEL = "QCD MC, Z+jets region\n(GenJets)"
QCD_ZpJ_QFLAV_LABEL = "QCD MC, Z+jets region\n(uds-matched)"
QCD_ZpJ_GFLAV_LABEL = "QCD MC, Z+jets region\n(g-matched)"

Dijet_LABEL = "Dijet region"
Dijet_GROOMED_LABEL = "Dijet region (groomed)"
QCD_Dijet_LABEL = "QCD MC, Dijet region"
QCD_Dijet_GEN_LABEL = "QCD MC, Dijet region\n(GenJets)"
QCD_Dijet_QFLAV_LABEL = "QCD MC, Dijet region\n(uds-matched)"
QCD_Dijet_GFLAV_LABEL = "QCD MC, Dijet region\n(g-matched)"

Dijet_CEN_LABEL = "Dijet region (central jet)"
QCD_Dijet_CEN_LABEL = "QCD MC, Dijet region (central jet)"
QCD_Dijet_CEN_GROOMED_LABEL = "QCD MC, Dijet region (central jet, groomed)"
QCD_Dijet_CEN_GEN_LABEL = "QCD MC, Dijet region (central jet)\n(GenJets)"
QCD_Dijet_CEN_QFLAV_LABEL = "QCD MC, Dijet region (central jet)\n(uds-matched)"
QCD_Dijet_CEN_GFLAV_LABEL = "QCD MC, Dijet region (central jet)\n(g-matched)"

Dijet_FWD_LABEL = "Dijet region (forward jet)"
QCD_Dijet_FWD_LABEL = "QCD MC, Dijet region (forward jet)"
QCD_Dijet_FWD_GROOMED_LABEL = "QCD MC, Dijet region (forward jet, groomed)"
QCD_Dijet_FWD_GEN_LABEL = "QCD MC, Dijet region (forward jet)\n(GenJets)"
QCD_Dijet_FWD_QFLAV_LABEL = "QCD MC, Dijet region (forward jet)\n(uds-matched)"
QCD_Dijet_FWD_GFLAV_LABEL = "QCD MC, Dijet region (forward jet)\n(g-matched)"

SINGLE_MU_LABEL = "SingleMu, Z+jets region"
JETHT_LABEL = "JetHT, Dijet region"
ZB_LABEL = "ZeroBias, Dijet region"
JETHT_ZB_LABEL = "JetHT+ZeroBias, Dijet region"

JETHT_CEN_LABEL = "JetHT, Dijet region (central jet)"
ZB_CEN_LABEL = "ZeroBias, Dijet region (central jet)"
JETHT_ZB_CEN_LABEL = "JetHT+ZeroBias, Dijet region (central jet)"

JETHT_FWD_LABEL = "JetHT, Dijet region (forward jet)"
ZB_FWD_LABEL = "ZeroBias, Dijet region (forward jet)"
JETHT_ZB_FWD_LABEL = "JetHT+ZeroBias, Dijet region (forward jet)"

# Dirs in ROOT files
ZPJ_RECOJET_RDIR = "ZPlusJets_QG"
DJ_RECOJET_RDIR = "Dijet_QG_tighter"
ZPJ_GENJET_RDIR = "ZPlusJets_genjet"
DJ_GENJET_RDIR = "Dijet_genjet"

# Common filenames
DY_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root"
QCD_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_QCD.root"
QCD_PYTHIA_ONLY_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_PYTHIA-QCD.root"
ZPJ_ALL_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_MGPYTHIA_DYJetsToLL_M-50_HT-All.root"

QCD_HERWIG_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD.root"
QCD_HERWIG_PTREWEIGHT_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD_PtReweight.root"
DY_HERWIG_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_HERWIG_DYJetsToLL.root"
DY_MG_HERWIG_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_MG_HERWIG_DYJetsToLL.root"
DY_AMCATNLO_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_AMCATNLO_DYJetsToLL.root"

SINGLE_MU_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_SingleMu.root"
JETHT_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT.root"
ZB_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_ZeroBias.root"
JETHT_ZB_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_ZeroBias.root"

TTBAR_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_TTBAR.root"
WW_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_WW.root"
WZ_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_WZ.root"
ZZ_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_ZZ.root"

JETHT_RUNB_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunB.root"
JETHT_RUNC_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunC.root"
JETHT_RUND_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunD.root"
JETHT_RUNE_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunE.root"
JETHT_RUNF_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunF.root"
JETHT_RUNG_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunG.root"
JETHT_RUNH_FILENAME = "uhh2.AnalysisModuleRunner.DATA.Data_JetHT_RunH.root"


def get_parser():
    """Return a parser to loop over several input dirs"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workdirs',
                        nargs='+',
                        help='Workdir(s) with ROOT files to process. '
                        'Each directory must have ROOT files '
                        '"uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" '
                        'and "uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" in it. '
                        'Several dirs can be specified here, separated by a space.')
    parser.add_argument("-o", "--output",
                        help="Directory to put output plot dirs into",
                        default=".")
    return parser


def extract_jet_config(workdir):
    """Return jet algo/PUS from dir name"""
    import re
    res = re.search(r"(ak|ca)([0-9]+)(chs|puppi)?", workdir, flags=re.IGNORECASE)
    if res:
        algo = res.groups()[0]
        r = res.groups()[1]
        pus = ""
        if len(res.groups()) == 3:
            pus = res.groups()[2]
        return "%s%s %s" % (algo.upper(), r, pus.upper())
    else:
        print("Couldn't decypher jet algo")
        return ""
