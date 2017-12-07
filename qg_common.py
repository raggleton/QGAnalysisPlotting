"""Common vars for doing QG plots"""

import argparse
from collections import namedtuple


# Purple & blue look ~ same to red-green colourblind people
# (V.Dark blue & light pruple do look a bit diff)
# Red and blue does look different
DY_COLOUR =  632 # kRed 
QCD_COLOUR = 867 # kAzure+7

DY_COLOUR =  880 # kViolet
# QCD_COLOUR = 600 # kBlue

# You shoudl use markers as well for colourblindness
DY_MARKER = 20
QCD_MARKER = 22

# When comparing e.g. PU bins
# FIXME: update for colourblindness
DY_COLOURS = [880, 881, 884, 871]
QCD_COLOURS = [867, 600, 853, 425]

PT_BINS = [(80, 100), (100, 200), (400, 500), (1000, 2000), (80, 2000)]
THEORY_PT_BINS = [(100, 200), (400, 500), (1000, 2000), (80, 2000)]
THEORY_PT_BINS = [(100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 800), (800, 1000), (1000, 2000), (80, 2000)]
# THEORY_PT_BINS = [(20, 40), (40, 60), (60, 80), (80, 100), (100, 120),
THEORY_PT_BINS = [(20, 25), (40, 45), (50, 55), (60, 70), (70, 80), (80, 90), (90, 100), (100, 120),
                  (120, 160), (160, 200), (200, 220), (280, 300), (300, 400), (400, 410),
                  (400, 500), (500, 600), (500, 550), (600, 800), (800, 1000),
                  (1000, 1400), (1400, 2000)]

Angle = namedtuple("Angle", ['var', 'kappa', 'beta', 'name', "lambda_str", "colour"])
COMMON_VARS_WITH_FLAV = [
    Angle("jet_puppiMultiplicity", 0, 0, "PUPPI Multiplicity", "#lambda_{0}^{0} (PUPPI)", 2),
    Angle('jet_multiplicity', 0, 0, "Multiplicity", "#lambda_{0}^{0}", 2),
    Angle('jet_pTD', 2, 0, "(p_{T}^{D})^{2}", "#lambda_{0}^{2}", 418),
    Angle('jet_LHA', 1, 0.5, "LHA", "#lambda_{0.5}^{1}", 600),
    Angle('jet_width', 1, 1, "Width", "#lambda_{1}^{1}", 861),
    Angle('jet_thrust', 1, 2, "Thrust", "#lambda_{2}^{1}", 617),
    Angle('jet_flavour', 0, 0, "Flavour", "PDGID", 7),
    Angle("jet_genParton_flavour", 0, 0, "Flavour", "PDGID", 8)
]

COMMON_VARS = COMMON_VARS_WITH_FLAV[:-2]

DY_ZpJ_LABEL = "DY+j, Z+jets selection"
DY_ZpJ_GEN_LABEL = "DY+j, Z+jets selection\n(GenJets)"
DY_ZpJ_QFLAV_LABEL = "DY+j, Z+jets selection\n(uds-matched)"
DY_ZpJ_GFLAV_LABEL = "DY+j, Z+jets selection\n(g-matched)"

DY_Dijet_LABEL = "DY+j, Dijet selection"
DY_Dijet_GEN_LABEL = "DY+j, Dijet selection\n(GenJets)"
DY_Dijet_QFLAV_LABEL = "DY+j, Dijet selection\n(uds-matched)"
DY_Dijet_GFLAV_LABEL = "DY+j, Dijet selection\n(g-matched)"

QCD_ZpJ_LABEL = "QCD, Z+jets selection"
QCD_ZpJ_GEN_LABEL = "QCD, Z+jets selection\n(GenJets)"
QCD_ZpJ_QFLAV_LABEL = "QCD, Z+jets selection\n(uds-matched)"
QCD_ZpJ_GFLAV_LABEL = "QCD, Z+jets selection\n(g-matched)"

QCD_Dijet_LABEL = "QCD, Dijet selection"
QCD_Dijet_GEN_LABEL = "QCD, Dijet selection\n(GenJets)"
QCD_Dijet_QFLAV_LABEL = "QCD, Dijet selection\n(uds-matched)"
QCD_Dijet_GFLAV_LABEL = "QCD, Dijet selection\n(g-matched)"

# Dirs in ROOT files
ZPJ_RECOJET_RDIR = "ZPlusJets_QG"
DJ_RECOJET_RDIR = "Dijet_QG"
ZPJ_GENJET_RDIR = "ZPlusJets_genjet"
DJ_GENJET_RDIR = "Dijet_genjet"

# Common filenames
DY_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root"
QCD_FILENAME = "uhh2.AnalysisModuleRunner.MC.MC_QCD_.root"


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
    parser.add_argument("-o", "--output", help="Directory to put output plot dirs into", default=None)
    return parser
