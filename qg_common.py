"""Common vars for doing QG plots"""


from collections import namedtuple


DY_COLOUR = 880 # kViolet
QCD_COLOUR = 867 #kAzure+7

DY_COLOURS = [880, 871, 882, 884]
QCD_COLOURS = [867, 868, 864, 853]

PT_BINS = [(80, 100), (100, 200), (400, 500), (1000, 2000), (80, 2000)]
THEORY_PT_BINS = [(100, 200), (400, 500), (1000, 2000), (80, 2000)]
THEORY_PT_BINS = [(100, 200), (200, 300), (300, 400), (400, 500), (500, 600), (600, 800), (800, 1000), (1000, 2000), (80, 2000)]

Angle = namedtuple("Angle", ['var', 'kappa', 'beta', 'name', "lambda_str", "colour"])
COMMON_VARS_WITH_FLAV = [
    Angle('jet_multiplicity', 0, 0, "Multiplicity", "#lambda_{0}^{0}", 2),
    Angle('jet_pTD', 2, 0, "(p_{T}^{D})^{2}", "#lambda_{0}^{2}", 418),
    Angle('jet_LHA', 1, 0.5, "LHA", "#lambda_{0.5}^{1}", 600),
    Angle('jet_width', 1, 1, "Width", "#lambda_{1}^{1}", 861),
    Angle('jet_thrust', 1, 2, "Thrust", "#lambda_{2}^{1}", 617),
    Angle('jet_flavour', 0, 0, "Flavour", "PDGID", 7),
    Angle("jet_genParton_flavour", 0, 0, "Flavour", "PDGID", 8)
]

COMMON_VARS = COMMON_VARS_WITH_FLAV[:]

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

# Dirs in ROOT files
ZPJ_RECOJET_RDIR = "ZPlusJets_QG"
DJ_RECOJET_RDIR = "Dijet_QG"
ZPJ_GENJET_RDIR = "ZPlusJets_genjet"
DJ_GENJET_RDIR = "Dijet_genjet"
