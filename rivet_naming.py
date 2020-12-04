"""
Common functions, etc for RIVET/YODA naming scheme for plots
"""

from collections import namedtuple


# PLOT NAMING SCHEME
# --------------------------------------------------------------------------
# d = channel (ak4/8 radii (10, 20) + {dijet cen / fwd, or Z+Jet} * groomed versions[1..4])
# x = lambda variable; neutral+charged & charged-only are treated separately
# y = pT bin


# The user should iterate over these collections, to be able to use get_plot_name()

# Jet radii
Radius = namedtuple("Radius", ['name'])

JET_RADII = [
    Radius("AK4"),
    Radius("AK8"),
]

# Rivet "directory" name for hists
DIJET_PATH = "CMS_2018_PAS_SMP_18_QGX_DIJET"

# Selection regions
# format: (ROOT TDirectory name, plot label)
# Note that both dijet and Z+J start from the same index
# (e.g. d11 for both AK4 dijet central, and AK4 Z+Jet)
# since the two schemes shouldn't overlap in reality
Region = namedtuple("Region", ["name", "label", "is_groomed"])

DIJET_REGIONS = [
    Region(name="Dijet_central", label="Central dijet", is_groomed=False),
    Region(name="Dijet_central_groomed", label="Central dijet", is_groomed=True),
    Region(name="Dijet_forward", label="Forward dijet", is_groomed=False),
    Region(name="Dijet_forward_groomed", label="Forward dijet", is_groomed=True),
]


ZPJ_PATH = "CMS_2018_PAS_SMP_18_QGX_ZPJ"

ZPJ_REGIONS = [
    Region(name="ZPlusJets", label="Z+jet", is_groomed=False),
    Region(name="ZPlusJets_groomed", label="Z+jet", is_groomed=True),
]


# Lambda (angle) variables
LambdaVar = namedtuple("LambdaVar", ["hist_name", "label", "lambda_label", "xlimit"])

LAMBDA_VARS = [
    LambdaVar(hist_name='jet_puppiMultiplicity',
              label="Multiplicity $\\lambda_{0}^{0}$",
              lambda_label="\\lambda_{0}^{0}",
              xlimit=(0, 150)),
    LambdaVar(hist_name='jet_pTD',
              label="$(p_{T}^{D})^{2}$ $\\lambda_{0}^{2}$",
              lambda_label="\\lambda_{0}^{2}",
              xlimit=(0, 1)),
    LambdaVar(hist_name='jet_LHA',
              label="LHA $\\lambda_{0.5}^{1}$",
              lambda_label="\\lambda_{0.5}^{1}",
              xlimit=(0, 1)),
    LambdaVar(hist_name='jet_width',
              label="Width $\\lambda_{1}^{1}$",
              lambda_label="\\lambda_{1}^{1}",
              xlimit=(0, 1)),
    LambdaVar(hist_name='jet_thrust',
              label="Thrust $\\lambda_{2}^{1}$",
              lambda_label="\\lambda_{2}^{1}",
              xlimit=(0, 0.5)),

    LambdaVar(hist_name='jet_puppiMultiplicity_charged',
              label="Multiplicity $\\lambda_{0}^{0}$ (charged only)",
              lambda_label="\\lambda_{0}^{0}",
              xlimit=(0, 150)),
    LambdaVar(hist_name='jet_pTD_charged',
              label="$(p_{T}^{D})^{2}$ $\\lambda_{0}^{2}$ (charged only)",
              lambda_label="\\lambda_{0}^{2}",
              xlimit=(0, 1)),
    LambdaVar(hist_name='jet_LHA_charged',
              label="LHA $\\lambda_{0.5}^{1}$ (charged only)",
              lambda_label="\\lambda_{0.5}^{1}",
              xlimit=(0, 1)),
    LambdaVar(hist_name='jet_width_charged',
              label="Width $\\lambda_{1}^{1}$ (charged only)",
              lambda_label="\\lambda_{1}^{1}",
              xlimit=(0, 1)),
    LambdaVar(hist_name='jet_thrust_charged',
              label="Thrust $\\lambda_{2}^{1}$ (charged only)",
              lambda_label="\\lambda_{2}^{1}",
              xlimit=(0, 0.5)),
]


PT_BINS_EDGES_DIJET = [
    50, 65, 88, 120, 150, 186, 254, 326, 408, 481, 614, 800, 1000, 4000
]

PT_BINS_EDGES_ZPJ = [
    50, 65, 88, 120, 150, 186, 254, 326, 408, 1500
]

# convert to pairs of (lower, upper) bins edges
PT_BINS_DIJET = list(zip(PT_BINS_EDGES_DIJET[:-1], PT_BINS_EDGES_DIJET[1:]))

PT_BINS_ZPJ = list(zip(PT_BINS_EDGES_ZPJ[:-1], PT_BINS_EDGES_ZPJ[1:]))


def get_plot_name(jet_radius, region, lambda_variable, pt_bin):
    """Get base RIVET plot name, given objects that specify
    jet radius, signal region, lambda var, pt bin.

    Parameters
    ----------
    jet_radius : str
        One of the JET_RADII entries
    region : tuple(str, str)
        One of the DIJET_REGIONS or ZPJ_REGIONS entries
    lambda_variable : LambdaVar
        One of the LAMBDA_VARS entrues
    pt_bin : tuple(float, float)
        One of the PT_BINS_DIJET or PT_BINS_ZPJ entries

    Returns
    -------
    str
        base RIVET plot name
    """
    def _get_index(list_thing, item, name):
        try:
            # +1 as everything in RIVET is 1-indexed
            return list_thing.index(item)+1
        except ValueError as e:
            raise ValueError("Cannot find this bin ('%s') in list of %s: %s" % (item, name, list_thing))

    radius_ind = _get_index(JET_RADII, jet_radius, "jet radii")

    is_zpj = False
    try:
        region_ind = _get_index(DIJET_REGIONS, region, "dijet regions")
    except ValueError:
        try:
            region_ind = _get_index(ZPJ_REGIONS, region, "Z+jet regions")
            is_zpj = True
        except ValueError:
            raise ValueError("Cannot find region ('%s') in either list of dijet or Z+J regions" % (region))

    lambda_ind = _get_index(LAMBDA_VARS, lambda_variable, "lambda variables")

    pt_ind = _get_index(PT_BINS_ZPJ if is_zpj else PT_BINS_DIJET, pt_bin, "pt bins")

    channel_ind = (10*radius_ind) + region_ind
    name = "d{:0>2}-x{:0>2}-y{:0>2}".format(channel_ind, lambda_ind, pt_ind)
    return name
