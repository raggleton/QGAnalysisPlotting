"""Common region dict generating function for unfolding"""


import os

import ROOT
from MyStyle import My_Style
My_Style.cd()

import qg_common as qgc

ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat(".3f")
ROOT.gStyle.SetHistTopMargin(0.)


# FOR DIJET:
def get_dijet_config(source_dir, central=True, groomed=False):
    source_dir_systs = os.path.join(source_dir, "systematics_files")

    input_mc_qcd_mgpythia_tfile = os.path.join(source_dir, qgc.QCD_FILENAME)
    input_mc_qcd_pythia_tfile = os.path.join(source_dir, qgc.QCD_PYTHIA_ONLY_FILENAME)
    input_mc_qcd_herwig_tfile = os.path.join(source_dir, qgc.QCD_HERWIG_FILENAME)
    input_mc_qcd_herwig_tfile_reweight = os.path.join(source_dir, qgc.QCD_HERWIG_PTREWEIGHT_FILENAME)

    input_jetht_tfile = os.path.join(source_dir, qgc.JETHT_ZB_FILENAME)

    dijet_region_dict_template = {
        "name": "Dijet",
        "dirname": "Dijet_QG_Unfold_central_tighter",
        "label": "Dijet",
        "data_tfile": input_jetht_tfile,
        "mc_tfile": input_mc_qcd_mgpythia_tfile,
        "mc_label": "MG+Pythia8",
        # "mc_tfile": input_mc_qcd_herwig_tfile,
        # "mc_label": "Herwig++",
        # "mc_tfile": input_mc_qcd_pythia_tfile,
        # "mc_label": "Pythia8",
        "alt_mc_tfile": input_mc_qcd_herwig_tfile,
        "alt_mc_label": "Herwig++",
        # "alt_mc_tfile": input_mc_qcd_pythia_tfile,
        # "alt_mc_label": "Pythia8",
        # "alt_mc_tfile": input_mc_qcd_herwig_tfile_reweight,
        # "alt_mc_label": "Herwig++ (p_{T} reweight)",
        "tau_limits": None,  # user should set this
        "unreg_unfolder": None,  # set later if regularisation used
        "experimental_systematics": [
            # {
            #     "label": "Charged hadron up",
            #     "tfile": os.path.join(source_dir_systs, 'chargedHadronShiftUp', qgc.QCD_FILENAME),
            #     "colour": ROOT.kAzure+1,
            # },
            # {
            #     "label": "Charged hadron down",
            #     "tfile": os.path.join(source_dir_systs, 'chargedHadronShiftDown', qgc.QCD_FILENAME),
            #     "colour": ROOT.kAzure+1,
            #     "linestyle": 2,
            # },
            # {
            #     "label": "Neutral hadron up",
            #     "tfile": os.path.join(source_dir_systs, 'neutralHadronShiftUp', qgc.QCD_FILENAME),
            #     "colour": ROOT.kOrange-4,
            # },
            # {
            #     "label": "Neutral hadron down",
            #     "tfile": os.path.join(source_dir_systs, 'neutralHadronShiftDown', qgc.QCD_FILENAME),
            #     "colour": ROOT.kOrange-4,
            #     "linestyle": 2,
            # },
            # {
            #     "label": "Photon up",
            #     "tfile": os.path.join(source_dir_systs, 'photonShiftUp', qgc.QCD_FILENAME),
            #     "colour": ROOT.kMagenta-3,
            # },
            # {
            #     "label": "Photon down",
            #     "tfile": os.path.join(source_dir_systs, 'photonShiftDown', qgc.QCD_FILENAME),
            #     "colour": ROOT.kMagenta-3,
            #     "linestyle": 2,
            # },
            # {
            #     "label": "JEC up",
            #     "tfile": os.path.join(source_dir_systs, 'jecsmear_directionUp', qgc.QCD_FILENAME),
            #     "colour": ROOT.kGreen+2,
            # },
            # {
            #     "label": "JEC down",
            #     "tfile": os.path.join(source_dir_systs, 'jecsmear_directionDown', qgc.QCD_FILENAME),
            #     "colour": ROOT.kGreen+2,
            #     "linestyle": 2,
            # },
            # {
            #     "label": "JER up",
            #     "tfile": os.path.join(source_dir_systs, 'jersmear_directionUp', qgc.QCD_FILENAME),
            #     "colour": ROOT.kOrange+3,
            # },
            # {
            #     "label": "JER down",
            #     "tfile": os.path.join(source_dir_systs, 'jersmear_directionDown', qgc.QCD_FILENAME),
            #     "colour": ROOT.kOrange+3,
            #     "linestyle": 2,
            # },
            # {
            #     "label": "Pileup up",
            #     "tfile": os.path.join(source_dir_systs, 'pileup_directionUp', qgc.QCD_FILENAME),
            #     "colour": ROOT.kBlue-4,
            # },
            # {
            #     "label": "Pileup down",
            #     "tfile": os.path.join(source_dir_systs, 'pileup_directionDown', qgc.QCD_FILENAME),
            #     "colour": ROOT.kBlue-4,
            #     "linestyle": 2,
            # },
            {
                "label": "Tracking up",
                "tfile": os.path.join(source_dir_systs, 'track_directionUp', qgc.QCD_FILENAME),
                "colour": ROOT.kMagenta+3,
            },
            {
                "label": "Tracking down",
                "tfile": os.path.join(source_dir_systs, 'track_directionDown', qgc.QCD_FILENAME),
                "colour": ROOT.kMagenta+3,
                "linestyle": 2,
            },
            {
                "label": "Herwig++",
                "tfile": input_mc_qcd_herwig_tfile,
                "colour": ROOT.kOrange-3,
            },
            {
                "label": "Herwig++ (p_{T} reweight)",
                "tfile": input_mc_qcd_herwig_tfile_reweight,
                "colour": ROOT.kOrange+4,
            },

        ],
        "model_systematics": [
            # {
            #     "label": "muR up, muF nominal",
            #     "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNominal', qgc.QCD_FILENAME),
            #     "colour": ROOT.kAzure,
            #     "unfolder": None,
            # },
            # {
            #     "label": "muR down, muF nominal",
            #     "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFNominal', qgc.QCD_FILENAME),
            #     "colour": ROOT.kAzure+1,
            #     "unfolder": None,
            # },
            # {
            #     "label": "muR nominal, muF up",
            #     "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFUp', qgc.QCD_FILENAME),
            #     "colour": ROOT.kAzure+2,
            #     "unfolder": None,
            # },
            # {
            #     "label": "muR nominal, muF down",
            #     "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRNominal_ScaleVariationMuFDown', qgc.QCD_FILENAME),
            #     "colour": ROOT.kAzure+3,
            #     "unfolder": None,
            # },
            # {
            #     "label": "muR down, muF down",
            #     "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFDown', qgc.QCD_FILENAME),
            #     "colour": ROOT.kAzure+4,
            #     "unfolder": None,
            # },
            # {
            #     "label": "muR up, muF up",
            #     "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFUp', qgc.QCD_FILENAME),
            #     "colour": ROOT.kAzure+5,
            #     "unfolder": None,
            # },
            {
                "label": "Herwig++",
                "tfile": input_mc_qcd_herwig_tfile,
                "colour": ROOT.kOrange-3,
                "unfolder": None,
            },
            {
                "label": "Herwig++ (p_{T} reweight)",
                "tfile": input_mc_qcd_herwig_tfile_reweight,
                "colour": ROOT.kOrange+4,
                "unfolder": None,
            },
            # {
            #     "label": "Pythia8",
            #     "tfile": input_mc_qcd_pythia_tfile,
            #     "colour": ROOT.kGreen-3,
            #     "unfolder": None,
            # },
        ],
        "pdf_systematics": [
            {
                "label": "PDF",  # this is a tempalte entry, used for future
                "tfile": os.path.join(source_dir_systs, 'PDFvariationsTrue', qgc.QCD_FILENAME),
                "colour": ROOT.kCyan+2,
                "unfolder": None,
                "variations": range(100),  # list of all the variation #s to be used
            },
        ]
    }

    if central and not groomed:
        this_dict = dijet_region_dict_template.copy()
        this_dict['dirname'] = 'Dijet_QG_Unfold_central_tighter'
        this_dict['label'] = 'Dijet central'
        this_dict['name'] = 'Dijet_central'
        return this_dict

    elif not central and not groomed:
        this_dict = dijet_region_dict_template.copy()
        this_dict['dirname'] = 'Dijet_QG_Unfold_forward_tighter'
        this_dict['label'] = 'Dijet forward'
        this_dict['name'] = 'Dijet_forward'
        return this_dict

    elif central and groomed:
        this_dict = dijet_region_dict_template.copy()
        this_dict['dirname'] = 'Dijet_QG_Unfold_central_tighter_groomed'
        this_dict['label'] = 'Dijet central'
        this_dict['name'] = 'Dijet_central_groomed'
        return this_dict

    elif not central and groomed:
        this_dict = dijet_region_dict_template.copy()
        this_dict['dirname'] = 'Dijet_QG_Unfold_forward_tighter_groomed'
        this_dict['label'] = 'Dijet forward'
        this_dict['name'] = 'Dijet_forward_groomed'
        return this_dict


def get_zpj_config(source_dir, groomed=False):
    source_dir_systs = os.path.join(source_dir, "systematics_files")

    input_mc_dy_mgpythia_tfile = os.path.join(source_dir, qgc.DY_FILENAME)
    input_mc_dy_mgherwig_tfile = os.path.join(source_dir, qgc.DY_MG_HERWIG_FILENAME)
    input_mc_dy_herwig_tfile = os.path.join(source_dir, qgc.DY_HERWIG_FILENAME)

    input_singlemu_tfile = os.path.join(source_dir, qgc.SINGLE_MU_FILENAME)

    zpj_region_dict = {
        "name": "ZPlusJets",
        "dirname": "ZPlusJets_QG_Unfold",
        "label": "Z+jets",
        "data_tfile": input_singlemu_tfile,
        "mc_tfile": input_mc_dy_mgpythia_tfile,
        "mc_label": "MG+Pythia8",
        "alt_mc_tfile": input_mc_dy_mgherwig_tfile,
        "alt_mc_label": "MG+Herwig++",
        "tau_limits": None,
        "backgrounds": [
            {
                "name": "t#bar{t}",
                "tfile": os.path.join(source_dir, qgc.TTBAR_FILENAME),
                "rate_unc": 1.,
            },
            {
                "name": "WW",
                "tfile": os.path.join(source_dir, qgc.WW_FILENAME),
                "rate_unc": 1.
            },
            {
                "name": "WZ",
                "tfile": os.path.join(source_dir, qgc.WZ_FILENAME),
                "rate_unc": 1.
            },
            {
                "name": "ZZ",
                "tfile": os.path.join(source_dir, qgc.ZZ_FILENAME),
                "rate_unc": 1.
            },
        ],
        "experimental_systematics": [
            # {
            #     "label": "Luminosity up",
            #     "tfile": None,
            #     "factor": 0.025,
            #     "colour": ROOT.kCyan,
            # },
            # {
            #     "label": "Luminosity down",
            #     "tfile": None,
            #     "factor": -0.025,
            #     "colour": ROOT.kCyan,
            #     'linestyle': 2,
            # },
            # {
            #     "label": "Charged hadron up",
            #     "tfile": os.path.join(source_dir_systs, 'chargedHadronShiftUp', qgc.DY_FILENAME),
            #     "colour": ROOT.kAzure+1,
            # },
            # {
            #     "label": "Charged hadron down",
            #     "tfile": os.path.join(source_dir_systs, 'chargedHadronShiftDown', qgc.DY_FILENAME),
            #     "colour": ROOT.kAzure+1,
            #     "linestyle": 2,
            # },
            # {
            #     "label": "Neutral hadron up",
            #     "tfile": os.path.join(source_dir_systs, 'neutralHadronShiftUp', qgc.DY_FILENAME),
            #     "colour": ROOT.kOrange-3,
            # },
            # {
            #     "label": "Neutral hadron down",
            #     "tfile": os.path.join(source_dir_systs, 'neutralHadronShiftDown', qgc.DY_FILENAME),
            #     "colour": ROOT.kOrange-3,
            #     "linestyle": 2,
            # },
            # {
            #     "label": "Photon up",
            #     "tfile": os.path.join(source_dir_systs, 'photonShiftUp', qgc.DY_FILENAME),
            #     "colour": ROOT.kMagenta-3,
            # },
            # {
            #     "label": "Photon down",
            #     "tfile": os.path.join(source_dir_systs, 'photonShiftDown', qgc.DY_FILENAME),
            #     "colour": ROOT.kMagenta-3,
            #     "linestyle": 2,
            # },
            # {
            #     "label": "JEC up",
            #     "tfile": os.path.join(source_dir_systs, 'jecsmear_directionUp', qgc.DY_FILENAME),
            #     "colour": ROOT.kGreen+2,
            # },
            # {
            #     "label": "JEC down",
            #     "tfile": os.path.join(source_dir_systs, 'jecsmear_directionDown', qgc.DY_FILENAME),
            #     "colour": ROOT.kGreen+2,
            #     "linestyle": 2,
            # },
            # {
            #     "label": "JER up",
            #     "tfile": os.path.join(source_dir_systs, 'jersmear_directionUp', qgc.DY_FILENAME),
            #     "colour": ROOT.kOrange+3,
            # },
            # {
            #     "label": "JER down",
            #     "tfile": os.path.join(source_dir_systs, 'jersmear_directionDown', qgc.DY_FILENAME),
            #     "colour": ROOT.kOrange+3,
            #     "linestyle": 2,
            # },
            # {
            #     "label": "Pileup up",
            #     "tfile": os.path.join(source_dir_systs, 'pileup_directionUp', qgc.DY_FILENAME),
            #     "colour": ROOT.kBlue-4,
            # },
            # {
            #     "label": "Pileup down",
            #     "tfile": os.path.join(source_dir_systs, 'pileup_directionDown', qgc.DY_FILENAME),
            #     "colour": ROOT.kBlue-4,
            #     "linestyle": 2,
            # },
            {
                "label": "Tracking up",
                "tfile": os.path.join(source_dir_systs, 'track_directionUp', qgc.DY_FILENAME),
                "colour": ROOT.kMagenta+3,
            },
            {
                "label": "Tracking down",
                "tfile": os.path.join(source_dir_systs, 'track_directionDown', qgc.DY_FILENAME),
                "colour": ROOT.kMagenta+3,
                "linestyle": 2,
            },
            # {
            #     "label": "MG+Herwig++",
            #     "tfile": input_mc_dy_mgherwig_tfile,
            #     "colour": ROOT.kOrange-3,
            # },
            # {
            #     "label": "Herwig++",
            #     "tfile": input_mc_dy_herwig_tfile,
            #     "colour": ROOT.kBlue-3,
            # },
        ],
        "model_systematics": [
            {
                "label": "muR up, muF nominal",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFNom', qgc.DY_FILENAME),
                "colour": ROOT.kAzure,
                "unfolder": None,
            },
            {
                "label": "muR down, muF nominal",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFNom', qgc.DY_FILENAME),
                "colour": ROOT.kAzure+1,
                "unfolder": None,
            },
            {
                "label": "muR nominal, muF up",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRNom_ScaleVariationMuFUp', qgc.DY_FILENAME),
                "colour": ROOT.kAzure+2,
                "unfolder": None,
            },
            {
                "label": "muR nominal, muF down",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRNom_ScaleVariationMuFDown', qgc.DY_FILENAME),
                "colour": ROOT.kAzure+3,
                "unfolder": None,
            },
            {
                "label": "muR down, muF down",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRDown_ScaleVariationMuFDown', qgc.DY_FILENAME),
                "colour": ROOT.kAzure+4,
                "unfolder": None,
            },
            {
                "label": "muR up, muF up",
                "tfile": os.path.join(source_dir_systs, 'ScaleVariationMuRUp_ScaleVariationMuFUp', qgc.DY_FILENAME),
                "colour": ROOT.kAzure+5,
                "unfolder": None,
            },
            # {
            #     "label": "MG+Herwig++",
            #     "tfile": input_mc_dy_mgherwig_tfile,
            #     "colour": ROOT.kOrange-3,
            #     "unfolder": None,
            # },
            # {
            #     "label": "Herwig++",
            #     "tfile": input_mc_dy_herwig_tfile,
            #     "colour": ROOT.kBlue-3,
            #     "unfolder": None,
            # },
        ],
        "pdf_systematics": [
            {
                "label": "PDF",
                "tfile": os.path.join(source_dir_systs, 'PDFvariationsTrue', qgc.DY_FILENAME),
                "colour": ROOT.kCyan+2,
                "unfolder": None,
                "variations": range(100),
            },
        ]
    }

    if not groomed:
        return zpj_region_dict.copy()
    else:
        this_dict = zpj_region_dict.copy()
        this_dict['dirname'] = 'ZPlusJets_QG_Unfold_groomed'
        this_dict['name'] = 'ZPlusJets_groomed'
        this_dict['label'] = 'Z+jets'
        return this_dict
