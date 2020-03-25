#!/usr/bin/env python


"""Script/module for unfolding logistics e.g. arg parser, output directory.

Useful for e.g. batch jobs
"""


from __future__ import print_function, division

import os
os.nice(10)
import sys
import argparse
import qg_common as qgc
import distutils
from distutils import util
import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(1)
ROOT.TH1.SetDefaultSumw2()


AREA_OPT_DICT = {'Area': ROOT.TUnfold.kEConstraintArea, 'None': ROOT.TUnfold.kEConstraintNone}
AREA_OPT_INV_DICT = {v: k for k, v in AREA_OPT_DICT.items()}


def get_unfolding_argparser(description='', parser=None):
    """Create ArgumentParser to handle various unfolding options

    To be used in unfolding.py

    If parser is None, creates one using description,
    otherwise adds all arguments etc to the object passed in
    """
    parser = parser or argparse.ArgumentParser(description=description)
    parser.add_argument("source",
                        help="Source directory with ROOT files")

    parser.add_argument("--angles",
                        choices=list(qgc.VAR_UNFOLD_DICT.keys()) + ["all"],
                        nargs='+',
                        help="Lambda angles to unfold, or 'all' for all of them")

    standard_bool_description = (" True values are 'y', 'yes', 't', 'true', 'on' and '1'."
                                 " False values are 'n', 'no', 'f', 'false', 'off' and '0'.")

    parser.add_argument("--doSummaryPlot",
                        type=lambda x:bool(distutils.util.strtobool(x)),
                        default=False,
                        help=('Do summary plot.'
                              + standard_bool_description))

    parser.add_argument("--noBinnedPlots",
                        action='store_true',
                        default=False,
                        help="Don't do binned plots")

    parser.add_argument("--outputDir",
                        default='',
                        help='Output directory')

    parser.add_argument("--areaConstraint",
                        default='None',
                        type=lambda x: AREA_OPT_DICT[x],
                        choices=AREA_OPT_DICT.keys(),
                        help='Area constraint.')

    # SIGNAL REGION OPTIONS
    # --------------------------------------------------------------------------
    region_group = parser.add_argument_group('Region selection')
    region_group.add_argument("--doAllRegions",
                              action='store_true',
                              help='Do unfolding for all regions (dijet, Z+J, groomed, ungroomed)')
    region_group.add_argument("--doDijetCentral",
                              action='store_true',
                              help='Do unfolding for dijet (central) jets')
    region_group.add_argument("--doDijetForward",
                              action='store_true',
                              help='Do unfolding for dijet (forward) jets')
    region_group.add_argument("--doDijetCentralGroomed",
                              action='store_true',
                              help='Do unfolding for groomed dijet (central) jets')
    region_group.add_argument("--doDijetForwardGroomed",
                              action='store_true',
                              help='Do unfolding for groomed dijet (forward) jets')
    region_group.add_argument("--doZPJ",
                              action='store_true',
                              help='Do unfolding for Z+jet jets')
    region_group.add_argument("--doZPJGroomed",
                              action='store_true',
                              help='Do unfolding for groomed Z+jet jets')

    # REGULARISATION OPTIONS
    # --------------------------------------------------------------------------
    regularization_group = parser.add_argument_group('Regularization options')
    regularization_group.add_argument("--regularize",
                                      choices=['None', 'tau', 'L'],
                                      default='None',
                                      help='Regularization scheme')
    regularization_group.add_argument("--regularizeAxis",
                                      choices=['both', 'pt', 'angle'],
                                      default='both',
                                      help='Axis to regularize')
    regularization_group.add_argument("--nScan",
                                      type=int,
                                      default=100,
                                      help='Number of scan points for regularization')
    regularization_group.add_argument("--biasFactor",
                                      type=float,
                                      default=0,
                                      help='Bias factor for regularization')

    # MC INPUT OPTIONS
    # --------------------------------------------------------------------------
    mc_group = parser.add_argument_group('MC input options')
    mc_group.add_argument("--MCinput",
                          type=lambda x: bool(distutils.util.strtobool(x)),
                          default=False,
                          help=('Unfold MC instead of data.'
                                + standard_bool_description))

    mc_group.add_argument("--MCsplit",
                          type=lambda x: bool(distutils.util.strtobool(x)),
                          default=False,
                          help=('Split MC between response & 1D reco, good for testing procedure.'
                                + standard_bool_description))

    # BACKGROUNDS OPTIONS
    # --------------------------------------------------------------------------
    bg_group = parser.add_argument_group("Backgrounds options")
    bg_group.add_argument("--subtractBackgrounds",
                          type=lambda x: bool(distutils.util.strtobool(x)),
                          default=False,
                          help=('Subtract true backgrounds (e.g. ttbar).'
                                + standard_bool_description))

    # EXPERIMENTAL SYST OPTIONS
    # --------------------------------------------------------------------------
    syst_group = parser.add_argument_group('Systematics options')
    syst_group.add_argument("--doExperimentalSysts",
                            type=lambda x: bool(distutils.util.strtobool(x)),
                            default=False,
                            help=('Do experimental systematics (i.e. those that modify response matrix).'
                                  + standard_bool_description))

    syst_group.add_argument("--doExperimentalSystsOnlyHerwig",
                            type=lambda x: bool(distutils.util.strtobool(x)),
                            default=False,
                            help=('Do only Herwig experimental systematics (i.e. those that modify response matrix).'
                                  + standard_bool_description))

    syst_group.add_argument("--doExperimentalSystsFromFile",
                            default=None,
                            help='Do experimental systematics (i.e. those that modify response matrix) ' \
                                 'but get shifts from previous unfolding. This should be a directory ' \
                                 'made by a previous running of unfolding.py ' \
                                 'that covers the different signal regions & variables')

    # MODEL SYST OPTIONS
    # --------------------------------------------------------------------------
    syst_group.add_argument("--doModelSysts",
                            type=lambda x: bool(distutils.util.strtobool(x)),
                            default=False,
                            help=('Do model systematics (i.e. those that modify input to be unfolded).'
                                   + standard_bool_description))

    syst_group.add_argument("--doModelSystsOnlyHerwig",
                            type=lambda x: bool(distutils.util.strtobool(x)),
                            default=False,
                            help=('Do only Herwig model systematics (i.e. those that modify input to be unfolded).'
                                   + standard_bool_description))

    syst_group.add_argument("--doModelSystsOnlyScale",
                            type=lambda x: bool(distutils.util.strtobool(x)),
                            default=False,
                            help=('Do only scale model systematics (i.e. those that modify input to be unfolded).'
                                   + standard_bool_description))

    syst_group.add_argument("--doModelSystsFromFile",
                            default=None,
                            help='Get model systematics from file. This should be a directory ' \
                                 'made by a previous running of unfolding.py ' \
                                 'that covers the different signal regions & variables')

    # PDF SYST OPTIONS
    # --------------------------------------------------------------------------
    syst_group.add_argument("--doPDFSysts",
                            type=lambda x: bool(distutils.util.strtobool(x)),
                            default=False,
                            help=('Do pdf systematics (may be slow!).'
                                   + standard_bool_description))

    syst_group.add_argument("--doPDFSystsFromFile",
                            default=None,
                            help='Get PDF systematics from file. This should be a directory ' \
                                 'made by a previous running of unfolding.py ' \
                                 'that covers the different signal regions & variables')

    # ALT RESPONSE
    # --------------------------------------------------------------------------
    syst_group.add_argument("--useAltResponse",
                            type=lambda x: bool(distutils.util.strtobool(x)),
                            default=False,
                            help=('Use alternate response matrix to unfold.'
                                   + standard_bool_description))
    return parser


def sanitise_args(args):
    """Do various post-parsing checks"""
    if args.doAllRegions:
        for x in ['doDijetCentral', 'doDijetForward', 'doDijetCentralGroomed', 'doDijetForwardGroomed', 'doZPJ', 'doZPJGroomed']:
            setattr(args, x, True)

    if not any([args.doDijetCentral, args.doDijetForward, args.doDijetCentralGroomed, args.doDijetForwardGroomed, args.doZPJ, args.doZPJGroomed]):
        raise RuntimeError("You need to specify at least one signal region e.g. --doDijetCentral")

    if not args.MCinput and args.doModelSysts:
        raise RuntimeError("You cannot do both model systs and run on data")

    if args.MCinput and args.subtractBackgrounds:
        print("")
        print("!!!! Cannot subtract backgrounds while using MC input, ignoring for now")
        print("")
        args.subtractBackgrounds = False

    if args.doPDFSysts and not args.MCinput:
        raise RuntimeError("Cannot do PDF systs and run over data")

    if (args.doModelSysts or args.doModelSystsOnlyHerwig or args.doModelSystsOnlyScale) and not args.MCinput:
        raise RuntimeError("Cannot do model systs and run over data")

    if args.doModelSystsOnlyHerwig and args.doModelSystsOnlyScale:
        raise RuntimeError("Cannot do both --doModelSystsOnlyHerwig and --doModelSystsOnlyScale")

    if (args.doExperimentalSysts or args.doExperimentalSystsOnlyHerwig) and args.doExperimentalSystsFromFile:
        args.doExperimentalSysts = False
        args.doExperimentalSystsOnlyHerwig = False
        print("Warning: will use experimental systs from --doExperimentalSystsFromFile option only, "
              "ignoring --doExperimentalSysts and --doExperimentalSystsOnlyHerwig")

    if (args.doModelSystsOnlyScale or args.doModelSystsOnlyHerwig or args.doModelSysts) and args.doModelSystsFromFile:
        args.doModelSysts = False
        args.doModelSystsOnlyScale = False
        args.doModelSystsOnlyHerwig = False
        print("Warning: will use model systs from --doModelSystsFromFile option only, "
              "ignoring --doModelSysts, --doModelSystsOnlyHerwig, -- doModelSystsOnlyScale")

    if args.doPDFSystsFromFile and args.doPDFSysts:
        args.doPDFSysts = False
        print("Warning: will use PDF systs from --doPDFSystsFromFile option only, ignoring --doPDFSysts")

    # if args.useAltResponse and args.doExperimentalSysts:
    #     args.doExperimentalSysts = False
    #     print("You cannot use both --useAltResponse and --doExperimentalSysts: disabling doExperimentalSysts")

    # # TODO handle both?
    # if args.useAltResponse and args.doModelSysts:
    #     args.doExperimentalSysts = False
    #     print("You cannot use both --useAltResponse and --doModelSysts: disabling doModelSysts")


def get_unfolding_output_dir(args):
    """Get name of outputdir for a given set of options

    args should be the parsed result from get_unfolding_argparser()
    """
    # Run with MC input instead of data
    mc_append = "_MC" if args.MCinput else "_DATA"

    if args.MCinput:
        mc_append += "_split" if args.MCsplit else "_all"

    SUBTRACT_FAKES = True  # this should alwys be True
    sub_append = "_subFakes" if SUBTRACT_FAKES else ""

    append = ""

    if args.subtractBackgrounds:
        append += "_subBkg"

    if args.doExperimentalSysts:
        append += "_experimentalSyst"

    if args.doExperimentalSystsOnlyHerwig:
        args.doExperimentalSysts = True
        append += "_experimentalSystOnlyHerwig"

    if args.doExperimentalSystsFromFile:
        append += "_experimentalSystFromFile"

    if args.doModelSysts:
        append += "_modelSyst"
        if not args.doExperimentalSysts:
            append += "NoExperimentalSyst"

    if args.doModelSystsOnlyHerwig:
        args.doModelSysts = True
        append += "_modelSystOnlyHerwig"
        if not args.doExperimentalSysts:
            append += "NoExperimentalSyst"

    elif args.doModelSystsOnlyScale:
        args.doModelSysts = True
        append += "_modelSystOnlyScale"
        if not args.doExperimentalSysts:
            append += "NoExperimentalSyst"

    if args.doModelSystsFromFile:
        append += "_modelSystFromFile"

    if args.doPDFSysts:
        append += "_pdfSyst"
        if not args.doExperimentalSysts:
            append += "NoExperimentalSyst"

    if args.doPDFSystsFromFile:
        append += "_pdfSystFromFile"

    if args.useAltResponse:
        append += "_altResponse"

    bias_str = ""
    if args.biasFactor != 0:
        bias_str = "_biasFactor%g" % args.biasFactor
        bias_str = bias_str.replace(".", "p")

    reg_axis_str = ""
    if args.regularize != "None":
        if args.regularizeAxis == 'pt':
            reg_axis_str = '_onlyRegPt'
        elif args.regularizeAxis == 'angle':
            reg_axis_str = '_onlyRegAngle'
        # reg_axis_str += "_onlyBinFactors"
        # reg_axis_str += "_invTruthMoreBins"
        reg_axis_str += "_invTruth"
        # reg_axis_str += "_invTruthUseUnfolded"

    regularize_str = "regularize%s%s%s" % (str(args.regularize).capitalize(), bias_str, reg_axis_str)

    str_parts = dict(
        regularize_str=regularize_str,
        mc_append=mc_append,
        area=AREA_OPT_INV_DICT[args.areaConstraint],
        append=append,
        sub_append=sub_append,
    )
    output_dir = "unfolding_{regularize_str}{mc_append}{sub_append}_densityModeBinWidth_constraint{area}{append}_signalRegionOnly".format(**str_parts)
    # Default to putting things into args.source, otherwise in wherever the user says
    if args.outputDir:
        output_dir = os.path.join(args.outputDir, output_dir)
    else:
        output_dir = os.path.join(args.source, output_dir)
    return output_dir


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    subparser = parser.add_subparsers(dest='out',
                                      help="Various sub commands. " \
                                           "All other arguments are passed to the main unfolding parser")
    unfolding_parser = subparser.add_parser('out', help='Get output directory name.')
    get_unfolding_argparser(description='', parser=unfolding_parser)
    args = parser.parse_args()
    sanitise_args(args)
    # print(args)
    if args.out:
        print(get_unfolding_output_dir(args))




