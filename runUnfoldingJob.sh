#!/bin/bash -e

# Run all the unfolding jobs for a given signal region, for a given angle
# Script gets given the following parameters:
# parameters   = CHANNEL ANGLE

# Need this as not empty at logon
export LD_LIBRARY_PATH=""

# Setup CMSSW etc
CMSSW_AREA=/nfs/dust/cms/user/aggleton/QG/102X/CMSSW_10_2_10/src

# To get conda
# Can't use $HOME because grid-control manges my env vars
source /afs/desy.de/user/a/aggleton/.bashrc
conda activate unfolding_py3

# Actually run all the unfolding jobs
set -x
# INPUTDIR and OUTPUTDIR come from the grid-control config
ODIR="${INPUTDIR}/${OUTPUTDIR}"

COMMON_ARGS="$INPUTDIR --outputDir=$ODIR --do$CHANNEL --angles $ANGLE --noBinnedPlots"

# MC input, all MC
# ------------------------------------------------------------------------------
#python unfolding.py ${COMMON_ARGS} --MCinput=True --MCsplit=False

# MC input, split (closure test), and alt response
# ------------------------------------------------------------------------------
#python unfolding.py ${COMMON_ARGS} --MCinput=True --MCsplit=True --useAltResponse=True

# MC input, not split - use to derive model, exp, pdf systs for data
# ------------------------------------------------------------------------------
# Do some per-variable/channel settings
EXP_SYST_AS_ALT_RSP=0
MERGE_LAST_PT_BIN=0
ZERO_BAD_RSP_BIN=0
if [[ "$CHANNEL" == *"Dijet"* ]]; then
    EXP_SYST_AS_ALT_RSP=1
else
    if [[ "ANGLE" == *"Multiplicity"* ]]; then
        MERGE_LAST_PT_BIN=1
        ZERO_BAD_RSP_BIN=1
    fi
fi

ARGS="${COMMON_ARGS} --MCinput=True "\
"--doExperimentalSysts=True --doExperimentalSystsAsAltResponse ${EXP_SYST_AS_ALT_RSP} "\
"--doScaleSysts=True --doPDFSysts=True "\
"--mergeLastPtBin ${MERGE_LAST_PT_BIN} --zeroBadResponseBins ${ZERO_BAD_RSP_BIN}"

# stash location of output file for use later with data
REFFILE=$(python unfolding_logistics.py out $ARGS)
python unfolding.py $ARGS

# Data input, use previous MC file to get systs if necessary
# ------------------------------------------------------------------------------
EXP_SYSTS_ARGS="--doExperimentalSystsFromFile=$REFFILE"
if [[ "$CHANNEL" == *"Dijet"* ]]; then
    EXP_SYSTS_ARGS="--doExperimentalSysts=True --doExperimentalSystsAsAltResponse=${EXP_SYST_AS_ALT_RSP}"
fi

ARGS="${COMMON_ARGS} "\
"${EXP_SYSTS_ARGS} "\
"--doScaleSystsFromFile=$REFFILE "\
"--doPDFSystsFromFile=$REFFILE "\
"--mergeLastPtBin ${MERGE_LAST_PT_BIN} --zeroBadResponseBins ${ZERO_BAD_RSP_BIN}"

python unfolding.py $ARGS

# Data input, use alt response matrix as cross-check
# ------------------------------------------------------------------------------
#python unfolding.py ${COMMON_ARGS} --useAltResponse=True

set +x
set +u
