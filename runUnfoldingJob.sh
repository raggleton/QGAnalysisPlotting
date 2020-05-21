#!/bin/bash -e

# Run all the unfolding jobs for a given signal region, for a given angle
# Script gets given the following parameters:
# parameters   = CHANNEL ANGLE

# Need this as not empty at logon
export LD_LIBRARY_PATH=""

# Setup CMSSW etc
CMSSW_AREA=/nfs/dust/cms/user/aggleton/QG/102X/CMSSW_10_2_10/src
# cmssw setup triggers this err so turn it off for now
#set +u
#source /cvmfs/cms.cern.ch/cmsset_default.sh
#cd ${CMSSW_AREA}
#eval `scramv1 runtime -sh`
#cd ${CMSSW_AREA}/UHH2/QGAnalysisPlotting/
#printenv | sort
#set -u
# To get conda
# Can't use $HOME because grid-control manges my env vars
source /afs/desy.de/user/a/aggleton/.bashrc
conda activate unfolding_py3

# Actually run all the unfolding jobs
set -x
ODIR="${INPUTDIR}/unfolding_output"

# MC input, split (closure test), and alt response
python unfolding.py $INPUTDIR --do$CHANNEL --regularize=None --MCinput=True --MCsplit=True --useAltResponse=True --angles $ANGLE --noBinnedPlots --outputDir=$ODIR

# MC input, not split - use to derive model, exp, pdf systs for data
ARGS="$INPUTDIR --outputDir=$ODIR --do$CHANNEL --regularize=None --MCinput=True --MCsplit=False --useAltResponse=False --doExperimentalSysts=True --doScaleSysts=True --doPDFSysts=True --angles $ANGLE --noBinnedPlots"
# stash location of output file for use later with data
REFFILE=$(python unfolding_logistics.py out $ARGS)
python unfolding.py $ARGS

# Data input, use previous MC file to get systs
python unfolding.py $INPUTDIR --outputDir=$ODIR --do$CHANNEL --regularize=None --MCinput=False --MCsplit=False --useAltResponse=False --doExperimentalSystsFromFile=$REFFILE --doScaleSystsFromFile=$REFFILE --doPDFSystsFromFile=$REFFILE --angles $ANGLE --noBinnedPlots

# Data input, use alt response matrix as cross-check
python unfolding.py --outputDir=$ODIR $INPUTDIR --do$CHANNEL --regularize=None --MCinput=False --MCsplit=False --useAltResponse=True --angles $ANGLE --noBinnedPlots

# Zip up output files
# Not needed?
#cd $ODIR && tar cvzf unfolding_${NICKNAME}.tar.gz ${OUTPUTDIR1#$ODIR/} ${OUTPUTDIR2#$ODIR/} ${OUTPUTDIR3#$ODIR/} ${OUTPUTDIR4#$ODIR/} && sleep 10 && cd ..
#cd $ODIR && tar cvzf unfolding_${CHANNEL}_${MCINPUT}_${MCSPLIT}_${ANGLE}.tar.gz $ODIR/* && sleep 10 && cd ..
set +x
set +u

