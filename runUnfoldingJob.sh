#!/bin/bash -e

function generateNickname {
    CHANNEL=$1
    MCINPUT=$2
    MCSPLIT=$3
    SYST=$4
    ANGLE=$5

    MCSTR=""
    if [[ $MCINPUT == "True" ]]
    then 
        MCSTR="_MC"
        if [[ $MCSPLIT == "True" ]]
        then 
            MCSTR="${MCSTR}_split" 
        else
            MCSTR="${MCSTR}_all"
        fi
    fi
    SYSTSTR=""
    if [[ $SYST == "True" ]]; then SYSTSTR="_syst"; fi

    echo "${CHANNEL}${MCSTR}${SYSTSTR}_densityModeBinWidth_constraintArea_${ANGLE}"
}

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /nfs/dust/cms/user/aggleton/QG/CMSSW_8_0_24_patch1/src
eval `scramv1 runtime -sh`
cd /nfs/dust/cms/user/aggleton/QG/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/Analysis/QGAnalysisPlotting
printenv | sort
# Run both unregularized and regularized
set -x
ODIR="unfolding_output"
# run with systematics
SYST="True"
NICKNAME=$(generateNickname $CHANNEL $MCINPUT $MCSPLIT $SYST $ANGLE)

OUTPUTDIR1="${ODIR}/unfolding_regNone_${NICKNAME}"
python unfolding.py $INPUTDIR --do$CHANNEL --regularize=None --MCinput=$MCINPUT --MCsplit=$MCSPLIT --doSyst=$SYST --outputDir=$OUTPUTDIR1 --angles $ANGLE

OUTPUTDIR2="${ODIR}/unfolding_regTau_${NICKNAME}"
python unfolding.py $INPUTDIR --do$CHANNEL --regularize=tau --MCinput=$MCINPUT --MCsplit=$MCSPLIT --doSyst=$SYST --outputDir=$OUTPUTDIR2 --angles $ANGLE

# repeat without any systematics
SYST="False"
NICKNAME=$(generateNickname $CHANNEL $MCINPUT $MCSPLIT $SYST $ANGLE)

OUTPUTDIR3="${ODIR}/unfolding_regNone_${NICKNAME}"
python unfolding.py $INPUTDIR --do$CHANNEL --regularize=None --MCinput=$MCINPUT --MCsplit=$MCSPLIT --doSyst=$SYST --outputDir=$OUTPUTDIR3 --angles $ANGLE

OUTPUTDIR4="${ODIR}/unfolding_regTau_${NICKNAME}"
python unfolding.py $INPUTDIR --do$CHANNEL --regularize=tau --MCinput=$MCINPUT --MCsplit=$MCSPLIT --doSyst=$SYST --outputDir=$OUTPUTDIR4 --angles $ANGLE

# Zip up output files
cd $ODIR && tar cvzf unfolding_${NICKNAME}.tar.gz ${OUTPUTDIR1#$ODIR/} ${OUTPUTDIR2#$ODIR/} ${OUTPUTDIR3#$ODIR/} ${OUTPUTDIR4#$ODIR/} && sleep 10 && cd ..
set +x
