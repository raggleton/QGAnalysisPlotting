#!/bin/bash -e

# Sync hadded files from NAF to here

function sync {
    # WORKDIR="workdir_ak4puppi_herwig"
    WORKDIR="$1"

    if [ ! -d "$WORKDIR" ]; then 
        echo "Creating $WORKDIR"
        mkdir "$WORKDIR"
    fi

    GEN=""
    if [[ "$WORKDIR" == *"_mgpythia"* ]]; then
        GEN="MGPythia"
    elif [[ "$WORKDIR" == *"_herwig"* ]]; then
        GEN="Herwig"
    elif [[ "$WORKDIR" == *"_pythiaOnlyFlat"* ]]; then
        GEN="PythiaOnlyFlat"
    elif [[ "$WORKDIR" == *"_powheg"* ]]; then
        GEN="Powheg"
    fi
    
    if [[ "$WORKDIR" != *"_pythia"* ]] && [[ "$WORKDIR" != *"_powheg"* ]]; then
        rsync -avzP NAF:"/nfs/dust/cms/user/aggleton/QG/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/Selection/$GEN/$WORKDIR/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" "$WORKDIR/"
    fi
    rsync -avzP NAF:"/nfs/dust/cms/user/aggleton/QG/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/Selection/$GEN/$WORKDIR/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" "$WORKDIR/"
}

for wdir in "$@";
do
    echo "Syncing $wdir"
    sync "$wdir"
done
