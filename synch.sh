#!/bin/bash -e

# Sync hadded files from NAF to here

function sync {
    # WORKDIR="workdir_ak4puppi_herwig"
    WORKDIR="$1"

    if [ ! -d "$WORKDIR" ]; then 
        echo "Creating $WORKDIR"
        mkdir "$WORKDIR"
    fi

    # Remember to change this
    GEN="MGPythia"
    GEN="Herwig"
    # GEN="PythiaOnlyFlat"

    rsync -avzP NAF:"/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/Selection/$GEN/$WORKDIR/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL_.root" "$WORKDIR/"
    rsync -avzP NAF:"/nfs/dust/cms/user/aggleton/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/Selection/$GEN/$WORKDIR/uhh2.AnalysisModuleRunner.MC.MC_QCD_.root" "$WORKDIR/"
}

for wdir in "$@";
do
    echo "Syncing $wdir"
    sync "$wdir"
done
