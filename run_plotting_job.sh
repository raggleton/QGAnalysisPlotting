#!/bin/bash


# Need this as not empty at logon
export LD_LIBRARY_PATH=""

# To get conda
# Can't use $HOME because grid-control manges my env vars
source /afs/desy.de/user/a/aggleton/.bashrc
conda activate unfolding_py3

./do_unfolding_plots.py $INPUTDIR --doBinnedPlotsGenPt --do$CHANNEL --angles all

