#!/bin/bash 
#
# Do all the various combos of response comparison: PUS, jet R, generator

CMD="./print_jet_response_qg_plots.py"

# PUS
for R in 4 8;
do
    $CMD workdir_ak${R}puppi_mgpythia workdir_ak${R}chs_mgpythia --output workdir_ak${R}puppi_mgpythia/response_plots_comparePUS
done

# jet R
for PUS in chs puppi;
do
    $CMD workdir_ak4${PUS}_mgpythia workdir_ak8${PUS}_mgpythia --output workdir_ak4${PUS}_mgpythia/response_plots_compareR
done

# generator
for R in 4 8; 
do
    for PUS in chs puppi;
    do
        $CMD workdir_ak${R}${PUS}_mgpythia workdir_ak${R}${PUS}_herwig_reweight --output workdir_ak${R}${PUS}_herwig_reweight/response_plots_compareGenerator
    done
done