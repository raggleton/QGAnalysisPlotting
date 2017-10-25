#!/bin/bash 
#
# Do all the various combos of response comparison: PUS, jet R, generator

CMD="./print_jet_response_qg_plots.py"

function doComparisons {
    # PUS
    for R in 4 8;
    do
        $CMD --comparison workdir_ak${R}puppi_mgpythia workdir_ak${R}chs_mgpythia --output workdir_ak${R}puppi_mgpythia/response_plots_comparePUS
        $CMD --comparison workdir_ak${R}puppi_herwig_reweight workdir_ak${R}chs_herwig_reweight --output workdir_ak${R}puppi_herwig_reweight/response_plots_comparePUS
    done

    # jet R
    for PUS in chs puppi;
    do
        $CMD --comparison workdir_ak4${PUS}_mgpythia workdir_ak8${PUS}_mgpythia --output workdir_ak4${PUS}_mgpythia/response_plots_compareR
        $CMD --comparison workdir_ak4${PUS}_herwig_reweight workdir_ak8${PUS}_herwig_reweight --output workdir_ak4${PUS}_herwig_reweight/response_plots_compareR
    done

    # generator
    for R in 4 8; 
    do
        for PUS in chs puppi;
        do
            $CMD --comparison workdir_ak${R}${PUS}_mgpythia workdir_ak${R}${PUS}_herwig_reweight --output workdir_ak${R}${PUS}_herwig_reweight/response_plots_compareGenerator
        done
    done
}

function doIndividuals {
    # individual ones
    for R in 4 8; 
    do
        for PUS in chs puppi;
        do
            for GEN in mgpythia herwig_reweight;
            do
                $CMD workdir_ak${R}${PUS}_${GEN}
            done
        done
    done
}

# doComparisons;
doIndividuals;