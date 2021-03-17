#!/bin/bash

# Create pileup profile per trigger
# Thanks to Guillelmo

# See https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData
# https://twiki.cern.ch/twiki/bin/view/CMS/BrilcalcQuickStart
export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH

CERTFILE="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"

triggers=( "ZeroBias" "PFJet40" "PFJet60" "PFJet80" "PFJet140" "PFJet200" "PFJet260" "PFJet320" "PFJet400" "PFJet450" "AK8PFJet40" "AK8PFJet60" "AK8PFJet80" "AK8PFJet140" "AK8PFJet200" "AK8PFJet260" "AK8PFJet320" "AK8PFJet400" "AK8PFJet450" )
for TRIG in "${triggers[@]}"
do
    echo "Doing $TRIG..."
    brilcalc lumi --byls --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json  -i  "$CERTFILE" --hltpath HLT_${TRIG}_v* -o output_${TRIG}.csv
    /cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_2_16/bin/slc6_amd64_gcc700/pileupReCalc_HLTpaths.py -i output_${TRIG}.csv --inputLumiJSON pileup_latest.txt -o HLT_corrected_PileupJSON_${TRIG}.txt --runperiod Run2
    /cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_2_16/bin/slc6_amd64_gcc700/pileupCalc.py -i "$CERTFILE" --inputLumiJSON HLT_corrected_PileupJSON_${TRIG}.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100  MyDataPileupHistogram_${TRIG}.root
done
