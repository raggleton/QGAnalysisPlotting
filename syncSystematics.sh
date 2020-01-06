#!/usr/bin/env bash -e
#
# Sync systematics files from NAF

VARIATIONS=("Up" "Down")
SYSTS=("chargedHadronShift" "neutralHadronShift" "photonShift" "jecsmear_direction" "jersmear_direction" "pileup_direction")
# SYSTS=("pileup_direction")
# SYSTS=("jecsmear_direction" "jersmear_direction")

OUTPUTDIR="workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet/"
OUTPUTDIR="workdir_ak4puppi_data_jetAsymCut_pt1RecoConstituents_V11JEC_JER_target0p6_wta_groomed_fwdcenDijet_ptAveBinning/"
OUTPUTDIR="workdir_ak4puppi_data_jetAsymCut_pt1RecoConstituents_V11JEC_JER_target0p5_wta_groomed_fwdcenDijet_ptAveBinning/"

OUTPUTDIR="workdir_ak4puppi_data_jetAsymCut_pt1RecoConstituents_V11JEC_JER_target0p5_wta_groomed_fwdcenDijet_ZReweight/"
OUTPUTDIR="workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet_Zreweight_noZjet2Cut_zPt30_zJetAyms0p4/"
OUTPUTDIR="workdir_ak4puppi_data_trigBinningBetter2_jetAsymCut_pt1RecoConstituents_V11JEC_JER_tUnfoldBetter_target0p5_wta_groomed_fwdcenDijet_Zreweight_noZjet2Cut_newBinning"

NAFSTEM="/nfs/dust/cms/user/aggleton/QG/CMSSW_8_0_24_patch1/src/UHH2/QGAnalysis/Selection/"
STEMSRCDIR="${NAFSTEM}/MGPythia/workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p6_noZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_ptAveBinning"
STEMSRCDIR="${NAFSTEM}/MGPythia/workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_noZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_ptAveBinning"
STEMSRCDIR="${NAFSTEM}/MGPythia/workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut_zPt30"
STEMSRCDIR="${NAFSTEM}/MGPythia/workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut_zPt30_zJetAsym0p4"
STEMSRCDIR="${NAFSTEM}/MGPythia/workdir_ak4puppi_mgpythia_newFlav_jetAsymCut_chargedVars_pt1RecoConstituents_V11JEC_JER_target0p5_ZReweight_wta_groomed_fwdcenDijet_betterLargeWeightVeto_noZjet2Cut_zPt30_newBinning"

CMD="rsync -ahvzP"

for syst in ${SYSTS[@]}; do
    for var in ${VARIATIONS[@]}; do
        echo "------------------------------------------------------------------"
        echo ">>> ${syst}${var}"
        echo "------------------------------------------------------------------"
        ODIR="${OUTPUTDIR}/systematics_files/${syst}${var}"
        mkdir -p "${ODIR}"
        $CMD NAF:"${STEMSRCDIR}_${syst}${var}/uhh2.AnalysisModuleRunner.MC.MC_[QD]*.root" "${ODIR}/"
        # $CMD NAF:"${STEMSRCDIR}_${syst}${var}/uhh2.AnalysisModuleRunner.MC.MC_QCD.root" "${ODIR}/"
        # $CMD NAF:"${STEMSRCDIR}_${syst}${var}/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root" "${ODIR}/"
    done
done

# Just list these ones explicitly, cba to write a loop
SYSTS=(
    "ScaleVariationMuRNominal_ScaleVariationMuFUp"
    "ScaleVariationMuRNominal_ScaleVariationMuFDown"
    "ScaleVariationMuRUp_ScaleVariationMuFNominal"
    "ScaleVariationMuRDown_ScaleVariationMuFNominal"
    "ScaleVariationMuRUp_ScaleVariationMuFUp"
    "ScaleVariationMuRDown_ScaleVariationMuFDown"
)

for syst in ${SYSTS[@]}; do
    echo "------------------------------------------------------------------"
    echo ">>> ${syst}"
    echo "------------------------------------------------------------------"
    ODIR="${OUTPUTDIR}/systematics_files/${syst}"
    mkdir -p "${ODIR}"
    $CMD NAF:"${STEMSRCDIR}_${syst}/uhh2.AnalysisModuleRunner.MC.MC_[QD]*.root" "${ODIR}/"
    # $CMD NAF:"${STEMSRCDIR}_${syst}/uhh2.AnalysisModuleRunner.MC.MC_DYJetsToLL.root" "${ODIR}/"
done