# Nutples & unfolded files:


Ntuple dirs:

```
# AK4:
NTUPLEAK4=/nfs/dust/cms/user/aggleton/QG/102X/CMSSW_10_2_16/src/UHH2/QGAnalysis/Selection/workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning6_PUDownX3
# AK8
NTUPLEAK8/nfs/dust/cms/user/aggleton/QG/102X/CMSSW_10_2_16/src/UHH2/QGAnalysis/Selection/workdir_102X_v3data_v2mc_ak8puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning6_PUDownX3
```

NB mostly symlinks to 
```
/nfs/dust/cms/user/aggleton/QG/102X/CMSSW_10_2_16/src/UHH2/QGAnalysis/Selection/workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning6
/nfs/dust/cms/user/aggleton/QG/102X/CMSSW_10_2_16/src/UHH2/QGAnalysis/Selection/workdir_102X_v3data_v2mc_ak8puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning6
```

Unfolded files:

```
# AK4:
UNFOLDAK4=/nfs/dust/cms/user/aggleton/QG/102X/CMSSW_10_2_16/src/UHH2/QGAnalysis/Selection/workdir_102X_v3data_v2mc_ak4puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning6_PUDownX3/unfolding_output/unfolding_regularizeNone_DATA_subFakes_densityModeBinWidth_constraintNone_dijet_experimentalSystAsAltResponse_scaleSyst_pdfSyst_zpj_ExperimentalSystFromFile_scalSystFromFile_pdfSystFromFile_mergeLastPtBin_zeroBadResponseBins_onlyZPJmulti_signalRegionOnly
# AK8:
UNFOLDAK8=/nfs/dust/cms/user/aggleton/QG/102X/CMSSW_10_2_16/src/UHH2/QGAnalysis/Selection/workdir_102X_v3data_v2mc_ak8puppi_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning6_PUDownX3/unfolding_output/unfolding_regularizeNone_DATA_subFakes_densityModeBinWidth_constraintNone_dijet_experimentalSystAsAltResponse_scaleSyst_pdfSyst_zpj_ExperimentalSystFromFile_scalSystFromFile_pdfSystFromFile_mergeLastPtBin_zeroBadResponseBins_onlyZPJmulti_signalRegionOnly
```

## RIVET files

I am assuming you have cloned `SubstructureProfessor` inside this dir. Change paths as necessary.

# Plots:

## Jet pT:

`./print_data_mc_jet_pt.py --paper ${NTUPLEAK4}`

Plots will be in: `${NTUPLEAK4}/data_mc_jet_pt`

## Flav fractions:

Main paper plot:

`./print_flav_fraction_comparison_plots.py --dj --zpj --gen --paper --dir ${NTUPLEAK4} --dirLabel "MG+Pythia8"`

Produces: `${NTUPLEAK4}/flav_fractions_gen/g_flav_fraction_compare_jet1_paper.pdf`


To make the supplementary plots with RIVET inputs, first run `yoda2root` on the yoda files: 

`find SubstructureProfessor/Julian/yoda_files/ -name "*Mar16.yoda" -execdir yoda2root {} \;`

Supplementary plot with RIVET inputs:

```
./print_flav_fraction_comparison_plots.py --dir ${NTUPLEAK4} --dirLabel "Nominal" \
--dj --zpj --gen  --paper \
--rivetdj SubstructureProfessor/Julian/yoda_files/QCD_Pythia8_CP5_Mar16.root \
--rivetdjLabel "Pythia8 CP5" \
--rivetzpj SubstructureProfessor/Julian/yoda_files/DY_Pythia8_CP5_Mar16.root \
--rivetzpjLabel "Pythia8 CP5" \
--rivetdj SubstructureProfessor/Julian/yoda_files/QCD_Herwig7_CH3_Mar16.root \
--rivetdjLabel "Herwig7 CH3" \
--rivetzpj SubstructureProfessor/Julian/yoda_files/DY_Herwig7_CH3_Mar16.root \
--rivetzpjLabel "Herwig7 CH3" \
--rivetdj SubstructureProfessor/Julian/yoda_files/QCD_Sherpa_LO_Mar16.root \
--rivetdjLabel "Sherpa LO"  \
--rivetzpj SubstructureProfessor/Julian/yoda_files/DY_Sherpa_LO_Mar16.root \
--rivetzpjLabel "Sherpa LO"  \
--rivetdj SubstructureProfessor/Julian/yoda_files/QCD_Sherpa_LOM_Mar16.root \
--rivetdjLabel "Sherpa LO+jet" \
--rivetzpj SubstructureProfessor/Julian/yoda_files/DY_Sherpa_LOM_Mar16.root \
--rivetzpjLabel "Sherpa LO+jet" \
--rivetdj SubstructureProfessor/Julian/yoda_files/QCD_Sherpa_NLO_Mar16.root \
--rivetdjLabel "Sherpa NLO" \
--rivetzpj SubstructureProfessor/Julian/yoda_files/DY_Sherpa_NLO_Mar16.root \
--rivetzpjLabel "Sherpa NLO"
```

Produces:

## Unfolded vs detector plots

```
./print_detector_unfolded_plots.py $UNFOLDAK4 --paper --angles jet_LHA jet_pTD
```

Produces plots in `${UNFOLDAK4}/detector_unfolded_dijet_zpj/`

## Unfolded & systematics plots

```
./do_unfolding_plots.py $UNFOLDAK4 --doBinnedPlotsGenPt --doBinnedPlots --paper --onlyPaperPlots --doDijetCentral --doZPJ --doDijetCentralGroomed --doZPJGroomed --angles all 
```

## Summary plots

Make the HDF5 files:
```
SUMMARYDIR=/nfs/dust/cms/user/aggleton/QG/102X/CMSSW_10_2_16/src/UHH2/QGAnalysis/Selection/summary_102X_v3data_v2mc_fixSelCutOrder_puppiJER_tightJetId_constitPt0MultPt1_WeightCuts_zjAsym_genjetGhostFlav_noBCpref_genJetNoMu_fixCharged_jetIDSel_newBinning6_PUDownX3_AsAltReponseDijet_MergeLastPtBinZeroBadResponseBinsZPJMultiplicity
mkdir $SUMMARYDIR
UNFOLDH5="${SUMMARYDIR}/unfolding.h5"
./extract_unfolding_summary_stats.py --ak4source $UNFOLDAK4 --ak8source $UNFOLDAK8 --h5output $UNFOLDH5
```

> NB if you can't install YODA, do this command under CMSSW_10_6_X to use YODA. Once you have the output H5 file, you can work outside CMSSW again

```
RIVETH5="${SUMMARYDIR}/rivet_new.h5"
./extract_rivet_summary_stats.py \
--yodaInputDijet SubstructureProfessor/Julian/yoda_files/QCD_Pythia8_CP5_Mar16.yoda  \
--yodaInputZPJ SubstructureProfessor/Julian/yoda_files/DY_Pythia8_CP5_Mar16.yoda  \
--yodaLabel "Pythia8 CP5"  \
--yodaInputDijet SubstructureProfessor/Julian/yoda_files/QCD_Pythia8_CP2_Mar16.yoda  \
--yodaInputZPJ SubstructureProfessor/Julian/yoda_files/DY_Pythia8_CP2_Mar16.yoda  \
--yodaLabel "Pythia8 CP2"  \
--yodaInputDijet SubstructureProfessor/Julian/yoda_files/QCD_Herwig7_CH3_Mar16.yoda  \
--yodaInputZPJ SubstructureProfessor/Julian/yoda_files/DY_Herwig7_CH3_Mar16.yoda  \
--yodaLabel "Herwig7 CH3"  \
--yodaInputDijet SubstructureProfessor/Julian/yoda_files/QCD_Sherpa_LO_Mar16.yoda   \
--yodaInputZPJ SubstructureProfessor/Julian/yoda_files/DY_Sherpa_LO_Mar16.yoda  \
--yodaLabel "Sherpa LO"  \
--yodaInputDijet SubstructureProfessor/Julian/yoda_files/QCD_Sherpa_LOM_Mar16.yoda   \
--yodaInputZPJ SubstructureProfessor/Julian/yoda_files/DY_Sherpa_LOM_Mar16.yoda  \
--yodaLabel "Sherpa LO+jet"  \
--h5output "$RIVETH5"
```

Now make the plots:

```
./do_summary_unfolding_plots.py --h5input "$UNFOLDH5" --h5inputRivet "$RIVETH5" --onlyYodaData --doMetricVsPt --paper
./do_summary_unfolding_plots.py --h5input "$UNFOLDH5" --h5inputRivet "$RIVETH5" --onlyYodaData --doMetricVsPt
./do_summary_unfolding_plots.py --h5input "$UNFOLDH5" --h5inputRivet "$RIVETH5" --onlyDataOldMC --doMetricVsPt --paper
./do_summary_unfolding_plots.py --h5input "$UNFOLDH5" --h5inputRivet "$RIVETH5" --onlyDataOldMC --doMetricVsPt
./do_summary_unfolding_plots.py --h5input "$UNFOLDH5" --h5inputRivet "$RIVETH5" --doMetricVsPt --paper --supplementary
```

Plots are in `$SUMMARYDIR`

# Copying plots to AN/paper

Use script `Figures/copyPlots.sh` in the AN repo, which copies plots to AN/Figures. Then use `copyFigureFilesToPaperRepo.sh` to copy them to the paper repo (check paths in script first!)

NB `copyPlots.sh` doesn't copy flavour fraction plots, you have to do these manually, sorry.
