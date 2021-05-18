# Quark-gluon analysis & plotting code

This code does further analysis & plotting on UHH2 Ntuples produced by https://github.com/raggleton/QGAnalysis

## Contents

  * [Ntuple organisation](#ntuple-organisation)
  * [Installation](#installation)
  * [Running](#running)
    + [Constructing binning schemes & plots](#constructing-binning-schemes---plots)
    + [Unfolding](#unfolding)
      - [Running unfolding in batch jobs on BIRD](#running-unfolding-in-batch-jobs-on-bird)
    + [Plotting the results of unfolding](#plotting-the-results-of-unfolding)
    + [Doing bottom-line tests](#doing-bottom-line-tests)
    + [Calculation of summary stats for unfolding & RIVET files](#calculation-of-summary-stats-for-unfolding---rivet-files)
    + [Summary plots](#summary-plots)
  * [Other scripts especially pertaining to paper plots](#other-scripts-especially-pertaining-to-paper-plots)
    + [Jet pT plot](#jet-pt-plot)
    + [Flavour fraction plot](#flavour-fraction-plot)
    + [Detector-unfolded plots](#detector-unfolded-plots)



## Ntuple organisation

All scripts assume the following structure containing the UHH2 ntuples:

- Directory for AK4 ntuples
    - `uhh2.AnalysisModuleRunner.DATA.Data_JetHT_ZeroBias.root`
    - `uhh2.AnalysisModuleRunner.DATA.Data_SingleMu.root`
    - `uhh2.AnalysisModuleRunner.MC.MC_QCD.root`
    - `uhh2.AnalysisModuleRunner.MC.MC_DY.root`
    - ... other data & MC nominal ntuples ...
    - `systematics_files`
        - `jes_directionUp`
            - `uhh2.AnalysisModuleRunner.MC.MC_QCD.root`
            - `uhh2.AnalysisModuleRunner.MC.MC_DY.root`
        - ... other systematics dirs
- Directory for AK8 ntuples
    - ... similar structure as for AK4 dir

Throughout the scripts, the ntuple directories might be referred to as "workdirs"

## Installation

Scripts are written for python 3, but may work on python2 by accident.
Most of this code relies on `ROOT` (especially `PyROOT` & `TUnfold`). Note that from version 6.22, changes to PyROOT were made that are incompatible with the code here.

:warning: **You must use ROOT < 6.22!** :warning:

See https://root.cern/manual/python/#new-pyroot-backwards-incompatible-changes for details (this code uses `ROOT.Double`, and there may be other changes I am unaware of).

I used ROOT 6.20.2, which uses TUnfold 17.6 - analysis unverified using newer TUnfold versions!

It also requires several python packages: `numpy`, `scipy`, `pandas` - these can generally be found in any vaguely recent CMSSW setup (>= 10_2_X).
Some of it also relies on less common packages, including the `uproot` (the older uproot3, not the more recent 4, but **not** via the `uproot3` package), `uncertainties`, `jax`, and `rootpy` packages.

If you need python 3 and/or ROOT, [miniconda](https://docs.conda.io/en/latest/miniconda.html) is my choice (especially on the NAF), since it requires no special permissions.

You can create a copy of my conda environment called `qganalysis` using:

```conda env create -f conda_environment.yml```

You can then activate with the usual `conda activate qganalysis`, etc.

However, YODA cannot be installed easily via conda, and there are API changes that clash with python3 (1.7.7 has different API to 1.8.X, e.g. `.bins` vs `.bins()`, but the former doesn't compile under Python3.7 TODO: update my usage of API). 
For the scripts that need YODA, I just do them within CMSSW 10_6_X (since they often do not need other packages outside of CMSSW), to convert them into a more portable format.

If you really want to install it:

To install YODA manually, into a directory `PREFIX` (must be full filepath) of your choosing:


```
conda activate qganalysis
pip install -U Cython
wget --no-check-certificate http://www.hepforge.org/archive/yoda/YODA-1.7.7.tar.gz
tar xzf YODA-1.7.7.tar.gz
cd YODA-1.7.7
cd pyext/yoda
cython ./core.pyx --cplus -I . -I ./include  -I . -I ./include -o core.cpp
cython ./util.pyx --cplus -I . -I ./include  -I . -I ./include -o util.cpp
cython ./rootcompat.pyx --cplus -I . -I ./include  -I . -I ./include -o rootcompat.cpp
./configure --prefix=$PREFIX
make -j4
make -j4 install
cp yodaenv.sh $PREFIX/yodaenv.sh
cp yodaenv.sh ${CONDA_PREFIX}/etc/conda/activate.d/
```

By copying to the `activate.d`, `yodaenv.sh` will be called whenever the conda environment is activated.

If you already have python3 & ROOT and just want the other packages, you can instead use the `requirements.txt` to install these via `pip install -r requirements.txt`.

## Running

Most of the scripts are relatively standalone, operating only on the Ntuples directly.

Scripts with options will give full info about them by using the `-h` or `--help` flags.

### Constructing binning schemes & plots

Use [`determine_lambda_binning.py uhh2.AnalysisModuelRunner.MC.MC_QCD.root`](determine_lambda_binning.py)

`--target XX` to set a desired purity & stability fraction (e.g. 0.5)

Also options to rebin other files using the derived binning scheme.

Binning scheme saved as txt file. To use it, it requires copying manually into `VAR_UNFOLD_DICT_TARGET0p5` in [`qg_common.py`](qg_common.py). It also requires copying into the vectors in the `Binning` namespace in [`QGAddModules.h`](https://github.com/raggleton/QGAnalysis/blob/master/include/QGAddModules.h#L651)

### Unfolding

The ntuples already contain the binned response matrices and 1D flattened distributions.
If you change the binning scheme, you must reprocess the ntuples first. 

Unfolding is done by the [`unfolding.py`](unfolding.py) script.
This relies on TUnfold, specifically v17.6 or newer.

```
unfolding.py <ak4 or ak8 workdir> <options> <regions> --angles <angles>
```

`<regions>` are one or more of: `--doDijetCentral`, `--doDijetCentralGroomed`, `--doDijetForward`, `--doDijetForwardGroomed`, `--doZPJ`, `--doZPJGroomed`. To do all regions, you can instead use `--doAllRegions`.

`<angles>` are one or more of: `jet_LHA`, `jet_width`, `jet_thrust`, `jet_puppiMultiplicity`, `jet_pTD`. Each also has a `*_charged` version (e.g. `jet_LHA_charged`). To do all angles, you can instead specify `all`.

There are a wide variety of options when unfolding, the most important being:

- `--MCinput 0/1` : 0 for unfolding data, 1 for unfolding MC
- `--MCsplit 1`: 1 for unfolding split sample MC (20% being unfolded, 80% in the response matrix)
- `--doExperimentalSysts 1`: for doing experimental systs, including shower & hadronization (i.e. Herwig++) syst
- `--doScalesSysts 1`: for doing scale variation uncertainty
- `--doPDFSysts 1`: for doing PDF uncertainty - warning SLOW as it unfolds 100 times

By default, the script will also produce lots of plots. This can be slow, so to disable it, use `--noBinnedPlots`. Plots can instead be made using the script below.

The output from unfolding is organised as follows:
- There is a "top-level" dir, `unfolding_regularizeNone_DATA_subFakes_densityModeBinWidth_constraintNone_experimentalSystFromFile_scalSystFromFile_pdfSystFromFile_signalRegionOnly`. This is by default put within the workdir specified when running the script. If you use the `--outputDir` option, it will instead be places within that directory. 
- This top-level directory then has separate directories for each region, and within that, separate directories for each angle
- The top-level directory is designed to reflect most of the options used when running. 
- Within each (region,angle) directory, you will find:
    - `args.json` which has a copy of all the options used when running
    - `unfolding_slim.root`, which only has the main normalised distribution hists
    - `unfolding_result.pkl`, a python pickle file, which has many individual distributions, e.g. the response matrix & unfolded results for each PDF variation
    - Many plots
    - Systematic variations will often have subdirectories with corresponding names

#### Running unfolding in batch jobs on BIRD

Doing the unfolding for all angles/regions/systematics will take a very long time. Instead, do it in parallel jobs on the BIRD cluster.
Jobs are handled using [grid-control](https://grid-control.github.io/).
This assumes you are using a conda environment.

The scripts [`gc_submitUnfoldingJobs_ak4.conf`](gc_submitUnfoldingJobs_ak4.conf) and [gc_submitUnfoldingJobs_ak8.conf](`gc_submitUnfoldingJobs_ak8.conf`) are the job submit files for AK4 and AK8 jobs, respectively. Each runs the commands in [`runUnfoldingJob.sh`](runUnfoldingJob.sh)
You should update `runUnfoldingJob.sh` to match your conda environment name, and any other dir locations. 

Note that you cannot run grid-control with python3, so for doing batch unfolding plots, you'll need to `conda deactivate` first.

### Plotting the results of unfolding

Use [`do_unfolding_plots.py <unfolding dir>`](do_unfolding_plots.py).

`<unfolding dir>` here is the top-level output directory from `unfolding.py` (e.g. in the example above, `<workdir>/unfolding_regularizeNone_DATA_subFakes_densityModeBinWidth_constraintNone_experimentalSystFromFile_scalSystFromFile_pdfSystFromFile_signalRegionOnly`)

Again, this uses the same region & angle options as the unfolding script.
There are also various options about which set(s) of plots to produce:

- `--doBinnedPlotsGenPt`: plots per generator pt bin
- `--doBinnedPlotsGenLambda`: plots per generator lambda bin
- `--doBinnedPlotsRecoPt`: plots per detector pt bin
- `--doBigNormed1DPlots`: do plots with normalised plots, across all pt bins
- `--doBigAbs1DPlots`: do plots with absolute plots, across all pt bins

Note that there is an optional `--paper` argument for formatting for the final publication (i.e. removing "Preliminary").

### Doing bottom-line tests

Use [`print_bottom_line_test.py`](print_bottom_line_test.py).
Again, this uses the same region & angle options as the unfolding script.

You might also want the `--noNullBins` options if the VV^-1 plots look bad (i.e. not unity).

### Calculation of summary stats for unfolding & RIVET files

Before we can plot the grand summary plots (e.g. mean vs pT), we first extract the summary stats (mean, RMS, etc), from the required inputs (unfolding results, and YODA files from RIVET), and save them as pandas dataframes, using the HDF5 output format.

**To convert unfolding results:**

`extract_unfolding_summary_stats.py --ak4source <ak4 unfolding dir> --ak8source <ak8 unfolding dir> --h5output <output H5 filename>`

**To convert RIVET results:**

This is more complicated, since there are many RIVET inputs corresponding to different generators. For each generator, we specify:

- the dijet YODA file (`--yodaInputDijet`)
- the Z+Jet YODA file (`--yodaInputZPJ`)
- the label (`--yodaInputLabel`)

NB the label here must be one of those keys specified in `SAMPLE_STYLE_DICTS` in `do_summary_unfolding_plots.py`.

Note that we can also specify the unfolded results, since we can then calculate the delta metric.

```
extract_rivet_summary_stats.py --ak4source <ak4 unfolding dir> --ak8source <ak8 unfolding dir> \
--h5output <output H5 filename> \
# for each generator, you need entries:
--yodaInputDijet <.../QCD.yoda> --yodaInputZPJ <.../ZPJ.yoda> --yodaLabel "MC GEN TUNE1"
```

### Summary plots

Use `do_summary_unfolding_plots.py` with arguments:

- `--h5input <unfolding.h5>`: the HDF5 file made by `extract_unfolding_summary_stats.py`
- `--h5inputRivet <rivet.h5>`: the HDF5 file made by `extract_rivet_summary_stats.py`
- `--outputDir`: directory to put plots

By default, the script will plot all inputs (data, nominal MCs, RIVET files). You can also specify which sets of inputs to plot:

- `--onlyYodaData`: only data & RIVET files
- `--onlyDataOldMC`: only plot data, MG+Pythia8, Herwig++

The script makes 2 types of summary plots: metric (e.g. mean) vs pt plots across all regions/angles; and the 5-column 5-point grand summary plots with choice bins shown.

- `--doMetricVsPt`: do metric vs pt plots
- `--doSummaryBins`: do grand summary plots

Note that there are also an optional `--paper` argument for formatting for the final publication (i.e. removing "Preliminary"), and an optional `--supplementary` argument for supplementary plots.

## Other scripts especially pertaining to paper plots

Note that these have an optional `--paper` argument for formatting for the final publication (i.e. removing "Preliminary").

### Jet pT plot

Use [`print_data_mc_jet_pt.py <workdir>`](print_data_mc_jet_pt.py), where `<workdir>` is a directory containing the UHH2 ntuples. By default, it includes systematics combined into one band.

Plots produced in `<workdir>/data_mc_jet_pt`

### Flavour fraction plot

Use [`print_flav_fraction_comparison_plots.py`](print_flav_fraction_comparison_plots.py), with arguments:

- `--dj` to do dijet region plots
- `--zpj` to do Z+Jet region plots
- `--gen` to use on genjet flavour (otherwise uses reco jet flavour)

This script can plot multiple UHH2 ntuples, as well as multiple RIVET files.

For each UHH2 ntuple, you must specify:

- `--dir`: the dir that has the ntuples - it will only use the files `uhh2.AnalysisModuleRunner.MC.MC_QCD.root`, `uhh2.AnalysisModuleRunner.MC.MC_DY.root`, and `uhh2.AnalysisModuleRunner.MC.MC_HERWIG_QCD.root`
- `--dirLabel`: the label to be attached to those ntuples (e.g. "MG+PYTHIA8")

For each RIVET input, you must specify:

- `--rivetdj`: dijet YODA file
- `--rivetdjLabel`: label to go with the dijet file
- `--rivetzpj`: Z+Jet YODA file
- `--rivetzpjLabel`: label to go with the Z+Jet file

### Detector & unfolded plots

Use [`print_detector_unfolded_plots.py <unfolding_dir>`](`print_detector_unfolded_plots.py`), where `<unfolding_dir>` is the top-level directory from `unfolding.py` (see description above).

You can specify which angles using the same options as for `unfolding.py` (e.g. "jet_LHA jet_pTD_charged", or "all")
