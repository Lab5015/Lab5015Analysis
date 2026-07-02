# Lab5015Analysis
collection of programs for the analysis of BTL sensor modules measurements with TOFHIR readout 

# Installation

## Login to a machine
```sh
ssh username@hostname
```

## Fresh installation of the analysis package
```sh
export MYNAME=putYourNameHere  #use your name for development
mkdir $MYNAME
cd $MYNAME
# please make sure to include the relevant branches, if needed
git clone --recursive https://github.com/Lab5015/Lab5015Analysis
cd Lab5015Analysis
source scripts/setup.sh
make
make exe
```

# Package structure
This package is structured as follows
- Utilities are organized in the `interface` and `src` folders:
    - `interface`: contains the class or functions headers
    - `src`: contains the class or functions implementation
- `main`: contains the main analysis code that make use of the Utilities
- `cfg`: contains the config files which are used to pass parameters to the executables
- `macros`: contains a collection of python codes that are used to further analyze the outcomes of the main codes

After the compilation, each step of the analysis can be executed from the main folder `Lab5015Analysis` with a command like:
`./bin/executable cfg/configfile`


*Every time you login* remember to source the setup script:
```
cd Lab5015Analysis
source scripts/setup.sh
```

*Before running the analysis*

The output filename and output location of each analysis step are defined in the cfg files. Before running the analysis make sure you have checked and if needed updated the relevant output paths in order not to overwrite the work of others.



# Run the legacy analysis (module characterization with time difference)
## Main workflow
The analysis of the collected data is structured in three steps
1. `moduleCharacterization_step1.cpp`:
   This step performs the first loop over the events and fills base histograms such as energy and ToT for each bar, each threshold and each over-voltage. A skim of the events according to the selections defined in the config file is also performed. The output is in the form of TTrees and histograms.

1. `moduleCharacterization_step2.cpp`:
   This step performs several loops over the events:
    1. loop over the histos filled in step1 and define the energy ranges
    1. loop over the skimmed TTrees and fill raw distribution of energy, ToT, energy ratio, ToT ratio according to the predefined bins
    1. loop over the skimmed TTrees and compute time walk corrections
    1. loop over the skimmed TTrees and apply the time walk corrections
    
   The output of this step are histograms.

1. `moduleCharacterizationSummaryPlots.py`:
   This step takes the output of step2 as an input and displays summary plots in a website. Loop over the events doesn't belong here.
   example:
   ```sh
   python moduleCharacterizationSummaryPlots.py -m 2 -i run6055 -o run6055
   ```

An additional code which is useful to plot the pulse shape for a given channel is `drawPulseShapeTB.exe`. 
example:
```sh
./bin/drawPulseShapeTB.exe cfg/drawPulseShapeTB.cfg
```

Remember to change the paths of input and output folders/files in all the configurations files as needed.


## Scripts to prepare configuration files
A set of scripts is available under the `scripts` folder to create configuration files (both for the moduleCharacterization executables and the drawPulseShapeTB.exe). 
The starting point are the base configuration files under the `cfg` folder
Edit the base configuration files (e.g. `moduleCharacterization_base.cfg`, `drawPulseShapeTB_base.cfg`, `minEnergies_base.txt`) and the `create_config.py` and `launch_create_config.sh` scripts, as needed.
To run the scripts:
```sh
cd scripts/
source launch_create_config.sh 
```

## Submit the analysis in parallel on lxplus Condor 
A set of scripts under the `scripts` folder allows the submission of multiple jobs (e.g. a set of runs corresponding to an overvoltage/threshold scan) in parallel singin Condor on lxplus. 
From `scripts` folder, edit `moduleCharacterizationCondor.sh`, `submit_moduleCharacterization_Condor.sub`, `list_cfg_moduleCharacterization.txt` as needed.

The script is used as follows:
```sh
condor_submit submit_moduleCharacterization_Condor.sub
```

# Integration of new methods for September 2025 TB
A new set of analysis methods has been introduced within the context of the analysis of September 2025 TB data. Please find the relevant reference material and the corresponding instructions in this twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MTDTB202509CERN




