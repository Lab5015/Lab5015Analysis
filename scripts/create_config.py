#!/usr/bin/env python
import shutil
from pathlib import Path
import argparse
import sys
from channelMapping import *

# ------ MODIFY YOUR PATHS HERE -------
Lab5015_path = "/eos/home-s/spalluot/MTD/TB_CERN_Sep25/Lab5015Analysis"
plot_path = "/eos/home-s/spalluot/www/MTD/MTDTB_CERN_Sep25"
# -------------------------------------


# parser
# ----------------------------
parser = argparse.ArgumentParser(description='This script creates moduleCharacterization , drawPulseShape, and minEnergy cfg')
parser.add_argument("-ml", "--modulelabel", required=True, type=str, help="module label")
parser.add_argument("-r", "--runs", required=True, type=str, help="comma-separated list of runs to be processed")
parser.add_argument("-t", "--temperature", required=True, type=str, help="temperature")
parser.add_argument("-ov", "--Vov", required=True, type=str, help="overvoltage")
parser.add_argument("-c", "--config", required=True, type=str, help="config number")
parser.add_argument("-th", "--threshold", required=True, type=str, help="threshold used for the scan")
parser.add_argument("-vth1", "--thresholdT1", required=False, default="-1", type=str, help="specific threshold for T1")
parser.add_argument("-vth2", "--thresholdT2", required=False, default="-1", type=str, help="specific threshold for T2")
parser.add_argument("-vthe", "--thresholdE",  required=False, default="-1", type=str, help="specific threshold for E")
parser.add_argument("-e", "--extraLabel", required=False,  type=str,  default="", help="eg: angle or check or whatever")
parser.add_argument("--whichEnergyIntercalib", required=False,  default="", type=str, help="i.e. TOFHIR, TOFHIR_LO or empty for None")
parser.add_argument("--dutASIC", required=False, default=7, type=int, help="the DUT ASIC position (default set to 7)")
parser.add_argument("--refASIC", required=False, default=4, type=int, help="the REF ASIC position (default set to 4)")
parser.add_argument("--refBar",  required=False, type=int, help="the REF bar on which ask for coincidence (default set to 7 in the code)")
parser.add_argument("--calibBar",  required=False, type=int, help="the REF bar on which coincidence was required to derive calib factors (default = refBar)")
parser.add_argument("--saveRefInfoFlag", required=False, type=int, default=0, help="0 or 1: flag to set if reference info should be saved")
parser.add_argument("--refCalibPath",    required=False, type=str, default="DM_348_Vov3.00_T18C_conf87_refBar9", help="REF module label used to get the energy intercalibration factors (default: DM_348_Vov3.00_T18C_conf87_refBar9)")
parser.add_argument("--refAmpWalkPerBar", required=False, default=0, type=bool, help="do you want to compute REF amp walk correction per DUT bar? (default=0)")
args = parser.parse_args()

if args.whichEnergyIntercalib != "":
    print("\nThe assumption is that the energy intercalibration file is of the kind <DM_label><refBar>, so the <extra_label> is not used to retrieve the calibration files")

# checks on ASICS : DUT e REF ASICs cannot be the same
# -----------------------------------------------------
if args.refASIC == args.dutASIC:
    print("[ERROR] DUT and REF ASICs cannot be the same")
    sys.exit(1)

# default ref bar
# ------------------------------
if args.refBar is None:
    args.refBar = 7
else:
    # - if the reference bar is not the default one, add an extra label
    if args.extraLabel:
        args.extraLabel += "_refBar"+str(args.refBar)
    else:
        args.extraLabel = "refBar"+str(args.refBar)

if args.calibBar is None:
    args.calibBar = args.refBar
    
# base label
# ----------------------------
Vov_float = float(args.Vov)
base_label = f"{args.modulelabel}_Vov{Vov_float:.2f}_T{args.temperature}C"

# channels
# -------------------------
chL = args.refASIC*32 + map_bar_LR[args.refBar][0]
chR = args.refASIC*32 + map_bar_LR[args.refBar][1]
print("\n-- Ref ch: \tchannel left ", chL, "   channel right ", chR)


# intercalib labels
# -------------------------
# the intercalib label is the one of the energy intercalibration csv files
intercalib_label = f"{base_label}_refBar{str(args.calibBar)}"

# when applying the energy intercalibrations, use the TOFHIR or TOFHIR_LO label to keep the info in the name of the file
if args.whichEnergyIntercalib != "":
    if args.extraLabel == "":
        args.extraLabel = f"{args.whichEnergyIntercalib}calib"
    else:
        args.extraLabel += f"_{args.whichEnergyIntercalib}calib"
    intercalib_path = "Plot_repo_path/energy_intercalibration/intercalibLabel/whichCalibration_calibration_factors.csv"
    intercalib_path_ref = f"Plot_repo_path/energy_intercalibration/{args.refCalibPath}/whichCalibration_calibration_factors.csv"
else:
    intercalib_path = "0"
    intercalib_path_ref = "0"

# general label
# ------------------------- 
if args.extraLabel:
    generalLabel = f"{base_label}_{args.extraLabel}"
else:
    generalLabel = base_label
print(f"-- Config: \t{args.config}")

    
# path
# ----------------------------
cfg_path = f"{Lab5015_path}/cfg"
cfgFolder = Path(cfg_path)
temp_min = cfgFolder / f"minEnergies_{args.modulelabel}.txt"


# copy minEnergy file if it does not exist
# ----------------------------
if not temp_min.exists():
    base_min = cfgFolder / "minEnergies_base.txt"
    shutil.copy(base_min, temp_min)

    
# function to write config file
# -------------------------------
def write_cfg(base_path: Path, out_path: Path, replacements: dict, check_Vov: bool = False):
    with open(base_path) as baseCfg, open(out_path, 'w') as newCfg:
        for line in baseCfg:
            original_line = line
            if check_Vov and line.startswith('Vov') and args.Vov not in line:
                print(f"ERROR: missing ov in {base_path.name} file")
                sys.exit(1)

            for key, val in replacements.items():
                if key in line:
                    line = line.replace(key, str(val))
            newCfg.write(line)


# ModuleCharacterization cfg
# ----------------------------
base_module_cfg = cfgFolder / "moduleCharacterization_base.cfg"
out_module_cfg = cfgFolder / f"moduleCharacterization_{generalLabel}.cfg"
print(f"-- Writing \t{out_module_cfg.name}")

# replace the words find in the cfg_base (key) with the item of this dict
replacements_module = {
    "interCalibPath" : intercalib_path,
    "interCalibpath_ref" : intercalib_path_ref,
    "Lab5015_repo_path": Lab5015_path,
    "Plot_repo_path": plot_path,
    "runNumbers": args.runs,
    "generalLabel": generalLabel,
    "moduleLabel": args.modulelabel,
    "confNumber": args.config,
    "vovLabel": args.Vov,
    "scanVth" : args.threshold,
    "setVth1" : args.thresholdT1,
    "setVth2" : args.thresholdT2,
    "setVthe" : args.thresholdE,
    "channelLeft" : chL,
    "channelRight": chR,
    "refASIC_ch" : args.refASIC,
    "dutASIC_ch" : args.dutASIC,
    "saveReferenceModuleInfoFlag" : args.saveRefInfoFlag,
    "channelMapping_list_from_bar0_to_bar15" : get_TOFHIR_channel_mapping(),
    "whichCalibration" : args.whichEnergyIntercalib,
    "intercalibLabel" : intercalib_label,
    "useAmpWalkPerBarFlag" : args.refAmpWalkPerBar
}

write_cfg(base_module_cfg, out_module_cfg, replacements_module, check_Vov=True)
