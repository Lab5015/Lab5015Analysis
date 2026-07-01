import json
import math

# path to ALDO correction values
json_file = "/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2025/SMs_QAQC/ALDOcor.json"
with open(json_file, "r") as f:
    data = json.load(f)

# ------- 
# production HPK SiPM : 25 um ES2
def PDE(ov):
    return 0.638 * (1.0 - math.exp(-0.651 * ov))

def Gain(ov):
    return 7.044E04 + 2.895E05*ov
# -------

# --- Compute correction factors for left/right sensor sides ---
def getCorrectionFactors(SM_ID, OV):
    # Find module configuration entry
    row = None
    for r in data:
        if r["id"] == SM_ID and float(r["ovSet"]) == float(OV):
            row = r
            break
    if row is None:
        raise ValueError(f"[ERROR] No match for id={SM_ID}, ov={OV}")
    ov_left  = row["ov0corr"]
    ov_right = row["ov1corr"]
    target = PDE(OV) * Gain(OV)

    # real values
    left_val  = PDE(ov_left)  * Gain(ov_left)
    right_val = PDE(ov_right) * Gain(ov_right)

    # correction factors
    left_factor  = target / left_val
    right_factor = target / right_val

    return { "left_factor": left_factor, "right_factor": right_factor }

# --- PDE normalization factors with corrected OV per side ---
def getPDEfactor(SM_ID, Vov):
    # find the row containing the SM_ID and the Vov
    row = None
    for r in data:
        if r["id"] == SM_ID and float(r["ovSet"]) == float(Vov):
            row = r
            break
    if row is None:
        raise ValueError(f"[ERROR] No match for id={SM_ID}, ov={Vov}")
    pde_val = PDE(Vov)
    return [ pde_val / row["ov0corr"], pde_val / row["ov1corr"] ]

# --- Relative change between two overvoltage points ---
def pde_gain_ratio(vov1, vov2):
    pde1 = PDE(vov1)
    pde2 = PDE(vov2)
    pde_ratio = pde1 / pde2
    pde = int((1-pde_ratio)*100)
    
    gain1 = Gain(vov1)
    gain2 = Gain(vov2)
    gain_ratio = gain1/ gain2
    gain = int((1-gain_ratio)*100)
    
    return [pde, gain]
