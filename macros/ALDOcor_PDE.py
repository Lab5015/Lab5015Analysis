import json
import math

# production HPK SiPM : 25 um ES2
def PDE(ov):
    return 0.638 * (1.0 - math.exp(-0.651 * ov))

def Gain(ov):
    return 7.044E04 + 2.895E05*ov

def getCorrectionFactors(SM_ID, OV):
    json_file = "/eos/home-s/spalluot/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/ALDOcor.json"
    with open(json_file, "r") as f:
        data = json.load(f)
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

def getPDEfactor(SM_ID, Vov):
    json_file = "/eos/home-s/spalluot/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/ALDOcor.json"
    with open(json_file, "r") as f:
        data = json.load(f)
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

# print(pde_gain_ratio(3.17, 2.82))
# 15 perc variations between elft right sides
# ovs = [(3.55, 3.0), (3.0, 2.55), (1.38,1.20), (1.20, 1.02), (1.03, 0.9), (0.9, 0.76)]
# print(" \t OV1  \t / \t OV1 \t\t [PDE impact, Gain impact]")
# for ovup, ovlow in ovs:
#     print(f" \t {ovup}  \t / \t {ovlow}  \t\t {pde_gain_ratio(ovup,ovlow)}")
