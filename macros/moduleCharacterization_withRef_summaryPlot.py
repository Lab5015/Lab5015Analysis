#! /usr/bin/env python
from reference_utils import *
# ---- YOUR PATH -----
eos_path = "/eos/home-s/spalluot/MTD/TB_CERN_Sep25/Lab5015Analysis/"
outdir  = '/eos/home-s/spalluot/www/MTD/MTDTB_CERN_Sep25/ModuleCharacterization/'
# ---------------
# --- arguments ---
parser = argparse.ArgumentParser(description='Module characterization summary plots with ref info')
parser.add_argument("-i",  "--inputLabel",   required=True, type=str, help="input label")
parser.add_argument("-o",  "--outFolder",     required=False, type=str, help="out folder, default=inputLabel")
parser.add_argument("-dm",  "--dm_id",     required=False, type=str, help="DM ID for moduleChar comparison")
args = parser.parse_args()

# - thRef at which moduleChar results were extracted
thRef = 10
# - time resolution measured with moduleChar tDiff vs bar for a few DMs
dm = {
    "9001" : [33.7608, 31.6392, 30.4032, 27.884, 29.1516, 29.2742, 29.3106, 27.1104, 27.7983, 28.4243, 30.0891, 29.3959, 32.2497, 31.2093, 30.9873, 34.3772],
    "9000" : [33.2457, 31.5417, 30.7123, 30.687, 29.0857, 28.9459, 29.9012, 29.074, 29.9042, 28.3576, 31.5272, 30.3022, 30.4894, 30.6531, 31.781, 33.1243],
    "9005" : [30.4445, 31.9715, 30.8398, 30.693, 30.1605, 29.8004, 28.5805, 30.3734, 30.6305, 28.886, 30.1636, 31.2088, 31.475, 30.2317, 31.6581, 32.567],
    "FE_4587" : [31.3033, 31.5171, 29.8797, 28.6355, 26.7881, 27.5626, 27.4332, 27.683, 27.872, 28.4059, 28.9689, 29.5674, 30.1509, 29.7557, 30.6913, 33.0048]
}



# ----- MAIN ------
# - Define paths
input_label = args.inputLabel
print(" - Input label: ", input_label)
if args.outFolder == None:
    args.outFolder = args.inputLabel        
inputdir = f'{eos_path}/plots/'
summarydir = inputdir
outdir = f"{outdir}/{args.outFolder}/summaryPlots_withRef/"
os.makedirs(outdir, exist_ok=True)
outFileName = summarydir+'/summaryPlots_withRef_'+args.outFolder+'.root'
print(' - Saving root file ', outFileName)
print(' - Saving plots in ', outdir)
outfile = ROOT.TFile(outFileName, 'RECREATE' )
f = ROOT.TFile.Open(f"{inputdir}/moduleCharacterization_step2_{input_label}.root")

# - Get histos from input                    
pattern = re.compile( r"h1_deltaT_(.+)_bar(\d+)_Vov3\.00_th([0-9\.]+)" )
data = defaultdict(lambda: defaultdict(dict))

# - Get functions from the input file
tf1_global = {}
for key in f.GetListOfKeys():
    obj = key.ReadObj()
    if obj.InheritsFrom("TF1"):
        tf1_global[obj.GetName()] = obj

# - Get histograms from the input file
for key in f.GetListOfKeys():
    obj = key.ReadObj()
    name = obj.GetName()
    if not obj.InheritsFrom("TH1"):
        continue
    m = pattern.match(name)
    if not m:
        continue
    quantity = m.group(1)
    bar = int(m.group(2))
    thr = int(m.group(3))
    h = obj
    # -- get the histo functions
    f1 = get_tf1_from_hist(h)
    if f1 is None:
        f1 = find_tf1_by_name(tf1_global, name)
    if f1 is None:
        print(f"[WARNING] No TF1 found for {name}")
        continue

    # -- get mean and sigma of the gaussian fit
    mean = f1.GetParameter(1)
    sigma = f1.GetParameter(2)
    emean  = f1.GetParError(1)
    esigma = f1.GetParError(2)
    data[quantity][thr][bar] = { "mean": mean, "emean": emean, "sigma": sigma, "esigma": esigma }
    
# - make plots of mean and sigma parameters from the fit vs bar
for quantity, thr_dict in data.items():
    for thr, bars in thr_dict.items():
        plot_vs_bar( outdir, {"" : data[quantity][thr] }, title=f"thr {thr}", ykey="sigma", ekey="esigma", plotlabel=f"sigma_vs_bar_{quantity}_th{thr:02d}", plotdir=quantity)
        plot_vs_bar( outdir, {"" : data[quantity][thr] }, title=f"thr {thr}", ykey="mean", ekey="emean", plotlabel=f"mean_vs_bar_{quantity}_th{thr:02d}", ylim=None, plotdir=quantity, ylabel=r"$\mu$ [ps]")
        
# - compute triplet sigmas as a cross check
derived = defaultdict(lambda: defaultdict(dict))
q_tLRef = "tLCor_tAveRefCor_phaseCor"
q_tRRef = "tRCor_tAveRefCor_phaseCor"
q_tLtR  = "tL_tR_eRatioCor_phaseMeanCor"
common_thrs = sorted( set(data[q_tLRef].keys()) & set(data[q_tRRef].keys()) & set(data[q_tLtR].keys()) )
for thr in common_thrs:
    common_bars = sorted( set(data[q_tLRef][thr].keys()) & set(data[q_tRRef][thr].keys()) & set(data[q_tLtR][thr].keys()) )
    for bar in common_bars:
        sigma12 = data[q_tLRef][thr][bar]["sigma"]
        sigma13 = data[q_tRRef][thr][bar]["sigma"]
        sigma23 = data[q_tLtR][thr][bar]["sigma"]
        esigma12 = data[q_tLRef][thr][bar]["esigma"]
        esigma13 = data[q_tRRef][thr][bar]["esigma"]
        esigma23 = data[q_tLtR][thr][bar]["esigma"]
        vals = sigma_triplet( sigma12=sigma12, sigma23=sigma23, sigma13=sigma13, err12=esigma12, err23=esigma23, err13=esigma13)
        derived["sRef"][thr][bar] = {"sigma": vals["sRef"][0], "esigma": vals["sRef"][1]}
        derived["sL"][thr][bar]   = {"sigma": vals["sL"][0],   "esigma": vals["sL"][1]  }
        derived["sR"][thr][bar]   = {"sigma": vals["sR"][0],   "esigma": vals["sR"][1]  }
        derived["sDiff"][thr][bar]= {"sigma": sigma23/2.0,     "esigma": esigma23/2.0   }
for thr in common_thrs:
    plot_vs_bar( outdir,  {"REF avg" : derived["sRef"][thr], "DUT L": derived["sL"][thr], "DUT R": derived["sR"][thr]}, title=f"thr {thr}", ykey="sigma", ekey="esigma", plotlabel=f"triplet_vs_bar_th{thr:02d}", plotdir="triplet_check")
    
# - compute average DUT bar time resolution 
derived["sAve"] = defaultdict(lambda: defaultdict(dict))
q_dut = "tAve_tAveRefCor_eCor_phaseCor"
for thr in data[q_dut].keys():
    common_bars = set(data[q_dut][thr].keys()) & set(derived["sRef"][thr].keys())
    for bar in common_bars:
        val, err = compute_sigma_differences(data,  q_dut, derived, "sRef", thr, bar)
        derived["sAve"][thr][bar] = {"sigma": val, "esigma": err}

# - tDiff vs tAvg comparison
plotname = f"tDiff_tAvg_vs_bar_th{thRef:02d}"
if args.dm_id and args.dm_id in dm.keys():
    dm_key = args.dm_id
    derived["sDiff_moduleChar"] = defaultdict(lambda: defaultdict(dict))
    bars_sorted = sorted(common_bars)
    for i, bar in enumerate(bars_sorted):
        derived["sDiff_moduleChar"][thRef][bar] = { "sigma": dm[dm_key][i], "esigma": 0 }
    plot_vs_bar( outdir,  {"tDiff": derived["sDiff"][thRef], "tAvg": derived["sAve"][thRef], "tDiff moduleChar": derived["sDiff_moduleChar"][thRef]}, title="", ykey="sigma", ekey="esigma", plotlabel=plotname)
else:
    if dm_key not in dm.keys():
        print("[WARNING] Missing time resolution values from moduleCharacterization tDiff for the specified DM ID \n")
    plot_vs_bar( outdir, {"tDiff": derived["sDiff"][thRef], "tAvg": derived["sAve"][thRef]}, title="", ykey="sigma", ekey="esigma", plotlabel=plotname)
    
print("\n - Time resolution values: \n")
for bar in derived["sDiff"][thRef]:
    sdiff = derived["sDiff"][thRef][bar]
    smod  = derived["sDiff_moduleChar"][thRef][bar]
    save  = derived["sAve"][thRef][bar]
    print( f"Bar {bar:2d} | "
           f"time difference referenceChar : {sdiff['sigma']:.0f} ± {sdiff['esigma']:.0f} ps | "
           f"time difference moduleChar    : {smod['sigma']:.0f} ± {smod['esigma']:.0f} ps | "
           f"time average referenceChar    : {save['sigma']:.0f} ± {save['esigma']:.0f} ps" )
