from utils import *
# - Each histogram has already been fitted with a gaussian in the step2
#    - get functions from the histo function
def get_tf1_from_hist(h):
    """Try to get TF1 attached to histogram"""
    funcs = h.GetListOfFunctions()
    if funcs:
        for obj in funcs:
            if obj.InheritsFrom("TF1"):
                return obj
    return None
def find_tf1_by_name(tf1_map, hist_name):
    """Fallback: match TF1 in global list"""
    for name, f in tf1_map.items():
        if hist_name in name:
            return f
    return None

def safe_sqrt(x):
    return math.sqrt(x) if x > 0 else 999

# - il triangolo no legacy
def err(s, A, sA, B, sB, C, sC):
    invS = 0.5 * 1/s
    A_piece = (invS * A * sA)**2
    B_piece = (invS * B * sB)**2
    C_piece = (invS * C * sC)**2
    return safe_sqrt(A_piece + B_piece + C_piece)

def sigma_triplet(sigma12, sigma23, sigma13, err12=0., err23=0., err13=0., sigmaCorr=0):
    argRef = 0.5 * (sigma12**2 + sigma13**2 - sigma23**2 - sigmaCorr**2)
    argL   = 0.5 * (sigma12**2 + sigma23**2 - sigma13**2 - sigmaCorr**2)
    argR   = 0.5 * (sigma13**2 + sigma23**2 - sigma12**2) - sigmaCorr**2

    sRef = safe_sqrt(argRef)
    sL   = safe_sqrt(argL)
    sR   = safe_sqrt(argR)

    err_sRef = err(sRef, sigma12, err12, sigma23, err23, sigma13, err13)
    err_sL   = err(sL,   sigma12, err12, sigma23, err23, sigma13, err13)
    err_sR   = err(sR,   sigma12, err12, sigma23, err23, sigma13, err13)    
    return { "sRef": (sRef, err_sRef), "sL": (sL, err_sL), "sR": (sR, err_sR)}

def compute_sigma_differences(data1, qt1, data2, qt2, thr, bar):
    sigma_dut_ref = data1[qt1][thr][bar]["sigma"]
    sigma_ref     = data2[qt2][thr][bar]["sigma"]
    err_dut_ref   = data1[qt1][thr][bar]["esigma"]
    err_ref       = data2[qt2][thr][bar]["esigma"]

    val = safe_sqrt(sigma_dut_ref**2 - sigma_ref**2)
    piece_dut_ref = (1/val * sigma_dut_ref * err_dut_ref)**2
    piece_ref     = (1/val * sigma_ref * err_ref)**2 
    err = safe_sqrt(piece_dut_ref + piece_ref)
    return (val, err)

def plot_vs_bar(outdir, series_dict, thr, title, plotlabel, ykey, ekey, ylabel=r"$\sigma$ [ps]", ylim=(0,100), plotdir=None):
    if plotdir:
        plot_dir = os.path.join(outdir, plotdir)
        os.makedirs(plot_dir, exist_ok=True)
    else:
        plot_dir = outdir
    fig, ax = plt.subplots(figsize=(12,10))
    for label, values_dict in series_dict.items():
        bars = sorted(values_dict.keys())
        y = np.array([values_dict[b][ykey] for b in bars])
        yerr = np.array([values_dict[b][ekey] for b in bars])
        ax.errorbar( bars, y, yerr=yerr, marker='o', linestyle='-', capsize=3, label=label)
    ax.set_xlabel("Bar")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid()
    if ylim:
        ax.set_ylim(*ylim)
    ax.legend()
    # Stats 
    text_lines = []
    for label, values_dict in series_dict.items():
        yvals = np.array([values_dict[b][ykey] for b in sorted(values_dict.keys())])
        mean = np.nanmean(yvals)
        rms  = np.nanstd(yvals)
        text_lines.append(f"{label}: mean={mean:.0f}, RMS={rms:.0f}")
    text = "\n\n".join(text_lines)        
    ax.text( 0.55, 0.65, text, transform=ax.transAxes, fontsize=20, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none') )
    plt.savefig(os.path.join(plot_dir, f"{plotlabel}_th{thr:02d}.png"))
    plt.close()
