import ROOT
import argparse
from utils import *

plotdir = "/eos/home-s/spalluot/www/MTD/MTDTB_CERN_Sep25/check_occupancy/"

# argparse
parser = argparse.ArgumentParser(description="2D Ref vs DUT with energy cuts")
parser.add_argument("--run", type=str, default="4140")
parser.add_argument("--emin_ref", type=float, default=200)
parser.add_argument("--emax_ref", type=float, default=1200)
parser.add_argument("--thr_down", type=float, default=100)
args = parser.parse_args()

run_number = args.run
emin_ref = args.emin_ref
emax_ref = args.emax_ref
thr_down = args.thr_down

# dataframe including all spills of the specified run
df = ROOT.RDataFrame("data",f"/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Sep2025/TOFHIR/reco/{run_number}/*.root")

# define channel mapping
ROOT.gInterpreter.Declare("""
// Upstream
std::map<unsigned int,int> chL = {
  {132,0},{129,1},{128,2},{131,3},{130,4},{134,5},{135,6},{137,7},
  {133,8},{139,9},{136,10},{140,11},{138,12},{142,13},{143,14},{141,15}
};
std::map<unsigned int,int> chR = {
  {157,0},{153,1},{152,2},{158,3},{159,4},{156,5},{155,6},{150,7},
  {154,8},{149,9},{151,10},{148,11},{147,12},{146,13},{144,14},{145,15}
};
// Downstream
std::map<unsigned int,int> chDL = {
  {228,0},{225,1},{224,2},{227,3},{226,4},{230,5},{231,6},{233,7},
  {229,8},{235,9},{232,10},{236,11},{234,12},{238,13},{239,14},{237,15}
};
std::map<unsigned int,int> chDR = {
  {253,0},{249,1},{248,2},{254,3},{255,4},{252,5},{251,6},{246,7},
  {250,8},{245,9},{247,10},{244,11},{243,12},{242,13},{240,14},{241,15}
};

""")

# define ref energy L+R
df2 = df.Define("refEnergy",
    """
    std::map<int,float> barE;
    for (size_t i=0;i<channelID.size();++i){
        auto ch = channelID[i];
        if (chL.count(ch)) barE[chL[ch]] += energy[i];
        if (chR.count(ch)) barE[chR[ch]] += energy[i];
    }
    return barE;
    """
)

# ref bar scan
df3 = df2.Define("refBar",
    f"""
    int selected = -1;
    int n = 0;
    for (auto const &p : refEnergy){{
        if (p.second > {emin_ref} && p.second < {emax_ref}){{
            selected = p.first;
            n++;
        }}
    }}
    if (n == 1) return selected;
    else return -1;
    """
).Filter("refBar >= 0")

# DUT counts with a conservative cut on emin
df4 = df3.Define("barD",
    f"""
    std::set<int> bars;
    for (size_t i=0;i<channelID.size();++i){{
        auto ch = channelID[i];
        if (chDL.count(ch) && energy[i] > {thr_down}) bars.insert(chDL[ch]);
        if (chDR.count(ch) && energy[i] > {thr_down}) bars.insert(chDR[ch]);
    }}
    return std::vector<int>(bars.begin(), bars.end());
    """
)

df5 = df4.Define("refVec", "return std::vector<int>(barD.size(), refBar);")
h2 = df5.Histo2D(("h2", ";Ref bar;DUT bar", 16,0,16,16,0,16),"refVec","barD")
c = ROOT.TCanvas("c", "", 700,600)
ROOT.gStyle.SetOptStat(0)
h2.Draw("COLZ")
c.SaveAs(f"{plotdir}/beam_profile_modules_occupancy_run{run_number}.png")
