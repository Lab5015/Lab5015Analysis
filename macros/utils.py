#! /usr/bin/env python
import numpy as np
import argparse
import re
from collections import defaultdict
from scipy.interpolate import interp1d
from ctypes import c_double, c_float
import CMS_lumi, tdrstyle

from slewRate import *
from SiPM import *
from moduleDict import *
from calibration_utils import *

cms_colors = [
    ROOT.TColor.GetColor("#3f90da"),
    ROOT.TColor.GetColor("#ffa90e"),
    ROOT.TColor.GetColor("#bd1f01"),
    ROOT.TColor.GetColor("#94a4a2"),
    ROOT.TColor.GetColor("#832db6"),
    ROOT.TColor.GetColor("#a96b59"),
    ROOT.TColor.GetColor("#e76300"),
    ROOT.TColor.GetColor("#b9ac70"),
    ROOT.TColor.GetColor("#717581"),
    ROOT.TColor.GetColor("#92dadd")
]

# ----- Set TDR style ------
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.055,'X')
ROOT.gStyle.SetLabelSize(0.055,'Y')
ROOT.gStyle.SetTitleSize(0.06,'X')
ROOT.gStyle.SetTitleSize(0.06,'Y')
ROOT.gStyle.SetTitleOffset(1.05,'X')
ROOT.gStyle.SetTitleOffset(1.12,'Y')
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.045)
ROOT.gStyle.SetPadTopMargin(0.07)
ROOT.gStyle.SetPadRightMargin(0.1)
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(111)

# ------- Utility functions ------
def intersection(lst1,lst2):
    lstIntersection = [value for value in lst1 if value in lst2] 
    return lstIntersection

def union(lst1,lst2):
    return lst1 + list(set(lst2) - set(lst1))

def remove_useless_dirs(path):
    """Remove empty directories or those containing only 'index.php'."""
    for root, dirs, files in os.walk(path, topdown=False):
        for d in dirs:
            full = os.path.join(root, d)            
            file_list = [f for f in os.listdir(full) if os.path.isfile(os.path.join(full, f))]
            subdirs = [f for f in os.listdir(full) if os.path.isdir(os.path.join(full, f))]
            if len(subdirs) == 0 and (len(file_list) == 0 or file_list == ["index.php"]):
                print(f"Removing directory: {full}")
                if "index.php" in file_list:
                    os.remove(os.path.join(full, "index.php"))
                os.rmdir(full)

# -------- Latex labels -------
def draw_logo():
    logo_x = 0.16
    logo_right = 0.78
    logo = ROOT.TLatex()
    logo.SetNDC()
    logo.SetTextSize(0.06)
    logo.SetTextFont(62)
    logo.DrawText(logo_x,0.95,'CMS')
    logo.SetTextSize(0.045)
    logo.SetTextFont(52)
    logo.DrawText(logo_x+0.11, 0.95, 'Preliminary')
    logo.DrawText(logo_right, 0.95, 'Phase II')
    return logo

def latex_vov(overv):
    latex_tmp = ROOT.TLatex(0.19,0.88,'Vov%.2f'%overv)
    latex_tmp.SetNDC()
    latex_tmp.SetTextSize(0.035)
    latex_tmp.SetTextFont(42)
    return latex_tmp

def latex_sipm(sip_):
    latex_s = ROOT.TLatex(0.17,0.83,'%s'%label_(sip_))
    latex_s.SetNDC()
    latex_s.SetTextSize(0.035)
    latex_s.SetTextFont(42)
    return latex_s


def latex_bar(bar_):
    latex_b = ROOT.TLatex(0.19,0.65,'bar%02d'%bar_)
    latex_b.SetNDC()
    latex_b.SetTextSize(0.035)
    latex_b.SetTextFont(42)
    return latex_b

# -------- ROOT functions ----------
def remove_points_beyond_rms(g, rms_thr):
    """Return a TGraphErrors with points outside rms threshold removed"""
    if not g:
        return
    n      = g.GetN()
    x_vals = g.GetX()
    y_vals = g.GetY()
    x_errs = g.GetEX()
    y_errs = g.GetEY()
    mean   = g.GetMean(2)
    rms    = g.GetRMS(2)
    new_g  = ROOT.TGraphErrors()
    for i in range(g.GetN()):
        if abs(y_vals[i] - mean) <= rms_thr * rms:
            new_g.SetPoint(new_g.GetN(), x_vals[i], y_vals[i])
            new_g.SetPointError(new_g.GetN()-1, x_errs[i], y_errs[i])
        else:
            print(" ~~~~~~~~~~~   point at ", x_vals[i], ' removed')
    return new_g

def interpolate_error(graph, x_val,verbose=False):
    """Interpolate the y-error of a TGraphErrors at a given x value."""
    n_points = graph.GetN()
    x_values = np.zeros(n_points)
    y_values = np.zeros(n_points)
    y_errors = np.zeros(n_points)
    for i in range(n_points):
        x_values[i] = graph.GetX()[i]
        y_values[i] = graph.GetY()[i]
        y_errors[i] = graph.GetEY()[i]
        error_interpolator = interp1d(x_values, y_errors, bounds_error=False, fill_value="extrapolate")
    if verbose:    
        print("\ny values: ", y_values, "\ty err : ", y_errors)
        print("x values: ", x_values)
    return error_interpolator(x_val)


# -------- Langaus function --------------
"""
 --- ROOT tutorial: https://root.cern/doc/v636/langaus_8C.html
 Fit parameters:
   - par[0] = Width (scale) parameter of Landau density
   - par[1] = Most Probable (MP, location) parameter of Landau density
   - par[2] = Total area (integral -inf to inf, normalization constant)
   - par[3] = Width (sigma) of convoluted Gaussian function
 Note:
  In the Landau distribution (represented by the CERNLIB approximation),
  the maximum is located at x=-0.22278298 with the location parameter=0.
  This shift is corrected within this function, so that the actual
  maximum is identical to the MP parameter.
 Usage:
  f_lg = ROOT.TF1("f_langau", getattr(ROOT,"langaufun"), 0, 1000, 4)
"""
ROOT.gInterpreter.Declare("""
double langaufun(double *x, double *par) {
   double invsq2pi = 0.3989422804014;
   double mpshift  = -0.22278298;
   double np = 50.0;
   double sc = 5.0;
   double xx;
   double mpc;
   double fland;
   double sum = 0.0;
   double xlow,xupp;
   double step;
   mpc = par[1] - mpshift * par[0];
   xlow = x[0] - sc * par[3];
   xupp = x[0] + sc * par[3];
   step = (xupp-xlow) / np;
   for(int i=1; i<=np/2; i++) {
      xx = xlow + (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);
      xx = xupp - (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);
   }
   return par[2] * step * sum * invsq2pi / par[3];
}
""")


def correct_histogram_x(h_orig, p0=0.0, p1=1.0, new_name=None):
    """
    Return a histogram with linearly corrected x-axis (smoothed)
    Developed for correcting energy histogram with calibration factors per bin
    """
    if new_name is None:
        new_name = h_orig.GetName() + "_calib"
    h_corr = h_orig.Clone(new_name)
    h_corr.Reset()
    for ibin in range(1, h_orig.GetNbinsX() + 1):
        E = h_orig.GetBinCenter(ibin)
        N = h_orig.GetBinContent(ibin)
        if N == 0:
            continue
        E_corr = p0 + p1*E
        # find adjacent bins
        bin_left = max(1, h_corr.GetXaxis().FindBin(E_corr) - 1)
        bin_right = bin_left + 1
        x_left = h_corr.GetBinCenter(bin_left)
        x_right = h_corr.GetBinCenter(bin_right)
        # linear weights to adjacent bins
        if x_right == x_left:
            w_left = 1.0
            w_right = 0.0
        else:
            w_right = (E_corr - x_left) / (x_right - x_left)
            w_left = 1.0 - w_right
        # proportional content
        h_corr.SetBinContent(bin_left, h_corr.GetBinContent(bin_left) + N*w_left)
        h_corr.SetBinContent(bin_right, h_corr.GetBinContent(bin_right) + N*w_right)
    return h_corr


def fit_landau_langaus(h, emin, emax, landau_only=False):
    """
    Given an histogram and the x range, it provides two fit functions and a dict:
      - landau fit function
      - langaus fit function
      - a dictionary containing the fit parameters
    Usage: 
      f_landau, f_lg, result = fit_landau_langaus(h, emin, emax)
    """
    integral = h.Integral()
    rebin = get_rebin_factor(integral)
    if rebin > 1:
        h.Rebin(rebin)
    bin_min = h.GetXaxis().FindBin(emin)
    bin_max = h.GetXaxis().FindBin(emax)
    max_bin = max(range(bin_min, bin_max+1), key=lambda b: h.GetBinContent(b))
    max_content = h.GetBinContent(max_bin)
    max_energy  = h.GetBinCenter(max_bin)
    f_landau = ROOT.TF1("f_landau", "[0]*TMath::Landau(x,[1],[2])", emin, emax)
    f_landau.SetParLimits(1, 0, 900)
    f_landau.SetParLimits(2, 0, 100)
    f_landau.SetParameters(max_content, max_energy, max_energy*0.15)
    h.Fit(f_landau, "QRS")
    mpv_landau = f_landau.GetParameter(1)
    w_landau = f_landau.GetParameter(2)
    mpv_landau_err = f_landau.GetParError(1)
    w_landau_err = f_landau.GetParError(2)    
    if landau_only:
        return f_landau, {
            "landau_mpv" : mpv_landau,
            "landau_width" : w_landau,
            "landau_mpv_err" : mpv_landau_err,
            "landau_width_err" : w_landau_err
        }
    else:
        f_lg = ROOT.TF1("f_lg", getattr(ROOT,"langaufun"), mpv_landau*0.8, mpv_landau*1.6, 4)
        f_lg.SetParLimits(0, 0.5*w_landau, 1.5*w_landau)
        f_lg.SetParLimits(1, mpv_landau*0.5, mpv_landau*1.5)    # MPV
        f_lg.SetParLimits(2, 0.5*integral, 25*integral)
        f_lg.SetParLimits(3, 0.01*mpv_landau, 0.15*mpv_landau)   # Gaussian sigma
        f_lg.SetParameters(w_landau, mpv_landau, integral, 0.04*mpv_landau)
        h.Fit(f_lg, "QRS")
        mpv_lang = f_lg.GetParameter(1)
        sigma    = f_lg.GetParameter(3)
        mpv_lang_err = f_lg.GetParError(1)
        sigma_err    = f_lg.GetParError(3)
        return f_landau, f_lg, {
            "landau_mpv" : mpv_landau,
            "landau_width" : w_landau,
            "landau_mpv_err" : mpv_landau_err,
            "landau_width_err" : w_landau_err,
            "langaus_width" : f_lg.GetParameter(0),
            "langaus_mpv": mpv_lang,
            "langaus_sigma": sigma,
            "langaus_mpv_err": mpv_lang_err,
            "langaus_sigma_err": sigma_err
        }
