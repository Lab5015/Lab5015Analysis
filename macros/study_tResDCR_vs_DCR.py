#! /usr/bin/env python
import os 
import shutil
import glob
import math
import array
import sys
import time
import argparse                                                                                                                                                          
import json
import numpy as np

from SiPM import *

import ROOT                                                                                                                                                              
import CMS_lumi, tdrstyle                                                                                                                                                                          
#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetFitFormat('3.2g')
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.04,'X')
ROOT.gStyle.SetLabelSize(0.04,'Y')
ROOT.gStyle.SetTitleSize(0.04,'X')
ROOT.gStyle.SetTitleSize(0.04,'Y')
ROOT.gStyle.SetTitleOffset(1.1,'X')
ROOT.gStyle.SetTitleOffset(1.4,'Y')
ROOT.gStyle.SetMarkerSize(1)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

from VovsEff import *



def computeVovEff(vov_set, Iarray, Ich):
    vov_eff =  vov_set - 1E-03*(Iarray*10. + Ich*240.)  
    return(vov_eff)


def getNpe(Ipho, laserFreq, gain, xtalk):
    npe = Ipho /(1.6E-19 * laserFreq * gain) / xtalk 
    return npe

# --------------------


setVovs = [0.80, 1.00, 1.20, 1.50, 2.00, 3.50]
laserTune = '77.5'
leds = ['8', '8.5', '9', '9.5', '10', '10.5', '11', '11.5','12'] 
sipm_type = 'HPK-PIT-C25-ES2'
laserFreq = 50000

fIpho  = 0.75 # 
fIpho_err = 0.10 # errore relativo
#fIdcrL = 1./16
#fIdcrR = 1./16
fIdcrL = 0.049
fIdcrR = 0.055
fIdcr_err = 0.15 # errore relativo

#outdir 
outdir = '/var/www/html/TOFHIR2C/904/study_DCR_LaserTune%.1f_v2/'%float(laserTune)

# input files  ( Vov, LED) : run )
runs_dict_DCR = {}
firstRun = 181
if (laserTune == 77.5): firstRun = 235
print('First run: %d'%firstRun)
it = 0
for led in leds:
    for ov in [1.0, 0.8, 1.20, 1.50, 2.00, 3.50]:
        runs_dict_DCR[ (ov, led)] =  firstRun +  it
        it = it + 1

print runs_dict_DCR


runs_dict_noDCR = {}
firstRun = 175
if (laserTune == 77.5): firstRun = 163    
print('First run: %d'%firstRun)
it = 0
for ov in [1.0, 0.8, 1.20, 1.50, 2.00, 3.50]:
    runs_dict_noDCR[ (ov, '0')] =  firstRun +  it
    it = it + 1

print runs_dict_noDCR
 

# IV with and without LED
g_IV_dark_L = {}
g_IV_dark_R = {}

g_IV_laser_L = {}
g_IV_laser_R = {}

g_IV_laser_led_L = {}
g_IV_laser_led_R = {}

Vbd_L=38.08
Vbd_R=38.15

f = {}

# dark: laser and LED off
fNameA = glob.glob('/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_gen23/logs_41.00/logIV_LYSO828_dark_ASIC0_ALDOA_ch2_time_2023-06-15_10:47:37.root')
fNameB = glob.glob('/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_gen23/logs_41.00/logIV_LYSO828_dark_ASIC0_ALDOB_ch1_time_2023-06-15_10:54:13.root')
print(fNameA)
print(fNameB)
f_L = ROOT.TFile.Open(fNameA[0])
f_R = ROOT.TFile.Open(fNameB[0])
g_IV_dark_L = f_L.Get('g_IV')
g_IV_dark_R = f_R.Get('g_IV')

# laser ON
fNameA = glob.glob('/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_gen23/logs_41.00/logIV_LYSO828_laser_50000_%s_ASIC0_ALDOA_ch2_*root'%laserTune)
fNameB = glob.glob('/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_gen23/logs_41.00/logIV_LYSO828_laser_50000_%s_ASIC0_ALDOB_ch1_*root'%laserTune)
print(fNameA)
print(fNameB)
f_L = ROOT.TFile.Open(fNameA[0])
f_R = ROOT.TFile.Open(fNameB[0])
g_IV_laser_L = f_L.Get('g_IV')
g_IV_laser_R = f_R.Get('g_IV')

# laser + LED ON
for led in leds:
    fNameA = glob.glob('/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_gen23/logs_41.00/logIV_LYSO828_laser_50000_%s_led_%s_ASIC0_ALDOA_ch2_*.root'%(laserTune,led))
    fNameB = glob.glob('/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_gen23/logs_41.00/logIV_LYSO828_laser_50000_%s_led_%s_ASIC0_ALDOB_ch1_*.root'%(laserTune,led))
    print(fNameA)
    print(fNameB)
    f_L = ROOT.TFile.Open(fNameA[0])
    f_R = ROOT.TFile.Open(fNameB[0])
    g_IV_laser_led_L[led] = f_L.Get('g_IV')
    g_IV_laser_led_R[led] = f_R.Get('g_IV')


# get energy, time resolution  vs Vov - NO DCR
g_tRes_vs_Vov_noDCR = ROOT.TGraphErrors()     
g_tRes_vs_Npe_noDCR = ROOT.TGraphErrors()     

for vov_set in setVovs:        
    fName = '../plots/summaryPlots_run%d.root'% runs_dict_noDCR[vov_set, '0']
    ff = ROOT.TFile.Open(fName)
    gRes = ff.Get('g_deltaT_totRatioCorr_bestTh_vs_vov_enBin01_average')        
    
    IphoL = g_IV_laser_L.Eval(Vbd_L+vov_set) - g_IV_dark_L.Eval(Vbd_L+vov_set)
    IphoR = g_IV_laser_R.Eval(Vbd_R+vov_set) - g_IV_dark_R.Eval(Vbd_R+vov_set)
    gain = Gain(sipm_type, vov_set)
    npe = 0.5 * ( getNpe(IphoL*1E-06, laserFreq, gain, 1./fIpho)  + getNpe(IphoR*1E-06, laserFreq, gain, 1./fIpho) )
    g_tRes_vs_Npe_noDCR.SetPoint(g_tRes_vs_Npe_noDCR.GetN(), npe, gRes.Eval(vov_set))
    g_tRes_vs_Vov_noDCR.SetPoint(g_tRes_vs_Vov_noDCR.GetN(), vov_set, gRes.Eval(vov_set))


# get energy, time resolution with DCR + estimate VovEff from energy
g_Ich_L_vs_VovEff = {}
g_Ich_R_vs_VovEff = {}
g_Ich_vs_VovEff = {}

g_DCR_L_vs_VovEff = {}
g_DCR_R_vs_VovEff = {}
g_DCR_vs_VovEff = {}

g_DCR_L_vs_PDE = {}
g_DCR_R_vs_PDE = {}
g_DCR_vs_PDE = {}

g_tRes_vs_VovEff_withDCR = {} 


vovEffL = {}
vovEffR = {}
vovEff  = {}
vovEffL_err  = {}
vovEffR_err  = {}
vovEff_err  = {}
npe = {}

for led in leds:
    g_Ich_L_vs_VovEff[led] = ROOT.TGraphErrors()
    g_Ich_R_vs_VovEff[led] = ROOT.TGraphErrors()
    g_Ich_vs_VovEff[led] = ROOT.TGraphErrors()

    g_DCR_L_vs_VovEff[led] = ROOT.TGraphErrors()
    g_DCR_R_vs_VovEff[led] = ROOT.TGraphErrors()
    g_DCR_vs_VovEff[led] = ROOT.TGraphErrors()

    g_DCR_L_vs_PDE[led] = ROOT.TGraphErrors()
    g_DCR_R_vs_PDE[led] = ROOT.TGraphErrors()                                                                                                               
    g_DCR_vs_PDE[led] = ROOT.TGraphErrors()

    g_tRes_vs_VovEff_withDCR[led] = ROOT.TGraphErrors()    

    for vov_set in setVovs:

        IarrayL = g_IV_laser_led_L[led].Eval(Vbd_L+vov_set)/1000 # uA->mA
        IarrayR = g_IV_laser_led_R[led].Eval(Vbd_R+vov_set)/1000 #uA->mA
        Iarray  = 0.5*(IarrayL+IarrayR)

        # approx. (it should be evaluated in OV_eff)
        IphoL   = (g_IV_laser_L.Eval(Vbd_L+vov_set) - g_IV_dark_L.Eval(Vbd_L+vov_set))/1000. # uA -> mA
        IphoR   = (g_IV_laser_R.Eval(Vbd_R+vov_set) - g_IV_dark_R.Eval(Vbd_R+vov_set))/1000. # uA -> mA

        IchL = IphoL*fIpho + (IarrayL - IphoL) * fIdcrL
        IchR = IphoR*fIpho + (IarrayR - IphoR) * fIdcrR
        Ich  = 0.5*(IchL+IchR)
                        
        vovEffL[led, vov_set] = computeVovEff(vov_set, IarrayL, IchL)
        vovEffR[led, vov_set] = computeVovEff(vov_set, IarrayR, IchR)
        vovEff[led, vov_set] = 0.5*(vovEffL[led, vov_set]+vovEffR[led, vov_set])    
        vovEffL_err[led, vov_set] = 0.5*abs(computeVovEff(vov_set, IarrayL, IchL*(1-fIdcr_err)) - computeVovEff(vov_set, IarrayL, IchL*(1+fIdcr_err))) 
        vovEffR_err[led, vov_set] = 0.5*abs(computeVovEff(vov_set, IarrayR, IchR*(1-fIdcr_err)) - computeVovEff(vov_set, IarrayR, IchR*(1+fIdcr_err))) 
        vovEff_err[led, vov_set] = 0.5 * (vovEffL_err[led, vov_set]+vovEffR_err[led, vov_set])

        gain = Gain(sipm_type, vovEff[led, vov_set])
        npe[led, vov_set] = 0.5 * ( getNpe(IphoL*1E-03, laserFreq, gain, 1./fIpho)  + getNpe(IphoR*1E-03, laserFreq, gain, 1./fIpho) ) 
        
        g_Ich_L_vs_VovEff[led].SetPoint(g_Ich_L_vs_VovEff[led].GetN(), vovEffL[led, vov_set], IchL)
        g_Ich_R_vs_VovEff[led].SetPoint(g_Ich_R_vs_VovEff[led].GetN(), vovEffR[led, vov_set], IchR)
        g_Ich_vs_VovEff[led].SetPoint(g_Ich_vs_VovEff[led].GetN(), vovEff[led, vov_set], Ich)

        dcrL = IchL*1E-03 / (1.602E-19 * Gain(sipm_type, vovEffL[led, vov_set])) * 1E-09 
        dcrR = IchR*1E-03 / (1.602E-19 * Gain(sipm_type, vovEffR[led, vov_set])) * 1E-09 
        dcr  = Ich*1E-03 / (1.602E-19 * Gain(sipm_type, vovEff[led, vov_set])) * 1E-09
        g_DCR_L_vs_VovEff[led].SetPoint(g_DCR_L_vs_VovEff[led].GetN(), vovEffL[led, vov_set], dcrL)
        g_DCR_R_vs_VovEff[led].SetPoint(g_DCR_R_vs_VovEff[led].GetN(), vovEffR[led, vov_set], dcrR)
        g_DCR_vs_VovEff[led].SetPoint(g_DCR_vs_VovEff[led].GetN(), vovEff[led, vov_set], dcr)
        
        g_DCR_L_vs_PDE[led].SetPoint(g_DCR_L_vs_PDE[led].GetN(), PDE(sipm_type, vovEffL[led, vov_set]), dcrL)
        g_DCR_R_vs_PDE[led].SetPoint(g_DCR_R_vs_PDE[led].GetN(), PDE(sipm_type, vovEffR[led, vov_set]), dcrR)
        g_DCR_vs_PDE[led].SetPoint(g_DCR_vs_PDE[led].GetN(), PDE(sipm_type, vovEff[led, vov_set]), dcr)


        fName = '../plots/summaryPlots_run%d.root'% runs_dict_DCR[vov_set,led]
        ff = ROOT.TFile.Open(fName)
        gRes = ff.Get('g_deltaT_totRatioCorr_bestTh_vs_vov_enBin01_average')        
        g_tRes_vs_VovEff_withDCR[led].SetPoint( g_tRes_vs_VovEff_withDCR[led].GetN(), vovEff[led,vov_set], gRes.Eval(vov_set)  )


leg2 = ROOT.TLegend(0.60, 0.70, 0.89, 0.89)
for i,led in enumerate(leds):
    c = ROOT.TCanvas('c_timeResolution_vs_VovEff_LED_%.1fV'%float(led), '', 700, 600)
    hPad = ROOT.gPad.DrawFrame(0.,0.,4.0,80.)
    hPad.SetTitle(";V_{OV}^{eff} [V]; #sigma_{t} [ps]")
    hPad.Draw()
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    g_tRes_vs_VovEff_withDCR[led].Draw('plsame')
    g_tRes_vs_Vov_noDCR.SetMarkerStyle(25)
    g_tRes_vs_Vov_noDCR.Draw('plsame')
    if (i==0):
        leg2.AddEntry(g_tRes_vs_VovEff_withDCR[led],'with DCR', 'PL')
        leg2.AddEntry(g_tRes_vs_Vov_noDCR,'without DCR', 'PL')
    leg2.Draw('same')
    c.SaveAs(outdir+c.GetName()+'.png')  
    c.SaveAs(outdir+c.GetName()+'.pdf')  

    c = ROOT.TCanvas('c_Ich_vs_VovEff_LED_%.1fV'%float(led), '', 700, 600)
    hPad = ROOT.gPad.DrawFrame(0.0,0.,4.0,2.)
    hPad.SetTitle(";V_{OV}^{eff} [V]; I_{channel} [mA]")
    hPad.Draw()
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    g_Ich_R_vs_VovEff[led].SetMarkerStyle(24)
    g_Ich_L_vs_VovEff[led].Draw('plsame')
    g_Ich_R_vs_VovEff[led].Draw('plsame')
    g_Ich_vs_VovEff[led].SetLineStyle(2)
    g_Ich_vs_VovEff[led].Draw('lsame')
    c.SaveAs(outdir+c.GetName()+'.png')
    c.SaveAs(outdir+c.GetName()+'.pdf')

    c = ROOT.TCanvas('c_DCR_vs_VovEff_LED_%.1fV'%float(led), '', 700, 600)
    hPad = ROOT.gPad.DrawFrame(0.0,0.,4.0,20.)
    hPad.SetTitle(";V_{OV}^{eff} [V]; DCR [GHz]")
    hPad.Draw()
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    g_DCR_R_vs_VovEff[led].SetMarkerStyle(24)
    g_DCR_L_vs_VovEff[led].Draw('plsame')
    g_DCR_R_vs_VovEff[led].Draw('plsame')
    g_DCR_vs_VovEff[led].SetLineStyle(2)
    g_DCR_vs_VovEff[led].Draw('lsame')
    c.SaveAs(outdir+c.GetName()+'.png')
    c.SaveAs(outdir+c.GetName()+'.pdf')

    c = ROOT.TCanvas('c_DCR_vs_PDE_LED_%.1fV'%float(led), '', 700, 600)
    hPad = ROOT.gPad.DrawFrame(0.0,0.,0.7,20.)
    hPad.SetTitle(";PDE; DCR [GHz]")
    hPad.Draw()
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    g_DCR_R_vs_PDE[led].SetMarkerStyle(24)
    g_DCR_L_vs_PDE[led].Draw('plsame')
    g_DCR_R_vs_PDE[led].Draw('plsame')
    g_DCR_vs_PDE[led].SetLineStyle(2)
    g_DCR_vs_PDE[led].Draw('lsame')
    c.SaveAs(outdir+c.GetName()+'.png')
    c.SaveAs(outdir+c.GetName()+'.pdf')



# now computre DCR contribution to the time resolution
g_tRes_DCR_all = ROOT.TGraphErrors()
for vov_set in setVovs:
    for led in leds:
        vov_eff = vovEff[led,vov_set]
        vov_eff_err = vovEff_err[led,vov_set]
        if (vov_eff< 0.80): continue # no DCR tRes in [0.80, 3.50]
        tRes0 = g_tRes_vs_Vov_noDCR.Eval(vov_eff)
        tRes  = g_tRes_vs_VovEff_withDCR[led].Eval(vov_eff)
        tRes_DCR = 0
        if (tRes > tRes0): 
            tRes_DCR = math.sqrt( tRes*tRes - tRes0*tRes0)
        else:
            continue
        dcr =  g_DCR_vs_VovEff[led].Eval(vov_eff)
        err_dcr = 0.5*abs( g_DCR_vs_VovEff[led].Eval(vov_eff+vov_eff_err) -g_DCR_vs_VovEff[led].Eval(vov_eff-vov_eff_err) )
        print led, vovEff[led, vov_set], npe[led, vov_set], tRes, tRes_DCR
        g_tRes_DCR_all.SetPoint(g_tRes_DCR_all.GetN(), dcr, npe[led,vov_set]/6000*tRes_DCR)
        g_tRes_DCR_all.SetPointError(g_tRes_DCR_all.GetN()-1, fIdcr_err*dcr, npe[led,vov_set]*fIpho_err/6000*tRes_DCR) 

            
c = ROOT.TCanvas('c_timeResolutionDCR_all', 'c_timeResolutionDCR_all', 700, 600)
hPad = ROOT.gPad.DrawFrame(0.,0.,40.,80.)
hPad.SetTitle(";DCR[GHz]; (Npe/6000) x #sigma_{t,DCR} [ps]")
hPad.Draw()
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
g_tRes_DCR_all.Draw('psame')
fitFun_all = ROOT.TF1('fitFun_all', '[1]*pow(x/30.,[0])' ,0, 100)
fitFun_all.SetParameters(0.4, 30.)
g_tRes_DCR_all.Fit(fitFun_all, 'QRS')
c.SaveAs(outdir+c.GetName()+'.png')
c.SaveAs(outdir+c.GetName()+'.pdf')

