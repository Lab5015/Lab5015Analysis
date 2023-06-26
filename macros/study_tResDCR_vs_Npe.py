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
laserTunes = ['74', '74.5','75','75.5','76','76.5','77.1','77.6','78.1']
sipm_type = 'HPK-PIT-C25-ES2'
laserFreq = 50000

fIpho  = 0.75
fIpho_err = 0.10 # errore relativo
#fIdcrL = 1./16
#fIdcrR = 1./16
fIdcrL = 0.055
fIdcrR = 0.055
fIdcr_err = 0.20 # errore relativo

thRef = 11

#outdir 
outdir = '/var/www/html/TOFHIR2C/904/study_DCR_vs_Npe_LED_12V/'

# input files  ( Vov, laserTune) : run )
runs_dict_DCR = {}
firstRun = 289
print('First run: %d'%firstRun)
it = 0
for laser in laserTunes:
    for ov in [1.0, 0.8, 1.20, 1.50, 2.00, 3.50]:
        runs_dict_DCR[ (ov, laser)] =  firstRun +  it
        it = it + 1

print runs_dict_DCR


runs_dict_noDCR = {}
firstRun = 439
print('First run: %d'%firstRun)
it = 0
for laser in laserTunes:
    for ov in [1.0, 0.8, 1.20, 1.50, 2.00, 3.50]:
        runs_dict_noDCR[ (ov, laser)] =  firstRun +  it
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
for laser in laserTunes:

    # dark: laser and LED off
    fNameA = glob.glob('/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_gen23/logs_41.00/logIV_LYSO828_dark_ASIC0_ALDOA_ch2_time_2023-06-15_10:47:37.root')[0]
    fNameB = glob.glob('/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_gen23/logs_41.00/logIV_LYSO828_dark_ASIC0_ALDOB_ch1_time_2023-06-15_10:54:13.root')[0]
    print(fNameA)
    print(fNameB)
    f_L = ROOT.TFile.Open(fNameA)
    f_R = ROOT.TFile.Open(fNameB)
    g_IV_dark_L[laser] = f_L.Get('g_IV')
    g_IV_dark_R[laser] = f_R.Get('g_IV')

    # laser ON
    fNameA = glob.glob('/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_gen23/logs_41.00/logIV_LYSO828_laser_50000_%s_led_0_ASIC0_ALDOA_ch2_*.root'%laser)[0]
    fNameB = glob.glob('/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_gen23/logs_41.00/logIV_LYSO828_laser_50000_%s_led_0_ASIC0_ALDOB_ch1_*.root'%laser)[0]
    print(fNameA)
    print(fNameB)
    f_L = ROOT.TFile.Open(fNameA)
    f_R = ROOT.TFile.Open(fNameB)
    g_IV_laser_L[laser] = f_L.Get('g_IV')
    g_IV_laser_R[laser] = f_R.Get('g_IV')

    # laser + LED ON
    fNameA = glob.glob('/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_gen23/logs_41.00/logIV_LYSO828_laser_50000_%s_led_12_ASIC0_ALDOA_ch2_*.root'%laser)[0]
    fNameB = glob.glob('/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_gen23/logs_41.00/logIV_LYSO828_laser_50000_%s_led_12_ASIC0_ALDOB_ch1_*.root'%laser)[0]
    print(fNameA)
    print(fNameB)
    f_L = ROOT.TFile.Open(fNameA)
    f_R = ROOT.TFile.Open(fNameB)
    g_IV_laser_led_L[laser] = f_L.Get('g_IV')
    g_IV_laser_led_R[laser] = f_R.Get('g_IV')


# get energy, time resolution  vs Vov - NO DCR
g_energyL_vs_Vov = {}
g_energyR_vs_Vov = {}
g_Vov_vs_energyL = {}
g_Vov_vs_energyR = {}


Npe = {}
g_Npe_vs_laserTune = {}
g_Npe_vs_energy = {}
g_tRes_vs_Npe_noDCR = {}
g_tRes_vs_Vov_noDCR = {}
for vov_set in setVovs: 
    g_Npe_vs_laserTune[vov_set]  = ROOT.TGraphErrors()
    g_tRes_vs_Npe_noDCR[vov_set] = ROOT.TGraphErrors()
    g_Npe_vs_energy[vov_set] = ROOT.TGraphErrors()


for laser in laserTunes:
    g_energyL_vs_Vov[laser] = ROOT.TGraphErrors()
    g_energyR_vs_Vov[laser] = ROOT.TGraphErrors()
    g_Vov_vs_energyL[laser] = ROOT.TGraphErrors()
    g_Vov_vs_energyR[laser] = ROOT.TGraphErrors()
    g_tRes_vs_Vov_noDCR[laser] = ROOT.TGraphErrors() 

    for vov_set in setVovs:        
        fName = '../plots/summaryPlots_run%d.root'% runs_dict_noDCR[vov_set,laser]
        ff = ROOT.TFile.Open(fName)
        gRes = ff.Get('g_deltaT_totRatioCorr_bestTh_vs_vov_enBin01_average')        
        
        gL = ff.Get('g_energyL_vs_th_bar00_Vov%.2f'%vov_set)
        gR = ff.Get('g_energyR_vs_th_bar00_Vov%.2f'%vov_set)

        IphoL = g_IV_laser_L[laser].Eval(Vbd_L+vov_set) - g_IV_dark_L[laser].Eval(Vbd_L+vov_set)
        IphoR = g_IV_laser_R[laser].Eval(Vbd_R+vov_set) - g_IV_dark_R[laser].Eval(Vbd_R+vov_set)
        gain = Gain(sipm_type, vov_set)
        npe = 0.5 * ( getNpe(IphoL*1E-06, laserFreq, gain, 1./fIpho)  + getNpe(IphoR*1E-06, laserFreq, gain, 1./fIpho) )
        g_Npe_vs_laserTune[vov_set].SetPoint(g_Npe_vs_laserTune[vov_set].GetN(), float(laser), npe)
        g_tRes_vs_Npe_noDCR[vov_set].SetPoint(g_tRes_vs_Npe_noDCR[vov_set].GetN(), npe, gRes.Eval(vov_set))
        g_tRes_vs_Vov_noDCR[laser].SetPoint(g_tRes_vs_Vov_noDCR[laser].GetN(), vov_set, gRes.Eval(vov_set))

        
        if (vov_set > 1.5): continue # different delayE used
        g_Npe_vs_energy[vov_set].SetPoint(g_Npe_vs_energy[vov_set].GetN(), gL.Eval(thRef) , npe)   

        g_energyL_vs_Vov[laser].SetPoint(g_energyL_vs_Vov[laser].GetN(), vov_set, gL.Eval(thRef))
        g_energyL_vs_Vov[laser].SetPointError(g_energyL_vs_Vov[laser].GetN()-1, 0, gL.GetErrorY(3))
    
        g_energyR_vs_Vov[laser].SetPoint(g_energyR_vs_Vov[laser].GetN(), vov_set, gR.Eval(thRef))
        g_energyR_vs_Vov[laser].SetPointError(g_energyR_vs_Vov[laser].GetN()-1, 0, gR.GetErrorY(3))

        g_Vov_vs_energyL[laser].SetPoint(g_Vov_vs_energyL[laser].GetN(), gL.Eval(thRef), vov_set)
        g_Vov_vs_energyL[laser].SetPointError(g_Vov_vs_energyL[laser].GetN()-1, gL.GetErrorY(3), 0)
    
        g_Vov_vs_energyR[laser].SetPoint(g_Vov_vs_energyR[laser].GetN(), gR.Eval(thRef), vov_set)
        g_Vov_vs_energyR[laser].SetPointError(g_Vov_vs_energyR[laser].GetN()-1, gR.GetErrorY(3), 0)


c = ROOT.TCanvas('c_energy_vs_Vov_noDCR','c_energy_vs_Vov_noDCR', 700, 600)
hPad = ROOT.gPad.DrawFrame(0.,0.,4.,1000.)
hPad.SetTitle(";V_{OV}^{eff} [V]; energy [ADC]")
hPad.Draw()
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
for i,laser in enumerate(laserTunes):
    g_energyL_vs_Vov[laser].SetMarkerStyle(20)
    g_energyR_vs_Vov[laser].SetMarkerStyle(24)
    g_energyL_vs_Vov[laser].SetMarkerColor(51+i*8)
    g_energyR_vs_Vov[laser].SetMarkerColor(51+i*8)
    g_energyL_vs_Vov[laser].SetLineColor(51+i*8)
    g_energyR_vs_Vov[laser].SetLineColor(51+i*8)
    g_energyL_vs_Vov[laser].Draw('plsame')
    g_energyR_vs_Vov[laser].Draw('plsame')
    c.SaveAs(outdir+c.GetName()+'.png')


c = ROOT.TCanvas('c_Npe_vs_laserTune_noDCR','c_Npe_vs_laserTune_noDCR', 700, 600)
hPad = ROOT.gPad.DrawFrame(73.,0.,80.,30000.)
hPad.SetTitle('; laser tune; Npe')
hPad.Draw()
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
for i,vov_set in enumerate(setVovs):
    g_Npe_vs_laserTune[vov_set].SetMarkerColor(51+i*5)
    g_Npe_vs_laserTune[vov_set].SetLineColor(51+i*5)
    g_Npe_vs_laserTune[vov_set].Draw('plsame')
c.SaveAs(outdir+c.GetName()+'.png')


c = ROOT.TCanvas('c_Npe_vs_energy_noDCR','c_Npe_vs_energy_noDCR', 700, 600)
hPad = ROOT.gPad.DrawFrame(0.,0.,1024.,30000.)
hPad.SetTitle('; energy[ADC]; Npe')
hPad.Draw()
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
for i,vov_set in enumerate(setVovs):
    g_Npe_vs_energy[vov_set].SetMarkerColor(51+i*5)
    g_Npe_vs_energy[vov_set].SetLineColor(51+i*5)
    g_Npe_vs_energy[vov_set].Draw('plsame')
c.SaveAs(outdir+c.GetName()+'.png')



# get energy, time resolution with DCR + estimate VovEff from energy
g_Npe_vs_VovEff_withDCR = {} 
g_tRes_vs_Npe_withDCR = {} 
for vov_set in setVovs:
    g_tRes_vs_Npe_withDCR[vov_set] = ROOT.TGraphErrors()



vovEff  = {}
vovEff_err  = {}

for laser in laserTunes:
    g_Npe_vs_VovEff_withDCR[laser] = ROOT.TGraphErrors()
    for vov_set in setVovs:

        fName = '../plots/summaryPlots_run%d.root'% runs_dict_DCR[vov_set,laser]
        ff = ROOT.TFile.Open(fName)
        gRes = ff.Get('g_deltaT_totRatioCorr_bestTh_vs_vov_enBin01_average')        

        #gL = ff.Get('g_energyL_vs_th_bar00_Vov%.2f'%vov_set)
        #gR = ff.Get('g_energyR_vs_th_bar00_Vov%.2f'%vov_set)
        
        #energyL = gL.Eval(thRef)
        #energyR = gR.Eval(thRef)
        #vovEffL =  g_Vov_vs_energyL[laser].Eval( energyL )
        #vovEffR =  g_Vov_vs_energyR[laser].Eval( energyR )
        #

        #if (energyL < min(g_Vov_vs_energyL[laser].GetX()) or energyR < min(g_Vov_vs_energyR[laser].GetX()) ): continue
        #if (energyL > max(g_Vov_vs_energyL[laser].GetX()) or energyR > max(g_Vov_vs_energyR[laser].GetX()) ): continue
        
        IarrayL = g_IV_laser_led_L[laser].Eval(Vbd_L+vov_set)/1000 # uA->mA
        IarrayR = g_IV_laser_led_R[laser].Eval(Vbd_R+vov_set)/1000 #uA->mA
        Iarray  = 0.5*(IarrayL+IarrayR)

        # approx. (it shoulb be evaluated in OV_eff)
        IphoL   = (g_IV_laser_L[laser].Eval(Vbd_L+vov_set) - g_IV_dark_L[laser].Eval(Vbd_L+vov_set))/1000. # uA -> mA
        IphoR   = (g_IV_laser_R[laser].Eval(Vbd_R+vov_set) - g_IV_dark_R[laser].Eval(Vbd_R+vov_set))/1000. # uA -> mA

        IchL = IphoL*fIpho + (IarrayL - IphoL) * fIdcrL
        IchR = IphoR*fIpho + (IarrayR - IphoR) * fIdcrR
        Ich  = 0.5*(IchL+IchR)
                        
        vovEffL = computeVovEff(vov_set, IarrayL, IchL)
        vovEffR = computeVovEff(vov_set, IarrayR, IchR)
        vovEff[vov_set,laser] = 0.5*(vovEffL+vovEffR)    
        vovEffL_err = 0.5*abs(computeVovEff(vov_set, IarrayL, IchL*(1-fIdcr_err)) - computeVovEff(vov_set, IarrayL, IchL*(1+fIdcr_err))) 
        vovEffR_err = 0.5*abs(computeVovEff(vov_set, IarrayR, IchR*(1-fIdcr_err)) - computeVovEff(vov_set, IarrayR, IchR*(1+fIdcr_err))) 
        vovEff_err[vov_set,laser] = 0.5 * (vovEffL_err+vovEffR_err)

        '''
        for it in range(0,10):
            IphoL   = (g_IV_laser_L[laser].Eval(Vbd_L+vovEffL) - g_IV_dark_L[laser].Eval(Vbd_L+vovEffL))/1000. # uA -> mA
            IphoR   = (g_IV_laser_R[laser].Eval(Vbd_R+vovEffR) - g_IV_dark_R[laser].Eval(Vbd_R+vovEffR))/1000. # uA -> mA
        
            IchL = IphoL*fIpho + (IarrayL - IphoL) * fIdcrL
            IchR = IphoR*fIpho + (IarrayR - IphoR) * fIdcrR
            Ich  = 0.5*(IchL+IchR)
            
            vovEffL = computeVovEff(vov_set, IarrayL, IchL)
            vovEffR = computeVovEff(vov_set, IarrayR, IchR)
            vovEff[vov_set] = 0.5*(vovEffL+vovEffR)    
        
            print(it, vov_set, vovEffL, vovEffR)
        '''

        gain = Gain(sipm_type, vovEff[vov_set, laser])
        npe = 0.5 * ( getNpe(IphoL*1E-03, laserFreq, gain, 1./fIpho)  + getNpe(IphoR*1E-03, laserFreq, gain, 1./fIpho) ) 

        g_tRes_vs_Npe_withDCR[vov_set].SetPoint( g_tRes_vs_Npe_withDCR[vov_set].GetN(), npe,  gRes.Eval(vov_set))
        g_Npe_vs_VovEff_withDCR[laser].SetPoint( g_Npe_vs_VovEff_withDCR[laser].GetN(), vovEff[vov_set, laser], npe  )



# now computre DCR contribution to the time resolution
g_tRes_DCR = {}
fitFun = {}

#for vov_set in setVovs:

for vov_eff in [0.80, 1.00, 1.20]:
    g_tRes_DCR[vov_eff] = ROOT.TGraphErrors()
    for laser in laserTunes:
        npe      = g_Npe_vs_VovEff_withDCR[laser].Eval(vov_eff)
        npe_up   = g_Npe_vs_VovEff_withDCR[laser].Eval(vov_eff+vovEff_err[vov_set])
        npe_down = g_Npe_vs_VovEff_withDCR[laser].Eval(vov_eff-vovEff_err[vov_set])
        npe_err  =  0.5*abs(npe_up-npe_down)
        npe_err  = math.sqrt( pow(npe_err,2) + pow(npe*fIpho_err,2))
        # print(laser, vov_set, vov_eff, npe, npe_err/npe)
        tRes0      = g_tRes_vs_Vov_noDCR[laser].Eval(vov_eff)
        tRes0_up   = g_tRes_vs_Vov_noDCR[laser].Eval(vov_eff-vovEff_err[vov_set] )
        tRes0_down = g_tRes_vs_Vov_noDCR[laser].Eval(vov_eff+vovEff_err[vov_set] )
        tRes0_err  = 0.5*abs(tRes0_up-tRes0_down) 
        tRes      = g_tRes_vs_Npe_withDCR[vov_set].Eval(npe)
        tRes_up   = g_tRes_vs_Npe_withDCR[vov_set].Eval(npe_down)
        tRes_down = g_tRes_vs_Npe_withDCR[vov_set].Eval(npe_up)
        tRes_err  = 0.5*(tRes_up-tRes_down)
        if (tRes < tRes0 or tRes_up < tRes0_up or tRes_down < tRes0_down) : continue
        tRes_DCR = math.sqrt( tRes*tRes - tRes0*tRes0)    
        #tRes_DCR_up = math.sqrt( tRes_up*tRes_up - tRes0_up*tRes0_up)    
        #tRes_DCR_down = math.sqrt( tRes_down*tRes_down - tRes0_down*tRes0_down)
        #tRes_DCR_err = 0.5*abs(tRes_DCR_up-tRes_DCR_down)
        tRes_DCR_err = 1./tRes_DCR* math.sqrt( pow(tRes*tRes_err, 2) + pow(tRes0*tRes0_err, 2)) 
        #print(tRes_DCR_err, 0.5*abs(tRes_DCR_up-tRes_DCR_down) )
        g_tRes_DCR[vov_set].SetPoint(g_tRes_DCR[vov_set].GetN(), npe, tRes_DCR)
        g_tRes_DCR[vov_set].SetPointError(g_tRes_DCR[vov_set].GetN()-1, npe_err, tRes_DCR_err)     
        

    c = ROOT.TCanvas('c_timeResolutionDCR_vs_Npe_VovEff%.02f'%vov_eff, 'c_timeResolutionDCR_vs_Npe_Vov%.02f'%vov_eff, 700, 600)
    hPad = ROOT.gPad.DrawFrame(0.,0.,25000.,80.)
    hPad.SetTitle("; Npe; #sigma_{t,DCR} [ps]")
    hPad.Draw()
    ROOT.gPad.SetGridx()
    ROOT.gPad.SetGridy()
    g_tRes_DCR[vov_set].Draw('psame')
    fitFun_old = ROOT.TF1('fitFun_%.2f'%vov_set, '[1]*pow(x/6000.,[0])', 3000, 10000) 
    fitFun_old.SetLineColor(ROOT.kRed-10)
    fitFun_old.SetLineStyle(7)
    fitFun_old.SetParameters(-1.0, 40.)
    fitFun[vov_set] = ROOT.TF1('fitFun_%.2f'%vov_set, '[1]*pow(x/6000.,[0])', 0, 100000)
    fitFun[vov_set].SetLineColor(2)
    fitFun[vov_set].SetParameters(-1.0, 40.)
    #fitFun[vov_set].FixParameter(0,-1)
    g_tRes_DCR[vov_set].Fit(fitFun[vov_set], 'QRS')
    fitFun_old.Draw('same')
    leg = ROOT.TLegend(0.6, 0.6, 0.89, 0.8)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.AddEntry(g_tRes_DCR[vov_set],'data - 25#mum','PL')
    leg.AddEntry(fitFun[vov_set],'fit ', 'L')
    leg.AddEntry(fitFun_old,'old param - 15#mum', 'L')
    leg.Draw()
    c.SaveAs(outdir+c.GetName()+'.png')
    c.SaveAs(outdir+c.GetName()+'.pdf')
