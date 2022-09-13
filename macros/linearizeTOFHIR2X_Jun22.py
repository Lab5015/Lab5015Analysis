#! /usr/bin/env python
import os
import shutil
import glob
import math
import array
import sys
import time
import argparse
import ROOT
import CMS_lumi, tdrstyle
from collections import OrderedDict

#parser = argparse.ArgumentParser(description='Module characterization summary plots')
##parser.add_argument("-r",  "--runs",          required=True, type=str, help="comma-separated list of runs to be processed")
#parser.add_argument("-i",  "--inputLabels",   required=True, type=str, help="comma-separated list of input labels")
#parser.add_argument("-m",  "--resMode",       required=True, type=int, help="resolution mode: 2 - tDiff, 1 - tAve")
#parser.add_argument("-o",  "--outFolder",     required=True, type=str, help="out folder")
#args = parser.parse_args()

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.04,'X')
ROOT.gStyle.SetLabelSize(0.04,'Y')
ROOT.gStyle.SetTitleSize(0.04,'X')
ROOT.gStyle.SetTitleSize(0.04,'Y')
ROOT.gStyle.SetTitleOffset(1.1,'X')
ROOT.gStyle.SetTitleOffset(1.2,'Y')
ROOT.gROOT.SetBatch(True)

#def Gain(ov, sipm, irr='0'):
#    k = 1.
#    if (irr == '2E14' and 'HPK' in sipm): k = 0.92 # gain reduction for HPK 2E14 irradiated SiPMs                                                                                                                                                                                                                                                                                                                                   ## 
#    if ('HPK' in sipm):
#        return k*(36890. + 97602.*ov) # HPK

def Gain(ov, sipm, irr='0'):
    k = 1.
    if (irr == '2E14' and 'HPK' in sipm): k = 0.92 # gain reduction for HPK 2E14 irradiated SiPMs
    
    if ('HPK' in sipm):
        return k*(ov+0.25)*(12.940-0.2*ov+3.2)/1.602/0.0001

def PDE(ov, sipm, irr='0'):
    k = 1.
    if (irr == '2E14' and 'HPK' in sipm): k = 0.78 # 22% PDE reduction for HPK SiPMs irradiated 2E14   
    if (irr == '1E14' and 'HPK' in sipm): k = 0.89 # 11% PDE reduction for HPK SiPMs irradiated 1E14 ?(assume that for 1E14 is half of 2E14) 
    if ('HPK' in sipm):
        return k  * 1.0228  * 0.384 * ( 1. - math.exp(-1.*0.583*ov) ) # 1.0228 factor to account for LYSO emission spectrum
    # FBK-MS
    #if ('FBK' in sipm):
    #    return k  0.8847*0.466  ( 1. - math.exp(-1.*0.314*ov) ) # 0.8847 factor to account for LYSO emission spectrum
    #FBK W4C
    if ('FBK' in sipm):
        return k  * 0.490  * ( 1. - math.exp(-1.*0.225*ov) )/1.071 # 1.071 factor to account for bech calib, convolution PDE with LYSO already accounted for

def ECF(ov):
    return 1. + 1.78409e-04*pow(ov,3.68613)
    

def GetSaturationCorrection(Ncells, Edep, LO):
    Npe    = Edep * LO
    Nfired = Ncells * (1 - math.exp(-Npe/Ncells))
    k = Nfired / Npe
    return k


thRef = 7
list_bars = [it for it in range(0, 16)]
list_Vovs = [1.5, 3.5]
list2_Vovs = [1.5, 2.5, 3.5, 5.0]
list_sides = ['L', 'R',]

goodBars = ['00L', '02L', '03R', '04R', '05L', '05R', '06L', '07L', '08R', '09R', '10R', '11R', '12L', '12R', '13L', '14L', '14R']

colors = { 1.5:ROOT.kRed, 2.5:ROOT.kOrange, 3.5:ROOT.kBlue, 5.0:ROOT.kGreen+1, 7.5:ROOT.kTeal }

infiles = OrderedDict()
infiles[(1.5, 0)]  = ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/simona_TB_CERN/Lab5015Analysis/plots/moduleCharacterization_step2_HPK_nonIrr_LYSO528_T10C_Vov1.50_angle_0.root')
infiles[(1.5, 25)] = ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/simona_TB_CERN/Lab5015Analysis/plots/moduleCharacterization_step2_HPK_nonIrr_LYSO528_T10C_Vov1.50_angle_25.root')
infiles[(1.5, 52)] = ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/simona_TB_CERN/Lab5015Analysis/plots/moduleCharacterization_step2_HPK_nonIrr_LYSO528_T10C_Vov1.50_ith2_8.root')
infiles[(3.5, 0)]  = ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/simona_TB_CERN/Lab5015Analysis/plots/moduleCharacterization_step2_HPK_nonIrr_LYSO528_T10C_Vov3.50_angle_0.root')
infiles[(3.5, 25)] = ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/simona_TB_CERN/Lab5015Analysis/plots/moduleCharacterization_step2_HPK_nonIrr_LYSO528_T10C_Vov3.50_angle_25.root')
infiles[(3.5, 52)] = ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/simona_TB_CERN/Lab5015Analysis/plots/moduleCharacterization_step2_HPK_nonIrr_LYSO528_T10C_Vov3.50_ith2_8.root')

infiles2 = OrderedDict()
infiles2[(1.5, 52)] = ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/martina_TB_CERN_June22/Lab5015Analysis/plots/moduleCharacterization_step2_HPK_nonIrr_LYSO528_T10C_Vov1.50.root','READ')
infiles2[(2.5, 52)] = ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/martina_TB_CERN_June22/Lab5015Analysis/plots/moduleCharacterization_step2_HPK_nonIrr_LYSO528_T10C_Vov2.50.root','READ')
infiles2[(3.5, 52)] = ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/martina_TB_CERN_June22/Lab5015Analysis/plots/moduleCharacterization_step2_HPK_nonIrr_LYSO528_T10C_Vov3.50.root','READ')
infiles2[(5.0, 52)] = ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/martina_TB_CERN_June22/Lab5015Analysis/plots/moduleCharacterization_step2_HPK_nonIrr_LYSO528_T10C_Vov5.00.root','READ')

c1 = {}
g_linearization_angle = {}
f_linearization_angle = {}
g_linearization2_angle = {}
f_linearization2_angle = {}

outfile = ROOT.TFile('linearizeTOFHIR2X_Jun22.root','RECREATE')

for bar in list_bars:
    for side in list_sides:
        
        barLabel = '%02d%s'%(bar,side)
        if barLabel not in goodBars:
            continue
        
    	leg = ROOT.TLegend(0.50, 0.70, 0.89, 0.90)
    	leg.SetBorderSize(0)
    	leg.SetFillStyle(0)
    	c1[bar] = ROOT.TCanvas('c1','',1400,600)
    	c1[bar].Divide(2,1)
    	c1[bar].cd(1)
    	ROOT.gPad.SetLogy()
    	hPad = ROOT.gPad.DrawFrame(-24.,1.,1000.,100000.)
    	hPad.GetXaxis().SetTitle('energy [ADC]')
    	hPad.GetYaxis().SetTitle('entries')
    	hPad.Draw()
    	lat = ROOT.TLatex( 0.17, 0.96, 'bar %02d %s'%(bar,side))
    	lat.SetNDC()
    	lat.SetTextSize(0.040)
    	lat.Draw('same')
        
        for Vov in list_Vovs:
            g_linearization_angle[bar,Vov] = ROOT.TGraphErrors()
        g_linearization2_angle[bar] = ROOT.TGraphErrors()
        
    	for key in infiles.keys():
    	    Vov = key[0]
    	    angle = key[1]
    	    energy = 0.86 * 3. / math.cos(angle*3.14159/180.) # 0.86 MeV/mm * 3 mm / cos(theta)
            
            
    	    histo = infiles[key].Get('h1_energy_bar%02d%s_Vov%.2f_th%02d'%(bar,side,Vov,thRef))
    	    if not histo:
    	        continue
    	    histo.SetLineColor(colors[Vov])
    	    histo.Draw('same')
    	    leg.AddEntry(histo,'V_{OV} = %.2f V, angle = %d deg'%(Vov,angle),'L')
    	    
    	    fitFunc = histo.GetFunction('f_landau_bar%02d%s_Vov%.2f_vth1_%02.0f'%(bar,side,Vov,thRef))
    	    if not fitFunc:
    	        continue
    	    fitFunc.SetLineColor(colors[Vov])
    	    
    	    g_linearization_angle[bar,Vov].SetPoint(g_linearization_angle[bar,Vov].GetN(),fitFunc.GetParameter(1),energy*Gain(Vov,'HPK','0E14')*PDE(Vov,'HPK','0E14')*ECF(Vov)*GetSaturationCorrection(40000.,energy,1250*PDE(Vov,'HPK','0E14')/PDE(1.5,'HPK','0E14'))/(Gain(1.5,'HPK','0E14')*PDE(1.5,'HPK','0E14')*ECF(1.5)))
    	    g_linearization_angle[bar,Vov].SetPointError(g_linearization_angle[bar,Vov].GetN()-1,fitFunc.GetParError(1),0.)
            
        
    	for key in infiles2.keys():
    	    Vov = key[0]
    	    angle = key[1]
    	    energy = 0.86 * 3. / math.cos(angle*3.14159/180.) # 0.86 MeV/mm * 3 mm / cos(theta)
            
            
    	    histo = infiles2[key].Get('h1_energy_bar%02d%s_Vov%.2f_th%02d'%(bar,side,Vov,thRef))
    	    if not histo:
    	        continue
    	    histo.SetLineColor(colors[Vov])
    	    histo.Draw('same')
    	    leg.AddEntry(histo,'V_{OV} = %.2f V, angle = %d deg'%(Vov,angle),'L')
    	    
    	    fitFunc = histo.GetFunction('f_landau_bar%02d%s_Vov%.2f_vth1_%02.0f'%(bar,side,Vov,thRef))
    	    if not fitFunc:
    	        continue
    	    fitFunc.SetLineColor(colors[Vov])
    	    
    	    g_linearization2_angle[bar].SetPoint(g_linearization2_angle[bar].GetN(),fitFunc.GetParameter(1),energy*Gain(Vov,'HPK','0E14')*PDE(Vov,'HPK','0E14')*ECF(Vov)*GetSaturationCorrection(40000.,energy,1250*PDE(Vov,'HPK','0E14')/PDE(1.5,'HPK','0E14'))/(Gain(1.5,'HPK','0E14')*PDE(1.5,'HPK','0E14')*ECF(1.5)))
    	    g_linearization2_angle[bar].SetPointError(g_linearization2_angle[bar].GetN()-1,fitFunc.GetParError(1),0.)
            
    	leg.Draw('same')
    	
    	c1[bar].cd(2)
    	hPad = ROOT.gPad.DrawFrame(0.,0.,800.,25.)
    	hPad.GetXaxis().SetTitle('energy [ADC]')
    	hPad.GetYaxis().SetTitle('E_{MeV} #times G #times PDE #times ECF #times k_{sat} / [G #times PDE #times ECF]_{3.5 V}')
    	hPad.Draw()
    	lat.Draw('same')
        
        for Vov in list_Vovs:
            g_linearization_angle[bar,Vov].SetMarkerStyle(20)
            g_linearization_angle[bar,Vov].SetMarkerSize(1.2)
            g_linearization_angle[bar,Vov].SetMarkerColor(colors[Vov])
            g_linearization_angle[bar,Vov].Draw('P,same')
            f_linearization_angle[bar,Vov] = ROOT.TF1('f_linearization_angle_bar%02d%s_Vov%.02f'%(bar,side,Vov),"pol1",0.,1000.)
            f_linearization_angle[bar,Vov].SetParameters(1,0.002)
            g_linearization_angle[bar,Vov].Fit(f_linearization_angle[bar,Vov],"QNS","",20.,1000.)
            f_linearization_angle[bar,Vov].SetLineColor(colors[Vov])
            f_linearization_angle[bar,Vov].SetLineStyle(2)
            f_linearization_angle[bar,Vov].SetLineWidth(1)
            f_linearization_angle[bar,Vov].Draw('same')
            outfile.cd()
            f_linearization_angle[bar,Vov].Write()

        g_linearization2_angle[bar].SetMarkerStyle(20)
        g_linearization2_angle[bar].SetMarkerSize(1.2)
        g_linearization2_angle[bar].SetMarkerColor(ROOT.kBlack)
        g_linearization2_angle[bar].Draw('P,same')
        f_linearization2_angle[bar] = ROOT.TF1('f_linearization_Vov_bar%02d%s'%(bar,side),"pol1",0.,1000.)
        f_linearization2_angle[bar].SetParameters(1,0.002)
        g_linearization2_angle[bar].Fit(f_linearization2_angle[bar],"QNS","",20.,1000.)
        f_linearization2_angle[bar].SetLineColor(ROOT.kBlack)
        f_linearization2_angle[bar].SetLineStyle(2)
        f_linearization2_angle[bar].SetLineWidth(1)
        f_linearization2_angle[bar].Draw('same')
        f_linearization2_angle[bar].Write()
        
    	#g_linearization[bar].SetMarkerStyle(20)
    	#g_linearization[bar].SetMarkerSize(1.2)
    	#g_linearization[bar].Draw('P,same')
    	#f_linearization[bar] = ROOT.TF1('f_linearization_bar%02d%s'%(bar,side),"pol1",0.,1000.)
    	#f_linearization[bar].SetParameters(1,0.002)
    	#g_linearization[bar].Fit(f_linearization[bar],"QNRS")
    	#f_linearization[bar].SetLineColor(ROOT.kBlack)
    	#f_linearization[bar].SetLineStyle(2)
    	#f_linearization[bar].SetLineWidth(1)
    	#f_linearization[bar].Draw('same')
        
    	leg2 = ROOT.TLegend(0.25, 0.75, 0.50, 0.85)
    	leg2.SetBorderSize(0)
    	leg2.SetFillStyle(0)
    	#leg2.AddEntry(g_linearization[bar],'MIP Landau MPV','PE')
    	leg2.Draw()
    	c1[bar].Print('c1_linearization_bar%02d%s.png'%(bar,side))
