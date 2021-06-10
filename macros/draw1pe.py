#! /usr/bin/env python
import os
import shutil
import glob
import math
import array
import sys
import time


import ROOT
import tdrstyle

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetTitleOffset(1.05,'Y')
ROOT.gStyle.SetLabelSize(0.04)
ROOT.gErrorIgnoreLevel = ROOT.kWarning;
ROOT.gROOT.SetBatch(True)

plotDir = '/var/www/html/TOFHIR2A/singlePhotoelectronMeasurement/'

dac_to_mV = {'vth2': 8.,
             'vth1_4': 4.,
             'vth1_3': 2.,
             'vth1_1': 1.,
             'vth1_0': 0.5}

runs_dict = { 964 : [ 0, 'vth2'  , 613.95],
              965 : [ 0, 'vth1_3', 613.95],
              966 : [ 1, 'vth2'  , 613.95],
              967 : [ 1, 'vth1_3', 613.95],

              968 : [ 0, 'vth2'  , 301.45],
              969 : [ 0, 'vth1_3', 301.45],
              970 : [ 1, 'vth2'  , 301.45],
              971 : [ 1, 'vth1_3', 301.45],

              972 : [ 0, 'vth2'  , 158.04],
              973 : [ 0, 'vth1_3', 158.04],
              974 : [ 1, 'vth2'  , 158.04],
              975 : [ 1, 'vth1_3', 158.04],

              976 : [ 0, 'vth2'  , 477.38],
              977 : [ 0, 'vth1_3', 477.38],
              978 : [ 1, 'vth2'  , 477.38],
              979 : [ 1, 'vth1_3', 477.38],

              986 : [ 0, 'vth2'  , 67.64],
              987 : [ 0, 'vth1_3', 67.64],
              988 : [ 1, 'vth2'  , 67.64],
              989 : [ 1, 'vth1_3', 67.64],
} 
  
g_1pe_vs_npe = {}
g_amp_vs_npe = {}
for gain in [0,1]:
    for th in ['vth2','vth1_3']:
        g_1pe_vs_npe[gain, th] = ROOT.TGraphErrors()
        g_amp_vs_npe[gain, th] = ROOT.TGraphErrors()

for run in runs_dict:
    print run
    gain = runs_dict[run][0]
    th   = runs_dict[run][1]
    npe  = runs_dict[run][2]

    f = ROOT.TFile.Open('/data/Lab5015Analysis/pulseShapes/pulseShape_run%04d.root'%run)
    g = f.Get('g_N_ch1_Vov3.5')
    if (g == None): continue

    f_sigmoid = ROOT.TF1("f_sigmoid","[0]*(1-0.5*(1.+TMath::Erf((x-[1])/[2])))",0.,64.)
    f_sigmoid.SetNpx(10000)
    index = 0
    for j in range(0, g.GetN()):
        if (g.GetY()[j] < g.GetY()[0]/2):
            index = j
            break
    f_sigmoid.SetParameters(g.GetY()[0], g.GetX()[index],1.)
    if (th == 'vth1_3'):
        f_sigmoid.SetRange(3,63)
    if (run  == 989):
        f_sigmoid.SetRange(3,63)
        f_sigmoid.SetParameters(g.GetY()[3], g.GetX()[index],1.)

    g.Fit(f_sigmoid,"QRS");
    g.Fit(f_sigmoid,"QRS");


    ctemp = ROOT.TCanvas('c_N_vs_th_run%04d'%run)
    g.GetXaxis().SetTitle('th [DAC]')
    g.GetYaxis().SetTitle('N')
    g.Draw('ap')
    ctemp.SaveAs(plotDir+ctemp.GetName()+'.png')

    if (f_sigmoid.GetNDF()==0): continue
    if (f_sigmoid.GetChisquare()/f_sigmoid.GetNDF() < 1): continue
    if (f_sigmoid.GetChisquare()/f_sigmoid.GetNDF() > 999): continue

    print f_sigmoid.GetChisquare(), f_sigmoid.GetNDF()
    
    #if ( f_sigmoid.GetChisquare()/f_sigmoid.GetNDF() > 200): continue

    amp = f_sigmoid.GetParameter(1)*dac_to_mV[th]
    amp_err = f_sigmoid.GetParError(1)*dac_to_mV[th]
    
    print gain,th, npe, amp, amp/npe
    g_1pe_vs_npe[gain,th].SetPoint(g_1pe_vs_npe[gain,th].GetN(), npe,  amp/npe)
    g_1pe_vs_npe[gain,th].SetPointError(g_1pe_vs_npe[gain,th].GetN()-1, 0.02*npe,  amp_err/npe)

    g_amp_vs_npe[gain,th].SetPoint(g_amp_vs_npe[gain,th].GetN(), npe,  amp)
    g_amp_vs_npe[gain,th].SetPointError(g_amp_vs_npe[gain,th].GetN()-1, 0.02*npe,  amp_err)
    


g_amp_vs_npe_Tahereh = ROOT.TGraph('Npe_Tahereh.txt')

legend = ROOT.TLegend(0.18, 0.65, 0.45, 0.92)


c1 = ROOT.TCanvas('c_1pe_vs_Npe','c_1pe_vs_Npe')
c1.SetGridx()
c1.SetGridy()
hdummy1 = ROOT.TH2F('hdummy','', 100, 0, 700, 100, 0., 0.75)
hdummy1.GetXaxis().SetTitle('N_{pe}')
hdummy1.GetYaxis().SetTitle('1 pe amplitude (mV/pe)')
hdummy1.Draw()

for gain in [0,1]:
    for th in ['vth2','vth1_3']:
        g_1pe_vs_npe[gain, th].SetLineColor(4-gain)
        g_1pe_vs_npe[gain, th].SetMarkerColor(4-gain)
        if (th == 'vth1_3'): g_1pe_vs_npe[gain, th].SetMarkerStyle(24)
        g_1pe_vs_npe[gain, th].Draw('p')
        legend.AddEntry(g_1pe_vs_npe[gain, th], 'preAmpGain%d - %s'%(gain, th))
legend.Draw('same')

c2 = ROOT.TCanvas('c_amp_vs_Npe','c_amp_vs_Npe')
c2.SetGridx()
c2.SetGridy()
hdummy2 = ROOT.TH2F('hdummy2','', 100, 0, 700, 100, 0., 400)
hdummy2.GetXaxis().SetTitle('N_{pe}')
hdummy2.GetYaxis().SetTitle('amp (mV)')
hdummy2.Draw()

for gain in [0,1]:
    for th in ['vth2','vth1_3']:
        g_amp_vs_npe[gain, th].SetLineColor(4-gain)
        g_amp_vs_npe[gain, th].SetMarkerColor(4-gain)
        if (th == 'vth1_3'): g_amp_vs_npe[gain, th].SetMarkerStyle(24)
        g_amp_vs_npe[gain, th].Draw('p')
g_amp_vs_npe_Tahereh.SetMarkerSize(1)
g_amp_vs_npe_Tahereh.SetMarkerStyle(22)
g_amp_vs_npe_Tahereh.SetLineColor(ROOT.kGray+2)
g_amp_vs_npe_Tahereh.SetLineWidth(1)
g_amp_vs_npe_Tahereh.Draw('plsame')
legend.Draw('same')

for c in c1, c2:
    c.SaveAs(plotDir+c.GetName()+'.png')
    c.SaveAs(plotDir+c.GetName()+'.pdf')
    c.SaveAs(plotDir+c.GetName()+'.C')

raw_input('ciao')
