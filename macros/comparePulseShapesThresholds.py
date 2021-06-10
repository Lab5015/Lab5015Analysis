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
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetTitleOffset(1.05,'Y')
ROOT.gErrorIgnoreLevel = ROOT.kWarning;

gain = int(sys.argv[1])

#files
fname = {}

# preAmpGain1
runs = {'vth2' : 966, 'vth1_3' : 967}

# preAmpGain0
if (gain == 0 ):
    runs = {'vth2' : 964, 'vth1_3' : 965}

for t in runs:
    fname[t]  = '/data/Lab5015Analysis/pulseShapes/pulseShape_run0%d.root'%runs[t] 
    
# graphs
f = {}
g_amp = {} # g[th]

for t in runs:
    f[t] = ROOT.TFile.Open(fname[t])
    g_amp[t] = f[t].Get('g_ps_totSel_ch1_Vov3.5')
    
fitSR = {}
tt = {}
fitSRmax = {}
npoints = 3

c1 = ROOT.TCanvas('c_pulseShape_preAmpGain%d_vth1-vth2', 'c_pulseShape_preAmpGain%d_vth1-vth2'%gain, 600, 600)
c1.SetGridx()
c1.SetGridy()
hdummy = ROOT.TH2F('hdummy','', 100, g_amp['vth2'].GetX()[0]-0.5, g_amp['vth2'].GetX()[0]+19.5, 100, -30., 500)
hdummy.GetXaxis().SetTitle('time (ns)')
hdummy.GetYaxis().SetTitle('amplitude (mV)')
hdummy.Draw()

for i, t in enumerate(runs):
    g_amp[t].SetMarkerStyle(20)
    g_amp[t].SetMarkerSize(0.7)
    g_amp[t].SetMarkerColor(51+i*8)
    g_amp[t].SetLineColor(51+i*8)
    g_amp[t].Draw('psame')
    #leg.AddEntry(g_amp[t] , 'n_{TAPS} = %d'%t, 'PL')
    
    
    #fitSR[t] = ROOT.TF1('fitSR_%s'%t, 'pol1', 0, 10)
    #fitSR[t].SetLineColor(g_amp[t].GetLineColor())
    #fitSR[t].SetLineWidth(3)
    ##fitSR[t].SetRange(g_amp[t].GetX()[12], g_amp[t].GetX()[18] )
    #fitSR[t].SetRange(g_amp[t].GetX()[8], g_amp[t].GetX()[12] )
    #fitSR[t].SetParameters(0, 50)
    #g_amp[t].Fit(fitSR[t], 'QRS+')
    #print 'N_taps %s : slew rate = %f'%(t,fitSR[t].GetParameter(1) )

    # get max slew rate
    # -- SR max
    fitSRmax[t] = ROOT.TF1('fitSRmax_%s'%t, 'pol1', 0, 10)
    fitSRmax[t].SetLineColor(g_amp[t].GetLineColor())
    fitSRmax[t].SetLineWidth(3)
    fitSRmax[t].SetParameters(0, 50)
    maxSR = 0
    jmax = 0
    for j in range(npoints, g_amp[t].GetN()/2-(npoints+1)):
        fitSRmax[t].SetRange(g_amp[t].GetX()[j-npoints], g_amp[t].GetX()[j+npoints] )
        g_amp[t].Fit(fitSRmax[t], 'QRN')
        #print '%s : slew rate = %f'%(d,fitSRmax[t].GetParameter(1) )
        if ( fitSRmax[t].GetParameter(1) > maxSR ):
            maxSR = fitSRmax[t].GetParameter(1)
            jmax = j
    fitSRmax[t].SetRange(g_amp[t].GetX()[jmax-npoints], g_amp[t].GetX()[jmax+npoints] )
    g_amp[t].Fit(fitSRmax[t], 'QRS+')
    print 'Threshold %s : slew rate = %f'%(t,fitSRmax[t].GetParameter(1) )
    
    tt[t] = ROOT.TLatex()
    tt[t].SetNDC()
    tt[t].SetTextColor( g_amp[t].GetLineColor() )
    tt[t].SetTextSize( 0.030 )
    tt[t].DrawLatex( 0.20, 0.90-i*0.03, '%s:  max slew rate = %0.1f mV/ns'%(t,fitSRmax[t].GetParameter(1)) )

for c in [c1]:
    c.SaveAs(c.GetName()+'.pdf')
    c.SaveAs(c.GetName()+'.png')

raw_input('ok?')


