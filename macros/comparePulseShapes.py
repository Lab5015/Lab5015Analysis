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

if (len(sys.argv) < 2):
    sys.exit('Usage:  comparePulseShape.py gain[0/1]')

gain = int(sys.argv[1])
print 'preAmpGain = ', gain

#files
fname = {}
fname['sim']  = '/data/Lab5015Analysis/pulseShapes/pulseShape_sim.root' # from Tahereh, 9500 pe, preAmpGain1
fname['data'] = '/data/Lab5015Analysis/pulseShapes//pulseShape_run0762.root' # preAmpGain1_vth2, 9500 pe, 3.5 OV
#fname['data'] = '/data/Lab5015Analysis/pulseShapes//pulseShape_run0814.root' # preAmpGain1_vth2
#fname['data'] = '/data/Lab5015Analysis/pulseShapes//pulseShape_run0887.root' # preAmpGain1_vth2
#fname['data'] = '/data/Lab5015Analysis/pulseShapes//pulseShape_run0913.root' # preAmpGain1_vth2
fname['data'] = '/data/Lab5015Analysis/pulseShapes//pulseShape_run0768.root' # preAmpGain1_vth2, 1.2 OV
fname['psi']  = '/data/Lab5015Analysis/pulseShapes/pulseShape_run2898_bar08.root' # preAmpGain1_vth2

if (gain == 0):
    fname['data'] = '/data/Lab5015Analysis/pulseShapes//pulseShape_run0895.root' # preAmpGain0_vth2, 9500 pe , 3.5 OV
    #fname['data'] = '/data/Lab5015Analysis/pulseShapes//pulseShape_run0910.root' # preAmpGain0_vth2, ~ 18000 pe, 3.5 OV
    #fname['data'] = '/data/Lab5015Analysis/pulseShapes//pulseShape_run0912.root' # preAmpGain0_vth2, ~ 9500 pe, 6.5 OV

gname = {'sim'  : 'g_pulseShape_preAmpGain%d'%gain,
         'data' : 'g_ps_totSel_ch1_Vov3.5',
         #'data' : 'g_ps_totSel_ch1_Vov1.5',
         'psi'  : 'g_ps_totSel_deltaT_ch2_Vov3.5'}


# graphs
f = {}
g_amp = {} # g[sim/data]
t_ref = {}
a_ref = {}

for d in ['sim', 'data', 'psi']:
    f[d] = ROOT.TFile.Open(fname[d])
    g_amp[d] = f[d].Get(gname[d])
    #ftemp = g_amp[d].GetFunction('fitFun')
    #if (ftemp): ftemp.Delete()
    
    if (d == 'sim'):
        #baseline = g_amp['sim'].GetY()[0]
        baseline = 56.
        for i in range(0, g_amp['sim'].GetN()):
            g_amp['sim'].SetPoint(i, g_amp['sim'].GetX()[i], g_amp['sim'].GetY()[i] - baseline ) # subtract baseline
    if (d == 'data'):
        #baseline = g_amp['sim'].GetY()[0]
        baseline = 0.
        for i in range(0, g_amp['data'].GetN()):
            g_amp['data'].SetPoint(i, g_amp['data'].GetX()[i], g_amp['data'].GetY()[i] - baseline ) # subtract baseline
            print g_amp['data'].GetX()[i], g_amp['data'].GetY()[i]/8.
    if (d == 'psi'):
        #baseline = g_amp['sim'].GetY()[0]
        baseline = 0.
        for i in range(0, g_amp['psi'].GetN()):
            g_amp['psi'].SetPoint(i, g_amp['psi'].GetX()[i], g_amp['psi'].GetY()[i] - baseline ) # subtract baseline


t_offset = g_amp['sim'].GetX()[0]-g_amp['data'].GetX()[0]
#a_offset = g_amp['sim'].GetY()[0]-g_amp['data'].GetY()[0]
a_offset = 0
for i in range(0, g_amp['sim'].GetN()):
        g_amp['sim'].SetPoint(i, g_amp['sim'].GetX()[i]-t_offset, g_amp['sim'].GetY()[i] - a_offset)

t_offset = g_amp['psi'].GetX()[0]-g_amp['data'].GetX()[0]
#a_offset = g_amp['sim'].GetY()[0]-g_amp['data'].GetY()[0]
a_offset = 0
for i in range(0, g_amp['psi'].GetN()):
        g_amp['psi'].SetPoint(i, g_amp['psi'].GetX()[i]-t_offset, g_amp['psi'].GetY()[i] - a_offset)


c = ROOT.TCanvas('c_pulseShape_preAmpGain%d'%gain, 'c_pulseShape_preAmpGain%d'%gain, 600, 600)
c.SetGridx()
c.SetGridy()
#hdummy = ROOT.TH2F('hdummy','', 100, g_amp['data'].GetX()[0]-0.5, g_amp['data'].GetX()[0]+19.5, 100, -30., 600)
hdummy = ROOT.TH2F('hdummy','', 100, g_amp['data'].GetX()[0]-0.5, g_amp['data'].GetX()[0]+19.5, 100, -30., 200)
hdummy.GetXaxis().SetTitle('time (ns)')
hdummy.GetYaxis().SetTitle('amplitude (mV)')
hdummy.Draw()
g_amp['sim'].SetMarkerStyle(24)
g_amp['sim'].SetMarkerSize(0.8)
g_amp['sim'].SetMarkerColor(1)
#g_amp['sim'].Draw('psame')
g_amp['data'].SetMarkerStyle(24)
g_amp['data'].SetMarkerSize(0.8)
g_amp['data'].Draw('psame')
g_amp['psi'].SetMarkerColor(3)
g_amp['psi'].SetMarkerStyle(24)
g_amp['psi'].SetMarkerSize(0.7)
#g_amp['psi'].Draw('psame')

leg = ROOT.TLegend(0.68, 0.82, 0.93, 0.94)
leg.AddEntry(g_amp['sim'] , 'simulation', 'PL')
leg.AddEntry(g_amp['data'], 'measurement', 'PL')
leg.Draw('same')

# --- slew rate
# -- at beginning of the pulse 
fitSR = {}
t = {}

fitSRmax = {}
tmax = {}
npoints = 1

fitSR1 = {}
t1 = {}

#for i, d in enumerate(['sim','data']):
for i, d in enumerate(['data']):
    # -- SR at beginning of the pulse  
    fitSR[d] = ROOT.TF1('fitSR_%s'%d, 'pol1', 0, 30)
    fitSR[d].SetLineColor(g_amp[d].GetLineColor())
    fitSR[d].SetLineWidth(3)
    fitSR[d].SetRange(g_amp[d].GetX()[0], g_amp[d].GetX()[4] )
    fitSR[d].SetParameters(0, 50)
    g_amp[d].Fit(fitSR[d], 'QRS+')
    print '%s : slew rate at timing thr. = %f'%(d,fitSR[d].GetParameter(1) )
    t[d] = ROOT.TLatex()
    t[d].SetNDC()
    t[d].SetTextColor( g_amp[d].GetLineColor() )
    t[d].SetTextSize( 0.028 )
    t[d].DrawLatex( 0.20, 0.91-i*0.03, 'slew rate at the timing thr. = %0.1f mV/ns'%fitSR[d].GetParameter(1) )

    # -- SR max
    fitSRmax[d] = ROOT.TF1('fitSRmax_%s'%d, 'pol1', 0, 30)
    fitSRmax[d].SetLineColor(g_amp[d].GetLineColor())
    fitSRmax[d].SetLineWidth(3)
    fitSRmax[d].SetParameters(0, 50)
    maxSR = 0
    jmax = 0
    for j in range(npoints, g_amp[d].GetN()/2-(npoints+1)):
        fitSRmax[d].SetRange(g_amp[d].GetX()[j-npoints], g_amp[d].GetX()[j+npoints] )
        g_amp[d].Fit(fitSRmax[d], 'QRN')
        #print '%s : slew rate = %f'%(d,fitSRmax[d].GetParameter(1) )
        if ( fitSRmax[d].GetParameter(1) > maxSR ):
            maxSR = fitSRmax[d].GetParameter(1)
            jmax = j
    fitSRmax[d].SetRange(g_amp[d].GetX()[jmax-npoints], g_amp[d].GetX()[jmax+npoints] )
    g_amp[d].Fit(fitSRmax[d], 'QRS+')
    print '%s : slew rate = %f'%(d,fitSRmax[d].GetParameter(1) )
    tmax[d] = ROOT.TLatex()
    tmax[d].SetNDC()
    tmax[d].SetTextColor( g_amp[d].GetLineColor() )
    tmax[d].SetTextSize( 0.028 )
    tmax[d].DrawLatex( 0.20, 0.83-i*0.03, 'max slew rate = %0.1f mV/ns'%maxSR)


    # -- SR after 1 ns
    fitSR1[d] = ROOT.TF1('fitSR1_%s'%d, 'pol1', 0, 30)
    fitSR1[d].SetLineColor(g_amp[d].GetLineColor())
    fitSR1[d].SetLineWidth(3)
    fitSR1[d].SetRange(0.5, 1.5)
    #fitSR1[d].SetRange(1.5, 2.5)
    fitSR1[d].SetParameters(0, 50)
    g_amp[d].Fit(fitSR1[d], 'QRN')
    print '%s : slew rate at 1 ns = %f'%(d,fitSR1[d].GetParameter(1) )
    t1[d] = ROOT.TLatex()
    t1[d].SetNDC()
    t1[d].SetTextColor( g_amp[d].GetLineColor() )
    t1[d].SetTextSize( 0.025 )
    #t1[d].DrawLatex( 0.20, 0.72-i*0.04, 'slew rate at 1 ns = %0.1f mV/ns'%fitSR1[d].GetParameter(1) )


 

c.SaveAs(c.GetName()+'.pdf')
c.SaveAs(c.GetName()+'.png')
raw_input('ok?')


