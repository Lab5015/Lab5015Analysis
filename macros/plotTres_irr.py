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
ROOT.gROOT.SetBatch(True)
#ROOT.gROOT.SetBatch(False)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0111)

from SiPM import *
from VovsEff import *

def getSlewRateFromPulseShape(g1, timingThreshold, npoints, gtemp, canvas=None):
    #if ( g1.GetN()/2 < npoints): return (-1, -1)
    if ( g1.GetN()/2 < npoints): npoints = g1.GetN()/2 
    # find index at the timing threshold
    itiming = 0
    for i in range(0,g1.GetN()):
        if (round(g1.GetY()[i]/0.313) == timingThreshold):
            itiming = i
            break

    ifirst = ROOT.TMath.LocMin(g1.GetN(), g1.GetX())
    imin = max(0, itiming-2)
    if ( imin >= 0 and g1.GetX()[imin+1] < g1.GetX()[imin] ): imin = ifirst
    tmin = g1.GetX()[imin]
    tmax = min(g1.GetX()[imin+npoints],3.)
    for i in range(imin, imin+npoints+1):
        gtemp.SetPoint(gtemp.GetN(), g1.GetX()[i], g1.GetY()[i])
        gtemp.SetPointError(gtemp.GetN()-1, g1.GetErrorX(i), g1.GetErrorY(i))
    fitSR = ROOT.TF1('fitSR', 'pol1', tmin, tmax)
    fitSR.SetLineColor(g1.GetMarkerColor()+1)
    fitSR.SetRange(tmin,tmax)
    fitSR.SetParameters(0, 10)
    fitStatus = int(gtemp.Fit(fitSR, 'QRS+'))
    sr = fitSR.Derivative( g1.GetX()[itiming])
    err_sr = fitSR.GetParError(1)
    if (canvas!=None):
        canvas.cd()
        gtemp.SetMarkerStyle(g1.GetMarkerStyle())
        gtemp.SetMarkerColor(g1.GetMarkerColor())
        gtemp.Draw('psames')
        g1.Draw('psames')
        fitSR.Draw('same')
        canvas.Update()
        #ps = g1.FindObject("stats")
        ps = gtemp.FindObject("stats")
        ps.SetTextColor(g1.GetMarkerColor())
        if ('L' in g1.GetName()):
            ps.SetY1NDC(0.85) # new y start position
            ps.SetY2NDC(0.95)# new y end position
        if ('R' in g1.GetName()):
            ps.SetY1NDC(0.73) # new y start position
            ps.SetY2NDC(0.83)# new y end position

    return(sr,err_sr)


def findTimingThreshold(g2):
    xmin = 0
    ymin = 9999
    for i in range(0, g2.GetN()):
        y = g2.GetY()[i]
        x = g2.GetX()[i]
        if ( y < ymin):
            ymin = y
            xmin = x 
    return xmin

# =====================================

outdir = '/eos/user/m/malberti/www/MTD/TOFHIR2X/MTDTB_CERN_Oct21/timeResolution_vs_Vov_2E14/'
#outdir = '/eos/user/m/malberti/www/MTD/TOFHIR2X/MTDTB_CERN_Oct21/timeResolution_vs_Vov_FBK_2E14/'
#outdir = '/eos/user/m/malberti/www/MTD/TOFHIR2X/MTDTB_CERN_Oct21/timeResolution_vs_Vov_HPK_1E13_LYSOtype1_58deg/'
#outdir = '/eos/user/m/malberti/www/MTD/TOFHIR2X/MTDTB_CERN_Oct21/timeResolution_vs_Vov_1E13/'

if (os.path.exists(outdir)==False):
    os.mkdir(outdir)
if (os.path.exists(outdir+'/plotsSR')==False):
    os.mkdir(outdir+'/plotsSR/')


outfile = ROOT.TFile.Open('plots_2E14.root','recreate')
#outfile = ROOT.TFile.Open('plots_1E13_type1.root','recreate')
#outfile = ROOT.TFile.Open('plots_1E13.root','recreate')

sipmTypes = ['HPK_2E14_T-40C','FBK_2E14_T-40C']
#sipmTypes = ['FBK_2E14_T-40C', 'FBK_2E14_T-32C', 'FBK_2E14_T-22C']
#sipmTypes = ['HPK_1E13_LYSOtype1_T0C', 'HPK_1E13_LYSOtype1_T-6C', 'HPK_1E13_LYSOtype1_T-20C', 'HPK_1E13_LYSOtype1_T-40C']
#sipmTypes = ['HPK_1E13_T0C','FBK_1E13_T0C']

fnames = {'HPK_2E14_T-40C' : '../plots/HPK_2E14_52deg_T-40C_summary.root',
          'FBK_2E14_T-40C' : '../plots/FBK_2E14_52deg_T-40C_summary.root',
          'FBK_2E14_T-32C' : '../plots/FBK_2E14_52deg_T-32C_summary.root',
          'FBK_2E14_T-22C' : '../plots/FBK_2E14_52deg_T-22C_summary.root',
          'HPK_1E13_LYSOtype1_T0C'   : '../plots/HPK_1E13_LYSOtype1_58deg_T0C_summary.root',
          'HPK_1E13_LYSOtype1_T-6C'  : '../plots/HPK_1E13_LYSOtype1_58deg_T-6C_summary.root',
          'HPK_1E13_LYSOtype1_T-20C' : '../plots/HPK_1E13_LYSOtype1_58deg_T-20C_summary.root',
          'HPK_1E13_LYSOtype1_T-40C' : '../plots/HPK_1E13_LYSOtype1_58deg_T-40C_summary.root',
          'HPK_1E13_T0C'             : '../plots/HPK_1E13_52deg_T0C_summary.root',
          'FBK_1E13_T0C'             : '../plots/FBK_1E13_52deg_T0C_summary.root'}

LO = { 'HPK_2E14_T-40C': 1100., #?
       'FBK_2E14_T-40C': 1100., #?
       'FBK_2E14_T-32C': 1100., #?
       'FBK_2E14_T-22C': 1100.,
       'HPK_1E13_T0C'  : 1100, 
       'FBK_1E13_T0C'  : 1100, 
       'HPK_1E13_LYSOtype1_T0C'  : 1100, 
       'HPK_1E13_LYSOtype1_T-6C' : 1100, 
       'HPK_1E13_LYSOtype1_T-20C': 1100, 
       'HPK_1E13_LYSOtype1_T-40C': 1100   } #?

np = 3
errSRsyst  = 0.10 # error on the slew rate

g = {}
g_Noise_vs_Vov = {}
g_Stoch_vs_Vov = {}
g_DCR_vs_Vov = {}
g_Tot_vs_Vov   = {}

g_Stoch_vs_Npe = {}
g_DCR_vs_Npe = {}

g_bestTh_vs_Vov = {}
g_SR_vs_Vov = {}

g_bestTh_vs_bar = {}
g_SR_vs_bar = {}
g_Noise_vs_bar = {}
g_Stoch_vs_bar = {}
g_DCR_vs_bar = {}

g_Tot_vs_SR   = {}


g_DCR_vs_DCRNpeSR = {}

g_DCR_vs_DCRNpeSR_all = {}

bars = {}
Vovs = {}
Npe = {}

sigma_stoch_ref = {'HPK_2E14_T-40C' : 40.,
                   'FBK_2E14_T-40C' : 47.,
                   'FBK_2E14_T-32C' : 47.,
                   'FBK_2E14_T-22C' : 47.,
                   'HPK_1E13_T0C'   : 40.,
                   'FBK_1E13_T0C'   : 47.,
                   'HPK_1E13_LYSOtype1_T0C' : 40.,
                   'HPK_1E13_LYSOtype1_T-6C' : 40.,
                   'HPK_1E13_LYSOtype1_T-20C' : 40.,
                   'HPK_1E13_LYSOtype1_T-40C' : 40.}

err_sigma_stoch_ref = 2.
ov_ref = 1.50

for sipm in sipmTypes:
    f = ROOT.TFile.Open(fnames[sipm])
    print sipm, fnames[sipm]
    listOfKeys = [key.GetName().replace('g_deltaT_energyRatioCorr_bestTh_vs_vov_','') for key in ROOT.gDirectory.GetListOfKeys() if ( 'g_deltaT_energyRatioCorr_bestTh_vs_vov_bar' in key.GetName())]
    bars[sipm] = []
    for k in listOfKeys:
        bars[sipm].append( int(k[3:5]) )

    listOfKeys2 = [key.GetName().replace('g_deltaT_energyRatioCorr_bestTh_vs_bar_','') for key in ROOT.gDirectory.GetListOfKeys() if key.GetName().startswith('g_deltaT_energyRatioCorr_bestTh_vs_bar_')]
    Vovs[sipm] = []
    for k in listOfKeys2:
        Vovs[sipm].append( float(k[3:7]) )

    print bars[sipm]
    print Vovs[sipm]
    

fPS = {}
f   = {}
for sipm in sipmTypes:
    Npe[sipm] = {}
    g[sipm] = {}
    g_Noise_vs_Vov[sipm] = {}
    g_Stoch_vs_Vov[sipm] = {}
    g_DCR_vs_Vov[sipm] = {}
    g_Tot_vs_Vov[sipm] = {}

    g_Tot_vs_SR[sipm] = {}

    g_Stoch_vs_Npe[sipm] = {}
    g_DCR_vs_Npe[sipm] = {}

    g_DCR_vs_DCRNpeSR[sipm] = {}
    g_DCR_vs_DCRNpeSR_all[sipm] = ROOT.TGraphErrors()

    g_SR_vs_Vov[sipm] = {}
    g_bestTh_vs_Vov[sipm] = {}

    g_SR_vs_bar[sipm] = {}
    g_bestTh_vs_bar[sipm] = {}
    g_Noise_vs_bar[sipm] = {}
    g_Stoch_vs_bar[sipm] = {}
    g_DCR_vs_bar[sipm] = {}

    fPS[sipm] = {}
    for ov in Vovs[sipm]:
        if (sipm == 'HPK_2E14_T-40C'): fPS[sipm][ov] = ROOT.TFile.Open('../plots/pulseShape_HPK_2E14_52deg_T-40C_Vov%.2f_ith1.root'%ov)
        if (sipm == 'FBK_2E14_T-40C'): fPS[sipm][ov] = ROOT.TFile.Open('../plots/pulseShape_FBK_2E14_52deg_T-40C_Vov%.2f_ith1.root'%ov)
        if (sipm == 'FBK_2E14_T-32C'): fPS[sipm][ov] = ROOT.TFile.Open('../plots/pulseShape_FBK_2E14_52deg_T-32C_Vov%.2f_ith1.root'%ov)
        if (sipm == 'FBK_2E14_T-22C'): fPS[sipm][ov] = ROOT.TFile.Open('../plots/pulseShape_FBK_2E14_52deg_T-22C_Vov%.2f_ith1.root'%ov)
        if (sipm == 'HPK_1E13_T0C')  : fPS[sipm][ov] = ROOT.TFile.Open('../plots/pulseShape_HPK_1E13_52deg_T0C_Vov%.2f_ith1.root'%ov)
        if (sipm == 'FBK_1E13_T0C')  : fPS[sipm][ov] = ROOT.TFile.Open('../plots/pulseShape_FBK_1E13_52deg_T0C_Vov%.2f_ith1.root'%ov)
        if (sipm == 'HPK_1E13_LYSOtype1_T0C'): fPS[sipm][ov] = ROOT.TFile.Open('../plots/pulseShape_HPK_1E13_LYSOtype1_58deg_T0C_Vov%.2f_ith1.root'%ov)
        if (sipm == 'HPK_1E13_LYSOtype1_T-6C'): fPS[sipm][ov] = ROOT.TFile.Open('../plots/pulseShape_HPK_1E13_LYSOtype1_58deg_T-6C_Vov%.2f_ith1.root'%ov)
        if (sipm == 'HPK_1E13_LYSOtype1_T-20C'): fPS[sipm][ov] = ROOT.TFile.Open('../plots/pulseShape_HPK_1E13_LYSOtype1_58deg_T-20C_Vov%.2f_ith1.root'%ov)
        if (sipm == 'HPK_1E13_LYSOtype1_T-40C'): fPS[sipm][ov] = ROOT.TFile.Open('../plots/pulseShape_HPK_1E13_LYSOtype1_58deg_T-40C_Vov%.2f_ith1.root'%ov)

        g_SR_vs_bar[sipm][ov] = ROOT.TGraphErrors()
        g_bestTh_vs_bar[sipm][ov] = ROOT.TGraphErrors()
        g_Noise_vs_bar[sipm][ov] = ROOT.TGraphErrors()
        g_Stoch_vs_bar[sipm][ov] = ROOT.TGraphErrors()
        g_DCR_vs_bar[sipm][ov] = ROOT.TGraphErrors()
        g_Tot_vs_SR[sipm][ov] = ROOT.TGraphErrors()

        
    f[sipm] = ROOT.TFile.Open(fnames[sipm])

    for bar in bars[sipm]:
        g[sipm][bar] = f[sipm].Get('g_deltaT_energyRatioCorr_bestTh_vs_vov_bar%02d_enBin01;1'%bar)
        if (g[sipm][bar].GetN()==0): 
            print 'No data for bar ', bar
            continue
        g_Noise_vs_Vov[sipm][bar] = ROOT.TGraphErrors()
        g_Stoch_vs_Vov[sipm][bar] = ROOT.TGraphErrors()
        g_DCR_vs_Vov[sipm][bar] = ROOT.TGraphErrors()
        g_Tot_vs_Vov[sipm][bar] = ROOT.TGraphErrors()

        g_Stoch_vs_Npe[sipm][bar] = ROOT.TGraphErrors()
        g_DCR_vs_Npe[sipm][bar] = ROOT.TGraphErrors()
        g_DCR_vs_DCRNpeSR[sipm][bar] = ROOT.TGraphErrors()

        g_SR_vs_Vov[sipm][bar] = ROOT.TGraphErrors()
        g_bestTh_vs_Vov[sipm][bar] = ROOT.TGraphErrors()
                   
        for ov in Vovs[sipm]:
            #if ( VovsEff[sipm][ov] < g[sipm][bar].GetX()[0] or VovsEff[sipm][ov] > g[sipm][bar].GetX()[g[sipm][bar].GetN()-1]): continue
            # get measured time resolution
            sigma_meas = g[sipm][bar].Eval(VovsEff[sipm][ov])
            err_sigma_meas = 0.
            for i in range(0,g[sipm][bar].GetN()):
                if g[sipm][bar].GetX()[i] == g[sipm][bar].Eval(VovsEff[sipm][ov]):
                    err_sigma_meas = g[sipm][bar].GetErrorY[i]
            # Npe and Gain at this OVeff
            irr = '2E14'
            if ('1E13' in sipm): irr = '1E13'
            Npe[sipm][ov]  = LO[sipm]*4.2*PDE(VovsEff[sipm][ov],sipm,irr)/PDE(3.50,sipm,'0') #LO is referred to 3.50 V OV
            gain = Gain(ov, sipm, irr)
            if ('type1' in sipm):
                Npe[sipm][ov]  = 0.86*3.75/math.cos(58*math.pi/180)*LO[sipm]*PDE(VovsEff[sipm][ov],sipm,'1E13')/PDE(3.50,sipm,'0') # type1 SiPM 
                gain = Gain(ov, sipm, '1E13')
            # get pulse shapes
            g_psL = fPS[sipm][ov].Get('g_pulseShapeL_bar%02d_Vov%.2f'%(bar,ov))
            g_psR = fPS[sipm][ov].Get('g_pulseShapeR_bar%02d_Vov%.2f'%(bar,ov))
            if (g_psL!=None): g_psL.SetName('g_pulseShapeL_bar%02d_Vov%.2f_%s'%(bar,ov,sipm))
            if (g_psR!=None): g_psR.SetName('g_pulseShapeR_bar%02d_Vov%.2f_%s'%(bar,ov,sipm))
            timingThreshold = findTimingThreshold(f[sipm].Get('g_deltaT_energyRatioCorr_vs_th_bar%02d_Vov%.2f_enBin01'%(bar,ov)))
            srL = -1
            srR = -1
            sr = -1
            err_srL = -1
            err_srR = -1
            c = ROOT.TCanvas('c_%s'%(g_psL.GetName().replace('g_pulseShapeL','pulseShape').replace('Vov%.2f'%ov,'VovEff%.2f'%VovsEff[sipm][ov])),'',600,600)  
            hdummy = ROOT.TH2F('hdummy','', 100, min(g_psR.GetX())-1., 30., 100, 0., 15.)
            #hdummy = ROOT.TH2F('hdummy','', 100, min(g_psL.GetX())-1., min(g_psL.GetX())+3, 100, 0., 15.)
            hdummy.GetXaxis().SetTitle('time [ns]')
            hdummy.GetYaxis().SetTitle('amplitude [#muA]')
            hdummy.Draw()
            gtempL = ROOT.TGraphErrors()
            gtempR = ROOT.TGraphErrors()
            if (g_psL!=None): srL,err_srL = getSlewRateFromPulseShape(g_psL, timingThreshold, np, gtempL, c)
            if (g_psR!=None): srR,err_srR = getSlewRateFromPulseShape(g_psR, timingThreshold, np, gtempR, c) 
            line = ROOT.TLine(min(g_psL.GetX())-1., timingThreshold*0.313, 30., timingThreshold*0.313)
            line.SetLineStyle(7)
            line.SetLineWidth(2)
            line.SetLineColor(ROOT.kOrange+1)        
            line.Draw('same')
            c.SaveAs(outdir+'/plotsSR/'+c.GetName()+'.png')   
            hdummy.Delete()
            if (srL>0 and srR>0):
                # weighted average
                sr =  ( (srL/(err_srL*err_srL) + srR/(err_srR*err_srR) ) / (1./(err_srL*err_srL) + 1./(err_srR*err_srR) ) )
                errSR = 1./math.sqrt( 1./(err_srL*err_srL)  +  1./(err_srR*err_srR) )
                errSR = errSR
            if (srL>0 and srR<0):
                sr = srL
                errSR = err_srL
            if (srL<0 and srR>0):
                sr = srR
                errSR = err_srR
            if (srL<0 and srR<0): continue
            errSR = math.sqrt(errSR*errSR+errSRsyst*errSRsyst*sr*sr) 
            ovEff = VovsEff[sipm][ov]
            #print sipm, ov, ovEff, gain, Npe[sipm][ov], srL, srR, sr, errSR
            g_SR_vs_Vov[sipm][bar].SetPoint( g_SR_vs_Vov[sipm][bar].GetN(), ovEff, sr )
            g_SR_vs_Vov[sipm][bar].SetPointError( g_SR_vs_Vov[sipm][bar].GetN()-1, 0, errSR )
            
            g_bestTh_vs_Vov[sipm][bar].SetPoint( g_bestTh_vs_Vov[sipm][bar].GetN(), ovEff, timingThreshold )
            g_bestTh_vs_Vov[sipm][bar].SetPointError( g_bestTh_vs_Vov[sipm][bar].GetN()-1, 0, 0 )
            
            g_SR_vs_bar[sipm][ov].SetPoint( g_SR_vs_bar[sipm][ov].GetN(), bar, sr )
            g_SR_vs_bar[sipm][ov].SetPointError( g_SR_vs_bar[sipm][ov].GetN()-1, 0, errSR )
            
            g_bestTh_vs_bar[sipm][ov].SetPoint( g_bestTh_vs_bar[sipm][ov].GetN(), bar, timingThreshold )
            g_bestTh_vs_bar[sipm][ov].SetPointError( g_bestTh_vs_bar[sipm][ov].GetN()-1, 0, 0)
            
            g_Noise_vs_bar[sipm][ov].SetPoint( g_Noise_vs_bar[sipm][ov].GetN(), bar, sigma_noise(sr) )
            g_Noise_vs_bar[sipm][ov].SetPointError( g_Noise_vs_bar[sipm][ov].GetN()-1, 0,  0.5*(sigma_noise(sr*(1-errSR/sr))-sigma_noise(sr*(1+errSR/sr))) )
            
            g_Noise_vs_Vov[sipm][bar].SetPoint(g_Noise_vs_Vov[sipm][bar].GetN(), ovEff, sigma_noise(sr))
            g_Noise_vs_Vov[sipm][bar].SetPointError(g_Noise_vs_Vov[sipm][bar].GetN()-1, 0, 0.5*(sigma_noise(sr*(1-errSR/sr))-sigma_noise(sr*(1+errSR/sr))))
            
            g_Tot_vs_SR[sipm][ov].SetPoint(g_Tot_vs_SR[sipm][ov].GetN(), sr, sigma_meas)
            g_Tot_vs_SR[sipm][ov].SetPointError(g_Tot_vs_SR[sipm][ov].GetN()-1, errSR, err_sigma_meas)
            
            # compute s_stoch by scaling the stochastic term measured for non-irradiated SiPMs (40 ps HPK, 45 ps FBK) for sqrt(PDE) 
            sigma_stoch = sigma_stoch_ref[sipm]/pow( PDE(ovEff,sipm,'2E14')/PDE(ov_ref,sipm,'0'), 0.70  )
            err_sigma_stoch = err_sigma_stoch_ref/math.sqrt( PDE(ovEff,sipm,'2E14')/PDE(ov_ref,sipm,'0') )
            if ('1E13' in sipm):
                sigma_stoch = sigma_stoch_ref[sipm]/pow( PDE(ovEff,sipm,'1E13')/PDE(ov_ref,sipm,'0'), 0.70  )
                err_sigma_stoch = err_sigma_stoch_ref/math.sqrt( PDE(ovEff,sipm,'1E13')/PDE(ov_ref,sipm,'0') )
            g_Stoch_vs_Vov[sipm][bar].SetPoint(g_Stoch_vs_Vov[sipm][bar].GetN(), ovEff, sigma_stoch)
            g_Stoch_vs_Vov[sipm][bar].SetPointError(g_Stoch_vs_Vov[sipm][bar].GetN()-1, 0, err_sigma_stoch)
            g_Stoch_vs_bar[sipm][ov].SetPoint( g_Stoch_vs_bar[sipm][ov].GetN(), bar, sigma_stoch )
            g_Stoch_vs_bar[sipm][ov].SetPointError( g_Stoch_vs_bar[sipm][ov].GetN()-1, 0,  err_sigma_stoch)
            # compute sigma_DCR as difference in quadrature between measure tRes and noise, stoch
            if ( sigma_meas*sigma_meas - sigma_stoch*sigma_stoch - sigma_noise(sr)*sigma_noise(sr) > 0):
                sigma_dcr = math.sqrt( sigma_meas*sigma_meas - sigma_stoch*sigma_stoch - sigma_noise(sr)*sigma_noise(sr) ) 
                err_sigma_dcr = 1./sigma_dcr * math.sqrt( pow( err_sigma_meas*sigma_meas,2) + pow(sigma_noise(sr)*g_Noise_vs_Vov[sipm][bar].GetErrorY(g_Noise_vs_Vov[sipm][bar].GetN()-1),2))
                g_DCR_vs_Vov[sipm][bar].SetPoint(g_DCR_vs_Vov[sipm][bar].GetN(), ovEff, sigma_dcr)
                g_DCR_vs_Vov[sipm][bar].SetPointError(g_DCR_vs_Vov[sipm][bar].GetN()-1, 0, err_sigma_dcr)
                g_DCR_vs_bar[sipm][ov].SetPoint( g_DCR_vs_bar[sipm][ov].GetN(), bar, sigma_dcr )
                g_DCR_vs_bar[sipm][ov].SetPointError( g_DCR_vs_bar[sipm][ov].GetN()-1, 0,  err_sigma_dcr)
                
                dcr = DCR[sipm][ov]
                if ('2E14' in sipm): dcr = 0.7*DCR[sipm][ov]
                
                g_DCR_vs_Npe[sipm][bar].SetPoint( g_DCR_vs_Npe[sipm][bar].GetN(), math.sqrt(dcr)/Npe[sipm][ov]/(math.sqrt(30.)/3000.), sigma_dcr )
                g_DCR_vs_Npe[sipm][bar].SetPointError( g_DCR_vs_Npe[sipm][bar].GetN()-1, 0,  err_sigma_dcr)

                beta = 0.5
                x = pow( dcr, 0.5 ) / Npe[sipm][ov] / pow( sr/Npe[sipm][ov], beta)
                x = x / (pow( 30., 0.5 ) / 3000. / pow( 5./3000., beta ))
                g_DCR_vs_DCRNpeSR[sipm][bar].SetPoint( g_DCR_vs_DCRNpeSR[sipm][bar].GetN(), x, sigma_dcr )
                g_DCR_vs_DCRNpeSR[sipm][bar].SetPointError( g_DCR_vs_DCRNpeSR[sipm][bar].GetN()-1, 0,  err_sigma_dcr)

                g_DCR_vs_DCRNpeSR_all[sipm].SetPoint( g_DCR_vs_DCRNpeSR_all[sipm].GetN(), x, sigma_dcr )
                g_DCR_vs_DCRNpeSR_all[sipm].SetPointError( g_DCR_vs_DCRNpeSR_all[sipm].GetN()-1, 0,  err_sigma_dcr)

            #
            sigma_tot = math.sqrt( sigma_stoch*sigma_stoch + sigma_noise(sr)*sigma_noise(sr) )
            err_sigma_tot = 1./sigma_tot * math.sqrt( pow( err_sigma_stoch*sigma_stoch,2) + pow(sigma_noise(sr)*g_Noise_vs_Vov[sipm][bar].GetErrorY(g_Noise_vs_Vov[sipm][bar].GetN()-1),2))

            g_Tot_vs_Vov[sipm][bar].SetPoint(g_Tot_vs_Vov[sipm][bar].GetN(), ovEff, sigma_tot)
            g_Tot_vs_Vov[sipm][bar].SetPointError(g_Tot_vs_Vov[sipm][bar].GetN()-1, 0, err_sigma_tot)

            #print sipm,' OV = %.2f  gain = %d  Npe = %d  bar = %02d  thr = %02d  SR = %.1f   noise = %.1f    stoch = %.1f   tot = %.1f'%(ov, gain, Npe, bar, timingThreshold, sr, sigma_noise(sr), sigma_stoch, sigma_tot)
            

g_DCR_vs_DCR = {}
g1_DCR_vs_Vov = {}
g_DCR_vs_DCRNpe_average = {}

for sipm in sipmTypes:
    g1_DCR_vs_Vov[sipm] = ROOT.TGraphErrors() 
    g_DCR_vs_DCRNpe_average[sipm] = ROOT.TGraphErrors()
    for ov in Vovs[sipm]:
        dcr = DCR[sipm][ov]
        if ('2E14' in sipm): dcr = 0.7*DCR[sipm][ov]
        g1_DCR_vs_Vov[sipm].SetPoint(g1_DCR_vs_Vov[sipm].GetN(), VovsEff[sipm][ov], DCR[sipm][ov])

        x = math.sqrt(DCR[sipm][ov])/Npe[sipm][ov]/ (math.sqrt(30.)/3000)
        if (g_DCR_vs_bar[sipm][ov].GetN()==0):continue
        g_DCR_vs_DCRNpe_average[sipm].SetPoint( g_DCR_vs_DCRNpe_average[sipm].GetN(), x,  g_DCR_vs_bar[sipm][ov].GetMean(2))
        g_DCR_vs_DCRNpe_average[sipm].SetPointError( g_DCR_vs_DCRNpe_average[sipm].GetN()-1, 0,  g_DCR_vs_bar[sipm][ov].GetRMS(2))


for bar in range(0,16):
    g_DCR_vs_DCR[bar] = ROOT.TGraphErrors()
    for sipm in sipmTypes:
        if (bar not in g_DCR_vs_Vov[sipm].keys()): continue
        ovEff = 1.60
        g_DCR_vs_DCR[bar].SetPoint( g_DCR_vs_DCR[bar].GetN(), g1_DCR_vs_Vov[sipm].Eval(ovEff), g_DCR_vs_Vov[sipm][bar].Eval(ovEff))
        g_DCR_vs_DCR[bar].SetPointError( g_DCR_vs_DCR[bar].GetN()-1, 0, g_DCR_vs_Vov[sipm][bar].GetErrorY(0))


# draw
c1 = {}
c2 = {}
c3 = {}
c4 = {}
hdummy1 = {}
hdummy2 = {}
hdummy3 = {}
hdummy4 = {}
leg = {}

for sipm in sipmTypes:
    c1[sipm] = {}
    hdummy1[sipm] = {}
    leg[sipm] = ROOT.TLegend(0.55,0.70,0.89,0.89)
    leg[sipm].SetBorderSize(0)
    leg[sipm].SetFillStyle(0)
    for i,bar in enumerate(bars[sipm]):
        if (bar not in g[sipm].keys()): continue
        if (g[sipm][bar].GetN()==0): continue
        c1[sipm][bar] =  ROOT.TCanvas('c_timeResolution_vs_Vov_%s_bar%02d'%(sipm,bar),'c_timeResolution_vs_Vov_%s_bar%02d'%(sipm,bar),600,600)
        c1[sipm][bar].SetGridy()
        c1[sipm][bar].cd()
        xmin = g[sipm][bar].GetX()[0]-0.2
        xmax = g[sipm][bar].GetX()[g[sipm][bar].GetN()-1]+0.2
        #xmin = 1.1
        #xmax = 2.5
        #if ('1E13' in sipm):
        #    xmin = 1.1
        #    xmax = 2.5
        #if ('1E13' in sipm and 'type1' in sipm):
        #    xmin = 1.1
        #    xmax = 5.0
        hdummy1[sipm][bar] = ROOT.TH2F('hdummy1_%s_%d'%(sipm,bar),'',100,xmin,xmax,180,0,180)
        hdummy1[sipm][bar].GetXaxis().SetTitle('V_{OV}^{eff} [V]')
        hdummy1[sipm][bar].GetYaxis().SetTitle('#sigma_{t} [ps]')
        hdummy1[sipm][bar].Draw()
        g[sipm][bar].SetMarkerStyle(20)
        g[sipm][bar].SetMarkerSize(1)
        g[sipm][bar].SetMarkerColor(1)
        g[sipm][bar].SetLineColor(1)
        g[sipm][bar].SetLineWidth(2)
        g[sipm][bar].Draw('plsame')
        if (bar not in g_Noise_vs_Vov[sipm].keys()): continue
        g_Noise_vs_Vov[sipm][bar].SetLineWidth(2)
        g_Noise_vs_Vov[sipm][bar].SetLineColor(ROOT.kBlue)
        g_Noise_vs_Vov[sipm][bar].SetFillColor(ROOT.kBlue)
        g_Noise_vs_Vov[sipm][bar].SetFillColorAlpha(ROOT.kBlue,0.2)
        g_Noise_vs_Vov[sipm][bar].SetFillStyle(3004)
        g_Noise_vs_Vov[sipm][bar].Draw('E3lsame')
        g_Stoch_vs_Vov[sipm][bar].SetLineWidth(2)
        g_Stoch_vs_Vov[sipm][bar].SetLineColor(ROOT.kGreen+2)
        g_Stoch_vs_Vov[sipm][bar].SetFillColor(ROOT.kGreen+2)
        g_Stoch_vs_Vov[sipm][bar].SetFillStyle(3001)
        g_Stoch_vs_Vov[sipm][bar].SetFillColorAlpha(ROOT.kGreen+2,0.2)
        g_Stoch_vs_Vov[sipm][bar].Draw('E3lsame')
        g_DCR_vs_Vov[sipm][bar].SetLineWidth(2)
        g_DCR_vs_Vov[sipm][bar].SetLineColor(ROOT.kOrange+2)
        g_DCR_vs_Vov[sipm][bar].SetFillColor(ROOT.kOrange+2)
        g_DCR_vs_Vov[sipm][bar].SetFillColorAlpha(ROOT.kOrange+2,0.2)
        g_DCR_vs_Vov[sipm][bar].SetFillStyle(3001)
        g_DCR_vs_Vov[sipm][bar].Draw('E3lsame')
        g_Tot_vs_Vov[sipm][bar].SetLineWidth(2)
        g_Tot_vs_Vov[sipm][bar].SetLineColor(ROOT.kRed+1)
        g_Tot_vs_Vov[sipm][bar].SetFillColor(ROOT.kRed+1)
        g_Tot_vs_Vov[sipm][bar].SetFillColorAlpha(ROOT.kRed+1,0.2)
        g_Tot_vs_Vov[sipm][bar].SetFillStyle(3001)
        #g_Tot_vs_Vov[sipm][bar].Draw('E3lsame')
        if (i==0):
            leg[sipm].AddEntry(g[sipm][bar], 'data', 'PL')
            leg[sipm].AddEntry(g_Noise_vs_Vov[sipm][bar], 'noise', 'PL')
            leg[sipm].AddEntry(g_Stoch_vs_Vov[sipm][bar], 'stoch', 'PL')
            leg[sipm].AddEntry(g_DCR_vs_Vov[sipm][bar], 'DCR', 'PL')
            #leg[sipm].AddEntry(g_Tot_vs_Vov[sipm][bar], 'stoch (+) noise', 'PL')
        leg[sipm].Draw('same')
        c1[sipm][bar].SaveAs(outdir+'/'+c1[sipm][bar].GetName()+'.png')
        hdummy1[sipm][bar].Delete()
        c1[sipm][bar].Delete()


# total time resolution vs SR
for sipm in sipmTypes:
    for ov in Vovs[sipm]:
        if (ov not in g_Tot_vs_SR[sipm].keys()): continue
        c =  ROOT.TCanvas('c_timeResolution_vs_SR_%s_Vov%.02f'%(sipm,VovsEff[sipm][ov]),'c_timeResolution_vs_SR_%s_Vov%.02f'%(sipm,VovsEff[sipm][ov]),600,600)
        c.SetGridy()
        c.cd()
        xmin = 0.
        xmax = g_Tot_vs_SR[sipm][ov].GetMean() + 5.
        ymin = g_Tot_vs_SR[sipm][ov].GetMean(2)-40.
        ymax = g_Tot_vs_SR[sipm][ov].GetMean(2)+40.
        hdummy = ROOT.TH2F('hdummy_%s_%d'%(sipm,ov),'',100, xmin, xmax, 100, ymin, ymax)
        hdummy.GetXaxis().SetTitle('slew rate at the timing thr. [#muA/ns]')
        hdummy.GetYaxis().SetTitle('#sigma_{t} [ps]')
        hdummy.Draw()
        g_Tot_vs_SR[sipm][ov].SetMarkerStyle(20)
        g_Tot_vs_SR[sipm][ov].SetMarkerSize(1)
        g_Tot_vs_SR[sipm][ov].SetMarkerColor(1)
        g_Tot_vs_SR[sipm][ov].SetLineColor(1)
        g_Tot_vs_SR[sipm][ov].SetLineWidth(2)
        g_Tot_vs_SR[sipm][ov].Draw('psame')
        c.SaveAs(outdir+'/'+c.GetName()+'.png')
        hdummy.Delete()
        c.Delete()



markers = { 'HPK_2E14_T-40C' : 20 ,
            'FBK_2E14_T-40C' : 21 ,
            'FBK_2E14_T-32C' : 21 ,
            'FBK_2E14_T-22C' : 21 ,
            'HPK_1E13_T0C' : 20,
            'FBK_1E13_T0C' : 21,
            'HPK_1E13_LYSOtype1_T0C' : 20,
            'HPK_1E13_LYSOtype1_T-6C' : 20,
            'HPK_1E13_LYSOtype1_T-20C' : 20,
            'HPK_1E13_LYSOtype1_T-40C' : 20}


cols = { 'HPK_2E14_T-40C' : ROOT.kBlack ,
         'FBK_2E14_T-40C' : ROOT.kRed ,
         'FBK_2E14_T-32C' : ROOT.kOrange+1 ,
         'FBK_2E14_T-22C' : ROOT.kMagenta,
         'HPK_1E13_T0C'   : ROOT.kBlack,
         'FBK_1E13_T0C'   : ROOT.kRed,
         'HPK_1E13_LYSOtype1_T0C'   : ROOT.kGreen,
         'HPK_1E13_LYSOtype1_T-6C'  : ROOT.kOrange+1,
         'HPK_1E13_LYSOtype1_T-20C' : ROOT.kMagenta,
         'HPK_1E13_LYSOtype1_T-40C' : ROOT.kRed}

# vs npe
for bar in range(0,16):      
    c2 =  ROOT.TCanvas('c_timeResolutionDCR_vs_Npe_bar%02d'%(bar),'c_timeResolutionDCR_vs_Npe_bar%02d'%(bar),600,600)
    c2.SetGridx()
    c2.SetGridy()
    c2.cd()
    xmax = 2
    ymax = 180
    #if ('1E13' in sipm):        
    #    xmax = 2
    #    ymax = 100
    hdummy2 = ROOT.TH2F('hdummy2_%d'%(bar),'',100,0,xmax,100,0,ymax)
    hdummy2.GetXaxis().SetTitle('#sqrt{DCR}/Npe')
    hdummy2.GetYaxis().SetTitle('#sigma_{t}^{DCR} [ps]')
    hdummy2.Draw()
    for sipm in sipmTypes:    
        if (bar not in g_DCR_vs_Npe[sipm].keys()): continue
        g_DCR_vs_Npe[sipm][bar].SetMarkerStyle(markers[sipm])
        g_DCR_vs_Npe[sipm][bar].SetMarkerColor(cols[sipm])
        g_DCR_vs_Npe[sipm][bar].SetLineWidth(1)
        g_DCR_vs_Npe[sipm][bar].SetLineColor(cols[sipm])
        g_DCR_vs_Npe[sipm][bar].Draw('psame')
        #fitFun = ROOT.TF1('fitFun_%s_%2d'%(sipm,bar),'pol1',g_DCR_vs_Npe[sipm][bar].GetX()[0]-0.05, g_DCR_vs_Npe[sipm][bar].GetX()[g_DCR_vs_Npe[sipm][bar].GetN()-1]+0.05)
        #fitFun.SetLineColor(cols[sipm])
        #fitFun.SetLineStyle(2)
        #g_DCR_vs_Npe[sipm][bar].Fit(fitFun,'QRS')
    c2.SaveAs(outdir+'/'+c2.GetName()+'.png')
    hdummy2.Delete()
    c2.Delete()

    c2 =  ROOT.TCanvas('c_timeResolutionDCR_vs_DCRNpeSR_bar%02d'%(bar),'c_timeResolutionDCR_vs_DCRNpeSR_bar%02d'%(bar),600,600)
    c2.SetGridx()
    c2.SetGridy()
    c2.cd()    
    hdummy2 = ROOT.TH2F('hdummy2_%d'%(bar),'',100,0,2,100,0,ymax)
    hdummy2.GetXaxis().SetTitle('DCR^{#alpha}/Npe/(SR/Npe)^{#beta}')
    hdummy2.GetYaxis().SetTitle('#sigma_{t}^{DCR} [ps]')
    hdummy2.Draw()
    for sipm in sipmTypes:
        if (bar not in g_DCR_vs_DCRNpeSR[sipm].keys()): continue
        g_DCR_vs_DCRNpeSR[sipm][bar].SetMarkerStyle(markers[sipm])
        g_DCR_vs_DCRNpeSR[sipm][bar].SetMarkerColor(cols[sipm])
        g_DCR_vs_DCRNpeSR[sipm][bar].SetLineWidth(1)
        g_DCR_vs_DCRNpeSR[sipm][bar].SetLineColor(cols[sipm])
        g_DCR_vs_DCRNpeSR[sipm][bar].Draw('psame')
        #fitFun = ROOT.TF1('fitFun_%s_%2d'%(sipm,bar),'pol1',g_DCR_vs_DCRNpeSR[sipm][bar].GetX()[0]-0.05, g_DCR_vs_DCRNpeSR[sipm][bar].GetX()[g_DCR_vs_DCRNpeSR[sipm][bar].GetN()-1]+0.05)
        #fitFun.SetLineColor(cols[sipm])
        #fitFun.SetLineStyle(2)
        #g_DCR_vs_DCRNpeSR[sipm][bar].Fit(fitFun,'QRS')
    c2.SaveAs(outdir+'/'+c2.GetName()+'.png')
    hdummy2.Delete()
    c2.Delete()


c2 =  ROOT.TCanvas('c_timeResolutionDCR_vs_DCRNpe_average','c_timeResolutionDCR_vs_DCRNpe_average',600,600)
c2.SetGridx()
c2.SetGridy()
c2.cd()    
hdummy2 = ROOT.TH2F('hdummy2_%d'%(bar),'',100,0,2.0,100,0,ymax)
hdummy2.GetXaxis().SetTitle('#sqrt{DCR}/Npe')
hdummy2.GetYaxis().SetTitle('#sigma_{t}^{DCR} [ps]')
hdummy2.Draw()
for sipm in sipmTypes:
    g_DCR_vs_DCRNpe_average[sipm].SetMarkerStyle(markers[sipm])
    g_DCR_vs_DCRNpe_average[sipm].SetMarkerColor(cols[sipm])
    g_DCR_vs_DCRNpe_average[sipm].SetLineWidth(1)
    g_DCR_vs_DCRNpe_average[sipm].SetLineColor(cols[sipm])
    g_DCR_vs_DCRNpe_average[sipm].Draw('psame')
    outfile.cd()
    g_DCR_vs_DCRNpe_average[sipm].Write('g_DCR_vs_DCRNpe_average_%s'%sipm)
c2.SaveAs(outdir+'/'+c2.GetName()+'.png')
hdummy2.Delete()
c2.Delete()


c2 =  ROOT.TCanvas('c_timeResolutionDCR_vs_DCRNpeSR_allBars','c_timeResolutionDCR_vs_DCRNpeSR_allBars',600,600)
c2.SetGridx()
c2.SetGridy()
c2.cd()    
hdummy2 = ROOT.TH2F('hdummy2_%d'%(bar),'',100,0,2,100,0,ymax)
hdummy2.GetXaxis().SetTitle('DCR^{#alpha}/Npe/(SR/Npe)^{#beta}')
hdummy2.GetYaxis().SetTitle('#sigma_{t}^{DCR} [ps]')
hdummy2.Draw()
for sipm in sipmTypes:
    g_DCR_vs_DCRNpeSR_all[sipm].SetMarkerStyle(markers[sipm])
    g_DCR_vs_DCRNpeSR_all[sipm].SetMarkerColor(cols[sipm])
    g_DCR_vs_DCRNpeSR_all[sipm].SetLineWidth(1)
    g_DCR_vs_DCRNpeSR_all[sipm].SetLineColor(cols[sipm])
    g_DCR_vs_DCRNpeSR_all[sipm].Draw('psame')
c2.SaveAs(outdir+'/'+c2.GetName()+'.png')
hdummy2.Delete()
c2.Delete()



# sigma_DCR vs DCR @ reference OV
for bar in range(0,16):      
    c2 =  ROOT.TCanvas('c_timeResolutionDCR_vs_DCR_bar%02d'%(bar),'c_timeResolutionDCR_vs_DCR_bar%02d'%(bar),600,600)
    c2.SetGridx()
    c2.SetGridy()
    c2.cd()
    hdummy2 = ROOT.TH2F('hdummy2_%d'%(bar),'',100,0,100,100,0,180)
    hdummy2.GetXaxis().SetTitle('DCR [GHz]')
    hdummy2.GetYaxis().SetTitle('#sigma_{t}^{DCR} [ps]')
    hdummy2.Draw()
    g_DCR_vs_DCR[bar].SetMarkerStyle(20)
    g_DCR_vs_DCR[bar].SetMarkerSize(1)
    g_DCR_vs_DCR[bar].Draw('plsame')
    fitFun = ROOT.TF1('fitFun_%2d'%(bar),'[0]*pow(x,[1])',0,200)
    fitFun.SetParameters(20.,0.5)
    g_DCR_vs_DCR[bar].Fit(fitFun,'QRS')
    c2.SaveAs(outdir+'/'+c2.GetName()+'.png')
    fitFun.Delete()
    hdummy2.Delete()
    c2.Delete()
    


# SR and best threshold vs Vov
leg2 = ROOT.TLegend(0.15,0.70,0.45,0.89)
leg2.SetBorderSize(0)
leg2.SetFillStyle(0)
for i,bar in enumerate(bars[sipm]):
    if (bar not in g[sipm].keys()): continue
    c3[bar] = ROOT.TCanvas('c_slewRate_vs_Vov_bar%02d'%(bar),'c_slewRate_vs_Vov_bar%02d'%(bar),600,600)
    c3[bar].SetGridy()
    c3[bar].cd()
    xmin = 0.9
    xmax = 2.5
    ymax = 15 
    if ('1E13' in sipm):
        xmin = 1.1
        xmax = 5.0
        ymax = 35.
    hdummy3[bar] = ROOT.TH2F('hdummy3_%d'%(bar),'',100,xmin,xmax,100,0,ymax)
    hdummy3[bar].GetXaxis().SetTitle('V_{OV}^{eff} [V]')
    hdummy3[bar].GetYaxis().SetTitle('slew rate at the timing thr. [#muA/ns]')
    hdummy3[bar].Draw()
    for j,sipm in enumerate(sipmTypes):
        if (bar not in g_SR_vs_Vov[sipm].keys()): continue
        if (i==0):
            leg2.AddEntry(g_SR_vs_Vov[sipm][bar], '%s'%sipm, 'PL')
        g_SR_vs_Vov[sipm][bar].SetMarkerStyle(markers[sipm])
        g_SR_vs_Vov[sipm][bar].SetMarkerColor(cols[sipm])
        g_SR_vs_Vov[sipm][bar].SetLineColor(cols[sipm])
        g_SR_vs_Vov[sipm][bar].Draw('plsame')
    leg2.Draw()
    c3[bar].SaveAs(outdir+'/'+c3[bar].GetName()+'.png')
    
    c4[bar] = ROOT.TCanvas('c_bestTh_vs_Vov_bar%02d'%(bar),'c_bestTh_vs_Vov_bar%02d'%(bar),600,600)
    c4[bar].SetGridy()
    c4[bar].cd()
    hdummy4[bar] = ROOT.TH2F('hdummy4_%d'%(bar),'',100,xmin,xmax,100,0,20)
    hdummy4[bar].GetXaxis().SetTitle('V_{OV}^{eff} [V]')
    hdummy4[bar].GetYaxis().SetTitle('best threshold [DAC]')
    hdummy4[bar].Draw()
    for j,sipm in enumerate(sipmTypes):
        if (bar not in g_bestTh_vs_Vov[sipm].keys()): continue
        g_bestTh_vs_Vov[sipm][bar].SetMarkerStyle(markers[sipm])
        g_bestTh_vs_Vov[sipm][bar].SetMarkerColor(cols[sipm])
        g_bestTh_vs_Vov[sipm][bar].SetLineColor(cols[sipm])
        g_bestTh_vs_Vov[sipm][bar].Draw('plsame')
    c4[bar].SaveAs(outdir+'/'+c4[bar].GetName()+'.png')



    
# SR and best threshold vs bar
c5 = {}
hdummy5 = {}
c6 = {}
hdummy6 = {}
c7 = {}
hdummy7 = {}
c8 = {}
hdummy8 = {}
for j,sipm in enumerate(sipmTypes):
    for ov in Vovs[sipm]:
        if (ov not in g_SR_vs_bar[sipm].keys()): continue
        c5[ov] = ROOT.TCanvas('c_slewRate_vs_bar_%s_Vov%.2f'%(sipm,VovsEff[sipm][ov]),'c_slewRate_vs_bar_%s_Vov%.2f'%(sipm,VovsEff[sipm][ov]),600,600)
        c5[ov].SetGridy()
        c5[ov].cd()
        #hdummy5[ov] = ROOT.TH2F('hdummy5_%d'%(ov),'',100,-0.5,15.5,100,0,15)
        hdummy5[ov] = ROOT.TH2F('hdummy5_%d'%(ov),'',100,-0.5,15.5,100,0,35)
        hdummy5[ov].GetXaxis().SetTitle('bar')
        hdummy5[ov].GetYaxis().SetTitle('slew rate at the timing thr. [#muA/ns]')
        hdummy5[ov].Draw()
        g_SR_vs_bar[sipm][ov].SetMarkerStyle(markers[sipm])
        g_SR_vs_bar[sipm][ov].SetMarkerColor(cols[sipm])
        g_SR_vs_bar[sipm][ov].SetLineColor(cols[sipm])
        g_SR_vs_bar[sipm][ov].Draw('psame')
        leg2.Draw()
        c5[ov].SaveAs(outdir+'/'+c5[ov].GetName()+'.png')

        c6[ov] = ROOT.TCanvas('c_bestTh_vs_bar_%s_Vov%.2f'%(sipm,VovsEff[sipm][ov]),'c_bestTh_vs_bar_%s_Vov%.2f'%(sipm,VovsEff[sipm][ov]),600,600)
        c6[ov].SetGridy()
        c6[ov].cd()
        hdummy6[ov] = ROOT.TH2F('hdummy6_%d'%(ov),'',100,-0.5,15.5,100,0,20)
        hdummy6[ov].GetXaxis().SetTitle('bar')
        hdummy6[ov].GetYaxis().SetTitle('timing threshold [DAC]')
        hdummy6[ov].Draw()
        g_bestTh_vs_bar[sipm][ov].SetMarkerStyle(markers[sipm])
        g_bestTh_vs_bar[sipm][ov].SetMarkerColor(cols[sipm])
        g_bestTh_vs_bar[sipm][ov].SetLineColor(cols[sipm])
        g_bestTh_vs_bar[sipm][ov].Draw('plsame')
        leg2.Draw()        
        c6[ov].SaveAs(outdir+'/'+c6[ov].GetName()+'.png')

        c7[ov] = ROOT.TCanvas('c_noise_vs_bar_%s_Vov%.2f'%(sipm,VovsEff[sipm][ov]),'c_noise_vs_bar_%s_Vov%.2f'%(sipm,VovsEff[sipm][ov]),600,600)
        c7[ov].SetGridy()
        c7[ov].cd()
        hdummy7[ov] = ROOT.TH2F('hdummy7_%d'%(ov),'',100,-0.5,15.5,100,0,80)
        hdummy7[ov].GetXaxis().SetTitle('bar')
        hdummy7[ov].GetYaxis().SetTitle('#sigma_{t, noise} [ps]')
        hdummy7[ov].Draw()
        g_Noise_vs_bar[sipm][ov].SetMarkerStyle( markers[sipm] )
        g_Noise_vs_bar[sipm][ov].SetMarkerColor(cols[sipm])
        g_Noise_vs_bar[sipm][ov].SetLineColor(cols[sipm])
        g_Noise_vs_bar[sipm][ov].Draw('psame')
        leg2.Draw()
        c7[ov].SaveAs(outdir+'/'+c7[ov].GetName()+'.png')

        c8[ov] = ROOT.TCanvas('c_stoch_vs_bar_%s_Vov%.2f'%(sipm,VovsEff[sipm][ov]),'c_stoch_vs_bar_%s_Vov%.2f'%(sipm,VovsEff[sipm][ov]),600,600)
        c8[ov].SetGridy()
        c8[ov].cd()
        hdummy8[ov] = ROOT.TH2F('hdummy8_%d'%(ov),'',100,-0.5,15.5,100,0,80)
        hdummy8[ov].GetXaxis().SetTitle('bar')
        hdummy8[ov].GetYaxis().SetTitle('#sigma_{t, stoch} [ps]')
        hdummy8[ov].Draw()
        g_Stoch_vs_bar[sipm][ov].SetMarkerStyle( markers[sipm] )
        g_Stoch_vs_bar[sipm][ov].SetMarkerColor(cols[sipm])
        g_Stoch_vs_bar[sipm][ov].SetLineColor(cols[sipm])
        g_Stoch_vs_bar[sipm][ov].Draw('psame')
        leg2.Draw()
        c8[ov].SaveAs(outdir+'/'+c8[ov].GetName()+'.png')

raw_input('OK?')
