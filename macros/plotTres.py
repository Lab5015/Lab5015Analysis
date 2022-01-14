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
ROOT.gStyle.SetOptFit(1111)


outdir = '/eos/user/m/malberti/www/MTD/TOFHIR2X/MTDTB_CERN_Oct21/timeResolution_vs_Vov_unirr/'
if (os.path.exists(outdir)==False):
    os.mkdir(outdir)
if (os.path.exists(outdir+'/plotsSR')==False):
    os.mkdir(outdir+'/plotsSR/')

def Gain(ov, sipm):
    if ('HPK' in sipm):
        return 36890. + 97602.*ov # HPK
    if ('FBK' in sipm):
        return 50739. + 95149.*ov # FBK

def PDE(ov, sipm):
        if ('HPK' in sipm):
            return 0.384 * ( 1. - math.exp(-1.*0.583*ov) )
        if ('FBK' in sipm):
            return 0.466 * ( 1. - math.exp(-1.*0.314*ov) )

def sigma_noise(sr):
    noise_single = math.sqrt( pow(420./sr,2) + 16.7*16.7)
    return noise_single / math.sqrt(2)


def getSlewRateFromPulseShape(g1, timingThreshold, npoints, canvas=None):
    if ( g1.GetN()/2 <= npoints): return (-1, -1)
    # find index at the timing threshold
    itiming = 0
    for i in range(0,g1.GetN()):
        #print g1.GetY()[i], round(g1.GetY()[i]/0.313)
        if (round(g1.GetY()[i]/0.313) == timingThreshold):
            itiming = i
            break
    #fitSR = ROOT.TF1('fitSR', 'expo', 0, 10)
    fitSR = ROOT.TF1('fitSR', 'pol1', 0, 10)
    fitSR.SetLineColor(g1.GetMarkerColor()+1)
    tmin = min(g1.GetX())
    tmax = g1.GetX()[npoints]
    fitSR.SetRange(tmin,tmax)
    fitSR.SetParameters(0, 10)
    fitStatus = int(g1.Fit(fitSR, 'QRS+'))
    #fitStatus = int(g1.Fit(fitSR, 'QRN'))
    #redefine errors to get chi2/NDF = 1
    #for i in range(0,g1.GetN()):
    #    if (fitSR.GetNDF()>0): 
    #        err = g1.GetErrorX(i)
    #        g1.SetPointError(i, err*math.sqrt( fitSR.GetChisquare()/fitSR.GetNDF() ) , 0)
    #fitStatus = int(g1.Fit(fitSR, 'QRS+'))
    sr = fitSR.Derivative( g1.GetX()[itiming])
    err_sr = fitSR.GetParError(1)
    #print g1.GetName(), 'SR = ', sr,'+/-', fitSR.GetParError(1) , '   chi2 = ', fitSR.GetChisquare(), '   NDF = ', fitSR.GetNDF(), '  fitStatus = ', fitStatus
    #if (fitSR.GetNDF()==0 or (fitSR.GetNDF()>0 and fitSR.GetChisquare()/fitSR.GetNDF()>10)):
    #if (fitSR.GetNDF()==0):
    #    sr = -1
    #    err_sr = -1
    tmin = g1.GetX()[0]-1.0
    tmax = g1.GetX()[g1.GetN()-1]+1.0
    if (canvas!=None):
        canvas.cd()
        g1.Draw('psames')
        canvas.Update()
        ps = g1.FindObject("stats")
        ps.SetTextColor(g1.GetMarkerColor())
        if ('R' in g1.GetName()):
            ps.SetY1NDC(0.5)# new y start position
            ps.SetY2NDC(0.7)#new y end position

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

#sipmTypes = ['HPK_unirr_LYSO528','FBK_unirr_LYSO422']
sipmTypes = ['HPK_unirr_LYSO528','FBK_unirr_LYSO422','HPK_unirr_LYSOwithSlit']
#sipmTypes = ['HPK_unirr_LYSOwithSlit','FBK_unirr_LYSO422']
#sipmTypes = ['HPK_unirr_LYSO528','HPK_unirr_LYSOwithSlit']

fnames = {'HPK_unirr_LYSO528'      : '../plots/HPK528_unirr_52deg_T10C_summary.root',
          'FBK_unirr_LYSO422'      : '../plots/FBK_unirr_52deg_T10C_summary.root',
          'HPK_unirr_LYSOwithSlit' : '../plots/HPK_unirr_LYSOwithSlit_52deg_T18C_new_summary.root',
          }

LO = { 'HPK_unirr_LYSO528' : 1450.,
       'FBK_unirr_LYSO422' :1040.,
       'HPK_unirr_LYSOwithSlit' : 1100.}

#fitFunc_SR = (ROOT.TFile.Open('SR_vs_amp.root')).Get('fitFunc_SR_vs_amp')

np = 3
errSRsyst  = 0.10 # error on the slew rate

g = {}
gNoise = {}
gStoch = {}
gTot   = {}
gStoch_vs_npe = {}

g_SR_vs_bar = {}
g_SR_vs_Vov = {}

g_bestTh_vs_bar = {}
g_bestTh_vs_Vov = {}

g_Noise_vs_bar = {}
g_Stoch_vs_bar = {}

bars = {}
Vovs = {}
for sipm in sipmTypes:
    f = ROOT.TFile.Open(fnames[sipm])
    print sipm, fnames[sipm]

    #listOfKeys = [key.GetName().replace('g_deltaT_energyRatioCorr_vs_vov_','') for key in ROOT.gDirectory.GetListOfKeys() if ( 'g_deltaT_energyRatioCorr_vs_vov_bar' in key.GetName() and 'th10' in key.GetName() )]
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
for sipm in sipmTypes:
    f = ROOT.TFile.Open(fnames[sipm])

    gNoise[sipm] = {}
    gStoch[sipm] = {}
    gTot[sipm] = {}
    g[sipm] = {}
    gStoch_vs_npe[sipm] = {}

    g_SR_vs_bar[sipm] = {}
    g_SR_vs_Vov[sipm] = {}
    g_bestTh_vs_bar[sipm] = {}
    g_bestTh_vs_Vov[sipm] = {}
    g_Noise_vs_bar[sipm] = {}
    g_Stoch_vs_bar[sipm] = {}

    fPS[sipm] = {}
    for ov in Vovs[sipm]:
        if (sipm == 'HPK_unirr_LYSO528'): fPS[sipm][ov] = ROOT.TFile.Open('../plots/pulseShape_HPK528_unirr_52deg_T10C_Vov%.2f_ith1.root'%ov)
        if (sipm == 'FBK_unirr_LYSO422'): fPS[sipm][ov] = ROOT.TFile.Open('../plots/pulseShape_FBK_unirr_52deg_T10C_Vov%.2f_ith1.root'%ov)
        if (sipm == 'HPK_unirr_LYSOwithSlit'): fPS[sipm][ov] = ROOT.TFile.Open('../plots/pulseShape_HPK_unirr_LYSOwithSlit_52deg_T18C_Vov%.2f_ith1.root'%ov)
        g_SR_vs_bar[sipm][ov] = ROOT.TGraphErrors()
        g_bestTh_vs_bar[sipm][ov] = ROOT.TGraphErrors()
        g_Noise_vs_bar[sipm][ov] = ROOT.TGraphErrors()
        g_Stoch_vs_bar[sipm][ov] = ROOT.TGraphErrors()
        
    for bar in bars[sipm]:
        #g[sipm][bar] = f.Get('g_deltaT_energyRatioCorr_vs_vov_bar%02d_th10_enBin01'%bar)
        g[sipm][bar] = f.Get('g_deltaT_energyRatioCorr_bestTh_vs_vov_bar%02d_enBin01'%bar)
        gNoise[sipm][bar] = ROOT.TGraphErrors()
        gStoch[sipm][bar] = ROOT.TGraphErrors()
        gTot[sipm][bar] = ROOT.TGraphErrors()
        gStoch_vs_npe[sipm][bar] = ROOT.TGraphErrors()

        g_SR_vs_Vov[sipm][bar] = ROOT.TGraphErrors()
        g_bestTh_vs_Vov[sipm][bar] = ROOT.TGraphErrors()
        
        sigma_stoch_ref = 0
        err_sigma_stoch_ref = 0
        ov_ref = 3.5
    
        for i in range(0,g[sipm][bar].GetN()):
            ov = g[sipm][bar].GetX()[i]
            if (ov != ov_ref): continue
            sigma_tot = g[sipm][bar].GetY()[i]
            err_sigma_tot = g[sipm][bar].GetErrorY(i)
            #Npe  = LO[sipm]*4.2*PDE(ov,sipm)/PDE(3.5,sipm)
            #gain = Gain(ov, sipm)
            #sr = fitFunc_SR.Eval(gain*Npe/(Gain(3.5, sipm)*9500)) 
            g_psL = fPS[sipm][ov].Get('g_pulseShapeL_bar%02d_Vov%.2f'%(bar,ov))
            g_psR = fPS[sipm][ov].Get('g_pulseShapeR_bar%02d_Vov%.2f'%(bar,ov))
            if (g_psL==None): continue
            if (g_psR==None): continue
            g_psL.SetName('g_pulseShapeL_bar%02d_Vov%.2f'%(bar,ov))
            g_psR.SetName('g_pulseShapeR_bar%02d_Vov%.2f'%(bar,ov))
            timingThreshold = findTimingThreshold(f.Get('g_deltaT_energyRatioCorr_vs_th_bar%02d_Vov%.2f_enBin01'%(bar,ov)))
            srL,err_srL = getSlewRateFromPulseShape(g_psL, timingThreshold, np)
            srR,err_srR = getSlewRateFromPulseShape(g_psR, timingThreshold, np)
            if (srL>0 and srR>0):
                # weighted average
                sr =  ( (srL/(err_srL*err_srL) + srR/(err_srR*err_srR) ) / (1./(err_srL*err_srL) + 1./(err_srR*err_srR) ) )
                errSR = 1./math.sqrt( 1./(err_srL*err_srL)  +  1./(err_srR*err_srR) )
                errSR = errSR/sr
            if (srL>0 and srR<0):
                sr = srL
                errSR = err_srL/sr
            if (srL<0 and srR>0):
                sr = srR
                errSR = err_srR/sr
            errSR = math.sqrt(errSR*errSR+errSRsyst*errSRsyst) 
            if (sigma_tot>=sigma_noise(sr)):
                sigma_stoch_ref = math.sqrt(sigma_tot*sigma_tot - sigma_noise(sr)*sigma_noise(sr))
                err_sigma_noise = 0.5*(sigma_noise(sr*(1-errSR))-sigma_noise(sr*(1+errSR)))
                err_sigma_stoch_ref = 1./sigma_stoch_ref*math.sqrt( pow(err_sigma_tot*sigma_tot,2)+pow( sigma_noise(sr)*err_sigma_noise ,2) )
            else:
                print 'skipping bar%02d:  %.1f   %.1f'%(bar, sigma_tot,sigma_noise(sr))
                continue
                
        for i in range(0,g[sipm][bar].GetN()):
            sigma_meas = g[sipm][bar].GetY()[i]
            err   = g[sipm][bar].GetErrorY(i)
            ov = g[sipm][bar].GetX()[i]
            if (ov not in Vovs[sipm]): continue
            Npe  = LO[sipm]*4.2*PDE(ov,sipm)/PDE(3.5,sipm)
            gain = Gain(ov, sipm)
            #sr = fitFunc_SR.Eval(gain*Npe/(Gain(3.5, sipm)*9500))
            g_psL = fPS[sipm][ov].Get('g_pulseShapeL_bar%02d_Vov%.2f'%(bar,ov))
            g_psR = fPS[sipm][ov].Get('g_pulseShapeR_bar%02d_Vov%.2f'%(bar,ov))
            if (g_psL!=None): g_psL.SetName('g_pulseShapeL_bar%02d_Vov%.2f_%s'%(bar,ov,sipm))
            if (g_psR!=None): g_psR.SetName('g_pulseShapeR_bar%02d_Vov%.2f_%s'%(bar,ov,sipm))
            timingThreshold = findTimingThreshold(f.Get('g_deltaT_energyRatioCorr_vs_th_bar%02d_Vov%.2f_enBin01'%(bar,ov)))
            srL = -1
            srR = -1
            sr = -1
            err_srL = -1
            err_srR = -1
            c = ROOT.TCanvas('c_%s'%(g_psL.GetName().replace('g_pulseShapeL','pulseShape')))  
            hdummy = ROOT.TH2F('hdummy','', 100, min(g_psL.GetX())-1., 30., 100, 0., 15.)
            #hdummy = ROOT.TH2F('hdummy','', 100, min(g_psL.GetX())-1., min(g_psL.GetX())+2, 100, 0., 15.)
            hdummy.GetXaxis().SetTitle('time [ns]')
            hdummy.GetYaxis().SetTitle('amplitude [#muA]')
            hdummy.Draw()
            if (g_psL!=None): srL,err_srL = getSlewRateFromPulseShape(g_psL, timingThreshold, np, c)
            if (g_psR!=None): srR,err_srR = getSlewRateFromPulseShape(g_psR, timingThreshold, np, c) 
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
            if (srL>0 and srR<0):
                sr = srL
                errSR = err_srL
            if (srL<0 and srR>0):
                sr = srR
                errSR = err_srR
            if (srL<0 and srR<0): continue
            errSR = math.sqrt(errSR*errSR+errSRsyst*errSRsyst*sr*sr) 
            #print ov, gain, Npe, srL, srR, sr, errSR
            g_SR_vs_Vov[sipm][bar].SetPoint( g_SR_vs_Vov[sipm][bar].GetN(), ov, sr )
            g_SR_vs_Vov[sipm][bar].SetPointError( g_SR_vs_Vov[sipm][bar].GetN()-1, 0, errSR )
            g_bestTh_vs_Vov[sipm][bar].SetPoint( g_bestTh_vs_Vov[sipm][bar].GetN(), ov, timingThreshold )
            g_bestTh_vs_Vov[sipm][bar].SetPointError( g_bestTh_vs_Vov[sipm][bar].GetN()-1, 0, 0 )
            g_SR_vs_bar[sipm][ov].SetPoint( g_SR_vs_bar[sipm][ov].GetN(), bar, sr )
            g_SR_vs_bar[sipm][ov].SetPointError( g_SR_vs_bar[sipm][ov].GetN()-1, 0, errSR )
            g_bestTh_vs_bar[sipm][ov].SetPoint( g_bestTh_vs_bar[sipm][ov].GetN(), bar, timingThreshold )
            g_bestTh_vs_bar[sipm][ov].SetPointError( g_bestTh_vs_bar[sipm][ov].GetN()-1, 0, 0)
            g_Noise_vs_bar[sipm][ov].SetPoint( g_Noise_vs_bar[sipm][ov].GetN(), bar, sigma_noise(sr) )
            g_Noise_vs_bar[sipm][ov].SetPointError( g_Noise_vs_bar[sipm][ov].GetN()-1, 0,  0.5*(sigma_noise(sr*(1-errSR/sr))-sigma_noise(sr*(1+errSR/sr))) )
            gNoise[sipm][bar].SetPoint(gNoise[sipm][bar].GetN(), ov, sigma_noise(sr))
            gNoise[sipm][bar].SetPointError(gNoise[sipm][bar].GetN()-1, 0, 0.5*(sigma_noise(sr*(1-errSR/sr))-sigma_noise(sr*(1+errSR/sr))))
            # compute s_stoch as diff in quadrature between measured tRes and noise term
            if ( sigma_meas > sigma_noise(sr) ):
                s = math.sqrt(sigma_meas*sigma_meas-sigma_noise(sr)*sigma_noise(sr))
                es = 1./s * math.sqrt( pow(sigma_meas*g[sipm][bar].GetErrorY(i),2) + pow( sigma_noise(sr)*gNoise[sipm][bar].GetErrorY(gNoise[sipm][bar].GetN()-1),2) )
                gStoch_vs_npe[sipm][bar].SetPoint(gStoch_vs_npe[sipm][bar].GetN(), Npe, s)
                gStoch_vs_npe[sipm][bar].SetPointError(gStoch_vs_npe[sipm][bar].GetN()-1, 0.05*Npe, es )
            # compute stoch by scaling from 3.5 V OV
            sigma_stoch = sigma_stoch_ref/math.sqrt(  PDE(ov,sipm)/PDE(ov_ref,sipm)  )
            err_sigma_stoch = err_sigma_stoch_ref/math.sqrt( PDE(ov,sipm)/PDE(ov_ref,sipm) )
            gStoch[sipm][bar].SetPoint(gStoch[sipm][bar].GetN(), ov, sigma_stoch)
            gStoch[sipm][bar].SetPointError(gStoch[sipm][bar].GetN()-1, 0, err_sigma_stoch)
            g_Stoch_vs_bar[sipm][ov].SetPoint( g_Stoch_vs_bar[sipm][ov].GetN(), bar, sigma_stoch )
            g_Stoch_vs_bar[sipm][ov].SetPointError( g_Stoch_vs_bar[sipm][ov].GetN()-1, 0,  err_sigma_stoch)
            # tot resolution summing noise + stochastic in quadrature
            sigma_tot = math.sqrt( sigma_stoch*sigma_stoch + sigma_noise(sr)*sigma_noise(sr) )
            err_sigma_tot = 1./sigma_tot * math.sqrt( pow( err_sigma_stoch*sigma_stoch,2) + pow(sigma_noise(sr)*gNoise[sipm][bar].GetErrorY(gNoise[sipm][bar].GetN()-1),2))
            gTot[sipm][bar].SetPoint(gTot[sipm][bar].GetN(), ov, sigma_tot)
            gTot[sipm][bar].SetPointError(gTot[sipm][bar].GetN()-1, 0, err_sigma_tot)
            #print sipm,' OV = %.2f  gain = %d  Npe = %d  bar = %02d  thr = %02d  SR = %.1f   noise = %.1f    stoch = %.1f   tot = %.1f'%(ov, gain, Npe, bar, timingThreshold, sr, sigma_noise(sr), sigma_stoch, sigma_tot)

# ratio of stochatic terms at 3.5 OV
g_ratio_stoch = ROOT.TGraphErrors()
for bar in range(0,16):
    if (bar not in bars[sipmTypes[1]]): continue
    if (bar not in bars[sipmTypes[0]]): continue
    #if (bar not in [6,7,8,9,10,11,12]): continue # only bars with good SR fit
    if (gStoch[sipmTypes[0]][bar].Eval(3.5)<=0): continue
    ratio_stoch =  gStoch[sipmTypes[1]][bar].Eval(3.5)/gStoch[sipmTypes[0]][bar].Eval(3.5)
    err1 = [  gStoch[sipmTypes[1]][bar].GetErrorY(i) for i in range(0, gStoch[sipmTypes[1]][bar].GetN()) if gStoch[sipmTypes[1]][bar].GetX()[i] == 3.50]
    err0 = [  gStoch[sipmTypes[0]][bar].GetErrorY(i) for i in range(0, gStoch[sipmTypes[0]][bar].GetN()) if gStoch[sipmTypes[0]][bar].GetX()[i] == 3.50]
    if (err1 == [] or err0 == []): continue
    err_ratio_stoch = ratio_stoch * math.sqrt( pow(err1[0]/gStoch[sipmTypes[1]][bar].Eval(3.5),2) + pow(err0[0]/gStoch[sipmTypes[0]][bar].Eval(3.5),2) ) 
    print 'ratio stochastic term at 3.5 V OV = ', ratio_stoch
    g_ratio_stoch.SetPoint(g_ratio_stoch.GetN(), bar, ratio_stoch)
    g_ratio_stoch.SetPointError(g_ratio_stoch.GetN()-1, 0, err_ratio_stoch)
        
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
    c2[sipm] = {}
    hdummy2[sipm] = {}
    leg[sipm] = ROOT.TLegend(0.55,0.70,0.89,0.89)
    leg[sipm].SetBorderSize(0)
    for i,bar in enumerate(bars[sipm]):
        c1[sipm][bar] =  ROOT.TCanvas('c_timeResolution_vs_Vov_%s_bar%02d'%(sipm,bar),'c_timeResolution_vs_Vov_%s_bar%02d'%(sipm,bar),600,600)
        c1[sipm][bar].SetGridy()
        c1[sipm][bar].cd()
        hdummy1[sipm][bar] = ROOT.TH2F('hdummy1_%s_%d'%(sipm,bar),'',100,0,6,100,0,100)
        hdummy1[sipm][bar].GetXaxis().SetTitle('V_{OV} [V]')
        hdummy1[sipm][bar].GetYaxis().SetTitle('#sigma_{t} [ps]')
        hdummy1[sipm][bar].Draw()
        g[sipm][bar].SetMarkerStyle(20)
        g[sipm][bar].SetMarkerSize(1)
        g[sipm][bar].SetMarkerColor(1)
        g[sipm][bar].SetLineColor(1)
        g[sipm][bar].SetLineWidth(2)
        g[sipm][bar].Draw('plsame')
        gNoise[sipm][bar].SetLineWidth(2)
        gNoise[sipm][bar].SetLineColor(ROOT.kBlue)
        gNoise[sipm][bar].SetFillColor(ROOT.kBlue)
        gNoise[sipm][bar].SetFillColorAlpha(ROOT.kBlue,0.2)
        gNoise[sipm][bar].SetFillStyle(3004)
        gNoise[sipm][bar].Draw('E3lsame')
        gStoch[sipm][bar].SetLineWidth(2)
        gStoch[sipm][bar].SetLineColor(ROOT.kGreen+2)
        gStoch[sipm][bar].SetFillColor(ROOT.kGreen+2)
        gStoch[sipm][bar].SetFillStyle(3001)
        gStoch[sipm][bar].SetFillColorAlpha(ROOT.kGreen+2,0.2)
        gStoch[sipm][bar].Draw('E3lsame')
        #gStoch[sipm][bar].SetLineColor(ROOT.kGreen+2)
        #gStoch[sipm][bar].Draw('lsame')
        gTot[sipm][bar].SetLineWidth(2)
        gTot[sipm][bar].SetLineColor(ROOT.kRed+1)
        #gTot[sipm][bar].SetLineStyle(2)
        #gTot[sipm][bar].SetLineWidth(2)
        gTot[sipm][bar].SetFillColor(ROOT.kRed+1)
        gTot[sipm][bar].SetFillColorAlpha(ROOT.kRed+1,0.2)
        gTot[sipm][bar].SetFillStyle(3001)
        gTot[sipm][bar].Draw('E3lsame')
        if (i==0):
            leg[sipm].AddEntry(g[sipm][bar], 'data', 'PL')
            leg[sipm].AddEntry(gNoise[sipm][bar], 'noise', 'PL')
            leg[sipm].AddEntry(gStoch[sipm][bar], 'stoch', 'PL')
            leg[sipm].AddEntry(gTot[sipm][bar], 'stoch (+) noise', 'PL')
        leg[sipm].Draw('same')
        c1[sipm][bar].SaveAs(outdir+'/'+c1[sipm][bar].GetName()+'.png')

        # vs npe
        c2[sipm][bar] =  ROOT.TCanvas('c_timeResolution_vs_Npe_%s_bar%02d'%(sipm,bar),'c_timeResolution_vs_Npe_%s_bar%02d'%(sipm,bar),600,600)
        c2[sipm][bar].SetGridy()
        c2[sipm][bar].cd()
        hdummy2[sipm][bar] = ROOT.TH2F('hdummy2_%s_%d'%(sipm,bar),'',10000,1000,10000,100,0,100)
        hdummy2[sipm][bar].GetXaxis().SetTitle('Npe')
        hdummy2[sipm][bar].GetYaxis().SetTitle('#sigma_{t} [ps]')
        hdummy2[sipm][bar].Draw()
        gStoch_vs_npe[sipm][bar].SetMarkerStyle(20)
        gStoch_vs_npe[sipm][bar].SetMarkerSize(0.8)
        gStoch_vs_npe[sipm][bar].SetMarkerColor(ROOT.kGreen+2)
        gStoch_vs_npe[sipm][bar].SetLineWidth(1)
        gStoch_vs_npe[sipm][bar].SetLineColor(ROOT.kGreen+2)
        gStoch_vs_npe[sipm][bar].Draw('plsame')
        fitFun = ROOT.TF1('fitFun_%s_%.2d'%(sipm,bar),'[0]*pow(x,[1])',2000,9500)
        fitFun.SetLineColor(ROOT.kGreen+4)
        fitFun.SetParameters(30,-0.5)
        gStoch_vs_npe[sipm][bar].Fit(fitFun,'QRS')
        c2[sipm][bar].SaveAs(outdir+'/'+c2[sipm][bar].GetName()+'.png')

# SR and best threshold vs Vov
markers = { 'HPK_unirr_LYSO528' : 20 ,
            'HPK_unirr_LYSOwithSlit' : 24 ,
            'FBK_unirr_LYSO422' : 24 }

cols = { 'HPK_unirr_LYSO528' : ROOT.kBlack ,
         'HPK_unirr_LYSOwithSlit' : ROOT.kBlue ,
         'FBK_unirr_LYSO422' : ROOT.kRed }

leg2 = ROOT.TLegend(0.15,0.75,0.40,0.89)
leg2.SetBorderSize(0)
for i,bar in enumerate(bars[sipm]):
    c3[bar] = ROOT.TCanvas('c_slewRate_vs_Vov_%s_%s_bar%02d'%(sipmTypes[0], sipmTypes[1],bar),'c_slewRate_vs_Vov_%s_%s_bar%02d'%(sipmTypes[0], sipmTypes[1],bar),600,600)
    c3[bar].SetGridy()
    c3[bar].cd()
    hdummy3[bar] = ROOT.TH2F('hdummy3_%d'%(bar),'',100,0,6,100,0,35)
    hdummy3[bar].GetXaxis().SetTitle('V_{OV} [V]')
    hdummy3[bar].GetYaxis().SetTitle('slew rate at the timing thr. [#muA/ns]')
    hdummy3[bar].Draw()
    for j,sipm in enumerate(sipmTypes):
        if (i==0):
            leg2.AddEntry(g_SR_vs_Vov[sipm][bar], '%s'%sipm, 'PL')
        g_SR_vs_Vov[sipm][bar].SetMarkerStyle( markers[sipm] )
        g_SR_vs_Vov[sipm][bar].SetMarkerColor(cols[sipm])
        g_SR_vs_Vov[sipm][bar].SetLineColor(cols[sipm])
        g_SR_vs_Vov[sipm][bar].Draw('psame')
    leg2.Draw()
    c3[bar].SaveAs(outdir+'/'+c3[bar].GetName()+'.png')
    
    c4[bar] = ROOT.TCanvas('c_bestTh_vs_Vov_%s_%s_bar%02d'%(sipmTypes[0], sipmTypes[1],bar),'c_bestTh_vs_Vov_%s_%s_bar%02d'%(sipmTypes[0], sipmTypes[1],bar),600,600)
    c4[bar].SetGridy()
    c4[bar].cd()
    hdummy4[bar] = ROOT.TH2F('hdummy4_%d'%(bar),'',100,0,6,100,0,20)
    hdummy4[bar].GetXaxis().SetTitle('V_{OV} [V]')
    hdummy4[bar].GetYaxis().SetTitle('best threshold [DAC]')
    hdummy4[bar].Draw()
    for j,sipm in enumerate(sipmTypes):
        g_bestTh_vs_Vov[sipm][bar].SetMarkerStyle( markers[sipm] )
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
for ov in Vovs[sipm]:
    c5[ov] = ROOT.TCanvas('c_slewRate_vs_bar_%s_%s_Vov%.2f'%(sipmTypes[0], sipmTypes[1],ov),'c_slewRate_vs_bar_%s_%s_Vov%.2f'%(sipmTypes[0], sipmTypes[1],ov),600,600)
    c5[ov].SetGridy()
    c5[ov].cd()
    hdummy5[ov] = ROOT.TH2F('hdummy5_%d'%(ov),'',100,-0.5,15.5,100,0,35)
    hdummy5[ov].GetXaxis().SetTitle('bar')
    hdummy5[ov].GetYaxis().SetTitle('slew rate at the timing thr. [#muA/ns]')
    hdummy5[ov].Draw()
    for j,sipm in enumerate(sipmTypes):
        if (ov not in g_SR_vs_bar[sipm].keys()): continue
        g_SR_vs_bar[sipm][ov].SetMarkerStyle( markers[sipm] )
        g_SR_vs_bar[sipm][ov].SetMarkerColor(cols[sipm])
        g_SR_vs_bar[sipm][ov].SetLineColor(cols[sipm])
        g_SR_vs_bar[sipm][ov].Draw('psame')
    leg2.Draw()
    c5[ov].SaveAs(outdir+'/'+c5[ov].GetName()+'.png')

    c6[ov] = ROOT.TCanvas('c_bestTh_vs_bar_%s_%s_Vov%.2f'%(sipmTypes[0], sipmTypes[1],ov),'c_bestTh_vs_bar_%s_%s_Vov%.2f'%(sipmTypes[0], sipmTypes[1],ov),600,600)
    c6[ov].SetGridy()
    c6[ov].cd()
    hdummy6[ov] = ROOT.TH2F('hdummy6_%d'%(ov),'',100,-0.5,15.5,100,0,20)
    hdummy6[ov].GetXaxis().SetTitle('bar')
    hdummy6[ov].GetYaxis().SetTitle('timing threshold [DAC]')
    hdummy6[ov].Draw()
    for j,sipm in enumerate(sipmTypes):
        if (ov not in g_bestTh_vs_bar[sipm].keys()): continue
        g_bestTh_vs_bar[sipm][ov].SetMarkerStyle( markers[sipm] )
        g_bestTh_vs_bar[sipm][ov].SetMarkerColor(cols[sipm])
        g_bestTh_vs_bar[sipm][ov].SetLineColor(cols[sipm])
        g_bestTh_vs_bar[sipm][ov].Draw('plsame')
    leg2.Draw()        
    c6[ov].SaveAs(outdir+'/'+c6[ov].GetName()+'.png')

    c7[ov] = ROOT.TCanvas('c_noise_vs_bar_%s_%s_Vov%.2f'%(sipmTypes[0], sipmTypes[1],ov),'c_noise_vs_bar_%s_%s_Vov%.2f'%(sipmTypes[0], sipmTypes[1],ov),600,600)
    c7[ov].SetGridy()
    c7[ov].cd()
    hdummy7[ov] = ROOT.TH2F('hdummy7_%d'%(ov),'',100,-0.5,15.5,100,0,80)
    hdummy7[ov].GetXaxis().SetTitle('bar')
    hdummy7[ov].GetYaxis().SetTitle('#sigma_{t, noise} [ps]')
    hdummy7[ov].Draw()
    for j,sipm in enumerate(sipmTypes):
        if (ov not in g_Noise_vs_bar[sipm].keys()): continue
        g_Noise_vs_bar[sipm][ov].SetMarkerStyle( markers[sipm] )
        g_Noise_vs_bar[sipm][ov].SetMarkerColor(cols[sipm])
        g_Noise_vs_bar[sipm][ov].SetLineColor(cols[sipm])
        g_Noise_vs_bar[sipm][ov].Draw('psame')
    leg2.Draw()        
    c7[ov].SaveAs(outdir+'/'+c7[ov].GetName()+'.png')

    c8[ov] = ROOT.TCanvas('c_stoch_vs_bar_%s_%s_Vov%.2f'%(sipmTypes[0], sipmTypes[1],ov),'c_stoch_vs_bar_%s_%s_Vov%.2f'%(sipmTypes[0], sipmTypes[1],ov),600,600)
    c8[ov].SetGridy()
    c8[ov].cd()
    hdummy8[ov] = ROOT.TH2F('hdummy8_%d'%(ov),'',100,-0.5,15.5,100,0,80)
    hdummy8[ov].GetXaxis().SetTitle('bar')
    hdummy8[ov].GetYaxis().SetTitle('#sigma_{t, stoch} [ps]')
    hdummy8[ov].Draw()
    for j,sipm in enumerate(sipmTypes):
        if (ov not in g_Stoch_vs_bar[sipm].keys()): continue
        g_Stoch_vs_bar[sipm][ov].SetMarkerStyle( markers[sipm] )
        g_Stoch_vs_bar[sipm][ov].SetMarkerColor(cols[sipm])
        g_Stoch_vs_bar[sipm][ov].SetLineColor(cols[sipm])
        g_Stoch_vs_bar[sipm][ov].Draw('psame')
    leg2.Draw()        
    c8[ov].SaveAs(outdir+'/'+c8[ov].GetName()+'.png')
    
cc = ROOT.TCanvas('c_ratioStoch_vs_bar_%s_%s_'%(sipmTypes[0], sipmTypes[1]),'c_ratioStoch_vs_bar_%s_%s_'%(sipmTypes[0], sipmTypes[1]))
cc.cd()
hdummy = ROOT.TH2F('hdummy2','',16,-0.5,15.5,100,0.8,1.5)
hdummy.GetXaxis().SetTitle('bar')
hdummy.GetYaxis().SetTitle('ratio stochastic term')
hdummy.Draw('')
g_ratio_stoch.SetMarkerStyle(20)
g_ratio_stoch.Draw('psame')
print 'ratio of stoch. terms expected from LO  = ', math.sqrt(LO[sipmTypes[0]]/LO[sipmTypes[1]])
print 'ratio of stoch. terms measured at 3.5 V = ', g_ratio_stoch.GetMean(2)
ll = ROOT.TLine(0, math.sqrt(LO[sipmTypes[0]]/LO[sipmTypes[1]]), 15, math.sqrt(LO[sipmTypes[0]]/LO[sipmTypes[1]]))
ll.SetLineStyle(7)
ll.SetLineWidth(2)
ll.SetLineColor(ROOT.kOrange+1)
ll.Draw('same')
#g_ratio_stoch.Fit('pol0')
leg2 = ROOT.TLegend(0.55,0.70,0.89,0.89)
leg2.AddEntry(g_ratio_stoch,'ratio of stoch. term','PL')
leg2.AddEntry(ll,'sqrt(LO_{HPK}/LO_{FBK})','PL')
leg2.Draw('same')
cc.SaveAs(outdir+'/'+cc.GetName()+'.png')                        


for sipm in sipmTypes:
    print sipm, '  average stoch. term = ', g_Stoch_vs_bar[sipm][3.50].GetMean(2),' ps' 
    print sipm, '  average noise  term = ', g_Noise_vs_bar[sipm][3.50].GetMean(2),' ps'
    


raw_input('OK?')
