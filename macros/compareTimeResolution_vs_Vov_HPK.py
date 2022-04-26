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
#ROOT.gROOT.SetBatch(True)
ROOT.gROOT.SetBatch(False)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gStyle.SetOptStat(0)

outdir = '/eos/user/m/malberti/www/MTD/TOFHIR2X/MTDTB_CERN_Oct21/'

sipmTypes = ['HPK_1E13_LYSOtype2_T0C','HPK_1E13_LYSOtype1_T0C','HPK_1E13_LYSOtype1_T-6C']
fnames = {}
fnames = {'HPK_1E13_LYSOtype2_T0C' : '../plots/HPK_1E13_52deg_T0C_summary.root',
          'HPK_1E13_LYSOtype1_T0C' : '../plots/HPK_1E13_LYSOtype1_58deg_T0C_summary.root',
          'HPK_1E13_LYSOtype1_T-6C' : '../plots/HPK_1E13_LYSOtype1_58deg_T-6C_summary.root'}


VovsEff = {}
VovsEff['HPK_1E13_LYSOtype2_T0C'] = { 1.50 : 1.21 ,
                                      1.65 : 1.30 ,
                                      1.80 : 1.40 ,
                                      2.00 : 1.53 ,
                                      2.30 : 1.72 ,
                                      2.80 : 2.01 ,
                                      3.20 : 2.19 }   

VovsEff['HPK_1E13_LYSOtype1_T-6C']  = { 1.60 : 1.33 ,
                                        1.80 : 1.45 ,
                                        2.00 : 1.57 ,
                                        2.40 : 1.81 ,
                                        2.60 : 1.92 ,
                                        2.80 : 2.04 ,
                                        3.00 : 2.15 }


VovsEff['HPK_1E13_LYSOtype1_T0C'] = { 1.60 : 1.20 ,
                                      1.80 : 1.32 ,
                                      2.00 : 1.45 ,
                                      2.20 : 1.57 ,
                                      2.40 : 1.66 ,
                                      2.60 : 1.76 ,
                                      2.80 : 1.85 ,
                                      3.00 : 1.94 }
 

g = {}
gMax = {}
gMin = {}
g = {}
Vovs = {}
for sipm in sipmTypes:
    f = ROOT.TFile.Open(fnames[sipm])
    g[sipm] = ROOT.TGraphErrors()
    gMax[sipm] = ROOT.TGraphErrors()
    gMin[sipm] = ROOT.TGraphErrors()
    Vovs[sipm] = []
    listOfKeys = [key.GetName().replace('g_deltaT_energyRatioCorr_bestTh_vs_bar_','') for key in ROOT.gDirectory.GetListOfKeys() if key.GetName().startswith('g_deltaT_energyRatioCorr_bestTh_vs_bar_')]
    for k in listOfKeys:
        Vovs[sipm].append( float (k[3:6]) )
    Vovs[sipm].sort()    
    print sipm, Vovs[sipm]
    for i,vov in enumerate(Vovs[sipm]):
        gg = f.Get('g_deltaT_energyRatioCorr_bestTh_vs_bar_Vov%.02f_enBin01'%(vov))
        fitFun = ROOT.TF1('fitFun','pol0',0,16)
        #fitFun.SetRange(3,12)
        gg.Fit(fitFun,'QR')
        g[sipm].SetPoint(g[sipm].GetN(), VovsEff[sipm][vov], fitFun.GetParameter(0))
        #g[sipm].SetPointError(g[sipm].GetN()-1, 0, fitFun.GetParError(0))
        #print VovsEff[sipm][vov], fitFun.GetParameter(0)
        g[sipm].SetPointError( g[sipm].GetN()-1, 0, gg.GetRMS(2) )
        #g[sipm].SetPoint(g[sipm].GetN(), VovsEff[sipm][vov], gg.GetMean(2))
        #g[sipm].SetPointError(g[sipm].GetN()-1, 0, gg.GetRMS(2)/math.sqrt(gg.GetN()))
        #gMax[sipm].SetPoint(gMax[sipm].GetN(), VovsEff[sipm][vov], max(gg.GetY()))
        #gMin[sipm].SetPoint(gMin[sipm].GetN(), VovsEff[sipm][vov], min(gg.GetY()))
        
c1 =  ROOT.TCanvas('c_timeResolution_bestTh_vs_Vov','c_timeResolution_bestTh_vs_Vov',600,600)
c1.SetTickx()
c1.SetTicky()
c1.SetGridy()
c1.cd()
n = g[sipmTypes[1]].GetN()
xmax = g[sipmTypes[1]].GetX()[n-1] + 0.5
xmin = g[sipmTypes[1]].GetX()[0]   - 0.5 
ymin = 40
ymax = 140
hdummy = ROOT.TH2F('hdummy','',100,xmin,xmax,100,ymin,ymax)
hdummy.GetXaxis().SetTitle('V_{OV}^{eff} [V]')
hdummy.GetYaxis().SetTitle('#sigma_{t} [ps]')
hdummy.Draw()
leg = ROOT.TLegend(0.15,0.75,0.6,0.89)
leg.SetBorderSize(0)
for i,sipm in enumerate(sipmTypes):
    g[sipm].SetMarkerStyle(20+i)
    g[sipm].SetMarkerColor(2+i*2)
    g[sipm].SetLineColor(2+i*2)
    g[sipm].SetLineStyle(1)
    g[sipm].Draw('plsame')
    gMax[sipm].SetLineColor(2+i*2)
    gMax[sipm].SetLineStyle(2)
    gMax[sipm].Draw('lsame')
    gMin[sipm].SetLineColor(2+i*2)
    gMin[sipm].SetLineStyle(2)
    gMin[sipm].Draw('lsame')
    leg.AddEntry( g[sipm], sipm.replace('_',' ').replace('T','T=').replace('C','#circC') , 'PL')
leg.Draw('same')

#latex = ROOT.TLatex(0.6,0.85,'%s 1 MeV n_{eq}')
#latex.SetNDC()
#latex.SetTextSize(0.035)
#latex.SetTextFont(42)
#latex.Draw('same')

for c in c1,:
    c.SaveAs(outdir+c.GetName()+'_HPK_LYSOtype2_LYSOtype1.png')
    c.SaveAs(outdir+c.GetName()+'_HPK_LYSOtype2_LYSOtype1.pdf')


raw_input('OK?')
