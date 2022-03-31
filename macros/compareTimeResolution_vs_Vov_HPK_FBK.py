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


outdir = '/eos/user/m/malberti/www/MTD/TOFHIR2X/MTDTB_CERN_Oct21/'


#irr = 'unirr'
#irr = '1E13'
irr = '2E14'
sipmTypes = ['HPK_%s_T0'%irr,'FBK_%s_T0'%irr]
if (irr == '2E14'):
    #sipmTypes = ['HPK_%s_T-40'%irr,'FBK_%s_T-40'%irr, 'FBK_%s_T-32'%irr, 'FBK_%s_T-22'%irr]
    sipmTypes = ['HPK_%s_T-40'%irr,'FBK_%s_T-40'%irr]
if (irr == 'unirr'):
    sipmTypes = ['HPK_unirr_LYSO528','FBK_unirr_LYSO422', 'HPK_unirr_LYSOwithSlit']
#fnames = {}
fnames = {'HPK_unirr_LYSO528' : '../plots/HPK528_unirr_52deg_T10C_summary.root',
          'FBK_unirr_LYSO422' : '../plots/FBK_unirr_52deg_T10C_summary.root',
          'HPK_unirr_LYSOwithSlit' : '../plots/HPK_unirr_LYSOwithSlit_52deg_T18C_new_summary.root',
          'HPK_1E13_T0' : '../plots/HPK_1E13_52deg_T0C_summary.root',
          'FBK_1E13_T0' : '../plots/FBK_1E13_52deg_T0C_summary.root',
          'HPK_2E14_T-40' : '../plots/HPK_2E14_52deg_T-40C_summary.root',
          'FBK_2E14_T-32' : '../plots/FBK_2E14_52deg_T-32C_summary.root',
          'FBK_2E14_T-40' : '../plots/FBK_2E14_52deg_T-40C_summary.root',
          'FBK_2E14_T-22' : '../plots/FBK_2E14_52deg_T-22C_summary.root' }


VovsEff = {}
VovsEff['HPK_unirr_LYSO528'] = { 1.50 : 1.50 ,
                                 1.75 : 1.75 ,
                                 2.00 : 2.00 ,
                                 2.50 : 2.50 ,
                                 3.50 : 3.50 ,
                                 5.00 : 5.00 }

VovsEff['FBK_unirr_LYSO422'] = { 1.50 : 1.50 ,
                                 1.75 : 1.75 ,
                                 2.00 : 2.00 ,
                                 2.50 : 2.50 ,
                                 3.50 : 3.50 ,
                                 5.00 : 5.00 }

VovsEff['HPK_unirr_LYSOwithSlit'] = { 1.50 : 1.50 ,
                                      2.50 : 2.50 ,
                                      3.50 : 3.50 }

VovsEff['HPK_1E13_T0'] = { 1.50 : 1.21 ,
                       1.65 : 1.30 ,
                       1.80 : 1.40 ,
                       2.00 : 1.53 ,
                       2.30 : 1.72 ,
                       2.80 : 2.01 ,
                       3.20 : 2.19 }   

VovsEff['FBK_1E13_T0'] = { 1.50 : 1.19 ,
                       1.70 : 1.39 ,
                       2.00 : 1.57 ,
                       2.50 : 1.88 ,
                       3.00 : 2.17 }

VovsEff['HPK_2E14_T-40'] = { 1.60 : 1.40,
                       1.70 : 1.46,
                       1.80 : 1.50,
                       1.90 : 1.55,
                       2.00 : 1.59}

VovsEff['FBK_2E14_T-32'] = { 1.40 : 1.30 ,
                       1.70 : 1.52 ,
                       2.10 : 1.76 ,
                       2.80 : 2.07 ,
                       3.70 : 2.33 }

dV = -0.6 # Vbd shift at TB
VovsEff['FBK_2E14_T-40'] = { 1.70 : 1.57+dV,
                            2.00  : 1.78+dV ,
                            2.50  : 2.06+dV ,
                            3.00  : 2.27+dV ,
                            3.50  : 2.40+dV }
print VovsEff['FBK_2E14_T-40']

VovsEff['FBK_2E14_T-22'] = { 1.50 : 1.23 ,
                             2.00 : 1.50 ,
                             2.50 : 1.69 ,
                             2.70 : 1.76 }

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
        Vovs[sipm].append( float (k[3:7]) )
    Vovs[sipm].sort()    
    print sipm, Vovs[sipm]
    if (sipm=='HPK1E13'):Vovs[sipm].remove(1.50)
    for i,vov in enumerate(Vovs[sipm]):
        gg = f.Get('g_deltaT_energyRatioCorr_bestTh_vs_bar_Vov%.02f_enBin01'%(vov))
        fitFun = ROOT.TF1('fitFun','pol0',0,16)
        #fitFun.SetRange(3,12)
        gg.Fit(fitFun,'QR')
        print sipm, VovsEff[sipm][vov], fitFun.GetParameter(0)
        g[sipm].SetPoint(g[sipm].GetN(), VovsEff[sipm][vov], fitFun.GetParameter(0))
        g[sipm].SetPointError(g[sipm].GetN()-1, 0, fitFun.GetParError(0))
        #g[sipm].SetPointError( g[sipm].GetN()-1, 0, gg.GetRMS(2) )
        #g[sipm].SetPoint(g[sipm].GetN(), VovsEff[sipm][vov], gg.GetMean(2))
        #g[sipm].SetPointError(g[sipm].GetN()-1, 0, gg.GetRMS(2)/math.sqrt(gg.GetN()))
        gMax[sipm].SetPoint(gMax[sipm].GetN(), VovsEff[sipm][vov], max(gg.GetY()))
        gMin[sipm].SetPoint(gMin[sipm].GetN(), VovsEff[sipm][vov], min(gg.GetY()))
        
c1 =  ROOT.TCanvas('c_timeResolution_bestTh_vs_Vov','c_timeResolution_bestTh_vs_Vov',600,600)
c1.SetGridy()
c1.cd()
n = g[sipmTypes[1]].GetN()
xmax = g[sipmTypes[1]].GetX()[n-1] + 0.5
xmin = g[sipmTypes[1]].GetX()[0]   - 0.5 
ymin = 70
ymax = 180
if (irr=='unirr'):
    ymin = 20
    ymax = 100
hdummy = ROOT.TH2F('hdummy','',100,xmin,xmax,100,ymin,ymax)
hdummy.GetXaxis().SetTitle('V_{OV}^{eff} [V]')
hdummy.GetYaxis().SetTitle('#sigma_{t} [ps]')
hdummy.Draw()
#leg = ROOT.TLegend(0.15,0.60,0.45,0.89)
leg = ROOT.TLegend(0.55,0.60,0.89,0.89)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
for i,sipm in enumerate(sipmTypes):
    g[sipm].SetMarkerStyle(20+i)
    g[sipm].SetMarkerSize(1)
    g[sipm].SetMarkerColor(2+i*2)
    g[sipm].SetLineColor(2+i*2)
    g[sipm].SetLineStyle(1)
    g[sipm].SetLineWidth(2)
    g[sipm].Draw('plsame')
    gMax[sipm].SetLineColor(2+i*2)
    gMax[sipm].SetLineStyle(2)
    gMin[sipm].SetLineColor(2+i*2)
    gMin[sipm].SetLineStyle(2)
    if (irr != 'unirr'):
        #gMax[sipm].Draw('lsame')
        #gMin[sipm].Draw('lsame')
        leg.AddEntry( g[sipm], sipm.replace(irr,'').replace('_',' ').replace('T','T=')+'^{o}C' , 'PL')
    else:
        leg.AddEntry( g[sipm], sipm.replace(irr,''), 'PL')
leg.Draw('same')

#latex = ROOT.TLatex(0.6,0.85,'%s 1 MeV n_{eq}'%irr)
latex = ROOT.TLatex(0.15,0.82,'%s'%irr)
if (irr == 'unirr'):
    latex = ROOT.TLatex(0.15,0.82,'non-irradiated')
latex.SetNDC()
latex.SetTextSize(0.045)
latex.SetTextFont(42)
latex.Draw('same')

for c in [c1]:
    c.SaveAs(outdir+c.GetName()+'_%s.png'%irr)
    c.SaveAs(outdir+c.GetName()+'_%s.pdf'%irr)


#outfile = ROOT.TFile('timeResolution_averaged_vs_VOV_%s.root'%irr,'recreate')

    
#raw_input('OK?')
