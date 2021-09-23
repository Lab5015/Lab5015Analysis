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
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning



labels = ['1E13_52deg_T0C',
          'unirr_52deg_T5C']

thresholds = [5,7,10,15,20,30,40,50]


# get average tRes vs threshold

g = {}
f = {}


for i,label in enumerate(labels):
    g[label] = ROOT.TGraphErrors()
    g[label].SetMarkerSize(1)
    g[label].SetMarkerColor(51+10*i)
    g[label].SetLineColor(51+10*i)

    f[label] = ROOT.TFile.Open('../plots/HPK_%s_summaryPlots.root'%label) 
    
    for thr in thresholds:
        htemp = ROOT.TH1F('htemp','',1000,0,1000)
        for bar in range(0,15):
            gg = f[label].Get('g_deltaT_energyRatioCorr_vs_th_bar%02d_Vov1.50_enBin01'%bar)
            if (gg == None): continue
            htemp.Fill( gg.Eval(thr))
        g[label].SetPoint(g[label].GetN(), thr, htemp.GetMean())
        g[label].SetPointError(g[label].GetN()-1, 0, htemp.GetRMS()/math.sqrt(16))
        htemp.Delete()

gdiff = ROOT.TGraphErrors()
gdiff.SetMarkerStyle(24)        
gdiff.SetLineStyle(2)        
for thr in thresholds:
    res0 = g[labels[0]].Eval(thr)
    res1 = g[labels[1]].Eval(thr)
    diff = math.sqrt(res0*res0-res1*res1 )
    gdiff.SetPoint(gdiff.GetN(), thr, diff)

leg = ROOT.TLegend(0.20,0.70, 0.55, 0.93)
leg.SetBorderSize(0)
leg.AddEntry(g[labels[0]], '#sigma_{bar} - HPK 1E13  ','PL')
leg.AddEntry(g[labels[1]], '#sigma_{bar} - HPK unirr.', 'PL')
leg.AddEntry(gdiff, 'diff. in quadrature','PL')


c = ROOT.TCanvas('c','c',600,600)
g[labels[0]].GetXaxis().SetTitle('threshold [DAC]')
g[labels[0]].GetYaxis().SetTitle('#sigma_{t} [ps]')
g[labels[0]].GetYaxis().SetRangeUser(0,300.)
g[labels[0]].Draw('apl')
g[labels[1]].Draw('plsame')
gdiff.Draw('plsame')
leg.Draw()


c.SaveAs('/var/www/html/TOFHIR2A/MTDTB_CERN_Jul21/ModuleCharacterization/c_tRes_vs_thr_Vov1.50_HPK.png')
c.SaveAs('/var/www/html/TOFHIR2A/MTDTB_CERN_Jul21/ModuleCharacterization/c_tRes_vs_thr_Vov1.50_HPK.pdf')

#raw_input('ok?')
