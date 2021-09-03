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
ROOT.gStyle.SetTitleOffset(1.2,'Y')
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


Vov = 1.25

file_dict = {'../plots/HPK_unirr_52deg_T5C_summaryPlots.root': '',
             '../plots/HPK_1E13_52deg_T-40C_summaryPlots.root': '../plots/DCR_vs_Vov_HPK_1E13_52deg_T-40C.root',
             '../plots/HPK_1E13_52deg_T-25C_summaryPlots.root': '../plots/DCR_vs_Vov_HPK_1E13_52deg_T-25C.root',
             '../plots/HPK_1E13_52deg_T0C_summaryPlots.root':'../plots/DCR_vs_Vov_HPK_1E13_52deg_T0C.root'} 


tRes = {}
err_dcr = {}
for f1,f2 in file_dict.items():
    print f1, f2
    if ('unirr' not in f1):
        fDCR = ROOT.TFile.Open(f2)
        dcr = 0.5 * ( (fDCR.Get('g_DCR_vs_Vov_ch3')).Eval(Vov) + (fDCR.Get('g_DCR_vs_Vov_ch4')).Eval(Vov))
        err_dcr[dcr] = 0.5 * abs(( (fDCR.Get('g_DCR_vs_Vov_ch3')).Eval(Vov) - (fDCR.Get('g_DCR_vs_Vov_ch4')).Eval(Vov)) )
    else:
        dcr = 0
        err_dcr[dcr] = 0.0
    fRes = ROOT.TFile.Open(f1)         
    #tRes[dcr] = {}
    tRes[dcr] = [(fRes.Get('g_deltaT_energyRatioCorr_bestTh_vs_bar_Vov%.02f_enBin01'%Vov)).GetMean(2) ,
                 (fRes.Get('g_deltaT_energyRatioCorr_bestTh_vs_bar_Vov%.02f_enBin01'%Vov)).GetRMS(2)/math.sqrt(16) ]

#sorted(tRes.items(), key=lambda x: x[1])
print tRes


g = ROOT.TGraphErrors()

for dcr in sorted(tRes):
    print dcr, err_dcr[dcr],tRes[dcr]
    tres = math.sqrt(tRes[dcr][0]*tRes[dcr][0] - tRes[0.0][0]*tRes[0.0][0])
    g.SetPoint(g.GetN(), dcr, tres)
#    g.SetPointError(g.GetN()-1, 0.12*dcr, tRes[dcr][1])
    g.SetPointError(g.GetN()-1, err_dcr[dcr], tRes[dcr][1])


fitFun = ROOT.TF1('fitFun','[0] * pow(x,[1])',0,100)
fitFun.SetLineStyle(2)
fitFun.SetLineColor(2)
fitFun.SetParameter(0,20)
fitFun.SetParameter(1,0.5)
g.Fit(fitFun)


c = ROOT.TCanvas('c','c',600,600)
g.SetMarkerStyle(20)
g.SetMarkerSize(1)
g.GetXaxis().SetTitle('DCR [GHz]')
g.GetYaxis().SetTitle('#sigma_{DCR} [ps]')
g.GetYaxis().SetRangeUser(0,150)
g.Draw('ap')
latex = ROOT.TLatex()
latex.SetNDC()
latex.DrawLatex(0.2,0.75,'Vov = %.02f V'%Vov)

c.SaveAs('c_tResDCR_vs_DCR_Vov%0.02f.png'%Vov)
c.SaveAs('c_tResDCR_vs_DCR_Vov%0.02f.pdf'%Vov)


raw_input('ok?')
