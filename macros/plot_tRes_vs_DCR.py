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


Vov = 1.50

file_dict = {'../plots/HPK_unirr_52deg_T5C_summaryPlots.root' : '',
             '../plots/HPK_1E13_52deg_T-40C_summaryPlots.root': '../plots/DCR_vs_Vov_HPK_1E13_52deg_T-40C.root',
             '../plots/HPK_1E13_52deg_T-25C_summaryPlots.root': '../plots/DCR_vs_Vov_HPK_1E13_52deg_T-25C.root',
             '../plots/HPK_1E13_52deg_T0C_summaryPlots.root'  : '../plots/DCR_vs_Vov_HPK_1E13_52deg_T0C.root'} 


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
    g = fRes.Get('g_deltaT_energyRatioCorr_bestTh_vs_bar_Vov%.02f_enBin01'%Vov)
    if (g == None): continue
    tRes[dcr] = [g.GetMean(2) , g.GetRMS(2)/math.sqrt(16) ]

#sorted(tRes.items(), key=lambda x: x[1])
print tRes


g = ROOT.TGraphErrors()

for dcr in sorted(tRes):
    print dcr, err_dcr[dcr],tRes[dcr]
    tres = math.sqrt(tRes[dcr][0]*tRes[dcr][0] - tRes[0.0][0]*tRes[0.0][0])
    tres_err = 0
    if tres > 0:
        tres_err = 1./tres * math.sqrt(pow(tRes[dcr][0]*tRes[dcr][1],2)+pow(tRes[0.0][0]*tRes[0.0][1],2))
    else:
        tres_err = tRes[dcr][1]
    g.SetPoint(g.GetN(), dcr, tres)
#    g.SetPointError(g.GetN()-1, err_dcr[dcr], tRes[dcr][1])
    g.SetPointError(g.GetN()-1, err_dcr[dcr], tres_err)


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
g.GetYaxis().SetRangeUser(0,180)
g.Draw('ap')
c.Update()
st = g.FindObject("stats")
st.SetX1NDC(0.7) 
st.SetX2NDC(0.93) 
st.SetY1NDC(0.8) #new x end position
st.SetY2NDC(0.93) #new x end position
c.Modified()

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.035)
latex.DrawLatex(0.2,0.88,'V_{ov} = %.02f V'%Vov)
latex.DrawLatex(0.2,0.83,'HPK arrays - 1e13 n_{eq}/cm^{2}')

c.SaveAs('/var/www/html/TOFHIR2A/MTDTB_CERN_Jul21/ModuleCharacterization/c_tResDCR_vs_DCR_Vov%0.02f_HPK_1E13.png'%Vov)
c.SaveAs('/var/www/html/TOFHIR2A/MTDTB_CERN_Jul21/ModuleCharacterization/c_tResDCR_vs_DCR_Vov%0.02f_HPK_1E13.pdf'%Vov)


raw_input('ok?')
