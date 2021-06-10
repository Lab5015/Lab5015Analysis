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
ROOT.gROOT.SetBatch(True)

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetTitleOffset(1.05,'Y')
ROOT.gStyle.SetLabelSize(0.04)
ROOT.gErrorIgnoreLevel = ROOT.kWarning;

plotDir = '/var/www/html/TOFHIR2A/ModuleCharacterization/scan_NPE/'

runs_dict = { 0 : [918, 920, 922, 924, 926, 928, 930, 932, 934, 936],
              1 : [919, 921, 923, 925, 927, 929, 931, 933, 935, 937] }

#runs_dict = { 0 : [918, 920, 922, 924, 928, 930, 932, 934, 936],
#              1 : [919, 921, 923, 925, 929, 931, 933, 935, 937] }

useBestThreshold = True
refTh = int(sys.argv[1])

c = ROOT.TCanvas('c_timeResolution_vs_Npe','c_timeResolution_vs_Npe',1000,800)
c.cd()
hPad = ROOT.gPad.DrawFrame(0.,0.,20000., 300.)
hPad.SetTitle(";N_{pe};#sigma_{t_{diff}} / 2 [ps]")
hPad.Draw()
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()


leg = ROOT.TLegend(0.20,0.80,0.35,0.92)
leg.SetBorderSize(0)
leg.SetFillStyle(0)


g = {}
g_bestTh = {}
fitFun = {}
fitFun2 = {}
ps = {}
ps2 = {}
for gain,runs in runs_dict.items():
        #        print gain, runs
        g[gain] = ROOT.TGraphErrors()
        for run in runs:
                inFile = ROOT.TFile.Open('/data/Lab5015Analysis/moduleCharacterization/moduleCharacterization_step2_run'+str(run)+'.root')
                
                # find best th
                bestRes = 9999
                bestTh = 10
                for th in [5, 7, 10, 15, 20, 25, 30, 40]:
                        histoRes = inFile.Get('h1_deltaT_energyRatioPhaseCorr_bar00L-R_Vov3.5_th%.02d_energyBin01'%th)
                        #histoRes = inFile.Get('h1_deltaT_energyRatioCorr_bar00L-R_Vov3.5_th%.02d_energyBin01'%th)
                        if (histoRes.GetEntries()<100 ): continue
                        fitRes = ROOT.TF1('fitRes_'+str(run), "gaus", histoRes.GetMean() -2 * histoRes.GetRMS(), histoRes.GetMean() +2 * histoRes.GetRMS())
                        fitRes.SetParameters(histoRes.GetMaximum(), histoRes.GetMean(), histoRes.GetRMS())
                        histoRes.Fit(fitRes, 'NQR')
                        histoRes.Fit(fitRes, 'NQR','', fitRes.GetParameter(1)-2*fitRes.GetParameter(2), fitRes.GetParameter(1)+2*fitRes.GetParameter(2))
                        if ( fitRes.GetParameter(2) < bestRes): 
                                bestRes = fitRes.GetParameter(2)
                                bestTh  = th

                if (useBestThreshold):
                        refTh = bestTh
                
                histoRes = inFile.Get('h1_deltaT_energyRatioPhaseCorr_bar00L-R_Vov3.5_th%.02d_energyBin01'%refTh)
                #histoRes = inFile.Get('h1_deltaT_energyRatioCorr_bar00L-R_Vov3.5_th%.02d_energyBin01'%refTh)
                histoEnergy = inFile.Get('h1_energy_bar00R_Vov3.5_th10') # ch111 used for Npe calibration
		
                fitRes = ROOT.TF1('fitRes_'+str(run), "gaus", histoRes.GetMean() -2 * histoRes.GetRMS(), histoRes.GetMean() +2 * histoRes.GetRMS())
                fitRes.SetParameters(histoRes.GetMaximum(), histoRes.GetMean(), histoRes.GetRMS())
                histoRes.Fit(fitRes, 'NQR')
                histoRes.Fit(fitRes, 'NQR','', fitRes.GetParameter(1)-2*fitRes.GetParameter(2), fitRes.GetParameter(1)+2*fitRes.GetParameter(2))
                
                fitEnergy = ROOT.TF1('fitEnergy_'+str(run), "gaus", histoEnergy.GetMean() -2 * histoEnergy.GetRMS(), histoEnergy.GetMean() +2 * histoEnergy.GetRMS())
                fitEnergy.SetParameters(histoEnergy.GetMaximum(), histoEnergy.GetMean(), histoEnergy.GetRMS())
                histoEnergy.Fit(fitEnergy, 'NQR')
                histoEnergy.Fit(fitEnergy, 'NQR','', fitEnergy.GetParameter(1)-2*fitEnergy.GetParameter(2), fitEnergy.GetParameter(1)+2*fitEnergy.GetParameter(2))

                print fitEnergy.GetParameter(1), fitEnergy.GetParameter(1)/250.*9500., bestTh
                        
                energyErr = fitEnergy.GetParameter(1)/250*9500 * 0.03
                #energyErr = fitEnergy.GetParError(1)/250*9500
                
                g[gain].SetPoint( g[gain].GetN(), fitEnergy.GetParameter(1)/250.*9500., fitRes.GetParameter(2)/2.)
                g[gain].SetPointError(g[gain].GetN()-1, energyErr, fitRes.GetParError(2)/2.)


                        



        # fit tRes vs Npe
        gg = g[gain].Clone()
        
        #fitFun[gain] = ROOT.TF1('fitFun_%d'%gain,'sqrt( [0]*[0] + [1]*[1]/x + [2]*[2]/x/x )', 2000, 20000)        
        fitFun[gain] = ROOT.TF1('fitFun_%d'%gain,'sqrt( [0]*[0] + [1]*[1]/(x/9500) + [2]*[2]/pow(x/9500,2))', 2000, 20000)        
        fitFun[gain].SetLineWidth(1)       
        fitFun[gain].SetLineColor(gain+1)       
        fitFun[gain].SetParName(0,'c')
        fitFun[gain].SetParName(1,'s')
        fitFun[gain].SetParName(2,'n')
        fitFun[gain].SetParameters(10, 30, 50)
        fitFun[gain].FixParameter(0,12./math.sqrt(2))
        g[gain].Fit(fitFun[gain],'QSR')
        g[gain].Draw('PSAME')

        # Retrieve the stat box
        c.Update()
        ps[gain] = g[gain].FindObject("stats");
        ps[gain].SetTextColor(gain+1)
        ps[gain].SetLineColor(gain+1)
        ps[gain].SetX1NDC(0.75)
        ps[gain].SetX2NDC(0.92)
        ps[gain].SetY1NDC(0.63  + gain*0.16)
        ps[gain].SetY2NDC(0.78 + gain*0.16)


        # fit with power law + const term
        fitFun2[gain] = ROOT.TF1('fitFun2_%d'%gain,'sqrt( [0]*[0] + [1]*[1]/(pow(x,[2]) * pow(x,[2])) )', 2000, 20000)        
        fitFun2[gain].SetLineWidth(1)       
        fitFun2[gain].SetLineStyle(2)       
        fitFun2[gain].SetLineColor(gain+1)       
        fitFun2[gain].SetParName(0,'c')
        fitFun2[gain].SetParName(1,'a')
        fitFun2[gain].SetParName(2,'#alpha')
        fitFun2[gain].SetParameters(10, 30, 0.5)
        fitFun2[gain].FixParameter(0, 12./math.sqrt(2))
        gg.Fit(fitFun2[gain],'SR+')  
        gg.Draw('PSAME')   
        
        # Retrieve the stat box
        c.Update()
        ps2[gain] = gg.FindObject("stats");
        ps2[gain].SetTextColor(gain+1)
        ps2[gain].SetLineColor(gain+1)
        ps2[gain].SetLineStyle(2)
        ps2[gain].SetX1NDC(0.75)
        ps2[gain].SetX2NDC(0.92)
        ps2[gain].SetY1NDC(0.30  + gain*0.16)
        ps2[gain].SetY2NDC(0.45 + gain*0.16)

        # draw graphs
        g[gain].SetMarkerStyle(20)
        g[gain].SetMarkerColor(gain+1)
        g[gain].SetLineColor(gain+1)
        g[gain].Draw('PSAME')

        leg.AddEntry(g[gain],"PreAmpGain%d"%gain, "pl")

leg.Draw("SAME")

if (useBestThreshold == False):
        c.Print(str(plotDir)+'c_tRes_vs_Npe_Vov3.5_th%.02d.png'%refTh)
        c.Print(str(plotDir)+'c_tRes_vs_Npe_Vov3.5_th%.02d.pdf'%refTh)
else:
        c.Print(str(plotDir)+'c_tRes_vs_Npe_Vov3.5_bestTh.png')
        c.Print(str(plotDir)+'c_tRes_vs_Npe_Vov3.5_bestTh.pdf')
        
		
			
