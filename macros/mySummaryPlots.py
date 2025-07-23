import os
import shutil
import glob
import math
import array
import sys
import time
import argparse
import json
import numpy as np
import ctypes
from array import array

parser = argparse.ArgumentParser(description='Module characterization summary plots')
parser.add_argument("-i",  "--inputLabels",   required=True, type=str, help="comma-separated list of input labels")
parser.add_argument("-m",  "--resMode",       required=True, type=int, help="resolution mode: 2 - tDiff, 1 - tAve")
parser.add_argument("-o",  "--outFolder",     required=True, type=str, help="out folder")
args = parser.parse_args()


import ROOT
import CMS_lumi, tdrstyle

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.055,'X')
ROOT.gStyle.SetLabelSize(0.055,'Y')
ROOT.gStyle.SetTitleSize(0.07,'X')
ROOT.gStyle.SetTitleSize(0.07,'Y')
ROOT.gStyle.SetTitleOffset(1.05,'X')
ROOT.gStyle.SetTitleOffset(1.1,'Y')
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.045)
ROOT.gStyle.SetPadTopMargin(0.07)
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

def getTimeResolution(h1_deltaT,hit):
   
   tRes = [-1,-1]

   h1_deltaT.GetListOfFunctions().Clear()

   h1_deltaT.GetXaxis().SetRangeUser(h1_deltaT.GetMean() - 5*h1_deltaT.GetRMS(), h1_deltaT.GetMean() + 5*h1_deltaT.GetRMS())  #commentato da me
                    
   fitFunc = ROOT.TF1('fitFunc','gaus',-10000, 10000)
   fitFunc.SetLineColor(ROOT.kBlue+3)
   fitFunc.SetLineWidth(2)
   fitFunc.SetParameters(h1_deltaT.GetMaximum(),h1_deltaT.GetMean(), h1_deltaT.GetRMS())
   
   fitXMin = h1_deltaT.GetBinCenter(h1_deltaT.GetMaximumBin()) - 200 #prima era -200
   fitXMax = h1_deltaT.GetBinCenter(h1_deltaT.GetMaximumBin()) + 200.
   #fitXMin = h1_deltaT.GetMean() - 3*h1_deltaT.GetRMS()
   #fitXMax = h1_deltaT.GetMean() + 3*h1_deltaT.GetRMS()
   fitFunc.SetRange(fitXMin, fitXMax)
   h1_deltaT.Fit('fitFunc','QNRL','', fitXMin, fitXMax)
   #fitFunc.SetRange(fitFunc.GetParameter(1) - 3.0*fitFunc.GetParameter(2), fitFunc.GetParameter(1) + 3.0*fitFunc.GetParameter(2))

   fitFunc.SetRange(fitFunc.GetParameter(1) - 1.0*fitFunc.GetParameter(2), fitFunc.GetParameter(1) + 1.0*fitFunc.GetParameter(2)) #commentato da me
   h1_deltaT.Fit('fitFunc','QNRL')

   fitFunc.SetRange(fitFunc.GetParameter(1) - 2.5*fitFunc.GetParameter(2), fitFunc.GetParameter(1) + 2.5*fitFunc.GetParameter(2)) #commentato da me
   h1_deltaT.Fit('fitFunc','QRSL+')

   if(False):
    legend = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)

    name= h1_deltaT.GetName()
    c = ROOT.TCanvas("c_"+name, "", 800, 600)
    h1_deltaT.Draw()
    fitFunc.Draw("same")
    c.SetGrid()
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    legend.AddEntry(fitFunc, f"single gauss", "l")
    legend.Draw()
    c.Update()
    c.SaveAs(f"/eos/home-f/ftonetto/www/MTD_TB_CERN_Sep23/run6289_timeDiff_bis/summaryPlots/energy/A_"+name+".png")
    c.SaveAs(f"/eos/home-f/ftonetto/www/MTD_TB_CERN_Sep23/run6289_timeDiff_bis/summaryPlots/energy/A_"+name+".pdf")
    c.Close()
   
   #if (fitFunc==None): continue                    
   #if (fitFunc.GetParameter(2) > 1000): continue
   #if (fitFunc.GetParameter(2) < 20): continue
   #if (fitFunc.GetParError(2) > 200): continue
   tRes = [ fitFunc.GetParameter(2),fitFunc.GetParError(2)]
   return tRes



def getTimeResolution_S(h1_deltaT):
    tRes=[-1,-1,-1,-1]

    h1_deltaT.GetListOfFunctions().Clear()
    h1_deltaT.GetXaxis().SetRangeUser(h1_deltaT.GetMean() - 5*h1_deltaT.GetRMS(), h1_deltaT.GetMean() + 5*h1_deltaT.GetRMS())

    double_gauss = "[0]*TMath::Gaus(x, [1], [2], false) + [3]*TMath::Gaus(x, [1], [4], false)"
    fitFunc = ROOT.TF1('fitFunc', double_gauss, -10000, 10000)
    fitFunc.SetLineColor(ROOT.kRed+3)
    fitFunc.SetLineWidth(2)
    fitFunc.SetParameters(h1_deltaT.GetMaximum(), h1_deltaT.GetMean(), 120,h1_deltaT.GetMaximum()/100,200)
    #fitFunc.SetParLimits(2,80,160)
    #fitFunc.SetParLimits(4,160,240)

    fitFunc.SetParName(0, "A_{1}")
    fitFunc.SetParName(1, "#mu")
    fitFunc.SetParName(2, "#sigma_{1}")
    fitFunc.SetParName(3, "A_{2}")
    fitFunc.SetParName(4, "#sigma_{2}")

    
    #h1_deltaT.SetStats(False)

    fitXMin= h1_deltaT.GetBinCenter(h1_deltaT.GetMaximumBin()) - 400
    fitXMax= h1_deltaT.GetBinCenter(h1_deltaT.GetMaximumBin()) + 400.
    fitFunc.SetRange(fitXMin, fitXMax)
    h1_deltaT.Fit('fitFunc','QRL','', fitXMin, fitXMax)

    #fitFunc.SetRange(fitFunc.GetParameter(1) - 1.0*fitFunc.GetParameter(2), fitFunc.GetParameter(1) + 1.0*fitFunc.GetParameter(2)) 
    #h1_deltaT.Fit('fitFunc','QNRL')

    #fitFunc.SetRange(fitFunc.GetParameter(1) - 2.5*fitFunc.GetParameter(2), fitFunc.GetParameter(1) + 2.5*fitFunc.GetParameter(2)) 
    #h1_deltaT.Fit('fitFunc','QRSL+')

    
    legend = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)
    
    name= h1_deltaT.GetName()
    c = ROOT.TCanvas("c_"+name, "", 800, 600)
    h1_deltaT.Draw()
    fitFunc.Draw("same")
    c.SetGrid()
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    legend.AddEntry(fitFunc, f"double gauss", "l")
    legend.Draw()
    c.Update()
    c.SaveAs(f"/eos/home-f/ftonetto/www/MTD_TB_CERN_Sep23/run6289_timeDiff_bis/summaryPlots/energy/"+name+".png")
    c.SaveAs(f"/eos/home-f/ftonetto/www/MTD_TB_CERN_Sep23/run6289_timeDiff_bis/summaryPlots/energy/"+name+".pdf")
    c.Close()

    tRes = [ fitFunc.GetParameter(2),fitFunc.GetParError(2), fitFunc.GetParameter(4),fitFunc.GetParError(4)]
    return tRes




# INPUT
inputdir ='/eos/home-f/ftonetto/Lab5015Analysis/A_output'
source = 'TB'

# OUTPUT
outdir  = '/eos/home-f/ftonetto/www/MTD_TB_CERN_Sep23/'
outdir=outdir+args.outFolder
outFileName = inputdir+'/summaryPlots_'+args.outFolder+'.root'
print('Saving root file ', outFileName)
print('Saving plots in ', outdir)
outfile = ROOT.TFile(outFileName, 'RECREATE' )

# --- colors
cols = { 0.50 : 51,
         0.60 : 51+8,
         0.80 : 51+16,
         1.00 : 51+24,
         1.25 : 51+48, #51+32 
         1.50 : 51+40,
         2.00 : 51+44,
         2.50 : 51+46,
         3.00 : 51+47,
         3.50 : 51+32 #51+48
}

# create files list
label_list = (args.inputLabels.split(','))
print(label_list)

# resolution mode : 0 : /1;  1: /sqrt(2) if CTR, 2: /2 if TDiff
kscale = 2.
if (args.resMode == 0): kscale = 1
if (args.resMode == 1): kscale = math.sqrt(2)

peaks=[0]
enBins=[1]

# --- prepare output dir 
if (os.path.isdir(outdir) == False):
    os.system('mkdir %s'%outdir)
    os.system('cp /eos/home-f/ftonetto/www/index.php %s'%outdir) #copia file index.php per visualizzazione web
    os.system('mkdir %s/summaryPlots/'%outdir)
    os.system('cp /eos/home-f/ftonetto/www/index.php %s/summaryPlots/'%outdir)
    os.system('mkdir %s/summaryPlots/tot/'%outdir)
    os.system('cp /eos/home-f/ftonetto/www/index.php %s/summaryPlots/tot/'%outdir)
    os.system('mkdir %s/summaryPlots/energy/'%outdir)
    os.system('cp /eos/home-f/ftonetto/www/index.php %s/summaryPlots/energy/'%outdir)

    #os.system('mkdir %s/summaryPlots/energy/singleGauss/'%outdir)
    #os.system('cp /eos/home-f/ftonetto/www/index.php %s/summaryPlots/energy/singleGauss/'%outdir)

    os.system('mkdir %s/summaryPlots/timeResolution/'%outdir)
    os.system('cp /eos/home-f/ftonetto/www/index.php %s/summaryPlots/timeResolution/'%outdir)
    os.system('mkdir %s/summaryPlots/timeResolution/fits/'%outdir)
    os.system('cp /eos/home-f/ftonetto/www/index.php %s/summaryPlots/timeResolution/fits/'%outdir)

# -- ref threhsold
thRef = 11

# -- get list of bars, Vovs, thresholds to be analyzed
bars = []
thresholds = []
Vovs = [] 
Hits = [1,2,3]
for label in label_list:
    inputFile = None
    inputFile = ROOT.TFile.Open(inputdir+'/moduleCharacterization_step2_%s.root'%label)
    listOfKeys = [key.GetName().replace('h1_deltaT_totRatioCorr_','') for key in ROOT.gDirectory.GetListOfKeys() if key.GetName().startswith('h1_deltaT_totRatioCorr_bar')]
    for k in listOfKeys:
        barNum = int (k.split('_')[0][3:5])
        bars.append(barNum)
        vov = float (k.split('_')[1][3:7])
        Vovs.append(vov)
        thr = int (k.split('_')[2][2:4])
        thresholds.append(thr)
# remove duplicates
bars = [i for n, i in enumerate(bars) if i not in bars[:n]]
Vovs = [i for n, i in enumerate(Vovs) if i not in Vovs[:n]]
thresholds = [i for n, i in enumerate(thresholds) if i not in thresholds[:n]]
bars.sort()
Vovs.sort()
thresholds.sort()

goodBars = {}
VovsEff = {}
plots_label = ''

for vov in Vovs:
    VovsEff[vov] = vov
    goodBars[vov] = bars

print('Vovs:',Vovs)
print('Bars:',bars)
print('good bars:', goodBars)
print('Vovs:',Vovs)
print('thresholds:', thresholds)


# --- Summary graphs
g_deltaT_corr_vs_bars={} #g[Vov,th,energyBin,hits]
g_deltaT_bestth_vs_bars={}
g_deltaT_single_vs_bars={}
new_selection = True
# --- time resolution
tRes_spread={}
tRes_best={}
tRes_all={}
if(new_selection):
    for bar in bars:
        for vov in Vovs:
            for enBin in enBins:
                for hit in Hits:
                    tRes_best[bar, vov, enBin,hit] = [9999, 9999]
else:
    for bar in bars:
        for vov in Vovs:
            for enBin in enBins:
                tRes_best[bar, vov, enBin] = [9999, 9999]


if(new_selection):
    #inizializzo i plot
    for vov in Vovs:
        for thr in thresholds:
            for enBin in enBins:
                g_deltaT_single_vs_bars[vov,thr,enBin,1]=ROOT.TGraphErrors()
                g_deltaT_single_vs_bars[vov,thr,enBin,2]=ROOT.TGraphErrors()
                for hit in Hits:
                    g_deltaT_corr_vs_bars[vov,thr,enBin,hit]=ROOT.TGraphErrors()
                    g_deltaT_bestth_vs_bars[vov,enBin,hit]=ROOT.TGraphErrors()

    # --- Read the histograms from moduleCharacterization_step2 file
    for label in label_list:
        print(label)
        inputFile == None
        inputFile = ROOT.TFile.Open(inputdir+'/moduleCharacterization_step2_%s.root'%label)
        print("drawing histo and fitting...")
        for bar in bars:
            for vov in Vovs:
                for thr in thresholds:
                    tRes_corr={}
                    for enBin in enBins:
                        for hit in Hits:
                            #h1_deltaT_Corr = inputFile.Get('h1_deltaT__bars%02d-%02d_Vov%.2f_th%02d_energyBin%02d_%02dHits'%(bar,bar+1, vov, thr, enBin,hit))
                            #h1_deltaT_Corr =inputFile.Get('h1_deltaT_energyRatioPhaseCorr_bars%02d-%02d_Vov%.2f_th%02d_energyBin%02d_%02dHits'%(bar,bar+1, vov, thr, enBin,hit))
                            #h1_deltaT_Corr =inputFile.Get('h1_deltaT_energyRatiCorr_ave_bars%02d-%02d_Vov%.2f_th%02d_energyBin%02d_%02dHits'%(bar,bar+1, vov, thr, enBin,hit))
                             
                            #h1_deltaT_w = inputFile.Get('h1_deltaT_w_bars%02d-%02d_Vov%.2f_th%02d_energyBin%02d_%02dHits'%(bar,bar+1, vov, thr, enBin,hit))
                            h1_deltaT_Corr = inputFile.Get('h1_deltaT_energyRatioPhaseCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d_%02dHits'%(bar, vov, thr, enBin,hit))
                            #h1_deltaT_Corr = inputFile.Get('h1_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d_%02dHits'%(bar, vov, thr, enBin,hit))
                            
                            if not h1_deltaT_Corr:
                                #if hit==1 or hit==3: continue
                                print(f"Warning: histogram not found for bar{bar:02d}, Vov={vov:.2f}, thr={thr}, enBin={enBin}, hit={hit}")
                                continue 

                            #if(hit==1):
                                #tRes_S, s_tRes_S,tRes_D,s_tRes_D = getTimeResolution_S(h1_deltaT_Corr)
                                #index=g_deltaT_single_vs_bars[vov,thr,enBin,1].GetN()
                                #g_deltaT_single_vs_bars[vov,thr,enBin,1].SetPoint(index,bar,tRes_S/kscale)
                                #g_deltaT_single_vs_bars[vov,thr,enBin,1].SetPointError(index,0,s_tRes_S/kscale)
                                #index=g_deltaT_single_vs_bars[vov,thr,enBin,2].GetN()
                                #g_deltaT_single_vs_bars[vov,thr,enBin,2].SetPoint(index,bar,tRes_D/kscale)
                                #g_deltaT_single_vs_bars[vov,thr,enBin,2].SetPointError(index,0,s_tRes_D/kscale)

                            tRes_corr[enBin,hit]=getTimeResolution(h1_deltaT_Corr,hit)
                            tRes, s_tRes= tRes_corr[enBin,hit]

                            tRes_all[bar,vov,thr,enBin,hit] = [0.5*tRes,0.5*s_tRes]

                            index=g_deltaT_corr_vs_bars[vov,thr,enBin,hit].GetN()
                            g_deltaT_corr_vs_bars[vov,thr,enBin,hit].SetPoint(index,bar,tRes/kscale)
                            g_deltaT_corr_vs_bars[vov,thr,enBin,hit].SetPointError(index,0,s_tRes/kscale)
                            if(tRes_corr[enBin,hit][0]<tRes_best[bar,vov,enBin,hit][0]):
                                tRes_best[bar,vov,enBin,hit]=tRes_corr[enBin,hit]
                            
                                
        for bar in bars:
            for vov in Vovs:
                for enBin in enBins:
                    for hit in Hits:
                        if(tRes_best[bar,vov,enBin,hit][0]!=9999):
                            index=g_deltaT_bestth_vs_bars[vov,enBin,hit].GetN()
                            g_deltaT_bestth_vs_bars[vov,enBin,hit].SetPoint(index,bar,tRes_best[bar,vov,enBin,hit][0]/kscale)
                            g_deltaT_bestth_vs_bars[vov,enBin,hit].SetPointError(index,0,tRes_best[bar,vov,enBin,hit][1]/kscale)
    
    print("Average time resolution for: ")
    for key, graph in g_deltaT_corr_vs_bars.items():
        vov, thr, enBin, hit = key
        label=f"Vov{vov}_th{thr:02d}_enBin{enBin}_{hit}Hits"

        c=ROOT.TCanvas(f"c_timeRes_corr_"+label,"",800,600)
        c.SetGrid()
        c.SetLeftMargin(0.15)
        c.SetBottomMargin(0.15)
        graph.SetTitle("")
        graph.GetXaxis().SetTitle("Bar")
        graph.GetYaxis().SetTitle("#sigma_{t} [ps]")
        graph.GetYaxis().SetRangeUser(0,140)
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(1.2)
        graph.SetLineWidth(2)
        graph.SetLineColor(ROOT.kRed)
        graph.SetMarkerColor(ROOT.kRed)
        graph.Draw("AP")
        fitRes=ROOT.TF1("fitRes",'pol0',1,13)
        graph.Fit(fitRes,"QRN")#quite range no-store
        tRes_ave=fitRes.GetParameter(0)
        
        c.Update()
        c.SaveAs(outdir+"/summaryPlots/timeResolution/c_timeRes_corr_"+label+".png")
        c.SaveAs(outdir+"/summaryPlots/timeResolution/c_timeRes_corr_"+label+".pdf")
        
        c2=ROOT.TCanvas(f"c_tRes_spread_"+label,"",800,600)
        c2.SetGrid()
        c2.SetLeftMargin(0.15)
        c2.SetBottomMargin(0.15)
        histo=ROOT.TH1F("tRes_spread"+label,"",10,-0.5,0.5)
        for i in range(1,13):
            x = ctypes.c_double()
            tRes = ctypes.c_double()
            graph.GetPoint(i, x, tRes)
            tRes=tRes.value
            if(tRes==0):
                spread=0
            else:
                spread=(tRes-tRes_ave)/tRes
            histo.Fill(spread)
        histo.GetXaxis().SetTitle("#frac{#sigma_{t}-<#sigma_{t}>}{#sigma_{t}}")
        histo.Draw()
        latex = ROOT.TLatex()
        latex.SetNDC()  # coordinate normalizzate (0-1)
        latex.SetTextSize(0.04)
        latex.SetTextColor(ROOT.kRed)
        latex.DrawLatex(0.15, 0.85, f"RMS: {histo.GetRMS():0.4f}")#(in alto a sinistra)
        c2.Update()
        c2.SaveAs(outdir+"/summaryPlots/timeResolution/c_timeRes_spread_"+label+".png")
        c2.SaveAs(outdir+"/summaryPlots/timeResolution/c_timeRes_spread_"+label+".pdf")
        if(tRes_ave==0 ): continue
        else:
            print(f"Vov{vov} th{thr} enBin{enBin} hits{hit}: {tRes_ave:.01f} [ps] spread (RMS) of {100*histo.GetRMS():0.2f}%")
    
    for vov in Vovs:
        for thr in thresholds:
            for enBin in enBins:

                label=f"Vov{vov}_th{thr:02d}_enBin{enBin}"
                c3=ROOT.TCanvas(f"c_timeRes_singles_"+label,"",800,600)
                c3.SetGrid()
                c3.SetLeftMargin(0.15)
                c3.SetBottomMargin(0.15)

                legend = ROOT.TLegend(0.65, 0.15, 0.88, 0.33)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.SetTextSize(0.03)
                colors01= [ROOT.kRed, ROOT.kBlue]
                for i in range(1,3):
                    key = (vov, thr, enBin, i)
                    key2= (vov,thr,enBin,1)
                    graph01=g_deltaT_single_vs_bars.get(key)
                    if not graph01 :
                        continue
                    graph01.SetMarkerStyle(20)
                    graph01.SetMarkerSize(1.1)
                    graph01.SetLineWidth(2)
                    graph01.SetLineColor(colors01[i-1])
                    graph01.SetMarkerColor(colors01[i-1])
                    graph02=g_deltaT_corr_vs_bars.get(key2)
                    if not graph02: 
                        continue
                    graph02.SetMarkerStyle(20)
                    graph02.SetMarkerSize(1.1)
                    graph02.SetLineWidth(2)
                    graph02.SetLineColor(ROOT.kBlack)
                    graph02.SetMarkerColor(ROOT.kBlack)
                    if i == 1:
                        legend.AddEntry(graph02, f"single", "ep")
                        graph01.Draw("AP")
                        graph01.GetXaxis().SetTitle("Bar")
                        graph01.GetYaxis().SetTitle("#sigma_{t} [ps]")
                        graph01.GetYaxis().SetRangeUser(0, 150)
                        graph01.SetTitle("")
                    else:
                        graph01.Draw("P")

                    graph02.Draw("P")
                    
                    legend.AddEntry(graph01, f"maybe {i}Hits", "ep")
                    
                
                legend.Draw()
                c3.Update()
                c3.SaveAs(outdir + "/summaryPlots/timeResolution/c_timeRes_singles_" + label + ".png")
                c3.SaveAs(outdir + "/summaryPlots/timeResolution/c_timeRes_singles_" + label + ".pdf")    

    for vov in Vovs:
        for enBin in enBins:
            for hit in Hits:
                c01=ROOT.TCanvas(f"c_timeRes_corr_Vov{vov}_ALLth_enBin{enBin}_{hit}Hits","",800,600)
                c01.SetGrid()
                c01.SetLeftMargin(0.15)
                c01.SetBottomMargin(0.15)

                legend = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.SetTextSize(0.03)
                colors01= [ROOT.kRed, ROOT.kBlue, ROOT.kBlack]
                for i,thr in enumerate(thresholds):
                    key=(vov,thr,enBin,hit)
                    graph01=g_deltaT_corr_vs_bars.get(key)
                    if not graph01:
                        continue
                    graph01.SetMarkerStyle(20)
                    graph01.SetMarkerSize(1.1)
                    graph01.SetLineWidth(2)
                    graph01.SetLineColor(colors01[i])
                    graph01.SetMarkerColor(colors01[i])

                    if i == 0:
                        graph01.Draw("AP")
                        graph01.GetXaxis().SetTitle("Bar")
                        graph01.GetYaxis().SetTitle("#sigma_{t} [ps]")
                        graph01.GetYaxis().SetRangeUser(0, 140)
                        graph01.SetTitle("")
                    else:
                        graph01.Draw("P")

                    legend.AddEntry(graph01, f"thr {thr}", "ep")
                
                legend.Draw()
                label = f"c_timeRes_corr_ALLthr_Vov{vov}_enBin{enBin}_{hit}Hits"
                c01.Update()
                c01.SaveAs(outdir + "/summaryPlots/timeResolution/" + label + ".png")
                c01.SaveAs(outdir + "/summaryPlots/timeResolution/" + label + ".pdf")
            
        
    print("Average time resolution at best th for: ")
    for key,graph in g_deltaT_bestth_vs_bars.items():
        vov, enBin, hit = key
        leg = ROOT.TLegend(0.20, 0.90, 0.60, 0.70)
        c2=ROOT.TCanvas(f"c_timeRes_bestTh_corr_Vov{vov}_enBin{enBin}_{hit}Hits","",800,600)
        c2.SetGrid()
        c2.SetLeftMargin(0.15)
        c2.SetBottomMargin(0.15)
        graph.SetTitle("")
        graph.GetXaxis().SetTitle("Bar")
        graph.GetYaxis().SetTitle("#sigma_{t} [ps]")
        graph.GetYaxis().SetRangeUser(0,140)
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(1.2)
        graph.SetLineWidth(2)
        graph.SetLineColor(ROOT.kGreen)
        graph.SetMarkerColor(ROOT.kGreen)

        graph.Draw("AP")
        fitRes=ROOT.TF1("fitRes",'pol0',1,13)
        graph.Fit(fitRes,"QRN")#quite range no-store
        #print(f"Vov{vov} enBin{enBin} hits{hit}: {fitRes.GetParameter(0):.01f} [ps]")
        tRes_ave=fitRes.GetParameter(0)

        line = ROOT.TLine(1,tRes_ave ,13,tRes_ave)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(2)  # linea tratteggiata
        line.Draw("same")
        label=f"c_timeRes_corr_Vov{vov}_best_th_enBin{enBin}_{hit}Hits"
        c2.SaveAs(outdir+"/summaryPlots/timeResolution/"+label+".png")
        c2.SaveAs(outdir+"/summaryPlots/timeResolution/"+label+".pdf")
        label=f"Vov{vov}_bestTh_enBin{enBin}_{hit}Hits"
        c2=ROOT.TCanvas(f"c_tRes_spread_"+label,"",800,600)
        c2.SetGrid()
        c2.SetLeftMargin(0.15)
        c2.SetBottomMargin(0.15)
        histo=ROOT.TH1F("tRes_spread"+label,"",20,-0.5,0.5)
        for i in range(1,14):
            x = ctypes.c_double()
            tRes = ctypes.c_double()
            graph.GetPoint(i, x, tRes)
            tRes=tRes.value
            if(tRes==0):
                spread=0
            else:
                spread=(tRes-tRes_ave)/tRes
            histo.Fill(spread)
        histo.GetXaxis().SetTitle("#frac{#sigma_{t}-<#sigma_{t}>}{#sigma_{t}}")
        histo.Draw()
        latex = ROOT.TLatex()
        latex.SetNDC()  # coordinate normalizzate (0-1)
        latex.SetTextSize(0.04)
        latex.SetTextColor(ROOT.kRed)
        latex.DrawLatex(0.15, 0.85, f"RMS: {histo.GetRMS():0.4f}")#(in alto a sinistra)
        c2.Update()
        c2.SaveAs(outdir+"/summaryPlots/timeResolution/c_timeRes_spread_"+label+".png")
        c2.SaveAs(outdir+"/summaryPlots/timeResolution/c_timeRes_spread_"+label+".pdf")
        if(tRes_ave==0 ): continue
        else:
            print(f"Vov{vov} enBin{enBin} hits{hit}: {tRes_ave:.01f} [ps] spread (RMS) of {100*histo.GetRMS():0.2f}%")

    # --- print time resolution for central bars
    print("Time resolution for central bars(6,7,8): ")
    for bar in [6,7,8]:
        for vov in Vovs:
            for th in thresholds:
                for enBin in enBins:
                    for hit in Hits:
                        if (tRes_all[bar,vov,th,enBin,hit][0]==0): continue
                        else:
                            pass
                            #print(f" Bar {bar:02d} Vov{vov:.2f} th{th:02d} enBin{enBin} {hit}Hits: ({tRes_all[bar,vov,th,enBin,hit][0]:.2f} +/- {tRes_all[bar,vov,th,enBin,hit][1]:.2f}) [ps]")


# ---------- old selection
else:
    for vov in Vovs:
        for thr in thresholds:
            for enBin in enBins:
                g_deltaT_corr_vs_bars[vov,thr,enBin]=ROOT.TGraphErrors()
                g_deltaT_bestth_vs_bars[vov,enBin]=ROOT.TGraphErrors()
    # --- Read the histograms from moduleCharacterization_step2 file
    for label in label_list:
        print(label)
        inputFile == None
        inputFile = ROOT.TFile.Open(inputdir+'/moduleCharacterization_step2_%s.root'%label)
        for bar in bars:
            for vov in Vovs:
                for thr in thresholds:
                    tRes_corr={}
                    for enBin in enBins:
                            h1_deltaT_Corr = inputFile.Get('h1_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d'%(bar, vov, thr, enBin))
                            #print('h1_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d_%02dHits'%(bar, vov, thr, enBin,hit))
                            if not h1_deltaT_Corr:
                                print(f"Warning: histogram not found for bar{bar:02d}, Vov={vov:.2f}, thr={thr}, enBin={enBin}")
                                continue
                            tRes_corr[enBin]=getTimeResolution(h1_deltaT_Corr)
                            tRes, s_tRes= tRes_corr[enBin]
                            index=g_deltaT_corr_vs_bars[vov,thr,enBin].GetN()
                            g_deltaT_corr_vs_bars[vov,thr,enBin].SetPoint(index,bar,tRes/kscale)
                            g_deltaT_corr_vs_bars[vov,thr,enBin].SetPointError(index,0,s_tRes/kscale)
                            if(tRes_corr[enBin][0]<tRes_best[bar,vov,enBin][0]):
                                tRes_best[bar,vov,enBin]=tRes_corr[enBin]
        
        for bar in bars:
            for vov in Vovs:
                for enBin in enBins:
                        if(tRes_best[bar,vov,enBin][0]!=9999):
                            index=g_deltaT_bestth_vs_bars[vov,enBin].GetN()
                            g_deltaT_bestth_vs_bars[vov,enBin].SetPoint(index,bar,tRes_best[bar,vov,enBin][0]/kscale)
                            g_deltaT_bestth_vs_bars[vov,enBin].SetPointError(index,0,tRes_best[bar,vov,enBin][1]/kscale)

    print("Average time resolution for: ")
    for key, graph in g_deltaT_corr_vs_bars.items():
        vov, thr, enBin = key
        label=f"Vov{vov}_th{thr:02d}_enBin{enBin}"
        c=ROOT.TCanvas(f"c_timeRes_corr_Vov{vov}_thr{thr}_enBin{enBin}","",800,600)
        c.SetGrid()
        c.SetLeftMargin(0.15)
        c.SetBottomMargin(0.15)
        graph.SetTitle("")
        graph.GetXaxis().SetTitle("Bar")
        graph.GetYaxis().SetTitle("#sigma_{t} [ps]")
        graph.GetYaxis().SetRangeUser(0,140)
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(1.2)
        graph.SetLineWidth(2)
        graph.SetLineColor(ROOT.kBlue)
        graph.SetMarkerColor(ROOT.kBlue)

        graph.Draw("AP")
        fitRes=ROOT.TF1("fitRes",'pol0',1,13)
        graph.Fit(fitRes,"QRN")#quite range no-store
        tRes_ave=fitRes.GetParameter(0)
        #print(f"Vov{vov} th{thr} enBin{enBin}: {fitRes.GetParameter(0):.01f}[ps]")
        label=f"c_timeRes_corr_Vov{vov}_thr{thr:02d}_enBin{enBin}"
        c.SaveAs(outdir+"/summaryPlots/timeResolution/"+label+".png")
        c.SaveAs(outdir+"/summaryPlots/timeResolution/"+label+".pdf")

        c2=ROOT.TCanvas(f"c_tRes_spread_"+label,"",800,600)
        c2.SetGrid()
        c2.SetLeftMargin(0.15)
        c2.SetBottomMargin(0.15)
        histo=ROOT.TH1F("tRes_spread"+label,"",10,-0.5,0.5)
        for i in range(graph.GetN()):
            x = ctypes.c_double()
            tRes = ctypes.c_double()
            graph.GetPoint(i, x, tRes)
            tRes=tRes.value
            spread=(tRes-tRes_ave)/tRes
            histo.Fill(spread)
        histo.GetXaxis().SetTitle("#frac{#sigma_{t}-<#sigma_{t}>}{#sigma_{t}}")
        histo.Draw()
        latex = ROOT.TLatex()
        latex.SetNDC()  # coordinate normalizzate (0-1)
        latex.SetTextSize(0.04)
        latex.SetTextColor(ROOT.kRed)
        latex.DrawLatex(0.15, 0.85, f"RMS: {histo.GetRMS():0.4f}")#(in alto a sinistra)
        c2.Update()
        c2.SaveAs(outdir+"/summaryPlots/timeResolution/c_timeRes_spread_"+label+".png")
        c2.SaveAs(outdir+"/summaryPlots/timeResolution/c_timeRes_spread_"+label+".pdf")

        print(f"Vov{vov} th{thr} enBin{enBin}: {tRes_ave:.01f} [ps] spread (RMS) of {100*histo.GetRMS():0.2f}%")
    
    for vov in Vovs:
        for enBin in enBins:
                c01=ROOT.TCanvas(f"c_timeRes_corr_Vov{vov}_ALLth_enBin{enBin}","",800,600)
                c01.SetGrid()
                c01.SetLeftMargin(0.15)
                c01.SetBottomMargin(0.15)

                legend = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.SetTextSize(0.03)
                colors01= [ROOT.kRed, ROOT.kBlue, ROOT.kBlack]
                for i,thr in enumerate(thresholds):
                    key = (vov, thr, enBin)
                    graph01=g_deltaT_corr_vs_bars.get(key)
                    if not graph01:
                        continue
                    graph01.SetMarkerStyle(20)
                    graph01.SetMarkerSize(1.1)
                    graph01.SetLineWidth(2)
                    graph01.SetLineColor(colors01[i])
                    graph01.SetMarkerColor(colors01[i])

                    if i == 0:
                        graph01.Draw("AP")
                        graph01.GetXaxis().SetTitle("Bar")
                        graph01.GetYaxis().SetTitle("#sigma_{t} [ps]")
                        graph01.GetYaxis().SetRangeUser(0, 140)
                        graph01.SetTitle("")
                    else:
                        graph01.Draw("P")

                    legend.AddEntry(graph01, f"thr {thr}", "ep")
                
                legend.Draw()
                label = f"c_timeRes_corr_ALLthr_Vov{vov}_enBin{enBin}"
                c01.Update()
                c01.SaveAs(outdir + "/summaryPlots/timeResolution/" + label + ".png")
                c01.SaveAs(outdir + "/summaryPlots/timeResolution/" + label + ".pdf")
    
    print("Average time resolution at best th for: ")
    for key,graph in g_deltaT_bestth_vs_bars.items():
        vov, enBin = key
        leg = ROOT.TLegend(0.20, 0.90, 0.60, 0.70)
        c=ROOT.TCanvas(f"c_timeRes_bestTh_corr_Vov{vov}_enBin{enBin}","",800,600)
        c.SetGrid()
        c.SetLeftMargin(0.15)
        c.SetBottomMargin(0.15)
        graph.SetTitle("")
        graph.GetXaxis().SetTitle("Bar")
        graph.GetYaxis().SetTitle("#sigma_{t} [ps]")
        graph.GetYaxis().SetRangeUser(0,140)
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(1.2)
        graph.SetLineWidth(2)
        graph.SetLineColor(ROOT.kGreen)
        graph.SetMarkerColor(ROOT.kGreen)

        graph.Draw("AP")
        fitRes=ROOT.TF1("fitRes",'pol0',1,13)
        graph.Fit(fitRes,"QRN")#quite range no-store
        line = ROOT.TLine(1,tRes_ave ,13,tRes_ave)
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(2)  # linea tratteggiata
        line.Draw("same")
        label=f"c_timeRes_corr_Vov{vov}_best_th_enBin{enBin}"
        c.SaveAs(outdir+"/summaryPlots/timeResolution/"+label+".png")
        c.SaveAs(outdir+"/summaryPlots/timeResolution/"+label+".pdf")
        tRes_ave=fitRes.GetParameter(0)
        #print(f"Vov{vov} enBin{enBin}: {fitRes.GetParameter(0):.01f} [ps]")
        label=f"Vov{vov}_th{thr:02d}_enBin{enBin}"
        c2=ROOT.TCanvas(f"c_tRes_spread_"+label,"",800,600)
        c2.SetGrid()
        c2.SetLeftMargin(0.15)
        c2.SetBottomMargin(0.18)
        histo=ROOT.TH1F("tRes_spread"+label,"",10,-0.5,0.5)
        for i in range(graph.GetN()):
            x = ctypes.c_double()
            tRes = ctypes.c_double()
            graph.GetPoint(i, x, tRes)
            tRes=tRes.value
            spread=(tRes-tRes_ave)/tRes
            histo.Fill(spread)
        histo.GetXaxis().SetTitle("#frac{#sigma_{t}-<#sigma_{t}>}{#sigma_{t}}")
        histo.Draw()
        latex = ROOT.TLatex()
        latex.SetNDC()  # coordinate normalizzate (0-1)
        latex.SetTextSize(0.04)
        latex.SetTextColor(ROOT.kRed)
        latex.DrawLatex(0.15, 0.85, f"RMS: {histo.GetRMS():0.4f}")#(in alto a sinistra)
        c2.Update()
        c2.SaveAs(outdir+"/summaryPlots/timeResolution/c_timeRes_spread_"+label+".png")
        c2.SaveAs(outdir+"/summaryPlots/timeResolution/c_timeRes_spread_"+label+".pdf")

        print(f"Vov{vov} th{thr} enBin{enBin}: {tRes_ave:.01f} [ps] spread (RMS) of {100*histo.GetRMS():0.2f}%")






