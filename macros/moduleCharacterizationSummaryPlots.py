#! /usr/bin/env python
import os
import shutil
import glob
import math
import array
import sys
import time
import argparse

parser = argparse.ArgumentParser(description='Module characterization summary plots')
#parser.add_argument("-r",  "--runs",          required=True, type=str, help="comma-separated list of runs to be processed")
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
ROOT.gStyle.SetLabelSize(0.04,'X')
ROOT.gStyle.SetLabelSize(0.04,'Y')
ROOT.gStyle.SetTitleSize(0.04,'X')
ROOT.gStyle.SetTitleSize(0.04,'Y')
ROOT.gStyle.SetTitleOffset(1.1,'X')
ROOT.gStyle.SetTitleOffset(1.2,'Y')
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

source = 'TB'
tResMax = 150

# create files list
label_list = (args.inputLabels.split(','))
'''
print label_list
for run in label_list:
    if ('-' in run ): 
        all_runs = [i for i in range(int(run.split('-')[0]), int(run.split('-')[1]))]
        label_list.remove(run)
        label_list=label_list+all_runs
label_list = [int(r) for r in label_list]
label_list.sort()
'''
print label_list

# resolution mode : 0 : /1;  1: /sqrt(2) if CTR, 2: /2 if TDiff
kscale = 2.
if (args.resMode == 0): kscale = 1
if (args.resMode == 1): kscale = math.sqrt(2)

# output
outdir = '/var/www/html/TOFHIR2X/MTDTB_CERN_Oct21/ModuleCharacterization/'+args.outFolder
print 'Saving plots in ', outdir
#outFileName = '/home/cmsdaq/Lab5015Analysis_new/TB_CERN_Oct21/Lab5015Analysis/plots/'+args.outFolder+'.root'
outFileName = '/home/cmsdaq/Lab5015Analysis_new/martina_TB_CERN_Oct21/Lab5015Analysis/plots/'+args.outFolder+'.root'
outfile = ROOT.TFile(outFileName, 'RECREATE' )



peaks = []
enBins = []
refPeak = 0
if (source == 'Na22'):
    peaks = [511, 1275]
    enBins = [5, 7] ##??? non me li ricordo
    refPeak = 511
if (source == 'Laser'):
    peaks = [0]
    enBins = [1]
if (source == 'TB'):
    peaks = [0]
    enBins = [1]


# --- colors
cols = { 1.00 : 49,  
         1.25 : 50,  
         1.40 : 50, 
         1.50 : 51, 
         1.75 : 51+5,
         1.70 : 51+5,
         1.80 : 51+5,
         2.00 : 51 + 10, 
         2.10 : 51 + 12, 
         2.25 : 51 + 15,
         2.50  : 51 + 20,
         2.80  : 51 + 22,
         3.00  : 51 + 25,
         3.50  : 51 + 40,
         3.70  : 51 + 40,
         4.00  : 51 + 45,
         5.00  : 51 + 50,
         3.83  : 4}


# --- prepare output dir
if (os.path.isdir(outdir) == False): 
    os.system('mkdir %s'%outdir)
os.system('mkdir %s/summaryPlots/'%outdir)
os.system('mkdir %s/summaryPlots/tot/'%outdir)
os.system('mkdir %s/summaryPlots/energy/'%outdir)
os.system('mkdir %s/summaryPlots/timeResolution/'%outdir)
    

# -- ref threhsold
thRef = 10

# -- get list of bars, Vovs, thresholds to be analyzed
bars = []
thresholds = []
Vovs = [] 
for label in label_list:
    inputFile = ROOT.TFile.Open('/home/cmsdaq/Lab5015Analysis_new/martina_TB_CERN_Oct21/Lab5015Analysis/plots/moduleCharacterization_step2_%s.root'%label)
    listOfKeys = [key.GetName().replace('h1_deltaT_energyRatioCorr_','') for key in ROOT.gDirectory.GetListOfKeys() if key.GetName().startswith('h1_deltaT_energyRatioCorr')]
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

bars.remove(14)
bars.remove(15)

print 'bars:', bars
print 'Vovs:',Vovs
print 'thresholds:', thresholds

#if ('HPK528' in args.outFolder):
#bars = [4,5,6,7,8,9,10,11,12]

#if ('FBK' in args.outFolder or '4355' in args.outFolder):
#    bars = [5,6,7,8,9,10,11,12]
#
#if ('thick' in args.outFolder):
#    bars = [4,5,6,7,9,10,11,12] 

#raw_input('continue?')

# --- Summary graphs
g_tot_vs_th  = {} # g [bar, l, vov] 
g_tot_vs_vov  = {} # g [bar, l, th] 
g_tot_vs_bar = {} # g [l, vov] 

g_energy_vs_th  = {} # g [bar, l, vov, peak] 
g_energy_vs_vov  = {} # g [bar, l, th, peak] 
g_energy_vs_bar = {} # g [l, vov, peak] 

g_deltaT_energyRatioCorr_vs_th  = {} # g [bar, vov, energyBin] 
g_deltaT_energyRatioCorr_vs_vov  = {} # g [bar, th, energyBin] 
g_deltaT_energyRatioCorr_vs_bar = {} # g [viv, thr, energyBin] 
g_deltaT_energyRatioCorr_bestTh_vs_vov = {} # g [bar, energyBin] 
g_deltaT_energyRatioCorr_bestTh_vs_bar = {} # g [vov, energyBin] 

for bar in bars:
    for l in ['L','R','L-R']:
        for vov in Vovs:
            g_tot_vs_th[bar, l, vov] = ROOT.TGraphErrors()
            for peak in peaks: g_energy_vs_th[bar, l, vov, peak] = ROOT.TGraphErrors()
            for enBin in enBins: g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin] = ROOT.TGraphErrors()
        for thr in thresholds:
            g_tot_vs_vov[bar, l, thr] = ROOT.TGraphErrors()
            for peak in peaks: g_energy_vs_vov[bar, l, thr, peak] = ROOT.TGraphErrors()
            for enBin in enBins: g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin] = ROOT.TGraphErrors()
                        
for l in ['L','R','L-R']:
    for vov in Vovs:
        for thr in thresholds:
            g_tot_vs_bar[l, vov, thr] = ROOT.TGraphErrors()
            for peak in peaks: g_energy_vs_bar[l, vov, thr, peak] = ROOT.TGraphErrors()
            if (l=='L-R'):
                for enBin in enBins: 
                    g_deltaT_energyRatioCorr_vs_bar[vov, thr, enBin] = ROOT.TGraphErrors()

for bar in bars:
    for enBin in enBins: 
        g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin] = ROOT.TGraphErrors()

for vov in Vovs:
    for enBin in enBins: 
        g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin] = ROOT.TGraphErrors()


# --- Read the histograms from moduleCharacterization_step2 file
for label in label_list:
    inputFile = ROOT.TFile.Open('/home/cmsdaq/Lab5015Analysis_new/martina_TB_CERN_Oct21/Lab5015Analysis/plots/moduleCharacterization_step2_%s.root'%label)

    for bar in bars:
        for l in ['L','R','L-R']:
            for vov in Vovs:
                for thr in thresholds:
                    # -- tot vs thr, Vov, bar
                    if ( inputFile.GetListOfKeys().Contains('h1_tot_bar%02d%s_Vov%.02f_th%02d'%(bar, l, vov, thr)) ): 
                        h1_tot = inputFile.Get('h1_tot_bar%02d%s_Vov%.02f_th%02d'%(bar, l, vov, thr))
                        if (h1_tot==None): continue
                        
                        max1 = h1_tot.GetBinCenter(h1_tot.GetMaximumBin())
                        h1_tot.GetXaxis().SetRangeUser(0.25*max1,2.*max1)
                        fitFunc = ROOT.TF1('fitFunc','gaus',max1-0.05*max1,max1+0.05*max1)
                        h1_tot.Fit('fitFunc','QNRS+')
                        h1_tot.Fit('fitFunc','QNS+','', fitFunc.GetParameter(1)-fitFunc.GetParameter(2), fitFunc.GetParameter(1)+fitFunc.GetParameter(2))
                        
                        g_tot_vs_th[bar, l, vov].SetPoint(g_tot_vs_th[bar, l, vov].GetN(), thr, fitFunc.GetParameter(1) )    
                        g_tot_vs_th[bar, l, vov].SetPointError(g_tot_vs_th[bar, l, vov].GetN()-1, 0, fitFunc.GetParError(1) )    
                        
                        g_tot_vs_vov[bar, l, thr].SetPoint(g_tot_vs_vov[bar, l, thr].GetN(), vov, fitFunc.GetParameter(1) )    
                        g_tot_vs_vov[bar, l, thr].SetPointError(g_tot_vs_vov[bar, l, thr].GetN()-1, 0, fitFunc.GetParError(1) )    
                        
                        if (thr == thRef):
                            g_tot_vs_bar[l, vov, thr].SetPoint(g_tot_vs_bar[l, vov, thr].GetN(), bar, fitFunc.GetParameter(1) )    
                            g_tot_vs_bar[l, vov, thr].SetPointError(g_tot_vs_bar[l, vov, thr].GetN()-1, 0, fitFunc.GetParError(1) )    
                    
                    # -- energy vs thr, Vov, bar
                    if ( inputFile.GetListOfKeys().Contains('h1_energy_bar%02d%s_Vov%.02f_th%02d'%(bar, l, vov, thr)) ): 
                        h1_energy = inputFile.Get('h1_energy_bar%02d%s_Vov%.02f_th%02d'%(bar, l, vov, thr))
                        if ( h1_energy == None): continue
                        energyPeak = {}
                        if (source == 'Na22'):
                            for peak in peaks:
                                energyPeak[peak] = [h1_energy.GetFunction('fitFunc_%d'%peak).GetParameter(1), h1_energy.GetFunction('fitFunc_%d'%peak).GetParError(1)]
                        elif ( source == 'Laser'):
                            fitFunc = ROOT.TF1('fitFunc','gaus',0, 1000)
                            h1_energy.Fit('fitFunc','QNRS+')
                            h1_energy.Fit('fitFunc','QNS+','', fitFunc.GetParameter(1)-fitFunc.GetParameter(2), fitFunc.GetParameter(1)+fitFunc.GetParameter(2))
                            for peak in peaks:
                                energyPeak[peak] = [ fitFunc.GetParameter(1), fitFunc.GetParError(1)]
                        elif ( source == 'TB'):
                            #fitFunc = h1_energy.GetFunction('f_landau_bar%02d%s_Vov%.02f_vth1_%02d'%(bar, l, vov, thr))
                            #if (fitFunc==None): continue
                            fitFunc = ROOT.TF1('f_landau_bar%02d%s_Vov%.02f_vth1_%02d'%(bar, l, vov, thr), '[0]*TMath::Landau(x,[1],[2])', 0,1000.)
                            h1_energy.GetXaxis().SetRangeUser(50,800)
                            emax = h1_energy.GetBinCenter(h1_energy.GetMaximumBin())
                            fitFunc.SetParameters(10, emax, 30)
                            fitFunc.SetRange(0.8*emax, 1.5*emax)
                            h1_energy.Fit(fitFunc,'QR')
                            print bar, vov, thr, emax, fitFunc.GetParameter(1), fitFunc.GetParError(1)
                            for peak in peaks:
                                energyPeak[peak] = [fitFunc.GetParameter(1), fitFunc.GetParError(1)] 
                                

                        for peak in peaks:
                            g_energy_vs_th[bar, l, vov, peak].SetPoint(g_energy_vs_th[bar, l, vov, peak].GetN(), thr, energyPeak[peak][0] )
                            g_energy_vs_th[bar, l, vov, peak].SetPointError(g_energy_vs_th[bar, l, vov, peak].GetN()-1, 0, energyPeak[peak][1])
                            
                            g_energy_vs_vov[bar, l, thr, peak].SetPoint(g_energy_vs_vov[bar, l, thr, peak].GetN(), vov, energyPeak[peak][0] )
                            g_energy_vs_vov[bar, l, thr, peak].SetPointError(g_energy_vs_vov[bar, l, thr, peak].GetN()-1, 0, energyPeak[peak][1] )
                    
                            if (thr == thRef):
                                g_energy_vs_bar[l, vov, thr, peak].SetPoint(g_energy_vs_bar[l, vov, thr, peak].GetN(), bar, energyPeak[peak][0] )
                                g_energy_vs_bar[l, vov, thr, peak].SetPointError(g_energy_vs_bar[l, vov, thr, peak].GetN()-1, 0, energyPeak[peak][1] )
                        
    # -- tRes vs Vov, thr, bar
    bestRes = {}
    bestTh = {}

    for bar in bars:
        for vov in Vovs:
            for enBin in enBins:    
                bestRes[bar, vov, enBin] = [9999, 9999]
                bestTh[bar, vov, enBin]  = 10 
    
    for bar in bars:
        for vov in Vovs:
            for thr in thresholds: 
                tRes = {}
                for enBin in enBins:
                    #h1_deltaT = inputFile.Get('h1_deltaT_energyRatioCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d'%(bar, vov, thr, enBin))
                    #h1_deltaT = inputFile.Get('h1_deltaT_energyRatioPhaseCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d'%(bar, vov, thr, enBin))
                    h1_deltaT = inputFile.Get('h1_deltaT_totRatioPhaseCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d'%(bar, vov, thr, enBin))
                    if (h1_deltaT == None): continue
                    if (h1_deltaT.GetEntries() < 10 ): continue

                    h1_deltaT.GetXaxis().SetRangeUser(h1_deltaT.GetMean() - 5*h1_deltaT.GetRMS(), h1_deltaT.GetMean() + 5*h1_deltaT.GetRMS())

                    
                    fitFunc = ROOT.TF1('fitFunc','gaus',-10000, 10000)
                    fitFunc.SetLineColor(ROOT.kGreen+3)
                    fitFunc.SetLineWidth(2)
                    fitFunc.SetParameters(h1_deltaT.GetMaximum(),h1_deltaT.GetMean(), h1_deltaT.GetRMS())
                    #fitXMin = h1_deltaT.GetBinCenter(h1_deltaT.GetMaximumBin()) - 200
                    #fitXMax = h1_deltaT.GetBinCenter(h1_deltaT.GetMaximumBin()) + 200.

                    fitXMin = h1_deltaT.GetMean() - 3*h1_deltaT.GetRMS()
                    fitXMax = h1_deltaT.GetMean() + 3*h1_deltaT.GetRMS()
                    fitFunc.SetRange(fitXMin, fitXMax)
                    h1_deltaT.Fit('fitFunc','QNRS+','', fitXMin, fitXMax)
                    fitFunc.SetRange(fitFunc.GetParameter(1) - 3.0*fitFunc.GetParameter(2), fitFunc.GetParameter(1) + 3.0*fitFunc.GetParameter(2))
                    h1_deltaT.Fit('fitFunc','QNRS+')
                    fitFunc.SetRange(fitFunc.GetParameter(1) - 2.5*fitFunc.GetParameter(2), fitFunc.GetParameter(1) + 2.5*fitFunc.GetParameter(2))
                    h1_deltaT.Fit('fitFunc','QRS+')

                    
                    if (fitFunc==None): continue
                    if (fitFunc.GetParameter(2) > 1000): continue
                    if (fitFunc.GetParameter(2) < 20): continue
                    if (fitFunc.GetParError(2) > 200): continue
                    tRes[enBin] = [ fitFunc.GetParameter(2),fitFunc.GetParError(2)]
                    print bar, vov, thr, tRes[enBin]
                    #print 'best res, par2 , err2 : ', (bestRes[bar, vov, enBin], fitFunc.GetParameter(2), fitFunc.GetParError(2))
                    if (fitFunc.GetParameter(2) < bestRes[bar, vov, enBin][0]):
                        bestRes[bar, vov, enBin] = [fitFunc.GetParameter(2),fitFunc.GetParError(2)]
                        bestTh[bar, vov, enBin]  = thr
                        
                    ctemp = ROOT.TCanvas()
                    h1_deltaT.GetYaxis().SetRangeUser(0, h1_deltaT.GetBinContent(h1_deltaT.GetMaximumBin())*1.2)                
                    h1_deltaT.Draw()                
                    ctemp.SaveAs(outdir+'/summaryPlots/timeResolution/'+'/c_h1_deltaT_energyRatioCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d.png'%(bar, vov, thr, enBin))

                    g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin].SetPoint(g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin].GetN(), thr, tRes[enBin][0]/kscale )
                    g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin].SetPointError(g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin].GetN()-1, 0, tRes[enBin][1]/kscale)
                        
                    g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin].SetPoint(g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin].GetN(), vov, tRes[enBin][0]/kscale )
                    g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin].SetPointError(g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin].GetN()-1, 0, tRes[enBin][1]/kscale)
                        
                    if (thr == thRef):
                        g_deltaT_energyRatioCorr_vs_bar[vov, thr, enBin].SetPoint(g_deltaT_energyRatioCorr_vs_bar[vov, thr, enBin].GetN(), bar, tRes[enBin][0]/kscale )
                        g_deltaT_energyRatioCorr_vs_bar[vov, thr, enBin].SetPointError(g_deltaT_energyRatioCorr_vs_bar[vov, thr, enBin].GetN()-1, 0, tRes[enBin][1]/kscale)
                
            for enBin in enBins:
                if (bestRes[bar, vov, enBin][0]==9999): continue

                #print bar, vov, enBin, bestRes[bar, vov, enBin]
                # -- tRes  vs Vov at the best threshold    
                g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].SetPoint(g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].GetN(), vov, (bestRes[bar, vov, enBin][0])/kscale )
                g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].SetPointError(g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].GetN()-1, 0, (bestRes[bar, vov, enBin][1])/kscale)

                # -- tRes  vs bar at the best threshold
                g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin].SetPoint(g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin].GetN(), bar, (bestRes[bar, vov, enBin][0])/kscale )
                g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin].SetPointError(g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin].GetN()-1, 0, (bestRes[bar, vov, enBin][1])/kscale)            


# -- Draw
for bar in bars:
    
    # -- tot vs threshold
    ctot1 = ROOT.TCanvas('c_tot_vs_th_bar%.02d'%bar)
    hPad1 = ROOT.TH2F('hPad1','', 100, -1., 64.,40, 0.,40.)
    hPad1.SetTitle(";threshold [DAC];ToT [ns]")
    hPad1.Draw()
    ctot1.SetGridy()
    leg = ROOT.TLegend(0.70, 0.70, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for i, vov in enumerate(Vovs):
        for l in ['L','R']:
            g_tot_vs_th[bar, l, vov].SetMarkerStyle(20)
            if (l == 'R'): g_tot_vs_th[bar, l, vov].SetMarkerStyle(24)
            g_tot_vs_th[bar, l, vov].SetMarkerColor(i+1)
            g_tot_vs_th[bar, l, vov].SetLineColor(i+1)
            g_tot_vs_th[bar, l, vov].Draw('plsame')
        leg.AddEntry(g_tot_vs_th[bar,'L', vov], 'V_{OV} = %.02f'%vov, 'PL')
    leg.Draw()    
    ctot1.SaveAs(outdir+'/summaryPlots/tot/'+ctot1.GetName()+'.png')
    ctot1.SaveAs(outdir+'/summaryPlots/tot/'+ctot1.GetName()+'.pdf')
    hPad1.Delete()

    # -- tot vs Vov
    ctot2 = ROOT.TCanvas('c_tot_vs_Vov_bar%.02d'%bar)
    hPad2 = ROOT.TH2F('hPad2','', 10, 0.,8.,40, 0.,40.)
    hPad2.SetTitle(";V_{OV} [V];ToT [ns]")
    hPad2.Draw()
    ctot2.SetGridy()
    leg = ROOT.TLegend(0.70, 0.70, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for i, thr in enumerate(thresholds):
        for l in ['L','R']:
            g_tot_vs_vov[bar,l, thr].SetMarkerStyle(20)
            if (l == 'R'): g_tot_vs_vov[bar, l, thr].SetMarkerStyle(24)
            g_tot_vs_vov[bar,l, thr].SetMarkerColor(i+1)
            g_tot_vs_vov[bar,l, thr].SetLineColor(i+1)
            g_tot_vs_vov[bar,l, thr].Draw('plsame')
        leg.AddEntry(g_tot_vs_vov[bar,'L', thr], 'th. = %d'%thr, 'PL')
    leg.Draw()    
    ctot2.SaveAs(outdir+'/summaryPlots/tot/'+ctot2.GetName()+'.png')
    ctot2.SaveAs(outdir+'/summaryPlots/tot/'+ctot2.GetName()+'.pdf')
    hPad2.Delete()
    
    # -- energy vs threshold
    cen1 = ROOT.TCanvas('c_energy_vs_th_bar%.02d'%bar)
    hPadEn1 = ROOT.TH2F('hPadEn1','', 100, -1., 64.,500, 0.,1000.)
    hPadEn1.SetTitle(";threshold [DAC]; energy")
    hPadEn1.Draw()
    cen1.SetGridy()
    leg = ROOT.TLegend(0.70, 0.50, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for i, vov in enumerate(Vovs):
        for l in ['L','R']:
            g_energy_vs_th[bar, l, vov, refPeak].SetMarkerStyle(20)
            if (l == 'R'): g_energy_vs_th[bar, l, vov, refPeak].SetMarkerStyle(24)
            g_energy_vs_th[bar, l, vov, refPeak].SetMarkerColor(i+1)
            g_energy_vs_th[bar, l, vov, refPeak].SetLineColor(i+1)
            g_energy_vs_th[bar, l, vov, refPeak].Draw('plsame')
        leg.AddEntry(g_energy_vs_th[bar,'L', vov, refPeak], 'V_{OV} = %.02f'%vov, 'PL')
    leg.Draw()
    cen1.SaveAs(outdir+'/summaryPlots/energy/'+cen1.GetName()+'.png')
    cen1.SaveAs(outdir+'/summaryPlots/energy/'+cen1.GetName()+'.pdf')
    hPadEn1.Delete()

    # -- energy vs Vov
    cen2 = ROOT.TCanvas('c_energy_vs_Vov_bar%.02d'%bar)
    hPadEn2 = ROOT.TH2F('hPadEn2','', 10, 0.,8.,50, 0.,1000.)
    hPadEn2.SetTitle(";V_{OV} [V]; energy")
    hPadEn2.Draw()
    cen2.SetGridy()
    leg = ROOT.TLegend(0.70, 0.50, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for i, thr in enumerate(thresholds):
        for l in ['L','R']:
            g_energy_vs_vov[bar,l, thr, refPeak].SetMarkerStyle(20)
            if (l == 'R'): g_energy_vs_vov[bar, l, thr, refPeak].SetMarkerStyle(24)
            g_energy_vs_vov[bar,l, thr, refPeak].SetMarkerColor(i+1)
            g_energy_vs_vov[bar,l, thr, refPeak].SetLineColor(i+1)
            g_energy_vs_vov[bar,l, thr, refPeak].Draw('plsame')
        leg.AddEntry(g_energy_vs_vov[bar,'L', thr, refPeak], 'th. = %d'%thr, 'PL')
    leg.Draw()
    cen2.SaveAs(outdir+'/summaryPlots/energy/'+cen2.GetName()+'.png')
    cen2.SaveAs(outdir+'/summaryPlots/energy/'+cen2.GetName()+'.pdf') 
    hPadEn2.Delete()   


    # -- time resolution vs threshold
    for enBin in enBins:
        ctres1 = ROOT.TCanvas('c_tRes_energyRatioCorr_vs_th_bar%.02d_enBin%02d'%(bar,enBin))
        #hPad = ROOT.gPad()
        #hPad.DrawFrame(-1.,0.,64.,500.)
        #hPad.SetTitle(";threshold [DAC];Energy [ns]")
        #hPad.Draw()
        #hPad.SetGridy()
        hPadT1 = ROOT.TH2F('hPadT1','', 100, -1., 64.,100, 0.,tResMax)
        hPadT1.SetTitle(";threshold [DAC]; #sigma_{t} [ns]")
        hPadT1.Draw()
        ctres1.SetGridy()
        leg = ROOT.TLegend(0.70, 0.50, 0.89, 0.89)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        for i, vov in enumerate(Vovs):
            #print 'Vov, entries = ', vov, g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin].GetN()
            g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin].SetMarkerStyle(20)
            g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin].SetMarkerColor(cols[vov])
            g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin].SetLineColor(cols[vov])
            g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin].Draw('plsame')
            leg.AddEntry(g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin], 'V_{OV} = %.02f'%vov, 'PL')
            outfile.cd()
            g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin].Write('g_deltaT_energyRatioCorr_vs_th_bar%02d_Vov%.02f_enBin%02d'%(bar, vov, enBin))
        leg.Draw()
        ctres1.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres1.GetName()+'.png')
        ctres1.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres1.GetName()+'.pdf')
        hPadT1.Delete()   

        # -- time resolution vs Vov
        ctres2 = ROOT.TCanvas('c_tRes_energyRatioCorr_vs_Vov_bar%.02d_enBin%02d'%(bar,enBin))
        hPadT2 = ROOT.TH2F('hPadT2','', 10, 0.,10.,10, 0.,tResMax)
        hPadT2.SetTitle(";V_{OV} [V];#sigma_{t} [ns]")
        hPadT2.Draw()
        ctres2.SetGridy()
        leg = ROOT.TLegend(0.70, 0.50, 0.89, 0.89)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        for i, thr in enumerate(thresholds):
            g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin].SetMarkerStyle(20)
            g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin].SetMarkerColor(i+1)
            g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin].SetLineColor(i+1)
            g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin].Draw('plsame')
            leg.AddEntry(g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin], 'th. = %d'%thr, 'PL')
        leg.Draw()
        ctres2.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres2.GetName()+'.png')
        ctres2.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres2.GetName()+'.pdf')
        hPadT2.Delete()   

        # -- time resolution vs Vov at the best Th
        ctres2 = ROOT.TCanvas('c_tRes_energyRatioCorr_bestTh_vs_Vov_bar%.02d_enBin%02d'%(bar,enBin))
        hPadT2 = ROOT.TH2F('hPadT2','', 10, 0.,10.,10, 0.,tResMax)
        hPadT2.SetTitle(";V_{OV} [V];#sigma_{t} [ns]")
        hPadT2.Draw()
        ctres2.SetGridy()
        g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].SetMarkerStyle(20)
        g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].SetMarkerColor(i+1)
        g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].SetLineColor(i+1)
        g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].Draw('plsame')
        outfile.cd() 
        g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].Write('g_deltaT_energyRatioCorr_bestTh_vs_vov_bar%02d_enBin%02d'%(bar, enBin))
        ctres2.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres2.GetName()+'.png')
        ctres2.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres2.GetName()+'.pdf')
        hPadT2.Delete()   

# -- plots vs bar

# -- tot vs bar 
for i, vov in enumerate(Vovs):
    ctot3 = ROOT.TCanvas('c_tot_vs_bar_Vov%.02f'%vov)
    hPad3 = ROOT.TH2F('hPad3','', 100, -0.5, 15.5,40, 0.,40.)
    hPad3.SetTitle("; bar; ToT [ns]")
    hPad3.Draw()
    ctot3.SetGridy()
    leg = ROOT.TLegend(0.70, 0.50, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for i, vov in enumerate(Vovs):
        for l in ['L','R']:
            g_tot_vs_bar[l, vov, thRef].SetMarkerStyle(20)
            if (l == 'R'): g_tot_vs_bar[l, vov, thRef].SetMarkerStyle(24)
            g_tot_vs_bar[l, vov, thRef].SetMarkerColor(cols[vov])
            g_tot_vs_bar[l, vov, thRef].SetLineColor(cols[vov])
            g_tot_vs_bar[l, vov, thRef].Draw('plsame')
        leg.AddEntry(g_tot_vs_bar['L', vov, thRef], 'V_{OV} = %.02f'%vov, 'PL')
    leg.Draw()
    ctot3.SaveAs(outdir+'/summaryPlots/tot/'+ctot3.GetName()+'.png')
    ctot3.SaveAs(outdir+'/summaryPlots/tot/'+ctot3.GetName()+'.pdf')    
    hPad3.Delete()

# -- energy vs bar
for i, vov in enumerate(Vovs):
    cen3 = ROOT.TCanvas('c_energy_vs_bar_Vov%.02f'%vov)
    #hPadEn3 = ROOT.TH2F('hPadEn3','', 100, -0.5, 15.5, 40, 0.,1000.)
    hPadEn3 = ROOT.TH2F('hPadEn3','', 100, -0.5, 15.5, 40, 0.,200.)
    hPadEn3.SetTitle("; bar; energy")
    hPadEn3.Draw()
    cen3.SetGridy()
    leg = ROOT.TLegend(0.70, 0.50, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for i, vov in enumerate(Vovs):
        for l in ['L','R']:
            g_energy_vs_bar[l, vov, thRef, refPeak].SetMarkerStyle(20)
            if (l == 'R'): g_energy_vs_bar[l, vov, thRef, refPeak].SetMarkerStyle(24)
            g_energy_vs_bar[l, vov, thRef, refPeak].SetMarkerColor(cols[vov])
            g_energy_vs_bar[l, vov, thRef, refPeak].SetLineColor(cols[vov])
            g_energy_vs_bar[l, vov, thRef, refPeak].Draw('plsame')
            print 'Vov = %.02f, average  energy %s = %.02f'%(vov, l, g_energy_vs_bar[l, vov, thRef, refPeak].GetMean(2))
            outfile.cd()  
            g_energy_vs_bar[l, vov, thRef, refPeak].Write('g_energy%s_vs_bar_Vov%.02f_th%02d'%(l,vov,thRef))
        leg.AddEntry(g_energy_vs_bar['L', vov, thRef, refPeak], 'V_{OV} = %.02f'%vov, 'PL')
    leg.Draw()
    cen3.SaveAs(outdir+'/summaryPlots/energy/'+cen3.GetName()+'.png')
    cen3.SaveAs(outdir+'/summaryPlots/energy/'+cen3.GetName()+'.pdf')    
    hPadEn3.Delete()

# -- time resolution vs bar at the ref threshold
for enBin in enBins:
    ctres3 = ROOT.TCanvas('c_tRes_energyRatioCorr_vs_bar_Vov%.02f_enBin%02d'%(vov,enBin))
    hPadT3 = ROOT.TH2F('hPadT3','', 100, -0.5, 15.5,100, 0.,tResMax)
    hPadT3.SetTitle("; bar; #sigma_{t}[ns]")
    hPadT3.Draw()
    ctres3.SetGridy()
    leg = ROOT.TLegend(0.70, 0.50, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for i, vov in enumerate(Vovs):
        g_deltaT_energyRatioCorr_vs_bar[vov, thRef, enBin].SetMarkerStyle(20)
        g_deltaT_energyRatioCorr_vs_bar[vov, thRef, enBin].SetMarkerColor(cols[vov])
        g_deltaT_energyRatioCorr_vs_bar[vov, thRef, enBin].SetLineColor(cols[vov])
        g_deltaT_energyRatioCorr_vs_bar[vov, thRef, enBin].Draw('plsame')
        outfile.cd() 
        g_deltaT_energyRatioCorr_vs_bar[vov, thRef, enBin].Write('g_deltaT_energyRatioCorr_vs_bar__Vov%.02f_th%02d'%(vov,thRef))
        leg.AddEntry(g_deltaT_energyRatioCorr_vs_bar[vov, thRef, enBin], 'V_{OV} = %.02f'%vov, 'PL')
    leg.Draw()
    ctres3.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres3.GetName()+'.png')
    ctres3.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres3.GetName()+'.pdf')    
    hPadT3.Delete()  

# -- time resolution vs bar at the best threshold
for enBin in enBins:
    ctres3 = ROOT.TCanvas('c_tRes_energyRatioCorr_bestTh_vs_bar_enBin%02d'%(enBin))
    hPadT3 = ROOT.TH2F('hPadT3','', 100, -0.5, 15.5,100, 0.,tResMax)
    hPadT3.SetTitle("; bar; #sigma_{t} [ns]")
    hPadT3.Draw()
    ctres3.SetGridy()
    leg = ROOT.TLegend(0.70, 0.50, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for i, vov in enumerate(Vovs):
        #fitRes = ROOT.TF1('fitRes','pol0',0,16)
        #g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin].Fit(fitRes,'QRS')
        #if (g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin].GetMean(2) > 0  ):
        g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin].SetMarkerStyle(20)
        g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin].SetMarkerColor(cols[vov])
        g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin].SetLineColor(cols[vov])
        g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin].Draw('plsame')
        print 'Vov = %0.02f --> Spread of tRes = %.03f'%(vov, g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin].GetRMS(2)/g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin].GetMean(2))
        leg.AddEntry(g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin], 'V_{OV} = %.02f'%vov, 'PL')
        outfile.cd()
        g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin].Write('g_deltaT_energyRatioCorr_bestTh_vs_bar_Vov%.02f_enBin%02d'%(vov, enBin))
    leg.Draw()
    ctres3.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres3.GetName()+'.png')
    ctres3.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres3.GetName()+'.pdf')    
    hPadT3.Delete()  

outfile.Close()

raw_input('OK?')

