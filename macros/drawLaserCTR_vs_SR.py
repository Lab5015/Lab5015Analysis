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

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetTitleOffset(1.05,'Y')
ROOT.gStyle.SetLabelSize(0.04)
ROOT.gErrorIgnoreLevel = ROOT.kWarning;
ROOT.gROOT.SetBatch(True)
#ROOT.gROOT.SetBatch(False)


plotDir = '/var/www/html/TOFHIR2X/MTDST_CERN_Oct21/'
outFile = ROOT.TFile("/home/petsys/TOFHiR2X/sw_analysis/Lab5015Analysis/plots_fede/drawLaserCTR_vs_SR.root","RECREATE");


dac_to_uA = { 
    'ith2_3': 1.250,
    'ith2_2': 0.940,
    'ith2_1': 0.630,
    'ith2_0': 0.313,
    'ith1_3': 0.630,
    'ith1_2': 0.470,
    'ith1_1': 0.313,
    'ith1_0': 0.156
}

VovList = [1.5, 2.5, 3.5]

runs_dict = { 

    84 : [ 'ith2_0', 0.60, [2.5, 3.5]],
    87 : [ 'ith2_0', 0.50, [2.5, 3.5]],
    89 : [ 'ith2_0', 0.40, [2.5, 3.5]],
    91 : [ 'ith2_0', 0.30, [2.5, 3.5]],
    93 : [ 'ith2_0', 0.20, [2.5, 3.5]],

    98 : [ 'ith2_3', 0.80, [2.5, 3.5]],
    99 : [ 'ith2_3', 0.70, [1.5, 2.5, 3.5]],
    100 : [ 'ith2_3', 0.60, [1.5, 2.5, 3.5]],
    101 : [ 'ith2_3', 0.50, [1.5, 2.5, 3.5]],
    102 : [ 'ith2_3', 0.40, [1.5, 2.5, 3.5]],
    103 : [ 'ith2_3', 0.30, [1.5, 2.5, 3.5]],
    104 : [ 'ith2_3', 0.20, [1.5, 2.5, 3.5]],
    105 : [ 'ith2_3', 0.10, [1.5, 2.5, 3.5]]

}

g_tRes_vs_SR = {}
g_tRes_vs_SR_all = ROOT.TGraphErrors()
for Vov in VovList:
    g_tRes_vs_SR[Vov] = ROOT.TGraphErrors()


#---------------------------
# -- compute SR vs threshold
npoints = 7

for run in sorted(runs_dict):
    print("===>>> Processing run: ", run)

    ithMode = runs_dict[run][0]
    laserTune = runs_dict[run][1]
    Vovs  = runs_dict[run][2]

    for Vov in Vovs:
        inFile = ROOT.TFile.Open('/home/petsys/TOFHiR2X/sw_analysis/Lab5015Analysis/plots_fede/pulseShape_run%04d.root'%run)
    
        SR = 0.
        thresh = 0.
        channels = ['ch1', 'ch2']

        graph1 = inFile.Get('g_ps_totSel_%s_Vov%.01f'%('ch1',Vov))
        graph2 = inFile.Get('g_ps_totSel_%s_Vov%.01f'%('ch2',Vov))
        if graph1 == None or graph2 == None:
            continue
        
        index_cen = int(min(graph1.GetN()/4,graph2.GetN()/4))
        index_min = max(0,index_cen-int(npoints/2))
        index_max = min(index_cen+int(npoints/2),int(min(graph1.GetN()/2,graph2.GetN()/2)))

        for ch in channels:
            graph = inFile.Get('g_ps_totSel_%s_Vov%.01f'%(ch,Vov))
            fitSR = ROOT.TF1('fitSR', 'pol1', -10., 30.)
            fitSR.SetRange( graph.GetPointX(index_min)-0.001, graph.GetPointX(index_max)+0.001 )
            fitSR.SetParameters(0,20)
            graph.Fit(fitSR,'QRS')
            outFile.cd()
            graph.Write('Run%s_g_ps_totSel_%s_Vov%.01f'%(run,ch,Vov))
        
            #print index_cen, index_min, index_max, graph.GetY()[index_cen], fitSR.GetParameter(1)
            SR += fitSR.GetParameter(1)
            thresh += int(round(graph.GetPointY(index_cen)/dac_to_uA[ithMode]))
        
        SR /= 2.
        thresh /= 2.
    
        print(SR,thresh)
    
        inFile = ROOT.TFile.Open('/home/petsys/TOFHiR2X/sw_analysis/Lab5015Analysis/plots_fede/moduleCharacterization_step2_run%.04d.root'%run)
        
        print('NAME: h1_deltaT_totRatioCorr_bar02L-R_Vov%.02f_th%.02d_energyBin01'%(Vov,thresh))
        histo = inFile.Get('h1_deltaT_totRatioCorr_bar02L-R_Vov%.02f_th%.02d_energyBin01'%(Vov,thresh))
        if histo == None:
            continue
        if histo.GetEntries() < 100:
            continue
        if histo.GetRMS() <= 0.:
            continue


        fitFunc = ROOT.TF1('fitFunc','gaus',-10000, 10000)
        fitFunc.SetRange(histo.GetMean()-2.0*histo.GetRMS(),histo.GetMean()+2.0*histo.GetRMS())
        histo.Fit('fitFunc','QNRS+')
        fitFunc.SetRange(fitFunc.GetParameter(1)-3.0*fitFunc.GetParameter(2),fitFunc.GetParameter(1)+3.0*fitFunc.GetParameter(2))
        histo.Fit('fitFunc','QNS+','', fitFunc.GetParameter(1)-3.0*fitFunc.GetParameter(2),fitFunc.GetParameter(1)+3.0*fitFunc.GetParameter(2))
        res = [fitFunc.GetParameter(2), fitFunc.GetParError(2)]
        
        print run, 'bestTh = ', thresh, '   time res = ',  res[0]/math.sqrt(2), '  SR =', SR
        g_tRes_vs_SR[Vov].SetPoint(g_tRes_vs_SR[Vov].GetN(), SR, res[0]/math.sqrt(2))
        g_tRes_vs_SR[Vov].SetPointError(g_tRes_vs_SR[Vov].GetN()-1, 0., math.sqrt(pow(res[1]/math.sqrt(2),2) + 2.*2.))
        
        g_tRes_vs_SR_all.SetPoint(g_tRes_vs_SR_all.GetN(), SR, res[0]/math.sqrt(2))
        g_tRes_vs_SR_all.SetPointError(g_tRes_vs_SR_all.GetN()-1, 0., math.sqrt(pow(res[1]/math.sqrt(2),2) + 2.*2.))
        print("   ")
        


c = ROOT.TCanvas('c_tRes_vs_SR','c_tRes_vs_SR',1200,700)
hPad = ROOT.gPad.DrawFrame(0.,0.,300.,150.)
hPad.SetTitle(";slew rate [#muA/ns];#sigma_{t}^{single} [ps]");
hPad.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();


legend = ROOT.TLegend(0.60, 0.70, 0.89, 0.89)
legend.SetBorderSize(0)

fit = {}
it = 1
for Vov in VovList:
    g_tRes_vs_SR[Vov].SetLineColor(51+8*it)
    g_tRes_vs_SR[Vov].SetMarkerColor(51+8*it)
    g_tRes_vs_SR[Vov].SetMarkerSize(1)
    g_tRes_vs_SR[Vov].Draw('psame')
    legend.AddEntry(g_tRes_vs_SR[Vov], 'V_{OV} = %.1f V'%Vov,'PL')
    fit[Vov] = ROOT.TF1('fit%d'%laserTune,'sqrt([0]*[0] + [1]*[1]/x/x )', 0, 300)
    fit[Vov].SetLineColor(51+8*it)
    fit[Vov].SetLineWidth(1)
    fit[Vov].SetParameters(12, 500)
    g_tRes_vs_SR[Vov].Fit(fit[Vov],'QRNS')
    fit[Vov].Draw('same')
    it += 1

legend.Draw('same')


fitAll = ROOT.TF1('fitAll','sqrt([0]*[0] + [1]*[1]/x/x )', 0, 1000)
fitAll.SetLineColor(ROOT.kBlack)
fitAll.SetLineStyle(2)
fitAll.SetParameters(12, 500)
g_tRes_vs_SR_all.Fit(fitAll,'SR')
fitAll.Draw('same')

latex = ROOT.TLatex(0.62,0.60,'#sigma_{t} = %.02f #muA / (dV/dt) #oplus %.01f ps'%(fitAll.GetParameter(1)/1000.,fitAll.GetParameter(0)))
latex.SetNDC()
latex.SetTextFont(42)
latex.SetTextSize(0.04)
latex.SetTextColor(ROOT.kBlack)
latex.Draw('same')

c.Print("c_tRes_vs_SR.png")

