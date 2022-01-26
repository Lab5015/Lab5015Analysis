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


plotDir = '/var/www/html/MTDST_CERN_Oct21/CCv2//'
outFile = ROOT.TFile("/home/data/mtd/RUptoro2/Lab5015Analysis_MTDST_CERN_Oct21/plots_fede/drawLaserCTR_vs_SR.root","RECREATE");
outDir = "/var/www/html/MTDST_CERN_Oct21/CCv2/slewRate/"


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

VovList = [1.0, 1.5, 2.0, 2.5, 3.0]

runs_dict = { 

    "2318-2348" : [ 'ith2_0', 0.50, [1.0]],
    "2349-2379" : [ 'ith2_0', 0.60, [1.0]],
    "2380-2410" : [ 'ith2_0', 0.70, [1.0]],
    "2411-2423" : [ 'ith2_0', 0.80, [1.0]],
    "2426-2456" : [ 'ith2_0', 0.75, [1.0]],
    "2457-2487" : [ 'ith2_0', 0.65, [1.0]],
    "2488-2518" : [ 'ith2_0', 0.40, [1.0]],
    "2519-2549" : [ 'ith2_0', 0.30, [1.0]],
    "2550-2578" : [ 'ith2_0', 0.20, [1.0]],
    "2581-2611" : [ 'ith2_0', 0.05, [1.0]],
    "2612-2642" : [ 'ith2_0', 0.50, [1.5]],
    "2643-2673" : [ 'ith2_0', 0.30, [1.5]],
    #"2674-2704" : [ 'ith2_0', 0.10, [1.5]],
    "2705-2765" : [ 'ith2_0', 0.10, [2.5,3.0]]

}

g_tRes_vs_SR = {}
g_tRes_vs_SR_all = ROOT.TGraphErrors()
for Vov in VovList:
    g_tRes_vs_SR[Vov] = ROOT.TGraphErrors()



def get_list(input_list):
    if input_list is None:
        return None
    vals = []
    comma_list = input_list.split(',')
    for item in comma_list:
        hyphen_list = item.split('-')
        if len(hyphen_list) > 1:
            for i in range(int(hyphen_list[0]), int(hyphen_list[1])+1):
                vals.append(i)
        else:
            vals.append(hyphen_list[0])
    return vals


#---------------------------
# -- compute SR vs threshold
npoints = 16
nintervals = 4
pngLabel = ""

for run in sorted(runs_dict):
        print("===>>> Processing run range: ", run)

        fullPath = outDir+"/"+run
        if (os.path.isdir(fullPath) == False):
            os.system('mkdir -p %s'%fullPath)

    
        pngLabel += str(run)+"-"

        ithMode = runs_dict[run][0]
        laserTune = runs_dict[run][1]
        Vovs  = runs_dict[run][2]

        for Vov in Vovs:
            inFile = ROOT.TFile.Open('/home/data/mtd/RUptoro2/Lab5015Analysis_MTDST_CERN_Oct21/plots_fede/pulseShape_run%s.root'%run)
    
            SR = 0.
            thresh = 0.
            channels = ['ch1', 'ch2']

            graph1 = inFile.Get('g_ps_totSel_%s_Vov%.01f'%('ch1',Vov))
            graph2 = inFile.Get('g_ps_totSel_%s_Vov%.01f'%('ch2',Vov))
            if graph1 == None or graph2 == None:
                continue
        
            index_cen = int( min(graph1.GetN()/nintervals, graph2.GetN()/nintervals)) -1
            index_min = max(1, index_cen - int(npoints/2)) -1
            index_max = min(index_cen + int(npoints/2), int(min(graph1.GetN()/2, graph2.GetN()/2))) -1

            for ch in channels:
                canvas = ROOT.TCanvas()
                graph = inFile.Get('g_ps_totSel_%s_Vov%.01f'%(ch,Vov))
                graph.Draw("AP")
                fitSR = ROOT.TF1('fitSR', 'pol1', -10., 30.)
                fitSR.SetRange( graph.GetPointX(index_min)-0.001, graph.GetPointX(index_max)+0.001 )
                fitSR.SetParameters(0,100)
                graph.Fit(fitSR,'QRS')
                canvas.Print(fullPath+'/ps_totSel_%s_Vov%.01f.png'%(ch,Vov))
                outFile.cd()
                graph.Write('Run%s_g_ps_totSel_%s_Vov%.01f'%(run,ch,Vov))
                
                print(ch, index_cen, index_min, index_max, graph.GetX()[index_cen], graph.GetX()[index_min], graph.GetX()[index_max], fitSR.GetParameter(1))
                SR += fitSR.GetParameter(1)
                thresh += int(round(graph.GetPointY(index_cen)/dac_to_uA[ithMode]))
        
            SR /= 2.
            thresh /= 2.
    
    
            inFile = ROOT.TFile.Open('/home/data/mtd/RUptoro2/Lab5015Analysis_MTDST_CERN_Oct21/plots_fede/moduleCharacterization_step2_run%s.root'%run)
        
            print('NAME: h1_deltaT_totRatioCorr_bar00L-R_Vov%.02f_th%.02d_energyBin01'%(Vov,thresh))
            histo = inFile.Get('h1_deltaT_totRatioCorr_bar00L-R_Vov%.02f_th%.02d_energyBin01'%(Vov,thresh))
            if histo == None:
                print("MISSING HISTO")
                continue
            if histo.GetEntries() < 100:
                continue
            if histo.GetRMS() <= 0.:
                continue

            canvas = ROOT.TCanvas()
            histo.Draw()
            fitFunc = ROOT.TF1('fitFunc','gaus',-10000, 10000)
            fitFunc.SetRange(histo.GetMean()-2.0*histo.GetRMS(),histo.GetMean()+2.0*histo.GetRMS())
            histo.Fit('fitFunc','QNRS+')
            fitFunc.SetRange(fitFunc.GetParameter(1)-3.0*fitFunc.GetParameter(2),fitFunc.GetParameter(1)+3.0*fitFunc.GetParameter(2))
            histo.Fit('fitFunc','QNS+','', fitFunc.GetParameter(1)-3.0*fitFunc.GetParameter(2),fitFunc.GetParameter(1)+3.0*fitFunc.GetParameter(2))
            res = [fitFunc.GetParameter(2), fitFunc.GetParError(2)]
            canvas.Print(fullPath+'/tRes_%s_thr%s_Vov%.01f.png'%(ch,thresh,Vov))
        
            print 'runRange = ', run, 'bestTh = ', thresh, '   time res = ',  res[0]/math.sqrt(2), '  SR =', SR
            g_tRes_vs_SR[Vov].SetPoint(g_tRes_vs_SR[Vov].GetN(), SR, res[0]/math.sqrt(2))
            g_tRes_vs_SR[Vov].SetPointError(g_tRes_vs_SR[Vov].GetN()-1, 0., math.sqrt(pow(res[1]/math.sqrt(2),2) + 2.*2.))
            
            g_tRes_vs_SR_all.SetPoint(g_tRes_vs_SR_all.GetN(), SR, res[0]/math.sqrt(2))
            g_tRes_vs_SR_all.SetPointError(g_tRes_vs_SR_all.GetN()-1, 0., math.sqrt(pow(res[1]/math.sqrt(2),2) + 2.*2.))
            print("   ")
        


c = ROOT.TCanvas('c_tRes_vs_SR','c_tRes_vs_SR',1200,700)
hPad = ROOT.gPad.DrawFrame(0.,0.,600.,100.)
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
    fit[Vov] = ROOT.TF1('fit%d'%laserTune,'sqrt([0]*[0] + [1]*[1]/x/x )', 0, 600)
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


c.Print(outDir+"/c_tRes_vs_SR_"+pngLabel+".png")

