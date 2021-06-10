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
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetTitleOffset(1.05,'Y')
ROOT.gErrorIgnoreLevel = ROOT.kWarning;

mV_to_pe = 0.15

fnames = {}

vovs = []
vovs.append(5.0)
vovs.append(3.5)
vovs.append(3.0)
vovs.append(2.5)
vovs.append(2.0)
vovs.append(1.6)
vovs.append(1.2)

for vov in vovs:
    fnames[vov] = {}

vths = []
vths.append("vth2")
vths.append("vth1_3")

fnames[5.0]["vth2"] = "/data/Lab5015Analysis/moduleCharacterization/moduleCharacterization_step3_run1037.root"
fnames[5.0]["vth1_3"] = "/data/Lab5015Analysis/moduleCharacterization/moduleCharacterization_step3_run1038.root"
fnames[3.5]["vth2"] = "/data/Lab5015Analysis/moduleCharacterization/moduleCharacterization_step3_run1041.root"
fnames[3.5]["vth1_3"] = "/data/Lab5015Analysis/moduleCharacterization/moduleCharacterization_step3_run1042.root"
fnames[3.0]["vth2"] = "/data/Lab5015Analysis/moduleCharacterization/moduleCharacterization_step3_run1045.root"
fnames[3.0]["vth1_3"] = "/data/Lab5015Analysis/moduleCharacterization/moduleCharacterization_step3_run1046.root"
fnames[2.5]["vth2"] = "/data/Lab5015Analysis/moduleCharacterization/moduleCharacterization_step3_run1049.root"
fnames[2.5]["vth1_3"] = "/data/Lab5015Analysis/moduleCharacterization/moduleCharacterization_step3_run1050.root"
fnames[2.0]["vth2"] = "/data/Lab5015Analysis/moduleCharacterization/moduleCharacterization_step3_run1053.root"
fnames[2.0]["vth1_3"] = "/data/Lab5015Analysis/moduleCharacterization/moduleCharacterization_step3_run1054.root"
fnames[1.6]["vth2"] = "/data/Lab5015Analysis/moduleCharacterization/moduleCharacterization_step3_run1057.root"
fnames[1.6]["vth1_3"] = "/data/Lab5015Analysis/moduleCharacterization/moduleCharacterization_step3_run1058.root"
fnames[1.2]["vth2"] = "/data/Lab5015Analysis/moduleCharacterization/moduleCharacterization_step3_run1061.root"
fnames[1.2]["vth1_3"] = "/data/Lab5015Analysis/moduleCharacterization/moduleCharacterization_step3_run1062.root"

#gtypes = ['timeRes', 'timeResToT']
gtypes = ['timeResToT']

colors = {}
colors["vth2"] = ROOT.kRed-1
colors["vth1_3"] = ROOT.kOrange-3
colors["vth1_1"] = ROOT.kSpring+5
colors["vth1_0"] = ROOT.kAzure-3

vth_to_mV = {}
vth_to_mV["vth2"]   = 8.
vth_to_mV["vth1_4"] = 4.
vth_to_mV["vth1_3"] = 2.
vth_to_mV["vth1_1"] = 1.
vth_to_mV["vth1_0"] = 0.5

plotDir = '/var/www/html/TOFHIR2A/ModuleCharacterization/summary_runs_1035-1062_preAmpGain0/'



graphs_final = {}
graphs_final_pe = {}

for gtype in gtypes:
    
    for vov in vovs:
        c = ROOT.TCanvas('c_Vov%.1f'%vov,'c_Vov%.1f'%vov,1200,600)
        c.Divide(2,1)
        c.cd(1)
        hPad1 = ROOT.gPad.DrawFrame(0.,0.,400.,250.)
        hPad1.SetTitle(";threshold [mV];#sigma_{t_{diff}} / 2 [ps]");
        hPad1.Draw();
        ROOT.gPad.SetGridx();
        ROOT.gPad.SetGridy();
        c.cd(2)
        hPad2 = ROOT.gPad.DrawFrame(0.,20.,60.,120.)
        hPad2.SetTitle(";threshold [mV];#sigma_{t_{diff}} / 2 [ps]");
        hPad2.Draw();
        ROOT.gPad.SetGridx();
        ROOT.gPad.SetGridy();

        c2 = ROOT.TCanvas('c2_Vov%.1f'%vov,'c2_Vov%.1f'%vov,1200,600)
        c2.Divide(2,1)
        c2.cd(1)
        hPad11 = ROOT.gPad.DrawFrame(0.,0.,400./mV_to_pe,250.)
        hPad11.SetTitle(";threshold [p.e.];#sigma_{t_{diff}} / 2 [ps]");
        hPad11.Draw();
        ROOT.gPad.SetGridx();
        ROOT.gPad.SetGridy();
        c2.cd(2)
        hPad22 = ROOT.gPad.DrawFrame(0.,20.,60./mV_to_pe,120.)
        hPad22.SetTitle(";threshold [p.e.];#sigma_{t_{diff}} / 2 [ps]");
        hPad22.Draw();
        ROOT.gPad.SetGridx();
        ROOT.gPad.SetGridy();
        
        legend = ROOT.TLegend(0.50,0.90-0.04*len(fnames[vov]),0.80,0.90)
        legend.SetTextFont(82) 
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        
        
        for label in fnames[vov]:
            inFile = ROOT.TFile.Open(fnames[vov][label])
            graph = inFile.Get('g_%s_vs_th_Vov%.1f_bar00_enBin1'%(gtype,vov))
            
            graphs_final['Vov%.1f'%vov+label+gtype] = ROOT.TGraphErrors()
            graph_final = graphs_final['Vov%.1f'%vov+label+gtype]

            graphs_final_pe['Vov%.1f'%vov+label+gtype] = ROOT.TGraphErrors()
            graph_final_pe = graphs_final_pe['Vov%.1f'%vov+label+gtype]
            
            for point in range(graph.GetN()):
                #if label =='vth2' and point == (graph.GetN() -1):
                #continue
                graph_final.SetPoint(graph_final.GetN(),graph.GetPointX(point)*vth_to_mV[label],graph.GetPointY(point))
                graph_final.SetPointError(graph_final.GetN()-1,0.,graph.GetErrorY(point))
                graph_final_pe.SetPoint(graph_final_pe.GetN(),graph.GetPointX(point)*vth_to_mV[label]/mV_to_pe,graph.GetPointY(point))
                graph_final_pe.SetPointError(graph_final_pe.GetN()-1,0.,graph.GetErrorY(point))
            
            graph_final.SetMarkerColor(colors[label])
            graph_final.SetLineColor(colors[label])
            graph_final_pe.SetMarkerColor(colors[label])
            graph_final_pe.SetLineColor(colors[label])
            
            c.cd(1)
            graph_final.Draw("PL,same")
            c.cd(2)
            graph_final.Draw("PL,same")

            c2.cd(1)
            graph_final_pe.Draw("PL,same")
            c2.cd(2)
            graph_final_pe.Draw("PL,same")
            
            legend.AddEntry(graph_final,"%s"%label,'PL')
          
        c.cd(1)
        legend.Draw("same")
        c.cd(2)
        legend.Draw("same")
        c.Print("%stRes_vs_vth_Vov%.1f.png"%(plotDir,vov))
        
        c2.cd(1)
        legend.Draw("same")
        c2.cd(2)
        legend.Draw("same")
        c2.Print("%stRes_vs_vth_pe_Vov%.1f.png"%(plotDir,vov))




for gtype in gtypes:
    
    for vth in vths:
        
        c = ROOT.TCanvas('c_%s'%vth,'c_%s'%vth,1200,600)
        c.Divide(2,1)
        c.cd(1)
        hPad1 = ROOT.gPad.DrawFrame(0.,0.,400.,250.)
        hPad1.SetTitle(";threshold [mV];#sigma_{t_{diff}} / 2 [ps]");
        hPad1.Draw();
        ROOT.gPad.SetGridx();
        ROOT.gPad.SetGridy();
        c.cd(2)
        hPad2 = ROOT.gPad.DrawFrame(0.,20.,60.,120.)
        hPad2.SetTitle(";threshold [mV];#sigma_{t_{diff}} / 2 [ps]");
        hPad2.Draw();
        ROOT.gPad.SetGridx();
        ROOT.gPad.SetGridy();
        
        c2 = ROOT.TCanvas('c2_%s'%vth,'c2_%s'%vth,1200,600)
        c2.Divide(2,1)
        c2.cd(1)
        hPad11 = ROOT.gPad.DrawFrame(0.,0.,400./mV_to_pe,250.)
        hPad11.SetTitle(";threshold [p.e.];#sigma_{t_{diff}} / 2 [ps]");
        hPad11.Draw();
        ROOT.gPad.SetGridx();
        ROOT.gPad.SetGridy();
        c2.cd(2)
        hPad22 = ROOT.gPad.DrawFrame(0.,20.,100./mV_to_pe,120.)
        hPad22.SetTitle(";threshold [p.e.];#sigma_{t_{diff}} / 2 [ps]");
        hPad22.Draw();
        ROOT.gPad.SetGridx();
        ROOT.gPad.SetGridy();
        
        legend = ROOT.TLegend(0.50,0.90-0.04*len(vovs),0.80,0.90)
        legend.SetTextFont(82) 
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        
        it = 0
        for vov in vovs:
            graph_final = graphs_final['Vov%.1f'%vov+vth+gtype]
            graph_final_pe = graphs_final_pe['Vov%.1f'%vov+vth+gtype]
            
            graph_final.SetMarkerColor(51+it*int(50/len(vovs)))
            graph_final.SetLineColor(51+it*int(50/len(vovs)))
            graph_final_pe.SetMarkerColor(51+it*int(50/len(vovs)))
            graph_final_pe.SetLineColor(51+it*int(50/len(vovs)))
            
            c.cd(1)
            graph_final.Draw("PL,same")
            c.cd(2)
            graph_final.Draw("PL,same")

            c2.cd(1)
            graph_final_pe.Draw("PL,same")
            c2.cd(2)
            graph_final_pe.Draw("PL,same")
            
            legend.AddEntry(graph_final,"V_{OV} = %.1f V"%vov,'PL')
            
            it += 1
        
        c.cd(1)
        legend.Draw("same")
        c.cd(2)
        legend.Draw("same")
        c.Print("%stRes_vs_vth_%s.png"%(plotDir,vth))
        
        c2.cd(1)
        legend.Draw("same")
        c2.cd(2)
        legend.Draw("same")
        c2.Print("%stRes_vs_vth_pe_%s.png"%(plotDir,vth))
