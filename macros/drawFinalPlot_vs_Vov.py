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
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetTitleOffset(1.05,'Y')
ROOT.gErrorIgnoreLevel = ROOT.kWarning;


fnames = {}
#fnames[" 9500 pe (1^{st})"] = "/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run772-799_vth2_step3.root"

# files for preAmpGain1
fnames[" 9500 pe"] = "/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run847-874_vth2_step3.root"
fnames["20000 pe"] = "/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run815-842_vth2_step3.root"

#gtypes = ['timeRes', 'timeResToT']
gtypes = ['timeRes']

colors = {}
#colors[" 9500 pe (1^{st})"] = ROOT.kOrange
colors[" 9500 pe"] = ROOT.kRed
colors["20000 pe"] = ROOT.kBlue

markerStyles = {}
markerStyles["timeRes"] = 20
markerStyles["timeResToT"] = 22

graphs_final = {}


c = ROOT.TCanvas('c','c',600,600)
hPad = ROOT.gPad.DrawFrame(0.,0.,6.,160.)
hPad.SetTitle(";V_{OV} [V];#sigma_{t_{diff}} / 2 [ps]");
hPad.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();

legend = ROOT.TLegend(0.50,0.90-0.04*len(fnames),0.80,0.90)
legend.SetTextFont(82)
legend.SetFillStyle(0)
legend.SetBorderSize(0)

for label in fnames:
    for gtype in gtypes:
        vovs = []
        inFile = ROOT.TFile.Open(fnames[label])
        objs = inFile.GetListOfKeys()
        for obj in objs:
            name = obj.GetName()
            if 'g_%s_vs_th'%gtype in name:
                vovs.append( float(name[12+int(len(gtype)):12+int(len(gtype))+3]) )
            vovs = list(set(vovs))
        vovs.sort()
        
        graphs_final[label+gtype] = ROOT.TGraphErrors()
        graph_final = graphs_final[label+gtype]
        for vov in vovs:
            graph = inFile.Get('g_%s_vs_th_Vov%.1f_bar00_enBin1'%(gtype,vov))
            tResBest = 999999.
            tResBestErr = 0.
            for point in range(graph.GetN()):
                tRes = graph.GetPointY(point)
                if tRes < tResBest and tRes > 10.:
                    tResBest = tRes
                    tResBestErr = graph.GetErrorY(point)
            graph_final.SetPoint(graph_final.GetN(),vov,tResBest)
            graph_final.SetPointError(graph_final.GetN()-1,0.,tResBestErr)
        
        c.cd()
        graph_final.SetMarkerColor(colors[label])
        graph_final.SetLineColor(colors[label])
        graph_final.SetMarkerStyle(markerStyles[gtype])
        graph_final.Draw("PL,same")
        
        legend.AddEntry(graph_final,"%s"%label,'PL')
        
legend.Draw('same')
c.Print("plots/tRes_vs_Vov.png")
