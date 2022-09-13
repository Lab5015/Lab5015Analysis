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
import json
import CMS_lumi, tdrstyle
from collections import OrderedDict

#parser = argparse.ArgumentParser(description='Module characterization summary plots')
##parser.add_argument("-r",  "--runs",          required=True, type=str, help="comma-separated list of runs to be processed")
#parser.add_argument("-i",  "--inputLabels",   required=True, type=str, help="comma-separated list of input labels")
#parser.add_argument("-m",  "--resMode",       required=True, type=int, help="resolution mode: 2 - tDiff, 1 - tAve")
#parser.add_argument("-o",  "--outFolder",     required=True, type=str, help="out folder")
#args = parser.parse_args()

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.04,'X')
ROOT.gStyle.SetLabelSize(0.04,'Y')
ROOT.gStyle.SetTitleSize(0.04,'X')
ROOT.gStyle.SetTitleSize(0.04,'Y')
ROOT.gStyle.SetTitleOffset(1.1,'X')
ROOT.gStyle.SetTitleOffset(1.2,'Y')
ROOT.gROOT.SetBatch(True)



thRef = 7
list_bars = [it for it in range(0, 16)]
list_Vovs = [1.5, 3.5]
list2_Vovs = [1.5, 2.5, 3.5, 5.0]
list_sides = ['L', 'R',]

goodBars = ['00L', '03R', '07L', '10R', '11R', '13L', '14R']

colors = [ ROOT.kRed, ROOT.kOrange, ROOT.kBlue, ROOT.kTeal ] 

infile_linearization = ROOT.TFile('linearizeTOFHIR2X_Jun22.root','READ')

infiles = OrderedDict()
infiles['HPK_nonIrr_T10C'] = [
    ( ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/simona_TB_CERN/Lab5015Analysis/plots/moduleCharacterization_step2_HPK_nonIrr_LYSO528_T10C_Vov1.50_ith2_8.root'), 1.50 )
]
infiles['HPK_2E14_T-40C']  = [
    ( ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/martina_TB_CERN_June22/Lab5015Analysis/plots/moduleCharacterization_step2_HPK_2E14_LYSO796_T-40C_Vov1.70.root'), 1.70 ),
    ( ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/martina_TB_CERN_June22/Lab5015Analysis/plots/moduleCharacterization_step2_HPK_2E14_LYSO796_T-40C_Vov1.90.root'), 1.90 )
]
infiles['HPK_1E14_T-40C']  = [
    ( ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/martina_TB_CERN_June22/Lab5015Analysis/plots/moduleCharacterization_step2_HPK_1E14_LYSO802_T-40C_Vov1.60.root'), 1.60 ),
    ( ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/martina_TB_CERN_June22/Lab5015Analysis/plots/moduleCharacterization_step2_HPK_1E14_LYSO802_T-40C_Vov1.80.root'), 1.80 )
]
infiles['FBK_2E14_T-40C']  = [
    ( ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/martina_TB_CERN_June22/Lab5015Analysis/plots/moduleCharacterization_step2_FBK_2E14_LYSO797_T-40C_Vov1.60.root'), 1.60 ),
    ( ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/martina_TB_CERN_June22/Lab5015Analysis/plots/moduleCharacterization_step2_FBK_2E14_LYSO797_T-40C_Vov1.80.root'), 1.80 )
]
infiles['FBK_1E14_T-40C']  = [
    ( ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/martina_TB_CERN_June22/Lab5015Analysis/plots/moduleCharacterization_step2_FBK_1E14_LYSO803_T-40C_Vov1.40.root'), 1.40 ),
    ( ROOT.TFile('/data1/cmsdaq/Lab5015Analysis_new/martina_TB_CERN_June22/Lab5015Analysis/plots/moduleCharacterization_step2_FBK_1E14_LYSO803_T-40C_Vov1.60.root'), 1.60 )
]
modules = infiles.keys()
module_ref = 'HPK_nonIrr_T10C'

outfile = ROOT.TFile('studyModuleLO_Jun22.root','RECREATE')

with open('/data1/html/TOFHIR2X/MTDTB_CERN_June22/Currents/VovsEff.json') as json_file:
    VovsEff_dict = json.load(json_file)

c1 = ROOT.TCanvas('c1','',700,600)
hPad = ROOT.gPad.DrawFrame(-1.,0.,16.,1.)
hPad.GetXaxis().SetTitle('bar ID')
hPad.GetYaxis().SetTitle('LO / LO_{ref}')
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
hPad.Draw()

leg = ROOT.TLegend(0.50, 0.70, 0.89, 0.90)
leg.SetBorderSize(0)
leg.SetFillStyle(0)

graphs = {}
fitFuncs = {}
for module in modules:
    graphs[module] = ROOT.TGraphErrors()
    fitFuncs[module] = ROOT.TF1('fitFunc_%s'%module,'pol0',-1.,16.)

for module in modules:
    for bar in list_bars:
        for side in list_sides:
            
            aldo = 'A'
            if side == 'R':
                aldo = 'B'
            
            barLabel = '%02d%s'%(bar,side)
            if barLabel not in goodBars:
                continue
            
            graph_temp = ROOT.TGraph()
            
            for infile in infiles[module]:
                
                Vov = infile[1]
                VovEff = Vov
                label = module+'_'+aldo
                if label in VovsEff_dict.keys():
                    print label,'%.2f'%Vov
                    VovEff = float(VovsEff_dict[label]['%.2f'%Vov][0])
                print VovEff
                angle = 52.
                
                histo = infile[0].Get('h1_energy_bar%02d%s_Vov%.2f_th%02d'%(bar,side,Vov,thRef))
                if not histo:
                    continue
                
                fitFunc = histo.GetFunction('f_landau_bar%02d%s_Vov%.2f_vth1_%02.0f'%(bar,side,Vov,thRef))
                if not fitFunc:
                    continue
                
                f_linearization = infile_linearization.Get('f_linearization_angle_bar%02d%s_Vov1.50'%(bar,side))
                energy_linearized = f_linearization.Eval(fitFunc.GetParameter(1))
                print module,bar,side,fitFunc.GetParameter(1),f_linearization.Eval(fitFunc.GetParameter(1))
            
                graph_temp.SetPoint(graph_temp.GetN(),VovEff,energy_linearized)
            
            graphs[module].SetPoint(graphs[module].GetN(),int(bar),graph_temp.Eval(1.5))

for module in modules:
    if module == module_ref:
        continue
    for point in range(graphs[module].GetN()):
        x = graphs[module].GetPointX(point)
        y = graphs[module].GetPointY(point)
        y_ref = graphs[module_ref].GetPointY(point)
        graphs[module].SetPoint(point,x,y/y_ref)
        print module,side,x,y/y_ref
#leg2 = ROOT.TLegend(0.25, 0.75, 0.50, 0.85)
#leg2.SetBorderSize(0)
#leg2.SetFillStyle(0)
#leg2.AddEntry(g_linearization[bar],'MIP Landau MPV','PE')
#leg2.Draw()
it = 0
latexs = {}
for module in modules:
    if module == module_ref:
        continue
    graphs[module].SetMarkerColor(colors[it])
    graphs[module].SetMarkerStyle(20)
    graphs[module].Draw('P,same')
    graphs[module].Fit(fitFuncs[module],'QNRS')
    fitFuncs[module].SetLineColor(colors[it])
    fitFuncs[module].SetLineStyle(2)
    fitFuncs[module].SetLineWidth(2)
    fitFuncs[module].Draw('same')
    latexs[module] = ROOT.TLatex( 0.20, 0.50-0.05*it, '%s:   #LT LO #GT / #LT LO #GT_{ref} = %.2f'%(module,fitFuncs[module].GetParameter(0)))
    latexs[module].SetNDC()
    latexs[module].SetTextSize(0.040)
    latexs[module].SetTextColor(colors[it])
    latexs[module].Draw('same')
    it += 1

c1.Print('c1_moduleLO_Jun22.png')
