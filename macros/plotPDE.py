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
#ROOT.gROOT.SetBatch(True)
ROOT.gROOT.SetBatch(False)
ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gStyle.SetOptStat(0)

from SiPM import *

#def PDE(ov, sipm):
#    if (sipm=='HPK'):
#        return 1.0228*0.384 * ( 1. - math.exp(-1.*0.583*ov) ) # 1.0228 factor to account for LYSO emission spectrum 
#    if (sipm=='FBK'):
#        return 0.8847*0.466 * ( 1. - math.exp(-1.*0.314*ov) ) # 0.8847 factor to account for LYSO emission spectrum   


g = {}
g['HPK'] = ROOT.TGraph()
g['FBK'] = ROOT.TGraph()
gRatio = ROOT.TGraph()

for i in range(10,800):
    ov = i*0.01
    for sipm in ['HPK','FBK']:
        g[sipm].SetPoint(g[sipm].GetN(), ov, PDE(ov, sipm))

    gRatio.SetPoint(gRatio.GetN(), ov, PDE(ov,'FBK')/PDE(ov,'HPK') )


legend = ROOT.TLegend(0.2,0.7,0.5,0.89)
legend.SetBorderSize(0)

c = ROOT.TCanvas('c_PDE','c_PDE')
c.SetGridx()
c.SetGridy()
c.SetTickx()
c.SetTicky()
hdummy = ROOT.TH2F('hdummy','',100,0,6,40,0,0.45)
hdummy.GetXaxis().SetTitle('V_{OV} [V]')
hdummy.GetYaxis().SetTitle('PDE')
hdummy.Draw()
for i,sipm in enumerate(['HPK','FBK']):
    g[sipm].SetMarkerStyle(20)
    g[sipm].SetMarkerSize(0.2)
    g[sipm].SetMarkerColor(i+1)
    g[sipm].SetLineColor(i+1)
    g[sipm].Draw('plsame')
    legend.AddEntry(g[sipm],'%s'%sipm,'PL')
legend.Draw('same')

c2 = ROOT.TCanvas('c_ratioPDE','c_ratioPDE')
c2.SetGridx()
c2.SetGridy()
c2.SetTickx()
c2.SetTicky()
gRatio.GetXaxis().SetRangeUser(0,6)
gRatio.GetXaxis().SetTitle('V_{OV} [V]')
gRatio.GetYaxis().SetTitle('PDE(FBK)/PDE(HPK)')
gRatio.SetMarkerStyle(20)
gRatio.SetMarkerSize(0.2)
gRatio.Draw('ap')

c.SaveAs(c.GetName()+'.pdf')
c.SaveAs(c.GetName()+'.png')
c2.SaveAs(c2.GetName()+'.pdf')

raw_input('ciao')
