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
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.05,'X')
ROOT.gStyle.SetLabelSize(0.05,'Y')
ROOT.gStyle.SetTitleSize(0.05,'X')
ROOT.gStyle.SetTitleSize(0.05,'Y')
ROOT.gStyle.SetTitleOffset(1.0,'X')
ROOT.gStyle.SetTitleOffset(1.2,'Y')                                                                                                                
ROOT.gROOT.SetBatch(True)


def Gain(vov, irr):
    gain = (50738.5+95149*vov)  
    if irr == '2E14':
        gain = gain/1.3
    return gain

def DCR(I, gain):
    return (I*0.001/16/(1.602E-19*gain)/1000000000)


irr  = '2E13'
temp = '-25'

Vbd = {'ch3': 36.69, 
       'ch4': 36.76}
inputDir = '/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2_pedSubtraction_tb/log_IV/ASIC0_HPK_non-irr__ASIC2_HPK_1E13__Tsipm-25C__52deg/'
if (temp == '-40'):
    Vbd = {'ch3': 36.08, 
           'ch4': 36.12}
    inputDir = '/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2_pedSubtraction_tb/log_IV/ASIC0_HPK_non-irr__ASIC2_HPK_1E13__Tsipm-40C__52deg/'

outDir = '/var/www/html/TOFHIR2A/MTDTB_CERN_Jul21/ModuleCharacterization/HPK_1E13_52deg_T%sC_summaryPlots/DCR/'%(temp)
outFileName = '../plots/DCR_vs_Vov_HPK_1E13_52deg_T%sC.root'%temp
outfile = ROOT.TFile(outFileName, 'RECREATE' )

    

if (os.path.isdir(outDir) == False):
    os.system('mkdir %s'%outDir)  

fileList = {} 
fileList['ch3'] = [f for f in os.listdir(inputDir)  if ('ASIC0' not in f and 'root' in f and 'ch3' in f) ]
fileList['ch4'] = [f for f in os.listdir(inputDir)  if ('ASIC0' not in f and 'root' in f and 'ch4' in f) ]
print fileList

g_I_vs_Vbias = {}
g_I_vs_Vov = {}
g_DCR_vs_Vov = {}

for ch in ['ch3','ch4']:
    g_I_vs_Vov[ch] = ROOT.TGraphErrors()
    g_DCR_vs_Vov[ch] = ROOT.TGraphErrors()
    if (ch == 'ch3'): 
        g_I_vs_Vov[ch].SetMarkerStyle(20)
        g_DCR_vs_Vov[ch].SetMarkerStyle(20)
    if (ch == 'ch4'): 
        g_I_vs_Vov[ch].SetMarkerStyle(24)
        g_DCR_vs_Vov[ch].SetMarkerStyle(24)
        
    f = ROOT.TFile.Open(inputDir+'/'+fileList[ch][-1])
    g_I_vs_Vbias[ch] = f.Get('g_IV')

    ctemp = ROOT.TCanvas('ctemp','ctemp')
    g_I_vs_Vbias[ch].SetMarkerStyle(20)
    g_I_vs_Vbias[ch].GetYaxis().SetRangeUser(0, 2000)
    g_I_vs_Vbias[ch].GetYaxis().SetTitle('I_array [#muA]')
    g_I_vs_Vbias[ch].GetXaxis().SetTitle('V_{bias} [V]')
    g_I_vs_Vbias[ch].Draw('ap')
    ctemp.SaveAs(outDir+'/Iarray_vs_Vbias_%s.png'%(ch))

    for i in range(0,g_I_vs_Vbias[ch].GetN()):
        vbias = g_I_vs_Vbias[ch].GetX()[i]
        vov   = vbias - Vbd[ch] 
        I_array = g_I_vs_Vbias[ch].Eval(vbias)/1000.
        gain = Gain(vov, irr)
        dcr  = DCR(I_array, gain)
        #dcr_err = 0.5*(DCR(I_array+I_array_err, gain)-DCR(I_array-I_array_err, gain))
        dcr_err = 0
        #print vov, dcr
        g_I_vs_Vov[ch].SetPoint(g_I_vs_Vov[ch].GetN(), vov, I_array)
        g_I_vs_Vov[ch].SetPointError(g_I_vs_Vov[ch].GetN()-1, 0, 0)

        g_DCR_vs_Vov[ch].SetPoint(g_DCR_vs_Vov[ch].GetN(), vov, dcr)
        g_DCR_vs_Vov[ch].SetPoint(g_DCR_vs_Vov[ch].GetN(), 0, dcr_err)



c1 = ROOT.TCanvas()
c1.SetGridx()
c1.SetGridy()
max = g_I_vs_Vov['ch3'].GetY()[g_I_vs_Vov['ch3'].GetN()-1]* 1.5
hdummy = ROOT.TH2F('hdummy','',100,0,2.0,100,0,max)
hdummy.GetXaxis().SetTitle('V_{OV} [V]')
hdummy.GetYaxis().SetTitle('I_{array} (mA)')
hdummy.Draw()
g_I_vs_Vov['ch3'].Draw('psame')
g_I_vs_Vov['ch4'].Draw('psame')


c2 = ROOT.TCanvas()
c2.SetGridx()
c2.SetGridy()
max = g_DCR_vs_Vov['ch3'].GetY()[g_DCR_vs_Vov['ch3'].GetN()-2]* 1.5
hdummy2 = ROOT.TH2F('hdummy2','',100,0,2.0,100,0,max)
hdummy2.GetXaxis().SetTitle('V_{OV} [V]')
hdummy2.GetYaxis().SetTitle('DCR [GHz]')
hdummy2.Draw()
g_DCR_vs_Vov['ch3'].Draw('psame')
g_DCR_vs_Vov['ch4'].Draw('psame')

outfile.cd()
g_DCR_vs_Vov['ch3'].Write('g_DCR_vs_Vov_ch3')
g_DCR_vs_Vov['ch4'].Write('g_DCR_vs_Vov_ch4')  
outfile.Close()

c1.SaveAs(outDir+'/Iarray_vs_Vov.png')
c1.SaveAs(outDir+'/Iarray_vs_Vov.pdf')

c2.SaveAs(outDir+'/DCR_vs_Vov.png')
c2.SaveAs(outDir+'/DCR_vs_Vov.pdf')

raw_input('OK?')

