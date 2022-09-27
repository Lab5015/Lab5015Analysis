#! /usr/bin/env python
import os
import shutil
import glob
import math
import array
import sys
import time
import ROOT
import json
from collections import OrderedDict
#set the tdr style
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
ROOT.gROOT.SetBatch(True)
#ROOT.gROOT.SetBatch(False)



elec = 1.602E-19

def TGraphRatio(g1,g2,scale=1.):
    g = ROOT.TGraphErrors()
    xmin = max(g1.GetPointX(0),g2.GetPointX(0))
    xmax = min(g1.GetPointX(g1.GetN()-1),g2.GetPointX(g2.GetN()-1))
    npoints = 100
    for point in range(npoints):
        Vov = xmin + 1.*(xmax-xmin)/npoints*point
        ratio = g2.Eval(Vov)/g1.Eval(Vov)/scale
        g.SetPoint(point,Vov,ratio)
    return g


def myfunc_DCR(x, par):
    xx = x[0]
    x0 = par[4]
    A = par[0]
    C = par[1]
    D = par[2]
    E = par[3]
    B = D*E*math.exp(E*x0) - 2.*C*x0
    F = A + B*x0 + C*x0*x0 - D*math.exp(E*x0)
    if xx < x0 :
        return par[5] * ( A + B*xx + C*xx*xx )
    else :
        return par[5] * ( D*math.exp(E*xx) + F )


def Gain(ov, sipm, irr='0'):
    k = 1.
    if (irr == '2E14' and 'HPK' in sipm): k = 0.92 # gain reduction for HPK 2E14 irradiated SiPMs 
    if (irr == '1E14' and 'HPK' in sipm): k = 0.96 # gain reduction for HPK 1E14 irradiated SiPMs (assume that for 1E14 is half of 2E14)
    
    # HPK - simple linear fit to data points from some slides
    if ('HPK' in sipm):
        return k*(36890. + 97602.*ov)
    
    #if ('HPK' in sipm):
    #    return k*(ov+0.25)*(12.940-0.2*ov+3.2)/1.602/0.0001      # Arjan defininition
    
    # FBK-MS
    #if ('FBK' in sipm):
    #    return k*(50739. + 95149.*ov) # FBK-MS
    # FBK-W4C 
    if ('FBK' in sipm):
        return k*91541.7*(ov+0.408182) # FBK-W4C




##################
### define infiles

outdir = '/var/www/html/TOFHIR2B/MTDTB_CERN_June22/Currents/'


infilenames_IV = OrderedDict()

infilenames_IV[('HPK','2E14',-40,'A')] = '/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_jun22/logs_14.03/logIV_ASIC0_ALDOA_ch5_time_2022-07-16_16:30:01.root'
infilenames_IV[('HPK','2E14',-40,'B')] = '/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_jun22/logs_14.03/logIV_ASIC0_ALDOB_ch6_time_2022-07-16_16:37:07.root'
infilenames_IV[('HPK','2E14',-35,'A')] = '/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_jun22/logs_14.01/logIV_ASIC0_ALDOA_ch5_time_2022-07-16_10:58:57.root'
infilenames_IV[('HPK','2E14',-35,'B')] = '/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_jun22/logs_14.01/logIV_ASIC0_ALDOB_ch6_time_2022-07-16_10:48:54.root'

infilenames_IV[('HPK','1E14',-35,'A')] = '/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_jun22/logs_15.00/logIV_ASIC0_ALDOA_ch5_time_2022-07-17_13:19:01.root'
infilenames_IV[('HPK','1E14',-35,'B')] = '/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_jun22/logs_15.00/logIV_ASIC0_ALDOB_ch6_time_2022-07-17_13:12:16.root'
#infilenames_IV[('HPK','1E14',-40,'A')] = ''
#infilenames_IV[('HPK','1E14',-40,'B')] = ''

infilenames_IV[('FBK','1E14',-35,'A')] = '/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_jun22/logs_16.00/logIV_ASIC0_ALDOA_ch5_time_2022-07-17_18:50:08.root'
infilenames_IV[('FBK','1E14',-35,'B')] = '/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_jun22/logs_16.00/logIV_ASIC0_ALDOB_ch6_time_2022-07-17_18:44:57.root'
#infilenames_IV[('FBK','1E14',-40,'A')] = ''
#infilenames_IV[('FBK','1E14',-40,'B')] = ''

infilenames_IV[('FBK','2E14',-35,'A')] = '/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_jun22/logs_17.00/logIV_ASIC0_ALDOA_ch5_time_2022-07-18_13:24:57.root'
infilenames_IV[('FBK','2E14',-35,'B')] = '/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_jun22/logs_17.00/logIV_ASIC0_ALDOB_ch6_time_2022-07-18_13:29:06.root'
#infilenames_IV[('FBK','2E14',-40,'A')] = ''
#infilenames_IV[('FBK','2E14',-40,'B')] = ''

VovsSet = OrderedDict()
VovsSet[('HPK','2E14',-40)] = [ 1.40, 1.60, 2.00 ]
VovsSet[('HPK','2E14',-35)] = [ 1.40, 1.60, 2.00 , 2.40, 2.70 ]
VovsSet[('HPK','1E14',-40)] = [ ]
VovsSet[('HPK','1E14',-35)] = [ 1.40, 1.60, 1.80, 2.00 ]

VovsSet[('FBK','2E14',-40)] = [ ]
VovsSet[('FBK','2E14',-35)] = [ 1.40, 1.60, 1.70, 1.80 ]
VovsSet[('FBK','1E14',-40)] = [ ]
VovsSet[('FBK','1E14',-35)] = [ 1.60, 1.70, 1.80, 2.00 ]
                              
Vbrs = OrderedDict()

Vbrs[('HPK','2E14',-40,'A')] = 38.39
Vbrs[('HPK','2E14',-40,'B')] = 38.52
Vbrs[('HPK','2E14',-35,'A')] = 38.59
Vbrs[('HPK','2E14',-35,'B')] = 38.69

Vbrs[('HPK','1E14',-35,'A')] = 37.59
Vbrs[('HPK','1E14',-35,'B')] = 37.57
Vbrs[('HPK','1E14',-40,'A')] = 37.35
Vbrs[('HPK','1E14',-40,'B')] = 37.32

Vbrs[('FBK','1E14',-35,'A')] = 32.64
Vbrs[('FBK','1E14',-35,'B')] = 33.00
Vbrs[('FBK','1E14',-40,'A')] = 32.54
Vbrs[('FBK','1E14',-40,'B')] = 32.79

Vbrs[('FBK','2E14',-35,'A')] = 33.83
Vbrs[('FBK','2E14',-35,'B')] = 33.06
Vbrs[('FBK','2E14',-40,'A')] = 33.53
Vbrs[('FBK','2E14',-40,'B')] = 32.77

VbrsUsed = OrderedDict()

VbrsUsed[('HPK','2E14',-40,'A')] = 38.39
VbrsUsed[('HPK','2E14',-40,'B')] = 38.52
VbrsUsed[('HPK','2E14',-35,'A')] = 38.59
VbrsUsed[('HPK','2E14',-35,'B')] = 38.69

VbrsUsed[('HPK','1E14',-35,'A')] = 37.59
VbrsUsed[('HPK','1E14',-35,'B')] = 37.57
VbrsUsed[('HPK','1E14',-40,'A')] = 37.35
VbrsUsed[('HPK','1E14',-40,'B')] = 37.30

VbrsUsed[('FBK','1E14',-35,'A')] = 32.64
VbrsUsed[('FBK','1E14',-35,'B')] = 33.00
VbrsUsed[('FBK','1E14',-40,'A')] = 32.54
VbrsUsed[('FBK','1E14',-40,'B')] = 32.79

VbrsUsed[('FBK','2E14',-35,'A')] = 33.83
VbrsUsed[('FBK','2E14',-35,'B')] = 33.06
VbrsUsed[('FBK','2E14',-40,'A')] = 33.53
VbrsUsed[('FBK','2E14',-40,'B')] = 32.77

vendors = []
fluences = []
temps = []
aldos = []
keys = infilenames_IV.keys()
for key in keys:
    vendors.append(key[0])
    fluences.append(key[1])
    temps.append(key[2])
    aldos.append(key[3])
vendors = list(dict.fromkeys(vendors))
fluences = list(dict.fromkeys(fluences))
temps = list(dict.fromkeys(temps))
aldos = list(dict.fromkeys(aldos))
aldos2 = aldos + ['AB']



############################
### define output dictionary

output_dict = OrderedDict()




##########################
### define and fill graphs

g_I_Vset = {}
g_I_Vset_avgLast = {}
g_I_OVeff_avgAll = {}
g_I_OVeff_avgLast = {}
g_DCR = {}

for key in keys:
    vendor = key[0]
    fluence = key[1]
    temp = key[2]
    aldo = key[3]
    Vbr = Vbrs[key]
    #print(vendor,fluence,temp,aldo,Vbr)
    
    infile = ROOT.TFile(infilenames_IV[key],'READ')
    g_I_Vset[key] = infile.Get("g_IV")
    g_I_Vset_avgLast[key] = ROOT.TGraphErrors()
    g_I_OVeff_avgAll[key] = ROOT.TGraphErrors()
    g_I_OVeff_avgLast[key] = ROOT.TGraphErrors()
    g_DCR[key] = ROOT.TGraphErrors()

    output_dict['%s_%s_T%dC_%s'%(vendor,fluence,temp,aldo)] = {}
    for point in range(g_I_Vset[key].GetN()):
        Vset = g_I_Vset[key].GetPointX(point)
        I = g_I_Vset[key].GetPointY(point)
        Veff = Vset-Vbr-15*I*1E-06 # JulyTB only 240 ohm / 16
        g_I_OVeff_avgAll[key].SetPoint(point,Veff,I*1E-03/16.)
        #print(point,Veff,I*1E-03)
        
        g_I = infile.Get('g_I_bv%.2f'%Vset)
        nmeas = g_I.GetN()
        xmin = g_I.GetPointX(int(round(2./3.*nmeas))-1)-0.5
        xmax = g_I.GetPointX(nmeas-1)+0.5
        fitFunc = ROOT.TF1('fitFunc','pol0',xmin,xmax)
        g_I.Fit(fitFunc,'QNRS')
        current = fitFunc.GetParameter(0)*1E-03/16.
        dcr = current*1E-03/elec/Gain(Veff,'HPK',fluence)/1E09
        g_I_Vset_avgLast[key].SetPoint(point,Vset,current)
        g_I_OVeff_avgLast[key].SetPoint(point,Veff,current)
        g_I_OVeff_avgLast[key].SetPointError(point,0.,fitFunc.GetParError(0)*1E-03/16.)
        g_DCR[key].SetPoint(point,Veff,dcr)


for key in keys:
    vendor = key[0]
    fluence = key[1]
    temp = key[2]
    aldo = key[3]
    Vbr = Vbrs[key]
    VbrUsed = VbrsUsed[key]

    output_dict['%s_%s_T%dC_%s'%(vendor,fluence,temp,aldo)] = OrderedDict()
    
    for VovSet in VovsSet[(vendor,fluence,temp)]:

        current = g_I_Vset_avgLast[key].Eval(VovSet+VbrUsed)
        Veff = VovSet+VbrUsed-Vbr-15*current*16/1E03 # Jul22 TB --> Vdrop only on 240 ohm/16 resistor
        dcr = g_DCR[key].Eval(Veff)
        output_dict['%s_%s_T%dC_%s'%(vendor,fluence,temp,aldo)]['%.2f'%(VovSet)] = ['%.2f'%round(Veff,2),'%.1f'%round(dcr,1)]
        print vendor,fluence,temp,VovSet,VbrUsed,current,Veff
        
with open(outdir+'VovsEff.json', 'w') as fp:
    json.dump(output_dict, fp, indent=2)




### ALDOA / ALDOB average
for vendor in vendors:
    for fluence in fluences:
        for temp in temps:
            
            keyA = (vendor,fluence,temp,'A')
            keyB = (vendor,fluence,temp,'B')
            keyAB = (vendor,fluence,temp,'AB')
            
            if keyA not in infilenames_IV.keys():
                continue
            if keyB not in infilenames_IV.keys():
                continue
            
            g_I_OVeff_avgLast[keyAB] = ROOT.TGraphErrors()
            g_DCR[keyAB] = ROOT.TGraphErrors()
            xmin = max(g_I_OVeff_avgLast[keyA].GetPointX(0),g_I_OVeff_avgLast[keyB].GetPointX(0))
            xmax = min(g_I_OVeff_avgLast[keyA].GetPointX(g_I_OVeff_avgLast[keyA].GetN()-1),g_I_OVeff_avgLast[keyB].GetPointX(g_I_OVeff_avgLast[keyB].GetN()-1))
            npoints = 100
            for point in range(npoints):
                x = xmin+1.*(xmax-xmin)/npoints*point
                g_I_OVeff_avgLast[keyAB].SetPoint(point,x,0.5*(g_I_OVeff_avgLast[keyA].Eval(x)+g_I_OVeff_avgLast[keyB].Eval(x)))
                g_DCR[keyAB].SetPoint(point,x,0.5*(g_DCR[keyA].Eval(x)+g_DCR[keyB].Eval(x)))




##################
### draw IV curves

canvasXmin = -2.
canvasXmax = 3.5
canvasYmin = 0.
canvasYmax = 3.5
canvasYminlog = 0.00001
canvasYmaxlog = 10.
for vendor in vendors:
    for fluence in fluences:
        for temp in temps:
            
            keyA = (vendor,fluence,temp,'A')
            keyB = (vendor,fluence,temp,'B')
            keyAB = (vendor,fluence,temp,'AB')
            if keyA not in infilenames_IV.keys():
                continue
            if keyB not in infilenames_IV.keys():
                continue
            
            c = ROOT.TCanvas('c_IV_%s_%s_%dC'%(vendor,fluence,temp),'c_IV_%s_%s_%dC'%(vendor,fluence,temp),1600,700)
            c.Divide(2,1)
            c.cd(1)
            hPad1 = ROOT.gPad.DrawFrame(canvasXmin,canvasYmin,canvasXmax,canvasYmax)
            hPad1.SetTitle(";V_{ov}^{eff} [V]; I_{array} / 16 [mA]")
            hPad1.Draw()
            ROOT.gPad.SetGridx()
            ROOT.gPad.SetGridy()
            c.cd(2)
            hPad2 = ROOT.gPad.DrawFrame(canvasXmin,canvasYminlog,canvasXmax,canvasYmaxlog)
            hPad2.SetTitle(";V_{ov}^{eff} [V]; I_{array} / 16 [mA]");
            hPad2.Draw()
            ROOT.gPad.SetLogy()
            ROOT.gPad.SetGridx()
            ROOT.gPad.SetGridy()
            
            label = ROOT.TLatex(0.17,0.94,'%s %s - T = %d#circ C'%(vendor,fluence,temp))
            label.SetNDC()
            label.SetTextSize(0.04)
            label.SetTextFont(42)
            
            latex = {}
            it = 1
            for aldo in aldos2:
                key = (vendor,fluence,temp,aldo)
                print key
                c.cd(1)
                if aldo != 'AB':
                    g_I_OVeff_avgLast[key].SetMarkerColor(1+it)
                    g_I_OVeff_avgLast[key].SetLineColor(1+it)
                    g_I_OVeff_avgLast[key].SetMarkerSize(1.)
                    g_I_OVeff_avgLast[key].SetMarkerStyle(20)
                    g_I_OVeff_avgLast[key].Draw('PEL,same')
                else:
                    g_I_OVeff_avgLast[key].SetLineColor(ROOT.kBlack)
                    g_I_OVeff_avgLast[key].Draw('L,same')
                    linex = ROOT.TLine(canvasXmin,g_I_OVeff_avgLast[key].Eval(1.5),1.5,g_I_OVeff_avgLast[key].Eval(1.5))
                    linex.SetLineStyle(7)
                    linex.Draw('same')
                    liney = ROOT.TLine(1.5,canvasYmin,1.5,g_I_OVeff_avgLast[key].Eval(1.5))
                    liney.SetLineStyle(7)
                    liney.Draw('same')
                c.cd(2)
                if aldo != 'AB':
                    g_I_OVeff_avgLast[key].Draw('PEL,same')
                else:
                    g_I_OVeff_avgLast[key].Draw('L,same')                    
                    linex.Draw('same')
                    lineylog = ROOT.TLine(1.5,canvasYminlog,1.5,g_I_OVeff_avgLast[key].Eval(1.5))
                    lineylog.SetLineStyle(7)
                    lineylog.Draw('same')
                    
                latex[aldo] = ROOT.TLatex(0.20,0.75-0.05*it,'ALDO %s: I_{ch.} = %.2f mA'%(aldo,g_I_OVeff_avgLast[key].Eval(1.5)))
                latex[aldo].SetNDC()
                latex[aldo].SetTextSize(0.04)
                latex[aldo].SetTextFont(42)
                latex[aldo].SetTextColor(1+it)
                if aldo == 'AB':
                    latex[aldo].SetTextColor(ROOT.kBlack)
                c.cd(1)
                latex[aldo].Draw('same')
                c.cd(2)
                latex[aldo].Draw('same')
                
                it += 1
            c.cd(1)
            label.Draw('same')
            c.cd(2)                
            label.Draw('same')
            c.Print(outdir+'%s.png'%c.GetName())
            
            outfile = ROOT.TFile('%s.root'%c.GetName(),'RECREATE')
            g_I_OVeff_avgLast[keyA].Write('g_IV_aldoA')
            g_I_OVeff_avgLast[keyB].Write('g_IV_aldoB')
            g_I_OVeff_avgLast[keyAB].Write('g_IV_aldoAB')




###################
### draw DCR curves

canvasXmin = 0.
canvasXmax = 3.5
canvasYmin = 0.
canvasYmax = 70.
canvasYminlog = 0.1
canvasYmaxlog = 100.
for vendor in vendors:
    for fluence in fluences:
        for temp in temps:
            
            keyA = (vendor,fluence,temp,'A')
            keyB = (vendor,fluence,temp,'B')
            keyAB = (vendor,fluence,temp,'AB')
            if keyA not in infilenames_IV.keys():
                continue
            if keyB not in infilenames_IV.keys():
                continue
            
            c = ROOT.TCanvas('c_DCR_%s_%s_%dC'%(vendor,fluence,temp),'c_DCR_%s_%s_%dC'%(vendor,fluence,temp),1600,700)
            c.Divide(2,1)
            c.cd(1)
            hPad1 = ROOT.gPad.DrawFrame(canvasXmin,canvasYmin,canvasXmax,canvasYmax)
            hPad1.SetTitle(";V_{ov}^{eff} [V]; DCR [GHz]")
            hPad1.Draw()
            ROOT.gPad.SetGridx()
            ROOT.gPad.SetGridy()
            c.cd(2)
            hPad2 = ROOT.gPad.DrawFrame(canvasXmin,canvasYminlog,canvasXmax,canvasYmaxlog)
            hPad2.SetTitle(";V_{ov}^{eff} [V]; DCR [GHz]");
            hPad2.Draw()
            ROOT.gPad.SetLogy()
            ROOT.gPad.SetGridx()
            ROOT.gPad.SetGridy()
            
            label = ROOT.TLatex(0.17,0.94,'%s %s - T = %d#circ C'%(vendor,fluence,temp))
            label.SetNDC()
            label.SetTextSize(0.04)
            label.SetTextFont(42)
            
            latex = {}
            it = 1
            for aldo in aldos2:
                key = (vendor,fluence,temp,aldo)
                print key
                c.cd(1)
                if aldo != 'AB':
                    g_DCR[key].SetMarkerColor(1+it)
                    g_DCR[key].SetLineColor(1+it)
                    g_DCR[key].SetMarkerSize(1.)
                    g_DCR[key].SetMarkerStyle(20)
                    g_DCR[key].Draw('PEL,same')
                else:
                    g_DCR[key].SetLineColor(ROOT.kBlack)
                    g_DCR[key].Draw('L,same')
                    linex = ROOT.TLine(canvasXmin,g_DCR[key].Eval(1.5),1.5,g_DCR[key].Eval(1.5))
                    linex.SetLineStyle(7)
                    linex.Draw('same')
                    liney = ROOT.TLine(1.5,canvasYmin,1.5,g_DCR[key].Eval(1.5))
                    liney.SetLineStyle(7)
                    liney.Draw('same')
                c.cd(2)
                if aldo != 'AB':
                    g_DCR[key].Draw('PEL,same')
                else:
                    g_DCR[key].Draw('L,same')                    
                    linex.Draw('same')
                    lineylog = ROOT.TLine(1.5,canvasYminlog,1.5,g_DCR[key].Eval(1.5))
                    lineylog.SetLineStyle(7)
                    lineylog.Draw('same')
                    
                latex[aldo] = ROOT.TLatex(0.50,0.40-0.05*it,'ALDO %s: DCR = %.1f GHz'%(aldo,g_DCR[key].Eval(1.5)))
                latex[aldo].SetNDC()
                latex[aldo].SetTextSize(0.04)
                latex[aldo].SetTextFont(42)
                latex[aldo].SetTextColor(1+it)
                if aldo == 'AB':
                    latex[aldo].SetTextColor(ROOT.kBlack)
                c.cd(1)
                latex[aldo].Draw('same')
                c.cd(2)
                latex[aldo].Draw('same')
                
                it += 1
            c.cd(1)
            label.Draw('same')
            c.cd(2)                
            label.Draw('same')
            c.Print(outdir+'%s.png'%c.GetName())

            outfile = ROOT.TFile('%s.root'%c.GetName(),'RECREATE')
            g_DCR[keyA].Write('g_DCR_aldoA')
            g_DCR[keyB].Write('g_DCR_aldoB')
            g_DCR[keyAB].Write('g_DCR_aldoAB')




##################
# make comparisons

canvasXmin = 0.
canvasXmax = 3.5
canvasYmin = 0.
canvasYmax = 2.
canvasYmin2 = 0.
canvasYmax2 = 3.
for vendor in vendors:
    for fluence in fluences:
        
        keyA1 = (vendor,fluence,-40,'A')
        keyB1 = (vendor,fluence,-40,'B')
        if keyA1 not in infilenames_IV.keys():
            continue
        if keyB1 not in infilenames_IV.keys():
            continue
        
        keyA2 = (vendor,fluence,-35,'A')
        keyB2 = (vendor,fluence,-35,'B')
        if keyA2 not in infilenames_IV.keys():
            continue
        if keyB2 not in infilenames_IV.keys():
            continue
        
        graphs = {}
        graphs['A'] = TGraphRatio(g_I_OVeff_avgLast[keyA1],g_I_OVeff_avgLast[keyA2])
        graphs['B'] = TGraphRatio(g_I_OVeff_avgLast[keyB1],g_I_OVeff_avgLast[keyB2])
        
        graphs2 = {}
        for aldo in aldos:
            graphs2[aldo] = ROOT.TGraphErrors()
            for point in range(graphs[aldo].GetN()):
                graphs2[aldo].SetPoint(point,graphs[aldo].GetPointX(point),pow(graphs[aldo].GetPointY(point),2))
        
        fits = {}
        for aldo in aldos:
            fits[aldo] = ROOT.TF1('fit_%s'%aldo,'[0]^(0.5)',0.5,10.)
            fits[aldo].SetParameter(0,2.)
        
        c = ROOT.TCanvas('c_IV_%s_%s_tempRatio'%(vendor,fluence),'c_IV_%s_%s_tempRatio'%(vendor,fluence),1600,700)
        c.Divide(2,1)
        c.cd(1)
        hPad1 = ROOT.gPad.DrawFrame(canvasXmin,canvasYmin,canvasXmax,canvasYmax)
        hPad1.SetTitle(";V_{ov}^{eff} [V]; I_{T = -35#circ C} / I_{T = -40#circ C}")
        hPad1.Draw()
        ROOT.gPad.SetGridx()
        ROOT.gPad.SetGridy()
        c.cd(2)
        hPad2 = ROOT.gPad.DrawFrame(canvasXmin,canvasYmin2,canvasXmax,canvasYmax2)
        hPad2.SetTitle(";V_{ov}^{eff} [V]; T_{coeff}")
        hPad2.Draw()
        ROOT.gPad.SetGridx()
        ROOT.gPad.SetGridy()
        
        label = ROOT.TLatex(0.17,0.94,'%s %s'%(vendor,fluence))
        label.SetNDC()
        label.SetTextSize(0.04)
        label.SetTextFont(42)
        
        c.cd(1)
        graphs['A'].SetLineWidth(2)
        graphs['A'].SetLineColor(ROOT.kRed)
        graphs['A'].Draw('L,same')
        graphs['B'].SetLineWidth(2)
        graphs['B'].SetLineColor(ROOT.kGreen)
        graphs['B'].Draw('L,same')
        label.Draw('same')

        c.cd(2)
        graphs2['A'].SetLineWidth(2)
        graphs2['A'].SetLineColor(ROOT.kRed)
        graphs2['A'].Draw('L,same')
        graphs2['B'].SetLineWidth(2)
        graphs2['B'].SetLineColor(ROOT.kGreen)
        graphs2['B'].Draw('L,same')
        label.Draw('same')
        
        c.Print(outdir+'%s.png'%c.GetName())




canvasXmin = 0.
canvasXmax = 3.5
canvasYmin = 0.
canvasYmax = 3.5
canvasYmin2 = 0.
canvasYmax2 = 3.
for vendor in vendors:
    for temp in temps:
        
        keyA1 = (vendor,'1E14',temp,'A')
        keyB1 = (vendor,'1E14',temp,'B')
        if keyA1 not in infilenames_IV.keys():
            continue
        if keyB1 not in infilenames_IV.keys():
            continue
        
        keyA2 = (vendor,'2E14',temp,'A')
        keyB2 = (vendor,'2E14',temp,'B')
        if keyA2 not in infilenames_IV.keys():
            continue
        if keyB2 not in infilenames_IV.keys():
            continue
        
        graphs = {}
        graphs['A'] = TGraphRatio(g_I_OVeff_avgLast[keyA1],g_I_OVeff_avgLast[keyA2],2.)
        graphs['B'] = TGraphRatio(g_I_OVeff_avgLast[keyB1],g_I_OVeff_avgLast[keyB2],2.)
        
        c = ROOT.TCanvas('c_IV_%s_%dC_fluenceRatio'%(vendor,temp),'c_IV_%s_%dC_fluenceRatio'%(vendor,temp),1600,700)
        c.Divide(2,1)
        c.cd(1)
        hPad1 = ROOT.gPad.DrawFrame(canvasXmin,canvasYmin,canvasXmax,canvasYmax)
        hPad1.SetTitle(";V_{ov}^{eff} [V]; I_{array} / 16 [mA]")
        hPad1.Draw()
        ROOT.gPad.SetGridx()
        ROOT.gPad.SetGridy()
        c.cd(2)
        hPad2 = ROOT.gPad.DrawFrame(canvasXmin,canvasYmin2,canvasXmax,canvasYmax2)
        hPad2.SetTitle(";V_{ov}^{eff} [V]; I_{2E14} / ( 2 #times I_{1E14} )")
        hPad2.Draw()
        ROOT.gPad.SetGridx()
        ROOT.gPad.SetGridy()
        
        label = ROOT.TLatex(0.17,0.94,'%s %dC'%(vendor,temp))
        label.SetNDC()
        label.SetTextSize(0.04)
        label.SetTextFont(42)
        
        c.cd(1)
        g_I_OVeff_avgLast[keyA1].SetMarkerColor(ROOT.kRed)
        g_I_OVeff_avgLast[keyA1].SetLineColor(ROOT.kRed)
        g_I_OVeff_avgLast[keyA1].SetMarkerSize(1.)
        g_I_OVeff_avgLast[keyA1].SetMarkerStyle(20)
        g_I_OVeff_avgLast[keyA1].Draw('PEL,same')
        g_I_OVeff_avgLast[keyB1].SetMarkerColor(ROOT.kGreen)
        g_I_OVeff_avgLast[keyB1].SetLineColor(ROOT.kGreen)
        g_I_OVeff_avgLast[keyB1].SetMarkerSize(1.)
        g_I_OVeff_avgLast[keyB1].SetMarkerStyle(20)
        g_I_OVeff_avgLast[keyB1].Draw('PEL,same')
        g_I_OVeff_avgLast[keyA2].SetMarkerColor(ROOT.kRed)
        g_I_OVeff_avgLast[keyA2].SetLineColor(ROOT.kRed)
        g_I_OVeff_avgLast[keyA2].SetMarkerSize(1.)
        g_I_OVeff_avgLast[keyA2].SetMarkerStyle(24)
        g_I_OVeff_avgLast[keyA2].Draw('PEL,same')
        g_I_OVeff_avgLast[keyB2].SetMarkerColor(ROOT.kGreen)
        g_I_OVeff_avgLast[keyB2].SetLineColor(ROOT.kGreen)
        g_I_OVeff_avgLast[keyB2].SetMarkerSize(1.)
        g_I_OVeff_avgLast[keyB2].SetMarkerStyle(24)
        g_I_OVeff_avgLast[keyB2].Draw('PEL,same')
        label.Draw('same')

        c.cd(2)
        graphs['A'].SetLineWidth(4)
        graphs['A'].SetLineColor(ROOT.kRed)
        graphs['A'].Draw('L,same')
        graphs['B'].SetLineWidth(4)
        graphs['B'].SetLineColor(ROOT.kGreen)
        graphs['B'].Draw('L,same')
        label.Draw('same')
        
        c.Print(outdir+'%s.png'%c.GetName())




canvasXmin = 0.
canvasXmax = 3.5
canvasYmin = 0.
canvasYmax = 3.5
canvasYmin2 = 0.
canvasYmax2 = 3.
for fluence in fluences:
    for temp in temps:
        
        keyA1 = ('HPK',fluence,temp,'A')
        keyB1 = ('HPK',fluence,temp,'B')
        if keyA1 not in infilenames_IV.keys():
            continue
        if keyB1 not in infilenames_IV.keys():
            continue
        
        keyA2 = ('FBK',fluence,temp,'A')
        keyB2 = ('FBK',fluence,temp,'B')
        if keyA2 not in infilenames_IV.keys():
            continue
        if keyB2 not in infilenames_IV.keys():
            continue
        
        graphs = {}
        graphs['A'] = TGraphRatio(g_I_OVeff_avgLast[keyA1],g_I_OVeff_avgLast[keyA2])
        graphs['B'] = TGraphRatio(g_I_OVeff_avgLast[keyB1],g_I_OVeff_avgLast[keyB2])
        
        c = ROOT.TCanvas('c_IV_%s_%dC_vendorRatio'%(fluence,temp),'c_IV_%s_%dC_vendorRatio'%(fluence,temp),1600,700)
        c.Divide(2,1)
        c.cd(1)
        hPad1 = ROOT.gPad.DrawFrame(canvasXmin,canvasYmin,canvasXmax,canvasYmax)
        hPad1.SetTitle(";V_{ov}^{eff} [V]; I_{array} / 16 [mA]")
        hPad1.Draw()
        ROOT.gPad.SetGridx()
        ROOT.gPad.SetGridy()
        c.cd(2)
        hPad2 = ROOT.gPad.DrawFrame(canvasXmin,canvasYmin2,canvasXmax,canvasYmax2)
        hPad2.SetTitle(";V_{ov}^{eff} [V]; I_{FBK} / I_{HPK}")
        hPad2.Draw()
        ROOT.gPad.SetGridx()
        ROOT.gPad.SetGridy()
        
        label = ROOT.TLatex(0.17,0.94,'%s %dC'%(fluence,temp))
        label.SetNDC()
        label.SetTextSize(0.04)
        label.SetTextFont(42)
        
        c.cd(1)
        g_I_OVeff_avgLast[keyA1].SetMarkerColor(ROOT.kRed)
        g_I_OVeff_avgLast[keyA1].SetLineColor(ROOT.kRed)
        g_I_OVeff_avgLast[keyA1].SetMarkerSize(1.)
        g_I_OVeff_avgLast[keyA1].SetMarkerStyle(20)
        g_I_OVeff_avgLast[keyA1].Draw('PEL,same')
        g_I_OVeff_avgLast[keyB1].SetMarkerColor(ROOT.kGreen)
        g_I_OVeff_avgLast[keyB1].SetLineColor(ROOT.kGreen)
        g_I_OVeff_avgLast[keyB1].SetMarkerSize(1.)
        g_I_OVeff_avgLast[keyB1].SetMarkerStyle(20)
        g_I_OVeff_avgLast[keyB1].Draw('PEL,same')
        g_I_OVeff_avgLast[keyA2].SetMarkerColor(ROOT.kRed)
        g_I_OVeff_avgLast[keyA2].SetLineColor(ROOT.kRed)
        g_I_OVeff_avgLast[keyA2].SetMarkerSize(1.)
        g_I_OVeff_avgLast[keyA2].SetMarkerStyle(24)
        g_I_OVeff_avgLast[keyA2].Draw('PEL,same')
        g_I_OVeff_avgLast[keyB2].SetMarkerColor(ROOT.kGreen)
        g_I_OVeff_avgLast[keyB2].SetLineColor(ROOT.kGreen)
        g_I_OVeff_avgLast[keyB2].SetMarkerSize(1.)
        g_I_OVeff_avgLast[keyB2].SetMarkerStyle(24)
        g_I_OVeff_avgLast[keyB2].Draw('PEL,same')
        label.Draw('same')

        c.cd(2)
        graphs['A'].SetLineWidth(4)
        graphs['A'].SetLineColor(ROOT.kRed)
        graphs['A'].Draw('L,same')
        graphs['B'].SetLineWidth(4)
        graphs['B'].SetLineColor(ROOT.kGreen)
        graphs['B'].Draw('L,same')
        label.Draw('same')
        
        c.Print(outdir+'%s.png'%c.GetName())
