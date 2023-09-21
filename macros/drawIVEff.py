#! /usr/bin/python

import ROOT
import glob
import math
import argparse
import json
import os
import numpy
from collections import OrderedDict 

from SiPM import *

aldos = ['A', 'B']

channelMap = {}
channelMap[(0,'A')] = 2
channelMap[(0,'B')] = 1
channelMap[(2,'A')] = 6
channelMap[(2,'B')] = 5

inputfolder = '/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2b_gen23/'



parser = argparse.ArgumentParser(description='draw IV scans from stress test')
parser.add_argument("--asic", type=int, required=True, help="asic") 
args = parser.parse_args()


################################
### edit here
################################

confs = [ 
#        '24.00', '24.03', '24.01', '24.04',
#        '25.00', '25.01', '25.02',
#        '27.00', '27.01', '27.02', '27.03',
#        '28.00', '28.01', '28.02', '28.03',
#        '29.00', '29.01',
        '36.00', '36.01', #TOFHIR2C
        '37.00',          #TOFHIR2C
        '38.00',          #TOFHIR2C
        '39.00',          #TOFHIR2C
]

temps = {}
temps['24.04'] = -30.
temps['24.03'] = -35.
temps['24.00'] = -40.
temps['24.01'] = -45.
temps['25.00'] = -40.
temps['25.01'] = -35.
temps['25.02'] = -30.
temps['27.00'] = -32.
temps['27.01'] = -37.
temps['27.02'] = -27.
temps['27.03'] = -22.
temps['28.00'] = -19.
temps['28.01'] = -32.
temps['28.02'] = 0.
temps['28.03'] = 12.
temps['29.00'] = -32.
temps['29.01'] = -22.
temps['36.00'] = -35.
temps['36.01'] = -30.
temps['37.00'] = -35.
temps['38.00'] = -30.
temps['39.00'] = -30.

sipmTypes = {}
sipmTypes['24.04'] = 'HPK-PIT-C25-ES2'
sipmTypes['24.03'] = 'HPK-PIT-C25-ES2'
sipmTypes['24.00'] = 'HPK-PIT-C25-ES2'
sipmTypes['24.01'] = 'HPK-PIT-C25-ES2'
sipmTypes['25.00'] = 'HPK-PIT-C20-ES2'
sipmTypes['25.01'] = 'HPK-PIT-C20-ES2'
sipmTypes['25.02'] = 'HPK-PIT-C20-ES2'
sipmTypes['27.00'] = 'HPK-PIT-C25-ES2'
sipmTypes['27.01'] = 'HPK-PIT-C25-ES2'
sipmTypes['27.02'] = 'HPK-PIT-C25-ES2'
sipmTypes['27.03'] = 'HPK-PIT-C25-ES2'
sipmTypes['28.00'] = 'HPK-PIT-C25-ES2'
sipmTypes['28.01'] = 'HPK-PIT-C25-ES2'
sipmTypes['28.02'] = 'HPK-PIT-C25-ES2'
sipmTypes['28.03'] = 'HPK-PIT-C25-ES2'
sipmTypes['29.00'] = 'HPK-PIT-C25-ES2'
sipmTypes['29.01'] = 'HPK-PIT-C25-ES2'
sipmTypes['36.00'] = 'HPK-PIT-C25-ES2'
sipmTypes['36.01'] = 'HPK-PIT-C25-ES2'
sipmTypes['37.00'] = 'HPK-PIT-C20-ES2'
sipmTypes['38.00'] = 'HPK-MS' #15 um 
sipmTypes['39.00'] = 'HPK-PIT-C25-ES2'

labels = {}
labels['24.04'] = 'HPK_2E14_LYSO815'
labels['24.03'] = 'HPK_2E14_LYSO815'
labels['24.00'] = 'HPK_2E14_LYSO815'
labels['24.01'] = 'HPK_2E14_LYSO815'
labels['25.00'] = 'HPK_2E14_LYSO825'
labels['25.01'] = 'HPK_2E14_LYSO825'
labels['25.02'] = 'HPK_2E14_LYSO825'
labels['27.00'] = 'HPK_1E14_LYSO819'
labels['27.01'] = 'HPK_1E14_LYSO819'
labels['27.02'] = 'HPK_1E14_LYSO819'
labels['27.03'] = 'HPK_1E14_LYSO819'
labels['28.00'] = 'HPK_1E13_LYSO829'
labels['28.01'] = 'HPK_1E13_LYSO829'
labels['28.02'] = 'HPK_1E13_LYSO829'
labels['28.03'] = 'HPK_1E13_LYSO829'
labels['29.00'] = 'HPK_1E14_LYSO817'
labels['29.01'] = 'HPK_1E14_LYSO817'  
labels['36.00'] = 'HPK_2E14_LYSO815'
labels['36.01'] = 'HPK_2E14_LYSO815'
labels['37.00'] = 'HPK_2E14_LYSO825'
labels['38.00'] = 'HPK_1E14_LYSO844'
labels['39.00'] = 'HPK_1E14_LYSO819'


gainDrops = {}
gainDrops['24.04'] = 0.08
gainDrops['24.03'] = 0.08
gainDrops['24.00'] = 0.08
gainDrops['24.01'] = 0.08
gainDrops['25.00'] = 0.08
gainDrops['25.01'] = 0.08
gainDrops['25.02'] = 0.08
gainDrops['27.00'] = 0.04
gainDrops['27.01'] = 0.04
gainDrops['27.02'] = 0.04
gainDrops['27.03'] = 0.04
gainDrops['28.00'] = 0.
gainDrops['28.01'] = 0.
gainDrops['28.02'] = 0.
gainDrops['28.03'] = 0.
gainDrops['29.00'] = 0.04
gainDrops['29.01'] = 0.04
gainDrops['36.00'] = 0.08
gainDrops['36.01'] = 0.08
gainDrops['37.00'] = 0.08
gainDrops['38.00'] = 0.04
gainDrops['39.00'] = 0.04




################################


graphs_IVarray = {}
graphs_IVch = {}
graphs_DCR = {}
for aldo in aldos:
        graphs_IVarray[(args.asic,aldo)] = ROOT.TGraph()
        graphs_IVch[(args.asic,aldo)] = ROOT.TGraph()
        graphs_DCR[(args.asic,aldo)] = ROOT.TGraph()


#outfile = ROOT.TFile('./logIVEff.root','RECREATE')
outfile = ROOT.TFile('./logIVEff_TOFHIR2C.root','RECREATE')

output_dict = OrderedDict()

for conf in confs:
        print('\n'+labels[conf]+' - '+'T = %s'%temps[conf])
        for aldo in aldos:
                graphs_IVarray[(conf,args.asic,aldo)] = ROOT.TGraph()
                graphs_IVch[(conf,args.asic,aldo)] = ROOT.TGraph()
                graphs_DCR[(conf,args.asic,aldo)] = ROOT.TGraph()
        
        for aldo in aldos:
            print('>>> ALDO'+aldo)
        
            infilenames = glob.glob('%s/logs_%s/logIV_*_ASIC%d_ALDO%s_ch%d_*.root'%(inputfolder,conf,args.asic,aldo,channelMap[(args.asic,aldo)]))
            print infilenames
            
            VbiasList = []
            
            for infilename in infilenames:
                # get Vbr from bias_settings_aldo.tsv
                infilealdo = open('%s/config_%s/bias_settings_aldo.tsv'%(inputfolder,conf),'r')
                Vbr = 0
                for line in infilealdo:
                        tokens = line.split()
                        if tokens[2] == str(args.asic) and tokens[3] == aldo:
                                Vbr = tokens[4]
                print("Vbr: "+str(Vbr))
                
                infile = ROOT.TFile(infilename,"READ")
                graph = infile.Get('g_IV')
                for point in range(graph.GetN()):
                        x = graph.GetPointX(point)
                        y = graph.GetPointY(point)
                        
                        VovEff = x - float(Vbr) - y*1E-06*25.
                        DCR = y*1E-6 / (1.602E-19 * Gain(sipmTypes[conf],VovEff) * (1.-gainDrops[conf]) ) * 1E-09 / 16.
                        #print(round(x,2),round(x-float(Vbr),2),round(VovEff,2),round(y,2),round(y/16.,2),round(Gain(sipmTypes[conf],VovEff),0),round(DCR,1))
                        
                        graphs_IVarray[(conf,args.asic,aldo)].SetPoint(graphs_IVarray[(conf,args.asic,aldo)].GetN(),VovEff,y)
                        graphs_IVch[(conf,args.asic,aldo)].SetPoint(graphs_IVch[(conf,args.asic,aldo)].GetN(),VovEff,y/16.)
                        graphs_DCR[(conf,args.asic,aldo)].SetPoint(graphs_DCR[(conf,args.asic,aldo)].GetN(),VovEff,DCR)
                
                
                output_dict['%s_T%dC_%s'%(labels[conf],temps[conf],aldo)] = OrderedDict()
                for VovSet in numpy.arange(0.5, 3.1, 0.05):
                        current = graph.Eval(float(Vbr)+VovSet)
                        VovEff = VovSet - current*1E-06*25.
                        DCR = current*1E-6 / (1.602E-19 * Gain(sipmTypes[conf],VovEff) * (1.-gainDrops[conf]) ) * 1E-09 / 16.
                        #output_dict['%s_T%dC_%s'%(labels[conf],temps[conf],aldo)]['%.2f'%(VovSet)] = ['%.2f'%round(VovEff,2),'%.1f'%round(DCR,1)]
                        output_dict['%s_T%dC_%s'%(labels[conf],temps[conf],aldo)]['%.2f'%(VovSet)] = ['%.2f'%round(VovEff,2),'%.1f'%round(DCR,1), '%.2f'%round(current/16*1E-3,2)]

        
        for aldo in aldos:
                outfile.cd()

                graphs_IVarray[(conf,args.asic,aldo)].Sort()
                graphs_IVch[(conf,args.asic,aldo)].Sort()
                graphs_DCR[(conf,args.asic,aldo)].Sort()

                graphs_IVarray[(conf,args.asic,aldo)].Write('g_IVEff_array_%s_T%.0f_ALDO%s'%(labels[conf],temps[conf],aldo))
                graphs_IVch[(conf,args.asic,aldo)].Write('g_IVEff_ch_%s_T%.0f_ALDO%s'%(labels[conf],temps[conf],aldo))
                graphs_DCR[(conf,args.asic,aldo)].Write('g_DCR_%s_T%.0f_ALDO%s'%(labels[conf],temps[conf],aldo))
                
                
outfile.Close()

#with open('VovsEff.json', 'w') as fp:
with open('VovsEff_TOFHIR2C.json', 'w') as fp:
        json.dump(output_dict, fp, indent=2)


'''
    c1 = ROOT.TCanvas('c1_ASIC%d_ALDO%s'%(args.asic,aldo),'c1_ASIC%d_ALDO%s'%(args.asic,aldo),1400,700)
    c1.Divide(2,1)
    c1.cd(1)
    hPad1 = ROOT.gPad.DrawFrame(xMin-0.1*(xMax-xMin),0.,xMax+0.1*(xMax-xMin),1.1*yMax)
    hPad1.SetTitle(";V_{bias} [V]; I [#muA]")
    hPad1.Draw() 
    graph_ave.SetMarkerStyle(20)
    graph_ave.SetMarkerSize(0.7)
    graph_ave.Draw("PL,same")
    
    c1.cd(2)
    hPad2 = ROOT.gPad.DrawFrame(xMin-0.1*(xMax-xMin),0.,xMax+0.1*(xMax-xMin),15.)
    hPad2.SetTitle(";V_{bias} [V]; #DeltalogI/#DeltaV [#muA/V]")
    hPad2.Draw()
    
    graph_dlogIdV = ROOT.TGraph()
    for point in range(1,graph_ave.GetN()):
        x1 = graph.GetPointX(point-1)
        y1 = graph.GetPointY(point-1)
        x2 = graph.GetPointX(point)
        y2 = graph.GetPointY(point)
        print(x1,y1,x2,y2)
        graph_dlogIdV.SetPoint(graph_dlogIdV.GetN(),0.5*(x1+x2),(math.log(y2)-math.log(y1))/(x2-x1))
    graph_dlogIdV.SetMarkerStyle(20)
    graph_dlogIdV.SetMarkerSize(0.7)
    graph_dlogIdV.Draw("PL,same")
    
    graph_dlogIdV_ave = ROOT.TGraph()
    for point in range(1,graph_dlogIdV.GetN()-1):
        x = graph_dlogIdV.GetPointX(point)
        y1 = graph_dlogIdV.GetPointY(point-1)
        y2 = graph_dlogIdV.GetPointY(point)
        y3 = graph_dlogIdV.GetPointY(point+1)
        graph_dlogIdV_ave.SetPoint(graph_dlogIdV_ave.GetN(),x,(y1+y2+y3)/3.)
    graph_dlogIdV_ave.SetLineColor(ROOT.kTeal)
    graph_dlogIdV_ave.SetLineWidth(2)
    graph_dlogIdV_ave.Draw("L,same")
    
    maximum = FindMaximumPoint(graph_dlogIdV_ave)
    fitFunc = ROOT.TF1("fitFunc","gaus(0)",xMin,xMax)
    fitFunc.SetParameter(1,maximum)
    graph_dlogIdV_ave.Fit(fitFunc,"QNRS+")
    
    xMin = maximum - fitFunc.GetParameter(2)
    xMax = maximum + fitFunc.GetParameter(2)
    fitFunc2 = ROOT.TF1("fitFunc2","gaus(0)",xMin,xMax)
    graph_dlogIdV_ave.Fit(fitFunc2,"QNRS+")
    fitFunc2.SetLineColor(ROOT.kRed)
    fitFunc2.SetLineWidth(1)
    fitFunc2.Draw("same")
    
    Vbr = fitFunc2.GetParameter(1)
    latex = ROOT.TLatex(Vbr,fitFunc2.Eval(Vbr),'V_{br.} = %.2f V'%Vbr)
    latex.SetTextFont(42)
    latex.SetTextSize(0.04)
    latex.SetTextColor(ROOT.kRed)
    latex.Draw("same")
    
    print('ASIC %d, ALDO %s:   Vbr = %.2f'%(args.asic,aldo,Vbr))
    c1.Print('../logs/IV_%s_ASIC%d_ALDO%s.png'%(args.label, args.asic,aldo))
    
    key='0       0       %d       %s'%(args.asic,aldo)
    VbrString='%.2f'%Vbr
    command = 'sed -i \"s%^'+key+'.*$%'+key+'        '+VbrString+'           5.00%\"'+' ../config/bias_settings_aldo.tsv'
    #print(command)
    os.system(command)
'''
