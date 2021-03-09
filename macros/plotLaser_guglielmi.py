#! /usr/bin/env python
import os
import shutil
import glob
import math
import array
import sys
import time
import numpy as np


import ROOT
from ROOT import *
import CMS_lumi, tdrstyle

from ROOT import  TLatex

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
ROOT.gStyle.SetImageScaling(2.5);



def GetLinearizedEnergy(x, par):
  xx = x[0]
  if( xx < par[0] ): return ( par[1]*(math.exp(par[2]*xx)-1) + par[3] );
  else:              return ( par[1]*(math.exp(par[2]*par[0])-1)+par[3])/(math.exp(par[4]*par[0]))*math.exp(par[4]*xx );


def GetSaturationCorrection(Ncells, Edep, LY, PDE, LCE):
  Npe    = Edep * LY * PDE * LCE
  Nfired = Ncells * (1 - math.exp(-Npe/Ncells))
  k = Nfired/Npe
  #print Edep, LY, PDE, LCE, ECF
  #print 'Ncells = %d, Npe = %f, Nfired = %f,  k = %f'%(Ncells, Npe, Nfired, k)
  return k





#to run python plotLaser_guglielmi.py /home/cmsdaq/guglielmi/Lab5015Analysis/ /data/ALIO/plots/ /home/cmsdaq/guglielmi/Lab5015Analysis/plots/

###
outdir = '/var/www/html/ModuleCharacterization/guglielmi/Laser/new_funcTot/'
path = os.path.join(outdir)
#os.mkdir(path)
plotdir = sys.argv[1] + 'plots/'
step1dir = sys.argv[2]
otherdir = sys.argv[3] 

Vov = 3.0
LY = 40000. # gamma/MeV
LCE = 0.15
Ncells = 40000
PDE = 0.351



#laserTunes = [10,20,30,40,50,60,70,80,83,85,86,87,88,89,90]
laserTunes = [20,30,40,50,60,70,75,80,85,87,90,91,92,93,94]
run =[2665, 2664, 2663, 2662, 2661, 2660, 2659, 2658, 2655, 2657, 2654, 2653, 2656, 2652, 2651]

f_3_20_Na22 = ROOT.TFile.Open('%smoduleCharacterization_step1_run2611.root'%step1dir)
f_5_20_Na22 = ROOT.TFile.Open('%smoduleCharacterization_step1_run2620.root'%step1dir)
f_7_20_Na22 = ROOT.TFile.Open('%smoduleCharacterization_step1_run2629.root'%step1dir)

h_en_Vov3_th20_Na22 = f_3_20_Na22.Get('h1_energy_bar00L-R_Vov3.0_th20')
h_en_Vov5_th20_Na22 = f_5_20_Na22.Get('h1_energy_bar00L-R_Vov5.0_th20')
h_en_Vov7_th20_Na22 = f_7_20_Na22.Get('h1_energy_bar00L-R_Vov7.0_th20')

f_3_20_Co60 = ROOT.TFile.Open('%smoduleCharacterization_step1_run2670.root'%step1dir)
f_5_20_Co60 = ROOT.TFile.Open('%smoduleCharacterization_step1_run2671.root'%step1dir)
f_7_20_Co60 = ROOT.TFile.Open('%smoduleCharacterization_step1_run2672.root'%step1dir)

h_en_Vov3_th20_Co60 = f_3_20_Co60.Get('h1_energy_bar00L-R_Vov3.0_th20')
h_en_Vov5_th20_Co60 = f_5_20_Co60.Get('h1_energy_bar00L-R_Vov5.0_th20')
h_en_Vov7_th20_Co60 = f_7_20_Co60.Get('h1_energy_bar00L-R_Vov7.0_th20')



energy = {}
linearizedEnergy = {}
f = {}
fstep2 = {}
h1_energy  = {}
g_tRes_vs_threshold = {}
g_enPeak_vs_bar_Vov3_th20 = {}

g_tRes_vs_energy = ROOT.TGraphErrors()
g_tRes_vs_linearizedEnergy = ROOT.TGraphErrors()
g_tRes_vs_energy_bestTh = ROOT.TGraphErrors()
g_tRes_vs_linearizedEnergy_bestTh = ROOT.TGraphErrors()
g_energy_vs_laserTune = ROOT.TGraphErrors()

tRes_err_sys = 2.
energy_err_sys = 0.015 # 1.5%


# outfile root
outFileName = ('%splotLaser.root' %plotdir)
outFile = ROOT.TFile.Open(str(outFileName),"RECREATE");




# from Andrea
funcLinearizedEnergy = ROOT.TF1('funcLinearizedEnergy',GetLinearizedEnergy, 0, 100,5)
funcLinearizedEnergy.SetParameter(0,1.79617e+01)
funcLinearizedEnergy.SetParameter(1,7.07013e+00)
funcLinearizedEnergy.SetParameter(2,1.87292e-02)
funcLinearizedEnergy.SetParameter(3,-1.70555e-02)
funcLinearizedEnergy.SetParameter(4,6.59613e-02)

funcLinearizedEnergyCo60 = ROOT.TF1('funcLinearizedEnergyCo60',GetLinearizedEnergy, 0, 100,5)
funcLinearizedEnergyCo60.SetParameter(0,1.79617e+01)
funcLinearizedEnergyCo60.SetParameter(1,7.07013e+00)
funcLinearizedEnergyCo60.SetParameter(2,1.87292e-02)
funcLinearizedEnergyCo60.SetParameter(3,-1.70555e-02)
funcLinearizedEnergyCo60.SetParameter(4,6.59613e-02)






#Linearization
f1 = ROOT.TFile.Open('%smoduleCharacterization_single_HPK_HDR2_step3.root'%plotdir)
f2 = ROOT.TFile.Open('%smoduleCharacterization_single_HPK_HDR2_Co60_step3.root'%plotdir)

'''funcLinearizedEnergy = TF1('f_linearity_Vov%.01f' %Vov, '[0]*(exp([1]*x)-1) + [2]' ,0.,200.)
funcLinearizedEnergy.SetParameters(6.92,0.0158,-0.0639)
funcLinearizedEnergyCo60 = TF1('f_linearity_Vov%.01f_Co60' %Vov, '[0]*(exp([1]*x)-1) + [2]' ,0.,200.)
funcLinearizedEnergyCo60.SetParameters(6.92,0.0158,-0.0639)'''



g_Linearization_Na22 = ROOT.TGraphErrors()
g_Linearization_Na22 = f1.Get('g_linearity_Vov%.01f' %Vov)
g_Linearization_Co60 = ROOT.TGraphErrors()
g_Linearization_Co60 = f2.Get('g_linearity_Vov%.01f' %Vov)
g_Vov3_Na22 = ROOT.TGraphErrors()
g_Vov3_Na22 = f1.Get('g_ov3_Vov%.01f' %Vov)
g_Vov5_Na22 = ROOT.TGraphErrors()
g_Vov5_Na22 = f1.Get('g_ov5_Vov%.01f' %Vov)
g_Vov7_Na22 = ROOT.TGraphErrors()
g_Vov7_Na22 = f1.Get('g_ov7_Vov%.01f' %Vov)
g_Vov3_Co60 = ROOT.TGraphErrors()
g_Vov3_Co60 = f2.Get('g_ov3_Vov%.01f' %Vov)
g_Vov5_Co60 = ROOT.TGraphErrors()
g_Vov5_Co60 = f2.Get('g_ov5_Vov%.01f' %Vov)
g_Vov7_Co60 = ROOT.TGraphErrors()
g_Vov7_Co60 = f2.Get('g_ov7_Vov%.01f' %Vov)



g_Linearization_Co60_rescaled = ROOT.TGraphErrors()
g_Vov3_Co60_rescaled = ROOT.TGraphErrors()
g_Vov5_Co60_rescaled = ROOT.TGraphErrors()
g_Vov7_Co60_rescaled = ROOT.TGraphErrors()
for i in range(g_Linearization_Co60.GetN()):
	g_Linearization_Co60_rescaled.SetPoint(i,g_Linearization_Co60.GetPointX(i)/0.93, g_Linearization_Co60.GetPointY(i))
	#print(g_Linearization_Co60.GetPointX(i)/0.93)
for i in range(g_Vov3_Co60.GetN()):
	g_Vov3_Co60_rescaled.SetPoint(i,g_Vov3_Co60.GetPointX(i)/0.93, g_Vov3_Co60.GetPointY(i))
	g_Vov5_Co60_rescaled.SetPoint(i,g_Vov5_Co60.GetPointX(i)/0.93, g_Vov5_Co60.GetPointY(i))
	g_Vov7_Co60_rescaled.SetPoint(i,g_Vov7_Co60.GetPointX(i)/0.93, g_Vov7_Co60.GetPointY(i))

'''g_Vov3_Co60_rescaled.SetPoint(0,7.57231769017414, g_Vov3_Co60.GetPointY(0))	
g_Vov3_Co60_rescaled.SetPoint(1,8.59874438475017, g_Vov3_Co60.GetPointY(1))
g_Vov3_Co60_rescaled.SetPoint(2,16.1710620749243, g_Vov3_Co60.GetPointY(2))
#g_Vov3_Co60_rescaled.SetPoint(2,g_Vov3_Co60.GetPointX(2), g_Vov3_Co60.GetPointY(2))
g_Linearization_Co60_rescaled.SetPoint(0,7.57231769017414, g_Vov3_Co60.GetPointY(0))
g_Linearization_Co60_rescaled.SetPoint(1,8.59874438475017, g_Vov3_Co60.GetPointY(1))
g_Linearization_Co60_rescaled.SetPoint(2,16.1710620749243, g_Vov3_Co60.GetPointY(2))
#g_Linearization_Co60_rescaled.SetPoint(2,g_Vov3_Co60.GetPointX(2), g_Vov3_Co60.GetPointY(2))

g_Vov5_Co60_rescaled.SetPoint(0,13.2141438415333, g_Vov5_Co60.GetPointY(0))
g_Vov5_Co60_rescaled.SetPoint(1,15.0053193494649, g_Vov5_Co60.GetPointY(1))
g_Vov5_Co60_rescaled.SetPoint(2,28.2194631909982, g_Vov5_Co60.GetPointY(2))
#g_Vov5_Co60_rescaled.SetPoint(2,g_Vov5_Co60.GetPointX(2), g_Vov5_Co60.GetPointY(2))
g_Linearization_Co60_rescaled.SetPoint(3,13.2141438415333, g_Vov5_Co60.GetPointY(0))
g_Linearization_Co60_rescaled.SetPoint(4,15.0053193494649, g_Vov5_Co60.GetPointY(1))
g_Linearization_Co60_rescaled.SetPoint(5,28.2194631909982, g_Vov5_Co60.GetPointY(2))
g_Linearization_Co60_rescaled.SetPoint(5,g_Vov5_Co60.GetPointX(2), g_Vov5_Co60.GetPointY(2))

g_Vov7_Co60_rescaled.SetPoint(0,18.033908680328, g_Vov7_Co60.GetPointY(0))
g_Vov7_Co60_rescaled.SetPoint(1,20.4784026958201, g_Vov7_Co60.GetPointY(1))
g_Vov7_Co60_rescaled.SetPoint(2,38.5123113761481, g_Vov7_Co60.GetPointY(2))
#g_Vov7_Co60_rescaled.SetPoint(2,g_Vov7_Co60.GetPointX(2), g_Vov7_Co60.GetPointY(2))
g_Linearization_Co60_rescaled.SetPoint(6,18.033908680328, g_Vov7_Co60.GetPointY(0))
g_Linearization_Co60_rescaled.SetPoint(7,20.4784026958201, g_Vov7_Co60.GetPointY(1))
g_Linearization_Co60_rescaled.SetPoint(8,38.5123113761481, g_Vov7_Co60.GetPointY(2))
#g_Linearization_Co60_rescaled.SetPoint(8,g_Vov7_Co60.GetPointX(2), g_Vov7_Co60.GetPointY(2))'''


c_colour_rescaled = ROOT.TCanvas('Lin coloured scaled','Lin coloured scaled',700,500)
hPad_lin10 = ROOT.gPad.DrawFrame(0.,0.,30.,10.)
hPad_lin10.GetYaxis().SetTitleSize(0.04)
hPad_lin10.GetXaxis().SetTitle('TOFPET2 amp (a.u.)')
hPad_lin10.GetYaxis().SetTitle('E_{MeV} #times G #times PDE #times ECF #times k_{sat}/[G #times PDE #times ECF]_{%.1f V}' %Vov)
hPad_lin10.Draw()
gPad.SetGridx()
gPad.SetGridy()
g_Vov3_Na22.SetMarkerStyle(20)
g_Vov3_Na22.SetMarkerSize(1.)
g_Vov3_Na22.SetMarkerColor(kGreen)
g_Vov3_Na22.Draw("P")
g_Vov5_Na22.SetMarkerStyle(20)
g_Vov5_Na22.SetMarkerSize(1.)
g_Vov5_Na22.SetMarkerColor(kRed)
g_Vov5_Na22.Draw("Psame")
g_Vov7_Na22.SetMarkerStyle(20)
g_Vov7_Na22.SetMarkerSize(1.)
g_Vov7_Na22.SetMarkerColor(kBlue)
g_Vov7_Na22.Draw("Psame")
g_Vov3_Co60_rescaled.SetMarkerStyle(20)
g_Vov3_Co60_rescaled.SetMarkerSize(1.)
g_Vov3_Co60_rescaled.SetMarkerColor(kYellow)
g_Vov3_Co60_rescaled.Draw("Psame")
g_Vov5_Co60_rescaled.SetMarkerStyle(20)
g_Vov5_Co60_rescaled.SetMarkerSize(1.)
g_Vov5_Co60_rescaled.SetMarkerColor(kViolet)
g_Vov5_Co60_rescaled.Draw("Psame")
g_Vov7_Co60_rescaled.SetMarkerStyle(20)
g_Vov7_Co60_rescaled.SetMarkerSize(1.)
g_Vov7_Co60_rescaled.SetMarkerColor(kBlack)
g_Vov7_Co60_rescaled.Draw("Psame")
'''g_Linearization_Co60_rescaled.SetMarkerStyle(20)
g_Linearization_Co60_rescaled.SetMarkerSize(1.)
g_Linearization_Co60_rescaled.SetMarkerColor(kYellow)
g_Linearization_Co60_rescaled.Draw("Psame")'''
leg_rescaled = ROOT.TLegend(0.33,0.73,0.53,0.85)
leg_rescaled.SetFillStyle(0)
leg_rescaled.SetBorderSize(0)
leg_rescaled.AddEntry(g_Vov3_Na22,'Na22 Vov3','PL')
leg_rescaled.AddEntry(g_Vov5_Na22,'Na22 Vov5','PL')
leg_rescaled.AddEntry(g_Vov7_Na22,'Na22 Vov7','PL')
leg_rescaled.AddEntry(g_Vov3_Co60_rescaled,'Co60 Vov3','PL')
leg_rescaled.AddEntry(g_Vov5_Co60_rescaled,'Co60 Vov5','PL')
leg_rescaled.AddEntry(g_Vov7_Co60_rescaled,'Co60 Vov7','PL')
leg_rescaled.Draw('same')




x = []
y = []
n = g_Linearization_Na22.GetN() + g_Linearization_Co60.GetN() 
for i in range( n ):
	if i<g_Linearization_Na22.GetN():
		x.append( g_Linearization_Na22.GetPointX(i))
		y.append( g_Linearization_Na22.GetPointY(i))
	else:
		x.append( g_Linearization_Co60_rescaled.GetPointX(i-g_Linearization_Na22.GetN()))
		y.append( g_Linearization_Co60_rescaled.GetPointY(i-g_Linearization_Na22.GetN()))
		'''print (x)
		print (y)
		print (i)
		print (i-g_Linearization_Na22.GetN())'''


g_Linearization = ROOT.TGraph( n, np.array(x), np.array(y) )
	
c_lin = ROOT.TCanvas('Linearization','Linearization',1300,600)
c_lin.Divide(2,2)	
'''c_lin.cd(1)
hPad_lin = ROOT.gPad.DrawFrame(0.,0.,35.,10.)
hPad_lin.GetYaxis().SetTitleSize(0.04)
hPad_lin.GetXaxis().SetTitle('TOFPET2 amp (a.u.)')
hPad_lin.GetYaxis().SetTitle('E_{MeV} #times G #times PDE #times ECF #times k_{sat}/[G #times PDE #times ECF]_{%.1f V}' %Vov)
hPad_lin.Draw()
gPad.SetGridx()
gPad.SetGridy()
g_Linearization_Na22.SetMarkerStyle(20)
g_Linearization_Na22.SetMarkerSize(1.)
g_Linearization_Na22.SetMarkerColor(kRed)
g_Linearization_Na22.Draw("P")
g_Linearization_Na22.Fit(funcLinearizedEnergy, "", "", 0, 35)
c_lin.cd(2)
hPad_lin1 = ROOT.gPad.DrawFrame(0.,0.,35.,10.)
hPad_lin1.GetYaxis().SetTitleSize(0.04)
hPad_lin1.GetXaxis().SetTitle('TOFPET2 amp (a.u.)')
hPad_lin1.GetYaxis().SetTitle('E_{MeV} #times G #times PDE #times ECF #times k_{sat}/[G #times PDE #times ECF]_{%.1f V}' %Vov)
hPad_lin1.Draw()
gPad.SetGridx()
gPad.SetGridy()
g_Linearization_Co60_rescaled.SetMarkerStyle(20)
g_Linearization_Co60_rescaled.SetMarkerSize(1.)
g_Linearization_Co60_rescaled.SetMarkerColor(kYellow)
g_Linearization_Co60_rescaled.Draw("P")
g_Linearization_Co60_rescaled.Fit(funcLinearizedEnergy, "", "", 0, 35)

c_lin.cd(3)
hPad_lin2 = ROOT.gPad.DrawFrame(0.,0.,35.,10.)
hPad_lin2.GetYaxis().SetTitleSize(0.04)
hPad_lin2.GetXaxis().SetTitle('TOFPET2 amp (a.u.)')
hPad_lin2.GetYaxis().SetTitle('E_{MeV} #times G #times PDE #times ECF #times k_{sat}/[G #times PDE #times ECF]_{%.1f V}' %Vov)
hPad_lin2.Draw()
gPad.SetGridx()
gPad.SetGridy()
g_Linearization_Na22.SetMarkerStyle(20)
g_Linearization_Na22.SetMarkerSize(1.)
g_Linearization_Na22.SetMarkerColor(kRed)
g_Linearization_Na22.Draw("P")
g_Linearization_Co60_rescaled.SetMarkerStyle(20)
g_Linearization_Co60_rescaled.SetMarkerSize(1.)
g_Linearization_Co60_rescaled.SetMarkerColor(kYellow)
g_Linearization_Co60_rescaled.Draw("P,same")'''
c_lin.cd(4)
hPad_lin3 = ROOT.gPad.DrawFrame(0.,0.,35.,10.)
hPad_lin3.GetYaxis().SetTitleSize(0.04)
hPad_lin3.GetXaxis().SetTitle('TOFPET2 amp (a.u.)')
hPad_lin3.GetYaxis().SetTitle('E_{MeV} #times G #times PDE #times ECF #times k_{sat}/[G #times PDE #times ECF]_{%.1f V}' %Vov)
hPad_lin3.Draw()
gPad.SetGridx()
gPad.SetGridy()
g_Linearization.SetMarkerStyle(20)
g_Linearization.SetMarkerSize(1.)
g_Linearization.Draw("P")
funcLinearizedEnergy.SetLineColor(kRed)
g_Linearization.Fit(funcLinearizedEnergy, "", "", 0, 35)


	

c_colour = ROOT.TCanvas('Lin coloured','Lin coloured',700,500)
hPad_lin9 = ROOT.gPad.DrawFrame(0.,0.,30.,10.)
hPad_lin9.GetYaxis().SetTitleSize(0.04)
hPad_lin9.GetXaxis().SetTitle('TOFPET2 amp (a.u.)')
hPad_lin9.GetYaxis().SetTitle('E_{MeV} #times G #times PDE #times ECF #times k_{sat}/[G #times PDE #times ECF]_{%.1f V}' %Vov)
hPad_lin9.Draw()
gPad.SetGridx()
gPad.SetGridy()
g_Vov3_Na22.SetMarkerStyle(20)
g_Vov3_Na22.SetMarkerSize(1.)
g_Vov3_Na22.SetMarkerColor(kGreen)
g_Vov3_Na22.Draw("P")
g_Vov5_Na22.SetMarkerStyle(20)
g_Vov5_Na22.SetMarkerSize(1.)
g_Vov5_Na22.SetMarkerColor(kRed)
g_Vov5_Na22.Draw("P")
g_Vov7_Na22.SetMarkerStyle(20)
g_Vov7_Na22.SetMarkerSize(1.)
g_Vov7_Na22.SetMarkerColor(kBlue)
g_Vov7_Na22.Draw("P")
g_Vov3_Co60.SetMarkerStyle(20)
g_Vov3_Co60.SetMarkerSize(1.)
g_Vov3_Co60.SetMarkerColor(kYellow)
g_Vov3_Co60.Draw("P,same")
g_Vov5_Co60.SetMarkerStyle(20)
g_Vov5_Co60.SetMarkerSize(1.)
g_Vov5_Co60.SetMarkerColor(kViolet)
g_Vov5_Co60.Draw("P,same")
g_Vov7_Co60.SetMarkerStyle(20)
g_Vov7_Co60.SetMarkerSize(1.)
g_Vov7_Co60.SetMarkerColor(kBlack)
g_Vov7_Co60.Draw("P,same")
leg_colour = ROOT.TLegend(0.33,0.73,0.53,0.85)
leg_colour.SetFillStyle(0)
leg_colour.SetBorderSize(0)
leg_colour.AddEntry(g_Vov3_Na22,'Na22 Vov3','PL')
leg_colour.AddEntry(g_Vov5_Na22,'Na22 Vov5','PL')
leg_colour.AddEntry(g_Vov7_Na22,'Na22 Vov7','PL')
leg_colour.AddEntry(g_Vov3_Co60,'Co60 Vov3','PL')
leg_colour.AddEntry(g_Vov5_Co60,'Co60 Vov5','PL')
leg_colour.AddEntry(g_Vov7_Co60,'Co60 Vov7','PL')
leg_colour.Draw('same')

for i in range(g_Vov3_Na22.GetN()):
	print('Na22Vov3  '+str(g_Vov3_Na22.GetPointX(i)))

for i in range(g_Vov5_Na22.GetN()):
	print('Na22Vov5  '+str(g_Vov5_Na22.GetPointX(i)))

for i in range(g_Vov7_Na22.GetN()):
	print('Na22Vov7  '+str(g_Vov7_Na22.GetPointX(i)))

for i in range(g_Vov3_Co60.GetN()):
	print('Co60Vov3  '+str(g_Vov3_Co60.GetPointX(i)))

for i in range(g_Vov5_Co60.GetN()):
	print('Co60Vov5  '+str(g_Vov5_Co60.GetPointX(i)))

for i in range(g_Vov7_Co60.GetN()):
	print('Co60Vov7  '+str(g_Vov7_Co60.GetPointX(i)))







funcLinearizedEnergy.SetParameter(0, funcLinearizedEnergy.GetParameter(0))
funcLinearizedEnergy.SetParameter(1, funcLinearizedEnergy.GetParameter(1))
funcLinearizedEnergy.SetParameter(2, funcLinearizedEnergy.GetParameter(2))
funcLinearizedEnergy.SetParameter(3, funcLinearizedEnergy.GetParameter(3))
funcLinearizedEnergy.SetParameter(4, funcLinearizedEnergy.GetParameter(4))




### LASER th20
i = 0;
for tune in laserTunes:
	#/home/cmsdaq/guglielmi/Lab5015Analysis/plots
	
	#fname  = '%sanalyzeTOFPET2_wirelessBarInArray_HDR2_UVlaser_tune%.2d_plots.root'%(otherdir,tune)
	fname  = '%sModuleCharTOFPET2_wirelessBarInArray_HDR2_UVlaser_tune%.2d_step3.root'%(otherdir,tune)
	fname2  = '%sModuleCharTOFPET2_wirelessBarInArray_HDR2_UVlaser_tune%.2d_step2_run%d.root'%(otherdir,tune,run[i])
	
	f[tune] = ROOT.TFile.Open(fname)
	fstep2[tune] = ROOT.TFile.Open(fname2)
	
	c1 = ROOT.TCanvas('c1%d' %tune, 'c1%d' %tune)
	h1_energy[tune]  = fstep2[tune].Get('h1_energy_bar00L-R_Vov%.01f_th20'%Vov)

	g_tRes_vs_threshold[tune] = f[tune].Get('g_timeRes_vs_th_Vov%.1f_bar00_enBin1'%Vov) 
	energy[tune] = f[tune].Get('g_en511_vs_bar_Vov%.01f_th20' %Vov).GetPointY(0)
	energy_err = energy_err_sys * energy[tune]
  
	
	#ALLA SOGLIA 20!!!!
	g_energy_vs_laserTune.SetPoint(g_energy_vs_laserTune.GetN(), tune, energy[tune])
	g_energy_vs_laserTune.SetPointError(g_energy_vs_laserTune.GetN()-1, 0, energy_err)




	g = f[tune].Get('g_timeRes_vs_enBin_bar00_Vov%.01f_th20' %Vov)
	tRes = g.GetPointY(0)
	tRes_err = g.GetErrorY(0)
	tRes_err = pow((tRes_err*tRes_err + tRes_err_sys*tRes_err_sys), 0.5)
	g_tRes_vs_energy.SetPoint(g_tRes_vs_energy.GetN(), energy[tune], tRes)
	g_tRes_vs_energy.SetPointError(g_tRes_vs_energy.GetN()-1, energy_err, tRes_err)


	linearizedEnergy[tune] = funcLinearizedEnergy.Eval(energy[tune])
	k = GetSaturationCorrection(Ncells,linearizedEnergy[tune], LY, PDE, LCE)
	linearizedEnergy[tune] = linearizedEnergy[tune]/k
	g_tRes_vs_linearizedEnergy.SetPoint(g_tRes_vs_linearizedEnergy.GetN(), linearizedEnergy[tune], tRes)
	linEnergy_err = 0.5 *(funcLinearizedEnergy.Eval(energy[tune]*(1+energy_err_sys)) - funcLinearizedEnergy.Eval(energy[tune]*(1-energy_err_sys)) )/k
	g_tRes_vs_linearizedEnergy.SetPointError(g_tRes_vs_linearizedEnergy.GetN()-1, linEnergy_err, tRes_err)

	outFile.cd()
	g_tRes_vs_linearizedEnergy.Write('g_tRes_vs_linearizedEnergy')

	g_bestTh = f[tune].Get('g_timeRes_vs_enBin_bar00_Vov%.01f_bestTh' %Vov)
	tRes_bestTh = g_bestTh.GetPointY(0)
	tRes_err_bestTh = g_bestTh.GetErrorY(0)
	tRes_err_bestTh = pow((tRes_err_bestTh*tRes_err_bestTh + tRes_err_sys*tRes_err_sys), 0.5)
	g_tRes_vs_energy_bestTh.SetPoint(g_tRes_vs_energy_bestTh.GetN(), energy[tune], tRes_bestTh)
	g_tRes_vs_energy_bestTh.SetPointError(g_tRes_vs_energy_bestTh.GetN()-1, energy_err, tRes_err_bestTh)

	g_tRes_vs_linearizedEnergy_bestTh.SetPoint(g_tRes_vs_linearizedEnergy_bestTh.GetN(), linearizedEnergy[tune], tRes_bestTh)
	g_tRes_vs_linearizedEnergy_bestTh.SetPointError(g_tRes_vs_linearizedEnergy_bestTh.GetN()-1, linEnergy_err, tRes_err_bestTh)
	i = i +1


### Na22
g_tRes_vs_energy_Na22 = ROOT.TGraphErrors()
g_tRes_vs_energy_Na22_bestTh = ROOT.TGraphErrors()
g_tRes_vs_energy_th20_Na22 = ROOT.TGraphErrors()
g_tRes_vs_linearizedEnergy_Na22 = ROOT.TGraphErrors()
g_tRes_vs_linearizedEnergy_Na22_bestTh = ROOT.TGraphErrors()
g_diff_linearizedEnergy_vs_peak = ROOT.TGraphErrors()

#f1 = ROOT.TFile.Open('%smoduleCharacterization_single_HPK_HDR2_step3.root'%plotdir)



g_tRes_vs_energy_Na22=f1.Get('g_timeRes_vs_enBin_bar00_Vov%.01f_thRef20' %Vov)
g_tRes_vs_energy_Na22_bestTh=f1.Get('g_timeRes_vs_enBin_bar00_Vov%.01f_bestTh' %Vov)
g_tRes_vs_energy_th20_Na22=f1.Get('g_timeRes_vs_enBin_bar00_Vov%.01f_th20' %Vov)

for eBin in range(2,13):
	g_tRes_vs_threshold[eBin] = f1.Get('g_timeRes_vs_th_Vov%.1f_bar00_enBin%d' %(Vov,eBin))
	

	c1 = ROOT.TCanvas('c1%d' %eBin, 'c1%d' %eBin)
	linearizedEnergy[eBin] = funcLinearizedEnergy.Eval(g_tRes_vs_energy_th20_Na22.GetPointX(eBin-1))
	#print('Na22')
	k = GetSaturationCorrection(Ncells,linearizedEnergy[eBin], LY, PDE, LCE)
	linearizedEnergy[eBin] = linearizedEnergy[eBin]/k
	#print(linearizedEnergy[eBin])
	g_tRes_vs_linearizedEnergy_Na22.SetPoint(g_tRes_vs_linearizedEnergy_Na22.GetN(), linearizedEnergy[eBin], g_tRes_vs_energy_Na22.GetPointY(eBin-1))
	linEnergy_err = 0.5 *(funcLinearizedEnergy.Eval(g_tRes_vs_energy_th20_Na22.GetPointX(eBin-1)*(1+energy_err_sys)) - funcLinearizedEnergy.Eval(g_tRes_vs_energy_th20_Na22.GetPointX(eBin-1)*(1-energy_err_sys)) )/k
	tRes_err = g_tRes_vs_energy_Na22.GetErrorY(eBin-1)
	tRes_err = pow((tRes_err*tRes_err + tRes_err_sys*tRes_err_sys), 0.5)
	g_tRes_vs_linearizedEnergy_Na22.SetPointError(g_tRes_vs_linearizedEnergy_Na22.GetN()-1, linEnergy_err, tRes_err)

	outFile.cd()
	g_tRes_vs_linearizedEnergy_Na22.Write('g_tRes_vs_linearizedEnergy_Na22')

	tRes_err = g_tRes_vs_energy_Na22_bestTh.GetErrorY(eBin-1)
	tRes_err = pow((tRes_err*tRes_err + tRes_err_sys*tRes_err_sys), 0.5)
	g_tRes_vs_linearizedEnergy_Na22_bestTh.SetPoint(g_tRes_vs_linearizedEnergy_Na22_bestTh.GetN(), linearizedEnergy[eBin], g_tRes_vs_energy_Na22_bestTh.GetPointY(eBin-1))
	g_tRes_vs_linearizedEnergy_Na22_bestTh.SetPointError(g_tRes_vs_linearizedEnergy_Na22_bestTh.GetN()-1, linEnergy_err, tRes_err)
	

### Na22 coinc
g_tRes_vs_energy_Na22_coinc = ROOT.TGraphErrors()
g_tRes_vs_energy_Na22_coinc_bestTh = ROOT.TGraphErrors()
g_tRes_vs_energy_th20_Na22_coinc = ROOT.TGraphErrors()
g_tRes_vs_linearizedEnergy_Na22_coinc = ROOT.TGraphErrors()
g_tRes_vs_linearizedEnergy_Na22_coinc_bestTh = ROOT.TGraphErrors()
g_diff_linearizedEnergy_vs_peak_coinc = ROOT.TGraphErrors()

f7 = ROOT.TFile.Open('%smoduleCharacterization_HPK_HDR2_externalBarCoincidence_step3.root'%plotdir)
g_tRes_vs_energy_Na22_coinc_bestTh=f7.Get('g_timeRes_vs_enBin_bar00_Vov%.01f_bestTh' %Vov)
g_tRes_vs_energy_Na22_coinc=f7.Get('g_timeRes_vs_enBin_bar00_Vov%.01f_thRef20' %Vov)
g_tRes_vs_energy_th20_Na22_coinc=f7.Get('g_timeRes_vs_enBin_bar00_Vov%.01f_th20' %Vov)
g_tRes_vs_threshold_Na22_coinc = f7.Get('g_timeRes_vs_th_Vov%.1f_bar00_enBin%d' %(Vov,eBin))

for eBin in range(1,2):
	linearizedEnergy[eBin] = funcLinearizedEnergy.Eval(g_tRes_vs_energy_th20_Na22_coinc.GetPointX(eBin-1))
	#print('Coinc:')
	k = GetSaturationCorrection(Ncells,linearizedEnergy[eBin], LY, PDE, LCE)
	linearizedEnergy[eBin] = linearizedEnergy[eBin]/k
	#print(linearizedEnergy[eBin])
	g_tRes_vs_linearizedEnergy_Na22_coinc.SetPoint(g_tRes_vs_linearizedEnergy_Na22_coinc.GetN(), linearizedEnergy[eBin], g_tRes_vs_energy_Na22_coinc.GetPointY(eBin-1))
	linEnergy_err = 0.5 *(funcLinearizedEnergy.Eval(g_tRes_vs_energy_th20_Na22_coinc.GetPointX(eBin-1)*(1+energy_err_sys)) - funcLinearizedEnergy.Eval(g_tRes_vs_energy_th20_Na22_coinc.GetPointX(eBin-1)*(1-energy_err_sys)) )/k
	tRes_err = g_tRes_vs_energy_Na22_coinc.GetErrorY(eBin-1)
	tRes_err = pow((tRes_err*tRes_err + tRes_err_sys*tRes_err_sys), 0.5)
	g_tRes_vs_linearizedEnergy_Na22_coinc.SetPointError(g_tRes_vs_linearizedEnergy_Na22_coinc.GetN()-1, linEnergy_err, tRes_err)
	
	tRes_err = g_tRes_vs_energy_Na22_coinc_bestTh.GetErrorY(eBin-1)
	tRes_err = pow((tRes_err*tRes_err + tRes_err_sys*tRes_err_sys), 0.5)
	g_tRes_vs_linearizedEnergy_Na22_coinc_bestTh.SetPoint(g_tRes_vs_linearizedEnergy_Na22_coinc_bestTh.GetN(), linearizedEnergy[eBin], g_tRes_vs_energy_Na22_coinc_bestTh.GetPointY(eBin-1))
	g_tRes_vs_linearizedEnergy_Na22_coinc_bestTh.SetPointError(g_tRes_vs_linearizedEnergy_Na22_coinc_bestTh.GetN()-1, linEnergy_err, tRes_err)	


	
### Co60
g_tRes_vs_linearizedEnergy_Co60 = ROOT.TGraphErrors()
g_tRes_vs_energy_Co60 = ROOT.TGraphErrors()
g_tRes_vs_energy_th20_Co60 = ROOT.TGraphErrors()
g_tRes_vs_threshold_Co60 = ROOT.TGraphErrors()

f_Co60 = ROOT.TFile.Open('%smoduleCharacterization_single_HPK_HDR2_Co60_Res_step3.root'%plotdir)
#g_tRes_vs_linearizedEnergy_Co60=f_Co60.Get('g_timeRes_vs_linearizedEnergy_Vov%.01f_bestTh' %Vov)
g_tRes_vs_energy_Co60 = f_Co60.Get('g_timeRes_vs_enBin_bar00_Vov%.01f_th20' %Vov)
g_tRes_vs_energy_th20_Co60 = f_Co60.Get('g_timeRes_vs_enBin_bar00_Vov%.01f_th20' %Vov)
g_tRes_vs_threshold_Co60 = f_Co60.Get('g_timeRes_vs_th_Vov%.1f_bar00_enBin1' %Vov)


for eBin in range(1,2):
	linearizedEnergy[eBin] = funcLinearizedEnergy.Eval(g_tRes_vs_energy_th20_Co60.GetPointX(eBin-1)/0.93)
	#k = GetSaturationCorrection(Ncells, linearizedEnergy[eBin], LY, PDE, LCE, ECF)
	#print ('Co60')
	
	k = GetSaturationCorrection(Ncells,linearizedEnergy[eBin], LY, PDE, LCE)
	linearizedEnergy[eBin] = linearizedEnergy[eBin]/k
	#print(linearizedEnergy[eBin])
	g_tRes_vs_linearizedEnergy_Co60.SetPoint(g_tRes_vs_linearizedEnergy_Co60.GetN(), linearizedEnergy[eBin], g_tRes_vs_energy_Co60.GetPointY(eBin-1))
	linEnergy_err = 0.5 *(funcLinearizedEnergy.Eval(g_tRes_vs_energy_th20_Co60.GetPointX(eBin-1)*(1+energy_err_sys)) - funcLinearizedEnergy.Eval(g_tRes_vs_energy_th20_Co60.GetPointX(eBin-1)*(1-energy_err_sys)) )/k
	tRes_err = g_tRes_vs_energy_Co60.GetErrorY(eBin-1)
	tRes_err = pow((tRes_err*tRes_err + tRes_err_sys*tRes_err_sys), 0.5)
	g_tRes_vs_linearizedEnergy_Co60.SetPointError(g_tRes_vs_linearizedEnergy_Co60.GetN()-1, linEnergy_err, tRes_err)

	outFile.cd()
	g_tRes_vs_linearizedEnergy_Co60.Write('g_tRes_vs_linearizedEnergy_Co60')
  
#find laser tune that gives energy closest to 511,1275 Na22
laserTune = {}
refPeaks = {5, 10}
refPeaks_Co60 = {1}
for eBin in refPeaks:
#for eBin in range(1, 13):
  d  = 999999
  for i,tune in enumerate(laserTunes):
    if (abs(g_tRes_vs_energy_Na22.GetPointX(eBin-1)-energy[tune])  < d):
      laserTune[eBin]= tune
      d = abs(g_tRes_vs_energy_Na22.GetPointX(eBin-1)-energy[tune]) 
  print 'laser tune for %d a.u. Na22 --> %d'%(g_tRes_vs_energy_Na22.GetPointX(eBin-1), laserTune[eBin])

#find laser tune that gives energy closest to sumPeak (2.5MeV) Co60
for eBin in refPeaks_Co60:
#for eBin in range(1, 13):
  d  = 999999
  for i,tune in enumerate(laserTunes):
    if (abs(g_tRes_vs_energy_Co60.GetPointX(eBin-1)-energy[tune])  < d):
      laserTune[eBin]= tune
      d = abs(g_tRes_vs_energy_Co60.GetPointX(eBin-1)-energy[tune]) 
  print 'laser tune for %d a.u. Co60 --> %d'%(g_tRes_vs_energy_Co60.GetPointX(eBin-1), laserTune[eBin])  

# check diff in quadrature of Na22 - laser vs threshold
#g_diff_vs_threshold = {}
#for eBin in refPeaks:
#  g_diff_vs_threshold[eBin] = ROOT.TGraphErrors()
#  for i in range(0, g_tRes_vs_threshold[eBin].GetN()):
#    resNa22 = ROOT.Double(0)
#    th = ROOT.Double(0)
#    g_tRes_vs_threshold[eBin].GetPoint(i, th, resNa22)
#    resLaser = g_tRes_vs_threshold[laserTune[eBin]].Eval(th)
#    diff = -1
#    if (resNa22 > resLaser):
#      diff = math.sqrt(resNa22*resNa22-resLaser*resLaser)
#      g_diff_vs_threshold[eBin].SetPoint(g_diff_vs_threshold[eBin].GetN(), th, diff)
  
### PLOTS
  
tl = ROOT.TLatex( 0.16, 0.93, 'LYSO:Ce 3x3x57 mm^{2} - HDR2 - TOFPET2' )
tl.SetNDC()
tl.SetTextSize(0.030)

 
# plot laser energy
canvas1 = ROOT.TCanvas('c_energy_laser','c_energy_laser')
canvas1.SetGridy()
canvas1.SetGridx()
hdummy1 = ROOT.TH2F('hdummy1','', 100, 0, 30, 100, 0.1, 4000)
hdummy1.GetXaxis().SetTitle('TOFPET2 amp (a.u.)')
hdummy1.GetYaxis().SetTitle('counts')
hdummy1.GetYaxis().SetTitleOffset(1.5)
hdummy1.Draw()
for i,tune in enumerate(laserTunes):
	h1_energy[tune].SetLineWidth(2)
	h1_energy[tune].SetLineColor(51 + i*3)
	h1_energy[tune].Draw('histo same')
tl.Draw()
 


# plot laser tune vs energy at th20
canvas4 = ROOT.TCanvas('c_energy_vs_laserTune','c_energy_vs_laserTune')
canvas4.SetGridy()
canvas4.SetGridx()
hdummy4 = ROOT.TH2F('hdummy4','', 100, 0, 100, 100, 0, 30)
hdummy4.GetYaxis().SetTitle('energy')
hdummy4.GetXaxis().SetTitle('laser tune')
hdummy4.Draw()
g_energy_vs_laserTune.Draw('plsame')
tl.Draw()

# plot time resolution vs threshold for laser and Na22 and Co60
canvas3 = ROOT.TCanvas('c_tRes_vs_threshold','c_tRes_vs_threshold')
canvas3.Divide(2)
#canvas3.SetGridy()
#canvas3.SetGridx()
hdummy3 = ROOT.TH2F('hdummy3','', 100, 0, 120, 250, 0, 700)
hdummy3.GetXaxis().SetTitle('threshold (DAC)')
hdummy3.GetYaxis().SetTitle('#sigma_{bar} (ps)')

legend1 = ROOT.TLegend(0.75,0.2,0.95,0.9)
legend1.SetBorderSize(0)
legend2 = ROOT.TLegend(0.75,0.2,0.95,0.9)
legend2.SetBorderSize(0)
canvas3.cd(1)
hPad1 = ROOT.gPad.DrawFrame(0.,0.,30.,7.);
hPad1.Draw();
gPad.SetGridx();
gPad.SetGridy();
hdummy3.Draw()
for i,tune in enumerate(laserTunes):
	g_tRes_vs_threshold[tune].SetMarkerColor(i+1)
	g_tRes_vs_threshold[tune].SetLineColor(i+1)
	if i == 9:
		g_tRes_vs_threshold[tune].SetMarkerColor(len(laserTunes)+1)
		g_tRes_vs_threshold[tune].SetLineColor(len(laserTunes)+1)
	g_tRes_vs_threshold[tune].SetMarkerStyle(26)
	g_tRes_vs_threshold[tune].Draw('plsame')
	legend1.AddEntry(g_tRes_vs_threshold[tune], 'laser tune %d'%tune, 'pl')
legend1.Draw()
tl.Draw()
canvas3.cd(2)
hPad2 = ROOT.gPad.DrawFrame(0.,0.,2.,500.);
hPad2.Draw();
gPad.SetGridx();
gPad.SetGridy();
hdummy3.Draw()
#g_tRes_vs_threshold[5].SetMarkerColor(eBin)
#g_tRes_vs_threshold[5].SetLineColor(eBin)
#g_tRes_vs_threshold[5].SetMarkerStyle(24)
#g_tRes_vs_threshold[5].Draw('plsame')
for eBin in range(2,13):
	g_tRes_vs_threshold[eBin].SetMarkerColor(eBin)
	g_tRes_vs_threshold[eBin].SetLineColor(eBin)
	g_tRes_vs_threshold[eBin].SetMarkerStyle(24)
	g_tRes_vs_threshold[eBin].Draw('plsame')
	legend2.AddEntry(g_tRes_vs_threshold[eBin], 'Na22 enBin %d'%eBin, 'pl')
g_tRes_vs_threshold_Co60.SetMarkerColor(kOrange)
g_tRes_vs_threshold_Co60.SetLineColor(kOrange)
g_tRes_vs_threshold_Co60.SetMarkerStyle(24)
g_tRes_vs_threshold_Co60.Draw('plsame')
legend2.AddEntry(g_tRes_vs_threshold_Co60, 'Co60 enBin SumPeak', 'pl')

'''g_tRes_vs_threshold_Na22_coinc.SetMarkerColor(kViolet)
g_tRes_vs_threshold_Na22_coinc.SetLineColor(kViolet)
g_tRes_vs_threshold_Na22_coinc.SetMarkerStyle(24)
g_tRes_vs_threshold_Na22_coinc.Draw('plsame')
legend2.AddEntry(g_tRes_vs_threshold_Na22_coinc, 'Na22 enBin 5 Coinc', 'pl')'''

legend2.Draw()
tl.Draw()

# plot time resolution vs threshold for laser and Na22
#canvas7 = ROOT.TCanvas('c_diff_vs_threshold','c_diff_vs_threshold')
#canvas7.SetGridy()
#canvas7.SetGridx()
#hdummy7 = ROOT.TH2F('hdummy7','', 100, 0, 64, 250, 0, 250)
#hdummy7.GetXaxis().SetTitle('threshold (DAC)')
#hdummy7.GetYaxis().SetTitle('diff (ps)')
#hdummy7.Draw()
#for i,peak in enumerate(refPeaks):
#  g_diff_vs_threshold[eBin].SetMarkerColor(i+1)
#  g_diff_vs_threshold[eBin].SetLineColor(i+1)
#  g_diff_vs_threshold[eBin].SetMarkerStyle(21)
#  g_diff_vs_threshold[eBin].Draw('plsame')
#tl.Draw()


# plot time resolution vs energy th20
canvas2 = ROOT.TCanvas('c_tRes_vs_energy','c_tRes_vs_energy')
canvas2.SetGridy()
canvas2.SetGridx()
hdummy2 = ROOT.TH2F('hdummy2','', 100, 0, 30, 250, 0, 500)
hdummy2.GetXaxis().SetTitle('TOFPET2 amp (a.u.)')
hdummy2.GetYaxis().SetTitle('#sigma_{bar} (ps)')
hdummy2.Draw()
g_tRes_vs_energy.SetLineWidth(2)
g_tRes_vs_energy.Draw('plsame')
g_tRes_vs_energy_Na22.SetMarkerStyle(21)
g_tRes_vs_energy_Na22.SetMarkerSize(0.7)
g_tRes_vs_energy_Na22.SetMarkerColor(ROOT.kRed)
g_tRes_vs_energy_Na22.SetLineColor(ROOT.kRed)
g_tRes_vs_energy_Na22.Draw('psame')
g_tRes_vs_energy_Na22_coinc.SetMarkerStyle(21)
g_tRes_vs_energy_Na22_coinc.SetMarkerSize(0.7)
g_tRes_vs_energy_Na22_coinc.SetMarkerColor(ROOT.kBlue)
g_tRes_vs_energy_Na22_coinc.SetLineColor(ROOT.kBlue)
g_tRes_vs_energy_Na22_coinc.Draw('psame')
g_tRes_vs_energy_Co60.SetMarkerStyle(21)
g_tRes_vs_energy_Co60.SetMarkerSize(0.7)
g_tRes_vs_energy_Co60.SetMarkerColor(ROOT.kGreen)
g_tRes_vs_energy_Co60.SetLineColor(ROOT.kGreen)
g_tRes_vs_energy_Co60.Draw('psame')

#leg = ROOT.TLegend(0.15,0.13,0.33,0.28)
leg = ROOT.TLegend(0.33,0.73,0.53,0.85)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.AddEntry(g_tRes_vs_energy,'UV laser','PL')
leg.AddEntry(g_tRes_vs_energy_Na22,'^{22}Na','PL')
leg.AddEntry(g_tRes_vs_energy_Na22_coinc,'^{22}Na coincidence','PL')
leg.AddEntry(g_tRes_vs_energy_Co60,'^{60}Co','PL')
leg.Draw('same')
tl.Draw()

#plot function for linearized energy
canvas6 = ROOT.TCanvas('c_function_linearity','c_function_linearity')
canvas6.SetGridy()
canvas6.SetGridx()
hdummy6 = ROOT.TH2F('hdummy6','', 30, 0, 30, 6, 0, 6)
hdummy6.GetXaxis().SetTitle('TOFPET2 amp (a.u.)')
hdummy6.GetYaxis().SetTitle('energy (MeV)')
hdummy6.Draw()
funcLinearizedEnergy.SetLineWidth(2)
funcLinearizedEnergy.Draw('same')

# plot time resolution vs linearized energy at th20
canvas5 = ROOT.TCanvas('c_tRes_vs_linearizedEnergy','c_tRes_vs_linearizedEnergy')
canvas5.SetGridy()
canvas5.SetGridx()
hdummy5 = ROOT.TH2F('hdummy5','', 100, 0, 7.0, 250, 0, 500)
hdummy5.GetXaxis().SetTitle('energy (MeV)')
hdummy5.GetYaxis().SetTitle('#sigma_{bar} (ps)')
hdummy5.Draw()
g_tRes_vs_linearizedEnergy.SetLineWidth(2)
g_tRes_vs_linearizedEnergy.SetMarkerSize(0.7)
g_tRes_vs_linearizedEnergy.Draw('psame')
'''fitFunc = ROOT.TF1('fitFunc','sqrt([0]*[0] + [1]*[1]/pow(x,[2])/pow(x,[2]))', 0.0,10.0)
fitFunc.SetRange(0.2,6.7)
fitFunc.SetLineWidth(2)
fitFunc.SetParName(0,'c')
fitFunc.SetParName(1,'s')
fitFunc.SetParName(2,'#alpha')
fitFunc.SetLineColor(1)
fitFunc.SetParameter(0, 15)
fitFunc.SetParameter(1, 50)
fitFunc.SetParameter(2, 0.5)
g_tRes_vs_linearizedEnergy.Fit('fitFunc','QRS')'''
fitFunc = ROOT.TF1('fitFunc','sqrt(pow([0]/x,2)+pow([1]/sqrt(x),2)+pow([2],2))',0.0,10.0)
fitFunc.SetRange(0.1, 7.0)
fitFunc.SetParName(0,'n')
fitFunc.SetParName(1,'s')
fitFunc.SetParName(2,'c')
fitFunc.SetLineColor(kBlack)
fitFunc.SetParameter(0, 50)
fitFunc.SetParameter(1, 100)
fitFunc.FixParameter(2, 15)
fitFunc.SetParLimits(0, 0.00001, 1000)
fitFunc.SetParLimits(1, 0.00001, 1000)

fitFunc2 = ROOT.TF1('fitFunc2','sqrt(pow([0]/x,2)+pow([1]/sqrt(x),2)+pow([2],2))',0.0,10.0)
fitFunc2.SetRange(0.1, 7.0)
fitFunc2.SetParName(0,'n')
fitFunc2.SetParName(1,'s')
fitFunc2.SetParName(2,'c')
fitFunc2.SetLineColor(kBlack)
fitFunc2.SetParameter(0, 50)
fitFunc2.SetParameter(1, 100)
fitFunc2.FixParameter(2, 15)
fitFunc2.SetParLimits(0, 0.00001, 1000)
fitFunc2.SetParLimits(1, 0.00001, 1000)
#fitFunc.SetParLimits(2, 0.00001, 1000)
g_tRes_vs_linearizedEnergy.Fit('fitFunc','QRS')
fitFunc2.SetLineColor(kRed)
g_tRes_vs_linearizedEnergy_Na22.Fit('fitFunc2','QRS+')
canvas5.Modified()
canvas5.Update()
gPad.Update()
st = g_tRes_vs_linearizedEnergy.GetListOfFunctions().FindObject('stats')
st.SetX1NDC(0.62)
st.SetX2NDC(0.98)
st.SetY1NDC(0.59)
st.SetY2NDC(0.79)
g_tRes_vs_linearizedEnergy_Na22.SetMarkerStyle(21)
g_tRes_vs_linearizedEnergy_Na22.SetMarkerSize(0.7)
g_tRes_vs_linearizedEnergy_Na22.SetMarkerColor(ROOT.kRed)
g_tRes_vs_linearizedEnergy_Na22.SetLineColor(ROOT.kRed)
g_tRes_vs_linearizedEnergy_Na22.SetLineWidth(2)
g_tRes_vs_linearizedEnergy_Na22.Draw('psame')
g_tRes_vs_linearizedEnergy_Na22_coinc.SetMarkerStyle(21)
g_tRes_vs_linearizedEnergy_Na22_coinc.SetMarkerSize(0.7)
g_tRes_vs_linearizedEnergy_Na22_coinc.SetMarkerColor(ROOT.kBlue)
g_tRes_vs_linearizedEnergy_Na22_coinc.SetLineColor(ROOT.kBlue)
g_tRes_vs_linearizedEnergy_Na22_coinc.SetLineWidth(2)
g_tRes_vs_linearizedEnergy_Na22_coinc.Draw('psame')
g_tRes_vs_linearizedEnergy_Co60.SetLineWidth(2)
g_tRes_vs_linearizedEnergy_Co60.SetMarkerSize(0.7)
g_tRes_vs_linearizedEnergy_Co60.SetMarkerStyle(21)
g_tRes_vs_linearizedEnergy_Co60.SetMarkerColor(kGreen)
g_tRes_vs_linearizedEnergy_Co60.Draw('psame')
canvas5.Modified()
canvas5.Update()
leg.Draw('same')
tl.Draw()

print 'Time resolution at 4.2 MeV = %.01f ps'%fitFunc.Eval(4.2)
vline = ROOT.TLine(4.2,0,4.2,fitFunc.Eval(4.2))
vline.SetLineStyle(2)
vline.Draw('same')

hline = ROOT.TLine(0,fitFunc.Eval(4.2),4.2,fitFunc.Eval(4.2))
hline.SetLineStyle(2)
hline.Draw('same')





# print difference in quadrature laser vs Na22
g_diff_vs_linearizedEnergy = ROOT.TGraphErrors()

for ipeak,eBin in enumerate(refPeaks):
  diff = -1
  diff = math.sqrt(fabs(g_tRes_vs_energy_Na22.Eval(g_tRes_vs_energy_Na22.GetPointX(eBin-1))*g_tRes_vs_energy_Na22.Eval(g_tRes_vs_energy_Na22.GetPointX(eBin-1)) - g_tRes_vs_energy.Eval(g_tRes_vs_energy_Na22.GetPointX(eBin-1))* g_tRes_vs_energy.Eval(g_tRes_vs_energy_Na22.GetPointX(eBin-1)) ))
  print '%d a.u.  -->   laser = %.02f ps    Na22 = %.02f ps     diff = %.02f ps '%(g_tRes_vs_energy_Na22.GetPointX(eBin-1), g_tRes_vs_energy.Eval(g_tRes_vs_energy_Na22.GetPointX(eBin-1)), g_tRes_vs_energy_Na22.Eval(g_tRes_vs_energy_Na22.GetPointX(eBin-1)), diff)
  g_diff_vs_linearizedEnergy.SetPoint(g_diff_vs_linearizedEnergy.GetN(), g_tRes_vs_linearizedEnergy_Na22.GetPointX(eBin-1), diff)
  diff_err = g_tRes_vs_energy_Na22.GetEY()[ipeak]
  g_diff_vs_linearizedEnergy.SetPointError(g_diff_vs_linearizedEnergy.GetN()-1, 0, diff_err)

# print difference in quadrature Na22 vs Na22coincidence
g_diff_vs_linearizedEnergy_Na22Coinc = ROOT.TGraphErrors()

refPeaksCoin = [5]

for ipeak,eBin in enumerate(refPeaksCoin):
  diff = -1
  diff = math.sqrt(fabs(g_tRes_vs_energy_Na22.Eval(g_tRes_vs_energy_Na22.GetPointX(eBin-1))*g_tRes_vs_energy_Na22.Eval(g_tRes_vs_energy_Na22.GetPointX(eBin-1)) - g_tRes_vs_energy_Na22_coinc.Eval(g_tRes_vs_energy_Na22_coinc.GetPointX(eBin-1))* g_tRes_vs_energy_Na22_coinc.Eval(g_tRes_vs_energy_Na22_coinc.GetPointX(eBin-1)) ))
  print '%d a.u.  -->   Na22Coinc = %.02f ps    Na22 = %.02f ps     diff = %.02f ps '%(g_tRes_vs_energy_Na22.GetPointX(eBin-1), g_tRes_vs_energy_Na22_coinc.Eval(g_tRes_vs_energy_Na22_coinc.GetPointX(eBin-1)), g_tRes_vs_energy_Na22.Eval(g_tRes_vs_energy_Na22.GetPointX(eBin-1)), diff)
  


  
# print difference in quadrature laser vs Co60
g_diff_vs_linearizedEnergy_Co60 = ROOT.TGraphErrors()
for ipeak,eBin in enumerate(refPeaks_Co60):
#for eBin in refPeaks:
  diff = -1
  #if ( g_tRes_vs_energy_Co60.Eval(g_tRes_vs_energy_Co60.GetPointX(eBin-1)) > g_tRes_vs_energy.Eval(energy[laserTune[eBin]])):
  diff = math.sqrt(fabs(g_tRes_vs_energy_Co60.Eval(g_tRes_vs_energy_Co60.GetPointX(eBin-1))*g_tRes_vs_energy_Co60.Eval(g_tRes_vs_energy_Co60.GetPointX(eBin-1)) - g_tRes_vs_energy.Eval(g_tRes_vs_energy_Co60.GetPointX(eBin-1))* g_tRes_vs_energy.Eval(g_tRes_vs_energy_Co60.GetPointX(eBin-1)) ))
  print '%d a.u.  -->   laser = %.02f ps    Co60 = %.02f ps     diff = %.02f ps '%(g_tRes_vs_energy_Co60.GetPointX(eBin-1), g_tRes_vs_energy.Eval(g_tRes_vs_energy_Co60.GetPointX(eBin-1)), g_tRes_vs_energy_Co60.Eval(g_tRes_vs_energy_Co60.GetPointX(eBin-1)), diff)
  g_diff_vs_linearizedEnergy_Co60.SetPoint(g_diff_vs_linearizedEnergy_Co60.GetN(), g_tRes_vs_linearizedEnergy_Co60.GetPointX(eBin-1), diff)
  diff_err = g_tRes_vs_energy_Co60.GetEY()[ipeak]
  g_diff_vs_linearizedEnergy_Co60.SetPointError(g_diff_vs_linearizedEnergy_Co60.GetN()-1, 0, diff_err)


# plot time resolution diff vs linearized energy
canvas8 = ROOT.TCanvas('c_diff_vs_linearizedEnergy','c_diff_vs_linearizedEnergy')
canvas8.SetGridy()
canvas8.SetGridx()
hdummy8 = ROOT.TH2F('hdummy8','', 100, 0, 6.5, 250, 0, 150)
hdummy8.GetXaxis().SetTitle('energy (MeV)')
hdummy8.GetYaxis().SetTitle('#sigma_{bar}^{source} (-) #sigma_{bar}^{laser} (ps)')
hdummy8.Draw()
g_diff_vs_linearizedEnergy.SetLineWidth(2)
g_diff_vs_linearizedEnergy.SetMarkerSize(0.8)
g_diff_vs_linearizedEnergy.SetMarkerStyle(20)
g_diff_vs_linearizedEnergy.SetLineColor(kRed)
g_diff_vs_linearizedEnergy.Draw('psame')
g_diff_vs_linearizedEnergy_Co60.SetLineWidth(2)
g_diff_vs_linearizedEnergy_Co60.SetMarkerSize(0.8)
g_diff_vs_linearizedEnergy_Co60.SetMarkerStyle(20)
g_diff_vs_linearizedEnergy_Co60.SetLineColor(kGreen)
g_diff_vs_linearizedEnergy_Co60.Draw('psame')


# plot time resolution vs energy at bestTh
canvas9 = ROOT.TCanvas('c_tRes_vs_energy_bestTh','c_tRes_vs_energy_bestTh')
canvas9.SetGridy()
canvas9.SetGridx()
hdummy9 = ROOT.TH2F('hdummy2','', 100, 0, 30, 250, 0, 500)
hdummy9.GetXaxis().SetTitle('TOFPET2 amp (a.u.)')
hdummy9.GetYaxis().SetTitle('#sigma_{bar} (ps)')
hdummy9.Draw()
g_tRes_vs_energy_bestTh.SetLineWidth(2)
g_tRes_vs_energy_bestTh.Draw('plsame')
g_tRes_vs_energy_Na22_bestTh.SetMarkerStyle(21)
g_tRes_vs_energy_Na22_bestTh.SetMarkerSize(0.7)
g_tRes_vs_energy_Na22_bestTh.SetMarkerColor(ROOT.kRed)
g_tRes_vs_energy_Na22_bestTh.SetLineColor(ROOT.kRed)
g_tRes_vs_energy_Na22_bestTh.Draw('psame')
g_tRes_vs_energy_Na22_coinc_bestTh.SetMarkerStyle(21)
g_tRes_vs_energy_Na22_coinc_bestTh.SetMarkerSize(0.7)
g_tRes_vs_energy_Na22_coinc_bestTh.SetMarkerColor(ROOT.kBlue)
g_tRes_vs_energy_Na22_coinc_bestTh.SetLineColor(ROOT.kBlue)
g_tRes_vs_energy_Na22_coinc_bestTh.Draw('psame')
g_tRes_vs_energy_Co60.SetMarkerStyle(21)
g_tRes_vs_energy_Co60.SetMarkerSize(0.7)
g_tRes_vs_energy_Co60.SetMarkerColor(ROOT.kGreen)
g_tRes_vs_energy_Co60.SetLineColor(ROOT.kGreen)
g_tRes_vs_energy_Co60.Draw('psame')

#leg = ROOT.TLegend(0.15,0.13,0.33,0.28)
leg9 = ROOT.TLegend(0.33,0.73,0.53,0.85)
leg9.SetFillStyle(0)
leg9.SetBorderSize(0)
leg9.AddEntry(g_tRes_vs_energy_bestTh,'UV laser','PL')
leg9.AddEntry(g_tRes_vs_energy_Na22_bestTh,'^{22}Na','PL')
leg9.AddEntry(g_tRes_vs_energy_Na22_coinc_bestTh,'^{22}Na coincidence','PL')
leg9.AddEntry(g_tRes_vs_energy_Co60,'^{60}Co','PL')
leg9.Draw('same')
tl.Draw()

# plot time resolution vs linearized energy at bestTh
canvas10 = ROOT.TCanvas('c_tRes_vs_linearizedEnergy_bestTh','c_tRes_vs_linearizedEnergy_bestTh')
canvas10.SetGridy()
canvas10.SetGridx()
hdummy10 = ROOT.TH2F('hdummy10','', 100, 0, 7.0, 250, 0, 500)
hdummy10.GetXaxis().SetTitle('energy (MeV)')
hdummy10.GetYaxis().SetTitle('#sigma_{bar} (ps)')
hdummy10.Draw()
g_tRes_vs_linearizedEnergy_bestTh.SetLineWidth(2)
g_tRes_vs_linearizedEnergy_bestTh.SetMarkerSize(0.7)
g_tRes_vs_linearizedEnergy_bestTh.Draw('psame')
fitFunc = ROOT.TF1('fitFunc','sqrt(pow([0]/x,2)+pow([1]/sqrt(x),2)+pow([2],2))',0.0,10.0)
fitFunc.SetRange(0.25,6.7)
fitFunc.SetParName(0,'n')
fitFunc.SetParName(1,'s')
fitFunc.SetParName(2,'c')
fitFunc.SetLineColor(kBlack)
fitFunc.SetParameter(0, 50)
fitFunc.SetParameter(1, 100)
fitFunc.SetParameter(2, 50)
g_tRes_vs_linearizedEnergy_bestTh.Fit('fitFunc','QRS')
'''fitFunc1 = ROOT.TF1('fitFunc1','sqrt(pow([0]/x,2)+pow([1]/sqrt(x),2)+pow([2],2))',0.0,10.0)
fitFunc1.SetRange(0.25,6.7)
fitFunc1.SetParName(0,'n')
fitFunc1.SetParName(1,'s')
fitFunc1.SetParName(2,'c')
fitFunc1.SetLineColor(kBlack)
fitFunc1.SetParameter(0, 50)
fitFunc1.SetParameter(1, 100)
fitFunc1.SetParameter(2, 50)
g_tRes_vs_linearizedEnergy_Na22_bestTh.Fit('fitFunc1','QRS')'''
canvas10.Modified()
canvas10.Update()
st1 = g_tRes_vs_linearizedEnergy_bestTh.GetListOfFunctions().FindObject('stats')
st1.SetX1NDC(0.6)
st1.SetX2NDC(0.90)
st1.SetY1NDC(0.7)
st1.SetY2NDC(0.89)
g_tRes_vs_linearizedEnergy_Na22_bestTh.SetMarkerStyle(21)
g_tRes_vs_linearizedEnergy_Na22_bestTh.SetMarkerSize(0.7)
g_tRes_vs_linearizedEnergy_Na22_bestTh.SetMarkerColor(ROOT.kRed)
g_tRes_vs_linearizedEnergy_Na22_bestTh.SetLineColor(ROOT.kRed)
g_tRes_vs_linearizedEnergy_Na22_bestTh.SetLineWidth(2)
g_tRes_vs_linearizedEnergy_Na22_bestTh.Draw('psame')
g_tRes_vs_linearizedEnergy_Na22_coinc_bestTh.SetMarkerStyle(21)
g_tRes_vs_linearizedEnergy_Na22_coinc_bestTh.SetMarkerSize(0.7)
g_tRes_vs_linearizedEnergy_Na22_coinc_bestTh.SetMarkerColor(ROOT.kBlue)
g_tRes_vs_linearizedEnergy_Na22_coinc_bestTh.SetLineColor(ROOT.kBlue)
g_tRes_vs_linearizedEnergy_Na22_coinc_bestTh.SetLineWidth(2)
g_tRes_vs_linearizedEnergy_Na22_coinc_bestTh.Draw('psame')
g_tRes_vs_linearizedEnergy_Co60.SetLineWidth(2)
g_tRes_vs_linearizedEnergy_Co60.SetMarkerSize(0.7)
g_tRes_vs_linearizedEnergy_Co60.SetMarkerStyle(21)
g_tRes_vs_linearizedEnergy_Co60.SetMarkerColor(kGreen)
g_tRes_vs_linearizedEnergy_Co60.Draw('psame')
canvas10.Modified()
canvas10.Update()
leg9.Draw('same')
tl.Draw()

#energy Na22, Co60 at th20
c_en_Na22 = ROOT.TCanvas('c_en_Na22','c_en_Na22')
c_en_Na22.SetGridy()
c_en_Na22.SetGridx()
hdummy_en_Na22 = ROOT.TH2F('hdummy_en_Na22','', 100, 0, 30, 100, 1, 140000)
hdummy_en_Na22.GetXaxis().SetTitle('TOFPET2 amp (a.u.)')
hdummy_en_Na22.GetYaxis().SetTitle('entries')
hdummy_en_Na22.Draw()
h_en_Vov3_th20_Na22.SetLineColor(kRed)
h_en_Vov5_th20_Na22.SetLineColor(kBlue)
h_en_Vov7_th20_Na22.SetLineColor(kGreen)
h_en_Vov3_th20_Na22.Draw('same')
h_en_Vov5_th20_Na22.Draw('same')
h_en_Vov7_th20_Na22.Draw('same')
c_en_Na22.Modified()
c_en_Na22.Update()
tl.Draw()

c_en_Co60 = ROOT.TCanvas('c_en_Co60','c_en_Co60')
c_en_Co60.SetGridy()
c_en_Co60.SetGridx()
c_en_Co60.cd()
hdummy_en_Co60 = ROOT.TH2F('hdummy_en_Co60','', 100, 0, 30, 100, 1, 200000)
hdummy_en_Co60.GetXaxis().SetTitle('TOFPET2 amp (a.u.)')
hdummy_en_Co60.GetYaxis().SetTitle('entries')
hdummy_en_Co60.Draw()
h_en_Vov3_th20_Co60.SetLineColor(kRed)
h_en_Vov5_th20_Co60.SetLineColor(kBlue)
h_en_Vov7_th20_Co60.SetLineColor(kGreen)
h_en_Vov3_th20_Co60.Draw('same')
h_en_Vov5_th20_Co60.Draw('same')
h_en_Vov7_th20_Co60.Draw('same')
c_en_Co60.Modified()
c_en_Co60.Update()
tl.Draw()

# plot source energy
canvas11 = ROOT.TCanvas('c_energy_sources','c_energy_sources')
canvas11.SetGridy()
canvas11.SetGridx()
canvas11.SetLogy()
hdummy11 = ROOT.TH2F('hdummy1','', 100, 0, 30, 100, 1, 500000)
hdummy11.GetXaxis().SetTitle('TOFPET2 amp (a.u.)')
hdummy11.GetYaxis().SetTitle('counts')
hdummy11.Draw()
h_en_Vov3_th20_Na22.SetLineColor(kRed)
h_en_Vov3_th20_Na22.SetMarkerColor(kRed)
h_en_Vov3_th20_Na22.SetLineWidth(2)
h_en_Vov3_th20_Na22.Draw('histo same')
#h1_energy[1173].Scale(300./600.)
h_en_Vov3_th20_Co60.SetLineColor(kBlack)
h_en_Vov3_th20_Co60.SetLineWidth(2)
h_en_Vov3_th20_Co60.Draw('histo same')
leg_source = ROOT.TLegend(0.73,0.73,0.85,0.85)
leg_source.SetFillStyle(0)
leg_source.SetBorderSize(0)
leg_source.AddEntry(h_en_Vov3_th20_Na22,'Na22','L')
leg_source.AddEntry(h_en_Vov3_th20_Co60,'Co60','L')
leg_source.Draw('same')
tl.Draw()




# save plots  
print 'Saving plots in %s'%outdir
os.system('mkdir %s'%outdir)
shutil.copy('/var/www/html/index.php', '%s'%outdir)


for c in [canvas1, canvas2, canvas3, canvas4, canvas5, canvas6, canvas8, canvas9, canvas10, canvas11, c_lin, c_colour, c_colour_rescaled, c_en_Na22, c_en_Co60]:
  c.SaveAs(outdir+'/'+c.GetName()+'.png')
  c.SaveAs(outdir+'/'+c.GetName()+'.pdf')
  c.SaveAs(outdir+'/'+c.GetName()+'.C')

outFile.Close()
