#! /usr/bin/env python
import os
import shutil
import glob
import math
import array
import sys
import time


import ROOT
import CMS_lumi, tdrstyle

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

def GetEnergyPeak(h, xmin):
  fGaus = ROOT.TF1('fGaus','gaus',0,50)
  fGaus.SetLineColor(1)
  binmin = h.FindBin(xmin)
  binmax = -999
  ymax = -999
  for ibin in range(binmin,h.GetNbinsX()):
    if (h.GetBinContent(ibin) > ymax): 
      ymax   = h.GetBinContent(ibin) 
      binmax = ibin
  xmax = h.GetBinCenter(int(binmax))
  fGaus.SetRange(xmax - 0.10 * xmax,  xmax + 0.10 * xmax)
  h.Fit('fGaus', 'QSR')
  fGaus.SetRange(xmax - fGaus.GetParameter(2),  xmax + fGaus.GetParameter(2) )
  h.Fit('fGaus', 'QSR')
  if (fGaus.GetParameter(1)<xmin):
    fGaus.SetRange(xmin, xmax*1.05)
    h.Fit('fGaus', 'QSR')
    fGaus.SetRange(xmax - fGaus.GetParameter(2),  xmax + fGaus.GetParameter(2) )
    h.Fit('fGaus', 'QSR')
  return fGaus.GetParameter(1), fGaus.GetParError(1)


def GetLinearizedEnergy(x, par):
  xx = x[0]
  if( xx < par[0] ): return ( par[1]*(math.exp(par[2]*xx)-1) + par[3] );
  else:              return ( par[1]*(math.exp(par[2]*par[0])-1)+par[3])/(math.exp(par[4]*par[0]))*math.exp(par[4]*xx );



def GetSaturationCorrection(Ncells, Edep, LY, PDE, LCE, ECF):
  Npe    = Edep * LY * PDE * LCE * ECF
  Nfired = Ncells * (1 - math.exp(-Npe/Ncells))
  k = Nfired/Npe
  #print Edep, LY, PDE, LCE, ECF
  #print 'Ncells = %d, Npe = %f, Nfired = %f,  k = %f'%(Ncells, Npe, Nfired, k)
  return k


###
outdir = sys.argv[1]

Vov = 3.0
laserTunes = [10,20,30,40,50,60,70,80,83,85,86,87,88,89,90]

LY = 40000. # gamma/MeV
LCE = 0.15
Ncells = 40000
PDE = 0.351
ECF = 1.023

tRes = {}
energy = {}
linearizedEnergy = {}
f = {}
h1_energy  = {}
g_tRes_vs_threshold = {}

g_tRes_vs_energy = ROOT.TGraphErrors()
g_tRes_vs_linearizedEnergy = ROOT.TGraphErrors()
g_energy_vs_laserTune = ROOT.TGraphErrors()

tRes_err_sys = 2.
energy_err_sys = 0.015 # 1.5%


# from Andrea
funcLinearizedEnergy = ROOT.TF1('funcLinearizedEnergy',GetLinearizedEnergy, 0, 100,5)
funcLinearizedEnergy.SetParameter(0,1.79617e+01)
funcLinearizedEnergy.SetParameter(1,7.07013e+00)
funcLinearizedEnergy.SetParameter(2,1.87292e-02)
funcLinearizedEnergy.SetParameter(3,-1.70555e-02)
funcLinearizedEnergy.SetParameter(4,6.59613e-02)


### LASER
for tune in laserTunes:
  fname  = '../plots/analyzeTOFPET2_wirelessBar_HDR2_UVlaser_tune%.2d_plots.root'%tune
  #print 'Reading ',fname
  f[tune] = ROOT.TFile.Open(fname)
  h1_energy[tune]  = f[tune].Get('h1_energy_barL-barR_Vov%.01f_th20'%Vov)
  g_tRes_vs_threshold[tune] = f[tune].Get('g_tRes_energyCorr_gaus_vs_th_barL-barR_Vov%.01f'%Vov)
  energy[tune] = GetEnergyPeak(h1_energy[tune], 0.1)[0]
  energy_err = energy_err_sys * energy[tune]
  
  g_energy_vs_laserTune.SetPoint(g_energy_vs_laserTune.GetN(), tune, energy[tune])
  g_energy_vs_laserTune.SetPointError(g_energy_vs_laserTune.GetN()-1, 0, energy_err)

  g = f[tune].Get('g_tRes_energyCorr_gaus_bestTh_vs_Vov_barL-barR_bestTh')
  tRes = g.Eval(Vov)
  tRes_err = f[tune].Get('h1_deltaT_energyCorr_barL-barR_Vov%.01f_th20'%Vov).GetFunction('fitFunc2_energyCorr_barL-barR_Vov%.01f_th20'%Vov).GetParError(2)
  tRes_err = math.sqrt(tRes_err*tRes_err + tRes_err_sys*tRes_err_sys)
  g_tRes_vs_energy.SetPoint(g_tRes_vs_energy.GetN(), energy[tune], tRes)
  g_tRes_vs_energy.SetPointError(g_tRes_vs_energy.GetN()-1, energy_err, tRes_err)

  linearizedEnergy[tune] = funcLinearizedEnergy.Eval(energy[tune])
  k = GetSaturationCorrection(Ncells,linearizedEnergy[tune], LY, PDE, LCE, ECF)
  linearizedEnergy[tune] = linearizedEnergy[tune]/k
  g_tRes_vs_linearizedEnergy.SetPoint(g_tRes_vs_linearizedEnergy.GetN(), linearizedEnergy[tune], tRes)
  linEnergy_err = 0.5 *(funcLinearizedEnergy.Eval(energy[tune]*(1+energy_err_sys)) - funcLinearizedEnergy.Eval(energy[tune]*(1-energy_err_sys)) )/k
  g_tRes_vs_linearizedEnergy.SetPointError(g_tRes_vs_linearizedEnergy.GetN()-1, linEnergy_err, tRes_err)




### Na22
g_tRes_vs_energy_Na22 = ROOT.TGraphErrors()
g_tRes_vs_linearizedEnergy_Na22 = ROOT.TGraphErrors()
g_diff_linearizedEnergy_vs_peak = ROOT.TGraphErrors()

refPeaks = [511, 1275, 1789, 2505]
for peak in refPeaks:
  if (peak != 2505): f[peak] = ROOT.TFile.Open('../plots/analyzeTOFPET2_wirelessBar_HDR2_Na22_%dkeV_plots.root'%peak)
  else: f[peak] = ROOT.TFile.Open('../plots/analyzeTOFPET2_wirelessBar_HDR2_Co60_%dkeV_plots.root'%peak)
  g_tRes_vs_threshold[peak] = f[peak].Get('g_tRes_energyCorr_gaus_vs_th_barL-barR_Vov%.01f'%Vov)
  h1_energy[peak] = f[peak].Get('h1_energy_barL-barR_Vov%.01f_th20'%Vov)
  xmin = 0.5*(f[peak].Get('g_energy_vs_Vov_barL_th20').Eval(Vov)+f[peak].Get('g_energy_vs_Vov_barR_th20').Eval(Vov))
  if (peak == 2505): xmin = 16.
  ctemp = ROOT.TCanvas()
  energy[peak] = GetEnergyPeak(h1_energy[peak], 0.95*xmin )[0]
  ctemp.SetLogy()
  ctemp.SaveAs(outdir+'/c_peak_%d.pdf'%peak)  
  energy_err = energy_err_sys * energy[peak]
  tRes = f[peak].Get('g_tRes_energyCorr_gaus_bestTh_vs_Vov_barL-barR_bestTh').Eval(Vov)
  tRes_err = f[peak].Get('h1_deltaT_energyCorr_barL-barR_Vov%.01f_th20'%Vov).GetFunction('fitFunc2_energyCorr_barL-barR_Vov%.01f_th20'%Vov).GetParError(2)
  tRes_err = math.sqrt(tRes_err*tRes_err + tRes_err_sys*tRes_err_sys)

  g_tRes_vs_energy_Na22.SetPoint(g_tRes_vs_energy_Na22.GetN(), energy[peak], tRes)
  g_tRes_vs_energy_Na22.SetPointError(g_tRes_vs_energy_Na22.GetN()-1, energy_err, tRes_err)

  linearizedEnergy[peak] = funcLinearizedEnergy.Eval(energy[peak])
  k = GetSaturationCorrection(Ncells, linearizedEnergy[peak], LY, PDE, LCE, ECF)
  linearizedEnergy[peak] = linearizedEnergy[peak]/k
  g_tRes_vs_linearizedEnergy_Na22.SetPoint(g_tRes_vs_linearizedEnergy_Na22.GetN(), linearizedEnergy[peak], tRes)
  linEnergy_err = 0.5 *(funcLinearizedEnergy.Eval(energy[peak]*(1+energy_err_sys)) - funcLinearizedEnergy.Eval(energy[peak]*(1-energy_err_sys)) )/k
  g_tRes_vs_linearizedEnergy_Na22.SetPointError(g_tRes_vs_linearizedEnergy_Na22.GetN()-1, linEnergy_err, tRes_err)
  #raw_input('ok?')

#cc = ROOT.TCanvas()
#g_diff_linearizedEnergy_vs_peak.Draw('ap*')
#raw_input('ok?')
  
#find laser tune that gives energy closest to 511,1275,1789 Na22
laserTune = {}
for peak in refPeaks:
  d  = 999999
  for i,tune in enumerate(laserTunes):
    if (abs(energy[peak]-energy[tune])  < d):
      laserTune[peak]= tune
      d = abs(energy[peak]-energy[tune]) 
  print 'laser tune for %d keV --> %d'%(peak, laserTune[peak])


# check diff in quadrature of Na22 - laser vs threshold
#g_diff_vs_threshold = {}
#for peak in refPeaks:
#  g_diff_vs_threshold[peak] = ROOT.TGraphErrors()
#  for i in range(0, g_tRes_vs_threshold[peak].GetN()):
#    resNa22 = ROOT.Double(0)
#    th = ROOT.Double(0)
#    g_tRes_vs_threshold[peak].GetPoint(i, th, resNa22)
#    resLaser = g_tRes_vs_threshold[laserTune[peak]].Eval(th)
#    diff = -1
#    if (resNa22 > resLaser):
#      diff = math.sqrt(resNa22*resNa22-resLaser*resLaser)
#      g_diff_vs_threshold[peak].SetPoint(g_diff_vs_threshold[peak].GetN(), th, diff)
  
### PLOTS
  
tl = ROOT.TLatex( 0.15, 0.92, 'LYSO:Ce 3x3x57 mm^{2} - HDR2 - TOFPET2' )
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

canvas11 = ROOT.TCanvas('c_energy_sources','c_energy_sources')
canvas11.SetGridy()
canvas11.SetGridx()
canvas11.SetLogy()
hdummy11 = ROOT.TH2F('hdummy1','', 100, 0, 30, 100, 1, 500000)
hdummy11.GetXaxis().SetTitle('TOFPET2 amp (a.u.)')
hdummy11.GetYaxis().SetTitle('counts')
hdummy11.Draw()
h1_energy[511].SetLineColor(1)
h1_energy[511].SetLineWidth(2)
h1_energy[511].Draw('histo same')
f[1173] = ROOT.TFile.Open('../plots/analyzeTOFPET2_wirelessBar_HDR2_Co60_1173keV_plots.root')
h1_energy[1173] = f[1173].Get('h1_energy_barL-barR_Vov%.01f_th20'%Vov)
#h1_energy[1173].Scale(300./600.)
h1_energy[1173].SetLineWidth(2)
h1_energy[1173].Draw('histo same')
tl.Draw()

# plot laser tune vs energy
canvas4 = ROOT.TCanvas('c_energy_vs_laserTune','c_energy_vs_laserTune')
canvas4.SetGridy()
canvas4.SetGridx()
hdummy4 = ROOT.TH2F('hdummy4','', 100, 0, 100, 100, 0, 30)
hdummy4.GetYaxis().SetTitle('energy')
hdummy4.GetXaxis().SetTitle('laser tune')
hdummy4.Draw()
g_energy_vs_laserTune.Draw('plsame')
tl.Draw()

# plot time resolution vs threshold for laser and Na22
canvas3 = ROOT.TCanvas('c_tRes_vs_threshold','c_tRes_vs_threshold')
canvas3.SetGridy()
canvas3.SetGridx()
hdummy3 = ROOT.TH2F('hdummy3','', 100, 0, 64, 250, 0, 350)
hdummy3.GetXaxis().SetTitle('threshold (DAC)')
hdummy3.GetYaxis().SetTitle('#sigma_{bar} (ps)')
hdummy3.Draw()
for i,peak in enumerate(refPeaks):
  g_tRes_vs_threshold[laserTune[peak]].SetMarkerColor(i+1)
  g_tRes_vs_threshold[laserTune[peak]].SetLineColor(i+1)
  g_tRes_vs_threshold[laserTune[peak]].Draw('plsame')
  g_tRes_vs_threshold[peak].SetMarkerColor(i+1)
  g_tRes_vs_threshold[peak].SetLineColor(i+1)
  g_tRes_vs_threshold[peak].SetMarkerStyle(21)
  g_tRes_vs_threshold[peak].Draw('plsame')
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
#  g_diff_vs_threshold[peak].SetMarkerColor(i+1)
#  g_diff_vs_threshold[peak].SetLineColor(i+1)
#  g_diff_vs_threshold[peak].SetMarkerStyle(21)
#  g_diff_vs_threshold[peak].Draw('plsame')
#tl.Draw()


# plot time resolution vs energy
canvas2 = ROOT.TCanvas('c_tRes_vs_energy','c_tRes_vs_energy')
canvas2.SetGridy()
canvas2.SetGridx()
hdummy2 = ROOT.TH2F('hdummy2','', 100, 0, 30, 250, 0, 250)
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
#leg = ROOT.TLegend(0.15,0.13,0.33,0.28)
leg = ROOT.TLegend(0.33,0.73,0.53,0.85)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.AddEntry(g_tRes_vs_energy,'UV laser','PL')
leg.AddEntry(g_tRes_vs_energy_Na22,'^{22}Na/^{60}Co','PL')
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

# plot time resolution vs linearized energy
canvas5 = ROOT.TCanvas('c_tRes_vs_linearizedEnergy','c_tRes_vs_linearizedEnergy')
canvas5.SetGridy()
canvas5.SetGridx()
hdummy5 = ROOT.TH2F('hdummy5','', 100, 0, 7.0, 250, 0, 250)
hdummy5.GetXaxis().SetTitle('energy (MeV)')
hdummy5.GetYaxis().SetTitle('#sigma_{bar} (ps)')
hdummy5.Draw()
g_tRes_vs_linearizedEnergy.SetLineWidth(2)
g_tRes_vs_linearizedEnergy.SetMarkerSize(0.7)
g_tRes_vs_linearizedEnergy.Draw('psame')
fitFunc = ROOT.TF1('fitFunc','sqrt([0]*[0] + [1]*[1]/pow(x,[2])/pow(x,[2]))', 0.0,10.0)
fitFunc.SetRange(0.2,6.7)
fitFunc.SetLineWidth(2)
fitFunc.SetParName(0,'c')
fitFunc.SetParName(1,'s')
fitFunc.SetParName(2,'#alpha')
fitFunc.SetLineColor(1)
fitFunc.SetParameter(0, 15)
fitFunc.SetParameter(1, 50)
fitFunc.SetParameter(2, 0.5)
g_tRes_vs_linearizedEnergy.Fit('fitFunc','QRS')
#fitFunc2 = ROOT.TF1('fitFunc2','sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/x/x)', 0.0,10)
#fitFunc2.SetRange(0.25,6.7)
#fitFunc2.SetParName(0,'c')
#fitFunc2.SetParName(1,'s')
#fitFunc2.SetParName(2,'n')
#fitFunc2.SetLineColor(3)
#fitFunc2.SetParameter(0, 15)
#fitFunc2.SetParameter(1, 50)
#fitFunc2.SetParameter(2, 50)
#g_tRes_vs_linearizedEnergy.Fit('fitFunc2','QRS')
canvas5.Modified()
canvas5.Update()
st = g_tRes_vs_linearizedEnergy.GetListOfFunctions().FindObject('stats')
st.SetX1NDC(0.6)
st.SetX2NDC(0.90)
st.SetY1NDC(0.7)
st.SetY2NDC(0.89)
g_tRes_vs_linearizedEnergy_Na22.SetMarkerStyle(21)
g_tRes_vs_linearizedEnergy_Na22.SetMarkerSize(0.7)
g_tRes_vs_linearizedEnergy_Na22.SetMarkerColor(ROOT.kRed)
g_tRes_vs_linearizedEnergy_Na22.SetLineColor(ROOT.kRed)
g_tRes_vs_linearizedEnergy_Na22.SetLineWidth(2)
g_tRes_vs_linearizedEnergy_Na22.Draw('psame')
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
for ipeak,peak in enumerate(refPeaks):
  diff = -1
  if ( g_tRes_vs_energy_Na22.Eval(energy[peak]) > g_tRes_vs_energy.Eval(energy[peak])):
    diff = math.sqrt(g_tRes_vs_energy_Na22.Eval(energy[peak])*g_tRes_vs_energy_Na22.Eval(energy[peak]) - g_tRes_vs_energy.Eval(energy[peak])* g_tRes_vs_energy.Eval(energy[peak]) )
  print '%d keV  -->   laser = %.02f ps    Na22 = %.02f ps     diff = %.02f ps '%(peak, g_tRes_vs_energy.Eval(energy[peak]), g_tRes_vs_energy_Na22.Eval(energy[peak]), diff)
  g_diff_vs_linearizedEnergy.SetPoint(g_diff_vs_linearizedEnergy.GetN(), peak/1000., diff)
  diff_err = g_tRes_vs_energy_Na22.GetEY()[ipeak]
  g_diff_vs_linearizedEnergy.SetPointError(g_diff_vs_linearizedEnergy.GetN()-1, 0, diff_err)
  

# plot time resolution diff vs linearized energy
canvas8 = ROOT.TCanvas('c_diff_vs_linearizedEnergy','c_diff_vs_linearizedEnergy')
canvas8.SetGridy()
canvas8.SetGridx()
hdummy8 = ROOT.TH2F('hdummy8','', 100, 0, 6.5, 250, 0, 100)
hdummy8.GetXaxis().SetTitle('energy (MeV)')
hdummy8.GetYaxis().SetTitle('#sigma_{bar}^{source} (-) #sigma_{bar}^{laser} (ps)')
hdummy8.Draw()
g_diff_vs_linearizedEnergy.SetLineWidth(2)
g_diff_vs_linearizedEnergy.SetMarkerSize(0.8)
g_diff_vs_linearizedEnergy.SetMarkerStyle(20)
g_diff_vs_linearizedEnergy.Draw('psame')


# save plots  
print 'Saving plots in %s'%outdir
os.system('mkdir %s'%outdir)
shutil.copy('/var/www/html/index.php', '%s'%outdir)


for c in [canvas1, canvas2, canvas3, canvas4, canvas5, canvas6, canvas8,canvas11]:
  c.SaveAs(outdir+'/'+c.GetName()+'.png')
  c.SaveAs(outdir+'/'+c.GetName()+'.pdf')
  c.SaveAs(outdir+'/'+c.GetName()+'.C')

###raw_input('ok?')
