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
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.04,'X')
ROOT.gStyle.SetLabelSize(0.04,'Y')
ROOT.gStyle.SetTitleSize(0.04,'X')
ROOT.gStyle.SetTitleSize(0.04,'Y')
ROOT.gStyle.SetTitleOffset(1.1,'X')
ROOT.gStyle.SetTitleOffset(1.2,'Y')

outdir = sys.argv[2]
print 'Saving plots in %s'%outdir
os.system('mkdir %s'%outdir)
shutil.copy('index.php', '%s'%outdir)

Vov = int(sys.argv[1])
energyMin = { 2 : 1,
              3 : 2,
              4 : 3,
              6 : 5}

arrays = ['Milano', 'Caltech']

colors = {'Caltech': ROOT.kRed, 'Milano' : ROOT.kBlue}

file = {}
for array in arrays:
  file[array] = ROOT.TFile.Open('../plotsFiles/v2/analyzeTOFPET_step2_conf9.1_%sMatrix.root'%array,"READ")


g_energyL_vs_bar = {}
g_energyR_vs_bar = {}
g_energy_vs_bar = {}
g_bestTh_vs_bar = {} 
g_tRes_tDiff_vs_bar   = {}
g_tRes_barL_vs_bar   = {} 
g_tRes_barR_vs_bar   = {} 
g_tRes_barL_noiseSubtracted_vs_bar   = {} 
g_tRes_barR_noiseSubtracted_vs_bar   = {} 
g_dVdt_vs_bar   = {} 
g_noise_vs_bar  = {} 
g_tRes_tDiff_vs_dVdt   = {}
g_tRes_tDiff_vs_noise   = {}
g_tRes_tDiff_vs_energy   = {}
g_tRes_barL_vs_noise   = {}
g_tRes_barR_vs_noise   = {}
g_tRes_barL_noiseSubtracted_vs_noise   = {}
g_tRes_barR_noiseSubtracted_vs_noise   = {}

# find bets threshold
bestTh ={} # bestTh[array][bar]
for array in arrays:
    bestTh[array] = {}
    for ibar in range(0,15):
        g = file[array].Get('g_tRes_energyCorr_gaus_vs_th_bar%dL-bar%dR_Vov%.01f'%(ibar, ibar,Vov))
        if ( g == None) :
            bestTh[array][ibar] = -1
            continue
        tRes_best = 9999.
        for i in range(0, g.GetN()):
            if (g.GetY()[i] < tRes_best):
                bestTh[array][ibar] = g.GetX()[i]
                tRes_best = g.GetY()[i]
        
# fill graphs vs bar
for array in arrays:
    g_energyL_vs_bar[array] = ROOT.TGraphErrors()
    g_energyR_vs_bar[array] = ROOT.TGraphErrors()
    g_energy_vs_bar[array] = ROOT.TGraphErrors()
    g_bestTh_vs_bar[array] = ROOT.TGraphErrors()
    g_tRes_tDiff_vs_bar[array]   = ROOT.TGraphErrors()
    g_tRes_barL_vs_bar[array]   = ROOT.TGraphErrors()
    g_tRes_barL_noiseSubtracted_vs_bar[array]   = ROOT.TGraphErrors()
    g_tRes_barR_vs_bar[array]   = ROOT.TGraphErrors()
    g_tRes_barR_noiseSubtracted_vs_bar[array]   = ROOT.TGraphErrors()
    g_dVdt_vs_bar[array]   = ROOT.TGraphErrors()
    g_noise_vs_bar[array]  = ROOT.TGraphErrors()
    g_tRes_tDiff_vs_dVdt[array]   = ROOT.TGraphErrors()
    g_tRes_tDiff_vs_noise[array]   = ROOT.TGraphErrors()
    g_tRes_tDiff_vs_energy[array]   = ROOT.TGraphErrors()
    g_tRes_barR_vs_noise[array]   = ROOT.TGraphErrors()
    g_tRes_barL_vs_noise[array]   = ROOT.TGraphErrors()
    g_tRes_barR_noiseSubtracted_vs_noise[array]   = ROOT.TGraphErrors()
    g_tRes_barL_noiseSubtracted_vs_noise[array]   = ROOT.TGraphErrors()
    for ibar in range(0,15):
          if (bestTh[array][ibar] != -1):
            g_bestTh_vs_bar[array].SetPoint(g_bestTh_vs_bar[array].GetN(), ibar, bestTh[array][ibar] )
           
            tRes_tDiff_bestTh  = file[array].Get('g_tRes_energyCorr_gaus_vs_th_bar%dL-bar%dR_Vov%.01f'%(ibar, ibar,Vov)).Eval(bestTh[array][ibar])
            tRes_barL_bestTh  = file[array].Get('g_tRes_energyCorr_gaus_vs_th_bar%dL-bar_Vov%.01f'%(ibar, Vov)).Eval(bestTh[array][ibar])
            tRes_barR_bestTh  = file[array].Get('g_tRes_energyCorr_gaus_vs_th_bar%dR-bar_Vov%.01f'%(ibar, Vov)).Eval(bestTh[array][ibar])
            dVdt_bestTh  = file[array].Get('g_dVdt_vs_th_bar%dL-bar_Vov%.01f'%(ibar, Vov)).Eval(bestTh[array][ibar])
            noiseL_bestTh = file[array].Get('g_tRes_noise_vs_th_bar%dL-bar_Vov%.01f'%(ibar, Vov)).Eval(bestTh[array][ibar])
            noiseR_bestTh = file[array].Get('g_tRes_noise_vs_th_bar%dR-bar_Vov%.01f'%(ibar, Vov)).Eval(bestTh[array][ibar])
            #noise_bestTh = file[array].Get('g_tRes_noise_vs_th_bar%dL-bar_Vov%.01f'%(ibar, Vov)).Eval(bestTh[array][ibar])
            noise_bestTh = 0.5 * (noiseL_bestTh+noiseR_bestTh)
            
            #energyL = file[array].Get('g_energy_vs_th_bar%dL_Vov%.01f'%(ibar, Vov)).Eval(bestTh[array][ibar])
            #energyR = file[array].Get('g_energy_vs_th_bar%dR_Vov%.01f'%(ibar, Vov)).Eval(bestTh[array][ibar])
            #g_energyL_vs_bar[array].SetPoint(g_energyL_vs_bar[array].GetN(), ibar, energyL)
            #g_energyR_vs_bar[array].SetPoint(g_energyR_vs_bar[array].GetN(), ibar, energyR)
            #g_energy_vs_bar[array].SetPoint(g_energy_vs_bar[array].GetN(), ibar, 0.5*(energyL+energyR))
            
            hL = file[array].Get('h1_energy_bar%dL_Vov%.01f_th%02d'%(ibar, Vov,bestTh[array][ibar]))
            hL.Rebin(2)
            hL.GetXaxis().SetRangeUser(energyMin[Vov],100)
            mpvbin = hL.GetMaximumBin()
            mpv = hL.GetBinCenter(mpvbin)
            #print mpv
            if (mpv > 19.99 and mpv < 20.01):
              hL.GetXaxis().SetRangeUser(energyMin[Vov],19.5)
              mpvbin = hL.GetMaximumBin()
              mpv = hL.GetBinCenter(mpvbin)
            fitL = ROOT.TF1('fitL','landau',0,100)
            fitL.SetLineColor(1)
            fitL.SetRange(mpv*0.85, mpv*1.5)
            hL.Fit(fitL,'QR')
            energyL =  fitL.GetParameter(1)

            #raw_input('continue?')
            
            hR = file[array].Get('h1_energy_bar%dR_Vov%.01f_th%02d'%(ibar, Vov,bestTh[array][ibar]))
            hR.Rebin(2)
            hR.GetXaxis().SetRangeUser(energyMin[Vov],100)
            mpvbin = hR.GetMaximumBin()
            mpv = hR.GetBinCenter(mpvbin)
            if (mpv > 19.99 and mpv < 20.01):
              hR.GetXaxis().SetRangeUser(energyMin[Vov],19.5)
              mpvbin = hR.GetMaximumBin()
              mpv = hR.GetBinCenter(mpvbin)
            fitR = ROOT.TF1('fitR','landau',0,100)
            fitR.SetLineColor(1)
            fitR.SetRange(mpv*0.85, mpv*1.5)
            hR.Fit(fitR,'QR')
            energyR =  fitR.GetParameter(1)

            #raw_input('continue?')

            # invert bar number for Caltech
            if ( array == 'Caltech'): ibar = 15 - ibar
            
            g_energyL_vs_bar[array].SetPoint(g_energyL_vs_bar[array].GetN(), ibar, fitL.GetParameter(1))
            g_energyL_vs_bar[array].SetPointError(g_energyL_vs_bar[array].GetN()-1, 0, fitL.GetParError(1))

            g_energyR_vs_bar[array].SetPoint(g_energyR_vs_bar[array].GetN(), ibar, fitR.GetParameter(1))
            g_energyR_vs_bar[array].SetPointError(g_energyR_vs_bar[array].GetN()-1, 0, fitR.GetParError(1))

            print array,ibar, energyL, energyR
            
            #raw_input('continue?')
            
            g_energy_vs_bar[array].SetPoint(g_energy_vs_bar[array].GetN(), ibar, 0.5*(energyL+energyR))

            
            #raw_input('continue?')
            #mybin = file[array].Get('g_energy_vs_th_bar%dL_Vov%.01f'%(ibar, Vov)).GetHistogram().GetXaxis().FindBin(bestTh[array][ibar])
            #print bestTh[array][ibar], mybin
            #errL   = file[array].Get('g_energy_vs_th_bar%dL_Vov%.01f'%(ibar, Vov)).GetErrorY[mybin]
            #mybin = file[array].Get('g_energy_vs_th_bar%dR_Vov%.01f'%(ibar, Vov)).GetXaxis().FindBin(bestTh[array][ibar])
            #errR  = file[array].Get('g_energy_vs_th_bar%dR_Vov%.01f'%(ibar, Vov)).GetEY()[mybin]
            #g_energyL_vs_bar[array].SetPointError(g_energyL_vs_bar[array].GetN()-1, 0, errL)
            #g_energyR_vs_bar[array].SetPointError(g_energyL_vs_bar[array].GetN()-1, 0, errR)
            #g_energy_vs_bar[array].SetPointError(g_energy_vs_bar[array].GetN()-1, 0, 0.5*math.sqrt(errL*errL+errR*errR))
            
            g_tRes_tDiff_vs_bar[array].SetPoint(g_tRes_tDiff_vs_bar[array].GetN(), ibar, tRes_tDiff_bestTh )
            
            g_tRes_barL_noiseSubtracted_vs_bar[array].SetPoint(g_tRes_barL_noiseSubtracted_vs_bar[array].GetN(), ibar, math.sqrt(tRes_barL_bestTh*tRes_barL_bestTh-noiseL_bestTh*noiseL_bestTh))
            g_tRes_barL_vs_bar[array].SetPoint(g_tRes_barL_vs_bar[array].GetN(), ibar, tRes_barL_bestTh )

            g_tRes_barR_noiseSubtracted_vs_bar[array].SetPoint(g_tRes_barR_noiseSubtracted_vs_bar[array].GetN(), ibar, math.sqrt(tRes_barR_bestTh*tRes_barR_bestTh-noiseR_bestTh*noiseR_bestTh))
            g_tRes_barR_vs_bar[array].SetPoint(g_tRes_barR_vs_bar[array].GetN(), ibar, tRes_barR_bestTh )

            g_dVdt_vs_bar[array].SetPoint(g_dVdt_vs_bar[array].GetN(), ibar, dVdt_bestTh )
           
            g_noise_vs_bar[array].SetPoint(g_noise_vs_bar[array].GetN(), ibar, noise_bestTh )
        
            g_tRes_tDiff_vs_dVdt[array].SetPoint(g_tRes_tDiff_vs_dVdt[array].GetN(), dVdt_bestTh, tRes_tDiff_bestTh )
            g_tRes_tDiff_vs_noise[array].SetPoint(g_tRes_tDiff_vs_noise[array].GetN(), noise_bestTh, tRes_tDiff_bestTh )
            g_tRes_tDiff_vs_energy[array].SetPoint(g_tRes_tDiff_vs_energy[array].GetN(), 0.5*(energyL+energyR), tRes_tDiff_bestTh )

            g_tRes_barL_vs_noise[array].SetPoint(g_tRes_barL_vs_noise[array].GetN(), noiseL_bestTh, tRes_barL_bestTh )
            g_tRes_barR_vs_noise[array].SetPoint(g_tRes_barR_vs_noise[array].GetN(), noiseR_bestTh, tRes_barR_bestTh )

            g_tRes_barL_noiseSubtracted_vs_noise[array].SetPoint(g_tRes_barL_noiseSubtracted_vs_bar[array].GetN(), noiseL_bestTh, math.sqrt(tRes_barL_bestTh*tRes_barL_bestTh-noiseL_bestTh*noiseL_bestTh))
            g_tRes_barR_noiseSubtracted_vs_noise[array].SetPoint(g_tRes_barR_noiseSubtracted_vs_bar[array].GetN(), noiseR_bestTh, math.sqrt(tRes_barR_bestTh*tRes_barR_bestTh-noiseR_bestTh*noiseR_bestTh))

            

#for array in arrays:
#  g_energyL_vs_bar[array].Fit('pol1')
#  g_energyL_vs_bar[array].Draw()
#  raw_input('continue?')
#  g_energyR_vs_bar[array].Fit('pol1')
#  g_energyR_vs_bar[array].Draw()
#  #raw_input('continue?')

            
leg = ROOT.TLegend(0.65, 0.75, 0.89, 0.89)
leg.SetBorderSize(0)

tl1 = {}
tl2 = {}
tl3 = {}
tl4 = {}
tl5 = {}
tl6 = {}
tl5r = {}
tl6r = {}
tl7 = {}
tl8 = {}
tl9l = {}
tl9r = {}
tl9a = {}
tl10l = {}
tl10r = {}
tl10a = {}
tl11 = {}

canvas1 = ROOT.TCanvas('c_bestTh','c_bestTh')
canvas1.SetGridy()
hdummy1 = ROOT.TH2F('hdummy1','', 16, -0.5, 15.5, 100, 0, 63)
hdummy1.GetXaxis().SetTitle('bar')
hdummy1.GetYaxis().SetTitle('threshold (DAC)')
hdummy1.Draw()
for ii, array in enumerate(arrays):
    g_bestTh_vs_bar[array].SetLineColor(colors[array])
    g_bestTh_vs_bar[array].SetMarkerColor(colors[array])
    g_bestTh_vs_bar[array].SetMarkerStyle(20)
    g_bestTh_vs_bar[array].Draw('pl same')
    tl1[array] = ROOT.TLatex( 0.17, 0.85-0.05*ii, 'mean: %.1f, RMS: %.1f'%(g_bestTh_vs_bar[array].GetMean(2), g_bestTh_vs_bar[array].GetRMS(2) ) )
    tl1[array].SetNDC()
    tl1[array].SetTextSize(0.030)
    tl1[array].SetTextColor(colors[array])
    tl1[array].Draw()
    leg.AddEntry(g_bestTh_vs_bar[array],'%s'%array, 'L')   
leg.Draw()




canvas2 = ROOT.TCanvas('c_tRes_tDiff','c_tRes_tDiff')
canvas2.SetGridy()
hdummy2 = ROOT.TH2F('hdummy2','', 16, -0.5, 15.5, 100, 40, 350)
hdummy2.GetXaxis().SetTitle('bar')
hdummy2.GetYaxis().SetTitle('#sigma(tDiff) [ps]')
hdummy2.Draw()
for ii,array in enumerate(arrays):
    g_tRes_tDiff_vs_bar[array].SetLineColor(colors[array])
    g_tRes_tDiff_vs_bar[array].SetMarkerColor(colors[array])
    g_tRes_tDiff_vs_bar[array].SetMarkerStyle(20)
    g_tRes_tDiff_vs_bar[array].Draw('p same')
    tl2[array] = ROOT.TLatex( 0.17, 0.85-0.05*ii, 'mean: %.1f ps, RMS: %.0f %%'%(g_tRes_tDiff_vs_bar[array].GetMean(2), 100*g_tRes_tDiff_vs_bar[array].GetRMS(2)/g_tRes_tDiff_vs_bar[array].GetMean(2) ) )
    tl2[array].SetNDC()
    tl2[array].SetTextSize(0.030)
    tl2[array].SetTextColor(colors[array])
    tl2[array].Draw()
leg.Draw()


canvas5 = ROOT.TCanvas('c_tRes_barL-bar','c_tRes_barL-bar')
canvas5.SetGridy()
hdummy5 = ROOT.TH2F('hdummy5','', 16, -0.5, 15.5, 100, 40, 350)
hdummy5.GetXaxis().SetTitle('bar')
hdummy5.GetYaxis().SetTitle('#sigma(L-bar) [ps]')
hdummy5.Draw()
for ii,array in enumerate(arrays):
    g_tRes_barL_vs_bar[array].SetLineColor(colors[array])
    g_tRes_barL_vs_bar[array].SetMarkerColor(colors[array])
    g_tRes_barL_vs_bar[array].SetMarkerStyle(20)
    g_tRes_barL_vs_bar[array].Draw('p same')
    tl5[array] = ROOT.TLatex( 0.17, 0.85-0.05*ii, 'mean: %.1f ps, RMS: %.0f %%'%(g_tRes_barL_vs_bar[array].GetMean(2), 100*g_tRes_barL_vs_bar[array].GetRMS(2)/g_tRes_barL_vs_bar[array].GetMean(2) ) )
    tl5[array].SetNDC()
    tl5[array].SetTextSize(0.030)
    tl5[array].SetTextColor(colors[array])
    tl5[array].Draw()
leg.Draw()


canvas5r = ROOT.TCanvas('c_tRes_barR-bar','c_tRes_barR-bar')
canvas5r.SetGridy()
hdummy5r = ROOT.TH2F('hdummy5r','', 16, -0.5, 15.5, 100, 40, 350)
hdummy5r.GetXaxis().SetTitle('bar')
hdummy5r.GetYaxis().SetTitle('#sigma(R-bar) [ps]')
hdummy5r.Draw()
for ii,array in enumerate(arrays):
    g_tRes_barR_vs_bar[array].SetLineColor(colors[array])
    g_tRes_barR_vs_bar[array].SetMarkerColor(colors[array])
    g_tRes_barR_vs_bar[array].SetMarkerStyle(20)
    g_tRes_barR_vs_bar[array].Draw('p same')
    tl5r[array] = ROOT.TLatex( 0.17, 0.85-0.05*ii, 'mean: %.1f ps, RMS: %.0f %%'%(g_tRes_barR_vs_bar[array].GetMean(2), 100*g_tRes_barR_vs_bar[array].GetRMS(2)/g_tRes_barR_vs_bar[array].GetMean(2) ) )
    tl5r[array].SetNDC()
    tl5r[array].SetTextSize(0.030)
    tl5r[array].SetTextColor(colors[array])
    tl5r[array].Draw()
leg.Draw()


canvas6 = ROOT.TCanvas('c_tRes_noiseSubtracted_barL-bar','c_tRes_noiseSubtracted_barL-bar')
canvas6.SetGridy()
hdummy6 = ROOT.TH2F('hdummy6','', 16, -0.5, 15.5, 100, 0, 350)
hdummy6.GetXaxis().SetTitle('bar')
hdummy6.GetYaxis().SetTitle('#sigma(L-bar) [ps]')
hdummy6.Draw()
for ii,array in enumerate(arrays):
    g_tRes_barL_noiseSubtracted_vs_bar[array].SetLineColor(colors[array])
    g_tRes_barL_noiseSubtracted_vs_bar[array].SetMarkerColor(colors[array])
    g_tRes_barL_noiseSubtracted_vs_bar[array].SetMarkerStyle(20)
    g_tRes_barL_noiseSubtracted_vs_bar[array].Draw('p same')
    tl6[array] = ROOT.TLatex( 0.17, 0.85-0.05*ii, 'mean: %.1f ps, RMS: %.0f %%'%(g_tRes_barL_noiseSubtracted_vs_bar[array].GetMean(2), 100*g_tRes_barL_noiseSubtracted_vs_bar[array].GetRMS(2)/g_tRes_barL_noiseSubtracted_vs_bar[array].GetMean(2) ) )
    tl6[array].SetNDC()
    tl6[array].SetTextSize(0.030)
    tl6[array].SetTextColor(colors[array])
    tl6[array].Draw()
leg.Draw()



canvas6r = ROOT.TCanvas('c_tRes_noiseSubtracted_barR-bar','c_tRes_noiseSubtracted_barR-bar')
canvas6r.SetGridy()
hdummy6r = ROOT.TH2F('hdummy6r','', 16, -0.5, 15.5, 100, 0, 350)
hdummy6r.GetXaxis().SetTitle('bar')
hdummy6r.GetYaxis().SetTitle('#sigma(R-bar) [ps]')
hdummy6r.Draw()
for ii,array in enumerate(arrays):
    g_tRes_barR_noiseSubtracted_vs_bar[array].SetLineColor(colors[array])
    g_tRes_barR_noiseSubtracted_vs_bar[array].SetMarkerColor(colors[array])
    g_tRes_barR_noiseSubtracted_vs_bar[array].SetMarkerStyle(20)
    g_tRes_barR_noiseSubtracted_vs_bar[array].Draw('p same')
    tl6r[array] = ROOT.TLatex( 0.17, 0.85-0.05*ii, 'mean: %.1f ps, RMS: %.0f %%'%(g_tRes_barR_noiseSubtracted_vs_bar[array].GetMean(2), 100*g_tRes_barR_noiseSubtracted_vs_bar[array].GetRMS(2)/g_tRes_barR_noiseSubtracted_vs_bar[array].GetMean(2) ) )
    tl6r[array].SetNDC()
    tl6r[array].SetTextSize(0.030)
    tl6r[array].SetTextColor(colors[array])
    tl6r[array].Draw()
leg.Draw()

canvas3 = ROOT.TCanvas('c_dVdt','c_dVdt')
canvas3.SetGridy()
hdummy3 = ROOT.TH2F('hdummy3','', 16, -0.5, 15.5, 100, 0, 150)
hdummy3.GetXaxis().SetTitle('bar')
hdummy3.GetYaxis().SetTitle('dV/dt (DAC/ns)')
hdummy3.Draw()
for ii, array in enumerate(arrays):
    g_dVdt_vs_bar[array].SetLineColor(colors[array])
    g_dVdt_vs_bar[array].SetMarkerColor(colors[array])
    g_dVdt_vs_bar[array].SetMarkerStyle(20)
    g_dVdt_vs_bar[array].Draw('p same')
    tl3[array] = ROOT.TLatex( 0.17, 0.85-0.05*ii, 'mean: %.1f DAC/ns, RMS: %.0f %%'%(g_dVdt_vs_bar[array].GetMean(2), 100*g_dVdt_vs_bar[array].GetRMS(2)/g_dVdt_vs_bar[array].GetMean(2) ) )
    tl3[array].SetNDC()
    tl3[array].SetTextSize(0.030)
    tl3[array].SetTextColor(colors[array])
    tl3[array].Draw()
leg.Draw()


canvas4 = ROOT.TCanvas('c_noise','c_noise')
canvas4.SetGridy()
hdummy4 = ROOT.TH2F('hdummy4','', 16, -0.5, 15.5, 100, 0, 250)
hdummy4.GetXaxis().SetTitle('bar')
hdummy4.GetYaxis().SetTitle('#sigma_{noise} (ps)')
hdummy4.Draw()
for ii, array in enumerate(arrays):
    g_noise_vs_bar[array].SetLineColor(colors[array])
    g_noise_vs_bar[array].SetMarkerColor(colors[array])
    g_noise_vs_bar[array].SetMarkerStyle(20)
    g_noise_vs_bar[array].Draw('p same')
    print 'mean: %f DAC/ns, RMS: %f DAC/ns'%(g_noise_vs_bar[array].GetMean(2), g_noise_vs_bar[array].GetRMS(2) ) 
    tl4[array] = ROOT.TLatex( 0.17, 0.85-0.05*ii, 'mean: %.1f ps, RMS: %.0f %%'%(g_noise_vs_bar[array].GetMean(2), 100*g_noise_vs_bar[array].GetRMS(2)/g_noise_vs_bar[array].GetMean(2) ) )
    tl4[array].SetNDC()
    tl4[array].SetTextSize(0.030)
    tl4[array].SetTextColor(colors[array])
    tl4[array].Draw()
leg.Draw()


canvas7 = ROOT.TCanvas('c_tRes_tDiff_vs_dVdt','c_tRes_tDiff_vs_dVdt')
canvas7.SetGridy()
hdummy7 = ROOT.TH2F('hdummy7','', 16, 0, 80, 100, 50, 250)
hdummy7.GetYaxis().SetTitle('#sigma_{tDiff} (ps)')
hdummy7.GetXaxis().SetTitle('dV/dt (DAC/ns)')
hdummy7.Draw()
for ii, array in enumerate(arrays):
    g_tRes_tDiff_vs_dVdt[array].SetLineColor(colors[array])
    g_tRes_tDiff_vs_dVdt[array].SetMarkerColor(colors[array])
    g_tRes_tDiff_vs_dVdt[array].SetMarkerStyle(20)
    g_tRes_tDiff_vs_dVdt[array].Draw('p same')
leg.Draw()

canvas8 = ROOT.TCanvas('c_tRes_tDiff_vs_noise','c_tRes_tDiff_vs_noise')
canvas8.SetGridy()
hdummy8 = ROOT.TH2F('hdummy8','', 16, 0, 120, 120, 0, 300)
hdummy8.GetYaxis().SetTitle('#sigma_{tDiff} (ps)')
hdummy8.GetXaxis().SetTitle('#sigma_{noise} (ps)')
hdummy8.Draw()
for ii, array in enumerate(arrays):
    g_tRes_tDiff_vs_noise[array].SetLineColor(colors[array])
    g_tRes_tDiff_vs_noise[array].SetMarkerColor(colors[array])
    g_tRes_tDiff_vs_noise[array].SetMarkerStyle(20)
    g_tRes_tDiff_vs_noise[array].Draw('p same')
leg.Draw()



canvas9 = ROOT.TCanvas('c_energy_Milano','c_energy_Milano')
canvas9.SetGridy()
hdummy9 = ROOT.TH2F('hdummy9','', 16, -0.5, 15.5, 100, g_energy_vs_bar['Caltech'].GetMean(2)-6, g_energy_vs_bar['Caltech'].GetMean(2)+6)
hdummy9.GetXaxis().SetTitle('bar')
hdummy9.GetYaxis().SetTitle('energy')
hdummy9.Draw()
array = 'Milano'
g_energyL_vs_bar[array].SetLineColor(colors[array])
g_energyL_vs_bar[array].SetMarkerColor(colors[array])
g_energyL_vs_bar[array].SetMarkerStyle(25)
g_energyL_vs_bar[array].Draw('pl same')
g_energyR_vs_bar[array].SetLineColor(colors[array])
g_energyR_vs_bar[array].SetMarkerColor(colors[array])
g_energyR_vs_bar[array].SetMarkerStyle(26)
g_energyR_vs_bar[array].Draw('pl same')
g_energy_vs_bar[array].SetLineColor(colors[array])
g_energy_vs_bar[array].SetMarkerColor(colors[array])
g_energy_vs_bar[array].SetMarkerStyle(20)
g_energy_vs_bar[array].Draw('pl same')
tl9l[array] = ROOT.TLatex( 0.17, 0.80, 'Left:   mean = %.1f, RMS = %.0f%% '%(g_energyL_vs_bar[array].GetMean(2), 100*g_energyL_vs_bar[array].GetRMS(2)/g_energyL_vs_bar[array].GetMean(2) ) )
tl9l[array].SetNDC()
tl9l[array].SetTextSize(0.030)
tl9l[array].SetTextColor(colors[array])
tl9l[array].Draw()
tl9r[array] = ROOT.TLatex( 0.17, 0.75, 'Right:   mean = %.1f, RMS = %.0f%% '%(g_energyR_vs_bar[array].GetMean(2), 100*g_energyR_vs_bar[array].GetRMS(2)/g_energyR_vs_bar[array].GetMean(2) ) )
tl9r[array].SetNDC()
tl9r[array].SetTextSize(0.030)
tl9r[array].SetTextColor(colors[array])
tl9r[array].Draw()
tl9a[array] = ROOT.TLatex( 0.17, 0.70, 'Average:   mean = %.1f, RMS = %.0f%% '%(g_energy_vs_bar[array].GetMean(2), 100*g_energy_vs_bar[array].GetRMS(2)/g_energy_vs_bar[array].GetMean(2) ) )
tl9a[array].SetNDC()
tl9a[array].SetTextSize(0.030)
tl9a[array].SetTextColor(colors[array])
tl9a[array].Draw()

canvas10 = ROOT.TCanvas('c_energy_Caltech','c_energy_Caltech')
canvas10.SetGridy()
hdummy10 = ROOT.TH2F('hdummy10','', 16, -0.5, 15.5, 100, g_energy_vs_bar['Caltech'].GetMean(2)-6, g_energy_vs_bar['Caltech'].GetMean(2)+6)
hdummy10.GetXaxis().SetTitle('bar')
hdummy10.GetYaxis().SetTitle('energy')
hdummy10.Draw()
array = 'Caltech'
g_energyL_vs_bar[array].SetLineColor(colors[array])
g_energyL_vs_bar[array].SetMarkerColor(colors[array])
g_energyL_vs_bar[array].SetMarkerStyle(25)
g_energyL_vs_bar[array].Draw('pl same')
g_energyR_vs_bar[array].SetLineColor(colors[array])
g_energyR_vs_bar[array].SetMarkerColor(colors[array])
g_energyR_vs_bar[array].SetMarkerStyle(26)
g_energyR_vs_bar[array].Draw('pl same')
g_energy_vs_bar[array].SetLineColor(colors[array])
g_energy_vs_bar[array].SetMarkerColor(colors[array])
g_energy_vs_bar[array].SetMarkerStyle(20)
g_energy_vs_bar[array].Draw('pl same')
tl10l[array] = ROOT.TLatex( 0.17, 0.80, 'Left:   mean = %.1f, RMS = %.0f%% '%(g_energyL_vs_bar[array].GetMean(2), 100*g_energyL_vs_bar[array].GetRMS(2)/g_energyL_vs_bar[array].GetMean(2) ) )
tl10l[array].SetNDC()
tl10l[array].SetTextSize(0.030)
tl10l[array].SetTextColor(colors[array])
tl10l[array].Draw()
tl10r[array] = ROOT.TLatex( 0.17, 0.75, 'Right:   mean = %.1f, RMS = %.0f%% '%(g_energyR_vs_bar[array].GetMean(2), 100*g_energyR_vs_bar[array].GetRMS(2)/g_energyR_vs_bar[array].GetMean(2) ) )
tl10r[array].SetNDC()
tl10r[array].SetTextSize(0.030)
tl10r[array].SetTextColor(colors[array])
tl10r[array].Draw()
tl10a[array] = ROOT.TLatex( 0.17, 0.70, 'Average:   mean = %.1f, RMS = %.0f%% '%(g_energy_vs_bar[array].GetMean(2), 100*g_energy_vs_bar[array].GetRMS(2)/g_energy_vs_bar[array].GetMean(2) ) )
tl10a[array].SetNDC()
tl10a[array].SetTextSize(0.030)
tl10a[array].SetTextColor(colors[array])
tl10a[array].Draw()



canvas11 = ROOT.TCanvas('c_tRes_tDiff_vs_energy','c_tRes_tDiff_vs_energy')
canvas11.SetGridy()
hdummy11 = ROOT.TH2F('hdummy11','', 16, 4, 10, 100, 50, 250)
hdummy11.GetYaxis().SetTitle('#sigma_{tDiff} (ps)')
hdummy11.GetXaxis().SetTitle('energy')
hdummy11.Draw()
for ii, array in enumerate(arrays):
    g_tRes_tDiff_vs_energy[array].SetLineColor(colors[array])
    g_tRes_tDiff_vs_energy[array].SetMarkerColor(colors[array])
    g_tRes_tDiff_vs_energy[array].SetMarkerStyle(20)
    g_tRes_tDiff_vs_energy[array].Draw('p same')
leg.Draw()





canvas12 = ROOT.TCanvas('c_tRes_barL_vs_noise','c_tRes_barL_vs_noise')
canvas12.SetGridy()
hdummy12 = ROOT.TH2F('hdummy12','', 16, 0, 120, 120, 0, 300)
hdummy12.GetYaxis().SetTitle('#sigma_{L-ref} (ps)')
hdummy12.GetXaxis().SetTitle('#sigma_{noise} (ps)')
hdummy12.Draw()
for ii, array in enumerate(arrays):
    g_tRes_barL_vs_noise[array].SetLineColor(colors[array])
    g_tRes_barL_vs_noise[array].SetMarkerColor(colors[array])
    g_tRes_barL_vs_noise[array].SetMarkerStyle(20)
    g_tRes_barL_vs_noise[array].Draw('p same')
leg.Draw()


canvas13 = ROOT.TCanvas('c_tRes_barR_vs_noise','c_tRes_barR_vs_noise')
canvas13.SetGridy()
hdummy13 = ROOT.TH2F('hdummy13','', 16, 0, 120, 120, 0, 300)
hdummy13.GetYaxis().SetTitle('#sigma_{R-ref} (ps)')
hdummy13.GetXaxis().SetTitle('#sigma_{noise} (ps)')
hdummy13.Draw()
for ii, array in enumerate(arrays):
    g_tRes_barR_vs_noise[array].SetLineColor(colors[array])
    g_tRes_barR_vs_noise[array].SetMarkerColor(colors[array])
    g_tRes_barR_vs_noise[array].SetMarkerStyle(20)
    g_tRes_barR_vs_noise[array].Draw('p same')
leg.Draw()


canvas14 = ROOT.TCanvas('c_tRes_barL_noiseSubtracted_vs_noise','c_tRes_barL_noiseSubtracted_vs_noise')
canvas14.SetGridy()
hdummy14 = ROOT.TH2F('hdummy14','', 16, 0, 120, 120, 0, 300)
hdummy14.GetYaxis().SetTitle('#sigma_{L-ref} (ps)')
hdummy14.GetXaxis().SetTitle('#sigma_{noise} (ps)')
hdummy14.Draw()
for ii, array in enumerate(arrays):
    g_tRes_barL_noiseSubtracted_vs_noise[array].SetLineColor(colors[array])
    g_tRes_barL_noiseSubtracted_vs_noise[array].SetMarkerColor(colors[array])
    g_tRes_barL_noiseSubtracted_vs_noise[array].SetMarkerStyle(20)
    g_tRes_barL_noiseSubtracted_vs_noise[array].Draw('p same')
leg.Draw()


canvas15 = ROOT.TCanvas('c_tRes_barR_noiseSubtracted_vs_noise','c_tRes_barR_noiseSubtracted_vs_noise')
canvas15.SetGridy()
hdummy15 = ROOT.TH2F('hdummy15','', 16, 0, 120, 100, 0, 300)
hdummy15.GetYaxis().SetTitle('#sigma_{R-ref} (ps)')
hdummy15.GetXaxis().SetTitle('#sigma_{noise} (ps)')
hdummy15.Draw()
for ii, array in enumerate(arrays):
    g_tRes_barR_noiseSubtracted_vs_noise[array].SetLineColor(colors[array])
    g_tRes_barR_noiseSubtracted_vs_noise[array].SetMarkerColor(colors[array])
    g_tRes_barR_noiseSubtracted_vs_noise[array].SetMarkerStyle(20)
    g_tRes_barR_noiseSubtracted_vs_noise[array].Draw('p same')
leg.Draw()


for c in canvas1, canvas2, canvas3, canvas4, canvas5, canvas5r, canvas6,canvas6r, canvas7, canvas8, canvas9, canvas10, canvas11, canvas12, canvas13, canvas14, canvas15:
    c.SaveAs(outdir+'/'+c.GetName()+'.png')
    c.SaveAs(outdir+'/'+c.GetName()+'.pdf')
    c.SaveAs(outdir+'/'+c.GetName()+'.C')
    
#raw_input('ok?')
