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
import tdrstyle

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetTitleOffset(1.05,'Y')
ROOT.gStyle.SetLabelSize(0.04)
ROOT.gErrorIgnoreLevel = ROOT.kWarning;
ROOT.gROOT.SetBatch(True)
#ROOT.gROOT.SetBatch(False)


parser = argparse.ArgumentParser(description='This script creates the final noise summary plot starting from the moduleCharacterization and drawPulseShape analyses')
parser.add_argument("-r",         "--rangeFile", required=True, type=str, help="input file containing the run ranges to be processed")
parser.add_argument("--pngName",  "--pngName", required=False, default="", type=str, help="name of the output picture")
args = parser.parse_args()


baseFolder    = "/home/data/mtd/RUptoro2/Lab5015Analysis_MTDST_CERN_Oct21"
outFile = ROOT.TFile("/home/data/mtd/RUptoro2/Lab5015Analysis_MTDST_CERN_Oct21/plots_fede/drawLaserCTR_vs_SR.root","RECREATE");
outDir = "/var/www/html/MTDST_CERN_Oct21/CCv2/slewRate/"


dac_to_uA = { 
    'ith2_3': 1.250,
    'ith2_2': 0.940,
    'ith2_1': 0.630,
    'ith2_0': 0.313,
    'ith1_3': 0.630,
    'ith1_2': 0.470,
    'ith1_1': 0.313,
    'ith1_0': 0.156
}

runs_dict = {}
cfgFile = open(baseFolder+'/'+args.rangeFile)
lines = cfgFile.readlines()

for line in lines:
   if "#" in line:
      continue

   line = line.strip()
   if line == "":
      continue

   line = line.split()

   if "ch1" in line:
      ch1 = line[1]
      continue
   if "ch2" in line:
      ch2 = line[1]
      continue
   if "pngName" in line:
      pngName = line[1]
      continue


   runs_dict[line[0]] = [line[1], float(line[2]), list(map(float, line[3].split(",")))]


VovList = []
for run_range, value in runs_dict.items():
    VovList += value[2]

VovList = sorted(set(VovList))


print(VovList)

g_tRes_vs_SR = {}
g_tRes_vs_SR_all = ROOT.TGraphErrors()
for Vov in VovList:
    g_tRes_vs_SR[Vov] = ROOT.TGraphErrors()


g_SRch2_vs_SRch1 = ROOT.TGraphErrors()



def get_list(input_list):
    if input_list is None:
        return None
    vals = []
    comma_list = input_list.split(',')
    for item in comma_list:
        hyphen_list = item.split('-')
        if len(hyphen_list) > 1:
            for i in range(int(hyphen_list[0]), int(hyphen_list[1])+1):
                vals.append(i)
        else:
            vals.append(hyphen_list[0])
    return vals


#---------------------------
# -- compute SR vs threshold
npoints = 15
nintervals = 4
pngLabel = ""

for it,run in enumerate(sorted(runs_dict)):
        print("===>>> Processing run range: ", run)

        fullPath = outDir+"/"+run
        if (os.path.isdir(fullPath) == False):
            os.system('mkdir -p %s'%fullPath)

    
        pngLabel += str(run)+"_"

        ithMode = runs_dict[run][0]
        laserTune = runs_dict[run][1]
        Vovs  = runs_dict[run][2]

        for Vov in Vovs:
            inFile = ROOT.TFile.Open('/home/data/mtd/RUptoro2/Lab5015Analysis_MTDST_CERN_Oct21/plots_fede/pulseShape_run%s.root'%run)
    
            SR = 0.
            SR_errSum = 0.
            SRsingleCh = []
            thresh = 0.
            channels = ['ch1', 'ch2']

            graph1 = inFile.Get('g_ps_totSel_%s_Vov%.01f'%('ch1',Vov))
            graph2 = inFile.Get('g_ps_totSel_%s_Vov%.01f'%('ch2',Vov))
            if graph1 == None or graph2 == None:
                continue
        
            index_cen = int( min(graph1.GetN()/nintervals, graph2.GetN()/nintervals)) -1
            index_min = max(1, index_cen - int(npoints/2)) -1
            index_max = min(index_cen + int(npoints/2), int(min(graph1.GetN()/2, graph2.GetN()/2))) -1

            for ch in channels:
                canvas = ROOT.TCanvas()
                canvas.cd()
                graph = inFile.Get('g_ps_totSel_%s_Vov%.01f'%(ch,Vov))
                graph.Draw("AP")

                fitSR_inverted = ROOT.TF1('fitSR_inverted', 'pol1')
                fitSR = ROOT.TF1('fitSR', 'pol1',-10.,10.)
                fitSR.SetNpx(10000)

                invertedGraph = ROOT.TGraphErrors()
                partialGraph = ROOT.TGraphErrors()
                for iPoint in range(index_min,index_max+1):
                    invertedGraph.SetPoint(iPoint-index_min, graph.GetPointY(iPoint), graph.GetPointX(iPoint))
                    invertedGraph.SetPointError(iPoint-index_min, 0., 0.015)
                    partialGraph.SetPoint(iPoint-index_min, graph.GetPointX(iPoint), graph.GetPointY(iPoint))

                partialGraph.Draw("P,sames")

                fitSR_inverted.SetParameters(0., 1e-3)
                invertedGraph.Fit(fitSR_inverted,'QS')

                fitSR.SetParameters( -fitSR_inverted.GetParameter(0)/fitSR_inverted.GetParameter(1), 1./fitSR_inverted.GetParameter(1) )
                errX = fitSR_inverted.GetParError(1) / fitSR_inverted.GetParameter(1) / fitSR_inverted.GetParameter(1)
                #print("ERROR: ",fitSR_inverted.GetParError(1), errX)
                fitSR.SetParError(1, errX)
                fitSR.Draw("sames")

                canvas.Print(fullPath+'/ps_totSel_%s_Vov%.01f.png'%(ch,Vov))
                outFile.cd()
                graph.Write('Run%s_g_ps_totSel_%s_Vov%.01f'%(run,ch,Vov))
                invertedGraph.Write('Run%s_g_ps_inverted_totSel_%s_Vov%.01f'%(run,ch,Vov))
                
                print(ch, index_cen, index_min, index_max, graph.GetX()[index_cen], graph.GetX()[index_min], graph.GetX()[index_max], fitSR.GetParameter(1))
                SR += fitSR.GetParameter(1) * 1./errX/errX
                SR_errSum += 1./errX/errX

                SRsingleCh.append(SR)
                thresh += int(round(graph.GetPointY(index_cen)/dac_to_uA[ithMode]))

            g_SRch2_vs_SRch1.SetPoint(it, SRsingleCh[0], SRsingleCh[1])
        
            SR /= SR_errSum
            thresh /= 2.
    
    
            inFile = ROOT.TFile.Open('/home/data/mtd/RUptoro2/Lab5015Analysis_MTDST_CERN_Oct21/plots_fede/moduleCharacterization_step2_run%s.root'%run)
        
            print('NAME: h1_deltaT_totRatioCorr_bar00L-R_Vov%.02f_th%.02d_energyBin01'%(Vov,thresh))
            histo = inFile.Get('h1_deltaT_totRatioCorr_bar00L-R_Vov%.02f_th%.02d_energyBin01'%(Vov,thresh))
            if histo == None:
                print("MISSING HISTO")
                continue
            if histo.GetEntries() < 100:
                continue
            if histo.GetRMS() <= 0.:
                continue

            canvas = ROOT.TCanvas()
            histo.Draw()
            fitFunc = ROOT.TF1('fitFunc','gaus',-10000, 10000)
            fitFunc.SetRange(histo.GetMean()-2.0*histo.GetRMS(),histo.GetMean()+2.0*histo.GetRMS())
            histo.Fit('fitFunc','QNRS+')
            fitFunc.SetRange(fitFunc.GetParameter(1)-3.0*fitFunc.GetParameter(2),fitFunc.GetParameter(1)+3.0*fitFunc.GetParameter(2))
            histo.Fit('fitFunc','QNS+','', fitFunc.GetParameter(1)-3.0*fitFunc.GetParameter(2),fitFunc.GetParameter(1)+3.0*fitFunc.GetParameter(2))
            res = [fitFunc.GetParameter(2), fitFunc.GetParError(2)]
            canvas.Print(fullPath+'/tRes_%s_thr%s_Vov%.01f.png'%(ch,thresh,Vov))

            print 'runRange = ', run, 'bestTh = ', thresh, '   time res = ',  res[0]/math.sqrt(2), '  SR =', SR
            g_tRes_vs_SR[Vov].SetPoint(g_tRes_vs_SR[Vov].GetN(), SR, res[0]/math.sqrt(2))
            g_tRes_vs_SR[Vov].SetPointError(g_tRes_vs_SR[Vov].GetN()-1, 1./math.sqrt(SR_errSum)/2., res[1]/math.sqrt(2))
            
            g_tRes_vs_SR_all.SetPoint(g_tRes_vs_SR_all.GetN(), SR, res[0]/math.sqrt(2))
            g_tRes_vs_SR_all.SetPointError(g_tRes_vs_SR_all.GetN()-1, 1./math.sqrt(SR_errSum)/2., res[1]/math.sqrt(2))
            print("   ")
        




b = ROOT.TCanvas('c_SRch2_vs_SRch1','c_SRch2_vs_SRch1',1200,700)
bPad = ROOT.gPad.DrawFrame(0.,0.,500.,500.)
bPad.SetTitle(";slew rate ch1 [#muA/ns];slew rate ch2 [#muA/ns]");
bPad.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();
g_SRch2_vs_SRch1.Draw("psame")



c = ROOT.TCanvas('c_tRes_vs_SR','c_tRes_vs_SR',1200,700)
hPad = ROOT.gPad.DrawFrame(0.,0.,300.,100.)
hPad.SetTitle(";slew rate [#muA/ns];#sigma_{t}^{single} [ps]");
hPad.Draw();
ROOT.gPad.SetGridx();
ROOT.gPad.SetGridy();

legLow = 0.94 - len(VovList)*0.03
legend = ROOT.TLegend(0.60, legLow, 0.89, 0.94)
legend.SetBorderSize(0)

fit = {}
it = 1
colors = {1:ROOT.kRed-4, 2:ROOT.kMagenta-4, 3:ROOT.kBlue-4, 4:ROOT.kCyan, 5:ROOT.kGreen+1, 6:ROOT.kOrange, 7:ROOT.kYellow-3, 8:ROOT.kBlack}
latexes = {}
for Vov in VovList:
    #g_tRes_vs_SR[Vov].SetLineColor(ROOT.kRainBow+3*it)
    #g_tRes_vs_SR[Vov].SetMarkerColor(ROOT.kRainBow+3*it)
    g_tRes_vs_SR[Vov].SetLineColor(colors[it])
    g_tRes_vs_SR[Vov].SetMarkerColor(colors[it])
    g_tRes_vs_SR[Vov].SetMarkerSize(1)
    g_tRes_vs_SR[Vov].Draw('psame')
    legend.AddEntry(g_tRes_vs_SR[Vov], 'V_{OV} = %.1f V'%Vov,'PL')
    fit[Vov] = ROOT.TF1('fit%d'%laserTune,'sqrt([0]*[0] + [1]*[1]/x/x )', 0, 1000)
    fit[Vov].SetLineColor(colors[it])
    fit[Vov].SetLineWidth(1)
    fit[Vov].SetParameters(12, 500)
    g_tRes_vs_SR[Vov].Fit(fit[Vov],'QRNS')
    fit[Vov].Draw('same')

    legY = 0.89 - it*0.04
    latexes[Vov] = ROOT.TLatex(0.55,legY,'V_{OV} = %.1f V:   #sigma_{t} = %.02f #muA / (dV/dt) #oplus %.01f ps'%(Vov,fit[Vov].GetParameter(1)/1000.,fit[Vov].GetParameter(0)))
    latexes[Vov].SetNDC()
    latexes[Vov].SetTextFont(42)
    latexes[Vov].SetTextSize(0.035)
    latexes[Vov].SetTextColor(colors[it])
    latexes[Vov].Draw('same')

    it += 1

#legend.Draw('same')


fitAll = ROOT.TF1('fitAll','sqrt([0]*[0] + [1]*[1]/x/x )', 0, 1000)
fitAll.SetLineColor(ROOT.kBlack)
fitAll.SetLineStyle(2)
fitAll.SetParameters(12, 500)
g_tRes_vs_SR_all.Fit(fitAll,'SR')
fitAll.Draw('same')

latex = ROOT.TLatex(0.55,0.89,'#sigma_{t} = %.02f #pm %.02f #muA / (dV/dt) #oplus %.01f #pm %.01f ps'%(fitAll.GetParameter(1)/1000., fitAll.GetParError(1)/1000., fitAll.GetParameter(0), fitAll.GetParError(0)))
latex.SetNDC()
latex.SetTextFont(42)
latex.SetTextSize(0.04)
latex.SetTextColor(ROOT.kBlack)
latex.Draw('same')

coordinates = pngName.strip(".png").split("_")
latex2 = ROOT.TLatex(0.63,0.965,coordinates[0]+"   "+coordinates[1]+"   "+coordinates[2])
latex2.SetNDC()
latex2.SetTextFont(42)
latex2.SetTextSize(0.03)
latex2.SetTextColor(ROOT.kBlack)
latex2.Draw('same')


if args.pngName != "":
    pngLabel = args.pngName
    c.Print(outDir+"/tRes_vs_SR_"+pngLabel)
    b.Print(outDir+"/SRch2_vs_SRc1_"+pngLabel)


else:
    pngLabel = pngLabel[:-1]
    c.Print(outDir+"/tRes_vs_SR_"+pngLabel+".png")
    b.Print(outDir+"/SRch2_vs_SRc1_"+pngLabel+".png")

