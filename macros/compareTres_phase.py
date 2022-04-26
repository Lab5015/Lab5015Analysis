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
ROOT.gROOT.SetBatch(False)
ROOT.gErrorIgnoreLevel = ROOT.kWarning




phases = ['300-350',
          '350-400',
          '400-450',
          '450-500',
          '500-550',
          '550-600',
          '600-650',
          '650-700',
          '700-750',
          '750-800',
          '300-800'
]


c = ROOT.TCanvas('c','c',600,600)
c.SetGridy()
c.SetGridx()
hdummy = ROOT.TH2F('hdummy','',100,0,40, 150, 30, 180)
hdummy.GetXaxis().SetTitle('threshold [DAC]')
hdummy.GetYaxis().SetTitle('#sigma_{t} [ps]')
hdummy.Draw('')

leg = ROOT.TLegend(0.20,0.65,0.55,0.89)

for i,phase in enumerate(phases):
    f = ROOT.TFile.Open('/afs/cern.ch/work/m/malberti/MTD/TBatH8Oct2021/Lab5015Analysis/plots/HPK528_unirr_52deg_T10C_Vov1.50_testPhase%s.root'%phase)
    g = f.Get('g_deltaT_energyRatioCorr_vs_th_bar06_Vov1.50_enBin01')
    if (g == None): continue
    g.SetMarkerSize(1)
    if (phase == '300-800'):
        g.SetMarkerStyle(21)
        g.SetMarkerColor(1)
        g.SetLineColor(1)
    else:
        g.SetMarkerStyle(20)
        g.SetMarkerColor(51+i*5)
        g.SetLineColor(51+i*5)
    leg.AddEntry(g,'%s'%phase,'PL')
    g.Draw('plsame')

leg.Draw('same')
c.SaveAs('c.png')

raw_input('ok')
