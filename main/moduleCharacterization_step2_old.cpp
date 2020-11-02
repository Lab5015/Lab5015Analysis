#include "interface/AnalysisUtils.h"
#include "interface/Na22SpectrumAnalyzer.h"
#include "interface/FitUtils.h"
#include "interface/SetTDRStyle.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <time.h>
#include <stdio.h>
#include <sys/stat.h>

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TSpectrum.h"





int main(int argc, char** argv)
{
  setTDRStyle();
  
  
  if( argc < 2 )
  {
    std::cout << ">>> moduleCharacterization_step2::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }
  
  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  
  
  //--- get parameters
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  system(Form("mkdir -p %s",plotDir.c_str()));
  system(Form("mkdir -p %s/qfine/",plotDir.c_str()));
  system(Form("mkdir -p %s/tot/",plotDir.c_str()));
  system(Form("mkdir -p %s/energy/",plotDir.c_str()));
  system(Form("mkdir -p %s/energyRatio/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR_energyRatioCorr/",plotDir.c_str()));
  
  system(Form("cp %s/../index.php %s",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/qfine/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/tot/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/energy/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/energyRatio/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/CTR/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/CTR_energyRatioCorr/",plotDir.c_str(),plotDir.c_str()));
  
  
  std::vector<std::string> LRLabels;
  LRLabels.push_back("L");
  LRLabels.push_back("R");
  
  
  // //--- get cuts per bar / Vov
  // std::map<unsigned int,std::map<float,float> > cut_qfineAcc;
  // std::map<unsigned int,std::map<float,float> > cut_totAcc;
  // std::map<unsigned int,std::map<float,float> > cut_energyAcc;
  // std::map<unsigned int,std::map<float,float> > cut_energyFitMin;
  // std::map<unsigned int,std::map<float,float> > cut_energyFitMax;
  // std::map<unsigned int, float > noise;
  // for(auto ch :  channels)
  // {
  //   int chID = opts.GetOpt<int>(Form("%s.chID",ch.c_str()));
  //   std::vector<float> Vovs          = opts.GetOpt<std::vector<float> >(Form("%s.Vovs",ch.c_str()));
  //   std::vector<float> qfineMins     = opts.GetOpt<std::vector<float> >(Form("%s.qfineMins",ch.c_str()));
  //   std::vector<float> totMins       = opts.GetOpt<std::vector<float> >(Form("%s.totMins",ch.c_str()));
  //   std::vector<float> energyMins    = opts.GetOpt<std::vector<float> >(Form("%s.energyMins",ch.c_str()));
  //   std::vector<float> energyFitMins = opts.GetOpt<std::vector<float> >(Form("%s.energyFitMins",ch.c_str()));
  //   std::vector<float> energyFitMaxs = opts.GetOpt<std::vector<float> >(Form("%s.energyFitMaxs",ch.c_str()));
  //   noise[chID] = opts.GetOpt<float>(Form("%s.noise",ch.c_str()));
  //   int iter = 0;
  //   for(auto Vov : Vovs)
  //   {
  //     cut_qfineAcc[chID][Vov]     = qfineMins.at(iter);
  //     cut_totAcc[chID][Vov]       = totMins.at(iter);
  //     cut_energyAcc[chID][Vov]    = energyMins.at(iter);
  //     cut_energyFitMin[chID][Vov] = energyFitMins.at(iter);
  //     cut_energyFitMax[chID][Vov] = energyFitMaxs.at(iter);
  //     ++iter;
  //   }
  // }
  // std::map<std::string,float> cut_energyMin;
  // std::map<std::string,float> cut_energyMax;
  
  
  //--- open files
  std::string step1FileName= opts.GetOpt<std::string>("Input.step1FileName");
  TFile* inFile = TFile::Open(step1FileName.c_str(),"READ");
  
  std::map<std::string,TTree*> trees;
  
  std::map<std::string,int> VovLabels;
  std::map<std::string,int> thLabels;
  std::vector<std::string> stepLabels;
  std::map<std::string,float> map_Vovs;
  std::map<std::string,float> map_ths;  
  
  TList* list = inFile -> GetListOfKeys();
  TIter next(list);
  TObject* object = 0;
  while( (object = next()) )
  {
    std::string name(object->GetName());
    std::vector<std::string> tokens = GetTokens(name,'_');
    std::cout << name << std::endl;
    std::size_t found;
    
    found = name.find("data_");
    if( found!=std::string::npos )
    {
      std::string label(Form("%s_%s_%s",tokens[1].c_str(),tokens[2].c_str(),tokens[3].c_str()));
      trees[label] = (TTree*)( inFile->Get(name.c_str()) );
    }
    
    found = name.find("h1_energy");
    if( found!=std::string::npos )
    {
      std::string stepLabel = tokens[3]+"_"+tokens[4];
      VovLabels[tokens[3]] += 1;
      thLabels[tokens[4]] += 1;
      stepLabels.push_back(stepLabel);
      
      std::string string_Vov = tokens[3];
      string_Vov.erase(0,3);
      map_Vovs[stepLabel] = atof(string_Vov.c_str());
      
      std::string string_th = tokens[4];
      string_th.erase(0,2);
      map_ths[stepLabel] = atof(string_th.c_str());
    }
  }
  std::sort(stepLabels.begin(),stepLabels.end());
  stepLabels.erase(std::unique(stepLabels.begin(),stepLabels.end()),stepLabels.end());
  
  
  //--- define histograms
  std::string plotFileName = opts.GetOpt<std::string>("Output.plotFileName");
  TFile* outFile = TFile::Open(plotFileName.c_str(),"RECREATE");
  
  std::map<std::string,TH1F*> h1_energyRatio;
  
  std::map<std::string,TH1F*> h1_deltaT_raw;
  std::map<std::string,TH1F*> h1_deltaT;
  std::map<std::string,TProfile*> p1_deltaT_vs_energyRatio;

  std::map<std::string,TH1F*> h1_deltaT_energyRatioCorr;
  
  std::map<std::string,std::vector<float>*> ranges;
  std::map<std::string,std::map<std::string,std::pair<float,float> > > peaks;
  
  
  
  //--- get plot settings
  float energyMin = opts.GetOpt<float>("Plots.energyMin");
  float energyMax = opts.GetOpt<float>("Plots.energyMax");
  float tResMin = opts.GetOpt<float>("Plots.tResMin");
  float tResMax = opts.GetOpt<float>("Plots.tResMax");
  int tResMode = opts.GetOpt<int>("Plots.tResMode");
  float refVov = opts.GetOpt<float>("Plots.refVov");
  float refTh = opts.GetOpt<float>("Plots.refTh");
  
  float corr = 1.;
  if( tResMode == 1 ) corr = 1./sqrt(2.);
  if( tResMode == 2 ) corr = 0.5;
  
  
  
  TCanvas* c;
  float* vals = new float[6];
  TLatex* latex;
  TH1F* histo;
  TProfile* prof;
  TLegend* legend;
  TLegend* legend2;
  TH1F* hPad;
  
  TLatex* latex_ref = new TLatex(0.20,0.85,Form("#splitline{V_{OV} = %.1f V}{th. = %.0f DAC}",refVov,refTh));
  latex_ref -> SetNDC();
  latex_ref -> SetTextFont(42);
  latex_ref -> SetTextSize(0.04);
  
  
  
  
  
  
  //------------------
  //--- draw 1st plots
  std::map<std::string,TGraphErrors*> g_energy_vs_th;
  std::map<std::string,TGraphErrors*> g_energy_vs_Vov;
  
  for(auto stepLabel : stepLabels)
  {
    float Vov = map_Vovs[stepLabel];
    float th = map_ths[stepLabel];
    std::string VovLabel(Form("Vov%.1f",Vov));
    std::string thLabel(Form("th%02.0f",th));
    
    
    //--------------------------------------------------------
    
    
    for(int iBar = 0; iBar < 16; ++iBar)
    {
      for(auto LRLabel : LRLabels )
      {
        std::string label(Form("bar%02d%s_%s",iBar,LRLabel.c_str(),stepLabel.c_str()));
        std::string label_vs_th(Form("bar%02d%s_%s",iBar,LRLabel.c_str(),VovLabel.c_str()));
        std::string label_vs_Vov(Form("bar%02d%s_%s",iBar,LRLabel.c_str(),thLabel.c_str()));
        
        histo = (TH1F*)( inFile->Get(Form("h1_qfine_%s",label.c_str())) );
        if( !histo ) continue;
        if( histo->GetEntries() < 100 ) continue;
        
        
        latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d%s}{V_{OV} = %.1f V, th. = %d DAC}",iBar,LRLabel.c_str(),Vov,int(th)));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kRed);
        
        
        c = new TCanvas(Form("c_qfine_%s",label.c_str()),Form("c_qfine_%s",label.c_str()));
        gPad -> SetLogy();
        
        histo = (TH1F*)( inFile->Get(Form("h1_qfine_%s",label.c_str())) );
        histo -> SetTitle(";Q_{fine} [ADC];entries");
        histo -> SetLineColor(kRed);
        histo -> Draw();
        histo -> GetXaxis() -> SetRangeUser(0.,600.);
        // TLine* line_qfineAcc1 = new TLine(cut_qfineAcc[chID][Vov],histo->GetMinimum(),cut_qfineAcc[chID][Vov],histo->GetMaximum());
        // line_qfineAcc1 -> SetLineColor(kBlack);
        // line_qfineAcc1 -> Draw("same");
        latex -> Draw("same");      
        histo -> Write();
        c -> Print(Form("%s/qfine/c_qfine__%s.png",plotDir.c_str(),label.c_str()));
        c -> Print(Form("%s/qfine/c_qfine__%s.pdf",plotDir.c_str(),label.c_str()));
        delete c;
        
        
        c = new TCanvas(Form("c_tot_%s",label.c_str()),Form("c_tot_%s",label.c_str()));
        // gPad -> SetLogy();
        
        histo = (TH1F*)( inFile->Get(Form("h1_tot_%s",label.c_str())) );
        histo -> SetTitle(";ToT [ns];entries");
        histo -> SetLineColor(kRed);
        histo -> Draw();
        // TLine* line_totAcc1 = new TLine(cut_totAcc[chID][Vov],histo->GetMinimum(),cut_totAcc[chID][Vov],histo->GetMaximum());
        // line_totAcc1 -> SetLineColor(kBlack);
        // line_totAcc1 -> Draw("same");
        latex -> Draw("same");      
        histo -> Write();
        c -> Print(Form("%s/tot/c_tot__%s.png",plotDir.c_str(),label.c_str()));
        c -> Print(Form("%s/tot/c_tot__%s.pdf",plotDir.c_str(),label.c_str()));
        delete c;
        
        
        c = new TCanvas(Form("c_energy_%s",label.c_str()),Form("c_energy_%s",label.c_str()));
        gPad -> SetLogy();
        
        histo = (TH1F*)( inFile->Get(Form("h1_energy_%s",label.c_str())) );      
        histo -> SetTitle(";energy [a.u.];entries");
        histo -> SetLineColor(kRed);
        histo -> SetLineWidth(2);
        histo -> Draw();
        
        peaks[label] = Na22SpectrumAnalyzer(histo);
        
        for(auto peak : peaks[label] )
        {
          histo -> GetXaxis() -> SetRangeUser(0.,5.*peak.second.first);
          break;
        }
        histo -> GetYaxis() -> SetRangeUser(3.,5.*histo->GetMaximum());
        
        latex -> Draw("same");      
        histo -> Write();
        c -> Print(Form("%s/energy/c_energy__%s.png",plotDir.c_str(),label.c_str()));
        c -> Print(Form("%s/energy/c_energy__%s.pdf",plotDir.c_str(),label.c_str()));
        delete c;
        
        
        for(auto peak : peaks[label])
        {
          std::string peakStr(Form("_%s",peak.first.c_str()));
          
          if( g_energy_vs_th[label_vs_th+peakStr] == NULL ) g_energy_vs_th[label_vs_th+peakStr] = new TGraphErrors();
          if( g_energy_vs_Vov[label_vs_Vov+peakStr] == NULL ) g_energy_vs_Vov[label_vs_Vov+peakStr] = new TGraphErrors();
          
          g_energy_vs_th[label_vs_th+peakStr] -> SetPoint(g_energy_vs_th[label_vs_th+peakStr]->GetN(),th,peak.second.first);
          g_energy_vs_th[label_vs_th+peakStr] -> SetPointError(g_energy_vs_th[label_vs_th+peakStr]->GetN()-1,0.,0.);
          
          g_energy_vs_Vov[label_vs_Vov+peakStr] -> SetPoint(g_energy_vs_Vov[label_vs_Vov+peakStr]->GetN(),Vov,peak.second.first);
          g_energy_vs_Vov[label_vs_Vov+peakStr] -> SetPointError(g_energy_vs_Vov[label_vs_Vov+peakStr]->GetN()-1,0.,0.);
        }
      }
      
      
      {
        std::string label(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
        
        latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.1f V, th. = %d DAC}",iBar,Vov,int(th)));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kRed);
        
        
        c = new TCanvas(Form("c_energy_%s",label.c_str()),Form("c_energy_%s",label.c_str()));
        gPad -> SetLogy();
        
        histo = (TH1F*)( inFile->Get(Form("h1_energy_%s",label.c_str())) );
        if( !histo ) continue;
        if( histo->GetEntries() < 100 ) continue;
        
        histo -> SetTitle(";energy [a.u.];entries");
        histo -> SetLineWidth(2);
        histo -> SetLineColor(kRed);
        histo -> Draw();
        
        ranges[label] = new std::vector<float>;
        peaks[label] = Na22SpectrumAnalyzer(histo,ranges[label]);
        
        for(auto peak : peaks[label] )
        {
          histo -> GetXaxis() -> SetRangeUser(0.,5.*peak.second.first);
          break;
        }
        histo -> GetYaxis() -> SetRangeUser(3.,5.*histo->GetMaximum());
        
        latex -> Draw("same");
        histo -> Write();
        c -> Print(Form("%s/energy/c_energy__%s.png",plotDir.c_str(),label.c_str()));
        c -> Print(Form("%s/energy/c_energy__%s.pdf",plotDir.c_str(),label.c_str()));
        delete c;
      }
    }
  }
  
  
  for(int iBar = 0; iBar < 16; ++iBar)
  {
    c = new TCanvas(Form("c_energy_vs_th_bar%02d_L-R",iBar),Form("c_energy_vs_th_bar%02d_L-R",iBar));
    
    hPad = (TH1F*)( gPad->DrawFrame(-1.,energyMin,64.,energyMax) );
    hPad -> SetTitle(";threshold [DAC];energy [a.u.]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    int iter = 0;
    bool skip = true;
    for(auto mapIt : VovLabels)
    {
      int nPeaks = peaks[Form("bar%02dL_%s_th%02.0f",iBar,mapIt.first.c_str(),refTh)].size();
      unsigned int size = VovLabels.size()*nPeaks;
      
      int iPeak = 0;
      for( auto peak : peaks[Form("bar%02dL_%s_th%02.0f",iBar,mapIt.first.c_str(),refTh)])
      {
        std::string peakStr(Form("_%s",peak.first.c_str()));
        
        std::string label1(Form("bar%02dL_%s%s",iBar,mapIt.first.c_str(),peakStr.c_str()));
        std::string label2(Form("bar%02dR_%s%s",iBar,mapIt.first.c_str(),peakStr.c_str()));
        TGraph* g_energyL = g_energy_vs_th[label1];
        TGraph* g_energyR = g_energy_vs_th[label2];
        if( !g_energyL || !g_energyR ) continue;
        if( g_energyL->GetN() == 0 || g_energyR->GetN() == 0 ) continue;
        
        g_energyL -> SetLineColor(51+iter*(int(50/size)));
        g_energyL -> SetMarkerColor(51+iter*(int(50/size)));
        g_energyL -> SetMarkerStyle(20+iPeak);
        g_energyL -> Draw("PL,same");
        
        g_energyR -> SetLineColor(51+iter*(int(50/size)));
        g_energyR -> SetMarkerColor(51+iter*(int(50/size)));
        g_energyR -> SetMarkerStyle(24+iPeak);
        g_energyR -> SetLineStyle(7);
        g_energyR -> Draw("PL,same");
        
        latex = new TLatex(0.45,0.85-0.04*iter,Form("%s - %s",peak.first.c_str(),mapIt.first.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(51+iter*(int(50/size)));
        latex -> Draw("same");
        
        skip = false;
        ++iPeak;
        ++iter;
      }
    }
    
    if( skip ) continue;
    // legend -> Draw("same");
    
    latex = new TLatex(0.25,0.85,Form("bar %02d",iBar));
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    latex -> Draw("same");
    
    c -> Print(Form("%s/c_energy_vs_th__bar%02d_L-R.png",plotDir.c_str(),iBar));
    c -> Print(Form("%s/c_energy_vs_th__bar%02d_L-R.pdf",plotDir.c_str(),iBar));
    delete c;
    
    
    c = new TCanvas(Form("c_energy_vs_Vov_bar%02d_L-R",iBar),Form("c_energy_vs_Vov_bar%02d_L-R",iBar));
    
    hPad = (TH1F*)( gPad->DrawFrame(0.,energyMin,10.,energyMax) );
    hPad -> SetTitle(";V_{OV} [V];energy [a.u.]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    iter = 0;
    skip = true;
    for(auto mapIt : thLabels)
    {
      unsigned int nPeaks = peaks[Form("bar%02dL_Vov%.1f_%s",iBar,refVov,mapIt.first.c_str())].size();
      unsigned int size = thLabels.size()*nPeaks;
      
      int iPeak = 0;
      for( auto peak : peaks[Form("bar%02dL_Vov%.1f_%s",iBar,refVov,mapIt.first.c_str())])
      {
        std::string peakStr(Form("_%s",peak.first.c_str()));
        std::string label1(Form("bar%02dL_%s%s",iBar,mapIt.first.c_str(),peakStr.c_str()));
        std::string label2(Form("bar%02dR_%s%s",iBar,mapIt.first.c_str(),peakStr.c_str()));
        TGraph* g_energyL = g_energy_vs_Vov[label1];
        TGraph* g_energyR = g_energy_vs_Vov[label2];
        if( !g_energyL || !g_energyR ) continue;
        if( g_energyL->GetN() == 0 || g_energyR->GetN() == 0 ) continue;
        
        g_energyL -> SetLineColor(51+iter*(int(50/size)));
        g_energyL -> SetMarkerColor(51+iter*(int(50/size)));
        g_energyL -> SetMarkerStyle(20+iPeak);
        g_energyL -> Draw("PL,same");
        
        g_energyR -> SetLineColor(51+iter*(int(50/size)));
        g_energyR -> SetMarkerColor(51+iter*(int(50/size)));
        g_energyR -> SetMarkerStyle(24+iPeak);
        g_energyR -> SetLineStyle(7);
        g_energyR -> Draw("PL,same");
        
        latex = new TLatex(0.45,0.85-0.04*iter,Form("%s - %s",peak.first.c_str(),mapIt.first.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(51+iter*(int(50/size)));
        latex -> Draw("same");
        
        skip = false;
        ++iPeak;
        ++iter;
      }
    }
    
    if( skip ) continue;
    // legend -> Draw("same");
    
    c -> Print(Form("%s/c_energy_vs_Vov__bar%02d_L-R.png",plotDir.c_str(),iBar));
    c -> Print(Form("%s/c_energy_vs_Vov__bar%02d_L-R.pdf",plotDir.c_str(),iBar));
    delete c;
  }
  
  
  c = new TCanvas(Form("c_energy_vs_bar"),Form("c_energy_vs_bar"));
  
  hPad = (TH1F*)( gPad->DrawFrame(-1.,energyMin,16.,energyMax) );
  hPad -> SetTitle(";bar ID;energy [a.u.]");
  hPad -> Draw();
  gPad -> SetGridy();
  
  unsigned int nPeaks = peaks[Form("bar%02dL_Vov%.1f_th%02.0f",7,refVov,refTh)].size();
  int iter = 0;
  
  legend = new TLegend(0.45,0.90-0.04*nPeaks*2,0.90,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  
  int iPeak = 0;
  for(auto peak : peaks[Form("bar%02dL_Vov%.1f_th%02.0f",7,refVov,refTh)])
  {
    TGraph* g_1 = new TGraph();
    TGraph* g_2 = new TGraph();
    
    for(int iBar = 0; iBar < 16; ++iBar)
    {
      std::string peakStr(Form("_%s",peak.first.c_str()));
      
      std::string label1(Form("bar%02dL_Vov%.1f%s",iBar,refVov,peakStr.c_str()));
      std::string label2(Form("bar%02dR_Vov%.1f%s",iBar,refVov,peakStr.c_str()));
      TGraph* g_energyL = g_energy_vs_th[label1];
      TGraph* g_energyR = g_energy_vs_th[label2];
      if( !g_energyL || !g_energyR ) continue;
      if( g_energyL->GetN() == 0 || g_energyR->GetN() == 0 ) continue;
      
      double x,y;
      
      for(int point = 0; point < g_energyL->GetN(); ++point)
      {
        g_energyL -> GetPoint(point,x,y);
        if( x == refTh ) break;
      }
      g_1 -> SetPoint(g_1->GetN(),iBar,y);
      
      for(int point = 0; point < g_energyR->GetN(); ++point)
      {
        g_energyR -> GetPoint(point,x,y);
        if( x == refTh ) break;
      }
      g_2 -> SetPoint(g_2->GetN(),iBar,y);
    }
    
    g_1 -> SetLineColor(51+iter*(int(50/nPeaks)));
    g_1 -> SetMarkerColor(51+iter*(int(50/nPeaks)));
    g_1 -> SetMarkerStyle(20+iPeak);
    g_1 -> Draw("PL,same");
    
    g_2 -> SetLineColor(51+iter*(int(50/nPeaks)));
    g_2 -> SetMarkerColor(51+iter*(int(50/nPeaks)));
    g_2 -> SetMarkerStyle(24+iPeak);
    g_2 -> SetLineStyle(7);
    g_2 -> Draw("PL,same");
    
    legend -> AddEntry(g_1,Form("%s - left",peak.first.c_str()),"PL");
    legend -> AddEntry(g_2,Form("%s - right",peak.first.c_str()),"PL");
    
    ++iPeak;
    ++iter;
  }
  
  latex_ref -> Draw("same");
  legend -> Draw("same");
  
  c -> Print(Form("%s/c_energy_vs_bar.png",plotDir.c_str()));
  c -> Print(Form("%s/c_energy_vs_bar.pdf",plotDir.c_str()));
  delete c;
  
  // end 1st plots
  
  
  
  
  
  
  
  //------------------------
  //--- 2nd loop over events
  std::map<std::string,std::map<int,bool> > accept;
  
  for(auto mapIt : trees)
  {
    std::string label = mapIt.first;
    
    ModuleEventClass* anEvent = new ModuleEventClass();
    mapIt.second -> SetBranchAddress("event",&anEvent);
    
    int nEntries = mapIt.second->GetEntries();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 2nd loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
      mapIt.second -> GetEntry(entry);
      
      accept[label][entry] = false;
      
      if( !ranges[label] ) continue;
      int energyBin = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges[label]);
      if( energyBin < 0 ) continue;
      
      accept[label][entry] = true;
      
      std::string labelLR(Form("bar%02dL-R_Vov%.1f_th%02d",anEvent->barID,anEvent->Vov,anEvent->vth1));
      std::string label_energyBin(Form("%s_energyBin%d",labelLR.c_str(),energyBin));
      
      if( h1_energyRatio[label_energyBin] == NULL )
      {
        h1_energyRatio[label_energyBin] = new TH1F(Form("h1_energyRatio_%s",label_energyBin.c_str()),"",1000,0.,5.);
        h1_deltaT_raw[label_energyBin] = new TH1F(Form("h1_deltaT_raw_%s",label_energyBin.c_str()),"",1250,-5000.,5000.);
      }
      
      h1_energyRatio[label_energyBin] -> Fill( anEvent->energyR / anEvent->energyL );
      h1_deltaT_raw[label_energyBin] -> Fill( anEvent->timeR-anEvent->timeL );
    }
    std::cout << std::endl;
  }
  
  
  
  //------------------
  //--- draw 2nd plots
  std::map<std::string,float> CTRMeans;
  std::map<std::string,float> CTRSigmas;
  
  std::map<std::string,TF1*> fitFunc_energyRatio;
  
  for(auto mapIt : h1_deltaT_raw)
  {
    FindSmallestInterval(vals,h1_deltaT_raw[mapIt.first],0.68);
    float mean = vals[0];
    float min = vals[4];
    float max = vals[5];
    float delta = max-min;
    float sigma = 0.5*delta;
    float effSigma = sigma;
    CTRMeans[mapIt.first] = mean;
    CTRSigmas[mapIt.first] = effSigma;
  }
  
  
  for(auto stepLabel : stepLabels)
  {
    float Vov = map_Vovs[stepLabel];
    float th = map_ths[stepLabel];
    
    for(int iBar = 0; iBar < 16; ++iBar)
    {
      std::string label(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
      
      if( !ranges[label] ) continue;
      int nEnergyBins = ranges[label]->size()-1;
      
      for(int iEnergyBin = 0; iEnergyBin < nEnergyBins; ++iEnergyBin)
      {
        std::string label_energyBin(Form("%s_energyBin%d",label.c_str(),iEnergyBin));
        
        c = new TCanvas(Form("c_energyRatio_%s",label_energyBin.c_str()),Form("c_energyRatio_%s",label_energyBin.c_str()));
        
        histo = h1_energyRatio[label_energyBin];
        histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
        histo -> SetTitle(Form(";energy_{right} / energy_{left};entries"));
        histo -> SetLineColor(kRed);
        histo -> SetLineWidth(2);
        histo -> Draw();
        histo -> Write();
        
        fitFunc_energyRatio[label_energyBin] = new TF1(Form("fitFunc_energyRatio_%s",label_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
        histo -> Fit(fitFunc_energyRatio[label_energyBin],"QRLS+");
        fitFunc_energyRatio[label_energyBin] -> SetLineColor(kBlack);
        fitFunc_energyRatio[label_energyBin] -> SetLineWidth(2);
        fitFunc_energyRatio[label_energyBin] -> Draw("same");
        
        latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.1f V, th. = %d DAC}",iBar,Vov,int(th)));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kRed);
        latex -> Draw("same");
        
        c -> Print(Form("%s/energyRatio/c_energyRatio__%s.png",plotDir.c_str(),label_energyBin.c_str()));
        c -> Print(Form("%s/energyRatio/c_energyRatio__%s.pdf",plotDir.c_str(),label_energyBin.c_str()));
        delete c;
      }
    }
  }
  
  
  
  
  
  
  //------------------------
  //--- 3rd loop over events
  for(auto mapIt : trees)
  {
    std::string label = mapIt.first;
    
    std::cout << ">>> 3rd loop: " << label << std::endl;
    
    ModuleEventClass* anEvent = new ModuleEventClass();
    mapIt.second -> SetBranchAddress("event",&anEvent);
    
    int nEntries = mapIt.second->GetEntries();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 3rd loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
      mapIt.second -> GetEntry(entry);
      
      if( !accept[label][entry] ) continue;
      
      std::string labelLR(Form("bar%02dL-R_Vov%.1f_th%02d",anEvent->barID,anEvent->Vov,anEvent->vth1));      
      int energyBin = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges[label]);
      std::string label_energyBin(Form("%s_energyBin%d",labelLR.c_str(),energyBin));
      
      float energyRatioMean = fitFunc_energyRatio[label_energyBin]->GetParameter(1);
      float energyRatioSigma = fitFunc_energyRatio[label_energyBin]->GetParameter(2);
      
      if( fabs(anEvent->energyR/anEvent->energyL-energyRatioMean) > 2.*energyRatioSigma )
      {
        accept[label][entry] = false;
        continue;
      }
      
      
      long long deltaT = anEvent->timeR - anEvent->timeL;
      
      if( h1_deltaT[label_energyBin] == NULL )
      {
        h1_deltaT[label_energyBin] = new TH1F(Form("h1_deltaT_%s",label_energyBin.c_str()),"",1000,-5000,5000.);
      }
      
      h1_deltaT[label_energyBin] -> Fill( deltaT );
      
      
      float timeLow = CTRMeans[label_energyBin] - 1.* CTRSigmas[label_energyBin];
      float timeHig = CTRMeans[label_energyBin] + 1.* CTRSigmas[label_energyBin];
      
      if( !p1_deltaT_vs_energyRatio[label_energyBin] )
      {
        p1_deltaT_vs_energyRatio[label_energyBin] = new TProfile(Form("p1_deltaT_vs_energyRatio_%s",label_energyBin.c_str()),"",50,energyRatioMean-2.*energyRatioSigma,energyRatioMean+2.*energyRatioSigma);
      }
      
      if( ( deltaT > timeLow ) && ( deltaT < timeHig ) )
        p1_deltaT_vs_energyRatio[label_energyBin] -> Fill( anEvent->energyR/anEvent->energyL,deltaT );
    }
    std::cout << std::endl;
  }
  
  
  
  //------------------
  //--- draw 3rd plots
  std::map<std::string,TF1*> fitFunc_energyRatioCorr;
  
  for(auto stepLabel : stepLabels)
  {
    float Vov = map_Vovs[stepLabel];
    float th = map_ths[stepLabel];
    
    for(int iBar = 0; iBar < 16; ++iBar)
    {
      std::string label(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
      
      if( !ranges[label] ) continue;
      int nEnergyBins = ranges[label]->size()-1;
      
      for(int iEnergyBin = 0; iEnergyBin < nEnergyBins; ++iEnergyBin)
      {
        std::string label_energyBin(Form("%s_energyBin%d",label.c_str(),iEnergyBin));
        
        c = new TCanvas(Form("c_deltaT_vs_energyRatio_%s",label_energyBin.c_str()),Form("c_deltaT_vs_energyRatio_%s",label_energyBin.c_str()));
        
        prof = p1_deltaT_vs_energyRatio[label_energyBin];
        prof -> SetTitle(Form("energy_{right} / energy_{left};#Deltat [ps]"));
        prof -> Draw("");
        
        latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.1f V, th. = %d DAC}",iBar,Vov,int(th)));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kRed);
        latex -> Draw("same");
        
        fitFunc_energyRatioCorr[label_energyBin] = new TF1(Form("fitFunc_energyRatioCorr_%s",label_energyBin.c_str()),"pol1",0.,5.);
        p1_deltaT_vs_energyRatio[label_energyBin] -> Fit(fitFunc_energyRatioCorr[label_energyBin],"QRS+");
        fitFunc_energyRatioCorr[label_energyBin] -> SetLineColor(kRed);
        fitFunc_energyRatioCorr[label_energyBin] -> SetLineWidth(2);
        fitFunc_energyRatioCorr[label_energyBin] -> Draw("same");
        
        c -> Print(Form("%s/CTR/c_deltaT_vs_energyRatio__%s.png",plotDir.c_str(),label_energyBin.c_str()));
        c -> Print(Form("%s/CTR/c_deltaT_vs_energyRatio__%s.pdf",plotDir.c_str(),label_energyBin.c_str()));
        delete c;
      }
    }
  }
  
  
  
  
  
  
  //------------------------
  //--- 4th loop over events
  for(auto mapIt : trees)
  {
    std::string label = mapIt.first;
    
    std::cout << ">>> 4th loop: " << label << std::endl;
    
    ModuleEventClass* anEvent = new ModuleEventClass();
    mapIt.second -> SetBranchAddress("event",&anEvent);
    
    int nEntries = mapIt.second->GetEntries();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 4th loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
      mapIt.second -> GetEntry(entry);
      
      if( !accept[label][entry] ) continue;
      
      
      int energyBin = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges[label]);
      std::string labelLR(Form("bar%02dL-R_Vov%.1f_th%02d",anEvent->barID,anEvent->Vov,anEvent->vth1));
      std::string label_energyBin(Form("%s_energyBin%d",labelLR.c_str(),energyBin));
      
      long long deltaT = anEvent->timeR - anEvent->timeL;
      
      float energyRatioCorr = 1.;
      if( fitFunc_energyRatioCorr[label_energyBin] )
        energyRatioCorr = fitFunc_energyRatioCorr[label_energyBin]->Eval(anEvent->energyR/anEvent->energyL) -
                          fitFunc_energyRatioCorr[label_energyBin]->Eval(fitFunc_energyRatio[label_energyBin]->GetParameter(1));
      
      if( h1_deltaT_energyRatioCorr[label_energyBin] == NULL )
      {
        h1_deltaT_energyRatioCorr[label_energyBin] = new TH1F(Form("h1_deltaT_energyRatioCorr_%s",label_energyBin.c_str()),"",1000,-5000.,5000.);
      }
      
      h1_deltaT_energyRatioCorr[label_energyBin] -> Fill( deltaT - energyRatioCorr );
    }
    std::cout << std::endl;
  }
  
  
  
  //------------------
  //--- draw 4th plots
  std::map<std::string,TGraphErrors*> g_tRes_gaus_vs_th;
  std::map<std::string,TGraphErrors*> g_tRes_gaus_bestTh_vs_Vov;
  
  std::map<std::string,TGraphErrors*> g_tRes_energyRatioCorr_gaus_vs_th;
  std::map<std::string,TGraphErrors*> g_tRes_energyRatioCorr_gaus_bestTh_vs_Vov;
  
  std::map<std::string,std::map<float,float> > tRes_gaus_bestTh;
  std::map<std::string,std::map<float,float> > tResErr_gaus_bestTh;
  std::map<std::string,std::map<float,float> > tRes_energyRatioCorr_gaus_bestTh;
  std::map<std::string,std::map<float,float> > tResErr_energyRatioCorr_gaus_bestTh;
  
  for(auto stepLabel : stepLabels)
  {
    float Vov = map_Vovs[stepLabel];
    float th = map_ths[stepLabel];
    
    for(int iBar = 0; iBar < 16; ++iBar)
    {
      std::string label(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
      
      if( !ranges[label] ) continue;
      int nEnergyBins = ranges[label]->size()-1;
      
      for(int iEnergyBin = 0; iEnergyBin < nEnergyBins; ++iEnergyBin)
      {
        std::string label_energyBin(Form("%s_energyBin%d",label.c_str(),iEnergyBin));
        
        
        c = new TCanvas(Form("c_deltaT_energyRatioCorr_%s",label_energyBin.c_str()),Form("c_deltaT_energyRatioCorr_%s",label_energyBin.c_str()));
        
	// -- energy corr deltaT
        histo = h1_deltaT_energyRatioCorr[label_energyBin];
        histo -> SetTitle(Form(";energy-corrected #Deltat [ps];entries"));
        histo -> SetLineWidth(2);
        histo -> SetLineColor(kBlue);
        histo -> SetMarkerColor(kBlue);
        histo -> Rebin(2);
        histo -> Draw("");
        
        float fitXMin = CTRMeans[label_energyBin] - 3.*CTRSigmas[label_energyBin];
        float fitXMax = CTRMeans[label_energyBin] + 3.*CTRSigmas[label_energyBin];
        TF1* fitFunc = new TF1(Form("fitFunc_energyRatioCorr_%s",label_energyBin.c_str()),"gaus",fitXMin,fitXMax);
        fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
        histo -> Fit(fitFunc,"QNRSL+","");
        histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
        histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-2.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.*fitFunc->GetParameter(2));
        histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
        
        float fitXMin2 = fitFunc->GetParameter(1)-fabs(2.0*fitFunc->GetParameter(2));
        float fitXMax2 = fitFunc->GetParameter(1)+fabs(2.0*fitFunc->GetParameter(2));
        TF1* fitFunc2 = new TF1(Form("fitFunc2_energyRatioCorr_%s",label_energyBin.c_str()),"gaus",fitXMin2,fitXMax2);
        fitFunc2 -> SetParameter(0,fitFunc->GetParameter(0));
        fitFunc2 -> SetParameter(1,fitFunc->GetParameter(1));
        fitFunc2 -> SetParameter(2,fitFunc->GetParameter(2));
        histo -> Fit(fitFunc2,"QRSL+");
        
        fitFunc2 -> SetLineColor(kBlue+1);
        fitFunc2 -> SetLineWidth(3);
        fitFunc2 -> Draw("same");
        
        FindSmallestInterval(vals,histo,0.68);
        float mean = vals[0];
        float min = vals[4];
        float max = vals[5];
        float delta = max-min;
        float sigma = 0.5*delta;
        float effSigma = sigma;
        
        latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corr.}^{eff} = %.0f ps}{#sigma_{corr.}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc2->GetParameter(2))));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kBlue);
        latex -> Draw("same");
        
        std::string label_vs_th(Form("bar%02dL-R_Vov%.1f_energyBin%d",iBar,Vov,iEnergyBin));
        if( g_tRes_gaus_vs_th[label_vs_th] == NULL )
        {
          g_tRes_gaus_vs_th[label_vs_th] = new TGraphErrors();
          g_tRes_energyRatioCorr_gaus_vs_th[label_vs_th] = new TGraphErrors();
        }
        
        std::string label_bestTh_vs_Vov(Form("bar%02dL-R_bestTh_energyBin%d",iBar,iEnergyBin));
        if( g_tRes_energyRatioCorr_gaus_bestTh_vs_Vov[label_bestTh_vs_Vov] == NULL )
        {
          g_tRes_gaus_bestTh_vs_Vov[label_bestTh_vs_Vov] = new TGraphErrors();
          g_tRes_energyRatioCorr_gaus_bestTh_vs_Vov[label_bestTh_vs_Vov] = new TGraphErrors();
        } 
        
        g_tRes_energyRatioCorr_gaus_vs_th[label_vs_th] -> SetPoint(g_tRes_energyRatioCorr_gaus_vs_th[label_vs_th]->GetN(),th,fabs(fitFunc2->GetParameter(2))*corr);
        g_tRes_energyRatioCorr_gaus_vs_th[label_vs_th] -> SetPointError(g_tRes_energyRatioCorr_gaus_vs_th[label_vs_th]->GetN()-1,0.,fitFunc2->GetParError(2)*corr);
        
        if( ( tRes_energyRatioCorr_gaus_bestTh[label_bestTh_vs_Vov][Vov] != 0. && fabs(fitFunc2->GetParameter(2))*corr < tRes_energyRatioCorr_gaus_bestTh[label_bestTh_vs_Vov][Vov] ) ||
            tRes_energyRatioCorr_gaus_bestTh[label_bestTh_vs_Vov][Vov] == 0 )
        {
          tRes_energyRatioCorr_gaus_bestTh[label_bestTh_vs_Vov][Vov] = fabs(fitFunc2->GetParameter(2))*corr;
          tResErr_energyRatioCorr_gaus_bestTh[label_bestTh_vs_Vov][Vov] = fitFunc2->GetParError(2)*corr;
        }
        
        
        
	// -- raw delta T
	c->cd();
        histo = h1_deltaT[label_energyBin];
        histo -> SetLineWidth(2);
        histo -> SetLineColor(kRed);
        histo -> SetMarkerColor(kRed);
        histo -> Rebin(2);
        histo -> Draw("same");
        
        fitFunc = new TF1(Form("fitFunc_%s",label_energyBin.c_str()),"gaus",fitXMin,fitXMax);
        fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
        histo -> Fit(fitFunc,"QNRSL+","");
        histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
        histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-2.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.*fitFunc->GetParameter(2));
        histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
        
        fitXMin2 = fitFunc->GetParameter(1)-fabs(2.*fitFunc->GetParameter(2));
        fitXMax2 = fitFunc->GetParameter(1)+fabs(2.*fitFunc->GetParameter(2));
        fitFunc2 = new TF1(Form("fitFunc2_energyRatioCorr_%s",label_energyBin.c_str()),"gaus",fitXMin2,fitXMax2);
        fitFunc2 -> SetParameter(0,fitFunc->GetParameter(0));
        fitFunc2 -> SetParameter(1,fitFunc->GetParameter(1));
        fitFunc2 -> SetParameter(2,fitFunc->GetParameter(2));
        histo -> Fit(fitFunc2,"QRSL+");
        
        fitFunc2 -> SetLineColor(kRed+1);
        fitFunc2 -> SetLineWidth(3);
        fitFunc2 -> Draw("same");
        
        FindSmallestInterval(vals,histo,0.68);
        mean = vals[0];
        min = vals[4];
        max = vals[5];
        delta = max-min;
        sigma = 0.5*delta;
        effSigma = sigma;
        
        histo -> GetXaxis() -> SetRangeUser(mean-5.*sigma,mean+5.*sigma);
        
        latex = new TLatex(0.20,0.85,Form("#splitline{#sigma_{raw}^{eff} = %.0f ps}{#sigma_{raw}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc2->GetParameter(2))));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kRed);
        latex -> Draw("same");
        
        g_tRes_gaus_vs_th[label_vs_th] -> SetPoint(g_tRes_gaus_vs_th[label_vs_th]->GetN(),th,fabs(fitFunc2->GetParameter(2))*corr);
        g_tRes_gaus_vs_th[label_vs_th] -> SetPointError(g_tRes_gaus_vs_th[label_vs_th]->GetN()-1,0.,fitFunc2->GetParError(2)*corr);
        
        if( ( tRes_gaus_bestTh[label_bestTh_vs_Vov][Vov] != 0. && fabs(fitFunc2->GetParameter(2))*corr < tRes_gaus_bestTh[label_bestTh_vs_Vov][Vov] ) ||
            tRes_gaus_bestTh[label_bestTh_vs_Vov][Vov] == 0. )
        {
          tRes_gaus_bestTh[label_bestTh_vs_Vov][Vov] = fabs(fitFunc2->GetParameter(2))*corr;
          tResErr_gaus_bestTh[label_bestTh_vs_Vov][Vov] = fitFunc2->GetParError(2)*corr;
        }
        
        c -> Print(Form("%s/CTR_energyRatioCorr/c_deltaT_energyRatioCorr__%s.png",plotDir.c_str(),label_energyBin.c_str()));
        c -> Print(Form("%s/CTR_energyRatioCorr/c_deltaT_energyRatioCorr__%s.pdf",plotDir.c_str(),label_energyBin.c_str()));
        delete c;
      }
    }
  }
  
  
  
  for(int iBar = 0; iBar < 16; ++iBar)
  {
    if( !ranges[Form("bar%02dL-R_Vov%.1f_th%02.0f",iBar,refVov,refTh)] ) continue;
    int nEnergyBins = ranges[Form("bar%02dL-R_Vov%.1f_th%02.0f",iBar,refVov,refTh)]->size()-1;
    
    for(int iEnergyBin = 0; iEnergyBin < nEnergyBins; ++iEnergyBin)
    {
      c = new TCanvas(Form("c_tRes_vs_th"),Form("c_tRes_vs_th"));
      
      TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,tResMin,64.,tResMax) );
      if( tResMode == 0 ) hPad -> SetTitle(";threshold [DAC];#sigma_{#Deltat} [ps]");
      if( tResMode == 1 ) hPad -> SetTitle(";threshold [DAC];#sigma_{#Deltat} / #sqrt{2} [ps]");
      if( tResMode == 2 ) hPad -> SetTitle(";threshold [DAC];#sigma_{#Deltat} / 2 [ps]");
      hPad -> Draw();
      gPad -> SetGridy();
      
      int iter = 0;
      for(auto mapIt : VovLabels)
      {
        std::string label(Form("bar%02dL-R_%s_energyBin%d",iBar,mapIt.first.c_str(),iEnergyBin));
        TGraph* g_gaus = g_tRes_gaus_vs_th[label];
        TGraph* g_energyRatioCorr_gaus = g_tRes_energyRatioCorr_gaus_vs_th[label];
        
        g_gaus -> SetLineColor(1+iter);
        g_gaus -> SetMarkerColor(1+iter);
        g_gaus -> SetMarkerStyle(25);
        g_gaus -> Draw("PL,same");
        
        g_energyRatioCorr_gaus -> SetLineColor(1+iter);
        g_energyRatioCorr_gaus -> SetMarkerColor(1+iter);
        g_energyRatioCorr_gaus -> SetMarkerStyle(20);
        g_energyRatioCorr_gaus -> Draw("PL,same");
        
        latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kBlack+iter);
        latex -> Draw("same");
        
        ++iter;
      }
      
      c -> Print(Form("%s/c_tRes_vs_th__bar%02d_energyBin%d.png",plotDir.c_str(),iBar,iEnergyBin));
      c -> Print(Form("%s/c_tRes_vs_th__bar%02d_energyBin%d.pdf",plotDir.c_str(),iBar,iEnergyBin));
      delete c;
    }
  }
  
  
  
  for(int iBar = 0; iBar < 16; ++iBar)
  {
    if( !ranges[Form("bar%02dL-R_Vov%.1f_th%02.0f",iBar,refVov,refTh)] ) continue;
    int nEnergyBins = ranges[Form("bar%02dL-R_Vov%.1f_th%02.0f",iBar,refVov,refTh)]->size()-1;
    
    for(int iEnergyBin = 0; iEnergyBin < nEnergyBins; ++iEnergyBin)
    {
      c = new TCanvas(Form("c_tRes_bestTh_vs_Vov"),Form("c_tRes_bestTh_vs_Vov"));
      
      legend2 = new TLegend(0.20,0.82,0.50,0.90);
      legend2 -> SetFillColor(0);
      legend2 -> SetFillStyle(1000);
      legend2 -> SetTextFont(42);
      
      TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,tResMin,10.,tResMax) );
      if( tResMode == 0 ) hPad -> SetTitle(";V_{OV} [V];#sigma_{#Deltat} [ps]");
      if( tResMode == 1 ) hPad -> SetTitle(";V_{OV} [V];#sigma_{#Deltat} / #sqrt{2} [ps]");
      if( tResMode == 2 ) hPad -> SetTitle(";V_{OV} [V];#sigma_{#Deltat} / 2 [ps]");
      hPad -> Draw();
      gPad -> SetGridy();
      
      std::string label(Form("bar%02dL-R_bestTh_energyBin%d",iBar,iEnergyBin));
      TGraph* g_gaus = g_tRes_gaus_bestTh_vs_Vov[label];
      TGraph* g_energyRatioCorr_gaus = g_tRes_energyRatioCorr_gaus_bestTh_vs_Vov[label];
      
      g_gaus -> SetLineColor(1);
      g_gaus -> SetMarkerColor(1);
      g_gaus -> SetMarkerStyle(25);
      g_gaus -> Draw("PL,same");
      
      g_energyRatioCorr_gaus -> SetLineColor(1);
      g_energyRatioCorr_gaus -> SetMarkerColor(1);
      g_energyRatioCorr_gaus -> SetMarkerStyle(20);
      g_energyRatioCorr_gaus -> Draw("PL,same");
      
      latex = new TLatex(0.55,0.85,Form("best th."));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack);
      latex -> Draw("same");
      
      legend2->AddEntry(g_gaus,"raw","PL");
      legend2->AddEntry(g_energyRatioCorr_gaus,"corrected","PL");
      legend2->Draw("same");

      c -> Print(Form("%s/c_tRes_bestTh_vs_Vov__bar%02d_l-R_energyBin%d.png",plotDir.c_str(),iBar,iEnergyBin));
      
      delete c;
    }

  }
  
  
  
  
  
  
  // outFile -> cd();
  
  // for(auto mapIt: g_energy_vs_th)  mapIt.second -> Write(Form("g_energy_vs_th_%s", mapIt.first.c_str()));
  // for(auto mapIt: g_energy_vs_Vov) mapIt.second -> Write(Form("g_energy_vs_Vov_%s",mapIt.first.c_str()));
  
  /*
  for(auto mapIt: g_tRes_gaus_vs_th)             mapIt.second -> Write(Form("g_tRes_gaus_vs_th_%s",             mapIt.first.c_str()));
  for(auto mapIt: g_tRes_gaus_vs_Vov)            mapIt.second -> Write(Form("g_tRes_gaus_vs_Vov_%s",            mapIt.first.c_str()));
  for(auto mapIt: g_tRes_gaus_bestTh_vs_Vov)     mapIt.second -> Write(Form("g_tRes_gaus_bestTh_vs_Vov_%s",     mapIt.first.c_str()));
  
  for(auto mapIt: g_tRes_energyRatioCorr_gaus_vs_th)             mapIt.second -> Write(Form("g_tRes_energyRatioCorr_gaus_vs_th_%s",             mapIt.first.c_str()));
  for(auto mapIt: g_tRes_energyRatioCorr_gaus_vs_Vov)            mapIt.second -> Write(Form("g_tRes_energyRatioCorr_gaus_vs_Vov_%s",            mapIt.first.c_str()));
  for(auto mapIt: g_tRes_energyRatioCorr_gaus_bestTh_vs_Vov)     mapIt.second -> Write(Form("g_tRes_energyRatioCorr_gaus_bestTh_vs_Vov_%s",     mapIt.first.c_str()));
  */
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
