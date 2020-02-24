#include "interface/AnalysisUtils.h"
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
#include "TError.h"





int main(int argc, char** argv)
{
  setTDRStyle();
  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(0);
  gErrorIgnoreLevel = kWarning; // to suppress Info messages
  
  if( argc < 2 )
  {
    std::cout << ">>> drawTOFPET2::usage:   " << argv[0] << " configFile.cfg" << std::endl;
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
  system(Form("mkdir -p %s/totRatio/",plotDir.c_str()));
  system(Form("mkdir -p %s/energy/",plotDir.c_str()));
  system(Form("mkdir -p %s/energyRatio/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR_energyCorr/",plotDir.c_str()));
  
  
  //--- open files
  std::string plotFileName= opts.GetOpt<std::string>("Input.plotFileName");
  TFile* inFile = TFile::Open(plotFileName.c_str(),"READ");
 
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
    std::size_t found = name.find("h1_energy_cut");
    if( found!=std::string::npos )
    {
      std::stringstream ss(name);
      std::string token;
      std::vector<std::string> tokens;
      while( std::getline(ss, token, '_') )
      {
        tokens.push_back(token);
      }
      std::string stepLabel = tokens[4]+"_"+tokens[5];
      VovLabels[tokens[4]] += 1;
      thLabels[tokens[5]] += 1;
      stepLabels.push_back(stepLabel);
      
      std::string string_Vov = tokens[4];
      string_Vov.erase(0,3);
      map_Vovs[stepLabel] = atof(string_Vov.c_str());
      
      std::string string_th = tokens[5];
      string_th.erase(0,2);
      map_ths[stepLabel] = atof(string_th.c_str());
    }
  }
  std::sort(stepLabels.begin(),stepLabels.end());
  stepLabels.erase(std::unique(stepLabels.begin(),stepLabels.end()),stepLabels.end());
  for(auto label : stepLabels)
    std::cout << label << std::endl;
  
  
  
  //--- define channels
  std::vector<std::string> channels = opts.GetOpt<std::vector<std::string> >("Channels.channels");
  
  std::vector<std::string> pairs = opts.GetOpt<std::vector<std::string> >("Channels.pairs");
  std::vector<std::pair<std::string,std::string> > pairsVec;
  for(unsigned int ii = 0; ii < pairs.size()/2; ++ii)
    pairsVec.push_back(std::make_pair(pairs.at(0+ii*2),pairs.at(1+ii*2)));
  std::vector<int> pairsMode = opts.GetOpt<std::vector<int> >("Channels.pairsMode");
  
  
  //--- get cuts per bar / Vov
  std::map<unsigned int,std::map<float,float> > cut_qfineAcc;
  std::map<unsigned int,std::map<float,float> > cut_totAcc;
  std::map<unsigned int,std::map<float,float> > cut_energyAcc;
  std::map<unsigned int,std::map<float,float> > cut_energyFitMin;
  std::map<unsigned int,std::map<float,float> > cut_energyFitMax;
  for(auto ch :  channels)
  {
    int chID = opts.GetOpt<int>(Form("%s.chID",ch.c_str()));
    std::vector<float> Vovs          = opts.GetOpt<std::vector<float> >(Form("%s.Vovs",ch.c_str()));
    std::vector<float> qfineMins     = opts.GetOpt<std::vector<float> >(Form("%s.qfineMins",ch.c_str()));
    std::vector<float> totMins       = opts.GetOpt<std::vector<float> >(Form("%s.totMins",ch.c_str()));
    std::vector<float> energyMins    = opts.GetOpt<std::vector<float> >(Form("%s.energyMins",ch.c_str()));
    std::vector<float> energyFitMins = opts.GetOpt<std::vector<float> >(Form("%s.energyFitMins",ch.c_str()));
    std::vector<float> energyFitMaxs = opts.GetOpt<std::vector<float> >(Form("%s.energyFitMaxs",ch.c_str()));
    int iter = 0;
    for(auto Vov : Vovs)
    {
      cut_qfineAcc[chID][Vov]     = qfineMins.at(iter);
      cut_totAcc[chID][Vov]       = totMins.at(iter);
      cut_energyAcc[chID][Vov]    = energyMins.at(iter);
      cut_energyFitMin[chID][Vov] = energyFitMins.at(iter);
      cut_energyFitMax[chID][Vov] = energyFitMaxs.at(iter);
      ++iter;
    }
  }
  std::map<std::string,float> cut_energyMin;
  std::map<std::string,float> cut_energyMax;
  
  
  //--- get plot settings
  std::vector<int> plots = opts.GetOpt<std::vector<int> >("Plots.plots");
  float tResMin = opts.GetOpt<float>("Plots.tResMin");
  float tResMax = opts.GetOpt<float>("Plots.tResMax");
  
  
  
  TCanvas* c;
  TCanvas *cc;
  TCanvas *ccc;
  float* vals = new float[6];
  TLatex* latex;
  TH1F* histo;
  TH2F* histo2;
  TProfile* prof;
  TLegend* legend;
  TLegend* legend2;
  

  TH1F* hamp;
  TH1F* hampL;
  TH1F* hampR;
  TF1 *fLandau1;
  TF1 *fLandau2;
  float norm;

  
  //------------------
  //--- draw 1st plots
  if( std::find(plots.begin(),plots.end(),1) != plots.end() )
  {
    std::map<std::string,TGraphErrors*> g_tot_vs_th;
    std::map<std::string,TGraphErrors*> g_tot_vs_Vov;
    std::map<std::string,TGraphErrors*> g_energy_vs_th;
    std::map<std::string,TGraphErrors*> g_energy_vs_Vov;
    
    for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float th = map_ths[stepLabel];
      std::string VovLabel(Form("Vov%.1f",Vov));
      std::string thLabel(Form("th%02.0f",th));
      
      
      //--------------------------------------------------------
      
      
      for(auto ch : channels)
      {
        int chID = opts.GetOpt<int>(Form("%s.chID",ch.c_str()));  
        std::string label(Form("%s_%s",ch.c_str(),stepLabel.c_str()));
        
        
        c = new TCanvas(Form("c_qfine_%s",label.c_str()),Form("c_qfine_%s",label.c_str()));
        gPad -> SetLogy();
        
        histo = (TH1F*)( inFile->Get(Form("h1_qfine_%s",label.c_str())) );
        histo -> SetTitle(";Q_{fine} [ADC];entries");
        histo -> SetLineColor(kRed);
        histo -> Draw();
        histo -> GetXaxis() -> SetRangeUser(12,513);
        TLine* line_qfineAcc1 = new TLine(cut_qfineAcc[chID][Vov],histo->GetMinimum(),cut_qfineAcc[chID][Vov],histo->GetMaximum());
        line_qfineAcc1 -> SetLineColor(kBlack);
        line_qfineAcc1 -> Draw("same");
        latex = new TLatex(0.65,0.85,Form("%s",ch.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kRed);
        latex -> Draw("same");      
        c -> Print(Form("%s/qfine/c_qfine__%s.png",plotDir.c_str(),label.c_str()));
        c -> Print(Form("%s/qfine/c_qfine__%s.pdf",plotDir.c_str(),label.c_str()));
        delete c;
        
        
        histo = (TH1F*)( inFile->Get(Form("h1_qfine_%s",label.c_str())) );      
        c = new TCanvas(Form("c_qfine_vs_tot_%s",label.c_str()),Form("c_qfine_vs_tot_%s",label.c_str()));
        gPad -> SetLogz();
        
        
        histo2 = (TH2F*)( inFile->Get(Form("h2_qfine_vs_tot_%s",label.c_str())) );            
        histo2 -> SetTitle(Form(";%s ToT [ns];%s Q_{fine} [ADC]",ch.c_str(),ch.c_str()));
        histo2 -> Draw("colz");
        
        c -> Print(Form("%s/qfine/c_qfine_vs_tot__%s.png",plotDir.c_str(),label.c_str()));
        c -> Print(Form("%s/qfine/c_qfine_vs_tot__%s.pdf",plotDir.c_str(),label.c_str()));
        delete c;
        
        
        c = new TCanvas(Form("c_tot_%s",label.c_str()),Form("c_tot_%s",label.c_str()));
        // gPad -> SetLogy();
        
        histo = (TH1F*)( inFile->Get(Form("h1_tot_%s",label.c_str())) );
        histo -> SetTitle(";ToT [ns];entries");
        histo -> SetLineColor(kRed);
        histo -> Draw();
        float max1 = FindXMaximum(histo,cut_totAcc[chID][Vov],1000.);
        histo -> GetXaxis() -> SetRangeUser(0.25*max1,2.*max1);
        TF1* fitFunc1 = new TF1("fitFunc1","gaus",max1-0.05*max1,max1+0.05*max1);
        histo -> Fit(fitFunc1,"QNRS+");
        fitFunc1 -> SetLineColor(kBlack);
        fitFunc1 -> SetLineWidth(3);
        fitFunc1 -> Draw("same");
        TLine* line_totAcc1 = new TLine(cut_totAcc[chID][Vov],histo->GetMinimum(),cut_totAcc[chID][Vov],histo->GetMaximum());
        line_totAcc1 -> SetLineColor(kBlack);
        line_totAcc1 -> Draw("same");
        latex = new TLatex(0.65,0.85,Form("%s",ch.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kRed);
        latex -> Draw("same");      
        c -> Print(Form("%s/tot/c_tot__%s.png",plotDir.c_str(),label.c_str()));
        c -> Print(Form("%s/tot/c_tot__%s.pdf",plotDir.c_str(),label.c_str()));
        delete c;
        
        
        if( g_tot_vs_th[Form("%s_%s",ch.c_str(),VovLabel.c_str())] == NULL )
          g_tot_vs_th[Form("%s_%s",ch.c_str(),VovLabel.c_str())] = new TGraphErrors();
        
        if( g_tot_vs_Vov[Form("%s_%s",ch.c_str(),thLabel.c_str())] == NULL )
          g_tot_vs_Vov[Form("%s_%s",ch.c_str(),thLabel.c_str())] = new TGraphErrors();
        
        g_tot_vs_th[Form("%s_%s",ch.c_str(),VovLabel.c_str())] -> SetPoint(g_tot_vs_th[Form("%s_%s",ch.c_str(),VovLabel.c_str())]->GetN(),th,fitFunc1->GetMaximumX());
        g_tot_vs_th[Form("%s_%s",ch.c_str(),VovLabel.c_str())] -> SetPointError(g_tot_vs_th[Form("%s_%s",ch.c_str(),VovLabel.c_str())]->GetN()-1,0.,0.);
        
        g_tot_vs_Vov[Form("%s_%s",ch.c_str(),thLabel.c_str())] -> SetPoint(g_tot_vs_Vov[Form("%s_%s",ch.c_str(),thLabel.c_str())]->GetN(),Vov,fitFunc1->GetMaximumX());
        g_tot_vs_Vov[Form("%s_%s",ch.c_str(),thLabel.c_str())] -> SetPointError(g_tot_vs_Vov[Form("%s_%s",ch.c_str(),thLabel.c_str())]->GetN()-1,0.,0.);
        
        
        c = new TCanvas(Form("c_energy_%s",label.c_str()),Form("c_energy_%s",label.c_str()));
        // gPad -> SetLogy();
        
        histo = (TH1F*)( inFile->Get(Form("h1_energy_%s",label.c_str())) );      
        histo -> SetTitle(";energy [a.u.];entries");
        histo -> SetLineColor(kRed);
        histo -> Draw();
        std::cout << label << std::endl;
        max1 = FindXMaximum(histo,cut_energyAcc[chID][Vov],100.);
        // histo -> GetXaxis() -> SetRangeUser(0.1*max1,2.5*max1);
        histo -> GetXaxis() -> SetRangeUser(0.,50.);
        fitFunc1 = new TF1("fitFunc1","gaus",max1-cut_energyFitMin[chID][Vov]*max1,max1+cut_energyFitMax[chID][Vov]*max1);
        histo -> Fit(fitFunc1,"QNRS+");
        fitFunc1 -> SetLineColor(kBlack);
        fitFunc1 -> SetLineWidth(3);
        fitFunc1 -> Draw("same");
        // cut_energyMin[Form("%s_%s",ch.c_str(),stepLabel.c_str())] = fitFunc1->GetMaximumX()-cut_energyFitMin[chID][Vov]*fitFunc1->GetMaximumX();
        // cut_energyMax[Form("%s_%s",ch.c_str(),stepLabel.c_str())] = fitFunc1->GetMaximumX()+cut_energyFitMax[chID][Vov]*fitFunc1->GetMaximumX();
        cut_energyMin[Form("%s_%s",ch.c_str(),stepLabel.c_str())] = cut_energyAcc[chID][Vov];
        cut_energyMax[Form("%s_%s",ch.c_str(),stepLabel.c_str())] = 50.;
        
        TLine* line_energyMin1 = new TLine(cut_energyMin[Form("%s_%s",ch.c_str(),stepLabel.c_str())],histo->GetMinimum(),
                                           cut_energyMin[Form("%s_%s",ch.c_str(),stepLabel.c_str())],histo->GetMaximum());
        line_energyMin1 -> SetLineColor(kBlack);
        line_energyMin1 -> Draw("same");
        TLine* line_energyMax1 = new TLine(cut_energyMax[Form("%s_%s",ch.c_str(),stepLabel.c_str())],histo->GetMinimum(),
                                           cut_energyMax[Form("%s_%s",ch.c_str(),stepLabel.c_str())],histo->GetMaximum());
        line_energyMax1 -> SetLineColor(kBlack);
        line_energyMax1 -> Draw("same");
        latex = new TLatex(0.65,0.85,Form("%s",ch.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kRed);
        latex -> Draw("same");      
        c -> Print(Form("%s/energy/c_energy__%s.png",plotDir.c_str(),label.c_str()));
        c -> Print(Form("%s/energy/c_energy__%s.pdf",plotDir.c_str(),label.c_str()));
        delete c;
      
        
        if( g_energy_vs_th[Form("%s_%s",ch.c_str(),VovLabel.c_str())] == NULL )
          g_energy_vs_th[Form("%s_%s",ch.c_str(),VovLabel.c_str())] = new TGraphErrors();
        
        if( g_energy_vs_Vov[Form("%s_%s",ch.c_str(),thLabel.c_str())] == NULL )
          g_energy_vs_Vov[Form("%s_%s",ch.c_str(),thLabel.c_str())] = new TGraphErrors();
        
        g_energy_vs_th[Form("%s_%s",ch.c_str(),VovLabel.c_str())] -> SetPoint(g_energy_vs_th[Form("%s_%s",ch.c_str(),VovLabel.c_str())]->GetN(),th,fitFunc1->GetMaximumX());
        g_energy_vs_th[Form("%s_%s",ch.c_str(),VovLabel.c_str())] -> SetPointError(g_energy_vs_th[Form("%s_%s",ch.c_str(),VovLabel.c_str())]->GetN()-1,0.,0.);
        
        g_energy_vs_Vov[Form("%s_%s",ch.c_str(),thLabel.c_str())] -> SetPoint(g_energy_vs_Vov[Form("%s_%s",ch.c_str(),thLabel.c_str())]->GetN(),Vov,fitFunc1->GetMaximumX());
        g_energy_vs_Vov[Form("%s_%s",ch.c_str(),thLabel.c_str())] -> SetPointError(g_energy_vs_Vov[Form("%s_%s",ch.c_str(),thLabel.c_str())]->GetN()-1,0.,0.);
      }
      
      
      //--------------------------------------------------------
      
      
      // for(auto pair : pairsVec)
      // {  
      //   std::string ch1 = pair.first;
      //   std::string ch2 = pair.second;
      //   std::string label12 = Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),stepLabel.c_str());
      
      //   c = new TCanvas(Form("c_tot_corr_%s",label12.c_str()),Form("c_tot_corr_%s",label12.c_str()));
      //   gPad -> SetLogz();
      
      //   histo2 = (TH2F*)( inFile->Get(Form("h2_tot_corr_%s",label12.c_str())) );
      //   histo2 -> SetTitle(Form(";%s ToT [ns];%s ToT [ns]",ch1.c_str(),ch2.c_str()));
      //   histo2 -> Draw("colz");
      
      //   c -> Print(Form("%s/tot/c_tot_corr__%s.png",plotDir.c_str(),label12.c_str()));
      //   c -> Print(Form("%s/tot/c_tot_corr__%s.pdf",plotDir.c_str(),label12.c_str()));
      //   delete c;
      
      
      //   c = new TCanvas(Form("c_energy_corr_%s",label12.c_str()),Form("c_energy_corr_%s",label12.c_str()));
      //   gPad -> SetLogz();
      
      //   histo2 = (TH2F*)( inFile->Get(Form("h2_energy_corr_%s",label12.c_str())) );
      //   histo2 -> SetTitle(Form(";%s energy [a.u.];%s energy [a.u.]",ch1.c_str(),ch2.c_str()));
      //   histo2 -> Draw("colz");
      
      //   c -> Print(Form("%s/energy/c_energy_corr__%s.png",plotDir.c_str(),label12.c_str()));
      //   c -> Print(Form("%s/energy/c_energy_corr__%s.pdf",plotDir.c_str(),label12.c_str()));
      //   delete c;
      // }
      
      
    }  
    
    //--------------------------------------------------------
    
    
    
    
    for(auto pair : pairsVec)
    {
      std::string ch1 = pair.first;
      std::string ch2 = pair.second;
      int isBar1 = opts.GetOpt<int>(Form("%s.isBar",ch1.c_str()));
      int isBar2 = opts.GetOpt<int>(Form("%s.isBar",ch2.c_str()));
      if( isBar1 || isBar2 ) continue;
      
      
      legend = new TLegend(0.20,0.82,0.50,0.90);
      legend -> SetFillColor(0);
      legend -> SetFillStyle(1000);
      legend -> SetTextFont(42);
      legend -> Draw("same");
      
      
      c = new TCanvas(Form("c_tot_vs_th_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_tot_vs_th_%s-%s",ch1.c_str(),ch2.c_str()));
      // gPad -> SetLogy();
      
      TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,500.) );
      hPad -> SetTitle(";threshold [DAC];ToT [ns]");
      hPad -> Draw();
      gPad -> SetGridy();
      
      int iter = 0;
      for(auto mapIt : VovLabels) 
      {
        std::string label1(Form("%s_%s",ch1.c_str(),mapIt.first.c_str()));
        std::string label2(Form("%s_%s",ch2.c_str(),mapIt.first.c_str()));
        TGraph* g_tot1 = g_tot_vs_th[label1];
        TGraph* g_tot2 = g_tot_vs_th[label2];
        
        g_tot1 -> SetLineColor(1+iter);
        g_tot1 -> SetMarkerColor(1+iter);
        g_tot1 -> SetMarkerStyle(20);
        g_tot1 -> Draw("PL,same");
        
        g_tot2 -> SetLineColor(1+iter);
        g_tot2 -> SetMarkerColor(1+iter);
        g_tot2 -> SetMarkerStyle(25);
        g_tot2 -> Draw("PL,same");
        
        latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kBlack+iter);
        latex -> Draw("same");
        
        if( iter == 0 )
        {
          legend -> AddEntry(g_tot1,ch1.c_str(),"PL");
          legend -> AddEntry(g_tot2,ch2.c_str(),"PL");
        }
        ++iter;
      }
      
      legend -> Draw("same");
      
      c -> Print(Form("%s/c_tot_vs_th__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c -> Print(Form("%s/c_tot_vs_th__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      delete c;
      
      
      c = new TCanvas(Form("c_tot_vs_Vov_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_tot_vs_Vov_%s-%s",ch1.c_str(),ch2.c_str()));
      // gPad -> SetLogy();
      
      hPad = (TH1F*)( gPad->DrawFrame(0.,0.,10.,500.) );
      hPad -> SetTitle(";V_{ov} [V];ToT [ns]");
      hPad -> Draw();
      gPad -> SetGridy();
      
      iter = 0;
      for(auto mapIt : thLabels)
      {
        std::string label1(Form("%s_%s",ch1.c_str(),mapIt.first.c_str()));
        std::string label2(Form("%s_%s",ch2.c_str(),mapIt.first.c_str()));
        TGraph* g_tot1 = g_tot_vs_Vov[label1];
        TGraph* g_tot2 = g_tot_vs_Vov[label2];
        
        g_tot1 -> SetLineColor(1+iter);
        g_tot1 -> SetMarkerColor(1+iter);
        g_tot1 -> SetMarkerStyle(20);
        g_tot1 -> Draw("PL,same");
        
        g_tot2 -> SetLineColor(1+iter);
        g_tot2 -> SetMarkerColor(1+iter);
        g_tot2 -> SetMarkerStyle(25);
        g_tot2 -> Draw("PL,same");
        
        latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kBlack+iter);
        latex -> Draw("same");
        
        ++iter;
      }
      
      legend -> Draw("same");
      
      c -> Print(Form("%s/c_tot_vs_Vov__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c -> Print(Form("%s/c_tot_vs_Vov__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      delete c;
      
      
      c = new TCanvas(Form("c_energy_vs_th_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_energy_vs_th_%s-%s",ch1.c_str(),ch2.c_str()));
      // gPad -> SetLogy();
      
      hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,50.) );
      hPad -> SetTitle(";threshold [DAC];energy [a.u.]");
      hPad -> Draw();
      gPad -> SetGridy();
      
      iter = 0;
      for(auto mapIt : VovLabels)
      {
        std::string label1(Form("%s_%s",ch1.c_str(),mapIt.first.c_str()));
        std::string label2(Form("%s_%s",ch2.c_str(),mapIt.first.c_str()));
        TGraph* g_energy1 = g_energy_vs_th[label1];
        TGraph* g_energy2 = g_energy_vs_th[label2];
        
        g_energy1 -> SetLineColor(1+iter);
        g_energy1 -> SetMarkerColor(1+iter);
        g_energy1 -> SetMarkerStyle(20);
        g_energy1 -> Draw("PL,same");
        
        g_energy2 -> SetLineColor(1+iter);
        g_energy2 -> SetMarkerColor(1+iter);
        g_energy2 -> SetMarkerStyle(25);
        g_energy2 -> Draw("PL,same");
        
        latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kBlack+iter);
        latex -> Draw("same");
        
        ++iter;
      }
      
      legend -> Draw("same");
      
      c -> Print(Form("%s/c_energy_vs_th__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c -> Print(Form("%s/c_energy_vs_th__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      delete c;
      
      
      c = new TCanvas(Form("c_energy_vs_Vov_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_energy_vs_Vov_%s-%s",ch1.c_str(),ch2.c_str()));
      // gPad -> SetLogy();
      
      hPad = (TH1F*)( gPad->DrawFrame(0.,0.,10.,50.) );
      hPad -> SetTitle(";V_{ov} [V];energy [a.u.]");
      hPad -> Draw();
      gPad -> SetGridy();
      
      iter = 0;
      for(auto mapIt : thLabels)
      {
        std::string label1(Form("%s_%s",ch1.c_str(),mapIt.first.c_str()));
        std::string label2(Form("%s_%s",ch2.c_str(),mapIt.first.c_str()));
        TGraph* g_energy1 = g_energy_vs_Vov[label1];
        TGraph* g_energy2 = g_energy_vs_Vov[label2];
        
        g_energy1 -> SetLineColor(1+iter);
        g_energy1 -> SetMarkerColor(1+iter);
        g_energy1 -> SetMarkerStyle(20);
        g_energy1 -> Draw("PL,same");
        
        g_energy2 -> SetLineColor(1+iter);
        g_energy2 -> SetMarkerColor(1+iter);
        g_energy2 -> SetMarkerStyle(25);
        g_energy2 -> Draw("PL,same");
        
        latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kBlack+iter);
        latex -> Draw("same");
        
        ++iter;
      }
      
      legend -> Draw("same");
      
      c -> Print(Form("%s/c_energy_vs_Vov__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c -> Print(Form("%s/c_energy_vs_Vov__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      delete c;
    }
  }
  
  
  
  
  
  
  //------------------
  //--- draw 2nd plots
  std::map<std::string,float> CTRMeans;
  std::map<std::string,float> CTRSigmas;
  
  for(auto stepLabel : stepLabels)
  {
    for(auto pair : pairsVec)
    {
      std::string ch1 = pair.first;
      std::string ch2 = pair.second;
      std::string label1(Form("%s_%s",ch1.c_str(),stepLabel.c_str()));
      std::string label2(Form("%s_%s",ch2.c_str(),stepLabel.c_str()));
      std::string label12 = Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),stepLabel.c_str());
      
      
      //--------------------------------------------------------
      
      
      histo = (TH1F*)( inFile->Get(Form("h1_deltaT_raw_%s",label12.c_str())) );
      FindSmallestInterval(vals,histo,0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      CTRMeans[label12] = mean;
      CTRSigmas[label12] = effSigma;
      
      
      //--------------------------------------------------------
      
      
    //   c = new TCanvas(Form("c_totRatio_%s",label12.c_str()),Form("c_totRatio_%s",label12.c_str()));
    //   gPad -> SetLogy();
      
    //   histo = ( TH1F* )( inFile->Get(Form("h1_totRatio_%s",label12.c_str())) );
    //   histo -> SetTitle(Form(";%s ToT / %s ToT;entries",ch1.c_str(),ch2.c_str()));
    //   histo -> Draw("colz");
      
    //   TF1* fitFunc = new TF1(Form("fitFunc_totRatio_%s",label12.c_str()),"gaus",
    //                          histo->GetMean()-1.5*histo->GetRMS(),
    //                          histo->GetMean()+1.5*histo->GetRMS());
    //   fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
    //   histo -> Fit(fitFunc,"QRSL+");
    //   fitFunc -> SetLineColor(kRed);
    //   fitFunc -> SetLineWidth(1);
    //   fitFunc -> Draw("same");
      
    //   latex = new TLatex(0.55,0.85,Form("#sigma = %.1f %%",100.*fitFunc->GetParameter(2)));
    //   latex -> SetNDC();
    //   latex -> SetTextFont(42);
    //   latex -> SetTextSize(0.04);
    //   latex -> SetTextColor(kRed);
    //   latex -> Draw("same");
      
    //   c -> Print(Form("%s/totRatio/c_totRatio__%s.png",plotDir.c_str(),label12.c_str()));
    //   c -> Print(Form("%s/totRatio/c_totRatio__%s.pdf",plotDir.c_str(),label12.c_str()));
    //   delete c;
      
      
    //   c = new TCanvas(Form("c_energyRatio_%s",label12.c_str()),Form("c_energyRatio_%s",label12.c_str()));
    //   gPad -> SetLogy();
      
    //   histo = ( TH1F* )( inFile->Get(Form("h1_energyRatio_%s",label12.c_str())) );
    //   histo -> SetTitle(Form(";%s energy / %s energy;entries",ch1.c_str(),ch2.c_str()));
    //   histo -> Draw("colz");
      
    //   fitFunc = new TF1(Form("fitFunc_energyRatio_%s",label12.c_str()),"gaus",
    //                     histo->GetMean()-1.5*histo->GetRMS(),
    //                     histo->GetMean()+1.5*histo->GetRMS());
    //   fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
    //   histo -> Fit(fitFunc,"QRSL+");
    //   fitFunc -> SetLineColor(kRed);
    //   fitFunc -> SetLineWidth(1);
    //   fitFunc -> Draw("same");
      
    //   latex = new TLatex(0.55,0.85,Form("#sigma = %.1f %%",100.*fitFunc->GetParameter(2)));
    //   latex -> SetNDC();
    //   latex -> SetTextFont(42);
    //   latex -> SetTextSize(0.04);
    //   latex -> SetTextColor(kRed);
    //   latex -> Draw("same");
      
    //   c -> Print(Form("%s/energyRatio/c_energyRatio__%s.png",plotDir.c_str(),label12.c_str()));
    //   c -> Print(Form("%s/energyRatio/c_energyRatio__%s.pdf",plotDir.c_str(),label12.c_str()));
    //   delete c;
    }
  }
  
  
  
  
  
  
  //------------------
  //--- draw 3rd plots
  if( std::find(plots.begin(),plots.end(),3) != plots.end() )
  {
    for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float th = map_ths[stepLabel];
      std::string VovLabel(Form("Vov%.1f",Vov));
      std::string thLabel(Form("th%02.0f",th));
      
      for(auto pair : pairsVec)
      {
        std::string ch1 = pair.first;
        std::string ch2 = pair.second;
        std::string label1(Form("%s_%s",ch1.c_str(),stepLabel.c_str()));
        std::string label2(Form("%s_%s",ch2.c_str(),stepLabel.c_str()));
        std::string label12 = Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),stepLabel.c_str());
        
        
        //--------------------------------------------------------
        
        
        c = new TCanvas(Form("c_deltaT_vs_energyRatio_%s",label12.c_str()),Form("c_deltaT_vs_energyRatio_%s",label12.c_str()));
        
        histo = (TH1F*)( inFile->Get(Form("h1_energyRatio_%s",label12.c_str())) );
        prof = (TProfile*)( inFile->Get(Form("p1_deltaT_vs_energyRatio_%s",label12.c_str())) );
        prof -> SetTitle(Form(";%s energy / %s energy;#Deltat [ps]",ch2.c_str(),ch1.c_str()));
        prof -> GetXaxis() -> SetRangeUser(histo->GetMean()-3.*histo->GetRMS(),
                                           histo->GetMean()+3.*histo->GetRMS());
        prof -> Draw("");
        
        // float fitXMin = histo->GetMean() - 2.*histo->GetRMS();
        // float fitXMax = histo->GetMean() + 2.*histo->GetRMS();
        // fitFunc_energyCorr[label12] = new TF1(Form("fitFunc_energyCorr_%s",label12.c_str()),"pol4",fitXMin,fitXMax);
        // prof -> Fit(fitFunc_energyCorr[label12],"QNRS+");
        // fitFunc_energyCorr[label12] -> SetLineColor(kRed);
        // fitFunc_energyCorr[label12] -> SetLineWidth(2);
        // fitFunc_energyCorr[label12] -> Draw("same");
        
        c -> Print(Form("%s/CTR/c_deltaT_vs_energyRatio__%s.png",plotDir.c_str(),label12.c_str()));
        c -> Print(Form("%s/CTR/c_deltaT_vs_energyRatio__%s.pdf",plotDir.c_str(),label12.c_str()));
        delete c;
      }
    }
  }
  
  
  
  
  
  
  //------------------
  //--- draw 4th plots
  std::map<std::string,TGraphErrors*> g_tRes_effSigma_vs_th;
  std::map<std::string,TGraphErrors*> g_tRes_gaus_vs_th;
  std::map<std::string,TGraphErrors*> g_tRes_effSigma_vs_Vov;
  std::map<std::string,TGraphErrors*> g_tRes_gaus_vs_Vov;
  std::map<std::string,TGraphErrors*> g_tRes_effSigma_bestTh_vs_Vov;
  std::map<std::string,TGraphErrors*> g_tRes_gaus_bestTh_vs_Vov;
  
  std::map<std::string,TGraphErrors*> g_tRes_energyCorr_effSigma_vs_th;
  std::map<std::string,TGraphErrors*> g_tRes_energyCorr_gaus_vs_th;
  std::map<std::string,TGraphErrors*> g_tRes_energyCorr_effSigma_vs_Vov;
  std::map<std::string,TGraphErrors*> g_tRes_energyCorr_gaus_vs_Vov;
  std::map<std::string,TGraphErrors*> g_tRes_energyCorr_effSigma_bestTh_vs_Vov;
  std::map<std::string,TGraphErrors*> g_tRes_energyCorr_gaus_bestTh_vs_Vov;
  
  std::map<std::string,TGraphErrors*> g_slewRate_vs_th;
  std::map<std::string,TGraphErrors*> g_slewRateNormalized_vs_th;

  std::map<std::string,std::map<float,float> > tRes_gaus_bestTh;
  std::map<std::string,std::map<float,float> > tResErr_gaus_bestTh;
  std::map<std::string,std::map<float,float> > tRes_effSigma_bestTh;
  std::map<std::string,std::map<float,float> > tRes_energyCorr_gaus_bestTh;
  std::map<std::string,std::map<float,float> > tResErr_energyCorr_gaus_bestTh;
  std::map<std::string,std::map<float,float> > tRes_energyCorr_effSigma_bestTh;
  
  if( std::find(plots.begin(),plots.end(),4) != plots.end() )
  {
    for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float th = map_ths[stepLabel];
      std::string VovLabel(Form("Vov%.1f",Vov));
      std::string thLabel(Form("th%02.0f",th));
      
      // t_Vov = Vov;
      // t_th = th;
      
      int pairsIt = 0;
      for(auto pair : pairsVec)
      {
        std::string ch1 = pair.first;
        std::string ch2 = pair.second;
        std::string label1(Form("%s_%s",ch1.c_str(),stepLabel.c_str()));
        std::string label2(Form("%s_%s",ch2.c_str(),stepLabel.c_str()));
        std::string label12 = Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),stepLabel.c_str());
        
        // t_arrayID1 = opts.GetOpt<unsigned int>(Form("%s.arrayID",ch1));
        // t_barID1   = opts.GetOpt<unsigned int>(Form("%s.barID",ch1));
        // t_lrID1    = opts.GetOpt<unsigned int>(Form("%s.lrID",ch1));
        // t_chID1    = ch1;
        // t_arrayID2 = opts.GetOpt<unsigned int>(Form("%s.arrayID",ch2.c_str()));
        // t_barID2   = opts.GetOpt<unsigned int>(Form("%s.barID",ch2.c_str()));
        // t_lrID2    = opts.GetOpt<unsigned int>(Form("%s.lrID",ch2.c_str()));
        // t_chID2    = ch2;
        
        
        //--------------------------------------------------------
        
        c = new TCanvas(Form("c_deltaT_energyCorr_%s",label12.c_str()),Form("c_deltaT_energyCorr_%s",label12.c_str()));
        
        histo = (TH1F*)( inFile->Get(Form("h1_deltaT_energyCorr_%s",label12.c_str())) );
        histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-3.*histo->GetRMS(),
                                            histo->GetMean()+3.*histo->GetRMS());          
        histo -> SetTitle(Form(";energy-corrected #Deltat [ps];entries"));
        histo -> SetLineWidth(2);
        histo -> SetLineColor(kBlue);
        histo -> SetMarkerColor(kBlue);
        histo -> Rebin(2);
        histo -> Draw("");
        
        float fitXMin = CTRMeans[label12] - 1.5*CTRSigmas[label12];
        float fitXMax = CTRMeans[label12] + 1.5*CTRSigmas[label12];
        TF1* fitFunc = new TF1(Form("fitFunc_energyCorr_%s",label12.c_str()),"gaus",fitXMin,fitXMax);
        fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
        histo -> Fit(fitFunc,"QNRSL+","");
        histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
        histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-2.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.*fitFunc->GetParameter(2));
        histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-1.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+1.5*fitFunc->GetParameter(2));
        float fitXMin2 = CTRMeans[label12] - 5.*CTRSigmas[label12];
        float fitXMax2 = CTRMeans[label12] + 5.*CTRSigmas[label12];
        TF1* fitFunc2 = new TF1(Form("fitFunc2_energyCorr_%s",label12.c_str()),"gaus(0)+gaus(3)+gaus(6)",fitXMin2,fitXMax2);
        fitFunc2 -> SetParameter(0,fitFunc->GetParameter(0));
        fitFunc2 -> SetParameter(1,fitFunc->GetParameter(1));
        fitFunc2 -> SetParameter(2,fitFunc->GetParameter(2));
        fitFunc2 -> SetParameter(3,fitFunc->GetParameter(0)/3.);
        fitFunc2 -> FixParameter(4,fitFunc->GetParameter(1)-330.);
        fitFunc2 -> SetParameter(5,fitFunc->GetParameter(2));
        fitFunc2 -> SetParameter(6,fitFunc->GetParameter(0)/3.);
        fitFunc2 -> FixParameter(7,fitFunc->GetParameter(1)+360.);
        fitFunc2 -> SetParameter(8,fitFunc->GetParameter(2));
        histo -> Fit(fitFunc2,"QNRSL+");
        
        fitFunc2 -> SetLineColor(kBlue+1);
        fitFunc2 -> SetLineWidth(3);
        fitFunc2 -> Draw("same");
        
        // t_tRes    = fitFunc -> GetParameter(2);
        // t_tResErr = fitFunc -> GetParError(2);
        // outTree -> Fill();
        
        FindSmallestInterval(vals,histo,0.68);
        float mean = vals[0];
        float min = vals[4];
        float max = vals[5];
        float delta = max-min;
        float sigma = 0.5*delta;
        float effSigma = sigma;
        
        histo -> GetXaxis() -> SetRangeUser(mean-5.*sigma,mean+5.*sigma);
        
        latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corr.}^{eff} = %.0f ps}{#sigma_{corr.}^{gaus} = %.0f ps}",effSigma,fitFunc->GetParameter(2)));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kBlue);
        latex -> Draw("same");
        
        std::string label_vs_th(Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),VovLabel.c_str()));
        if( g_tRes_effSigma_vs_th[label_vs_th] == NULL )
        {
          g_tRes_effSigma_vs_th[label_vs_th] = new TGraphErrors();
          g_tRes_gaus_vs_th[label_vs_th] = new TGraphErrors();
          g_tRes_energyCorr_effSigma_vs_th[label_vs_th] = new TGraphErrors();
          g_tRes_energyCorr_gaus_vs_th[label_vs_th] = new TGraphErrors();
          
          g_slewRate_vs_th[label_vs_th] = new TGraphErrors();
          g_slewRateNormalized_vs_th[label_vs_th] = new TGraphErrors();
        }
        
        std::string label_vs_Vov(Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),thLabel.c_str()));
        if( g_tRes_effSigma_vs_Vov[label_vs_Vov] == NULL )
        {
          g_tRes_effSigma_vs_Vov[label_vs_Vov] = new TGraphErrors();
          g_tRes_gaus_vs_Vov[label_vs_Vov] = new TGraphErrors();
          g_tRes_energyCorr_effSigma_vs_Vov[label_vs_Vov] = new TGraphErrors();
          g_tRes_energyCorr_gaus_vs_Vov[label_vs_Vov] = new TGraphErrors();
        } 
        
        std::string label_bestTh_vs_Vov(Form("%s-%s_bestTh",ch1.c_str(),ch2.c_str()));
        if( g_tRes_energyCorr_effSigma_bestTh_vs_Vov[label_bestTh_vs_Vov] == NULL )
        {
          g_tRes_effSigma_bestTh_vs_Vov[label_bestTh_vs_Vov] = new TGraphErrors();
          g_tRes_gaus_bestTh_vs_Vov[label_bestTh_vs_Vov] = new TGraphErrors();
          g_tRes_energyCorr_effSigma_bestTh_vs_Vov[label_bestTh_vs_Vov] = new TGraphErrors();
          g_tRes_energyCorr_gaus_bestTh_vs_Vov[label_bestTh_vs_Vov] = new TGraphErrors();
        } 
        
        float corr = 1.;
        if( pairsMode.at(pairsIt) == 1 ) corr = 1./sqrt(2.);
        if( pairsMode.at(pairsIt) == 2 ) corr = 0.5;
        
        g_tRes_effSigma_vs_th[label_vs_th] -> SetPoint(g_tRes_effSigma_vs_th[label_vs_th]->GetN(),th,effSigma*corr);
        g_tRes_effSigma_vs_th[label_vs_th] -> SetPointError(g_tRes_effSigma_vs_th[label_vs_th]->GetN()-1,0.,5.);
        g_tRes_gaus_vs_th[label_vs_th] -> SetPoint(g_tRes_gaus_vs_th[label_vs_th]->GetN(),th,fitFunc->GetParameter(2)*corr);
        g_tRes_gaus_vs_th[label_vs_th] -> SetPointError(g_tRes_gaus_vs_th[label_vs_th]->GetN()-1,0.,fitFunc->GetParError(2)*corr);
        
        g_tRes_effSigma_vs_Vov[label_vs_Vov] -> SetPoint(g_tRes_effSigma_vs_Vov[label_vs_Vov]->GetN(),Vov,effSigma*corr);
        g_tRes_effSigma_vs_Vov[label_vs_Vov] -> SetPointError(g_tRes_effSigma_vs_Vov[label_vs_Vov]->GetN()-1,0.,5.);
        g_tRes_gaus_vs_Vov[label_vs_Vov] -> SetPoint(g_tRes_gaus_vs_Vov[label_vs_Vov]->GetN(),Vov,fitFunc->GetParameter(2)*corr);
        g_tRes_gaus_vs_Vov[label_vs_Vov] -> SetPointError(g_tRes_gaus_vs_Vov[label_vs_Vov]->GetN()-1,0.,fitFunc->GetParError(2)*corr);
        
        if( ( tRes_energyCorr_gaus_bestTh[label_bestTh_vs_Vov][Vov] != 0. && fitFunc->GetParameter(2)*corr < tRes_energyCorr_gaus_bestTh[label_bestTh_vs_Vov][Vov] ) ||
            tRes_energyCorr_gaus_bestTh[label_bestTh_vs_Vov][Vov] == 0 )
        {
          tRes_energyCorr_gaus_bestTh[label_bestTh_vs_Vov][Vov] = fitFunc->GetParameter(2)*corr;
          tResErr_energyCorr_gaus_bestTh[label_bestTh_vs_Vov][Vov] = fitFunc->GetParError(2)*corr;
        }
        if( ( tRes_energyCorr_effSigma_bestTh[label_bestTh_vs_Vov][Vov] != 0. && effSigma*corr < tRes_energyCorr_effSigma_bestTh[label_bestTh_vs_Vov][Vov] ) ||
            tRes_energyCorr_effSigma_bestTh[label_bestTh_vs_Vov][Vov] == 0. )
        {
          tRes_energyCorr_effSigma_bestTh[label_bestTh_vs_Vov][Vov] = effSigma*corr;
        }
        
        
        histo = (TH1F*)( inFile->Get(Form("h1_deltaT_%s",label12.c_str())) );
        histo -> SetLineWidth(2);
        histo -> SetLineColor(kRed);
        histo -> SetMarkerColor(kRed);
        histo -> Rebin(2);
        histo -> Draw("same");
        
        fitXMin = CTRMeans[label12] - 1.5*CTRSigmas[label12];
        fitXMax = CTRMeans[label12] + 1.5*CTRSigmas[label12];
        fitFunc = new TF1(Form("fitFunc_%s",label12.c_str()),"gaus",fitXMin,fitXMax);
        fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
        histo -> Fit(fitFunc,"QNRSL+","");
        histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
        histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-2.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.*fitFunc->GetParameter(2));
        histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-1.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+1.5*fitFunc->GetParameter(2));
        fitFunc2 = new TF1(Form("fitFunc2_energyCorr_%s",label12.c_str()),"gaus(0)+gaus(3)+gaus(6)",fitXMin2,fitXMax2);
        fitFunc2 -> SetParameter(0,fitFunc->GetParameter(0));
        fitFunc2 -> SetParameter(1,fitFunc->GetParameter(1));
        fitFunc2 -> SetParameter(2,fitFunc->GetParameter(2));
        fitFunc2 -> SetParameter(3,fitFunc->GetParameter(0)/3.);
        fitFunc2 -> FixParameter(4,fitFunc->GetParameter(1)-330.);
        fitFunc2 -> SetParameter(5,fitFunc->GetParameter(2));
        fitFunc2 -> SetParameter(6,fitFunc->GetParameter(0)/3.);
        fitFunc2 -> FixParameter(7,fitFunc->GetParameter(1)+360.);
        fitFunc2 -> SetParameter(8,fitFunc->GetParameter(2));
        histo -> Fit(fitFunc2,"QNRSL+");
        
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
        
        latex = new TLatex(0.20,0.85,Form("#splitline{#sigma_{raw}^{eff} = %.0f ps}{#sigma_{raw}^{gaus} = %.0f ps}",effSigma,fitFunc->GetParameter(2)));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kRed);
        latex -> Draw("same");
        
        g_tRes_energyCorr_effSigma_vs_th[label_vs_th] -> SetPoint(g_tRes_energyCorr_effSigma_vs_th[label_vs_th]->GetN(),th,effSigma*corr);
        g_tRes_energyCorr_effSigma_vs_th[label_vs_th] -> SetPointError(g_tRes_energyCorr_effSigma_vs_th[label_vs_th]->GetN()-1,0.,5.);
        g_tRes_energyCorr_gaus_vs_th[label_vs_th] -> SetPoint(g_tRes_energyCorr_gaus_vs_th[label_vs_th]->GetN(),th,fitFunc->GetParameter(2)*corr);
        g_tRes_energyCorr_gaus_vs_th[label_vs_th] -> SetPointError(g_tRes_energyCorr_gaus_vs_th[label_vs_th]->GetN()-1,0.,fitFunc->GetParError(2)*corr);
        
        g_slewRate_vs_th[label_vs_th] -> SetPoint(g_slewRate_vs_th[label_vs_th]->GetN(),fitFunc->GetParameter(1),th);
        g_slewRate_vs_th[label_vs_th] -> SetPointError(g_slewRate_vs_th[label_vs_th]->GetN()-1,fitFunc->GetParError(1),0.);
	
	// normalize pulse shape to mip amplitude
	fLandau1 = new TF1("fLandau1","landau", 0, 50);
	fLandau2 = new TF1("fLandau2","landau", 0, 50);
	int isBar1 = opts.GetOpt<int>(Form("%s.isBar",ch1.c_str()));
	int isBar2 = opts.GetOpt<int>(Form("%s.isBar",ch2.c_str()));
	if (!isBar1 && !isBar2){
	  hamp = (TH1F*)( inFile->Get(Form("h1_energy_%s",label1.c_str())) );
	  int chID = opts.GetOpt<int>(Form("%s.chID",ch1.c_str()));
	  float max1 = FindXMaximum(hamp,cut_energyAcc[chID][Vov],100.);
	  fLandau1->SetRange(max1*0.1, max1*1.5);
	  hamp->Fit(fLandau1,"QR");
	  norm = fLandau1->GetParameter(1);
	}
	else if (isBar1 || isBar2){
	  std::string channelL, channelR;
	  if (isBar1){
	    channelL= opts.GetOpt<std::string>(Form("%s.channelL",ch1.c_str()));
	    channelR = opts.GetOpt<std::string>(Form("%s.channelR",ch1.c_str()));
	  }
	  if (isBar2){
	    channelL = opts.GetOpt<std::string>(Form("%s.channelL",ch2.c_str()));
	    channelR = opts.GetOpt<std::string>(Form("%s.channelR",ch2.c_str()));
	  }
	  int chIDL = opts.GetOpt<int>(Form("%s.chID",channelL.c_str()));
	  int chIDR = opts.GetOpt<int>(Form("%s.chID",channelR.c_str()));
	  
	  hampL = (TH1F*)( inFile->Get(Form("h1_energy_%s_%s_%s",channelL.c_str(),VovLabel.c_str(), thLabel.c_str())) );
	  hampR = (TH1F*)( inFile->Get(Form("h1_energy_%s_%s_%s",channelR.c_str(), VovLabel.c_str(), thLabel.c_str())) );
	  float max1 = FindXMaximum(hampL,cut_energyAcc[chIDL][Vov],100.);
          fLandau1->SetRange(max1*0.1, max1*1.5);
	  hampL->Fit(fLandau1,"QR");
	  float max2 = FindXMaximum(hampR,cut_energyAcc[chIDR][Vov],100.);
          fLandau2->SetRange(max2*0.1, max2*1.5);
          hampR->Fit(fLandau2,"QR");
	  norm = 0.5 * (fLandau1->GetParameter(1) + fLandau2->GetParameter(1) );
	}

	g_slewRateNormalized_vs_th[label_vs_th] -> SetPoint(g_slewRateNormalized_vs_th[label_vs_th]->GetN(),fitFunc->GetParameter(1),th/norm);
        g_slewRateNormalized_vs_th[label_vs_th] -> SetPointError(g_slewRateNormalized_vs_th[label_vs_th]->GetN()-1,fitFunc->GetParError(1),0.);
	


        g_tRes_energyCorr_effSigma_vs_Vov[label_vs_Vov] -> SetPoint(g_tRes_energyCorr_effSigma_vs_Vov[label_vs_Vov]->GetN(),Vov,effSigma*corr);
        g_tRes_energyCorr_effSigma_vs_Vov[label_vs_Vov] -> SetPointError(g_tRes_energyCorr_effSigma_vs_Vov[label_vs_Vov]->GetN()-1,0.,5.);
        g_tRes_energyCorr_gaus_vs_Vov[label_vs_Vov] -> SetPoint(g_tRes_energyCorr_gaus_vs_Vov[label_vs_Vov]->GetN(),Vov,fitFunc->GetParameter(2)*corr);
        g_tRes_energyCorr_gaus_vs_Vov[label_vs_Vov] -> SetPointError(g_tRes_energyCorr_gaus_vs_Vov[label_vs_Vov]->GetN()-1,0.,fitFunc->GetParError(2)*corr);
        
        if( ( tRes_gaus_bestTh[label_bestTh_vs_Vov][Vov] != 0. && fitFunc->GetParameter(2)*corr < tRes_gaus_bestTh[label_bestTh_vs_Vov][Vov] ) ||
            tRes_gaus_bestTh[label_bestTh_vs_Vov][Vov] == 0. )
        {
          tRes_gaus_bestTh[label_bestTh_vs_Vov][Vov] = fitFunc->GetParameter(2)*corr;
          tResErr_gaus_bestTh[label_bestTh_vs_Vov][Vov] = fitFunc->GetParError(2)*corr;
        }
        if( ( tRes_effSigma_bestTh[label_bestTh_vs_Vov][Vov] != 0. && effSigma*corr < tRes_effSigma_bestTh[label_bestTh_vs_Vov][Vov] ) ||
            tRes_effSigma_bestTh[label_bestTh_vs_Vov][Vov] == 0. )
        {
          tRes_effSigma_bestTh[label_bestTh_vs_Vov][Vov] = effSigma*corr;
        }
        
        c -> Print(Form("%s/CTR_energyCorr/c_deltaT_energyCorr__%s.png",plotDir.c_str(),label12.c_str()));
        c -> Print(Form("%s/CTR_energyCorr/c_deltaT_energyCorr__%s.pdf",plotDir.c_str(),label12.c_str()));
        delete c;
        
        ++pairsIt;
      }
    }
    
    for(auto pair : pairsVec)
    {
      std::string ch1 = pair.first;
      std::string ch2 = pair.second;
      std::string label(Form("%s-%s_bestTh",ch1.c_str(),ch2.c_str()));
      
      for(auto mapIt : tRes_effSigma_bestTh[label])
      {
        if( mapIt.first == 0 ) continue;
        
        g_tRes_effSigma_bestTh_vs_Vov[label] -> SetPoint(g_tRes_effSigma_bestTh_vs_Vov[label]->GetN(),mapIt.first,tRes_effSigma_bestTh[label][mapIt.first]);
        g_tRes_energyCorr_effSigma_bestTh_vs_Vov[label] -> SetPoint(g_tRes_energyCorr_effSigma_bestTh_vs_Vov[label]->GetN(),mapIt.first,tRes_energyCorr_effSigma_bestTh[label][mapIt.first]);
        
        g_tRes_gaus_bestTh_vs_Vov[label] -> SetPoint(g_tRes_gaus_bestTh_vs_Vov[label]->GetN(),mapIt.first,tRes_gaus_bestTh[label][mapIt.first]);
        g_tRes_energyCorr_gaus_bestTh_vs_Vov[label] -> SetPoint(g_tRes_energyCorr_gaus_bestTh_vs_Vov[label]->GetN(),mapIt.first,tRes_energyCorr_gaus_bestTh[label][mapIt.first]);
      }
    }
    
    
    //--------------------------------------------------------
    
    
    int pairsIt = 0;
    for(auto pair : pairsVec)
    {
      std::string ch1 = pair.first;
      std::string ch2 = pair.second;
      
      c = new TCanvas(Form("c_tRes_vs_th_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_tRes_vs_th_%s-%s",ch1.c_str(),ch2.c_str()));
      // gPad -> SetLogy();
      
      TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,tResMin,64.,tResMax) );
      if( pairsMode.at(pairsIt) == 0 )
        hPad -> SetTitle(";threshold [DAC];#sigma_{t_{diff}} [ps]");
      if( pairsMode.at(pairsIt) == 1 )
        hPad -> SetTitle(";threshold [DAC];#sigma_{t_{diff}} / #sqrt{2} [ps]");
      if( pairsMode.at(pairsIt) == 2 )
        hPad -> SetTitle(";threshold [DAC];#sigma_{t_{diff}} / 2 [ps]");
      hPad -> Draw();
      gPad -> SetGridy();
      
      int iter = 0;
      for(auto mapIt : VovLabels)
      {
        std::string label(Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),mapIt.first.c_str()));
        TGraph* g_effSigma = g_tRes_effSigma_vs_th[label];
        TGraph* g_gaus = g_tRes_gaus_vs_th[label];
        TGraph* g_energyCorr_effSigma = g_tRes_energyCorr_effSigma_vs_th[label];
        TGraph* g_energyCorr_gaus = g_tRes_energyCorr_gaus_vs_th[label];
        
        g_effSigma -> SetLineColor(1+iter);
        g_effSigma -> SetMarkerColor(1+iter);
        g_effSigma -> SetMarkerStyle(20);
        // g_effSigma -> Draw("PL,same");
        
        g_gaus -> SetLineColor(1+iter);
        g_gaus -> SetMarkerColor(1+iter);
        g_gaus -> SetMarkerStyle(20);
        g_gaus -> Draw("PL,same");
        
        g_energyCorr_effSigma -> SetLineColor(1+iter);
        g_energyCorr_effSigma -> SetMarkerColor(1+iter);
        g_energyCorr_effSigma -> SetMarkerStyle(25);
        // g_energyCorr_effSigma -> Draw("PL,same");
        
        g_energyCorr_gaus -> SetLineColor(1+iter);
        g_energyCorr_gaus -> SetMarkerColor(1+iter);
        g_energyCorr_gaus -> SetMarkerStyle(25);
        g_energyCorr_gaus -> Draw("PL,same");
        
        latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kBlack+iter);
        latex -> Draw("same");
        
        ++iter;
      }
      
      c -> Print(Form("%s/c_tRes_vs_th__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c -> Print(Form("%s/c_tRes_vs_th__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      delete c;
      
      
      c = new TCanvas(Form("c_slewRate_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_slewRate_%s-%s",ch1.c_str(),ch2.c_str()));
      // gPad -> SetLogy();
      
      hPad = (TH1F*)( gPad->DrawFrame(-0.5,0.,2.5,64.) );
      hPad -> SetTitle(";#LT t_{diff} #GT [ns];threshold [DAC]");
      hPad -> Draw();
      gPad -> SetGridy();


      cc = new TCanvas(Form("c_slewRateNormalized_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_slewRateNormalized_%s-%s",ch1.c_str(),ch2.c_str()));
      TH1F *hPad2 = (TH1F*)( gPad->DrawFrame(-0.5,0.,2.5,15.) );
      hPad2 -> SetTitle(";#LT t_{diff} #GT [ns]; normalized threshold");
      hPad2 -> Draw();
      gPad -> SetGridy();

      ccc = new TCanvas(Form("c_dVdt_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_dVdt_%s-%s",ch1.c_str(),ch2.c_str()));
      TH1F *hPad3 = (TH1F*)( gPad->DrawFrame(-0.5,0.,64,500.) );
      hPad3 -> SetTitle("; threshold [DAC]; dV/dt [a.u.]");
      hPad3 -> Draw();
      gPad -> SetGridy();
      
      iter = 0;
      for(auto mapIt : VovLabels)
      {
        std::string label(Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),mapIt.first.c_str()));
        TGraphErrors* g_slewRate_final = new TGraphErrors();
        TGraph* g_slewRate = g_slewRate_vs_th[label];
        double x0,y0;
        g_slewRate -> GetPoint(0,x0,y0);
        for(int point = 0; point < g_slewRate->GetN(); ++point)
        {
          double x,y;
          g_slewRate -> GetPoint(point,x,y);
          g_slewRate_final -> SetPoint(point,fabs(x-x0)/1000.,y);
          std::cout << "x: " << x << "   y: " << y << "   y0: " << y0 << "   val: " << fabs(y-y0)/1000. << std::endl;
        }
        
	c->cd();
        g_slewRate_final -> SetLineColor(1+iter);
        g_slewRate_final -> SetMarkerColor(1+iter);
        g_slewRate_final -> SetMarkerStyle(20);
        g_slewRate_final -> Draw("PL,same");

        latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kBlack+iter);
        latex -> Draw("same");


	// normalized to amp
        TGraphErrors* g_slewRateNormalized_final = new TGraphErrors();
        TGraph* g_slewRateNormalized = g_slewRateNormalized_vs_th[label];
        g_slewRateNormalized -> GetPoint(0,x0,y0);
        for(int point = 0; point < g_slewRateNormalized->GetN(); ++point)
        {
          double x,y;
          g_slewRateNormalized -> GetPoint(point,x,y);
          g_slewRateNormalized_final -> SetPoint(point,fabs(x-x0)/1000.,y);
          //std::cout << "x: " << x << "   y: " << y << "   y0: " << y0 << "   val: " << fabs(y-y0)/1000. << std::endl;
        }
        
	cc->cd();
        g_slewRateNormalized_final -> SetLineColor(1+iter);
        g_slewRateNormalized_final -> SetMarkerColor(1+iter);
        g_slewRateNormalized_final -> SetMarkerStyle(20);
        g_slewRateNormalized_final -> Draw("PL,same");

	latex -> Draw("same");
        

	// derivative of the pulse shape vs threshold
	TGraphErrors* g_dVdt = new TGraphErrors();
	for(int point = 0; point < g_slewRate->GetN(); ++point)
	  {
	    double x,y;
	    g_slewRate -> GetPoint(point,x,y);
	    float delta = 50.; // ps
 	    float dVdt = ( g_slewRate->Eval(x-delta) - g_slewRate->Eval(x+delta))/(2*delta);
	    g_dVdt -> SetPoint(point, y, dVdt*1000.);
	    //std::cout << "x: " << x << "   y: " << y << "   y0: " << y0 << "   val: " << fabs(y-y0)/1000. << std::endl; 
	  }
	
	ccc->cd();
        g_dVdt -> SetLineColor(1+iter);
        g_dVdt -> SetMarkerColor(1+iter);
        g_dVdt -> SetMarkerStyle(20);
        g_dVdt -> Draw("PL,same");

        latex -> Draw("same");

        ++iter;
      }
      
      c -> Print(Form("%s/c_slewRate__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c -> Print(Form("%s/c_slewRate__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c -> Print(Form("%s/c_slewRate__%s-%s.C",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      cc -> Print(Form("%s/c_slewRateNormalized__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      cc -> Print(Form("%s/c_slewRateNormalized__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      ccc -> Print(Form("%s/c_dVdt__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      ccc -> Print(Form("%s/c_dVdt__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      ccc -> Print(Form("%s/c_dVdt__%s-%s.C",plotDir.c_str(),ch1.c_str(),ch2.c_str()));

      delete c;
      delete cc;
      delete ccc;
      
      
      ++pairsIt;
    }
    
    
    
    //--------------------------------------------------------
    
    
    pairsIt = 0;
    for(auto pair : pairsVec)
    {
      std::string ch1 = pair.first;
      std::string ch2 = pair.second;
      
      c = new TCanvas(Form("c_tRes_vs_Vov_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_tRes_vs_Vov_%s-%s",ch1.c_str(),ch2.c_str()));
      // gPad -> SetLogy();
      
      TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,tResMin,10.,tResMax) );
      if( pairsMode.at(pairsIt) == 0 )
        hPad -> SetTitle(";V_{OV} [V];#sigma_{t_{diff}} [ps]");
      if( pairsMode.at(pairsIt) == 1 )
        hPad -> SetTitle(";V_{OV} [V];#sigma_{t_{diff}} / #sqrt{2} [ps]");
      if( pairsMode.at(pairsIt) == 2 )
        hPad -> SetTitle(";V_{OV} [V];#sigma_{t_{diff}} / 2 [ps]");
      hPad -> Draw();
      gPad -> SetGridy();
      
      int iter = 0;
      for(auto mapIt : thLabels)
      {
        std::string label(Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),mapIt.first.c_str()));
        TGraph* g_effSigma = g_tRes_effSigma_vs_Vov[label];
        TGraph* g_gaus = g_tRes_gaus_vs_Vov[label];
        TGraph* g_energyCorr_effSigma = g_tRes_energyCorr_effSigma_vs_Vov[label];
        TGraph* g_energyCorr_gaus = g_tRes_energyCorr_gaus_vs_Vov[label];
        
        g_effSigma -> SetLineColor(1+iter);
        g_effSigma -> SetMarkerColor(1+iter);
        g_effSigma -> SetMarkerStyle(20);
        // g_effSigma -> Draw("PL,same");
        
        g_gaus -> SetLineColor(1+iter);
        g_gaus -> SetMarkerColor(1+iter);
        g_gaus -> SetMarkerStyle(20);
        g_gaus -> Draw("PL,same");
        
        g_energyCorr_effSigma -> SetLineColor(1+iter);
        g_energyCorr_effSigma -> SetMarkerColor(1+iter);
        g_energyCorr_effSigma -> SetMarkerStyle(25);
        // g_energyCorr_effSigma -> Draw("PL,same");
        
        g_energyCorr_gaus -> SetLineColor(1+iter);
        g_energyCorr_gaus -> SetMarkerColor(1+iter);
        g_energyCorr_gaus -> SetMarkerStyle(25);
        g_energyCorr_gaus -> Draw("PL,same");
        
        latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kBlack+iter);
        latex -> Draw("same");
        
        ++iter;
      }
      
      c -> Print(Form("%s/c_tRes_vs_Vov__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c -> Print(Form("%s/c_tRes_vs_Vov__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      delete c;
      
      ++pairsIt;
    }
    
    
    //--------------------------------------------------------
    
    
    pairsIt = 0;
    for(auto pair : pairsVec)
    {
      std::string ch1 = pair.first;
      std::string ch2 = pair.second;
      
      c = new TCanvas(Form("c_tRes_bestTh_vs_Vov_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_tRes_bestTh_vs_Vov_%s-%s",ch1.c_str(),ch2.c_str()));
      // gPad -> SetLogy();

      legend2 = new TLegend(0.20,0.82,0.50,0.90);
      legend2 -> SetFillColor(0);
      legend2 -> SetFillStyle(1000);
      legend2 -> SetTextFont(42);
      
      TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,tResMin,10.,tResMax) );
      if( pairsMode.at(pairsIt) == 0 )
        hPad -> SetTitle(";V_{OV} [V];#sigma_{t_{diff}} [ps]");
      if( pairsMode.at(pairsIt) == 1 )
        hPad -> SetTitle(";V_{OV} [V];#sigma_{t_{diff}} / #sqrt{2} [ps]");
      if( pairsMode.at(pairsIt) == 2 )
        hPad -> SetTitle(";V_{OV} [V];#sigma_{t_{diff}} / 2 [ps]");
      hPad -> Draw();
      gPad -> SetGridy();
      
      std::string label(Form("%s-%s_bestTh",ch1.c_str(),ch2.c_str()));
      TGraph* g_effSigma = g_tRes_effSigma_bestTh_vs_Vov[label];
      TGraph* g_gaus = g_tRes_gaus_bestTh_vs_Vov[label];
      TGraph* g_energyCorr_effSigma = g_tRes_energyCorr_effSigma_bestTh_vs_Vov[label];
      TGraph* g_energyCorr_gaus = g_tRes_energyCorr_gaus_bestTh_vs_Vov[label];
      
      g_effSigma -> SetLineColor(1);
      g_effSigma -> SetMarkerColor(1);
      g_effSigma -> SetMarkerStyle(25);
      // g_effSigma -> Draw("PL,same");
      
      g_gaus -> SetLineColor(1);
      g_gaus -> SetMarkerColor(1);
      g_gaus -> SetMarkerStyle(25);
      g_gaus -> Draw("PL,same");
      
      g_energyCorr_effSigma -> SetLineColor(1);
      g_energyCorr_effSigma -> SetMarkerColor(1);
      g_energyCorr_effSigma -> SetMarkerStyle(20);
      // g_energyCorr_effSigma -> Draw("PL,same");
      
      g_energyCorr_gaus -> SetLineColor(1);
      g_energyCorr_gaus -> SetMarkerColor(1);
      g_energyCorr_gaus -> SetMarkerStyle(20);
      g_energyCorr_gaus -> Draw("PL,same");
      
      latex = new TLatex(0.55,0.85,Form("best th."));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack);
      latex -> Draw("same");
      
      legend2->AddEntry(g_gaus,"raw","PL");
      legend2->AddEntry(g_energyCorr_gaus,"corrected","PL");
      
      legend2->Draw("same");

      c -> Print(Form("%s/c_tRes_bestTh_vs_Vov__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c -> Print(Form("%s/c_tRes_bestTh_vs_Vov__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      delete c;
      
      ++pairsIt;
    }
  }
}
