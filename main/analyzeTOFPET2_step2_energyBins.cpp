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






int main(int argc, char** argv)
{
  setTDRStyle();
  
  
  if( argc < 2 )
  {
    std::cout << ">>> analyzeBars::usage:   " << argv[0] << " configFile.cfg" << std::endl;
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
  system(Form("mkdir -p %s/tracker/",plotDir.c_str()));
  system(Form("mkdir -p %s/qfine/",plotDir.c_str()));
  system(Form("mkdir -p %s/tot/",plotDir.c_str()));
  system(Form("mkdir -p %s/totRatio/",plotDir.c_str()));
  system(Form("mkdir -p %s/energy/",plotDir.c_str()));
  system(Form("mkdir -p %s/energyRatio/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR_energyCorr/",plotDir.c_str()));
  
  system(Form("cp %s/../index.php %s",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/tracker/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/qfine/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/tot/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/totRatio/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/energy/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/energyRatio/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/CTR/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/CTR_energyCorr/",plotDir.c_str(),plotDir.c_str()));
  
  //--- open files and make the tree chain
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  std::string fileBaseName = opts.GetOpt<std::string>("Input.fileBaseName");
  std::string runs = opts.GetOpt<std::string>("Input.runs");
  
  
  //--- track position cuts
  int doTracks = opts.GetOpt<int>("TrackCuts.doTracks");
  float cut_Xmin = -99;
  float cut_Xmax = 99;
  float cut_Ymin = -99;
  float cut_Ymax = 99;
  if( doTracks )
  {
    cut_Xmin = opts.GetOpt<float>("TrackCuts.Xmin");
    cut_Xmax = opts.GetOpt<float>("TrackCuts.Xmax");
    cut_Ymin = opts.GetOpt<float>("TrackCuts.Ymin");
    cut_Ymax = opts.GetOpt<float>("TrackCuts.Ymax");
  }
  
  
  //--- define channels
  std::vector<std::string> channels = opts.GetOpt<std::vector<std::string> >("Channels.channels");
  
  std::vector<std::string> pairs = opts.GetOpt<std::vector<std::string> >("Channels.pairs");
  std::vector<std::pair<std::string,std::string> > pairsVec;
  for(unsigned int ii = 0; ii < pairs.size()/2; ++ii)
  {
    pairsVec.push_back(std::make_pair(pairs.at(0+ii*2),pairs.at(1+ii*2)));
  }
  std::vector<int> pairsMode = opts.GetOpt<std::vector<int> >("Channels.pairsMode");
  
  
  //--- get cuts per bar / Vov
  std::map<unsigned int,std::map<float,float> > cut_qfineAcc;
  std::map<unsigned int,std::map<float,float> > cut_totAcc;
  std::map<unsigned int,std::map<float,float> > cut_energyAcc1;
  std::map<unsigned int,std::map<float,float> > cut_energyAcc2;
  std::map<unsigned int,std::map<float,float> > cut_energyFitMin;
  std::map<unsigned int,std::map<float,float> > cut_energyFitMax;
  std::map<unsigned int, float > noise;
  for(auto ch :  channels)
  {
    int isBar = opts.GetOpt<int>(Form("%s.isBar",ch.c_str()));

    if( isBar == 0 )
    {
      int chID = opts.GetOpt<int>(Form("%s.chID",ch.c_str()));
      std::vector<float> Vovs          = opts.GetOpt<std::vector<float> >(Form("%s.Vovs",ch.c_str()));
      std::vector<float> qfineMins     = opts.GetOpt<std::vector<float> >(Form("%s.qfineMins",ch.c_str()));
      std::vector<float> totMins       = opts.GetOpt<std::vector<float> >(Form("%s.totMins",ch.c_str()));
      std::vector<float> energyMins    = opts.GetOpt<std::vector<float> >(Form("%s.energyMins",ch.c_str()));
      std::vector<float> energyMaxs    = opts.GetOpt<std::vector<float> >(Form("%s.energyMaxs",ch.c_str()));
      std::vector<float> energyFitMins = opts.GetOpt<std::vector<float> >(Form("%s.energyFitMins",ch.c_str()));
      std::vector<float> energyFitMaxs = opts.GetOpt<std::vector<float> >(Form("%s.energyFitMaxs",ch.c_str()));
      noise[chID] = opts.GetOpt<float>(Form("%s.noise",ch.c_str()));
      int iter = 0;
      for(auto Vov : Vovs)
      {
        cut_qfineAcc[chID][Vov]     = qfineMins.at(iter);
        cut_totAcc[chID][Vov]       = totMins.at(iter);
        cut_energyAcc1[chID][Vov]   = energyMins.at(iter);
        cut_energyAcc2[chID][Vov]   = energyMaxs.at(iter);
        cut_energyFitMin[chID][Vov] = energyFitMins.at(iter);
        cut_energyFitMax[chID][Vov] = energyFitMaxs.at(iter);
      ++iter;
      }
    }
  }
  std::map<std::string,float> cut_energyMin;
  std::map<std::string,float> cut_energyMax;
  
  
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
  
    std::size_t found;
    
    found = name.find("data_");
    if( found!=std::string::npos )
    {
      std::string label(Form("%s_%s_%s",tokens[1].c_str(),tokens[2].c_str(),tokens[3].c_str()));
      trees[label] = (TTree*)( inFile->Get(name.c_str()) );
    }
    
    found = name.find("h1_energy_cut");
    if( found!=std::string::npos )
    {
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
  
  
  //--- define histograms
  std::string plotFileName = opts.GetOpt<std::string>("Output.plotFileName");
  TFile* outFile = TFile::Open(plotFileName.c_str(),"RECREATE");
  
  std::map<std::string,TH1F*> h1_energy_allTh;
  
  std::map<std::string,TH2F*> h2_tot_corr;
  std::map<std::string,TH1F*> h1_totRatio;
  
  std::map<std::string,TH2F*> h2_energy_corr;
  std::map<std::string,TH1F*> h1_energyRatio;

  std::map<std::string,TH1F*> h1_deltaT_raw;
  std::map<std::string,TH1F*> h1_deltaT;
  std::map<std::string,TProfile*> p1_deltaT_vs_energyRatio;

  std::map<std::string,TH1F*> h1_deltaT_energyCorr;
  std::map<std::string,TProfile*> p1_deltaT_energyCorr_vs_pos;

  std::map<std::string,TH1F*> h1_deltaT_energyCorr_posCorr;
  
  
  
  
  //--- get plot settings
  std::vector<int> plots = opts.GetOpt<std::vector<int> >("Plots.plots");
  int photopeakSelection = opts.GetOpt<int>("Plots.photopeakSelection");
  float energyMin = opts.GetOpt<float>("Plots.energyMin");
  float energyMax = opts.GetOpt<float>("Plots.energyMax");
  float tResMin = opts.GetOpt<float>("Plots.tResMin");
  float tResMax = opts.GetOpt<float>("Plots.tResMax");
  
  TCanvas* c;
  TCanvas* c2;
  TCanvas* c3;
  TCanvas *c4;
  float* vals = new float[6];
  TLatex* latex;
  TH1F* histo;
  TH2F* histo2;
  TProfile* prof;
  TProfile2D* prof2;
  TLegend* legend;
  TLegend* legend2;
  
  
  
  
  //-------------------
  // define energy bins
  int nEnergyBins = opts.GetOpt<float>("Plots.nEnergyBins");
  std::vector<float> energyBins = opts.GetOpt<std::vector<float> >("Plots.energyBins");
  
  std::map<std::string,float*> energyBinEdges;
  std::map<std::string,int> nEnergyBinsEff;
  unsigned int it = 0;
  while( it < energyBins.size() )
  {
    float ov = energyBins.at(it);
    std::string VovLabel(Form("Vov%.1f",ov));
    
    ++it;
    nEnergyBinsEff[VovLabel] = energyBins.at(it);
    ++it;
    
    energyBinEdges[VovLabel] = new float[1+nEnergyBinsEff[VovLabel]];
    
    
    std::cout << "ov: " << ov << "   nBins: " << nEnergyBinsEff[VovLabel] << "   ";
    for(int jj = 0; jj < nEnergyBinsEff[VovLabel]+1; ++jj)
    {
      std::cout << energyBins.at(it) << ",";
      energyBinEdges[VovLabel][jj] = energyBins.at(it);
      ++it;
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  
  
  //------------------
  //--- draw 1st plots
  std::map<std::string,TGraphErrors*> g_tot_vs_th;
  std::map<std::string,TGraphErrors*> g_tot_vs_Vov;
  std::map<std::string,TGraphErrors*> g_energy_vs_th;
  std::map<std::string,TGraphErrors*> g_energy_vs_Vov;
  
  std::map<std::string,TH1F*> h1_energyBins;
  std::map<std::string,std::map<int,TH1F*> > h1_energy_bin;
  
  if( std::find(plots.begin(),plots.end(),1) != plots.end() )
  {
    for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float th = map_ths[stepLabel];
      std::string VovLabel(Form("Vov%.1f",Vov));
      std::string thLabel(Form("th%02.0f",th));
      
      
      //--------------------------------------------------------
      
      
      for(auto ch : channels)
      {
        int isBar = opts.GetOpt<int>(Form("%s.isBar",ch.c_str()));
        int chID = -1;
        
        if( isBar == 1 )
        {
          std::string channelL = opts.GetOpt<std::string>(Form("%s.channelL",ch.c_str()));
          chID = opts.GetOpt<int>(Form("%s.chID",channelL.c_str()));
        }
        else
          int chID = opts.GetOpt<int>(Form("%s.chID",ch.c_str()));  
        std::string label(Form("%s_%s",ch.c_str(),stepLabel.c_str()));
        
        
        if( doTracks )
        {
          c = new TCanvas(Form("c_eff_vs_XY_%s",label.c_str()),Form("c_eff_vs_XY_%s",label.c_str()));
          
          prof2 = (TProfile2D*)( inFile->Get(Form("p2_eff_vs_XY_%s",label.c_str())) );
          prof2 -> SetTitle(";tracker x [mm];tracker y [mm]");
          prof2 -> Draw("colz");
          prof2 -> GetXaxis() -> SetRangeUser(-10.,40.);
          prof2 -> GetYaxis() -> SetRangeUser(5.,40.);
          TLine* line1 = new TLine(cut_Xmin,cut_Ymin,cut_Xmax,cut_Ymin);
          line1 -> SetLineColor(kMagenta);
          line1 -> Draw("same");
          TLine* line2 = new TLine(cut_Xmin,cut_Ymax,cut_Xmax,cut_Ymax);
          line2 -> SetLineColor(kMagenta);
          line2 -> Draw("same");
          TLine* line3 = new TLine(cut_Xmin,cut_Ymin,cut_Xmin,cut_Ymax);
          line3 -> SetLineColor(kMagenta);
          line3 -> Draw("same");
          TLine* line4 = new TLine(cut_Xmax,cut_Ymin,cut_Xmax,cut_Ymax);
          line4 -> SetLineColor(kMagenta);
          line4 -> Draw("same");
          latex = new TLatex(0.65,0.85,Form("%s",ch.c_str()));
          latex -> SetNDC();
          latex -> SetTextFont(42);
          latex -> SetTextSize(0.04);
          latex -> Draw("same");      
          prof2 -> Write();
          c -> Print(Form("%s/tracker/c_eff_vs_XY__%s.png",plotDir.c_str(),label.c_str()));
          c -> Print(Form("%s/tracker/c_eff_vs_XY__%s.pdf",plotDir.c_str(),label.c_str()));
          delete c;
        }
        
        c = new TCanvas(Form("c_qfine_%s",label.c_str()),Form("c_qfine_%s",label.c_str()));
        gPad -> SetLogy();
        
        histo = (TH1F*)( inFile->Get(Form("h1_qfine_%s",label.c_str())) );
        histo -> SetTitle(";Q_{fine} [ADC];entries");
        histo -> SetLineColor(kRed);
        histo -> Draw();
        histo -> GetXaxis() -> SetRangeUser(0.,600.);
        TLine* line_qfineAcc1 = new TLine(cut_qfineAcc[chID][Vov],histo->GetMinimum(),cut_qfineAcc[chID][Vov],histo->GetMaximum());
        line_qfineAcc1 -> SetLineColor(kBlack);
        line_qfineAcc1 -> Draw("same");
        latex = new TLatex(0.65,0.85,Form("%s",ch.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kRed);
        latex -> Draw("same");      
        histo -> Write();
        c -> Print(Form("%s/qfine/c_qfine__%s.png",plotDir.c_str(),label.c_str()));
        c -> Print(Form("%s/qfine/c_qfine__%s.pdf",plotDir.c_str(),label.c_str()));
        delete c;
        
        
        histo = (TH1F*)( inFile->Get(Form("h1_qfine_%s",label.c_str())) );      
        c = new TCanvas(Form("c_qfine_vs_tot_%s",label.c_str()),Form("c_qfine_vs_tot_%s",label.c_str()));
        gPad -> SetLogz();
        
        
        histo2 = (TH2F*)( inFile->Get(Form("h2_qfine_vs_tot_%s",label.c_str())) );            
        histo2 -> SetTitle(Form(";%s ToT [ns];%s Q_{fine} [ADC]",ch.c_str(),ch.c_str()));
        histo2 -> Draw("colz");
        
        histo -> Write();
        histo2 -> Write();
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
        histo -> Write();
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
        
        
        
        if( h1_energy_allTh[ch+"_"+VovLabel] == NULL )
          h1_energy_allTh[ch+"_"+VovLabel] = new TH1F(Form("h1_energy_allTh_%s_%s",ch.c_str(),VovLabel.c_str()),"",400,0.,400.);
        histo = (TH1F*)( inFile->Get(Form("h1_energy_%s",label.c_str())) );
        h1_energy_allTh[ch+"_"+VovLabel] -> Add(histo);
        
        
        
        c = new TCanvas(Form("c_energy_%s",label.c_str()),Form("c_energy_%s",label.c_str()));
        gPad -> SetLogy();
        
        histo = (TH1F*)( inFile->Get(Form("h1_energy_%s",label.c_str())) );      
        histo -> SetTitle(";energy [a.u.];entries");
        histo -> SetLineColor(kRed);
        histo -> Draw();
        max1 = FindXMaximum(histo,cut_energyAcc1[chID][Vov],200.);
        histo -> GetXaxis() -> SetRangeUser(energyMin,energyMax);
        fitFunc1 = new TF1("fitFunc1","gaus",max1-cut_energyFitMin[chID][Vov]*max1,max1+cut_energyFitMax[chID][Vov]*max1);
        histo -> Fit(fitFunc1,"QNRS+");
        fitFunc1 -> SetLineColor(kBlack);
        fitFunc1 -> SetLineWidth(3);
        fitFunc1 -> Draw("same");
        if( photopeakSelection )
        {
          cut_energyMin[Form("%s_%s",ch.c_str(),stepLabel.c_str())] = fitFunc1->GetMaximumX()-cut_energyFitMin[chID][Vov]*fitFunc1->GetMaximumX();
          cut_energyMax[Form("%s_%s",ch.c_str(),stepLabel.c_str())] = fitFunc1->GetMaximumX()+cut_energyFitMax[chID][Vov]*fitFunc1->GetMaximumX();
        }
        else
        {
          cut_energyMin[Form("%s_%s",ch.c_str(),stepLabel.c_str())] = cut_energyAcc1[chID][Vov];
          cut_energyMax[Form("%s_%s",ch.c_str(),stepLabel.c_str())] = cut_energyAcc2[chID][Vov];
        }
        
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
        histo -> Write();
        
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
      int size = VovLabels.size();
      for(auto mapIt : VovLabels) 
      {
        std::string label1(Form("%s_%s",ch1.c_str(),mapIt.first.c_str()));
        std::string label2(Form("%s_%s",ch2.c_str(),mapIt.first.c_str()));
        TGraph* g_tot1 = g_tot_vs_th[label1];
        TGraph* g_tot2 = g_tot_vs_th[label2];
        
        g_tot1 -> SetLineColor(51+iter*(int(50/size)));
        g_tot1 -> SetMarkerColor(51+iter*(int(50/size)));
        g_tot1 -> SetMarkerStyle(20);
        g_tot1 -> Draw("PL,same");
        
        g_tot2 -> SetLineColor(51+iter*(int(50/size)));
        g_tot2 -> SetMarkerColor(51+iter*(int(50/size)));
        g_tot2 -> SetMarkerStyle(25);
        g_tot2 -> Draw("PL,same");
        
        latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(51+iter*(int(50/size)));
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
      size = thLabels.size();
      for(auto mapIt : thLabels)
      {
        std::string label1(Form("%s_%s",ch1.c_str(),mapIt.first.c_str()));
        std::string label2(Form("%s_%s",ch2.c_str(),mapIt.first.c_str()));
        TGraph* g_tot1 = g_tot_vs_Vov[label1];
        TGraph* g_tot2 = g_tot_vs_Vov[label2];
        
        g_tot1 -> SetLineColor(51+iter*(int(50/size)));
        g_tot1 -> SetMarkerColor(51+iter*(int(50/size)));
        g_tot1 -> SetMarkerStyle(20);
        g_tot1 -> Draw("PL,same");
        
        g_tot2 -> SetLineColor(51+iter*(int(50/size)));
        g_tot2 -> SetMarkerColor(51+iter*(int(50/size)));
        g_tot2 -> SetMarkerStyle(25);
        g_tot2 -> Draw("PL,same");
        
        latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(51+iter*(int(50/size)));
        latex -> Draw("same");
        
        ++iter;
      }
      
      legend -> Draw("same");
      
      c -> Print(Form("%s/c_tot_vs_Vov__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c -> Print(Form("%s/c_tot_vs_Vov__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      delete c;
      
      
      c = new TCanvas(Form("c_energy_vs_th_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_energy_vs_th_%s-%s",ch1.c_str(),ch2.c_str()));
      // gPad -> SetLogy();
      
      hPad = (TH1F*)( gPad->DrawFrame(-1.,energyMin,64.,energyMax) );
      hPad -> SetTitle(";threshold [DAC];energy [a.u.]");
      hPad -> Draw();
      gPad -> SetGridy();
      
      iter = 0;
      size = VovLabels.size();
      for(auto mapIt : VovLabels)
      {
        std::string label1(Form("%s_%s",ch1.c_str(),mapIt.first.c_str()));
        std::string label2(Form("%s_%s",ch2.c_str(),mapIt.first.c_str()));
        TGraph* g_energy1 = g_energy_vs_th[label1];
        TGraph* g_energy2 = g_energy_vs_th[label2];
        
        g_energy1 -> SetLineColor(51+iter*(int(50/size)));
        g_energy1 -> SetMarkerColor(51+iter*(int(50/size)));
        g_energy1 -> SetMarkerStyle(20);
        g_energy1 -> Draw("PL,same");
        
        g_energy2 -> SetLineColor(51+iter*(int(50/size)));
        g_energy2 -> SetMarkerColor(51+iter*(int(50/size)));
        g_energy2 -> SetMarkerStyle(25);
        g_energy2 -> Draw("PL,same");
        
        latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(51+iter*(int(50/size)));
        latex -> Draw("same");
        
        ++iter;
      }
      
      legend -> Draw("same");
      
      c -> Print(Form("%s/c_energy_vs_th__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c -> Print(Form("%s/c_energy_vs_th__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      delete c;
      
      
      c = new TCanvas(Form("c_energy_vs_Vov_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_energy_vs_Vov_%s-%s",ch1.c_str(),ch2.c_str()));
      // gPad -> SetLogy();
      
      hPad = (TH1F*)( gPad->DrawFrame(0.,energyMin,10.,energyMax) );
      hPad -> SetTitle(";V_{ov} [V];energy [a.u.]");
      hPad -> Draw();
      gPad -> SetGridy();
      
      iter = 0;
      size = thLabels.size();
      for(auto mapIt : thLabels)
      {
        std::string label1(Form("%s_%s",ch1.c_str(),mapIt.first.c_str()));
        std::string label2(Form("%s_%s",ch2.c_str(),mapIt.first.c_str()));
        TGraph* g_energy1 = g_energy_vs_Vov[label1];
        TGraph* g_energy2 = g_energy_vs_Vov[label2];
        
        g_energy1 -> SetLineColor(51+iter*(int(50/size)));
        g_energy1 -> SetMarkerColor(51+iter*(int(50/size)));
        g_energy1 -> SetMarkerStyle(20);
        g_energy1 -> Draw("PL,same");
        
        g_energy2 -> SetLineColor(51+iter*(int(50/size)));
        g_energy2 -> SetMarkerColor(51+iter*(int(50/size)));
        g_energy2 -> SetMarkerStyle(25);
        g_energy2 -> Draw("PL,same");
        
        latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(51+iter*(int(50/size)));
        latex -> Draw("same");
        
        ++iter;
      }
      
      legend -> Draw("same");
      
      c -> Print(Form("%s/c_energy_vs_Vov__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c -> Print(Form("%s/c_energy_vs_Vov__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      delete c;
    }
    
    
    for(auto ch : channels)
    {
      for(auto mapIt : VovLabels)
      {
        //int chID = opts.GetOpt<int>(Form("%s.chID",ch.c_str()));  
        
        std::string VovLabel = mapIt.first;
        std::string string_Vov = VovLabel;
        string_Vov.erase(0,3);
        float Vov = atof(string_Vov.c_str());
      
        c = new TCanvas(Form("c_energy_allTh_%s_%s",ch.c_str(),VovLabel.c_str()),Form("c_energy_allTh_%s_%s",ch.c_str(),VovLabel.c_str()));
        gPad -> SetLogy();
        
        histo = h1_energy_allTh[ch+"_"+VovLabel];
        
        histo -> SetTitle(";energy [a.u.];entries");
        histo -> SetLineColor(kRed);
        histo -> Draw();
        
        // max1 = FindXMaximum(histo,cut_energyA<cc1[chID][Vov],200.);
        histo -> GetXaxis() -> SetRangeUser(energyMin,energyMax);
        // fitFunc1 = new TF1("fitFunc1","gaus",max1-cut_energyFitMin[chID][Vov]*max1,max1+cut_energyFitMax[chID][Vov]*max1);
        // histo -> Fit(fitFunc1,"QNRS+");
        // fitFunc1 -> SetLineColor(kBlack);
        // fitFunc1 -> SetLineWidth(3);
        // fitFunc1 -> Draw("same");
        // if( photopeakSelection )
        // {
        //   cut_energyMin[Form("%s_%s",ch.c_str(),stepLabel.c_str())] = fitFunc1->GetMaximumX()-cut_energyFitMin[chID][Vov]*fitFunc1->GetMaximumX();
        //   cut_energyMax[Form("%s_%s",ch.c_str(),stepLabel.c_str())] = fitFunc1->GetMaximumX()+cut_energyFitMax[chID][Vov]*fitFunc1->GetMaximumX();
        // }
        // else
        // {
        //   cut_energyMin[Form("%s_%s",ch.c_str(),stepLabel.c_str())] = cut_energyAcc1[chID][Vov];
        //   cut_energyMax[Form("%s_%s",ch.c_str(),stepLabel.c_str())] = cut_energyAcc2[chID][Vov];
        // }
        
        
        
        // // find energy bins
        // int firstBin = histo->FindBin(cut_energyAcc1[chID][Vov]);
        // int lastBin = histo->FindBin(cut_energyAcc2[chID][Vov]);
        
        // int nEventsTot = histo -> Integral(firstBin,lastBin);
        // int nEventsBin = nEventsTot/nEnergyBins;
        
        // float* energyBinEdges = new float[1+nEnergyBins];
        
        // energyBinEdges[0] = histo->GetBinLowEdge(firstBin);
        // nEnergyBinsEff[ch+"_"+VovLabel] = 0;
        
        // float integral = 0;
        // int startBin = firstBin;
        // for(int ii = 0; ii < nEnergyBins; ++ii)
        // {
        //   for(int bin = startBin; bin <= lastBin; ++bin)
        //   {
        //     integral += histo->GetBinContent(bin);
            
        //     if( integral >= nEventsBin || bin == lastBin )
        //     {
        //       startBin = bin+1;
        //       energyBinEdges[ii+1] = histo->GetBinLowEdge(bin) + histo->GetBinWidth(bin);
        //       std::cout << ">>>>>> ii: " << ii << "   bin: " << bin << "   energyBinEdges: [" << energyBinEdges[ii] << "," << energyBinEdges[ii+1] << "]   integral: " << integral << std::endl;
        //       nEnergyBinsEff[ch+"_"+VovLabel] += 1;
        //       integral = 0;
              
        //       TLine* line_bin = new TLine(energyBinEdges[ii+1],histo->GetMinimum(),energyBinEdges[ii+1],histo->GetMaximum());
        //       line_bin -> SetLineColor(kBlack);
        //       line_bin -> SetLineStyle(2);
        //       line_bin -> Draw("same");
        //       break;
        //     }
        //   }
        // }
        // std::cout << ch + "_" + VovLabel << ":   " << (nEnergyBinsEff[ch+"_"+VovLabel]) << std::endl;
        // energyBinEdges[(nEnergyBinsEff[ch+"_"+VovLabel])] = histo->GetBinLowEdge(lastBin) + histo->GetBinWidth(lastBin);
        
        // h1_energyBins[Form("%s_%s",ch.c_str(),VovLabel.c_str())] = new TH1F(Form("h1_energyBins_%s_%s",ch.c_str(),VovLabel.c_str()),"",nEnergyBinsEff[ch+"_"+VovLabel],energyBinEdges);
        // for(int jj = 1; jj <= nEnergyBinsEff[ch+"_"+VovLabel]; ++jj)
        //   (h1_energy_bin[Form("%s_%s",ch.c_str(),VovLabel.c_str())])[jj] = new TH1F(Form("h1_energy_%s_%s_bin%d",ch.c_str(),VovLabel.c_str(),jj),"",1000,0.,1000.);
        
        for(int ii = 0; ii < nEnergyBinsEff[VovLabel]+1; ++ii)
        {
          TLine* line_bin = new TLine(energyBinEdges[VovLabel][ii],histo->GetMinimum(),energyBinEdges[VovLabel][ii],histo->GetMaximum());
          line_bin -> SetLineColor(kBlack);
          line_bin -> SetLineStyle(2);
          line_bin -> Draw("same");
        }
        
        h1_energyBins[Form("%s_%s",ch.c_str(),VovLabel.c_str())] = new TH1F(Form("h1_energyBins_%s_%s",ch.c_str(),VovLabel.c_str()),"",nEnergyBinsEff[VovLabel],energyBinEdges[VovLabel]);
        for(int jj = 1; jj <= nEnergyBinsEff[VovLabel]; ++jj)
          (h1_energy_bin[Form("%s_%s",ch.c_str(),VovLabel.c_str())])[jj] = new TH1F(Form("h1_energy_%s_%s_bin%d",ch.c_str(),VovLabel.c_str(),jj),"",1000,0.,1000.);
        
        
        c -> Print(Form("%s/energy/c_energy_allTh__%s_%s.png",plotDir.c_str(),ch.c_str(),VovLabel.c_str()));
        c -> Print(Form("%s/energy/c_energy_allTh__%s_%s.pdf",plotDir.c_str(),ch.c_str(),VovLabel.c_str()));
      }
    }
  }
  
  
  
  
  //------------------------
  //--- 2nd loop over events
  std::map<std::string,std::map<int,bool> > accept;
  
  for(auto mapIt : trees)
  {
    std::string label = mapIt.first;
    std::cout << ">>> 2nd loop: " << label << std::endl;
    
    std::vector<std::string> tokens = GetTokens(label,'_');
    std::string stepLabel(Form("%s_%s",tokens[1].c_str(),tokens[2].c_str()));
    
    float Vov = map_Vovs[stepLabel];
    std::string VovLabel(Form("Vov%.1f",Vov));
    
    EventClass* anEvent = new EventClass();
    mapIt.second -> SetBranchAddress("event",&anEvent);
    
    int nEntries = mapIt.second->GetEntries();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      // if( entry%1000 == 0 ) std::cout << ">>> 2nd loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
      mapIt.second -> GetEntry(entry);
      
      accept[label][entry] = false;
      
      // selection on track position
      if( doTracks )
      {
        if( (cut_Xmin != -999 && cut_Xmax != -999 ) && (anEvent->x < cut_Xmin || anEvent->x > cut_Xmax) ) continue;
        if( (cut_Ymin != -999 && cut_Ymax != -999 ) && (anEvent->y < cut_Ymin || anEvent->y > cut_Ymax) ) continue;
      }
      
      
      if( anEvent->isBar1 == 0 )
      {
        int chID1 = opts.GetOpt<int>(Form("%s.chID",anEvent->ch1.c_str()));
        
        if( anEvent->qfine1 < cut_qfineAcc[chID1][Vov] ) continue;
        if( anEvent->tot1 < cut_totAcc[chID1][Vov] ) continue;
        if( anEvent->energy1 < cut_energyMin[Form("%s_%s",anEvent->ch1.c_str(),anEvent->stepLabel.c_str())] ) continue;
        if( anEvent->energy1 > cut_energyMax[Form("%s_%s",anEvent->ch1.c_str(),anEvent->stepLabel.c_str())] ) continue;
      }
      else
      {
        std::string channelL = opts.GetOpt<std::string>(Form("%s.channelL",anEvent->ch1.c_str()));
        int chIDL = opts.GetOpt<int>(Form("%s.chID",channelL.c_str()));
        std::string channelR = opts.GetOpt<std::string>(Form("%s.channelR",anEvent->ch1.c_str()));
        int chIDR = opts.GetOpt<int>(Form("%s.chID",channelR.c_str()));
        
        if( anEvent->qfine1L < cut_qfineAcc[chIDL][Vov] ) continue;
        if( anEvent->tot1L < cut_totAcc[chIDL][Vov] ) continue;
        if( anEvent->energy1L < cut_energyMin[Form("%s_%s",channelL.c_str(),anEvent->stepLabel.c_str())] ) continue;
        if( anEvent->energy1L > cut_energyMax[Form("%s_%s",channelL.c_str(),anEvent->stepLabel.c_str())] ) continue;
        
        if( anEvent->qfine1R < cut_qfineAcc[chIDR][Vov] ) continue;
        if( anEvent->tot1R < cut_totAcc[chIDR][Vov] ) continue;
        if( anEvent->energy1R < cut_energyMin[Form("%s_%s",channelR.c_str(),anEvent->stepLabel.c_str())] ) continue;
        if( anEvent->energy1R > cut_energyMax[Form("%s_%s",channelR.c_str(),anEvent->stepLabel.c_str())] ) continue;
      }
      
      
      if( anEvent->isBar2 == 0 )
      {
        int chID2 = opts.GetOpt<int>(Form("%s.chID",anEvent->ch2.c_str()));
        
        if( anEvent->qfine2 < cut_qfineAcc[chID2][Vov] ) continue;
        if( anEvent->tot2 < cut_totAcc[chID2][Vov] ) continue;
        if( anEvent->energy2 < cut_energyMin[Form("%s_%s",anEvent->ch2.c_str(),anEvent->stepLabel.c_str())] ) continue;
        if( anEvent->energy2 > cut_energyMax[Form("%s_%s",anEvent->ch2.c_str(),anEvent->stepLabel.c_str())] ) continue;
      }
      else
      {
        std::string channelL = opts.GetOpt<std::string>(Form("%s.channelL",anEvent->ch2.c_str()));
        int chIDL = opts.GetOpt<int>(Form("%s.chID",channelL.c_str()));
        std::string channelR = opts.GetOpt<std::string>(Form("%s.channelR",anEvent->ch2.c_str()));
        int chIDR = opts.GetOpt<int>(Form("%s.chID",channelR.c_str()));
        
        if( anEvent->qfine2L < cut_qfineAcc[chIDL][Vov] ) continue;
        if( anEvent->tot2L < cut_totAcc[chIDL][Vov] ) continue;
        if( anEvent->energy2L < cut_energyMin[Form("%s_%s",channelL.c_str(),anEvent->stepLabel.c_str())] ) continue;
        if( anEvent->energy2L > cut_energyMax[Form("%s_%s",channelL.c_str(),anEvent->stepLabel.c_str())] ) continue;

        if( anEvent->qfine2R < cut_qfineAcc[chIDR][Vov] ) continue;
        if( anEvent->tot2R < cut_totAcc[chIDR][Vov] ) continue;
        if( anEvent->energy2R < cut_energyMin[Form("%s_%s",channelR.c_str(),anEvent->stepLabel.c_str())] ) continue;
        if( anEvent->energy2R > cut_energyMax[Form("%s_%s",channelR.c_str(),anEvent->stepLabel.c_str())] ) continue;
      }
      
      int energyBin = h1_energyBins[Form("%s_%s",anEvent->ch1.c_str(),VovLabel.c_str())] -> Fill( anEvent->energy1 );
      if( energyBin < 1 || energyBin > nEnergyBinsEff[VovLabel] )
        continue;
      
      accept[label][entry] = true;
      
      (h1_energy_bin[Form("%s_%s",anEvent->ch1.c_str(),VovLabel.c_str())])[energyBin] -> Fill( anEvent->energy1 );
      std::string label12(Form("%s_energyBin%d",anEvent->label12.c_str(),energyBin));
      
      if( h1_totRatio[label12] == NULL )
      {
        h1_totRatio[label12] = new TH1F(Form("h1_totRatio_%s",label12.c_str()),"",1000,0.,5.);
        h1_energyRatio[label12] = new TH1F(Form("h1_energyRatio_%s",label12.c_str()),"",125,0.,5.);
        h1_deltaT_raw[label12] = new TH1F(Form("h1_deltaT_raw_%s",label12.c_str()),"",1250,-5000.,5000.);
      }
      
      h1_totRatio[label12] -> Fill( anEvent->tot2 / anEvent->tot1 );
      h1_energyRatio[label12] -> Fill( anEvent->energy2 / anEvent->energy1 );
      h1_deltaT_raw[label12] -> Fill( anEvent->time2-anEvent->time1 );
    }
    std::cout << std::endl;
  }
  
  
  
  //------------------
  //--- draw 2nd plots
  std::map<std::string,float> CTRMeans;
  std::map<std::string,float> CTRSigmas;
  
  for(auto stepLabel : stepLabels)
  {
    float Vov = map_Vovs[stepLabel];
    std::string VovLabel(Form("Vov%.1f",Vov));
    
    for(auto pair : pairsVec)
    {
      std::string ch1 = pair.first;
      std::string ch2 = pair.second;
      std::string label1(Form("%s_%s",ch1.c_str(),stepLabel.c_str()));
      std::string label2(Form("%s_%s",ch2.c_str(),stepLabel.c_str()));
      
      for(int ii = 1; ii <= nEnergyBinsEff[VovLabel]; ++ii)
      {
        std::string label12 = Form("%s-%s_%s_energyBin%d",ch1.c_str(),ch2.c_str(),stepLabel.c_str(),ii);
        
        std::cout << "find: " << label12 << std::endl;
        FindSmallestInterval(vals,h1_deltaT_raw[label12],0.68);
        float mean = vals[0];
        float min = vals[4];
        float max = vals[5];
        float delta = max-min;
        float sigma = 0.5*delta;
        float effSigma = sigma;
        CTRMeans[label12] = mean;
        CTRSigmas[label12] = effSigma;
      
        
        c = new TCanvas(Form("c_energyRatio_%s",label12.c_str()),Form("c_energyRatio_%s",label12.c_str()));
        gPad -> SetLogy();
        
        histo = h1_energyRatio[label12];
        histo -> SetTitle(Form(";%s energy / %s energy;entries",ch2.c_str(),ch1.c_str()));
        histo -> SetLineColor(kRed);
        histo -> Draw();
        
        float maxX = FindXMaximum(h1_energyRatio[label12],0.,5.);
        TF1* fitFunc = new TF1("fitFunc","gaus(0)",maxX-0.75*h1_energyRatio[label12]->GetRMS(),maxX+0.75*h1_energyRatio[label12]->GetRMS());
        h1_energyRatio[label12] -> Fit(fitFunc,"QRS+");
        
        histo -> Write();
        c -> Print(Form("%s/energyRatio/c_energyRatio__%s.png",plotDir.c_str(),label12.c_str()));
        c -> Print(Form("%s/energyRatio/c_energyRatio__%s.pdf",plotDir.c_str(),label12.c_str()));
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
    
    std::vector<std::string> tokens = GetTokens(label,'_');
    std::string stepLabel(Form("%s_%s",tokens[1].c_str(),tokens[2].c_str()));
    
    float Vov = map_Vovs[stepLabel];
    std::string VovLabel(Form("Vov%.1f",Vov));
    
    EventClass* anEvent = new EventClass();
    mapIt.second -> SetBranchAddress("event",&anEvent);
    
    int nEntries = mapIt.second->GetEntries();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      // if( entry%1000 == 0 ) std::cout << ">>> 3rd loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
      mapIt.second -> GetEntry(entry);
      
      if( !accept[label][entry] ) continue;
      
      int energyBin = h1_energyBins[Form("%s_%s",anEvent->ch1.c_str(),VovLabel.c_str())] -> Fill( anEvent->energy1 );
      std::string label12(Form("%s_energyBin%d",anEvent->label12.c_str(),energyBin));
      
      float timeLow = CTRMeans[label12] - 1.* CTRSigmas[label12];
      float timeHig = CTRMeans[label12] + 1.* CTRSigmas[label12];
      long long deltaT = anEvent->time2 - anEvent->time1;
      
      if( h1_deltaT[label12] == NULL )
      {
        h1_deltaT[label12] = new TH1F(Form("h1_deltaT_%s",label12.c_str()),"",1000,-5000,5000.);
      }
      
      h1_deltaT[label12] -> Fill( deltaT );

      if( !p1_deltaT_vs_energyRatio[label12] )
      {
        // float maxX = FindXMaximum(h1_energyRatio[label12],0.,5.);
        float maxX = h1_energyRatio[label12]->GetMean();
        TF1* fitFunc = new TF1("fitFunc","gaus(0)",maxX-2.*h1_energyRatio[label12]->GetRMS(),maxX+2.*h1_energyRatio[label12]->GetRMS());
        h1_energyRatio[label12] -> Fit(fitFunc,"QRS+");
        
        float xMin = maxX-5.00*h1_energyRatio[label12]->GetRMS();
        float xMax = maxX+5.00*h1_energyRatio[label12]->GetRMS();
        p1_deltaT_vs_energyRatio[label12] = new TProfile(Form("p1_deltaT_vs_energyRatio_%s",label12.c_str()),"",100,xMin,xMax);
      }
      
      if( ( deltaT > timeLow ) && ( deltaT < timeHig ) )
        p1_deltaT_vs_energyRatio[label12] -> Fill( anEvent->energy2/anEvent->energy1,anEvent->time2-anEvent->time1 );
    }
    std::cout << std::endl;
  }
  
  
  
  
  std::map<std::string,TF1*> fitFunc_energyCorr;
  
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
        
      for(int ii = 1; ii <= nEnergyBinsEff[VovLabel]; ++ii)
      {
        std::string label12 = Form("%s-%s_%s_energyBin%d",ch1.c_str(),ch2.c_str(),stepLabel.c_str(),ii);
        
        // TF1* fitFunc = (TF1*)( h1_energyRatio[label12]->GetFunction("fitFunc"))
        // float fitXMin = fitFunc->GetParameter(1) - 2.*fitFunc->GetParameter(2);
        // float fitXMax = fitFunc->GetParameter(1) + 2.*fitFunc->GetParameter(2);
        float fitXMin = h1_energyRatio[label12]->GetMean() - 3.*h1_energyRatio[label12]->GetRMS();
        float fitXMax = h1_energyRatio[label12]->GetMean() + 3.*h1_energyRatio[label12]->GetRMS();
        
        // fitFunc_energyCorr[label12] = new TF1(Form("fitFunc_energyCorr_%s",label12.c_str()),"pol1",fitXMin,fitXMax);
        fitFunc_energyCorr[label12] = new TF1(Form("fitFunc_energyCorr_%s",label12.c_str()),"pol4",fitXMin,fitXMax);
        fitFunc_energyCorr[label12] -> SetNpx(10000);
        fitFunc_energyCorr[label12] -> SetLineWidth(3);
        p1_deltaT_vs_energyRatio[label12] -> Fit(fitFunc_energyCorr[label12],"QRS+");
      }
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
        
        for(int ii = 1; ii <= nEnergyBinsEff[VovLabel]; ++ii)
        {
          std::string label12 = Form("%s-%s_%s_energyBin%d",ch1.c_str(),ch2.c_str(),stepLabel.c_str(),ii);
          
          
          //--------------------------------------------------------
          
          
          c = new TCanvas(Form("c_deltaT_vs_energyRatio_%s",label12.c_str()),Form("c_deltaT_vs_energyRatio_%s",label12.c_str()));
          
          prof = (TProfile*)( outFile->Get(Form("p1_deltaT_vs_energyRatio_%s",label12.c_str())) );
          prof -> SetTitle(Form(";%s energy / %s energy;#Deltat [ps]",ch2.c_str(),ch1.c_str()));
          prof -> Draw("");
          
          c -> Print(Form("%s/CTR/c_deltaT_vs_energyRatio__%s.png",plotDir.c_str(),label12.c_str()));
          c -> Print(Form("%s/CTR/c_deltaT_vs_energyRatio__%s.pdf",plotDir.c_str(),label12.c_str()));
          delete c;
        }
      }
    }
  }
  
  
  
  
  //------------------------
  //--- 4th loop over events
  for(auto mapIt : trees)
  {
    std::string label = mapIt.first;
    std::cout << ">>> 4th loop: " << label << std::endl;
    
    std::vector<std::string> tokens = GetTokens(label,'_');
    std::string stepLabel(Form("%s_%s",tokens[1].c_str(),tokens[2].c_str()));
    
    float Vov = map_Vovs[stepLabel];
    std::string VovLabel(Form("Vov%.1f",Vov));
    
    EventClass* anEvent = new EventClass();
    mapIt.second -> SetBranchAddress("event",&anEvent);
    
    int nEntries = mapIt.second->GetEntries();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      // if( entry%1000 == 0 ) std::cout << ">>> 4th loop (" << label << "): reading entry " << entry << " / " << nEntries << "\r" << std::flush;
      mapIt.second -> GetEntry(entry);
      
      if( !accept[label][entry] ) continue;
      
      int energyBin = h1_energyBins[Form("%s_%s",anEvent->ch1.c_str(),VovLabel.c_str())] -> Fill( anEvent->energy1 );
      std::string label12(Form("%s_energyBin%d",anEvent->label12.c_str(),energyBin));
      
      float timeLow = CTRMeans[label12] - 2.* CTRSigmas[label12];
      float timeHig = CTRMeans[label12] + 2.* CTRSigmas[label12];
      
      long long deltaT = anEvent->time2 - anEvent->time1;
      
      float energyCorr = fitFunc_energyCorr[label12]->Eval(anEvent->energy2/anEvent->energy1) -
                         fitFunc_energyCorr[label12]->Eval(h1_energyRatio[label12]->GetMean());
      
      if( h1_deltaT_energyCorr[label12] == NULL )
      {
        h1_deltaT_energyCorr[label12] = new TH1F(Form("h1_deltaT_energyCorr_%s",label12.c_str()),"",1000,-5000.,5000.);
      }
      
      h1_deltaT_energyCorr[label12] -> Fill( deltaT - energyCorr );
      
      if( !p1_deltaT_energyCorr_vs_pos[label12] )
      {
        p1_deltaT_energyCorr_vs_pos[label12] =  new TProfile(Form("p1_deltaT_energyCorr_vs_pos_%s",label12.c_str()),"",200,0.,50.);
      }
      
      if( ( (deltaT-energyCorr) > timeLow ) && ( (deltaT-energyCorr) < timeHig ) )
      {
        if( (anEvent->isHorizontal1!=-1 && anEvent->isHorizontal1==1) || (anEvent->isHorizontal2!=-1 && anEvent->isHorizontal2==1) ) 
          p1_deltaT_energyCorr_vs_pos[label12] -> Fill( anEvent->x,deltaT-energyCorr );
        if( (anEvent->isHorizontal1!=-1 && anEvent->isHorizontal1==0) || (anEvent->isHorizontal2!=-1 && anEvent->isHorizontal2==0) ) 
          p1_deltaT_energyCorr_vs_pos[label12] -> Fill( anEvent->y,deltaT-energyCorr );
      }
    }
    std::cout << std::endl;
  }




  //------------------
  //--- draw 4th plots
  TH1F* hamp;
  TH1F* hampL;
  TH1F* hampR;
  TF1 *fLandau1;
  TF1 *fLandau2;
  float norm;
  
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

  std::map<std::string,TGraphErrors*> g_tRes_gaus_vs_dVdt;
  std::map<std::string,TGraphErrors*> g_tRes_energyCorr_gaus_vs_dVdt;
  
  std::map<std::string,TGraphErrors*> g_slewRate_vs_th;
  std::map<std::string,TGraphErrors*> g_slewRateNormalized_vs_th;
  std::map<std::string,TGraphErrors*> g_dVdt_vs_th;
  std::map<std::string,TGraphErrors*> g_tRes_noise_vs_th;
  std::map<std::string,TGraphErrors*> g_tRes_noise_corr_vs_th;

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
        
        for(int ii = 1; ii <= nEnergyBinsEff[VovLabel]; ++ii)
        {
          std::string label12 = Form("%s-%s_%s_energyBin%d",ch1.c_str(),ch2.c_str(),stepLabel.c_str(),ii);
          
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
          
          // -- energy corr deltaT
          histo = h1_deltaT_energyCorr[label12];
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
          histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
          // float fitXMin2 = CTRMeans[label12] - 5.*CTRSigmas[label12];
          // float fitXMax2 = CTRMeans[label12] + 5.*CTRSigmas[label12];
          // TF1* fitFunc2 = new TF1(Form("fitFunc2_energyCorr_%s",label12.c_str()),"gaus(0)+gaus(3)+gaus(6)",fitXMin2,fitXMax2);
          // fitFunc2 -> SetParameter(0,fitFunc->GetParameter(0));
          // fitFunc2 -> SetParameter(1,fitFunc->GetParameter(1));
          // fitFunc2 -> SetParameter(2,fitFunc->GetParameter(2));
          // fitFunc2 -> SetParameter(3,fitFunc->GetParameter(0)/3.);
          // //fitFunc2 -> FixParameter(4,fitFunc->GetParameter(1)-330.);
          // fitFunc2 -> FixParameter(4,fitFunc->GetParameter(1)-350.);
          // fitFunc2 -> SetParameter(5,fitFunc->GetParameter(2));
          // fitFunc2 -> SetParameter(6,fitFunc->GetParameter(0)/3.);
          // //fitFunc2 -> FixParameter(7,fitFunc->GetParameter(1)+360.);
          // fitFunc2 -> FixParameter(7,fitFunc->GetParameter(1)+350.);
          // fitFunc2 -> SetParameter(8,fitFunc->GetParameter(2));
          float fitXMin2 = fitFunc->GetParameter(1)-fabs(2.0*fitFunc->GetParameter(2));
          float fitXMax2 = fitFunc->GetParameter(1)+fabs(2.0*fitFunc->GetParameter(2));
          TF1* fitFunc2 = new TF1(Form("fitFunc2_energyCorr_%s",label12.c_str()),"gaus(0)",fitXMin2,fitXMax2);
          fitFunc2 -> SetParameter(0,fitFunc->GetParameter(0));
          fitFunc2 -> SetParameter(1,fitFunc->GetParameter(1));
          fitFunc2 -> SetParameter(2,fitFunc->GetParameter(2));
          histo -> Fit(fitFunc2,"QRSL+");
          
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
          
          latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corr.}^{eff} = %.0f ps}{#sigma_{corr.}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc2->GetParameter(2))));
          latex -> SetNDC();
          latex -> SetTextFont(42);
          latex -> SetTextSize(0.04);
          latex -> SetTextColor(kBlue);
          latex -> Draw("same");
          
          /*
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
          
          g_tRes_energyCorr_effSigma_vs_th[label_vs_th] -> SetPoint(g_tRes_energyCorr_effSigma_vs_th[label_vs_th]->GetN(),th,effSigma*corr);
          g_tRes_energyCorr_effSigma_vs_th[label_vs_th] -> SetPointError(g_tRes_energyCorr_effSigma_vs_th[label_vs_th]->GetN()-1,0.,5.);
          g_tRes_energyCorr_gaus_vs_th[label_vs_th] -> SetPoint(g_tRes_energyCorr_gaus_vs_th[label_vs_th]->GetN(),th,fabs(fitFunc2->GetParameter(2))*corr);
          g_tRes_energyCorr_gaus_vs_th[label_vs_th] -> SetPointError(g_tRes_energyCorr_gaus_vs_th[label_vs_th]->GetN()-1,0.,fitFunc2->GetParError(2)*corr);
          
          g_tRes_energyCorr_effSigma_vs_Vov[label_vs_Vov] -> SetPoint(g_tRes_energyCorr_effSigma_vs_Vov[label_vs_Vov]->GetN(),Vov,effSigma*corr);
          g_tRes_energyCorr_effSigma_vs_Vov[label_vs_Vov] -> SetPointError(g_tRes_energyCorr_effSigma_vs_Vov[label_vs_Vov]->GetN()-1,0.,5.);
          g_tRes_energyCorr_gaus_vs_Vov[label_vs_Vov] -> SetPoint(g_tRes_energyCorr_gaus_vs_Vov[label_vs_Vov]->GetN(),Vov,fabs(fitFunc2->GetParameter(2))*corr);
          g_tRes_energyCorr_gaus_vs_Vov[label_vs_Vov] -> SetPointError(g_tRes_energyCorr_gaus_vs_Vov[label_vs_Vov]->GetN()-1,0.,fitFunc2->GetParError(2)*corr);
          
          if( ( tRes_energyCorr_gaus_bestTh[label_bestTh_vs_Vov][Vov] != 0. && fabs(fitFunc2->GetParameter(2))*corr < tRes_energyCorr_gaus_bestTh[label_bestTh_vs_Vov][Vov] ) ||
              tRes_energyCorr_gaus_bestTh[label_bestTh_vs_Vov][Vov] == 0 )
          {
            tRes_energyCorr_gaus_bestTh[label_bestTh_vs_Vov][Vov] = fabs(fitFunc2->GetParameter(2))*corr;
            tResErr_energyCorr_gaus_bestTh[label_bestTh_vs_Vov][Vov] = fitFunc2->GetParError(2)*corr;
          }
          if( ( tRes_energyCorr_effSigma_bestTh[label_bestTh_vs_Vov][Vov] != 0. && effSigma*corr < tRes_energyCorr_effSigma_bestTh[label_bestTh_vs_Vov][Vov] ) ||
              tRes_energyCorr_effSigma_bestTh[label_bestTh_vs_Vov][Vov] == 0. )
          {
            tRes_energyCorr_effSigma_bestTh[label_bestTh_vs_Vov][Vov] = effSigma*corr;
          }
          */
          
          
          // -- raw delta T
          c->cd();
          histo = h1_deltaT[label12];
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
          histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
          // fitFunc2 = new TF1(Form("fitFunc2_energyCorr_%s",label12.c_str()),"gaus(0)+gaus(3)+gaus(6)",fitXMin2,fitXMax2);
          // fitFunc2 -> SetParameter(0,fitFunc->GetParameter(0));
          // fitFunc2 -> SetParameter(1,fitFunc->GetParameter(1));
          // fitFunc2 -> SetParameter(2,fitFunc->GetParameter(2));
          // fitFunc2 -> SetParameter(3,fitFunc->GetParameter(0)/3.);
          // fitFunc2 -> FixParameter(4,fitFunc->GetParameter(1)-350.);
          // fitFunc2 -> SetParameter(5,fitFunc->GetParameter(2));
          // fitFunc2 -> SetParameter(6,fitFunc->GetParameter(0)/3.);
          // fitFunc2 -> FixParameter(7,fitFunc->GetParameter(1)+350.);
          // fitFunc2 -> SetParameter(8,fitFunc->GetParameter(2));
          fitXMin2 = fitFunc->GetParameter(1)-fabs(2.*fitFunc->GetParameter(2));
          fitXMax2 = fitFunc->GetParameter(1)+fabs(2.*fitFunc->GetParameter(2));
          fitFunc2 = new TF1(Form("fitFunc2_energyCorr_%s",label12.c_str()),"gaus(0)",fitXMin2,fitXMax2);
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
          
          /*
          g_tRes_effSigma_vs_th[label_vs_th] -> SetPoint(g_tRes_effSigma_vs_th[label_vs_th]->GetN(),th,effSigma*corr);
          g_tRes_effSigma_vs_th[label_vs_th] -> SetPointError(g_tRes_effSigma_vs_th[label_vs_th]->GetN()-1,0.,5.);
          g_tRes_gaus_vs_th[label_vs_th] -> SetPoint(g_tRes_gaus_vs_th[label_vs_th]->GetN(),th,fabs(fitFunc2->GetParameter(2))*corr);
          g_tRes_gaus_vs_th[label_vs_th] -> SetPointError(g_tRes_gaus_vs_th[label_vs_th]->GetN()-1,0.,fitFunc2->GetParError(2)*corr);
          
          g_tRes_effSigma_vs_Vov[label_vs_Vov] -> SetPoint(g_tRes_effSigma_vs_Vov[label_vs_Vov]->GetN(),Vov,effSigma*corr);
          g_tRes_effSigma_vs_Vov[label_vs_Vov] -> SetPointError(g_tRes_effSigma_vs_Vov[label_vs_Vov]->GetN()-1,0.,5.);
          g_tRes_gaus_vs_Vov[label_vs_Vov] -> SetPoint(g_tRes_gaus_vs_Vov[label_vs_Vov]->GetN(),Vov,fabs(fitFunc2->GetParameter(2))*corr);
          g_tRes_gaus_vs_Vov[label_vs_Vov] -> SetPointError(g_tRes_gaus_vs_Vov[label_vs_Vov]->GetN()-1,0.,fitFunc2->GetParError(2)*corr);
          
          if( ( tRes_gaus_bestTh[label_bestTh_vs_Vov][Vov] != 0. && fabs(fitFunc2->GetParameter(2))*corr < tRes_gaus_bestTh[label_bestTh_vs_Vov][Vov] ) ||
              tRes_gaus_bestTh[label_bestTh_vs_Vov][Vov] == 0. )
          {
            tRes_gaus_bestTh[label_bestTh_vs_Vov][Vov] = fabs(fitFunc2->GetParameter(2))*corr;
            tResErr_gaus_bestTh[label_bestTh_vs_Vov][Vov] = fitFunc2->GetParError(2)*corr;
          }
          if( ( tRes_effSigma_bestTh[label_bestTh_vs_Vov][Vov] != 0. && effSigma*corr < tRes_effSigma_bestTh[label_bestTh_vs_Vov][Vov] ) ||
              tRes_effSigma_bestTh[label_bestTh_vs_Vov][Vov] == 0. )
          {
            tRes_effSigma_bestTh[label_bestTh_vs_Vov][Vov] = effSigma*corr;
          }
          */
          
          
          c -> Print(Form("%s/CTR_energyCorr/c_deltaT_energyCorr__%s.png",plotDir.c_str(),label12.c_str()));
          c -> Print(Form("%s/CTR_energyCorr/c_deltaT_energyCorr__%s.pdf",plotDir.c_str(),label12.c_str()));
          delete c;
          
          ++pairsIt;
        }
      }
    }
  }
  
  
  /*
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
      int isBar1 = opts.GetOpt<int>(Form("%s.isBar",ch1.c_str()));
      int isBar2 = opts.GetOpt<int>(Form("%s.isBar",ch2.c_str()));

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
        g_effSigma -> SetMarkerStyle(25);
        // g_effSigma -> Draw("PL,same");
        
        g_gaus -> SetLineColor(1+iter);
        g_gaus -> SetMarkerColor(1+iter);
        g_gaus -> SetMarkerStyle(25);
        g_gaus -> Draw("PL,same");
        
        g_energyCorr_effSigma -> SetLineColor(1+iter);
        g_energyCorr_effSigma -> SetMarkerColor(1+iter);
        g_energyCorr_effSigma -> SetMarkerStyle(20);
        // g_energyCorr_effSigma -> Draw("PL,same");
        
        g_energyCorr_gaus -> SetLineColor(1+iter);
        g_energyCorr_gaus -> SetMarkerColor(1+iter);
        g_energyCorr_gaus -> SetMarkerStyle(20);
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


      c2 = new TCanvas(Form("c_slewRateNormalized_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_slewRateNormalized_%s-%s",ch1.c_str(),ch2.c_str()));
      TH1F *hPad2 = (TH1F*)( gPad->DrawFrame(-0.5,0.,2.5,15.) );
      hPad2 -> SetTitle(";#LT t_{diff} #GT [ns]; normalized threshold");
      hPad2 -> Draw();
      gPad -> SetGridy();

      c3 = new TCanvas(Form("c_dVdt_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_dVdt_%s-%s",ch1.c_str(),ch2.c_str()));
      TH1F *hPad3 = (TH1F*)( gPad->DrawFrame(-0.5,0.,64,400.) );
      hPad3 -> SetTitle("; threshold [DAC]; dV/dt [a.u.]");
      hPad3 -> Draw();
      gPad -> SetGridy();
      
      c4 = new TCanvas(Form("c_tRes_noise_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_tRes_noise_%s-%s",ch1.c_str(),ch2.c_str()));
      TH1F *hPad4 = (TH1F*)( gPad->DrawFrame(-0.5,0.,64,150.) );
      hPad4 -> SetTitle("; threshold [DAC]; #sigma_{noise} (ps)");
      hPad4 -> Draw();
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
	  //std::cout << "x: " << x << "   y: " << y << "   y0: " << y0 << "   val: " << fabs(y-y0)/1000. << std::endl;
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
        
	c2->cd();
        g_slewRateNormalized_final -> SetLineColor(1+iter);
        g_slewRateNormalized_final -> SetMarkerColor(1+iter);
        g_slewRateNormalized_final -> SetMarkerStyle(20);
        g_slewRateNormalized_final -> Draw("PL,same");

	latex -> Draw("same");
        

	// derivative of the pulse shape vs threshold
	g_dVdt_vs_th[label] = new TGraphErrors();
	g_tRes_noise_vs_th[label] = new TGraphErrors();
	float sigmaV;
	if (isBar1 == 0) {
	  int chID = opts.GetOpt<int>(Form("%s.chID",ch1.c_str()));
	  sigmaV = noise[chID];
	}
	else{
	  std::string channelL = opts.GetOpt<std::string>(Form("%s.channelL",ch1.c_str()));
	  std::string channelR = opts.GetOpt<std::string>(Form("%s.channelR",ch1.c_str()));
	  int chIDL = opts.GetOpt<int>(Form("%s.chID",channelL.c_str()));
	  int chIDR = opts.GetOpt<int>(Form("%s.chID",channelR.c_str()));
	  sigmaV = 0.5*(noise[chIDL] + noise[chIDR]);
	}
	for(int point = 0; point < g_slewRate->GetN(); ++point)
	  {
	    double x,y;
	    g_slewRate -> GetPoint(point,x,y);
            float delta = 50.; // ps
            float dVdt = ( g_slewRate->Eval(x-delta) - g_slewRate->Eval(x+delta))/(2*delta);
            g_dVdt_vs_th[label] -> SetPoint(point, y, dVdt*1000.);
            g_tRes_noise_vs_th[label] -> SetPoint(point, y, sigmaV/dVdt);
	    //std::cout << "x: " << x << "   y: " << y << "   y0: " << y0 << "   val: " << fabs(y-y0)/1000. << std::endl; 
	  }
	
	c3->cd();
        g_dVdt_vs_th[label] -> SetLineColor(1+iter);
        g_dVdt_vs_th[label] -> SetMarkerColor(1+iter);
        g_dVdt_vs_th[label] -> SetMarkerStyle(20);
        g_dVdt_vs_th[label] -> Draw("PL,same");

        latex -> Draw("same");


        c4->cd();
        g_tRes_noise_vs_th[label] -> SetLineColor(1+iter);                                                                                                                      
        g_tRes_noise_vs_th[label] -> SetMarkerColor(1+iter);                                                                                                                    
        //g_tRes_noise_vs_th[label] -> SetLineColor(color[vov]);
        //g_tRes_noise_vs_th[label] -> SetMarkerColor(color[vov]);
        g_tRes_noise_vs_th[label] -> SetMarkerStyle(20);
        g_tRes_noise_vs_th[label] -> Draw("PL,same");

        latex -> Draw("same");        
	++iter;
      }
      
      c -> Print(Form("%s/c_slewRate__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c -> Print(Form("%s/c_slewRate__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c2 -> Print(Form("%s/c_slewRateNormalized__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c2 -> Print(Form("%s/c_slewRateNormalized__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c3 -> Print(Form("%s/c_dVdt__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c3 -> Print(Form("%s/c_dVdt__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c4 -> Print(Form("%s/c_tRes_noise_%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
      c4 -> Print(Form("%s/c_tRes_noise_%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));

      delete c;
      delete c2;
      delete c3;
      delete c4;
      
      ++pairsIt;
    }
    
    
    //--------------------------------------------------------
    pairsIt = 0;
    for(auto pair : pairsVec)
      {
	std::string ch1 = pair.first;
	std::string ch2 = pair.second;

	c = new TCanvas(Form("c_tRes_vs_dVdt_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_tRes_vs_dVdt_%s-%s",ch1.c_str(),ch2.c_str()));
	// gPad -> SetLogy();                                                                                                     

	TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,tResMin,400,tResMax) );
	if( pairsMode.at(pairsIt) == 0 )
	  hPad -> SetTitle("; dV/dt [DAC/ns];#sigma_{t_{diff}} [ps]");
	if( pairsMode.at(pairsIt) == 1 )
	  hPad -> SetTitle(";dV/dt [DAC/ns];#sigma_{t_{diff}} / #sqrt{2} [ps]");
	if( pairsMode.at(pairsIt) == 2 )
	  hPad -> SetTitle(";dV/dt [DAC/ns];#sigma_{t_{diff}} / 2 [ps]");
	hPad -> Draw();
	gPad -> SetGridy();

	int iter = 0;
	for(auto mapIt : VovLabels)
	  {
	    std::string label(Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),mapIt.first.c_str()));
	    g_tRes_gaus_vs_dVdt[label] = new TGraphErrors();
	    g_tRes_energyCorr_gaus_vs_dVdt[label] = new TGraphErrors();
	    for(int point = 0; point < g_tRes_gaus_vs_th[label]->GetN(); ++point)
	      {
		double x,y;
		g_tRes_gaus_vs_th[label] -> GetPoint(point,x,y);
		g_tRes_gaus_vs_dVdt[label]->SetPoint(point,g_dVdt_vs_th[label]->Eval(x),y);
		g_tRes_gaus_vs_dVdt[label]->SetPointError(point,0,g_tRes_gaus_vs_th[label]->GetEY()[point]);
		g_tRes_energyCorr_gaus_vs_th[label] -> GetPoint(point,x,y);
		g_tRes_energyCorr_gaus_vs_dVdt[label]->SetPoint(point,g_dVdt_vs_th[label]->Eval(x),y);
		g_tRes_energyCorr_gaus_vs_dVdt[label]->SetPointError(point,0,g_tRes_energyCorr_gaus_vs_th[label]->GetEY()[point]);
	      }
	  
	
	    g_tRes_gaus_vs_dVdt[label] -> SetLineColor(1+iter);
	    g_tRes_gaus_vs_dVdt[label] -> SetMarkerColor(1+iter);
	    g_tRes_gaus_vs_dVdt[label] -> SetMarkerStyle(25);
	    //g_tRes_gaus_vs_dVdt[label] -> Draw("PL,same");
	    
	    g_tRes_energyCorr_gaus_vs_dVdt[label] -> SetLineColor(1+iter);
	    g_tRes_energyCorr_gaus_vs_dVdt[label] -> SetMarkerColor(1+iter);
	    g_tRes_energyCorr_gaus_vs_dVdt[label] -> SetMarkerStyle(20);
	    g_tRes_energyCorr_gaus_vs_dVdt[label] -> Draw("PL,same");
	    
	    latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
	    latex -> SetNDC();
	    latex -> SetTextFont(42);
	    latex -> SetTextSize(0.04);
	    latex -> SetTextColor(kBlack+iter);
	    latex -> Draw("same");
	    
	    ++iter;
	  }

	c -> Print(Form("%s/c_tRes_vs_dVdt__%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
	c -> Print(Form("%s/c_tRes_vs_dVdt__%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
	delete c;

	++pairsIt;      
      }


    
    //--------------------------------------------------------
    // correlated noise term
    pairsIt = 0;
    for(auto pair : pairsVec)
      {
	std::string ch1 = pair.first;
	std::string ch2 = pair.second;
        int isBar1 = opts.GetOpt<int>(Form("%s.isBar",ch1.c_str()));
        int isBar2 = opts.GetOpt<int>(Form("%s.isBar",ch2.c_str()));

        // consider only bar vs pixel as ref
        if ( !(  ( isBar1 && ch2.find("pixel") != std::string::npos ) ||
                 ( isBar2 && ch1.find("pixel") != std::string::npos )  ) ) continue;

        c = new TCanvas(Form("c_tRes_corr_vs_th_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_tRes_corr_vs_th_%s-%s",ch1.c_str(),ch2.c_str()));
        // gPad -> SetLogy();

        TH1F* hPad = (TH1F*)( gPad->DrawFrame(-0.5,tResMin,64,tResMax) );
        if( pairsMode.at(pairsIt) == 0 )
          hPad -> SetTitle("; threshold [DAC];#sigma_{t_{diff}} [ps]");
        if( pairsMode.at(pairsIt) == 1 )
          hPad -> SetTitle(";threshold [DAC];#sigma_{t_{diff}} / #sqrt{2} [ps]");
        if( pairsMode.at(pairsIt) == 2 )
          hPad -> SetTitle(";threhsold [DAC];#sigma_{t_{diff}} / 2 [ps]");
        hPad -> Draw();
        gPad -> SetGridy();

	int iter = 0;
        for(auto mapIt : VovLabels)
          {
	    std::string vov = mapIt.first.c_str();
	    std::string label(Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),mapIt.first.c_str()));
	    std::string label2(Form("%sL-%sR_%s",ch1.c_str(),ch1.c_str(),mapIt.first.c_str()));
	    std::string label3(Form("pixel3x3-pixel2x2_%s",mapIt.first.c_str()));
            if (isBar2){
              label2 = Form("%sL-%sR_%s",ch2.c_str(),ch2.c_str(),mapIt.first.c_str());
            }
	    std::cout << label.c_str() <<  "  " << label2.c_str()<< "  " <<  label3.c_str() << std::endl;
            g_tRes_noise_corr_vs_th[label] = new TGraphErrors();
            for(int point = 0; point < g_tRes_energyCorr_gaus_vs_th[label]->GetN(); ++point)
              {
                double y1, y2, y3, th;
                g_tRes_energyCorr_gaus_vs_th[label] -> GetPoint(point,th,y1);// bar -pixel                                                                                        
                g_tRes_energyCorr_gaus_vs_th[label2] -> GetPoint(point,th,y2); // diff                                                                                            
                g_tRes_energyCorr_gaus_vs_th[label3] -> GetPoint(point,th,y3);// pixel                                                                                            
                double y = sqrt(y1*y1 - (y2*y2/4) - (y3*y3/2));
                g_tRes_noise_corr_vs_th[label]->SetPoint(point,th,y);
              }

            //g_tRes_noise_corr_vs_th[label] -> SetLineColor(color[vov]);
            //g_tRes_noise_corr_vs_th[label] -> SetMarkerColor(color[vov]);
	    g_tRes_noise_corr_vs_th[label] -> SetLineColor(iter+1);
	    g_tRes_noise_corr_vs_th[label] -> SetMarkerColor(iter+1);
            g_tRes_noise_corr_vs_th[label] -> SetMarkerStyle(20);
            g_tRes_noise_corr_vs_th[label] -> Draw("PL,same");

            latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            //latex -> SetTextColor(color[vov]);
	    latex -> SetTextColor(iter+1);
            latex -> Draw("same");

            ++iter;
          }

        if (isBar1){
          c -> Print(Form("%s/c_tRes_noise_corr_%sL-%sR.png",plotDir.c_str(),ch1.c_str(),ch1.c_str()));
          c -> Print(Form("%s/c_tRes_noise_corr_%sL-%sR.pdf",plotDir.c_str(),ch1.c_str(),ch1.c_str()));
        }
        if (isBar2){
          c -> Print(Form("%s/c_tRes_noise_corr_%sL-%sR.png",plotDir.c_str(),ch2.c_str(),ch2.c_str()));
          c -> Print(Form("%s/c_tRes_noise_corr_%sL-%sR.pdf",plotDir.c_str(),ch2.c_str(),ch2.c_str()));
        }

        delete c;

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
        g_effSigma -> SetMarkerStyle(25);
        // g_effSigma -> Draw("PL,same");
        
        g_gaus -> SetLineColor(1+iter);
        g_gaus -> SetMarkerColor(1+iter);
        g_gaus -> SetMarkerStyle(25);
        g_gaus -> Draw("PL,same");
        
        g_energyCorr_effSigma -> SetLineColor(1+iter);
        g_energyCorr_effSigma -> SetMarkerColor(1+iter);
        g_energyCorr_effSigma -> SetMarkerStyle(20);
        // g_energyCorr_effSigma -> Draw("PL,same");
        
        g_energyCorr_gaus -> SetLineColor(1+iter);
        g_energyCorr_gaus -> SetMarkerColor(1+iter);
        g_energyCorr_gaus -> SetMarkerStyle(20);
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
          */
  
  
  
  outFile -> cd();
  
  // for(auto mapIt: g_tot_vs_th)  mapIt.second -> Write(Form("g_tot_vs_th_%s", mapIt.first.c_str()));
  // for(auto mapIt: g_tot_vs_Vov) mapIt.second -> Write(Form("g_tot_vs_Vov_%s",mapIt.first.c_str()));
  // for(auto mapIt: g_energy_vs_th)  mapIt.second -> Write(Form("g_energy_vs_th_%s", mapIt.first.c_str()));
  // for(auto mapIt: g_energy_vs_Vov) mapIt.second -> Write(Form("g_energy_vs_Vov_%s",mapIt.first.c_str()));
  
  // for(auto mapIt: g_tRes_effSigma_vs_th)         mapIt.second -> Write(Form("g_tRes_effSigma_vs_th_%s",         mapIt.first.c_str()));
  // for(auto mapIt: g_tRes_gaus_vs_th)             mapIt.second -> Write(Form("g_tRes_gaus_vs_th_%s",             mapIt.first.c_str()));
  // for(auto mapIt: g_tRes_effSigma_vs_Vov)        mapIt.second -> Write(Form("g_tRes_effSigma_vs_Vov_%s",        mapIt.first.c_str()));
  // for(auto mapIt: g_tRes_gaus_vs_Vov)            mapIt.second -> Write(Form("g_tRes_gaus_vs_Vov_%s",            mapIt.first.c_str()));
  // for(auto mapIt: g_tRes_effSigma_bestTh_vs_Vov) mapIt.second -> Write(Form("g_tRes_effSigma_bestTh_vs_Vov_%s", mapIt.first.c_str()));
  // for(auto mapIt: g_tRes_gaus_bestTh_vs_Vov)     mapIt.second -> Write(Form("g_tRes_gaus_bestTh_vs_Vov_%s",     mapIt.first.c_str()));
  
  // for(auto mapIt: g_tRes_energyCorr_effSigma_vs_th)         mapIt.second -> Write(Form("g_tRes_energyCorr_effSigma_vs_th_%s",         mapIt.first.c_str()));
  // for(auto mapIt: g_tRes_energyCorr_gaus_vs_th)             mapIt.second -> Write(Form("g_tRes_energyCorr_gaus_vs_th_%s",             mapIt.first.c_str()));
  // for(auto mapIt: g_tRes_energyCorr_effSigma_vs_Vov)        mapIt.second -> Write(Form("g_tRes_energyCorr_effSigma_vs_Vov_%s",        mapIt.first.c_str()));
  // for(auto mapIt: g_tRes_energyCorr_gaus_vs_Vov)            mapIt.second -> Write(Form("g_tRes_energyCorr_gaus_vs_Vov_%s",            mapIt.first.c_str()));
  // for(auto mapIt: g_tRes_energyCorr_effSigma_bestTh_vs_Vov) mapIt.second -> Write(Form("g_tRes_energyCorr_effSigma_bestTh_vs_Vov_%s", mapIt.first.c_str()));
  // for(auto mapIt: g_tRes_energyCorr_gaus_bestTh_vs_Vov)     mapIt.second -> Write(Form("g_tRes_energyCorr_gaus_bestTh_vs_Vov_%s",     mapIt.first.c_str()));
  
  // for(auto mapIt: g_tRes_gaus_vs_dVdt)            mapIt.second -> Write(Form("g_tRes_gaus_vs_dVdt_%s",           mapIt.first.c_str()));
  // for(auto mapIt: g_tRes_energyCorr_gaus_vs_dVdt) mapIt.second -> Write(Form("g_tRes_energyCorr_gaus_vs_dVdt_%s",mapIt.first.c_str()));
  
  // for(auto mapIt: g_slewRate_vs_th)           mapIt.second -> Write(Form("g_slewRate_vs_th_%s",           mapIt.first.c_str()));
  // for(auto mapIt: g_slewRateNormalized_vs_th) mapIt.second -> Write(Form("g_slewRateNormalized_vs_th_%s", mapIt.first.c_str()));
  // for(auto mapIt:  g_dVdt_vs_th)              mapIt.second -> Write(Form(" g_dVdt_vs_th_%s",              mapIt.first.c_str()));

  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
