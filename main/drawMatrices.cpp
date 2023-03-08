#include "interface/AnalysisUtils.h"
#include "interface/FitUtils.h"
#include "interface/SetTDRStyle.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <vector>
#include <time.h>
#include <stdio.h>
#include <sys/stat.h>
#include <algorithm>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"



struct channel
{
  int barID;
  int lrID;
};



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
  
  
  
  //--- open files and make the tree chain
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  std::vector<std::string> fileBaseNames = opts.GetOpt<std::vector<std::string> >("Input.fileBaseNames");
  int maxEntries = opts.GetOpt<int>("Input.maxEntries");
  TChain* tree = new TChain("data","data");
  
  time_t timesec;
  for(auto fileBaseName : fileBaseNames)
  {
    //std::string fileName = Form("%s/%s_ped_e.root",inputDir.c_str(),fileBaseName.c_str());
    std::string fileName = Form("%s/%s_e.root",inputDir.c_str(),fileBaseName.c_str());
    std::cout << ">>> Adding flle " << fileName << std::endl;
    tree -> Add(fileName.c_str());
    
    struct stat t_stat;
    //stat(Form("/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Feb2020/TOFHIR/RawData/%s.rawf",fileBaseName.c_str()), &t_stat);
    stat(Form("/storage/TOFHIR2/raw/%s.rawf",fileBaseName.c_str()), &t_stat);
    struct tm * timeinfo = localtime(&t_stat.st_mtime);
    timesec = mktime(timeinfo);
    std::cout << "Time and date of raw file of run " << fileBaseName << ": " << asctime(timeinfo);
  }
  
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  system(Form("mkdir -p %s",plotDir.c_str()));
  
  
  
  std::vector<unsigned int> channelMapping = opts.GetOpt<std::vector<unsigned int> >("Channels.channelMapping");
  

  
  float qfineMin = opts.GetOpt<float>("Cuts.qfineMin");
  int nEnergyBins = opts.GetOpt<int>("Cuts.nEnergyBins");
  float energyMin = opts.GetOpt<float>("Cuts.energyMin");
  float energyMax = opts.GetOpt<float>("Cuts.energyMax");
  
  
  
  //--- define branches
  float step1;
  int channelIdx[256];
  std::vector<unsigned short> *qfine = 0;
  std::vector<float> *tot = 0;
  std::vector<float> *energy = 0;

  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("step1",     1); tree -> SetBranchAddress("step1",    &step1);
  tree -> SetBranchStatus("channelIdx",1); tree -> SetBranchAddress("channelIdx",channelIdx);
  tree -> SetBranchStatus("qfine",  1); tree -> SetBranchAddress("qfine",   &qfine);  
  tree -> SetBranchStatus("tot",    1); tree -> SetBranchAddress("tot",       &tot);
  tree -> SetBranchStatus("energy", 1); tree -> SetBranchAddress("energy", &energy);
  
  
  
  //--- define histograms
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileName");
  TFile* outFile = TFile::Open(Form("%s",outFileName.c_str()),"RECREATE");
  outFile -> cd();
  
  std::map<std::string,std::map<int,std::map<int,TH1F*> > > h1_energy;
  std::vector<std::string> histoLabels;
  histoLabels.push_back("L");
  histoLabels.push_back("R");
  histoLabels.push_back("LR");
  
  
  
  //--- loop over events
  std::map<float,int> VovMap;
  
  int nEntries = tree->GetEntries();
  if( maxEntries > nEntries ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    tree -> GetEntry(entry);
    if( entry%100000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    
    VovMap[step1] += 1;
    
    
    float energyLRMax = -1;
    int barIDMax = -1;
    for(int barIt = 0; barIt < 16; ++barIt)
    {
      if (channelIdx[channelMapping[barIt*2+0]] < 0 || channelIdx[channelMapping[barIt*2+1]] < 0 ) continue;

      float energyL = (*qfine)[channelIdx[channelMapping[barIt*2+0]]] > qfineMin ? (*energy)[channelIdx[channelMapping[barIt*2+0]]] : 0.;
      float energyR = (*qfine)[channelIdx[channelMapping[barIt*2+1]]] > qfineMin ? (*energy)[channelIdx[channelMapping[barIt*2+1]]] : 0.;
      float energyLR = 0.5*(energyL+energyR);
      
      if( energyLR > energyLRMax )
      {
        energyLRMax = energyLR;
        barIDMax = barIt;
      }
    }
    
    if (channelIdx[channelMapping[barIDMax*2+0]] < 0 || channelIdx[channelMapping[barIDMax*2+1]] < 0 ) continue;

    if( (*qfine)[channelIdx[channelMapping[barIDMax*2+0]]] > qfineMin || (*qfine)[channelIdx[channelMapping[barIDMax*2+1]]] > qfineMin )
    {
      if( h1_energy["LR"][step1][barIDMax] == NULL )
      {
        std::string label_L(Form("h1_energy_%s_bar%02d_L_Vov%.1f","prova",barIDMax,step1));
        std::string label_R(Form("h1_energy_%s_bar%02d_R_Vov%.1f","prova",barIDMax,step1));
        std::string label_LR(Form("h_energy_%s_bar%02d_LR_Vov%.1f","prova",barIDMax,step1));
        
        h1_energy["L"][step1][barIDMax] = new TH1F(label_L.c_str(),"",nEnergyBins,energyMin,energyMax);
        h1_energy["R"][step1][barIDMax] = new TH1F(label_R.c_str(),"",nEnergyBins,energyMin,energyMax);
        h1_energy["LR"][step1][barIDMax] = new TH1F(label_LR.c_str(),"",nEnergyBins,energyMin,energyMax);
      }
      
      if( (*qfine)[channelIdx[channelMapping[barIDMax*2+0]]] > qfineMin ) h1_energy["L"][step1][barIDMax] -> Fill( (*energy)[channelIdx[channelMapping[barIDMax*2+0]]] );
      if( (*qfine)[channelIdx[channelMapping[barIDMax*2+1]]] > qfineMin ) h1_energy["R"][step1][barIDMax] -> Fill( (*energy)[channelIdx[channelMapping[barIDMax*2+1]]] );
      if( (*qfine)[channelIdx[channelMapping[barIDMax*2+1]]] > qfineMin ) h1_energy["LR"][step1][barIDMax] -> Fill( 0.5 * ((*energy)[channelIdx[channelMapping[barIDMax*2+0]]] + (*energy)[channelIdx[channelMapping[barIDMax*2+1]]] ) );
    }
    
  } // loop over entries
  std::cout << std::endl;
  
  
  
  //--- draw plots
  for( auto mapIt : VovMap)
  {
    float Vov = mapIt.first;
    
    std::map<std::string,TGraphErrors*> g_511keV;
    for(auto histoLabel : histoLabels)
      g_511keV[histoLabel] = new TGraphErrors();
    
    float maxN = -999.;
    std::map<std::string,TGraphErrors*> g_N;
    for(auto histoLabel : histoLabels)
      g_N[histoLabel] = new TGraphErrors();
    
    
    for(int barIt = 0; barIt < 16; ++barIt)
      {
        TCanvas* c1 = new TCanvas(Form("c1_bar%d",barIt),Form("c1_bar%d",barIt));
        //c1 -> SetLogy();
	c1 -> cd(1);
        
        std::string label_L = Form("h1_energy_bar%02d_L_Vov%.1f", barIt,Vov);
        std::string label_R = Form("h1_energy_bar%02d_R_Vov%.1f", barIt,Vov);
        std::string label_LR = Form("h1_energy_bar%02d_LR_Vov%.1f", barIt,Vov);
        TH1F* histo_L = h1_energy["L"][Vov][barIt];
        TH1F* histo_R = h1_energy["R"][Vov][barIt];
        TH1F* histo_LR = h1_energy["LR"][Vov][barIt];
        
        if( histo_L == NULL || histo_R == NULL ) continue;
        
        // if( histo_L->Integral() > 0. ) histo_L->Scale(1./histo_L->Integral());
        // if( histo_R->Integral() > 0. ) histo_R->Scale(1./histo_R->Integral());
        
        float maxL = histo_L->GetMaximum();
        float maxR = histo_R->GetMaximum();
        
        histo_L -> SetMaximum(1.2*std::max(maxL,maxR));
        histo_R -> SetMaximum(1.2*std::max(maxL,maxR));
        histo_L -> GetXaxis() -> SetRangeUser(energyMin,energyMax);
        histo_R -> GetXaxis() -> SetRangeUser(energyMin,energyMax);
        
        histo_L -> SetTitle(";energy [a.u.];events");
        histo_L -> SetLineColor(kRed);
        histo_L -> Draw("hist");
        histo_R -> SetLineColor(kBlue);
        histo_R -> Draw("hist,sames");
        histo_LR -> SetLineColor(kBlack);
        histo_LR -> SetLineWidth(1);
        histo_LR -> Draw("hist,sames");
        
	for(auto histoLabel : histoLabels)
	  {
	    std::cout << "barIt: " << barIt << "   " << histoLabel << std::endl;
	    
	    int nPeaks = 5;
	    TSpectrum* spectrum = new TSpectrum(nPeaks);
	    int bin = h1_energy[histoLabel][Vov][barIt]->FindLastBinAbove(0);
	    float emin = h1_energy[histoLabel][Vov][barIt]->GetBinCenter(bin)*0.10;
	    h1_energy[histoLabel][Vov][barIt] -> GetXaxis()->SetRangeUser(emin,std::min(double(energyMax),900.));
	    int nFound = spectrum -> Search(h1_energy[histoLabel][Vov][barIt], 5., "nodraw", 0.1);
	    h1_energy[histoLabel][Vov][barIt] -> GetXaxis()->SetRangeUser(0,energyMax);
	    double* peaks = spectrum -> GetPositionX();
	    
	    // -- fit 511 keV peak only
	    //float xMax = FindXMaximum(histo_L,5.,20.);
	    float xMax = 0.;
	    for(int jj = 0; jj < nFound; ++jj)
	      if( peaks[jj] > xMax) xMax = peaks[jj];
	    
	    TF1* f_gaus = new TF1("f_gaus","gaus(0)",xMax-0.12*xMax,xMax+0.12*xMax);
	    h1_energy[histoLabel][Vov][barIt] -> Fit(f_gaus,"QNRS+");
	    f_gaus -> SetLineColor(h1_energy[histoLabel][Vov][barIt]->GetLineColor());
	    f_gaus -> SetLineWidth(5);
	    f_gaus -> Draw("same");
	    
	    g_511keV[histoLabel] -> SetPoint(g_511keV[histoLabel]->GetN(),barIt,f_gaus -> GetParameter(1));
	    g_511keV[histoLabel] -> SetPointError(g_511keV[histoLabel]->GetN()-1,0,f_gaus -> GetParError(1));
	    
	    float integral = f_gaus->Integral(f_gaus->GetParameter(1)-f_gaus->GetParameter(2),
					      f_gaus->GetParameter(1)+f_gaus->GetParameter(2));
	    
	    g_N[histoLabel] -> SetPoint(g_N[histoLabel]->GetN(),barIt,integral);
	    g_N[histoLabel] -> SetPointError(g_N[histoLabel]->GetN()-1,0,0.);
	    if( integral > maxN ) maxN = integral;
	  }
	
        c1 -> Print(Form("%s/c1_energy__bar%02d__Vov%.1f.png",plotDir.c_str(),barIt,Vov));
        delete c1;
      }
    
    
    TCanvas* c2 = new TCanvas("c2","c2");
    
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,energyMin,17.,energyMax) );
    hPad -> SetTitle(";bar ID;photopeak energy [a.u.]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    int jj = 0;
    for(auto histoLabel : histoLabels)
      {
	g_511keV[histoLabel] -> SetMarkerStyle(20+jj);
	g_511keV[histoLabel] -> SetMarkerColor(h1_energy[histoLabel][Vov][7]->GetLineColor());
	if( histoLabel == "LR" ) g_511keV[histoLabel] -> Draw("PL,same");
	else 	                 g_511keV[histoLabel] -> Draw("P,same");
	
	outFile -> cd();
	g_511keV[histoLabel] -> Write(Form("g_511keV_%s__Vov%.1f",histoLabel.c_str(),Vov));
	
	++jj;
      }
    
    c2 -> Print(Form("%s/c2_energy__Vov%.1f.png",plotDir.c_str(),Vov));
    
    
    TCanvas* c3 = new TCanvas("c3","c3");
    
    hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,17.,1.2*maxN) );
    hPad -> SetTitle(";bar ID;number of events");
    hPad -> Draw();
    gPad -> SetGridy();
    
    jj = 0;
    for(auto histoLabel : histoLabels)
      {
	g_N[histoLabel] -> SetMarkerStyle(20+jj);
	g_N[histoLabel] -> SetMarkerColor(h1_energy[histoLabel][Vov][7]->GetLineColor());
	if( histoLabel == "LR" ) g_N[histoLabel] -> Draw("PL,same");
	else 	                 g_N[histoLabel] -> Draw("P,same");
	
	outFile -> cd();
	g_N[histoLabel] -> Write(Form("g_N_%s__Vov%.1f",histoLabel.c_str(),Vov));
	
	++jj;
      }
    
    c3 -> Print(Form("%s/c3_N__Vov%.1f.png",plotDir.c_str(),Vov));
  }
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  
  
  return 0;
}
