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
  std::string array;
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
  
  
  
  std::vector<unsigned int> channelMapping1 = opts.GetOpt<std::vector<unsigned int> >("Channels.channelMapping1");
  std::vector<unsigned int> channelMapping2 = opts.GetOpt<std::vector<unsigned int> >("Channels.channelMapping2");
  std::vector<unsigned int> channelMapping;
  
  std::vector<std::string> arrays = opts.GetOpt<std::vector<std::string> >("Channels.arrays");
  for(auto array : arrays)
  {
    unsigned int FEBDPort = opts.GetOpt<unsigned int>(Form("%s.FEBDPort",array.c_str()));
    std::string connectorID = opts.GetOpt<std::string>(Form("%s.connectorID",array.c_str()));
    int offset = 128*FEBDPort;
    
    if( connectorID == "S1_1" || connectorID == "S2_1" )
    {
      channelMapping = channelMapping1;
    }
    else if( connectorID == "S1_2" || connectorID == "S2_2" )
    {
      channelMapping = channelMapping2;
    }
    else continue;
    
    if( connectorID == "S2_1" || connectorID == "S2_2" )
    {
      offset += 64;
    }
    
    for(unsigned int ii = 0; ii < channelMapping.size(); ++ii)
      channelMapping[ii] += offset;
  }
  
  
  
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
  
  std::map<int,std::map<int,TH1F*> > h1_energy_L;
  std::map<int,std::map<int,TH1F*> > h1_energy_R;
  std::map<int,std::map<int,TH1F*> > h1_energy_LR;
  // std::map<int,TH2F*> h2_energy_LRCorr;
  
  
  
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
      if( h1_energy_L[step1][barIDMax] == NULL )
      {
        std::string label_L(Form("h1_energy_%s_bar%02d_L_Vov%.1f","prova",barIDMax,step1));
        std::string label_R(Form("h1_energy_%s_bar%02d_R_Vov%.1f","prova",barIDMax,step1));
        std::string label_LR(Form("h_energy_%s_bar%02d_LR_Vov%.1f","prova",barIDMax,step1));

        
        h1_energy_L[step1][barIDMax] = new TH1F(label_L.c_str(),"",nEnergyBins,energyMin,energyMax);
        h1_energy_R[step1][barIDMax] = new TH1F(label_R.c_str(),"",nEnergyBins,energyMin,energyMax);
        h1_energy_LR[step1][barIDMax] = new TH1F(label_LR.c_str(),"",nEnergyBins,energyMin,energyMax);
        //h1_energy_L[step1][barID] = new TH1F(label_L.c_str(),"",150,0.,150.);
        //h1_energy_R[step1][barID] = new TH1F(label_R.c_str(),"",150,0.,150.);
        // h2_energy_LRCorr[barID] = new TH2F(label_LR.c_str(),"",150,0.,150.,150,0.,150.);
      }
      
      if( (*qfine)[channelIdx[channelMapping[barIDMax*2+0]]] > qfineMin ) h1_energy_L[step1][barIDMax] -> Fill( (*energy)[channelIdx[channelMapping[barIDMax*2+0]]] );
      if( (*qfine)[channelIdx[channelMapping[barIDMax*2+1]]] > qfineMin ) h1_energy_R[step1][barIDMax] -> Fill( (*energy)[channelIdx[channelMapping[barIDMax*2+1]]] );
      if( (*qfine)[channelIdx[channelMapping[barIDMax*2+1]]] > qfineMin ) h1_energy_LR[step1][barIDMax] -> Fill( 0.5 * ((*energy)[channelIdx[channelMapping[barIDMax*2+0]]] + (*energy)[channelIdx[channelMapping[barIDMax*2+1]]] ) );
    }
    
    
    /*
    for(auto mapIt: channels_key)
    {
      int chID    = mapIt.first;
      std::string array = mapIt.second.array;
      int barID   = mapIt.second.barID;
      int lrID    = mapIt.second.lrID;
      
      if( qfine[chID] > qfineMin )
      {
        if( h1_energy_L[step1][barID] == NULL )
        {
          std::string label_L(Form("h1_energy_%s_bar%02d_L_Vov%.1f",array.c_str(),barID,step1));
          std::string label_R(Form("h1_energy_%s_bar%02d_R_Vov%.1f",array.c_str(),barID,step1));
          // std::string label_LR(Form("h2_energy_%s_bar%02d_LRCorr_Vov%.1f",array.c_str(),barID,step1));
          
          h1_energy_L[step1][barID] = new TH1F(label_L.c_str(),"",400,0.,40.);
          h1_energy_R[step1][barID] = new TH1F(label_R.c_str(),"",400,0.,40.);
          //h1_energy_L[step1][barID] = new TH1F(label_L.c_str(),"",150,0.,150.);
          //h1_energy_R[step1][barID] = new TH1F(label_R.c_str(),"",150,0.,150.);
          // h2_energy_LRCorr[barID] = new TH2F(label_LR.c_str(),"",150,0.,150.,150,0.,150.);
        }
        
        if( lrID == 0 ) h1_energy_L[step1][barID] -> Fill( energy[chID] );
        else            h1_energy_R[step1][barID] -> Fill( energy[chID] );
        
        // if( lrID == 0 )
        // {
        //   for(auto mapIt2: channels_key)
        //   {
        //         mapIt2.second.barID == barID &&
        //         mapIt2.second.lrID == 1 )
        //     {
        //       int chID2 = mapIt2.first;
              
        //       h2_energy_LRCorr[barID] -> Fill( energy[chID],energy[chID2] );
        //     }
        //   }
        // }
      }
    } */
  }
  std::cout << std::endl;
  
  
  
  //--- draw plots
  for( auto mapIt : VovMap)
  {
    float Vov = mapIt.first;
    
    for(auto array : arrays)
    {
      TGraphErrors* g_511keV_L = new TGraphErrors();
      TGraphErrors* g_511keV_R = new TGraphErrors();
      TGraphErrors* g_511keV_LR = new TGraphErrors();
      
      float maxN = -999.;
      TGraphErrors* g_N_L = new TGraphErrors();
      TGraphErrors* g_N_R = new TGraphErrors();
      TGraphErrors* g_N_LR = new TGraphErrors();
      
      
      for(int barIt = 0; barIt < 16; ++barIt)
      {
        // TCanvas* c1 = new TCanvas(Form("c1_bar%d",barIt),Form("c1_bar%d",barIt),1200,600);
        // c1 -> Divide(2,1);
        TCanvas* c1 = new TCanvas(Form("c1_bar%d",barIt),Form("c1_bar%d",barIt));
        c1 -> SetLogy();
        // c1 -> Divide(2,1);
        // c1 -> cd(1);
        
        std::string label_L = Form("h1_energy_%s_bar%02d_L_Vov%.1f", array.c_str(),barIt,Vov);
        std::string label_R = Form("h1_energy_%s_bar%02d_R_Vov%.1f", array.c_str(),barIt,Vov);
        std::string label_LR = Form("h1_energy_%s_bar%02d_LR_Vov%.1f", array.c_str(),barIt,Vov);
        TH1F* histo_L = h1_energy_L[Vov][barIt];
        TH1F* histo_R = h1_energy_R[Vov][barIt];
        TH1F* histo_LR = h1_energy_LR[Vov][barIt];
        
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
        histo_LR -> SetLineWidth(2);
        histo_LR -> Draw("hist,sames");
        
        if( histo_L->Integral(histo_L->FindBin(7.),histo_L->FindBin(20.)) > 100 )
        {
          std::cout << "barIt: " << barIt << "   left " << std::endl;
          
          int nPeaks = 5;
          TSpectrum* spectrum = new TSpectrum(nPeaks);
          int bin = histo_L->FindLastBinAbove(0);
          float emin = histo_L->GetBinCenter(bin)*0.25;
          histo_L-> GetXaxis()->SetRangeUser(emin,energyMax);
          int nFound = spectrum -> Search(histo_L, 5.0, "nodraw", 0.01);
          histo_L-> GetXaxis()->SetRangeUser(0,energyMax);
          double* peaks = spectrum -> GetPositionX();
          
          // -- fit 511 keV peak only
          //float xMax = FindXMaximum(histo_L,5.,20.);
          float xMax = 0.;
          for(int jj = 0; jj < nFound; ++jj)
            if( peaks[jj] > xMax) xMax = peaks[jj];
	  
          TF1* f_gaus = new TF1("f_gaus","gaus(0)",xMax-0.12*xMax,xMax+0.12*xMax);
          histo_L -> Fit(f_gaus,"QNRS+");
          f_gaus -> SetLineColor(kRed);
          f_gaus -> SetLineWidth(3);
          f_gaus -> Draw("same");
          
          g_511keV_L -> SetPoint(g_511keV_L->GetN(),barIt,f_gaus -> GetParameter(1));
          g_511keV_L -> SetPointError(g_511keV_L->GetN()-1,0,f_gaus -> GetParError(1));
          
          float integral = f_gaus->Integral(f_gaus->GetParameter(1)-f_gaus->GetParameter(2),
                                            f_gaus->GetParameter(1)+f_gaus->GetParameter(2));
          
          g_N_L -> SetPoint(g_N_L->GetN(),barIt,integral);
          g_N_L -> SetPointError(g_N_L->GetN()-1,0,0.);
          if( integral > maxN ) maxN = integral;
        }
        
        if( histo_R->Integral(histo_R->FindBin(5.),histo_R->FindBin(20.)) > 100 )
        {
          std::cout << "barIt: " << barIt << "   right " << std::endl;
          
          int nPeaks = 5;
          TSpectrum* spectrum = new TSpectrum(nPeaks);
          int bin = histo_R->FindLastBinAbove(0);
          float emin = histo_R->GetBinCenter(bin)*0.25;
          histo_R-> GetXaxis()->SetRangeUser(emin,energyMax);
          int nFound = spectrum -> Search(histo_R, 5., "nodraw", 0.01);
          histo_R-> GetXaxis()->SetRangeUser(0,energyMax);
          double* peaks = spectrum -> GetPositionX();
          

          // -- fit 511 keV peak only
          // float xMax = FindXMaximum(histo_R,5.,20.);
          float xMax = 0.;
          for(int jj = 0; jj < nFound; ++jj)
            if( peaks[jj] > xMax) xMax = peaks[jj];
          
          TF1* f_gaus = new TF1("f_gaus","gaus(0)",xMax-0.12*xMax,xMax+0.12*xMax);
          histo_R -> Fit(f_gaus,"QNRS+");
          f_gaus -> SetLineColor(kBlue);
          f_gaus -> SetLineWidth(3);
          f_gaus -> Draw("same");
          
          g_511keV_R -> SetPoint(g_511keV_R->GetN(),barIt,f_gaus -> GetParameter(1));
          g_511keV_R -> SetPointError(g_511keV_R->GetN()-1,0,f_gaus -> GetParError(1));
          
          float integral = f_gaus->Integral(f_gaus->GetParameter(1)-f_gaus->GetParameter(2),
                                            f_gaus->GetParameter(1)+f_gaus->GetParameter(2));
          
          g_N_R -> SetPoint(g_N_R->GetN(),barIt,integral);
          g_N_R -> SetPointError(g_N_R->GetN()-1,0,0.);
          if( integral > maxN ) maxN = integral;
        }
        
        if( histo_LR->Integral(histo_LR->FindBin(3.),histo_LR->FindBin(20.)) > 100 )
        {
          std::cout << "barIt: " << barIt << "   left-right " << std::endl;
          
          int nPeaks = 5;
          TSpectrum* spectrum = new TSpectrum(nPeaks);
          int bin = histo_LR->FindLastBinAbove(0);
          float emin = histo_LR->GetBinCenter(bin)*0.25;
          histo_LR-> GetXaxis()->SetRangeUser(emin,energyMax);
          int nFound = spectrum -> Search(histo_LR, 5., "nodraw", 0.01);
          double* peaks = spectrum -> GetPositionX();
	  histo_LR-> GetXaxis()->SetRangeUser(energyMin,energyMax);
		    
          // -- fit 511 keV peak only
          // float xMax = FindXMaximum(histo_R,5.,20.);
          float xMax = 0.;
          for(int jj = 0; jj < nFound; ++jj)
            if( peaks[jj] > xMax) xMax = peaks[jj];
          TF1* f_gaus = new TF1("f_gaus","gaus(0)",xMax-0.12*xMax,xMax+0.12*xMax);
          histo_LR -> Fit(f_gaus,"QNRS+");
          f_gaus -> SetLineColor(kBlack);
          f_gaus -> SetLineWidth(3);
          f_gaus -> Draw("same");
          
          g_511keV_LR -> SetPoint(g_511keV_LR->GetN(),barIt,f_gaus -> GetParameter(1));
          g_511keV_LR -> SetPointError(g_511keV_LR->GetN()-1,0,f_gaus -> GetParError(1));
          
          float integral = f_gaus->Integral(f_gaus->GetParameter(1)-f_gaus->GetParameter(2),
                                            f_gaus->GetParameter(1)+f_gaus->GetParameter(2));
          
          g_N_LR -> SetPoint(g_N_LR->GetN(),barIt,integral);
          g_N_LR -> SetPointError(g_N_LR->GetN()-1,0,0.);
          if( integral > maxN ) maxN = integral;
        }
        
        
        // c1 -> cd(2);
        
        // std::string label_LR = Form("h2_energy_%s_bar%02d_LRCorr_Vov%.1f", array.c_str(),barIt,Vov);
        // h2_energy_LRCorr[label_LR] -> Draw("COLZ");
        
        c1 -> Print(Form("%s/c1_energy__%s__bar%02d__Vov%.1f.png",plotDir.c_str(),array.c_str(),barIt,Vov));
        delete c1;
      }
      
      
      TCanvas* c2 = new TCanvas("c2","c2");
      
      TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,energyMin,17.,energyMax) );
      hPad -> SetTitle(";bar ID;photopeak energy [a.u.]");
      hPad -> Draw();
      gPad -> SetGridy();
      
      g_511keV_L -> SetMarkerStyle(20);
      g_511keV_L -> SetMarkerColor(kRed);
      g_511keV_L -> Draw("P,same");
      g_511keV_R -> SetMarkerColor(kBlue);
      g_511keV_R -> SetMarkerStyle(21);
      g_511keV_R -> Draw("P,same");
      g_511keV_LR -> SetLineColor(kBlack);
      g_511keV_LR -> SetMarkerColor(kBlack);
      g_511keV_LR -> SetMarkerStyle(22);
      g_511keV_LR -> Draw("PL,same");
      
      outFile -> cd();
      g_511keV_L -> Write(Form("g_511keV_L__%s__Vov%.1f",array.c_str(),Vov));
      g_511keV_R -> Write(Form("g_511keV_R__%s__Vov%.1f",array.c_str(),Vov));
      g_511keV_LR -> Write(Form("g_511keV_LR__%s__Vov%.1f",array.c_str(),Vov));
      
      c2 -> Print(Form("%s/c2_energy__%s__Vov%.1f.png",plotDir.c_str(),array.c_str(),Vov));
      
      
      TCanvas* c3 = new TCanvas("c3","c3");
      
      hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,17.,1.2*maxN) );
      hPad -> SetTitle(";bar ID;number of events");
      hPad -> Draw();
      gPad -> SetGridy();
      
      g_N_L -> SetMarkerStyle(20);
      g_N_L -> SetMarkerColor(kRed);
      g_N_L -> Draw("P,same");
      g_N_R -> SetMarkerColor(kBlue);
      g_N_R -> SetMarkerStyle(21);
      g_N_R -> Draw("P,same");
      g_N_LR -> SetMarkerColor(kBlack);
      g_N_LR -> SetMarkerStyle(22);
      g_N_LR -> Draw("P,same");
      
      outFile -> cd();
      g_N_L -> Write(Form("g_N_L__%s__Vov%.1f",array.c_str(),Vov));
      g_N_R -> Write(Form("g_N_R__%s__Vov%.1f",array.c_str(),Vov));
      g_N_LR -> Write(Form("g_N_LR__%s__Vov%.1f",array.c_str(),Vov));
      
      c3 -> Print(Form("%s/c3_N__%s__Vov%.1f.png",plotDir.c_str(),array.c_str(),Vov));
    }
  }
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  
  
  return 0;
}
