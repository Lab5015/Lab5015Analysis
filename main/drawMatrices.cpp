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



struct channel
{
  std::string array;
  int arrayID;
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
    std::string fileName = Form("%s/%s_events.root",inputDir.c_str(),fileBaseName.c_str());
    std::cout << ">>> Adding flle " << fileName << std::endl;
    tree -> Add(fileName.c_str());
    
    struct stat t_stat;
    stat(Form("/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Feb2020/TOFHIR/RawData/%s.rawf",fileBaseName.c_str()), &t_stat);
    struct tm * timeinfo = localtime(&t_stat.st_mtime);
    timesec = mktime(timeinfo);
    std::cout << "Time and date of raw file of run " << fileBaseName << ": " << asctime(timeinfo);
  }
  
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  system(Form("mkdir -p %s",plotDir.c_str()));
  
  
  
  std::vector<unsigned int> channelMapping1 = opts.GetOpt<std::vector<unsigned int> >("Channels.channelMapping1");
  std::vector<unsigned int> channelMapping2 = opts.GetOpt<std::vector<unsigned int> >("Channels.channelMapping2");
  std::vector<unsigned int> channelMapping;
  
  std::map<int,channel> channels_key;
  std::vector<std::string> arrays = opts.GetOpt<std::vector<std::string> >("Channels.arrays");
  for(auto array : arrays)
  {
    unsigned int FEBDPort = opts.GetOpt<unsigned int>(Form("%s.FEBDPort",array.c_str()));
    std::string connectorID = opts.GetOpt<std::string>(Form("%s.connectorID",array.c_str()));
    int arrayID = opts.GetOpt<unsigned int>(Form("%s.arrayID",array.c_str()));
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
    
    for(int ii = 0; ii < 32; ++ii)
    {
      int chID = channelMapping[ii] + offset;
      channels_key[chID] = {array,arrayID,(ii/2),(ii%2)};
    }
  }
  
  
  //--- define branches
  float energy[256];
  float qfine[256];
  float tot[256];
  
  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("energy",      1); tree -> SetBranchAddress("energy",      energy);
  tree -> SetBranchStatus("qfine",       1); tree -> SetBranchAddress("qfine",       qfine);
  tree -> SetBranchStatus("tot",         1); tree -> SetBranchAddress("tot",         tot);
  
  
  
  //--- define histograms
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileName");
  TFile* outFile = TFile::Open(Form("%s",outFileName.c_str()),"RECREATE");
  outFile -> cd();
  
  
  // TTree* outTree = new TTree("data","data");
  // unsigned long long int t_time = timesec;
  // unsigned int t_arrayID;
  // unsigned int t_barID;
  // unsigned int t_lrID;
  // float t_peak1;
  // float t_peak1Err;
  // float t_sigma1;
  // float t_sigma1Err;
  // outTree -> Branch("time",     &t_time,          "time/l");
  // outTree -> Branch("arrayID",  &t_arrayID,    "arrayID/i");
  // outTree -> Branch("barID",    &t_barID,        "barID/i");
  // outTree -> Branch("lrID",     &t_lrID,          "lrID/i");
  // outTree -> Branch("peak1",    &t_peak1,        "peak1/F");
  // outTree -> Branch("peak1Err", &t_peak1Err,  "peak1Err/F");
  // outTree -> Branch("sigma1",   &t_sigma1,      "sigma1/F");
  // outTree -> Branch("sigma1Err",&t_sigma1Err,"sigma1Err/F");
  
  
  std::map<std::string,TH1F*> h1_energy;
  std::map<std::string,TH2F*> h2_energy_LRCorr;
  
  for(auto mapIt: channels_key)
  {
    if( mapIt.second.lrID != 0 ) continue;
    
    std::string array = mapIt.second.array;
    int arrayID = mapIt.second.arrayID;
    int barID = mapIt.second.barID;
    
    std::string label_L = Form("h1_energy_%s_bar%02d_L",array.c_str(),barID);
    std::string label_R = Form("h1_energy_%s_bar%02d_R",array.c_str(),barID);
    
    h1_energy[label_L] = new TH1F(label_L.c_str(),"",500,0.,50.);
    h1_energy[label_R] = new TH1F(label_R.c_str(),"",500,0.,50.);
    
    std::string label_LR = Form("h2_energy_%s_bar%02d_LRCorr",array.c_str(),barID);
    h2_energy_LRCorr[label_LR] = new TH2F(label_LR.c_str(),"",500,0.,50.,500,0.,50.);
  }
  
  
  
  //--- loop over events
  int nEntries = tree->GetEntries();
  if( maxEntries > nEntries ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    tree -> GetEntry(entry);
    if( entry%100000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    
    for(auto mapIt: channels_key)
    {
      int chID    = mapIt.first;
      std::string array = mapIt.second.array;
      int arrayID = mapIt.second.arrayID;
      int barID   = mapIt.second.barID;
      int lrID    = mapIt.second.lrID;
      
      if( qfine[chID] > 13 )
      {
        if( lrID == 0 ) h1_energy[Form("h1_energy_%s_bar%02d_L",array.c_str(),barID)] -> Fill( energy[chID] );
        else            h1_energy[Form("h1_energy_%s_bar%02d_R",array.c_str(),barID)] -> Fill( energy[chID] );
        
        if( lrID == 0 )
        {
          for(auto mapIt2: channels_key)
          {
            if( mapIt2.second.arrayID == arrayID &&
                mapIt2.second.barID == barID &&
                mapIt2.second.lrID == 1 )
            {
              int chID2 = mapIt2.first;
              
              h2_energy_LRCorr[Form("h2_energy_%s_bar%02d_LRCorr",array.c_str(),barID)] -> Fill( energy[chID],energy[chID2] );
            }
          }
        }
        
      }
    } 
  }
  std::cout << std::endl;
  
  
  
  //--- draw plots
  for(auto array : arrays)
  {
    TGraphErrors* g_511keV_L = new TGraphErrors();
    TGraphErrors* g_511keV_R = new TGraphErrors();
    
    float maxN = -999.;
    TGraphErrors* g_N_L = new TGraphErrors();
    TGraphErrors* g_N_R = new TGraphErrors();
    
    
    for(int barIt = 0; barIt < 16; ++barIt)
    {
      // t_arrayID = arrayIt;
      // t_barID = barIt;
      
      TCanvas* c1 = new TCanvas(Form("c1_bar%d",barIt),Form("c1_bar%d",barIt),1200,600);
      c1 -> Divide(2,1);
      
      c1 -> cd(1);
      
      std::string label_L = Form("h1_energy_%s_bar%02d_L", array.c_str(),barIt);
      std::string label_R = Form("h1_energy_%s_bar%02d_R", array.c_str(),barIt);
      
      if( h1_energy[label_L] == NULL || h1_energy[label_R] == NULL ) continue;
      
      if( h1_energy[label_L]->Integral() > 0. ) h1_energy[label_L]->Scale(1./h1_energy[label_L]->Integral());
      if( h1_energy[label_R]->Integral() > 0. ) h1_energy[label_R]->Scale(1./h1_energy[label_R]->Integral());
      
      float max = h1_energy[label_L] -> GetMaximum();
      if( h1_energy[label_R]->GetMaximum() > max ) max = h1_energy[label_R]->GetMaximum();
      h1_energy[label_L]->SetMaximum(max);
      
      h1_energy[label_L] -> SetTitle(";energy [a.u.];fraction of events");
      h1_energy[label_L] -> SetLineColor(kRed);
      h1_energy[label_L] -> Draw("hist");
      h1_energy[label_R] -> SetLineColor(kBlue);
      h1_energy[label_R] -> Draw("hist,sames");
      
      g_N_L -> SetPoint(g_N_L->GetN(),barIt,h1_energy[label_L]->GetEntries());
      g_N_L -> SetPointError(g_N_L->GetN()-1,0,sqrt(h1_energy[label_L]->GetEntries()));
      if( h1_energy[label_L]->GetEntries() > maxN ) maxN = h1_energy[label_L]->GetEntries();
      g_N_R -> SetPoint(g_N_R->GetN(),barIt,h1_energy[label_R]->GetEntries());
      g_N_R -> SetPointError(g_N_R->GetN()-1,0,sqrt(h1_energy[label_R]->GetEntries()));
      if( h1_energy[label_R]->GetEntries() > maxN ) maxN = h1_energy[label_R]->GetEntries();
      
      
      std::vector<std::string> labels;
      labels.push_back(label_L);
      labels.push_back(label_R);
      
      for(auto label: labels)
      {
        // if( label == label_L ) t_lrID = 0;
        // if( label == label_R ) t_lrID = 1;
        
        
        TH1F* histo = h1_energy[label];
        if( histo->GetEntries() < 100 ) continue;
        
        
        // -- fit 511 keV peak only
        float xMax = FindXMaximum(histo,5.,20.);
        TF1* f_gaus = new TF1("f_gaus","gaus(0)",xMax-0.05*xMax,xMax+0.05*xMax);
        histo -> Fit(f_gaus,"QNRS+");
        f_gaus -> SetLineColor(kBlack);
        f_gaus -> SetLineWidth(3);
        f_gaus -> Draw("same");
        
        if( label == label_L )
        {
          g_511keV_L -> SetPoint(g_511keV_L->GetN(),barIt,f_gaus -> GetParameter(1));
          g_511keV_L -> SetPointError(g_511keV_L->GetN()-1,0,f_gaus -> GetParError(1));
        }
        if( label == label_R )
        {
          g_511keV_R -> SetPoint(g_511keV_R->GetN(),barIt,f_gaus -> GetParameter(1));
          g_511keV_R -> SetPointError(g_511keV_R->GetN()-1,0,f_gaus -> GetParError(1));
        }
      }
      
      
      c1 -> cd(2);
      
      std::string label_LR = Form("h2_energy_%s_bar%02d_LRCorr", array.c_str(),barIt);
      h2_energy_LRCorr[label_LR] -> Draw("COLZ");
      
      c1 -> Print(Form("%s/c1_energy__%s__bar%02d.png",plotDir.c_str(),array.c_str(),barIt));
      delete c1;
    }
    
    TCanvas* c2 = new TCanvas("c2","c2");
    
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,17.,30.) );
    hPad -> SetTitle(";bar ID;photopeak energy [a.u.]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    g_511keV_L -> SetMarkerStyle(20);
    g_511keV_L -> SetMarkerColor(kRed);
    g_511keV_L -> Draw("P,same");
    g_511keV_R -> SetMarkerColor(kBlue);
    g_511keV_R -> SetMarkerStyle(21);
    g_511keV_R -> Draw("P,same");
    
    c2 -> Print(Form("%s/c2_energy__%s.png",plotDir.c_str(),array.c_str()));
    
    
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
    
    c3 -> Print(Form("%s/c3_N__%s.png",plotDir.c_str(),array.c_str()));
  }
  
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  
  
  return 0;
}
