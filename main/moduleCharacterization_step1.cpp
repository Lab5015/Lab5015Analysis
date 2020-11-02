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






int main(int argc, char** argv)
{
  setTDRStyle();
  float cpu[2]{0}, mem[2]={0}, vsz[2]={0}, rss[2]={0};  
  
  
  if( argc < 2 )
  {
    std::cout << ">>> moduleCharacterization_step1::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }
  
  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  
  
  
  //--- open files and make the tree chain
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  std::string fileBaseName = opts.GetOpt<std::string>("Input.fileBaseName");
  std::string runs = opts.GetOpt<std::string>("Input.runs");
  int maxEntries = opts.GetOpt<int>("Input.maxEntries");
  int usePedestals = opts.GetOpt<int>("Input.usePedestals");
  TChain* tree = new TChain("data","data");
  
  std::stringstream ss(runs);
  std::string token;
  while( std::getline(ss,token,',') )
  {
    std::stringstream ss2(token);
    std::string token2;
    int runMin = -1;
    int runMax = -1;
    while( std::getline(ss2,token2,'-') )
    {
      if( runMin != -1 && runMax == -1 ) runMax = atoi(token2.c_str());
      if( runMin == -1 ) runMin = atoi(token2.c_str());
    }
    if( runMax == -1 ) runMax = runMin;
    
    for(int run = runMin; run <= runMax; ++run)
    {
      std::string fileName;
      if( !usePedestals ) fileName = Form("%s/%s%04d_t.root",inputDir.c_str(),fileBaseName.c_str(),run);
      else                fileName = Form("%s/%s%04d_ped_t.root",inputDir.c_str(),fileBaseName.c_str(),run);
      std::cout << ">>> Adding file " << fileName << std::endl;
      tree -> Add(fileName.c_str());
      
      struct stat t_stat;
      stat(Form("/data/TOFPET2/raw/run%04d.rawf",run), &t_stat);
      struct tm * timeinfo = localtime(&t_stat.st_mtime);
      std::cout << "Time and date of raw file of run" << run << ": " << asctime(timeinfo);
    }
  }
  
  
  
  //--- define channels
  int chsL[16];
  int chsR[16];
  chsL[0] = 249;
  chsR[0] = 4;
  chsL[1] = 242;
  chsR[1] = 6;
  chsL[2] = 251;
  chsR[2] = 15;
  chsL[3] = 253;
  chsR[3] = 8;
  chsL[4] = 250;
  chsR[4] = 13;
  chsL[5] = 254;
  chsR[5] = 2;
  chsL[6] = 201;
  chsR[6] = 11;
  chsL[7] = 230;
  chsR[7] = 27;
  chsL[8] = 226;
  chsR[8] = 32;
  chsL[9] = 225;
  chsR[9] = 31;
  chsL[10] = 228;
  chsR[10] = 30;
  chsL[11] = 227;
  chsR[11] = 29;
  chsL[12] = 229;
  chsR[12] = 28;
  chsL[13] = 231;
  chsR[13] = 26;
  chsL[14] = 233;
  chsR[14] = 24;
  chsL[15] = 232;
  chsR[15] = 25;
  
  
  
  //--- define branches
  float step1, step2;
  float tot[256];
  float qfine[256];
  float energy[256];
  long long time[256];
  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("step1",  1); tree -> SetBranchAddress("step1",  &step1);
  tree -> SetBranchStatus("step2",  1); tree -> SetBranchAddress("step2",  &step2);
  tree -> SetBranchStatus("qfine",  1); tree -> SetBranchAddress("qfine",   qfine);
  tree -> SetBranchStatus("tot",    1); tree -> SetBranchAddress("tot",       tot);
  tree -> SetBranchStatus("energy", 1); tree -> SetBranchAddress("energy", energy);
  tree -> SetBranchStatus("time",   1); tree -> SetBranchAddress("time",     time);
  
  
  
  //--- get plot settings
  std::vector<int> Vov = opts.GetOpt<std::vector<int> >("Plots.Vov");
  std::vector<int> energyBins = opts.GetOpt<std::vector<int> >("Plots.energyBins");
  std::map<int,int> map_energyBins;
  for(unsigned int ii = 0; ii < Vov.size(); ++ii)
    map_energyBins[Vov[ii]] = energyBins[ii];
  
  float energyMin = opts.GetOpt<float>("Plots.energyMin");
  float energyMax = opts.GetOpt<float>("Plots.energyMax");
  
  
  
  //--- define histograms
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileNameStep1");
  TFile* outFile = TFile::Open(Form("%s",outFileName.c_str()),"RECREATE");
  outFile -> cd();
  
  std::map<int,TTree*> outTrees;
  
  std::map<int,TH1F*> h1_qfineL;
  std::map<int,TH1F*> h1_qfineR;
  std::map<int,TH1F*> h1_totL;
  std::map<int,TH1F*> h1_totR;
  std::map<int,TH1F*> h1_energyL;
  std::map<int,TH1F*> h1_energyR;
  std::map<int,TH1F*> h1_energyLR;
  
  
  
  //------------------------
  //--- 1st loop over events
  ModuleEventClass anEvent;
  
  int nEntries = tree->GetEntries();
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    tree -> GetEntry(entry);
    if( entry%200000 == 0 )
    {
      std::cout << "\n>>> 1st loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << std::endl;
      TrackProcess(cpu, mem, vsz, rss);
    }
    
    float Vov = step1;
    float vth1 = float(int(step2/10000)-1);;
    // float vth2 = float(int((step2-10000*vth1)/100)-1);
    // float vthe = float(int((step2-10000*vth1-step2-100*vth2)/1)-1);
    
    
    for(int iBar = 0; iBar < 16; ++iBar)
    {
      int index( (10000*int(Vov*100.)) + (100*vth1) + iBar );
      
      
      int chL(chsL[iBar]);
      int chR(chsR[iBar]);
              
      float qfineL(qfine[chL]);
      float qfineR(qfine[chR]);
      float totL(0.001*tot[chL]);
      float totR(0.001*tot[chR]);
      
      if( totL <= 0. || totR <= 0. ) continue;
      if( qfineL < 13. || qfineR < 13. ) continue;
      
      float energyL(energy[chL]);
      float energyR(energy[chR]);
      long long timeL(time[chL]);
      long long timeR(time[chR]);
      
      
      //--- create histograms, if needed
      if( h1_totL[index] == NULL )
      {
        h1_qfineL[index] = new TH1F(Form("h1_qfine_bar%02dL_Vov%.1f_th%02.0f",iBar,Vov,vth1),"",512,-0.5,511.5);
        h1_qfineR[index] = new TH1F(Form("h1_qfine_bar%02dR_Vov%.1f_th%02.0f",iBar,Vov,vth1),"",512,-0.5,511.5);
        
        h1_totL[index] = new TH1F(Form("h1_tot_bar%02dL_Vov%.1f_th%02.0f",iBar,Vov,vth1),"",2000,0.,1000.);
        h1_totR[index] = new TH1F(Form("h1_tot_bar%02dR_Vov%.1f_th%02.0f",iBar,Vov,vth1),"",2000,0.,1000.);
        
        h1_energyL[index] = new TH1F(Form("h1_energy_bar%02dL_Vov%.1f_th%02.0f",iBar,Vov,vth1),"",map_energyBins[Vov],energyMin,energyMax);
        h1_energyR[index] = new TH1F(Form("h1_energy_bar%02dR_Vov%.1f_th%02.0f",iBar,Vov,vth1),"",map_energyBins[Vov],energyMin,energyMax);
        
        outTrees[index] = new TTree(Form("data_bar%02dL-R_Vov%.1f_th%02.0f",iBar,Vov,vth1),Form("data_bar%02dL-R_Vov%.1f_th%02.0f",iBar,Vov,vth1));
        outTrees[index] -> Branch("event",&anEvent);
        
        h1_energyLR[index] = new TH1F(Form("h1_energy_bar%02dL-R_Vov%.1f_th%02.0f",iBar,Vov,vth1),"",map_energyBins[Vov],energyMin,energyMax);
      }
      
      
      //--- fill histograms
      h1_qfineL[index] -> Fill( qfineL );
      h1_totL[index] -> Fill( totL );
      h1_energyL[index] -> Fill( energyL );
      
      h1_qfineR[index] -> Fill( qfineR );
      h1_totR[index] -> Fill( totR );
      h1_energyR[index] -> Fill( energyR );
      
      h1_energyLR[index] -> Fill(0.5*(energyL+energyR));
      
      anEvent.barID = iBar;
      anEvent.Vov = Vov;
      anEvent.vth1 = vth1;
      anEvent.energyL = energyL;
      anEvent.energyR = energyR;
      anEvent.timeL = timeL;
      anEvent.timeR = timeR;
      outTrees[index] -> Fill();
    }
  }
  std::cout << std::endl;
  
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
