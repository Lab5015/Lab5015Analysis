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
  std::string fileBaseName = opts.GetOpt<std::string>("Input.fileBaseName");
  std::string runs = opts.GetOpt<std::string>("Input.runs");
  int maxEntries = opts.GetOpt<int>("Input.maxEntries");
  TChain* tree = new TChain("data","data");
  
  std::stringstream ss(runs);
  std::string token;
  time_t timesec;
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
      //std::string fileName = Form("%s/*%04d*.root",inputDir.c_str(),run);
      std::string fileName = Form("%s/%s%04d_events.root",inputDir.c_str(),fileBaseName.c_str(),run);
      std::cout << ">>> Adding file " << fileName << std::endl;
      tree -> Add(fileName.c_str());
      
      struct stat t_stat;
      stat(Form("/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_FNAL_Feb2020/TOFHIR/RawData/run%04d.rawf",run), &t_stat);
      struct tm * timeinfo = localtime(&t_stat.st_mtime);
      timesec = mktime(timeinfo);
      std::cout << "Time and date of raw file of run" << run << ": " << asctime(timeinfo);
    }
  }
  
  
  //--- define channels
  std::vector<std::string> channels = opts.GetOpt<std::vector<std::string> >("Channels.channels");
  
  std::vector<std::string> pairs = opts.GetOpt<std::vector<std::string> >("Channels.pairs");
  std::vector<std::pair<std::string,std::string> > pairsVec;
  for(unsigned int ii = 0; ii < pairs.size()/2; ++ii)
  {
    pairsVec.push_back(std::make_pair(pairs.at(0+ii*2),pairs.at(1+ii*2)));
  }
  
  
  
  //--- get cuts per bar / Vov
  std::map<unsigned int,std::map<float,float> > cut_qfineAcc;
  std::map<unsigned int,std::map<float,float> > cut_qfineMax;
  std::map<unsigned int,std::map<float,float> > cut_totAcc;
  std::map<unsigned int,std::map<float,float> > cut_energyAcc;
  std::map<unsigned int,std::map<float,float> > cut_energyFitMin;
  std::map<unsigned int,std::map<float,float> > cut_energyFitMax;
  for(auto ch :  channels)
  {
    int chID = opts.GetOpt<int>(Form("%s.chID",ch.c_str()));
    std::vector<float> Vovs          = opts.GetOpt<std::vector<float> >(Form("%s.Vovs",ch.c_str()));
    std::vector<float> qfineMins     = opts.GetOpt<std::vector<float> >(Form("%s.qfineMins",ch.c_str()));
    std::vector<float> qfineMaxs     = opts.GetOpt<std::vector<float> >(Form("%s.qfineMaxs",ch.c_str()));
    std::vector<float> totMins       = opts.GetOpt<std::vector<float> >(Form("%s.totMins",ch.c_str()));
    std::vector<float> energyMins    = opts.GetOpt<std::vector<float> >(Form("%s.energyMins",ch.c_str()));
    std::vector<float> energyFitMins = opts.GetOpt<std::vector<float> >(Form("%s.energyFitMins",ch.c_str()));
    std::vector<float> energyFitMaxs = opts.GetOpt<std::vector<float> >(Form("%s.energyFitMaxs",ch.c_str()));
    int iter = 0;
    for(auto Vov : Vovs)
    {
      cut_qfineAcc[chID][Vov]     = qfineMins.at(iter);
      cut_qfineMax[chID][Vov]     = qfineMaxs.at(iter);
      cut_totAcc[chID][Vov]       = totMins.at(iter);
      cut_energyAcc[chID][Vov]    = energyMins.at(iter);
      cut_energyFitMin[chID][Vov] = energyFitMins.at(iter);
      cut_energyFitMax[chID][Vov] = energyFitMaxs.at(iter);
      ++iter;
    }
  }
  std::map<std::string,float> cut_energyMin;
  std::map<std::string,float> cut_energyMax;
  
  
  //--- define branches
  float step1, step2;
  float tot[256];
  float qfine[256];
  float energy[256];
  long long time[256];
  float xIntercept;
  float yIntercept;
  int ntracks;
  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("step1",       1); tree -> SetBranchAddress("step1",       &step1);
  tree -> SetBranchStatus("step2",       1); tree -> SetBranchAddress("step2",       &step2);
  tree -> SetBranchStatus("step2",       1); tree -> SetBranchAddress("step2",       &step2);
  tree -> SetBranchStatus("qfine",       1); tree -> SetBranchAddress("qfine",        qfine);
  tree -> SetBranchStatus("tot",         1); tree -> SetBranchAddress("tot",          tot);
  tree -> SetBranchStatus("energy",      1); tree -> SetBranchAddress("energy",       energy);
  tree -> SetBranchStatus("time",        1); tree -> SetBranchAddress("time",         time);
  tree -> SetBranchStatus("xIntercept",  1); tree -> SetBranchAddress("xIntercept",  &xIntercept);
  tree -> SetBranchStatus("yIntercept",  1); tree -> SetBranchAddress("yIntercept",  &yIntercept);
  tree -> SetBranchStatus("ntracks",     1); tree -> SetBranchAddress("ntracks",     &ntracks);
  
  
  //--- define histograms
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileName");
  TFile* outFile = TFile::Open(Form("%s.root",outFileName.c_str()),"RECREATE");
  outFile -> cd();
  
  std::map<std::string,TTree*> outTrees;
  
  std::map<std::string,int> VovLabels;
  std::map<std::string,int> thLabels;
  std::vector<std::string> stepLabels;
  std::map<std::string,float> map_Vovs;
  std::map<std::string,float> map_ths;
  
  std::map<std::string,TProfile2D*> p2_eff_vs_XY;
  
  std::map<std::string,TH1F*> h1_qfine;
  std::map<std::string,TH2F*> h2_qfine_vs_tot;
  
  std::map<std::string,TH1F*> h1_tot;
  std::map<std::string,TH2F*> h2_tot_corr;
  std::map<std::string,TH1F*> h1_totRatio;
  
  std::map<std::string,TH1F*> h1_energy;
  std::map<std::string,TH1F*> h1_energy_cut;
  std::map<std::string,TH2F*> h2_energy_corr;
  std::map<std::string,TH1F*> h1_energyRatio;
  
  std::map<std::string,TH1F*> h1_deltaT_raw;
  std::map<std::string,TH1F*> h1_deltaT;
  std::map<std::string,TProfile*> p1_deltaT_vs_energyRatio;
  
  std::map<std::string,TH1F*> h1_deltaT_energyCorr;
  std::map<std::string,TProfile*> p1_deltaT_energyCorr_vs_pos;
  
  std::map<std::string,TH1F*> h1_deltaT_energyCorr_posCorr;
  
  
  //------------------------
  //--- 1st loop over events
  int nEntries = tree->GetEntries();
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    tree -> GetEntry(entry);
    if( entry%10000 == 0 )
    {
      std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << std::endl;
      TrackProcess(cpu, mem, vsz, rss);
    }
    
    float vth1 = float(int(step2/10000)-1);;
    // float vth2 = float(int((step2-10000*vth1)/100)-1);
    // float vthe = float(int((step2-10000*vth1-step2-100*vth2)/1)-1);
    
    std::string VovLabel = Form("Vov%.1f",step1);
    std::string thLabel = Form("th%02.0f",vth1);
    std::string stepLabel = Form("Vov%.1f_th%02.0f",step1,vth1);
    
    EventClass anEvent;
    
    
    //--- create histograms, if needed
    for(auto ch : channels)
    {
      std::string label = Form("%s_%s",ch.c_str(),stepLabel.c_str());
      
      if( h1_tot[label] == NULL )
      {
        p2_eff_vs_XY[label] = new TProfile2D(Form("p2_eff_vs_XY_%s",label.c_str()),"",200,-10.,40.,200,0.,50.);
  
        h1_qfine[label] = new TH1F(Form("h1_qfine_%s",label.c_str()),"",512,-0.5,511.5);
        h2_qfine_vs_tot[label] = new TH2F(Form("h2_qfine_vs_tot_%s",label.c_str()),"",100,0.,500,512,-0.5,511.5);
        
        h1_tot[label] = new TH1F(Form("h1_tot_%s",label.c_str()),"",2000,0.,1000.);
        
        h1_energy[label] = new TH1F(Form("h1_energy_%s",label.c_str()),"",1000,0.,50.);
        h1_energy_cut[label] = new TH1F(Form("h1_energy_cut_%s",label.c_str()),"",1000,0.,50.);
        
        VovLabels[VovLabel] += 1;
        thLabels[thLabel] += 1;
        stepLabels.push_back(stepLabel);
        map_Vovs[stepLabel] = step1;
        map_ths[stepLabel] = vth1;
      }
    }
    for(auto pair : pairsVec)
    {
      std::string ch1 = pair.first;
      std::string ch2 = pair.second;
      std::string label12 = Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),stepLabel.c_str());
      
      if( h2_tot_corr[label12] == NULL )
      {
        outTrees[label12] = new TTree(Form("data_%s",label12.c_str()),Form("data_%s",label12.c_str()));
        outTrees[label12] -> Branch("event",&anEvent);
        
        h2_tot_corr[label12] = new TH2F(Form("h2_tot_corr_%s",label12.c_str()),"",100,0.,500.,100,0.,500.);
        h1_totRatio[label12] = new TH1F(Form("h1_totRatio_%s",label12.c_str()),"",1000,0.,5.);
        
        h2_energy_corr[label12] = new TH2F(Form("h2_energy_corr_%s",label12.c_str()),"",200,0.,50.,200,0.,50.);
        h1_energyRatio[label12] = new TH1F(Form("h1_energyRatio_%s",label12.c_str()),"",1000,0.,5.);
        
        h1_deltaT_raw[label12] = new TH1F(Form("h1_deltaT_raw_%s",label12.c_str()),"",1250,-5000.,5000.);
        h1_deltaT[label12] = new TH1F(Form("h1_deltaT_%s",label12.c_str()),"",1250,-5000,5000.);
        
        h1_deltaT_energyCorr[label12] = new TH1F(Form("h1_deltaT_energyCorr_%s",label12.c_str()),"",1250,-5000.,5000.);
        h1_deltaT_energyCorr_posCorr[label12] = new TH1F(Form("h1_deltaT_energyCorr_posCorr_%s",label12.c_str()),"",1250,-5000.,5000.);
      }
    }
    
    
    //--- fill histograms
    for(auto ch : channels)
    {
      int chID = opts.GetOpt<int>(Form("%s.chID",ch.c_str()));
      std::string label = Form("%s_%s",ch.c_str(),stepLabel.c_str());
      
      float qfine1  = qfine[int(chID)];
      float tot1    = tot[int(chID)]/1000.;
      float energy1 = energy[int(chID)];
      
      //p2_eff_vs_XY[label] -> Fill( xIntercept,yIntercept,energy1>cut_energyAcc[chID][step1] );
      int weight = (energy1>cut_energyAcc[chID][step1]) && (qfine1 > cut_qfineAcc[chID][step1]) && (tot1 > cut_totAcc[chID][step1]);
      p2_eff_vs_XY[label] -> Fill( xIntercept,yIntercept,weight );
      h1_qfine[label] -> Fill( qfine1 );
      h1_tot[label] -> Fill( tot1 );
      h2_qfine_vs_tot[label] -> Fill( tot1,qfine1 );
      
      if( qfine1 < cut_qfineAcc[chID][step1] ) continue;
      if( qfine1 > cut_qfineMax[chID][step1] ) continue;
      if( tot1 < cut_totAcc[chID][step1] ) continue;
      
      h1_energy[label] -> Fill( energy1 );
    }
    for(auto pair : pairsVec)
    {
      std::string ch1 = pair.first;
      std::string ch2 = pair.second;
      int isBar1 = opts.GetOpt<int>(Form("%s.isBar",ch1.c_str()));
      int isBar2 = opts.GetOpt<int>(Form("%s.isBar",ch2.c_str()));
      int isHorizontal1 = opts.GetOpt<int>(Form("%s.isHorizontal",ch1.c_str()));
      int isHorizontal2 = opts.GetOpt<int>(Form("%s.isHorizontal",ch2.c_str()));
      
      std::string label1 = Form("%s_%s",ch1.c_str(),stepLabel.c_str());
      std::string label2 = Form("%s_%s",ch2.c_str(),stepLabel.c_str());
      std::string label12 = Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),stepLabel.c_str());
      
      
      float qfine1, qfine1L, qfine1R;
      float tot1, tot1L, tot1R;
      float energy1, energy1L, energy1R;
      long long time1, time1L, time1R;
      if( isBar1 == 0 )
      {
        int chID1 = opts.GetOpt<int>(Form("%s.chID",ch1.c_str()));        
        qfine1    = qfine[int(chID1)];
        tot1      = tot[int(chID1)]/1000.;
        energy1   = energy[int(chID1)];
        time1     = time[int(chID1)];
      }
      else
      {
        std::string channelL = opts.GetOpt<std::string>(Form("%s.channelL",ch1.c_str()));
        int chIDL = opts.GetOpt<int>(Form("%s.chID",channelL.c_str()));
        std::string channelR = opts.GetOpt<std::string>(Form("%s.channelR",ch1.c_str()));
        int chIDR = opts.GetOpt<int>(Form("%s.chID",channelR.c_str()));
        
        qfine1L  = qfine[int(chIDL)];
        tot1L    = tot[int(chIDL)]/1000.;
        energy1L = energy[int(chIDL)];
        time1L   = time[int(chIDL)];
        
        qfine1R  = qfine[int(chIDR)];
        tot1R    = tot[int(chIDR)]/1000.;
        energy1R = energy[int(chIDR)];
        time1R   = time[int(chIDR)];
        
        qfine1    = 0.5*(qfine1L+qfine1R);
        tot1      = 0.5*(tot1L+tot1R);
        energy1   = 0.5*(energy1L+energy1R);
        time1     = 0.5*(time1L+time1R);
      }
      
      
      float qfine2, qfine2L, qfine2R;
      float tot2, tot2L, tot2R;
      float energy2, energy2L, energy2R;
      long long time2, time2L, time2R;
      if( isBar2 == 0 )
      {
        int chID2 = opts.GetOpt<int>(Form("%s.chID",ch2.c_str()));        
        qfine2    = qfine[int(chID2)];
        tot2      = tot[int(chID2)]/1000.;
        energy2   = energy[int(chID2)];
        time2     = time[int(chID2)];
      }
      else
      {
        std::string channelL = opts.GetOpt<std::string>(Form("%s.channelL",ch2.c_str()));
        int chIDL = opts.GetOpt<int>(Form("%s.chID",channelL.c_str()));
        std::string channelR = opts.GetOpt<std::string>(Form("%s.channelR",ch2.c_str()));
        int chIDR = opts.GetOpt<int>(Form("%s.chID",channelR.c_str()));        
        
        qfine2L  = qfine[int(chIDL)];
        tot2L    = tot[int(chIDL)]/1000.;
        energy2L = energy[int(chIDL)];
        time2L   = time[int(chIDL)];
        
        qfine2R  = qfine[int(chIDR)];
        tot2R    = tot[int(chIDR)]/1000.;
        energy2R = energy[int(chIDR)];
        time2R   = time[int(chIDR)];
        
        qfine2    = 0.5*(qfine2L+qfine2R);
        tot2      = 0.5*(tot2L+tot2R);
        energy2   = 0.5*(energy2L+energy2R);
        time2     = 0.5*(time2L+time2R);
      }
      
      
      if( tot1 < -1. || tot2 < -1. ) continue;
      
      h2_tot_corr[label12] -> Fill( tot1,tot2 );
      h2_tot_corr[label12] -> Fill( tot1,tot2 );
      h2_energy_corr[label12] -> Fill( energy1,energy2 );
      h2_energy_corr[label12] -> Fill( energy1,energy2 );
      
      
      bool accept = true;
      std::vector<std::string> coincidenceChs1 = opts.GetOpt<std::vector<std::string> >(Form("%s.coincidenceCh",ch1.c_str()));
      for(auto coincidenceCh : coincidenceChs1)
      {
        if( coincidenceCh == "NULL" ) continue;
        
        int coincidenceID = opts.GetOpt<int>(Form("%s.chID",coincidenceCh.c_str()));
        float coincidenceEnergy = energy[int(coincidenceID)];
        if( coincidenceEnergy < cut_energyAcc[coincidenceID][step1] ) accept = false;
      }
      if( !accept ) continue;
      
      
      anEvent.stepLabel = stepLabel;
      anEvent.ch1 = ch1;
      anEvent.ch2 = ch2;
      anEvent.label1 = label1;
      anEvent.label2 = label2;
      anEvent.label12 = label12;
      anEvent.x = xIntercept;
      anEvent.y = yIntercept;
      anEvent.isBar1 = isBar1;
      anEvent.isBar2 = isBar2;
      anEvent.isHorizontal1 = isHorizontal1;
      anEvent.isHorizontal2 = isHorizontal2;
      anEvent.qfine1 = qfine1;
      anEvent.qfine1L = qfine1L;
      anEvent.qfine1R = qfine1R;
      anEvent.qfine2 = qfine2;
      anEvent.qfine2L = qfine2L;
      anEvent.qfine2R = qfine2R;
      anEvent.tot1 = tot1;
      anEvent.tot1L = tot1L;
      anEvent.tot1R = tot1R;
      anEvent.tot2 = tot2;
      anEvent.tot2L = tot2L;
      anEvent.tot2R = tot2R;
      anEvent.energy1 = energy1;
      anEvent.energy1L = energy1L;
      anEvent.energy1R = energy1R;
      anEvent.energy2 = energy2;
      anEvent.energy2L = energy2L;
      anEvent.energy2R = energy2R;
      anEvent.time1 = time1;
      anEvent.time2 = time2;
      outTrees[label12] -> Fill();
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
