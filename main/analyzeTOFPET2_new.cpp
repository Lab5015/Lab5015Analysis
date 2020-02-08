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



struct Event
{
  std::string stepLabel;
  std::string ch1;
  std::string ch2;
  std::string label1;
  std::string label2;
  std::string label12;
  int isBar1;
  int isBar2;
  float qfine1;
  float qfine1L;
  float qfine1R;
  float qfine2;
  float qfine2L;
  float qfine2R;
  float tot1;
  float tot1L;
  float tot1R;
  float tot2;
  float tot2L;
  float tot2R;
  float energy1;
  float energy1L;
  float energy1R;
  float energy2;
  float energy2L;
  float energy2R;
  long long time1;
  long long time2;
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
      for(int run = runMin; run <= runMax; ++run)
	{
	  //std::string fileName = Form("%s/*%04d*.root",inputDir.c_str(),run);
	  std::string fileName = Form("%s/%s%04d_t.root",inputDir.c_str(),fileBaseName.c_str(),run);
	  std::cout << ">>> Adding flle " << fileName << std::endl;
	  tree -> Add(fileName.c_str());
          
          struct stat t_stat;
          stat(Form("/data/TOFPET2/raw/run%04d.rawf",run), &t_stat);
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
  
  
  //--- define branches
  float step1, step2;
  unsigned int channelCount[256];
  float tot[256];
  float qfine[256];
  float energy[256];
  long long time[256];
  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("step1",       1); tree -> SetBranchAddress("step1",       &step1);
  tree -> SetBranchStatus("step2",       1); tree -> SetBranchAddress("step2",       &step2);
  tree -> SetBranchStatus("channelCount",1); tree -> SetBranchAddress("channelCount", channelCount);
  tree -> SetBranchStatus("qfine",       1); tree -> SetBranchAddress("qfine",        qfine);
  tree -> SetBranchStatus("tot",         1); tree -> SetBranchAddress("tot",          tot);
  tree -> SetBranchStatus("energy",      1); tree -> SetBranchAddress("energy",       energy);
  tree -> SetBranchStatus("time",        1); tree -> SetBranchAddress("time",         time);
  
  
  //--- define histograms
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileName");
  TFile* outFile = TFile::Open(Form("%s.root",outFileName.c_str()),"RECREATE");
  outFile -> cd();
  
  
  TTree* outTree = new TTree("data","data");
  unsigned long long int t_time = timesec;
  unsigned int t_arrayID1;
  unsigned int t_barID1;
  unsigned int t_lrID1;
  unsigned int t_chID1;
  unsigned int t_arrayID2;
  unsigned int t_barID2;
  unsigned int t_lrID2;
  unsigned int t_chID2;
  float t_Vov;
  float t_th;
  float t_tRes;
  float t_tResErr;
  outTree -> Branch("time",      &t_time,        "time/l");
  outTree -> Branch("arrayID1",  &t_arrayID1,"arrayID1/i");
  outTree -> Branch("barID1",    &t_barID1,    "barID1/i");
  outTree -> Branch("lrID1",     &t_lrID1,      "lrID1/i");
  outTree -> Branch("chID1",     &t_chID1,      "chID1/i");
  outTree -> Branch("arrayID2",  &t_arrayID2,"arrayID2/i");
  outTree -> Branch("barID2",    &t_barID2,    "barID2/i");
  outTree -> Branch("lrID2",     &t_lrID2,      "lrID2/i");
  outTree -> Branch("chID2",     &t_chID2,      "chID2/i");
  outTree -> Branch("tRes",      &t_tRes,        "tRes/F");
  outTree -> Branch("tResErr",   &t_tResErr,  "tResErr/F");
  
  
  std::map<std::string,int> VovLabels;
  std::map<std::string,int> thLabels;
  std::vector<std::string> stepLabels;
  std::map<std::string,float> map_Vovs;
  std::map<std::string,float> map_ths;
  
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
  
  
  
  
  //------------------------
  //--- 1st loop over events
  std::map<std::string,std::vector<Event> > events;
  
  int nEntries = tree->GetEntries();
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    tree -> GetEntry(entry);
    if( entry%10000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    
    float vth1 = float(int(step2/10000)-1);;
    // float vth2 = float(int((step2-10000*vth1)/100)-1);
    // float vthe = float(int((step2-10000*vth1-step2-100*vth2)/1)-1);
    
    std::string VovLabel = Form("Vov%.1f",step1);
    std::string thLabel = Form("th%02.0f",vth1);
    std::string stepLabel = Form("Vov%.1f_th%02.0f",step1,vth1);
    
        
    //--- create histograms, if needed
    for(auto ch : channels)
    {
      std::string label = Form("%s_%s",ch.c_str(),stepLabel.c_str());
      
      if( h1_tot[label] == NULL )
      {
        h1_qfine[label] = new TH1F(Form("h1_qfine_%s",label.c_str()),"",512,-0.5,511.5);
        h2_qfine_vs_tot[label] = new TH2F(Form("h2_qfine_vs_tot_%s",label.c_str()),"",100,0.,500,256,-0.5,255.5);
        
        h1_tot[label] = new TH1F(Form("h1_tot_%s",label.c_str()),"",2000,0.,1000.);
        
        h1_energy[label] = new TH1F(Form("h1_energy_%s",label.c_str()),"",1000,-10.,40.);
        h1_energy_cut[label] = new TH1F(Form("h1_energy_cut_%s",label.c_str()),"",1000,-10.,40.);
        
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
        h2_tot_corr[label12] = new TH2F(Form("h2_tot_corr_%s",label12.c_str()),"",100,0.,500.,100,0.,500.);
        h1_totRatio[label12] = new TH1F(Form("h1_totRatio_%s",label12.c_str()),"",1000,0.,5.);
        
        h2_energy_corr[label12] = new TH2F(Form("h2_energy_corr_%s",label12.c_str()),"",200,0.,50.,200,0.,50.);
        h1_energyRatio[label12] = new TH1F(Form("h1_energyRatio_%s",label12.c_str()),"",1000,0.,5.);
        
        h1_deltaT_raw[label12] = new TH1F(Form("h1_deltaT_raw_%s",label12.c_str()),"",1250,-5000.,5000.);
        h1_deltaT[label12] = new TH1F(Form("h1_deltaT_%s",label12.c_str()),"",1250,-5000,5000.);
        p1_deltaT_vs_energyRatio[label12] = new TProfile(Form("p1_deltaT_vs_energyRatio_%s",label12.c_str()),"",1000,0.,5.);
        
        h1_deltaT_energyCorr[label12] = new TH1F(Form("h1_deltaT_energyCorr_%s",label12.c_str()),"",1250,-5000.,5000.);
      }
    }
    
    
    //--- fill histograms
    for(auto ch : channels)
    {
      int chID = opts.GetOpt<int>(Form("%s.chID",ch.c_str()));
      std::string label = Form("%s_%s",ch.c_str(),stepLabel.c_str());
      
      float qfine1  = channelCount[int(chID)] == 1 ? qfine[int(chID)]     : -1.;
      float tot1    = channelCount[int(chID)] == 1 ? tot[int(chID)]/1000. : -1.;
      float energy1 = channelCount[int(chID)] == 1 ? energy[int(chID)]    : -1.;
      
      if( channelCount[int(chID)] != 1 ) continue;

      h1_qfine[label] -> Fill( qfine1 );
      h1_tot[label] -> Fill( tot1 );
      h2_qfine_vs_tot[label] -> Fill( tot1,qfine1 );
      
      if( qfine1 < cut_qfineAcc[chID][step1] ) continue;
      if( tot1 < cut_totAcc[chID][step1] ) continue;
      
      h1_energy[label] -> Fill( energy1 );
    }
    for(auto pair : pairsVec)
    {
      std::string ch1 = pair.first;
      std::string ch2 = pair.second;
      int isBar1 = opts.GetOpt<int>(Form("%s.isBar",ch1.c_str()));
      int isBar2 = opts.GetOpt<int>(Form("%s.isBar",ch2.c_str()));
      
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
        qfine1    = channelCount[int(chID1)] == 1 ? qfine[int(chID1)]     : -1.;
        tot1      = channelCount[int(chID1)] == 1 ? tot[int(chID1)]/1000. : -1.;
        energy1   = channelCount[int(chID1)] == 1 ? energy[int(chID1)]    : -1.;
        time1     = channelCount[int(chID1)] == 1 ? time[int(chID1)]      : -1.;
      }
      else
      {
        std::string channelL = opts.GetOpt<std::string>(Form("%s.channelL",ch1.c_str()));
        int chIDL = opts.GetOpt<int>(Form("%s.chID",channelL.c_str()));
        std::string channelR = opts.GetOpt<std::string>(Form("%s.channelR",ch1.c_str()));
        int chIDR = opts.GetOpt<int>(Form("%s.chID",channelR.c_str()));
        
        qfine1L  = channelCount[int(chIDL)] == 1 ? qfine[int(chIDL)]     : -1.;
        tot1L    = channelCount[int(chIDL)] == 1 ? tot[int(chIDL)]/1000. : -1.;
        energy1L = channelCount[int(chIDL)] == 1 ? energy[int(chIDL)]    : -1.;
        time1L   = channelCount[int(chIDL)] == 1 ? time[int(chIDL)]      : -1.;
        
        qfine1R  = channelCount[int(chIDR)] == 1 ? qfine[int(chIDR)]     : -1.;
        tot1R    = channelCount[int(chIDR)] == 1 ? tot[int(chIDR)]/1000. : -1.;
        energy1R = channelCount[int(chIDR)] == 1 ? energy[int(chIDR)]    : -1.;
        time1R   = channelCount[int(chIDR)] == 1 ? time[int(chIDR)]      : -1.;
        
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
        qfine2    = channelCount[int(chID2)] == 1 ? qfine[int(chID2)]     : -1.;
        tot2      = channelCount[int(chID2)] == 1 ? tot[int(chID2)]/1000. : -1.;
        energy2   = channelCount[int(chID2)] == 1 ? energy[int(chID2)]    : -1.;
        time2     = channelCount[int(chID2)] == 1 ? time[int(chID2)]      : -1.;
      }
      else
      {
        std::string channelL = opts.GetOpt<std::string>(Form("%s.channelL",ch2.c_str()));
        int chIDL = opts.GetOpt<int>(Form("%s.chID",channelL.c_str()));
        std::string channelR = opts.GetOpt<std::string>(Form("%s.channelR",ch2.c_str()));
        int chIDR = opts.GetOpt<int>(Form("%s.chID",channelR.c_str()));        
        
        qfine2L  = channelCount[int(chIDL)] == 1 ? qfine[int(chIDL)]     : -1.;
        tot2L    = channelCount[int(chIDL)] == 1 ? tot[int(chIDL)]/1000. : -1.;
        energy2L = channelCount[int(chIDL)] == 1 ? energy[int(chIDL)]    : -1.;
        time2L   = channelCount[int(chIDL)] == 1 ? time[int(chIDL)]      : -1.;
        
        qfine2R  = channelCount[int(chIDR)] == 1 ? qfine[int(chIDR)]     : -1.;
        tot2R    = channelCount[int(chIDR)] == 1 ? tot[int(chIDR)]/1000. : -1.;
        energy2R = channelCount[int(chIDR)] == 1 ? energy[int(chIDR)]    : -1.;
        time2R   = channelCount[int(chIDR)] == 1 ? time[int(chIDR)]      : -1.;
        
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
        float coincidenceEnergy = channelCount[int(coincidenceID)] == 1 ? energy[int(coincidenceID)] : -1.;
        if( coincidenceEnergy < cut_energyAcc[coincidenceID][step1] ) accept = false;
      }
      if( !accept ) continue;
      
      
      Event anEvent;
      anEvent.stepLabel = stepLabel;
      anEvent.ch1 = ch1;
      anEvent.ch2 = ch2;
      anEvent.label1 = label1;
      anEvent.label2 = label2;
      anEvent.label12 = label12;
      anEvent.isBar1 = isBar1;
      anEvent.isBar2 = isBar2;
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
      events[label12].push_back(anEvent);
    }
  }
  std::cout << std::endl;
  
  std::vector<std::string>::iterator iter;
  iter = std::unique(stepLabels.begin(),stepLabels.end());
  stepLabels.resize( std::distance(stepLabels.begin(),iter) );  
  
  std::sort(stepLabels.begin(),stepLabels.end());
  
  
  
  
  for(auto stepLabel : stepLabels)
  {
    float Vov = map_Vovs[stepLabel];
    float th = map_ths[stepLabel];
    std::string VovLabel(Form("Vov%.1f",Vov));
    std::string thLabel(Form("th%02.0f",th));
    
    
    for(auto ch : channels)
    {
      int chID = opts.GetOpt<int>(Form("%s.chID",ch.c_str()));
      std::string label(Form("%s_%s",ch.c_str(),stepLabel.c_str()));
      
      float max1 = FindXMaximum(h1_energy[label],cut_energyAcc[chID][Vov],50.);
      TF1* fitFunc = new TF1("fitFunc","gaus",max1-cut_energyFitMin[chID][Vov]*max1,max1+cut_energyFitMax[chID][Vov]*max1);
      h1_energy[label] -> Fit(fitFunc,"QRS+");
      cut_energyMin[Form("%s_%s",ch.c_str(),stepLabel.c_str())] = fitFunc->GetMaximumX()-cut_energyFitMin[chID][Vov]*fitFunc->GetMaximumX();
      cut_energyMax[Form("%s_%s",ch.c_str(),stepLabel.c_str())] = fitFunc->GetMaximumX()+cut_energyFitMax[chID][Vov]*fitFunc->GetMaximumX();
    }
  }
  
  
  
  
  //------------------------
  //--- 2nd loop over events
  std::map<std::string,std::vector<Event> > events2;
  
  for(auto mapIt : events)
  {
    std::string label = mapIt.first;
    float Vov = map_Vovs[label];
        
    nEntries = mapIt.second.size();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 2nd loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
      Event anEvent = mapIt.second.at(entry);
      
      
      if( anEvent.isBar1 == 0 )
      {
        int chID1 = opts.GetOpt<int>(Form("%s.chID",anEvent.ch1.c_str()));
        
        if( anEvent.qfine1 < cut_qfineAcc[chID1][Vov] ) continue;
        if( anEvent.tot1 < cut_totAcc[chID1][Vov] ) continue;
        if( anEvent.energy1 < cut_energyMin[Form("%s_%s",anEvent.ch1.c_str(),anEvent.stepLabel.c_str())] ) continue;
        if( anEvent.energy1 > cut_energyMax[Form("%s_%s",anEvent.ch1.c_str(),anEvent.stepLabel.c_str())] ) continue;
      }
      else
      {
        std::string channelL = opts.GetOpt<std::string>(Form("%s.channelL",anEvent.ch1.c_str()));
        int chIDL = opts.GetOpt<int>(Form("%s.chID",channelL.c_str()));
        std::string channelR = opts.GetOpt<std::string>(Form("%s.channelR",anEvent.ch1.c_str()));
        int chIDR = opts.GetOpt<int>(Form("%s.chID",channelR.c_str()));
        
        if( anEvent.qfine1L < cut_qfineAcc[chIDL][Vov] ) continue;
        if( anEvent.tot1L < cut_totAcc[chIDL][Vov] ) continue;
        if( anEvent.energy1L < cut_energyMin[Form("%s_%s",channelL.c_str(),anEvent.stepLabel.c_str())] ) continue;
        if( anEvent.energy1L > cut_energyMax[Form("%s_%s",channelL.c_str(),anEvent.stepLabel.c_str())] ) continue;
        
        if( anEvent.qfine1R < cut_qfineAcc[chIDR][Vov] ) continue;
        if( anEvent.tot1R < cut_totAcc[chIDR][Vov] ) continue;
        if( anEvent.energy1R < cut_energyMin[Form("%s_%s",channelR.c_str(),anEvent.stepLabel.c_str())] ) continue;
        if( anEvent.energy1R > cut_energyMax[Form("%s_%s",channelR.c_str(),anEvent.stepLabel.c_str())] ) continue;
      }
      
      
      if( anEvent.isBar2 == 0 )
      {
        int chID2 = opts.GetOpt<int>(Form("%s.chID",anEvent.ch2.c_str()));
        
        if( anEvent.qfine2 < cut_qfineAcc[chID2][Vov] ) continue;
        if( anEvent.tot2 < cut_totAcc[chID2][Vov] ) continue;
        if( anEvent.energy2 < cut_energyMin[Form("%s_%s",anEvent.ch2.c_str(),anEvent.stepLabel.c_str())] ) continue;
        if( anEvent.energy2 > cut_energyMax[Form("%s_%s",anEvent.ch2.c_str(),anEvent.stepLabel.c_str())] ) continue;
      }
      else
      {
        std::string channelL = opts.GetOpt<std::string>(Form("%s.channelL",anEvent.ch2.c_str()));
        int chIDL = opts.GetOpt<int>(Form("%s.chID",channelL.c_str()));
        std::string channelR = opts.GetOpt<std::string>(Form("%s.channelR",anEvent.ch2.c_str()));
        int chIDR = opts.GetOpt<int>(Form("%s.chID",channelR.c_str()));
        
        if( anEvent.qfine2L < cut_qfineAcc[chIDL][Vov] ) continue;
        if( anEvent.tot2L < cut_totAcc[chIDL][Vov] ) continue;
        if( anEvent.energy2L < cut_energyMin[Form("%s_%s",channelL.c_str(),anEvent.stepLabel.c_str())] ) continue;
        if( anEvent.energy2L > cut_energyMax[Form("%s_%s",channelL.c_str(),anEvent.stepLabel.c_str())] ) continue;

        if( anEvent.qfine2R < cut_qfineAcc[chIDR][Vov] ) continue;
        if( anEvent.tot2R < cut_totAcc[chIDR][Vov] ) continue;
        if( anEvent.energy2R < cut_energyMin[Form("%s_%s",channelR.c_str(),anEvent.stepLabel.c_str())] ) continue;
        if( anEvent.energy2R > cut_energyMax[Form("%s_%s",channelR.c_str(),anEvent.stepLabel.c_str())] ) continue;
      }
      
      
      h1_totRatio[anEvent.label12] -> Fill( anEvent.tot2 / anEvent.tot1 );
      h1_energyRatio[anEvent.label12] -> Fill( anEvent.energy2 / anEvent.energy1 );
      
      h1_deltaT_raw[anEvent.label12] -> Fill( anEvent.time2-anEvent.time1 );
      
      
      Event anEvent2;
      anEvent2.stepLabel = anEvent.stepLabel;
      anEvent2.ch1 = anEvent.ch1;
      anEvent2.ch2 = anEvent.ch2;
      anEvent2.label1 = anEvent.label1;
      anEvent2.label2 = anEvent.label2;
      anEvent2.label12 = anEvent.label12;
      anEvent2.isBar1 = anEvent.isBar1;
      anEvent2.isBar2 = anEvent.isBar2;
      anEvent2.qfine1 = anEvent.qfine1;
      anEvent2.qfine1L = anEvent.qfine1L;
      anEvent2.qfine1R = anEvent.qfine1R;
      anEvent2.qfine2 = anEvent.qfine2;
      anEvent2.qfine2L = anEvent.qfine2L;
      anEvent2.qfine2R = anEvent.qfine2R;
      anEvent2.tot1 = anEvent.tot1;
      anEvent2.tot1L = anEvent.tot1L;
      anEvent2.tot1R = anEvent.tot1R;
      anEvent2.tot2 = anEvent.tot2;
      anEvent2.tot2L = anEvent.tot2L;
      anEvent2.tot2R = anEvent.tot2R;
      anEvent2.energy1 = anEvent.energy1;
      anEvent2.energy1L = anEvent.energy1L;
      anEvent2.energy1R = anEvent.energy1R;
      anEvent2.energy2 = anEvent.energy2;
      anEvent2.energy2L = anEvent.energy2L;
      anEvent2.energy2R = anEvent.energy2R;
      anEvent2.time1 = anEvent.time1;
      anEvent2.time2 = anEvent.time2;
      events2[anEvent.label12].push_back(anEvent2);
    }
    std::cout << std::endl;
  }
  
  events.clear();
  
  
  
  
  std::map<std::string,float> CTRMeans;
  std::map<std::string,float> CTRSigmas;
  
  float* vals = new float[6];
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
      
      
      FindSmallestInterval(vals,h1_deltaT_raw[label12],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      CTRMeans[label12] = mean;
      CTRSigmas[label12] = effSigma;
    }
  }
  
  
  
  
  //------------------------
  //--- 3rd loop over events
  for(auto mapIt : events2)
  {
    std::string label = mapIt.first;
    
    nEntries = mapIt.second.size();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 3rd loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
      Event anEvent = mapIt.second.at(entry);
      
      float timeLow = CTRMeans[anEvent.label12] - 2.* CTRSigmas[anEvent.label12];
      float timeHig = CTRMeans[anEvent.label12] + 2.* CTRSigmas[anEvent.label12];
      long long deltaT = anEvent.time2 - anEvent.time1;
      
      h1_deltaT[anEvent.label12] -> Fill( deltaT );
      if( ( deltaT > timeLow ) &&
          ( deltaT < timeHig ) )
        p1_deltaT_vs_energyRatio[anEvent.label12] -> Fill( anEvent.energy2/anEvent.energy1,anEvent.time2-anEvent.time1 );    
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
      std::string label2(Form("%s_%s",ch2.c_str(),stepLabel.c_str()));
      std::string label12 = Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),stepLabel.c_str());
      
      float fitXMin = h1_energyRatio[label12]->GetMean() - 2.*h1_energyRatio[label12]->GetRMS();
      float fitXMax = h1_energyRatio[label12]->GetMean() + 2.*h1_energyRatio[label12]->GetRMS();
      fitFunc_energyCorr[label12] = new TF1(Form("fitFunc_energyCorr_%s",label12.c_str()),"pol4",fitXMin,fitXMax);
      p1_deltaT_vs_energyRatio[label12] -> Fit(fitFunc_energyCorr[label12],"QRS+");
    }
  }
  
  
  
  
  //------------------------
  //--- 4th loop over events
  for(auto mapIt : events2)
  {
    std::string label = mapIt.first;
    
    nEntries = mapIt.second.size();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 4th loop (" << label << "): reading entry " << entry << " / " << nEntries << "\r" << std::flush;
      Event anEvent = mapIt.second.at(entry);
      
      long long deltaT = anEvent.time2-anEvent.time1;
      
      float energyCorr = fitFunc_energyCorr[label]->Eval(anEvent.energy2/anEvent.energy1) -
                         fitFunc_energyCorr[label]->Eval(h1_energyRatio[anEvent.label12]->GetMean());
      h1_deltaT_energyCorr[label] -> Fill( deltaT - energyCorr );
    }
    std::cout << std::endl;
  }
  
  
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
