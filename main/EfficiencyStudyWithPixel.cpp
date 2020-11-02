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


  int mystep1 = 6; // VoV
  int mystep2 = 211606; // threshold configuration
  int chRef   = 249; // channel for coincidence (conf9.2)

  int denominatorType = atoi(argv[2]);
  //int denominatorType = 0; // only track
  //int denominatorType = 1; // track + qfine
  //int denominatorType = 2; // track + qfine + pixel energy
  //int denominatorType = 3; // track + qfine + pixel energy + pixel position
  //int denominatorType = 4; // track + pixel energy
  
  
  if( argc < 2 )
    {
      std::cout << ">>> analyzeBars::usage:   " << argv[0] << " configFile.cfg    denominatorType [0/1/2/3]" << std::endl;
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


  //--- define branches
  int run;
  float step1, step2;
  float tot[256];
  float qfine[256];
  float energy[256];
  long long time[256];
  float xIntercept;
  float yIntercept;
  int ntracks;
  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("run",       1); tree -> SetBranchAddress("run",       &run);
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
  std::string suffix = "preselTracks";
  if (denominatorType == 1 ) suffix = "preselTracksQfine";
  if (denominatorType == 2 ) suffix = "preselTracksQfinePixelEnergy";
  if (denominatorType == 3 ) suffix = "preselTracksQfinePixelEnergyPosition";
  if (denominatorType == 4 ) suffix = "preselTracksPixelEnergy";

  TFile* outFile = TFile::Open(Form("%s_%s.root",outFileName.c_str(),suffix.c_str()),"RECREATE");
  std::cout << "Output file name : " << outFile->GetName() <<std::endl;
  outFile -> cd();

  std::map<std::string,TProfile*> p_eff_vs_X_effCutEnergy;  
  std::map<std::string,TProfile2D*> p2_eff_vs_XY_effCutEnergy;

  std::map<std::string,TProfile*> p_eff_vs_X_effCutEnergyQfine;  
  std::map<std::string,TProfile2D*> p2_eff_vs_XY_effCutEnergyQfine;

  std::map<std::string,TProfile*> p_eff_vs_X_effCutEnergyQfineTot;  
  std::map<std::string,TProfile2D*> p2_eff_vs_XY_effCutEnergyQfineTot;

  std::map<std::string,TH2F*> h2_qfine_vs_tot;
  std::map<std::string,TH2F*> h2_qfine_vs_energy;
  std::map<std::string,TH2F*> h2_qfine1_vs_qfine2;
  std::map<std::string,TH2F*> h2_energy1_vs_energy2;

  //------------------------
  //--- 1st loop over events
  int nEntries = tree->GetEntries();
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
    {
      tree -> GetEntry(entry);
      if( entry%50000 == 0 )
	{
	  std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << std::endl;
	  std::cout << "Run " << run <<std::endl;
	}

      // -- channel ref = 34 for conf 11.3
      if (run>=26828) chRef = 34;

      //-- select step1, step2
      if (step1!=mystep1) continue;
      if (step2!=mystep2) continue;
      
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
      
	  if( p2_eff_vs_XY_effCutEnergy[label] == NULL )
	    {

	      p2_eff_vs_XY_effCutEnergy[label] = new TProfile2D(Form("p2_eff_vs_XY_effCutEnergy_%s",label.c_str()),"",200,-10.,40.,200,0.,50.);
	      p_eff_vs_X_effCutEnergy[label] = new TProfile(Form("p_eff_vs_X_effCutEnergy_%s",label.c_str()),"",200,-10.,40.);

	      p2_eff_vs_XY_effCutEnergyQfine[label] = new TProfile2D(Form("p2_eff_vs_XY_effCutEnergyQfine_%s",label.c_str()),"",200,-10.,40.,200,0.,50.);
	      p_eff_vs_X_effCutEnergyQfine[label] = new TProfile(Form("p_eff_vs_X_effCutEnergyQfine_%s",label.c_str()),"",200,-10.,40.);

	      p2_eff_vs_XY_effCutEnergyQfineTot[label] = new TProfile2D(Form("p2_eff_vs_XY_effCutEnergyQfineTot_%s",label.c_str()),"",200,-10.,40.,200,0.,50.);
	      p_eff_vs_X_effCutEnergyQfineTot[label] = new TProfile(Form("p_eff_vs_X_effCutEnergyQfineTot_%s",label.c_str()),"",200,-10.,40.);

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

	  // preselection cuts on tracks and qfine
	  bool accept = 0;
	  if (denominatorType == 0) accept = (ntracks == 1 &&  xIntercept > 0);

	  if (denominatorType == 1) accept = (ntracks == 1 &&  xIntercept > 0 
					    && qfine1 > 15 && qfine1 < 500);

	  if (denominatorType == 2) accept = (ntracks == 1 &&  xIntercept > 0 
					    && qfine1 > 15 && qfine1 < 500 
					    &&  energy[int(chRef)] > 21 && qfine[int(chRef)]>15 && qfine[int(chRef)]<500);

	  if (denominatorType == 3) accept = (ntracks == 1 &&  xIntercept > 0 && qfine1 > 15 && qfine1 < 500 
					    &&  energy[int(chRef)] > 21 && qfine[int(chRef)]>15 && qfine[int(chRef)]<500
					    && xIntercept > 16 && xIntercept < 20);

	  if (denominatorType == 4) accept = (ntracks == 1 &&  xIntercept > 0 
					    &&  energy[int(chRef)] > 21 && qfine[int(chRef)]>15 && qfine[int(chRef)]<500);
	  if (!accept) continue;


	  // eff of only energy cut
	  int weight = (energy1>cut_energyAcc[chID][step1]);
	  p_eff_vs_X_effCutEnergy[label] -> Fill( xIntercept,weight );
	  p2_eff_vs_XY_effCutEnergy[label] -> Fill( xIntercept,yIntercept,weight );

	  // eff of energy+qfine cut
          weight = (energy1>cut_energyAcc[chID][step1] && qfine1 > 15 && qfine1 < 500);
          p_eff_vs_X_effCutEnergyQfine[label] -> Fill( xIntercept,weight );
          p2_eff_vs_XY_effCutEnergyQfine[label] -> Fill( xIntercept,yIntercept,weight );

	  // eff of energy+qfine+tot cut
          weight = (energy1>cut_energyAcc[chID][step1] && qfine1 > 15 && qfine1 < 500 && tot1>180);
          p_eff_vs_X_effCutEnergyQfineTot[label] -> Fill( xIntercept,weight );
          p2_eff_vs_XY_effCutEnergyQfineTot[label] -> Fill( xIntercept,yIntercept,weight );
	  

	}

    }

  std::cout << std::endl;
    
  outFile -> Write();
 
}
