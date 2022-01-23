#include "interface/TOFHIRThresholdZero.h"
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
  

  if( argc < 2 )
    {
      std::cout << ">>> drawPulseShape::usage:   " << argv[0] << " configFile.cfg" << std::endl;
      return -1;
    }
  
  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  
  //  int run = opts.GetOpt<int>("Input.run");
  std::string runs = opts.GetOpt<std::string>("Input.runs");
  std::string tofhirVersion = opts.GetOpt<std::string>("Input.tofhirVersion");
  std::cout << tofhirVersion <<std::endl;
  std::string ithMode = opts.GetOpt<std::string>("Input.ithMode");
  float frequency = opts.GetOpt<float>("Input.frequency");

  std::string discCalibrationFile = opts.GetOpt<std::string>("Input.discCalibration");
  TOFHIRThresholdZero thrZero(discCalibrationFile,0);

  int ch1   = opts.GetOpt<float>("Input.ch1");
  int ch2   = opts.GetOpt<float>("Input.ch2");
  int chRef = opts.GetOpt<float>("Input.chRef");
  
  std::vector<int> channels;
  channels.push_back(ch1);
  channels.push_back(ch2);


  float energyMin = opts.GetOpt<float>("Cuts.energyMin");
  float energyMax = opts.GetOpt<float>("Cuts.energyMax");
  float energyMinRef = opts.GetOpt<float>("Cuts.energyMinRef");
  float energyMaxRef = opts.GetOpt<float>("Cuts.energyMaxRef");

  std::string coincidence = opts.GetOpt<std::string>("Cuts.coincidence");
  int ch1Ext   = opts.GetOpt<float>("Cuts.ch1Ext");
  int ch2Ext   = opts.GetOpt<float>("Cuts.ch2Ext");           
  float energyMinExt = opts.GetOpt<float>("Cuts.energyMinExt");
  float energyMaxExt = opts.GetOpt<float>("Cuts.energyMaxExt");
    

  //------------------------------
  // open file and define branches
  //std::string inFileName(Form("/data/TOFHIR2/reco/run%04d_ped_e.root",run));
  //std::string inFileName(Form("/data/TOFHIR2/reco/run%04d_e.root",run));
  //std::cout << "Opening file " << inFileName << std::endl;
  //TFile* inFile = TFile::Open(inFileName.c_str(),"READ");
  //TTree* data = (TTree*)( inFile->Get("data") );

  TChain* data = new TChain("data","data");
  //std::string inFileName = Form("/data/tofhir2/h8/reco/%04d/*ped_e.root",run);
  //data -> Add(inFileName.c_str());

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
    
      for(int run = runMin; run <= runMax; ++run) {
	//std::string inFileName = Form("/data/tofhir2/reco/run%04d*_e.root",run); 
	std::string inFileName = Form("/home/data/mtd/data/tofhir2/reco/run%04d_e.root",run); 
	std::cout << ">>> Adding file " << inFileName << std::endl;
	data -> Add(inFileName.c_str());
      }
    }
  



  float step1, step2;
  int channelIdx[2048];
  std::vector<float>* tot = 0;
  std::vector<float>* energy = 0;
  std::vector<long long>* time = 0;
  data -> SetBranchStatus("*",0);
  data -> SetBranchStatus("step1",     1); data -> SetBranchAddress("step1",          &step1);
  data -> SetBranchStatus("step2",     1); data -> SetBranchAddress("step2",          &step2);
  data -> SetBranchStatus("channelIdx",1); data -> SetBranchAddress("channelIdx", channelIdx);
  data -> SetBranchStatus("tot",       1); data -> SetBranchAddress("tot",              &tot);
  data -> SetBranchStatus("energy",    1); data -> SetBranchAddress("energy",        &energy);
  data -> SetBranchStatus("time",      1); data -> SetBranchAddress("time",            &time);
  
  
  //---------------
  // define outfile
  //TFile* outFile = new TFile(Form("./plots/pulseShape_run%04d.root",run),"RECREATE");
  TFile* outFile = new TFile(Form("./plots_fede/pulseShape_run%s.root",runs.c_str()),"RECREATE");
  //  TFile* outFile = new TFile(Form("/data/Lab5015Analysis/pulseShapes/pulseShape_run%04d.root",run),"RECREATE");
  
  
  //------------------
  // define histograms
  std::map<float, std::map<int,TH1F*> > h1_tot_ch1;
  std::map<float, std::map<int,TH1F*> > h1_energy_ch1;
  std::map<float, std::map<int,TH1F*> > h1_tot_ch2;
  std::map<float, std::map<int,TH1F*> > h1_energy_ch2;
  
  std::map<float, std::map<int,TH1F*> > h1_tot_totSel_ch1;
  std::map<float, std::map<int,TH1F*> > h1_energy_totSel_ch1;
  std::map<float, std::map<int,TH1F*> > h1_time1_wide_ch1;
  std::map<float, std::map<int,TH1F*> > h1_time1_totSel_ch1;
  std::map<float, std::map<int,TH1F*> > h1_time2_totSel_ch1;
  std::map<float, std::map<int,TH1F*> > h1_tot_totSel_ch2;
  std::map<float, std::map<int,TH1F*> > h1_energy_totSel_ch2;
  std::map<float, std::map<int,TH1F*> > h1_time1_wide_ch2;
  std::map<float, std::map<int,TH1F*> > h1_time1_totSel_ch2;
  std::map<float, std::map<int,TH1F*> > h1_time2_totSel_ch2;
  
  std::map<float, std::map<int,TH1F*> > h1_deltaT1_totSel;
  std::map<float, std::map<int,TH1F*> > h1_deltaT2_totSel;
  
  
  //-----------------
  // pre-loop over events  
  int nEntries = data -> GetEntries();
  for(int entry = 0; entry < nEntries; ++entry)
  {
    data -> GetEntry(entry);
    if( entry%10000 == 0 )
    {
      std::cout << ">>>reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    }      
    
    float Vov = step1;
    int ith1 = int(step2/10000.)-1;
    int ith2 = int((step2-10000*(ith1+1))/100.)-1;
    
    int ith = -1;
    if( ithMode.find("ith1") != std::string::npos ) ith = ith1;
    if( ithMode.find("ith2") != std::string::npos ) ith = ith2;
    if( ithMode.find("vth1") != std::string::npos ) ith = ith1;
    if( ithMode.find("vth2") != std::string::npos ) ith = ith2;
    
    // -- coincidence with external bar
    if (coincidence.find("yes") != std::string::npos){
      if( channelIdx[ch1Ext] < 0 ) continue; 
      if( channelIdx[ch2Ext] < 0 ) continue; 
      if( (*tot)[channelIdx[ch1Ext]]/1000. < 0. || (*tot)[channelIdx[ch1Ext]]/1000. > 100. ) continue;
      if( (*tot)[channelIdx[ch2Ext]]/1000. < 0. || (*tot)[channelIdx[ch2Ext]]/1000. > 100. ) continue;
      float energyExt = 0.5 * (  (*energy)[channelIdx[ch1Ext]] + (*energy)[channelIdx[ch2Ext]] );
      if ( energyExt < energyMinExt  || energyExt > energyMaxExt) continue;
    }


    for(int ch :  channels)
    {
      if( channelIdx[ch] < 0 ) continue;
      if( ( thrZero.GetThresholdZero(ch,ithMode) + ith + 1) > 63. ) continue;

      if( ch == ch1 && !h1_time1_wide_ch1[Vov][ith] )
	{
	  h1_time1_wide_ch1[Vov][ith]   = new TH1F(Form("h1_time1_wide_ch1_Vov%.1f_ith%02d",Vov,ith),"",10000,-10000.,10000.);
	}
      if( ch == ch2 && !h1_time1_wide_ch2[Vov][ith] )
      {
        h1_time1_wide_ch2[Vov][ith]   = new TH1F(Form("h1_time1_wide_ch2_Vov%.1f_ith%02d",Vov,ith),"",10000,-10000.,10000.);
      }
      
      if( (*energy)[channelIdx[ch]] < energyMin || (*energy)[channelIdx[ch]] > energyMax ) continue;
      
      if( (*tot)[channelIdx[ch]]/1000. < 0. || (*tot)[channelIdx[ch]]/1000. > 100. ) continue;


      // --  laser
      if (frequency > -1){
      
	long int scale = 1000000000./(frequency);
	
	if( ch == ch1 )
	  {
	    h1_time1_wide_ch1[Vov][ith] -> Fill( ((*time)[channelIdx[ch1]]%scale)/1000. );
	  }
	if( ch == ch2 )
	  {
	    h1_time1_wide_ch2[Vov][ith] -> Fill( ((*time)[channelIdx[ch2]]%scale)/1000. );
	  }
      }


      // -- ref channel
      else {

	if( channelIdx[chRef] < 0 ) continue; 
	if( (*energy)[channelIdx[chRef]] < energyMinRef || (*energy)[channelIdx[chRef]] > energyMaxRef ) continue;     

	if( ch == ch1 ) {
	  h1_time1_wide_ch1[Vov][ith] -> Fill( ((*time)[channelIdx[ch1]]-(*time)[channelIdx[chRef]])/1000.  );
	}
        if( ch == ch2 ) {
	  h1_time1_wide_ch2[Vov][ith] -> Fill( ((*time)[channelIdx[ch2]]-(*time)[channelIdx[chRef]])/1000. );
	}
      }

    }
  }
  
    
  
  std::cout << std::endl;
  
  std::map<float,std::map<int,int> > lowestThr;
  std::map<float,std::map<int,float> > timeOffset;
  for(auto mapIt : h1_time1_wide_ch1)
    {
      float Vov = mapIt.first;
      std::map<int,TH1F*> histos = mapIt.second;
      
      for(auto mapIt2 : histos)
	{
	  int ith = mapIt2.first;
	  if(ith < 10) continue;
	  TH1F* histo = mapIt2.second;
	  histo->Write();

	  std::cout << "===>>> " << Vov << " " << ith << " " << histo->GetMean() << std::endl;
	  if( (lowestThr[Vov][ch1] != 0 && ith < lowestThr[Vov][ch1]) || (lowestThr[Vov][ch1] == 0) )
	    {
	      std::cout << "=========>>> " << Vov << " " << ith << " " << histo->GetMean() << std::endl;

	      timeOffset[Vov][ch1] = histo->GetBinCenter(histo->GetMaximumBin());
	      std::cout << "timeOffset = " << timeOffset[Vov][ch1] <<std::endl;
	      //timeOffset[Vov][ch1] = 230;
	      lowestThr[Vov][ch1] = ith;
	    }
	}
    }
  for(auto mapIt : h1_time1_wide_ch2)
    {
      float Vov = mapIt.first;
      std::map<int,TH1F*> histos = mapIt.second;
      
      for(auto mapIt2 : histos)
	{
	  int ith = mapIt2.first;
	  if(ith < 10) continue;
	  TH1F* histo = mapIt2.second;
	  histo->Write();
	  
	  if( (lowestThr[Vov][ch2] != 0 && ith < lowestThr[Vov][ch2]) || (lowestThr[Vov][ch2] == 0) )
	    {
	      timeOffset[Vov][ch2] = histo->GetBinCenter(histo->GetMaximumBin());
	      //timeOffset[Vov][ch2] = 230;
	      lowestThr[Vov][ch2] = ith;
	    }
	}
    }
  
  
  //-----------------
  // loop over events  
  for(int entry = 0; entry < nEntries; ++entry)
  {
    data -> GetEntry(entry);
    if( entry%10000 == 0 )
    {
      std::cout << ">>>reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    }      
    
    float Vov = step1;
    int ith1 = int(step2/10000.)-1;
    int ith2 = int((step2-10000*(ith1+1))/100.)-1;
    
    int ith = -1;
    if( ithMode.find("ith1") != std::string::npos ) ith = ith1;
    if( ithMode.find("ith2") != std::string::npos ) ith = ith2;
    
    if( ithMode.find("vth1") != std::string::npos ) ith = ith1;
    if( ithMode.find("vth2") != std::string::npos ) ith = ith2;
    

 
    // -- coincidence with external bar
    if (coincidence.find("yes") != std::string::npos){
      if( channelIdx[ch1Ext] < 0 ) continue;
      if( channelIdx[ch2Ext] < 0 ) continue;
      if( (*tot)[channelIdx[ch1Ext]]/1000. < 0. || (*tot)[channelIdx[ch1Ext]]/1000. > 100. ) continue;
      if( (*tot)[channelIdx[ch2Ext]]/1000. < 0. || (*tot)[channelIdx[ch2Ext]]/1000. > 100. ) continue;
      float energyExt = 0.5 * (  (*energy)[channelIdx[ch1Ext]] + (*energy)[channelIdx[ch2Ext]] );
      if ( energyExt < energyMinExt  || energyExt > energyMaxExt) continue;
    }

    
    for(int ch :  channels)
    {
      if( channelIdx[ch] < 0 ) continue;
      
      if( ( thrZero.GetThresholdZero(ch,ithMode) + ith + 1) > 63. ) continue;

      if( ch == ch1 && !h1_tot_ch1[Vov][ith] )
	{
        h1_tot_ch1[Vov][ith]    = new TH1F(Form("h1_tot_ch1_Vov%.1f_ith%02d",Vov,ith),"",15000,-50000.,100000.);
        h1_energy_ch1[Vov][ith] = new TH1F(Form("h1_energy_ch1_Vov%.1f_ith%02d",Vov,ith),"",3000,-0.5,2999.5);
        
        h1_tot_totSel_ch1[Vov][ith]     = new TH1F(Form("h1_tot_totSel_ch1_Vov%.1f_ith%02d",Vov,ith),"",1000,0.,100.);
        h1_energy_totSel_ch1[Vov][ith]  = new TH1F(Form("h1_energy_totSel_ch1_Vov%.1f_ith%02d",Vov,ith),"",3000,-0.5,2999.5);
        h1_time1_totSel_ch1[Vov][ith]   = new TH1F(Form("h1_time1_totSel_ch1_Vov%.1f_ith%02d",Vov,ith),"",5000,timeOffset[Vov][ch1]-50.,timeOffset[Vov][ch1]+50.);
        h1_time2_totSel_ch1[Vov][ith]   = new TH1F(Form("h1_time2_totSel_ch1_Vov%.1f_ith%02d",Vov,ith),"",5000,timeOffset[Vov][ch1]-50.,timeOffset[Vov][ch1]+50.);
        
        h1_deltaT1_totSel[Vov][ith] = new TH1F(Form("h1_deltaT1_totSel_Vov%.1f_ith%02d",Vov,ith),"",10000,-100.,100.);
        h1_deltaT2_totSel[Vov][ith] = new TH1F(Form("h1_deltaT2_totSel_Vov%.1f_ith%02d",Vov,ith),"",10000,-100.,100.);
      }
      if( ch == ch2 && !h1_tot_ch2[Vov][ith] )
      {
        h1_tot_ch2[Vov][ith]    = new TH1F(Form("h1_tot_ch2_Vov%.1f_ith%02d",Vov,ith),"",15000,-50000.,100000.);
        h1_energy_ch2[Vov][ith] = new TH1F(Form("h1_energy_ch2_Vov%.1f_ith%02d",Vov,ith),"",3000,-0.5,2999.5);
        
        h1_tot_totSel_ch2[Vov][ith]     = new TH1F(Form("h1_tot_totSel_ch2_Vov%.1f_ith%02d",Vov,ith),"",1000,0.,100.);
        h1_energy_totSel_ch2[Vov][ith]  = new TH1F(Form("h1_energy_totSel_ch2_Vov%.1f_ith%02d",Vov,ith),"",3000,-0.5,2999.5);
        h1_time1_totSel_ch2[Vov][ith]   = new TH1F(Form("h1_time1_totSel_ch2_Vov%.1f_ith%02d",Vov,ith),"",500,timeOffset[Vov][ch2]-50.,timeOffset[Vov][ch2]+50.);
        h1_time2_totSel_ch2[Vov][ith]   = new TH1F(Form("h1_time2_totSel_ch2_Vov%.1f_ith%02d",Vov,ith),"",500,timeOffset[Vov][ch2]-50.,timeOffset[Vov][ch2]+50.);
      }
      
      if( (*energy)[channelIdx[ch]] < energyMin || (*energy)[channelIdx[ch]] > energyMax ) continue;
      
      if( ch == ch1 )
      {
        h1_tot_ch1[Vov][ith] -> Fill( (*tot)[channelIdx[ch1]]/1000. );
        h1_energy_ch1[Vov][ith] -> Fill( (*energy)[channelIdx[ch1]] );
      }
      if( ch == ch2 )
      {
        h1_tot_ch2[Vov][ith] -> Fill( (*tot)[channelIdx[ch2]]/1000. );
        h1_energy_ch2[Vov][ith] -> Fill( (*energy)[channelIdx[ch2]] );
      }
      
      if( (*tot)[channelIdx[ch]]/1000. < 0. || (*tot)[channelIdx[ch]]/1000. > 100. ) continue;
      

      // -- laser 
      if (frequency > -1) {

	long int scale = 1000000000/(frequency);
      
	if( ch == ch1 )
	  {
	    h1_tot_totSel_ch1[Vov][ith] -> Fill( (*tot)[channelIdx[ch1]]/1000. );
	    h1_energy_totSel_ch1[Vov][ith] -> Fill( (*energy)[channelIdx[ch1]] );
	    h1_time1_totSel_ch1[Vov][ith] -> Fill( ((*time)[channelIdx[ch1]]%scale)/1000. );
	    h1_time2_totSel_ch1[Vov][ith] -> Fill( ((*time)[channelIdx[ch1]]%scale)/1000. + (*tot)[channelIdx[ch1]]/1000. );
	  }
	if( ch == ch2 )
	  {
	    h1_tot_totSel_ch2[Vov][ith] -> Fill( (*tot)[channelIdx[ch2]]/1000. );
	    h1_energy_totSel_ch2[Vov][ith] -> Fill( (*energy)[channelIdx[ch2]] );
	    h1_time1_totSel_ch2[Vov][ith] -> Fill( ((*time)[channelIdx[ch2]]%scale)/1000. );
	    h1_time2_totSel_ch2[Vov][ith] -> Fill( ((*time)[channelIdx[ch2]]%scale)/1000. + (*tot)[channelIdx[ch2]]/1000. );
	  }
      }
      
      // -- ch ref
      else {
       if( channelIdx[chRef] < 0 ) continue;
       if( (*energy)[channelIdx[chRef]] < energyMinRef || (*energy)[channelIdx[chRef]] > energyMaxRef ) continue;     

       if( ch == ch1 )
	 {
	   h1_tot_totSel_ch1[Vov][ith] -> Fill( (*tot)[channelIdx[ch1]]/1000. );
	   h1_energy_totSel_ch1[Vov][ith] -> Fill( (*energy)[channelIdx[ch1]] );
	   h1_time1_totSel_ch1[Vov][ith] -> Fill( ((*time)[channelIdx[ch1]]-(*time)[channelIdx[chRef]])/1000. );
	   h1_time2_totSel_ch1[Vov][ith] -> Fill( ((*time)[channelIdx[ch1]]-(*time)[channelIdx[chRef]])/1000. + (*tot)[channelIdx[ch1]]/1000. );
	 }

       if( ch == ch2 )
	 {
	   h1_tot_totSel_ch2[Vov][ith] -> Fill( (*tot)[channelIdx[ch2]]/1000. );
	   h1_energy_totSel_ch2[Vov][ith] -> Fill( (*energy)[channelIdx[ch2]] );
	   h1_time1_totSel_ch2[Vov][ith] -> Fill( ((*time)[channelIdx[ch2]]-(*time)[channelIdx[chRef]])/1000. );
	   h1_time2_totSel_ch2[Vov][ith] -> Fill( ((*time)[channelIdx[ch2]]-(*time)[channelIdx[chRef]])/1000. + (*tot)[channelIdx[ch2]]/1000. );
	 }

     }
      //std::cout << time[ch1] << " - " << time[ch2] << " - " << scale << " - " << time[ch1] << " - " << (time[ch1]%scale)/1000. << " - " << time[ch2] << " - " << (time[ch2]%scale)/1000. << std::endl;
      
      // h1_deltaT1_totSel[Vov][ith] -> Fill( ((*time)[channelIdx[ch1]]-(*time)[channelIdx[ch2]])/1000. );
      // h1_deltaT2_totSel[Vov][ith] -> Fill( (float((*time)[channelIdx[ch1]]-(*time)[channelIdx[ch2]])+(*tot)[channelIdx[ch1]])/1000. );
    }
  }
  std::cout << std::endl;
  
  
  
  //-----------------
  // draw pulse shape
  float dac_to_uA = -1.;
  if( ithMode == "ith2_3" ) dac_to_uA = 1.250;
  if( ithMode == "ith2_2" ) dac_to_uA = 0.940;
  if( ithMode == "ith2_1" ) dac_to_uA = 0.630;
  if( ithMode == "ith2_0" ) dac_to_uA = 0.313;
  if( ithMode == "ith1_3" ) dac_to_uA = 0.630;
  if( ithMode == "ith1_2" ) dac_to_uA = 0.470;
  if( ithMode == "ith1_1" ) dac_to_uA = 0.313;
  if( ithMode == "ith1_0" ) dac_to_uA = 0.156;
  
  if (tofhirVersion.find("2A")!= std::string::npos ){
    if( ithMode == "vth2" ) dac_to_uA = 8.;
    if( ithMode == "vth1_4" ) dac_to_uA = 4.;
    if( ithMode == "vth1_3" ) dac_to_uA = 2.;
    if( ithMode == "vth1_1" ) dac_to_uA = 1.;
    if( ithMode == "vth1_0" ) dac_to_uA = 0.5;
  }

  std::cout << " vth mode " << ithMode << "   dac_to_mV = " << dac_to_uA <<std::endl;

  
  std::map<float, TGraphErrors*> g_N_ch1;
  std::map<float, TGraphErrors*> g_tot_ch1;
  std::map<float, TGraphErrors*> g_energy_ch1;
  std::map<float, TGraphErrors*> g_N_ch2;
  std::map<float, TGraphErrors*> g_tot_ch2;
  std::map<float, TGraphErrors*> g_energy_ch2;
  
  std::map<float, TGraphErrors*> g_N_totSel_ch1;
  std::map<float, TGraphErrors*> g_tot_totSel_ch1;
  std::map<float, TGraphErrors*> g_energy_totSel_ch1;
  std::map<float, TGraphErrors*> g_ps_totSel_ch1;
  std::map<float, TGraphErrors*> g_N_totSel_ch2;
  std::map<float, TGraphErrors*> g_tot_totSel_ch2;
  std::map<float, TGraphErrors*> g_energy_totSel_ch2;
  std::map<float, TGraphErrors*> g_ps_totSel_ch2;
  
  std::map<float, TGraphErrors*> g_ps_totSel_deltaT;
  
  for(auto mapIt : h1_tot_ch1)
    {
      float Vov = mapIt.first;
      for(auto mapIt2 : mapIt.second)
	{
	  int ith = mapIt2.first;
	  TH1F* histo = mapIt2.second;
	  
	  if( !g_N_ch1[Vov] ) g_N_ch1[Vov] = new TGraphErrors();
	  g_N_ch1[Vov] -> SetPoint(g_N_ch1[Vov]->GetN(),ith,histo->Integral());
          g_N_ch1[Vov] -> SetPointError(g_N_ch1[Vov]->GetN()-1,0,sqrt(histo->Integral()));
	  
	  if( !g_tot_ch1[Vov] ) g_tot_ch1[Vov] = new TGraphErrors();
	  g_tot_ch1[Vov] -> SetPoint(g_tot_ch1[Vov]->GetN(),ith,histo->GetMean());
	  g_tot_ch1[Vov] -> SetPointError(g_tot_ch1[Vov]->GetN()-1,0.,histo->GetRMS());
	  
	  histo -> Write();
	  
	  histo = h1_tot_totSel_ch1[Vov][ith];
	  if( histo->Integral() <= 0. ) continue;
	  
	  if( !g_N_totSel_ch1[Vov] ) g_N_totSel_ch1[Vov] = new TGraphErrors();
	  g_N_totSel_ch1[Vov] -> SetPoint(g_N_totSel_ch1[Vov]->GetN(),ith,histo->Integral());
	  
	  if( !g_tot_totSel_ch1[Vov] ) g_tot_totSel_ch1[Vov] = new TGraphErrors();
	  g_tot_totSel_ch1[Vov] -> SetPoint(g_tot_totSel_ch1[Vov]->GetN(),ith,histo->GetMean());
	  g_tot_totSel_ch1[Vov] -> SetPointError(g_tot_totSel_ch1[Vov]->GetN()-1,0.,histo->GetRMS());
	  
	  histo -> Write();
	}
    }
  for(auto mapIt : h1_tot_ch2)
    {
      float Vov = mapIt.first;
      for(auto mapIt2 : mapIt.second)
	{
	  int ith = mapIt2.first;
	  TH1F* histo = mapIt2.second;
	  
	  if( !g_N_ch2[Vov] ) g_N_ch2[Vov] = new TGraphErrors();
	  g_N_ch2[Vov] -> SetPoint(g_N_ch2[Vov]->GetN(),ith,histo->Integral());
	  
	  if( !g_tot_ch2[Vov] ) g_tot_ch2[Vov] = new TGraphErrors();
	  g_tot_ch2[Vov] -> SetPoint(g_tot_ch2[Vov]->GetN(),ith,histo->GetMean());
	  g_tot_ch2[Vov] -> SetPointError(g_tot_ch2[Vov]->GetN()-1,0.,histo->GetRMS());
	  
	  histo -> Write();
	  
	  histo = h1_tot_totSel_ch2[Vov][ith];
	  if( histo->Integral() <= 0. ) continue;
	  
	  if( !g_N_totSel_ch2[Vov] ) g_N_totSel_ch2[Vov] = new TGraphErrors();
	  g_N_totSel_ch2[Vov] -> SetPoint(g_N_totSel_ch2[Vov]->GetN(),ith,histo->Integral());
	  
	  if( !g_tot_totSel_ch2[Vov] ) g_tot_totSel_ch2[Vov] = new TGraphErrors();
	  g_tot_totSel_ch2[Vov] -> SetPoint(g_tot_totSel_ch2[Vov]->GetN(),ith,histo->GetMean());
	  g_tot_totSel_ch2[Vov] -> SetPointError(g_tot_totSel_ch2[Vov]->GetN()-1,0.,histo->GetRMS());
	  
	  histo -> Write();
	}
    }
  
  for(auto mapIt : h1_energy_ch1)
    {
      float Vov = mapIt.first;
      for(auto mapIt2 : mapIt.second)
	{
	  int ith = mapIt2.first;
	  TH1F* histo = mapIt2.second;
	  if( histo->Integral() <= 0. ) continue;
	  
	  if( !g_energy_ch1[Vov] ) g_energy_ch1[Vov] = new TGraphErrors();
	  g_energy_ch1[Vov] -> SetPoint(g_energy_ch1[Vov]->GetN(),ith,histo->GetMean());
	  g_energy_ch1[Vov] -> SetPointError(g_energy_ch1[Vov]->GetN()-1,0.,histo->GetRMS());
	  
	  histo -> Write();
	  
	  histo = h1_energy_totSel_ch1[Vov][ith];
	  if( histo->Integral() <= 0. ) continue;
      
	  if( !g_energy_totSel_ch1[Vov] ) g_energy_totSel_ch1[Vov] = new TGraphErrors();
	  g_energy_totSel_ch1[Vov] -> SetPoint(g_energy_totSel_ch1[Vov]->GetN(),ith,histo->GetMean());
	  g_energy_totSel_ch1[Vov] -> SetPointError(g_energy_totSel_ch1[Vov]->GetN()-1,0.,histo->GetRMS());
	  
	  histo -> Write();
	}
    }
  for(auto mapIt : h1_energy_ch2)
    {
      float Vov = mapIt.first;
      for(auto mapIt2 : mapIt.second)
	{
	  int ith = mapIt2.first;
	  TH1F* histo = mapIt2.second;
	  if( histo->Integral() <= 0. ) continue;
      
	  if( !g_energy_ch2[Vov] ) g_energy_ch2[Vov] = new TGraphErrors();
	  g_energy_ch2[Vov] -> SetPoint(g_energy_ch2[Vov]->GetN(),ith,histo->GetMean());
	  g_energy_ch2[Vov] -> SetPointError(g_energy_ch2[Vov]->GetN()-1,0.,histo->GetRMS());
	  
	  histo -> Write();
	  
	  histo = h1_energy_totSel_ch2[Vov][ith];
	  if( histo->Integral() <= 0. ) continue;
      
	  if( !g_energy_totSel_ch2[Vov] ) g_energy_totSel_ch2[Vov] = new TGraphErrors();
	  g_energy_totSel_ch2[Vov] -> SetPoint(g_energy_totSel_ch2[Vov]->GetN(),ith,histo->GetMean());
	  g_energy_totSel_ch2[Vov] -> SetPointError(g_energy_totSel_ch2[Vov]->GetN()-1,0.,histo->GetRMS());
	}
    }
  
  // for(auto mapIt : h1_deltaT1_totSel)
  //   {
  //     float Vov = mapIt.first;
  //     for(auto mapIt2 : mapIt.second)
  // 	{
  // 	  int ith = mapIt2.first;
  // 	  TH1F* histo = mapIt2.second;
  // 	  if( histo->Integral() <= 0. ) continue;
	  
  // 	  if( !g_ps_totSel_deltaT[Vov] ) g_ps_totSel_deltaT[Vov] = new TGraphErrors();
  // 	  g_ps_totSel_deltaT[Vov] -> SetPoint(g_ps_totSel_deltaT[Vov]->GetN(),histo->GetMean(),ith*dac_to_uA);
  // 	  g_ps_totSel_deltaT[Vov] -> SetPointError(g_ps_totSel_deltaT[Vov]->GetN()-1,histo->GetMeanError(),0.);
	  
  // 	  histo -> Write();
	  
  // 	  histo = h1_deltaT2_totSel[Vov][ith];
  // 	  if( histo->Integral() <= 0. ) continue;
	  
  // 	  g_ps_totSel_deltaT[Vov] -> SetPoint(g_ps_totSel_deltaT[Vov]->GetN(),histo->GetMean(),ith*dac_to_uA);
  // 	  g_ps_totSel_deltaT[Vov] -> SetPointError(g_ps_totSel_deltaT[Vov]->GetN()-1,histo->GetMeanError(),0.);
	  
  // 	  histo -> Write();
  // 	}
  //   }

  for(auto mapIt : h1_time1_totSel_ch1)
  {
    float Vov = mapIt.first;
    for(auto mapIt2 : mapIt.second)
    {
      int ith = mapIt2.first;
      TH1F* histo = mapIt2.second;
      histo -> Write();

      std::cout << "+++++>>>> " << ith << " " << histo->Integral() << std::endl;
      //if( histo->Integral() <= 0.8*h1_time1_totSel_ch1[Vov][lowestThr[Vov][ch1]]->Integral() ) continue;
      //if( histo->Integral() < 10) continue;

      if( !g_ps_totSel_ch1[Vov] ) g_ps_totSel_ch1[Vov] = new TGraphErrors();
      g_ps_totSel_ch1[Vov] -> SetPoint(g_ps_totSel_ch1[Vov]->GetN(),histo->GetMean()-timeOffset[Vov][ch1],ith*dac_to_uA);
      g_ps_totSel_ch1[Vov] -> SetPointError(g_ps_totSel_ch1[Vov]->GetN()-1,histo->GetMeanError(),0.);
      
    }
    for(auto mapIt2 : mapIt.second)
    {
      int ith = mapIt2.first;
      TH1F* histo = h1_time2_totSel_ch1[Vov][ith];

      std::cout << "+++++++++++>>>> " << ith << " " << histo->Integral() << std::endl;
      //if( histo->Integral() <= 0.8*h1_time2_totSel_ch1[Vov][lowestThr[Vov][ch1]]->Integral() ) continue;
      //if( histo->Integral() < 10) continue;  

      g_ps_totSel_ch1[Vov] -> SetPoint(g_ps_totSel_ch1[Vov]->GetN(),histo->GetMean()-timeOffset[Vov][ch1],ith*dac_to_uA);
      g_ps_totSel_ch1[Vov] -> SetPointError(g_ps_totSel_ch1[Vov]->GetN()-1,histo->GetMeanError(),0.);
      
      histo -> Write();
    }
  }
  for(auto mapIt : h1_time1_totSel_ch2)
  {
    float Vov = mapIt.first;
    for(auto mapIt2 : mapIt.second)
    {
      int ith = mapIt2.first;
      TH1F* histo = mapIt2.second;
      histo -> Write();
      //if( histo->Integral() <= 0.9*h1_time1_totSel_ch2[Vov][lowestThr[Vov][ch2]]->Integral() ) continue;
    //  if( histo->Integral() < 10 ) continue;

      if( !g_ps_totSel_ch2[Vov] ) g_ps_totSel_ch2[Vov] = new TGraphErrors();
      g_ps_totSel_ch2[Vov] -> SetPoint(g_ps_totSel_ch2[Vov]->GetN(),histo->GetMean()-timeOffset[Vov][ch2],ith*dac_to_uA);
      printf("time = %0.20f\t amplitude = %0.20f\n", histo->GetMean()-timeOffset[Vov][ch2],ith*dac_to_uA);
      g_ps_totSel_ch2[Vov] -> SetPointError(g_ps_totSel_ch2[Vov]->GetN()-1,histo->GetMeanError(),0.);
    }
    for(auto mapIt2 : mapIt.second)
    {
      int ith = mapIt2.first;
      TH1F* histo = h1_time2_totSel_ch2[Vov][ith];
      //if( histo->Integral() <= 0.9*h1_time2_totSel_ch2[Vov][lowestThr[Vov][ch2]]->Integral() ) continue;
      //if( histo->Integral() < 10 ) continue; 

      g_ps_totSel_ch2[Vov] -> SetPoint(g_ps_totSel_ch2[Vov]->GetN(),histo->GetMean()-timeOffset[Vov][ch2],ith*dac_to_uA);
      g_ps_totSel_ch2[Vov] -> SetPointError(g_ps_totSel_ch2[Vov]->GetN()-1,histo->GetMeanError(),0.);
      
      histo -> Write();
    }
  }
  
  
  //-----------
  // draw plots
  std::cout << "+++ DRAW PLOTS +++" << std::endl;
  std::string plotDir(Form("/var/www/html/MTDST_CERN_Oct21/CCv2/pulseShapes/run%s/",runs.c_str()));
  system(Form("mkdir -p %s",plotDir.c_str()));
  
  TCanvas* c;
  TH1F* hPad;
  
  
  TF1* f_sigmoid_ch1;
  TF1* f_sigmoid_ch2;
  
  for(auto mapIt : h1_tot_ch1)
    {  
      float Vov = mapIt.first;
      
      c = new TCanvas("c","c");
      //hPad = (TH1F*)( gPad->DrawFrame(-0.5,0.,63.5,12000*frequency/10.) );
      hPad = (TH1F*)( gPad->DrawFrame(-0.5,0.,63.5, g_N_totSel_ch1[Vov] -> GetY()[10]*1.5) );
      hPad -> SetTitle(Form(";%s [DAC]; number of hits",ithMode.c_str()));
      hPad -> Draw();
      g_N_totSel_ch1[Vov] -> SetMarkerColor(kRed);
      g_N_totSel_ch1[Vov] -> SetMarkerStyle(26);
      g_N_totSel_ch1[Vov] -> Draw("P,same");
      if( g_N_totSel_ch2[Vov] ) g_N_totSel_ch2[Vov] -> SetMarkerColor(kBlue);
      if( g_N_totSel_ch2[Vov] ) g_N_totSel_ch2[Vov] -> SetMarkerStyle(32);
      if( g_N_totSel_ch2[Vov] ) g_N_totSel_ch2[Vov] -> Draw("P,same");
      
      f_sigmoid_ch1 = new TF1("f_sigmoid_ch1","[0]*(1-0.5*(1.+TMath::Erf((x-[1])/[2])))",0.,64.);
      f_sigmoid_ch1 -> SetNpx(10000);
      f_sigmoid_ch1 -> SetLineWidth(2);
      f_sigmoid_ch1 -> SetLineColor(kRed);
      int index = 0;
      for (int j = 0; j < g_N_ch1[Vov] -> GetN(); j++ ){
        if ( g_N_ch1[Vov]->GetPointY(j) < g_N_ch1[Vov]->GetPointY(0)/2){
          index = j;
          break;
        } 
      }
      //f_sigmoid_ch1 -> SetParameters(g_N_ch1[Vov]->GetPointY(0),12.,3.);
      f_sigmoid_ch1 -> SetParameters(g_N_ch1[Vov]->GetPointY(0), g_N_ch1[Vov]->GetPointX(index),3.);
      g_N_ch1[Vov] -> Fit(f_sigmoid_ch1,"QRS");
      f_sigmoid_ch1->Draw("same");
      
      TLatex* latex_ch1 = new TLatex(0.40,0.90,Form("amplitude = %.1f #muA",dac_to_uA*f_sigmoid_ch1->GetParameter(1)));
      latex_ch1 -> SetNDC();
      latex_ch1 -> SetTextFont(82);
      latex_ch1 -> SetTextSize(0.04);
      latex_ch1 -> SetTextAlign(11);
      latex_ch1 -> SetTextColor(kRed);
      latex_ch1 -> Draw("same");
      
      f_sigmoid_ch2 = new TF1("f_sigmoid_ch2","[0]*(1-0.5*(1.+TMath::Erf((x-[1])/[2])))",0.,64.);
      f_sigmoid_ch2 -> SetNpx(10000);
      f_sigmoid_ch2 -> SetLineWidth(2);
      f_sigmoid_ch2 -> SetLineColor(kBlue);
      index = 0;

      for (int j = 0; j < g_N_ch2[Vov] -> GetN(); j++ ){
        if ( g_N_ch2[Vov]->GetPointY(j) < g_N_ch2[Vov]->GetPointY(0)/2){
          index = j;
          break;
        } 
      }
      //f_sigmoid_ch2 -> SetParameters(g_N_ch2[Vov]->GetPointY(0),12.,3.);
      f_sigmoid_ch2 -> SetParameters(g_N_ch2[Vov]->GetPointY(0), g_N_ch2[Vov]->GetPointX(index),3.);
      g_N_ch2[Vov] -> Fit(f_sigmoid_ch2,"QRS");
      f_sigmoid_ch2->Draw("same");
      
      TLatex* latex_ch2 = new TLatex(0.40,0.85,Form("amplitude = %.1f #muA",dac_to_uA*f_sigmoid_ch2->GetParameter(1)));
      latex_ch2 -> SetNDC();
      latex_ch2 -> SetTextFont(82);
      latex_ch2 -> SetTextSize(0.04);
      latex_ch2 -> SetTextAlign(11);
      latex_ch2 -> SetTextColor(kBlue);
      latex_ch2 -> Draw("same");
      
      c -> Print(Form("%s/g_N_Vov%.1f_ch1_%d_ch2_%d.png",plotDir.c_str(),Vov, ch1, ch2));
      
      delete c;
    }
  for(auto mapIt : h1_tot_ch1)
  {  
    float Vov = mapIt.first;
    
    c = new TCanvas("c","c");
    hPad = (TH1F*)( gPad->DrawFrame(-0.5,0.,63.5,10.) );
    hPad -> SetTitle(Form(";%s [DAC]; ToT [ns]",ithMode.c_str()));
    hPad -> Draw();
    g_tot_ch1[Vov] -> SetMarkerColor(kRed);
    g_tot_ch1[Vov] -> SetMarkerStyle(22);
    g_tot_ch1[Vov] -> Draw("PL,same");
    if( g_tot_ch2[Vov] ) g_tot_ch2[Vov] -> SetMarkerColor(kBlue);
    if( g_tot_ch2[Vov] ) g_tot_ch2[Vov] -> SetMarkerStyle(23);
    if( g_tot_ch2[Vov] ) g_tot_ch2[Vov] -> Draw("PL,same");
    g_tot_totSel_ch1[Vov] -> SetMarkerColor(kRed-4);
    g_tot_totSel_ch1[Vov] -> SetMarkerStyle(26);
    g_tot_totSel_ch1[Vov] -> Draw("PL,same");
    if( g_tot_totSel_ch2[Vov] ) g_tot_totSel_ch2[Vov] -> SetMarkerColor(kBlue-4);
    if( g_tot_totSel_ch2[Vov] ) g_tot_totSel_ch2[Vov] -> SetMarkerStyle(32);
    if( g_tot_totSel_ch2[Vov] ) g_tot_totSel_ch2[Vov] -> Draw("PL,same");
    c -> Print(Form("%s/g_tot_Vov%.1f_ch1_%d_ch2_%d.png",plotDir.c_str(),Vov, ch1, ch2));
    delete c;
  }
  
  for(auto mapIt : h1_energy_ch1)
  {  
    float Vov = mapIt.first;
    
    c = new TCanvas("c","c");
    hPad = (TH1F*)( gPad->DrawFrame(-0.5,0.,63.5,1023.5) );
    hPad -> SetTitle(Form(";%s [DAC]; energy [ADC]",ithMode.c_str()));
    hPad -> Draw();
    g_energy_ch1[Vov] -> SetMarkerColor(kRed);
    g_energy_ch1[Vov] -> SetMarkerStyle(22);
    g_energy_ch1[Vov] -> Draw("PL,same");
    if( g_energy_ch2[Vov] ) g_energy_ch2[Vov] -> SetMarkerColor(kBlue);
    if( g_energy_ch2[Vov] )g_energy_ch2[Vov] -> SetMarkerStyle(23);
    if( g_energy_ch2[Vov] )g_energy_ch2[Vov] -> Draw("PL,same");
    g_energy_totSel_ch1[Vov] -> SetMarkerColor(kRed-4);
    g_energy_totSel_ch1[Vov] -> SetMarkerStyle(26);
    g_energy_totSel_ch1[Vov] -> Draw("PL,same");
    if( g_energy_totSel_ch2[Vov] ) g_energy_totSel_ch2[Vov] -> SetMarkerColor(kBlue-4);
    if( g_energy_totSel_ch2[Vov] ) g_energy_totSel_ch2[Vov] -> SetMarkerStyle(32);
    if( g_energy_totSel_ch2[Vov] ) g_energy_totSel_ch2[Vov] -> Draw("PL,same");
    c -> Print(Form("%s/g_energy_Vov%.1f_ch1_%d_ch2_%d.png",plotDir.c_str(),Vov, ch1, ch2));
    delete c;
  }
  
  
  // for(auto mapIt : h1_tot_ch1)
  //   {  
  //     float Vov = mapIt.first;
      
  //     c = new TCanvas("c","c");
  //     //hPad = (TH1F*)( gPad->DrawFrame(-2.,0.,20.,40.) );
  //     //hPad = (TH1F*)( gPad->DrawFrame(-2.,0.,20.,100.) );
  //     hPad = (TH1F*)( gPad->DrawFrame(-2.,0.,20.,200.) );
  //     hPad -> SetTitle(Form(";time [ns]; pulse shape [mV]"));
  //     hPad -> Draw();
  //     g_ps_totSel_deltaT[Vov] -> SetMarkerColor(kGray+2);
  //     g_ps_totSel_deltaT[Vov] -> SetLineColor(kGray+2);
  //     g_ps_totSel_deltaT[Vov] -> SetMarkerStyle(22);
  //     g_ps_totSel_deltaT[Vov] -> Draw("P,same");
  //     c -> Print(Form("%s/g_ps_deltaT_Vov%.1f.png",plotDir.c_str(),Vov));
  //     delete c;
  //   }
  
  
  for(auto mapIt : h1_time1_totSel_ch1)
  {  

    float Vov = mapIt.first;
    
    float fitXMin = 0.;
    float fitXMax = 999.;
    
    c = new TCanvas("c","c");
    //hPad = (TH1F*)( gPad->DrawFrame(-2.,0.,20.,40.) );
    //hPad = (TH1F*)( gPad->DrawFrame(-2.,0.,20.,100.) );
    hPad = (TH1F*)( gPad->DrawFrame(-2.,0.,10.,65*dac_to_uA) );
    hPad -> SetTitle(Form(";time [ns]; pulse shape [#muA]"));
    if (tofhirVersion.find("2A")!= std::string::npos ){ hPad -> SetTitle(Form(";time [ns]; pulse shape [mV]"));  }
    hPad -> Draw();
    g_ps_totSel_ch1[Vov] -> SetLineColor(kRed-4);
    g_ps_totSel_ch1[Vov] -> SetMarkerColor(kRed-4);
    g_ps_totSel_ch1[Vov] -> SetMarkerStyle(26);
    g_ps_totSel_ch1[Vov] -> Draw("P,same");
    if( g_ps_totSel_ch2[Vov] ) g_ps_totSel_ch2[Vov] -> SetMarkerColor(kBlue-4);
    if( g_ps_totSel_ch2[Vov] ) g_ps_totSel_ch2[Vov] -> SetLineColor(kBlue-4);
    if( g_ps_totSel_ch2[Vov] ) g_ps_totSel_ch2[Vov] -> SetMarkerStyle(32);
    if( g_ps_totSel_ch2[Vov] ) g_ps_totSel_ch2[Vov] -> Draw("P,same");

    float slewRate = 0.;
    TF1* fitFunc_ch1 = new TF1("fitFunc_ch1","pol1",-5.,10.);
    //-- slew rate max 
    for(int point1 = 0; point1 < g_ps_totSel_ch1[Vov]->GetN()-10; ++point1)
    {
      TGraph* g_temp = new TGraph();
      for(int point2 = point1; point2 < point1+10; ++point2)
      {
        g_temp -> SetPoint(g_temp->GetN(),g_ps_totSel_ch1[Vov]->GetPointX(point2),g_ps_totSel_ch1[Vov]->GetPointY(point2));
      }
      
      TF1* f_temp = new TF1("f_temp","pol1",-10.,100.);
      g_temp -> Fit(f_temp,"QNRS");
      
      if( f_temp->GetParameter(1) > slewRate )
      {
        slewRate = f_temp->GetParameter(1);
        fitFunc_ch1 -> SetParameters(f_temp->GetParameter(0),f_temp->GetParameter(1));
      }
    }

    //-- slew rate at low threshold
    TF1* fitFuncLow_ch1 = new TF1("fitFuncLow_ch1","pol1",-5.,10.);
    TGraph* g_temp = new TGraph();
    for(int point1 = 0; point1 < 4; ++point1){
      g_temp -> SetPoint(g_temp->GetN(),g_ps_totSel_ch1[Vov]->GetPointX(point1),g_ps_totSel_ch1[Vov]->GetPointY(point1)); 
    }
    TF1* f_temp = new TF1("f_temp","pol1", g_temp ->GetPointX(0), g_temp ->GetPointX(3));
    f_temp->SetParameter(1,10.);
    g_temp -> Fit(f_temp,"QNRS"); 
    fitFuncLow_ch1 -> SetParameters(f_temp->GetParameter(0),f_temp->GetParameter(1));
    float slewRate_low = fitFuncLow_ch1->GetParameter(1); 
    std::cout << "Slew rate at low threshold = " << slewRate_low << std::endl;
    for(int point = 0; point < g_ps_totSel_ch1[Vov]->GetN(); ++point)
      if( g_ps_totSel_ch1[Vov]->GetPointY(point) > 4. )
      {
        fitXMin = g_ps_totSel_ch1[Vov]->GetPointX(point);
        break;
      }
    for(int point = 0; point < g_ps_totSel_ch1[Vov]->GetN(); ++point)
      if( g_ps_totSel_ch1[Vov]->GetPointY(point) > 15. )
      {
        fitXMax = g_ps_totSel_ch1[Vov]->GetPointX(point);
        std::cout << fitXMax << std::endl;
        break;
      }
    // TF1* fitFunc_ch1 = new TF1("fitFunc_ch1","pol1",0.,7.);
    fitFunc_ch1 -> SetParameters(0.,250.);
    g_ps_totSel_ch1[Vov] -> Fit(fitFunc_ch1,"QNS+","",fitXMin,fitXMax);
    fitFunc_ch1 -> SetLineColor(kRed-4);
    fitFunc_ch1 -> Draw("same");
    TLatex* latex_ch1 = new TLatex(0.40,0.80,Form("slew rate = %.1f #muA/ns",fitFunc_ch1->GetParameter(1)));
    if (tofhirVersion.find("2A")!= std::string::npos ){
      latex_ch1 = new TLatex(0.40,0.80,Form("slew rate = %.1f mV/ns",fitFunc_ch1->GetParameter(1)));     
    }
    latex_ch1 -> SetNDC();
    latex_ch1 -> SetTextFont(82);
    latex_ch1 -> SetTextSize(0.04);
    latex_ch1 -> SetTextAlign(11);
    latex_ch1 -> SetTextColor(kRed-4);
    latex_ch1 -> Draw("same");

    // TLatex* latex_ch1_low = new TLatex(0.40,0.70,Form("slew rate = %.1f #muA/ns",fitFuncLow_ch1->GetParameter(1)));
    // if (tofhirVersion.find("2A")!= std::string::npos ){
    //   latex_ch1_low = new TLatex(0.40,0.70,Form("slew rate = %.1f mV/ns",fitFuncLow_ch1->GetParameter(1)));     
    // }
    // latex_ch1_low -> SetNDC();
    // latex_ch1_low -> SetTextFont(82);
    // latex_ch1_low -> SetTextSize(0.04);
    // latex_ch1_low -> SetTextAlign(11);
    // latex_ch1_low -> SetTextColor(kRed-4);
    // latex_ch1_low -> Draw("same");
    
    slewRate = 0.;
    TF1* fitFunc_ch2 = new TF1("fitFunc_ch","pol1",-10.,100.);
    for(int point1 = 0; point1 < g_ps_totSel_ch2[Vov]->GetN()-10; ++point1)
    {
      TGraph* g_temp = new TGraph();
      for(int point2 = point1; point2 < point1+10; ++point2)
      {
        g_temp -> SetPoint(g_temp->GetN(),g_ps_totSel_ch2[Vov]->GetPointX(point2),g_ps_totSel_ch2[Vov]->GetPointY(point2));
      }
      
      TF1* f_temp = new TF1("f_temp","pol1",-10.,100.);
      g_temp -> Fit(f_temp,"QNRS");
      
      if( f_temp->GetParameter(1) > slewRate )
      {
        slewRate = f_temp->GetParameter(1);
        fitFunc_ch2 -> SetParameters(f_temp->GetParameter(0),f_temp->GetParameter(1));
      }
    }
    if( g_ps_totSel_ch2[Vov] )
    {
      fitXMin = 0.;
      fitXMax = 999.;
      for(int point = 0; point < g_ps_totSel_ch2[Vov]->GetN(); ++point)
        if( g_ps_totSel_ch2[Vov]->GetPointY(point) > 4. )
        {
          fitXMin = g_ps_totSel_ch2[Vov]->GetPointX(point);
          break;
        }
      for(int point = 0; point < g_ps_totSel_ch2[Vov]->GetN(); ++point)
        if( g_ps_totSel_ch2[Vov]->GetPointY(point) > 15. )
        {
          fitXMax = g_ps_totSel_ch2[Vov]->GetPointX(point);
          break;
        }
    }
    // TF1* fitFunc_ch2 = new TF1("fitFunc_ch2","pol1",0.,7.);
    fitFunc_ch2 -> SetParameters(0.,250.);
    g_ps_totSel_ch2[Vov] -> Fit(fitFunc_ch2,"QNS+","",fitXMin,fitXMax);
    fitFunc_ch2 -> SetLineColor(kBlue-4);
    fitFunc_ch2 -> Draw("same");
    TLatex* latex_ch2 = new TLatex(0.40,0.76,Form("slew rate = %.1f #muA/ns",fitFunc_ch2->GetParameter(1)));
    if (tofhirVersion.find("2A")!= std::string::npos ){
      latex_ch2 = new TLatex(0.40,0.76,Form("slew rate = %.1f mV/ns",fitFunc_ch2->GetParameter(1)));     
    }
    latex_ch2 -> SetNDC();
    latex_ch2 -> SetTextFont(82);
    latex_ch2 -> SetTextSize(0.04);
    latex_ch2 -> SetTextAlign(11);
    latex_ch2 -> SetTextColor(kBlue-4);
    latex_ch2 -> Draw("same");
    
    c -> Print(Form("%s/g_ps_Vov%.1f_ch1_%d_ch2_%d.png", plotDir.c_str(),Vov, ch1, ch2));
    delete c;
  }
  
  
  //-----------
  // save plots
  
  for(auto mapIt : g_ps_totSel_ch1)
  {
    mapIt.second -> Write(Form("g_ps_totSel_ch1_Vov%.1f",mapIt.first));
  }
  for(auto mapIt : g_ps_totSel_ch2)
  {
    mapIt.second -> Write(Form("g_ps_totSel_ch2_Vov%.1f",mapIt.first));
  }

  for (auto mapIt:  g_N_ch1)
  {
    mapIt.second -> Write(Form("g_N_ch1_Vov%.1f",mapIt.first));
  } 

  for (auto mapIt:  g_N_ch2)
  {
    mapIt.second -> Write(Form("g_N_ch2_Vov%.1f",mapIt.first));
  } 
  
  
  //gApplication->Terminate(); 
  
  outFile -> Close();
}
