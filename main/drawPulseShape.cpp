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
  
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  std::string outputDir = opts.GetOpt<std::string>("Output.outputDir");
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  
  float totMin = opts.GetOpt<float>("Cuts.totMin");
  float totMax = opts.GetOpt<float>("Cuts.totMax");
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
  TChain* data = new TChain("data","data");
  
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
	std::string inFileName = Form("%s/run%04d*_e.root",inputDir.c_str(),run); 
	std::cout << ">>> Adding file " << inFileName << std::endl;
	data -> Add(inFileName.c_str());
      }
    }
  



  float step1, step2;
  int channelIdx[128];
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
  TFile* outFile = new TFile(Form("%s/pulseShape_run%s.root",outputDir.c_str(),runs.c_str()),"RECREATE");
  
  
  //------------------
  // define histograms
  std::map<float,int> VovMap;
  
  std::map<int, std::map<float, std::map<int,TH1F*> > > h1_tot;
  std::map<int, std::map<float, std::map<int,TH1F*> > > h1_energy;
  std::map<int, std::map<float, std::map<int,TH1F*> > > h1_tot_totSel;
  std::map<int, std::map<float, std::map<int,TH1F*> > > h1_energy_totSel;
  std::map<int, std::map<float, std::map<int,TH1F*> > > h1_time1_wide;
  std::map<int, std::map<float, std::map<int,TH1F*> > > h1_time1_totSel;
  std::map<int, std::map<float, std::map<int,TH1F*> > > h1_time2_totSel;
  
  
  //---------------------
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
    if (coincidence.find("yes") != std::string::npos)
      {
	if( channelIdx[ch1Ext] < 0 ) continue; 
	if( channelIdx[ch2Ext] < 0 ) continue; 
	if( (*tot)[channelIdx[ch1Ext]]/1000. < totMin || (*tot)[channelIdx[ch1Ext]]/1000. > totMax ) continue;
	if( (*tot)[channelIdx[ch2Ext]]/1000. < totMin || (*tot)[channelIdx[ch2Ext]]/1000. > totMax ) continue;
	float energyExt = 0.5 * (  (*energy)[channelIdx[ch1Ext]] + (*energy)[channelIdx[ch2Ext]] );
	if ( energyExt < energyMinExt  || energyExt > energyMaxExt) continue;
      }
    
    VovMap[Vov] += 1;
    
    int chIt = 0;
    for(int ch :  channels)
      {
	++chIt;
	if( channelIdx[ch] < 0 ) continue;
	if( ( thrZero.GetThresholdZero(ch,ithMode) + ith + 1) > 63. ) continue;
	
	if( !h1_time1_wide[ch][Vov][ith] )
	  {
	    h1_time1_wide[ch][Vov][ith] = new TH1F(Form("h1_time1_wide_ch%d_Vov%.02f_ith%02d",chIt,Vov,ith),"",10000,-10000.,10000.);
	  }
	
	if( (*energy)[channelIdx[ch]] < energyMin || (*energy)[channelIdx[ch]] > energyMax ) continue;
	
	if( (*tot)[channelIdx[ch]]/1000. < totMin || (*tot)[channelIdx[ch]]/1000. > totMax ) continue;
	
	// --  laser
	if( frequency > -1 )
	  {
	    long int scale = 1000000000/(frequency);
	    h1_time1_wide[ch][Vov][ith] -> Fill( ((*time)[channelIdx[ch]]%scale)/1000. );
	  }
	
	
	// -- ref channel
	else
	  {
	    if( channelIdx[chRef] < 0 ) continue; 
	    if( (*energy)[channelIdx[chRef]] < energyMinRef || (*energy)[channelIdx[chRef]] > energyMaxRef ) continue;     
	    
	    h1_time1_wide[ch][Vov][ith] -> Fill( ((*time)[channelIdx[ch]]-(*time)[channelIdx[chRef]])/1000.  );
	  }
	
      } // loop over channels
    
  } // pre-loop over entries
  std::cout << std::endl;
  
  
  //-------------
  // find offsets
  std::map<float,std::map<int,int> > lowestThr;
  std::map<float,std::map<int,float> > timeOffset;
  for(int ch :  channels)
    for(auto mapIt : h1_time1_wide[ch])
      {
	float Vov = mapIt.first;
	std::map<int,TH1F*> histos = mapIt.second;
	
	for(auto mapIt2 : histos)
	  {
	    int ith = mapIt2.first;
	    TH1F* histo = mapIt2.second;
	    histo->Write();
	    
	    std::cout << "===>>> ch " << ch << "   Vov " << Vov << "   ith " << ith << "   t " << histo->GetMean() << std::endl;
	    if( (lowestThr[Vov][ch] != 0 && ith < lowestThr[Vov][ch]) || (lowestThr[Vov][ch] == 0) )
	      {
		std::cout << "===>>>>>> " << ith << " " << histo->GetMean();
		
		timeOffset[Vov][ch] = histo->GetBinCenter(histo->GetMaximumBin());
		lowestThr[Vov][ch] = ith;
		std::cout << "   timeOffset = " << timeOffset[Vov][ch] << std::endl;
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
    if( coincidence.find("yes") != std::string::npos )
      {
	if( channelIdx[ch1Ext] < 0 ) continue;
	if( channelIdx[ch2Ext] < 0 ) continue;
	if( (*tot)[channelIdx[ch1Ext]]/1000. < totMin || (*tot)[channelIdx[ch1Ext]]/1000. > totMax ) continue;
	if( (*tot)[channelIdx[ch2Ext]]/1000. < totMin || (*tot)[channelIdx[ch2Ext]]/1000. > totMax ) continue;
	float energyExt = 0.5 * ( (*energy)[channelIdx[ch1Ext]] + (*energy)[channelIdx[ch2Ext]] );
	if( energyExt < energyMinExt || energyExt > energyMaxExt ) continue;
      }
    
    
    int chIt = 0;
    for(int ch : channels)
      {
	++chIt;
	
	if( channelIdx[ch] < 0 ) continue;
	
	if( ( thrZero.GetThresholdZero(ch,ithMode) + ith + 1) > 63. ) continue;
	
	if( !h1_tot[ch][Vov][ith] )
	  {
	    h1_tot[ch][Vov][ith]    = new TH1F(Form("h1_tot_ch%d_Vov%.02f_ith%02d",chIt,Vov,ith),"",15000,-50000.,100000.);
	    h1_energy[ch][Vov][ith] = new TH1F(Form("h1_energy_ch%d_Vov%.02f_ith%02d",chIt,Vov,ith),"",1000,-0.5,999.5);
	    
	    h1_tot_totSel[ch][Vov][ith]    = new TH1F(Form("h1_tot_totSel_ch%d_Vov%.02f_ith%02d",chIt,Vov,ith),"",1000,0.,100.);
	    h1_energy_totSel[ch][Vov][ith] = new TH1F(Form("h1_energy_totSel_ch%d_Vov%.02f_ith%02d",chIt,Vov,ith),"",1000,-0.5,999.5);
	    h1_time1_totSel[ch][Vov][ith]  = new TH1F(Form("h1_time1_totSel_ch%d_Vov%.02f_ith%02d",chIt,Vov,ith),"",5000,timeOffset[Vov][ch]-50.,timeOffset[Vov][ch]+50.);
	    h1_time2_totSel[ch][Vov][ith]  = new TH1F(Form("h1_time2_totSel_ch%d_Vov%.02f_ith%02d",chIt,Vov,ith),"",5000,timeOffset[Vov][ch]-50.,timeOffset[Vov][ch]+50.);
	  }
	
	if( (*energy)[channelIdx[ch]] < energyMin || (*energy)[channelIdx[ch]] > energyMax ) continue;
	
	h1_tot[ch][Vov][ith] -> Fill( (*tot)[channelIdx[ch]]/1000. );
	h1_energy[ch][Vov][ith] -> Fill( (*energy)[channelIdx[ch]] );
	
	if( (*tot)[channelIdx[ch]]/1000. < 0. || (*tot)[channelIdx[ch]]/1000. > 100. ) continue;
	
	
	// -- laser 
	if( frequency > -1 )
	  {
	    long int scale = 1000000000/(frequency);
	    
	    int bin1 = h1_time1_totSel[ch][Vov][ith] -> Fill( ((*time)[channelIdx[ch]]%scale)/1000. );
	    int bin2 = h1_time2_totSel[ch][Vov][ith] -> Fill( ((*time)[channelIdx[ch]]%scale)/1000. + (*tot)[channelIdx[ch]]/1000. );
	    
	    if( bin1 >= 0 && bin1 < h1_time1_totSel[ch][Vov][ith]->GetNbinsX() && 
		bin2 >= 0 && bin2 < h1_time2_totSel[ch][Vov][ith]->GetNbinsX() )
	      {
		h1_tot_totSel[ch][Vov][ith] -> Fill( (*tot)[channelIdx[ch]]/1000. );
		h1_energy_totSel[ch][Vov][ith] -> Fill( (*energy)[channelIdx[ch]] );
	      }
	  }
	
	// -- ch ref
	else
	  {
	    if( channelIdx[chRef] < 0 ) continue;
	    if( (*energy)[channelIdx[chRef]] < energyMinRef || (*energy)[channelIdx[chRef]] > energyMaxRef ) continue;     
	    
	    int bin1 = h1_time1_totSel[ch][Vov][ith] -> Fill( ((*time)[channelIdx[ch]]-(*time)[channelIdx[chRef]])/1000. );
	    int bin2 = h1_time2_totSel[ch][Vov][ith] -> Fill( ((*time)[channelIdx[ch]]-(*time)[channelIdx[chRef]])/1000. + (*tot)[channelIdx[ch]]/1000. );
	    
	    if( bin1 >= 0 && bin1 < h1_time1_totSel[ch][Vov][ith]->GetNbinsX() && 
		bin2 >= 0 && bin2 < h1_time2_totSel[ch][Vov][ith]->GetNbinsX() )
	      {
		h1_tot_totSel[ch][Vov][ith] -> Fill( (*tot)[channelIdx[ch]]/1000. );
		h1_energy_totSel[ch][Vov][ith] -> Fill( (*energy)[channelIdx[ch]] );
	      }
	  }
	
	//std::cout << time[ch1] << " - " << time[ch2] << " - " << scale << " - " << time[ch1] << " - " << (time[ch1]%scale)/1000. << " - " << time[ch2] << " - " << (time[ch2]%scale)/1000. << std::endl;
	
      } // loop over channels
    
  } // loop over entries
  std::cout << std::endl;
  
  
  //------------
  // fill graphs
  float dac_to_uA = -1.;
  if( ithMode == "ith2_3" ) dac_to_uA = 1.250;
  if( ithMode == "ith2_2" ) dac_to_uA = 0.940;
  if( ithMode == "ith2_1" ) dac_to_uA = 0.630;
  if( ithMode == "ith2_0" ) dac_to_uA = 0.313;
  if( ithMode == "ith1_3" ) dac_to_uA = 0.630;
  if( ithMode == "ith1_2" ) dac_to_uA = 0.470;
  if( ithMode == "ith1_1" ) dac_to_uA = 0.313;
  if( ithMode == "ith1_0" ) dac_to_uA = 0.156;
  
  if( tofhirVersion.find("2A")!= std::string::npos )
    {
      if( ithMode == "vth2" ) dac_to_uA = 8.;
      if( ithMode == "vth1_4" ) dac_to_uA = 4.;
      if( ithMode == "vth1_3" ) dac_to_uA = 2.;
      if( ithMode == "vth1_1" ) dac_to_uA = 1.;
      if( ithMode == "vth1_0" ) dac_to_uA = 0.5;
    }
  
  std::cout << " vth mode " << ithMode << "   dac_to_mV = " << dac_to_uA <<std::endl;
  
  
  std::map<int, std::map<float, TGraphErrors*> > g_tot;
  std::map<int, std::map<float, TGraphErrors*> > g_energy;
  std::map<int, std::map<float, TGraphErrors*> > g_N_totSel;
  std::map<int, std::map<float, TGraphErrors*> > g_tot_totSel;
  std::map<int, std::map<float, TGraphErrors*> > g_energy_totSel;
  std::map<int, std::map<float, TGraphErrors*> > g_ps_totSel;
  std::map<int, std::map<float, TGraphErrors*> > g_SR_totSel;
  
  for(int ch : channels)
    for(auto mapIt : h1_time1_totSel[ch])
      {
	float Vov = mapIt.first;
	for(auto mapIt2 : mapIt.second)
	  {
	    int ith = mapIt2.first;
	    TH1F* histo = mapIt2.second;
	    if( histo->Integral() <= 0. ) continue;
	    
	    if( !g_N_totSel[ch][Vov] ) g_N_totSel[ch][Vov] = new TGraphErrors();
	    g_N_totSel[ch][Vov] -> SetPoint(g_N_totSel[ch][Vov]->GetN(),ith,histo->Integral());
	    
	    histo = h1_tot[ch][Vov][ith];
	    if( histo->Integral() <= 0. ) continue;
	    
	    if( !g_tot[ch][Vov] ) g_tot[ch][Vov] = new TGraphErrors();
	    g_tot[ch][Vov] -> SetPoint(g_tot[ch][Vov]->GetN(),ith,histo->GetMean());
	    g_tot[ch][Vov] -> SetPointError(g_tot[ch][Vov]->GetN()-1,0.,histo->GetRMS());
	    
	    histo -> Write();
	    
	    histo = h1_tot_totSel[ch][Vov][ith];
	    if( histo->Integral() <= 0. ) continue;
	    
	    if( !g_tot_totSel[ch][Vov] ) g_tot_totSel[ch][Vov] = new TGraphErrors();
	    g_tot_totSel[ch][Vov] -> SetPoint(g_tot_totSel[ch][Vov]->GetN(),ith,histo->GetMean());
	    g_tot_totSel[ch][Vov] -> SetPointError(g_tot_totSel[ch][Vov]->GetN()-1,0.,histo->GetRMS());
	    
	    histo -> Write();
	  }
      }
  
  for(int ch : channels)
    for(auto mapIt : h1_energy[ch])
      {
	float Vov = mapIt.first;
	for(auto mapIt2 : mapIt.second)
	  {
	    int ith = mapIt2.first;
	    TH1F* histo = mapIt2.second;
	    if( histo->Integral() <= 0. ) continue;
	    
	    if( !g_energy[ch][Vov] ) g_energy[ch][Vov] = new TGraphErrors();
	    g_energy[ch][Vov] -> SetPoint(g_energy[ch][Vov]->GetN(),ith,histo->GetMean());
	    g_energy[ch][Vov] -> SetPointError(g_energy[ch][Vov]->GetN()-1,0.,histo->GetRMS());
	    
	    histo -> Write();
	    
	    histo = h1_energy_totSel[ch][Vov][ith];
	    if( histo->Integral() <= 0. ) continue;
	    
	    if( !g_energy_totSel[ch][Vov] ) g_energy_totSel[ch][Vov] = new TGraphErrors();
	    g_energy_totSel[ch][Vov] -> SetPoint(g_energy_totSel[ch][Vov]->GetN(),ith,histo->GetMean());
	    g_energy_totSel[ch][Vov] -> SetPointError(g_energy_totSel[ch][Vov]->GetN()-1,0.,histo->GetRMS());
	    
	    histo -> Write();
	  }
      }
  
  for(int ch : channels)
    for(auto mapIt : h1_time1_totSel[ch])
      {
	float Vov = mapIt.first;
	for(auto mapIt2 : mapIt.second)
	  {
	    int ith = mapIt2.first;
	    TH1F* histo = mapIt2.second;
	    histo -> Write();
	    
	    std::cout << "+++>>> ch: " << ch << "   ith: " << ith << "    time 1 integral: " << histo->Integral() << std::endl;
	    if( histo->Integral() <= 0.8*h1_time1_totSel[ch][Vov][lowestThr[Vov][ch]]->Integral() ) continue;
	    if( histo->Integral() < 10 ) continue;
	    
	    if( !g_ps_totSel[ch][Vov] ) g_ps_totSel[ch][Vov] = new TGraphErrors();
	    g_ps_totSel[ch][Vov] -> SetPoint(g_ps_totSel[ch][Vov]->GetN(),histo->GetMean()-timeOffset[Vov][ch],ith*dac_to_uA);
	    g_ps_totSel[ch][Vov] -> SetPointError(g_ps_totSel[ch][Vov]->GetN()-1,histo->GetMeanError(),0.);
	    
	  }
	for(auto mapIt2 : mapIt.second)
	  {
	    int ith = mapIt2.first;
	    TH1F* histo = h1_time2_totSel[ch][Vov][ith];
	    
	    std::cout << "+++>>> ch: " << ch << "   ith: " << ith << "   time 2 integral: " << histo->Integral() << std::endl;
	    if( histo->Integral() <= 0.8*h1_time2_totSel[ch][Vov][lowestThr[Vov][ch]]->Integral() ) continue;
	    if( histo->Integral() < 10 ) continue;  
	    
	    g_ps_totSel[ch][Vov] -> SetPoint(g_ps_totSel[ch][Vov]->GetN(),histo->GetMean()-timeOffset[Vov][ch],ith*dac_to_uA);
	    g_ps_totSel[ch][Vov] -> SetPointError(g_ps_totSel[ch][Vov]->GetN()-1,histo->GetMeanError(),0.);
	    
	    histo -> Write();
	  }
      }
  
  for(int ch : channels)
    for(auto mapIt : VovMap)
      {
	float Vov = mapIt.first;
	
	if( !g_SR_totSel[ch][Vov] ) g_SR_totSel[ch][Vov] = new TGraphErrors();

	for(int point1 = 0; point1 < 0.5*g_ps_totSel[ch][Vov]->GetN()-3; ++point1)
	  {
	    int point2 = point1+3;
	    float SR =
	      ( g_ps_totSel[ch][Vov]->GetPointY(point2) - g_ps_totSel[ch][Vov]->GetPointY(point1) ) / 
	      ( g_ps_totSel[ch][Vov]->GetPointX(point2) - g_ps_totSel[ch][Vov]->GetPointX(point1) );
	    g_SR_totSel[ch][Vov] -> SetPoint(g_SR_totSel[ch][Vov]->GetN(),g_ps_totSel[ch][Vov]->GetPointY(point1)/dac_to_uA,SR);
	  }
      }
  
  
  //-----------
  // draw plots
  plotDir += Form("/run%s/",runs.c_str());
  system(Form("mkdir -p %s",plotDir.c_str()));
  
  TCanvas* c;
  TH1F* hPad;
  TH1F* hPad2;
  
  std::vector<int> colors;
  colors.push_back(kBlack);
  colors.push_back(kRed);
  colors.push_back(kBlue);
  
  for(auto mapIt : VovMap)
    {
      float Vov = mapIt.first;
      
      c = new TCanvas("c","c");
      float yMax = std::max(g_N_totSel[ch1][Vov]->GetPointY(0),g_N_totSel[ch2][Vov]->GetPointY(0));
      hPad = (TH1F*)( gPad->DrawFrame(-0.5,0.,63.5,1.5*yMax) );
      hPad -> SetTitle(Form(";%s [DAC]; number of hits",ithMode.c_str()));
      hPad -> Draw();
      
      int chIt = 0;
      for(int ch : channels)
	{
	  ++chIt;
	  
	  g_N_totSel[ch][Vov] -> SetMarkerColor(colors[chIt]);
	  g_N_totSel[ch][Vov] -> SetLineColor(colors[chIt]);
	  g_N_totSel[ch][Vov] -> SetMarkerStyle(20);
	  g_N_totSel[ch][Vov] -> SetMarkerSize(0.7);
	  g_N_totSel[ch][Vov] -> Draw("P,same");
	  
	  TF1* f_sigmoid = new TF1(Form("f_sigmoid_ch%d",chIt),"[0]*(1-0.5*(1.+TMath::Erf((x-[1])/[2])))",0.,64.);
	  f_sigmoid -> SetNpx(10000);
	  f_sigmoid -> SetLineWidth(2);
	  f_sigmoid -> SetLineColor(colors[chIt]);
	  
	  // find middle point of the sigmoid
	  int index = 0;
	  for(int j = 0; j < g_N_totSel[ch][Vov] -> GetN(); ++j)
	    {
	      if( g_N_totSel[ch][Vov]->GetPointY(j) < g_N_totSel[ch][Vov]->GetPointY(0)/2.)
		{
		  index = j;
		  break;
		} 
	    }
	  f_sigmoid -> SetParameters(g_N_totSel[ch][Vov]->GetPointY(0),g_N_totSel[ch][Vov]->GetPointX(index),3.);
	  g_N_totSel[ch][Vov] -> Fit(f_sigmoid,"QRS");
	  f_sigmoid -> Draw("same");
	  
	  TLatex* latex = new TLatex(0.40,0.90-0.05*chIt,Form("amplitude = %.1f #muA",dac_to_uA*f_sigmoid->GetParameter(1)));
	  latex -> SetNDC();
	  latex -> SetTextFont(82);
	  latex -> SetTextSize(0.04);
	  latex -> SetTextAlign(11);
	  latex -> SetTextColor(colors[chIt]);
	  latex -> Draw("same");
	}
      
      c -> Print(Form("%s/g_N_totSel_Vov%.02f.png",plotDir.c_str(),Vov));
      delete c;
    }
  
  
  for(auto mapIt : VovMap)
    {  
      float Vov = mapIt.first;
      
      c = new TCanvas("c","c");
      hPad = (TH1F*)( gPad->DrawFrame(-0.5,0.,63.5,20.) );
      hPad -> SetTitle(Form(";%s [DAC]; ToT [ns]",ithMode.c_str()));
      hPad -> Draw();
      
      int chIt = 0;
      for(int ch : channels)
	{
	  ++chIt;
	  
	  g_tot[ch][Vov] -> SetMarkerColor(colors[chIt]);
	  g_tot[ch][Vov] -> SetLineColor(colors[chIt]);
	  g_tot[ch][Vov] -> SetLineStyle(2);
	  g_tot[ch][Vov] -> SetMarkerStyle(24);
	  g_tot[ch][Vov] -> SetMarkerSize(0.7);
	  g_tot[ch][Vov] -> Draw("PL,same");
	  
	  g_tot_totSel[ch][Vov] -> SetMarkerColor(colors[chIt]);
	  g_tot_totSel[ch][Vov] -> SetLineColor(colors[chIt]);
	  g_tot_totSel[ch][Vov] -> SetMarkerStyle(20);
	  g_tot_totSel[ch][Vov] -> SetMarkerSize(0.7);
	  g_tot_totSel[ch][Vov] -> Draw("PL,same");
	}
      
      c -> Print(Form("%s/g_tot_Vov%.02f.png",plotDir.c_str(),Vov));
      delete c;
    }
  
  
  for(auto mapIt : VovMap)
    {  
      float Vov = mapIt.first;  
      
      c = new TCanvas("c","c");
      hPad = (TH1F*)( gPad->DrawFrame(-0.5,-10.,63.5,1023.5) );
      hPad -> SetTitle(Form(";%s [DAC]; energy [ADC]",ithMode.c_str()));
      hPad -> Draw();
      
      int chIt = 0;
      for(int ch : channels)
	{
	  ++chIt;
	  
	  g_energy[ch][Vov] -> SetMarkerColor(colors[chIt]);
	  g_energy[ch][Vov] -> SetLineColor(colors[chIt]);
	  g_energy[ch][Vov] -> SetLineStyle(2);
	  g_energy[ch][Vov] -> SetMarkerStyle(24);
	  g_energy[ch][Vov] -> SetMarkerSize(0.7);
	  g_energy[ch][Vov] -> Draw("PL,same");
	  
	  g_energy_totSel[ch][Vov] -> SetMarkerColor(colors[chIt]);
	  g_energy_totSel[ch][Vov] -> SetLineColor(colors[chIt]);
	  g_energy_totSel[ch][Vov] -> SetMarkerStyle(20);
	  g_energy_totSel[ch][Vov] -> SetMarkerSize(0.7);
	  g_energy_totSel[ch][Vov] -> Draw("PL,same");
	}
      
      c -> Print(Form("%s/g_energy_Vov%.02f.png",plotDir.c_str(),Vov));
      delete c;
    }
  
  
  std::map<int,float> slewRateMax;
  for(auto mapIt : VovMap)
    {  
      float Vov = mapIt.first;
      
      float fitXMin = 0.;
      float fitXMax = 999.;
      
      c = new TCanvas("c","c",1400,700);
      c -> Divide(2,1);
      c -> cd(1);
      hPad = (TH1F*)( gPad->DrawFrame(-2.,0.,20.,65*dac_to_uA) );
      hPad -> SetTitle(Form(";time [ns]; pulse shape [#muA]"));
      if( tofhirVersion.find("2A")!= std::string::npos ) hPad -> SetTitle(Form(";time [ns]; pulse shape [mV]"));
      hPad -> Draw();
      c -> cd(2);
      hPad2 = (TH1F*)( gPad->DrawFrame(-2.,0.,2.,65*dac_to_uA) );
      hPad2 -> SetTitle(Form(";time [ns]; pulse shape [#muA]"));
      if( tofhirVersion.find("2A")!= std::string::npos ) hPad2 -> SetTitle(Form(";time [ns]; pulse shape [mV]"));
      hPad2 -> Draw();
      
      int chIt = 0;
      for(int ch : channels)
	{
	  ++chIt;
	  
	  g_ps_totSel[ch][Vov] -> SetLineColor(colors[chIt]);
	  g_ps_totSel[ch][Vov] -> SetMarkerColor(colors[chIt]);
	  g_ps_totSel[ch][Vov] -> SetMarkerStyle(20);
	  g_ps_totSel[ch][Vov] -> SetMarkerSize(0.7);
	  c -> cd(1);
	  g_ps_totSel[ch][Vov] -> Draw("P,same");
	  c -> cd(2);
	  g_ps_totSel[ch][Vov] -> Draw("P,same");
	  
	  // slew rate max 
	  int nPoints_SRMax = 5;
	  TF1* fitFunc_SRMax = new TF1(Form("fitFunc_SRMax_ch%d",chIt),"pol1",-5.,10.);
	  for(int point1 = 0; point1 < g_ps_totSel[ch][Vov]->GetN()-nPoints_SRMax; ++point1)
	    {
	      TGraph* g_temp = new TGraph();
	      for(int point2 = point1; point2 < point1+nPoints_SRMax; ++point2)
		{
		  g_temp -> SetPoint(g_temp->GetN(),g_ps_totSel[ch][Vov]->GetPointX(point2),g_ps_totSel[ch][Vov]->GetPointY(point2));
		}
	      
	      TF1* f_temp = new TF1("f_temp","pol1",-10.,100.);
	      g_temp -> Fit(f_temp,"QNRS");
	      
	      if( f_temp->GetParameter(1) > slewRateMax[ch] )
		{
		  slewRateMax[ch] = f_temp->GetParameter(1);
		  fitFunc_SRMax -> SetParameters(f_temp->GetParameter(0),f_temp->GetParameter(1));
		  //fitFunc_SRMax -> SetRange(g_temp->GetPointX(0),g_temp->GetPointX(g_temp->GetN()-1));
		}
	      delete g_temp;
	      delete f_temp;
	    }
	  
	  // slew rate at low threshold
	  float lowTh = opts.GetOpt<float>("Input.lowTh");
	  int nPoints_SRLow = opts.GetOpt<int>("Input.nPointsLow");
	  int pointLow = -1;
	  float slewRateLow = 0.;
	  TF1* fitFunc_SRLow = new TF1(Form("fitFunc_SRLow_ch%d",chIt),"pol1",-5.,10.);
	  
	  for(int point1 = 0; point1 < 0.5*g_ps_totSel[ch][Vov]->GetN(); ++point1)
	    {
	      if( g_ps_totSel[ch][Vov]->GetPointY(point1) > lowTh )
		{
		  pointLow = point1;
	   	  break;
	   	}
	    }
	  TGraph* g_temp = new TGraph();
	  if( pointLow >= 0 )
	    {
	      for(int point1 = pointLow-0.5*nPoints_SRLow; point1 < pointLow+0.5*nPoints_SRLow; ++point1)
		{
		  g_temp -> SetPoint(g_temp->GetN(),g_ps_totSel[ch][Vov]->GetPointX(point1),g_ps_totSel[ch][Vov]->GetPointY(point1)); 
		}
	      g_temp -> Fit(fitFunc_SRLow,"QNRS"); 
	      fitFunc_SRLow -> SetRange(g_temp->GetPointX(0),g_temp->GetPointX(g_temp->GetN()-1));
	      slewRateLow = fitFunc_SRLow->GetParameter(1); 
	    }
	  
	  // draw
	  fitFunc_SRMax -> SetLineColor(colors[chIt]);
	  fitFunc_SRMax -> SetLineWidth(1);
	  fitFunc_SRMax -> SetLineStyle(2);
	  c -> cd(1);
	  fitFunc_SRMax -> Draw("same");
	  
	  TLatex* latex_SRMax = new TLatex(0.30,0.90-0.05*chIt,Form("slope_{max.} = %.1f #muA/ns",fitFunc_SRMax->GetParameter(1)));
	  if( tofhirVersion.find("2A") != std::string::npos )
	    latex_SRMax = new TLatex(0.30,0.90-0.05*chIt,Form("slope_{max.} = %.1f mV/ns",fitFunc_SRMax->GetParameter(1)));
	  latex_SRMax -> SetNDC();
	  latex_SRMax -> SetTextFont(82);
	  latex_SRMax -> SetTextSize(0.04);
	  latex_SRMax -> SetTextAlign(11);
	  latex_SRMax -> SetTextColor(colors[chIt]);
	  c -> cd(1);
	  latex_SRMax -> Draw("same");
	  
	  fitFunc_SRLow -> SetLineColor(colors[chIt]-2);
	  fitFunc_SRLow -> SetLineWidth(1);
	  fitFunc_SRLow -> SetLineWidth(4);
	  c -> cd(2);
	  fitFunc_SRLow -> Draw("same");
	  
	  TLatex* latex_SRLow = new TLatex(0.30,0.90-0.05*chIt,Form("slope_{%.1f uA} = %.1f #muA/ns",lowTh,fitFunc_SRLow->GetParameter(1)));
	  if( tofhirVersion.find("2A") != std::string::npos )
	    latex_SRLow = new TLatex(0.30,0.90-0.05*chIt,Form("slope_{%.1f mV} = %.1f mV/ns",lowTh,fitFunc_SRLow->GetParameter(1)));
	  latex_SRLow -> SetNDC();
	  latex_SRLow -> SetTextFont(82);
	  latex_SRLow -> SetTextSize(0.04);
	  latex_SRLow -> SetTextAlign(11);
	  latex_SRLow -> SetTextColor(colors[chIt]-2);
	  c -> cd(2);
	  latex_SRLow -> Draw("same");
	}
      
      c -> Print(Form("%s/g_ps_ch1_ch2_Vov%.02f.png",plotDir.c_str(),Vov));
      delete c;
    }

  for(auto mapIt : VovMap)
    {  
      float Vov = mapIt.first;  
      
      c = new TCanvas("c","c");
      hPad = (TH1F*)( gPad->DrawFrame(-0.5,0.,63.5,1.5*std::max(slewRateMax[channels[0]],slewRateMax[channels[1]])) );
      hPad -> SetTitle(Form(";%s [DAC]; slope [#muA/ns]",ithMode.c_str()));
      hPad -> Draw();
      
      int chIt = 0;
      for(int ch : channels)
	{
	  ++chIt;

	  g_SR_totSel[ch][Vov] -> SetMarkerColor(colors[chIt]);
	  g_SR_totSel[ch][Vov] -> SetLineColor(colors[chIt]);
	  g_SR_totSel[ch][Vov] -> SetMarkerStyle(20);
	  g_SR_totSel[ch][Vov] -> SetMarkerSize(0.7);
	  g_SR_totSel[ch][Vov] -> Draw("PL,same");
	}
      
      c -> Print(Form("%s/g_SR_Vov%.02f.png",plotDir.c_str(),Vov));
      delete c;
    }
  
  
  //-----------
  // save plots
  
  int chIt = 0;
  for(auto mapIt : g_ps_totSel)
  {
    ++chIt;
    
    for(auto mapIt2 : g_ps_totSel[mapIt.first])
      mapIt2.second -> Write(Form("g_ps_totSel_ch%d_Vov%.02f",chIt,mapIt2.first));
  }

  chIt = 0;
  for(auto mapIt : g_N_totSel)
  {
    ++chIt;
    
    for(auto mapIt2 : g_N_totSel[mapIt.first])
      mapIt2.second -> Write(Form("g_N_totSel_ch%d_Vov%.02f",chIt,mapIt2.first));
  }



  chIt = 0;
  for(auto mapIt : g_SR_totSel)
  {
    ++chIt;
    
    for(auto mapIt2 : g_SR_totSel[mapIt.first])
      mapIt2.second -> Write(Form("g_SR_totSel_ch%d_Vov%.02f",chIt,mapIt2.first));
  }
  
  outFile -> Close();
}
