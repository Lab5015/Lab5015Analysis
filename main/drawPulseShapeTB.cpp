#include "interface/TOFHIRThresholdZero.h"
#include "interface/SetTDRStyle.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/FitUtils.h"

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
  
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  std::string runs = opts.GetOpt<std::string>("Input.runs");
  std::string ithMode = opts.GetOpt<std::string>("Input.ithMode");

  std::string discCalibrationFile = opts.GetOpt<std::string>("Input.discCalibration");
  TOFHIRThresholdZero thrZero(discCalibrationFile,0);

  int chRef = opts.GetOpt<float>("Input.chRef");

  std::string outName = opts.GetOpt<std::string>("Output.outName");
  std::string outDir  = opts.GetOpt<std::string>("Output.outDir");
  std::string plotDirName  = opts.GetOpt<std::string>("Output.plotDir");

  float energyMinRef = opts.GetOpt<float>("Cuts.energyMinRef");
  float energyMaxRef = opts.GetOpt<float>("Cuts.energyMaxRef");

  // -- read minimum energy for each bar from file
  std::string minEnergiesFileName = opts.GetOpt<std::string>("Cuts.minEnergiesFileName");
  std::map < std::pair<int, float>, float> minE; 
  std::cout << minEnergiesFileName <<std::endl;
  if (minEnergiesFileName != "") {
    std::ifstream minEnergiesFile;
    minEnergiesFile.open(minEnergiesFileName);
    std::string line;
    int bar;
    float ov;
    float value;
    while ( minEnergiesFile.good() ){
      getline(minEnergiesFile, line);
      std::istringstream ss(line);
      ss >> bar >> ov >> value; 
      minE[std::make_pair(bar,ov)] = value; 
      std::cout<< bar <<  "   " << ov << "  " << minE[std::make_pair(bar,ov)] <<std::endl;
    }
  }
  //else{
    //    for(unsigned int iBar = 0; iBar < 16; ++iBar){
    //  for(unsigned int ii = 0; ii < Vov.size(); ++ii){
    //	minE[std::make_pair(iBar, Vov[ii])] =   ;
    //  }
    //}
  //}
  

  std::string coincidence = opts.GetOpt<std::string>("Cuts.coincidence");
  int ch1Ext   = opts.GetOpt<float>("Cuts.ch1Ext");
  int ch2Ext   = opts.GetOpt<float>("Cuts.ch2Ext");           
  float energyMinExt = opts.GetOpt<float>("Cuts.energyMinExt");
  float energyMaxExt = opts.GetOpt<float>("Cuts.energyMaxExt");
    
  //--- define channels (read mapping from the configuration file)
  std::vector<unsigned int> channelMapping = opts.GetOpt<std::vector<unsigned int> >("Channels.channelMapping");

  int chL[16];
  int chR[16];
  
  for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar){
    if(opts.GetOpt<int>("Channels.array")==0){
      chL[iBar] = channelMapping[iBar*2+0];
      chR[iBar] = channelMapping[iBar*2+1];
    }
    if(opts.GetOpt<int>("Channels.array")==1){
      chL[iBar] = channelMapping[iBar*2+0]+64;
      chR[iBar] = channelMapping[iBar*2+1]+64;
    }
    std::cout << "Bar: " << iBar << "   chL: "<< chL[iBar] << "    chR: " <<chR[iBar] <<std::endl;
  }

  //------------------------------
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
	//std::string inFileName = Form("/data/tofhir2/h8/reco/%04d/*_e.root",run); 
	//std::string inFileName = Form("/data1/cmsdaq/tofhir2/h8/reco/%04d/*_e.root",run);
	//std::string inFileName = Form("/eos/cms/store/group/dpg_mtd/comm_mtd/TB/MTDTB_H8_Oct2021/TOFHIR2/h8/reco/%04d/*_e.root",run); 
	//std::string inFileName = Form("/afs/cern.ch/work/m/malberti/MTD/TBatFNALMar2023/Lab5015Analysis/data/run%05d_e.root",run); 
	std::string inFileName = Form("/eos/uscms/store/group/cmstestbeam/2023_03_cmstiming_BTL/TOFHIR/RecoData/run%05d_e.root",run); 
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
  TFile* outFile = new TFile(Form("%s/pulseShape_%s.root", outDir.c_str(), outName.c_str()),"RECREATE");
  
  
  //------------------
  // define histograms

  std::map<int,TH1F*> h1_energyLR;
  std::map<int,TH1F*> h1_energyL;
  std::map<int,TH1F*> h1_energyR;

  std::map<int,TH1F*> h1_totL;
  std::map<int,TH1F*> h1_totR;
  
  std::map<int,TH1F*> h1_time1_wide_chL;
  std::map<int,TH1F*> h1_time1_wide_chR;

  std::map<int,TH1F*> h1_time1_totSel_chL;
  std::map<int,TH1F*> h1_time1_totSel_chR;

  std::map<int,TH1F*> h1_time2_totSel_chL;
  std::map<int,TH1F*> h1_time2_totSel_chR;

  std::map<int,TH2F*> h2_time1_vs_energy_totSel_chL;
  std::map<int,TH2F*> h2_time1_vs_tot_totSel_chL;
  
  std::map<int,bool> acceptEvent; 

  //float fract = 0.90;
  float fract = 0.95;
  int minEntries = 100;

  //-----------------
  // pre-loop over events  
  int nEntries = data -> GetEntries();
  //int nEntries = 1000000;
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
      
      acceptEvent[entry] = true;

      // -- remove showering events on array0,1
      int nActiveBars0 = 0;
      int nActiveBars1 = 0;
      for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar){          
        if ( channelIdx[chL[iBar]] > 0  &&  channelIdx[chR[iBar]] > 0 && (*energy)[channelIdx[chL[iBar]]] > 0 && (*energy)[channelIdx[chR[iBar]]] > 0 )
	  nActiveBars0+=1;
        if ( channelIdx[chL[iBar]+64] > 0  &&  channelIdx[chR[iBar]+64] > 0 && (*energy)[channelIdx[chL[iBar]+64]] > 0 && (*energy)[channelIdx[chR[iBar]+64]] > 0 )
	  nActiveBars1+=1;
      }
      
      int maxActiveBars = 3;
      //if (Vov>4.0) maxActiveBars = 5;
      if (nActiveBars0 > maxActiveBars || nActiveBars1 > maxActiveBars){
	acceptEvent[entry] = false;
	continue;
      }

      // --- keep only events with coincidence with an external bar
      if (coincidence.find("yes") != std::string::npos){
	if( channelIdx[ch1Ext] < 0 ) continue; 
	if( channelIdx[ch2Ext] < 0 ) continue; 
	if( (*tot)[channelIdx[ch1Ext]]/1000. < -10. || (*tot)[channelIdx[ch1Ext]]/1000. > 100. ) continue;
	if( (*tot)[channelIdx[ch2Ext]]/1000. < -10. || (*tot)[channelIdx[ch2Ext]]/1000. > 100. ) continue;
	float energyExt = 0.5 * (  (*energy)[channelIdx[ch1Ext]] + (*energy)[channelIdx[ch2Ext]] );
	if ( energyExt < energyMinExt  || energyExt > energyMaxExt) continue;
      }
      	
      float energyL[16];
      float energyR[16];
      
      float totL[16];
      float totR[16];
      
      long long timeL[16];
      long long timeR[16];
      
      // -- loop over bars in the module
      for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar){
	
	if ( channelIdx[chL[iBar]] < 0  ||  channelIdx[chR[iBar]] < 0 ) continue;

	energyL[iBar]=(*energy)[channelIdx[chL[iBar]]];
	energyR[iBar]=(*energy)[channelIdx[chR[iBar]]];
	totL[iBar]=(*tot)[channelIdx[chL[iBar]]];
	totR[iBar]=(*tot)[channelIdx[chR[iBar]]];
	timeL[iBar]=(*time)[channelIdx[chL[iBar]]];
	timeR[iBar]=(*time)[channelIdx[chR[iBar]]];
	
	if( totL[iBar]/1000 <= -10. || totR[iBar]/1000 <= -10. ) continue;
	if( totL[iBar]/1000 >= 100. ||  totR[iBar]/1000 >= 100.) continue;
	if( ( thrZero.GetThresholdZero(chL[iBar],ithMode) + ith) > 63. ) continue;
	if( ( thrZero.GetThresholdZero(chR[iBar],ithMode) + ith) > 63. ) continue;

	int index( (10000*round(Vov*100.)) + (100*ith) + iBar );

	// -- book histograms if needed
	if (!h1_energyLR[index]){
	  //h1_totL[index] = new TH1F(Form("h1_tot_bar%02dL_Vov%.2f_ith%02d",iBar,Vov,ith),"",1100,-10.,100.);
	  //h1_totR[index] = new TH1F(Form("h1_tot_bar%02dR_Vov%.2f_ith%02d",iBar,Vov,ith),"",1100,-10.,100.);
	  h1_totL[index] = new TH1F(Form("h1_tot_bar%02dL_Vov%.2f_ith%02d",iBar,Vov,ith),"",2000,-10.,10.);
	  h1_totR[index] = new TH1F(Form("h1_tot_bar%02dR_Vov%.2f_ith%02d",iBar,Vov,ith),"",2000,-10.,10.);
	  h1_energyL[index] = new TH1F(Form("h1_energy_bar%02dL_Vov%.2f_ith%02d",iBar,Vov,ith),"",512,0.,1024.);
	  h1_energyR[index] = new TH1F(Form("h1_energy_bar%02dR_Vov%.2f_ith%02d",iBar,Vov,ith),"",512,0.,1024.);
	  h1_energyLR[index] = new TH1F(Form("h1_energy_bar%02dL-R_Vov%.2f_ith%02d",iBar,Vov,ith),"",1024,0.,1024.);
	  h1_time1_wide_chL[index] = new TH1F(Form("h1_time1_wide_bar%02dL_Vov%.2f_ith%02d",iBar,Vov,ith),"",10000,-10000.,10000.);
	  h1_time1_wide_chR[index] = new TH1F(Form("h1_time1_wide_bar%02dR_Vov%.2f_ith%02d",iBar,Vov,ith),"",10000,-10000.,10000.);
	}

	h1_energyL[index]  -> Fill( energyL[iBar] );
	h1_energyR[index]  -> Fill( energyR[iBar] );
	h1_energyLR[index] -> Fill( 0.5*(energyL[iBar] + energyR[iBar]) );
	//h1_totL[index]  -> Fill( totL[iBar]/1000. );
	//h1_totR[index]  -> Fill( totR[iBar]/1000. );
	
	if ( 0.5*(energyL[iBar] + energyR[iBar]) < minE[std::make_pair(iBar,Vov)] ||  0.5*(energyL[iBar] + energyR[iBar]) > 940. ) continue;
    
	// -- ref channel
	if( channelIdx[chRef] < 0 ) continue; 
	if( (*energy)[channelIdx[chRef]] < energyMinRef || (*energy)[channelIdx[chRef]] > energyMaxRef ) continue;     
	
	h1_time1_wide_chL[index] -> Fill( (timeL[iBar] - (*time)[channelIdx[chRef]])/1000.  );
	h1_time1_wide_chR[index] -> Fill( (timeR[iBar] - (*time)[channelIdx[chRef]])/1000.  );
      
      }// end loop over bars
    }// end pre-loop over events

  
  std::cout << std::endl;

  // -- find energy ranges
  std::map<int,float> energyMins;
  std::map<int,float> energyMaxs;
  for(auto mapIt : h1_energyLR){
    int index = mapIt.first;
    TH1F* histo = mapIt.second;
    
    float Vov = int(index/10000)/100.;
    int ith = int ((index - int(10000*100*Vov))/100.);
    int iBar = index - (Vov*10000*100) - (100*ith); 

    histo->GetXaxis()->SetRangeUser(minE[std::make_pair(iBar,Vov)], 940);
    int max =  histo->GetBinCenter(histo->GetMaximumBin());
    TF1 *f_landau = new TF1(Form("f_landau_bar%02d_Vov%.2f_th%02d",iBar,Vov,ith),"[0]*TMath::Landau(x,[1],[2])", 0,1000.);                 
    float xmin = max * 0.65;                                                                                                                                              
    float xmax = std::min(max*2.5, 940.);
    f_landau -> SetRange(xmin,xmax);                                                                                                                               
    f_landau -> SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, max, 0.1*max);                                                       
    f_landau -> SetParLimits(1,0,9999);                                                                                                                            
    f_landau -> SetParLimits(2,0,9999);                                                                                                                            
    histo -> Fit(f_landau,"QRS");                                                                                                                                  
    if ( f_landau->GetParameter(1) > 0 ){                                                                                                                          
      xmin = f_landau->GetParameter(1) - 2.0 * std::abs(f_landau->GetParameter(2));                                                                           
      if (xmin < minE[std::make_pair(iBar, Vov)]) xmin = minE[std::make_pair(iBar, Vov)] ;                                                                                
      xmax = std::min(f_landau->GetParameter(1) * 2.5, 940.);                                                                                                      
      f_landau -> SetRange(xmin, xmax);                                                                                                                            
      f_landau -> SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, f_landau->GetParameter(1), 0.1*f_landau->GetParameter(1));
    }
    histo -> Fit(f_landau,"QRS");  
    histo->GetXaxis()->SetRangeUser(0,1024);
    
    if (  f_landau->GetParameter(1) > minE[std::make_pair(iBar, Vov)] &&                                                                                              
	  (f_landau->GetParameter(1) - 2.0 * std::abs(f_landau->GetParameter(2))) >= minE[std::make_pair(iBar, Vov)] &&
	  (f_landau->GetParameter(1) - 2.0 * std::abs(f_landau->GetParameter(2))) < 940) {     
	  energyMins[index] =  f_landau->GetParameter(1) - 2.0 * std::abs(f_landau->GetParameter(2));
	  //(f_landau->GetParameter(1)*0.80) >= minE[std::make_pair(iBar, Vov)] && (f_landau->GetParameter(1) * 0.80) < 940) {
	  //energyMins[index] =  f_landau->GetParameter(1)*0.80;
	  // energy max
	  energyMaxs[index] = 940; // take full mip spectrum
	  //energyMaxs[index] = std::min(f_landau->GetParameter(1)*2.0, 940.); // select around mip spectrum
    }
    else {
      energyMins[index] = minE[std::make_pair(iBar,Vov)];
      energyMaxs[index] = 940;// take full mip spectrum
    }
    std::cout << Vov << "  th = " << ith  << "   bar = "  << iBar <<  "   minEnergy = " << energyMins[index] <<  "  " <<  f_landau->GetParameter(1)  << "  " << minE[std::make_pair(iBar,Vov)] <<std::endl;
    histo->Write();      
  }


  // -- find time offset
  std::map<int,int>  lowestThrL;
  std::map<int,float>  timeOffsetL;
  std::map<int,int>  lowestThrR;
  std::map<int,float>  timeOffsetR;
  for(auto mapIt : h1_time1_wide_chL)
    {
      int index = mapIt.first;
      TH1F* histo = mapIt.second;
      
      float Vov =  int(index/10000)/100.;
      int ith = int ((index - int(10000*100*Vov))/100.);

      int index2 = index - ith*100;

      histo->Write();
      if ( histo ->GetEntries() < 1000 ) continue;

      //std::cout << "===>>> " << Vov << " " << ith << " " << histo->GetMean() << std::endl;

      if( (lowestThrL[index2] != 0 && ith < lowestThrL[index2]) || (lowestThrL[index2] == 0) )
	{
	  //std::cout << "=========>>> " << Vov << " " << ith << " " << histo->GetMean() << std::endl;
	  timeOffsetL[index2] = histo->GetBinCenter(histo->GetMaximumBin());
	  //std::cout << "chL timeOffset = " << timeOffsetL[index2] <<std::endl;
	  lowestThrL[index2] = ith;
	}
    }


  for(auto mapIt : h1_time1_wide_chR)
    {
      int index = mapIt.first;
      TH1F* histo = mapIt.second;
      
      float Vov =  int(index/10000)/100.;
      int ith = int ((index - int(10000*100*Vov))/100.);
      
      int index2 = index - ith*100;

      histo->Write();
      if ( histo ->GetEntries() < 1000 ) continue;

      //std::cout << "===>>> " << Vov << " " << ith << " " << histo->GetMean() << std::endl;

      if( (lowestThrR[index2] != 0 && ith < lowestThrR[index2]) || (lowestThrR[index2] == 0) )
	{
	  //std::cout << "=========>>> " << Vov << " " << ith << " " << histo->GetMean() << std::endl;
	  timeOffsetR[index2] = histo->GetBinCenter(histo->GetMaximumBin());
	  //std::cout << "chR timeOffset = " << timeOffsetR[index2] <<std::endl;
	  lowestThrR[index2] = ith;
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
    
    if (!acceptEvent[entry]) continue;

    float Vov = roundf(step1*100)/100;
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
      if( (*tot)[channelIdx[ch1Ext]]/1000. < -10. || (*tot)[channelIdx[ch1Ext]]/1000. > 100. ) continue;
      if( (*tot)[channelIdx[ch2Ext]]/1000. < -10. || (*tot)[channelIdx[ch2Ext]]/1000. > 100. ) continue;
      float energyExt = 0.5 * (  (*energy)[channelIdx[ch1Ext]] + (*energy)[channelIdx[ch2Ext]] );
      if ( energyExt < energyMinExt  || energyExt > energyMaxExt) continue;
    }

    float energyL[16];
    float energyR[16];
    
    float totL[16];
    float totR[16];
    
    long long timeL[16];
    long long timeR[16];
    
    // -- loop over bars in the module
    for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar){
      
      if ( channelIdx[chL[iBar]] < 0  ||  channelIdx[chR[iBar]] < 0 ) continue;
      
      //int index( (10000*int(Vov*100.)) + (100*ith) + iBar );
      int index( (10000*round(Vov*100.)) + (100*ith) + iBar );
      int index2 = index - ith*100;

      energyL[iBar]=(*energy)[channelIdx[chL[iBar]]];
      energyR[iBar]=(*energy)[channelIdx[chR[iBar]]];
      totL[iBar]=(*tot)[channelIdx[chL[iBar]]];
      totR[iBar]=(*tot)[channelIdx[chR[iBar]]];
      timeL[iBar]=(*time)[channelIdx[chL[iBar]]];
      timeR[iBar]=(*time)[channelIdx[chR[iBar]]];
      
      if( totL[iBar]/1000 <= -10. || totR[iBar]/1000 <= -10. ) continue;
      if( totL[iBar]/1000 >= 100. ||  totR[iBar]/1000 >= 100.) continue;
    
      if( ( thrZero.GetThresholdZero(chL[iBar],ithMode) + ith) > 63. ) continue;
      if( ( thrZero.GetThresholdZero(chR[iBar],ithMode) + ith) > 63. ) continue;

      if ( 0.5*(energyL[iBar] + energyR[iBar]) < energyMins[index] ||  0.5*(energyL[iBar] + energyR[iBar]) > energyMaxs[index] ) continue;
      	
      // -- book histograms if needed
      if (!h1_time1_totSel_chL[index]){
	h1_time1_totSel_chL[index] = new TH1F(Form("h1_time1_totSel_bar%02dL_Vov%.2f_ith%02d",iBar,Vov,ith),"",5000,timeOffsetL[index2]-50.,timeOffsetL[index2]+50.);
	h1_time1_totSel_chR[index] = new TH1F(Form("h1_time1_totSel_bar%02dR_Vov%.2f_ith%02d",iBar,Vov,ith),"",5000,timeOffsetR[index2]-50.,timeOffsetR[index2]+50.);
	h1_time2_totSel_chL[index] = new TH1F(Form("h1_time2_totSel_bar%02dL_Vov%.2f_ith%02d",iBar,Vov,ith),"",5000,timeOffsetL[index2]-50.,timeOffsetL[index2]+50.);
	h1_time2_totSel_chR[index] = new TH1F(Form("h1_time2_totSel_bar%02dR_Vov%.2f_ith%02d",iBar,Vov,ith),"",5000,timeOffsetR[index2]-50.,timeOffsetR[index2]+50.);
	h2_time1_vs_energy_totSel_chL[index] = new TH2F(Form("h2_time1_vs_energy_totSel_bar%02dL_Vov%.2f_ith%02d",iBar,Vov,ith),"",1024,0,1024,5000,timeOffsetL[index2]-50.,timeOffsetL[index2]+50.);
	h2_time1_vs_tot_totSel_chL[index] = new TH2F(Form("h2_time1_vs_tot_totSel_bar%02dL_Vov%.2f_ith%02d",iBar,Vov,ith),"",100,0,100,5000,timeOffsetL[index2]-50.,timeOffsetL[index2]+50.);
      }
    
      // -- ref channel
      if( channelIdx[chRef] < 0 ) continue; 
      if( (*energy)[channelIdx[chRef]] < energyMinRef || (*energy)[channelIdx[chRef]] > energyMaxRef ) continue;     
           

      h1_time1_totSel_chL[index] -> Fill( (timeL[iBar] - (*time)[channelIdx[chRef]])/1000.  );
      h1_time1_totSel_chR[index] -> Fill( (timeR[iBar] - (*time)[channelIdx[chRef]])/1000.  );
      h1_time2_totSel_chL[index] -> Fill( (timeL[iBar] - (*time)[channelIdx[chRef]])/1000.  + totL[iBar]/1000. );
      h1_time2_totSel_chR[index] -> Fill( (timeR[iBar] - (*time)[channelIdx[chRef]])/1000.  + totR[iBar]/1000. );
      h2_time1_vs_energy_totSel_chL[index]->Fill(0.5*(energyL[iBar] + energyR[iBar]) , (timeL[iBar] - (*time)[channelIdx[chRef]])/1000. );
      h2_time1_vs_tot_totSel_chL[index]->Fill( totL[iBar]/1000. , (timeL[iBar] - (*time)[channelIdx[chRef]])/1000. );

      h1_totL[index]  -> Fill( totL[iBar]/1000. );
      h1_totR[index]  -> Fill( totR[iBar]/1000. );

    }// end loop over bars
  }// end loop over events

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
  
  std::cout << " vth mode " << ithMode << "   dac_to_uA = " << dac_to_uA <<std::endl;
  
  std::map<int, TGraphErrors*> g_N_L;
  std::map<int, TGraphErrors*> g_totL;
  std::map<int, TGraphErrors*> g_energyL;
  std::map<int, TGraphErrors*> g_N_R;
  std::map<int, TGraphErrors*> g_totR;
  std::map<int, TGraphErrors*> g_energyR;
  std::map<int, TGraphErrors*> g_energyLR;
  
  std::map<int, TGraphErrors*> g_pulseShapeL;
  std::map<int, TGraphErrors*> g_pulseShapeR;

  
  for(auto mapIt : h1_totL)
    {
      int index = mapIt.first;
      TH1F* histo = mapIt.second;
    
      float Vov =  int(index/10000)/100.;
      int ith = int ((index - int(10000*100*Vov))/100.);
      
      int index2 = index - ith*100;

      if( !g_N_L[index2] ) g_N_L[index2] = new TGraphErrors();
      g_N_L[index2] -> SetPoint(g_N_L[index2]->GetN(),ith,histo->Integral());
      g_N_L[index2] -> SetPointError(g_N_L[index2]->GetN()-1,0,sqrt(histo->Integral()));
	  
      if( !g_totL[index2] ) g_totL[index2] = new TGraphErrors();
      g_totL[index2] -> SetPoint(g_totL[index2]->GetN(),ith,histo->GetMean());
      g_totL[index2] -> SetPointError(g_totL[index2]->GetN()-1,0.,histo->GetRMS());
	  
      histo -> Write();
    }


  
  for(auto mapIt : h1_totR)
    {
      int index = mapIt.first;
      TH1F* histo = mapIt.second;
        
      float Vov =  int(index/10000)/100.;
      int ith = int ((index - int(10000*100*Vov))/100.);
      int index2 = index - ith*100;

      if( !g_N_R[index2] ) g_N_R[index2] = new TGraphErrors();
      g_N_R[index2] -> SetPoint(g_N_R[index2]->GetN(),ith,histo->Integral());
      g_N_R[index2] -> SetPointError(g_N_R[index2]->GetN()-1,0,sqrt(histo->Integral()));
	  
      if( !g_totR[index2] ) g_totR[index2] = new TGraphErrors();
      g_totR[index2] -> SetPoint(g_totR[index2]->GetN(),ith,histo->GetMean());
      g_totR[index2] -> SetPointError(g_totR[index2]->GetN()-1,0.,histo->GetRMS());
	  
      histo -> Write();
    }


  for(auto mapIt : h1_energyL)
    {
      int index = mapIt.first;
      TH1F* histo = mapIt.second;
      
      float Vov =  int(index/10000)/100.;
      int ith = int ((index - int(10000*100*Vov))/100.);
      int iBar = index - (Vov*10000*100) - (100*ith);  
      int index2 = index - ith*100;

      int max =  histo->GetBinCenter(histo->GetMaximumBin());
      TF1 *f_landau = new TF1(Form("f_landauL_bar%02dL_Vov%.2f_ith%02d",iBar,Vov,ith),"[0]*TMath::Landau(x,[1],[2])", 0,1000.);                 
      f_landau -> SetRange(max*0.75,max*2.5);
      f_landau -> SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, max, 0.1*max);                                                   
      f_landau -> SetParLimits(1,0,9999);                                                                                                                        
      f_landau -> SetParLimits(2,0,9999);                                                                                                                        
      histo -> Fit(f_landau,"QRS+");
            
      if( !g_energyL[index2] ) g_energyL[index2] = new TGraphErrors();
      g_energyL[index2] -> SetPoint(g_energyL[index2]->GetN(),ith, f_landau->GetParameter(1));
      g_energyL[index2] -> SetPointError(g_energyL[index2]->GetN()-1,0., f_landau->GetParError(1));
	  
      histo -> Write();
    }



  for(auto mapIt : h1_energyR)
    {
      int index = mapIt.first;
      TH1F* histo = mapIt.second;
      
      float Vov =  int(index/10000)/100.;
      int ith = int ((index - int(10000*100*Vov))/100.);
      int iBar = index - (Vov*10000*100) - (100*ith);  
      int index2 = index - ith*100;
      
      int max =  histo->GetBinCenter(histo->GetMaximumBin());
      TF1 *f_landau = new TF1(Form("f_landauL_bar%02dR_Vov%.2f_ith%02d",iBar,Vov,ith),"[0]*TMath::Landau(x,[1],[2])", 0,1000.);                 
      f_landau -> SetRange(max*0.75,max*2.5);
      f_landau -> SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, max, 0.1*max);                                                   
      f_landau -> SetParLimits(1,0,9999);                                                                                                                        
      f_landau -> SetParLimits(2,0,9999);                                                                                                                        
      histo -> Fit(f_landau,"QRS+");
            
      if( !g_energyR[index2] ) g_energyR[index2] = new TGraphErrors();
      g_energyR[index2] -> SetPoint(g_energyR[index2]->GetN(),ith, f_landau->GetParameter(1));
      g_energyR[index2] -> SetPointError(g_energyR[index2]->GetN()-1,0., f_landau->GetParError(1));
	  
      histo -> Write();
    }

  // -- pulse shape chL
  float* vals = new float[6]; 

  for(auto mapIt : h1_time1_totSel_chL)
    {
      int index = mapIt.first;
      TH1F* histo = mapIt.second;
      
      float Vov =  int(index/10000)/100.;
      int ith = int ((index - int(10000*100*Vov))/100.);
      int iBar = index - (Vov*10000*100) - (100*ith);  
      int index2 = index - ith*100;

      if( histo->Integral() < minEntries) continue;
      
      if( !g_pulseShapeL[index2] ) g_pulseShapeL[index2] = new TGraphErrors();

      // -- find smallest interval containing 90% of the events
      FindSmallestInterval(vals,histo,fract);
      float mean = vals[0];
      float meanErr = vals[1];
      histo->GetXaxis()->SetRangeUser(vals[4],vals[5]);
      g_pulseShapeL[index2] -> SetPoint(g_pulseShapeL[index2]->GetN(),mean-timeOffsetL[index2],ith*dac_to_uA);
      g_pulseShapeL[index2] -> SetPointError(g_pulseShapeL[index2]->GetN()-1,meanErr,0.);
            
      histo -> Write();
      h2_time1_vs_energy_totSel_chL[index]->Write();
      h2_time1_vs_tot_totSel_chL[index]->Write();
    }

  for(auto mapIt : h1_time2_totSel_chL)
    {
      int index = mapIt.first;
      TH1F* histo = mapIt.second;
      
      float Vov =  int(index/10000)/100.;
      int ith = int ((index - int(10000*100*Vov))/100.);
      int iBar = index - (Vov*10000*100) - (100*ith);
      int index2 = index - ith*100;
 
      if( histo->Integral() < minEntries) continue;

      // -- find smallest interval containing 90% of the events
      FindSmallestInterval(vals,histo,fract);
      float mean = vals[0];
      float meanErr = vals[1];
      histo->GetXaxis()->SetRangeUser(vals[4],vals[5]);
      //g_pulseShapeL[index2] -> SetPoint(g_pulseShapeL[index2]->GetN(),mean-timeOffsetL[index2],ith*dac_to_uA);
      //g_pulseShapeL[index2] -> SetPointError(g_pulseShapeL[index2]->GetN()-1,meanErr,0.);
            
      histo -> Write();
    }



  // -- pulse shape chR
  //
  for(auto mapIt : h1_time1_totSel_chR)
    {
      int index = mapIt.first;
      TH1F* histo = mapIt.second;
      
      float Vov =  int(index/10000)/100.;
      int ith = int ((index - int(10000*100*Vov))/100.);
      int iBar = index - (Vov*10000*100) - (100*ith);
      int index2 = index - ith*100;
      
      if( histo->Integral() < minEntries) continue;
            
      if( !g_pulseShapeR[index2] ) g_pulseShapeR[index2] = new TGraphErrors();
      
      // -- find smallest interval containing 90% of the events
      FindSmallestInterval(vals,histo,fract);
      float mean = vals[0];
      float meanErr = vals[1];
      histo->GetXaxis()->SetRangeUser(vals[4],vals[5]); 
      g_pulseShapeR[index2] -> SetPoint(g_pulseShapeR[index2]->GetN(),mean-timeOffsetR[index2],ith*dac_to_uA);
      g_pulseShapeR[index2] -> SetPointError(g_pulseShapeR[index2]->GetN()-1,meanErr,0.);

      histo -> Write();
    }

  for(auto mapIt : h1_time2_totSel_chR)
    {
      int index = mapIt.first;
      TH1F* histo = mapIt.second;
      
      float Vov =  int(index/10000)/100.;
      int ith = int ((index - int(10000*100*Vov))/100.);
      int iBar = index - (Vov*10000*100) - (100*ith);
      int index2 = index - ith*100; 

      if( histo->Integral() < minEntries) continue;
      
      // -- find mean in the smallest interval containing 90% of the events
      FindSmallestInterval(vals,histo,fract);
      float mean = vals[0];
      float meanErr = vals[1];
      histo->GetXaxis()->SetRangeUser(vals[4],vals[5]);
      //g_pulseShapeR[index2] -> SetPoint(g_pulseShapeR[index2]->GetN(),mean-timeOffsetR[index2],ith*dac_to_uA);
      //g_pulseShapeR[index2] -> SetPointError(g_pulseShapeR[index2]->GetN()-1,meanErr,0.);      
            
      histo -> Write();
    }

  //-----------
  // draw plots
  std::cout << "Plotting..."<<std::endl;
  //std::string plotDir(Form("/var/www/html/TOFHIR2X/MTDTB_CERN_June22/pulseShapes/%s",outName.c_str()));
  //std::string plotDir(Form("/var/www/html/TOFHIR2B/MTDTB_CERN_June22/pulseShapes/%s",outName.c_str()));
  //std::string plotDir(Form("/eos/user/m/malberti/www/MTD/TOFHIR2X/MTDTB_FNAL_Mar23/pulseShapes/%s",outName.c_str()));  
  std::string plotDir(Form("%s/%s", plotDirName.c_str(),outName.c_str()));  
  system(Form("mkdir -p %s",plotDir.c_str()));
  


  TCanvas* c;
  TH1F* hPad;

  for(auto mapIt : h1_energyLR)
    {
      int index = mapIt.first;
      TH1F* histo = mapIt.second;
      
      float Vov =  int(index/10000)/100.;
      int ith = int ((index - int(10000*100*Vov))/100.);
      int iBar = index - (Vov*10000*100) - (100*ith);  
      int index2 = index - ith*100;

      TF1 *f_landau = histo->GetFunction(Form("f_landau_bar%02d_Vov%.2f_th%02d",iBar,Vov,ith));
      if (!f_landau) continue;
            
      if( !g_energyLR[index2] ) g_energyLR[index2] = new TGraphErrors();
      g_energyLR[index2] -> SetPoint(g_energyLR[index2]->GetN(),ith, f_landau->GetParameter(1));
      g_energyLR[index2] -> SetPointError(g_energyLR[index2]->GetN()-1,0., f_landau->GetParError(1));

      TCanvas *c = new TCanvas("c","c");
      c->SetLogy();
      TH1F *hPad = (TH1F*)( gPad->DrawFrame(0.,0.1, 1024,histo->GetMaximum()*10.) );
      hPad -> GetXaxis()->SetTitle("energy [ADC]");
      hPad -> GetYaxis()->SetTitle("events");
      hPad -> Draw();
      histo->Draw("same");
      float yval = std::max(10., histo->GetBinContent(histo->FindBin(energyMins[index])));
      TLine* line1 = new TLine(energyMins[index],0.,energyMins[index], yval);
      line1 -> SetLineWidth(1);
      line1 -> SetLineStyle(7);
      line1 -> Draw("same");
      TLine* line2 = new TLine(energyMaxs[index],0.,energyMaxs[index], yval);
      line2 -> SetLineWidth(1);
      line2 -> SetLineStyle(7);
      line2 -> Draw("same");
      c -> Print(Form("%s/c_energy_bar%02dLR_Vov%.2f_th%02d.png",plotDir.c_str(),iBar,Vov,ith));
      c -> Print(Form("%s/c_energy_bar%02dLR_Vov%.2f_th%02d.pdf",plotDir.c_str(),iBar,Vov,ith));
      delete c;
    }
  

  
  std::map<int,TF1*> f_sigmoid_L;
  std::map<int,TF1*> f_sigmoid_R;
  
  for(auto mapIt : g_N_L)
    {  
      int index2 = mapIt.first;
      float Vov =  int(index2/10000)/100.;
      int iBar = index2 - (Vov*10000*100);
            
      c = new TCanvas("c","c");
      hPad = (TH1F*)( gPad->DrawFrame(-0.5,0.,63.5, g_N_L[index2] -> GetY()[0]*1.5) );
      hPad -> SetTitle(Form(";%s [DAC]; number of hits",ithMode.c_str()));
      hPad -> Draw();
      g_N_L[index2] -> SetMarkerColor(kRed);
      g_N_L[index2] -> SetMarkerStyle(20);
      g_N_L[index2] -> Draw("P,same");
      if( g_N_R[index2] ) g_N_R[index2] -> SetMarkerColor(kBlue);
      if( g_N_R[index2] ) g_N_R[index2] -> SetMarkerStyle(20);
      if( g_N_R[index2] ) g_N_R[index2] -> Draw("P,same");
            
      f_sigmoid_L[index2] = new TF1("f_sigmoid_L","[0]*(1-0.5*(1.+TMath::Erf((x-[1])/[2])))",0.,64.);
      f_sigmoid_L[index2] -> SetNpx(10000);
      f_sigmoid_L[index2] -> SetLineWidth(2);
      f_sigmoid_L[index2] -> SetLineColor(kRed);
      int myindex = 0;
      for (int j = 0; j < g_N_L[index2] -> GetN(); j++ ){
        if ( g_N_L[index2]->GetPointY(j) < g_N_L[index2]->GetPointY(0)/2){
          myindex = j;
          break;
        } 
      }
      f_sigmoid_L[index2] -> SetParameters(g_N_L[index2]->GetPointY(0), g_N_L[index2]->GetPointX(myindex),3.);
      g_N_L[index2] -> Fit(f_sigmoid_L[index2],"QRS");
      f_sigmoid_L[index2]->Draw("same");
      
      TLatex* latexL = new TLatex(0.40,0.90,Form("amplitude = %.1f #muA",dac_to_uA*f_sigmoid_L[index2]->GetParameter(1)));
      latexL -> SetNDC();
      latexL -> SetTextFont(82);
      latexL -> SetTextSize(0.04);
      latexL -> SetTextAlign(11);
      latexL -> SetTextColor(kRed);
      latexL -> Draw("same");
      
      
      f_sigmoid_R[index2] = new TF1("f_sigmoid_R","[0]*(1-0.5*(1.+TMath::Erf((x-[1])/[2])))",0.,64.);
      f_sigmoid_R[index2] -> SetNpx(10000);
      f_sigmoid_R[index2] -> SetLineWidth(2);
      f_sigmoid_R[index2] -> SetLineColor(kBlue);
      myindex = 0;
      for (int j = 0; j < g_N_R[index2] -> GetN(); j++ ){
        if ( g_N_R[index2]->GetPointY(j) < g_N_R[index2]->GetPointY(0)/2){
          myindex = j;
          break;
        } 
      }
      f_sigmoid_R[index2] -> SetParameters(g_N_R[index2]->GetPointY(0), g_N_R[index2]->GetPointX(myindex),3.);
      g_N_R[index2] -> Fit(f_sigmoid_R[index2],"QRS");
      f_sigmoid_R[index2]->Draw("same");
      
      TLatex* latexR = new TLatex(0.40,0.85,Form("amplitude = %.1f #muA",dac_to_uA*f_sigmoid_R[index2]->GetParameter(1)));
      latexR -> SetNDC();
      latexR -> SetTextFont(82);
      latexR -> SetTextSize(0.04);
      latexR -> SetTextAlign(11);
      latexR -> SetTextColor(kBlue);
      latexR -> Draw("same");
      
      c -> Print(Form("%s/g_N_bar%02d_Vov%.2f.png",plotDir.c_str(),iBar,Vov));
      c -> Print(Form("%s/g_N_bar%02d_Vov%.2f.pdf",plotDir.c_str(),iBar,Vov));
      
      delete c;
    }


  
  for(auto mapIt : g_totL)
  {  
    int index2 = mapIt.first;
    float Vov =  int(index2/10000)/100.;
    int iBar = index2 - (Vov*10000*100);

    c = new TCanvas("c","c");
    //    hPad = (TH1F*)( gPad->DrawFrame(-0.5,0.,63.5,50.) );
    hPad = (TH1F*)( gPad->DrawFrame(-0.5,0.,63.5,2.) );
    hPad -> SetTitle(Form(";%s [DAC]; ToT [ns]",ithMode.c_str()));
    hPad -> Draw();
    g_totL[index2] -> SetMarkerColor(kRed);
    g_totL[index2] -> SetMarkerStyle(20);
    g_totL[index2] -> Draw("PL,same");
    if( g_totR[index2] ) g_totR[index2] -> SetMarkerColor(kBlue);
    if( g_totR[index2] ) g_totR[index2] -> SetMarkerStyle(20);
    if( g_totR[index2] ) g_totR[index2] -> Draw("PL,same");
    c -> Print(Form("%s/g_tot_bar%02d_Vov%.2f.png",plotDir.c_str(),iBar,Vov));
    c -> Print(Form("%s/g_tot_bar%02d_Vov%.2f.pdf",plotDir.c_str(),iBar,Vov));
    delete c;

    g_totL[index2] -> Write(Form("g_totL_bar%02d_Vov%.2f",iBar,Vov));
    g_totR[index2] -> Write(Form("g_totR_bar%02d_Vov%.2f",iBar,Vov));
  }
  
  for(auto mapIt : g_energyLR)
  {  
    int index2 = mapIt.first;
    float Vov =  int(index2/10000)/100.;
    int iBar = index2 - (Vov*10000*100);

    c = new TCanvas("c","c");
    hPad = (TH1F*)( gPad->DrawFrame(-0.5,0.,63.5,g_energyLR[index2] -> GetMean(2)*1.5) );
    hPad -> SetTitle(Form(";%s [DAC]; energy [ADC]",ithMode.c_str()));
    hPad -> Draw();
    g_energyLR[index2] -> SetMarkerColor(kRed);
    g_energyLR[index2] -> SetMarkerStyle(20);
    g_energyLR[index2] -> Draw("PL,same");
    c -> Print(Form("%s/g_energy_bar%02d_Vov%.2f.png",plotDir.c_str(),iBar,Vov));
    c -> Print(Form("%s/g_energy_bar%02d_Vov%.2f.pdf",plotDir.c_str(),iBar,Vov));
    delete c;
  }
    
  for(auto mapIt : g_pulseShapeL)
  {  
    int index2 = mapIt.first;
    float Vov =  int(index2/10000)/100.;
    int iBar = index2 - (Vov*10000*100);

    c = new TCanvas("c","c");
    hPad = (TH1F*)( gPad->DrawFrame(-2.,0.,25.,65*dac_to_uA) );
    hPad -> SetTitle(Form(";time [ns]; pulse shape [#muA]"));
    hPad -> Draw();
    g_pulseShapeL[index2] -> SetLineColor(kRed-4);
    g_pulseShapeL[index2] -> SetMarkerColor(kRed-4);
    g_pulseShapeL[index2] -> SetMarkerStyle(20);
    g_pulseShapeL[index2] -> Draw("P,same");
    if( g_pulseShapeR[index2] ) g_pulseShapeR[index2] -> SetMarkerColor(kBlue-4);
    if( g_pulseShapeR[index2] ) g_pulseShapeR[index2] -> SetLineColor(kBlue-4);
    if( g_pulseShapeR[index2] ) g_pulseShapeR[index2] -> SetMarkerStyle(20);
    if( g_pulseShapeR[index2] ) g_pulseShapeR[index2] -> Draw("P,same");

    float slewRate = 0.;
    TF1* fitFuncL = new TF1("fitFuncL","pol1",-5.,7.);
    //-- slew rate max 
    int npoints = 4;
    for(int point1 = 0; point1 <  g_pulseShapeL[index2]->GetN()-npoints; ++point1)
    {
      TGraph* g_temp = new TGraph();
      for(int point2 = point1; point2 < point1+npoints; ++point2)
      {
        g_temp -> SetPoint(g_temp->GetN(), g_pulseShapeL[index2]->GetPointX(point2), g_pulseShapeL[index2]->GetPointY(point2));
      }
      
      TF1* f_temp = new TF1("f_temp","pol1",-10.,100.);
      g_temp -> Fit(f_temp,"QNRS");
      
      if( f_temp->GetParameter(1) > slewRate )
      {
        slewRate = f_temp->GetParameter(1);
        fitFuncL -> SetParameters(f_temp->GetParameter(0),f_temp->GetParameter(1));
	fitFuncL -> SetRange(g_temp->GetPointX(0), g_temp->GetPointX(g_temp->GetN()-1));
      }
      delete g_temp;
    }

    //-- slew rate at low threshold
    TF1* fitFuncLowL = new TF1("fitFuncLowL","pol1",-5.,7.);
    TGraph* g_temp = new TGraph();
    for(int point1 = 0; point1 < npoints; ++point1){
      g_temp -> SetPoint(g_temp->GetN(), g_pulseShapeL[index2]->GetPointX(point1), g_pulseShapeL[index2]->GetPointY(point1)); 
    }
    TF1* f_temp = new TF1("f_temp","pol1", g_temp ->GetPointX(0), g_temp ->GetPointX(npoints-1));
    f_temp->SetParameter(1,10.);
    g_temp -> Fit(f_temp,"QNRS"); 
    fitFuncLowL -> SetParameters(f_temp->GetParameter(0),f_temp->GetParameter(1));
    fitFuncLowL -> SetRange(g_temp ->GetPointX(0), g_temp ->GetPointX(npoints-1));
    float slewRate_low = fitFuncLowL->GetParameter(1); 
    std::cout << "ch1 - Slew rate max              = " << slewRate << std::endl;
    std::cout << "ch1 - Slew rate at low threshold = " << slewRate_low << std::endl;
    delete g_temp;

    // -- draw
    fitFuncL -> SetLineColor(kRed-4);
    fitFuncL -> Draw("same");
    TLatex* latexL = new TLatex(0.30,0.80,Form("slew rate max = %.1f #muA/ns",fitFuncL->GetParameter(1)));
    latexL -> SetNDC();
    latexL -> SetTextFont(82);
    latexL -> SetTextSize(0.04);
    latexL -> SetTextAlign(11);
    latexL -> SetTextColor(kRed-4);
    latexL -> Draw("same");

    fitFuncLowL -> SetLineColor(kRed+2);
    fitFuncLowL -> Draw("same");
    TLatex* latexL_low = new TLatex(0.30,0.70,Form("slew rate timing th. = %.1f #muA/ns",fitFuncLowL->GetParameter(1)));
    latexL_low -> SetNDC();
    latexL_low -> SetTextFont(82);
    latexL_low -> SetTextSize(0.04);
    latexL_low -> SetTextAlign(11);
    latexL_low -> SetTextColor(kRed+2);
    latexL_low -> Draw("same");
    
    g_pulseShapeL[index2] ->Write(Form("g_pulseShapeL_bar%02d_Vov%.2f",iBar,Vov));

    // --- channel R
    if( !g_pulseShapeR[index2] ) continue;
    slewRate = 0.;
    TF1* fitFuncR = new TF1("fitFuncR","pol1",-5.,7.);
    for(int point1 = 0; point1 <  g_pulseShapeR[index2]->GetN()-npoints; ++point1)
    {
      TGraph* g_temp = new TGraph();
      for(int point2 = point1; point2 < point1+npoints; ++point2)
      {
        g_temp -> SetPoint(g_temp->GetN(), g_pulseShapeR[index2]->GetPointX(point2), g_pulseShapeR[index2]->GetPointY(point2));
      }
      
      TF1* f_temp = new TF1("f_temp","pol1",-10.,100.);
      g_temp -> Fit(f_temp,"QNRS");
      
      if( f_temp->GetParameter(1) > slewRate )
      {
        slewRate = f_temp->GetParameter(1);
        fitFuncR -> SetParameters(f_temp->GetParameter(0),f_temp->GetParameter(1));
	fitFuncR -> SetRange(g_temp->GetPointX(0), g_temp->GetPointX(g_temp->GetN()-1));
      }
      delete g_temp;
    }
    
    TF1* fitFuncLow_ch2 = new TF1("fitFuncLow_ch2","pol1",-5.,7.);
    //-- slew rate at low threshold
    g_temp = new TGraph();
    for(int point1 = 0; point1 < npoints; ++point1){
      g_temp -> SetPoint(g_temp->GetN(), g_pulseShapeR[index2]->GetPointX(point1), g_pulseShapeR[index2]->GetPointY(point1)); 
    }
    f_temp = new TF1("f_temp","pol1", g_temp ->GetPointX(0), g_temp ->GetPointX(npoints-1));
    f_temp->SetParameter(1,10.);
    g_temp -> Fit(f_temp,"QNRS"); 
    fitFuncLow_ch2 -> SetParameters(f_temp->GetParameter(0),f_temp->GetParameter(1));
    fitFuncLow_ch2 -> SetRange(g_temp ->GetPointX(0), g_temp ->GetPointX(npoints-1));
    slewRate_low = fitFuncLow_ch2->GetParameter(1); 
    std::cout << "ch2 - Slew rate max              = " << slewRate << std::endl;
    std::cout << "ch2 - Slew rate at low threshold = " << slewRate_low << std::endl;
    delete g_temp;



    fitFuncR -> SetLineColor(kBlue-4);
    fitFuncR -> Draw("same");
    TLatex* latexR = new TLatex(0.30,0.76,Form("slew rate max = %.1f #muA/ns",fitFuncR->GetParameter(1)));
    latexR -> SetNDC();
    latexR -> SetTextFont(82);
    latexR -> SetTextSize(0.04);
    latexR -> SetTextAlign(11);
    latexR -> SetTextColor(kBlue-4);
    latexR -> Draw("same");


    fitFuncLow_ch2 -> SetLineColor(kBlue+2);
    fitFuncLow_ch2 -> Draw("same");
    TLatex* latexR_low = new TLatex(0.30,0.66,Form("slew rate timing th. = %.1f #muA/ns",fitFuncLow_ch2->GetParameter(1)));
    latexR_low -> SetNDC();
    latexR_low -> SetTextFont(82);
    latexR_low -> SetTextSize(0.04);
    latexR_low -> SetTextAlign(11);
    latexR_low -> SetTextColor(kBlue+2);
    latexR_low -> Draw("same");
    
    c -> Print(Form("%s/g_ps_bar%02d_Vov%.2f.png",plotDir.c_str(),iBar,Vov));
    c -> Print(Form("%s/g_ps_bar%02d_Vov%.2f.pdf",plotDir.c_str(),iBar,Vov));

    delete c;

    g_pulseShapeR[index2] ->Write(Form("g_pulseShapeR_bar%02d_Vov%.2f",iBar,Vov));
  }
  
    
  outFile -> Close();
}
