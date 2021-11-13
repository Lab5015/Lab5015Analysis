#include "interface/TOFHIRThresholdZero.h"
#include "interface/AnalysisUtils.h"
#include "interface/FitUtils.h"
#include "interface/SetTDRStyle.h"
#include "interface/Na22SpectrumAnalyzer.h"
#include "interface/Na22SpectrumAnalyzerSingleBar.h"
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
  
  
  //--- open files and make the tree chain
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  std::string fileBaseName = opts.GetOpt<std::string>("Input.fileBaseName");
  std::string runs = opts.GetOpt<std::string>("Input.runs");
  int maxEntries = opts.GetOpt<int>("Input.maxEntries");
  int usePedestals = opts.GetOpt<int>("Input.usePedestals");
  std::string source = opts.GetOpt<std::string>("Input.sourceName");
  int useTrackInfo = opts.GetOpt<int>("Input.useTrackInfo");
  
  std::string discCalibrationFile = opts.GetOpt<std::string>("Input.discCalibration");
  TOFHIRThresholdZero thrZero(discCalibrationFile,1);

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
    
    for(int run = runMin; run <= runMax; ++run) {
      std::string fileName;
      if( !usePedestals ) fileName = Form("%s/%s%04d_*e*.root",inputDir.c_str(),fileBaseName.c_str(),run);
      //else                fileName = Form("%s/%s%04d_*ped_e.root",inputDir.c_str(),fileBaseName.c_str(),run);
      else                fileName = Form("%s/%04d/*ped_e.root",inputDir.c_str(),run);
      std::cout << ">>> Adding file " << fileName << std::endl;
      tree -> Add(fileName.c_str());
      
      struct stat t_stat;
      //stat(Form("/data/TOFHIR2/raw/run%04d.rawf",run), &t_stat);
      stat(Form("/data/tofhir2/h8/raw/%04d/",run), &t_stat);
      struct tm * timeinfo = localtime(&t_stat.st_mtime);
      std::cout << "Time and date of raw file of run" << run << ": " << asctime(timeinfo);
    }
  }
  
     
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
    std::cout<<chL[iBar]<<"    "<<chR[iBar]<<std::endl;
  }
  
  //--- define branches
  float step1, step2;
  int channelIdx[256];
  std::vector<unsigned short> *qfine = 0;
  std::vector<float> *tot = 0;
  std::vector<float> *energy = 0;
  std::vector<long long> *time = 0;
  std::vector<float>* qT1 = 0;
  std::vector<unsigned short>* t1fine = 0;
  
  int nhits;
  float x, y;
  
  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("step1",  1); tree -> SetBranchAddress("step1",  &step1);
  tree -> SetBranchStatus("step2",  1); tree -> SetBranchAddress("step2",  &step2);
  
  tree -> SetBranchStatus("channelIdx",  1); tree -> SetBranchAddress("channelIdx",  channelIdx);
  
  tree -> SetBranchStatus("qfine",  1); tree -> SetBranchAddress("qfine",   &qfine);  
  tree -> SetBranchStatus("tot",    1); tree -> SetBranchAddress("tot",       &tot);
  tree -> SetBranchStatus("energy", 1); tree -> SetBranchAddress("energy", &energy);
  tree -> SetBranchStatus("time",   1); tree -> SetBranchAddress("time",     &time);
  
  tree -> SetBranchStatus("qT1",       1); tree -> SetBranchAddress("qT1",      &qT1);
  tree -> SetBranchStatus("t1fine",    1); tree -> SetBranchAddress("t1fine",   &t1fine);
  
  if ( !opts.GetOpt<std::string>("Input.sourceName").compare("TB") &&  useTrackInfo ){
    tree -> SetBranchStatus("nhits_WC", 1); tree -> SetBranchAddress("nhits_WC",  &nhits);
    tree -> SetBranchStatus("x_WC", 1);     tree -> SetBranchAddress("x_WC",          &x);
    tree -> SetBranchStatus("y_WC", 1);     tree -> SetBranchAddress("y_WC",          &y);
  }
  

  //--- get plot settings
  std::vector<float> Vov = opts.GetOpt<std::vector<float> >("Plots.Vov");
  std::vector<int> energyBins = opts.GetOpt<std::vector<int> >("Plots.energyBins");
  std::vector<int> energyMins = opts.GetOpt<std::vector<int> >("Plots.energyMins");
  std::vector<int> energyMaxs = opts.GetOpt<std::vector<int> >("Plots.energyMaxs");
  
  std::map<float,int> map_energyBins;
  std::map<float,int> map_energyMins;
  std::map<float,int> map_energyMaxs;
  for(unsigned int ii = 0; ii < Vov.size(); ++ii)
    {
      map_energyBins[Vov[ii]] = energyBins[ii];
      map_energyMins[Vov[ii]] = energyMins[ii];
      map_energyMaxs[Vov[ii]] = energyMaxs[ii];
    }
  
  
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
  std::map<int,TH1F*> h1_energyLR_ext;
  std::map<int,TCanvas*> c;
  std::map<int,std::vector<float>*> rangesLR;
  std::map<int,std::map<std::string,std::pair<float,float> > > peaksLR;
  std::map<int,bool> acceptEvent;

  // -- Coincidence pre loop
  if( !opts.GetOpt<std::string>("Coincidence.status").compare("yes") &&
      ( !opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar") ||
	!opts.GetOpt<std::string>("Input.sourceName").compare("Na22") ||
	!opts.GetOpt<std::string>("Input.sourceName").compare("TB") ||
	!opts.GetOpt<std::string>("Input.sourceName").compare("keepAll") ) )
    {
      float energyL_ext;
      float energyR_ext;
      int chL_ext = opts.GetOpt<float>("Coincidence.chL");//NB: gli passo direttamente da cfg il ch barra 8 +64
      int chR_ext = opts.GetOpt<float>("Coincidence.chR");//NB: gli passo direttamente da cfg il ch barra 8 +64
      
      int nEntries = tree->GetEntries();
      if( maxEntries > 0 ) nEntries = maxEntries;
      for(int entry = 0; entry < nEntries; ++entry){
	tree -> GetEntry(entry);
	acceptEvent[entry] = false;
	if( entry%200000 == 0 ){
	  std::cout << "\n>>> external bar loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << std::endl;
	  TrackProcess(cpu, mem, vsz, rss);
	}
      
	float Vov = step1;
        float vth1 = float(int(step2/10000)-1);
        float vth2 = int((step2-10000*(vth1+1))/100.)-1;
        float vth = 0;
	if(!opts.GetOpt<std::string>("Input.vth").compare("vth1"))  { vth = vth1;}
	if(!opts.GetOpt<std::string>("Input.vth").compare("vth2"))  { vth = vth2;}
       	
	// --- calculate energy sum for module - useful to remove showering events
	if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB")){

	    float energySumArray = 0.;
	    int   nBarsArray = 0;
	    for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar) {                             
	      int chL_iext = channelMapping[iBar*2+0];// array0 for coincidence is hard coded... - to be fixed
	      int chR_iext = channelMapping[iBar*2+1];// array0 for coincidence is hard coded... - to be fixed
	      float energyL_iext = (*energy)[channelIdx[chL_iext]];              
	      float energyR_iext = (*energy)[channelIdx[chR_iext]]; 
	      float totL_iext    = 0.001*(*tot)[channelIdx[chL_iext]];              
	      float totR_iext    = 0.001*(*tot)[channelIdx[chR_iext]]; 
	      if ( totL_iext > 0 && totL_iext < 20 && totR_iext > 0 && totR_iext < 20   ){
		float energyMean=(energyL_iext+energyR_iext)/2;
		if (energyMean>0){
		  energySumArray+=energyMean;
		  nBarsArray+=1;
		}
	      }
	    }
	    if (energySumArray > 800. ) continue;
	    if (nBarsArray > 5 ) continue;
	}

	if (channelIdx[chL_ext] <0 || channelIdx[chR_ext] <0) continue;
	
	energyL_ext = (*energy)[channelIdx[chL_ext]];
	energyR_ext = (*energy)[channelIdx[chR_ext]];
	
	int index( (10000*int(Vov*100.)) + (100*vth) + 99 );
	
	//--- create histograms, if needed
	if( h1_energyLR_ext[index] == NULL ){
	  c[index] = new TCanvas(Form("c1_Vov%.2f_th%02.0f",Vov,vth), Form("c1_Vov%.2f_th%02.0f",Vov,vth));
	  c[index] -> cd();
	  h1_energyLR_ext[index] = new TH1F(Form("h1_energy_external_barL-R_Vov%.2f_th%02.0f",Vov,vth),"",map_energyBins[Vov],map_energyMins[Vov],map_energyMaxs[Vov]);
	}

	acceptEvent[entry] = true;
	
	h1_energyLR_ext[index] -> Fill(0.5*(energyL_ext + energyR_ext));    
      }
      
      std::cout << std::endl;
      
      if (!opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar") || !opts.GetOpt<std::string>("Input.sourceName").compare("Na22")){		
	for( auto index : h1_energyLR_ext){
	  rangesLR[index.first] = new std::vector<float>;
	  if(!opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar")) peaksLR[index.first] = Na22SpectrumAnalyzerSingleBar(index.second,rangesLR[index.first]);
	  if(!opts.GetOpt<std::string>("Input.sourceName").compare("Na22")) peaksLR[index.first] = Na22SpectrumAnalyzer(index.second,rangesLR[index.first]);
	}
      }
      
      if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB")){
	for( auto index : h1_energyLR_ext){
	  rangesLR[index.first] = new std::vector<float>;
	  
	  float Vov = float ((int(index.first /10000))/100.);
	  float vth1 = float(int((index.first-Vov*100000*100)/100.));
	  float vth2 = int((step2-10000*(vth1+1))/100.)-1;
          float vth = 0;
          if(!opts.GetOpt<std::string>("Input.vth").compare("vth1"))  { vth = vth1;}
          if(!opts.GetOpt<std::string>("Input.vth").compare("vth2"))  { vth = vth2;}
	  
	  if( opts.GetOpt<int>("Channels.array") == 0){
	    index.second->GetXaxis()->SetRangeUser(Vov*20.,700);
	  }
	  if( opts.GetOpt<int>("Channels.array") == 1){
	    TF1 *ftemp = new TF1("ftemp","gaus",0,1000);
	    index.second->Fit(ftemp,"QRS");
	    //index.second->GetXaxis()->SetRangeUser(Vov*20.,700);
	    index.second->GetXaxis()->SetRangeUser(ftemp->GetParameter(1),700);
	    
	    // July21 TB: runs 4184 and 4212 are at fixed OV and thr!
	    if ( (opts.GetOpt<std::string>("Input.runs")).find("4184") != std::string::npos ||
		 (opts.GetOpt<std::string>("Input.runs")).find("4212") != std::string::npos )
	      { index.second->GetXaxis()->SetRangeUser(200,700);}
	  }

	  float max = index.second->GetBinCenter(index.second->GetMaximumBin());
	  index.second->GetXaxis()->SetRangeUser(0,1000);
	  
	  TF1* f_pre = new TF1(Form("fit_energy_coincBar_Vov%.2f_vth1_%02.0f",Vov,vth), "[0]*TMath::Landau(x,[1],[2])", 0, 1000.); 
	  f_pre -> SetRange(max*0.8, max*1.5);
	  f_pre -> SetLineColor(kBlack);
          f_pre -> SetLineWidth(2);
          f_pre -> SetParameters(index.second->Integral(index.second->GetMaximumBin(), index.second->GetNbinsX())/10, max, 5*Vov);
	  if ( (opts.GetOpt<std::string>("Input.runs")).find("4184") != std::string::npos ||  (opts.GetOpt<std::string>("Input.runs")).find("4212") != std::string::npos ) 
	    f_pre -> SetParameters(index.second->Integral(index.second->GetMaximumBin(), index.second->GetNbinsX())/10, max, 25);
	  index.second->Fit(f_pre, "QRS");      
	  
	  if (f_pre->GetParameter(1)>10)
	    rangesLR[index.first] -> push_back( 0.80*f_pre->GetParameter(1));
	  else
	    rangesLR[index.first] -> push_back( 20 );
	  rangesLR[index.first] -> push_back( 700. );
	  
	  std::cout << "Vov = " << Vov << "  vth1 = " << vth1 << "   vth2 = " << vth2 
		    << "    Coincidence bar - energy range:  " << rangesLR[index.first]->at(0) << " - " << rangesLR[index.first]->at(1)<< std::endl;
	}
      }
    }
  
  
  
  //------------------------
  //--- 1st loop over events
  
  ModuleEventClass anEvent;

  unsigned short qfineL[16];
  unsigned short qfineR[16];    
  float totL[16];
  float totR[16];
  long long timeL[16];
  long long timeR[16];
  unsigned short t1fineL[16]; 
  unsigned short t1fineR[16]; 
  float energyL[16];
  float energyR[16];
  
  int nEntries = tree->GetEntries();
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry) {
    tree -> GetEntry(entry);
    if( entry%200000 == 0 )
      {
	std::cout << "\n>>> 1st loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << std::endl;
	TrackProcess(cpu, mem, vsz, rss);
      }
    
    if (useTrackInfo && nhits > 0 &&  (x < -100 || y < -100 ) ) continue;

    float Vov = step1;
    float vth1 = float(int(step2/10000)-1);
    float vth2 = int((step2-10000*(vth1+1))/100.)-1;
    float vth = 0;
    std::string vthMode = opts.GetOpt<std::string>("Input.vth");
    if(!opts.GetOpt<std::string>("Input.vth").compare("vth1"))  { vth = vth1;}
    if(!opts.GetOpt<std::string>("Input.vth").compare("vth2"))  { vth = vth2;}
    // float vthe = float(int((step2-10000*vth1-step2-100*vth2)/1)-1);
    
    
    // --- check coincidence with another channel 
    if(!opts.GetOpt<std::string>("Coincidence.status").compare("yes"))
      {
	if(!acceptEvent[entry] ) continue;
	
	int chL_ext = opts.GetOpt<int>("Coincidence.chL");
	int chR_ext = opts.GetOpt<int>("Coincidence.chR");
	float energyL_ext = (*energy)[channelIdx[chL_ext]];
	float energyR_ext = (*energy)[channelIdx[chR_ext]];
	
	int label = (10000*int(Vov*100.)) + (100*vth) + 99;
	int eBin = opts.GetOpt<int>("Coincidence.peak511eBin");
	float avEn = 0.5 * ( energyL_ext + energyR_ext);
	if ( (!opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar") || !opts.GetOpt<std::string>("Input.sourceName").compare("Na22")) &&  avEn > rangesLR[label]-> at(eBin)) {
	  continue;
	}
	
	if ( (!opts.GetOpt<std::string>("Input.sourceName").compare("TB")) && ( avEn < rangesLR[label]-> at(0) || avEn > rangesLR[label]-> at(1) ) ) {
	  continue;
	}
      }
    
    
    for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar)
      {
	if (channelIdx[chL[iBar]] >=0 && channelIdx[chR[iBar]] >=0)
	  {
	    qfineL[iBar]=(*qfine)[channelIdx[chL[iBar]]];
	    qfineR[iBar]=(*qfine)[channelIdx[chR[iBar]]];
	    totL[iBar]=0.001*(*tot)[channelIdx[chL[iBar]]];
	    totR[iBar]=0.001*(*tot)[channelIdx[chR[iBar]]];
	    energyL[iBar]=(*energy)[channelIdx[chL[iBar]]];
	    energyR[iBar]=(*energy)[channelIdx[chR[iBar]]];
	    timeL[iBar]=(*time)[channelIdx[chL[iBar]]];
	    timeR[iBar]=(*time)[channelIdx[chR[iBar]]];
	    t1fineL[iBar]=(*t1fine)[channelIdx[chL[iBar]]];
	    t1fineR[iBar]=(*t1fine)[channelIdx[chR[iBar]]];
	}
	else
	  {
	    qfineL[iBar]=-10;
	    qfineR[iBar]=-10;
	    totL[iBar]=-10;
	    totR[iBar]=-10;
	    energyL[iBar]=-10;
	    energyR[iBar]=-10;
	    timeL[iBar]=-10;
	    timeR[iBar]=-10;
	    t1fineL[iBar]=-10;
	    t1fineR[iBar]=-10;
	  }     
      }
    
    int maxEn=0;
    int maxBar=0;
    
    float energySumArray = 0;
    int   nBarsArray = 0;
    
    for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar)
      {
	if (totL[iBar]>0 && totR[iBar]>0 && totL[iBar]<100 && totR[iBar]<100){
	  float energyMean=(energyL[iBar]+energyR[iBar])/2;
	  energySumArray+=energyMean;
	  nBarsArray+=1;
	  if(energyMean>maxEn){
	    maxEn = energyMean;
	    maxBar = iBar;
	  }
	}
    
	int index( (10000*int(Vov*100.)) + (100*vth) + iBar );
	
	
	//--- create histograms, if needed
	if( h1_totL[index] == NULL )
	{
	  h1_qfineL[index] = new TH1F(Form("h1_qfine_bar%02dL_Vov%.2f_th%02.0f",iBar,Vov,vth),"",512,-0.5,511.5);
	  h1_qfineR[index] = new TH1F(Form("h1_qfine_bar%02dR_Vov%.2f_th%02.0f",iBar,Vov,vth),"",512,-0.5,511.5);
	  
	  h1_totL[index] = new TH1F(Form("h1_tot_bar%02dL_Vov%.2f_th%02.0f",iBar,Vov,vth),"",500,0.,100.);
	  h1_totR[index] = new TH1F(Form("h1_tot_bar%02dR_Vov%.2f_th%02.0f",iBar,Vov,vth),"",500,0.,100.);
	  
	  h1_energyL[index] = new TH1F(Form("h1_energy_bar%02dL_Vov%.2f_th%02.0f",iBar,Vov,vth),"",map_energyBins[Vov],map_energyMins[Vov],map_energyMaxs[Vov]);
	  h1_energyR[index] = new TH1F(Form("h1_energy_bar%02dR_Vov%.2f_th%02.0f",iBar,Vov,vth),"",map_energyBins[Vov],map_energyMins[Vov],map_energyMaxs[Vov]);
	  
	  outTrees[index] = new TTree(Form("data_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth),Form("data_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth));
	  outTrees[index] -> Branch("event",&anEvent);
	  
	  h1_energyLR[index] = new TH1F(Form("h1_energy_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth),"",map_energyBins[Vov],map_energyMins[Vov],map_energyMaxs[Vov]);
	}
      
	
	//--- fill histograms for each bar for Co60 & TB analysis
	if( !opts.GetOpt<std::string>("Input.sourceName").compare("Co60") ||
	    !opts.GetOpt<std::string>("Input.sourceName").compare("Co60SumPeak") ||
	    !opts.GetOpt<std::string>("Input.sourceName").compare("TB") ||
	    !opts.GetOpt<std::string>("Input.sourceName").compare("keepAll")
	    )
	{
	  
	  if( totL[iBar] <= 0. || totR[iBar] <= 0. ) continue;
	  if( totL[iBar] >= 50. ||  totR[iBar] >= 50.) continue;
	  if( ( thrZero.GetThresholdZero(chL[iBar],vthMode) + vth) > 63. ) continue;
          if( ( thrZero.GetThresholdZero(chR[iBar],vthMode) + vth) > 63. ) continue;

	  if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB") && (energySumArray > 800 || nBarsArray > 5)) continue; // to remove showering events


	  h1_qfineL[index] -> Fill( qfineL[iBar] );
	  h1_totL[index] -> Fill( totL[iBar]  );
	  h1_energyL[index] -> Fill( energyL[iBar] );
	  
	  h1_qfineR[index] -> Fill( qfineR[iBar] );
	  h1_totR[index] -> Fill( totR[iBar] );
	  h1_energyR[index] -> Fill( energyR[iBar] );
	  
	  h1_energyLR[index] -> Fill(0.5*(energyL[iBar]+energyR[iBar]));
	  
	  anEvent.barID = iBar;
	  anEvent.Vov = Vov;
	  anEvent.vth1 = vth;
	  anEvent.energyL = energyL[iBar];
	  anEvent.energyR = energyR[iBar];
	  anEvent.totL = totL[iBar];
	  anEvent.totR = totR[iBar];
	  anEvent.timeL = timeL[iBar];
	  anEvent.timeR = timeR[iBar];
	  anEvent.t1fineL = t1fineL[iBar];
	  anEvent.t1fineR = t1fineR[iBar];
	  if(useTrackInfo){
	    anEvent.nhits = nhits;
	    anEvent.x = x;
	    anEvent.y = y;
	  }
	  else{
	    anEvent.nhits = -1;
	    anEvent.x = -999.;
	    anEvent.y = -999.;
	  }
	  
	  
	  outTrees[index] -> Fill();
	  
	}
      }// -- end loop over bars
    
    // --- for Na22 or Laser analysis use only the bar with max energy to remove cross-talk between adjacent bars
    if( !opts.GetOpt<std::string>("Input.sourceName").compare("Na22") ||
	!opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar") ||
	!opts.GetOpt<std::string>("Input.sourceName").compare("Laser"))
      {
	int index( (10000*int(Vov*100.)) + (100*vth) + maxBar );
	
	if( totL[maxBar] <= 0. || totR[maxBar] <= 0. ) continue;
	if( totL[maxBar] >= 50. ||  totR[maxBar] >= 50.) continue;
	if( ( thrZero.GetThresholdZero(chL[maxBar],vthMode) + vth) > 63. ) continue;
        if( ( thrZero.GetThresholdZero(chR[maxBar],vthMode) + vth) > 63. ) continue;

	//--- fill histograms
	h1_qfineL[index] -> Fill( qfineL[maxBar] );
	h1_totL[index] -> Fill( totL[maxBar] );
	h1_energyL[index] -> Fill( energyL[maxBar] );
	
	h1_qfineR[index] -> Fill( qfineR[maxBar] );
	h1_totR[index] -> Fill( totR[maxBar] );
	h1_energyR[index] -> Fill( energyR[maxBar] );
	  
	h1_energyLR[index] -> Fill(0.5*(energyL[maxBar]+energyR[maxBar]));
	
	anEvent.barID = maxBar;
	anEvent.Vov = Vov;
	anEvent.vth1 = vth;
	anEvent.energyL = energyL[maxBar];
	anEvent.energyR = energyR[maxBar];
	anEvent.totL = totL[maxBar];
	anEvent.totR = totR[maxBar];
	anEvent.timeL = timeL[maxBar];
	anEvent.timeR = timeR[maxBar];
	anEvent.t1fineL = t1fineL[maxBar];
	anEvent.t1fineR = t1fineR[maxBar];
	if(useTrackInfo){
	  anEvent.nhits = nhits;
	  anEvent.x = x;
	  anEvent.y = y;
	}
	else{
	  anEvent.nhits = -1;
	  anEvent.x = -999.;
	  anEvent.y = -999.;
	}
	outTrees[index] -> Fill();
      }
  } // --- end loop over events
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  
}
