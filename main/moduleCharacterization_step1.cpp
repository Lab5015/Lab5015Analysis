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
#include <sstream>
#include <vector>
#include <map>
#include <time.h>
#include <stdio.h>
#include <sys/stat.h>
#include <algorithm>
#include <iterator>
#include <dirent.h>

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
  
  // - parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  // - open files and make the tree chain
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  std::string fileBaseName = opts.GetOpt<std::string>("Input.fileBaseName");
  std::string runs = opts.GetOpt<std::string>("Input.runs");
  int maxEntries = opts.GetOpt<int>("Input.maxEntries");
  std::string source = opts.GetOpt<std::string>("Input.sourceName");
  int useTrackInfo = opts.GetOpt<int>("Input.useTrackInfo");
  int saveRefInfo = opts.GetOpt<int>("Input.saveRefInfo");
  float my_step1 = opts.GetOpt<float>("Input.vov") ;
  int DUTasic = opts.GetOpt<int>("Channels.DUTasic");
  int REFasic = opts.GetOpt<int>("Channels.REFasic");
  std::string discCalibrationFile = opts.GetOpt<std::string>("Input.discCalibration");
  std::string interCalibrationFile = opts.GetOpt<std::string>("Input.interCalibration");
  std::string interCalibrationFile_ref = opts.GetOpt<std::string>("Input.interCalibration_ref");
  std::vector<unsigned int> channelMapping = opts.GetOpt<std::vector<unsigned int> >("Channels.channelMapping");
  // - Note: the channels corresponding to the external reference bar are directly set in the cfg. The coincidence is required only if it is specified in the cfg
  int chL_ext = opts.GetOpt<float>("Coincidence.chL");
  int chR_ext = opts.GetOpt<float>("Coincidence.chR");
  // - cross check that channelMapping is behaving as expected and get reference bar
  std::string sideL_ext;
  std::string sideR_ext;
  // - Note: the current logic is valid for reference module coincidence requirement, it may vary for other conditions (e.g. lab ref on one channel)
  int barL_ext = findBarAndSide(chL_ext, REFasic, channelMapping, sideL_ext);
  int barR_ext = findBarAndSide(chR_ext, REFasic, channelMapping, sideR_ext);
  std::cout<<"[INFO] Coincidence required on channels corresponding to bar "<<barL_ext <<std::endl;
  if (sideL_ext != "L" || sideR_ext != "R" || barL_ext != barR_ext){
    std::cout <<"[ERROR] channelMapping has something wrong "<<std::endl;
    return 1;
  }
  // - [WARNING]: hard-coded, but always valid for Sep25 TB
  float Vov_ext = 3.0;
  float th_ext = 10.0;
  
  std::vector<std::string> zombieFiles;
  TOFHIRThresholdZero thrZero(discCalibrationFile,0);
  int maxActiveBars = 3;
  TChain* tree = new TChain("data","data");

  // - energy intercalibration map (if specified in the config) for the DUT module
  std::map<std::tuple<int,std::string,int,int>, float> calibMap;
  if (interCalibrationFile != "0") {
    calibMap = loadCalibrationMap(interCalibrationFile);
  }
  // - energy intercalibration map for the REF module
  std::map<std::tuple<int,std::string,int,int>, float> calibMap_ref;
  if (interCalibrationFile_ref != "0" && saveRefInfo == 1) {
    calibMap_ref = loadCalibrationMap(interCalibrationFile_ref);
  }
  
  // - determine run numbers
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
      // --- analyze only spills at a chosen OV to speed up analysis
      // --- list of files in run folder
      DIR *dir_ptr;
      struct dirent *diread;
      std::vector<std::string> filenames;
      std::string directory_path = Form("%s/%04d/",inputDir.c_str(),run);
      std::cout << directory_path.c_str()<<std::endl;
      if ((dir_ptr = opendir(directory_path.c_str())) != nullptr) {
        while ((diread = readdir(dir_ptr)) != nullptr) {
	  std::string fname(diread->d_name);
	  filenames.push_back(fname);
	}
	closedir(dir_ptr);
      }

      // --- loop on filenames to add events to tree
      for (auto fname: filenames) {
	if (fname == ".") continue;
	if (fname == "..") continue;

	// ---- check if Vov selected
	bool addFile = true;
	TFile *f = TFile::Open((directory_path+fname).c_str());
	std::string fullPath = directory_path + fname;
	if (!f || f->IsZombie()) {
	  std::cerr << "[ERROR] Cannot open file: " << fullPath << std::endl;
	  zombieFiles.push_back(fullPath);
	  if(f) f->Close();
	  continue;
	}
	TTree *tmpTree = f->Get<TTree>("data");
	if (!tmpTree) {
	  std::cerr << "[ERROR] TTree 'data' not found in file: " << fullPath << std::endl;
	  zombieFiles.push_back(fullPath);
	  f->Close();
	  continue;
	}
	float step1;
	tmpTree->SetBranchAddress("step1",&step1);
	tmpTree->GetEntry(0);
	if (my_step1 > 0 && step1!=my_step1) addFile = false;
	delete tmpTree;
	f->Close();	
	if (addFile){
	  std::cout << ">>> step1 = " << step1 << " --> Adding file: " << fname.c_str()<< "\r" << std::flush;
	  tree->Add((directory_path+fname).c_str());
	}
      } // end of loop on filenames

      // --- print stats for raw files
      std::string rawDir = inputDir;
      if (auto pos = rawDir.find("reco"); pos != std::string::npos)
	rawDir.replace(pos, 4, "raw");
      struct stat t_stat;
      std::string fullPath = Form("%s/%04d/", rawDir.c_str(), run);
      stat(fullPath.c_str(), &t_stat);
      struct tm * timeinfo = localtime(&t_stat.st_mtime);
      std::cout << "Time and date of raw file of run" << run << ": " << asctime(timeinfo);
    } // end of loop on runs
  } // end of loop on token lines

  // - define channels (read mapping from the configuration file)
  int chL[16];
  int chR[16];
  for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar){
      chL[iBar] = channelMapping[iBar*2+0]+32*DUTasic;
      chR[iBar] = channelMapping[iBar*2+1]+32*DUTasic;
  }
  
  // - define branches
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
  
  // - get plot settings from the cfg file
  std::vector<float> Vov = opts.GetOpt<std::vector<float> >("Plots.Vov");
  std::vector<int> energyBins = opts.GetOpt<std::vector<int> >("Plots.energyBins");
  std::vector<int> energyMins = opts.GetOpt<std::vector<int> >("Plots.energyMins");
  std::vector<int> energyMaxs = opts.GetOpt<std::vector<int> >("Plots.energyMaxs");
  std::map<float,int> map_energyBins;
  std::map<float,int> map_energyMins;
  std::map<float,int> map_energyMaxs;
  for(unsigned int ii = 0; ii < Vov.size(); ++ii){
    map_energyBins[Vov[ii]] = energyBins[ii];
    map_energyMins[Vov[ii]] = energyMins[ii];
    map_energyMaxs[Vov[ii]] = energyMaxs[ii];
  }
  
  // - read minimum energy for each bar from the minEnergies cfg file
  int vetoOtherBars = opts.GetOpt<int>("Cuts.vetoOtherBars");
  std::string minEnergiesFileName = opts.GetOpt<std::string>("Cuts.minEnergiesFileName");
  std::map < std::pair<int, float>, float> minE;
  if (minEnergiesFileName != "") {
    std::ifstream minEnergiesFile;
    minEnergiesFile.open(minEnergiesFileName);
    std::string line;
    int bar;
    float ov;
    float value;
    while (getline(minEnergiesFile, line)) {
      if (line.empty()) continue;
      std::istringstream ss(line);
      ss >> bar >> ov >> value;
      minE[std::make_pair(bar,ov)] = value;
    }
    for (unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar) {
      auto key = std::make_pair(iBar, my_step1);
      if (minE.find(key) == minE.end()) {
	std::cerr << "[WARNING] Missing values in minimum energy values file for bar  "<<iBar << "  Vov  " <<my_step1 <<" --> setting to default"<< std::endl;
	minE[key] = map_energyMins[my_step1]; }
    }
  }
  else{
    for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar){
      minE[std::make_pair(iBar, my_step1)] = map_energyMins[my_step1];
    }
  }
  
  // - define histograms
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
  std::map<int,bool> acceptEvent;

  // - Coincidence pre loop
  if( !opts.GetOpt<std::string>("Coincidence.status").compare("yes") &&
      ( !opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar") ||
	!opts.GetOpt<std::string>("Input.sourceName").compare("Na22") ||
	!opts.GetOpt<std::string>("Input.sourceName").compare("TB") ||
	!opts.GetOpt<std::string>("Input.sourceName").compare("keepAll") ) )
    {
      float energyL_ext;
      float energyR_ext;
      int nEntries = tree->GetEntries();
      if( maxEntries > 0 ) nEntries = maxEntries;

      // -- loop on entries 
      for(int entry = 0; entry < nEntries; ++entry){
	tree -> GetEntry(entry);
	acceptEvent[entry] = false;
	if( entry%500000 == 0 ){
	  std::cout << "\n>>> external bar loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << std::endl;
	  TrackProcess(cpu, mem, vsz, rss);
	}
	float Vov = step1;
	// --- Decode TOFHIR thresholds from the step2 value
	//     - step2 encodes 3 thresh values
	//     - step2 = 10000*(vth1+1) + 100*(vth2+1) + (vthE+1)
        float vth1 = float(int(step2/10000)-1);
        float vth2 = float(int((step2-10000*(vth1+1))/100.)-1);
        float vthE = float(int((step2-10000*(vth1+1)) - 100*(vth2+1))-1);
        float vth = 0.;
	if( entry%500000 == 0 ){
	  std::cout << step2 << " ith1: " << vth1 << " ith2: " << vth2 << " E: " << vthE << std::endl;
	}
	if(!opts.GetOpt<std::string>("Input.vth").compare("vth1"))  { vth = vth1;}
	if(!opts.GetOpt<std::string>("Input.vth").compare("vth2"))  { vth = vth2;}
	if(!opts.GetOpt<std::string>("Input.vth").compare("vthE"))  { vth = vthE;}

	// --- select only one OV
	if (my_step1 > 0  && my_step1 != step1) continue;
	if (channelIdx[chL_ext] <0 || channelIdx[chR_ext] <0) continue;
	
	// --- remove showering and cross talk events
	if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB")){
	  float energySumArray = 0.;
	  int   nActiveBarsArray = 0;	  
	  for(int iBar = 0; iBar < int(channelMapping.size())/2; ++iBar) {
	    int chL_iext = channelMapping[iBar*2+0] + REFasic*32;
	    int chR_iext = channelMapping[iBar*2+1] + REFasic*32;
	    float energyL_iext = (*energy)[channelIdx[chL_iext]];              
	    float energyR_iext = (*energy)[channelIdx[chR_iext]]; 
	    float totL_iext    = 0.001*(*tot)[channelIdx[chL_iext]];              
	    float totR_iext    = 0.001*(*tot)[channelIdx[chR_iext]]; 
	    if ( totL_iext > -10 && totL_iext < 100 && totR_iext > -10 && totR_iext < 100   ){
	      float energyMean=(energyL_iext+energyR_iext)/2;
	      if (energyMean>0){
		energySumArray+=energyMean;
		nActiveBarsArray+=1;
	      }
	    }
	  }
	  if (nActiveBarsArray > maxActiveBars ) continue;
	}
	energyL_ext = (*energy)[channelIdx[chL_ext]];
	energyR_ext = (*energy)[channelIdx[chR_ext]];
	if (interCalibrationFile_ref != "0")
	  {
	    applyCalibration(energyL_ext, barL_ext, "L", Vov_ext, th_ext, calibMap_ref);
	    applyCalibration(energyR_ext, barR_ext, "R", Vov_ext, th_ext, calibMap_ref);
	  }

	// --- encode conf parameters into an index
	//     - Vov scaled by 100 and stored in the 10^4 position
	//     - vth stored in the 10^2 position
	//     - constant offset (99) to reserve last two digits
	int index( (10000*int(Vov*100.)) + (100*vth) + 99 );
	
	// --- create histograms for the external bar used for requiring coincidence with the DUT
	if( h1_energyLR_ext[index] == NULL ){
	  c[index] = new TCanvas(Form("c1_Vov%.2f_th%02.0f",Vov,vth), Form("c1_Vov%.2f_th%02.0f",Vov,vth));
	  c[index] -> cd();
	  c[index] ->SetLogy();
	  h1_energyLR_ext[index] = new TH1F(Form("h1_energy_external_barL-R_Vov%.2f_th%02.0f",Vov,vth),"",map_energyBins[Vov],map_energyMins[Vov],map_energyMaxs[Vov]);
	}
	acceptEvent[entry] = true;
	h1_energyLR_ext[index] -> Fill(0.5*(energyL_ext + energyR_ext));
      }
      std::cout << std::endl;

      // -- analysis for Na22 source
      if (!opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar") || !opts.GetOpt<std::string>("Input.sourceName").compare("Na22")){		
	for( auto index : h1_energyLR_ext){
	  rangesLR[index.first] = new std::vector<float>;
	  rangesLR[index.first]->push_back(30);
	  rangesLR[index.first]->push_back(950);	  
	  if(!opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar")) Na22SpectrumAnalyzerSingleBar(index.second,rangesLR[index.first]);
	  if(!opts.GetOpt<std::string>("Input.sourceName").compare("Na22")) Na22SpectrumAnalyzer(index.second,rangesLR[index.first]);
	}
      }

      // -- energy selections for the TB DUT spectra (still inside the coincidence loop)
      if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB")){
	for( auto index : h1_energyLR_ext){
	  // --- rangesLR is a map having the index encoding thresholds as key and the energy ranges as items
	  rangesLR[index.first] = new std::vector<float>;

	  // --- same encoding as above
	  float Vov = float ((int(index.first /10000))/100.);
	  float vth = float(int((index.first-Vov*10000*100)/100.));

	  // --- define energy ranges 
	  // FIXME: is it worth it to use minE ?
	  index.second->GetXaxis()->SetRangeUser(200,900);
	  float max = index.second->GetBinCenter(index.second->GetMaximumBin());
	  index.second->GetXaxis()->SetRangeUser(0,1024);

	  // --- fit the energy spectrum with a Landau function
	  TF1* f_pre = new TF1(Form("fit_energy_coincBar_Vov%.2f_vth_%02.0f",Vov,vth), "[0]*TMath::Landau(x,[1],[2])", 0, 1000.);
	  f_pre -> SetRange(max*0.70, max*1.5);
	  f_pre -> SetLineColor(kBlack);
          f_pre -> SetLineWidth(2);
          f_pre -> SetParameters(index.second->Integral(index.second->GetMaximumBin(), index.second->GetNbinsX())/10, max, 0.1*max);
	  f_pre -> SetParLimits(1, 0, 9999);
	  f_pre -> SetParLimits(2, 0, 9999);
	  index.second->Fit(f_pre, "QRS+");
	  if (f_pre->GetParameter(1)>20)
	    rangesLR[index.first] -> push_back( 0.70*f_pre->GetParameter(1));
	  else
	    rangesLR[index.first] -> push_back( 20 );
	  rangesLR[index.first] -> push_back( 950 );	  
	  std::cout << "Vov = " << Vov << "  vth = " << vth
		    << "    Coincidence bar - energy range:  " << rangesLR[index.first]->at(0) << " - " << rangesLR[index.first]->at(1)<< std::endl;
	}
      }
    } 

  // ------------------------
  // - 1st loop over events  
  ModuleEventWithRefClass anEvent;

  // - reference module info (stored only if saveRefInfo 1 in config)
  unsigned short qfineL_ref = -10;
  unsigned short qfineR_ref = -10;
  float totL_ref = -10;
  float totR_ref = -10;
  long long timeL_ref = -10;
  long long timeR_ref = -10;
  long long timeL_ref_cor = -10;
  long long timeR_ref_cor = -10;
  unsigned short t1fineL_ref = -10;
  unsigned short t1fineR_ref = -10; 
  float energyL_ref = -10;
  float energyR_ref = -10;

  // - module under test info
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

  // - loop over entries
  for(int entry = 0; entry < nEntries; ++entry) {
    tree -> GetEntry(entry);
    if( entry%500000 == 0 )
      {
	std::cout << "\n>>> 1st loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << std::endl;
	TrackProcess(cpu, mem, vsz, rss);
      }

    // -- if useTrackInfo, apply position selections on wire chambers info
    if (useTrackInfo && nhits > 0 &&  (x < -100 || y < -100 ) ) continue;

    // -- decode step2 values
    float Vov = step1;
    float vth1 = float(int(step2/10000)-1);
    float vth2 = float(int((step2-10000*(vth1+1))/100.)-1);
    float vthE = float(int((step2-10000*(vth1+1)) - 100*(vth2+1))-1);
    float vth = 0;
    std::string vthMode = opts.GetOpt<std::string>("Input.vth");
    if(!opts.GetOpt<std::string>("Input.vth").compare("vth1"))  { vth = vth1;}
    if(!opts.GetOpt<std::string>("Input.vth").compare("vth2"))  { vth = vth2;}
    if(!opts.GetOpt<std::string>("Input.vth").compare("vthE"))  { vth = vthE;}
    if( entry%500000 == 0 ){
      std::cout << step2 << " ith1: " << vth1 << " ith2: " << vth2 << " E: " << vthE << std::endl;
    }
    
    // -- select only one OV
    if (my_step1 > 0  && my_step1 != step1) continue;
    
    // -- check coincidence with another channel 
    if(!opts.GetOpt<std::string>("Coincidence.status").compare("yes"))
      {
	if(!acceptEvent[entry] ) continue;	
	float energyL_ext = (*energy)[channelIdx[chL_ext]];
	float energyR_ext = (*energy)[channelIdx[chR_ext]];
	if (interCalibrationFile_ref != "0")
	  {
	    applyCalibration(energyL_ext, barL_ext, "L", Vov_ext, th_ext, calibMap_ref);
	    applyCalibration(energyR_ext, barR_ext, "R", Vov_ext, th_ext, calibMap_ref);
	  }	
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
    
    // -- store info for each bar
    for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar){
      if (channelIdx[chL[iBar]] >=0 && channelIdx[chR[iBar]] >=0){	  
        qfineL[iBar]=(*qfine)[channelIdx[chL[iBar]]];
	qfineR[iBar]=(*qfine)[channelIdx[chR[iBar]]];
	totL[iBar]=0.001*(*tot)[channelIdx[chL[iBar]]];
	totR[iBar]=0.001*(*tot)[channelIdx[chR[iBar]]];
	energyL[iBar]=(*energy)[channelIdx[chL[iBar]]];
	energyR[iBar]=(*energy)[channelIdx[chR[iBar]]];
	
	// --- if intercalibraion file specified in the config, apply corrections to the DUT energies
	if (interCalibrationFile != "0")
	  {
	    applyCalibration(energyL[iBar], iBar, "L", Vov, vth, calibMap);
	    applyCalibration(energyR[iBar], iBar, "R", Vov, vth, calibMap);
	  }
	timeL[iBar]=(*time)[channelIdx[chL[iBar]]];
	timeR[iBar]=(*time)[channelIdx[chR[iBar]]];
	t1fineL[iBar]=(*t1fine)[channelIdx[chL[iBar]]];
	t1fineR[iBar]=(*t1fine)[channelIdx[chR[iBar]]];
	if (saveRefInfo)
	  {
	    qfineL_ref=(*qfine)[channelIdx[chL_ext]];
	    qfineR_ref=(*qfine)[channelIdx[chR_ext]];
	    totL_ref=0.001*(*tot)[channelIdx[chL_ext]];
	    totR_ref=0.001*(*tot)[channelIdx[chR_ext]];
	    energyL_ref=(*energy)[channelIdx[chL_ext]];
	    energyR_ref=(*energy)[channelIdx[chR_ext]];

	    // --- if intercalibraion file specified in the config, apply corrections to the REF energies
	    if (interCalibrationFile_ref != "0")
	      {
		applyCalibration(energyL_ref, barL_ext, "L", Vov_ext, th_ext, calibMap_ref);
		applyCalibration(energyR_ref, barR_ext, "R", Vov_ext, th_ext, calibMap_ref);
	      }
	    
	    timeL_ref=(*time)[channelIdx[chL_ext]];
	    timeR_ref=(*time)[channelIdx[chR_ext]];
	    t1fineL_ref=(*t1fine)[channelIdx[chL_ext]];
	    t1fineR_ref=(*t1fine)[channelIdx[chR_ext]];
	  }
	}
      else
	{
	  qfineL[iBar]=-10;
	  qfineR[iBar]=-10;
	  totL[iBar]=-9999;
	  totR[iBar]=-9999;
	  energyL[iBar]=-10;
	  energyR[iBar]=-10;
	  timeL[iBar]=-10;
	  timeR[iBar]=-10;
	  t1fineL[iBar]=-10;
	  t1fineR[iBar]=-10;
	} 
    } // end loop over bars

    int maxEn=0;
    int maxBar=0;
    float energySumArray = 0;
    int nActiveBarsArray = 0;
    int nBarsVeto[16];

    // -- determine DUT active bars
    for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar) {
      nBarsVeto[iBar] = 0;
      if (totL[iBar]>-10 && totR[iBar]>-10 && totL[iBar]<100 && totR[iBar]<100) {
	float energyMean=(energyL[iBar]+energyR[iBar])/2;
	if (energyL[iBar]>0 && energyR[iBar]>0 && energyMean > 0){
	  energySumArray+=energyMean;
	  nActiveBarsArray+=1;
	}	    

	// --- check energy in adjacent bars -----> this is currently unused but might be useful for specific studies
	for (int jBar = int(iBar) - 2; jBar < int(iBar) + 3; ++jBar){
	  if (jBar == int(iBar)) continue;
	  if (jBar < 0 || jBar > 15 ) continue;
	  if (totL[jBar]<-10 || totL[jBar]>100) continue;
	  if (totR[jBar]<-10 || totR[jBar]>100) continue;
	  float en = (energyL[jBar]+energyR[jBar])/2;
	  if ( en > minE[std::make_pair(jBar, Vov)] && minE[std::make_pair(jBar, Vov)]>1 && en<1024 ){
	    nBarsVeto[iBar]+=1;
	  }
	}
	// --- find bar having maximum average energy LR
	if(energyMean>maxEn){
	  maxEn = energyMean;
	  maxBar = iBar;
	}
      }
    } // end loop over bars

    // -- fill histograms and branch info
    for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar){
      if (totL[iBar]>-10 && totR[iBar]>-10 && totL[iBar]<100 && totR[iBar]<100){  
	int index( (10000*int(Vov*100.)) + (100*vth) + iBar );	
	
	// --- create histograms
	if( h1_totL[index] == NULL )
	{
	  h1_qfineL[index] = new TH1F(Form("h1_qfine_bar%02dL_Vov%.2f_th%02.0f",iBar,Vov,vth),"",512,-0.5,511.5);
	  h1_qfineR[index] = new TH1F(Form("h1_qfine_bar%02dR_Vov%.2f_th%02.0f",iBar,Vov,vth),"",512,-0.5,511.5);
	  h1_totL[index] = new TH1F(Form("h1_tot_bar%02dL_Vov%.2f_th%02.0f",iBar,Vov,vth),"",400,-5.,35.);
	  h1_totR[index] = new TH1F(Form("h1_tot_bar%02dR_Vov%.2f_th%02.0f",iBar,Vov,vth),"",400,-5.,35.);
	  h1_energyL[index] = new TH1F(Form("h1_energy_bar%02dL_Vov%.2f_th%02.0f",iBar,Vov,vth),"",map_energyBins[Vov],map_energyMins[Vov],map_energyMaxs[Vov]);
	  h1_energyR[index] = new TH1F(Form("h1_energy_bar%02dR_Vov%.2f_th%02.0f",iBar,Vov,vth),"",map_energyBins[Vov],map_energyMins[Vov],map_energyMaxs[Vov]);
	  outTrees[index] = new TTree(Form("data_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth),Form("data_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth));
	  outTrees[index] -> Branch("event",&anEvent);
	  h1_energyLR[index] = new TH1F(Form("h1_energy_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth),"",map_energyBins[Vov],map_energyMins[Vov],map_energyMaxs[Vov]);
	}

	// --- fill histograms for each bar for Co60 & TB analysis
	if( !opts.GetOpt<std::string>("Input.sourceName").compare("Co60") ||
	    !opts.GetOpt<std::string>("Input.sourceName").compare("Co60SumPeak") ||
	    !opts.GetOpt<std::string>("Input.sourceName").compare("TB")
	    )
	  {
	    if( totL[iBar] <= -10. || totR[iBar] <= -10. ) continue;
	    if( totL[iBar] >= 50. ||  totR[iBar] >= 50.) continue;
	    if( ( thrZero.GetThresholdZero(chL[iBar],vthMode) + vth) > 63. ) continue;
	    if( ( thrZero.GetThresholdZero(chR[iBar],vthMode) + vth) > 63. ) continue;  

	    // ---- remove showering events
	    if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB") && (vetoOtherBars && nActiveBarsArray > maxActiveBars)) continue;	  
	    h1_qfineL[index] -> Fill( qfineL[iBar] );
	    h1_totL[index] -> Fill( totL[iBar]  );
	    h1_energyL[index] -> Fill( energyL[iBar] );
	    h1_qfineR[index] -> Fill( qfineR[iBar] );
	    h1_totR[index] -> Fill( totR[iBar] );
	    h1_energyR[index] -> Fill( energyR[iBar] );
	    h1_energyLR[index] -> Fill(0.5*(energyL[iBar]+energyR[iBar]));
	    anEvent.barID = iBar;
	    anEvent.Vov = Vov;
	    anEvent.vth = vth;
	    anEvent.energyL = energyL[iBar];
	    anEvent.energyR = energyR[iBar];
	    anEvent.totL = totL[iBar];
	    anEvent.totR = totR[iBar];
	    anEvent.timeL = timeL[iBar];
	    anEvent.timeR = timeR[iBar];
	    anEvent.t1fineL = t1fineL[iBar];
	    anEvent.t1fineR = t1fineR[iBar];

	    if(saveRefInfo)
	      {
		anEvent.energyL_ref = energyL_ref;
		anEvent.energyR_ref = energyR_ref;
		anEvent.totL_ref = totL_ref;
		anEvent.totR_ref = totR_ref;
		anEvent.timeL_ref = timeL_ref;
		anEvent.timeR_ref = timeR_ref;
		anEvent.timeL_ref_cor = timeL_ref_cor;
		anEvent.timeR_ref_cor = timeR_ref_cor;
		anEvent.t1fineL_ref = t1fineL_ref;
		anEvent.t1fineR_ref = t1fineR_ref;
	      }
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
      }
    } // end loop over bars
    // -- for Na22 or Laser analysis use only the bar with max energy to remove cross-talk between adjacent bars
    if( !opts.GetOpt<std::string>("Input.sourceName").compare("Na22") ||
	!opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar") ||
	!opts.GetOpt<std::string>("Input.sourceName").compare("Laser") ||
	!opts.GetOpt<std::string>("Input.sourceName").compare("keepAll") )
      {
	int index( (10000*int(Vov*100.)) + (100*vth) + maxBar );
	if( totL[maxBar] <= -10. || totR[maxBar] <= -10. ) continue;
	if( totL[maxBar] >= 50. ||  totR[maxBar] >= 50.) continue;
	if( ( thrZero.GetThresholdZero(chL[maxBar],vthMode) + vth) > 63. ) continue;
        if( ( thrZero.GetThresholdZero(chR[maxBar],vthMode) + vth) > 63. ) continue;
	
	// --- fill histograms
	h1_qfineL[index] -> Fill( qfineL[maxBar] );
	h1_totL[index] -> Fill( totL[maxBar] );
	h1_energyL[index] -> Fill( energyL[maxBar] );
	h1_qfineR[index] -> Fill( qfineR[maxBar] );
	h1_totR[index] -> Fill( totR[maxBar] );
	h1_energyR[index] -> Fill( energyR[maxBar] );
	h1_energyLR[index] -> Fill(0.5*(energyL[maxBar]+energyR[maxBar]));
	anEvent.barID = maxBar;
	anEvent.Vov = Vov;
	anEvent.vth = vth;
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
  } // end loop over events

  // summary of zombie files
  if (!zombieFiles.empty()) {
    std::cout << "\n=== ZOMBIE FILES ===" << std::endl;
    for (const auto &f : zombieFiles) {
      std::cout << f << std::endl;
    }
  }
  if (interCalibrationFile_ref != "0")
    std::cout <<"\n[WARNING] REF module calibrations are hard-coded to take Vov 3V and threshold 10. In Sep25 TB no other configurations were acquired." <<std::endl ;
  // summary output sizes
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
