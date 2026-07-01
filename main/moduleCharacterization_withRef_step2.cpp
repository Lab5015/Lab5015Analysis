#include "interface/AnalysisUtils.h"
#include "interface/Na22SpectrumAnalyzer.h"
#include "interface/Na22SpectrumAnalyzerSingleBar_TOFHIR2.h"
#include "interface/Co60SpectrumAnalyzer_2Peaks.h"
#include "interface/FitUtils.h"
#include "interface/SetTDRStyle.h"
#include "interface/referenceCharacterization_helper.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <time.h>
#include <stdio.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <algorithm>
#include <iterator>
#include <filesystem>
#include <limits>

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
#include "TLegend.h"
#include "TSpectrum.h"


// ------ MAIN --------
// --------------------
int main(int argc, char** argv)
{
  setTDRStyle();
  gStyle->SetOptStat(1111);
  gErrorIgnoreLevel = kError;
  typedef std::numeric_limits<double> dbl;
  std::cout.precision(dbl::max_digits10);
  if( argc < 2 )
  {
    std::cout << ">>> code usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }

  // - fixed range values for visualization ---- to be FIXED, temporary patch
  float minEnergyRatio = 0.5;
  float maxEnergyRatio = 1.5;
  int nbins = 400;
  int nbins2d = 100;
  int nbins2d_phase = 70;
  int nbinsProf = 50;
  float minT = 500.;
  float maxT = 2500;
  float minPhase = 100;
  float maxPhase = 1000;
  float minDeltaT = -1000;
  float maxDeltaT = 1000;
  
  // - parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  std::vector<float> Vov = opts.GetOpt<std::vector<float> >("Plots.Vov");
  std::vector<int> energyMins = opts.GetOpt<std::vector<int> >("Plots.energyMins");
  std::vector<int> energyMaxs = opts.GetOpt<std::vector<int> >("Plots.energyMaxs");
  std::string minEnergiesFileName = opts.GetOpt<std::string>("Cuts.minEnergiesFileName");
  std::string refFileName= opts.GetOpt<std::string>("Output.refFileName");
  std::vector<int> barList = opts.GetOpt<std::vector<int> >("Plots.barList");
  
  std::filesystem::create_directories(plotDir);
  const std::vector<std::string> LRLabels = {"L","R","L-R"};
  std::map<float,int> map_energyMins;
  std::map<float,int> map_energyMaxs;
  for(unsigned int ii = 0; ii < Vov.size(); ++ii) {
    map_energyMins[Vov[ii]] = energyMins[ii];
    map_energyMaxs[Vov[ii]] = energyMaxs[ii];
  }
  
  // - read minimum energy for each bar from the minEnergies config file
  std::map < std::pair<int, float>, float> minE;
  std::cout << "> Reading minimum energy from :" <<minEnergiesFileName << std::endl;
  if( minEnergiesFileName != "" )
    {
      std::ifstream minEnergiesFile;
      minEnergiesFile.open(minEnergiesFileName);
      if(!minEnergiesFile.is_open()) {
	std::cerr << "[ERROR] Cannot open file: " << minEnergiesFileName << std::endl;
	return -1; }
      std::string line;
      int bar;
      float ov;
      float value;
      while (getline(minEnergiesFile, line)) {
	if (line.empty()) continue;
	std::istringstream ss(line);
	if( !(ss >> bar >> ov >> value) ) {
	  std::cout <<"[ERROR] Corrupted line in min energy file" <<std::endl;
	  continue;
	}       
	minE[std::make_pair(bar,ov)] = value;
      }
      // -- fill missing entries
      for (unsigned int iBar = 0; iBar < 16; ++iBar) {
	for (unsigned int ii = 0; ii < Vov.size(); ++ii) {
	  auto key = std::make_pair(iBar, Vov[ii]);
	  if (minE.find(key) == minE.end()) {
	    minE[key] = map_energyMins[Vov[ii]]; }
	}
      }
    }
  else
    {
      for(unsigned int iBar = 0; iBar < 16; ++iBar)
	for(unsigned int ii = 0; ii < Vov.size(); ++ii)
	  minE[std::make_pair(iBar, Vov[ii])] = map_energyMins[Vov[ii]];
    }

  // - open step1 file
  TFile* inFile = TFile::Open(refFileName.c_str(),"READ");
  if(!inFile || inFile->IsZombie()) {
    std::cerr << "[ERROR] Cannot open input file: " << refFileName << std::endl;
    return -1; }

  // - trees and labels
  std::map<std::string,TTree*> trees;
  std::vector<std::string> stepLabels;
  std::map<std::string,float> map_Vovs;
  std::map<std::string,float> map_ths;

  // - read ranges from the ref root file - dedicated tree
  TTree* rangeTree = dynamic_cast<TTree*>(inFile->Get("ranges"));
  std::map<std::string, std::map<int,std::vector<float>*>> ranges;
  if(!rangeTree) {
    std::cerr << "[ERROR] No ranges tree found\n";
    return -1; }
  int bar;
  int vth;
  float VovR;
  char LR[16];
  float rmin;
  float rmax;
  rangeTree->SetBranchAddress("bar",&bar);
  rangeTree->SetBranchAddress("vth",&vth);
  rangeTree->SetBranchAddress("Vov",&VovR);
  rangeTree->SetBranchAddress("LR",LR);
  rangeTree->SetBranchAddress("rmin",&rmin);
  rangeTree->SetBranchAddress("rmax",&rmax);
  for(Long64_t i=0; i<rangeTree->GetEntries(); ++i) {
    rangeTree->GetEntry(i);     
    int index = (10000 * int(VovR*100.)) + (100 * vth) + bar;      
    ranges[LR][index] = new std::vector<float>{rmin,rmax}; }
  
  // - loop on ROOT keys
  TIter next(inFile->GetListOfKeys());
  TObject* object = nullptr;  
  while((object = next()))
    {
      std::string name = object->GetName();
      // -- trees
      if(name.find("data_") != std::string::npos)
	{
	  TTree* tree = dynamic_cast<TTree*>(inFile->Get(name.c_str()));
	  if(!tree) continue;
	  std::vector<std::string> tokens = GetTokens(name,'_');	  
	  if(tokens.size() < 4) continue;
	  std::string label = Form("%s_%s_%s", tokens[1].c_str(), tokens[2].c_str(), tokens[3].c_str());
	  trees[label] = tree;
	}
      // -- histograms used only to recover Vov/th labels
      if(name.find("h1_energy_b") != std::string::npos)
	{
	  std::vector<std::string> tokens = GetTokens(name,'_');
	  if(tokens.size() < 5) continue;
	  std::string stepLabel = tokens[3] + "_" + tokens[4];
	  stepLabels.push_back(stepLabel);
	  std::string string_Vov = tokens[3];
	  string_Vov.erase(0,3);
	  std::string string_th = tokens[4];
	  string_th.erase(0,2);
	  map_Vovs[stepLabel] = atof(string_Vov.c_str());
	  map_ths[stepLabel] = atof(string_th.c_str());
	}
    }
  
  // - sort and remove duplicates
  std::sort(stepLabels.begin(),stepLabels.end());
  stepLabels.erase(std::unique(stepLabels.begin(),stepLabels.end()),stepLabels.end());
  // - define output files
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileNameStep2");
  TFile* outFile = TFile::Open(outFileName.c_str(),"UPDATE");
  outFile->cd();
  
  // ========================
  // 1st loop 
  //  - build corrected reference time
  //  - apply energy correctons to the reference module time
  // =========================
  // - vars
  TF1* f;
  
  // 0a. DUT energy
  std::map<double,TH1F*> h1_eL;
  std::map<double,TH1F*> h1_eR;  
  std::map<double,TF1*> f_eL;
  std::map<double,TF1*> f_eR;

  // 0b. tdiff objects
  std::map<double,TH1F*>     h1_deltaT_tL_tR;
  std::map<double,TH1F*>     h1_eRatio;
  std::map<double,TF1*>      f_eRatio;
  std::map<double,TF1*>      f_LR;
  std::map<double,TProfile*> p1_deltaT_tL_tR_vs_eRatio;

  // 0c. phase histo
  std::map<double,TH1F*> h1_phaseAve;
  std::map<double,TH1F*> h1_phaseR;
  std::map<double,TH1F*> h1_phaseL;

  // 1. DUT - REF raw 
  // -    DUT module infos - tAverage REF 
  std::map<double,TH1F*> h1_deltaT_tL_tAveRef;
  std::map<double,TH1F*> h1_deltaT_tR_tAveRef;
  std::map<double,TH2F*> h2_deltaT_tL_tAveRef_vs_eL;
  std::map<double,TH2F*> h2_deltaT_tL_tAveRef_vs_eAveRef;

  // 2. DUT - REF_cor 
  // -    DUT module infos - tAverage REF corrected for TW
  std::map<double,TH1F*> h1_deltaT_tL_tAveRefCor;
  std::map<double,TH1F*> h1_deltaT_tR_tAveRefCor;
  std::map<double,TH2F*> h2_deltaT_tL_tAveRefCor_vs_eL;
  std::map<double,TH2F*> h2_deltaT_tL_tAveRefCor_vs_eAveRef;

  // 3. get DUT TW
  std::map<double,TProfile*> p1_deltaT_tL_tAveRefCor_vs_eL;
  std::map<double,TProfile*> p1_deltaT_tR_tAveRefCor_vs_eR;
  std::map<double,TF1*> f_LAve_eL;
  std::map<double,TF1*> f_RAve_eR;

  std::unordered_set<int> barSet(barList.begin(), barList.end());
  std::map<int,std::map<int,bool> > accept;
  for(auto& mapIt : trees)
    {
      ModuleEventWithRefClass* anEvent = new ModuleEventWithRefClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);
      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
	{
	  if( entry%100000 == 0 ) {
	    std::cout << ">>> 1st loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
	  }
	  mapIt.second -> GetEntry(entry);      
	  if (!barSet.count(anEvent->barID)) continue;
	  int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );
	  accept[index1][entry] = false;	  
	  int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;
	  double  index2( 10000000*energyBinAverage+index1 );

	  // ========== SELECTIONS ===========
	  double deltaTL        = anEvent->timeL_ref - (0.5*(anEvent->timeL     + anEvent->timeR));
	  double deltaTR        = anEvent->timeR_ref - (0.5*(anEvent->timeL     + anEvent->timeR));
	  double tAve_ref       = 0.5*( anEvent->timeL_ref     +  anEvent->timeR_ref     );
	  double deltaTL_AveRef = anEvent->timeL - tAve_ref;
	  double deltaTR_AveRef = anEvent->timeR - tAve_ref;
	  bool pass = PassSelection(anEvent, deltaTL, deltaTR, deltaTL_AveRef, deltaTR_AveRef, ranges, index1, anEvent->energyL, anEvent->energyR);
	  if (pass==false) continue;
	  accept[index1][entry] = true;
	  // =================================
	  
	  if( !h1_eL.count(index2) )
	    {
	      std::string labelLR_energyBin(Form("bar%02d_Vov%.02f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth,energyBinAverage));

	      // --- 0a. DUT energy
	      h1_eL[index2] = new TH1F(Form("h1_eL_%s",labelLR_energyBin.c_str()),";energy_{L};Events",nbins, 0.,1024.);
	      h1_eR[index2] = new TH1F(Form("h1_eR_%s",labelLR_energyBin.c_str()),";energy_{R};Events",nbins, 0.,1024.);
	      // --- 0a. DUT time difference
	      h1_eRatio[index2] = new TH1F(Form("h1_eRatio_%s",labelLR_energyBin.c_str()),";energy_{L}/energy_{R};Events",nbins, minEnergyRatio, maxEnergyRatio);
	      h1_deltaT_tL_tR[index2] = new TH1F(Form("h1_deltaT_tL_tR_%s",labelLR_energyBin.c_str()),";t_{L} - t_{R}; Events",nbins, minDeltaT, maxDeltaT);
	      p1_deltaT_tL_tR_vs_eRatio[index2] = new TProfile(Form("p1_deltaT_tL_tR_vs_eRatio_%s",labelLR_energyBin.c_str()),";energy_{L}/energy_{R};t_{L} - t_{R}", nbinsProf, minEnergyRatio, maxEnergyRatio);	      
	      // --- 0c. DUT phase
	      h1_phaseL[index2]   = new TH1F(Form("h1_phaseL_%s",     labelLR_energyBin.c_str()), ";phase_{L};Events"  ,nbins, minPhase, maxPhase);
	      h1_phaseR[index2]   = new TH1F(Form("h1_phaseR_%s",     labelLR_energyBin.c_str()), ";phase_{R};Events"  ,nbins, minPhase, maxPhase);
	      h1_phaseAve[index2] = new TH1F(Form("h1_phaseAve_%s",labelLR_energyBin.c_str()),    ";phase_{avg};Events",nbins, minPhase, maxPhase);

	      // --- 1. DUT - REF raw
	      h1_deltaT_tL_tAveRef[index2] = new TH1F(Form("h1_deltaT_tL_tAveRef_%s",labelLR_energyBin.c_str()),";t_{L} - t_{avg}^{REF}; Events",nbins, minT, maxT);
	      h1_deltaT_tR_tAveRef[index2] = new TH1F(Form("h1_deltaT_tR_tAveRef_%s",labelLR_energyBin.c_str()),";t_{R} - t_{avg}^{REF}; Events",nbins, minT, maxT);
	      h2_deltaT_tL_tAveRef_vs_eL[index2] = new TH2F(Form("h2_deltaT_tL_tAveRef_vs_eL_%s",labelLR_energyBin.c_str()), ";energy_{L};t_{L} - t_{avg}^{REF};Events", nbins2d, 0, 1024, nbins2d, minT, maxT);
	      h2_deltaT_tL_tAveRef_vs_eAveRef[index2] = new TH2F(Form("h2_deltaT_tL_tAveRef_vs_eAveRef_%s",labelLR_energyBin.c_str()), ";energy_{avg}^{REF};t_{L} - t_{avg}^{REF};Events", nbins2d, 0, 1024, nbins2d, minT,maxT);

	      // --- 2. DUT - REF_cor 
	      h1_deltaT_tL_tAveRefCor[index2] = new TH1F(Form("h1_deltaT_tL_tAveRefCor_%s",labelLR_energyBin.c_str()),";t_{L} - t_{avg}^{REF, cor}; Events",nbins, minT, maxT);
	      h1_deltaT_tR_tAveRefCor[index2] = new TH1F(Form("h1_deltaT_tR_tAveRefCor_%s",labelLR_energyBin.c_str()),";t_{R} - t_{avg}^{REF, cor}; Events",nbins, minT, maxT);
	      h2_deltaT_tL_tAveRefCor_vs_eL[index2] = new TH2F(Form("h2_deltaT_tL_tAveRefCor_vs_eL_%s",labelLR_energyBin.c_str()), ";energy_{L};t_{L} - t_{avg}^{REF,cor};Events", nbins2d, 0, 1024, nbins2d, minT, maxT);
	      h2_deltaT_tL_tAveRefCor_vs_eAveRef[index2] = new TH2F(Form("h2_deltaT_tL_tAveRefCor_vs_eAveRef_%s",labelLR_energyBin.c_str()), ";energy_{avg}^{REF};t_{L} - t_{avg}^{REF,cor};Events", nbins2d, 0, 1024, nbins2d, minT, maxT);

	      // --- 3. get DUT TW 	      
	      p1_deltaT_tL_tAveRefCor_vs_eL[index2] = new TProfile(Form("p1_deltaT_tL_tAveRefCor_vs_eL_%s",labelLR_energyBin.c_str()),";energy_{L};t_{L} - t_{avg}^{REF, cor}",nbinsProf,0,1024);
	      p1_deltaT_tR_tAveRefCor_vs_eR[index2] = new TProfile(Form("p1_deltaT_tR_tAveRefCor_vs_eR_%s",labelLR_energyBin.c_str()),";energy_{R};t_{R} - t_{avg}^{REF, cor}",nbinsProf,0,1024);
	    }

	  // -- average quantitites
	  double tAve_ref_cor    = 0.5*( anEvent->timeL_ref_cor +  anEvent->timeR_ref_cor );
	  double energyRefAve    = 0.5*( anEvent->energyL_ref   +  anEvent->energyR_ref   );  

	  // -- 2. DUT time - REFavg_cor
	  double deltaTL_AveRefCor = anEvent->timeL - tAve_ref_cor;
	  double deltaTR_AveRefCor = anEvent->timeR - tAve_ref_cor;

	  // -- fill histograms --
	  // -- 0a. DUT energy
	  h1_eL[index2] -> Fill( anEvent->energyL );
	  h1_eR[index2] -> Fill( anEvent->energyR );
	  // -- 0b. DUT time difference
	  h1_eRatio[index2] -> Fill(anEvent->energyL/anEvent->energyR);
	  h1_deltaT_tL_tR[index2] -> Fill(anEvent->timeL - anEvent->timeR);
	  p1_deltaT_tL_tR_vs_eRatio[index2] -> Fill(anEvent->energyL/anEvent->energyR, anEvent->timeL - anEvent->timeR);
	  // -- 0c. phase
	  h1_phaseL[index2]   -> Fill( anEvent->t1fineL );
	  h1_phaseR[index2]   -> Fill( anEvent->t1fineR );
	  h1_phaseAve[index2] -> Fill( 0.5*(anEvent->t1fineL + anEvent->t1fineR) );

	  // -- 1. DUT time - REFavg_raw
	  h1_deltaT_tL_tAveRef[index2] -> Fill( deltaTL_AveRef );
	  h1_deltaT_tR_tAveRef[index2] -> Fill( deltaTR_AveRef );
	  h2_deltaT_tL_tAveRef_vs_eL[index2]       -> Fill( anEvent->energyL,  anEvent->timeL - tAve_ref);
	  h2_deltaT_tL_tAveRef_vs_eAveRef[index2]  -> Fill( energyRefAve,      anEvent->timeL - tAve_ref);

	  // -- 2. DUT time - REFavg_cor
	  h1_deltaT_tL_tAveRefCor[index2] -> Fill( deltaTL_AveRefCor );
	  h1_deltaT_tR_tAveRefCor[index2] -> Fill( deltaTR_AveRefCor );	 
	  h2_deltaT_tL_tAveRefCor_vs_eL[index2]       -> Fill( anEvent->energyL,  deltaTL_AveRefCor);
	  h2_deltaT_tL_tAveRefCor_vs_eAveRef[index2]  -> Fill( energyRefAve,      deltaTL_AveRefCor);

	  // --- 3. get DUT TW 
	  p1_deltaT_tL_tAveRefCor_vs_eL[index2] -> Fill( anEvent->energyL, deltaTL_AveRefCor );
	  p1_deltaT_tR_tAveRefCor_vs_eR[index2] -> Fill( anEvent->energyR, deltaTR_AveRefCor );
	} // end loop over entries
    }
  
  // - draw 1st loop plots
  for(auto& it : h1_eL)
    {      
      double index2 = it.first;
      outFile->cd();
      // -- 0a. DUT energy
      f_eL[index2] = FitAndSaveHisto(outFile, h1_eL[index2], plotDir,1.0,0.5, 1, kLandau);
      f_eR[index2] = FitAndSaveHisto(outFile, h1_eR[index2], plotDir,1.0,0.5, 1, kLandau);
      // -- 0b. DUT time difference
      f_eRatio[index2] = FitAndSaveHisto(outFile, h1_eRatio[index2], plotDir);
      f = FitAndSaveHisto(outFile,h1_deltaT_tL_tR[index2], plotDir);
      f_LR[index2] = FitAndSaveProfile(outFile, p1_deltaT_tL_tR_vs_eRatio[index2], plotDir);
      // -- 0c. phase
      SaveHistoToCanvas(outFile, h1_phaseL[index2], plotDir);
      SaveHistoToCanvas(outFile, h1_phaseR[index2], plotDir);
      SaveHistoToCanvas(outFile, h1_phaseAve[index2], plotDir);

      // -- 1. DUT time - REFavg_raw
      f = FitAndSaveHisto(outFile, h1_deltaT_tL_tAveRef[index2], plotDir);
      f = FitAndSaveHisto(outFile, h1_deltaT_tR_tAveRef[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRef_vs_eL[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRef_vs_eAveRef[index2], plotDir);
      
      // -- 2. DUT time - REFavg_cor
      f = FitAndSaveHisto(outFile, h1_deltaT_tL_tAveRefCor[index2], plotDir);
      f = FitAndSaveHisto(outFile, h1_deltaT_tR_tAveRefCor[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRefCor_vs_eL[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRefCor_vs_eAveRef[index2], plotDir);
      
      // --- 3. get DUT TW
      p1_deltaT_tL_tAveRefCor_vs_eL[index2]->GetYaxis()->SetRangeUser(minT, maxT);
      p1_deltaT_tR_tAveRefCor_vs_eR[index2]->GetYaxis()->SetRangeUser(minT, maxT);
      f_LAve_eL[index2] = FitAndSaveProfile(outFile, p1_deltaT_tL_tAveRefCor_vs_eL[index2], plotDir);
      f_RAve_eR[index2]	= FitAndSaveProfile(outFile, p1_deltaT_tR_tAveRefCor_vs_eR[index2], plotDir);
    }
  
  // ========================
  // 2nd loop 
  //  - build corrected DUT time
  //  - derive phase corrections
  // =========================
  // 1. DUT - corrected REF
  std::map<double,TH1F*> h1_deltaT_tL_tAveRefCor_eCor;
  std::map<double,TH1F*> h1_deltaT_tR_tAveRefCor_eCor;

  // 2. check residual dependence on energy DUT and energy REF
  std::map<double,TH2F*> h2_deltaT_tL_tAveRefCor_eCor_vs_eL;
  std::map<double,TH2F*> h2_deltaT_tR_tAveRefCor_eCor_vs_eR;
  std::map<double,TH2F*> h2_deltaT_tL_tAveRefCor_eCor_vs_eAveRef;

  // 3. derive phase corrections 
  std::map<double,TProfile*> p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL;
  std::map<double,TProfile*> p1_deltaT_tL_tAveRefCor_eCor_vs_phaseRefAve;
  std::map<double,TProfile*> p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR;  
  std::map<double,TH2F*> h2_deltaT_tL_tAveRefCor_eCor_vs_phaseL;
  std::map<double,TH2F*> h2_deltaT_tL_tAveRefCor_eCor_vs_phaseRefAve;
  std::map<double,TH2F*> h2_deltaT_tR_tAveRefCor_eCor_vs_phaseR;  

  // 4. time difference TW corrected and derive phase corrections
  std::map<double,TH1F*>     h1_deltaT_tL_tR_eRatioCor;
  std::map<double,TH2F*>     h2_deltaT_tL_tR_eRatioCor_vs_eRatio;
  std::map<double,TProfile*> p1_deltaT_tL_tR_eRatioCor_vs_phaseMean;  
  std::map<double,TH2F*>     h2_deltaT_tL_tR_eRatioCor_vs_phaseMean;  

  for(auto mapIt : trees)
    {
      ModuleEventWithRefClass* anEvent = new ModuleEventWithRefClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);
      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
	{
	  if( entry%100000 == 0 ) {
	    std::cout << ">>> 2nd loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
	  }
	  mapIt.second -> GetEntry(entry);      
	  if (!barSet.count(anEvent->barID)) continue;
	  int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );
	  if( !accept[index1][entry] ) continue ;
	  int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;
	  double  index2( 10000000*energyBinAverage+index1 );
	  if(!f_LAve_eL[index2] || !f_RAve_eR[index2]){
	    accept[index1][entry] = false;
	    continue;
	  }	  
	  if(h1_deltaT_tL_tAveRefCor_eCor[index2] == NULL)
	    {
	      std::string labelLR_energyBin(Form("bar%02d_Vov%.02f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth,energyBinAverage));
	      // --- 1. DUT - corrected REF
	      h1_deltaT_tL_tAveRefCor_eCor[index2] = new TH1F(Form("h1_deltaT_tL_tAveRefCor_eCor_%s",labelLR_energyBin.c_str()), ";t_{L}^{cor} - t_{avg}^{REF,cor}; Events", nbins, minT, maxT);
	      h1_deltaT_tR_tAveRefCor_eCor[index2] = new TH1F(Form("h1_deltaT_tR_tAveRefCor_eCor_%s",labelLR_energyBin.c_str()), ";t_{R}^{cor} - t_{avg}^{REF,cor}; Events", nbins, minT, maxT);

	      // --- 2. check residual dependence on energy DUT and energy REF
	      h2_deltaT_tL_tAveRefCor_eCor_vs_eL[index2] = new TH2F(Form("h2_deltaT_tL_tAveRefCor_eCor_vs_eL_%s",labelLR_energyBin.c_str()), ";energy_{L};t_{L}^{cor} - t_{avg}^{REF,cor};Events", nbins2d, 0, 1024, nbins2d, minT, maxT);
	      h2_deltaT_tR_tAveRefCor_eCor_vs_eR[index2] = new TH2F(Form("h2_deltaT_tR_tAveRefCor_eCor_vs_eR_%s",labelLR_energyBin.c_str()), ";energy_{R};t_{R}^{cor} - t_{avg}^{REF,cor};Events", nbins2d, 0, 1024, nbins2d, minT, maxT);
	      h2_deltaT_tL_tAveRefCor_eCor_vs_eAveRef[index2] = new TH2F(Form("h2_deltaT_tL_tAveRefCor_eCor_vs_eAveRef_%s",labelLR_energyBin.c_str()), ";energy_{avg}^{REF};t_{L}^{cor} - t_{avg}^{REF,cor};Events", nbins2d, 0, 1024, nbins2d, minT, maxT);

	      // --- 3. derive phase corrections
	      p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2]       = new TProfile(Form("p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL_%s",labelLR_energyBin.c_str()),";phase_{L};t_{L}^{cor} - t_{avg}^{REF, cor}",nbinsProf,minPhase, maxPhase);
	      p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2]       = new TProfile(Form("p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR_%s",labelLR_energyBin.c_str()),";phase_{R};t_{R}^{cor} - t_{avg}^{REF, cor}",nbinsProf,minPhase, maxPhase);
	      p1_deltaT_tL_tAveRefCor_eCor_vs_phaseRefAve[index2]  = new TProfile(Form("p1_deltaT_tL_tAveRefCor_eCor_vs_phaseRefAve_%s",labelLR_energyBin.c_str()), ";phase_{avg}^{REF};t_{L}^{cor} - t_{avg}^{REF, cor}",nbinsProf,minPhase, maxPhase);	      
	      h2_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2]      = new TH2F(Form("h2_deltaT_tL_tAveRefCor_eCor_vs_phaseL_%s",labelLR_energyBin.c_str()),";phase_{L};t_{L}^{cor} - t_{avg}^{REF, cor};Events",nbins2d_phase,minPhase, maxPhase, nbins2d_phase, minT, maxT);
	      h2_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2]      = new TH2F(Form("h2_deltaT_tR_tAveRefCor_eCor_vs_phaseR_%s",labelLR_energyBin.c_str()),";phase_{R};t_{R}^{cor} - t_{avg}^{REF, cor};Events",nbins2d_phase,minPhase, maxPhase, nbins2d_phase, minT, maxT);
	      h2_deltaT_tL_tAveRefCor_eCor_vs_phaseRefAve[index2] = new TH2F(Form("h2_deltaT_tL_tAveRefCor_eCor_vs_phaseRefAve_%s",labelLR_energyBin.c_str()), ";phase_{avg}^{REF};t_{L}^{cor} - t_{avg}^{REF, cor}",nbins2d_phase,minPhase, maxPhase, nbins2d_phase, minT, maxT);
	      // --- 4. time difference TW corrected and derive phase corrections 
	      h1_deltaT_tL_tR_eRatioCor[index2] = new TH1F(Form("h1_deltaT_tL_tR_eRatioCor_%s",labelLR_energyBin.c_str()),";(t_{L} - t_{R})^{cor}; Events",nbins, -1000, 1000);
	      h2_deltaT_tL_tR_eRatioCor_vs_eRatio[index2] = new TH2F(Form("h2_deltaT_tL_tR_eRatioCor_vs_eRatio_%s",labelLR_energyBin.c_str()), ";energy_{L}/energy_{R};t_{L}^{cor} - t_{R}^{cor}",nbins2d, minEnergyRatio, maxEnergyRatio, nbins2d, minDeltaT, maxDeltaT);
	      p1_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2] = new TProfile(Form("p1_deltaT_tL_tR_eRatioCor_vs_phaseMean_%s",labelLR_energyBin.c_str()), ";phase_{avg};t_{L}^{cor} - t_{R}^{cor}",nbinsProf, minPhase, maxPhase);
	      h2_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2] = new TH2F(Form("h2_deltaT_tL_tR_eRatioCor_vs_phaseMean_%s",labelLR_energyBin.c_str()), ";phase_{avg};t_{L}^{cor} - t_{R}^{cor}",nbins2d_phase, minPhase, maxPhase, nbins2d_phase,minDeltaT, maxDeltaT);
	    }
	  // -- define quantities
	  double tAve_ref_cor   = 0.5*( anEvent->timeL_ref_cor +  anEvent->timeR_ref_cor );
	  double energyRefAve   = 0.5*( anEvent->energyL_ref   + anEvent->energyR_ref    );
	  double phaseRefAve    = 0.5*( anEvent->t1fineL_ref   + anEvent->t1fineR_ref    );	  	  
	  // -- apply DUT TW on DUT-REFavg
	  double eLCor         = f_LAve_eL[index2] -> Eval(anEvent->energyL) - f_LAve_eL[index2] -> Eval(f_eL[index2]->GetParameter(1));
	  double eRCor         = f_RAve_eR[index2] -> Eval(anEvent->energyR) - f_RAve_eR[index2] -> Eval(f_eR[index2]->GetParameter(1));	  
	  double deltaTL_eLCor = anEvent->timeL - tAve_ref_cor - eLCor;
	  double deltaTR_eRCor = anEvent->timeR - tAve_ref_cor - eRCor;
	  // -- apply DUT TW on DUT tDiff
	  double eRatioCor          = f_LR[index2] -> Eval( anEvent->energyL/anEvent->energyR ) - f_LR[index2] -> Eval( f_eRatio[index2]->GetParameter(1) );
	  double deltaTLR_eRatioCor = anEvent->timeL - anEvent->timeR - eRatioCor;
	  
	  // --- fill histograms
	  // 1. DUT - corrected REF
	  h1_deltaT_tL_tAveRefCor_eCor[index2] -> Fill( deltaTL_eLCor );
	  h1_deltaT_tR_tAveRefCor_eCor[index2] -> Fill( deltaTR_eRCor );

	  // 2. check residual dependence on energy DUT and energy REF
	  h2_deltaT_tL_tAveRefCor_eCor_vs_eL[index2]      -> Fill( anEvent->energyL ,  deltaTL_eLCor );
	  h2_deltaT_tR_tAveRefCor_eCor_vs_eR[index2]      -> Fill( anEvent->energyR ,  deltaTR_eRCor );
	  h2_deltaT_tL_tAveRefCor_eCor_vs_eAveRef[index2] -> Fill( energyRefAve ,      deltaTL_eLCor );

	  // 3. derive phase corrections
	  p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2] -> Fill( anEvent->t1fineL, deltaTL_eLCor );
	  p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2] -> Fill( anEvent->t1fineR, deltaTR_eRCor );
	  p1_deltaT_tL_tAveRefCor_eCor_vs_phaseRefAve[index2] -> Fill( phaseRefAve, deltaTL_eLCor );
	  h2_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2] -> Fill( anEvent->t1fineL, deltaTL_eLCor );
	  h2_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2] -> Fill( anEvent->t1fineR, deltaTR_eRCor );
	  h2_deltaT_tL_tAveRefCor_eCor_vs_phaseRefAve[index2] -> Fill( phaseRefAve, deltaTL_eLCor );

	  // 4. time difference TW corrected and derive phase corrections 
	  h1_deltaT_tL_tR_eRatioCor[index2]      -> Fill(deltaTLR_eRatioCor);
	  h2_deltaT_tL_tR_eRatioCor_vs_eRatio[index2] -> Fill(anEvent->energyL/anEvent->energyR , deltaTLR_eRatioCor);
	  p1_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2] -> Fill( 0.5*(anEvent->t1fineL + anEvent->t1fineR), deltaTLR_eRatioCor);
	  h2_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2] -> Fill( 0.5*(anEvent->t1fineL + anEvent->t1fineR), deltaTLR_eRatioCor);
	}
    }
  // - save and draw objects  2nd loop
  for(auto& it : h1_deltaT_tL_tAveRefCor_eCor)
    {
      double index2 = it.first;
      outFile->cd();
      // 1. DUT - corrected REF
      f = FitAndSaveHisto(outFile, h1_deltaT_tL_tAveRefCor_eCor[index2], plotDir);
      f = FitAndSaveHisto(outFile, h1_deltaT_tR_tAveRefCor_eCor[index2], plotDir);

      // 2. check residual dependence on energy DUT and energy REF
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRefCor_eCor_vs_eL[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tR_tAveRefCor_eCor_vs_eR[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRefCor_eCor_vs_eAveRef[index2], plotDir);

      // 3. derive phase corrections
      p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2]->GetYaxis()->SetRangeUser(minT, maxT);
      p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2]->GetYaxis()->SetRangeUser(minT, maxT);
      p1_deltaT_tL_tAveRefCor_eCor_vs_phaseRefAve[index2]->GetYaxis()->SetRangeUser(minT, maxT);
      SaveProfileToCanvas(outFile, p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2], plotDir);
      SaveProfileToCanvas(outFile, p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2], plotDir);
      SaveProfileToCanvas(outFile, p1_deltaT_tL_tAveRefCor_eCor_vs_phaseRefAve[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRefCor_eCor_vs_phaseRefAve[index2], plotDir);

      // 4. time difference TW corrected and derive phase corrections
      p1_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2]->GetYaxis()->SetRangeUser(minDeltaT, maxDeltaT);
      f = FitAndSaveHisto(outFile, h1_deltaT_tL_tR_eRatioCor[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tR_eRatioCor_vs_eRatio[index2], plotDir);
      SaveProfileToCanvas(outFile, p1_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2] , plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2], plotDir);
    }

  // ========================
  // 3rd loop
  // - apply phase corrections
  // =================

  // 1. apply phase corrections
  std::map<double,TH1F*> h1_deltaT_tL_tAveRefCor_eCor_phaseCor;
  std::map<double,TH1F*> h1_deltaT_tR_tAveRefCor_eCor_phaseCor;
  std::map<double,TH1F*> h1_deltaT_tAve_tAveRefCor_eCor_phaseCor;
  std::map<double,TH1F*> h1_deltaT_tL_tR_eRatioCor_phaseMeanCor;

  // 2. check residual dependence on phase
  std::map<double,TH2F*> h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseL;
  std::map<double,TH2F*> h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseAveRef;
  
  for(auto mapIt : trees)
    {
      ModuleEventWithRefClass* anEvent = new ModuleEventWithRefClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);
      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
	{
	  if( entry%100000 == 0 ) {
	    std::cout << ">>> 3rd loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
	  }
	  mapIt.second -> GetEntry(entry);      
	  if (!barSet.count(anEvent->barID)) continue;
	  int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );
	  if( !accept[index1][entry] ) continue ;
	  int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;
	  double  index2( 10000000*energyBinAverage+index1 );
	  if(h1_deltaT_tL_tAveRefCor_eCor_phaseCor[index2] == NULL)
	    {
	      std::string labelLR_energyBin(Form("bar%02d_Vov%.02f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth,energyBinAverage));
	      // --- 1. apply phase corrections
	      h1_deltaT_tL_tAveRefCor_eCor_phaseCor[index2] = new TH1F(Form("h1_deltaT_tL_tAveRefCor_eCor_phaseCor_%s",labelLR_energyBin.c_str()), ";t_{L}^{cor, ph} - t_{avg}^{REF,cor}; Events", nbins, minT, maxT);
	      h1_deltaT_tR_tAveRefCor_eCor_phaseCor[index2] = new TH1F(Form("h1_deltaT_tR_tAveRefCor_eCor_phaseCor_%s",labelLR_energyBin.c_str()), ";t_{R}^{cor, ph} - t_{avg}^{REF,cor}; Events", nbins, minT, maxT);
	      h1_deltaT_tAve_tAveRefCor_eCor_phaseCor[index2] = new TH1F(Form("h1_deltaT_tAve_tAveRefCor_eCor_phaseCor_%s",labelLR_energyBin.c_str()), ";t_{avg}^{cor, ph} - t_{avg}^{REF,cor}; Events", nbins, minT, maxT);
	      h1_deltaT_tL_tR_eRatioCor_phaseMeanCor[index2] = new TH1F(Form("h1_deltaT_tL_tR_eRatioCor_phaseMeanCor_%s",labelLR_energyBin.c_str()), ";t_{L}^{cor, ph} - t_{R}^{cor, ph}; Events", nbins, minDeltaT, maxDeltaT);

	      // --- 2. check residual dependence on phase
	      h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseAveRef[index2] = new TH2F(Form("h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseAveRef_%s",labelLR_energyBin.c_str()), ";phase_{avg}^{REF}; t_{L}^{cor, ph} - t_{avg}^{REF,cor}; Events", nbins2d_phase, minPhase, maxPhase, nbins2d_phase,minT,maxT);
	      h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseL[index2] = new TH2F(Form("h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseL_%s",labelLR_energyBin.c_str()), ";phase_{L}; t_{L}^{cor, ph} - t_{avg}^{REF,cor}; Events", nbins2d_phase, minPhase, maxPhase, nbins2d_phase, minT, maxT);
	    }
	  // -- define quantities
	  double tAve_ref_cor   = 0.5*( anEvent->timeL_ref_cor +  anEvent->timeR_ref_cor );
	  double phaseRefAve    = 0.5*( anEvent->t1fineL_ref   + anEvent->t1fineR_ref    );	  	  
	  // -- apply DUT TW on DUT-REFavg
	  double eLCor         = f_LAve_eL[index2] -> Eval(anEvent->energyL) - f_LAve_eL[index2] -> Eval(f_eL[index2]->GetParameter(1));
	  double eRCor         = f_RAve_eR[index2] -> Eval(anEvent->energyR) - f_RAve_eR[index2] -> Eval(f_eR[index2]->GetParameter(1));	  
	  // -- apply DUT TW on DUT tDiff
	  double eRatioCor          = f_LR[index2] -> Eval( anEvent->energyL/anEvent->energyR ) - f_LR[index2] -> Eval( f_eRatio[index2]->GetParameter(1) );

	  // -- correction DUT vs phase
	  int phaseBinL      = p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2]->FindBin( anEvent->t1fineL );
	  int phaseBinL_offs = p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2]->FindBin( h1_phaseL[index2]->GetMean());
	  double phaseLCor   = p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2]->GetBinContent(phaseBinL) - p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2]->GetBinContent(phaseBinL_offs);	  
	  int phaseBinR      = p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2]->FindBin( anEvent->t1fineR );
	  int phaseBinR_offs = p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2]->FindBin(  h1_phaseR[index2]->GetMean() );
	  double phaseRCor   = p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2]->GetBinContent(phaseBinR) - p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2]->GetBinContent(phaseBinR_offs);

	  double phaseMean      = 0.5*(anEvent->t1fineL + anEvent->t1fineR );
	  int phaseBinLR        = p1_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2]->FindBin( phaseMean );
	  int phaseBinLR_offs   = p1_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2]->FindBin( h1_phaseAve[index2]->GetMean() );
	  double phaseMeanCor   = p1_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2]->GetBinContent(phaseBinLR) - p1_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2]->GetBinContent(phaseBinLR_offs);

	  // -- DUT time L-R - TW cor - phase cor
	  double deltaTLR_eRatioCor_phaseMeanCor = anEvent->timeL - anEvent->timeR - eRatioCor - phaseMeanCor;
	  
	  // -- DUT time - REFavg - eRefCor - eCor - phaseCor
	  double deltaTL_eLCor_phaseLCor = anEvent->timeL - tAve_ref_cor - eLCor - phaseLCor;
	  double deltaTR_eRCor_phaseRCor = anEvent->timeR - tAve_ref_cor - eRCor - phaseRCor;

	  // -- (DUT avg time - DUT avg eCor - DUT avg phaseCor) - (REF avg time - REF avg eCor)
	  double deltaT_tAve_tAveRefCor_eCor_phaseCor = (0.5* (anEvent->timeL+anEvent->timeR) - 0.5*(eLCor+eRCor)  - 0.5*(phaseLCor+phaseRCor)) - tAve_ref_cor;

	  // -- fill histograms
	  // 1. apply phase corrections
	  h1_deltaT_tL_tAveRefCor_eCor_phaseCor[index2] -> Fill( deltaTL_eLCor_phaseLCor );
	  h1_deltaT_tR_tAveRefCor_eCor_phaseCor[index2] -> Fill( deltaTR_eRCor_phaseRCor );	  
	  h1_deltaT_tAve_tAveRefCor_eCor_phaseCor[index2] -> Fill( deltaT_tAve_tAveRefCor_eCor_phaseCor );
	  h1_deltaT_tL_tR_eRatioCor_phaseMeanCor[index2]  -> Fill( deltaTLR_eRatioCor_phaseMeanCor);
	  
	  // --- 2. check residual dependence on phase
	  h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseL[index2] -> Fill( anEvent->t1fineL,    deltaTL_eLCor_phaseLCor);
	  h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseAveRef[index2] -> Fill( phaseRefAve ,    deltaTL_eLCor_phaseLCor);
	} // end loop entries
    } // end loop tree

  // draw 3rd loop 
  for(auto& it : h1_deltaT_tL_tAveRefCor_eCor_phaseCor)
    {
      double index2 = it.first;
      outFile->cd();
      // 1. apply phase corrections
      f = FitAndSaveHisto(outFile, h1_deltaT_tL_tAveRefCor_eCor_phaseCor[index2], plotDir);      
      f = FitAndSaveHisto(outFile, h1_deltaT_tR_tAveRefCor_eCor_phaseCor[index2], plotDir);
      f = FitAndSaveHisto(outFile, h1_deltaT_tAve_tAveRefCor_eCor_phaseCor[index2], plotDir);
      f = FitAndSaveHisto(outFile, h1_deltaT_tL_tR_eRatioCor_phaseMeanCor[index2], plotDir);
      
      // 2. check residual dependence on phase  
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseAveRef[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseL[index2], plotDir);
    }

  int bytes = outFile -> Write();
  SortPlotsByQuantity(plotDir);
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
