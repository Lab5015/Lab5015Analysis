
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




// Questo codice mi e sfuggito di mano, ma ormai va cosi

// --------------------
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
    std::cout << ">>> moduleCharacterization_step2::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }

  // hard cut on energy
  float minEnergyRatio = 0.5;
  float maxEnergyRatio = 1.5;
  
  // - parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  // - get parameters
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  std::filesystem::create_directories(plotDir);
  std::vector<std::string> LRLabels;
  LRLabels.push_back("L");
  LRLabels.push_back("R");
  LRLabels.push_back("L-R");
  std::vector<float> Vov = opts.GetOpt<std::vector<float> >("Plots.Vov");
  std::vector<int> energyMins = opts.GetOpt<std::vector<int> >("Plots.energyMins");
  std::vector<int> energyMaxs = opts.GetOpt<std::vector<int> >("Plots.energyMaxs");
  std::map<float,int> map_energyMins;
  std::map<float,int> map_energyMaxs;
  for(unsigned int ii = 0; ii < Vov.size(); ++ii) {
    map_energyMins[Vov[ii]] = energyMins[ii];
    map_energyMaxs[Vov[ii]] = energyMaxs[ii];
  }
  
  // - read minimum energy for each bar from the minEnergies config file
  std::string minEnergiesFileName = opts.GetOpt<std::string>("Cuts.minEnergiesFileName");
  std::map < std::pair<int, float>, float> minE;
  std::cout << "> Reading minimum energy from :" <<minEnergiesFileName << std::endl;
  if( minEnergiesFileName != "" )
    {
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
  std::string step1FileName= opts.GetOpt<std::string>("Input.step1FileName");
  TFile* inFile = TFile::Open(step1FileName.c_str(),"READ");
  std::map<std::string,TTree*> trees;
  std::map<std::string,int> VovLabels;
  std::map<std::string,int> thLabels;
  std::vector<std::string> stepLabels;
  std::map<std::string,float> map_Vovs;
  std::map<std::string,float> map_ths;
  TList* list = inFile -> GetListOfKeys();
  TIter next(list);
  TObject* object = 0;
  
  // - loop over all the objects inside the file
  while( (object = next()) )
  {
    std::string name(object->GetName());
    std::vector<std::string> tokens = GetTokens(name,'_');
    std::size_t found;

    // -- if a tree is found, it stores the tree in the trees map
    found = name.find("data_");
    if( found!=std::string::npos )
    {
      std::string label(Form("%s_%s_%s",tokens[1].c_str(),tokens[2].c_str(),tokens[3].c_str()));
      trees[label] = (TTree*)( inFile->Get(name.c_str()) );
    }
    found = name.find("h1_energy_b");
    if( found!=std::string::npos )
    {
     // --- extract Vov e th from tokens
      std::string stepLabel = tokens[3]+"_"+tokens[4]; // Vov and threshold label
      VovLabels[tokens[3]] += 1;
      thLabels[tokens[4]] += 1;
      stepLabels.push_back(stepLabel);
      std::string string_Vov = tokens[3];
      string_Vov.erase(0,3);
      map_Vovs[stepLabel] = atof(string_Vov.c_str());
      std::string string_th = tokens[4];
      string_th.erase(0,2);
      map_ths[stepLabel] = atof(string_th.c_str());
    }
  }

  // - sort and remove duplicates
  std::sort(stepLabels.begin(),stepLabels.end());
  stepLabels.erase(std::unique(stepLabels.begin(),stepLabels.end()),stepLabels.end());
  // - define output files
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileNameStep2");
  std::string csvFileName = outFileName;
  // - replace outfilename ext
  size_t pos = csvFileName.find_last_of(".");
  if (pos != std::string::npos) {
    csvFileName = csvFileName.substr(0, pos) + ".csv";
  } else {
    csvFileName += ".csv";
  }
  std::ofstream out(csvFileName);
  out << "index,mpvL,mpvR,mpvLR,slopeL_ref,slopeR_ref,offsetL_ref,offsetR_ref,slopeL,slopeR,offsetL,offsetR,sigmaLRef,sigmaLR,sigmaRRef,sigmaRef,sigmaL,sigmaR,sigmaDutRef\n";
  TFile* outFile = TFile::Open(outFileName.c_str(),"RECREATE");
  outFile->cd();

  std::map<int, EntryData> dataMap;
  // - define histograms and TProfiles

  // -    REF module infos - tAverage DUT : time left/right REF - time average DUT
  // 0. checks and offset
  // 0a. DUT REF correlation
  std::map<double,TH2F*> h2_eL_eLRef;
  std::map<double,TH2F*> h2_eR_eRRef;

  // 0b. reference energy
  std::map<double,TH1F*> h1_eLRef;
  std::map<double,TH1F*> h1_eRRef;  
  std::map<double,TF1*> f_eLRef;
  std::map<double,TF1*> f_eRRef;

  // 1. RAW tRef - tDUT
  // 1.a distrib
  std::map<double,TH1F*> h1_deltaT_tLRef_tAve;
  std::map<double,TH1F*> h1_deltaT_tRRef_tAve;  

  // 1b. get corrections vs energy REF
  std::map<double,TProfile*> p1_deltaT_tLRef_tAve_vs_eLRef;
  std::map<double,TProfile*> p1_deltaT_tRRef_tAve_vs_eRRef;
  std::map<double,TF1*> f_LAve_eLRef;
  std::map<double,TF1*> f_RAve_eRRef;

  // 1c. corrected tRef - tDUT
  std::map<double,TH1F*> h1_deltaT_tLRef_tAve_eLRefCor;
  std::map<double,TH1F*> h1_deltaT_tRRef_tAve_eRRefCor;

  // 1d. check residual dependence on energy DUT and energy REF
  std::map<double,TH2F*> h2_deltaT_tLRef_tAve_eRefCor_vs_eLRef;
  std::map<double,TH2F*> h2_deltaT_tLRef_tAve_eRefCor_vs_eAve;
  std::map<double,TH2F*> h2_deltaT_tRRef_tAve_eRefCor_vs_eAve;
  
  // - other maps
  std::map<std::string, std::map<int, std::vector<float>*> > ranges; //ranges[LRlabel][index]
  std::map<std::string, std::map<int, std::vector<float>*> > narrow_ranges; //ranges[LRlabel][index]
  std::map<std::string, std::map<int, float> > mpv; // mpv[LRlabel][index] 
  std::map<std::string, std::map<int, std::map<std::string,std::pair<float,float> > > > peaks;	//peaks[LRlabel][index][energyPeak]
  std::map<std::string, std::map<int, std::map<int,float> > > energyBin; // energyBin[LRlabel][index]
  std::map<int,TF1*>  f_landau; 

  // - vars
  TH1F* histo;
  TF1* f;

  int nbins = 400;
  int nbins2d = 100;
  int nbinsProf = 50;
  float minT = -3.0e3;
  float maxT = 0.;
  
  // ====================================
  // - 1st LOOP -
  // - - analyze the energy spectra: different analysis
  //     ddepending on the source.
  //     It results with the identification of the energy
  //     ranges. 
  // ====================================
  std::string source = opts.GetOpt<std::string>("Input.sourceName");
  std::string TB = "TB";
  std::vector<int> barList = opts.GetOpt<std::vector<int> >("Plots.barList");// list of bars to be analyzed read from cfg
  // - loop over the Vov_th labels
  for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float vth = map_ths[stepLabel];

      // ================ SPEED UP : only one threshold ===================
      if(vth!=10.) continue;
      // ==================================================================
      std::string VovLabel(Form("Vov%.2f",Vov));
      std::string thLabel(Form("th%02.0f",vth));
      std::cout << ">>> 1st loop: Vov " <<VovLabel <<" th " <<thLabel << std::flush;

      // selecting external bar core energy
      histo = (TH1F*)( inFile->Get(Form("h1_energy_external_barL-R_Vov%.2f_th%02.0f", Vov, vth)));
      if (!histo) {
	std::cerr << "[ERROR] missing external histogram for Vov=" << Vov << " th=" << vth << std::endl;
	continue;
      }

      // -- loop over DUT bars
      for(int iBar = 0; iBar < 16; ++iBar) {
	std::cout <<"bar "<< iBar<<std::endl;
	bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
	if (!barFound) continue;
	int index( (10000*int(Vov*100.)) + (100*vth) + iBar );

	// --- loop over L, R, LR
	for(auto LRLabel : LRLabels ) {	  
	  std::string label(Form("bar%02d%s_%s",iBar,LRLabel.c_str(),stepLabel.c_str()));	  
	  histo = (TH1F*)( inFile->Get(Form("h1_energy_%s",label.c_str())) );
	  if( !histo ) continue;

	  ranges[LRLabel][index] = new std::vector<float>;
	  narrow_ranges[LRLabel][index] = new std::vector<float>;
	  mpv[LRLabel][index] = 0.0;
	  
	  // ---- if test beam (MIP peak), peaks are identified with a landau fit
	  if(!source.compare(TB)){ 
	    float max = histo->GetBinCenter(histo->GetMaximumBin());
	    histo->GetXaxis()->SetRangeUser(minE[std::make_pair(iBar, Vov)], 950); // minE is set in the minEnergies config file to avoid fitting noise	    

	    // ----- landau fit to determine optimal range for event selections
	    f_landau[index] = new TF1(Form("f_landau_bar%02d%s_Vov%.2f_vth_%02.0f", iBar,LRLabel.c_str(),Vov,vth),"[0]*TMath::Landau(x,[1],[2])", 0,1000.);
	    float xmin = max * 0.65;
	    float xmax = std::min(max*2.5, 940.);
	    f_landau[index] -> SetRange(xmin,xmax);
	    f_landau[index] -> SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, max, 0.1*max);
	    f_landau[index] -> SetParLimits(1,0,9999);
	    f_landau[index] -> SetParLimits(2,0,9999);
	    histo -> Fit(f_landau[index],"QRS");
	    if ( f_landau[index]->GetParameter(1) > 0 ){
	      xmin = f_landau[index]->GetParameter(1) - 2 * std::abs(f_landau[index]->GetParameter(2));
	      if (xmin < minE[std::make_pair(iBar, Vov)]) xmin = minE[std::make_pair(iBar, Vov)] ;
	      xmax = std::min(f_landau[index]->GetParameter(1) * 2.5, 940.);
	      f_landau[index] -> SetRange(xmin, xmax);
	      f_landau[index] -> SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, f_landau[index]->GetParameter(1), 0.1*f_landau[index]->GetParameter(1));
	    }
	    histo -> Fit(f_landau[index],"QRS");
	    f_landau[index] -> SetLineColor(kBlack);
	    f_landau[index] -> SetLineWidth(2);
	    f_landau[index] -> Draw("same");
	    if ( f_landau[index]->GetNDF() >0 && f_landau[index]->GetParameter(1) > minE[std::make_pair(iBar, Vov)] &&
		 (f_landau[index]->GetParameter(1) - 2.0 * std::abs(f_landau[index]->GetParameter(2))) >=  minE[std::make_pair(iBar, Vov)] &&
		 (f_landau[index]->GetParameter(1) - 2.0 * std::abs(f_landau[index]->GetParameter(2))) < 950) {
	      ranges[LRLabel][index] -> push_back( f_landau[index]->GetParameter(1) - 2.0 * std::abs(f_landau[index]->GetParameter(2)));
	      mpv[LRLabel][index] = f_landau[index]->GetParameter(1);
	      narrow_ranges[LRLabel][index] -> push_back( f_landau[index]->GetParameter(1) - 1.0 * std::abs(f_landau[index]->GetParameter(2)));
	    }
	    else
	      ranges[LRLabel][index] -> push_back( minE[std::make_pair(iBar, Vov)] );
	    
	    // ----- set the energy maximum value to maximum ADC 
	    ranges[LRLabel][index] -> push_back( 940 );
	    narrow_ranges[LRLabel][index] -> push_back(f_landau[index]->GetParameter(1) + 1.0 * std::abs(f_landau[index]->GetParameter(2)));
	    // ----- store energy mean of each bin in the energyBin map
	    GetEnergyBins(histo, ranges[LRLabel][index], energyBin[LRLabel][index]);
	  }// end MIP (TB)
	  
	  // ---- draw energy plots
	  histo->GetXaxis()->SetRangeUser(0,1024);
	  outFile -> cd();
	  histo->Write();
	}// ---- end loop over L, R, L-R labels	
      }// --- end loop over bars
    } // -- end loop over stepLabels
  // -  end 1st plots  

  // ====================================
  // - 2nd LOOP -
  //   fill ref histos
  //   selecting energy average DUT
  // ====================================
  std::unordered_set<int> barSet(barList.begin(), barList.end());
  std::map<int,std::map<int,bool> > accept;
  
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

	  // --- only bars in the barList specified in the config are accepted
	  if (!barSet.count(anEvent->barID)) continue;
	  int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );
	  accept[index1][entry] = false;

	  // --- define deltaT
	  //       REF - DUTavg
	  double deltaTL        = anEvent->timeL_ref - (0.5*(anEvent->timeL     + anEvent->timeR));
	  double deltaTR        = anEvent->timeR_ref - (0.5*(anEvent->timeL     + anEvent->timeR));
	  //       DUT - REFavg
	  double deltaTL_AveRef = anEvent->timeL     - (0.5*(anEvent->timeL_ref + anEvent->timeR_ref));
	  double deltaTR_AveRef = anEvent->timeR     - (0.5*(anEvent->timeL_ref + anEvent->timeR_ref));
	  float energyAve = 0.5*(anEvent->energyL + anEvent->energyR);

	  // ========== SELECTIONS ===========
	  // ================ SPEED UP : only one threshold ===================
	  if(anEvent->vth!=10.) continue;
	  // ==================================================================	 
	  if( std::abs(deltaTL)>5000 || std::abs(deltaTR)>5000 || std::abs(deltaTL_AveRef)>5000 || std::abs(deltaTR_AveRef)>5000) continue;
	  if( anEvent->totL < 0 || anEvent->totL > 25 || anEvent->totR < 0 || anEvent->totR > 25 || anEvent->totL_ref < 0 || anEvent->totL_ref > 25 || anEvent->totR_ref < 0 || anEvent->totR_ref > 25) continue;
	  if(!ranges["L-R"][index1] || !ranges["L"][index1] || !ranges["R"][index1]) continue;
	  int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;
	  if( energyBinAverage < 1 ) continue;
	  if( energyAve<ranges["L-R"][index1]->at(0) || energyAve>ranges["L-R"][index1]->at(1)) continue;
	  accept[index1][entry] = true;
	  // --- please note that for this first loop we require the DUT events to have kinda fixed energy so that tDUT dependency on amplitude is reduced
	  if (energyAve < narrow_ranges["L-R"][index1]->at(0) || energyAve> narrow_ranges["L-R"][index1]->at(1)) continue;	  
	  // =================================
	  double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );
	  if( h1_deltaT_tLRef_tAve[index2] == NULL )
	    {
	      std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth,energyBinAverage));	      
	      h1_deltaT_tLRef_tAve[index2] = new TH1F(Form("h1_deltaT_tLRef_tAve_%s",labelLR_energyBin.c_str()),";t_{L}^{REF} - t_{avg}; Events",nbins,minT,maxT);
	      h1_deltaT_tRRef_tAve[index2] = new TH1F(Form("h1_deltaT_tRRef_tAve_%s",labelLR_energyBin.c_str()),";t_{R}^{REF} - t_{avg}; Events",nbins,minT,maxT);

	      h2_eL_eLRef[index2] = new TH2F(Form("h2_eL_eLRef_%s",labelLR_energyBin.c_str()), ";energy_{L}^{REF};energy_{L};Events",   nbins2d,0,1024, nbins2d,0,1024);
	      h2_eR_eRRef[index2] = new TH2F(Form("h2_eR_eRRef_%s",labelLR_energyBin.c_str()), ";energy_{R}^{REF};energy_{R};Events",   nbins2d,0,1024, nbins2d,0,1024);

	      h1_eLRef[index2] = new TH1F(Form("h1_eLRef_%s",labelLR_energyBin.c_str()),";energy_{L}^{REF};Events",nbins, 0.,1024.);
	      h1_eRRef[index2] = new TH1F(Form("h1_eRRef_%s",labelLR_energyBin.c_str()),";energy_{R}^{REF};Events",nbins, 0.,1024.);

	      p1_deltaT_tLRef_tAve_vs_eLRef[index2] = new TProfile(Form("p1_deltaT_tLRef_tAve_vs_eLRef_%s",labelLR_energyBin.c_str()),";energy_{L}^{REF};t_{L}^{REF} - t_{avg}",nbinsProf,0,1024);
	      p1_deltaT_tRRef_tAve_vs_eRRef[index2] = new TProfile(Form("p1_deltaT_tRRef_tAve_vs_eRRef_%s",labelLR_energyBin.c_str()),";energy_{R}^{REF};t_{R}^{REF} - t_{avg}",nbinsProf,0,1024);	      
	    }

	  // REF energy
	  h1_eLRef[index2]   -> Fill( anEvent->energyL_ref );
	  h1_eRRef[index2]   -> Fill( anEvent->energyR_ref );

	  // REF vs DUT quantities
	  h2_eL_eLRef[index2] -> Fill(anEvent->energyL_ref, anEvent->energyL );
	  h2_eR_eRRef[index2] -> Fill(anEvent->energyR_ref, anEvent->energyR );

	  // deltaT REF - DUTavg 
	  h1_deltaT_tLRef_tAve[index2] -> Fill( deltaTL );
	  h1_deltaT_tRRef_tAve[index2] -> Fill( deltaTR );
	  
	  // DEPENDENCIES
	  //    vs energy REF
	  p1_deltaT_tLRef_tAve_vs_eLRef[index2] -> Fill( anEvent->energyL_ref, deltaTL );
	  p1_deltaT_tRRef_tAve_vs_eRRef[index2] -> Fill( anEvent->energyR_ref, deltaTR );
	} // end loop over entries
    }      
  
  // - draw 2nd loop plots
  for(auto& it : h1_deltaT_tLRef_tAve)
    {      
      double index2 = it.first;      

      // plt TH2F
      SaveHisto2ToCanvas(outFile, h2_eL_eLRef[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_eR_eRRef[index2], plotDir);

      // fit gaus and plot energy
      f_eLRef[index2]    = FitAndSaveHisto(outFile, h1_eLRef[index2], plotDir,1.0,0.5, 1, kLandau);
      f_eRRef[index2]    = FitAndSaveHisto(outFile, h1_eRRef[index2], plotDir,1.0,0.5, 1, kLandau);

      // fit gaus and plot deltaT
      f = FitAndSaveHisto(outFile, h1_deltaT_tLRef_tAve[index2], plotDir);
      f = FitAndSaveHisto(outFile, h1_deltaT_tRRef_tAve[index2], plotDir);

      // fit pol 1 and plot
      //     vs energy REF
      f_LAve_eLRef[index2]    = FitAndSaveProfile(outFile, p1_deltaT_tLRef_tAve_vs_eLRef[index2], plotDir);
      f_RAve_eRRef[index2]    = FitAndSaveProfile(outFile, p1_deltaT_tRRef_tAve_vs_eRRef[index2], plotDir);

      int index1 = (int)index2 % 10000000;
      EntryData &e = dataMap[index2];     
      e.index = index2;
      e.mpvL = mpv["L"][index1];
      e.mpvR = mpv["R"][index1];
      e.mpvLR = mpv["L-R"][index1];	
      e.slopeL_ref  = f_LAve_eLRef[index2]->GetParameter(1);
      e.slopeR_ref  = f_RAve_eRRef[index2]->GetParameter(1);
      e.offsetL_ref = f_LAve_eLRef[index2]->GetParameter(0);
      e.offsetR_ref = f_RAve_eRRef[index2]->GetParameter(0);      
    }

  // =========================
  // 3rd loop
  //  - build corrected reference time
  //  - apply energy correctons to the reference module time
  // =========================

  // 0c. DUT energy
  std::map<double,TH1F*> h1_eL;
  std::map<double,TH1F*> h1_eR;  
  std::map<double,TF1*> f_eL;
  std::map<double,TF1*> f_eR;

  // 0d. DUT phase   
  std::map<double,TH2F*> h2_phaseL_vs_phaseAveRef;
  std::map<double,TH2F*> h2_phaseR_vs_phaseAveRef;

  // 2. DUT - REF eCor 
  // -    DUT module infos - tAverage REF corrected for TW: time left/right DUT - time average REF corrected
  // 2a. non-corrected REF
  std::map<double,TH1F*> h1_deltaT_tL_tAveRef;
  std::map<double,TH1F*> h1_deltaT_tR_tAveRef;

  std::map<double,TH2F*> h2_deltaT_tL_tAveRef_vs_eL;
  std::map<double,TH2F*> h2_deltaT_tL_tAveRef_vs_eAveRef;

  // 2b. corrected REF
  std::map<double,TH1F*> h1_deltaT_tL_tAveRefCor;
  std::map<double,TH1F*> h1_deltaT_tR_tAveRefCor;

  std::map<double,TH2F*> h2_deltaT_tL_tAveRefCor_vs_eL;
  std::map<double,TH2F*> h2_deltaT_tL_tAveRefCor_vs_eAveRef;
  
  // 2c. get DUT energy correction
  std::map<double,TProfile*> p1_deltaT_tL_tAveRefCor_vs_eL;
  std::map<double,TProfile*> p1_deltaT_tR_tAveRefCor_vs_eR;

  std::map<double,TF1*> f_LAve_eL;
  std::map<double,TF1*> f_RAve_eR;

  // 2d. tdiff objects
  std::map<double,TH1F*>     h1_deltaT_tL_tR;
  std::map<double,TH1F*>     h1_eRatio;
  std::map<double,TF1*>      f_eRatio;
  std::map<double,TF1*>      f_LR;
  std::map<double,TProfile*> p1_deltaT_tL_tR_vs_eRatio;
  std::map<double,TH1F*>     h1_deltaT_tL_tR_eRatioCor;

  std::map<double,TH1F*> h1_phaseAve;
  std::map<double,TH1F*> h1_phaseR;
  std::map<double,TH1F*> h1_phaseL;
  
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
	  if( !accept[index1][entry] ) continue;
	  int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;
	  double  index2( 10000000*energyBinAverage+index1 );
	  if (!f_LAve_eLRef[index2] || !f_RAve_eRRef[index2] ) {
	    accept[index1][entry] = false;
	    continue;
	  }
	  
	  if( h1_deltaT_tLRef_tAve_eLRefCor[index2] == NULL )
	    {
	      std::string labelLR_energyBin(Form("bar%02d_Vov%.02f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth,energyBinAverage));

	      // ---- REF cor time - DUTavg 
	      h1_deltaT_tLRef_tAve_eLRefCor[index2] = new TH1F(Form("h1_deltaT_tLRef_tAve_eLRefCor_%s",labelLR_energyBin.c_str()),";t_{L}^{REF, cor} - t_{avg}; Events",nbins,minT,maxT);
	      h1_deltaT_tRRef_tAve_eRRefCor[index2] = new TH1F(Form("h1_deltaT_tRRef_tAve_eRRefCor_%s",labelLR_energyBin.c_str()),";t_{R}^{REF, cor} - t_{avg}; Events",nbins,minT,maxT);

	      h2_deltaT_tLRef_tAve_eRefCor_vs_eLRef[index2] = new TH2F(Form("h2_deltaT_tLRef_tAve_eRefCor_vs_eLRef_%s",labelLR_energyBin.c_str()),";energy_{L}^{REF};t_{L}^{REF,cor} - t_{avg}",nbins2d, 0,1024., nbins2d,minT,maxT);
	      h2_deltaT_tLRef_tAve_eRefCor_vs_eAve[index2] = new TH2F(Form("h2_deltaT_tLRef_tAve_eRefCor_vs_eAve_%s",labelLR_energyBin.c_str()),";energy_{avg};t_{L}^{REF,cor} - t_{avg}",nbins2d, 0,1024., nbins2d,minT,maxT);
	      h2_deltaT_tRRef_tAve_eRefCor_vs_eAve[index2] = new TH2F(Form("h2_deltaT_tRRef_tAve_eRefCor_vs_eAve_%s",labelLR_energyBin.c_str()),";energy_{avg};t_{R}^{REF,cor} - t_{avg}",nbins2d, 0,1024., nbins2d,minT,maxT);
	      
	      // ---- DUT time - REFavg
	      // phase
	      h1_phaseL[index2]   = new TH1F(Form("h1_phaseL_%s",     labelLR_energyBin.c_str()), ";phase_{L};Events"  ,nbins, 100.,1000.);
	      h1_phaseR[index2]   = new TH1F(Form("h1_phaseR_%s",     labelLR_energyBin.c_str()), ";phase_{R};Events"  ,nbins, 100.,1000.);
	      h1_phaseAve[index2] = new TH1F(Form("h1_phaseAve_%s",labelLR_energyBin.c_str()),    ";phase_{avg};Events",nbins, 100.,1000.);
	      
	      h2_phaseL_vs_phaseAveRef[index2] = new TH2F(Form("h1_phaseL_vs_phaseAveRef_%s",labelLR_energyBin.c_str()),";phase_{avg}^{REF};phase_{L};Events",nbins2d, 100.,1000.,nbins2d, 100.,1000.);
	      h2_phaseR_vs_phaseAveRef[index2] = new TH2F(Form("h1_phaseR_vs_phaseAveRef_%s",labelLR_energyBin.c_str()),";phase_{avg}^{REF};phase_{R};Events",nbins2d, 100.,1000.,nbins2d, 100.,1000.);

	      // energy
	      h1_eL[index2] = new TH1F(Form("h1_eL_%s",labelLR_energyBin.c_str()),";energy_{L};Events",nbins, 0.,1024.);
	      h1_eR[index2] = new TH1F(Form("h1_eR_%s",labelLR_energyBin.c_str()),";energy_{R};Events",nbins, 0.,1024.);

	      // deltaT
	      h1_deltaT_tL_tAveRef[index2] = new TH1F(Form("h1_deltaT_tL_tAveRef_%s",labelLR_energyBin.c_str()),";t_{L} - t_{avg}^{REF}; Events",nbins, maxT,-minT);
	      h1_deltaT_tR_tAveRef[index2] = new TH1F(Form("h1_deltaT_tR_tAveRef_%s",labelLR_energyBin.c_str()),";t_{R} - t_{avg}^{REF}; Events",nbins, maxT,-minT);

	      h2_deltaT_tL_tAveRef_vs_eL[index2] = new TH2F(Form("h2_deltaT_tL_tAveRef_vs_eL_%s",labelLR_energyBin.c_str()), ";energy_{L};t_{L} - t_{avg}^{REF};Events", nbins2d, 0, 1024, nbins2d,maxT,-minT);
	      h2_deltaT_tL_tAveRef_vs_eAveRef[index2] = new TH2F(Form("h2_deltaT_tL_tAveRef_vs_eAveRef_%s",labelLR_energyBin.c_str()), ";energy_{avg}^{REF};t_{L} - t_{avg}^{REF};Events", nbins2d, 0, 1024, nbins2d,maxT,-minT);

	      // ---- DUT time - REFavg - eRefCorAvg
	      h1_deltaT_tL_tAveRefCor[index2] = new TH1F(Form("h1_deltaT_tL_tAveRefCor_%s",labelLR_energyBin.c_str()),";t_{L} - t_{avg}^{REF, cor}; Events",nbins, maxT,-minT);
	      h1_deltaT_tR_tAveRefCor[index2] = new TH1F(Form("h1_deltaT_tR_tAveRefCor_%s",labelLR_energyBin.c_str()),";t_{R} - t_{avg}^{REF, cor}; Events",nbins, maxT,-minT);

	      h2_deltaT_tL_tAveRefCor_vs_eL[index2] = new TH2F(Form("h2_deltaT_tL_tAveRefCor_vs_eL_%s",labelLR_energyBin.c_str()), ";energy_{L};t_{L} - t_{avg}^{REF,cor};Events", nbins2d, 0, 1024, nbins2d,maxT,-minT);
	      h2_deltaT_tL_tAveRefCor_vs_eAveRef[index2] = new TH2F(Form("h2_deltaT_tL_tAveRefCor_vs_eAveRef_%s",labelLR_energyBin.c_str()), ";energy_{avg}^{REF};t_{L} - t_{avg}^{REF,cor};Events", nbins2d, 0, 1024, nbins2d,maxT,-minT);

	      p1_deltaT_tL_tAveRefCor_vs_eL[index2] = new TProfile(Form("p1_deltaT_tL_tAveRefCor_vs_eL_%s",labelLR_energyBin.c_str()),";energy_{L};t_{L} - t_{avg}^{REF, cor}",nbinsProf,0,1024);
	      p1_deltaT_tR_tAveRefCor_vs_eR[index2] = new TProfile(Form("p1_deltaT_tR_tAveRefCor_vs_eR_%s",labelLR_energyBin.c_str()),";energy_{R};t_{R} - t_{avg}^{REF, cor}",nbinsProf,0,1024);

	      // ---- DUT time difference
	      h1_eRatio[index2] = new TH1F(Form("h1_eRatio_%s",labelLR_energyBin.c_str()),";energy_{L}/energy_{R};Events",nbins, minEnergyRatio, maxEnergyRatio);
	      h1_deltaT_tL_tR[index2] = new TH1F(Form("h1_deltaT_tL_tR_%s",labelLR_energyBin.c_str()),";t_{L} - t_{R}; Events",nbins, -1000, 1000);
	      p1_deltaT_tL_tR_vs_eRatio[index2] = new TProfile(Form("p1_deltaT_tL_tR_vs_eRatio_%s",labelLR_energyBin.c_str()),";energy_{L}/energy_{R};t_{L} - t_{R}", nbinsProf, minEnergyRatio, maxEnergyRatio);
	    }

	  // ---------- Define quantitities --------
	  double energyAve   = 0.5*( anEvent->energyL   + anEvent->energyR   );
	  long long tAve_ref = 0.5*( anEvent->timeL_ref + anEvent->timeR_ref );
	  double eLRefCor    = f_LAve_eLRef[index2] -> Eval( anEvent->energyL_ref ) - f_LAve_eLRef[index2] -> Eval( f_eLRef[index2]->GetParameter(1) );
	  double eRRefCor    = f_RAve_eRRef[index2] -> Eval( anEvent->energyR_ref ) - f_RAve_eRRef[index2] -> Eval( f_eRRef[index2]->GetParameter(1) );
	  
	  // --- REF time - DUTavg - eRefCor
	  double deltaTL_eRefCor   = anEvent->timeL_ref - ( 0.5*( anEvent->timeL + anEvent->timeR) ) - eLRefCor;
	  double deltaTR_eRefCor   = anEvent->timeR_ref - ( 0.5*( anEvent->timeL + anEvent->timeR) ) - eRRefCor;

	  // --- DUT time - REFavg
	  double deltaTL_AveRef = anEvent->timeL - tAve_ref;
	  double deltaTR_AveRef = anEvent->timeR - tAve_ref;

	  // --- DUT time - REFavg - eRefCorAvg
	  double deltaTL_AveRefCor = anEvent->timeL - (tAve_ref - 0.5*( eLRefCor + eRRefCor ));
	  double deltaTR_AveRefCor = anEvent->timeR - (tAve_ref - 0.5*( eLRefCor + eRRefCor ));

	  double energyRefAve = 0.5*( anEvent->energyL_ref + anEvent->energyR_ref );
	  double phaseRefAve = 0.5*( anEvent->t1fineL_ref + anEvent->t1fineR_ref );

	  // ---------- Fill histograms -------
	  // --- REF time - DUTavg - eRefCor
	  if (energyAve > 250 && energyAve<350) {	  
	    h1_deltaT_tLRef_tAve_eLRefCor[index2] -> Fill( deltaTL_eRefCor );
	    h1_deltaT_tRRef_tAve_eRRefCor[index2] -> Fill( deltaTR_eRefCor );

	    h2_deltaT_tLRef_tAve_eRefCor_vs_eLRef[index2] -> Fill( anEvent->energyL_ref, deltaTL_eRefCor );
	    h2_deltaT_tLRef_tAve_eRefCor_vs_eAve[index2] -> Fill( energyAve, deltaTL_eRefCor );
	    h2_deltaT_tRRef_tAve_eRefCor_vs_eAve[index2] -> Fill( energyAve, deltaTR_eRefCor );

	    h2_phaseL_vs_phaseAveRef[index2] -> Fill( phaseRefAve, anEvent->t1fineL);
	    h2_phaseR_vs_phaseAveRef[index2] -> Fill( phaseRefAve, anEvent->t1fineR);
	  }

	  // --- DUT energy (filled here not to account for energy cut in the loop above)
	  h1_eL[index2] -> Fill( anEvent->energyL );
	  h1_eR[index2] -> Fill( anEvent->energyR );
	  
	  h1_phaseL[index2]   -> Fill( anEvent->t1fineL );
	  h1_phaseR[index2]   -> Fill( anEvent->t1fineR );
	  h1_phaseAve[index2] -> Fill( 0.5*(anEvent->t1fineL + anEvent->t1fineR) );

	  // --- DUT time - REFavg
	  h1_deltaT_tL_tAveRef[index2] -> Fill( deltaTL_AveRef );
	  h1_deltaT_tR_tAveRef[index2] -> Fill( deltaTR_AveRef );

	  h2_deltaT_tL_tAveRef_vs_eL[index2]       -> Fill( anEvent->energyL,  anEvent->timeL - tAve_ref);
	  h2_deltaT_tL_tAveRef_vs_eAveRef[index2]  -> Fill( energyRefAve,      anEvent->timeL - tAve_ref);

	  // --- DUT time - REFavg - eRefCorAvg
	  h1_deltaT_tL_tAveRefCor[index2] -> Fill( deltaTL_AveRefCor );
	  h1_deltaT_tR_tAveRefCor[index2] -> Fill( deltaTR_AveRefCor );
	  
	  p1_deltaT_tL_tAveRefCor_vs_eL[index2] -> Fill( anEvent->energyL, deltaTL_AveRefCor );
	  p1_deltaT_tR_tAveRefCor_vs_eR[index2] -> Fill( anEvent->energyR, deltaTR_AveRefCor );

	  h2_deltaT_tL_tAveRefCor_vs_eL[index2]       -> Fill( anEvent->energyL,  deltaTL_AveRefCor);
	  h2_deltaT_tL_tAveRefCor_vs_eAveRef[index2]  -> Fill( energyRefAve,      deltaTL_AveRefCor);
	  
	  // --- DUT time L - DUT time R
	  if (anEvent->energyL/anEvent->energyR >minEnergyRatio && anEvent->energyL/anEvent->energyR < maxEnergyRatio)
	    {
	      h1_eRatio[index2] -> Fill(anEvent->energyL/anEvent->energyR);
	      h1_deltaT_tL_tR[index2] -> Fill(anEvent->timeL - anEvent->timeR);
	      p1_deltaT_tL_tR_vs_eRatio[index2] -> Fill(anEvent->energyL/anEvent->energyR, anEvent->timeL - anEvent->timeR);
	    }
	} // end loop over entries
    }      
  // - draw 3rd loop plots
  for(auto& it : h1_deltaT_tLRef_tAve_eLRefCor)
    {      
      double index2 = it.first;
      outFile->cd();
      f = FitAndSaveHisto(outFile, h1_deltaT_tLRef_tAve_eLRefCor[index2], plotDir);
      f = FitAndSaveHisto(outFile, h1_deltaT_tRRef_tAve_eRRefCor[index2], plotDir);

      SaveHisto2ToCanvas(outFile, h2_phaseL_vs_phaseAveRef[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_phaseR_vs_phaseAveRef[index2], plotDir);
      
      SaveHisto2ToCanvas(outFile, h2_deltaT_tLRef_tAve_eRefCor_vs_eLRef[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tLRef_tAve_eRefCor_vs_eAve[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tRRef_tAve_eRefCor_vs_eAve[index2], plotDir);

      f = FitAndSaveHisto(outFile, h1_deltaT_tL_tAveRef[index2], plotDir);
      f = FitAndSaveHisto(outFile, h1_deltaT_tR_tAveRef[index2], plotDir);

      f = FitAndSaveHisto(outFile, h1_deltaT_tL_tAveRefCor[index2], plotDir);
      f = FitAndSaveHisto(outFile, h1_deltaT_tR_tAveRefCor[index2], plotDir);

      f_eL[index2] = FitAndSaveHisto(outFile, h1_eL[index2], plotDir,1.0,0.5, 1, kLandau);
      f_eR[index2] = FitAndSaveHisto(outFile, h1_eR[index2], plotDir,1.0,0.5, 1, kLandau);
	
      f_LAve_eL[index2] = FitAndSaveProfile(outFile, p1_deltaT_tL_tAveRefCor_vs_eL[index2], plotDir);
      f_RAve_eR[index2]	= FitAndSaveProfile(outFile, p1_deltaT_tR_tAveRefCor_vs_eR[index2], plotDir);

      EntryData &e = dataMap[index2];     
      e.index = index2;
      e.slopeL  = f_LAve_eL[index2]->GetParameter(1);
      e.slopeR  = f_RAve_eR[index2]->GetParameter(1);
      e.offsetL = f_LAve_eL[index2]->GetParameter(0);
      e.offsetR = f_RAve_eR[index2]->GetParameter(0);

      // deltaT LR
      f_eRatio[index2] = FitAndSaveHisto(outFile, h1_eRatio[index2], plotDir);
      f = FitAndSaveHisto(outFile,h1_deltaT_tL_tR[index2], plotDir);
      f_LR[index2] = FitAndSaveProfile(outFile, p1_deltaT_tL_tR_vs_eRatio[index2], plotDir);

      SaveHistoToCanvas(outFile, h1_phaseL[index2], plotDir);
      SaveHistoToCanvas(outFile, h1_phaseR[index2], plotDir);
      SaveHistoToCanvas(outFile, h1_phaseAve[index2], plotDir);
    }
  

  // =========================
  // 4th loop
  //  - build corrected DUT time
  // =========================
  // 2d. corrected REF and DUT time 
  std::map<double,TH1F*> h1_deltaT_tL_tAveRefCor_eCor;
  std::map<double,TH1F*> h1_deltaT_tR_tAveRefCor_eCor;

  // 2e. check residual dependence on energy DUT and energy REF
  std::map<double,TH2F*> h2_deltaT_tL_tAveRefCor_eCor_vs_eL;
  std::map<double,TH2F*> h2_deltaT_tL_tAveRefCor_eCor_vs_eAveRef;

  // JOLLY. triangolo variables
  std::map<double,double> s_tLRef;
  std::map<double,double> s_tRRef;
  std::map<double,double> s_tLR;
  std::map<double,TProfile*> p1_deltaT_tL_tR_eRatioCor_vs_phaseMean;
  
  // 3. derive phase corrections 
  std::map<double,TProfile*> p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL;
  std::map<double,TProfile*> p1_deltaT_tL_tAveRefCor_eCor_vs_phaseRefAve;
  std::map<double,TProfile*> p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR;  

  for(auto mapIt : trees)
    {
      ModuleEventWithRefClass* anEvent = new ModuleEventWithRefClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);
      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
	{
	  if( entry%100000 == 0 ) {
	    std::cout << ">>> 4th loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
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
	      h1_deltaT_tL_tAveRefCor_eCor[index2] = new TH1F(Form("h1_deltaT_tL_tAveRefCor_eCor_%s",labelLR_energyBin.c_str()), ";t_{L}^{cor} - t_{avg}^{REF,cor}; Events", nbins, maxT, -minT);
	      h1_deltaT_tR_tAveRefCor_eCor[index2] = new TH1F(Form("h1_deltaT_tR_tAveRefCor_eCor_%s",labelLR_energyBin.c_str()), ";t_{R}^{cor} - t_{avg}^{REF,cor}; Events", nbins, maxT, -minT);

	      h2_deltaT_tL_tAveRefCor_eCor_vs_eL[index2] = new TH2F(Form("h2_deltaT_tL_tAveRefCor_eCor_vs_eL_%s",labelLR_energyBin.c_str()), ";energy_{L};t_{L}^{cor} - t_{avg}^{REF,cor};Events", nbins2d, 0, 1024, nbins2d,maxT,-minT);
	      h2_deltaT_tL_tAveRefCor_eCor_vs_eAveRef[index2] = new TH2F(Form("h2_deltaT_tL_tAveRefCor_eCor_vs_eAveRef_%s",labelLR_energyBin.c_str()), ";energy_{avg}^{REF};t_{L}^{cor} - t_{avg}^{REF,cor};Events", nbins2d, 0, 1024, nbins2d,maxT,-minT);
	      
	      p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2]    = new TProfile(Form("p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL_%s",labelLR_energyBin.c_str()),";phase_{L};t_{L}^{cor} - t_{avg}^{REF, cor}",nbinsProf,100,1000);
	      p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2]    = new TProfile(Form("p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR_%s",labelLR_energyBin.c_str()),";phase_{R};t_{R}^{cor} - t_{avg}^{REF, cor}",nbinsProf,100,1000);
	      p1_deltaT_tL_tAveRefCor_eCor_vs_phaseRefAve[index2] = new TProfile(Form("p1_deltaT_tL_tAveRefCor_eCor_vs_phaseRefAve_%s",labelLR_energyBin.c_str()), ";phase_{avg}^{REF};t_{L}^{cor} - t_{avg}^{REF, cor}",nbinsProf,100,1000);	      

	      h1_deltaT_tL_tR_eRatioCor[index2] = new TH1F(Form("h1_deltaT_tL_tR_eRatioCor_%s",labelLR_energyBin.c_str()),";(t_{L} - t_{R})^{cor}; Events",nbins, -1000, 1000);
	      p1_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2] = new TProfile(Form("p1_deltaT_tL_tR_eRatioCor_vs_phaseMean_%s",labelLR_energyBin.c_str()), ";phase_{avg};t_{L}^{cor} - t_{R}^{cor}",nbinsProf, 100, 1000);
	    }
	  long long tAve_ref  = 0.5*( anEvent->timeL_ref   + anEvent->timeR_ref   );
	  double energyRefAve = 0.5*( anEvent->energyL_ref + anEvent->energyR_ref );
	  double phaseRefAve = 0.5*( anEvent->t1fineL_ref + anEvent->t1fineR_ref );	  

	  // --- correction REF vs energyRef
	  float eLRefCor   = f_LAve_eLRef[index2] -> Eval( anEvent->energyL_ref ) - f_LAve_eLRef[index2] -> Eval( f_eLRef[index2]->GetParameter(1) );
	  float eRRefCor   = f_RAve_eRRef[index2] -> Eval( anEvent->energyR_ref ) - f_RAve_eRRef[index2] -> Eval( f_eRRef[index2]->GetParameter(1) );
	  
	  // --- correction DUT vs energy
	  double eLCor = f_LAve_eL[index2] -> Eval(anEvent->energyL) - f_LAve_eL[index2] -> Eval(f_eL[index2]->GetParameter(1));
	  double eRCor = f_RAve_eR[index2] -> Eval(anEvent->energyR) - f_RAve_eR[index2] -> Eval(f_eR[index2]->GetParameter(1));	  

	  // --- DUT time - REFavg - eRefCor - eCor
	  double deltaTL_eLCor = anEvent->timeL - (tAve_ref - 0.5*(eLRefCor+eRRefCor)) - eLCor;
	  double deltaTR_eRCor = anEvent->timeR - (tAve_ref - 0.5*(eLRefCor+eRRefCor)) - eRCor;

	  h1_deltaT_tL_tAveRefCor_eCor[index2] -> Fill( deltaTL_eLCor );
	  h1_deltaT_tR_tAveRefCor_eCor[index2] -> Fill( deltaTR_eRCor );

	  h2_deltaT_tL_tAveRefCor_eCor_vs_eL[index2]      -> Fill( anEvent->energyL ,  deltaTL_eLCor );
	  h2_deltaT_tL_tAveRefCor_eCor_vs_eAveRef[index2] -> Fill( energyRefAve ,      deltaTL_eLCor );

	  p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2] -> Fill( anEvent->t1fineL, deltaTL_eLCor );
	  p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2] -> Fill( anEvent->t1fineR, deltaTR_eRCor );
	  p1_deltaT_tL_tAveRefCor_eCor_vs_phaseRefAve[index2] -> Fill( phaseRefAve, deltaTL_eLCor );
	  
	  // --- DUT time difference corrected
	  double eRatioCor = f_LR[index2] -> Eval( anEvent->energyL/anEvent->energyR ) - f_LR[index2] -> Eval( f_eRatio[index2]->GetParameter(1) );
	  double deltaTLR_eRatioCor = anEvent->timeL - anEvent->timeR - eRatioCor;
	  h1_deltaT_tL_tR_eRatioCor[index2]      -> Fill(deltaTLR_eRatioCor);
	  p1_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2] -> Fill( 0.5*(anEvent->t1fineL + anEvent->t1fineR), deltaTLR_eRatioCor);
	}
    }
  for(auto& it : h1_deltaT_tL_tAveRefCor_eCor)
    {
      double index2 = it.first;
      outFile->cd();
      // triangolo
      f = FitAndSaveHisto(outFile, h1_deltaT_tL_tAveRefCor_eCor[index2], plotDir);
      f = FitAndSaveHisto(outFile, h1_deltaT_tR_tAveRefCor_eCor[index2], plotDir);
      f = FitAndSaveHisto(outFile, h1_deltaT_tL_tR_eRatioCor[index2], plotDir);
      f = FitAndSaveProfile(outFile, p1_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2] , plotDir);
      
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRef_vs_eL[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRef_vs_eAveRef[index2], plotDir);      

      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRefCor_vs_eL[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRefCor_vs_eAveRef[index2], plotDir);      

      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRefCor_eCor_vs_eL[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRefCor_eCor_vs_eAveRef[index2], plotDir);

      f = FitAndSaveProfile(outFile, p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2], plotDir);
      f = FitAndSaveProfile(outFile, p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2], plotDir);
      f = FitAndSaveProfile(outFile, p1_deltaT_tL_tAveRefCor_eCor_vs_phaseRefAve[index2], plotDir);
    }

  // =================
  // 5th loop
  // =================

  // 3b. apply phase corrections
  std::map<double,TH1F*> h1_deltaT_tL_tAveRefCor_eCor_phaseCor;
  std::map<double,TH1F*> h1_deltaT_tR_tAveRefCor_eCor_phaseCor;

  std::map<double,TH1F*> h1_deltaT_tAve_tAveRefCor_eCor_phaseCor;

  // 3c. check residual dependence on phase
  std::map<double,TH2F*> h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseL;
  std::map<double,TH2F*> h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseAveRef;

  // JOLLY: triangolo
  std::map<double,TH1F*> h1_deltaT_tL_tR_eRatioCor_phaseMeanCor;
  
  minT = 0;
  maxT = 3000;
  for(auto mapIt : trees)
    {
      ModuleEventWithRefClass* anEvent = new ModuleEventWithRefClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);
      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
	{
	  if( entry%100000 == 0 ) {
	    std::cout << ">>> 5th loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
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
	      h1_deltaT_tL_tAveRefCor_eCor_phaseCor[index2] = new TH1F(Form("h1_deltaT_tL_tAveRefCor_eCor_phaseCor_%s",labelLR_energyBin.c_str()), ";t_{L}^{cor, ph} - t_{avg}^{REF,cor}; Events", nbins, minT, maxT);
	      h1_deltaT_tR_tAveRefCor_eCor_phaseCor[index2] = new TH1F(Form("h1_deltaT_tR_tAveRefCor_eCor_phaseCor_%s",labelLR_energyBin.c_str()), ";t_{R}^{cor, ph} - t_{avg}^{REF,cor}; Events", nbins, minT, maxT);
	      h1_deltaT_tAve_tAveRefCor_eCor_phaseCor[index2] = new TH1F(Form("h1_deltaT_tAve_tAveRefCor_eCor_phaseCor_%s",labelLR_energyBin.c_str()), ";t_{avg}^{cor, ph} - t_{avg}^{REF,cor}; Events", nbins, minT, maxT);

	      h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseAveRef[index2] = new TH2F(Form("h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseAveRef_%s",labelLR_energyBin.c_str()), ";t1fine_{avg}^{REF}; t_{L}^{cor, ph} - t_{avg}^{REF,cor}; Events", nbins2d, 100,1000, nbins2d,minT,maxT);
	      h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseL[index2] = new TH2F(Form("h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseL_%s",labelLR_energyBin.c_str()), ";t1fine_{L}; t_{L}^{cor, ph} - t_{avg}^{REF,cor}; Events", nbins2d, 100,1000, nbins2d, minT, maxT);


	      h1_deltaT_tL_tR_eRatioCor_phaseMeanCor[index2] = new TH1F(Form("h1_deltaT_tL_tR_eRatioCor_phaseMeanCor_%s",labelLR_energyBin.c_str()), ";t_{L}^{cor, ph} - t_{R}^{cor, ph}; Events", nbins, -1000,1000);
	    }
	  long long tAve_ref  = 0.5*( anEvent->timeL_ref   + anEvent->timeR_ref   );
	  double phaseRefAve = 0.5*( anEvent->t1fineL_ref + anEvent->t1fineR_ref );
	  
	  // --- correction REF vs energyRef
	  float eLRefCor   = f_LAve_eLRef[index2] -> Eval( anEvent->energyL_ref ) - f_LAve_eLRef[index2] -> Eval( f_eLRef[index2]->GetParameter(1) );
	  float eRRefCor   = f_RAve_eRRef[index2] -> Eval( anEvent->energyR_ref ) - f_RAve_eRRef[index2] -> Eval( f_eRRef[index2]->GetParameter(1) );
	  
	  // --- correction DUT vs energy
	  double eLCor = f_LAve_eL[index2] -> Eval(anEvent->energyL) - f_LAve_eL[index2] -> Eval(f_eL[index2]->GetParameter(1));
	  double eRCor = f_RAve_eR[index2] -> Eval(anEvent->energyR) - f_RAve_eR[index2] -> Eval(f_eR[index2]->GetParameter(1));	  

	  // --- correction DUT vs phase
	  int phaseBinL      = p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2]->FindBin( anEvent->t1fineL );
	  int phaseBinL_offs =  p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2]->FindBin( h1_phaseL[index2]->GetMean());
	  double phaseLCor   = p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2]->GetBinContent(phaseBinL) - p1_deltaT_tL_tAveRefCor_eCor_vs_phaseL[index2]->GetBinContent(phaseBinL_offs);
	  
	  int phaseBinR      = p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2]->FindBin( anEvent->t1fineR );
	  int phaseBinR_offs = p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2]->FindBin(  h1_phaseR[index2]->GetMean() );
	  double phaseRCor   = p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2]->GetBinContent(phaseBinR) - p1_deltaT_tR_tAveRefCor_eCor_vs_phaseR[index2]->GetBinContent(phaseBinR_offs);
	  
	  // --- DUT time - REFavg - eRefCor - eCor - phaseCor
	  double deltaTL_eLCor_phaseLCor = anEvent->timeL - (tAve_ref - 0.5*(eLRefCor+eRRefCor)) - eLCor - phaseLCor;
	  double deltaTR_eRCor_phaseRCor = anEvent->timeR - (tAve_ref - 0.5*(eLRefCor+eRRefCor)) - eRCor - phaseRCor;

	  // --- (DUT avg time - DUT avg eCor - DUT avg phaseCor) - (REF avg time - REF avg eCor)
	  double deltaT_tAve_tAveRefCor_eCor_phaseCor = (0.5* (anEvent->timeL+anEvent->timeR) - 0.5*(eLCor+eRCor)  - 0.5*(phaseLCor+phaseRCor)) - (tAve_ref - 0.5*(eLRefCor+eRRefCor));

	  h1_deltaT_tL_tAveRefCor_eCor_phaseCor[index2] -> Fill( deltaTL_eLCor_phaseLCor );
	  h1_deltaT_tR_tAveRefCor_eCor_phaseCor[index2] -> Fill( deltaTR_eRCor_phaseRCor );
	  
	  h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseAveRef[index2] -> Fill( phaseRefAve ,    deltaTL_eLCor_phaseLCor);
	  h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseL[index2] -> Fill( anEvent->t1fineL,    deltaTL_eLCor_phaseLCor);

	  h1_deltaT_tAve_tAveRefCor_eCor_phaseCor[index2] -> Fill( deltaT_tAve_tAveRefCor_eCor_phaseCor );
	  
	  // triangolo deltaT LR
	  // --- correction DUT vs phase
	  double phaseMean    = 0.5*(anEvent->t1fineL + anEvent->t1fineR );
	  int phaseBinLR      =  p1_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2]->FindBin( phaseMean );
	  int phaseBinLR_offs = p1_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2]->FindBin( h1_phaseAve[index2]->GetMean() );
	  double phaseMeanCor   = p1_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2]->GetBinContent(phaseBinLR) - p1_deltaT_tL_tR_eRatioCor_vs_phaseMean[index2]->GetBinContent(phaseBinLR_offs);

	  double eRatioCor = f_LR[index2] -> Eval( anEvent->energyL/anEvent->energyR ) - f_LR[index2] -> Eval( f_eRatio[index2]->GetParameter(1) );
	  double deltaTLR_eRatioCor_phaseMeanCor = anEvent->timeL - anEvent->timeR - eRatioCor - phaseMeanCor;

	  h1_deltaT_tL_tR_eRatioCor_phaseMeanCor[index2]  -> Fill( deltaTLR_eRatioCor_phaseMeanCor);	  
	}
    }
  for(auto& it : h1_deltaT_tL_tAveRefCor_eCor_phaseCor)
    {
      double index2 = it.first;
      outFile->cd();
      
      f = FitAndSaveHisto(outFile, h1_deltaT_tL_tAveRefCor_eCor_phaseCor[index2], plotDir);
      s_tLRef[index2] = f->GetParameter(2);
      
      f = FitAndSaveHisto(outFile, h1_deltaT_tR_tAveRefCor_eCor_phaseCor[index2], plotDir);
      s_tRRef[index2] =	f->GetParameter(2);
      
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseAveRef[index2], plotDir);
      SaveHisto2ToCanvas(outFile, h2_deltaT_tL_tAveRefCor_eCor_phaseCor_vs_phaseL[index2], plotDir);

      f = FitAndSaveHisto(outFile, h1_deltaT_tL_tR_eRatioCor_phaseMeanCor[index2], plotDir);
      s_tLR[index2] = f->GetParameter(2);

      f = FitAndSaveHisto(outFile, h1_deltaT_tAve_tAveRefCor_eCor_phaseCor[index2], plotDir);
      
      EntryData &e  = dataMap[index2];
      e.index       = index2;
      e.sigmaLRef   = s_tLRef[index2];
      e.sigmaLR     = s_tLR[index2];
      e.sigmaRRef   = s_tRRef[index2];
      e.sigmaDutRef = f->GetParameter(2);
      SigmaResult res = ilTriangoloNo(e.sigmaLRef, e.sigmaLR, e.sigmaRRef);
      e.sigmaRef    = res.sRef;
      e.sigmaL      = res.sL;
      e.sigmaR      = res.sR;
      std::cout <<"Triangolando " <<std::endl;
      std::cout <<"\t " <<e.sigmaLRef <<" \t "	<<e.sigmaLR <<" \t "	<<e.sigmaRRef <<" \t "    <<e.sigmaRef <<" \t "    <<e.sigmaL <<" \t "    <<e.sigmaR <<std::endl;  
    }

  std::cout <<"[INFO] Save output csv" <<std::endl;
  for (const auto &pair : dataMap) {
    const EntryData &e = pair.second;
    out << e.index << ","
        << e.mpvL << ","
	<< e.mpvR << ","
	<< e.mpvLR << ","
        << e.slopeL_ref << ","
        << e.slopeR_ref << ","
        << e.offsetL_ref << ","
        << e.offsetR_ref << ","
        << e.slopeL << ","
        << e.slopeR << ","
        << e.offsetL << ","
        << e.offsetR << ","
        << e.sigmaLRef << ","
        << e.sigmaLR << ","
        << e.sigmaRRef << ","
        << e.sigmaRef << ","
        << e.sigmaL << ","
        << e.sigmaR << ","
        << e.sigmaDutRef << "\n";
  }
  out.close();  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
