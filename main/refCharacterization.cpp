#include "interface/AnalysisUtils.h"
#include "interface/FitUtils.h"
#include "interface/SetTDRStyle.h"
#include "interface/referenceCharacterization_helper.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <filesystem>

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TNamed.h"

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
    std::cout << ">>> refCharacterization::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }

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
  // - options defining which amp walk functions to use: either each DUT bar computing its own for the same REF bar, or the refDUTbar for all DUT bars
  bool usePerBarAmpWalk = opts.GetOpt<bool>("Reference.useAmpWalkPerBar");
  int fixedAmpBar = opts.GetOpt<int>("Reference.ampWalkFixedDUTBar");
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
  if( minEnergiesFileName != "" ) {
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
      minE[std::make_pair(bar,ov)] = value; }
    for (unsigned int iBar = 0; iBar < 16; ++iBar) {
      for (unsigned int ii = 0; ii < Vov.size(); ++ii) {
	auto key = std::make_pair(iBar, Vov[ii]);
	if (minE.find(key) == minE.end()) {
	  minE[key] = map_energyMins[Vov[ii]]; } } } }
  else {
    for(unsigned int iBar = 0; iBar < 16; ++iBar)
      for(unsigned int ii = 0; ii < Vov.size(); ++ii)
	minE[std::make_pair(iBar, Vov[ii])] = map_energyMins[Vov[ii]]; }
  
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
  std::string outTreeFileName = opts.GetOpt<std::string>("Output.refFileName");
  TFile* outFile = TFile::Open(outTreeFileName.c_str(), "RECREATE");
  std::string step2FileName= opts.GetOpt<std::string>("Output.outFileNameStep2");
  TFile* step2File = TFile::Open(step2FileName.c_str(), "RECREATE");
  
  // - define histograms and TProfiles
  // 0b. reference energy
  std::map<double,TH1F*> h1_eLRef;
  std::map<double,TH1F*> h1_eRRef;  
  std::map<double,TF1*> f_eLRef;
  std::map<double,TF1*> f_eRRef;
  std::map<double,TProfile*> p1_deltaT_tLRef_tAve_vs_eLRef;
  std::map<double,TProfile*> p1_deltaT_tRRef_tAve_vs_eRRef;
  std::map<double,TF1*> f_LAve_eLRef;
  std::map<double,TF1*> f_RAve_eRRef;

  // - other maps
  std::map<std::string, std::map<int, std::vector<float>*> > ranges; //ranges[LRlabel][index]
  std::map<std::string, std::map<int, std::vector<float>*> > narrow_ranges; //ranges[LRlabel][index]
  std::map<std::string, std::map<int, std::map<int,float> > > energyBin; // energyBin[LRlabel][index]
  std::map<int,TF1*>  f_landau; 

  // - energy ranges (window around the MPV) are stored in a dedicated tree
  outFile->cd();
  TTree* rangeTree = new TTree("ranges","ranges");
  int tree_bar;
  int tree_vth;
  float tree_Vov;
  char tree_LR[16];
  float tree_rmin;
  float tree_rmax;
  float tree_rmin_narrow;
  float tree_rmax_narrow;
  rangeTree->Branch("bar",&tree_bar,"bar/I");
  rangeTree->Branch("vth",&tree_vth,"vth/I");
  rangeTree->Branch("Vov",&tree_Vov,"Vov/F");
  rangeTree->Branch("LR",tree_LR,"LR/C");
  rangeTree->Branch("rmin",&tree_rmin,"rmin/F");
  rangeTree->Branch("rmax",&tree_rmax,"rmax/F");
  rangeTree->Branch("rmin_narrow",&tree_rmin_narrow,"rmin_narrow/F");
  rangeTree->Branch("rmax_narrow",&tree_rmax_narrow,"rmax_narrow/F");
  
  // - vars
  TH1F* histo;
  int nbins = 400;
  int nbinsProf = 50;
  
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
      std::string VovLabel(Form("Vov%.2f",Vov));
      std::string thLabel(Form("th%02.0f",vth));
      std::cout << ">>> 1st loop: Vov " <<VovLabel <<" th " <<thLabel << std::endl;
      histo = (TH1F*)( inFile->Get(Form("h1_energy_external_barL-R_Vov%.2f_th%02.0f", Vov, vth)));
      if (!histo) {
	std::cerr << "[ERROR] missing external histogram for Vov=" << Vov << " th=" << vth << std::endl;
	continue;
      }
      // -- selecting external bar core energy, loop over DUT bars
      for(int iBar = 0; iBar < 16; ++iBar) {
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
	      narrow_ranges[LRLabel][index] -> push_back( f_landau[index]->GetParameter(1) - 1.0 * std::abs(f_landau[index]->GetParameter(2)));
	    }
	    else {
	      narrow_ranges[LRLabel][index] -> push_back( minE[std::make_pair(iBar, Vov)] );
	      ranges[LRLabel][index] -> push_back( minE[std::make_pair(iBar, Vov)] ); }
	    
	    // ----- set the energy maximum value to maximum ADC 
	    ranges[LRLabel][index] -> push_back( 940 );
	    narrow_ranges[LRLabel][index] -> push_back(f_landau[index]->GetParameter(1) + 1.0 * std::abs(f_landau[index]->GetParameter(2)));
	    // ----- store energy mean of each bin in the energyBin map
	    GetEnergyBins(histo, ranges[LRLabel][index], energyBin[LRLabel][index]);

	    // ----- fill ranges tree
	    tree_bar = iBar;
	    tree_vth = vth;
	    tree_Vov = Vov;
	    snprintf(tree_LR, sizeof(tree_LR), "%s", LRLabel.c_str());
	    tree_rmin = ranges[LRLabel][index]->at(0);
	    tree_rmax = ranges[LRLabel][index]->at(1);	    
	    tree_rmin_narrow = narrow_ranges[LRLabel][index]->at(0);
	    tree_rmax_narrow = narrow_ranges[LRLabel][index]->at(1);
	    rangeTree->Fill();
	  }// end MIP (TB)	  
	  histo->GetXaxis()->SetRangeUser(0,1024);
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
	  	  
	  // ========== SELECTIONS ===========
	  double deltaTL        = anEvent->timeL_ref - (0.5*(anEvent->timeL     + anEvent->timeR));
	  double deltaTR        = anEvent->timeR_ref - (0.5*(anEvent->timeL     + anEvent->timeR));
	  double deltaTL_AveRef = anEvent->timeL     - (0.5*(anEvent->timeL_ref + anEvent->timeR_ref));
	  double deltaTR_AveRef = anEvent->timeR     - (0.5*(anEvent->timeL_ref + anEvent->timeR_ref));	  
	  bool pass = PassSelection(anEvent, deltaTL, deltaTR, deltaTL_AveRef, deltaTR_AveRef, ranges, index1, anEvent->energyL, anEvent->energyR);
	  if (pass==false) continue;	   
	  accept[index1][entry] = true;
	  int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

	  // --- for the first loop we require the DUT events to have kinda fixed energy so that tDUT dependency on amplitude is reduced
	  float energyAve = 0.5*(anEvent->energyL + anEvent->energyR);
	  if (energyAve < narrow_ranges["L-R"][index1]->at(0) || energyAve> narrow_ranges["L-R"][index1]->at(1)) continue;	  
	  // =================================
	  double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );
	  if( p1_deltaT_tLRef_tAve_vs_eLRef[index2] == NULL )
	    {
	      std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth,energyBinAverage));

	      h1_eLRef[index2] = new TH1F(Form("h1_eLRef_%s",labelLR_energyBin.c_str()),";energy_{L}^{REF};Events",nbins, 0.,1024.);
	      h1_eRRef[index2] = new TH1F(Form("h1_eRRef_%s",labelLR_energyBin.c_str()),";energy_{R}^{REF};Events",nbins, 0.,1024.);

	      p1_deltaT_tLRef_tAve_vs_eLRef[index2] = new TProfile(Form("p1_deltaT_tLRef_tAve_vs_eLRef_%s",labelLR_energyBin.c_str()),";energy_{L}^{REF};t_{L}^{REF} - t_{avg}",nbinsProf,0,1024);
	      p1_deltaT_tRRef_tAve_vs_eRRef[index2] = new TProfile(Form("p1_deltaT_tRRef_tAve_vs_eRRef_%s",labelLR_energyBin.c_str()),";energy_{R}^{REF};t_{R}^{REF} - t_{avg}",nbinsProf,0,1024);
	    }
	  h1_eLRef[index2]   -> Fill( anEvent->energyL_ref );
	  h1_eRRef[index2]   -> Fill( anEvent->energyR_ref );
	  
	  p1_deltaT_tLRef_tAve_vs_eLRef[index2] -> Fill( anEvent->energyL_ref, deltaTL );
	  p1_deltaT_tRRef_tAve_vs_eRRef[index2] -> Fill( anEvent->energyR_ref, deltaTR );
	} // end loop over entries
    }      
  
  // - draw 2nd loop plots
  for(auto& it : p1_deltaT_tLRef_tAve_vs_eLRef)
    {      
      double index2 = it.first;
      outFile->cd();
      f_eLRef[index2]    = FitAndSaveHisto(outFile, h1_eLRef[index2], plotDir,1.0,0.5, 0);
      f_eRRef[index2]    = FitAndSaveHisto(outFile, h1_eRRef[index2], plotDir,1.0,0.5, 0);

      step2File->cd();
      f_LAve_eLRef[index2]    = FitAndSaveProfile(step2File, p1_deltaT_tLRef_tAve_vs_eLRef[index2], plotDir);
      f_RAve_eRRef[index2]    = FitAndSaveProfile(step2File, p1_deltaT_tRRef_tAve_vs_eRRef[index2], plotDir);
    }
  
  // ====================================
  // - 3rd LOOP -
  //   correct ref time and store to tree
  // ====================================
  for(auto mapIt : trees)
    {
      TTree* tree = mapIt.second;      
      ModuleEventWithRefClass* anEvent = nullptr;
      tree->SetBranchAddress("event",&anEvent);
      outFile->cd();
      TTree* outTree = tree->CloneTree(0);
      int nEntries = tree->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
	{
	  tree->GetEntry(entry);
	  if( entry%100000 == 0 ) {
	    std::cout << ">>> 3rd loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
	  }
	  if (!barSet.count(anEvent->barID)) continue;
	  int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );
	  if( !accept[index1][entry] ) continue;
	  int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;
	  int barForAmpWalk = usePerBarAmpWalk ? anEvent->barID : fixedAmpBar;
	  double index2 = 10000000*energyBinAverage + 10000*int(anEvent->Vov*100.) + 100*anEvent->vth + barForAmpWalk;
	  if (!f_LAve_eLRef[index2]) {
	    std::cout << "NULL f_LAve at index2=" << index2 << std::endl;
	    continue;
	  }
	  if (!f_LAve_eLRef[index2] || !f_eLRef[index2]) continue;
	  double eLRefCor    = f_LAve_eLRef[index2] -> Eval( anEvent->energyL_ref ) - f_LAve_eLRef[index2] -> Eval( f_eLRef[index2]->GetParameter(1) );
	  double eRRefCor    = f_RAve_eRRef[index2] -> Eval( anEvent->energyR_ref ) - f_RAve_eRRef[index2] -> Eval( f_eRRef[index2]->GetParameter(1) );
	  anEvent->timeL_ref_cor = anEvent->timeL_ref - eLRefCor;
	  anEvent->timeR_ref_cor = anEvent->timeR_ref - eRRefCor;
	  outTree->Fill();
	}
      outFile->cd();
      outTree->Write();
      tree->ResetBranchAddresses();
    }
  outFile->cd();
  rangeTree->Write();
  outFile->Close();
  step2File->Close();
  std::cout << "============================================"  << std::endl;
  std::cout << "============================================"  << std::endl;
}
