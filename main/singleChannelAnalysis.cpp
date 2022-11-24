#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
 
#include "interface/SetTDRStyle.h"  

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TLatex.h"
#include "TMath.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TApplication.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>


using namespace std;

void ilTriangoloNo(float s12, float s23, float s13, float err_s12, float err_s23, float err_s13, std::vector<float> &result, std::vector<float> &err) {
  
  float s1 = sqrt( 0.5 * ( s12*s12 + s13*s13 - s23*s23) );
  float s2 = sqrt( 0.5 * ( s12*s12 + s23*s23 - s13*s13) );
  float s3 = sqrt( 0.5 * ( s13*s13 + s23*s23 - s12*s12) );

  float err_s1 = 1./2/s1 * sqrt( pow(s12*err_s12,2) + pow(s13*err_s13,2) + pow(s23*err_s23,2) );
  float err_s2 = 1./2/s2 * sqrt( pow(s12*err_s12,2) + pow(s23*err_s23,2) + pow(s13*err_s13,2) );
  float err_s3 = 1./2/s3 * sqrt( pow(s13*err_s13,2) + pow(s23*err_s23,2) + pow(s12*err_s12,2) );

  result.push_back(s1);
  result.push_back(s2);
  result.push_back(s3);

  err.push_back(err_s1);
  err.push_back(err_s2);
  err.push_back(err_s3);
}



void getTimeResolution(TH1F *histo, TF1 *& fitFunc){

  float fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
  float fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;

  fitFunc -> SetParameters(histo->GetMaximum(), histo->GetBinCenter(histo->GetMaximumBin()), histo->GetRMS());
  fitFunc -> SetRange(fitXMin, fitXMax);
  histo -> Fit(fitFunc,"QRNL");
  fitFunc -> SetRange(fitFunc->GetParameter(1)-1.0*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+1.0*fitFunc->GetParameter(2));
  histo -> Fit(fitFunc,"QRNL");
  fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
  histo -> Fit(fitFunc,"QRSL+");
  histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
  histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-7.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+7.*fitFunc->GetParameter(2));
 }


// ================================================
int main(int argc, char** argv){

  setTDRStyle();
  gErrorIgnoreLevel = kError;  

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
  
  std::string runs = opts.GetOpt<std::string>("Input.runs");
  int chRef1 = opts.GetOpt<int>("Input.chRef1");
  int chRef2 = opts.GetOpt<int>("Input.chRef2");
  int useTimeAverage = opts.GetOpt<int>("Input.useTimeAverage");
  std::string source = opts.GetOpt<std::string>("Input.source");
 
  std::cout<< "Reference channels :  " <<  chRef1 << "  " << chRef2 <<std::endl;
  if ( useTimeAverage ) std::cout << "Using tAverage of the front module as reference"<<std::endl; 
  
  std::vector<unsigned int> channelMapping = opts.GetOpt<std::vector<unsigned int> >("Channels.channelMapping");
  int array = opts.GetOpt<int>("Channels.array"); 

  std::vector<std::string> labelLR = {"L","R"};
  map<std::string, int> chID;  
  for(int iBar = 0; iBar < 16; ++iBar){
    for (auto label : labelLR ){
      if ( label == "L") chID[ Form("bar%02d%s", iBar, label.c_str()) ] = channelMapping[iBar*2+0]+array*64;
      if ( label == "R") chID[ Form("bar%02d%s", iBar, label.c_str()) ] = channelMapping[iBar*2+1]+array*64;
    }
  }

  
  for(int iBar = 0; iBar < 16; ++iBar){
    for (auto label : labelLR ){
      std::string chLabel = Form("bar%02d%s", iBar, label.c_str()); 
      std::cout << iBar << "  " << label << "  " << chID[chLabel] <<   std::endl;
    }
  }

  int maxActiveChannels =  opts.GetOpt<int>("Cuts.maxActiveChannels");
  float minEnergy = opts.GetOpt<int>("Cuts.minEnergy");
  float maxEnergy = 950;
  //int mystep2 = opts.GetOpt<int>("Cuts.step2"); 

  std::vector<float> Vovs;
  Vovs.push_back(1.50);
  Vovs.push_back(2.50);
  Vovs.push_back(3.50);
  Vovs.push_back(4.00);
  Vovs.push_back(5.00);

  std::vector<int> thresholds;
  for (int i = 0; i < 64; i++){
    thresholds.push_back(i);
  }
  
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
      //std::cout<< bar <<  "   " << ov << "  " << minE[std::make_pair(bar,ov)] <<std::endl;
    }
  }
  else{
    for(unsigned int iBar = 0; iBar < 16; ++iBar){
      for(unsigned int ii = 0; ii < Vovs.size(); ++ii){
        minE[std::make_pair(iBar, Vovs[ii])] = minEnergy;
      }
    }
  }
          
  // --- reading tree
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
	std::string inFileName = Form("/data/TOFHIR2/reco/run%04d_e.root",run);
	std::cout << ">>> Adding file " << inFileName << std::endl;
	data -> Add(inFileName.c_str());
      }
    }

  //--- define branches
  float step1, step2;
  int channelIdx[256];
  std::vector<float> *tot = 0;
  std::vector<float> *energy = 0;
  std::vector<long long> *time = 0;
  std::vector<unsigned short>* t1fine = 0;
  
  data -> SetBranchStatus("*",0);
  data -> SetBranchStatus("step1",  1); data -> SetBranchAddress("step1",  &step1);
  data -> SetBranchStatus("step2",  1); data -> SetBranchAddress("step2",  &step2);
  data -> SetBranchStatus("channelIdx",  1); data -> SetBranchAddress("channelIdx",  channelIdx);
  data -> SetBranchStatus("tot",    1); data -> SetBranchAddress("tot",       &tot);
  data -> SetBranchStatus("energy", 1); data -> SetBranchAddress("energy", &energy);
  data -> SetBranchStatus("time",   1); data -> SetBranchAddress("time",     &time);
  data -> SetBranchStatus("t1fine",   1); data -> SetBranchAddress("t1fine",     &t1fine);

  int nEntries = data->GetEntries();
  cout << "Number of entries = " << nEntries << endl;
  //int maxEntries = 400000;
  int maxEntries = nEntries;
  


  // -- book histograms 
  std::cout<<"Booking histograms..."<<std::endl;
  map<int, TH1F*> h_nActiveChannels0;
  map<int, TH1F*> h_nActiveChannels1;
  map<int, TH1F*> h_energy_chRef;
  map<int, TH1F*> h_energyRatio_LR_chRef;
  map<int, TH1F*> h_deltaT_LR_chRef;
  map<int, TH1F*> h_deltaT_LR_energyRatioCorr_chRef;
  map<int, TH1F*> h_deltaT_LR_energyRatioCorr_phaseCorr_chRef;
  map<int, TProfile*> p_deltaT_LR_vs_energyRatio_chRef;
  map<int, TProfile*> p_deltaT_LR_energyRatioCorr_vs_t1fine_chRef;
   
  map<std::string, map<int,TH1F*> >      h_energy;
  map<std::string, map<int,TH1F*> >      h_energyRatio;
  map<std::string, map<int,TH2F*> >      h2_energyRatio_vs_totRatio;
  map<std::string, map<int,TH1F*> >      h_deltaT;
  map<std::string, map<int,TH1F*> >      h_deltaT_energyRatioCorr;
  map<std::string, map<int,TH1F*> >      h_deltaT_energyRatioCorr_phaseCorr;
  map<std::string, map<int,TProfile*> >  p_deltaT_vs_energyRatio;
  map<std::string, map<int,TH2F*> >      h2_deltaT_vs_energyRatio;
  map<std::string, map<int,TProfile*> >  p_deltaT_energyRatioCorr_vs_t1fine;
  map<std::string, map<int,TH2F*> >      h2_deltaT_energyRatioCorr_vs_t1fine;
  map<std::string, map<int,TProfile*> >  p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef;
  map<std::string, map<int,TH2F*> >      h2_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef;

  map<int, map<int,TH1F*> >      h_energy_LR;
  map<int, map<int,TH1F*> >      h_energyRatio_LR;
  map<int, map<int,TH1F*> >      h_deltaT_LR;
  map<int, map<int,TH1F*> >      h_deltaT_LR_energyRatioCorr;
  map<int, map<int,TH1F*> >      h_deltaT_LR_energyRatioCorr_phaseCorr;
  map<int, map<int,TProfile*> >  p_deltaT_LR_vs_energyRatio;
  map<int, map<int,TProfile*> >  p_deltaT_LR_energyRatioCorr_vs_t1fine;


  for (auto & th : thresholds){
        
    h_nActiveChannels0[th] = new TH1F(Form("h_nActiveChannels0_th%02d",th),"h_nActiveChannels0",32,-0.5,31.5);
    h_nActiveChannels1[th] = new TH1F(Form("h_nActiveChannels1_th%02d",th),"h_nActiveChannels1",32,-0.5,31.5);
    h_energy_chRef[th]     = new TH1F(Form("h_energy_chRef_th%02d",th),"h_energy_chRef",512,0,1024);
    h_energyRatio_LR_chRef[th]= new TH1F(Form("h_energyRatio_LR_chRef_th%02d",th),"h_energyRatio_LR_chRef",100,0,3);
    h_deltaT_LR_chRef[th]  = new TH1F(Form("h_deltaT_LR_chRef_th%02d",th),"h_deltaT_LR_chRef",  4000, -24000, 24000); 
    h_deltaT_LR_energyRatioCorr_chRef[th]  = new TH1F(Form("h_deltaT_LR_energyRatioCorr_chRef_th%02d",th),"h_deltaT_LR_energyRatioCorr_chRef",  4000, -24000, 24000); 
    h_deltaT_LR_energyRatioCorr_phaseCorr_chRef[th]  = new TH1F(Form("h_deltaT_LR_energyRatioCorr_phaseCorr_chRef_th%02d",th),"h_deltaT_LR_energyRatioCorr_phaseCorr_chRef",  4000, -24000, 24000); 
    p_deltaT_LR_vs_energyRatio_chRef[th]  = new TProfile(Form("p_deltaT_LR_vs_energyRatio_chRef_th%02d",th),"p_deltaT_LR_vs_energyRatio_chRef",  180, 0, 3); 
    p_deltaT_LR_energyRatioCorr_vs_t1fine_chRef[th]  = new TProfile(Form("p_deltaT_LR_energyRatioCorr_vs_t1fine_chRef_th%02d",th),"p_deltaT_LR_energyRatioCorr_vs_t1fine_chRef", 50, 0, 1000); 
        
    for(int iBar = 0; iBar < 16; ++iBar){        
      
      // L-R
      h_energy_LR[iBar][th] = new TH1F(Form("h_energy_LR_bar%02d_th%02d", iBar, th), Form("h_energy_LR_bar%02d_th%02d", iBar, th), 512, 0, 1024);
      h_energyRatio_LR[iBar][th] = new TH1F(Form("h_energyRatio_LR_bar%02d_th%02d", iBar, th), Form("h_energyRatio_LR_bar%02d_th%02d", iBar, th), 100, 0, 3);
      h_deltaT_LR[iBar][th] = new TH1F(Form("h_deltaT_LR_bar%02d_th%02d", iBar, th), Form("h_deltaT_LR_bar%02d_th%02d", iBar, th), 4000, -24000, 24000);
      h_deltaT_LR_energyRatioCorr[iBar][th] = new TH1F(Form("h_deltaT_LR_energyRatioCorr_bar%02d_th%02d", iBar, th), Form("h_deltaT_LR_energyRatioCorr_bar%02d_th%02d", iBar, th), 4000, -24000, 24000);
      h_deltaT_LR_energyRatioCorr_phaseCorr[iBar][th] = new TH1F(Form("h_deltaT_LR_energyRatioCorr_phaseCorr_bar%02d_th%02d", iBar, th), Form("h_deltaT_LR_energyRatioCorr_phaseCorr_bar%02d_th%02d", iBar, th), 4000, -24000, 24000);
      p_deltaT_LR_vs_energyRatio[iBar][th] = new TProfile(Form("p_deltaT_LR_vs_energyRatio_bar%02d_th%02d", iBar, th), Form("p_deltaT_LR_vs_energyRatio_bar%02d_th%02d", iBar, th), 60, 0, 3);
      p_deltaT_LR_energyRatioCorr_vs_t1fine[iBar][th] = new TProfile(Form("p_deltaT_LR_energyRatioCorr_vs_t1fine_bar%02d_th%02d", iBar, th), Form("p_deltaT_LR_energyRatioCorr_vs_t1fine_bar%02d_th%02d", iBar, th), 50, 0, 1000);

      for (auto label : labelLR ){
	std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
	
	h_energy[chLabel][th] = new TH1F(Form("h_energy_%s_th%02d", chLabel.c_str(), th) , Form("h_energy_%s_th%02d", chLabel.c_str(), th), 512, 0, 1024);
	h_energyRatio[chLabel][th] = new TH1F(Form("h_energyRatio_%s_th%02d", chLabel.c_str(), th) , Form("h_energyRatio_%s_th%02d", chLabel.c_str(), th), 100, 0, 3);
	h2_energyRatio_vs_totRatio[chLabel][th] = new TH2F(Form("h2_energyRatio_vs_totRatio_%s_th%02d", chLabel.c_str(), th) , Form("h_energyRatio_vs_totRatio_%s_th%02d", chLabel.c_str(), th), 100, 0, 3,100,0,3);
	h_deltaT[chLabel][th] = new TH1F(Form("h_deltaT_%s_th%02d", chLabel.c_str(), th) , Form("h_deltaT_%s_th%02d", chLabel.c_str(), th), 4000, -24000, 24000);
	h_deltaT_energyRatioCorr[chLabel][th] = new TH1F(Form("h_deltaT_energyRatioCorr_%s_th%02d", chLabel.c_str(), th) , Form("h_deltaT_energyRatioCorr_%s_th%02d", chLabel.c_str(), th), 4000, -24000, 24000);
	h_deltaT_energyRatioCorr_phaseCorr[chLabel][th] = new TH1F(Form("h_deltaT_energyRatioCorr_phaseCorr_%s_th%02d", chLabel.c_str(), th) , Form("h_deltaT_energyRatioCorr_phaseCorr_%s_th%02d", chLabel.c_str(), th), 4000, -24000, 24000);
	p_deltaT_energyRatioCorr_vs_t1fine[chLabel][th] = new TProfile(Form("p_deltaT_energyRatioCorr_vs_t1fine_%s_th%02d", chLabel.c_str(), th) , Form("p_deltaT_energyRatioCorr_vs_t1fine_%s_th%02d", chLabel.c_str(), th), 50, 0, 1000);
	h2_deltaT_energyRatioCorr_vs_t1fine[chLabel][th] = new TH2F(Form("h2_deltaT_energyRatioCorr_vs_t1fine_%s_th%02d", chLabel.c_str(), th) , Form("h2_deltaT_energyRatioCorr_vs_t1fine_%s_th%02d", chLabel.c_str(), th), 50, 0, 1000, 4000, -24000, 24000);
	p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel][th] = new TProfile(Form("p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef_%s_th%02d", chLabel.c_str(), th) , Form("p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef_%s_th%02d", chLabel.c_str(), th), 50, 0, 1000);
	h2_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel][th] = new TH2F(Form("h2_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef_%s_th%02d", chLabel.c_str(), th) , Form("h2_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef_%s_th%02d", chLabel.c_str(), th), 50, 0, 1000, 4000, -24000, 24000);
      }
    } 
  }

  map< int , map <int, map<int, bool> > > acceptEvent; //acceptEvent[entry][bar][thr]
  map< int, map<int, bool> >  acceptEvent_chRef;// acceptEvent_chRef[entry][thr] 

  map< int, map<int, int> > nActiveChannels0;
  map< int, map<int, int> > nActiveChannels1;

  float maxDeltaT = 20000;
  
  // -- first loop over events
  cout << "First loop over events to find the mip peak" <<endl;
  float Vov = 0;
  for (int entry = 0; entry < maxEntries; ++entry){
    
    data->GetEntry(entry);
    Vov = step1;

    if( entry%1000 == 0 ) std::cout << ">>> Reading entry " << entry << " / " << nEntries << "\r" << std::flush;

    //if (step2 != mystep2) continue;
    int th = int(step2/10000)-1;

    
    // -- count active channels in the two modules
    nActiveChannels0[entry][th] = 0;
    for (int ch = 0; ch < 32; ch++) {
      if ( channelIdx[ch] < 0 ) continue;
      if ( (*tot)[channelIdx[ch]]/1000 < -10. || (*tot)[channelIdx[ch]]/1000 > 50. ) continue;
      if ((*energy)[channelIdx[ch]] > 0)
	nActiveChannels0[entry][th]+=1;
    }    
    
    nActiveChannels1[entry][th] = 0;
    for (int ch = 64; ch < 96 ; ch++) {
      if ( channelIdx[ch] < 0 ) continue;
      if ( (*tot)[channelIdx[ch]]/1000 < -10. || (*tot)[channelIdx[ch]]/1000 > 50. ) continue;
      if ((*energy)[channelIdx[ch]] > 0)
	nActiveChannels1[entry][th]+=1;
    }
    
    h_nActiveChannels0[th] ->Fill(nActiveChannels0[entry][th]);
    h_nActiveChannels1[th] ->Fill(nActiveChannels1[entry][th]);


    if (nActiveChannels0[entry][th] >  maxActiveChannels) continue;
    if (nActiveChannels1[entry][th] >  maxActiveChannels) continue;

 

    //-- ref channel
    if ( channelIdx[chRef1] < 0 ) continue;
    if ( channelIdx[chRef2] < 0 ) continue;
    if ( (*tot)[channelIdx[chRef1]]/1000 < -10. || (*tot)[channelIdx[chRef1]]/1000 > 50. ) continue;      
    if ( (*tot)[channelIdx[chRef2]]/1000 < -10. || (*tot)[channelIdx[chRef2]]/1000 > 50. ) continue;      

    float energyRef = 0.5* ((*energy)[channelIdx[chRef1]]+(*energy)[channelIdx[chRef2]] );
    if ( !useTimeAverage) energyRef  = (*energy)[channelIdx[chRef1]];

    h_energy_chRef[th] -> Fill( energyRef );
    
    for(int iBar = 0; iBar < 16; ++iBar){        
      float energyR = 0;
      float energyL = 0;
      for (auto label : labelLR ){
	std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
	int ch = chID[chLabel]; 
	
	if ( channelIdx[ch] < 0 ) continue;
	if ( (*tot)[channelIdx[ch]]/1000 < -10. || (*tot)[channelIdx[ch]]/1000 > 50. ) continue;

	if (label == "R") energyR = (*energy)[channelIdx[ch]];
	if (label == "L") energyL = (*energy)[channelIdx[ch]];
      }

      if (energyL > 0 && energyR > 0 ) h_energy_LR[iBar][th] -> Fill( 0.5*(energyL+energyR) );

    }// end loop over bars

  }// -- end first loop over entries


  // -- find min energy
  map < int, float > energyMin_chRef;

  for (auto & th: thresholds){
    energyMin_chRef[th] = 0;
    
    if ( h_energy_chRef[th]->GetEntries() == 0) continue; 

    if (!source.compare("Laser")){
      TF1 *fitEnergy_chRef = new TF1(Form("fitGaus_chRef_th%02d", th) ,"gaus", 0, 1000);
      h_energy_chRef[th]->Fit(fitEnergy_chRef,"QR");
      energyMin_chRef[th] = fitEnergy_chRef->GetParameter(1) - 5.0*fitEnergy_chRef->GetParameter(2);
    }
    else{
      TF1 *fitEnergy_chRef = new TF1(Form("fitLandau_chRef_th%02d", th),"landau", 0, 1000);
      h_energy_chRef[th]->GetXaxis()->SetRangeUser(200,800);
      int maxbin = h_energy_chRef[th]->GetMaximumBin();
      float peak = h_energy_chRef[th]->GetBinCenter(maxbin);
      fitEnergy_chRef->SetRange(peak*0.8, peak*1.2);
      fitEnergy_chRef->SetParameter(1,peak);
      fitEnergy_chRef->SetParameter(2,0.1*peak);
      h_energy_chRef[th]->Fit(fitEnergy_chRef,"QR");
      energyMin_chRef[th] = fitEnergy_chRef->GetParameter(1)*0.8;
      h_energy_chRef[th]->GetXaxis()->SetRangeUser(0,1000);
    }
    std::cout << " th = " << th <<    " minEnergy chRef = " << energyMin_chRef[th] << std::endl;
  }


  map<std::string, map<int, float> > energyMin;
  map<std::string, map<int, TF1*> > fitLandau;
  map<int, map<int, float> > energyMinLR;
  map<int, map<int, TF1*> > fitEnergyLR;

  for(int iBar = 0; iBar < 16; ++iBar){        
    // -- fit mip peak single channels
    /*for (auto label : labelLR ){
      std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
      fitLandau[chLabel] = new TF1(Form("fitLandau_%s",chLabel.c_str()),"landau", 0, 1000);    
      h_energy[chLabel]->GetXaxis()->SetRangeUser(minE[std::make_pair(iBar, Vov)],800);
      if (chLabel == "bar13L") h_energy[chLabel]->GetXaxis()->SetRangeUser(minE[std::make_pair(iBar, Vov)]/5,800);
      if (chLabel == "bar13L" && runs == "5352") h_energy[chLabel]->GetXaxis()->SetRangeUser(60,800);
      if (chLabel == "bar06L" && runs == "5299") h_energy[chLabel]->GetXaxis()->SetRangeUser(40,800);
      int maxbin = h_energy[chLabel]->GetMaximumBin();
      float peak = h_energy[chLabel]->GetBinCenter(maxbin);
      fitLandau[chLabel]->SetRange(peak*0.85, peak*1.2);
      fitLandau[chLabel]->SetParameter(1,peak);
      fitLandau[chLabel]->SetParameter(2,0.1*peak);
      h_energy[chLabel]->Fit(fitLandau[chLabel],"QR");
      fitLandau[chLabel]->SetRange(fitLandau[chLabel]->GetParameter(1)*0.85, fitLandau[chLabel]->GetParameter(1)*1.2);  
      energyMin[chLabel] = fitLandau[chLabel]->GetParameter(1)*0.70;
      //energyMin[chLabel] = fitLandau[chLabel]->GetParameter(1)-2*fitLandau[chLabel]->GetParameter(2);
      h_energy[chLabel]->GetXaxis()->SetRangeUser(0,1000);
      }*/

    for (auto & th: thresholds){
      energyMin_chRef[th] = 0;
      
      if (h_energy_LR[iBar][th]->GetEntries() == 0) continue;
      
      // -- fit average L/R
      if (!source.compare("Laser")){
	fitEnergyLR[iBar][th] = new TF1(Form("fitGaus_bar%02d_th%02d", iBar, th),"gaus", 0, 1000);    
	h_energy_LR[iBar][th]->Fit(fitEnergyLR[iBar][th],"QR");
	energyMinLR[iBar][th] = fitEnergyLR[iBar][th]->GetParameter(1) - 5.0 * fitEnergyLR[iBar][th]->GetParameter(2);
      }
      else{
	fitEnergyLR[iBar][th] = new TF1(Form("fitLandau_bar%02d_th%02d", iBar, th),"landau", 0, 1000);    
	h_energy_LR[iBar][th]->GetXaxis()->SetRangeUser(minE[std::make_pair(iBar, Vov)],800);
	int maxbin = h_energy_LR[iBar][th]->GetMaximumBin();
	float peak = h_energy_LR[iBar][th]->GetBinCenter(maxbin);
	fitEnergyLR[iBar][th]->SetRange(peak*0.85, peak*1.2);
	fitEnergyLR[iBar][th]->SetParameter(1,peak);
	fitEnergyLR[iBar][th]->SetParameter(2,0.1*peak);
	h_energy_LR[iBar][th]->Fit(fitEnergyLR[iBar][th],"QR");
	fitEnergyLR[iBar][th]->SetRange(fitEnergyLR[iBar][th]->GetParameter(1)*0.85, fitEnergyLR[iBar][th]->GetParameter(1)*1.2);  
	energyMinLR[iBar][th] = fitEnergyLR[iBar][th]->GetParameter(1)*0.70;
      }
    
      //std::cout << "bar " iBar << "   th = " << th  << "   minEnergy = " << energyMinLR[iBar][th]<<std::endl; 
      
      h_energy_LR[iBar][th]->GetXaxis()->SetRangeUser(0,1000);
    }
  }
  
  // -- second loop over events 
  cout << "Second loop over events" <<endl;
  for (int entry = 0; entry < maxEntries; entry++){
    
    if( entry%1000 == 0 ) std::cout << ">>> Reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    data->GetEntry(entry);
    
    int th = int(step2/10000)-1;

    acceptEvent_chRef[entry][th] = false;
    for(int iBar = 0; iBar < 16; ++iBar){
      acceptEvent[entry][iBar][th] = false; 
    }
    
    // -- remove showering events
    if (nActiveChannels0[entry][th] > maxActiveChannels) continue;
    if (nActiveChannels1[entry][th] > maxActiveChannels) continue;
    
    //-- ref channel
    if ( channelIdx[chRef1] < 0 ) continue;
    if ( channelIdx[chRef2] < 0 ) continue;
    if ( (*tot)[channelIdx[chRef1]]/1000 < -10. || (*tot)[channelIdx[chRef1]]/1000 > 50. ) continue;      
    if ( (*tot)[channelIdx[chRef2]]/1000 < -10. || (*tot)[channelIdx[chRef2]]/1000 > 50. ) continue;      
    float energyRef = 0.5*((*energy)[channelIdx[chRef1]]+(*energy)[channelIdx[chRef2]]);
    if ( !useTimeAverage) energyRef  = (*energy)[channelIdx[chRef1]];
    if ( energyRef < energyMin_chRef[th] || energyRef > maxEnergy) continue;
    
    acceptEvent_chRef[entry][th] = true;
    
    // -- single channels
    for(int iBar = 0; iBar < 16; ++iBar){        
      
      int chL = chID[Form("bar%02dL", iBar)];
      int chR = chID[Form("bar%02dR", iBar)];
      
      if ( channelIdx[chL] < 0 ) continue;  
      if ( channelIdx[chR] < 0 ) continue;  
      
      if ( (*tot)[channelIdx[chL]]/1000 < -10. || (*tot)[channelIdx[chL]]/1000 > 50. ) continue;
      if ( (*tot)[channelIdx[chR]]/1000 < -10. || (*tot)[channelIdx[chR]]/1000 > 50. ) continue;

      float energyLR = 0.5 * ( (*energy)[channelIdx[chL]] + (*energy)[channelIdx[chR]] );

      if ( energyLR < energyMinLR[iBar][th] || energyLR > maxEnergy ) continue;
      
      acceptEvent[entry][iBar][th] = true;
      
      for (auto label : labelLR ){
	std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
	int ch = chID[chLabel]; 
	h_energy[chLabel][th] -> Fill( (*energy)[channelIdx[ch]] );	
	h_energyRatio[chLabel][th] -> Fill( (*energy)[channelIdx[ch]]/energyRef );	
      }

    }// end loop over bars
  }// end loop over events


  // -- second loop over events to get amp walk corrections
  cout << "Second loop over events to get amp walk corrections" <<endl;
  for (int entry = 0; entry < maxEntries; entry++){
    
    if( entry%1000 == 0 ) std::cout << ">>> Reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    data->GetEntry(entry);

    int th = int(step2/10000)-1;

    if (!acceptEvent_chRef[entry][th]) continue;

    //-- ref channel
    float energyRef = 0.5*((*energy)[channelIdx[chRef1]]+(*energy)[channelIdx[chRef2]]);
    if ( !useTimeAverage) energyRef  = (*energy)[channelIdx[chRef1]];
    long long tRef = 0.5*((*time)[channelIdx[chRef1]]+(*time)[channelIdx[chRef2]]);
    if ( !useTimeAverage) tRef  = (*time)[channelIdx[chRef1]];

    // -- single channels
    for(int iBar = 0; iBar < 16; ++iBar){        

      if (!acceptEvent[entry][iBar][th]) continue;

      int chL = chID[Form("bar%02dL", iBar)];
      int chR = chID[Form("bar%02dR", iBar)];

      for (auto label : labelLR ){
	std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
	int ch = chID[chLabel]; 

	long long deltaT = (*time)[channelIdx[ch]] - tRef;
      
	float energyRatio = (*energy)[channelIdx[ch]]/energyRef ;
	float totRatio = (*tot)[channelIdx[ch]]/(*tot)[channelIdx[chRef1]] ;
	
	if ( fabs(deltaT) > maxDeltaT) continue;
	
	h_deltaT[chLabel][th]  -> Fill( deltaT );	
	h2_energyRatio_vs_totRatio[chLabel][th] -> Fill( totRatio, energyRatio );	


	if (p_deltaT_vs_energyRatio[chLabel][th] == NULL) {
	  float ratioMin = std::max(h_energyRatio[chLabel][th]->GetMean()-5.*h_energyRatio[chLabel][th]->GetRMS(), 0.);
	  float ratioMax = h_energyRatio[chLabel][th]->GetMean()+5.*h_energyRatio[chLabel][th]->GetRMS();
	  p_deltaT_vs_energyRatio[chLabel][th] = new TProfile(Form("p_deltaT_vs_energyRatio_%s_th%02d", chLabel.c_str(), th) , Form("p_deltaT_vs_energyRatio_%s_th%02d", chLabel.c_str(), th), 100, ratioMin, ratioMax);
	  h2_deltaT_vs_energyRatio[chLabel][th] = new TH2F(Form("h2_deltaT_vs_energyRatio_%s_th%02d", chLabel.c_str(), th) , Form("h2_deltaT_vs_energyRatio_%s_th%02d", chLabel.c_str(), th), 100, ratioMin, ratioMax , 1000, -12000, 12000);
	}
	
	p_deltaT_vs_energyRatio[chLabel][th] -> Fill( energyRatio , deltaT );	
	h2_deltaT_vs_energyRatio[chLabel][th] -> Fill( energyRatio , deltaT );	
      }

      //-- tDiff
      if ( acceptEvent[entry][iBar][th] && acceptEvent_chRef[entry][th])    {      
	h_energyRatio_LR[iBar][th]-> Fill( (*energy)[channelIdx[chL]]/(*energy)[channelIdx[chR]] );
	long long deltaT = (*time)[channelIdx[chL]] - (*time)[channelIdx[chR]];
	if ( fabs(deltaT) < maxDeltaT  ){
	  h_deltaT_LR[iBar][th]-> Fill( deltaT );
	  p_deltaT_LR_vs_energyRatio[iBar][th]-> Fill((*energy)[channelIdx[chL]]/(*energy)[channelIdx[chR]], deltaT );
	}
      }
      
    }// end loop over bars

    // fill deltaT L-R for ref bar
    if ( acceptEvent[entry][7][th] && acceptEvent_chRef[entry][th])    {      
      float energyRatio = (*energy)[channelIdx[chRef1]]/(*energy)[channelIdx[chRef2]];
      h_energyRatio_LR_chRef[th] -> Fill( energyRatio);
      h_deltaT_LR_chRef[th] ->Fill(  (*time)[channelIdx[chRef1]] - (*time)[channelIdx[chRef2]] );
      p_deltaT_LR_vs_energyRatio_chRef[th] ->Fill( energyRatio, (*time)[channelIdx[chRef1]] - (*time)[channelIdx[chRef2]] );
    }
    
  }// -- end second loop over entries
  
  
  // ---  amp walk corr
  cout << "Fitting amplitude walk corrections..."<<endl;
  map<std::string, map<int,TF1*> > fitFun_energyRatio;
  map<std::string, map<int,TF1*> > fitFun_energyRatioCorr;
 
  map<int, map<int,TF1*> > fitFun_energyRatio_LR;
  map<int, map<int,TF1*> > fitFun_energyRatioCorr_LR;

  map<int,TF1*>  fitFun_energyRatio_LR_chRef;
  map<int,TF1*>  fitFun_energyRatioCorr_LR_chRef;

  
  for (auto & th: thresholds){
    
    // -- L-R ref bar
    fitFun_energyRatio_LR_chRef[th] = new TF1(Form("fitFun_energyRatio_LR_chRef_th%02d", th), "gaus", 0,10);  
    h_energyRatio_LR_chRef[th] -> Fit(fitFun_energyRatio_LR_chRef[th],"QR");

    fitFun_energyRatioCorr_LR_chRef[th] = new TF1(Form("fitFun_energyRatioCorr_LR_chRef_th%02d", th), "pol4", 0,10);
    fitFun_energyRatioCorr_LR_chRef[th]->SetRange( fitFun_energyRatio_LR_chRef[th]->GetParameter(1)-5*fitFun_energyRatio_LR_chRef[th]->GetParameter(2), fitFun_energyRatio_LR_chRef[th]->GetParameter(1)+5*fitFun_energyRatio_LR_chRef[th]->GetParameter(2));
    p_deltaT_LR_vs_energyRatio_chRef[th] -> Fit(fitFun_energyRatioCorr_LR_chRef[th],"QRS"); 
    
    for(int iBar = 0; iBar < 16; ++iBar){

      // -- L-R
      fitFun_energyRatio_LR[iBar][th] = new TF1(Form("fitFun_energyRatio_LR_bar%02d_th%02d", iBar, th), "gaus", 0,10);  
      h_energyRatio_LR[iBar][th] -> Fit(fitFun_energyRatio_LR[iBar][th],"QR");
      
      fitFun_energyRatioCorr_LR[iBar][th] = new TF1(Form("fitFun_energyRatioCorr_LR_bar%02d_th%02d", iBar, th), "pol4", 0,10);
      fitFun_energyRatioCorr_LR[iBar][th]->SetRange( fitFun_energyRatio_LR[iBar][th]->GetParameter(1)-5*fitFun_energyRatio_LR[iBar][th]->GetParameter(2), fitFun_energyRatio_LR[iBar][th]->GetParameter(1)+5*fitFun_energyRatio_LR[iBar][th]->GetParameter(2));
      p_deltaT_LR_vs_energyRatio[iBar][th] -> Fit(fitFun_energyRatioCorr_LR[iBar][th],"QRS"); 

      // -- single channels 
      for (auto label : labelLR ){
	std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
	int ch = chID[chLabel]; 
	
	if ( p_deltaT_vs_energyRatio[chLabel][th] == NULL) continue;
	
	fitFun_energyRatio[chLabel][th] = new TF1(Form("fitFun_energyRatio_%s_th%02d",chLabel.c_str(), th), "gaus", 0,10);  
	if ( h_energyRatio[chLabel][th] -> GetEntries() ==  0) continue;
	h_energyRatio[chLabel][th] -> Fit(fitFun_energyRatio[chLabel][th],"QR");
	
	//fitFun_energyRatioCorr[chLabel][th] = new TF1(Form("fitFun_energyRatioCorr_ch%02d_th%02d",ch, th), "[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*pow(log(x),3)", 0,10);
	fitFun_energyRatioCorr[chLabel][th] = new TF1(Form("fitFun_energyRatioCorr_ch%02d_th%02d",ch, th), "[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*pow(log(x),3)+[4]*pow(log(x),4)", 0,10);
	fitFun_energyRatioCorr[chLabel][th]->SetParameters( p_deltaT_vs_energyRatio[chLabel][th] -> GetMean(2), -100, -50, -50 );
	
	p_deltaT_vs_energyRatio[chLabel][th] -> Fit(fitFun_energyRatioCorr[chLabel][th],"QRS");
      }
    }
  }
  
  // -- third loop over events to apply amp walk corrections
  cout << "Third loop over events to apply amp walk corrections" <<endl;
  for (int entry = 0; entry < maxEntries; entry++){
    
    if( entry%1000 == 0 ) std::cout << ">>> Reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    data->GetEntry(entry);

    int th = int(step2/10000)-1;

    if ( !acceptEvent_chRef[entry][th]) continue;

    float energyRef = 0.5*((*energy)[channelIdx[chRef1]]+(*energy)[channelIdx[chRef2]]);
    long long tRef = 0.5*((*time)[channelIdx[chRef1]]+(*time)[channelIdx[chRef2]]);
    if ( !useTimeAverage) {
      energyRef  = (*energy)[channelIdx[chRef1]];
      tRef       = (*time)[channelIdx[chRef1]];
    }

    // -- single channels
    for(int iBar = 0; iBar < 16; ++iBar){        
      
      if (!acceptEvent[entry][iBar][th]) continue;
      
      for (auto label : labelLR ){
	std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
	int ch = chID[chLabel]; 
	
	if ( fitFun_energyRatioCorr[chLabel][th] == NULL) continue;
	
	float energyRatio = (*energy)[channelIdx[ch]]/energyRef;
	float energyRatioCorr = fitFun_energyRatioCorr[chLabel][th] -> Eval( energyRatio ) - fitFun_energyRatioCorr[chLabel][th] -> Eval( h_energyRatio[chLabel][th] -> GetMean() ); 
	long long deltaT = (*time)[channelIdx[ch]] - tRef;

	if ( fabs(deltaT) < maxDeltaT && fabs(deltaT-energyRatioCorr) < maxDeltaT){   
	  h_deltaT_energyRatioCorr[chLabel][th] -> Fill( deltaT - energyRatioCorr);
	  p_deltaT_energyRatioCorr_vs_t1fine[chLabel][th]  -> Fill( (*t1fine)[channelIdx[ch]] , deltaT - energyRatioCorr );
	  h2_deltaT_energyRatioCorr_vs_t1fine[chLabel][th] -> Fill( (*t1fine)[channelIdx[ch]] , deltaT - energyRatioCorr );
	}
      } // end loop L,R


      // tDiff L-R
      int chL = chID[Form("bar%02dL", iBar)];
      int chR = chID[Form("bar%02dR", iBar)];

      if ( acceptEvent[entry][iBar][th] && acceptEvent_chRef[entry][th])    {      
	float energyRatio = (*energy)[channelIdx[chL]]/(*energy)[channelIdx[chR]];
	float energyRatioCorr = fitFun_energyRatioCorr_LR[iBar][th] -> Eval( energyRatio ) - fitFun_energyRatioCorr_LR[iBar][th] -> Eval( h_energyRatio_LR[iBar][th]->GetMean() ); 
	long long deltaT = (*time)[channelIdx[chL]] - (*time)[channelIdx[chR]];
	
	if ( fabs(deltaT) < maxDeltaT && fabs(deltaT - energyRatioCorr) < maxDeltaT){
	  h_deltaT_LR_energyRatioCorr[iBar][th]-> Fill( deltaT - energyRatioCorr);
	  p_deltaT_LR_energyRatioCorr_vs_t1fine[iBar][th]-> Fill( 0.5*( (*t1fine)[channelIdx[chL]]+(*t1fine)[channelIdx[chR]]), deltaT - energyRatioCorr);
	}
      }

    }// end loop over bars


    // fill deltaT L-R for ref bar
    if ( acceptEvent[entry][7][th] && acceptEvent_chRef[entry][th])    {      
      float energyRatio = (*energy)[channelIdx[chRef1]] / (*energy)[channelIdx[chRef2]] ;
      float energyRatioCorr = fitFun_energyRatioCorr_LR_chRef[th] -> Eval( energyRatio) - fitFun_energyRatioCorr_LR_chRef[th] -> Eval( h_energyRatio_LR_chRef[th]->GetMean() );
      h_deltaT_LR_energyRatioCorr_chRef[th] ->Fill(  (*time)[channelIdx[chRef1]] - (*time)[channelIdx[chRef2]] - energyRatioCorr );
      p_deltaT_LR_energyRatioCorr_vs_t1fine_chRef[th] ->Fill( 0.5*( (*t1fine)[channelIdx[chRef1]]+(*t1fine)[channelIdx[chRef2]]  )    , (*time)[channelIdx[chRef1]] - (*time)[channelIdx[chRef2]] - energyRatioCorr);
    }

  }// end loop over events


  // -- fourth loop over events to apply phase corrections
  cout << "Fourth loop over events to apply phase corrections "  << endl;
  for (int entry = 0; entry < maxEntries; entry++){

    if( entry%1000 == 0 ) std::cout << ">>> Reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    data->GetEntry(entry); 

    int th = int(step2/10000)-1;

    if ( !acceptEvent_chRef[entry][th]) continue;
    
    float energyRef = 0.5*((*energy)[channelIdx[chRef1]]+(*energy)[channelIdx[chRef2]]);
    long long tRef = 0.5*((*time)[channelIdx[chRef1]]+(*time)[channelIdx[chRef2]]);
    if ( !useTimeAverage) {
      energyRef  = (*energy)[channelIdx[chRef1]];
      tRef       = (*time)[channelIdx[chRef1]];
    }

    // -- single channels    
    for(int iBar = 0; iBar < 16; ++iBar){        
      
      if ( !acceptEvent[entry][iBar][th] ) continue;
      
      for (auto label : labelLR ){
	std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
	int ch = chID[chLabel]; 
	
	long long deltaT = (*time)[channelIdx[ch]] - tRef;
	
	if ( fitFun_energyRatioCorr[chLabel][th] == NULL) continue;

	//float energyRatio = (*tot)[channelIdx[ch]]/(*tot)[channelIdx[chRef]] ;
	float energyRatio = (*energy)[channelIdx[ch]]/energyRef;
	float energyRatioCorr = fitFun_energyRatioCorr[chLabel][th] -> Eval( energyRatio ) - fitFun_energyRatioCorr[chLabel][th] -> Eval( h_energyRatio[chLabel][th] -> GetMean() );
	int bin1  = p_deltaT_energyRatioCorr_vs_t1fine[chLabel][th]->FindBin( (*t1fine)[channelIdx[ch]]) ; 
	int bin2 = p_deltaT_energyRatioCorr_vs_t1fine[chLabel][th]->FindBin( p_deltaT_energyRatioCorr_vs_t1fine[chLabel][th]->GetMean() );
	float phaseCorr = p_deltaT_energyRatioCorr_vs_t1fine[chLabel][th] -> GetBinContent(bin1) -  p_deltaT_energyRatioCorr_vs_t1fine[chLabel][th] -> GetBinContent(bin2);
	
	if ( fabs(deltaT) < maxDeltaT && fabs(deltaT-energyRatioCorr) < maxDeltaT){   
	  h_deltaT_energyRatioCorr_phaseCorr[chLabel][th] -> Fill( deltaT - energyRatioCorr - phaseCorr);
	  p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel][th]  -> Fill( (*t1fine)[channelIdx[chRef2]] , deltaT - energyRatioCorr - phaseCorr);
	  h2_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel][th] -> Fill( (*t1fine)[channelIdx[chRef2]] , deltaT - energyRatioCorr - phaseCorr);
	}
      } // end loop L,R
      
      
      // tDiff L-R
      int chL = chID[Form("bar%02dL", iBar)];
      int chR = chID[Form("bar%02dR", iBar)];
      
      if ( acceptEvent[entry][iBar][th] && acceptEvent_chRef[entry][th])    {      
	float energyRatio = (*energy)[channelIdx[chL]]/(*energy)[channelIdx[chR]];
	float energyRatioCorr = fitFun_energyRatioCorr_LR[iBar][th] -> Eval( energyRatio ) - fitFun_energyRatioCorr_LR[iBar][th] -> Eval( h_energyRatio_LR[iBar][th]->GetMean() ); 
	long long deltaT = (*time)[channelIdx[chL]] - (*time)[channelIdx[chR]];
	int bin1  = p_deltaT_LR_energyRatioCorr_vs_t1fine[iBar][th]->FindBin( 0.5 * ( (*t1fine)[channelIdx[chL]]+(*t1fine)[channelIdx[chR]]) ); 
	int bin2 = p_deltaT_LR_energyRatioCorr_vs_t1fine[iBar][th]->FindBin( p_deltaT_LR_energyRatioCorr_vs_t1fine[iBar][th]->GetMean() );
	float phaseCorr = p_deltaT_LR_energyRatioCorr_vs_t1fine[iBar][th] -> GetBinContent(bin1) -  p_deltaT_LR_energyRatioCorr_vs_t1fine[iBar][th] -> GetBinContent(bin2);
	
	if ( fabs(deltaT)<maxDeltaT && fabs(deltaT-energyRatioCorr)<maxDeltaT){
	  h_deltaT_LR_energyRatioCorr_phaseCorr[iBar][th] -> Fill( deltaT - energyRatioCorr - phaseCorr);
	}
      }
    }    // end loop over bars


    // fill deltaT L-R for ref bar
    if ( acceptEvent[entry][7][th] && acceptEvent_chRef[entry][th])    {      
      float energyRatio = (*energy)[channelIdx[chRef1]] / (*energy)[channelIdx[chRef2]] ;
      float energyRatioCorr = fitFun_energyRatioCorr_LR_chRef[th] -> Eval( energyRatio ) - fitFun_energyRatioCorr_LR_chRef[th] -> Eval( h_energyRatio_LR_chRef[th]->GetMean() );
      int bin1 = p_deltaT_LR_energyRatioCorr_vs_t1fine_chRef[th]->FindBin( 0.5 * ( (*t1fine)[channelIdx[chRef1]]+(*t1fine)[channelIdx[chRef2]]) ); 
      int bin2 = p_deltaT_LR_energyRatioCorr_vs_t1fine_chRef[th]->FindBin( p_deltaT_LR_energyRatioCorr_vs_t1fine_chRef[th]->GetMean() );
      float phaseCorr = p_deltaT_LR_energyRatioCorr_vs_t1fine_chRef[th] -> GetBinContent(bin1) -  p_deltaT_LR_energyRatioCorr_vs_t1fine_chRef[th] -> GetBinContent(bin2);      
      h_deltaT_LR_energyRatioCorr_phaseCorr_chRef[th] ->Fill(  (*time)[channelIdx[chRef1]] - (*time)[channelIdx[chRef2]] - energyRatioCorr - phaseCorr);
    }
  }
  

  // -- gaus fit deltaT for each channel
  map<int, TGraphErrors*> g_tRes_LR_vs_bar; 
  map<int, TGraphErrors*> g_tRes_energyRatioCorr_LR_vs_bar;
  map<int, TGraphErrors*> g_tRes_energyRatioCorr_phaseCorr_LR_vs_bar;

  map<int, map<int,TF1*> > fitGaus_LR;
  map<int, map<int,TF1*> > fitGaus_energyRatioCorr_LR;
  map<int, map<int,TF1*> > fitGaus_energyRatioCorr_phaseCorr_LR;
  

  map<std::string, map<int,TGraphErrors*> > g_tRes_vs_bar;
  map<std::string, map<int,TGraphErrors*> > g_tRes_energyRatioCorr_vs_bar;
  map<std::string, map<int,TGraphErrors*> > g_tRes_energyRatioCorr_phaseCorr_vs_bar;

  map<std::string, map<int,TGraphErrors*> > g_tRes_vs_th;
  map<std::string, map<int,TGraphErrors*> > g_tRes_energyRatioCorr_vs_th;
  map<std::string, map<int,TGraphErrors*> > g_tRes_energyRatioCorr_phaseCorr_vs_th;

  map<std::string, map<int,TGraphErrors*> > g_tRes_vs_th_final; // after subtracting sigma_ref
  map<std::string, map<int,TGraphErrors*> > g_tRes_energyRatioCorr_vs_th_final;
  map<std::string, map<int,TGraphErrors*> > g_tRes_energyRatioCorr_phaseCorr_vs_th_final;

  map<std::string, map<int,TF1*> > fitGaus;
  map<std::string, map<int,TF1*> > fitGaus_energyRatioCorr;
  map<std::string, map<int,TF1*> > fitGaus_energyRatioCorr_phaseCorr;


  map<int,TF1*> fitGaus_LR_chRef;
  map<int,TF1*> fitGaus_energyRatioCorr_LR_chRef;
  map<int,TF1*> fitGaus_energyRatioCorr_phaseCorr_LR_chRef;

  // delta T L-R
  std::cout << "Time resolution from tDiff(L-R)"<<std::endl;

  for (auto & th : thresholds){
    
    //L-R chRef
    fitGaus_LR_chRef[th] = new TF1(Form("fitGaus_LR_chRef_th%02d", th), "gaus",-maxDeltaT,maxDeltaT);
    getTimeResolution(h_deltaT_LR_chRef[th], fitGaus_LR_chRef[th]);    

    fitGaus_energyRatioCorr_LR_chRef[th] = new TF1(Form("fitGaus_energyRatioCorr_LR_chRef_th%02d", th), "gaus",-maxDeltaT,maxDeltaT);
    getTimeResolution(h_deltaT_LR_energyRatioCorr_chRef[th], fitGaus_energyRatioCorr_LR_chRef[th]);    
    
    fitGaus_energyRatioCorr_phaseCorr_LR_chRef[th] = new TF1(Form("fitGaus_energyRatioCorr_phaseCorr_LR_chRef_th%02d", th), "gaus",-maxDeltaT,maxDeltaT);
    getTimeResolution(h_deltaT_LR_energyRatioCorr_phaseCorr_chRef[th], fitGaus_energyRatioCorr_phaseCorr_LR_chRef[th]);    

    /// L-R
    g_tRes_LR_vs_bar[th] = new TGraphErrors();
    g_tRes_energyRatioCorr_LR_vs_bar[th] = new TGraphErrors();
    g_tRes_energyRatioCorr_phaseCorr_LR_vs_bar[th] = new TGraphErrors();
    
    for(int iBar = 0; iBar < 16; ++iBar){
      
      if ( h_deltaT_LR[iBar][th] -> GetEntries() == 0) continue;
      
      // -- no corr                                                                                                                                                         
      fitGaus_LR[iBar][th] = new TF1(Form("fitGaus_LR_bar%02d_th%02d", iBar, th), "gaus",-maxDeltaT,maxDeltaT);                                                                                getTimeResolution(h_deltaT_LR[iBar][th], fitGaus_LR[iBar][th]);
      g_tRes_LR_vs_bar[th] -> SetPoint( g_tRes_LR_vs_bar[th]->GetN(), iBar, fitGaus_LR[iBar][th]->GetParameter(2));
      g_tRes_LR_vs_bar[th] -> SetPointError( g_tRes_LR_vs_bar[th]->GetN()-1, 0, fitGaus_LR[iBar][th]->GetParError(2));

      // -- energy corr
      fitGaus_energyRatioCorr_LR[iBar][th] = new TF1(Form("fitGaus_energyRatioCorr_LR_bar%02d_th%02d", iBar, th), "gaus",-maxDeltaT,maxDeltaT);
      getTimeResolution(h_deltaT_LR_energyRatioCorr[iBar][th], fitGaus_energyRatioCorr_LR[iBar][th]); 
      g_tRes_energyRatioCorr_LR_vs_bar[th] -> SetPoint( g_tRes_energyRatioCorr_LR_vs_bar[th]->GetN(), iBar, fitGaus_energyRatioCorr_LR[iBar][th]->GetParameter(2));
      g_tRes_energyRatioCorr_LR_vs_bar[th] -> SetPointError( g_tRes_energyRatioCorr_LR_vs_bar[th]->GetN()-1, 0, fitGaus_energyRatioCorr_LR[iBar][th]->GetParError(2));
      
      // -- energy + phase corr
      fitGaus_energyRatioCorr_phaseCorr_LR[iBar][th] = new TF1(Form("fitGaus_energyRatioCorr_phaseCorr_LR_bar%02d_th%02d", iBar, th), "gaus",-maxDeltaT,maxDeltaT);
      getTimeResolution(h_deltaT_LR_energyRatioCorr_phaseCorr[iBar][th], fitGaus_energyRatioCorr_phaseCorr_LR[iBar][th]); 
      g_tRes_energyRatioCorr_phaseCorr_LR_vs_bar[th]-> SetPoint( g_tRes_energyRatioCorr_phaseCorr_LR_vs_bar[th]->GetN(), iBar, fitGaus_energyRatioCorr_phaseCorr_LR[iBar][th]->GetParameter(2));
      g_tRes_energyRatioCorr_phaseCorr_LR_vs_bar[th]-> SetPointError( g_tRes_energyRatioCorr_phaseCorr_LR_vs_bar[th]->GetN()-1, 0, fitGaus_energyRatioCorr_phaseCorr_LR[iBar][th]->GetParError(2));
      
    }
  }



  // deltaT ch - Ref
  std::cout << "Time resolution from tDiff(ch-Ref)"<<std::endl;

  for (auto label : labelLR ){
    for (auto & th : thresholds) {
      g_tRes_vs_bar[label][th] = new TGraphErrors();
      g_tRes_energyRatioCorr_vs_bar[label][th] = new TGraphErrors();
      g_tRes_energyRatioCorr_phaseCorr_vs_bar[label][th] = new TGraphErrors();  
    }
	   
    for(int iBar = 0; iBar < 16; ++iBar){        
      g_tRes_vs_th[label][iBar] = new TGraphErrors();
      g_tRes_energyRatioCorr_vs_th[label][iBar] = new TGraphErrors();
      g_tRes_energyRatioCorr_phaseCorr_vs_th[label][iBar] = new TGraphErrors();  

      g_tRes_vs_th_final[label][iBar] = new TGraphErrors();
      g_tRes_energyRatioCorr_vs_th_final[label][iBar] = new TGraphErrors();
      g_tRes_energyRatioCorr_phaseCorr_vs_th_final[label][iBar] = new TGraphErrors();  
    }
  }
  
  for (auto & th : thresholds) {
    for (auto label : labelLR ){
      for(int iBar = 0; iBar < 16; ++iBar){        
	
	std::string chLabel = Form("bar%02d%s", iBar, label.c_str());

	if ( h_deltaT[chLabel][th] -> GetEntries() == 0) continue;
	
	// -- no corr
	fitGaus[chLabel][th] = new TF1(Form("fitGaus_%s_th%02d",chLabel.c_str(), th), "gaus",-maxDeltaT,maxDeltaT);
	getTimeResolution(h_deltaT[chLabel][th], fitGaus[chLabel][th]);
	g_tRes_vs_bar[label][th]  -> SetPoint( g_tRes_vs_bar[label][th]->GetN(), iBar, fitGaus[chLabel][th]->GetParameter(2));
	g_tRes_vs_bar[label][th]  -> SetPointError( g_tRes_vs_bar[label][th]->GetN()-1, 0, fitGaus[chLabel][th]->GetParError(2));
	g_tRes_vs_th[label][iBar] -> SetPoint( g_tRes_vs_th[label][iBar]->GetN(), th, fitGaus[chLabel][th]->GetParameter(2));
	g_tRes_vs_th[label][iBar] -> SetPointError( g_tRes_vs_th[label][iBar]->GetN()-1, 0, fitGaus[chLabel][th]->GetParError(2));
	float corr_sigma = 0;
	if ( fitGaus[chLabel][th]->GetParameter(2) > (0.5*fitGaus_LR_chRef[th]->GetParameter(2))) {
	  corr_sigma = sqrt( pow(fitGaus[chLabel][th]->GetParameter(2),2) - pow(0.5*fitGaus_LR_chRef[th]->GetParameter(2),2));  
	  g_tRes_vs_th_final[label][iBar] -> SetPoint( g_tRes_vs_th_final[label][iBar]->GetN(), th, corr_sigma);	  
	}

	// -- energy corr
	fitGaus_energyRatioCorr[chLabel][th] = new TF1(Form("fitGaus_energyRatioCorr_%s_th%02d",chLabel.c_str(),th), "gaus",-maxDeltaT,maxDeltaT);
	getTimeResolution(h_deltaT_energyRatioCorr[chLabel][th], fitGaus_energyRatioCorr[chLabel][th]); 
	g_tRes_energyRatioCorr_vs_bar[label][th]  -> SetPoint( g_tRes_energyRatioCorr_vs_bar[label][th]->GetN(), iBar, fitGaus_energyRatioCorr[chLabel][th]->GetParameter(2));
	g_tRes_energyRatioCorr_vs_bar[label][th]  -> SetPointError( g_tRes_energyRatioCorr_vs_bar[label][th]->GetN()-1, 0, fitGaus_energyRatioCorr[chLabel][th]->GetParError(2));
	g_tRes_energyRatioCorr_vs_th[label][iBar] -> SetPoint( g_tRes_energyRatioCorr_vs_th[label][iBar]->GetN(), th, fitGaus_energyRatioCorr[chLabel][th]->GetParameter(2));
	g_tRes_energyRatioCorr_vs_th[label][iBar] -> SetPointError( g_tRes_energyRatioCorr_vs_th[label][iBar]->GetN()-1, 0, fitGaus_energyRatioCorr[chLabel][th]->GetParError(2));
	if ( fitGaus_energyRatioCorr[chLabel][th]->GetParameter(2) > (0.5*fitGaus_energyRatioCorr_LR_chRef[th]->GetParameter(2))) {
	  corr_sigma = sqrt( pow(fitGaus_energyRatioCorr[chLabel][th]->GetParameter(2),2) - pow(0.5*fitGaus_energyRatioCorr_LR_chRef[th]->GetParameter(2),2));  
	  g_tRes_energyRatioCorr_vs_th_final[label][iBar] -> SetPoint( g_tRes_energyRatioCorr_vs_th_final[label][iBar]->GetN(), th, corr_sigma);	  
	}
	
	// -- energy + phase corr
	fitGaus_energyRatioCorr_phaseCorr[chLabel][th] = new TF1(Form("fitGaus_energyRatioCorr_phaseCorr_%s_th%02d",chLabel.c_str(),th), "gaus",-maxDeltaT,maxDeltaT);
	getTimeResolution(h_deltaT_energyRatioCorr_phaseCorr[chLabel][th], fitGaus_energyRatioCorr_phaseCorr[chLabel][th]);
	g_tRes_energyRatioCorr_phaseCorr_vs_bar[label][th]  -> SetPoint( g_tRes_energyRatioCorr_phaseCorr_vs_bar[label][th]->GetN(), iBar, fitGaus_energyRatioCorr_phaseCorr[chLabel][th]->GetParameter(2));
	g_tRes_energyRatioCorr_phaseCorr_vs_bar[label][th]  -> SetPointError( g_tRes_energyRatioCorr_phaseCorr_vs_bar[label][th]->GetN()-1, 0, fitGaus_energyRatioCorr_phaseCorr[chLabel][th]->GetParError(2));
	g_tRes_energyRatioCorr_phaseCorr_vs_th[label][iBar] -> SetPoint( g_tRes_energyRatioCorr_phaseCorr_vs_th[label][iBar]->GetN(), th, fitGaus_energyRatioCorr_phaseCorr[chLabel][th]->GetParameter(2));
	g_tRes_energyRatioCorr_phaseCorr_vs_th[label][iBar] -> SetPointError( g_tRes_energyRatioCorr_phaseCorr_vs_th[label][iBar]->GetN()-1, 0, fitGaus_energyRatioCorr_phaseCorr[chLabel][th]->GetParError(2));
	if ( fitGaus_energyRatioCorr_phaseCorr[chLabel][th]->GetParameter(2) > (0.5*fitGaus_energyRatioCorr_phaseCorr_LR_chRef[th]->GetParameter(2))) {
	  corr_sigma = sqrt( pow(fitGaus_energyRatioCorr_phaseCorr[chLabel][th]->GetParameter(2),2) - pow(0.5*fitGaus_energyRatioCorr_phaseCorr_LR_chRef[th]->GetParameter(2),2));  
	  g_tRes_energyRatioCorr_phaseCorr_vs_th_final[label][iBar] -> SetPoint( g_tRes_energyRatioCorr_phaseCorr_vs_th_final[label][iBar]->GetN(), th, corr_sigma);	  
	}	
      }
    }
  }
  
  
  // --- Triangulation L, R, chRef
  std::cout << "Triangulation ..." << std::endl;
  map<int, TGraphErrors*> g_tRes_L_vs_bar_fromTriangulation;
  map<int, TGraphErrors*> g_tRes_R_vs_bar_fromTriangulation;
  map<int, TGraphErrors*> g_tRes_chRef_vs_bar_fromTriangulation;

  for (auto & th : thresholds){
    g_tRes_L_vs_bar_fromTriangulation[th] = new TGraphErrors();
    g_tRes_R_vs_bar_fromTriangulation[th] = new TGraphErrors();
    g_tRes_chRef_vs_bar_fromTriangulation[th] = new TGraphErrors();
  }

  map<int, TGraphErrors*> g_tRes_L_vs_th_fromTriangulation;
  map<int, TGraphErrors*> g_tRes_R_vs_th_fromTriangulation;
  map<int, TGraphErrors*> g_tRes_chRef_vs_th_fromTriangulation;

  for (int iBar = 0 ; iBar < 16; iBar++){
    g_tRes_L_vs_th_fromTriangulation[iBar] = new TGraphErrors();
    g_tRes_R_vs_th_fromTriangulation[iBar] = new TGraphErrors();
    g_tRes_chRef_vs_th_fromTriangulation[iBar] = new TGraphErrors();
  }
 
  for (auto & th : thresholds){

    for (int iBar = 0 ; iBar < 16; iBar++){
      
      if (h_deltaT_LR_energyRatioCorr_phaseCorr[iBar][th] -> GetEntries() == 0) continue;
      if (h_deltaT_energyRatioCorr_phaseCorr[ Form("bar%02dR", iBar)][th] -> GetEntries() == 0) continue;
      if (h_deltaT_energyRatioCorr_phaseCorr[ Form("bar%02dL", iBar)][th] -> GetEntries() == 0) continue;
      
      float tRes_L_R = 0;
      float tRes_L_chRef = 0;
      float tRes_R_chRef = 0;
      
      float err_tRes_L_R = 0;
      float err_tRes_L_chRef = 0;
      float err_tRes_R_chRef = 0;
    
      for ( int i = 0; i < g_tRes_energyRatioCorr_phaseCorr_LR_vs_bar[th]-> GetN(); i++){
	if (  int(g_tRes_energyRatioCorr_phaseCorr_LR_vs_bar[th] -> GetPointX(i)) == iBar) {
	  tRes_L_R = g_tRes_energyRatioCorr_phaseCorr_LR_vs_bar[th] -> GetPointY(i);
	  err_tRes_L_R = g_tRes_energyRatioCorr_phaseCorr_LR_vs_bar[th] -> GetErrorY(i);
	  break;
	}
      }

      for ( int i = 0; i < g_tRes_energyRatioCorr_phaseCorr_vs_bar["L"][th]-> GetN(); i++){
	if (  g_tRes_energyRatioCorr_phaseCorr_vs_bar["L"][th]-> GetPointX(i) == iBar) {
	  tRes_L_chRef = g_tRes_energyRatioCorr_phaseCorr_vs_bar["L"][th]-> GetPointY(i);
	  err_tRes_L_chRef = g_tRes_energyRatioCorr_phaseCorr_vs_bar["L"][th]-> GetErrorY(i);
	  break;
	}
      }
      
      for ( int i = 0; i < g_tRes_energyRatioCorr_phaseCorr_vs_bar["R"][th]-> GetN(); i++){
	if (  g_tRes_energyRatioCorr_phaseCorr_vs_bar["R"][th]-> GetPointX(i) == iBar) {
	  tRes_R_chRef = g_tRes_energyRatioCorr_phaseCorr_vs_bar["R"][th]-> GetPointY(i);
	  err_tRes_R_chRef = g_tRes_energyRatioCorr_phaseCorr_vs_bar["R"][th]-> GetErrorY(i);
	  break;
	}
      }
      
      cout <<  "  ***   bar " << iBar << "  th = " << th <<  "  sigma(L-R) = " << tRes_L_R << "  sigma(R-ref) = " << tRes_R_chRef  << "  sigma(L-ref) = " << tRes_L_chRef <<endl;
      if ( tRes_L_R==0 || tRes_R_chRef==0 || tRes_L_chRef==0) continue;    
      
      std::vector<float> tRes, err_tRes ;
      tRes.clear();
      err_tRes.clear();
      
      ilTriangoloNo(tRes_L_R,tRes_R_chRef,tRes_L_chRef, err_tRes_L_R, err_tRes_R_chRef, err_tRes_L_chRef, tRes, err_tRes);
      
      // vs bar
      if ( !isnan(tRes[0]) ) g_tRes_L_vs_bar_fromTriangulation[th] -> SetPoint( g_tRes_L_vs_bar_fromTriangulation[th]->GetN(), iBar, tRes[0]);
      if ( !isnan(tRes[1]) ) g_tRes_R_vs_bar_fromTriangulation[th] -> SetPoint( g_tRes_R_vs_bar_fromTriangulation[th]->GetN(), iBar, tRes[1]);
      if ( !isnan(tRes[2]) ) g_tRes_chRef_vs_bar_fromTriangulation[th] -> SetPoint( g_tRes_chRef_vs_bar_fromTriangulation[th]->GetN(), iBar, tRes[2]);
      
      if ( err_tRes[0]>=0 ) g_tRes_L_vs_bar_fromTriangulation[th] -> SetPointError( g_tRes_L_vs_bar_fromTriangulation[th]->GetN()-1, 0, err_tRes[0]);
      if ( err_tRes[1]>=0 ) g_tRes_R_vs_bar_fromTriangulation[th] -> SetPointError( g_tRes_R_vs_bar_fromTriangulation[th]->GetN()-1, 0, err_tRes[1]);
      if ( err_tRes[2]>=0 ) g_tRes_chRef_vs_bar_fromTriangulation[th] -> SetPointError( g_tRes_chRef_vs_bar_fromTriangulation[th]->GetN()-1, 0, err_tRes[2]);


      // vs th
      if ( !isnan(tRes[0]) ) g_tRes_L_vs_th_fromTriangulation[iBar] -> SetPoint( g_tRes_L_vs_th_fromTriangulation[iBar]->GetN(), th, tRes[0]);
      if ( !isnan(tRes[1]) ) g_tRes_R_vs_th_fromTriangulation[iBar] -> SetPoint( g_tRes_R_vs_th_fromTriangulation[iBar]->GetN(), th, tRes[1]);
      if ( !isnan(tRes[2]) ) g_tRes_chRef_vs_th_fromTriangulation[iBar] -> SetPoint( g_tRes_chRef_vs_th_fromTriangulation[iBar]->GetN(), th, tRes[2]);
      
      if ( err_tRes[0]>=0 ) g_tRes_L_vs_th_fromTriangulation[iBar] -> SetPointError( g_tRes_L_vs_th_fromTriangulation[iBar]->GetN()-1, 0, err_tRes[0]);
      if ( err_tRes[1]>=0 ) g_tRes_R_vs_th_fromTriangulation[iBar] -> SetPointError( g_tRes_R_vs_th_fromTriangulation[iBar]->GetN()-1, 0, err_tRes[1]);
      if ( err_tRes[2]>=0 ) g_tRes_chRef_vs_th_fromTriangulation[iBar] -> SetPointError( g_tRes_chRef_vs_th_fromTriangulation[iBar]->GetN()-1, 0, err_tRes[2]);
    }
  }




  // ======  save histograms in a file
  string foutName = Form("plots/analysisSingleChannel_runs%s.root",runs.c_str());
  if (useTimeAverage) foutName = Form("plots/analysisSingleChannel_runs%s_timeAverage.root",runs.c_str());
  std::cout<< "Saving histograms on file " << foutName.c_str() << std::endl;

  
  TFile *fout = new TFile(foutName.c_str(),"recreate");

  for (auto & th : thresholds ){
    for (auto label : labelLR ){
      for(int iBar = 0; iBar < 16; ++iBar){        
	std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
	
	if (h_energy[chLabel][th]-> GetEntries() == 0) continue;

	h_energyRatio[chLabel][th]->Write();
	h2_energyRatio_vs_totRatio[chLabel][th]->Write();
	h_energy[chLabel][th]->Write();
	h_deltaT[chLabel][th]->Write();
	h_deltaT_energyRatioCorr[chLabel][th]->Write();
	h_deltaT_energyRatioCorr_phaseCorr[chLabel][th]->Write();
	if ( p_deltaT_vs_energyRatio[chLabel][th] ) {
	  p_deltaT_vs_energyRatio[chLabel][th]->Write();
	  h2_deltaT_vs_energyRatio[chLabel][th]->Write();
	}
      }
      if (g_tRes_vs_bar[label][th]->GetN() == 0) continue;
      g_tRes_vs_bar[label][th]->Write(Form("g_tRes_vs_bar_%s_th%02d",label.c_str(), th));
      g_tRes_energyRatioCorr_vs_bar[label][th]->Write(Form("g_tRes_energyRatioCorr_%s_vs_bar_th%02d",label.c_str(), th));
      g_tRes_energyRatioCorr_phaseCorr_vs_bar[label][th]->Write( Form("g_tRes_energyRatioCorr_phaseCorr_%s_vs_bar_th%02d",label.c_str(), th));
    }
  }
  
  for (auto label : labelLR ){
    for(int iBar = 0; iBar < 16; ++iBar){        
      std::string chLabel = Form("bar%02d%s", iBar, label.c_str());
      if (g_tRes_vs_th[label][iBar]->GetN() == 0) continue;
      g_tRes_vs_th[label][iBar]->Write(Form("g_tRes_vs_th_%s",chLabel.c_str()));
      g_tRes_energyRatioCorr_vs_th[label][iBar]->Write(Form("g_tRes_energyRatioCorr_vs_th_%s",chLabel.c_str()));
      g_tRes_energyRatioCorr_phaseCorr_vs_th[label][iBar]->Write( Form("g_tRes_energyRatioCorr_phaseCorr_vs_th_%s",chLabel.c_str()));

      g_tRes_vs_th_final[label][iBar]->Write(Form("g_tRes_vs_th_%s_final",chLabel.c_str()));
      g_tRes_energyRatioCorr_vs_th_final[label][iBar]->Write(Form("g_tRes_energyRatioCorr_vs_th_%s_final",chLabel.c_str()));
      g_tRes_energyRatioCorr_phaseCorr_vs_th_final[label][iBar]->Write( Form("g_tRes_energyRatioCorr_phaseCorr_vs_th_%s_final",chLabel.c_str()));
    }
  }

  for (auto & th : thresholds ){
    if ( g_tRes_L_vs_bar_fromTriangulation[th] -> GetN() == 0) continue;
    g_tRes_L_vs_bar_fromTriangulation[th] -> Write(Form("g_tRes_L_vs_bar_fromTriangulation_th%02d",th));
    g_tRes_R_vs_bar_fromTriangulation[th] -> Write(Form("g_tRes_R_vs_bar_fromTriangulation_th%02d",th));
    g_tRes_chRef_vs_bar_fromTriangulation[th] -> Write(Form("g_tRes_chRef_vs_bar_fromTriangulation_th%02d",th));
  }

  for(int iBar = 0; iBar < 16; ++iBar){
    if ( g_tRes_L_vs_th_fromTriangulation[iBar] -> GetN() == 0) continue;
    g_tRes_L_vs_th_fromTriangulation[iBar] -> Write(Form("g_tRes_L_vs_th_fromTriangulation_bar%02d",iBar));
    g_tRes_R_vs_th_fromTriangulation[iBar] -> Write(Form("g_tRes_R_vs_th_fromTriangulation_bar%02d",iBar));
    g_tRes_chRef_vs_th_fromTriangulation[iBar] -> Write(Form("g_tRes_chRef_vs_th_fromTriangulation_bar%02d",iBar));
  }

  fout->Close();
  
  
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);

  // ======== PLOT 
  //std::string plotDir(Form("/var/www/html/TOFHIR2B/MTDTB_CERN_June22/analysisSingleChannel/%s/",runs.c_str() ));
  //if (useTimeAverage) plotDir = Form("/var/www/html/TOFHIR2B/MTDTB_CERN_June22/analysisSingleChannel/%s_timeAverage/",runs.c_str() );
  //std::string plotDir(Form("/var/www/html/TOFHIR2X/MTDTB_CERN_June22/analysisSingleChannel/%s/",runs.c_str() ));
  //if (useTimeAverage) plotDir = Form("/var/www/html/TOFHIR2X/MTDTB_CERN_June22/analysisSingleChannel/%s_timeAverage/",runs.c_str() );
  std::string plotDir(Form("/var/www/html/TOFHIR2B/analysisSingleChannel/%s/",runs.c_str() ));
  if (useTimeAverage) plotDir = Form("/var/www/html/TOFHIR2B/analysisSingleChannel/%s_timeAverage/",runs.c_str() );

  system(Form("mkdir -p %s",plotDir.c_str()));  

  cout<< "Printing plots in "<< plotDir.c_str() <<endl;

  int th = 10;
  
  TCanvas *c;

  c = new TCanvas("c","c", 800, 600);
  h_nActiveChannels0[th] -> GetXaxis()->SetTitle("nActiveChannels");
  h_nActiveChannels0[th] -> Draw(); 
  TLine* line = new TLine(maxActiveChannels,0.,maxActiveChannels+0.5,h_nActiveChannels0[th] ->GetMaximum());
  line -> SetLineWidth(1);
  line -> SetLineStyle(7);
  line -> Draw("same");     
  c->Print(Form("%s/c_nActiveChannels_0.png",plotDir.c_str()));
  c->Print(Form("%s/c_nActiveChannels_0.pdf",plotDir.c_str()));
  delete c;

  c = new TCanvas("c","c", 800, 600);
  h_nActiveChannels1[th] -> GetXaxis()->SetTitle("nActiveChannels");
  h_nActiveChannels1[th] -> Draw();
  line -> Draw("same");           
  c->Print(Form("%s/c_nActiveChannels_1.png",plotDir.c_str()));
  c->Print(Form("%s/c_nActiveChannels_1.pdf",plotDir.c_str()));
  delete c;


  c = new TCanvas("c","c", 800, 600);
  //c->SetLogy();
  h_energy_chRef[th]-> GetXaxis()->SetTitle("energy [ADC]");
  h_energy_chRef[th]->Draw();      
  c->Print(Form("%s/c_energy_chRef.png",plotDir.c_str()));
  c->Print(Form("%s/c_energy_chRef.pdf",plotDir.c_str()));
  delete c;

  c = new TCanvas("c","c", 800, 600);
  //c->SetLogy();
  p_deltaT_LR_vs_energyRatio_chRef[th]-> GetXaxis()->SetTitle("energy ratio");
  p_deltaT_LR_vs_energyRatio_chRef[th]->GetXaxis()->SetRangeUser( p_deltaT_LR_vs_energyRatio_chRef[th]->GetMean()-5*p_deltaT_LR_vs_energyRatio_chRef[th]->GetRMS(), p_deltaT_LR_vs_energyRatio_chRef[th]->GetMean()+5*p_deltaT_LR_vs_energyRatio_chRef[th]->GetRMS());
  p_deltaT_LR_vs_energyRatio_chRef[th]->Draw();
  c->Print(Form("%s/c_deltaT_LR_vs_energyRatio_chRef.png",plotDir.c_str()));
  c->Print(Form("%s/c_deltaT_LR_vs_energyRatio_chRef.pdf",plotDir.c_str()));
  delete c;



  c = new TCanvas("c","c", 800, 600);
  h_deltaT_LR_energyRatioCorr_phaseCorr_chRef[th]-> GetXaxis()-> SetTitle("#Deltat [ps] ");
  TF1 *ff = new TF1("ff","gaus",-maxDeltaT,maxDeltaT);
  h_deltaT_LR_energyRatioCorr_phaseCorr_chRef[th]-> Fit(ff,"QRS");
  ff->SetRange( ff->GetParameter(1)-2*ff->GetParameter(2), ff->GetParameter(1)+2*ff->GetParameter(2));
  h_deltaT_LR_energyRatioCorr_phaseCorr_chRef[th]-> Fit(ff,"QRS");
  h_deltaT_LR_energyRatioCorr_phaseCorr_chRef[th]-> SetMarkerStyle(20);
  h_deltaT_LR_energyRatioCorr_phaseCorr_chRef[th]-> GetXaxis()-> SetRangeUser( ff->GetParameter(1)-7*ff->GetParameter(2),ff->GetParameter(1)+7*ff->GetParameter(2) );
  h_deltaT_LR_energyRatioCorr_phaseCorr_chRef[th]-> Draw("e");
  h_deltaT_LR_chRef[th]-> SetMarkerStyle(24);
  h_deltaT_LR_chRef[th]-> Draw("same");
  c->Print(Form("%s/c_deltaT_LR_chRef.png",plotDir.c_str()));
  c->Print(Form("%s/c_deltaT_LR_chRef.pdf",plotDir.c_str()));
  delete c;



  for(int iBar = 0; iBar < 16; ++iBar){
    
    if ( h_energy_LR[iBar][th]->GetEntries() == 0) continue;

    c = new TCanvas("c","c", 800, 600);
    h_energy_LR[iBar][th]-> GetXaxis()-> SetTitle("energy [ADC] ");
    h_energy_LR[iBar][th]-> SetMarkerStyle(20);
    h_energy_LR[iBar][th]-> Draw("");
    TLine* line = new TLine(energyMinLR[iBar][th],0.,energyMinLR[iBar][th],h_energy_LR[iBar][th]->GetMaximum());
    line -> SetLineWidth(1);
    line -> SetLineStyle(7);
    line -> Draw("same");   
    c->Print(Form("%s/c_energy_LR_bar%02d.png",plotDir.c_str(), iBar) );
    c->Print(Form("%s/c_energy_LR_bar%02d.pdf",plotDir.c_str(), iBar) );
    delete c;

    c = new TCanvas("c","c", 800, 600);
    h_deltaT_LR_energyRatioCorr_phaseCorr[iBar][th]-> GetXaxis()-> SetTitle("#Deltat [ps] ");
    h_deltaT_LR_energyRatioCorr_phaseCorr[iBar][th]-> SetMarkerStyle(20);
    h_deltaT_LR_energyRatioCorr_phaseCorr[iBar][th]-> Draw("e");
    c->Print(Form("%s/c_deltaT_LR_energyRatioCorr_phaseCorr_bar%02d.png",plotDir.c_str(), iBar) );
    c->Print(Form("%s/c_deltaT_LR_energyRatioCorr_phaseCorr_bar%02d.pdf",plotDir.c_str(), iBar) );
    delete c;


  }
  
  
  
  for (auto label : labelLR ){
    for(int iBar = 0; iBar < 16; ++iBar){        
      std::string chLabel = Form("bar%02d%s", iBar, label.c_str());

      if (h_energy[chLabel][th]->GetEntries() == 0) continue;

      if ( p_deltaT_vs_energyRatio[chLabel][th] == NULL) continue;

      // -- energy
      c = new TCanvas("c","c", 800, 600);
      //c->SetLogy();
      h_energy[chLabel][th]-> GetXaxis()->SetTitle("energy [ADC]");
      h_energy[chLabel][th]->Draw();
      /*TLine* line = new TLine(energyMin[chLabel],0.,energyMin[chLabel],h_energy[chLabel]->GetMaximum());
      line -> SetLineWidth(1);
      line -> SetLineStyle(7);
      line -> Draw("same");
      TLine* line2 = new TLine(maxEnergy,0.,maxEnergy,h_energy[chLabel]->GetMaximum());
      line2 -> SetLineWidth(1);
      line2 -> SetLineStyle(7);
      line2 -> Draw("same");
      */
      c->Print(Form("%s/c_energy_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_energy_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;
      
      // -- deltaT
      c = new TCanvas("c","c", 800, 600);
      h_deltaT[chLabel][th]-> GetXaxis()->SetTitle("#Deltat [ps]");
      h_deltaT[chLabel][th]->SetMarkerStyle(20);
      h_deltaT[chLabel][th]->Draw("e");
      c->Print(Form("%s/c_deltaT_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_deltaT_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;
      
      // -- tot Ratio
      c = new TCanvas("c","c", 800, 600);
      h_energyRatio[chLabel][th]-> GetXaxis()->SetTitle("E_{ch}/E_{chRef}");
      h_energyRatio[chLabel][th]->Draw();
      c->Print(Form("%s/c_energyRatio_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_energyRatio_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;
      
      // -- deltaT vs energyRatio 
      c = new TCanvas("c","c", 800, 600);
      //gStyle->SetOptFit(0);  
      p_deltaT_vs_energyRatio[chLabel][th]-> SetMarkerStyle(20);
      p_deltaT_vs_energyRatio[chLabel][th]-> SetMarkerSize(1);
      p_deltaT_vs_energyRatio[chLabel][th]-> GetXaxis()-> SetRangeUser( p_deltaT_vs_energyRatio[chLabel][th]->GetMean(1) - 7*p_deltaT_vs_energyRatio[chLabel][th]->GetRMS(1), p_deltaT_vs_energyRatio[chLabel][th]->GetMean(1) + 7*p_deltaT_vs_energyRatio[chLabel][th]->GetRMS(1));
      p_deltaT_vs_energyRatio[chLabel][th]-> GetYaxis()-> SetRangeUser(p_deltaT_vs_energyRatio[chLabel][th]->GetMean(2) - 3*p_deltaT_vs_energyRatio[chLabel][th]->GetRMS(2), p_deltaT_vs_energyRatio[chLabel][th]->GetMean(2) + 3*p_deltaT_vs_energyRatio[chLabel][th]->GetRMS(2));
      p_deltaT_vs_energyRatio[chLabel][th]-> GetXaxis()->SetTitle("E_{ch}/E_{chRef}");
      p_deltaT_vs_energyRatio[chLabel][th]-> GetYaxis()->SetTitle("#DeltaT [ps]");
      p_deltaT_vs_energyRatio[chLabel][th]->Draw();
      h2_deltaT_vs_energyRatio[chLabel][th]->Draw("colz same");
      p_deltaT_vs_energyRatio[chLabel][th]->Draw("same");
      fitFun_energyRatioCorr[chLabel][th]->Draw("same");
      c->Print(Form("%s/c_deltaT_vs_energyRatio_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_deltaT_vs_energyRatio_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;
      gStyle->SetOptFit(1111);  
      
      // -- deltaT corr
      c = new TCanvas("c","c", 800, 600);
      h_deltaT_energyRatioCorr[chLabel][th]-> GetXaxis()->SetTitle("#Deltat [ps]");
      h_deltaT_energyRatioCorr[chLabel][th]-> SetMarkerStyle(20);
      h_deltaT_energyRatioCorr[chLabel][th]-> Draw("e");
      c->Print(Form("%s/c_deltaT_energyRatioCorr_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_deltaT_energyRatioCorr_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;
      
      // -- deltaT energy+phase corr
      c = new TCanvas("c","c", 800, 600);
      h_deltaT_energyRatioCorr_phaseCorr[chLabel][th]-> GetXaxis()->SetTitle("#Deltat [ps]");
      h_deltaT_energyRatioCorr_phaseCorr[chLabel][th]-> SetMarkerStyle(20);
      h_deltaT_energyRatioCorr_phaseCorr[chLabel][th]-> Draw("e");
      c->Print(Form("%s/c_deltaT_energyRatioCorr_phaseCorr_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_deltaT_energyRatioCorr_phaseCorr_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;

      // -- deltaT corr vs t1fine
      c = new TCanvas("c","c", 800, 600);
      p_deltaT_energyRatioCorr_vs_t1fine[chLabel][th]->SetMarkerStyle(20);
      p_deltaT_energyRatioCorr_vs_t1fine[chLabel][th]->SetMarkerSize(1);
      p_deltaT_energyRatioCorr_vs_t1fine[chLabel][th]-> GetXaxis()->SetTitle("t1fine");
      p_deltaT_energyRatioCorr_vs_t1fine[chLabel][th]-> GetYaxis()->SetTitle("#Deltat [ps]");
      p_deltaT_energyRatioCorr_vs_t1fine[chLabel][th]-> GetYaxis()->SetRangeUser( h_deltaT_energyRatioCorr[chLabel][th]->GetMean() - 600, h_deltaT_energyRatioCorr[chLabel][th]->GetMean() + 600);
      p_deltaT_energyRatioCorr_vs_t1fine[chLabel][th]->Draw();
      h2_deltaT_energyRatioCorr_vs_t1fine[chLabel][th]->Draw("colz same");
      p_deltaT_energyRatioCorr_vs_t1fine[chLabel][th]->Draw("same");
      
      c->Print(Form("%s/c_deltaT_energyRatioCorr_vs_t1fine_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_deltaT_energyRatioCorr_vs_t1fine_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;
      
      
      // -- deltaT corr vs t1fine ref channel
      c = new TCanvas("c","c", 800, 600);
      p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel][th]->SetMarkerStyle(20);
      p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel][th]->SetMarkerSize(1);
      p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel][th]-> GetXaxis()->SetTitle("t1fine_{chRef}");
      p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel][th]-> GetYaxis()->SetRangeUser( h2_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel][th]->GetMean(2) - 600, h2_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel][th]->GetMean(2) + 600);
      p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel][th]->Draw();
      h2_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel][th]->Draw("colz same");
      p_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef[chLabel][th]->Draw("same");
      
      c->Print(Form("%s/c_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef_%s.png",plotDir.c_str(),chLabel.c_str()));
      c->Print(Form("%s/c_deltaT_energyRatioCorr_phaseCorr_vs_t1fineRef_%s.pdf",plotDir.c_str(),chLabel.c_str()));
      delete c;


      c = new TCanvas("c","c", 800, 600);
      g_tRes_energyRatioCorr_phaseCorr_vs_th_final[label][iBar]-> GetXaxis()->SetTitle("ith1 [DAC]");
      g_tRes_energyRatioCorr_phaseCorr_vs_th_final[label][iBar]-> GetYaxis()->SetTitle("#sigma_{t} [ps]");
      g_tRes_energyRatioCorr_phaseCorr_vs_th_final[label][iBar]-> Draw("ap");
      c->Print(Form("%s/c_tRes_energyRatioCorr_phaseCorr_vs_th_%s.png",plotDir.c_str(), chLabel.c_str()) );
    }
  }
 
  gStyle->SetOptStat(0);
  
  // -- no corr
  float ymin = 30;
  float ymax = 200;
  if (step1 == 1.50) ymax = 240;

  c = new TCanvas("c","c", 800, 600);
  c->SetGridx();
  c->SetGridy();
  TH2F* hdummy1 = new TH2F("hdummy1","",16,-0.5,15.5,100,ymin,ymax);
  hdummy1->GetXaxis()->SetTitle("bar");
  hdummy1->GetYaxis()->SetTitle("#sigma(t_{ch} - t_{ref}) [ps]");
  hdummy1->Draw();
  TLegend *leg0 = new TLegend(0.20, 0.20, 0.35, 0.35);
  leg0->SetBorderSize(0);
  leg0->SetFillStyle(0);
  for (auto label : labelLR ){   
    g_tRes_vs_bar[label][th]->SetMarkerStyle(20);
    if (label == "R") g_tRes_vs_bar[label][th]->SetMarkerStyle(24);
    g_tRes_vs_bar[label][th]->SetMarkerSize(1);
    g_tRes_vs_bar[label][th]->Draw("psame") ;
    g_tRes_energyRatioCorr_vs_bar[label][th]->SetMarkerColor(2);
    g_tRes_energyRatioCorr_vs_bar[label][th]->SetLineColor(2);
    g_tRes_energyRatioCorr_vs_bar[label][th]->SetMarkerStyle(20);
    if (label == "R") g_tRes_energyRatioCorr_vs_bar[label][th]->SetMarkerStyle(24);
    g_tRes_energyRatioCorr_vs_bar[label][th]->SetMarkerSize(1);
    g_tRes_energyRatioCorr_vs_bar[label][th]->Draw("psame") ;
    leg0->AddEntry(g_tRes_vs_bar[label][th], label.c_str(), "PL");
  }
  leg0->Draw("same");
  c->Print(Form("%s/c_tRes_vs_bar.png",plotDir.c_str()));                                                                                                     
  c->Print(Form("%s/c_tRes_vs_bar.pdf",plotDir.c_str()));                                                                                                     
  delete c;
  
  
  // -- energy corr
  c = new TCanvas("c","c", 800, 600);
  c->SetGridx();
  c->SetGridy();
  TH2F* hdummy2 = new TH2F("hdummy2","",16,-0.5,15.5,100,ymin,ymax);
  hdummy2->GetXaxis()->SetTitle("bar");
  hdummy2->GetYaxis()->SetTitle("#sigma(t_{ch} - t_{ref}) [ps]");
  hdummy2->Draw();
  for (auto label : labelLR ){   
    g_tRes_energyRatioCorr_vs_bar[label][th]->SetMarkerColor(2);
    g_tRes_energyRatioCorr_vs_bar[label][th]->SetLineColor(2);
    g_tRes_energyRatioCorr_vs_bar[label][th]->SetMarkerStyle(20);
    if (label == "R") g_tRes_energyRatioCorr_vs_bar[label][th]->SetMarkerStyle(24);
    g_tRes_energyRatioCorr_vs_bar[label][th]->SetMarkerSize(1);
    g_tRes_energyRatioCorr_vs_bar[label][th]->Draw("psame") ;
  }
  c->Print(Form("%s/c_tRes_energyRatioCorr_vs_bar.png",plotDir.c_str()));
  c->Print(Form("%s/c_tRes_energyRatioCorr_vs_bar.pdf",plotDir.c_str()));
  delete c;
  
    
  // -- energy+phase corr
  c = new TCanvas("c","c", 800, 600);
  c->SetGridx();
  c->SetGridy();
  TH2F* hdummy3 = new TH2F("hdummy3","",16,-0.5,15.5,100,ymin,ymax);
  hdummy3->GetXaxis()->SetTitle("bar");
  hdummy3->GetYaxis()->SetTitle("#sigma(t_{ch} - t_{ref}) [ps]");
  hdummy3->Draw();
  for (auto label : labelLR ){   
    g_tRes_energyRatioCorr_phaseCorr_vs_bar[label][th]->SetMarkerColor(4);
    g_tRes_energyRatioCorr_phaseCorr_vs_bar[label][th]->SetLineColor(4);
    g_tRes_energyRatioCorr_phaseCorr_vs_bar[label][th]->SetMarkerStyle(20);
    if (label == "R") g_tRes_energyRatioCorr_phaseCorr_vs_bar[label][th]->SetMarkerStyle(24);
    g_tRes_energyRatioCorr_phaseCorr_vs_bar[label][th]->SetMarkerSize(1);
    g_tRes_energyRatioCorr_phaseCorr_vs_bar[label][th]->Draw("psame") ;
  }
  c->Print(Form("%s/c_tRes_energyRatioCorr_phaseCorr_vs_bar.png",plotDir.c_str()));
  c->Print(Form("%s/c_tRes_energyRatioCorr_phaseCorr_vs_bar.pdf",plotDir.c_str()));
  delete c;

  // -- tRes of L, R, chRef after triangulation
  c = new TCanvas("c","c", 800, 600);
  c->SetGridx();
  c->SetGridy();
  ymin = 0;
  ymax = 120;
  if ( step1 == 1.50){
    ymax = 160;
  }
  TH2F* hdummy4 = new TH2F("hdummy4","",16,-0.5,15.5,100,ymin,ymax);
  hdummy4->GetXaxis()->SetTitle("bar");
  hdummy4->GetYaxis()->SetTitle("#sigma_{t} [ps]");
  hdummy4->Draw();
  g_tRes_L_vs_bar_fromTriangulation[th] ->SetMarkerStyle(20);
  g_tRes_L_vs_bar_fromTriangulation[th] -> SetLineColor(51);
  g_tRes_L_vs_bar_fromTriangulation[th] -> SetMarkerColor(51);
  g_tRes_R_vs_bar_fromTriangulation[th] -> SetMarkerStyle(24);
  g_tRes_R_vs_bar_fromTriangulation[th] -> SetLineColor(51);
  g_tRes_R_vs_bar_fromTriangulation[th] -> SetMarkerColor(51);
  g_tRes_chRef_vs_bar_fromTriangulation[th] ->SetMarkerStyle(21);
  g_tRes_chRef_vs_bar_fromTriangulation[th] ->SetMarkerColor(12);
  g_tRes_chRef_vs_bar_fromTriangulation[th] ->SetLineColor(12);
  g_tRes_L_vs_bar_fromTriangulation[th]->Draw("psame") ;
  g_tRes_R_vs_bar_fromTriangulation[th]->Draw("psame") ;
  g_tRes_chRef_vs_bar_fromTriangulation[th]->Draw("psame") ;
  
  TLegend *leg = new TLegend(0.20, 0.20, 0.35, 0.35);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(g_tRes_L_vs_bar_fromTriangulation[th], "L", "P");
  leg->AddEntry(g_tRes_R_vs_bar_fromTriangulation[th], "R", "P");
  if (useTimeAverage) leg->AddEntry(g_tRes_chRef_vs_bar_fromTriangulation[th], "Ref (average)", "P");
  else leg->AddEntry(g_tRes_chRef_vs_bar_fromTriangulation[th], "Ref (single)", "P");
  leg->Draw("same");


  TF1 *myfitL = new TF1("myfitL","pol0", 0, 100);
  myfitL -> SetLineColor(51);
  myfitL -> SetParameter(0, g_tRes_L_vs_bar_fromTriangulation[th] -> GetMean(2));
  g_tRes_L_vs_bar_fromTriangulation[th] ->Fit("myfitL");
  TLatex *latexL = new TLatex(0.50,0.32,Form("<#sigma_{t}> = %.1f #pm %.1f ps", myfitL->GetParameter(0), myfitL->GetParError(0)));
  latexL -> SetNDC();
  latexL -> SetTextFont(42);
  latexL -> SetTextSize(0.04);
  latexL->Draw("same");
  TF1 *myfitR = new TF1("myfitR","pol0", 0, 100);
  myfitR -> SetLineColor(51);
  myfitR -> SetLineStyle(2);
  myfitR -> SetParameter(0, g_tRes_R_vs_bar_fromTriangulation[th] -> GetMean(2));
  g_tRes_R_vs_bar_fromTriangulation[th] ->Fit("myfitR");
  TLatex *latexR = new TLatex(0.50,0.27,Form("<#sigma_{t}> = %.1f #pm %.1f ps", myfitR->GetParameter(0), myfitR->GetParError(0)));
  latexR -> SetNDC();
  latexR -> SetTextFont(42);
  latexR -> SetTextSize(0.04);
  latexR->Draw("same");
  TF1 *myfitChRef = new TF1("myfitChRef","pol0", 0, 100);
  myfitChRef -> SetLineColor(12);
  myfitChRef -> SetParameter(0,g_tRes_chRef_vs_bar_fromTriangulation[th] -> GetMean(2));
  g_tRes_chRef_vs_bar_fromTriangulation[th] ->Fit("myfitChRef");
  TLatex *latexRef = new TLatex(0.50,0.22,Form("<#sigma_{t}> = %.1f #pm %.1f ps", myfitChRef->GetParameter(0), myfitChRef->GetParError(0)));
  latexRef -> SetNDC();
  latexRef -> SetTextFont(42);
  latexRef -> SetTextSize(0.04);
  latexRef->Draw("same");
  
  std::cout << "L    : " <<  "   <sigma_t> = " << myfitL->GetParameter(0) << "+/-" << myfitL->GetParError(0) <<std::endl;
  std::cout << "R    : " <<  "   <sigma_t> = " << myfitR->GetParameter(0) << "+/-" << myfitR->GetParError(0) <<std::endl;
  std::cout << "ChRef: " <<  "   <sigma_t> = " << myfitChRef->GetParameter(0) << "+/-" << myfitChRef->GetParError(0) <<std::endl;
    
  c->Print(Form("%s/c_tRes.png",plotDir.c_str()));
  c->Print(Form("%s/c_tRes.pdf",plotDir.c_str()));
  delete c;

}
  
