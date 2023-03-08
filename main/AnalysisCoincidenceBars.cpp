#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
 
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
//#include "TPaveStats.h"
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
#include "TLine.h"
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



int main(int argc, char** argv){
  
  map<int, std::vector<int>> bars; 
  bars[0].push_back(8);
  bars[0].push_back(9);
  bars[0].push_back(10);

  bars[1].push_back(8);
  bars[1].push_back(9);
  bars[1].push_back(10);

  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);

  int run  = opts.GetOpt<int>("Inputs.run");
  int pedestals  = opts.GetOpt<int>("Inputs.pedestals");

  //  float enMin = opts.GetOpt<float>("Plots.energyMin");
  //float Vov = opts.GetOpt<float>("Plots.Vov");
  std::vector<float> Vovs = opts.GetOpt<std::vector<float> >("Plots.Vov");
  std::vector<float> energyMins = opts.GetOpt<std::vector<float> >("Plots.energyMin");
  map<float, float> enMins;
  for (unsigned int ivov = 0; ivov < Vovs.size(); ivov++){
    enMins[Vovs[ivov]] = energyMins[ivov];
  }
  

  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  system(Form("mkdir -p %s",plotDir.c_str()));


  TChain* tree = new TChain("data","data");
  if (pedestals)
    tree->Add(Form("/data/tofhir2/reco/run%d_ped_t2.root", run));
  else
    tree->Add(Form("/data/tofhir2/reco_v2/run%d_t2.root", run));

  //--- define branches
  float step1, step2;
  int channelIdx[128];
  std::vector<float> *qfine = 0;
  std::vector<float> *tot = 0;
  std::vector<float> *energy = 0;
  std::vector<long long> *time = 0;
  
  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("step1",  1); tree -> SetBranchAddress("step1",  &step1);
  tree -> SetBranchStatus("step2",  1); tree -> SetBranchAddress("step2",  &step2);
  tree -> SetBranchStatus("channelIdx",  1); tree -> SetBranchAddress("channelIdx",  channelIdx);
  tree -> SetBranchStatus("qfine",  1); tree -> SetBranchAddress("qfine",   &qfine);  
  tree -> SetBranchStatus("tot",    1); tree -> SetBranchAddress("tot",       &tot);
  tree -> SetBranchStatus("energy", 1); tree -> SetBranchAddress("energy", &energy);
  tree -> SetBranchStatus("time",   1); tree -> SetBranchAddress("time",     &time);

  int nEntries = tree->GetEntries();
  cout << "Number of entries = " << nEntries << endl;
  //int maxEntries = 12000000;
  int maxEntries = nEntries;
  
  // -- channel mapping
  //--- define channels (read mapping from the configuration file)
  std::vector<unsigned int> channelMapping = opts.GetOpt<std::vector<unsigned int> >("Channels.channelMapping");
  
  int chL[2][16];
  int chR[2][16];
  
  for(int iBar = 0; iBar < 16; ++iBar){
      chL[0][iBar] = channelMapping[iBar*2+0];
      chR[0][iBar] = channelMapping[iBar*2+1];
      chL[1][iBar] = channelMapping[iBar*2+0]+64;
      chR[1][iBar] = channelMapping[iBar*2+1]+64;
  }
      
  std::cout<< "Analyzing channels: " <<std::endl;
  for (int iMod = 0; iMod < 2; iMod++){
    for (int iBar = 0; iBar < 16; iBar++){ 
      bool found = std::find(std::begin(bars[iMod]), std::end(bars[iMod]), iBar) != std::end(bars[iMod]);
      if (!found) continue;
      std::cout<<iMod << "  bar = "<< iBar<< "   chL = " << chL[iMod][iBar] << "   chR = "<<chR[iMod][iBar]<<std::endl;
    }
  }

  // -- book histograms 

  // histograms for single bars
  map<float, map<int, map<int,TH1F*> > > h_energyL;
  map<float, map<int, map<int,TH1F*> > > h_energyR;
  map<float, map<int, map<int,TH1F*> > > h_energyLR;

  map<float, map<int, map<int,TH1F*> > > h_energyRatio;
  map<float, map<int, map<int,TProfile*> > > p_dtLR_vs_energyRatio;
  map<float, map<int, map<int,TH2F*> > > h2_dtLR_vs_energyRatio;
  map<float, map<int, map<int,TH1F*> > > h_dtLR;
  map<float, map<int, map<int,TH1F*> > > h_dtLR_energyRatioCorr;


  
  for (int iMod = 0; iMod < 2; iMod++){
    for (int iBar = 0; iBar < 16; iBar++){
  
      bool found = std::find(std::begin(bars[iMod]), std::end(bars[iMod]), iBar) != std::end(bars[iMod]);
      if (!found) continue;    
      
      for (unsigned int ivov = 0; ivov < Vovs.size(); ivov++){
	float Vov = Vovs[ivov];
	h_energyL[Vov][iMod][iBar] = new TH1F(Form("h_energyL_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), Form("h_energyL_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), 1000, 0, 1000);
	h_energyR[Vov][iMod][iBar] = new TH1F(Form("h_energyR_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), Form("h_energyR_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), 1000, 0, 1000);
	h_energyLR[Vov][iMod][iBar] = new TH1F(Form("h_energyLR_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), Form("h_energyLR_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), 1000, 0, 1000);	
	h_energyRatio[Vov][iMod][iBar] = new TH1F(Form("h_energyRatio_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), Form("h_energyRatio_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), 200, 0, 2);
	p_dtLR_vs_energyRatio[Vov][iMod][iBar] = new TProfile(Form("p_dtLR_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), Form("p_dtLR_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), 200, 0, 2, -5000, 5000);
	h2_dtLR_vs_energyRatio[Vov][iMod][iBar] = new TH2F(Form("h2_dtLR_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), Form("h2_dtLR_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), 200, 0, 2, 1000, -5000, 5000);
	h_dtLR[Vov][iMod][iBar] = new TH1F(Form("h_dtLR_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), Form("h_dtLR_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), 1000, -5000, 5000);
	h_dtLR_energyRatioCorr[Vov][iMod][iBar] = new TH1F(Form("h_dtLR_energyRatioCorr_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), Form("h_dtLR_energyRatioCorr_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), 1000, -5000, 5000);
      }
    }
  }

  // compare two bars
  map<float, map<int, map<int, TProfile*> > > p_dtAve_vs_energyRatio; 
  map<float, map<int, map<int, TH2F*> > > h2_dtAve_vs_energyRatio;
  map<float, map<int, map<int, TH1F*> > > h_energyAveRatio;
  map<float, map<int, map<int, TH1F*> > > h_dtAve;
  map<float, map<int, map<int, TH1F*> > > h_dtAve_energyRatioCorr;
  map<float, map<int, map<int, TH1F*> > > h_tmp;

  for (int iBar0 = 0; iBar0 < 16; iBar0++){
    bool found0 = std::find(std::begin(bars[0]), std::end(bars[0]), iBar0) != std::end(bars[0]);    
    if (!found0) continue;
    for (int iBar1 = 0; iBar1 < 16; iBar1++){   
      bool found1 = std::find(std::begin(bars[1]), std::end(bars[1]), iBar1) != std::end(bars[1]);    
      if (!found1) continue;
      for (unsigned int ivov = 0; ivov < Vovs.size(); ivov++){
        float Vov = Vovs[ivov];  
	p_dtAve_vs_energyRatio[Vov][iBar0][iBar1] = new TProfile(Form("p_dtAve_vs_energyRatio_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov), Form("p_dtAve_vs_energyRatio_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov), 200, 0, 2, -5000, 5000);
	h2_dtAve_vs_energyRatio[Vov][iBar0][iBar1] = new TH2F(Form("h2_dtAve_vs_energyRatio_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov), Form("h2_dtAve_vs_energyRatio_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov), 200, 0, 2, 1000,-5000, 5000);
	h_energyAveRatio[Vov][iBar0][iBar1] = new TH1F(Form("h_energyAveRatio_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov), Form("h_energyAveRatio_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov), 400, 0, 2);
	h_dtAve[Vov][iBar0][iBar1] = new TH1F(Form("h_dtAve_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov),Form("h_dtAve_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov), 1000, -5000, 5000);
	h_dtAve_energyRatioCorr[Vov][iBar0][iBar1] = new TH1F(Form("h_dtAve_energyRatioCorr_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov), Form("h_dtAve_energyRatioCorr_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov), 1000, -5000, 5000 );
	h_tmp[Vov][iBar0][iBar1] = new TH1F(Form("h_tmp_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov), Form("h_tmp_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov), 100000, -2000000, 2000000);
      }
    }
  }
  

  //acceptEvent[Vov][Imod][ibar][entry]
  map<float, map<int, map <int, map <int, bool> > > > acceptEvent;

  // -- first loop over events
  cout << "First loop over events to find the mip peak" <<endl;
  for (int entry = 0; entry < maxEntries; entry++){
    
    tree->GetEntry(entry);

    if( entry%1000 == 0 ) std::cout << ">>> Reading entry " << entry << " / " << nEntries << "\r" << std::flush;

    if (step2 != 110406 ) continue;

    float Vov = step1;
    bool goodVov = std::find(std::begin(Vovs), std::end(Vovs), Vov) != std::end(Vovs);   
    if (!goodVov) continue;
    
    
    float aveEnergy[2][16];  
    long long aveTime[2][16];   
      
    for (int iMod = 0; iMod < 2; iMod++){
      for (int iBar = 0; iBar < 16; iBar++){
	
	aveEnergy[iMod][iBar] = -999;  
	aveTime[iMod][iBar] = -999;   
	
	bool found = std::find(std::begin(bars[iMod]), std::end(bars[iMod]), iBar) != std::end(bars[iMod]);
	if (!found) continue;          
	
	if (channelIdx[chL[iMod][iBar]] < 0 || channelIdx[chR[iMod][iBar]] < 0);
	
	float totL =(*tot)[channelIdx[chL[iMod][iBar]]];
	float totR =(*tot)[channelIdx[chR[iMod][iBar]]];
	
	if (totL < 0 ) continue;
	if (totR < 0 ) continue;
	
	float enL =(*energy)[channelIdx[chL[iMod][iBar]]];
	float enR =(*energy)[channelIdx[chR[iMod][iBar]]];
	
	h_energyL[Vov][iMod][iBar] -> Fill(enL);	
	h_energyR[Vov][iMod][iBar]  -> Fill(enR);	
	h_energyLR[Vov][iMod][iBar]  -> Fill(0.5*(enL+enR));	
	
	if (0.5 * (enL+enR) < enMins[Vov]) continue;                                                                                                                
	long long tL = (*time)[channelIdx[chL[iMod][iBar]]];
	long long tR = (*time)[channelIdx[chR[iMod][iBar]]];
	aveEnergy[iMod][iBar] = 0.5 * (enL+enR) ;  
	aveTime[iMod][iBar] = 0.5 * (tL+tR) ;
	h_energyRatio[Vov][iMod][iBar]->Fill(enL/enR);  
      }
    }
    
    
    for (int iBar0 = 0; iBar0 < 16; iBar0++){
      bool found0 = std::find(std::begin(bars[0]), std::end(bars[0]), iBar0) != std::end(bars[0]);
      if (!found0) continue;
      for (int iBar1 = 0; iBar1 < 16; iBar1++){
	bool found1 = std::find(std::begin(bars[1]), std::end(bars[1]), iBar1) != std::end(bars[1]);
	if (!found1) continue;
	
	if ( aveEnergy[0][iBar0] > enMins[Vov] &&  aveEnergy[1][iBar1] > enMins[Vov]){
	  h_tmp[Vov][iBar0][iBar1] -> Fill(aveTime[0][iBar0]  - aveTime[1][iBar1] );
	  float eRatio = aveEnergy[0][iBar0]/aveEnergy[1][iBar1] ; 
	  h_energyAveRatio[Vov][iBar0][iBar1]->Fill( eRatio );
	}
      }
    }
  }// -- end first loop over entries


  // offset[Vov][iBar1][iBar2]
  map<float, map<int, map<int,long long > > > offset;
  for (int iBar0 = 0; iBar0 < 16; iBar0++){                                                                                                                                 
    bool found0 = std::find(std::begin(bars[0]), std::end(bars[0]), iBar0) != std::end(bars[0]);                                                                            
    if (!found0) continue;                                                                                                                                                  
    for (int iBar1 = 0; iBar1 < 16; iBar1++){                                                                                                                               
      bool found1 = std::find(std::begin(bars[1]), std::end(bars[1]), iBar1) != std::end(bars[1]);                                                                          
      if (!found1) continue;              
      for (unsigned int ivov =0 ; ivov < Vovs.size(); ivov++){
	float Vov = Vovs[ivov];
	offset[Vov][iBar0][iBar1] = h_tmp[Vov][iBar0][iBar1] ->GetBinCenter(h_tmp[Vov][iBar0][iBar1]->GetMaximumBin());
      }
    }
  }
  

  // -- find min energy
  map<float, map<int, map<int, float> > > energyMin;
  
  for (int iMod = 0; iMod < 2; iMod++){
    for (int iBar = 0; iBar < 16; iBar++){
      
      bool found = std::find(std::begin(bars[iMod]), std::end(bars[iMod]), iBar) != std::end(bars[iMod]);
      if (!found) continue;          
      
      for (unsigned int ivov =0 ; ivov < Vovs.size(); ivov++){
        float Vov = Vovs[ivov];       
	//TF1 *fun = new TF1("fun","landau", 0, 1000);
	TF1 *fun = new TF1("fun","gaus", 0, 1000);
	fun->SetLineColor(1);
	h_energyLR[Vov][iMod][iBar]->GetXaxis()->SetRangeUser(enMins[Vov],1000);
	int maxbin = h_energyLR[Vov][iMod][iBar]->GetMaximumBin();
	float peak = h_energyLR[Vov][iMod][iBar]->GetBinCenter(maxbin);
	fun->SetRange(peak*0.9, peak*1.1);
	fun->SetParameter(1,peak);
	fun->SetParameter(2,0.1*peak);
	h_energyLR[Vov][iMod][iBar]->Fit("fun","QR");
	energyMin[Vov][iMod][iBar] = fun->GetParameter(1)*0.8;
	if (energyMin[Vov][iMod][iBar] < enMins[Vov]) energyMin[Vov][iMod][iBar] = enMins[Vov];
	delete fun;
	cout << "Vov = " << Vov << "   MOD "<<iMod<<"   BAR "<< iBar << "   Energy min = " << energyMin[Vov][iMod][iBar]<<endl;
	h_energyLR[Vov][iMod][iBar]->GetXaxis()->SetRangeUser(0,1000);
      }
    }
  }



  // -- energy ratio range
  map<float, map<int, map<int, float> > > enRatioMin;  
  map<float, map<int, map<int, float> > > enRatioMax;  

  for (int iMod = 0; iMod < 2; iMod++){
    for (int iBar=0; iBar < 16; iBar++){

      bool found = std::find(std::begin(bars[iMod]), std::end(bars[iMod]), iBar) != std::end(bars[iMod]);
      if (!found) continue;          

      for (unsigned int ivov =0 ; ivov < Vovs.size(); ivov++){
	float Vov = Vovs[ivov];     
	TF1 *fitfun = new TF1("fitfun", "gaus", 0, 2);
	fitfun->SetRange(h_energyRatio[Vov][iMod][iBar]->GetMean()-3*h_energyRatio[Vov][iMod][iBar]->GetRMS(), 
			 h_energyRatio[Vov][iMod][iBar]->GetMean()+3*h_energyRatio[Vov][iMod][iBar]->GetRMS());
	h_energyRatio[Vov][iMod][iBar]-> Fit(fitfun, "QRS");
	enRatioMin[Vov][iMod][iBar] = fitfun->GetParameter(1) - 3.5*fitfun->GetParameter(2);
	enRatioMax[Vov][iMod][iBar] = fitfun->GetParameter(1) + 3.5*fitfun->GetParameter(2);
	delete fitfun;
      }
    }
  }

  map<float, map<int, map<int, float> > > enRatioAveMin;  
  map<float, map<int, map<int, float> > > enRatioAveMax;  
  for (int iBar0 = 0; iBar0 < 16; iBar0++){
    bool found0 = std::find(std::begin(bars[0]), std::end(bars[0]), iBar0) != std::end(bars[0]);
    if (!found0) continue;
    for (int iBar1 = 0; iBar1 < 16; iBar1++){
      bool found1 = std::find(std::begin(bars[1]), std::end(bars[1]), iBar1) != std::end(bars[1]);
      if (!found1) continue; 
      
      for (unsigned int ivov =0 ; ivov < Vovs.size(); ivov++){
	float Vov = Vovs[ivov]; 
	TF1 *fitfun = new TF1("fitfun", "gaus", 0, 2);
	fitfun->SetRange(h_energyAveRatio[Vov][iBar0][iBar1]->GetMean()-3*h_energyAveRatio[Vov][iBar0][iBar1]->GetRMS(), 
			 h_energyAveRatio[Vov][iBar0][iBar1]->GetMean()+3*h_energyAveRatio[Vov][iBar0][iBar1]->GetRMS());
	h_energyAveRatio[Vov][iBar0][iBar1]-> Fit(fitfun, "QRS");
	enRatioAveMin[Vov][iBar0][iBar1] = fitfun->GetParameter(1) - 3.5*fitfun->GetParameter(2);
	enRatioAveMax[Vov][iBar0][iBar1] = fitfun->GetParameter(1) + 3.5*fitfun->GetParameter(2); 
      }
    }
  }


  // -- second loop over events to get amp walk corrections
  cout << "Second loop over events to get amp walk corrections" <<endl;
  for (int entry = 0; entry < maxEntries; entry++){
    
    if( entry%1000 == 0 ) std::cout << ">>> Reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    tree->GetEntry(entry);

    if (step2 != 110406 ) continue;                                                                                                                                           

    float Vov = step1;
    bool goodVov = std::find(std::begin(Vovs), std::end(Vovs), Vov) != std::end(Vovs);
    if (!goodVov) continue; 

    
    float aveEnergy[2][16];
    long long aveTime[2][16];
    
    // -- single bars
    for (int iMod = 0; iMod < 2; iMod++){
      for (int iBar=0; iBar < 16; iBar++){
	
	acceptEvent[Vov][iMod][iBar][entry] = false;
	  
	aveEnergy[iMod][iBar] = -999;
	aveTime[iMod][iBar] = -999;
	
	bool found = std::find(std::begin(bars[iMod]), std::end(bars[iMod]), iBar) != std::end(bars[iMod]);
	if (!found) continue;          
	
	if (channelIdx[chL[iMod][iBar]] < 0 || channelIdx[chR[iMod][iBar]] < 0);
	
	float totL =(*tot)[channelIdx[chL[iMod][iBar]]];
	float totR =(*tot)[channelIdx[chR[iMod][iBar]]];
	
	if (totL < 0 ) continue;
	if (totR < 0 ) continue;
	
	float enL =(*energy)[channelIdx[chL[iMod][iBar]]];
	float enR =(*energy)[channelIdx[chR[iMod][iBar]]];
	
	aveEnergy[iMod][iBar] = 0.5 * (enL+enR) ;
	
	if (0.5 * (enL+enR) < energyMin[Vov][iMod][iBar]) continue;
	if (  (enL/enR) < enRatioMin[Vov][iMod][iBar] || (enL/enR) > enRatioMax[Vov][iMod][iBar] ) continue;
	
	acceptEvent[Vov][iMod][iBar][entry] = true;	
	
	long long tL = (*time)[channelIdx[chL[iMod][iBar]]];
	long long tR = (*time)[channelIdx[chR[iMod][iBar]]];
	aveTime[iMod][iBar] = 0.5 * (tL+tR) ;
	
	float dt = tL - tR;
	h_dtLR[Vov][iMod][iBar]->Fill(dt);
	p_dtLR_vs_energyRatio[Vov][iMod][iBar]->Fill(enL/enR,dt);
	h2_dtLR_vs_energyRatio[Vov][iMod][iBar]->Fill(enL/enR,dt);
      }
    }// end loop over modules
    
    // -- 2 bars
    for (int iBar0 = 0; iBar0 < 16; iBar0++){
      bool found0 = std::find(std::begin(bars[0]), std::end(bars[0]), iBar0) != std::end(bars[0]);
      if (!found0) continue;
      for (int iBar1 = 0; iBar1 < 16; iBar1++){
	bool found1 = std::find(std::begin(bars[1]), std::end(bars[1]), iBar1) != std::end(bars[1]);
	if (!found1) continue;
	
	if( !acceptEvent[Vov][0][iBar0][entry] || !acceptEvent[Vov][1][iBar1][entry]) continue;
	  
	if ( aveEnergy[0][iBar0] > energyMin[0][iBar0][Vov] &&  aveEnergy[1][iBar1] > energyMin[1][iBar1][Vov]){
	  float eRatio = aveEnergy[0][iBar0]/aveEnergy[1][iBar1] ;
	  float dt = aveTime[0][iBar0]  - aveTime[1][iBar1] ;
	  dt = dt - offset[Vov][iBar0][iBar1];
	  h_dtAve[Vov][iBar0][iBar1]->Fill(dt);
	  p_dtAve_vs_energyRatio[Vov][iBar0][iBar1] ->Fill(eRatio, dt);
	    h2_dtAve_vs_energyRatio[Vov][iBar0][iBar1] ->Fill(eRatio, dt);
	}
      }
    }
  }// -- end second loop over entries
  
  
  // ---  amp walk corr
  map<float, map<int, map <int,TF1*> > >fitFunCorr;

  for (int iMod = 0; iMod < 2; iMod++){
    for (int iBar=0; iBar < 16; iBar++){

      bool found = std::find(std::begin(bars[iMod]), std::end(bars[iMod]), iBar) != std::end(bars[iMod]);                                                                   
      if (!found) continue;  

      for (unsigned int ivov = 0; ivov < Vovs.size(); ivov++){
	float Vov = Vovs[ivov];
	fitFunCorr[Vov][iMod][iBar] = new TF1(Form("fitFunCorr_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov), "pol4",0, 2 );
	fitFunCorr[Vov][iMod][iBar]->SetRange(enRatioMin[Vov][iMod][iBar], enRatioMax[Vov][iMod][iBar]);
	p_dtLR_vs_energyRatio[Vov][iMod][iBar]-> Fit(fitFunCorr[Vov][iMod][iBar],"QRS");
      }
    }
  }
  

  map<float, map<int, map <int,TF1*> > > fitFunCorr2;
  for (int iBar0 = 0; iBar0 < 16; iBar0++){                                                                                                                                 
    bool found0 = std::find(std::begin(bars[0]), std::end(bars[0]), iBar0) != std::end(bars[0]);                                                                            
    if (!found0) continue;                                                                                                                                                  
    for (int iBar1 = 0; iBar1 < 16; iBar1++){                                                                                                                               
      bool found1 = std::find(std::begin(bars[1]), std::end(bars[1]), iBar1) != std::end(bars[1]);                                                                          
      if (!found1) continue;  
      for (unsigned int ivov = 0; ivov < Vovs.size(); ivov++){
        float Vov = Vovs[ivov]; 
	fitFunCorr2[Vov][iBar0][iBar1] = new TF1(Form("fitFunCorr2_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov), "pol4",0, 2 );
	fitFunCorr2[Vov][iBar0][iBar1]->SetRange(enRatioAveMin[Vov][iBar0][iBar1], enRatioAveMax[Vov][iBar0][iBar1]);
	p_dtAve_vs_energyRatio[Vov][iBar0][iBar1]-> Fit(fitFunCorr2[Vov][iBar0][iBar1], "QRS");
      }
    }
  }


  // -- third loop over events to apply amp walk corrections
  cout << "Third loop over events to apply amp walk corrections" <<endl;
  for (int entry = 0; entry < maxEntries; entry++){
    
    if( entry%1000 == 0 ) std::cout << ">>> Reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    tree->GetEntry(entry);

    if (step2 != 110406 ) continue;  
    
    float Vov = step1;
    bool goodVov = std::find(std::begin(Vovs), std::end(Vovs), Vov) != std::end(Vovs);
    if (!goodVov) continue;                                                                                                                                       


    float aveEnergy[2][16];
    long long aveTime[2][16];
    
    // -- single bars
    for (int iMod = 0; iMod < 2; iMod++){
      for (int iBar=0; iBar < 16; iBar++){

	if( !acceptEvent[Vov][iMod][iBar][entry]) continue;

	aveEnergy[iMod][iBar] = -999;
	aveTime[iMod][iBar] = -999;

	bool found = std::find(std::begin(bars[iMod]), std::end(bars[iMod]), iBar) != std::end(bars[iMod]);
        if (!found) continue;   

        float enL =(*energy)[channelIdx[chL[iMod][iBar]]];
        float enR =(*energy)[channelIdx[chR[iMod][iBar]]];

        aveEnergy[iMod][iBar] = 0.5 * (enL+enR) ;

        long long tL = (*time)[channelIdx[chL[iMod][iBar]]];
        long long tR = (*time)[channelIdx[chR[iMod][iBar]]];
        aveTime[iMod][iBar] = 0.5 * (tL+tR) ;

        float dt = tL - tR;
	float dtCorr = dt - fitFunCorr[Vov][iMod][iBar]->Eval(enL/enR) + fitFunCorr[Vov][iMod][iBar]->Eval( h_energyRatio[Vov][iMod][iBar]->GetMean());
        h_dtLR_energyRatioCorr[Vov][iMod][iBar]->Fill(dtCorr);
      }
    }
    

    // -- 2 bars
    for (int iBar0 = 0; iBar0 < 16; iBar0++){
      bool found0 = std::find(std::begin(bars[0]), std::end(bars[0]), iBar0) != std::end(bars[0]);
      if (!found0) continue;
      for (int iBar1 = 0; iBar1 < 16; iBar1++){
        bool found1 = std::find(std::begin(bars[1]), std::end(bars[1]), iBar1) != std::end(bars[1]);
        if (!found1) continue;
        if( !acceptEvent[Vov][0][iBar0][entry] || !acceptEvent[Vov][1][iBar1][entry]) continue;
        if ( aveEnergy[0][iBar0] > energyMin[0][iBar0][Vov] &&  aveEnergy[1][iBar1] > energyMin[1][iBar1][Vov]){    
	  float eRatio = aveEnergy[0][iBar0]/aveEnergy[1][iBar1] ;
	  float dt = aveTime[0][iBar0]  - aveTime[1][iBar1] ;
	  dt = dt - offset[Vov][iBar0][iBar1];
	  float dtCorr = dt - fitFunCorr2[Vov][iBar0][iBar1]->Eval(eRatio);
	  h_dtAve_energyRatioCorr[Vov][iBar0][iBar1]->Fill(dtCorr);
	}
      }
    }
  }// - end third loop over events
  
  

  // -- single bar
  map<float, map<int, map<int,TF1*> > > fGaus;
  map<float, map<int, map<int,TF1*> > > fGausCorr;

  for (int iMod = 0; iMod < 2; iMod++){
    for (int iBar=0; iBar < 16; iBar++){

      bool found = std::find(std::begin(bars[iMod]), std::end(bars[iMod]), iBar) != std::end(bars[iMod]);                                                                   
      if (!found) continue;  
      
      for (unsigned int ivov = 0; ivov < Vovs.size(); ivov++){
        float Vov = Vovs[ivov];     

	fGaus[Vov][iMod][iBar] = new TF1(Form("fGaus_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov),"gaus",-5000, 5000);
	fGaus[Vov][iMod][iBar]->SetRange(h_dtLR[Vov][iMod][iBar]->GetMean()-2*h_dtLR[Vov][iMod][iBar]->GetRMS(), 
					 h_dtLR[Vov][iMod][iBar]->GetMean()+2*h_dtLR[Vov][iMod][iBar]->GetRMS());
	h_dtLR[Vov][iMod][iBar]->Fit(fGaus[Vov][iMod][iBar],"QSR");
	fGaus[Vov][iMod][iBar]->SetRange(fGaus[Vov][iMod][iBar]->GetParameter(1)-2*fGaus[Vov][iMod][iBar]->GetParameter(2),
					 fGaus[Vov][iMod][iBar]->GetParameter(1)+2*fGaus[Vov][iMod][iBar]->GetParameter(2));
	h_dtLR[Vov][iMod][iBar]->Fit(fGaus[Vov][iMod][iBar],"QSR");
	
	fGausCorr[Vov][iMod][iBar] = new TF1(Form("fGausCorr_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov),"gaus",-5000, 5000);
	fGausCorr[Vov][iMod][iBar]->SetRange(h_dtLR_energyRatioCorr[Vov][iMod][iBar]->GetMean()-2*h_dtLR_energyRatioCorr[Vov][iMod][iBar]->GetRMS(), 
					     h_dtLR_energyRatioCorr[Vov][iMod][iBar]->GetMean()+2*h_dtLR_energyRatioCorr[Vov][iMod][iBar]->GetRMS());
	h_dtLR_energyRatioCorr[Vov][iMod][iBar]->Fit(fGausCorr[Vov][iMod][iBar],"QSR");
	fGausCorr[Vov][iMod][iBar]->SetRange(fGausCorr[Vov][iMod][iBar]->GetParameter(1)-2*fGausCorr[Vov][iMod][iBar]->GetParameter(2),
					     fGausCorr[Vov][iMod][iBar]->GetParameter(1)+2*fGausCorr[Vov][iMod][iBar]->GetParameter(2));
	h_dtLR_energyRatioCorr[Vov][iMod][iBar]->Fit(fGausCorr[Vov][iMod][iBar],"QSR");
	
	cout << "Vov = "<<Vov<<"   M0D " << iMod <<  "     bar " << iBar << 
	  "     sigma(L-R)/2 = " << 0.5*(h_dtLR_energyRatioCorr[Vov][iMod][iBar]->GetFunction(Form("fGausCorr_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov))->GetParameter(2))<<endl;
      }
    }
  }

  // -- 2 bars
  map<float, map<int, map<int,TF1*> > > fGaus2;
  map<float, map<int, map<int,TF1*> > > fGausCorr2;

  for (int iBar0 = 0; iBar0 < 16; iBar0++){                                                                                                                                 
    bool found0 = std::find(std::begin(bars[0]), std::end(bars[0]), iBar0) != std::end(bars[0]);                                                                            
    if (!found0) continue;                                                                                                                                                  
    for (int iBar1 = 0; iBar1 < 16; iBar1++){                                                                                                                               
      bool found1 = std::find(std::begin(bars[1]), std::end(bars[1]), iBar1) != std::end(bars[1]);                                                                          
      if (!found1) continue;
    
      for (unsigned int ivov = 0; ivov < Vovs.size(); ivov++){
        float Vov = Vovs[ivov];
	
	fGaus2[Vov][iBar0][iBar1] = new TF1(Form("fGaus2_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov),"gaus", -5000, 5000);
	fGaus2[Vov][iBar0][iBar1]->SetRange(h_dtAve[Vov][iBar0][iBar1]->GetMean()-2*h_dtAve[Vov][iBar0][iBar1]->GetRMS(),
					    h_dtAve[Vov][iBar0][iBar1]->GetMean()+2*h_dtAve[Vov][iBar0][iBar1]->GetRMS());
	h_dtAve[Vov][iBar0][iBar1]->Fit(fGaus2[Vov][iBar0][iBar1],"QSR");
	fGaus2[Vov][iBar0][iBar1]->SetRange(fGaus2[Vov][iBar0][iBar1]->GetParameter(1)-2*fGaus2[Vov][iBar0][iBar1]->GetParameter(2), 
					    fGaus2[Vov][iBar0][iBar1]->GetParameter(1)+2*fGaus2[Vov][iBar0][iBar1]->GetParameter(2));
	h_dtAve[Vov][iBar0][iBar1]->Fit(fGaus2[Vov][iBar0][iBar1],"QSR");
	
	fGausCorr2[Vov][iBar0][iBar1] = new TF1(Form("fGausCorr2_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov),"gaus", -5000, 5000);  
	fGausCorr2[Vov][iBar0][iBar1]->SetRange(h_dtAve_energyRatioCorr[Vov][iBar0][iBar1]->GetMean()-2*h_dtAve_energyRatioCorr[Vov][iBar0][iBar1]->GetRMS(), 
						h_dtAve_energyRatioCorr[Vov][iBar0][iBar1]->GetMean()+2*h_dtAve_energyRatioCorr[Vov][iBar0][iBar1]->GetRMS());
	h_dtAve_energyRatioCorr[Vov][iBar0][iBar1]->Fit(fGausCorr2[Vov][iBar0][iBar1],"QSR");
	fGausCorr2[Vov][iBar0][iBar1]->SetRange(fGausCorr2[Vov][iBar0][iBar1]->GetParameter(1)-2*fGausCorr2[Vov][iBar0][iBar1]->GetParameter(2), 
						fGausCorr2[Vov][iBar0][iBar1]->GetParameter(1)+2*fGausCorr2[Vov][iBar0][iBar1]->GetParameter(2));
	h_dtAve_energyRatioCorr[Vov][iBar0][iBar1]->Fit(fGausCorr2[Vov][iBar0][iBar1],"QSR");
	cout << "Vov = "<< Vov << "   BAR "<<iBar0 << "   BAR "<<iBar1<<"   sigma(tAve1-tAve2)/sqrt(2) = " << fGausCorr2[Vov][iBar0][iBar1]->GetParameter(2) / sqrt(2) <<endl;
      }
    }
  }

    


  // ======  save histograms in a file
  //string foutName = Form("plots/analysisCoincidenceBars_run%d_Vov%.01f.root", run, Vov);
  //if (!pedestals) foutName = Form("plots/analysisCoincidenceBars_run%d_Vov%.01f_noPedSub.root", run, Vov);
  string foutName = Form("plots/analysisCoincidenceBars_run%d.root",run);
  if (!pedestals) foutName = Form("plots/analysisCoincidenceBars_run%d_noPedSub.root",run);
  TFile *fout = new TFile(foutName.c_str(),"recreate");

  for (int iMod = 0; iMod < 2; iMod++){
    for (int iBar=0; iBar < 16; iBar++){

      bool found = std::find(std::begin(bars[iMod]), std::end(bars[iMod]), iBar) != std::end(bars[iMod]);
      if (!found) continue;    

      for (unsigned int ivov = 0; ivov < Vovs.size(); ivov++){
	float Vov = Vovs[ivov];
	h_energyL[Vov][iMod][iBar]->Write();
	h_energyR[Vov][iMod][iBar]->Write();
	h_energyLR[Vov][iMod][iBar]->Write();
	h_dtLR[Vov][iMod][iBar]->Write();
	h_dtLR_energyRatioCorr[Vov][iMod][iBar]->Write();
	h_energyRatio[Vov][iMod][iBar]->Write();
	p_dtLR_vs_energyRatio[Vov][iMod][iBar]->Write();
      }
    }
  }


  for (int iBar0 = 0; iBar0 < 16; iBar0++){
    bool found0 = std::find(std::begin(bars[0]), std::end(bars[0]), iBar0) != std::end(bars[0]);
    if (!found0) continue;
    for (int iBar1 = 0; iBar1 < 16; iBar1++){
      bool found1 = std::find(std::begin(bars[1]), std::end(bars[1]), iBar1) != std::end(bars[1]);
      if (!found1) continue; 
      for (unsigned int ivov = 0; ivov < Vovs.size(); ivov++){
        float Vov = Vovs[ivov];
	h_dtAve[Vov][iBar0][iBar1]->Write();
	h_dtAve_energyRatioCorr[Vov][iBar0][iBar1]->Write();
	h_energyAveRatio[Vov][iBar0][iBar1]->Write();
	p_dtAve_vs_energyRatio[Vov][iBar0][iBar1]->Write();
	h_tmp[Vov][iBar0][iBar1]->Write();
      }
    }
  }
  
  fout->Close();


  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);


  // ======== PLOT 
  
  // -- deltat(L-R)
  cout<< "Printing dt(L-R) plots"<<endl;

  for (int iMod = 0; iMod < 2; iMod++){
    for (int iBar=0; iBar < 16; iBar++){
      
      bool found = std::find(std::begin(bars[iMod]), std::end(bars[iMod]), iBar) != std::end(bars[iMod]);
      if (!found) continue;

      for (unsigned int ivov = 0; ivov < Vovs.size(); ivov++){
        float Vov = Vovs[ivov];
      
	cout << "Vov = " << Vov <<endl;
	
	TCanvas *c = new TCanvas("c","c", 700, 600);
	h_dtLR_energyRatioCorr[Vov][iMod][iBar]->GetXaxis()->SetTitle("t_L - t_R (ps)");
	h_dtLR_energyRatioCorr[Vov][iMod][iBar]->GetXaxis()->SetRangeUser(-1000, 1000);
	h_dtLR_energyRatioCorr[Vov][iMod][iBar]->GetYaxis()->SetRangeUser( 0, h_dtLR_energyRatioCorr[Vov][iMod][iBar]->GetMaximum()*1.3);
	h_dtLR_energyRatioCorr[Vov][iMod][iBar]->GetFunction(Form("fGausCorr_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov))->SetLineColor(4);
	h_dtLR_energyRatioCorr[Vov][iMod][iBar]->SetLineColor(4);
	h_dtLR_energyRatioCorr[Vov][iMod][iBar]->Draw();
	h_dtLR[Vov][iMod][iBar]->SetLineColor(2);
	h_dtLR[Vov][iMod][iBar]->GetFunction(Form("fGaus_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov))->SetLineColor(2);
	h_dtLR[Vov][iMod][iBar]->Draw("same");
	float sigmaCorr = h_dtLR_energyRatioCorr[Vov][iMod][iBar]->GetFunction(Form("fGausCorr_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov))->GetParameter(2);
	float sigma  = h_dtLR[Vov][iMod][iBar]->GetFunction(Form("fGaus_MOD%02d_BAR%02d_Vov%.01f",iMod,iBar,Vov))->GetParameter(2);
	TLatex *lt  = new TLatex(0.15, 0.85, Form("raw: #sigma_{t} = %.01f ps", sigma));
	TLatex *ltc = new TLatex(0.15, 0.80, Form("corr.: #sigma_{t} = %.01f ps", sigmaCorr));
	lt->SetNDC();
	ltc->SetNDC();
	lt->SetTextColor(2);
	ltc->SetTextColor(4);
	lt->Draw();
	ltc->Draw();
	c->Print(Form("%s/c_dtLR_MOD%02d_BAR%02d_Vov%.01f.png",plotDir.c_str(),iMod,iBar,Vov));
	c->Print(Form("%s/c_dtLR_MOD%02d_BAR%02d_Vov%.01f.pdf",plotDir.c_str(),iMod,iBar,Vov));
		
	TCanvas *c2 = new TCanvas("c2","c2", 700, 600);
	h2_dtLR_vs_energyRatio[Vov][iMod][iBar]->GetXaxis()->SetTitle("energy ratio");
	h2_dtLR_vs_energyRatio[Vov][iMod][iBar]->GetYaxis()->SetTitle("t_L - t_R (ps)");
	h2_dtLR_vs_energyRatio[Vov][iMod][iBar]->GetXaxis()->SetRangeUser(0.5, 1.5);
	h2_dtLR_vs_energyRatio[Vov][iMod][iBar]->GetYaxis()->SetRangeUser(-1000, 1000);
	h2_dtLR_vs_energyRatio[Vov][iMod][iBar]->Draw("colz");
	p_dtLR_vs_energyRatio[Vov][iMod][iBar]->SetMarkerStyle(20);
	p_dtLR_vs_energyRatio[Vov][iMod][iBar]->SetMarkerSize(0.6);
	p_dtLR_vs_energyRatio[Vov][iMod][iBar]->Draw("same");
	c2->Print(Form("%s/c_dtLR_vs_energyRatio_MOD%02d_BAR%02d_Vov%.01f.png",plotDir.c_str(),iMod,iBar,Vov));
	c2->Print(Form("%s/c_dtLR_vs_energyRatio_MOD%02d_BAR%02d_Vov%.01f.pdf",plotDir.c_str(),iMod,iBar,Vov));

	cout << "Deleting canvas..." <<endl;

	delete c;
	delete c2;
      }
    }
  }

  cout << "Printing coincidence plots ... " << endl;

  for (int iBar0 = 0; iBar0 < 16; iBar0++){
    bool found0 = std::find(std::begin(bars[0]), std::end(bars[0]), iBar0) != std::end(bars[0]);
    if (!found0) continue;
    for (int iBar1 = 0; iBar1 < 16; iBar1++){
      bool found1 = std::find(std::begin(bars[1]), std::end(bars[1]), iBar1) != std::end(bars[1]);
      if (!found1) continue;

      for (unsigned int ivov = 0; ivov < Vovs.size(); ivov++){
        float Vov = Vovs[ivov];

	// energy
	TCanvas *c_energy = new TCanvas("c_energy","c_energy", 700, 600);
	c_energy->SetLogy();
	h_energyLR[Vov][0][iBar0]->SetLineColor(2);
	h_energyLR[Vov][1][iBar1]->SetLineColor(4);
	h_energyLR[Vov][0][iBar0]-> GetXaxis()->SetTitle("energy (ADC)");
	h_energyLR[Vov][0][iBar0]->Draw("");
	h_energyLR[Vov][1][iBar1]->Draw("same");
	TLine *line0 = new TLine(energyMin[Vov][0][iBar0], 0, energyMin[Vov][0][iBar0], h_energyLR[Vov][0][iBar0]->GetMaximum());
	line0->SetLineColor(2);
	line0->SetLineStyle(2);
	TLine *line1 = new TLine(energyMin[Vov][1][iBar1], 0, energyMin[Vov][1][iBar1], h_energyLR[Vov][1][iBar1]->GetMaximum());
	line1->SetLineColor(4);
	line1->SetLineStyle(2);
	line0->Draw("same");
	line1->Draw("same");
	c_energy->Print(Form("%s/c_energy_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f.png",plotDir.c_str(),iBar0,iBar1,Vov));
	c_energy->Print(Form("%s/c_energy_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f.pdf",plotDir.c_str(),iBar1,iBar1,Vov));
	
	// delta tAverage 
	TCanvas *c_dtAverage = new TCanvas("c_dtAverage","c_dtAverage", 700, 600);
	h_dtAve_energyRatioCorr[Vov][iBar0][iBar1]->GetXaxis()->SetTitle("t_{Ave1} - t_{Ave2} (ps)");
	h_dtAve_energyRatioCorr[Vov][iBar0][iBar1]->GetXaxis()->SetRangeUser(-1000, 1000);
	h_dtAve_energyRatioCorr[Vov][iBar0][iBar1]->GetYaxis()->SetRangeUser(0, h_dtAve_energyRatioCorr[Vov][iBar0][iBar1]->GetMaximum()*1.3);
	h_dtAve_energyRatioCorr[Vov][iBar0][iBar1]->SetLineColor(4);
	h_dtAve_energyRatioCorr[Vov][iBar0][iBar1]->GetFunction(Form("fGausCorr2_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov))->SetLineColor(4);
	h_dtAve_energyRatioCorr[Vov][iBar0][iBar1]->Draw();
	h_dtAve[Vov][iBar0][iBar1]->SetLineColor(2);
	h_dtAve[Vov][iBar0][iBar1]->GetFunction(Form("fGaus2_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov))->SetLineColor(2);
	h_dtAve[Vov][iBar0][iBar1]->Draw("same");
	float sigma2 = h_dtAve[Vov][iBar0][iBar1]->GetFunction(Form("fGaus2_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov))->GetParameter(2);
	float sigmac2 = h_dtAve_energyRatioCorr[Vov][iBar0][iBar1]->GetFunction(Form("fGausCorr2_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f",iBar0,iBar1,Vov))->GetParameter(2);
	TLatex *lt2 = new TLatex(0.15, 0.85, Form("raw: #sigma_{t} = %.01f ps", sigma2));
	TLatex *ltc2 = new TLatex(0.15, 0.80, Form("corr.: #sigma_{t} = %.01f ps", sigmac2));
	lt2->SetNDC();
	ltc2->SetNDC();
	lt2->SetTextColor(2);
	ltc2->SetTextColor(4);
	lt2->Draw();
	ltc2->Draw();
	c_dtAverage->Print(Form("%s/c_dtAverage_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f.png",plotDir.c_str(),iBar0,iBar1,Vov));
	c_dtAverage->Print(Form("%s/c_dtAverage_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f.pdf",plotDir.c_str(),iBar0,iBar1,Vov));
	
	
	TCanvas *c3 = new TCanvas("c3","c3", 700, 600);
	h2_dtAve_vs_energyRatio[Vov][iBar0][iBar1]->GetXaxis()->SetTitle("energy ratio");
	h2_dtAve_vs_energyRatio[Vov][iBar0][iBar1]->GetYaxis()->SetTitle("t_{Ave1} - t_{Ave2} (ps)");
	h2_dtAve_vs_energyRatio[Vov][iBar0][iBar1]->GetYaxis()->SetRangeUser(-1000, 1000);
	h2_dtAve_vs_energyRatio[Vov][iBar0][iBar1]->Draw("colz");
	p_dtAve_vs_energyRatio[Vov][iBar0][iBar1]->SetMarkerStyle(20);
	p_dtAve_vs_energyRatio[Vov][iBar0][iBar1]->SetMarkerSize(0.7);
	p_dtAve_vs_energyRatio[Vov][iBar0][iBar1]->Draw("same");
	c3->Print(Form("%s/c_dtAve_vs_energyRatio_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f.png",plotDir.c_str(),iBar0,iBar1,Vov));
	c3->Print(Form("%s/c_dtAve_vs_energyRatio_MOD00_BAR%02d_MOD01_BAR%02d_Vov%.01f.pdf",plotDir.c_str(),iBar0,iBar1,Vov));
	
	delete c_energy;
	delete c_dtAverage;
	delete c3;
      }
    }
  }
 
}
