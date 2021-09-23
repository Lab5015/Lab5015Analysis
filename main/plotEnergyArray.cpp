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
#include "TPaveStats.h"
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
#include "TLine.h"

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

  gStyle->SetPalette(kRainBow);


  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);

  int run  = opts.GetOpt<int>("Inputs.run");

  std::vector<float> Vovs = opts.GetOpt<std::vector<float> >("Plots.Vov");
  std::vector<float> energyMins = opts.GetOpt<std::vector<float> >("Plots.energyMin");
  std::vector<float> energyMaxs = opts.GetOpt<std::vector<float> >("Plots.energyMax");
  map<float, float> enMins;
  map<float, float> enMaxs;
  for (unsigned int ivov = 0; ivov < Vovs.size(); ivov++){
    enMins[Vovs[ivov]] = energyMins[ivov];
    enMaxs[Vovs[ivov]] = energyMaxs[ivov];
  }
  

  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  system(Form("mkdir -p %s",plotDir.c_str()));


  // -- choose vth1
  //int mystep2 = 211102;
  int myvth1  = 20;

  // -- max energySum
  float maxEnergySum = 800;
  int maxNbars = 3;

  TChain* tree = new TChain("data","data");
  //tree->Add(Form("/data/TOFHIR2/MTDTB_CERN_Jul21/reco/%d/*_ped_e.root", run));
  tree->Add(Form("/data/tofhir2/h8/reco/%d/*_ped_e.root", run));

  //--- define branches
  float step1, step2;
  int channelIdx[128];
  std::vector<float> *tot = 0;
  std::vector<float> *energy = 0;
  
  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("step1",  1); tree -> SetBranchAddress("step1",  &step1);
  tree -> SetBranchStatus("step2",  1); tree -> SetBranchAddress("step2",  &step2);
  tree -> SetBranchStatus("channelIdx",  1); tree -> SetBranchAddress("channelIdx",  channelIdx);
  tree -> SetBranchStatus("tot",    1); tree -> SetBranchAddress("tot",       &tot);
  tree -> SetBranchStatus("energy", 1); tree -> SetBranchAddress("energy", &energy);

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
      std::cout<<iMod << "  bar = "<< iBar<< "   chL = " << chL[iMod][iBar] << "   chR = "<<chR[iMod][iBar]<<std::endl;
    }
  }

  // -- book histograms 

  map<float, map<int, map<int,TH1F*> > > h_energyLR;
  map<float, map<int, map<int,TH1F*> > > h_energyLR_selEnergySum;
  map<float, map<int, map<int,TH1F*> > > h_energyLR_selEnergyFraction;
  map<float, map<int, map<int,TH1F*> > > h_energyLR_selEnergySumNbars;
  map<float, map<int, map<int,TH1F*> > > h_energyLR_selNbars;
  map<float, map<int, TH1F*> >  h_energySum;
  map<float, map<int, TH1F*> >  h_nBars;
  map<float, map<int, map<int,TH2F*> > > h_energyLR_vs_energySum;
  map<float, map<int, map<int,TH2F*> > > h_energyFraction_vs_energySum;
  map<float, map<int, map<int,TH2F*> > > h_energyFraction_vs_energyLR;
  map<float, map<int, map<int,TH2F*> > > h_energyLR_vs_nBars;
  map<float, map<int, TH2F*> >  h_energySum_vs_nBars;
  map<float, map<int, map<int,TH1F*> > > h_energyMap_Shower;
  map<float, map<int, map<int,TH1F*> > > h_energyMap_MIP;
    
  for (unsigned int ivov = 0; ivov < Vovs.size(); ivov++){
    float Vov = Vovs[ivov];
    for (int iMod = 0; iMod < 2; iMod++){
      for (int iBar = 0; iBar < 16; iBar++){
	h_energyLR[Vov][iMod][iBar] = new TH1F(Form("h_energyLR_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), Form("h_energyLR_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), 500, 0, 1000);	
	h_energyLR_selEnergySum[Vov][iMod][iBar] = new TH1F(Form("h_energyLR_selEnergySum_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), Form("h_energyLR_selEnergySum_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), 500, 0, 1000);	
	h_energyLR_selEnergySumNbars[Vov][iMod][iBar] = new TH1F(Form("h_energyLR_selEnergySumNbars_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), Form("h_energyLR_selEnergySumNbars_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), 500, 0, 1000);	
	h_energyLR_selNbars[Vov][iMod][iBar] = new TH1F(Form("h_energyLR_selNbars_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), Form("h_energyLR_selNbars_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), 500, 0, 1000);	
	h_energyLR_selEnergyFraction[Vov][iMod][iBar] = new TH1F(Form("h_energyLR_selEnergyFraction_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), Form("h_energyLR_selEnergyFraction_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), 500, 0, 1000);	
	h_energyLR_vs_energySum[Vov][iMod][iBar] = new TH2F(Form("h_energyLR_vs_energySum_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), Form("h_energyLR_vs_energySum_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), 500, 0, 10000, 500, 0, 1000);	
	h_energyFraction_vs_energySum[Vov][iMod][iBar] = new TH2F(Form("h_energyFraction_vs_energySum_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), Form("h_energyFraction_vs_energySum_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), 500, 0, 10000, 110, 0, 1.1);	
	h_energyFraction_vs_energyLR[Vov][iMod][iBar] = new TH2F(Form("h_energyFraction_vs_energyLR_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), Form("h_energyFraction_vs_energyLR_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), 500, 0, 1000, 110, 0, 1.1);	
	h_energyLR_vs_nBars[Vov][iMod][iBar] = new TH2F(Form("h_energyLR_vs_nBars_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), Form("h_energyLR_vs_nBars_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), 17, -0.5, 16.5, 500, 0, 1000);	
	h_energyMap_MIP[Vov][iMod][iBar] = new TH1F(Form("h_energyMap_MIP_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), Form("h_energyMap_MIP_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov),17,-0.5,16.5);
	h_energyMap_Shower[Vov][iMod][iBar] = new TH1F(Form("h_energyMap_Shower_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov), Form("h_energyMap_Shower_array%02d_bar%02d_Vov%.01f",iMod,iBar,Vov),17,-0.5,16.5);
      }
      h_energySum[Vov][iMod] = new TH1F(Form("h_energySum_array%02d_Vov%.01f",iMod,Vov), Form("h_energySum_array%02d_Vov%.01f",iMod,Vov),5000,0,10000);
      h_nBars[Vov][iMod] = new TH1F(Form("h_nBars_array%02d_Vov%.01f",iMod,Vov), Form("h_nBars_array%02d_Vov%.01f",iMod,Vov),17,-0.5,16.5);
      h_energySum_vs_nBars[Vov][iMod] = new TH2F(Form("h_energySum_vs_nBars_array%02d_Vov%.01f",iMod,Vov), Form("h_energySum_vs_nBars_array%02d_Vov%.01f",iMod,Vov), 17, -0.5, 16.5, 500, 0, 10000);	
    }
  }


  // -- first loop over events
  cout << "First loop over events to find the mip peak" <<endl;
  for (int entry = 0; entry < maxEntries; entry++){
    
    tree->GetEntry(entry);

    if( entry%1000 == 0 ) std::cout << ">>> Reading entry " << entry << " / " << nEntries << "\r" << std::flush;

    //if (step2 != mystep2 ) continue;
    int vth1 = step2/10000-1;
    if (vth1 != myvth1 ) continue;
    
    float Vov = step1;
    bool goodVov = std::find(std::begin(Vovs), std::end(Vovs), Vov) != std::end(Vovs);   
    if (!goodVov) continue;
   
    float aveEnergy[2][16];
    float energySum[2];
    int n[2];
    
    for (int iMod = 0; iMod < 2; iMod++){
      
      energySum[iMod] =  0;
      n[iMod] = 0;

      for (int iBar = 0; iBar < 16; iBar++){
	
	aveEnergy[iMod][iBar] = -999.;

	if (channelIdx[chL[iMod][iBar]] < 0 || channelIdx[chR[iMod][iBar]] < 0);
	
	float totL =(*tot)[channelIdx[chL[iMod][iBar]]] / 1000. ;
	float totR =(*tot)[channelIdx[chR[iMod][iBar]]] / 1000. ;
	
	if (totL < 0 || totL > 20.) continue;
	if (totR < 0 || totR > 20.) continue;
	
	float enL =(*energy)[channelIdx[chL[iMod][iBar]]];
	float enR =(*energy)[channelIdx[chR[iMod][iBar]]];
	aveEnergy[iMod][iBar] = 0.5 * (enL+enR);
	
	if (aveEnergy[iMod][iBar]>0) {
	  energySum[iMod]+=aveEnergy[iMod][iBar]; 
	  n[iMod]+=1;
 	  h_energyLR[Vov][iMod][iBar]->Fill(aveEnergy[iMod][iBar]);
	}
      }//end loop over bars in iMod

      h_energySum[Vov][iMod]->Fill(energySum[iMod]);
      h_nBars[Vov][iMod]->Fill(n[iMod]);
      h_energySum_vs_nBars[Vov][iMod] ->Fill(n[iMod], energySum[iMod]);

      for (int iBar = 0; iBar < 16; iBar++){
	if ( aveEnergy[iMod][iBar] > 0) {
	  
	  h_energyLR_vs_energySum[Vov][iMod][iBar]->Fill(energySum[iMod], aveEnergy[iMod][iBar]);
	  //if (n[iMod] <= maxNbars) h_energyLR_vs_energySum[Vov][iMod][iBar]->Fill(energySum[iMod], aveEnergy[iMod][iBar]);
	  h_energyFraction_vs_energySum[Vov][iMod][iBar]->Fill(energySum[iMod], aveEnergy[iMod][iBar]/energySum[iMod]);
	  
	  h_energyFraction_vs_energyLR[Vov][iMod][iBar]->Fill(aveEnergy[iMod][iBar], aveEnergy[iMod][iBar]/energySum[iMod]);
	  h_energyLR_vs_nBars[Vov][iMod][iBar]->Fill(n[iMod], aveEnergy[iMod][iBar]);
	  //if (energySum[iMod] < maxEnergySum ) h_energyLR_vs_nBars[Vov][iMod][iBar]->Fill(n[iMod], aveEnergy[iMod][iBar]);
 
	  if (energySum[iMod] < maxEnergySum ) h_energyLR_selEnergySum[Vov][iMod][iBar]->Fill(aveEnergy[iMod][iBar]);
	  if (energySum[iMod] < maxEnergySum && n[iMod] <= maxNbars) h_energyLR_selEnergySumNbars[Vov][iMod][iBar]->Fill(aveEnergy[iMod][iBar]);
	  if (n[iMod] <= maxNbars) h_energyLR_selNbars[Vov][iMod][iBar]->Fill(aveEnergy[iMod][iBar]);
	  if (aveEnergy[iMod][iBar]/energySum[iMod] > 0.9 ) h_energyLR_selEnergyFraction[Vov][iMod][iBar]->Fill(aveEnergy[iMod][iBar]);
	  
	  // -- showering events
	  if (energySum[iMod] > maxEnergySum && n[iMod] > maxNbars) {
	    for (int jBar = 0; jBar < 16; jBar++){
	      h_energyMap_Shower[Vov][iMod][iBar]->Fill(jBar, aveEnergy[iMod][jBar]/energySum[iMod]);
	    }
	  }

	  // -- MIP events
	  //if (energySum[iMod] < maxEnergySum && n[iMod] < 6 && aveEnergy[iMod][iBar] > 400 && aveEnergy[iMod][iBar] < 600) {
	  if (energySum[iMod] < maxEnergySum && n[iMod] <= maxNbars) {
	    for (int jBar = 0; jBar < 16; jBar++){
	      if ( Vov == 5 && aveEnergy[iMod][iBar]> 400 && aveEnergy[iMod][iBar] < maxEnergySum) 
		h_energyMap_MIP[Vov][iMod][iBar]->Fill(jBar, aveEnergy[iMod][jBar]/energySum[iMod]);
	      if ( Vov == 1.5 && aveEnergy[iMod][iBar]> 50 && aveEnergy[iMod][iBar] < maxEnergySum) 
		h_energyMap_MIP[Vov][iMod][iBar]->Fill(jBar, aveEnergy[iMod][jBar]/energySum[iMod]);
	    }
	  }
	}
      }


    }// end loop over module
  }// end loop over events


  // ======  save histograms in a file
  string foutName = Form("plots/plotEnergyArray_run%d.root",run);

  TFile *fout = new TFile(foutName.c_str(),"recreate");

  for (unsigned int ivov = 0; ivov < Vovs.size(); ivov++){
    float Vov = Vovs[ivov];
    for (int iMod = 0; iMod < 2; iMod++){
      for (int iBar=0; iBar < 16; iBar++){
	if (h_energyLR[Vov][iMod][iBar]->GetEntries()==0) continue;
	h_energyLR[Vov][iMod][iBar]->Write();
	h_energyLR_selEnergySum[Vov][iMod][iBar]->Write();
	h_energyLR_selEnergySumNbars[Vov][iMod][iBar]->Write();
	h_energyLR_selNbars[Vov][iMod][iBar]->Write();
	h_energyLR_selEnergyFraction[Vov][iMod][iBar]->Write();
	h_energyFraction_vs_energySum[Vov][iMod][iBar]->Write();
	h_energyFraction_vs_energyLR[Vov][iMod][iBar]->Write();
	h_energyLR_vs_nBars[Vov][iMod][iBar]->Write();
	h_energyMap_MIP[Vov][iMod][iBar]->Write();
	h_energyMap_Shower[Vov][iMod][iBar]->Write();
      }
      h_energySum[Vov][iMod]->Write();
      h_nBars[Vov][iMod]->Write();
    }
  }

  fout->Close();


  //gStyle->SetOptFit(1111);
  //gStyle->SetOptTitle(0);


  // ======== PLOT 

  for (unsigned int ivov = 0; ivov < Vovs.size(); ivov++){
    float Vov = Vovs[ivov];
    for (int iMod = 0; iMod < 2; iMod++){
      for (int iBar=0; iBar < 16; iBar++){

	if (h_energyLR[Vov][iMod][iBar]->GetEntries()==0) continue;

	TCanvas *c = new TCanvas("c","c", 700, 600);

	c->cd();
	c->SetLogy();
	h_energyLR[Vov][iMod][iBar]->GetXaxis()->SetTitle("energy [ADC]");
	h_energyLR_selEnergySum[Vov][iMod][iBar]->SetLineColor(kRed);
	h_energyLR_selEnergySumNbars[Vov][iMod][iBar]->SetLineColor(kMagenta);
	h_energyLR_selEnergySumNbars[Vov][iMod][iBar]->SetLineStyle(3);
	h_energyLR_selNbars[Vov][iMod][iBar]->SetLineColor(kGreen+1);
	h_energyLR_selNbars[Vov][iMod][iBar]->SetLineStyle(1);
	h_energyLR_selEnergyFraction[Vov][iMod][iBar]->SetLineColor(kOrange);
	h_energyLR_selEnergyFraction[Vov][iMod][iBar]->SetLineStyle(2);
	h_energyLR[Vov][iMod][iBar]->Draw();
	h_energyLR_selEnergySum[Vov][iMod][iBar]->Draw("same");
	h_energyLR_selEnergySumNbars[Vov][iMod][iBar]->Draw("same");
	h_energyLR_selNbars[Vov][iMod][iBar]->Draw("same");
	h_energyLR_selEnergyFraction[Vov][iMod][iBar]->Draw("same");
	TLegend *leg = new TLegend();
	leg->AddEntry(h_energyLR[Vov][iMod][iBar],"all","L");
	leg->AddEntry(h_energyLR_selEnergySum[Vov][iMod][iBar],"E_{array} < 800","L");
	leg->AddEntry(h_energyLR_selEnergySumNbars[Vov][iMod][iBar],"E_{array} < 800 && Nbars < 5","L");
	leg->AddEntry(h_energyLR_selNbars[Vov][iMod][iBar],"Nbars < 5","L");
	leg->AddEntry(h_energyLR_selEnergyFraction[Vov][iMod][iBar],"E_{bar}/E_{array} > 0.9","L");
	leg->Draw("same");
	c->Print(Form("%s/c_energy_array%02d_bar%02d_Vov%.01f.png",plotDir.c_str(),iMod,iBar,Vov));
	c->Print(Form("%s/c_energy_array%02d_bar%02d_Vov%.01f.pdf",plotDir.c_str(),iMod,iBar,Vov));


	c->cd();
	c->SetLogy(0);
	c->SetLogz();
	h_energyLR_vs_energySum[Vov][iMod][iBar]->GetXaxis()->SetRangeUser(0, 5000);
	h_energyLR_vs_energySum[Vov][iMod][iBar]->GetXaxis()->SetTitle("energy sum [ADC]");
	h_energyLR_vs_energySum[Vov][iMod][iBar]->GetYaxis()->SetTitle("energy [ADC]");
	h_energyLR_vs_energySum[Vov][iMod][iBar]->Draw("colz");
	c->Print(Form("%s/c_energy_energySum_array%02d_bar%02d_Vov%.01f.png",plotDir.c_str(),iMod,iBar,Vov));
	c->Print(Form("%s/c_energy_energySum_array%02d_bar%02d_Vov%.01f.pdf",plotDir.c_str(),iMod,iBar,Vov));

	c->cd();
	c->SetLogz();
	h_energyFraction_vs_energySum[Vov][iMod][iBar]->GetXaxis()->SetRangeUser(0, 5000);
	h_energyFraction_vs_energySum[Vov][iMod][iBar]->GetXaxis()->SetTitle("energy sum [ADC]");
	h_energyFraction_vs_energySum[Vov][iMod][iBar]->GetYaxis()->SetTitle("energy fraction");
	h_energyFraction_vs_energySum[Vov][iMod][iBar]->Draw("colz");
	c->Print(Form("%s/c_energyFraction_energySum_array%02d_bar%02d_Vov%.01f.png",plotDir.c_str(),iMod,iBar,Vov));
	c->Print(Form("%s/c_energyFraction_energySum_array%02d_bar%02d_Vov%.01f.pdf",plotDir.c_str(),iMod,iBar,Vov));


	c->cd();
	c->SetLogz();
	h_energyFraction_vs_energyLR[Vov][iMod][iBar]->GetXaxis()->SetRangeUser(0, 1000);
	h_energyFraction_vs_energyLR[Vov][iMod][iBar]->GetXaxis()->SetTitle("energy [ADC]");
	h_energyFraction_vs_energyLR[Vov][iMod][iBar]->GetYaxis()->SetTitle("energy fraction");
	h_energyFraction_vs_energyLR[Vov][iMod][iBar]->Draw("colz");
	c->Print(Form("%s/c_energyFraction_energyLR_array%02d_bar%02d_Vov%.01f.png",plotDir.c_str(),iMod,iBar,Vov));
	c->Print(Form("%s/c_energyFraction_energyLR_array%02d_bar%02d_Vov%.01f.pdf",plotDir.c_str(),iMod,iBar,Vov));

	c->cd();
	c->SetLogz();
	h_energyLR_vs_nBars[Vov][iMod][iBar]->GetXaxis()->SetTitle("n bars");
	h_energyLR_vs_nBars[Vov][iMod][iBar]->GetYaxis()->SetTitle("energy [ADC]");
	h_energyLR_vs_nBars[Vov][iMod][iBar]->Draw("colz");
	c->Print(Form("%s/c_energy_nBars_array%02d_bar%02d_Vov%.01f.png",plotDir.c_str(),iMod,iBar,Vov));
	c->Print(Form("%s/c_energy_nBars_array%02d_bar%02d_Vov%.01f.pdf",plotDir.c_str(),iMod,iBar,Vov));


	c->cd();
	h_energyMap_MIP[Vov][iMod][iBar]->GetXaxis()->SetTitle("bar");
	h_energyMap_MIP[Vov][iMod][iBar]->GetYaxis()->SetTitle("energy fraction");
	h_energyMap_MIP[Vov][iMod][iBar]->DrawNormalized("histo");
	c->Print(Form("%s/c_energyMap_MIP_array%02d_bar%02d_Vov%.01f.png",plotDir.c_str(),iMod,iBar,Vov));
	c->Print(Form("%s/c_energyMap_MIP_array%02d_bar%02d_Vov%.01f.pdf",plotDir.c_str(),iMod,iBar,Vov));

	c->cd();
	h_energyMap_Shower[Vov][iMod][iBar]->GetXaxis()->SetTitle("bar");
	h_energyMap_Shower[Vov][iMod][iBar]->GetYaxis()->SetTitle("energy fraction");
	h_energyMap_Shower[Vov][iMod][iBar]->DrawNormalized("histo");
	c->Print(Form("%s/c_energyMap_Shower_array%02d_bar%02d_Vov%.01f.png",plotDir.c_str(),iMod,iBar,Vov));
	c->Print(Form("%s/c_energyMap_Shower_array%02d_bar%02d_Vov%.01f.pdf",plotDir.c_str(),iMod,iBar,Vov));

      }
      
      TCanvas *cc = new TCanvas("cc","cc", 700, 600);
      cc->cd();
      cc->SetLogy();
      h_energySum[Vov][iMod]->GetXaxis()->SetTitle("energy sum [ADC]");
      h_energySum[Vov][iMod]->Draw();
      cc->Print(Form("%s/c_energySum_array%02d_Vov%.01f.png",plotDir.c_str(),iMod,Vov));
      cc->Print(Form("%s/c_energySum_array%02d_Vov%.01f.pdf",plotDir.c_str(),iMod,Vov));

      cc->cd();
      h_nBars[Vov][iMod]->GetXaxis()->SetTitle("N bars");
      h_nBars[Vov][iMod]->Draw();
      cc->Print(Form("%s/c_nBars_array%02d_Vov%.01f.png",plotDir.c_str(),iMod,Vov));
      cc->Print(Form("%s/c_nBars_array%02d_Vov%.01f.pdf",plotDir.c_str(),iMod,Vov));


      cc->cd();
      cc->SetLogy(0);
      cc->SetLogz(1);
      h_energySum_vs_nBars[Vov][iMod]->GetXaxis()->SetTitle("N bars");
      h_energySum_vs_nBars[Vov][iMod]->GetYaxis()->SetTitle("energy sum [ADC]");
      h_energySum_vs_nBars[Vov][iMod]->Draw("colz");
      cc->Print(Form("%s/c_energySum_vs_nBars_array%02d_Vov%.01f.png",plotDir.c_str(),iMod,Vov));
      cc->Print(Form("%s/c_energySum_vs_nBars_array%02d_Vov%.01f.pdf",plotDir.c_str(),iMod,Vov));
      
    }
  }

 
}
