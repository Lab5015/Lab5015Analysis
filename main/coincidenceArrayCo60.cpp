#include "interface/AnalysisUtils.h"
#include "interface/Na22SpectrumAnalyzer.h"
#include "interface/FitUtils.h"
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
#include "TLegend.h"
#include "TSpectrum.h"
#include "TMultiGraph.h"



int main(int argc, char** argv)
{
  setTDRStyle();


  if( argc < 2 )
  {
    std::cout << ">>> coincideceArrayCo60::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }


  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);

  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);


  typedef std::numeric_limits<double> dbl;
  std::cout.precision(dbl::max_digits10);
  
  //--- get parameters
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  //system(Form("rm -r %s/coincidenceResults", plotDir.c_str()));
  system(Form("mkdir -p %s",plotDir.c_str()));
  system(Form("mkdir -p %s/coincidenceResults/",plotDir.c_str()));
  //system(Form("mkdir -p %s/coincidenceResults/tot/",plotDir.c_str()));
  system(Form("mkdir -p %s/coincidenceResults/energy/",plotDir.c_str()));
	system(Form("mkdir -p %s/coincidenceResults/energy/runs",plotDir.c_str()));
	system(Form("mkdir -p %s/coincidenceResults/energy/total",plotDir.c_str()));
  system(Form("mkdir -p %s/coincidenceResults/timeResolution/",plotDir.c_str()));
	system(Form("mkdir -p %s/coincidenceResults/timeResolution/runs",plotDir.c_str()));
	system(Form("mkdir -p %s/coincidenceResults/timeResolution/total",plotDir.c_str()));
	system(Form("mkdir -p %s/coincidenceResults/timeResolution/ratioEnergyCorrection",plotDir.c_str()));
	system(Form("mkdir -p %s/coincidenceResults/timeResolution/timeResolutionComparison",plotDir.c_str()));
	

  //copia il file .php che si trova in una cartella fuori plotDir in plotDir definita nel config file
  system(Form("cp %s/../index.php %s",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/coincidenceResults/",plotDir.c_str(),plotDir.c_str()));
  //system(Form("cp %s/../index.php %s/coincidenceResults/tot/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/coincidenceResults/energy/",plotDir.c_str(),plotDir.c_str()));
	system(Form("cp %s/../index.php %s/coincidenceResults/energy/runs",plotDir.c_str(),plotDir.c_str()));
	system(Form("cp %s/../index.php %s/coincidenceResults/energy/total",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/coincidenceResults/timeResolution/",plotDir.c_str(),plotDir.c_str()));
	system(Form("cp %s/../index.php %s/coincidenceResults/timeResolution/runs",plotDir.c_str(),plotDir.c_str()));
	system(Form("cp %s/../index.php %s/coincidenceResults/timeResolution/total",plotDir.c_str(),plotDir.c_str()));
	system(Form("cp %s/../index.php %s/coincidenceResults/timeResolution/ratioEnergyCorrection",plotDir.c_str(),plotDir.c_str()));
	system(Form("cp %s/../index.php %s/coincidenceResults/timeResolution/timeResolutionComparison",plotDir.c_str(),plotDir.c_str()));

  //--- open files and make the tree chain
  //--- Prima di eseguire cambia il numero del run nel config file
  std::string inputDirStep1 = opts.GetOpt<std::string>("Input.inputDirStep1");
  std::string inputDirStep2 = opts.GetOpt<std::string>("Input.inputDirStep2");
  std::string fileBaseName1 = opts.GetOpt<std::string>("Input.fileBaseName1");
  std::string fileBaseName2 = opts.GetOpt<std::string>("Input.fileBaseName2");
	std::string inputDirData = opts.GetOpt<std::string>("Input.inputDirData");
  std::string fileBaseNameData = opts.GetOpt<std::string>("Input.fileBaseNameData");
  std::string runs = opts.GetOpt<std::string>("Input.runs");


  //Define variables
  std::vector<std::string> listStep1;
  std::vector<std::string> listStep2;
	std::vector<std::string> listData;
	std::vector<int> listRun;

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
    
    for(int run = runMin; run <= runMax; ++run)
    {
      std::string fileNameStep1;
      std::string fileNameStep2;
			std::string fileNameData;
      fileNameStep1 = Form("%s%s1_%s%04d.root",inputDirStep1.c_str(),fileBaseName1.c_str(),fileBaseName2.c_str(),run);
      fileNameStep2 = Form("%s%s2_%s%04d.root",inputDirStep2.c_str(),fileBaseName1.c_str(),fileBaseName2.c_str(),run);
			fileNameData =  Form("%s/%s%04d_t.root",inputDirData.c_str(),fileBaseNameData.c_str(),run);
      std::cout << ">>> Adding file " << fileNameStep1 << std::endl;
      std::cout << ">>> Adding file " << fileNameStep2 << std::endl;
			std::cout << ">>> Adding file " << fileNameData << std::endl;

      listStep1.push_back(fileNameStep1);
      listStep2.push_back(fileNameStep2); 
			listData.push_back(fileNameData);
			listRun.push_back(run);      
    }
  }
  
  //--- define output root file
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileName");
  TFile* outFile = TFile::Open(outFileName.c_str(),"RECREATE");

  
  std::map<std::string,TTree*> trees;
	std::map<std::string,TChain*> dataTrees;
  std::map<std::string,TTree*> trees2;
  std::vector<std::string> stepLabels;
	std::map< std::string, std::vector <float>> energyRanges;
	std::map< std::string, std::map<int,TF1*>> tProfileFitFunc;
	std::map< std::string, std::map<int,TF1*>> energyRatioFitFunc;
	std::map< std::string, std::map<int,TCanvas*>> c_timeRes_summary;
	std::map< std::string, std::map<int,TLegend*>> l_timeRes_summary;
	std::map< std::string, std::map<int,std::vector <TGraphErrors*>>> gr_timeRes_summary;
	
  

	//open Step2 root file and data root file
	int id = 0;

	for(auto file: listStep2){
	  std::string label;
		std::string labelData;
		std::string label2;
	  TFile* inFile = TFile::Open(file.c_str(),"READ");
	  
	  TList* list = inFile -> GetListOfKeys();
	  TIter next(list);
	  TObject* object = 0;
		while( (object = next()) ){
	    std::string name(object->GetName());
	    std::vector<std::string> tokens = GetTokens(name,'_');
	    std::size_t found; 
			std::size_t found2; 
	    found = name.find("data_");
	    found2 = name.find("dataRes_");
	    if( found!=std::string::npos ){
	      label = Form("%s_%s_%s_%04d", tokens[1].c_str(), tokens[2].c_str(), tokens[3].c_str(), listRun[id]);
				labelData = Form("%s_%s_%04d", tokens[2].c_str(), tokens[3].c_str(), listRun[id]);
	      trees[label] = (TTree*)( inFile->Get(name.c_str()) );
	    }
			 if( found2!=std::string::npos ){
	      label2 = Form("%s_%s_%s_%s_%04d", tokens[1].c_str(), tokens[2].c_str(), tokens[3].c_str(), tokens[4].c_str(), listRun[id]);
	      trees2[label2] = (TTree*)( inFile->Get(name.c_str()) );
	    }
		}				
		dataTrees[labelData] = new TChain("data","data");
		dataTrees[labelData] -> Add(listData[id].c_str());
	  id ++;
	}
	
	
	

	//taking 7-8-9 energyBin from filled Trees
	

	for (auto mapIt : trees){
		
    float  energyBinVal[3];
    float theIndex = -1;
    mapIt.second -> SetBranchStatus("*",0);
    mapIt.second -> SetBranchStatus("energyBinValCo60",  1); mapIt.second -> SetBranchAddress("energyBinValCo60",  energyBinVal);
    mapIt.second -> SetBranchStatus("indexID",  1); mapIt.second -> SetBranchAddress("indexID",   &theIndex);
		int nEntries = mapIt.second->GetEntries();
    for(int entry = 0; entry < nEntries; ++entry){
      mapIt.second -> GetEntry(entry);
			
			std::vector <float> vec;
			if (theIndex > 0){
				for ( int i = 0; i < 3; i++){
					vec.push_back(energyBinVal[i]);
				}
				std::sort(vec.begin(), vec.end());
				energyRanges[mapIt.first] = vec;   
			}
		
		}	
		
	}
	
	for ( auto mapIt : trees2){
		int iBin;
		float tRes;
		int run;
		mapIt.second -> SetBranchStatus("*",0);
		mapIt.second -> SetBranchStatus("timeResolution",  1); mapIt.second -> SetBranchAddress("timeResolution",  &tRes);
		int nEntries = mapIt.second->GetEntries();
		std::vector<std::string> tokens = GetTokens(mapIt.first,'_');
		std::string label_summary = Form( "%s_%s_%s",tokens[0].c_str(), tokens[1].c_str(), tokens[2].c_str());
		std::vector<std::string> tokens2 = GetTokens(tokens[3].c_str(),'n');
		iBin = int( atof(tokens2[2].c_str()));
		run = int( atof(tokens[4].c_str()));
    for(int entry = 0; entry < nEntries; ++entry){
      mapIt.second -> GetEntry(entry);
			
			if ( iBin == 8 || iBin == 9){		
				if (c_timeRes_summary[label_summary][iBin] == NULL){
					c_timeRes_summary[label_summary][iBin] = new TCanvas ( Form("c_time_resolution_comparison%s_energyBin%d", label_summary.c_str(), iBin), Form("c_time_resolution_comparison%s_energyBin%d", label_summary.c_str(), iBin)); 
					TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
					hPad -> Draw();
					gPad -> SetGridy();
					l_timeRes_summary[label_summary][iBin] = new TLegend (0.70,0.40,0.95,0.95);
					l_timeRes_summary[label_summary][iBin] -> SetBorderSize(0);
				}
					
				TGraphErrors* graph = new TGraphErrors();
				graph -> SetPoint (0, gr_timeRes_summary[label_summary][iBin].size()+1, tRes);
				graph -> SetPointError(0, 0, 0);
				graph -> SetMarkerColor(1 + gr_timeRes_summary[label_summary][iBin].size());
				graph -> SetMarkerStyle(25);
				gr_timeRes_summary[label_summary][iBin].push_back((TGraphErrors*) graph);
				l_timeRes_summary[label_summary][iBin] -> AddEntry ( graph, Form("run %4d", run), "p");	
			}
		}	
	}
	

	//bars of interest
	std::map <std::string, std::vector<int>> barMap;
	for( auto index :energyRanges){
			std::vector<std::string> tokens = GetTokens(index.first,'_');
			std::string label = Form("%s_%s_%s", tokens[1].c_str(), tokens[2].c_str(), tokens[3].c_str()); 
			std::vector<std::string> tokens2 = GetTokens(tokens[0],'r');
			std::vector<std::string> tokens3 = GetTokens(tokens2[1],'L');
			int iBar = int(atof(tokens3[0].c_str()));
			barMap[label].push_back(iBar);
	}
	
	//taking TProfiles  and energyRatio for 8-9 energyBin and bars of interest
	id = 0;
	std::map < std::string, TCanvas*> c_step2_TProfile_run;
	for(auto file: listStep2){
	  std::string label;
		std::string labelTProfile;
		std::string labelEnergyRatio;
	  TFile* inFile = TFile::Open(file.c_str(),"READ");
	
		for ( auto mapIt : energyRanges){
			label = mapIt.first;
			std::vector<std::string> tokens = GetTokens(mapIt.first,'_');
			if ( atof(tokens[3].c_str()) != listRun[id]) continue;
			for( int iBin = 8; iBin < 10; iBin++){
				std::string labelEnFile = Form("h1_energyRatio_%s_%s_%s_energyBin%d", tokens[0].c_str(), tokens[1].c_str(), tokens[2].c_str(), iBin);
				labelEnergyRatio = Form("h1_energyRatio_%s_%s_%s_energyBin%d_%04d", tokens[0].c_str(), tokens[1].c_str(), tokens[2].c_str(), iBin, listRun[id]);
				TCanvas * canvas = new TCanvas ( Form("canvas_%s",labelEnergyRatio.c_str() ), Form("canvas_%s",labelEnergyRatio.c_str() ));
				canvas -> cd();
				
				TH1F * histo = (TH1F*) ( inFile -> Get( Form("%s", labelEnFile.c_str())));
				if( !histo) continue;
				energyRatioFitFunc[label][iBin] = new TF1(Form("energyRatioFitFunc%s",Form("bin%d_%s", iBin, label.c_str())),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
    		histo -> Fit(energyRatioFitFunc[label][iBin],"QNRS");
				histo -> Fit(energyRatioFitFunc[label][iBin],"QS+","",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
				outFile -> cd();
				histo -> Write();

				

				labelTProfile = Form("p1_deltaT_vs_energyRatio_%s_%s_%s_energyBin%d", tokens[0].c_str(), tokens[1].c_str(), tokens[2].c_str(), iBin);
				std::string labelCanvas = Form( "%s_%s_%s_energyBin%d_%04d",tokens[0].c_str(), tokens[1].c_str(), tokens[2].c_str(), iBin, listRun[id]);
				
				TProfile* prof = (TProfile*)( inFile->Get(Form("%s",labelTProfile.c_str())));
				if( !prof) continue;
				c_step2_TProfile_run[labelCanvas] = new TCanvas(Form("c_deltaT_vs_energyRatio_%s", labelCanvas.c_str()), Form("c_deltaT_vs_energyRatio_%s", labelCanvas.c_str()));
				
				tProfileFitFunc[label][iBin] = new TF1(Form("tProfileFitFunc%s",Form("bin%d_%s", iBin, label.c_str())),"pol1",-0.88,1.12);
    		prof -> Fit(tProfileFitFunc[label][iBin],"QRS+");
				
				
	      TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
				hPad -> Draw();
	 			gPad -> SetGridy();
				prof -> SetTitle(Form(";energy_{right} / energy_{left};#Deltat [ps]"));
				prof -> Draw("");
				tProfileFitFunc[label][iBin] -> SetLineColor(kRed);
        tProfileFitFunc[label][iBin] -> SetLineWidth(2);
        tProfileFitFunc[label][iBin] -> Draw("same");
				outFile -> cd();
				prof -> Write();
				
				
			}
		}
		id++;
	}
	


					
	
	int chsL[16];
  int chsR[16];
  chsL[0] = 249;
  chsR[0] = 4;
  chsL[1] = 242;
  chsR[1] = 6;
  chsL[2] = 251;
  chsR[2] = 15;
  chsR[3] = 8;
  chsL[4] = 250;
  chsR[4] = 13;
  chsL[5] = 254;
  chsR[5] = 2;
  chsL[6] = 201;
  chsR[6] = 11;
  chsL[7] = 230;
  chsR[7] = 27;
  chsL[8] = 226;
  chsR[8] = 32;
  chsL[9] = 225;
  chsR[9] = 31;
  chsL[10] = 228;
  chsR[10] = 30;
  chsL[11] = 227;
  chsR[11] = 29;
  chsL[12] = 229;
  chsR[12] = 28;
  chsL[13] = 231;
  chsR[13] = 26;
  chsL[14] = 233;
  chsR[14] = 24;
  chsL[15] = 232;
  chsR[15] = 25;

	std::map < std::string, TCanvas*> c_energy_bar_run;
	std::map < std::string, TCanvas*> c_energy_bar_tot;
	std::map < std::string, TCanvas*> c_energy_1173_bar_run;
	std::map < std::string, TCanvas*> c_energy_1173_bar_tot;
	std::map < std::string, TCanvas*> c_energy_1332_bar_run;
	std::map < std::string, TCanvas*> c_energy_1332_bar_tot;
	std::map < std::string, TCanvas*> c_energy_coincidence_run;
	std::map < std::string, TCanvas*> c_energy_coincidence_tot;
	std::map < std::string, TCanvas*> c_timeRes_coincidence_run;
	std::map < std::string, TCanvas*> c_timeRes_coincidence_tot;
	std::map < std::string, std::map <int,TCanvas*>>	c_TProfile_energyRatio_corr;

	std::map < std::string, TH1F*> h_energy_bar_run;
	std::map < std::string, TH1F*> h_energy_bar_tot;
	std::map < std::string, TH1F*> h_energy_1173_bar_run;
	std::map < std::string, TH1F*> h_energy_1173_bar_tot;
	std::map < std::string, TH1F*> h_energy_1332_bar_run;
	std::map < std::string, TH1F*> h_energy_1332_bar_tot;
	std::map < std::string, TH1F*> h_energy_coincidence_run;
	std::map < std::string, TH1F*> h_energy_coincidence_tot;
	std::map < std::string, TH1F*> h_timeRes_coincidence_run;
	std::map < std::string, TH1F*> h_timeRes_coincidence_tot;
	std::map < std::string, TH1F*> h_timeRes_corr_coincidence_run;
	std::map < std::string, TH1F*> h_timeRes_corr_coincidence_tot;
	std::map < std::string, TH1F*> h_timeRes_corr_data_coincidence_run;
	std::map < std::string, TH1F*> h_timeRes_corr_data_coincidence_tot;
	std::map < std::string, std::map <int,TProfile*>>	p1_TProfile_energyRatio_corr;
	std::map < std::string, std::map <int,TH1F*>>	h1_energyRatio;
	
	for (auto mapIt : dataTrees){
		for( int iBar : barMap[mapIt.first]){
				
			std::string label_c = Form( "bar%02dL-R_%s", iBar, mapIt.first.c_str());
			std::vector<std::string> tokens = GetTokens(mapIt.first,'_');
			std::string labelTot_c = Form( "bar%02dL-R_%s_%s", iBar, tokens[0].c_str(), tokens[1].c_str());

			
			
			//all spectrum histo and canvan
			if(c_energy_bar_run[label_c]==NULL){
	      c_energy_bar_run[label_c] = new TCanvas(Form("c_energy_%s",label_c.c_str()), Form("c_energy_%s",label_c.c_str()));
				c_energy_bar_run[label_c] -> cd();
	      TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
	      hPad -> SetTitle(";entries;energy [a.u.]");
				hPad -> Draw();
	 			gPad -> SetGridy();
				h_energy_bar_run[label_c] = new TH1F ( Form("h_energy_%s",label_c.c_str()), Form("h_energy_%s",label_c.c_str()), 1000, 0, 50);
			}
			
			if(c_energy_bar_tot[labelTot_c]==NULL){
	      c_energy_bar_tot[labelTot_c] = new TCanvas(Form("c_energy_%s",labelTot_c.c_str()), Form("c_energy_%s",labelTot_c.c_str()));
	      TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
	      hPad -> SetTitle(";entries;energy [a.u.]");
				hPad -> Draw();
	 			gPad -> SetGridy();
				h_energy_bar_tot[labelTot_c] = new TH1F ( Form("h_energy_%s",labelTot_c.c_str()), Form("h_energy_%s",labelTot_c.c_str()), 1000, 0, 50);
			}
			//1173 peak histo and canvan
			if(c_energy_1173_bar_run[label_c]==NULL){
	      c_energy_1173_bar_run[label_c] = new TCanvas(Form("c_1173_peak_%s",label_c.c_str()), Form("c_1173_peak_%s",label_c.c_str()));
	      TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
	      hPad -> SetTitle(";entries;energy [a.u.]");
				hPad -> Draw();
	 			gPad -> SetGridy();
				h_energy_1173_bar_run[label_c] = new TH1F ( Form("h_1173_peak_%s",label_c.c_str()), Form("h_1173_peak_%s",label_c.c_str()), 1000, 0, 50);
			}
				
			if(c_energy_1173_bar_tot[labelTot_c]==NULL){
	      c_energy_1173_bar_tot[labelTot_c] = new TCanvas(Form("c_1173_peak_%s",labelTot_c.c_str()), Form("c_1173_peak_%s",labelTot_c.c_str()));
	      TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
	      hPad -> SetTitle(";entries;energy [a.u.]");
				hPad -> Draw();
	 			gPad -> SetGridy();
				h_energy_1173_bar_tot[labelTot_c] = new TH1F ( Form("h_1173_peak_%s",labelTot_c.c_str()), Form("h_1173_peak_%s",labelTot_c.c_str()), 1000, 0, 50);
			}
			//1332 peak histo and canvan
			if(c_energy_1332_bar_run[label_c]==NULL){
	      c_energy_1332_bar_run[label_c] = new TCanvas(Form("c_1332_peak_%s",label_c.c_str()), Form("c_1332_peak_%s",label_c.c_str()));
	      TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
	      hPad -> SetTitle(";entries;energy [a.u.]");
				hPad -> Draw();
	 			gPad -> SetGridy();
				h_energy_1332_bar_run[label_c] = new TH1F ( Form("h_1332_peak_%s",label_c.c_str()), Form("h_1332_peak_%s",label_c.c_str()), 1000, 0, 50);
			}
				
			if(c_energy_1332_bar_tot[labelTot_c]==NULL){
	      c_energy_1332_bar_tot[labelTot_c] = new TCanvas(Form("c_1332_peak_%s",labelTot_c.c_str()), Form("c_1332_peak_%s",labelTot_c.c_str()));
	      TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
	      hPad -> SetTitle(";entries;energy [a.u.]");
				hPad -> Draw();
	 			gPad -> SetGridy();
				h_energy_1332_bar_tot[labelTot_c] = new TH1F ( Form("h_1332_peak_%s",labelTot_c.c_str()), Form("h_1332_peak_%s",labelTot_c.c_str()), 1000, 0, 50);
			}
			//peaks -> iBar coincidence with the others (energy and time resolution)
			for (int coinBar : barMap[mapIt.first]){
				if(coinBar != iBar){
					std::string label_coin =  Form( "1173_bar%02dL-R_vs_1332_bar%02dL-R_%s", iBar, coinBar, mapIt.first.c_str());
					std::string labelTot_coin = Form( "1173_bar%02dL-R_vs_1332_bar%02dL-R_%s_%s", iBar, coinBar, tokens[0].c_str(), tokens[1].c_str());
					std::string label_coin2 =  Form( "1332_bar%02dL-R_vs_1173_bar%02dL-R_%s", iBar, coinBar, mapIt.first.c_str());
					std::string labelTot_coin2 = Form( "1332_bar%02dL-R_vs_1173_bar%02dL-R_%s_%s", iBar, coinBar, tokens[0].c_str(), tokens[1].c_str());
				
					if(c_energy_coincidence_run[label_coin]==NULL){
					  c_energy_coincidence_run[label_coin] = new TCanvas(Form("c_energy_coincidence_%s",label_coin.c_str()), Form("c_energy_coincidence_%s",label_coin.c_str()));
					  TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
					  hPad -> SetTitle(";entries;energy [a.u.]");
						hPad -> Draw();
			 			gPad -> SetGridy();
						h_energy_coincidence_run[label_coin] = new TH1F ( Form("h_energy_coincidence_%s",label_coin.c_str()), Form("h_energy_coincidence_%s",label_coin.c_str()), 1000, 0, 50);
					}
						
					if(c_energy_coincidence_tot[labelTot_coin]==NULL){
					  c_energy_coincidence_tot[labelTot_coin] = new TCanvas(Form("c_energy_coincidence_%s",labelTot_coin.c_str()), Form("c_energy_coincidence_%s",labelTot_coin.c_str()));
					  TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
					  hPad -> SetTitle(";entries;energy [a.u.]");
						hPad -> Draw();
			 			gPad -> SetGridy();
						h_energy_coincidence_tot[labelTot_coin] = new TH1F ( Form("h_energy_coincidence_%s",labelTot_coin.c_str()), Form("h_energy_coincidence_%s",labelTot_coin.c_str()), 1000, 0, 50);
					}
				
					if(c_energy_coincidence_run[label_coin2]==NULL){
					  c_energy_coincidence_run[label_coin2] = new TCanvas(Form("c_energy_coincidence_%s",label_coin2.c_str()), Form("c_energy_coincidence_%s",label_coin2.c_str()));
					  TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
					  hPad -> SetTitle(";entries;energy [a.u.]");
						hPad -> Draw();
			 			gPad -> SetGridy();
						h_energy_coincidence_run[label_coin2] = new TH1F ( Form("h_energy_coincidence_%s",label_coin2.c_str()), Form("h_energy_coincidence_%s",label_coin2.c_str()), 1000, 0, 50);
					}
						
					if(c_energy_coincidence_tot[labelTot_coin2]==NULL){
					  c_energy_coincidence_tot[labelTot_coin2] = new TCanvas(Form("c_energy_coincidence_%s",labelTot_coin2.c_str()), Form("c_energy_coincidence_%s",labelTot_coin2.c_str()));
					  TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
					  hPad -> SetTitle(";entries;energy [a.u.]");
						hPad -> Draw();
			 			gPad -> SetGridy();
						h_energy_coincidence_tot[labelTot_coin2] = new TH1F ( Form("h_energy_coincidence_%s",labelTot_coin2.c_str()), Form("h_energy_coincidence_%s",labelTot_coin2.c_str()), 1000, 0, 50);
					}
					if(c_timeRes_coincidence_run[label_coin]==NULL){
					  c_timeRes_coincidence_run[label_coin] = new TCanvas(Form("c_timeRes_coincidence_%s",label_coin.c_str()), Form("c_timeRes_coincidence_%s",label_coin.c_str()));
					  TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
					  hPad -> SetTitle("energy-corrected #Deltat [ps];entries");
						hPad -> Draw();
			 			gPad -> SetGridy();
						h_timeRes_coincidence_run[label_coin] = new TH1F ( Form("h_timeRes_coincidence_%s",label_coin.c_str()), Form("h_timeRes_coincidence_%s",label_coin.c_str()), 200, -4000, 4000);
						h_timeRes_corr_coincidence_run[label_coin] = new TH1F ( Form("h_timeRes_corr_coincidence_%s",label_coin.c_str()), Form("h_timeRes_corr_coincidence_%s",label_coin.c_str()), 200, -4000, 4000);
						h_timeRes_corr_data_coincidence_run[label_coin] = new TH1F ( Form("h_timeRes_corr_data_coincidence_%s",label_coin.c_str()), Form("h_timeRes_corr_data_coincidence_%s",label_coin.c_str()), 200, -4000, 4000);
					}
						
					if(c_timeRes_coincidence_tot[labelTot_coin]==NULL){
					  c_timeRes_coincidence_tot[labelTot_coin] = new TCanvas(Form("c_timeRes_coincidence_%s",labelTot_coin.c_str()), Form("c_timeRes_coincidence_%s",labelTot_coin.c_str()));
					  TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
					  hPad -> SetTitle(";energy-corrected #Deltat [ps];entries");
						hPad -> Draw();
			 			gPad -> SetGridy();
						h_timeRes_coincidence_tot[labelTot_coin] = new TH1F ( Form("h_timeRes_coincidence_%s",labelTot_coin.c_str()), Form("h_timeRes_coincidence_%s",labelTot_coin.c_str()), 200, -4000, 4000);
						h_timeRes_corr_coincidence_tot[labelTot_coin] = new TH1F ( Form("h_timeRes_corr_coincidence_%s",labelTot_coin.c_str()), Form("h_timeRes_corr_coincidence_%s",labelTot_coin.c_str()), 200, -4000, 4000);
						h_timeRes_corr_data_coincidence_tot[labelTot_coin] = new TH1F ( Form("h_timeRes_corr_data_coincidence_%s",labelTot_coin.c_str()), Form("h_timeRes_corr_data_coincidence_%s",labelTot_coin.c_str()), 200, -4000, 4000);
					}
			
					if(c_timeRes_coincidence_run[label_coin2]==NULL){
					  c_timeRes_coincidence_run[label_coin2] = new TCanvas(Form("c_timeRes_coincidence_%s",label_coin2.c_str()), Form("c_timeRes_coincidence_%s",label_coin2.c_str()));
					  TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
					  hPad -> SetTitle("energy-corrected #Deltat [ps];entries");
						hPad -> Draw();
			 			gPad -> SetGridy();
						h_timeRes_coincidence_run[label_coin2] = new TH1F ( Form("h_timeRes_coincidence_%s",label_coin2.c_str()), Form("h_timeRes_coincidence_%s",label_coin2.c_str()), 200, -4000, 4000);
						h_timeRes_corr_coincidence_run[label_coin2] = new TH1F ( Form("h_timeRes_corr_coincidence_%s",label_coin2.c_str()), Form("h_timeRes_corr_coincidence_%s",label_coin2.c_str()), 200, -4000, 4000);
						h_timeRes_corr_data_coincidence_run[label_coin2] = new TH1F ( Form("h_timeRes_corr_data_coincidence_%s",label_coin2.c_str()), Form("h_timeRes_corr_data_coincidence_%s",label_coin2.c_str()), 200, -4000, 4000);
					}
						
					if(c_timeRes_coincidence_tot[labelTot_coin2]==NULL){
					  c_timeRes_coincidence_tot[labelTot_coin2] = new TCanvas(Form("c_timeRes_coincidence_%s",labelTot_coin2.c_str()), Form("c_timeRes_coincidence_%s",labelTot_coin2.c_str()));
					  TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
					  hPad -> SetTitle(";energy-corrected #Deltat [ps];entries");
						hPad -> Draw();
			 			gPad -> SetGridy();
						h_timeRes_coincidence_tot[labelTot_coin2] = new TH1F ( Form("h_timeRes_coincidence_%s",labelTot_coin2.c_str()), Form("h_timeRes_coincidence_%s",labelTot_coin2.c_str()), 200, -4000, 4000);
						h_timeRes_corr_coincidence_tot[labelTot_coin2] = new TH1F ( Form("h_timeRes_corr_coincidence_%s",labelTot_coin2.c_str()), Form("h_timeRes_corr_coincidence_%s",labelTot_coin2.c_str()), 200, -4000, 4000);
						h_timeRes_corr_data_coincidence_tot[labelTot_coin2] = new TH1F ( Form("h_timeRes_corr_data_coincidence_%s",labelTot_coin2.c_str()), Form("h_timeRes_corr_data_coincidence_%s",labelTot_coin2.c_str()), 200, -4000, 4000);
					}
				}
			}
			for(int iBin = 8; iBin < 10; iBin ++){
				if ( c_TProfile_energyRatio_corr[labelTot_c][iBin] == NULL){
					c_TProfile_energyRatio_corr[labelTot_c][iBin] = new TCanvas( Form("c_deltaT_vs_energyRatio_%s_energyBin%d", labelTot_c.c_str(), iBin), Form("c_deltaT_vs_energyRatio_%s_energyBin%d", labelTot_c.c_str(),iBin));
					TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
					hPad -> Draw();
					gPad -> SetGridy();
					p1_TProfile_energyRatio_corr[labelTot_c][iBin] = new TProfile(Form("p1_deltaT_vs_energyRatio_%s_energyBin%d",labelTot_c.c_str(), iBin),Form("p1_deltaT_vs_energyRatio_%s_energyBin%d",labelTot_c.c_str(), iBin),100, 0, 2, -4000,4000);
					h1_energyRatio[labelTot_c][iBin] = new TH1F(Form("h1_energyRatio_%s_energyBin%d",labelTot_c.c_str(),iBin),"",1000,0.,5.);
				}
			}
		}
	}

	std::map<std::string, std::map < int ,bool>> acceptEvent;
		
	for (auto mapIt : dataTrees){
		std::cout << mapIt.first << std::endl;
		int nEvents = mapIt.second -> GetEntries();

	  float qfine[256];
	  float energy[256];
	  long long time[256];
	  mapIt.second -> SetBranchStatus("*",0);
	  mapIt.second -> SetBranchStatus("qfine",       1); mapIt.second -> SetBranchAddress("qfine",        qfine);
	  mapIt.second -> SetBranchStatus("energy",      1); mapIt.second -> SetBranchAddress("energy",       energy);
	  mapIt.second -> SetBranchStatus("time",        1); mapIt.second -> SetBranchAddress("time",         time);

		for (int entry = 0; entry < nEvents; entry++){
    	mapIt.second -> GetEntry(entry);
			if (entry%500000 == 0) std::cout << "Analyzing entry " << entry << "/" << nEvents <<  std::endl;
			acceptEvent[mapIt.first][entry] = false;
			for( int iBar : barMap[mapIt.first]){
				
				std::string label_c = Form( "bar%02dL-R_%s", iBar, mapIt.first.c_str());
				std::string label_en = Form( "bar%02dL-R_%s", iBar, mapIt.first.c_str());
				std::vector<std::string> tokens = GetTokens(mapIt.first,'_');
				std::string labelTot_c = Form( "bar%02dL-R_%s_%s", iBar, tokens[0].c_str(), tokens[1].c_str());

				int chL(chsL[iBar]);
				int chR(chsR[iBar]);
				          
				float qfineL(qfine[chL]);
				float qfineR(qfine[chR]);

				if( qfineL < 13. || qfineR < 13. ) continue;
				if( qfineL > 500. || qfineR > 500. ) continue;
				  
				float energyL(energy[chL]);
				float energyR(energy[chR]);
				long long timeL(time[chL]);
				long long timeR(time[chR]);
				float avEnergy = 0.5*(energyL+energyR);
				//energy spectrum
				h_energy_bar_run[label_c] -> Fill( avEnergy);
				h_energy_bar_tot[labelTot_c] -> Fill( avEnergy);
				
				if( avEnergy > energyRanges[label_en][0] && 	avEnergy < energyRanges[label_en][1] ){	
	
					//1173 peak
					h_energy_1173_bar_run[label_c] -> Fill( avEnergy);
					h_energy_1173_bar_tot[labelTot_c] -> Fill( avEnergy);

					//coincidence 1173 vs 1332 (no energy in other bars)
					int coin_chL[barMap[mapIt.first].size()];
					int coin_chR[barMap[mapIt.first].size()];
					float coin_qfineL[barMap[mapIt.first].size()];
					float coin_qfineR[barMap[mapIt.first].size()];
					
					int counter = 0;
					for ( int i = 0; i < barMap[mapIt.first].size(); i++){
						int coinBar = barMap[mapIt.first][i];
						coin_chL[i] = chsL[coinBar];
						coin_chR[i] = chsR[coinBar];
						coin_qfineL[i] = qfine[coin_chL[i]];
						coin_qfineR[i] = qfine[coin_chR[i]];
						if ( coin_qfineL[i] < 13. && coin_qfineR[i] < 13.){
							counter ++;
						}			
					}
					if(counter == barMap[mapIt.first].size() -2){
						for ( int i = 0; i < barMap[mapIt.first].size(); i++){
							int coinBar = barMap[mapIt.first][i];
							if (coinBar != iBar && coin_qfineL[i] > 13. && coin_qfineR[i] > 13.){ 
								float coinEnergy = 0.5 *( energy[coin_chL[i]] + energy[coin_chR[i]]);
								if( coinEnergy > energyRanges[label_en][1] && 	coinEnergy < energyRanges[label_en][2] ){
									std::string label_coin =  Form( "1173_bar%02dL-R_vs_1332_bar%02dL-R_%s", iBar, coinBar, mapIt.first.c_str());
									std::string labelTot_coin = Form( "1173_bar%02dL-R_vs_1332_bar%02dL-R_%s_%s", iBar, coinBar, tokens[0].c_str(), tokens[1].c_str());
									h_energy_coincidence_run[label_coin] -> Fill(avEnergy);
									h_energy_coincidence_tot[labelTot_coin] -> Fill(avEnergy);
									//time resolution
									h_timeRes_coincidence_run[label_coin] -> Fill(timeR - timeL);
									h_timeRes_coincidence_tot[labelTot_coin] -> Fill(timeR -timeL);
									acceptEvent[mapIt.first][entry] = true;
									
									float energyRatioCorr = tProfileFitFunc[label_en][9]->Eval(energyR/energyL) - tProfileFitFunc[label_en][9]->Eval(energyRatioFitFunc[label_en][9]->GetParameter(1));
									h_timeRes_corr_coincidence_run[label_coin] -> Fill(timeR -timeL -energyRatioCorr);
									h_timeRes_corr_coincidence_tot[labelTot_coin] -> Fill(timeR -timeL -energyRatioCorr);
									
									h1_energyRatio[labelTot_c][8] -> Fill (energyR/energyL);
									p1_TProfile_energyRatio_corr[labelTot_c][8] -> Fill( energyR/energyL, timeR - timeL);
								}	
							}
						}
					}
				}
				
				if( avEnergy > energyRanges[label_en][1] && 	avEnergy < energyRanges[label_en][2] ){

					//1332 peak
					h_energy_1332_bar_run[label_c] -> Fill( avEnergy);
					h_energy_1332_bar_tot[labelTot_c] -> Fill( avEnergy);
			
					//coincidence 1332 vs 1173 (no energy in other bars)
					int coin_chL[barMap[mapIt.first].size()];
					int coin_chR[barMap[mapIt.first].size()];
					float coin_qfineL[barMap[mapIt.first].size()];
					float coin_qfineR[barMap[mapIt.first].size()];

					int counter = 0;
					for ( int i = 0; i < barMap[mapIt.first].size(); i++){
						int coinBar = barMap[mapIt.first][i];
						coin_chL[i] = chsL[coinBar];
						coin_chR[i] = chsR[coinBar];
						coin_qfineL[i] = qfine[coin_chL[i]];
						coin_qfineR[i] = qfine[coin_chR[i]];
						if ( coin_qfineL[i] < 13. && coin_qfineR[i] < 13.){
							counter ++;
						}			
					}
					if(counter == barMap[mapIt.first].size() -2){
						for ( int i = 0; i < barMap[mapIt.first].size(); i++){
							int coinBar = barMap[mapIt.first][i];
							if (coinBar != iBar && coin_qfineL[i] > 13. && coin_qfineR[i] > 13.){ 
								float coinEnergy = 0.5 *( energy[coin_chL[i]] + energy[coin_chR[i]]);
								if( coinEnergy > energyRanges[label_en][0] && 	coinEnergy < energyRanges[label_en][1] ){
									std::string label_coin =  Form( "1332_bar%02dL-R_vs_1173_bar%02dL-R_%s", iBar, coinBar, mapIt.first.c_str());
									std::string labelTot_coin = Form( "1332_bar%02dL-R_vs_1173_bar%02dL-R_%s_%s", iBar, coinBar, tokens[0].c_str(), tokens[1].c_str());
									h_energy_coincidence_run[label_coin] -> Fill(avEnergy);
									h_energy_coincidence_tot[labelTot_coin] -> Fill(avEnergy);

									//time resolution
									h_timeRes_coincidence_run[label_coin] -> Fill(timeR - timeL);
									h_timeRes_coincidence_tot[labelTot_coin] -> Fill(timeR -timeL);
									acceptEvent[mapIt.first][entry] = true;

									float energyRatioCorr = tProfileFitFunc[label_en][9]->Eval(energyR/energyL) - tProfileFitFunc[label_en][9]->Eval(energyRatioFitFunc[label_en][9]->GetParameter(1));
									h_timeRes_corr_coincidence_run[label_coin] -> Fill(timeR -timeL -energyRatioCorr);
									h_timeRes_corr_coincidence_tot[labelTot_coin] -> Fill(timeR -timeL -energyRatioCorr);
									h1_energyRatio[labelTot_c][9] -> Fill (energyR/energyL);
									p1_TProfile_energyRatio_corr[labelTot_c][9] -> Fill( energyR/energyL, timeR - timeL);
								}	
							}
						}
					}
				}
			}
		}
	}



	std::map < std::string, std::map <int,TF1*>>	p1_TProfile_energyRatio_corr_FitFunc;
	std::map < std::string, std::map <int,TF1*>>	h1_energyRatio_FitFunc;
	
	for (auto mapIt :c_TProfile_energyRatio_corr){
		for(int iBin = 8; iBin < 10; iBin++){
			
			c_TProfile_energyRatio_corr[mapIt.first][iBin] -> cd();
			TProfile *prof = p1_TProfile_energyRatio_corr[mapIt.first][iBin]; 
			p1_TProfile_energyRatio_corr_FitFunc[mapIt.first][iBin] = new TF1(Form("p1_TProfile_energyRatio_corr_FitFunc%s",Form("bin%d_%s", iBin, mapIt.first.c_str())),"pol4",0,10);
			float xmin = p1_TProfile_energyRatio_corr[mapIt.first][iBin] -> GetMean() - 2.5*p1_TProfile_energyRatio_corr[mapIt.first][iBin] -> GetRMS();
    	float xmax = p1_TProfile_energyRatio_corr[mapIt.first][iBin] -> GetMean() + 2.5*p1_TProfile_energyRatio_corr[mapIt.first][iBin] -> GetRMS();
   	  p1_TProfile_energyRatio_corr_FitFunc[mapIt.first][iBin] -> SetRange(xmin, xmax);
    	prof -> Fit(p1_TProfile_energyRatio_corr_FitFunc[mapIt.first][iBin],"QRS+");
			prof -> SetTitle(Form(";energy_{right} / energy_{left};#Deltat [ps]"));
			prof -> Draw("");
			p1_TProfile_energyRatio_corr_FitFunc[mapIt.first][iBin] -> SetLineColor(kRed);
      p1_TProfile_energyRatio_corr_FitFunc[mapIt.first][iBin] -> SetLineWidth(2);
      p1_TProfile_energyRatio_corr_FitFunc[mapIt.first][iBin] -> Draw("same");
		
			outFile -> cd();
			prof -> Write();		
		}
	}
	for (auto mapIt : h1_energyRatio){
		for(int iBin = 8; iBin < 10; iBin++){
			
			TCanvas * canvas = new TCanvas (Form("canvas_energyRatio%s",Form("bin%d_%s", iBin, mapIt.first.c_str())), Form("canvas_energyRatio%s",Form("bin%d_%s", iBin, mapIt.first.c_str())));
			canvas -> cd();
			TH1F* histo = h1_energyRatio[mapIt.first][iBin];
			h1_energyRatio_FitFunc[mapIt.first][iBin] = new TF1(Form("h1_energyRatio_FitFunc%s",Form("bin%d_%s", iBin, mapIt.first.c_str())),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
    	histo -> Fit(h1_energyRatio_FitFunc[mapIt.first][iBin],"QNRS");
			histo -> Fit(h1_energyRatio_FitFunc[mapIt.first][iBin],"QS+","",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());

			outFile -> cd();		
			histo -> Write();
		}
	}



	std::cout << " Second loop over events " << std::endl;
	for (auto mapIt : dataTrees){
		std::cout << mapIt.first << std::endl;
		int nEvents = mapIt.second -> GetEntries();
		
		float qfine[256];
	  float energy[256];
	  long long time[256];
	  
		for (int entry = 0; entry < nEvents; entry++){
    	mapIt.second -> GetEntry(entry);
			if (entry%500000 == 0) std::cout << "Analyzing entry " << entry << "/" << nEvents <<  std::endl;
			if(	!acceptEvent[mapIt.first][entry] ) continue;
			for( int iBar : barMap[mapIt.first]){
				
				std::string label_c = Form( "bar%02dL-R_%s", iBar, mapIt.first.c_str());
				std::vector<std::string> tokens = GetTokens(mapIt.first,'_');
				std::string labelTot_c = Form( "bar%02dL-R_%s_%s", iBar, tokens[0].c_str(), tokens[1].c_str());

				int chL(chsL[iBar]);
				int chR(chsR[iBar]);
				          
				float qfineL(qfine[chL]);
				float qfineR(qfine[chR]);

				if( qfineL < 13. || qfineR < 13. ) continue;
				if( qfineL > 500. || qfineR > 500. ) continue;
				  
				float energyL(energy[chL]);
				float energyR(energy[chR]);
				long long timeL(time[chL]);
				long long timeR(time[chR]);
				float avEnergy = 0.5*(energyL+energyR);
				
				
				if( avEnergy > energyRanges[label_c][0] && 	avEnergy < energyRanges[label_c][1] ){	
					

					//coincidence 1173 vs 1332 (no energy in other bars)
					int coin_chL[barMap[mapIt.first].size()];
					int coin_chR[barMap[mapIt.first].size()];
					float coin_qfineL[barMap[mapIt.first].size()];
					float coin_qfineR[barMap[mapIt.first].size()];
					
					int counter = 0;
					for ( int i = 0; i < barMap[mapIt.first].size(); i++){
						int coinBar = barMap[mapIt.first][i];
						coin_chL[i] = chsL[coinBar];
						coin_chR[i] = chsR[coinBar];
						coin_qfineL[i] = qfine[coin_chL[i]];
						coin_qfineR[i] = qfine[coin_chR[i]];
						if ( coin_qfineL[i] < 13. && coin_qfineR[i] < 13.){
							counter ++;
						}			
					}
					if(counter == barMap[mapIt.first].size() -2){
						for ( int i = 0; i < barMap[mapIt.first].size(); i++){
							int coinBar = barMap[mapIt.first][i];
							if (coinBar != iBar && coin_qfineL[i] > 13. && coin_qfineR[i] > 13.){ 
								float coinEnergy = 0.5 *( energy[coin_chL[i]] + energy[coin_chR[i]]);
								if( coinEnergy > energyRanges[label_c][1] && 	coinEnergy < energyRanges[label_c][2] ){
									std::string label_coin =  Form( "1173_bar%02dL-R_vs_1332_bar%02dL-R_%s", iBar, coinBar, mapIt.first.c_str());
									std::string labelTot_coin = Form( "1173_bar%02dL-R_vs_1332_bar%02dL-R_%s_%s", iBar, coinBar, tokens[0].c_str(), tokens[1].c_str());
									std::string label_TProf = Form("bar%02dL-R_%s_%s",iBar, tokens[0].c_str(), tokens[1].c_str());

									float energyRatioCorr = p1_TProfile_energyRatio_corr_FitFunc[label_TProf][9]->Eval(energyR/energyL) - p1_TProfile_energyRatio_corr_FitFunc[label_TProf][9]->Eval(h1_energyRatio_FitFunc[label_TProf][9]->GetParameter(1));
									h_timeRes_corr_data_coincidence_run[label_coin] -> Fill(timeR -timeL -energyRatioCorr);
									h_timeRes_corr_data_coincidence_tot[labelTot_coin] -> Fill(timeR -timeL -energyRatioCorr);
									
								}	
							}
						}
					}
				}
				
				if( avEnergy > energyRanges[label_c][1] && 	avEnergy < energyRanges[label_c][2] ){
					//coincidence 1332 vs 1173 (no energy in other bars)
					int coin_chL[barMap[mapIt.first].size()];
					int coin_chR[barMap[mapIt.first].size()];
					float coin_qfineL[barMap[mapIt.first].size()];
					float coin_qfineR[barMap[mapIt.first].size()];

					int counter = 0;
					for ( int i = 0; i < barMap[mapIt.first].size(); i++){
						int coinBar = barMap[mapIt.first][i];
						coin_chL[i] = chsL[coinBar];
						coin_chR[i] = chsR[coinBar];
						coin_qfineL[i] = qfine[coin_chL[i]];
						coin_qfineR[i] = qfine[coin_chR[i]];
						if ( coin_qfineL[i] < 13. && coin_qfineR[i] < 13.){
							counter ++;
						}			
					}
					if(counter == barMap[mapIt.first].size() -2){
						for ( int i = 0; i < barMap[mapIt.first].size(); i++){
							int coinBar = barMap[mapIt.first][i];
							if (coinBar != iBar && coin_qfineL[i] > 13. && coin_qfineR[i] > 13.){ 
								float coinEnergy = 0.5 *( energy[coin_chL[i]] + energy[coin_chR[i]]);
								if( coinEnergy > energyRanges[label_c][0] && 	coinEnergy < energyRanges[label_c][1] ){
									std::string label_coin =  Form( "1332_bar%02dL-R_vs_1173_bar%02dL-R_%s", iBar, coinBar, mapIt.first.c_str());
									std::string labelTot_coin = Form( "1332_bar%02dL-R_vs_1173_bar%02dL-R_%s_%s", iBar, coinBar, tokens[0].c_str(), tokens[1].c_str());
									std::string label_TProf = Form("bar%02dL-R_%s_%s",iBar, tokens[0].c_str(), tokens[1].c_str());
									float energyRatioCorr = p1_TProfile_energyRatio_corr_FitFunc[label_TProf][8]->Eval(energyR/energyL) - p1_TProfile_energyRatio_corr_FitFunc[label_TProf][8]->Eval(h1_energyRatio_FitFunc[label_TProf][8]->GetParameter(1));
									h_timeRes_corr_data_coincidence_run[label_coin] -> Fill(timeR -timeL -energyRatioCorr);
									h_timeRes_corr_data_coincidence_tot[labelTot_coin] -> Fill(timeR -timeL -energyRatioCorr);
								}	
							}
						}
					}
				}
			}
		}
	}

	//Draw and save graphs
	for ( auto label : c_energy_bar_run){
		label.second -> cd();
		h_energy_bar_run[label.first] ->  SetTitle(";energy [a.u.];entries");
		h_energy_bar_run[label.first] -> Draw();
		label.second-> Print(Form("%s/coincidenceResults/energy/runs/c_energy_%s.png",plotDir.c_str(),label.first.c_str()));
    label.second-> Print(Form("%s/coincidenceResults/energy/runs/c_energy_%s.pdf",plotDir.c_str(),label.first.c_str()));
		outFile -> cd();
		h_energy_bar_run[label.first] -> Write();	
	}
	for ( auto label : c_energy_bar_tot){
		label.second -> cd();
		h_energy_bar_tot[label.first] ->  SetTitle(";energy [a.u.];entries");
		h_energy_bar_tot[label.first] -> Draw();
		label.second-> Print(Form("%s/coincidenceResults/energy/total/c_energy_%s.png",plotDir.c_str(),label.first.c_str()));
    label.second-> Print(Form("%s/coincidenceResults/energy/total/c_energy_%s.pdf",plotDir.c_str(),label.first.c_str()));
		outFile -> cd();
		h_energy_bar_tot[label.first] -> Write();	
	}
	for ( auto label : c_energy_1173_bar_run){
		label.second -> cd();
		h_energy_1173_bar_run[label.first] -> SetTitle(";energy [a.u.];entries");
		h_energy_1173_bar_run[label.first] -> Draw();
		label.second-> Print(Form("%s/coincidenceResults/energy/runs/c_1173_peak_%s.png",plotDir.c_str(),label.first.c_str()));
    label.second-> Print(Form("%s/coincidenceResults/energy/runs/c_1173_peak_%s.pdf",plotDir.c_str(),label.first.c_str()));
		outFile -> cd();
		h_energy_1173_bar_run[label.first] -> Write();	
	}
	for ( auto label : c_energy_1173_bar_tot){
		label.second -> cd();
		h_energy_1173_bar_tot[label.first] -> SetTitle(";energy [a.u.];entries");
		h_energy_1173_bar_tot[label.first] -> Draw();
		label.second-> Print(Form("%s/coincidenceResults/energy/total/c_1173_peak_%s.png",plotDir.c_str(),label.first.c_str()));
    label.second-> Print(Form("%s/coincidenceResults/energy/total/c_1173_peak_%s.pdf",plotDir.c_str(),label.first.c_str()));
		outFile -> cd();
		h_energy_1173_bar_tot[label.first] -> Write();	
	}
	for ( auto label : c_energy_1332_bar_run){
		label.second -> cd();
		h_energy_1332_bar_run[label.first] -> SetTitle(";energy [a.u.];entries");
		h_energy_1332_bar_run[label.first] -> Draw();
		label.second-> Print(Form("%s/coincidenceResults/energy/runs/c_1332_peak_%s.png",plotDir.c_str(),label.first.c_str()));
    label.second-> Print(Form("%s/coincidenceResults/energy/runs/c_1332_peak_%s.pdf",plotDir.c_str(),label.first.c_str()));
		outFile -> cd();
		h_energy_1332_bar_run[label.first] -> Write();	
	}
	for ( auto label : c_energy_1332_bar_tot){
		label.second -> cd();
		h_energy_1332_bar_tot[label.first] -> SetTitle(";energy [a.u.];entries");
		h_energy_1332_bar_tot[label.first] -> Draw();
		label.second-> Print(Form("%s/coincidenceResults/energy/total/c_1332_peak_%s.png",plotDir.c_str(),label.first.c_str()));
    label.second-> Print(Form("%s/coincidenceResults/energy/total/c_1332_peak_%s.pdf",plotDir.c_str(),label.first.c_str()));
		outFile -> cd();
		h_energy_1332_bar_tot[label.first] -> Write();	
	}
	for ( auto label : c_energy_coincidence_run){
		label.second -> cd();
		h_energy_coincidence_run[label.first] -> SetTitle(";energy [a.u.];entries");
		h_energy_coincidence_run[label.first] -> Draw();
		label.second-> Print(Form("%s/coincidenceResults/energy/runs/c_energy_coincidence_%s.png",plotDir.c_str(),label.first.c_str()));
    label.second-> Print(Form("%s/coincidenceResults/energy/runs/c_energy_coincidence_%s.pdf",plotDir.c_str(),label.first.c_str()));
		outFile -> cd();
		h_energy_coincidence_run[label.first] -> Write();	
	}
	for ( auto label : c_energy_coincidence_tot){
		label.second -> cd();
		h_energy_coincidence_tot[label.first] -> SetTitle(";energy [a.u.];entries");
		//h_energy_coincidence_tot[label.first] -> GetYaxis() -> SetRangeUser ( 0, h_energy_coincidence_tot[label.first] -> GetMaximum() *1.1);
		h_energy_coincidence_tot[label.first] -> Draw();
		label.second-> Print(Form("%s/coincidenceResults/energy/total/c_energy_coincidence_%s.png",plotDir.c_str(),label.first.c_str()));
    label.second-> Print(Form("%s/coincidenceResults/energy/total/c_energy_coincidence_%s.pdf",plotDir.c_str(),label.first.c_str()));
		outFile -> cd();
		h_energy_coincidence_tot[label.first] -> Write();	
	}

	for ( auto label : c_timeRes_coincidence_run){
		label.second -> cd();
		h_timeRes_coincidence_run[label.first] -> SetTitle(";energy-corrected #Deltat [ps];entries");
		h_timeRes_coincidence_run[label.first] -> Draw();
		label.second-> Print(Form("%s/coincidenceResults/timeResolution/runs/c_energy_coincidence_%s.png",plotDir.c_str(),label.first.c_str()));
    label.second-> Print(Form("%s/coincidenceResults/timeResolution/runs/c_energy_coincidence_%s.pdf",plotDir.c_str(),label.first.c_str()));
		outFile -> cd();
		h_timeRes_coincidence_run[label.first] -> Write();	
	}
	
	for ( auto label : c_timeRes_coincidence_tot){
		label.second -> cd();
		//raw time resolution
		h_timeRes_coincidence_tot[label.first] -> SetTitle(";energy-corrected #Deltat [ps];entries");
		h_timeRes_coincidence_tot[label.first] -> SetLineColor(kRed+1);
		float fitXMin = h_timeRes_coincidence_tot[label.first] -> GetBinCenter(h_timeRes_coincidence_tot[label.first]-> GetMaximumBin()) - 200.;
    float fitXMax = h_timeRes_coincidence_tot[label.first] -> GetBinCenter(h_timeRes_coincidence_tot[label.first] -> GetMaximumBin()) + 200.;
    TF1* fitFunc = new TF1(Form("fitFunc_%s",label.first.c_str()),"gaus",fitXMin,fitXMax);
    fitFunc -> SetParameters(1, h_timeRes_coincidence_tot[label.first] -> GetMean(), h_timeRes_coincidence_tot[label.first] -> GetRMS());
    h_timeRes_coincidence_tot[label.first] -> Fit(fitFunc,"QNRSL+","");
    h_timeRes_coincidence_tot[label.first] -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
		h_timeRes_coincidence_tot[label.first] -> Draw();

		float fitXMin2 = fitFunc->GetParameter(1)-fabs(1.5*fitFunc->GetParameter(2));
    float fitXMax2 = fitFunc->GetParameter(1)+fabs(1.5*fitFunc->GetParameter(2));
    TF1* fitFunc2 = new TF1(Form("fitFunc2_energy_%s",label.first.c_str()),"gaus(0)",fitXMin2,fitXMax2);
    fitFunc2 -> SetParameter(0,fitFunc->GetParameter(0));
    fitFunc2 -> SetParameter(1,fitFunc->GetParameter(1));
    fitFunc2 -> SetParameter(2,fitFunc->GetParameter(2));
    h_timeRes_coincidence_tot[label.first] -> Fit(fitFunc2,"QRSL+");
    h_timeRes_coincidence_tot[label.first] -> Fit(fitFunc2,"QRSL+");
        
    fitFunc2 -> SetLineColor(kRed+1);
    fitFunc2 -> SetLineWidth(3);
    fitFunc2 -> Draw("same");
       
    TLatex* latex = new TLatex(0.20,0.85,Form("{#sigma_{raw}^{gaus} = %.0f ps}",fabs(fitFunc2->GetParameter(2))));
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    latex -> SetTextColor(kRed);
    latex -> Draw("same");

		//step 2 corrected time resolution
		h_timeRes_corr_coincidence_tot[label.first] -> SetLineColor(kBlue+1);
		fitXMin = h_timeRes_corr_coincidence_tot[label.first] -> GetBinCenter(h_timeRes_corr_coincidence_tot[label.first]-> GetMaximumBin()) - 200.;
    fitXMax = h_timeRes_corr_coincidence_tot[label.first] -> GetBinCenter(h_timeRes_corr_coincidence_tot[label.first] -> GetMaximumBin()) + 200.;
    TF1* fitFunc3 = new TF1(Form("fitFunc_%s",label.first.c_str()),"gaus",fitXMin,fitXMax);
    fitFunc3 -> SetParameters(1, h_timeRes_corr_coincidence_tot[label.first] -> GetMean(), h_timeRes_corr_coincidence_tot[label.first] -> GetRMS());
    h_timeRes_corr_coincidence_tot[label.first] -> Fit(fitFunc3,"QNRSL+","");
    h_timeRes_corr_coincidence_tot[label.first] -> Fit(fitFunc3,"QNRSL+","",fitFunc3->GetParameter(1)-fitFunc3->GetParameter(2),fitFunc3->GetParameter(1)+fitFunc3->GetParameter(2));
		h_timeRes_corr_coincidence_tot[label.first] -> Draw("same");

		fitXMin2 = fitFunc3->GetParameter(1)-fabs(1.5*fitFunc3->GetParameter(2));
    fitXMax2 = fitFunc3->GetParameter(1)+fabs(1.5*fitFunc3->GetParameter(2));
    TF1* fitFunc4 = new TF1(Form("fitFunc2_energyCorr_%s",label.first.c_str()),"gaus(0)",fitXMin2,fitXMax2);
    fitFunc4 -> SetParameter(0,fitFunc3->GetParameter(0));
    fitFunc4 -> SetParameter(1,fitFunc3->GetParameter(1));
    fitFunc4 -> SetParameter(2,fitFunc3->GetParameter(2));
    h_timeRes_corr_coincidence_tot[label.first] -> Fit(fitFunc4,"QRSL+");
    h_timeRes_corr_coincidence_tot[label.first] -> Fit(fitFunc4,"QRSL+");
        
    fitFunc4 -> SetLineColor(kBlue+1);
    fitFunc4 -> SetLineWidth(3);
    fitFunc4 -> Draw("same");
       
    latex = new TLatex(0.55,0.85,Form("{#sigma_{corr. step2}^{gaus} = %.0f ps}",fabs(fitFunc4->GetParameter(2))));
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    latex -> SetTextColor(kBlue);
    latex -> Draw("same");
		
		
		

		//data energy ratio corrected time resolution
		h_timeRes_corr_data_coincidence_tot[label.first] -> SetLineColor(kBlue+1);
		fitXMin = h_timeRes_corr_data_coincidence_tot[label.first] -> GetBinCenter(h_timeRes_corr_data_coincidence_tot[label.first]-> GetMaximumBin()) - 200.;
    fitXMax = h_timeRes_corr_data_coincidence_tot[label.first] -> GetBinCenter(h_timeRes_corr_data_coincidence_tot[label.first] -> GetMaximumBin()) + 200.;
    TF1* fitFunc5 = new TF1(Form("fitFunc_%s",label.first.c_str()),"gaus",fitXMin,fitXMax);
    fitFunc5 -> SetParameters(1, h_timeRes_corr_data_coincidence_tot[label.first] -> GetMean(), h_timeRes_corr_data_coincidence_tot[label.first] -> GetRMS());
    h_timeRes_corr_data_coincidence_tot[label.first] -> Fit(fitFunc5,"QNRSL+","");
    h_timeRes_corr_data_coincidence_tot[label.first] -> Fit(fitFunc5,"QNRSL+","",fitFunc5->GetParameter(1)-fitFunc5->GetParameter(2),fitFunc5->GetParameter(1)+fitFunc5->GetParameter(2));
		h_timeRes_corr_data_coincidence_tot[label.first] -> Draw("same");

		fitXMin2 = fitFunc5->GetParameter(1)-fabs(1.5*fitFunc5->GetParameter(2));
    fitXMax2 = fitFunc5->GetParameter(1)+fabs(1.5*fitFunc5->GetParameter(2));
    TF1* fitFunc6 = new TF1(Form("fitFunc2_energyCorr_%s",label.first.c_str()),"gaus(0)",fitXMin2,fitXMax2);
    fitFunc6 -> SetParameter(0,fitFunc5->GetParameter(0));
    fitFunc6 -> SetParameter(1,fitFunc5->GetParameter(1));
    fitFunc6 -> SetParameter(2,fitFunc5->GetParameter(2));
    h_timeRes_corr_data_coincidence_tot[label.first] -> Fit(fitFunc6,"QRSL+");
    h_timeRes_corr_data_coincidence_tot[label.first] -> Fit(fitFunc6,"QRSL+");
        
    fitFunc6 -> SetLineColor(kGreen+3);
    fitFunc6 -> SetLineWidth(3);
    fitFunc6 -> Draw("same");
       
    latex = new TLatex(0.55,0.70,Form("{#sigma_{corr. data}^{gaus} = %.0f ps}",fabs(fitFunc6->GetParameter(2))));
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    latex -> SetTextColor(kGreen+3);
    latex -> Draw("same");

		
		label.second-> Print(Form("%s/coincidenceResults/timeResolution/total/c_energy_coincidence_%s.png",plotDir.c_str(),label.first.c_str()));
    label.second-> Print(Form("%s/coincidenceResults/timeResolution/total/c_energy_coincidence_%s.pdf",plotDir.c_str(),label.first.c_str()));
		outFile -> cd();
		h_timeRes_coincidence_tot[label.first] -> Write();	
		h_timeRes_corr_coincidence_tot[label.first] -> Write();



		std::vector<std::string> tokens = GetTokens(label.first,'_');
		int iBin = 0;
		std::string label_graph = Form("%s_%s_%s", tokens[1].c_str(), tokens[5].c_str(), tokens[6].c_str());
		if( atof(tokens[0].c_str()) == 1173){
			iBin = 8;
		}	
		if( atof(tokens[0].c_str()) == 1332){
			iBin = 9;
		}	
		TGraphErrors* graph = new TGraphErrors();		
		graph -> SetPoint( 0, gr_timeRes_summary[label_graph][iBin].size() +1, fitFunc4->GetParameter(2));
		graph -> SetPointError( 0, 0, 0);
		graph -> SetMarkerColor (2* (gr_timeRes_summary[label_graph][iBin].size() +1));
		graph -> SetMarkerStyle(26); 
		gr_timeRes_summary[label_graph][iBin].push_back( (TGraphErrors*) graph);
		l_timeRes_summary[label_graph][iBin] -> AddEntry( graph, Form("%s_coincidence_step2Corr",tokens[4].c_str()), "p");

		TGraphErrors* graph2 = new TGraphErrors();
		graph2 -> SetPoint( 0, gr_timeRes_summary[label_graph][iBin].size() +1, fitFunc6->GetParameter(2));
		graph2 -> SetPointError( 0, 0, 0);
		graph2 -> SetMarkerColor (2*(gr_timeRes_summary[label_graph][iBin].size() +1)-1);
		graph2 -> SetMarkerStyle(24); 
		gr_timeRes_summary[label_graph][iBin].push_back( (TGraphErrors*) graph2);
		l_timeRes_summary[label_graph][iBin] -> AddEntry( graph2, Form("%s_coincidence_dataCorr",tokens[4].c_str()), "p");
	
	}

	for( auto label: c_step2_TProfile_run){
		label.second -> cd();
		label.second -> Print(Form("%s/coincidenceResults/timeResolution/ratioEnergyCorrection/c_deltaT_vs_energyRatio_%s.png",plotDir.c_str(),label.first.c_str()));
   	label.second -> Print(Form("%s/coincidenceResults/timeResolution/ratioEnergyCorrection/c_deltaT_vs_energyRatio_%s.pdf",plotDir.c_str(),label.first.c_str()));
	}
		
	for( auto label: c_TProfile_energyRatio_corr){
		for( int iBin = 8; iBin < 10; iBin++){
		
			label.second[iBin] -> cd();
			label.second[iBin] -> Print(Form("%s/coincidenceResults/timeResolution/ratioEnergyCorrection/c_deltaT_vs_energyRatio_%s_energyBin%d_coincidenceData.png",plotDir.c_str(),label.first.c_str(),iBin));
	   	label.second[iBin] -> Print(Form("%s/coincidenceResults/timeResolution/ratioEnergyCorrection/c_deltaT_vs_energyRatio_%s_energyBin%d_coincidenceData.pdf",plotDir.c_str(),label.first.c_str(),iBin));
		}
	}


	for ( auto label: c_timeRes_summary){
		for( int iBin = 8; iBin < 10; iBin++){
			label.second[iBin] -> cd();
			gr_timeRes_summary[label.first][iBin][gr_timeRes_summary[label.first][iBin].size()-1] -> SetTitle(";;energy-corrected #Deltat [ps]");
			gr_timeRes_summary[label.first][iBin][gr_timeRes_summary[label.first][iBin].size()-1] -> GetHistogram() -> GetXaxis() -> SetLimits(0,gr_timeRes_summary[label.first][iBin].size() + 10);
			gr_timeRes_summary[label.first][iBin][gr_timeRes_summary[label.first][iBin].size()-1] -> GetYaxis() -> SetRangeUser(0, 600);
			gr_timeRes_summary[label.first][iBin][gr_timeRes_summary[label.first][iBin].size()-1] -> Draw("AP");
			for ( int i=0; i < gr_timeRes_summary[label.first][iBin].size()-1; i++){
				gr_timeRes_summary[label.first][iBin][i] -> Draw("psame");
			}
			l_timeRes_summary[label.first][iBin]-> Draw("same");
			label.second[iBin] -> Print(Form("%s/coincidenceResults/timeResolution/timeResolutionComparison/c_timeResolutioComparison_%s_energyBin%d.png",plotDir.c_str(),label.first.c_str(),iBin));
	   	label.second[iBin] -> Print(Form("%s/coincidenceResults/timeResolution/timeResolutionComparison/c_timeResolutioComparison_%s_energyBin%d.pdf",plotDir.c_str(),label.first.c_str(),iBin));
		}
	}
}


