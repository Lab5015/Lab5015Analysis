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



int main(int argc, char** argv)
{
  setTDRStyle();


  if( argc < 2 )
  {
    std::cout << ">>> moduleCharacterization_step3::usage:   " << argv[0] << " configFile.cfg" << std::endl;
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
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDirStep3");
  //system(Form("rm -r %s", plotDir.c_str()));
  system(Form("mkdir -p %s",plotDir.c_str()));
  system(Form("mkdir -p %s/summaryPlots/",plotDir.c_str()));
  system(Form("mkdir -p %s/summaryPlots/tot/",plotDir.c_str()));
  system(Form("mkdir -p %s/summaryPlots/energy/",plotDir.c_str()));
  system(Form("mkdir -p %s/summaryPlots/timeResolution/",plotDir.c_str()));

  //copia il file .php che si trova in una cartella fuori plotDir in plotDir definita nel config file
  system(Form("cp %s/../index.php %s",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/summaryPlots/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/summaryPlots/tot/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/summaryPlots/energy/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/summaryPlots/timeResolution/",plotDir.c_str(),plotDir.c_str()));


  //--- open files and make the tree chain
  //--- Prima di eseguire cambia il numero del run nel config file
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  std::string fileBaseName1 = opts.GetOpt<std::string>("Input.fileBaseName1");
  std::string fileBaseName2 = opts.GetOpt<std::string>("Input.fileBaseName2");
  std::string runs = opts.GetOpt<std::string>("Input.runs");


  //Define variables
  std::vector<std::string> listStep1;
  std::vector<std::string> listStep2;

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
      fileNameStep1 = Form("%s%s1_%s%04d.root",inputDir.c_str(),fileBaseName1.c_str(),fileBaseName2.c_str(),run);
      fileNameStep2 = Form("%s%s2_%s%04d.root",inputDir.c_str(),fileBaseName1.c_str(),fileBaseName2.c_str(),run);
      std::cout << ">>> Adding file " << fileNameStep1 << std::endl;
      std::cout << ">>> Adding file " << fileNameStep2 << std::endl;

      listStep1.push_back(fileNameStep1);
      listStep2.push_back(fileNameStep2);      
    }
  }
  
  //--- define output root file
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileNameStep3");
  TFile* outFile = TFile::Open(outFileName.c_str(),"RECREATE");
  

  //--- get plot settings
  TCanvas* c;
  TCanvas* c1;
  TLatex* latex;
  TH1F* histo;
  //define maps
  std::vector<std::string> LRLabels;
  LRLabels.push_back("L");
  LRLabels.push_back("R");
  
  float refTh = opts.GetOpt<float>("Plots.refTh");
  
  std::map<std::string,TTree*> trees;
  std::map<std::string,TTree*> trees2;
  std::map<std::string,int> VovLabels;
  std::map<std::string,int> thLabels;
  std::map<std::string,int> EnBinLabels;
  std::vector<std::string> stepLabels;
  std::map<std::string,float> map_Vovs;
  std::map<std::string,float> map_ths;
  std::map<std::string,float> map_EnBin; 

  std::vector<float> vec_th1;
 
  //define TGraph
  std::map<std::string,TGraphErrors*> g_tot_vs_th;
  std::map<std::string,TGraphErrors*> g_tot_vs_Vov;
  std::map<std::string,TGraphErrors*> g_tot_vs_bar;

  for(auto file: listStep1){
	  
	  //--- open files

	  
	  TFile* inFile = TFile::Open(file.c_str(),"READ");
	  
	  TList* list = inFile -> GetListOfKeys();
	  TIter next(list);
	  TObject* object = 0;
	  while( (object = next()) )
	  {
	    std::string name(object->GetName());
	    std::vector<std::string> tokens = GetTokens(name,'_');
	    std::size_t found;
	    
	    found = name.find("h1_energy");
	    if( found!=std::string::npos )
	    {
	     //Vov e th
	      std::string stepLabel = tokens[3]+"_"+tokens[4];
	      VovLabels[tokens[3]] += 1;
	      thLabels[tokens[4]] += 1;
	      stepLabels.push_back(stepLabel);
	      std::string string_Vov = tokens[3];
	      string_Vov.erase(0,3);
	      map_Vovs[stepLabel] = atof(string_Vov.c_str());
	      std::string string_th = tokens[4];
	      string_th.erase(0,2);
	      map_ths[stepLabel] = atof(string_th.c_str());
	      float vth1 = map_ths[stepLabel];
              vec_th1.push_back(vth1);
	    }
	  }
	  std::sort(stepLabels.begin(),stepLabels.end());
	  stepLabels.erase(std::unique(stepLabels.begin(),stepLabels.end()),stepLabels.end());
  	  std::sort(vec_th1.begin(),vec_th1.end());
          vec_th1.erase(std::unique(vec_th1.begin(),vec_th1.end()),vec_th1.end());  



	  
	  
	  for(auto stepLabel : stepLabels)
	  {
	    float Vov = map_Vovs[stepLabel];
	    float vth1 = map_ths[stepLabel];
	    std::string VovLabel(Form("Vov%.1f",Vov));
	    std::string thLabel(Form("th%02.0f",vth1));
	    
	    
	    //--------------------------------------------------------
	    
	    //loop sulle barre 
	    for(int iBar = 0; iBar < 16; ++iBar)
	    {
	      int index( (10000*int(Vov*100.)) + (100*vth1) + iBar );
	      
	      for(auto LRLabel : LRLabels )
	      {

		//label histo
		std::string label(Form("bar%02d%s_%s",iBar,LRLabel.c_str(),stepLabel.c_str()));
		
		histo = (TH1F*)( inFile->Get(Form("h1_tot_%s",label.c_str())) );
		if( !histo ) continue;
		if( histo->GetEntries() < 100 ) continue;
		

		latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d%s}{V_{OV} = %.1f V, th. = %d DAC}",iBar,LRLabel.c_str(),Vov,int(vth1)));
		latex -> SetNDC();
		latex -> SetTextFont(42);
		latex -> SetTextSize(0.04);
		latex -> SetTextColor(kRed);
		
		c = new TCanvas(Form("c_tot_%s",label.c_str()),Form("c_tot_%s",label.c_str()));
		
		histo = (TH1F*)( inFile->Get(Form("h1_tot_%s",label.c_str())) );
		histo -> SetTitle(";ToT [ns];entries");
		histo -> SetLineColor(kRed);
		histo -> Draw();
		latex -> Draw("same");
		float max1 = FindXMaximum(histo, 0., 500.);
		histo -> GetXaxis() -> SetRangeUser(0.25*max1,2.*max1);
		TF1* fitFunc1 = new TF1("fitFunc1","gaus",max1-0.05*max1,max1+0.05*max1);
		//TF1* fitFunc1 = new TF1("fitFunc1","gaus",0,500);
		histo -> Fit(fitFunc1,"QNRS+");
		histo -> Fit(fitFunc1,"QNS+", " ", fitFunc1->GetParameter(1)-fitFunc1->GetParameter(2), fitFunc1->GetParameter(1)+fitFunc1->GetParameter(2));
		fitFunc1 -> SetLineColor(kBlack);
		fitFunc1 -> SetLineWidth(3);
		fitFunc1 -> Draw("same");
		//outFile -> cd();     
		//histo -> Write();
		
		/*c -> Print(Form("%s/summaryPlots/tot/c_tot__%s.png",plotDir.c_str(),label.c_str()));
		c -> Print(Form("%s/summaryPlots/tot/c_tot__%s.pdf",plotDir.c_str(),label.c_str()));
		delete c;*/

		
		if(LRLabel=="L"){
		
		  //g_tot_vs_th
		
		  if( g_tot_vs_th[Form("bar%02dL_%s",iBar,VovLabel.c_str())] == NULL ){
		    g_tot_vs_th[Form("bar%02dL_%s",iBar,VovLabel.c_str())] = new TGraphErrors();
		  }
		  g_tot_vs_th[Form("bar%02dL_%s",iBar,VovLabel.c_str())] -> SetPoint(g_tot_vs_th[Form("bar%02dL_%s",iBar,VovLabel.c_str())]->GetN(),vth1,fitFunc1->GetMaximumX());
		  g_tot_vs_th[Form("bar%02dL_%s",iBar,VovLabel.c_str())] -> SetPointError(g_tot_vs_th[Form("bar%02dL_%s",iBar,VovLabel.c_str())]->GetN()-1,0,0);

		  //g_tot_vs_Vov

		  if( g_tot_vs_Vov[Form("bar%02dL_%s",iBar,thLabel.c_str())] == NULL ){
		    g_tot_vs_Vov[Form("bar%02dL_%s",iBar,thLabel.c_str())] = new TGraphErrors();
		  }
		  g_tot_vs_Vov[Form("bar%02dL_%s",iBar,thLabel.c_str())] -> SetPoint(g_tot_vs_Vov[Form("bar%02dL_%s",iBar,thLabel.c_str())]->GetN(),Vov,fitFunc1->GetMaximumX());
		  g_tot_vs_Vov[Form("bar%02dL_%s",iBar,thLabel.c_str())] -> SetPointError(g_tot_vs_Vov[Form("bar%02dL_%s",iBar,thLabel.c_str())]->GetN()-1,0,0);

		  //g_tot_vs_bar

		  if( g_tot_vs_bar[Form("L_%s_%s",VovLabel.c_str(),thLabel.c_str())] == NULL ){
		    g_tot_vs_bar[Form("L_%s_%s",VovLabel.c_str(),thLabel.c_str())] = new TGraphErrors();
		  }
		  g_tot_vs_bar[Form("L_%s_%s",VovLabel.c_str(),thLabel.c_str())] -> SetPoint(g_tot_vs_bar[Form("L_%s_%s",VovLabel.c_str(),thLabel.c_str())]->GetN(),iBar,fitFunc1->GetMaximumX());
		  g_tot_vs_bar[Form("L_%s_%s",VovLabel.c_str(),thLabel.c_str())] -> SetPointError(g_tot_vs_bar[Form("L_%s_%s",VovLabel.c_str(),thLabel.c_str())]->GetN()-1,0,0);

		}
		if(LRLabel=="R"){

		  //g_tot_vs_th

		  if( g_tot_vs_th[Form("bar%02dR_%s",iBar,VovLabel.c_str())] == NULL ){	
		    g_tot_vs_th[Form("bar%02dR_%s",iBar,VovLabel.c_str())] = new TGraphErrors();
		  }
		  g_tot_vs_th[Form("bar%02dR_%s",iBar,VovLabel.c_str())] -> SetPoint(g_tot_vs_th[Form("bar%02dR_%s",iBar,VovLabel.c_str())]->GetN(),vth1,fitFunc1->GetMaximumX());
		  g_tot_vs_th[Form("bar%02dR_%s",iBar,VovLabel.c_str())] -> SetPointError(g_tot_vs_th[Form("bar%02dR_%s",iBar,VovLabel.c_str())]->GetN()-1,0,0);	

		  //g_tot_vs_Vov

		  if( g_tot_vs_Vov[Form("bar%02dR_%s",iBar,thLabel.c_str())] == NULL ){
		    g_tot_vs_Vov[Form("bar%02dR_%s",iBar,thLabel.c_str())] = new TGraphErrors();
		  }
		  g_tot_vs_Vov[Form("bar%02dR_%s",iBar,thLabel.c_str())] -> SetPoint(g_tot_vs_Vov[Form("bar%02dR_%s",iBar,thLabel.c_str())]->GetN(),Vov,fitFunc1->GetMaximumX());
		  g_tot_vs_Vov[Form("bar%02dR_%s",iBar,thLabel.c_str())] -> SetPointError(g_tot_vs_Vov[Form("bar%02dR_%s",iBar,thLabel.c_str())]->GetN()-1,0,0);

		  //g_tot_vs_bar

		  if( g_tot_vs_bar[Form("R_%s_%s",VovLabel.c_str(),thLabel.c_str())] == NULL ){
		    g_tot_vs_bar[Form("R_%s_%s",VovLabel.c_str(),thLabel.c_str())] = new TGraphErrors();
		  }
		  g_tot_vs_bar[Form("R_%s_%s",VovLabel.c_str(),thLabel.c_str())] -> SetPoint(g_tot_vs_bar[Form("R_%s_%s",VovLabel.c_str(),thLabel.c_str())]->GetN(),iBar,fitFunc1->GetMaximumX());
		  g_tot_vs_bar[Form("R_%s_%s",VovLabel.c_str(),thLabel.c_str())] -> SetPointError(g_tot_vs_bar[Form("R_%s_%s",VovLabel.c_str(),thLabel.c_str())]->GetN()-1,0,0);
		}  
	      }
	    }
	  }
  }
   
    
	  //------------------
	  //--- draw summary plots

	  for(int iBar=0; iBar<16; iBar++){
	    
	    c1 = new TCanvas(Form("c_tot_vs_th_bar%02d",iBar),Form("c_tot_vs_th_bar%02d",iBar));
	    
	    TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,500.) );
	    hPad -> SetTitle(";threshold [DAC];ToT [ns]");
	    hPad -> Draw();
	    gPad -> SetGridy();

	    int iter = 0;
	    for(auto mapIt : VovLabels)
	    {
	      std::string label1(Form("bar%02dL_%s",iBar,mapIt.first.c_str()));
	      std::string label2(Form("bar%02dR_%s",iBar,mapIt.first.c_str()));

	      if(! g_tot_vs_th[Form("bar%02dL_%s",iBar,mapIt.first.c_str())]){
		std::cout<<"NON ESISTE barL "<<iBar<<"map "<<mapIt.first.c_str()<<std::endl;
		continue;
	      }

	      if(! g_tot_vs_th[Form("bar%02dR_%s",iBar,mapIt.first.c_str())]){ 
		std::cout<<"NON ESISTE barR "<<iBar<<"map "<<mapIt.first.c_str()<<std::endl;
		continue;
	      }

	      TGraphErrors* g_totL = g_tot_vs_th[label1];
	      TGraphErrors* g_totR = g_tot_vs_th[label2];
	       
	      g_totL -> SetLineColor(1+iter);
	      g_totL -> SetMarkerColor(1+iter);
	      g_totL -> SetMarkerStyle(20);
	      g_totL -> Draw("PL,same");
	      outFile -> cd();
	      g_totL -> Write(Form("g_tot_vs_th_bar%02dL_%s",iBar,mapIt.first.c_str()));
	      
	      g_totR -> SetLineColor(1+iter);
	      g_totR -> SetMarkerColor(1+iter);
	      g_totR -> SetMarkerStyle(25);
	      g_totR -> Draw("PL,same");
	      outFile -> cd();
	      g_totR -> Write(Form("g_tot_vs_th_bar%02dR_%s",iBar,mapIt.first.c_str()));
	      
	      latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kBlack+iter);
	      latex -> Draw("same");

	      ++iter;
	      //std::cout<<"iter "<<iter<<std::endl;
	    }
		  
	    c1 -> Print(Form("%s/summaryPlots/tot/c_tot_vs_th_bar%02d.png",plotDir.c_str(),iBar));
	    c1 -> Print(Form("%s/summaryPlots/tot/c_tot_vs_th_bar%02d.pdf",plotDir.c_str(),iBar));



	    c1 = new TCanvas(Form("c_tot_vs_Vov_bar%02d",iBar),Form("c_tot_vs_Vov_bar%02d",iBar));
	    
	    hPad = (TH1F*)( gPad->DrawFrame(0.,0.,10.,500.) );
	    hPad -> SetTitle(";V_{ov} [V];ToT [ns]");
	    hPad -> Draw();
	    gPad -> SetGridy();

	    iter = 0;
	    for(auto mapIt : thLabels)
	    {
	      std::string label1(Form("bar%02dL_%s",iBar,mapIt.first.c_str()));
	      std::string label2(Form("bar%02dR_%s",iBar,mapIt.first.c_str()));

	      if(! g_tot_vs_Vov[Form("bar%02dL_%s",iBar,mapIt.first.c_str())]){
		std::cout<<"NON ESISTE barL "<<iBar<<"map "<<mapIt.first.c_str()<<std::endl;
		continue;
	      }

	      if(! g_tot_vs_Vov[Form("bar%02dR_%s",iBar,mapIt.first.c_str())]){ 
		std::cout<<"NON ESISTE barR "<<iBar<<"map "<<mapIt.first.c_str()<<std::endl;
		continue;
	      }

	      TGraphErrors* g_totL = g_tot_vs_Vov[label1];
	      TGraphErrors* g_totR = g_tot_vs_Vov[label2];
	       
	      g_totL -> SetLineColor(1+iter);
	      g_totL -> SetMarkerColor(1+iter);
	      g_totL -> SetMarkerStyle(20);
	      g_totL -> Draw("PL,same");
	      outFile -> cd();
	      g_totL -> Write(Form("g_tot_vs_Vov_bar%02dL_%s",iBar,mapIt.first.c_str()));
	      
	      g_totR -> SetLineColor(1+iter);
	      g_totR -> SetMarkerColor(1+iter);
	      g_totR -> SetMarkerStyle(25);
	      g_totR -> Draw("PL,same");
	      outFile -> cd();
	      g_totR -> Write(Form("g_tot_vs_Vov_bar%02dR_%s",iBar,mapIt.first.c_str()));
	      
	      latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kBlack+iter);
	      latex -> Draw("same");

	      ++iter;
	      //std::cout<<"iter "<<iter<<std::endl;
	    }
		  
	    c1 -> Print(Form("%s/summaryPlots/tot/c_tot_vs_Vov_bar%02d.png",plotDir.c_str(),iBar));
	    c1 -> Print(Form("%s/summaryPlots/tot/c_tot_vs_Vov_bar%02d.pdf",plotDir.c_str(),iBar));

	  }


	  for(auto mapIt1 : thLabels){
	    //std::cout<<"ref"<<refTh<<std::endl;
            //std::cout<<"quiiii"<<Form("th%d",int(refTh))<<std::endl;
	    if( mapIt1.first != Form("th%d",int(refTh))) continue;
	    //std::cout<<"quiiiiiii"<<std::endl;
	    for(auto mapIt2 : VovLabels){
	     
	      c1 = new TCanvas(Form("c_tot_vs_bar_%s_%s",mapIt2.first.c_str(),mapIt1.first.c_str()),Form("c_tot_vs_bar_%s_%s",mapIt2.first.c_str(),mapIt1.first.c_str()));
	    
	      TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,200.,17.,300.) );
	      hPad -> SetTitle(";ID bar;ToT [ns]");
	      hPad -> Draw();
	      gPad -> SetGridy();        
		
	      std::string label1(Form("L_%s_%s",mapIt2.first.c_str(),mapIt1.first.c_str()));
	      std::string label2(Form("R_%s_%s",mapIt2.first.c_str(),mapIt1.first.c_str()));

	      if( g_tot_vs_bar[Form("L_%s_%s",mapIt2.first.c_str(),mapIt1.first.c_str())]){
		TGraphErrors* g_totL = g_tot_vs_bar[label1];
		if(g_totL->GetN()>0){
		  g_totL -> SetLineColor(kBlue);
	          g_totL -> SetMarkerColor(kBlue);
		  g_totL -> SetMarkerStyle(20);
		  g_totL -> Draw("P");
		  outFile -> cd();
	          g_totL -> Write(Form("g_tot_vs_barL_%s_%s",mapIt2.first.c_str(),mapIt1.first.c_str()));

                  /*std::string label_prova = "Blue: L	Red: R";
		  latex = new TLatex(0.55,0.85-0.04,label_prova.c_str());
	          latex -> SetNDC();
	          latex -> SetTextFont(42);
	          latex -> SetTextSize(0.04);
	      	  latex -> SetTextColor(kBlack);
	      	  latex -> Draw("same");*/
		  TLegend* legenda = new TLegend(0.55,0.80,0.70,0.85);
		  legenda->AddEntry(g_totL, "Left", "P");
		  legenda -> SetBorderSize(0);
		  legenda->Draw();


	      	  c1 -> Print(Form("%s/summaryPlots/tot/c_tot_vs_bar_%s_%s.png",plotDir.c_str(),mapIt2.first.c_str(),mapIt1.first.c_str()));
	      	  c1 -> Print(Form("%s/summaryPlots/tot/c_tot_vs_bar_%s_%s.pdf",plotDir.c_str(),mapIt2.first.c_str(),mapIt1.first.c_str()));
                }
	      }

	      if(g_tot_vs_bar[Form("R_%s_%s",mapIt2.first.c_str(),mapIt1.first.c_str())]){
		TGraphErrors* g_totR = g_tot_vs_bar[label2];
		if(g_totR->GetN()>0){
		  g_totR -> SetLineColor(kRed);
		  g_totR -> SetMarkerColor(kRed);
		  g_totR -> SetMarkerStyle(20);
		  g_totR -> Draw("P,same");
		  outFile -> cd();
	          g_totR -> Write(Form("g_tot_vs_barR_%s_%s",mapIt2.first.c_str(),mapIt1.first.c_str()));
	
		  TLegend* legenda = new TLegend(0.55,0.85,0.70,0.89);
		  legenda->AddEntry(g_totR, "Right", "P");
		  legenda -> SetBorderSize(0);
		  legenda->Draw();
		

		  c1 -> Print(Form("%s/summaryPlots/tot/c_tot_vs_bar_%s_%s.png",plotDir.c_str(),mapIt2.first.c_str(),mapIt1.first.c_str()));
	          c1 -> Print(Form("%s/summaryPlots/tot/c_tot_vs_bar_%s_%s.pdf",plotDir.c_str(),mapIt2.first.c_str(),mapIt1.first.c_str()));
		}
	      }

	      

	    }
	  }



  //int d=0;
  std::vector<float> vec_Vov;
  std::vector<float> vec_th;	
  for(auto file: listStep2){
	  //int i=0;
	  TFile* inFile = TFile::Open(file.c_str(),"READ");
	  
	  TList* list = inFile -> GetListOfKeys();
	  TIter next(list);
	  TObject* object = 0;
	  while( (object = next()) )
	  {
	    std::string name(object->GetName());
	    std::vector<std::string> tokens = GetTokens(name,'_');
	    std::size_t found;
            std::size_t found1;  
	    found = name.find("data_");
            found1 = name.find("dataRes_");
	    //tree
	    if( found!=std::string::npos )
	    {
	      std::string label(Form("%s_%s_%s",tokens[1].c_str(),tokens[2].c_str(),tokens[3].c_str()));
	      trees[label] = (TTree*)( inFile->Get(name.c_str()) );
	      std::string stepLabel = tokens[2]+"_"+tokens[3];
	      VovLabels[tokens[2]] += 1;
	      thLabels[tokens[3]] += 1;
	      stepLabels.push_back(stepLabel);
	      std::string string_Vov = tokens[2];
	      string_Vov.erase(0,3);
	      map_Vovs[stepLabel] = atof(string_Vov.c_str());
	      std::string string_th = tokens[3];
	      string_th.erase(0,2);
	      map_ths[stepLabel] = atof(string_th.c_str());
	      //i++;
	    }
	  
	  std::sort(stepLabels.begin(),stepLabels.end());
	  stepLabels.erase(std::unique(stepLabels.begin(),stepLabels.end()),stepLabels.end());

	    if( found1!=std::string::npos )
	    {
	      std::string label2(Form("%s_%s_%s_%s",tokens[1].c_str(),tokens[2].c_str(),tokens[3].c_str(),tokens[4].c_str()));
	      trees2[label2] = (TTree*)( inFile->Get(name.c_str()) );
	      std::string stepLabel = tokens[2]+"_"+tokens[3];
	      VovLabels[tokens[2]] += 1;
	      thLabels[tokens[3]] += 1;
	      stepLabels.push_back(stepLabel);
	      std::string string_Vov = tokens[2];
	      string_Vov.erase(0,3);
	      map_Vovs[stepLabel] = atof(string_Vov.c_str());
	      std::string string_th = tokens[3];
	      string_th.erase(0,2);
	      map_ths[stepLabel] = atof(string_th.c_str());
              float Vov = map_Vovs[stepLabel];
	      float vth1 = map_ths[stepLabel];
              vec_Vov.push_back(Vov);
	      vec_th.push_back(vth1);
	      EnBinLabels[tokens[4]] += 1;
	      std::string string_EnBin = tokens[4];
	      string_EnBin.erase(0,9);
	      map_EnBin[stepLabel] = atof(string_EnBin.c_str());
	      //d++;
	    }
	  
	  }
	  //std::cout<<"i="<<i<<std::endl;
	  //std::cout<<"d="<<d<<std::endl;

  }

  std::sort(vec_Vov.begin(),vec_Vov.end());
  std::sort(vec_th.begin(),vec_th.end());
  vec_Vov.erase(std::unique(vec_Vov.begin(),vec_Vov.end()),vec_Vov.end());
  vec_th.erase(std::unique(vec_th.begin(),vec_th.end()),vec_th.end()); 


  std::map<float,float> map_en511;
  std::map<float,float> map_en1275;
  std::map<float,float> map_enRatio;
 
  std::map<float,float> map_en511R;
  std::map<float,float> map_en1275R;
  std::map<float,float> map_enRatioR;

  std::map<float,float> map_en511L;
  std::map<float,float> map_en1275L;
  std::map<float,float> map_enRatioL;

  std::map<float,TGraphErrors*> g_en511_vs_th;
  std::map<float,TGraphErrors*> g_en511_vs_Vov;
  std::map<float,TGraphErrors*> g_en511_vs_bar;
  std::map<float,TGraphErrors*> g_en1275_vs_th;
  std::map<float,TGraphErrors*> g_en1275_vs_Vov;
  std::map<float,TGraphErrors*> g_en1275_vs_bar;
  std::map<float,TGraphErrors*> g_enRatio_vs_th;
  std::map<float,TGraphErrors*> g_enRatio_vs_Vov;
  std::map<float,TGraphErrors*> g_enRatio_vs_bar;
  std::map<float,TGraphErrors*> g_en511_vs_th_L;
  std::map<float,TGraphErrors*> g_en511_vs_Vov_L;
  std::map<float,TGraphErrors*> g_en511_vs_bar_L;
  std::map<float,TGraphErrors*> g_en1275_vs_th_L;
  std::map<float,TGraphErrors*> g_en1275_vs_Vov_L;
  std::map<float,TGraphErrors*> g_en1275_vs_bar_L;
  std::map<float,TGraphErrors*> g_enRatio_vs_th_L;
  std::map<float,TGraphErrors*> g_enRatio_vs_Vov_L;
  std::map<float,TGraphErrors*> g_enRatio_vs_bar_L;
  std::map<float,TGraphErrors*> g_en511_vs_th_R;
  std::map<float,TGraphErrors*> g_en511_vs_Vov_R;
  std::map<float,TGraphErrors*> g_en511_vs_bar_R;
  std::map<float,TGraphErrors*> g_en1275_vs_th_R;
  std::map<float,TGraphErrors*> g_en1275_vs_Vov_R;
  std::map<float,TGraphErrors*> g_en1275_vs_bar_R;
  std::map<float,TGraphErrors*> g_enRatio_vs_th_R;
  std::map<float,TGraphErrors*> g_enRatio_vs_Vov_R;
  std::map<float,TGraphErrors*> g_enRatio_vs_bar_R;
  
  //int index( (10000*int(Vov*100.)) + (100*vth1) + iBar );


  for (auto mapIt : trees){
    float energy511 = -1;
    float energy1275 = -1;
    float energy511L = -1;
    float energy1275L = -1;
    float energy511R = -1;
    float energy1275R = -1;
    float theIndex = -1;
    mapIt.second -> SetBranchStatus("*",0);
    mapIt.second -> SetBranchStatus("energyPeak511LR",  1); mapIt.second -> SetBranchAddress("energyPeak511LR",  &energy511);
    mapIt.second -> SetBranchStatus("energyPeak1275LR",  1); mapIt.second -> SetBranchAddress("energyPeak1275LR",  &energy1275);
    mapIt.second -> SetBranchStatus("energyPeak511L",  1); mapIt.second -> SetBranchAddress("energyPeak511L",  &energy511L);
    mapIt.second -> SetBranchStatus("energyPeak1275L",  1); mapIt.second -> SetBranchAddress("energyPeak1275L",  &energy1275L);
    mapIt.second -> SetBranchStatus("energyPeak511R",  1); mapIt.second -> SetBranchAddress("energyPeak511R",  &energy511R);
    mapIt.second -> SetBranchStatus("energyPeak1275R",  1); mapIt.second -> SetBranchAddress("energyPeak1275R",  &energy1275R);
    mapIt.second -> SetBranchStatus("indexID",  1); mapIt.second -> SetBranchAddress("indexID",   &theIndex);
    int nEntries = mapIt.second->GetEntries();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      mapIt.second -> GetEntry(entry);
    }
    if(energy511 > 0.1){  
      map_en511[theIndex] = energy511;
    }
    if(energy1275 > 0.1){  
      map_en1275[theIndex] = energy1275;
    }
    if(energy511 > 0.1 && energy1275 > 0.1){
      map_enRatio[theIndex] = energy1275/energy511;
    }

    if(energy511L > 0.1){  
      map_en511L[theIndex] = energy511L;
    }
    if(energy1275L > 0.1){  
      map_en1275L[theIndex] = energy1275L;
    }
    if(energy511L > 0.1 && energy1275L > 0.1){
      map_enRatioL[theIndex] = energy1275L/energy511L;
    }

    if(energy511R > 0.1){  
      map_en511R[theIndex] = energy511R;
    }
    if(energy1275R > 0.1){  
      map_en1275R[theIndex] = energy1275R;
    }
    if(energy511R > 0.1 && energy1275R > 0.1){
      map_enRatioR[theIndex] = energy1275R/energy511R;
    }


  }

  //energy Peak 511 vs th and vs Vov and vs iBar (LR, L, R)
  for(std::map<float,float>::iterator index = map_en511.begin(); index != map_en511.end(); index++){
    //std::cout<<"index"<<index->first<<std::endl;
    float Vov;
    int th;
    int iBar;
    Vov = float((int(index->first/10000))/100);
    th = int((index->first - Vov*1000000)/100);
    iBar = int(index->first - Vov*1000000 - th*100);
    int Vov_iBar_ID;
    int th_iBar_ID;
    int Vov_th_ID;
    Vov_iBar_ID = 10000*int(Vov*100.) + iBar;
    th_iBar_ID = 100*th + iBar;
    Vov_th_ID =(10000*int(Vov*100.)) + (100*th);

    if( g_en511_vs_th[Vov_iBar_ID] == NULL ){	
      g_en511_vs_th[Vov_iBar_ID] = new TGraphErrors();
    }
    if( g_en511_vs_Vov[th_iBar_ID] == NULL ){	
      g_en511_vs_Vov[th_iBar_ID] = new TGraphErrors();
    }
    if( g_en511_vs_bar[Vov_th_ID] == NULL ){	
      g_en511_vs_bar[Vov_th_ID] = new TGraphErrors();
    }

    if( g_en511_vs_bar_L[Vov_th_ID] == NULL ){	
      g_en511_vs_bar_L[Vov_th_ID] = new TGraphErrors();
    }
    
    if( g_en511_vs_bar_R[Vov_th_ID] == NULL ){	
      g_en511_vs_bar_R[Vov_th_ID] = new TGraphErrors();
    }
    
    if (index->second>0.1){
      //511 Peak LR
      g_en511_vs_th[Vov_iBar_ID]->SetPoint(g_en511_vs_th[Vov_iBar_ID]->GetN(), th, index->second);
      //std::cout<<"N"<<g_en511_vs_th[Vov_iBar_ID]->GetN()<<"th"<<th<<"en"<<index->second<<"VovIbarID"<<Vov_iBar_ID<<"Vov"<<Vov<<std::endl;
      g_en511_vs_th[Vov_iBar_ID]->SetPointError(g_en511_vs_th[Vov_iBar_ID]->GetN()-1, 0., 0.);

      g_en511_vs_Vov[th_iBar_ID]->SetPoint(g_en511_vs_Vov[th_iBar_ID]->GetN(), Vov, index->second);
      //std::cout<<"N"<<g_en511_vs_Vov[th_iBar_ID]->GetN()<<"Vov"<<Vov<<"en"<<index->second<<"thIbarID"<<th_iBar_ID<<"th"<<th<<std::endl;
      g_en511_vs_Vov[th_iBar_ID]->SetPointError(g_en511_vs_Vov[th_iBar_ID]->GetN()-1, 0., 0.);

      g_en511_vs_bar[Vov_th_ID]->SetPoint(g_en511_vs_bar[Vov_th_ID]->GetN(), iBar, index->second);
      //std::cout<<"N"<<g_en511_vs_bar[Vov_th_ID]->GetN()<<"Vov"<<Vov<<"en"<<index->second<<"thVovID"<<Vov_th_ID<<"th"<<th<<"bar"<<iBar<<std::endl;
      g_en511_vs_bar[Vov_th_ID]->SetPointError(g_en511_vs_bar[Vov_th_ID]->GetN()-1, 0., 0.);
      
    }

    if (map_en511L[index->first]>0.1){
      //511 Peak Left
    /*  g_en511_vs_th_L[Vov_iBar_ID]->SetPoint(g_en511_vs_th_L[Vov_iBar_ID]->GetN(), th, map_en511L[index->first]);
      g_en511_vs_th_L[Vov_iBar_ID]->SetPointError(g_en511_vs_th_L[Vov_iBar_ID]->GetN()-1, 0., 0.);

      g_en511_vs_Vov_L[th_iBar_ID]->SetPoint(g_en511_vs_Vov_L[th_iBar_ID]->GetN(), Vov, map_en511L[index->first]);
      g_en511_vs_Vov_L[th_iBar_ID]->SetPointError(g_en511_vs_Vov_L[th_iBar_ID]->GetN()-1, 0., 0.);*/

      g_en511_vs_bar_L[Vov_th_ID]->SetPoint(g_en511_vs_bar_L[Vov_th_ID]->GetN(), iBar, map_en511L[index->first]);
      g_en511_vs_bar_L[Vov_th_ID]->SetPointError(g_en511_vs_bar_L[Vov_th_ID]->GetN()-1, 0., 0.);
      
    }

    if (map_en511R[index->first]>0.1){
      //511 Peak Right
    /*  g_en511_vs_th_R[Vov_iBar_ID]->SetPoint(g_en511_vs_th_R[Vov_iBar_ID]->GetN(), th, map_en511R[index->first]);
      g_en511_vs_th_R[Vov_iBar_ID]->SetPointError(g_en511_vs_th_R[Vov_iBar_ID]->GetN()-1, 0., 0.);

      g_en511_vs_Vov_R[th_iBar_ID]->SetPoint(g_en511_vs_Vov_R[th_iBar_ID]->GetN(), Vov, map_en511R[index->first]);
      g_en511_vs_Vov_R[th_iBar_ID]->SetPointError(g_en511_vs_Vov_R[th_iBar_ID]->GetN()-1, 0., 0.);*/

      g_en511_vs_bar_R[Vov_th_ID]->SetPoint(g_en511_vs_bar_R[Vov_th_ID]->GetN(), iBar, map_en511R[index->first]);
      g_en511_vs_bar_R[Vov_th_ID]->SetPointError(g_en511_vs_bar_R[Vov_th_ID]->GetN()-1, 0., 0.);
      
    }

  }


  //energy Peak 1275 vs th and vs Vov and vs iBar (LR, L, R)
  for(std::map<float,float>::iterator index = map_en1275.begin(); index != map_en1275.end(); index++){
    //std::cout<<"index"<<index->first<<std::endl;
    float Vov;
    int th;
    int iBar;
    Vov = float((int(index->first/10000))/100);
    th = int((index->first - Vov*1000000)/100);
    iBar = int(index->first - Vov*1000000 - th*100);
    int Vov_iBar_ID;
    int th_iBar_ID;
    int Vov_th_ID;
    Vov_iBar_ID = 10000*int(Vov*100.) + iBar;
    th_iBar_ID = 100*th + iBar;
    Vov_th_ID = 10000*int(Vov*100.) + 100*th;

    if( g_en1275_vs_th[Vov_iBar_ID] == NULL ){	
      g_en1275_vs_th[Vov_iBar_ID] = new TGraphErrors();
    }
    if( g_en1275_vs_Vov[th_iBar_ID] == NULL ){	
      g_en1275_vs_Vov[th_iBar_ID] = new TGraphErrors();
    }
    if( g_en1275_vs_bar[Vov_th_ID] == NULL ){	
      g_en1275_vs_bar[Vov_th_ID] = new TGraphErrors();
    }
    
    if( g_en1275_vs_bar_L[Vov_th_ID] == NULL ){	
      g_en1275_vs_bar_L[Vov_th_ID] = new TGraphErrors();
    }
    
    if( g_en1275_vs_bar_R[Vov_th_ID] == NULL ){	
      g_en1275_vs_bar_R[Vov_th_ID] = new TGraphErrors();
    }
    
    if (index->second>0.1){
      //1275 Peak LR
      g_en1275_vs_th[Vov_iBar_ID]->SetPoint(g_en1275_vs_th[Vov_iBar_ID]->GetN(), th, index->second);
      g_en1275_vs_th[Vov_iBar_ID]->SetPointError(g_en1275_vs_th[Vov_iBar_ID]->GetN()-1, 0., 0.);

      g_en1275_vs_Vov[th_iBar_ID]->SetPoint(g_en1275_vs_Vov[th_iBar_ID]->GetN(), Vov, index->second);
      g_en1275_vs_Vov[th_iBar_ID]->SetPointError(g_en1275_vs_Vov[th_iBar_ID]->GetN()-1, 0., 0.);

      g_en1275_vs_bar[Vov_th_ID]->SetPoint(g_en1275_vs_bar[Vov_th_ID]->GetN(), iBar, index->second);
      g_en1275_vs_bar[Vov_th_ID]->SetPointError(g_en1275_vs_bar[Vov_th_ID]->GetN()-1, 0., 0.);
    }
	
    if (map_en1275R[index->first]>0.1){
      //1275 Peak Right
    /*  g_en1275_vs_th_R[Vov_iBar_ID]->SetPoint(g_en1275_vs_th_R[Vov_iBar_ID]->GetN(), th, map_en1275R[index->first]);
      g_en1275_vs_th_R[Vov_iBar_ID]->SetPointError(g_en1275_vs_th_R[Vov_iBar_ID]->GetN()-1, 0., 0.);

      g_en1275_vs_Vov_R[th_iBar_ID]->SetPoint(g_en1275_vs_Vov_R[th_iBar_ID]->GetN(), Vov, map_en1275R[index->first]);
      g_en1275_vs_Vov_R[th_iBar_ID]->SetPointError(g_en1275_vs_Vov_R[th_iBar_ID]->GetN()-1, 0., 0.);*/

      g_en1275_vs_bar_R[Vov_th_ID]->SetPoint(g_en1275_vs_bar_R[Vov_th_ID]->GetN(), iBar, map_en1275R[index->first]);
      g_en1275_vs_bar_R[Vov_th_ID]->SetPointError(g_en1275_vs_bar_R[Vov_th_ID]->GetN()-1, 0., 0.);
    }
		
    if (map_en1275L[index->first]>0.1){
      //1275 Peak Left
     /* g_en1275_vs_th_L[Vov_iBar_ID]->SetPoint(g_en1275_vs_th_L[Vov_iBar_ID]->GetN(), th, map_en1275L[index->first]);
      g_en1275_vs_th_L[Vov_iBar_ID]->SetPointError(g_en1275_vs_th_L[Vov_iBar_ID]->GetN()-1, 0., 0.);

      g_en1275_vs_Vov_L[th_iBar_ID]->SetPoint(g_en1275_vs_Vov_L[th_iBar_ID]->GetN(), Vov, map_en1275L[index->first]);
      g_en1275_vs_Vov_L[th_iBar_ID]->SetPointError(g_en1275_vs_Vov_L[th_iBar_ID]->GetN()-1, 0., 0.);*/

      g_en1275_vs_bar_L[Vov_th_ID]->SetPoint(g_en1275_vs_bar_L[Vov_th_ID]->GetN(), iBar, map_en1275L[index->first]);
      g_en1275_vs_bar_L[Vov_th_ID]->SetPointError(g_en1275_vs_bar_L[Vov_th_ID]->GetN()-1, 0., 0.);
    }
  }


  //energy Peak Ratio vs th and vs Vov vs iBar (LR, L, R)
  for(std::map<float,float>::iterator index = map_enRatio.begin(); index != map_enRatio.end(); index++){
    //std::cout<<"index"<<index->first<<std::endl;
    float Vov;
    int th;
    int iBar;
    Vov = float((int(index->first/10000))/100);
    th = int((index->first - Vov*1000000)/100);
    iBar = int(index->first - Vov*1000000 - th*100);
    int Vov_iBar_ID;
    int th_iBar_ID;
    int Vov_th_ID;
    Vov_iBar_ID = 10000*int(Vov*100.) + iBar;
    th_iBar_ID = 100*th + iBar;
    Vov_th_ID = 10000*int(Vov*100.) + 100*th; 

    if( g_enRatio_vs_th[Vov_iBar_ID] == NULL ){	
      g_enRatio_vs_th[Vov_iBar_ID] = new TGraphErrors();
    }
    if( g_enRatio_vs_Vov[th_iBar_ID] == NULL ){	
      g_enRatio_vs_Vov[th_iBar_ID] = new TGraphErrors();
    }

    if( g_enRatio_vs_bar[Vov_th_ID] == NULL ){	
      g_enRatio_vs_bar[Vov_th_ID] = new TGraphErrors();
    }
		
    if( g_enRatio_vs_bar_L[Vov_th_ID] == NULL ){	
      g_enRatio_vs_bar_L[Vov_th_ID] = new TGraphErrors();
    }
    
    if( g_enRatio_vs_bar_R[Vov_th_ID] == NULL ){	
      g_enRatio_vs_bar_R[Vov_th_ID] = new TGraphErrors();
    }
    
    if (index->second>0.1){
      //Ratio Peak LR
      g_enRatio_vs_th[Vov_iBar_ID]->SetPoint(g_enRatio_vs_th[Vov_iBar_ID]->GetN(), th, index->second);
      g_enRatio_vs_th[Vov_iBar_ID]->SetPointError(g_enRatio_vs_th[Vov_iBar_ID]->GetN()-1, 0., 0.);

      g_enRatio_vs_Vov[th_iBar_ID]->SetPoint(g_enRatio_vs_Vov[th_iBar_ID]->GetN(), Vov, index->second);
      g_enRatio_vs_Vov[th_iBar_ID]->SetPointError(g_enRatio_vs_Vov[th_iBar_ID]->GetN()-1, 0., 0.);

      g_enRatio_vs_bar[Vov_th_ID]->SetPoint(g_enRatio_vs_bar[Vov_th_ID]->GetN(), iBar, index->second);
      g_enRatio_vs_bar[Vov_th_ID]->SetPointError(g_enRatio_vs_bar[Vov_th_ID]->GetN()-1, 0., 0.);
    }

    if (map_enRatioR[index->first]>0.1){
      //Ratio Peak R
     /* g_enRatio_vs_th_R[Vov_iBar_ID]->SetPoint(g_enRatio_vs_th_R[Vov_iBar_ID]->GetN(), th, map_enRatioR[index->first]);
      g_enRatio_vs_th_R[Vov_iBar_ID]->SetPointError(g_enRatio_vs_th_R[Vov_iBar_ID]->GetN()-1, 0., 0.);

      g_enRatio_vs_Vov_R[th_iBar_ID]->SetPoint(g_enRatio_vs_Vov_R[th_iBar_ID]->GetN(), Vov, map_enRatioR[index->first]);
      g_enRatio_vs_Vov_R[th_iBar_ID]->SetPointError(g_enRatio_vs_Vov_R[th_iBar_ID]->GetN()-1, 0., 0.);*/

      g_enRatio_vs_bar_R[Vov_th_ID]->SetPoint(g_enRatio_vs_bar_R[Vov_th_ID]->GetN(), iBar, map_enRatioR[index->first]);
      g_enRatio_vs_bar_R[Vov_th_ID]->SetPointError(g_enRatio_vs_bar_R[Vov_th_ID]->GetN()-1, 0., 0.);
    }

    if (map_enRatioL[index->first]>0.1){
      //Ratio Peak L
     /* g_enRatio_vs_th_L[Vov_iBar_ID]->SetPoint(g_enRatio_vs_th_L[Vov_iBar_ID]->GetN(), th, map_enRatioL[index->first]);
      g_enRatio_vs_th_L[Vov_iBar_ID]->SetPointError(g_enRatio_vs_th_L[Vov_iBar_ID]->GetN()-1, 0., 0.);

      g_enRatio_vs_Vov_L[th_iBar_ID]->SetPoint(g_enRatio_vs_Vov_L[th_iBar_ID]->GetN(), Vov, map_enRatioL[index->first]);
      g_enRatio_vs_Vov_L[th_iBar_ID]->SetPointError(g_enRatio_vs_Vov_L[th_iBar_ID]->GetN()-1, 0., 0.);*/

      g_enRatio_vs_bar_L[Vov_th_ID]->SetPoint(g_enRatio_vs_bar_L[Vov_th_ID]->GetN(), iBar, map_enRatioL[index->first]);
      g_enRatio_vs_bar_L[Vov_th_ID]->SetPointError(g_enRatio_vs_bar_L[Vov_th_ID]->GetN()-1, 0., 0.);
    }
  }
  
  
  
  
  std::map<int, TCanvas*> c_en511_vs_th;
  std::map<int, TCanvas*> c_en511_vs_Vov;
  std::map<int, TCanvas*> c_en511_vs_bar;
  std::map<int, TCanvas*> c_en1275_vs_th;
  std::map<int, TCanvas*> c_en1275_vs_Vov;
  std::map<int, TCanvas*> c_en1275_vs_bar;
  std::map<int, TCanvas*> c_enRatio_vs_th;
  std::map<int, TCanvas*> c_enRatio_vs_Vov;
  std::map<int, TCanvas*> c_enRatio_vs_bar;

  std::map<int, TCanvas*> c_en511_vs_th_L;
  std::map<int, TCanvas*> c_en511_vs_Vov_L;
  std::map<int, TCanvas*> c_en511_vs_bar_L;
  std::map<int, TCanvas*> c_en1275_vs_th_L;
  std::map<int, TCanvas*> c_en1275_vs_Vov_L;
  std::map<int, TCanvas*> c_en1275_vs_bar_L;
  std::map<int, TCanvas*> c_enRatio_vs_th_L;
  std::map<int, TCanvas*> c_enRatio_vs_Vov_L;
  std::map<int, TCanvas*> c_enRatio_vs_bar_L;

  std::map<int, TCanvas*> c_en511_vs_th_R;
  std::map<int, TCanvas*> c_en511_vs_Vov_R;
  std::map<int, TCanvas*> c_en511_vs_bar_R;
  std::map<int, TCanvas*> c_en1275_vs_th_R;
  std::map<int, TCanvas*> c_en1275_vs_Vov_R;
  std::map<int, TCanvas*> c_en1275_vs_bar_R;
  std::map<int, TCanvas*> c_enRatio_vs_th_R;
  std::map<int, TCanvas*> c_enRatio_vs_Vov_R;
  std::map<int, TCanvas*> c_enRatio_vs_bar_R;

  std::map<int, int> iter;
  


  //summary plots Energy Peak ( 511, 1275, Ratio, L-R-LR) vs th
  for(std::map<float,TGraphErrors*>::iterator index = g_en511_vs_th.begin(); index != g_en511_vs_th.end(); index++){
    //std::cout<<"indexfirst= "<<index->first<<std::endl;
    int Vov_iBar_ID;
    float Vov;
    int iBar;
    Vov = float((int(index->first/10000))/100);
    iBar = int(index->first - Vov*1000000);
    Vov_iBar_ID = 10000*int(Vov*100.) + iBar; 
    for(int bar=0; bar<16; bar++){
      if(bar==iBar){
        if(c_en511_vs_th[bar]==NULL){
          c_en511_vs_th[bar] = new TCanvas(Form("c_en511_vs_th_bar%02d",bar),Form("c_en511_vs_th_bar%02d",bar));
          iter[bar] = 0;
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
          hPad -> SetTitle(";threshold [DAC];energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();
	}
				
	/*if(c_en511_vs_th_L[bar]==NULL){
          c_en511_vs_th_L[bar] = new TCanvas(Form("c_en511_vs_th_bar%02dL",bar),Form("c_en511_vs_th_bar%02dL",bar));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
          hPad -> SetTitle(";threshold [DAC];energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();
	}
	
	if(c_en511_vs_th_R[bar]==NULL){
          c_en511_vs_th_R[bar] = new TCanvas(Form("c_en511_vs_th_bar%02dR",bar),Form("c_en511_vs_th_bar%02dR",bar));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,25.) );
          hPad -> SetTitle(";threshold [DAC];energy [a.u.]");
	  hPad -> Draw();
          gPad -> SetGridy();
	}*/

	if(c_en1275_vs_th[bar]==NULL){
          c_en1275_vs_th[bar] = new TCanvas(Form("c_en1275_vs_th_bar%02d",bar),Form("c_en1275_vs_th_bar%02d",bar));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,35.) );
          hPad -> SetTitle(";threshold [DAC];energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();
        }

	/*if(c_en1275_vs_th_L[bar]==NULL){
          c_en1275_vs_th_L[bar] = new TCanvas(Form("c_en1275_vs_th_bar%02dL",bar),Form("c_en1275_vs_th_bar%02dL",bar));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,35.) );
          hPad -> SetTitle(";threshold [DAC];energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();
	}

	if(c_en1275_vs_th_R[bar]==NULL){
          c_en1275_vs_th_R[bar] = new TCanvas(Form("c_en1275_vs_th_bar%02dR",bar),Form("c_en1275_vs_th_bar%02dR",bar));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,35.) );
          hPad -> SetTitle(";threshold [DAC];energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();
	}*/

        if(c_enRatio_vs_th[bar]==NULL){
          c_enRatio_vs_th[bar] = new TCanvas(Form("c_enRatio_vs_th_bar%02d",bar),Form("c_enRatio_vs_th_bar%02d",bar));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,1.,64.,3.) );
          hPad -> SetTitle(";threshold [DAC];energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();
	}
      
        /*if(c_enRatio_vs_th_L[bar]==NULL){
          c_enRatio_vs_th_L[bar] = new TCanvas(Form("c_enRatio_vs_th_bar%02dL",bar),Form("c_enRatio_vs_th_bar%02dL",bar));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,1.,64.,3.) );
          hPad -> SetTitle(";threshold [DAC];energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();
	}

        if(c_enRatio_vs_th_R[bar]==NULL){
          c_enRatio_vs_th_R[bar] = new TCanvas(Form("c_enRatio_vs_th_bar%02dR",bar),Form("c_enRatio_vs_th_bar%02dR",bar));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,1.,64.,3.) );
          hPad -> SetTitle(";threshold [DAC];energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();
	}*/
            

        
        c_en511_vs_th[bar]->cd();
        TGraph* g_energy = g_en511_vs_th[Vov_iBar_ID];
         
        g_energy -> SetLineColor(1+iter[bar]);
        g_energy -> SetMarkerColor(1+iter[bar]);
        g_energy -> SetMarkerStyle(20);
        g_energy -> Draw("PL,same");
	outFile -> cd();
	g_energy -> Write(Form("g_en511_vs_th_Vov%.01f_bar%02d",Vov,bar));

        std::string VovLabel = Form("Vov%.01f", Vov);
        latex = new TLatex(0.55,0.85-0.04*iter[bar],VovLabel.c_str());
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kBlack+iter[bar]);
        latex -> Draw("same");

	/*c_en511_vs_th_L[bar]->cd();
        TGraph* g_energy_R = g_en511_vs_th_L[Vov_iBar_ID];
         
        g_energy_R -> SetLineColor(1+iter[bar]);
        g_energy_R -> SetMarkerColor(1+iter[bar]);
        g_energy_R -> SetMarkerStyle(20);
        g_energy_R -> Draw("PL,same");

       
        latex -> Draw("same");

	c_en511_vs_th_R[bar]->cd();
        TGraph* g_energy_L = g_en511_vs_th_R[Vov_iBar_ID];
         
        g_energy_L -> SetLineColor(1+iter[bar]);
        g_energy_L -> SetMarkerColor(1+iter[bar]);
        g_energy_L -> SetMarkerStyle(20);
        g_energy_L -> Draw("PL,same");

        latex -> Draw("same");*/
				
	c_en1275_vs_th[bar]->cd();
        TGraph* g_energy1275 = g_en1275_vs_th[Vov_iBar_ID];
         
        g_energy1275 -> SetLineColor(1+iter[bar]);
        g_energy1275 -> SetMarkerColor(1+iter[bar]);
        g_energy1275 -> SetMarkerStyle(20);
        g_energy1275 -> Draw("PL,same");
        outFile -> cd();
	g_energy1275 -> Write(Form("g_en1275_vs_th_Vov%.01f_bar%02d",Vov,bar));

        latex -> Draw("same");
	
	/*c_en1275_vs_th_R[bar]->cd();
        TGraph* g_energy1275_R = g_en1275_vs_th_R[Vov_iBar_ID];
         
        g_energy1275_R -> SetLineColor(1+iter[bar]);
        g_energy1275_R -> SetMarkerColor(1+iter[bar]);
        g_energy1275_R -> SetMarkerStyle(20);
        g_energy1275_R -> Draw("PL,same");
        
        latex -> Draw("same");

	c_en1275_vs_th_L[bar]->cd();
        TGraph* g_energy1275_L = g_en1275_vs_th_L[Vov_iBar_ID];
         
        g_energy1275_L -> SetLineColor(1+iter[bar]);
        g_energy1275_L -> SetMarkerColor(1+iter[bar]);
        g_energy1275_L -> SetMarkerStyle(20);
        g_energy1275_L -> Draw("PL,same");
        
        latex -> Draw("same");*/

        c_enRatio_vs_th[bar]->cd();
        TGraph* g_energyRatio = g_enRatio_vs_th[Vov_iBar_ID];
         
        g_energyRatio -> SetLineColor(1+iter[bar]);
        g_energyRatio -> SetMarkerColor(1+iter[bar]);
        g_energyRatio -> SetMarkerStyle(20);
        g_energyRatio -> Draw("PL,same");
	outFile -> cd();
	g_energyRatio -> Write(Form("g_enRatio_vs_th_Vov%.01f_bar%02d",Vov,bar));
	

        latex -> Draw("same");

        /*c_enRatio_vs_th_L[bar]->cd();
        TGraph* g_energyRatio_L = g_enRatio_vs_th_L[Vov_iBar_ID];
         
        g_energyRatio_L -> SetLineColor(1+iter[bar]);
        g_energyRatio_L -> SetMarkerColor(1+iter[bar]);
        g_energyRatio_L -> SetMarkerStyle(20);
        g_energyRatio_L -> Draw("PL,same");

        latex -> Draw("same");

        c_enRatio_vs_th_R[bar]->cd();
        TGraph* g_energyRatio_R = g_enRatio_vs_th_R[Vov_iBar_ID];
         
        g_energyRatio_R -> SetLineColor(1+iter[bar]);
        g_energyRatio_R -> SetMarkerColor(1+iter[bar]);
        g_energyRatio_R -> SetMarkerStyle(20);
        g_energyRatio_R -> Draw("PL,same");

        latex -> Draw("same");*/
      
       
      
      
        ++iter[bar];
      }
    }
    
  }

  for(std::map<int,TCanvas*>::iterator index = c_en511_vs_th.begin(); index != c_en511_vs_th.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energy511_vs_th_bar%02d.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energy511_vs_th_bar%02d.pdf",plotDir.c_str(),index->first));
  }

  /*for(std::map<int,TCanvas*>::iterator index = c_en511_vs_th_R.begin(); index != c_en511_vs_th_R.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energy511_vs_th_bar%02dR.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energy511_vs_th_bar%02dR.pdf",plotDir.c_str(),index->first));
  }
  
  for(std::map<int,TCanvas*>::iterator index = c_en511_vs_th_L.begin(); index != c_en511_vs_th_L.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energy511_vs_th_bar%02dL.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energy511_vs_th_bar%02dL.pdf",plotDir.c_str(),index->first));
  }*/
	
  for(std::map<int,TCanvas*>::iterator index = c_en1275_vs_th.begin(); index != c_en1275_vs_th.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energy1275_vs_th_bar%02d.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energy1275_vs_th_bar%02d.pdf",plotDir.c_str(),index->first));
  }
	
 /* for(std::map<int,TCanvas*>::iterator index = c_en1275_vs_th_L.begin(); index != c_en1275_vs_th_L.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energy1275_vs_th_bar%02dL.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energy1275_vs_th_bar%02dL.pdf",plotDir.c_str(),index->first));
  }
	
  for(std::map<int,TCanvas*>::iterator index = c_en1275_vs_th_R.begin(); index != c_en1275_vs_th_R.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energy1275_vs_th_bar%02dR.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energy1275_vs_th_bar%02dR.pdf",plotDir.c_str(),index->first));
  }*/

  for(std::map<int,TCanvas*>::iterator index = c_enRatio_vs_th.begin(); index != c_enRatio_vs_th.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energyRatio_vs_th_bar%02d.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energyRatio_vs_th_bar%02d.pdf",plotDir.c_str(),index->first));
  }

 /* for(std::map<int,TCanvas*>::iterator index = c_enRatio_vs_th_L.begin(); index != c_enRatio_vs_th_L.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energyRatio_vs_th_bar%02dL.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energyRatio_vs_th_bar%02dL.pdf",plotDir.c_str(),index->first));
  }

  for(std::map<int,TCanvas*>::iterator index = c_enRatio_vs_th_R.begin(); index != c_enRatio_vs_th_R.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energyRatio_vs_th_bar%02dR.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energyRatio_vs_th_bar%02dR.pdf",plotDir.c_str(),index->first));
  }*/
      
  //summary plots Energy Peak 511, 1275, Ratio vs Vov (LR, L, R)
  for(std::map<float,TGraphErrors*>::iterator index = g_en511_vs_Vov.begin(); index != g_en511_vs_Vov.end(); index++){
    //std::cout<<"indexfirst= "<<index->first<<std::endl;
    int th;
    int iBar;
    int th_iBar_ID;
    th = int((index->first)/100);
    iBar = int(index->first - th*100);
    th_iBar_ID = 100*th + iBar;
    
    for(int bar=0; bar<16; bar++){
      if(bar==iBar){
        if(c_en511_vs_Vov[bar]==NULL){
          c_en511_vs_Vov[bar] = new TCanvas(Form("c_en511_vs_Vov_bar%02d",bar),Form("c_en511_vs_Vov_bar%02d",bar));
          iter[bar] = 0;
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,10.,25.) );
          hPad -> SetTitle(";V_{ov} [V];energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();       
	}

	/*if(c_en511_vs_Vov_L[bar]==NULL){
          c_en511_vs_Vov_L[bar] = new TCanvas(Form("c_en511_vs_Vov_bar%02dL",bar),Form("c_en511_vs_Vov_bar%02dL",bar));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,10.,25.) );
          hPad -> SetTitle(";V_{ov} [V];energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();
	}
	
	if(c_en511_vs_Vov_R[bar]==NULL){
          c_en511_vs_Vov_R[bar] = new TCanvas(Form("c_en511_vs_Vov_bar%02dR",bar),Form("c_en511_vs_Vov_bar%02dR",bar));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,10.,25.) );
          hPad -> SetTitle(";V_{ov} [V];energy [a.u.]");
	  hPad -> Draw();
          gPad -> SetGridy();
	}*/

	if(c_en1275_vs_Vov[bar]==NULL){
          c_en1275_vs_Vov[bar] = new TCanvas(Form("c_en1275_vs_Vov_bar%02d",bar),Form("c_en1275_vs_Vov_bar%02d",bar));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,10.,35.) );
          hPad -> SetTitle(";V_{ov} [V];energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();
        }

	/*if(c_en1275_vs_Vov_L[bar]==NULL){
          c_en1275_vs_Vov_L[bar] = new TCanvas(Form("c_en1275_vs_Vov_bar%02dL",bar),Form("c_en1275_vs_Vov_bar%02dL",bar));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,10.,35.) );
          hPad -> SetTitle(";V_{ov} [V];energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();
	}

	if(c_en1275_vs_Vov_R[bar]==NULL){
          c_en1275_vs_Vov_R[bar] = new TCanvas(Form("c_en1275_vs_Vov_bar%02dR",bar),Form("c_en1275_vs_Vov_bar%02dR",bar));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,10.,35.) );
          hPad -> SetTitle(";V_{ov} [V];energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();
	}*/

	if(c_enRatio_vs_Vov[bar]==NULL){
          c_enRatio_vs_Vov[bar] = new TCanvas(Form("c_enRatio_vs_Vov_bar%02d",bar),Form("c_enRatio_vs_Vov_bar%02d",bar));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,1.,10.,3.) );
          hPad -> SetTitle(";V_{ov} [V];energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();       
	}
	
	/*if(c_enRatio_vs_Vov_R[bar]==NULL){
          c_enRatio_vs_Vov_R[bar] = new TCanvas(Form("c_enRatio_vs_Vov_bar%02dR",bar),Form("c_enRatio_vs_Vov_bar%02dR",bar));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,1.,10.,3.) );
          hPad -> SetTitle(";V_{ov} [V];energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();       
	}
	
	if(c_enRatio_vs_Vov_L[bar]==NULL){
          c_enRatio_vs_Vov_L[bar] = new TCanvas(Form("c_enRatio_vs_Vov_bar%02dL",bar),Form("c_enRatio_vs_Vov_bar%02dL",bar));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,1.,10.,3.) );
          hPad -> SetTitle(";V_{ov} [V];energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();       
	}*/


        c_en511_vs_Vov[bar]->cd();
        TGraph* g_energy = g_en511_vs_Vov[th_iBar_ID];
         
        g_energy -> SetLineColor(1+iter[bar]);
        g_energy -> SetMarkerColor(1+iter[bar]);
        g_energy -> SetMarkerStyle(20);
        g_energy -> Draw("PL,same");
	outFile -> cd();
	g_energy -> Write(Form("g_en511_vs_Vov_th%d_bar%02d",th,bar));

        std::string thLabel = Form("th%d", th);
        latex = new TLatex(0.55,0.85-0.04*iter[bar],thLabel.c_str());
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kBlack+iter[bar]);
        latex -> Draw("same");

        /*c_en511_vs_Vov_R[bar]->cd();
        TGraph* g_energy_R = g_en511_vs_Vov_R[th_iBar_ID];
         
        g_energy_R -> SetLineColor(1+iter[bar]);
        g_energy_R -> SetMarkerColor(1+iter[bar]);
        g_energy_R -> SetMarkerStyle(20);
        g_energy_R -> Draw("PL,same");

       
        latex -> Draw("same");

	c_en511_vs_Vov_L[bar]->cd();
        TGraph* g_energy_L = g_en511_vs_Vov_L[th_iBar_ID];
         
        g_energy_L -> SetLineColor(1+iter[bar]);
        g_energy_L -> SetMarkerColor(1+iter[bar]);
        g_energy_L -> SetMarkerStyle(20);
        g_energy_L -> Draw("PL,same");

        latex -> Draw("same");*/
				
	c_en1275_vs_Vov[bar]->cd();
        TGraph* g_energy1275 = g_en1275_vs_Vov[th_iBar_ID];
         
        g_energy1275 -> SetLineColor(1+iter[bar]);
        g_energy1275 -> SetMarkerColor(1+iter[bar]);
        g_energy1275 -> SetMarkerStyle(20);
        g_energy1275 -> Draw("PL,same");
	outFile -> cd();
	g_energy1275 -> Write(Form("g_en1275_vs_Vov_th%d_bar%02d",th,bar));

        latex -> Draw("same");
	
	/*c_en1275_vs_Vov_R[bar]->cd();
        TGraph* g_energy1275_R = g_en1275_vs_Vov_R[th_iBar_ID];
         
        g_energy1275_R -> SetLineColor(1+iter[bar]);
        g_energy1275_R -> SetMarkerColor(1+iter[bar]);
        g_energy1275_R -> SetMarkerStyle(20);
        g_energy1275_R -> Draw("PL,same");
        
        latex -> Draw("same");

	c_en1275_vs_Vov_L[bar]->cd();
        TGraph* g_energy1275_L = g_en1275_vs_Vov_L[th_iBar_ID];
         
        g_energy1275_L -> SetLineColor(1+iter[bar]);
        g_energy1275_L -> SetMarkerColor(1+iter[bar]);
        g_energy1275_L -> SetMarkerStyle(20);
        g_energy1275_L -> Draw("PL,same");
        
        latex -> Draw("same");*/
      
	c_enRatio_vs_Vov[bar]->cd();
        TGraph* g_energyRatio = g_enRatio_vs_Vov[th_iBar_ID];
         
        g_energyRatio -> SetLineColor(1+iter[bar]);
        g_energyRatio -> SetMarkerColor(1+iter[bar]);
        g_energyRatio -> SetMarkerStyle(20);
        g_energyRatio -> Draw("PL,same");
	outFile -> cd();
	g_energyRatio -> Write(Form("g_enRatio_vs_Vov_th%d_bar%02d",th,bar));

       
        latex -> Draw("same");

	/*c_enRatio_vs_Vov_L[bar]->cd();
        TGraph* g_energyRatio_L = g_enRatio_vs_Vov_L[th_iBar_ID];
         
        g_energyRatio_L -> SetLineColor(1+iter[bar]);
        g_energyRatio_L -> SetMarkerColor(1+iter[bar]);
        g_energyRatio_L -> SetMarkerStyle(20);
        g_energyRatio_L -> Draw("PL,same");

       
        latex -> Draw("same");

	c_enRatio_vs_Vov_R[bar]->cd();
        TGraph* g_energyRatio_R = g_enRatio_vs_Vov_R[th_iBar_ID];
         
        g_energyRatio_R -> SetLineColor(1+iter[bar]);
        g_energyRatio_R -> SetMarkerColor(1+iter[bar]);
        g_energyRatio_R -> SetMarkerStyle(20);
        g_energyRatio_R -> Draw("PL,same");

       
        latex -> Draw("same");*/

        ++iter[bar];
      }
    }
    
  }
 
  for(std::map<int,TCanvas*>::iterator index = c_en511_vs_Vov.begin(); index != c_en511_vs_Vov.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energy511_vs_Vov_bar%02d.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energy511_vs_Vov_bar%02d.pdf",plotDir.c_str(),index->first));
  }

 /* for(std::map<int,TCanvas*>::iterator index = c_en511_vs_Vov_R.begin(); index != c_en511_vs_Vov_R.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energy511_vs_Vov_bar%02dR.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energy511_vs_Vov_bar%02dR.pdf",plotDir.c_str(),index->first));
  }
  
  for(std::map<int,TCanvas*>::iterator index = c_en511_vs_Vov_L.begin(); index != c_en511_vs_Vov_L.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energy511_vs_Vov_bar%02dL.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energy511_vs_Vov_bar%02dL.pdf",plotDir.c_str(),index->first));
  }*/
	
  for(std::map<int,TCanvas*>::iterator index = c_en1275_vs_Vov.begin(); index != c_en1275_vs_Vov.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energy1275_vs_Vov_bar%02d.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energy1275_vs_Vov_bar%02d.pdf",plotDir.c_str(),index->first));
  }
	
  /*for(std::map<int,TCanvas*>::iterator index = c_en1275_vs_Vov_L.begin(); index != c_en1275_vs_Vov_L.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energy1275_vs_Vov_bar%02dL.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energy1275_vs_Vov_bar%02dL.pdf",plotDir.c_str(),index->first));
  }
	
  for(std::map<int,TCanvas*>::iterator index = c_en1275_vs_Vov_R.begin(); index != c_en1275_vs_Vov_R.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energy1275_vs_Vov_bar%02dR.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energy1275_vs_Vov_bar%02dR.pdf",plotDir.c_str(),index->first));
  } */  

  for(std::map<int,TCanvas*>::iterator index = c_enRatio_vs_Vov.begin(); index != c_enRatio_vs_Vov.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energyRatio_vs_Vov_bar%02d.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energyRatio_vs_Vov_bar%02d.pdf",plotDir.c_str(),index->first));
  }
  /*for(std::map<int,TCanvas*>::iterator index = c_enRatio_vs_Vov_R.begin(); index != c_enRatio_vs_Vov_R.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energyRatio_vs_Vov_bar%02dR.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energyRatio_vs_Vov_bar%02dR.pdf",plotDir.c_str(),index->first));
  }
  for(std::map<int,TCanvas*>::iterator index = c_enRatio_vs_Vov_L.begin(); index != c_enRatio_vs_Vov_L.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energyRatio_vs_Vov_bar%02dL.png",plotDir.c_str(),index->first));
    index->second -> Print(Form("%s/summaryPlots/energy/c_energyRatio_vs_Vov_bar%02dL.pdf",plotDir.c_str(),index->first));
  } */

  
  //summary plots Energy Peak 511, 1275, Ratio vs iBar (LR, L, R)
  for(std::map<float,TGraphErrors*>::iterator index = g_en511_vs_bar.begin(); index != g_en511_vs_bar.end(); index++){
    
    int th;
    float Vov;
    int Vov_th_ID;
    Vov = float((int(index->first/10000))/100);
    th = int((index->first - Vov*1000000)/100);
    Vov_th_ID = 10000*int(Vov*100.) + 100*th;
    
    if( th != refTh ) continue;
    
        if(c_en511_vs_bar[Vov_th_ID]==NULL){
          c_en511_vs_bar[Vov_th_ID] = new TCanvas(Form("c_en511_vs_bar_Vov%.01f_th%d",Vov,th),Form("c_en511_vs_bar_Vov%.01f_th%d",Vov,th));
          iter[Vov_th_ID] = 0;
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,17.,25.) );
          hPad -> SetTitle(";ID bar;energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();       
	}
        
	 if(c_en1275_vs_bar[Vov_th_ID]==NULL){
          c_en1275_vs_bar[Vov_th_ID] = new TCanvas(Form("c_en1275_vs_bar_Vov%.01f_th%d",Vov,th),Form("c_en1275_vs_bar_Vov%.01f_th%d",Vov,th));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,17.,35.) );
          hPad -> SetTitle(";ID bar;energy [a.u.]");
	  hPad -> Draw();
	  gPad -> SetGridy();       
	}
         
	if(c_enRatio_vs_bar[Vov_th_ID]==NULL){
          c_enRatio_vs_bar[Vov_th_ID] = new TCanvas(Form("c_enRatio_vs_bar_Vov%.01f_th%d",Vov,th),Form("c_enRatio_vs_bar_Vov%.01f_th%d",Vov,th));
          TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,1.,17.,3.) );
          hPad -> SetTitle(";ID bar;energy_{1275 keV} / energy_{511 keV}");
	  hPad -> Draw();
	  gPad -> SetGridy();       
	}
        
        c_en511_vs_bar[Vov_th_ID]->cd();
        TGraph* g_energy = g_en511_vs_bar[Vov_th_ID];
        TGraph* g_energy_L = g_en511_vs_bar_L[Vov_th_ID];
        TGraph* g_energy_R = g_en511_vs_bar_R[Vov_th_ID];

        g_energy -> SetLineColor(1+iter[Vov_th_ID]);
        g_energy -> SetMarkerColor(1+iter[Vov_th_ID]);
        g_energy -> SetMarkerStyle(20);
        g_energy -> Draw("PL,same");
	outFile -> cd();
	g_energy -> Write(Form("g_en511_vs_bar_Vov%.01f_th%d",Vov,th));
	g_energy_L -> SetLineStyle(9);
	g_energy_L -> SetMarkerStyle(24);
        g_energy_L -> Draw("PL,same");
	outFile -> cd();
	g_energy_L -> Write(Form("g_en511_vs_barL_Vov%.01f_th%d",Vov,th));
	g_energy_R -> SetLineStyle(9);
	g_energy_R -> SetMarkerStyle(25);
        g_energy_R -> Draw("PL,same");
	outFile -> cd();
	g_energy_R -> Write(Form("g_en511_vs_barR_Vov%.01f_th%d",Vov,th));
        
        /*std::string VovthLabel = Form("Vov%.01f th%d", Vov,th);
        latex = new TLatex(0.55,0.85-0.04*iter[Vov_th_ID],VovthLabel.c_str());
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kBlack+iter[Vov_th_ID]);
        latex -> Draw("same");*/
        
        latex = new TLatex(0.25,0.87, Form("Average = %.1f, RMS = %.1f %%", g_energy ->GetMean(2), g_energy->GetRMS(2)/g_energy ->GetMean(2)*100));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.03);
        latex ->Draw("same");
        
        latex = new TLatex(0.25,0.83, Form("Left = %.1f, RMS = %.1f %%", g_energy_L ->GetMean(2), g_energy_L->GetRMS(2)/g_energy_L ->GetMean(2)*100 ));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.03);
        latex ->Draw("same");

        latex = new TLatex(0.25,0.79, Form("Right = %.1f, RMS = %.1f %%", g_energy_R ->GetMean(2), g_energy_R->GetRMS(2)/g_energy_R ->GetMean(2)*100));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.03);
        latex ->Draw("same");
        
        TLegend* legenda = new TLegend(0.20,0.78,0.25,0.90);
	legenda->AddEntry(g_energy, "", "P");
	legenda->AddEntry(g_energy_L, "", "P");
	legenda->AddEntry(g_energy_R, "", "P");
	legenda -> SetBorderSize(0);
	legenda->Draw();
	
	c_en1275_vs_bar[Vov_th_ID]->cd();
        TGraph* g_energy1275 = g_en1275_vs_bar[Vov_th_ID];
        TGraph* g_energy1275_L = g_en1275_vs_bar_L[Vov_th_ID];
        TGraph* g_energy1275_R = g_en1275_vs_bar_R[Vov_th_ID];
         
        g_energy1275 -> SetLineColor(1+iter[Vov_th_ID]);
        g_energy1275 -> SetMarkerColor(1+iter[Vov_th_ID]);
        g_energy1275 -> SetMarkerStyle(20);
        g_energy1275 -> Draw("PL,same");
	outFile -> cd();
	g_energy1275 -> Write(Form("g_en1275_vs_bar_Vov%.01f_th%d",Vov,th));
	g_energy1275_L -> SetLineStyle(9);
	g_energy1275_L -> SetMarkerStyle(24);
        g_energy1275_L -> Draw("PL,same");
	outFile -> cd();
	g_energy1275_L-> Write(Form("g_en1275_vs_barL_Vov%.01f_th%d",Vov,th));
	g_energy1275_R -> SetLineStyle(9);
	g_energy1275_R -> SetMarkerStyle(25);
        g_energy1275_R -> Draw("PL,same");
	outFile -> cd();
	g_energy1275_R -> Write(Form("g_en1275_vs_barR_Vov%.01f_th%d",Vov,th));
        
        latex = new TLatex(0.25,0.87, Form("Average = %.1f, RMS = %.1f %%", g_energy1275 ->GetMean(2), g_energy1275->GetRMS(2)/g_energy1275 ->GetMean(2)*100));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.03);
        latex ->Draw("same");
        
        latex = new TLatex(0.25,0.83, Form("Left = %.1f, RMS = %.1f %%", g_energy1275_L ->GetMean(2), g_energy1275_L->GetRMS(2)/g_energy1275_L ->GetMean(2)*100 ));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.03);
        latex ->Draw("same");

        latex = new TLatex(0.25,0.79, Form("Right = %.1f, RMS = %.1f %%", g_energy1275_R ->GetMean(2), g_energy1275_R->GetRMS(2)/g_energy1275_R ->GetMean(2)*100));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.03);
        latex ->Draw("same");
        
	TLegend* legenda2 = new TLegend(0.20,0.78,0.25,0.90);
	legenda2->AddEntry(g_energy1275, "", "P");
	legenda2->AddEntry(g_energy1275_L, "", "P");
	legenda2->AddEntry(g_energy1275_R, "", "P");
	legenda2 -> SetBorderSize(0);
	legenda2->Draw();
        
	c_enRatio_vs_bar[Vov_th_ID]->cd();
        TGraph* g_energyRatio = g_enRatio_vs_bar[Vov_th_ID];
         
        g_energyRatio -> SetLineColor(1+iter[Vov_th_ID]);
        g_energyRatio -> SetMarkerColor(1+iter[Vov_th_ID]);
        g_energyRatio -> SetMarkerStyle(20);
        g_energyRatio -> Draw("PL,same");
	outFile -> cd();
	g_energyRatio -> Write(Form("g_enRatio_vs_bar_Vov%.01f_th%d",Vov,th));
        
        
        ++iter[Vov_th_ID];
  }
  
  for(std::map<int,TCanvas*>::iterator index = c_en511_vs_bar.begin(); index != c_en511_vs_bar.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energy511_vs_bar_Vov%.01f_th%d.png",plotDir.c_str(),float((int(index->first/10000))/100),int((index->first - (float((int(index->first/10000))/100))*1000000)/100)));
    index->second-> Print(Form("%s/summaryPlots/energy/c_energy511_vs_bar_Vov%.01f_th%d.pdf",plotDir.c_str(),float((int(index->first/10000))/100),int((index->first - (float((int(index->first/10000))/100))*1000000)/100)));
  }
  
  for(std::map<int,TCanvas*>::iterator index = c_en1275_vs_bar.begin(); index != c_en1275_vs_bar.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energy1275_vs_bar_Vov%.01f_th%d.png",plotDir.c_str(),float((int(index->first/10000))/100),int((index->first - (float((int(index->first/10000))/100))*1000000)/100)));
    index->second-> Print(Form("%s/summaryPlots/energy/c_energy1275_vs_bar_Vov%.01f_th%d.pdf",plotDir.c_str(),float((int(index->first/10000))/100),int((index->first - (float((int(index->first/10000))/100))*1000000)/100)));
  }
  
  for(std::map<int,TCanvas*>::iterator index = c_enRatio_vs_bar.begin(); index != c_enRatio_vs_bar.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/energy/c_energyRatio_vs_bar_Vov%.01f_th%d.png",plotDir.c_str(),float((int(index->first/10000))/100),int((index->first - (float((int(index->first/10000))/100))*1000000)/100)));
    index->second-> Print(Form("%s/summaryPlots/energy/c_energyRatio_vs_bar_Vov%.01f_th%d.pdf",plotDir.c_str(),float((int(index->first/10000))/100),int((index->first - (float((int(index->first/10000))/100))*1000000)/100)));
  }
  
  
  
  
  // Time Resolution
  std::map<double,float> map_timeRes;
  std::map<double,float> map_timeRes511;
  std::map<double,float> map_timeRes1275;
  std::map<double,float> map_enBin;
  std::map<double,float> map_enBin511;
  std::map<double,float> map_enBin1275;

  std::map<double,TGraphErrors*> g_timeRes_vs_th;
  std::map<double,TGraphErrors*> g_timeRes_vs_Vov;
  std::map<double,TGraphErrors*> g_timeRes511_vs_bar;
  std::map<double,TGraphErrors*> g_timeRes1275_vs_bar;
  std::map<double,TGraphErrors*> g_timeRes_vs_enBin;

  for (auto mapIt : trees2){
    float timeRes = -1;
    double theIndex2 = -1;
    float enBin = -1;
    mapIt.second -> SetBranchStatus("*",0);
    mapIt.second -> SetBranchStatus("timeResolution",  1); mapIt.second -> SetBranchAddress("timeResolution",  &timeRes);
    mapIt.second -> SetBranchStatus("indexID2",  1); mapIt.second -> SetBranchAddress("indexID2",   &theIndex2);
    mapIt.second -> SetBranchStatus("energyBin",  1); mapIt.second -> SetBranchAddress("energyBin",   &enBin);
    int nEntries = mapIt.second->GetEntries();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      mapIt.second -> GetEntry(entry);
    }
    //std::cout<<"Time Res = "<<timeRes<<" indexID"<<theIndex2<<"enBin"<<enBin<<std::endl;
    if(timeRes > 0){  
      map_timeRes[theIndex2] = timeRes;
      map_enBin[theIndex2] = enBin;
      if (int(theIndex2/10000000) == 3){
        map_enBin511[theIndex2] = enBin;
	map_timeRes511[theIndex2] = timeRes;
      }
      if (int(theIndex2/10000000) == 8){
        map_enBin1275[theIndex2] = enBin;
	map_timeRes1275[theIndex2] = timeRes;
      }
      
    }
  }


  //TimeRes vs th, vs Vov, vs enBin and vs iBar
 for(std::map<double,float>::iterator index = map_timeRes.begin(); index != map_timeRes.end(); index++){
    //std::cout<<"index"<<index->first<<std::endl;
    // int index2( (1000000*energyBin+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );  
    float Vov;
    int th;
    int iBar;
    int enBin;
    enBin = int(index->first/10000000.);
    Vov = float(int(double(index->first-enBin*10000000)/10000.)/100);
    th = int(double(index->first - enBin*10000000 - Vov*1000000)/100);
    iBar = int(double(index->first - enBin*10000000 - Vov*1000000 - th*100));
    double Vov_iBar_enBin_ID;
    Vov_iBar_enBin_ID = double(10000000*enBin) + double(Vov*1000000) + double(iBar);
    double th_iBar_enBin_ID;
    th_iBar_enBin_ID = double(10000000*enBin) + double(th*100) + double(iBar);
    double Vov_th_iBar_ID;
    Vov_th_iBar_ID = double(Vov*1000000) + double(th*100) + double(iBar);
    double Vov_th_ID;
    Vov_th_ID = double(Vov*1000000) + double(th*100);
    double index2;
    index2 = double(10000000*enBin) + double(Vov*1000000) + double(th*100) + double(iBar);

    if( g_timeRes_vs_th[Vov_iBar_enBin_ID] == NULL ){
      g_timeRes_vs_th[Vov_iBar_enBin_ID] = new TGraphErrors();
    }
    if( g_timeRes_vs_Vov[th_iBar_enBin_ID] == NULL ){	
      g_timeRes_vs_Vov[th_iBar_enBin_ID] = new TGraphErrors();
    }
    if( g_timeRes_vs_enBin[Vov_th_iBar_ID] == NULL ){	
      g_timeRes_vs_enBin[Vov_th_iBar_ID] = new TGraphErrors();
    }
    
    if (enBin ==3){
      if( g_timeRes511_vs_bar[Vov_th_ID] == NULL ){	
      g_timeRes511_vs_bar[Vov_th_ID] = new TGraphErrors();
      }
    }

    if (enBin ==8){
      if( g_timeRes1275_vs_bar[Vov_th_ID] == NULL ){	
      g_timeRes1275_vs_bar[Vov_th_ID] = new TGraphErrors();
      }
    }


    if (index->second>0){
      g_timeRes_vs_th[Vov_iBar_enBin_ID]->SetPoint(g_timeRes_vs_th[Vov_iBar_enBin_ID]->GetN(), th, (index->second)/2.);
      //std::cout<<"N"<<g_timeRes_vs_th[Vov_iBar_enBin_ID]->GetN()<<"th"<<th<<"RES "<<(index->second)/2<<"   VovIbarenBinID"<<Vov_iBar_enBin_ID<<"Vov"<<Vov<<std::endl;
      g_timeRes_vs_th[Vov_iBar_enBin_ID]->SetPointError(g_timeRes_vs_th[Vov_iBar_enBin_ID]->GetN()-1, 0., 0.);

      g_timeRes_vs_Vov[th_iBar_enBin_ID]->SetPoint(g_timeRes_vs_Vov[th_iBar_enBin_ID]->GetN(), Vov, (index->second)/2.);
      //std::cout<<"N"<<g_timeRes_vs_Vov[th_iBar_enBin_ID]->GetN()<<"th"<<th<<"RES "<<(index->second)/2<<"   thIbarenBinID"<<th_iBar_enBin_ID<<"Vov"<<Vov<<std::endl;
      g_timeRes_vs_Vov[th_iBar_enBin_ID]->SetPointError(g_timeRes_vs_Vov[th_iBar_enBin_ID]->GetN()-1, 0., 0.);

      g_timeRes_vs_enBin[Vov_th_iBar_ID]->SetPoint(g_timeRes_vs_enBin[Vov_th_iBar_ID]->GetN(), map_enBin[index->first], (index->second)/2.);
      g_timeRes_vs_enBin[Vov_th_iBar_ID]->SetPointError(g_timeRes_vs_enBin[Vov_th_iBar_ID]->GetN()-1, 0., 0.);

    }

    if (enBin ==3){
      g_timeRes511_vs_bar[Vov_th_ID]->SetPoint(g_timeRes511_vs_bar[Vov_th_ID]->GetN(), iBar, (index->second)/2);
      //std::cout<<"N"<<g_timeRes511_vs_bar[Vov_th_ID]->GetN()<<"iBar"<<iBar<<"RES "<<(index->second)/2<<std::endl;
      g_timeRes511_vs_bar[Vov_th_ID]->SetPointError(g_timeRes511_vs_bar[Vov_th_ID]->GetN()-1, 0., 0.);
    }

    if (enBin ==8){
      g_timeRes1275_vs_bar[Vov_th_ID]->SetPoint(g_timeRes1275_vs_bar[Vov_th_ID]->GetN(), iBar, (index->second)/2);
      //std::cout<<"N"<<g_timeRes1275_vs_bar[Vov_th_ID]->GetN()<<"iBar"<<iBar<<"RES "<<(index->second)/2<<std::endl;
      g_timeRes1275_vs_bar[Vov_th_ID]->SetPointError(g_timeRes1275_vs_bar[Vov_th_ID]->GetN()-1, 0., 0.);
    }

    
      
  }
  
  
  //float tResMin = opts.GetOpt<float>("Plots.tResMin");
  //float tResMax = opts.GetOpt<float>("Plots.tResMax");
  std::map<double, TCanvas*> c_timeRes_vs_th;
  std::map<double, TCanvas*> c_timeRes_vs_Vov;
  std::map<double, TCanvas*> c_timeRes_vs_enBin;
  std::map<double, TCanvas*> c_timeRes_vs_bar;
  std::map<double, TCanvas*> c_timeRes_vs_enBin_thBest;
  std::map<double, float> tBest;
  double Vov_iBar_enBin_ID;
  double th_iBar_enBin_ID;
  double iBar_enBin_ID;
  double Vov_th_ID;
  float Vov;
  int iBar;
  int enBin;
  int th;

  //summary plot time resolution vs th
  for(std::map<double,TGraphErrors*>::iterator index = g_timeRes_vs_th.begin(); index != g_timeRes_vs_th.end(); index++){
 
    enBin = int(double(index->first/10000000.));
    Vov = float(int(double(index->first-enBin*10000000.)/10000.)/100.);
    iBar = int(double(index->first - enBin*10000000 - Vov*1000000));
    Vov_iBar_enBin_ID = double(10000000*enBin) + double(Vov*1000000) + double(iBar);
    iBar_enBin_ID = 10000000*enBin + iBar;

    for(int bar=0; bar<16; bar++){
      if(bar==iBar){
        for(int energyBin=1; energyBin<9; energyBin++){
          if(energyBin==enBin){
            if(enBin==3 || enBin==8){  
              if(c_timeRes_vs_th[iBar_enBin_ID]==NULL){
                c_timeRes_vs_th[iBar_enBin_ID] = new TCanvas(Form("c_timeRes_vs_th_bar%02d_enBin%d",bar,energyBin),Form("c_timeRes_vs_th_bar%02d_enBin%d",bar,energyBin));
                iter[Vov_iBar_enBin_ID] = 0;
                TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0,64.,500) );
                hPad -> SetTitle(";threshold [DAC];#sigma_{t_{Diff}}/2 [ps]");
	        hPad -> Draw();
	        gPad -> SetGridy();
              }

	      c_timeRes_vs_th[iBar_enBin_ID]->cd();
              TGraph* g_tRes = g_timeRes_vs_th[Vov_iBar_enBin_ID];
            
              g_tRes -> SetLineColor(1+iter[iBar_enBin_ID]);
              g_tRes -> SetMarkerColor(1+iter[iBar_enBin_ID]);
              g_tRes -> SetMarkerStyle(20);
              g_tRes -> Draw("PL,same");
	      outFile -> cd();
	      g_tRes -> Write(Form("g_timeRes_vs_th_Vov%.01f_bar%02d_enBin%d",Vov,bar,energyBin));
	  
              std::string VovLabel = Form("Vov%.01f", Vov);
              latex = new TLatex(0.55,0.85-0.04*iter[iBar_enBin_ID],VovLabel.c_str());
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kBlack+iter[iBar_enBin_ID]);
              latex -> Draw("same");
      
              ++iter[iBar_enBin_ID];
	      }
	    for (float over : vec_Vov){
	      if (over == Vov){
	        double val = 9999;
	        int xVal = 0;
		for ( int nPoint = 0; nPoint < g_timeRes_vs_th[Vov_iBar_enBin_ID]->GetN(); nPoint++){
		  if ( g_timeRes_vs_th[Vov_iBar_enBin_ID]->GetPointY(nPoint)< val){
		    val = g_timeRes_vs_th[Vov_iBar_enBin_ID]->GetPointY(nPoint);
		    xVal = nPoint;
		  }
		  tBest[Vov_iBar_enBin_ID] = g_timeRes_vs_th[Vov_iBar_enBin_ID]->GetPointX(xVal);		 
                }
	      }	  		  
	    }
	  }
        }
      }
    }
  }  
  

  for(std::map<double,TCanvas*>::iterator index = c_timeRes_vs_th.begin(); index != c_timeRes_vs_th.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/timeResolution/c_timeRes_vs_th_bar%02d_enBin%d.png",plotDir.c_str(),int(index->first-double((int(index->first/10000000))*10000000)), int(double(index->first/10000000.))));
    index->second-> Print(Form("%s/summaryPlots/timeResolution/c_timeRes_vs_th_bar%02d_enBin%d.pdf",plotDir.c_str(),int(index->first-double((int(index->first/10000000))*10000000)), int(double(index->first/10000000.))));
 }

    
  //summary plot time resolution vs Vov
  for(std::map<double,TGraphErrors*>::iterator index = g_timeRes_vs_Vov.begin(); index != g_timeRes_vs_Vov.end(); index++){
    //std::cout<<"indexfirst= "<<index->first<<std::endl;
    enBin = int(double(index->first/10000000.));
    th = int(double(index->first - enBin*10000000)/100);
    iBar = int(double(index->first - enBin*10000000- th*100));
    th_iBar_enBin_ID = 10000000*enBin + th*100 + iBar;
    iBar_enBin_ID = 10000000*enBin + iBar;
    for(int bar=0; bar<16; bar++){
      if(bar==iBar){
        for(int energyBin=1; energyBin<9; energyBin++){
          if(energyBin==enBin){
            if(enBin==3 || enBin==8){
	      for ( int nPoint = 0; nPoint < g_timeRes_vs_Vov[th_iBar_enBin_ID]->GetN(); nPoint++){
	      	float over = g_timeRes_vs_Vov[th_iBar_enBin_ID]->GetPointX(nPoint);
	     	if(th == tBest[double(10000000*enBin) + double(over*1000000) + double(iBar)]){
		  if(c_timeRes_vs_Vov[iBar_enBin_ID]==NULL){
		    c_timeRes_vs_Vov[iBar_enBin_ID] = new TCanvas(Form("c_timeRes_vs_Vov_bar%02d_enBin%d",bar,energyBin),Form("c_timeRes_vs_Vov_bar%02d_enBin%d",bar,energyBin));
		    iter[iBar_enBin_ID] = 0;
		    TH1F* hPad = (TH1F*)( gPad->DrawFrame(0,0,10.,500) );
		    hPad -> SetTitle(";V_{ov} [V];#sigma_{t_{Diff}}/2 [ps]");
	            hPad -> Draw();
		    gPad -> SetGridy();
		  }

		  c_timeRes_vs_Vov[iBar_enBin_ID]->cd();
		  TGraph* g_tRes = g_timeRes_vs_Vov[th_iBar_enBin_ID];

		  g_tRes -> SetLineColor(1+iter[iBar_enBin_ID]);
		  g_tRes -> SetMarkerColor(1+iter[iBar_enBin_ID]);
		  g_tRes -> SetMarkerStyle(20);
		  g_tRes -> Draw("PL,same");
		  outFile -> cd();
	          g_tRes -> Write(Form("g_timeRes_vs_Vov_th%d_bar%02d_enBin%d",th,bar,energyBin));
		  
		  
		  std::string thLabel = Form("th%d", th);
		  latex = new TLatex(0.55,0.85-0.04*iter[iBar_enBin_ID],thLabel.c_str());
	          latex -> SetNDC();
	          latex -> SetTextFont(42);
	          latex -> SetTextSize(0.04);
	          latex -> SetTextColor(kBlack+iter[iBar_enBin_ID]);
	          latex -> Draw("same");
	      
		  ++iter[iBar_enBin_ID];
	        }
	      }
 	    }
	  }
        }
      }
    }
  }  




  for(std::map<double,TCanvas*>::iterator index = c_timeRes_vs_Vov.begin(); index != c_timeRes_vs_Vov.end(); index++){
    index->second-> Print(Form("%s/summaryPlots/timeResolution/c_timeRes_vs_Vov_bar%02d_enBin%d.png",plotDir.c_str(),int(index->first-double((int(index->first/10000000))*10000000)), int(double(index->first/10000000.))));
    index->second-> Print(Form("%s/summaryPlots/timeResolution/c_timeRes_vs_Vov_bar%02d_enBin%d.pdf",plotDir.c_str(),int(index->first-double((int(index->first/10000000))*10000000)), int(double(index->first/10000000.))));
  }

  std::map<double, std::vector<std::pair<float,float>>>vec_tRes;
  std::map<double, std::vector<std::pair<float,float>>>vec_tRes_error;
  //summary plot time resolution vs energy bin
  for(std::map<double,TGraphErrors*>::iterator index = g_timeRes_vs_enBin.begin(); index != g_timeRes_vs_enBin.end(); index++){
  
    Vov = float(int(index->first/10000.)/100.);
    th = int(double(index->first - Vov*1000000)/100.);
    iBar = int(double(index->first - Vov*1000000 - th*100));
    double Vov_th_iBar_ID;
    double iBar_Vov_ID; 
    Vov_th_iBar_ID = Vov*1000000 + th*100 + iBar;
    iBar_Vov_ID = Vov*1000000 + iBar;

    for(int bar=0; bar<16; bar++){
      if(bar==iBar){
        for(int i=0; i<vec_Vov.size(); i++){
          if(vec_Vov[i]==Vov){
	    for(int i=0; i<vec_th.size(); i++){
              if(vec_th[i]==th){
                if(c_timeRes_vs_enBin[Vov_th_iBar_ID]==NULL){
                  c_timeRes_vs_enBin[Vov_th_iBar_ID] = new TCanvas(Form("c_timeRes_vs_enBin_bar%02d_Vov%.01f_th%d",bar,Vov,th),Form("c_timeRes_vs_enBin_bar%02d_Vov%.01f_th%d",bar,Vov,th));
                  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0,0,30., 1.1*g_timeRes_vs_enBin[Vov_th_iBar_ID]->GetHistogram()->GetMaximum()) );
                  hPad -> SetTitle(";energy [a.u.];#sigma_{t_{Diff}}/2 [ps]");
	          hPad -> Draw();
	          gPad -> SetGridy();
		}
		if(c_timeRes_vs_enBin_thBest[iBar_Vov_ID]==NULL){
                  c_timeRes_vs_enBin_thBest[iBar_Vov_ID] = new TCanvas(Form("c_timeRes_vs_enBin_bar%02d_Vov%.01f_bestTh",bar,Vov),Form("c_timeRes_vs_enBin_bar%02d_Vov%.01f_bestTh",bar,Vov));
                  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0,0,30., 1.1*g_timeRes_vs_enBin[Vov_th_iBar_ID]->GetHistogram()->GetMaximum()) );
                  hPad -> SetTitle(";energy [a.u.];#sigma_{t_{Diff}}/2 [ps]");
	          hPad -> Draw();
	          gPad -> SetGridy();
		}
		

	        c_timeRes_vs_enBin[Vov_th_iBar_ID]->cd();
                TGraph* g_tRes_en = g_timeRes_vs_enBin[Vov_th_iBar_ID];
            
                g_tRes_en -> SetLineColor(kBlack);
                g_tRes_en -> SetMarkerColor(kBlack);
                g_tRes_en -> SetMarkerStyle(20);
                g_tRes_en -> Draw("PL");
		outFile -> cd();
	        g_tRes_en -> Write(Form("g_timeRes_vs_enBin_bar%02d_Vov%.01f_th%d",iBar,Vov, th));
		

	      

	        for ( int nPoint = 0; nPoint < g_timeRes_vs_enBin[Vov_th_iBar_ID]->GetN(); nPoint++){
	          float eBin = g_timeRes_vs_enBin[Vov_th_iBar_ID]->GetPointX(nPoint);
		  
                  if(vec_th[i] == tBest[double(10000000*(nPoint+1)) + double(Vov*1000000) + double(bar)]){ 
	            vec_tRes[iBar_Vov_ID].push_back(std::make_pair ( eBin, g_timeRes_vs_enBin[Vov_th_iBar_ID]->GetPointY(nPoint)));
		    vec_tRes_error[iBar_Vov_ID].push_back(std::make_pair ( g_timeRes_vs_enBin[Vov_th_iBar_ID]->GetErrorX(nPoint), g_timeRes_vs_enBin[Vov_th_iBar_ID]->GetErrorY(nPoint)));             
			
	          }
		}

	      }
            }
          }
        }
      }
    }
  } 

 

  for(std::map<double,TCanvas*>::iterator index = c_timeRes_vs_enBin.begin(); index != c_timeRes_vs_enBin.end(); index++){
    Vov = float(int(index->first/10000.)/100.);
    th = int(double(index->first - Vov*1000000)/100.);
    iBar = int(double(index->first - Vov*1000000 - th*100));
    index->second-> Print(Form("%s/summaryPlots/timeResolution/c_timeRes_vs_enBin_bar%02d_Vov%.01f_th%d.png",plotDir.c_str(),iBar,Vov, th));
     index->second-> Print(Form("%s/summaryPlots/timeResolution/c_timeRes_vs_enBin_bar%02d_Vov%.01f_th%d.pdf",plotDir.c_str(),iBar,Vov, th));
  }
  
   for(std::map<double,TCanvas*>::iterator index = c_timeRes_vs_enBin_thBest.begin(); index != c_timeRes_vs_enBin_thBest.end(); index++){
     index->second -> cd();
     Vov = float(int(index->first/10000.)/100.);
     iBar = int(double(index->first - Vov*1000000));
     TGraphErrors * graph = new TGraphErrors();
     int up = vec_tRes[index->first].size();
     for (int i = 0; i < up; i++){
       float min = 9999;
       int minIndex = 0;
       for (int j = 0; j < vec_tRes[index->first].size(); j++){
         if (vec_tRes[index->first][j].first < min){
	   min = vec_tRes[index->first][j].first;
	   minIndex = j;
	 }
       }
       graph -> SetPoint( i, vec_tRes[index->first][minIndex].first, vec_tRes[index->first][minIndex].second);
       graph -> SetPointError( i, vec_tRes_error[index->first][minIndex].first, vec_tRes_error[index->first][minIndex].second);
       vec_tRes[index->first].erase(vec_tRes[index->first].begin() + minIndex);
       vec_tRes_error[index->first].erase(vec_tRes_error[index->first].begin() + minIndex);
    }
    graph-> SetLineColor(kBlack);
    graph-> SetMarkerColor(kBlack);
    graph -> SetMarkerStyle(20);
    graph-> Draw("PL");
    outFile -> cd();
    graph -> Write(Form("g_timeRes_vs_enBin_bar%02d_Vov%.01f_bestTh",iBar,Vov));
    
    index->second-> Print(Form("%s/summaryPlots/timeResolution/c_timeRes_vs_enBin_bar%02d_Vov%.01f_bestTh.png",plotDir.c_str(),iBar,Vov));
    index->second-> Print(Form("%s/summaryPlots/timeResolution/c_timeRes_vs_enBin_bar%02d_Vov%.01f_bestTh.pdf",plotDir.c_str(),iBar,Vov));
  }


  //summary plot time resolution 511 peak and 1275 peak vs iBar
  std::map<double, std::vector<std::pair<float,float>>>vec_tRes_511;
  std::map<double, std::vector<std::pair<float,float>>>vec_tRes_error_511;
  std::map<double, std::vector<std::pair<float,float>>>vec_tRes_error_1275;
  std::map<double, std::vector<std::pair<float,float>>>vec_tRes_1275;
  for(std::map<double,TGraphErrors*>::iterator index = g_timeRes511_vs_bar.begin(); index != g_timeRes511_vs_bar.end(); index++){

    Vov = float(int(index->first/10000)/100.);
    th = int((index->first - Vov*1000000)/100);

    Vov_th_ID = Vov*1000000 + th*100;
    	
    
      if(c_timeRes_vs_bar[Vov]==NULL){
        c_timeRes_vs_bar[Vov] = new TCanvas(Form("c_timeRes_vs_bar_Vov%.01f",Vov),Form("c_timeRes_vs_bar_Vov%.01f",Vov));
        c_timeRes_vs_bar[Vov]->cd();
        TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,17.,500.) );
        hPad -> SetTitle(";ID bar;#sigma_{t_{Diff}}/2 [ps]");
        hPad -> Draw();
        gPad -> SetGridy();    
      }

    
    TGraph* g_timeRes511 = g_timeRes511_vs_bar[Vov_th_ID];
    TGraph* g_timeRes1275 = g_timeRes1275_vs_bar[Vov_th_ID];
        
    for ( int nPoint = 0; nPoint < g_timeRes511->GetN(); nPoint++){
      if (th == tBest[ double(10000000*3) + double(Vov*1000000) + double(int(g_timeRes511->GetPointX(nPoint)))]){	
	vec_tRes_511[Vov].push_back(std::make_pair ( int(g_timeRes511->GetPointX(nPoint)), g_timeRes511->GetPointY(nPoint)));
        vec_tRes_error_511[Vov].push_back(std::make_pair ( g_timeRes511->GetErrorX(nPoint), g_timeRes511->GetErrorY(nPoint)));             
      }
    }
		
    
         
        
    for ( int nPoint = 0; nPoint < g_timeRes1275->GetN(); nPoint++){
      if (th == tBest[ double(10000000*8) + double(Vov*1000000) + double(int(g_timeRes1275->GetPointX(nPoint)))]){
	vec_tRes_1275[Vov].push_back(std::make_pair ( int(g_timeRes1275->GetPointX(nPoint)), g_timeRes1275->GetPointY(nPoint)));
        vec_tRes_error_1275[Vov].push_back(std::make_pair ( g_timeRes1275->GetErrorX(nPoint), g_timeRes1275->GetErrorY(nPoint)));   
      }
    }   
  } 

  for(std::map<double,TCanvas*>::iterator index = c_timeRes_vs_bar.begin(); index != c_timeRes_vs_bar.end(); index++){
   
    index->second -> cd();
	//Vov = index->first;
    TGraphErrors * graph1 = new TGraphErrors();
    TGraphErrors * graph2 = new TGraphErrors();
    int up = vec_tRes_511[index->first].size();
    for (int i = 0; i < up; i++){
      float min = 9999;
      int minIndex = 0;
      for (int j = 0; j < vec_tRes_511[index->first].size(); j++){
        if (vec_tRes_511[index->first][j].first < min){
	  min = vec_tRes_511[index->first][j].first;
	  minIndex = j;
	}
      }
      graph1 -> SetPoint( i, vec_tRes_511[index->first][minIndex].first, vec_tRes_511[index->first][minIndex].second);
      graph1 -> SetPointError( i, vec_tRes_error_511[index->first][minIndex].first, vec_tRes_error_511[index->first][minIndex].second);
      vec_tRes_511[index->first].erase(vec_tRes_511[index->first].begin() + minIndex);
      vec_tRes_error_511[index->first].erase(vec_tRes_error_511[index->first].begin() + minIndex);
    }
     
    graph1-> SetLineColor(kRed);
    graph1-> SetMarkerColor(kRed);
    graph1 -> SetMarkerStyle(20);
    graph1-> Draw("PL,same");
    outFile -> cd();
    graph1 -> Write(Form("g_timeRes511_vs_bar_Vov%.01f_bestTh",Vov));
	
	
    int up2 = vec_tRes_1275[index->first].size();
    for (int k = 0; k < up2; k++){
      float min = 9999;
      int minIndex = 0;
      for (int l = 0; l < vec_tRes_1275[index->first].size(); l++){
        if (vec_tRes_1275[index->first][l].first < min){
	  min = vec_tRes_1275[index->first][l].first;
	  minIndex = l;
	}
      }
      graph2 -> SetPoint( k, vec_tRes_1275[index->first][minIndex].first, vec_tRes_1275[index->first][minIndex].second);
      graph2 -> SetPointError( k, vec_tRes_error_1275[index->first][minIndex].first, vec_tRes_error_1275[index->first][minIndex].second);
      vec_tRes_1275[index->first].erase(vec_tRes_1275[index->first].begin() + minIndex);
      vec_tRes_error_1275[index->first].erase(vec_tRes_error_1275[index->first].begin() + minIndex);
    }
    graph2-> SetLineColor(kBlue);
    graph2-> SetMarkerColor(kBlue);
    graph2 -> SetMarkerStyle(20);
    graph2-> Draw("PL, same");
    outFile -> cd();
    graph2 -> Write(Form("g_timeRes1275_vs_bar_Vov%.01f_bestTh",Vov));
    

    std::string Label1 = "511 keV Peak";
    TLatex* latex1 = new TLatex(0.55,0.85,Label1.c_str());
    latex1 -> SetNDC();
    latex1 -> SetTextFont(42);
    latex1 -> SetTextSize(0.04);
    latex1 -> SetTextColor(kRed);
    latex1 -> Draw("same"); 
    std::string Label2 = "1275 keV Peak";
    TLatex* latex2 = new TLatex(0.55,0.85-0.04,Label2.c_str());
    latex2 -> SetNDC();
    latex2 -> SetTextFont(42);
    latex2 -> SetTextSize(0.04);
    latex2 -> SetTextColor(kBlue);
    latex2 -> Draw("same");
   
    index->second-> Print(Form("%s/summaryPlots/timeResolution/c_timeRes_vs_bar_Vov%.01f.png",plotDir.c_str(),float(index->first)));
    index->second-> Print(Form("%s/summaryPlots/timeResolution/c_timeRes_vs_bar_Vov%.01f.pdf",plotDir.c_str(),float(index->first)));
  }



int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;

outFile -> Close();

}
                                                                                                                                                                              












