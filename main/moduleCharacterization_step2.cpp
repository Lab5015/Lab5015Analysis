#include "interface/AnalysisUtils.h"
#include "interface/Na22SpectrumAnalyzer.h"
//#include "interface/Na22SpectrumAnalyzerSingleBar.h"
#include "interface/Na22SpectrumAnalyzerSingleBar_TOFHIR2.h"
//#include "interface/Na22SpectrumAnalyzerModule_TOFHIR2.h"
#include "interface/Co60SpectrumAnalyzer_2Peaks.h"
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
#include <stdlib.h>
#include <algorithm>
#include <iterator>

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



Double_t langaufun(Double_t *x, Double_t *par){
  
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
  
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0];
  
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  
  return (par[2] * step * sum * invsq2pi / par[3]);
}






//---- find energy bins
void GetEnergyBins(TH1F *h, std::vector<float> *r, std::map<int, float> & b){

  for(unsigned int i = 1; i < r->size(); i++){
    TH1F *binHisto = new TH1F ( "binHisto", "binHisto", h -> FindBin(r->at(i)) - h->FindBin(r-> at(i-1)), r-> at(i-1), r->at(i));
    int j = 1;
    for (int bin = h->FindBin(r->at(i-1)) ; bin < h -> FindBin(r->at(i))+1 ; bin++){
      binHisto -> SetBinContent( j, h->GetBinContent(bin));
      j++;
    }
    b[i] = binHisto -> GetMean();
    binHisto -> Delete();
  }
}


// ---- Draw  DeltaT plots                                                                                                                                                   
void drawDeltaT(TCanvas *& c, TH1F *histo, TF1 *& fitFunc, std::string xaxis_label, std::string latex_label, std::string drawSame ){

  c->cd();

  // -- first histo
  histo -> SetTitle(Form("; %s #Deltat [ps];entries", xaxis_label.c_str()));
  histo -> Draw(drawSame.c_str());

  float* vals = new float[6];
  FindSmallestInterval(vals,histo,0.68);
  float min = vals[4];
  float max = vals[5];
  float delta = max-min;
  float sigma = 0.5*delta;
  float effSigma = sigma;

  float fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
  float fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;

  fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
  fitFunc -> SetRange(fitXMin, fitXMax);
  histo -> Fit(fitFunc,"QNRSL");
  fitFunc -> SetRange(fitFunc->GetParameter(1)-1.0*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+1.0*fitFunc->GetParameter(2));
  histo -> Fit(fitFunc,"QNRSL");
  fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
  histo -> Fit(fitFunc,"QRSL+");

  fitFunc -> SetLineColor( histo -> GetLineColor() + 1 );
  fitFunc -> SetLineWidth(3);
  fitFunc -> Draw("same");

  histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
  histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-7.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+7.*fitFunc->GetParameter(2));
  //histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-10.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+10.*fitFunc->GetParameter(2));

  TLatex* latex = new TLatex(0.55,0.75,Form("#splitline{#sigma_{%s}^{eff} = %.0f ps}{#sigma_{%s}^{gaus} = %.0f ps}",latex_label.c_str(),effSigma, latex_label.c_str(),fabs(fitFunc->GetParameter(2))));
  if (drawSame == "same")
    latex = new TLatex(0.20,0.85,Form("#splitline{#sigma_{%s}^{eff} = %.0f ps}{#sigma_{%s}^{gaus} = %.0f ps}",latex_label.c_str(),effSigma, latex_label.c_str(),fabs(fitFunc->GetParameter(2))));
  latex -> SetNDC();
  latex -> SetTextFont(42);
  latex -> SetTextSize(0.04);
  latex -> SetTextColor( histo -> GetLineColor() );
  latex -> Draw("same");

}




// ============  ********* MAIN  ************ =============== //
int main(int argc, char** argv) {
    setTDRStyle();
    float cpu[2]{0}, mem[2]={0}, vsz[2]={0}, rss[2]={0};

    gErrorIgnoreLevel = kError;
  
    typedef std::numeric_limits<double> dbl;
    std::cout.precision(dbl::max_digits10);
    if( argc < 2 )
    {
        std::cout << ">>> moduleCharacterization_step2::usage:   " << argv[0] << " configFile.cfg" << std::endl;
        return -1;
    }
  

  
    //--- parse the config file
    CfgManager opts;
    opts.ParseConfigFile(argv[1]);
  
    int debugMode = 0;
    if( argc > 2 ) debugMode = atoi(argv[2]);
  
  
    //--- get parameters
    std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
	plotDir.append("/energySharing");
    //system(Form("rm -r %s", plotDir.c_str())); // questo non va bene se stiamo lavorando in parallelo
    system(Form("mkdir -p %s",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s",plotDir.c_str()));
    system(Form("mkdir -p %s/tot/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/tot/",plotDir.c_str()));
    system(Form("mkdir -p %s/totRatio/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/totRatio/",plotDir.c_str()));
    system(Form("mkdir -p %s/energy/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/energy/",plotDir.c_str()));

	system(Form("mkdir -p %s/energy/singleBars",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/energy/singleBars",plotDir.c_str()));

	system(Form("mkdir -p %s/raw_deltaT",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/raw_deltaT",plotDir.c_str()));
	system(Form("mkdir -p %s/energyCorrelation",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/energyCorrelation",plotDir.c_str()));
	system(Form("mkdir -p %s/timeCorrelation",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/timeCorrelation",plotDir.c_str()));
	system(Form("mkdir -p %s/timeCorrelation/w_REF",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/timeCorrelation/w_REF",plotDir.c_str()));
	system(Form("mkdir -p %s/timeCorrelation/vs_eREF",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/timeCorrelation/vs_eREF",plotDir.c_str()));
	system(Form("mkdir -p %s/timeCorrelation/tREF_vs_eREF",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/timeCorrelation/tREF_vs_eREF",plotDir.c_str()));

	system(Form("mkdir -p %s/CTR_REF",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/CTR_REF",plotDir.c_str()));
	system(Form("mkdir -p %s/CTR_REF/CTR_REF_raw",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/CTR_REF/CTR_REF_raw",plotDir.c_str()));
	system(Form("mkdir -p %s/CTR_REF/CTR_REF_1Bar",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/CTR_REF/CTR_REF_1Bar",plotDir.c_str()));
	system(Form("mkdir -p %s/CTR_REF/CTR_REF_2LR",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/CTR_REF/CTR_REF_2LR",plotDir.c_str()));
	system(Form("mkdir -p %s/CTR_REF/CTR_REF_3LR",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/CTR_REF/CTR_REF_3LR",plotDir.c_str()));

	system(Form("mkdir -p %s/modulePosition",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/modulePosition",plotDir.c_str()));

    system(Form("mkdir -p %s/externalBar",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/externalBar",plotDir.c_str()));

    system(Form("mkdir -p %s/energyRatio/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/energyRatio/",plotDir.c_str()));
    system(Form("mkdir -p %s/t1fine/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/t1fine/",plotDir.c_str()));
    system(Form("mkdir -p %s/qT1/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/qT1/",plotDir.c_str()));
    system(Form("mkdir -p %s/energyRatioCorr/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/energyRatioCorr/",plotDir.c_str()));
    system(Form("mkdir -p %s/energyRatioCorr_totRatioCorr/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/energyRatioCorr_totRatioCorr/",plotDir.c_str()));
    system(Form("mkdir -p %s/totRatioCorr/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/totRatioCorr/",plotDir.c_str()));
    system(Form("mkdir -p %s/phaseCorr/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/phaseCorr/",plotDir.c_str()));
    system(Form("mkdir -p %s/positionCorr/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/positionCorr/",plotDir.c_str()));
    system(Form("mkdir -p %s/CTR_energyRatioCorr/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/CTR_energyRatioCorr/",plotDir.c_str()));
    system(Form("mkdir -p %s/CTR_totRatioCorr/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/CTR_totRatioCorr/",plotDir.c_str()));
    system(Form("mkdir -p %s/CTR_energyRatioCorr_totRatioCorr/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/CTR_energyRatioCorr_totRatioCorr/",plotDir.c_str()));  
    system(Form("mkdir -p %s/CTR_energyRatioCorr_phaseCorr/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/CTR_energyRatioCorr_phaseCorr/",plotDir.c_str()));
    system(Form("mkdir -p %s/CTR_totRatioCorr_phaseCorr/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/CTR_totRatioCorr_phaseCorr/",plotDir.c_str()));
    system(Form("mkdir -p %s/CTR_energyRatioCorr_totRatioCorr_phaseCorr/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/CTR_energyRatioCorr_totRatioCorr_phaseCorr/",plotDir.c_str()));
    system(Form("mkdir -p %s/CTR_energyRatioCorr_phaseCorr_posCorr/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/CTR_energyRatioCorr_phaseCorr_posCorr/",plotDir.c_str()));
    system(Form("mkdir -p %s/CTR_totRatioCorr_phaseCorr_posCorr/",plotDir.c_str()));
    system(Form("cp -n /eos/home-f/ftonetto/www/index.php %s/CTR_totRatioCorr_phaseCorr_posCorr/",plotDir.c_str()));  
  

  
    std::vector<std::string> LRLabels;
    LRLabels.push_back("L");
    LRLabels.push_back("R");
    LRLabels.push_back("L-R");
  
    std::vector<std::string> PPLabels;
    PPLabels.push_back("PRE");
    PPLabels.push_back("POST");
    PPLabels.push_back("PREPOST");
    

    std::vector<float> Vov = opts.GetOpt<std::vector<float> >("Plots.Vov");
    std::vector<int> energyMins = opts.GetOpt<std::vector<int> >("Plots.energyMins");
    std::vector<int> energyMaxs = opts.GetOpt<std::vector<int> >("Plots.energyMaxs");
  
    std::map<float,int> map_energyMins;
    std::map<float,int> map_energyMaxs;
    for(unsigned int ii = 0; ii < Vov.size(); ++ii)
    {
        map_energyMins[Vov[ii]] = energyMins[ii];
        map_energyMaxs[Vov[ii]] = energyMaxs[ii];
    }
  
    int useTrackInfo = opts.GetOpt<int>("Input.useTrackInfo");
  
    // -- read minimum energy for each bar from file
    std::string minEnergiesFileName = opts.GetOpt<std::string>("Cuts.minEnergiesFileName");
    std::map < std::pair<int, float>, float> minE; 
    std::cout << minEnergiesFileName << std::endl;
    if( minEnergiesFileName != "" ) {
        std::ifstream minEnergiesFile;
        minEnergiesFile.open(minEnergiesFileName);
        std::string line;
        int bar;
        float ov;
        float value;
        while(minEnergiesFile.good() ){
	        getline(minEnergiesFile, line);
	        std::istringstream ss(line);
	        ss >> bar >> ov >> value; 
	        minE[std::make_pair(bar,ov)] = value; 
	        //std::cout << "minEnergies:   bar " << bar <<  "   Vov " << ov << "   E " << minE[std::make_pair(bar,ov)] <<std::endl;
	    }
    }
    else {
      for(unsigned int iBar = 0; iBar < 16; ++iBar)
	        for(unsigned int ii = 0; ii < Vov.size(); ++ii)
	            minE[std::make_pair(iBar, Vov[ii])] = map_energyMins[Vov[ii]];
    }



    //--- open files
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
    while( (object = next()) ) {
        std::string name(object->GetName());
        std::vector<std::string> tokens = GetTokens(name,'_');
        std::size_t found;
     
        found = name.find("data_");
        //tree
        if( found!=std::string::npos ){
            std::string label(Form("%s_%s_%s",tokens[1].c_str(),tokens[2].c_str(),tokens[3].c_str()));
            trees[label] = (TTree*)( inFile->Get(name.c_str()) );
        }
        found = name.find("h1_energy_b");
        if( found!=std::string::npos ) {
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
        }
    }
    std::sort(stepLabels.begin(),stepLabels.end());
    stepLabels.erase(std::unique(stepLabels.begin(),stepLabels.end()),stepLabels.end());
  
  
    //--- define histograms
    std::string outFileName = opts.GetOpt<std::string>("Output.outFileNameStep2");
    TFile* outFile = TFile::Open(outFileName.c_str(),"RECREATE");
    outFile->cd();

    std::map<double,TH1F*> h1_energyRatio;
    std::map<double,TH1F*> h1_totRatio;
    std::map<double,TH1F*> h1_t1fineMean;
    std::map<double,TH1F*> h1_qT1Mean;
      
    // -- raw deltaT
    std::map<double,TH1F*> h1_deltaT_raw;
    std::map<double,TH1F*> h1_deltaT;

    // ---  deltaT of bars, i.e. deltaT = t_ave(i) - t_ave(i+1)
	std::map<double,TH1F*> h1_deltaT_bars;
	std::map<double,TH1F*> h1_deltaT_bars_w;
	std::map<double,TH1F*> h1_deltaT_bars_w_raw;
	std::map<double,TH1F*> h1_energyRatio_bars_F;
	std::map<double,TH1F*> h1_energyRatio_bars_B;
	std::map<double,TH2F*> h2_deltaT_energyRatio_bars;
	std::map<double,TProfile*> p1_deltaT_energyRatio_bars;
	std::map<double,TH2F*> h2_deltaT_energyRatio_bars_w;
	std::map<double,TProfile*> p1_deltaT_energyRatio_bars_w;
	

	std::map<double,TProfile*> p1_tL_eL;
	std::map<double,TProfile*> p1_tR_eR;

    // -- energy correlation plots
	std::map<double, TH2F*> h2_energy_correlation;
	std::map<double, TProfile*> p1_energy_correlation;
	// -- time correlation
	std::map<double, TH2F*> h2_timeDiff_correlation;
	std::map<double, TProfile*> p1_timeDiff_correlation;

	std::map<double, TH2F*> h2_tL_REFvsEL;
	std::map<double, TH2F*> h2_tR_REFvsER;

	std::map<double, TProfile*> p1_tL_REFvsEL;
	std::map<double, TProfile*> p1_tR_REFvsER;

	std::map<double, TH2F*> h2_tL_wREFvsEL;
	std::map<double, TH2F*> h2_tR_wREFvsER;

	std::map<double, TProfile*> p1_tL_wREFvsEL;
	std::map<double, TProfile*> p1_tR_wREFvsER;

	std::map<double, TH2F*> h2_tL_REFvsE_REF;
	std::map<double, TH2F*> h2_tR_REFvsE_REF;

    // -- energy and/or tot corr  
    std::map<double,TProfile*> p1_deltaT_vs_energyRatio;
    std::map<double,TProfile*> p1_deltaT_vs_totRatio;
    std::map<double,TProfile*> p1_deltaT_energyRatioCorr_vs_totRatio;
  
    // -- phase corr
    std::map<double,TProfile*> p1_deltaT_energyRatioCorr_vs_t1fineMean;
    std::map<double,TProfile*> p1_deltaT_totRatioCorr_vs_t1fineMean;
    std::map<double,TProfile*> p1_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean;

    std::map<double,TH2F*> h2_deltaT_energyRatioCorr_vs_t1fineMean;
    std::map<double,TH2F*> h2_deltaT_totRatioCorr_vs_t1fineMean;
    std::map<double,TH2F*> h2_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean;
 
    // -- position corr 
    std::map<double,TProfile*> p1_deltaT_energyRatioCorr_vs_posX;
    std::map<double,TProfile*> p1_deltaT_totRatioCorr_vs_posX;
    std::map<double,TProfile*> p1_deltaT_energyRatioCorr_totRatioCorr_vs_posX;

    std::map<double,TH2F*> h2_deltaT_vs_totRatio;
    std::map<double,TProfile*> p1_deltaT_totRatioCorr_vs_totRatio;
    std::map<double,TH2F*> h2_deltaT_totRatioCorr_vs_totRatio;
    std::map<double,TH2F*> h2_deltaT_energyRatioCorr_vs_totRatio;

    // -- corrected deltaT
    std::map<double,TH1F*> h1_deltaT_energyRatioCorr;
    std::map<double,TH1F*> h1_deltaT_totRatioCorr;
    std::map<double,TH1F*> h1_deltaT_energyRatioCorr_totRatioCorr;
    std::map<double,TH1F*> h1_deltaT_energyRatioPhaseCorr;
	std::map<double,TH1F*> h1_deltaT_energyRatioPhaseCorr_bars;
    std::map<double,TH1F*> h1_deltaT_totRatioPhaseCorr;
    std::map<double,TH1F*> h1_deltaT_energyRatioCorr_totRatioCorr_phaseCorr;
    std::map<double,TH1F*> h1_deltaT_energyRatioPhasePosCorr;
    std::map<double,TH1F*> h1_deltaT_totRatioPhasePosCorr;

    std::map<double,TProfile*> p1_deltaT_totRatioPhaseCorr_vs_totRatio;
    std::map<double,TH2F*> h2_deltaT_totRatioPhaseCorr_vs_totRatio;

	std::map<double, TProfile*> p1_deltaTave_vs_energyRatio;
	std::map<double,TH2F*> h2_deltaTave_vs_energyRatio;

	std::map<double, TProfile*> p1_enL_enR_ave;
	std::map<double, TH2F*> h2_enL_enR_ave;

	// --- deltaT DUT - REF
	std::map<double, TH1F*> h1_deltaT_raw_REF;
	std::map<double, TH1F*> h1_deltaT_w_REF;
	std::map<double, TH1F*> h1_deltaT_corr_REF;
	std::map<double, TH1F*> h1_deltaT_LR_corr_REF;
	std::map<double, TH1F*> h1_deltaT_LR_corr_REF_w;

	// --- external bar
    std::map<int, TH1F*> h1_energyRatio_REF;
	std::map<int, TF1*>  fitFunc_energyRatio_REF;
	std::map<int, TH2F*> h2_deltaT_vs_energyRatio_REF;
	std::map<int, TProfile*> p1_deltaT_vs_energyRatio_REF;
	std::map<int, TF1*> fitFunc_energyRatioCorrection_REF;
	std::map<int, TH1F*> h1_deltaT_eneryRatioCorr_REF;
    
	std::map<int,TH1F*> h1_t1fineMean_REF;
    std::map<int,TH2F*> h2_deltaT_energyRatioCorr_vs_t1fineMean_REF;
    std::map<int,TProfile*> p1_deltaT_energyRatioCorr_vs_t1fineMean_REF;
	std::map<int, TH1F*> h1_deltaT_eneryRatioCorr_pasheCorr_REF;

    // this ranges are used for single hits selection
    std::map<std::string, std::map<int, std::vector<float>*> > ranges; //ranges[LRlabel][index]

	std::map<double ,std::vector<double>* > energyRatio_ranges;

	//std::map<int,std::vector<float> >  my_ranges;

	//this ranges are used for double hits selection
	std::map<int,std::vector<float>*> ranges_doubleHits; 

    std::map<std::string, std::map<int, std::map<std::string,std::pair<float,float> > > > peaks;	//peaks[LRlabel][index][energyPeak]
    std::map<std::string, std::map<int, std::map<int,float> > > energyBin; // energyBin[LRlabel][index]

    std::map<int,TF1*>  f_langaus; // f_langaus[index]
    std::map<int,TF1*>  f_gaus; // f_gaus[index]
    std::map<int,TF1*>  f_landau; // f_landau[index]
    std::map<int,TF1*>  f_landau1; 	
    std::map<float,int>  Vov_LandauMin; //Vov_LandauMin[Vov]
    Vov_LandauMin[1.0] = 0;
    Vov_LandauMin[1.5] = 50;
    Vov_LandauMin[2.0] = 100;
    Vov_LandauMin[2.5] = 180;
    Vov_LandauMin[3.0] = 250;
    Vov_LandauMin[3.5] = 300;
    Vov_LandauMin[6.0] = 300;
  
  
    //--- get plot settings
    TCanvas* c;
    TCanvas* c2;
    float* vals = new float[6];
	float* my_vals = new float[6];
    TLatex* latex;
    TH1F* histo;
    TProfile* prof;
    TH2F* h2;
  

	std::map<double, double> my_cut; //my_cut[index2]

    //------------------
    //--- draw 1st plots
    std::string source = opts.GetOpt<std::string>("Input.sourceName");
    std::string Na22 = "Na22";
    std::string Na22SingleBar = "Na22SingleBar";
    std::string Co60 = "Co60";
    std::string Co60SumPeak = "Co60SumPeak";	
    std::string Laser = "Laser";	
    std::string TB = "TB";
    std::string keepAll = "keepAll";
    std::vector<int> barList = opts.GetOpt<std::vector<int> >("Plots.barList");// list of bars to be analyzed read from cfg

    for(auto stepLabel : stepLabels) {   
        std::cout<<"##########\nstepLabel: "<<stepLabel<<"\n"<<std::endl;
		//TrackProcess(cpu, mem, vsz, rss);
        
        float Vov = map_Vovs[stepLabel];
        float vth1 = map_ths[stepLabel];
        std::string VovLabel(Form("Vov%.2f",Vov));
        std::string thLabel(Form("th%02.0f",vth1));
      
        // -----------
        // --- external bar
        // -- draw energy
        c = new TCanvas(Form("c_energy_external_Vov%.2f_th%02.0f",Vov,vth1), Form("c_energy_external_Vov%.2f_th%02.0f",Vov,vth1));
        gPad -> SetLogy();
        histo = (TH1F*)( inFile->Get(Form("h1_energy_external_barL-R_Vov%.2f_th%02.0f", Vov, vth1)));
        if ( histo ) {
	        histo -> SetTitle(";energy [a.u.];entries");
	        histo -> SetLineColor(kRed);                                                                                                                                              
	        histo -> SetLineWidth(2);                                                                                                                                                 
	        histo -> Draw();  
	  
	        TF1 *ftemp = histo->GetFunction( ((histo->GetListOfFunctions()->FirstLink())->GetObject())->GetName() );
	        if (ftemp != NULL){
	            ftemp->SetNpx(1000);
	            ftemp->SetLineWidth(2);
	            ftemp->SetLineColor(1);
	            ftemp->Draw("same");
				float max_ext = ftemp->GetParameter(1);
				TLine* line_ext = new TLine(0.85*max_ext, 0,0.85*max_ext,histo->GetMaximum());
				line_ext->SetLineColor(kBlue);
				line_ext->SetLineWidth(2);
				line_ext->SetLineStyle(2);
				line_ext->Draw("same");
				line_ext = new TLine(400, 0,400,histo->GetMaximum());
				line_ext->SetLineColor(kBlue);
				line_ext->SetLineWidth(2);
				line_ext->SetLineStyle(2);
				line_ext->Draw("same");
				
	        }
	        c -> Print(Form("%s/energy/c_energy_external__Vov%.2f_th%02.0f.png",plotDir.c_str(), Vov, vth1));
	        c -> Print(Form("%s/energy/c_energy_external__Vov%.2f_th%02.0f.pdf",plotDir.c_str(), Vov, vth1));
	        delete c;
	        delete ftemp;
	    }

		// --- timeREF vs energyREF (left and right ch )
		c = new TCanvas(Form("c_tREFvseREF_extBar_Vov%.2f_th%02.0f_D",Vov,vth1),Form("c_tREFvseREF_extBar_Vov%.2f_th%02.0f_D",Vov,vth1));
		h2 = (TH2F*)( inFile->Get(Form("h2_tREF_L_vs_eREF_Vov%.2f_th%02.0f", Vov, vth1)));
		h2 -> GetYaxis()->SetRangeUser(-500, 500);
        h2 -> SetTitle(Form(";energy_{REF} [a.u.]; time_{REF} - time_{REF}^{L} [ps]"));
		h2->Draw("colz");
                        
		prof = (TProfile*)( inFile->Get(Form("p1_tREF_L_vs_eREF_Vov%.2f_th%02.0f", Vov, vth1)));
		prof -> SetMarkerSize(0.4);
		prof -> Draw("psame");
                    
		latex = new TLatex(0.40,0.8,Form("#splitline{ext. bar}{V_{OV} = %.2f V, th. = %d DAC}",Vov,int(vth1)));
		latex -> SetNDC();
		latex -> SetTextFont(42);
		latex -> SetTextSize(0.04);
		latex -> SetTextColor(kRed);
		latex -> Draw("same");
		c->Print(Form("%s/timeCorrelation/tREF_vs_eREF/c_tREF_vs_eREF_extBar_%s_D.png",plotDir.c_str(),stepLabel.c_str()));
		c->Print(Form("%s/timeCorrelation/tREF_vs_eREF/c_tREF_vs_eREF_extBar_%s_D.pdf",plotDir.c_str(),stepLabel.c_str()));
		delete c;

		c = new TCanvas(Form("c_tREFvseREFr_extBar_Vov%.2f_th%02.0f_D",Vov,vth1),Form("c_tREFvseREFr_extBar_Vov%.2f_th%02.0f_D",Vov,vth1));
		h2 = (TH2F*)( inFile->Get(Form("h2_tREF_R_vs_eREF_Vov%.2f_th%02.0f", Vov, vth1)));
		h2 -> GetYaxis()->SetRangeUser(-500, 500);
        h2 -> SetTitle(Form(";energy_{REF} [a.u.]; time_{REF} - time_{REF}^{R} [ps]"));
		h2->Draw("colz");
                        
		prof = (TProfile*)( inFile->Get(Form("p1_tREF_R_vs_eREF_Vov%.2f_th%02.0f", Vov, vth1)));
		prof -> SetMarkerSize(0.4);
		prof -> Draw("psame");
                    
		latex = new TLatex(0.40,0.8,Form("#splitline{ext. bar}{V_{OV} = %.2f V, th. = %d DAC}",Vov,int(vth1)));
		latex -> SetNDC();
		latex -> SetTextFont(42);
		latex -> SetTextSize(0.04);
		latex -> SetTextColor(kRed);
		latex -> Draw("same");
		c->Print(Form("%s/timeCorrelation/tREF_vs_eREF/c_tREFr_vs_eREF_extBar_%s_D.png",plotDir.c_str(),stepLabel.c_str()));
		c->Print(Form("%s/timeCorrelation/tREF_vs_eREF/c_tREFr_vs_eREF_extBar_%s_D.pdf",plotDir.c_str(),stepLabel.c_str()));
		delete c;

		// trackInfo
		if(useTrackInfo && !opts.GetOpt<std::string>("Coincidence.status").compare("no") ){
		    c = new TCanvas(Form("c_trackInfo_noCoin_bar06_Vov%.2f_th%02.0f_D",Vov,vth1),Form("c_trackInfo_noCoin_bar06_Vov%.2f_th%02.0f_D",Vov,vth1));
			h2 = (TH2F*)(inFile->Get(Form("h2_xy_double_noCoin_bar06_Vov%.2f_th%02.0f",Vov,vth1))); 
			if(!h2){
				std::cout<<Form("h2_xy_double_noCoin_bar06_Vov%.2f_th%02.0f",Vov,vth1) <<" not found"<<std::endl;
			    continue;
			}
			h2 -> GetXaxis() -> SetRangeUser(-50,50);
			h2 -> GetYaxis() -> SetRangeUser(-50,50);
			h2 -> SetTitle(Form(";x [mm];y [mm]"));
			h2 -> Draw("colz");

			latex = new TLatex(0.40,0.85,Form("#splitline{bars 06-07}{V_{OV} = %.2f V, th. = %d DAC}",Vov,int(vth1)));
			latex -> SetNDC();
			latex -> SetTextFont(42);
			latex -> SetTextSize(0.04);
			latex -> SetTextColor(kRed);
			latex -> Draw("same");
			c -> Print(Form("%s/modulePosition/c_trackInfo_noCoin_bar06_D_%s.png",plotDir.c_str(),stepLabel.c_str()));
			c -> Print(Form("%s/modulePosition/c_trackInfo_noCoin_bar06_D_%s.pdf",plotDir.c_str(),stepLabel.c_str()));
			delete c;
			delete latex;
		}
       if(useTrackInfo && !opts.GetOpt<std::string>("Coincidence.status").compare("yes") ){
		    c = new TCanvas(Form("c_trackInfo_bar06_Vov%.2f_th%02.0f_D",Vov,vth1),Form("c_trackInfo_bar06_Vov%.2f_th%02.0f_D",Vov,vth1));
			h2 = (TH2F*)(inFile->Get(Form("h2_xy_double_bar06_Vov%.2f_th%02.0f",Vov,vth1))); 
			if(!h2){
				std::cout<<Form("h2_xy_double_bar06_Vov%.2f_th%02.0f",Vov,vth1) <<" not found"<<std::endl;
			    continue;
			}
			h2 -> GetXaxis() -> SetRangeUser(-50,50);
			h2 -> GetYaxis() -> SetRangeUser(-50,50);
			h2 -> SetTitle(Form(";x [mm];y [mm]"));
			h2 -> Draw("colz");

			latex = new TLatex(0.40,0.85,Form("#splitline{bars 06-07}{V_{OV} = %.2f V, th. = %d DAC}",Vov,int(vth1)));
			latex -> SetNDC();
			latex -> SetTextFont(42);
			latex -> SetTextSize(0.04);
			latex -> SetTextColor(kRed);
			latex -> Draw("same");
			c -> Print(Form("%s/modulePosition/c_trackInfo_bar06_D_%s.png",plotDir.c_str(),stepLabel.c_str()));
			c -> Print(Form("%s/modulePosition/c_trackInfo_bar06_D_%s.pdf",plotDir.c_str(),stepLabel.c_str()));
			delete c;
			delete latex;

			c = new TCanvas(Form("c_trackInfo_externalBar_Vov%.2f_th%02.0f_D",Vov,vth1),Form("c_trackInfo_externalBar_Vov%.2f_th%02.0f_D",Vov,vth1));
			h2 = (TH2F*)(inFile->Get(Form("h2_trackInfo_REF_Vov%.2f_th%02.0f",Vov,vth1))); 
			if(!h2){
				std::cout<<Form("h2_trackInfo_REF_%.2f_th%02.0f",Vov,vth1) <<" not found"<<std::endl;
			    continue;
			}
			h2 -> GetXaxis() -> SetRangeUser(-50,50);
			h2 -> GetYaxis() -> SetRangeUser(-50,50);
			h2 -> SetTitle(Form(";x [mm];y [mm]"));
			h2 -> Draw("colz");

			latex = new TLatex(0.40,0.85,Form("#splitline{external bar}{V_{OV} = %.2f V, th. = %d DAC}",Vov,int(vth1)));
			latex -> SetNDC();
			latex -> SetTextFont(42);
			latex -> SetTextSize(0.04);
			latex -> SetTextColor(kRed);
			latex -> Draw("same");
			c -> Print(Form("%s/modulePosition/c_trackInfo_externalBar_%s.png",plotDir.c_str(),stepLabel.c_str()));
			c -> Print(Form("%s/modulePosition/c_trackInfo_externalBar_%s.pdf",plotDir.c_str(),stepLabel.c_str()));
			delete c;
			delete latex;

			c = new TCanvas(Form("c_trackInfo_externalBar_cuts_Vov%.2f_th%02.0f_D",Vov,vth1),Form("c_trackInfo_externalBar_cuts_Vov%.2f_th%02.0f_D",Vov,vth1));
			h2 = (TH2F*)(inFile->Get(Form("h2_trackInfo_REF_cuts_Vov%.2f_th%02.0f",Vov,vth1))); 
			if(!h2){
				std::cout<<Form("h2_trackInfo_REF_cuts_%.2f_th%02.0f",Vov,vth1) <<" not found"<<std::endl;
			    continue;
			}
			h2 -> GetXaxis() -> SetRangeUser(-50,50);
			h2 -> GetYaxis() -> SetRangeUser(-50,50);
			h2 -> SetTitle(Form(";x [mm];y [mm]"));
			h2 -> Draw("colz");

			latex = new TLatex(0.40,0.85,Form("#splitline{external bar}{V_{OV} = %.2f V, th. = %d DAC}",Vov,int(vth1)));
			latex -> SetNDC();
			latex -> SetTextFont(42);
			latex -> SetTextSize(0.04);
			latex -> SetTextColor(kRed);
			latex -> Draw("same");
			c -> Print(Form("%s/modulePosition/c_trackInfo_externalBar_cuts_%s.png",plotDir.c_str(),stepLabel.c_str()));
			c -> Print(Form("%s/modulePosition/c_trackInfo_externalBar_cuts_%s.pdf",plotDir.c_str(),stepLabel.c_str()));
			delete c;
			delete latex;
		}
        //--------------------------------------------------------
        // --- loop over bars
        for(int iBar = 0; iBar < 16; ++iBar) {
	
	        bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
	        if (!barFound) continue;
	
	        int index( (10000*int(Vov*100.)) + (100*vth1) + iBar );
	
	        // -- loop over L, R, LR
	        for(auto LRLabel : LRLabels ) {
	  
	            //label histo
	            std::string label(Form("bar%02d%s_%s",iBar,LRLabel.c_str(),stepLabel.c_str()));
	            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d%s}{V_{OV} = %.2f V, th. = %d DAC}",iBar,LRLabel.c_str(),Vov,int(vth1)));
	            if (LRLabel == "L-R") { 
	                latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	            }
	            latex -> SetNDC();
	            latex -> SetTextFont(42);
	            latex -> SetTextSize(0.04);
	            latex -> SetTextColor(kRed);
          
	  
	            // -- qfine and tot only per L, R
	            if (LRLabel == "R" || LRLabel == "L") {
	                // -- draw ToT
	                c = new TCanvas(Form("c_tot_%s",label.c_str()),Form("c_tot_%s",label.c_str()));
	                gPad -> SetLogy();
	      
	                histo = (TH1F*)( inFile->Get(Form("h1_tot_%s",label.c_str())) );
	                if (histo){
		                histo -> SetTitle(";ToT [ns];entries");
		                histo -> SetLineColor(kRed);
		                histo -> Draw();
		                // TLine* line_totAcc1 = new TLine(cut_totAcc[chID][Vov],histo->GetMinimum(),cut_totAcc[chID][Vov],histo->GetMaximum());
		                // line_totAcc1 -> SetLineColor(kBlack);
		                // line_totAcc1 -> Draw("same");
		                latex -> Draw("same");      
		                histo -> Write();
		                c -> Print(Form("%s/tot/c_tot__%s.png",plotDir.c_str(),label.c_str()));
		                c -> Print(Form("%s/tot/c_tot__%s.pdf",plotDir.c_str(),label.c_str()));
		                delete c;
		                }
	            }
	  
	            // -- draw energy
	            c = new TCanvas(Form("c_energy_%s",label.c_str()),Form("c_energy_%s",label.c_str()));
	            gPad -> SetLogy();
	  
	            histo = (TH1F*)( inFile->Get(Form("h1_energy_%s",label.c_str())) );      
	            if( !histo ) continue;
	            histo -> SetTitle(";energy [a.u.];entries");
	            histo -> SetLineColor(kRed);
	            histo -> SetLineWidth(2);
	            histo -> Draw();
	  
	  
	            // --- look for peaks and define energy ranges
	            if( source.compare(Na22) && source.compare(Na22SingleBar) && source.compare(Co60) && source.compare(Co60SumPeak) && source.compare(Laser) && source.compare(TB) && source.compare(keepAll) )
	            {
	                std::cout << " Source not found !!! " << std::endl;
	                return(0);
	            }
	  
	            ranges[LRLabel][index] = new std::vector<float>;
	  
	            // --- Na22 or Co60 spectrum
	            if(!source.compare(Na22)  || !source.compare(Na22SingleBar) ||  !source.compare(Co60) ){
	                std::string firstPeak = "";
	                if (!source.compare(Na22)){
		                //hardcoding histo range for central bar
		                ranges[LRLabel][index]->push_back(30);
		                ranges[LRLabel][index]->push_back(950);
		                if (iBar == 8){
		                    ranges[LRLabel][index]->clear();
		                    ranges[LRLabel][index]->push_back(150);
		                    ranges[LRLabel][index]->push_back(950);
		                }	  
		                peaks[LRLabel][index] = Na22SpectrumAnalyzer(histo,ranges[LRLabel][index]);
		                firstPeak = "0.511 MeV";
		            }
	                if (!source.compare(Na22SingleBar)) {
		                //peaks[LRLabel][index] = Na22SpectrumAnalyzerSingleBar(histo,ranges[LRLabel][index]);//
		                // !!!!! fare in modo di avere un singolo codice Na22SpectrumAnalyzer, eventualmente con delle flag...
		                peaks[LRLabel][index] = Na22SpectrumAnalyzerSingleBar_TOFHIR2(histo,ranges[LRLabel][index]);//run330, barInArray
		                //peaks[LRLabel][index] = Na22SpectrumAnalyzerModule_TOFHIR2(histo,ranges[LRLabel][index]);//run555, bar read out by two modules, "wirelessbar"
		                firstPeak = "0.511 MeV";	
	                }
	                if (!source.compare(Co60)) {
		                peaks[LRLabel][index] = Co60SpectrumAnalyzer_2Peaks(histo,ranges[LRLabel][index]);
		                firstPeak = "1.173 MeV";
	                }
	      
	                if (peaks[LRLabel][index][firstPeak].first > 0){
		                //histo -> GetXaxis() -> SetRangeUser(0.,5.*peaks[LRLabel][index][firstPeak].first);
		                histo -> GetXaxis() -> SetRangeUser(0.,1000.);
	                }
	      
	                if (peaks[LRLabel][index][firstPeak].first== -9999){
		                histo -> GetXaxis() -> SetRangeUser(0., map_energyMaxs[Vov]); 
		                peaks[LRLabel].erase(index);
		                ranges[LRLabel].erase(index);
		                if (!source.compare(Na22) || !source.compare(Na22SingleBar) ) peaks[LRLabel][index]["0.511 MeV"].first = -10;
		                if (!source.compare(Na22) || !source.compare(Na22SingleBar) ) peaks[LRLabel][index]["1.275 MeV"].first = -10;
		                if (!source.compare(Na22SingleBar))                           peaks[LRLabel][index]["1.786 MeV"].first = -10;
		                if (!source.compare(Co60))                                    peaks[LRLabel][index]["1.173 MeV"].first = -10;
		                if (!source.compare(Co60))                                    peaks[LRLabel][index]["1.332 MeV"].first = -10;
		                if (!source.compare(Co60))                                    peaks[LRLabel][index]["2.505 MeV"].first = -10;
	                }
	      
	                if(peaks[LRLabel][index][firstPeak].first != -10){
		                GetEnergyBins(histo, ranges[LRLabel][index], energyBin[LRLabel][index]);
	                }          
	            } // end if Na22/Co60
                
				// -- if Co60 sum peak or laser, we don't use the spectrum analyzers - just gaussian fit to the energy peak.
	            if( !source.compare(Co60SumPeak) || !source.compare(Laser) ){ 
	                histo -> GetXaxis() -> SetRangeUser(map_energyMins[Vov],map_energyMaxs[Vov]);
	                float max = FindXMaximum(histo,map_energyMins[Vov],map_energyMaxs[Vov]);
	      
	                TF1* fitFunc = new TF1 ( Form("func_%d",index), "gaus(0)", max-2*histo->GetRMS(), max+2*histo->GetRMS() );
	                fitFunc -> SetLineColor(kBlack);
	                fitFunc -> SetLineWidth(2);
	                histo -> Fit( fitFunc, "NQR");
	                fitFunc -> SetRange(fitFunc->GetParameter(1) - fitFunc->GetParameter(2)*5, fitFunc->GetParameter(1) + fitFunc->GetParameter(2)*5 );
	                histo -> Fit( fitFunc, "QRS+");
	      
	                // -- questa condizione er solo per il laser... puo' andare bene tenerla anche per 60Co?
	                //if ((fitFunc->GetParameter(1) - fitFunc->GetParameter(2)*5)>0) 
	      
	                //ranges[LRLabel][index] -> push_back( fitFunc->GetParameter(1) - fitFunc->GetParameter(2)*5);
	                ranges[LRLabel][index] -> push_back( fitFunc->GetParameter(1) - std::max(fitFunc->GetParameter(2)*5, histo->GetBinWidth(1)));
	                // else 
	                //   ranges[LRLabel][index] -> push_back(0);
	                //ranges[LRLabel][index] -> push_back( fitFunc->GetParameter(1) + fitFunc->GetParameter(2)*5);
	                ranges[LRLabel][index] -> push_back( fitFunc->GetParameter(1) + std::max(fitFunc->GetParameter(2)*5, histo->GetBinWidth(1)));
	      
	                for(auto range: (*ranges[LRLabel][index])){
		                float yval = std::max(10., histo->GetBinContent(histo->FindBin(range)));
		                TLine* line = new TLine(range,3.,range, yval);
		                line -> SetLineWidth(1);
		                line -> SetLineStyle(7);
		                line -> Draw("same");
	                }
	      
	                GetEnergyBins(histo, ranges[LRLabel][index], energyBin[LRLabel][index]);
	            }// end Co60SumPeak or laser
	  
	        	
	  
	            // -- keep all events within energyMin and energyMax
	            if( !source.compare(keepAll) ) { 
	                ranges[LRLabel][index] -> push_back( map_energyMins[Vov] );
	                //ranges[LRLabel][index] -> push_back( map_energyMins[Vov] );
	                ranges[LRLabel][index] -> push_back( map_energyMaxs[Vov] );
	      
	                for(auto range: (*ranges[LRLabel][index])){
		                float yval = std::max(10., histo->GetBinContent(histo->FindBin(range)));
		                TLine* line = new TLine(range,3.,range, yval);
		                line -> SetLineWidth(1);
		                line -> SetLineStyle(7);
		                line -> Draw("same");
		            }
	      
	                GetEnergyBins(histo, ranges[LRLabel][index], energyBin[LRLabel][index]);
	            }// end keepAll
	  
	  
	            // -- if MIP peak, we don't use the spectrum analyzers - just langaus fit to the energy peak.
	            if(!source.compare(TB)){ 
	                /*
	                //TF1* f_langaus = new TF1("f_langaus", langaufun, 300.,1000.,4);
	                f_langaus[index] = new TF1(Form("fit_energy_bar%02d_Vov%.2f_vth1_%02.0f",iBar,Vov,vth1),langaufun,Vov_LandauMin[Vov],1000.,4);
	                f_langaus[index] -> SetParameters(f_landau->GetParameter(2),f_landau->GetParameter(1),histo->Integral(histo->FindBin(300.),histo->FindBin(1000.))*histo->GetBinWidth(1),10.);
	                f_langaus[index] -> SetLineColor(kBlack);
	                f_langaus[index] -> SetLineWidth(2);
	                histo -> Fit(f_langaus[index],"QNRS+");
	                f_langaus[index] -> Draw("same");
	    
	                ranges[LRLabel][index] -> push_back( 0.80*f_langaus[index]->GetParameter(1));
	                ranges[LRLabel][index] -> push_back( histo -> GetBinCenter(500) );
	                */
	    
	                if( opts.GetOpt<int>("Channels.array") == 1){
	                    histo->GetXaxis()->SetRangeUser(minE[std::make_pair(iBar, Vov)], 950);
	                }
	                if( opts.GetOpt<int>("Channels.array") == 0){
	                    histo->GetXaxis()->SetRangeUser(minE[std::make_pair(iBar, Vov)], 950);
	                }
	                float max = histo->GetBinCenter(histo->GetMaximumBin());
	                histo->GetXaxis()->SetRangeUser(0,1024);
	    
	                f_gaus[index] = new TF1(Form("fit_energy_bar%02d%s_Vov%.2f_vth1_%02.0f",iBar,LRLabel.c_str(),Vov,vth1), "gaus", max-50, max+50);
	                f_gaus[index]->SetParameters(histo->GetMaximumBin(), max, 70);
	                histo->Fit(f_gaus[index], "QRS");
	                f_gaus[index]->SetRange(f_gaus[index]->GetParameter(1)-f_gaus[index]->GetParameter(2), f_gaus[index]->GetParameter(1)+f_gaus[index]->GetParameter(2));
	                histo->Fit(f_gaus[index], "QRS");
	                f_gaus[index] -> SetLineColor(kBlue);
	                f_gaus[index] -> SetLineWidth(2);
	                f_gaus[index] -> SetLineStyle(2);
	                //f_gaus[index] -> Draw("same");
	    
	                //ranges[LRLabel][index] -> push_back( 0.80*f_gaus[index]->GetParameter(1));
	                //ranges[LRLabel][index] -> push_back( histo -> GetBinCenter(700) ); // to avoid sturation
	    
	                f_landau[index] = new TF1(Form("f_landau_bar%02d%s_Vov%.2f_vth1_%02.0f", iBar,LRLabel.c_str(),Vov,vth1),"[0]*TMath::Landau(x,[1],[2])", 0,1000.);
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
	    
	                if( f_landau[index]->GetNDF() >0 && f_landau[index]->GetParameter(1) > minE[std::make_pair(iBar, Vov)] &&  
		                (f_landau[index]->GetParameter(1) - 2.0 * std::abs(f_landau[index]->GetParameter(2))) >=  minE[std::make_pair(iBar, Vov)] &&
		                (f_landau[index]->GetParameter(1) - 2.0 * std::abs(f_landau[index]->GetParameter(2))) < 950) {
	                    //ranges[LRLabel][index] -> push_back( 0.75*f_landau[index]->GetParameter(1));
	                    //ranges[LRLabel][index] -> push_back( 0.60*f_landau[index]->GetParameter(1));
	                    ranges[LRLabel][index] -> push_back( f_landau[index]->GetParameter(1) - 2.0 * std::abs(f_landau[index]->GetParameter(2)));
	                }
	                else  ranges[LRLabel][index] -> push_back( minE[std::make_pair(iBar, Vov)] ); // 
	    
	                if ( LRLabel=="L-R" && int(vth1)==10) std::cout << iBar << "  " << Vov  << "  " << ranges[LRLabel][index] ->at(0) <<std::endl;
	    
	    
	                //ranges[LRLabel][index] -> push_back( std::min(f_landau[index]->GetParameter(1)*2.0, 940.)); // tight selection around the MIP peak
	                ranges[LRLabel][index] -> push_back( 940 ); // use the entire mip spectrum
	    
	                for(auto range: (*ranges[LRLabel][index])){
	                    TLine* line = new TLine(range,0.,range, histo->GetMaximum());
	                    line -> SetLineWidth(2);
	                    line -> SetLineStyle(7);
	                    line -> Draw("same");
	                }
	    
	                GetEnergyBins(histo, ranges[LRLabel][index], energyBin[LRLabel][index]);
	            }// end MIP (TB)	
	  
	  
	  
	            // -- draw and print energy plots 
	            latex -> Draw("same"); 
	            outFile -> cd();
	            histo->Write();     
	            c -> Print(Form("%s/energy/c_energy__%s.png",plotDir.c_str(),label.c_str()));
	            c -> Print(Form("%s/energy/c_energy__%s.pdf",plotDir.c_str(),label.c_str()));
	            delete c;
	            delete latex;
	  
	        }// end loop over L, R, L-R labels
	        
			ranges_doubleHits[index] = new std::vector<float>;
			// -- loop over PRE ,POST, PREPOST labels
			for(auto PPLabel : PPLabels){
				std::string label1(Form("%s_bar%02dL-R_%s",PPLabel.c_str(),iBar,stepLabel.c_str()));
				std::string label2(Form("%s_ONLY_bar%02dL-R_%s",PPLabel.c_str(),iBar,stepLabel.c_str()));
				
                double index2( (100000000*2)+ 10000000*1 + index ); // index2 per my_cut consistente a quello in 2nd loop 

				//double index2( (100000000*anEvent->nClusters)+(10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );  ex di index2 in 2nd loop

				latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d%s}{V_{OV} = %.2f V, th. = %d DAC}",iBar,PPLabel.c_str(),Vov,int(vth1)));
				latex -> SetNDC();
	            latex -> SetTextFont(42);
	            latex -> SetTextSize(0.04);
	            latex -> SetTextColor(kRed);

				// -- draw energy
	            c = new TCanvas(Form("c_energy_%s",label1.c_str()),Form("c_energy_%s",label1.c_str()));
	            gPad -> SetLogy();
	  
	            histo = (TH1F*)( inFile->Get(Form("h1_energy_%s",label1.c_str())) );      
	            if( !histo ) continue;
	            histo -> SetTitle(";energy [a.u.];entries");
	            histo -> SetLineColor(kRed);
	            histo -> SetLineWidth(2);
	            histo -> Draw();
				if(PPLabel == "POST"){ //-----start fit-----
					float max = histo->GetBinCenter(histo->GetMaximumBin());
					histo->GetXaxis()->SetRangeUser(0, 2024);
					f_landau1[index] = new TF1(Form("f_landau_POST_bar%02dL-R_Vov%.2f_vth1_%02.0f", iBar,Vov,vth1),"[0]*TMath::Landau(x,[1],[2])", 0,2000.);
					float xmin = max * 0.65;
					float xmax = std::min(max*2.5, 1400.);
					f_landau1[index]->SetRange(xmin,xmax);
					//setting dei parametri; 0 ampiezza, 1 mpv, 2 width
					f_landau1[index] -> SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, max, 0.1*max);
					f_landau1[index] -> SetParLimits(1,0,9999);
					f_landau1[index] -> SetParLimits(2,0,9999);
					histo->Fit(f_landau1[index],"QRS");
					if ( f_landau1[index]->GetParameter(1) > 0 ){
						xmin = f_landau1[index]->GetParameter(1) - 2 * std::abs(f_landau1[index]->GetParameter(2));
						if (xmin < minE[std::make_pair(iBar, Vov)]) xmin = minE[std::make_pair(iBar, Vov)] ;
						xmax = std::min(f_landau1[index]->GetParameter(1) * 2.5, 1400.);
						f_landau1[index] -> SetRange(xmin, xmax);
						f_landau1[index] -> SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, f_landau1[index]->GetParameter(1),0.1*f_landau1[index]->GetParameter(1));
					}
					histo->Fit(f_landau1[index],"QRS");
					f_landau1[index] -> SetLineColor(kBlack);
					f_landau1[index] -> SetLineWidth(2);
					f_landau1[index] -> Draw("same");
					//std::cout<<f_landau1[index]->GetParameter(0)<<"  "<<f_landau1[index]->GetParameter(1)<<"  "<<std::abs(f_landau1[index]->GetParameter(2))<<std::endl;
					ranges_doubleHits[index]->push_back( f_landau1[index]->GetParameter(1) - 2.0 * std::abs(f_landau1[index]->GetParameter(2)));
					ranges_doubleHits[index]->push_back(std::min(f_landau1[index]->GetParameter(1)*2.5, 1400.));
					//std::cout<<" "<<std::endl;
					for(auto range: (*ranges_doubleHits[index])){
						TLine* line = new TLine(range,0.,range, histo->GetMaximum());
						line -> SetLineWidth(2);
						line -> SetLineStyle(7);
						line -> Draw("same");
					}
				}//-----end fit-----

				latex -> Draw();
				outFile -> cd();
				histo -> Write();
				c -> Print(Form("%s/energy/c_energy__%s.png",plotDir.c_str(),label1.c_str()));
	            c -> Print(Form("%s/energy/c_energy__%s.pdf",plotDir.c_str(),label1.c_str()));
	            delete c;
	            delete latex;

				latex = new TLatex(0.40,0.85,Form("#splitline{only bar %02d%s}{V_{OV} = %.2f V, th. = %d DAC}",iBar,PPLabel.c_str(),Vov,int(vth1)));
				latex -> SetNDC();
	            latex -> SetTextFont(42);
	            latex -> SetTextSize(0.04);
	            latex -> SetTextColor(kRed);
				
				c = new TCanvas(Form("c_energy_%s",label2.c_str()),Form("c_energy_%s",label2.c_str()));
	            gPad -> SetLogy();
	  
	            histo = (TH1F*)( inFile->Get(Form("h1_energy_%s",label2.c_str())) );      
	            if( !histo ) continue;
	            histo -> SetTitle(";energy [a.u.];entries");
				histo -> GetXaxis()->SetRangeUser(0,1024);
	            histo -> SetLineColor(kRed);
	            histo -> SetLineWidth(2);
	            histo -> Draw();

				latex -> Draw();
				outFile -> cd();
				histo -> Write();
				c -> Print(Form("%s/energy/c_energy__%s.png",plotDir.c_str(),label2.c_str()));
	            c -> Print(Form("%s/energy/c_energy__%s.pdf",plotDir.c_str(),label2.c_str()));
	            delete c;
	            delete latex;

				if(PPLabel == "POST") {
					// --- linear plots for energy of double hits
					c = new TCanvas(Form("c_energy_%s_linear",label1.c_str()),Form("c_energy_%s_linear",label1.c_str()));
					
	                histo = (TH1F*)( inFile->Get(Form("h1_energy_%s",label1.c_str())) );      
	                if( !histo ) continue;
					histo -> GetXaxis()->SetRangeUser(200,1024);
	                histo -> SetTitle(";energy [a.u.];entries");
	                histo -> SetLineColor(kRed);
	                histo -> SetLineWidth(2);
	                histo -> Draw();
					double max = histo->GetBinCenter(histo->GetMaximumBin());
					
					TF1* fitFunc_cut = new TF1(Form("fitFunc_cut_%s",label1.c_str()),"gaus", max-50,max+50);
					histo->Fit(fitFunc_cut, "QRS");
					fitFunc_cut -> SetLineColor(kBlack);
	                fitFunc_cut -> SetLineWidth(2);
	                fitFunc_cut -> Draw("same");
					my_cut[index2] = fitFunc_cut->GetParameter(1) + 2*fitFunc_cut->GetParameter(2); // cut at mean + 2sigma

					//std::cout<<"index: "<<index<<"\t my_cut["<<index2<<"]: "<<my_cut[index2]<<std::endl;

					TLine* line_cut = new TLine(my_cut[index2],0.,my_cut[index2],histo->GetMaximum());
					line_cut -> SetLineColor(kBlack);
					line_cut -> SetLineWidth(2);
					line_cut -> SetLineStyle(3);
					line_cut -> Draw("same");
                    
					latex = new TLatex(0.55,0.75,Form("#splitline{bar %02d%s}{V_{OV} = %.2f V, th. = %d DAC}",iBar,PPLabel.c_str(),Vov,int(vth1)));
				    latex -> SetNDC();
	                latex -> SetTextFont(42);
	                latex -> SetTextSize(0.04);
	                latex -> SetTextColor(kRed);
					latex -> Draw("same");		 

					c -> Print(Form("%s/energy/c_energy__%s_L.png",plotDir.c_str(),label1.c_str()));
					c -> Print(Form("%s/energy/c_energy__%s_L.pdf",plotDir.c_str(),label1.c_str()));
					
					delete c;
					delete latex;
					delete fitFunc_cut;
					delete line_cut;

				}
			    

                if(PPLabel != "PREPOST") {
                    std::string label3(Form("LR_%s_triple_bar%02d_%s",PPLabel.c_str(),iBar,stepLabel.c_str()));
				    latex = new TLatex(0.40,0.85,Form("#splitline{triple hits, %s bar %02d )}{V_{OV} = %.2f V, th. = %d DAC}",PPLabel.c_str(),iBar,Vov,int(vth1)));
					latex -> SetNDC();
	                latex -> SetTextFont(42);
	                latex -> SetTextSize(0.04);
	                latex -> SetTextColor(kRed);

					c = new TCanvas(Form("c_energy_%s",label3.c_str()),Form("c_energy_%s",label3.c_str()));
	                gPad -> SetLogy();

					histo = (TH1F*)( inFile->Get(Form("h1_energy_%s",label3.c_str())) );      
	                if( !histo ) continue;
	                histo -> SetTitle(";energy [a.u.];entries");
				    histo -> GetXaxis()->SetRangeUser(0,1024);
	                histo -> SetLineColor(kRed);
	                histo -> SetLineWidth(2);
	                histo -> Draw();

					latex -> Draw();
				    c -> Print(Form("%s/energy/c_energy__%s.png",plotDir.c_str(),label3.c_str()));
	                c -> Print(Form("%s/energy/c_energy__%s.pdf",plotDir.c_str(),label3.c_str()));
	                delete c;
	                delete latex;

				}

			}//end loop over PP labels

			//histo sum for energy sharing 
			/*
			TH1F* histo1;
			TH1F* histo2;
			TH1F* histo3;
			TH1F* histo4;

			std::string label(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
			std::cout<<label<<std::endl;
			latex = new TLatex(0.40,0.85,Form("#splitline{ bar %02dL-R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
			latex -> SetNDC();
	        latex -> SetTextFont(42);
	        latex -> SetTextSize(0.04);
	        latex -> SetTextColor(kRed);

			histo1 = (TH1F*)(inFile->Get( Form("h1_energy_%s",label.c_str()) ) );
			if(!histo1){std::cout<<Form("h1_energy_%s",label.c_str()) <<" not found"<<std::endl;
			            continue;}
			else{std::cout<<Form("h1_energy_%s",label.c_str()) <<" found"<<std::endl;}
			histo1->GetXaxis()->SetRangeUser(0,2024);
			

			histo2 = (TH1F*)(inFile->Get( Form("h1_energy_PRE_%s",label.c_str()) ) );
			if(!histo2){std::cout<<Form("h1_energy_PRE_%s",label.c_str()) <<" not found"<<std::endl;}
			else {std::cout<<Form("h1_energy_PRE_%s",label.c_str()) <<" found"<<std::endl;}

			histo3 = (TH1F*)(inFile->Get( Form("h1_energy_POST_%s",label.c_str()) ) );
			if(!histo3){std::cout<<Form("h1_energy_POST_%s",label.c_str()) <<" not found"<<std::endl;}
			else{std::cout<<Form("h1_energy_POST_%s",label.c_str()) <<" found"<<std::endl;}

			histo4 = (TH1F*)(inFile->Get( Form("h1_energy_PREPOST_%s",label.c_str()) ) );
			if(!histo4){std::cout<<Form("h1_energy_PREPOST_%s",label.c_str()) <<" not found"<<std::endl;}
			else{std::cout<<Form("h1_energy_PREPOST_%s",label.c_str()) <<" found"<<std::endl;}

			histo1->Add(histo2);
			histo1->Add(histo3);
			histo1->Add(histo4);

			c = new TCanvas(Form("c_energy_sum_bar%02dL-R_%s",iBar,stepLabel.c_str()),Form("c_energy_sum_bar%02dL-R_%s",iBar,stepLabel.c_str()));
			gPad -> SetLogy();
			histo1 -> SetTitle(";energy [a.u.];entries");
	        histo1 -> SetLineColor(kRed);
	        histo1 -> SetLineWidth(2);
	        histo1 -> Draw();
			histo1->GetXaxis()->SetRangeUser(minE[std::make_pair(iBar, Vov)], 1600);
			//-----start fit-----
			float max = histo1->GetBinCenter(histo1->GetMaximumBin());
			histo1->GetXaxis()->SetRangeUser(0, 2024);
            f_landau1[index] = new TF1(Form("f_landau_bar%02dL-R_Vov%.2f_vth1_%02.0f", iBar,Vov,vth1),"[0]*TMath::Landau(x,[1],[2])", 0,2000.);
			float xmin = max * 0.65;
			float xmax = std::min(max*2.5, 1600.);
			f_landau1[index]->SetRange(xmin,xmax);
            //setting dei parametri da chiedere; 0 ampiezza, 1 mpv, 2 width
			f_landau1[index] -> SetParameters(histo1->Integral(histo1->GetMaximumBin(), histo1->GetNbinsX())/10, max, 0.1*max);
	        f_landau1[index] -> SetParLimits(1,0,9999);
	        f_landau1[index] -> SetParLimits(2,0,9999);
            histo1->Fit(f_landau1[index],"QRS");
			if ( f_landau1[index]->GetParameter(1) > 0 ){
				xmin = f_landau1[index]->GetParameter(1) - 2 * std::abs(f_landau1[index]->GetParameter(2));
				if (xmin < minE[std::make_pair(iBar, Vov)]) xmin = minE[std::make_pair(iBar, Vov)] ;
				xmax = std::min(f_landau1[index]->GetParameter(1) * 2.5, 1600.);
				f_landau1[index] -> SetRange(xmin, xmax);
				f_landau1[index] -> SetParameters(histo1->Integral(histo1->GetMaximumBin(), histo1->GetNbinsX())/10, f_landau1[index]->GetParameter(1), 0.1*f_landau1[index]->GetParameter(1));
			}
			histo1->Fit(f_landau1[index],"QRS");
			f_landau1[index] -> SetLineColor(kBlack);
	        f_landau1[index] -> SetLineWidth(2);
	        f_landau1[index] -> Draw("same");
			std::cout<<f_landau1[index]->GetParameter(0)<<"  "<<f_landau1[index]->GetParameter(1)<<"  "<<std::abs(f_landau1[index]->GetParameter(2))<<std::endl;
			my_ranges[index].push_back( f_landau1[index]->GetParameter(1) - 2.0 * std::abs(f_landau1[index]->GetParameter(2)));
			my_ranges[index].push_back(std::min(f_landau1[index]->GetParameter(1)*2.5, 1600.));
            std::cout<<" "<<std::endl;
			for(auto range: (my_ranges[index])){
				TLine* line = new TLine(range,0.,range, histo1->GetMaximum());
				line -> SetLineWidth(2);
				line -> SetLineStyle(7);
				line -> Draw("same");
			}
            //-----end fit-----
			latex -> Draw("same");
			outFile -> cd();
			histo1 -> Write();
			c->Update();
			c -> Print(Form("%s/energy/singleBars/c_energy_sum_%s.png",plotDir.c_str(),label.c_str()));
	        c -> Print(Form("%s/energy/singleBars/c_energy_sum_%s.pdf",plotDir.c_str(),label.c_str()));
	        delete c;
	        delete latex;
			*/
        }// -- end loop over bars

		//draw counter
		std::string label(Form("%s",stepLabel.c_str()));
        c = new TCanvas(Form("c_counter_%s",stepLabel.c_str()),Form("c_counter_%s",stepLabel.c_str()));
		histo = (TH1F*)( inFile->Get(Form("h1_counter_%s",label.c_str())) );      
	    if( !histo ) continue;

		if (histo->Integral() > 0) {
			histo->Scale(1.0 / histo->Integral());
		}
        histo->GetYaxis()->SetRangeUser(0, 1);
	    histo -> SetTitle(";nHits;frequency");
	    histo -> SetLineColor(kBlue);
	    histo -> SetLineWidth(2);
	    histo -> Draw();

		latex = new TLatex(0.40,0.85,Form("#splitline{events distribution}{V_{OV} = %.2f V, th. = %d DAC}",Vov,int(vth1)));
		latex -> SetNDC();
	    latex -> SetTextFont(42);
	    latex -> SetTextSize(0.04);
	    latex -> SetTextColor(kRed);

		latex -> Draw();
		outFile -> cd();
		histo -> Write();
		
		c -> Print(Form("%s/energy/c_counter__%s.png",plotDir.c_str(),label.c_str()));
	    c -> Print(Form("%s/energy/c_counter__%s.pdf",plotDir.c_str(),label.c_str()));
	    delete c;
	    delete latex;
		 
    } // -- end loop over stepLabels
  
    // ---  end 1st plots
  
    
    //------------------------
    //--- 2nd loop over events
    std::map<int,std::map<int,bool> > accept;
    int accepted2=0;
	std::map<double,int> g_index;
    for(auto mapIt : trees){
        ModuleEventClass* anEvent = new ModuleEventClass();
      
        mapIt.second -> SetBranchAddress("event",&anEvent);
      
        int nEntries = mapIt.second->GetEntries();
		//std::cout<<mapIt.first<<std::endl;
        for(int entry = 0; entry < nEntries; ++entry){
	        if( entry%100000 == 0 ) {
	            std::cout << ">>> 2nd loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
	            //TrackProcess(cpu, mem, vsz, rss);
	        }
	  
	        mapIt.second -> GetEntry(entry);
	        
	        bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
	        if (!barFound) continue;
	  
	        int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
			int index3( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + 99 );
	  
	        accept[index1][entry] = false;
			int energyBinAverage=0;

	        if(anEvent->nClusters==1) {
				if(!ranges["L-R"][index1] ) continue;
				energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;
				if(energyBinAverage < 1 ) continue; // selezione eventi singoli
			}
			else if(anEvent->nClusters==2 && anEvent->energyL_post>0 && anEvent->energyR_post>0){
				if(!ranges_doubleHits[index1] ) continue;
				energyBinAverage = FindBin(anEvent->energySum,ranges_doubleHits[index1])+1;
				if(energyBinAverage < 1 ) continue; // selezione eventi doppi "POST"
			}
			else if(anEvent->nClusters==2 && anEvent->energyL_post<0 && anEvent->energyR_post<0){
				energyBinAverage = -1;
				continue; // escludo eventi "PRE"
			}
			else if(anEvent->nClusters==3){
				energyBinAverage=1; 
				if(0.5*(anEvent->energyL + anEvent->energyR) > 940) continue; // seleziono eventi tripli con energia totale < 940
			}
			
            double index2( (100000000*anEvent->nClusters)+(10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
	  
	        accept[index1][entry] = true;
	         
			

	        if( h1_energyRatio[index2] == NULL ){// h1 per energy ratio(L/R),tot(L/R) ,fase , carica e deltaT_raw per una baretta 
	            std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
	      
                h1_energyRatio[index2] = new TH1F(Form("h1_energyRatio_%s",labelLR_energyBin.c_str()),"",1000,0.,5.);
	            h1_totRatio[index2] = new TH1F(Form("h1_totRatio_%s",labelLR_energyBin.c_str()),"",2000,0.,5.);
	            h1_t1fineMean[index2] = new TH1F(Form("h1_t1fineMean_%s",labelLR_energyBin.c_str()),"",1000,0.,1000.);
	            h1_qT1Mean[index2] = new TH1F(Form("h1_qT1Mean_%s",labelLR_energyBin.c_str()),"",250,0.5,1.5);
	            h1_deltaT_raw[index2] = new TH1F(Form("h1_deltaT_raw_%s",labelLR_energyBin.c_str()),"",2000,-24000.,24000.);
	        }

			if(h1_deltaT_bars[index2]==NULL){// deltaT tra tempi medi di due barrette, con media aritmetica e pesata & energy ratio tra due barrette
				std::string labelLR_energyBin(Form("bars%02d-%02d_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->barID+1,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
				h1_deltaT_bars[index2] = new TH1F(Form("h1_deltaT__%s",labelLR_energyBin.c_str()),"",8000,-24000.,24000.);//2000 bins prima
				h1_deltaT_bars_w[index2] = new TH1F(Form("h1_deltaT_w_%s",labelLR_energyBin.c_str()),"",8000,-24000.,24000.);
				h1_deltaT_bars_w_raw[index2] = new TH1F(Form("h1_deltaT_w_raw%s",labelLR_energyBin.c_str()),"",8000,-24000.,24000.);
                h1_energyRatio_bars_F[index2] = new TH1F(Form("h1_energyRatio_F_%s",labelLR_energyBin.c_str()),"",1000,0.,5.);
				h1_energyRatio_bars_B[index2] = new TH1F(Form("h1_energyRatio_B_%s",labelLR_energyBin.c_str()),"",1000,0.,5.);
	            h2_deltaT_energyRatio_bars[index2] = new TH2F(Form("h2_detaT_energyRatio_%s",labelLR_energyBin.c_str()),"",2000,0.,5.,2000,-24000.,24000.);
				h2_deltaT_energyRatio_bars_w[index2] = new TH2F(Form("h2_detaT_energyRatio_w_%s",labelLR_energyBin.c_str()),"",2000,0.,5.,2000,-24000.,24000.); //2000 bins prima su asse y
	            p1_deltaT_energyRatio_bars[index2] = new TProfile(Form("p1_deltaT_energyRatio_%s",labelLR_energyBin.c_str()),"",80,0.,5.);
				p1_deltaT_energyRatio_bars_w[index2] = new TProfile(Form("p1_deltaT_energyRatio_w_%s",labelLR_energyBin.c_str()),"",80,0.,5.);
			}

			if(h2_energy_correlation[index2]==NULL){// correlazione di energia tra le due barrette da cui calcolo deltaT
				std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
				
				h2_energy_correlation[index2]= new TH2F(Form("h2_energy_correlation_%s",labelLR_energyBin.c_str()),"",500,0,1000,500,0,1000);
				p1_energy_correlation[index2]= new TProfile(Form("p1_energy_correleation_%s",labelLR_energyBin.c_str()),"",100,0,1000);
			}

			if(p1_tL_eL[index2]==NULL){
				std::string labelL(Form("bar%02dL_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
				std::string labelR(Form("bar%02dR_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
				std::string label(Form("bars%02d-%02d_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->barID+1,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
				// time L vs energy L & time R vs energy R di singola barretta
				p1_tL_eL[index2] = new TProfile(Form("p1_tL_eL%s",labelL.c_str()),"",100,0,1000);
				p1_tR_eR[index2] = new TProfile(Form("p1_tR_eR%s",labelR.c_str()),"",100,0,1000);
                // energyL media vs energyR media di due barrette da cui clacolo deltaT
				p1_enL_enR_ave[index2] = new TProfile(Form("p1_enL_enR_ave_%s",label.c_str()),"",100,0,1000);
				h2_enL_enR_ave[index2] = new TH2F(Form("h2_enL_enR_ave_%s",label.c_str()),"",512,0,1024,512,0,1024);
			}

			if(h1_energyRatio_REF[index3]==NULL){
				h1_energyRatio_REF[index3] = new TH1F(Form("h1_energyRatio_externalBar_Vov%.2f_th%02d",anEvent->Vov,anEvent->vth1),"",1000,0.,5.);
				h1_t1fineMean_REF[index3] = new TH1F(Form("h1_t1fineMean_externalBar_Vov%.2f_th%02d",anEvent->Vov,anEvent->vth1),"",1000,0.,1000.);
			}
            // reference module 
			if (fabs(anEvent->timeR_ext-anEvent->timeL_ext)<10000 && anEvent->energyR_ext>0 && anEvent->energyL_ext>0){
				h1_energyRatio_REF[index3] -> Fill( anEvent->energyR_ext / anEvent->energyL_ext );
				h1_t1fineMean_REF[index3] -> Fill(0.5*(anEvent->t1fineR_ext + anEvent->t1fineL_ext));
			}

			
	        if (fabs(anEvent->timeR-anEvent->timeL)<10000){
				accepted2++;
	            //if ((anEvent->energyR / anEvent->energyL >0) & (anEvent->energyR / anEvent->energyL <5)){
	            if ( ((anEvent->energyR / anEvent->energyL > -999) & (anEvent->energyR / anEvent->energyL <9999999))){
		            h1_energyRatio[index2] -> Fill( anEvent->energyR / anEvent->energyL );						     
		            h1_totRatio[index2] -> Fill( anEvent->totR / anEvent->totL );
		            h1_deltaT_raw[index2] -> Fill( anEvent->timeR-anEvent->timeL );

		            h1_t1fineMean[index2] -> Fill( 0.5 * (anEvent->t1fineR + anEvent->t1fineL) );
		            h1_qT1Mean[index2] -> Fill( 0.5 * (anEvent->qT1R + anEvent->qT1L) );
	            }
	        }

			if(anEvent->nClusters==2){
				p1_tL_eL[index2] -> Fill(anEvent->energyL,anEvent->timeL);
				p1_tR_eR[index2] -> Fill(anEvent->energyR,anEvent->timeR);
			}

           	// --- double hits (post) 
			if(anEvent->nClusters==2 && anEvent->energyL_post>0 && anEvent->energyR_post>0){
				if(anEvent->barID==15) continue;

                // --- variables definition
				double E_mean=0.5*(anEvent->energyL+anEvent->energyR);
				double E_mean_post=0.5*(anEvent->energyL_post+anEvent->energyR_post);
				float enL_mean = 0.5*(anEvent->energyL+anEvent->energyL_post);
				float enR_mean = 0.5*(anEvent->energyR+anEvent->energyR_post);
				  // -- arithmetic mean 
				double t_mean=0.5*(anEvent->timeL+anEvent->timeR);
				double t_mean_post=0.5*(anEvent->timeL_post+anEvent->timeR_post);
				  // -- weighted mean
				double t_mean_w=((static_cast<double>(anEvent->energyL)*static_cast<double>(anEvent->timeL))+(static_cast<double>(anEvent->energyR)*static_cast<double>(anEvent->timeR)))/(static_cast<double>(anEvent->energyL)+static_cast<double>(anEvent->energyR));
				double t_mean_post_w=((static_cast<double>(anEvent->energyL_post)*static_cast<double>(anEvent->timeL_post))+(static_cast<double>(anEvent->energyR_post)*static_cast<double>(anEvent->timeR_post)))/(static_cast<double>(anEvent->energyL_post)+static_cast<double>(anEvent->energyR_post));

				if(fabs(t_mean-t_mean_post)<10000){
					// --- energy correlation between bar(i) and bar(i+1)
					h2_energy_correlation[index2]->Fill(E_mean_post,E_mean);
					p1_energy_correlation[index2]->Fill(E_mean_post,E_mean);
					// --- energy left average vs energy right average
					h2_enL_enR_ave[index2]->Fill(enL_mean,enR_mean);
					p1_enL_enR_ave[index2]->Fill(enL_mean,enR_mean);
					if((E_mean+E_mean_post)<=my_cut[index2]){
						// --- energy correlation between bar(i) and bar(i+1) and vice versa
						h1_energyRatio_bars_F[index2]->Fill(E_mean/E_mean_post);
					    h1_energyRatio_bars_B[index2]->Fill(E_mean_post/E_mean);
                        // --- deltaT vs energyRatio for bar(i) and bar(i+1)
						h2_deltaT_energyRatio_bars[index2]->Fill(E_mean/E_mean_post,t_mean-t_mean_post);
					    p1_deltaT_energyRatio_bars[index2]->Fill(E_mean/E_mean_post,t_mean-t_mean_post);
					}	
				}
				if(fabs(t_mean_w-t_mean_post_w)<10000){
					h1_deltaT_bars_w_raw[index2] -> Fill(t_mean_w-t_mean_post_w);
					/*if((E_mean+E_mean_post)<=my_cut[index2]){
						// --- deltaT vs energyRatio for bar(i) and bar(i+1) WEIGHTED MEAN
						h2_deltaT_energyRatio_bars_w[index2]->Fill(E_mean/E_mean_post,t_mean_w-t_mean_post_w);
					    p1_deltaT_energyRatio_bars_w[index2]->Fill(E_mean/E_mean_post,t_mean_w-t_mean_post_w);
					}*/
				}
                																 
			} // --- double hits (post)
 
	    } // end loop over entries
    }

    int totEntries=0;
	std::cout<<"\n"<<std::endl;
	for(auto mapIt: trees){
		//std::cout<<"data tree: "<<mapIt.first<<" entries: "<<mapIt.second->GetEntries()<<std::endl;
		totEntries+=mapIt.second->GetEntries();
	}
	std::cout<<"\ntotal events at start for new selection: "<<totEntries<<std::endl;
      
  
    //------------------
    //--- draw 2nd plots
	std::cout<<"2nd loop accepted events: "<<accepted2<<std::endl;
    std::map<double,float> CTRMeans;
    std::map<double,float> CTRSigmas;
	std::map<double,float> my_CTRMeans;
    std::map<double,float> my_CTRSigmas;
  
    std::map<double,TF1*> fitFunc_energyRatio;
    std::map<double,TF1*> fitFunc_totRatio;

	//std::map<double,TF1*> fitFunc_energyRatio_bars;
    std::map<double, TF1*> fitFunc_timeL;
	std::map<double, TF1*> fitFunc_timeR;
	std::map<double, TF1*> fitFunc_timeL_2;
	std::map<double, TF1*> fitFunc_timeR_2;
	std::map<double,TF1*> fitFunc_energyRatioCorr_bars;
	std::map<double,TF1*> fitFunc_energyRatioCorr_bars_w;
  
    for(auto mapIt : h1_deltaT_raw){
        double index = mapIt.first;
      
        FindSmallestInterval(vals,h1_deltaT_raw[index],0.68);
        float mean = vals[0];
        float min = vals[4];
        float max = vals[5];
        float delta = max-min;
        float sigma = 0.5*delta;
        float effSigma = sigma;
        CTRMeans[index] = mean;
        CTRSigmas[index] = effSigma;
    }
    
	// ho provato questa selezione solo per deltaT tra barre con t_ave(i) calcolato con media pesata, peggiora la situazione e non l'ho implementatato per gli altri metodi 
	for(auto mapIt : h1_deltaT_bars_w_raw){
        double index = mapIt.first;
      
        FindSmallestInterval(my_vals,h1_deltaT_bars_w_raw[index],0.68);
        float mean = my_vals[0];
        float min = my_vals[4];
        float max = my_vals[5];
        float delta = max-min;
        float sigma = 0.5*delta;
        float effSigma = sigma;
        my_CTRMeans[index] = mean;
        my_CTRSigmas[index] = effSigma;
    }

    std::cout<<"Drawing 2nd plots ... "<<std::endl;
    for(auto stepLabel : stepLabels){
        float Vov = map_Vovs[stepLabel];
        float vth1 = map_ths[stepLabel];

        std::string extLabel(Form("externalBar_L-R_%s",stepLabel.c_str())); 
		int index3( (10000*int(Vov*100.)) + (100*vth1) + 99 );

		// --- draw energy ratio for external bar
		c = new TCanvas(Form("c_energyRatio_%s",extLabel.c_str()),Form("c_energyRatio_%s",extLabel.c_str()));
		histo = h1_energyRatio_REF[index3];
		histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
		histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
		histo -> SetTitle(Form(";energy_{right} / energy_{left};entries"));
		histo -> SetLineColor(kRed);
		histo -> SetLineWidth(2);
		histo -> Draw();
		histo -> Write();

		fitFunc_energyRatio_REF[index3] = new TF1(Form("fitFunc_energyRatio_%s",extLabel.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
		histo -> Fit(fitFunc_energyRatio_REF[index3],"QNRS");
		histo -> Fit(fitFunc_energyRatio_REF[index3],"QSR+","",fitFunc_energyRatio_REF[index3]->GetParameter(1)-2.*fitFunc_energyRatio_REF[index3]->GetParameter(2),fitFunc_energyRatio_REF[index3]->GetParameter(1)+2.*fitFunc_energyRatio_REF[index3]->GetParameter(2));
		histo -> Fit(fitFunc_energyRatio_REF[index3],"QSR+","",fitFunc_energyRatio_REF[index3]->GetParameter(1)-2.*fitFunc_energyRatio_REF[index3]->GetParameter(2),fitFunc_energyRatio_REF[index3]->GetParameter(1)+2.*fitFunc_energyRatio_REF[index3]->GetParameter(2));
		  
		fitFunc_energyRatio_REF[index3] -> SetLineColor(kBlack);
		fitFunc_energyRatio_REF[index3] -> SetLineWidth(2);
		fitFunc_energyRatio_REF[index3] -> Draw("same");

		//FIXME
		//fitFunc_energyRatio[index2] -> SetParameter(1,histo->GetMean());
		//fitFunc_energyRatio[index2] -> SetParameter(2,histo->GetRMS());
			  
		latex = new TLatex(0.40,0.85,Form("#splitline{external bar}{V_{OV} = %.2f V, th. = %d DAC}",Vov,int(vth1)));
		latex -> SetNDC();
		latex -> SetTextFont(42);
		latex -> SetTextSize(0.04);
		latex -> SetTextColor(kRed);
		latex -> Draw("same");
			  
		c -> Print(Form("%s/externalBar/c_energyRatio__%s.png",plotDir.c_str(),extLabel.c_str()));
		c -> Print(Form("%s/externalBar/c_energyRatio__%s.pdf",plotDir.c_str(),extLabel.c_str()));
		delete c;

		// -- draw t1fine (average between left and right) for external bar
		c = new TCanvas(Form("c_t1fineMean_%s",extLabel.c_str()),Form("c_t1fineMean_%s",extLabel.c_str()));
		histo = h1_t1fineMean_REF[index3];
		histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
		histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
		histo -> SetTitle(Form(";(t1fine_{right}+t1fine_{left})/2;entries"));
		histo -> SetLineColor(kRed);
		histo -> SetLineWidth(2);
		histo -> Draw();
		histo -> Write();
			  
		latex -> Draw("same");
			  
		c -> Print(Form("%s/externalBar/c_t1fineMean__%s.png",plotDir.c_str(),extLabel.c_str()));
		c -> Print(Form("%s/externalBar/c_t1fineMean__%s.pdf",plotDir.c_str(),extLabel.c_str()));
		delete latex;
		delete c;

        for(int iBar = 0; iBar < 16; ++iBar) {
	        bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
            if (!barFound) continue;      
	  
	        std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));;          
			std::string labelBars(Form("bars%02d-%02d_%s",iBar,iBar+1,stepLabel.c_str()));
	  
	        int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );
	  
	        for(int i=1; i<=3;i++){
				int nEnergyBins;
				if(i==1){
				    if( !ranges["L-R"][index1]) continue;
					nEnergyBins = ranges["L-R"][index1]->size()-1;
				}
				else if(i==2){
				    if( !ranges_doubleHits[index1]) continue;
					nEnergyBins = ranges_doubleHits[index1]->size()-1;
				}
				else{ nEnergyBins = 1;}

				for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin){
					//if (ranges["L-R"][index1]->at(iEnergyBin)<0) continue;
					double index2( 100000000*i+(10000000*iEnergyBin+index1) );
			  
					if (!h1_energyRatio[index2]) continue;
					std::string labelLR_energyBin;
					std::string labelBars_energyBin;
					if(i==1){labelLR_energyBin = Form("%s_energyBin%02d_S",labelLR.c_str(),iEnergyBin);}
					else if(i==2){
						labelLR_energyBin = Form("%s_energyBin%02d_D",labelLR.c_str(),iEnergyBin);
						labelBars_energyBin = Form("%s_energyBin%02d_D",labelBars.c_str(),iEnergyBin);
					}
					else if(i==3){labelLR_energyBin = Form("%s_energyBin%02d_T",labelLR.c_str(),iEnergyBin);}
			  
			  
					// -- draw energy ratio 
					c = new TCanvas(Form("c_energyRatio_%s",labelLR_energyBin.c_str()),Form("c_energyRatio_%s",labelLR_energyBin.c_str()));
					histo = h1_energyRatio[index2];
					histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
					histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
					histo -> SetTitle(Form(";energy_{right} / energy_{left};entries"));
					histo -> SetLineColor(kRed);
					histo -> SetLineWidth(2);
					histo -> Draw();
					histo -> Write();
			        
					energyRatio_ranges[index2] = new std::vector<double>;
					fitFunc_energyRatio[index2] = new TF1(Form("fitFunc_energyRatio_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
					histo -> Fit(fitFunc_energyRatio[index2],"QNRS");
					histo -> Fit(fitFunc_energyRatio[index2],"QSR+","",fitFunc_energyRatio[index2]->GetParameter(1)-2.*fitFunc_energyRatio[index2]->GetParameter(2),fitFunc_energyRatio[index2]->GetParameter(1)+2.*fitFunc_energyRatio[index2]->GetParameter(2));
					histo -> Fit(fitFunc_energyRatio[index2],"QSR+","",fitFunc_energyRatio[index2]->GetParameter(1)-2.*fitFunc_energyRatio[index2]->GetParameter(2),fitFunc_energyRatio[index2]->GetParameter(1)+2.*fitFunc_energyRatio[index2]->GetParameter(2));
			  
					fitFunc_energyRatio[index2] -> SetLineColor(kBlack);
					fitFunc_energyRatio[index2] -> SetLineWidth(2);
					fitFunc_energyRatio[index2] -> Draw("same");
					energyRatio_ranges[index2]->push_back(fitFunc_energyRatio[index2]->GetParameter(1) - std::abs(fitFunc_energyRatio[index2]->GetParameter(2)));
					energyRatio_ranges[index2]->push_back(fitFunc_energyRatio[index2]->GetParameter(1) + std::abs(fitFunc_energyRatio[index2]->GetParameter(2)));
					if(i==2){
						for(auto range: (*energyRatio_ranges[index2])){
	                    	TLine* line = new TLine(range,0.,range, histo->GetMaximum());
	                    	line -> SetLineWidth(2);
	                    	line -> SetLineStyle(7);
	                    	line -> Draw("same");
	                	}
				    }

			  
					//FIXME
					fitFunc_energyRatio[index2] -> SetParameter(1,histo->GetMean());
					fitFunc_energyRatio[index2] -> SetParameter(2,histo->GetRMS());
			  
					latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					latex -> SetNDC();
					latex -> SetTextFont(42);
					latex -> SetTextSize(0.04);
					latex -> SetTextColor(kRed);
					latex -> Draw("same");
			  
					c -> Print(Form("%s/energyRatio/c_energyRatio__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/energyRatio/c_energyRatio__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete latex;
					delete c;
			  
			  
					// -- draw tot ratio 
					c = new TCanvas(Form("c_totRatio_%s",labelLR_energyBin.c_str()),Form("c_totRatio_%s",labelLR_energyBin.c_str()));
					histo = h1_totRatio[index2];
					histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
					histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
					histo -> SetTitle(Form(";tot_{right} / tot_{left};entries"));
					histo -> SetLineColor(kRed);
					histo -> SetLineWidth(2);
					histo -> Draw();
					histo -> Write();
			  
					fitFunc_totRatio[index2] = new TF1(Form("fitFunc_totRatio_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
					histo -> Fit(fitFunc_totRatio[index2],"QNRS");
					histo -> Fit(fitFunc_totRatio[index2],"QSR+","",fitFunc_totRatio[index2]->GetParameter(1)-2.*fitFunc_totRatio[index2]->GetParameter(2),fitFunc_totRatio[index2]->GetParameter(1)+2.*fitFunc_totRatio[index2]->GetParameter(2));
					histo -> Fit(fitFunc_totRatio[index2],"QSR+","",fitFunc_totRatio[index2]->GetParameter(1)-2.*fitFunc_totRatio[index2]->GetParameter(2),fitFunc_totRatio[index2]->GetParameter(1)+2.*fitFunc_totRatio[index2]->GetParameter(2));
			  
					fitFunc_totRatio[index2] -> SetLineColor(kBlack);
					fitFunc_totRatio[index2] -> SetLineWidth(2);
					fitFunc_totRatio[index2] -> Draw("same");
			  
					latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					latex -> SetNDC();
					latex -> SetTextFont(42);
					latex -> SetTextSize(0.04);
					latex -> SetTextColor(kRed);
					latex -> Draw("same");
			  
					c -> Print(Form("%s/totRatio/c_totRatio__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/totRatio/c_totRatio__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete latex;
					delete c;
			  
			  
					// -- draw t1fine (average between left and right)
					c = new TCanvas(Form("c_t1fineMean_%s",labelLR_energyBin.c_str()),Form("c_t1fineMean_%s",labelLR_energyBin.c_str()));
					histo = h1_t1fineMean[index2];
					histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
					histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
					histo -> SetTitle(Form(";(t1fine_{right}+t1fine_{left})/2;entries"));
					histo -> SetLineColor(kRed);
					histo -> SetLineWidth(2);
					histo -> Draw();
					histo -> Write();
			  
					latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					latex -> SetNDC();
					latex -> SetTextFont(42);
					latex -> SetTextSize(0.04);
					latex -> SetTextColor(kRed);
					latex -> Draw("same");
			  
					c -> Print(Form("%s/t1fine/c_t1fineMean__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/t1fine/c_t1fineMean__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete latex;
					delete c;
			  
			  
			  
					delete h1_energyRatio[index2];
					delete h1_totRatio[index2];
			  
					// -- draw qT1 (average between left and right)
					c = new TCanvas(Form("c_qT1Mean_%s",labelLR_energyBin.c_str()),Form("c_qT1Mean_%s",labelLR_energyBin.c_str()));
					histo = h1_qT1Mean[index2];
					histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
					histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
					histo -> SetTitle(Form(";(qT1_{right}+qT1_{left})/2;entries"));
					histo -> SetLineColor(kRed);
					histo -> SetLineWidth(2);
					histo -> Draw();
					histo -> Write();
			  
					latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					latex -> SetNDC();
					latex -> SetTextFont(42);
					latex -> SetTextSize(0.04);
					latex -> SetTextColor(kRed);
					latex -> Draw("same");
			  
					c -> Print(Form("%s/qT1/c_qT1Mean__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/qT1/c_qT1Mean__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete c;
					delete latex;

                    // ---  energy ratio + energy scatter plots
                    if(i==2){
						if(iBar==15) continue;
						/*
						// --- draw tL vs eL
						c = new TCanvas(Form("c_tL_vs_eL_bar%02dL_%s_D",iBar,stepLabel.c_str()),Form("c_tL_vs_eL_bar%02dL_%s_D",iBar,stepLabel.c_str()));
						prof = p1_tL_eL[index2];
						prof -> SetTitle(Form(";energy_{L} [a.u.];time_{L} [ps]"));
						prof -> SetMarkerSize(0.4);
						prof -> Draw("");

						latex = new TLatex(0.40,0.5,Form("#splitline{bar %02dL}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kRed);
					    latex -> Draw("same");

						c->Print(Form("%s/energyCorrelation/c_tL_vs_eL_bar%02dL_%s_D.png",plotDir.c_str(),iBar,stepLabel.c_str()));
						c->Print(Form("%s/energyCorrelation/c_tL_vs_eL_bar%02dL_%s_D.pdf",plotDir.c_str(),iBar,stepLabel.c_str()));
						delete latex;
						delete c;
                        
						// --- draw tR vs eR
						c = new TCanvas(Form("c_tR_vs_eR_bar%02dR_%s_D",iBar,stepLabel.c_str()),Form("c_tR_vs_eR_bar%02dR_%s_D",iBar,stepLabel.c_str()));
						prof = p1_tL_eL[index2];
						prof -> SetTitle(Form(";energy_{R} [a.u.];time_{R} [ps]"));
						prof -> SetMarkerSize(0.4);
						prof -> Draw("");

						latex = new TLatex(0.40,0.85,Form("#splitline{bar %02dR}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kRed);
					    latex -> Draw("same");

						c->Print(Form("%s/energyCorrelation/c_tR_vs_eR_bar%02dR_%s_D.png",plotDir.c_str(),iBar,stepLabel.c_str()));
						c->Print(Form("%s/energyCorrelation/c_tR_vs_eR_bar%02dR_%s_D.pdf",plotDir.c_str(),iBar,stepLabel.c_str()));
						delete latex;
						delete c;
						*/

						// --- draw energyL-ave vs energyR-ave
                        c = new TCanvas(Form("c_enL_vs_enR_ave_%s",labelBars_energyBin.c_str()),Form("c_enL_vs_enR_ave_%s",labelBars_energyBin.c_str()));
						if(!p1_enL_enR_ave[index2]) continue;
					    c -> SetGridy();
			  
					    h2 = h2_enL_enR_ave[index2]; 
					    h2 -> SetTitle(Form(";energyR [a.u.];energyL [a.u.]"));
					    h2 -> Draw("colz");
					    prof = p1_enL_enR_ave[index2];
						prof -> SetMarkerSize(0.5);
					    prof -> Draw("psame");
                        
						
					    latex = new TLatex(0.40,0.85,Form("#splitline{bars %02d-%02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,iBar+1,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kRed);
					    latex -> Draw("same");
			  
					    c -> Print(Form("%s/energy/c_enL_vs_enR_ave_%s.png",plotDir.c_str(),labelBars_energyBin.c_str()));
						c -> Print(Form("%s/energy/c_enL_vs_enR_ave_%s.pdf",plotDir.c_str(),labelBars_energyBin.c_str()));
					    delete c;
					    delete latex;
						// -----------------------------------------------------------------------------------

						// -- draw energy ratio for bars : Forward(bar/post_bar)
						c = new TCanvas(Form("c_energyRatio_F_%s",labelBars_energyBin.c_str()),Form("c_energyRatio_F_%s",labelBars_energyBin.c_str()));
						histo = h1_energyRatio_bars_F[index2];
					    histo -> SetTitle(Form(";energy_{%02d} / energy_{%02d};entries",iBar,iBar+1));
					    histo -> SetLineColor(kRed);
					    histo -> SetLineWidth(2);
					    histo -> Draw();
					    histo -> GetXaxis() -> SetRangeUser(0,4.5);
						

						/*fitFunc_energyRatio_bars[index2] = new TF1(Form("fitFunc_energyRatio_%s",labelBars_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
					    
						histo -> Fit(fitFunc_energyRatio_bars[index2],"QNRS");
					    histo -> Fit(fitFunc_energyRatio_bars[index2],"QSR+","",fitFunc_energyRatio_bars[index2]->GetParameter(1)-2.*fitFunc_energyRatio_bars[index2]->GetParameter(2),fitFunc_energyRatio_bars[index2]->GetParameter(1)+2.*fitFunc_energyRatio_bars[index2]->GetParameter(2));
					    //histo -> Fit(fitFunc_energyRatio[index2],"QSR+","",fitFunc_energyRatio[index2]->GetParameter(1)-2.*fitFunc_energyRatio[index2]->GetParameter(2),fitFunc_energyRatio[index2]->GetParameter(1)+2.*fitFunc_energyRatio[index2]->GetParameter(2));
			  
					    fitFunc_energyRatio_bars[index2] -> SetLineColor(kBlack);
					    fitFunc_energyRatio_bars[index2] -> SetLineWidth(2);
					    fitFunc_energyRatio_bars[index2] -> Draw("same");
			  
					    fitFunc_energyRatio_bars[index2] -> SetParameter(1,histo->GetMean());
					    fitFunc_energyRatio_bars[index2] -> SetParameter(2,histo->GetRMS());
			            */
					    latex = new TLatex(0.40,0.85,Form("#splitline{bars %02d-%02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,iBar+1,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.03);
					    latex -> SetTextColor(kRed);
					    latex -> Draw("same");
                        double maxVal = histo->GetBinCenter(histo->GetMaximumBin());
						latex = new TLatex(0.40,0.78,Form("max = %.2f ",maxVal));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kBlue);
					    latex -> Draw("same");
                        TLine* line = new TLine(maxVal,0,maxVal,histo ->GetBinContent(histo -> GetMaximumBin()));
						line -> SetLineColor(kBlue);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						line -> Draw("same");
			            histo -> Write();
					    c -> Print(Form("%s/energyRatio/c_energyRatio__F_%s.png",plotDir.c_str(),labelBars_energyBin.c_str()));
					    c -> Print(Form("%s/energyRatio/c_energyRatio__F_%s.pdf",plotDir.c_str(),labelBars_energyBin.c_str()));
						delete line;
					    delete latex;
					    delete c;

                        /*
						// -- draw energy ratio for bars : Backward(post_bar/bar)
                        c = new TCanvas(Form("c_energyRatio_B_%s",labelBars_energyBin.c_str()),Form("c_energyRatio_B_%s",labelBars_energyBin.c_str()));
						histo = h1_energyRatio_bars_B[index2];
					    histo -> SetTitle(Form(";energy_{%02d} / energy_{%02d};entries",iBar+1,iBar));
					    histo -> SetLineColor(kRed);
					    histo -> SetLineWidth(2);
					    // histo -> Draw();
					    histo -> GetXaxis() -> SetRangeUser(0,4.5);

						latex = new TLatex(0.40,0.85,Form("#splitline{bars %02d-%02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,iBar+1,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.03);
					    latex -> SetTextColor(kRed);
					   // latex -> Draw("same");
                        double maxVal_B = histo->GetBinCenter(histo->GetMaximumBin());
						latex = new TLatex(0.40,0.78,Form("max = %.2f ",maxVal_B));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kBlue);
					    //latex -> Draw("same");
                        line = new TLine(maxVal,0,maxVal_B,histo ->GetBinContent(histo -> GetMaximumBin()));
						line -> SetLineColor(kBlue);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						//line -> Draw("same");
			            histo -> Write();
					    //c -> Print(Form("%s/energyRatio/c_energyRatio__B_%s.png",plotDir.c_str(),labelBars_energyBin.c_str()));
					    //c -> Print(Form("%s/energyRatio/c_energyRatio__B_%s.pdf",plotDir.c_str(),labelBars_energyBin.c_str()));
						delete line;
					    delete latex;
					    delete c;
						*/

						c = new TCanvas(Form("c_energyRatio_FBratio_%s",labelBars_energyBin.c_str()),Form("c_energyRatio_FBratio_%s",labelBars_energyBin.c_str()));
						TH1F* h1_energyRatio_FBratio = (TH1F*)h1_energyRatio_bars_F[index2]->Clone(Form("h1_energyRatio_FBratio_%s",labelBars_energyBin.c_str()));
						h1_energyRatio_FBratio -> Divide(h1_energyRatio_bars_B[index2]);
						h1_energyRatio_FBratio -> SetTitle(Form(";energy_{%02d} / energy_{%02d};ratio",iBar,iBar+1));
						h1_energyRatio_FBratio -> SetLineColor(kRed);
					    h1_energyRatio_FBratio -> SetLineWidth(2); 
					    h1_energyRatio_FBratio -> Draw();
						latex = new TLatex(0.40,0.85,Form("#splitline{bars %02d-%02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,iBar+1,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.03);
					    latex -> SetTextColor(kRed);
						latex -> Draw("same");
						
						c -> Print(Form("%s/energyRatio/c_energyRatio__FBratio_%s.png",plotDir.c_str(),labelBars_energyBin.c_str()));
					    c -> Print(Form("%s/energyRatio/c_energyRatio__FBratio_%s.pdf",plotDir.c_str(),labelBars_energyBin.c_str()));
						delete c;


						// --draw deltaT_bars vs energy ratio
						c = new TCanvas(Form("c_deltaT_energyRatio_%s",labelBars_energyBin.c_str()),Form("c_deltaT_energyRatio_%s",labelBars_energyBin.c_str()));
						h2 = h2_deltaT_energyRatio_bars[index2];
						h2 -> GetYaxis()->SetRangeUser(-5000, +5000);
					    h2 -> SetTitle(Form(";energy_{%02d} / energy_{%02d};#DeltaT_{%02d-%02d} [ps]",iBar,iBar+1,iBar,iBar+1));
					    h2 -> Draw("colz");
					    outFile -> cd();
					    h2 -> Write();

						prof = p1_deltaT_energyRatio_bars[index2];
						prof -> SetMarkerSize(0.4);
						prof -> Draw("psame");

						fitFunc_energyRatioCorr_bars[index2] = new TF1(Form("fitFunc_energyRatioCorr_%s",labelLR_energyBin.c_str()),"[0]/x + [1]/(x*x) + [2]/(x*x*x)- [3]*x - [4]*x*x -[5]*x*x*x",0.2,5);
						//fitFunc_energyRatioCorr_bars[index2] = new TF1(Form("fitFunc_energyRatioCorr_%s",labelLR_energyBin.c_str()),"[0]/x - [1]*x - [2]*x*x -[3]*x*x*x",0.2,5);
					    prof -> Fit(fitFunc_energyRatioCorr_bars[index2],"QRS");
					    fitFunc_energyRatioCorr_bars[index2] -> SetLineColor(kRed);
					    fitFunc_energyRatioCorr_bars[index2] -> SetLineWidth(2);
					    fitFunc_energyRatioCorr_bars[index2] -> Draw("same");

			            gPad->Update();

			            line = new TLine(0.2,-5000,0.2,5000);
						line -> SetLineColor(kBlack);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						line -> Draw("same");
						
						line = new TLine(0,0,5,0);
						line -> SetLineColor(kRed);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						line -> Draw("same");
						
						line = new TLine(1,-5000,1,5000);
						line -> SetLineColor(kRed);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						line -> Draw("same");

					    latex = new TLatex(0.40,0.85,Form("#splitline{bars %02d-%02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,iBar+1,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kRed);
					    latex -> Draw("same");
			  
					    c -> Print(Form("%s/energyRatioCorr/c_deltaT_energyRatio__%s.png",plotDir.c_str(),labelBars_energyBin.c_str()));
					    c -> Print(Form("%s/energyRatioCorr/c_deltaT_energyRatio__%s.pdf",plotDir.c_str(),labelBars_energyBin.c_str())); 
                        delete line;
						delete latex;
					    delete c;
                        
						c = new TCanvas(Form("c_deltaT_energyRatio_prof_%s",labelBars_energyBin.c_str()),Form("c_deltaT_energyRatio_prof_%s",labelBars_energyBin.c_str()));
						prof -> GetYaxis()->SetRangeUser(-5000, +5000);
					    prof -> SetTitle(Form(";energy_{%02d} / energy_{%02d};#DeltaT_{%02d-%02d} [ps]",iBar,iBar+1,iBar,iBar+1));
						prof -> SetMarkerSize(0.4);
						prof -> Draw("psame");

						gStyle->SetOptFit(1111);
                        prof->SetStats(kTRUE); 

					    fitFunc_energyRatioCorr_bars[index2] -> Draw("same");

			            line = new TLine(0.2,-5000,0.2,5000);
						line -> SetLineColor(kBlack);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						line -> Draw("same");

					    latex = new TLatex(0.40,0.85,Form("#splitline{bars %02d-%02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,iBar+1,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kRed);
					    latex -> Draw("same");
			  
					    c -> Print(Form("%s/energyRatioCorr/c_deltaT_energyRatio_prof_%s.png",plotDir.c_str(),labelBars_energyBin.c_str()));
					    c -> Print(Form("%s/energyRatioCorr/c_deltaT_energyRatio_prof_%s.pdf",plotDir.c_str(),labelBars_energyBin.c_str())); 
                       
						delete latex;
					    delete c;

						// --draw deltaT_bars vs energy ratio for weigthed mean
						/*
						c = new TCanvas(Form("c_deltaT_energyRatio_w_%s",labelBars_energyBin.c_str()),Form("c_deltaT_energyRatio_w_%s",labelBars_energyBin.c_str()));
						h2 = h2_deltaT_energyRatio_bars_w[index2];
						h2 -> GetYaxis()->SetRangeUser(-5000, +5000);
					    h2 -> SetTitle(Form(";energy_{%02d} / energy_{%02d};#DeltaT_{%02d-%02d} [ps]",iBar,iBar+1,iBar,iBar+1));
					    h2 -> Draw("colz");
					    outFile -> cd();
					    h2 -> Write();

						prof = p1_deltaT_energyRatio_bars_w[index2];
						prof -> SetMarkerSize(0.4);
						prof -> Draw("psame");

						fitFunc_energyRatioCorr_bars_w[index2] = new TF1(Form("fitFunc_energyRatioCorr_%s",labelLR_energyBin.c_str()),"[0]/x + [1]/(x*x) + [2]/(x*x*x)- [3]*x - [4]*x*x -[5]*x*x*x",0.2,5);
						//fitFunc_energyRatioCorr_bars_w[index2] = new TF1(Form("fitFunc_energyRatioCorr_%s",labelLR_energyBin.c_str()),"[0]/x - [1]*x - [2]*x*x -[3]*x*x*x",0.2,5);
					    prof -> Fit(fitFunc_energyRatioCorr_bars_w[index2],"QRS");
					    fitFunc_energyRatioCorr_bars_w[index2] -> SetLineColor(kRed);
					    fitFunc_energyRatioCorr_bars_w[index2] -> SetLineWidth(2);
					    fitFunc_energyRatioCorr_bars_w[index2] -> Draw("same");

			            gPad->Update();

			            line = new TLine(0.2,-5000,0.2,5000);
						line -> SetLineColor(kBlack);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						line -> Draw("same");
						
						line = new TLine(0,0,5,0);
						line -> SetLineColor(kRed);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						line -> Draw("same");
						
						line = new TLine(1,-5000,1,5000);
						line -> SetLineColor(kRed);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						line -> Draw("same");

					    latex = new TLatex(0.40,0.85,Form("#splitline{(w) bars %02d-%02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,iBar+1,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kRed);
					    latex -> Draw("same");
			  
					    c -> Print(Form("%s/energyRatioCorr/c_deltaT_energyRatio_w__%s.png",plotDir.c_str(),labelBars_energyBin.c_str()));
					    c -> Print(Form("%s/energyRatioCorr/c_deltaT_energyRatio_w__%s.pdf",plotDir.c_str(),labelBars_energyBin.c_str())); 
                        delete line;
						delete latex;
					    delete c;

						c = new TCanvas(Form("c_deltaT_energyRatio_prof_w_%s",labelBars_energyBin.c_str()),Form("c_deltaT_energyRatio_prof_w_%s",labelBars_energyBin.c_str()));
						prof -> GetYaxis()->SetRangeUser(-5000, +5000);
					    prof -> SetTitle(Form(";energy_{%02d} / energy_{%02d};#DeltaT_{%02d-%02d} [ps]",iBar,iBar+1,iBar,iBar+1));
						prof -> SetMarkerSize(0.4);
						prof -> Draw("psame");

						gStyle->SetOptFit(1111);
                        prof->SetStats(kTRUE); 

					    fitFunc_energyRatioCorr_bars_w[index2] -> Draw("same");

                        line = new TLine(0.2,-5000,0.2,5000);
						line -> SetLineColor(kBlack);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						line -> Draw("same");
			            

					    latex = new TLatex(0.40,0.85,Form("#splitline{(w) bars %02d-%02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,iBar+1,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kRed);
					    latex -> Draw("same");
			  
					    c -> Print(Form("%s/energyRatioCorr/c_deltaT_energyRatio_prof_w_%s.png",plotDir.c_str(),labelBars_energyBin.c_str()));
					    c -> Print(Form("%s/energyRatioCorr/c_deltaT_energyRatio_prof_w_%s.pdf",plotDir.c_str(),labelBars_energyBin.c_str())); 
                       
						delete latex;
					    delete c;
                        */
					    // --- draw energy correlation
						c = new TCanvas(Form("c_energyCorrelation_%s",labelLR_energyBin.c_str()),Form("c_energyCorrelation_%s",labelLR_energyBin.c_str()));
						if(!p1_energy_correlation[index2]) continue;
					    c -> SetGridy();
			  
					    h2 = h2_energy_correlation[index2]; 
					    h2 -> SetTitle(Form(";bar %02d energy [a.u.];bar %02d energy [a.u.]",iBar+1,iBar));
					    h2 -> Draw("colz");
					    prof = p1_energy_correlation[index2];
						prof -> SetMarkerSize(0.5);
					    prof -> Draw("psame");
                        //std::cout<<"my_cut["<<index2<<"]: "<<my_cut[index2]<<std::endl;
						line = new TLine(0,my_cut[index2],my_cut[index2],0);
						line -> SetLineColor(kRed);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						line -> Draw("same");
			           
					    latex = new TLatex(0.40,0.85,Form("V_{OV} = %.2f V, th. = %d DAC",Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kRed);
					    latex -> Draw("same");
			  
					    c -> Print(Form("%s/energyCorrelation/c_energyCorrelation_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
						c -> Print(Form("%s/energyCorrelation/c_energyCorrelation_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					    delete c;
					    delete latex;	
					}
			  
				} // --- end loop over energy bins

		    }// end loop over single double triple
	  
	    } // --- end loop ober bars
      
    } // --- end loop over stepLabels
  

  
    //------------------------
    //--- 3rd loop over events
	int accepted3=0;
	std::map<double,double> xMax_h2_tL_REFvseL;
	std::map<double,double> xMax_h2_tR_REFvseR;
	std::map<double, std::vector<int> > cut_ev_counterL; 
	std::map<double, std::vector<int> > cut_ev_counterR;
    for(auto mapIt : trees){
        ModuleEventClass* anEvent = new ModuleEventClass();
        mapIt.second -> SetBranchAddress("event",&anEvent);
      
        int nEntries = mapIt.second->GetEntries();
        for(int entry = 0; entry < nEntries; ++entry){
	        if( entry%100000 == 0 ){
	            std::cout << ">>> 3rd loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
	            //TrackProcess(cpu, mem, vsz, rss);
	        }
	        mapIt.second -> GetEntry(entry);
	  
	        bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
            if (!barFound) continue;
	  
	        int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
			int index3( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + 99 );
	  
	        if( !accept[index1][entry] ) continue;
	  
	        int energyBinAverage=0;
			if(anEvent->nClusters==1) {energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;	}
			else if(anEvent->nClusters==2 && anEvent->energyL_post>0 && anEvent->energyR_post>0){energyBinAverage = FindBin(anEvent->energySum,ranges_doubleHits[index1])+1;}
			else if(anEvent->nClusters==2 && anEvent->energyL_post<0 && anEvent->energyR_post<0){
				energyBinAverage = -1;
				continue; // escludo eventi "PRE"
			}
			else if(anEvent->nClusters==3){	energyBinAverage=1;}

	        double index2( (100000000*anEvent->nClusters)+(10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

			if(!fitFunc_energyRatio[index2]) continue;
	        float energyRatioMean = fitFunc_energyRatio[index2]->GetParameter(1);
	        float energyRatioSigma = fitFunc_energyRatio[index2]->GetParameter(2);
	        if(!fitFunc_totRatio[index2]) continue;
            float totRatioMean = fitFunc_totRatio[index2]->GetParameter(1);
	        float totRatioSigma = fitFunc_totRatio[index2]->GetParameter(2);
	  
	        if( fabs(anEvent->totR/anEvent->totL-totRatioMean) > 3.*totRatioSigma  ||  (anEvent->totR/anEvent->totL)>5 || (anEvent->totR/anEvent->totL)<0 ) {
	            accept[index1][entry] = false;
	            continue;
	        }
	        if(anEvent->nClusters==1){
				float energyMean = 0.5*(anEvent->energyR + anEvent->energyL );			
	        	if(!source.compare(TB) && energyMean < ranges["L-R"][index1]->at(0) ) {
	            	accept[index1][entry] = false;
	            	continue;
	        	}
			}
	        
	  
	        if( h1_deltaT[index2] == NULL ) {
	            std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
	      
	            h1_deltaT[index2] = new TH1F(Form("h1_deltaT_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
	            p1_deltaT_vs_energyRatio[index2] = new TProfile(Form("p1_deltaT_vs_energyRatio_%s",labelLR_energyBin.c_str()),"",50,energyRatioMean-3.*energyRatioSigma,energyRatioMean+3.*energyRatioSigma);
	            // p1_deltaT_vs_totRatio[index2] = new TProfile(Form("p1_deltaT_vs_totRatio_%s",labelLR_energyBin.c_str()),"",50,totRatioMean-3.*totRatioSigma, totRatioMean+3.*totRatioSigma);
	            p1_deltaT_vs_totRatio[index2] = new TProfile(Form("p1_deltaT_vs_totRatio_%s",labelLR_energyBin.c_str()),"",50,totRatioMean-5.*totRatioSigma, totRatioMean+5.*totRatioSigma);
	            h2_deltaT_vs_totRatio[index2] = new TH2F(Form("h2_deltaT_vs_totRatio_%s",labelLR_energyBin.c_str()),"",50,totRatioMean-3.*totRatioSigma, totRatioMean+3.*totRatioSigma, 2000, -12000., 12000.);
	        }

			if(h1_deltaT_raw_REF[index2] == NULL){ // delta T with REF module
				std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
				h1_deltaT_raw_REF[index2] = new TH1F(Form("h1_deltaT_raw_REF_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
				h1_deltaT_w_REF[index2] = new TH1F(Form("h1_deltaT_w_REF_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000);
			}

			if(h2_tL_REFvsEL[index2]==NULL){ // time L/R - time REF vs energy
				std::string labelL(Form("bar%02dL_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
				std::string labelR(Form("bar%02dR_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));

				h2_tL_REFvsEL[index2] = new TH2F(Form("h2_tL_REFvsEL_%s",labelL.c_str()),"",500,0,1000,2000,-24000,24000);
				p1_tL_REFvsEL[index2] = new TProfile(Form("p1_tL_REFvsEL_%s",labelL.c_str()),"",80,0,1000);
				h2_tR_REFvsER[index2] = new TH2F(Form("h2_tR_REFvsER_%s",labelR.c_str()),"",500,0,1000,2000,-24000,24000);
				p1_tR_REFvsER[index2] = new TProfile(Form("p1_tR_REFvsER_%s",labelR.c_str()),"",80,0,1000);

				h2_tL_wREFvsEL[index2] = new TH2F(Form("h2_tL_wREFvsEL_%s",labelL.c_str()),"",500,0,1000,2000,-24000,24000);
				p1_tL_wREFvsEL[index2] = new TProfile(Form("p1_tL_wREFvsEL_%s",labelL.c_str()),"",80,0,1000);
				h2_tR_wREFvsER[index2] = new TH2F(Form("h2_tR_wREFvsER_%s",labelR.c_str()),"",500,0,1000,2000,-24000,24000);
				p1_tR_wREFvsER[index2] = new TProfile(Form("p1_tR_wREFvsER_%s",labelR.c_str()),"",80,0,1000);

				h2_tL_REFvsE_REF[index2] = new TH2F(Form("h2_tL_REFvsE_REF_%s",labelL.c_str()),"",500,0,450,2000,-24000,24000);
				h2_tR_REFvsE_REF[index2] = new TH2F(Form("h2_tR_REFvsE_REF_%s",labelR.c_str()),"",500,0,450,2000,-24000,24000);
			}

			if(h2_deltaT_vs_energyRatio_REF[index3]==NULL){ // delta T  vs energy Ratio of reference module 
				h2_deltaT_vs_energyRatio_REF[index3] = new TH2F(Form("h2_deltaT_vs_energyRatio_REF_Vov%.2f_th%02d",anEvent->Vov,anEvent->vth1),"",50,0,5,2000,-12000,12000);
				p1_deltaT_vs_energyRatio_REF[index3] = new TProfile(Form("p1_deltaT_vs_energyRatio_REF_Vov%.2f_th%02d",anEvent->Vov,anEvent->vth1),"",50,0,5);
			}
            
			long long deltaT = anEvent->timeR - anEvent->timeL;


			if(fabs(anEvent->timeL_ext - anEvent->timeR_ext)<10000 && anEvent->energyR_ext>0 && anEvent->energyL_ext>0){
				// --- deltaT vs energy ratio of reference module
				h2_deltaT_vs_energyRatio_REF[index3] -> Fill(anEvent->energyR_ext/anEvent->energyL_ext, anEvent->timeL_ext - anEvent->timeR_ext);
				p1_deltaT_vs_energyRatio_REF[index3] -> Fill(anEvent->energyR_ext/anEvent->energyL_ext, anEvent->timeL_ext - anEvent->timeR_ext);
			}

	        if(fabs(deltaT)>10000) continue;
	        accepted3++;
	        h1_deltaT[index2] -> Fill( deltaT );    
            float timeLow = CTRMeans[index2] - 3.* CTRSigmas[index2];
	        float timeHig = CTRMeans[index2] + 3.* CTRSigmas[index2];
	  
	        if( ( deltaT > timeLow ) && ( deltaT < timeHig ) ){
	            p1_deltaT_vs_energyRatio[index2] -> Fill( anEvent->energyR/anEvent->energyL,deltaT );
	            p1_deltaT_vs_totRatio[index2] -> Fill( anEvent->totR/anEvent->totL,deltaT );
	            h2_deltaT_vs_totRatio[index2] -> Fill( anEvent->totR/anEvent->totL,deltaT );
	        }

            if (cut_ev_counterL.count(index2) == 0) {
                cut_ev_counterL[index2] = std::vector<int>(2, 0);  // [0, 0]
				cut_ev_counterR[index2] = std::vector<int>(2, 0);  // [0, 0]
            }


			if(anEvent->nClusters==2 && anEvent->energyL_post>0 && anEvent->energyR_post>0){
				if(anEvent->barID==15) continue;
				double E_mean=0.5*(anEvent->energyL+anEvent->energyR);
				double E_mean_post=0.5*(anEvent->energyL_post+anEvent->energyR_post);
				double E_mean_REF=0.5*(anEvent->energyL_ext+anEvent->energyR_ext);// energy REF module

				double t_mean_a=0.5*(static_cast<double>(anEvent->timeL)+static_cast<double>(anEvent->timeR));//average, arithmetic mean
				double t_mean_post_a=0.5*(static_cast<double>(anEvent->timeL_post)+static_cast<double>(anEvent->timeR_post));

				double t_mean_w=((static_cast<double>(anEvent->energyL)*static_cast<double>(anEvent->timeL))+(static_cast<double>(anEvent->energyR)*static_cast<double>(anEvent->timeR)))/(static_cast<double>(anEvent->energyL)+static_cast<double>(anEvent->energyR));// weighted mean
				double t_mean_post_w=((static_cast<double>(anEvent->energyL_post)*static_cast<double>(anEvent->timeL_post))+(static_cast<double>(anEvent->energyR_post)*static_cast<double>(anEvent->timeR_post)))/(static_cast<double>(anEvent->energyL_post)+static_cast<double>(anEvent->energyR_post));
				
                if(fabs(t_mean_w-t_mean_post_w)<10000){
					float my_timeLow = my_CTRMeans[index2] - 3.* my_CTRSigmas[index2];
	                float my_timeHig = my_CTRMeans[index2] + 3.* my_CTRSigmas[index2];
					if( ( t_mean_w-t_mean_post_w > my_timeLow ) && ( t_mean_w-t_mean_post_w < my_timeHig ) ){
						h2_deltaT_energyRatio_bars_w[index2]->Fill(E_mean/E_mean_post,t_mean_w-t_mean_post_w);
                        p1_deltaT_energyRatio_bars_w[index2]->Fill(E_mean/E_mean_post,t_mean_w-t_mean_post_w);
					}
				}

				if(fabs(t_mean_a-t_mean_post_a)<10000 ){
					if((E_mean+E_mean_post)<=my_cut[index2]){
                        if(!energyRatio_ranges[index2] || !energyRatio_ranges[index2+1]) continue;
                        if(!fitFunc_energyRatioCorr_bars[index2] ) continue;
                        if(!h1_deltaT_bars[index2]) continue;

						
						if((anEvent->energyR/anEvent->energyL > energyRatio_ranges[index2]->at(0)) && (anEvent->energyR/anEvent->energyL < energyRatio_ranges[index2]->at(1) ) &&
					       (anEvent->energyR_post/anEvent->energyL_post)>energyRatio_ranges[index2+1]->at(0) && (anEvent->energyR_post/anEvent->energyL_post < energyRatio_ranges[index2+1]->at(1))){

                            double energyRatioCorr_bars = fitFunc_energyRatioCorr_bars[index2]->Eval(E_mean/E_mean_post) - fitFunc_energyRatioCorr_bars[index2]->Eval(1.);
							//double energyRatioCorr_bars_w = fitFunc_energyRatioCorr_bars_w[index2]->Eval(E_mean/E_mean_post) - fitFunc_energyRatioCorr_bars_w[index2]->Eval(1.);
                            h1_deltaT_bars[index2]->Fill(t_mean_a - t_mean_post_a - energyRatioCorr_bars);
							//h1_deltaT_bars_w[index2]->Fill(t_mean_w - t_mean_post_w - energyRatioCorr_bars_w);
						}
					}			
				}

				double time_ave = 0.5*(t_mean_a+t_mean_post_a);
				double time_ave_w = (t_mean_a*E_mean + t_mean_post_a*E_mean_post)/(E_mean+E_mean_post);
				double time_ave_REF = 0.5*(static_cast<double>(anEvent->timeL_ext)+static_cast<double>(anEvent->timeR_ext));// time REF module 
				double time_ave_wREF= (static_cast<double>(anEvent->energyL_ext)*static_cast<double>(anEvent->timeL_ext)+static_cast<double>(anEvent->energyR_ext)*static_cast<double>(anEvent->timeR_ext))/(static_cast<double>(anEvent->energyL_ext)+static_cast<double>(anEvent->energyR_ext));// weighted mean REF module
				
				if(anEvent->energyR_ext>0 && anEvent->energyL_ext>0){
					h1_deltaT_raw_REF[index2] -> Fill(time_ave - time_ave_REF);
				    h1_deltaT_w_REF[index2] -> Fill( time_ave_w - time_ave_REF);

				    h2_tL_REFvsE_REF[index2] -> Fill(E_mean_REF, anEvent->timeL - time_ave_REF);
				    h2_tR_REFvsE_REF[index2] -> Fill(E_mean_REF, anEvent->timeR - time_ave_REF);

				    h2_tL_wREFvsEL[index2] ->Fill(anEvent->energyL, anEvent->timeL - time_ave_wREF);
				    p1_tL_wREFvsEL[index2] ->Fill(anEvent->energyL, anEvent->timeL - time_ave_wREF);

				    h2_tR_wREFvsER[index2] ->Fill(anEvent->energyR, anEvent->timeR - time_ave_wREF);
				    p1_tR_wREFvsER[index2] ->Fill(anEvent->energyR, anEvent->timeR - time_ave_wREF);

                    if((anEvent->timeL - time_ave_REF) >3500){ 
					    h2_tL_REFvsEL[index2] -> Fill(anEvent->energyL, anEvent->timeL - time_ave_REF);
				        p1_tL_REFvsEL[index2] -> Fill(anEvent->energyL, anEvent->timeL - time_ave_REF);
					    cut_ev_counterL[index2][0]++;
				    }
				    else if((anEvent->timeL - time_ave_REF) <= 3500) cut_ev_counterL[index2][1]++;
                    if((anEvent->timeR - time_ave_REF) >3500){
				        h2_tR_REFvsER[index2] -> Fill(anEvent->energyR, anEvent->timeR - time_ave_REF);
				        p1_tR_REFvsER[index2] -> Fill(anEvent->energyR, anEvent->timeR - time_ave_REF);
					    cut_ev_counterR[index2][0]++;
			        }
				    else if((anEvent->timeR - time_ave_REF) <= 3500) cut_ev_counterR[index2][1]++;
				}
				
			}
	    }
      
        std::cout << std::endl;
    }
  
  
    //------------------
    //--- draw 3rd plots
    std::map<double,TF1*> fitFunc_energyRatioCorr;
    std::map<double,TF1*> fitFunc_totRatioCorr;
    std::map<double,TF1*> fitFunc_energyRatioCorr_totRatioCorr;
	std::cout<<"3rd loop accepted events: "<<accepted3<<std::endl;
    std::cout<<"\nDrawing 3rd plots ... "<<std::endl;
    for(auto stepLabel : stepLabels){
        float Vov = map_Vovs[stepLabel];
        float vth1 = map_ths[stepLabel];

		std::string extLabel(Form("externalBar_L-R_%s",stepLabel.c_str())); 
		int index3( (10000*int(Vov*100.)) + (100*vth1) + 99 );
      

		// -- draw deltaT vs energy ratio for external bar
		c = new TCanvas(Form("c_deltaT_vs_energyRatio_%s",extLabel.c_str()),Form("c_deltaT_vs_energyRatio_%s",extLabel.c_str()));
		h2 = h2_deltaT_vs_energyRatio_REF[index3];
		h2 -> SetTitle(Form(";energy_{right} / energy_{left};#Deltat [ps]"));
		h2 -> GetYaxis() -> SetRangeUser(-1000,1000);
		h2 -> Draw("colz");

		prof = p1_deltaT_vs_energyRatio_REF[index3];
		prof -> SetMarkerSize(0.4);
		prof -> Draw("psame");
			
		latex = new TLatex(0.40,0.85,Form("#splitline{external bar}{V_{OV} = %.2f V, th. = %d DAC}",Vov,int(vth1)));
		latex -> SetNDC();
		latex -> SetTextFont(42);
		latex -> SetTextSize(0.04);
		latex -> SetTextColor(kRed);
		latex -> Draw("same");
			
		float fitXMin = fitFunc_energyRatio_REF[index3]->GetParameter(1) - 1.*fitFunc_energyRatio_REF[index3]->GetParameter(2);
		float fitXMax = fitFunc_energyRatio_REF[index3]->GetParameter(1) + 1.*fitFunc_energyRatio_REF[index3]->GetParameter(2);
			
		fitFunc_energyRatioCorrection_REF[index3] = new TF1(Form("fitFunc_energyRatioCorr_%s",extLabel.c_str()),"pol3",fitXMin,fitXMax);
		prof -> Fit(fitFunc_energyRatioCorrection_REF[index3],"QRS+");
		fitFunc_energyRatioCorrection_REF[index3] -> SetLineColor(kRed);
		fitFunc_energyRatioCorrection_REF[index3] -> SetLineWidth(2);
		fitFunc_energyRatioCorrection_REF[index3] -> Draw("same");
			
		c -> Print(Form("%s/externalBar/c_deltaT_vs_energyRatio__%s.png",plotDir.c_str(),extLabel.c_str()));
		c -> Print(Form("%s/externalBar/c_deltaT_vs_energyRatio__%s.pdf",plotDir.c_str(),extLabel.c_str()));
		delete c;

		c = new TCanvas(Form("c_deltaT_vs_energyRatio_p_%s",extLabel.c_str()),Form("c_deltaT_vs_energyRatio_p_%s",extLabel.c_str()));
        prof = p1_deltaT_vs_energyRatio_REF[index3];
		prof -> SetMarkerSize(0.4);
		prof -> SetTitle(Form(";energy_{right} / energy_{left};#Deltat [ps]"));
		prof -> GetYaxis() -> SetRangeUser(-500,500);
		prof -> Draw("psame");
		latex -> Draw("same");
		fitFunc_energyRatioCorrection_REF[index3] -> Draw("same");
		TLine *line = new TLine(fitXMin,-500,fitXMin,500);
		line -> SetLineColor(kBlue);
		line -> SetLineWidth(2);
		line -> SetLineStyle(2);
		line -> Draw("same");
		line = new TLine(fitXMax,-500,fitXMax,500);
		line -> SetLineColor(kBlue);
		line -> SetLineWidth(2);
		line -> SetLineStyle(2);
		line -> Draw("same");
		c -> Print(Form("%s/externalBar/c_deltaT_vs_energyRatio_p_%s.png",plotDir.c_str(),extLabel.c_str()));
		c -> Print(Form("%s/externalBar/c_deltaT_vs_energyRatio_p_%s.pdf",plotDir.c_str(),extLabel.c_str()));
		delete c;

        for(int iBar = 0; iBar < 16; ++iBar){
	
            bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;                                                                                       
	        if (!barFound) continue;  
	
	        std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
			std::string labelBars(Form("bars%02d-%02d_%s",iBar,iBar+1,stepLabel.c_str()));
	        int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );	
	        if( !ranges["L-R"][index1] ) continue;
	
	        int nEnergyBins = ranges["L-R"][index1]->size()-1;
	
	        for(int i=1;i<=3;i++){//loop over single double triple 
				for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin) {
					//if (ranges["L-R"][index1]->at(iEnergyBin)<0) continue;
					double  index2((100000000*i)+10000000*iEnergyBin+index1 );
					if(!p1_deltaT_vs_energyRatio[index2]) continue;

					std::string labelLR_energyBin;
					std::string labelBars_energyBin;
					if(i==1){labelLR_energyBin = Form("%s_energyBin%02d_S",labelLR.c_str(),iEnergyBin);}
					else if(i==2){
						labelLR_energyBin = Form("%s_energyBin%02d_D",labelLR.c_str(),iEnergyBin);
						labelBars_energyBin = Form("%s_energyBin%02d_D",labelBars.c_str(),iEnergyBin);
					}
					else if(i==3){labelLR_energyBin = Form("%s_energyBin%02d_T",labelLR.c_str(),iEnergyBin);}
			
			
					// -- draw deltaT vs energy ratio
					c = new TCanvas(Form("c_deltaT_vs_energyRatio_%s",labelLR_energyBin.c_str()),Form("c_deltaT_vs_energyRatio_%s",labelLR_energyBin.c_str()));
			
					prof = p1_deltaT_vs_energyRatio[index2];
					prof -> SetTitle(Form(";energy_{right} / energy_{left};#Deltat [ps]"));
					prof -> GetYaxis() -> SetRangeUser(CTRMeans[index2]-3.*CTRSigmas[index2],CTRMeans[index2]+3.*CTRSigmas[index2]);
					prof -> Draw("");
			
					latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					latex -> SetNDC();
					latex -> SetTextFont(42);
					latex -> SetTextSize(0.04);
					latex -> SetTextColor(kRed);
					latex -> Draw("same");
			
					float fitXMin = fitFunc_energyRatio[index2]->GetParameter(1) - 3.*fitFunc_energyRatio[index2]->GetParameter(2);
					float fitXMax = fitFunc_energyRatio[index2]->GetParameter(1) + 3.*fitFunc_energyRatio[index2]->GetParameter(2);
			
					fitFunc_energyRatioCorr[index2] = new TF1(Form("fitFunc_energyRatioCorr_%s",labelLR_energyBin.c_str()),"pol3",fitXMin,fitXMax);
					prof -> Fit(fitFunc_energyRatioCorr[index2],"QRS+");
					fitFunc_energyRatioCorr[index2] -> SetLineColor(kRed);
					fitFunc_energyRatioCorr[index2] -> SetLineWidth(2);
					fitFunc_energyRatioCorr[index2] -> Draw("same");
			
					c -> Print(Form("%s/energyRatioCorr/c_deltaT_vs_energyRatio__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/energyRatioCorr/c_deltaT_vs_energyRatio__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete latex;
					delete c;
			
			
					// -- draw deltaT vs tot ratio
					if(!p1_deltaT_vs_totRatio[index2]) continue;
			
					c = new TCanvas(Form("c_deltaT_vs_totRatio_%s",labelLR_energyBin.c_str()),Form("c_deltaT_vs_totRatio_%s",labelLR_energyBin.c_str()));
			
					prof = p1_deltaT_vs_totRatio[index2];
					prof -> SetTitle(Form(";ToT_{right} / ToT_{left};#Deltat [ps]"));
					prof -> GetYaxis() -> SetRangeUser(CTRMeans[index2]-3.*CTRSigmas[index2],CTRMeans[index2]+3.*CTRSigmas[index2]);
					//prof -> GetXaxis() -> SetRangeUser(-1.2,1.2);
					prof -> Draw("");
			
					latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					latex -> SetNDC();
					latex -> SetTextFont(42);
					latex -> SetTextSize(0.04);
					latex -> SetTextColor(kRed);
					latex -> Draw("same");
			
					fitXMin = fitFunc_totRatio[index2]->GetParameter(1) - 3.*fitFunc_totRatio[index2]->GetParameter(2);
					fitXMax = fitFunc_totRatio[index2]->GetParameter(1) + 3.*fitFunc_totRatio[index2]->GetParameter(2);
			
					fitFunc_totRatioCorr[index2] = new TF1(Form("fitFunc_totRatioCorr_%s",labelLR_energyBin.c_str()),"pol3",fitXMin,fitXMax);
					prof -> Fit(fitFunc_totRatioCorr[index2],"QRS+");
					fitFunc_totRatioCorr[index2] -> SetLineColor(kRed);
					fitFunc_totRatioCorr[index2] -> SetLineWidth(2);
					fitFunc_totRatioCorr[index2] -> Draw("same");
			
					c -> Print(Form("%s/totRatioCorr/c_deltaT_vs_totRatio__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/totRatioCorr/c_deltaT_vs_totRatio__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete c;
					delete latex;

                    if(i==2){
						// --- draw deltaT_bars vs energy ratio for weigthed mean : start
						c = new TCanvas(Form("c_deltaT_energyRatio_w_%s",labelBars_energyBin.c_str()),Form("c_deltaT_energyRatio_w_%s",labelBars_energyBin.c_str()));
						h2 = h2_deltaT_energyRatio_bars_w[index2];
						h2 -> GetYaxis()->SetRangeUser(-5000, +5000);
						h2 -> SetTitle(Form(";energy_{%02d} / energy_{%02d};#DeltaT_{%02d-%02d} [ps]",iBar,iBar+1,iBar,iBar+1));
						h2 -> Draw("colz");
						outFile -> cd();
						h2 -> Write();
						prof = p1_deltaT_energyRatio_bars_w[index2];
						prof -> SetMarkerSize(0.4);
						prof -> Draw("psame");
						fitFunc_energyRatioCorr_bars_w[index2] = new TF1(Form("fitFunc_energyRatioCorr_%s",labelLR_energyBin.c_str()),"[0]/x + [1]/(x*x) + [2]/(x*x*x)- [3]*x - [4]*x*x -[5]*x*x*x",0.2,5);
						//fitFunc_energyRatioCorr_bars_w[index2] = new TF1(Form("fitFunc_energyRatioCorr_%s",labelLR_energyBin.c_str()),"[0]/x - [1]*x - [2]*x*x -[3]*x*x*x",0.2,5);
						prof -> Fit(fitFunc_energyRatioCorr_bars_w[index2],"QRS");
						fitFunc_energyRatioCorr_bars_w[index2] -> SetLineColor(kRed);
						fitFunc_energyRatioCorr_bars_w[index2] -> SetLineWidth(2);
						fitFunc_energyRatioCorr_bars_w[index2] -> Draw("same");
						gPad->Update();

						line = new TLine(0.2,-5000,0.2,5000);
						line -> SetLineColor(kBlack);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						line -> Draw("same");
						
						line = new TLine(0,0,5,0);
						line -> SetLineColor(kRed);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						line -> Draw("same");
						
						line = new TLine(1,-5000,1,5000);
						line -> SetLineColor(kRed);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						line -> Draw("same");

						latex = new TLatex(0.40,0.85,Form("#splitline{(w) bars %02d-%02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,iBar+1,Vov,int(vth1)));
						latex -> SetNDC();
						latex -> SetTextFont(42);
						latex -> SetTextSize(0.04);
						latex -> SetTextColor(kRed);
						latex -> Draw("same");
			  
						c -> Print(Form("%s/energyRatioCorr/c_deltaT_energyRatio_w__%s.png",plotDir.c_str(),labelBars_energyBin.c_str()));
						c -> Print(Form("%s/energyRatioCorr/c_deltaT_energyRatio_w__%s.pdf",plotDir.c_str(),labelBars_energyBin.c_str())); 
						delete line;
						delete latex;
						delete c;

						c = new TCanvas(Form("c_deltaT_energyRatio_prof_w_%s",labelBars_energyBin.c_str()),Form("c_deltaT_energyRatio_prof_w_%s",labelBars_energyBin.c_str()));
						prof -> GetYaxis() -> SetRangeUser(my_CTRMeans[index2]-3.*my_CTRSigmas[index2],my_CTRMeans[index2]+3.*my_CTRSigmas[index2]);
						prof -> SetTitle(Form(";energy_{%02d} / energy_{%02d};#DeltaT_{%02d-%02d} [ps]",iBar,iBar+1,iBar,iBar+1));
						prof -> SetMarkerSize(0.4);
						prof -> Draw("psame");

						gStyle->SetOptFit(1111);
						prof->SetStats(kTRUE); 

						fitFunc_energyRatioCorr_bars_w[index2] -> Draw("same");

						line = new TLine(0.2,-5000,0.2,5000);
						line -> SetLineColor(kBlack);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						line -> Draw("same");
			            

						latex = new TLatex(0.40,0.85,Form("#splitline{(w) bars %02d-%02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,iBar+1,Vov,int(vth1)));
						latex -> SetNDC();
						latex -> SetTextFont(42);
						latex -> SetTextSize(0.04);
						latex -> SetTextColor(kRed);
						latex -> Draw("same");
			  
						c -> Print(Form("%s/energyRatioCorr/c_deltaT_energyRatio_prof_w_%s.png",plotDir.c_str(),labelBars_energyBin.c_str()));
						c -> Print(Form("%s/energyRatioCorr/c_deltaT_energyRatio_prof_w_%s.pdf",plotDir.c_str(),labelBars_energyBin.c_str())); 
                       
						delete latex;
						delete c;
                        // --- draw deltaT_bars vs energy ratio for weigthed mean : end

						c = new TCanvas(Form("c_deltaT_%s",labelBars_energyBin.c_str()),Form("c_deltaT_%s",labelBars_energyBin.c_str()));
						histo = h1_deltaT_bars[index2];
					    histo -> SetLineWidth(2);
					    histo -> SetLineColor(kBlue);
					    histo -> SetMarkerColor(kBlue);
			            TF1* fitFunc = new TF1(Form("fitFunc_deltaT_%s",labelBars_energyBin.c_str()),"gaus",-10000, 10000);
					    drawDeltaT(c, histo, fitFunc, Form("%02d - %02d diff ",iBar,iBar+1), "enCorr","");
					    outFile -> cd();
					    histo -> Write();

						c -> Print(Form("%s/CTR_energyRatioCorr/c_deltaT__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
				    	c -> Print(Form("%s/CTR_energyRatioCorr/c_deltaT__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
						delete c;
                        /*
						c = new TCanvas(Form("c_deltaT_w_%s",labelBars_energyBin.c_str()),Form("c_deltaT_w_%s",labelBars_energyBin.c_str()));
						histo = h1_deltaT_bars_w[index2];
					    histo -> SetLineWidth(2);
					    histo -> SetLineColor(kGreen);
					    histo -> SetMarkerColor(kGreen);
					    drawDeltaT(c, histo, fitFunc, Form("%02d - %02d diff ",iBar,iBar+1), "enCorr_w","");
					    outFile -> cd();
					    histo -> Write();

						c -> Print(Form("%s/CTR_energyRatioCorr/c_deltaT__w_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
				    	c -> Print(Form("%s/CTR_energyRatioCorr/c_deltaT__w_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
						delete fitFunc;
						delete c;
						*/

						c = new TCanvas(Form("c_deltaT_REF_raw_%s",labelBars_energyBin.c_str()),Form("c_deltaT_REF_raw_%s",labelBars_energyBin.c_str()));
						histo = h1_deltaT_raw_REF[index2];
					    histo -> SetLineWidth(2);
					    histo -> SetLineColor(kAzure);
					    histo -> SetMarkerColor(kAzure);
			            fitFunc = new TF1(Form("fitFunc_deltaT_raw_REF_%s",labelBars_energyBin.c_str()),"gaus",-10000, 10000);
					    drawDeltaT(c, histo, fitFunc, Form("[(%02d + %02d) - REF]  ",iBar,iBar+1), "raw","");
					    outFile -> cd();
					    histo -> Write();

						c -> Print(Form("%s/CTR_REF/CTR_REF_raw/c_deltaT_REF_raw_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
				    	c -> Print(Form("%s/CTR_REF/CTR_REF_raw/c_deltaT_REF_raw_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
						delete c;
                        
						c = new TCanvas(Form("c_deltaT_REF_raw_w_%s",labelBars_energyBin.c_str()),Form("c_deltaT_REF_raw_%s",labelBars_energyBin.c_str()));
						histo = h1_deltaT_w_REF[index2];
					    histo -> SetLineWidth(2);
					    histo -> SetLineColor(kRed+2);
					    histo -> SetMarkerColor(kRed+2);
			            fitFunc = new TF1(Form("fitFunc_deltaT_raw_REF_w_%s",labelBars_energyBin.c_str()),"gaus",-10000, 10000);
					    drawDeltaT(c, histo, fitFunc, Form("[(%02d + %02d) - REF]  ",iBar,iBar+1), "raw(w)","");
					    outFile -> cd();
					    histo -> Write();

						c -> Print(Form("%s/CTR_REF/CTR_REF_raw/c_deltaT_REF_w_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
				    	c -> Print(Form("%s/CTR_REF/CTR_REF_raw/c_deltaT_REF_w_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
						delete c;
                      
                        // --- draw timeL-time_REF vs energyL
						c = new TCanvas(Form("c_tL_REFvsEL_bar%02dL_%s_D",iBar,stepLabel.c_str()),Form("c_tL_REFvsEL_bar%02dL_%s_D",iBar,stepLabel.c_str()));
						h2 = h2_tL_REFvsEL[index2];
						int bx,by,bz;
						h2->GetMaximumBin(bx,by,bz);
						double xMax = h2->GetXaxis()->GetBinCenter(bx);
						xMax_h2_tL_REFvseL[index2] = xMax;
						h2 -> GetYaxis()->SetRangeUser(0, 10000);
                        h2 -> SetTitle(Form(";energy_{L} [a.u.]; time_{L} - time_{REF} [ps]"));
						h2->Draw("colz");
						prof = p1_tL_REFvsEL[index2];
						prof -> SetMarkerSize(0.4);
						prof -> Draw("psame");
						fitFunc_timeL[index2] = new TF1(Form("fitFunc_timeL_bar%02d_%s",iBar,stepLabel.c_str()),"[0]/x + [1]/(x*x) - [2]*x - [3]*x*x",80,800);
						prof->Fit(fitFunc_timeL[index2],"QRS");
						fitFunc_timeL[index2] -> SetLineColor(kRed);
						fitFunc_timeL[index2] -> SetLineWidth(2);
						fitFunc_timeL[index2] -> Draw("same");

                        double ev_below_cut = static_cast<double>(cut_ev_counterL[index2][1])/(cut_ev_counterL[index2][0]+cut_ev_counterL[index2][1]) *100. ; 
						latex = new TLatex(0.40,0.73,Form("below cut: %.2f%",ev_below_cut));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.03);
					    latex -> SetTextColor(kBlue+2);
					    latex -> Draw("same");

						latex = new TLatex(0.40,0.8,Form("#splitline{bar %02dL}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kRed);
					    latex -> Draw("same");

						c->Print(Form("%s/timeCorrelation/c_tL_REFvsEL_bar%02dL_%s_D.png",plotDir.c_str(),iBar,stepLabel.c_str()));
						c->Print(Form("%s/timeCorrelation/c_tL_REFvsEL_bar%02dL_%s_D.pdf",plotDir.c_str(),iBar,stepLabel.c_str()));
						delete c;
                        // --- profile tL_eL
						c = new TCanvas(Form("c_tL_REFvsEL_p_bar%02dL_%s_D",iBar,stepLabel.c_str()),Form("c_tL_REFvsEL_p_bar%02dL_%s_D",iBar,stepLabel.c_str()));
						prof -> GetYaxis()->SetRangeUser(2500, 10000);
						prof -> SetTitle(Form(";energy_{L} [a.u.]; time_{L} - time_{REF} [ps]"));
						prof -> SetMarkerSize(0.4);
						prof -> Draw("psame");

						gStyle->SetOptFit(1111);
                        prof->SetStats(kTRUE); 
						fitFunc_timeL[index2] -> Draw("same");
						latex -> Draw("same");
						
						TLine* line_h = new TLine(0,fitFunc->GetParameter(1),1000,fitFunc->GetParameter(1));
						line_h -> SetLineColor(kBlue);
						line_h -> SetLineWidth(2);
						line_h -> SetLineStyle(2);
						line_h -> Draw("same");

						TLine* line_v = new TLine(xMax,2500,xMax,10000);
						line_v -> SetLineColor(kBlue);
						line_v -> SetLineWidth(2);
						line_v -> SetLineStyle(2);
						line_v -> Draw("same");

						TLine* line_80 = new TLine(80,2500,80,10000);
						line_80 -> SetLineColor(kBlack);
						line_80 -> SetLineWidth(2);
						line_80 -> SetLineStyle(2);
						line_80 -> Draw("same");

						TLine* line_800 = new TLine(800,2500,800,10000);
						line_800 -> SetLineColor(kBlack);
						line_800 -> SetLineWidth(2);
						line_800 -> SetLineStyle(2);
						line_800 -> Draw("same");
						

						c->Print(Form("%s/timeCorrelation/c_tL_REFvsEL_p_bar%02dL_%s_D.png",plotDir.c_str(),iBar,stepLabel.c_str()));
						c->Print(Form("%s/timeCorrelation/c_tL_REFvsEL_p_bar%02dL_%s_D.pdf",plotDir.c_str(),iBar,stepLabel.c_str()));
						
						delete latex;
						delete c;

                        
						// --- draw timeR-time_REF vs energyR
						c = new TCanvas(Form("c_tR_REFvsER_bar%02dR_%s_D",iBar,stepLabel.c_str()),Form("c_tR_REFvsER_bar%02dR_%s_D",iBar,stepLabel.c_str()));
						h2 = h2_tR_REFvsER[index2];
						h2->GetMaximumBin(bx,by,bz);
						xMax = h2->GetXaxis()->GetBinCenter(bx);
						xMax_h2_tR_REFvseR[index2] = xMax;
						h2 -> GetYaxis()->SetRangeUser(0, 10000);
                        h2 -> SetTitle(Form(";energy_{R} [a.u.];time_{R} - time_{REF} [ps]"));
						h2->Draw("colz");
						prof = p1_tR_REFvsER[index2];
						prof -> SetMarkerSize(0.4);
						prof -> Draw("psame");

						fitFunc_timeR[index2] = new TF1(Form("fitFunc_timeR_bar%02dR_%s",iBar,stepLabel.c_str()),"[0]/x + [1]/(x*x) - [2]*x - [3]*x*x",80,600);
						prof->Fit(fitFunc_timeR[index2],"QRS");
						fitFunc_timeR[index2] -> SetLineColor(kRed);
						fitFunc_timeR[index2] -> SetLineWidth(2);
						fitFunc_timeR[index2] -> Draw("same");

						ev_below_cut = static_cast<double>(cut_ev_counterR[index2][1])/(cut_ev_counterR[index2][0]+cut_ev_counterR[index2][1]) *100. ; 
						latex = new TLatex(0.40,0.73,Form("below cut: %.2f",ev_below_cut));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.03);
					    latex -> SetTextColor(kBlue+2);
					    latex -> Draw("same");

						latex = new TLatex(0.40,0.8,Form("#splitline{bar %02dR}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kRed);
					    latex -> Draw("same");

						c->Print(Form("%s/timeCorrelation/c_tR_REFvsER_bar%02dR_%s_D.png",plotDir.c_str(),iBar,stepLabel.c_str()));
						c->Print(Form("%s/timeCorrelation/c_tR_REFvsER_bar%02dR_%s_D.pdf",plotDir.c_str(),iBar,stepLabel.c_str()));
						delete c;
                        
						// --- profile tR_eR
						c = new TCanvas(Form("c_tR_REFvsER_p_bar%02dR_%s_D",iBar,stepLabel.c_str()),Form("c_tR_REFvsER_p_bar%02dR_%s_D",iBar,stepLabel.c_str()));
						prof -> GetYaxis()->SetRangeUser(2500, 10000);
						prof -> SetTitle(Form(";energy_{R} [a.u.]; time_{R} - time_{REF} [ps]"));
						prof -> SetMarkerSize(0.4);
						prof -> Draw("psame");

						gStyle -> SetOptFit(1111);
                        prof -> SetStats(kTRUE); 
						fitFunc_timeR[index2] -> Draw("same");
						latex -> Draw("same");
						line_h -> Draw("same");
						line_v = new TLine(xMax,2500,xMax,10000);
						line_v -> SetLineColor(kBlue);
						line_v -> SetLineWidth(2);
						line_v -> SetLineStyle(2);
						line_v -> Draw("same");
						line_80 -> Draw("same");
						line_800 = new TLine(600,2500,600,10000);
						line_800 -> SetLineColor(kBlack);
						line_800 -> SetLineWidth(2);
						line_800 -> SetLineStyle(2);
						line_800 -> Draw("same");


						c->Print(Form("%s/timeCorrelation/c_tR_REFvsER_p_bar%02dR_%s_D.png",plotDir.c_str(),iBar,stepLabel.c_str()));
						c->Print(Form("%s/timeCorrelation/c_tR_REFvsER_p_bar%02dR_%s_D.pdf",plotDir.c_str(),iBar,stepLabel.c_str()));
						delete line_h;
						delete line_v;
						delete line_80;
						delete line_800;
						delete latex;
						delete c;

						// --- tL-REF(W) vs EL
						c = new TCanvas(Form("c_tL_wREFvsEL_bar%02dL_%s_D",iBar,stepLabel.c_str()),Form("c_tL_wREFvsEL_bar%02dL_%s_D",iBar,stepLabel.c_str()));
						h2 = h2_tL_wREFvsEL[index2];
						h2 -> GetYaxis()->SetRangeUser(0, 10000);
                        h2 -> SetTitle(Form(";energy_{L} [a.u.]; time_{L} - time_{REF(w)} [ps]"));
						h2->Draw("colz");
						prof = p1_tL_wREFvsEL[index2];
						prof -> SetMarkerSize(0.4);
						prof -> Draw("psame");

						latex = new TLatex(0.40,0.8,Form("#splitline{bar %02dL}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kRed);
					    latex -> Draw("same");

						c->Print(Form("%s/timeCorrelation/w_REF/c_tL_wREFvsEL_bar%02dL_%s_D.png",plotDir.c_str(),iBar,stepLabel.c_str()));
						c->Print(Form("%s/timeCorrelation/w_REF/c_tL_wREFvsEL_bar%02dL_%s_D.pdf",plotDir.c_str(),iBar,stepLabel.c_str()));
						delete c;

						// --- tR-REF(W) vs ER
						c = new TCanvas(Form("c_tR_wREFvsER_bar%02dR_%s_D",iBar,stepLabel.c_str()),Form("c_tR_wREFvsER_bar%02dR_%s_D",iBar,stepLabel.c_str()));
						h2 = h2_tR_wREFvsER[index2];
						h2 -> GetYaxis()->SetRangeUser(0, 10000);
                        h2 -> SetTitle(Form(";energy_{R} [a.u.]; time_{R} - time_{REF(w)} [ps]"));
						h2->Draw("colz");
						prof = p1_tR_wREFvsER[index2];
						prof -> SetMarkerSize(0.4);
						prof -> Draw("psame");

						latex = new TLatex(0.40,0.8,Form("#splitline{bar %02dR}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kRed);
					    latex -> Draw("same");

						c->Print(Form("%s/timeCorrelation/w_REF/c_tR_wREFvsER_bar%02dR_%s_D.png",plotDir.c_str(),iBar,stepLabel.c_str()));
						c->Print(Form("%s/timeCorrelation/w_REF/c_tR_wREFvsER_bar%02dR_%s_D.pdf",plotDir.c_str(),iBar,stepLabel.c_str()));
						delete c;

						// --- tL-REF vs eREF
						c = new TCanvas(Form("c_tL_REFvseREF_bar%02dL_%s_D",iBar,stepLabel.c_str()),Form("c_tL_REFvseREF_bar%02dL_%s_D",iBar,stepLabel.c_str()));
						h2 = h2_tL_REFvsE_REF[index2];
						h2 -> GetYaxis()->SetRangeUser(0, 10000);
                        h2 -> SetTitle(Form(";energy_{REF} [a.u.]; time_{L} - time_{REF} [ps]"));
						h2->Draw("colz");

						latex = new TLatex(0.40,0.8,Form("#splitline{bar %02dL}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kRed);
					    latex -> Draw("same");

						c->Print(Form("%s/timeCorrelation/vs_eREF/c_tL_wREFvsEL_bar%02dL_%s_D.png",plotDir.c_str(),iBar,stepLabel.c_str()));
						c->Print(Form("%s/timeCorrelation/vs_eREF/c_tL_wREFvsEL_bar%02dL_%s_D.pdf",plotDir.c_str(),iBar,stepLabel.c_str()));
						delete c;

                        // --- tR-REF vs eREF
						c = new TCanvas(Form("c_tR_REFvseREF_bar%02dR_%s_D",iBar,stepLabel.c_str()),Form("c_tR_REFvseREF_bar%02dR_%s_D",iBar,stepLabel.c_str()));
						h2 = h2_tR_REFvsE_REF[index2];
						h2 -> GetYaxis()->SetRangeUser(0, 10000);
                        h2 -> SetTitle(Form(";energy_{REF} [a.u.]; time_{R} - time_{REF} [ps]"));
						h2->Draw("colz");

						latex = new TLatex(0.40,0.8,Form("#splitline{bar %02dR}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kRed);
					    latex -> Draw("same");

						c->Print(Form("%s/timeCorrelation/vs_eREF/c_tR_wREFvsER_bar%02dR_%s_D.png",plotDir.c_str(),iBar,stepLabel.c_str()));
						c->Print(Form("%s/timeCorrelation/vs_eREF/c_tR_wREFvsER_bar%02dR_%s_D.pdf",plotDir.c_str(),iBar,stepLabel.c_str()));
						delete c;
					}
				}
			}
        }
    }
  
  
  
  
    //------------------------
    //--- 4th loop over events
    int accepted4=0;
    gStyle->SetOptFit(1111);
    for(auto mapIt : trees){
        ModuleEventClass* anEvent = new ModuleEventClass();
        mapIt.second -> SetBranchAddress("event",&anEvent);
        int nEntries = mapIt.second->GetEntries();
        for(int entry = 0; entry < nEntries; ++entry){

	        if( entry%100000 == 0 ) std::cout << ">>> 4th loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
	        
			mapIt.second -> GetEntry(entry);
	  
	        bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
            if (!barFound) continue;
	  
	        int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
	    	int index3( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + 99 );
	  
	        if( !accept[index1][entry] ) continue;
	  
	        int energyBinAverage=0;
			if(anEvent->nClusters==1) {energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;	}
			else if(anEvent->nClusters==2 && anEvent->energyL_post>0 && anEvent->energyR_post>0){energyBinAverage = FindBin(anEvent->energySum,ranges_doubleHits[index1])+1;}
			else if(anEvent->nClusters==2 && anEvent->energyL_post<0 && anEvent->energyR_post<0){
				energyBinAverage = -1;
				continue; // escludo eventi "PRE"
			}
			else if(anEvent->nClusters==3){	energyBinAverage=1;}

	        double  index2((100000000*anEvent->nClusters)+10000000*energyBinAverage+index1 );     
	        long long deltaT = anEvent->timeR - anEvent->timeL;
	  
	        float t1fineMean = 0.5 * ( anEvent->t1fineR + anEvent->t1fineL );
	  
	        if( !fitFunc_energyRatioCorr[index2] ) continue;
	  
	        float energyRatioCorr = fitFunc_energyRatioCorr[index2]->Eval(anEvent->energyR/anEvent->energyL) - fitFunc_energyRatioCorr[index2]->Eval(fitFunc_energyRatio[index2]->GetParameter(1));
	  
	        if( !fitFunc_totRatioCorr[index2] ) continue;
	  
	        float totRatioCorr = fitFunc_totRatioCorr[index2]->Eval(anEvent->totR/anEvent->totL) - fitFunc_totRatioCorr[index2]->Eval(fitFunc_totRatio[index2]->GetParameter(1));
	  
	  
			if( h1_deltaT_energyRatioCorr[index2] == NULL ){
				std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.02f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
				h1_deltaT_energyRatioCorr[index2] = new TH1F(Form("h1_deltaT_energyRatioCorr_%s",labelLR_energyBin.c_str()),"",2000,-12000.,12000.);
				h1_deltaT_totRatioCorr[index2]    = new TH1F(Form("h1_deltaT_totRatioCorr_%s",labelLR_energyBin.c_str()),"",2000,-12000.,12000.);
				p1_deltaT_totRatioCorr_vs_totRatio[index2] = new TProfile(Form("p1_deltaT_totRatioCorr_vs_totRatio_%s",labelLR_energyBin.c_str()),"",50,fitFunc_totRatio[index2]->GetParameter(1)-3.*fitFunc_totRatio[index2]->GetParameter(2), fitFunc_totRatio[index2]->GetParameter(1)+3.*fitFunc_totRatio[index2]->GetParameter(2));
				h2_deltaT_totRatioCorr_vs_totRatio[index2] = new TH2F(Form("h2_deltaT_totRatioCorr_vs_totRatio_%s",labelLR_energyBin.c_str()),"",50,fitFunc_totRatio[index2]->GetParameter(1)-3.*fitFunc_totRatio[index2]->GetParameter(2), fitFunc_totRatio[index2]->GetParameter(1)+3.*fitFunc_totRatio[index2]->GetParameter(2), 2000, -12000., 12000. );
				p1_deltaT_energyRatioCorr_vs_totRatio[index2] = new TProfile(Form("p1_deltaT_energyRatioCorr_vs_totRatio_%s",labelLR_energyBin.c_str()),"",50,fitFunc_totRatio[index2]->GetParameter(1)-3.*fitFunc_totRatio[index2]->GetParameter(2), fitFunc_totRatio[index2]->GetParameter(1)+3.*fitFunc_totRatio[index2]->GetParameter(2));
				h2_deltaT_energyRatioCorr_vs_totRatio[index2] = new TH2F(Form("h2_deltaT_energyRatioCorr_vs_totRatio_%s",labelLR_energyBin.c_str()),"",50,fitFunc_totRatio[index2]->GetParameter(1)-3.*fitFunc_totRatio[index2]->GetParameter(2), fitFunc_totRatio[index2]->GetParameter(1)+3.*fitFunc_totRatio[index2]->GetParameter(2), 2000, -12000., 12000. );
				p1_deltaT_energyRatioCorr_vs_t1fineMean[index2] = new TProfile(Form("p1_deltaT_energyRatioCorr_vs_t1fineMean_%s",labelLR_energyBin.c_str()),"",50,0,1000);
				p1_deltaT_totRatioCorr_vs_t1fineMean[index2] = new TProfile(Form("p1_deltaT_totRatioCorr_vs_t1fineMean_%s",labelLR_energyBin.c_str()),"",50,0,1000);
				p1_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean[index2] = new TProfile(Form("p1_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean_%s",labelLR_energyBin.c_str()),"",50,0,1000);  
				h2_deltaT_energyRatioCorr_vs_t1fineMean[index2] = new TH2F(Form("h2_deltaT_energyRatioCorr_vs_t1fineMean_%s",labelLR_energyBin.c_str()),"",50,0,1000,2000, -12000., 12000.);
				h2_deltaT_totRatioCorr_vs_t1fineMean[index2] = new TH2F(Form("h2_deltaT_totRatioCorr_vs_t1fineMean_%s",labelLR_energyBin.c_str()),"",50,0,1000, 2000, -12000., 12000.);
				h2_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean[index2] = new TH2F(Form("h2_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean_%s",labelLR_energyBin.c_str()),"",50,0,1000, 2000, -12000., 12000.);
			}
			if(h1_deltaT_corr_REF[index2]==NULL){
				std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
                h1_deltaT_corr_REF[index2] = new TH1F(Form("h1_deltaT_corr_REF_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000);
                h1_deltaT_LR_corr_REF[index2] = new TH1F(Form("h1_deltaT_LR_corr_REF_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000);
                h1_deltaT_LR_corr_REF_w[index2] = new TH1F(Form("h1_deltaT_LR_corr_REF_w_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000);
			}
			if(h1_deltaT_eneryRatioCorr_REF[index3]==NULL){
				h1_deltaT_eneryRatioCorr_REF[index3] = new TH1F(Form("h1_deltaT_energyRatioCorr_externalBar_Vov%.2f_th%02d",anEvent->Vov,anEvent->vth1),"",2000,-12000,12000);

				h2_deltaT_energyRatioCorr_vs_t1fineMean_REF[index3] = new TH2F(Form("h2_deltaT_energyRatioCorr_vs_t1fineMean_REF_Vov%.2f_th%02d",anEvent->Vov,anEvent->vth1),"",100,0,1000,2000,-12000,12000);
				p1_deltaT_energyRatioCorr_vs_t1fineMean_REF[index3] = new TProfile(Form("p1_deltaT_energyRatioCorr_vs_t1fineMean_REF_Vov%.2f_th%02d",anEvent->Vov,anEvent->vth1),"",50,0,1000);
			}
	  
            float enRatioCorr_REF = fitFunc_energyRatioCorrection_REF[index3]->Eval(anEvent->energyR_ext/anEvent->energyL_ext) - fitFunc_energyRatioCorrection_REF[index3]->Eval(fitFunc_energyRatio_REF[index3]->GetParameter(1)); 
			long long deltaT_REF = anEvent->timeL_ext - anEvent->timeR_ext;
			float t1fineMean_REF = 0.5 * ( anEvent->t1fineR_ext + anEvent->t1fineL_ext );

            if(fabs(deltaT_REF - enRatioCorr_REF)<10000 ){
				float enMin = fitFunc_energyRatio_REF[index3]->GetParameter(1) - 1.*fitFunc_energyRatio_REF[index3]->GetParameter(2);
				float enMax = fitFunc_energyRatio_REF[index3]->GetParameter(1) + 1.*fitFunc_energyRatio_REF[index3]->GetParameter(2);
				if(anEvent->energyR_ext/anEvent->energyL_ext > enMin && anEvent->energyR_ext/anEvent->energyL_ext < enMax){
					h1_deltaT_eneryRatioCorr_REF[index3] -> Fill( deltaT_REF - enRatioCorr_REF );
				    h2_deltaT_energyRatioCorr_vs_t1fineMean_REF[index3] -> Fill( t1fineMean_REF, deltaT_REF - enRatioCorr_REF );
				    p1_deltaT_energyRatioCorr_vs_t1fineMean_REF[index3] -> Fill( t1fineMean_REF, deltaT_REF - enRatioCorr_REF );
				}	
			}

			if (fabs(deltaT - energyRatioCorr)<10000 ) {
				accepted4++;
				h1_deltaT_energyRatioCorr[index2] -> Fill( deltaT  - energyRatioCorr );
				//if (fabs(deltaT - energyRatioCorr - h1_deltaT[index2]->GetMean() )< 3*h1_deltaT[index2]->GetRMS() ){
				p1_deltaT_energyRatioCorr_vs_t1fineMean[index2] -> Fill( t1fineMean, deltaT - energyRatioCorr );
				h2_deltaT_energyRatioCorr_vs_t1fineMean[index2] -> Fill( t1fineMean, deltaT - energyRatioCorr );
				p1_deltaT_energyRatioCorr_vs_totRatio[index2] -> Fill( anEvent->totR/anEvent->totL, deltaT - energyRatioCorr );
				h2_deltaT_energyRatioCorr_vs_totRatio[index2] -> Fill( anEvent->totR/anEvent->totL, deltaT - energyRatioCorr );
				//}
			}
	  
			if (fabs(deltaT - totRatioCorr)<10000 ) {
				h1_deltaT_totRatioCorr[index2] -> Fill( deltaT  - totRatioCorr );
				h2_deltaT_totRatioCorr_vs_totRatio[index2] -> Fill( anEvent->totR/anEvent->totL, deltaT  - totRatioCorr );
				//if (fabs(deltaT - totRatioCorr - h1_deltaT[index2]->GetMean()) < 3.*h1_deltaT[index2]->GetRMS() ){ 
				p1_deltaT_totRatioCorr_vs_t1fineMean[index2] -> Fill( t1fineMean, deltaT - totRatioCorr );
				h2_deltaT_totRatioCorr_vs_t1fineMean[index2] -> Fill( t1fineMean, deltaT - totRatioCorr );   
				//}
			}
            
			if(anEvent->nClusters==2 && anEvent->energyL_post>0 && anEvent->energyR_post>0){
				if(anEvent->barID==15) continue;
				double E_mean=0.5*(anEvent->energyL+anEvent->energyR);
				double E_mean_post=0.5*(anEvent->energyL_post+anEvent->energyR_post);
				double t_mean_w=((static_cast<double>(anEvent->energyL)*static_cast<double>(anEvent->timeL))+(static_cast<double>(anEvent->energyR)*static_cast<double>(anEvent->timeR)))/(static_cast<double>(anEvent->energyL)+static_cast<double>(anEvent->energyR));// weighted mean
				double t_mean_post_w=((static_cast<double>(anEvent->energyL_post)*static_cast<double>(anEvent->timeL_post))+(static_cast<double>(anEvent->energyR_post)*static_cast<double>(anEvent->timeR_post)))/(static_cast<double>(anEvent->energyL_post)+static_cast<double>(anEvent->energyR_post));
                if(!fitFunc_energyRatioCorr_bars_w[index2]) continue;
				double energyRatioCorr_bars_w = fitFunc_energyRatioCorr_bars_w[index2]->Eval(E_mean/E_mean_post) - fitFunc_energyRatioCorr_bars_w[index2]->Eval(1.);
				if (fabs((t_mean_w-t_mean_post_w) - energyRatioCorr_bars_w)<10000 ){
					float my_timeLow = my_CTRMeans[index2] - 3.* my_CTRSigmas[index2];
	                float my_timeHig = my_CTRMeans[index2] + 3.* my_CTRSigmas[index2];
					if( ( t_mean_w-t_mean_post_w > my_timeLow ) && ( t_mean_w-t_mean_post_w < my_timeHig ) ){
						h1_deltaT_bars_w[index2]->Fill(t_mean_w - t_mean_post_w - energyRatioCorr_bars_w);
					}
				}

			}

			if(anEvent->nClusters==2 && anEvent->energyL_post>0 && anEvent->energyR_post>0 && anEvent->energyL_ext>0 && anEvent->energyR_ext>0){
                double tL_corr, tL_post_corr, tR_corr, tR_post_corr;
			    double time_ave_REF = 0.5*(static_cast<double>(anEvent->timeL_ext)+static_cast<double>(anEvent->timeR_ext));
				double time_ave_REF_w = (static_cast<double>(anEvent->timeL_ext)*static_cast<double>(anEvent->energyL_ext)+static_cast<double>(anEvent->timeR_ext)*static_cast<double>(anEvent->energyR_ext))/(static_cast<double>(anEvent->energyL_ext)+static_cast<double>(anEvent->energyR_ext));
			    double xMax;
			    double corr;
			    double t_ave, t_ave_post;
			    if((anEvent->timeL-time_ave_REF) > 3500 &&(anEvent->timeL_post-time_ave_REF) > 3500 && (anEvent->timeR-time_ave_REF) > 3500 && (anEvent->timeR_post-time_ave_REF) > 3500){
				    if(anEvent->barID==15) continue;
					if(!fitFunc_timeL[index2+1] || !fitFunc_timeR[index2+1]) continue;
				    xMax = xMax_h2_tL_REFvseL[index2];
				    corr = fitFunc_timeL[index2]->Eval(anEvent->energyL) - fitFunc_timeL[index2]->Eval(xMax);
				    tL_corr = static_cast<double>(anEvent->timeL) - corr;

					xMax = xMax_h2_tR_REFvseR[index2];
				    corr = fitFunc_timeR[index2]->Eval(anEvent->energyR) - fitFunc_timeR[index2]->Eval(xMax);
				    tR_corr = static_cast<double>(anEvent->timeR) - corr;

					xMax = xMax_h2_tL_REFvseL[index2+1];
				    corr = fitFunc_timeL[index2+1]->Eval(anEvent->energyL_post) - fitFunc_timeL[index2+1]->Eval(xMax);
				    tL_post_corr =static_cast<double>(anEvent->timeL_post) - corr;

					xMax = xMax_h2_tR_REFvseR[index2+1];
				    corr = fitFunc_timeR[index2+1]->Eval(anEvent->energyR_post) - fitFunc_timeR[index2+1]->Eval(xMax);
				    tR_post_corr =static_cast<double>(anEvent->timeR_post) - corr;

                    t_ave = 0.5*(tL_corr + tR_corr);
				    t_ave_post = 0.5*(tL_post_corr + tR_post_corr);
				    double energy = 0.5*(anEvent->energyL + anEvent->energyR);
				    double energy_post = 0.5*(anEvent->energyL_post + anEvent->energyR_post);
                    double tL_ave = (static_cast<double>(anEvent->energyL)*tL_corr + static_cast<double>(anEvent->energyL_post)*tL_post_corr)/(static_cast<double>(anEvent->energyL) + static_cast<double>(anEvent->energyL_post));
                    double tR_ave = (static_cast<double>(anEvent->energyR)*tR_corr + static_cast<double>(anEvent->energyR_post)*tR_post_corr)/(static_cast<double>(anEvent->energyR) + static_cast<double>(anEvent->energyR_post));
                    
					if((anEvent->energyL>80 && anEvent->energyR>80 && anEvent->energyL_post>80 && anEvent->energyR_post>80) &&
					   (anEvent->energyL<800 && anEvent->energyR<600 && anEvent->energyL_post<800 && anEvent->energyR_post<600)){
			           
						h1_deltaT_corr_REF[index2] -> Fill((energy*t_ave + energy_post*t_ave_post)/(energy + energy_post) - time_ave_REF);
						
						h1_deltaT_LR_corr_REF[index2] -> Fill(0.5*(tL_ave + tR_ave) - time_ave_REF);
						h1_deltaT_LR_corr_REF_w[index2] -> Fill(0.5*(tL_ave + tR_ave) - time_ave_REF_w);
						/*if(entry%100000==0){
							std::cout<<"\nentry: "<<entry<<"  index2: "<<index2<<std::endl;
							std::cout<<"tL_corr: "<<tL_corr<<"  tR_corr: "<<tR_corr<<std::endl;
							std::cout<<"tL_post_corr: "<<tL_post_corr<<"  tR_post_corr: "<<tR_post_corr<<std::endl;
							std::cout<<"energyL: "<<anEvent->energyL<<"  energyR: "<<anEvent->energyR<<std::endl;
							std::cout<<"energyL_post: "<<anEvent->energyL_post<<"  energyR_post: "<<anEvent->energyR_post<<std::endl;
							std::cout<<"energyL_ref: "<<anEvent->energyL_ext<<"  energyR_ref: "<<anEvent->energyR_ext<<std::endl;
							std::cout<<"tL_ave: "<<tL_ave<<"  tR_ave: "<<tR_ave<<std::endl;
							std::cout<<"time_ave_REF: "<<time_ave_REF<<"  time_ave_REF_w: "<<time_ave_REF_w<<std::endl;
							std::cout<<"h1_deltaT_LR_corr_REF[index2] filled with: "<<0.5*(tL_ave + tR_ave) - time_ave_REF<<std::endl;
						}*/
                        
					}
			    }
		    }
	    }
        std::cout << std::endl;
    }
  
  
    //------------------
    //--- draw 4th plots
	std::cout<<"4th loop accepted events: "<<accepted4<<std::endl;
	std::cout<<"Draw 4th plots ... "<<std::endl;
    for(auto stepLabel : stepLabels) {
        float Vov = map_Vovs[stepLabel];
        float vth1 = map_ths[stepLabel];

		std::string extLabel(Form("externalBar_L-R_%s",stepLabel.c_str())); 
		int index3( (10000*int(Vov*100.)) + (100*vth1) + 99 );
        // -- energy corr deltaT
		c = new TCanvas(Form("c_deltaT_energyRatioCorr_%s",extLabel.c_str()),Form("c_deltaT_energyRatioCorr_%s",extLabel.c_str()));
			  
		histo = h1_deltaT_eneryRatioCorr_REF[index3];
		histo -> SetLineWidth(2);
		histo -> SetLineColor(kBlue+2);
		histo -> SetMarkerColor(kBlue+2);
			  
		TF1* fitFunc = new TF1(Form("fitFunc_energyCorr_%s",extLabel.c_str()),"gaus",-10000, 10000);
		drawDeltaT(c, histo, fitFunc, "energy-corrected", "en.Corr","");
			  
		outFile -> cd();
		histo -> Write();

		c -> Print(Form("%s/externalBar/c_deltaT_energyRatioCorr__%s.pdf",plotDir.c_str(),extLabel.c_str()));
		c -> Print(Form("%s/externalBar/c_deltaT_energyRatioCorr__%s.png",plotDir.c_str(),extLabel.c_str()));
		delete c;

		c = new TCanvas(Form("c_deltaT_energyRatioCorr__vs_t1fineMean_%s",extLabel.c_str()),Form("c_deltaT_energyRatioCorr_vs_t1fineMean_%s",extLabel.c_str()));			  
		h2 = h2_deltaT_energyRatioCorr_vs_t1fineMean_REF[index3];
		h2 -> GetYaxis()->SetRangeUser(h2 -> GetMean(2) - 600., h2 -> GetMean(2)+ 600);
		h2 -> SetTitle(Form(";t1fineMean;#Deltat [ps]"));
		h2 -> Draw("colz");
		prof = p1_deltaT_energyRatioCorr_vs_t1fineMean_REF[index3];
		prof -> SetTitle(Form(";t1fineMean;#Deltat [ps]"));
		prof -> Draw("plsame");
			  
		latex = new TLatex(0.40,0.85,Form("#splitline{external bar}{V_{OV} = %.2f V, th. = %d DAC}",Vov,int(vth1)));
		latex -> SetNDC();
		latex -> SetTextFont(42);
		latex -> SetTextSize(0.04);
		latex -> SetTextColor(kRed);
		latex -> Draw("same");
		 
		c -> Print(Form("%s/externalBar/c_deltaT_energyRatioCorr_vs_t1fineMean__%s.png",plotDir.c_str(),extLabel.c_str()));
		c -> Print(Form("%s/externalBar/c_deltaT_energyRatioCorr_vs_t1fineMean__%s.pdf",plotDir.c_str(),extLabel.c_str()));
		delete c;
		delete latex;

        
		for(int iBar = 0; iBar < 16; ++iBar) {
	        bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
	        if (!barFound) continue;
	  
	        int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );
      	    if( !ranges["L-R"][index1] ) continue;
	  
	        std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
	  
	        int nEnergyBins = ranges["L-R"][index1]->size()-1;
	  
	        for(int i=1;i<=3;i++){//loop over single double triple 
				for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin){
					double  index2( (100000000*i)+10000000*iEnergyBin + index1 );
			  
					if(!h1_deltaT_energyRatioCorr[index2]) continue;
			        
					std::string labelBars(Form("bars%02d-%02d_%s",iBar,iBar+1,stepLabel.c_str()));
					std::string labelLR_energyBin;
					std::string labelBars_energyBin;
					if(i==1){labelLR_energyBin = Form("%s_energyBin%02d_S",labelLR.c_str(),iEnergyBin);}
					else if(i==2){
						labelLR_energyBin = Form("%s_energyBin%02d_D",labelLR.c_str(),iEnergyBin);
						labelBars_energyBin = Form("%s_energyBin%02d_D",labelBars.c_str(),iEnergyBin);
					}
					else if(i==3){labelLR_energyBin = Form("%s_energyBin%02d_T",labelLR.c_str(),iEnergyBin);}
			  
					// -- energy corr deltaT
					c = new TCanvas(Form("c_deltaT_energyRatioCorr_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioCorr_%s",labelLR_energyBin.c_str()));
			  
					histo = h1_deltaT_energyRatioCorr[index2];
					histo -> SetLineWidth(2);
					histo -> SetLineColor(kBlue);
					histo -> SetMarkerColor(kBlue);
			  
					TF1* fitFunc = new TF1(Form("fitFunc_energyCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					drawDeltaT(c, histo, fitFunc, "energy-corrected", "en.Corr","");
			  
					outFile -> cd();
					histo -> Write();
			  
			  
					// -- totRatio corr deltaT
					c2 = new TCanvas(Form("c_deltaT_totRatioCorr_%s",labelLR_energyBin.c_str()),Form("c_deltaT_totRatioCorr_%s",labelLR_energyBin.c_str()));
			  
					histo = h1_deltaT_totRatioCorr[index2];
					histo -> SetLineWidth(2);
					histo -> SetLineColor(kBlue);
					histo -> SetMarkerColor(kBlue);
			  
					fitFunc = new TF1(Form("fitFunc_totCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					drawDeltaT(c2, histo, fitFunc, "ToT-corrected", "totCorr", "");
			  
					outFile -> cd();
					histo -> Write();
			  
			  
					// -- raw delta T
					histo = h1_deltaT[index2];
					histo -> SetLineWidth(2);
					histo -> SetLineColor(kRed);
					histo -> SetMarkerColor(kRed);
			  
					fitFunc = new TF1(Form("fitFunc_%s",labelLR_energyBin.c_str()),"gaus",-10000,10000);
					drawDeltaT(c, histo, fitFunc, "", "raw", "same");
					drawDeltaT(c2, histo, fitFunc, "", "raw", "same");
			  
					outFile -> cd();
					histo -> Write();
			  
					c -> Print(Form("%s/CTR_energyRatioCorr/c_deltaT_energyRatioCorr__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/CTR_energyRatioCorr/c_deltaT_energyRatioCorr__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete c;
			  
					c2 -> Print(Form("%s/CTR_totRatioCorr/c_deltaT_energyRatioCorr__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					c2 -> Print(Form("%s/CTR_totRatioCorr/c_deltaT_energyRatioCorr__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete c2;
			  
			  
			  
					// -- draw deltaT energyRatioCorr vs totRatio
					if(!p1_deltaT_energyRatioCorr_vs_totRatio[index2]) continue;
					c = new TCanvas(Form("c_deltaT_energyRatioCorr_vs_totRatio_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioCorr_vs_totRatio_%s",labelLR_energyBin.c_str()));                  
			  
					prof = p1_deltaT_energyRatioCorr_vs_totRatio[index2];
					prof -> SetTitle(Form(";ToT_{right} / ToT_{left};#Deltat [ps]"));
					prof -> GetYaxis() -> SetRangeUser(CTRMeans[index2]-3.*CTRSigmas[index2],CTRMeans[index2]+3.*CTRSigmas[index2]);
					prof -> Draw("pl");
			  
					latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					latex -> SetNDC();
					latex -> SetTextFont(42);
					latex -> SetTextSize(0.04);
					latex -> SetTextColor(kRed);
					latex -> Draw("same");                                                                                                                                        
			  
					float fitXMin = fitFunc_totRatio[index2]->GetParameter(1) - 3.*fitFunc_totRatio[index2]->GetParameter(2);
					float fitXMax = fitFunc_totRatio[index2]->GetParameter(1) + 3.*fitFunc_totRatio[index2]->GetParameter(2);
					fitFunc_energyRatioCorr_totRatioCorr[index2] = new TF1(Form("fitFunc_energyRatioCorr_totRatioCorr_%s",labelLR_energyBin.c_str()),"pol3",fitXMin,fitXMax);
					prof -> Fit(fitFunc_energyRatioCorr_totRatioCorr[index2],"QRS+");
					fitFunc_energyRatioCorr_totRatioCorr[index2] -> SetLineColor(kRed);
					fitFunc_energyRatioCorr_totRatioCorr[index2] -> SetLineWidth(2);
					fitFunc_energyRatioCorr_totRatioCorr[index2] -> Draw("same");
			  
					c -> Print(Form("%s/energyRatioCorr_totRatioCorr/c_deltaT_energyRatioCorr_vs_totRatio__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/energyRatioCorr_totRatioCorr/c_deltaT_energyRatioCorr_vs_totRatio__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete c;                       
			  
			  
					// -- draw deltaT vs t1fine
					if(!p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]) continue;
					c = new TCanvas(Form("c_deltaT_energyRatioCorr__vs_t1fineMean_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioCorr_vs_t1fineMean_%s",labelLR_energyBin.c_str()));
			  
					h2 = h2_deltaT_energyRatioCorr_vs_t1fineMean[index2];
					h2 -> GetYaxis()->SetRangeUser(h2 -> GetMean(2) - 600., h2 -> GetMean(2)+ 600);
					h2 -> SetTitle(Form(";t1fineMean;#Deltat [ps]"));
					h2 -> Draw("colz");
					prof = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2];
					prof -> SetTitle(Form(";t1fineMean;#Deltat [ps]"));
					prof -> Draw("plsame");
			  
					latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					latex -> SetNDC();
					latex -> SetTextFont(42);
					latex -> SetTextSize(0.04);
					latex -> SetTextColor(kRed);
					latex -> Draw("same");
			  
					c -> Print(Form("%s/phaseCorr/c_deltaT_energyRatioCorr_vs_t1fineMean__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/phaseCorr/c_deltaT_energyRatioCorr_vs_t1fineMean__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete c;
					delete latex;
			  
			  
					// -- draw deltaT vs t1fine
					if(!p1_deltaT_totRatioCorr_vs_t1fineMean[index2]) continue;
					c = new TCanvas(Form("c_deltaT_totRatioCorr__vs_t1fineMean_%s",labelLR_energyBin.c_str()),Form("c_deltaT_totRatioCorr_vs_t1fineMean_%s",labelLR_energyBin.c_str()));
					c -> SetGridy();
			  
					h2 = h2_deltaT_totRatioCorr_vs_t1fineMean[index2]; 
					h2 -> GetYaxis()->SetRangeUser(h2 -> GetMean(2) -600., h2 -> GetMean(2)+600);
					h2 -> SetTitle(Form(";t1fineMean;#Deltat [ps]"));
					h2 -> Draw("colz");
					prof = p1_deltaT_totRatioCorr_vs_t1fineMean[index2];
					//prof -> Draw("pl");
					prof -> Draw("plsame");
			  
					latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					latex -> SetNDC();
					latex -> SetTextFont(42);
					latex -> SetTextSize(0.04);
					latex -> SetTextColor(kRed);
					latex -> Draw("same");
			  
					c -> Print(Form("%s/phaseCorr/c_deltaT_totRatioCorr_vs_t1fineMean__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/phaseCorr/c_deltaT_totRatioCorr_vs_t1fineMean__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete c;
					delete latex; 
					if(i==2){
						// --- draw deltaT from bars difference with weigthed mean : start
						c = new TCanvas(Form("c_deltaT_w_%s",labelBars_energyBin.c_str()),Form("c_deltaT_w_%s",labelBars_energyBin.c_str()));
						histo = h1_deltaT_bars_w[index2];
					    histo -> SetLineWidth(2);
					    histo -> SetLineColor(kGreen);
					    histo -> SetMarkerColor(kGreen);
					    drawDeltaT(c, histo, fitFunc, Form("%02d - %02d diff ",iBar,iBar+1), "enCorr_w","");
					    outFile -> cd();
					    histo -> Write();

						c -> Print(Form("%s/CTR_energyRatioCorr/c_deltaT__w_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
				    	c -> Print(Form("%s/CTR_energyRatioCorr/c_deltaT__w_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
						delete fitFunc;
						delete c;
						// --- draw deltaT from bars difference with weigthed mean : end

					    c = new TCanvas(Form("c_deltaT_REF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_REF_%s",labelLR_energyBin.c_str()));
						histo = h1_deltaT_corr_REF[index2];
					    histo -> SetLineWidth(2);
					    histo -> SetLineColor(kRed+5);
					    histo -> SetMarkerColor(kRed+5);
			            fitFunc = new TF1(Form("fitFunc_deltaT_REF_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					    drawDeltaT(c, histo, fitFunc, Form("[(%02d + %02d) - REF]  ",iBar,iBar+1), "en.Corr","");
					    outFile -> cd();
					    histo -> Write();

						c -> Print(Form("%s/CTR_REF/CTR_REF_1Bar/c_deltaT_REF_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
				    	c -> Print(Form("%s/CTR_REF/CTR_REF_1Bar/c_deltaT_REF_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
						delete c;

						c = new TCanvas(Form("c_deltaT_LR_REF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_LR_REF_%s",labelLR_energyBin.c_str()));
						histo = h1_deltaT_LR_corr_REF[index2];
					    histo -> SetLineWidth(2);
					    histo -> SetLineColor(kGreen+3);
					    histo -> SetMarkerColor(kGreen+3);
			            fitFunc = new TF1(Form("fitFunc_deltaT_LR_REF_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					    drawDeltaT(c, histo, fitFunc, Form("[(%02d + %02d) - REF]  ",iBar,iBar+1), "en.Corr","");
					    outFile -> cd();
					    histo -> Write();

						c -> Print(Form("%s/CTR_REF/CTR_REF_2LR/c_deltaT_LRcorr_REF_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
				    	c -> Print(Form("%s/CTR_REF/CTR_REF_2LR/c_deltaT_LRcorr_REF_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
						delete c;

						c = new TCanvas(Form("c_deltaT_LR_REFw_%s",labelLR_energyBin.c_str()),Form("c_deltaT_LR_REFw_%s",labelLR_energyBin.c_str()));
						histo = h1_deltaT_LR_corr_REF_w[index2];
					    histo -> SetLineWidth(2);
					    histo -> SetLineColor(kGreen+2);
					    histo -> SetMarkerColor(kGreen+2);
			            fitFunc = new TF1(Form("fitFunc_deltaT_LR_REFw_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					    drawDeltaT(c, histo, fitFunc, Form("[(%02d + %02d) - REF]  ",iBar,iBar+1), "en.Corr","");
					    outFile -> cd();
					    histo -> Write();

						c -> Print(Form("%s/CTR_REF/CTR_REF_3LR/c_deltaT_LRcorr_REFw_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
				    	c -> Print(Form("%s/CTR_REF/CTR_REF_3LR/c_deltaT_LRcorr_REFw_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
						delete c;
					}
				}
			}
	    }
    }
  
  
  
  
    //------------------------
    //--- 5th loop over events
	int accepted5=0;
    for(auto mapIt : trees) {
        ModuleEventClass* anEvent = new ModuleEventClass();
        mapIt.second -> SetBranchAddress("event",&anEvent);
      
        int nEntries = mapIt.second->GetEntries();

        for(int entry = 0; entry < nEntries; ++entry) {
	        if( entry%100000 == 0 ) std::cout << ">>> 5th loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
	        mapIt.second -> GetEntry(entry);
	  
	        bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;                                                                             
	        if (!barFound) continue;
	        
	        int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID ); 
		    int index3( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + 99 );

	        if( !accept[index1][entry] ) continue;
	  
	        int energyBinAverage=0;
			if(anEvent->nClusters==1) {energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;	}
			else if(anEvent->nClusters==2 && anEvent->energyL_post>0 && anEvent->energyR_post>0){energyBinAverage = FindBin(anEvent->energySum,ranges_doubleHits[index1])+1;}
			else if(anEvent->nClusters==2 && anEvent->energyL_post<0 && anEvent->energyR_post<0){
				energyBinAverage = -1;
				continue; // escludo eventi "PRE"
			}
			else if(anEvent->nClusters==3){	energyBinAverage=1;}

	        double  index2( (100000000*anEvent->nClusters)+10000000*energyBinAverage+index1 );     
	  
	        long long deltaT = anEvent->timeR - anEvent->timeL;
	  
	        if( !fitFunc_energyRatioCorr[index2] )	continue;
	        if( !p1_deltaT_energyRatioCorr_vs_t1fineMean[index2] )	continue;
	  
	        float energyRatioCorr = fitFunc_energyRatioCorr[index2]->Eval(anEvent->energyR/anEvent->energyL) -
	        fitFunc_energyRatioCorr[index2]->Eval(fitFunc_energyRatio[index2]->GetParameter(1));
	 
	        if( !fitFunc_energyRatioCorr_totRatioCorr[index2] )	continue;
	        float energyRatioCorr_totRatioCorr = fitFunc_energyRatioCorr_totRatioCorr[index2]->Eval(anEvent->totR/anEvent->totL) - fitFunc_energyRatioCorr_totRatioCorr[index2]->Eval(fitFunc_totRatio[index2]->GetParameter(1));                                 
	  
	        float t1fineMean = 0.5* ( anEvent->t1fineR + anEvent->t1fineL );
	        int t1fineBin = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]->FindBin(t1fineMean);
	        int t1fineBin2 = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]->FindBin(h1_t1fineMean[index2]->GetMean());
	        float t1fineCorr = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]->GetBinContent(t1fineBin) - p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]->GetBinContent( t1fineBin2 );

			float deltaT_corr = deltaT  - energyRatioCorr - t1fineCorr;
            
			float t1fineMean_REF = 0.5* ( anEvent->t1fineR_ext + anEvent->t1fineL_ext );
	        int t1fineBin_REF = p1_deltaT_energyRatioCorr_vs_t1fineMean_REF[index3]->FindBin(t1fineMean_REF);
	        int t1fineBin2_REF = p1_deltaT_energyRatioCorr_vs_t1fineMean_REF[index3]->FindBin(h1_t1fineMean_REF[index3]->GetMean());
	        float t1fineCorr_REF = p1_deltaT_energyRatioCorr_vs_t1fineMean_REF[index3]->GetBinContent(t1fineBin_REF) - p1_deltaT_energyRatioCorr_vs_t1fineMean_REF[index3]->GetBinContent( t1fineBin2_REF );
			float enRatioCorr_REF = fitFunc_energyRatioCorrection_REF[index3]->Eval(anEvent->energyR_ext/anEvent->energyL_ext) - fitFunc_energyRatioCorrection_REF[index3]->Eval(fitFunc_energyRatio_REF[index3]->GetParameter(1)); 
			long long deltaT_REF = anEvent->timeL_ext - anEvent->timeR_ext; 
			float deltaT_corr_REF = deltaT_REF - enRatioCorr_REF - t1fineCorr_REF;
              
			long long deltaT_post=0;
			float energyRatioCorr_post=0;
			float t1fineMean_post=0;
			int t1fineBin_post=0;
			int t1fineBin2_post=0;
			float t1fineCorr_post=0;
			float deltaT_post_corr=0;
			if(anEvent->nClusters==2 && anEvent->energyL_post>0 && anEvent->energyR_post>0){
				deltaT_post = anEvent->timeR_post - anEvent->timeL_post;
				if( !fitFunc_energyRatioCorr[index2+1] )	continue;
	            if( !p1_deltaT_energyRatioCorr_vs_t1fineMean[index2+1] )continue;
				energyRatioCorr_post = fitFunc_energyRatioCorr[index2+1]->Eval(anEvent->energyR_post/anEvent->energyL_post) - fitFunc_energyRatioCorr[index2+1]->Eval(fitFunc_energyRatio[index2+1]->GetParameter(1));
				t1fineMean_post = 0.5* ( anEvent->t1fineR_post + anEvent->t1fineL_post );
	            t1fineBin_post = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2+1]->FindBin(t1fineMean_post);
	            t1fineBin2_post = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2+1]->FindBin(h1_t1fineMean[index2+1]->GetMean());
	            t1fineCorr_post = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2+1]->GetBinContent(t1fineBin_post) - p1_deltaT_energyRatioCorr_vs_t1fineMean[index2+1]->GetBinContent( t1fineBin2_post );
				deltaT_post_corr = deltaT_post - energyRatioCorr_post - t1fineCorr_post;
			}

	  
	        if( h1_deltaT_energyRatioPhaseCorr[index2] == NULL ) {
	            std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
				std::string labelLR_energyBin_bars(Form("bars%02d-%02d_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->barID+1,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));

	            h1_deltaT_energyRatioPhaseCorr[index2] = new TH1F(Form("h1_deltaT_energyRatioPhaseCorr_%s",labelLR_energyBin.c_str()),"",2000,-12000.,12000.);
				h1_deltaT_energyRatioPhaseCorr_bars[index2] = new TH1F(Form("h1_deltaT_energyRatioPhaseCorr_%s",labelLR_energyBin_bars.c_str()),"",2000,-12000.,12000.);
	            h1_deltaT_energyRatioCorr_totRatioCorr[index2] = new TH1F(Form("h1_deltaT_energyRatioCorr_totRatioCorr_%s",labelLR_energyBin.c_str()),"",2000,-12000.,12000.);
	            p1_deltaT_energyRatioCorr_vs_posX[index2] = new TProfile(Form("p1_deltaT_energyRatioCorr_vs_posX_%s",labelLR_energyBin.c_str()),"",100,-50,50);
	        }

            if(h2_deltaTave_vs_energyRatio[index2] == NULL){
			    std::string label(Form("bars%02d-%02d_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->barID+1,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
				p1_deltaTave_vs_energyRatio[index2] = new TProfile(Form("p1_deltaT_vs_energyRatio_ave_%s",label.c_str()),"",80,0,5.);
				h2_deltaTave_vs_energyRatio[index2] = new TH2F(Form("h2_deltaT_vs_energyRatio_ave_%s",label.c_str()),"",100,0.,5.,2000,-24000.,24000.);
			}
			
	        // -- energyRatio corr
	        if (fabs(deltaT - energyRatioCorr)<10000 ){
				accepted5++;
	            h1_deltaT_energyRatioPhaseCorr[index2] -> Fill( deltaT  - energyRatioCorr - t1fineCorr );
	            h1_deltaT_energyRatioCorr_totRatioCorr[index2] -> Fill( deltaT  - energyRatioCorr - energyRatioCorr_totRatioCorr );
	            p1_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean[index2] -> Fill( t1fineMean, deltaT - energyRatioCorr - energyRatioCorr_totRatioCorr );  
	            h2_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean[index2] -> Fill( t1fineMean, deltaT - energyRatioCorr - energyRatioCorr_totRatioCorr );  
	            if (useTrackInfo && anEvent->nhits>0 && anEvent->x>-100) {
	                p1_deltaT_energyRatioCorr_vs_posX[index2] ->Fill( anEvent->x, deltaT  - energyRatioCorr - t1fineCorr);
	            }
	        }
			if(h1_deltaT_eneryRatioCorr_pasheCorr_REF[index3]==NULL){
				h1_deltaT_eneryRatioCorr_pasheCorr_REF[index3] = new TH1F(Form("h1_deltaT_eneryRatioCorr_pasheCorr_REF_Vov%.2f_th%02d",anEvent->Vov,anEvent->vth1),"",2000,-12000.,12000.);
			}

			if(fabs(deltaT_REF - enRatioCorr_REF)<10000 && anEvent->energyL_ext>0 && anEvent->energyR_ext>0){
                float enMin = fitFunc_energyRatio_REF[index3]->GetParameter(1) - 1.*fitFunc_energyRatio_REF[index3]->GetParameter(2);
                float enMax = fitFunc_energyRatio_REF[index3]->GetParameter(1) + 1.*fitFunc_energyRatio_REF[index3]->GetParameter(2);
                if(anEvent->energyR_ext/anEvent->energyL_ext > enMin && anEvent->energyR_ext/anEvent->energyL_ext < enMax){
                    h1_deltaT_eneryRatioCorr_pasheCorr_REF[index3] ->Fill(deltaT_corr_REF);
                }	
            }


			if(anEvent->nClusters==2 && anEvent->energyL_post>0 && anEvent->energyR_post>0){
				    float E_mean = 0.5*(anEvent->energyL+anEvent->energyR);
					float E_mean_post = 0.5 * (anEvent->energyL_post + anEvent->energyR_post);
					h2_deltaTave_vs_energyRatio[index2] -> Fill((anEvent->energyL+anEvent->energyR)/(anEvent->energyL_post+anEvent->energyR_post), 0.5*(deltaT_corr + deltaT_post_corr));
				    p1_deltaTave_vs_energyRatio[index2] -> Fill((anEvent->energyL+anEvent->energyR)/(anEvent->energyL_post+anEvent->energyR_post), 0.5*(deltaT_corr + deltaT_post_corr));
					
					h1_deltaT_energyRatioPhaseCorr_bars[index2] -> Fill(0.5*(deltaT_corr + deltaT_post_corr));   
			}
	  
	  
	        // -- totRatio corr
	        if( !fitFunc_totRatioCorr[index2] )	continue;
	        if( !p1_deltaT_totRatioCorr_vs_t1fineMean[index2] )	continue;
	  
	        float totRatioCorr = fitFunc_totRatioCorr[index2]->Eval(anEvent->totR/anEvent->totL) - fitFunc_totRatioCorr[index2]->Eval(fitFunc_totRatio[index2]->GetParameter(1));
	  
	        t1fineBin = p1_deltaT_totRatioCorr_vs_t1fineMean[index2]->FindBin(t1fineMean);
	        t1fineBin2 = p1_deltaT_totRatioCorr_vs_t1fineMean[index2]->FindBin(h1_t1fineMean[index2]->GetMean());
	        t1fineCorr = p1_deltaT_totRatioCorr_vs_t1fineMean[index2]->GetBinContent(t1fineBin) - p1_deltaT_totRatioCorr_vs_t1fineMean[index2]->GetBinContent( t1fineBin2 );
	  
	        if( h1_deltaT_totRatioPhaseCorr[index2] == NULL ) {
	            std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
	      
	            h1_deltaT_totRatioPhaseCorr[index2] = new TH1F(Form("h1_deltaT_totRatioPhaseCorr_%s",labelLR_energyBin.c_str()),"",2000,-12000.,12000.);
	            p1_deltaT_totRatioCorr_vs_posX[index2] = new TProfile(Form("p1_deltaT_totRatioCorr_vs_posX_%s",labelLR_energyBin.c_str()),"",100,-50,50);
	      
	            p1_deltaT_totRatioPhaseCorr_vs_totRatio[index2] = new TProfile(Form("p1_deltaT_totRatioPhaseCorr_vs_totRatio_%s",labelLR_energyBin.c_str()),"",50,fitFunc_totRatio[index2]->GetParameter(1)-3.*fitFunc_totRatio[index2]->GetParameter(2), fitFunc_totRatio[index2]->GetParameter(1)+3.*fitFunc_totRatio[index2]->GetParameter(2));
	            h2_deltaT_totRatioPhaseCorr_vs_totRatio[index2] = new TH2F(Form("h2_deltaT_totPhaseRatioCorr_vs_totRatio_%s",labelLR_energyBin.c_str()),"",50,fitFunc_totRatio[index2]->GetParameter(1)-3.*fitFunc_totRatio[index2]->GetParameter(2), fitFunc_totRatio[index2]->GetParameter(1)+3.*fitFunc_totRatio[index2]->GetParameter(2), 2000, -12000., 12000. );
	        }
	  
	        //if (fabs(deltaT - totRatioCorr)>10000 ) continue;
	        if (fabs(deltaT - totRatioCorr)<10000 ) {
	            h1_deltaT_totRatioPhaseCorr[index2] -> Fill( deltaT  - totRatioCorr - t1fineCorr );
	            p1_deltaT_totRatioPhaseCorr_vs_totRatio[index2] -> Fill( anEvent->totR/anEvent->totL, deltaT  - totRatioCorr - t1fineCorr ); 
	            h2_deltaT_totRatioPhaseCorr_vs_totRatio[index2] -> Fill( anEvent->totR/anEvent->totL, deltaT  - totRatioCorr - t1fineCorr ); 
	            if (useTrackInfo && anEvent->nhits>0 && anEvent->x>-100) p1_deltaT_totRatioCorr_vs_posX[index2] ->Fill( anEvent->x, deltaT  - totRatioCorr - t1fineCorr);
	        }
	    }
        std::cout << std::endl;
    }
  
  
    //------------------
    //--- draw 5th plots
    std::map<double,TF1*> fitFunc1_posCorr;
    std::map<double,TF1*> fitFunc2_posCorr;
	std::map<double,TF1*> fitFunc_energyRatioCorr_bars_ave;
	std::map<double,TF1*> fitFunc_delatT_bars;
	std::cout<<"5th loop accepted events: "<<accepted5<<std::endl;
    std::cout<<"Drawing 5th plots ... "<<std::endl;
    for(auto stepLabel : stepLabels) {
        float Vov = map_Vovs[stepLabel];
        float vth1 = map_ths[stepLabel];

        std::string extLabel(Form("externalBar_L-R_%s",stepLabel.c_str())); 
		int index3( (10000*int(Vov*100.)) + (100*vth1) + 99 );

		// -- energy corr phase corr deltaT
		c = new TCanvas(Form("c_deltaT_energyRatioPhaseCorr_%s",extLabel.c_str()),Form("c_deltaT_energyRatioPhaseCorr_%s",extLabel.c_str()));
			  
		histo = h1_deltaT_eneryRatioCorr_pasheCorr_REF[index3];
		histo -> SetLineWidth(2);
		histo -> SetLineColor(kGreen+2);
		histo -> SetMarkerColor(kGreen+2);
			  
		TF1* fitFunc = new TF1(Form("fitFunc_energyCorr_%s",extLabel.c_str()),"gaus",-10000, 10000);
		drawDeltaT(c, histo, fitFunc, "phase-corrected", "en.Corr+phaseCorr","");
			  
		outFile -> cd();
		histo -> Write();

		c -> Print(Form("%s/externalBar/c_deltaT_energyRatioPhaseCorr__%s.pdf",plotDir.c_str(),extLabel.c_str()));
		c -> Print(Form("%s/externalBar/c_deltaT_energyRatioPhaseCorr__%s.png",plotDir.c_str(),extLabel.c_str()));
		delete c;

        for(int iBar = 0; iBar < 16; ++iBar) {
	        bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
	        if (!barFound) continue; 
	  
	        std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
	  
	        int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );
	        if( !ranges["L-R"][index1] ) continue;
	  
	        int nEnergyBins = ranges["L-R"][index1]->size()-1;
	  
	        for(int i=1; i<=3;i++){
				for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin) {
					double  index2( (100000000*i)+10000000*iEnergyBin + index1 );
			  
					if(!h1_deltaT_energyRatioPhaseCorr[index2]) continue;
					std::string labelLR_energyBin;
					if(i==1){labelLR_energyBin = Form("%s_energyBin%02d_S",labelLR.c_str(),iEnergyBin);}
					else if(i==2){labelLR_energyBin = Form("%s_energyBin%02d_D",labelLR.c_str(),iEnergyBin);}
					else if(i==3){labelLR_energyBin = Form("%s_energyBin%02d_T",labelLR.c_str(),iEnergyBin);}
			  

			  
					std::cout << labelLR_energyBin.c_str()<<std::endl;
			  
			  
					// -- energy and phase corr deltaT
					c = new TCanvas(Form("c_deltaT_energyRatioPhaseCorr_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioPhaseCorr_%s",labelLR_energyBin.c_str()));
					histo = h1_deltaT_energyRatioPhaseCorr[index2];
					histo -> SetLineWidth(2);
					histo -> SetLineColor(kGreen+1);
					histo -> SetMarkerColor(kGreen+1);
			  
					TF1* fitFunc = new TF1(Form("fitFunc_energyPhaseCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					drawDeltaT(c, histo, fitFunc, "phase-corrected", "ph.Corr","");
			  
					outFile -> cd();
					histo -> Write();
			  
					// -- energy corr deltaT
					histo = h1_deltaT_energyRatioCorr[index2];
					histo -> SetLineWidth(2);
					histo -> SetLineColor(kBlue);
					histo -> SetMarkerColor(kBlue);
			  
					fitFunc = new TF1(Form("fitFunc_energyCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					drawDeltaT(c, histo, fitFunc, "energy-corrected", "en.Corr", "same");
			  
					c -> Print(Form("%s/CTR_energyRatioCorr_phaseCorr/c_deltaT_energyRatioPhaseCorr__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/CTR_energyRatioCorr_phaseCorr/c_deltaT_energyRatioPhaseCorr__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete c;
            
			  
					// -- tot and phase corr deltaT
					if (!h1_deltaT_totRatioPhaseCorr[index2]) continue;
			  
					c = new TCanvas(Form("c_deltaT_totRatioPhaseCorr_%s",labelLR_energyBin.c_str()),Form("c_deltaT_totRatioPhaseCorr_%s",labelLR_energyBin.c_str()));
			  
					histo = h1_deltaT_totRatioPhaseCorr[index2];
					histo -> SetLineWidth(2);
					histo -> SetLineColor(kGreen+1);
					histo -> SetMarkerColor(kGreen+1);
			  
					fitFunc = new TF1(Form("fitFunc_totPhaseCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					drawDeltaT(c, histo, fitFunc, "phase-corrected", "ph.Corr", "");      	      
			  
					outFile -> cd();
					histo -> Write();
			  
					// -- tot corr deltaT
					histo = h1_deltaT_totRatioCorr[index2];
					histo -> SetLineWidth(2);
					histo -> SetLineColor(kBlue);
					histo -> SetMarkerColor(kBlue);
			  
					fitFunc = new TF1(Form("fitFunc_totCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					drawDeltaT(c, histo, fitFunc, "ToT-corrected", "totCorr", "same"); 
			  
			  
					c -> Print(Form("%s/CTR_totRatioCorr_phaseCorr/c_deltaT_totRatioPhaseCorr__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/CTR_totRatioCorr_phaseCorr/c_deltaT_totRatioPhaseCorr__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete c;
			  
			  
					// -- energy and tot corr deltaT
					if (!h1_deltaT_energyRatioCorr_totRatioCorr[index2]) continue;
					c = new TCanvas(Form("c_deltaT_energyRatioCorr_totRatioCorr_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioCorr_totRatioCorr_%s",labelLR_energyBin.c_str()));                   
					histo = h1_deltaT_energyRatioCorr_totRatioCorr[index2];
					histo -> SetLineWidth(2);
					histo -> SetLineColor(kAzure);
					histo -> SetMarkerColor(kAzure);                                                                                                                            
			  
					fitFunc = new TF1(Form("fitFunc_energyRatioCorr_totRatioCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					drawDeltaT(c, histo, fitFunc, "energy+Tot corrected", "en+ToT.corr", "");
			  
					outFile -> cd();
					histo -> Write();
			  
					// -- energy corrected deltaT
					histo = h1_deltaT_energyRatioCorr[index2];
					histo -> SetLineWidth(2);
					histo -> SetLineColor(kBlue);
					histo -> SetMarkerColor(kBlue);
			  
					fitFunc = new TF1(Form("fitFunc_energyCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					drawDeltaT(c, histo, fitFunc, "energy-corrected", "en.Corr", "same");
			  
					c -> Print(Form("%s/CTR_energyRatioCorr_totRatioCorr/c_deltaT_energyRatioCorr_totRatioCorr__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/CTR_energyRatioCorr_totRatioCorr/c_deltaT_energyRatioCorr_totRatioCorr__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete c;
			  
			  
					// -- draw deltaT (energy+tot corr) vs t1Fine 
					if(!p1_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean[index2]) continue;
			  
					c = new TCanvas(Form("c_deltaT_energyRatioCorr_totRatioCorr__vs_t1fineMean_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean_%s",labelLR_energyBin.c_str()));
					c -> SetGridy();
			  
					h2 = h2_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean[index2];
					h2 -> GetYaxis()->SetRangeUser(h2 -> GetMean(2) -600., h2 -> GetMean(2)+600);
					h2 -> SetTitle(Form(";t1fineMean;#Deltat [ps]"));
					h2 -> Draw("colz");
					prof = p1_deltaT_totRatioCorr_vs_t1fineMean[index2];
					prof -> Draw("plsame");
			  
					latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					latex -> SetNDC();
					latex -> SetTextFont(42);
					latex -> SetTextSize(0.04);
					latex -> SetTextColor(kRed);
					latex -> Draw("same");
			  
					c -> Print(Form("%s/phaseCorr/c_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/phaseCorr/c_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete c;
					delete latex;

					// -- energy and phase corr deltaT between bars
					if(i==2){

						if(!h1_deltaT_energyRatioPhaseCorr_bars[index2]) continue;
						c = new TCanvas(Form("c_deltaT_energyRatioPhaseCorr_bars_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioPhaseCorr_bars_%s",labelLR_energyBin.c_str()));
					    histo = h1_deltaT_energyRatioPhaseCorr_bars[index2];
					    histo -> SetLineWidth(2);
					    histo -> SetLineColor(kGreen+1);
					    histo -> SetMarkerColor(kGreen+1);
			  
					    fitFunc_delatT_bars[index2] = new TF1(Form("fitFunc_energyPhaseCorr_bars_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					    drawDeltaT(c, histo, fitFunc_delatT_bars[index2], Form("%02d - %02d average ",iBar,iBar+1), "enCorr+phaseCorr","");
			  
					    outFile -> cd();
					    histo -> Write();
			        			  
					    c -> Print(Form("%s/CTR_energyRatioCorr_phaseCorr/c_deltaT_bars_energyRatioPhaseCorr__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					    c -> Print(Form("%s/CTR_energyRatioCorr_phaseCorr/c_deltaT_bars_energyRatioPhaseCorr__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					    delete c;

						c = new TCanvas(Form("c_deltaTave_vs_enRatio_%s",labelLR_energyBin.c_str()),Form("c_deltaTave_vs_enRatio_%s",labelLR_energyBin.c_str()));
						if(!p1_deltaTave_vs_energyRatio[index2]) continue;
						h2 = h2_deltaTave_vs_energyRatio[index2];
						h2 -> GetYaxis()->SetRangeUser(-2000, +2000);
					    h2 -> SetTitle(Form(";energy_{%02d} / energy_{%02d};#DeltaT_{average} [ps]",iBar,iBar+1));
					    h2 -> Draw("colz");
					    //outFile -> cd();
					    //h2 -> Write();

						prof = p1_deltaTave_vs_energyRatio[index2];
						//prof -> GetYaxis()->SetRangeUser(-500, +500);
						//prof -> SetTitle(Form(";energy_{%02d} / energy_{%02d};#DeltaT_{average} [ps]",iBar,iBar+1));
						prof -> SetMarkerSize(0.4);
						prof -> Draw("psame");

						latex = new TLatex(0.40,0.85,Form("#splitline{bars %02d-%02d average }{V_{OV} = %.2f V, th. = %d DAC}",iBar,iBar+1,Vov,int(vth1)));
						latex -> SetNDC();
						latex -> SetTextFont(42);
						latex -> SetTextSize(0.04);
						latex -> SetTextColor(kRed);
						latex -> Draw("same");

						c -> Print(Form("%s/energyRatioCorr/c_delatTave_energyRatio__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
						c -> Print(Form("%s/energyRatioCorr/c_delatTave_energyRatio__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
						delete c;

						c = new TCanvas(Form("c_deltaTave_vs_enRatio_p_%s",labelLR_energyBin.c_str()),Form("c_deltaTave_vs_enRatio_p_%s",labelLR_energyBin.c_str()));
						prof -> GetYaxis()->SetRangeUser(-500, +500);
						prof -> SetTitle(Form(";energy_{%02d} / energy_{%02d};#DeltaT_{average} [ps]",iBar,iBar+1));
						prof -> SetMarkerSize(0.4);
						prof -> Draw();
						fitFunc_energyRatioCorr_bars_ave[index2] = new TF1(Form("fitFunc_deltaT_enRatio_%s",labelLR_energyBin.c_str()),"pol3",0.2,5.);
					
						prof -> Fit(fitFunc_energyRatioCorr_bars_ave[index2],"QRS+");
						fitFunc_energyRatioCorr_bars_ave[index2] -> SetLineColor(kRed);
						fitFunc_energyRatioCorr_bars_ave[index2] -> SetLineWidth(2);
						fitFunc_energyRatioCorr_bars_ave[index2] -> Draw("same");
						latex -> Draw("same");

						TLine* line = new TLine(0,fitFunc_delatT_bars[index2]->GetParameter(1),5,fitFunc_delatT_bars[index2]->GetParameter(1));
						line -> SetLineColor(kBlue);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						line -> Draw("same");

						line = new TLine(1,-500,1,500);
						line -> SetLineColor(kBlue);
						line -> SetLineWidth(2);
						line -> SetLineStyle(2);
						line -> Draw("same");

						c -> Print(Form("%s/energyRatioCorr/c_delatTave_energyRatio_p__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
						c -> Print(Form("%s/energyRatioCorr/c_delatTave_energyRatio_p__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
						delete latex;
						delete c;

					}
			  
			  
					// -- draw deltaT vs position
					if (useTrackInfo) {
						c = new TCanvas(Form("c_deltaT_energyRatioCorr_vs_posX_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioCorr_vs_posX_%s",labelLR_energyBin.c_str()));
			
						prof = p1_deltaT_energyRatioCorr_vs_posX[index2];
						prof -> SetTitle(Form("; x [mm] ;#Deltat [ps]"));
						prof -> GetYaxis() -> SetRangeUser(prof->GetMean(2)-3*prof->GetRMS(2), prof->GetMean(2)+3*prof->GetRMS(2));
						prof -> Draw("");
				
						latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
						latex -> SetNDC();
						latex -> SetTextFont(42);
						latex -> SetTextSize(0.04);
						latex -> SetTextColor(kRed);
						latex -> Draw("same");
				
						fitFunc1_posCorr[index2] = new TF1(Form("fitFunc1_posCorr_%s",labelLR_energyBin.c_str()),"pol1",-50,50);
						fitFunc1_posCorr[index2] -> SetRange( prof -> GetMean()-3*prof->GetRMS(), prof -> GetMean()+3*prof->GetRMS());
						prof -> Fit(fitFunc1_posCorr[index2],"QRS+");
						fitFunc1_posCorr[index2] -> SetLineColor(kRed);
						fitFunc1_posCorr[index2] -> SetLineWidth(2);
						fitFunc1_posCorr[index2] -> Draw("same");
				
						c -> Print(Form("%s/positionCorr/c_deltaT_energyRatioCorr_vs_posX__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
						c -> Print(Form("%s/positionCorr/c_deltaT_energyRatioCorr_vs_posX_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
						delete c;
						delete latex;
				
				
						c = new TCanvas(Form("c_deltaT_totRatioCorr_vs_posX_%s",labelLR_energyBin.c_str()),Form("c_deltaT_totRatioCorr_vs_posX_%s",labelLR_energyBin.c_str()));
				
						prof = p1_deltaT_totRatioCorr_vs_posX[index2];
						prof -> SetTitle(Form("; x [mm] ;#Deltat [ps]"));
						prof -> GetYaxis() -> SetRangeUser(prof->GetMean(2)-3*prof->GetRMS(2), prof->GetMean(2)+3*prof->GetRMS(2));
						prof -> Draw("");
				
						latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
						latex -> SetNDC();
						latex -> SetTextFont(42);
						latex -> SetTextSize(0.04);
						latex -> SetTextColor(kRed);
						latex -> Draw("same");
				
						fitFunc2_posCorr[index2] = new TF1(Form("fitFunc2_posCorr_%s",labelLR_energyBin.c_str()),"pol1",-50,50);
						fitFunc2_posCorr[index2] -> SetRange( prof -> GetMean()-3*prof->GetRMS(), prof -> GetMean()+3*prof->GetRMS());
						prof -> Fit(fitFunc2_posCorr[index2],"QRS+");
						fitFunc2_posCorr[index2] -> SetLineColor(kRed);
						fitFunc2_posCorr[index2] -> SetLineWidth(2);
						fitFunc2_posCorr[index2] -> Draw("same");
				
						c -> Print(Form("%s/positionCorr/c_deltaT_totRatioCorr_vs_posX__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
						c -> Print(Form("%s/positionCorr/c_deltaT_totRatioCorr_vs_posX_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
						delete c;
						delete latex;
	
					}
				}
			}
	    }
    }// end draw 5th plots
  
  
  
  
    //------------------------
    //--- 6th loop over events
    //  if (useTrackInfo){
	std::map<double,TH1F*> h1_deltaT_energyRatioCorr_ave_bars;
    for(auto mapIt : trees) {
        ModuleEventClass* anEvent = new ModuleEventClass();
        mapIt.second -> SetBranchAddress("event",&anEvent);
      
        int nEntries = mapIt.second->GetEntries();
        for(int entry = 0; entry < nEntries; ++entry) {
	  
	        if( entry%100000 == 0 ) std::cout << ">>> 6th loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
	        mapIt.second -> GetEntry(entry);
	  
	        bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;                                                                             
	        if (!barFound) continue;
	  
	        int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
	        if( !accept[index1][entry] ) continue;
	  
	       int energyBinAverage=0;
			if(anEvent->nClusters==1) {energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;	}
			else if(anEvent->nClusters==2 && anEvent->energyL_post>0 && anEvent->energyR_post>0){energyBinAverage = FindBin(anEvent->energySum,ranges_doubleHits[index1])+1;}
			else if(anEvent->nClusters==2 && anEvent->energyL_post<0 && anEvent->energyR_post<0){
				energyBinAverage = -1;
				continue; // escludo eventi "PRE"
			}
			else if(anEvent->nClusters==3){	energyBinAverage=1;}
			

	        double  index2((100000000*anEvent->nClusters)+10000000*energyBinAverage+index1 );     
	  
	        long long deltaT = anEvent->timeR - anEvent->timeL;
	  
	  
	        // --- energy ratio + tot + phase corr
	        if( !fitFunc_energyRatioCorr[index2] )    continue;
	        if( !fitFunc_energyRatioCorr_totRatioCorr[index2] )    continue;
	        if( !p1_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean[index2] )      continue;     
	  
	        float energyRatioCorr = fitFunc_energyRatioCorr[index2]->Eval(anEvent->energyR/anEvent->energyL) - fitFunc_energyRatioCorr[index2]->Eval(fitFunc_energyRatio[index2]->GetParameter(1));
	  
	        float energyRatioCorr_totRatioCorr = fitFunc_energyRatioCorr_totRatioCorr[index2]->Eval(anEvent->totR/anEvent->totL) - 
	                                       fitFunc_energyRatioCorr_totRatioCorr[index2]->Eval(fitFunc_totRatio[index2]->GetParameter(1));
	  
	        float t1fineMean = 0.5* ( anEvent->t1fineR + anEvent->t1fineL );
	        int t1fineBin = p1_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean[index2]->FindBin(t1fineMean);
	        int t1fineBin2 = p1_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean[index2]->FindBin(h1_t1fineMean[index2]->GetMean());
	        float t1fineCorr = p1_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean[index2]->GetBinContent(t1fineBin) - p1_deltaT_energyRatioCorr_totRatioCorr_vs_t1fineMean[index2]->GetBinContent( t1fineBin2 );
	  
	        if( h1_deltaT_energyRatioCorr_totRatioCorr_phaseCorr[index2] == NULL ){
	            std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
	            h1_deltaT_energyRatioCorr_totRatioCorr_phaseCorr[index2] = new TH1F(Form("h1_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_%s",labelLR_energyBin.c_str()),"",2000,-12000.,12000.);         
	        } 

			if(h2_timeDiff_correlation[index2]==NULL){// deltaT correlation
				std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
				
                h2_timeDiff_correlation[index2] = new TH2F(Form("h2_timeDiff_correlation_%s",labelLR_energyBin.c_str()),"",500,-2000,2000,500,-2000,2000);
				p1_timeDiff_correlation[index2] = new TProfile(Form("p1_timeDiff_correlation_%s",labelLR_energyBin.c_str()),"",80,-2000,2000);
			}
	  
	        if (fabs(deltaT - energyRatioCorr)<10000){
	            h1_deltaT_energyRatioCorr_totRatioCorr_phaseCorr[index2] -> Fill( deltaT  - energyRatioCorr - energyRatioCorr_totRatioCorr - t1fineCorr);
	            //h1_deltaT_energyRatioCorr_totRatioCorr_phaseCorr[index2] -> Fill( deltaT  - energyRatioCorr - energyRatioCorr_totRatioCorr);
	        }    
	  
	        if(h1_deltaT_energyRatioCorr_ave_bars[index2]==NULL){
				std::string labelLR_energyBin_bars(Form("bars%02d-%02d_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->barID+1,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
                h1_deltaT_energyRatioCorr_ave_bars[index2] = new TH1F(Form("h1_deltaT_energyRatioCorr_ave_%s",labelLR_energyBin_bars.c_str()),"",2000,-12000.,12000.);
			}
			float deltaT_corr = deltaT  - energyRatioCorr - t1fineCorr;
             
			long long deltaT_post=0;
			float energyRatioCorr_post=0;
			float t1fineMean_post=0;
			int t1fineBin_post=0;
			int t1fineBin2_post=0;
			float t1fineCorr_post=0;
			float deltaT_post_corr=0;
			if(anEvent->nClusters==2 && anEvent->energyL_post>0 && anEvent->energyR_post>0){
				deltaT_post = anEvent->timeR_post - anEvent->timeL_post;
				if( !fitFunc_energyRatioCorr[index2+1] )	continue;
	            if( !p1_deltaT_energyRatioCorr_vs_t1fineMean[index2+1] )continue;
				energyRatioCorr_post = fitFunc_energyRatioCorr[index2+1]->Eval(anEvent->energyR_post/anEvent->energyL_post) - fitFunc_energyRatioCorr[index2+1]->Eval(fitFunc_energyRatio[index2+1]->GetParameter(1));
				t1fineMean_post = 0.5* ( anEvent->t1fineR_post + anEvent->t1fineL_post );
	            t1fineBin_post = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2+1]->FindBin(t1fineMean_post);
	            t1fineBin2_post = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2+1]->FindBin(h1_t1fineMean[index2+1]->GetMean());
	            t1fineCorr_post = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2+1]->GetBinContent(t1fineBin_post) - p1_deltaT_energyRatioCorr_vs_t1fineMean[index2+1]->GetBinContent( t1fineBin2_post );
				deltaT_post_corr = deltaT_post - energyRatioCorr_post - t1fineCorr_post;
			}


			if(anEvent->nClusters==2 && anEvent->energyL_post>0 && anEvent->energyR_post>0){
				    float E_mean = 0.5*(anEvent->energyL+anEvent->energyR);
					float E_mean_post = 0.5 * (anEvent->energyL_post + anEvent->energyR_post);

					h2_timeDiff_correlation[index2]->Fill(deltaT_post_corr,deltaT_corr);
					p1_timeDiff_correlation[index2]->Fill(deltaT_post_corr,deltaT_corr);

					if(!fitFunc_energyRatioCorr_bars_ave[index2]) continue;
					float energyRatioCorr_bars_ave = fitFunc_energyRatioCorr_bars_ave[index2]->Eval(E_mean/E_mean_post) - fitFunc_energyRatioCorr_bars_ave[index2]->Eval(1.);
					
					h1_deltaT_energyRatioCorr_ave_bars[index2] -> Fill(0.5*(deltaT_corr + deltaT_post_corr) - energyRatioCorr_bars_ave);   
			}


	        // --- positon corrections
	        if (!useTrackInfo) continue;
	  
	  
	        // --- energy ratio + phase + positon corrections
	        if( !fitFunc_energyRatioCorr[index2] )	continue;
	        if( !p1_deltaT_energyRatioCorr_vs_t1fineMean[index2] )	continue;
	        if( !fitFunc1_posCorr[index2] )	continue;
	  
	  
	        t1fineMean = 0.5* ( anEvent->t1fineR + anEvent->t1fineL );
	        t1fineBin = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]->FindBin(t1fineMean);
	        t1fineBin2 = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]->FindBin(h1_t1fineMean[index2]->GetMean());
	        t1fineCorr = p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]->GetBinContent(t1fineBin) - p1_deltaT_energyRatioCorr_vs_t1fineMean[index2]->GetBinContent( t1fineBin2 );
	  
	        float posCorr = fitFunc1_posCorr[index2]->Eval(anEvent->x) - fitFunc1_posCorr[index2]->Eval( p1_deltaT_energyRatioCorr_vs_posX[index2]->GetMean());
	  
	        if( h1_deltaT_energyRatioPhasePosCorr[index2] == NULL ) {
	            std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
	            h1_deltaT_energyRatioPhasePosCorr[index2] = new TH1F(Form("h1_deltaT_energyRatioPhasePosCorr_%s",labelLR_energyBin.c_str()),"",2000,-12000.,12000.);
	        }
	  
	        if (fabs(deltaT - energyRatioCorr)<10000 ){
	            h1_deltaT_energyRatioPhasePosCorr[index2] -> Fill( deltaT  - energyRatioCorr - t1fineCorr - posCorr);
	        }
	  
	        // --- tot ratio + phase + positon corrections
	        if( !fitFunc_totRatioCorr[index2] )	continue;
	        if( !p1_deltaT_totRatioCorr_vs_t1fineMean[index2] )	continue;
	        if( !fitFunc2_posCorr[index2] )     continue;
	  
	        float totRatioCorr = fitFunc_totRatioCorr[index2]->Eval(anEvent->totR/anEvent->totL) -
	        fitFunc_totRatioCorr[index2]->Eval(fitFunc_totRatio[index2]->GetParameter(1));
	  
	  
	        t1fineBin = p1_deltaT_totRatioCorr_vs_t1fineMean[index2]->FindBin(t1fineMean);
	        t1fineBin2 = p1_deltaT_totRatioCorr_vs_t1fineMean[index2]->FindBin(h1_t1fineMean[index2]->GetMean());
	        t1fineCorr = p1_deltaT_totRatioCorr_vs_t1fineMean[index2]->GetBinContent(t1fineBin) - p1_deltaT_totRatioCorr_vs_t1fineMean[index2]->GetBinContent( t1fineBin2 );
	  
	        posCorr = fitFunc2_posCorr[index2]->Eval(anEvent->x) - fitFunc2_posCorr[index2]->Eval( p1_deltaT_totRatioCorr_vs_posX[index2]->GetMean());
	  
	        if( h1_deltaT_totRatioPhasePosCorr[index2] == NULL ) {
	            std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d_%02dHits",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage,anEvent->nClusters));
	            h1_deltaT_totRatioPhasePosCorr[index2] = new TH1F(Form("h1_deltaT_totRatioPhasePosCorr_%s",labelLR_energyBin.c_str()),"",2000,-12000.,12000.);
	        }
	  
	        if (fabs(deltaT - totRatioCorr)>10000 ) continue;
	        h1_deltaT_totRatioPhasePosCorr[index2] -> Fill( deltaT  - totRatioCorr - t1fineCorr - posCorr);
	    }
        std::cout << std::endl;
    }
  
  
    //------------------
    // draw 6th plots
  
    for(auto stepLabel : stepLabels) {
        float Vov = map_Vovs[stepLabel];
        float vth1 = map_ths[stepLabel];
      
        for(int iBar = 0; iBar < 16; ++iBar) {
	        bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
	        if (!barFound) continue;
	  
	        std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
	  
	        int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );
	        if( !ranges["L-R"][index1] ) continue;
	  
	        int nEnergyBins = ranges["L-R"][index1]->size()-1;
	  
	        for(int i=1;i<=3;i++){
				for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin) {
					double  index2((100000000*i)+10000000*iEnergyBin + index1 );
					std::string labelLR_energyBin;
					if(i==1){labelLR_energyBin = Form("%s_energyBin%02d_S",labelLR.c_str(),iEnergyBin);}
					else if(i==2){labelLR_energyBin = Form("%s_energyBin%02d_D",labelLR.c_str(),iEnergyBin);}
					else if(i==3){labelLR_energyBin = Form("%s_energyBin%02d_T",labelLR.c_str(),iEnergyBin);}
			  

					std::cout << labelLR_energyBin.c_str()<<std::endl;
			  
					// -- energy+tot+phase corr
					if ( !h1_deltaT_energyRatioCorr_totRatioCorr_phaseCorr[index2] ) continue;
					c = new TCanvas(Form("c_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_%s",labelLR_energyBin.c_str()));		
					histo = h1_deltaT_energyRatioCorr_totRatioCorr_phaseCorr[index2];
					histo -> SetLineWidth(2);
					histo -> SetLineColor(kGreen+1);
					histo -> SetMarkerColor(kGreen+1);
			  
					TF1* fitFunc = new TF1(Form("fitFunc_energyTotPhaseCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					drawDeltaT(c, histo, fitFunc, "corrected", "allCorr","");
			  
					outFile -> cd();
					histo -> Write();
			        if(!h1_deltaT_energyRatioCorr_totRatioCorr[index2]) continue;
					histo = h1_deltaT_energyRatioCorr_totRatioCorr[index2];
					histo -> SetLineWidth(2);
					histo -> SetLineColor(kAzure);
					histo -> SetMarkerColor(kAzure);
			  
					fitFunc = new TF1(Form("fitFunc_energyRatioCorr_totRatioCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);                                            
					drawDeltaT(c, histo, fitFunc, "corrected", "en+totCorr","same");
			  
					c -> Print(Form("%s/CTR_energyRatioCorr_totRatioCorr_phaseCorr/c_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/CTR_energyRatioCorr_totRatioCorr_phaseCorr/c_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete c;

					if(i==2){
                        // --- draw deltaT correlation plots
						c = new TCanvas(Form("c_timeDiffCorrelation_%s",labelLR_energyBin.c_str()),Form("c_timeDiffCorrelation_%s",labelLR_energyBin.c_str()));
						if(!p1_timeDiff_correlation[index2]) continue;
					    c -> SetGridy();

						h2 = h2_timeDiff_correlation[index2]; 
					    h2 -> SetTitle(Form(";bar %02d #DeltaT [ps];bar %02d #DeltaT [ps]",iBar+1,iBar));
					    h2 -> Draw("colz");

					    prof = p1_timeDiff_correlation[index2];
						prof -> SetMarkerSize(0.5);
					    prof -> Draw("psame");

						latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
					    latex -> SetNDC();
					    latex -> SetTextFont(42);
					    latex -> SetTextSize(0.04);
					    latex -> SetTextColor(kRed);
					    latex -> Draw("same");

						c -> Print(Form("%s/timeCorrelation/c_timeDiffCorrelation_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
						c -> Print(Form("%s/timeCorrelation/c_timeDiffCorrelation_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					    delete c;
					    delete latex;

                        // --- draw deltaT_ave plots
						if(!h1_deltaT_energyRatioCorr_ave_bars[index2]) continue;
						c = new TCanvas(Form("c_deltaT_energyRatioCorr_ave_bars_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioCorr_ave_bars_%s",labelLR_energyBin.c_str()));
					    histo = h1_deltaT_energyRatioCorr_ave_bars[index2];
					    histo -> SetLineWidth(2);
					    histo -> SetLineColor(kGreen+1);
					    histo -> SetMarkerColor(kGreen+1);
			  
					    fitFunc = new TF1(Form("fitFunc_energyCorr_ave_bars_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					    drawDeltaT(c, histo, fitFunc, Form("%02d - %02d average ",iBar,iBar+1), "allCorr","");
			  
					    outFile -> cd();
					    histo -> Write();
			        			  
					    c -> Print(Form("%s/CTR_energyRatioCorr_phaseCorr/c_deltaT_bars_energyRatioCorr_ave__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					    c -> Print(Form("%s/CTR_energyRatioCorr_phaseCorr/c_deltaT_bars_energyRatioCorr_ave__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					    delete c;
					}
			  
			  
					// position corrections
					if (!useTrackInfo) continue;
					if(!h1_deltaT_energyRatioPhasePosCorr[index2]) continue;
			  
					// -- energy, phase, pos corr deltaT
					c = new TCanvas(Form("c_deltaT_energyRatioPhasePosCorr_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioPhasePosCorr_%s",labelLR_energyBin.c_str()));		
					histo = h1_deltaT_energyRatioPhasePosCorr[index2];
					histo -> SetLineWidth(2);
					histo -> SetLineColor(kMagenta);
					histo -> SetMarkerColor(kMagenta);
			  
					fitFunc = new TF1(Form("fitFunc_energyPhasePosCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					drawDeltaT(c, histo, fitFunc, "corrected", "pos.Corr","");
			  
					outFile -> cd();
					histo -> Write();
			  
			  
					// -- energy, phase corr deltaT    
					histo = h1_deltaT_energyRatioPhaseCorr[index2];
					histo -> SetLineWidth(2);
					histo -> SetLineColor(kGreen+1);
					histo -> SetMarkerColor(kGreen+1);
			  
					fitFunc = new TF1(Form("fitFunc_energyPhaseCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					drawDeltaT(c, histo, fitFunc, "corrected", "ph.Corr","same");
			  
					outFile -> cd();
					histo -> Write();
			  
					c -> Print(Form("%s/CTR_energyRatioCorr_phaseCorr_posCorr/c_deltaT_energyRatioPhasePosCorr_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/CTR_energyRatioCorr_phaseCorr_posCorr/c_deltaT_energyRatioPhasePosCorr_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete c;
			  
			  
					// -- tot, phase, pos corr deltaT
					c = new TCanvas(Form("c_deltaT_totRatioPhasePosCorr_%s",labelLR_energyBin.c_str()),Form("c_deltaT_totRatioPhasePosCorr_%s",labelLR_energyBin.c_str()));
					histo = h1_deltaT_totRatioPhasePosCorr[index2];
					histo -> SetLineWidth(2);
					histo -> SetLineColor(kMagenta);
					histo -> SetMarkerColor(kMagenta);
			  
					fitFunc = new TF1(Form("fitFunc_totRatioPhasePosCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					drawDeltaT(c, histo, fitFunc, "corrected", "pos.Corr",""); 
			  
					outFile -> cd();
					histo -> Write();
			  
			  
					// -- tot, phase corr deltaT    
					histo = h1_deltaT_totRatioPhaseCorr[index2];
					histo -> SetLineWidth(2);
					histo -> SetLineColor(kGreen+1);
					histo -> SetMarkerColor(kGreen+1);
			  
					fitFunc = new TF1(Form("fitFunc_totPhaseCorr_%s",labelLR_energyBin.c_str()),"gaus",-10000, 10000);
					drawDeltaT(c, histo, fitFunc, "corrected", "totCorr","same");
			  
					outFile -> cd();
					histo -> Write();
			  
					c -> Print(Form("%s/CTR_totRatioCorr_phaseCorr_posCorr/c_deltaT_totRatioPhasePosCorr_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
					c -> Print(Form("%s/CTR_totRatioCorr_phaseCorr_posCorr/c_deltaT_totRatioPhasePosCorr_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
					delete c;
				}
			}
	    }
    }
  std::cout<<"\ntotal events at start for new selection: "<<totEntries<<std::endl;
  std::cout<<"2nd loop -> accepted "<<accepted2<<" events"<<std::endl;
  std::cout<<"3rd loop -> accepted "<<accepted3<<" events"<<std::endl;
  std::cout<<"4th loop -> accepted "<<accepted4<<" events"<<std::endl;
  std::cout<<"5th loop -> accepted "<<accepted5<<" events"<<std::endl;
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
