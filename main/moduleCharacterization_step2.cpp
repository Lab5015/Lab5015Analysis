#include "interface/AnalysisUtils.h"
#include "interface/Na22SpectrumAnalyzer.h"
#include "interface/Na22SpectrumAnalyzerSingleBar.h"
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


// find energy bins
void GetEnergyBins(TH1F *h, std::vector<float> *r, std::map<int, float> b){

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



// ============  ********* MAIN  ************ =============== //
int main(int argc, char** argv)
{
  setTDRStyle();
  
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
  //system(Form("rm -r %s", plotDir.c_str())); // questo non va bene se stiamo lavorando in parallelo
  system(Form("mkdir -p %s",plotDir.c_str()));
  system(Form("mkdir -p %s/qfine/",plotDir.c_str()));
  system(Form("mkdir -p %s/tot/",plotDir.c_str()));
  system(Form("mkdir -p %s/energy/",plotDir.c_str()));
  system(Form("mkdir -p %s/energyRatio/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR_energyRatioCorr/",plotDir.c_str()));
  
  
  //copia il file .php che si trova in una cartella fuori plotDir in plotDir definita nel config file
  system(Form("cp %s/../index.php %s",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/qfine/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/tot/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/energy/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/energyRatio/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/CTR/",plotDir.c_str(),plotDir.c_str()));
  system(Form("cp %s/../index.php %s/CTR_energyRatioCorr/",plotDir.c_str(),plotDir.c_str()));
  
  
  std::vector<std::string> LRLabels;
  LRLabels.push_back("L");
  LRLabels.push_back("R");
  LRLabels.push_back("L-R");
  
  
  // //--- get cuts per bar / Vov
  // std::map<unsigned int,std::map<float,float> > cut_qfineAcc;
  // std::map<unsigned int,std::map<float,float> > cut_totAcc;
  // std::map<unsigned int,std::map<float,float> > cut_energyAcc;
  // std::map<unsigned int,std::map<float,float> > cut_energyFitMin;
  // std::map<unsigned int,std::map<float,float> > cut_energyFitMax;
  // std::map<unsigned int, float > noise;
  // for(auto ch :  channels)
  // {
  //   int chID = opts.GetOpt<int>(Form("%s.chID",ch.c_str()));
  //   std::vector<float> Vovs          = opts.GetOpt<std::vector<float> >(Form("%s.Vovs",ch.c_str()));
  //   std::vector<float> qfineMins     = opts.GetOpt<std::vector<float> >(Form("%s.qfineMins",ch.c_str()));
  //   std::vector<float> totMins       = opts.GetOpt<std::vector<float> >(Form("%s.totMins",ch.c_str()));
  //   std::vector<float> energyMins    = opts.GetOpt<std::vector<float> >(Form("%s.energyMins",ch.c_str()));
  //   std::vector<float> energyFitMins = opts.GetOpt<std::vector<float> >(Form("%s.energyFitMins",ch.c_str()));
  //   std::vector<float> energyFitMaxs = opts.GetOpt<std::vector<float> >(Form("%s.energyFitMaxs",ch.c_str()));
  //   noise[chID] = opts.GetOpt<float>(Form("%s.noise",ch.c_str()));
  //   int iter = 0;
  //   for(auto Vov : Vovs)
  //   {
  //     cut_qfineAcc[chID][Vov]     = qfineMins.at(iter);
  //     cut_totAcc[chID][Vov]       = totMins.at(iter);
  //     cut_energyAcc[chID][Vov]    = energyMins.at(iter);
  //     cut_energyFitMin[chID][Vov] = energyFitMins.at(iter);
  //     cut_energyFitMax[chID][Vov] = energyFitMaxs.at(iter);
  //     ++iter;
  //   }
  // }
  // std::map<std::string,float> cut_energyMin;
  // std::map<std::string,float> cut_energyMax;
  

  float energyMax = opts.GetOpt<float>("Plots.energyMax");
  
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
  while( (object = next()) )
  {
    std::string name(object->GetName());
    std::vector<std::string> tokens = GetTokens(name,'_');
    std::size_t found;
     
    found = name.find("data_");
    //tree
    if( found!=std::string::npos )
    {
      std::string label(Form("%s_%s_%s",tokens[1].c_str(),tokens[2].c_str(),tokens[3].c_str()));
      trees[label] = (TTree*)( inFile->Get(name.c_str()) );
    }
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
    }
  }
  std::sort(stepLabels.begin(),stepLabels.end());
  stepLabels.erase(std::unique(stepLabels.begin(),stepLabels.end()),stepLabels.end());
  
  
  //--- define histograms
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileNameStep2");
  TFile* outFile = TFile::Open(outFileName.c_str(),"RECREATE");
  outFile->cd();

  std::map<int,TTree*> outTrees;
  std::map<double,TTree*> outTrees2;
  float energy511LR;
  float energy1275LR;
  float energy1786LR;
  float energy511R;
  float energy1275R;
  float energy1786R;
  float energy511L;
  float energy1275L;
  float energy1786L;
  float theIndex;
  float timeRes;


  std::map<double,TH1F*> h1_energyRatio;
  
  std::map<double,TH1F*> h1_deltaT_raw;
  std::map<double,TH1F*> h1_deltaT;
  std::map<double,TProfile*> p1_deltaT_vs_energyRatio;

  std::map<double,TH1F*> h1_deltaT_energyRatioCorr;

  std::map<std::string, std::map<int, std::vector<float>*> > ranges; //ranges[LRlabel][index]
  std::map<std::string, std::map<int, std::map<std::string,std::pair<float,float> > > > peaks;	//peaks[LRlabel][index][energyPeak]
  std::map<std::string, std::map<int, std::map<int,float> > > energyBin; // energyBin[LRlabel][index]
  
  
  //--- get plot settings
  TCanvas* c;
  float* vals = new float[6];
  TLatex* latex;
  TH1F* histo;
  TProfile* prof;

  
  //------------------
  //--- draw 1st plots
  std::string source = opts.GetOpt<std::string>("Input.sourceName");
  std::string Na22 = "Na22";
  std::string Na22SingleBar = "Na22SingleBar";
  std::string Co60 = "Co60";
  std::string Co60SumPeak = "Co60SumPeak";	
  std::string Laser = "Laser";	
  
  for(auto stepLabel : stepLabels)
  {
    float Vov = map_Vovs[stepLabel];
    float vth1 = map_ths[stepLabel];
    std::string VovLabel(Form("Vov%.1f",Vov));
    std::string thLabel(Form("th%02.0f",vth1));
    
    
    //--------------------------------------------------------
    
    // --- loop over bars
    for(int iBar = 0; iBar < 16; ++iBar)
    {
      int index( (10000*int(Vov*100.)) + (100*vth1) + iBar );
      
      outTrees[index] = new TTree(Form("data_bar%02dL-R_Vov%.1f_th%02.0f",iBar,Vov,vth1),Form("data_bar%02dL-R_Vov%.1f_th%02.0f",iBar,Vov,vth1));
      outTrees[index] -> Branch("energyPeak511LR",&energy511LR);
      outTrees[index] -> Branch("energyPeak1275LR",&energy1275LR);
      outTrees[index] -> Branch("energyPeak1786LR",&energy1786LR);
      outTrees[index] -> Branch("energyPeak511L",&energy511L);
      outTrees[index] -> Branch("energyPeak1275L",&energy1275L);
      outTrees[index] -> Branch("energyPeak1786L",&energy1786L);
      outTrees[index] -> Branch("energyPeak511R",&energy511R);
      outTrees[index] -> Branch("energyPeak1275R",&energy1275R);
      outTrees[index] -> Branch("energyPeak1786R",&energy1786R);
      outTrees[index] -> Branch("indexID",&theIndex);
      
      
      // -- loop over L, R, LR
      for(auto LRLabel : LRLabels )
      {
	//label histo
        std::string label(Form("bar%02d%s_%s",iBar,LRLabel.c_str(),stepLabel.c_str()));

        latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d%s}{V_{OV} = %.1f V, th. = %d DAC}",iBar,LRLabel.c_str(),Vov,int(vth1)));
        if (LRLabel == "L-R") 
          latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.1f V, th. = %d DAC}",iBar,Vov,int(vth1)));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kRed);
          

        // -- qfine and tot only per L, R
        if (LRLabel == "R" || LRLabel == "L") {

          histo = (TH1F*)( inFile->Get(Form("h1_qfine_%s",label.c_str())) );
          if( !histo ) continue;
          if( histo->GetEntries() < 100 ) continue;
          
          c = new TCanvas(Form("c_qfine_%s",label.c_str()),Form("c_qfine_%s",label.c_str()));
          gPad -> SetLogy();
                  
          histo -> SetTitle(";Q_{fine} [ADC];entries");
          histo -> SetLineColor(kRed);
          histo -> Draw();
          histo -> GetXaxis() -> SetRangeUser(0.,600.);
          // TLine* line_qfineAcc1 = new TLine(cut_qfineAcc[chID][Vov],histo->GetMinimum(),cut_qfineAcc[chID][Vov],histo->GetMaximum());
          // line_qfineAcc1 -> SetLineColor(kBlack);
          // line_qfineAcc1 -> Draw("same");
          latex -> Draw("same");      
          histo -> Write();
          c -> Print(Form("%s/qfine/c_qfine__%s.png",plotDir.c_str(),label.c_str()));
          c -> Print(Form("%s/qfine/c_qfine__%s.pdf",plotDir.c_str(),label.c_str()));
          delete c;
        
        
          c = new TCanvas(Form("c_tot_%s",label.c_str()),Form("c_tot_%s",label.c_str()));
          // gPad -> SetLogy();
        
          histo = (TH1F*)( inFile->Get(Form("h1_tot_%s",label.c_str())) );
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
        
        // -- draw energy
        c = new TCanvas(Form("c_energy_%s",label.c_str()),Form("c_energy_%s",label.c_str()));
        gPad -> SetLogy();
        
        histo = (TH1F*)( inFile->Get(Form("h1_energy_%s",label.c_str())) );      
        histo -> SetTitle(";energy [a.u.];entries");
        histo -> SetLineColor(kRed);
        histo -> SetLineWidth(2);
        histo -> Draw();
	


        // --- look for peaks and define energy ranges
        if( source.compare(Na22) && source.compare(Na22SingleBar) && source.compare(Co60) && source.compare(Co60SumPeak) && source.compare(Laser)){
          std::cout << " Missspelled radioactive source " << std::endl;
          return(0);
        }
	
        ranges[LRLabel][index] = new std::vector<float>;

        // Na22 or Co60 spectrum
        if(!source.compare(Na22)  || !source.compare(Na22SingleBar) ||  !source.compare(Co60) ){        
          std::string firstPeak = "";
          if (!source.compare(Na22))          {
            peaks[LRLabel][index] = Na22SpectrumAnalyzer(histo,ranges[LRLabel][index]);
            firstPeak = "0.511 MeV";
          }
          if (!source.compare(Na22SingleBar)) {
            peaks[LRLabel][index] = Na22SpectrumAnalyzerSingleBar(histo,ranges[LRLabel][index]);
            firstPeak = "0.511 MeV";
          }
          if (!source.compare(Co60))          {
            peaks[LRLabel][index] = Co60SpectrumAnalyzer_2Peaks(histo,ranges[LRLabel][index]);
            firstPeak = "1.173 MeV";
          }
          
          
          if (peaks[LRLabel][index][firstPeak].first > 0){
            histo -> GetXaxis() -> SetRangeUser(0.,5.*peaks[LRLabel][index][firstPeak].first);
          }
            
          if (peaks[LRLabel][index][firstPeak].first== -9999){
            histo -> GetXaxis() -> SetRangeUser(0., energyMax); 
            peaks[LRLabel].erase(index);
            ranges[LRLabel].erase(index);
            if (!source.compare(Na22) || !source.compare(Na22SingleBar) ) peaks[LRLabel][index]["0.511 MeV"].first = -10;
            if (!source.compare(Na22) || !source.compare(Na22SingleBar) ) peaks[LRLabel][index]["1.275 MeV"].first = -10;
            if (!source.compare(Na22SingleBar))                           peaks[LRLabel][index]["1.786 MeV"].first = -10;
            if (!source.compare(Co60))                                    peaks[LRLabel][index]["1.173 MeV"].first = -10;
            if (!source.compare(Co60))                                    peaks[LRLabel][index]["1.332 MeV"].first = -10;
            if (!source.compare(Co60))                                    peaks[LRLabel][index]["2.505 MeV"].first = -10;
          }
          
          histo -> GetYaxis() -> SetRangeUser(3.,5.*histo->GetMaximum());
          
          if(peaks[LRLabel][index][firstPeak].first != -10){
            GetEnergyBins(histo, ranges[LRLabel][index], energyBin[LRLabel][index]);
          }          
        } // end if Na22/Co60
          
        
        // -- if Co60 sum peak or laser, we don't use the spectrum analyzers - just gaussian fit to the energy peak.
        if(!source.compare(Co60SumPeak) ||  !source.compare(Laser)){ 
            
          TF1* fitFunc = new TF1 ( Form("funcL_%d",index), "gaus(0)", histo -> GetBinCenter(histo->GetMaximumBin()) - 2*histo->GetRMS(), histo -> GetBinCenter(histo->GetMaximumBin()) + 2*histo->GetRMS());
        
          histo -> Fit( fitFunc, "NQR");
          fitFunc -> SetParameter(0, fitFunc->GetParameter(0));
          fitFunc -> SetParameter(1, fitFunc->GetParameter(1));
          fitFunc -> SetParameter(2, fitFunc->GetParameter(2));
          fitFunc -> SetLineColor(kBlack);
          fitFunc -> SetLineWidth(2);
          
          histo -> Fit ( fitFunc, "QRS+", "",  fitFunc->GetParameter(1) - fitFunc->GetParameter(2) *5, fitFunc->GetParameter(1) + fitFunc->GetParameter(2) *5);
          
          // -- questa condizione er solo per il laser... puo' andare bene tenerla anche per 60Co?
          if ((fitFunc->GetParameter(1) - fitFunc->GetParameter(2) *5)>0) 
            ranges[LRLabel][index] -> push_back( fitFunc->GetParameter(1) - fitFunc->GetParameter(2) *5);
          else 
            ranges[LRLabel][index] -> push_back(0);
          ranges[LRLabel][index] -> push_back( fitFunc->GetParameter(1) + fitFunc->GetParameter(2) *5);

          // use energy511XX branch to save the peak value 
          if (LRLabel == "L") energy511L = fitFunc->GetParameter(1);					
          if (LRLabel == "R") energy511R = fitFunc->GetParameter(1);					
          if (LRLabel == "L-R") energy511LR = fitFunc->GetParameter(1);					

          for(auto range: (*ranges[LRLabel][index])){
            float yval = std::max(10., histo->GetBinContent(histo->FindBin(range)));
            TLine* line = new TLine(range,3.,range, yval);
            line -> SetLineWidth(1);
            line -> SetLineStyle(7);
            line -> Draw("same");
           }

           GetEnergyBins(histo, ranges[LRLabel][index], energyBin[LRLabel][index]);
         }// end Co60SumPeak or laser	



         // -- draw and print energy plots 
         latex -> Draw("same"); 
         outFile -> cd();
         histo->Write();     
         c -> Print(Form("%s/energy/c_energy__%s.png",plotDir.c_str(),label.c_str()));
         c -> Print(Form("%s/energy/c_energy__%s.pdf",plotDir.c_str(),label.c_str()));
         delete c;

       }// end loop over L, R, L-R labels



       // --- now fill the output tree.  
       // *********** MM: Ma serve ? non basterebbe leggere i valori di posizione dei picchi e numero di eventi dagli istogrammi che vengono comunque salvati nel file di output? 

       theIndex = index;

       if(!source.compare(Na22)  || !source.compare(Na22SingleBar) ){
         energy511R  = peaks["R"][index]["0.511 MeV"].first;		
         energy511L  = peaks["L"][index]["0.511 MeV"].first;		
         energy511LR = peaks["L-R"][index]["0.511 MeV"].first;

         energy1275R = peaks["R"][index]["1.275 MeV"].first;
         energy1275L = peaks["L"][index]["1.275 MeV"].first;
         energy1275LR = peaks["L-R"][index]["1.275 MeV"].first;

         if( !source.compare(Na22SingleBar) ) {
           energy1786R = peaks["R"][index]["1.786 MeV"].first;
           energy1786L = peaks["L"][index]["1.786 MeV"].first;
           energy1786LR = peaks["L-R"][index]["1.786 MeV"].first;
         }       
       }

       if(!source.compare(Co60)){
         energy511R  = peaks["R"][index]["1.173 MeV"].first;
         energy511L  = peaks["L"][index]["1.173 MeV"].first;
         energy511LR = peaks["L-R"][index]["1.173 MeV"].first;

         energy1275R = peaks["R"][index]["1.332 MeV"].first;
         energy1275L = peaks["L"][index]["1.332 MeV"].first;
         energy1275LR = peaks["L-R"][index]["1.332 MeV"].first;

         energy1786R = peaks["R"][index]["2.505 MeV"].first;
         energy1786L = peaks["L"][index]["2.505 MeV"].first;
         energy1786LR = peaks["L-R"][index]["2.505 MeV"].first;
       }


       outTrees[index] -> Fill();


     }// -- end loop over bars

   } // -- end loop over stepLabels


   // ---  end 1st plots




   //------------------------
   //--- 2nd loop over events
   std::map<int,std::map<int,bool> > accept;

   for(auto mapIt : trees)
   {
     ModuleEventClass* anEvent = new ModuleEventClass();
     mapIt.second -> SetBranchAddress("event",&anEvent);

     int nEntries = mapIt.second->GetEntries();
     for(int entry = 0; entry < nEntries; ++entry)
     {
       if( entry%100000 == 0 ) std::cout << ">>> 2nd loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
       mapIt.second -> GetEntry(entry);

       int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

       accept[index1][entry] = false;

       if(!ranges["L-R"][index1] ) continue;

       int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

       if( energyBinAverage < 1 ) continue;

       double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

       accept[index1][entry] = true;

       if( h1_energyRatio[index2] == NULL )
       {
         std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.1f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));

         h1_energyRatio[index2] = new TH1F(Form("h1_energyRatio_%s",labelLR_energyBin.c_str()),"",1000,0.,5.);

         h1_deltaT_raw[index2] = new TH1F(Form("h1_deltaT_raw_%s",labelLR_energyBin.c_str()),"",1250,-5000.,5000.);
       }

       if ((anEvent->energyR / anEvent->energyL >0) && (anEvent->energyR / anEvent->energyL <5)){
         h1_energyRatio[index2] -> Fill( anEvent->energyR / anEvent->energyL );						     
         h1_deltaT_raw[index2] -> Fill( anEvent->timeR-anEvent->timeL );
       }
     }
     std::cout << std::endl;
   }


   //------------------
   //--- draw 2nd plots
   std::map<double,float> CTRMeans;
   std::map<double,float> CTRSigmas;

   std::map<double,TF1*> fitFunc_energyRatio;

   for(auto mapIt : h1_deltaT_raw)
   {
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

   for(auto stepLabel : stepLabels)
   {
     float Vov = map_Vovs[stepLabel];
     float vth1 = map_ths[stepLabel];

     for(int iBar = 0; iBar < 16; ++iBar)
     {
       int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );
       std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));

       if( !ranges["L-R"][index1] ) continue;
       int nEnergyBins = ranges["L-R"][index1]->size()-1;

       for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
       {
         //if (ranges["L-R"][index1]->at(iEnergyBin)<0) continue;
         double index2( 10000000*iEnergyBin+index1 );

         if (!h1_energyRatio[index2]) continue;

         std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));


         c = new TCanvas(Form("c_energyRatio_%s",labelLR_energyBin.c_str()),Form("c_energyRatio_%s",labelLR_energyBin.c_str()));
         histo = h1_energyRatio[index2];
         histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
         histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
         histo -> SetTitle(Form(";energy_{right} / energy_{left};entries"));
         histo -> SetLineColor(kRed);
         histo -> SetLineWidth(2);
         histo -> Draw();
         histo -> Write();

         fitFunc_energyRatio[index2] = new TF1(Form("fitFunc_energyRatio_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
         histo -> Fit(fitFunc_energyRatio[index2],"QNRS");
         histo -> Fit(fitFunc_energyRatio[index2],"QSR+","",fitFunc_energyRatio[index2]->GetParameter(1)-2.*fitFunc_energyRatio[index2]->GetParameter(2),fitFunc_energyRatio[index2]->GetParameter(1)+2.*fitFunc_energyRatio[index2]->GetParameter(2));
         histo -> Fit(fitFunc_energyRatio[index2],"QSR+","",fitFunc_energyRatio[index2]->GetParameter(1)-2.*fitFunc_energyRatio[index2]->GetParameter(2),fitFunc_energyRatio[index2]->GetParameter(1)+2.*fitFunc_energyRatio[index2]->GetParameter(2));

         fitFunc_energyRatio[index2] -> SetLineColor(kBlack);
         fitFunc_energyRatio[index2] -> SetLineWidth(2);
         fitFunc_energyRatio[index2] -> Draw("same");

         latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.1f V, th. = %d DAC}",iBar,Vov,int(vth1)));
         latex -> SetNDC();
         latex -> SetTextFont(42);
         latex -> SetTextSize(0.04);
         latex -> SetTextColor(kRed);
         latex -> Draw("same");

         c -> Print(Form("%s/energyRatio/c_energyRatio__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
         c -> Print(Form("%s/energyRatio/c_energyRatio__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
         delete c;


       } // --- end loop over energy bins
     }// --- end loop ober bars
   }// --- end loop over stepLabels



   //------------------------
   //--- 3rd loop over events
   for(auto mapIt : trees)
   {
     ModuleEventClass* anEvent = new ModuleEventClass();
     mapIt.second -> SetBranchAddress("event",&anEvent);

     int nEntries = mapIt.second->GetEntries();
     for(int entry = 0; entry < nEntries; ++entry)
     {
       if( entry%100000 == 0 ) std::cout << ">>> 3rd loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
       mapIt.second -> GetEntry(entry);

       int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

       if( !accept[index1][entry] ) continue;

       int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

       double  index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

       float energyRatioMean = fitFunc_energyRatio[index2]->GetParameter(1);
       float energyRatioSigma = fitFunc_energyRatio[index2]->GetParameter(2);


       if( fabs(anEvent->energyR/anEvent->energyL-energyRatioMean) > 2.*energyRatioSigma )
       {
         accept[index1][entry] = false;
         continue;
       }


       if( (anEvent->energyR/anEvent->energyL)>5 || (anEvent->energyR/anEvent->energyL)<0)
       {
         accept[index1][entry] = false;
         continue;
       }



       long long deltaT = anEvent->timeR - anEvent->timeL;

       if( h1_deltaT[index2] == NULL )
       {
         std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.1f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));

         h1_deltaT[index2] = new TH1F(Form("h1_deltaT_%s",labelLR_energyBin.c_str()),"",1000,-5000,5000.);

         p1_deltaT_vs_energyRatio[index2] = new TProfile(Form("p1_deltaT_vs_energyRatio_%s",labelLR_energyBin.c_str()),"",50,energyRatioMean-3.*energyRatioSigma,energyRatioMean+3.*energyRatioSigma);
       }

       h1_deltaT[index2] -> Fill( deltaT );    

       float timeLow = CTRMeans[index2] - 1.* CTRSigmas[index2];
       float timeHig = CTRMeans[index2] + 1.* CTRSigmas[index2];

       if( ( deltaT > timeLow ) && ( deltaT < timeHig ) )
         p1_deltaT_vs_energyRatio[index2] -> Fill( anEvent->energyR/anEvent->energyL,deltaT );
     }
     std::cout << std::endl;
   }


   //------------------
   //--- draw 3rd plots
   std::map<double,TF1*> fitFunc_energyRatioCorr;

   for(auto stepLabel : stepLabels)
   {
     float Vov = map_Vovs[stepLabel];
     float vth1 = map_ths[stepLabel];

     for(int iBar = 0; iBar < 16; ++iBar)
     {
       int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );

       std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));

       if( !ranges["L-R"][index1] ) continue;

       int nEnergyBins = ranges["L-R"][index1]->size()-1;

       for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
       {
         //if (ranges["L-R"][index1]->at(iEnergyBin)<0) continue;
         double  index2( 10000000*iEnergyBin+index1 );
         if(!p1_deltaT_vs_energyRatio[index2]) continue;

         std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));


         c = new TCanvas(Form("c_deltaT_vs_energyRatio_%s",labelLR_energyBin.c_str()),Form("c_deltaT_vs_energyRatio_%s",labelLR_energyBin.c_str()));

         prof = p1_deltaT_vs_energyRatio[index2];
         prof -> SetTitle(Form(";energy_{right} / energy_{left};#Deltat [ps]"));
         prof -> GetYaxis() -> SetRangeUser(CTRMeans[index2]-2.*CTRSigmas[index2],CTRMeans[index2]+2.*CTRSigmas[index2]);
         prof -> Draw("");

         latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.1f V, th. = %d DAC}",iBar,Vov,int(vth1)));
         latex -> SetNDC();
         latex -> SetTextFont(42);
         latex -> SetTextSize(0.04);
         latex -> SetTextColor(kRed);
         latex -> Draw("same");

         float fitXMin = fitFunc_energyRatio[index2]->GetParameter(1) - 2.*fitFunc_energyRatio[index2]->GetParameter(2);
         float fitXMax = fitFunc_energyRatio[index2]->GetParameter(1) + 2.*fitFunc_energyRatio[index2]->GetParameter(2);

         fitFunc_energyRatioCorr[index2] = new TF1(Form("fitFunc_energyRatioCorr_%s",labelLR_energyBin.c_str()),"pol1",fitXMin,fitXMax);
         prof -> Fit(fitFunc_energyRatioCorr[index2],"QRS+");
         fitFunc_energyRatioCorr[index2] -> SetLineColor(kRed);
         fitFunc_energyRatioCorr[index2] -> SetLineWidth(2);
         fitFunc_energyRatioCorr[index2] -> Draw("same");

         c -> Print(Form("%s/CTR/c_deltaT_vs_energyRatio__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
         c -> Print(Form("%s/CTR/c_deltaT_vs_energyRatio__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
         delete c;
       }
     }
   }



   //------------------------
   //--- 4th loop over events

   gStyle->SetOptFit(1111);

   for(auto mapIt : trees)
   {
     ModuleEventClass* anEvent = new ModuleEventClass();
     mapIt.second -> SetBranchAddress("event",&anEvent);

     int nEntries = mapIt.second->GetEntries();
     for(int entry = 0; entry < nEntries; ++entry)
     {

       if( entry%100000 == 0 ) std::cout << ">>> 4th loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
       mapIt.second -> GetEntry(entry);

       int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

       if( !accept[index1][entry] ) continue;

       int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;
       double  index2( 10000000*energyBinAverage+index1 );     

       long long deltaT = anEvent->timeR - anEvent->timeL;

       float energyRatioCorr = 1.;

       if( !fitFunc_energyRatioCorr[index2] )	continue;


       energyRatioCorr = fitFunc_energyRatioCorr[index2]->Eval(anEvent->energyR/anEvent->energyL) -
         fitFunc_energyRatioCorr[index2]->Eval(fitFunc_energyRatio[index2]->GetParameter(1));


       if( h1_deltaT_energyRatioCorr[index2] == NULL )
       {

         std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.1f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));

         h1_deltaT_energyRatioCorr[index2] = new TH1F(Form("h1_deltaT_energyRatioCorr_%s",labelLR_energyBin.c_str()),"",1000,-5000.,5000.);

       }

       h1_deltaT_energyRatioCorr[index2] -> Fill( deltaT - energyRatioCorr );

     }
     std::cout << std::endl;
   }


   double  theIndex2;
   float enBin;
   float theDeltaTMean;
   float theEntriesCTR;
   float errTimeRes;
   float timeResEffSigma;


   //------------------
   //--- draw 4th plots
   for(auto stepLabel : stepLabels)
   {
     float Vov = map_Vovs[stepLabel];
     float vth1 = map_ths[stepLabel];
     
     for(int iBar = 0; iBar < 16; ++iBar)
     {
       int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );
       
       std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
       
       if( !ranges["L-R"][index1] ) continue;
       
       int nEnergyBins = ranges["L-R"][index1]->size()-1;
       
       
       for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
       {
         double  index2( 10000000*iEnergyBin + index1 );
         
         outTrees2[index2] = new TTree(Form("dataRes_bar%02dL-R_Vov%.1f_th%02.0f_enBin%02d",iBar,Vov,vth1,iEnergyBin),Form("dataRes_bar%02dL-R_Vov%.1f_th%02.0f_enBin%02d",iBar,Vov,vth1,iEnergyBin));
         outTrees2[index2] -> Branch("energyBin", &enBin);
         outTrees2[index2] -> Branch("timeResolution",&timeRes);
         outTrees2[index2] -> Branch("effSigma",&timeResEffSigma);
         outTrees2[index2] -> Branch("errTimeResolution",&errTimeRes);
         outTrees2[index2] -> Branch("indexID2",&theIndex2);
         outTrees2[index2] -> Branch("DeltaTMean",&theDeltaTMean);
         outTrees2[index2] -> Branch("EntriesCTR",&theEntriesCTR);


         if(!h1_deltaT_energyRatioCorr[index2]) continue;

         std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));


         c = new TCanvas(Form("c_deltaT_energyRatioCorr_%s",labelLR_energyBin.c_str()),Form("c_deltaT_energyRatioCorr_%s",labelLR_energyBin.c_str()));

         // -- energy corr deltaT
         histo = h1_deltaT_energyRatioCorr[index2];

         histo -> SetTitle(Form(";energy-corrected #Deltat [ps];entries"));
         histo -> SetLineWidth(2);
         histo -> SetLineColor(kBlue);
         histo -> SetMarkerColor(kBlue);
         histo -> RebinX(2);

         histo -> Draw("");


         FindSmallestInterval(vals,histo,0.68);
         float min = vals[4];
         float max = vals[5];
         float delta = max-min;
         float sigma = 0.5*delta;
         float effSigma = sigma;


         float fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
         float fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;

         TF1* fitFunc = new TF1(Form("fitFunc_energyCorr_%s",labelLR_energyBin.c_str()),"gaus",-5000, 5000);
         fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
         // non mi prende il range di fit se non faccio SetRange a mano...
         fitFunc -> SetRange(fitXMin, fitXMax);
         histo -> Fit(fitFunc,"QNRSL+","", fitXMin, fitXMax);
         fitFunc -> SetRange(fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
         histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
         fitFunc -> SetRange(fitFunc->GetParameter(1)-1.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+1.5*fitFunc->GetParameter(2));
         histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-1.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+1.5*fitFunc->GetParameter(2));

         fitFunc -> SetLineColor(kBlue+1);
         fitFunc -> SetLineWidth(3);
         fitFunc -> Draw("same");         
 
         outFile -> cd();
         histo -> Write();
 	
         enBin = energyBin["L-R"][index1][iEnergyBin];
         theIndex2 = index2;
         timeRes = fabs(fitFunc->GetParameter(2));
         timeResEffSigma = effSigma;
	
         theDeltaTMean = fitFunc->GetParameter(1);
         theEntriesCTR = fitFunc->GetParameter(0);
         errTimeRes = fabs(fitFunc->GetParError(2));
	
         int contr = 0;
         if(!source.compare("Na22") && fitFunc->GetParameter(0)<20){
           timeRes = 0;
           effSigma = 0;
           contr = 1;
         }
                
	outTrees2[index2]->Fill();
        
        if (!source.compare("Laser")) histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-10.*fitFunc->GetParameter(2),
                                                                          fitFunc->GetParameter(1)+10.*fitFunc->GetParameter(2));
	
	histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
	
        
        
        latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corr.}^{eff} = %.0f ps}{#sigma_{corr.}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc->GetParameter(2))));
				if(contr==1) latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{corr.}^{eff} = %.0f ps}{#sigma_{corr.}^{gaus} = %.0f ps}",effSigma,timeRes));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kBlue);
        latex -> Draw("same");
        
	// -- raw delta T
        c->cd();
        histo = h1_deltaT[index2];
        histo -> SetLineWidth(2);
        histo -> SetLineColor(kRed);
        histo -> SetMarkerColor(kRed);
        histo -> RebinX(2);

        histo -> Draw("same");


        FindSmallestInterval(vals,histo,0.68);
        min = vals[4];
        max = vals[5];
        delta = max-min;
        sigma = 0.5*delta;
        effSigma = sigma;
        
        fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
        fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;
                
        fitFunc = new TF1(Form("fitFunc_%s",labelLR_energyBin.c_str()),"gaus",-5000,5000);
        fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
        fitFunc -> SetRange(fitXMin, fitXMax);
        histo -> Fit(fitFunc,"QNRSL+","", fitXMin, fitXMax);
        fitFunc -> SetRange(fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
        histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-fitFunc->GetParameter(2),fitFunc->GetParameter(1)+fitFunc->GetParameter(2));
        fitFunc -> SetRange(fitFunc->GetParameter(1)-1.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+1.5*fitFunc->GetParameter(2));
        histo -> Fit(fitFunc,"QNRSL+","",fitFunc->GetParameter(1)-1.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+1.5*fitFunc->GetParameter(2));
                
        fitFunc -> SetLineColor(kRed+1);
        fitFunc -> SetLineWidth(3);
        fitFunc -> Draw("same");
        
        latex = new TLatex(0.20,0.85,Form("#splitline{#sigma_{raw}^{eff} = %.0f ps}{#sigma_{raw}^{gaus} = %.0f ps}",effSigma,fabs(fitFunc->GetParameter(2))));
        latex -> SetNDC();
        latex -> SetTextFont(42);
        latex -> SetTextSize(0.04);
        latex -> SetTextColor(kRed);
        latex -> Draw("same");
        outFile -> cd();
        histo -> Write();
        c -> Print(Form("%s/CTR_energyRatioCorr/c_deltaT_energyRatioCorr__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
        c -> Print(Form("%s/CTR_energyRatioCorr/c_deltaT_energyRatioCorr__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
        delete c;
      }
    }
  }
    
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
