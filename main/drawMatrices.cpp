#include "interface/AnalysisUtils.h"
#include "interface/FitUtils.h"
#include "interface/SetTDRStyle.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <vector>
#include <time.h>
#include <stdio.h>
#include <sys/stat.h>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"



struct channel
{
  int arrayID;
  int barID;
  int lrID;
};



int main(int argc, char** argv)
{
  setTDRStyle();
  
  
  if( argc < 2 )
  {
    std::cout << ">>> analyzeBars::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }
  
  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  
  
  
  //--- open files and make the tree chain
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  std::vector<std::string> fileBaseNames = opts.GetOpt<std::vector<std::string> >("Input.fileBaseNames");
  int maxEntries = opts.GetOpt<int>("Input.maxEntries");
  TChain* tree = new TChain("data","data");
  
  time_t timesec;
  for(auto fileBaseName : fileBaseNames)
  {
    std::string fileName = Form("%s/%s_t.root",inputDir.c_str(),fileBaseName.c_str());
    std::cout << ">>> Adding flle " << fileName << std::endl;
    tree -> Add(fileName.c_str());
    
    struct stat t_stat;
    stat(Form("/data/TOFPET2/raw/%s.rawf",fileBaseName.c_str()), &t_stat);
    struct tm * timeinfo = localtime(&t_stat.st_mtime);
    timesec = mktime(timeinfo);
    std::cout << "Time and date of raw file of run " << fileBaseName << ": " << asctime(timeinfo);
  }
  
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  system(Form("mkdir -p %s",plotDir.c_str()));
  system(Form("cp /var/www/html/index.php %s",plotDir.c_str()));
  
  
  
  std::vector<unsigned int> channelMapping1 = opts.GetOpt<std::vector<unsigned int> >("Channels.channelMapping1");
  std::vector<unsigned int> channelMapping2 = opts.GetOpt<std::vector<unsigned int> >("Channels.channelMapping2");
  std::vector<unsigned int> channelMapping;
  
  std::map<int,channel> channels_key;
  std::vector<std::string> arrays = opts.GetOpt<std::vector<std::string> >("Channels.arrays");
  for(auto array : arrays)
  {
    unsigned int FEBDPort = opts.GetOpt<unsigned int>(Form("%s.FEBDPort",array.c_str()));
    std::string connectorID = opts.GetOpt<std::string>(Form("%s.connectorID",array.c_str()));
    int arrayID = opts.GetOpt<unsigned int>(Form("%s.arrayID",array.c_str()));
    int offset = 128*FEBDPort;
    
    if( connectorID == "S1_1" || connectorID == "S2_1" )
    {
      channelMapping = channelMapping1;
    }
    else if( connectorID == "S1_2" || connectorID == "S2_2" )
    {
      channelMapping = channelMapping2;
    }
    else continue;
    
    if( connectorID == "S2_1" || connectorID == "S2_2" )
    {
      offset += 64;
    }
    
    for(int ii = 0; ii < 32; ++ii)
    {
      int chID = channelMapping[ii] + offset;
      channels_key[chID] = {arrayID,(ii/2),(ii%2)};
    }
  }
  
  
  //--- define branches
  unsigned int channelID[256];
  unsigned int channelCount[256];
  float energy[256];
  float qfine[256];
  float tot[256];
  
  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("channelID",   1); tree -> SetBranchAddress("channelID",   channelID);
  tree -> SetBranchStatus("channelCount",1); tree -> SetBranchAddress("channelCount",channelCount);
  tree -> SetBranchStatus("energy",      1); tree -> SetBranchAddress("energy",      energy);
  tree -> SetBranchStatus("qfine",       1); tree -> SetBranchAddress("qfine",       qfine);
  tree -> SetBranchStatus("tot",         1); tree -> SetBranchAddress("tot",         tot);
  
  
  
  //--- define histograms
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileName");
  TFile* outFile = TFile::Open(Form("%s",outFileName.c_str()),"RECREATE");
  outFile -> cd();
  
  
  // TTree* outTree = new TTree("data","data");
  // unsigned long long int t_time = timesec;
  // unsigned int t_arrayID;
  // unsigned int t_barID;
  // unsigned int t_lrID;
  // float t_peak1;
  // float t_peak1Err;
  // float t_sigma1;
  // float t_sigma1Err;
  // outTree -> Branch("time",     &t_time,          "time/l");
  // outTree -> Branch("arrayID",  &t_arrayID,    "arrayID/i");
  // outTree -> Branch("barID",    &t_barID,        "barID/i");
  // outTree -> Branch("lrID",     &t_lrID,          "lrID/i");
  // outTree -> Branch("peak1",    &t_peak1,        "peak1/F");
  // outTree -> Branch("peak1Err", &t_peak1Err,  "peak1Err/F");
  // outTree -> Branch("sigma1",   &t_sigma1,      "sigma1/F");
  // outTree -> Branch("sigma1Err",&t_sigma1Err,"sigma1Err/F");
  
  
  std::map<std::string,TH1F*> h1_energy;
  std::map<std::string,TH2F*> h2_energy_LRCorr;
  
  for(auto mapIt: channels_key)
  {
    if( mapIt.second.lrID != 0 ) continue;
    
    int arrayID = mapIt.second.arrayID;
    int barID = mapIt.second.barID;
    
    std::string label_L = Form("h1_energy_array%d_bar%02d_L",arrayID,barID);
    std::string label_R = Form("h1_energy_array%d_bar%02d_R",arrayID,barID);
    
    h1_energy[label_L] = new TH1F(label_L.c_str(),"",300,0.,30.);
    h1_energy[label_R] = new TH1F(label_R.c_str(),"",300,0.,30.);
    
    std::string label_LR = Form("h2_energy_array%d_bar%02d_LRCorr",arrayID,barID);
    h2_energy_LRCorr[label_LR] = new TH2F(label_LR.c_str(),"",300,0.,30.,300,0.,30.);
  }
  
  
  
  //--- loop over events
  int nEntries = tree->GetEntries();
  if( maxEntries > nEntries ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    tree -> GetEntry(entry);
    if( entry%100000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    
    for(auto mapIt: channels_key)
    {
      int chID    = mapIt.first;
      int arrayID = mapIt.second.arrayID;
      int barID   = mapIt.second.barID;
      int lrID    = mapIt.second.lrID;
      
      if( channelCount[chID] == 1 && qfine[chID] > 13 )
      {
        if( lrID == 0 ) h1_energy[Form("h1_energy_array%d_bar%02d_L",arrayID,barID)] -> Fill( energy[chID] );
        else            h1_energy[Form("h1_energy_array%d_bar%02d_R",arrayID,barID)] -> Fill( energy[chID] );
        
        if( lrID == 0 )
        {
          for(auto mapIt2: channels_key)
          {
            if( mapIt2.second.arrayID == arrayID &&
                mapIt2.second.barID == barID &&
                mapIt2.second.lrID == 1 )
            {
              int chID2 = mapIt2.first;
              
              h2_energy_LRCorr[Form("h2_energy_array%d_bar%02d_LRCorr",arrayID,barID)] -> Fill( energy[chID],energy[chID2] );
            }
          }
        }
        
      }
    } 
  }
  std::cout << std::endl;
  
  
  
  //--- draw plots
  for(auto array : arrays)
  {
    TGraphErrors* g_511keV_L = new TGraphErrors();
    TGraphErrors* g_511keV_R = new TGraphErrors();
    
    
    for(int barIt = 0; barIt < 16; ++barIt)
    {
      // t_arrayID = arrayIt;
      // t_barID = barIt;
      
      TCanvas* c1 = new TCanvas("c1","c1",1200,600);
      c1 -> Divide(2,1);
      
      c1 -> cd(1);
      
      std::string label_L = Form("h1_energy_%s_bar%02d_L", array.c_str(),barIt);
      std::string label_R = Form("h1_energy_%s_bar%02d_R", array.c_str(),barIt);
      
      if( h1_energy[label_L] == NULL || h1_energy[label_R] == NULL ) continue;
      
      if( h1_energy[label_L]->Integral() > 0. ) h1_energy[label_L]->Scale(1./h1_energy[label_L]->Integral());
      if( h1_energy[label_R]->Integral() > 0. ) h1_energy[label_R]->Scale(1./h1_energy[label_R]->Integral());
      
      float max = h1_energy[label_L] -> GetMaximum();
      if( h1_energy[label_R]->GetMaximum() > max ) max = h1_energy[label_R]->GetMaximum();
      h1_energy[label_L]->SetMaximum(max);
      
      h1_energy[label_L] -> SetTitle(";energy [a.u.];fraction of events");
      h1_energy[label_L] -> SetLineColor(kRed);
      h1_energy[label_L] -> Draw("hist");
      h1_energy[label_R] -> SetLineColor(kBlue);
      h1_energy[label_R] -> Draw("hist,sames");
      
      
      std::vector<std::string> labels;
      labels.push_back(label_L);
      labels.push_back(label_R);
      
      for(auto label: labels)
      {
        // if( label == label_L ) t_lrID = 0;
        // if( label == label_R ) t_lrID = 1;
        
        
        TH1F* histo = h1_energy[label];
        if( histo->GetEntries() < 100 ) continue;
        
        
        // -- fit 511 keV peak only
        float xMax = FindXMaximum(histo,5.,20.);
        TF1* f_gaus = new TF1("f_gaus","gaus(0)",xMax-0.05*xMax,xMax+0.05*xMax);
        histo -> Fit(f_gaus,"QNRS+");
        f_gaus -> SetLineColor(kBlack);
        f_gaus -> SetLineWidth(3);
        f_gaus -> Draw("same");
        
        if( label == label_L )
        {
          g_511keV_L -> SetPoint(g_511keV_L->GetN(),barIt,f_gaus -> GetParameter(1));
          g_511keV_L -> SetPointError(g_511keV_L->GetN()-1,0,f_gaus -> GetParError(1));
        }
        if( label == label_R )
        {
          g_511keV_R -> SetPoint(g_511keV_R->GetN(),barIt,f_gaus -> GetParameter(1));
          g_511keV_R -> SetPointError(g_511keV_R->GetN()-1,0,f_gaus -> GetParError(1));
        }
        
        // t_peak1 = f_gaus -> GetParameter(1);
        // t_peak1Err = f_gaus -> GetParError(1);
        // t_sigma1 = f_gaus -> GetParameter(2);
        // t_sigma1Err = f_gaus -> GetParError(2);
        
        // outTree -> Fill();
        
        // // -- fit whole 22Na spectrum
        // std::string comptonScattering2 = " 1./(1 + exp( [1] * (x -[2]) ) )";
        // std::string emissionPeak2 = "[3] * exp(-pow(x-[4],2)/(2*[5]*[5]))";
        // std::string backScattering2 = "[15] * exp(-pow(x-[16],2)/(2*[17]*[17]))";
        
        // std::string comptonScattering1 = "[6]/(1 + exp( [7] * (x -[8]) ) )";
        // std::string emissionPeak1 = "[9] * exp(-pow(x-[10],2)/(2*[11]*[11]))";
        // std::string backScattering1 = "[12] * exp(-pow(x-[13],2)/(2*[14]*[14]))";
        
        // std::string megafit ="[0]*("
        //   + comptonScattering2 + " + "
        //   + emissionPeak2      + " + "
        //   + backScattering2    + " + "
        //   + comptonScattering1 + " + "
        //   + emissionPeak1      + " + "
        //   + backScattering1    + ")" ;
        
        // TF1* f_megafit = new TF1("f_megafit", megafit.c_str(),5.,30.);
        
        
        // // *** compton scattering
        // //f_megafit->SetParameter(7, 0.1);
        // f_megafit -> SetParameter(8, peak1);
        
        // //f_megafit->SetParameter(1, 0.1);
        // f_megafit -> SetParameter(2, peak1*2.5*0.7);
        
        
        // // *** emission peak
        // f_megafit->FixParameter(10, peak1);
        // f_megafit->FixParameter(11, peak1sigma);
        
        // f_megafit->SetParameter(4, peak1*2.5*0.7);
        // f_megafit->SetParameter(5, peak1sigma/sqrt(2.5));
        
        
        // // *** backscatter peak
        // f_megafit->FixParameter(12, 0);
        // f_megafit->SetParameter(13, peak1*0.35);
        // f_megafit->SetParameter(14, peak1sigma/sqrt(0.35));
        
        // //f_megafit->SetParameter(15, 0);
        // f_megafit->SetParameter(16, peak1*0.35);
        // f_megafit->SetParameter(17, peak1sigma/sqrt(0.35));
        
        // histo -> Fit(f_megafit,"QNRS+");
        // f_megafit -> SetLineColor(kGreen);
        // f_megafit -> SetLineWidth(3);
        // f_megafit -> Draw("same");
      }
      
      
      c1 -> cd(2);
      
      std::string label_LR = Form("h2_energy_%s_bar%02d_LRCorr", array.c_str(),barIt);
      h2_energy_LRCorr[label_LR] -> Draw("COLZ");
      
      
      c1 -> Print(Form("%s/c1_energy__%s__bar%02d.png",plotDir.c_str(),array.c_str(),barIt));
      delete c1;
    }
    
    TCanvas* c2 = new TCanvas("c2","c2");
    
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,17.,30.) );
    hPad -> SetTitle(";bar ID;photopeak energy [a.u.]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    g_511keV_L -> SetMarkerStyle(20);
    g_511keV_L -> SetMarkerColor(kRed);
    g_511keV_L -> Draw("P,same");
    g_511keV_R -> SetMarkerColor(kBlue);
    g_511keV_R -> SetMarkerStyle(21);
    g_511keV_R -> Draw("P,same");
    
    c2 -> Print(Form("%s/c2_energy__%s.png",plotDir.c_str(),array.c_str()));
  }
  
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  
  
  return 0;
}
