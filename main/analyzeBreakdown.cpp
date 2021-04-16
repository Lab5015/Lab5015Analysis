#include "interface/AnalysisUtils.h"
#include "interface/FitUtils.h"
#include "interface/SetTDRStyle.h"
#include "interface/SiPM_HDR2.h"
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



std::string exec(std::string command)
{
  char buffer[128];
  std::string result = "";
  
  // Open pipe to file
  FILE* pipe = popen(command.c_str(), "r");
  if( !pipe ) 
    return "popen failed!";
  
  // read till end of process:
  while( !feof(pipe) )
  {
    // use buffer to read and add to result
    if( fgets(buffer,128,pipe) != NULL )
      result += buffer;
  }
  
  pclose(pipe);
  return result;
}



double amp_vs_ov(double* x, double* par)
{
  double xx = x[0];
  return par[0] * PDE_vs_OV(xx-par[1]) * Gain_vs_OV(xx-par[1]);
}



int main(int argc, char** argv)
{
  setTDRStyle();
  
  
  if( argc < 2 )
  {
    std::cout << ">>> analyzeBreakdown::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }
  
  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  
  
  //--- get plot settings
  float energyBins = opts.GetOpt<float>("Plots.energyBins");
  float energyMin = opts.GetOpt<float>("Plots.energyMin");
  float energyMax = opts.GetOpt<float>("Plots.energyMax");
  float energyRangeMin = opts.GetOpt<float>("Plots.energyRangeMin");
  float energyRangeMax = opts.GetOpt<float>("Plots.energyRangeMax");
  float fitMinOV = opts.GetOpt<float>("Plots.fitMinOV");
  float fitMaxOV = opts.GetOpt<float>("Plots.fitMaxOV");
  
  //--- open files and make the tree chain
  std::string inputRecoDir = opts.GetOpt<std::string>("Input.inputRecoDir");
  std::string inputRawDir = opts.GetOpt<std::string>("Input.inputRawDir");
  std::vector<std::string> recoFileNames = opts.GetOpt<std::vector<std::string> >("Input.recoFileNames");
  std::vector<std::string> rawFileNames  = opts.GetOpt<std::vector<std::string> >("Input.rawFileNames");
  int maxEntries = opts.GetOpt<int>("Input.maxEntries");
  TChain* tree = new TChain("data","data");
  
  float t_temp;
  time_t timesec;
  int ii = 0;
  for(auto recoFileName : recoFileNames)
  {
    std::cout << ">>> Adding flle " << recoFileName << std::endl;
    tree -> Add((inputRecoDir+"/"+recoFileName).c_str());
    
    struct stat t_stat;
    stat((inputRawDir+"/"+rawFileNames[ii]).c_str(), &t_stat);
    struct tm * timeinfo = localtime(&t_stat.st_mtime);
    timesec = mktime(timeinfo);
    std::cout << "Time and date of raw file of run " << rawFileNames[ii] << ": " << asctime(timeinfo) << std::endl;
    
    std::string datetime(Form("%d-%d-%d %02d:%02d:%02d",timeinfo->tm_year+1900,timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_hour,timeinfo->tm_min,timeinfo->tm_sec));
    std::string temp = exec(("ssh pi@100.100.100.5 \'./getTemperature.py "+datetime+"\'"));
    t_temp = atof(temp.c_str());
    
    ++ii;
  }
  
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  system(Form("mkdir -p %s",plotDir.c_str()));
  system(Form("cp %s/../index.php %s",plotDir.c_str(),plotDir.c_str()));  
  
  
  
  std::vector<unsigned int> channels = opts.GetOpt<std::vector<unsigned int> >("Channels.channels");
  int combineChannels = opts.GetOpt<int>("Channels.combineChannels");
  
  

  
  
  //--- define branches
  float step1, step2;
  int channelIdx[128];
  std::vector<float> *qfine = 0;
  std::vector<float> *tot = 0;
  std::vector<float> *energy = 0;
  
  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("step1",  1); tree -> SetBranchAddress("step1",  &step1);
  tree -> SetBranchStatus("step2",  1); tree -> SetBranchAddress("step2",  &step2);
  
  tree -> SetBranchStatus("channelIdx",  1); tree -> SetBranchAddress("channelIdx",  channelIdx);
  
  tree -> SetBranchStatus("qfine",  1); tree -> SetBranchAddress("qfine",   &qfine);  
  tree -> SetBranchStatus("tot",    1); tree -> SetBranchAddress("tot",       &tot);
  tree -> SetBranchStatus("energy", 1); tree -> SetBranchAddress("energy", &energy);

  
  //--- define histograms
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileName");
  TFile* outFile = TFile::Open(Form("%s",outFileName.c_str()),"RECREATE");
  outFile -> cd();
  
  
  std::map<std::string,TH1F*> h1_energy;
  
  
  
  //--- loop over events
  int nEntries = tree->GetEntries();
  if( maxEntries > nEntries ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    tree -> GetEntry(entry);
    if( entry%100000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    
    
    if( combineChannels )
    {
      for(int it = 0; it < channels.size()/2; ++it)
      {
        int ch1 = channels.at(0+it*2);
        int ch2 = channels.at(1+it*2);

        if (channelIdx[ch1] <0 || channelIdx[ch2] <0) continue;
  
 
       std::string label(Form("ch%d-ch%d_ov%.1f",ch1,ch2,step1));
        if( !h1_energy[label] )
          h1_energy[label] = new TH1F(label.c_str(),"",energyBins,energyMin,energyMax);
        h1_energy[label] -> Fill( 0.5*((*energy)[channelIdx[ch1]]+(*energy)[channelIdx[ch2]]) );
      }
    }
    
    else
    {
      for(auto ch: channels)
      {
        std::string label(Form("ch%d_ov%.1f",ch,step1));
        if( !h1_energy[label] )
          h1_energy[label] = new TH1F(label.c_str(),"",energyBins,energyMin,energyMax);
        
        h1_energy[label] -> Fill( (*energy)[channelIdx[ch]] );
       } 
    }
  }
  std::cout << std::endl;
  
  
  
  //--- draw plots
  std::map<std::string,TGraphErrors*> g_peak511;

  if( combineChannels )
  {
    for(int it = 0; it < channels.size()/2; ++it)
    {
      int ch1 = channels.at(0+it*2);
      int ch2 = channels.at(1+it*2);
      std::string label(Form("ch%d-ch%d",ch1,ch2));
      g_peak511[label] = new TGraphErrors();      
    }
  }
  else
  {
    for(auto ch : channels)
    {
      std::string label(Form("ch%d",ch));
      g_peak511[label] = new TGraphErrors();
    }
  }
  
  for(auto mapIt : h1_energy)
  {
    std::string label = mapIt.first;

    std::size_t pos = label.find("_");
    std::string ch_str = label.substr(0,pos);
    
    pos = label.find("ov");
    std::string ov_str = label.substr(pos+2);
    float ov = atof(ov_str.c_str());
    
    
    TCanvas* c1 = new TCanvas(Form("%s",label.c_str()),Form("%s",label.c_str()));
    c1 -> cd(1);
    
    TH1F* histo = h1_energy[label];
    
    if( histo == NULL ) continue;
    if( histo->GetEntries() < 100 ) continue;
    
    histo -> SetTitle(";energy [a.u.];entries");
    histo -> SetLineColor(kRed);
    histo -> SetLineWidth(2);
    histo -> Draw("hist");
    
    
    
    // -- fit 511 keV peak only
    float xMax = FindXMaximum(histo,energyRangeMin,energyRangeMax);
    TF1* f_gaus = new TF1("f_gaus","gaus(0)",xMax-0.20*xMax,xMax+0.30*xMax);
    histo -> Fit(f_gaus,"QNRS");
    f_gaus -> SetLineColor(kBlue);
    f_gaus -> SetLineWidth(1);
    f_gaus -> Draw("same");
    float mean = f_gaus->GetParameter(1);
    float sigma = f_gaus->GetParameter(2);
    TF1* f_gaus2 = new TF1("f_gaus2","gaus(0)",mean-0.5*sigma,mean+sigma);
    histo -> Fit(f_gaus2,"QNRS+");
    f_gaus2 -> SetLineColor(kBlack);
    f_gaus2 -> SetLineWidth(2);
    f_gaus2 -> Draw("same");
    
    if( !g_peak511[ch_str] ) g_peak511[ch_str] = new TGraphErrors();
    g_peak511[ch_str] -> SetPoint(g_peak511[ch_str]->GetN(),ov,f_gaus2->GetParameter(1));
    g_peak511[ch_str] -> SetPointError(g_peak511[ch_str]->GetN()-1,0.,f_gaus2->GetParError(1));
    
    
    c1 -> Print(Form("%s/c1_energy_%s.png",plotDir.c_str(),label.c_str()));
    delete c1;
  }
  
  
  for(auto mapIt : g_peak511)
  {
    std::string ch = mapIt.first;
    
    TCanvas* c1 = new TCanvas(Form("%s",ch.c_str()),Form("%s",ch.c_str()));
    c1 -> cd(1);
    
    TH1F* hpad = (TH1F*)( gPad->DrawFrame(-1.,0.,10.,energyMax) );
    hpad -> SetTitle(";over-voltage [V];511 keV peak amplitude [a.u.]");
    hpad -> Draw();
    
    g_peak511[ch] -> Draw("P,same");
    
    TF1* fitFunc = new TF1("fitFunc",amp_vs_ov,-1.,10.,2);
    fitFunc -> SetNpx(10000);
    fitFunc -> SetParameter(0,20./Gain_vs_OV(3.)/PDE_vs_OV(3.));
    fitFunc -> SetParameter(1,0.);
    g_peak511[ch] -> Fit(fitFunc,"RNS+","",fitMinOV,fitMaxOV);
    fitFunc -> Draw("same");
    
    TLatex* latex = new TLatex(0.25,0.85,Form("#Delta_{OV} = (%.2f #pm %.2f) V",fitFunc->GetParameter(1),fitFunc->GetParError(1)));
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    latex -> SetTextColor(kRed);
    latex -> Draw("same");
    
    c1 -> Print(Form("%s/c1_scan_vs_ov_%s.png",plotDir.c_str(),ch.c_str()));
    delete c1;
    
    
    TGraph* g_nonLinearity = new TGraph();
    for(float jj = 0; jj < 10; jj+=0.1)
    {
      g_nonLinearity -> SetPoint(g_nonLinearity->GetN(),g_peak511[ch]->Eval(jj),fitFunc->Eval(jj)/g_peak511[ch]->Eval(jj));
    }
    outFile -> cd();
    g_nonLinearity -> Write("g_nonLinearity");
  }
  
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  
  
  return 0;
}
