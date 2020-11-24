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
  return par[0] * PDE_vs_OV_HDR2(xx-par[1]) * Gain_vs_OV_HDR2(xx-par[1]);
}



int main(int argc, char** argv)
{
  setTDRStyle();
  
  
  if( argc < 2 )
  {
    std::cout << ">>> analyzeSlewRate::usage:   " << argv[0] << " configFile.cfg" << std::endl;
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
  float timeBins = opts.GetOpt<float>("Plots.timeBins");
  float timeMin = opts.GetOpt<float>("Plots.timeMin");
  float timeMax = opts.GetOpt<float>("Plots.timeMax");
  
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
  std::vector<unsigned int> srChannels = opts.GetOpt<std::vector<unsigned int> >("Channels.srChannels");
  
  
  
  //--- define branches
  float step1;
  float step2;
  float energy[256];
  float qfine[256];
  float tot[256];
  long long int time[256];
  
  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("step1", 1); tree -> SetBranchAddress("step1", &step1);
  tree -> SetBranchStatus("step2", 1); tree -> SetBranchAddress("step2", &step2);
  tree -> SetBranchStatus("energy",1); tree -> SetBranchAddress("energy",energy);
  tree -> SetBranchStatus("qfine", 1); tree -> SetBranchAddress("qfine",  qfine);
  tree -> SetBranchStatus("tot",   1); tree -> SetBranchAddress("tot",      tot);
  tree -> SetBranchStatus("time",  1); tree -> SetBranchAddress("time",    time);
  
  
  
  //--- define histograms
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileName");
  TFile* outFile = TFile::Open(Form("%s",outFileName.c_str()),"RECREATE");
  outFile -> cd();
  
  
  std::map<std::string,TH1F*> h1_energy;
  std::map<std::string,TH1F*> h1_deltaT;
  
  
  
  //--- loop over events
  int nEntries = tree->GetEntries();
  if( maxEntries > nEntries ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    tree -> GetEntry(entry);
    if( entry%100000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    
    float Vov = step1;
    float vth1 = float(int(step2/10000)-1);;
    
    for(auto ch: channels)
    {
      std::string label(Form("ch%d__Vov_%.1f",ch,Vov));
      if( !h1_energy[label] )
        h1_energy[label] = new TH1F(Form("h1_energy_%s",label.c_str()),"",energyBins,energyMin,energyMax);
      
      if( qfine[ch] > 13 )
      {
        h1_energy[label] -> Fill( energy[ch] );
      }
    } 
    
    for(auto ch: srChannels)
    {
      int refCh = opts.GetOpt<int>(Form("ch%d.refChannel",ch));
      
      std::string label(Form("ch%d__Vov_%.1f__vth1_%02d",ch,Vov,int(vth1)));
      if( !h1_deltaT[label] )
        h1_deltaT[label] = new TH1F(Form("h1_deltaT_%s",label.c_str()),"",timeBins,timeMin,timeMax);
      
      if( qfine[ch] > 13 )
      {
        h1_deltaT[label] -> Fill( time[ch]-time[refCh] );
      }
    } 
  }
  std::cout << std::endl;
  
  
  
  //--- draw plots
  for(auto mapIt : h1_energy)
  {
    std::string label = mapIt.first;
    
    std::size_t pos = label.find("_");
    std::string ch_str = label.substr(0,pos);
    
    
    TCanvas* c1 = new TCanvas(Form("%s",label.c_str()),Form("%s",label.c_str()));
    c1 -> cd(1);
    
    TH1F* histo = h1_energy[label];
    
    if( histo == NULL ) continue;
    if( histo->GetEntries() < 100 ) continue;
    
    histo -> SetTitle(";energy [a.u.];entries");
    histo -> SetLineColor(kRed);
    histo -> SetLineWidth(2);
    histo -> Draw("hist");
    
    
    c1 -> Print(Form("%s/c1_energy_%s.png",plotDir.c_str(),label.c_str()));
    delete c1;
  }
  
  std::map<int,TGraphErrors*> g_slewrate;
  
  for(auto mapIt : h1_deltaT)
  {
    std::string label = mapIt.first;
    
    std::size_t pos = label.find("_");
    std::string ch_str = label.substr(2,pos);
    int ch = atoi(ch_str.c_str());
    
    pos = label.find("vth1_");
    std::string vth1_str = label.substr(pos+5,pos+7);
    
    if( !g_slewrate[ch] )
      g_slewrate[ch] = new TGraphErrors();
    
    
    TCanvas* c1 = new TCanvas(Form("%s",label.c_str()),Form("%s",label.c_str()));
    c1 -> cd(1);
    
    TH1F* histo = h1_deltaT[label];
    
    if( histo == NULL ) continue;
    if( histo->GetEntries() < 100 ) continue;
    
    histo -> SetTitle(";t-t_{ref} [ps];entries");
    histo -> SetLineColor(kRed);
    histo -> SetLineWidth(2);
    histo -> Draw("hist");
    
    TF1* fitFunc = new TF1(Form("fitFunc_%s",label.c_str()),"gaus(0)",timeMin,timeMax);
    fitFunc -> SetParameters(histo->GetMaximum(),histo->GetMean(),histo->GetRMS());
    histo -> Fit(fitFunc,"QNRSM");
    fitFunc -> Draw("same");
    
    g_slewrate[ch] -> SetPoint(g_slewrate[ch]->GetN(),fitFunc->GetParameter(1),atoi(vth1_str.c_str()));
    g_slewrate[ch] -> SetPointError(g_slewrate[ch]->GetN()-1,fitFunc->GetParError(1),0.);
    
    c1 -> Print(Form("%s/c1_deltaT_%s.png",plotDir.c_str(),label.c_str()));
    delete c1;
  }
  
  
  for(auto mapIt : g_slewrate )
  {
    int ch = mapIt.first;
    TGraphErrors* graph = mapIt.second;
    
    TCanvas* c2 = new TCanvas(Form("c_slewrate_ch%d",ch),Form("c_slewrate_ch%d",ch));
    c2 -> cd();
    
    graph -> Draw("APL");
    
    c2 -> Print(Form("%s/c1_slewrate_ch%d.png",plotDir.c_str(),ch));
  }
  
  
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  
  
  return 0;
}
