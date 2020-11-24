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



int main(int argc, char** argv)
{
  setTDRStyle();
  
  
  if( argc < 2 )
  {
    std::cout << ">>> analyzeStability::usage:   " << argv[0] << " configFile.cfg" << std::endl;
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
    std::cout << ("ssh pi@100.100.100.5 \'./getTemperature.py "+datetime+"\'") << std::endl;
    t_temp = atof(temp.c_str());
    
    ++ii;
  }
  
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  system(Form("mkdir -p %s",plotDir.c_str()));
  system(Form("cp %s/../index.php %s",plotDir.c_str(),plotDir.c_str()));  
  
  
  
  std::vector<unsigned int> channels = opts.GetOpt<std::vector<unsigned int> >("Channels.channels");
  
  
  
  //--- define branches
  float energy[256];
  float qfine[256];
  float tot[256];
  
  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("energy",1); tree -> SetBranchAddress("energy",energy);
  tree -> SetBranchStatus("qfine", 1); tree -> SetBranchAddress("qfine",  qfine);
  tree -> SetBranchStatus("tot",   1); tree -> SetBranchAddress("tot",      tot);
  
  
  
  //--- define histograms
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileName");
  TFile* outFile = TFile::Open(Form("%s",outFileName.c_str()),"RECREATE");
  outFile -> cd();
  
  
  TTree* outTree = new TTree("data","data");
  unsigned long long int t_time = timesec;
  unsigned int t_channelID;
  float t_peak1;
  float t_peak1Err;
  float t_sigma1;
  float t_sigma1Err;
  outTree -> Branch("time",     &t_time,          "time/l");
  outTree -> Branch("temp",     &t_temp,          "temp/F");
  outTree -> Branch("channelID",&t_channelID,"channelID/i");
  outTree -> Branch("peak1",    &t_peak1,        "peak1/F");
  outTree -> Branch("peak1Err", &t_peak1Err,  "peak1Err/F");
  outTree -> Branch("sigma1",   &t_sigma1,      "sigma1/F");
  outTree -> Branch("sigma1Err",&t_sigma1Err,"sigma1Err/F");
  
  
  std::map<std::string,TH1F*> h1_energy;
  
  for(auto ch: channels)
  {
    std::string label = Form("h1_energy_ch%d",ch);
    
    h1_energy[label] = new TH1F(label.c_str(),"",energyBins,energyMin,energyMax);
  }
  
  
  
  //--- loop over events
  int nEntries = tree->GetEntries();
  if( maxEntries > nEntries ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    tree -> GetEntry(entry);
    if( entry%100000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    
    for(auto ch: channels)
    {
      if( qfine[ch] > 13 )
      {
        h1_energy[Form("h1_energy_ch%d",ch)] -> Fill( energy[ch] );
      }
    } 
  }
  std::cout << std::endl;
  
  
  
  //--- draw plots
  for(auto ch: channels)
  {
    TCanvas* c1 = new TCanvas(Form("c1_ch%d",ch),Form("c1_ch%d",ch));
    c1 -> cd(1);
    
    std::string label = Form("h1_energy_ch%d",ch);
    TH1F* histo = h1_energy[label];
    
    if( histo == NULL ) continue;
    if( histo->GetEntries() < 100 ) continue;
    
    histo -> SetTitle(";energy [a.u.];entries");
    histo -> SetLineColor(kRed);
    histo -> SetLineWidth(2);
    histo -> Draw("hist");
    
    
    
    // -- fit 511 keV peak only
    float xMax = FindXMaximum(histo,energyRangeMin,energyRangeMax);
    TF1* f_gaus = new TF1("f_gaus","gaus(0)",xMax-0.1*xMax,xMax+0.1*xMax);
    histo -> Fit(f_gaus,"QNRS");
    float mean = f_gaus->GetParameter(1);
    float sigma = f_gaus->GetParameter(2);
    TF1* f_gaus2 = new TF1("f_gaus2","gaus(0)",mean-1.5*sigma,mean+2.*sigma);
    histo -> Fit(f_gaus2,"QNRS+");
    f_gaus2 -> SetLineColor(kBlack);
    f_gaus2 -> SetLineWidth(2);
    f_gaus2 -> Draw("same");
    
    
    t_channelID = ch;
    t_peak1 = f_gaus -> GetParameter(1);
    t_peak1Err = f_gaus -> GetParError(1);
    t_sigma1 = f_gaus -> GetParameter(2);
    t_sigma1Err = f_gaus -> GetParError(2);
    outTree -> Fill();
    
    
    c1 -> Print(Form("%s/c1_energy_ch%d.png",plotDir.c_str(),ch));
    delete c1;
  }
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  
  
  return 0;
}
