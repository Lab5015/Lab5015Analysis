#include "interface/FitUtils.h"
#include "interface/TreeUtils.h"
#include "interface/SetTDRStyle.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

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






float FindXMaximum(TH1F* histo, const float& xMin, const float& xMax)
{
  float max = -999999999.;
  int binMax = -1;
  for(int bin = 1; bin <= histo->GetNbinsX(); ++bin)
  {
    if( histo->GetBinCenter(bin) < xMin ) continue;
    if( histo->GetBinCenter(bin) > xMax ) continue;
    if( histo->GetBinContent(bin) > max ) { max = histo->GetBinContent(bin); binMax = bin; };
  }
  return histo->GetBinCenter(binMax);
}

struct EventSingle
{
  std::string stepLabel;
  unsigned int ch1;
  unsigned int ch2;
  std::string label1;
  std::string label2;
  std::string label12;
  float qfine1;
  float qfine2;
  float tot1;
  float tot2;
  float energy1;
  float energy2;
  long long time1;
  long long time2;
};

struct EventCoinc
{
  std::string stepLabel;
  unsigned int ch1;
  unsigned int ch2;
  unsigned int ch3;
  unsigned int ch4;
  std::string label1;
  std::string label2;
  std::string label3;
  std::string label4;
  std::string label1234;
  float qfine1;
  float qfine2;
  float qfine3;
  float qfine4;
  float tot1;
  float tot2;
  float tot3;
  float tot4;
  float energy1;
  float energy2;
  float energy3;
  float energy4;
  long long time1;
  long long time2;
  long long time3;
  long long time4;
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
  

  
  //--- get parameters
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  system(Form("mkdir -p %s",plotDir.c_str()));
  system(Form("mkdir -p %s/qfine/",plotDir.c_str()));
  system(Form("mkdir -p %s/tot/",plotDir.c_str()));
  system(Form("mkdir -p %s/totRatio/",plotDir.c_str()));
  system(Form("mkdir -p %s/energy/",plotDir.c_str()));
  system(Form("mkdir -p %s/energyRatio/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR/",plotDir.c_str()));
  system(Form("mkdir -p %s/CTR_energyCorr/",plotDir.c_str()));
  system(Form("cp /afs/cern.ch/user/a/abenagli/public/index.php %s",plotDir.c_str()));
  system(Form("cp /afs/cern.ch/user/a/abenagli/public/index.php %s/qfine/",plotDir.c_str()));
  system(Form("cp /afs/cern.ch/user/a/abenagli/public/index.php %s/tot/",plotDir.c_str()));
  system(Form("cp /afs/cern.ch/user/a/abenagli/public/index.php %s/totRatio/",plotDir.c_str()));
  system(Form("cp /afs/cern.ch/user/a/abenagli/public/index.php %s/energy/",plotDir.c_str()));
  system(Form("cp /afs/cern.ch/user/a/abenagli/public/index.php %s/energyRatio/",plotDir.c_str()));
  system(Form("cp /afs/cern.ch/user/a/abenagli/public/index.php %s/CTR/",plotDir.c_str()));
  system(Form("cp /afs/cern.ch/user/a/abenagli/public/index.php %s/CTR_energyCorr/",plotDir.c_str()));
  
  
  //--- open files and make the tree chain
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  std::string fileBaseName = opts.GetOpt<std::string>("Input.fileBaseName");
  std::string runs = opts.GetOpt<std::string>("Input.runs");
  int maxEntries = opts.GetOpt<int>("Input.maxEntries");
  TChain* tree = new TChain("data","data");
  
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
      for(int run = runMin; run <= runMax; ++run)
	{
	  //std::string fileName = Form("%s/*%04d*.root",inputDir.c_str(),run);
	  std::string fileName = Form("%s/*%s*%d*.root",inputDir.c_str(),fileBaseName.c_str(),run);
	  std::cout << ">>> Adding flle " << fileName << std::endl;
	  tree -> Add(fileName.c_str());
	}
    }
  
  
  //--- define channels
  std::vector<unsigned int> channels = opts.GetOpt<std::vector<unsigned int> >("Channels.channels");
  
  std::vector<int> pairsMode = opts.GetOpt<std::vector<int> >("Channels.pairsMode");
  std::vector<unsigned int> pairs = opts.GetOpt<std::vector<unsigned int> >("Channels.pairs");
  std::vector<std::pair<unsigned int,unsigned int> > pairsVec;
  for(unsigned int ii = 0; ii < pairs.size()/2; ++ii)
  {
    pairsVec.push_back(std::make_pair(pairs.at(0+ii*2),pairs.at(1+ii*2)));
  }
  
  std::vector<unsigned int> bars = opts.GetOpt<std::vector<unsigned int> >("Channels.bars");
  std::vector<std::pair<std::pair<int,int>,std::pair<unsigned int,unsigned int> > > barsVec;
  for(unsigned int ii = 0; ii < bars.size()/4; ++ii)
  {
    barsVec.push_back( std::make_pair(std::make_pair(bars.at(0+ii*4),bars.at(1+ii*4)),std::make_pair(bars.at(2+ii*4),bars.at(3+ii*4))) );
  }
  
  
  //--- get cuts per bar / Vov
  std::map<unsigned int,std::map<float,float> > cut_qfineAcc;
  std::map<unsigned int,std::map<float,float> > cut_totAcc;
  std::map<unsigned int,std::map<float,float> > cut_energyAcc;
  std::map<unsigned int,std::map<float,float> > cut_energyFitMin;
  std::map<unsigned int,std::map<float,float> > cut_energyFitMax;
  for(auto ch :  channels)
  {
    std::vector<float> Vovs          = opts.GetOpt<std::vector<float> >(Form("ch%d.Vovs",int(ch)));
    std::vector<float> qfineMins     = opts.GetOpt<std::vector<float> >(Form("ch%d.qfineMins",int(ch)));
    std::vector<float> totMins       = opts.GetOpt<std::vector<float> >(Form("ch%d.totMins",int(ch)));
    std::vector<float> energyMins    = opts.GetOpt<std::vector<float> >(Form("ch%d.energyMins",int(ch)));
    std::vector<float> energyFitMins = opts.GetOpt<std::vector<float> >(Form("ch%d.energyFitMins",int(ch)));
    std::vector<float> energyFitMaxs = opts.GetOpt<std::vector<float> >(Form("ch%d.energyFitMaxs",int(ch)));
    int iter = 0;
    for(auto Vov : Vovs)
    {
      cut_qfineAcc[ch][Vov]  = qfineMins.at(iter);
      cut_totAcc[ch][Vov]    = totMins.at(iter);
      cut_energyAcc[ch][Vov] = energyMins.at(iter);
      cut_energyFitMin[ch][Vov] = energyFitMins.at(iter);
      cut_energyFitMax[ch][Vov] = energyFitMaxs.at(iter);
      ++iter;
    }
  }
  std::map<std::string,float> cut_energyMin;
  std::map<std::string,float> cut_energyMax;
  
  
  //--- get plot settings
  float tResMin = opts.GetOpt<float>("Plots.tResMin");
  float tResMax = opts.GetOpt<float>("Plots.tResMax");
  
  
  //--- define branches
  float step1, step2;
  unsigned int channelCount[256];
  float tot[256];
  float qfine[256];
  float energy[256];
  long long time[256];
  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("step1",       1); tree -> SetBranchAddress("step1",       &step1);
  tree -> SetBranchStatus("step2",       1); tree -> SetBranchAddress("step2",       &step2);
  tree -> SetBranchStatus("channelCount",1); tree -> SetBranchAddress("channelCount", channelCount);
  tree -> SetBranchStatus("qfine",       1); tree -> SetBranchAddress("qfine",        qfine);
  tree -> SetBranchStatus("tot",         1); tree -> SetBranchAddress("tot",          tot);
  tree -> SetBranchStatus("energy",      1); tree -> SetBranchAddress("energy",       energy);
  tree -> SetBranchStatus("time",        1); tree -> SetBranchAddress("time",         time);
  
  
  //--- define histograms
  std::map<std::string,int> VovLabels;
  std::map<std::string,int> thLabels;
  std::vector<std::string> stepLabels;
  std::map<std::string,float> map_Vovs;
  std::map<std::string,float> map_ths;
  
  std::map<std::string,TH1F*> h1_qfine;
  std::map<std::string,TH2F*> h2_qfine_vs_tot;
  
  std::map<std::string,TH1F*> h1_tot;
  std::map<std::string,TH2F*> h2_tot_corr;
  std::map<std::string,TH1F*> h1_totRatio;
  
  std::map<std::string,TH1F*> h1_energy;
  std::map<std::string,TH1F*> h1_energy_cut;
  std::map<std::string,TH2F*> h2_energy_corr;
  std::map<std::string,TH1F*> h1_energyRatio;
  
  std::map<std::string,TH1F*> h1_deltaT_raw;
  std::map<std::string,TH1F*> h1_deltaT;
  std::map<std::string,TProfile*> p1_deltaT_vs_energyRatio;

  std::map<std::string,TH1F*> h1_deltaT_energyCorr;  
  
  
  
  
  //------------------------
  //--- 1st loop over events
  std::map<std::string,std::vector<EventSingle> > eventsSingle;
  std::map<std::string,std::vector<EventCoinc> > eventsCoinc;
  
  int nEntries = tree->GetEntries();
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    tree -> GetEntry(entry);
    if( entry%10000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    
    float vth1 = float(int(step2/10000)-1);;
    // float vth2 = float(int((step2-10000*vth1)/100)-1);
    // float vthe = float(int((step2-10000*vth1-step2-100*vth2)/1)-1);
    
    std::string VovLabel = Form("Vov%.1f",step1);
    std::string thLabel = Form("th%02.0f",vth1);
    std::string stepLabel = Form("Vov%.1f_th%02.0f",step1,vth1);
    
        
    //--- create histograms, if needed
    for(auto ch : channels)
    {
      std::string label = Form("ch%d_%s",ch,stepLabel.c_str());
      
      if( h1_tot[label] == NULL )
      {
        h1_qfine[label] = new TH1F(Form("h1_qfine_%s",label.c_str()),"",512,-0.5,511.5);
        h2_qfine_vs_tot[label] = new TH2F(Form("h2_qfine_vs_tot_%s",label.c_str()),"",100,0.,500,256,-0.5,255.5);
        
        h1_tot[label] = new TH1F(Form("h1_tot_%s",label.c_str()),"",2000,0.,1000.);
        
        h1_energy[label] = new TH1F(Form("h1_energy_%s",label.c_str()),"",1000,-10.,40.);
        h1_energy_cut[label] = new TH1F(Form("h1_energy_cut_%s",label.c_str()),"",1000,-10.,40.);
        
        VovLabels[VovLabel] += 1;
        thLabels[thLabel] += 1;
        stepLabels.push_back(stepLabel);
        map_Vovs[stepLabel] = step1;
        map_ths[stepLabel] = vth1;
      }
    }
    for(auto pair : pairsVec)
    {
      unsigned int ch1 = pair.first;
      unsigned int ch2 = pair.second;
      std::string label12 = Form("ch%d-ch%d_%s",ch1,ch2,stepLabel.c_str());
      
      if( h2_tot_corr[label12] == NULL )
      {
        h2_tot_corr[label12] = new TH2F(Form("h2_tot_corr_%s",label12.c_str()),"",100,0.,500.,100,0.,500.);
        h1_totRatio[label12] = new TH1F(Form("h1_totRatio_%s",label12.c_str()),"",1000,0.,5.);
        
        h2_energy_corr[label12] = new TH2F(Form("h2_energy_corr_%s",label12.c_str()),"",200,0.,50.,200,0.,50.);
        h1_energyRatio[label12] = new TH1F(Form("h1_energyRatio_%s",label12.c_str()),"",1000,0.,5.);
        
        h1_deltaT_raw[label12] = new TH1F(Form("h1_deltaT_raw_%s",label12.c_str()),"",1250,-2500.,2500.);
        h1_deltaT[label12] = new TH1F(Form("h1_deltaT_%s",label12.c_str()),"",1250,-2500,2500.);
        p1_deltaT_vs_energyRatio[label12] = new TProfile(Form("p1_deltaT_vs_energyRatio_%s",label12.c_str()),"",1000,0.,5.);
        
        h1_deltaT_energyCorr[label12] = new TH1F(Form("h1_deltaT_energyCorr_%s",label12.c_str()),"",1250,-2500.,2500.);
      }
    }
    for(auto bar : barsVec)
    {
      unsigned int ch1 = bar.first.first;
      unsigned int ch2 = bar.first.second;
      unsigned int ch3 = bar.second.first;
      unsigned int ch4 = bar.second.second;
      
      std::string label1234 = Form("ch%d+ch%d-ch%d+ch%d_%s",ch1,ch2,ch3,ch4,stepLabel.c_str());
      
      if( h1_deltaT[label1234] == NULL )
      { 
        h1_totRatio[label1234] = new TH1F(Form("h1_totRatio_%s",label1234.c_str()),"",1000,0.,5.);
        h2_tot_corr[label1234] = new TH2F(Form("h2_tot_corr_%s",label1234.c_str()),"",100,0.,500.,100,0.,500.);
        
        h1_energyRatio[label1234] = new TH1F(Form("h1_energyRatio_%s",label1234.c_str()),"",1000,0.,5.);
        h2_energy_corr[label1234] = new TH2F(Form("h2_energy_corr_%s",label1234.c_str()),"",200,0.,50.,200,0.,50.);
        
        h1_deltaT_raw[label1234] = new TH1F(Form("h1_deltaT_raw_%s",label1234.c_str()),"",1000,-5000.,5000.);
        h1_deltaT[label1234] = new TH1F(Form("h1_deltaT_%s",label1234.c_str()),"",250,-5000.,5000.);
        p1_deltaT_vs_energyRatio[label1234] = new TProfile(Form("p1_deltaT_vs_energyRatio_%s",label1234.c_str()),"",1000,0.,5.);
        
        h1_deltaT_energyCorr[label1234] = new TH1F(Form("h1_deltaT_energyCorr_%s",label1234.c_str()),"",250,-5000.,5000.);
      }
    }
    
    
    //--- fill histograms
    for(auto ch1 : channels)
    {
      std::string label1 = Form("ch%d_%s",ch1,stepLabel.c_str());
      
      float qfine1  = channelCount[int(ch1)] == 1 ? qfine[int(ch1)]     : -1.;
      float tot1    = channelCount[int(ch1)] == 1 ? tot[int(ch1)]/1000. : -1.;
      float energy1 = channelCount[int(ch1)] == 1 ? energy[int(ch1)]    : -1.;
      
      if( channelCount[int(ch1)] == 1 )
      {
        h1_qfine[label1] -> Fill( qfine1 );
        h2_qfine_vs_tot[label1] -> Fill( tot1,qfine1 );
        h1_tot[label1] -> Fill( tot1 );
        h1_energy[label1] -> Fill( energy1 );
      }
    }
    for(auto pair : pairsVec)
    {
      unsigned int ch1 = pair.first;
      unsigned int ch2 = pair.second;
      std::string label1 = Form("ch%d_%s",ch1,stepLabel.c_str());
      std::string label2 = Form("ch%d_%s",ch2,stepLabel.c_str());
      std::string label12 = Form("ch%d-ch%d_%s",ch1,ch2,stepLabel.c_str());
      
      float qfine1    = channelCount[int(ch1)] == 1 ? qfine[int(ch1)]     : -1.;
      float qfine2    = channelCount[int(ch2)] == 1 ? qfine[int(ch2)]     : -1.;
      float tot1      = channelCount[int(ch1)] == 1 ? tot[int(ch1)]/1000. : -1.;
      float tot2      = channelCount[int(ch2)] == 1 ? tot[int(ch2)]/1000. : -1.;
      float energy1   = channelCount[int(ch1)] == 1 ? energy[int(ch1)]    : -1.;
      float energy2   = channelCount[int(ch2)] == 1 ? energy[int(ch2)]    : -1.;
      long long time1 = channelCount[int(ch1)] == 1 ? time[int(ch1)]      : -1.;
      long long time2 = channelCount[int(ch2)] == 1 ? time[int(ch2)]      : -1.;
      
      if( channelCount[int(ch1)] == 1 && channelCount[int(ch2)] == 1 )
      {
        h2_tot_corr[label12] -> Fill( tot1,tot2 );
        h2_tot_corr[label12] -> Fill( tot1,tot2 );
        h2_energy_corr[label12] -> Fill( energy1,energy2 );
        h2_energy_corr[label12] -> Fill( energy1,energy2 );
        
        EventSingle anEvent;
        anEvent.stepLabel = stepLabel;
        anEvent.ch1 = ch1;
        anEvent.ch2 = ch2;
        anEvent.label1 = label1;
        anEvent.label2 = label2;
        anEvent.label12 = label12;
        anEvent.qfine1 = qfine1;
        anEvent.qfine2 = qfine2;
        anEvent.tot1 = tot1;
        anEvent.tot2 = tot2;
        anEvent.energy1 = energy1;
        anEvent.energy2 = energy2;
        anEvent.time1 = time1;
        anEvent.time2 = time2;
        eventsSingle[label12].push_back(anEvent);
      }
    }
    for(auto bar : barsVec)
    {
      unsigned int ch1 = bar.first.first;
      unsigned int ch2 = bar.first.second;
      std::string label1 = Form("ch%d_%s",ch1,stepLabel.c_str());
      std::string label2 = Form("ch%d_%s",ch2,stepLabel.c_str());
      
      unsigned int ch3 = bar.second.first;
      unsigned int ch4 = bar.second.second;
      std::string label3 = Form("ch%d_%s",ch3,stepLabel.c_str());
      std::string label4 = Form("ch%d_%s",ch4,stepLabel.c_str());
      
      std::string label1234 = Form("ch%d+ch%d-ch%d+ch%d_%s",ch1,ch2,ch3,ch4,stepLabel.c_str());
      
      float qfine1 = channelCount[int(ch1)] == 1 ? qfine[int(ch1)] : -1.;
      float qfine2 = channelCount[int(ch2)] == 1 ? qfine[int(ch2)] : -1.;
      float qfine3 = channelCount[int(ch3)] == 1 ? qfine[int(ch3)] : -1.;
      float qfine4 = channelCount[int(ch4)] == 1 ? qfine[int(ch4)] : -1.;
      float tot1 = channelCount[int(ch1)] == 1 ? tot[int(ch1)]/1000. : -1.;
      float tot2 = channelCount[int(ch2)] == 1 ? tot[int(ch2)]/1000. : -1.;
      float tot3 = channelCount[int(ch3)] == 1 ? tot[int(ch3)]/1000. : -1.;
      float tot4 = channelCount[int(ch4)] == 1 ? tot[int(ch4)]/1000. : -1.;
      float energy1 = channelCount[int(ch1)] == 1 ? energy[int(ch1)] : -1.;
      float energy2 = channelCount[int(ch2)] == 1 ? energy[int(ch2)] : -1.;
      float energy3 = channelCount[int(ch3)] == 1 ? energy[int(ch3)] : -1.;
      float energy4 = channelCount[int(ch4)] == 1 ? energy[int(ch4)] : -1.;
      long long time1 = channelCount[int(ch1)] == 1 ? time[int(ch1)] : -1.;
      long long time2 = channelCount[int(ch2)] == 1 ? time[int(ch2)] : -1.;
      long long time3 = channelCount[int(ch3)] == 1 ? time[int(ch3)] : -1.;
      long long time4 = channelCount[int(ch4)] == 1 ? time[int(ch4)] : -1.;
      
      if( channelCount[int(ch1)] == 1 && channelCount[int(ch2)] == 1 && channelCount[int(ch3)] == 1 && channelCount[int(ch4)] == 1)
      {
        h2_tot_corr[label1234] -> Fill( 0.5*(tot1+tot2),0.5*(tot3+tot4) );
        h2_energy_corr[label1234] -> Fill( 0.5*(energy1+energy2),0.5*(energy3+energy4) );
        
        EventCoinc anEvent;
        anEvent.stepLabel = stepLabel;
        anEvent.ch1 = ch1;
        anEvent.ch2 = ch2;
        anEvent.ch3 = ch3;
        anEvent.ch4 = ch4;
        anEvent.label1 = label1;
        anEvent.label2 = label2;
        anEvent.label3 = label3;
        anEvent.label4 = label4;
        anEvent.label1234 = label1234;
        anEvent.qfine1 = qfine1;
        anEvent.qfine2 = qfine2;
        anEvent.qfine3 = qfine3;
        anEvent.qfine4 = qfine4;
        anEvent.tot1 = tot1;
        anEvent.tot2 = tot2;
        anEvent.tot3 = tot3;
        anEvent.tot4 = tot4;
        anEvent.energy1 = energy1;
        anEvent.energy2 = energy2;
        anEvent.energy3 = energy3;
        anEvent.energy4 = energy4;
        anEvent.time1 = time1;
        anEvent.time2 = time2;
        anEvent.time3 = time3;
        anEvent.time4 = time4;
        eventsCoinc[label1234].push_back(anEvent);
      }
    }
  }
  std::cout << std::endl;
  
  std::vector<std::string>::iterator iter;
  iter = std::unique(stepLabels.begin(),stepLabels.end());
  stepLabels.resize( std::distance(stepLabels.begin(),iter) );  
  
  std::sort(stepLabels.begin(),stepLabels.end());
  
  
  
  
  //------------------
  //--- draw 1st plots
  TCanvas* c;
  float* vals = new float[6];
  TLatex* latex;
  
  std::map<std::string,TGraphErrors*> g_tot_vs_th;
  std::map<std::string,TGraphErrors*> g_tot_vs_Vov;
  std::map<std::string,TGraphErrors*> g_energy_vs_th;
  std::map<std::string,TGraphErrors*> g_energy_vs_Vov;
  
  for(auto stepLabel : stepLabels)
  {
    float Vov = map_Vovs[stepLabel];
    float th = map_ths[stepLabel];
    std::string VovLabel(Form("Vov%.1f",Vov));
    std::string thLabel(Form("th%02.0f",th));
    
    
    //--------------------------------------------------------
    
    
    for(auto ch : channels)
    {
      std::string label(Form("ch%d_%s",int(ch),stepLabel.c_str()));
      
      
      c = new TCanvas(Form("c_qfine_%s",label.c_str()),Form("c_qfine_%s",label.c_str()));
      // gPad -> SetLogy();
      
      h1_qfine[label] -> SetTitle(";Q_{fine} [ADC];entries");
      h1_qfine[label] -> SetLineColor(kRed);
      h1_qfine[label] -> Draw();
      h1_qfine[label] -> GetXaxis() -> SetRangeUser(12,200);
      TLine* line_qfineAcc1 = new TLine(cut_qfineAcc[ch][Vov],h1_qfine[label]->GetMinimum(),cut_qfineAcc[ch][Vov],h1_qfine[label]->GetMaximum());
      line_qfineAcc1 -> SetLineColor(kBlack);
      line_qfineAcc1 -> Draw("same");
      latex = new TLatex(0.65,0.85,Form("ch%d",ch));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
      latex -> Draw("same");      
      c -> Print(Form("%s/qfine/c_qfine_%s.png",plotDir.c_str(),label.c_str()));
      c -> Print(Form("%s/qfine/c_qfine_%s.pdf",plotDir.c_str(),label.c_str()));
      
      
      c = new TCanvas(Form("c_qfine_vs_tot_%s",label.c_str()),Form("c_qfine_vs_tot_%s",label.c_str()));
      gPad -> SetLogz();
      
      h2_qfine_vs_tot[label] -> SetTitle(Form(";ch%d ToT [ns];ch%d Q_{fine} [ADC]",ch,ch));
      h2_qfine_vs_tot[label] -> Draw("colz");
      
      c -> Print(Form("%s/qfine/c_qfine_vs_tot_%s.png",plotDir.c_str(),label.c_str()));
      c -> Print(Form("%s/qfine/c_qfine_vs_tot_%s.pdf",plotDir.c_str(),label.c_str()));
      
      
      c = new TCanvas(Form("c_tot_%s",label.c_str()),Form("c_tot_%s",label.c_str()));
      // gPad -> SetLogy();
      
      h1_tot[label] -> SetTitle(";ToT [ns];entries");
      h1_tot[label] -> SetLineColor(kRed);
      h1_tot[label] -> Draw();
      float max1 = FindXMaximum(h1_tot[label],cut_totAcc[ch][Vov],1000.);
      h1_tot[label] -> GetXaxis() -> SetRangeUser(0.25*max1,2.*max1);
      TF1* fitFunc1 = new TF1("fitFunc1","gaus",max1-0.05*max1,max1+0.05*max1);
      h1_tot[label] -> Fit(fitFunc1,"QNRS+");
      fitFunc1 -> SetLineColor(kBlack);
      fitFunc1 -> SetLineWidth(3);
      fitFunc1 -> Draw("same");
      TLine* line_totAcc1 = new TLine(cut_totAcc[ch][Vov],h1_tot[label]->GetMinimum(),cut_totAcc[ch][Vov],h1_tot[label]->GetMaximum());
      line_totAcc1 -> SetLineColor(kBlack);
      line_totAcc1 -> Draw("same");
      latex = new TLatex(0.65,0.85,Form("ch%d",ch));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
      latex -> Draw("same");      
      c -> Print(Form("%s/tot/c_tot_%s.png",plotDir.c_str(),label.c_str()));
      c -> Print(Form("%s/tot/c_tot_%s.pdf",plotDir.c_str(),label.c_str()));
      
      
      if( g_tot_vs_th[Form("ch%d_%s",ch,VovLabel.c_str())] == NULL )
        g_tot_vs_th[Form("ch%d_%s",ch,VovLabel.c_str())] = new TGraphErrors();
      
      if( g_tot_vs_Vov[Form("ch%d_%s",ch,thLabel.c_str())] == NULL )
        g_tot_vs_Vov[Form("ch%d_%s",ch,thLabel.c_str())] = new TGraphErrors();
      
      g_tot_vs_th[Form("ch%d_%s",ch,VovLabel.c_str())] -> SetPoint(g_tot_vs_th[Form("ch%d_%s",ch,VovLabel.c_str())]->GetN(),th,fitFunc1->GetMaximumX());
      g_tot_vs_th[Form("ch%d_%s",ch,VovLabel.c_str())] -> SetPointError(g_tot_vs_th[Form("ch%d_%s",ch,VovLabel.c_str())]->GetN()-1,0.,0.);
      
      g_tot_vs_Vov[Form("ch%d_%s",ch,thLabel.c_str())] -> SetPoint(g_tot_vs_Vov[Form("ch%d_%s",ch,thLabel.c_str())]->GetN(),Vov,fitFunc1->GetMaximumX());
      g_tot_vs_Vov[Form("ch%d_%s",ch,thLabel.c_str())] -> SetPointError(g_tot_vs_Vov[Form("ch%d_%s",ch,thLabel.c_str())]->GetN()-1,0.,0.);
      
      
      c = new TCanvas(Form("c_energy_%s",label.c_str()),Form("c_energy_%s",label.c_str()));
      // gPad -> SetLogy();
      
      h1_energy[label] -> SetTitle(";energy [a.u.];entries");
      h1_energy[label] -> SetLineColor(kRed);
      h1_energy[label] -> Draw();
      max1 = FindXMaximum(h1_energy[label],cut_energyAcc[ch][Vov],100.);
      h1_energy[label] -> GetXaxis() -> SetRangeUser(0.1*max1,2.5*max1);
      fitFunc1 = new TF1("fitFunc1","gaus",max1-cut_energyFitMin[ch][Vov]*max1,max1+cut_energyFitMax[ch][Vov]*max1);
      h1_energy[label] -> Fit(fitFunc1,"QNRS+");
      fitFunc1 -> SetLineColor(kBlack);
      fitFunc1 -> SetLineWidth(3);
      fitFunc1 -> Draw("same");
      cut_energyMin[Form("ch%d_%s",ch,stepLabel.c_str())] = fitFunc1->GetMaximumX()-cut_energyFitMin[ch][Vov]*fitFunc1->GetMaximumX();
      cut_energyMax[Form("ch%d_%s",ch,stepLabel.c_str())] = fitFunc1->GetMaximumX()+cut_energyFitMax[ch][Vov]*fitFunc1->GetMaximumX();
      TLine* line_energyMin1 = new TLine(cut_energyMin[Form("ch%d_%s",ch,stepLabel.c_str())],h1_energy[label]->GetMinimum(),cut_energyMin[Form("ch%d_%s",ch,stepLabel.c_str())],h1_energy[label]->GetMaximum());
      line_energyMin1 -> SetLineColor(kBlack);
      line_energyMin1 -> Draw("same");
      TLine* line_energyMax1 = new TLine(cut_energyMax[Form("ch%d_%s",ch,stepLabel.c_str())],h1_energy[label]->GetMinimum(),cut_energyMax[Form("ch%d_%s",ch,stepLabel.c_str())],h1_energy[label]->GetMaximum());
      line_energyMax1 -> SetLineColor(kBlack);
      line_energyMax1 -> Draw("same");
      latex = new TLatex(0.65,0.85,Form("ch%d",ch));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
      latex -> Draw("same");      
      c -> Print(Form("%s/energy/c_energy_%s.png",plotDir.c_str(),label.c_str()));
      c -> Print(Form("%s/energy/c_energy_%s.pdf",plotDir.c_str(),label.c_str()));
      
      
      if( g_energy_vs_th[Form("ch%d_%s",ch,VovLabel.c_str())] == NULL )
        g_energy_vs_th[Form("ch%d_%s",ch,VovLabel.c_str())] = new TGraphErrors();
      
      if( g_energy_vs_Vov[Form("ch%d_%s",ch,thLabel.c_str())] == NULL )
        g_energy_vs_Vov[Form("ch%d_%s",ch,thLabel.c_str())] = new TGraphErrors();
      
      g_energy_vs_th[Form("ch%d_%s",ch,VovLabel.c_str())] -> SetPoint(g_energy_vs_th[Form("ch%d_%s",ch,VovLabel.c_str())]->GetN(),th,fitFunc1->GetMaximumX());
      g_energy_vs_th[Form("ch%d_%s",ch,VovLabel.c_str())] -> SetPointError(g_energy_vs_th[Form("ch%d_%s",ch,VovLabel.c_str())]->GetN()-1,0.,0.);
      
      g_energy_vs_Vov[Form("ch%d_%s",ch,thLabel.c_str())] -> SetPoint(g_energy_vs_Vov[Form("ch%d_%s",ch,thLabel.c_str())]->GetN(),Vov,fitFunc1->GetMaximumX());
      g_energy_vs_Vov[Form("ch%d_%s",ch,thLabel.c_str())] -> SetPointError(g_energy_vs_Vov[Form("ch%d_%s",ch,thLabel.c_str())]->GetN()-1,0.,0.);
    }
    
    
    //--------------------------------------------------------
    
    
    for(auto pair : pairsVec)
    {  
      unsigned int ch1 = pair.first;
      unsigned int ch2 = pair.second;
      std::string label12 = Form("ch%d-ch%d_%s",ch1,ch2,stepLabel.c_str());
      
      c = new TCanvas(Form("c_tot_corr_%s",label12.c_str()),Form("c_tot_corr_%s",label12.c_str()));
      gPad -> SetLogz();
      
      h2_tot_corr[label12] -> SetTitle(Form(";ch%d ToT [ns];ch%d ToT [ns]",ch1,ch2));
      h2_tot_corr[label12] -> Draw("colz");
      
      c -> Print(Form("%s/tot/c_tot_corr_%s.png",plotDir.c_str(),label12.c_str()));
      c -> Print(Form("%s/tot/c_tot_corr_%s.pdf",plotDir.c_str(),label12.c_str()));
      
      
      c = new TCanvas(Form("c_energy_corr_%s",label12.c_str()),Form("c_energy_corr_%s",label12.c_str()));
      gPad -> SetLogz();
      
      h2_energy_corr[label12] -> SetTitle(Form(";ch%d energy [a.u.];ch%d energy [a.u.]",ch1,ch2));
      h2_energy_corr[label12] -> Draw("colz");
      
      c -> Print(Form("%s/energy/c_energy_corr_%s.png",plotDir.c_str(),label12.c_str()));
      c -> Print(Form("%s/energy/c_energy_corr_%s.pdf",plotDir.c_str(),label12.c_str()));
    }
    
    
    //--------------------------------------------------------
    
    
    for(auto bar : barsVec)
    {
      unsigned int ch1 = bar.first.first;
      unsigned int ch2 = bar.first.second;
      std::string label1 = Form("ch%d_%s",ch1,stepLabel.c_str());
      std::string label2 = Form("ch%d_%s",ch2,stepLabel.c_str());
      
      unsigned int ch3 = bar.second.first;
      unsigned int ch4 = bar.second.second;
      std::string label3 = Form("ch%d_%s",ch3,stepLabel.c_str());
      std::string label4 = Form("ch%d_%s",ch4,stepLabel.c_str());
      
      std::string label1234 = Form("ch%d+ch%d-ch%d+ch%d_%s",ch1,ch2,ch3,ch4,stepLabel.c_str());
      
      c = new TCanvas(Form("c_tot_corr_%s",label1234.c_str()),Form("c_tot_corr_%s",label1234.c_str()));
      gPad -> SetLogz();
      
      h2_tot_corr[label1234] -> SetTitle(Form(";ch%d+ch%d ToT [ns];ch%d+ch%d ToT [ns]",ch1,ch2,ch3,ch4));
      h2_tot_corr[label1234] -> Draw("colz");
      
      c -> Print(Form("%s/tot/c_tot_corr_%s.png",plotDir.c_str(),label1234.c_str()));
      c -> Print(Form("%s/tot/c_tot_corr_%s.pdf",plotDir.c_str(),label1234.c_str()));
      
      
      c = new TCanvas(Form("c_energy_corr_%s",label1234.c_str()),Form("c_energy_corr_%s",label1234.c_str()));
      gPad -> SetLogz();
      
      h2_energy_corr[label1234] -> SetTitle(Form(";ch%d+ch%d energy [a.u.];ch%d+ch%d energy [a.u.]",ch1,ch2,ch3,ch4));
      h2_energy_corr[label1234] -> Draw("colz");
      
      c -> Print(Form("%s/energy/c_energy_corr_%s.png",plotDir.c_str(),label1234.c_str()));
      c -> Print(Form("%s/energy/c_energy_corr_%s.pdf",plotDir.c_str(),label1234.c_str()));
    }
  }
  
  
  //--------------------------------------------------------
  
  
  for(auto pair : pairsVec)
  {
    unsigned int ch1 = pair.first;
    unsigned int ch2 = pair.second;
    
    c = new TCanvas(Form("c_tot_vs_th_ch%d-ch%d",ch1,ch2),Form("c_tot_vs_th_ch%d-ch%d",ch1,ch2));
    // gPad -> SetLogy();
    
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,500.) );
    hPad -> SetTitle(";threshold [DAC];ToT [ns]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    int iter = 0;
    for(auto mapIt : VovLabels)
    {
      std::string label1(Form("ch%d_%s",ch1,mapIt.first.c_str()));
      std::string label2(Form("ch%d_%s",ch2,mapIt.first.c_str()));
      TGraph* g_tot1 = g_tot_vs_th[label1];
      TGraph* g_tot2 = g_tot_vs_th[label2];
      
      g_tot1 -> SetLineColor(1+iter);
      g_tot1 -> SetMarkerColor(1+iter);
      g_tot1 -> SetMarkerStyle(20);
      g_tot1 -> Draw("PL,same");
      
      g_tot2 -> SetLineColor(1+iter);
      g_tot2 -> SetMarkerColor(1+iter);
      g_tot2 -> SetMarkerStyle(25);
      g_tot2 -> Draw("PL,same");
      
      latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack+iter);
      latex -> Draw("same");
      
      ++iter;
    }
    
    c -> Print(Form("%s/c_tot_vs_th_ch%d-ch%d.png",plotDir.c_str(),ch1,ch2));
    c -> Print(Form("%s/c_tot_vs_th_ch%d-ch%d.pdf",plotDir.c_str(),ch1,ch2));
    
    
    c = new TCanvas(Form("c_tot_vs_Vov_ch%d-ch%d",ch1,ch2),Form("c_tot_vs_Vov_ch%d-ch%d",ch1,ch2));
    // gPad -> SetLogy();
    
    hPad = (TH1F*)( gPad->DrawFrame(0.,0.,10.,500.) );
    hPad -> SetTitle(";V_{ov} [V];ToT [ns]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    iter = 0;
    for(auto mapIt : thLabels)
    {
      std::string label1(Form("ch%d_%s",ch1,mapIt.first.c_str()));
      std::string label2(Form("ch%d_%s",ch2,mapIt.first.c_str()));
      
      TGraph* g_tot1 = g_tot_vs_Vov[label1];
      TGraph* g_tot2 = g_tot_vs_Vov[label2];
      
      g_tot1 -> SetLineColor(1+iter);
      g_tot1 -> SetMarkerColor(1+iter);
      g_tot1 -> SetMarkerStyle(20);
      g_tot1 -> Draw("PL,same");
      
      g_tot2 -> SetLineColor(1+iter);
      g_tot2 -> SetMarkerColor(1+iter);
      g_tot2 -> SetMarkerStyle(25);
      g_tot2 -> Draw("PL,same");
      
      latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack+iter);
      latex -> Draw("same");
      
      ++iter;
    }
    
    c -> Print(Form("%s/c_tot_vs_Vov_ch%d-ch%d.png",plotDir.c_str(),ch1,ch2));
    c -> Print(Form("%s/c_tot_vs_Vov_ch%d-ch%d.pdf",plotDir.c_str(),ch1,ch2));
    
    
    
    c = new TCanvas(Form("c_energy_vs_th_ch%d-ch%d",ch1,ch2),Form("c_energy_vs_th_ch%d-ch%d",ch1,ch2));
    // gPad -> SetLogy();
    
    hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,50.) );
    hPad -> SetTitle(";threshold [DAC];energy [a.u.]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    iter = 0;
    for(auto mapIt : VovLabels)
    {
      std::string label1(Form("ch%d_%s",ch1,mapIt.first.c_str()));
      std::string label2(Form("ch%d_%s",ch2,mapIt.first.c_str()));
      TGraph* g_energy1 = g_energy_vs_th[label1];
      TGraph* g_energy2 = g_energy_vs_th[label2];
      
      g_energy1 -> SetLineColor(1+iter);
      g_energy1 -> SetMarkerColor(1+iter);
      g_energy1 -> SetMarkerStyle(20);
      g_energy1 -> Draw("PL,same");
      
      g_energy2 -> SetLineColor(1+iter);
      g_energy2 -> SetMarkerColor(1+iter);
      g_energy2 -> SetMarkerStyle(25);
      g_energy2 -> Draw("PL,same");
      
      latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack+iter);
      latex -> Draw("same");
      
      ++iter;
    }
    
    c -> Print(Form("%s/c_energy_vs_th_ch%d-ch%d.png",plotDir.c_str(),ch1,ch2));
    c -> Print(Form("%s/c_energy_vs_th_ch%d-ch%d.pdf",plotDir.c_str(),ch1,ch2));
    
    
    c = new TCanvas(Form("c_energy_vs_Vov_ch%d-ch%d",ch1,ch2),Form("c_energy_vs_Vov_ch%d-ch%d",ch1,ch2));
    // gPad -> SetLogy();
    
    hPad = (TH1F*)( gPad->DrawFrame(0.,0.,10.,50.) );
    hPad -> SetTitle(";V_{ov} [V];energy [a.u.]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    iter = 0;
    for(auto mapIt : thLabels)
    {
      std::string label1(Form("ch%d_%s",ch1,mapIt.first.c_str()));
      std::string label2(Form("ch%d_%s",ch2,mapIt.first.c_str()));
      TGraph* g_energy1 = g_energy_vs_Vov[label1];
      TGraph* g_energy2 = g_energy_vs_Vov[label2];
      
      g_energy1 -> SetLineColor(1+iter);
      g_energy1 -> SetMarkerColor(1+iter);
      g_energy1 -> SetMarkerStyle(20);
      g_energy1 -> Draw("PL,same");
      
      g_energy2 -> SetLineColor(1+iter);
      g_energy2 -> SetMarkerColor(1+iter);
      g_energy2 -> SetMarkerStyle(25);
      g_energy2 -> Draw("PL,same");
      
      latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack+iter);
      latex -> Draw("same");
      
      ++iter;
    }
    
    c -> Print(Form("%s/c_energy_vs_Vov_ch%d-ch%d.png",plotDir.c_str(),ch1,ch2));
    c -> Print(Form("%s/c_energy_vs_Vov_ch%d-ch%d.pdf",plotDir.c_str(),ch1,ch2));
  }
  
  
  
  
  //------------------------
  //--- 2nd loop over events
  std::map<std::string,std::vector<EventSingle> > eventsSingle2;
  std::map<std::string,std::vector<EventCoinc> > eventsCoinc2;
  
  for(auto mapIt : eventsSingle)
  {
    std::string label = mapIt.first;
    float Vov = map_Vovs[label];
        
    nEntries = mapIt.second.size();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 2nd loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
      EventSingle anEvent = mapIt.second.at(entry);
      
      if( anEvent.qfine1 < cut_qfineAcc[anEvent.ch1][Vov] ) continue;
      if( anEvent.qfine2 < cut_qfineAcc[anEvent.ch2][Vov] ) continue;
      if( anEvent.tot1 < cut_totAcc[anEvent.ch1][Vov] ) continue;
      if( anEvent.tot2 < cut_totAcc[anEvent.ch2][Vov] ) continue;
      
      if( anEvent.energy1 > cut_energyMin[Form("ch%d_%s",anEvent.ch1,anEvent.stepLabel.c_str())] &&
          anEvent.energy1 < cut_energyMax[Form("ch%d_%s",anEvent.ch1,anEvent.stepLabel.c_str())] )
        h1_energy_cut[anEvent.label1] -> Fill( anEvent.energy1 );
      if( anEvent.energy2 > cut_energyMin[Form("ch%d_%s",anEvent.ch2,anEvent.stepLabel.c_str())] &&
          anEvent.energy2 < cut_energyMax[Form("ch%d_%s",anEvent.ch2,anEvent.stepLabel.c_str())] )
        h1_energy_cut[anEvent.label2] -> Fill( anEvent.energy2 );
      
      if( anEvent.energy1 < cut_energyMin[Form("ch%d_%s",anEvent.ch1,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy2 < cut_energyMin[Form("ch%d_%s",anEvent.ch2,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy1 > cut_energyMax[Form("ch%d_%s",anEvent.ch1,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy2 > cut_energyMax[Form("ch%d_%s",anEvent.ch2,anEvent.stepLabel.c_str())] ) continue;
      
      h1_totRatio[anEvent.label12] -> Fill( anEvent.tot2 / anEvent.tot1 );
      h1_energyRatio[anEvent.label12] -> Fill( anEvent.energy2 / anEvent.energy1 );
      
      h1_deltaT_raw[anEvent.label12] -> Fill( anEvent.time2-anEvent.time1 );
      
      EventSingle anEvent2;
      anEvent2.stepLabel = anEvent.stepLabel;
      anEvent2.ch1 = anEvent.ch1;
      anEvent2.ch2 = anEvent.ch2;
      anEvent2.label1 = anEvent.label1;
      anEvent2.label2 = anEvent.label2;
      anEvent2.label12 = anEvent.label12;
      anEvent2.qfine1 = anEvent.qfine1;
      anEvent2.qfine2 = anEvent.qfine2;
      anEvent2.tot1 = anEvent.tot1;
      anEvent2.tot2 = anEvent.tot2;
      anEvent2.energy1 = anEvent.energy1;
      anEvent2.energy2 = anEvent.energy2;
      anEvent2.time1 = anEvent.time1;
      anEvent2.time2 = anEvent.time2;
      eventsSingle2[anEvent.label12].push_back(anEvent2);
    }
    std::cout << std::endl;
  }
  
  for(auto mapIt : eventsCoinc)
  {
    std::string label = mapIt.first;
    
    nEntries = mapIt.second.size();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 2nd loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
      EventCoinc anEvent = mapIt.second.at(entry);
      
      if( anEvent.energy1 < cut_energyMin[Form("ch%d_%s",anEvent.ch1,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy2 < cut_energyMin[Form("ch%d_%s",anEvent.ch2,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy3 < cut_energyMin[Form("ch%d_%s",anEvent.ch3,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy4 < cut_energyMin[Form("ch%d_%s",anEvent.ch4,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy1 > cut_energyMax[Form("ch%d_%s",anEvent.ch1,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy2 > cut_energyMax[Form("ch%d_%s",anEvent.ch2,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy3 > cut_energyMax[Form("ch%d_%s",anEvent.ch3,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy4 > cut_energyMax[Form("ch%d_%s",anEvent.ch4,anEvent.stepLabel.c_str())] ) continue;
      
      h1_totRatio[anEvent.label1234] -> Fill( (anEvent.tot3+anEvent.tot4) / (anEvent.tot1+anEvent.tot2) );
      
      h1_deltaT_raw[anEvent.label1234] -> Fill( 0.5*(anEvent.time3+anEvent.time4) - 0.5*(anEvent.time1+anEvent.time2) );
      
      EventCoinc anEvent2;
      anEvent2.stepLabel = anEvent.stepLabel;
      anEvent2.ch1 = anEvent.ch1;
      anEvent2.ch2 = anEvent.ch2;
      anEvent2.ch3 = anEvent.ch3;
      anEvent2.ch4 = anEvent.ch4;
      anEvent2.label1 = anEvent.label1;
      anEvent2.label2 = anEvent.label2;
      anEvent2.label3 = anEvent.label3;
      anEvent2.label4 = anEvent.label4;
      anEvent2.label1234 = anEvent.label1234;
      anEvent2.qfine1 = anEvent.qfine1;
      anEvent2.qfine2 = anEvent.qfine2;
      anEvent2.qfine3 = anEvent.qfine3;
      anEvent2.qfine4 = anEvent.qfine4;
      anEvent2.tot1 = anEvent.tot1;
      anEvent2.tot2 = anEvent.tot2;
      anEvent2.tot3 = anEvent.tot3;
      anEvent2.tot4 = anEvent.tot4;
      anEvent2.energy1 = anEvent.energy1;
      anEvent2.energy2 = anEvent.energy2;
      anEvent2.energy3 = anEvent.energy3;
      anEvent2.energy4 = anEvent.energy4;
      anEvent2.time1 = anEvent.time1;
      anEvent2.time2 = anEvent.time2;
      anEvent2.time3 = anEvent.time3;
      anEvent2.time4 = anEvent.time4;
      eventsCoinc2[anEvent.label1234].push_back(anEvent2);
    }
    std::cout << std::endl;
  }
  
  
  
  
  //------------------
  //--- draw 2nd plots
  std::map<std::string,float> CTRMeans;
  std::map<std::string,float> CTRSigmas;
  
  for(auto stepLabel : stepLabels)
  {
    for(auto pair : pairsVec)
    {
      unsigned int ch1 = pair.first;
      unsigned int ch2 = pair.second;
      std::string label1(Form("ch%d_%s",ch1,stepLabel.c_str()));
      std::string label2(Form("ch%d_%s",ch2,stepLabel.c_str()));
      std::string label12 = Form("ch%d-ch%d_%s",ch1,ch2,stepLabel.c_str());
      
      
      //--------------------------------------------------------
      
      
      FindSmallestInterval(vals,h1_deltaT_raw[label12],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      CTRMeans[label12] = mean;
      CTRSigmas[label12] = effSigma;
      
      
      //--------------------------------------------------------
      
      
      c = new TCanvas(Form("c_totRatio_%s",label12.c_str()),Form("c_totRatio_%s",label12.c_str()));
      // gPad -> SetLogy();
      
      h1_totRatio[label12] -> SetTitle(Form(";%d ToT / %d ToT;entries",ch1,ch2));
      h1_totRatio[label12] -> Draw("colz");
      
      TF1* fitFunc = new TF1(Form("fitFunc_totRatio_%s",label12.c_str()),"gaus",
                             h1_totRatio[label12]->GetMean()-1.5*h1_totRatio[label12]->GetRMS(),
                             h1_totRatio[label12]->GetMean()+1.5*h1_totRatio[label12]->GetRMS());
      fitFunc -> SetParameters(1,h1_totRatio[label12]->GetMean(),h1_totRatio[label12]->GetRMS());
      h1_totRatio[label12] -> Fit(fitFunc,"QNRSL+");
      fitFunc -> SetLineColor(kRed);
      fitFunc -> SetLineWidth(1);
      fitFunc -> Draw("same");
      
      latex = new TLatex(0.55,0.85,Form("#sigma = %.1f %%",100.*fitFunc->GetParameter(2)));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
      latex -> Draw("same");
      
      c -> Print(Form("%s/totRatio/c_totRatio_%s.png",plotDir.c_str(),label12.c_str()));
      c -> Print(Form("%s/totRatio/c_totRatio_%s.pdf",plotDir.c_str(),label12.c_str()));
      
      
      c = new TCanvas(Form("c_energyRatio_%s",label12.c_str()),Form("c_energyRatio_%s",label12.c_str()));
      // gPad -> SetLogy();
      
      h1_energyRatio[label12] -> SetTitle(Form(";ch%d energy / ch%d energy;entries",ch1,ch2));
      h1_energyRatio[label12] -> Draw("colz");
      
      fitFunc = new TF1(Form("fitFunc_energyRatio_%s",label12.c_str()),"gaus",
                        h1_energyRatio[label12]->GetMean()-1.5*h1_energyRatio[label12]->GetRMS(),
                        h1_energyRatio[label12]->GetMean()+1.5*h1_energyRatio[label12]->GetRMS());
      fitFunc -> SetParameters(1,h1_energyRatio[label12]->GetMean(),h1_energyRatio[label12]->GetRMS());
      h1_energyRatio[label12] -> Fit(fitFunc,"QNRSL+");
      fitFunc -> SetLineColor(kRed);
      fitFunc -> SetLineWidth(1);
      fitFunc -> Draw("same");
      
      latex = new TLatex(0.55,0.85,Form("#sigma = %.1f %%",100.*fitFunc->GetParameter(2)));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
      latex -> Draw("same");
      
      c -> Print(Form("%s/energyRatio/c_energyRatio_%s.png",plotDir.c_str(),label12.c_str()));
      c -> Print(Form("%s/energyRatio/c_energyRatio_%s.pdf",plotDir.c_str(),label12.c_str()));
    }
    
    
    for(auto bar : barsVec)
    {
      unsigned int ch1 = bar.first.first;
      unsigned int ch2 = bar.first.second;
      std::string label1 = Form("%d_%s",ch1,stepLabel.c_str());
      std::string label2 = Form("%d_%s",ch2,stepLabel.c_str());

      unsigned int ch3 = bar.second.first;
      unsigned int ch4 = bar.second.second;
      std::string label3 = Form("%d_%s",ch3,stepLabel.c_str());
      std::string label4 = Form("%d_%s",ch4,stepLabel.c_str());

      std::string label1234(Form("ch%d+ch%d-ch%d+ch%d_%s",ch1,ch2,ch3,ch4,stepLabel.c_str()));
      
      
      //--------------------------------------------------------
      
      
      FindSmallestInterval(vals,h1_deltaT_raw[label1234],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      CTRMeans[label1234] = mean;
      CTRSigmas[label1234] = effSigma;
      
      
      //--------------------------------------------------------
      
      
      c = new TCanvas(Form("c_totRatio_%s",label1234.c_str()),Form("c_totRatio_%s",label1234.c_str()));
      gPad -> SetLogy();
      
      h1_totRatio[label1234] -> SetTitle(Form(";%d+%d ToT / %d+%d ToT;entries",ch1,ch2,ch3,ch4));
      h1_totRatio[label1234] -> Draw("colz");
      
      c -> Print(Form("%s/totRatio/c_totRatio_%s.png",plotDir.c_str(),label1234.c_str()));
      c -> Print(Form("%s/totRatio/c_totRatio_%s.pdf",plotDir.c_str(),label1234.c_str()));
    }
  }
  
  
  
  
  //------------------------
  //--- 3rd loop over events
  for(auto mapIt : eventsSingle)
  {
    std::string label = mapIt.first;
    
    nEntries = mapIt.second.size();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 3rd loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
      EventSingle anEvent = mapIt.second.at(entry);
      
      if( anEvent.energy1 < cut_energyMin[Form("ch%d_%s",anEvent.ch1,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy2 < cut_energyMin[Form("ch%d_%s",anEvent.ch2,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy1 > cut_energyMax[Form("ch%d_%s",anEvent.ch1,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy2 > cut_energyMax[Form("ch%d_%s",anEvent.ch2,anEvent.stepLabel.c_str())] ) continue;
      
      float timeLow = CTRMeans[anEvent.label12] - 2.* CTRSigmas[anEvent.label12];
      float timeHig = CTRMeans[anEvent.label12] + 2.* CTRSigmas[anEvent.label12];
      long long deltaT = anEvent.time2 - anEvent.time1;
      
      h1_deltaT[anEvent.label12] -> Fill( deltaT );
      if( ( deltaT > timeLow ) &&
          ( deltaT < timeHig ) )
        p1_deltaT_vs_energyRatio[anEvent.label12] -> Fill( anEvent.energy2/anEvent.energy1,anEvent.time2-anEvent.time1 );    
    }
    std::cout << std::endl;
  }
  
  for(auto mapIt : eventsCoinc)
  {
    std::string label = mapIt.first;
    
    nEntries = mapIt.second.size();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 3rd loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
      EventCoinc anEvent = mapIt.second.at(entry);
      
      if( anEvent.energy1 < cut_energyMin[Form("ch%d_%s",anEvent.ch1,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy2 < cut_energyMin[Form("ch%d_%s",anEvent.ch2,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy3 < cut_energyMin[Form("ch%d_%s",anEvent.ch3,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy4 < cut_energyMin[Form("ch%d_%s",anEvent.ch4,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy1 > cut_energyMax[Form("ch%d_%s",anEvent.ch1,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy2 > cut_energyMax[Form("ch%d_%s",anEvent.ch2,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy3 > cut_energyMax[Form("ch%d_%s",anEvent.ch3,anEvent.stepLabel.c_str())] ) continue;
      if( anEvent.energy4 > cut_energyMax[Form("ch%d_%s",anEvent.ch4,anEvent.stepLabel.c_str())] ) continue;      
      
      float timeLow = CTRMeans[anEvent.label1234] - 2.* CTRSigmas[anEvent.label1234];
      float timeHig = CTRMeans[anEvent.label1234] + 2.* CTRSigmas[anEvent.label1234];
      
      long long timeComb = 0.5*(anEvent.time3+anEvent.time4) - 0.5*(anEvent.time1+anEvent.time2);
      h1_deltaT[anEvent.label1234] -> Fill( timeComb );
      if( ( timeComb > timeLow ) &&
          ( timeComb < timeHig ) )
        p1_deltaT_vs_energyRatio[anEvent.label1234] -> Fill( (anEvent.energy3+anEvent.energy4)/(anEvent.energy1+anEvent.energy2),timeComb );
    }
    std::cout << std::endl;
  }
  
  
  
  
  //------------------
  //--- draw 3rd plots
  std::map<std::string,TF1*> fitFunc_energyCorr;
  
  for(auto stepLabel : stepLabels)
  {
    float Vov = map_Vovs[stepLabel];
    float th = map_ths[stepLabel];
    std::string VovLabel(Form("Vov%.1f",Vov));
    std::string thLabel(Form("th%02.0f",th));
    
    for(auto pair : pairsVec)
    {
      unsigned int ch1 = pair.first;
      unsigned int ch2 = pair.second;
      std::string label1(Form("ch%d_%s",ch1,stepLabel.c_str()));
      std::string label2(Form("ch%d_%s",ch2,stepLabel.c_str()));
      std::string label12 = Form("ch%d-ch%d_%s",ch1,ch2,stepLabel.c_str());
      
      
      //--------------------------------------------------------
      
      
      c = new TCanvas(Form("c_deltaT_vs_energyRatio_%s",label12.c_str()),Form("c_deltaT_vs_energyRatio_%s",label12.c_str()));
      
      p1_deltaT_vs_energyRatio[label12] -> SetTitle(Form(";ch%d energy / ch%d energy;#Deltat [ps]",ch2,ch1));
      p1_deltaT_vs_energyRatio[label12] -> GetXaxis() -> SetRangeUser(h1_energyRatio[label12]->GetMean()-3.*h1_energyRatio[label12]->GetRMS(),
                                                                      h1_energyRatio[label12]->GetMean()+3.*h1_energyRatio[label12]->GetRMS());
      p1_deltaT_vs_energyRatio[label12] -> Draw("");
      
      float fitXMin = h1_energyRatio[label12]->GetMean() - 2.*h1_energyRatio[label12]->GetRMS();
      float fitXMax = h1_energyRatio[label12]->GetMean() + 2.*h1_energyRatio[label12]->GetRMS();
      fitFunc_energyCorr[label12] = new TF1(Form("fitFunc_energyCorr_%s",label12.c_str()),"pol4",fitXMin,fitXMax);
      p1_deltaT_vs_energyRatio[label12] -> Fit(fitFunc_energyCorr[label12],"QNRS+");
      fitFunc_energyCorr[label12] -> SetLineColor(kRed);
      fitFunc_energyCorr[label12] -> SetLineWidth(2);
      fitFunc_energyCorr[label12] -> Draw("same");
      
      c -> Print(Form("%s/CTR/c_deltaT_vs_energyRatio_%s.png",plotDir.c_str(),label12.c_str()));
      c -> Print(Form("%s/CTR/c_deltaT_vs_energyRatio_%s.pdf",plotDir.c_str(),label12.c_str()));
    }
    
    
    for(auto bar : barsVec)
    {
      unsigned int ch1 = bar.first.first;
      unsigned int ch2 = bar.first.second;
      std::string label1(Form("ch%d_%s",ch1,stepLabel.c_str()));
      std::string label2(Form("ch%d_%s",ch2,stepLabel.c_str()));

      unsigned int ch3 = bar.second.first;
      unsigned int ch4 = bar.second.second;
      std::string label3(Form("ch%d_%s",ch3,stepLabel.c_str()));
      std::string label4(Form("ch%d_%s",ch4,stepLabel.c_str()));

      std::string label1234(Form("ch%d+ch%d-ch%d+ch%d_%s",ch1,ch2,ch3,ch4,stepLabel.c_str()));
      
      
      //--------------------------------------------------------
      
      
      c = new TCanvas(Form("c_deltaT_vs_energyRatio_%s",label1234.c_str()),Form("c_deltaT_vs_energyRatio_%s",label1234.c_str()));
      
      p1_deltaT_vs_energyRatio[label1234] -> SetTitle(Form(";ch%d+ch%d energy / ch%d+ch%d energy;#Deltat [ps]",ch3,ch4,ch1,ch2));
      p1_deltaT_vs_energyRatio[label1234] -> GetXaxis() -> SetRangeUser(h1_energyRatio[label1234]->GetMean()-3.*h1_energyRatio[label1234]->GetRMS(),
                                                                        h1_energyRatio[label1234]->GetMean()+3.*h1_energyRatio[label1234]->GetRMS());
      p1_deltaT_vs_energyRatio[label1234] -> Draw("");
      
      float fitXMin = h1_energyRatio[label1234]->GetMean() - 2.*h1_energyRatio[label1234]->GetRMS();
      float fitXMax = h1_energyRatio[label1234]->GetMean() + 2.*h1_energyRatio[label1234]->GetRMS();
      fitFunc_energyCorr[label1234] = new TF1(Form("fitFunc_energyCorr_%s",label1234.c_str()),"pol4",fitXMin,fitXMax);
      p1_deltaT_vs_energyRatio[label1234] -> Fit(fitFunc_energyCorr[label1234],"QNRS+");
      fitFunc_energyCorr[label1234] -> SetLineColor(kRed);
      fitFunc_energyCorr[label1234] -> SetLineWidth(2);
      fitFunc_energyCorr[label1234] -> Draw("same");
      
      c -> Print(Form("%s/CTR/c_deltaT_vs_energyRatio_%s.png",plotDir.c_str(),label1234.c_str()));
      c -> Print(Form("%s/CTR/c_deltaT_vs_energyRatio_%s.pdf",plotDir.c_str(),label1234.c_str()));
    }
  }
  
  
  
  
  //------------------------
  //--- 4th loop over events
  for(auto mapIt : eventsSingle2)
  {
    std::string label = mapIt.first;
    
    nEntries = mapIt.second.size();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 4th loop (" << label << "): reading entry " << entry << " / " << nEntries << "\r" << std::flush;
      EventSingle anEvent = mapIt.second.at(entry);
      
      long long deltaT = anEvent.time2-anEvent.time1;
      float energyCorr = fitFunc_energyCorr[label]->Eval(anEvent.energy2/anEvent.energy1) - fitFunc_energyCorr[label]->Eval(h1_energy_cut[anEvent.label2]->GetMean()/h1_energy_cut[anEvent.label1]->GetMean());
      h1_deltaT_energyCorr[label] -> Fill( deltaT - energyCorr );
    }
    std::cout << std::endl;
  }
  
  for(auto mapIt : eventsCoinc2)
  {
    std::string label = mapIt.first;
    
    nEntries = mapIt.second.size();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 4th loop (" << label << "): reading entry " << entry << " / " << nEntries << "\r" << std::flush;
      EventCoinc anEvent = mapIt.second.at(entry);
      
      long long deltaT = 0.5*(anEvent.time3+anEvent.time4) - 0.5*(anEvent.time1+anEvent.time2);
      float energyCorr = fitFunc_energyCorr[label]->Eval((anEvent.energy3+anEvent.energy4)/(anEvent.energy1+anEvent.energy2)) - fitFunc_energyCorr[label]->Eval((h1_energy_cut[anEvent.label3]->GetMean()+h1_energy_cut[anEvent.label4]->GetMean())/(h1_energy_cut[anEvent.label1]->GetMean()+h1_energy_cut[anEvent.label2]->GetMean()));
      h1_deltaT_energyCorr[label] -> Fill( deltaT - energyCorr );
    }
    std::cout << std::endl;
  }
  
  
  
  //------------------
  //--- draw 4th plots
  std::map<std::string,TGraphErrors*> g_tRes_effSigma_vs_th;
  std::map<std::string,TGraphErrors*> g_tRes_gaus_vs_th;
  std::map<std::string,TGraphErrors*> g_tRes_effSigma_vs_Vov;
  std::map<std::string,TGraphErrors*> g_tRes_gaus_vs_Vov;
  
  std::map<std::string,TGraphErrors*> g_tRes_energyCorr_effSigma_vs_th;
  std::map<std::string,TGraphErrors*> g_tRes_energyCorr_gaus_vs_th;
  std::map<std::string,TGraphErrors*> g_tRes_energyCorr_effSigma_vs_Vov;
  std::map<std::string,TGraphErrors*> g_tRes_energyCorr_gaus_vs_Vov;
  
  for(auto stepLabel : stepLabels)
  {
    float Vov = map_Vovs[stepLabel];
    float th = map_ths[stepLabel];
    std::string VovLabel(Form("Vov%.1f",Vov));
    std::string thLabel(Form("th%02.0f",th));
    
    int pairsIt = 0;
    for(auto pair : pairsVec)
    {
      unsigned int ch1 = pair.first;
      unsigned int ch2 = pair.second;
      std::string label1(Form("ch%d_%s",ch1,stepLabel.c_str()));
      std::string label2(Form("ch%d_%s",ch2,stepLabel.c_str()));
      std::string label12 = Form("ch%d-ch%d_%s",ch1,ch2,stepLabel.c_str());
      
      
      //--------------------------------------------------------
      
      c = new TCanvas(Form("c_deltaT_energyCorr_%s",label12.c_str()),Form("c_deltaT_energyCorr_%s",label12.c_str()));
      
      h1_deltaT_energyCorr[label12] -> GetXaxis() -> SetRangeUser(h1_deltaT_energyCorr[label12]->GetMean()-2.*h1_deltaT_energyCorr[label12]->GetRMS(),
                                                                  h1_deltaT_energyCorr[label12]->GetMean()+2.*h1_deltaT_energyCorr[label12]->GetRMS());          
      h1_deltaT_energyCorr[label12] -> SetTitle(Form(";energy-corrected #Deltat [ps];entries"));
      h1_deltaT_energyCorr[label12] -> SetLineWidth(2);
      h1_deltaT_energyCorr[label12] -> SetLineColor(kBlue);
      h1_deltaT_energyCorr[label12] -> SetMarkerColor(kBlue);
      h1_deltaT_energyCorr[label12] -> Draw("");
      
      float fitXMin = CTRMeans[label12] - 1.5*CTRSigmas[label12];
      float fitXMax = CTRMeans[label12] + 1.5*CTRSigmas[label12];
      TF1* fitFunc = new TF1(Form("fitFunc_energyCorr_%s",label12.c_str()),"gaus",fitXMin,fitXMax);
      fitFunc -> SetParameters(1,h1_deltaT_energyCorr[label12]->GetMean(),h1_deltaT_energyCorr[label12]->GetRMS());
      h1_deltaT_energyCorr[label12] -> Fit(fitFunc,"QNRSL+");
      fitFunc -> SetLineColor(kBlue);
      fitFunc -> SetLineWidth(2);
      fitFunc -> Draw("same");
      
      FindSmallestInterval(vals,h1_deltaT_energyCorr[label12],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      
      latex = new TLatex(0.55,0.85,Form("time walk corr. #splitline{#sigma_{CTR}^{eff} = %.1f ps}{#sigma_{CTR}^{gaus} = %.1f ps}",effSigma,fitFunc->GetParameter(2)));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlue);
      latex -> Draw("same");
      
      if( g_tRes_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] == NULL )
      {
        g_tRes_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] = new TGraphErrors();
        g_tRes_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] = new TGraphErrors();
        g_tRes_energyCorr_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] = new TGraphErrors();
        g_tRes_energyCorr_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] = new TGraphErrors();        
      }
      if( g_tRes_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] == NULL )
      {
        g_tRes_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] = new TGraphErrors();
        g_tRes_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] = new TGraphErrors();
        g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] = new TGraphErrors();
        g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] = new TGraphErrors();
      } 
      
      if( pairsMode.at(pairsIt) == 0 )
      {
        g_tRes_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPoint(g_tRes_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN(),th,effSigma);
        g_tRes_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPointError(g_tRes_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN()-1,0.,5.);
        g_tRes_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPoint(g_tRes_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN(),th,fitFunc->GetParameter(2));
        g_tRes_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPointError(g_tRes_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN()-1,0.,fitFunc->GetParError(2));
        
        g_tRes_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPoint(g_tRes_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN(),Vov,effSigma);
        g_tRes_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPointError(g_tRes_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN()-1,0.,5.);
        g_tRes_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPoint(g_tRes_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN(),Vov,fitFunc->GetParameter(2));
        g_tRes_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPointError(g_tRes_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN()-1,0.,fitFunc->GetParError(2));
      }
      if( pairsMode.at(pairsIt) == 1 )
      {
        g_tRes_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPoint(g_tRes_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN(),th,effSigma/sqrt(2));
        g_tRes_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPointError(g_tRes_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN()-1,0.,5.);
        g_tRes_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPoint(g_tRes_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN(),th,fitFunc->GetParameter(2)/sqrt(2));
        g_tRes_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPointError(g_tRes_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN()-1,0.,fitFunc->GetParError(2)/sqrt(2));
        
        g_tRes_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPoint(g_tRes_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN(),Vov,effSigma/sqrt(2));
        g_tRes_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPointError(g_tRes_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN()-1,0.,5.);
        g_tRes_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPoint(g_tRes_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN(),Vov,fitFunc->GetParameter(2)/sqrt(2));
        g_tRes_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPointError(g_tRes_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN()-1,0.,fitFunc->GetParError(2)/sqrt(2));
      }
      if( pairsMode.at(pairsIt) == 2 )
      {
        g_tRes_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPoint(g_tRes_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN(),th,effSigma/2);
        g_tRes_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPointError(g_tRes_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN()-1,0.,5.);
        g_tRes_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPoint(g_tRes_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN(),th,fitFunc->GetParameter(2)/2);
        g_tRes_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPointError(g_tRes_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN()-1,0.,fitFunc->GetParError(2)/2.);
        
        g_tRes_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPoint(g_tRes_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN(),Vov,effSigma/2);
        g_tRes_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPointError(g_tRes_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN()-1,0.,5.);
        g_tRes_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPoint(g_tRes_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN(),Vov,fitFunc->GetParameter(2)/2);
        g_tRes_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPointError(g_tRes_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN()-1,0.,fitFunc->GetParError(2)/2.);
      }
      
      
      
      h1_deltaT[label12] -> SetLineWidth(2);
      h1_deltaT[label12] -> SetLineColor(kRed);
      h1_deltaT[label12] -> SetMarkerColor(kRed);
      h1_deltaT[label12] -> Draw("same");
      
      fitXMin = CTRMeans[label12] - 1.5*CTRSigmas[label12];
      fitXMax = CTRMeans[label12] + 1.5*CTRSigmas[label12];
      fitFunc = new TF1(Form("fitFunc_%s",label12.c_str()),"gaus",fitXMin,fitXMax);
      fitFunc -> SetParameters(1,h1_deltaT[label12]->GetMean(),h1_deltaT[label12]->GetRMS());
      h1_deltaT[label12] -> Fit(fitFunc,"QNRSL+");
      fitFunc -> SetLineColor(kRed);
      fitFunc -> SetLineWidth(2);
      fitFunc -> Draw("same");
      
      FindSmallestInterval(vals,h1_deltaT[label12],0.68);
      mean = vals[0];
      min = vals[4];
      max = vals[5];
      delta = max-min;
      sigma = 0.5*delta;
      effSigma = sigma;
      
      h1_deltaT_energyCorr[label12] -> GetXaxis() -> SetRangeUser(mean-5.*sigma,mean+5.*sigma);
      
      latex = new TLatex(0.55,0.65,Form("#splitline{raw #sigma_{CTR}^{eff} = %.1f ps}{raw #sigma_{CTR}^{gaus} = %.1f ps}",effSigma,fitFunc->GetParameter(2)));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
      latex -> Draw("same");
      
      if( pairsMode.at(pairsIt) == 0 )
      {
        g_tRes_energyCorr_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPoint(g_tRes_energyCorr_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN(),th,effSigma);
        g_tRes_energyCorr_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPointError(g_tRes_energyCorr_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN()-1,0.,5.);
        g_tRes_energyCorr_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPoint(g_tRes_energyCorr_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN(),th,fitFunc->GetParameter(2));
        g_tRes_energyCorr_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPointError(g_tRes_energyCorr_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN()-1,0.,fitFunc->GetParError(2));
        
        g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPoint(g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN(),Vov,effSigma);
        g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPointError(g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN()-1,0.,5.);
        g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPoint(g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN(),Vov,fitFunc->GetParameter(2));
        g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPointError(g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN()-1,0.,fitFunc->GetParError(2));
      }
      if( pairsMode.at(pairsIt) == 1 )
      {
        g_tRes_energyCorr_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPoint(g_tRes_energyCorr_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN(),th,effSigma/sqrt(2));
        g_tRes_energyCorr_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPointError(g_tRes_energyCorr_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN()-1,0.,5.);
        g_tRes_energyCorr_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPoint(g_tRes_energyCorr_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN(),th,fitFunc->GetParameter(2)/sqrt(2));
        g_tRes_energyCorr_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPointError(g_tRes_energyCorr_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN()-1,0.,fitFunc->GetParError(2)/sqrt(2));
        
        g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPoint(g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN(),Vov,effSigma/sqrt(2));
        g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPointError(g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN()-1,0.,5.);
        g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPoint(g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN(),Vov,fitFunc->GetParameter(2)/sqrt(2));
        g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPointError(g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN()-1,0.,fitFunc->GetParError(2)/sqrt(2));
      }
      if( pairsMode.at(pairsIt) == 2 )
      {
        g_tRes_energyCorr_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPoint(g_tRes_energyCorr_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN(),th,effSigma/2.);
        g_tRes_energyCorr_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPointError(g_tRes_energyCorr_effSigma_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN()-1,0.,5.);
        g_tRes_energyCorr_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPoint(g_tRes_energyCorr_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN(),th,fitFunc->GetParameter(2)/2.);
        g_tRes_energyCorr_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())] -> SetPointError(g_tRes_energyCorr_gaus_vs_th[Form("ch%d-ch%d_%s",ch1,ch2,VovLabel.c_str())]->GetN()-1,0.,fitFunc->GetParError(2)/2.);
        
        g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPoint(g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN(),Vov,effSigma/2.);
        g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPointError(g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN()-1,0.,5.);
        g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPoint(g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN(),Vov,fitFunc->GetParameter(2)/2.);
        g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())] -> SetPointError(g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d-ch%d_%s",ch1,ch2,thLabel.c_str())]->GetN()-1,0.,fitFunc->GetParError(2)/2.);
      }
      
      c -> Print(Form("%s/CTR_energyCorr/c_deltaT_energyCorr_%s.png",plotDir.c_str(),label12.c_str()));
      c -> Print(Form("%s/CTR_energyCorr/c_deltaT_energyCorr_%s.pdf",plotDir.c_str(),label12.c_str()));
      
      ++pairsIt;
    }
    
    
    for(auto bar : barsVec)
    {
      unsigned int ch1 = bar.first.first;
      unsigned int ch2 = bar.first.second;
      std::string label1(Form("ch%d_%s",ch1,stepLabel.c_str()));
      std::string label2(Form("ch%d_%s",ch2,stepLabel.c_str()));
      
      unsigned int ch3 = bar.second.first;
      unsigned int ch4 = bar.second.second;
      std::string label3(Form("ch%d_%s",ch3,stepLabel.c_str()));
      std::string label4(Form("ch%d_%s",ch4,stepLabel.c_str()));
      
      std::string label1234(Form("ch%d+ch%d-ch%d+ch%d_%s",ch1,ch2,ch3,ch4,stepLabel.c_str()));
      
      
      //--------------------------------------------------------
      
      
      c = new TCanvas(Form("c_deltaT_energyCorr_%s",label1234.c_str()),Form("c_deltaT_energyCorr_%s",label1234.c_str()));
      // gPad -> SetLogy();
      
      h1_deltaT_energyCorr[label1234] -> SetTitle(Form(";energy-corrected #Deltat [ps];entries"));
      h1_deltaT_energyCorr[label1234] -> SetMarkerColor(kBlue);
      h1_deltaT_energyCorr[label1234] -> SetLineColor(kBlue);
      h1_deltaT_energyCorr[label1234] -> Draw("");
      
      float fitXMin = CTRMeans[label1234] - 1.*CTRSigmas[label1234];
      float fitXMax = CTRMeans[label1234] + 1.*CTRSigmas[label1234];
      TF1* fitFunc = new TF1(Form("fitFunc_energyCorr_%s",label1234.c_str()),"gaus",fitXMin,fitXMax);
      fitFunc -> SetParameters(1,h1_deltaT_energyCorr[label1234]->GetMean(),h1_deltaT_energyCorr[label1234]->GetRMS());
      h1_deltaT_energyCorr[label1234] -> Fit(fitFunc,"QNRSL+");
      fitFunc -> SetLineColor(kBlue-1);
      fitFunc -> SetLineWidth(2);
      fitFunc -> Draw("same");
      
      FindSmallestInterval(vals,h1_deltaT_energyCorr[label1234],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      
      latex = new TLatex(0.55,0.85,Form("#splitline{time walk corr. #sigma_{CTR}^{eff} = %.1f ps}{time walk corr. #sigma_{CTR}^{gaus} = %.1f ps}",effSigma,fitFunc->GetParameter(2)));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlue);
      latex -> Draw("same");
      
      if( g_tRes_effSigma_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())] == NULL )
      {
        g_tRes_effSigma_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())] = new TGraphErrors();
        g_tRes_gaus_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())] = new TGraphErrors();
        g_tRes_energyCorr_effSigma_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())] = new TGraphErrors();
        g_tRes_energyCorr_gaus_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())] = new TGraphErrors();
      }
      if( g_tRes_effSigma_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())] == NULL )
      {          
        g_tRes_effSigma_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())] = new TGraphErrors();
        g_tRes_gaus_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())] = new TGraphErrors();
        g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())] = new TGraphErrors();
        g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())] = new TGraphErrors();
      }
      
      g_tRes_effSigma_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())] -> SetPoint(g_tRes_effSigma_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())]->GetN(),th,effSigma/sqrt(2.));
      g_tRes_effSigma_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())] -> SetPointError(g_tRes_effSigma_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())]->GetN()-1,0.,5.);
      g_tRes_gaus_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())] -> SetPoint(g_tRes_gaus_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())]->GetN(),th,fitFunc->GetParameter(2)/sqrt(2.));
      g_tRes_gaus_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())] -> SetPointError(g_tRes_gaus_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())]->GetN()-1,0.,fitFunc->GetParError(2));
      
      g_tRes_effSigma_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())] -> SetPoint(g_tRes_effSigma_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())]->GetN(),Vov,effSigma/sqrt(2.));
      g_tRes_effSigma_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())] -> SetPointError(g_tRes_effSigma_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())]->GetN()-1,0.,5.);
      g_tRes_gaus_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())] -> SetPoint(g_tRes_gaus_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())]->GetN(),Vov,fitFunc->GetParameter(2)/sqrt(2.));
      g_tRes_gaus_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())] -> SetPointError(g_tRes_gaus_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())]->GetN()-1,0.,fitFunc->GetParError(2));
      
      
      h1_deltaT[label1234] -> SetMarkerColor(kRed);
      h1_deltaT[label1234] -> SetLineColor(kRed);
      h1_deltaT[label1234] -> Draw("same");
      
      fitXMin = CTRMeans[label1234] - 1.*CTRSigmas[label1234];
      fitXMax = CTRMeans[label1234] + 1.*CTRSigmas[label1234];
      fitFunc = new TF1(Form("fitFunc_%s",label1234.c_str()),"gaus",fitXMin,fitXMax);
      fitFunc -> SetParameters(1,h1_deltaT[label1234]->GetMean(),h1_deltaT[label1234]->GetRMS());
      h1_deltaT[label1234] -> Fit(fitFunc,"QNRSL+");
      fitFunc -> SetLineColor(kRed-1);
      fitFunc -> SetLineWidth(2);
      fitFunc -> Draw("same");
      
      FindSmallestInterval(vals,h1_deltaT[label1234],0.68);
      mean = vals[0];
      min = vals[4];
      max = vals[5];
      delta = max-min;
      sigma = 0.5*delta;
      effSigma = sigma;
      
      latex = new TLatex(0.55,0.65,Form("#splitline{raw #sigma_{CTR}^{eff} = %.1f ps}{raw #sigma_{CTR}^{gaus} = %.1f ps}",effSigma,fitFunc->GetParameter(2)));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
      latex -> Draw("same");
      
      g_tRes_energyCorr_effSigma_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())] -> SetPoint(g_tRes_energyCorr_effSigma_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())]->GetN(),th,effSigma/sqrt(2.));
      g_tRes_energyCorr_effSigma_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())] -> SetPointError(g_tRes_energyCorr_effSigma_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())]->GetN()-1,0.,5.);
      g_tRes_energyCorr_gaus_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())] -> SetPoint(g_tRes_energyCorr_gaus_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())]->GetN(),th,fitFunc->GetParameter(2)/sqrt(2.));
      g_tRes_energyCorr_gaus_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())] -> SetPointError(g_tRes_energyCorr_gaus_vs_th[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+VovLabel.c_str())]->GetN()-1,0.,fitFunc->GetParError(2));
      
      g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())] -> SetPoint(g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())]->GetN(),Vov,effSigma/sqrt(2.));
      g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())] -> SetPointError(g_tRes_energyCorr_effSigma_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())]->GetN()-1,0.,5.);
      g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())] -> SetPoint(g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())]->GetN(),Vov,fitFunc->GetParameter(2)/sqrt(2.));
      g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())] -> SetPointError(g_tRes_energyCorr_gaus_vs_Vov[Form("ch%d+ch%d-ch%d+ch%d+_%s",ch1,ch2,ch3,ch4,+thLabel.c_str())]->GetN()-1,0.,fitFunc->GetParError(2));
      
      c -> Print(Form("%s/CTR_energyCorr/c_deltaT_energyCorr_%s.png",plotDir.c_str(),label1234.c_str()));
      c -> Print(Form("%s/CTR_energyCorr/c_deltaT_energyCorr_%s.pdf",plotDir.c_str(),label1234.c_str()));
    }
  }
  
  
  //--------------------------------------------------------
  
  
  int pairsIt = 0;
  for(auto pair : pairsVec)
  {
    unsigned int ch1 = pair.first;
    unsigned int ch2 = pair.second;
    
    c = new TCanvas(Form("c_tRes_vs_th_ch%d-ch%d",ch1,ch2),Form("c_tRes_vs_th_ch%d-ch%d",ch1,ch2));
    // gPad -> SetLogy();
    
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,tResMin,64.,tResMax) );
    if( pairsMode.at(pairsIt) == 0 )
      hPad -> SetTitle(";threshold [DAC];#sigma_{t_{diff}} [ps]");
    if( pairsMode.at(pairsIt) == 1 )
      hPad -> SetTitle(";threshold [DAC];#sigma_{t_{diff}} / #sqrt{2} [ps]");
    if( pairsMode.at(pairsIt) == 2 )
      hPad -> SetTitle(";threshold [DAC];#sigma_{t_{diff}} / 2 [ps]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    int iter = 0;
    for(auto mapIt : VovLabels)
    {
      std::string label(Form("ch%d-ch%d_%s",ch1,ch2,mapIt.first.c_str()));
      TGraph* g_effSigma = g_tRes_effSigma_vs_th[label];
      TGraph* g_gaus = g_tRes_gaus_vs_th[label];
      TGraph* g_energyCorr_effSigma = g_tRes_energyCorr_effSigma_vs_th[label];
      TGraph* g_energyCorr_gaus = g_tRes_energyCorr_gaus_vs_th[label];
      
      g_effSigma -> SetLineColor(1+iter);
      g_effSigma -> SetMarkerColor(1+iter);
      g_effSigma -> SetMarkerStyle(20);
      // g_effSigma -> Draw("PL,same");
      
      g_gaus -> SetLineColor(1+iter);
      g_gaus -> SetMarkerColor(1+iter);
      g_gaus -> SetMarkerStyle(20);
      g_gaus -> Draw("PL,same");
      
      g_energyCorr_effSigma -> SetLineColor(1+iter);
      g_energyCorr_effSigma -> SetMarkerColor(1+iter);
      g_energyCorr_effSigma -> SetMarkerStyle(21);
      // g_energyCorr_effSigma -> Draw("PL,same");
      
      g_energyCorr_gaus -> SetLineColor(1+iter);
      g_energyCorr_gaus -> SetMarkerColor(1+iter);
      g_energyCorr_gaus -> SetMarkerStyle(25);
      g_energyCorr_gaus -> Draw("PL,same");
      
      latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack+iter);
      latex -> Draw("same");
      
      ++iter;
    }
    
    c -> Print(Form("%s/c_tRes_vs_th_ch%d-ch%d.png",plotDir.c_str(),ch1,ch2));
    c -> Print(Form("%s/c_tRes_vs_th_ch%d-ch%d.pdf",plotDir.c_str(),ch1,ch2));
    
    ++pairsIt;
  }
  
  for(auto bar : barsVec)
  {
    unsigned int ch1 = bar.first.first;
    unsigned int ch2 = bar.first.second;
    unsigned int ch3 = bar.second.first;
    unsigned int ch4 = bar.second.second;
    
    c = new TCanvas(Form("c_tRes_vs_th_ch%d+ch%d-ch%d+ch%d",ch1,ch2,ch3,ch4),Form("c_tRes_vs_th_ch%d+ch%d-ch%d+ch%d",ch1,ch2,ch3,ch4));
    // gPad -> SetLogy();
    
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,tResMin,64.,tResMax) );
    hPad -> SetTitle(";threshold [DAC];#sigma_{CTR} / #sqrt{2} [ps]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    int iter = 0;
    for(auto mapIt : VovLabels)
    {
      std::string label(Form("ch%d+ch%d-ch%d+ch%d_%s",ch1,ch2,ch3,ch4,mapIt.first.c_str()));
      TGraph* g_effSigma = g_tRes_effSigma_vs_th[label];
      TGraph* g_gaus = g_tRes_gaus_vs_th[label];
      TGraph* g_energyCorr_effSigma = g_tRes_energyCorr_effSigma_vs_th[label];
      TGraph* g_energyCorr_gaus = g_tRes_energyCorr_gaus_vs_th[label];
      
      g_effSigma -> SetLineColor(1+iter);
      g_effSigma -> SetMarkerColor(1+iter);
      g_effSigma -> SetMarkerStyle(20);
      // g_effSigma -> Draw("PL,same");
      
      g_gaus -> SetLineColor(1+iter);
      g_gaus -> SetMarkerColor(1+iter);
      g_gaus -> SetMarkerStyle(20);
      g_gaus -> Draw("PL,same");
      
      g_energyCorr_effSigma -> SetLineColor(1+iter);
      g_energyCorr_effSigma -> SetMarkerColor(1+iter);
      g_energyCorr_effSigma -> SetMarkerStyle(21);
      // g_energyCorr_effSigma -> Draw("PL,same");
      
      g_energyCorr_gaus -> SetLineColor(1+iter);
      g_energyCorr_gaus -> SetMarkerColor(1+iter);
      g_energyCorr_gaus -> SetMarkerStyle(25);
      g_energyCorr_gaus -> Draw("PL,same");
      
      latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack+iter);
      latex -> Draw("same");
      
      ++iter;
    }
    
    c -> Print(Form("%s/c_tRes_vs_th_ch%d+ch%d_ch%d+ch%d.png",plotDir.c_str(),ch1,ch2,ch3,ch4));
    c -> Print(Form("%s/c_tRes_vs_th_ch%d+ch%d_ch%d+ch%d.pdf",plotDir.c_str(),ch1,ch2,ch3,ch4));
  }
  
  
  //--------------------------------------------------------
  
  
  for(auto pair : pairsVec)
  {
    unsigned int ch1 = pair.first;
    unsigned int ch2 = pair.second;
    
    c = new TCanvas(Form("c_tRes_vs_Vov_ch%d-ch%d",ch1,ch2),Form("c_tRes_vs_Vov_ch%d-ch%d",ch1,ch2));
    // gPad -> SetLogy();
    
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,tResMin,10.,tResMax) );
    hPad -> SetTitle(";V_{ov} [V];#sigma_{t_{diff}} / 2 [ps]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    int iter = 0;
    for(auto mapIt : thLabels)
    {
      std::string label(Form("ch%d-ch%d_%s",ch1,ch2,mapIt.first.c_str()));
      TGraph* g_effSigma = g_tRes_effSigma_vs_Vov[label];
      TGraph* g_gaus = g_tRes_gaus_vs_Vov[label];
      TGraph* g_energyCorr_effSigma = g_tRes_energyCorr_effSigma_vs_Vov[label];
      TGraph* g_energyCorr_gaus = g_tRes_energyCorr_gaus_vs_Vov[label];
      
      g_effSigma -> SetLineColor(1+iter);
      g_effSigma -> SetMarkerColor(1+iter);
      g_effSigma -> SetMarkerStyle(20);
      // g_effSigma -> Draw("PL,same");
      
      g_gaus -> SetLineColor(1+iter);
      g_gaus -> SetMarkerColor(1+iter);
      g_gaus -> SetMarkerStyle(20);
      g_gaus -> Draw("PL,same");
      
      g_energyCorr_effSigma -> SetLineColor(1+iter);
      g_energyCorr_effSigma -> SetMarkerColor(1+iter);
      g_energyCorr_effSigma -> SetMarkerStyle(21);
      // g_energyCorr_effSigma -> Draw("PL,same");
      
      g_energyCorr_gaus -> SetLineColor(1+iter);
      g_energyCorr_gaus -> SetMarkerColor(1+iter);
      g_energyCorr_gaus -> SetMarkerStyle(25);
      g_energyCorr_gaus -> Draw("PL,same");
      
      latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack+iter);
      latex -> Draw("same");
      
      ++iter;
    }
    
    c -> Print(Form("%s/c_tRes_vs_Vov_ch%d-ch%d.png",plotDir.c_str(),ch1,ch2));
    c -> Print(Form("%s/c_tRes_vs_Vov_ch%d-ch%d.pdf",plotDir.c_str(),ch1,ch2));
  }
  
  for(auto bar : barsVec)
  {
    unsigned int ch1 = bar.first.first;
    unsigned int ch2 = bar.first.second;
    unsigned int ch3 = bar.second.first;
    unsigned int ch4 = bar.second.second;
    
    c = new TCanvas(Form("c_tRes_vs_Vov_ch%d+ch%d-ch%d+ch%d",ch1,ch2,ch3,ch4),Form("c_tRes_vs_Vov_ch%d+ch%d-ch%d+ch%d",ch1,ch2,ch3,ch4));
    // gPad -> SetLogy();
    
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,tResMin,10.,tResMax) );
    hPad -> SetTitle(";V_{ov} [V];#sigma_{CTR} / #sqrt{2} [ps]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    int iter = 0;
    for(auto mapIt : thLabels)
    {
      std::string label(Form("ch%d+ch%d-ch%d+ch%d_%s",ch1,ch2,ch3,ch4,mapIt.first.c_str()));
      TGraph* g_effSigma = g_tRes_effSigma_vs_Vov[label];
      TGraph* g_gaus = g_tRes_gaus_vs_Vov[label];
      TGraph* g_energyCorr_effSigma = g_tRes_energyCorr_effSigma_vs_Vov[label];
      TGraph* g_energyCorr_gaus = g_tRes_energyCorr_gaus_vs_Vov[label];
      
      g_effSigma -> SetLineColor(1+iter);
      g_effSigma -> SetMarkerColor(1+iter);
      g_effSigma -> SetMarkerStyle(20);
      // g_effSigma -> Draw("PL,same");
      
      g_gaus -> SetLineColor(1+iter);
      g_gaus -> SetMarkerColor(1+iter);
      g_gaus -> SetMarkerStyle(20);
      g_gaus -> Draw("PL,same");
      
      g_energyCorr_effSigma -> SetLineColor(1+iter);
      g_energyCorr_effSigma -> SetMarkerColor(1+iter);
      g_energyCorr_effSigma -> SetMarkerStyle(21);
      // g_energyCorr_effSigma -> Draw("PL,same");
      
      g_energyCorr_gaus -> SetLineColor(1+iter);
      g_energyCorr_gaus -> SetMarkerColor(1+iter);
      g_energyCorr_gaus -> SetMarkerStyle(25);
      g_energyCorr_gaus -> Draw("PL,same");
      
      latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack+iter);
      latex -> Draw("same");
      
      ++iter;
    }
    
    c -> Print(Form("%s/c_tRes_vs_Vov_ch%d+ch%d_ch%d+ch%d.png",plotDir.c_str(),ch1,ch2,ch3,ch4));
    c -> Print(Form("%s/c_tRes_vs_Vov_ch%d+ch%d_ch%d+ch%d.pdf",plotDir.c_str(),ch1,ch2,ch3,ch4));
  }
  
  
}
