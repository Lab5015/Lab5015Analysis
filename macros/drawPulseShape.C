//#include "PulseShape.h"

void drawPulseShape(const int& run, const std::string& vthMode, const float frequency = 10.)
{
  int ch1 = 111;
  int ch2 = 47;
  std::vector<int> channels;
  channels.push_back(ch1);
  channels.push_back(ch2);
  
  cout << ch1 << endl;
  cout << ch2 << endl;
  
  //------------------------
  // define global variables
  // float timeOffset = 237.;
  float timeOffset[2];
  // timeOffset[ch1] = 236.;
  // timeOffset[ch2] = 233.;
  timeOffset[ch1] = 240.;
  timeOffset[ch2] = 237.;
  
  
  //------------------------------
  // open file and define branches
  TFile* inFile = TFile::Open(Form("/data/TOFHIR2/reco/run%04d_ped_e.root",run),"READ");
  //TFile* inFile = TFile::Open(Form("/data/TOFHIR2/reco/run%04d_e.root",run),"READ");
  TTree* data = (TTree*)( inFile->Get("data") );
  
  float step1, step2;
  int channelIdx[128];
  std::vector<float>* tot = 0;
  std::vector<float>* energy = 0;
  std::vector<long long>* time = 0;
  data -> SetBranchStatus("*",0);
  data -> SetBranchStatus("step1",     1); data -> SetBranchAddress("step1",          &step1);
  data -> SetBranchStatus("step2",     1); data -> SetBranchAddress("step2",          &step2);
  data -> SetBranchStatus("channelIdx",1); data -> SetBranchAddress("channelIdx", channelIdx);
  data -> SetBranchStatus("tot",       1); data -> SetBranchAddress("tot",              &tot);
  data -> SetBranchStatus("energy",    1); data -> SetBranchAddress("energy",        &energy);
  data -> SetBranchStatus("time",      1); data -> SetBranchAddress("time",            &time);
  
  
  //------------------
  // define histograms
  //---------------
  // define outfile
  TFile* outFile = new TFile(Form("/data/Lab5015Analysis/pulseShapes/pulseShape_run%04d.root",run),"RECREATE");
  
  std::map<float, std::map<int,TH1F*> > h1_tot_ch1;
  std::map<float, std::map<int,TH1F*> > h1_energy_ch1;
  std::map<float, std::map<int,TH1F*> > h1_tot_ch2;
  std::map<float, std::map<int,TH1F*> > h1_energy_ch2;
  
  std::map<float, std::map<int,TH1F*> > h1_tot_totSel_ch1;
  std::map<float, std::map<int,TH1F*> > h1_energy_totSel_ch1;
  std::map<float, std::map<int,TH1F*> > h1_time1_totSel_ch1;
  std::map<float, std::map<int,TH1F*> > h1_time2_totSel_ch1;
  std::map<float, std::map<int,TH1F*> > h1_tot_totSel_ch2;
  std::map<float, std::map<int,TH1F*> > h1_energy_totSel_ch2;
  std::map<float, std::map<int,TH1F*> > h1_time1_totSel_ch2;
  std::map<float, std::map<int,TH1F*> > h1_time2_totSel_ch2;
  
  std::map<float, std::map<int,TH1F*> > h1_deltaT1_totSel;
  std::map<float, std::map<int,TH1F*> > h1_deltaT2_totSel;
  
  
  //-----------------
  // loop over events  
  int nEntries = data -> GetEntries();
  for(int entry = 0; entry < nEntries; ++entry)
  {
    data -> GetEntry(entry);
    if( entry%10000 == 0 )
    {
      std::cout << ">>>reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    }      
    
    float Vov = step1;
    int vth1 = int(step2/10000.)-1;
    int vth2 = int((step2-10000*(vth1+1))/100.)-1;
    
    int vth = -1;
    if( vthMode.find("vth1") != std::string::npos ) vth = vth1;
    if( vthMode == "vth2" ) vth = vth2;
    
    for(int ch :  channels)
    {
      // if( Vov != 7 ) continue;
      if( channelIdx[ch] < 0 ) continue;
      
      if( ch == ch1 && !h1_tot_ch1[Vov][vth] )
      {
        h1_tot_ch1[Vov][vth]    = new TH1F(Form("h1_tot_ch1_Vov%.1f_vth%02d",Vov,vth),"",15000,-50000.,100000.);
        h1_energy_ch1[Vov][vth] = new TH1F(Form("h1_energy_ch1_Vov%.1f_vth%02d",Vov,vth),"",1000,-0.5,999.5);
        
        h1_tot_totSel_ch1[Vov][vth]     = new TH1F(Form("h1_tot_totSel_ch1_Vov%.1f_vth%02d",Vov,vth),"",1000,0.,100.);
        h1_energy_totSel_ch1[Vov][vth]  = new TH1F(Form("h1_energy_totSel_ch1_Vov%.1f_vth%02d",Vov,vth),"",1000,-0.5,999.5);
        h1_time1_totSel_ch1[Vov][vth]   = new TH1F(Form("h1_time1_totSel_ch1_Vov%.1f_vth%02d",Vov,vth),"",5000,timeOffset[ch1]-50.,timeOffset[ch1]+50.);
        h1_time2_totSel_ch1[Vov][vth]   = new TH1F(Form("h1_time2_totSel_ch1_Vov%.1f_vth%02d",Vov,vth),"",5000,timeOffset[ch1]-50.,timeOffset[ch1]+50.);
        
        h1_deltaT1_totSel[Vov][vth] = new TH1F(Form("h1_deltaT1_totSel_Vov%.1f_vth%02d",Vov,vth),"",10000,-100.,100.);
        h1_deltaT2_totSel[Vov][vth] = new TH1F(Form("h1_deltaT2_totSel_Vov%.1f_vth%02d",Vov,vth),"",10000,-100.,100.);
      }
      if( ch == ch2 && !h1_tot_ch2[Vov][vth] )
      {
        h1_tot_ch2[Vov][vth]    = new TH1F(Form("h1_tot_ch2_Vov%.1f_vth%02d",Vov,vth),"",15000,-50000.,100000.);
        h1_energy_ch2[Vov][vth] = new TH1F(Form("h1_energy_ch2_Vov%.1f_vth%02d",Vov,vth),"",1000,-0.5,999.5);
        
        h1_tot_totSel_ch2[Vov][vth]     = new TH1F(Form("h1_tot_totSel_ch2_Vov%.1f_vth%02d",Vov,vth),"",1000,0.,100.);
        h1_energy_totSel_ch2[Vov][vth]  = new TH1F(Form("h1_energy_totSel_ch2_Vov%.1f_vth%02d",Vov,vth),"",1000,-0.5,999.5);
        h1_time1_totSel_ch2[Vov][vth]   = new TH1F(Form("h1_time1_totSel_ch2_Vov%.1f_vth%02d",Vov,vth),"",1000,timeOffset[ch1]-50.,timeOffset[ch1]+50.);
        h1_time2_totSel_ch2[Vov][vth]   = new TH1F(Form("h1_time2_totSel_ch2_Vov%.1f_vth%02d",Vov,vth),"",1000,timeOffset[ch1]-50.,timeOffset[ch1]+50.);
      }
      
      if( ch == ch1 )
      {
        h1_tot_ch1[Vov][vth] -> Fill( (*tot)[channelIdx[ch1]]/1000. );
        h1_energy_ch1[Vov][vth] -> Fill( (*energy)[channelIdx[ch1]] );
      }
      if( ch == ch2 )
      {
        h1_tot_ch2[Vov][vth] -> Fill( (*tot)[channelIdx[ch2]]/1000. );
        h1_energy_ch2[Vov][vth] -> Fill( (*energy)[channelIdx[ch2]] );
      }
      
      if( (*tot)[channelIdx[ch]]/1000. < 0. || (*tot)[channelIdx[ch]]/1000. > 100. ) continue;
      
      long int scale = 1000000000/(frequency);
      
      if( ch == ch1 )
      {
        h1_tot_totSel_ch1[Vov][vth] -> Fill( (*tot)[channelIdx[ch1]]/1000. );
        h1_energy_totSel_ch1[Vov][vth] -> Fill( (*energy)[channelIdx[ch1]] );
        h1_time1_totSel_ch1[Vov][vth] -> Fill( ((*time)[channelIdx[ch1]]%scale)/1000. );
        h1_time2_totSel_ch1[Vov][vth] -> Fill( ((*time)[channelIdx[ch1]]%scale)/1000. + (*tot)[channelIdx[ch1]]/1000. );
      }
      if( ch == ch2 )
      {
        h1_tot_totSel_ch2[Vov][vth] -> Fill( (*tot)[channelIdx[ch2]]/1000. );
        h1_energy_totSel_ch2[Vov][vth] -> Fill( (*energy)[channelIdx[ch2]] );
        h1_time1_totSel_ch2[Vov][vth] -> Fill( ((*time)[channelIdx[ch2]]%scale)/1000. );
        h1_time2_totSel_ch2[Vov][vth] -> Fill( ((*time)[channelIdx[ch2]]%scale)/1000. + (*tot)[channelIdx[ch2]]/1000. );
      }
      
      //std::cout << time[ch1] << " - " << time[ch2] << " - " << scale << " - " << time[ch1] << " - " << (time[ch1]%scale)/1000. << " - " << time[ch2] << " - " << (time[ch2]%scale)/1000. << std::endl;
      
      // h1_deltaT1_totSel[Vov][vth] -> Fill( ((*time)[channelIdx[ch1]]-(*time)[channelIdx[ch2]])/1000. );
      // h1_deltaT2_totSel[Vov][vth] -> Fill( (float((*time)[channelIdx[ch1]]-(*time)[channelIdx[ch2]])+(*tot)[channelIdx[ch1]])/1000. );
    }
  }
  std::cout << std::endl;
  
  
  
  //-----------------
  // draw pulse shape
  float dac_to_mV = -1.;
  if( vthMode == "vth2"   ) dac_to_mV = 8.;
  if( vthMode == "vth1_4" ) dac_to_mV = 4.;
  if( vthMode == "vth1_3" ) dac_to_mV = 2.;
  if( vthMode == "vth1_1" ) dac_to_mV = 1.;
  if( vthMode == "vth1_0" ) dac_to_mV = 0.5;
  
  
  std::map<float, TGraphErrors*> g_N_ch1;
  std::map<float, TGraphErrors*> g_tot_ch1;
  std::map<float, TGraphErrors*> g_energy_ch1;
  std::map<float, TGraphErrors*> g_N_ch2;
  std::map<float, TGraphErrors*> g_tot_ch2;
  std::map<float, TGraphErrors*> g_energy_ch2;
  
  std::map<float, TGraphErrors*> g_N_totSel_ch1;
  std::map<float, TGraphErrors*> g_tot_totSel_ch1;
  std::map<float, TGraphErrors*> g_energy_totSel_ch1;
  std::map<float, TGraphErrors*> g_ps_totSel_ch1;
  std::map<float, TGraphErrors*> g_N_totSel_ch2;
  std::map<float, TGraphErrors*> g_tot_totSel_ch2;
  std::map<float, TGraphErrors*> g_energy_totSel_ch2;
  std::map<float, TGraphErrors*> g_ps_totSel_ch2;
  
  std::map<float, TGraphErrors*> g_ps_totSel_deltaT;
  
  for(auto mapIt : h1_tot_ch1)
    {
      float Vov = mapIt.first;
      for(auto mapIt2 : mapIt.second)
	{
	  int vth = mapIt2.first;
	  TH1F* histo = mapIt2.second;
	  
	  if( !g_N_ch1[Vov] ) g_N_ch1[Vov] = new TGraphErrors();
	  g_N_ch1[Vov] -> SetPoint(g_N_ch1[Vov]->GetN(),vth,histo->Integral());
          g_N_ch1[Vov] -> SetPointError(g_N_ch1[Vov]->GetN()-1,0,sqrt(histo->Integral()));
	  
	  if( !g_tot_ch1[Vov] ) g_tot_ch1[Vov] = new TGraphErrors();
	  g_tot_ch1[Vov] -> SetPoint(g_tot_ch1[Vov]->GetN(),vth,histo->GetMean());
	  g_tot_ch1[Vov] -> SetPointError(g_tot_ch1[Vov]->GetN()-1,0.,histo->GetRMS());
	  
	  histo -> Write();
	  
	  histo = h1_tot_totSel_ch1[Vov][vth];
	  if( histo->Integral() <= 0. ) continue;
	  
	  if( !g_N_totSel_ch1[Vov] ) g_N_totSel_ch1[Vov] = new TGraphErrors();
	  g_N_totSel_ch1[Vov] -> SetPoint(g_N_totSel_ch1[Vov]->GetN(),vth,histo->Integral());
	  
	  if( !g_tot_totSel_ch1[Vov] ) g_tot_totSel_ch1[Vov] = new TGraphErrors();
	  g_tot_totSel_ch1[Vov] -> SetPoint(g_tot_totSel_ch1[Vov]->GetN(),vth,histo->GetMean());
	  g_tot_totSel_ch1[Vov] -> SetPointError(g_tot_totSel_ch1[Vov]->GetN()-1,0.,histo->GetRMS());
	  
	  histo -> Write();
	}
    }
  for(auto mapIt : h1_tot_ch2)
    {
      float Vov = mapIt.first;
      for(auto mapIt2 : mapIt.second)
	{
	  int vth = mapIt2.first;
	  TH1F* histo = mapIt2.second;
	  
	  if( !g_N_ch2[Vov] ) g_N_ch2[Vov] = new TGraphErrors();
	  g_N_ch2[Vov] -> SetPoint(g_N_ch2[Vov]->GetN(),vth,histo->Integral());
	  
	  if( !g_tot_ch2[Vov] ) g_tot_ch2[Vov] = new TGraphErrors();
	  g_tot_ch2[Vov] -> SetPoint(g_tot_ch2[Vov]->GetN(),vth,histo->GetMean());
	  g_tot_ch2[Vov] -> SetPointError(g_tot_ch2[Vov]->GetN()-1,0.,histo->GetRMS());
	  
	  histo -> Write();
	  
	  histo = h1_tot_totSel_ch2[Vov][vth];
	  if( histo->Integral() <= 0. ) continue;
	  
	  if( !g_N_totSel_ch2[Vov] ) g_N_totSel_ch2[Vov] = new TGraphErrors();
	  g_N_totSel_ch2[Vov] -> SetPoint(g_N_totSel_ch2[Vov]->GetN(),vth,histo->Integral());
	  
	  if( !g_tot_totSel_ch2[Vov] ) g_tot_totSel_ch2[Vov] = new TGraphErrors();
	  g_tot_totSel_ch2[Vov] -> SetPoint(g_tot_totSel_ch2[Vov]->GetN(),vth,histo->GetMean());
	  g_tot_totSel_ch2[Vov] -> SetPointError(g_tot_totSel_ch2[Vov]->GetN()-1,0.,histo->GetRMS());
	  
	  histo -> Write();
	}
    }
  
  for(auto mapIt : h1_energy_ch1)
    {
      float Vov = mapIt.first;
      for(auto mapIt2 : mapIt.second)
	{
	  int vth = mapIt2.first;
	  TH1F* histo = mapIt2.second;
	  if( histo->Integral() <= 0. ) continue;
	  
	  if( !g_energy_ch1[Vov] ) g_energy_ch1[Vov] = new TGraphErrors();
	  g_energy_ch1[Vov] -> SetPoint(g_energy_ch1[Vov]->GetN(),vth,histo->GetMean());
	  g_energy_ch1[Vov] -> SetPointError(g_energy_ch1[Vov]->GetN()-1,0.,histo->GetRMS());
	  
	  histo -> Write();
	  
	  histo = h1_energy_totSel_ch1[Vov][vth];
	  if( histo->Integral() <= 0. ) continue;
      
	  if( !g_energy_totSel_ch1[Vov] ) g_energy_totSel_ch1[Vov] = new TGraphErrors();
	  g_energy_totSel_ch1[Vov] -> SetPoint(g_energy_totSel_ch1[Vov]->GetN(),vth,histo->GetMean());
	  g_energy_totSel_ch1[Vov] -> SetPointError(g_energy_totSel_ch1[Vov]->GetN()-1,0.,histo->GetRMS());
	  
	  histo -> Write();
	}
    }
  for(auto mapIt : h1_energy_ch2)
    {
      float Vov = mapIt.first;
      for(auto mapIt2 : mapIt.second)
	{
	  int vth = mapIt2.first;
	  TH1F* histo = mapIt2.second;
	  if( histo->Integral() <= 0. ) continue;
      
	  if( !g_energy_ch2[Vov] ) g_energy_ch2[Vov] = new TGraphErrors();
	  g_energy_ch2[Vov] -> SetPoint(g_energy_ch2[Vov]->GetN(),vth,histo->GetMean());
	  g_energy_ch2[Vov] -> SetPointError(g_energy_ch2[Vov]->GetN()-1,0.,histo->GetRMS());
	  
	  histo -> Write();
	  
	  histo = h1_energy_totSel_ch2[Vov][vth];
	  if( histo->Integral() <= 0. ) continue;
      
	  if( !g_energy_totSel_ch2[Vov] ) g_energy_totSel_ch2[Vov] = new TGraphErrors();
	  g_energy_totSel_ch2[Vov] -> SetPoint(g_energy_totSel_ch2[Vov]->GetN(),vth,histo->GetMean());
	  g_energy_totSel_ch2[Vov] -> SetPointError(g_energy_totSel_ch2[Vov]->GetN()-1,0.,histo->GetRMS());
	}
    }
  
  /*
  for(auto mapIt : h1_deltaT1_totSel)
    {
      float Vov = mapIt.first;
      for(auto mapIt2 : mapIt.second)
	{
	  int vth = mapIt2.first;
	  TH1F* histo = mapIt2.second;
	  if( histo->Integral() <= 0. ) continue;
	  
	  if( !g_ps_totSel_deltaT[Vov] ) g_ps_totSel_deltaT[Vov] = new TGraphErrors();
	  g_ps_totSel_deltaT[Vov] -> SetPoint(g_ps_totSel_deltaT[Vov]->GetN(),histo->GetMean(),vth*dac_to_mV);
	  g_ps_totSel_deltaT[Vov] -> SetPointError(g_ps_totSel_deltaT[Vov]->GetN()-1,histo->GetMeanError(),0.);
	  
	  histo -> Write();
	  
	  histo = h1_deltaT2_totSel[Vov][vth];
	  if( histo->Integral() <= 0. ) continue;
	  
	  g_ps_totSel_deltaT[Vov] -> SetPoint(g_ps_totSel_deltaT[Vov]->GetN(),histo->GetMean(),vth*dac_to_mV);
	  g_ps_totSel_deltaT[Vov] -> SetPointError(g_ps_totSel_deltaT[Vov]->GetN()-1,histo->GetMeanError(),0.);
	  
	  histo -> Write();
	}
    }
  */
  
  for(auto mapIt : h1_time1_totSel_ch1)
  {
    float Vov = mapIt.first;
    for(auto mapIt2 : mapIt.second)
    {
      int vth = mapIt2.first;
      TH1F* histo = mapIt2.second;
      if( histo->Integral() <= 1000. ) continue;
      
      if( !g_ps_totSel_ch1[Vov] ) g_ps_totSel_ch1[Vov] = new TGraphErrors();
      g_ps_totSel_ch1[Vov] -> SetPoint(g_ps_totSel_ch1[Vov]->GetN(),histo->GetMean()-timeOffset[ch1],vth*dac_to_mV);
      g_ps_totSel_ch1[Vov] -> SetPointError(g_ps_totSel_ch1[Vov]->GetN()-1,histo->GetMeanError(),0.);
      
      histo -> Write();
    }
    for(auto mapIt2 : mapIt.second)
    {
      int vth = mapIt2.first;
      TH1F* histo = h1_time2_totSel_ch1[Vov][vth];
      if( histo->Integral() <= 1000. ) continue;
      
      g_ps_totSel_ch1[Vov] -> SetPoint(g_ps_totSel_ch1[Vov]->GetN(),histo->GetMean()-timeOffset[ch1],vth*dac_to_mV);
      g_ps_totSel_ch1[Vov] -> SetPointError(g_ps_totSel_ch1[Vov]->GetN()-1,histo->GetMeanError(),0.);
      
      histo -> Write();
    }
  }
  for(auto mapIt : h1_time1_totSel_ch2)
  {
    float Vov = mapIt.first;
    for(auto mapIt2 : mapIt.second)
    {
      int vth = mapIt2.first;
      TH1F* histo = mapIt2.second;
      if( histo->Integral() <= 1000. ) continue;
      
      if( !g_ps_totSel_ch2[Vov] ) g_ps_totSel_ch2[Vov] = new TGraphErrors();
      g_ps_totSel_ch2[Vov] -> SetPoint(g_ps_totSel_ch2[Vov]->GetN(),histo->GetMean()-timeOffset[ch2],vth*dac_to_mV);
      g_ps_totSel_ch2[Vov] -> SetPointError(g_ps_totSel_ch2[Vov]->GetN()-1,histo->GetMeanError(),0.);
      
      histo -> Write();
    }
    for(auto mapIt2 : mapIt.second)
    {
      int vth = mapIt2.first;
      TH1F* histo = h1_time2_totSel_ch2[Vov][vth];
      if( histo->Integral() <= 1000. ) continue;
      
      g_ps_totSel_ch2[Vov] -> SetPoint(g_ps_totSel_ch2[Vov]->GetN(),histo->GetMean()-timeOffset[ch2],vth*dac_to_mV);
      g_ps_totSel_ch2[Vov] -> SetPointError(g_ps_totSel_ch2[Vov]->GetN()-1,histo->GetMeanError(),0.);
      
      histo -> Write();
    }
  }
  
  
  //-----------
  // draw plots
  std::string plotDir(Form("/var/www/html/TOFHIR2A/debug/run%04d/",run));
  system(Form("mkdir -p %s",plotDir.c_str()));
  
  TCanvas* c;
  TH1F* hPad;
  
  
  for(auto mapIt : h1_tot_ch1)
    {  
      float Vov = mapIt.first;
      
      c = new TCanvas("c","c");
      hPad = (TH1F*)( gPad->DrawFrame(-0.5,0.,63.5,12000*frequency/10.) );
      hPad -> SetTitle(Form(";%s [DAC]; number of hits in ch111",vthMode.c_str()));
      hPad -> Draw();
      g_N_ch1[Vov] -> SetMarkerColor(kRed);
      g_N_ch1[Vov] -> SetMarkerStyle(22);
      g_N_ch1[Vov] -> Draw("P,same");
      if( g_N_ch2[Vov] ) g_N_ch2[Vov] -> SetMarkerColor(kBlue);
      if( g_N_ch2[Vov] ) g_N_ch2[Vov] -> SetMarkerStyle(23);
      if( g_N_ch2[Vov] ) g_N_ch2[Vov] -> Draw("PL,same");
      g_N_totSel_ch1[Vov] -> SetMarkerColor(kRed-4);
      g_N_totSel_ch1[Vov] -> SetMarkerStyle(26);
      g_N_totSel_ch1[Vov] -> Draw("P,same");
      if( g_N_totSel_ch2[Vov] ) g_N_totSel_ch2[Vov] -> SetMarkerColor(kBlue-4);
      if( g_N_totSel_ch2[Vov] ) g_N_totSel_ch2[Vov] -> SetMarkerStyle(32);
      if( g_N_totSel_ch2[Vov] ) g_N_totSel_ch2[Vov] -> Draw("PL,same");
      
      TF1* f_sigmoid_ch1 = new TF1("f_sigmoid_ch1","[0]*(1-0.5*(1.+TMath::Erf((x-[1])/[2])))",0.,64.);
      f_sigmoid_ch1 -> SetNpx(10000);
      f_sigmoid_ch1 -> SetLineWidth(1);
      f_sigmoid_ch1 -> SetLineColor(kRed-4);
      int index = 0;
      for (int j = 0; j < g_N_ch1[Vov] -> GetN(); j++ ){
        if ( g_N_ch1[Vov]->GetPointY(j) < g_N_ch1[Vov]->GetPointY(0)/2){
          index = j;
          break;
        } 
      }
      //f_sigmoid_ch1 -> SetParameters(g_N_ch1[Vov]->GetPointY(0),12.,3.);
      f_sigmoid_ch1 -> SetParameters(g_N_ch1[Vov]->GetPointY(0), g_N_ch1[Vov]->GetPointX(index),3.);
      g_N_ch1[Vov] -> Fit(f_sigmoid_ch1,"QRS");
      f_sigmoid_ch1->Draw("same");
      
      TLatex* latex_ch1 = new TLatex(0.40,0.80,Form("amplitude = %.1f mV",dac_to_mV*f_sigmoid_ch1->GetParameter(1)));
      latex_ch1 -> SetNDC();
      latex_ch1 -> SetTextFont(82);
      latex_ch1 -> SetTextSize(0.04);
      latex_ch1 -> SetTextAlign(11);
      latex_ch1 -> SetTextColor(kRed-4);
      latex_ch1 -> Draw("same");
      
      c -> Print(Form("%s/g_N_Vov%.1f.png",plotDir.c_str(),Vov));
      
      delete c;
    }
  
  for(auto mapIt : h1_tot_ch1)
  {  
    float Vov = mapIt.first;
    
    c = new TCanvas("c","c");
    hPad = (TH1F*)( gPad->DrawFrame(-0.5,0.,63.5,20.) );
    hPad -> SetTitle(Form(";%s [DAC]; ToT [ns]",vthMode.c_str()));
    hPad -> Draw();
    g_tot_ch1[Vov] -> SetMarkerColor(kRed);
    g_tot_ch1[Vov] -> SetMarkerStyle(22);
    g_tot_ch1[Vov] -> Draw("PL,same");
    if( g_tot_ch2[Vov] ) g_tot_ch2[Vov] -> SetMarkerColor(kBlue);
    if( g_tot_ch2[Vov] ) g_tot_ch2[Vov] -> SetMarkerStyle(23);
    if( g_tot_ch2[Vov] ) g_tot_ch2[Vov] -> Draw("PL,same");
    g_tot_totSel_ch1[Vov] -> SetMarkerColor(kRed-4);
    g_tot_totSel_ch1[Vov] -> SetMarkerStyle(26);
    g_tot_totSel_ch1[Vov] -> Draw("PL,same");
    if( g_tot_totSel_ch2[Vov] ) g_tot_totSel_ch2[Vov] -> SetMarkerColor(kBlue-4);
    if( g_tot_totSel_ch2[Vov] ) g_tot_totSel_ch2[Vov] -> SetMarkerStyle(32);
    if( g_tot_totSel_ch2[Vov] ) g_tot_totSel_ch2[Vov] -> Draw("PL,same");
    c -> Print(Form("%s/g_tot_Vov%.1f.png",plotDir.c_str(),Vov));
    delete c;
  }
  
  for(auto mapIt : h1_energy_ch1)
  {  
    float Vov = mapIt.first;
    
    c = new TCanvas("c","c");
    hPad = (TH1F*)( gPad->DrawFrame(-0.5,0.,63.5,1023.5) );
    hPad -> SetTitle(Form(";%s [DAC]; energy [ADC]",vthMode.c_str()));
    hPad -> Draw();
    g_energy_ch1[Vov] -> SetMarkerColor(kRed);
    g_energy_ch1[Vov] -> SetMarkerStyle(22);
    g_energy_ch1[Vov] -> Draw("PL,same");
    if( g_energy_ch2[Vov] ) g_energy_ch2[Vov] -> SetMarkerColor(kBlue);
    if( g_energy_ch2[Vov] )g_energy_ch2[Vov] -> SetMarkerStyle(23);
    if( g_energy_ch2[Vov] )g_energy_ch2[Vov] -> Draw("PL,same");
    g_energy_totSel_ch1[Vov] -> SetMarkerColor(kRed-4);
    g_energy_totSel_ch1[Vov] -> SetMarkerStyle(26);
    g_energy_totSel_ch1[Vov] -> Draw("PL,same");
    if( g_energy_totSel_ch2[Vov] ) g_energy_totSel_ch2[Vov] -> SetMarkerColor(kBlue-4);
    if( g_energy_totSel_ch2[Vov] ) g_energy_totSel_ch2[Vov] -> SetMarkerStyle(32);
    if( g_energy_totSel_ch2[Vov] ) g_energy_totSel_ch2[Vov] -> Draw("PL,same");
    c -> Print(Form("%s/g_energy_Vov%.1f.png",plotDir.c_str(),Vov));
    delete c;
  }
  
  
  /*
  for(auto mapIt : h1_tot_ch1)
    {  
      float Vov = mapIt.first;
      
      c = new TCanvas("c","c");
      //hPad = (TH1F*)( gPad->DrawFrame(-2.,0.,20.,40.) );
      //hPad = (TH1F*)( gPad->DrawFrame(-2.,0.,20.,100.) );
      hPad = (TH1F*)( gPad->DrawFrame(-2.,0.,20.,200.) );
      hPad -> SetTitle(Form(";time [ns]; pulse shape [mV]"));
      hPad -> Draw();
      g_ps_totSel_deltaT[Vov] -> SetMarkerColor(kGray+2);
      g_ps_totSel_deltaT[Vov] -> SetLineColor(kGray+2);
      g_ps_totSel_deltaT[Vov] -> SetMarkerStyle(22);
      g_ps_totSel_deltaT[Vov] -> Draw("P,same");
      c -> Print(Form("%s/g_ps_deltaT_Vov%.1f.png",plotDir.c_str(),Vov));
      delete c;
    }
  */
  
  for(auto mapIt : h1_time1_totSel_ch1)
  {  
    float Vov = mapIt.first;
    
    float fitXMin = 0.;
    float fitXMax = 999.;
    
    c = new TCanvas("c","c");
    //hPad = (TH1F*)( gPad->DrawFrame(-2.,0.,20.,40.) );
    //hPad = (TH1F*)( gPad->DrawFrame(-2.,0.,20.,100.) );
    hPad = (TH1F*)( gPad->DrawFrame(-2.,0.,20.,65*dac_to_mV) );
    hPad -> SetTitle(Form(";time [ns]; pulse shape [mV]"));
    hPad -> Draw();
    g_ps_totSel_ch1[Vov] -> SetLineColor(kRed-4);
    g_ps_totSel_ch1[Vov] -> SetMarkerColor(kRed-4);
    g_ps_totSel_ch1[Vov] -> SetMarkerStyle(26);
    g_ps_totSel_ch1[Vov] -> Draw("P,same");
    if( g_ps_totSel_ch2[Vov] ) g_ps_totSel_ch2[Vov] -> SetMarkerColor(kBlue-4);
    if( g_ps_totSel_ch2[Vov] ) g_ps_totSel_ch2[Vov] -> SetLineColor(kBlue-4);
    if( g_ps_totSel_ch2[Vov] ) g_ps_totSel_ch2[Vov] -> SetMarkerStyle(32);
    if( g_ps_totSel_ch2[Vov] ) g_ps_totSel_ch2[Vov] -> Draw("P,same");
    
    float slewRate = 0.;
    TF1* fitFunc_ch1 = new TF1("fitFunc_ch1","pol1",0.,10.);
    for(int point1 = 0; point1 < g_ps_totSel_ch1[Vov]->GetN()-4; ++point1)
    {
      TGraph* g_temp = new TGraph();
      for(int point2 = point1; point2 < point1+4; ++point2)
      {
        g_temp -> SetPoint(g_temp->GetN(),g_ps_totSel_ch1[Vov]->GetPointX(point2),g_ps_totSel_ch1[Vov]->GetPointY(point2));
      }
      
      TF1* f_temp = new TF1("f_temp","pol1",0.,100.);
      g_temp -> Fit(f_temp,"QNRS");
      
      if( f_temp->GetParameter(1) > slewRate )
      {
        slewRate = f_temp->GetParameter(1);
        fitFunc_ch1 -> SetParameters(f_temp->GetParameter(0),f_temp->GetParameter(1));
      }
    }
    // for(int point = 0; point < g_ps_totSel_ch1[Vov]->GetN(); ++point)
    //   if( g_ps_totSel_ch1[Vov]->GetPointY(point) > 30. )
    //   {
    //     fitXMin = g_ps_totSel_ch1[Vov]->GetPointX(point);
    //     break;
    //   }
    // for(int point = 0; point < g_ps_totSel_ch1[Vov]->GetN(); ++point)
    //   if( g_ps_totSel_ch1[Vov]->GetPointY(point) > 100. )
    //   {
    //     fitXMax = g_ps_totSel_ch1[Vov]->GetPointX(point);
    //     std::cout << fitXMax << std::endl;
    //     break;
    //   }
    // TF1* fitFunc_ch1 = new TF1("fitFunc_ch1","pol1",0.,7.);
    // fitFunc_ch1 -> SetParameters(0.,250.);
    // g_ps_totSel_ch1[Vov] -> Fit(fitFunc_ch1,"QNS+","",fitXMin,fitXMax);
    fitFunc_ch1 -> SetLineColor(kRed-4);
    // fitFunc_ch1 -> Draw("same");
    TLatex* latex_ch1 = new TLatex(0.40,0.80,Form("slew rate = %.1f mV/ns",fitFunc_ch1->GetParameter(1)));
    latex_ch1 -> SetNDC();
    latex_ch1 -> SetTextFont(82);
    latex_ch1 -> SetTextSize(0.04);
    latex_ch1 -> SetTextAlign(11);
    latex_ch1 -> SetTextColor(kRed-4);
    // latex_ch1 -> Draw("same");
    
    if( g_ps_totSel_ch2[Vov] )
    {
      fitXMin = 0.;
      fitXMax = 999.;
      for(int point = 0; point < g_ps_totSel_ch2[Vov]->GetN(); ++point)
        if( g_ps_totSel_ch2[Vov]->GetPointY(point) > 30. )
        {
          fitXMin = g_ps_totSel_ch2[Vov]->GetPointX(point);
          break;
        }
      for(int point = 0; point < g_ps_totSel_ch2[Vov]->GetN(); ++point)
        if( g_ps_totSel_ch2[Vov]->GetPointY(point) > 100. )
        {
          fitXMax = g_ps_totSel_ch2[Vov]->GetPointX(point);
          break;
        }
      TF1* fitFunc_ch2 = new TF1("fitFunc_ch2","pol1",0.,7.);
      fitFunc_ch2 -> SetParameters(0.,250.);
      g_ps_totSel_ch2[Vov] -> Fit(fitFunc_ch2,"QNS+","",fitXMin,fitXMax);
      fitFunc_ch2 -> SetLineColor(kBlue-4);
      //fitFunc_ch2 -> Draw("same");
      TLatex* latex_ch2 = new TLatex(0.40,0.76,Form("slew rate = %.1f mV/ns",fitFunc_ch2->GetParameter(1)));
      latex_ch2 -> SetNDC();
      latex_ch2 -> SetTextFont(82);
      latex_ch2 -> SetTextSize(0.04);
      latex_ch2 -> SetTextAlign(11);
      latex_ch2 -> SetTextColor(kBlue-4);
      //latex_ch2 -> Draw("same");
    }
    
    c -> Print(Form("%s/g_ps_ch1_ch2_Vov%.1f.png",plotDir.c_str(),Vov));
    delete c;
  }
  
  //-----------
  // save plots
  
  for(auto mapIt : g_ps_totSel_ch1)
  {
    mapIt.second -> Write(Form("g_ps_totSel_ch1_Vov%.1f",mapIt.first));
  }
  for(auto mapIt : g_ps_totSel_ch2)
  {
    mapIt.second -> Write(Form("g_ps_totSel_ch2_Vov%.1f",mapIt.first));
  }

  for (auto mapIt:  g_N_ch1)
  {
    mapIt.second -> Write(Form("g_N_ch1_Vov%.1f",mapIt.first));
  } 

  gApplication->Terminate(); 
  
  /*
  p_tot_vs_qT -> Write();
  
  g_N -> Write("g_N");
  g_tot -> Write("g_tot");
  g_energy -> Write("g_energy");
  
  g_N_totSel -> Write("g_N_totSel");
  g_tot_totSel -> Write("g_tot_totSel");
  g_energy_totSel -> Write("g_energy_totSel");
  g_ps_totSel -> Write("g_ps_totSel");
  */

  outFile -> Close();
}
