#include <ctime>



void drawStability_energy(const int& firstRun, const int& lastRun)
{
  gStyle->SetLabelSize(0.03, "X");
  gStyle->SetLabelSize(0.04, "Y");
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.07);
  gStyle->SetTitleOffset(0.7, "Y");
  
  unsigned long long int tMin = 1591963200;
  unsigned long long int tMax = 1592481600;
  
  float yMin = 0.90;
  float yMax = 1.50;
  
  float tempMin = 19.;
  float tempMax = 33.;
  
  
  
  //--- draw temp graph
  TGraph* g_temp = new TGraph;
  
  // std::ifstream inFile("temp.txt",std::ios::in);
  // std::string line;
  // while(1)
  // {
  //   getline(inFile,line,'\n');
  //   if( !inFile.good() ) break;

  //   std::stringstream ss(line);
  //   std::string date, hour, temp, dummy;
  //   inFile >> date >> hour >> dummy >> temp >> dummy;
    
  //   struct tm timeinfo;
  //   strptime((date+" "+hour).c_str(),"%F %H:%M:%S",&timeinfo);
  //   time_t timesec =  mktime(&timeinfo);
  //   g_temp -> SetPoint(g_temp->GetN(),timesec,yMin+(yMax-yMin)/(tempMax-tempMin)*(atof(temp.c_str())-tempMin));
  // }


  
  std::map<int,TGraph*> graphs;
  std::map<int,TGraph*> graphs_ped;
  
  std::map<int,TH1F*> histos;
  std::map<int,float> means;
  std::map<int,int> entries;
  
  std::vector<int> channels;
  channels.push_back(49);
  channels.push_back(249);
  for(auto ch : channels)
  {
    graphs[ch] = new TGraph();
    graphs_ped[ch] = new TGraph();
    histos[ch] = new TH1F(Form("histo_ch%d",ch),"",100,0.96,1.04);
  }
  
  
  
  for(int run = firstRun; run <= lastRun; ++run)
  {
    TFile* inFile = TFile::Open(Form("plots/analyzeStability/run%04d_stability.root",run),"READ");
    if( !inFile ) continue;
    
    TTree* tree = (TTree*)( inFile->Get("data") );
    
    unsigned long long int t_time;
    float t_temp;
    unsigned int t_channelID;
    float t_peak1;
    float t_peak1Err;
    float t_sigma1;
    float t_sigma1Err;
    tree -> SetBranchAddress("time",      &t_time);
    tree -> SetBranchAddress("temp",      &t_temp);
    tree -> SetBranchAddress("channelID", &t_channelID);
    tree -> SetBranchAddress("peak1",     &t_peak1);
    tree -> SetBranchAddress("peak1Err",  &t_peak1Err);
    tree -> SetBranchAddress("sigma1",    &t_sigma1);
    tree -> SetBranchAddress("sigma1Err", &t_sigma1Err);
    
    for(int entry = 0; entry < tree->GetEntries(); ++entry)
    {
      tree -> GetEntry(entry);
      graphs[t_channelID] -> SetPoint(graphs[t_channelID]->GetN(),t_time,t_peak1);
    }
    
    inFile -> Close();
    
    
    inFile = TFile::Open(Form("plots/analyzeStability/run%04d_ped_stability.root",run),"READ");
    if( !inFile ) continue;
    
    tree = (TTree*)( inFile->Get("data") );
    tree -> SetBranchAddress("time",      &t_time);
    tree -> SetBranchAddress("temp",      &t_temp);
    tree -> SetBranchAddress("channelID", &t_channelID);
    tree -> SetBranchAddress("peak1",     &t_peak1);
    tree -> SetBranchAddress("peak1Err",  &t_peak1Err);
    tree -> SetBranchAddress("sigma1",    &t_sigma1);
    tree -> SetBranchAddress("sigma1Err", &t_sigma1Err);
    
    for(int entry = 0; entry < tree->GetEntries(); ++entry)
    {
      tree -> GetEntry(entry);
      graphs_ped[t_channelID] -> SetPoint(graphs_ped[t_channelID]->GetN(),t_time,t_peak1);
      
      if( t_time > 1592316000 && t_time < 1592377200 )
      {
        entries[t_channelID] += 1;
        means[t_channelID] += t_peak1;
      }
      
      if( entry == 0 )
        g_temp -> SetPoint(g_temp->GetN(),t_time,yMin+(yMax-yMin)/(tempMax-tempMin)*(t_temp-tempMin));
    }
    
    inFile -> Close();
  }
  
  
  
  //--- normalize graphs to first point
  std::map<int,double> y0;
  for(auto ch : channels)
  {
    double x0;
    double x,y,err;
    
    TGraph* g = graphs[ch];
    g -> GetPoint(0,x0,y0[ch]);
    for(int point = 0; point < g->GetN(); ++point)
    {
      g -> GetPoint(point,x,y);
      g -> SetPoint(point,x,y/y0[ch]);
      g -> GetErrorY(point);
    }
    
    g = graphs_ped[ch];
    g -> GetPoint(0,x0,y0[ch]);
    for(int point = 0; point < g->GetN(); ++point)
    {
      g -> GetPoint(point,x,y);
      g -> SetPoint(point,x,y/y0[ch]);
      g -> GetErrorY(point);
    }
  }
  
  
  
  for(int run = firstRun; run <= lastRun; ++run)
  {
    unsigned long long int t_time;
    float t_temp;
    unsigned int t_channelID;
    float t_peak1;
    float t_peak1Err;
    float t_sigma1;
    float t_sigma1Err;
    
    TFile* inFile = TFile::Open(Form("plots/analyzeStability/run%04d_ped_stability.root",run),"READ");
    if( !inFile ) continue;
    
    TTree* tree = (TTree*)( inFile->Get("data") );
    tree -> SetBranchAddress("time",      &t_time);
    tree -> SetBranchAddress("temp",      &t_temp);
    tree -> SetBranchAddress("channelID", &t_channelID);
    tree -> SetBranchAddress("peak1",     &t_peak1);
    tree -> SetBranchAddress("peak1Err",  &t_peak1Err);
    tree -> SetBranchAddress("sigma1",    &t_sigma1);
    tree -> SetBranchAddress("sigma1Err", &t_sigma1Err);
    
    for(int entry = 0; entry < tree->GetEntries(); ++entry)
    {
      tree -> GetEntry(entry);
      
      if( t_time > 1592316000 && t_time < 1592377200 )
        histos[t_channelID] -> Fill(t_peak1/(means[t_channelID]/entries[t_channelID]));
    }
  }  
  

  
  //--- draw graphs
  TCanvas* c1 = new TCanvas("c","c",1400,600);
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(tMin,yMin,tMax,yMax) );
  hPad -> GetXaxis() -> SetTimeFormat("%b%d %Hh%m%F1970-01-01 00:00:00");
  hPad -> GetXaxis() -> SetTimeDisplay(1);
  hPad -> SetTitle(";time;stability");
  hPad -> Draw();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  std::vector<int> colors;
  colors.push_back(kRed);
  colors.push_back(kGreen);
  int chIt = 0;
  for(auto ch : channels)
  {
    TGraph* g = graphs[ch];
    g -> SetMarkerColor(colors.at(chIt));
    g -> SetMarkerSize(0.5);
    g -> SetLineColor(colors.at(chIt));
    g -> SetLineWidth(1);
    g -> Draw("Lsame");
    
    g = graphs_ped[ch];
    g -> SetMarkerColor(colors.at(chIt)+2);
    g -> SetMarkerSize(0.5);
    g -> SetLineColor(colors.at(chIt)+2);
    g -> SetLineWidth(2);
    g -> Draw("Lsame");
    
    // draw an axis on the right side
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
                              gPad->GetUxmax(),gPad->GetUymax(),tempMin,tempMax,510,"+L");
    axis->SetLineColor(kBlack);
    axis->SetLabelColor(kBlack);
    axis->SetTitle("temperature [#circC]");
    axis->Draw("same");
    
    g_temp -> SetLineWidth(2);
    g_temp -> Draw("L,same");
    
    // TLine* line = new TLine(1592316000,means[ch]/entries[ch]/y0[ch],tMax,means[ch]/entries[ch]/y0[ch]);
    // line -> SetLineStyle(7);
    // line -> Draw("same");
    // line -> SetLineColor(colors.at(chIt)+2);
    
    ++chIt;
  }
  
  TLine* lineMin = new TLine(1592316000,yMin,1592316000,yMax);
  lineMin -> SetLineStyle(7);
  lineMin -> Draw("same");
  TLine* lineMax = new TLine(1592377200,yMin,1592377200,yMax);
  lineMax -> SetLineStyle(7);
  lineMax -> Draw("same");
  
  
  
  TCanvas* c2 = new TCanvas("c2","c2");
  
  int it = 0;
  for(auto ch : channels)
  {
    histos[ch] -> SetLineColor(colors.at(it)+2);
    if( it == 0 ) histos[ch] -> Draw();
    else          histos[ch] -> Draw("same");
    
    ++it;
  }
}
