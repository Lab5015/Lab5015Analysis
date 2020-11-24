#include <ctime>



void drawStability_energy(const int& firstRun, const int& lastRun)
{
  gStyle->SetLabelSize(0.03, "X");
  gStyle->SetLabelSize(0.04, "Y");
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.07);
  gStyle->SetTitleOffset(0.7, "Y");
  
  unsigned long long int tMin = 1579262400;
  unsigned long long int tMax = 1579608000;
  
  float yMin = 0.8;
  float yMax = 1.4;
  
  float tempMin = 15.;
  float tempMax = 35.;
  
  
  
  //--- draw temp graph
  TGraph* g_temp = new TGraph;
  
  std::ifstream inFile("temp.txt",std::ios::in);
  std::string line;
  while(1)
  {
    getline(inFile,line,'\n');
    if( !inFile.good() ) break;

    std::stringstream ss(line);
    std::string date, hour, temp, dummy;
    inFile >> date >> hour >> dummy >> temp >> dummy;
    
    struct tm timeinfo;
    strptime((date+" "+hour).c_str(),"%F %H:%M:%S",&timeinfo);
    time_t timesec =  mktime(&timeinfo);
    g_temp -> SetPoint(g_temp->GetN(),timesec,yMin+(yMax-yMin)/(tempMax-tempMin)*(atof(temp.c_str())-tempMin));
  }


  
  std::map<std::pair<int,std::pair<int,int> >,TGraphErrors*> graphs;
  
  for(int arrayIt = 0; arrayIt < 2; ++arrayIt)
    for(int barIt = 0; barIt < 16; ++barIt)
    {
      graphs[std::make_pair(arrayIt,std::make_pair(barIt,0))] = new TGraphErrors();
      graphs[std::make_pair(arrayIt,std::make_pair(barIt,1))] = new TGraphErrors();
    }
  
  
  
  for(int run = firstRun; run <= lastRun; ++run)
  {
    TFile* inFile = TFile::Open(Form("../plots/matricesPlot_run%04d.root",run),"READ");
    if( !inFile ) continue;
    
    TTree* tree = (TTree*)( inFile->Get("data") );
    
    unsigned long long int t_time;
    unsigned int t_arrayID;
    unsigned int t_barID;
    unsigned int t_lrID;
    float t_peak1;
    float t_peak1Err;
    float t_sigma1;
    float t_sigma1Err;
    tree -> SetBranchAddress("time",      &t_time);
    tree -> SetBranchAddress("arrayID",   &t_arrayID);
    tree -> SetBranchAddress("barID",     &t_barID);
    tree -> SetBranchAddress("lrID",      &t_lrID);
    tree -> SetBranchAddress("peak1",     &t_peak1);
    tree -> SetBranchAddress("peak1Err",  &t_peak1Err);
    tree -> SetBranchAddress("sigma1",    &t_sigma1);
    tree -> SetBranchAddress("sigma1Err", &t_sigma1Err);
    
    for(int entry = 0; entry < tree->GetEntries(); ++entry)
    {
      tree -> GetEntry(entry);
      
      TGraphErrors* g = graphs[std::make_pair(t_arrayID,std::make_pair(t_barID,t_lrID))];
      
      g -> SetPoint(g->GetN(),t_time,t_peak1);
      g -> SetPointError(g->GetN()-1,0.,t_peak1Err);
    }
  }
  
  
  
  //--- normalize graphs to first point
  for(int arrayIt = 0; arrayIt < 2; ++arrayIt)
    for(int barIt = 0; barIt < 16; ++barIt)
    {
      double x0,y0;
      double x,y,err;
      
      TGraphErrors* g_L = graphs[std::make_pair(arrayIt,std::make_pair(barIt,0))];
      g_L -> GetPoint(125,x0,y0);
      for(int point = 0; point < g_L->GetN(); ++point)
      {
        g_L -> GetPoint(point,x,y);
        g_L -> SetPoint(point,x,y/y0);
        g_L -> GetErrorY(point);
        g_L -> SetPointError(point,0,err/y0);
      }
      
      TGraphErrors* g_R = graphs[std::make_pair(arrayIt,std::make_pair(barIt,1))];
      g_R -> GetPoint(125,x0,y0);
      for(int point = 0; point < g_R->GetN(); ++point)
      {
        g_R -> GetPoint(point,x,y);
        g_R -> SetPoint(point,x,y/y0);
        g_R -> GetErrorY(point);
        g_R -> SetPointError(point,0,err/y0);
      }
    }
  
  
  
  //--- draw graphs
  for(int arrayIt = 0; arrayIt < 2; ++arrayIt)
  {
    TCanvas* c1 = new TCanvas(Form("array%d",arrayIt),Form("array%d",arrayIt),1200,600);
    
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(tMin,yMin,tMax,yMax) );
    hPad -> GetXaxis() -> SetTimeFormat("%b%d %Hh%m%F1970-01-01 00:00:00");
    hPad -> GetXaxis() -> SetTimeDisplay(1);
    hPad -> SetTitle(";time;stability");
    hPad -> Draw();
    gPad -> SetGridx();
    gPad -> SetGridy();
    
    for(int barIt = 0; barIt < 16; ++barIt)
    {
      if( barIt < 4 ) continue;
      
      TGraphErrors* g_L = graphs[std::make_pair(arrayIt,std::make_pair(barIt,0))];
      g_L -> SetMarkerColor(51+barIt*3);
      g_L -> SetMarkerSize(0.5);
      g_L -> SetLineColor(51+barIt*3);
      g_L -> SetLineStyle(1);
      g_L -> Draw("PLsame");
      
      TGraphErrors* g_R = graphs[std::make_pair(arrayIt,std::make_pair(barIt,1))];
      g_R -> SetMarkerColor(51+barIt*3);
      g_R -> SetMarkerSize(0.5);
      g_R -> SetLineColor(51+barIt*3);
      g_R -> SetLineStyle(2);
      g_R -> Draw("PLsame");
      
      // draw an axis on the right side
      TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
                                gPad->GetUxmax(), gPad->GetUymax(),tempMin,tempMax,510,"+L");
      axis->SetLineColor(kBlack);
      axis->SetLabelColor(kBlack);
      axis->SetTitle("temperature [#circC]");
      axis->Draw("same");
      
      g_temp -> Draw("L,same");
    }
  }
  
  
  
}




void drawStability_tRes(const int& firstRun, const int& lastRun)
{
  gStyle->SetLabelSize(0.03, "X");
  gStyle->SetLabelSize(0.04, "Y");
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.07);
  gStyle->SetTitleOffset(0.7, "Y");
  
  unsigned long long int tMin = 1579262400;
  unsigned long long int tMax = 1579608000;
  
  float yMin = 0.;
  float yMax = 200.;
  
  float tempMin = 15.;
  float tempMax = 35.;  
  
  
  
  //--- draw temp graph
  TGraph* g_temp = new TGraph;
  
  std::ifstream inFile("temp.txt",std::ios::in);
  std::string line;
  while(1)
  {
    getline(inFile,line,'\n');
    if( !inFile.good() ) break;

    std::stringstream ss(line);
    std::string date, hour, temp, dummy;
    inFile >> date >> hour >> dummy >> temp >> dummy;
    
    struct tm timeinfo;
    strptime((date+" "+hour).c_str(),"%F %H:%M:%S",&timeinfo);
    time_t timesec =  mktime(&timeinfo);
    g_temp -> SetPoint(g_temp->GetN(),timesec,yMin+(yMax-yMin)/(tempMax-tempMin)*(atof(temp.c_str())-tempMin));
  }
  
  
  
  std::map<std::pair<int,std::pair<int,int> >,TGraphErrors*> graphs;
  
  for(int arrayIt = 1; arrayIt < 2; ++arrayIt)
    for(int barIt = 0; barIt < 16; ++barIt)
    {
      graphs[std::make_pair(arrayIt,std::make_pair(barIt,0))] = new TGraphErrors();
      graphs[std::make_pair(arrayIt,std::make_pair(barIt,1))] = new TGraphErrors();
    }
  
  
  
  for(int run = firstRun; run <= lastRun; ++run)
  {
    TFile* inFile = TFile::Open(Form("../plots/analyzeTOFPET2_run%04d.root",run),"READ");
    if( !inFile ) continue;
    
    TTree* tree = (TTree*)( inFile->Get("data") );
    
    unsigned long long int t_time;
    unsigned int t_arrayID;
    unsigned int t_barID;
    unsigned int t_lrID;
    float t_tRes;
    float t_tResErr;
    tree -> SetBranchAddress("time",    &t_time);
    tree -> SetBranchAddress("arrayID1",&t_arrayID);
    tree -> SetBranchAddress("barID1",  &t_barID);
    tree -> SetBranchAddress("lrID1",   &t_lrID);
    tree -> SetBranchAddress("tRes",    &t_tRes);
    tree -> SetBranchAddress("tResErr", &t_tResErr);
    
    for(int entry = 0; entry < tree->GetEntries(); ++entry)
    {
      tree -> GetEntry(entry);
      
      TGraphErrors* g = graphs[std::make_pair(t_arrayID,std::make_pair(t_barID,t_lrID))];

      if( t_tRes > 10. && t_tResErr < 20. )
      {
        g -> SetPoint(g->GetN(),t_time,t_tRes/2.);
        g -> SetPointError(g->GetN()-1,0.,t_tResErr/2.);
      }
    }
  }
  
  
  
  //--- draw graphs
  for(int arrayIt = 1; arrayIt < 2; ++arrayIt)
  {
    TCanvas* c1 = new TCanvas(Form("array%d",arrayIt),Form("array%d",arrayIt),1200,600);
    
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(tMin,0.,tMax,200.) );
    hPad -> GetXaxis() -> SetTimeFormat("%b%d %Hh%m%F1970-01-01 00:00:00");
    hPad -> GetXaxis() -> SetTimeDisplay(1);
    hPad -> SetTitle(";time;#sigma_{t_{diff}} / 2 [ps]");
    hPad -> Draw();
    gPad -> SetGridx();
    gPad -> SetGridy();
    
    for(int barIt = 0; barIt < 16; ++barIt)
    {
      if( barIt < 4 ) continue;
      
      TGraphErrors* g_L = graphs[std::make_pair(arrayIt,std::make_pair(barIt,0))];
      g_L -> SetMarkerColor(51+barIt*3);
      g_L -> SetMarkerSize(0.5);
      g_L -> SetLineColor(51+barIt*3);
      g_L -> SetLineStyle(1);
      g_L -> Draw("PLsame");
      
      
      // draw an axis on the right side
      TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
                                gPad->GetUxmax(), gPad->GetUymax(),tempMin,tempMax,510,"+L");
      axis->SetLineColor(kBlack);
      axis->SetLabelColor(kBlack);
      axis->SetTitle("temperature [#circC]");
      axis->Draw("same");
      
      g_temp -> Draw("L,same");
      
      // TGraphErrors* g_R = graphs[std::make_pair(arrayIt,std::make_pair(barIt,1))];
      // g_R -> SetMarkerColor(51+barIt*3);
      // g_R -> SetMarkerSize(0.5);
      // g_R -> SetLineColor(51+barIt*3);
      // g_R -> SetLineStyle(2);
      // g_R -> Draw("PLsame");
    }
  }
}
