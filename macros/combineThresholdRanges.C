void CalibrateThresholdDAC(TGraph* g, const float& scale)
{
  for(int point = 0; point < g->GetN(); ++point)
    g -> SetPointX(point,g->GetPointX(point)*scale);
}


void combineThresholdRanges()
{
  float Vov = 3.5;
  
  std::vector<float> DAC_to_mV;
  DAC_to_mV.push_back(8.0);
  DAC_to_mV.push_back(2.0);
  // DAC_to_mV.push_back(1.0);
  // DAC_to_mV.push_back(0.5);
  
  std::vector<int> colors;
  colors.push_back(kRed-1);
  colors.push_back(kOrange-3);
  // colors.push_back(kSpring+5);
  // colors.push_back(kAzure-3);
  
  // std::vector<std::tuple<std::string,std::string,std::string,std::string> > dict;
  std::vector<std::tuple<std::string,std::string> > dict;
  
  dict.push_back(
    std::make_tuple<std::string,std::string>(
      //"/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run726-750_vth2_step3.root",
      //"/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run726-750_vth1_3_step3.root",
      //"/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run726-750_vth1_1_step3.root",
      //"/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run726-750_vth1_0_step3.root"
      "/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run772-799_vth2_step3.root",
      "/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run772-799_vth1_3_step3.root"
      // "/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run772-799_vth1_1_step3.root",
      // "/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run772-799_vth1_0_step3.root"
      )
    );
  
  dict.push_back(
    std::make_tuple<std::string,std::string>(
      "/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run847-874_vth2_step3.root",
      "/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run847-874_vth1_3_step3.root"
      // "/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run847-874_vth1_1_step3.root",
      // "/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run847-874_vth1_0_step3.root"
      )
    );
  
  dict.push_back(
    std::make_tuple<std::string,std::string>(
      "/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run815-842_vth2_step3.root",
      "/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run815-842_vth1_3_step3.root"
      // "/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run815-842_vth1_1_step3.root",
      // "/data/guglielmi/TOFHIR2/moduleCharacterization_single_HPK_HDR2_run815-842_vth1_0_step3.root"
      )
    );
  
  // std::map<float,std::tuple<int,int,int> > dict;
  // dict[Vov] = std::make_tuple<int,int,int>(217,218,219);
  // dict[20.] = std::make_tuple<int,int,int>(225,226,227);
  // dict[30.] = std::make_tuple<int,int,int>(233,234,235);
  // dict[40.] = std::make_tuple<int,int,int>(241,242,243);
  // dict[50.] = std::make_tuple<int,int,int>(249,250,251);
  // dict[60.] = std::make_tuple<int,int,int>(257,258,259);
  // dict[70.] = std::make_tuple<int,int,int>(265,266,267);
  // dict[75.] = std::make_tuple<int,int,int>(273,274,275);
  // dict[80.] = std::make_tuple<int,int,int>(281,282,283);
  // dict[85.] = std::make_tuple<int,int,int>(289,290,291);
  
  for(auto runs : dict)
  {
    // float laserTune = mapIt.first;
    // std::tuple<int,int,int> runs = mapIt.second;
    
    std::vector<TFile*> inFiles;
    inFiles.push_back( TFile::Open(Form("%s",(std::get<0>(runs)).c_str())) );
    inFiles.push_back( TFile::Open(Form("%s",(std::get<1>(runs)).c_str())) );
    // inFiles.push_back( TFile::Open(Form("%s",(std::get<2>(runs)).c_str())) );
    // inFiles.push_back( TFile::Open(Form("%s",(std::get<3>(runs)).c_str())) );
    
    std::vector<TGraphErrors*> g_tRes_energyCorr;
    for(int ii = 0; ii < 2; ++ii)
    {
      g_tRes_energyCorr.push_back( (TGraphErrors*)( inFiles.at(ii)->Get(Form("g_timeRes_vs_th_Vov%.1f_bar00_enBin1",Vov))) );
      CalibrateThresholdDAC(g_tRes_energyCorr.at(ii),DAC_to_mV.at(ii));
    }
    
    TGraphErrors* g_tRes_energyCorr_all = new TGraphErrors();
    for(int ii = 0; ii < 2; ++ii)
    {
      for(int point = 0; point < g_tRes_energyCorr.at(ii)->GetN(); ++point)
        g_tRes_energyCorr_all -> SetPoint(g_tRes_energyCorr_all->GetN(),g_tRes_energyCorr.at(ii)->GetPointX(point),g_tRes_energyCorr.at(ii)->GetPointY(point));
    }
    TF1* fitFunc = new TF1("fitFunc","pol3",0.,200.);
    g_tRes_energyCorr_all -> Fit(fitFunc,"QNRS+");
    fitFunc -> SetLineColor(kRed);
    fitFunc -> SetLineWidth(3);
    
    TCanvas* c = new TCanvas();
    gPad -> SetGridx();
    gPad -> SetGridy();
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,350.,200.) );
    hPad -> SetTitle(";threshold [mV];#sigma_{t}^{diff} / 2 [ps]");
    hPad -> Draw();
    
    for(int ii = 0; ii < 2; ++ii)
    {
      if( ii != 0 ) g_tRes_energyCorr.at(ii) -> SetMarkerSize(0.7);
      g_tRes_energyCorr.at(ii) -> SetMarkerColor(colors.at(ii));
      g_tRes_energyCorr.at(ii) -> SetLineColor(colors.at(ii));
      g_tRes_energyCorr.at(ii) -> Draw("PL,same");
    }
    //fitFunc -> Draw("same");
  }
  
}
