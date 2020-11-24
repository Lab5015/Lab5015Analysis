void studyLinearity()
{
  TFile* sipmParams = TFile::Open("sipm_spec_input_HDR2-015-v2-1e13.root","READ");
  TF1* f_PDE = new TF1("func","[0]*(1-exp(-1.*[1]*x))",0.,10.);
  f_PDE -> SetParameters(39.4321,0.738063);
  TF1* f_gain = (TF1*)( sipmParams->Get("fGain_vs_OV") );
  TF1* f_ENF = (TF1*)( sipmParams->Get("fENF_vs_OV") );
  
  
  TFile* inFile = TFile::Open("histograms.root","READ");
  
  std::map<int,TH1F*> histos;
  histos[3] = (TH1F*)( inFile->Get("h1_energy_ov3_ch49") );
  histos[5] = (TH1F*)( inFile->Get("h1_energy_ov5_ch49") );
  histos[7] = (TH1F*)( inFile->Get("h1_energy_ov7_ch49") );
                       
  std::map<int,TF1*> fitFuncs_511;
  fitFuncs_511[3] = new TF1("fitFunc_511_ov3","gaus(0)",15.,23.);
  fitFuncs_511[3] -> SetParameters(100000,20.,2.);
  fitFuncs_511[5] = new TF1("fitFunc_511_ov5","gaus(0)",25.,35.);
  fitFuncs_511[5] -> SetParameters(100000,30.,2.);
  fitFuncs_511[7] = new TF1("fitFunc_511_ov7","gaus(0)",7.5,9.5);
  fitFuncs_511[7] -> SetParameters(100000,8.5,0.7);
  
  std::map<int,TF1*> fitFuncs_1275;
  fitFuncs_1275[3] = new TF1("fitFunc_1275_ov3","gaus(0)+gaus(3)+gaus(6)",13.,20.);
  fitFuncs_1275[3] -> SetParameters(20000,7.5,0.7,1000,10.,2.,500.,11.,1.);
  fitFuncs_1275[5] = new TF1("fitFunc_1275_ov5","gaus(0)+gaus(3)+gaus(6)",13.,20.);
  fitFuncs_1275[5] -> SetParameters(20000,13.5,0.7,1000,16.,2.,500.,18.,1.);
  fitFuncs_1275[7] = new TF1("fitFunc_1275_ov7","gaus(0)+gaus(3)+gaus(6)",18.,26.);
  fitFuncs_1275[7] -> SetParameters(20000,19.,0.7,1000,21.,2.,500.,25.,1.);
  
  TCanvas* c1 = new TCanvas("c1","c1");
  c1 -> SetLogy();
  
  int it = 1;
  for(auto mapIt : histos )
  {
    int ov = mapIt.first;
    
    histos[ov] -> SetLineColor(it);
    if( it == 1 ) histos[ov] -> Draw();
    else          histos[ov] -> Draw("same");
    
    histos[ov] -> Fit(fitFuncs_511[ov],"NRSM+");
    fitFuncs_511[ov] -> SetLineColor(it);
    fitFuncs_511[ov] -> Draw("same");
    
    histos[ov] -> Fit(fitFuncs_1275[ov],"NRSM+");
    fitFuncs_1275[ov] -> SetLineColor(it);
    fitFuncs_1275[ov] -> Draw("same");
    
    ++it;
  }
  
  /*
  std::map<int,TGraphErrors*> g_linearity;
  std::map<int,TGraphErrors*> g_linearity_rel;
  
  it = 1;
  for(auto mapIt : datas )
  {
    int ov = mapIt.first;
    
    
    g_linearity[ov] = new TGraphErrors();
    
    g_linearity[ov] -> SetPoint(g_linearity[ov]->GetN(),511,fitFuncs_511[ov]->GetParameter(1));
    g_linearity[ov] -> SetPointError(g_linearity[ov]->GetN()-1,0.,fitFuncs_511[ov]->GetParError(1));
    
    g_linearity[ov] -> SetPoint(g_linearity[ov]->GetN(),1275,fitFuncs_1275[ov]->GetParameter(1));
    g_linearity[ov] -> SetPointError(g_linearity[ov]->GetN()-1,0.,fitFuncs_1275[ov]->GetParError(1));
    
    g_linearity[ov] -> SetPoint(g_linearity[ov]->GetN(),1786,fitFuncs_1275[ov]->GetParameter(4));
    g_linearity[ov] -> SetPointError(g_linearity[ov]->GetN()-1,0.,fitFuncs_1275[ov]->GetParError(7));
    
    g_linearity[ov] -> SetMarkerStyle(20);
    g_linearity[ov] -> SetMarkerSize(0.7);
    g_linearity[ov] -> SetMarkerColor(it);
    
    
    g_linearity_rel[ov] = new TGraphErrors();
    
    g_linearity_rel[ov] -> SetPoint(g_linearity_rel[ov]->GetN(),511,fitFuncs_511[ov]->GetParameter(1)/fitFuncs_511[ov]->GetParameter(1));
    g_linearity_rel[ov] -> SetPointError(g_linearity_rel[ov]->GetN()-1,0.,fitFuncs_511[ov]->GetParError(1)/fitFuncs_511[ov]->GetParameter(1));
    
    g_linearity_rel[ov] -> SetPoint(g_linearity_rel[ov]->GetN(),1275,fitFuncs_1275[ov]->GetParameter(1)/fitFuncs_511[ov]->GetParameter(1));
    g_linearity_rel[ov] -> SetPointError(g_linearity_rel[ov]->GetN()-1,0.,fitFuncs_1275[ov]->GetParError(1)/fitFuncs_511[ov]->GetParameter(1));
    
    g_linearity_rel[ov] -> SetPoint(g_linearity_rel[ov]->GetN(),1786,fitFuncs_1275[ov]->GetParameter(4)/fitFuncs_511[ov]->GetParameter(1));
    g_linearity_rel[ov] -> SetPointError(g_linearity_rel[ov]->GetN()-1,0.,fitFuncs_1275[ov]->GetParError(4)/fitFuncs_511[ov]->GetParameter(1));
    
    g_linearity_rel[ov] -> SetMarkerStyle(20);
    g_linearity_rel[ov] -> SetMarkerSize(0.7);
    g_linearity_rel[ov] -> SetMarkerColor(it);
    
    
    ++it;
  }
  
  TCanvas* c2 = new TCanvas("c2","c2",1200,500);
  c2 -> Divide(2,1);
  
  c2 -> cd(1);
  
  TH1F* hpad1 = (TH1F*)( gPad->DrawFrame(0.,0.,2000.,25.) );
  hpad1 -> SetTitle(";energy [keV];amplitude [a.u.]");
  hpad1 -> Draw();
  
  for(auto mapIt : datas )
  {
    int ov = mapIt.first;
    g_linearity[ov] -> Draw("P,same");
  }
  
  c2 -> cd(2);
  
  TH1F* hpad2 = (TH1F*)( gPad->DrawFrame(0.,0.,2000.,5.) );
  hpad2 -> SetTitle(";energy [keV];relative amplitude [a.u.]");
  hpad2 -> Draw();
  
  it = 1;
  for(auto mapIt : datas )
  {
    int ov = mapIt.first;
    g_linearity_rel[ov] -> Draw("P,same");
    
    TF1* func = new TF1(Form("func_%d",ov),"[0]*(1-exp(-x/[1]))");
    func -> SetLineColor(it);
    g_linearity_rel[ov] -> Fit(func);
    
    ++it;
  }
  
  
  
  TGraphErrors* linearity_pe = new TGraphErrors();
  
  for(auto mapIt : datas )
  {
    int ov = mapIt.first;
    
    linearity_pe -> SetPoint(linearity_pe->GetN(),fitFuncs_511[ov]->GetParameter(1),40000*0.511*0.15*f_PDE->Eval(ov)/100.*f_gain->Eval(ov)*f_ENF->Eval(ov));
    linearity_pe -> SetPointError(linearity_pe->GetN()-1,fitFuncs_511[ov]->GetParError(1),0.);
    
    linearity_pe -> SetPoint(linearity_pe->GetN(),fitFuncs_1275[ov]->GetParameter(1),40000*1.275*0.15*f_PDE->Eval(ov)/100.*f_gain->Eval(ov)*f_ENF->Eval(ov));
    linearity_pe -> SetPointError(linearity_pe->GetN()-1,fitFuncs_1275[ov]->GetParError(1),0.);
    
    linearity_pe -> SetPoint(linearity_pe->GetN(),fitFuncs_1275[ov]->GetParameter(4),40000*1.786*0.15*f_PDE->Eval(ov)/100.*f_gain->Eval(ov)*f_ENF->Eval(ov));
    linearity_pe -> SetPointError(linearity_pe->GetN()-1,fitFuncs_1275[ov]->GetParError(4),0.);
  }
  
  TCanvas* c3 = new TCanvas("c3","c3");
  
  linearity_pe -> SetMarkerStyle(20);
  linearity_pe -> SetMarkerSize(0.7);
  linearity_pe -> Draw("AP");
  */
  
  // TF1* func = new TF1("func","[0]*(1-exp(-x/[1]))",0.,10000.);
  // func -> SetParameters(170.,40000.);
  // linearity_pe->Fit(func,"RSM+");
  // func -> Draw("same");
  
  // TFile* inFile_511  = TFile::Open("../plots/analyzeTOFPET2_wirelessBar2_511keV_plots.root", "READ");
  // TFile* inFile_1275 = TFile::Open("../plots/analyzeTOFPET2_wirelessBar2_1275keV_plots.root","READ");
  
  
  // std::map<int,TGraphErrors*> g_energy_vs_Vov;
  // g_energy_vs_Vov[511]  = (TGraphErrors*)( inFile_511 ->Get("g_energy_vs_Vov_barR_th30") );
  // g_energy_vs_Vov[1275] = (TGraphErrors*)( inFile_1275->Get("g_energy_vs_Vov_barR_th30") );
  
  // TGraphErrors* linearity = new TGraphErrors();
  
  // for(auto mapIt : g_energy_vs_Vov)
  // {
  //   double ov,amp;
  //   for(int point = 0; point < g_energy_vs_Vov[mapIt.first]->GetN(); ++point)
  //   {
  //     g_energy_vs_Vov[mapIt.first] -> GetPoint(point,ov,amp);
      
  //     double Eeff = mapIt.first * f_PDE->Eval(ov)/f_PDE->Eval(7.) * f_gain->Eval(ov)/f_gain->Eval(7.);
  //     double Npe = 40000. * mapIt.first/1000. * 0.15 * f_PDE->Eval(ov)/100.;
  //     // linearity -> SetPoint(linearity->GetN(),Eeff,amp);
  //     linearity -> SetPoint(linearity->GetN(),Npe,amp/f_gain->Eval(ov)*f_gain->Eval(2.));
  //   }
  // }
  
  // linearity -> SetMarkerStyle(20);
  // linearity -> SetMarkerSize(1);
  // linearity -> Draw("AP");
}
