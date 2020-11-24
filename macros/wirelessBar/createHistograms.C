void createHistograms()
{
  TFile* outFile = TFile::Open("histograms_new.root","RECREATE");
  
  TFile* inFile;
  TTree* data;
  
  std::vector<int> channels;
  channels.push_back(49);
  channels.push_back(249);
  
  
  inFile = TFile::Open("/data/TOFPET2/reco/run1170_ped_t.root","READ");
  data = (TTree*)( inFile->Get("data") );
  
  for(auto ch : channels)
  {
    TH1F* histo = new TH1F(Form("h1_energy_ov3_ch%d",ch),"",200,0.,200.);
    data -> Draw(Form("energy[49]>>h1_energy_ov3_ch%d",ch),"","goff");
    
    outFile -> cd();
    histo -> Write();
  }
  
  
  inFile = TFile::Open("/data/TOFPET2/reco/run1178_ped_t.root","READ");
  data = (TTree*)( inFile->Get("data") );
  
  for(auto ch : channels)
  {
    TH1F* histo = new TH1F(Form("h1_energy_ov5_ch%d",ch),"",200,0.,200.);
    data -> Draw(Form("energy[49]>>h1_energy_ov5_ch%d",ch),"","goff");
    
    outFile -> cd();
    histo -> Write();
  }
  
  
  inFile = TFile::Open("/data/TOFPET2/reco/run1179_ped_t.root","READ");
  data = (TTree*)( inFile->Get("data") );
  
  for(auto ch : channels)
  {
    TH1F* histo = new TH1F(Form("h1_energy_ov7_ch%d",ch),"",200,0.,200.);
    data -> Draw(Form("energy[49]>>h1_energy_ov7_ch%d",ch),"","goff");
    
    outFile -> cd();
    histo -> Write();
  }
}
