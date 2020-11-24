void createHistograms()
{
  TFile* outFile = TFile::Open("histograms.root","RECREATE");
  
  TFile* inFile;
  TTree* data;
  
  std::vector<int> channels;
  channels.push_back(195);
  channels.push_back(217);
  
  
  inFile = TFile::Open("/data/TOFPET2/reco/run1496_ped_t.root","READ");
  data = (TTree*)( inFile->Get("data") );
  
  for(auto ch : channels)
  {
    TH1F* histo = new TH1F(Form("h1_energy_ov2_ch%d",ch),"",200,0.,200.);
    data -> Draw(Form("energy[%d]>>h1_energy_ov2_ch%d",ch,ch),"","goff");
    
    outFile -> cd();
    histo -> Write();
  }
  
  
  inFile = TFile::Open("/data/TOFPET2/reco/run1497_ped_t.root","READ");
  data = (TTree*)( inFile->Get("data") );
  
  for(auto ch : channels)
  {
    TH1F* histo = new TH1F(Form("h1_energy_ov4_ch%d",ch),"",200,0.,200.);
    data -> Draw(Form("energy[%d]>>h1_energy_ov4_ch%d",ch,ch),"","goff");
    
    outFile -> cd();
    histo -> Write();
  }
  
  
  inFile = TFile::Open("/data/TOFPET2/reco/run1498_ped_t.root","READ");
  data = (TTree*)( inFile->Get("data") );
  
  for(auto ch : channels)
  {
    TH1F* histo = new TH1F(Form("h1_energy_ov7_ch%d",ch),"",200,0.,200.);
    data -> Draw(Form("energy[%d]>>h1_energy_ov7_ch%d",ch,ch),"","goff");
    
    outFile -> cd();
    histo -> Write();
  }
}
