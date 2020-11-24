void createHistograms()
{
  TFile* outFile = TFile::Open("histograms_new.root","RECREATE");
  
  TFile* inFile;
  TTree* data;
  
  std::vector<int> channels;
  channels.push_back(49);
  channels.push_back(249);
  
  
  inFile = TFile::Open("/data/TOFPET2/reco/run1540_ped_t.root","READ");
  data = (TTree*)( inFile->Get("data") );
  
  for(auto ch : channels)
  {
    TH1F* histo = new TH1F(Form("h1_energy_ov2_ch%d",ch),"",200,0.,200.);
    data -> Draw(Form("energy[%d]>>h1_energy_ov2_ch%d",ch,ch),"","goff");
    
    outFile -> cd();
    histo -> Write();
  }
  
  
  inFile = TFile::Open("/data/TOFPET2/reco/run1541_ped_t.root","READ");
  data = (TTree*)( inFile->Get("data") );
  
  for(auto ch : channels)
  {
    TH1F* histo = new TH1F(Form("h1_energy_ov4_ch%d",ch),"",200,0.,200.);
    data -> Draw(Form("energy[%d]>>h1_energy_ov4_ch%d",ch,ch),"","goff");
    
    outFile -> cd();
    histo -> Write();
  }
  
  
  inFile = TFile::Open("/data/TOFPET2/reco/run1542_ped_t.root","READ");
  data = (TTree*)( inFile->Get("data") );
  
  for(auto ch : channels)
  {
    TH1F* histo = new TH1F(Form("h1_energy_ov7_ch%d",ch),"",200,0.,200.);
    data -> Draw(Form("energy[%d]>>h1_energy_ov7_ch%d",ch,ch),"","goff");
    
    outFile -> cd();
    histo -> Write();
  }
}
