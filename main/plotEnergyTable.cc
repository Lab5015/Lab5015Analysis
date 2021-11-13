using namespace std;

map<float,int> ovmap ;// = { {1,0}, {,1}, {,2}, {,3}, {,4}, {,5}, {,6}, {,7}, {,8} };

void plotEnergyTable(int run, int iCh1, int iCh2)
{

  gStyle->SetOptStat(0);

  int nCh [2] = {iCh1,iCh2};

  TCanvas* c[2][9];
  TH2F* h[2][9];
  
  auto fin = TFile::Open(Form("~/TOFHiR2X/sw_daq_tofhir2_Milan/sw_daq_tofhir2/data/reco/run%04i_e.root",run));
  auto data = (TTree*)fin->Get("data");

  vector<uint>* channelID;
  vector<float>* energy;
  float step1, step2;

  data->SetBranchAddress("channelID",&channelID);
  data->SetBranchAddress("energy",&energy);
  data->SetBranchAddress("step1",&step1);
  data->SetBranchAddress("step2",&step2);
  
  for(int i=0; i<9; ++i) {
    ovmap[1+i*0.5]=i;
    h[0][i] = new TH2F(Form("h%i%i",iCh1,i),Form("Events with valid energy - OV %.1f - ch %i;Att. gain;Delay E setting",1+i*0.5,iCh1),8,0,8,8,0,8);
    h[1][i] = new TH2F(Form("h%i%i",iCh2,i),Form("Events with valid energy - OV %.1f - ch %i;Att. gain;Delay E setting",1+i*0.5,iCh2),8,0,8,8,0,8);
  }

  for(int iEv=0; iEv<data->GetEntries(); ++iEv) {
    data->GetEntry(iEv);

    for (int iCh=0; iCh<energy->size(); ++iCh) {
      if (energy->at(iCh)<30 || energy->at(iCh)>950) continue;
      if (channelID->at(iCh)%32 == iCh1)
	h[0][ovmap[step1]]->Fill(step2/8,((int)step2)%8);
      if (channelID->at(iCh)%32 == iCh2)
	h[1][ovmap[step1]]->Fill(step2/8,((int)step2)%8);
    }

  }

  string plotDir = Form("/var/www/html/TOFHIR2X/MTDST_CERN_Oct21/EnergyCharacterization/run%04i/attVsDelayE",run);
  system(Form("mkdir -p %s",plotDir.c_str()));

  for(int i=0; i<9; ++i)
    for(int iCh=0; iCh<2; ++iCh) {
      c[iCh][i]= new TCanvas();
      h[iCh][i]->Draw("colz");
      c[iCh][i]->SaveAs(Form("%s/ov%.1f-ch%i.png",plotDir.c_str(),1+i*0.5,nCh[iCh]));
  }
}
