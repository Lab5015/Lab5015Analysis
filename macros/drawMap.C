
std::vector<int> unpack_range(std::string range)
{
  std::vector<int> unpacked_list;

  std::stringstream ss(range);
  std::string token;
  while( std::getline(ss,token,',') )
  {
    std::stringstream ss2(token);
    std::string token2;
    int min = -1;
    int max = -1;
    while( std::getline(ss2,token2,'-') )
    {
      if( min != -1 && max == -1 ) max = atoi(token2.c_str());
      if( min == -1 ) min = atoi(token2.c_str());
    }
    if( max == -1 ) max = min;

    for(int ii=min; ii<=max; ++ii)
    {
      unpacked_list.push_back(ii);
    }
  }
  return unpacked_list;
}


void drawMap(const std::string& fileName, int chipID, std::string chList, float vov=-1, int ith1=-1, int ith2=-1)
{
  TFile* inFile = TFile::Open(fileName.c_str(),"READ");
  TTree* tree = (TTree*)( inFile->Get("data") );


  std::map<int, std::pair<int, int>> channelMapping; //x,y 
  int firstCh = chipID*32;


  //left
  channelMapping[firstCh+14] = std::make_pair(2,0);
  channelMapping[firstCh+12] = std::make_pair(2,1);
  channelMapping[firstCh+10] = std::make_pair(2,2);
  channelMapping[firstCh+8]  = std::make_pair(2,3);
  channelMapping[firstCh+6]  = std::make_pair(2,4);
  channelMapping[firstCh+4]  = std::make_pair(2,5);
  channelMapping[firstCh+2]  = std::make_pair(2,6);
  channelMapping[firstCh+0]  = std::make_pair(2,7);
  channelMapping[firstCh+1]  = std::make_pair(2,8);
  channelMapping[firstCh+3]  = std::make_pair(2,9);
  channelMapping[firstCh+5]  = std::make_pair(2,10);
  channelMapping[firstCh+7]  = std::make_pair(2,11);
  channelMapping[firstCh+9]  = std::make_pair(2,12);
  channelMapping[firstCh+11] = std::make_pair(2,13);
  channelMapping[firstCh+13] = std::make_pair(2,14);
  channelMapping[firstCh+15] = std::make_pair(2,15);
  //right
  channelMapping[firstCh+17] = std::make_pair(0,0);
  channelMapping[firstCh+19] = std::make_pair(0,1);
  channelMapping[firstCh+21] = std::make_pair(0,2);
  channelMapping[firstCh+23]  = std::make_pair(0,3);
  channelMapping[firstCh+25]  = std::make_pair(0,4);
  channelMapping[firstCh+27]  = std::make_pair(0,5);
  channelMapping[firstCh+29]  = std::make_pair(0,6);
  channelMapping[firstCh+31]  = std::make_pair(0,7);
  channelMapping[firstCh+30]  = std::make_pair(0,8);
  channelMapping[firstCh+28]  = std::make_pair(0,9);
  channelMapping[firstCh+26]  = std::make_pair(0,10);
  channelMapping[firstCh+24]  = std::make_pair(0,11);
  channelMapping[firstCh+22]  = std::make_pair(0,12);
  channelMapping[firstCh+20] = std::make_pair(0,13);
  channelMapping[firstCh+18] = std::make_pair(0,14);
  channelMapping[firstCh+16] = std::make_pair(0,15);
  //left
  channelMapping[firstCh+46] = std::make_pair(2,17);
  channelMapping[firstCh+44] = std::make_pair(2,18);
  channelMapping[firstCh+42] = std::make_pair(2,19);
  channelMapping[firstCh+40]  = std::make_pair(2,20);
  channelMapping[firstCh+38]  = std::make_pair(2,21);
  channelMapping[firstCh+36]  = std::make_pair(2,22);
  channelMapping[firstCh+34]  = std::make_pair(2,23);
  channelMapping[firstCh+32]  = std::make_pair(2,24);
  channelMapping[firstCh+33]  = std::make_pair(2,25);
  channelMapping[firstCh+35]  = std::make_pair(2,26);
  channelMapping[firstCh+37]  = std::make_pair(2,27);
  channelMapping[firstCh+39]  = std::make_pair(2,28);
  channelMapping[firstCh+41]  = std::make_pair(2,29);
  channelMapping[firstCh+43] = std::make_pair(2,30);
  channelMapping[firstCh+45] = std::make_pair(2,31);
  channelMapping[firstCh+47] = std::make_pair(2,32);
  //right
  channelMapping[firstCh+49] = std::make_pair(0,17);
  channelMapping[firstCh+51] = std::make_pair(0,18);
  channelMapping[firstCh+53] = std::make_pair(0,19);
  channelMapping[firstCh+55]  = std::make_pair(0,20);
  channelMapping[firstCh+57]  = std::make_pair(0,21);
  channelMapping[firstCh+59]  = std::make_pair(0,22);
  channelMapping[firstCh+61]  = std::make_pair(0,23);
  channelMapping[firstCh+63]  = std::make_pair(0,24);
  channelMapping[firstCh+62]  = std::make_pair(0,25);
  channelMapping[firstCh+60]  = std::make_pair(0,26);
  channelMapping[firstCh+58]  = std::make_pair(0,27);
  channelMapping[firstCh+56]  = std::make_pair(0,28);
  channelMapping[firstCh+54]  = std::make_pair(0,29);
  channelMapping[firstCh+52] = std::make_pair(0,30);
  channelMapping[firstCh+50] = std::make_pair(0,31);
  channelMapping[firstCh+48] = std::make_pair(0,32);



  TProfile2D* totProfArray = new TProfile2D("totProfArray","totProfArray; x; y; <tot>", 3,0,3, 33,0,33);
  std::vector<int> chVec   = unpack_range(chList);

  std::map<int, TH1F*> totMap;
  std::map<int, TH1F*> totHisto;

  int colorItr = 1;
  for(auto ch : chVec)
  {
    char* hName = Form("tot_ch_%d", ch);
    totMap[ch] = new TH1F(hName, hName,1000,0,30);
    totMap[ch]->SetLineColor(kRainBow+5*colorItr);
    totMap[ch]->SetLineWidth(3);

    ++colorItr;
  }

  float tot, step2, step1;
  unsigned int channelID;

  tree->SetBranchAddress("tot", &tot);
  tree->SetBranchAddress("channelID", &channelID);
  tree->SetBranchAddress("step2", &step2);
  tree->SetBranchAddress("step1", &step1);

  
  for(int entry = 0; entry < tree->GetEntries(); ++entry)
  {
    tree -> GetEntry(entry);
    if(tot/1000 > 300 || tot/1000 < 0)
      continue;

    if(vov != -1 && vov != step1)
      continue;

    if(ith1 != -1 && ith1 != float(int(step2/10000)-1))
      continue;

    if(ith2 != -1 && ith2 != int((step2-10000*(ith1+1))/100.)-1)
      continue;

    totProfArray->Fill(channelMapping[channelID].first, channelMapping[channelID].second , tot/1000);
    
    for(auto ch : chVec)
    {
      if(channelID == ch)
      	totMap[channelID]->Fill(tot/1000);
    }
  }


  TCanvas* c0 = new TCanvas();
  totProfArray->Draw("COLZTEXT");

  TCanvas* c2 = new TCanvas();
  auto it = totMap.begin();
  it->second->DrawNormalized("HISTO");
  for(it = std::next(totMap.begin()); it!=totMap.end(); ++it)
  {
    it->second->DrawNormalized("HISTO,sames");
  }

  c2->BuildLegend(0.48,0.68,0.88,0.88);

  return;
}
