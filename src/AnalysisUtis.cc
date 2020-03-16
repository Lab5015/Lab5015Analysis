#include "interface/AnalysisUtils.h"



void TrackProcess(float* cpu, float* mem, float* vsz, float* rss)
{
  std::string dummy1, dummy2, dummy3, dummy4, dummy5;
  std::string time;
  
  //---get cpu/mem info                                                                                                                                                                                                                     
  int pid = getpid();
  std::string ps_command = "ps up "+std::string(Form("%d",pid))+" >.proc.tmp";
  system(ps_command.c_str());
  std::ifstream proc_tmp(".proc.tmp", std::ios::in);
  getline(proc_tmp, dummy1);
  proc_tmp >> dummy1 >> dummy2 >> cpu[0] >> mem[0] >> vsz[0] >> rss[0]
           >> dummy3 >> dummy4 >> dummy5 >> time;
  vsz[0] = vsz[0]/1000;
  rss[0] = rss[0]/1000;
  proc_tmp.close();
  if(cpu[0]>cpu[1])
    cpu[1] = cpu[0];
  if(mem[0]>mem[1])
    mem[1] = mem[0];
  if(vsz[0]>vsz[1])
    vsz[1] = vsz[0];
  if(rss[0]>rss[1])
    rss[1] = rss[0];
  
  //---print statistics                                                                                                                                                                                                                     
  std::cout << "\033[0m""-----Machine stats---current/max-----" << std::endl
            << "CPU(%): " << cpu[0] << "/" << cpu[1] << std::endl
            << "MEM(%): " << mem[0] << "/" << mem[1] << std::endl
            << "VSZ(M): " << vsz[0] << "/" << vsz[1] << std::endl
            << "RSS(M): " << rss[0] << "/" << rss[1] << std::endl
            << "time lasted: " << time << std::endl;
}



std::vector<std::string> GetTokens(const std::string& input, const char& sep)
{
  std::stringstream ss(input);
  std::string token;
  std::vector<std::string> tokens;
  
  while( std::getline(ss,token,sep) )
  {
    tokens.push_back(token);
  }
  
  return tokens;
}




float DeltaEta(const float& eta1, const float& eta2)
{
  return fabs( eta1 - eta2 );
}

float DeltaPhi(const float& phi1, const float& phi2)
{
  float dphi = fabs( phi1 - phi2 );
  if( dphi > PI ) dphi = 2*PI - dphi;
  return dphi;
}

float DeltaR(const float& eta1, const float& phi1,
             const float& eta2, const float& phi2)
{
  return sqrt( DeltaEta(eta1,eta2)*DeltaEta(eta1,eta2) +
               DeltaPhi(phi1,phi2)*DeltaPhi(phi1,phi2) );
}

float FindXMaximum(TH1F* histo, const float& xMin, const float& xMax)
{
  float max = -999999999.;
  int binMax = -1;
  for(int bin = 2; bin <= histo->GetNbinsX()-1; ++bin)
  {
    if( histo->GetBinCenter(bin) < xMin ) continue;
    if( histo->GetBinCenter(bin) > xMax ) continue;
    if( histo->GetBinContent(bin) > max )
    {
      float delta_min = ( histo->GetBinContent(bin-1) - histo->GetBinContent(bin) ) / histo->GetBinContent(bin);
      float delta_max = ( histo->GetBinContent(bin+1) - histo->GetBinContent(bin) ) / histo->GetBinContent(bin);
      if( fabs(delta_min) < 0.2 && fabs(delta_max) < 0.2 )
      {
        max = histo->GetBinContent(bin); binMax = bin;
      }
    };
  }
  return histo->GetBinCenter(binMax);
}
