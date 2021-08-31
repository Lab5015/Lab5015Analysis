#include "interface/TOFHIRThresholdZero.h"


TOFHIRThresholdZero::TOFHIRThresholdZero(const std::string& cfgFileName, const int& verbosity)
{
  std::cout << cfgFileName << std::endl;
  
  std::ifstream cfgFile(cfgFileName.c_str(),std::ios::in);

  if( !cfgFile.is_open() )
    {
      std::cout << "ERROR: cannot open file " << cfgFileName << std::endl;
      return;
    }
  
  std::string line;
  std::string token;
  while( getline(cfgFile,line) )
    {
      if( line.at(0) == '#' )
	continue;
      
      std::stringstream ss(line);
      
      std::vector<std::string> tokens;
      while( std::getline(ss,token,'\t') )
	{
	  tokens.push_back(token);
	  if( verbosity ) std::cout << token << " - ";
	}
      if( verbosity ) std::cout << std::endl;
      
      int asic = atoi(tokens[2].c_str());
      int ch = atoi(tokens[3].c_str());
      float zero_th1 = atof(tokens[6].c_str());
      float zero_th2 = atof(tokens[7].c_str());
      float zero_thE = atof(tokens[8].c_str());
      if( verbosity ) std::cout << asic << " " << ch << " " << zero_th1 << " " << zero_th2 << " " << zero_thE << std::endl;
      
      std::tuple<float,float,float> zero(zero_th1,zero_th2,zero_thE);
      _zeros[ch+asic*32] = zero;
    }
}



float TOFHIRThresholdZero::GetThresholdZero(const int& ch, const std::string& thr)
{
  float value = 63.;
  
  std::tuple<float,float,float> aTuple = _zeros[ch];
  
  if( thr.find("ith1") != std::string::npos || thr.find("vth1") != std::string::npos )
    value = std::get<0>(aTuple);
  
  if( thr.find("ith2") != std::string::npos || thr.find("vth2") != std::string::npos )
    value = std::get<1>(aTuple);
  
  if( thr.find("ithE") != std::string::npos || thr.find("vthE") != std::string::npos )
    value = std::get<2>(aTuple);
  
  return value;
}
