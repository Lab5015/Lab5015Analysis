#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>



class TOFHIRThresholdZero {

 public:
  TOFHIRThresholdZero(const std::string& cfgFileName, const int& verbosity = 0);
  ~TOFHIRThresholdZero() {};
  
  float GetThresholdZero(const int& ch, const std::string& thr);
  
 private:
  std::map<int,std::tuple<float,float,float> > _zeros;
};
