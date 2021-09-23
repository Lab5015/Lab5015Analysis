#ifndef ANALYSIS_UTILS_H
#define ANALYSIS_UTILS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include "TH1F.h"

#define PI 3.14159265359



class EventClass : public TObject {

public:
  std::string stepLabel;
  std::string ch1;
  std::string ch2;
  std::string label1;
  std::string label2;
  std::string label12;
  float x;
  float y;
  int isBar1;
  int isBar2;
  int isBarSide1;
  int isBarSide2;
  int isHorizontal1;
  int isHorizontal2;
  float qfine1;
  float qfine1L;
  float qfine1R;
  float qfine2;
  float qfine2L;
  float qfine2R;
  float tot1;
  float tot1L;
  float tot1R;
  float tot2;
  float tot2L;
  float tot2R;
  float energy1;
  float energy1L;
  float energy1R;
  float energy2;
  float energy2L;
  float energy2R;
  long long time1;
  long long time2;
  float t1fine1;
  float t1fine1L;
  float t1fine1R;
  float t1fine2;
  float t1fine2L;
  float t1fine2R;
  
  ClassDef(EventClass,1);
};



class ModuleEventClass : public TObject {

public:
  int barID;
  float Vov;
  int vth1;
  float energyL;
  float energyR;
  float totL;
  float totR;
  long long timeL;
  long long timeR;
  unsigned short t1fineL;
  unsigned short t1fineR;
  int nhits;
  float x;
  float y;
  ClassDef(ModuleEventClass,1);
};


struct Event
{
  std::string stepLabel;
  std::string ch1;
  std::string ch2;
  std::string label1;
  std::string label2;
  std::string label12;
  float x;
  float y;
  int isBar1;
  int isBar2;
  int isHorizontal1;
  int isHorizontal2;
  float qfine1;
  float qfine1L;
  float qfine1R;
  float qfine2;
  float qfine2L;
  float qfine2R;
  float tot1;
  float tot1L;
  float tot1R;
  float tot2;
  float tot2L;
  float tot2R;
  float energy1;
  float energy1L;
  float energy1R;
  float energy2;
  float energy2L;
  float energy2R;
  long long time1;
  long long time2;
  float t1fine1;
  float t1fine1L;
  float t1fine1R;
  float t1fine2;
  float t1fine2L;
  float t1fine2R;
};



void TrackProcess(float* cpu, float* mem, float* vsz, float* rss);

std::vector<std::string> GetTokens(const std::string& input, const char& sep);

float DeltaEta(const float& eta1, const float& eta2);
float DeltaPhi(const float& phi1, const float& phi2);
float DeltaR(const float& eta1, const float& phi1,
             const float& eta2, const float& phi2);

float FindXMaximum(TH1F* histo, const float& xMin, const float& xMax, const bool& checkDerivative = false);

int FindBin(const float& val, const std::vector<float>* ranges);

#endif
