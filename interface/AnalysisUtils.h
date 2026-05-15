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
#include <tuple>

#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"

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

class ModuleEventWithRefClass : public TObject {
public:
  int barID;
  float Vov;
  int vth;
  // DUT bar
  float energyL;
  float energyR;
  float totL;
  float totR;
  long long timeL;
  long long timeR;
  unsigned short t1fineL;
  unsigned short t1fineR;
  float qT1L;
  float qT1R;
  // reference bar
  float energyL_ref;
  float energyR_ref;
  float totL_ref;
  float totR_ref;
  long long timeL_ref;
  long long timeR_ref;
  long long timeL_ref_cor;
  long long timeR_ref_cor;
  unsigned short t1fineL_ref;
  unsigned short t1fineR_ref;
  float qT1L_ref;
  float qT1R_ref;
  // tracking info
  int nhits;
  float x;
  float y;
  ClassDef(ModuleEventWithRefClass,1);
};

class ModuleEventClass : public TObject {
public:
  int barID;
  float Vov;
  int vth;
  float energyL;
  float energyR;
  float totL;
  float totR;
  long long timeL;
  long long timeR;
  unsigned short t1fineL;
  unsigned short t1fineR;
  float qT1L;
  float qT1R;
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
float DeltaR(const float& eta1, const float& phi1, const float& eta2, const float& phi2);
float FindXMaximum(TH1F* histo, const float& xMin, const float& xMax, const bool& checkDerivative = false);
int FindBin(const float& val, const std::vector<float>* ranges);
std::map<std::tuple<int,std::string,int,int>, float> loadCalibrationMap(const std::string& filename);
Double_t langaufun(Double_t *x, Double_t *par);
void GetEnergyBins(TH1F *h, std::vector<float> *r, std::map<int, float> & b);
void drawDeltaT(TCanvas *& c, TH1F *histo, TF1 *& fitFunc, std::string xaxis_label, std::string latex_label, std::string drawSame);
int findBarAndSide(int ch, int asic, std::vector<unsigned int>& channelMapping, std::string& side);
using CalibKey = std::tuple<int,std::string,int,int>;
using CalibMap = std::map<CalibKey, float>;
void applyCalibration(float& energy, int bar, const std::string& side, float Vov, float vth, const CalibMap& calibMap);
#endif
