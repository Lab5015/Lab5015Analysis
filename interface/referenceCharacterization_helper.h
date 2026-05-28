#ifndef REFERENCECHARACTERIZATION_HELPER_H
#define REFERENCECHARACTERIZATION_HELPER_H

#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <regex>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// ROOT
#include "TH1F.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"

#include "interface/AnalysisUtils.h"

struct EntryData {
  int index;
  float mpvL;
  float mpvR;
  float mpvLR;
  float slopeL;
  float slopeR;
  float offsetL;
  float offsetR;
  float slopeL_ref;
  float slopeR_ref;
  float offsetL_ref;
  float offsetR_ref;
  float sigmaLRef;
  float sigmaLR;
  float sigmaRRef;
  float sigmaRef;
  float sigmaL;
  float sigmaR;
  float sigmaDutRef;
};
struct SigmaResult {
  float sRef;
  float sL;
  float sR;
};
enum FitType { kGaussian, kLandau };

SigmaResult ilTriangoloNo(float sigma12, float sigma23, float sigma13, float sigmaCorr = 0);

TH1F* AutoRebinAndRange(TH1F* h, double nSigma = 5.);
void  AutoRangeRMS(TH1F* h, double nSigma = 5.);

void SaveHistoToCanvas(TFile* outFile, TH1F* histo, const std::string& plotdirectory);
void SaveHisto2ToCanvas(TFile* outFile, TH2* histo, const std::string& plotdirectory);
void SaveProfileToCanvas(TFile* outFile, TProfile* prof, const std::string& plotdirectory);

TF1* FitAndSaveHisto(TFile* outFile, TH1F* histo, const std::string& plotdirectory, double nSigmaLow = 1., double nSigmaUp = 1., int saveFlag=1, FitType fitType = kGaussian);
TF1* FitAndSaveProfile(TFile* outFile, TProfile* prof, const std::string& plotdirectory, int saveFlag=1);

bool PassSelection(ModuleEventWithRefClass* anEvent, double deltaTL, double deltaTR, double deltaTL_AveRef, double deltaTR_AveRef, std::map<std::string, std::map<int, std::vector<float>*> > ranges, int index1, double energyL, double energyR);

std::string GetQuantityFromFilename(const std::string& name);
void SortPlotsByQuantity(const std::string& plotDir);

#endif // REFERENCECHARACTERIZATION_HELPER_H
