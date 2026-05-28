#include "interface/referenceCharacterization_helper.h"

SigmaResult ilTriangoloNo(float sigma12, float sigma23, float sigma13, float sigmaCorr) {
  SigmaResult res;
  res.sRef = sqrt(0.5 * (sigma12*sigma12 + sigma13*sigma13 - sigma23*sigma23 - sigmaCorr*sigmaCorr));
  res.sL = sqrt(0.5 * (sigma12*sigma12 + sigma23*sigma23 - sigma13*sigma13 - sigmaCorr*sigmaCorr));
  res.sR = sqrt(0.5 * (sigma13*sigma13 + sigma23*sigma23 - sigma12*sigma12) - sigmaCorr*sigmaCorr);
  return res;
}

TH1F* AutoRebinAndRange(TH1F* h, double nSigma)
{
    int nBinsOrig = h->GetNbinsX();
    double mean = h->GetMean();
    double rms  = h->GetRMS();
    double xMin = mean - nSigma * rms;
    double xMax = mean + nSigma * rms;
    double targetBinWidth = rms / 5.;
    double currentBinWidth = (h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin()) / nBinsOrig;
    int rebinFactor = std::round(targetBinWidth / currentBinWidth);
    if (rebinFactor < 1) rebinFactor = 1;
    while (nBinsOrig % rebinFactor != 0)
        rebinFactor--;
    if (rebinFactor < 1) rebinFactor = 1;
    TH1F* hRebinned = (TH1F*)h->Rebin(rebinFactor, Form("%s_rebinned_tmp", h->GetName()));
    // new histo with reduced range
    int nBinsNew = (xMax - xMin) / (rebinFactor * currentBinWidth);
    if (nBinsNew < 10) nBinsNew = 10;
    TH1F* hFinal = new TH1F(Form("%s_final", h->GetName()), h->GetTitle(), nBinsNew, xMin, xMax);
    for (int i = 1; i <= hRebinned->GetNbinsX(); i++)
    {
        double x = hRebinned->GetBinCenter(i);
        double content = hRebinned->GetBinContent(i);
        if (x >= xMin && x <= xMax)
            hFinal->Fill(x, content);
    }
    return hFinal;
}

void AutoRangeRMS(TH1F* h, double nSigma)
{
    double mean = h->GetMean();
    double rms  = h->GetRMS();
    h->GetXaxis()->SetRangeUser(mean - nSigma*rms, mean + nSigma*rms);
}

void SaveHistoToCanvas(TFile* outFile, TH1F* histo, const std::string& plotdirectory) {
  if (!histo){
    std::cerr << "[SKIP] Histo null: " << std::endl;
    return;
  }
  outFile->cd();
  AutoRangeRMS(histo);
  TCanvas* c = new TCanvas(Form("c_%s", histo->GetName()), histo->GetName(), 600, 500);
  histo->Draw();    
  histo->Write();
  c->SaveAs(Form("%s/%s.png", plotdirectory.c_str(), c->GetName()));
  delete c;
}
void SaveHisto2ToCanvas(TFile* outFile, TH2* histo, const std::string& plotdirectory) {
  if (!histo){
    std::cerr << "[SKIP] Histo null: " << std::endl;
    return;
  }
  outFile->cd();
  TCanvas* c = new TCanvas(Form("c_%s", histo->GetName()), histo->GetName(), 600, 500);
  histo->Draw("COLZ");    
  histo->Write();
  c->SaveAs(Form("%s/%s.png", plotdirectory.c_str(), c->GetName()));
  delete c;
}
void SaveProfileToCanvas(TFile* outFile, TProfile* prof, const std::string& plotdirectory) {
  if (!prof) {
    std::cerr << "[SKIP] Profile null: " << std::endl;
    return;
  }
  if (prof->GetEntries() < 10) {
    std::cerr << "[SKIP] Few entries (" << prof->GetEntries() << "): " <<prof->GetName() << std::endl;
    return;
  }
  // - compute global Y mean/RMS to identify bulk region
  double meanY = prof->GetMean(2);
  double rmsY  = prof->GetRMS(2);
  double nsigma = 2.0;
  int firstFitBin = -1;
  int lastFitBin  = -1;
  int goodBins    = 0;
  for (int i = 1; i <= prof->GetNbinsX(); ++i) {
    if (prof->GetBinEntries(i) < 5)
      continue;
    double y = prof->GetBinContent(i);
    if (fabs(y - meanY) > nsigma * rmsY)
      continue;
    if (firstFitBin == -1)
      firstFitBin = i;
    lastFitBin = i;
    goodBins++;
  }
  if (goodBins < 3) {
    std::cerr << "[SKIP] Too few valid bins after selection: " << prof->GetName() << std::endl;
    return;
  }
  double xmin = prof->GetXaxis()->GetBinLowEdge(firstFitBin);
  double xmax = prof->GetXaxis()->GetBinUpEdge(lastFitBin);
  if (xmin >= xmax) {
    std::cerr << "[SKIP] Bad fit range for " << prof->GetName() << std::endl;
    return;
  }
  TCanvas* c = new TCanvas(Form("c_%s", prof->GetName()), prof->GetName(), 600, 500);
  prof->Draw();
  c->SaveAs(Form("%s/%s.png", plotdirectory.c_str(), c->GetName()));
  outFile->cd();
  prof->Write();
  delete c;
}

TF1* FitAndSaveHisto(TFile* outFile, TH1F* histo, const std::string& plotdirectory, double nSigmaLow, double nSigmaUp, int saveFlag, FitType fitType) {
  if (!histo){
    std::cerr << "[SKIP] Histo null: " << std::endl;
    return nullptr;
  }
  if (histo->GetEntries()<10)
    {
      std::cerr << "[SKIP] Few entries " << histo->GetName() << std::endl;
      return nullptr;
    }
  outFile->cd();
  AutoRangeRMS(histo);
  
  int peakBin = histo->GetMaximumBin();   
  double peakX = histo->GetXaxis()->GetBinCenter(peakBin);
  double sigmaGuess = histo->GetRMS();
  double xmin = peakX - nSigmaLow * sigmaGuess;
  double xmax = peakX + nSigmaUp * sigmaGuess;
  // safety bounds
  if (xmin < histo->GetXaxis()->GetXmin())
    xmin = histo->GetXaxis()->GetXmin();
  if (xmax > histo->GetXaxis()->GetXmax())
    xmax = histo->GetXaxis()->GetXmax();
  if (xmin >= xmax) {
    std::cerr << "[SKIP] Bad range: " << histo->GetName() << std::endl;
    return nullptr;
  }
  const char* funcName = nullptr;
  const char* funcExpr = nullptr;
  switch (fitType) {
  case kGaussian:
    funcName = "gaus";
    funcExpr = "gaus";
    break;
  case kLandau:
    funcName = "landau";
    funcExpr = "landau";
    break;
  }
  TCanvas* c = new TCanvas(Form("c_%s", histo->GetName()), histo->GetName(), 600, 500);
  TF1* f = new TF1(Form("f_%s_%s", histo->GetName(), funcName), funcExpr, xmin, xmax);
  f->SetParameters(histo->GetMaximum(), peakX, 20);
  histo->Fit(f, "QRS");
  histo->Draw();
  c->SaveAs(Form("%s/%s.png", plotdirectory.c_str(), c->GetName()));
  if (saveFlag == 1) {
    outFile->cd();
    histo->Write();
    f->Write();
  }
  delete c;
  return f;
}

TF1* FitAndSaveProfile(TFile* outFile, TProfile* prof, const std::string& plotdirectory, int saveFlag) {
  if (!prof) {
    std::cerr << "[SKIP] Profile null: " << std::endl;
    return nullptr;
  }
  if (prof->GetEntries() < 10) {
    std::cerr << "[SKIP] Few entries (" << prof->GetEntries() << "): " <<prof->GetName() << std::endl;
    return nullptr;
  }
  // - compute global Y mean/RMS to identify bulk region
  double meanY = prof->GetMean(2);
  double rmsY  = prof->GetRMS(2);
  double nsigma = 2.0;
  int firstFitBin = -1;
  int lastFitBin  = -1;
  int goodBins    = 0;
  for (int i = 1; i <= prof->GetNbinsX(); ++i) {
    if (prof->GetBinEntries(i) < 5)
      continue;
    double y = prof->GetBinContent(i);
    if (fabs(y - meanY) > nsigma * rmsY)
      continue;
    if (firstFitBin == -1)
      firstFitBin = i;
    lastFitBin = i;
    goodBins++;
  }
  if (goodBins < 3) {
    std::cerr << "[SKIP] Too few valid bins after selection: " << prof->GetName() << std::endl;
    return nullptr;
  }
  double xmin = prof->GetXaxis()->GetBinLowEdge(firstFitBin);
  double xmax = prof->GetXaxis()->GetBinUpEdge(lastFitBin);
  if (xmin >= xmax) {
    std::cerr << "[SKIP] Bad fit range for " << prof->GetName() << std::endl;
    return nullptr;
  }
  TCanvas* c = new TCanvas(Form("c_%s", prof->GetName()), prof->GetName(), 600, 500);
  TF1* f = new TF1(Form("f_%s", prof->GetName()), "pol1", xmin, xmax);
  prof->Fit(f, "QRS");
  prof->Draw();
  c->SaveAs(Form("%s/%s.png", plotdirectory.c_str(), c->GetName()));
  if (saveFlag == 1) {
    outFile->cd();
    prof->Write();
    f->Write();
  }
  delete c;
  return f;
}

bool PassSelection(ModuleEventWithRefClass* anEvent, double deltaTL, double deltaTR, double deltaTL_AveRef, double deltaTR_AveRef, std::map<std::string, std::map<int, std::vector<float>*> > ranges, int index1, double energyL, double energyR)
{
    if (std::abs(deltaTL) > 5000 || std::abs(deltaTR) > 5000 || std::abs(deltaTL_AveRef) > 5000 || std::abs(deltaTR_AveRef) > 5000)
        return false;
    if (anEvent->totL < 0 || anEvent->totL > 25 || anEvent->totR < 0 || anEvent->totR > 25 || anEvent->totL_ref < 0 || anEvent->totL_ref > 25 || anEvent->totR_ref < 0 || anEvent->totR_ref > 25)
        return false;
    if (!ranges.at("L-R")[index1] || !ranges.at("L")[index1] || !ranges.at("R")[index1])
        return false;
    auto it = ranges.at("L-R").find(index1);
    if (it == ranges.at("L-R").end() || !it->second || it->second->size() < 2)
        return false;
    int energyBinAverage = FindBin(0.5 * (anEvent->energyL + anEvent->energyR), ranges.at("L-R")[index1]) + 1;
    if (energyBinAverage < 1)
        return false;
    double energyAve = 0.5*(energyL+energyR);   
    if (energyAve < ranges.at("L-R")[index1]->at(0) || energyAve > ranges.at("L-R")[index1]->at(1))
        return false;
    return true;
}

std::string GetQuantityFromFilename(const std::string& name)
{
    // name of the kind of c_<obj>_<quantity>_barXX...
    std::regex re("^c_[^_]+_(.+?)_bar");
    std::smatch match;
    if (!std::regex_search(name, match, re))
        return "UNKNOWN";
    std::string quantity = match[1];
    // remove vs_...
    size_t pos = quantity.find("vs_");
    if (pos != std::string::npos)
        quantity = quantity.substr(0, pos);
    // remove underscores
    while (!quantity.empty() && quantity.back() == '_')
        quantity.pop_back();
    return quantity;
}

void SortPlotsByQuantity(const std::string& plotDir)
{
    DIR* dir = opendir(plotDir.c_str());
    if (!dir) {
      std::cerr << "Cannot open directory: " << plotDir << std::endl;
      return; }
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr)
    {
        std::string fileName = entry->d_name;
        if (fileName == "." || fileName == "..")
            continue;
        if (fileName.find(".png") == std::string::npos)
            continue;
        std::string quantity = GetQuantityFromFilename(fileName);
        std::string targetDir = plotDir + "/" + quantity;
        mkdir(targetDir.c_str(), 0755);
        std::string oldPath = plotDir + "/" + fileName;
        std::string newPath = targetDir + "/" + fileName;
        rename(oldPath.c_str(), newPath.c_str());
    }
    closedir(dir);
}
