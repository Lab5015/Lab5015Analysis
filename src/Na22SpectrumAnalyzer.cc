#include "interface/Na22SpectrumAnalyzer.h"



std::map<std::string,std::pair<float,float> > Na22SpectrumAnalyzer(TH1F* histo,
                                                                   std::vector<float>* ranges)
{
  gStyle -> SetOptFit(0);
  
  
  std::map<std::string,std::pair<float,float> > res;
  
  /*
  //--- find first peak (compton + threshold)
  float firstPeakVal = -1;
  float firstPeakPos = -1;
  int firstPeakCount = 0;
  bool firstPeak = false;
  for(int bin = 1; bin <= histo->GetNbinsX(); ++bin)
  {
    float binContent = histo->GetBinContent(bin);
    float binCenter = histo->GetBinCenter(bin);
    
    if( !firstPeak )
    {
      if( binContent >= firstPeakVal && binContent >= 1000. )
      {
        firstPeakVal = binContent;
        firstPeakPos = binCenter;
        firstPeakCount = 0;
      }
      else if( binContent < firstPeakVal )
        ++firstPeakCount;
    }
    if( firstPeakCount == 20 )
      firstPeak = true;
  }
  // TLine* line = new TLine(firstPeakPos,0.,firstPeakPos,histo->GetMaximum());
  // line -> Draw("same");
  
  
  //--- find first dip (just before 511 peak)
  float firstDipVal = 999999999;
  float firstDipPos = -1;
  int firstDipCount = 0;
  bool firstDip = false;
  for(int bin = 1; bin <= histo->GetNbinsX(); ++bin)
  {
    float binContent = histo->GetBinContent(bin);
    float binCenter = histo->GetBinCenter(bin);
    if( binCenter <= firstPeakPos ) continue;
    
    if( !firstDip )
    {
      if( binContent <= firstDipVal && binContent >= 100. )
      {
        firstDipVal = binContent;
        firstDipPos = binCenter;
        firstDipCount = 0;
      }
      else if( binContent > firstDipVal )
        ++firstDipCount;
    }
    if( firstDipCount == 10 )
      firstDip = true;
  }
  // TLine* line2 = new TLine(firstDipPos,0.,firstDipPos,histo->GetMaximum());
  // line2 -> Draw("same");
  
  */
  histo -> GetXaxis() -> SetRangeUser(5.,40.);
  
  int nPeaks = 3;                                                                                                                                                                                                                                                                                                                                                                                                          
  TSpectrum* spectrum = new TSpectrum(nPeaks);                                                                                                                                                                                                                                                                                                                                                                             
  int nFound = spectrum -> Search(histo, 0.5, "", 0.001);
  double* peaks = spectrum -> GetPositionX();
  
  
  //--- find and fit 511 keV peak
  TF1* fitFunc_511 = new TF1("fitFunc_511","gaus",0.90*peaks[0],1.10*peaks[0]);
  fitFunc_511 -> SetLineColor(kBlack);
  fitFunc_511 -> SetLineWidth(2);
  histo -> Fit(fitFunc_511,"QRS+");
  
  res["0.511 MeV"] = std::make_pair(fitFunc_511->GetParameter(1),fitFunc_511->GetParameter(2));
  histo -> GetXaxis() -> SetRangeUser(0.,3.*fitFunc_511->GetParameter(1));
  
  
  //--- find and fit 1275 keV peak
  int iPeak_1275;
  for(int iPeak = 1; iPeak < nFound; ++iPeak)
  {
    if( peaks[iPeak]/peaks[0] > 1.85 )
    {
      iPeak_1275 = iPeak;
      break;
    }
  }
  TF1* fitFunc_1275 = new TF1("fitFunc_1275","gaus",0.95*peaks[iPeak_1275],1.05*peaks[iPeak_1275]);
  fitFunc_1275 -> SetLineColor(kBlack);
  fitFunc_1275 -> SetLineWidth(2);
  histo -> Fit(fitFunc_1275,"QRS+");
  
  res["1.275 MeV"] = std::make_pair(fitFunc_1275->GetParameter(1),fitFunc_1275->GetParameter(2));
  
  
  if( ranges )
  {
    ranges->push_back(res["0.511 MeV"].first-8.*res["0.511 MeV"].second);
    ranges->push_back(res["0.511 MeV"].first-5.*res["0.511 MeV"].second);
    ranges->push_back(res["0.511 MeV"].first-2.*res["0.511 MeV"].second);
    ranges->push_back(res["0.511 MeV"].first+2.*res["0.511 MeV"].second);
    ranges->push_back(res["0.511 MeV"].first+5.*res["0.511 MeV"].second);
    ranges->push_back(0.5*(res["0.511 MeV"].first+5.*res["0.511 MeV"].second+res["1.275 MeV"].first-5.*res["1.275 MeV"].second));
    ranges->push_back(res["1.275 MeV"].first-5.*res["1.275 MeV"].second);
    ranges->push_back(res["1.275 MeV"].first-2.*res["1.275 MeV"].second);
    ranges->push_back(res["1.275 MeV"].first+2.*res["1.275 MeV"].second);
    
    for(auto range: (*ranges))
    {
      TLine* line = new TLine(range,3.,range,histo->GetBinContent(histo->FindBin(range)));
      line -> SetLineWidth(1);
      line -> SetLineStyle(7);
      line -> Draw("same");
    }
  }
  
  
  return res;
}
