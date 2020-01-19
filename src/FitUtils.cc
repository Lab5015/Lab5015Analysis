#include "interface/FitUtils.h"



/*** find effective sigma ***/
void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction)
{
  float integralMax = fraction * histo->Integral();
  
  // find first and last bin with non-null content
  int N = histo -> GetNbinsX();
  int M1 = 1;
  int M2 = 1;
  for(int bin1 = 1; bin1 <= N; ++bin1)
  {
    if( histo->GetBinContent(bin1) > 0. && M1 == 1 ) M1 = bin1;
    if( histo->GetBinContent(bin1) > 0. )            M2 = bin1;
  }

  std::map<int,float> binCenters;
  std::map<int,float> binContents;
  std::map<int,float> binIntegrals;
  for(int bin1 = M1; bin1 <= M2; ++bin1)
  {
    if( histo->GetBinContent(bin1) == 0 ) continue;

    binCenters[bin1] = histo->GetBinCenter(bin1);
    binContents[bin1] = histo->GetBinContent(bin1);

    for(int bin2 = M1; bin2 <= bin1; ++bin2)
      binIntegrals[bin1] += histo->GetBinContent(bin2);
  }
  
  float min = 0.;
  float max = 0.;
  float delta = 999999.;
  for(std::map<int,float>::const_iterator mapIt1 = binIntegrals.begin(); mapIt1 != binIntegrals.end(); ++mapIt1)
  {
    for(std::map<int,float>::const_iterator mapIt2 = ++binIntegrals.begin(); mapIt2 != binIntegrals.end(); ++mapIt2)
    {
      if( (mapIt2->second-mapIt1->second) < integralMax ) continue;

      float tmpMin = binCenters[mapIt1->first];
      float tmpMax = binCenters[mapIt2->first];

      if( (tmpMax-tmpMin) < delta )
      {
        delta = (tmpMax - tmpMin);
        min = tmpMin;
        max = tmpMax;
      }

      break;
    }
  }
  
  TH1F* smallHisto = (TH1F*)( histo->Clone("smallHisto") );
  for(int bin = 1; bin <= smallHisto->GetNbinsX(); ++bin)
  {
    if( smallHisto->GetBinCenter(bin) < min )
      smallHisto -> SetBinContent(bin,0);

    if( smallHisto->GetBinCenter(bin) > max )
      smallHisto -> SetBinContent(bin,0);
  }
  smallHisto -> SetFillColor(kYellow);

  float mean = smallHisto -> GetMean();
  float meanErr = smallHisto -> GetMeanError();
  float RMS = smallHisto -> GetRMS();
  float RMSErr = smallHisto -> GetRMSError();

  ret[0] = mean;
  ret[1] = meanErr;
  ret[2] = RMS;
  ret[3] = RMSErr;
  ret[4] = min;
  ret[5] = max;

  // std::cout << "mean: " << mean << "   min: " << min << "   max: " << max << std::endl;
}
