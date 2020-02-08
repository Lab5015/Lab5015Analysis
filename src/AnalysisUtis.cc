#include "interface/AnalysisUtils.h"



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
  for(int bin = 1; bin <= histo->GetNbinsX(); ++bin)
  {
    if( histo->GetBinCenter(bin) < xMin ) continue;
    if( histo->GetBinCenter(bin) > xMax ) continue;
    if( histo->GetBinContent(bin) > max ) { max = histo->GetBinContent(bin); binMax = bin; };
  }
  return histo->GetBinCenter(binMax);
}
