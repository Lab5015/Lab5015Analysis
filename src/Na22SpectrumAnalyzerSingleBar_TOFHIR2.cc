#include "interface/Na22SpectrumAnalyzerSingleBar_TOFHIR2.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"


std::map<std::string,std::pair<float,float> > Na22SpectrumAnalyzerSingleBar_TOFHIR2(TH1F* histo,
                                                                   std::vector<float>* ranges)
{
  gStyle -> SetOptFit(0);
  
  std::map<std::string,std::pair<float,float> > res;
  
  
  // use TSpectrum to localize peaks
  histo -> GetXaxis() -> SetRangeUser(1.,500.);
  
  int nPeaks = 4;
  TSpectrum* spectrum = new TSpectrum(nPeaks);
  int nFound = spectrum -> Search(histo, 20., "", 0.001);		
  double* peaks = spectrum -> GetPositionX();				
 
  if(histo->GetNbinsX() == 2048){
	histo->RebinX(2);
  }
	
  // totalPeaks -> all found Peaks, nBins -> bin number of found Peaks, realPeaks -> 511 and 1275 
  std::vector<double> totalPeaks;
  std::vector<double> realPeaks;
  std::vector<int> nBins;
  double maxValue = 0;
  double secondMaxValue = 0;
  int maxIndex = 0;
  int secondMaxIndex = 0;
  double maxSigma = 0;  	
  double secondMaxSigma = 0;
  int endBin = 0;
  int nFoundAdd = 0;

  if(nFound > 1){
  	for(int i = 0; i < nFound; ++i)
		{
  		nBins.push_back( histo->FindBin(peaks[i]) );
  		realPeaks.push_back( peaks[i] );
				//std::cout << peaks[i] << "    aaa" << std::endl;
		}
		sort(nBins.begin(),nBins.end());
		sort(realPeaks.begin(),realPeaks.end());
	}

		/*endBin = nBins[0]; 
		while ( histo -> GetBinContent(endBin) > 10 && histo -> GetBinContent(endBin-1) > 10 && histo -> GetBinContent(endBin+1) > 10 ){
			endBin++;
		}
	if (fabs(histo -> GetBinCenter(histo -> GetMaximumBin()) - totalPeaks[0]) > 10 && fabs(histo -> GetBinCenter(histo -> GetMaximumBin()) - totalPeaks[1]) > 10){ 
		histo -> GetXaxis() -> SetRangeUser(0, totalPeaks[0]-10);
 		int nPeaks3 = 2;
  		TSpectrum* spectrum3 = new TSpectrum(nPeaks3);
  		int nFound3 = spectrum3 -> Search(histo, 0.1, "", 0.001);		
  		double* peaks3 = spectrum3 -> GetPositionX();	
		if(nFound3 ==2 ){
			totalPeaks.push_back(peaks3[1]);
		nBins.push_back( histo->FindBin(peaks3[1]) );
		}
		if(nFound3 != 0){
		totalPeaks.push_back(peaks3[0]);
		nBins.push_back( histo->FindBin(peaks3[0]) );
		sort(nBins.begin(),nBins.end());
		sort(totalPeaks.begin(),totalPeaks.end());
		}
	}	
	
	if (histo->GetBinContent(nBins[0])/histo->GetBinContent(nBins[1]) < histo->GetBinContent(nBins[1])/histo->GetBinContent(nBins[2])){
		realPeaks.push_back(totalPeaks[1]);
	}
	if (histo->GetBinContent(nBins[0])/histo->GetBinContent(nBins[1]) > histo->GetBinContent(nBins[1])/histo->GetBinContent(nBins[2])){
		realPeaks.push_back(totalPeaks[0]);
	}*/
			 
    //Fitting 511 keV peak
  	if( realPeaks.size() >= 1 )
  	{
			realPeaks[0] = histo->GetBinCenter(histo->GetMaximumBin());
  	  TF1*fitFunc_511 = new TF1("fitFunc_511","gaus",realPeaks[0] *0.90,realPeaks[0] *1.2);
  	  fitFunc_511 -> SetLineColor(kBlack);
  	  fitFunc_511 -> SetLineWidth(2);
  	  histo -> Fit(fitFunc_511,"QRS+");	
  	  res["0.511 MeV"] = std::make_pair(fitFunc_511->GetParameter(1),fitFunc_511->GetParameter(2));
			
  	  histo -> GetXaxis() -> SetRangeUser(0.,3.*fitFunc_511->GetParameter(1));			
  	}
  /* float med1 = 0;
	 float med2 = 0;
	 int co = 0;
	 for( int i = nBins[totalPeaks.size() -2]; i < nBins[totalPeaks.size() -1]; i++){
		med1 += histo->GetBinContent(i);
		co ++;
	 }
	 med1 /= float(co);
	 co = 0;
	 for( int i = nBins[totalPeaks.size() -3]; i < nBins[totalPeaks.size() -2]; i++){
		med2 += histo->GetBinContent(i);
		co ++;
	 }
	 med2 /= float(co);
	
  
  
  	if( realPeaks.size() ==1 ){
  		histo -> GetXaxis() -> SetRangeUser(totalPeaks[totalPeaks.size() -2] +20 ,500.);
 		  int nPeaks2 = 2;
  		TSpectrum* spectrum2 = new TSpectrum(nPeaks2);
  		int nFound2 = spectrum2 -> Search(histo, 0.5, "", 0.001);		
  		double* peaks2 = spectrum2 -> GetPositionX();	
			if (nFound2 != 0 & peaks2[0] < peaks2[1]){
				realPeaks.push_back(peaks2[0]);
			}
			if (nFound2 != 0 & peaks2[0] > peaks2[1]){
				realPeaks.push_back(peaks2[1]);
			}
			if (nFound2 > 1 & peaks2[0] > peaks2[1]){
				realPeaks.push_back(peaks2[0]);
			}
			if (nFound2 > 1 & peaks2[0] < peaks2[1]){
				realPeaks.push_back(peaks2[1]);
			}
 		}*/
	/*if( realPeaks.size() ==2 ){
  		histo -> GetXaxis() -> SetRangeUser(realPeaks[1] +20 ,500.);
 		int nPeaks4 = 1;
  		TSpectrum* spectrum4 = new TSpectrum(nPeaks4);
  		int nFound4 = spectrum4 -> Search(histo, 0.5, "", 0.001);		
  		double* peaks4 = spectrum4 -> GetPositionX();	
		if (nFound4 != 0){
			realPeaks.push_back(peaks4[0]);
		}
 	}*/
  
  //}

  histo -> GetXaxis() -> SetRangeUser(0.  ,500.);
  
  //Fitting 1275 keV peak
  if( realPeaks.size() > 1 )
  {  
    TF1* fitFunc_1275 = new TF1("fitFunc_1275","gaus",realPeaks[1]-5, realPeaks[1]+5);
    fitFunc_1275 -> SetParameter(1, realPeaks[1]);
    fitFunc_1275 -> SetParameter(2, res["0.511 Mev"].second/2.);
    histo -> Fit(fitFunc_1275,"QRS+");
    
    fitFunc_1275 -> SetParameter(0, fitFunc_1275 -> GetParameter(0));
    fitFunc_1275 -> SetParameter(1, fitFunc_1275 -> GetParameter(1));
    fitFunc_1275 -> SetParameter(2, fitFunc_1275 -> GetParameter(2));    
    fitFunc_1275 -> SetLineColor(kBlack);
    fitFunc_1275 -> SetLineWidth(2);
    histo -> Fit(fitFunc_1275,"QRS+", "", fitFunc_1275->GetParameter(1)-fitFunc_1275->GetParameter(2),fitFunc_1275->GetParameter(1)+1.5*fitFunc_1275->GetParameter(2));		
    res["1.275 MeV"] = std::make_pair(fitFunc_1275->GetParameter(1),fitFunc_1275->GetParameter(2));
  }

   //Fitting 1786 keV peak
  /*int c = 0;
  if(realPeaks.size()==3){
  	TF1* fitFunc_1786 = new TF1("fitFunc_1786","gaus",realPeaks[2]-2, realPeaks[2]+2);
		//fitFunc_1786 -> SetLineColor(kBlack);
    //fitFunc_1786 -> SetLineWidth(2);
	  fitFunc_1786 -> SetParameter(1, realPeaks[2]);
		fitFunc_1786 -> SetParameter(2, res["1275 Mev"].second/2.);
    histo -> Fit(fitFunc_1786,"QRS+");
    
    fitFunc_1786 -> SetParameter(0, fitFunc_1786 -> GetParameter(0));
    fitFunc_1786 -> SetParameter(1, fitFunc_1786 -> GetParameter(1));
    fitFunc_1786 -> SetParameter(2, fitFunc_1786 -> GetParameter(2));
    fitFunc_1786 -> SetLineColor(kBlack);
    fitFunc_1786 -> SetLineWidth(2);
    histo -> Fit(fitFunc_1786,"QRS+", "", fitFunc_1786->GetParameter(1)-fitFunc_1786->GetParameter(2),fitFunc_1786->GetParameter(1)+fitFunc_1786->GetParameter(2));
		histo -> Fit(fitFunc_1786,"QRS+", "", fitFunc_1786->GetParameter(1)-fitFunc_1786->GetParameter(2),fitFunc_1786->GetParameter(1)+fitFunc_1786->GetParameter(2));		
    res["1.786 MeV"] = std::make_pair(fitFunc_1786->GetParameter(1),fitFunc_1786->GetParameter(2));
		c = 1;
  } */

	



  
  //Not Na22 spectrum controll
  if( realPeaks.size() == 0 || realPeaks.size() > 3  )
  {
    std::cout << "Errore" << std::endl;
    res["0.511 MeV"] = std::make_pair(-9999,0);		
  }
  
	int ultimo_bin=0;
  for(int i=histo->GetNbinsX(); i>0; i--){
   if(histo->GetBinContent(i)>2){ 
     ultimo_bin=i;
		 break;
	  }
  }





//single Bar
					
if ( realPeaks.size() > 0 && realPeaks.size() < 4 ){
  	  //ranges->push_back(0.5*(histo->GetBinCenter(unemptyBin)+res["0.511 MeV"].first-0.6*res["0.511 MeV"].first));
	  //ranges->push_back(res["0.511 MeV"].first-0.6*res["0.511 MeV"].first);
  	  //ranges->push_back(res["0.511 MeV"].first-0.4*res["0.511 MeV"].first);
	  //ranges->push_back((res["0.511 MeV"].first-0.4*res["0.511 MeV"].first+res["0.511 MeV"].first-0.1*res["0.511 MeV"].first)/2.);
	  ranges->push_back(res["0.511 MeV"].first-0.1*res["0.511 MeV"].first);
	  ranges->push_back(res["0.511 MeV"].first+0.1*res["0.511 MeV"].first);
	  //ranges->push_back(0.5*(res["0.511 MeV"].first+0.1*res["0.511 MeV"].first+(2*0.25*(res["0.511 MeV"].first+0.1*res["0.511 MeV"].first+res["1.275 MeV"].first-0.1*res["1.275 MeV"].first))));
		//ranges->push_back(res["0.511 MeV"].first+0.4*res["0.511 MeV"].first);	
	
}	


  if ( realPeaks.size() == 2  || realPeaks.size() == 3){
		//ranges->push_back(0.25*(res["0.511 MeV"].first+0.1*res["0.511 MeV"].first+res["1.275 MeV"].first-0.1*res["1.275 MeV"].first));
		//ranges->push_back(2*0.25*(res["0.511 MeV"].first+0.1*res["0.511 MeV"].first+res["1.275 MeV"].first-0.1*res["1.275 MeV"].first));
		
		//ranges->push_back(0.5*(res["1.275 MeV"].first-1*res["1.275 MeV"].second+2*0.25*(res["0.511 MeV"].first+0.1*res["0.511 MeV"].first+res["1.275 MeV"].first-0.1*res["1.275 MeV"].first)));
		ranges->push_back(res["1.275 MeV"].first-1*res["1.275 MeV"].second);
		ranges->push_back(res["1.275 MeV"].first+1*res["1.275 MeV"].second);
		//ranges->push_back(res["1.275 MeV"].first+2*res["1.275 MeV"].second);
		//ranges->push_back(0.5*(res["1.275 MeV"].first+0.01*res["1.275 MeV"].first+res["1.786 MeV"].first-0.06*res["1.786 MeV"].first));
		//ranges->push_back(res["1.786 MeV"].first-0.06*res["1.786 MeV"].first);
		//ranges->push_back(res["1.786 MeV"].first+0.06*res["1.786 MeV"].first);
		
		//ranges->push_back((histo->GetBinCenter(ultimo_bin)-(res["1.275 MeV"].first+0.05*res["1.275 MeV"].first))/4+res["1.275 MeV"].first+0.05*res["1.275 MeV"].first);
		//ranges->push_back(histo->GetBinCenter(ultimo_bin));
   }




  if ( realPeaks.size() == 1 ){
		float second = res["0.511 MeV"].first *2.5;
		//ranges->push_back(0.5*(res["0.511 MeV"].first+0.4*res["0.511 MeV"].second-0.2*second));
		//ranges->push_back(second-0.2*second);
		ranges->push_back(second-0.1*second);
		ranges->push_back(second+0.1*second);
		//ranges->push_back((histo->GetBinCenter(ultimo_bin)-(res["1.275 MeV"].first+0.1*res["1.275 MeV"].first))/3+res["1.275 MeV"].first+0.1*res["1.275 MeV"].first);
		//ranges->push_back((res["1.275 MeV"].first)*1.25+0.2*(res["1.275 MeV"].first)*1.25);
		//ranges->push_back(histo->GetBinCenter(ultimo_bin));
   }
  
  for(auto range: (*ranges)){
		
		if ( range  < 0  && realPeaks.size() == 1) { //error controll
			std::cout << "Errore" << std::endl;
    	res["0.511 MeV"] = std::make_pair(-9999,0);	
		}
	  float yval = std::max(10., histo->GetBinContent(histo->FindBin(range)));
		TLine* line = new TLine(range,0.,range, yval);
		line -> SetLineWidth(1);
		line -> SetLineStyle(7);
		line -> Draw("same");
    
  }
  
  
  return res;
}

