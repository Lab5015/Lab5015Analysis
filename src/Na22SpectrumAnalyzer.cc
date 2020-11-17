#include "interface/Na22SpectrumAnalyzer.h"



std::map<std::string,std::pair<float,float> > Na22SpectrumAnalyzer(TH1F* histo,
                                                                   std::vector<float>* ranges)
{
  gStyle -> SetOptFit(0);
  
  
  std::map<std::string,std::pair<float,float> > res;
  
  int nPeaks = 5;
  TSpectrum* spectrum = new TSpectrum(nPeaks);

  histo -> GetXaxis() -> SetRangeUser(0.,40.);
                                                                                                                                                                                                                                                                                                                                                                             
  int nFound = spectrum -> Search(histo, 0.5, "", 0.001);		
  double* peaks = spectrum -> GetPositionX();				
 
  //totalPeaks -> all found Peaks, nBins -> bin number of found Peaks, realPeaks -> 511 e 1275 
  std::vector<double> totalPeaks;
  std::vector<double> realPeaks;
  std::vector<int> nBins;
 
  if(nFound >2){
	
	for( int i = 0; i < nFound; i++){
  		nBins.push_back(histo -> FindBin(peaks[i]));
		totalPeaks.push_back(peaks[i]);
 	}
	sort(nBins.begin(), nBins.end());
  	sort(totalPeaks.begin(), totalPeaks.end());

	//Derivative calculation
	std::vector<double> derivative;
	derivative.push_back(0);

 	int unemptyBin = 2;
  	while ( histo -> GetBinContent(unemptyBin) < 100 && histo -> GetBinContent(unemptyBin-1) < 100 ){
		unemptyBin ++;
  	}
	for( int bin = 1; bin < nBins[nBins.size()-1] +20; bin++){
		if ( bin < unemptyBin){
			derivative.push_back(0);
		}
		if( bin >= unemptyBin){	
			double difference = histo->GetBinContent(bin+1) - histo->GetBinContent(bin);
			double distance =  histo->GetBinCenter(bin+1) - histo->GetBinCenter(bin);		
			derivative.push_back(difference/distance);		
		}
	  }



	//Finding  first and second maximum in totalPeaks
	//Fitting them with gaussian and getting their sigma
	double maxValue = 0;
	double secondMaxValue = 0;
  	int maxIndex = 0;
	int secondMaxIndex = 0;
	double maxSigma =0;  	
	double secondMaxSigma =0;

	//First max
  	for (int i = 0; i < nFound; i++){
		if ( histo -> GetBinContent(nBins[i]) > maxValue){
			maxValue = histo -> GetBinContent(nBins[i]);
			maxIndex = i;
		}
	}
	TF1* fitFunc_gaus = new TF1("fitFunc_gaus","gaus",0.95*double(totalPeaks[maxIndex]),1.05*double(totalPeaks[maxIndex]));
	histo -> Fit(fitFunc_gaus,"QR");
	maxSigma = fitFunc_gaus->GetParameter(2);
	
	//Second Max
	for (int i = 0; i < nFound; i++){
		if(i != maxIndex){
			if ( histo -> GetBinContent(nBins[i]) > secondMaxValue){
				secondMaxValue = histo -> GetBinContent(nBins[i]);
				secondMaxIndex = i;
			}
		}
		
	}
	TF1* fitFunc_gaus2 = new TF1("fitFunc_gaus2","gaus",0.95*double(totalPeaks[secondMaxIndex]),1.05*double(totalPeaks[secondMaxIndex]));
	histo -> Fit(fitFunc_gaus2,"QR");
	secondMaxSigma = fitFunc_gaus2->GetParameter(2);


	//Finding compton edge, looking for the range (width 10 bins) with smallest derivative not included in [first max - second max]
	int sigmaBin = std::min(maxSigma, secondMaxSigma) / (histo->GetBinWidth(10));
	int windowBin = 10;	
	double minPendence = 1000;
	int compton;

	for ( int j = 1; j < nFound-1; j++){

		int nRanges = 0;
		if ( double(nBins[j+1] -sigmaBin * 2 - nBins[j])/windowBin - int((nBins[j+1] -sigmaBin * 2 - nBins[j])/windowBin) < 0.5){
			nRanges = int((nBins[j+1] -sigmaBin * 2 - nBins[j])/windowBin);
		}
		if ( double(nBins[j+1] -sigmaBin * 2 - nBins[j])/windowBin - int((nBins[j+1] -sigmaBin * 2. - nBins[j])/windowBin) >= 0.5){
			nRanges = int((nBins[j+1] -sigmaBin * 2 - nBins[j])/windowBin) +1;
		}
		
		for(int r = 0; r < nRanges; r++){
			double pendence = 0;	
			for ( int i = nBins[j]+ sigmaBin + windowBin*r; i < nBins[j] + sigmaBin +windowBin*(r+1); i++){
				pendence += derivative[i];
			}
			
			if(std::abs(pendence/windowBin) <= minPendence && std::abs(pendence/windowBin)>0.1 && j+1 > maxIndex && j+1 > secondMaxIndex){
				minPendence = std::abs(pendence/windowBin);
				compton = j+1;
			}
		}
  	}
	
			
   	
   	//Compton is the i-th found peaks after the compton edge
	//Based on the number of peaks before and after the compton edge, 511 and 1275 Peaks are selected
	if( compton > maxIndex && compton > secondMaxIndex){
		int firstPeak = 0;
		int afterCompton = nFound -1 -compton;
		
		do{ firstPeak++;}
				while(double(histo->GetBinContent(nBins[compton-firstPeak]))/double(maxValue) < 0.1);

		realPeaks.push_back(totalPeaks[compton-firstPeak]);
		if(afterCompton == 1 || afterCompton == 0){
			realPeaks.push_back(totalPeaks[compton+afterCompton]);
		}
		
		if(afterCompton > 1){
			double peaksRatio = histo->GetBinContent(nBins[nFound-2]) / histo->GetBinContent(nBins[nFound-1]);
			if ( peaksRatio > 3 ) { 
				realPeaks.push_back(totalPeaks[nFound-2]);
			}
			if ( peaksRatio < 3 && peaksRatio > 0.3333){				
				realPeaks.push_back(totalPeaks[nFound-1]);	 
			}
		}
	}
  }

  	
	
  //Fitting 511 keV peak
  if (realPeaks.size() == 1 || realPeaks.size() == 2 ){
  	TF1* fitFunc_511 = new TF1("fitFunc_511","gaus",0.90*realPeaks[0],1.10*realPeaks[0]);
	fitFunc_511 -> SetLineColor(kBlack);
	fitFunc_511 -> SetLineWidth(2);
	histo -> Fit(fitFunc_511,"QRS+");	
	res["0.511 MeV"] = std::make_pair(fitFunc_511->GetParameter(1),fitFunc_511->GetParameter(2));
	histo -> GetXaxis() -> SetRangeUser(0.,3.*fitFunc_511->GetParameter(1));			
  }
  
  //Fitting 1275 keV peak
  if (realPeaks.size() == 2){  
	TF1* fitFunc_1275 = new TF1("fitFunc_1275","gaus",0.95*realPeaks[1],1.05*realPeaks[1]);
	fitFunc_1275 -> SetLineColor(kBlack);
	fitFunc_1275 -> SetLineWidth(2);
	histo -> Fit(fitFunc_1275,"QRS+");		
	res["1.275 MeV"] = std::make_pair(fitFunc_1275->GetParameter(1),fitFunc_1275->GetParameter(2));
  }
 
  //Not Na22 spectrum controll
  if (realPeaks.size() == 0 || realPeaks.size() > 2 || nFound < 3){
	std::cout << "Errore" << std::endl;
	res["0.511 MeV"] = std::make_pair(-9999,0);		
  }

  //Locating and drawing ranges of interest
  if ( realPeaks.size() > 0 && realPeaks.size() < 3 ){  
  	ranges->push_back(res["0.511 MeV"].first-8.*res["0.511 MeV"].second);
  	ranges->push_back(res["0.511 MeV"].first-5.*res["0.511 MeV"].second);
	ranges->push_back(res["0.511 MeV"].first-2.*res["0.511 MeV"].second);
	ranges->push_back(res["0.511 MeV"].first+2.*res["0.511 MeV"].second);
	ranges->push_back(res["0.511 MeV"].first+5.*res["0.511 MeV"].second);
   }	

  if ( realPeaks.size() == 2 ){
	ranges->push_back(0.5*(res["0.511 MeV"].first+5.*res["0.511 MeV"].second+res["1.275 MeV"].first-5.*res["1.275 MeV"].second));
	ranges->push_back(res["1.275 MeV"].first-5.*res["1.275 MeV"].second);
	ranges->push_back(res["1.275 MeV"].first-2.*res["1.275 MeV"].second);
	ranges->push_back(res["1.275 MeV"].first+2.*res["1.275 MeV"].second);
   }
    
  for(auto range: (*ranges)){
	TLine* line = new TLine(range,3.,range,histo->GetBinContent(histo->FindBin(range)));
	line -> SetLineWidth(1);
	line -> SetLineStyle(7);
	line -> Draw("same");
    
  }
  
  
  return res;
}

