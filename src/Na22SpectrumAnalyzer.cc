#include "interface/Na22SpectrumAnalyzer.h"



std::map<std::string,std::pair<float,float> > Na22SpectrumAnalyzer(TH1F* histo,
                                                                   std::vector<float>* ranges)
{
  gStyle -> SetOptFit(0);
  
  std::map<std::string,std::pair<float,float> > res;
  
  
  // use TSpectrum to localize peaks
  histo -> GetXaxis() -> SetRangeUser(1.,40.);
  
  int nPeaks = 5;
  TSpectrum* spectrum = new TSpectrum(nPeaks);
  int nFound = spectrum -> Search(histo, 0.5, "", 0.001);		
  double* peaks = spectrum -> GetPositionX();				
  
  
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
		    totalPeaks.push_back( peaks[i] );
		  }
		  sort(nBins.begin(),nBins.end());
		  sort(totalPeaks.begin(),totalPeaks.end());

		endBin = nBins[0]; 
		while ( histo -> GetBinContent(endBin) > 10 && histo -> GetBinContent(endBin-1) > 10 && histo -> GetBinContent(endBin+1) > 10 ){
			endBin++;
		}
	//looking for one other peak after the found ones
		if(totalPeaks[totalPeaks.size()-1]+2 < histo->GetBinCenter(endBin)){
			histo -> GetXaxis() -> SetRangeUser(totalPeaks[totalPeaks.size()-1]*1.1 ,40.);
	 		int nPeaksAdd = 1;
			TSpectrum* spectrumAdd = new TSpectrum(nPeaksAdd);
			nFoundAdd = spectrumAdd -> Search(histo, 0.5, "", 0.001);		
			double* peaksAdd = spectrumAdd -> GetPositionX();	
			if (nFoundAdd != 0 && peaksAdd[0] != totalPeaks[totalPeaks.size()-1] && histo->FindBin(peaksAdd[0]) < endBin ){
				totalPeaks.push_back(peaksAdd[0]);
				nBins.push_back(histo->FindBin(peaksAdd[0]));
			}
		}
	}
  
  histo -> GetXaxis() -> SetRangeUser(0.  ,40.);

 
		
 
  if( nFound + nFoundAdd > 2 )
  {
    
    
    //Derivative calculation
    std::vector<double> derivative;
    derivative.push_back(0);
    
    int unemptyBin = 2;
    while( histo -> GetBinContent(unemptyBin) < 100 && histo -> GetBinContent(unemptyBin-1) < 100 )
    {
      ++unemptyBin;
    }
    for(int bin = 1; bin < nBins[nBins.size()-1]+20; ++bin)
    {
      if ( bin < unemptyBin) derivative.push_back(0);
      else
      {	
        double difference = histo->GetBinContent(bin+1) - histo->GetBinContent(bin);
        double distance =  histo->GetBinCenter(bin+1) - histo->GetBinCenter(bin);		
        derivative.push_back(difference/distance);		
      }
    }
    
    // Finding  first and second maximum in totalPeaks
    // Fitting them with gaussian and getting their sigma
    
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
		if(nFound < 4){
			secondMaxValue = maxValue;
			secondMaxIndex = maxIndex;
		}
		if(nFound > 3){
		  for (int i = 0; i < nFound; i++){
		    if(i != maxIndex){
		      if ( histo -> GetBinContent(nBins[i]) > secondMaxValue){
		        secondMaxValue = histo -> GetBinContent(nBins[i]);
		        secondMaxIndex = i;
		      }
		    }
		  }
		}
    TF1* fitFunc_gaus2 = new TF1("fitFunc_gaus2","gaus",0.95*double(totalPeaks[secondMaxIndex]),1.05*double(totalPeaks[secondMaxIndex]));
    histo -> Fit(fitFunc_gaus2,"QR");
    secondMaxSigma = fitFunc_gaus2->GetParameter(2);
    
		 //Finding compton edge between 511 e 1275
		 int sigmaBin = std::min(maxSigma, secondMaxSigma) / (histo->GetBinWidth(10));
		 int windowBin = 5;	
		 double minPendence = 1000;
		 int compton = totalPeaks.size()-1;
		 
		 for ( int j = 0; j < nFound-1; j++)
		 {
		   int nRanges = 0;
		   if ( double(nBins[j+1] -sigmaBin * 2 - nBins[j])/windowBin - int((nBins[j+1] -sigmaBin * 2 - nBins[j])/windowBin) < 0.5){
		     nRanges = int((nBins[j+1] -sigmaBin * 2 - nBins[j])/windowBin);
		   }
		   if ( double(nBins[j+1] -sigmaBin * 2 - nBins[j])/windowBin - int((nBins[j+1] -sigmaBin * 2. - nBins[j])/windowBin) >= 0.5){
		   	nRanges = int((nBins[j+1] -sigmaBin * 2 - nBins[j])/windowBin) +1;
		   }
		   
		   for(int r = 0; r < nRanges; r++){
		    double pendence = 0;
				int counter = 0;	
		    for ( int i = nBins[j]+ sigmaBin + windowBin*r; i < nBins[j] + sigmaBin +windowBin*(r+1); i++){
		    	pendence += derivative[i];
			  	counter ++;
		    }
		
		    if(std::abs(pendence/counter) <= minPendence &&  j+1 > maxIndex && j+1 > secondMaxIndex && histo->GetBinContent(nBins[j+1]) > 50 ){ //std::abs(pendence/counter)>0.1 &&
		      minPendence = std::abs(pendence/counter);
		      compton = j+1;
		    }
			}
		}
	  int afterCompton = nFound + nFoundAdd - 1 - compton;
	  // Compton is the i-th found peaks after the compton edge
	  // Looking for 511 peak valuating the amplitude of peaks before compton
		  
	  int firstPeak = 0;
		int found1275  = 0; 
		while ( found1275 == 0){ 
	  	firstPeak++; 
	  	if ( double(histo->GetBinContent(nBins[compton-firstPeak]))/double(secondMaxValue) > 0.8 || double(histo->GetBinContent(nBins[compton-firstPeak]))/double(maxValue) > 0.2  ){ 		
				realPeaks.push_back(totalPeaks[compton-firstPeak]);
				found1275 = 1;
			}
		}
		 
    //Fitting 511 keV peak
  	if( realPeaks.size() >= 1 )
  	{
  	  TF1*fitFunc_511 = new TF1("fitFunc_511","gaus",realPeaks[0] *0.90,realPeaks[0] *1.1);
  	  fitFunc_511 -> SetLineColor(kBlack);
  	  fitFunc_511 -> SetLineWidth(2);
  	  histo -> Fit(fitFunc_511,"QRS+");	
  	  res["0.511 MeV"] = std::make_pair(fitFunc_511->GetParameter(1),fitFunc_511->GetParameter(2));
  	  histo -> GetXaxis() -> SetRangeUser(0.,3.*fitFunc_511->GetParameter(1));			
  	}
   
		//Looking for 1275 keV peak in a window [1.8* 511 position, 4* 511 position]
		int secondFound = 0;
		std::vector <int> secondPeaksIndex;
				
		for( int i = 0 ; i < totalPeaks.size(); i++){
			if ( totalPeaks[i] > res["0.511 MeV"].first * 1.8 && totalPeaks[i] < res["0.511 MeV"].first * 4){
				secondFound ++;
				secondPeaksIndex.push_back(i);
			}
		}
		if( secondFound == 1){
			realPeaks.push_back(totalPeaks[secondPeaksIndex[0]]);
		}
   
		if( secondPeaksIndex.size() == 2){
			double peaksRatio = histo->GetBinContent(nBins[secondPeaksIndex[0]]) / histo->GetBinContent(nBins[secondPeaksIndex[1]]);
    	if ( peaksRatio > 3 )                       realPeaks.push_back(totalPeaks[secondPeaksIndex[0]]);
    	if ( peaksRatio < 3 && peaksRatio > 0.3333) realPeaks.push_back(totalPeaks[secondPeaksIndex[1]]);
		}
			
		if( secondPeaksIndex.size() == 3){
			realPeaks.push_back(totalPeaks[secondPeaksIndex[1]]);
		}
		if( secondPeaksIndex.size() == 4){
			realPeaks.push_back(totalPeaks[secondPeaksIndex[2]]);
		}
  }
  
  if( realPeaks.size() ==1 ){
		histo -> GetXaxis() -> SetRangeUser(totalPeaks[totalPeaks.size()-1]+2 ,40.);
 		int nPeaks2 = 1;
  	TSpectrum* spectrum2 = new TSpectrum(nPeaks2);
  	int nFound2 = spectrum2 -> Search(histo, 0.5, "", 0.001);		
  	double* peaks2 = spectrum2 -> GetPositionX();	
		if (nFound2 != 0){
			realPeaks.push_back(peaks2[0]);
		}
	}
  
  histo -> GetXaxis() -> SetRangeUser(0.  ,40.);
  
  //Fitting 1275 keV peak
  if( realPeaks.size() == 2 )
  {  
    TF1* fitFunc_1275 = new TF1("fitFunc_1275","gaus",realPeaks[1] *0.90, realPeaks[1] *1.10);
		fitFunc_1275 -> SetParameter(1, realPeaks[1]);
		fitFunc_1275 -> SetParameter(2, res["0.511 Mev"].second/2.);
    fitFunc_1275 -> SetLineColor(kBlack);
    fitFunc_1275 -> SetLineWidth(2);
    histo -> Fit(fitFunc_1275,"QRS+");		
    res["1.275 MeV"] = std::make_pair(fitFunc_1275->GetParameter(1),fitFunc_1275->GetParameter(2));
  }
  
  //Not Na22 spectrum controll
  if( realPeaks.size() == 0 || realPeaks.size() > 2  )
  {
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
		
		if ( range  < 0  && realPeaks.size() == 1) { //error controll
			std::cout << "Errore" << std::endl;
    	res["0.511 MeV"] = std::make_pair(-9999,0);	
		}
		TLine* line = new TLine(range,3.,range,histo->GetBinContent(histo->FindBin(range)));
		line -> SetLineWidth(1);
		line -> SetLineStyle(7);
		line -> Draw("same");
    
  }
  
  
  return res;
}

