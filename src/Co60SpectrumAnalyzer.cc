#include "interface/Co60SpectrumAnalyzer.h"



std::map<std::string,std::pair<float,float> > Co60SpectrumAnalyzer(TH1F* histo,
                                                                   std::vector<float>* ranges)
{
  gStyle -> SetOptFit(0);
  
  
  std::map<std::string,std::pair<float,float> > res;
  

  int startBin = 2;
  int endBin;
  
  while ( histo -> GetBinContent(startBin) < 100 && histo -> GetBinContent(startBin-1) < 100 ){
	startBin ++;
  }
 
  endBin = startBin;
  while ( histo -> GetBinContent(endBin) > 2 && histo -> GetBinContent(endBin-1) > 2 && histo -> GetBinContent(endBin+1) > 2){
		endBin++;
  }

  //Derivative calculation	
  std::vector<double> derivative;
  derivative.push_back(0);
  for( int bin = 1; bin < endBin; bin++){
  	if ( bin < startBin){
		derivative.push_back(0);
	
	}
	if( bin >= startBin){	
	
		double difference = histo->GetBinContent(bin+1) - histo->GetBinContent(bin);
		double distance =  histo->GetBinCenter(bin+1) - histo->GetBinCenter(bin);		
		derivative.push_back(difference/distance);		
	}
  }

  //Looking for the right side of the 1332 Peak
  std::vector <double> rangesDerivative;
  double drop = startBin+1;
  
  int counter = 0;
  double sum = 0;
  for ( int j = startBin; j < endBin; j++){
		sum += derivative[j];
		counter ++;
		
		if ( counter == 5 ){
			rangesDerivative.push_back(double(sum)/5.);
			counter = 0;
			sum = 0;	
		}
  }
  
  int check = 0;
  
  while ( check != 99){
	double min = *std::min_element( rangesDerivative.begin(), rangesDerivative.end());
	int minIndex = 5;
	do{minIndex ++;}
		while(rangesDerivative[minIndex] != min);
		
	check =  99;

	for ( int i = minIndex +1;  i < minIndex + 5; i++){
		if( rangesDerivative[i] > 0){
			check = 0;
		}
	}
	if (check == 0){
		rangesDerivative.erase(rangesDerivative.begin() + minIndex);
	}
	drop = histo -> GetBinCenter(startBin + (minIndex+5) *5);
	
 }
	

  //Searching 1173 peak
  int nPeaks = 4;
  TSpectrum* spectrum = new TSpectrum(nPeaks);

  histo -> GetXaxis() -> SetRangeUser(0., drop);
                                                                                                                                                                                                                                                                                                                                                                             
  int nFound = spectrum -> Search(histo, 0.01, "", 0.001);		
  double* peaks = spectrum -> GetPositionX();			
 
  std::vector<double> totalPeaks;
  std::vector<double> totalPeaks2;
  std::vector<double> realPeaks;


  if(nFound >1){
	for( int i = 0; i < nFound; i++){
		totalPeaks.push_back(peaks[i]);
 	}	
  
  	sort(totalPeaks.begin(), totalPeaks.end());
	
	for( int j = 0; j < totalPeaks.size(); j++){
	}	
  }

  if(totalPeaks.size() == 4){
  	realPeaks.push_back(totalPeaks[totalPeaks.size()-2]);
	realPeaks.push_back(totalPeaks[totalPeaks.size()-1]);
  }
  if(totalPeaks.size() <4){
	double d1 = drop - histo->GetBinWidth(10)*20 -totalPeaks[totalPeaks.size()-1];
	TF1* fitFunc_gaus = new TF1("fitFunc_gaus","gaus",0.95*totalPeaks[2],1.05*totalPeaks[2]);
	histo -> Fit(fitFunc_gaus,"QR"); 
	double sigma = fitFunc_gaus -> GetParameter(2);
	if(d1 > sigma*1.5){
		realPeaks.push_back(totalPeaks[totalPeaks.size()-1]);
	}
	if(d1 < sigma*1.5){
		realPeaks.push_back(totalPeaks[totalPeaks.size()-2]);
		realPeaks.push_back(totalPeaks[totalPeaks.size()-1]);
	}
  }

 

  //Fitting 1173 keV peak
  if (realPeaks.size() == 1 || realPeaks.size() == 2 ){
  	TF1* fitFunc_1173 = new TF1("fitFunc_1173","gaus",0.95*realPeaks[0],1.05*realPeaks[0]);
	fitFunc_1173 -> SetLineColor(kBlack);
	fitFunc_1173 -> SetLineWidth(2);
	histo -> Fit(fitFunc_1173,"QRS+"); 
	  
	res["1.173 MeV"] = std::make_pair(fitFunc_1173->GetParameter(1),fitFunc_1173->GetParameter(2));
	//histo -> GetXaxis() -> SetRangeUser(0.,3.*fitFunc_1173->GetParameter(1));			
  }

  //Searching 1332 peaks if previously failed
  if(realPeaks.size() == 1){
	  histo -> GetXaxis() -> SetRangeUser ( res["1.173 MeV"].first +res["1.173 MeV"].second, 40.);
	  nPeaks = 2;
	  TSpectrum* spectrum2 = new TSpectrum(nPeaks);
	  nFound = spectrum2 -> Search(histo, 0.01, "", 0.001);		
	  double* peaks2 = spectrum2 -> GetPositionX();			

	  if(nFound >0){
		for( int i = 0; i < nFound; i++){
			totalPeaks2.push_back(peaks2[i]);
	 	}	

	  	sort(totalPeaks2.begin(), totalPeaks2.end());	
	
	  	realPeaks.push_back(totalPeaks2[0]);
	  }
  }
  
  histo -> GetXaxis() -> SetRangeUser(0., 40.);
  
  //Fitting 1332 keV peak
  if (realPeaks.size() == 2){  

	TF1* fitFunc_1332 = new TF1("fitFunc_1332","gaus",0.95*realPeaks[1],1.05*realPeaks[1]);
	fitFunc_1332 -> SetLineColor(kBlack);
	fitFunc_1332 -> SetLineWidth(2);
	histo -> Fit(fitFunc_1332,"QRS+");		
	
	res["1.332 MeV"] = std::make_pair(fitFunc_1332->GetParameter(1),fitFunc_1332->GetParameter(2));
  
  }
 
 //Not Co60 spectrum controll
  if (realPeaks.size() == 0 || realPeaks.size() > 2 ){
	std::cout << "Errore" << std::endl;
	res["1.173 MeV"] = std::make_pair(-9999,0);
  }

  //Locating and drawing ranges of interest
  if ( realPeaks.size() > 0 && realPeaks.size() < 3 ){
  	ranges->push_back(res["1.173 MeV"].first-8.*res["1.173 MeV"].second);
  	ranges->push_back(res["1.173 MeV"].first-5.*res["1.173 MeV"].second);
	ranges->push_back(res["1.173 MeV"].first-2.*res["1.173 MeV"].second);
	ranges->push_back(res["1.173 MeV"].first+2.*res["1.173 MeV"].second);
	ranges->push_back(res["1.173 MeV"].first+5.*res["1.173 MeV"].second);
   }	
  if ( realPeaks.size() == 2 ){
	ranges->push_back(0.5*(res["1.173 MeV"].first+5.*res["1.173 MeV"].second+res["1.332 MeV"].first-5.*res["1.332 MeV"].second));
	ranges->push_back(res["1.332 MeV"].first-5.*res["1.332 MeV"].second);
	ranges->push_back(res["1.332 MeV"].first-2.*res["1.332 MeV"].second);
	ranges->push_back(res["1.332 MeV"].first+2.*res["1.332 MeV"].second);
   }
    
  
  for(auto range: (*ranges)){
	TLine* line = new TLine(range,3.,range,histo->GetBinContent(histo->FindBin(range)));
	line -> SetLineWidth(1);
	line -> SetLineStyle(7);
	line -> Draw("same");


    
  }
  
  return res;
}
