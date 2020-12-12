#include "interface/Co60SpectrumAnalyzer.h"



std::map<std::string,std::pair<float,float> > Co60SpectrumAnalyzer(TH1F* histo,
                                                                   std::vector<float>* ranges)
{
  gStyle -> SetOptFit(0);
  histo->GetXaxis()->SetRangeUser(0,50);
  std::map<std::string,std::pair<float,float> > res;
  

  int startBin = 2;
  int endBin;
  
  while ( histo -> GetBinContent(startBin) < 100 && histo -> GetBinContent(startBin-1) < 100 ){
	startBin ++;
  }

  endBin =  histo->GetNbinsX() /2.;
  while ( histo -> GetBinContent(endBin) < 50 && histo -> GetBinContent(endBin-1) < 50 && histo -> GetBinContent(endBin+1) < 50 && endBin < histo->GetNbinsX() -50){
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
  std::vector <int> erasedIndex;
	erasedIndex.push_back(0);
	erasedIndex.push_back(rangesDerivative.size());
	int minIndex = 5;
	double min;
  while ( check != 99){
		min = 10000;
		for( int j = 0; j < erasedIndex.size()-1; j++){
			for ( int i = erasedIndex[j] +1; i < erasedIndex[j+1] ; i++){
				if(rangesDerivative[i] < min){
					min = rangesDerivative[i];
					minIndex = i;
				}
			}
		}
		check =  99;
		int risalite = 0;
		int esclusioni = 0;
		for ( int i = minIndex -1;  i < minIndex + 4; i++){
			if( rangesDerivative[i] > 0){
				risalite ++;
				if( rangesDerivative[i] > 160){
					esclusioni ++;
				}
			}
		}
		if( risalite > 1 || esclusioni > 0){
			check = 0;
		}
		if (check == 0){
			erasedIndex.push_back(minIndex);
		}
		sort(erasedIndex.begin(), erasedIndex.end());	
		drop = histo -> GetBinCenter(startBin + (minIndex+5) *5);
	}
	std::vector<double> totalPeaks;
	std::vector<double> realPeaks;
	
	TF1* fitFunc;
	

	if (drop > histo->GetBinCenter((endBin-startBin)/3.) ){
		//Searching 1173 and 1332 peaks
		int nPeaks = 4;
		int nFound = 0;
		TSpectrum* spectrum = new TSpectrum(nPeaks);
		float start =  histo->GetBinCenter(startBin);
		float stop =   histo->GetBinCenter(endBin);
		histo -> GetXaxis() -> SetRangeUser(start, stop);
	                                                                                                                                                                                                                                                                                                                                                                           
		nFound = spectrum -> Search(histo, 1,"" , 0.001);
		double* peaks = spectrum -> GetPositionX();
		int index = -1;
		if (nFound > 1){
			for( int i = 0; i < nFound; i++){
				totalPeaks.push_back(peaks[i]);
			}
			std::sort(totalPeaks.begin(),totalPeaks.end());
			for( int i=0; i< totalPeaks.size(); i++){
				if( totalPeaks[i] < drop){
					index ++;
				}
			}
			int a;
			int b;
			if( histo -> GetBinContent( histo -> FindBin(totalPeaks[index])) > histo -> GetBinContent( histo -> FindBin(totalPeaks[index-1]))){
				a = index;
				b = index -1;
			}
			if( histo -> GetBinContent( histo -> FindBin(totalPeaks[index])) < histo -> GetBinContent( histo -> FindBin(totalPeaks[index-1]))){
				b = index;
				a = index -1;
			}
 
			if( histo -> GetBinContent( histo -> FindBin(totalPeaks[a])) / histo -> GetBinContent( histo -> FindBin(totalPeaks[b])) < 10){
			realPeaks.push_back(totalPeaks[index-1]);
			realPeaks.push_back(totalPeaks[index]);
			}
			if( histo -> GetBinContent( histo -> FindBin(totalPeaks[a])) / histo -> GetBinContent( histo -> FindBin(totalPeaks[b])) > 10){
			realPeaks.push_back(totalPeaks[index-2]);
			realPeaks.push_back(totalPeaks[index-1]);
			}
			
		}
	}
	histo -> GetXaxis() -> SetRangeUser(0, 50);
	if (realPeaks.size() == 1 || realPeaks.size() == 2 ){
			fitFunc = new TF1( "fitFunc", "gaus", realPeaks[0] -0.5, realPeaks[0] +0.5);
			histo-> Fit(fitFunc, "NQR");
			TF1* fitFunc_1173 = new TF1("fitFunc_1173","gaus",fitFunc->GetParameter(1)-fitFunc->GetParameter(2)*0.75, fitFunc->GetParameter(1)+fitFunc->GetParameter(2)*0.5);
			fitFunc_1173 -> SetParameter(0, fitFunc->GetParameter(0));
			fitFunc_1173 -> SetParameter(1, fitFunc->GetParameter(1));
			fitFunc_1173 -> SetParameter(2, fitFunc->GetParameter(2));
			fitFunc_1173 -> SetLineColor(kBlack);
			fitFunc_1173 -> SetLineWidth(2);
			histo -> Fit(fitFunc_1173,"QRS+"); 
			res["1.173 MeV"] = std::make_pair(fitFunc_1173->GetParameter(1),fitFunc_1173->GetParameter(2));
			//histo -> GetXaxis() -> SetRangeUser(0.,3.*fitFunc_1173->GetParameter(1));			
		}
  //Fitting 1332 keV peak
  if (realPeaks.size() == 2){  
		TF1* fitFunc2 = new TF1( "fitFunc2", "gaus", realPeaks[1] -0.5, realPeaks[1] +0.5);
		histo-> Fit(fitFunc2, "NQR");
		TF1* fitFunc_1332 = new TF1("fitFunc_1332","gaus",fitFunc2->GetParameter(1)-fitFunc2->GetParameter(2)*0.75, fitFunc2->GetParameter(1)+fitFunc2->GetParameter(2)*0.75);
			fitFunc_1332 -> SetParameter(0, fitFunc2->GetParameter(0));
			fitFunc_1332 -> SetParameter(1, fitFunc2->GetParameter(1));
			fitFunc_1332 -> SetParameter(2, fitFunc2->GetParameter(2));
			fitFunc_1332 -> SetLineColor(kBlack);
			fitFunc_1332 -> SetLineWidth(2);
			histo -> Fit(fitFunc_1332,"QRS+");		
	
		res["1.332 MeV"] = std::make_pair(fitFunc_1332->GetParameter(1),fitFunc_1332->GetParameter(2));
  
  }


  //Locating and drawing ranges of interest
  if ( realPeaks.size() > 0 && realPeaks.size() < 3 ){
		ranges->push_back(res["1.173 MeV"].first-12.*res["1.173 MeV"].second);
		ranges->push_back(res["1.173 MeV"].first-10.5*res["1.173 MeV"].second);
		ranges->push_back(res["1.173 MeV"].first-9.*res["1.173 MeV"].second);
		ranges->push_back(res["1.173 MeV"].first-7.5*res["1.173 MeV"].second);
		ranges->push_back(res["1.173 MeV"].first-6.*res["1.173 MeV"].second);
  	ranges->push_back(res["1.173 MeV"].first-4.5*res["1.173 MeV"].second);
  	ranges->push_back(res["1.173 MeV"].first-3.*res["1.173 MeV"].second);
		ranges->push_back(res["1.173 MeV"].first-1.5*res["1.173 MeV"].second);
   }	
  if ( realPeaks.size() == 2 ){
		ranges->push_back(0.5*(res["1.173 MeV"].first+4.*res["1.173 MeV"].second+res["1.332 MeV"].first-4.*res["1.332 MeV"].second));
		ranges->push_back(res["1.332 MeV"].first+1.5*res["1.332 MeV"].second);
		ranges->push_back(res["1.332 MeV"].first+3.*res["1.332 MeV"].second);
  }
  //Not Co60 spectrum controll
  if (realPeaks.size() == 0 || realPeaks.size() > 2 ){
		std::cout << "Errore" << std::endl;
		res["1.173 MeV"] = std::make_pair(-9999,0);
  }
	if( realPeaks.size() == 1 || realPeaks.size() == 2){
		for(auto range: (*ranges)){
			TLine* line = new TLine(range,3.,range,histo->GetBinContent(histo->FindBin(range)));
			line -> SetLineWidth(1);
			line -> SetLineStyle(7);
			line -> Draw("same");    
		}
	}
  return res;
}
