#include "interface/Co60SpectrumAnalyzer_2Peaks.h"



std::map<std::string,std::pair<float,float> > Co60SpectrumAnalyzer_2Peaks(TH1F* histo,
                                                                   std::vector<float>* ranges)
{
  gStyle -> SetOptFit(0);
  
  
  std::map<std::string,std::pair<float,float> > res;
  


	
	std::vector<double> totalPeaks;
	std::vector<double> totalPeaks2;
	std::vector<double> realPeaks;
  
  int startBin = 2;
  int endBin;
  
  while ( histo -> GetBinContent(startBin) < 100 && histo -> GetBinContent(startBin-1) < 100 ){
	startBin ++;
  }

  endBin =  histo->GetNbinsX() /6.;
  while ( histo -> GetBinContent(endBin) > 50 && histo -> GetBinContent(endBin-1) > 50 && histo -> GetBinContent(endBin+1) > 50 && endBin < histo->GetNbinsX() -50){
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
  double drop = histo->GetBinCenter(startBin+1);
  
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
 
  
  if( rangesDerivative.size() > 10){
		
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
			drop = histo -> GetBinCenter(startBin + (minIndex+5) *5);
			if(drop < float((histo->GetBinCenter(endBin)-histo->GetBinCenter(startBin)))/3.){
				check = 0;
			}
			if (check == 0){
				erasedIndex.push_back(minIndex);
			}
			sort(erasedIndex.begin(), erasedIndex.end());	
			
			
		}
		
		
		TF1* fitFunc;
		

		if (drop > histo->GetBinCenter((endBin-startBin)/3.) ){
			//Searching 1173 and 1332 peaks
			int nPeaks = 4;
			int nFound = 0;
			TSpectrum* spectrum = new TSpectrum(nPeaks);
			float start =  histo->GetBinCenter(startBin);
			float stop =   histo->GetBinCenter(endBin);
			histo -> GetXaxis() -> SetRangeUser(start, stop);
			std::cout<<"start"<<start<<"stop"<<stop<<std::endl;			
                                                                                                                                                                                                                                                                                                                                                                       
			nFound = spectrum -> Search(histo, 1,"" , 0.001);
			//std::cout<<"nFound"<<nFound<<std::endl;
			double* peaks = spectrum -> GetPositionX();
			int index = -1;
			if (nFound > 1){
				for( int i = 0; i < nFound; i++){
					totalPeaks.push_back(peaks[i]);
					//std::cout<<"totalPeaks "<<i<<" "<<totalPeaks[i]<<std::endl;
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
		histo -> GetXaxis() -> SetRangeUser(0, 35);
		if (realPeaks.size() == 1 || realPeaks.size() == 2 ){
				fitFunc = new TF1( "fitFunc", "gaus", realPeaks[0] -0.5, realPeaks[0] +0.5);
				histo-> Fit(fitFunc, "NQR");
				TF1* fitFunc_1173 = new TF1("fitFunc_1173","gaus",fitFunc->GetParameter(1)-fitFunc->GetParameter(2)*0.75, fitFunc->GetParameter(1)+fitFunc->GetParameter(2)*0.9);
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
			TF1* fitFunc_1332 = new TF1("fitFunc_1332","gaus",fitFunc2->GetParameter(1)-fitFunc2->GetParameter(2)*0.75, fitFunc2->GetParameter(1)+fitFunc2->GetParameter(2)*0.9);
				fitFunc_1332 -> SetParameter(0, fitFunc2->GetParameter(0));
				fitFunc_1332 -> SetParameter(1, fitFunc2->GetParameter(1));
				fitFunc_1332 -> SetParameter(2, fitFunc2->GetParameter(2));
				fitFunc_1332 -> SetLineColor(kBlack);
				fitFunc_1332 -> SetLineWidth(2);
				histo -> Fit(fitFunc_1332,"QRS+");		
		
			res["1.332 MeV"] = std::make_pair(fitFunc_1332->GetParameter(1),fitFunc_1332->GetParameter(2));
		
		}


		int nPeaks2 = 2;
			int nFound2 = 0;
			TSpectrum* spectrum2 = new TSpectrum(nPeaks2);
			float start2 =  histo->GetBinCenter(endBin);
			float stop2 =  histo->GetBinCenter(endBin)+100; 
			std::cout<<"start2"<<start2<<"end2"<<stop2<<std::endl;
			histo -> GetXaxis() -> SetRangeUser(start2, stop2);
			                                                                                                                                                                                                                                                                                                                                                                         
			nFound2 = spectrum2 -> Search(histo, 1,"" , 0.001);
			std::cout<<"nFound2"<<nFound2<<std::endl;
			double* peaks2 = spectrum2 -> GetPositionX();
			int index2 = -1;
			if (nFound2 > 0){
				for( int i = 0; i < nFound2; i++){
					totalPeaks2.push_back(peaks2[i]);
					std::cout<<"totalPeaks2 "<<i<<" "<<totalPeaks2[i]<<std::endl;
				}
			}

		//Fitting 2.5 MeV peak
		if (totalPeaks2.size() == 1){  
			TF1* fitFunc3 = new TF1( "fitFunc3", "gaus", totalPeaks2[0]-0.5, totalPeaks2[0]+2.);
			fitFunc3 -> SetLineColor(kBlack);
			fitFunc3 -> SetLineWidth(2);
			histo-> Fit(fitFunc3, "NQRS+");
			TF1* fitFunc_2505 = new TF1("fitFunc_2505","gaus",fitFunc3->GetParameter(1)-0.75*fitFunc3->GetParameter(2), totalPeaks2[0]+2.);
				fitFunc_2505 -> SetParameter(0, fitFunc3->GetParameter(0));
				fitFunc_2505 -> SetParameter(1, fitFunc3->GetParameter(1));
				fitFunc_2505 -> SetParameter(2, fitFunc3->GetParameter(2));
				fitFunc_2505 -> SetLineColor(kBlack);
				fitFunc_2505 -> SetLineWidth(2);

				histo -> Fit(fitFunc_2505,"QRS+");

		
			res["2.505 MeV"] = std::make_pair(fitFunc_2505->GetParameter(1),fitFunc_2505->GetParameter(2));
		}	

		double mean = 0;
		if (totalPeaks2.size() == 2){
			if (totalPeaks2[0]>totalPeaks2[1]) mean = totalPeaks2[0];
			else mean = totalPeaks2[1]; 
			TF1* fitFunc3 = new TF1( "fitFunc3", "gaus", mean-0.5, mean+2.);
			fitFunc3 -> SetLineColor(kBlack);
			fitFunc3 -> SetLineWidth(2);
			histo-> Fit(fitFunc3, "NQRS+");
			TF1* fitFunc_2505 = new TF1("fitFunc_2505","gaus",fitFunc3->GetParameter(1)-0.75*fitFunc3->GetParameter(2), totalPeaks2[0]+2.);
				fitFunc_2505 -> SetParameter(0, fitFunc3->GetParameter(0));
				fitFunc_2505 -> SetParameter(1, fitFunc3->GetParameter(1));
				fitFunc_2505 -> SetParameter(2, fitFunc3->GetParameter(2));
				fitFunc_2505 -> SetLineColor(kBlack);
				fitFunc_2505 -> SetLineWidth(2);

				histo -> Fit(fitFunc_2505,"QRS+");

		
			res["2.505 MeV"] = std::make_pair(fitFunc_2505->GetParameter(1),fitFunc_2505->GetParameter(2));
		}
		
	}



			



 //Not Co60 spectrum controll
  if (realPeaks.size() == 0 || realPeaks.size() > 2 ){
	std::cout << "Errore" << std::endl;
	res["1.173 MeV"] = std::make_pair(-9999,0);
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
		//ranges->push_back(0.5*(res["1.173 MeV"].first+4.*res["1.173 MeV"].second+res["1.332 MeV"].first-4.*res["1.332 MeV"].second));
		ranges->push_back(0.5*(res["1.173 MeV"].first-1.5*res["1.173 MeV"].second+res["1.332 MeV"].first+1.5*res["1.332 MeV"].second));
		ranges->push_back(res["1.332 MeV"].first+1.5*res["1.332 MeV"].second);
		ranges->push_back(res["2.505 MeV"].first-1.*res["2.505 MeV"].second);
		ranges->push_back(res["2.505 MeV"].first+1.*res["2.505 MeV"].second);
  }
  if (realPeaks.size() == 0){
		ranges -> push_back(-1);
	}

  if ( realPeaks.size() == 1 ){
		float second = res["1.173 MeV"].first *1.3;
		ranges->push_back(0.5*(res["1.173 MeV"].first+0.4*res["1.173 MeV"].second-0.2*second));
		ranges->push_back(second-0.2*second);
		ranges->push_back(second-0.1*second);
		ranges->push_back(second+0.1*second);
		//ranges->push_back((histo->GetBinCenter(ultimo_bin)-(res["1.275 MeV"].first+0.1*res["1.275 MeV"].first))/3+res["1.275 MeV"].first+0.1*res["1.275 MeV"].first);
		//ranges->push_back((res["1.275 MeV"].first)*1.25+0.2*(res["1.275 MeV"].first)*1.25);
		//ranges->push_back(histo->GetBinCenter(ultimo_bin));
   }

    
  for(auto range: (*ranges)){
		
		if ( range  < 0  && realPeaks.size() != 2) { //error controll
			std::cout << "Errore" << std::endl;
    	res["1.173 MeV"] = std::make_pair(-9999,0);	
		}
	  float yval = std::max(10., histo->GetBinContent(histo->FindBin(range)));
		TLine* line = new TLine(range,3.,range, yval);
		line -> SetLineWidth(1);
		line -> SetLineStyle(7);
		line -> Draw("same");
    
  }
  


    
  
  
  return res;
}
