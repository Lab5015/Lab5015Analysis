#ifndef NA22SPECTRUMANALYZER_H
#define NA22SPECTRUMANALYZER_H

#include <iostream>
#include <vector>

#include "TH1F.h"
#include "TF1.h"
#include "TLine.h"
#include "TSpectrum.h"
#include "TStyle.h"



std::map<std::string,std::pair<float,float> > Na22SpectrumAnalyzer(TH1F* histo,
                                                                   std::vector<float>* ranges = NULL);

#endif
