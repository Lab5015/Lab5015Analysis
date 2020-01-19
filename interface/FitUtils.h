#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <string>
#include <map>

#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TPad.h"
#include "TLatex.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLine.h"


/*** find effective sigma ***/
void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction);

/*** draw time resolution plot ***/
struct CTRResult
{
  float effSigma;
  float gausSigma;
  float gausSigmaErr;
};
