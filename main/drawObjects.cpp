#include "interface/SetTDRStyle.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <map>
#include <vector>
#include <filesystem>

#include "TFile.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TClass.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"

int main(int argc, char** argv)
{
  setTDRStyle();
  gErrorIgnoreLevel = kError;
  if (argc < 2)
  {
    std::cout << "Usage: " << argv[0] << " config.cfg" << std::endl;
    return -1;
  }
  // - get config options
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  std::string inputFile = opts.GetOpt<std::string>("Output.outFileNameStep2");
  std::filesystem::create_directories(plotDir);

  // - open file
  TFile* f = TFile::Open(inputFile.c_str(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error: cannot open file " << inputFile << std::endl;
    return -1;
  }

  // - create objects for storage
  std::map<std::string, TH1*> h1;
  std::map<std::string, TH2*> h2;
  std::map<std::string, TProfile*> prof;
  std::multimap<std::string, TF1*> funcs;

  // - get objects from keys
  TIter next(f->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next()))
    {
      TObject* obj = key->ReadObj();
      std::string name = obj->GetName();     
      if (obj->InheritsFrom(TProfile::Class()))
	{
	  prof[name] = (TProfile*)obj;
	}
      else if (obj->InheritsFrom(TH2::Class()))
	{
	  h2[name] = (TH2*)obj;
	}
      else if (obj->InheritsFrom(TH1::Class()))
	{
	  h1[name] = (TH1*)obj;
	}
      else if (obj->InheritsFrom(TF1::Class()))
	{
	  funcs.insert({name, (TF1*)obj});
	}
    }

  // - create output dirs
  auto makeDir = [&](const std::string& s)
  {
    std::string d = plotDir + "/" + s;
    std::filesystem::create_directories(d);
    return d;
  };
  std::string dirH1 = makeDir("TH1");
  std::string dirH2 = makeDir("TH2");
  std::string dirProf = makeDir("TProfile");

  // - draw TH1 + function
  for (auto& [name, h] : h1)
  {
    TCanvas c("c", "c", 800, 600);
    h->Draw("hist");
    for (auto& [fname, f] : funcs)
    {
      if (fname.find(name) != std::string::npos)
      {
        f->SetLineColor(kRed);
        f->SetLineWidth(2);
        f->Draw("SAME");
      }
    }
    c.SaveAs((dirH1 + "/" + name + ".png").c_str());
  }

  // - draw TH2
  for (auto& [name, h] : h2)
  {
    TCanvas c("c", "c", 800, 600);
    h->Draw("colz");
    c.SaveAs((dirH2 + "/" + name + ".png").c_str());
  }

  // - draw TProfile + function
  for (auto& [name, p] : prof)
  {
    TCanvas c("c", "c", 800, 600);
    p->Draw("E1");
    for (auto& [fname, f] : funcs)
    {
      if (fname.find(name) != std::string::npos)
      {
        f->SetLineColor(kRed);
        f->SetLineWidth(2);
        f->Draw("SAME");
      }
    }
    c.SaveAs((dirProf + "/" + name + ".png").c_str());
  }
  f->Close();
  return 0;
}
