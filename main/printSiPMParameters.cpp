#include "interface/SiPM_HDR2.h"

#include <iostream>
#include <string>


int main(int argc, char** argv)
{
  if( argc < 3 )
  {
    std::cout << ">>> printSiPMParameters::usage:   " << argv[0] << " SiPM type   ov" << std::endl;
    return -1;
  }
  
  
  //--- parse the parameters
  std::string SiPMType(argv[1]);
  float Vov(atof(argv[2]));
  
  
  std::cout << "  G(" << Vov << ") = " << Gain_vs_OV(Vov,SiPMType)/1E05 << "E05" << std::endl;
  std::cout << "PDE(" << Vov << ") = " << PDE_vs_OV(Vov,SiPMType) << std::endl;
  std::cout << "ECF(" << Vov << ") = " << ECF_vs_OV(Vov,SiPMType) << std::endl;
  std::cout << "GxPDExECF(" << Vov << ") = " << Gain_vs_OV(Vov,SiPMType)*PDE_vs_OV(Vov,SiPMType)*ECF_vs_OV(Vov,SiPMType) << std::endl;
}
