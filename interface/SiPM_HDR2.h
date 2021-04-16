#include <string>
#include <cmath>

float PDE_vs_OV(const float& ov, const std::string& type = "HDR2")
{
  if( type == "HDR2" )
    return 0.394321 * ( 1. - exp(-1.*0.738063*ov) );
  
  if( type == "FBK_W7S" )
    return 0.391748 * ( 1. - exp(-1.*0.288262*ov) );
  
  return 0.;
}

float Gain_vs_OV(const float& ov, const std::string& type = "HDR2")
{
  if( type == "HDR2" )
    return 36890.225 + 97602.904*ov;
    
  if( type == "FBK_W7S" )
    //return 29131.588 + 104567.57 * ov;
    //return 54174.9 + 95169.2 * ov;
    return 50738.5 + 95149 * ov;
  
  return 0.;
}

float ECF_vs_OV(const float& ov, const std::string& type = "HDR2")
{
  if( type == "HDR2" )
    return 1. + 0.0030812903 * ov + 0.0015623938 * ov * ov;
  
  if( type == "FBK_W7S" )
    //return 1. - 0.006356446 * ov + 0.0057041163 * ov * ov;
    return 1.02366 - 0.00168754 * ov + 0.00284079 * ov * ov;
  
  return 0.;
}
