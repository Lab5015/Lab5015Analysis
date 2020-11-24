// HDR2 SiPM parameters
float PDE_vs_OV_HDR2(const float& ov)
{
  return 39.4321 * ( 1. - exp(-1.*0.738063*ov) );
}

float Gain_vs_OV_HDR2(const float& ov)
{
  return 36890.225 + 97602.904*ov;
}

float ECF_vs_OV_HDR2(const float& ov){
  return 1. + 0.0030812903 * ov + 0.0015623938 * ov * ov;
}





// FBK SiPM parameters
/*float PDE_vs_OV(const float& ov)
{

  //return ( -1.5389053 + 22.812361*ov  - 10.091946*ov*ov + 3.6867230 *ov*ov*ov -0.73282017*ov*ov*ov*ov + 0.056861423*ov*ov*ov*ov*ov ) ;

  return 45.5275 * ( 1. - exp(-1.*0.36229*ov) );  
}

float Gain_vs_OV(const float& ov)
{
  return 29131.588 + 104567.57 * ov;

}*/
