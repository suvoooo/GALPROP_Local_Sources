
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * D_pp.cc *                                     galprop package * 2/13/2003 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// diffusion in momentum space
// formula from Seo and Ptuskin ApJ 431, 705 (1994); Ptuskin-2003 (Berezinskii et al. book)
// with w=1 (since can be subsumed in v_alfven)
// v_alfven in km s-1
// NB Dpp_constant was defined but not used in galprop v25 and earlier

using namespace std;//AWS20050624
#include<iostream>
#include<cmath>
#include"galprop_classes.h"

double Galprop::D_pp(double p,double a,double v_alfven,double D_xx)
{
   double Dpp_constant;
   
   if(galdef.diff_reacc==1) Dpp_constant=2.*2./(3.*(2.-a)*(4.-a)*(2.+a)*a); // Seo & Ptuskin
   else                     Dpp_constant=1./(a*(4.-a)*(4.-a*a));            // IMOS20030213 Ptuskin-2003
   double D_pp_=Dpp_constant*pow(p,2.0) * pow(v_alfven*1.e5,2) /D_xx;

// cout<<"D_pp "<<D_pp_<<" "<<a<<" "<<Dpp_constant<<endl;  
   return D_pp_;
}


