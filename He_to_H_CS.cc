
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * He_to_H_CS.cc *                               galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//                                    *** I.Moskalenko, ver. 12 Jun., 2000 ***
// For any given primary (IZI,IAI) and secondary (IZF,IAF) nuclei calculates
// ratios of the total cross sections (He+AI)/(H+AI) and (He+AI->AF)/(H+AI->AF).
//    INPUT:
// E1      - energy of the primary (GeV/nucleon);
// IZI,IAI - primary charge and atomic number;
// IZF,IAF - secondary charge and atomic number.
//    OUTPUT:
// CSratio     - ratio of the cross sections (He+AI->AF)/(H+AI->AF);
// CStot_ratio - ratio of the total cross sections (He+AI)/(H+AI).
//    REFERENCES:
// Ferrando P. et al. 1988, PRC 37, 1490
//    MODIFICATIONS OF THE ORIGINAL SCHEME:
// - AMU,DELTA are found from linear extrapolation for E1<E(1);
// - for linear inter-/extra-polation of DELTA the log scale in E is applied;
// - FZI is found from linear inter-/extra-polation on log scale in Z.
//=="====!===="====!===="====!===="====!===="====!===="====!===="====!===="====!

using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"

//int He_to_H_CS(double,int,int,int,int,double*,double*); main()

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int test_He_to_H_CS()
{
   cout<<" >>>>test_He_to_H_CS "<<endl;
// this can be compared to output from F77 original in He_to_H_CS.f

   int IZI=26, IAI=56, IAF=55;
   double E=2.10, CSratio,CStot_ratio,factor;

   for(int IZF = 4; IZF<26; IZF++)
   {
      He_to_H_CS(E,IZI,IAI,IZF,IAF,&CSratio,&CStot_ratio);
      factor =pow( pow(4,3./8.)+pow(IAI,3./8.)-1. ,2)/pow(IAI,3./4.);
      cout<<E<<" "<<IZI<<" "<<IAI<<" "<<IZF<<" "<<IAF<<" "
	  <<CSratio<<" "<<CStot_ratio<<" "<<endl;
   }
   cout<<" <<<<test_He_to_H_CS "<<endl;
   return 0;
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

double FI(double X,double X1,double X2,double F1,double F2) // linear interpolation
{
   return ((F1-F2)*X+X1*F2-X2*F1)/(X1-X2);
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int He_to_H_CS(double E1,int IZI,int IAI,int IZF,int IAF,double* CSratio,double* CStot_ratio)
{
   double X,Y,a,b;
   double AMU,DELTA,FZI;
   double  E[]= {-999,  0.43,0.73,1.51},  // zeroth element dummy to conform to F77 original
           d[]= {-999,  2.45,3.45,4.40},
	  am[]= {-999,  0.082,0.047,0.031},
	   Z[]= {-999,  6.,8.,26.},
           F[]= {-999, -0.76,-0.41,+1.00};

   *CSratio    = 0.;
   *CStot_ratio= 2.10/pow(IAI,0.055);
   if(IAI <= IAF) return 0;

   if(E1 < E[2])          // interpolation/extrapolation
   {
      AMU  = FI(E1,E[1],E[2],am[1],am[2]);
      DELTA= FI( log(E1), log(E[1]), log(E[2]),d[1],d[2]);  //log scale
   } else
   {
      AMU  = FI(E1,E[2],E[3],am[2],am[3]);
      DELTA= FI( log(E1), log(E[2]), log(E[3]),d[2],d[3]);  //log scale
   }    
   if(E1 < E[1])          // asymptotics E1->0
   {
      AMU  = am[1];
      DELTA=  d[1];
   }             
   if(E1 > E[3])          // asymptotics E1->inf
   {
      AMU  = am[3];
      DELTA=  d[3];
   }             

   if(IZI < Z[2]) FZI = FI( log(IZI*1.), log(Z[1]), log(Z[2]),F[1],F[2]); // inter-/extra-polation on log scale
   else           FZI = FI( log(IZI*1.), log(Z[2]), log(Z[3]),F[2],F[3]);   

   if(IZI-IZF <= IZI/2) *CSratio = exp(AMU*pow(fabs(IZI-IZF-FZI*DELTA),1.43));
   else                 // if (IZI-IZF) > IZI/2, then another approximation (requires continuity in Z)
   {
         X =   exp(AMU*pow(fabs(IZI/2-  FZI*DELTA), 1.43));      //value
	 Y = X-exp(AMU*pow(fabs(IZI/2-1-FZI*DELTA), 1.43));      // derivative in dZ
         b = IZI/2*Y/X;
         a = X/pow(IZI/2,b);
         *CSratio = a*pow(IZI-IZF,b);
   }
   return 0;
}
     
