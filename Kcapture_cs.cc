//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Kcapture_cs.cc *                              galprop package * 2001/08/16
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// Kcapture_cs                              ### I.Moskalenko ###  8/16/2001 ###
// calculation of the electron attachment and stripping cross sections.  
//    INPUT:
// Ek -kinetic energy per nucleon of the beam momentum particles (MeV);//AWS20010829
// Zp,Zt -the charge of the beam and target nuclei;
//    OUTPUT: 
// attach  [mbarn] is the attachment cross section;
// strip   [mbarn] is the stripping cross section;
//    REFERENCES:
// Crawford H.J. 1979, PhD Thesis UC at Berkeley
// Pratt R.H. et al. 1973, Rev.Mod.Phys. 45, 273
// Wilson L.W. 1978, PhD Thesis UC at Berkeley

using namespace std;
#include<cmath>
#include"constants.h"

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void Kcapture_cs(double Ek, int Zp, int Zt, double *attach, double *strip)
{
   double gam = 1.+Ek*1.e-3/amu, sT =1.e27 * 8./3.*Pi*pow(Rele,2), //AWS20010829
          beta,Mbet,Nbet,a,fcor;

   beta = sqrt(1.-pow(gam,-2));
   Mbet = 4./3.+gam*(gam-2.)/(gam+1.)*(1.-log((1.+beta)/(1.-beta))/2./beta/pow(gam,2));
   Nbet = pow(beta,-3) *((-4.*gam +34. -63./gam +25./pow(gam,2) +8./pow(gam,3))/15.
	   -(gam-2.)*(gam-1.)/2./beta/pow(gam,3) *log((1.+beta)/(1.-beta)));
   a = Zp*ALPHAf;
   fcor = pow(a,2.*(sqrt(1.-a*a)-1.)) *exp(-2.*a/beta*acos(a)) *(1.+Pi*a*Nbet/Mbet);
//    factor 1.202 accounts for contribution of states higher than 1s
   *attach = 1.202 *fcor *3./2.*sT*pow(1.*Zp,5)*Zt*pow(ALPHAf,4) *beta*gam *pow(gam-1.,-3)*Mbet;
   *strip = 3./2.*sT*pow(Zp*beta*ALPHAf,-2) 
            *0.285*(2.*log(2./0.219*beta*gam/Zp/ALPHAf)-pow(beta,2))*Zt*(Zt+1.);
   return;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
/*
#include<stdio.h>
#include<iostream.h>
void Kcapture_cs(double,int,int,double*,double*);
main() //AWS20010829 this routine not yet modified for GeV->MeV
{
//   FILE *f1;
   double Ek,attach,strip, ttt;
   int i,Zp[]={4,18,20,22,23,24,25,26,27,28};

//   f1=fopen("nucleon.out","w");                     // print data to a file
//   fprintf(f1,"%3s%11s%11s%11s%11s%11s%11s%11s\n", "# A",
//      "Ek, GeV","pp_inel","pA_inel","app_non","apA_non","app_ann","apA_ann");
   printf(  "\n%3s%13s%13s%13s", "\\Ek"," 1 GeV"," 3 GeV"," 5 GeV");
   printf(  "\n%3s%13s%13s%13s", "A \\","cs_att","cs_att","cs_att");

   for(i=0; i<10.; i++)
   {
      printf("\n%3d",Zp[i]);
      for(Ek=1.e3; Ek<6.e3; Ek+=2.)
      {
         Kcapture_cs(Ek,Zp[i],1,&attach,&strip);
// factor 1.148 accounts for contribution of heavier than H atoms in the ISM
// (with ISM abundances from Cameron 1970, SSR 15, 121):
// factor = SUM(i) Z_i (n_i/n_H).
         printf("%13.3e",attach *1.148);
      }
   }
   printf("\n");
   printf("%.4f   %.4f\n",1.+2.*0.0695,1.+2.*0.0695+5.*0.001); 

   Kcapture_cs(140.,6,13,&attach,&strip);
   printf("C(140 MeV/nucleon) + Al: %14.3e bn\n",strip *1.e-27);
   Kcapture_cs(250.,6,6,&attach,&strip);
   printf("C(250 MeV/nucleon) + C:  %14.3e bn\n",strip *1.e-27);
   Kcapture_cs(250.,6,13,&attach,&strip);
   printf("C(250 MeV/nucleon) + Al: %14.3e bn\n",strip *1.e-27);
   Kcapture_cs(1050.,10,6,&attach,&strip);
   printf("Ne(1.05 GeV/nucleon)+C:  %14.3e bn\n",strip *1.e-27);
   Kcapture_cs(1050.,10,1,&attach,&strip);
   printf("Ne(1.05 GeV/nucleon)+H:  %14.3e bn\n",strip *1.e-27);
}
*/





