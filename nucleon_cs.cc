//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * nucleon_cs.cc *                               galprop package * 2001/05/11 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// parametrization of the pp-, pA-, AA-, (anti_p)p-, and (anti_p)A total
// inelastic cross sections, as well as (anti_p) annihilation cross sect.
// A sample driver program is in the bottom of this file.
//    INPUT:
// option -pA total inelastic cross section 
//       =0 -[L83]; 
//       =1 -[WA96] for Zt>5 and [BP01] for Zt<=5;
//       =2 -[BP01];
// Ek -kinetic energy per nucleon of the beam momentum particles (GeV);
// Zp=+/-1 is the charge of beam and At is the atomic number of target nuclei
// (Zp=At=1 for pp-collisions, Zp= -1 for antiprotons);
//    OUTPUT: 
// PP_inel  [mbarn] is the proton-proton inelastic cross sect.;
// PA_inel  [mbarn] is the proton-nucleus inelastic cross sect. (p+A2);
// aPP_non  [mbarn] is the antiproton-proton inelastic cross sect.(nonannihil.);
// aPA_non  [mbarn] put =PA_inel, is the antiproton-nucleus inelastic cross sect.(nonan.);
// aPP_ann  [mbarn] is the antiproton-proton annihilation cross sect.,
// aPA_ann  [mbarn] is the antiproton-nucleus annihilation cross sect.,
//    REFERENCES:
// [W79] Westfall et al. 1979, PRC, 19, 1309
// [L83] Letaw et al. 1983, ApJS, 51, 271
// [TN83] Tan & Ng 1983, J.Phys.G:Nucl.Phys., 9, 227
// [PDG00] D.E.Groom et al. 2000, Europ. Phys. J. C15, 1 
// [MO97] Moiseev & Ormes 1997, Astropart. Phys.6, 379
// [WA96] Wellisch H.P., Axen D. 1996, PRC 54, 1329; Wellish 1999, private comm. 
//    (typo corrections & code)
// [BP01] V.S.Barashenkov, A.Polanski code used for calc. of the pA total inelastic 
//    cross sections

using namespace std;//AWS20050624
#include<cmath>
#include"fort_interface.h"
#include"constants.h"

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void nucleon_cs(int option, double Ek, int Zp, int Zt, int At,
   double *PP_inel,double *PA_inel,double *aPP_non,double *aPA_non,double *aPP_ann,double *aPA_ann)
{

   *PP_inel= *PA_inel= *aPP_non= *aPA_non= *aPP_ann= *aPA_ann=0.;
   if(Ek <= 0.) return;

   double U,Cp,C1,s2,Z,A, aPP_inel=0.,aPA_inel=0., Emev=Ek, b0,Fcorr,rN,s0,p1,p2,p3,p4,p5,x,f1,f2;
   double PZ=fabs(1.*Zp),Em=1000.*Ek, TZ=Zt, TA=At;
   int ISS=2;

//## Proton-Proton INELASTIC cross section, mb [TN83,p.234]
   if(Ek > 0.3)
   {
      U = log((Ek+Mp)/200.);
      *PP_inel = 32.2*(1.+0.0273*U);
      if(U >= 0.) *PP_inel += 32.2*0.01*U*U;
      if(Ek < 3.) *PP_inel /= 1.+0.00262/pow(Ek,17.9+13.8*log(Ek)+4.41*pow(log(Ek),2));
   }
   if(Zp*At == 1) return;

//## Proton-Nucleus INELASTIC cross section, mb
   switch(option)
   {
   case 0:                                        // [L83]
      C1 = (At == 4) ? 0.8 : 1.;                  // a correction for He
      if(At == 9) C1 = 1.+0.75*exp(-Ek/0.075);    // a correction for Be
      *PA_inel = C1 *45. *pow(TA,0.7) *(1.+0.016*sin(5.3-2.63*log(TA)));
      if(Ek < 5.) *PA_inel *= 1.-0.62 *exp(-Ek/0.2) *sin(10.9/pow(Ek*1.e3,0.28));
      if(At == 4) *PA_inel = (Ek > 0.01) ?        // pHe, my fit
         111.*(1.-(1.-sin(9.72*pow(log10(Ek*1000.),0.319)-4.14))*exp(-3.84*(Ek-0.1))) : 0.;
      break;

   case 1:                                        // Zt>5 [WA96], Zt<=5 [BP01]
      if(Zt>5)
      {
         b0 = 2.247-0.915*(1.-pow(TA,-1./3.));
         Fcorr = (1.-0.15*exp(-Emev))/(1.-0.0007*At);   // high-energy correction
         rN = (At-Zt>1.5) ? log(TA-Zt) : 1.;
         s0 = Pi*10.*pow(1.36,2.)*Fcorr*rN*(1.+pow(TA,1./3.)-b0*(1.-pow(TA,-1./3.)));
         p1 = 8.-8./At-0.008*At;
         p2 = 2.*(1.17-2.7/At-0.0014*At);
         p3 = 0.8+18./At-0.002*At;
         p4 = 5.6-0.016*At;
         p5 = 1.37*(1.+1./At);
         x = log10(Emev);
         f1 = 1./(1.+exp(-p1*(x+p2)));                 // low-energy return to zero
         f2 = 1. +p3 *( 1. -1./(1.+exp(-p4*(x+p5))) ); // low-energy rise
            *PA_inel = f1*f2*s0;
	 break;
      }

   case 2:                                        // [BP01]
   default:
      if (Em<14.)  Em=14.;
      if (Em>1.e6) Em=1.e6;
      *PA_inel = sighad_cc(ISS,PZ,PZ,TA,TZ,Em);   // IMOS20020502
   }
   if(Zp*At >= 1) return;

//## AntiProton-Proton ANNIHILATION cross section [TN83]
   if(Ek < 10.) *aPP_ann = 661.*(1.+0.0115/pow(Ek,0.774)-0.948*pow(Ek,0.0151)); // 0.1GeV<Ek<12GeV
   else
   {
// assuming aPP_ann = aPP_tot -PP_tot (i.e., aPP_elast = PP_elast); (aPP_tot-PP_tot) from [PDG00]
      s2 = 2.*Mp*(Ek+2*Mp);                   // square of the total CMS energy
      *aPP_ann = 2*35.43/pow(s2,0.560);
   }

//## AntiProton-Proton TOTAL INELASTIC cross section
   aPP_inel = *PP_inel + *aPP_ann;
   if(Ek <= 14.)
   { 
      aPP_inel = 24.7*(1.+0.584/pow(Ek,0.115)+0.856/pow(Ek,0.566));
      if(*aPP_ann > aPP_inel) *aPP_ann = aPP_inel;
   }

//## AntiProton-Proton TOTAL INELASTIC NON-ANNIHILATION cross section
   *aPP_non = aPP_inel - *aPP_ann;
   if(*aPP_non < 0.) *aPP_non = 0.;

//## AntiProton-NUCLEUS cross sections
   if(At > 1)
   {
//## AntiProton-NUCLEUS TOTAL INELASTIC NON-ANNIHILATION cross section
      *aPA_non = *PA_inel;

//## AntiProton-NUCLEUS ANNIHILATION cross section on 12C-nucleus [mb] (0.4<Pp<300) [MO97]
      A = At;
      Z = Zt;                               // Z = 0.59*pow(A,.927);  for Z > 2 nuclei
      if(At == 4) { Z = 2.; A = 3.30; }     // Modified to agree with HE p-He cs / imos
      *aPA_ann = pow(A,2./3.)               // Scaling to other nuclei
//         *(48.2 +19./pow(Ek-0.02,0.55) 
         *(48.2 +19./pow(Ek,0.55)           // modified to agree w. He@<100 MeV / imos
         +(0.1-0.18/pow(Ek,1.2))*Z +0.0012/pow(Ek,1.5)*Z*Z)  - *aPA_non;
      if(*aPA_ann < 0.) *aPA_ann = 0.;
      if(*aPA_ann < *aPP_ann) *aPA_ann = *aPP_ann;
      if(At == 4 && Ek > 5.)  *aPA_ann = *aPP_ann;
/*
//## my fit to AntiProton-NUCLEUS total cross section on 12C-nucleus [mb] (0.4<Pp<300)
         double Pp =sqrt(pow(Ek+Mp,2)-Mp*Mp); // GeV,kin. momentum per nucleon
         *aPA_ann = (Pp > 40.) ? 236.*(1.+6.9e-2*exp(-Pp/100.)) :
            209.7*pow(.542/Pp,.565)+29.6*cos(log10(Pp/9.29)*5.11)+257.9;
         *aPA_ann *= pow(At/12.,2./3.);                             // scaling to other nuclei
*/
   }
   return;
}


//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

// A SAMPLE TEST PROGRAM
/*
#include<stdio.h>
#include<iostream.h>
extern "C" void sigtap2_(int*);
extern "C" double sighad_(int*,double*,double*,double*,double*,double*);
void nucleon_cs(int,double,int,int,int,double*,double*,double*,double*,double*,double*);
main()
{
   FILE *f1,*f2;
   int i,j,Zp=-1,Zt=2,At=4,  ISS;
   double Emev,cs,PA=1.,PZ=1.,TA,TZ;
   double Ek,dE=1.2, PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann;

   ISS = -1; sigtap2_cc(&ISS); ISS = 2;

   f1=fopen("nucleon.out1","w");                     // output1
   fprintf(f1,"%3s%11s%11s%11s%11s%11s%11s%11s\n", "# A",
      "Ek, GeV","PP_inel","PA_inel","aPP_non","aPA_non","aPP_ann","aPA_ann");
   printf(  "\n%3s%11s%11s%11s%11s%11s%11s%11s\n", "# A",
      "Ek, GeV","PP_inel","PA_inel","aPP_non","aPA_non","aPP_ann","aPA_ann");

   for(i=At; i<At+1; i++)  // for plots
   {
      for(Ek=0.02; Ek<1.e4; Ek*=dE)
      {
         nucleon_cs(2,Ek,Zp,Zt,i,&PP_inel,&PA_inel,&aPP_non,&aPA_non,&aPP_ann,&aPA_ann);

         Emev=1000*Ek; TA = i; TZ = Zt;

// for plots
         fprintf(f1, "   %11.3e%11.3e%11.3e%11.3e%11.3e%11.3e%11.3e%11.3e%11.3e\n",
	    Ek,PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann,aPP_non+aPP_ann,aPA_ann+aPA_non);
// for tests
//         fprintf(f1,"%3d%3d%11.2E%11.2E%11.2E%11.2E%11.2E%11.2E%11.2E\n", Zt,
//	    i,Ek,PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann);
         printf(    "%3d%3d%11.3E%11.3E%11.3E%11.3E%11.3E%11.3E%11.3E\n", Zt,
	    i,Ek,PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann);
      }
   }
   fclose(f1);

   f2=fopen("nucleon.out2","w");                     // output2
   fprintf(f2,  "#%5s%11s%11s%11s%11s\n", "Z  A","Ek, GeV","PA_inel0","PA_inel1","PA_inel2");
   printf (   "\n#%5s%11s%11s%11s%11s\n", "Z  A","Ek, GeV","PA_inel0","PA_inel1","PA_inel2");

   int incr = 2;
   for(i=4; i<64; i+=incr)
   {
      Zt = (int) (0.59*pow(i,.927)+0.4);
      fprintf(f2,"\n#%3d%3d\n", Zt,i);
      printf (   "\n#%3d%3d\n", Zt,i);
      for(Ek=0.02; Ek<2000.; Ek*=1.1, fprintf(f2,"\n"), printf("\n"))
      {
         fprintf(f2,"      %11.2E",Ek);
         printf (   "      %11.2E",Ek);
         for(j=0; j<3; j++)
	 {
            nucleon_cs(j,Ek,Zp,Zt,i,&PP_inel,&PA_inel,&aPP_non,&aPA_non,&aPP_ann,&aPA_ann);
            fprintf(f2,"%11.2E",PA_inel);
            printf (   "%11.2E",PA_inel);
         }
      } if(i%10==0) incr+=2;
   }
   fclose(f2);
   exit(1);

   f1=fopen("nucleon.out1","w");                    // print data to a file

   Zt = 26; At = 56;
   for(Ek=0.005; Ek<2000.; Ek*=dE)
   {
      nucleon_cs(2,Ek,Zp,Zt,At,&PP_inel,&PA_inel,&aPP_non,&aPA_non,&aPP_ann,&aPA_ann);
      fprintf(f1,"%3d%3d%11.2E%11.2E%11.2E%11.2E%11.2E%11.2E%11.2E\n", Zt,
         At,Ek,PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann);
      printf(    "%3d%3d%11.3E%11.3E%11.3E%11.3E%11.3E%11.3E%11.3E\n", Zt,
         At,Ek,PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann);
   }
}
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

// Barashenkov & Polanski pA total cross section  IMOS20020502
double sighad_cc(int IS, double PA, double PZ, double TA, double TZ, double E)
{ 
   return( sighad_(&IS, &PA, &PZ, &TA, &TZ, &E) );
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

// initialization of the Barashenkov & Polanski cross section code
void sigtap_cc(int ISS)
{ 
   sigtap2_(&ISS);
}

*/
