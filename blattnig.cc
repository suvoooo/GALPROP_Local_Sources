//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * blattnig.cc *                                 galprop package * 2005/02/15 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include "constants.h"
#define max(a,b) (((a) > (b)) ? (a) : (b))
using namespace std;
                                                                                
extern "C" void sim1_(double*,double*,double*,double*,double*,double(double*),double*);
extern "C" double fepi_(double*);
double blattnig_gamma(double, double, int, int, int);
double blattnig_pi0(double, double, int, int, int);
int kinematic(int,int,char*,double&,double&,double&,double&,double&,double&,int);

double Pp1;
int NA1,NA2;

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// differential cross section for pp->pi0,barn/GeV (Blattnig etal. PRD 62,094030)
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
double blattnig_pi0(double Epi, double Pp1, int NA1, int NA2, int key)
{ 
  double se,Ekin,Pp=Pp1/NA1;
  int test=0;

  if(Pp<=Pth0) return (0.); // the lowest threshold momentum for pp->pi0+X
  Ekin=sqrt(Pp*Pp+Mp*Mp)-Mp;
// differential cross section for pp->pi0, barn/GeV
  se=1.e-3 *exp(-5.8 -1.82*pow(Ekin,-0.4) +13.5*pow(Epi-Mpi0,-0.2) -4.5*pow(Epi-Mpi0,-0.4));
  key=0;
  return (se);
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// differential cross section for pp->pi0->gamma,barn/GeV
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
double blattnig_gamma(double Esec, double Pp2, int NA11, int NA21, int key)
{
  double h,reps=0.01,aeps=1.e-40,ai,Pp,s,Mx,gam_cms,betgam_cms,
    gam_pi_max,betgam_pi_max,EL,EU,se;

  Pp1=Pp2; NA1=NA11; NA2=NA21; 
  Pp = Pp1/NA1;                         // momentum per nucleon
  if(Pp<=Pth0) return (0.);             // the lowest threshold momentum for pp->pi0+X
  s = 2.*Mp*(sqrt(Pp*Pp+Mp*Mp)+Mp);     // Total pp-energy in CMS square
  gam_cms = sqrt(s)/2./Mp;              // CMS Lorentz factor (Lf)
  betgam_cms = sqrt(pow(gam_cms,2)-1.); // CMS beta*gamma
  Mx = 2.*Mp;                           // The mass in the channel X (=2Mp)
  gam_pi_max = (s-Mx*Mx+Mpi0*Mpi0)/2./Mpi0/sqrt(s); // max pion Lf in CMS
  if(gam_pi_max <= 1.) gam_pi_max = 1.;
  betgam_pi_max = sqrt(pow(gam_pi_max,2)-1.);       // max beta*gamma -"-
  EL = max(Mpi0,Esec+Mpi0*Mpi0/4./Esec);            // Lower limit for pion LS Energy
  EU = Mpi0*(gam_cms*gam_pi_max+betgam_cms*betgam_pi_max); // Upper limit
  if(EL > EU) return(0.);
  h=EL/10.;                                          // initial integration step
  sim1_(&EU,&EL,&h,&reps,&aeps,&fepi_,&ai);            // integration
  se=-ai*pow((pow(NA1,3./8.) +pow(NA2,3./8.)-1.),2); // a correction for A+A interaction
  key=0;
  return(se);
}
 
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// A routine to integrate by dEpi (x=Epi -total pion energy, GeV)
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
extern "C" double fepi_(double* x)
{
// (2/sqrt(..)) comes from the distr. of gammas from pi0-decay
  return (blattnig_pi0(*x,Pp1,NA1,NA2,0) *2./sqrt(pow(*x,2)-Mpi0*Mpi0));
}
 
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// Test routine, requires integr.f, constants.h
// g++ -c blattnig.cc
// g++ *.o -L/usr/lib/gcc-lib/i386-redhat-linux/3.3.2 
// g++ *.o -L/usr/lib/gcc-lib/i386-redhat-linux/2.96/collect2
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

/*
int main()
{
  cout<<"<<test_blattnig>>"<<endl;

  double Esec,sum,Pp, Ekin0, Ep, deltaEs=1.4, deltaEp=1.2, Ekin;
  int i,j,key=0;

  Pp1=20.; NA1=1; NA2=1;
  Pp=Pp1/NA1;

  for(cout<<"# Neutral pions"<<endl, i=1;i<20;i++)
    {
      Esec=i;
      sum=blattnig_pi0(Esec,Pp1,NA1,NA2,key);
      cout<<"# "<<Pp<<" "<<Esec<<" "<<sum<<endl;
    }

  for(cout<<"# Gammas"<<endl, i=1;i<20;i++)
    {
      Esec=i;
      sum=blattnig_gamma(Esec,Pp1,NA1,NA2,key);
      cout<<"  "<<Pp<<" "<<Esec<<" "<<sum<<endl;
    }

// calculation of a test spectrum using CR protons ~2.2 E^-2.75
  for(cout<<"# Gammas: test spectrum"<<endl, Ekin0=0.3, Esec=1.e-2, i=1;i<40;i++, Esec*=deltaEs)
    {
      for(Ekin=1.e6, sum=0.; Ekin>Esec && Ekin>Ekin0;)
	{

	  Pp1=sqrt(Ekin*Ekin+2.*Ekin*Mp);
          Ep=Ekin+Mp;
	  sum+=pow(Ep,-2.75) *(Ekin*log(deltaEp)) *blattnig_gamma(Esec,Pp1,NA1,NA2,key);
//	  cout<<"  "<<Pp1<<" "<<Esec<<" "<<blattnig_gamma(Esec,Pp1,NA1,NA2,key)<<endl;
	  Ekin/=deltaEp;
	} // factor 1.45 to take into account heavier nuclei; 1.e-24 -conversion to cm^-2
      cout<<"  "<<Esec<<" "<<1.45 *2.2 *4.*Pi *sum*1.e-24
                      <<" "<<1.45 *2.2 *4.*Pi *sum*1.e-24 *pow(Esec,2.75)<<endl;
    }
}


*/


















