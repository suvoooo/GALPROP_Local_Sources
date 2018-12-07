using namespace std;
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include"galprop_classes.h"
#include"galprop_internal.h"
#define sign(a,b) (((b) > 0.) ? (a) : (-a))

double bessj(int, double);
double hyp1F1(double, double, double);
double fu(double);
//double gamma(const double);
double sim(double, double, double, double, double, double(*)(double));

//** see a sample main routine at the end **
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// Propagation of electrons in the Galaxy cylindrical geometry
// INPUT:
// Ee         = electron kinetic energy, GeV
// R, Z       = (R,Z) point coordinates, kpc
// Rg, hg     = radius and 1/2 thickness of the area filled with sources, kpc
// hh         = halo size, kpc
// elossconst = energy loss constant, 1/(GeV s); can not be =0
// elossconst=32./9.*Pi*pow(Rele/Mele,2)*C*Ugevcc; Ugevcc = ISRF energy density
// gamma_e    = electron spectral index 
// Dxx, g1    = diffusion coefficient normalization @ 1 GeV (kpc^2/s), and index
// requires file j0zero500.dat with zeros of Bessel J0 function
// OUTPUT:
// Electron density @ (R,Z), 1/(cc GeV) 
// -calculated for a uniform distribution of sources (within R<Rg, |Z|<hg), 
// the source normalization is 1/(cc GeV) at 1 GeV; 
// the energy losses dE/dt= -elossconst*E^2
// REFERENCE: Bulanov, S.V., Dogiel, V.A. 1974, Astrophys. Spa. Sci. 29,305,
//            their eq.(8) has to be divided by 2pi (-perhaps an error)
//                                                  I.V.Moskalenko  11/01/2006
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double eprop(double Ee, double R, double Z, 
	     double Rg, double hg, double hh, double elossconst,
	     double gamma_e, double Dxx, double g1)
{
  static int readarray=0;
  static double J0zero[500];
  int i, j, m=100, n=200;
  double x,S1,S2,S3,Pi,eps=1.e-16;

  Pi=acos(-1.);
  if(readarray==0)
    {
      ifstream data;
      data.open("j0zero500.dat");                    // open file if exists
      if(data.fail())
	{
	  cerr<<" >>eprop>> Error opening file "<<"j0zero500.dat"<<endl;
	  exit(1);
	}
      for(i=0; i<500; data >> J0zero[i++]);
      data.close();
      readarray=1;
    }

  for(S2=1.,S3=0., i=0; i<m || S2>eps*S3; i++)
    {
      for(S1=0., j=0;j<n;j++)
	{
	  x=-(pow(Pi*(i+0.5),2)+pow(J0zero[j]*hh/Rg,2))
	    *Dxx*pow(Ee,g1-1)/(hh*hh*elossconst)/(1.-g1);
	  S1+= 
	    bessj(0,J0zero[j]*R/Rg) / bessj(1,J0zero[j]) / J0zero[j]*
	    hyp1F1(1.,(gamma_e-g1)/(1.-g1),x);	  	  
	}
      S2=S1*sin(Pi*hg/hh*(i+0.5))*cos(Pi*Z/hh*(i+0.5))/(i+0.5);
      S3+=S2;
    }
  S3*=4.*pow(Ee,-gamma_e-1)/(Pi*(gamma_e-1.)*elossconst);
  return S3*1000.; // units 1/(cm3 GeV)
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// Bessel J0, J1 functions: adaptations of CERNLIB C312 fortran routines 
//                                                  I.V.Moskalenko  11/01/2006
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double bessj(int order, double x)
{
  double A,B,F,P,Q,V;

  switch(order)
    {
    case 0: // bessj0
      V=fabs(x);
      if(V<8.)
	{
	  F=0.0625*x*x-2.;
	  A =           - 0.0000000000000008;
	  B = F * A     + 0.0000000000000413;
	  A = F * B - A - 0.0000000000019438;
	  B = F * A - B + 0.0000000000784870;
	  A = F * B - A - 0.0000000026792535;
	  B = F * A - B + 0.0000000760816359;
	  A = F * B - A - 0.0000017619469078;
	  B = F * A - B + 0.0000324603288210;
	  A = F * B - A - 0.0004606261662063;
	  B = F * A - B + 0.0048191800694676;
	  A = F * B - A - 0.0348937694114089;
	  B = F * A - B + 0.1580671023320973;
	  A = F * B - A - 0.3700949938726498;
	  B = F * A - B + 0.2651786132033368;
	  A = F * B - A - 0.0087234423528522;
	  A = F * A - B + 0.3154559429497802;
	  return 0.5*(A-B);
	}
      else
	{
	  F=256./(x*x)-2.;
	  B =           + 0.0000000000000007;
	  A = F * B     - 0.0000000000000051;
	  B = F * A - B + 0.0000000000000433;
	  A = F * B - A - 0.0000000000004305;
	  B = F * A - B + 0.0000000000051683;
	  A = F * B - A - 0.0000000000786409;
	  B = F * A - B + 0.0000000016306465;
	  A = F * B - A - 0.0000000517059454;
	  B = F * A - B + 0.0000030751847875;
	  A = F * B - A - 0.0005365220468132;
	  A = F * A - B + 1.9989206986950373;
	  P=A-B;
	  B =           - 0.0000000000000006;
	  A = F * B     + 0.0000000000000043;
	  B = F * A - B - 0.0000000000000334;
	  A = F * B - A + 0.0000000000003006;
	  B = F * A - B - 0.0000000000032067;
	  A = F * B - A + 0.0000000000422012;
	  B = F * A - B - 0.0000000007271916;
	  A = F * B - A + 0.0000000179724572;
	  B = F * A - B - 0.0000007414498411;
	  A = F * B - A + 0.0000683851994261;
	  A = F * A - B - 0.0311117092106740;
	  Q=8.0*(A-B)/V;
	  F=V-0.785398163397448;
	  A=cos(F);
	  B=sin(F);
	  F=0.398942280401432/sqrt(V);
	  return F*(P*A-Q*B);
	}

    case 1: // bessj1
      V=fabs(x);
      if(V<8.)
	{
	  F=0.0625*x*x-2.;
	  B =           + 0.0000000000000114;
	  A = F * B     - 0.0000000000005777;
	  B = F * A - B + 0.0000000000252812;
	  A = F * B - A - 0.0000000009424213;
	  B = F * A - B + 0.0000000294970701;
	  A = F * B - A - 0.0000007617587805;
	  B = F * A - B + 0.0000158870192399;
	  A = F * B - A - 0.0002604443893486;
	  B = F * A - B + 0.0032402701826839;
	  A = F * B - A - 0.0291755248061542;
	  B = F * A - B + 0.1777091172397283;
	  A = F * B - A - 0.6614439341345433;
	  B = F * A - B + 1.2879940988576776;
	  A = F * B - A - 1.1918011605412169;
	  A = F * A - B + 1.2967175412105298;
	  return 0.0625*(A-B)*x;
	}
      else
	{
	  F=256./(x*x)-2.;
	  B =           - 0.0000000000000007;
	  A = F * B     + 0.0000000000000055;
	  B = F * A - B - 0.0000000000000468;
	  A = F * B - A + 0.0000000000004699;
	  B = F * A - B - 0.0000000000057049;
	  A = F * B - A + 0.0000000000881690;
	  B = F * A - B - 0.0000000018718907;
	  A = F * B - A + 0.0000000617763396;
	  B = F * A - B - 0.0000039872843005;
	  A = F * B - A + 0.0008989898330859;
	  A = F * A - B + 2.0018060817200274;
	  P=A-B;
	  B =           + 0.0000000000000007;
	  A = F * B     - 0.0000000000000046;
	  B = F * A - B + 0.0000000000000360;
	  A = F * B - A - 0.0000000000003264;
	  B = F * A - B + 0.0000000000035152;
	  A = F * B - A - 0.0000000000468636;
	  B = F * A - B + 0.0000000008229193;
	  A = F * B - A - 0.0000000209597814;
	  B = F * A - B + 0.0000009138615258;
	  A = F * B - A - 0.0000962772354916;
	  A = F * A - B + 0.0935555741390707;
	  Q=8.*(A-B)/V;
	  F=V-2.356194490192345;
	  A=cos(F);
	  B=sin(F);
	  F=0.398942280401432/sqrt(V);
	  return( x>0. ? F*(P*A-Q*B): -F*(P*A-Q*B));
	}
    default:
      cout<<" >>bessj>> routine has been called with order="<<order<<endl
	  <<" >>bessj>> routine calculates Bessel function of orders 0 and 1"<<endl;
      exit(1);
    }
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// Hypergeometrical function 1F1(a,b,z)
//                                                  I.V.Moskalenko  11/01/2006
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double a, b, z;
double hyp1F1(double a1, double b1, double z1)
{
  double ga, gb, gba;
  a=a1; b=b1; z=z1;

  if (a  > 0.) ga = (a  <1. ? log(gamma(a  +1.))-log(a)   : log(gamma(a)));
  if (b  > 0.) gb = (b  <1. ? log(gamma(b  +1.))-log(b)   : log(gamma(b)));
  if (b-a> 0.) gba= (b-a<1. ? log(gamma(b-a+1.))-log(b-a) : log(gamma(b-a)));

  return( exp(gb-ga-gba)*sim(0.,1.,1.e-4,1.e-8,1.e-20,&fu) );
}

double fu(double t)
{
    return(exp(z*t)*pow(t,a-1.)*pow(1.-t,b-a-1.));
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// Gamma function: adaptation of CERNLIB C305 fortran routine 
//                                                  I.V.Moskalenko  11/01/2006
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

/*double gamma(const double x)
{
  double G,F,Z;
  double C[13]={
    0.000539698958808, 0.002619307282746, 0.020449630823590,
    0.073094836414370, 0.279643691578538, 0.553387692385769,
    0.999999999999998,-0.000832724708684, 0.004698658079622,
    0.022523834747260,-0.170447932874746,-0.056810335086194,
    1.130603357286556};
  
  if(x==0.) return 1.;
  if(x == -(int) fabs(x))
    {
      cout<<" >>gamma>> argument is non-positive integer = "<<x<<endl;
      exit(1);
    }
  Z=x;
  if(x<0.) Z=1.-Z;
  F=1./Z;
  if(Z>1.)
    {
      for(F=1.;Z>=2.;F*=--Z);
      Z--;
    }
  G= F*((((((C[0]*Z+C[1])*Z+C[2])*Z+C[3])*Z+C[4])*Z+C[5])*Z+C[6])/
    ((((((C[7]*Z+C[8])*Z+C[9])*Z+C[10])*Z+C[11])*Z+C[12])*Z+1.);
  if(x>0.) return G;
  return 3.141592653589793/(sin(3.141592653589793*x)*G);
}
*/
/***********************************************************************
c ### The routine has been rewritten from FORTRAN source code ###
c calculation the definite integral by Simpson's method with the automatic
c choice of the integration step 
C INPUT: A1,B1 - the limits of integration; H1 - the initial step;
C REPS,AEPS - the relative and absolute precision; FU - the name of the 
C user-defined function f(x); OUTPUT: sim - the value of the integral;
c other values that are calculated in parallel:
C AIH - the value of integral with one more step of integration;
C AIABS - the value of the integral for module of the integrand;
C # NOTE # the subprogram returns the value of integral as one of the
C precise conditions (AEPS1,EPS1) are reached; when AEPS1=EPS1=0, 
c then it is calculated with the constant step H1. See appended test case.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
***********************************************************************/

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
/*
double sim(double A, double B, double H, double REPS, double AEPS, double (*fu)(double))
{
  int K;
  double F[8],P[6],S,C,X,X0,AI,AIH,AIABS,DI1,DI2,DI3,EPS,DELTA;
  H=sign(H,B-A);
  S=sign(1.,H);
  AI=AIH=AIABS=0.;
  P[2]=P[4]=4.;
  P[3]=2.;
  P[5]=1.;
  if(B-A==0.) return(AI);
  REPS=abs(REPS);
  AEPS=abs(AEPS);
  for(K=1;K<8;F[K++]=1.e20);
  X=A;
  C=0.;
  F[1]=fu(X)/3.;

 L4:
  X0=X;
  if((X0+4.*H-B)*S>0.)
    {
      H=(B-X0)/4.;
      if(H==0.) return(AI);
      for(K=2;K<8;F[K++]=1.e20);
      C=1.;
    }

 L5:
  DI2=F[1];
  DI3=abs(F[1]);
  for(K=2;K<6;K++)
    {
      X+=H;
      if((X-B)*S>=0.) X=B;
      if(F[K]-1.e20==0.) F[K]=fu(X)/3.;
      DI2+=P[K]*F[K];
      DI3+=P[K]*abs(F[K]);
    }
  DI1=(F[1]+4.*F[3]+F[5])*2.*H;
  DI2*=H;
  DI3*=H;
  if(REPS==0.&& AEPS==0.) goto L14;
  EPS=abs((AIABS+DI3)*REPS);
  if(EPS-AEPS<0) EPS=AEPS;
  DELTA=abs(DI2-DI1);
  if(DELTA-EPS<0.) { if(DELTA-EPS/8.>=0.) goto L14; }
  else goto L21;
  H*=2.;
  F[1]=F[5];
  F[2]=F[6];
  F[3]=F[7];
  for(K=4;K<8;F[K++]=1.e20);
  goto L18;

 L14:  
  F[1]=F[5];
  F[3]=F[6];
  F[5]=F[7];
  F[2]=F[4]=F[6]=F[7]=1.e20;

 L18:  
  DI1=DI2+(DI2-DI1)/15.;
  AI+=DI1;
  AIH+=DI2;
  AIABS+=DI3;
  goto L22;

 L21:
  H/=2.;
  F[7]=F[5];
  F[6]=F[4];
  F[5]=F[3];
  F[3]=F[2];
  F[2]=F[4]=1.e20;
  X=X0;
  C=0.;
  goto L5;

 L22:
  if(C==0) goto L4;
  return(AI);
}
*/
/*
int main()
#include <iomanip>
{
  double 
    Rg=30.,      //kpc, R Galaxy
    hg=0.2,      //kpc, disk (with the sources) half-thickness 
    hh=4.,       //kpc, halo size
    Ugevcc=1.e-9,//GeV/cc, photon field energy density (IC energy losses)
    Ke=1.,       //electron spectrum normalization at 1 GeV
    gamma_e=2.4, //electron spectrum injection index
    Dxx=1.e28,   //cm2/s, normalization of the diffusion coefficient @ 1 GeV
    g1=0.6,      //diffusion coefficient energy-dependence index
    kpc=3.085677e21; //cm, =1 kpc
  double 
    Rele=2.8179409238e-13,     // cm, =e^2/mc^2 class. electron radius
    C=2.99792458e10,           // cm/s, =c speed of light;
    Mele = 0.5109990615e-3;    // GeV/c^2, electron rest mass

  double Pi=acos(-1.);
  double elossconst=32./9.*Pi*pow(Rele/Mele,2)*C*Ugevcc; // 1/(GeV s), energy loss const
  Dxx*=pow(kpc,-2);
  double Ee, R=8.5, Z=0.01;

// output = spectrum: Ee [MeV], (c/4pi)*Flux [1/(cm2 s sr MeV)]
  for(Ee=0.01;Ee<1.e4;Ee*=2.)
    cout<<Ee *1.e3<<" "<<Ke*C/(4.*Pi)/1000. //Flux electrons/(cm2 s sr MeV); 1/GeV->1/MeV
      *eprop(Ee, R, Z, Rg, hg, hh, elossconst, gamma_e, Dxx, g1)
      <<" "<<Ke*C/(4.*Pi) *hh*hg*pow(Ee,-gamma_e-g1)/Dxx*(1.-Z/hh)//1D solution w/o losses
      <<endl;
// output = Z-profile
  for(Ee=1., Z=0.;Z<=hh;Z+=0.1)
    cout<<Z<<" "<<Ke*C/(4.*Pi)/1000. //Flux electrons/(cm2 s sr MeV); 1/GeV->1/MeV
      *eprop(Ee, R, Z, Rg, hg, hh, elossconst, gamma_e, Dxx, g1)
      <<" "<<Ke*C/(4.*Pi) *hh*hg*pow(Ee,-gamma_e-g1)/Dxx*(1.-Z/hh)//1D solution w/o losses
      <<endl;

  exit(0);
//*************** TESTS of special functions *******************
  double x,b,z;

  for (x=-10.;x<11.;x+=2.)  // Bessel functions
    cout<<x<<" "<<" "<<bessj(0,x)<<" "<<bessj(1,x)<<endl;

  for(z=-1.,b=3;b<10;b++)   // hypergeometrical function
    cout<<hyp1F1(1.,b,z)<<endl;

  for (x=1.;x<1.2;x+=0.005) // gamma-function
    cout<<x<<" "<<setprecision(17)<<" "<<gamma(x)<<endl;
}

*/

