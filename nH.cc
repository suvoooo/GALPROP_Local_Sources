
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * nH.cc *                                       galprop package * 2/20/2001 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// The formalism is described in Moskalenko I.V. et al. 2002, ApJ 565, 280
//
//1 The routine nH2 calculates H2 number density in mol/cm3
//1 nH2(R,Z) = epsilon0(R) *X /kpc2cm *exp(-ln2 (Z-Z0)^2/Zh^2), 
//1 where epsilon0(R) (K km s^-1 kpc^-1) -CO volume emissivity, and
//1 Z0(R), Zh(R) are taken from [B88]/Table 3/Cols 4,7,10;
//1 Additional data from W90 and F07.  Only average radial dependence of the
//1 surface density of F07 model is used.  The scale height is assumed to be
//1 70 pc and the X factor 5e19.  Only innermost 1.5 kpc is used.
//1 X = nH2/nCO =1.9e20 is the conversion factor taken from [SM96].
//1 [B88] Bronfman L. et al. 1988, ApJ 324, 248
//1 [W90] Wouterloot J. G. A. et al. 1990, A&A 230, 21
//1 [F07] Ferriere K. et al.  2007, A&A 467, 611
//1 [SM96] Strong & Mattox 1996, A&A 308, L21
//
//2 The routine nHI calculates HI number density in cm^-3
//2 nHI(R,Z) = Y(R)/nGB *fZ(Z), 
//2 where Y(R) (cm^-3) - nHI number density taken from [GB76],
//2 nGB = 0.33 cm^-3 their average density of the disk at 4-8 kpc,
//2 fZ(Z) - Z-distribution: 
//2 0<R<8 [DL90], R>10 kpc [C86], with linear interpolation in 8-10 kpc.
//2 The total distribution is thus renormalized to the average surface density
//2 =6.2e20 cm^-2 at 4-8 kpc [DL90].
//2 References:
//2 [C86]  Cox, Kruegel, Mezger 1986, A&A 155, 380
//2 [DL90] Dickey & Lockman 1990, Ann. Rev. Astron. Astrophys. 28, 215
//2 [GB76] Gordon & Burton 1976, ApJ 208, 346/Table 1
//
//3 The routine nHII calculates HII number density in cm^-3
//3 Cordes J.M., et al. 1991, Nature 354, 121 / Equation (6) and Table 1
//3 Annular Gaussian preferred model for narrow component
//3 nHII model can now be chosen:
//3 nHII_model=1   original Cordes et al 1991.
//3 nHII_model=2 : NE2001 from Cordes and Fazio 2002 , 	arXiv:astro-ph/0207156v3
//3 nHII_model=3 : Gaensler et al. 2009 Publications of the Astronomical Society of Australia 25(4), 184 : wide component.  NE2001: narrow component.
//
//4 nH*_av routines make averaging over smaller steps.
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

using namespace std;//AWS20050624
#include<cmath>
#include<iostream>  //AWS20090813
#include"constants.h"
#include "galprop_internal.h"
#include <cstdlib>
#include <sstream>

#include "ErrorLogger.h"

//AWS20090814
// Parameters to allow choice of models
// Implemented via a namespace to keep them addressed only in this code unit and avoid changing the calling routines.
// Defined outside functions to preserve between function calls (i.e. they are global, but namespace avoids any conflict).
// Initialized to original models so will be upward compatible with any previous use in inside or outside galprop.
// Only nHII is affected in this version.


namespace nH_models {int nHI_model=1, nH2_model=1, nHII_model=1; int debug=0; } // intialize to original models
using namespace nH_models;

//////////////////////////////////////////////////////////////////////////////////////
int nH_set_model(int nHI_model_, int nH2_model_, int nHII_model_, int debug_)
{
  nHI_model = nHI_model_;
  nH2_model = nH2_model_;
  nHII_model= nHII_model_;
  debug=debug_;

  ostringstream buf;
  buf<<"nHI_model: "<<nHI_model <<" nH2_model: " <<nH2_model<<" nHII_model:" <<nHII_model<<" debug: "<<debug;
  INFO(buf.str());
  return 0;
}
//////////////////////////////////////////////////////////////////////////////////////

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double nH2(double Rkpc, double Zkpc)
{
   double H2toCO = 1.9e20;
   return nH2(Rkpc,Zkpc,H2toCO);
   /*
   int i;
   double nH2_ = 0.0, fR,fZ0,fZh;                                              // [B88]/Table 3
   double R[18] ={ 0.00, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75,       // kpc, col.1
                         6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75,10.25},
          Y[18] ={ 0.00,  1.5,  3.3,  5.8,  5.5,  8.4,  9.0,  9.6,  8.6,       // CO, K km s^-1
		          9.1,  7.9,  9.2,  7.7,  5.0,  3.6,  4.8,  1.7,  0.0},// (col.4)
	  Z0[18]={0.039,0.039,0.036,0.000,-.008,0.001,-.010,-.001,-.004,       // kpc, col.7
                        -.019,-.022,-.014,-.009,-.004,0.013,-.004,-.020,-.020},
	  Zh[18]={0.077,0.077,0.080,0.061,0.065,0.071,0.072,0.082,0.083,       // kpc, col.10
                        0.073,0.063,0.058,0.072,0.080,0.066,0.023,0.147,0.147};
   double H2toCO = 1.9e20;              // [SM96]

   if(Rkpc > R[17]) return nH2_;
   for (i=0; i<17; i++)  if(R[i] <= Rkpc && Rkpc <= R[i+1])  break;

   fR =  Y[i] + ( Y[i+1] - Y[i])/(R[i+1] - R[i])*(Rkpc - R[i]); 
   fZ0= Z0[i] + (Z0[i+1] -Z0[i])/(R[i+1] - R[i])*(Rkpc - R[i]);
   fZh= Zh[i] + (Zh[i+1] -Zh[i])/(R[i+1] - R[i])*(Rkpc - R[i]);
   nH2_ =  fR * exp( -log(2.)*pow( (Zkpc-fZ0)/fZh, 2 ) )  *H2toCO/kpc2cm;
   return nH2_< 0. ? 0.: nH2_;
   */
}
/////////////////////////////////////////////////////////////////////
// version with H2toCO input
/////////////////////////////////////////////////////////////////////
double nH2(double Rkpc, double Zkpc,double H2toCO) // AWS20090616
{
   int i;
   double nH2_ = 0.0, fR,fZ0,fZh;                                              // [B88]/Table 3
   double R[44] ={ 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80,
                         0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50,             // F07 model
                         2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75,       // kpc, col.1
                         6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75,
                        10   ,11   ,12   ,13   ,14   ,15   ,16   ,17   ,
                        18   ,19   ,20   ,21   },                                 //Add data from [W90] Gulli 20100218
          Y[44] ={ 43.7, 24.5, 10.7,  1.6,  1.4,  1.5,  1.1,  0.9, 0.8,
	                  0.7,  0.6,  0.5,  0.4,  0.3,  0.2,  0.1,             // F07 model converted into average surface density and the CO using average scale heighe
	                  1.5,  3.3,  5.8,  5.5,  8.4,  9.0,  9.6,  8.6,       // CO, K km s^-1 
		          9.1,  7.9,  9.2,  7.7,  5.0,  3.6,  4.8,  1.7,
	                  0.4,  1.2,  1.0,  0.7,  0.3, 0.15, 0.09, 0.07,
	                 0.05,0.005,0.008,0.004},// (col.4)
	  Z0[44]={0.000,    0,    0,    0,    0,    0,    0,    0,    0,       //Ferriere's model is symmetric around z (on average at least)
	                    0,    0,    0,    0,    0,    0,    0,
	                0.039,0.036,0.000,-.008,0.001,-.010,-.001,-.004,       // kpc, col.7
                        -.019,-.022,-.014,-.009,-.004,0.013,-.004,-.020,
	                0    ,    0,    0,    0,    0,    0,    0,    0,
	                0    ,    0,    0,    0},
	  Zh[44]={0.070,0.070,0.070,0.070,0.070,0.070,0.070,0.070,0.070,
	                0.070,0.070,0.070,0.070,0.070,0.070,0.070,             //The scale height should be handled better, but the Ferriere model does not fit nicely into this functional form.  So this fixed approximation is made
	                0.077,0.080,0.061,0.065,0.071,0.072,0.082,0.083,       // kpc, col.10
                        0.073,0.063,0.058,0.072,0.080,0.066,0.023,0.147,
	                0.111,0.136,0.147,0.160,0.223,0.257,0.220,0.200,
	                0.188,0.200,0.200,0.200};
	  //   double H2toCO = 1.9e20;              // [SM96]

   if(Rkpc > R[43]) return nH2_;
   for (i=0; i<43; i++)  if(R[i] <= Rkpc && Rkpc <= R[i+1])  break;

   fR =  Y[i] + ( Y[i+1] - Y[i])/(R[i+1] - R[i])*(Rkpc - R[i]); 
   fZ0= Z0[i] + (Z0[i+1] -Z0[i])/(R[i+1] - R[i])*(Rkpc - R[i]);
   fZh= Zh[i] + (Zh[i+1] -Zh[i])/(R[i+1] - R[i])*(Rkpc - R[i]);
   nH2_ =  fR * exp( -log(2.)*pow( (Zkpc-fZ0)/fZh, 2 ) )  *H2toCO/kpc2cm;
   return nH2_< 0. ? 0.: nH2_;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double nHI(double Rkpc, double Zkpc)
{
   int i;                                                             // Table 1 [GB76]
   double R[30] ={ 0.0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0,  // kpc, col.1
                   6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,10.0,10.5,11.0,
                  11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,15.5,16.0},
          Y[30] ={ .10, .13, .14, .16, .19, .25, .30, .33, .32, .31,  // nHI, cm^-3
		   .30, .37, .38, .36, .32, .29, .38, .40, .25, .23,  // (col.3)
                   .32, .36, .32, .25, .16, .10, .09, .08, .06, .00};
   double fR, fZ,fZ1=0.,fZ2=0., R1,R2=R[29], Y1,Y2=Y[29];
   double nGB =0.33, nDL =0.57;       // cm^-3, disk density @ 4-8 kpc; [GB76], [DL90]
   double A1=0.395,     z1=0.212/2.,  // cm^-3, kpc; Z-distribution parameters from [DL90]
          A2=0.107,     z2=0.530/2.,
          B =0.064,     zh=0.403;

   for (i=0; i<28; i++)  if(R[i] <= Rkpc && Rkpc <= R[i+1])  break;  //Gulli20070810  i=28 if condition never met 

   R1 = (R[i]+R[i+1])/2;   Y1 = Y[i];
   if(Rkpc < R1)
   {  
      if(i> 0)    { R2 = (R[i-1]+R[i])/2;   Y2 = Y[i-1]; }
      else        { R2 = R[0];              Y2 = Y[0];   }
   }
   else  if(i<28) { R2 = (R[i+1]+R[i+2])/2; Y2 = Y[i+1]; }

   fR = Y1 +(Y2 -Y1)/(R2 -R1)*(Rkpc -R1);                             // interpolation in R

   R2 = (R[28] +R[29]) /2;
   if(Rkpc > R2) fR = Y[28]*exp(-(Rkpc-R2)/3);                        // extrapolation in R

// calculation of Z-dependence
   if(Rkpc <10.)                                                      // [DL90]
      fZ1 =A1*exp(-log(2.)*pow(Zkpc/z1,2))+A2*exp(-log(2.)*pow(Zkpc/z2,2))+B*exp(-fabs(Zkpc)/zh);
   if(Rkpc > 8.) fZ2=nDL*exp(-pow(Zkpc /(0.0523*exp(0.11*Rkpc)), 2)); // [C86] IMOS20010220

   if(Rkpc <= 8.) fZ = fZ1;
   else
   {   if(Rkpc >=10.) fZ = fZ2;
       else fZ = fZ1 +(fZ2 -fZ1)/2.*(Rkpc -8.);                       // interp. [DL90] & [C86]
   }
   return fZ *fR/nGB;
}

// original nHII routine, now replaced by multiple choice
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
/*
double nHII(double Rkpc,double Zkpc)
{
   double fne1=0.025, H1=1.00, A1=20.0;
   double fne2=0.200, H2=0.15, A2= 2.0;
   double R2=4.0;
   double ne1 = fne1 * exp(-fabs(Zkpc)/H1) * exp (-pow( Rkpc    /A1, 2));
   double ne2 = fne2 * exp(-fabs(Zkpc)/H2) * exp (-pow((Rkpc-R2)/A2, 2));
   return ne1+ne2;
}
*/
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double nHII(double Rkpc,double Zkpc) //AWS20090814
{
  // version with model choice 
  //  nH_set_model(1,1,3,0); just to test this function
 
  int model_found = 0;

  if(nHII_model==1) 
  { 
   // original galprop model from Cordes et al. 1991
   double fne1=0.025, H1=1.00, A1=20.0;
   double fne2=0.200, H2=0.15, A2= 2.0;
   double R2=4.0;
   double ne1 = fne1 * exp(-fabs(Zkpc)/H1) * exp (-pow( Rkpc    /A1, 2));
   double ne2 = fne2 * exp(-fabs(Zkpc)/H2) * exp (-pow((Rkpc-R2)/A2, 2));

  if (debug==1)cout<<"nHII: HII_model="<<nHII_model <<" Rkpc="<<Rkpc <<" Zkpc="<<Zkpc   <<" ne1="<<ne1<<" ne2="<<ne2 <<endl  ;
  model_found=1;

   return ne1+ne2;
  }

  if(nHII_model==2) 
  {
   // NE2001 from Cordes and Fazio 2002 , 	arXiv:astro-ph/0207156v3
   // a simplified version using only the 'smooth' components
   // Table 2 for formulae and Table 3 right column for values

    // wide component
    double H1 = 0.95 ;    
    double n1 = 0.033/H1 ; // since n1*H1 is given
    double A1 = 17.;

    // narrow component
    double H2 = 0.14;
    double n2 = 0.09;
    double A2 = 3.7;

    double Rsun = 8.5;

    double pi=acos(-1.0);

    // wide component
    double  g1  = cos(pi*Rkpc/(2.*A1)) / cos (pi*Rsun/(2.*A1));
    if(Rkpc >= A1) g1 = 0.0; // step function
    double  h1  = pow(cosh(Zkpc/H1),-2.0);  // sech^2
    double ne1  = n1 * g1 * h1;

    // narrow component
    // NB the formula for g2 in the NE2001 is wrong, does not agree with the code or the plots in the paper. It uses A2^2 in the denominator while the correct value is 1.8^2 
    //    it is given correctly in Taylor and Cordes 1993 ApJ 411, 674, and NE2001 uses the same function for g2.
    double  g2  = exp( -(pow(Rkpc-A2 ,2.))  / pow(1.8,2.) );
    double  h2 =  pow(cosh(Zkpc/H2),-2.0);  // sech^2
    double ne2 =  n2  * g2 * h2;
    
    if (debug==1)cout<<"nHII: HII_model="<<nHII_model <<" Rkpc="<<Rkpc <<" Zkpc="<<Zkpc   <<" ne1="<<ne1<<" ne2="<<ne2 <<endl  ;
    model_found=1;

    return ne1 + ne2;
  }

  if(nHII_model==3) 
  {
    // Gaensler et al. 2009  Publications of the Astronomical Society of Australia 25(4), 184
    // only z-dependence in this paper
    // use same notation as NE2001 for clarity

    double H1 = 1.8          ; //exponential scaleheight in kpc
    double n1 = 0.014        ; // density at z=0

    double ne1 = n1*exp(-fabs(Zkpc/H1));
    

   // narrow component from NE2001 (since not in Gaensler 2008 model).
    double H2 = 0.14;
    double n2 = 0.09;
    double A2 = 3.7;

   
    double  g2  = exp( -(pow(Rkpc-A2 ,2.))  / pow(1.8,2.) ); // see nHII==2 case for explanation
    double  h2 =  pow(cosh(Zkpc/H2),-2.0);  // sech^2
    double ne2 =  n2  * g2 * h2;

   if (debug==1)cout<<"nHII: HII_model="<<nHII_model <<" Rkpc="<<Rkpc <<" Zkpc="<<Zkpc   <<" ne1="<<ne1<<" ne2="<<ne2 <<endl  ;
   model_found=1;

   return ne1 + ne2;
  }

  if(model_found==0){ cout<<"nHII: invalid model:"<<nHII_model<<" exiting !  "<<endl; exit(0);}

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double nH_av(double x, double y, double z, double dz, double dzz, double(*nH)(double,double)) //IMOS20080114
{
  double R=sqrt(x*x + y*y), nH_av_=0;
  int nuse=0;

  for(double zz=z-dz/2.+dzz/2.; zz<z+dz/2.; zz+=dzz)
    {
      nH_av_+=nH(R,zz);
      nuse++;
    }
  return nH_av_/nuse;
}
/////////////////////////////////////////////////////////////////////
// version with H2toCO input
/////////////////////////////////////////////////////////////////////
double nH_av(double x, double y, double z, double dz, double dzz, double H2toCO,double(*nH)(double,double,double)) //AWS20090616
{
  double R=sqrt(x*x + y*y), nH_av_=0;
  int nuse=0;

  for(double zz=z-dz/2.+dzz/2.; zz<z+dz/2.; zz+=dzz)
    {
      nH_av_+=nH(R,zz,H2toCO);
      nuse++;
    }
  return nH_av_/nuse;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// (l,b) = galactic coordinates; db - size of b-grid; ddb -finer step to average
// d - the distance; dd - line-of-sight integration step; ddd -finer step to average
// make sure that ddb = db/N; ddd = dd/N, where N=integer (N~5-10) 
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double nH_av_lb(double l, double b, double db, double ddb, double d, double dd, double ddd, double(*nH)(double,double)) //IMOS20080114
{
  double dtr=acos(-1.)/180.; // conversion degrees to radians
  double sinl=sin(l*dtr);
  double cosl=cos(l*dtr);
  double nH_av_lb_=0;
  int nuse=0;
  
  for(double bb=b-db/2.+ddb/2.; bb<b+db/2.; bb+=ddb)
    {
      if (abs(bb) > 90.) continue;
      double sinb=sin(bb*dtr);
      double cosb=cos(bb*dtr);

      for(double d1=d-dd/2.+ddd/2.; d1<d+dd/2.; d1+=ddd)
	{
	  double z=d1*sinb;                                              // altitude of the current point 
	  double R=sqrt(Rsun*Rsun+pow(d1*cosb,2)-2.0*Rsun*d1*cosb*cosl); // Galactocentric distance of the current point
	  
	  nH_av_lb_+=nH(R,z);
	  nuse++;
	}
    }
  return nH_av_lb_/nuse;
}
/////////////////////////////////////////////////////////////////////
// version with H2toCO input
/////////////////////////////////////////////////////////////////////
double nH_av_lb(double l, double b, double db, double ddb, double d, double dd, double ddd, double H2toCO,double(*nH)(double,double,double)) //AWS20090622
{
  double dtr=acos(-1.)/180.; // conversion degrees to radians
  double sinl=sin(l*dtr);
  double cosl=cos(l*dtr);
  double nH_av_lb_=0;
  int nuse=0;
  
  for(double bb=b-db/2.+ddb/2.; bb<b+db/2.; bb+=ddb)
    {
      if (abs(bb) > 90.) continue;
      double sinb=sin(bb*dtr);
      double cosb=cos(bb*dtr);

      for(double d1=d-dd/2.+ddd/2.; d1<d+dd/2.; d1+=ddd)
	{
	  double z=d1*sinb;                                              // altitude of the current point 
	  double R=sqrt(Rsun*Rsun+pow(d1*cosb,2)-2.0*Rsun*d1*cosb*cosl); // Galactocentric distance of the current point
	  
	  nH_av_lb_+=nH(R,z,H2toCO);
	  nuse++;
	}
    }
  return nH_av_lb_/nuse;
}
/* //IMOS20080114

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double nH2_av(double x,double y,double z,double dz,double dzz)
{  
   double R=sqrt(x*x + y*y);
   double nH2_av_=0.0;
   int nuse=0;

   for (double zz=z-dz/2.; zz<=z+dz/2.; zz+=dzz)
   {
      nH2_av_+=nH2(R,zz);
      nuse++;
   }
   return nH2_av_/nuse;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double nHI_av(double x,double y,double z,double dz,double dzz)
{
  
   double R=sqrt(x*x + y*y);
   double nHI_av_=0.0;
   int nuse=0;

   for (double zz=z-dz/2.; zz<=z+dz/2.; zz+=dzz)
   {
      nHI_av_+=nHI(R,zz);
      nuse++;
   }
   return nHI_av_/nuse;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double nHII_av(double x,double y,double z,double dz,double dzz)
{  
   double R=sqrt(x*x + y*y);
   double nHII_av_=0.0;
   int nuse=0;

   for(double zz=z-dz/2.; zz<=z+dz/2.; zz+=dzz)
   {
      nHII_av_+=nHII(R,zz);
      nuse++;
   }
   return nHII_av_/nuse;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
*/

/*
#include<stdio.h>
main()
{
   for(double R=0.; R<20.; R+=0.05) printf("%f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n",
     R,2*nH2(R,0.),2*nH2(R,0.1),2*nH2(R,0.2),nHI(R,0.),nHI(R,0.1),nHI(R,0.2),nHII(R,0.),nHII(R,0.1),nHII(R,0.2));

   for(double R=0.; R<20.; R+=1) printf("%f  %f  %f  %f  %f\n",
     R,nH2_av(R,0.,0.1,0.1,0.01),nH_av(R,0.,0.1,0.1,0.01,nH2),nHI_av(R,0.,0.1,0.1,0.01),nH_av(R,0.,0.1,0.1,0.01,nHI));
}
*/
