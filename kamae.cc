//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * kamae.cc *                                    galprop package * 2006/11/07
//**.****!****.****!****.****!****.****!****.****!****.****!****.****!****.****!
//
// Implementation of parameterization of Kamae et al., ApJ, 647, 692
//
//**.****!****.****!****.****!****.****!****.****!****.****!****.****!****.****!

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <iostream>
#include <cmath>
#include <valarray>
#include "constants.h"
//#define Mp 0.938
#define max(a,b) (((a) > (b)) ? (a) : (b))
using namespace std;

// functions for component parameter calculation
valarray<double> kamae_gamma_param_nd(double);
valarray<double> kamae_gamma_param_diff(double);
valarray<double> kamae_gamma_param_delta(double);
valarray<double> kamae_gamma_param_res(double);
valarray<double> kamae_elec_param_nd(double);
valarray<double> kamae_elec_param_diff(double);
valarray<double> kamae_elec_param_delta(double);
valarray<double> kamae_elec_param_res(double);
valarray<double> kamae_posi_param_nd(double);
valarray<double> kamae_posi_param_diff(double);
valarray<double> kamae_posi_param_delta(double);
valarray<double> kamae_posi_param_res(double);
// functions for component differential cross section
double kamae_nd(double, double, valarray<double>, int);
double param_diff(double, double, valarray<double>);
double param_delta(double, double, valarray<double>);
double kamae_res(double, double, valarray<double>);
// functions for total inelastic differential cross sections
double kamae(double, double, int, int, int, valarray<double>, valarray<double>, valarray<double>, valarray<double>);

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//		calculate non-diff parameters for gamma-rays
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
valarray<double> kamae_gamma_param_nd(double Tp) {

  valarray<double> param_a(0., 9);

  double y, z;
  int i;
  
  y = log10(Tp*0.001);
  
  if ((Tp > 0.487) && (Tp < 512000.1)) {

    z = y + 3.3;
    param_a[0] = z*(-0.51187 + z*(7.6179 + z*(-2.1332 + 0.22184*z)));
    param_a[1] = -1.2592e-5 + 1.4439e-5*exp(-0.29360*(y + 3.4)) + 5.9363e-5/(y + 4.1485) + y*(2.2640e-6 - 3.3723e-7*y);
    param_a[2] = -174.83 + 152.78*log10(1.5682*(y + 3.4)) - 808.74/(y + 4.6157);
    param_a[3] = 0.81177 + y*(0.56385 + y*(0.0040031 + y*(-0.0057658 + 0.00012057*y)));
    z = y + 3.32;
    param_a[4] = z*(0.68631 + z*(10.145 + z*(-4.6176 + z*(0.86824 - 0.053741*z))));
    z = y + 4.7171;
    param_a[5] = 9.0466e-7 + 1.4539e-6*log10(0.015204*(y + 3.4)) + 0.00013253/(z*z) + y*(-4.1228e-7 + 2.2036e-7*y);
    param_a[6] = -339.45 + 618.73*log10(0.31595*(y + 3.9)) + 250.20/((y + 4.4395)*(y + 4.4395));
    param_a[7] = -35.105 + y*(36.167 + y*(-9.3575 + 0.33717*y));
    param_a[8] = 0.17554 + y*(0.37300 + y*(-0.014938 + y*(0.0032314 + 0.0025579*y)));
  } //else {
	
    //for (i = 0; i < 9; i++)
  //param_a[i] = 0.0;

  //}

  return param_a;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//		calculate diff. diss. parameters for gamma-rays
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
valarray<double> kamae_gamma_param_diff(double Tp) {

  valarray<double> param_b(0., 8);

  double y, z1, z2, pow;
  int i;
  
  y = log10(Tp*0.001);
  
  if ((Tp > 1.94) && (Tp < 512000.1)) {

    if (Tp > 5.51) {

      z1 = y + 0.59913;
      z2 = y + 9.4773;
      param_b[0] = 60.142*tanh(-0.37555*(y + 2.2)) - 5.9564*z1*z1 + 0.0060162*z2*z2*z2*z2;
      z1 = y + 369.13;
      param_b[1] = 35.322 + 3.8026*tanh(-2.4979*(y + 1.9)) - 0.00021870*z1*z1;
      z1 = y + 252.43;
      param_b[2] = -15.732 - 0.082064*tanh(-1.9621*(y + 2.1)) + 0.00023355*z1*z1;
      pow = (y + 1.0444)/(1.0 + 0.27437*(y + 1.0444));
      param_b[3] = -0.086827 + 0.37646*exp(-0.53053*pow*pow);
    
    } else {
	  
      param_b[0] = 0.0;
      param_b[1] = 0.0;
      param_b[2] = 0.0;
      param_b[3] = 0.0;
    
    }
    
    z1 = y + 2.95;
    pow = (y + 2.45) - 0.19717*(y + 2.45)*(y + 2.45);
    param_b[4] = 2.5982 + 0.39131*z1*z1 - 0.0049693*z1*z1*z1*z1 + 0.94131*exp(-24.347*pow*pow);
    z1 = (y - 0.83562)/(1.0 + 0.33933*(y - 0.83562));
    param_b[5] = 0.11198 + y*(-0.64582 + 0.16114*y) + 2.2853*exp(-0.0032432*z1*z1);
    param_b[6] = 1.7843 + y*(0.91914 + y*(0.050118 + y*(0.038096 + y*(-0.027334 + y*(-0.0035556 + 0.0025742*y)))));
    z1 = y + 1.8441;
    param_b[7] = -0.19870 + y*(-0.071003 + 0.019328*y) - 0.28321*exp(-6.0516*z1*z1);
  } //else {

    //for (i = 0; i < 8; i++)
  //param_b[i] = 0.0;
    
  //}

  return param_b;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//		calculate delta(1232) parameters for gamma-rays
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
valarray<double> kamae_gamma_param_delta(double Tp) {

  valarray<double> param_c(0., 5);

  double y, pow;
  int i;
  
  y = log10(Tp*0.001);
  
  if ((Tp < 0.488) || (Tp > 1.95))
    for (i = 0; i < 5; i++)
      param_c[i] = 0.0;
  else {
      
    pow = ((y + 3.1301)/(1.0 + 0.14921*(y + 3.1301)));
    param_c[0] = 2.4316*exp(-69.484*pow*pow) - (6.3003 + 9.5349/y - 0.38121*y*y);
    param_c[1] = 56.872 + y*(40.627 + 7.7528*y);
    param_c[2] = -5.4918 - 6.7872*tanh(4.7128*(y + 2.1)) + 0.68048*y;
    param_c[3] = -0.36414 + 0.039777*y;
    param_c[4] = -0.72807 + y*(-0.48828 - 0.092876*y);
  
  }

  return param_c;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//		calculate res(1600) parameters for gamma-rays
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
valarray<double> kamae_gamma_param_res(double Tp) {

  valarray<double> param_d(0., 5);

  double y, pow;
  int i;
  
  y = log10(Tp*0.001);
  
  /* 06/06/06: removed unneccessary use of pow() to increase performance
     also added use of pow = <expression> */
  if ((Tp < 0.69) || (Tp > 2.76)) {

    for (i = 0; i < 5; i++)
      param_d[i] = 0.0;
    
  } else {

    pow = ((y + 2.9507)/(1.0 + 1.2912*(y + 2.9507)));
    param_d[0] = 3.2433*exp(-57.133*pow*pow) - (1.0640 + 0.43925*y);
    param_d[1] = 16.901 + y*(5.9539 + y*(-2.1257 - 0.92057*y));
    param_d[2] = -6.6638 - 7.5010*tanh(30.322*(y + 2.1)) + 0.54662*y;
    param_d[3] = -1.50648 + y*(-0.87211 - 0.17097*y);
    param_d[4] = 0.42795 + y*(0.55136 + y*(0.20707 + 0.027552*y));
    
  }

  return param_d;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//		calculate non-diff parameters for electrons
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
valarray<double> kamae_elec_param_nd(double Tp) {

  valarray<double> param_a(0., 9);

  double y, z;
  int i;
  
  y = log10(Tp*0.001);
  
  if ((Tp > 0.487) && (Tp < 512000.1)) {

    z = y + 3.3;
    param_a[0] = z*(-0.018639 + z*(2.4315 + z*(-0.57719 + 0.063435*z)));
    param_a[1] = 7.1827e-6 + y*(-3.5067e-6 + y*(1.3264e-6 + y*(-3.3481e-7 + y*(2.3551e-8 + 3.4297e-9*y))));
    z = y + 7.9031;
    param_a[2] = 563.91 - 362.18*log10(2.7187*(y + 3.4)) - 2.8924e4/(z*z);
    param_a[3] = 0.52684 + y*(0.57717 + y*(0.0045336 - 0.0089066*y));
    z = y + 3.32;
    param_a[4] = z*(0.36108 + z*(1.6963 + z*(-0.074456 + z*(-0.071455 + 0.010473*z))));
    param_a[5] = 9.7387e-5 + 7.8573e-5*log10(0.0036055*(y + 4.3)) + 0.00024660/(y + 4.9390) - 3.8097e-7*y*y;
    param_a[6] = -273.00 - 106.22*log10(0.34100*(y + 3.4)) + 89.037*y - 12.546*y*y;
    z = y + 8.5518;
    param_a[7] = 432.53 - 883.99*log10(0.19737*(y + 3.9)) - 4.1938e4/(z*z);
    param_a[8] = -0.12756 + y*(0.43478 + y*(-0.0027797 - 0.0083074*y));
  
  } //else {

    //for (i = 0; i < 9; i++)
  //param_a[i] = 0.0;
    
  //}

  return param_a;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//		calculate diff. diss. parameters for electrons
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
valarray<double> kamae_elec_param_diff(double Tp) {

  valarray<double> param_b(0., 8);

  double y, z1, z2, pow;
  int i;
  
  y = log10(Tp*0.001);
  
  if ((Tp > 1.94) && (Tp < 512000.1)) {

    if (Tp > 5.51) {

      z1 = y + 1.6878;
      z2 = y + 9.6400;
      param_b[0] = 0.20463*tanh(-6.2370*(y + 2.2)) - 0.16362*z1*z1 + 3.5183e-4*z2*z2*z2*z2;
      pow = (y + 2.0154)/(1.0 + 0.62779*(y + 2.0154));
      param_b[1] = 1.6537 + 3.8530*exp(-3.2027*pow*pow);
      z1 = y + 256.63;
      param_b[2] = -10.722 + 0.082672*tanh(1.8879*(y + 2.1)) + 0.00014895*z1*z1;
      pow = (y + 1.9877)/(1.0 + 0.40300*(y + 1.988));
      param_b[3] = -0.023752 - 0.51734*exp(-3.3087*pow*pow);
    
    } else {

      param_b[0] = 0.0;
      param_b[1] = 0.0;
      param_b[2] = 0.0;
      param_b[3] = 0.0;

    }
    
    z1 = y + 2.9;
    param_b[4] = 0.94921 + 0.12280*z1*z1 - 7.1585e-4*z1*z1*z1*z1 + 0.52130*log10(z1);
    param_b[5] = -4.2295 - 1.0025*tanh(9.0733*(y + 1.9)) - 0.11452*(y - 62.382);
    param_b[6] = 1.4862 + y*(0.99544 + y*(-0.042763 + y*(-0.0040065 + 0.0057987*y)));
    z1 = y - 2.8542;
    param_b[7] = 6.2629 + 6.9517*tanh(-0.36480*(y + 2.1)) - 0.26033*z1*z1;
  
  } //else {

    //for (i = 0; i < 8; i++)
  //param_b[i] = 0.0;
  
  //}

  return param_b;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//		calculate delta(1232) parameters for electrons
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
void kamae_elec_param_delta(double Tp, double* param_c) {

  param_c[0] = param_c[1] = param_c[2] = param_c[3] = param_c[4] = 0.0;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//		calculate res(1600) parameters for electrons
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
valarray<double> kamae_elec_param_res(double Tp) {

  valarray<double> param_d(0., 5);

  double y, pow;
  int i;
  
  y = log10(Tp*0.001);
  
  if ((Tp < 0.69) || (Tp > 2.76)) {

    for (i = 0; i < 5; i++)
      param_d[i] = 0.0;
  
  } else {

    pow = (y + 2.9537)/(1.0 + 1.5221*(y + 2.9537));
    param_d[0] = 0.37790*exp(-56.826*pow*pow) - 0.059458 + 0.0096583*y*y;
    param_d[1] = -5.5135 - 3.3988*y;
    param_d[2] = -7.1209 - 7.1850*tanh(30.801*(y + 2.1)) + 0.35108*y;
    param_d[3] = -6.7841 - 4.8385*y - 0.91523*y*y;
    param_d[4] = -134.03 - 139.63*y - 48.316*y*y - 5.5526*y*y*y;
  
  }

  return param_d;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//		calculate non-diff parameters for positrons
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
valarray<double> kamae_posi_param_nd(double Tp) {

  valarray<double> param_a(0., 9);

  double y, z;
  int i;
  
  y = log10(Tp*0.001);
  
  if ((Tp > 0.487) && (Tp < 512000.1)) {

    z = y + 3.3;
    param_a[0] = z*(-0.79606 + z*(7.7496 + z*(-3.9326 + z*(0.80202 - 0.054994*z))));
    param_a[1] = 6.7943e-6 + y*(-3.5345e-6 + y*(6.0927e-7 + y*(2.0219e-7 + y*(5.1005e-8 - 4.2622e-8*y))));
    param_a[2] = 44.827 + 81.378*log10(0.027733*(y + 3.5)) - 1.3886e4/((y + 8.4417)*(y + 8.4417));
    param_a[3] = 0.52010 + y*(0.59336 + y*(0.012032 - 0.0064242*y));
    z = y + 3.32;
    param_a[4] = z*(2.1361 + z*(1.8514 + z*(-0.47872 + z*(0.0032043 + 0.0082955*z))));
    param_a[5] = 1.0845e-6 + 1.4336e-6*log10(0.0077255*(y + 4.3)) + 0.00013018/((y + 4.8188)*(y + 4.8188)) + 9.3601e-8*y;
    param_a[6] = -267.74 + 14.175*log10(0.35391*(y + 3.4)) + y*(64.669 - 7.7036*y);
    param_a[7] = 138.26 - 529.84*log10(0.12467*(y + 3.9)) - 1.9869e4/((y + 7.6884)*(y + 7.6884)) + 1.0675*y*y;
    param_a[8] = -0.14707 + y*(0.40135 + y*(0.0039899 - 0.0016602*y));
  
  } 
  //else {

  //for (i = 0; i < 9; i++)
  //  param_a[i] = 0.0;
  
  //}

  return param_a;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//		calculate diff. diss. parameters for positrons
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
valarray<double> kamae_posi_param_diff(double Tp) {

  valarray<double> param_b(0., 8);

  double y, z1, z2, pow;
  int i;
  
  y = log10(Tp*0.001);
  
  if ((Tp > 1.94) && (Tp < 512000.1)) {

    if (Tp > 11.0) {

      z1 = y + 0.67500;
      z2 = y + 9.0824;
      param_b[0] = 29.192*tanh(-0.37879*(y + 2.2)) - 3.2196*z1*z1 + 3.6687e-3*z2*z2*z2*z2;
      pow = (y + 1.8781)/(1.0 + 3.8389*(y + 1.8781));
      param_b[1] = -142.97 + 147.86*exp(-0.37194*pow*pow);
      z1 = y + 234.65;
      param_b[2] = -14.487 - 4.2223*tanh(-13.546*(y + 2.2)) + 0.00016988*z1*z1;
      pow = (y + 1.8194)/(1.0 + 0.99946*(y + 1.8194));
      param_b[3] = -0.0036974 - 0.41976*exp(-6.1527*pow*pow);
    
    } else {

      param_b[0] = 0.0;
      param_b[1] = 0.0;
      param_b[2] = 0.0;
      param_b[3] = 0.0;
    
    }
    
    z1 = y + 2.95;
    pow = y + 2.29 - 0.18967*(y + 2.29);
    param_b[4] = 1.8108 + z1*z1*(0.18545 - 2.0049e-3*z1*z1) + 0.85084*exp(-14.987*pow*pow);
    param_b[5] = 2.0404 - 0.51548*tanh(2.2758*(y + 1.9)) - 0.035009*(y - 6.6555);
    param_b[6] = 1.5258 + y*(1.0132 + y*(-0.064388 + y*(-0.0040209 + 0.0082772*y)));
    z1 = y - 2.7718;
    param_b[7] = 3.0551 + 3.5240*tanh(-0.36739*(y + 2.1)) - 0.13382*z1*z1;
  
  } //else {

    //for (i = 0; i < 8; i++)
  //param_b[i] = 0.0;
  
  //}

  return param_b;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//		calculate delta(1232) parameters for positrons
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
valarray<double> kamae_posi_param_delta(double Tp) {

  valarray<double> param_c(0., 5);

  double y, pow;
  int i;
  
  y = log10(Tp*0.001);
  
  if ((Tp < 0.488) || (Tp > 1.95))
    for (i = 0; i < 5; i++)
      param_c[i] = 0.0;
  else {
    
    pow = ((y + 3.1272)/(1.0 + 0.22831*(y + 3.1272)));
    param_c[0] = 2.9841*exp(-67.857*pow*pow) - (6.5855 + 9.6984/y - 0.41256*y*y);
    param_c[1] = 6.8276 + 5.2236*y + 1.4630*y*y;
    param_c[2] = -6.0291 - 6.4581*tanh(5.0830*(y + 2.1)) + 0.46352*y;
    param_c[3] = 0.59300 + 0.36093*y;
    param_c[4] = 0.77368 + 0.44776*y + 0.056409*y*y;
  
  }

  return param_c;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//		calculate res(1600) parameters for positrons
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
valarray<double> kamae_posi_param_res(double Tp) {

  valarray<double> param_d(0., 5);

  double y, pow;
  int i;
  
  y = log10(Tp*0.001);
  
  if ((Tp < 0.69) || (Tp >= 2.76)) {

    for (i = 0; i < 5; i++)
      param_d[i] = 0.0;
  
  } else {

    pow = ((y + 2.9485)/(1.0 + 1.2892*(y + 2.9485)));
    param_d[0] = 1.9186*exp(-56.544*pow*pow) - (0.23720 - 0.041315*y*y);
    param_d[1] = -4.9866 - 3.1435*y;
    param_d[2] = -7.0550 - 7.2165*tanh(31.033*(y + 2.1)) + 0.38541*y;
    param_d[3] = -2.8915 - 2.1495*y - 0.45006*y*y;
    param_d[4] = -1.2970 - 0.13947*y + 0.41197*y*y + 0.10641*y*y*y;
    
  }

  return param_d;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// non-diff differential cross section, mbarn/GeV
//
// E        energy of secondary particle (gamma, e-, e+) in GeV
// Tp       proton kinetic energy in GeV
// particle the particle type, i.e. gamma, e- or e+
//
// return differential cross section dsigma/dlogE in mbarn/GeV 
//**.****!****.****!****.****!****.****!****.****!****.****!****.****!****.****!
double kamae_nd(double E, double Tp, valarray<double> param_a, int particle) {

  double Wl,Wh, Lmin, Lmax;
  double x, xa3, xa8;
  double y;
  double pow1, pow2;
  double sigma;
  double factor, r_factor;
  
  // Table 2 of Kamae et al.
  double L_MAX[3] = {0.96, 0.96, 0.94};
  double W_NDL[3] = {15.0, 20.0, 15.0};
  double W_NDH[3] = {44.0, 45.0, 47.0};
  
  // init some variables, given in table 2
  Lmin = -2.6;
  Lmax = L_MAX[particle]*log10(Tp);
  Wl = W_NDL[particle];
  Wh = W_NDH[particle];
  
  // calculate log(E) and log(Tp)
  x = log10(E);
  y = log10(Tp*0.001);
  
  // calculate the flux due to non-diffractive process for given gamma-ray energy
  xa3 = x - param_a[3];
  pow1 = xa3*(1 + param_a[2]*xa3);
  xa8 = x - param_a[8];
  pow2 = xa8*(1 + xa8*(param_a[6] + param_a[7]*xa8));
  sigma = param_a[0]*exp(-param_a[1]*pow1*pow1) + param_a[4]*exp(-param_a[5]*pow2*pow2);
  
  // factor is the kinematic limit function as in the paper
  factor = (1.0/(1.0 + exp(Wl*(Lmin - x))))*(1.0/(1.0 + exp(Wh*(x - Lmax))));
  sigma = sigma*factor;
  
  if (sigma < 0.0)
    sigma = 0.0;
  
  // renormalization
  r_factor = 1.0;
  
  switch (particle) {

  case 0: // gamma-rays
    if (Tp <= 1.95) {

      pow1 = (y + 3.25)/(1.0 + 8.08*(y + 3.25));
      r_factor = 3.05*exp(-107.0*pow1*pow1);
    }
    else
      r_factor = 1.01;
    break;
  case (-1): // electrons
    if (Tp <= 15.6) {

      pow1 = (y + 3.26)/(1.0 + 9.21*(y + 3.26));
      r_factor = 3.63*exp(-106*pow1*pow1) + y*(-0.182 - 0.175*y);
    }
    else
      r_factor = 1.01;
    break;
    
  case 1: // positrons
    if (Tp <= 5.52) {

      pow1 = (y + 3.25)/(1.0 + 10.4*(y + 3.25));
      r_factor = 2.22*exp(-98.9*pow1*pow1);
    }
    break;
    
  }	
  
  sigma = sigma*r_factor;
  
  return sigma;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// diff. diss. differential cross section, mbarn/GeV
//
// E       energy of secondary particle (gamma, e-, e+) in GeV
// Tp      proton kinetic energy in GeV
//
// return differential cross section dsigma/dlogE in mbarn/GeV 
//**.****!****.****!****.****!****.****!****.****!****.****!****.****!****.****!
double kamae_diff(double E, double Tp, valarray<double> param_b) {

  double W=75.0, Lmax=log10(Tp);
  double x=log10(E);
  double pow1, pow2;
  double sigma;
  double factor;
  
  // calculate the sigma due to diffractive process for given gamma-ray energy
  pow1 = (x - param_b[2])/(1.0 + param_b[3]*(x - param_b[2]));
  pow2 = (x - param_b[6])/(1.0 + param_b[7]*(x - param_b[6]));
  sigma = param_b[0]*exp(-param_b[1]*pow1*pow1) + param_b[4]*exp(-param_b[5]*pow2*pow2);
  
  // factor is the kinematic limit function as in the paper
  factor = 1.0/(1.0 + exp(W*(x - Lmax)));
  sigma = sigma*factor;
  
  if (sigma < 0.0)
    sigma = 0.0;
  
  return sigma;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// delta(1232) differential cross section, mbarn/GeV
//
// E       energy of secondary particle (gamma, e-, e+) in GeV
// Tp      proton kinetic energy in GeV
//
// return differential cross section dsigma/dlogE in mbarn/GeV 
//**.****!****.****!****.****!****.****!****.****!****.****!****.****!****.****!
double kamae_delta(double E, double Tp, valarray<double> param_c) {

  double W=75.0, Lmax=log10(Tp);
  double x=log10(E), xc2;
  double pow;
  double sigma;
  double factor;
  
  // calculate the sigma due to resonance process for given gamma-ray energy
  xc2 = x - param_c[2];
  pow = xc2/(1.0 + xc2*(param_c[3] + param_c[4]*xc2));
  sigma = param_c[0]*exp(-param_c[1]*pow*pow);
  
  // factor is the kinematic limit function as in the paper
  factor = 1.0/(1.0 + exp(W*(x - Lmax)));
  sigma = sigma*factor;
  
  if (sigma < 0.0)
    sigma = 0.0;
  
  return sigma;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// res(1600) differential cross section, mbarn/GeV
//
// E       energy of secondary particle (gamma, e-, e+) in GeV
// Tp      proton kinetic energy in GeV
//
// return differential cross section dsigma/dlogE in mbarn/GeV 
//**.****!****.****!****.****!****.****!****.****!****.****!****.****!****.****!
double kamae_res(double E, double Tp, valarray<double> param_d) {

  return kamae_delta(E, Tp, param_d);

}


//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// differential cross section for pp->pi0->gamma,barn/GeV
//
// Esec    energy of secondary particle (gamma, e-, e+) in GeV
// Pp2     primary proton momentum in GeV/c
// NA11
// NA21
// key     specify type of secondary particle, 0: gamma, 1: e+ and -1: e-
// params  array of four pointers to arrays of doubles with the current
//         parameter values
//
// return  differential cross section dsigma/dE in barn/GeV
//**.****!****.****!****.****!****.****!****.****!****.****!****.****!****.****!
double 
kamae(double Esec, double Pp2, int NA11, int NA21, int key, 
      valarray<double> params0, valarray<double> params1, 
      valarray<double> params2, valarray<double> params3) {

  double s_nd,s_diff,s_delta,s_res;
  double Pp1=Pp2;
  int NA1=NA11;
  int NA2=NA21;
  double Pp=Pp1/NA1;                   // momentum per nucleon
  double Tp=sqrt(Pp*Pp+Mp*Mp)-Mp;
  
  // calculate dsigma/dlogE for each contributioon
  s_nd = 1.0e-3*kamae_nd(Esec,Tp,params0,key);
  s_diff = 1.0e-3*kamae_diff(Esec,Tp,params1);
  s_delta = 1.0e-3*kamae_delta(Esec,Tp,params2);
  s_res = 1.0e-3*kamae_res(Esec,Tp,params3);
  
  // 1.0/Esec because dsigma/dlogE = 1.0/E * dsigma/dE
  return (s_nd+s_diff+s_delta+s_res)/Esec;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// Test routine
//
// Calculates gamma-ray spectrum for power-law protons of index=2 as in
// figure 11 of Kamae et al., ApJ, 647, 692.
//
// Compile: g++ -o kamae kamae.cc
// Run: ./kamae > spectrum.dat
//**.****!****.****!****.****!****.****!****.****!****.****!****.****!****.****!

/*
int main()
{
  // the proton kinetic energies to calculate for
  double Tp[41] = {512.0e3, 362.0e3, 256.0e3, 181.0e3, 128.0e3, 90.5e3, 64.0e3, 45.3e3, 32.0e3, 22.6e3, 
		   16.0e3, 11.3e3, 8.0e3, 5.66e3, 4.0e3, 2.8e3, 2.0e3, 1.41e3, 1.0e3, 707.0, 500.0, 354.0,
		   250.0, 177.0, 125.0, 88.4, 62.5, 44.2, 31.3, 22.1, 15.6, 11.1, 7.81, 5.52, 3.91, 2.76,
		   1.95, 1.38, 0.98, 0.69, 0.488};
  
  // allocate memory for the spectrum (180 bins)
  double* spectrum = new double[180];
  
  // key = 0 for gamma-rays
  int key = 0;
  
  // loop over all Tp
  for (int i = 0; i < 41; i++) {
    double param_a[9];
    double param_b[8];
    double param_c[5];
    double param_d[5];
    
    kamae_gamma_param_nd(Tp[i], param_a);
    kamae_gamma_param_diff(Tp[i], param_b);
    kamae_gamma_param_delta(Tp[i], param_c);
    kamae_gamma_param_res(Tp[i], param_d);
    
    double* params[4] = {param_a, param_b, param_c, param_d};
    // index factor for normalization to proton power-law of index=2
    double ind2factor = 1.0/(Tp[i]*1.0e-3);
    
    // loop over all gamma-ray energy bins and calculate
    for (int j = 0; j < 180; j++) {
      // calculate the gamma-ray energy, ranging from 10^-3 to 10^6 GeV in
      // log bins 0.05 wide
      double E = pow(10.0, j*0.05 - 3.0);
      double Pp = sqrt(Tp[i]*Tp[i] + 2*Mp*Tp[i]);
      
      int NA11 = 1;
      int NA21 = 1;
      
      // calculate cross section and add to overall spectrum
      // multiply by 1.0e3 and E to get dsigma/dlogE in mbarn
      double s = kamae(E, Pp, NA11, NA21, key, params)*1.0e3*E;
      // multiply with correct index factor
      spectrum[j] += s*ind2factor;
    }
  }
  // this outputs the spectrum vs. logE in a format readable by hippodraw
  cout << "# gamma-ray spectrum of index=2 protons" << endl;
  cout << "logE" << "\t" << "s" << endl;
  for (int i = 0; i < 180; i++) {
    double logE = i*0.05 - 3.0;
    double s = log10(spectrum[i] + 1.0e-12) + logE;
    cout << logE << " " << s << endl;
  }
  
  delete [] spectrum;
}
*/
