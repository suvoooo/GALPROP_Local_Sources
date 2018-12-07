//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * galprop.h *                                   galprop package * 08/16/2001 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
#ifndef _galprop_internal_h_
#define _galprop_internal_h_

#include <constants.h>

#include <cmath>
#include <iostream>
#include <string>
#include <valarray>

using namespace std;                                             //AWS20050912

//#include "fort_interface.h"

// function prototypes

// only those functions which are not in Galprop class
// and hence stand-alone.

double eprop(double,double,double,double,double,double,double,double,double,double);//IMOS20061030
void Kcapture_cs(double,int,int,double*,double*);               // IMOS20010816
void nucleon_cs(int,double,int,int,int,double*,double*,double*,double*,double*,double*);// IMOS20010511
double isotope_cs(double,int,int,int,int,int,int*);
void read_nucdata(const string& path);
void cleanup_nucdata();                                         //IMOS20060420
double nucdata(int,int,int,int,int,int,int*,int*,double*);      // IMOS20010816
double nucleon_loss(int,int,double,double,double,double,        // nucleon
		    double*, double*);                                    // energy losses
double electron_loss(double,double,double,double,double,double, // electron
		     double*,double*,double*,double*,double*,double*);     // energy losses
double blattnig_gamma(double,double,int,int,int);               // gammas from pi0-decay, Blattnig etal. formalism

double kamae(double,double,int,int,int, valarray<double>, valarray<double>, valarray<double>, valarray<double>);//double**);              // gammas from pi0-decay, Kamae etal. formalism IMOS20061114
 // functions for component parameter calculation IMOS20061114
valarray<double> kamae_gamma_param_nd(double);//, double*);
valarray<double> kamae_gamma_param_diff(double);//, double*);
valarray<double> kamae_gamma_param_delta(double);//, double*);
valarray<double> kamae_gamma_param_res(double);//, double*);

double sigma_boron_dec_heinbach_simon(int, int, int, int, double);

//double sim(double,double,double,double,double,double(*)(double)); //integration
int kinematic(int, int, char*, double&, double&, double&, double&, double&, double&, int);

int tridag(float*, float*, float*, float*, float*, int);
int tridag(double*, double*, double*, double*, double*, int); //IMOS20030217
int tridag_sym(float*, float*, float*, float*, float*, int);
int tridag_sym(double*, double*, double*, double*, double*, int); //IMOS20030217
double sim(double, double, double, double, double, double(*)(double)); //integration

int He_to_H_CS(double, int, int, int, int, double*, double*);

double aic_cc(int, int, double, double, double, double, double, double, double, double, double);//IMOS20060420
double fjones_cc(double, double, double);          //IMOS20060420

double ionization_bethe(int Z, double beta);

double nHI(double, double);
double nH2(double, double);
double nH2(double, double, double);//AWS20090616
double nHII(double, double);
int    nH_set_model(int nHI_model_, int nH2_model_, int nHII_model_, int debug_); //AWS20090814

//double nHI_av (double,double,double,double,double); //IMOS20080114
//double nH2_av (double,double,double,double,double); //IMOS20080114
//double nHII_av(double,double,double,double,double); //IMOS20080114

double nH_av(double, double, double, double, double, double(*)(double, double)); //IMOS20080114
double nH_av(double, double, double, double, double, double, double(*)(double, double, double));                      //AWS20090616
double nH_av_lb(double, double, double, double, double, double, double, double(*)(double, double)); //IMOS20080114
double nH_av_lb(double, double, double, double, double, double, double, double,double(*)(double, double, double)); //AWS20090616


double B_field_model(double, double, int);
double B_field_model(double, double, double, int);

double gauss(double mean, double sigma);

int test_sigma_boron_dec_heinbach_simon();
int test_kinematic();
int test_He_to_H_CS();
int test_nH(); 
int test_Distribution();
int test_float_accuracy();


#endif
