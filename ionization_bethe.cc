
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * ionization_bethe.cc *                         galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include<cmath>
#include<iostream>

/*
 Ionization cross-section in cm^2 from Bethe (1933)
 as quoted by Spitzer and Tomasko ApJ 152, 971, (1968)

The factor 5/3 is included.
*/
double ionization_bethe(int Z, double beta){

double ionization_bethe_;
double beta2= beta*beta;

 ionization_bethe_= 1.23e-20 * Z * Z /beta2
                    * (6.20 + log10(beta2/(1.0-beta2)) -0.43*beta2)
                    * 5.0/3.0;




return ionization_bethe_;
}
