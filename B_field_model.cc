
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * B_field_model.cc *                            galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include<cmath>     //AWS20050624

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// Magnetic field in Tesla (=1e4 Gauss); r,z in kpc

double B_field_model(double r,double z,int model)
{
   float Bo, rscale, zscale;
   float b_field=0.0;
   float ro=8.5;      // Sun galactocentric distance

   if (model==1)
   {
      Bo=6.0e-10;
      rscale=20.;
      zscale=5.;
      b_field=Bo *exp(-(r-ro)/rscale) * exp(-fabs(z)/zscale);
   }

// exactly the model used in cal3prop:
   if (model==2)
   {
      Bo=6.0e-10;
      rscale=20.;
      zscale=5.;
      b_field=Bo *exp(-r/rscale) * exp(-fabs(z)/zscale);
   }

   if (model> 2)
   {
      Bo=model*1.0e-10;
      rscale=20.;
      zscale=5.;
      b_field=Bo *exp(-(r-ro)/rscale) * exp(-fabs(z)/zscale);
   }

// Bo rscale zscale encoded in 9-digit number: BBBrrrzzz in units of 0.1
// e.g. 123456789 : Bo= 12.3E-10 Tesla  rscale=45.6 kpc zscale=78.9 kpc 

   if (model > 1000)
   {
      Bo=           (model/1000000)                 * 0.1 *1.0e-10;        
      rscale=(model-(model/1000000)*1000000 )/1000  * 0.1         ;
      zscale=(model%1000)                           * 0.1         ;
      b_field=Bo *exp(-(r-ro)/rscale) * exp(-fabs(z)/zscale);
   }
   return b_field;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// 3D  interface

double B_field_model(double x, double y, double z, int model)
{
   return B_field_model(sqrt(x*x+y*y),z,model);
}
