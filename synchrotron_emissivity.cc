
#include<iostream>
#include <cmath>
#include <string>
#include <cstdlib>
using namespace std;
//C linkage since GSL is C library, assumed to be present as library compiled with C. 

/*
// these are the natural forms, but give error for large x so use error-handling versions instead
extern "C" double     gsl_sf_synchrotron_1                     (const double x);
extern "C" double     gsl_sf_synchrotron_2                     (const double x);
extern "C" double     gsl_sf_bessel_Knu_scaled(const double nu, const double x);
extern "C" double     gsl_sf_bessel_Knu       (const double nu, const double x);
*/

// error-handling versions are preferred

typedef struct{double val;double err;} gsl_sf_result; // structure used by GSL routines

extern "C" int        gsl_sf_synchrotron_1_e                     (const double x, gsl_sf_result *result);
extern "C" int        gsl_sf_synchrotron_2_e                     (const double x, gsl_sf_result *result);
extern "C" double     gsl_sf_bessel_Knu_scaled_e(const double nu, const double x, gsl_sf_result *result);
extern "C" double     gsl_sf_bessel_Knu_e       (const double nu, const double x, gsl_sf_result *result);

extern "C" void       gsl_set_error_handler_off();


//////////////////////////////////////////////////////////////////////////////////
double synchrotron_emissivity(double gamma, double nu,  double Bperp, double Brand,
                              double &synch_emissivity_regular,       double &synch_emissivity_parallel, double &synch_emissivity_perpendicular,
                              double &synch_emissivity_random,
                              int debug=0 )
{
  //  Computes synchrotron emissivity for isotropic distribution of electrons for regular and random magnetic fields. Includes polarization.
  //  Developed by Andy Strong  & Elena Orlando, MPE
  //  First version AWS 20071211
  //  It has been tested against the galprop fortran version and by integrating over angle to get the random from the regular 
  //  (see appended test program). The agreement is excellent.


  //  input  arguments
  //   double gamma  electron Lorentz factor;
  //   double nu     synchrotron frequency, Hz
  //   double Bperp  perpendicular (to observer line-of-sight) component of regular field, Gauss
  //   double Brand  total                                                  random  field, Gauss

  //  output arguments: synchrotron power radiated, integrated over the emission cone
  //   units: erg s-1 Hz^-1
  //   double synch_emissivity_parallel;          // parallel      polarization for regular field
  //   double synch_emissivity_perpendicular;     // perpendicular polarization for regular field
  //   double synch_emissivity_regular;           // total power                for regular field (sum of polarizations)
  //   double synch_emissivity_random;            // total power                for random  field

  //  returned by function:
  //   double synchrotron_emissivity              // total power regular and random         fields
 

 /*
 GSL (GNU Subroutine Library)
 http://www.gnu.org/software/gsl/

 GSL (GNU Subroutine Library) Manual:

 7.29 Synchrotron Functions

The functions described in this section are declared in the header file gsl_sf_synchrotron.h.
 Function: double gsl_sf_synchrotron_1 (double x)
 Function: int gsl_sf_synchrotron_1_e (double x, gsl_sf_result * result)

    These routines compute the first synchrotron function x \int_x^\infty dt K_{5/3}(t) for x >= 0. 

 Function: double gsl_sf_synchrotron_2 (double x)
 Function: int gsl_sf_synchrotron_2_e (double x, gsl_sf_result * result)

    These routines compute the second synchrotron function x K_{2/3}(x) for x >= 0.
*/


  // this is a stand-alone routine so define the constants 

   double pi  = acos(-1.0);
   double c   = 2.99792e10;     // speed of light, cm s^-1
   double m_e =  9.10939e-28 ;  // electron mass, g
   double e_e =  4.80321e-10 ;  // electron charge, esu
   
   // derived constants

   double emc     =      e_e   / ( m_e     * c   );     //  e  /mc
   double e3mc2   =  pow(e_e,3) /( m_e * pow(c,2));     //  e^3/mc^2


 
  synch_emissivity_regular         = 0.;
  synch_emissivity_parallel        = 0.;
  synch_emissivity_perpendicular   = 0.;
  synch_emissivity_random          = 0.;


  ////////////////////////////////////////////////////////////////////////////
  // regular field
  ////////////////////////////////////////////////////////////////////////////

  // Notes on various formulations:
  // sqrt(3) e^3/mc^2 Bperp  F(nu/nu_c) : Lang 2006 Vol 1 p. 30 equations  1.158,159,160) 
  // 1.158,1.159 have 1/2pi, 1.160 has no 1/pi, because 1.158,159 are per sr, and 1/2 from integration of sin theta.
  // Longair High Energy Astrophysics Second Edition Vol 2 Chapter 18.1 equation 18.31, 18.35 (NB uses angular freq omega, not freq)
  // omega=2pi nu  , d/dnu=2pi d/domega (see Longair eq 18.36)
  // Longair 18.36: total(omega)=sqrt(3)/(8pi^2 e0 m_e c) F(x) - SI units, e0=1/mu0 c^2 = 1/(4pi 10^-7 c^2) -> sqrt(3)/(2pi/c). 
  //-> total(nu)=sqrt(3) like Lang 1.160 and Ginberg&Syrovatskii 1964 eq. 4.20
  // NB the GSL functions  are exactly the  F(x), G(x) of Longair

  int status;
  gsl_sf_result result;

  gsl_set_error_handler_off(); // since we check the status and act on it

  if(Bperp>0.0)
  {

  double nu_c = 3./(4.*pi) * emc  * Bperp * gamma * gamma ; // critical frequency
  double    x = nu/nu_c;
  
  
  double F=0.0; status= gsl_sf_synchrotron_1_e(x,&result) ;if(status==0) F=result.val;
  double G=0.0; status= gsl_sf_synchrotron_2_e(x,&result) ;if(status==0) G=result.val;

  double constant_reg=sqrt(3.)/2. *  e3mc2 * Bperp ;
      
  synch_emissivity_parallel       = constant_reg *  (F - G); // parallel      polarization 
  synch_emissivity_perpendicular  = constant_reg *  (F + G); // perpendicular polarization

  synch_emissivity_regular = synch_emissivity_parallel + synch_emissivity_perpendicular; // total regular
                                                                 // = sqrt(3) e^3/mc^2 F(x) as in Ginzberg&Syrovatskii 1964 eq 4.20, Lang 1.160

  if(debug==1)
  {
   cout<<"gamma= "<<gamma<<" nu="<<nu<<" Bperp="<<Bperp<<" nu_c= "<<nu_c << " Hz"<<endl;
   cout <<" x="<<x<<" F(x)="<<F<<" G(x)="<<G
      <<" total reg synch emissivity=  "<< synch_emissivity_regular
      <<" polarized: parallel =  "<< synch_emissivity_parallel
      <<               " perp =  "<< synch_emissivity_perpendicular
      <<endl;
  }// debug

  }// if

  ////////////////////////////////////////////////////////////////////////////
  // random field
  ////////////////////////////////////////////////////////////////////////////

  // Strong et al. 2000 ApJ 537,763 Appendix B, quoting Ghisellini, Guilbert & Svenson 1988, ApJ 334, L5
  // total power erg s-1 Hz-1 = 4sqrt(3) pi r_e m_e c nuB x^2(K4/3 K 1/3 - 3/5 x (K4/3^2-K1/3^2))
  // where nuB=eB/(2pi m_e c), x=nu/(3 gamma^2 nuB) and  B is the total random field
  //  4sqrt(3) pi r_e m_e c nuB=  4sqrt(3) pi e^2/(m_e c^2) m_e c eB/(2pi m_e c) 
  // =  2sqrt(3) e^3/(m_e c^2) B      ( twice the formula in Lang 1.160 for reg field using Bperp)
  // x=nu/(3*gamma^2*nuB)-> nu_c=3*eB/(2pi m_e c)*gamma^2 = 3/2pi e/mc B gamma^2 (cf 3/4pi e/mc Bperp gamma^2 for regular B)




  if(Brand>0.) 
  {  
  
  double nu_c = 3./(2.*pi) * emc  * Brand * gamma * gamma ; // critical frequency
  double x=nu/nu_c;

  double K13=0.0; status = gsl_sf_bessel_Knu_e (1./3.,  x, &result); if(status==0) K13=result.val; // modified Bessel function of order 4/3
  double K43=0.0; status = gsl_sf_bessel_Knu_e (4./3.,  x, &result); if(status==0) K43=result.val; // modified Bessel function of order 4/3

  double  constant_random=2.*sqrt(3.) *  e3mc2 * Brand ;
  
  synch_emissivity_random = constant_random *  x*x* (K43*K13 -3./5. *x*(K43*K43 - K13*K13));

 


  if(debug==1)
  {
  cout<<"gamma= "<<gamma<<" nu="<<nu<<" Brand="<<Brand<<" nu_c= "<<nu_c << " Hz"<<endl;
  cout <<" x="<<x<<" K1/3="<<K13<<" K4/3="<<K43
      <<"  synch emissivity for random field=  "<< synch_emissivity_random
      <<endl;
  }// debug

 } //if

  double synch_emissivity_total =    synch_emissivity_regular + synch_emissivity_random;

  return  synch_emissivity_total;

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double synchrotron_cc_new(double gamma, double nu, double Brand)
{
  // this is a substitute for the original galprop fortran routine, to enable a direct replacement without changing galprop code
  // the fortran code is for random field only, using the same Ghisellini et al. formula as here

 double Bperp=0.e-6; //dummy value: this routine is just for the random field

 double synch_emissivity_total;
 double synch_emissivity_reg;
 double synch_emissivity_par;
 double synch_emissivity_perp;
 double synch_emissivity_random;

 int debug=0;

 synch_emissivity_total
     = synchrotron_emissivity( gamma, nu, Bperp, Brand,
                               synch_emissivity_reg, synch_emissivity_par, synch_emissivity_perp, synch_emissivity_random, debug );

  return    synch_emissivity_random;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// tests
/*

typical usage;

f95i synchrotron_galprop.f -c

GSL library is assumed to exist:

icpc synchrotron_emissivity.cc /afs/ipp-garching.mpg.de/home/a/aws/gsl/gsl-1.10/olga/lib/libgsl.a synchrotron_galprop.o -lifcore -limf -lm -L/afs/ipp/i386_linux24/soft/intel/compiler64.f90/lib

setenv LD_LIBRARY_PATH /afs/ipp/i386_linux24/soft/intel/compiler64.f90/lib:/afs/ipp/i386_linux24/soft/intel/compiler64.f90/lib:/usr/lib/qt3/lib
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// prefix  '//' on next line to test the routine
 /* 

extern "C" double synchrotron_(double*,double*,double*); // galprop fortran version for tests


int main()
{


  cout<<"testing synchrotron_emissivity.cc"<<endl;

  double x,F,G;
  double gamma,nu,Bperp,Brand;


 for(x=1e-6;x<10.;x*=1.1)
 {
  F= gsl_sf_synchrotron_1(x);
  G= gsl_sf_synchrotron_2(x);

  cout<<"x="<<x<<" F(x)="<<F<<" G(x)="<<G<<endl;
 }

 /////////////////////////////////////////////////////////////////////////
 double nuK=4./3.;
 x=1.51;
 cout<<"irreg modified Bessel function of second kind Knu nu="<<nuK<<" x="<<x<<" Knu="<<gsl_sf_bessel_Knu (nuK,  x)<<endl;
 nuK=1./3.;
 cout<<"irreg modified Bessel function of second kind Knu nu="<<nuK<<" x="<<x<<" Knu="<<gsl_sf_bessel_Knu (nuK,  x)<<endl;
 nuK=5./3.;
 cout<<"irreg modified Bessel function of second kind Knu nu="<<nuK<<" x="<<x<<" Knu="<<gsl_sf_bessel_Knu (nuK,  x)<<endl;
 nuK=2./3.;
 cout<<"irreg modified Bessel function of second kind Knu nu="<<nuK<<" x="<<x<<" Knu="<<gsl_sf_bessel_Knu (nuK,  x)<<endl;
 
 G= gsl_sf_synchrotron_2(x);
 nuK=2./3.;
 cout<<"irreg modified Bessel function of second kind Knu nuK="<<nuK<<" x="<<x<<" x*Knu="<<x*gsl_sf_bessel_Knu (nuK,  x)<<endl;
 cout<<" comparing with GSL synchrotron function x="<<x<<" G(x)="<<G<<endl;



 //////////////////////////////////////////////////////////////////////


 gamma=1e4/.511; // 10 GeV electrons = 1e4 MeV
 nu=1.e9; // 1000 MHz
 Bperp=6e-6;// 6 microgauss
 Brand=6e-6;

 double synch_emissivity_total;
 double synch_emissivity_reg;
 double synch_emissivity_par;
 double synch_emissivity_perp;
 double synch_emissivity_random;

 int debug=1;

 synch_emissivity_total
     = synchrotron_emissivity( gamma, nu, Bperp, Brand,
                         synch_emissivity_reg, synch_emissivity_par, synch_emissivity_perp, synch_emissivity_random, debug );

 cout<<"testing synchrotron routine"<<endl;
 cout<<endl;
 cout<<"gamma= "<<gamma<<" nu="<<nu<<" Bperp="<<Bperp   << " Brand="<<Brand   <<endl

  <<" synch emissivity random =  " << synch_emissivity_random <<endl
  <<" synch emissivity reg    = "  << synch_emissivity_reg    <<endl
  <<" synch emissivity parallel = "<< synch_emissivity_par    <<endl
  <<" synch emissivity perp = "    << synch_emissivity_perp   <<endl 
  <<" synch emissivity total =  "  << synch_emissivity_total  <<endl       
  <<endl;




 // check the Ghisellini et al. formula by integrating the regular field emission over all angles weighted by sin(theta)

 double theta;
 double dtr=acos(-1.)/180.; // degrees to radians
 double sum_reg     =0;
 double sum_sintheta=0;

 debug=0;

 for (theta=1;theta<179;theta++)
 {
   Bperp=Brand*sin(theta*dtr);
   synch_emissivity_total= synchrotron_emissivity( gamma, nu, Bperp, Brand,  
                                              synch_emissivity_reg, synch_emissivity_par, synch_emissivity_perp, synch_emissivity_random, debug );

   cout<<"gamma= "<<gamma<< " nu="<<nu<< " theta="<<theta<< " Bperp="<<Bperp  << " Brand="<<Brand 
       <<" synch emissivity random = "  << synch_emissivity_random 
       <<                 " reg = "     << synch_emissivity_reg    
       <<                 " parallel = "<< synch_emissivity_par     
       <<                 " perp = "    << synch_emissivity_perp 
       <<                 " total = "   << synch_emissivity_total           
       <<endl;

   sum_reg     += synch_emissivity_reg *  sin(theta*dtr); // sin(theta) weighted average
   sum_sintheta+=                         sin(theta*dtr);

   }

 double synch_emissivity_reg_B_averaged = sum_reg/sum_sintheta;

 cout<<" random field emissivity based on Ghisellini formula                 : "<<  synch_emissivity_random        <<endl;
 cout<<" sin(theta)-weighted average regular field emissivity => random field: "<<  synch_emissivity_reg_B_averaged<<endl;


 // compare with fortran version from galprop
 double synch_emissivity_fort=synchrotron_( &gamma, &nu, &Brand);

 cout<<"galprop fortran version: "
     <<" synch emissivity=  "<< synch_emissivity_fort<<" diff new version cf fortran="<<synch_emissivity_random - synch_emissivity_fort
     <<" frac diff="<< (synch_emissivity_random - synch_emissivity_fort )/synch_emissivity_random  <<endl;


 cout<<" new              version: "
     <<" synch emissivity random=  "<< synch_emissivity_random<<" diff cf sintheta av="<<synch_emissivity_random - synch_emissivity_reg_B_averaged
     <<" frac diff="<< (synch_emissivity_random - synch_emissivity_reg_B_averaged )/synch_emissivity_random  <<endl;



 cout<<"testing routine synchrotron_cc_new (fortran replacement): emissivity= "<<synchrotron_cc_new(gamma,nu,Brand)<<endl;

  return 0;
}

// prefix '//' on next line to test the routine
 */
