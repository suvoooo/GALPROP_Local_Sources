using namespace std;
#include<iostream>


#include "B_field_3D_model.h"
#include "synchrotron_emissivity.h"
#include "synchrotron_emissivity_B_field.h"

///////////////////////////////////////////////////////////////


double synchrotron_emissivity_B_field( double gamma, double nu, 
				       char *name, double *parameters,
				       double x, double y, double z,  
                                       int options,
                                       double x0, double y0, double z0,
                                       double     &synch_emissivity_reg, double &synch_emissivity_par, double &synch_emissivity_perp, double &synch_emissivity_random,
                                       int debug=0 )

{

  
 
  double  Breg,          Bregx,         Bregy,        Bregz;
  double  Bran,          Branx,         Brany,        Branz;

  

 B_field_3D_model ( name, parameters, x,  y,  z,  options,
                 Breg,  Bregx, Bregy, Bregz,
		 Bran,  Branx, Brany, Branz,
                 debug );

 if(debug==1)
 {

 cout<<"(x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  ";
 cout<<" B model= "<<name<<"  Breg="<<  Breg <<" Bregx="<<  Bregx<<" Bregy="<< Bregy<<" Bregz="<< Bregz
                          <<" Bran="<<  Bran <<" Branx="<<  Branx<<" Brany="<< Brany<<" Branz="<< Branz<<endl;

 cout<<"Bperp for regular field ="<<Bperp( x, y, z,  Bregx, Bregy, Bregz,  x0,y0,z0)<<endl;
 cout<<"Bperp for random  field ="<<Bperp( x, y, z,  Branx, Brany, Branz,  x0,y0,z0)<<endl;


 }

 // use synchrotron routine

 double synch_emissivity_total;
 


 double Brand; // NB names inconsistent
 Brand =Bran;
 double Bperp_reg=Bperp( x, y, z,  Bregx, Bregy, Bregz,  x0,y0,z0);
 double Bperp_ran=Bperp( x, y, z,  Branx, Brany, Branz,  x0,y0,z0);

 /* solved by error handling AWS20080307
 Bperp_reg += 1.e-9; // smaller field gives error in synchrotron.c
 Bperp_ran += 1.e-9; // smaller field gives error in synchrotron.c
 */

 if(debug==1)
 {
 cout<<endl;
 cout<<"using synchrotron routine"<<endl;
 cout<<endl;

 cout<<"regular field:"<<endl;
 cout<<"gamma= "<<gamma<<" nu="<<nu<<" Bperp_reg="<<Bperp_reg   << " Brand="<<Brand   <<endl;

 }


 synch_emissivity_total
     = synchrotron_emissivity( gamma, nu, Bperp_reg, Brand,
                         synch_emissivity_reg, synch_emissivity_par, synch_emissivity_perp, synch_emissivity_random, debug );


 if(debug==1)
 {
   cout
  <<" synch emissivity random   =  "  << synch_emissivity_random <<endl
  <<" synch emissivity reg      =  "  << synch_emissivity_reg    <<endl
  <<" synch emissivity parallel =  "  << synch_emissivity_par    <<endl
  <<" synch emissivity perp     =  "  << synch_emissivity_perp   <<endl 
  <<" synch emissivity total    =  "  << synch_emissivity_total  <<endl       
  <<endl;
 }


 /*
// this piece is redundant and in fact wrong since it double-counts the perp. component of the random field
// something like it may be useful if simulation a real random field
 if(debug==1)
 {
 cout<<"random  field:"<<endl;
 cout<<"gamma= "<<gamma<<" nu="<<nu<<" Bperp_ran="<<Bperp_ran   << " Brand="<<Brand   <<endl;
 }
//                                          V  here   V  
 synch_emissivity_total
     = synchrotron_emissivity( gamma, nu, Bperp_ran, Brand,
                         synch_emissivity_reg, synch_emissivity_par, synch_emissivity_perp, synch_emissivity_random, debug );


 if(debug==1)
 {
  cout
  <<" synch emissivity random   = "  << synch_emissivity_random <<endl
  <<" synch emissivity reg      = "  << synch_emissivity_reg    <<endl
  <<" synch emissivity parallel = "  << synch_emissivity_par    <<endl
  <<" synch emissivity perp     = "  << synch_emissivity_perp   <<endl 
  <<" synch emissivity total    = "  << synch_emissivity_total  <<endl       
  <<endl;
 }

 */ 

 return synch_emissivity_total ;

}


/////////////////////////////////////////////////////////////// version with Stokes AWS20100706

void   stokes ( double x , double y , double z , double x0, double y0, double z0,
                double Bx, double By, double Bz, double perp, double par,
                double &I, double&Q, double &U, int debug);


////////////////////////////////////////////////////////////////////////////////////////////////
double synchrotron_emissivity_B_field( double gamma, double nu, 
				       char *name, double *parameters,
				       double x, double y, double z,  
                                       int options,
                                       double x0, double y0, double z0,
                                       double     &synch_emissivity_reg, double &synch_emissivity_par, double &synch_emissivity_perp, double &synch_emissivity_random,
                                       double     &I, double     &Q, double &U,
                                       int debug=0 )

{

  
 
  double  Breg,          Bregx,         Bregy,        Bregz;
  double  Bran,          Branx,         Brany,        Branz;

  

 B_field_3D_model ( name, parameters, x,  y,  z,  options,
                 Breg,  Bregx, Bregy, Bregz,
		 Bran,  Branx, Brany, Branz,
                 debug );

 if(debug==1)
 {

 cout<<"(x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  ";
 cout<<" B model= "<<name<<"  Breg="<<  Breg <<" Bregx="<<  Bregx<<" Bregy="<< Bregy<<" Bregz="<< Bregz
                          <<" Bran="<<  Bran <<" Branx="<<  Branx<<" Brany="<< Brany<<" Branz="<< Branz<<endl;

 cout<<"Bperp for regular field ="<<Bperp( x, y, z,  Bregx, Bregy, Bregz,  x0,y0,z0)<<endl;
 cout<<"Bperp for random  field ="<<Bperp( x, y, z,  Branx, Brany, Branz,  x0,y0,z0)<<endl;


 }

 // use synchrotron routine

 double synch_emissivity_total;
 


 double Brand; // NB names inconsistent
 Brand =Bran;
 double Bperp_reg=Bperp( x, y, z,  Bregx, Bregy, Bregz,  x0,y0,z0);
 double Bperp_ran=Bperp( x, y, z,  Branx, Brany, Branz,  x0,y0,z0);

 /* solved by error handling AWS20080307
 Bperp_reg += 1.e-9; // smaller field gives error in synchrotron.c
 Bperp_ran += 1.e-9; // smaller field gives error in synchrotron.c
 */

 if(debug==1)
 {
 cout<<endl;
 cout<<"using synchrotron routine"<<endl;
 cout<<endl;

 cout<<"regular field:"<<endl;
 cout<<"gamma= "<<gamma<<" nu="<<nu<<" Bperp_reg="<<Bperp_reg   << " Brand="<<Brand   <<endl;

 }


 synch_emissivity_total
     = synchrotron_emissivity( gamma, nu, Bperp_reg, Brand,
                         synch_emissivity_reg, synch_emissivity_par, synch_emissivity_perp, synch_emissivity_random, debug );


 if(debug==1)
 {
   cout
  <<" synch emissivity random   =  "  << synch_emissivity_random <<endl
  <<" synch emissivity reg      =  "  << synch_emissivity_reg    <<endl
  <<" synch emissivity parallel =  "  << synch_emissivity_par    <<endl
  <<" synch emissivity perp     =  "  << synch_emissivity_perp   <<endl 
  <<" synch emissivity total    =  "  << synch_emissivity_total  <<endl       
  <<endl;
 }


 /*
// this piece is redundant and in fact wrong since it double-counts the perp. component of the random field
// something like it may be useful if simulation a real random field
 if(debug==1)
 {
 cout<<"random  field:"<<endl;
 cout<<"gamma= "<<gamma<<" nu="<<nu<<" Bperp_ran="<<Bperp_ran   << " Brand="<<Brand   <<endl;
 }
//                                          V  here   V  
 synch_emissivity_total
     = synchrotron_emissivity( gamma, nu, Bperp_ran, Brand,
                         synch_emissivity_reg, synch_emissivity_par, synch_emissivity_perp, synch_emissivity_random, debug );


 if(debug==1)
 {
  cout
  <<" synch emissivity random   = "  << synch_emissivity_random <<endl
  <<" synch emissivity reg      = "  << synch_emissivity_reg    <<endl
  <<" synch emissivity parallel = "  << synch_emissivity_par    <<endl
  <<" synch emissivity perp     = "  << synch_emissivity_perp   <<endl 
  <<" synch emissivity total    = "  << synch_emissivity_total  <<endl       
  <<endl;
 }

 */ 

 int stokes_debug=debug;
 

 stokes  (  x,y,z,x0,y0,z0 ,  Bregx, Bregy,  Bregz, synch_emissivity_perp, synch_emissivity_par, I,Q,U,stokes_debug);

 return synch_emissivity_total ;

}


////////////////////////////////////////////////////////////////////////////////// Stokes routines added AWS20100706



 #include<cmath>

/////////////////////////////////////////////////////////////////////////////////////////////////////////
void   stokes ( double x , double y , double z , double x0, double y0, double z0,
                double Bx, double By, double Bz, double perp, double par,
                double &I, double&Q, double &U, int debug=0)
{


  // computes Stokes parameters. The coordinate system is Galactic, with psi measured from Galactic N-S latitude meridians
 
  if(debug==1)  cout<<endl<<"========= stokes"<<endl;

  double l,b;
  double X[3],Y[3],Z[3]; // rotated coordinates: YZ is plane of sky perp to line of sight. Z is parallel to Galactic latitude pointing to Galactic pole.

  double dtr=acos(-1.0)/180.;

  // direction of point x,y,z seen from Solar position (x0,y0,z0). 
  //  +ve long=+ve y for Z=-X^Y (LH) system: see gen_synch_skymap.cc
  l=atan2(y-y0,x0-x);
  b=atan2(z-z0,sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)));


  double sinl=sin(l);
  double cosl=cos(l);
  double sinb=sin(b);
  double cosb=cos(b);

  // rotated system

  X[0] = cosb*cosl;
  X[1] = cosb*sinl;
  X[2] = sinb;


  Y[0] = -sinl;
  Y[1] =  cosl;
  Y[2] =  0.;


  Z[0] = -sinb*cosl;
  Z[1] = -sinb*sinl;
  Z[2] =  cosb;


  // check orthogonality of axes using dot products 
  double XY = X[0]*Y[0]+X[1]*Y[1]+X[2]*Y[2];
  double XZ = X[0]*Z[0]+X[1]*Z[1]+X[2]*Z[2];
  double YZ = Y[0]*Z[0]+Y[1]*Z[1]+Y[2]*Z[2];




  // B in rotated system (p=primed)
  double Bxp, Byp, Bzp;
  Bxp= Bx*X[0]+By*X[1]+Bz*X[2];  // B  . X'
  Byp= Bx*Y[0]+By*Y[1]+Bz*Y[2];  // B  . Y'
  Bzp= Bx*Z[0]+By*Z[1]+Bz*Z[2];  // B  . Z'



  double psi;

  psi=atan2(Byp,Bzp); 

  if(debug==1)
  {
  cout<<endl;
  cout << "l = "<<l<<" b= "<<b<<endl;
  cout<<"local coordinate system :"<<endl;
  cout<<"X = "; for (int i=0;i<3;i++) cout<<X[i]<<" "; cout<<endl;
  cout<<"Y = "; for (int i=0;i<3;i++) cout<<Y[i]<<" "; cout<<endl;
  cout<<"Z = "; for (int i=0;i<3;i++) cout<<Z[i]<<" "; cout<<endl;


  cout<<" X.Y="<<XY;
  cout<<" X.Z="<<XZ;
  cout<<" Y.Z="<<YZ<<endl;
  
  cout<<"Bx  ="<<Bx <<" By  ="<<By <<" Bz  ="<<Bz <<" Btot ="<<sqrt(Bx *Bx +By *By +Bz *Bz )<<endl;
  cout<<"Bxp ="<<Bxp<<" Byp ="<<Byp<<" Bzp ="<<Bzp<<" Btotp="<<sqrt(Bxp*Bxp+Byp*Byp+Bzp*Bzp)<<endl;
  }


  

  // polarized intensity = perp-par


  // as in Hammurabi
  Q =       (perp-par)*cos(2.*psi);
  U =       (perp-par)*sin(2.*psi);

  // based on definition of Stokes "Q=Ex^2-Ey^2" with B defining x-axis
  Q =       (par-perp)*cos(2.*psi);
  U =       (par-perp)*sin(2.*psi);

  
  // I here defined as polarized intensity only
  I =sqrt(Q*Q + U*U);




  if(debug==1||debug==2)
  {
  cout << "stokes: x = "<<x    <<" y= "<<y <<" z= "<<z;  
  cout <<        " x0= "<<x0   <<" y0="<<y0<<" z0="<<z0; 
  cout<<" Bx  ="<<Bx <<" By  ="<<By <<" Bz  ="<<Bz ;
  cout<<" perp="<<perp<<" par="<<par;
  cout << "  /  l = "<<l/dtr<<" b= "<<b/dtr;
  cout<<" Stokes I = "<<I<<" Q = "<<Q  <<" U = "<<U         ;     
  cout<<" psi (deg) = "<<psi/dtr;
  cout<<" polarization="<<(perp-par)/(perp+par);
  cout<<" total intensity = "<<perp+par;
  cout<<" unpol intensity = "<<2*par<<endl;
  }

  return;
}

//////////////////////////////////////////////////////////////
int stokes_test()
{
 double Bx,By,Bz;
 double perp,par;
 double x,y,z,x0,y0,z0;

 x0=8.5; y0=0;z0=0;
 perp=1.; par=.15; // roughly realistic for 75% polarized
 int debug=2;

 double I,Q,U; // output parameters

 cout<<"=============== stokes_test"<<endl;
 x=0;y=0;z=0; Bx=1; By=0;Bz=0; stokes  (  x,y,z,x0,y0,z0 ,  Bx, By,  Bz, perp, par, I,Q,U,debug);
 x=0;y=0;z=0; Bx=0; By=1;Bz=0; stokes  (  x,y,z,x0,y0,z0 ,  Bx, By,  Bz, perp, par, I,Q,U,debug);
 x=0;y=0;z=0; Bx=0; By=0;Bz=1; stokes  (  x,y,z,x0,y0,z0 ,  Bx, By,  Bz, perp, par, I,Q,U,debug);
 x=0;y=0;z=0; Bx=0; By=1;Bz=1; stokes  (  x,y,z,x0,y0,z0 ,  Bx, By,  Bz, perp, par, I,Q,U,debug);


x=x0;y=2;z=0; Bx=1; By=0;Bz=0; stokes  (  x,y,z,x0,y0,z0 ,  Bx, By,  Bz, perp, par, I,Q,U,debug);
 x=9;y=0;z=0; Bx=0; By=1;Bz=0; stokes  (  x,y,z,x0,y0,z0 ,  Bx, By,  Bz, perp, par, I,Q,U,debug);
 x=9;y=0;z=0; Bx=0; By=1;Bz=1; stokes  (  x,y,z,x0,y0,z0 ,  Bx, By,  Bz, perp, par, I,Q,U,debug);

  return 0;
}
/////////////////////////////////////////////////////////
/*
int main()
{stokes_test();return 0;}
*/
