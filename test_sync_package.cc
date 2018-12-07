using namespace std;
#include<iostream>


#include "B_field_3D_model.h"
#include "synchrotron_emissivity.h"
#include "synchrotron_emissivity_B_field.h"
#include "synchrotron_emissivity_aws.h"
#include <cstring>

extern "C" double synchrotron_(double*,double*,double*); // galprop fortran version for tests

////////////////////////////////////////////////////////////////////////
//                      test routine
////////////////////////////////////////////////////////////////////////
int test_sync_package()
{

  char name[20];
  double parameters[20];
  double x,y,z;
  int options;
  double x0,y0,z0;



 double gamma,nu;

 double synch_emissivity_total;
 double synch_emissivity_reg;
 double synch_emissivity_par;
 double synch_emissivity_perp;
 double synch_emissivity_random;

 int debug=1;

 double E_electron=1.0e4; // 10 GeV electrons = 1e4 MeV

 gamma=E_electron   /.511; 
 nu=1.e9; // 1000 MHz

  strcpy(name,"test");
  strcpy(name,"circular");
  strcpy(name,"circular2"); parameters[0]=1e-6; parameters[1]=2e-6; // parameters: Breg, Bran
  strcpy(name,"spiral");    parameters[0]=1e-6; parameters[1]=20.; parameters[2]=2e-6; // parameters: Breg, pitch_angle (deg), Bran
  strcpy(name,"galprop_original");    parameters[0]= 50100010; // as B_field_model.cc, galdef.B_field_model coded as parameter, avoid leading zeros, taken as hex!

  x=8;
  y=1;
  z=1;

  x0=8.5; // solar position
  y0=0.;
  z0=0.;

  options=0;

 for (x= -10.; x<+10.; x+=2)
 for (y= -10.; y<+10.; y+=2)
 for (z=  -1; z< +1.1;  z+=1.)
 {



 synch_emissivity_total  =
 synchrotron_emissivity_B_field(  gamma,  nu, 
				  name, parameters, x, y, z, options,
                                  x0, y0, z0,
                                  synch_emissivity_reg,synch_emissivity_par,synch_emissivity_perp,synch_emissivity_random,
			      	  debug );

 cout<<"------------------- testing synchrotron_emissivity_B_field at ";
 cout<<"(x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<endl;
 cout<<"B field name="<<name<<"  ";

 cout
  <<" electron energy= "<<E_electron<<" MeV  gamma= "<<gamma<<" nu="<<nu << " Hz" <<endl
  <<" synch emissivity random        =  "  << synch_emissivity_random <<endl
  <<" synch emissivity regular       =  "  << synch_emissivity_reg    <<endl
  <<" synch emissivity  pol parallel =   " << synch_emissivity_par    <<endl
  <<" synch emissivity  pol perp     =   " << synch_emissivity_perp   <<endl
  <<" synch emissivity total         =  "  << synch_emissivity_total  <<endl
  <<endl;


  } // for x y z





 ////////////////////////////////////////////////////////////////////////////////
 // interface to galprop
 ////////////////////////////////////////////////////////////////////////////////


 
 int galdef_B_field_model;

 debug=1;

 for (int icase=1;icase<=2;icase++)
 {




 if(icase==1) galdef_B_field_model= 50100010; // galprop original B_field_model.cc, tagged as negative NB avoid leading zeros, taken as hex!
 if(icase==2) galdef_B_field_model= 5;
 if(icase==3) galdef_B_field_model= 6;
 if(icase==4) galdef_B_field_model= 7;

 for (x= -10.; x<+10.; x+=2)
 for (y= -10.; y<+10.; y+=2)
 for (z=  -2; z< +2.1;  z+=1.)
 {



   synch_emissivity_total  = synchrotron_emissivity_aws(gamma,nu,x,y,z,galdef_B_field_model);

 cout<<"------------------- testing synchrotron_emissivity_aws at ";
 cout<<"(x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  " <<" synch emissivity total         =  "  << synch_emissivity_total  <<endl;
  

 }

 double R;
 

 for (R= 1.; R<+10.; R+=2.)
 for (z=  -2; z< +2.1;  z+=1.)
 {



   synch_emissivity_total  = synchrotron_emissivity_aws(gamma,nu,R,z,galdef_B_field_model,debug);

 cout<<"------------------- testing synchrotron_emissivity_aws at ";
 cout<<"(R, z) = ("<<R<<", "<<z<<")  " <<" synch emissivity total         =  "  << synch_emissivity_total  <<endl;
  

 }

}//icase



 cout<<endl<<"----------------------- comparing new routines with old fortran"<<endl;

 // galprop debug output:
 //gamma nu r z 1957.95 1e+09 8 0 galaxy.B_field (old)=4.0402e-10  sync_per_electron(old)= 3.014e-34 (new)= 4.56011e-34

 gamma=1957.95;
 nu=1e9;
 E_electron=gamma*.511;
 cout  <<" electron energy= "<<E_electron<<" MeV  gamma= "<<gamma<<" nu="<<nu << " Hz" <<endl;

 double Brand=4.0402e-6; // Gauss to compare with galprop debug
 //       Brand=4.0000e-6; // Gauss for exact comparison of routines

 double Brand_Tesla=Brand/1e4; // since fortran routine uses Tesla

 // compare with fortran version from galprop
 double synch_emissivity_fort=synchrotron_( &gamma, &nu, &Brand);





 x=8.0; // corresponds to grid point
 y=0.;
 z=0.;
 galdef_B_field_model= 40500030; // 4 microGauss at r=ro z=0
 debug=1;
 synch_emissivity_total  = synchrotron_emissivity_aws(gamma,nu,x,y,z,galdef_B_field_model,debug);

 cout<<endl;
 cout<<"emissivity fortran           ="<< synch_emissivity_fort<<endl;
 cout<<"synchrotron_emissivity_aws   ="<< synch_emissivity_total<<endl;
 cout<<"difference                   ="<< synch_emissivity_total- synch_emissivity_fort <<endl;
 cout<<"frac difference              ="<<( synch_emissivity_total- synch_emissivity_fort)/synch_emissivity_total    <<endl;
 cout<<endl<<endl;

 // test 2D for R=x
 synch_emissivity_total  = synchrotron_emissivity_aws(gamma,nu,x,  z,galdef_B_field_model,debug);
 cout<<"synchrotron_emissivity_aws Rz="<< synch_emissivity_total<<endl;


  strcpy(name,"galprop_original");    parameters[0]= 40500030; 
 synch_emissivity_total  =
 synchrotron_emissivity_B_field(  gamma,  nu, 
				  name, parameters, x, y, z, options,
                                  x0, y0, z0,
                                  synch_emissivity_reg,synch_emissivity_par,synch_emissivity_perp,synch_emissivity_random,
			      	  debug );

 cout<<"------------------- testing synchrotron_emissivity_B_field at ";
 cout<<"(x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<endl;
 cout<<"B field name="<<name<<"  ";

 cout
  <<" electron energy= "<<E_electron<<" MeV  gamma= "<<gamma<<" nu="<<nu << " Hz" <<endl
  <<" synch emissivity random        =  "  << synch_emissivity_random <<endl
  <<" synch emissivity regular       =  "  << synch_emissivity_reg    <<endl
  <<" synch emissivity  pol parallel =   " << synch_emissivity_par    <<endl
  <<" synch emissivity  pol perp     =   " << synch_emissivity_perp   <<endl
  <<" synch emissivity total         =  "  << synch_emissivity_total  <<endl
  <<endl;

  return 0;
}
