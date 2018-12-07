#include "B_field_3D_model.h"
#include "synchrotron_emissivity.h"
#include "synchrotron_emissivity_B_field.h"

#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>

using namespace std;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double synchrotron_emissivity_aws(double gamma, double nu, double x, double y, double z,  char* name, double *parameters, int debug=0)
{

 int options;
  double x0,y0,z0;

 double synch_emissivity_total;
 double synch_emissivity_reg;
 double synch_emissivity_par;
 double synch_emissivity_perp;
 double synch_emissivity_random;


  x0=8.5; // solar position
  y0=0.;
  z0=0.;

  options=0;


 synch_emissivity_total  =
 synchrotron_emissivity_B_field(  gamma,  nu, 
				  name, parameters, x, y, z, options,
                                  x0, y0, z0,
                                  synch_emissivity_reg,synch_emissivity_par,synch_emissivity_perp,synch_emissivity_random,
			      	  debug );

 if(debug==1)
 {
 cout<<"-------------------  synchrotron_emissivity_aws: synchrotron_emissivity_B_field at ";
 cout<<"(x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<endl;
 cout<<"B field name="<<name<<"  ";

 cout
  <<" gamma= "<<gamma<<" nu="<<nu << " Hz" <<endl
  <<" synch emissivity random        =  "  << synch_emissivity_random <<endl
  <<" synch emissivity regular       =  "  << synch_emissivity_reg    <<endl
  <<" synch emissivity  pol parallel =   " << synch_emissivity_par    <<endl
  <<" synch emissivity  pol perp     =   " << synch_emissivity_perp   <<endl
  <<" synch emissivity total         =  "  << synch_emissivity_total  <<endl
  <<endl;

}



  return synch_emissivity_total;
}


//////////////////////////////////////////////////////////////////////////////////////////////////

double synchrotron_emissivity_aws(double gamma, double nu, double R, double z, char* name, double *parameters, int debug=0)

{

 double synch_emissivity_total;
 double x,y;
 double phi;
 double pi=acos(-1.0);
 phi=pi/4.;             // arbitrary angle, should really take azimuthal average
 x=R*cos(phi); 
 y=R*sin(phi);
 synch_emissivity_total  =  synchrotron_emissivity_aws(gamma,nu,x,y,z,  name, parameters, debug);

 if(debug==1)
 {
 cout<<"-------------------  synchrotron_emissivity_aws: synchrotron_emissivity_B_field at ";
 cout<<"(R, x, y, z) = ("<<R<<","<<x<<", "<<y<<", "<<z<<")  "
     <<" synch emissivity total         =  "  << synch_emissivity_total  <<endl;
 }

 return synch_emissivity_total;
}



// Stokes parameters version AWS20100707
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double synchrotron_emissivity_aws(double gamma, double nu, double x, double y, double z,  char* name, double *parameters, double &I, double &Q, double &U, int debug=0)
{

 int options;
  double x0,y0,z0;

 double synch_emissivity_total;
 double synch_emissivity_reg;
 double synch_emissivity_par;
 double synch_emissivity_perp;
 double synch_emissivity_random;


  x0=8.5; // solar position
  y0=0.;
  z0=0.;

  options=0;


 synch_emissivity_total  =
 synchrotron_emissivity_B_field(  gamma,  nu, 
				  name, parameters, x, y, z, options,
                                  x0, y0, z0,
                                  synch_emissivity_reg,synch_emissivity_par,synch_emissivity_perp,synch_emissivity_random,
                                  I, Q, U,
			      	  debug );

 if(debug==1)
 {
 cout<<"-------------------  synchrotron_emissivity_aws: synchrotron_emissivity_B_field at ";
 cout<<"(x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<endl;
 cout<<"B field name="<<name<<"  ";

 cout
  <<" gamma= "<<gamma<<" nu="<<nu << " Hz" <<endl
  <<" synch emissivity random        =  "  << synch_emissivity_random <<endl
  <<" synch emissivity regular       =  "  << synch_emissivity_reg    <<endl
  <<" synch emissivity  pol parallel =   " << synch_emissivity_par    <<endl
  <<" synch emissivity  pol perp     =   " << synch_emissivity_perp   <<endl
  <<" synch emissivity Stokes I pol  =  "  << I                       <<endl
  <<" synch emissivity Q             =  "  << Q                       <<endl
  <<" synch emissivity U             =  "  << U                       <<endl
  <<endl;

}



  return synch_emissivity_total;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

double synchrotron_emissivity_aws(double gamma, double nu, double R, double z, char* name, double *parameters, double &I, double &Q, double &U, int debug=0)

{

 double synch_emissivity_total;
 double x,y;
 double phi;
 double pi=acos(-1.0);
 phi=pi/4.;             // arbitrary angle, should really take azimuthal average
 x=R*cos(phi); 
 y=R*sin(phi);
 synch_emissivity_total  =  synchrotron_emissivity_aws(gamma,nu,x,y,z,  name, parameters, I, Q, U,  debug);

 if(debug==1)
 {
 cout<<"-------------------  synchrotron_emissivity_aws: synchrotron_emissivity_B_field at ";
 cout<<"(R, x, y, z) = ("<<R<<","<<x<<", "<<y<<", "<<z<<")  "
     <<" synch emissivity total         =  "  << synch_emissivity_total  <<endl;
 }

 return synch_emissivity_total;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// following routines are obsolete
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double synchrotron_emissivity_aws(double gamma, double nu, double x, double y, double z,  int galdef_B_field_model, int debug=0)
{

 

  char name[20];
  double parameters[20];
 
  int options;
  double x0,y0,z0;



 

 double synch_emissivity_total;
 double synch_emissivity_reg;
 double synch_emissivity_par;
 double synch_emissivity_perp;
 double synch_emissivity_random;

 
 // extract information from galdef parameter, here provided as an argument
 // to make the package independent of galprop

 


 strcpy(name,"undefined");   //AWS20080310
 
 if  (galdef_B_field_model == 1)
  strcpy(name,"test");

 if  (galdef_B_field_model == 2)
  strcpy(name,"circular");

 if  (galdef_B_field_model == 3)
  {strcpy(name,"circular2"); parameters[0]=1e-6; parameters[1]=2e-6;} // parameters: Breg, Bran

 if  (galdef_B_field_model == 4)
  {strcpy(name,"spiral");    parameters[0]=1e-6; parameters[1]=20.; parameters[2]=2e-6;} // parameters: Breg, pitch_angle (deg), Bran

 /* model invalid
 if  (galdef_B_field_model == 5)
  {strcpy(name,"wmap_page");    parameters[0]=3.e-6; parameters[1]=25.;parameters[2]=35.;parameters[3]=0.9;parameters[4]=0.;} 
                             // parameters: Breg, chi0=z dependence,psi0=opening angle of spiral arms; psi1=radial dependence of opening angle,  Bran
			     */

 if  (galdef_B_field_model == 5)// Han 1994 model and parameter values (except z-dependence from  Miville-Deschenes 2008, undefined in Han)
  {strcpy(name,"han"    );    parameters[0]=1.8e-6; parameters[1]=8.0;   parameters[2]=-8.2;      parameters[3]=11.9; parameters[4]=0.;} 
                           // Breg (G)              chi0=z dependence      p=pitch angle (deg);   ro (kpc)            Bran (G)

 if  (galdef_B_field_model == 6)// Han 1994 model, parameters values from Miville-Deschenes 2008
  {strcpy(name,"han"    );    parameters[0]=3.0e-6; parameters[1]=8.0;   parameters[2]=-8.5;      parameters[3]=11.0; parameters[4]=0.;} 
                           // Breg (G)              chi0=z dependence    p=pitch angle (deg);   ro (kpc)            Bran (G)

 if  (galdef_B_field_model == 7)// Han 1994 model,  parameters values from Miville-Deschenes 2008, plus their random field = 0.57Breg                                 
  {strcpy(name,"han"    );    parameters[0]=3.0e-6; parameters[1]=8.0;   parameters[2]=-8.5;      parameters[3]=11.0; parameters[4]=1.71e-6;} 
                           // Breg (G)              chi0=z dependence      p=pitch angle (deg);   ro (kpc)            Bran (G)


 if  (galdef_B_field_model == 8)//  Tinyakov & Tkachev 2002    model, parameters value as in their paper              AWS20080311
                                //  parameters arranged as for Han model for easy reading (1 and 3 not used)
  {strcpy(name,"tinyakov"    );    parameters[0]=1.4e-6; parameters[1]=0.0;   parameters[2]=-8.0;      parameters[3]= 0.0; parameters[4]=0.;} 
                                // Breg (G)                                p=pitch angle (deg);                       Bran (G)

 if  (galdef_B_field_model == 9)// Han 1994 model,  parameters values from Miville-Deschenes 2008, increasing Breg to fit 408 MHz with regular only                 
  {strcpy(name,"han"    );    parameters[0]=14.0e-6;parameters[1]=8.0;   parameters[2]=-8.5;      parameters[3]=11.0; parameters[4]=0.00e-6;} 
                           // Breg (G)              chi0=z dependence      p=pitch angle (deg);   ro (kpc)            Bran (G)



 // original galprop exponential model

 if  (galdef_B_field_model > 1000) // original model from galprop B_field_model.cc, transferring model code via parameters[0]
   {strcpy(name,"galprop_original");    parameters[0]= galdef_B_field_model;}
  


 if(strcmp(name,"undefined")==0) //AWS20080310
 {
    cout<<"synchrotron_emissivity_aws: B_field_model="<<galdef_B_field_model<<" undefined !"<<endl;
    exit(1);
 }



  x0=8.5; // solar position
  y0=0.;
  z0=0.;

  options=0;

 



 synch_emissivity_total  =
 synchrotron_emissivity_B_field(  gamma,  nu, 
				  name, parameters, x, y, z, options,
                                  x0, y0, z0,
                                  synch_emissivity_reg,synch_emissivity_par,synch_emissivity_perp,synch_emissivity_random,
			      	  debug );

 if(debug==1)
 {
 cout<<"-------------------  synchrotron_emissivity_aws: synchrotron_emissivity_B_field at ";
 cout<<"(x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<endl;
 cout<<"B field name="<<name<<"  ";

 cout
  <<" gamma= "<<gamma<<" nu="<<nu << " Hz" <<endl
  <<" synch emissivity random        =  "  << synch_emissivity_random <<endl
  <<" synch emissivity regular       =  "  << synch_emissivity_reg    <<endl
  <<" synch emissivity  pol parallel =   " << synch_emissivity_par    <<endl
  <<" synch emissivity  pol perp     =   " << synch_emissivity_perp   <<endl
  <<" synch emissivity total         =  "  << synch_emissivity_total  <<endl
  <<endl;

 }
 


  return synch_emissivity_total;
}


//////////////////////////////////////////////////////////////////////////////////////////////////

double synchrotron_emissivity_aws(double gamma, double nu, double R, double z, int galdef_B_field_model, int debug=0)

{

 double synch_emissivity_total;
 double x,y;
 double phi;
 double pi=acos(-1.0);
 phi=pi/4.;
 x=R*cos(phi); 
 y=R*sin(phi);
 synch_emissivity_total  =  synchrotron_emissivity_aws(gamma,nu,x,y,z,  galdef_B_field_model,debug);

 if(debug==1)
 {
 cout<<"-------------------  synchrotron_emissivity_aws: synchrotron_emissivity_B_field at ";
 cout<<"(R, x, y, z) = ("<<R<<","<<x<<", "<<y<<", "<<z<<")  "
     <<" synch emissivity total         =  "  << synch_emissivity_total  <<endl;
 }

 return synch_emissivity_total;
}
