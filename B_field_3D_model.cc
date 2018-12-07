#include<iostream>
#include <cmath>
#include <cstring> // AWS20080711 r194: for compatibility with gcc4.3 

using namespace std;

void B_field_3D_model
(char* name, double *parameters,
 double x, double y, double z,
 int options,
 double &Breg,  double &Bregx, double &Bregy, double &Bregz,
 double &Bran,  double &Branx, double &Brany, double &Branz,
 int debug=0 )
{

  // General purpose interface for Galactic magnetic field models
  // including regular and random components

  // input  arguments

  // name: name of model
  // parameters: array of model parameters for model with this name
  // options:    for addional control
  // x,y,z  Galactic coordinates of point, kpc  (Sun at x=Rsun, y=z=0)


  // output arguments
  // magnetic field strengths in Gauss
  // Breg   regular field total
  // Bregx  regular field x-component
  // Bregy  regular field y-component
  // Bregz  regular field z-component

  // Bran   random  field total
  // Branx  random  field x-component (for a random direction)
  // Brany  random  field y-component (for a random direction) 
  // Branz  random  field z-component (for a random direction)

  // local variables
  double theta,phi; // spherical angles
  double dtr=acos(-1.)/180.; // degrees to radians
  double pi =acos(-1.);

  if(debug==1)cout<<"B_field_3D_model"<<endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  if(strcmp(name,"test")==0)
  {
    // uniform regular field strength and direction, and uniform random field
    Breg=6e-6;
    theta=89*dtr;
    phi= 100*dtr;
    Bregx=Breg*sin(theta)*cos(phi);
    Bregy=Breg*sin(theta)*sin(phi);
    Bregz=Breg*cos(theta);         

    Bran=5e-6;
    theta=30*dtr;
    phi= 200*dtr;
    Branx=Bran*sin(theta)*cos(phi);
    Brany=Bran*sin(theta)*sin(phi);
    Branz=Bran*cos(theta);  
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(strcmp(name,"circular")==0)
  {
    // uniform regular field strength parallel to plane, circular geometry, and uniform random field
    Breg=6e-6;
    theta=90.*dtr;  // zenith angle for in-plane field
    phi= atan2(y,x)+ pi/2; // this makes field direction circular
    Bregx=Breg*sin(theta)*cos(phi);
    Bregy=Breg*sin(theta)*sin(phi);
    Bregz=Breg*cos(theta);         

    Bran=5e-6;
    //   theta=30*dtr;  in future should assign a random direction
    //    phi= 200*dtr;
    Branx=Bran*sin(theta)*cos(phi);
    Brany=Bran*sin(theta)*sin(phi);
    Branz=Bran*cos(theta);  
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(strcmp(name,"circular2")==0)
  {
    // uniform regular field strength parallel to plane, circular geometry, and uniform random field
    // parameterized as below
    Breg=parameters[0];
    theta=90.*dtr;  // zenith angle for in-plane field
    phi= atan2(y,x)+ pi/2; // this makes field direction circular
    Bregx=Breg*sin(theta)*cos(phi);
    Bregy=Breg*sin(theta)*sin(phi);
    Bregz=Breg*cos(theta);         

    Bran=parameters[1];
    //   theta=30*dtr;  in future should assign a random direction
    //    phi= 200*dtr;
    Branx=Bran*sin(theta)*cos(phi);
    Brany=Bran*sin(theta)*sin(phi);
    Branz=Bran*cos(theta);  
  }



 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(strcmp(name,"spiral")==0)
  {
    // uniform regular field strength parallel to plane, spiral geometry, and uniform random field
    // parameterized as below
    Breg=parameters[0];
    double pitch_angle=parameters[1]*dtr; 
    theta=90.*dtr;  // zenith angle for in-plane field
    phi= atan2(y,x)+ pi/2 + pitch_angle; // pitch angle relative to circle
    Bregx=Breg*sin(theta)*cos(phi);
    Bregy=Breg*sin(theta)*sin(phi);
    Bregz=Breg*cos(theta);         

    Bran=parameters[2];
    //   theta=30*dtr;  in future should assign a random direction
    //    phi= 200*dtr;
    Branx=Bran*sin(theta)*cos(phi);
    Brany=Bran*sin(theta)*sin(phi);
    Branz=Bran*cos(theta);  
  }






  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  if(strcmp(name,"galprop_original")==0)
  {
    // same model as in galprop B_field_model.cc 
    // Bo rscale zscale encoded in 9-digit number: BBBrrrzzz in units of 0.1
    // e.g. 123456789 : Bo= 12.3E-10 Tesla  rscale=45.6 kpc zscale=78.9 kpc



    double Bo,rscale,zscale;
    double ro=8.5;
    double b_field;

    

    double r=sqrt(x*x+y*y);


    /* this was original test case, no longer used
    int model=int(parameters[0]+.001);// using this to represent galplotdef.B_field_model; parameters are double so convert to integer for galprop B_field_model.cc case
                                      // NB avoid leading zeros when defining B_field_model using integer (it is interpreted as hex!)
    // following 7-line  code segment is copied from  galprop B_field_model.cc, gives total B random in Tesla
   if (model > 1000)
   {
      Bo=           (model/1000000)                 * 0.1 *1.0e-10;
      rscale=(model-(model/1000000)*1000000 )/1000  * 0.1         ;
      zscale=(model%1000)                           * 0.1         ;
      b_field=Bo *exp(-(r-ro)/rscale) * exp(-fabs(z)/zscale);

   if(debug==1)
   cout<<"galprop original model="<<model<<" Bo="<<Bo<<" rscale="<<rscale<<" zscale="<<zscale
    <<" (x, y, z) = ("<<x<<", "<<y<<", "<<z<<")" << " r=" <<r 
    <<" b_field="<<b_field<<endl;

   }
    */

   // now take parameters as specified for 3D model in galdef file           AWS20080314
   Bo   = parameters[0]; // Gauss
   rscale=parameters[1]; // kpc
   zscale=parameters[2]; // kpc
   b_field=Bo *exp(-(r-ro)/rscale) * exp(-fabs(z)/zscale);

 

   // the regular field is zero
   Breg =0.0;
   Bregx=0.0;
   Bregy=0.0;
   Bregz=0.0;

   Bran=b_field;
   Branx=Bran;  // put all field in x-direction (has no significance)
   Brany=0.0;
   Branz=0.0;

   theta=0.0;
   phi  =0.0;

  } //galprop_original

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

    // WMAP regular field strength and direction (x,y,z); right handed system; Sun in x; B in Page et al. 2007 adapted to this coordinate system
    // from Elena Orlando: but she says this formulation is invalid: keep for reference
if(strcmp(name,"wmap_page")==0)
  {
  
       
  
   double chi0,chi,psi0,psi,psi1;
   double ro=8;

    Breg=parameters[0];
    phi=atan2(y,x);
    chi0=parameters[1]*dtr;//parameter[1]=25 in Page et al 2007; z dependence
    chi=chi0*tanh(z);
    psi0=parameters[2]*dtr;//perameter[2]=35; opening angle of spiral arms
    psi1=parameters[3]*dtr;//parameter[3]=0.9 radial dependence of opening angle
    psi=psi0+psi1*log(pow(x*x+y*y,0.5)/ro);
    Bregx=Breg*cos(chi)*sin(psi+phi);
    Bregy=-Breg*cos(chi)*cos(psi+phi);
    Bregz=Breg*sin(chi);  

   Bran=parameters[4]; 
   Branx=Bran;  // put all field in x-direction (has no significance)
   Brany=0.0;
   Branz=0.0;

   theta=0.0;


}

/////////////////////////////////////////////////////////////////////////////////////////////

// Han, J.L. 1994 A&A 288,759-772 formulation = WMAP Miville-Deschenes et al.  2008 arXiv:0802.3345
// select the parameters. bi-symmetric spiral; Sun on x-axis  , theta=clockwise 
// in the direction of positive y; reference system of galprop.
// A constant random component can be added via parameters[4]
// from Elena Orlando

if(strcmp(name,"han")==0)
 {

 
   //    Breg=parameters[0];           //=1.8 const in Han; =3 const in WMAP M.        //AWS20080331 Breg is output parameter!

    double B0=parameters[0];           //=1.8 const in Han; =3 const in WMAP M.        //AWS20080331
    theta=atan2(y,x);             // NB theta defined differently from above 
    double chi0=parameters[1]*dtr;// not defined in Han; parameter[1]=8 in WMAP Miville
    double chi=chi0*tanh(z);      //  z dependence, not defined in Han;
    double p=parameters[2]*dtr;   // pitch angle =-8.2 in Han;-8.5 in WMAP M.
    double psi=1./(tan(p));       // radial dependence of opening angle
    double ro=parameters[3];      // =11.9 in Han, =11 in WMAP M.     
    
    //    Bregx=Breg*cos(theta-psi*log(pow(x*x+y*y,0.5)/ro))*sin(p-theta)*cos(chi);    //AWS20080331
    //    Bregy=Breg*cos(theta-psi*log(pow(x*x+y*y,0.5)/ro))*cos(p-theta)*cos(chi);    //AWS20080331
    Bregx=B0*cos(theta-psi*log(pow(x*x+y*y,0.5)/ro))*sin(p-theta)*cos(chi);            //AWS20080331
    Bregy=B0*cos(theta-psi*log(pow(x*x+y*y,0.5)/ro))*cos(p-theta)*cos(chi);            //AWS20080331

    Bregz=0.0; // Han (and Miville-Deschenes ??) have zero z-component
    Bregz=parameters[5];                                                               //AWS20080711

    //Bregz=Breg*sin(chi); in case we need a z-component in future

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy +  Bregz*Bregz);                               //AWS20080331

    // constant random component can be added 
    Bran=parameters[4]; 
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

   
}
/////////////////////////////////////////////////////////////////////////////////////////////
// same as model "han" but random field is scaled to total regular at each point using parameters[4]
// contributed by AWS

  if(strcmp(name,"han_scaled_ran")==0) //AWS20080319
 {



    double B0=parameters[0];           //=1.8 const in Han; =3 const in WMAP M.        //AWS20080331
    theta=atan2(y,x);             // NB theta defined differently from above 
    double chi0=parameters[1]*dtr;// not defined in Han; parameter[1]=8 in WMAP Miville
    double chi=chi0*tanh(z/parameters[9]);      //  z dependence, not defined in Han;  //EO20081120
    //   double chi=tanh(z/parameters[9]);      //  z dependence, not defined in Han; error in Miville formula ??             //EO20081120
    double p=parameters[2]*dtr;   // pitch angle =-8.2 in Han;-8.5 in WMAP M.
    double psi=1./(tan(p));       // radial dependence of opening angle
    double ro=parameters[3];      // =11.9 in Han, =11 in WMAP M.     
    
  
    Bregx=B0*cos(theta-psi*log(pow(x*x+y*y,0.5)/ro))*sin(p-theta)*cos(chi);            //AWS20080331
    Bregy=B0*cos(theta-psi*log(pow(x*x+y*y,0.5)/ro))*cos(p-theta)*cos(chi);            //AWS20080331

    Bregz=0.0; // Han (and Miville-Deschenes ??) have zero z-component
    Bregz=parameters[5];                                                               //AWS20080711

    //Bregz=Breg*sin(chi); in case we need a z-component in future

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy +  Bregz*Bregz);                               //AWS20080331

    // random component proportional to total regular
    //    Bran=parameters[4]*sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);            //AWS20080331

    Bran=parameters[4]*Breg;                                                           //AWS20080331
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

   
}

/////////////////////////////////////////////////////////////////////////////////////////////
// same as model "han" but random field constant but with z dependence


  if(strcmp(name,"han_ran")==0) //EO20081127
 {



    double B0=parameters[0];           
    theta=atan2(y,x);             // NB theta defined differently from above 
    double chi0=parameters[1]*dtr;// not defined in Han; parameter[1]=8 in WMAP Miville
    double chi=chi0*tanh(z/parameters[9]);      //  z dependence, not defined in Han;  //EO20081120
    //   double chi=tanh(z/parameters[9]);      //  z dependence, not defined in Han; error in Miville formula ??             //EO20081120
    double p=parameters[2]*dtr;   // pitch angle =-8.2 in Han;-8.5 in WMAP M.
    double psi=1./(tan(p));       // radial dependence of opening angle
    double ro=parameters[3];      // =11.9 in Han, =11 in WMAP M.     
    
  
    Bregx=B0*cos(theta-psi*log(pow(x*x+y*y,0.5)/ro))*sin(p-theta)*cos(chi);            //AWS20080331
    Bregy=B0*cos(theta-psi*log(pow(x*x+y*y,0.5)/ro))*cos(p-theta)*cos(chi);            //AWS20080331

    Bregz=0.0; // Han (and Miville-Deschenes ??) have zero z-component
    Bregz=parameters[5];                                                              

    //Bregz=Breg*sin(chi); in case we need a z-component in future

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy +  Bregz*Bregz);                            

  

    Bran=parameters[4]*cos(chi);                                                         
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

   
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
// Han with exponential z dependence, new model
// from Elena Orlando

 if(strcmp(name,"Han_exp_scaled")==0)           
 {
    double B0=parameters[0];
                                                         
   double  ro=8.5;
   theta=atan2(y,x);
   double p=parameters[2]*dtr;//=-8 pitch angle
   double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   //   double z0=1.5;                                            
   double z0=parameters[6];                         
   //   double fz=(z/fabs(z))*exp(-fabs(z/z0));   
   double fz=exp(-fabs(z/z0));//                                    
   if(z<0.) fz=-fz; // to avoid problem at z=0
   double  psi=1./(tan(p));//radial dependence of opening angle   
   double r=pow(x*x+y*y,0.5);


    Bregx=B0  *cos(theta-psi*log(    r     /r0))*sin(p-theta)*fz;                        
    Bregy=B0  *cos(theta-psi*log(    r     /r0))*cos(p-theta)*fz;
    Bregz=0.;
   
    Bregz=parameters[5];                                                            


    Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                                 

    // random component can be added for testing
    Bran=parameters[4]*Breg;
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

 }

//////////////////////////////////////////////////////////////////////////////////////////////////
// Han with exponential z dependence, new model   EO20081112
// from Elena Orlando

 if(strcmp(name,"Han_exp")==0)
 {
    double B0=parameters[0];

   double  ro=8.5;
   theta=atan2(y,x);
   double p=parameters[2]*dtr;//=-8 pitch angle
   double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   //   double z0=1.5;
   double z0=parameters[6];
   //   double fz=(z/fabs(z))*exp(-fabs(z/z0));
   double fz=exp(-fabs(z/z0));//
   if(z<0.) fz=-fz; // to avoid problem at z=0
   double  psi=1./(tan(p));//radial dependence of opening angle
   double r=pow(x*x+y*y,0.5);

    Bregx=B0  *cos(theta-psi*log(    r     /r0))*sin(p-theta)*fz;
    Bregy=B0  *cos(theta-psi*log(    r     /r0))*cos(p-theta)*fz;
    Bregz=0.;

    Bregz=parameters[5];

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);

    // random component can be added for testing
    Bran=parameters[4];
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

 }
/////////////////////////////////////////////////////////////////////////////////////////////////

// Han with exponential z dependence, new model   EO20081112
// from Elena Orlando

 if(strcmp(name,"Han_ran_exp")==0)
 {
    double B0=parameters[0];

   double  ro=8.5;
   theta=atan2(y,x);
   double p=parameters[2]*dtr;//=-8 pitch angle
   double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   //   double z0=1.5;
   double z0=parameters[6];
   //   double fz=(z/fabs(z))*exp(-fabs(z/z0));
   double fz=exp(-fabs(z/z0));//
   if(z<0.) fz=-fz; // to avoid problem at z=0
   double  psi=1./(tan(p));//radial dependence of opening angle
   double r=pow(x*x+y*y,0.5);

    Bregx=B0  *cos(theta-psi*log(    r     /r0))*sin(p-theta)*fz;
    Bregy=B0  *cos(theta-psi*log(    r     /r0))*cos(p-theta)*fz;
    Bregz=0.;

    Bregz=parameters[5];

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);

    // random component can be added for testing
    Bran=parameters[4];
    Branx=Bran*fz;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

 }


//////////////////////////////////////////////////////////////////////////////////////////////////

// Tinyakov & Tkachev 2002,Astroparticle Physics 18, 165  BSS-A model, radial dependence of Bo, no Bz
// from Elena Orlando

 if(strcmp(name,"tinyakov")==0)            //AWS20080311
 {
   double B0;                                                                              //AWS20080331
   double  ro=8.5;
   theta=atan2(y,x);
   double p=parameters[2]*dtr;//=-8 pitch angle
   double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   //   double z0=1.5;                                               AWS20081107
   double z0=parameters[6];                          //              AWS20081107
   //   double fz=(z/fabs(z))*exp(-fabs(z/z0));      //              AWS20080609
   double fz=exp(-fabs(z/z0));//                                     AWS20080609
   if(z<0.) fz=-fz; // to avoid problem at z=0
   double  psi=1./(tan(p));//radial dependence of opening angle    //AWS20080609
   double r=pow(x*x+y*y,0.5);

   double rbreak = parameters[8];                                  //EO20081119

    if(r<=rbreak) B0  =parameters[0]*(ro/rbreak);//parameters[0]=1.4
    if(r> rbreak) B0  =parameters[0]*(ro/r);

    Bregx=B0  *cos(theta-psi*log(    r     /r0))*sin(p-theta)*fz;                          //AWS20080311 20081107
    Bregy=B0  *cos(theta-psi*log(    r     /r0))*cos(p-theta)*fz; //correction AWS20080320   AWS20080311 20081107
    Bregz=0.;
    Bregz=parameters[5];                                                               //AWS20080711

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy +  Bregz*Bregz);                                   //AWS20080331

    // random component can be added for testing
    Bran=parameters[4]; 
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

 }





  if(debug==1)
 cout<<"B model= "<<name
     <<" (x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<"(theta,phi)= ("<<theta<<","<<phi<<")"
     <<"  Breg="<< Breg <<" Bregx="<<  Bregx<<" Bregy="<< Bregy<<" Bregz="<< Bregz
                         <<" Bran="<< Bran <<" Branx="<<  Branx<<" Brany="<< Brany<<" Branz="<< Branz<<endl;



//////////////////////////////////////////////////////////////////////////////////////////////////

// Tinyakov & Tkachev 2002,Astroparticle Physics 18, 165  BSS-A model, radial dependence of Bo, no Bz
// from Elena Orlando
// same as model "tinyakov" but random field is scaled to total regular at each point using parameters[4]

 if(strcmp(name,"tinyakov_scaled_ran")==0)            //AWS20080311
 {
   double B0;                                                                              //AWS20080331
   double  ro=8.5;
   theta=atan2(y,x);
   double p=parameters[2]*dtr;//=-8 pitch angle
   double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   //   double z0=1.5;                                               AWS20081107
   double z0=parameters[6];                          //              AWS20081107
   //   double fz=(z/fabs(z))*exp(-fabs(z/z0));      //              AWS20080609
   double fz=exp(-fabs(z/z0));//                                     AWS20080609
   if(z<0.) fz=-fz; // to avoid problem at z=0
   double  psi=1./(tan(p));//radial dependence of opening angle    //AWS20080609
   double r=pow(x*x+y*y,0.5);

   double rbreak = parameters[8];                                  //EO20081119

    if(r<=rbreak) B0  =parameters[0]*(ro/rbreak);//parameters[0]=1.4
    if(r> rbreak) B0  =parameters[0]*(ro/r);

    Bregx=B0  *cos(theta-psi*log(    r     /r0))*sin(p-theta)*fz;                          //AWS20080311  20081107
    Bregy=B0  *cos(theta-psi*log(    r     /r0))*cos(p-theta)*fz; //correction AWS20080320   AWS20080311  20081107
    Bregz=0.;
    Bregz=parameters[5];                                                               //AWS20080711

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy +  Bregz*Bregz);                                   //AWS20080331

  


    Bran=parameters[4]*Breg;                                                           //AWS20080331
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;



 }

 //////////////////////////////////////////////////////////////////////////////////////
// Tinyakov & Tkachev 2002,Astroparticle Physics 18, 165  BSS-A model,
// radial dependence of Bo
// from Elena Orlando

 if(strcmp(name,"tinyakov_exp_ran")==0)            //EO20081118
 {
   double B0;                                                                             

   double  ro=8.5;
   theta=atan2(y,x);
   double p=parameters[2]*dtr;                      //=-8 pitch angle
   double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   
   double z0=parameters[6];                          
   double fz=exp(-fabs(z/z0));                                 
   if(z<0.) fz=-fz;                              // to avoid problem at z=0

   double  psi=1./(tan(p));//radial dependence of opening angle   

   double r=pow(x*x+y*y,0.5);
   double rbreak = parameters[8];                                  //EO20081119

    if(r<=rbreak) B0  =parameters[0]*(ro/rbreak);//parameters[0]=1.4
    if(r> rbreak) B0  =parameters[0]*(ro/r);

    Bregx=B0  *cos(theta-psi*log(r/r0))*sin(p-theta)*fz;                          
    Bregy=B0  *cos(theta-psi*log(r/r0))*cos(p-theta)*fz;
    Bregz=parameters[5];                                                              

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                                  

    // random component can be added for testing
    Bran=parameters[4];
    Branx=Bran*fz;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

 }


////////////////////////////////////////////////////////////////////////////////////////////
// halo dipole component from prouza 2003 formulae, same parameters

  if(strcmp(name,"halo_dipole")==0)
 {

        
    theta=atan2(y,x);             // NB theta defined differently from above
   double r=pow(x*x+y*y,0.5);
   double beta=90.*dtr-atan(z/r);

 if(fabs(z)<=0.3 && r<=0.1) //Btot=const =2mG from prouza
{double k=pow(10.,-3.);
Bregx=-(3./2.)*k*sin(2.*beta)*sin(theta);
Bregy=-(3./2.)*k*sin(2.*beta)*cos(theta);
Bregz=-k*(3.*cos(beta)*cos(beta)-1.);
}

 if(fabs(z)>0.3 && r<=0.1) //Btot=const =2mG from prouza 
{double k=2.*pow(10.,-7.);
Bregx=-(3./2.)*k*sin(2.*beta)*sin(theta)/pow(r,3.);
Bregy=-(3./2.)*k*sin(2.*beta)*cos(theta)/pow(r,3.);
Bregz=-k*(3.*cos(beta)*cos(beta)-1.)/pow(r,3.);
}

if( r>0.1 && r<=2.) // from prouza
{double k=2.*pow(10.,-7.);
Bregx=-(3./2.)*k*sin(2.*beta)*sin(theta)/pow(r,3.);
Bregy=-(3./2.)*k*sin(2.*beta)*cos(theta)/pow(r,3.);
Bregz=-k*(3.*cos(beta)*cos(beta)-1.)/pow(r,3.);
}

if( r>2. && r<=5.) // from prouza
{double k=5.*pow(10.,-7.);
Bregx=-(3./2.)*k*sin(2.*beta)*sin(theta);
Bregy=-(3./2.)*k*sin(2.*beta)*cos(theta);
Bregz=-k*(3.*cos(beta)*cos(beta)-1.);
}

if(  r>5.) // from prouza
{double k=pow(10.,-4.);
Bregx=-(3./2.)*k*sin(2.*beta)*sin(theta)/pow(r,3.);
Bregy=-(3./2.)*k*sin(2.*beta)*cos(theta)/pow(r,3.);
Bregz=-k*(3.*cos(beta)*cos(beta)-1.)/pow(r,3.);
}



    Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                              

// random component set to zero
   Bran =0.0;
   Branx=0.0;  
   Brany=0.0;
   Branz=0.0;


}


/////////////////////////////////////////////////////// EO20090701
  //Sun model for B in the plane. ASS-RING model

 if(strcmp(name,"sun")==0)            //EO20090207
 {
   double B0;                                                                             

   double  ro=8.5;
   theta=atan2(y,x);
   double p=parameters[2]*dtr;                      //=-12 pitch angle
   //double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   double r0=10.;   

   double z0=parameters[6];         //=1 kpc                 
   double fz=exp(-fabs(z/z0));                                 
   if(z<0.) fz=-fz;                              // to avoid problem at z=0

   // double  psi=1./(tan(p));//radial dependence of opening angle   

   double r=pow(x*x+y*y,0.5);
   double rbreak = parameters[8];//rc=5 kpc                                  //EO20081119
   double D1;
   double D2;
   B0=parameters[0];//2 microGauss

   if(r<=rbreak && fabs(z)<=1) D1  = B0;
   if(r<=rbreak && fabs(z)>1) D1  = 0.;//artificial cut not present in Sun formulation 
   if(r> rbreak) D1  =B0*exp(-(r-ro)/r0)*fz;

   if(r>7.5) D2=1;
 if(r<=7.5 & r>6.) D2=-1;
 if(r>5.&& r<=6) D2=1;
 if(r<=5) D2=-1;


    Bregx=D1*D2*sin(p-theta);                          
    Bregy=D1*D2*cos(p-theta);  
    Bregz=parameters[5];//=0 in Sun                                                              

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                                  

    // random component can be added for testing
    Bran=parameters[4];//3 microgauss
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

 }



  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 // from Elena Orlando                                                                    EO20090716
  if(strcmp(name,"random")==0)
  {
 


    double Bo,rscale,zscale;
    double ro=8.5;
    double b_field;

    

    double r=sqrt(x*x+y*y);


   //just to test the random compoennt since the contrubution of the models for the regular is negligible EO20090716
   Bo   = parameters[0]; // Gauss
   rscale=parameters[3]; // kpc
   zscale=parameters[2]; // kpc

   double chi0=parameters[1]*dtr;// not defined in Han; parameter[1]=8 in WMAP Miville
    double chi=chi0*tanh(z/parameters[9]);      //  z dependence, not defined in Han;  //EO20081120

    b_field=Bo *exp(-(r-ro)/rscale)*cos(chi);// * exp(-fabs(z)/zscale);

 

   // the regular field is zero
   Breg =0.0;
   Bregx=0.0;
   Bregy=0.0;
   Bregz=0.0;

   Bran=b_field;
   Branx=Bran;  // put all field in x-direction (has no significance)
   Brany=0.0;
   Branz=0.0;

   theta=0.0;
   phi  =0.0;

  } //only random component

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 // from Elena Orlando                                                                    EO20090716
  if(strcmp(name,"random1")==0)
  {
 


    double Bo,rscale,zscale;
    double ro=8.5;
    double b_field;

    

    double r=sqrt(x*x+y*y);


   //just to test the random compoennt since the contrubution of the models for the regular is negligible EO20090716
   Bo   = parameters[0]; // Gauss
   rscale=parameters[3]; // kpc
   zscale=parameters[2]; // kpc

   double chi0=parameters[1]*dtr;// not defined in Han; parameter[1]=8 in WMAP Miville
    double chi=chi0*tanh(z/parameters[9]);      //  z dependence, not defined in Han;  //EO20081120
    double rc=parameters[4];

    if(r<=rc) b_field=Bo *exp(-(r-ro)/rscale)*cos(chi);// * exp(-fabs(z)/zscale);
    if(r>rc) b_field=Bo *exp(-(rc-ro)/rscale)*cos(chi);// * exp(-fabs(z)/zscale);
 

   // the regular field is zero
   Breg =0.0;
   Bregx=0.0;
   Bregy=0.0;
   Bregz=0.0;

   Bran=b_field;
   Branx=Bran;  // put all field in x-direction (has no significance)
   Brany=0.0;
   Branz=0.0;

   theta=0.0;
   phi  =0.0;

  } //only random component


/////////////////////////////////////////////////////// EO20090701
  //Sun model for B in the plane. ASS-RING model

 if(strcmp(name,"jansson_disk")==0)            //EO20090207
 {
   double B0;                                                                             

   double  ro=8.5;
   theta=atan2(y,x);
   double p=-5*dtr;                      //=-12 pitch angle
   //double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   double r0=5.1;   

   double z0=parameters[6];         //=1 kpc                 
   double fz=exp(-fabs(z/z0));                                 
   if(z<0.) fz=-fz;                              // to avoid problem at z=0

   // double  psi=1./(tan(p));//radial dependence of opening angle   

   double r=pow(x*x+y*y,0.5);
   double rbreak = 5.7;//rc=5 kpc                                  //EO20081119
   double D1;
   double D2;
   double Bc=0.16e-6;
   B0=1.1e-6;//2 microGauss

   if(r<=rbreak) D1  = Bc*fz;
   //  if(r<=rbreak) D1  = 0.;//artificial cut not present in Sun formulation 
   if(r> rbreak) D1  =B0*exp(-(r-ro)/r0)*fz;

   if(r>7.5) D2=1;
 if(r<=7.5 & r>6.) D2=-1;
 if(r>5.&& r<=6) D2=1;
 if(r<=5) D2=-1;


    Bregx=D1*D2*sin(p-theta);                          
    Bregy=D1*D2*cos(p-theta);  
    Bregz=0.0;//=0 in Sun                                                              

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                                  

    // random component can be added for testing
    Bran=0.0;//3 microgauss
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

 }


/////////////////////////////////////////////////////// EO20090701
  //Sun model for B in the plane. ASS-RING model
 //to be finished!!!!
 if(strcmp(name,"jansson_halo")==0)            //EO20090207
 {
   double B0;                                                                             

   double  ro=8.5;
   theta=atan2(y,x);
   double p;                      //=-12 pitch angle
   //double r0=(ro-0.5)*exp((-pi/2)*tan(p));
   double r0=28.;   

   double z0=parameters[6];         //=1 kpc                 
   double fz=exp(-fabs(z/z0));                                 
   if(z<0.) fz=-fz;                              // to avoid problem at z=0

   // double  psi=1./(tan(p));//radial dependence of opening angle   

   double r=pow(x*x+y*y,0.5);
   double rbreak = 8.72;//rc=5 kpc                                  //EO20081119
   // double D1;
   // double D2;
   // double Bc=0.16;
   B0=2.3e-6;//2 microGauss

   if(r<=rbreak) p=-30.*dtr;//     D1  = Bc*fz;
   //  if(r<=rbreak) D1  = 0.;//artificial cut not present in Sun formulation 
   if(r> rbreak) p=-2.*dtr;  // D1  =B0*exp(-r/r0)*fz;

   //   if(r>7.5) D2=1;
   // if(r<=7.5 & r>6.) D2=-1;
   // if(r>5.&& r<=6) D2=1;
   // if(r<=5) D2=-1;


    Bregx=B0*exp(-r/r0)*fz*sin(p-theta);                          
    Bregy=B0*exp(-r/r0)*fz*cos(p-theta);  
    Bregz=0.0;//=0 in Sun                                                              

    Breg=sqrt(Bregx*Bregx + Bregy*Bregy + Bregz*Bregz);                                  

    // random component can be added for testing
    Bran=0.0;//3 microgauss
    Branx=Bran;  // put all field in x-direction (has no significance)
    Brany=0.0;
    Branz=0.0;

 }


 /////////////////////////////////////////////////////////////////////////////////////////////
 // from Elena Orlando                                                                    EO20081114
 //  Sun et al. 2008, A&A 477, 573; Prouza & Smida 2003 A&A 410, 1
 //  halo field, added to regular field

 // always added but controlled by parameters[7] which can be set to zero if halo not required
 // so this code should come at end

 double r=pow(x*x+y*y,0.5);

 double fi=atan2(r,z);
double B0t=parameters[7]; //2 microgauss
 double z1=0.;
 if (fabs(z) <1.5) z1=0.2;
 if (fabs(z)>=1.5) z1=0.4;

 double B1=1/(1+pow((fabs(z)-1.5)/z1, 2));
 double B2=exp(-(r-4)/4);
 double Bt=B0t*B1*B2*r/4;

 theta=atan2(y,x);
 double Bhalox = Bt*cos(theta);
 double Bhaloy =-Bt*sin(theta);
 double Bhaloz=0.;

 if (z<0.){Bhalox=-Bhalox; Bhaloy=-Bhaloy;}

 Bregx += Bhalox;
 Bregy += Bhaloy;
 Bregz += Bhaloz;

 Breg=sqrt(Bregx*Bregx + Bregy*Bregy +  Bregz*Bregz); 

//////////////////////////////////////////////////////////////// 

 


//============================================================================





  if(debug==1)
 cout<<"B model= "<<name
     <<" (x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  "<<"(theta,phi)= ("<<theta<<","<<phi<<")"
     <<"  Breg="<< Breg <<" Bregx="<<  Bregx<<" Bregy="<< Bregy<<" Bregz="<< Bregz
                         <<" Bran="<< Bran <<" Branx="<<  Branx<<" Brany="<< Brany<<" Branz="<< Branz<<endl;


}

////////////////////////////////////////////////////////////////////////////////////////////////////////
double Bperp(double x,double y,double z, double Bx,double By,double Bz, double x0,double y0,double z0)
{
  //returns the field (at x,y,z)  component perpendicular to the line-of sight to x0,y0,z0

  double Bperp_,Bperp2;

  // component parallel to line of sight by projection onto line of sight
  double Bpar=        ((x-x0)*Bx +     (y-y0)*By +    (z-z0)*Bz    ) 
                /sqrt ((x-x0)*(x-x0) + (y-y0)*(y-y0)+ (z-z0)*(z-z0)) ;

  // component perpendicular by Pythagoras: sqrt(B_total^2 - Bpar^2)
  double B_total2=(Bx*Bx + By*By + Bz*Bz);

  Bperp2 =      Bx*Bx + By*By + Bz*Bz - Bpar*Bpar ;
  Bperp_= 0.0;
  if(Bperp2 > 0.0)          // to avoid rounding error problem
  Bperp_ = sqrt(Bperp2);

  //  cout<<"Bperp routine: B_total^2="<<B_total2<< " Bpar="<<Bpar <<" Bperp2="<<Bperp2<<" Bperp_="<<Bperp_<<endl;

  

  return Bperp_;
}


//////////////////////////////////////////////////////////////////////////////////////////////

double B_field_3D_model_tot
(char* name, double *parameters,
 double x, double y, double z,
 int debug=0 )
{
  // return just total B, used e.g. for energy losses

  int options=0;
  double Breg,   Bregx,  Bregy,  Bregz;
  double Bran,   Branx,  Brany,  Branz;
  double B_tot;


  //  debug=1; // test AWS20080331

 B_field_3D_model
( name, parameters,
  x,  y,  z,
  options,
  Breg,  Bregx, Bregy, Bregz,
  Bran,  Branx, Brany, Branz,
  debug );

 B_tot=sqrt(Breg*Breg + Bran*Bran);

 return B_tot;

}

//////////////////////////////////////////////////////////////////////////////////////////////

double B_field_3D_model_tot
(char* name, double *parameters,
 double R, double z,
 int debug=0 )
{
  // return just total B, used e.g. for energy losses

  int options=0;
  double x,y;
  double Breg,   Bregx,  Bregy,  Bregz;
  double Bran,   Branx,  Brany,  Branz;
  double B_tot;

  // compute field on GC-Sun axis since x,y undefined (better would be to take azimuthal average)
  x = R;
  y = 0.0;

 B_field_3D_model
( name, parameters,
  x,  y,  z,
  options,
  Breg,  Bregx, Bregy, Bregz,
  Bran,  Branx, Brany, Branz,
  debug );

 B_tot=sqrt(Breg*Breg + Bran*Bran);

 return B_tot;

}

///////////////////////////////////////////////////////////////
//             test program
///////////////////////////////////////////////////////////////
/*

#include "synchrotron_emissivity.h"
int main()
{

  char name[20];
  double parameters[20];
  double x,y,z;
  int options;
  double  Breg,          Bregx,         Bregy,        Bregz;
  double  Bran,          Branx,         Brany,        Branz;
  int debug;

  strcpy(name,"test");
  strcpy(name,"circular");

  x=8;
  y=1;
  z=1;

  options=0;
  debug=1;

  for (x= -10.; x<+10.; x+=2) // avoid x0,y0,z0
  for (y= -10.; y<+10.; y+=2)
  for (z=  -1; z< +1.; z+=1.) 
  {

 B_field_3D_model
( name, parameters,
  x,  y,  z,
  options,
 Breg,  Bregx, Bregy, Bregz,
  Bran,  Branx, Brany, Branz,
  debug );

 cout<<"(x, y, z) = ("<<x<<", "<<y<<", "<<z<<")  ";
 cout<<" B model= "<<name<<"  Breg="<<  Breg <<" Bregx="<<  Bregx<<" Bregy="<< Bregy<<" Bregz="<< Bregz
                         <<" Bran="<<  Bran <<" Branx="<<  Branx<<" Brany="<< Brany<<" Branz="<< Branz<<endl;

 double x0,y0,z0;
 x0=8.5; // solar position
 y0=0;
 z0=0;

 cout<<"Bperp for regular field ="<<Bperp( x, y, z,  Bregx, Bregy, Bregz,  x0,y0,z0)<<endl;
 
 cout<<"Bperp for random  field ="<<Bperp( x, y, z,  Branx, Brany, Branz,  x0,y0,z0)<<endl;

 // test synchrotron routine

 double gamma,nu;
 double synch_emissivity_total;
 double synch_emissivity_reg;
 double synch_emissivity_par;
 double synch_emissivity_perp;
 double synch_emissivity_random;

 int debug=1;

 gamma=1e4/.511; // 10 GeV electrons = 1e4 MeV
 nu=1.e9; // 1000 MHz


 double Brand; // NB names inconsistent
 Brand=Bran;
 double Bperp_reg=Bperp( x, y, z,  Bregx, Bregy, Bregz,  x0,y0,z0);
 double Bperp_ran=Bperp( x, y, z,  Branx, Brany, Branz,  x0,y0,z0);


 Bperp_reg += 1.e-9; // smaller field gives error in synchrotron.c
 Bperp_ran += 1.e-9; // smaller field gives error in synchrotron.c

 cout<<endl;
 cout<<"testing synchrotron routine"<<endl;
 cout<<endl;

 cout<<"regular field:"<<endl;
 cout<<"gamma= "<<gamma<<" nu="<<nu<<" Bperp_reg="<<Bperp_reg   << " Brand="<<Brand   <<endl;

 synch_emissivity_total
     = synchrotron_emissivity( gamma, nu, Bperp_reg, Brand,
                         synch_emissivity_reg, synch_emissivity_par, synch_emissivity_perp, synch_emissivity_random, debug );


 cout<<"regular field:"<<endl;
 cout<<"gamma= "<<gamma<<" nu="<<nu<<" Bperp_reg="<<Bperp_reg   << " Brand="<<Brand   <<endl

  <<" synch emissivity random   =  "  << synch_emissivity_random <<endl
  <<" synch emissivity reg      =  "  << synch_emissivity_reg    <<endl
  <<" synch emissivity parallel =  "  << synch_emissivity_par    <<endl
  <<" synch emissivity perp     =  "  << synch_emissivity_perp   <<endl 
  <<" synch emissivity total    =  "  << synch_emissivity_total  <<endl       
  <<endl;


 cout<<"random  field:"<<endl;
 cout<<"gamma= "<<gamma<<" nu="<<nu<<" Bperp_ran="<<Bperp_ran   << " Brand="<<Brand   <<endl;

 synch_emissivity_total
     = synchrotron_emissivity( gamma, nu, Bperp_ran, Brand,
                         synch_emissivity_reg, synch_emissivity_par, synch_emissivity_perp, synch_emissivity_random, debug );


 cout<<endl;
 cout<<"random  field:"<<endl;
 cout<<"gamma= "<<gamma<<" nu="<<nu<<" Bperp_ran="<<Bperp_ran   << " Brand="<<Brand   <<endl

  <<" synch emissivity random   = "  << synch_emissivity_random <<endl
  <<" synch emissivity reg      = "  << synch_emissivity_reg    <<endl
  <<" synch emissivity parallel = "  << synch_emissivity_par    <<endl
  <<" synch emissivity perp     = "  << synch_emissivity_perp   <<endl 
  <<" synch emissivity total    = "  << synch_emissivity_total  <<endl       
  <<endl;


  }// for x y z

 return 0;

}
*/
