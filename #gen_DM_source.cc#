
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_DM_source.cc *                             galprop package * 9/09/2005 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|


//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// The routine gen_DM_source calculates the source functions of the products of the
// dark matter (DM) particle annihilation [cm^-3 s^-1 MeV^-1].
// The routine can be used to calculate source function of positrons, electrons,
// and antiprotons.
// Use gen_DM_emiss to define gamma-ray emissivity (cm^-3 s^-1 MeV^-1)
// in terms (dn/dEdt *c/4pi), where n is the number density, c is speed of light.
// The user must use the parameters DM_double0-9 and DM_int0-9 (galdef-file) to 
// specify the Galactic DM profile, branching, decay channels, and spectra (see 
// the template below). The DM profile is defined in the DM_profile routine.
// The profile is then averaged over the grid step (dR,dz) or (dx,dy,dz) with 
// a smaller step: normally 1/10 of the grid size.          IMOS20050912
//
// See example in Moskalenko I.V., Strong A.W. 1999, Phys. Rev. D 60, 063003
// and realization below.
//=="====!===="====!===="====!===="====!===="====!===="====!===="====!===="====!
using namespace std;
#include"galprop_classes.h"
#include"galprop_internal.h"

#include <fort_interface.h>
#include <fstream>
#include <string.h>
#include <iostream>

//extern "C" void RHO_DARKSUSY_F77(double*,double*,double*,double*); //IMOS20060901

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

int Galprop::gen_DM_source(Particle &particle)
{
   cout<<"gen_DM_source"<<endl;
   cout<<"generating "<<particle.name<<" source function for n_spatial_dimensions="
       <<gcr[0].n_spatial_dimensions<<endl;

   double DMwidth,DMbranching,   // annihilation product distribution
     DMsecondary_spectrum,       // spectrum of secondaries from DM annihilation
     DME0,                       // delta function energy used for Green's function
     DMmass  =galdef.DM_double2, // DM particle mass
     DMcs_v  =galdef.DM_double9, // DM <cross_sec*V> -thermally overaged, cm3 s-1 
     dzz=0.01;                   // kpc, gas averaging step
   int stat=0;

 // define the spectra of annihilation products: positrons, electrons, antiprotons

   if(strcmp(particle.name,"DM_positrons")==0)
     {
       DMwidth     =galdef.DM_double3;
       DMbranching =galdef.DM_double4;
     }

   if(strcmp(particle.name,"DM_electrons")==0)
     {
       DMwidth     =galdef.DM_double5;
       DMbranching =galdef.DM_double6;
     }

   if(strcmp(particle.name,"DM_antiprotons")==0)
     {
       DMwidth     =galdef.DM_double7;
       DMbranching =galdef.DM_double8;
     }

// assign the source function (2D)

   if(galaxy.n_spatial_dimensions==2)
     {
       for(int ir=0; ir<gcr[0].n_rgrid; ir++)
	 {
	   for(int iz=0; iz<gcr[0].n_zgrid; iz++)
	     {
	       for(int ip=0; ip<particle.n_pgrid; ip++)
		 {
// test of electron propagation vs analytical calculations IMOS20061030
// to run test, assign galdef.DM_int0=99, other parameters:
// galdef.DM_double6 - the half thickness of the disk source distribution (e.g. 0.1 kpc), the source 
//                     distribution is uniform within the disk; normalization =1 at the normalization energy
// galdef.DM_double7 - the photon field energy density (e.g. 1 eV/cc)
// galdef.DM_double8 - the normalization energy of the electron spectrum (e.g. 10^3 MeV)
// galdef.DM_double9 - the injection spectral index of electrons (e.g. 2.4)
		   if(abs(galdef.DM_int0)==99 && particle.A==0) 
		     {
		       if(strcmp(particle.name,"DM_electrons")==0) //numerical  solution "DM_electrons"
			 particle.secondary_source_function.d2[ir][iz].s[ip]= 
			   (galdef.DM_double6 <= fabs(galaxy.z[iz])) ? 0.:
			   C/4./Pi*pow(particle.Ekin[ip]/galdef.DM_double8,-galdef.DM_double9);
		       if(strcmp(particle.name,"DM_positrons")==0) //analytical solution "DM_positrons"
			 particle.secondary_source_function.d2[ir][iz].s[ip]=0.;
		       continue;		     
		     }
// end of the test area
		   if(galdef.DM_int1==9) // Green's function to work with DarkSUSY IMOS20060901
		     {
		       if(DME0<particle.Ekin[ip] || DME0/DMwidth>particle.Ekin[ip]) continue;
		       particle.secondary_source_function.d2[ir][iz].s[ip]
			 +=pow(DM_profile_av(galaxy.r[ir], galaxy.z[iz], galaxy.dr, galaxy.dz, dzz),2)
			 *DMsecondary_spectrum*DMbranching/4./Pi*C;
		       continue;
		     }
		   if(particle.Etot[ip]*1.e-3<=DMmass) 
		     particle.secondary_source_function.d2[ir][iz].s[ip]+= DMcs_v*
		       pow(DM_profile_av(galaxy.r[ir], galaxy.z[iz], galaxy.dr, galaxy.dz, dzz)/DMmass,2)
		       *C/4./Pi*DMbranching*exp(-pow((DMmass-particle.Etot[ip]*1.e-3)/DMwidth,2))/DMmass*1.e-3;
		 } // ip
	     }  //  iz
	 }  //  ir
     }  //  particle.n_spatial_dimensions==2
   
// assign the source function (3D)

   if(galaxy.n_spatial_dimensions==3)
     {
       for(int ix=0; ix<gcr[0].n_xgrid; ix++)
	 {
	   for(int iy=0; iy<gcr[0].n_ygrid; iy++)
	     {
	       for(int iz=0; iz<gcr[0].n_zgrid; iz++)
		 {
		   for(int ip=0; ip<particle.n_pgrid; ip++)
		     {
		       if(galdef.DM_int1==9) // Green's function to work with DarkSUSY IMOS20060901
			 {
			   if(DME0<particle.Ekin[ip] || DME0/DMwidth>particle.Ekin[ip]) continue;
			   particle.secondary_source_function.d3[ix][iy][iz].s[ip]
			     +=pow(DM_profile_av(galaxy.r[ix], galaxy.r[iy], galaxy.z[iz], galaxy.dx, galaxy.dy, galaxy.dz, dzz),2)
			     *DMsecondary_spectrum*DMbranching/4./Pi*C;
			   continue;
			 }
		       if(particle.Etot[ip]*1.e-3<=DMmass) 
			 particle.secondary_source_function.d3[ix][iy][iz].s[ip]+= DMcs_v*
			   pow(DM_profile_av(galaxy.x[ix], galaxy.y[ix], galaxy.z[iz], galaxy.dx, galaxy.dy, galaxy.dz, dzz)/DMmass,2)
		       *C/4./Pi*DMbranching*exp(-pow((DMmass-particle.Etot[ip]*1.e-3)/DMwidth,2))/DMmass*1.e-3;
		     } //ip
		 }  //  iz
	     }  //  iy
	 }  //  ix
     }  //  particle.n_spatial_dimensions==3
 
 // test printout

   if(galdef.verbose>=2)
     {
       cout<<"   particle.secondary_source_function for "<<particle.name<<endl;
       particle.secondary_source_function.print();
     }
   cout<<" <<<< gen_DM_source"<<endl;
   return stat;
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

int Galprop::gen_DM_emiss()
{
   cout<<"gen_DM_emiss"<<endl;
   double 
     DMmass      =galdef.DM_double2, // DM particle mass
     DMcs_v      =galdef.DM_double9, // DM <cross_sec*V> -thermally overaged, cm3 s-1 
     DMbranching =0.1,
     dzz=0.01;                       // kpc, gas averaging step
   int stat=0;

   galaxy.DM_emiss=0.;

// define the spectra of annihilation products: gammas
   
   if(galdef.n_spatial_dimensions==2)
     {
       cout<<"generating DM emissivity for n_spatial_dimensions="<<galdef.n_spatial_dimensions<<endl;
       for(int ir=0; ir<gcr[0].n_rgrid; ir++)
	 {
	   for(int iz=0; iz<gcr[0].n_zgrid; iz++)
	     {
               for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
		 {
		   if(galaxy.E_gamma[iEgamma]*1.e-3>DMmass) 
		     {
		       galaxy.DM_emiss.d2[ir][iz].s[iEgamma]=0;
		       continue;
		     }
		   galaxy.DM_emiss.d2[ir][iz].s[iEgamma]= DMcs_v *DMbranching/(4.*Pi)// sr^-1 IMOS20060420
		     *pow(DM_profile_av(galaxy.r[ir], galaxy.z[iz], galaxy.dr, galaxy.dz, dzz)/DMmass,2)
		     /galaxy.E_gamma[iEgamma];
		 }
	     }
	 }
     }
   if(galdef.n_spatial_dimensions==3)
     {
       cout<<"generating DM emissivity for n_spatial_dimensions="<<galdef.n_spatial_dimensions<<endl;
       for(int ix=0; ix<gcr[0].n_rgrid; ix++)
	 {
	   for(int iy=0; iy<gcr[0].n_rgrid; iy++)
	     {
	       for(int iz=0; iz<gcr[0].n_zgrid; iz++)
		 {
		   for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
		     {
		       if(galaxy.E_gamma[iEgamma]*1.e-3>DMmass) 
			 {
			   galaxy.DM_emiss.d3[ix][iy][iz].s[iEgamma]=0;
			   continue;
			 }
		       galaxy.DM_emiss.d3[ix][iy][iz].s[iEgamma]=  DMcs_v *DMbranching/(4.*Pi) // sr^-1 IMOS20060420
			 *pow(DM_profile_av(galaxy.x[ix], galaxy.y[ix], galaxy.z[iz], galaxy.dx, galaxy.dy, galaxy.dz, dzz)/DMmass,2)
			 /galaxy.E_gamma[iEgamma];
		     }
		 }
	     }
	 }
     }
   cout<<" <<<< gen_DM_emiss"<<endl;
   return(stat);
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

double Galprop::DM_profile(double Xkpc, double Ykpc, double Zkpc)
{
  double R=sqrt(Xkpc*Xkpc+Ykpc*Ykpc+Zkpc*Zkpc),
    Rsun =8.5,                     //kpc, galactocentric distance of the solar system 
    Rc         =galdef.DM_double0, //core radius
    rho0       =galdef.DM_double1; //local DM mass density
  int profile_key =galdef.DM_int0; //profile type
  
  switch(profile_key)
    {
    case 0:   //NFW profile
      return(rho0*Rc/R*pow(1.+R/Rc,-2));
      
    case 1:   //isothermal profile
      return(rho0*(pow(Rc,2)+pow(Rsun,2))/(pow(Rc,2)+pow(R,2)));
      
    case 2:   //Evans profile
      return(rho0*pow(pow(Rc,2)+pow(Rsun,2),2)/(3.*pow(Rc,2)+pow(Rsun,2))
	     *(3.*pow(Rc,2)+pow(R,2))/pow(pow(Rc,2)+pow(R,2),2));
      
    case 3:   //alternative profile
      return(rho0*pow(Rc+Rsun,2)/pow(Rc+R,2));
      
    case 9:   //DarkSUSY profile (use only if the DarkSUSY and GALPROP combined) IMOS20060901
      RHO_DARKSUSY_F77(&Xkpc,&Ykpc,&Zkpc,&rho0);

      if(rho0<0.)
	{
	  cout<<"gen_DM_source: rho_darksusy() function is not defined"<<endl;
	  exit(0);
	}
      return(rho0);

    default:
      return(rho0);
    }
}
  

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

 double Galprop::DM_profile_av(double r,double z,double dr,double dz,double dzz)
   {  
     double DM_profile_av_=0.0;
     int nuse=0;
     
     for (double zz=z-dz/2.; zz<=z+dz/2.; zz+=dzz)
       for (double rr=r-dr/2.; rr<=r+dr/2.; rr+=dr/10.)
	 { 
	   if (rr<0.) continue;
	   DM_profile_av_+=DM_profile(rr,0,zz);
	   nuse++; 
	 }
     return (DM_profile_av_/nuse);
   }
 
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
 
double Galprop::DM_profile_av(double x,double y,double z,double dx,double dy,double dz,double dzz)
   {  
     double DM_profile_av_=0.0;
     int nuse=0;
     
     for (double zz=z-dz/2.; zz<=z+dz/2.; zz+=dzz)
       for (double xx=x-dx/2.; xx<=x+dx/2.; xx+=dx/10.)
	 for (double yy=y-dy/2.; yy<=y+dy/2.; yy+=dy/10.)
	   {
	     DM_profile_av_+=DM_profile(xx,yy,zz);
	     nuse++;
	   }
     return DM_profile_av_/nuse;
   }
 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
int Galprop::gen_DM_annihilation(Particle &particle)
{
   cout<<"gen_DM_annihilation"<<endl;
   cout<<"generating "<<particle.name<<" source function for n_spatial_dimensions="
       <<gcr[0].n_spatial_dimensions<<endl;

   int stat=0;

   double ss = 1000;
   double temp;
   double mod_DM_mass;
   for(int ip=0; ip<particle.n_pgrid; ip++){
     temp = particle.p[ip]-galdef.DM_mass;
     if(temp<ss){
       ss = abs(particle.p[ip]-galdef.DM_mass);
       mod_DM_mass = particle.p[ip];
     }
   }
   cerr << " DM mass[MeV] : " << mod_DM_mass << endl;
     
 // define the spectra of annihilation products: positrons, electrons, antiprotons

   if(strcmp(particle.name,"primary_DM_positron")==0){
     cout<<"assign the primary positron source function of DM"<<endl;
     if(galaxy.n_spatial_dimensions==2){
       for(int ir=0; ir<gcr[0].n_rgrid; ir++){
	 for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	   for(int ip=0; ip<particle.n_pgrid; ip++){
	     particle.secondary_source_function.d2[ir][iz].s[ip] = 
	       galdef.DM_boostfactor*galdef.DM_sigmav*pow(DM_Distribution(particle.r[ir],particle.z[iz]),2)/(2*pow(mod_DM_mass,2))*func_DM_Annihilation(particle.p[ip],mod_DM_mass)*C/4/Pi;
	     if((particle.r[ir]==-20&&particle.z[iz]==-4)&&abs(particle.p[ip]-mod_DM_mass)<0.01) cerr << particle.secondary_source_function.d2[ir][iz].s[ip] << endl; 
	   } // ip
	 }  //  iz
       }  //  ir
     }  //  particle.n_spatial_dimensions==2
   }
   
   if(strcmp(particle.name,"primary_DM_electron")==0){
     cout<<"assign the primary DM electron source function of DM"<<endl;
     if(galaxy.n_spatial_dimensions==2){
       for(int ir=0; ir<gcr[0].n_rgrid; ir++){
	 for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	   for(int ip=0; ip<particle.n_pgrid; ip++){
	     particle.secondary_source_function.d2[ir][iz].s[ip] = 
	       galdef.DM_boostfactor*galdef.DM_sigmav*pow(DM_Distribution(particle.r[ir],particle.z[iz]),2)/(2*pow(mod_DM_mass,2))*func_DM_Annihilation(particle.p[ip],mod_DM_mass)*C/4/Pi;
	   } // ip
	 }  //  iz
       }  //  ir
     }  //  particle.n_spatial_dimensions==2
   }
   
   /*
   if(strcmp(particle.name,"primary_DM_antiproton")==0){
     ifstream fin;
     fin.open("antiproton_DM.dat",ios::in);
     cout<<"assign the primary antiprotons source function of DM"<<endl;
     if(galaxy.n_spatial_dimensions==2){
       for(int ip=0; ip<particle.n_pgrid; ip++){
	 for(int ir=0; ir<gcr[0].n_rgrid; ir++){
	   for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	     //if(particle.Etot[ip]*1.e-3<=DMmass) 
	     fin>>source;
	     particle.secondary_source_function.d2[ir][iz].s[ip] = source*C/4/Pi;
	   } // ip
	 }  //  iz
       }  //  ir
     }  //  particle.n_spatial_dimensions==2
   }
   */

   cout<<" <<<< gen_DM_annihilation"<<endl;
   return stat;
}
//Source spectrum in DM annihilation model
double Galprop::func_DM_Annihilation(double E_,double mod_DM_mass_){
  int model_num = galdef.DM_model;
  double mod_DM_mass = mod_DM_mass_;
  double E = E_;
  double source;
  switch(model_num){
  case 1: //mono
    if(abs(E-mod_DM_mass)<0.001) source = galdef.power_law_1;
    else source = 0;
    break;
  case 2: //flat
    if(E<=mod_DM_mass) source = 1/(mod_DM_mass-0);
    else source = 0;
    break;
  case 3: //double-peak
    if(E<=mod_DM_mass) source = 3/pow(mod_DM_mass,3)*(pow(E-mod_DM_mass/2,2)+pow(mod_DM_mass,2)/4);
    else source = 0;
    break;
  }
  
  return source;
}
double Galprop::DM_Distribution(double r_, double z_){
  double rho;
  double r = r_;
  double z = z_;
  const double r_sun   = 8.5;  //kpc
  const double r_c     = 3.5;  //kpc
  const double rho_sun = 0.3;  //GeV/cm^3
  rho = rho_sun*(pow(r_c,2)+pow(r_sun,2))/(pow(r_c,2)+pow(r,2)+pow(z,2));
  
  return rho;
}
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
int Galprop::gen_galactic_Pulsar(Particle &particle)
{
   cout<<"gen_galactic_Pulsar"<<endl;
   cout<<"generating "<<particle.name<<" source function for n_spatial_dimensions="
       <<gcr[0].n_spatial_dimensions<<endl;

   int stat=0;
   double E_cut = galdef.energy_cutoff;
   double index = galdef.power_law_1;

 // define the spectra of annihilation products: positrons, electrons, antiprotons

   if(strcmp(particle.name,"primary_DM_electron")==0){
     cout<<"assign the primary electron source function of Pulsar"<<endl;
     if(galaxy.n_spatial_dimensions==2){
       for(int ir=0; ir<gcr[0].n_rgrid; ir++){
	 for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	   for(int ip=0; ip<particle.n_pgrid; ip++){
	     particle.secondary_source_function.d2[ir][iz].s[ip] = 
	       Pulsar_Distribution(particle.r[ir],particle.z[iz])*func_galactic_Pulsar(particle.p[ip],E_cut,index);
	   } // ip
	 }  //  iz
       }  //  ir
     }  //  particle.n_spatial_dimensions==2
   }
   
   if(strcmp(particle.name,"primary_DM_positron")==0){
     cout<<"assign the primary DM positron source function of DM"<<endl;
     if(galaxy.n_spatial_dimensions==2){
       for(int ir=0; ir<gcr[0].n_rgrid; ir++){
	 for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	   for(int ip=0; ip<particle.n_pgrid; ip++){
	     particle.secondary_source_function.d2[ir][iz].s[ip] = 
	       Pulsar_Distribution(particle.r[ir],particle.z[iz])*func_galactic_Pulsar(particle.p[ip],E_cut,index);
	   } // ip
	 }  //  iz
       }  //  ir
     }  //  particle.n_spatial_dimensions==2
   }
   
   cout<<" <<<< gen_galactic_Pulsar"<<endl;
   return stat;
}
//Source spectrum in galactic Pulsars
double Galprop::func_galactic_Pulsar(double E,double E_cut,double a){
  double source;
  source = pow(E,-a)*exp((-E/E_cut));
  return source;
}
double Galprop::Pulsar_Distribution(double r, double z){
  double rho;
  double a = 2.35;
  double b = 5.56;  
  const double r_sun   = 8.5;  //kpc
  const double z_s     = 0.2;  //kpc 
  rho = pow(r/r_sun,a)*exp(-b*(r-r_sun)/r_sun)*exp(-abs(z)/z_s);  

  return rho;
}
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
int Galprop::gen_TD_injection(Particle &particle) //YO/20140220 -> 20141010
{
   int stat=0;

 // define the spectra of annihilation products: positrons, electrons
   int model_num = galdef.TD_SNR_model;
   switch(model_num){
     static int flag = 0;
   case 1: //SNR burst-like model
     if(flag==0) cout << "  <<< SNR burst-like model" << endl;
     flag = 1;
     SNR_burstlike(particle);
     break;
   case 2: //time-dependent pulsar decay model
     if(flag==0) cout << "  <<< pulsar time-dependent decay" << endl; 
     flag = 1;
     pulsar_TD_decay(particle);
     break;
   case 3: //exponential SNR decay model
     if(flag==0) cout << "  <<< SNR exponential decay model" << endl; 
     flag = 1;
     SNR_exp_decay(particle);
     break;
   case 4: //energy dependent TeV electron model
     if(flag==0) cout << "  <<< energy dependent TeV electron model" << endl; 
     flag = 1;
     TeVelectron_energy_dependent_decay(particle);
     break;
   }
   
   //cout <<" <<<< gen_TD_injection"<<endl;
   return stat;
}
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
int Galprop::SNR_burstlike(Particle &particle) //YO/20141010 -> 20141117
{
  //burst-like SNR model : T.Kobayashi, Apj:601,340-351, 2004.
  double source;
  double E_cut = galdef.energy_cutoff;
  double index;
  
  //Positron injection
  if(strcmp(particle.name,"primary_TD_positron")==0){
    
    if(galdef.check_num==325){
      cerr << " propagating time = " << particle.m_step_total_time << " dt = " << particle.m_dt << endl;
    }
    
    if((particle.m_final_time-particle.m_step_total_time)+particle.m_dt>particle.m_source_age/*galdef.source_age*/ &&
       particle.m_source_age/*galdef.source_age*/>=(particle.m_final_time-particle.m_step_total_time)){
      
      if(galdef.check_num==325){
	cerr << " SNR model = " << galdef.TD_SNR_model << endl;
	cerr << " total time - propagating time = " << particle.m_final_time-particle.m_step_total_time 
	     << " source_age = " << particle.m_source_age/*galdef.source_age*/ << endl;
	cerr << " => USE CR source" << endl;
      }
      
      //ifstream fin;
      //fin.open("source.dat",ios::in);
      if(galaxy.n_spatial_dimensions==3){
	for(int ix=0; ix<gcr[0].n_xgrid; ix++){
	  for(int iy=0; iy<gcr[0].n_ygrid; iy++){
	    for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	      for(int ip=0; ip<particle.n_pgrid; ip++){
		particle.secondary_source_function.d3[ix][iy][iz].s[ip] = 0;
		//fin>>source;
		//particle.secondary_source_function.d3[ix][iy][iz].s[ip] = source*C/4/Pi;
		/*if(particle.x[ix]==galdef.SNR_x && particle.y[iy]==galdef.SNR_y && particle.z[iz]==galdef.SNR_z)
		  particle.secondary_source_function.d3[ix][iy][iz].s[ip] = function(particle.p[ip],E_cut,index);
		else
		particle.secondary_source_function.d3[ix][iy][iz].s[ip] = 0;*/
	      }}}}
      }  // particle.n_spatial_dimensions==3
    }
    else{
      
      if(galdef.check_num==325){
	cerr << " SNR model = " << galdef.TD_SNR_model << endl;
	cerr << " total time - propagating time = " << particle.m_final_time-particle.m_step_total_time 
	     << "  source_age = " << particle.m_source_age/*galdef.source_age*/ << endl;
	cerr << " => NOT use CR source" << endl;
      }
      
      if(galaxy.n_spatial_dimensions==3){
	for(int ix=0; ix<gcr[0].n_xgrid; ix++){
	  for(int iy=0; iy<gcr[0].n_ygrid; iy++){
	    for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	      for(int ip=0; ip<particle.n_pgrid; ip++){
		particle.secondary_source_function.d3[ix][iy][iz].s[ip] = 0;
	      }}}}
      }  // particle.n_spatial_dimensions==3
    }
  } // end positron
  
  //Electron injection
  if(strcmp(particle.name,"primary_TD_electron")==0){
    
    if(galdef.check_num==325){
      cerr << " propagating time = " << particle.m_step_total_time << " dt = " << particle.m_dt << endl;
    }
    
    if((particle.m_final_time-particle.m_step_total_time)+particle.m_dt>particle.m_source_age/*galdef.source_age*/ &&
       particle.m_source_age/*galdef.source_age*/>=(particle.m_final_time-particle.m_step_total_time)){
      
      if(galdef.check_num==325){
	cerr << " SNR model = " << galdef.TD_SNR_model << endl;
	cerr << " total time - propagating time= " << particle.m_final_time-particle.m_step_total_time 
	     << " source_age = " << particle.m_source_age/*galdef.source_age*/ << endl;
	cerr << " => USE CR source" << endl;
      }
      
      //ifstream fin;
      //fin.open("electron_source.dat",ios::in);
      if(galaxy.n_spatial_dimensions==3){
	for(int ix=0; ix<gcr[0].n_xgrid; ix++){
	  for(int iy=0; iy<gcr[0].n_ygrid; iy++){
	    for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	      for(int ip=0; ip<particle.n_pgrid; ip++){
		//fin>>source;
		//particle.secondary_source_function.d3[ix][iy][iz].s[ip] = source*C/4/Pi;
		if(abs(particle.x[ix]-particle.m_sourceX/*galdef.SNR_x*/)<0.01 &&
		   abs(particle.y[iy]-particle.m_sourceY/*galdef.SNR_y*/)<0.01 &&
		   abs(particle.z[iz]-particle.m_sourceZ/*galdef.SNR_z*/)<0.01){
		  if(particle.p[ip]<galdef.power_law_br) index = galdef.power_law_1;
		  if(particle.p[ip]>=galdef.power_law_br)  index = galdef.power_law_2;
		  particle.secondary_source_function.d3[ix][iy][iz].s[ip] = particle.m_nor_fac*func_SNR_burstlike(particle.p[ip],E_cut,index)*C/4/Pi;
		  if(galdef.check_num==1117){
		    cerr << " source coordinate (" << particle.x[ix] << "," << particle.y[iy] << "," << particle.z[iz] << ")";
		    cerr << ", normalization factor = " << particle.m_nor_fac;
		    cerr << " : func_SNR_burstlike = " << func_SNR_burstlike(particle.p[ip],E_cut,index) << endl;
		  }
		}
		else particle.secondary_source_function.d3[ix][iy][iz].s[ip] = 0;
	      }}}}
      }  // particle.n_spatial_dimensions==3
    }
    else{
      
      if(galdef.check_num==325){
	cerr << " SNR model = " << galdef.TD_SNR_model << endl;
	cerr << " total time - propagating time = " << particle.m_final_time-particle.m_step_total_time 
	     << " source_age = " << particle.m_source_age/*galdef.source_age*/ << endl;
	cerr << " => NOT use CR source" << endl;
      }
      
      if(galaxy.n_spatial_dimensions==3){
	for(int ix=0; ix<gcr[0].n_xgrid; ix++){
	  for(int iy=0; iy<gcr[0].n_ygrid; iy++){
	    for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	      for(int ip=0; ip<particle.n_pgrid; ip++){
		particle.secondary_source_function.d3[ix][iy][iz].s[ip] = 0;
	      }}}}
      }  // particle.n_spatial_dimensions==3
    }
  } // end electron

  return 0;
  
}
//get normalization factor with Simpson method
int Galprop::nor_SNR_burstlike(Particle &particle){
  cerr << " >>>> calculating normalization factor at SNR in burst-like model";
  
  double E_cut = galdef.energy_cutoff;
  double index = galdef.power_law_1-1;
  double E_total = galdef.total_energy; //erg

  double E_min = galdef.Ekin_min;
  double E_max = galdef.Ekin_max*1; //or10^20eV?

  double div_num = (E_max-E_min)*10;
  double dx = (E_max-E_min)/(2*div_num);
  
  double integral,x,nor_fac;
  double dy;

  //必要な係数　index,E_cut,E_total => galdefから読み込み
  x = E_min;
  if(x>=galdef.power_law_br) index = galdef.power_law_2-1;
  integral = func_SNR_burstlike(x,E_cut,index);

  double l = 0;
  for(double i=1;i<div_num;i=i+1.0){
    if(x<galdef.power_law_br) index = galdef.power_law_1-1;
    if(x>=galdef.power_law_br)  index = galdef.power_law_2-1;
    dy = 4*func_SNR_burstlike((x+dx),E_cut,index)+2*func_SNR_burstlike((x+2*dx),E_cut,index);
    integral += dy;
    x += 2*dx;
    l=l+1.0;
  }

  integral += (4*func_SNR_burstlike((x+dx),E_cut,index)+func_SNR_burstlike((x+2*dx),E_cut,index));
  integral *= dx/3;

  nor_fac = E_total/MEV2ERG/integral/(galdef.start_timestep*year2sec)/(particle.dx*particle.dy*particle.dz*pow(kpc2cm,3));
  particle.m_nor_fac = nor_fac;
  
  cerr << " : i= " << l << ", integral = " << integral << ", nor_fac = " << particle.m_nor_fac << endl;
  
  return 0;
  
}
//Source spectrum in SNR burst-like model
double Galprop::func_SNR_burstlike(double E,double E_cut,double a){
  double result;
  result = pow(E,-a)*exp((-E/E_cut));
  //cerr << "E= " << E << " Ecut = " << E_cut << " index= " << a << endl;
  //cerr << "result = " << result << endl;
  return result;
}
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
int Galprop::pulsar_TD_decay(Particle &particle) //YO/20141010
{
  //time-dependent pulsar decay model : N.Kawanaka, arXiv:0903.3782, 2009.
  double source;
  double tau_0 = 1.e5;
  
  //Positron injection
  if(strcmp(particle.name,"primary_TD_positron")==0){
     
    if(galdef.check_num==325){
      cerr << "propagating time = " << particle.m_step_total_time << " dt = " << particle.m_dt << endl;
    }
    
    if((particle.m_final_time-particle.m_step_total_time)<=galdef.source_age){
      
      if(galdef.check_num==325){
	cerr << " SNR model = " << galdef.TD_SNR_model << endl;
	cerr << " total time - propagating time = " << particle.m_final_time-particle.m_step_total_time 
	     << " source_age = " << galdef.source_age << endl;
	cerr << " => USE CR source" << endl;
      }
      
      ifstream fin;
      fin.open("source.dat",ios::in);
      if(galaxy.n_spatial_dimensions==3){
	for(int ix=0; ix<gcr[0].n_xgrid; ix++){
	  for(int iy=0; iy<gcr[0].n_ygrid; iy++){
	    for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	      for(int ip=0; ip<particle.n_pgrid; ip++){
		fin>>source;
		particle.secondary_source_function.d3[ix][iy][iz].s[ip] =
		  source*C/4/Pi/pow((1+(galdef.source_age-(particle.m_final_time-particle.m_step_total_time)))/tau_0,2);
	      }}}}
      }  // particle.n_spatial_dimensions==3
    }
    else{
      
      if(galdef.check_num==325){
	cerr << " SNR model = " << galdef.TD_SNR_model << endl;
	cerr << " total time - propagating time = " << particle.m_final_time-particle.m_step_total_time 
	     << "  source_age = " << galdef.source_age << endl;
	cerr << " => NOT use CR source" << endl;
      }
      
      if(galaxy.n_spatial_dimensions==3){
	for(int ix=0; ix<gcr[0].n_xgrid; ix++){
	  for(int iy=0; iy<gcr[0].n_ygrid; iy++){
	    for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	      for(int ip=0; ip<particle.n_pgrid; ip++){
		particle.secondary_source_function.d3[ix][iy][iz].s[ip] = 0;
	      }}}}
      }  // particle.n_spatial_dimensions==3
    }
  } // end positron
  
  //Electron injection
  if(strcmp(particle.name,"primary_TD_electron")==0){
    
    if(galdef.check_num==325){
      cerr << "propagating time = " << particle.m_step_total_time << " dt = " << particle.m_dt << endl;
    }

    if((particle.m_final_time-particle.m_step_total_time)<=galdef.source_age){
      
      if(galdef.check_num==325){
	cerr << " SNR model = " << galdef.TD_SNR_model << endl;
	cerr << " total time - propagating time= " << particle.m_final_time-particle.m_step_total_time 
	     << " source_age = " << galdef.source_age << endl;
	cerr << " => USE CR source" << endl;
      }
      
      ifstream fin;
      fin.open("electron_source.dat",ios::in);
      cout<<" assign the primary electron source function of time dependent injection"<<endl;
      if(galaxy.n_spatial_dimensions==3){
	for(int ix=0; ix<gcr[0].n_xgrid; ix++){
	  for(int iy=0; iy<gcr[0].n_ygrid; iy++){
	    for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	      for(int ip=0; ip<particle.n_pgrid; ip++){
		fin>>source;
		particle.secondary_source_function.d3[ix][iy][iz].s[ip] =
		  source*C/4/Pi/pow((1+(galdef.source_age-(particle.m_final_time-particle.m_step_total_time)))/tau_0,2);
	      }}}}
      }  // particle.n_spatial_dimensions==3
    }
    else{
      
      if(galdef.check_num==325){
	cerr << " SNR model = " << galdef.TD_SNR_model << endl;
	cerr << " total time - propagating time = " << particle.m_final_time-particle.m_step_total_time 
	     << " source_age = " << galdef.source_age << endl;
	cerr << " => NOT use CR source" << endl;
      }
      
      if(galaxy.n_spatial_dimensions==3){
	for(int ix=0; ix<gcr[0].n_xgrid; ix++){
	  for(int iy=0; iy<gcr[0].n_ygrid; iy++){
	    for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	      for(int ip=0; ip<particle.n_pgrid; ip++){
		particle.secondary_source_function.d3[ix][iy][iz].s[ip] = 0;
	      }}}}
      }  // particle.n_spatial_dimensions==3
    }
  } // end electron
  
  return 0;

}
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
int Galprop::SNR_exp_decay(Particle &particle) //YO/20141010 -> 20141202
{
  //exponential pulsar decay model : N.Kawanaka, arXiv:0903.3782, 2009.
  double source;
  double index = galdef.power_law_1;
  
  //Positron injection
  if(strcmp(particle.name,"primary_TD_positron")==0){
     
    if(galdef.check_num==325){
      cerr << "propagating time = " << particle.m_step_total_time << " dt = " << particle.m_dt << endl;
    }
    
    if((particle.m_final_time-particle.m_step_total_time)<=galdef.source_age){
      
      if(galdef.check_num==325){
	cerr << " SNR model = " << galdef.TD_SNR_model << endl;
	cerr << " total time - propagating time = " << particle.m_final_time-particle.m_step_total_time 
	     << " source_age = " << galdef.source_age << endl;
	cerr << " => USE CR source" << endl;
      }
      
      //ifstream fin;
      //fin.open("source.dat",ios::in);
      if(galaxy.n_spatial_dimensions==3){
	for(int ix=0; ix<gcr[0].n_xgrid; ix++){
	  for(int iy=0; iy<gcr[0].n_ygrid; iy++){
	    for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	      for(int ip=0; ip<particle.n_pgrid; ip++){
		if(abs(particle.x[ix]-particle.m_sourceX)<0.01 &&
		   abs(particle.y[iy]-particle.m_sourceY)<0.01 &&
		   abs(particle.z[iz]-particle.m_sourceZ)<0.01){
		  particle.secondary_source_function.d3[ix][iy][iz].s[ip] = particle.m_nor_fac*func_SNR_exp_decay1(particle.p[ip],index)*func_SNR_exp_decay2((galdef.source_age-(particle.m_final_time-particle.m_step_total_time)))*C/4/Pi;
		  if(galdef.check_num==1117){
		    cerr << " source coordinate (" << particle.x[ix] << "," << particle.y[iy] << "," << particle.z[iz] << ")";
		    cerr << ", normalization factor = " << particle.m_nor_fac;
		    cerr << " : secondary_source_function = " << particle.secondary_source_function.d3[ix][iy][iz].s[ip] << endl;
		  }
		}
		else particle.secondary_source_function.d3[ix][iy][iz].s[ip] = 0;
	      }}}}
      }  // particle.n_spatial_dimensions==3
    }
    else{
      
      if(galdef.check_num==325){
	cerr << " SNR model = " << galdef.TD_SNR_model << endl;
	cerr << " total time - propagating time = " << particle.m_final_time-particle.m_step_total_time 
	     << "  source_age = " << galdef.source_age << endl;
	cerr << " => NOT use CR source" << endl;
      }
      
      if(galaxy.n_spatial_dimensions==3){
	for(int ix=0; ix<gcr[0].n_xgrid; ix++){
	  for(int iy=0; iy<gcr[0].n_ygrid; iy++){
	    for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	      for(int ip=0; ip<particle.n_pgrid; ip++){
		particle.secondary_source_function.d3[ix][iy][iz].s[ip] = 0;
	      }}}}
      }  // particle.n_spatial_dimensions==3
    }
  } // end positron
  
  //Electron injection
  if(strcmp(particle.name,"primary_TD_electron")==0){
    //particle.t_d += particle.m_dt;
    
    if(galdef.check_num==325){
      cerr << "propagating time = " << particle.m_step_total_time << " dt = " << particle.m_dt << endl;
    }

    if((particle.m_final_time-particle.m_step_total_time)<=galdef.source_age){
      
      if(galdef.check_num==325){
	cerr << " SNR model = " << galdef.TD_SNR_model << endl;
	cerr << " total time - propagating time= " << particle.m_final_time-particle.m_step_total_time 
	     << " source_age = " << galdef.source_age << endl;
	cerr << " => USE CR source" << endl;
      }
      
      //ifstream fin;
      //fin.open("electron_source.dat",ios::in);
      if(galaxy.n_spatial_dimensions==3){
	for(int ix=0; ix<gcr[0].n_xgrid; ix++){
	  for(int iy=0; iy<gcr[0].n_ygrid; iy++){
	    for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	      for(int ip=0; ip<particle.n_pgrid; ip++){
		//fin>>source;
		/*particle.secondary_source_function.d3[ix][iy][iz].s[ip] =
		  source*C/4/Pi*exp(-((galdef.source_age-((particle.m_final_time-particle.m_step_total_time)))*log(4))/tau_0);*/
		if(abs(particle.x[ix]-particle.m_sourceX)<0.01 &&
		   abs(particle.y[iy]-particle.m_sourceY)<0.01 &&
		   abs(particle.z[iz]-particle.m_sourceZ)<0.01){
		  particle.secondary_source_function.d3[ix][iy][iz].s[ip] = particle.m_nor_fac*func_SNR_exp_decay1(particle.p[ip],index)*func_SNR_exp_decay2((galdef.source_age-(particle.m_final_time-particle.m_step_total_time)))*C/4/Pi;
		  if(galdef.check_num==1117){
		    cerr << " source coordinate (" << particle.x[ix] << "," << particle.y[iy] << "," << particle.z[iz] << ")";
		    cerr << ", normalization factor = " << particle.m_nor_fac;
		    cerr << " : secondary_source_function = " << particle.secondary_source_function.d3[ix][iy][iz].s[ip] << endl;
		  }
		}
		else particle.secondary_source_function.d3[ix][iy][iz].s[ip] = 0;
	      }}}}
      }  // particle.n_spatial_dimensions==3
    }
    else{
      
      if(galdef.check_num==325){
	cerr << " SNR model = " << galdef.TD_SNR_model << endl;
	cerr << " total time - propagating time = " << particle.m_final_time-particle.m_step_total_time 
	     << " source_age = " << galdef.source_age << endl;
	cerr << " => NOT use CR source" << endl;
      }
      
      if(galaxy.n_spatial_dimensions==3){
	for(int ix=0; ix<gcr[0].n_xgrid; ix++){
	  for(int iy=0; iy<gcr[0].n_ygrid; iy++){
	    for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	      for(int ip=0; ip<particle.n_pgrid; ip++){
		particle.secondary_source_function.d3[ix][iy][iz].s[ip] = 0;
	      }}}}
      }  // particle.n_spatial_dimensions==3
    }
  } // end electron

  return 0;

}
//Get normalization factor with Simpson method
int Galprop::nor_SNR_exp_decay(Particle &particle){ //YO/20141203 未完
  cerr << " >>>> calculating normalization factor at SNR in exp-decay model";
  
  //normalization factor for func1
  double index = galdef.power_law_1-1;   //times E
  double E_total = galdef.total_energy;
  double E_min = galdef.Ekin_min;
  double E_max = galdef.Ekin_max;

  double div_num1 = (E_max-E_min)*10;
  double dx1 = (E_max-E_min)/(2*div_num1);

  double integral1,x1,nor_fac1,dy1;

  x1 = E_min;
  integral1 = func_SNR_exp_decay1(x1,index);
  
  for(double i=1;i<div_num1;i=i+1.){
    dy1 = 4*func_SNR_exp_decay1((x1+dx1),index)+2*func_SNR_exp_decay1((x1+2*dx1),index);
    integral1 += dy1;
    x1 += 2*dx1;
   }
  
  integral1 += (4*func_SNR_exp_decay1((x1+dx1),index)+func_SNR_exp_decay1((x1+2*dx1),index));
  integral1 *= dx1/3;

  //normalization factor for func2
  double source_age = galdef.source_age;
  double t_max = galdef.source_age;
  double t_min = 0;

  double div_num2 = (t_max-t_min)*10;
  double dx2 = (t_max-t_min)/(2*div_num2);

  double integral2,x2,nor_fac2,dy2;

  x2 = t_min;
  integral2 = func_SNR_exp_decay2(x2);
  
  for(double i=1;i<div_num2;i=i+1.){
    dy2 = 4*func_SNR_exp_decay2((x2+dx2))+2*func_SNR_exp_decay2((x2+2*dx2));
    integral2 += dy2;
    x2 += 2*dx2;
  }
  
  integral2 += (4*func_SNR_exp_decay2((x2+dx2))+func_SNR_exp_decay2((x2+2*dx2)));
  integral2 *= dx2/3;

  //func1 times func2
  double nor_fac;
  nor_fac = E_total/MEV2ERG/integral1/integral2/year2sec/(particle.dx*particle.dy*particle.dz*pow(kpc2cm,3));
  particle.m_nor_fac = nor_fac;

  cerr << " : nor_fac = " << particle.m_nor_fac << endl;

  return 0;

}
//Source spectrum in SNR exp-decay model
double Galprop::func_SNR_exp_decay1(double E,double a){ //YO/20141203
  double result1;
  result1 = pow(E,-a);
  return result1;
} 
double Galprop::func_SNR_exp_decay2(double t){ //YO/20141204
  double result2;
  double tau_0 = 1.e5;
  result2 = exp((-t*log(4)/tau_0));
  //result2 = exp((-(galdef.source_age-(particle.m_final_time-particle.m_step_total_time))*log(4)/tau_0));
  return result2;
} 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
int Galprop::TeVelectron_energy_dependent_decay(Particle &particle) //YO/20141021
{
  //TeV electron energy-dependent injection model : N.Kawanaka, Apj729:93, 2011.
  double source;
  double source1, source2; //escape1, escape2
  source1 = 0;
  source2 = 0;
  
  //Positron injection
  if(strcmp(particle.name,"primary_TD_positron")==0){
     
    if(galdef.check_num==325){
      cerr << "propagating time = " << particle.m_step_total_time << " dt = " << particle.m_dt << endl;
    }
    
    if((particle.m_final_time-particle.m_step_total_time)<=galdef.source_age){
      
      if(galdef.check_num==325){
	cerr << " SNR model = " << galdef.TD_SNR_model << endl;
	cerr << " total time - propagating time = " << particle.m_final_time-particle.m_step_total_time 
	     << " source_age = " << galdef.source_age << endl;
	cerr << " => USE CR source" << endl;
      }
      
      ifstream fin;
      fin.open("source.dat",ios::in);
      if(galaxy.n_spatial_dimensions==3){
	for(int ix=0; ix<gcr[0].n_xgrid; ix++){
	  for(int iy=0; iy<gcr[0].n_ygrid; iy++){
	    for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	      for(int ip=0; ip<particle.n_pgrid; ip++){
		fin>>source;
		particle.secondary_source_function.d3[ix][iy][iz].s[ip] = source*C/4/Pi;
	      }}}}
      }  // particle.n_spatial_dimensions==3
    }
    else{
      
      if(galdef.check_num==325){
	cerr << " SNR model = " << galdef.TD_SNR_model << endl;
	cerr << " total time - propagating time = " << particle.m_final_time-particle.m_step_total_time 
	     << "  source_age = " << galdef.source_age << endl;
	cerr << " => NOT use CR source" << endl;
      }
      
      if(galaxy.n_spatial_dimensions==3){
	for(int ix=0; ix<gcr[0].n_xgrid; ix++){
	  for(int iy=0; iy<gcr[0].n_ygrid; iy++){
	    for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	      for(int ip=0; ip<particle.n_pgrid; ip++){
		particle.secondary_source_function.d3[ix][iy][iz].s[ip] = 0;
	      }}}}
      }  // particle.n_spatial_dimensions==3
    }
  } // end positron
  
  //Electron injection
  if(strcmp(particle.name,"primary_TD_electron")==0){
    double t_sedov = 200; // the begging of Sedov phase [yr]
    double tau_0 = 10^4; // spin-down time scale [yr]
    double a = 2.0; // power law index
    
    if(galdef.check_num==325){
      cerr << "propagating time = " << particle.m_step_total_time << " dt = " << particle.m_dt << endl;
    }

    if((particle.m_final_time-particle.m_step_total_time)<=galdef.source_age){
      double E_esc;
      E_esc = pow(10,9.5)*pow((galdef.source_age-(particle.m_final_time-particle.m_step_total_time))/t_sedov,-2.6);
       // escape energy threshold [MeV]
      
      if(galdef.check_num==325){
	cerr << " SNR model = " << galdef.TD_SNR_model << endl;
	cerr << " total time - propagating time= " << particle.m_final_time-particle.m_step_total_time 
	     << " source_age = " << galdef.source_age << endl;
	cerr << " escape energy threshold [MeV] = " << E_esc << endl;
	cerr << " => USE CR source" << endl;
      }
      
      ifstream fin;
      fin.open("electron_source.dat",ios::in);
      cout<<" assign the primary electron source function of time dependent injection"<<endl;
      if(galaxy.n_spatial_dimensions==3){
	for(int ix=0; ix<gcr[0].n_xgrid; ix++){
	  for(int iy=0; iy<gcr[0].n_ygrid; iy++){
	    for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	      for(int ip=0; ip<particle.n_pgrid; ip++){
		fin>>source; //q_0*E^-a
		//source1 = flux_escape1(particle,ip,E_esc,tau_0,source); //escape term 1
		//cerr << " >>>>>>>>>>>>>>>>>>start escape2>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl; 
		source2 = flux_escape2(particle,ip,E_esc,a,t_sedov,source); //escape term 2
		//cerr << " >>>>>>>>>>>>>>>>>>end escape2>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl; 
		particle.secondary_source_function.d3[ix][iy][iz].s[ip] = (/*source1+*/source2)*C/4/Pi;
	      }}}}
      }  // particle.n_spatial_dimensions==3
    }
    else{
      
      if(galdef.check_num==325){
	cerr << " SNR model = " << galdef.TD_SNR_model << endl;
	cerr << " total time - propagating time = " << particle.m_final_time-particle.m_step_total_time 
	     << " source_age = " << galdef.source_age << endl;
	cerr << " => NOT use CR source" << endl;
      }
      
      if(galaxy.n_spatial_dimensions==3){
	for(int ix=0; ix<gcr[0].n_xgrid; ix++){
	  for(int iy=0; iy<gcr[0].n_ygrid; iy++){
	    for(int iz=0; iz<gcr[0].n_zgrid; iz++){
	      for(int ip=0; ip<particle.n_pgrid; ip++){
		particle.secondary_source_function.d3[ix][iy][iz].s[ip] = 0;
	      }}}}
      }  // particle.n_spatial_dimensions==3
    }
  } // end electron

  return 0;

}
//Calculate electron escape term
double Galprop::flux_escape1(Particle &particle,int ip,double E_esc,double tau_0,double source){ //YO/20141024
  double source1;
  
  if(E_esc <= particle.p[ip]){
    source1 = source*pow(1+(galdef.source_age-(particle.m_final_time-particle.m_step_total_time))/tau_0,-2);
  }
  if(E_esc > particle.p[ip]){
    source1 = 0;
  }
 
  double temp = particle.p[ip]/pow(10,6); //expの中身はあらかじめ計算しないとダメみたい
  source1 *= exp(-temp);
  
  return source1; //E_cut=10TeV
}
double Galprop::flux_escape2(Particle &particle,int ip,double E_esc,double a,double t_sedov,double source){ //YO/20141024
  double source2 = 0;
  double N_conf=0;
  double t,t_cre;

  t=galdef.source_age-(particle.m_final_time-particle.m_step_total_time); //SNから現在までの時間
  t_cre=pow(pow(t_sedov,-2.6)/pow(10,9.5)*particle.p[ip]*pow(t,2/5),-1/2.2); //しきい値を超えるときの時間

  if(t<t_cre) t_cre = t;

  /*if(t>=3000)
    cerr << "            >>>>>>>> ip = " << ip << "   t = " << t << "   t_cre = " << t_cre << endl;*/

  //get Nconf with integration
  if(t_sedov>t_cre) source2 = 0;
  if(t_sedov<=t_cre){
    for(int t_i=t_sedov;t_i<=t_cre;t_i=t_i+1/*particle.m_dt*/){
      double e = -particle.p[ip]*pow((t/t_i),2/5)/pow(10,6);
      double temp = source/pow((1+t/t_i),2)*pow((t/t_i),(2/5*(1-a)))*exp(e); //Ecut=10TeV
      /*if(ip==29){cerr << "        energy = " << particle.p[ip] << " t_i = " << t_i << " e = " << e << " exp(e) = " << exp(e) << endl;}*/
      N_conf += temp; //Nconf(E,t)
    }
  }
  cerr << "      Nconf(E,t) = " << N_conf << endl;
  //temp *= q_0*pow(particle.p[ip],-a); //Nconf(E,t)

  double temp2=pow(10,9.5)*(-2.6)*pow((t/t_sedov),-3.6)/t_sedov + 2*particle.p[ip]/5/t;
  if(particle.p[ip-1]<=E_esc && particle.p[ip]>E_esc)
    source2=-1*N_conf*temp2; //N'esc2(E,t)
  else source2=0;

  return source2;
}
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
int Galprop::mk_primary_electron_table(Particle &particle) //YO/20140616
{
  
  //primary_electronのfluxをテキストに保存する
  
  ofstream fout;
  fout.open("primary_electron.dat");
  
  for(int ix=0; ix<particle.n_xgrid; ix++){
    for(int iy=0; iy<particle.n_ygrid; iy++){
      for(int iz=0; iz<particle.n_zgrid; iz++){
	for(int ip=0; ip<particle.n_pgrid; ip++){
	  //cerr << particle.x[ix] << " " << particle.y[iy] <<" " << particle.z[iz] << " " << particle.primary_source_function.d3[ix][iy][iz].s[ip] << endl;
	  fout << particle.primary_source_function.d3[ix][iy][iz].s[ip] << endl;
	}//ip
      }//iz
    }//iy
  }//ix
	  
  return 0;
  
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
int Galprop::read_primary_electron_table(Particle &particle) //YO/20140616
{
  
  double flux;

  ifstream fin;
  fin.open("primary_electron.dat");

  for(int ix=0; ix<particle.n_xgrid; ix++){
    for(int iy=0; iy<particle.n_ygrid; iy++){
      for(int iz=0; iz<particle.n_zgrid; iz++){
	for(int ip=0; ip<particle.n_pgrid; ip++){
	  fin >> flux;
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] = flux;
	}//ip
      }//iz
    }//iy
  }//ix
  
  return 0;
  
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
int Galprop::mod_primary_electron_table(Particle &particle) //YO/20140616
{

  //保存したprimary_electronのテーブルを読み込み、近傍ソースの寄与だけ0にしたものを
  //primary_source_functionとして使用

  double flux;

  ifstream fin;
  fin.open("primary_electron.dat");

  for(int ix=0; ix<particle.n_xgrid; ix++){
    for(int iy=0; iy<particle.n_ygrid; iy++){
      for(int iz=0; iz<particle.n_zgrid; iz++){
	for(int ip=0; ip<particle.n_pgrid; ip++){
	  fin >> flux;
	  if(sqrt((particle.x[ix]-8.5)*(particle.x[ix]-8.5)+particle.y[iy]*particle.y[iy]+particle.z[iz]*particle.z[iz])<=galdef.near_distance){
	    particle.primary_source_function.d3[ix][iy][iz].s[ip] = 0;
	  }
	  else{
	    particle.primary_source_function.d3[ix][iy][iz].s[ip] = flux;
	  }
	}//ip
      }//iz
    }//iy
  }//ix

  return 0;

}
