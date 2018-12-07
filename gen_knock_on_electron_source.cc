
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_knock_on_electron_source.cc *             galprop package * 2/2/2006 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|


//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// The routine to calculate the knock-on electron source function IMOS20060504 
//
// CR density gcr.cr_density is in c/4pi * n(E) [cm s^-1 sr^-1 cm^-3 MeV^-1]
//
// The routine knock_on_cross_section calculates the knock-on electron production 
// spectrum vs. energy (cm^2/MeV). Lorentz factors of the projectile particle
// and electron, and the atomic number of the projectile.  
//
// The knock-on source function [cm^-2 s^-2 sr^-1 MeV^-1] as used in
// galprop is defined as following (c/4pi * q)  [q = cm^-3 s^-1 MeV^-1]:
//                ___      ___  
//         c      \        \    /                   c   d sigma_ij(p,p')
// q(p) * --- = c /__  n_i /__  \ dp' beta n_j(p') ---  --------------- ,
//        4pi    i=H,He     j   /                  4pi        dp  
// 
// where n_i is the gas density, d sigma_ij(p,p')/dp is
// the production cross section, n_j(p') is the CR species density, 
// and p' is the total momentum of a nucleus.
// Substitution of dp' with d(log Ekin) gives:
//                ___                          ___  
//       c        \        /                   \               c  d sigma_ij(p,Ekin)
// q(p)*--- = c A /__ n_i  \ d(log Ekin)  Ekin /__  n_j(Ekin) --- -----------------
//      4pi      i=H,He    /                    j             4pi       dp  
//                         ___     ___      ___
//                         \       \        \              c   d sigma_ij(p,Ekin)
//      = c A /\(log Ekin) /__ n_i /__ Ekin /__ n_j(Ekin) ---  -----------------,
//                        i=H,He   Ekin      j            4pi        dp             
// 
// where /\=Delta, and we used dp'=1/beta A Ekin d(log Ekin).
// Since positrons/electrons are assumed massless Etot=p.
//=="====!===="====!===="====!===="====!===="====!===="====!===="====!===="====!

using namespace std;
#include"galprop_classes.h"
#include"galprop_internal.h"
#include <cstring>

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

int Galprop::gen_knock_on_electron_source(Particle &particle)
{
  cout<<" >>>> gen_knock_on_electron_source"<<endl;
  cout<<" generating "<<particle.name<<" source function for n_spatial_dimensions="
      <<gcr[0].n_spatial_dimensions<<endl;
  
  int key1=-9999;
  if(strcmp(particle.name,"knock_on_electrons")==0) key1=0;
  if(key1==-9999) { cout<<" invalid particle "<<particle.name<<endl;  return 2; }
  
  int stat=0, iprotons=-1, iHelium =-1, i;
  Distribution protons;
  
// identify CR protons
  if(galdef.n_spatial_dimensions==2) protons.init(gcr[0].n_rgrid,                 gcr[0].n_zgrid, gcr[0].n_pgrid);
  if(galdef.n_spatial_dimensions==3) protons.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
  protons=0.;
  for(i=0; i<n_species; i++)  
    if(101==100*gcr[i].Z+gcr[i].A)
      {
	iprotons=i;
	protons+=gcr[iprotons].cr_density;
	cout<<"  CR protons found as species #"<<iprotons<<endl;
      }
  if(iprotons==-1) { cout<<"  CR protons not found!"<<endl;  return 1; }

// identify CR Helium
  for(i=0; i<n_species; i++)  if(204==100*gcr[i].Z+gcr[i].A)  iHelium =i;
  if(iHelium ==-1) cout<<"  CR Helium  not found!"<<endl;
  else cout<<"  CR Helium  found as species #"<<iHelium <<endl;

  if(galdef.knock_on_electrons==2) iHelium =-1;

  int Z[2]={1,2}, A[2]={1,4}; // p, He
  double cs[2][2];

  for(int ip_sec=0; ip_sec<particle.n_pgrid; ip_sec++)
    {
      for(int ip=0; ip<gcr[iprotons].n_pgrid; ip++)
	{
//  beam+target: p+H
	  cs[0][0] =pow(1.*Z[0],2)*Z[0]*knock_on_cross_section(gcr[iprotons].gamma[ip], particle.gamma[ip_sec], A[0]);
//  beam+target: p+He
	  cs[0][1] =pow(1.*Z[0],2)*Z[1]*knock_on_cross_section(gcr[iprotons].gamma[ip], particle.gamma[ip_sec], A[0]);

	  if(iHelium !=-1)
	    {
//  beam+target: He+H
	  cs[1][0] =pow(1.*Z[1],2)*Z[0]*knock_on_cross_section(gcr[iprotons].gamma[ip], particle.gamma[ip_sec], A[1]);
//  beam+target: He+He
	  cs[1][1] =pow(1.*Z[1],2)*Z[1]*knock_on_cross_section(gcr[iprotons].gamma[ip], particle.gamma[ip_sec], A[1]);
	    }
	 
//	  if(ip_sec==5 && ip==gcr[iprotons].n_pgrid/3) 
//	    cout<<" knock_on "<<" "<<ip_sec<<" "<<ip<<" "<<cs[0][0]<<" "<<cs[0][1]<<" "<<cs[1][0]<<" "<<cs[1][1]<<endl;

	  if(galaxy.n_spatial_dimensions==2)
            for(int ir=0; ir<gcr[iprotons].n_rgrid; ir++)
	      for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
		{
                 particle.secondary_source_function.d2[ir][iz].s[ip_sec ]+=
		    (galaxy.n_HI.d2[ir][iz].s[0]+2.0*galaxy.n_H2.d2[ir][iz].s[0])
                    *(cs[0][0]  +cs[0][1] *galdef.He_H_ratio) 
		    *protons.d2[ir][iz].s[ip] *gcr[iprotons].Ekin[ip]*gcr[iprotons].A;

                  if(iHelium !=-1) 
		    particle.secondary_source_function.d2[ir][iz].s[ip_sec ]+=
		      (galaxy.n_HI.d2[ir][iz].s[0]+2.0*galaxy.n_H2.d2[ir][iz].s[0])
		      *(cs[1][0] +cs[1][1] *galdef.He_H_ratio) 
		      *gcr[iHelium ].cr_density.d2[ir][iz].s[ip] *gcr[iHelium ].Ekin[ip]*gcr[iHelium].A;
		}  //  iz //  ir //  particle.n_spatial_dimensions==2
	  
	  if(galaxy.n_spatial_dimensions==3)
            for(int ix=0; ix<gcr[iprotons].n_xgrid; ix++)
	      for(int iy=0; iy<gcr[iprotons].n_ygrid; iy++)
		for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
		  {
		    particle.secondary_source_function.d3[ix][iy][iz].s[ip_sec ]+=
		      (galaxy.n_HI.d3[ix][iy][iz].s[0]+2.0*galaxy.n_H2.d3[ix][iy][iz].s[0])
		      *(cs[0][0]  +cs[0][1] *galdef.He_H_ratio) 
		      *protons.d3[ix][iy][iz].s[ip] *gcr[iprotons].Ekin[ip]*gcr[iprotons].A;
		    if(iHelium !=-1) 
		      particle.secondary_source_function.d3[ix][iy][iz].s[ip_sec ]+=
			(galaxy.n_HI.d3[ix][iy][iz].s[0]+2.0*galaxy.n_H2.d3[ix][iy][iz].s[0])
			*(cs[1][0] +cs[1][1]*galdef.He_H_ratio) 
			*gcr[iHelium ].cr_density.d3[ix][iy][iz].s[ip] *gcr[iHelium ].Ekin[ip]*gcr[iHelium].A;
                  }  //  iz //  iy  //  ix  //  particle.n_spatial_dimensions==3
	}  //  ip
    }  //  iEgamma

   double factor= C*log(galdef.Ekin_factor); // constant factor
   if(galdef.knock_on_electrons==2) factor*=1.75; //to account for all nuclei Z>1
   particle.secondary_source_function *= factor;
   protons.delete_array();

   if(galdef.verbose>=2)
   {
      cout<<"   particle.secondary_source_function for "<<particle.name<<endl;
      particle.secondary_source_function.print();
   }
   cout<<" <<<< gen_knock_on_electron_source"<<endl;
   return stat;
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// A routine to calculate knock-on electron cross section
// Refs: Abraham etal. 1966, PR 150, 1088; Berrington & Dermer 2003, ApJ 594,709
// Input: 
//   gam_p -Lorentz-factor of the projectile particle
//   gam_e -Lorentz-factor of the knock-on electron
//   A     -atomic number of the projectile nucleus 
// Output: cross section in cm^2/MeV
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
double Galprop::knock_on_cross_section(double gam_p, double gam_e, int Ap)
{
  double gam_max, gam_min, S, Ma=Ap*Mp*1.e3;  // Ma -projectile mass in MeV 

  gam_max =1. +2.*(pow(gam_p,2)-1.)*Ma*Ma/(2.*Mele*Ma*gam_p +Mele*Mele +Ma*Ma);
  if (gam_e > gam_max) return(0.);  // max Lorentz-factor of electron

  gam_min =Mele/Ma*(gam_e-1.)/2. 
    +sqrt((gam_e+1.)/2.+(pow(gam_e,2)-1.)*pow(Mele/Ma/2.,2));
  if (gam_p < gam_min) return(0.);  // min Lorentz-factor of projectile

  S =2.*pow(gam_e-1.,-2) 
    -(2.*Mele*Ma*gam_p +Mele*Mele +Ma*Ma)/pow(gam_p*Ma,2)/(gam_e-1.)
    +pow(Mele/(Ma*gam_p),2);

  return(Pi*pow(Rele,2)/Mele *S/(1.-pow(gam_p,-2)));
}
