
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_secondary_antiproton_source.cc *          galprop package * 2001/05/11
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// The routine to calculate the antiproton source function. 
//
// CR density gcr.cr_density is in c/4pi * n(E) [cm s^-1 sr^-1 cm^-3 MeV^-1]
//
// The routine ANTIPROTON written in FORTRAN-77 is designed to calculate
// the antiproton (+antineutron) production spectrum vs. momentum (barn/GeV). 
// Antiproton momentum and nucleus momentum (GeV) per nucleon are used as input
// parameters as well as beam and target nuclei atomic numbers.
//
// The antiproton source function [cm^-2 s^-2 sr^-1 MeV^-1] as used in galprop is 
// defined as following (c/4pi * q)  [q = cm^-3 s^-1 MeV^-1]:
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
// 
// To transfer to units cm^2/MeV we need a factor= 1.0e-24 *1.0e-3.
// Ref.: Moskalenko I.V. et al. 2002, ApJ 565, 280
//=="====!===="====!===="====!===="====!===="====!===="====!===="====!===="====!

using namespace std;//AWS20050624

#include"galprop_classes.h"
#include"galprop_internal.h"

#include <fort_interface.h>

#include <cstring>

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Galprop::gen_secondary_antiproton_source(Particle &particle)
{
   if(galdef.verbose>=1) cout<<"gen_secondary_antiproton_source"<<endl;
   if(galdef.verbose>=1) cout<<"generating "<<particle.name<<" source function for n_spatial_dimensions="
       <<gcr[0].n_spatial_dimensions<<endl;

   if(strcmp(particle.name,"secondary_antiprotons")!=0)
   {  cout<<"invalid particle "<<particle.name<<endl; return 2; }

   int stat=0, iprotons=-1, iHelium =-1,  Z1, A1, Z2, A2;
   float cs_p_HI, cs_p_He, cs_He_HI, cs_He_He;
   Distribution protons;                 // IMOS20000606.6

// identify CR protons                   // IMOS20000606.7
   if(galdef.n_spatial_dimensions==2) protons.init(gcr[0].n_rgrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
   if(galdef.n_spatial_dimensions==3) protons.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
   protons=0.;
   for(int i=0; i<n_species; i++)  
      if(101==100*gcr[i].Z+gcr[i].A)
      {
         iprotons=i;
 	 protons+=gcr[iprotons].cr_density;
         if(galdef.verbose>=1) cout<<"  CR protons found as species #"<<iprotons<<endl;
      }
   if(iprotons==-1) { cout<<"CR protons not found!"<<endl; return 1; }
 
// identify CR Helium
   for(int i=0; i<n_species; i++) if(204 == 100*gcr[i].Z+gcr[i].A) iHelium =i;
   if(iHelium ==-1) { cout<<"CR Helium  not found!"<<endl; return 1; }
   else if(galdef.verbose>=1) cout<<"  CR Helium  found as species #"<<iHelium <<endl;
//Gulli20070821 
#pragma omp parallel for schedule(dynamic) default(shared) private(Z1,Z2,A1,A2,cs_p_HI,cs_He_HI,cs_p_He,cs_He_He)
   for(int ip_sec=0; ip_sec<particle.n_pgrid; ip_sec++)
   {
      for(int ip=0; ip<gcr[iprotons].n_pgrid; ip++)
      {
         Z1=gcr[iprotons].Z;  A1=gcr[iprotons].A;  Z2=1;  A2=1;    // beam+target: p+HI
         cs_p_HI =antiproton_cc(galdef.total_cross_section,particle.p[ip_sec]*1.e-3, gcr[iprotons].p[ip]*1.e-3, Z1,A1,Z2,A2); // IMOS20010511 IMOS20000601

// secondary_antiprotons =1 uses scaling to calc.; =2 uses factors by Simon et al. 1998
         if(galdef.secondary_antiprotons == 2)                                                     // IMOS20000802.2
         {
            cs_p_HI*=0.12/pow(particle.Ekin[ip_sec]/1000,1.67)+1.78;
            cs_He_HI = cs_p_He = cs_He_He = 0.;
         }
         else
	 {
            Z1=gcr[iprotons].Z;  A1=gcr[iprotons].A;  Z2=2;  A2=4;    // beam+target: p+He
            cs_p_He =antiproton_cc(galdef.total_cross_section,particle.p[ip_sec]*1.e-3, gcr[iprotons].p[ip]*1.e-3, Z1,A1,Z2,A2); // IMOS20010511 IMOS20000601

            Z1=gcr[iHelium ].Z;  A1=gcr[iHelium ].A;  Z2=1;  A2=1;    // beam+target: He+HI
            cs_He_HI=antiproton_cc(galdef.total_cross_section,particle.p[ip_sec]*1.e-3, gcr[iHelium ].p[ip]*1.e-3, Z1,A1,Z2,A2); // IMOS20010511 IMOS20000601

            Z1=gcr[iHelium ].Z;  A1=gcr[iHelium ].A;  Z2=2;  A2=4;    // beam+target: He+He
            cs_He_He=antiproton_cc(galdef.total_cross_section,particle.p[ip_sec]*1.e-3, gcr[iHelium ].p[ip]*1.e-3, Z1,A1,Z2,A2); // IMOS20010511 IMOS20000601
	 }
         if(galaxy.n_spatial_dimensions==2)
         {
            for(int ir=0; ir<gcr[iprotons].n_rgrid; ir++)
            {
               for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
               {
                  particle.secondary_source_function.d2[ir][iz].s[ip_sec ]+=  
                     (galaxy.n_HI.d2[ir][iz].s[0]+2.0*galaxy.n_H2.d2[ir][iz].s[0]+galaxy.n_HII.d2[ir][iz].s[0])
                  *( (cs_p_HI  +cs_p_He *galdef.He_H_ratio) 
                     *protons.d2[ir][iz].s[ip] *gcr[iprotons].Ekin[ip]                  // IMOS20000606.8
                    +(cs_He_HI +cs_He_He*galdef.He_H_ratio) 
                     *gcr[iHelium ].cr_density.d2[ir][iz].s[ip] *gcr[iHelium ].Ekin[ip] *gcr[iHelium].A );
               }  //  iz
            }  //  ir
         }  //  particle.n_spatial_dimensions==2

         if(galaxy.n_spatial_dimensions==3)
         {
            for(int ix=0; ix<gcr[iprotons].n_xgrid; ix++)
            {
               for(int iy=0; iy<gcr[iprotons].n_ygrid; iy++)
               {
                  for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
                  {
                     particle.secondary_source_function.d3[ix][iy][iz].s[ip_sec ]+=
                        (galaxy.n_HI.d3[ix][iy][iz].s[0]+2.0*galaxy.n_H2.d3[ix][iy][iz].s[0]+galaxy.n_HII.d3[ix][iy][iz].s[0])
                     *( (cs_p_HI  +cs_p_He *galdef.He_H_ratio) 
                        *protons.d3[ix][iy][iz].s[ip] *gcr[iprotons].Ekin[ip]           // IMOS20000606.9
                       +(cs_He_HI +cs_He_He*galdef.He_H_ratio) 
                        *gcr[iHelium ].cr_density.d3[ix][iy][iz].s[ip] *gcr[iHelium ].Ekin[ip] *gcr[iHelium].A );        
                  }  //  iz
               }  //  iy
            }  //  ix
         }  //  particle.n_spatial_dimensions==3
      }  //  ip
   }  //  ip_sec

   double factor=1.e-24 *1.e-3 *C *log(galdef.Ekin_factor); // transformation to cm2/MeV and constant factors
   particle.secondary_source_function *= factor;

   protons.delete_array();                  // IMOS20000606.10

   if(galdef.verbose>=2)
   {
      cout<<"   particle.secondary_source_function for "<<particle.name<<endl;
      particle.secondary_source_function.print();
    }
   if(galdef.verbose>=1) cout<<" <<<< gen_secondary_antiproton_source"<<endl;
   return stat;
}
