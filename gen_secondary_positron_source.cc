
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_secondary_positron_source.cc *            galprop package * 6/02/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|


//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// The routine to calculate the secondary positron and electron source function. 
//
// CR density gcr.cr_density is in c/4pi * n(E) [cm s^-1 sr^-1 cm^-3 MeV^-1]
//
// The routine PP_MESON written in FORTRAN-77 is designed to calculate
// sec. positron (or sec. electron) production spectrum vs. energy (barn/GeV). 
// Positron/electron energy and total nucleus momentum (GeV) are used as input
// parameters as well as beam and target nuclei atomic numbers.
//
// The positron/electron source function [cm^-2 s^-2 sr^-1 MeV^-1] as used in
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
// 
// To transfer to units cm^2/MeV we need a factor= 1.0e-24 *1.0e-3.
// Ref.: Moskalenko I.V., Strong A.W. 1998, ApJ 493, 694
//=="====!===="====!===="====!===="====!===="====!===="====!===="====!===="====!
using namespace std;//AWS20050624

#include"galprop_classes.h"
#include"galprop_internal.h"

#include <fort_interface.h>

#include <string>
#include <cstring>
#include <sstream>

#include <ErrorLogger.h>

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

int Galprop::gen_secondary_positron_source(Particle& particle) {

  INFO("Entry");

  //cout<<"gen_secondary_positron_source"<<endl;
  //cout<<"generating "<<particle.name<<" source function for n_spatial_dimensions="
  //  <<gcr[0].n_spatial_dimensions<<endl;
  
  ostringstream buf;
  buf << "Generating " << particle.name << " source function for spatial dimensions " << gcr[0].n_spatial_dimensions;
  INFO(buf.str());

  //int key1=-9999;
  //   if(strcmp(particle.name,"secondary_positrons")==0) key1=+3;  // e+ from pi+ decay including kaons
  //   if(strcmp(particle.name,"secondary_electrons")==0) key1=-3;  // e- from pi- decay including kaons

  const string particleName = particle.name;

  const int key1 =  (particleName == "secondary_positrons" ?  4 :                     //TAP20091010 
		    (particleName == "secondary_electrons" ? -4 : -9999));            //TAP20091010

  // The following assigned the positrons incorrectly and was corrected TAP on 20091010 and documented by AWS on 20091012
  //if(strcmp(particle.name,"secondary_electrons")==0 || strcmp(particle.name,"secondary_positrons")==0) key1=4; //IMOS20060724 error fix (e+/- from muon decay)

  if (-9999 == key1) {

    ostringstream buf;
    buf << "Invalid particle " << particleName;
    INFO(buf.str());
    INFO("Exit");
    return 2;

  }

  //if(key1==-9999) { cout<<"invalid particle "<<particle.name<<endl;  return 2; }
  
  int stat=0, iprotons=-1, iHelium =-1,  NA1, NA2, i;
  double cs_p_HI, cs_p_He, cs_He_HI, cs_He_He;
  Distribution protons;                 // IMOS20000606.1
  
  // identify CR protons                   // IMOS20000606.2
  if (2 == galdef.n_spatial_dimensions) 
    protons.init(gcr[0].n_rgrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
  
  if (3 == galdef.n_spatial_dimensions) 
    protons.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);

  protons=0.;

  for(i=0; i<n_species; i++)  
    if(101==100*gcr[i].Z+gcr[i].A)
      {
	iprotons=i;
	protons+=gcr[iprotons].cr_density;
	buf.str("");
	buf<<"  CR protons found as species #"<<iprotons;
	INFO(buf.str());
      }
  if(iprotons==-1) { cout<<"  CR protons not found!"<<endl;  return 1; }
  
  // identify CR Helium
  for(i=0; i<n_species; i++)  if(204==100*gcr[i].Z+gcr[i].A)  iHelium =i;
  if(iHelium ==-1) WARNING("  CR Helium  not found!"); // IMOS20030217
  else {
     buf.str("");
     buf<<"  CR Helium  found as species #"<<iHelium;
     INFO(buf.str());
  }
  
  //Gulli20070810
#pragma omp parallel for schedule(dynamic) default(shared) private(NA1,NA2,cs_p_HI,cs_p_He,cs_He_HI,cs_He_He)
  for(int ip_sec=0; ip_sec<particle.n_pgrid; ip_sec++) {

    //    cout << ip_sec << endl;

    for(int ip=0; ip<gcr[iprotons].n_pgrid; ip++) {

      //ostringstream buf;
      //buf << "Generating for secondary energy " << ip_sec << " " << particle.Etot[ip_sec] << " primary energy " << ip << " " << gcr[iprotons].p[ip];
	//INFO(buf.str());

      NA1=1;   NA2=1;  //  beam+target: p+HI
      cs_p_HI =pp_meson_cc(particle.Etot[ip_sec]*1.e-3, gcr[iprotons].p[ip]*1.e-3, NA1, NA2, key1 );
      NA1=1;   NA2=4;  //  beam+target: p+He
      cs_p_He =pp_meson_cc(particle.Etot[ip_sec]*1.e-3, gcr[iprotons].p[ip]*1.e-3, NA1, NA2, key1 );
      if(iHelium !=-1)                                                                      // IMOS20030217
	{ 
	  NA1=4;   NA2=1;  //  beam+target: He+HI
	  cs_He_HI=pp_meson_cc(particle.Etot[ip_sec]*1.e-3, gcr[iHelium ].p[ip]*1.e-3, NA1, NA2, key1 );
	  NA1=4;   NA2=4;  //  beam+target: He+He
	  cs_He_He=pp_meson_cc(particle.Etot[ip_sec]*1.e-3, gcr[iHelium ].p[ip]*1.e-3, NA1, NA2, key1 );
	}
      
      if(galaxy.n_spatial_dimensions==2)
	for(int ir=0; ir<gcr[iprotons].n_rgrid; ir++)
	  for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
	    {
	      particle.secondary_source_function.d2[ir][iz].s[ip_sec ]+=                   // IMOS20030217
		(galaxy.n_HI.d2[ir][iz].s[0]+2.0*galaxy.n_H2.d2[ir][iz].s[0]+galaxy.n_HII.d2[ir][iz].s[0])
		*(cs_p_HI  +cs_p_He *galdef.He_H_ratio) 
                     *protons.d2[ir][iz].s[ip] *gcr[iprotons].Ekin[ip];
                  if(iHelium !=-1) particle.secondary_source_function.d2[ir][iz].s[ip_sec ]+=  // IMOS20030217
                     (galaxy.n_HI.d2[ir][iz].s[0]+2.0*galaxy.n_H2.d2[ir][iz].s[0]+galaxy.n_HII.d2[ir][iz].s[0])
                    *(cs_He_HI +cs_He_He*galdef.He_H_ratio) 
                     *gcr[iHelium ].cr_density.d2[ir][iz].s[ip] *gcr[iHelium ].Ekin[ip]*gcr[iHelium].A;
               }  //  iz //  ir //  particle.n_spatial_dimensions==2
 
         if(galaxy.n_spatial_dimensions==3)
            for(int ix=0; ix<gcr[iprotons].n_xgrid; ix++)
               for(int iy=0; iy<gcr[iprotons].n_ygrid; iy++)
                  for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
		  {
                     particle.secondary_source_function.d3[ix][iy][iz].s[ip_sec ]+=                  // IMOS20030217
                        (galaxy.n_HI.d3[ix][iy][iz].s[0]+2.0*galaxy.n_H2.d3[ix][iy][iz].s[0]+galaxy.n_HII.d3[ix][iy][iz].s[0])
                       *(cs_p_HI  +cs_p_He *galdef.He_H_ratio) 
                        *protons.d3[ix][iy][iz].s[ip] *gcr[iprotons].Ekin[ip];
                     if(iHelium !=-1) particle.secondary_source_function.d3[ix][iy][iz].s[ip_sec ]+= // IMOS20030217
                        (galaxy.n_HI.d3[ix][iy][iz].s[0]+2.0*galaxy.n_H2.d3[ix][iy][iz].s[0]+galaxy.n_HII.d3[ix][iy][iz].s[0])
                       *(cs_He_HI +cs_He_He*galdef.He_H_ratio) 
                        *gcr[iHelium ].cr_density.d3[ix][iy][iz].s[ip] *gcr[iHelium ].Ekin[ip]*gcr[iHelium].A;
                  }  //  iz //  iy  //  ix  //  particle.n_spatial_dimensions==3
      }  //  ip
   }  //  iEgamma
 
   double factor=1.e-24 *1.e-3 *C *log(galdef.Ekin_factor); // transformation to cm2/MeV and constant factors
   particle.secondary_source_function *= factor;
   protons.delete_array();                         // IMOS20000606.5

   if(galdef.verbose>=2)
   {
      cout<<"   particle.secondary_source_function for "<<particle.name<<endl;
      particle.secondary_source_function.print();
   }
   //cout<<" <<<< gen_secondary_positron_source"<<endl;
   
   INFO("Exit");
   return stat;
}
