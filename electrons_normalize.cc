
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * electrons_normalize.cc *                      galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"
#include <cstring>

#include "ErrorLogger.h"
#include <sstream>

//int Galprop::electrons_normalize()
int Galprop::electrons_normalize(Particle &particle) //YO/20151123
{
   INFO("Entry");

// identify the primary electrons 
  int ielectrons=-1;
  for(int i=0;i<n_species;i++) 
    if(strcmp(gcr[i].name,"primary_electrons")==0) ielectrons=i;
  if(ielectrons==-1){WARNING("primary electrons not found!"); return 1;}
  ostringstream buf;
  buf<<"  primary electrons found as species #"<<ielectrons;
  INFO(buf.str());

  double  r0=8.5; // solar Galactocentric radius, kpc
// deriving the normalization grid point  
  int ip=(int)(log(galdef.electron_norm_Ekin/galdef.Ekin_min)/log(galdef.Ekin_factor) + 0.5);//IMOS20060420
  int iz=(int)((0.-galdef.z_min)/galdef.dz + 0.5); // Z=0, Galactic plane                      IMOS20060420
  
  double v1,v2,v3,v4,v5,v6;
  if(galdef.n_spatial_dimensions==2)
    {
      int ir=(int)((r0-galdef.r_min)/galdef.dr + 0.5);//IMOS20060420
      buf.str("");
      buf<<"Grid point for normalization: ir r[ir] iz z[iz] ip Ekin[ip] "<<ir<<" " <<gcr[ielectrons].r[ir]<<" " <<iz <<" "<<gcr[ielectrons].z[iz]<<" "<<ip<<" "<< gcr[ielectrons].Ekin[ip];
      INFO(buf.str());
      
      v1=gcr[ielectrons].cr_density.d2[ir  ][iz].s[ip];
      v2=gcr[ielectrons].cr_density.d2[ir+1][iz].s[ip];
      v3=gcr[ielectrons].cr_density.d2[ir  ][iz].s[ip+1];
      v4=gcr[ielectrons].cr_density.d2[ir+1][iz].s[ip+1];
      v5=v1+(r0-gcr[ielectrons].r[ir])/galdef.dr*(v2-v1); // r0 ip
      v6=v3+(r0-gcr[ielectrons].r[ir])/galdef.dr*(v4-v3); // r0 ip+1
    }//n_spatial_dimensions==2
  
  if(galdef.n_spatial_dimensions==3)
    {
      int ix=(int)((r0-galdef.x_min)/galdef.dx + 0.5); //IMOS20060420
      int iy=(int)((0.-galdef.y_min)/galdef.dy + 0.5); //IMOS20060420
      buf.str("");
      buf<<"Grid point for normalization: ix x[ix] iy y[iy] iz z[iz] ip Ekin[ip] "<<ix<<" " <<gcr[ielectrons].x[ix]<<" "<<iy<<" "<<gcr[ielectrons].y[iy]<<" " <<iz <<" "<<gcr[ielectrons].z[iz]<<" "<<ip<<" "<< gcr[ielectrons].Ekin[ip];            //AWS20001121
      INFO(buf.str());
      
      v1=gcr[ielectrons].cr_density.d3[ix  ][iy][iz].s[ip];   //AWS20001121
      v2=gcr[ielectrons].cr_density.d3[ix+1][iy][iz].s[ip];   //AWS20001121
      v3=gcr[ielectrons].cr_density.d3[ix  ][iy][iz].s[ip+1]; //AWS20001121
      v4=gcr[ielectrons].cr_density.d3[ix+1][iy][iz].s[ip+1]; //AWS20001121
      v5=v1+(r0-gcr[ielectrons].x[ix])/galdef.dx*(v2-v1); // r0 ip
      v6=v3+(r0-gcr[ielectrons].x[ix])/galdef.dx*(v4-v3); // r0 ip+1
      
    }//n_spatial_dimensions==3
  
  double vnorm=exp( log(v5)+log(galdef.electron_norm_Ekin/gcr[ielectrons].Ekin[ip])/log(galdef.Ekin_factor)*log(v6/v5) );
  
  buf.str("");
  buf<<"v1 v2 v3 v4 v5 v6 vnorm  "<<v1<<" " <<v2<<" " <<v3 <<" "<<v4 <<" "<<v5<<" "<< v6<<" "<<vnorm;
  INFO(buf.str());

 // normalize electrons
  galdef.electron_source_normalization *= galdef.electron_norm_flux/vnorm; // IMOS20031016
  gcr[ielectrons].cr_density           *=(galdef.electron_norm_flux/vnorm);
  gcr[ielectrons].normalization_factor  =(galdef.electron_norm_flux/vnorm); //AWS20010121

  if(galdef.check_num==1123){
    //double burst_norm = nor_SNR_burstlike(particle);
    cerr << " >>>> BG electron normalization factor: " << gcr[ielectrons].normalization_factor << endl; 
  }
  
  if(galdef.verbose>=2){
    INFO("primary electrons after normalization:");
    gcr[ielectrons].cr_density.print();
  }  

  INFO("Exit");
return 0;
}
