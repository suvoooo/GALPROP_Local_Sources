
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_tertiary_antiproton_source.cc *           galprop package * 2001/05/11 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"
#include <cstring>

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// See more details in the header of routine gen_secondary_antiproton_source.cc
// Ref.: Moskalenko I.V. et al. 2002, ApJ 565, 280
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Galprop::gen_tertiary_antiproton_source(Particle &particle)
{
   cout<<"gen_tertiary_antiproton_source"<<endl;
   cout<<"generating "<<particle.name<<" source function for n_spatial_dimensions="
       <<gcr[0].n_spatial_dimensions<<endl;

   int stat=0, ipbars=-1, A1=-1, Z2=2, A2=4;
   double PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann;

   if(strcmp(particle.name,"tertiary_antiprotons")!=0)
   { cout<<"invalid particle "<<particle.name<<endl; return 2; }

// identify the secondary antiprotons
   for(int i=n_species-1; i>=0; i--)  if(-99 == 100*gcr[i].Z+gcr[i].A) {  ipbars=i;   break;  }
   if(ipbars==-1) { cout<<"secondary antiprotons not found!"<<endl; return 1; }
   if(strcmp(gcr[ipbars].name,"secondary_antiprotons")!=0) { cout<<"secondary antiprotons not found!"<<endl; return 1; }
   cout<<"  secondary antiprotons found as species #"<<ipbars<<endl;

   for(int ip_sec=0; ip_sec<particle.n_pgrid; ip_sec++)
   {
      for(int ip=0; ip<gcr[ipbars].n_pgrid; ip++)
      {  
         nucleon_cs(galdef.total_cross_section,gcr[ipbars].Ekin[ip]*1.e-3,A1,Z2,A2,&PP_inel,&PA_inel,&aPP_non,&aPA_non,&aPP_ann,&aPA_ann);  // IMOS20010511

         if(galaxy.n_spatial_dimensions==2)
         {
            for(int ir=0; ir<gcr[ipbars].n_rgrid; ir++)
            {
               for(int iz=0; iz<gcr[ipbars].n_zgrid; iz++)
               {                       // Pbar distribution after scattering ~1/Ekin' [TN83,p.235]
                  if(gcr[ipbars].Ekin[ip]<particle.Ekin[ip_sec]) continue; 
                  particle.secondary_source_function.d2[ir][iz].s[ip_sec]+= particle.beta[ip_sec]
		     *(galaxy.n_HI.d2[ir][iz].s[0]+2.0*galaxy.n_H2.d2[ir][iz].s[0]+galaxy.n_HII.d2[ir][iz].s[0])
                     *(aPP_non +galdef.He_H_ratio *aPA_non) *gcr[ipbars].cr_density.d2[ir][iz].s[ip];
               }  //  iz
            }  //  ir
         }

         if(galaxy.n_spatial_dimensions==3)
         {
            for(int ix=0; ix<gcr[ipbars].n_xgrid; ix++)
            {
               for(int iy=0; iy<gcr[ipbars].n_ygrid; iy++)
               {
                  for(int iz=0; iz<gcr[ipbars].n_zgrid; iz++)
                  {                    // Pbar distribution after scattering ~1/Ekin' [TN83,p.235]
                     if(gcr[ipbars].Ekin[ip]<particle.Ekin[ip_sec]) continue; 
                     particle.secondary_source_function.d3[ix][iy][iz].s[ip_sec]+= particle.beta[ip_sec]
			*(galaxy.n_HI.d3[ix][iy][iz].s[0]+2.0*galaxy.n_H2.d3[ix][iy][iz].s[0]+galaxy.n_HII.d3[ix][iy][iz].s[0])
                        *(aPP_non +galdef.He_H_ratio *aPA_non) *gcr[ipbars].cr_density.d3[ix][iy][iz].s[ip];
                  }  //  iz
               }  //  iy
            }  //  ix
         }
      }  //  ip
   }  //  iEgamma
 
   double factor=1.e-27 *C *log(galdef.Ekin_factor); // transformation mb -> cm2 and constant factors
   particle.secondary_source_function *= factor;

   if(galdef.verbose>=2)
   {
      cout<<"   particle.tertiary_source_function for "<<particle.name<<endl;
      particle.secondary_source_function.print();
    }
   cout<<" <<<< gen_tertiary_antiproton_source"<<endl;
   return stat;
}
