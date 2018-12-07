
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_secondary_proton_source.cc *              galprop package * 2001/05/11 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624

#include"galprop_classes.h"
#include"galprop_internal.h"
#include <cstring>

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// See more details in the header of routine gen_secondary_antiproton_source.cc
// Ref.: Moskalenko I.V. et al. 2002, ApJ 565, 280
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Galprop::gen_secondary_proton_source(Particle &particle)
{
   cout<<"gen_secondary_proton_source"<<endl;
   cout<<"generating "<<particle.name<<" source function for n_spatial_dimensions="
       <<gcr[0].n_spatial_dimensions<<endl;

   int stat=0, iprotons=-1, A1=1, Z2=2, A2=4;
   double PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann;

   if(strcmp(particle.name,"secondary_protons")!=0)
   { cout<<"invalid particle "<<particle.name<<endl; return 2; }

// identify CR protons
   for(int i=n_species-1; i>=0; i--) if(101 == 100*gcr[i].Z+gcr[i].A) {  iprotons=i;   break;  }
   if(iprotons==-1) { cout<<"CR protons not found!"<<endl; return 1; }
   if(strcmp(gcr[iprotons].name,"Hydrogen_1")!=0) { cout<<"CR protons not found!"<<endl; return 1; }
   cout<<"  CR protons found as species #"<<iprotons<<endl;

   //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann)
   for(int ip_sec=0; ip_sec<particle.n_pgrid; ip_sec++)
   {
      for(int ip=0; ip<gcr[iprotons].n_pgrid; ip++)
      {  
         nucleon_cs(galdef.total_cross_section,gcr[iprotons].Ekin[ip]*1.e-3,A1,Z2,A2,&PP_inel,&PA_inel,&aPP_non,&aPA_non,&aPP_ann,&aPA_ann);  // IMOS20010511

         if(galaxy.n_spatial_dimensions==2)
         {
            for(int ir=0; ir<gcr[iprotons].n_rgrid; ir++)
            {
               for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
               {                       // Pbar distribution after scattering ~1/Ep' [TN83,p.235]
                  if(gcr[iprotons].Ekin[ip]<particle.Ekin[ip_sec]) continue; 
                  particle.secondary_source_function.d2[ir][iz].s[ip_sec]+= particle.beta[ip_sec]
		     *(galaxy.n_HI.d2[ir][iz].s[0]+2.0*galaxy.n_H2.d2[ir][iz].s[0]+galaxy.n_HII.d2[ir][iz].s[0])
                     *(PP_inel +galdef.He_H_ratio *PA_inel) *gcr[iprotons].cr_density.d2[ir][iz].s[ip];
               }  //  iz
            }  //  ir
         }

         if(galaxy.n_spatial_dimensions==3)
         {
            for(int ix=0; ix<gcr[iprotons].n_xgrid; ix++)
            {
               for(int iy=0; iy<gcr[iprotons].n_ygrid; iy++)
               {
                  for(int iz=0; iz<gcr[iprotons].n_zgrid; iz++)
                  {                       // Pbar distribution after scattering ~1/Ep' [TN83,p.235]
                     if(gcr[iprotons].Ekin[ip]<particle.Ekin[ip_sec]) continue; 
                     particle.secondary_source_function.d3[ix][iy][iz].s[ip_sec]+= particle.beta[ip_sec]
			*(galaxy.n_HI.d3[ix][iy][iz].s[0]+2.0*galaxy.n_H2.d3[ix][iy][iz].s[0]+galaxy.n_HII.d3[ix][iy][iz].s[0])
                        *(PP_inel +galdef.He_H_ratio *PA_inel) *gcr[iprotons].cr_density.d3[ix][iy][iz].s[ip];
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
      cout<<"   particle.secondary_source_function for "<<particle.name<<endl;
      particle.secondary_source_function.print();
    }
   cout<<" <<<< gen_secondary_proton_source"<<endl;
   return stat;
}
