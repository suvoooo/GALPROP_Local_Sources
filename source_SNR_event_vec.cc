
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * source_SNR_event_vec.cc *                     galprop package * 10/03/2003 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// vectorizable version of source_SNR_event
// CC -c source_SNR_event_vec.C -pvctl,fullmsg -D_BUILTIN_
// need  -D_BUILTIN_ to vectorize fmod()

#define NMAX 1000000
using namespace std;//AWS20050624

#include"galprop_classes.h"
#include"galprop_internal.h"
#include<cstdlib> // exit()
#include <cstring>

int init=0;
int nkxyz;
int nkt;

float SNR_cell_time  [NMAX];
float SNR_cell_phase [NMAX];
double SNR_electron_dg[NMAX]; //AWS20010410
double SNR_nuc_dg     [NMAX]; //AWS20010410

float  p1[NMAX];
float  p2[NMAX];

int Galprop::source_SNR_event_vec(Particle &particle,double t,
                         float *total_source_function_x)
{
int ix,iy,iz,kk;
int ns,i,k;
int stat;
float SNR_rate;

if(galdef.verbose>=1)    cout<<" >>>> source_SNR_event_vec"<<endl;
stat=0;
 
 if(galdef.n_spatial_dimensions==2) return stat;  // not applicable in 2D

 if(galdef.n_spatial_dimensions==3){

   if(init==0){
     cout<<"source_SNR_event_vec: initializing"<<endl;
     nkxyz=particle.n_xgrid*particle.n_ygrid*particle.n_zgrid;
     if(nkxyz >NMAX){
        cout<< "work array too small, need "<<nkxyz
        <<" got only "<<NMAX<<endl; exit(0);
      }
    nkt=particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid;
    init=1;
    kk=0;
    SNR_rate=0.;

    for(    iz=0;iz<particle.n_zgrid;iz++){
    for(    iy=0;iy<particle.n_ygrid;iy++){
    for(    ix=0;ix<particle.n_xgrid;ix++){
  
     
     SNR_cell_time [kk]=galaxy.SNR_cell_time  .d3[ix][iy][iz].s[0];
     SNR_cell_phase[kk]=galaxy.SNR_cell_phase .d3[ix][iy][iz].s[0];

    if (galdef.use_symmetry==0)         SNR_rate+=1./SNR_cell_time [kk];
    //since symmetrical bin at z=0 need to avoid double counting.
    //SNR concentrated at z=0, x and y not critical here
    if (galdef.use_symmetry==1 && iz==0) SNR_rate+=4./SNR_cell_time [kk];
    if (galdef.use_symmetry==1 && iz >0) SNR_rate+=8./SNR_cell_time [kk];

    SNR_electron_dg[kk]=galaxy.SNR_electron_dg.d3[ix][iy][iz].s[0];//AWS20010410
    SNR_nuc_dg     [kk]=galaxy.SNR_nuc_dg.     d3[ix][iy][iz].s[0];//AWS20010410    


     kk++;
    }//iz
    }//iy
    }//ix



    SNR_rate*= 100;// convert from per year to per century
    cout<<"SNR_rate="<<SNR_rate<<" per century"<<endl;

    }// init

   double g_0=0., rigid_br0=0.,                                   // IMOS20031012
          rigid_br=galdef.nuc_rigid_br,                           // IMOS20000607
          g_1=galdef.nuc_g_1,
          g_2=galdef.nuc_g_2;
   if(strcmp(particle.name,"primary_electrons")==0)
   {
      g_0=galdef.electron_g_0;                                    // IMOS20031012
      rigid_br0=galdef.electron_rigid_br0;
      g_1=galdef.electron_g_1;
      rigid_br=galdef.electron_rigid_br;
      g_2=galdef.electron_g_2;
   }
//   cout<<particle.name<<" g_0="<<g_0<<"  rigid_br0= "<<rigid_br0  // IMOS20031012
//                      <<" g_1="<<g_1<<"  rigid_br= " <<rigid_br <<" g_2="<<g_2<<endl;
   double Rsun=8.5; // Galactocentre radius of Sun in kpc
   double cell_volume=galdef.dx * galdef.dy * galdef.dz;
   double factor= galdef.SNR_interval/galdef.SNR_livetime* source_distribution(Rsun,0.0,0.0, galdef.source_model, galdef.source_parameters)/cell_volume;
                  
   int n_SNR=0;

   //for(kk=0;kk<nkt;kk++)total_source_function_x[kk]=0.0; do in protri to vectorize
     
     for(kk=0;kk<nkxyz;kk++){
     
     p1[kk]=fmod (t                     /SNR_cell_time[kk],1.0);
     p2[kk]=fmod((t+galdef.SNR_livetime)/SNR_cell_time[kk],1.0);

     // cout<<" t kk  p1 p2 SNR_cell_time SNR_livetime "<<t<<" "<<kk<<" "<<p1[kk]<<" "<<p2[kk]<<" "<<SNR_cell_time[kk]<<" "<<galdef.SNR_livetime<<endl;
     }

  

     for(kk=0;kk<nkxyz;kk++){

       //      cout<<"kk  p1 p2 SNR_cell_phase "<<kk<<" "<<p1[kk]<<" "<<p2[kk]<<" "<<SNR_cell_phase[kk]<<endl;


     if(p1[kk]<       SNR_cell_phase[kk] && 
        p2[kk]>       SNR_cell_phase[kk]) {

  
          
          for(int ip=0;ip<particle.n_pgrid;ip++){

           double spec_shape;

	    if(particle.rigidity[ip]< rigid_br0)                                   // IMOS20031012
               spec_shape =pow(particle.rigidity[ip]/rigid_br0,-g_0) *pow(rigid_br0/rigid_br,-g_1);
            if(rigid_br0<= particle.rigidity[ip] && particle.rigidity[ip]< rigid_br)
               spec_shape =pow(particle.rigidity[ip]/rigid_br, -g_1);
            if(rigid_br <= particle.rigidity[ip])
               spec_shape =pow(particle.rigidity[ip]/rigid_br, -g_2);

           double spec_dg_ratio;                             //AWS20010410
           if(strcmp(particle.name,"primary_electrons")==0){ //AWS20010410
	     spec_dg_ratio=
	     pow(particle.rigidity[ip]/galdef.SNR_electron_dgpivot,SNR_electron_dg[kk]);

             if(galdef.verbose==-501)// selectable debug
	       cout<<"SNR_electron_dg="<<SNR_electron_dg[kk]
                <<" rigidity="<<particle.rigidity[ip]<<"  spec_dg_ratio="<<spec_dg_ratio<<endl;
	   }

           if(strcmp(particle.name,"primary_electrons")!=0){ //AWS20010410
	     spec_dg_ratio=
	     pow(particle.rigidity[ip]/galdef.SNR_nuc_dgpivot, SNR_nuc_dg[kk] );

             if(galdef.verbose==-501)// selectable debug
               cout<<"SNR_nuc_dg="<<SNR_nuc_dg[kk]
               <<" rigidity="<<particle.rigidity[ip]<<"  spec_dg_ratio="<<spec_dg_ratio<<endl;
	   }

           spec_shape*=spec_dg_ratio;




           total_source_function_x[ip*nkxyz+kk]= factor*particle.primary_abundance*spec_shape   ;

  if(strcmp(particle.name,"primary_electrons")!=0)total_source_function_x[ip*nkxyz+kk]*=( pow(particle.A, g_2-1)* pow(particle.Z,-g_2));

	  }//ip

    n_SNR+=1;

 if(galdef.verbose>=1)    cout<<" source_SNR_event at t="<<t<<" kk x y z= "<<kk<<" "<<galaxy.x[ix]<<" "<<galaxy.y[iy]<<" "<<galaxy.z[iz]<<" SNR_cell_time "<<SNR_cell_time[kk]<<" p1 SNR_cell_phase p2 "<<p1[kk]<<" "<<SNR_cell_phase[kk]<<" "<<p2[kk]<<" total_source_function_x="<<total_source_function_x[0*nkxyz+kk]<<endl;
    //    cout<<" source_SNR_event at t="<<t<<" kk="<<kk<<endl;

     }  
     }//kk

 if(galdef.verbose>=0)   cout<<"source_SNR_event: number of live SNR at this time = "<<n_SNR<<endl;

 }// 3D
   
 if(galdef.verbose>=2)     cout<<" <<<< source_SNR_event_vec"<<endl;
return stat;
}
