
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * create_SNR.cc *                               galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

using namespace std;//AWS20050624
#include<cstdlib>
#include"galprop_classes.h"
#include"galprop_internal.h"


int Galprop::create_SNR()
{
int stat;

  cout<<" >>>> create_SNR"<<endl;
stat=0;


 
 
 if(galdef.n_spatial_dimensions==2) return stat;  // not applicable in 2D

 if(galdef.n_spatial_dimensions==3){

 
  
   double Rsun = 8.5;
   double cell_volume=galdef.dx * galdef.dy * galdef.dz;

   unsigned seed;
   seed=1234;
   srand(seed); // eventually use a galdef parameter

   for(int ix=0;ix<galaxy.n_xgrid;ix++){
   for(int iy=0;iy<galaxy.n_ygrid;iy++){
   for(int iz=0;iz<galaxy.n_zgrid;iz++){
  
  
     galaxy.SNR_cell_time .d3[ix][iy][iz].s[0] =
       galdef.SNR_interval* 
       source_distribution(Rsun, 0, 0, galdef.source_model, galdef.source_parameters)/
       (source_distribution(galaxy.x[ix],
			    galaxy.y[iy],
			    galaxy.z[iz],
			    galdef.source_model,
                            galdef.source_parameters)+1.e-30)/cell_volume;

     
     //    cout<<" create_SNR "<<galaxy.x[ix]<<" "<<galaxy.y[iy]<<" "<<galaxy.z[iz]<<" "<<galdef.SNR_interval<<" "<<source_distribution(galaxy.x[ix],galaxy.y[iy],galaxy.z[iz])<<" "<<source_distribution(Rsun,               0.0,         0.0  )<<" "<<cell_volume<<" "<<galaxy.SNR_cell_time .d3[ix][iy][iz].s[0]<<endl;


     
     galaxy.SNR_cell_phase.d3[ix][iy][iz].s[0]=float(rand())/RAND_MAX;


     // Gaussian distributed source index delta                                       AWS20010410
     galaxy.SNR_electron_dg.d3[ix][iy][iz].s[0]=gauss(0.0,galdef.SNR_electron_sdg); //AWS20010410
     galaxy.SNR_nuc_dg     .d3[ix][iy][iz].s[0]=gauss(0.0,galdef.SNR_nuc_sdg);      //AWS20010410

     }
   }
   }
   }
   
 
 if(galdef.verbose>=2 || galdef.verbose==-401){ // selectable debug AWS20010410
  cout<<"galaxy.SNR_cell_time :"<<endl;galaxy.SNR_cell_time .print();
  cout<<"galaxy.SNR_cell_phase:"<<endl;galaxy.SNR_cell_phase.print();
                      
  cout<<"galaxy.SNR_electron_dg:"<<endl;galaxy.SNR_electron_dg.print(); //AWS20010410
  cout<<"galaxy.SNR_nuc_dg:     "<<endl;galaxy.SNR_nuc_dg     .print(); //AWS20010410
 }

  cout<<" <<<< create_SNR"<<endl;
return stat;
}
