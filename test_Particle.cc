
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * test_Particle.cc *                            galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"

int Galprop::test_Particle()
{
   cout<<" >>>> test_Particle"<<endl;

   int stat=0;
   Particle particle;  
   char name[100];
   int Z=30, A=60;
   int K_electron=1;                        // AWS20010731
   double t_half=12340.0;
 
   sprintf(name,  "test_particle_%d",A);
   particle.primary_abundance=99.99;

   if(galdef.n_spatial_dimensions==2)
   particle.init(name,Z,A,K_electron,t_half,// AWS20010731
                       galdef.r_min,  galdef.r_max, galdef.dr,  
                       galdef.z_min,  galdef.z_max, galdef.dz,
                       galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
                       galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
                       galdef.p_Ekin_grid);  


   if(galdef.n_spatial_dimensions==3)
   particle.init(name,Z,A,K_electron,t_half,// AWS20010731
                       galdef.x_min,  galdef.x_max, galdef.dx,  
                       galdef.y_min,  galdef.y_max, galdef.dy,
                       galdef.z_min,  galdef.z_max, galdef.dz,
                       galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
                       galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
                       galdef.p_Ekin_grid); 

   particle.cr_density=4321.; 
   particle.print();
 
   cout<<endl<<"   testing : particle_copy=particle"<<endl<<endl;
   Particle particle_copy;
   particle_copy.init();
   particle_copy=particle;
   particle_copy.print();
   particle_copy.cr_density.print();

   cout<<endl<<"   testing : particle_copy2=particle_copy"<<endl<<endl;
   Particle particle_copy2;
   particle_copy2.init();
   particle_copy2=particle_copy;
   particle_copy2.print();
   particle_copy2.cr_density.print();

   cout<<endl<<"   testing : particle_copy2=particle_copy again"<<endl<<endl;
   particle_copy2.cr_density=-1234;
   particle_copy2.cr_density.print();
   particle_copy2=particle_copy;
   particle_copy2.print();
   particle_copy2.cr_density.print();

   cout<<" <<<< test_Particle"<<endl;

   exit(0);
   return stat;
}
