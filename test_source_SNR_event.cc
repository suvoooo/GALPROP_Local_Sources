
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * test_source_SNR_event.cc *                    galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"


int Galprop::test_source_SNR_event()
{
int stat;

  cout<<" >>>> test_source_SNR_event"<<endl;
stat=0;

 Particle particle;
 particle.init();// to signal that arrays not yet created

 particle=gcr[0];
 particle.create_transport_arrays(); // creates primary source function array

 double t,dt,tmax;
 dt=10;
 tmax=1.e5;
 for (t=0.;t<tmax;t+=dt){
   source_SNR_event(particle,t);  
   cout<<"source_SNR_event: primary source function for particle "<<particle.name<<endl;
   particle.primary_source_function .print();
}
 
 

     
 

  cout<<" <<<< test_source_SNR_event"<<endl;
return stat;
}
