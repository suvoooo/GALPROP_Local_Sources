//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * source_SNR_event.cc *                         galprop package * 10/12/2003 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include"galprop_classes.h"
#include"galprop_internal.h"

#include <ErrorLogger.h>

#include <cstring>

using namespace std;//AWS20050624

int Galprop::source_SNR_event(Particle& particle, double t) {

  INFO("Entry");

  //if (galdef.verbose>=2)    cout<<" >>>> source_SNR_event"<<endl;
  int stat = 0;
 
  if (2 == galdef.n_spatial_dimensions) 
    return stat;  // not applicable in 2D

  //if (3 == galdef.n_spatial_dimensions) { 

  double p1, p2;
  
  particle.primary_source_function=0.0;
  
  double g_0=0., rigid_br0=0.,                                   // IMOS20031012
    rigid_br=galdef.nuc_rigid_br,                           // IMOS20000607
    g_1=galdef.nuc_g_1,
    g_2=galdef.nuc_g_2;
  
  if ("primary_electrons" == particle.name) { //strcmp(particle.name,"primary_electrons")==0)
    //{
      
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
  double factor= galdef.SNR_interval/galdef.SNR_livetime*source_distribution(Rsun,0,0, galdef.source_model,galdef.source_parameters)/cell_volume;
  
  int n_SNR=0;
  
  for (int ix=0;ix<particle.n_xgrid;ix++){
    for (int iy=0;iy<particle.n_ygrid;iy++){
      for (int iz=0;iz<particle.n_zgrid;iz++){
		
	p1 = fmod(t/galaxy.SNR_cell_time.d3[ix][iy][iz].s[0], 1.);
	p2 = fmod((t + galdef.SNR_livetime)/galaxy.SNR_cell_time.d3[ix][iy][iz].s[0], 1.);
		
	if (p1 < galaxy.SNR_cell_phase.d3[ix][iy][iz].s[0] && 
	    p2 > galaxy.SNR_cell_phase.d3[ix][iy][iz].s[0]) {
	  	  
          for (int ip = 0; ip < particle.n_pgrid; ++ip) {
	    
	    double spec_shape;
	    
	    if (particle.rigidity[ip]< rigid_br0) // IMOS20031012
	      spec_shape = pow(particle.rigidity[ip]/rigid_br0, -g_0)*pow(rigid_br0/rigid_br, -g_1);

            if (rigid_br0 <= particle.rigidity[ip] && 
		particle.rigidity[ip]< rigid_br)
               spec_shape = pow(particle.rigidity[ip]/rigid_br, -g_1);

            if (rigid_br <= particle.rigidity[ip])
	      spec_shape = pow(particle.rigidity[ip]/rigid_br, -g_2);

	    double spec_dg_ratio;                             //AWS20010410
	    
	    if ("primary_electrons" == particle.name) { //strcmp(particle.name,"primary_electrons")==0){ //AWS20010410
	    
	      spec_dg_ratio=
		pow(particle.rigidity[ip]/galdef.SNR_electron_dgpivot, 1.*galaxy.SNR_electron_dg.d3[ix][iy][iz].s[0]);

	      if(galdef.verbose==-501)// selectable debug
		cout<<"SNR_electron_dg="<<galaxy.SNR_electron_dg.d3[ix][iy][iz].s[0]
		    <<" rigidity="<<particle.rigidity[ip]<<"  spec_dg_ratio="<<spec_dg_ratio<<endl;
	    }

	    if ("primary_electrons" != particle.name) { //strcmp(particle.name,"primary_electrons")!=0){ //AWS20010410
	      spec_dg_ratio=
		pow(particle.rigidity[ip]/galdef.SNR_nuc_dgpivot, 1.*galaxy.SNR_nuc_dg.d3[ix][iy][iz].s[0]);

	      if(galdef.verbose==-501)// selectable debug
		cout<<"SNR_nuc_dg="<<galaxy.SNR_nuc_dg.d3[ix][iy][iz].s[0] 
		    <<" rigidity="<<particle.rigidity[ip]<<"  spec_dg_ratio="<<spec_dg_ratio<<endl;
	    }
	    
	    spec_shape*=spec_dg_ratio;
	    
	    particle.primary_source_function.d3[ix][iy][iz].s[ip]=factor*particle.primary_abundance*spec_shape   ;

	  }//ip

	  n_SNR += 1;

	  cout<<" source_SNR_event at t="<<t<<" x y z= "<<galaxy.x[ix]<<" "<<galaxy.y[iy]<<" "<<galaxy.z[iz]<<" SNR_cell_time "<<galaxy.SNR_cell_time .d3[ix][iy][iz].s[0]<<" p1 SNR_cell_phase p2 "<<p1<<" "<<galaxy.SNR_cell_phase.d3[ix][iy][iz].s[0]<<" "<<p2<<" primary_source_function="<<particle.primary_source_function.d3[ix][iy][iz].s[0]<<endl; //AWS20001121
     
	}
	
      }//iz
   
    }//iy
   
  }//ix

  if ("primary_electrons" != particle.name) //strcmp(particle.name,"primary_electrons")!=0)
    particle.primary_source_function *= (pow(particle.A, g_2-1)*pow(particle.Z,-g_2)); // (Note 1)

  cout<<"source_SNR_event: number of live SNR at this time = "<<n_SNR<<endl;
  
  //}// 3D
    
  if(galdef.verbose>=1){
    
    cout<<"source_SNR_event: primary source function for particle "<<particle.name<<endl;
    particle.primary_source_function .print();
    
  }

  INFO("Exit");

  return stat;

}
