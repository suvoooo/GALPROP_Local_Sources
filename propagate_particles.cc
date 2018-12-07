
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * propagate_particles.cc *                      galprop package * 1/29/2003 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"
#include <cstring>
#include <algorithm>

#include <ErrorLogger.h>

int Galprop::propagate_particles() { //AWS20050816

  //cout<<">>>>propagate_particles"<<endl;

  INFO("Entry");
  
  int i,net_iter;           //IMOS20032901
  Particle particle;
  
  if(galdef.warm_start==0) for(i=0; i<n_species; i++) gcr[i].cr_density=0.0; //AWS20010121
  if(galdef.warm_start==1) read_gcr();                                       //AWS20010121  
  
  particle.init();  // to signal that arrays not yet allocated
  
  particle=gcr[0];
  particle.create_transport_arrays();

  int network_iterations = max(galdef.network_iterations, galdef.network_iter_compl);
  network_iterations = max(network_iterations, galdef.network_iter_sec);
  
  for (int r_net_iter=network_iterations; r_net_iter>=1; --r_net_iter) { //IMOS20030129
   
     net_iter = network_iterations-r_net_iter+1;
    //cout<<"    Network iteration "<<net_iter<<endl;   //IMOS20030129
    
    //Gulli20100324
    //Network iterations only used to get damping correct and B/C ratio
    //(hopefully)
    //Limit secondary low A (<=1) generation to last iteration
    //If more than two iterations, only propagate protons except for last two
    //iterations.
    for (i = n_species-1; i>-1; i--) {
      
      particle=gcr[i];
      if ( r_net_iter > galdef.network_iter_compl && ! ( particle.Z == 1 && particle.A == 1 ) )
	 continue;
      
      create_transport_arrays(particle);

      // cout<<"particle.Dxx:"<<endl;particle.Dxx.print();
      
      //cout<<">>>>generating secondary source for particle "<<i<<endl;
      if ( r_net_iter <= galdef.network_iter_sec || ( particle.A >= 1 && particle.Z != -1 ) )
	 if (0 != gen_secondary_source(particle)) 
	    return 1;
      
      if(galdef.verbose==10) particle.secondary_source_function.print();
      //if(galdef.verbose>= 1) particle.print();
      if(galdef.verbose>= 0) {

	ostringstream lvl0Buf;
	lvl0Buf << "\n Network iteration "<<net_iter   //IMOS20030129
		<<" species "<<i<<" "<<gcr[i].name<<" (Z,A) = ("<<gcr[i].Z
		<<","<<gcr[i].A<<")"; 
	INFO(lvl0Buf.str());
      
      }

      if(propel(particle)!=0) return 1;
      // cout<<"============particle.cr_density:"<<endl;particle.cr_density.print();
      
      gcr[i]=particle;
      if (galdef.verbose >= 1) {

	ostringstream lvl1Buf;
	lvl1Buf <<"Network iteration "<<net_iter     //IMOS20030129
		<<" species "<<i<<"  "<<gcr[i].name;
	INFO(lvl1Buf.str());
	gcr[i].cr_density.print();
	particle.cr_density.print();
	particle.print();
      }
    } 
    
    // test of electron propagation vs analytical calculations (only runs for galdef.DM_int0>0)  IMOS20061030
    if(galdef.DM_int0==99)
      {
	cout<<" ***** Analytical test of electron propagation  *****"<<endl
	    <<" ***** results are stored in DM_positrons array *****"<<endl;
	
	int iDM_positrons=-1;// identify test DM_positrons array
	for(i=0; i<n_species; i++)  
	  if(strcmp(gcr[i].name,"DM_positrons")==0)
	    {
	      iDM_positrons=i;
	      cout<<"  DM_positrons found as species #"<<iDM_positrons<<endl;
	    }
	if(iDM_positrons==-1) { cout<<"  DM_positrons not found!"<<endl;  return 1; }
	
	// assigning analytical solution to DM_positrons array
	double elossconst=32./9.*Pi*pow(Rele/Mele*1.e3,2)*C*galdef.DM_double7*1.e-9;// 1/(GeV s)
	int iz1=0, iz2=gcr[0].n_zgrid-1;
	if(!galdef.output_gcr_full) iz1=iz2=(int)(1.e-6-galdef.z_min/galdef.dz);//z=0,Galactic plane
	for(int ir=0; ir<gcr[0].n_rgrid; cout<<">> r= "<<galaxy.r[ir]<<endl, ir++)
	  {
	    for(int iz=iz1; iz<iz2+1; cout<<" "<<galaxy.z[iz]<<":", iz++)
	      {
		for(int ip=0; ip<gcr[iDM_positrons].n_pgrid; ip++)
		  {
		    gcr[iDM_positrons].cr_density.d2[ir]    [iz].s[ip]= C/(4.*Pi)/1.e3
		      *eprop(gcr[iDM_positrons].Ekin[ip]/1.e3,galaxy.r[ir],galaxy.z[iz],
			     galdef.r_max,galdef.DM_double6,galdef.z_max,elossconst,
			     galdef.electron_g_0,galdef.D0_xx*pow(kpc2cm,-2),galdef.D_g_1);
		  }
	      }
	  }
	cout<<" *****   End of test of electron propagation    *****"<<endl;
      }
    
    // convert density per momentum to flux per (KE/nucleon) for output NOT done in store_gcr
    // only for nuclei
    // for(i=0; i<n_species; i++) if(gcr[i].A!=0) gcr[i].cr_density*=gcr[i].A;
    
    if(galdef.proton_norm_flux > 0) nuclei_normalize();                                                //IMOS20030129
    //if(galdef.primary_electrons==1 && r_net_iter <= galdef.network_iter_compl && galdef.electron_norm_flux > 0) electrons_normalize();             //IMOS20030129
    if(galdef.primary_electrons==1 && r_net_iter <= galdef.network_iter_compl && galdef.electron_norm_flux > 0) electrons_normalize(particle); //YO/20151123
    //cout << '\a';                                                      //IMOS20030129
  } //nuc_iter                                                          //IMOS20030129
  
  //Gulli20070810 clean up the gcr transport arrays
  ostringstream cleanBuf;
  cleanBuf << "Cleaning up for: "<<n_species<<" number of species.";
  INFO(cleanBuf.str());
  //for(i=0; i<n_species; i++)
  // cout<< "Deleting transport arrays:" << i << endl;
  // gcr[i].delete_transport_arrays();
  //particle.delete_arrays();
  
  INFO("Exit");
  //cout<<"<<<<propagate_particles"<<endl;
  return 0;
}


