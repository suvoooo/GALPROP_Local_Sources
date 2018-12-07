//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_secondary_source.cc *                     galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <valarray>

using namespace std;//AWS20050624

#include"galprop_classes.h"
#include"galprop_internal.h"

#include <ErrorLogger.h>

int Galprop::gen_secondary_source(Particle& particle) {

  INFO("Entry");

  //AWS20091009 fixed the changes introduced by TAP in r338: particle.name is a char array not a string

  //string particle_name;                                             //AWS20091009
  //particle_name = particle.name; // transfer from char to string    //AWS20091009

  const string particle_name = particle.name; 

  if (galdef.verbose >= 1) {

    ostringstream buf;
    buf << "Generating for (Z, A, K_electron) = (" 
	<< particle.Z << ", "
	<< particle.A << ", "
	<< particle.K_electron << ")";
    INFO(buf.str());

  }

  Distribution gas_times_cross_section;
  Distribution decay_source;                           //AWS20010806
  int i,ip, Z1,A1, galdef_network_par = 0;
  
  valarray<double> cross_section(0., gcr[0].n_pgrid);
  double t_half;

  //double* cross_section, t_half;
  //cross_section=new double[gcr[0].n_pgrid];

  if (2 == galdef.n_spatial_dimensions) {                   //AWS20010806
   
    gas_times_cross_section.init(gcr[0].n_rgrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
    decay_source.init(gcr[0].n_rgrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
  
  }
  
  if (3 == galdef.n_spatial_dimensions) {                   //AWS20010806
    
    gas_times_cross_section.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
    decay_source.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);          

  }

  particle.secondary_source_function = 0;

  if (particle_name == "primary_electrons") { // Primary

    INFO("Exit");
    return 0;

  }

  if (particle_name == "secondary_electrons" || 
      particle_name == "secondary_positrons") {

    gen_secondary_positron_source(particle);
    INFO("Exit");
    return 0;

  }

  if (particle_name == "knock_on_electrons") {

    gen_knock_on_electron_source(particle);
    INFO("Exit");
    return 0;

  }

  if (particle_name == "secondary_protons") {

    gen_secondary_proton_source(particle);
    INFO("Exit");
    return 0;

  }

  if (particle_name == "secondary_antiprotons") {

    gen_secondary_antiproton_source(particle);
    INFO("Exit");
    return 0;

  }

  if (particle_name == "tertiary_antiprotons") {

    gen_tertiary_antiproton_source(particle);
    INFO("Exit");
    return 0;

  }

  /*if(strcmp(particle.name,"primary_electrons"    )==0)       // since primary
    {if(galdef.verbose>=1) cout<<" <<<<gen_secondary_source"<<endl  ;  return 0;}

   if(strcmp(particle.name,"secondary_electrons"  )==0)       // secondary e-
      {gen_secondary_positron_source  (particle);  return 0;}

   if(strcmp(particle.name,"knock_on_electrons"  )==0)        // knock-on e-
      {gen_knock_on_electron_source  (particle);  return 0;}

   if(strcmp(particle.name,"secondary_positrons"  )==0)       // secondary e+
      {gen_secondary_positron_source  (particle);  return 0;}

   if(strcmp(particle.name,"secondary_antiprotons")==0)       // secondary p-
      {gen_secondary_antiproton_source(particle);  return 0;}

   if(strcmp(particle.name,"tertiary_antiprotons")==0)        // tertiary p- IMOS20000605
      {gen_tertiary_antiproton_source(particle);  return 0;}  //             IMOS20000605

   if(strcmp(particle.name,"secondary_protons")==0)           // secondary protons IMOS20000605
      {gen_secondary_proton_source(particle);  return 0;}     //                   IMOS20000605
  */

  // particles by DM annihilation BiXJ 2007/2/2
  if (particle_name == "primary_DM_positron" ||
      particle_name == "primary_DM_electron" ||
      particle_name == "primary_DM_antiproton") {

    int model_num = galdef.TD_SNR_model;
    switch(model_num){ //YO/20151224
    case -1:
      gen_DM_annihilation(particle);
      break;
    case -2:
      gen_galactic_Pulsar(particle);
      break;
    }
    
    INFO("Exit");
    return 0;

  }

  /*if(strcmp(particle.name,"primary_DM_positron"  )==0)              // DM e+  BiXJ 2007/2/2
    {gen_DM_annihilation(particle);  return 0;}
  
  if(strcmp(particle.name,"primary_DM_electron"  )==0)              // DM e- BiXJ 2007/2/2
    {gen_DM_annihilation(particle);  return 0;}
  
  if(strcmp(particle.name,"primary_DM_antiproton")==0)              // DM p- BiXJ 2002/2/2
    {gen_DM_annihilation(particle);  return 0;}
  */

// positron by time dependent injection YO/20140220
  /*if (particle_name == "TD_injection") {

    gen_TD_injection(particle);
    INFO("Exit");
    return 0;

    }*/
 
  // DM source

  if (particle_name == "DM_positrons" ||
      particle_name == "DM_electrons" ||
      particle_name == "DM_antiprotons") {

    gen_DM_source(particle);
    INFO("Exit");
    return 0;

  }

  /*if(strcmp(particle.name,"DM_positrons"  )==0)              // DM e+:  IMOS20050912
    {gen_DM_source(particle);  return 0;}
  
  if(strcmp(particle.name,"DM_electrons"  )==0)              // DM e-:  IMOS20050912
    {gen_DM_source(particle);  return 0;}

  if(strcmp(particle.name,"DM_antiprotons")==0)              // DM p-:  IMOS20050912
    {gen_DM_source(particle);  return 0;}
  */

  // FRAGMENTATION SOURCE
  
  for (i = 0; i < n_species; ++i) {

    cross_section = 0;
    
    //for (ip = 0; ip < gcr[0].n_pgrid; cross_section[ip++] = 0.);            // IMOS20010816
      
    if (gcr[i].A > particle.A && !particle.K_electron) {
	
      decayed_cross_sections(gcr[i].Z,gcr[i].A, particle.Z,particle.A, particle.Ekin, particle.n_pgrid, &cross_section[0]);
      
      if (galdef.verbose>=2) { 
	  
	ostringstream buf;
	buf << "Cross_section for Z, A = " <<gcr[i].Z << ", " <<gcr[i].A << "-> " << particle.Z << ", " <<particle.A << " :";
	for (ip = 0; ip < gcr[0].n_pgrid; buf << cross_section[ip++] << " "); 
	INFO(buf.str());
	
      }
	
      // apply factors as shown in galprop notes
      
      //for (ip = 0; ip < gcr[0].n_pgrid; ++ip) 
      //cross_section[ip] *= double(gcr[i].A)/particle.A;
      
      cross_section *= double(gcr[i].A)/particle.A;
      
      if (galdef.verbose >= 2) {
	
	ostringstream buf;
	buf << "Cross_section with Aprim/Asec factors: ";
	for (ip = 0; ip < gcr[0].n_pgrid; buf << cross_section[ip++] << " ");
	INFO(buf.str());
	
      }
      
      if (galdef.verbose == -601 || galdef.verbose == 10) { // selectable debug
	cout << "Progenitor density of " << gcr[i].name << " for " << particle.name << ":\n";
	gcr[i].cr_density.print();
      }
      
    }  //  gcr.A>particle.A
    
    if (1 == galdef.K_capture) { //AWS20010828
      
      // ELECTRON ATTACHMENT or STRIPPING SOURCE   IMOS20010816
      
      double attach_H, strip_H, attach_He, strip_He;
      int zH = 1, zHe = 2;
      
      // Source from isotopes quickly decaying after EC 
      if (gcr[i].A == particle.A && gcr[i].Z != particle.Z)
	if (nucdata(galdef_network_par,
		    gcr[i].Z,
		    gcr[i].A,
		    1,
		    particle.Z,
		    particle.A, 
		    &Z1,
		    &A1,
		    &t_half)) {
	  
	  if (t_half > 0. && t_half/year2sec < galdef.t_half_limit) {
	    
	    for (ip = 0; ip<gcr[0].n_pgrid; ++ip) {
	      
	      Kcapture_cs(gcr[i].Ekin[ip],gcr[i].Z,zH, &attach_H ,&strip_H );
	      Kcapture_cs(gcr[i].Ekin[ip],gcr[i].Z,zHe,&attach_He,&strip_He);
	      cross_section[ip] = attach_H + galdef.He_H_ratio*attach_He;
	      
	    }
	    
	    if (galdef.verbose==-602 ||galdef.verbose>=1) { // selectable debug                  //AWS20010828
	      
	      ostringstream buf;
	      buf << "gen_secondary_source: electron attachment and fast decay cross_section for Z,A ="
		  <<gcr[i].Z<<","<<gcr[i].A<<" K_electron="<<gcr[i].K_electron
		  <<"->>"
		  <<particle.Z<<","<<particle.A<<"  K_electron="<<particle.K_electron
		  <<" t_half="<<t_half/year2sec<<" yr : ";
	      for(ip=0; ip<gcr[0].n_pgrid; buf << cross_section[ip++]<<" "); 
	      INFO(buf.str());
	      
	    }
	    
	    if (galdef.verbose==-601 ||galdef.verbose==10) { // selectable debug
	      
	      cout<<" progenitor density of "<<gcr[i].name<<" for "<<particle.name<<":\n";
	      gcr[i].cr_density.print();
              
	    }
	    
	  }
	  
	}
      
      // Source from electron attachment or stripping (for slow EC decays) 
      if (100*gcr[i].Z+gcr[i].A == 100*particle.Z+particle.A)
	if (gcr[i].K_electron != particle.K_electron) {
	  
	  for (ip = 0; ip < gcr[0].n_pgrid; ++ip) {
	    
	    Kcapture_cs(gcr[i].Ekin[ip],gcr[i].Z,zH, &attach_H ,&strip_H );
	    Kcapture_cs(gcr[i].Ekin[ip],gcr[i].Z,zHe,&attach_He,&strip_He);
	    
	    if(particle.K_electron) 
	      cross_section[ip] = attach_H + galdef.He_H_ratio*attach_He;
	    else                    
	      cross_section[ip] = strip_H + galdef.He_H_ratio*strip_He ;
	    
	  }
	  
	  if (galdef.verbose==-602||galdef.verbose>=1) { // selectable debug                  //AWS20010828
	    
	    ostringstream buf;
	    
	    if (particle.K_electron) 
	      buf<<"gen_secondary_source: electron attachment for Z,A =";
	    else
	      buf<<"gen_secondary_source: stripping cross_section for Z,A =";
	    
	    buf <<gcr[i].Z<<","<<gcr[i].A<<" K_electron="<<gcr[i].K_electron<<"->"
		<<particle.Z<<","<<particle.A<<"  K_electron="<<particle.K_electron<<" :";
	    for (ip=0; ip<gcr[0].n_pgrid; buf<<cross_section[ip++]<<" "); 
	    
	    INFO(buf.str());
	    
	  }
	  
	  if(galdef.verbose==-601 ||galdef.verbose==10) {// selectable debug   
	    
	    cout<<" progenitor density of "<<gcr[i].name<<" for "<<particle.name<<":\n";
	    gcr[i].cr_density.print();
	    
	  }
	  
	} // gcr[i].K_electron != particle.K_electron
	
    } //galdef.K_capture==1
      
    cross_section *= 1e-27; // mb -> cm^2
    
    //for (ip = 0; ip < gcr[0].n_pgrid; cross_section[ip++] *= 1.e-27); // mb -> cm^2
    
    if (2 == galdef.n_spatial_dimensions)
      for(int ir=0; ir<gcr[0].n_rgrid; ir++)
	for(int iz=0; iz<gcr[0].n_zgrid; iz++)
	  for(int ip=0; ip<gcr[0].n_pgrid; ip++)
	    gas_times_cross_section.d2[ir][iz].s[ip] = 
	      cross_section[ip]*gcr[i].beta[ip]*C*(galaxy.n_HI.d2[ir][iz].s[0] + 2.*galaxy.n_H2.d2[ir][iz].s[0] + galaxy.n_HII.d2[ir][iz].s[0]);
    
    if (3 == galdef.n_spatial_dimensions)
      for(int ix=0; ix<gcr[0].n_xgrid; ix++)
	for(int iy=0; iy<gcr[0].n_ygrid; iy++)
	  for(int iz=0; iz<gcr[0].n_zgrid; iz++)
	    for(int ip=0; ip<gcr[0].n_pgrid; ip++)
	      gas_times_cross_section.d3[ix][iy][iz].s[ip] = 
		cross_section[ip]*gcr[i].beta[ip]*C*(galaxy.n_HI.d3[ix][iy][iz].s[0] + 2.*galaxy.n_H2.d3[ix][iy][iz].s[0] + galaxy.n_HII.d3[ix][iy][iz].s[0]);
    
    gas_times_cross_section *= gcr[i].cr_density;
    particle.secondary_source_function += gas_times_cross_section;
    
    if (galdef.verbose==-603 ||galdef.verbose==10) {// selectable debug   
      
      cout<<" gas times cross section times progenitor density of "<<gcr[i].name<<" for "<<particle.name<<":\n";
      gas_times_cross_section.print();
      
    }
    
    // RADIOACTIVE DECAY SOURCE (beta+/-, EC) AWS20010806 IMOS20010816
    
    //double branching_ratio;
    
    if (gcr[i].A >= particle.A) {
      
      if (100*gcr[i].Z+gcr[i].A != 100*particle.Z+particle.A) {
	
	const double branching_ratio = nucdata(galdef_network_par,gcr[i].Z,  
					       gcr[i].A, gcr[i].K_electron,
					       particle.Z,particle.A, 
					       &Z1, &A1, &t_half);
	
	if (branching_ratio != 0 && t_half != 0. && 100*Z1+A1 == 0) {
	  
	  if(galdef.verbose==-602||galdef.verbose>=1) { // selectable debug                  //AWS20010828
	    
	    ostringstream buf;
	    buf <<"gen_secondary_source: (Z,A,K_electron) "
		<< gcr[i].Z<<" "<<  gcr[i].A<<" "<<  gcr[i].K_electron<<"->"
		<<particle.Z<<" "<<particle.A
		<<" t_half="<<t_half/year2sec<<" yr  branching ratio="<<branching_ratio;
	    INFO(buf.str());
	    
	  }
	  
	  if (2 == galdef.n_spatial_dimensions)
	    for(int ir=0; ir<gcr[0].n_rgrid; ir++)
	      for(int iz=0; iz<gcr[0].n_zgrid; iz++)
		for(int ip=0; ip<gcr[0].n_pgrid; ip++)
		  decay_source.d2[ir] [iz].s[ip] = 
		    branching_ratio*gcr[i].cr_density.d2[ir][iz].s[ip]
		    /(gcr[i].gamma[ip]*t_half/log(2.0));
	  
	  if (3 == galdef.n_spatial_dimensions)
	    for(int ix=0; ix<gcr[0].n_xgrid; ix++)
	      for(int iy=0; iy<gcr[0].n_ygrid; iy++)
		for(int iz=0; iz<gcr[0].n_zgrid; iz++)
		  for(int ip=0; ip<gcr[0].n_pgrid; ip++)
		    decay_source.d3[ix][iy][iz].s[ip] = 
		      branching_ratio*gcr[i].cr_density.d3[ix][iy][iz].s[ip]
		      /(gcr[i].gamma[ip]*t_half/log(2.0));
	  
	  particle.secondary_source_function += decay_source;
          
	} // branching_ratio != 0 && t_half != 0. && 100*Z1+A1 == 0
        
      } // 100*gcr[i].Z+gcr[i].A != 100*particle.Z+particle.A
      
    } // gcr[i].A >= particle.A
    
  }  //  species
     
  if(galdef.verbose==-603 ||galdef.verbose==10) { // selectable debug
      
    cout<<" secondary source for "<<particle.name<<":\n";
    particle.secondary_source_function.print();
    
  }
  
  if(galdef.verbose==10) 
    particle.print();
  
  gas_times_cross_section.delete_array();
  decay_source.delete_array(); // AWS20010806
  //delete[] cross_section;
  
  INFO("Exit");
  
  return 0;
    
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
int Galprop::gen_TD_secondary_source(Particle& particle) { //YO/20140221

  //INFO("Entry");
  //AWS20091009 fixed the changes introduced by TAP in r338: particle.name is a char array not a string

  //string particle_name;                                             //AWS20091009
  //particle_name = particle.name; // transfer from char to string    //AWS20091009

  const string particle_name = particle.name; 

  if (galdef.verbose >= 1) {

    ostringstream buf;
    buf << "Generating for (Z, A, K_electron) = (" 
	<< particle.Z << ", "
	<< particle.A << ", "
	<< particle.K_electron << ")";
    INFO(buf.str());

  }

  Distribution gas_times_cross_section;
  Distribution decay_source;                           //AWS20010806
  int i,ip, Z1,A1, galdef_network_par = 0;
  
  valarray<double> cross_section(0., gcr[0].n_pgrid);
  double t_half;

  if (2 == galdef.n_spatial_dimensions) {                   //AWS20010806
   
    gas_times_cross_section.init(gcr[0].n_rgrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
    decay_source.init(gcr[0].n_rgrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
  
  }
  
  if (3 == galdef.n_spatial_dimensions) {                   //AWS20010806
    
    gas_times_cross_section.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
    decay_source.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);          

  }

  particle.secondary_source_function = 0;

  if (particle_name == "primary_electrons") { // Primary

    INFO("Exit");
    return 0;

  }

  if (particle_name == "secondary_electrons" || 
      particle_name == "secondary_positrons") {

    gen_secondary_positron_source(particle);
    INFO("Exit");
    return 0;

  }

  if (particle_name == "knock_on_electrons") {

    gen_knock_on_electron_source(particle);
    INFO("Exit");
    return 0;

  }

  if (particle_name == "secondary_protons") {

    gen_secondary_proton_source(particle);
    INFO("Exit");
    return 0;

  }

  if (particle_name == "secondary_antiprotons") {

    gen_secondary_antiproton_source(particle);
    INFO("Exit");
    return 0;

  }

  if (particle_name == "tertiary_antiprotons") {

    gen_tertiary_antiproton_source(particle);
    INFO("Exit");
    return 0;

  }

  // particles by DM annihilation BiXJ 2007/2/2
  if (particle_name == "primary_DM_positron" ||
      particle_name == "primary_DM_electron" ||
      particle_name == "primary_DM_antiproton") {

    gen_DM_annihilation(particle);
    INFO("Exit");
    return 0;

  }

// positron by time dependent injection YO/20140220 -> 20140408 changed
  if (particle_name == "primary_TD_positron" ||
      particle_name == "primary_TD_electron") {

    gen_TD_injection(particle);
    //INFO("Exit"); //YO/20140305
    return 0;

    }
 
// DM source

  if (particle_name == "DM_positrons" ||
      particle_name == "DM_electrons" ||
      particle_name == "DM_antiprotons") {

    gen_DM_source(particle);
    INFO("Exit");
    return 0;

  }

// FRAGMENTATION SOURCE
  
  for (i = 0; i < n_species; ++i) {

    cross_section = 0;
    
    if (gcr[i].A > particle.A && !particle.K_electron) {
	
      decayed_cross_sections(gcr[i].Z,gcr[i].A, particle.Z,particle.A, particle.Ekin, particle.n_pgrid, &cross_section[0]);
      
      if (galdef.verbose>=2) { 
	  
	ostringstream buf;
	buf << "Cross_section for Z, A = " <<gcr[i].Z << ", " <<gcr[i].A << "-> " << particle.Z << ", " <<particle.A << " :";
	for (ip = 0; ip < gcr[0].n_pgrid; buf << cross_section[ip++] << " "); 
	INFO(buf.str());
	
      }
	
      // apply factors as shown in galprop notes
      cross_section *= double(gcr[i].A)/particle.A;
      
      if (galdef.verbose >= 2) {
	
	ostringstream buf;
	buf << "Cross_section with Aprim/Asec factors: ";
	for (ip = 0; ip < gcr[0].n_pgrid; buf << cross_section[ip++] << " ");
	INFO(buf.str());
	
      }
      
      if (galdef.verbose == -601 || galdef.verbose == 10) { // selectable debug
	cout << "Progenitor density of " << gcr[i].name << " for " << particle.name << ":\n";
	gcr[i].cr_density.print();
      }
      
    }  //  gcr.A>particle.A
    
    if (1 == galdef.K_capture) { //AWS20010828
      
      // ELECTRON ATTACHMENT or STRIPPING SOURCE   IMOS20010816
      
      double attach_H, strip_H, attach_He, strip_He;
      int zH = 1, zHe = 2;
      
      // Source from isotopes quickly decaying after EC 
      if (gcr[i].A == particle.A && gcr[i].Z != particle.Z)
	if (nucdata(galdef_network_par,
		    gcr[i].Z,
		    gcr[i].A,
		    1,
		    particle.Z,
		    particle.A, 
		    &Z1,
		    &A1,
		    &t_half)) {
	  
	  if (t_half > 0. && t_half/year2sec < galdef.t_half_limit) {
	    
	    for (ip = 0; ip<gcr[0].n_pgrid; ++ip) {
	      
	      Kcapture_cs(gcr[i].Ekin[ip],gcr[i].Z,zH, &attach_H ,&strip_H );
	      Kcapture_cs(gcr[i].Ekin[ip],gcr[i].Z,zHe,&attach_He,&strip_He);
	      cross_section[ip] = attach_H + galdef.He_H_ratio*attach_He;
	      
	    }
	    
	    if (galdef.verbose==-602 ||galdef.verbose>=1) { // selectable debug                  //AWS20010828
	      
	      ostringstream buf;
	      buf << "gen_secondary_source: electron attachment and fast decay cross_section for Z,A ="
		  <<gcr[i].Z<<","<<gcr[i].A<<" K_electron="<<gcr[i].K_electron
		  <<"->>"
		  <<particle.Z<<","<<particle.A<<"  K_electron="<<particle.K_electron
		  <<" t_half="<<t_half/year2sec<<" yr : ";
	      for(ip=0; ip<gcr[0].n_pgrid; buf << cross_section[ip++]<<" "); 
	      INFO(buf.str());
	      
	    }
	    
	    if (galdef.verbose==-601 ||galdef.verbose==10) { // selectable debug
	      
	      cout<<" progenitor density of "<<gcr[i].name<<" for "<<particle.name<<":\n";
	      gcr[i].cr_density.print();
              
	    }
	    
	  }
	  
	}
      
      // Source from electron attachment or stripping (for slow EC decays) 
      if (100*gcr[i].Z+gcr[i].A == 100*particle.Z+particle.A)
	if (gcr[i].K_electron != particle.K_electron) {
	  
	  for (ip = 0; ip < gcr[0].n_pgrid; ++ip) {
	    
	    Kcapture_cs(gcr[i].Ekin[ip],gcr[i].Z,zH, &attach_H ,&strip_H );
	    Kcapture_cs(gcr[i].Ekin[ip],gcr[i].Z,zHe,&attach_He,&strip_He);
	    
	    if(particle.K_electron) 
	      cross_section[ip] = attach_H + galdef.He_H_ratio*attach_He;
	    else                    
	      cross_section[ip] = strip_H + galdef.He_H_ratio*strip_He ;
	    
	  }
	  
	  if (galdef.verbose==-602||galdef.verbose>=1) { // selectable debug                  //AWS20010828
	    
	    ostringstream buf;
	    
	    if (particle.K_electron) 
	      buf<<"gen_secondary_source: electron attachment for Z,A =";
	    else
	      buf<<"gen_secondary_source: stripping cross_section for Z,A =";
	    
	    buf <<gcr[i].Z<<","<<gcr[i].A<<" K_electron="<<gcr[i].K_electron<<"->"
		<<particle.Z<<","<<particle.A<<"  K_electron="<<particle.K_electron<<" :";
	    for (ip=0; ip<gcr[0].n_pgrid; buf<<cross_section[ip++]<<" "); 
	    
	    INFO(buf.str());
	    
	  }
	  
	  if(galdef.verbose==-601 ||galdef.verbose==10) {// selectable debug   
	    
	    cout<<" progenitor density of "<<gcr[i].name<<" for "<<particle.name<<":\n";
	    gcr[i].cr_density.print();
	    
	  }
	  
	} // gcr[i].K_electron != particle.K_electron
	
    } //galdef.K_capture==1
      
    cross_section *= 1e-27; // mb -> cm^2
    
    if (2 == galdef.n_spatial_dimensions)
      for(int ir=0; ir<gcr[0].n_rgrid; ir++)
	for(int iz=0; iz<gcr[0].n_zgrid; iz++)
	  for(int ip=0; ip<gcr[0].n_pgrid; ip++)
	    gas_times_cross_section.d2[ir][iz].s[ip] = 
	      cross_section[ip]*gcr[i].beta[ip]*C*(galaxy.n_HI.d2[ir][iz].s[0] + 2.*galaxy.n_H2.d2[ir][iz].s[0] + galaxy.n_HII.d2[ir][iz].s[0]);
    
    if (3 == galdef.n_spatial_dimensions)
      for(int ix=0; ix<gcr[0].n_xgrid; ix++)
	for(int iy=0; iy<gcr[0].n_ygrid; iy++)
	  for(int iz=0; iz<gcr[0].n_zgrid; iz++)
	    for(int ip=0; ip<gcr[0].n_pgrid; ip++)
	      gas_times_cross_section.d3[ix][iy][iz].s[ip] = 
		cross_section[ip]*gcr[i].beta[ip]*C*(galaxy.n_HI.d3[ix][iy][iz].s[0] + 2.*galaxy.n_H2.d3[ix][iy][iz].s[0] + galaxy.n_HII.d3[ix][iy][iz].s[0]);
    
    gas_times_cross_section *= gcr[i].cr_density;
    particle.secondary_source_function += gas_times_cross_section;
    
    if (galdef.verbose==-603 ||galdef.verbose==10) {// selectable debug   
      
      cout<<" gas times cross section times progenitor density of "<<gcr[i].name<<" for "<<particle.name<<":\n";
      gas_times_cross_section.print();
      
    }
    
    // RADIOACTIVE DECAY SOURCE (beta+/-, EC) AWS20010806 IMOS20010816
    if (gcr[i].A >= particle.A) {
      
      if (100*gcr[i].Z+gcr[i].A != 100*particle.Z+particle.A) {
	
	const double branching_ratio = nucdata(galdef_network_par,gcr[i].Z,  
					       gcr[i].A, gcr[i].K_electron,
					       particle.Z,particle.A, 
					       &Z1, &A1, &t_half);
	
	if (branching_ratio != 0 && t_half != 0. && 100*Z1+A1 == 0) {
	  
	  if(galdef.verbose==-602||galdef.verbose>=1) { // selectable debug                  //AWS20010828
	    
	    ostringstream buf;
	    buf <<"gen_secondary_source: (Z,A,K_electron) "
		<< gcr[i].Z<<" "<<  gcr[i].A<<" "<<  gcr[i].K_electron<<"->"
		<<particle.Z<<" "<<particle.A
		<<" t_half="<<t_half/year2sec<<" yr  branching ratio="<<branching_ratio;
	    INFO(buf.str());
	    
	  }
	  
	  if (2 == galdef.n_spatial_dimensions)
	    for(int ir=0; ir<gcr[0].n_rgrid; ir++)
	      for(int iz=0; iz<gcr[0].n_zgrid; iz++)
		for(int ip=0; ip<gcr[0].n_pgrid; ip++)
		  decay_source.d2[ir] [iz].s[ip] = 
		    branching_ratio*gcr[i].cr_density.d2[ir][iz].s[ip]
		    /(gcr[i].gamma[ip]*t_half/log(2.0));
	  
	  if (3 == galdef.n_spatial_dimensions)
	    for(int ix=0; ix<gcr[0].n_xgrid; ix++)
	      for(int iy=0; iy<gcr[0].n_ygrid; iy++)
		for(int iz=0; iz<gcr[0].n_zgrid; iz++)
		  for(int ip=0; ip<gcr[0].n_pgrid; ip++)
		    decay_source.d3[ix][iy][iz].s[ip] = 
		      branching_ratio*gcr[i].cr_density.d3[ix][iy][iz].s[ip]
		      /(gcr[i].gamma[ip]*t_half/log(2.0));
	  
	  particle.secondary_source_function += decay_source;
          
	} // branching_ratio != 0 && t_half != 0. && 100*Z1+A1 == 0
        
      } // 100*gcr[i].Z+gcr[i].A != 100*particle.Z+particle.A
      
    } // gcr[i].A >= particle.A
    
  }  //  species
     
  if(galdef.verbose==-603 ||galdef.verbose==10) { // selectable debug
      
    cout<<" secondary source for "<<particle.name<<":\n";
    particle.secondary_source_function.print();
    
  }
  
  if(galdef.verbose==10) 
    particle.print();
  
  gas_times_cross_section.delete_array();
  decay_source.delete_array(); // AWS20010806
  
  INFO("Exit");
  
  return 0;
    
  }


//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
/*int Galprop::gen_TD_secondary_source(Particle& particle) { //YO/20140220

  INFO("Entry");

  //AWS20091009 fixed the changes introduced by TAP in r338: particle.name is a char array not a string

  //string particle_name;                                             //AWS20091009
  //particle_name = particle.name; // transfer from char to string    //AWS20091009

  const string particle_name = particle.name; 

  if (galdef.verbose >= 1) {

    ostringstream buf;
    buf << "Generating for (Z, A, K_electron) = (" 
	<< particle.Z << ", "
	<< particle.A << ", "
	<< particle.K_electron << ")";
    INFO(buf.str());

  }

  Distribution gas_times_cross_section;
  Distribution decay_source;                           //AWS20010806
  int i,ip, Z1,A1, galdef_network_par = 0;
  
  valarray<double> cross_section(0., gcr[0].n_pgrid);
  double t_half;

  //double* cross_section, t_half;
  //cross_section=new double[gcr[0].n_pgrid];

  if (2 == galdef.n_spatial_dimensions) {                   //AWS20010806
   
    gas_times_cross_section.init(gcr[0].n_rgrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
    decay_source.init(gcr[0].n_rgrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
  
  }
  
  if (3 == galdef.n_spatial_dimensions) {                   //AWS20010806
    
    gas_times_cross_section.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
    decay_source.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);          

  }

  particle.secondary_source_function = 0;

 // positron by time dependent injection YO/20140220
  if (particle_name == "TD_injection") {

    gen_TD_injection(particle);
    //INFO("Exit");
    //return 0;

  }

  INFO("Exit");
  
  //return 0;

  }
*/
