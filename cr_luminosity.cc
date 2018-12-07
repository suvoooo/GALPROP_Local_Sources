
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * cr_luminosity.cc *                      galprop package * 10/03/2003 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

using namespace std;//AWS20050624


#include"galprop_classes.h"
#include"galprop_internal.h"
#include"ErrorLogger.h"
#include <cstring>
#include <sstream>

// compute cosmic-ray luminosity of galaxy

/*
 

CR density  gcr.cr_density is in c/4pi * n(p) [ cm s^-1  * cm^-3 MeV^-1]
primary source function is    in c/4pi        [ cm s^-1  * cm^-3 MeV^-1 s-1]
luminosity (cm^-3)= Integral n(Ekin). A.Ekin. Ekin dlog(Ekin)
since each nucleus has KE=A.Ekin


The particle spectra are assumed to be on equal kinetic energy per nucleon grids
which is the standard for galprop.
    UNITS OF c/4pi*n(p) = (1/A)flux(Ekin)
 =                (1/A)c/4pi*beta*n(Ekin)

so factor required = A.4pi/c/beta 
*/

int Galprop::cr_luminosity()
{
   INFO("Entry");
   ostringstream buf;
   buf<<"computing cr_luminosity for n_spatial_dimensions="<<gcr[0].n_spatial_dimensions;
   INFO(buf.str());

   int stat=0;
   double CR_luminosity=0;
   double CR_particle_luminosity; // in particles instead of energy, e.g. for positrons  AWS20100309

// identify CR protons
   int iprotons=-1;
   for(int i=0;i<n_species;i++) if(gcr[i].A==1 && gcr[i].Z==1) iprotons=i;
   if(iprotons==-1) { WARNING("CR protons not found!"); return 1;}
   else {
      buf.str("");
      buf<<"  CR protons found as species #"<<iprotons;
      INFO(buf.str());
   }
 
// identify CR Helium
   int iHelium =-1;
   for(int i=0;i<n_species;i++) if(gcr[i].A==4 && gcr[i].Z==2) iHelium =i;
   if(iHelium ==-1) { WARNING("CR Helium  not found!"); return 1; }
   else {
      buf.str("");
      buf<<"  CR Helium  found as species #"<<iHelium;
      INFO(buf.str());
   }

   Particle particle;
   particle.init();  // to signal that arrays not yet allocated

   particle=gcr[0];
   particle.create_transport_arrays();

   for (int iprHe = 1; iprHe<=2; iprHe++)
   {
     if(iprHe==1) particle=gcr[iprotons];
     if(iprHe==2) particle=gcr[iHelium ];
     create_transport_arrays(particle);
     for(int ip=0;ip<gcr[iprotons].n_pgrid;ip++)
     {
        if(galaxy.n_spatial_dimensions==2)
           for(int ir=0;ir<gcr[iprotons].n_rgrid;ir++)
              for(int iz=0;iz<gcr[iprotons].n_zgrid;iz++)
                 CR_luminosity+= particle.primary_source_function.d2[ir][iz].s[ip] 
		   //            / particle.beta[ip]* pow(particle.Ekin[ip],2) *particle.A  *particle.A 
                   / particle.beta[ip]* pow(particle.Ekin[ip],2) *particle.A// number of nucleons gives total energy but Ekin/nucleon is the reference AWS20100322
                    *2.0*Pi*particle.r[ir]*galdef.dr*galdef.dz;
 
        if(galaxy.n_spatial_dimensions==3)
           for(int ix=0;ix<gcr[iprotons].n_xgrid;ix++)
              for(int iy=0;iy<gcr[iprotons].n_ygrid;iy++)
                 for(int iz=0;iz<gcr[iprotons].n_zgrid;iz++)
                 {
// weights for non-symmetric case
                    double sym_weight=1.0;
// weights for fully symmetric case
                    if(galdef.use_symmetry==1 && iz >0) sym_weight= 8.;
                    if(galdef.use_symmetry==1 && iz==0) sym_weight= 4.;// to avoid double-counting at z=0
                    CR_luminosity+= particle.primary_source_function.d3[ix][iy][iz].s[ip] 
		      / particle.beta[ip]* pow(particle.Ekin[ip],2) *particle.A  // AWS20100322
                       *galdef.dx*galdef.dy*galdef.dz *sym_weight;
                 }//iz//iy//ix//particle.n_spatial_dimensions==3
//cout<<"cr_luminosity:  ip= "<<ip<<" CR_luminosity="<<  CR_luminosity<<endl;

     }//ip
   }//iprHe

   CR_luminosity *= 4.0*Pi/c *1.0e6 *eV_to_erg *pow(kpc2cm,3) 
                    *log(gcr[iprotons].Ekin[1]/ gcr[iprotons].Ekin[0]) ;

 //CR_luminosity *= gcr[iprotons].normalization_factor;//IMOS20030217

   cout<<endl;
   cout<<"====================================================="<<endl;
   cout<<" CR p+He luminosity="<<  CR_luminosity <<" erg s^-1" <<endl;
     //<<  CR_luminosity*gcr[iprotons].normalization_factor <<endl;
   cout<<"====================================================="<<endl;
   cout<<endl;


   // protons and Helium separately AWS20100322

   for (int iprHe = 1; iprHe<=2; iprHe++)
   {
     if(iprHe==1) particle=gcr[iprotons];
     if(iprHe==2) particle=gcr[iHelium ];


     create_transport_arrays(particle);

     CR_luminosity=0.;
     for(int ip=0;ip<gcr[iprotons].n_pgrid;ip++)
     {
        if(galaxy.n_spatial_dimensions==2)
           for(int ir=0;ir<gcr[iprotons].n_rgrid;ir++)
              for(int iz=0;iz<gcr[iprotons].n_zgrid;iz++)
                 CR_luminosity+= particle.primary_source_function.d2[ir][iz].s[ip] 
		   / particle.beta[ip]* pow(particle.Ekin[ip],2) *particle.A  // number of nucleons gives total energy but Ekin/nucleon is the reference
                    *2.0*Pi*particle.r[ir]*galdef.dr*galdef.dz;
 
        if(galaxy.n_spatial_dimensions==3)
           for(int ix=0;ix<gcr[iprotons].n_xgrid;ix++)
              for(int iy=0;iy<gcr[iprotons].n_ygrid;iy++)
                 for(int iz=0;iz<gcr[iprotons].n_zgrid;iz++)
                 {
// weights for non-symmetric case
                    double sym_weight=1.0;
// weights for fully symmetric case
                    if(galdef.use_symmetry==1 && iz >0) sym_weight= 8.;
                    if(galdef.use_symmetry==1 && iz==0) sym_weight= 4.;// to avoid double-counting at z=0
                    CR_luminosity+= particle.primary_source_function.d3[ix][iy][iz].s[ip] 
                       / particle.beta[ip]* pow(particle.Ekin[ip],2) *particle.A // number of nucleons gives total energy but Ekin/nucleon is the reference
                       *galdef.dx*galdef.dy*galdef.dz *sym_weight;
                 }//iz//iy//ix//particle.n_spatial_dimensions==3
//cout<<"cr_luminosity:  ip= "<<ip<<" CR_luminosity="<<  CR_luminosity<<endl;

     }//ip
 

   CR_luminosity *= 4.0*Pi/c *1.0e6 *eV_to_erg *pow(kpc2cm,3) 
                    *log(gcr[iprotons].Ekin[1]/ gcr[iprotons].Ekin[0]) ;



   cout<<endl;
   cout<<"====================================================="<<endl;

   if(iprHe==1)cout<<" CR proton luminosity="<<  CR_luminosity <<" erg s^-1" <<endl;
   if(iprHe==2)cout<<" CR Helium luminosity="<<  CR_luminosity <<" erg s^-1" <<endl;
   cout<<"====================================================="<<endl;
   cout<<endl;

  }//iprHe

   //----------------------------------------------------------------


   if (galdef.output_gcr_full)
   for (int iprHe = 1; iprHe<=2; iprHe++) //AWS20100311
   {
     if(iprHe==1) particle=gcr[iprotons];
     if(iprHe==2) particle=gcr[iHelium ];

     create_transport_arrays(particle);

     store_gcr_source_functions(particle);
   }

 //------------------------------------------------------------------------------


// identify primary electrons IMOS20070725
  int ielectrons=-1;
   for(int i=0;i<n_species;i++) 
       if(!strcmp(gcr[i].name,"primary_electrons")) ielectrons=i;
   if(ielectrons==-1) { WARNING("CR primary electrons not found!"); return stat; }
   else {
      buf.str("");
      buf<<"  CR primary electrons found as species #"<<ielectrons;
      INFO(buf.str());
   }

   CR_luminosity=0;
   particle=gcr[ielectrons];
 
   create_transport_arrays(particle);
   for(int ip=0;ip<gcr[ielectrons].n_pgrid;ip++)
     {
       if(galaxy.n_spatial_dimensions==2)
	 for(int ir=0;ir<gcr[ielectrons].n_rgrid;ir++)
	   for(int iz=0;iz<gcr[ielectrons].n_zgrid;iz++)
	     CR_luminosity+= particle.primary_source_function.d2[ir][iz].s[ip] 
	       / particle.beta[ip]* pow(particle.Ekin[ip],2)
	       *2.0*Pi*particle.r[ir]*galdef.dr*galdef.dz;
       
       if(galaxy.n_spatial_dimensions==3)
	 for(int ix=0;ix<gcr[ielectrons].n_xgrid;ix++)
	   for(int iy=0;iy<gcr[ielectrons].n_ygrid;iy++)
	     for(int iz=0;iz<gcr[ielectrons].n_zgrid;iz++)
	       {
// weights for non-symmetric case
		 double sym_weight=1.0;
// weights for fully symmetric case
		 if(galdef.use_symmetry==1 && iz >0) sym_weight= 8.;
		 if(galdef.use_symmetry==1 && iz==0) sym_weight= 4.;// to avoid double-counting at z=0
		 CR_luminosity+= particle.primary_source_function.d3[ix][iy][iz].s[ip] 
		   / particle.beta[ip]* pow(particle.Ekin[ip],2)
		   *galdef.dx*galdef.dy*galdef.dz *sym_weight;
	       }//iz//iy//ix//particle.n_spatial_dimensions==3
//cout<<"cr_luminosity:  ip= "<<ip<<" CR_luminosity="<<  CR_luminosity<<endl;

     }//ip

   CR_luminosity *= 4.0*Pi/c *1.0e6 *eV_to_erg *pow(kpc2cm,3) 
                    *log(gcr[ielectrons].Ekin[1]/ gcr[ielectrons].Ekin[0]) ;

   cout<<endl;
   cout<<"====================================================="<<endl;
   cout<<" CR primary electron luminosity="<< CR_luminosity <<" erg s^-1" <<endl;
     //<< CR_luminosity*gcr[ielectrons].normalization_factor<<endl;
   cout<<"====================================================="<<endl;
   cout<<endl;


   if (galdef.output_gcr_full)store_gcr_source_functions(particle); //AWS20100311

   //--------------------------------------------------------------------------------------------------------------------

   
// identify secondary positrons IMOS20070725
   ielectrons=-1;
   for(int i=0;i<n_species;i++) 
       if(!strcmp(gcr[i].name,"secondary_positrons")) ielectrons=i;
   if(ielectrons==-1) { WARNING("CR secondary positrons not found!"); return stat; }
   else{
      buf.str("");
      buf<<"  CR secondary positrons found as species #"<<ielectrons;
      INFO(buf.str());
   }

   CR_luminosity=0;
   CR_particle_luminosity=0; //AWS20100309

   particle=gcr[ielectrons];
 
   create_transport_arrays(particle);
   gen_secondary_source(particle);

   for(int ip=0;ip<gcr[ielectrons].n_pgrid;ip++)
     {
       if(galaxy.n_spatial_dimensions==2)
	 for(int ir=0;ir<gcr[ielectrons].n_rgrid;ir++)
	   for(int iz=0;iz<gcr[ielectrons].n_zgrid;iz++)
	   {
             //AWS20100309:
	     CR_luminosity+= particle.secondary_source_function.d2[ir][iz].s[ip] 
	       / particle.beta[ip] *pow(particle.Ekin[ip],2)
	       *2.0*Pi*particle.r[ir]*galdef.dr*galdef.dz;

             //AWS20100309:
	     CR_particle_luminosity+= particle.secondary_source_function.d2[ir][iz].s[ip] 
	       / particle.beta[ip] *pow(particle.Ekin[ip],1)         
	       *2.0*Pi*particle.r[ir]*galdef.dr*galdef.dz;
	   } 
       
       if(galaxy.n_spatial_dimensions==3)
	 for(int ix=0;ix<gcr[ielectrons].n_xgrid;ix++)
	   for(int iy=0;iy<gcr[ielectrons].n_ygrid;iy++)
	     for(int iz=0;iz<gcr[ielectrons].n_zgrid;iz++)
	       {
// weights for non-symmetric case
		 double sym_weight=1.0;
// weights for fully symmetric case
		 if(galdef.use_symmetry==1 && iz >0) sym_weight= 8.;
		 if(galdef.use_symmetry==1 && iz==0) sym_weight= 4.;// to avoid double-counting at z=0

		 //AWS20100309:
		 CR_luminosity+= particle.secondary_source_function.d3[ix][iy][iz].s[ip] 
		   / particle.beta[ip] *pow(particle.Ekin[ip],2)
		   *galdef.dx*galdef.dy*galdef.dz *sym_weight;

		 //AWS20100309:
		 CR_particle_luminosity+= particle.secondary_source_function.d3[ix][iy][iz].s[ip] 
		   / particle.beta[ip] *pow(particle.Ekin[ip],1)    
		   *galdef.dx*galdef.dy*galdef.dz *sym_weight;


	       }//iz//iy//ix//particle.n_spatial_dimensions==3


//cout<<"cr_luminosity:  ip= "<<ip<<" CR_luminosity="<<  CR_luminosity<<endl;

     }//ip


   CR_luminosity *= 4.0*Pi/c *pow(kpc2cm,3)           *1.0e6*eV_to_erg                 //AWS20100309
                    *log(gcr[ielectrons].Ekin[1]/ gcr[ielectrons].Ekin[0]) ;

   CR_particle_luminosity *= 4.0*Pi/c *pow(kpc2cm,3)                                   //AWS20100309
                    *log(gcr[ielectrons].Ekin[1]/ gcr[ielectrons].Ekin[0]) ;

   cout<<endl;
   cout<<"====================================================="<<endl;
   cout<<" CR secondary positron luminosity="<< CR_luminosity                <<" erg s^-1" <<endl;
   cout<<" CR secondary positron luminosity="<< CR_particle_luminosity <<" positrons s^-1" <<endl;  
   cout<<"====================================================="<<endl;
   cout<<endl;
   
   if (galdef.output_gcr_full)store_gcr_source_functions(particle); //AWS20100311

   //-----------------------------------------------------------------------------------


// identify secondary electrons AWS20100310
   ielectrons=-1;
   for(int i=0;i<n_species;i++) 
       if(!strcmp(gcr[i].name,"secondary_electrons")) ielectrons=i;
   if(ielectrons==-1) { WARNING("CR secondary electrons not found!"); return stat; }
   else {
      buf.str("");
      buf<<"  CR secondary electrons found as species #"<<ielectrons;
      INFO(buf.str());
   }

   CR_luminosity=0;
   CR_particle_luminosity=0; //AWS20100309

   particle=gcr[ielectrons];
 
   create_transport_arrays(particle);
   gen_secondary_source(particle);

   for(int ip=0;ip<gcr[ielectrons].n_pgrid;ip++)
     {
       if(galaxy.n_spatial_dimensions==2)
	 for(int ir=0;ir<gcr[ielectrons].n_rgrid;ir++)
	   for(int iz=0;iz<gcr[ielectrons].n_zgrid;iz++)
	   {
             //AWS20100309:
	     CR_luminosity+= particle.secondary_source_function.d2[ir][iz].s[ip] 
	       / particle.beta[ip] *pow(particle.Ekin[ip],2)
	       *2.0*Pi*particle.r[ir]*galdef.dr*galdef.dz;

             //AWS20100309:
	     CR_particle_luminosity+= particle.secondary_source_function.d2[ir][iz].s[ip] 
	       / particle.beta[ip] *pow(particle.Ekin[ip],1)         
	       *2.0*Pi*particle.r[ir]*galdef.dr*galdef.dz;
	   } 
       
       if(galaxy.n_spatial_dimensions==3)
	 for(int ix=0;ix<gcr[ielectrons].n_xgrid;ix++)
	   for(int iy=0;iy<gcr[ielectrons].n_ygrid;iy++)
	     for(int iz=0;iz<gcr[ielectrons].n_zgrid;iz++)
	       {
// weights for non-symmetric case
		 double sym_weight=1.0;
// weights for fully symmetric case
		 if(galdef.use_symmetry==1 && iz >0) sym_weight= 8.;
		 if(galdef.use_symmetry==1 && iz==0) sym_weight= 4.;// to avoid double-counting at z=0

		 //AWS20100309:
		 CR_luminosity+= particle.secondary_source_function.d3[ix][iy][iz].s[ip] 
		   / particle.beta[ip] *pow(particle.Ekin[ip],2)
		   *galdef.dx*galdef.dy*galdef.dz *sym_weight;

		 //AWS20100309:
		 CR_particle_luminosity+= particle.secondary_source_function.d3[ix][iy][iz].s[ip] 
		   / particle.beta[ip] *pow(particle.Ekin[ip],1)    
		   *galdef.dx*galdef.dy*galdef.dz *sym_weight;


	       }//iz//iy//ix//particle.n_spatial_dimensions==3


//cout<<"cr_luminosity:  ip= "<<ip<<" CR_luminosity="<<  CR_luminosity<<endl;

     }//ip


   CR_luminosity *= 4.0*Pi/c *pow(kpc2cm,3)           *1.0e6*eV_to_erg                 //AWS20100309
                    *log(gcr[ielectrons].Ekin[1]/ gcr[ielectrons].Ekin[0]) ;

   CR_particle_luminosity *= 4.0*Pi/c *pow(kpc2cm,3)                                   //AWS20100309
                    *log(gcr[ielectrons].Ekin[1]/ gcr[ielectrons].Ekin[0]) ;

   cout<<endl;
   cout<<"====================================================="<<endl;
   cout<<" CR secondary electron luminosity="<< CR_luminosity                <<" erg s^-1" <<endl;
   cout<<" CR secondary electron luminosity="<< CR_particle_luminosity <<" electrons s^-1" <<endl;  
   cout<<"====================================================="<<endl;
   cout<<endl;
   
   if (galdef.output_gcr_full)store_gcr_source_functions(particle); //AWS20100311


   INFO("Exit");
	 particle.delete_arrays();
   return stat;
}
