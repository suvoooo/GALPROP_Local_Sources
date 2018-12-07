
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * create_gcr.cc *                               galprop package * 08/16/2001 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"
#include <cstring>

#include <ErrorLogger.h>

int Galprop::create_gcr()//AWS20050816
{
  //cout<<" >>>> create_gcr"<<endl;
  INFO("Entry"); 
  if(galdef.verbose >= 1)galdef.print();
  
   int i=0,j,Z,A,Z1,A1, stat=0;
   char name[100];
   char *element[]=
      {
        "Hydrogen",    "Helium",   "Lithium", "Beryllium",     "Boron",
          "Carbon",  "Nitrogen",    "Oxygen",  "Fluorine",      "Neon",
	  "Sodium", "Magnesium",  "Aluminium",   "Silicon","Phosphorus",
          "Sulphur",  "Chlorine",     "Argon", "Potassium",   "Calcium",
        "Scandium",  "Titanium",  "Vanadium",  "Chromium", "Manganese",
            "Iron",    "Cobalt",    "Nickel",    "Copper",      "Zinc"
      };
   double t_half[2],branching_ratio;                                   // IMOS20010816
   int galdef_network_par=0;         // imos network, use everywhere
   int K_electron;                                                     // AWS20010731

// calculate the number of species
   if(!galdef.secondary_antiprotons) galdef.tertiary_antiprotons=0;    // IMOS20000802

   n_species=0;
   if(galdef.secondary_positrons)   n_species++;
   if(galdef.knock_on_electrons)    n_species++; //IMOS20060504
   if(galdef.secondary_electrons)   n_species++;
   if(galdef.primary_electrons)     n_species++;

//DM annihilation positron electron and antiprotons by BiXJ
   if(galdef.primary_DM_positron)   n_species++; // BiXJ 2007/2/1
   if(galdef.primary_DM_electron)   n_species++; // BiXJ 2007/2/1
   if(galdef.primary_DM_antip)   n_species++; // BiXJ 2007/2/1

//time dependent
   if(galdef.primary_TD_positron) n_species++; //YO/20140220 -> YO/20140408 changed
   if(galdef.primary_TD_electron) n_species++; //YO/20140408

// DM decay species: IMOS20050912
   if(galdef.DM_positrons)          n_species++; 
   if(galdef.DM_electrons)          n_species++;
   if(galdef.DM_antiprotons)        n_species++;

   if(galdef.tertiary_antiprotons)  n_species++;                       // IMOS20000605
   if(galdef.secondary_antiprotons) n_species++;
   if(galdef.secondary_protons)     n_species++;                       // IMOS20000605
   for(Z=1; Z<=galdef.max_Z; Z++)
   {
      if(!galdef.use_Z[Z]) continue;
      for(A=2*Z-2; A<2.5*Z+4.2; A++)                                   // IMOS20010816 whole loop
      {
 	 if(!nucdata(galdef_network_par,Z,A,0,Z,A,&Z1,&A1,&t_half[0])) continue;
 	 if(!nucdata(galdef_network_par,Z,A,1,Z,A,&Z1,&A1,&t_half[1])) continue;
	 for(K_electron=0;K_electron<=galdef.K_capture;K_electron++)
	 {
            if(t_half[K_electron]>0. && t_half[K_electron]/year2sec<galdef.t_half_limit) continue;
            if(K_electron) if(t_half[0] == t_half[1]) continue;
            n_species++;
         }
      }
   }

   if(galdef.verbose >= 1)cout<<endl<<"Number of species to create: n_species= "<<n_species<<endl<<endl;
   if(!n_species) {cout<<"create_gcr.cc: No particles specified -exit"<<endl; exit(1);}

// create a Particle array
   delete[] gcr;
   gcr=new Particle[n_species];
   K_electron=0; // for all secondaries and non-nuclei                    AWS20010731

// DM POSITRONS   IMOS20050912

   if(galdef.DM_positrons)
   {
      strcpy(name,"DM_positrons");
      Z=1;  A=0;  t_half[0]=0.0;                                                  // IMOS20010816

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                               // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                               // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print(); 
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// DM ELECTRONS   IMOS20050912

   if(galdef.DM_electrons)
   {
      strcpy(name,"DM_electrons");
      Z=-1;  A=0;  t_half[0]=0.0;                                                  // IMOS20010816

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                               // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                               // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print(); 
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// DM annihilated positron by BiXJ 2007/2/1

   if(galdef.primary_DM_positron)
   {
      strcpy(name,"primary_DM_positron");
      Z=1;  A=0;  t_half[0]=0.0;                                                 

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                            
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                             
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      gcr[i].print(); 
      cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// DM annihilated electron by BiXJ 2007/2/1

   if(galdef.primary_DM_electron)
   {
      strcpy(name,"primary_DM_electron");
      Z=-1;  A=0;  t_half[0]=0.0;                                                  

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                               
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                              
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      gcr[i].print(); 
      cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// DM annihilated antiproton  by BiXJ 2007/2/1

   if(galdef.primary_DM_antip)
   {
      strcpy(name,"primary_DM_antiproton");
      Z=-1;  A=1;  t_half[0]=0.0;                                         

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                       
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                       
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      gcr[i].print();
      cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

//time dependent positron YO/20140220 -> 20140408 changed

   if(galdef.primary_TD_positron)
   {
      strcpy(name,"primary_TD_positron");
      Z=1;  A=0;  t_half[0]=0.0;                                                 

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                            
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                             
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      gcr[i].print(); 
      cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
      }

//time dependent electron YO/20140408

   if(galdef.primary_TD_electron)
   {
      strcpy(name,"primary_TD_electron");
      Z=-1;  A=0;  t_half[0]=0.0;                                                 

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                            
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                             
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      gcr[i].print(); 
      cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
      }

// SECONDARY POSITRONS

   if(galdef.secondary_positrons)
   {
      strcpy(name,"secondary_positrons");
      Z=1;  A=0;  t_half[0]=0.0;                                                  // IMOS20010816

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                               // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                               // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print(); 
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// KNOCK-ON ELECTRONS

   if(galdef.knock_on_electrons)                                         // IMOS20060504
   {
      strcpy(name,"knock_on_electrons");
      Z=-1;  A=0;  t_half[0]=0.0;                                         // IMOS20010816

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                       // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                       // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print();
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// SECONDARY ELECTRONS

   if(galdef.secondary_electrons)
   {
      strcpy(name,"secondary_electrons");
      Z=-1;  A=0;  t_half[0]=0.0;                                         // IMOS20010816

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                       // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                       // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print();
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// PRIMARY ELECTRONS

   if(galdef.primary_electrons)
   {
      strcpy(name,"primary_electrons");
      Z=-1;  A=0;  t_half[0]=0.0;                                         // IMOS20010816

      gcr[i].primary_abundance=1.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                       // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                       // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx, 
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print(); 
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// DM ANTIPROTONS   IMOS20050912

   if(galdef.DM_antiprotons)
   {
      strcpy(name,"DM_antiprotons");
      Z=-1;  A=1;  t_half[0]=0.0;                                         // IMOS20010816

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                       // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                       // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print();
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// TERTIARY ANTIPROTONS

   if(galdef.tertiary_antiprotons)
   {
      strcpy(name,"tertiary_antiprotons");
      Z=-1;  A=1;  t_half[0]=0.0;                                         // IMOS20010816

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                       // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                       // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print();
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// SECONDARY ANTIPROTONS

   if(galdef.secondary_antiprotons)
   {
      strcpy(name,"secondary_antiprotons");
      Z=-1;  A=1;  t_half[0]=0.0;                                         // IMOS20010816

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                       // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                       // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print();
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// NUCLEONS: SECONDARY PROTONS

   if(galdef.secondary_protons)
   {
      strcpy(name,"secondary_protons");
      Z=1;  A=1;  t_half[0]=0.0;                                          // IMOS20010816

      gcr[i].primary_abundance=0.0;

      if(galdef.n_spatial_dimensions==2)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                       // IMOS20010816
            galdef.r_min,  galdef.r_max, galdef.dr,  
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid);  

      if(galdef.n_spatial_dimensions==3)
         gcr[i].init(name,Z,A,K_electron,t_half[0],                       // IMOS20010816
            galdef.x_min,  galdef.x_max, galdef.dx,  
            galdef.y_min,  galdef.y_max, galdef.dy,
            galdef.z_min,  galdef.z_max, galdef.dz,
            galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
            galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
            galdef.p_Ekin_grid); 

      if(galdef.verbose >= 1)gcr[i].print();
      if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
      i++;
   }

// OTHER NUCLEONS

   strcpy(name,"particle");

   for(Z=1; Z<=galdef.max_Z; Z++)
   {
      if(!galdef.use_Z[Z]) continue;
      for(A=2*Z-2; A<2.5*Z+4.2; A++)
      {
 	 if(!nucdata(galdef_network_par,Z,A,0,Z,A,&Z1,&A1,&t_half[0])) continue;          // IMOS20010816
 	 if(!nucdata(galdef_network_par,Z,A,1,Z,A,&Z1,&A1,&t_half[1])) continue;          // IMOS20010816
         for(K_electron=galdef.K_capture;K_electron>=0;K_electron--)                      // AWS20010731
	 {
            if(galdef.verbose >= 1)cout<<"nucleus being tested Z A K_electron "<<Z<<" "<<A<<" "<<K_electron<<endl;//AWS20010731
	                                                                                  // IMOS20010816 next 2 lines
            if(t_half[K_electron]>0. && t_half[K_electron]/year2sec<galdef.t_half_limit) continue;
            if(K_electron) if(t_half[0] == t_half[1]) continue;

            t_half[K_electron]/= year2sec;
            if(galdef.verbose >= 1)cout<<"nucleus used Z A t_half "<<Z<<" "<<A<<" "<<t_half[K_electron]<<endl;
            for(j=1; j<31; j++) if(Z == j) sprintf(name,  "%s_%d",element[j-1],A);

            gcr[i].primary_abundance=galdef.isotopic_abundance[Z][A];
            if(K_electron >0) gcr[i].primary_abundance= 0.0;   // no prim. K-capture nuclei AWS20010731
     
            if(K_electron==1) strcat(name,"*");                                         // AWS20010731
            strcat(name,"\0");                                                          // IMOS20010816

            if(galdef.n_spatial_dimensions==2)
               gcr[i].init(name,Z,A,K_electron,t_half[K_electron],                      // IMOS20010816
                  galdef.r_min,  galdef.r_max, galdef.dr,  
                  galdef.z_min,  galdef.z_max, galdef.dz,
                  galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
                  galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
                  galdef.p_Ekin_grid);  

            if(galdef.n_spatial_dimensions==3)
               gcr[i].init(name,Z,A,K_electron,t_half[K_electron],                      // IMOS20010816
                  galdef.x_min,  galdef.x_max, galdef.dx,  
                  galdef.y_min,  galdef.y_max, galdef.dy,
                  galdef.z_min,  galdef.z_max, galdef.dz,
                  galdef.   p_min,  galdef.   p_max, galdef.   p_factor,
                  galdef.Ekin_min,  galdef.Ekin_max, galdef.Ekin_factor,
                  galdef.p_Ekin_grid); 

            if(galdef.verbose >= 1)gcr[i].print();
            if(galdef.verbose >= 1)cout<<"============== completed creation of "<<gcr[i].name<<endl<<endl;
            i++;
	 } // K_electron
      } //A
   } //Z

   INFO("Exit");
   //cout<<" <<<< create_gcr"<<endl;
   return stat;
}




