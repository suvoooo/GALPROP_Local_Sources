
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Galaxy.cc *                                   galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include<iostream>
#include"Galaxy.h"
#include"ErrorLogger.h"

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

//using namespace rf;

Galaxy::Galaxy() {

   init_pointers();

  //fISRFModel = 0;

}

Galaxy::Galaxy(const Galaxy & old) {
   //Copy all the variables
   copy_variables(old);
   //Do a deep copy of all the pointers that are initialized in the old version
   copy_pointers_deep(old);
}

void Galaxy::init_pointers() {

  x = y = z = r = E_gamma = nu_synch = 0 ;
  
  R_bins = 0;

  ISRF = ISRF_energy_density = IC_iso_emiss = IC_aniso_emiss = IC_iso_skymap = IC_aniso_skymap = 0;

  IC_iso_hp_skymap = IC_aniso_hp_skymap = bremss_H2R_hp_skymap = bremss_HIR_hp_skymap = bremss_HII_hp_skymap = 0;

  pi0_decay_H2R_hp_skymap = pi0_decay_HIR_hp_skymap = pi0_decay_HII_hp_skymap = 0;
  
}

void Galaxy::copy_variables(const Galaxy & old)
{
   z_min = old.z_min;
   z_max = old.z_max;
   dz = old.dz;
   r_min = old.r_min;
   r_max = old.r_max;
   dr = old.dr;
   x_min = old.x_min;
   x_max = old.x_max;
   dx = old.dx;
   y_min = old.y_min;
   y_max = old.y_max;
   dy = old.dy;
   p_min = old.p_min;
   p_max = old.p_max;
   p_factor = old.p_factor;
   n_spatial_dimensions = old.n_spatial_dimensions;
   n_pgrid = old.n_pgrid;
   n_xgrid = old.n_xgrid;
   n_rgrid = old.n_rgrid;
   n_zgrid = old.n_zgrid;
   n_ygrid = old.n_ygrid;
   n_E_gammagrid = old.n_E_gammagrid;
   E_gamma_min = old.E_gamma_min;
   E_gamma_max = old.E_gamma_max;
   E_gamma_factor = old.E_gamma_factor;
   n_lat = old.n_lat;
   n_long = old.n_long;
   long_min = old.long_min;
   long_max = old.long_max;
   lat_min = old.lat_min;
   lat_max = old.lat_max;
   d_long = old.d_long;
   d_lat = old.d_lat;
   n_nu_synchgrid = old.n_nu_synchgrid;
   nu_synch_min = old.nu_synch_min;
   nu_synch_max = old.nu_synch_max;
   nu_synch_factor = old.nu_synch_factor;
   n_H2 = old.n_H2;
   n_HI = old.n_HI;
   n_HII = old.n_HII;
   HIR = old.HIR;
   COR = old.COR;
   n_Ring = old.n_Ring;
   hpCOR = old.hpCOR;
   hpHIR = old.hpHIR;
   B_field = old.B_field;
   n_ISRF_components = old.n_ISRF_components;
   fISRFFactors.resize(old.fISRFFactors.size());
   fISRFFactors = old.fISRFFactors;
   nu_ISRF.resize(old.nu_ISRF.size());
   nu_ISRF = old.nu_ISRF;
   SNR_cell_time = old.SNR_cell_time;
   SNR_cell_phase = old.SNR_cell_phase;
   SNR_electron_dg = old.SNR_electron_dg;
   SNR_nuc_dg = old.SNR_nuc_dg;
   bremss_emiss = old.bremss_emiss;
   bremss_ionized_emiss = old.bremss_ionized_emiss;
   pi0_decay_emiss = old.pi0_decay_emiss;
   bremss_skymap = old.bremss_skymap;
   bremss_ionized_skymap = old.bremss_ionized_skymap;
   pi0_decay_skymap = old.pi0_decay_skymap;
   bremss_H2R_skymap = old.bremss_H2R_skymap;
   bremss_HIR_skymap = old.bremss_HIR_skymap;
   bremss_HII_skymap = old.bremss_HII_skymap;
   pi0_decay_H2R_skymap = old.pi0_decay_H2R_skymap;
   pi0_decay_HIR_skymap = old.pi0_decay_HIR_skymap;
   pi0_decay_HII_skymap = old.pi0_decay_HII_skymap;
   bremss_hp_skymap = old.bremss_hp_skymap;
   bremss_ionized_hp_skymap = old.bremss_ionized_hp_skymap;
   pi0_decay_hp_skymap = old.pi0_decay_hp_skymap;
   synchrotron_emiss = old.synchrotron_emiss;
   synchrotron_skymap = old.synchrotron_skymap;
   ionization_rate = old.ionization_rate;
   DM_emiss = old.DM_emiss;
   DM_skymap = old.DM_skymap;
   DM_hp_skymap = old.DM_hp_skymap;
}

void Galaxy::copy_pointers_deep(const Galaxy &old) 
{
   //Does not free allocated memory so this can be used in copy constructor
   init_pointers(); //Don't want un-initialized pointers.
   if (old.x) {
      x = new double[old.n_xgrid];
      for (int i = 0; i < old.n_xgrid; ++i) {
	 x[i] = old.x[i];
      }
   }
   if (old.r) {
      r = new double[old.n_rgrid];
      for (int i = 0; i < old.n_rgrid; ++i) {
	 r[i] = old.r[i];
      }
   }
   if (old.y) {
      y = new double[old.n_ygrid];
      for (int i = 0; i < old.n_ygrid; ++i) {
	 y[i] = old.y[i];
      }
   }
   if (old.z) {
      z = new double[old.n_zgrid];
      for (int i = 0; i < old.n_zgrid; ++i) {
	 z[i] = old.z[i];
      }
   }
   if (old.E_gamma) {
      E_gamma = new double[old.n_E_gammagrid];
      for (int i = 0; i < old.n_E_gammagrid; ++i) {
	 E_gamma[i] = old.E_gamma[i];
      }
   }
   if (old.nu_synch) {
      nu_synch = new double[old.n_nu_synchgrid];
      for (int i = 0; i < old.n_nu_synchgrid; ++i) {
	 nu_synch[i] = old.nu_synch[i];
      }
   }
   if (old.R_bins) {
      R_bins = new float[old.n_Ring];
      for (int i = 0; i < old.n_Ring; ++i) {
	 R_bins[i] = old.R_bins[i];
      }
   }
   if (old.ISRF) {
      ISRF = new Distribution[old.n_ISRF_components];
      for (int i = 0; i < old.n_ISRF_components; ++i) {
	 ISRF[i] = old.ISRF[i];
      }
   }
   if (old.ISRF_energy_density) {
      ISRF_energy_density = new Distribution[old.n_ISRF_components];
      for (int i = 0; i < old.n_ISRF_components; ++i) {
	 ISRF_energy_density[i] = old.ISRF_energy_density[i];
      }
   }
   if (old.IC_iso_emiss) {
      IC_iso_emiss = new Distribution[old.n_ISRF_components];
      for (int i = 0; i < old.n_ISRF_components; ++i) {
	 IC_iso_emiss[i] = old.IC_iso_emiss[i];
      }
   }
   if (old.IC_aniso_emiss) {
      IC_aniso_emiss = new Distribution;
      *IC_aniso_emiss = *(old.IC_aniso_emiss);
   }
   if (old.IC_iso_skymap) {
      IC_iso_skymap = new Distribution[old.n_ISRF_components];
      for (int i = 0; i < old.n_ISRF_components; ++i) {
	 IC_iso_skymap[i] = old.IC_iso_skymap[i];
      }
   }
   if (old.IC_aniso_skymap) {
      IC_aniso_skymap = new Distribution[old.n_ISRF_components];
      for (int i = 0; i < old.n_ISRF_components; ++i) {
	 IC_aniso_skymap[i] = old.IC_aniso_skymap[i];
      }
   }
   if (old.IC_iso_hp_skymap) {
      IC_iso_hp_skymap = new Skymap<double>[old.n_ISRF_components];
      for (int i = 0; i < old.n_ISRF_components; ++i) {
	 IC_iso_hp_skymap[i] = old.IC_iso_hp_skymap[i];
      }
   }
   if (old.IC_aniso_hp_skymap) {
      IC_aniso_hp_skymap = new Skymap<double>[old.n_ISRF_components];
      for (int i = 0; i < old.n_ISRF_components; ++i) {
	 IC_aniso_hp_skymap[i] = old.IC_aniso_hp_skymap[i];
      }
   }
   if (old.bremss_H2R_hp_skymap){
      bremss_H2R_hp_skymap = new Skymap<double>[old.n_Ring];
      for (int i = 0; i < old.n_Ring; ++i) {
	 bremss_H2R_hp_skymap[i] = old.bremss_H2R_hp_skymap[i];
      }
   }
   if (old.bremss_HIR_hp_skymap){
      bremss_HIR_hp_skymap = new Skymap<double>[old.n_Ring];
      for (int i = 0; i < old.n_Ring; ++i) {
	 bremss_HIR_hp_skymap[i] = old.bremss_HIR_hp_skymap[i];
      }
   }
   if (old.bremss_HII_hp_skymap){
      bremss_HII_hp_skymap = new Skymap<double>[old.n_Ring];
      for (int i = 0; i < old.n_Ring; ++i) {
	 bremss_HII_hp_skymap[i] = old.bremss_HII_hp_skymap[i];
      }
   }
   if (old.pi0_decay_H2R_hp_skymap){
      pi0_decay_H2R_hp_skymap = new Skymap<double>[old.n_Ring];
      for (int i = 0; i < old.n_Ring; ++i) {
	 pi0_decay_H2R_hp_skymap[i] = old.pi0_decay_H2R_hp_skymap[i];
      }
   }
   if (old.pi0_decay_HIR_hp_skymap){
      pi0_decay_HIR_hp_skymap = new Skymap<double>[old.n_Ring];
      for (int i = 0; i < old.n_Ring; ++i) {
	 pi0_decay_HIR_hp_skymap[i] = old.pi0_decay_HIR_hp_skymap[i];
      }
   }
   if (old.pi0_decay_HII_hp_skymap){
      pi0_decay_HII_hp_skymap = new Skymap<double>[old.n_Ring];
      for (int i = 0; i < old.n_Ring; ++i) {
	 pi0_decay_HII_hp_skymap[i] = old.pi0_decay_HII_hp_skymap[i];
      }
   }

}

Galaxy::~Galaxy() {

   free_memory();
}

void Galaxy::free_memory() {
  delete[] x;
  delete[] y;
  delete[] z;
  delete[] r;

  delete[] E_gamma;
  delete[] nu_synch;

  delete[] R_bins;

  delete[] ISRF;
  delete[] ISRF_energy_density;
  
  delete[] IC_iso_emiss;
  delete IC_aniso_emiss;
  delete[] IC_iso_skymap;
  delete[] IC_aniso_skymap;
  delete[] IC_iso_hp_skymap;          // inverse Compton   isotropic intensity skymap for each ISRF component
  delete[] IC_aniso_hp_skymap;        // inverse Compton anisotropic intensity skymap for each ISRF component IMOS20060420
  delete[] bremss_H2R_hp_skymap;   // bremsstrahlung intensity skymap on CO  in rings                 AWS20041214
  delete[] bremss_HIR_hp_skymap;   // bremsstrahlung intensity skymap on HI  in rings                 AWS20041214
  delete[] bremss_HII_hp_skymap;   // bremsstrahlung intensity skymap on HII  in rings                IMOS20080114*
  delete[] pi0_decay_H2R_hp_skymap;   // pi0 decay      intensity skymap on H2  in rings                 AWS20041214
  delete[] pi0_decay_HIR_hp_skymap;   // pi0 decay      intensity skymap on HI  in rings                 AWS20041214
  delete[] pi0_decay_HII_hp_skymap;   // pi0 decay      intensity skymap on HII in rings                 IMOS20080114*
  
}

void Galaxy::init(double r_min_, double r_max_, double dr_, 
                  double z_min_, double z_max_, double dz_) {

  INFO("Initializing 2D");
  n_spatial_dimensions=2;       //2D    
  
  r_min=r_min_;
  r_max=r_max_;
  dr   =   dr_;
  z_min=z_min_;
  z_max=z_max_;
  dz   =   dz_;
  
  n_rgrid=(int)((r_max-r_min)/dr + 1.5);
  n_zgrid=(int)((z_max-z_min)/dz + 1.5);
  
  delete[] r;
  r=new double[n_rgrid];
  delete[] z;
  z=new double[n_zgrid];
  
  int ir,iz;    
  for(ir=0; ir<n_rgrid; ir++) r[ir]=r_min+ir*dr;
  for(iz=0; iz<n_zgrid; iz++) z[iz]=z_min+iz*dz;
  
  n_HI.init(n_rgrid, n_zgrid, 1);
  n_H2.init(n_rgrid, n_zgrid, 1);
  n_HII.init(n_rgrid, n_zgrid, 1);
  B_field.init(n_rgrid, n_zgrid, 1);
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void Galaxy::init( double x_min_, double x_max_, double dx_, 
                   double y_min_, double y_max_, double dy_, 
                   double z_min_, double z_max_, double dz_) {  

  INFO("Initializing 3D");
  n_spatial_dimensions=3;       //3D    
  
  x_min=x_min_;
  x_max=x_max_;
  dx   =   dx_;
  y_min=y_min_;
  y_max=y_max_;
  dy   =   dy_;
  z_min=z_min_;
  z_max=z_max_;
  dz   =   dz_;
  
  n_xgrid=(int)((x_max-x_min)/dx + 1.5);
  n_ygrid=(int)((y_max-y_min)/dy + 1.5);
  n_zgrid=(int)((z_max-z_min)/dz + 1.5);
  
  delete[] x;
  x=new double[n_xgrid];
  delete[] y;
  y=new double[n_ygrid];
  delete[] z;
  z=new double[n_zgrid];
  
  int ix,iy,iz;    
  for(ix=0; ix<n_xgrid; ix++) x[ix]=x_min+ix*dx; 
  for(iy=0; iy<n_ygrid; iy++) y[iy]=y_min+iy*dy; 
  for(iz=0; iz<n_zgrid; iz++) z[iz]=z_min+iz*dz; 
  
  n_HI.init(n_xgrid, n_ygrid, n_zgrid, 1);
  n_H2.init(n_xgrid, n_ygrid, n_zgrid, 1);
  n_HII.init(n_xgrid, n_ygrid, n_zgrid, 1);	
  B_field.init(n_xgrid, n_ygrid, n_zgrid, 1);
  
  /* Moved to create_galaxy          Gulli20070810
     SNR_cell_time  .init(n_xgrid, n_ygrid, n_zgrid,1);
     SNR_cell_phase .init(n_xgrid, n_ygrid, n_zgrid,1);
     
     SNR_electron_dg.init(n_xgrid, n_ygrid, n_zgrid,1); //AWS20010410
     SNR_nuc_dg     .init(n_xgrid, n_ygrid, n_zgrid,1); //AWS20010410
  */
  
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void Galaxy::print() {

  cout<<"Galaxy: n_spatial_dimensions="<<n_spatial_dimensions<<endl;
  cout<<"        n_ISRF_components   ="<<n_ISRF_components   <<endl;
  
  if(n_spatial_dimensions==2) {

    cout<<"r_min="<<r_min<<endl;
    cout<<"r_max="<<r_max<<endl;
    cout<<"dr   ="<<dr   <<endl;
    cout<<"z_min="<<z_min<<endl;
    cout<<"z_max="<<z_max<<endl;
    cout<<"dz   ="<<dz   <<endl;
    
    cout<<"n_rgrid="<<n_rgrid<<endl;
    cout<<"n_zgrid="<<n_zgrid<<endl;
    
    cout<<"r grid:"<<endl;
    for(int ir=0; ir<n_rgrid; cout<<r[ir++]<<" "); cout<<endl;
    cout<<"z grid:"<<endl;
    for(int iz=0; iz<n_zgrid; cout<<z[iz++]<<" "); cout<<endl;
    
    cout<<"n_HI[0][0].s[0]="<<n_HI.d2[0][0].s[0]<<endl;
  }//(n_spatial_dimensions==2

  if(n_spatial_dimensions==3) {

    cout<<"x_min="<<x_min<<endl;
    cout<<"x_max="<<x_max<<endl;
    cout<<"dx   ="<<dx   <<endl;
    cout<<"y_min="<<y_min<<endl;
    cout<<"y_max="<<y_max<<endl;
    cout<<"dy   ="<<dy   <<endl;
    cout<<"z_min="<<z_min<<endl;
    cout<<"z_max="<<z_max<<endl;
    cout<<"dz   ="<<dz   <<endl;
    
    cout<<"n_xgrid="<<n_xgrid<<endl;
    cout<<"n_ygrid="<<n_ygrid<<endl;
    cout<<"n_zgrid="<<n_zgrid<<endl;
    
    cout<<"x grid:"<<endl;
    for(int ix=0; ix<n_xgrid; cout<<x[ix++]<<" ");  cout<<endl;
    cout<<"y grid:"<<endl;
    for(int iy=0; iy<n_ygrid; cout<<y[iy++]<<" ");  cout<<endl;
    cout<<"z grid:"<<endl;
    for(int iz=0; iz<n_zgrid; cout<<z[iz++]<<" ");  cout<<endl;
 
    //AWS20090805 removed this block since gamma grid is undefined at this point in the present version
    /*   
    cout<<"n_E_gammagrid="<<n_E_gammagrid<<endl;
    cout<<"E_gamma grid:"<<endl;
    for(int iEgamma=0; iEgamma<n_E_gammagrid; cout<<E_gamma[iEgamma++]<<" ");  cout<<endl;
    */

    cout<<"n_HI[0][0][0].s[0]="<<n_HI.d3[0][0][0].s[0]<<endl;
    cout<<"n_H2[0][0][0].s[0]="<<n_H2.d3[0][0][0].s[0]<<endl;
  
  } //(n_spatial_dimensions==3

}

Galaxy & Galaxy::operator = (const Galaxy & old) {
   //Avoid self-assignment
   if (this != &old) {
      free_memory();
      copy_variables(old);
      copy_pointers_deep(old);
   }
   return (*this);
}
