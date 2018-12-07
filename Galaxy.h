
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Galaxy.h *                                    galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#ifndef Galaxy_h
#define Galaxy_h

#include"Distribution.h"
//#include <GalacticRadiationField.h>
#include "Skymap.h"

#include <valarray>

class Galaxy {

 public:

  Galaxy();
  Galaxy(const Galaxy&);
  ~Galaxy();

  double z_min, z_max, dz;                // for 1,2,3D    
  double r_min, r_max, dr;                // for 2D 
  double x_min, x_max, dx, y_min, y_max, dy; // for 3D 
  
  double p_min, p_max, p_factor; // momentum start, end, factor
   
  int n_spatial_dimensions;// 1,2,3D
  int n_pgrid;             // number of points in momentum
  int n_rgrid;             // number of points in radius (2D)   
  int n_zgrid;             // number of points in z (1D,2D)  
  int n_xgrid;             // number of points in x (3D)
  int n_ygrid;             // number of points in y (3D)    
  
  double* x;             // x grid
  double* y;             // y grid
  double* z;             // z grid 
  double* r;             // r grid 
  
  int n_E_gammagrid;     // number of points in gamma-ray energy
  double E_gamma_min,E_gamma_max,E_gamma_factor;// min,max, factor for gamma-ray energy
  double* E_gamma;       // gamma ray energy grid  
  int n_lat,n_long;      // dimensions of gamma-ray skymaps
  double long_min,long_max; // longitude range of gamma-ray skymaps
  double lat_min, lat_max; // latitude range of gamma-ray skymaps
  double d_long,d_lat;     // longitude, latitude binsize of  gamma-ray skymaps
  
  int n_nu_synchgrid;    // number of points in synchrotron frequency
  double nu_synch_min,nu_synch_max,nu_synch_factor; // min,max, factor for synchrotron frequency
  double* nu_synch;                     // synchrotron frequency grid
  
  Distribution n_H2;
  Distribution n_HI;
  Distribution n_HII;
  
  Distribution HIR; // HI column density map in Galactocentric rings AWS2001025
  Distribution COR; // CO column density map in Galactocentric rings AWS2001025

  int n_Ring;   // the number of gas rings IMOS20080114
  float* R_bins;// pointer to an array of ring radii IMOS20080114

  Skymap<double> hpCOR;  // CO column density healpix map in Galactocentric rings 
  Skymap<double> hpHIR;  // HI column density healpix map in Galactocentric rings 
 
  Distribution B_field;
  Distribution* ISRF;
  Distribution* ISRF_energy_density;
  int n_ISRF_components;

  std::valarray<double> fISRFFactors;

  //rf::GalacticRadiationField* fISRFModel;

  valarray<double> nu_ISRF;       // ISRF frequencies
  
  Distribution SNR_cell_time;           // time between SNR for each cell
  Distribution SNR_cell_phase;          // phase of SNR for each cell
  Distribution SNR_electron_dg;         // electron injection spectral index delta (Gaussian distributed) AWS20010410 
  Distribution SNR_nuc_dg;              // nucleus  injection spectral index delta (Gaussian distributed) AWS20010410
  
  Distribution  bremss_emiss;           // bremsstrahlung emissivity on neutral gas
  Distribution  bremss_ionized_emiss;   // bremsstrahlung emissivity on ionized gas
  Distribution* IC_iso_emiss;           // inverse Compton isotropic emissivity for each ISRF component
  Distribution* IC_aniso_emiss;
  Distribution  pi0_decay_emiss;        // pi0 decay emissivity
  
  Distribution  bremss_skymap;          // bremsstrahlung intensity skymap on neutral gas
  Distribution  bremss_ionized_skymap;  // bremsstrahlung intensity skymap on HII in rings 
  Distribution* IC_iso_skymap;          // inverse Compton   isotropic intensity skymap for each ISRF component
  Distribution* IC_aniso_skymap;        // inverse Compton anisotropic intensity skymap for each ISRF component IMOS20060420
  Distribution  pi0_decay_skymap;       // pi0 decay  intensity skymap
  
  Distribution     bremss_H2R_skymap;   // bremsstrahlung intensity skymap on CO  in rings                AWS20041214
  Distribution     bremss_HIR_skymap;   // bremsstrahlung intensity skymap on HI  in rings                AWS20041214
  Distribution     bremss_HII_skymap;   // bremsstrahlung intensity skymap on HII  in rings               IMOS20080114*
  Distribution  pi0_decay_H2R_skymap;   // pi0 decay      intensity skymap on H2  in rings                AWS20041214
  Distribution  pi0_decay_HIR_skymap;   // pi0 decay      intensity skymap on HI  in rings                AWS20041214
  Distribution  pi0_decay_HII_skymap;   // pi0 decay      intensity skymap on HII in rings                IMOS20080114*

  Skymap<double>  bremss_hp_skymap;          // bremsstrahlung intensity skymap on neutral gas
  Skymap<double>  bremss_ionized_hp_skymap;  // bremsstrahlung intensity skymap on ionized gas
  Skymap<double>* IC_iso_hp_skymap;          // inverse Compton   isotropic intensity skymap for each ISRF component
  Skymap<double>* IC_aniso_hp_skymap;        // inverse Compton anisotropic intensity skymap for each ISRF component IMOS20060420
  Skymap<double>  pi0_decay_hp_skymap;       // pi0 decay  intensity skymap
  
  Skymap<double>*     bremss_H2R_hp_skymap;   // bremsstrahlung intensity skymap on CO  in rings                 AWS20041214
  Skymap<double>*     bremss_HIR_hp_skymap;   // bremsstrahlung intensity skymap on HI  in rings                 AWS20041214
  Skymap<double>*     bremss_HII_hp_skymap;   // bremsstrahlung intensity skymap on HII  in rings                IMOS20080114*
  Skymap<double>*  pi0_decay_H2R_hp_skymap;   // pi0 decay      intensity skymap on H2  in rings                 AWS20041214
  Skymap<double>*  pi0_decay_HIR_hp_skymap;   // pi0 decay      intensity skymap on HI  in rings                 AWS20041214
  Skymap<double>*  pi0_decay_HII_hp_skymap;   // pi0 decay      intensity skymap on HII in rings                 IMOS20080114*
  
  
  Distribution synchrotron_emiss;         // synchrotron emissivity
  Distribution synchrotron_Q_emiss;       // synchrotron Stokes Q emissivity                                       AWS20100708
  Distribution synchrotron_U_emiss;       // synchrotron Stokes U emissivity                                       AWS20100708

  Distribution synchrotron_skymap;        // synchrotron total intensity     Cartesian skymap
  Distribution synchrotron_Q_skymap;      // synchrotron Stokes Q            Cartesian skymap                      AWS20100708
  Distribution synchrotron_U_skymap;      // synchrotron Stokes U            Cartesian skymap                      AWS20100708
  Distribution synchrotron_P_skymap;      // synchrotron polarized intensity Cartesian skymap                      AWS20100708

  Skymap<double> synchrotron_hp_skymap;   // synchrotron total intensity       healpix skymap
  Skymap<double> synchrotron_Q_hp_skymap; // synchrotron Stokes Q              healpix skymap                      AWS20100708
  Skymap<double> synchrotron_U_hp_skymap; // synchrotron Stokes U              healpix skymap                      AWS20100708
  Skymap<double> synchrotron_P_hp_skymap; // synchrotron polarized intensity   healpix skymap                      AWS20100708

  Distribution ionization_rate;         // cosmic-ray ionization rate 
  
  Distribution DM_emiss;                // DM gamma-ray emission distribution IMOS20050912
  Distribution DM_skymap;               // DM gamma-ray emission skymap       IMOS20050912
  Skymap<double> DM_hp_skymap;
  
  //interface functions prototypes
  void init( double r_min_, double r_max_, double dr_, 
	     double z_min_, double z_max_, double dz_); 
  
  void init( double x_min_, double x_max_, double dx_, 
	     double y_min_, double y_max_, double dy_, 
	     double z_min_, double z_max_, double dz_);
  
  void print();

  Galaxy & operator = (const Galaxy &);

 private:
  void copy_variables(const Galaxy &);
  void copy_pointers_deep(const Galaxy &);
  void free_memory();
  void init_pointers();

};

#endif

 
