
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Galdef.h *                                    galprop package * 10/12/2003 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#ifndef Galdef_h
#define Galdef_h

//using namespace std; //AWS20050624

#include <cstdio>     //AWS20050624
#include <string>     //AWS20050624
#include <vector>

class Galdef {
   
 private:
  
  bool deleteMemory; // Only delete the memory pointed by our variables if we allocate it within the class.

 public:
  
  std::string version;    // galprop version number
  std::string run_no;    // identifier of galprop run
  std::string galdef_ID; // full identifier e.g. 01_123456
  
  int n_spatial_dimensions;             // 1,2 or 3
  double z_min,z_max,dz;                // for 1,2,3D    
  double r_min,r_max,dr;                // for 2D 
  double x_min,x_max,dx,y_min,y_max,dy; // for 3D 
  double p_min,p_max,p_factor;          // momentum start, end, factor
  double Ekin_min,Ekin_max,Ekin_factor; // kinetic energy per nucleon start, end, factor
  
  char p_Ekin_grid[5];                    // "p"||"Ekin": construct grid in p or Ekin
  
  double E_gamma_min,E_gamma_max,E_gamma_factor; // gamma-ray energy (MeV) start, end, factor
  double long_min,long_max;                      // gamma-ray skymap longitude min,max (degrees)
  double  lat_min, lat_max;                      // gamma-ray skymap latitude  min,max (degrees)
  double d_long,d_lat;                           // gamma-ray skymap longitude,latitude binsize (degrees)       
  int integration_mode;                          // integr.over the particle spectrum: =1-old E*logE; !=1-power-law analyt.
  int healpix_order;                             // The order for healpix maps.
  int lat_substep_number;                        // latitude bin splitting (0,1=no split, 2=split in 2...) IMOS20080114
  double LoS_step;                               // kpc, Line of Sight (LoS) integration step              IMOS20080114
  int LoS_substep_number;                        // number of substeps per LoS integration step for gas number density averaging (0,1=no substeps) IMOS20080114

  double nu_synch_min,nu_synch_max,nu_synch_factor;// synchrotron frequency min,max (Hz), factor
  
  double D0_xx;                         // diffusion coefficient at break  rigidity   
  double D_rigid_br;                    // break     rigidity for diffusion coefficient in MV
  double D_g_1;                         // diffusion coefficient index below break  rigidity
  double D_g_2;                         // diffusion coefficient index above break  rigidity 
  
  int diff_reacc;                       // 1,2=incl.diff.reacc.; 11=Kolmogorov+damping; 12=Kraichnan+damping
  double v_Alfven;                      // Alfven speed in km s-1
  double damping_p0;                    // ~ 1.e5 MV -characteristic rigidity         IMOS20030129
  double damping_max_path_L;            // ~ 1.*kpc2cm; Lmax~1 kpc, max free path     IMOS20030129
  double damping_const_G;               // a const derived from fitting B/C           IMOS20030129
  
  int convection;                       // 1=include convection                       AWS20010323
  double   v0_conv;                     // v=v0_conv+dvdz_conv*abs(z)                 AWS20010323
  double dvdz_conv;                     // v=v0_conv+dvdz_conv*abs(z)                 AWS20010323
  
  double rigid_min;                     // minimum rigidity of spectra
  double rigid_max;                     // maximum rigidity of spectra
  double nuc_rigid_br;                  // break rigidity for nuclei injection  in MV
  double nuc_g_1;                       // spectral index below break  
  double nuc_g_2;                       // spectral index above break 
  
  char inj_spectrum_type[9];            // "rigidity"||"beta_rig"||"Etot": choose the nucleon injection spectrum IMOS20000613.1
  
  double electron_g_0;                  // index below electron_rigid_br0            IMOS20031012
  double electron_rigid_br0;            // break rigidity0 for electron injection in MV
  double electron_g_1;                  // spectral index between breaks
  double electron_rigid_br;             // break rigidity1 for electron injection in MV
  double electron_g_2;                  // spectral index above electron_rigid_br
  
  double He_H_ratio;                    // He/H of ISM, by number
  int     n_X_CO;                       // option for X_CO values IMOS20080114
  double  X_CO;                         // conversion factor CO integrated temperature -> H2 column density IMOS20080114
  double *X_CO_radius, *X_CO_values;    // for n_X_CO=1
  double X_CO_parameters[10];           // for n_X_CO=2
  int n_X_CO_values;                    //   --- || ---
  int    propagation_X_CO;              // option to control H2 from CO for propagation  AWS20090623
  int    nHI_model;                     // selection of model for HI  gas density        AWS20090814
  int    nH2_model;                     // selection of model for H2  gas density        AWS20090814
  int    nHII_model;                    // selection of model for HII gas density        AWS20090814

  char    COR_filename[300];            // CO -molecular gas file IMOS20080114
  char    HIR_filename[300];            // HI -atomic    gas file IMOS20080114

  char    GCR_data_filename[300];       // Cosmic-ray database GJ20100103

  int fragmentation;                    // 1=include fragmentation
  int momentum_losses;                  // 1=include momentum losses
  int radioactive_decay;                // 1=include radioactive_decay
  int K_capture;                        // 1=include K_capture                        AWS20010731
  int ionization_rate;                  // 1=compute ionization rate          IMOS20060420
  
  double start_timestep;                // start time step in years
  double   end_timestep;                //   end time step in years
  double       timestep_factor;         //   factor to multiply timestep
  int          timestep_repeat;         //   number of times to repeat for factor
  int          timestep_repeat2;        //   number of  times to repeat in timestep_mode=2
  int          timestep_print ;         //   number of timesteps between printings
  int          timestep_diagnostics;    //   number of timesteps between propel diagnostics
  int           control_diagnostics;    //   control details of propel diagnostics
  
  int  network_iterations;              //   number of iterations of the protons (Needed for damping)
  int  network_iter_compl;              //   number of iterations for the entire network (For accurate B/C we need 2)
  int  network_iter_sec;                //   number of iterations for low A (<1) secondaries (usually only need 1, saves cpu time)
                                        //   network_iterations > network_iter_compl > network_iter_sec
  
  int prop_r;                           // for 2D: 1=propagate in r;
  int prop_x,prop_y,prop_z;             // for 2D: 1=propagate in z;for 3D 1=propagate in x,y,z
  int prop_p;                           // propagate in momentum
  int use_symmetry;                     // xyz symmetry (3D)
  
  int vectorized;                       // 0=unvectorized   code, 1=vectorized   code
  
  int source_specification;             // 0,1,2
  double source_normalization;          // 1.                                         IMOS20030129
  double electron_source_normalization; // 1.                                         IMOS20031016
  
  int  max_Z;                           // maximum   nucleus  Z in galdef file
  int *use_Z;                           // 1=use this nucleus Z
  
  double **isotopic_abundance;          // isotopic abundances
  
  int total_cross_section;              // total inel. cross-sections option AWS20010620
  int cross_section_option;             // controls which cross-sections used
  
  double t_half_limit;                  // lower limit on radioactive half-life for explicit inclusion  AWS20010731
  
  int primary_electrons;                // 1=propagate primary electrons
  int secondary_positrons;              // 1=propagate secondary positrons
  int secondary_electrons;              // 1=propagate secondary electrons
  int knock_on_electrons;               // 1=propagate knock-on electrons     IMOS20060504
  int tertiary_antiprotons;             // 1=propagate tertiary antiprotons   IMOS20000605.13
  int secondary_antiprotons;            // 1=propagate secondary antiprotons
  int secondary_protons;                // 1=propagate secondary protons      IMOS20000605.14
  
  int gamma_rays;                       // 1=compute gamma-ray emission
  int pi0_decay;                        // 1 - Dermer 1986 formalism, 2 - Blattnig et al. 2000,PRD 62,094030  IMOS20050216
  int IC_isotropic;                     // 1=compute isotropic inverse Compton IMOS20060420
  int IC_anisotropic;                   // 1=compute anisotropic inverse Compton
  int bremss;                           // 1=compute bremsstrahlung            IMOS20060420
  int synchrotron;                      // 1=compute synchrotron emission
  
  // DM parameters IMOS20050912
  int DM_positrons;                     // 1=compute positrons from DM
  int DM_electrons;                     // 1=compute electrons from DM
  int DM_antiprotons;                   // 1=compute antiprotons from DM
  int DM_gammas;                        // 1=compute gamma rays from DM
  double                                // user-defined params of DM (double)
    DM_double0, DM_double1, DM_double2, DM_double3, DM_double4,
    DM_double5, DM_double6, DM_double7, DM_double8, DM_double9;
  int                                   // user-defined params of DM (int)
    DM_int0, DM_int1, DM_int2, DM_int3, DM_int4,
    DM_int5, DM_int6, DM_int7, DM_int8, DM_int9;
  
  
  int source_model;                     //1= 2= 3= 5= 6=S&Mattox with cutoff
  int source_model_electron;
  std::vector<double> source_parameters; // for source_model=1
  std::vector<double> source_parameters_electron;
  double *source_radius, *source_values;// for source_model=8
  int n_source_values;                  //   --- || ---
  
  int   n_cr_sources;                   // number of pointlike cosmic-ray sources
  double *cr_source_x;                  // source x positions
  double *cr_source_y;                  // source y positions
  double *cr_source_z;                  // source z positions
  double *cr_source_w;                  // source width sigma in kpc
  double *cr_source_L;                  // source luminosity in TBD units
  
  int    SNR_events;                    // handle stochastic SNR events
  double SNR_interval;                  // time in years between SNRs in 1 kpc^3
  double SNR_livetime;                  // CR-producing live-time in years of an SNR
  double SNR_electron_sdg;              // delta electron source index Gaussian sigma
  double SNR_nuc_sdg;                   // delta nucleus  source index Gaussian sigma
  double SNR_electron_dgpivot;          // delta electron source index pivot rigidity (MeV)   AWS20010410
  double SNR_nuc_dgpivot;               // delta nucleus  source index pivot rigidity (MeV)   AWS20010410 
  
//  int HI_survey;                        // HI survey : 8=original 8 rings+high-latitudes, 9= 9 rings all-sky  AWS20050913
//  int CO_survey;                        // CO survey : 8=original 8 rings                 9= 9 rings all-sky  AWS20050913
  
  int     B_field_model;                // >1000=parameterized model
  char    B_field_name[100];            // 3D B-field model name ("galprop_original" uses B_field_model)      AWS20080313 
  int   n_B_field_parameters;           // number of B-field parameters                                       AWS20080313
  double *B_field_parameters;           // parameters for 3D B-field models                                   AWS20080313

  char ISRF_file[100];   // ISRF input file AWS20050301
  int ISRF_filetype;        // 0 for old (<Jan07), 1 for new FITS, 2 for new HEALPix
  int ISRF_healpixOrder; // Healpix binning for ISRF skymaps
  double ISRF_factors[3];  // ISRF factors for inverse Compton calculation       AWS20050301          
  
  
  double   proton_norm_Ekin;            // proton   kinetic energy for normalization (MeV)
  double   proton_norm_flux;            // flux of protons   at normalization energy (cm^-2 sr^-1 s^-1 MeV^-1)
  double electron_norm_Ekin;            // electron kinetic energy for normalization (MeV)
  double electron_norm_flux;            // flux of electrons at normalization energy (cm^-2 sr^-1 s^-1 MeV^-1)

  int skymap_format;                    // Select the output format for the skymap
  
  int output_gcr_full;                  // output full 2D or 3D gcr array
  int warm_start;                       // read nuclei file and continue run   AWS20010121
  
  int verbose;                          // verbosity: 0=min 10=max 
  int test_suite;                       // run test suit instead of normal run

  int primary_DM_positron;              // 1=propagate positrons produced by DM annihilation BiXJ 2007/2/1
  int primary_DM_electron;              // 1=propagate electrons produced by DM annihilation BiXJ 2007/2/1
  int primary_DM_antip;                 // 1=propagate antiprotons produced by DM annihilation BiXJ 2007/2/1
  char DM_positron_filename[300];
  char DM_electron_filename[300];
  char DM_antip_filename[300];
  
  //interface functions prototypes
  Galdef();
  ~Galdef();
  int read(const std::string& version, 
	   const std::string& runNumber,
	   const std::string& galdefDirectory);
	   //char *version_,char  *run_no_, const std::string& galdef_directory);//char *galdef_directory);
  int read_galdef_parameter(const std::string& filename, char *parstring,double *value);
  int read_galdef_parameter(const std::string& filename, char *parstring,double *value, bool); // AV 2010
  int read_galdef_parameter(const std::string& filename, char *parstring,double *value, int size);
  int read_galdef_parameter(const std::string& filename, char *parstring,int    *value);
  int read_galdef_parameter(const std::string& filename, char *parstring,char *value);
  void print();
};

#endif










