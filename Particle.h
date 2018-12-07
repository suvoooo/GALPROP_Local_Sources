
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Particle.h *                                  galprop package * 4/14/2000
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#ifndef Particle_h
#define Particle_h

#include<iostream>   //AWS20050624
#include<cmath>      //AWS20050624

#include"Distribution.h"

using namespace std;
  
class Particle {

 public: 
 
  char name[100];
  int Z,A;
  int K_electron;           // number of K-electrons (0 or 1)      AWS20010731
  double mass;
  double t_half;            // half-life in years
  double primary_abundance; // primary isotopic abundance
  
  double z_min,z_max,dz;                // for 1,2,3D    
  double r_min,r_max,dr;                // for 2D 
  double x_min,x_max,dx,y_min,y_max,dy; // for 3D 
  
  double    p_min,    p_max,    p_factor; // momentum start, end, factor
  double Ekin_min, Ekin_max, Ekin_factor; // kinetic energy/nucleon start,end,factor
  
  int n_spatial_dimensions;// 1,2,3D
  int n_pgrid;             // number of points in momentum
  int n_rgrid;             // number of points in radius (2D)   
  int n_zgrid;             // number of points in z (1D,2D)  
  int n_xgrid;             // number of points in x (3D)
  int n_ygrid;             // number of points in y (3D)    
  
  double *x;             // x grid
  double *y;             // y grid
  double *z;             // z grid 
  double *r;             // r grid 
  
  double *p;             // total momentum of particle in MV
  double *Etot;          // total energy   of particle in MeV
  double *Ekin;          // kinetic energy per nucleon in MeV
  double *beta;          // velocity/c
  double *gamma;         // Etot/mass
  double *rigidity;      // rigidity p/Z in MV  

  double m_dt;               //YO/20140305  time_step
  double t_d;                //YO/20140310  propagating time
  double m_step_total_time;  //YO/20140611  propagating time
  double m_final_time;       //YO/20140325  total propagate time

  double m_sourceX,m_sourceY,m_sourceZ;  //YO/20141213 single source coordinate
  double m_source_age;                   //YO/20141224 single source age

  double m_nor_fac;          //YO/20141117  normalization factor in burst-like model

  char species[20];
  
  int arrays_assigned;
  
  double normalization_factor;// to normalize protons or electrons AWS20010121
  
  Distribution  cr_density;
  Distribution  primary_source_function;
  Distribution  secondary_source_function;
  Distribution  nearby_source_function;//SB20161712
  Distribution fragment;// fragmentation destruction rate
  Distribution decay   ;// radioactive decay rate        
  Distribution dpdt;    // momentum change rate       
  
  Distribution Dxx ;    // spatial  diffusion coefficient
  Distribution Dyy ;    // spatial  diffusion coefficient AWS20090625
  Distribution Dzz ;    // spatial  diffusion coefficient AWS20090625

  // Distribution Dxx_TeV ; // spatial  diffusion coefficient in TeV region //YO/20140422
  
  Distribution Dpp ;    // momentum diffusion coefficient
  Distribution v_conv;  // convection velocity in z-direction
  
  //interface functions prototypes
  Particle(){init();}
  Particle(const Particle &old);
  ~Particle();
  void init(char *name_,int Z_,int A_,int K_electron_,double t_half_, //AWS20010731
	    double r_min_, double r_max_, double dr_, 
	    double z_min_, double z_max_, double dz_,
	    double    p_min_,double    p_max_,double    p_factor_,
	    double Ekin_min_,double Ekin_max_,double Ekin_factor_,
	    char *p_Ekin_grid);
  void init(char *name_,int Z_,int A_,int K_electron_,double t_half_, //AWS20010731
	    double x_min_, double x_max_, double dx_, 
	    double y_min_, double y_max_, double dy_, 
	    double z_min_, double z_max_, double dz_,
	    double    p_min_,double    p_max_,double    p_factor_,
	    double Ekin_min_,double Ekin_max_,double Ekin_factor_,
	    char *p_Ekin_grid);
  void init();     // flags object as having arrays not yet assigned 
  int delete_transport_arrays();
  int create_transport_arrays();
  int delete_arrays();  //Gulli20070810
  void print();
  Particle &operator=(const Particle &particle);
};
  
#endif
