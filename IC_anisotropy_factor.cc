
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * IC_anisotropy_factor.cc *                     galprop package * 4/20/2006
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"

#include <fort_interface.h>
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// call a fortran routine cfactor to calculate the ratio anisoIC/isoIC
// whole routine has been changed IMOS20060420
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
double Galprop::IC_anisotropy_factor (double E, double rho,double xi, double z, int isrf_comp) 
{
//cout<<" >>>> IC_anisotropy_factor"<<endl;
//cout<<"IC_anisotropy_factor:isrf_comp="<<isrf_comp<<" rho="<<rho<<" xi="<<xi<<" z="<<z<<endl;

  if(isrf_comp==2) {return 1.0;}                  //      CMB,                  isrf_comp=2
  double           E0=0.01e-6;                    // MeV, dust      <E>=0.01 eV,isrf_comp=1
  if(isrf_comp==0) E0=0.90e-6;                    //      starlight <E>=0.90 eV,isrf_comp=0
  int kbg=1;                                      // transparent disk
  double PLindex=3;                               // power-law index of the electron spectrum
  double RG=15.;                                  // kpc,radius of Galactic disk
  double RS=8.5;                                  // kpc,distance of the sun from Galactic center  
  
  isrf_energy_density_i_comp=isrf_comp;           // as required by isrf_energy_density.cc
  double factor=cfactor_cc(kbg,E0/Mele,E/Mele,PLindex,RG,rho,xi,z,RS);
  
//cout<<"isrf_comp="<<isrf_comp<<" rho="<<rho<<" xi="<<xi<<" z="<<z<<" cfactor="<<factor<<endl;
//cout<<" <<<< IC_anisotropy_factor"<<endl;
  return factor;
}
/*
      real*8 function CFactor(kbg,E0,E,PLindex,RG,rho,xi,z,RS)
c***********************************************************************
c                            *** I.Moskalenko, version of 31.03.1998 ***
c calculation of the correction factor for
c INVERSE COMPTON scattering in an ANISOTROPIC photon field
c    INPUT:
c kbg =-1 for isotropic background photon field (constant); 
c     = 0 for opaque disk photon field (~cos(theta));
c     = 1 for transparent disk photon field (~1/cos(theta)); 
c E0 is the energy of a background photon (in electron rest mass units mc^2);
c E  is the energy of scattered photon (in electron rest mass units mc^2);
c PLindex is the power-law index of the electron spectrum (>0);
c RG is the Galactic radius (kpc);
c (rho,xi,z) are the cylindrical coordinates of the electron position
c   in respect to the Galactic center (rho,z in kpc, 0< xi (radians) <2Pi 
c   is polar angle between the electron position and observer direction);
c RS is the Solar distance from the Galactic center (kpc).
c    internal OUTPUT:
c SPEC - the differential energy spectrum of gamma-rays (phot/ccm/sec/energy)
c    per one electron in ccm divided by Pi*r0^2*c;
c DENS - number density of background photons at the electron position
c    (the emissivity of unit area of the Galactic plane is equal to 1).
c    REFERENCES:
c Moskalenko I.V., Strong A.W. 2000, ApJ 528, 357
c***********************************************************************
*/
