#include <ErrorLogger.h>

#include <galprop_classes.h>
#include <galprop_internal.h>

#include <config.h>

extern Galprop* gGalprop;

//IMOS20060420 transferred from fort_interface1.cc
// the routine is called by FORTRAN routine emiss(r) (in file cfactor.f)
// hence underscore is appended and extern "C" supplied

#define ISRF_ENERGY_DENSITY_F77 F77_FUNC_(isrf_energy_density,ISRF_ENERGY_DENSITY)
#ifdef __cplusplus
extern "C" 
#endif
void ISRF_ENERGY_DENSITY_F77(float* r, float* z, float* energy_density) {

  *energy_density = gGalprop->isrf_energy_density(*r,*z);
//   cout<<"energy_density = "<<*energy_density<<endl;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// IMOS20060901
// a dummy routine to be replaced by DarkSUSY when compiled with GALPROP 

#define RHO_DARKSUSY_F77 F77_FUNC_(rho_darksusy,RHO_DARKSUSY)
#ifdef __cplusplus
extern "C" 
#endif
void RHO_DARKSUSY_F77(double* Xkpc, double* Ykpc, double* Zkpc, double* rho0) {

  *rho0 = -1.;
  
}

