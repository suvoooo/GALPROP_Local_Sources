
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_bremss_emiss.cc *                       galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cassert>
#include <string>
#include <cstring>
#include <valarray>

using namespace std;//AWS20050624

#include"galprop_classes.h"
#include"galprop_internal.h"
#include "fitsio.h" 

#include <ErrorLogger.h>

void Galprop::store_IC_emiss(const string& type) {

  INFO("Entry");

  fitsfile* fptr = 0;

  const long nAxes = 4;

  long axes[nAxes];

  assert(galaxy.n_spatial_dimensions == 2 || galaxy.n_spatial_dimensions == 3);
  
  axes[0] = (galaxy.n_spatial_dimensions == 2 ? galaxy.n_rgrid : galaxy.n_xgrid);
  axes[1] = galaxy.n_zgrid;
  axes[2] = galaxy.n_E_gammagrid;
  axes[3] = galaxy.n_ISRF_components + 1;

  const long nElements = axes[0]*axes[1]*axes[2]*axes[3];

  valarray<float> array(0., nElements);

  const string filename = "!" + configure.fOutputDirectory + configure.fOutputPrefix + "IC_emiss_" + galdef.galdef_ID + ".gz";

  int status = 0;

  fits_create_file(&fptr, filename.c_str(), &status);

  fits_create_img(fptr, FLOAT_IMG, nAxes, axes, &status);

  // Write some keywords giving information about whether its an isotropic 
  // emissivity, anisotropic calculated for particular viewing location, etc.

  long isotropicIC = (type == "isotropic" ? 1 : 0);
  
  fits_update_key(fptr, TLONG, "ISOTROPIC", &isotropicIC, "Isotropic/Anisotropic cross section", &status);

  // for 3D case store x-dimension at y = 0 -- for now consistent with brem
  // emissivity storage etc. In future will export full 3D 
 
  for (unsigned int iComp = 0; iComp < axes[3]; ++iComp) {

    for (unsigned int iP = 0; iP < axes[2]; ++iP) {

      for (unsigned int iZ = 0; iZ < axes[1]; ++iZ) {

	for (unsigned int iR = 0; iR < axes[0]; ++iR) {

	  unsigned int index = iComp*axes[2]*axes[1]*axes[0] + iP*axes[1]*axes[0] + iZ*axes[0] + iR; // explicit encoding -- use it!!

	  if (iComp < axes[3] - 1) { // Individual components

	    if (type == "isotropic")
	      array[index] = (galaxy.n_spatial_dimensions == 2 ? galaxy.IC_iso_emiss[iComp].d2[iR][iZ].s[iP] : galaxy.IC_iso_emiss[iComp].d3[iR][0][iZ].s[iP]);
	    else
	      array[index] = (galaxy.n_spatial_dimensions == 2 ? galaxy.IC_aniso_emiss->d2[iR][iZ].s[iP] : galaxy.IC_aniso_emiss->d3[iR][0][iZ].s[iP]);

	  } else { // Total emissivity

	    double sum = 0;

	    for (unsigned int iC = 0; iC < axes[3]-1; ++iC)
	      if (type == "isotropic")
		sum += (galaxy.n_spatial_dimensions == 2 ? galaxy.IC_iso_emiss[iC].d2[iR][iZ].s[iP] : galaxy.IC_iso_emiss[iC].d3[iR][0][iZ].s[iP]);
	      else
		sum += (galaxy.n_spatial_dimensions == 2 ? galaxy.IC_aniso_emiss->d2[iR][iZ].s[iP] : galaxy.IC_aniso_emiss->d3[iR][0][iZ].s[iP]);

	    array[index] = sum;

	  }

	  array[index] *= galaxy.E_gamma[iP]*galaxy.E_gamma[iP];

	}

      }

    }

  }
  
  // Write the array of floats to the image 
  long fPixel = 1;
  fits_write_img(fptr, TFLOAT, fPixel, nElements, &array[0], &status);
  
  // write basic FITS keywords
  
  double crval1 = (galaxy.n_spatial_dimensions == 2 ? galaxy.r_min : galaxy.x_min);
  double crval2 = galaxy.z_min;
  double crval3 = log10(galaxy.E_gamma_min);
  double crval4 = 1;

  double cdelt1 = (galaxy.n_spatial_dimensions == 2 ? galaxy.dr : galaxy.dx);
  double cdelt2 = galaxy.dz;
  double cdelt3 = log10(galaxy.E_gamma_factor);
  double cdelt4 = 1;
  
  fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1,"Start of radial dimension (kpc)", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2,"Start of Z dimension (kpc)", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL3", &crval3,"Start of log10(energy grid/MeV)", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL4", &crval4,"Start of component numbering", &status);
  
  fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1,"Increment of radial dimension (kpc)", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2,"Increment of Z dimension (kpc)", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT3", &cdelt3,"Increment of log10(energy grid/MeV)", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT4", &cdelt4,"Increment of components", &status);
  
  fits_close_file(fptr, &status);   
  fits_report_error(stderr, status);

  INFO("Exit");

}
