//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_pi0_decay_emiss.cc *                    galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cassert>
#include <cstring>
#include <string>
#include <valarray>

using namespace std;

#include "galprop_classes.h"
#include "galprop_internal.h"

#include "fitsio.h" 

#include <ErrorLogger.h>

int Galprop::store_pi0_decay_emiss() {

  INFO("Entry");

  int naxis = 4;
  long naxes[naxis];

  if (2 == galaxy.n_spatial_dimensions)
    naxes[0]=galaxy.n_rgrid;
    
  if (3 == galaxy.n_spatial_dimensions)
    naxes[0]=galaxy.n_xgrid;
    
  naxes[1] = galaxy.n_zgrid;             
  naxes[2] = galaxy.n_E_gammagrid;
  naxes[3] = 1;

  long nElements = naxes[0]*naxes[1]*naxes[2]*naxes[3];

  valarray<float> array(0., nElements);

  // for 3D case store x-dimension at y=0
  
  int i = 0;
  
  for (int ip = 0; ip < naxes[2]; ++ip) {

    for (int iz = 0; iz < naxes[1]; ++iz) {
      
      for (int ir = 0; ir < naxes[0]; ++ir) {
	
	if (2 == galaxy.n_spatial_dimensions)
	  array[i] = galaxy.pi0_decay_emiss.d2[ir][iz].s[ip];

	if (3 == galaxy.n_spatial_dimensions)
	  array[i] = galaxy.pi0_decay_emiss.d3[ir][0][iz].s[ip];

	array[i] *= pow(galaxy.E_gamma[ip], 2.);
	++i;

      }

    }

  }
   
  const std::string outfile = "!" + configure.fOutputDirectory + configure.fOutputPrefix + "pion_decay_emiss_" + galdef.galdef_ID + ".gz";
  
  fitsfile* fptr = 0;
  int status = 0;    
     
  fits_create_file(&fptr, outfile.c_str(), &status);   /* create new file or overwrite existing one */
  
  /* Create the primary array image (32-bit float pixels */
  fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);

  /* Write a keyword; must pass the ADDRESS of the value */
  float exposure = 1500;
  fits_update_key(fptr, TLONG, "EXPOSURE", &exposure, "Total Exposure Time", &status);

  long fpixel = 1;

  /* Write the array of floats to the image */
  fits_write_img(fptr, TFLOAT, fpixel, nElements, &array[0], &status);

  // write basic FITS keywords

  double crval1, crval2, crval3, crval4;
  double cdelt1, cdelt2, cdelt3, cdelt4;
 
  if (2 == galaxy.n_spatial_dimensions)
    crval1 = galaxy.r_min;

  if (3 == galaxy.n_spatial_dimensions)
    crval1 = galaxy.x_min;

  crval2 = galaxy.z_min;
  crval3 = log10(galaxy.E_gamma_min);
  crval4 = 1;
 
  if (2 == galaxy.n_spatial_dimensions)
    cdelt1 = galaxy.dr;
  
  if (3 == galaxy.n_spatial_dimensions)
    cdelt1 = galaxy.dx;
             
  cdelt2 = galaxy.dz;
  cdelt3 = log10(galaxy.E_gamma_factor);
  cdelt4 = 1;

  fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1,"Start of axis 1", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2,"Start of axis 2", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL3", &crval3,"Start of axis 3", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL4", &crval4,"Start of axis 4", &status);

  fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1,"Increment of axis 1", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2,"Increment of axis 2", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT3", &cdelt3,"Increment of axis 3", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT4", &cdelt4,"Increment of axis 4", &status);

  fits_close_file(fptr, &status);     

  fits_report_error(stderr, status); 

  INFO("Exit");

  return status;

}
