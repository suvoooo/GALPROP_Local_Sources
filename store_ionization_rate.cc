
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_ionization_rate.cc *                    galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cassert>
#include <string>
#include <cstring>
#include <valarray>

using namespace std;//AWS20050624

#include "galprop_classes.h"
#include "galprop_internal.h"

#include "fitsio.h" 

#include <ErrorLogger.h>

int Galprop::store_ionization_rate() {

  INFO("Entry");

  int status = 0;

  fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
  long  fpixel = 1, naxis = 4, nelements, exposure;
  long naxes[4]  ; 
    
  if(galaxy.n_spatial_dimensions==2)naxes[0]=galaxy.n_rgrid;
  if(galaxy.n_spatial_dimensions==3)naxes[0]=galaxy.n_xgrid;
  naxes[1]=galaxy.n_zgrid;             
  naxes[2]=1;
  naxes[3]=1;
  
  //cout<<galaxy.n_spatial_dimensions<<" "<<naxes[0]<<" "<<naxes[1]<<"  "<<naxes[2]<<endl;

  nelements=naxes[0]*naxes[1]*naxes[2]*naxes[3];
  
  valarray<float> array(0., nelements);

  const std::string outfile = "!" + configure.fOutputDirectory + configure.fOutputPrefix + "ionization_rate_" + galdef.galdef_ID + ".gz";
  
  fits_create_file(&fptr, outfile.c_str(), &status);   /* create new file or overwrite existing one */
  
  /* Create the primary array image (32-bit float pixels */
  fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);
  
  /* Write a keyword; must pass the ADDRESS of the value */
  exposure = 1500;
  fits_update_key(fptr, TLONG, "EXPOSURE", &exposure, "Total Exposure Time", &status);
 
  // for 3D case store x-dimension at y=0
  
  int i=0;
    
  for (int iz       =0;        iz<naxes[1];       iz++){
    for (int ir       =0;        ir<naxes[0];       ir++){
      
      if(galaxy.n_spatial_dimensions==2)
	array[i]=galaxy.ionization_rate.d2[ir]   [iz].s[0];
      if(galaxy.n_spatial_dimensions==3)
	array[i]=galaxy.ionization_rate.d3[ir][0][iz].s[0];
      
      i++;
    }
  }
  
  /* Write the array of floats to the image */
  fits_write_img(fptr, TFLOAT, fpixel, nelements, &array[0], &status);
    
  // write basic FITS keywords
  
  double crval1,crval2,crval3,crval4;
  double cdelt1,cdelt2,cdelt3,cdelt4;
  
  if(galaxy.n_spatial_dimensions==2)crval1=galaxy.r_min;
  if(galaxy.n_spatial_dimensions==3)crval1=galaxy.x_min;
  crval2=galaxy.z_min;
  crval3=1;
  crval4=1;
  
  if(galaxy.n_spatial_dimensions==2)cdelt1=galaxy.dr;
  if(galaxy.n_spatial_dimensions==3)cdelt1=galaxy.dx;
  cdelt2=galaxy.dz;
  cdelt3=1;
  cdelt4=1;
  
  fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1,"Start of axis 1", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2,"Start of axis 2", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL3", &crval3,"Start of axis 3", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL4", &crval4,"Start of axis 4", &status);
  
  fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1,"Increment of axis 1", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2,"Increment of axis 2", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT3", &cdelt3,"Increment of axis 3", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT4", &cdelt4,"Increment of axis 4", &status);
  
  fits_close_file(fptr, &status);            /* close the file */
  
  fits_report_error(stderr, status);  /* print out any error messages */

  INFO("Exit");

  return status;

}
