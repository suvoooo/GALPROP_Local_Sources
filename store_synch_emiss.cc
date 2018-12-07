//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_synch_emiss.cc *                    galprop package *  
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cassert>
#include <string>
#include <cstring>
#include <valarray>

using namespace std;

#include "galprop_classes.h"
#include "galprop_internal.h"

#include "fitsio.h" 

#include <ErrorLogger.h>

int Galprop::store_synch_emiss() {                //AWS20080314

  INFO("Entry");
  
  int status; //AWS20100708 moved here

  int naxis = 4;
  long naxes[naxis]; 

  if (2 == galaxy.n_spatial_dimensions) { // all frequencies for 2D
    
    naxes[0] = galaxy.n_rgrid;
    naxes[1] = galaxy.n_zgrid;             
    naxes[2] = galaxy.n_nu_synchgrid;
    naxes[3] = 1;
    
  }

  if (3 == galaxy.n_spatial_dimensions) { // first  frequency only for 3D
   
    naxes[0] = galaxy.n_xgrid;
    naxes[1] = galaxy.n_ygrid;             
    naxes[2] = galaxy.n_zgrid;      
    naxes[3] = 1; // in case frequency required in future
    
  }

  const long nElements = naxes[0]*naxes[1]*naxes[2]*naxes[3];

  valarray<double> array(0., nElements);

  // for 3D case store x-dimension at y=0
  
  for(int stokes=0;stokes<=2;stokes++)//AWS20100708
  {

  int i = 0;
 
  if (2 == galaxy.n_spatial_dimensions) {

    for (int ip = 0; ip < naxes[2]; ++ip)
      for (int iz = 0; iz < naxes[1]; ++iz)
	for (int ir = 0; ir < naxes[0]; ++ir) {

	  if(stokes==0)array[i] = galaxy.synchrotron_emiss  .d2[ir][iz].s[ip];//AWS20100708
	  if(stokes==1)array[i] = galaxy.synchrotron_Q_emiss.d2[ir][iz].s[ip];//AWS20100708
	  if(stokes==2)array[i] = galaxy.synchrotron_U_emiss.d2[ir][iz].s[ip];//AWS20100708
	  ++i;

	}

  }
   
  if (3 == galaxy.n_spatial_dimensions) {

    for (int iz = 0; iz < naxes[2]; ++iz)
      for (int iy = 0; iy < naxes[1]; ++iy)
	for (int ix = 0; ix < naxes[0]; ++ix) {
 
	  if(stokes==0)array[i] = galaxy.synchrotron_emiss  .d3[ix][iy][iz].s[0];//AWS20100708
	  if(stokes==1)array[i] = galaxy.synchrotron_Q_emiss.d3[ix][iy][iz].s[0];//AWS20100708
	  if(stokes==2)array[i] = galaxy.synchrotron_U_emiss.d3[ix][iy][iz].s[0];//AWS20100708
	  ++i;

	}

  }

   string outfile; 
  if(stokes==0)                  outfile = "!" + configure.fOutputDirectory + configure.fOutputPrefix + "synchrotron_emiss_"   + galdef.galdef_ID + ".gz";//AWS20100708
  if(stokes==1)                  outfile = "!" + configure.fOutputDirectory + configure.fOutputPrefix + "synchrotron_Q_emiss_" + galdef.galdef_ID + ".gz";//AWS20100708
  if(stokes==2)                  outfile = "!" + configure.fOutputDirectory + configure.fOutputPrefix + "synchrotron_U_emiss_" + galdef.galdef_ID + ".gz";//AWS20100708

  status = 0; //AWS20100708
  
  fitsfile* fptr = 0;
 
  fits_create_file(&fptr, outfile.c_str(), &status);   /* create new file or overwrite existing one */

  /* Create the primary array image (64-bit double pixels */
  fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);

  long  fpixel = 1;

  fits_write_img(fptr, TDOUBLE, fpixel, nElements, &array[0], &status);

  // write basic FITS keywords

  double crval1, crval2, crval3, crval4;
  double cdelt1, cdelt2, cdelt3, cdelt4;
 
  if (2 == galaxy.n_spatial_dimensions) {

    crval1 = galaxy.r_min;
    crval2 = galaxy.z_min;
    crval3 = log10(galaxy.nu_synch_min);
    
    cdelt1 = galaxy.dr;
    cdelt2 = galaxy.dz;
    cdelt3 = log10(galaxy.nu_synch_factor);
    
  }
                                  
  if (3 == galaxy.n_spatial_dimensions) {

    crval1 = galaxy.x_min;
    crval2 = galaxy.y_min;
    crval3 = galaxy.z_min;
    
    cdelt1 = galaxy.dx;
    cdelt2 = galaxy.dy;
    cdelt3 = galaxy.dz;
    
  }

  fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1,"Start of axis 1", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2,"Start of axis 2", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL3", &crval3,"Start of axis 3", &status);
  //  fits_update_key(fptr, TDOUBLE, "CRVAL4", &crval4,"Start of axis 4", &status);
  
  fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1,"Increment of axis 1", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2,"Increment of axis 2", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT3", &cdelt3,"Increment of axis 3", &status);
    //  fits_update_key(fptr, TDOUBLE, "CDELT4", &cdelt4,"Increment of axis 4", &status);

  fits_close_file(fptr, &status);            /* close the file */

  fits_report_error(stderr, status);  /* print out any error messages */

  } //AWS20100708

  INFO("Exit");

  return status;

}
