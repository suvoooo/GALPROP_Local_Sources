//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_pi0_decay_skymap.cc *                    galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cassert>
#include <string>
#include <cstring>
#include <valarray>

using namespace std;

#include "galprop_classes.h"
#include "galprop_internal.h"
#include "fitsio.h"
#include "SkymapFitsio.h"

#include <ErrorLogger.h>

int Galprop::store_pi0_decay_skymap() {

  INFO("Entry");

  int status = 0;

  if (3 == galdef.skymap_format) {
    
    string filename = configure.fOutputDirectory;
    filename += configure.fOutputPrefix;
    filename += "pi0_decay_healpix_";
    filename += galdef.galdef_ID;
    filename += ".gz";
    SkymapToFits(galaxy.pi0_decay_hp_skymap, filename, "Energy", "MeV");

  } else {
    
    //fitsfile *fptr;       // pointer to the FITS file; defined in fitsio.h

    const int naxis = 4;
    long naxes[naxis];
    
    naxes[0] = galaxy.n_long;
    naxes[1] = galaxy.n_lat;             
    naxes[2] = galaxy.n_E_gammagrid;
    naxes[3] = 1;
    
    const long nElements = naxes[0]*naxes[1]*naxes[2]*naxes[3];
    
    valarray<float> array(0., nElements);

    int i = 0; 
    
    if (galdef.skymap_format != 1) { 

      double crval[naxis], cdelt[naxis];

      crval[0] = galaxy.long_min;
      crval[1] = galaxy. lat_min;
      crval[2] = log10(galaxy.E_gamma_min);
      crval[3] = 1;
      
      cdelt[0] = galaxy.d_long;
      cdelt[1] = galaxy.d_lat;
      cdelt[2] = log10(galaxy.E_gamma_factor);
      cdelt[3] = 1;
    
      i = 0;

      for (int ip = 0; ip < naxes[2]; ++ip)
	for (int ib = 0; ib < naxes[1]; ++ib)
	  for (int il = 0; il < naxes[0]; ++il) {
	  
	  array[i] = galaxy.pi0_decay_skymap.d2[il][ib].s[ip];
	  array[i] *= pow(galaxy.E_gamma[ip], 2.);
	  ++i;

	  }

      status = store_skymap(&array[0], naxes, "pion_decay_skymap_", crval, cdelt);

    }

    if ((galdef.skymap_format == 1) || 
	(galdef.skymap_format == 2)) { //!< Mapcube output compatable with Glast science tools
      
      i = 0;

      for (int ip = 0; ip < naxes[2]; ++ip)
        for (int ib = 0; ib < naxes[1]; ++ib)
	  for (int il = 0; il < naxes[0]; ++il) {
	  
	    array[i] = galaxy.pi0_decay_skymap.d2[il][ib].s[ip];
	    ++i;

	}

      status = store_mapcube_skymap(&array[0], galaxy.E_gamma, 1, galaxy.n_E_gammagrid, "pion_decay_mapcube_", true);
    }

  }
  
  INFO("Exit");
  
  return status;

}
