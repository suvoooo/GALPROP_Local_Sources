//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_synch_skymap.cc *                       galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cassert>
#include <cstring>
#include <string>
#include <sstream>
#include <valarray>

#include "galprop_classes.h"
#include "galprop_internal.h"

using namespace std;

#include "SkymapFitsio.h"

#include "fitsio.h" 

#include <ErrorLogger.h>

int Galprop::store_synch_skymap() { //AWS20050817

  INFO("Entry");

  int status = 0;

  if (3 == galdef.skymap_format) {

    string filename;                                                                                                                                //AWS20100709
    filename =       configure.fOutputDirectory + configure.fOutputPrefix + "synchrotron_healpix_" + galdef.galdef_ID + ".gz";                      //AWS20100107
    SkymapToFits(galaxy.synchrotron_hp_skymap, filename, "Frequency", "Hz");

    filename =       configure.fOutputDirectory + configure.fOutputPrefix + "synchrotron_Q_healpix_" + galdef.galdef_ID + ".gz";                    //AWS20100709
    SkymapToFits(galaxy.synchrotron_Q_hp_skymap, filename, "Frequency", "Hz");                                                                      //AWS20100709
    
    filename =       configure.fOutputDirectory + configure.fOutputPrefix + "synchrotron_U_healpix_" + galdef.galdef_ID + ".gz";             //AWS20100709
    SkymapToFits(galaxy.synchrotron_U_hp_skymap, filename, "Frequency", "Hz");                                                                      //AWS20100709


  } else {

    long naxes[4]; 
    double crval[4], cdelt[4];
        
    naxes[0] = galaxy.n_long;
    naxes[1] = galaxy.n_lat;             
    naxes[2] = galaxy.n_nu_synchgrid;
    naxes[3] = 1;
    
    const long nElements = naxes[0]*naxes[1]*naxes[2]*naxes[3];
    
    valarray<float> array(0., nElements);
    valarray<float> arrayQ(0., nElements);//AWS20100709
    valarray<float> arrayU(0., nElements);//AWS20100709

    //float *array;          
    //array=new float[nelements];
  
    int i = 0;
  
    for (int ip = 0; ip < naxes[2]; ++ip) {

      for (int ib = 0; ib < naxes[1]; ++ib) {
	
	for (int il = 0; il < naxes[0]; ++il) {
	  
	  array [i] = galaxy.synchrotron_skymap  .d2[il][ib].s[ip];
	  arrayQ[i] = galaxy.synchrotron_Q_skymap.d2[il][ib].s[ip];//AWS20100709
	  arrayU[i] = galaxy.synchrotron_U_skymap.d2[il][ib].s[ip];//AWS20100709

	  ++i;
	  
	}
	
      }
      
    }

    if (galdef.skymap_format != 1) {

      crval[0] = galaxy.long_min;
      crval[1] = galaxy. lat_min;
      crval[2] = log10(galaxy.nu_synch_min);
      crval[3] = 1;

      cdelt[0] = galaxy.d_long;
      cdelt[1] = galaxy.d_lat;
      cdelt[2] = log10(galaxy.nu_synch_factor);
      cdelt[3] = 1;

      //Use the standard method to store the skymap
      status = store_skymap(&array [0], naxes, "synchrotron_skymap_",   crval, cdelt);
      status = store_skymap(&arrayQ[0], naxes, "synchrotron_Q_skymap_", crval, cdelt);//AWS20100709
      status = store_skymap(&arrayU[0], naxes, "synchrotron_U_skymap_", crval, cdelt);//AWS20100709
    }

    if ((1 == galdef.skymap_format) || 
	(2 == galdef.skymap_format)) { //!< Mapcube output compatable with Glast science tools

      status = store_mapcube_skymap(&array [0], galaxy.nu_synch, 1, galaxy.n_nu_synchgrid, "synchrotron_mapcube_",   false);
      status = store_mapcube_skymap(&arrayQ[0], galaxy.nu_synch, 1, galaxy.n_nu_synchgrid, "synchrotron_Q_mapcube_", false);
      status = store_mapcube_skymap(&arrayU[0], galaxy.nu_synch, 1, galaxy.n_nu_synchgrid, "synchrotron_U_mapcube_", false);

    }

  }

  INFO("Exit");

  return status;

}
