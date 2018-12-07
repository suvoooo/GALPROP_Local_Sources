//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_bremss_ionized_skymap.cc *              galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cassert>
#include <string>
#include <cstring>
#include <valarray>

using namespace std;//AWS20050624

#include "galprop_classes.h"
#include "galprop_internal.h"
#include "SkymapFitsio.h"

#include "fitsio.h" 

#include <ErrorLogger.h>

int Galprop::store_bremss_ionized_skymap() {

  INFO("Entry");

  int status = 0;

  if (galdef.skymap_format == 3) {
    string filename = configure.fOutputDirectory + configure.fOutputPrefix;
    filename += "bremss_ionized_healpix_";
    filename += galdef.galdef_ID;
    filename += ".gz";
    SkymapToFits(galaxy.bremss_ionized_hp_skymap, filename, "Energy", "MeV");
    
  }else{
    //stat=0;
    //fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
    
    long nelements;
    long naxes[4]; 
    double crval[4], cdelt[4];
        
    naxes[0]=galaxy.n_long;
    naxes[1]=galaxy.n_lat;             
    naxes[2]=galaxy.n_E_gammagrid;
    naxes[3]=1;                   

    nelements=naxes[0]*naxes[1]*naxes[2]*naxes[3];

    valarray<float> array(0., nelements);

    if (galdef.skymap_format != 1){ //!< Old output
      int i=0;
      
      for (int ip       =0;        ip<naxes[2];       ip++){
	for (int ib       =0;        ib<naxes[1];       ib++){
	  for (int il       =0;        il<naxes[0];       il++){
	    array[i]=galaxy.bremss_ionized_skymap.d2[il][ib].s[ip] *pow(galaxy.E_gamma[ip],2);//IMOS20080114*
	    i++;
	  }
	}
      }
      
      crval[0]=galaxy.long_min;
      crval[1]=galaxy. lat_min;
      crval[2]=log10(galaxy.E_gamma_min);
      crval[3]=1;
      
      cdelt[0]=galaxy.d_long;
      cdelt[1]=galaxy.d_lat;
      cdelt[2]=log10(galaxy.E_gamma_factor);
      cdelt[3]=1;
      
      //Use the standard method to store the skymap
      status = store_skymap(&array[0], naxes, "bremss_ionized_skymap_", crval, cdelt);
    }
    if(galdef.skymap_format == 1 || galdef.skymap_format == 2){ //!< Mapcube output compatible with Glast science tools
      int i=0;
      
      for (int ip       =0;        ip<naxes[2];       ip++){
	for (int ib       =0;        ib<naxes[1];       ib++){
	  for (int il       =0;        il<naxes[0];       il++){
	    array[i]=galaxy.bremss_ionized_skymap.d2[il][ib].s[ip];//IMOS20080114*
	    i++;
	  }
	}
      }
      status = store_mapcube_skymap(&array[0], galaxy.E_gamma, 1, galaxy.n_E_gammagrid, "bremss_ionized_mapcube_", true);
    }
  
  }

  INFO("Exit");

  return status;

}
