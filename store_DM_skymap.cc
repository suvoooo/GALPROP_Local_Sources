
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_DM_skymap.cc *                           galprop package * 9/14/2005 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

//  IMOS20050912

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

int Galprop::store_DM_skymap() {
  
  INFO("Entry");

  int status = 0;

  if (3 == galdef.skymap_format) {

    string filename = configure.fOutputDirectory + configure.fOutputPrefix;
    filename += "DM_healpix_";
    filename += galdef.galdef_ID;
    filename += ".gz";
    SkymapToFits(galaxy.DM_hp_skymap, filename, "Energy", "MeV");
  
  }else{

    //fitsfile *fptr;       // pointer to the FITS file; defined in fitsio.h
    int  status, ii, jj;
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
      for (int ip       =0;        ip<naxes[2];       ip++)
	for (int ib       =0;        ib<naxes[1];       ib++)
	  for (int il       =0;        il<naxes[0];       il++)
	    {
	      array[i]=0.0;
	      array[i]+=galaxy.DM_skymap .d2[il][ib].s[ip];
	      array[i]*=pow(galaxy.E_gamma[ip],2);
	      i++;
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
      status = store_skymap(&array[0], naxes, "DM_skymap_", crval, cdelt);
      
    }

    if(galdef.skymap_format == 1 || galdef.skymap_format == 2){ //!< Mapcube output compatable with Glast science tools
      int i=0; 
      for (int ip       =0;        ip<naxes[2];       ip++)
	for (int ib       =0;        ib<naxes[1];       ib++)
	  for (int il       =0;        il<naxes[0];       il++)
	    {
	      array[i]=galaxy.DM_skymap .d2[il][ib].s[ip];
	      i++;
	    }
      status = store_mapcube_skymap(&array[0], galaxy.E_gamma, 1, galaxy.n_E_gammagrid, "DM_mapcube_", true);
    }
  
  }

  INFO("Exit");

  return status;

}
