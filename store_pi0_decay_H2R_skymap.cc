#include <cassert>
#include <string>
#include <cstring>
#include <valarray>

using namespace std;//AWS20050624

#include "galprop_classes.h"
#include "galprop_internal.h"
#include "fitsio.h"
#include "SkymapFitsio.h"

#include <ErrorLogger.h>

int Galprop::store_pi0_decay_H2R_skymap() {

  INFO("Entry");

  int status = 0;

  if (3 == galdef.skymap_format) {

    char index[3];
    string fileprefix = configure.fOutputDirectory + configure.fOutputPrefix;

    fileprefix += "pi0_decay_H2R_ring_";
    for (int i_Ring=0; i_Ring<galaxy.hpHIR.nSpectra(); i_Ring++){
      sprintf(index, "%d", i_Ring+1);
      string filename = fileprefix + index;
      filename += "_healpix_";
      filename += galdef.galdef_ID;
      filename += ".gz";
      SkymapToFits(galaxy.pi0_decay_H2R_hp_skymap[i_Ring], filename, "Energy", "MeV");
    }

  } else {
    
    //fitsfile *fptr;       // pointer to the FITS file; defined in fitsio.h
    int status, ii, jj;
    long nelements;
    long naxes[4]; 
    double crval[4], cdelt[4];
    
    naxes[0]=galaxy.n_long;
    naxes[1]=galaxy.n_lat;   
    naxes[2]=galaxy.pi0_decay_H2R_skymap.n_zgrid; // number of Galactocentric rings
    naxes[3]=galaxy.n_E_gammagrid;
    
    nelements=naxes[0]*naxes[1]*naxes[2]*naxes[3];
    
    // The line below causes a segfault with icpc/ifort 11.1.056 using output
    // option `0'. Change to declare, then resize and initialise the array, 
    // instead of initialisation at declaration, works.

    valarray<float> array(0., nelements);
    //valarray<float> array;
    //array.resize(nelements);
    //array = 0;

    if (galdef.skymap_format != 1){ //!< Old output
      int i=0; 
      for     (int ip       =0;            ip <naxes[3];        ip++)
	for    (int i_Ring   =0;        i_Ring <naxes[2];    i_Ring++)
	  for   (int ib       =0;            ib <naxes[1];        ib++)
	    for (int il       =0;            il <naxes[0];        il++)
	      {
		//array[i]=0.0;
		array[i] =galaxy.pi0_decay_H2R_skymap.d3[il][ib][i_Ring].s[ip] *pow(galaxy.E_gamma[ip],2);
		i++;
	      }
      
      crval[0]=galaxy.long_min;
      crval[1]=galaxy. lat_min;
      crval[2]=0; //IMOS20080114
      crval[3]=log10(galaxy.E_gamma_min);
      
      cdelt[0]=galaxy.d_long;
      cdelt[1]=galaxy.d_lat;
      cdelt[2]=1;
      cdelt[3]=log10(galaxy.E_gamma_factor);
      
      //Use the standard method to store the skymap

      const string name = "pion_decay_H2R_skymap_";

      status = store_skymap(&array[0], naxes, name, crval, cdelt);
      
    }

    if(galdef.skymap_format == 1 || galdef.skymap_format == 2){ //!< Mapcube output compatable with Glast science tools
      int i=0; 
      for    (int i_Ring   =0;        i_Ring <naxes[2];    i_Ring++)
        for (int ip =0; ip <naxes[3]; ip++) // IMOS20080114
	  for   (int ib       =0;            ib <naxes[1];        ib++)
	    for (int il       =0;            il <naxes[0];        il++)
	      {
		array[i]=galaxy.pi0_decay_H2R_skymap .d3[il][ib][i_Ring].s[ip];
		i++;
	      }
      
      status = store_mapcube_skymap(&array[0], galaxy.E_gamma, galaxy.pi0_decay_H2R_skymap.n_zgrid, galaxy.n_E_gammagrid, "pion_decay_H2R_mapcube_", true);
    }

  }

  INFO("Exit");

  return status;

}
