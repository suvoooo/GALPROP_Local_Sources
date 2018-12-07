//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_IC_skymap.cc *                          galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

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

int Galprop::store_IC_skymap(const std::string& IC_type) { //AWS20090415
 
  INFO("Exit");

  int status = 0;

  if (3 == galdef.skymap_format) {

    string filename = configure.fOutputDirectory + configure.fOutputPrefix;
    
    if (IC_type== "isotropic")  filename += "ics_isotropic_";   //AWS20090415
    if (IC_type=="anisotropic") filename += "ics_anisotropic_"; //AWS20090415
    filename += "healpix_";
    filename += galdef.galdef_ID;
    filename += ".gz";
    if (IC_type== "isotropic"){                             //AWS20090415
      //Add all the maps to an initial map
      Skymap<double> outmap(galaxy.IC_iso_hp_skymap[0]);
      for (int i_comp=1; i_comp<galaxy.n_ISRF_components; i_comp++){
	outmap += galaxy.IC_iso_hp_skymap[i_comp];
      }
      SkymapToFits(outmap, filename, "Energy", "MeV");
    }

    if(       IC_type=="anisotropic" ){                               //AWS20090415
      //Add all the maps to an initial map
      Skymap<double> outmap(galaxy.IC_aniso_hp_skymap[0]);
      for (int i_comp=1; i_comp<galaxy.n_ISRF_components; i_comp++)
	outmap += galaxy.IC_aniso_hp_skymap[i_comp];
      SkymapToFits(outmap, filename, "Energy", "MeV");
    }

  } else {
    //fitsfile *fptr;       // pointer to the FITS file; defined in fitsio.h
    
    long  nelements;
    long naxes[4];
    string filename;
    double crval[4], cdelt[4];
    
    naxes[0]=galaxy.n_long;
    naxes[1]=galaxy.n_lat;             
    naxes[2]=galaxy.n_E_gammagrid;
    naxes[3]=1;
    nelements=naxes[0]*naxes[1]*naxes[2]*naxes[3];
    
    valarray<float> array(0., nelements);

    if(       IC_type!= "isotropic"  && IC_type!="anisotropic")   //AWS20090415
      { cout<<" store_IC_skymap: invalid IC_type "<<IC_type<<endl; return 1; }
    if(       IC_type==  "isotropic"    ) filename = "ics_isotropic_";     //AWS20090415
    if(       IC_type=="anisotropic"    ) filename = "ics_anisotropic_";   //AWS20090415
    
    if (galdef.skymap_format != 1){ //!< Old output
      int i=0; 
      for (int ip=0; ip<naxes[2]; ip++)
	for (int ib=0; ib<naxes[1]; ib++)
	  for (int il=0; il<naxes[0]; il++)
	    {
	      array[i]=0.0;
	      for (int i_comp=0; i_comp<galaxy.n_ISRF_components; i_comp++)
		{ 
		  if(       IC_type==  "isotropic"    ) array[i]+=galaxy.IC_iso_skymap  [i_comp].d2[il][ib].s[ip]; //AWS20090415
		  if(       IC_type=="anisotropic"    ) array[i]+=galaxy.IC_aniso_skymap[i_comp].d2[il][ib].s[ip]; //AWS20090415
		}
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
      status = store_skymap(&array[0], naxes, filename+"skymap_", crval, cdelt);
      
    }
    if(galdef.skymap_format == 1 || galdef.skymap_format == 2){ //!< Mapcube output compatable with Glast science tools
      int i=0; 
      for (int ip=0; ip<naxes[2]; ip++)
	for (int ib=0; ib<naxes[1]; ib++)
	  for (int il=0; il<naxes[0]; il++)
	    {
	      array[i]=0.0;
	      for (int i_comp=0; i_comp<galaxy.n_ISRF_components; i_comp++)
		{ 
		  if(       IC_type==  "isotropic"    ) array[i]+=galaxy.IC_iso_skymap  [i_comp].d2[il][ib].s[ip]; //AWS20090415
		  if(       IC_type=="anisotropic"    ) array[i]+=galaxy.IC_aniso_skymap[i_comp].d2[il][ib].s[ip]; //AWS20090415
		}
	      i++;
	    }
      
      status = store_mapcube_skymap(&array[0], galaxy.E_gamma, 1, galaxy.n_E_gammagrid, filename+"mapcube_", true);
    }

  }

  INFO("Exit");

  return status;

}
