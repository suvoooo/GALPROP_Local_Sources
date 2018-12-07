//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_bremss_skymap.cc *                      galprop package * 4/14/2000 
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

int Galprop::store_bremss_skymap() {
  
  //int stat=0;
  //cout<<" >>>> store_bremss_skymap"<<endl;

  INFO("Entry");

  int status = 0;

  if (3 == galdef.skymap_format) {
    string filename = configure.fOutputDirectory + configure.fOutputPrefix;
    filename += "bremss_healpix_";
    filename += galdef.galdef_ID;
    filename += ".gz";
    SkymapToFits(galaxy.bremss_hp_skymap, filename, "Energy", "MeV");
    
  }else{
    
   //fitsfile *fptr;       // pointer to the FITS file; defined in fitsio.h
   int status, ii, jj;
   long nelements;
   long naxes[4]; 
   double crval[4], cdelt[4];

   naxes[0]=galaxy.n_long;
   naxes[1]=galaxy.n_lat;             
   naxes[2]=galaxy.n_E_gammagrid;
   naxes[3]=1;
   nelements=naxes[0]*naxes[1]*naxes[2]*naxes[3];

   valarray<float> array(0., nelements);

   //char outfile[100];
   //strcpy(outfile,"!"); // create new file or overwrite existing one
   //strcat(outfile,configure.fFITSDataDirectory);
   //strcat(outfile,"bremss_skymap_");
   //strcat(outfile,galdef.galdef_ID);
   //cout<<" storing bremss skymap total in file "<<outfile<<endl;

   //status = 0;         // initialize status before calling fitsio routines
   //fits_create_file(&fptr, outfile, &status);   // create new file or overwrite existing one
//cout<<"create_file"<<endl;
   //fits_report_error(stderr, status);  // print out any error messages

//Create the primary array image (32-bit float pixels)
   //fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);
//cout<<"create_img"<<endl;
   //fits_report_error(stderr, status);  // print out any error messages

   if (galdef.skymap_format != 1){ //!< Old output
     int i=0; 
     for (int ip       =0;        ip<naxes[2];       ip++)
       for (int ib       =0;        ib<naxes[1];       ib++)
         for (int il       =0;        il<naxes[0];       il++)
	   {
	     array[i]=0.0;
	     array[i]+=galaxy.bremss_skymap .d2[il][ib].s[ip] *pow(galaxy.E_gamma[ip],2);
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
     status = store_skymap(&array[0], naxes, "bremss_skymap_", crval, cdelt);
   }

   if(galdef.skymap_format == 1 || galdef.skymap_format == 2){ //!< Mapcube output compatable with Glast science tools
     int i=0; 
     for (int ip       =0;        ip<naxes[2];       ip++) // IMOS20080114
       for (int ib       =0;        ib<naxes[1];       ib++)
         for (int il       =0;        il<naxes[0];       il++)
	   {
	     array[i]=galaxy.bremss_skymap .d2[il][ib].s[ip];
	     i++;
	   }
     
     status = store_mapcube_skymap(&array[0], galaxy.E_gamma, 1, galaxy.n_E_gammagrid, "bremss_mapcube_", true);
   }

  }

  INFO("Exit");

  return status;

}
