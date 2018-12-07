//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_bremss_emiss.cc *                       galprop package * 4/14/2000 
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

int Galprop::store_bremss_emiss() {

  INFO("Entry");

  int stat;
  
  stat=0;
  fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
  int status, ii, jj;
  long  fpixel = 1, naxis = 4, nelements, exposure;
  long naxes[4]  ; 
  
  if(galaxy.n_spatial_dimensions==2) naxes[0]=galaxy.n_rgrid;
  if(galaxy.n_spatial_dimensions==3) naxes[0]=galaxy.n_xgrid;
  naxes[1]=galaxy.n_zgrid;             
  naxes[2]=galaxy.n_E_gammagrid;
  naxes[3]=1;
  
  //cout<<galaxy.n_spatial_dimensions<<" "<<naxes[0]<<" "<<naxes[1]<<"  "<<naxes[2]<<endl;
  
  nelements=naxes[0]*naxes[1]*naxes[2]*naxes[3];
  
  valarray<float> array(0., nelements);
  
  //char outfile[100];
  //strcpy(outfile,"!"); /* create new file or overwrite existing one */
  //strcat(outfile,configure.fFITSDataDirectory);
  //strcat(outfile,"bremss_emiss_");
  //strcat(outfile,galdef.galdef_ID);
  //strcat(outfile,".gz");
  const std::string outfile = "!" + configure.fOutputDirectory + configure.fOutputPrefix + "bremss_emiss_" + galdef.galdef_ID + ".gz";
  //cout<<" storing bremss_emiss in file "<<outfile<<endl;
  
  status = 0;         /* initialize status before calling fitsio routines */
  fits_create_file(&fptr, outfile.c_str(), &status);   /* create new file or overwrite existing one */
  
  /* Create the primary array image (32-bit float pixels */
  fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);
  
  /* Write a keyword; must pass the ADDRESS of the value */
  exposure = 1500;
  fits_update_key(fptr, TLONG, "EXPOSURE", &exposure, "Total Exposure Time", &status);
 
  // for 3D case store x-dimension at y=0
  
  int i=0;
  
  for (int ip       =0;        ip<naxes[2];       ip++)
    {
      for (int iz       =0;        iz<naxes[1];       iz++)
	{
	  for (int ir       =0;        ir<naxes[0];       ir++)
	    {
	      if(galaxy.n_spatial_dimensions==2)array[i]=galaxy.bremss_emiss.d2[ir]   [iz].s[ip];
	      if(galaxy.n_spatial_dimensions==3)array[i]=galaxy.bremss_emiss.d3[ir][0][iz].s[ip];
	      array[i]*=pow(galaxy.E_gamma[ip],2);
	      i++;
	    }
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
  crval3=log10(galaxy.E_gamma_min);
  crval4=1;
  
  if(galaxy.n_spatial_dimensions==2)cdelt1=galaxy.dr;
  if(galaxy.n_spatial_dimensions==3)cdelt1=galaxy.dx;
  cdelt2=galaxy.dz;
  cdelt3=log10(galaxy.E_gamma_factor);
  cdelt4=1;
  
  fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1,"Start of axis 1", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2,"Start of axis 2", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL3", &crval3,"Start of axis 3", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL4", &crval4,"Start of axis 4", &status);
  
  fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1,"Increment of axis 1", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2,"Increment of axis 2", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT3", &cdelt3,"Increment of axis 3", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT4", &cdelt4,"Increment of axis 4", &status);
  
  
  /*
  // write keywords describing nuclei
  char keyword[20];
  char comment[40];
  
  for(int i_nucleus=0;i_nucleus<n_species;i_nucleus++){
  sprintf(keyword,"NUCZ%03d",       i_nucleus+1 ); // e.g. NUCZ012
  sprintf(comment,"Z of nucleus %d",i_nucleus+1 ); 
  cout<<keyword<<" "<<gcr[i_nucleus].Z<<endl;
  fits_update_key(fptr, TINT   , keyword , &gcr[i_nucleus].Z,comment, &status);
  
  sprintf(keyword,"NUCA%03d",       i_nucleus+1 ); // e.g. NUCA012
  sprintf(comment,"A of nucleus %d",i_nucleus+1 ); 
  cout<<keyword<<" "<<gcr[i_nucleus].A<<endl;
  fits_update_key(fptr, TINT   , keyword , &gcr[i_nucleus].A,comment, &status);
  
  
  }
  */
  
  fits_close_file(fptr, &status);            /* close the file */
  
  fits_report_error(stderr, status);  /* print out any error messages */
  //delete[] array;       //Gulli20070810
  //return( status );
  
  //cout<<" <<<< store_bremss_emiss"<<endl;
  //return stat;

  INFO("Exit");

  return status;

}
