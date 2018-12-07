
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_gcr.cc *                                galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cassert>
#include <string>
#include <cstring>
#include <valarray>
#include <sstream>

using namespace std;//AWS20050624

#include "galprop_classes.h"
#include "galprop_internal.h"
#include "fitsio.h"

#include <ErrorLogger.h>

int Galprop::store_gcr() {

  INFO("Entry");

  fitsfile *fptr;                            // pointer to the FITS file; defined in fitsio.h
  int status = 0;
  long fpixel = 1, naxis = 4, nelements, exposure;
  long naxes[4]  ; 
  
  if(gcr[0].n_spatial_dimensions==2) naxes[0]=gcr[0].n_rgrid;
  if(gcr[0].n_spatial_dimensions==3) naxes[0]=gcr[0].n_xgrid;
  naxes[1]=1;                                // z = 0 only
  naxes[2]=gcr[0].n_pgrid;
  naxes[3]=n_species;
  
  nelements=naxes[0]*naxes[1]*naxes[2]*naxes[3];
  
  valarray<float> array(0., nelements);

  const std::string outfile = "!" + configure.fOutputDirectory + configure.fOutputPrefix + "nuclei_" + galdef.galdef_ID + ".gz";
  //cout<<" storing nuclei in file "<<outfile<<endl;

   status = 0;                                // initialize status before calling fitsio routines
   fits_create_file(&fptr, outfile.c_str(), &status); // create new file or overwrite existing one

  // Create the primary array image (32-bit float pixels)
   fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);

  // Write a keyword; must pass the ADDRESS of the value
   exposure = 1500;
   fits_update_key(fptr, TLONG, "EXPOSURE", &exposure,"Total Exposure Time", &status);
    
   int iz1=(int)(1.e-6 -galdef.z_min/galdef.dz);// must correspond to z=0
   ostringstream buf;
   buf<<" iz value for output = "<<iz1<<",   z[iz1] = "<<gcr[0].z[iz1];
   INFO(buf.str());

   if(gcr[0].n_spatial_dimensions==2)
   {
      int i=0;
      for (int i_species=0;i_species< naxes[3];i_species++)
	for (int ip=0;        ip<naxes[2];       ip++)
	  for (int iz=0;        iz<naxes[1];       iz++)
	    for (int ir=0;        ir<naxes[0];       ir++)
	      {
		if(gcr[i_species].A!=0)
		  array[i]=    gcr[i_species].cr_density.d2[ir]    [iz1+iz].s[ip] 
		    *    gcr[i_species].A
		    *pow(gcr[i_species].Ekin[ip],2);
		if(gcr[i_species].A==0)  // electrons, positrons
		  array[i]=    gcr[i_species].cr_density.d2[ir]    [iz1+iz].s[ip] 
		    *pow(gcr[i_species].Ekin[ip],2);
		i++;
	      }
    }
  
  if(gcr[0].n_spatial_dimensions==3)
    {
      int iy=(int)(1.e-6 -galdef.y_min/galdef.dy); // must correspond to y=0
      if(galdef.verbose>=1)cout<<" iy value for output = "<<iy<<",   y[iy] = "<<gcr[0].y[iy]<<endl;
      int i=0;
      for (int i_species=0;i_species< naxes[3];i_species++)
	for (int ip=0;        ip<naxes[2];       ip++)
	  for (int iz=0;        iz<naxes[1];       iz++)
	    for (int ix=0;        ix<naxes[0];       ix++)
	      {
		if(gcr[i_species].A!=0)
		  array[i]=    gcr[i_species].cr_density.d3[ix][iy][iz1+iz].s[ip] 
		    *    gcr[i_species].A
		    *pow(gcr[i_species].Ekin[ip],2);
		if(gcr[i_species].A==0)     // electrons, positrons
		  array[i]=    gcr[i_species].cr_density.d3[ix][iy][iz1+iz].s[ip]
		    *pow(gcr[i_species].Ekin[ip],2);
		i++;
	      }
    }
  
  // Write the array of floats to the image
  fits_write_img(fptr, TFLOAT, fpixel, nelements, &array[0], &status);
  
  // write basic FITS keywords
  double crval1,crval2,crval3,crval4;
  double cdelt1,cdelt2,cdelt3,cdelt4;
  
  if(gcr[0].n_spatial_dimensions==2) crval1=gcr[0].r[0];
  if(gcr[0].n_spatial_dimensions==3) crval1=gcr[0].x[0];
  crval2=gcr[0].z[iz1];
  crval3=log10(gcr[0].Ekin[0]);
  crval4=1;
  if(gcr[0].n_spatial_dimensions==2) cdelt1=gcr[0].dr;  
  if(gcr[0].n_spatial_dimensions==3) cdelt1=gcr[0].dx;  
  cdelt2=gcr[0].dz;    
  cdelt3=log10(gcr[0].Ekin_factor);
  cdelt4=1;
  
  fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1,"Start of axis 1", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2,"Start of axis 2", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL3", &crval3,"Start of axis 3", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL4", &crval4,"Start of axis 4", &status);
  
  fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1,"Increment of axis 1", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2,"Increment of axis 2", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT3", &cdelt3,"Increment of axis 3", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT4", &cdelt4,"Increment of axis 4", &status);
  
  // write keywords describing nuclei
  char keyword[20], comment[40];
  
  for(int i_nucleus=0;i_nucleus<n_species;i_nucleus++)
    {
      sprintf(keyword,"NUCZ%03d",       i_nucleus+1 ); // e.g. NUCZ012
      sprintf(comment,"Z of nucleus %d",i_nucleus+1 ); 
      buf.str("");
      buf<<keyword<<" "<<gcr[i_nucleus].Z;
      INFO(buf.str());
      fits_update_key(fptr, TINT   , keyword , &gcr[i_nucleus].Z,comment, &status);
      
      sprintf(keyword,"NUCA%03d",       i_nucleus+1 ); // e.g. NUCA012
      sprintf(comment,"A of nucleus %d",i_nucleus+1 ); 
      buf.str("");
      buf<<keyword<<" "<<gcr[i_nucleus].A;
      INFO(buf.str());
      fits_update_key(fptr, TINT   , keyword , &gcr[i_nucleus].A,comment, &status);
      
      sprintf(keyword,"NUCK%03d",       i_nucleus+1 ); // e.g. NUCK012                      //AWS20010731
      sprintf(comment,"K-electrons of nucleus %d",i_nucleus+1 );                            //
      buf.str("");
      buf<<keyword<<" "<<gcr[i_nucleus].K_electron;
      INFO(buf.str());
      fits_update_key(fptr, TINT   , keyword , &gcr[i_nucleus].K_electron,comment, &status);//
    }
  fits_close_file(fptr, &status);     // close the file
  fits_report_error(stderr, status);  // print out any error messages
  
  INFO("Exit");

  return status;

}
