//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * store_gcr_full.cc *                           galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cassert>
#include <cstring>
#include <string>
#include <sstream>

#include "galprop_classes.h"
#include "galprop_internal.h"

using namespace std;

#include "fitsio.h" 

#include <ErrorLogger.h>

int Galprop::store_gcr_full() {

  INFO("Entry");

  long nAxis = -1, nElements = -1;
  valarray<long> nAxes;

  if (2 == gcr[0].n_spatial_dimensions) {
  
    nAxis = 4;
    nAxes.resize(nAxis);

    nAxes[0] = gcr[0].n_rgrid;
    nAxes[1] = gcr[0].n_zgrid;
    nAxes[2] = gcr[0].n_pgrid;
    nAxes[3] = n_species;
    
    nElements = nAxes[0]*nAxes[1]*nAxes[2]*nAxes[3];
    
  }

  if (3 == gcr[0].n_spatial_dimensions) {
  
    nAxis = 5;
    nAxes.resize(nAxis);
    
    nAxes[0] = gcr[0].n_xgrid;
    nAxes[1] = gcr[0].n_ygrid;
    nAxes[2] = gcr[0].n_zgrid;
    nAxes[3] = gcr[0].n_pgrid;
    nAxes[4] = n_species;
   
    nElements = nAxes[0]*nAxes[1]*nAxes[2]*nAxes[3]*nAxes[4];
  
  }

  assert (nAxis > 0 && nElements > 0 && nAxes.size());

  valarray<float> array(0., nElements);
    
  //float* array = new float[nelements];

  //char outfile[100];
  //strcpy(outfile,"!"); /* create new file or overwrite existing one */
  //strcat(outfile,configure.fFITSDataDirectory);
  //strcat(outfile,"nuclei_full_");
  //strcat(outfile,galdef.galdef_ID);
  //strcat(outfile,".gz");
    
  const std::string outfile = "!" + configure.fOutputDirectory + configure.fOutputPrefix + "nuclei_full_" + galdef.galdef_ID + ".gz";

  ostringstream buf;
  buf << "Storing full nuclei array in " << outfile;
  INFO(buf.str());

  //cout<<" storing full nuclei array in file "<<outfile<<endl;

  int status = 0;

  fitsfile* fptr = 0; /* pointer to the FITS file; defined in fitsio.h */
    
  //int status, ii, jj;
  //long fpixel = 1, exposure;

  fits_create_file(&fptr, outfile.c_str(), &status);   /* create new file or overwrite existing one */

  /* Create the primary array image (32-bit float pixels */
  fits_create_img(fptr, FLOAT_IMG, nAxis, &nAxes[0], &status);

  /* Write a keyword; must pass the ADDRESS of the value */
  long exposure = 1500;
  fits_update_key(fptr, TLONG, "EXPOSURE", &exposure, "Total Exposure Time", &status);

  if (2 == gcr[0].n_spatial_dimensions) {

    int i = 0;

    for (int i_species = 0; i_species < nAxes[3]; ++i_species) {
 
      for (int ip = 0; ip < nAxes[2]; ++ip) {

	for (int iz = 0; iz < nAxes[1]; ++iz) {
	  
	  for (int ir = 0; ir < nAxes[0]; ++ir) {
	    
	    if (0 != gcr[i_species].A) // nuclei
	      array[i] = gcr[i_species].cr_density.d2[ir][iz].s[ip] 
		*gcr[i_species].A
		*pow(gcr[i_species].Ekin[ip], 2.);
	    
	    if (0 == gcr[i_species].A)  // electrons, positrons
	      array[i] = gcr[i_species].cr_density.d2[ir][iz].s[ip] 
		*pow(gcr[i_species].Ekin[ip], 2.);
	    
	    //cout << i_species << " " << ip << " " << iz << " " << ir << " " << i << " " << array[i] << " " << array[i]/pow(gcr[i_species].Ekin[ip], 2.) << endl;
	    ++i;

	  }//ir

	}//iz

      }//ip

    }//i_species

  }//n_spatial_dimensions==2
  
  if (3 == gcr[0].n_spatial_dimensions) {
    
    int i = 0;

    for (int i_species = 0; i_species < nAxes[4]; ++i_species) {
 
      for (int ip = 0; ip < nAxes[3]; ++ip) {

	for (int iz = 0; iz < nAxes[2]; ++iz) {

	  for (int iy = 0; iy < nAxes[1]; ++iy) {

	    for (int ix = 0; ix < nAxes[0]; ++ix) {
	      
	      if (0 != gcr[i_species].A) // nuclei
		array[i] =  gcr[i_species].cr_density.d3[ix][iy][iz].s[ip] 
		  *gcr[i_species].A
		  *pow(gcr[i_species].Ekin[ip], 2.);
	      	      
	      if (0 == gcr[i_species].A) // electrons, positrons
		array[i] = gcr[i_species].cr_density.d3[ix][iy][iz].s[ip]
		  *pow(gcr[i_species].Ekin[ip], 2.);
	      
	      ++i;

	    }//ix

	  }//iy
	
	}//iz
      
      }//ip
    
    }//i_species
  
  }//n_spatial_dimensions==3
  
  long fpixel = 1;

  /* Write the array of floats to the image */
  fits_write_img(fptr, TFLOAT, fpixel, nElements, &array[0], &status);
    
  // write basic FITS keywords
  
  double crval1, crval2, crval3, crval4, crval5;
  double cdelt1, cdelt2, cdelt3, cdelt4, cdelt5;
  
  if (2 == gcr[0].n_spatial_dimensions) {

    crval1 = gcr[0].r[0];
    crval2 = gcr[0].z[0];
    crval3 = log10(gcr[0].Ekin[0]);
    crval4 = 1;
    cdelt1 = gcr[0].dr;  
    cdelt2 = gcr[0].dz;    
    cdelt3 = log10(gcr[0].Ekin_factor);
    cdelt4 = 1;
  
  }

  if (3 == gcr[0].n_spatial_dimensions) {

    crval1 = gcr[0].x[0];
    crval2 = gcr[0].y[0];
    crval3 = gcr[0].z[0];
    crval4 = log10(gcr[0].Ekin[0]);
    crval5 = 1;
    cdelt1 = gcr[0].dx;
    cdelt2 = gcr[0].dy;    
    cdelt3 = gcr[0].dz;    
    cdelt4 = log10(gcr[0].Ekin_factor);
    cdelt5 = 1;
  
  }

  fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1,"Start of axis 1", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2,"Start of axis 2", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL3", &crval3,"Start of axis 3", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL4", &crval4,"Start of axis 4", &status);
  
  if (3 == gcr[0].n_spatial_dimensions)
    fits_update_key(fptr, TDOUBLE, "CRVAL5", &crval5,"Start of axis 5", &status);

  fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1,"Increment of axis 1", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2,"Increment of axis 2", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT3", &cdelt3,"Increment of axis 3", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT4", &cdelt4,"Increment of axis 4", &status);
  
  if (3 == gcr[0].n_spatial_dimensions)
    fits_update_key(fptr, TDOUBLE, "CDELT5", &cdelt5,"Increment of axis 5", &status);
  
  // write keywords describing nuclei
  char keyword[20];
  char comment[40];
  
  for (int i_nucleus = 0; i_nucleus < n_species; ++i_nucleus) {

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
    sprintf(comment,"K-electrons of nucleus %d",i_nucleus+1 );                            //AWS20010731
    buf.str("");
    buf<<keyword<<" "<<gcr[i_nucleus].K_electron;
    INFO(buf.str());
    fits_update_key(fptr, TINT   , keyword , &gcr[i_nucleus].K_electron,comment, &status);//AWS20010731
 
  }

  // write normalization keywords   AWS20010121
  for (int i_nucleus = 0; i_nucleus < n_species; ++i_nucleus) {

    if (gcr[i_nucleus].A==1 && gcr[i_nucleus].Z==1)
      fits_update_key(fptr, TDOUBLE, "NUCNORM", &galdef.source_normalization,"Nucleon  norm. factor", &status); // IMOS20030217
    //         &gcr[i_nucleus].normalization_factor,"Nucleon  norm. factor", &status);
    
    if (gcr[i_nucleus].A==0 && gcr[i_nucleus].Z==-1)
      fits_update_key(fptr, TDOUBLE, "ELENORM", &galdef.electron_source_normalization,"Electron norm. factor", &status); // IMOS20031016
    //         &gcr[i_nucleus].normalization_factor,"Electron norm. factor", &status);
    
  }
  
  fits_close_file(fptr, &status);            /* close the file */
  
  fits_report_error(stderr, status);  /* print out any error messages */
  //delete[] array;
  
  INFO("Exit");
  //cout<<" <<<< store_gcr_full"<<endl;
  return status;

}
