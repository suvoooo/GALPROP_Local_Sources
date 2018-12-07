

#include <cassert>
#include <cstring>
#include <string>
#include <sstream>

#include "galprop_classes.h"
#include "galprop_internal.h"

using namespace std;

#include "fitsio.h" 

#include <ErrorLogger.h>

int Galprop::store_gcr_source_functions(Particle &particle) {

  INFO("Entry");

  int status = 0;

  ostringstream buf;
  long nAxis = -1, nElements = -1;
  valarray<long> nAxes;

  if (2 == gcr[0].n_spatial_dimensions) {
  
    nAxis = 3;
    nAxes.resize(nAxis);

    nAxes[0] = gcr[0].n_rgrid;
    nAxes[1] = gcr[0].n_zgrid;
    nAxes[2] = gcr[0].n_pgrid;
    
    
    nElements = nAxes[0]*nAxes[1]*nAxes[2];
    
  }

  if (3 == gcr[0].n_spatial_dimensions) {
  
    nAxis = 4;
    nAxes.resize(nAxis);
    
    nAxes[0] = gcr[0].n_xgrid;
    nAxes[1] = gcr[0].n_ygrid;
    nAxes[2] = gcr[0].n_zgrid;
    nAxes[3] = gcr[0].n_pgrid;
   
   
    nElements = nAxes[0]*nAxes[1]*nAxes[2]*nAxes[3];
  
  }

  assert (nAxis > 0 && nElements > 0 && nAxes.size());

  valarray<float> array(0., nElements);


  string origin="origin_not_assigned";    
  string filename,particle_name;
  particle_name=particle.name;

  buf<<"store_gcr_source_functions: particle="<<particle_name<<endl;
  

  filename="filename_not_assigned";

  if(particle_name=="Hydrogen_1"       )  { filename="primary_protons_source_function";     origin="primary";  }
  if(particle_name=="Helium_4"         )  { filename="primary_Helium_source_function";      origin="primary";  }
  if(particle_name=="primary_electrons")  { filename="primary_electrons_source_function";   origin="primary";  }

  if(particle_name=="secondary_electrons"){ filename="secondary_electrons_source_function"; origin="secondary";}
  if(particle_name=="secondary_positrons"){ filename="secondary_positrons_source_function"; origin="secondary";}

  buf<<" "<<origin<<" "<<filename<<endl;

  if(origin=="origin_not_assigned"||filename=="filename_not_assigned") {  INFO(buf.str());INFO("Exit");status=1; return status;}

  const std::string outfile = "!" + configure.fOutputDirectory + configure.fOutputPrefix + filename + "_"+ galdef.galdef_ID + ".gz";


  buf << "Storing " <<particle_name   <<" as "<<origin<<" source function array in " << outfile;
  INFO(buf.str());

  
  





  fitsfile* fptr = 0; /* pointer to the FITS file; defined in fitsio.h */
    
 

  fits_create_file(&fptr, outfile.c_str(), &status);   /* create new file or overwrite existing one */

  /* Create the primary array image (32-bit float pixels */
  fits_create_img(fptr, FLOAT_IMG, nAxis, &nAxes[0], &status);

 

  if (2 == gcr[0].n_spatial_dimensions) {

    int i = 0;

  
 
      for (int ip = 0; ip < nAxes[2]; ++ip) {

	for (int iz = 0; iz < nAxes[1]; ++iz) {
	  
	  for (int ir = 0; ir < nAxes[0]; ++ir) {
	    
	
	    if(origin=="primary"  )array[i] = particle.  primary_source_function.d2[ir][iz].s[ip] ;
	    if(origin=="secondary")array[i] = particle.secondary_source_function.d2[ir][iz].s[ip] ;	       
	    
	    //	    cout << "ip= " << ip << " " << iz << " " << ir << " " << i << " " << array[i]  << endl;
	    ++i;

	  }//ir

	}//iz

      }//ip

   

  }//n_spatial_dimensions==2
  
  if (3 == gcr[0].n_spatial_dimensions) {
    
    int i = 0;

    
 
      for (int ip = 0; ip < nAxes[3]; ++ip) {

	for (int iz = 0; iz < nAxes[2]; ++iz) {

	  for (int iy = 0; iy < nAxes[1]; ++iy) {

	    for (int ix = 0; ix < nAxes[0]; ++ix) {
	      
	   
		 if(origin=="primary"  )array[i] =   particle.primary_source_function  .d3[ix][iy][iz].s[ip]; 
		 if(origin=="secondary")array[i] =   particle.secondary_source_function.d3[ix][iy][iz].s[ip]; 
	      
	      ++i;

	    }//ix

	  }//iy
	
	}//iz
      
      }//ip
    
    
  
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
    
  
  }

  fits_update_key(fptr, TDOUBLE, "CRVAL1", &crval1,"Start of axis 1", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL2", &crval2,"Start of axis 2", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL3", &crval3,"Start of axis 3", &status);
  fits_update_key(fptr, TDOUBLE, "CRVAL4", &crval4,"Start of axis 4", &status);
  


  fits_update_key(fptr, TDOUBLE, "CDELT1", &cdelt1,"Increment of axis 1", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT2", &cdelt2,"Increment of axis 2", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT3", &cdelt3,"Increment of axis 3", &status);
  fits_update_key(fptr, TDOUBLE, "CDELT4", &cdelt4,"Increment of axis 4", &status);
  

  
  // write keywords describing dataset
  char keyword[20];
  char comment[40];
  
  // To be done
  
  fits_close_file(fptr, &status);            /* close the file */
  
  fits_report_error(stderr, status);  /* print out any error messages */

  //  delete[] array;
  
  INFO("Exit");
  
  return status;

}
