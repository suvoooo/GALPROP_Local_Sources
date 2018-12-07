/**\brief Store a skymap in the conventional galprop format of 4 dimensions
 *
 * \param array has the elements in correct order for the fits writer.
 * \param naxes is a 4 dimensional array containing the dimension of the image.
 * \param name is a string representing the name of the output file.
 * \param crval and param \cdelt are arrays for corresponding fits keywords.
 *
 * \author Gudlaugur Johannesson
 */

#include "fitsio.h"
#include "galprop_classes.h"
#include "galprop_internal.h"

#include <string>
#include <cstdio>
#include <cstring>

int Galprop::store_skymap(float array[], 
			  long naxes[4], 
			  const std::string name, 
			  double crval[4], 
			  double cdelt[4]) {

  int status = 0;
  
  fitsfile* fptr = 0;
  const int naxis = 4;
  const int fpixel = 1;
  long nelements = naxes[0]*naxes[1]*naxes[2]*naxes[3];
  std::string outFile;
  
  //Set the output filename.  We add ! in front to delete already existing
  //files
  outFile = "!";
  outFile += configure.fOutputDirectory;
  outFile += configure.fOutputPrefix;
  outFile += name;
  outFile += galdef.galdef_ID;
  outFile += ".gz";
  
  //Create a file
  fits_create_file(&fptr, outFile.c_str(), &status);
  fits_report_error(stderr, status);  // print out any error messages
  
  //Create the primary array image (32-bit float pixels)	
  fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);
  //cout<<"create_img"<<endl;
  fits_report_error(stderr, status);  // print out any error messages
  
  //Write the array of floats to the image
  fits_write_img(fptr, TFLOAT, fpixel, nelements, array, &status);
  //cout<<"write_img"<<endl;
  fits_report_error(stderr, status);  // print out any error messages
  
  // write basic FITS keywords
  //char keyname[6];
  //char comment[50];
  for (int ii = 0; ii < naxis; ++ii){

    ostringstream keyname1, comment1, keyname2, comment2;
    
    keyname1 << "CRVAL" << ii+1;
    comment1 << "Start of axis " << ii+1;
    
    keyname2 << "CDELT" << ii+1;
    comment2 << "Increment of axis " << ii+1;

    fits_update_key(fptr, TDOUBLE, (char*)keyname1.str().c_str(), &crval[ii], (char*)comment1.str().c_str(), &status);
    
    fits_update_key(fptr, TDOUBLE, (char*)keyname2.str().c_str(), &cdelt[ii], (char*)comment2.str().c_str(), &status);

    //sprintf(keyname, "CRVAL%d", ii+1);
    //sprintf(comment, "Start of axis %d", ii+1);
    //fits_update_key(fptr, TDOUBLE, keyname, &crval[ii], comment, &status);
    //sprintf(keyname, "CDELT%d", ii+1);
    //sprintf(comment, "Increment of axis %d", ii+1);
    //fits_update_key(fptr, TDOUBLE, keyname, &cdelt[ii], comment, &status);
  }
  
  fits_report_error(stderr, status);  // print out any error messages
  fits_close_file(fptr, &status);     // close the file
  fits_report_error(stderr, status);  // print out any error messages
  return( status );
}
