/**\file store_GLAST_skymap.cc
 * \brief store skymaps in a mapcube format compatible with the GLAST science
 * tools
 *
 * Seperate components are output in different filenames.  The component number
 * is appended to the name.
 */
#include"galprop_classes.h"
#include"galprop_internal.h"

#include <ErrorLogger.h>

#include <PhysicalConstants.h>


#include <vector>
#include <string>
#include <sstream>
#include <cstring>

int 
Galprop::store_mapcube_skymap(float* array, 
			      double* energy, 
			      const int nComponents, 
			      const int nEnergies, 
			      const std::string& name, 
			      const bool MeV) const {

  int status = 0;
  fitsfile* fptr;
  
  const int naxis = 3;
  const int fpixel = 1;
  long int naxes[3];
  double crval[3], cdelt[3], crpix[3];
  //char keyname[16];
  //char keyval[30];
  
  naxes[0] = galaxy.n_long;
  naxes[1] = galaxy.n_lat;
  naxes[2] = nEnergies;
  const long int nelements = naxes[0]*naxes[1]*naxes[2];
  
  cdelt[0] = galaxy.d_long;
  cdelt[1] = galaxy.d_lat;
  cdelt[2] = log10(energy[1]/energy[0]);
  
  crval[0] = 180.;//galaxy.long_min+(naxes[0]-1)*cdelt[0]/2.0;
  crval[1] = 0.;//galaxy.lat_min+(naxes[1]-1)*cdelt[1]/2.0;
  crval[2] = log10(energy[0]);
  
  crpix[0] = (crval[0]-galaxy.long_min)/galaxy.d_long+1;//(naxes[0]+1)/2;
  crpix[1] = (crval[1]-galaxy.lat_min)/galaxy.d_lat+1;//(naxes[1]+1)/2;
  crpix[2] = 1;

  std::string outFileBase;
  std::string outFile;
  //Set the output filename.  We add ! in front to delete already existing
  //files
  outFileBase = "!";
  outFileBase += configure.fOutputDirectory;
  outFileBase += configure.fOutputPrefix;
  outFileBase += name;
  if (nComponents > 1){
    outFileBase += "comp_"; 
  }else{
    outFile = outFileBase + galdef.galdef_ID + ".gz";
  }
  
  //Loop all components creating a seperate file for each
  for (int i = 0; i<nComponents; ++i){
    if (nComponents > 1){
      ostringstream buf;
      buf << outFileBase << i;
      outFile = buf.str();
      outFile += "_";
      outFile += galdef.galdef_ID;
      outFile += ".gz";
    }

    // Obtain total number of photons in map. This will be written into the 
    // header so the user will not have to explicitly calculate it. Note
    // map is in units intensity (MeV^-1 cm^-2 s^-1 sr^-1) -> need to 
    // convert this to m^-2 s^-1 starting with the lowest energy in the map
    // -- TAP 11 Jul 2008
    // Use power law interpolation to do the integral, more accurate that way.

    double sum = 0;
    
    const unsigned long iindex = i*nelements;

    for (long ip = 0; ip < naxes[2]-1; ++ip) {

      const unsigned long pindex = ip*naxes[1]*naxes[0] + iindex;

      const unsigned long pindex2 = (ip+1)*naxes[1]*naxes[0] + iindex;

      for (long ib = 0; ib < naxes[1]; ++ib) { // ib

        const unsigned long bindex = ib*naxes[0];

	const double b = galaxy.lat_min + ib*galaxy.d_lat;

	const double bRad = utl::kPi/2 - b*utl::kConvertDegreesToRadians;

	const double bWeight = sin(bRad)*cdelt[1]*utl::kConvertDegreesToRadians;

	for (long il = 0; il < naxes[0]; ++il) { // il
	  
	  const unsigned long index = bindex + il;
	  
	  const double sWeight = bWeight*cdelt[0]*utl::kConvertDegreesToRadians;
	  
	  //Only do power law interpolation if both are positive
	  if (array[index+pindex] > 0 && array[index+pindex2] > 0) {
	     const double gamma = log(array[index+pindex]/array[index+pindex2])/log(energy[ip]/energy[ip+1]);
	     const double f0 = array[index+pindex];
	     if ( fabs(gamma+1) > 1e-10 ) {
	        sum += f0/(1+gamma)*energy[ip]*(pow(energy[ip+1]/energy[ip],gamma+1)-1)*sWeight*1e4; // convert cm^-2 s^-1 -> m^-2 s^-1
	     } else {
		sum += f0*(log(energy[ip+1]/energy[ip]))*sWeight*1e4;
	     }
	  } else {
             sum += 0.5*(array[index+pindex]+array[index+pindex2])*(energy[ip+1]-energy[ip])*sWeight*1e4;
	  }
	}
      }
    }

    ostringstream sBuf;
    sBuf << "MapCube: " << sum;
    INFO(sBuf.str());

    //Create a file
    fits_create_file(&fptr, outFile.c_str(), &status);
    fits_report_error(stderr, status);  // print out any error messages
    
    //Create the primary array image (32-bit float pixels)	
    fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);
    //cout<<"create_img"<<endl;
    fits_report_error(stderr, status);  // print out any error messages
    
    //Write the array of floats to the image
    fits_write_img(fptr, TFLOAT, fpixel, nelements, &array[i*nelements], &status);
    //cout<<"write_img"<<endl;
    fits_report_error(stderr, status);  // print out any error messages
    
    //sprintf(keyname,"FLUX");

    ostringstream fBuf;
    fBuf << "FLUX";
    fits_write_key (fptr, TDOUBLE, (char*)fBuf.str().c_str(), &sum, "photon flux (m^-2 s^-1)",  &status);
    fits_report_error(stderr, status);

    //Write the headers
    for (int ii=0; ii<3; ++ii) {
      //write all keys to new file
      //sprintf(keyname,"CRVAL%d",ii+1);
      ostringstream crBuf;
      crBuf << "CRVAL" << ii+1;
      fits_write_key (fptr, TDOUBLE, (char*)crBuf.str().c_str(), &crval[ii], "updated ref point",  &status);
      fits_report_error(stderr, status);
      
      //sprintf(keyname,"CDELT%d",ii+1);
      ostringstream cdBuf;
      cdBuf << "CDELT" << ii+1;
      fits_write_key (fptr, TDOUBLE, (char*)cdBuf.str().c_str(), &cdelt[ii], "step size",  &status);
      fits_report_error(stderr, status);
      //sprintf(keyname,"CRPIX%d",ii+1);
      ostringstream cpBuf;
      cpBuf << "CRPIX" << ii+1;
      fits_write_key (fptr, TDOUBLE, (char*)cpBuf.str().c_str(), &crpix[ii], "",  &status);
      fits_report_error(stderr, status);
    }//for ii keys
    
    //---add all remaining missing keys---

    INFO("Add remaining keys to header ...");

    ostringstream keyval;
    
    keyval << "GLON-CAR";
    //strcpy(keyval,"GLON-CAR");
    fits_write_key (fptr, TSTRING, "CTYPE1", (void*)keyval.str().c_str(), "",  &status);
    keyval.str("");
    keyval.clear();
    keyval << "GLAT-CAR";
    fits_write_key (fptr, TSTRING, "CTYPE2", (void*)keyval.str().c_str(), "",  &status);
    keyval.str("");
    keyval.clear();
    keyval << "deg";

    fits_write_key (fptr, TSTRING, "CUNIT1", (void*)keyval.str().c_str(), "",  &status);
    fits_write_key (fptr, TSTRING, "CUNIT2", (void*)keyval.str().c_str(), "",  &status);

    //strcpy(keyval,"GLAT-CAR");
    
    //fits_write_key (fptr, TSTRING, "CTYPE2", keyval, "",  &status);
    
    //strcpy(keyval,"deg");
    //fits_write_key (fptr, TSTRING, "CUNIT2", keyval, "",  &status);

    if (MeV){
   
      keyval.str("");
      keyval.clear();
      keyval << "photon energy";
      //strcpy(keyval,"photon energy");
      fits_write_key (fptr, TSTRING, "CTYPE3", (void*)keyval.str().c_str(), "",  &status);
      keyval.str("");
      keyval.clear();
      keyval << "log_MeV";
      //strcpy(keyval,"log_MeV");
      fits_write_key (fptr, TSTRING, "CUNIT3", (void*)keyval.str().c_str(), "",  &status);
    }else{
      keyval.str("");
      keyval.clear();
      keyval << "photon frequency";
      //strcpy(keyval,"photon frequency");
      fits_write_key (fptr, TSTRING, "CTYPE3", (void*)keyval.str().c_str(), "",  &status);
      keyval.str("");
      keyval.clear();
      keyval << "log_Hz";
      //strcpy(keyval,"log_Hz");
      fits_write_key (fptr, TSTRING, "CUNIT3", (void*)keyval.str().c_str(), "",  &status);
    }
    fits_report_error(stderr, status);
    
    //create output binary table
    char* ttype[] = {new char[15]};
    char* tform[] = {new char[10]};
    char* tunit[] = {new char[10]};
    if (MeV){
      strncpy(ttype[0],"Energy",14);
      ttype[0][14] = '\0';
      strncpy(tform[0],"1D",9);
      tform[0][9] = '\0';
      strncpy(tunit[0],"MeV",9);
      tunit[0][9] = '\0';
    }else{
      strncpy(ttype[0],"Frequency",14);
      ttype[0][14] = '\0';
      strncpy(tform[0],"1D",9);
      tform[0][9] = '\0';
      strncpy(tunit[0],"Hz",9);
      tunit[0][9] = '\0';
    }
    char* extname = "ENERGIES";
    fits_create_tbl (fptr, BINARY_TBL, naxes[2], 1, ttype, tform, tunit, extname, &status);
    fits_report_error(stderr, status);
    
    fits_write_col (fptr, TDOUBLE, 1, 1, 1, naxes[2], energy, &status);
    fits_report_error(stderr, status);
    
    fits_close_file(fptr, &status);     // close the file
    fits_report_error(stderr, status);  // print out any error messages

    delete[] ttype[0];
    delete[] tform[0];
    delete[] tunit[0];

  } //Components
  return (status);
}

