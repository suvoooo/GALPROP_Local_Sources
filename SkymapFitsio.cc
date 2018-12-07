#include "SkymapFitsio.h"
#include "fitsio.h"
#include <sstream>
#include <cmath>
#include <string.h>
#include <valarray>

/*
	int FitsToSkymap (Skymap<double> & skymap, const std::string &filename){
		int status = 0;
		fitsfile *fptr;

		//Open the file
		fits_open_file(&fptr, filename.c_str(), READONLY, &status);
		if (status){
			std::cerr<<"Error opening file '"<<filename<<"'."<<std::endl;
			return status;
		}

		//Get the SKYMAP HDU
		fits_movnam_hdu(fptr, BINARY_TBL, "SKYMAP", 0, &status);
		if (status == BAD_HDU_NUM){
			std::cerr<<"File '"<<filename<<"' not in correct format!"<<std::endl;
			return status;
		}
		//Read in the size of the thing
		int nSide, nRows, nSpectra;
		fits_read_key(fptr, TINT, "NAXIS2", &nRows, NULL, &status);
		fits_read_key(fptr, TINT, "NSIDE", &nSide, NULL, &status);
		fits_read_key(fptr, TINT, "NBRBINS", &nSpectra, NULL, &status);

		std::cout<<"nRows="<<nRows<<std::endl;
		std::cout<<"nSide="<<nSide<<std::endl;
		std::cout<<"nSpectra="<<nSpectra<<std::endl;

		//Read the actual data
		long nElements = nRows*nSpectra;
		double *array = new double[nElements];
		double nulval = 0;
		int anynul;
		fits_read_col(fptr, TDOUBLE, 1, 1L, 1L, nElements, 0, array, &anynul, &status);

		return(status);
	}
*/

	/**\brief Read skymap from a fits file*/
	int FitsToSkymap (Skymap<long> & skymap, const std::string &filename){
		int status = 0;
		int nRows, nSpectra, nElements, nSide;
		fitsfile *fptr;
		char extname[100];
		char comment[100];
		double nulval = 0;
		int anynul = 0;
		char ordering[10];
		
		//Open fits file in readonly mode
		fits_open_file(&fptr, filename.c_str(), READONLY, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		//Get the skymap HDU
		strcpy(extname,"SKYMAP");
		fits_movnam_hdu(fptr, BINARY_TBL, extname, 0, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		//Read specific keywords to tell datatype and length
		fits_read_key(fptr, TINT, "NAXIS2", &nRows, comment, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		fits_read_key(fptr, TINT, "NBRBINS", &nSpectra, comment, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		fits_read_key(fptr, TINT, "NSIDE", &nSide, comment, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		fits_read_key(fptr, TSTRING, "ORDERING", ordering, comment, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		//Read in the skymap data
		nElements = nRows*nSpectra;
		long *array = new long[nElements];
		long frow=1, felement=1;
		int cnum=1;
		fits_read_col(fptr,TLONG, cnum, frow, felement, nElements, &nulval, array, &anynul, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		//Assign the values to the skymap
		std::valarray<long> spectra(nSpectra);
		int order = int(log(nSide)/log(2)+0.1);
		if (! strcmp(ordering,"NEST") ){
			skymap.Resize(order, nSpectra, NEST);
		}else{
			skymap.Resize(order, nSpectra, RING);
		}
		for (int i=0; i<nRows; ++i){
			for (int j=0; j<nSpectra; ++j){
				int l = i*nSpectra+j;
				spectra[j] = array[l];
			}
			skymap[i] = spectra;
		}

		delete[] array;

		//Get the energies HDU
		strcpy(extname, "ENERGIES");
		fits_movnam_hdu(fptr, BINARY_TBL, extname, 0, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		//Read the spectra
		double *arr = new double[nSpectra];
		fits_read_col(fptr,TDOUBLE, 1, 1, 1, nSpectra, &nulval, arr, &anynul, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		skymap.setSpectra(arr, nSpectra);
		delete[] arr;
		fits_close_file(fptr, &status);
		if(status) fits_report_error(stderr, status);

		return status;
	}//FitsToSkymap

	int FitsToSkymap (Skymap<double> & skymap, const std::string &filename){
		int status = 0;
		int nRows, nSpectra, nElements, nSide;
		fitsfile *fptr;
		char extname[100];
		char comment[100];
		double nulval = 0;
		int anynul = 0;
		char ordering[10];
		
		//Open fits file in readonly mode
		fits_open_file(&fptr, filename.c_str(), READONLY, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		//Get the skymap HDU
		strcpy(extname,"SKYMAP");
		fits_movnam_hdu(fptr, BINARY_TBL, extname, 0, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		//Read specific keywords to tell datatype and length
		fits_read_key(fptr, TINT, "NAXIS2", &nRows, comment, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		fits_read_key(fptr, TINT, "NBRBINS", &nSpectra, comment, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		fits_read_key(fptr, TINT, "NSIDE", &nSide, comment, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		fits_read_key(fptr, TSTRING, "ORDERING", ordering, comment, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		//Read in the skymap data
		nElements = nRows*nSpectra;
		double *array = new double[nElements];
		long frow=1, felement=1;
		int cnum=1;
		fits_read_col(fptr,TDOUBLE, cnum, frow, felement, nElements, &nulval, array, &anynul, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		//Assign the values to the skymap
		std::valarray<double> spectra(nSpectra);
		int order = int(log(nSide)/log(2)+0.1);
		if (! strcmp(ordering,"NEST") ){
			skymap.Resize(order, nSpectra, NEST);
		}else{
			skymap.Resize(order, nSpectra, RING);
		}
		for (int i=0; i<nRows; ++i){
			for (int j=0; j<nSpectra; ++j){
				int l = i*nSpectra+j;
				spectra[j] = array[l];
			}
			skymap[i] = spectra;
		}

		delete[] array;

		//Get the energies HDU
		strcpy(extname, "ENERGIES");
		fits_movnam_hdu(fptr, BINARY_TBL, extname, 0, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		//Read the spectra
		array = new double[nSpectra];
		fits_read_col(fptr,TDOUBLE, 1, 1, 1, nSpectra, &nulval, array, &anynul, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		skymap.setSpectra(array, nSpectra);
		delete[] array;
		fits_close_file(fptr, &status);
		if(status) fits_report_error(stderr, status);

		return status;
	}//FitsToSkymap

	/**\brief Output skymap to a fits file*/
	int SkymapToFits (const Skymap<double> & skymap, const std::string & filename, const std::string & sUnit, const std::string & sType){
		int status = 0;
		int keyvalue;
		double keyv;
		unsigned int nelements;
		char *ttypes[] = {new char[15]};
		char *tunits[] = {new char[15]};
		char *tforms[] = {new char[10]};
		char keyval[30];
		char extname[100];
		strcpy(extname,"SKYMAP");
		fitsfile *fptr;
		
		//Prepend ! so the file gets overwritten
		std::string fname = "!" + filename;
		
		//Create the file and an empty header unit
		fits_create_file(&fptr, fname.c_str(), &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		fits_create_img(fptr, FLOAT_IMG, 0, NULL, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		//Create the data table to store the actual image
		std::ostringstream strForm;
		strForm << skymap.nSpectra() << "D";
		strcpy(tforms[0], strForm.str().c_str());
		strcpy(ttypes[0], "Spectra");
		strcpy(tunits[0], "Intensity");
		fits_create_tbl(fptr, BINARY_TBL, 0, 1, ttypes, tforms, tunits, extname, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		//Convert the data to an array
		double *array;
		array = skymap.IOarray( nelements );

		//Write that array to the correct column
		fits_write_col(fptr, TDOUBLE, 1, 1, 1, nelements, array, &status);
		if (status){
			std::cout<<"Error in Skymap output"<<std::endl;
			fits_report_error(stderr, status);
			return status;
		}

		delete[] array;
		
		//Get the spectra from the skymap
		std::valarray<double> spectra = skymap.getSpectra();

		//Write keywords
		//Healpix
		strcpy(keyval,"HEALPIX");
		fits_write_key(fptr, TSTRING, "PIXTYPE", keyval, "", &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		if (NEST == skymap.Scheme()){
			strcpy(keyval,"NEST");
		}else{
			strcpy(keyval,"RING");
		}
		fits_write_key(fptr, TSTRING, "ORDERING", keyval, "", &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		keyvalue = skymap.Nside();
		fits_write_key(fptr, TINT, "NSIDE", &keyvalue, "", &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		keyvalue = 0;
		fits_write_key(fptr, TINT, "FIRSTPIX", &keyvalue, "", &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		keyvalue = skymap.Npix()-1;
		fits_write_key(fptr, TINT, "LASTPIX", &keyvalue, "", &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		//Energy binning
		keyvalue = skymap.nSpectra();
		fits_write_key(fptr, TINT, "NBRBINS", &keyvalue, "Number of energy bins", &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		keyv = spectra[0];
		fits_write_key(fptr, TDOUBLE, "EMIN", &keyv, "Minimum energy", &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		keyv = 0;
		if (spectra.size() > 1){
			keyv = log(spectra[1]/spectra[0]);
		}
		fits_write_key(fptr, TDOUBLE, "DELTAE", &keyv, "Step in energy (log)", &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		//Create another table to store the energy
		strcpy(tforms[0], "1D");
		strcpy(ttypes[0], sType.c_str());
		strcpy(tunits[0], sUnit.c_str());
		strcpy(extname, "ENERGIES");

		fits_create_tbl(fptr, BINARY_TBL, 0, 1, ttypes, tforms, tunits, extname, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		fits_write_col(fptr, TDOUBLE, 1, 1, 1, skymap.nSpectra(), &spectra[0], &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		fits_close_file(fptr, &status);
		fits_report_error(stderr, status);
		return status;
	}//int SkymapToFits

	int SkymapToFits (const Skymap<long> & skymap, const std::string & filename, const std::string & sUnit, const std::string & sType){
		int status = 0;
		int keyvalue;
		double keyv;
		unsigned int nelements;
		char *ttypes[] = {new char[15]};
		char *tunits[] = {new char[15]};
		char *tforms[] = {new char[10]};
		char keyval[30];
		char extname[100];
		strcpy(extname,"SKYMAP");
		fitsfile *fptr;
		
		//Prepend ! so the file gets overwritten
		std::string fname = "!" + filename;
		
		//Create the file and an empty header unit
		fits_create_file(&fptr, fname.c_str(), &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		fits_create_img(fptr, FLOAT_IMG, 0, NULL, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		//Create the data table to store the actual image
		std::ostringstream strForm;
		strForm << skymap.nSpectra() << "K";
		strcpy(tforms[0], strForm.str().c_str());
		strcpy(ttypes[0], "Spectra");
		strcpy(tunits[0], "Intensity");
		fits_create_tbl(fptr, BINARY_TBL, 0, 1, ttypes, tforms, tunits, extname, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		//Convert the data to an array
		long *array;
		array = skymap.IOarray( nelements );
		std::cout<<"Nelements:"<<nelements<<", Actual number of elements:"<<skymap.Npix()*skymap.nSpectra()<<std::endl;

		//Write that array to the correct column
		fits_write_col(fptr, TLONGLONG, 1, 1, 1, nelements, array, &status);
		if (status){
			std::cout<<"Error in Skymap output"<<std::endl;
			fits_report_error(stderr, status);
			return status;
		}

		delete[] array;
		
		//Get the spectra from the skymap
		std::valarray<double> spectra = skymap.getSpectra();

		//Write keywords
		//Healpix
		strcpy(keyval,"HEALPIX");
		fits_write_key(fptr, TSTRING, "PIXTYPE", keyval, "", &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		if (NEST == skymap.Scheme()){
			strcpy(keyval,"NEST");
		}else{
			strcpy(keyval,"RING");
		}
		fits_write_key(fptr, TSTRING, "ORDERING", keyval, "", &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		keyvalue = skymap.Nside();
		fits_write_key(fptr, TINT, "NSIDE", &keyvalue, "", &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		keyvalue = 0;
		fits_write_key(fptr, TINT, "FIRSTPIX", &keyvalue, "", &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		keyvalue = skymap.Npix()-1;
		fits_write_key(fptr, TINT, "LASTPIX", &keyvalue, "", &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		//Energy binning
		keyvalue = skymap.nSpectra();
		fits_write_key(fptr, TINT, "NBRBINS", &keyvalue, "Number of energy bins", &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		keyv = spectra[0];
		fits_write_key(fptr, TDOUBLE, "EMIN", &keyv, "Minimum energy", &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		keyv = 0;
		if (spectra.size() > 1){
			keyv = log(spectra[1]/spectra[0]);
		}
		fits_write_key(fptr, TDOUBLE, "DELTAE", &keyv, "Step in energy (log)", &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		//Create another table to store the energy
		strcpy(tforms[0], "1D");
		strcpy(ttypes[0], sType.c_str());
		strcpy(tunits[0], sUnit.c_str());
		strcpy(extname, "ENERGIES");

		fits_create_tbl(fptr, BINARY_TBL, 0, 1, ttypes, tforms, tunits, extname, &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}

		fits_write_col(fptr, TDOUBLE, 1, 1, 1, skymap.nSpectra(), &spectra[0], &status);
		if (status){
			fits_report_error(stderr, status);
			return status;
		}
		fits_close_file(fptr, &status);
		fits_report_error(stderr, status);
		return status;
	}//int SkymapToFits
