/**\class Skymap
 * \brief Skymaps for gamma ray calculation in galprop
 *
 * The skymap class is used to store skymaps calculated in gardian
 * The information is stored in a healpix format.
 * Internally we use a std::valarray to store the spectra for each pixel
 * \author Gudlaugur Johannesson
 */

#ifndef SKYMAP_H
#define SKYMAP_H

#include "HealpixBaseExtended.h"
#include "Coordinate.h"
#include "PhysicalConstants.h"
#include "ArraySlice.h"

#include <healpix_map.h>
#include <arr.h>

#include <CCfits/CCfits>

#include <fitsio.h>

#include <iostream>
#include <vector>
#include <valarray>
#include <map>
#include <iterator>
#include <algorithm>
#include <functional>
#include <cmath>
#include <typeinfo>

//Coordinate transformation if we find libastro
#ifdef HAVE_ASTRO
extern "C"{
#include <astro.h>
}
//Do to stupid convention in astro, we must create wrapper functions
void eq_gal2(double mj, double DEC, double RA, double *b, double *l);
void eq_ecl2(double mj, double DEC, double RA, double *lt, double *lg);
void gal_eq2(double mj, double b, double l, double *DEC, double *RA);
void ecl_eq2(double mj, double lt, double lg, double *DEC, double *RA);
//Functions to convert from gal to ecl and wise versa
void gal_ecl2(double mj, double b, double l, double *lt, double *lg);
void ecl_gal2(double mj, double lt, double lg, double *b, double *l);
void empty(double mj, double x, double y, double *z, double *w);
#endif



/**\brief Store a full sky map in a healpix pixelation.
 *
 * This class is templated to allow for a flexible storage.  Inherits from Healpix_Base so all its methods
 * are available.  Can only be initialized with a whole number order, i.e. nside = 2^order.
 */
template <typename T>
class Skymap : public HealpixBaseExtended {

 public:
  enum CoordSys {EQ, GAL, ECL};
 private:
  std::valarray<T> fMap; //!< A valarray to store the values, we use slices to get single spectra
  std::valarray<double> fSpectra; //!< The values at which the spectra is evaluated, usually energy for photon maps.
  std::valarray<double> fSpecMin, fSpecMax;  //!< If the skymap is binned, we store maximum and minimum as well as mid value for compatibility.
  bool fBinned; //!< Tells if the skymap is spectrally binned or instantaneous.
  
		/** \brief Returns true if the maps have the same size and same fSpectra */
  bool equiv(const Skymap<T> & otherMap) const{
    //First we check for binning equivalance
    bool eq = (fBinned == otherMap.fBinned);
    bool spectra = (fSpectra.size() == otherMap.fSpectra.size());
    if (eq && spectra) {
      //Check for equivalence of spectra
      if (fBinned) {
	for (int i = 0; i < fSpectra.size(); ++i){
	  spectra &= (fabs((fSpecMin[i] - otherMap.fSpecMin[i])/fSpecMin[i]) < 1e-6);
	  spectra &= (fabs((fSpecMax[i] - otherMap.fSpecMax[i])/fSpecMax[i]) < 1e-6);
	}
      }else{
	for (int i = 0; i < fSpectra.size(); ++i){
	  spectra &= (fabs((fSpectra[i] - otherMap.fSpectra[i])/fSpectra[i]) < 1e-6);
	}
      }
    }
    return (scheme_ == otherMap.scheme_ && order_ == otherMap.order_ && spectra && eq);
  }
  
  /**\brief Convert a CAR map to healpix format
   *
   * Given arrays for CRVAL, CDELT and CRPIX (at the center of the pixel),
   * turn the mapcube into pixels.  The skymap has to be resized properly
   * before entering this method.  The \a image is assumed to be three
   * dimensional in the standard FITS way, longitude evolving fastest, then
   * latitude and last spectra.
   *
   * \param image is the valarray of image data, as from ccfits
   * \param crval is the CRVAL array from the fits header
   * \param crpix is the CRPIX array from the fits header
   * \param cdelt is the CDELT array from the fits header
   */
  void fillMapcube(const std::valarray<T> &image, 
		   const std::valarray<long> &axis, 
		   const std::valarray<double> &crval, 
		   const std::valarray<double> &cdelt, 
		   const std::valarray<double> &crpix, 
		   bool correctSA = false){
    //Divide each HEALPix pixel into 256 smaller ones and fill those with the
    //value from the map pixel directly beneeth its center.  The HEALPix
    //pixel value is the average of those 256 smaller ones.
    const int dOrder = std::min(5,13-order_); //Order can not exceed 13
    const int nPix = (1<<dOrder)*(1<<dOrder);
    //Create a nested healpix base for us to work with
    const Healpix_Base finer(order_+dOrder, NEST);
    //We must multiply with the solid angle if requested, since the method
    //does not work unless the value of the pixel is independent of its
    //solid angle
    double hpSolidAngle = 1;
    if (correctSA) {
      hpSolidAngle = solidAngle();
    }
#pragma omp parallel for default(shared) schedule(static)
    for (int p = 0; p<npix_; ++p){
      //Use pixel numbering from the nested scheme
      const int pp = (scheme_ == NEST) ? p : ring2nest(p);
      //To store the total value for all the pixels
      std::valarray<T> totSpectra(T(0), fSpectra.size());
      //Loop over all of the subpixels
      pointing pnt1 = pix2ang(p);
      for (int sp = pp*nPix; sp<(pp+1)*nPix; ++sp){
	const pointing pnt = finer.pix2ang(sp);
	
	//Find the correct index in the array, first the longitude which can
	//loop
	double dil = (pnt.phi*180/utl::kPi-crval[0])/cdelt[0] + crpix[0] + 0.5;
	//If the pixel value is negative, we must loop it over the 360
	//degrees
	if (dil < 1) { 
	  dil += fabs(360./cdelt[0]); 
	  //Pixel values above the boundaries could also indicate looping
	} else if (dil > axis[0]) {
	  dil -= fabs(360./cdelt[0]);
	}
	int il = int(dil - 1); //Since the pixels are 1 based and the array is 0 based
	
	//Then the latitude, which can't loop
	int ib = int( ( (90 - pnt.theta*180/utl::kPi)-crval[1])/cdelt[1] + crpix[1] + 0.5) - 1;
	
	
	//They must both be within bounds to do something useful
	if (ib >= 0 && ib < axis[1] && il >= 0 && il < axis[0]) {
	  const int ind1 = il + ib*axis[0];
	  //We must divide with the solid angle if requested, since the
	  //method does not work unless the value of the pixel is independent
	  //of the solid angle of the pixel
	  double lbSolidAngle = 1;
	  if (correctSA) {
	    double bMiddle = crval[1] + (ib+1 - crpix[1])*cdelt[1];
	    //Note that sin(x-pi/2) = -cos(x)
	    lbSolidAngle = cdelt[0]*utl::kPi/180.*(cos((bMiddle-cdelt[1]/2.)*utl::kPi/180.) - cos((bMiddle+cdelt[1]/2.)*utl::kPi/180.));
	  }
	  for (int is=0; is<fSpectra.size(); ++is){
	    totSpectra[is] += T(image[ind1+is*axis[0]*axis[1]]/lbSolidAngle);
	  }
	} else {
	  std::cerr<<"Pixels fall outside of boundaries in l and b conversion to healpix"<<std::endl;
	  std::cerr<<"ib: "<<ib<<", il: "<<il<<std::endl;
	  std::cerr<<"Phi: "<<pnt.phi*180/utl::kPi<<", theta: "<<90 - pnt.theta*180/utl::kPi<<std::endl;
	}
      }
      size_t first = p*fSpectra.size();
      for (size_t l = 0; l < nSpectra(); ++l){
	fMap[first+l] = T(totSpectra[l]/double(nPix)*hpSolidAngle);
      }
    }
  }
  
#ifdef HAVE_ASTRO
  /**\brief Helper class to create a conversion function for coordinate
   * transformation
   */
  class TransformFunction {
  private:
    void (*tfunc)(double,double,double,double*,double*);
  public:
    TransformFunction(CoordSys from, CoordSys to){
      if (from == to)
	tfunc = &empty;
      else if (from == GAL) {
	if (to == EQ)
	  tfunc = &gal_eq2;
	else
	  tfunc = &gal_ecl2;
      } else if (from == EQ) {
	if (to == ECL)
	  tfunc = &eq_ecl2;
	else
	  tfunc = &gal_eq2;
      } else {
	if (to == GAL)
	  tfunc = &ecl_gal2;
	else
	  tfunc = &ecl_eq2;
      }
    }
    void operator () (double mj, const pointing &from, pointing &to){
      double fl = from.phi;
      double fb = PI/2 - from.theta;
      double tl, tb;
      (*tfunc)(mj,fb,fl,&tb,&tl);
      to.phi = tl;
      to.theta = PI/2 - tb;
    }
  };
  
#endif
  
 public:
  /**\brief Default constructor
   *
   * Sets the size to 0 and scheme to RING
   */
  Skymap(){ Resize(0,0); }
  
  /**\brief Constructor that takes an order of the skymap, the value at which the spectra is evaluated
   * an ordering scheme and a default value for the map.
   *
   * \param order is the order of the healpix coordinate system.  The number
   * of sides per base pixel,  nside, is 2^order.
   * \param spectra is the values at which the spectra is evaluated in each pixel
   * \param scheme is the ordering scheme of the healpix array
   * \param default_value is the default value for the skymap
   */
  Skymap(const int order, 
	 const std::valarray<double> spectra, 
	 const Healpix_Ordering_Scheme scheme=RING, 
	 const T & default_value = T(0)){
    Resize(order,spectra,scheme,default_value);
  }
  
  /**\brief Constructor that takes an order of the skymap, the
   * boundaries of the spectral bins, an ordering scheme and a default value for the map.
   *
   * \param order is the order of the healpix coordinate system.  The number
   * of sides per base pixel,  nside, is 2^order.
   * \param specMin are the lower boundaries of the spectral bins
   * \param specMax are the upper boundaries of the spectral bins
   * \param scheme is the ordering scheme of the healpix array
   * \param default_value is the default value for the skymap
   */
  Skymap(const int order, 
	 const std::valarray<double> specMin, 
	 const std::valarray<double> specMax, 
	 const Healpix_Ordering_Scheme scheme=RING, 
	 const T & default_value = T(0)){
    Resize(order,specMin,specMax,scheme,default_value);
  }
  
  /**\brief Constructor that takes an order of the skymap, the size of the spectra,
   * an ordering scheme and a default value for the map.
   *
   * \param order is the order of the healpix coordinate system.  The number
   * of sides per base pixel,  nside, is 2^order.
   * \param nSpectra is the size of the spectra
   * \param scheme is the ordering scheme of the healpix array (defaults to RING)
   * \param default_value is the default value for the skymap (defaults to T(0))
   *
   * The values at which the spectra is located is set to 1 in all cases.  This is done so we can take the
   * logarithm of the values, but it is expected that the values at which the spectra is evaluated is in most
   * cases logarithmically distributed.
   */
  Skymap(const int order, const int nSpectra, const Healpix_Ordering_Scheme scheme=RING, const T & default_value = T(0)){
    Resize(order,nSpectra,scheme,default_value);
  }
  
  /**\brief Construct a skymap from file, either a healpix fits file or CAR
   * projected fits image.
   *
   * \param fileName is the name of the file to be opened
   */
  Skymap(const std::string & fileName, int order=-1) {
    load(fileName, order);
  }
  
  /**\brief Copy constructor */
  Skymap(const Skymap<T> & oldMap){
    if (oldMap.fBinned){
      Resize(oldMap.order_, oldMap.fSpecMin, oldMap.fSpecMax, oldMap.scheme_);
    }else{
      Resize(oldMap.order_, oldMap.fSpectra, oldMap.scheme_);
    }
    fMap = oldMap.fMap;
  }
  
  /**\brief Resize the map to order and size of spectra
   *
   * \param order is the order of the healpix coordinate system
   * \param nSpectra is the size of the spectral array
   * \param scheme is the ordering scheme of the healpix map (defaults to RING)
   * \param defaultValue is the default value for the map (defaults to T(0))
   *
   * This method destroys the data in the skymap.  The value at which the spectra is evaluated is set to 1.
   */
  void Resize(const int order, const int nSpectra, const Healpix_Ordering_Scheme scheme=RING, const T & defaultValue = T(0)){
    std::valarray<double> spectra(1.0,nSpectra);
    Resize(order, spectra, scheme, defaultValue);
  }
  
  /**\brief Resize the map to order and size of spectra
   *
   * \param order is the order of the healpix coordinate system
   * \param spectra is the values at which the spectra is evaluated
   * \param scheme is the ordering scheme of the healpix map (defaults to RING)
   * \param defaultValue is the default value for the map (defaults to T(0))
   *
   * This method destroys the data in the skymap.
   */
  void Resize(const int order, const std::valarray<double> &spectra, const Healpix_Ordering_Scheme scheme=RING, const T & defaultValue = T(0)){
    Set(order, scheme); // Call the construcor for healpix with correct arguments
    fSpectra.resize(spectra.size());
    fSpectra = spectra;
    fBinned = false;
    fMap.resize(Npix()*fSpectra.size(), defaultValue);
  }
  
  /**\brief Resize the map to order and size of spectra in a binned fashion
   *
   * \param order is the order of the healpix coordinate system
   * \param specMin are the lower boundaries of the spectral bins
   * \param specMax are the upper boundaries of the spectral bins
   * \param scheme is the ordering scheme of the healpix map (defaults to RING)
   * \param defaultValue is the default value for the map (defaults to T(0))
   *
   * This method destroys the data in the skymap.
   */
  void Resize(const int order, const std::valarray<double> &specMin, const std::valarray<double> &specMax, const Healpix_Ordering_Scheme scheme=RING, const T & defaultValue = T(0)){
    //The spectral sizes must be the same
    if (specMin.size() != specMax.size()){
      std::cerr<<"Spectral sizes not equal for boundary arrays"<<std::endl;
      std::cerr<<specMin.size()<<" != "<<specMax.size()<<std::endl;
      throw(1);
    }
    Set(order, scheme); // Call the construcor for healpix with correct arguments
    fSpectra.resize(specMin.size());
    fSpecMax.resize(specMin.size());
    fSpecMin.resize(specMin.size());
    fSpectra = 0.5*(specMin+specMax);
    fSpecMin = specMin;
    fSpecMax = specMax;
    fBinned = true;
    fMap.resize(Npix()*fSpectra.size(), defaultValue);
  }
  
  /**\brief Load a skymap from a file
   *
   * \param fileName is the name of the file to load.
   * \param order is the order of the skymap if we are rebinning
   * CAR maps.  If -1, it is automatically determined from the
   * CAR binning.  On return it is set to the order of the skymap
   *
   * Tries to be smart and first checks if it is a healpix file, and then
   * tries to load a standard FITS image.  This method is limited to CAR
   * projection for the moment and it does not check!!
   */
  void load(const std::string &fileName, int order=-1){
    //Create the ccfits object, read only, finding the
    //extension with the pixtype HEALPIX and emin keyword
    try {
      std::vector< std::string > keywords(1,"");
      std::vector< std::string > values(1,"");
      keywords[0] = "PIXTYPE";
      values[0] = "HEALPIX";
      CCfits::FITS fits(fileName, CCfits::Read, keywords, values);
      
      //A reference to the table containing the skymap data
      CCfits::ExtHDU &skymapTable = fits.currentExtension(); 
      
      //Read the keywords to set up the skymap object
      int nRows, nSpectra, nSide;
      std::string ordering;
      skymapTable.readKey("NAXIS2", nRows);
      skymapTable.readKey("NBRBINS", nSpectra);
      skymapTable.readKey("NSIDE", nSide);
      skymapTable.readKey("ORDERING", ordering);
      //Calculate the order
      int hporder = int(log(double(nSide))/log(2.0)+0.1);// original: int order = int(log(nSide)/log(2)+0.1); AWS20080519
      
      //Try to find the EBOUNDS or ENERGIES extensions, either must exist
      bool foundEnExt = false;
      //First try the energies table
      try {
	fits.read(std::vector<std::string> (1,"ENERGIES"));
	CCfits::ExtHDU &energyTable = fits.extension("ENERGIES");
	std::valarray<double> energies;
	energyTable.column(1).read(energies, 1, nSpectra);
	
	if (ordering == "NEST") {
	  if (Order() != hporder || nSpectra != fSpectra.size() || Scheme() != NEST ) {
	    Resize(hporder, energies, NEST);
	  } else {
	    setSpectra(energies);
	  }
	} else {
	  if (Order() != hporder || nSpectra != fSpectra.size() || Scheme() != RING ) {
	    Resize(hporder, energies, RING);
	  } else {
	    setSpectra(energies);
	  }
	}
	foundEnExt=true;
      } catch (CCfits::FITS::NoSuchHDU) {
	try { //Then EBOUNDS table
	  fits.read(std::vector<std::string> (1,"EBOUNDS"));
	  CCfits::ExtHDU &energyTable = fits.extension("EBOUNDS");
	  std::valarray<double> eMin, eMax;
	  energyTable.column(2).read(eMin, 1, nSpectra);
	  energyTable.column(3).read(eMax, 1, nSpectra);
	  
	  if (ordering == "NEST") {
	    if (Order() != hporder || nSpectra != fSpectra.size() || Scheme() != NEST ) {
	      Resize(hporder, eMin, eMax, NEST);
	    } else {
	      setSpectra(eMin, eMax);
	    }
	  } else {
	    if (Order() != hporder || nSpectra != fSpectra.size() || Scheme() != RING ) {
	      Resize(hporder, eMin, eMax, RING);
	    } else {
	      setSpectra(eMin, eMax);
	    }
	  }
	  foundEnExt=true;
	} catch (CCfits::FITS::NoSuchHDU) { }
      }
      
      //If we found an energy extension, read in the data
      if (foundEnExt) {
	//Read in the data from the table
	//Skip the unnecessary overhead of using CCfits this
	//time and use cfitsio routines
	//Create a map of typeid's to cfitsio datatype
	std::map<const char*, int> formatMap;
	formatMap[typeid(char).name()] = TSBYTE;
	formatMap[typeid(short).name()] = TSHORT;
	formatMap[typeid(int).name()] = TINT;
	formatMap[typeid(long).name()] = TLONG;
	formatMap[typeid(float).name()] = TFLOAT;
	formatMap[typeid(double).name()] = TDOUBLE;
	formatMap[typeid(unsigned char).name()] = TBYTE;
	formatMap[typeid(unsigned short).name()] = TUSHORT;
	formatMap[typeid(unsigned int).name()] = TUINT;
	formatMap[typeid(unsigned long).name()] = TULONG;
	//Select the appropriate datatype
	int dataType = formatMap[typeid(T).name()];
	//Point the fits file pointer to the correct extension
	skymapTable.makeThisCurrent();
	//Get the fits pointer
	fitsfile* ff = fits.fitsPointer();
	//Read the data
	int status(0), anynul(0);
	T null(0);
	fits_read_col(ff, dataType, 1, 1, 1, fMap.size(), &null, &fMap[0], &anynul, &status);
	if (status != 0) {
	  fits_report_error(stderr, status);
	  throw(1);
	}
      } else {
	std::cerr<<"Not a compatible fits file, did not find an EBOUNDS or ENERGIES extension"<<std::endl;
	throw(std::string("Not a compatible fits file, did not find an EBOUNDS or ENERGIES extension"));
      }
    } catch (CCfits::FITS::NoSuchHDU) {
      //Assume we have a CAR fits image
      CCfits::FITS fits(fileName);
      CCfits::PHDU &mapCube = fits.pHDU();
      
      //Read the number of axes
      long axes = mapCube.axes();
      
      //Throw an error if the number of axes is less than 2
      if (axes < 2) throw(std::string("Number of axis less than 2, cannot continue"));
      
      //We take at most 3 axes into account
      axes = std::min(long(3),axes);
      
      //Create a vector of keywords for the CR values
      std::stringstream ss;  //To write the numbers to
      std::valarray<double> crval(0.0, 3), crpix(1.0, 3), cdelt(1.0, 3);
      std::valarray<long> axis(3);
      for (int i = 0; i < axes; ++i) {
	axis[i] = mapCube.axis(i);
	//Seek to the beginning of the stringstream to overwrite old values
	ss.str("");
	ss << "CRPIX" << i+1;
	try {
	  mapCube.readKey(ss.str(), crpix[i]);
	} catch (CCfits::HDU::NoSuchKeyword) {} //Assume the value is 1 if undefined
	//Seek to the beginning of the stringstream to overwrite old values
	ss.str("");
	ss << "CDELT" << i+1;
	try {
	  mapCube.readKey(ss.str(), cdelt[i]);
	} catch (CCfits::HDU::NoSuchKeyword) {
	  //Assume whole sky maps and 1 for all
	  //other axis
	  if (i == 0) {
	    cdelt[i] = 360./axis[i];
	  } else if (i == 1) {
	    cdelt[i] = 180./axis[i];
	  }
	} 
	//Seek to the beginning of the stringstream to overwrite old values
	ss.str("");
	ss << "CRVAL" << i+1;
	try {
	  mapCube.readKey(ss.str(), crval[i]);
	} catch (CCfits::HDU::NoSuchKeyword) {
	  //Assume full sky maps and 0 for everything else
	  if (i == 0) {
	    crval[i] = 0 + cdelt[i]/2.;
	  } else if (i == 1) {
	    crval[i] = -90 + cdelt[i]/2.;
	  }
	} 
      }
      
      //Read the data and resize the skymap.  Let the skymap be of RING
      //structure, most favorable for convolution.
      //The resolution determines the order for the map, which is always
      //greater or equal to the resolution of the map
      double res = std::min(fabs(cdelt[0]),fabs(cdelt[1]));
      if (order == -1)
	order = int(log(sqrt(3./utl::kPi)*60/res)/log(2.0)) + 1; //log(nside)/log(2) + 1      log(2)->log(2.0) AWS20080519
      
      //If this is a proper mapcube, read in the energies value
      std::valarray<double> energies, emin, emax;
      try {
	//Try energies extension first
	CCfits::ExtHDU & energyTable = fits.extension("ENERGIES");
	int nSpectra;
	energyTable.readKey("NAXIS2", nSpectra);
	energyTable.column(1).read(energies, 1, nSpectra);
      } catch (CCfits::FITS::NoSuchHDU) {
	try{
	  //Then ebounds
	  CCfits::ExtHDU & energyTable = fits.extension("EBOUNDS");
	  int nSpectra;
	  energyTable.readKey("NAXIS2", nSpectra);
	  energyTable.column(2).read(emin, 1, nSpectra);
	  energyTable.column(3).read(emax, 1, nSpectra);
	  energies.resize(emin.size());
	} catch (CCfits::FITS::NoSuchHDU) {}
      }
      
      //If no energies extension is found, create the values from CRVAL and
      //CDELT values
      if (energies.size() == 0) {
	if (axes < 3) {
	  energies.resize(1);
	  energies[0] = 1;
	  axis[2] = 1;
	}else{
	  energies.resize(axis[2]);
	  //Assume CDELT represents logarithmic values
	  for (int i = 0; i < axis[2]; ++i) {
	    energies[i] = pow(10,crval[2] + i*cdelt[2]);
	  }
	}
      }
      
      //Now we can set up the map
      if (emin.size() == 0){
	if (Order() != order || energies.size() != fSpectra.size() || Scheme() != RING ) {
	  Resize(order, energies, RING);
	}else{
	  setSpectra(energies);
	}
      } else {
	if (Order() != order || emin.size() != fSpectra.size() || Scheme() != RING ) {
	  Resize(order, emin, emax, RING);
	}else{
	  setSpectra(emin, emax);
	}
      }
      
      //Read the image data and fill the skymap
      std::valarray<T> image;
      long npix(1);
      for (int i = 0; i < axes; ++i){
	npix *= axis[i];
      }
      mapCube.read(image,1,npix);
      fillMapcube(image, axis, crval, cdelt, crpix);
      std::cout<<"All done"<<std::endl;
    }
  }
  
  /**\brief Write the skymap to a file
   *
   * \param fileName is the name of the file to write to
   *
   * The file is overwritten without warning.
   */
  void write(const std::string & fileName) const {
    //Do nothing if there is no data
    if (fSpectra.size() == 0) return;
    //Prepend ! so the file gets overwritten
    std::string fname = fileName[0] != '!' ? "!" + fileName : fileName;
    //Append .gz to allow compression
#ifdef ENABLE_COMPRESSION
    if ( fname.substr(fname.size()-3).compare(".gz") )
			   fname += ".gz";
#endif
    
    //Create a CCFITS object with the filename
    CCfits::FITS fits(fname, CCfits::Write );
    
    //Create a map of typeid's to cfitsio format characters
    std::map<const char*, std::string> formatMap;
    formatMap[typeid(char).name()] = "S";
    formatMap[typeid(short).name()] = "I";
    formatMap[typeid(int).name()] = "J";
    formatMap[typeid(long).name()] = "K";
    formatMap[typeid(float).name()] = "E";
    formatMap[typeid(double).name()] = "D";
    formatMap[typeid(unsigned char).name()] = "B";
    formatMap[typeid(unsigned short).name()] = "U";
    formatMap[typeid(unsigned int).name()] = "V";
    formatMap[typeid(unsigned long).name()] = "K";
    //Create the data table to store the actual image.  The extension name is
    //SKYMAP
    std::ostringstream strForm;
    strForm << fSpectra.size() << formatMap[typeid(T).name()];
    std::vector<std::string> colNames(1,"Spectra"), colForm(1,strForm.str()), colUnit(1,"");
    CCfits::Table *imgTable = fits.addTable("SKYMAP", Npix(), colNames, colForm, colUnit);
    
    //Write the data
    if ( fSpectra.size() == 1 ) {
      imgTable->column("Spectra").write(fMap,1);
    }else{
      imgTable->column("Spectra").write(fMap,Npix(),1);
    }
    
    
    //Write keywords
    imgTable->addKey("PIXTYPE", "HEALPIX", "Healpix pixel scheme");
    std::string keyvalue = "RING";
    if (NEST == Order()) keyvalue = "NEST";
    imgTable->addKey("ORDERING", keyvalue, "Ring or nested ordering of pixels");
    imgTable->addKey("NSIDE", Nside(), "Number of sides in a base pixel");
    imgTable->addKey("FIRSTPIX", 0, "Number of first pixel");
    imgTable->addKey("LASTPIX", Npix()-1, "Number of last pixel");
    imgTable->addKey("NBRBINS", nSpectra(), "Number of energy bins in spectra");
    imgTable->addKey("EMIN", fSpectra[0], "Minimum energy of spectra (MeV)");
    double keyv = 0;
    if (fSpectra.size() > 1){
      keyv = log(fSpectra[1]/fSpectra[0]);
    }
    imgTable->addKey("EBIN", keyv, "Energy bin size (logarithmic)");
    imgTable->addKey("COORDTYPE", "GAL", "");
    
    //Create another table to store the energy
    if (fBinned){
      colNames.resize(3); colUnit.resize(3); colForm.resize(3);
      colNames[0] = "CHANNEL"; colNames[1] = "E_MIN"; colNames[2] = "E_MAX";
      colUnit[0] = ""; colUnit[1] = ""; colUnit[2] = "";
      colForm[0] = "I"; colForm[1] = "D"; colForm[2] = "D";
      CCfits::Table *eTable = fits.addTable("EBOUNDS", nSpectra(), colNames, colForm, colUnit);
      //Create the channel array
      std::valarray<int> channel(fSpectra.size());
      for (int i = 0; i < int(channel.size()); ++i){
	channel[i] = i+1;
      }
      eTable->column("CHANNEL").write(channel,1);
      eTable->column("E_MIN").write(fSpecMin,1);
      eTable->column("E_MAX").write(fSpecMax,1);
    }else{
      colNames[0] = "Energy";
      colUnit[0] = "";
      colForm[0] = "D";
      CCfits::Table *eTable = fits.addTable("ENERGIES", nSpectra(), colNames, colForm, colUnit);
      
      eTable->column("Energy").write(fSpectra,1);
    }
  }
  
  /**\brief Convert all the pixels to a different type and returns the new copy of the map
   *
   * \param dummy controls the output type
   */
  template <typename C>
    Skymap<C> convert(C dummy) const{
    //Create an output, the same size and same spectra
    Skymap<C> output;
    if (fBinned) {
      output.Resize(order_, fSpecMin, fSpecMax, scheme_);
    } else {
      output.Resize(order_, fSpectra, scheme_);
    }
    for (size_t i = 0; i < Npix(); ++i){
      for (size_t j = 0; j < nSpectra(); ++j){
	output[i][j] = C((*this)[i][j]);
      }
    }
    return output;
  }
  
  /** \brief Interpolate the map to a different order.
   *
   * \param order specifies the new order.  It makes little sense
   * to interpolate to lower order, but it is not forbidden.
   *
   * Internally we use the get_interpol method of Healpix_Base to
   * calculate the interpolation value at the centers of the new
   * pixels
   */
  Skymap<T> interpolate(int order) const{
    //If the order is the same, just return
    if (order == order_) return *this;
    
    //Create a new skymap with the new order, keeping the spectra and the
    //scheme_
    Skymap<T> newMap;
    if (fBinned) {
      newMap.Resize(order, fSpecMin, fSpecMax, scheme_);
    }else{
      newMap.Resize(order, fSpectra, scheme_);
    }
    
    //Loop the new map and calculate the interpolation
    fix_arr<int,4> pixels;
    fix_arr<double,4> weight;
    for (int i = 0; i < newMap.Npix(); ++i) {
      get_interpol(newMap.pix2ang(i),pixels,weight);
      for (int j = 0; j < 4; ++j){
	for (int k = 0; k < fSpectra.size(); ++k){
	  newMap[i][k] += weight[j]*(*this)[pixels[j]][k];
	}
      }
    }
    return newMap;
  }
  
  
  
  /** \brief Rebin the map to a different order.
   *
   * \param order specifies the new order.  If it is greater than the old
   * order, the map size is increased and the new map contains multiple
   * pixels with the same values.  If it is less than the old order, the new
   * pixel values will be the average of the old pixels contained in the new
   * set.
   *
   * It is not wise to rebin the map to a higher order, rather use the
   * coordinate selector to get the pixel values.  That is what this
   * rebinning does anyway.
   */
  Skymap<T> rebin(int order, bool SAcorrect = true) const{
    //If the order is the same, just return
    if (order == order_) return *this;
    
    //Create a new skymap with the new order, keeping the spectra and the
    //scheme_
    Skymap<T> newMap;
    if (fBinned) {
      newMap.Resize(order, fSpecMin, fSpecMax, scheme_);
    }else{
      newMap.Resize(order, fSpectra, scheme_);
    }
    
    //What we do now depends if the order, if it is greater, loop the new map
    //and set the pixel from their coordinate
    if (order > order_){
      int dOrder = order - order_;
      const int nPix = (1<<dOrder)*(1<<dOrder);
      for (Iterator it = newMap.begin(); it != newMap.end(); ++it){
	SM::Coordinate co = it.coord();
	*it = (*this)[co];
      }
      if (!SAcorrect)
	newMap /= double(nPix);
    }else{
      //If the order is less, we have to average over several pixels
      //The difference in the number of pixels
      int dOrder = order_ - order;
      const int nPix = (1<<dOrder)*(1<<dOrder);
      //Create a memory to work with and an
      //ArraySlice
#pragma omp parallel default(shared)
      {
	std::valarray<T> avStore(0.0, fSpectra.size());
	ArraySlice<T> average(&avStore[0], fSpectra.size());
#pragma omp for schedule(static)
	for (long p = 0; p < newMap.Npix(); ++p){
	  average = 0;
	  //Use pixel numbering from the nested scheme
	  const int pp = (newMap.Scheme() == NEST) ? p : newMap.ring2nest(p);
	  for (long po = pp*nPix; po < (pp+1)*nPix; ++po){
	    //Convert pixel back to correct scheme if needed
	    const int pop = (Scheme() == NEST) ? po : nest2ring(po);
	    average += (*this)[pop];
	  }
	  //Divide by the number of pixels to get the average
	  if (SAcorrect) average /= double(nPix);
	  newMap[p] = average;
	}
      }//End parallel
    }
    //Assign the newMap to us
    return newMap;
  }
  
  /**\brief The size of the spectra */
  int nSpectra() const{ return fSpectra.size();}

  /**\brief Set the values at which the spectra is evaluated 
   *
   * \param spectra are the new values, stored in a valarray
   *
   * \return if the size of the input values do not confirm with the spectral size, a false value is returned.
   *
   * This method changes the skymap to unbinned mode
   */
  bool setSpectra(const std::valarray<double> &spectra){
    if (spectra.size() == fSpectra.size()){
      fSpectra = spectra;
      fBinned = false;
      return true;
    }
    return false;
  }
  
  /**\brief Set the values at which the spectra is evaluated in
   * binned mode
   *
   * \param specMin are the new lower boundaries
   * \param specMax are the new upper boundaries
   *
   * \return if the size of the input values do not confirm with the spectral size, a false value is returned.
   *
   * This method changes the skymap to binned mode
   */
  bool setSpectra(const std::valarray<double> &specMin, const std::valarray<double> &specMax){
    if (specMin.size() == fSpectra.size() && specMax.size() == fSpectra.size()){
      fSpecMin.resize(fSpectra.size());
      fSpecMax.resize(fSpectra.size());
      fSpectra = 0.5*(specMin+specMax);
      fSpecMin = specMin;
      fSpecMax = specMax;
      fBinned = true;
      return true;
    }
    return false;
  }
  
  /**\brief Set the values at which the spectra is evaluated
   *
   * \param spectra are the new values in a C array
   * \param size is the size of the new array.
   *
   * \return if the size of the input values do not confirm with the spectral size, a false value is returned.
   */
  bool setSpectra(const double spectra[], const int size){
    if (size == fSpectra.size()){
      for (int i=0; i<fSpectra.size(); ++i){
	fSpectra[i] = spectra[i];
      }
      return true;
    }
    return false;
  }
  
  /**\brief Return the values at which the spectra is evaluated*/
  const std::valarray<double> & getSpectra() const{return fSpectra;}
  
  /**\brief Return the boundaries for binned maps.
   *
   * \return false if the map is not binned
   */
  bool getBoundaries(std::valarray<double> & specMin, std::valarray<double> & specMax) const{
    if (fBinned) {
      specMin.resize(fSpecMin.size());
      specMax.resize(fSpecMax.size());
      specMin = fSpecMin;
      specMax = fSpecMax;
    }
    return fBinned;
  }
  
  /**\brief Fill skymap with data from an l and b map
   *
   * \param map is an one dimensional array formatted in mapcube format (order l,b,spectra)
   * \param nl the number of l bins
   * \param nb the number of b bins
   * \param correctSA is a boolean which controls wether we correct for solid angle or not.
   * 
   * The spectral size is assumed to be the same as the skymap and no checking is performed.
   * This method splits the healpix map into finer bins and selects the undelying l and b pixel.
   * The input map must be full sky and the first pixel in the input map must have the edge at 
   * l = 0 and b = -90.  It is assumed that the pixel values are given at the center of the pixels,
   * so the pixel (ib, il) spans the area  l = [il*360/nl, (il+1)*360/nl] and b = [ib*180/nl-90, (ib+1)*180/nl-90].
   * Note that ib and il are numbered from 0.
   */
  void filllbCARarray( const T map[], const int nl, const int nb, const bool correctSA=false){
    //Set up the crval, cdelt, crpix and pass it on
    std::valarray<double> crval(0.0, 3), cdelt(1.0,3), crpix(1.0,3);
    std::valarray<long> axis(3);
    
    //Full skymaps
    cdelt[0] = 360./nl;
    cdelt[1] = 180./nb;
    crval[0] = cdelt[0]/2.;
    crval[1] = -90+cdelt[1]/2.;
    axis[0] = nl;
    axis[1] = nb;
    axis[2] = nSpectra();
    
    //The valarray for the image data
    std::valarray<T> image(map,axis[0]*axis[1]*axis[2]);

    fillMapcube(image,axis,crval,cdelt,crpix,correctSA);
  }
  
  /**\bried Convert to a mapcube in CAR projection in galactic
   * coordinates and write it to a file
   *
   * Spatial coordinates are given as in the fits standard, with
   * vectors of length 2: NAXIS, CRVAL, CDELT, CRPIX.  The first
   * element corresponds to longitude but the second to latitude.
   *
   * If NAXIS is 0, full sky is assumed, trying to take CDELT
   * and CRVAL into account.  If CDELT is also 0, the resolution
   * of the skymap is used.  We put a pixel at 0,0 by default.
   * 
   * \param fileName is the name of the output file, with full
   * path.
   * \param NAXIS is the number of pixels per axis (length 2)
   * \param CRVAL is the value of pixel CRPIX (length 2)
   * \param CDELT is the difference in value between two pixels
   * (length 2)
   * \param CRPIX is the pixel number that is at CRVAL (note that
   * pixels are numbered from 1 like when you count, not from 0
   * like in C/C++)
   * \param correctSA should be true if you are working with
   * counts, false if you are working with flux.
   *
   * \return a 1D valarray that contains all the pixels, where l
   * loops fastest, then b and last the energy.
   *
   * Note that if the length of NAXIS et al. is greater than 2,
   * the values are simply ignored.  There is no way to rebin in
   * energy.  On output, the arrays will be resized to length 3
   * and the corresponding value for the energy dimension added
   * to it.  It will be in log energy.
   *
   * This routine works by splitting the output pixels up into
   * 1024 smaller pixels and assigning it the value that is
   * directly under the center.  Then the final bigger pixels is
   * achieved by taking the average.
   */
  void writeMapcube(const std::string &fileName, 
		    std::vector<long> &NAXIS, 
		    std::vector<double> &CRVAL, 
		    std::vector<double> &CDELT, 
		    std::vector<double> &CRPIX, 
		    bool correctSA = false) const{
    //Check for consistency in the input parameters
    NAXIS.resize(2,0);
    CDELT.resize(2,0.0);
    CRVAL.resize(2,0.0);
    CRPIX.resize(2,1.0);
    double res = resolution();
    if (NAXIS[0] == 0) {
      if (CDELT[0] == 0) {
	NAXIS[0] = long(360/res);
	CDELT[0] = 360./double(NAXIS[0]);
      } else {
	NAXIS[0] = long(360/CDELT[0]);
      }
      CRVAL[0] = 0;
      CRPIX[0] = (NAXIS[0]+1)/2;
    }
    if (NAXIS[1] == 0) {
      if (CDELT[1] == 0) {
	NAXIS[1] = long(180/res);
      } else {
	NAXIS[1] = long(180/CDELT[1]);
      }
      if (NAXIS[1] % 2) {
	CDELT[1] = 180./double(NAXIS[1]-1);
      } else {
	CDELT[1] = 180./double(NAXIS[1]);
	++NAXIS[1];
      }
      CRVAL[1] = 0;
      CRPIX[1] = (NAXIS[1]+1)/2;
    }
    //Fix CDELT, if it is still 0
    if (CDELT[0] == 0) {
      CDELT[0] = 360./double(NAXIS[0]);
    }
    if (CDELT[1] == 0) {
      CDELT[1] = 180./double(NAXIS[1]);
    }
    
    //Do the conversion
    std::valarray<T> mapCube = toMapcube(NAXIS,CRVAL,CDELT,CRPIX,correctSA);
    
    //Create the fits file
    CCfits::FITS fits(fileName, DOUBLE_IMG, NAXIS.size(), &NAXIS[0]);
    
    fits.pHDU().write(1, NAXIS[0]*NAXIS[1]*NAXIS[2], mapCube);
    
    //Write the keywords
    fits.pHDU().addKey("CRVAL1", CRVAL[0], "Value of longitude in pixel CRPIX1");
    fits.pHDU().addKey("CDELT1", CDELT[0], "Step size in longitude");
    fits.pHDU().addKey("CRPIX1", CRPIX[0], "Pixel that has value CRVAL1");
    fits.pHDU().addKey("CTYPE1", "GLON-CAR", "The type of parameter 1 (Galactic longitude in CAR projection)");
    fits.pHDU().addKey("CUNIT1", "deg", "The unit of parameter 1");
    fits.pHDU().addKey("CRVAL2", CRVAL[1], "Value of latitude in pixel CRPIX2");
    fits.pHDU().addKey("CDELT2", CDELT[1], "Step size in latitude");
    fits.pHDU().addKey("CRPIX2", CRPIX[1], "Pixel that has value CRVAL2");
    fits.pHDU().addKey("CTYPE2", "GLAT-CAR", "The type of parameter 2 (Galactic latitude in CAR projection)");
    fits.pHDU().addKey("CUNIT2", "deg", "The unit of parameter 2");
    fits.pHDU().addKey("CRVAL3", CRVAL[2], "Energy of pixel CRPIX3");
    fits.pHDU().addKey("CDELT3", CDELT[2], "log10 of step size in energy (if it is logarithmically distributed, see ENERGIES/EBOUNDS extension)");
    fits.pHDU().addKey("CRPIX3", CRPIX[2], "Pixel that has value CRVAL3");
    fits.pHDU().addKey("CTYPE3", "Energy", "Axis 3 is the spectra");
    fits.pHDU().addKey("CUNIT3", "MeV", "The unit of axis 3");

    //Write the energy extension
    if (fBinned){
      std::vector<std::string> colNames(3), colForm(3), colUnit(3);
      colNames[0] = "CHANNEL"; colNames[1] = "E_MIN"; colNames[2] = "E_MAX";
      colUnit[0] = ""; colUnit[1] = ""; colUnit[2] = "";
      colForm[0] = "I"; colForm[1] = "D"; colForm[2] = "D";
      std::valarray<int> channel(nSpectra());
      for (int i = 0; i < channel.size(); ++i){
	channel[i] = i+1;
      }
      CCfits::Table *eTable = fits.addTable("EBOUNDS", nSpectra(), colNames, colForm, colUnit);
      eTable->column("CHANNEL").write(channel,1);
      eTable->column("E_MIN").write(fSpecMin,1);
      eTable->column("E_MAX").write(fSpecMax,1);
    }else{
      std::vector<std::string> colNames(1), colForm(1), colUnit(1);
      colNames[0] = "Energy";
      colUnit[0] = "";
      colForm[0] = "D";
      CCfits::Table *eTable = fits.addTable("ENERGIES", nSpectra(), colNames, colForm, colUnit);
      
      eTable->column("Energy").write(fSpectra,1);
    }
    
  }
  
  /**\bried Convert to a mapcube in CAR projection in galactic
   * coordinates
   *
   * Spatial coordinates are given as in the fits standard, with
   * vectors of length 2: NAXIS, CRVAL, CDELT, CRPIX.  The first
   * element corresponds to longitude but the second to latitude.
   * No effort is put into correcting for input errors.
   *
   * \param NAXIS is the number of pixels per axis (length 2)
   * \param CRVAL is the value of pixel CRPIX (length 2)
   * \param CDELT is the difference in value between two pixels
   * (length 2)
   * \param CRPIX is the pixel number that is at CRVAL (note that
   * pixels are numbered from 1 like when you count, not from 0
   * like in C/C++)
   * \param correctSA should be true if you are working with
   * counts, false if you are working with flux.
   *
   * \return a 1D valarray that contains all the pixels, where l
   * loops fastest, then b and last the energy.
   *
   * Note that if the length of NAXIS et al. is greater than 2,
   * the values are simply ignored.  There is no way to rebin in
   * energy.  On output, the arrays will be resized to length 3
   * and the corresponding value for the energy dimension added
   * to it.  It will be in log energy.
   *
   * This routine works by splitting the output pixels up into
   * 1024 smaller pixels and assigning it the value that is
   * directly under the center.  Then the final bigger pixels is
   * achieved by taking the average.
   */
  std::valarray<T> toMapcube(std::vector<long> &NAXIS, 
			     std::vector<double> &CRVAL, 
			     std::vector<double> &CDELT, 
			     std::vector<double> &CRPIX, 
			     bool correctSA = false) const{
    //Create the output array
    std::valarray<T> output(NAXIS[0]*NAXIS[1]*fSpectra.size());
    
    //Calculating the finer resolution, taking the resolution of
    //the current HEALPix map into account, not making it too
    //small
    const double res = resolution();
    const int nfiner = 32; //The nominal number of subdivision, unless the pixel size of the HEALPix map is small enough
    const int nlfiner = CDELT[0] < res ? int(CDELT[0]/res*nfiner)+1 : nfiner;
    const int nbfiner = CDELT[1] < res ? int(CDELT[1]/res*nfiner)+1 : nfiner;
    int npfiner = nlfiner*nbfiner;
    
    //The finer resolution
    double cdeltl = CDELT[0]/nlfiner;
    double cdeltb = CDELT[0]/nbfiner;
    
    //Add the spectral information to the vectors
    NAXIS.resize(3);
    CRVAL.resize(3);
    CDELT.resize(3);
    CRPIX.resize(3);
    NAXIS[2] = fSpectra.size();
    CRVAL[2] = fSpectra[0];
    CRPIX[2] = 1;
    CDELT[2] = fSpectra.size() > 1 ? log10(fSpectra[1]/fSpectra[0]) : 0;
    
    //Loop over the output pixels, dividing them up into the
    //finer grid
#pragma omp parallel for default(shared) schedule(static) private(npfiner)
    for (long ib = 0; ib < NAXIS[1]; ++ib){
      //Create storage to calculate the average
      std::valarray<T> totStore(T(0),fSpectra.size());
      ArraySlice<T> totSpectra(&totStore[0], fSpectra.size());
      const size_t indb = ib*NAXIS[0];
      //The value of b at the edge of the pixel
      const double b = CRVAL[1] + CDELT[1]*(ib+0.5-CRPIX[1]);
      //Correct for solid angle if necessary
      double lbSolidAngle = 1;
      if (correctSA){
	lbSolidAngle = CDELT[0]*utl::kPi/180. * (sin(utl::kPi*(b+CDELT[1])/180.) - sin(utl::kPi*b/180.));
      }
      for (long il = 0; il < NAXIS[0]; ++il){
	totSpectra = 0;
	const size_t indl = indb + il;
	//The value of l at the edge of the pixel
	const double l = CRVAL[0] + CDELT[0]*(il+0.5-CRPIX[0]);
	//Loop over the finer grid and calculate the sum
	npfiner = 0;
	for (int ibf = 0; ibf < nbfiner; ++ibf){
	  const double bf = b + (ibf+0.5)*cdeltb;
	  if (fabs(bf) >= 90) continue;
	  for (int ilf = 0; ilf < nlfiner; ++ilf){
	    double lf = l + (ilf+0.5)*cdeltl;
	    while (lf < 0) lf += 360;
	    while (lf >= 360) lf -= 360;
	    totSpectra += (*this)[SM::Coordinate(lf,bf)];
	    ++npfiner;
	  }
	}
	//Now assign the average to the pixels
	for (int j = 0; j < totSpectra.size(); ++j){
	  output[j*NAXIS[0]*NAXIS[1]+indl] = totSpectra[j]/double(npfiner)*lbSolidAngle;
	}
      }
    }
    
    //Correct for solid Angle if needed
    if (correctSA)
      output/=solidAngle();
    
    return output;
  }
  
  /**\brief Convert to a l and b map in CAR projection
   *
   * \param nl is the number of l values in the output array
   * \param nb is the number of b values in the output array
   * \param correctSA is a boolean which tells the routine weather to correct for solid angle or not.
   * (When converting counts map, correctSA should be true.)
   *
   * \return a pointer to a dynamically allocated one dimensional c array compatible with mapcubes.  The 
   * l values are looped fastest, then b values and spectra last.  Please deallocate the array once done using
   * it to prevent memory leaks.
   */
  T * lbCARarray( const int nl, const int nb, double CRVAL[3], int CRPIX[3], double CDELT[3], const bool correctSA = false) const{
    T *array;
    array = new T[nl*nb*fSpectra.size()];
    const int nfiner = 32;
    const double lres = 360./nl;
    const double bres = 180./nb;
    const double mapres = resolution();
    //The resolution of the conversion has to be 1/10th of the resolution of
    //the resulting map.  But it only has to have 1/10th of the resolution of
    //the healpix map.
    const int nlfiner = lres < mapres ? int(lres/mapres*nfiner)+1 : nfiner;
    const int nbfiner = bres < mapres ? int(bres/mapres*nfiner)+1 : nfiner;
    const int nboxes = nbfiner*nlfiner;
    //Populate the fits header arrays
    CRVAL[0] = 0; CRVAL[1] = -90+bres/2.; CRVAL[2] = fSpectra[0];
    CRPIX[0] = CRPIX[1] = CRPIX[2] = 1; //fits starts counting at 1
    CDELT[0] = lres; CDELT[1] = bres; CDELT[2] = fSpectra.size() > 1 ? fSpectra[1] - fSpectra[0] : 0;
    //Divide by the solid angle if needed to get the correct results
    double hpSolidAngle = 1;
    if (correctSA){
      hpSolidAngle = solidAngle();
    }
    //Create storage
    std::valarray<T> totStore(T(0),fSpectra.size());
    ArraySlice<T> totSpectra(&totStore[0], fSpectra.size());
#pragma omp parallel for default(shared) schedule(static)
    for (int j = 0; j<nb; ++j){
      const int ind1 = j*nl;
      double bb = bres * j - 90;
      for (int i = 0; i<nl; ++i){
	totSpectra = 0;
	int ind2 = ind1 + i;
	//First pixel's center is at 0
	double ll = lres*i - 0.5*lres;
	//Correct for solid angle if necessary
	double lbSolidAngle = 1;
	if (correctSA){
	  //Note that we use cos, since sin(x-pi/2) = -cos(x)
	  lbSolidAngle = lres * utl::kPi/180 * (cos(utl::kPi/double(nb)*j)-cos(utl::kPi/double(nb)*(j+1))); //AWS20081020
	}
	for (int jj = 0; jj<nbfiner; ++jj){
	  const double b = bb + (jj+0.5)*bres/nbfiner;
	  for (int ii = 0; ii<nlfiner; ++ii){
	    double l = ll + (ii+0.5)*lres/nlfiner;
	    if (l < 0) l += 360;
	    totSpectra += (*this)[SM::Coordinate(l,b)];
	  }
	}
	for (int ll = 0; ll<fSpectra.size(); ++ll){
	  int ind = ind2 +ll*nl*nb;
	  array[ind] = totSpectra[ll]/double(nboxes)/ hpSolidAngle * lbSolidAngle;
	}
      }
    }
    return array;
  }
  
  /**\brief Fill the skymap with data from a Healpix_Map
   *
   * \param iSpectra the point in spectra at which the map should be inserted
   * \param inmap is the Healpix_Map to insert
   */
  void fromHealpixMap(const int iSpectra, const Healpix_Map<T> & inmap){
    //Check for iSpectra bounds and the order of the map
    if (iSpectra < fSpectra.size() && order_ == inmap.Order() && scheme_ == inmap.Scheme()){
      for(int i=0; i<npix_; ++i){
	fMap[i*fSpectra.size()+iSpectra] = inmap[i];
      }
    } else {
      std::cerr<<"Failed to insert Healpix_Map"<<std::endl;
    }
  }
  
  /**\brief Convert the data to a healpix map
   *
   * \param iSpectra the point in spectra for the returned map
   * \return the Healpix map at the corresponding point in spectra
   */
  Healpix_Map<T> toHealpixMap(const int iSpectra) const{
    //Create the output map and fill it with zeros
    Healpix_Map<T> outmap(order_, scheme_);
    outmap.fill(0.0);
    
    //iSpectra must be within bounds, otherwise return a zero map
    if (iSpectra < fSpectra.size()){
      for(int i=0; i<npix_; ++i){
	outmap[i] = fMap[i*fSpectra.size()+iSpectra];
      }
    }
    
    return outmap;
  }
  
  /**\brief Reference to a spectra given a coordinate
   *
   * \param coordinate is an instance of the coordinate class \see SM::Coordinate
   * \return a valarray of the spectra in the pixel corresponding to coordinate.
   */
  ArraySlice<T> operator [] (const SM::Coordinate & coordinate){
    return (*this)[ang2pix(coordinate.healpixAng())];
		}
  /** \overload */
  const ArraySlice<T> operator [] (const SM::Coordinate & coordinate) const{
    return (*this)[ang2pix(coordinate.healpixAng())];
  }
  /**\brief Reference to a spectra given a pixel
   *
   * \param pixel is a pixel number for the skymap
   * \return a std::valarray of the spectra in the pixel
   */
  ArraySlice<T> operator [] (const int pixel){
    return ArraySlice<T>(&fMap[pixel*fSpectra.size()], fSpectra.size());
  }
  /** \overload */
  const ArraySlice<T> operator [] (const int pixel) const{
    return ArraySlice<T>(&fMap[pixel*fSpectra.size()], fSpectra.size());
  }
  
  
  /** \brief Bidirectional iterator for the skymap */
  class Iterator : public std::iterator<std::bidirectional_iterator_tag, std::vector<T> > {
  private:
    int m_index; //!< The index for the array
    Skymap<T> & m_map; //!< A reference to a skymap
    Iterator(){} //!< The default iterator constructor should not be allowed
  public:
    /**\brief Constructor */
  Iterator(const int index, Skymap<T> & map) : m_index(index), m_map(map) {};
    /**\brief Pre-increment operator */
    Iterator & operator ++ () { ++m_index; return *this; }
    /**\brief Post-increment operator */
    Iterator operator ++ (int i) { Iterator tmp(*this); ++m_index; return tmp; }
    /**\brief Pre-decrement operator */
    Iterator & operator -- () { --m_index; return *this; }
    /**\brief Post-decrement operator */
    Iterator operator -- (int i) { Iterator tmp(*this); --m_index; return tmp; }
    /**\brief Dereference operator */
    ArraySlice<T> operator * () { return m_map[m_index];}
    /**\brief Unequal operator */
    bool operator != (const Iterator compare) { return m_index != compare.m_index; }
    /**\brief Return the corresponding coordinate */
    SM::Coordinate coord () { SM::Coordinate co(m_map.pix2ang(m_index)); return co; }
  };
  /** \brief Constant bidirectional iterator for the skymap */
  class constIterator : public std::iterator<std::bidirectional_iterator_tag, std::vector<T> > {
  private:
    int m_index; //!< The index for the array
    const Skymap<T> & m_map; //!< A reference to a skymap
    constIterator(){} //!< The default iterator constructor should not be allowed
  public:
    /**\brief Constructor*/
  constIterator(const int index, const Skymap<T> & map) : m_index(index), m_map(map) {};
    /**\brief Pre-increment operator */
    constIterator & operator ++ () { ++m_index; return *this; }
    /**\brief Post-increment operator */
    constIterator operator ++ (int i) { Iterator tmp(*this); ++m_index; return tmp; }
    /**\brief Pre-decrement operator */
    constIterator & operator -- () { --m_index; return *this; }
				/**\brief Post-decrement operator */
    constIterator operator -- (int i) { Iterator tmp(*this); --m_index; return tmp; }
    /**\brief Dereference operator */
    const ArraySlice<T> operator * () { return m_map[m_index];}
    /**\brief Unequal operator */
    bool operator != (const constIterator compare) { return m_index != compare.m_index; }
    /**\brief Return the corresponding coordinate */
    SM::Coordinate coord () { SM::Coordinate co(m_map.pix2ang(m_index)); return co; }
  };
  
  
  /**\brief Returns an Iterator for the first pixel */
  Iterator begin() { return Iterator(0, *this);}
  /** \overload */
  constIterator begin() const { return constIterator(0, *this);}
  /**\brief Returns an Iterator to the last+1 pixel */
  Iterator end() { return Iterator(npix_, *this); }
  /** \overload */
  constIterator end() const { return constIterator(npix_, *this); }
  
#ifdef HAVE_ASTRO
  /**\brief Convert the skymap from one coordinate system to the
   * other
   *
   * \parameter from is the coordinate system the map is
   * currently in
   * \parameter to is the coordinate system to transfer to
   * \parmeter mj is the modified Julian date
   *
   * \return Skymap converted to the new coordinate system
   *
   * Uses similar methods as transformation to a flat CAR
   * projection, creating a finer grid for the transformation and
   * rebinning to the same grid size again.
   */
  Skymap<T> CoordTransform(CoordSys from, CoordSys to, double mj){
    //If from and to are the same, just return this
    if (from == to)
      return *this;
    //Create the transformation function.  We do it reversed,
    //since we take the to coordinate and transfer to from
    TransformFunction f(to,from);
    //Create the output map and a finer healpix grid to do the
    //transformation on
    Skymap<T> out(*this);
    const int dOrder = std::min(5,13-order_); //Order can not exceed 13
    const int nPix = (1<<dOrder)*(1<<dOrder);
    //Create a nested healpix base for us to work with
    const Healpix_Base finer(order_+dOrder, NEST);
    //Loop over the map, then the finer base and fill the map
#pragma omp parallel for default(shared) schedule(static)
    for (int p = 0; p<npix_; ++p){
      //Use pixel numbering from the nested scheme
      const int pp = (scheme_ == NEST) ? p : ring2nest(p);
      //To store the total value for all the pixels
      std::valarray<T> totSpectra(T(0), fSpectra.size());
      ArraySlice<T> totsl(&totSpectra[0], fSpectra.size());
      //Loop over all of the subpixels
      pointing pnt1 = pix2ang(p);
      for (int sp = pp*nPix; sp<(pp+1)*nPix; ++sp){
	const pointing pnt = finer.pix2ang(sp);
	pointing old;
	f(mj, pnt, old);
	totsl += (*this)[old];
      }
      totsl /= double(nPix);
      out[p] = totsl;
    }
    return out;
  }
#endif
  
  /**\brief Print the skymap
   *
   * \param os is the output stream to which to print the skymap
   *
   * First outputs the coordinate and then the spectra.  This routine is for debugging of
   * small skymaps.  The output gets very large very quickly.  The SkymapFitsio routines are
   * recommended.
   */
  void print (std::ostream & os){
    for (int i = 0; i < npix_; ++i){
      SM::Coordinate co(pix2ang(i));
      os << co;
      os << std::endl << "  "; //!< First output the coordinates
      for(int j = 0; j < fSpectra.size(); ++j){
	os << fMap[i*fSpectra.size()+j] << " ";
      }
      os<<std::endl;
    }
  }
  
  /** \brief Returns the minimum value of the map */
  T min () const {
    return fMap.min();
  }
  /** \brief Returns the minimum value of the map at a given point in the spectra 
   *
   * \param iSpectra is the point in spectra to evaluate the minimum value.
   */
  T min(int iSpectra) const {
    T m(0);
    if (iSpectra >= 0 && iSpectra < int(fSpectra.size())){
      m = fMap[iSpectra];
      for (int i = 1; i < npix_; ++i){
	m = std::min(fMap[i*fSpectra.size()+iSpectra], m);
      }
    }
    return m;
  }
  
  /** \brief Returns the maximum value of the map */
  T max () const {
    return fMap.max();
  }
  /** \brief Returns the maximum value of the map at a given point in the spectra 
   *
   * \param iSpectra is the point in spectra to evaluate the maximum value.
   */
  T max(int iSpectra) const {
    T m(0);
    if (iSpectra >= 0 && iSpectra < int(fSpectra.size())){
      m = fMap[iSpectra];
      for (int i = 1; i < npix_; ++i){
	m = std::max(fMap[i*fSpectra.size()+iSpectra], m);
      }
    }
    return m;
  }
  
  /**\brief Sum of all the pixels for a point in spectra
   *
   * \param iSpectra is the point in spectra to evaluate the sum
   */
  T sum (int iSpectra) const {
    T s(0);
    //Check if iSpectra is within bounds
    if (iSpectra >= 0 && iSpectra < fSpectra.size()){
      for (int i = 0; i < npix_; ++i){
	s += fMap[i*fSpectra.size()+iSpectra];
      }
    }
    return s;
  }
  
  /**\brief Sum all of the pixels in the map */
  T sum () const {
    return fMap.sum();
  }
  
  /**\brief Assigns a single value to the map */ 
  template <typename C>
    Skymap<T> & operator = (const C number){
    fMap = T(number);
    return(*this);
  }
  
  /**\brief Assignment for maps, basically return the same map */
  Skymap<T> & operator = (const Skymap<T> & oldMap){
    //Avoid self-assignment
    if (this != &oldMap){
      if (oldMap.fBinned){
	Resize(oldMap.order_, oldMap.fSpecMin, oldMap.fSpecMax, oldMap.scheme_);
      }else{
	Resize(oldMap.order_, oldMap.fSpectra, oldMap.scheme_);
      }
      fMap = oldMap.fMap;
    }
    return(*this);
  }
  
  //! Add a number to the current skymap
  template <typename C>
    Skymap<T> & operator += (C number) {
    fMap += T(number);
    return (*this);
  }
  //! Add a number to a skymap
  template <typename C>
    Skymap<T> operator + (C number) const {
    Skymap<T> returnMap(*this);
    returnMap += number;
    return returnMap;
  }
  
  //! Withdraw a number from the current skymap
  template <typename C>
    Skymap<T> & operator -= (C number) {
    fMap -= number;
    return (*this);
  }
  //! Withdraw a number from a skymap
  template <typename C>
    Skymap<T> operator - (C number) const {
    Skymap<T> returnMap(*this);
    returnMap -= number;
    return returnMap;
  }
  
  //! Multiply the current skymap with a number
  template <typename C>
    Skymap<T> & operator *= (C number) {
    fMap *= number;
    return (*this);
  }
  //! Multiply a skymap with a number
  template <typename C>
    Skymap<T> operator * (C number) const {
    Skymap<T> returnMap(*this);
    returnMap *= number;
    return returnMap;
  }
  
  //! Divide the current skymap with a number
  template <typename C>
    Skymap<T> & operator /= (C number) {
    fMap /= number;
    return (*this);
  }
  //! Divide a skymap with a number
  template <typename C>
    Skymap<T> operator / (C number) const {
    Skymap<T> returnMap(*this);
    returnMap /= number;
    return returnMap;
  }
  
  /**\brief Add a skymap to the current one 
   *
   * Only adds the skymap if the values at which the spectra is evaluated are the same and the skymaps are the
   * same size.
   */
  Skymap<T> & operator += (const Skymap<T> & otherMap){
    if (! equiv(otherMap) ){
      std::cerr<<"Skymaps are not equivalent in addition"<<std::endl;
      std::cerr<<"There sizes are not equal or their spectra is not the same"<<std::endl;
      return (*this);
    }
    fMap += otherMap.fMap;
    return (*this);
  }
  
  /**\brief Add a skymap to another skymap
   *
   * Only adds the skymaps if the values at which the spectra is evaluated are the same and the skymaps are the
   * same size.
   */
  Skymap<T> operator + (const Skymap<T> & otherMap) const{
    Skymap<T> returnMap(*this);
    returnMap += otherMap;
    return (returnMap);
  }
  
  /**\brief Subtract a skymap from the current one
   *
   * Only subtracts the skymap if the values at which the spectra is evaluated are the same and the skymaps are the
   * same size.
   */
  Skymap<T> & operator -= (const Skymap<T> & otherMap){
    if (! equiv(otherMap) ){
      std::cerr<<"The skymaps arent equivalent when subtracting them.  Check that their spectra and sizes are the same."<<std::endl;
      return (*this);
    }
    fMap -= otherMap.fMap;
    return (*this);
  }
  
  /**\brief Subtract a skymap from another skymap
   *
   * Only subtracts the skymaps if the values at which the spectra is evaluated are the same and the skymaps are the
   * same size.
   */
  Skymap<T> operator - (const Skymap<T> & otherMap) const{
    Skymap<T> returnMap(*this);
    returnMap -= otherMap;
    return (returnMap);
  }
  
  /**\brief Multiply a skymap with the current one
   *
   * Only multiplies the skymap if the values at which the spectra is evaluated are the same and the skymaps are the
   * same size.
   */
  Skymap<T> & operator *= (const Skymap<T> & otherMap){
    if (! equiv(otherMap) ){
      std::cerr<<"The skymaps arent equivalent when multiplying them.  Check that their spectra and sizes are the same."<<std::endl;
      return (*this);
    }
    fMap *= otherMap.fMap;
    return (*this);
  }
  
  /**\brief Multiply a skymap with another skymap
   *
   * Only multiplies the skymaps if the values at which the spectra is evaluated are the same and the skymaps are the
   * same size.
		 */
  Skymap<T> operator * (const Skymap<T> & otherMap) const{
    Skymap<T> returnMap(*this);
    returnMap *= otherMap;
    return (returnMap);
  }
  
  /**\brief Divide a skymap into the current one
   *
   * Only divides the skymap if the values at which the spectra is evaluated are the same and the skymaps are the
   * same size.
   */
  Skymap<T> & operator /= (const Skymap<T> & otherMap){
    if (! equiv(otherMap) ){
      std::cerr<<"The skymaps arent equivalent when dividing them.  Check that their spectra and sizes are the same."<<std::endl;
      return (*this);
    }
#pragma omp parallel for default(shared) schedule(static) if(fMap.size() > 1e5)
    for (int i = 0; i < fMap.size(); ++i){
      fMap[i] /= otherMap.fMap[i];
    }
    return (*this);
  }
  
  /**\brief Divide a skymap with another skymap
   *
   * Only divides the skymaps if the values at which the spectra is evaluated are the same and the skymaps are the
   * same size.
   */
  Skymap<T> operator / (const Skymap<T> & otherMap) const{
    Skymap<T> returnMap(*this);
    returnMap /= otherMap;
    return (returnMap);
  }
  
  
  /**\brief Equality operator, returns true if every pixel has the same value and the values at which the spectra
   * is evaluated is the same as well.
   */
  bool operator == (const Skymap<T> & otherMap) const{
    //Check for equivalence of the map
    bool map = equiv(otherMap);
    if (! map) return map;
    //And the map, if needed
    for (int i = 0; i < fMap.size(); ++i){
      map &= (fMap[i] == otherMap.fMap[i]);
      if (! map) {
	return map;
      }
    }
    return map;
  }
  
  /** \brief Non-equality operator, the inverse of the equality operator. */
  bool operator != (const Skymap<T> & otherMap) const{
    return ( ! ((*this) == otherMap));
  }
 
  T* IOarray(unsigned int& nElements) const {

    T* array;
    
    array = new T[npix_*fSpectra.size()];
    nElements = 0;
    
    for (int i = 0; i < npix_; ++i){
      for (int j = 0; j < fSpectra.size(); ++j){
	array[nElements] = fMap[i*fSpectra.size()+j];//][j];
	++nElements;
      }
    }
    return array;
  }
  
};

#endif
