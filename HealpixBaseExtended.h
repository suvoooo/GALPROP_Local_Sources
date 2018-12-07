/** \class HealpixBaseExtended
 * \brief An extension for the Healpix_Base class with local modifications
 *
 * A place to put local additions to the Healpix_Base class that are necessary
 * for the Skymap but not dependent on the template type
 */
#ifndef HealpixBaseExtended_h
#define HealpixBaseExtended_h

#include "healpix_base.h"
#include "Coordinate.h"
#include "PhysicalConstants.h"
#include "Region.h"
#include <cmath>
#include <vector>
#include <set>

class HealpixBaseExtended : public Healpix_Base {
   private:
      /**\brief A generator to insert simple ranges into containers */
      class counter {
	 public:
	    counter(int init=0) : n(init) {}
	    int operator () () { return ++n; }
	 private:
	    int n;
      };

   public:
      /**\brief Return the resolution in degrees.
       *
       * This is the square root of the solid angle.
       */
      double resolution() const{return sqrt(3/utl::kPi)*60/nside_;}
      /**\brief Solid angle of each pixel in radians */
      double solidAngle() const{return utl::kPi/(3*nside_*nside_);}
      //! Coordinate for a given pixel
      SM::Coordinate pix2coord(int pix) const { return SM::Coordinate(pix2ang(pix)); }
      //! Pixel for a given coordinate
      int coord2pix(const SM::Coordinate & coord) const { return ang2pix(coord.healpixAng()); }


      /** \brief Return the pixels for a given region of the sky
       *
       * \parameters region is the selected region
       *
       * \returns a set of pixels
       *
       * All pixels that have their centers within the region are selected.
       */
      std::set<int> regionToPixels(const SkyRegion & region) const;
		
      /**\brief Return the pixel numbers within a rectangular region
       * 
       * \param pointingll is the lower left corner of the rectangle
       * \param pointingur is the upper right corner of the rectangle
       * \param listpix is the list of pixels who's center lies within the region
       *
       * This routine assumes a doughnut shaped world, so if the lower limit 
       * of the rectangle is above (north, theta=0, is up here) the upper limit, 
       * the pixels within the region above the lower limit and below the upper
       * limit is selected.  This also applies to the left and right region 
       * which is not as unusual.
       */
      void query_rectangle(const pointing &pointingll, const pointing &pointingur, std::vector<int> &listpix) const;

      /** \brief Return the pixels under a given pixel in a different
       * HealpixBase
       *
       * \parameter hp is the other HEALPix base
       * \parameter pixel is the pixel number in the other base
       * \param listpix is the list of pixels who's center lies within the region
       *
       */
      void query_pixel(const Healpix_Base &hp, int pixel, std::vector<int> &listpix) const;
};

#endif
