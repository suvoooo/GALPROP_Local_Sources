#ifndef REGION_H
#define REGION_H

#include "Coordinate.h"
#include <iostream>

enum RegionType {
   DISC, //!< Region selected by query_disk method of Healpix
   RECTANGLE //!< A square region on a CAR projection
};

class SkyRegion {
   public:
      SkyRegion();
      /** \brief Create a region with the given parameters
       *  
       *  \parameter type is either DISC or RECTANGLE
       *  \parameter center is the center of the region
       *  \parameter radius is the radius of the DISC, ignored if the region is
       *  a RECTANGLE
       *  \parameter deltal is the longitude range of the RECTANGLE, ignored if the
       *  region is a DISC
       *  \parameter deltab is the latitude range of the RECTANGLE, ignored if the
       *  region is a DISC
       *
       *  The RECTANGLE region takes into account that we are selecting on a
       *  sphere and does not loop from the north pole to the south pole.  A
       *  rectangle that traverses either pole results in two triangles, both with
       *  a vertex at the pole (unless the latitude range covers the whole
       *  sphere).
       */
      SkyRegion(RegionType type, const SM::Coordinate &center, double radius, double deltal, double deltab);
      //! Methods to get information on the region
      const RegionType & type() const;
      const SM::Coordinate & center() const;
      const double & radius() const;
      const double & deltal() const;
      const double & deltab() const;
      /** \brief Read a region from a stream
       *
       * The format of the region is
       * <type> <lcenter> <bcenter> <radius/deltal> <deltab>
       * where <type> is either RECTANGLE or DISC.  DISC requires 3 arguments
       * (center and radius) while RECTANGLE requires 4 (center and width and
       * height).  All dimensions are in degrees.
       */
      friend std::istream & operator >> (std::istream &is, SkyRegion &region);
      /** \brief Write a region to a stream
       *
       * The format is the same as for reading
       */
      friend std::ostream & operator << (std::ostream &os, const SkyRegion &region);
   private:
      RegionType ftype;
      //Parameters of the region
      SM::Coordinate fcenter;
      double fradius, fdeltal, fdeltab;
};
#endif
