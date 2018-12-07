#ifndef SKYSELECTION_H
#define SKYSELECTION_H

#include "Region.h"
#include "HealpixBaseExtended.h"
#include <string>
#include <set>
#include <vector>
#include <map>

class SkySelection {
   public:
      /** \brief Initialize the object with regions from a file
       *
       * The format of the file is
       * <ADD/REMOVE> <Region>
       * where the format of <Region> is defined in Region.h
       */
      SkySelection(const std::string & filename);
      /** \brief Add to the selection from a file
       *
       * The format of the file is
       * <ADD/REMOVE> <Region>
       * where the format of <Region> is defined in Region.h
       */
      void addFromFile(const std::string & filename);
      //! Add a region to the selection
      void addRegion(const SkyRegion & region);
      //! Remove a region from the selection
      void removeRegion(const SkyRegion & region);
      //! Return the selected pixels given a HealpixBaseExtended object
      std::set<int> pixelsForHealpix(const HealpixBaseExtended & hp) const;
   private:
      enum Inclusion {ADD, REMOVE};
      std::vector<std::pair<Inclusion,SkyRegion> > fselection;
};

#endif
