#ifndef _rf_RadiationField_h_
#define _rf_RadiationField_h_

#include <map>
#include <string>
#include <valarray>
#include <vector>

#include <CLHEP/Vector/ThreeVector.h>

#include <Skymap.h>

namespace rf {

  // Explicit assumption made is that the supplied radiation field
  // is generated on a regular spatial grid. That is, there are n*m entries 
  // in the grid. For each entry in the grid, there is a range of validity 
  // (given by the cell size). Outside this range the radiation field is 
  // determined by interpolation within the grid. For points outside the grid
  // the radiation field is considered zero

  class RadiationField {

  public:

    enum STELLARCOMPONENT { TOTAL, DIRECT, SCATTERED, TRANSIENT, THERMAL };

#ifdef CLHEP_V1_8
    typedef Hep3Vector ThreeVector;
#else
    typedef CLHEP::Hep3Vector ThreeVector;
#endif

    RadiationField();  
    
    RadiationField(const std::string& filename,
		   const std::valarray<double>& freq,
		   int rebinnedSkymapOrder = 0);

    ~RadiationField();

    // Note: The method below is *very* expensive if the desired 
    // healpix order is different from the rebinned order supplied
    // in the constructor above. Use with caution!

    const Skymap<double> GetSkymap(const ThreeVector& pos, 
				   const STELLARCOMPONENT component,
				   const int healpixOrder);

    // This method is *much* faster if you just want the number density

    const valarray<double> GetNumberDensity(const ThreeVector& pos,
					    const STELLARCOMPONENT component); 

  private:

    const Skymap<double> GetSkymap(const ThreeVector& pos,
				   const STELLARCOMPONENT component);

    bool fCacheBuilt;

    string fPrefix;

    int fRebinnedSkymapOrder;

    unsigned int fNumberOfComponents;

    double fLuminosity, fDustMass;

    std::valarray<double> fFrequency, fWavelength, fStellarComponentLuminosity;

    std::vector<std::string> fStellarComponentName;

    std::map<ThreeVector, std::vector<std::string> > fFilenameData;
  
    std::vector< std::vector<ThreeVector> > fPositionData, fRangeData;

    std::valarray<double> fRData, fZData, fRRangeData, fZRangeData;

    std::vector< std::vector< std::vector<std::string>* > > fFilenameOrderedData;

    std::vector< std::vector< std::vector<Skymap<double>* > > > fSkymapOrderedData;

    std::vector< std::vector< std::vector< valarray<double>* > > > fEnergyDensity;

    void ReadRadiationField(const std::string& filename);
    
    void BuildSkymapCache();
    
    void FlushSkymapCache();
    
    void ClearData();
  
  };

}

#endif
