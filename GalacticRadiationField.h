#ifndef _rf_GalacticRadiationField_h_
#define _rf_GalacticRadiationField_h_

#include <fitsio.h>

#include <string>
//#include <vector>
//#include <cmath>
#include <valarray>

namespace rf {

  class GalacticRadiationField {
    
  public:

    enum Component { STELLAR, SCATTERED, INFRARED, CMB, TOTAL };
    
    GalacticRadiationField(const std::string& filename, 
			   bool useIsotropic = false);
    ~GalacticRadiationField();

    // Energy density returned in (eV cm^-3 micron^-1)*micron

    double GetEnergyDensity(const double wl, 
			    const double r, 
			    const double z, 
			    const Component comp) const;

    // Intensity returned in (eV cm^-2 s^-1 sr^-1 micron^-1)*micron

    double GetIntensity(const double wl,
			const double r,
			const double z,
			const double azimuth,
			const double cosZenith,
			const Component comp) const;

    const std::valarray<double>& GetWavelengthData() const { return fWavelengthData; }
    const std::valarray<double>& GetAzimuthData() const { return fAzimuthData; }
    const std::valarray<double>& GetCosZenithData() const { return fCosZenithData; }
    const std::valarray<double>& GetRData() const { return fRData; }
    const std::valarray<double>& GetZData() const { return fZData; }

  private:

    bool fInitialised;

    double fAzimuthDelta, fCosZenithDelta;

    std::valarray<double> fRData, fZData, fWavelengthData, fAzimuthData, fCosZenithData;

    //float* fEnergyDensity;
    //float* fAngularDistribution;

    std::valarray< std::valarray< std::valarray<double> > > fEnergyDensityStellar, fEnergyDensityScattered, fEnergyDensityInfrared;

    std::valarray< std::valarray< std::valarray< std::valarray< std::valarray<double> > > > > //fIntensityStellar, fIntensityScattered, fIntensityInfrared, 
      fIntensityTotal;

    GalacticRadiationField();

    void ExtractFITSData(const std::string& filename, const bool useIsotropic);
    void ExtractFITSEnergyDensity(fitsfile* fptr);
    void ExtractFITSAngularDistribution(fitsfile* fptr);

    double BlackBodyEnergyDensity(const double energy, 
				  const double kT) const;

  };

}

#endif
