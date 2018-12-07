#include <GalacticRadiationField.h>

//#include <Utilities.h>
#include <PhysicalConstants.h>

//#include <fitsio.h>

#include <iostream>
#include <fstream>

using namespace rf;
using namespace utl;
using namespace std;

const double gKT = kBoltzmann*2.735; // CMB

GalacticRadiationField::
GalacticRadiationField(const string& filename, bool useIsotropic) {
  
  fInitialised = false;

  fAzimuthDelta = fCosZenithDelta = 0;

  ExtractFITSData(filename, useIsotropic);

}

GalacticRadiationField::~GalacticRadiationField() {

  //delete[] fAngularDistribution;
  //delete[] fEnergyDensity;

}

void GalacticRadiationField::ExtractFITSData(const string& filename,
					     const bool useIsotropic) {

  fitsfile* fptrEnergyDensity = 0;
  fitsfile* fptrAngularDistribution = 0;

  int status = 0;


  fits_open_file(&fptrEnergyDensity, filename.c_str(), READONLY, &status);

  if (!fptrEnergyDensity) {

    cerr << "Error opening FITS file " << filename << endl;
    fits_report_error(stderr, status);
    exit(1);

  }

  ExtractFITSEnergyDensity(fptrEnergyDensity);

  const string strippedFilename = filename.substr(0, filename.find(".fits"));

  const string intensityFilename = strippedFilename + "AngularDistribution.fits";

  if (!useIsotropic) {

    fits_open_file(&fptrAngularDistribution, intensityFilename.c_str(), READONLY, &status);

    if (!fptrAngularDistribution) {

      cerr << "Error opening FITS file " << filename << endl;
      fits_report_error(stderr, status);
      
    } else {
      
      ExtractFITSAngularDistribution(fptrAngularDistribution);
      
    }

  }

  fInitialised = true;

}

void GalacticRadiationField::ExtractFITSEnergyDensity(fitsfile* fptr) {

  int naxis, naxis1, naxis2, naxis3, naxis4, status = 0;
  float crval1, crval2, crval3, crval4;
  float crdelt1, crdelt2, crdelt3, crdelt4;

  char comment[1024];

  // Energy density

  fits_read_key(fptr, TINT, "NAXIS", &naxis, comment, &status);
  fits_read_key(fptr, TINT, "NAXIS1", &naxis1, comment, &status);
  fits_read_key(fptr, TINT, "NAXIS2", &naxis2, comment, &status);
  fits_read_key(fptr, TINT, "NAXIS3", &naxis3, comment, &status);
  fits_read_key(fptr, TINT, "NAXIS4", &naxis4, comment, &status);

  fits_read_key(fptr, TFLOAT, "CRVAL1", &crval1, comment, &status);
  fits_read_key(fptr, TFLOAT, "CRVAL2", &crval2, comment, &status);
  fits_read_key(fptr, TFLOAT, "CRVAL3", &crval3, comment, &status);
  fits_read_key(fptr, TFLOAT, "CRVAL4", &crval4, comment, &status);
  
  fits_read_key(fptr, TFLOAT, "CDELT1", &crdelt1, comment, &status);
  fits_read_key(fptr, TFLOAT, "CDELT2", &crdelt2, comment, &status);
  fits_read_key(fptr, TFLOAT, "CDELT3", &crdelt3, comment, &status);
  fits_read_key(fptr, TFLOAT, "CDELT4", &crdelt4, comment, &status);

  long nElements = naxis1*naxis2*naxis3*naxis4, fElement = 1;
  float* energyDensity = new float[nElements];
  float nullValue = 0;
  int anyNull;

  fits_read_img(fptr, TFLOAT, fElement, nElements, &nullValue, energyDensity, 
		&anyNull, &status); 

  //fEnergyDensity = energyDensity;

  const unsigned long rBins = naxis1, zBins = naxis2, wlBins = naxis3, 
    compBins = naxis4;

  //const double rDelta = crdelt1;
  //const double zDelta = crdelt2;
  //const double wavelengthDelta = crdelt3;

  fEnergyDensityStellar.resize(rBins);
  fEnergyDensityScattered.resize(rBins);
  fEnergyDensityInfrared.resize(rBins);

  for (unsigned long i = 0; i < rBins; ++i) {

    fEnergyDensityStellar[i].resize(zBins);
    fEnergyDensityScattered[i].resize(zBins);
    fEnergyDensityInfrared[i].resize(zBins);

    for (unsigned long j = 0; j < zBins; ++j) {

      fEnergyDensityStellar[i][j].resize(wlBins);
      fEnergyDensityScattered[i][j].resize(wlBins);
      fEnergyDensityInfrared[i][j].resize(wlBins);

      fEnergyDensityStellar[i][j] = 0;
      fEnergyDensityScattered[i][j] = 0;
      fEnergyDensityInfrared[i][j] = 0;

      for (unsigned long k = 0; k < wlBins; ++k) {

	const unsigned indexStellar = 
	  0*wlBins*rBins*zBins + k*rBins*zBins + j*rBins + i;

	fEnergyDensityStellar[i][j][k] = energyDensity[indexStellar];

	const unsigned indexScattered = 
	  1*wlBins*rBins*zBins + k*rBins*zBins + j*rBins + i;

	fEnergyDensityScattered[i][j][k] = energyDensity[indexScattered];

	const unsigned indexInfrared = 
	  2*wlBins*rBins*zBins + k*rBins*zBins + j*rBins + i;

	fEnergyDensityInfrared[i][j][k] = energyDensity[indexInfrared];

	//	cout << i << " " << j << " " << k << " " << fEnergyDensityStellar[i][j][k] << endl;

      }

    }

  }
  
  fRData.resize(rBins);

  for (unsigned long i = 0; i < rBins; ++i) {

    const double r = crval1 + (i - 0.5)*crdelt1;// - 0.5*crdelt1;
    fRData[i] = r;

  }

  fZData.resize(zBins);

  for (unsigned long i = 0; i < zBins; ++i) {

    const double z = crval2 + i*crdelt2;
    fZData[i] = z;

  }

  fWavelengthData.resize(wlBins);

  for (unsigned long i = 0; i < wlBins; ++i) {

    const double wl = pow(10.0, double(crval3) + i*double(crdelt3));
    fWavelengthData[i] = wl*m/micron;

  }

  fits_close_file(fptr, &status);

  delete[] energyDensity;

}

void GalacticRadiationField::ExtractFITSAngularDistribution(fitsfile* fptr) {

  //return;

  int naxis, naxis1, naxis2, naxis3, naxis4, status = 0;
  float angval1, angval2, angval3, angval4;
  float angdelt1, angdelt2, angdelt3, angdelt4;

  char comment[1024];

  fits_read_key(fptr, TINT, "NAXIS", &naxis, comment, &status);
  fits_read_key(fptr, TINT, "NAXIS1", &naxis1, comment, &status);
  fits_read_key(fptr, TINT, "NAXIS2", &naxis2, comment, &status);
  fits_read_key(fptr, TINT, "NAXIS3", &naxis3, comment, &status);
  fits_read_key(fptr, TINT, "NAXIS4", &naxis4, comment, &status);
  
  fits_read_key(fptr, TFLOAT, "ANGVAL1", &angval1, comment, &status);
  fits_read_key(fptr, TFLOAT, "ANGVAL2", &angval2, comment, &status);
  fits_read_key(fptr, TFLOAT, "ANGVAL3", &angval3, comment, &status);
  fits_read_key(fptr, TFLOAT, "ANGVAL4", &angval4, comment, &status);
  
  fits_read_key(fptr, TFLOAT, "ANGDELT1", &angdelt1, comment, &status);
  fits_read_key(fptr, TFLOAT, "ANGDELT2", &angdelt2, comment, &status);
  fits_read_key(fptr, TFLOAT, "ANGDELT3", &angdelt3, comment, &status);
  fits_read_key(fptr, TFLOAT, "ANGDELT4", &angdelt4, comment, &status);
  
  long nElements = naxis1*naxis2*naxis3*naxis4, fElement = 1;
  float* angularDistribution = new float[nElements];
  float nullValue = 0;
  int anyNull;
    
  fits_read_img(fptr, TFLOAT, fElement, nElements, &nullValue, 
		angularDistribution, &anyNull, &status); 

  //fAngularDistribution = angularDistribution;

  const unsigned long azBins = naxis1, cZBins = naxis2;

  const unsigned long rBins = fRData.size(), zBins = fZData.size(), wlBins = fWavelengthData.size();

  fAzimuthData.resize(azBins);

  for (unsigned long i = 0; i < azBins; ++i)     
    fAzimuthData[i] = angval1 + i*angdelt1;
      
  fCosZenithData.resize(cZBins);

  for (unsigned long i = 0; i < cZBins; ++i)     
    fCosZenithData[i] = angval2 + i*angdelt2;

  const double cosZenithDelta = angdelt2, azimuthDelta = angdelt1;

  //cout << "ISRF Intensity" << " " << cZBins << " " << cosZenithDelta << " " << azimuthDelta << " " << azBins << endl;

  fIntensityTotal.resize(rBins);
  //fIntensityStellar.resize(rBins);
  //fIntensityScattered.resize(rBins);
  //fIntensityInfrared.resize(rBins);

  for (unsigned long i = 0; i < rBins; ++i) {

    fIntensityTotal[i].resize(zBins);
    //fIntensityStellar[i].resize(zBins);
    ///fIntensityScattered[i].resize(zBins);
    //fIntensityInfrared[i].resize(zBins);

    for (unsigned long j = 0; j < zBins; ++j) {

      fIntensityTotal[i][j].resize(wlBins);
      //fIntensityStellar[i][j].resize(wlBins);
      //fIntensityScattered[i][j].resize(wlBins);
      //fIntensityInfrared[i][j].resize(wlBins);

      for (unsigned long k = 0; k < wlBins; ++k) {

	fIntensityTotal[i][j][k].resize(azBins);
	//fIntensityStellar[i][j][k].resize(azBins);
	//fIntensityScattered[i][j][k].resize(azBins);
	//fIntensityInfrared[i][j][k].resize(azBins);

	for (unsigned long l = 0; l < azBins; ++l) {

	  fIntensityTotal[i][j][k][l].resize(cZBins);
	  //fIntensityStellar[i][j][k][l].resize(cZBins);
	  //fIntensityScattered[i][j][k][l].resize(cZBins);
	  //fIntensityInfrared[i][j][k][l].resize(cZBins);

	  fIntensityTotal[i][j][k][l] = 0;
	  //fIntensityStellar[i][j][k][l] = 0;
	  //fIntensityScattered[i][j][k][l] = 0;
	  //fIntensityInfrared[i][j][k][l] = 0;

	  for (unsigned int m = 0; m < cZBins; ++m) {

	    const double solidAngle = 
	      (m == 0 || m == cZBins-1 ? 0.5 : 1)*cosZenithDelta*azimuthDelta;

	    const unsigned long indexStellarIntensity = 
	      0*wlBins*rBins*zBins*azBins*cZBins + 
	      k*rBins*zBins*azBins*cZBins + 
	      j*rBins*azBins*cZBins + 
	      i*azBins*cZBins + l*cZBins + m;

	    //fIntensityStellar[i][j][k][l][m] = 
	    //angularDistribution[indexStellarIntensity]/solidAngle*
	    //fEnergyDensityStellar[i][j][k]*kSpeedOfLight_SI/cm;

	    const unsigned long indexScatteredIntensity = 
	      1*wlBins*rBins*zBins*azBins*cZBins + 
	      k*rBins*zBins*azBins*cZBins + 
	      j*rBins*azBins*cZBins + 
	      i*azBins*cZBins + l*cZBins + m;

	    //fIntensityScattered[i][j][k][l][m] = 
	    //angularDistribution[indexScatteredIntensity]/solidAngle*
	    //fEnergyDensityScattered[i][j][k]*kSpeedOfLight_SI/cm;

	    const unsigned long indexInfraredIntensity = 
	      2*wlBins*rBins*zBins*azBins*cZBins + 
	      k*rBins*zBins*azBins*cZBins + 
	      j*rBins*azBins*cZBins + 
	      i*azBins*cZBins + l*cZBins + m;
	    
	    //fIntensityInfrared[i][j][k][l][m] = 
	    //angularDistribution[indexInfraredIntensity]/solidAngle*
	    //fEnergyDensityInfrared[i][j][k]*kSpeedOfLight_SI/cm;

	    fIntensityTotal[i][j][k][l][m] = 
	      (angularDistribution[indexStellarIntensity]*
	       fEnergyDensityStellar[i][j][k] +
	       angularDistribution[indexScatteredIntensity]*
	       fEnergyDensityScattered[i][j][k] + 
	       angularDistribution[indexInfraredIntensity]*
	       fEnergyDensityInfrared[i][j][k])/solidAngle*kSpeedOfLight_SI/cm;
	      
	  }

	}

      }

    }

  }

  fits_close_file(fptr, &status);

  delete[] angularDistribution;
    
}

double GalacticRadiationField::
GetEnergyDensity(const double wl,
		 const double r,
		 const double z,
		 const Component comp) const {

  if (comp == CMB) {

    const double energy =
      kPlanck/s*kSpeedOfLight_SI/(wl*micron/m);

    return BlackBodyEnergyDensity(energy, gKT);
   
  }

  const unsigned long wlBins = fWavelengthData.size();
  const unsigned long rBins = fRData.size();
  const unsigned long zBins = fZData.size();

  //cout << wl << endl;

  if (wl < fWavelengthData[0] || wl > fWavelengthData[wlBins-1])
    return 0;
  
  //const double r = sqrt(x*x + y*y);

  const unsigned long rBin = 
    long(r/(fRData[rBins-1] - fRData[0])*rBins);

  const unsigned long zBin = 
    long((z - fZData[0])/(fZData[zBins-1] - fZData[0])*zBins);

  const unsigned long wlBin = 
    long((log(wl) - log(fWavelengthData[0]))/
	 (log(fWavelengthData[wlBins-1]) - log(fWavelengthData[0]))*wlBins);

  //  cout << r << " " << rBin << " " << z << " " << zBin << " " << rBins << " " << zBins << endl;

  // Bi-linear interpolation in radius and z except at boundaries

  double result = 0;

  if (rBin < rBins - 1 && zBin < zBins - 1 && wlBin < wlBins) {

    const double rF = (r - fRData[rBin])/(fRData[rBin+1] - fRData[rBin]);
    const double zF = (z - fZData[zBin])/(fZData[zBin+1] - fZData[zBin]);

    const valarray<double> &stellar1 = fEnergyDensityStellar[rBin][zBin];
    const valarray<double> &stellar2 = fEnergyDensityStellar[rBin+1][zBin];
    const valarray<double> &stellar3 = fEnergyDensityStellar[rBin][zBin+1];
    const valarray<double> &stellar4 = fEnergyDensityStellar[rBin+1][zBin+1];

    valarray<double> zStellar1(wlBins), zStellar2(wlBins);

    zStellar1 = stellar1 + (stellar3 - stellar1)*zF;
    zStellar2 = stellar2 + (stellar4 - stellar2)*zF;

    valarray<double> rStellar(wlBins);

    rStellar = zStellar1 + (zStellar2 - zStellar1)*rF;
    
    //    cout << rBin << " " << zBin << " " << wlBin << " " << fWavelengthData[0] << " " << fWavelengthData[wlBins-1] << " " << wl << " " << stellar1[wlBin] << endl;

    const valarray<double> &scattered1 = fEnergyDensityScattered[rBin][zBin];
    const valarray<double> &scattered2 = fEnergyDensityScattered[rBin+1][zBin];
    const valarray<double> &scattered3 = fEnergyDensityScattered[rBin][zBin+1];
    const valarray<double> &scattered4 = fEnergyDensityScattered[rBin+1][zBin+1];

    valarray<double> zScattered1(wlBins), zScattered2(wlBins);

    zScattered1 = scattered1 + (scattered3 - scattered1)*zF;
    zScattered2 = scattered2 + (scattered4 - scattered2)*zF;

    valarray<double> rScattered(wlBins);

    rScattered = zScattered1 + (zScattered2 - zScattered1)*rF;
    

    const valarray<double> &infrared1 = fEnergyDensityInfrared[rBin][zBin];
    const valarray<double> &infrared2 = fEnergyDensityInfrared[rBin+1][zBin];
    const valarray<double> &infrared3 = fEnergyDensityInfrared[rBin][zBin+1];
    const valarray<double> &infrared4 = fEnergyDensityInfrared[rBin+1][zBin+1];

    valarray<double> zInfrared1(wlBins), zInfrared2(wlBins);

    zInfrared1 = infrared1 + (infrared3 - infrared1)*zF;
    zInfrared2 = infrared2 + (infrared4 - infrared2)*zF;

    valarray<double> rInfrared(wlBins);

    rInfrared = zInfrared1 + (zInfrared2 - zInfrared1)*rF; 
    

    result = (comp == STELLAR ? rStellar[wlBin] : 
	      (comp == SCATTERED ? rScattered[wlBin] : 
	       (comp == INFRARED ? rInfrared[wlBin] : 
		(comp == TOTAL ? 
		 rStellar[wlBin] + rScattered[wlBin] + rInfrared[wlBin] : 0))));
       
  }

  return result;

}

double GalacticRadiationField::
GetIntensity(const double wl,
	     const double r,
	     const double z,
	     const double azimuth,
	     const double cosZenith,
	     const Component comp) const {

  if (comp == CMB) {

    const double energy =
      kPlanck/s*kSpeedOfLight_SI/(wl*micron/m);

    const double cmb = BlackBodyEnergyDensity(energy, gKT)*
      kSpeedOfLight_SI*m/cm*kOneOnFourPi;
  
    return cmb;

  }

  if (comp == STELLAR || comp == SCATTERED || comp == INFRARED)
    return 0.;

  const unsigned long wlBins = fWavelengthData.size();
  const unsigned long rBins = fRData.size();
  const unsigned long zBins = fZData.size();
  const unsigned long azBins = fAzimuthData.size();
  const unsigned long cZBins = fCosZenithData.size();

  //cout << wl << endl;

  if (wl < fWavelengthData[0] || wl > fWavelengthData[wlBins-1])
    return 0;
  
  //const double r = sqrt(x*x + y*y);

  const unsigned long rBin = 
    long(r/(fRData[rBins-1] - fRData[0])*rBins);

  const unsigned long zBin = 
    long((z - fZData[0])/(fZData[zBins-1] - fZData[0])*zBins);

  const unsigned long wlBin = 
    long((log(wl) - log(fWavelengthData[0]))/
	 (log(fWavelengthData[wlBins-1]) - log(fWavelengthData[0]))*wlBins);

  const unsigned long azBin = 
    long(fmod(azimuth, kTwoPi)/(fAzimuthData[azBins-1] - fAzimuthData[0])*(azBins-1));

  const unsigned long cZBin =
    (cosZenith <= -1 ? 0 :
     (cosZenith >= 1 ? cZBins - 1 :
      long((cosZenith + 1.)/(fCosZenithData[cZBins-1] - fCosZenithData[0])*cZBins)));
 
  // Bi-linear interpolation in radius and z except at boundaries

  const double cosZenithDelta = fabs(fCosZenithData[1] - fCosZenithData[0]);

  const double azimuthDelta = fabs(fAzimuthData[1] - fAzimuthData[0]);

  const double solidAngle = 
    (cZBin == 0 || cZBin == cZBins-1 ? 0.5 : 1)*cosZenithDelta*azimuthDelta;

  //cout << cosZenith << " " << cZBin << " " << azimuth << " " << azBin << endl;

  double result = 0;

  if (rBin < rBins - 1 && 
      zBin < zBins - 1 && 
      wlBin < wlBins && 
      azBin < azBins && 
      cZBin < cZBins) {

    const double rF = (r - fRData[rBin])/(fRData[rBin+1] - fRData[rBin]);
    const double zF = (z - fZData[zBin])/(fZData[zBin+1] - fZData[zBin]);

    //const double stellar1 = fIntensityStellar[rBin][zBin][wlBin][azBin][cZBin];
    //const double stellar2 = fIntensityStellar[rBin][zBin][wlBin][azBin+1][cZBin];
    //const double stellar3 = fIntensityStellar[rBin][zBin][wlBin][azBin][cZBin+1];
    //const double stellar4 = fIntensityStellar[rBin][zBin][wlBin][azBin+1][cZBin+1];

    //const double cStellar1 = stellar1;// + (stellar3 - stellar1)*cF;
    //const double cStellar2 = stellar2;// + (stellar4 - stellar2)*cF;

    //const double aStellar = stellar1;//cStellar1 + (cStellar2 - cStellar1)*aF;
    

    //const double scattered1 = fIntensityScattered[rBin][zBin][wlBin][azBin][cZBin];
    //const double scattered2 = fIntensityScattered[rBin][zBin][wlBin][azBin+1][cZBin];
    //const double scattered3 = fIntensityScattered[rBin][zBin][wlBin][azBin][cZBin+1];
    //const double scattered4 = fIntensityScattered[rBin][zBin][wlBin][azBin+1][cZBin+1];

    //const double cScattered1 = scattered1;// + (scattered3 - scattered1)*cF;
    //const double cScattered2 = scattered2;// + (scattered4 - scattered2)*cF;

    //const double aScattered = scattered1;//cScattered1 + (cScattered2 - cScattered1)*aF;
    
    
    //const double infrared1 = fIntensityInfrared[rBin][zBin][wlBin][azBin][cZBin];
    //const double infrared2 = fIntensityInfrared[rBin][zBin][wlBin][azBin+1][cZBin];
    //const double infrared3 = fIntensityInfrared[rBin][zBin][wlBin][azBin][cZBin+1];
    //const double infrared4 = fIntensityInfrared[rBin][zBin][wlBin][azBin+1][cZBin+1];

    //const double cInfrared1 = infrared1;// + (infrared3 - infrared1)*cF;
    //const double cInfrared2 = infrared2;// + (infrared4 - infrared2)*cF;

    //const double aInfrared = infrared1;//cInfrared1 + (cInfrared2 - cInfrared1)*aF;
    
    
    result = fIntensityTotal[rBin][zBin][wlBin][azBin][cZBin];
    //(comp == STELLAR ? aStellar : 
    //      (comp == SCATTERED ? aScattered : 
    //       (comp == INFRARED ? aInfrared : 
    //	(comp == TOTAL ? aStellar + aScattered + aInfrared : 0))));
	      
  }
  
  return result;

}

// Energy in eV, kT in units of eV. Result returned as eV cm^-3

double 
GalacticRadiationField::BlackBodyEnergyDensity(const double energy, 
					       const double kT) const {

  const double energyOnPi = energy/kPi;

  return pow(kPlanckReduced/s/eV*kSpeedOfLight_SI/cm, -3.0)*
    energyOnPi*energyOnPi/(exp(energy/kT) - 1.0)*energy*energy;

}
