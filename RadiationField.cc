#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <PhysicalConstants.h>
#include <ErrorLogger.h>

#include <RadiationField.h>

#include <CCfits/CCfits>

//#include <CLHEP/Vector/ThreeVector.h>

using namespace std;
using namespace rf;
using namespace utl;
//using namespace CLHEP;

RadiationField::RadiationField() : 
  fRebinnedSkymapOrder(1), 
  fNumberOfComponents(0), 
  fLuminosity(0), 
  fDustMass(0), 
  fCacheBuilt(false) {

}

RadiationField::RadiationField(const string& filename,
			       const valarray<double>& freq,
			       int rebinnedSkymapOrder) :
  fRebinnedSkymapOrder(rebinnedSkymapOrder),
  fNumberOfComponents(0),
  fLuminosity(0),
  fDustMass(0),
  fCacheBuilt(false) {

  assert(!filename.empty());

  fFrequency.resize(freq.size());
  fFrequency = freq;
  
  ReadRadiationField(filename);

}

RadiationField::~RadiationField() {

  ClearData();

}

const Skymap<double>
RadiationField::GetSkymap(const ThreeVector& pos,
			  const STELLARCOMPONENT component,
			  const int healpixOrder) {

  assert (component >= TOTAL && component <= THERMAL);

  if (fRebinnedSkymapOrder == healpixOrder)
    return GetSkymap(pos, component);
  else {

    //FlushSkymapCache(component);

    //const int rebinnedSkymapOrder = fRebinnedSkymapOrder;

    //fRebinnedSkymapOrder = healpixOrder;

    Skymap<double> skymapTmp = GetSkymap(pos, component);

    Skymap<double> skymap = skymapTmp.rebin(healpixOrder);

    //FlushSkymapCache(component);

    //fRebinnedSkymapOrder = rebinnedSkymapOrder;

    return skymap;

  }

}

const Skymap<double> 
RadiationField::GetSkymap(const ThreeVector& pos, 
			  const STELLARCOMPONENT component) {

  //cout << "Cache built: " << fCacheBuilt[component] << endl;

  if (!fCacheBuilt)//[component])
    BuildSkymapCache();//component);

  const double r = sqrt(pos.x()*pos.x() + pos.y()*pos.y()), z = pos.z(), absZ = fabs(z);

  if (r > fRData[fRData.size()-1] || absZ > fZData[fZData.size()-1]) {

    const Skymap<double>& skymap = *fSkymapOrderedData[0][0][component];

    return ((Skymap<double>() = skymap) = 0); 

  }

  // Bilinear interpolation to form skymap

  const double* rVal = std::lower_bound(&fRData[0], &fRData[fRData.size()-1], r);

  const int rIndex = rVal - &fRData[0];

  double rCoeff = 0, dR = 0;

  if (r <= fRData[rIndex] && r >= fRData[rIndex])
    rCoeff = 1;
  else if (r < fRData[fRData.size()-1]) {
    
    dR = fRData[rIndex] - fRData[rIndex-1];
    
    rCoeff = (fRData[rIndex] - r)/dR;
    
  }

  const double* zVal = std::lower_bound(&fZData[0], &fZData[fZData.size()-1], absZ);

  const int zIndex = zVal - &fZData[0];
 
  double zCoeff = 0, dZ = 0;

  if (absZ <= fZData[zIndex] && absZ >= fZData[zIndex]) 
    zCoeff = 1;
  else if (absZ < fZData[fZData.size()-1]) {
    
    dZ = fZData[zIndex] - fZData[zIndex-1];
    
    zCoeff = (fZData[zIndex] - absZ)/dZ;
    
  }
  
  // Form the bilinear interpolated skymap 

  Skymap<double> skymap;
  
  if (rCoeff >= 1. && zCoeff >= 1.) {

    if (z >= 0.)
      return *fSkymapOrderedData[rIndex][zIndex][component];
    else {

      const Skymap<double>& skymap1 = *fSkymapOrderedData[rIndex][zIndex][component];

      Skymap<double> skymap1Mirror = skymap1;

      for (unsigned long i = 0; i < skymap1.Npix(); ++i) {

	SM::Coordinate coord = skymap1.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap1Mirror[mirrorCoord] = skymap1[coord];

      }

      return skymap1Mirror;

    }

  } else if (rCoeff >= 1.) {

    const Skymap<double>& skymap1 = *fSkymapOrderedData[rIndex][zIndex-1][component];
    
    const Skymap<double>& skymap2 = *fSkymapOrderedData[rIndex][zIndex][component];

    if (z >= 0) 
      skymap = skymap1*zCoeff + skymap2*(1. - zCoeff);
    else {

      Skymap<double> skymap1Mirror = skymap1;
      Skymap<double> skymap2Mirror = skymap2;

      for (unsigned long i = 0; i < skymap1.Npix(); ++i) {

	SM::Coordinate coord = skymap1.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap1Mirror[mirrorCoord] = skymap1[coord];

      }

      for (unsigned long i = 0; i < skymap2.Npix(); ++i) {

	SM::Coordinate coord = skymap2.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap2Mirror[mirrorCoord] = skymap2[coord];

      }

      skymap = skymap1Mirror*zCoeff + skymap2Mirror*(1. - zCoeff); 

    }

  } else if (zCoeff >= 1.) {

    const Skymap<double>& skymap1 = *fSkymapOrderedData[rIndex-1][zIndex][component];

    const Skymap<double>& skymap2 = *fSkymapOrderedData[rIndex][zIndex][component];

    if (z >= 0.)
      skymap = skymap1*rCoeff + skymap2*(1. - rCoeff);
    else {

      Skymap<double> skymap1Mirror = skymap1;
      Skymap<double> skymap2Mirror = skymap2;

      for (unsigned long i = 0; i < skymap1.Npix(); ++i) {

	SM::Coordinate coord = skymap1.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap1Mirror[mirrorCoord] = skymap1[coord];

      }

      for (unsigned long i = 0; i < skymap2.Npix(); ++i) {

	SM::Coordinate coord = skymap2.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap2Mirror[mirrorCoord] = skymap2[coord];

      }

      skymap = skymap1Mirror*rCoeff + skymap2Mirror*(1. - rCoeff);

    }

  } else {

    const Skymap<double>& skymap1 = *fSkymapOrderedData[rIndex-1][zIndex-1][component];

    const Skymap<double>& skymap2 = *fSkymapOrderedData[rIndex-1][zIndex][component];

    const Skymap<double>& skymap3 = *fSkymapOrderedData[rIndex][zIndex-1][component];

    const Skymap<double>& skymap4 = *fSkymapOrderedData[rIndex][zIndex][component];

    if (z >= 0.) 
      skymap = skymap1*rCoeff*zCoeff + skymap2*(1. - zCoeff)*rCoeff + skymap3*zCoeff*(1. - rCoeff) + skymap4*(1. - rCoeff)*(1. - zCoeff);
    else {

      Skymap<double> skymap1Mirror = skymap1;
      Skymap<double> skymap2Mirror = skymap2;
      Skymap<double> skymap3Mirror = skymap3;
      Skymap<double> skymap4Mirror = skymap4;

      for (unsigned long i = 0; i < skymap1.Npix(); ++i) {

	SM::Coordinate coord = skymap1.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap1Mirror[mirrorCoord] = skymap1[coord];

      }

      for (unsigned long i = 0; i < skymap2.Npix(); ++i) {

	SM::Coordinate coord = skymap2.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap2Mirror[mirrorCoord] = skymap2[coord];

      }

      for (unsigned long i = 0; i < skymap3.Npix(); ++i) {

	SM::Coordinate coord = skymap3.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap3Mirror[mirrorCoord] = skymap3[coord];

      }

      for (unsigned long i = 0; i < skymap4.Npix(); ++i) {

	SM::Coordinate coord = skymap4.pix2coord(i);

	SM::Coordinate mirrorCoord(coord.l(), -coord.b());

	skymap4Mirror[mirrorCoord] = skymap4[coord];

      }
      
      skymap = skymap1Mirror*rCoeff*zCoeff + skymap2Mirror*(1. - zCoeff)*rCoeff + skymap3Mirror*zCoeff*(1. - rCoeff) + skymap4Mirror*(1. - rCoeff)*(1. - zCoeff); 

    }

  }

  return skymap;

}

const valarray<double>
RadiationField::GetNumberDensity(const ThreeVector& pos, 
				 const STELLARCOMPONENT component) {

  valarray<double> spec(0., fFrequency.size());

  const double r = sqrt(pos.x()*pos.x() + pos.y()*pos.y()), z = pos.z(), absZ = fabs(z);
  
  if (r > fRData[fRData.size()-1] || z > fZData[fZData.size()-1]) {

    return spec;

  }

  // Bilinear interpolation 

  const double* rVal = std::lower_bound(&fRData[0], &fRData[fRData.size()-1], r);

  const int rIndex = rVal - &fRData[0];

  double rCoeff = 0, dR = 0;

  if (r <= fRData[rIndex] && r >= fRData[rIndex])
    rCoeff = 1;
  else if (r < fRData[fRData.size()-1]) {
    
    dR = fRData[rIndex] - fRData[rIndex-1];
    
    rCoeff = (fRData[rIndex] - r)/dR;
    
  }

  const double* zVal = std::lower_bound(&fZData[0], &fZData[fZData.size()-1], absZ);

  const int zIndex = zVal - &fZData[0];
 
  double zCoeff = 0, dZ = 0;

  if (absZ <= fZData[zIndex] && absZ >= fZData[zIndex]) 
    zCoeff = 1;
  else if (absZ < fZData[fZData.size()-1]) {
    
    dZ = fZData[zIndex] - fZData[zIndex-1];
    
    zCoeff = (fZData[zIndex] - absZ)/dZ;
    
  }
  
  //cout << z << " " << zIndex << " " << dZ << " " << zCoeff << " " << fZData[zIndex] << endl;

  //exit(0);

  const valarray<double>& wl = fWavelength;

  valarray<double> freq(0., wl.size());
      
  for (unsigned int iWl = 0; iWl < wl.size(); ++iWl)
    freq[iWl] = kSpeedOfLight_SI*1./(wl[wl.size()-1-iWl]*micron/m);
  
  valarray<double> energyDensity(0., wl.size());  

  if (rCoeff >= 1. && zCoeff >= 1.) {

    energyDensity = *fEnergyDensity[rIndex][zIndex][component];

  } else if (rCoeff >= 1.) {

    const valarray<double>& energyDensity1 = *fEnergyDensity[rIndex][zIndex-1][component];
    
    const valarray<double>& energyDensity2 = *fEnergyDensity[rIndex][zIndex][component];

    energyDensity = energyDensity1*zCoeff + energyDensity2*(1. - zCoeff);

  } else if (zCoeff >= 1.) {

    const valarray<double>& energyDensity1 = *fEnergyDensity[rIndex-1][zIndex][component];

    const valarray<double>& energyDensity2 = *fEnergyDensity[rIndex][zIndex][component];

    energyDensity = energyDensity1*rCoeff + energyDensity2*(1. - rCoeff);

  } else {

    const valarray<double>& energyDensity1 = *fEnergyDensity[rIndex-1][zIndex-1][component];

    const valarray<double>& energyDensity2 = *fEnergyDensity[rIndex-1][zIndex][component];

    const valarray<double>& energyDensity3 = *fEnergyDensity[rIndex][zIndex-1][component];

    const valarray<double>& energyDensity4 = *fEnergyDensity[rIndex][zIndex][component];

    energyDensity = energyDensity1*rCoeff*zCoeff + energyDensity2*(1. - zCoeff)*rCoeff + energyDensity3*zCoeff*(1. - rCoeff) + energyDensity4*(1. - rCoeff)*(1. - zCoeff);

  }

  valarray<double> binCount(0., fFrequency.size()), energy(fFrequency.size());
  
  energy = kPlanck_SI/e_SI*fFrequency;
 
  for (unsigned int iFreq = 0; iFreq < freq.size(); ++iFreq) {
	  
    const int index = (log10(freq[iFreq]) - log10(fFrequency[0]))/(log10(fFrequency[fFrequency.size()-1]) - log10(fFrequency[0]))*fFrequency.size();
    
    if (index >= 0 && index < fFrequency.size()) {
	    
      const double specVal = energyDensity[freq.size()-1-iFreq];
	    
      spec[index] += specVal;
      binCount[index] += 1.;
	    
    }
	  
  }
	
  spec *= 1./binCount*1./energy*1./energy;

  return spec;

}

void RadiationField::ReadRadiationField(const string& filename) {

  //#pragma omp critical 
  {

    if (filename.empty()) {

      INFO("Empty radiation field filename passed to ReadRadiationField.");
      exit(-1);

    }
  
    ostringstream buf;
    buf << "Reading from " << filename;
    INFO(buf.str());
    
    ClearData();
    
    const string prefix = filename.substr(0, filename.find_last_of("/")+1);
    
    //cout << prefix << endl;
    
    fPrefix = prefix;
    
    ifstream rfFile(filename.c_str());

    unsigned long volumeElements, wlBins, numFilters;
    double luminosity, dustMass;
    
    rfFile >> volumeElements >> wlBins >> numFilters;
    rfFile >> luminosity >> dustMass;
    
    if (!rfFile) {
      buf.str("");
      buf << "Error reading from file " << filename;
      ERROR(buf.str());
      throw(std::invalid_argument(buf.str()));
    }
    
    fLuminosity = luminosity;
    fDustMass = dustMass;
    
    fWavelength.resize(wlBins);
    fWavelength = 0;
    
    unsigned long stellarComponents;
    
    rfFile >> stellarComponents;
    
    fStellarComponentLuminosity.resize(stellarComponents);
    
    for (unsigned long i = 0; i < stellarComponents; ++i) {
      
      double componentLuminosity;
      string componentName;
      
      rfFile >> componentName >> componentLuminosity;
      
      fStellarComponentLuminosity[i] = componentLuminosity;
      fStellarComponentName.push_back(componentName);
      
    }
    
    unsigned long geometry;
    double modelRegionData[6];
    
    rfFile >> geometry;
    rfFile >> modelRegionData[0] >> modelRegionData[1] >> modelRegionData[2] >> modelRegionData[3] >> modelRegionData[4] >> modelRegionData[5];
    
    vector<ThreeVector> posVec, rangeVec;
    
    for (unsigned long i = 0; i < volumeElements; ++i) {
      
      double index;
      double x, y, z, dX, dY, dZ;
      
      rfFile >> index;
      rfFile >> x >> y >> z;
      rfFile >> dX >> dY >> dZ;
      
      ThreeVector pos(x, y, z), range(dX, dY, dZ);
      
      //cout << pos << endl;
      
      if (!posVec.size()) {
	
	posVec.push_back(pos);
	rangeVec.push_back(range);
	
      } else if (fabs(posVec[0].x() - pos.x()) < dX) {
	
	posVec.push_back(pos);
	rangeVec.push_back(range);
	
      }
      
      if (fabs(posVec[0].x() - pos.x()) > dX || i == volumeElements - 1) {
	
	fPositionData.push_back(posVec);
	fRangeData.push_back(rangeVec);
	
	posVec.clear();
	posVec.push_back(pos);
	
	rangeVec.clear();
	rangeVec.push_back(range);
	
      }
      
      string filterFilename, filterCountFilename, directFilename, scatteredFilename, transientFilename, thermalFilename, totalFilename, filterFluxFilename, fluxFilename;
      
      rfFile >> filterFilename;
      rfFile >> filterCountFilename;
      rfFile >> directFilename;
      rfFile >> scatteredFilename;
      rfFile >> transientFilename;
      rfFile >> thermalFilename;
      rfFile >> totalFilename;
      rfFile >> filterFluxFilename;
      rfFile >> fluxFilename;
      
      fNumberOfComponents = 7;
      
      //fCacheBuilt.resize(fNumberOfComponents);
      //fCacheBuilt = false;
      
      fFilenameData[pos].resize(fNumberOfComponents+1);
      
      fFilenameData[pos][0] = totalFilename;
      fFilenameData[pos][1] = directFilename;
      fFilenameData[pos][2] = scatteredFilename;
      fFilenameData[pos][3] = transientFilename;
      fFilenameData[pos][4] = thermalFilename;
      fFilenameData[pos][5] = filterFilename;
      fFilenameData[pos][6] = filterCountFilename;
      fFilenameData[pos][7] = fluxFilename;
      
    }
    
    for (unsigned int i = 0; i < fPositionData.size(); ++i) {
      
      vector< vector<string>* > pStrData;
      
      vector< vector<Skymap<double>* > > pSkymapData;
      
      vector< vector< valarray<double>* > > pEnergyDensityData;
      
      for (unsigned int j = 0; j < fPositionData[i].size(); ++j) {
	
	const ThreeVector& pos = fPositionData[i][j];
	
	map<ThreeVector, vector<string> >::iterator fIt = fFilenameData.find(pos);
	
	pStrData.push_back(&fIt->second);
	
	vector< Skymap<double>* > skymaps;
	skymaps.resize(pStrData.back()->size()); 
	skymaps.assign(pStrData.back()->size(), 0);
	
	pSkymapData.push_back(skymaps);
	
	// Read in the flux file at this stage -- we're ready to go at exit from
	// this method if this is all that the user will be requesting.
	
	vector< valarray<double>* > energyDensity;
	
	for (int iC = 0; iC < 5; ++iC) {
	  
	  valarray<double>* vec = new valarray<double>(0., wlBins);
	  energyDensity.push_back(vec);
	  
	}
	
	const string fname = prefix + "/" + pStrData.back()->back();
	
	CCfits::FITS fluxFile(fname, CCfits::Read);
	
	CCfits::ExtHDU& table = fluxFile.currentExtension(); // It should be the only one in this file ...
	
	vector<double> wl, total, direct, scattered, transient, thermal;
	
	table.column(1).read(wl, 1, wlBins);
	table.column(2).read(total, 1, wlBins);
	table.column(3).read(direct, 1, wlBins);
	table.column(4).read(scattered, 1, wlBins);
	table.column(5).read(transient, 1, wlBins);
	table.column(6).read(thermal, 1, wlBins);
	
	for (int iWl = 0; iWl < wlBins; ++iWl) {
	  
	  fWavelength[iWl] = wl[iWl];
	  
	  (*(energyDensity[0]))[iWl] = total[iWl];
	  (*(energyDensity[1]))[iWl] = direct[iWl];
	  (*(energyDensity[2]))[iWl] = scattered[iWl];
	  (*(energyDensity[3]))[iWl] = transient[iWl];
	  (*(energyDensity[4]))[iWl] = thermal[iWl];
	  
	  //cout << i << " " << j << " " << iWl << " " << wl[iWl] << " " << (*(energyDensity[0]))[iWl] << " " << (*(energyDensity[4]))[iWl] << endl;
	  
	}
	
	//exit(0);
	
	pEnergyDensityData.push_back(energyDensity);
	
      }
      
      fFilenameOrderedData.push_back(pStrData);
      fSkymapOrderedData.push_back(pSkymapData);
      fEnergyDensity.push_back(pEnergyDensityData);
      
    }
    
    //exit(0);
    
    fRData.resize(fPositionData.size());
    
    for (unsigned int i = 0; i < fPositionData.size(); ++i) {
      
      fRData[i] = fPositionData[i][0].x();
            
    }
    
    fZData.resize(fPositionData[0].size());
    
    for (unsigned int i = 0; i < fPositionData[0].size(); ++i)
      fZData[i] = fPositionData[0][i].z();
    
    fRRangeData.resize(fRangeData.size());
        
    for (unsigned int i = 0; i < fRangeData.size(); ++i) {
      
      fRRangeData[i] = fRangeData[i][0].x();
            
    }
    
    fZRangeData.resize(fRangeData[0].size());
    
    for (unsigned int i = 0; i < fRangeData[0].size(); ++i)
      fZRangeData[i] = fRangeData[0][i].z();
    
  }
    
}

void RadiationField::BuildSkymapCache() { //const unsigned int component) {

  //#pragma omp critical
  {

    INFO("Entry");

    if (!fCacheBuilt) {
        
      ostringstream buf;
      buf << "Building skymap cache";// for component" << component;
      INFO(buf.str());
    
    //assert(component < fNumberOfComponents);
      
      valarray<double> energy(fFrequency.size());
      
      energy = kPlanck_SI/e_SI*fFrequency;
      
      for (int i = 0; i < fSkymapOrderedData.size(); ++i) {
	
	for (int j = 0; j < fSkymapOrderedData[i].size(); ++j) {
	  
	  for (int k = RadiationField::TOTAL; k <= RadiationField::THERMAL; ++k) {
	    
	    const string& filename = (*fFilenameOrderedData[i][j])[k];
	    
	    //cout << filename << endl;
	    
	    Skymap<double> skymap(fPrefix + filename);
	    
	    const valarray<double>& wl = skymap.getSpectra();
	    
	    valarray<double> freq(0., wl.size()), en(0., wl.size());
	    
	    for (unsigned int iWl = 0; iWl < wl.size(); ++iWl)
	      freq[iWl] = kSpeedOfLight_SI*1./(wl[wl.size()-1-iWl]*micron/m);
	    
	    en = kPlanck_SI/e_SI*freq;
	    
	    Skymap<double> skymapRebinned(fRebinnedSkymapOrder, freq);
	    
	    skymapRebinned = skymap.rebin(fRebinnedSkymapOrder);
	    
	    fSkymapOrderedData[i][j][k] = new Skymap<double>(fRebinnedSkymapOrder, fFrequency);	
	    
	    for (unsigned int nPix = 0; nPix < skymapRebinned.Npix(); ++nPix) {
	      
	      valarray<double> spec(0., fFrequency.size()), binCount(0., fFrequency.size());
	      
	      for (unsigned int iFreq = 0; iFreq < freq.size(); ++iFreq) {
		
		const int index = (log10(freq[iFreq]) - log10(fFrequency[0]))/(log10(fFrequency[fFrequency.size()-1]) - log10(fFrequency[0]))*fFrequency.size();
		
		if (index >= 0 && index < fFrequency.size()) {
		  
		  const double specVal = skymapRebinned[nPix][freq.size()-1-iFreq];//energyRaw[iFreq]/energyRaw[iFreq];
		  
		  spec[index] += specVal;
		  binCount[index] += 1.;
		  
		}
		
	      }
	      
	      spec *= 1./binCount*1./energy*1./energy*1./(kSpeedOfLight_SI*m/cm);
	      
	      (*fSkymapOrderedData[i][j][k])[nPix] = spec; // eV^-1 cm^-3 -- per pixel
	      
	    }
	    
	    /*for (size_t iFreq = 0; iFreq < freq.size(); ++iFreq) {
	      
	      size_t index = (log10(freq[iFreq]) - log10(fFrequency[0]))/(log10(fFrequency[fFrequency.size()-1]) - log10(fFrequency[0]))*fFrequency.size();
	      
	      cout << iFreq << " " 
	      << freq[iFreq] << " " 
	      << en[iFreq] << " " 
	      << index << " " 
	      << (index >= 0 && index < fFrequency.size() ? fFrequency[index] : 0) << " "
	      << (index >= 0 && index < fFrequency.size() ? energy[index] : 0) << " " 
	      << skymap.sum(freq.size()-1-iFreq)*1./(kSpeedOfLight_SI*m/cm) << " " 
	      << skymapRebinned.sum(freq.size()-1-iFreq)*1./(kSpeedOfLight_SI*m/cm) << " "
	      << (index >= 0 && index < fFrequency.size() ? fSkymapOrderedData[i][j][component]->sum(index)*energy[index]*energy[index] : 0) << " "
	      << fFrequency.size() << endl;
	      
	      }
	      
	      exit(0);
	    */
	  }
	  
	}     
	
      }
      
      fCacheBuilt = true;//[component] = true;

    } else {

      INFO("Cache already built.");

    }

    INFO("Exit");

  }

}
  
void RadiationField::FlushSkymapCache() {//const unsigned int component) {

  INFO("Entry");

  //assert(component < fNumberOfComponents);
  
  for (unsigned int i = 0; i < fSkymapOrderedData.size(); ++i) {
    
    for (unsigned int j = 0; j < fSkymapOrderedData[i].size(); ++j) {
      
      for (unsigned int k = 0; k < fNumberOfComponents; ++k) {

	delete fSkymapOrderedData[i][j][k];
	fSkymapOrderedData[i][j][k] = 0;
      
      }

    }
    
  }

  fCacheBuilt = false;//[component] = false;
  
  INFO("Exit");

}

/*void RadiationField::FlushSkymapCache() {

  for (unsigned int i = 0; i < fNumberOfComponents; ++i)
    FlushSkymapCache(i);
  
}
*/

void RadiationField::ClearData() {

  //#pragma omp critical
  {

    if (fCacheBuilt) {
      
      for (unsigned int i = 0; i < fPositionData.size(); ++i) {
	
	fPositionData[i].clear();
	fRangeData[i].clear();
	
      }
      
      for (unsigned int i = 0; i < fFilenameOrderedData.size(); ++i)
	fFilenameOrderedData[i].clear(); // Don't own the pointers
      
      fFilenameOrderedData.clear();
      
      for (unsigned int i = 0; i < fSkymapOrderedData.size(); ++i) {
	
	for (unsigned int j = 0; j < fSkymapOrderedData[i].size(); ++j) {
	  
	  for (unsigned int k = 0; k < fSkymapOrderedData[i][j].size(); ++k) 
	    delete fSkymapOrderedData[i][j][k];
	  
	  fSkymapOrderedData[i][j].clear();
	  
	  for (unsigned int k = 0; k < fEnergyDensity[i][j].size(); ++k) 
	    delete fEnergyDensity[i][j][k];
	  
	  fEnergyDensity[i][j].clear();
	  
	}
	
	fSkymapOrderedData[i].clear();
	fEnergyDensity[i].clear();
	
      }
      
      fSkymapOrderedData.clear();
      fEnergyDensity.clear();
      
      fCacheBuilt = false;
      
    } else {

      INFO("Call to clear cache when not built");

    }

  }

}
