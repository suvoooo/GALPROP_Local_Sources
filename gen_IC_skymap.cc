
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_IC_skymap.cc *                            galprop package * 4/20/2006
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// generate inverse Compton skymaps 

// list of selectable debugs:  IMOS20080114
// galdef.verbose==-457 is the return to the old method of integration to compare with older versions

#include"galprop_classes.h"
#include"galprop_internal.h"

#include <cassert>
#include <sstream>
#include <valarray>
#include <vector>

#include <ErrorLogger.h>
#include <PhysicalConstants.h>
#include <RadiationField.h>
#include <Skymap.h>

#include <config.h>

using namespace std;

// Inverse Compton source function assuming isotropic electrons and targets
// Units : photons MeV^-1 cm^2 electron^-1 photon^-1

static double IsotropicICSSourceFunction(const double epsilon2,
					 const double epsilon1,
					 const double gamma) {

  using namespace utl;

  double result = 0.;

  const double p = 4.*epsilon1*gamma;

  if (gamma > 0. &&
      gamma > epsilon2 &&
      epsilon2 <= p*gamma/(1. + p)) {

    const double q = epsilon2/p/(gamma - epsilon2);

    const double pq = p*q;

    const double F =
     2.*q*log(q) + (1. + 2.*q)*(1. - q) + 0.5*pq*pq/(1. + pq)*(1. - q);

    result =
      3./4.*(kThomsonCrossSection_SI*utl::m/utl::cm*utl::m/utl::cm)*1./kElectronMass*1./(utl::eV/utl::MeV)*F/gamma/gamma/epsilon1;

    //cout << "Iso: " << result*kElectronMass << " " << kElectronMass << " " << kThomsonCrossSection_SI*utl::m/utl::cm*utl::m/utl::cm << endl;

  }

  return result;

}

// General inverse Compton source function assuming isotropic electrons
// Units : photons MeV^-1 cm^2 electron^-1 photon^-1 

static double AnisotropicICSSourceFunction(const double epsilon2,
					   const double epsilon1,
					   const double gamma,
					   const double cosZeta) {
  
  using namespace utl;

  double result = 0.;
  
  if (gamma > 1. &&
      gamma > epsilon2) {
    
    const double oneOnGamma = 1./gamma;
    
    const double beta = sqrt(1. - oneOnGamma*oneOnGamma);
    
    const double epsilon1Prime = epsilon1*gamma*(1. + beta*cosZeta);
    
    if (epsilon2 <= 2.*gamma*epsilon1Prime/(1. + 2.*epsilon1Prime)) {
      
      const double epsilon2OnGamma = epsilon2*oneOnGamma;
      
      const double oneOnEpsilon1Prime = 1./epsilon1Prime;
      
      const double gammaMinusEpsilon2 = gamma - epsilon2;
      
      const double F = 2. - 2.*epsilon2OnGamma*(oneOnEpsilon1Prime + 2.) +
        epsilon2OnGamma*epsilon2OnGamma*(oneOnEpsilon1Prime*oneOnEpsilon1Prime +
                                         2*oneOnEpsilon1Prime + 3.) -
        epsilon2OnGamma*epsilon2OnGamma*epsilon2OnGamma;
      
      result = 3./8.*(kThomsonCrossSection_SI/utl::cm/utl::cm)*F/epsilon1/
        gammaMinusEpsilon2/gammaMinusEpsilon2/(kElectronMass*utl::eV/utl::MeV);
      
    }
    
  }

  return result;
    
}

static double BlackBodyNumberDensity(const double energy, const double kT) {   

  using namespace utl;

  const double energyOnPi = energy/kPi,
    constant = 1./(kPlanckReduced/s/eV*kSpeedOfLight_SI/cm);

  return constant*constant*constant*
    energyOnPi*energyOnPi/(exp(energy/kT) - 1.); // eV^-1 cm^-3

}

int Galprop::gen_IC_skymap() { //IMOS20080114

  INFO("Entry: gen_IC_skymap");
  
  int stat = 0;

  const unsigned int nComps = galaxy.n_ISRF_components;

  Skymap<double> anisoRatioOpt, anisoRatioIR, anisoRatioCMB;

  const unsigned int nEGammaBins = galaxy.n_E_gammagrid;

  bool anisoValid = false;

  if (galdef.IC_anisotropic) {

    bool haveOMP = false;

#ifdef HAVE_OPENMP
    haveOMP = true;
#endif

    if (galdef.ISRF_filetype < 3 || nComps < 3 || !haveOMP) {

      ostringstream buf;
      buf << "Anisotropic IC calculation only available for version 3 ISRF files and later with all 3 ISRF components (optical, infrared, CMB) with openmp enabled. Only isotropic IC will be calculated.";
      INFO(buf.str());

    } else {

      anisoValid = true;

      const unsigned int internalHealpixOrder = 2; // Set this to whatever required for the internal calculation skymaps

      const int cosThetaBins = 500; // Set this to the number for precomputing the anisotropic cross section. Accuracy is best for 500, 5000 (need to follow up once working).

      const double dCosTheta = 2./cosThetaBins;

      Distribution electrons;

      // identify the electrons/positrons IMOS20030217
      if (2 == galdef.n_spatial_dimensions) 
	electrons.init(gcr[0].n_rgrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
      
      if (3 == galdef.n_spatial_dimensions) 
	electrons.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
      
      electrons = 0.;
      
      int iE, ielectrons = -1;
      
      for (iE = 0, ielectrons = -1; iE < n_species; ++iE)  
	if (100 == 100*abs(gcr[iE].Z) + gcr[iE].A) {
	  
	  ielectrons = iE;
	  electrons += gcr[ielectrons].cr_density;
	  
	  ostringstream buf;
	  buf << "CR " << gcr[ielectrons].name << " found as species #" << ielectrons;
	  INFO(buf.str());
	  
	}
      
      if (-1 == ielectrons) {
	
	ostringstream buf;
	buf << "CR electrons/positrons not found.";
	INFO(buf.str());
	electrons.delete_array();
	return 1;
	
      }

      valarray<double> gammaE(0., nEGammaBins), eps2(0., nEGammaBins);

      for (unsigned int i = 0; i < nEGammaBins; ++i)
	gammaE[i] = galaxy.E_gamma[i];

      eps2 = gammaE/Mele;

      Skymap<double> anisoSkymapOpt(internalHealpixOrder, gammaE);
      Skymap<double> anisoSkymapIR(internalHealpixOrder, gammaE);
      Skymap<double> anisoSkymapCMB(internalHealpixOrder, gammaE);
      Skymap<double> isoSkymapOpt(internalHealpixOrder, gammaE);
      Skymap<double> isoSkymapIR(internalHealpixOrder, gammaE);
      Skymap<double> isoSkymapCMB(internalHealpixOrder, gammaE);
      
      anisoSkymapOpt = 0;
      isoSkymapOpt = 0;
      anisoSkymapIR = 0;
      isoSkymapIR = 0;
      anisoSkymapCMB = 0;
      isoSkymapCMB = 0;

      const unsigned int electronBins = gcr[ielectrons].n_pgrid;

      valarray<double> electronE(0., electronBins), gamma(0., electronBins);

      for (unsigned int i = 0; i < electronBins; ++i)
	electronE[i] = gcr[ielectrons].Ekin[i];

      gamma = electronE/Mele;

      const unsigned int targetBins = galaxy.ISRF[0].n_pgrid;

      valarray<double> targetE(0., targetBins), eps1(0., targetBins), targetFreq(0., targetBins);

      for (unsigned int i = 0; i < targetBins; ++i) {

	targetFreq[i] = galaxy.nu_ISRF[i];
	targetE[i] = h_planck*galaxy.nu_ISRF[i]*erg_to_eV;  // target photon energy in eV

      }

      eps1 = targetE*1e-6/Mele;

      INFO("Precalculating cross section data");
	
      valarray<double> isoCrossSection(0., nEGammaBins*electronBins*targetBins), anisoCrossSection(0., nEGammaBins*electronBins*targetBins*cosThetaBins);

      for (unsigned long i = 0; i < nEGammaBins; ++i)
	for (unsigned long j = 0; j < electronBins; ++j)
	  for (unsigned long k = 0; k < targetBins; ++k) {
	    
	    unsigned long isoIndex = i*electronBins*targetBins + j*targetBins + k;
	    
	    isoCrossSection[isoIndex] = IsotropicICSSourceFunction(eps2[i], eps1[k], gamma[j]);

	    for (unsigned long l = 0; l < cosThetaBins; ++l) {
	      
	      unsigned long anisoIndex = i*electronBins*targetBins*cosThetaBins + j*targetBins*cosThetaBins + k*cosThetaBins + l;
	      
	      const double cosTheta = -1. + (l + 0.5)*dCosTheta;
	      
	      anisoCrossSection[anisoIndex] = AnisotropicICSSourceFunction(eps2[i], eps1[k], gamma[j], cosTheta);
	      
	    }

	  }   

      // Construct the radiation field

      const std::string fitsDirectory = configure.fFITSDataDirectory;
      const std::string isrfFilename = galdef.ISRF_file;
      const std::string filename = fitsDirectory + isrfFilename;
      
      ostringstream buf1;
      buf1 << "Reading ISRF from " << filename;
      INFO(buf1.str());
      
      rf::RadiationField rf(filename, targetFreq, galdef.ISRF_healpixOrder);

      const unsigned int rBins3D = std::max(galaxy.n_ygrid, galaxy.n_xgrid)/2; 
      
      const unsigned int rBins = (2 == gcr[0].n_spatial_dimensions ? galaxy.n_rgrid : utl::kSqrtTwo*rBins3D + 1);//galaxy.n_ygrid/2*galaxy.n_ygrid/2 + galaxy.n_xgrid/2*galaxy.n_xgrid/2) + 1);
      const unsigned int zBins = galaxy.n_zgrid;
      
      assert (nComps <= 3);

      valarray<double> rGrid(0., rBins);

      const double rMax3D = std::max(galaxy.y_max, galaxy.x_max);

      const double rMax = (2 == gcr[0].n_spatial_dimensions ? galaxy.r_max : utl::kSqrtTwo*rMax3D); 
      //	   sqrt(galaxy.x_max*galaxy.x_max + galaxy.y_max*galaxy.y_max));
      
      const int targetSkymapOrder = galdef.ISRF_healpixOrder;
      
      ostringstream buf2;
      buf2 << "Constructing ISRF angular distribution for healpix order " << galdef.ISRF_healpixOrder;
      INFO(buf2.str());
      
      Skymap<double> cmbSkymap(galdef.ISRF_healpixOrder, targetFreq);
      
      valarray<double> cmbNumberDensity(0., targetFreq.size());
      
      for (unsigned int i = 0; i < targetE.size(); ++i) {
	
	cmbNumberDensity[i] = BlackBodyNumberDensity(targetE[i], 2.735*utl::kBoltzmann_SI/utl::e_SI); // eV^-1 cm^-3 
		
      }
      
      for (unsigned int i = 0; i < cmbSkymap.Npix(); ++i)
	for (unsigned int iT = 0; iT < targetE.size(); ++iT)
	  cmbSkymap[i][iT] = (cmbNumberDensity[iT]*1./utl::kFourPi*cmbSkymap.solidAngle()); // Convert to eV^-1 cm^-3 * dSA/4pi per pixel

      // The intensity maps
      
      valarray< Skymap<double> > optAngDist(rBins*zBins), irAngDist(rBins*zBins);

      // The corresponding energy densities

      valarray< valarray<double> > optISRF(rBins*zBins), irISRF(rBins*zBins);

      // Build the array of skymaps and energy densities -- only use 2D ISRF
      // for now due to computational limitations (it takes much more time, and 
      // a lot of memory, to calculate a full 3D ISRF and I don't think 
      // much is gained in the galprop aniso IC calculation even if this is 
      // done)

      //#pragma omp parallel for schedule(dynamic) default(shared)
      for (unsigned int iR = 0; iR < rBins; ++iR)
	for (unsigned int iZ = 0; iZ < zBins; ++iZ) {
	  
	  //cout << "Built skymap " << iR << " " << iZ << endl;

	  const double r = rGrid[iR] = (2 == gcr[0].n_spatial_dimensions ? galaxy.r[iR] : rMax/rBins*iR);
	  const double z = galaxy.z[iZ];
	  
	  rf::RadiationField::ThreeVector pos(r, 0, z);
	  
	  const unsigned int index = iR*zBins + iZ;
	  
	  const Skymap<double> directSkymap = rf.GetSkymap(pos, rf::RadiationField::DIRECT, galdef.ISRF_healpixOrder)*cmbSkymap.solidAngle();
	  
	  const Skymap<double> scatteredSkymap = rf.GetSkymap(pos, rf::RadiationField::SCATTERED, galdef.ISRF_healpixOrder)*cmbSkymap.solidAngle();
	  
	  const Skymap<double> opticalSkymap = (directSkymap + scatteredSkymap);

	  const Skymap<double> transientSkymap = rf.GetSkymap(pos, rf::RadiationField::TRANSIENT, galdef.ISRF_healpixOrder)*cmbSkymap.solidAngle();
	  
	  const Skymap<double> thermalSkymap = rf.GetSkymap(pos, rf::RadiationField::THERMAL, galdef.ISRF_healpixOrder)*cmbSkymap.solidAngle();

	  const Skymap<double> infraredSkymap = (transientSkymap + thermalSkymap);

	  optAngDist[index] = opticalSkymap;

	  irAngDist[index] = infraredSkymap;
	  
	  //optAngDist[index].Resize(galdef.ISRF_healpixOrder, directSkymap.getSpectra());

	  //irAngDist[index].Resize(galdef.ISRF_healpixOrder, transientSkymap.getSpectra());

	  //for (unsigned int iPix = 0; iPix < directSkymap.Npix(); ++iPix) {

	  //for (unsigned int iFreq = 0; iFreq < directSkymap.getSpectra().size(); ++iFreq) {
	    
	  //  optAngDist[index][iPix][iFreq] = directSkymap[iPix][iFreq] + scatteredSkymap[iPix][iFreq];

	  //  irAngDist[index][iPix][iFreq] = transientSkymap[iPix][iFreq] + thermalSkymap[iPix][iFreq];

	  //  }
	  // }
	  	  
	  optISRF[index].resize(targetBins);
	  irISRF[index].resize(targetBins);
	  
	  for (unsigned int iTarget = 0; iTarget < targetBins; ++iTarget) {
	    
	    optISRF[index][iTarget] = optAngDist[index].sum(iTarget);
	    irISRF[index][iTarget] = irAngDist[index].sum(iTarget);
	    	    
	  }
	  
	  /*for (int iFreq = 0; iFreq < targetBins; ++iFreq) {
	    
	    const Skymap<double>& opt = optAngDist[index];
	    const Skymap<double>& ir = irAngDist[index];
	    
	    const double fac = targetE[iFreq]*targetE[iFreq];
	    
	    cout << iFreq << " " 
		 << targetE[iFreq] << " " 
		 << opt.sum(iFreq)*fac << " "
		 << optISRF[index][iFreq]*fac << " "
		 << galaxy.ISRF[0].d2[iR][iZ].s[iFreq] << " "
		 << ir.sum(iFreq)*fac << " "
		 << irISRF[index][iFreq]*fac << " "
		 << galaxy.ISRF[1].d2[iR][iZ].s[iFreq] << " "
		 << cmbSkymap.sum(iFreq)*fac << " "
		 << cmbNumberDensity[iFreq]*fac << " "
		 << galaxy.ISRF[2].d2[iR][iZ].s[iFreq] << endl;
	    
		 }*/
	  
	}
      
      //exit(0);

      // Precompute the inner product for target photon and emitted 
      // photon directions
      
      valarray<vec3> dirTarget(vec3(0, 0, 0), cmbSkymap.Npix());
      
      for (unsigned int i = 0; i < dirTarget.size(); ++i)
	dirTarget[i] = cmbSkymap.pix2coord(i).healpixAng().to_vec3();
      
      valarray<vec3> dirGamma(vec3(0, 0, 0), anisoSkymapOpt.Npix());
      
      for (unsigned int i = 0; i < dirGamma.size(); ++i)
	dirGamma[i] = anisoSkymapOpt.pix2coord(i).healpixAng().to_vec3();
      
      valarray<double> cosZeta(0., dirTarget.size()*dirGamma.size());
      
      for (unsigned int i = 0; i < dirTarget.size(); ++i)
	for (unsigned int j = 0; j < dirGamma.size(); ++j) {
	  
	  const unsigned int index = i*dirGamma.size() + j;
	  
	  cosZeta[index] = dotprod(-dirTarget[i], dirGamma[j]);
	  
	}
      
      // Easier to modify location of camera in future ...

      const double cameraX = Rsun;
      const double cameraY = 0;
      const double cameraZ = 0;

      const double ds = 0.25; // kpc -- hard wired so get predictable integration time. Errors due to the coarseness of this step size are effectively cancelled since we calculate the aniso/iso ratio

      INFO("Calculating interpolation skymaps");

#pragma omp parallel for schedule(dynamic) default(shared)
      for (int iPix = 0; iPix < anisoSkymapOpt.Npix(); ++iPix) {

	//ostringstream buf;
	//buf << "Beginning pixel " << iPix << " of " << anisoSkymapOpt.Npix();
	//INFO(buf.str());
	
	SM::Coordinate coord(anisoSkymapOpt.pix2ang(iPix));
	
	const double l = coord.l(), b = coord.b();
	const double lRad = l*utl::kConvertDegreesToRadians, bRad = b*utl::kConvertDegreesToRadians;

	const double sinb = sin(bRad), cosb = cos(bRad), sinl = sin(lRad), cosl = cos(lRad);
	  
	double s = 0;

	bool complete = false;

	while (1) {

	  s += ds;
	  
	  const double dx = s*cosb*cosl;
	  const double dy = s*cosb*sinl;
	  const double dz = s*sinb;
	  
	  const double x = cameraX - dx;
	  const double y = cameraY - dy;
	  const double z = cameraZ + dz;
	  
	  const double r = sqrt(x*x + y*y);

	  // These are the loop break criteria. Much simpler encoding of 
	  // the break condition than in the earlier pixel routines.

	  if (2 == gcr[0].n_spatial_dimensions) {

	    complete = (r >= galaxy.r_max || 
			z <= galaxy.z_min || z >= galaxy.z_max);

	  }

	  if (3 == gcr[0].n_spatial_dimensions) {

	    complete = (x <= galaxy.x_min || x >= galaxy.x_max || 
			y <= galaxy.y_min || y >= galaxy.y_max || 
			z <= galaxy.z_min || z >= galaxy.z_max);

	  }

	  //cout << iPix << " " << x << " " << y << " " << z << " " << r << " " << complete << endl;

	  if (complete)
	    break;

	  // Determine angle of gamma rays toward observer from the 
	  // emission region. This is the direction we use for the inner product
	  // in the cross section calculation. Having found the direction we 
	  // just index into the precomputed cosZeta array. Recall also that
	  // we are doing it this way because the coordinate system for the 
	  // ISRF is centred toward the GC due to the azimuthal symmetry.

	  const double cosXi = z/s;
	  const double cos2Xi = cosXi*cosXi;
	  const double sin2Xi = (cos2Xi > 1. ? 0. : 1. - cos2Xi);
	  const double sinXi = sqrt(sin2Xi);
	  
	  const double theta = atan2(y, x);
	  
	  const double rho = sqrt(dx*dx + dy*dy);

	  const double sinChi = -cameraX/rho*sin(theta);
	  const double cosChi = -(r - cameraX*cos(theta))/rho;

	  const vec3 dir(cosChi*sinXi, sinChi*sinXi, cosXi);

	  pointing point(dir);
	  
	  // This is the index into the inner product array calculated earlier

	  const int gammaPixIndex = anisoSkymapOpt.ang2pix(point);

	  valarray<double> anisoEmissivityOpt(0., nEGammaBins), isoEmissivityOpt(0., nEGammaBins);
	  valarray<double> anisoEmissivityIR(0., nEGammaBins), isoEmissivityIR(0., nEGammaBins);
	  valarray<double> anisoEmissivityCMB(0., nEGammaBins), isoEmissivityCMB(0., nEGammaBins);
	  
	  // Essentially the same algorithm for 2D and 3D: if we're at the 
	  // spatial bin boundary, use the electron spectra, target intensities
	  // and energy densities, there. Otherwise, interpolate (bilinear for
	  // 2D, trilinear for 3D) amongst the 4(8) surrounding grid points
	  // to get the electron spectra, etc. Once that is done, it's the 
	  // usual calculation for the isotropic cross section. For the 
	  // anisotropic cross section, there is an inner loop over the 
	  // number of pixels in the target skymap that sums the anisotropic
	  // cross section weighted by the target photon intensity at each
	  // pixel. The anisotropic cross section is evaluated for the 
	  // cosTheta value corresponding to the emission angle toward 
	  // the camera.

	  if (2 == gcr[0].n_spatial_dimensions) {

	    const double dR = (r - galaxy.r_min)/(galaxy.r_max - galaxy.r_min);
	    const double dZ = (z - galaxy.z_min)/(galaxy.z_max - galaxy.z_min);

	    const unsigned int iR = (dR < 1. ? dR*(galaxy.n_rgrid-1) : galaxy.n_rgrid-1);
	    const unsigned int iZ = (dZ < 1. ? dZ*(galaxy.n_zgrid-1) : galaxy.n_zgrid-1);
      
	    if (dR >=1 || dZ >= 1) {

	      const int index = iR*galaxy.n_zgrid + iZ;

	      const Skymap<double>& opt = optAngDist[index];
	      const Skymap<double>& ir = irAngDist[index];

	      for (unsigned int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {

		for (unsigned int iElectron = 0; iElectron < electronBins; ++iElectron) {

		  double anisoSumOpt = 0, isoSumOpt = 0;
		  double anisoSumIR = 0, isoSumIR = 0;
		  double anisoSumCMB = 0, isoSumCMB = 0; // For checking

		  for (unsigned int iTarget = 0; iTarget < targetBins; ++iTarget) {

		    const unsigned int isoIndex = iEGamma*electronBins*targetBins + iElectron*targetBins + iTarget;

		    // This is the isotropic cross section summation
		    
		    isoSumOpt += isoCrossSection[isoIndex]*
		      targetE[iTarget]*log(targetE[1]/targetE[0])*
		      optISRF[index][iTarget];

		    isoSumIR += isoCrossSection[isoIndex]*
		      targetE[iTarget]*log(targetE[1]/targetE[0])*
		      irISRF[index][iTarget];

		    isoSumCMB += isoCrossSection[isoIndex]*
		      targetE[iTarget]*log(targetE[1]/targetE[0])*
		      cmbNumberDensity[iTarget];

		    // Inner loop over the target skymaps
		
		    double sumSkymapOpt = 0, sumSkymapIR = 0, sumSkymapCMB = 0;
    
		    for (unsigned int iSkymap = 0; iSkymap < opt.Npix(); ++iSkymap) {

		      const unsigned int cosZetaIndex = iSkymap*dirGamma.size() + gammaPixIndex;

		      const unsigned int cosThetaIndex = 
			((1. + cosZeta[cosZetaIndex]) < 2. ? int(0.5*(1. + cosZeta[cosZetaIndex])*cosThetaBins) : cosThetaBins-1); 

		      const unsigned int anisoIndex = iEGamma*electronBins*targetBins*cosThetaBins + iElectron*targetBins*cosThetaBins + iTarget*cosThetaBins + cosThetaIndex;
		      
		      sumSkymapOpt += anisoCrossSection[anisoIndex]*
			opt[iSkymap][iTarget];
		      
		      sumSkymapIR += anisoCrossSection[anisoIndex]*
			ir[iSkymap][iTarget];
		      
		      sumSkymapCMB += anisoCrossSection[anisoIndex]*
			cmbSkymap[iSkymap][iTarget];

		    }
		     
		    anisoSumOpt += sumSkymapOpt*targetE[iTarget]*log(targetE[1]/targetE[0]);
		    anisoSumIR += sumSkymapIR*targetE[iTarget]*log(targetE[1]/targetE[0]);
		    anisoSumCMB += sumSkymapCMB*targetE[iTarget]*log(targetE[1]/targetE[0]);

		  }

		  const double elecSum = electrons.d2[iR][iZ].s[iElectron]*electronE[iElectron]*log(electronE[1]/electronE[0]);
		  
		  isoEmissivityOpt[iEGamma] += isoSumOpt*elecSum;
		  isoEmissivityIR[iEGamma] += isoSumIR*elecSum;
		  isoEmissivityCMB[iEGamma] += isoSumCMB*elecSum;

		  anisoEmissivityOpt[iEGamma] += anisoSumOpt*elecSum;
		  anisoEmissivityIR[iEGamma] += anisoSumIR*elecSum;
		  anisoEmissivityCMB[iEGamma] += anisoSumCMB*elecSum;

		}

	      }

	    } else {

	      // Bilinear interpolation

	      const double rCoeff = (r - galaxy.r[iR])/(galaxy.r[iR+1] - galaxy.r[iR]);
	      const double zCoeff = (z - galaxy.z[iZ])/(galaxy.z[iZ+1] - galaxy.z[iZ]);

	      const unsigned int index1 = iR*galaxy.n_zgrid + iZ;
	      const unsigned int index2 = (iR+1)*galaxy.n_zgrid + iZ;
	      const unsigned int index3 = iR*galaxy.n_zgrid + (iZ+1);
	      const unsigned int index4 = (iR+1)*galaxy.n_zgrid + (iZ+1);

	      const Skymap<double>& opt1 = optAngDist[index1];
	      const Skymap<double>& ir1 = irAngDist[index1];

	      const Skymap<double>& opt2 = optAngDist[index2];
	      const Skymap<double>& ir2 = irAngDist[index2];

	      const Skymap<double>& opt3 = optAngDist[index3];
	      const Skymap<double>& ir3 = irAngDist[index3];

	      const Skymap<double>& opt4 = optAngDist[index4];
	      const Skymap<double>& ir4 = irAngDist[index4];

	      for (unsigned int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {

		for (unsigned int iElectron = 0; iElectron < electronBins; ++iElectron) {
		  
		  double anisoSumOpt = 0, isoSumOpt = 0;
		  double anisoSumIR = 0, isoSumIR = 0;
		  double anisoSumCMB = 0, isoSumCMB = 0; // For checking

		  for (unsigned int iTarget = 0; iTarget < targetBins; ++iTarget) {

		    const unsigned int isoIndex = iEGamma*electronBins*targetBins + iElectron*targetBins + iTarget;

		    // This is the isotropic cross section summation
		    
		    isoSumOpt += 
		      isoCrossSection[isoIndex]*
		      targetE[iTarget]*log(targetE[1]/targetE[0])*
		      (optISRF[index1][iTarget]*(1. - rCoeff)*(1. - zCoeff) +
		       optISRF[index2][iTarget]*rCoeff*(1. - zCoeff) +
		       optISRF[index3][iTarget]*(1. - rCoeff)*zCoeff +
		       optISRF[index4][iTarget]*rCoeff*zCoeff);
		   
		    isoSumIR += 
		      isoCrossSection[isoIndex]*
		      targetE[iTarget]*log(targetE[1]/targetE[0])*
		      (irISRF[index1][iTarget]*(1. - rCoeff)*(1. - zCoeff) +
		       irISRF[index2][iTarget]*rCoeff*(1. - zCoeff) +
		       irISRF[index3][iTarget]*(1. - rCoeff)*zCoeff +
		       irISRF[index4][iTarget]*rCoeff*zCoeff);
						
		    isoSumCMB += 
		      isoCrossSection[isoIndex]*
		      targetE[iTarget]*log(targetE[1]/targetE[0])*
		      cmbNumberDensity[iTarget];
		
		    // Inner loop over the target skymaps

		    double sumSkymapOpt = 0, sumSkymapIR = 0, sumSkymapCMB = 0;
    
		    for (unsigned int iSkymap = 0; iSkymap < opt1.Npix(); ++iSkymap) {

		      const unsigned int cosZetaIndex = iSkymap*dirGamma.size() + gammaPixIndex;

		      const unsigned int cosThetaIndex = 
			((1. + cosZeta[cosZetaIndex]) < 2. ? int(0.5*(1. + cosZeta[cosZetaIndex])*cosThetaBins) : cosThetaBins-1); 

		      const unsigned int anisoIndex = iEGamma*electronBins*targetBins*cosThetaBins + iElectron*targetBins*cosThetaBins + iTarget*cosThetaBins + cosThetaIndex;
		    
		      sumSkymapOpt += anisoCrossSection[anisoIndex]*
			(opt1[iSkymap][iTarget]*(1. - rCoeff)*(1. - zCoeff) +
			 opt2[iSkymap][iTarget]*rCoeff*(1. - zCoeff) +
			 opt3[iSkymap][iTarget]*(1. - rCoeff)*zCoeff +
			 opt4[iSkymap][iTarget]*rCoeff*zCoeff);

		      sumSkymapIR += anisoCrossSection[anisoIndex]*
			(ir1[iSkymap][iTarget]*(1. - rCoeff)*(1. - zCoeff) +
			 ir2[iSkymap][iTarget]*rCoeff*(1. - zCoeff) +
			 ir3[iSkymap][iTarget]*(1. - rCoeff)*zCoeff +
			 ir4[iSkymap][iTarget]*rCoeff*zCoeff);
		      
		      sumSkymapCMB += anisoCrossSection[anisoIndex]*
			cmbSkymap[iSkymap][iTarget];

		    }
		     
		    anisoSumOpt += sumSkymapOpt*targetE[iTarget]*log(targetE[1]/targetE[0]);
		    anisoSumIR += sumSkymapIR*targetE[iTarget]*log(targetE[1]/targetE[0]);
		    anisoSumCMB += sumSkymapCMB*targetE[iTarget]*log(targetE[1]/targetE[0]);

		  }

		  // All 4 grid points for the bilinear interpolation
		  
		  const double elecSum =
		    (electrons.d2[iR][iZ].s[iElectron]*(1. - rCoeff)*(1. - zCoeff) +
		     electrons.d2[iR+1][iZ].s[iElectron]*rCoeff*(1. - zCoeff) +
		     electrons.d2[iR][iZ+1].s[iElectron]*(1. - rCoeff)*zCoeff +
		     electrons.d2[iR+1][iZ+1].s[iElectron]*rCoeff*zCoeff)*
		    electronE[iElectron]*log(electronE[1]/electronE[0]);
		  
		  isoEmissivityOpt[iEGamma] += isoSumOpt*elecSum;
		  isoEmissivityIR[iEGamma] += isoSumIR*elecSum;
		  isoEmissivityCMB[iEGamma] += isoSumCMB*elecSum;

		  anisoEmissivityOpt[iEGamma] += anisoSumOpt*elecSum;
		  anisoEmissivityIR[iEGamma] += anisoSumIR*elecSum;
		  anisoEmissivityCMB[iEGamma] += anisoSumCMB*elecSum;

		  //cout << iEGamma << " " << iElectron << " " << isoSumCMB << endl;

		}

	      }

	    }

	  }

	  if (3 == gcr[0].n_spatial_dimensions) {

	    const double dX = (x - galaxy.x_min)/(galaxy.x_max - galaxy.x_min);
	    const double dY = (y - galaxy.y_min)/(galaxy.y_max - galaxy.y_min);
	    const double dZ = (z - galaxy.z_min)/(galaxy.z_max - galaxy.z_min);

	    const unsigned int iX = (dX < 1. ? dX*(galaxy.n_xgrid-1) : galaxy.n_xgrid-1);
	    const unsigned int iY = (dY < 1. ? dY*(galaxy.n_ygrid-1) : galaxy.n_ygrid-1);
	    const unsigned int iZ = (dZ < 1. ? dZ*(galaxy.n_zgrid-1) : galaxy.n_zgrid-1);

	    // Radiation field is sampled only on a R,z grid due to 
	    // computational limitations (for now).

	    const double dR = r/rMax;
	    
	    const unsigned int iR = (dR < 1. ? dR*(rBins-1) : rBins-1);	    
   
	    if (dX >= 1 || dY >= 1 || dZ >= 1 || iR >= 1) {

	      const int index = iR*galaxy.n_zgrid + iZ;

	      const Skymap<double>& opt = optAngDist[index];
	      const Skymap<double>& ir = irAngDist[index];

	      for (unsigned int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {

		for (unsigned int iElectron = 0; iElectron < electronBins; ++iElectron) {

		  double anisoSumOpt = 0, isoSumOpt = 0;
		  double anisoSumIR = 0, isoSumIR = 0;
		  double anisoSumCMB = 0, isoSumCMB = 0; // For checking

		  for (unsigned int iTarget = 0; iTarget < targetBins; ++iTarget) {

		    const unsigned int isoIndex = iEGamma*electronBins*targetBins + iElectron*targetBins + iTarget;

		    // This is the isotropic cross section summation
		    
		    isoSumOpt += isoCrossSection[isoIndex]*
		      targetE[iTarget]*log(targetE[1]/targetE[0])*
		      optISRF[index][iTarget];

		    isoSumIR += isoCrossSection[isoIndex]*
		      targetE[iTarget]*log(targetE[1]/targetE[0])*
		      irISRF[index][iTarget];

		    isoSumCMB += isoCrossSection[isoIndex]*
		      targetE[iTarget]*log(targetE[1]/targetE[0])*
		      cmbNumberDensity[iTarget];

		    // Inner loop over target skymaps
		
		    double sumSkymapOpt = 0, sumSkymapIR = 0, sumSkymapCMB = 0;
    
		    for (unsigned int iSkymap = 0; iSkymap < opt.Npix(); ++iSkymap) {

		      const unsigned int cosZetaIndex = iSkymap*dirGamma.size() + gammaPixIndex;

		      const unsigned int cosThetaIndex = 
			((1. + cosZeta[cosZetaIndex]) < 2. ? int(0.5*(1. + cosZeta[cosZetaIndex])*cosThetaBins) : cosThetaBins-1); 

		      const unsigned int anisoIndex = iEGamma*electronBins*targetBins*cosThetaBins + iElectron*targetBins*cosThetaBins + iTarget*cosThetaBins + cosThetaIndex;
		      
		      sumSkymapOpt += anisoCrossSection[anisoIndex]*
			opt[iSkymap][iTarget];
		      
		      sumSkymapIR += anisoCrossSection[anisoIndex]*
			ir[iSkymap][iTarget];
		      
		      sumSkymapCMB += anisoCrossSection[anisoIndex]*
			cmbSkymap[iSkymap][iTarget];

		    }
		     
		    anisoSumOpt += sumSkymapOpt*targetE[iTarget]*log(targetE[1]/targetE[0]);
		    anisoSumIR += sumSkymapIR*targetE[iTarget]*log(targetE[1]/targetE[0]);
		    anisoSumCMB += sumSkymapCMB*targetE[iTarget]*log(targetE[1]/targetE[0]);

		  }

		  const double elecSum = electrons.d3[iX][iY][iZ].s[iElectron]*electronE[iElectron]*log(electronE[1]/electronE[0]);
		  
		  isoEmissivityOpt[iEGamma] += isoSumOpt*elecSum;
		  isoEmissivityIR[iEGamma] += isoSumIR*elecSum;
		  isoEmissivityCMB[iEGamma] += isoSumCMB*elecSum;

		  anisoEmissivityOpt[iEGamma] += anisoSumOpt*elecSum;
		  anisoEmissivityIR[iEGamma] += anisoSumIR*elecSum;
		  anisoEmissivityCMB[iEGamma] += anisoSumCMB*elecSum;

		}

	      }

	    } else {

	      // Trilinear interpolation

	      const double xCoeff = (x - galaxy.x[iX])/(galaxy.x[iX+1] - galaxy.x[iX]);

	      const double yCoeff = (y - galaxy.y[iY])/(galaxy.y[iY+1] - galaxy.y[iY]);

	      const double zCoeff = (z - galaxy.z[iZ])/(galaxy.z[iZ+1] - galaxy.z[iZ]);

	      const double rCoeff = (r - rGrid[iR])/(rGrid[iR+1] - rGrid[iR]);
 
	      const unsigned int index1 = iR*galaxy.n_zgrid + iZ;
	      const unsigned int index2 = (iR+1)*galaxy.n_zgrid + iZ;
	      const unsigned int index3 = iR*galaxy.n_zgrid + (iZ+1);
	      const unsigned int index4 = (iR+1)*galaxy.n_zgrid + (iZ+1);

	      const Skymap<double>& opt1 = optAngDist[index1];
	      const Skymap<double>& ir1 = irAngDist[index1];

	      const Skymap<double>& opt2 = optAngDist[index2];
	      const Skymap<double>& ir2 = irAngDist[index2];

	      const Skymap<double>& opt3 = optAngDist[index3];
	      const Skymap<double>& ir3 = irAngDist[index3];

	      const Skymap<double>& opt4 = optAngDist[index4];
	      const Skymap<double>& ir4 = irAngDist[index4];

	      for (unsigned int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {

		for (unsigned int iElectron = 0; iElectron < electronBins; ++iElectron) {
		  
		  double anisoSumOpt = 0, isoSumOpt = 0;
		  double anisoSumIR = 0, isoSumIR = 0;
		  double anisoSumCMB = 0, isoSumCMB = 0; // For checking

		  for (unsigned int iTarget = 0; iTarget < targetBins; ++iTarget) {

		    const unsigned int isoIndex = iEGamma*electronBins*targetBins + iElectron*targetBins + iTarget;

		    // This is the isotropic cross section summation
		    
		    isoSumOpt += 
		      isoCrossSection[isoIndex]*
		      targetE[iTarget]*log(targetE[1]/targetE[0])*
		      (optISRF[index1][iTarget]*(1. - rCoeff)*(1. - zCoeff) +
		       optISRF[index2][iTarget]*rCoeff*(1. - zCoeff) +
		       optISRF[index3][iTarget]*(1. - rCoeff)*zCoeff +
		       optISRF[index4][iTarget]*rCoeff*zCoeff);
		   
		    isoSumIR += 
		      isoCrossSection[isoIndex]*
		      targetE[iTarget]*log(targetE[1]/targetE[0])*
		      (irISRF[index1][iTarget]*(1. - rCoeff)*(1. - zCoeff) +
		       irISRF[index2][iTarget]*rCoeff*(1. - zCoeff) +
		       irISRF[index3][iTarget]*(1. - rCoeff)*zCoeff +
		       irISRF[index4][iTarget]*rCoeff*zCoeff);
						
		    isoSumCMB += 
		      isoCrossSection[isoIndex]*
		      targetE[iTarget]*log(targetE[1]/targetE[0])*
		      cmbNumberDensity[iTarget];

		    // Inner loop over target skymaps
		
		    double sumSkymapOpt = 0, sumSkymapIR = 0, sumSkymapCMB = 0;
    
		    for (unsigned int iSkymap = 0; iSkymap < opt1.Npix(); ++iSkymap) {

		      const unsigned int cosZetaIndex = iSkymap*dirGamma.size() + gammaPixIndex;

		      const unsigned int cosThetaIndex = 
			((1. + cosZeta[cosZetaIndex]) < 2. ? int(0.5*(1. + cosZeta[cosZetaIndex])*cosThetaBins) : cosThetaBins-1); 

		      const unsigned int anisoIndex = iEGamma*electronBins*targetBins*cosThetaBins + iElectron*targetBins*cosThetaBins + iTarget*cosThetaBins + cosThetaIndex;
		    
		      sumSkymapOpt += anisoCrossSection[anisoIndex]*
			(opt1[iSkymap][iTarget]*(1. - rCoeff)*(1. - zCoeff) +
			 opt2[iSkymap][iTarget]*rCoeff*(1. - zCoeff) +
			 opt3[iSkymap][iTarget]*(1. - rCoeff)*zCoeff +
			 opt4[iSkymap][iTarget]*rCoeff*zCoeff);

		      sumSkymapIR += anisoCrossSection[anisoIndex]*
			(ir1[iSkymap][iTarget]*(1. - rCoeff)*(1. - zCoeff) +
			 ir2[iSkymap][iTarget]*rCoeff*(1. - zCoeff) +
			 ir3[iSkymap][iTarget]*(1. - rCoeff)*zCoeff +
			 ir4[iSkymap][iTarget]*rCoeff*zCoeff);
		      
		      sumSkymapCMB += anisoCrossSection[anisoIndex]*
			cmbSkymap[iSkymap][iTarget];

		    }
		     
		    anisoSumOpt += sumSkymapOpt*targetE[iTarget]*log(targetE[1]/targetE[0]);
		    anisoSumIR += sumSkymapIR*targetE[iTarget]*log(targetE[1]/targetE[0]);
		    anisoSumCMB += sumSkymapCMB*targetE[iTarget]*log(targetE[1]/targetE[0]);

		  }
		  
		  // All 8 grid points for the trilinear interpolation

		  const double elecSum =
		    (electrons.d3[iX][iY][iZ].s[iElectron]*(1. - xCoeff)*(1. - yCoeff)*(1. - zCoeff) +
		     electrons.d3[iX+1][iY][iZ].s[iElectron]*xCoeff*(1. - yCoeff)*(1. - zCoeff) +
		     electrons.d3[iX][iY+1][iZ].s[iElectron]*(1. - xCoeff)*yCoeff*(1. - zCoeff) +
		     electrons.d3[iX+1][iY+1][iZ].s[iElectron]*xCoeff*yCoeff*(1. - zCoeff) +
		     electrons.d3[iX][iY][iZ+1].s[iElectron]*(1. - xCoeff)*(1. - yCoeff)*zCoeff +
		     electrons.d3[iX+1][iY][iZ+1].s[iElectron]*xCoeff*(1. - yCoeff)*zCoeff +
		     electrons.d3[iX][iY+1][iZ+1].s[iElectron]*(1. - xCoeff)*yCoeff*zCoeff +
		     electrons.d3[iX+1][iY+1][iZ+1].s[iElectron]*xCoeff*yCoeff*zCoeff)*
		    electronE[iElectron]*log(electronE[1]/electronE[0]);
		  
		  isoEmissivityOpt[iEGamma] += isoSumOpt*elecSum;
		  isoEmissivityIR[iEGamma] += isoSumIR*elecSum;
		  isoEmissivityCMB[iEGamma] += isoSumCMB*elecSum;

		  anisoEmissivityOpt[iEGamma] += anisoSumOpt*elecSum;
		  anisoEmissivityIR[iEGamma] += anisoSumIR*elecSum;
		  anisoEmissivityCMB[iEGamma] += anisoSumCMB*elecSum;

		  //cout << iEGamma << " " << iElectron << " " << isoSumCMB << endl;

		}

	      }

	    }

	  }

	  for (int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {

	    anisoSkymapOpt[iPix][iEGamma] += anisoEmissivityOpt[iEGamma]*ds*kpc2cm;
	    anisoSkymapIR[iPix][iEGamma] += anisoEmissivityIR[iEGamma]*ds*kpc2cm;
	    anisoSkymapCMB[iPix][iEGamma] += anisoEmissivityCMB[iEGamma]*ds*kpc2cm;
	    
	    isoSkymapOpt[iPix][iEGamma] += isoEmissivityOpt[iEGamma]*ds*kpc2cm;
	    isoSkymapIR[iPix][iEGamma] += isoEmissivityIR[iEGamma]*ds*kpc2cm;
	    isoSkymapCMB[iPix][iEGamma] += isoEmissivityCMB[iEGamma]*ds*kpc2cm;
	    
	  }

	}

#pragma omp critical 
	{
	  ostringstream buf;
	  buf << "Finished pixel " << iPix << " of " << anisoSkymapOpt.Npix();
	  INFO(buf.str());
	}

	/*for (unsigned int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {

	  SM::Coordinate coord(anisoSkymapOpt.pix2ang(iPix));
	
	  const double l = coord.l(), b = coord.b();
	  
	  cout << "EGamma: " 
	       << iPix << " " 
	       << l << " " << b << " " 
	       << iEGamma << " " 
	       << isoSkymapCMB[iPix][iEGamma] << " "
	    //<< anisoSkymapCMB[iPix][iEGamma] << " "
	       << (isoSkymapCMB[iPix][iEGamma] > 0. ? anisoSkymapCMB[iPix][iEGamma]/isoSkymapCMB[iPix][iEGamma] : 0.) << " "
	       << isoSkymapIR[iPix][iEGamma] << " "
	    // << anisoSkymapIR[iPix][iEGamma] << " "
	       << (isoSkymapIR[iPix][iEGamma] > 0. ? anisoSkymapIR[iPix][iEGamma]/isoSkymapIR[iPix][iEGamma] : 0.) << " "
	       << isoSkymapOpt[iPix][iEGamma] << " "
	    // << anisoSkymapOpt[iPix][iEGamma] << " "
	       << (isoSkymapOpt[iPix][iEGamma] > 0. ? anisoSkymapOpt[iPix][iEGamma]/isoSkymapOpt[iPix][iEGamma] : 0.) << endl;

	       }*/

      }

      const Skymap<double> optRatio = anisoSkymapOpt/isoSkymapOpt;
      const Skymap<double> irRatio = anisoSkymapIR/isoSkymapIR;
      const Skymap<double> cmbRatio = anisoSkymapCMB/isoSkymapCMB;
 
      // The skymaps calculated above are most likely done for a healpix 
      // pixelisation that is considerably coarser than the pixelisation 
      // used for the skymaps to be output to disk. Here we interpolate 
      // the skymaps calculated above to the higher order pixelisation 
      // and obtain a ratio map as a function of energy. This can be applied 
      // to the isotropic IC skymap calculated below to obtain the full 
      // anisotropic skymap for the higher pixelisations (order 6, 7, ..) 
      // typically used in calculations.

      const unsigned int healpixOrder = galdef.healpix_order;
 
      anisoRatioOpt = optRatio.interpolate(healpixOrder);
      anisoRatioIR = irRatio.interpolate(healpixOrder);
      anisoRatioCMB = cmbRatio.interpolate(healpixOrder);

      const string fN = configure.fOutputDirectory + configure.fOutputPrefix;
    
      //isoSkymapCMB.write(fN + "iso_cmb_hp.fits.gz");
      //anisoSkymapCMB.write(fN + "aniso_cmb_hp.fits.gz");
      //cmbRatio.write(fN + "ratio_cmb_hp.fits.gz");

      //isoSkymapIR.write(fN + "iso_ir_hp.fits.gz");
      //anisoSkymapIR.write(fN + "aniso_ir_hp.fits.gz");
      //irRatio.write(fN + "ratio_ir_hp.fits.gz");

      //isoSkymapOpt.write(fN + "iso_opt_hp.fits.gz");
      //anisoSkymapOpt.write(fN + "aniso_opt_hp.fits.gz");
      //optRatio.write(fN + "ratio_opt_hp.fits.gz");

      //anisoRatioOpt.write(fN + "aniso_ratio_opt_inter_hp.fits.gz");
      //anisoRatioIR.write(fN + "aniso_ratio_ir_inter_hp.fits.gz");
      //anisoRatioCMB.write(fN + "aniso_ratio_cmb_inter_hp.fits.gz");

      electrons.delete_array();

    }
    
  } 

  if (3 == galdef.skymap_format || 4 == galdef.skymap_format) {

    //cout << "iPix: ";

#pragma omp parallel for schedule(dynamic) default(shared) 
    for (int iPix = 0; iPix < galaxy.IC_iso_hp_skymap[0].Npix(); ++iPix) {

      //cout << iPix << " ";

      SM::Coordinate co(galaxy.IC_iso_hp_skymap[0].pix2ang(iPix));
      const double l = co.l();
      const double b = co.b();
      vector< vector<double> > iso_IC;//, aniso_IC;
      
      //gen_IC_skymap_pixel(l, b, electrons, ielectrons, iso_IC, aniso_IC, Etarget, factor);
      
      gen_IC_skymap_pixel(l, b, iso_IC);

      for (int iComp = 0; iComp < nComps; ++iComp)
	for (int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {

	  galaxy.IC_iso_hp_skymap[iComp][co][iEGamma] = iso_IC[iComp][iEGamma]*galdef.ISRF_factors[iComp];

	  //cout << iComp << " " << iEGamma << " " << galaxy.IC_iso_hp_skymap[iComp][co][iEGamma];

	  //if(galdef.IC_anisotropic) galaxy.IC_aniso_hp_skymap[i_comp][co][iEgamma] = aniso_IC[i_comp][iEgamma]*galdef.ISRF_factors[i_comp];
	
	  //cout << endl;

	}
    }

    //cout << endl;

    if (galdef.IC_anisotropic && anisoValid) {

      for (int iPix = 0; iPix < galaxy.IC_iso_hp_skymap[0].Npix(); ++iPix) {

	for (int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {

	  galaxy.IC_aniso_hp_skymap[0][iPix][iEGamma] = double(anisoRatioOpt[iPix][iEGamma])*galaxy.IC_iso_hp_skymap[0][iPix][iEGamma];
	  
	  galaxy.IC_aniso_hp_skymap[1][iPix][iEGamma] = double(anisoRatioIR[iPix][iEGamma])*galaxy.IC_iso_hp_skymap[1][iPix][iEGamma];
	  
	  galaxy.IC_aniso_hp_skymap[2][iPix][iEGamma] = double(anisoRatioCMB[iPix][iEGamma])*galaxy.IC_iso_hp_skymap[2][iPix][iEGamma];

	}

      }
    
    }
	
    for (int iComp = 0; iComp < nComps; ++iComp) {

      galaxy.IC_iso_hp_skymap[iComp].setSpectra(galaxy.E_gamma, nEGammaBins);
      
      if (galdef.IC_anisotropic && anisoValid) 
	galaxy.IC_aniso_hp_skymap[iComp].setSpectra(galaxy.E_gamma, nEGammaBins);

    }

    if (galdef.verbose >= 2) {

      for (unsigned int iComp = 0; iComp < nComps; ++iComp) {

	ostringstream isoBuf;
	isoBuf << "Isotropic inverse Compton skymap for ISRF component " << iComp;
	INFO(isoBuf.str());

	galaxy.IC_iso_hp_skymap[iComp].print(cout);
	
	if (galdef.IC_anisotropic && anisoValid) {

	  ostringstream anisoBuf;
	  anisoBuf << "Anisotropic inverse Compton skymap for component " << iComp;
	  INFO(anisoBuf.str());

	  galaxy.IC_aniso_hp_skymap[iComp].print(cout);
      
	}
      
      }

    }
      
  } else {

    cout << "iLong: ";

#pragma omp parallel for schedule(dynamic) default(shared)
    for (int iLong = 0; iLong < galaxy.n_long; ++iLong) {

      cout << iLong << " ";
      
      for (int iLat = 0; iLat < galaxy.n_lat; ++iLat) {

	const double l = galaxy.long_min + iLong*galaxy.d_long;
	const double b = galaxy.lat_min + iLat*galaxy.d_lat;
	vector< vector<double> > iso_IC;//, aniso_IC;
	
	//gen_IC_skymap_pixel(l, b, electrons, ielectrons, iso_IC, aniso_IC, Etarget, factor);
	
	gen_IC_skymap_pixel(l, b, iso_IC);

	//Store and apply user defined factors to compontents

	for (int iComp = 0; iComp < nComps; ++iComp)
	  for (int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {
	    
	    galaxy.IC_iso_skymap[iComp].d2[iLong][iLat].s[iEGamma] = iso_IC[iComp][iEGamma]*galdef.ISRF_factors[iComp];

	    //if(galdef.IC_anisotropic) galaxy.IC_aniso_skymap[i_comp].d2[i_long][i_lat].s[iEgamma] = aniso_IC[i_comp][iEgamma]*galdef.ISRF_factors[i_comp];
	  }
    
      }//lat
    
    }//long

    cout << endl;

    if (galdef.IC_anisotropic && anisoValid) {

      for (int iLong = 0; iLong < galaxy.n_long; ++iLong) {
	
	for (int iLat = 0; iLat < galaxy.n_lat; ++iLat) {
	  
	  const double l = galaxy.long_min + iLong*galaxy.d_long;
	  const double b = galaxy.lat_min + iLat*galaxy.d_lat;
	  
	  SM::Coordinate coord(l, b);

	  const int iPix = anisoRatioOpt.ang2pix(coord.healpixAng());
	
	  for (int iEGamma = 0; iEGamma < nEGammaBins; ++iEGamma) {
	  
	    galaxy.IC_aniso_skymap[0].d2[iLong][iLat].s[iEGamma] = galaxy.IC_iso_skymap[0].d2[iLong][iLat].s[iEGamma]*double(anisoRatioOpt[iPix][iEGamma]);

	    galaxy.IC_aniso_skymap[1].d2[iLong][iLat].s[iEGamma] = galaxy.IC_iso_skymap[1].d2[iLong][iLat].s[iEGamma]*double(anisoRatioIR[iPix][iEGamma]);

	    galaxy.IC_aniso_skymap[2].d2[iLong][iLat].s[iEGamma] = galaxy.IC_iso_skymap[2].d2[iLong][iLat].s[iEGamma]*double(anisoRatioCMB[iPix][iEGamma]);

	  }

	}
    
      }

    }

    if (galdef.verbose >= 2) {

      for (int iComp = 0; iComp < nComps; ++iComp) {

	ostringstream isoBuf;
	isoBuf << "Isotropic inverse Compton skymap for ISRF component " << iComp;
	INFO(isoBuf.str());

	galaxy.IC_iso_skymap[iComp].print();
	
	if (galdef.IC_anisotropic && anisoValid) {

	  ostringstream anisoBuf;
	  anisoBuf << "Anisotropic inverse Compton skymap for component " << iComp;
	  INFO(anisoBuf.str());

	  galaxy.IC_aniso_skymap[iComp].print();
      
	}
      
      }
    
    } // galdef.verbose>=2
  
  }

  //if(galdef.IC_anisotropic==1) {

  //electrons.delete_array();  // IMOS20060420
  //}

  INFO("Exit: gen_IC_skymap");

  //exit(0);

  return stat;
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// calc. of a 2D array of g-ray emission for the given E_gammagrid for a
// particular pixel (l,b)   IMOS20080114
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
int Galprop::gen_IC_skymap_pixel(const double l, const double b, 
				 //const Distribution & electrons, 
				 //const int ielectrons, 
				 vector< vector<double> >& iso_IC) { //, 
				 //vector< vector<double> >& aniso_IC, 
				 //const double Etarget[], const double factor) {
  
  double dtr=acos(-1.)/180.;
  int ir,ix,iy,iz;
  
  //IMOS20060420 full calculation of anisoIC
  double RG=galaxy.r_max;                           // kpc,radius of Galactic disk
  double DENS;
  double ISRF_over_nu[galaxy.ISRF[0].n_pgrid];
  
  double sinb=sin(b*dtr);
  double cosb=cos(b*dtr);
  double sinl=sin(l*dtr);
  double cosl=cos(l*dtr);
  double d=0;
  int complete=0;

  double dd=galdef.LoS_step /galdef.LoS_substep_number /(fabs(sinb)+1.e-6);// variable step depends on b
  if(dd>galdef.LoS_step) dd=galdef.LoS_step;                               // max integration step in kpc 

  if(galdef.verbose==-457) {// selectable debug (-457 -old method)
    
    double dd=0.1/(fabs(sinb)+1.e-6);
    if(dd>0.5) dd=0.5;
  }
  
  //defining vectors
  iso_IC.resize(galaxy.n_ISRF_components);
  //aniso_IC.resize(galaxy.n_ISRF_components);
  
  //zeroing vectors
  for (int i=0; i<galaxy.n_ISRF_components; ++i) {

    iso_IC[i].resize(galaxy.n_E_gammagrid, 0);
    //aniso_IC[i].resize(galaxy.n_E_gammagrid, 0);
  }
  
  // integration along the line of sight
  
  while(complete==0) {

    d += dd;
    double zz=d*sinb;
    double RR=sqrt(Rsun*Rsun+pow(d*cosb,2)-2.0*Rsun*d*cosb*cosl); //IMOS20080114
    double xx,yy;                                                 //IMOS20080114
    
    // used for anisotropic IC
    //double costheta=(Rsun-d*cosb*cosl)/RR;
    //if(costheta> 1.0) costheta= 1.0;
    //if(costheta<-1.0) costheta=-1.0;
    //double theta=acos(costheta);
    
    // checks if we got to the Galactic boundary in 2D, if so stop integration
    if(gcr[0].n_spatial_dimensions==2) {

      // find the nearest grid points on the LEFT side of the current point
      ir=(int)((RR-galaxy.r_min)/galaxy.dr);
      iz=(int)((zz-galaxy.z_min)/galaxy.dz);
      
      if(galdef.verbose==-457) {// selectable debug (-457 -old method, the nearest grid point)
	
	ir=(int)((RR-galaxy.r_min)/galaxy.dr + 0.5);
	iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.5);
      }
      
      // checks if we got to the Galactic boundary in 2D, if so stop integration
      if(RR>galaxy.r_max)                    complete=1;
      if(zz<galaxy.z_min || zz>galaxy.z_max) complete=1;
      
      if(ir>galaxy.n_rgrid-1) {              complete=1; ir=galaxy.n_rgrid-1; }
      
      if(iz<0               ) {              complete=1; iz=0; } 
      if(iz>galaxy.n_zgrid-1) {              complete=1; iz=galaxy.n_zgrid-1; }
      // cout<<"d RR zz ir iz "<<d<<" "<<RR<<" "<<zz<<" "<<ir<<" "<<iz<<endl;
      
    } // particle.n_spatial_dimensions==2
    
    // checks if we got to the Galactic boundary in 3D, if so stop integration
    
    if (3 == gcr[0].n_spatial_dimensions) {
    
      xx = Rsun-d*cosb*cosl; // 3D: Sun on x axis at x=+Rsun
      yy = -d*cosb*sinl; // 3D: Sun at y=0; +ve long=-ve y since Z=X^Y system
      
      if (galdef.use_symmetry==1) {	    
	
	xx=fabs(xx);
	yy=fabs(yy);
	zz=fabs(zz);
      
      }
      
      // find the nearest grid points on the LEFT side of the current point
      ix=(int)((xx - galaxy.x_min)/galaxy.dx);
      iy=(int)((yy - galaxy.y_min)/galaxy.dy);
      iz=(int)((zz - galaxy.z_min)/galaxy.dz);
      
      // checks if we got to the Galactic boundary in 3D, if so stop integration
      if(ix<0               ) { complete=1; ix=0;                }
      if(iy<0               ) { complete=1; iy=0;                }  
      if(iz<0               ) { complete=1; iz=0;                } 
      if(ix>galaxy.n_xgrid-1) { complete=1; ix=galaxy.n_xgrid-1; }
      if(iy>galaxy.n_ygrid-1) { complete=1; iy=galaxy.n_ygrid-1; }
      if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }
      
      if (zz<galaxy.z_min || zz>galaxy.z_max) complete=1;
      //	  if(fabs(zz) > zzmax                  ) complete=1;
      //  cout<<"d RR xx yy zz ix iy iz "<<d<<" "<<RR<<" "<<xx<<" "<<yy<<" "<<zz<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
      
    } //particle.n_spatial_dimensions==3
    
    for (int i_comp = 0; i_comp < galaxy.n_ISRF_components; ++i_comp) {//IMOS20060420
      
      //cout<<"dd d R theta Z i_comp IC_aniso_factor "<<dd<<" "<<d<<" "<<RR<<" "
      //<<theta<<" "<<zz<<" "<<i_comp<<" "<<IC_aniso_factor<<endl;
      
      //full calculation of anisoIC - too complicated already, so use a simple (old) method
      /*if(galdef.IC_anisotropic==1) {

	if(gcr[0].n_spatial_dimensions==2) {

	  ir=(int)((RR-galaxy.r_min)/galaxy.dr + 0.5);
	  iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.5);
	  
	  for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++)
	    ISRF_over_nu[inu] =galaxy.ISRF[i_comp].d2[ir][iz]    .s[inu]/galaxy.nu_ISRF[inu];
	}
	if(gcr[0].n_spatial_dimensions==3) {
	  ix=(int)((xx-galaxy.x_min)/galaxy.dx + 0.5);
	  iy=(int)((yy-galaxy.y_min)/galaxy.dy + 0.5);
	  iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.5);
	  
	  for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++)
	    ISRF_over_nu[inu] =galaxy.ISRF[i_comp].d3[ix][iy][iz].s[inu]/galaxy.nu_ISRF[inu];
	}
	}*/
      
      for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++) {//IMOS20060420
	
	//double IC_aniso_emiss=0., IC_aniso_factor=1.;
	double delta, x[8][3],f[8],y[7];
	
	//full calculation of anisoIC - too complicated already, so use a simple (old) method
	/*if(galdef.IC_anisotropic==1 && i_comp<2) { // i_comp=2 -CMB
	  
	  for(int ip=0; ip<gcr[ielectrons].n_pgrid; ip++) {

	    double sum=0.0;
	    for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++) {
	      
	      if(galaxy.E_gamma[iEgamma] > 4.*Etarget[inu]*1.e-6 *pow(gcr[ielectrons].gamma[ip],2)  // IMOS20030217
		 /(1.+4.*Etarget[inu]*1.e-6/Mele*gcr[ielectrons].gamma[ip])) continue;
	      
	      aic_cc(   0,1,Etarget[inu]*1.e-6/Mele,galaxy.E_gamma[iEgamma]/Mele,gcr[ielectrons].gamma[ip],RG,RR,theta,fabs(zz),Rsun, DENS);
	      
	      sum+= ISRF_over_nu[inu] 
		*aic_cc(1,1,Etarget[inu]*1.e-6/Mele,galaxy.E_gamma[iEgamma]/Mele,gcr[ielectrons].gamma[ip],RG,RR,theta,fabs(zz),Rsun, DENS);
	    }
	    
	    //cout<<"ir iz i_comp Ekin E_gamma sum cr_density "<<ir<<" "<<iz<<" "<<i_comp<<" "<<gcr[ielectrons].Ekin[ip]
	    //<<" "<<galaxy.E_gamma[iEgamma]<<" "<<sum<<" "<<gcr[ielectrons].cr_density.d2[ir][iz].s[ip]<<endl;
	    
	    sum*=factor; // to avoid loss of accuracy
	    
	    if(gcr[0].n_spatial_dimensions==2) 
	      IC_aniso_emiss+=sum *electrons.d2[ir][iz]    .s[ip]*gcr[ielectrons].Ekin[ip];
	    
	    if(gcr[0].n_spatial_dimensions==3) 
	      IC_aniso_emiss+=sum *electrons.d3[ix][iy][iz].s[ip]*gcr[ielectrons].Ekin[ip];
	  }//ip
	  
	  IC_aniso_factor=1.;
	  
	  if(gcr[0].n_spatial_dimensions==2 && galaxy.IC_iso_emiss[i_comp].d2[ir][iz]    .s[iEgamma]!=0.) 
	    IC_aniso_factor=IC_aniso_emiss/galaxy.IC_iso_emiss[i_comp].d2[ir][iz]    .s[iEgamma];
	  
	  if(gcr[0].n_spatial_dimensions==3 && galaxy.IC_iso_emiss[i_comp].d3[ix][iy][iz].s[iEgamma]!=0.) 
	    IC_aniso_factor=IC_aniso_emiss/galaxy.IC_iso_emiss[i_comp].d3[ix][iy][iz].s[iEgamma];
	}//full anisoIC
	*/
	
	//if(galdef.IC_anisotropic==2) IC_aniso_factor=IC_anisotropy_factor(galaxy.E_gamma[iEgamma],RR,theta,fabs(zz),i_comp);//IMOS20060420
	//if(galdef.IC_anisotropic==3) IC_aniso_factor=1.0;            // anisotropic=isotropic skymap  IMOS20060420
	//if(galdef.IC_anisotropic==1 || galdef.IC_anisotropic==2) 
	//cout<<"ir iz i_comp E_gamma IC_aniso_factor-> "<<ir<<" "<<iz<<" "<<i_comp<<" "<<galaxy.E_gamma[iEgamma]<<" "<<IC_aniso_factor<<endl;
	
	if(gcr[0].n_spatial_dimensions==2) {

	  if (ir == galaxy.n_rgrid-1 || 
	      iz == galaxy.n_zgrid-1 || 
	      galdef.verbose == -457) { //-457 -old method
	    
	    delta = dd*kpc2cm*galaxy.IC_iso_emiss[i_comp].d2[ir][iz].s[iEgamma];
	  
	  } else {  // linear interpolation
	    
	    //  x[0]=(R0,z1), x[1]=(R1,z1);  y[0] -location of the grid points
	    //  x[2]=(R0,z0), x[3]=(R1,z0);  y[1]
	    x[0][0] = galaxy.r[ir]; 
	    x[0][1] = galaxy.z[iz+1];  
	    f[0] = galaxy.IC_iso_emiss[i_comp].d2[ir][iz+1].s[iEgamma];
	    
	    x[1][0] = galaxy.r[ir+1]; 
	    x[1][1] = galaxy.z[iz+1];  
	    f[1] = galaxy.IC_iso_emiss[i_comp].d2[ir+1][iz+1].s[iEgamma];
	    
	    x[2][0] = galaxy.r[ir]; 
	    x[2][1] = galaxy.z[iz];  
	    f[2] = galaxy.IC_iso_emiss[i_comp].d2[ir][iz].s[iEgamma];
	    
	    x[3][0] = galaxy.r[ir+1]; 
	    x[3][1] = galaxy.z[iz];  
	    f[3] = galaxy.IC_iso_emiss[i_comp].d2[ir+1][iz].s[iEgamma];
	    
	    y[0] = (f[0] - f[1])/(x[0][0] - x[1][0])*(RR - x[0][0]) + f[0]; // interpolation in R
	    y[1] = (f[2] - f[3])/(x[2][0] - x[3][0])*(RR - x[2][0]) + f[2];
	    
	    y[2]=(y[0]-y[1])/(x[0][1]-x[2][1])*(zz-x[0][1])+y[0];   // interpolation in z
	    
	    delta = dd*kpc2cm*y[2];

	  }

	}

	if (3 == gcr[0].n_spatial_dimensions) {
 	      
	  if (ix == galaxy.n_xgrid-1 || 
	      iy == galaxy.n_ygrid-1 || 
	      iz == galaxy.n_zgrid-1) {

	    delta = dd*kpc2cm*galaxy.IC_iso_emiss[i_comp].d3[ix][iy][iz].s[iEgamma];
	  
	  } else {  // linear interpolation
	    
	    //  x[0]=(x0,z1,y0), x[1]=(x1,z1,y0);   x[4]=(x0,z1,y1), x[5]=(x1,z1,y1);
	    //  x[2]=(x0,z0,y0), x[3]=(x1,z0,y0);   x[6]=(x0,z0,y1), x[7]=(x1,z0,y1);
	    x[0][0] = galaxy.x[ix]; 
	    x[0][1] = galaxy.z[iz+1]; 
	    x[0][2] = galaxy.y[iy];  
	    f[0] = galaxy.IC_iso_emiss[i_comp].d3[ix][iy][iz+1].s[iEgamma];

	    x[1][0] = galaxy.x[ix+1]; 
	    x[1][1] = galaxy.z[iz+1]; 
	    x[1][2] = galaxy.y[iy];  
	    f[1] = galaxy.IC_iso_emiss[i_comp].d3[ix+1][iy][iz+1].s[iEgamma];

	    x[2][0] = galaxy.x[ix]; 
	    x[2][1] = galaxy.z[iz]; 
	    x[2][2] = galaxy.y[iy];  
	    f[2] = galaxy.IC_iso_emiss[i_comp].d3[ix][iy][iz].s[iEgamma];

	    x[3][0] = galaxy.x[ix+1]; 
	    x[3][1] = galaxy.z[iz]; 
	    x[3][2] = galaxy.y[iy];  
	    f[3] = galaxy.IC_iso_emiss[i_comp].d3[ix+1][iy][iz].s[iEgamma];

	    x[4][0] = galaxy.x[ix]; 
	    x[4][1] = galaxy.z[iz+1]; 
	    x[4][2] = galaxy.y[iy+1];  
	    f[4] = galaxy.IC_iso_emiss[i_comp].d3[ix][iy+1][iz+1].s[iEgamma];

	    x[5][0] = galaxy.x[ix+1]; 
	    x[5][1] = galaxy.z[iz+1]; 
	    x[5][2] = galaxy.y[iy+1];  
	    f[5] = galaxy.IC_iso_emiss[i_comp].d3[ix+1][iy+1][iz+1].s[iEgamma];

	    x[6][0] = galaxy.x[ix]; 
	    x[6][1] = galaxy.z[iz]; 
	    x[6][2] = galaxy.y[iy+1];  
	    f[6] = galaxy.IC_iso_emiss[i_comp].d3[ix][iy+1][iz].s[iEgamma];

	    x[7][0] = galaxy.x[ix+1]; 
	    x[7][1] = galaxy.z[iz]; 
	    x[7][2] = galaxy.y[iy+1];  
	    f[7] = galaxy.IC_iso_emiss[i_comp].d3[ix+1][iy+1][iz].s[iEgamma];
	    
	    //f[0] = f[1] = f[2] = f[3] = f[4] = f[5] = f[6] = f[7] = 1;

	    y[0] = (f[0] - f[1])/(x[0][0] - x[1][0])*(xx - x[0][0]) + f[0]; // interpolation in x
	    y[1] = (f[2] - f[3])/(x[2][0] - x[3][0])*(xx - x[2][0]) + f[2];
	    y[2] = (f[4] - f[5])/(x[4][0] - x[5][0])*(xx - x[4][0]) + f[4];
	    y[3] = (f[6] - f[7])/(x[6][0] - x[7][0])*(xx - x[6][0]) + f[6];
	    
	    y[4] = (y[0] - y[1])/(x[0][1] - x[2][1])*(zz - x[0][1]) + y[0];   // interpolation in z
	    y[5] = (y[2] - y[3])/(x[4][1] - x[6][1])*(zz - x[4][1]) + y[2];
	    
	    y[6] = (y[4] - y[5])/(x[0][2] - x[4][2])*(yy - x[0][2]) + y[4];   // interpolation in y
	    
	    //		  cout<<" y= "; for(int j=0;j<7;j++) cout<<y[j]<<" "; cout<<" >> "<<xx+zz+yy<<endl; 
	    //		  for(int j=0;j<3;j++) {cout<<" x,f= "; for(int i=0;i<8;i++) cout<<x[i][j]<<" "<<f[i]<<"   "; cout<<endl;} exit(0);
	    
	    delta = dd*kpc2cm*y[6];

	  }

	}

	iso_IC[i_comp][iEgamma] += delta;
	//if(galdef.IC_anisotropic) aniso_IC[i_comp][iEgamma] +=delta*IC_aniso_factor;
	
	//cout<<"ir iz i_comp  E_gamma IC_iso_emiss "<<ir<<" "<<iz<<" "<<i_comp<<" "<<" "
	//<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.IC_iso_emiss[i_comp].d2[ir][iz].s[iEgamma]<<endl;     
      }//iEgamma
    
    }//ISRF_components
  
  }//complete==0
  
  return 0;

}
