//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * source_distribution.cc *                      galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include"galprop_classes.h"
#include"galprop_internal.h"

#include <cmath>
#include <vector>

using namespace std;

#include <StellarGeometryDistributions.h>

static vector<rf::ArmData> gArmData;

// cosmic ray source distribution: x,y,z in kpc

double Galprop::source_distribution (const double x, 
				     const double y, 
				     const double z, 
				     int srcModel, const std::vector<double> &parameters) {
   
  const double r=sqrt(x*x + y*y);
  double alpha,beta,rmax,ro,rc,zscale,Value;
  double rm,rs,rconst,roff;                               //AWS20091001,20091008

  //double xx = x; //YO/20140316

  const double rSun = 8.5; // kpc
  
  double result = 0;
  
  // test of electron propagation vs analytical calculations IMOS20061030
  if(abs(galdef.DM_int0)==99) 
    {
      if(galdef.DM_double6 < fabs(z)) return 0.;
      else return 1.;
    }
  // end of the test area
  
  if (1 == srcModel)  //arbitrary parametrization
    {         
      alpha= parameters[1];
      beta = parameters[2]; 
      rmax = parameters[3]; 
      rconst=parameters[4]; // for rconst<r<rmax, set to value at rconst.     in kpc from GC  AWS20091008
      roff  =parameters[5]; // replace r with x=r+roff and r0 with x0=r0+roff in the equations (Yusifov & Kucuc 2004).   in kpc  GJ20100223
      
      ro=8.5;
      zscale=0.200;  //nominal value
      double x      = r      + roff;
      double xo     = ro     + roff;
      double xconst = rconst + roff;
      
      //cerr << "x= " << xx << " y= " << y << " z= " << z << endl;
      result =pow(x     /xo,alpha) * exp(-beta*(x     -xo)/xo) * exp(-fabs(z)/zscale);
      if (r >= rconst)result =pow(xconst/xo,alpha) * exp(-beta*(xconst-xo)/xo) * exp(-fabs(z)/zscale);//AWS20091008
      
      if (r >= rmax)  result = 0.0;

      /*if (galdef.n_spatial_dimensions == 3 && galdef.rm_nearbysource == 1){ //YO/20140314
	if (sqrt((xx-8.5)*(xx-8.5) + y*y + z*z) <= galdef.near_distance){
	  result = 0.0;
	  if(galdef.check_num==314){
	    cerr << "distribution= " << result << " in x= " << xx << " y= " << y << " z= " << z << endl;
	  }
	}
	}*/
      
    }
  
  if (2 == srcModel)  //SNR:  Case and Bhattacharya 3rd Compton p437
    {
      alpha=1.69;
      beta=3.33;
      ro=8.5;
      zscale=0.200;  // nominal value;
      result =pow(r/ro,alpha) * exp(-beta*(r-ro)/ro) * exp(-fabs(z)/zscale);
    }
  
  if (3 == srcModel)  //pulsars: Taylor, Manchester and Lyne 1993 ApJSupp 88,259
    {
      rc=3.5;
      zscale=0.200; //nominal value, not for pulsars
      result =cosh(ro/rc)/cosh(r/rc)  * exp(-fabs(z)/zscale);
    }    
  
  
  if (5 == srcModel) //Strong and Mattox gamma-ray distribution E> 100 MeV
    {
      zscale=0.200;
      Value=19.3;
      if (r >  4.0)  Value=21.9;
      if (r >  8.0)  Value=15.8;
      if (r > 10.0)  Value=18.3;
      if (r > 12.0)  Value=13.3;
      if (r > 15.0)  Value= 7.4;
      result = Value * exp(-fabs(z)/zscale);
    }
  
  if (6 == srcModel)  //Strong and Mattox gamma-ray distribution E> 100 MeV R<15 kpc, zero beyond
    {
      zscale=0.200;
      Value=19.3;
      if (r >  4.0)  Value=21.9;
      if (r >  8.0)  Value=15.8;
      if (r > 10.0)  Value=18.3;
      if (r > 12.0)  Value=13.3;
      if (r > 15.0)  Value= 0.0;
      result =  Value  * exp(-fabs(z)/zscale);
    }         
  
  
  if (7 == srcModel)  // Gaussian                AWS20091001
    {         
      rm   = parameters[1]; // mean  of Gaussian in kpc from GC
      rs   = parameters[2]; // sigma of Gaussian in kpc
      rmax = parameters[3]; // cutoff radius     in kpc from GC
      rconst=parameters[4]; // for rconst<r<rmax, set to value at rconst.     in kpc from GC  AWS20091008
      
      
      zscale=0.200;  //nominal value
      
      result =   exp(-(r     -rm)*(r     -rm)/(2.*rs*rs)) * exp(-fabs(z)/zscale);
      if (r >= rconst)result =   exp(-(rconst-rm)*(rconst-rm)/(2.*rs*rs)) * exp(-fabs(z)/zscale); //AWS20091008
      if (r >= rmax  )result =   0.0;
    }
  
  if (8 == srcModel) // Linear interpolation from tabulated values
    {
      //Find correct bin for interpolation (extrapolation if needed)
      if (r <= galdef.source_radius[0]) {
	result = galdef.source_values[0];
      } else if ( r >= galdef.source_radius[galdef.n_source_values-1]) {
	result = galdef.source_values[galdef.n_source_values-1];
      } else {
	int i = 0;
	while (r > galdef.source_radius[i])
	  i++;
	result = galdef.source_values[i-1] + 
	  (galdef.source_values[i]-galdef.source_values[i-1])/(galdef.source_radius[i]-galdef.source_radius[i-1])*
	  (r-galdef.source_radius[i-1]);
      }
      zscale=0.200;  //nominal value
      result *= exp(-fabs(z)/zscale);
    }
  
  if (9 == srcModel) //Total gas distribution
    {
      //Simply sum up the analytial gas functions
      result = nHI(r,z) + 2*nH2(r,z,fX_CO(r));
    }
  
  if (10 == srcModel) //Only CO
    {
      result = nH2(r,z,fX_CO(r));
    }
  
  if (11 == srcModel) {
    
    result = nHII(r, z);
    
  }
  
  if (12 == srcModel) {

    // Corresponds to OB stars in arms + thin disc from the ISRF model

    const double rMin = 0., rMax = 15.;
    const double zMin = -2., zMax = 2.;
    const double phiMin = 0., phiMax = utl::kTwoPi;
 
    const double discDensity = parameters[2], armDensity = parameters[3];
    const double densityOStar = 1e1, zScaleOStar = 0.035;
    const double densityBStar = 1e2, zScaleBStar = 0.2;
    const double rScale = parameters[1];
    const double hole = 1.3, holeIndex = 1.7;

    const double discContribution = 
      discDensity*exp(-(r - rSun)/rScale)* 
      (1. - exp(-pow(r/hole, holeIndex)))*
      (//densityOStar*exp(-fabs(z)/zScaleOStar) +
       densityBStar*exp(-fabs(z)/zScaleBStar));

    static double armNormalisation = 0.;
    
    if (0 == gArmData.size()) {

      rf::ArmData arm1(4.45, 3., 15., -0.6, 6., 0.4);
      rf::ArmData arm2(4.45, 3., 15., 2.55, 6., 0.4);
      rf::ArmData arm3(4.45, 4., 15., -1.2, 6., 0.4);
      rf::ArmData arm4(4.45, 4., 15., 1.95, 6., 0.4);
      rf::ArmData arm5(3.2, 8.1, 12., 2.71, 0.55, 0.4);

      gArmData.push_back(arm1);
      gArmData.push_back(arm2);
      gArmData.push_back(arm3);
      gArmData.push_back(arm4);
      gArmData.push_back(arm5);

      const int rBins = 600, zBins = 80, phiBins = 72; 

      const double rDelta = fabs(rMax - rMin)/rBins;
      const double zDelta = fabs(zMax - zMin)/zBins;
      const double phiDelta = fabs(phiMax - phiMin)/phiBins;

      double discSum = 0, armSum = 0;

      // Get normalisation based on relative surface area of disc vs. arms

      for (int iR = 0; iR < rBins; ++iR) {

	const double r = (iR+0.5)*rDelta;

	for (int iPhi = 0; iPhi < phiBins; ++iPhi) {
	  
	  const double phi = (iPhi+0.5)*phiDelta;
	  
	  double arm = 0;
	  
	  for (vector<rf::ArmData>::const_iterator aIt = gArmData.begin();
	       aIt != gArmData.end(); ++aIt) 
	    arm += utl::CompositeArm(r, rScale, rSun, phi, 0., 1., *aIt);
	  
	  armSum += arm*r*rDelta*phiDelta;
	  
	  discSum += exp(-(r - rSun)/rScale)*(1. - exp(-pow(r/hole, holeIndex)))*r*rDelta*phiDelta;
	  
	}

      }
      
      armNormalisation = discSum/armSum;
     
    }
    
    double arm = 0;

    if (armNormalisation > 0.) {

      const double phi = atan2(y, x);
      
      for (vector<rf::ArmData>::const_iterator aIt = gArmData.begin();
	   aIt != gArmData.end(); ++aIt) 
	arm += (2 == galdef.n_spatial_dimensions ? 
		densityOStar*
		utl::CompositeArmAverage(r, rScale, rSun, z, zScaleOStar, *aIt,
					 20.*utl::kPi/180.) + 
		densityBStar*
		utl::CompositeArmAverage(r, rScale, rSun, z, zScaleBStar, *aIt,
					 20.*utl::kPi/180.) :
		densityOStar*
		utl::CompositeArm(r, rScale, rSun, phi, z, zScaleOStar, *aIt) +
		densityBStar*
		utl::CompositeArm(r, rScale, rSun, phi, z, zScaleBStar, *aIt));
      
    }

    const double armContribution = armDensity*arm*armNormalisation;
   
    result = discContribution;// + armContribution;
    
  }

  if (13 == srcModel) {

    // Corresponds to bulge + thin disc from the ISRF model
 
    const double discDensity = parameters[0];
    const double rScale = parameters[1];
    const double zScale = parameters[2];
    const double bulgeDensity = parameters[3];
    const double bulgeA = parameters[4];
    const double bulgeB = parameters[5];
    const double bulgeScaleLength = parameters[6];
    const double bulgeIndex = parameters[7];
    const double phiOffset = parameters[8]*utl::kConvertDegreesToRadians;
    
    const double rS = 8.5;
    const double hole = 1.3, holeIndex = 1.7;
    
    const double discContribution = 
      discDensity*exp(-(r - rSun)/rScale)* 
      (1. - exp(-pow(r/hole, holeIndex)))*exp(-fabs(z)/zScale);
  
    const double phi = atan2(y, x);

    const double bulgeContribution =
      (2 == galdef.n_spatial_dimensions ?
       utl::LopezCorredoiraBulgeAverage(r, z, phiOffset, 1., bulgeScaleLength, bulgeA, bulgeB, bulgeIndex, utl::kPi/8) : 
       utl::LopezCorredoiraBulge(r, phi, z, phiOffset, 1., bulgeScaleLength, bulgeA, bulgeB, bulgeIndex));
     
    result = discContribution + bulgeContribution;
    
    //cout << x << " " << y << " " << r << " " << z << " " << discContribution << " " << bulgeContribution << " " << result << endl;

  }

  if (14 == srcModel) {

    // Corresponds to parameterised model + ellipsoid centred on GC.

    const double bulgeDensity = parameters[0];
    const double bulgeA = parameters[1];
    const double bulgeB = parameters[2];
    const double bulgeScaleLength = parameters[3];
    const double bulgeIndex = 2.;
    const double phiOffset = 0.0*utl::kConvertDegreesToRadians;; // deg -> radians
    
    const double phi = atan2(y, x);

    const double bulge =
      bulgeDensity*(2 == galdef.n_spatial_dimensions ?
       utl::LopezCorredoiraBulgeAverage(r, z, phiOffset, 1., bulgeScaleLength, bulgeA, bulgeB, bulgeIndex, utl::kPi/8) : 
       utl::LopezCorredoiraBulge(r, phi, z, phiOffset, 1., bulgeScaleLength, bulgeA, bulgeB, bulgeIndex));
 
    const double alpha = parameters[4];
    const double beta = parameters[5];
    const double rMax = parameters[6];
    const double rConst = parameters[7];
    const double rOffset = parameters[8];

    const double zScale = parameters[9];

    const double x = r + rOffset;
    const double xS = rSun + rOffset;
    const double xConst = rConst + rOffset;
      
    double parameterised = pow(x/xS, alpha)*exp(-beta*(x - xS)/xS)*exp(-fabs(z)/zScale);

    if (r >= rConst)
      parameterised = pow(xConst/xS, alpha)*exp(-beta*(xConst - xS)/xS)*exp(-fabs(z)/zScale);

    if (r >= rMax)
      parameterised = 0;

    result = bulge + parameterised;

    //cout << r << " " << z << " " << bulge << " " << parameterised << endl;
    
  }

  if (3 == galdef.n_spatial_dimensions) {
     
    for (int i_cr_source=0; i_cr_source<galdef.n_cr_sources; i_cr_source++) {

      const double r2 = 
	pow(galdef.cr_source_x[i_cr_source] - x, 2.) + 
	pow(galdef.cr_source_y[i_cr_source] - y, 2.) +
	pow(galdef.cr_source_z[i_cr_source] - z, 2.);

      const double src = galdef.cr_source_L[i_cr_source]*exp(-r2/(2.*pow(galdef.cr_source_w[i_cr_source], 2.)));

      result += src;
// cout<<" source_distribution: x y z cr source r2 L s:"<<x<<" "<<y<<" "<<z<<" "<<i_cr_source<<" "<<r2<<" "<<s<<endl;
    }
  
  }   //  galdef.n_spatial_dimensions==3
// cout<<"source distribution R z model source "<< r <<" "<<z<<" "<<galdef.source_model<<" "<<result<<endl;

   return result;
}
