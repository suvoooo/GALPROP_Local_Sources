#ifndef _utl_Distributions_h_
#define _utl_Distributions_h_

#include <PhysicalConstants.h>
#include <ArmData.h>

#include <cmath>

namespace utl {

  static double ExpDiscWithHole(const double r, 
				const double rScale, 
				const double rCut, 
				const double rSun, 
				const double z, 
				const double zScale,
				const double rHole,
				const double holeIndex) {

    using namespace std;

    const double expFact = exp(-(r - rSun)/rScale - fabs(z)/zScale);//(r <= rCut ? rScale : 0.5) - 
    //fabs(z)/zScale);

    return expFact*(1. - exp(-pow(r/rHole, holeIndex)));

  }

  static double ExpDiscWithHoleComposite(const double r,
					 const double rScale,
					 const double rCut,
					 const double rSun,
					 const double z,
					 const double zScale,
					 const double rHole,
					 const double holeIndex) {

    using namespace std;

    const double expFact = exp(-(r - rSun)/(r <= rCut ? rScale : 0.5));

    return expFact*(1. - exp(-pow(r/rHole, holeIndex)))/
      pow(cosh(-fabs(z)/zScale), 2.);

  }

  static double ExpDisc(const double r,
			const double rScale,
			const double rCut,
			const double rSun,
			const double z,
			const double zScale) {

    using namespace std;

    return (r <= rCut ? exp(-(r - rSun)/rScale - fabs(z)/zScale) : 0);

  }

  static double WainscoatBulge(const double r,
			       const double rMin,
			       const double r1,
			       const double z,
			       const double k1) {

    using namespace std;

    const double rs = (r < rMin ? rMin : r), zs = k1*z,
      x = sqrt(rs*rs + zs*zs)/r1;

    return exp(-x*x*x)/pow(x, 1.8);

  }

  static double WainscoatHalo(const double r,
			      const double rE,
			      const double z,
			      const double zE) {

    using namespace std;

    const double zp = zE*z;

    const double alpha = sqrt(r*r + zp*zp)/rE;

    return pow(10., -3.3307*(pow(alpha, 0.25) - 1.));

  }

  static double SDSSHalo(const double r, const double rScale, 
			 const double z, const double alpha, 
			 const double qH) {

    using namespace std;

    const double zOnQH = z/qH, rp = (r < 0.1 ? 0.1 : r);

    return pow(sqrt(rp*rp + zOnQH*zOnQH)/rScale, -alpha);

  }

  static double WainscoatRing(const double r, 
			      const double rRing,
			      const double rSigma,
			      const double rSun,
			      const double z,
			      const double zScale) {

    using namespace std;

    const double r0 = rRing*rSun, rS = 0.15*rSun;//2.0*rSigma*rSun;

    return exp(-pow((r - r0)/rS, 2.) - fabs(z)/zScale);

  }

  static double LopezCorredoiraBulge(const double r, 
				     const double phi,
				     const double z,
				     const double phiOffset,
				     const double k,
				     const double scaleLength,
				     const double a,
				     const double b,
				     const double index) {
    
    using namespace std;

    // Rotate coordinates into 'bulge' frame

    const double x = r*cos(phi), y = r*sin(phi),
      cosPhiOffset = cos(phiOffset), sinPhiOffset = sin(phiOffset), 
      xp = cosPhiOffset*x + sinPhiOffset*y, 
      yp = -sinPhiOffset*x + cosPhiOffset*y;

    const double xpIndex = pow(xp, index), ypIndex = pow(yp/a, index), 
      zpIndex = pow(z/b, index);

    const double t = pow(xpIndex + ypIndex + zpIndex, 1./index);

    return k*exp(-t/scaleLength);

  }

  static double LopezCorredoiraBulgeAverage(const double r,
					    const double z,
					    const double phiOffset,
					    const double k,
					    const double scaleLength,
					    const double a,
					    const double b,
					    const double index,
					    const double dPhi) {

    using namespace std;

    const unsigned int phiSteps = int(utl::kPi/dPhi) + 1;

    double result = 0;

    for (unsigned int i = 0; i < phiSteps; ++i) {

      const double phi = i*dPhi;

      result += LopezCorredoiraBulge(r, phi, z, phiOffset, k, scaleLength, a, b, index);

    }

    return result/phiSteps;

  }

  // This should more properly be referred to as a form for a triaxial bulge
  // c.f. the L-C boxy bulge above

  static double FreudenreichBar(const double r,
				const double barX,
				const double barY,
				const double barPerp,
				const double barREnd,
				const double phi,
				const double z,
				const double barZ,
				const double barPara,
				const double barHEnd, 
				const double phiOffset) {

    using namespace std;

    // Rotate current x and y into bar frame

    const double x = r*cos(phi), y = r*sin(phi);

    const double cosPhiOffset = cos(phiOffset), sinPhiOffset = sin(phiOffset);

    const double xp = cosPhiOffset*x + sinPhiOffset*y;
    const double yp = -sinPhiOffset*x + cosPhiOffset*y;

    const double rPerp = pow(pow(fabs(xp)/barX, barPerp) + 
			     pow(fabs(yp)/barY, barPerp), 1./barPerp),
      rs = pow(pow(rPerp, barPara) + pow(fabs(z)/barZ, barPara), 1./barPara);

    return (r <= barREnd ? 1. : exp(-pow((rs - barREnd)/barHEnd, 2.)))/
      pow(cosh(rs), 2.);

  }

  static double FreudenreichBarAverage(const double r,
				       const double barX,
				       const double barY,
				       const double barPerp,
				       const double barREnd,
				       const double z,
				       const double barZ,
				       const double barPara,
				       const double barHEnd,
				       const double dPhi) {

    using namespace std;

    const unsigned int phiSteps = int(utl::kPi/dPhi) + 1;

    double result = 0.0;

    for (unsigned int i = 0; i < phiSteps; ++i) {

      const double phi = i*dPhi;

      result += FreudenreichBar(r, barX, barY, barPerp, barREnd, 
				phi, z, barZ, barPara, barHEnd, 0.0);

    }

    return result/phiSteps;

  }

  static double WainscoatArm(const double r,
			     const double rScale,
			     const double rSun,
			     const double phi, 
			     const double z,
			     const double zScale,
			     const rf::ArmData& arm) {

    using namespace std;

    double result = 0;

    if (r >= arm.fRMin && r <= arm.fRMax) {
 
      double armPhi = arm.Angle(r);

      if (armPhi <= arm.fPhiMin + arm.fPhiExtent) {
	
	armPhi += kPi;

	const double x = r*cos(phi), y = r*sin(phi), 
	  armX = r*cos(armPhi), armY = r*sin(armPhi), 
	  dX = armX - x, dY = armY - y, s = sqrt(dX*dX + dY*dY);

	result = 
	  (s <= arm.fWidth ? exp(-(r - rSun)/rScale - fabs(z)/zScale) : 0.);

      }

    }

    return result;

  }

  static double CompositeArm(const double r,
			     const double rScale,
			     const double rSun,
			     const double phi, 
			     const double z,
			     const double zScale,
			     const rf::ArmData& arm,
			     double armScale = 1.) {

    using namespace std;

    double result = 0;

    if (r >= arm.fRMin && r <= arm.fRMax) {
 
      const double armPhi = arm.Angle(r);

      if (armPhi <= arm.fPhiMin + arm.fPhiExtent) {
	
	const double x = r*cos(phi), y = r*sin(phi), 
	  armX = r*cos(armPhi), armY = r*sin(armPhi), 
	  dX = armX - x, dY = armY - y, s = sqrt(dX*dX + dY*dY),
	  armSigma = 0.5*arm.fWidth/utl::kSqrtTwo/utl::kLogTwo*armScale;

	result = exp(-0.5*pow(s/(utl::kSqrtTwo*armSigma), 2.))*exp(-(r - rSun)/rScale - fabs(z)/zScale);

      }

    }

    return result;

  }

  static double WainscoatArmAverage(const double r,
				    const double rScale,
				    const double rSun,
				    const double z,
				    const double zScale,
				    const rf::ArmData& arm,
				    const double dPhi) {

    const unsigned int phiSteps = int(utl::kTwoPi/dPhi) + 1;

    double result = 0.0;

    for (unsigned int i = 0; i < phiSteps; ++i) {
      
      const double phi = i*dPhi;

      result += WainscoatArm(r, rScale, rSun, phi, z, zScale, arm);

    }

    return result/phiSteps;

  }
    
  static double CompositeArmAverage(const double r,
				    const double rScale,
				    const double rSun,
				    const double z,
				    const double zScale,
				    const rf::ArmData& arm,
				    const double dPhi, 
				    double armScale = 1.) {

    const unsigned int phiSteps = int(utl::kTwoPi/dPhi) + 1;

    double result = 0.0;

    for (unsigned int i = 0; i < phiSteps; ++i) {
      
      const double phi = i*dPhi;

      result += CompositeArm(r, rScale, rSun, phi, z, zScale, arm);

    }

    return result/phiSteps;

  }
    
}

#endif
