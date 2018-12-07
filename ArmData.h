#ifndef _rf_SpiralArmData_h_
#define _rf_SpiralArmData_h_

#include <iostream>

namespace rf {

  class ArmData {

  public:

    ArmData() {}
    ArmData(const double a,
	    const double rMin,
	    const double rMax,
	    const double phiMin,
	    const double phiExtent,
	    const double width) :
      fA(a), fRMin(rMin), fRMax(rMax), fPhiMin(phiMin), fPhiExtent(phiExtent),
      fWidth(width) {}

      double fA, fRMin, fRMax, fPhiMin, fPhiExtent, fWidth;

      // return angle in radians

      inline const double Angle(const double r) const {

	return fA*std::log(r/fRMin) + fPhiMin;

      }

      friend
	std::ostream& operator<<(std::ostream& out, const ArmData& data) {

	return out << "("
		   << data.fA << " "
		   << data.fRMin << " "
		   << data.fRMax << " "
		   << data.fPhiMin << " "
		   << data.fPhiExtent << " "
		   << data.fWidth << ")";

      }

  };

}

#endif
