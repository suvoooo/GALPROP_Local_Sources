#include "Galprop.h"
#include "galprop_internal.h"
#include <stdexcept>
#include <cmath>

Galprop::GasFunction::GasFunction( const std::string &type, double rPLindex, Galprop& gp) :
   fgp(gp),
   frInd(rPLindex)
{
   if ( type == "HI" )
      ftype = HI;
   else if ( type == "H2" )
      ftype = H2;
   else if ( type == "CO" )
      ftype = CO;
   else if ( type == "HII" )
      ftype = HII;
   else
      throw(std::invalid_argument("Type can only be \"HI\", \"H2\", \"CO\", or \"HII\" in GasFunction"));
}

double Galprop::GasFunction::operator () (double x, double y, double z) const{
   const double r = sqrt(x*x + y*y);
   const double rScale = pow(r, frInd);
   switch (ftype) {
      case HI:
	 return nHI( r, z )*rScale;
      case CO:
	 return 2*nH2( r, z, 1.0 )*rScale;
      case H2:
	 return 2*nH2( r, z, fgp.fX_CO(r) )*rScale;
      case HII:
	 return nHII( r, z )*rScale;
      default:
	 return 0;
   }
}
