#ifndef _DistributionFunction_h_
#define _DistributionFunction_h_

#include "Distribution.h"
#include "los_integration.h"
#include <valarray>
#include <stdexcept>

class DistributionFunction : public SM::LOSfunction<std::valarray<double> > {
   private:
      const Distribution &fdist;
      const std::vector<double> fx, fy, fz, fr;
   public:
      DistributionFunction(const Distribution &dist, const std::vector<double> &z, const std::vector<double> &r);
      DistributionFunction(const Distribution &dist, const std::vector<double> &z, const std::vector<double> &x, const std::vector<double> &y);

      //Use bilinear interpolation to estimate value in x, y, z
      virtual std::valarray<double> operator () ( double x, double y, double z ) const;
};

#endif
