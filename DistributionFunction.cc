#include "DistributionFunction.h"
#include "Distribution.h"

#include <valarray>
#include <cmath>
#include <stdexcept>

DistributionFunction::DistributionFunction(const Distribution &dist, const std::vector<double> &z, const std::vector<double> &r) : 
   fdist(dist),
   fz(z),
   fr(r)
{
   //Check for dimensions in distribution
   if ( dist.n_spatial_dimensions != 2 )
      throw(std::invalid_argument ( "Distribution must have 2 spatial dimensions when DistributionFunction is initialized with radial boundaries" ));
   if ( dist.n_rgrid != r.size() )
      throw(std::invalid_argument ( "DistributionFunction: Size of r array must equal size of r grid in Distribution" ));
   if ( dist.n_zgrid != z.size() )
      throw(std::invalid_argument ( "DistributionFunction: Size of z array must equal size of z grid in Distribution" ));
}

DistributionFunction::DistributionFunction(const Distribution &dist, const std::vector<double> &z, const std::vector<double> &x, const std::vector<double> &y) : 
   fdist(dist),
   fz(z),
   fx(x),
   fy(y)
{
   //Check for dimensions in distribution
   if ( dist.n_spatial_dimensions != 3 )
      throw(std::invalid_argument ( "Distribution must have 3 spatial dimensions when DistributionFunction is initialized with Cartesian boundaries" ));
   if ( dist.n_xgrid != x.size() )
      throw(std::invalid_argument ( "DistributionFunction: Size of x array must equal size of x grid in Distribution" ));
   if ( dist.n_ygrid != y.size() )
      throw(std::invalid_argument ( "DistributionFunction: Size of y array must equal size of y grid in Distribution" ));
   if ( dist.n_zgrid != z.size() )
      throw(std::invalid_argument ( "DistributionFunction: Size of z array must equal size of z grid in Distribution" ));
}


std::valarray<double> DistributionFunction::operator () ( double x, double y, double z ) const {

   std::valarray<double> output(fdist.n_pgrid);

   //We always interpolate in z, 
   size_t iz(0);
   //-2 as we want iz+1 to be a valid index after the loop.
   for ( ; iz < fz.size()-2; ++iz )
      if ( fz[iz] <= z && z < fz[iz+1] )
	 break;
   double lzf = (fz[iz+1]-z)/(fz[iz+1]-fz[iz]);
   double uzf = 1-lzf;

   //Check for the dimensions of the distribution, to find how to interpolate
   if ( fdist.n_spatial_dimensions == 2 ) {
      double r = sqrt(x*x+y*y);
      size_t ir(0);
      for ( ; ir < fr.size()-2; ++ir )
	 if ( fr[ir] <= r && r < fr[ir+1] )
	    break;
      double lrf = (fr[ir+1]-r)/(fr[ir+1]-fr[ir]);
      double urf = 1-lrf;

      for ( size_t i = 0; i < fdist.n_pgrid; ++i ) 
	 output[i] = lzf*lrf*fdist.d2[ir  ][iz  ].s[i]
	           + lzf*urf*fdist.d2[ir+1][iz  ].s[i]
		   + uzf*lrf*fdist.d2[ir  ][iz+1].s[i]
		   + uzf*urf*fdist.d2[ir+1][iz+1].s[i];
   } else if ( fdist.n_spatial_dimensions == 3 ) {
      size_t ix(0);
      for ( ; ix < fx.size()-2; ++ix )
	 if ( fx[ix] <= x && x < fx[ix+1] )
	    break;
      double lxf = (fx[ix+1]-x)/(fx[ix+1]-fx[ix]);
      double uxf = 1-lxf;

      size_t iy(0);
      for ( ; iy < fy.size()-2; ++iy )
	 if ( fy[iy] <= y && y < fy[iy+1] )
	    break;
      double lyf = (fy[iy+1]-y)/(fy[iy+1]-fy[iy]);
      double uyf = 1-lyf;
      for ( size_t i = 0; i < fdist.n_pgrid; ++i ) 
	 output[i] = lxf*lyf*lzf*fdist.d3[ix  ][iy  ][iz  ].s[i]
	           + lxf*lyf*uzf*fdist.d3[ix  ][iy  ][iz+1].s[i]
		   + lxf*uyf*lzf*fdist.d3[ix  ][iy+1][iz  ].s[i]
		   + lxf*uyf*uzf*fdist.d3[ix  ][iy+1][iz+1].s[i]
	           + uxf*lyf*lzf*fdist.d3[ix+1][iy  ][iz  ].s[i]
	           + uxf*lyf*uzf*fdist.d3[ix+1][iy  ][iz+1].s[i]
		   + uxf*uyf*lzf*fdist.d3[ix+1][iy+1][iz  ].s[i]
		   + uxf*uyf*uzf*fdist.d3[ix+1][iy+1][iz+1].s[i];
   }

   return output;
}
