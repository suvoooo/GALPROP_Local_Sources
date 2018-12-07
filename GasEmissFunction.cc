#include "Galprop.h"
#include "DistributionFunction.h"
#include <valarray>
#include <stdexcept>

Galprop::GasEmissFunction::GasEmissFunction( const std::string &type, const std::string &gas_type, Galprop&gp ) :
   fgp(gp),
   fgf(gas_type, 0.0, gp),
   fdf(0)
{
   //Set up arrays, depending on dimension of galprop
   std::vector<double> z(fgp.galaxy.z, fgp.galaxy.z+fgp.galaxy.n_zgrid);
   if ( type == "BREMSS" || type == "Bremss" || type == "bremss" ) {
      ftype = BREMSS;
      if ( fgp.galaxy.n_spatial_dimensions == 2 ) {
	 std::vector<double> r(fgp.galaxy.r, fgp.galaxy.r+fgp.galaxy.n_rgrid);
	 fdf = new DistributionFunction(fgp.galaxy.bremss_emiss, z, r);
      } else if ( fgp.galaxy.n_spatial_dimensions == 3 ) {
	 std::vector<double> x(fgp.galaxy.x, fgp.galaxy.x+fgp.galaxy.n_xgrid);
	 std::vector<double> y(fgp.galaxy.y, fgp.galaxy.y+fgp.galaxy.n_ygrid);
	 fdf = new DistributionFunction(fgp.galaxy.bremss_emiss, z, x, y);
      }
   } else if ( type == "PION" || type == "Pion" || type == "pion" || type == "PI0" || type == "Pi0" || type == "pi0" ) {
      ftype = PION;
      if ( fgp.galaxy.n_spatial_dimensions == 2 ) {
	 std::vector<double> r(fgp.galaxy.r, fgp.galaxy.r+fgp.galaxy.n_rgrid);
	 fdf = new DistributionFunction(fgp.galaxy.pi0_decay_emiss, z, r);
      } else if ( fgp.galaxy.n_spatial_dimensions == 3 ) {
	 std::vector<double> x(fgp.galaxy.x, fgp.galaxy.x+fgp.galaxy.n_xgrid);
	 std::vector<double> y(fgp.galaxy.y, fgp.galaxy.y+fgp.galaxy.n_ygrid);
	 fdf = new DistributionFunction(fgp.galaxy.pi0_decay_emiss, z, x, y);
      }
   } else {
      throw(std::invalid_argument("Type can only be \"PION\" or \"BREMSS\" in GasEmissFunction"));
   }
}

Galprop::GasEmissFunction::~GasEmissFunction() {
   delete fdf;
}

std::valarray<double> Galprop::GasEmissFunction::operator () (double x, double y, double z) const {
   //Just multiply the gas function with the distribution function, easy peasy
   return (*fdf)(x,y,z)*fgf(x,y,z);
}
