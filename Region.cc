#include "Region.h"
#include "Coordinate.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cctype>

SkyRegion::SkyRegion() {
   ftype = DISC;
   fcenter = SM::Coordinate(0,0);
   fradius = 0;
}
SkyRegion::SkyRegion(RegionType type, const SM::Coordinate &center, double radius, double deltal, double deltab) :
   ftype(type),
   fcenter(center)
{
   if (radius >= 0 && radius <= 180) {
      fradius = radius;
   } else {
      std::cerr<<"Radius should be between 0 and 180"<<std::endl;
      if (radius < 0) {
	 std::cerr<<"Value set to 0"<<std::endl;
	 fradius = 0;
      } else {
	 std::cerr<<"Value set to 180"<<std::endl;
	 fradius = 180;
      }
   }
   if (deltal >= 0 && radius <= 360) {
      fdeltal = radius;
   } else {
      std::cerr<<"Deltal should be between 0 and 360"<<std::endl;
      if (deltal < 0) {
	 std::cerr<<"Value set to 0"<<std::endl;
	 fdeltal = 0;
      } else {
	 std::cerr<<"Value set to 360"<<std::endl;
	 fdeltal = 360;
      }
   }
   if (deltab >= 0 && radius <= 360) {
      fdeltab = radius;
   } else {
      std::cerr<<"Deltab should be between 0 and 360"<<std::endl;
      if (deltab < 0) {
	 std::cerr<<"Value set to 0"<<std::endl;
	 fdeltab = 0;
      } else {
	 std::cerr<<"Value set to 360"<<std::endl;
	 fdeltab = 360;
      }
   }
}

const RegionType & SkyRegion::type() const { return ftype; }
const SM::Coordinate & SkyRegion::center() const { return fcenter; }
const double & SkyRegion::radius() const { return fradius; }
const double & SkyRegion::deltal() const { return fdeltal; }
const double & SkyRegion::deltab() const { return fdeltab; }

std::istream & operator >> (std::istream &is, SkyRegion &region) {
   //We have skipws flag everywhere.  I should check the flag initially, but
   //this is easier
   std::string type;
   is >> std::skipws >> type;
   //Convert type to all uppercase and check DISC or RECTANGLE.  Otherwise set
   //the fail bit of the is
   std::transform(type.begin(), type.end(), type.begin(), toupper);
   if (type == "DISC") {
      region.ftype = DISC;
   } else if (type == "RECTANGLE") {
      region.ftype = RECTANGLE;
   } else {
      is.setstate(std::ios::failbit);
      std::cerr<<"Region type \""<<type<<"\" unknown!"<<std::endl;
      std::cerr<<"Allowed types:"<<std::endl;
      std::cerr<<"\tDISC"<<std::endl;
      std::cerr<<"\tRECTANGLE"<<std::endl;
      return is;
   }

   //Read the coordinates
   double l, b;
   is >> std::skipws >> l;
   if (! is.good() ) return is;
   is >> std::skipws >> b;
   if (! is.good() ) return is;
   //Make sure we are within range
   if (l < -360 || l > 360){
      is.setstate(std::ios::failbit);
      std::cerr<<"lcenter should be between -360 and 360"<<std::endl;
      return is;
   }
   if (b < -90 || b > 90){
      std::cerr<<"bcenter should be between -90 and 90"<<std::endl;
      is.setstate(std::ios::failbit);
      return is;
   }
   //Fix l to be in range 0 to 360
   if (l < 0) l += 360;
   region.fcenter = SM::Coordinate(l,b);

   switch (region.ftype) {
      case DISC:
	 is >> std::skipws >> region.fradius;
	 if (! is.good() ) return is;
	 if ( region.fradius < 0 ||  region.fradius > 180 ){
	    std::cerr<<"radius should be between 0 and 180"<<std::endl;
	    is.setstate(std::ios::failbit);
	    return is;
	 }
	 break;
      case RECTANGLE:
	 is >> std::skipws >> region.fdeltal;
	 if (! is.good() ) return is;
	 if ( region.fdeltal < 0 ||  region.fdeltal > 360 ){
	    std::cerr<<"deltal should be between 0 and 360"<<std::endl;
	    is.setstate(std::ios::failbit);
	    return is;
	 }
	 is >> std::skipws >> region.fdeltab;
	 if (! is.good() ) return is;
	 if ( region.fdeltab < 0 ||  region.fdeltab > 360 ){
	    std::cerr<<"deltab should be between 0 and 360"<<std::endl;
	    is.setstate(std::ios::failbit);
	    return is;
	 }
	 break;
   }
   return is;
}

std::ostream & operator << (std::ostream &os, const SkyRegion &region) {
   switch (region.ftype) {
      case DISC:
	 os << std::setw(15) << std::left << "DISC";
	 break;
      case RECTANGLE:
	 os << std::setw(15) << std::left << "RECTANGLE";
	 break;
   }
   os << std::setw(12) << std::setprecision(7) << std::left << region.fcenter.l();
   os << std::setw(11) << std::setprecision(7) << std::left << region.fcenter.b();
   switch (region.ftype) {
      case DISC:
	 os << std::setw(12) << std::setprecision(7) << std::left << region.fradius;
	 break;
      case RECTANGLE:
	 os << std::setw(12) << std::setprecision(7) << std::left << region.fdeltal;
	 os << std::setw(12) << std::setprecision(7) << std::left << region.fdeltab;
	 break;
   }
   return os;
}
