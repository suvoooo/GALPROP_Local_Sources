#include "Coordinate.h"
#include "pointing.h"
#include "PhysicalConstants.h"
#include <iostream>

SM::Coordinate::Coordinate(const double l, const double b) : m_l(l), m_b(b){
   //Normalize the longitude range to 0,360
   while (m_l < 0){
      m_l += 360;
   }
   while (m_l >= 360){
      m_l -= 360;
   }
   //Clip the latitude range at +-90
   if (fabs(m_b) > 90)
      m_b = (m_b<0) ? -90 : 90;
}


SM::Coordinate::Coordinate(const pointing &point) {
   m_b = 90 - point.theta/utl::kConvertDegreesToRadians;
   m_l = point.phi/utl::kConvertDegreesToRadians;
}

pointing SM::Coordinate::healpixAng() const{
   const double theta = utl::kPi/2 - utl::kConvertDegreesToRadians*m_b;
   const double phi = utl::kConvertDegreesToRadians*m_l;
   return pointing(theta,phi);
}

namespace SM {
std::ostream & operator << (std::ostream &os, const SM::Coordinate &coord) {
   os << "(" << coord.l() <<","<<coord.b()<<")";
   return os;
}
}
