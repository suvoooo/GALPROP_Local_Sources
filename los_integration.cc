#include "los_integration.h"
#include <valarray>
#include <iostream>

template <> void SM::LOSintegrator<std::valarray<double> >::output (std::valarray<double> t) const{
   for (int i = 0; i < t.size(); ++i) {
      std::cout<<t[i]<<", "<<std::endl;
   }
}

