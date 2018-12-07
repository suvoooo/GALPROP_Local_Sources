/*
 * Generic line of sight integration routine.
 * Provided a vector of functions (x,y,z), it returns an array (R) where the
 * output is binned into Galacto-centric rings.
 */

#ifndef _los_integration_h_
#define _los_integration_h_

#include "integ.h"
#include "constants.h"
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
using namespace std;

namespace SM {
//! The function to be integrated over the LOS
template <typename T>
class LOSfunction {
public:
   virtual ~LOSfunction() {}
    virtual T operator () ( double x, double y, double z ) const = 0;
};

//! Out of bounds error
class LOSboundError : public out_of_range {
public:
    explicit LOSboundError ( const string& arg ) :
        out_of_range ( "Out of bounds error in LOS integration: " + arg )
   {}
};

//! Invalid argument error
class LOSargumentError : public invalid_argument {
public:
   explicit LOSargumentError ( const string& arg ) :
      invalid_argument ( "Invalid argument in LOS integration: " + arg )
   {}
};

/** \brief The LOS integrator
 *
 * Assumes boundaries are in units of kpc but integrates in units of cm
 * The plane b = 0 is parallel to z = 0 and l = 0 points to x,y = (0,0) from
 * origin. Origin is limited to positioning in z and x.
 */
template <typename T>
class LOSintegrator {
    double fxmin, fxmax, fymin, fymax, fzmin, fzmax, fRmax;
    double fxorig, fyorig, fzorig;
    std::vector<double> fRbins;
    double fstepSize;
    double frelTolerance;
    double fabsTolerance;
    const unsigned int fdimensions;
    size_t fiLocalAnnulus;

    // Given a d,l,b returns a x,y,z
    void coordinates(double d, double cosl, double sinl, double cosb, double sinb, double &x, double &y, double &z);

    void findLocalAnnulus();
    size_t findAnnulus(double R) const;
    void checkBins() const;
    void output(T v) const;

    // TODO fix when variable step size integration is implemented.
    class simFun : public integ::vecFun<T> {
        //Create the function to integrate, needs to know the geometry
        const LOSintegrator &lint;
    };

public:
    LOSintegrator ( double Rmax, double zmin, double zmax, std::vector<double> Rbins, double Rorig, double zorig, double stepSize ) :
            fRmax ( Rmax ),
            fzmin ( zmin ),
            fzmax ( zmax ),
            fRbins ( Rbins ),
            fxorig ( Rorig ),
            fzorig ( zorig ),
            fstepSize ( stepSize ),
            frelTolerance ( 1e-3 ),
            fabsTolerance ( 0 ),
            fdimensions ( 2 ) {
        if ( fzmin >= 0 )
            throw ( LOSboundError ( "zmin must be < 0" ) );
        if ( fzmax <= 0 )
            throw ( LOSboundError ( "zmax must be > 0" ) );
        if ( fRmax <= 0 )
            throw ( LOSboundError ( "Rmax must be > 0" ) );
        if ( fzorig < fzmin || fzorig > fzmax )
            throw ( LOSboundError ( "we must have zmin <= zorig <= zmax" ) );
        if ( fxorig < 0 || fxorig > fRmax )
            throw ( LOSboundError ( "we must have 0 <= Rorig <= Rmax" ) );
	checkBins();
	findLocalAnnulus();
    }

    LOSintegrator ( double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, std::vector<double> Rbins, double xorig, double zorig, double stepSize ) :
            fxmin ( xmin ),
            fxmax ( xmax ),
            fymin ( ymin ),
            fymax ( ymax ),
            fzmin ( zmin ),
            fzmax ( zmax ),
            fzorig ( zorig ),
            fxorig ( xorig ),
            fyorig ( 0 ),
            fRbins ( Rbins ),
            fstepSize ( stepSize ),
            frelTolerance ( 1e-3 ),
            fabsTolerance ( 0 ),
            fdimensions ( 3 ) {
        if ( fzmin >= 0 )
            throw ( LOSboundError ( "zmin must be < 0" ) );
        if ( fzmax <= 0 )
            throw ( LOSboundError ( "zmax must be > 0" ) );
        if ( fxmin >= 0 )
            throw ( LOSboundError ( "xmin must be < 0" ) );
        if ( fxmax <= 0 )
            throw ( LOSboundError ( "xmax must be > 0" ) );
        if ( fymin >= 0 )
            throw ( LOSboundError ( "ymin must be < 0" ) );
        if ( fymax <= 0 )
            throw ( LOSboundError ( "ymax must be > 0" ) );
        if ( fzorig < fzmin || fzorig > fzmax )
            throw ( LOSboundError ( "we must have zmin <= zorig <= zmax" ) );
        if ( fxorig < fxmin || fxorig > fxmax )
            throw ( LOSboundError ( "we must have xmin <= xorig <= xmax" ) );
	checkBins();
	findLocalAnnulus();
    }

    std::vector<std::vector<T> > integrate ( const double l, const double b, const double dl, std::vector<LOSfunction<T>*> );
};

template <typename T>
void LOSintegrator<T>::output(T t) const {
   std::cout<<t<<std::endl;
}

template <> void SM::LOSintegrator<std::valarray<double> >::output (std::valarray<double> t) const;

template <typename T>
void LOSintegrator<T>::coordinates( double d, double cosl, double sinl, double cosb, double sinb, double &x, double &y, double &z ) {
   z = d*sinb + fzorig;
   if ( fxorig < 0 ) {
      x = fxorig + d*cosb*cosl;
      y = d*cosb*sinl;
   } else {
      x = fxorig - d*cosb*cosl;
      y = -d*cosb*sinl;
   }
}

template <typename T>
void LOSintegrator<T>::findLocalAnnulus() {
   fiLocalAnnulus = findAnnulus(abs(fxorig));
}
template <typename T>
size_t LOSintegrator<T>::findAnnulus(double R) const {
   size_t i(0);
   while (fRbins[i] < R ) 
      ++i;
   if ( i != 0 ) --i; //Since we want to return the lower boundary
   return i;
}

template <typename T>
void LOSintegrator<T>::checkBins() const {
   //Make sure they are monotonically increasing
   for ( size_t i = 0; i < fRbins.size()-1; ++i ) {
      if ( fRbins[i] >= fRbins[i+1] )
	 throw(LOSargumentError("Rbins must be monotonically increasing"));
   }
   //Make sure the outermost ring contains the entire integration region
   double Rmax;
   if (fdimensions == 2) {
      Rmax = fRmax;
   } else if ( fdimensions == 3 ) {
      double xmax = fxmax > abs(fxmin) ? fxmax : abs(fxmin);
      double ymax = fymax > abs(fymin) ? fymax : abs(fymin);
      Rmax = sqrt(xmax*xmax+ymax*ymax);
   }
   if ( fRbins[fRbins.size()-1] < Rmax ) {
      throw(LOSboundError("Integration region must be entirely bound by Rbins"));
   }
}



template <typename T>
std::vector< std::vector<T> > LOSintegrator<T>::integrate ( const double l, const double b, const double dl, std::vector<LOSfunction<T>*> funcs) {
    //Use a proper integrator with variable step size
    //Begin by finding the limits on the distance along the line of sight
    double dmin ( 0 ), dmax;
    const double dtr=acos ( -1. ) /180.;                           // conversion degrees to radians
    const double sinb=sin ( b*dtr );
    const double cosb=cos ( b*dtr );
    const double sinl=sin ( l*dtr );
    const double cosl=cos ( l*dtr );

    double zdmax ( 0 );
    if ( b != 0 ){
        zdmax = b < 0 ? ( fzmin-fzorig ) /sinb : ( fzmax - fzorig ) /sinb;
    } else {
       // make it 10 times the biggest possible dimension
       if ( fdimensions == 2 ) {
	  zdmax = 10*sqrt(4*fRmax*fRmax + (fzmax-fzmin)*(fzmax-fzmin));
       } else if ( fdimensions == 3 ) {
	  zdmax = 10*sqrt( (fxmax-fxmin)*(fxmax-fxmin) + (fymax-fymin)*(fymax-fymin) + (fzmax-fzmin)*(fzmax-fzmin) );
       }
    }

    if ( fdimensions == 2 ) {
        if ( cosb == 0 ) {
            dmax = zdmax;
        } else {
            double rdmax = ( fxorig*cosb*cosl + sqrt ( fxorig*fxorig*cosb*cosb*cosl + ( fRmax*fRmax - fxorig*fxorig ) *cosb*cosb ) ) / ( cosb*cosb );
            dmax = zdmax < rdmax ? zdmax : rdmax;
        }
    } else if ( fdimensions == 3 ) {
        if ( cosb == 0 ) {
            dmax = zdmax;
        } else {
            if ( fxorig < 0 ) {
                if ( l == 0 ) {
                    double xdmax = ( fxmax - fxorig ) /cosb;
                    dmax = xdmax < zdmax ? xdmax : zdmax;
                } else if ( l == 180 ) {
                    double xdmax = ( fxorig - fxmin ) /cosb;
                    dmax = xdmax < zdmax ? xdmax : zdmax;
                } else if ( l == 90 ) {
                    double ydmax = fymax/cosb;
                    dmax = ydmax < zdmax ? ydmax : zdmax;
                } else if ( l == 270 ) {
                    double ydmax = fymin/cosb;
                    dmax = ydmax < zdmax ? ydmax : zdmax;
                } else {
                    double xdmax = ( l > 270 || l < 90 ) ? ( fxmax-fxorig ) / ( cosb*cosl ) : ( fxmin-fxorig ) / ( cosb*cosl );
                    double ydmax = ( l < 180 ) ? fymax/ ( cosb*sinl ) : fymin/ ( cosb*sinl );
                    double max = xdmax < ydmax ? xdmax : ydmax;
                    dmax = max < zdmax ? max : zdmax;
                }
            } else {
                if ( l == 0 ) {
                    double xdmax = ( fxorig-fxmin ) /cosb;
                    dmax = xdmax < zdmax ? xdmax : zdmax;
                } else if ( l == 180 ) {
                    double xdmax = ( fxmax-fxorig ) /cosb;
                    dmax = xdmax < zdmax ? xdmax : zdmax;
                } else if ( l == 90 ) {
                    double ydmax = fymin/cosb;
                    dmax = ydmax < zdmax ? ydmax : zdmax;
		} else if ( l == 270 ) {
	    	   double ydmax = fymax/cosb;
       		   dmax = ydmax < zdmax ? ydmax : zdmax;
   		} else {
		   double xdmax = ( l > 270 || l < 90 ) ? ( fxorig-fxmin ) / ( cosb*cosl ) : ( fxorig-fxmax ) / ( cosb*cosl );
		   double ydmax = ( l < 180 ) ? -fymin/ ( cosb*sinl ) : -fymax/ ( cosb*sinl );
		   double max = xdmax < ydmax ? xdmax : ydmax;
		   dmax = max < zdmax ? max : zdmax;
		}
	    }
	}

    }


    // Integrate between dmin and dmax, but split up into rings
    // Should possibly make ring integration an option.

    // TODO Find boundaries for the rings (dual in the inner galaxy) and
    // integrate
    //
    // For now we just use a fixed step size and trapezoid rule for
    // integration.
    double d(0);
    double xx, yy, zz;

    // Calculate the first point and initialize the output vector.  As we
    // expect valarrays we have to initialize on copy constructor rather than
    // assignment.  This also initializes the integral in the local ring.
    coordinates(d,cosl,sinl,cosb,sinb,xx,yy,zz);
    std::vector< std::vector<T> > integ(funcs.size());
    for ( size_t i = 0; i < funcs.size(); ++i ) {
       T tmp((*funcs[i])(xx,yy,zz));
       integ[i].resize(fRbins.size()-1, tmp);
       integ[i][fiLocalAnnulus] *= fstepSize;
       for (size_t j = 0; j < fiLocalAnnulus; ++j)
	  integ[i][j] = 0;
       for (size_t j = fiLocalAnnulus+1; j < integ[i].size(); ++j)
	  integ[i][j] = 0;
    }

    size_t prevAnn(fiLocalAnnulus), currAnn;

    // Boolean variable to keep track of split pixels
    bool noSplit = true;

    // Now we step along the line of sight
    while ( d + fstepSize < dmax ) {
       d += fstepSize;
       coordinates(d,cosl,sinl,cosb,sinb,xx,yy,zz);
       double RR = sqrt(xx*xx+yy*yy);
       // TODO make the step size such that we end on a boundary
       // Now we readjust the boundaries to the step size.  Should be fine
       // for now.
       currAnn = findAnnulus(RR);

       // Perform check for a ring boundary if dl > 0.
       // Note that this splits the integration value between the annuli as if the pixels
       // where square shaped
       if ( dl > 0 && ( l < 90 || l > 270 ) ) {
	  double lmin = l-dl*0.49;
	  double lmax = l+dl*0.49;
	  double xmin, ymin, zmin, xmax, ymax, zmax;
	  coordinates(d, cos(lmin*dtr), sin(lmin*dtr), cosb, sinb, xmin, ymin, zmin);
	  coordinates(d, cos(lmax*dtr), sin(lmax*dtr), cosb, sinb, xmax, ymax, zmax);
	  double Rmin = sqrt(xmin*xmin+ymin*ymin);
	  double Rmax = sqrt(xmax*xmax+ymax*ymax);

      	  size_t minAnn, maxAnn;
	  minAnn = findAnnulus(Rmin);
	  maxAnn = findAnnulus(Rmax);

	  // Find the edge of the ring, to make sure we only apply this
	  // correction for tangential split of pixels
	  double ledge(0);
	  if ( l > 270. )
	     ledge=360 - asin ( fRbins[maxAnn+1]/abs(fxorig) ) / dtr;
	  else
    	     ledge=asin ( fRbins[minAnn+1]/abs(fxorig) ) / dtr;

	  //std::cout<<"Edge test: "<<ledge<<", "<<lmin<<", "<<lmax<<std::endl;
	  //std::cout<<fRbins[maxAnn+1]<<", "<<fxorig<<", "<<dtr<<std::endl;

      	  if ( lmin < ledge && ledge < lmax && minAnn != maxAnn ) {
	     if ( currAnn != minAnn ) {
		//We apply the weight so that when all annuli are summed
		//up, it would give very similar results as if no splitting
		//between annuli was performed.  This is of course not
		//entirely true as we change the integration path in the
		//split pixel.
		double weight = (ledge-lmin)/(lmax-lmin);
		if ( noSplit ) { // Last step we had no split
		   for ( size_t i = 0; i < funcs.size(); ++i ) {
		      T tmp((*funcs[i])(xx,yy,zz));
		      integ[i][minAnn] += (*funcs[i])(xmin,ymin,zmin)*fstepSize*weight;
		      integ[i][currAnn] += tmp*fstepSize*(1-weight);
		      integ[i][prevAnn] += tmp*fstepSize;
		   }
		} else {
		   for ( size_t i = 0; i < funcs.size(); ++i ) {
		      T tmp((*funcs[i])(xx,yy,zz));
		      integ[i][minAnn] += (*funcs[i])(xmin,ymin,zmin)*fstepSize*weight*2;
		      integ[i][currAnn] += tmp*fstepSize*(1-weight);
		      integ[i][prevAnn] += tmp*fstepSize*(1-weight);
		   }
		}
	     } else {
		double weight = (lmax-ledge)/(lmax-lmin);
		if ( noSplit ) { // Last step we had no split
		   for ( size_t i = 0; i < funcs.size(); ++i ) {
		      T tmp((*funcs[i])(xx,yy,zz));
		      integ[i][maxAnn] += (*funcs[i])(xmax,ymax,zmax)*fstepSize*weight;
		      integ[i][currAnn] += tmp*fstepSize*(1-weight);
		      integ[i][prevAnn] += tmp*fstepSize;
		   }
		} else {
		   for ( size_t i = 0; i < funcs.size(); ++i ) {
		      T tmp((*funcs[i])(xx,yy,zz));
		      integ[i][maxAnn] += (*funcs[i])(xmax,ymax,zmax)*fstepSize*weight*2;
		      integ[i][currAnn] += tmp*fstepSize*(1-weight);
		      integ[i][prevAnn] += tmp*fstepSize*(1-weight);
		   }
		}
	     }
	     noSplit = false;
	  } else {
	     //There should be a final addition to the split ring, but it
	     //is only a minor error.
	     noSplit = true;
	  }

       } //End split pixel test

       if ( noSplit ) {
	  // Perform the unsplit integration.  Since we are using the trapezoid
	  // rule, we add the value to both previous and current annuli.
	  // If the are the same, we get twice the value as we are supposed
	  // to.  Otherwise, we both end and start the integral.
	  for ( size_t i = 0; i < funcs.size(); ++i ) {
	     T tmp((*funcs[i])(xx,yy,zz));
	     integ[i][currAnn] += tmp*fstepSize;
	     integ[i][prevAnn] += tmp*fstepSize;
	  }
       }

       prevAnn = currAnn;
    }

    //Now we must add the last step, assuming it is not split, which is
    //safe in most cases
    coordinates(dmax,cosl,sinl,cosb,sinb,xx,yy,zz);
    double RR = sqrt(xx*xx+yy*yy);
    currAnn = findAnnulus(RR);
    for ( size_t i = 0; i < funcs.size(); ++i ) 
       integ[i][currAnn] +=(*funcs[i])(xx,yy,zz)*(dmax-d);

    //We must also multiply with kpc2cm and divide by 2 to get the correct
    //answer
    for ( size_t i = 0; i < funcs.size(); ++i )
       for ( size_t j = 0; j < integ[i].size(); ++j )
   	  integ[i][j] *= kpc2cm/2;

    return integ;
}

} //End namespace
#endif
