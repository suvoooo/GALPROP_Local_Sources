#ifndef _integ_h_
#define _integ_h_

#include<vector>
#include<valarray>

namespace integ {
template <typename T>
class vecFun {
public:
   virtual ~vecFun() {}
    virtual std::vector<T> operator () ( double x ) const = 0;
};

template <typename T>
double simvec ( double lbound, double ubound, double step, double relEps, double absEps, const vecFun<T> &fu );

inline double min ( double value ) { return value; }

inline double max ( const std::valarray<double> &value ) { return value.max(); }

template <typename T>
double simvec ( double lbound, double ubound, double step, double relEps, double absEps, const vecFun<T> &fu ) {
    step = ubound >= lbound ? abs ( step ) : -abs ( step );
    double sign = step >= 0 ? 1. : -1.;
    double x ( lbound );
    double C ( 0. );

    //Here is a little hack to initialize T in case it is a valarray
    std::vector<T> *f0 = new std::vector<T> ( fu ( x ) );
    std::vector<T> integ ( *f0 ), absInteg ( *f0 );
    for ( size_t i = 0; i < f0->size(); ++i ) {
        integ[i] = 0;
        absInteg[i] = 0;
    }

    double P[5];
    P[1]=P[3]=4./3.;
    P[2]=2./3.;
    P[4]=1./3.;

    std::vector< std::vector<T>* > fvalue ( 8,0 );
    fvalue[0] = f0;

    if ( ubound==lbound )
        return integ;

    relEps = abs ( relEps );
    absEps = abs ( absEps );
    while ( C == 0 ) {
        double x0 = x;
        if ( ( x0 + 4.*step - ubound ) *sign > 0 ) {
            step = ( ubound - x0 ) /4.;
            if ( step == 0 )
                return integ;
            for ( int k = 1; k < 7; ++k ) {
                delete fvalue[k] = 0;
                fvalue[k] = 0;
            }
            C=1.;
        }
        std::vector<T> DI2 ( *fvalue[0] );
        std::vector<T> DI3 ( *fvalue[0] );
        for ( size_t i = 0; i < DI3.size(); ++i )
            DI3[i] = abs ( DI3[i] );

        for ( int k = 1; k < 5; ++k ) {
            x += step;
            if ( ( x-ubound ) *sign >= 0 )
                x = ubound;
            if ( fvalue[k] == 0 )
                fvalue[k] = new std::vector<T> ( fu ( x ) );
            for ( size_t i = 0; i < DI2.size; ++i ) {
                DI2[i] = DI2[i] + P[k]* ( *fvalue[k] ) [i]*step;
                DI3[i] = DI3[i] + P[k]*abs ( ( *fvalue[k] ) [i] ) *step;
            }
        }
        std::vector<T> DI1 ( *fvalue[0] );
        for ( size_t i = 0; i < DI1.size; ++i )
            DI1[i] = ( ( *fvalue[0] ) [i] + 4.* ( *fvalue[2] ) [i] + ( *fvalue[4] ) [i] ) *2./3.*step;
        T tmp ( DI2[0]-DI1[0] );
        double DELTA = max ( tmp );
        for ( size_t i = 1; i < DI1.size; ++i ) {
            tmp = DI2[i]-DI1[i];
            double d = max ( tmp );
            DELTA = d > DELTA ? d : DELTA;
        }
        double EPS;
        if ( relEps != 0 || absEps != 0 ) {
            tmp = absInteg[0]+DI3[0];
            double maxVal = max ( tmp );
            for ( size_t i = 1; i < DI1.size; ++i ) {
                tmp = absInteg[i]+DI3[i];
                double d = max ( tmp );
                maxVal = d > maxVal ? d : maxVal;
            }

            EPS = maxVal*relEps;
            if ( EPS < absEps )
                EPS = absEps;
        } else {
            EPS = DELTA*2.;
        }
        if ( DELTA < EPS ) {
            if ( DELTA < EPS/8. ) {
                step = step*2.;
                fvalue[0] = fvalue[4];
                fvalue[1] = fvalue[5];
                fvalue[2] = fvalue[6];
                for ( int k = 3; k < 7; ++k ) {
                    delete fvalue[k];
                    fvalue[k] = 0;
                }
            } else {
                fvalue[0] = fvalue[4];
                fvalue[2] = fvalue[5];
                fvalue[4] = fvalue[6];
                delete fvalue[1];
                fvalue[1] = 0;
                delete fvalue[3];
                fvalue[3] = 0;
                delete fvalue[5];
                fvalue[5] = 0;
	       	delete fvalue[6];
                fvalue[6] = 0;
   	    }
            for ( size_t i = 0; i < DI1.size(); ++i ) {
                integ[i] += DI2[i] + ( DI2[i]-DI1[i] ) /15.;
                absInteg[i] += DI3[i];
            }
        } else {
            step = step/2.0;
            fvalue[6] = fvalue[4];
            fvalue[5] = fvalue[3];
            fvalue[4] = fvalue[2];
            fvalue[2] = fvalue[1];
            delete fvalue[1];
            fvalue[1] = 0;
            delete fvalue[3];
            fvalue[3] = 0;
            x = x0;
            C = 0;
        }
    }
    return integ;
}
}

#endif
