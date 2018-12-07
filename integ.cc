//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * sim.cc                                         galprop package * 3/21/2006
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

/***********************************************************************
c ### The routine has been rewritten from FORTRAN source code ###
c calculation the definite integral by Simpson's method with the automatic
c choice of the integration step
C INPUT: A1,B1 - the limits of integration; H1 - the initial step;
C REPS,AEPS - the relative and absolute precision; FU - the name of the
C user-defined function f(x); OUTPUT: sim - the value of the integral;
c other values that are calculated in parallel:
C AIH - the value of integral with one more step of integration;
C AIABS - the value of the integral for module of the integrand;
C # NOTE # the subprogram returns the value of integral as one of the
C precise conditions (AEPS1,EPS1) are reached; when AEPS1=EPS1=0,
c then it is calculated with the constant step H1. See appended test case.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
***********************************************************************/

using namespace std;
#include <cmath>
#include <stdio.h>
#include <iostream>
#define sign(a,b) (((b) > 0.) ? (a) : (-a))

double simnew ( double A, double B, double H, double REPS, double AEPS, double ( *fu ) ( double ) ) {
    H = B >= A ? abs ( H ) : -abs ( H );
    double S = H >= 0 ? 1. : -1.;
    double AI ( 0 ), AIH ( 0 ), AIABS ( 0 );
    double P[5], F[7];
    P[1]=P[3]=4.;
    P[2]=2.;
    P[4]=1.;
    if ( B==A )
        return AI;
    REPS = abs ( REPS );
    AEPS = abs ( AEPS );
    for ( int k=0; k < 7; ++k )
        F[k] = 1e20;
    double X ( A );
    double C ( 0. );
    F[0] = fu ( X ) /3.;
    while ( C == 0 ) {
        double X0 = X;
        if ( ( X0 + 4.*H - B ) *S > 0 ) {
            H = ( B - X0 ) /4.;
            if ( H == 0 )
                return AI;
            for ( int k = 1; k < 7; ++k )
                F[k] = 1e20;
            C=1.;
        }
        double DI2 = F[0];
        double DI3 = abs ( F[0] );
        for ( int k = 1; k < 5; ++k ) {
            X += H;
            if ( ( X-B ) *S >= 0 )
                X = B;
            if ( F[k] == 1e20 )
                F[k] = fu ( X ) /3.;
            DI2 = DI2 + P[k]*F[k];
            DI3 = DI3 + P[k]*abs ( F[k] );
        }
        double DI1 = ( F[0] + 4.*F[2] + F[4] ) *2.*H;
        DI2 = DI2*H;
        DI3 = DI3*H;
        double DELTA = abs ( DI2-DI1 );
        double EPS;
        if ( REPS != 0 || AEPS != 0 ) {
            EPS = abs ( ( AIABS+DI3 ) *REPS );
            if ( EPS < AEPS )
                EPS = AEPS;
        } else {
            EPS = DELTA*2.;
        }
        if ( DELTA < EPS ) {
            if ( DELTA < EPS/8. ) {
                H = H*2.;
                F[0] = F[4];
                F[1] = F[5];
                F[2] = F[6];
                for ( int k = 3; k < 7; ++k )
                    F[k] = 1e20;
            } else {
                F[0] = F[4];
                F[2] = F[5];
                F[4] = F[6];
                F[1] = 1e20;
                F[3] = 1e20;
                F[5] = 1e20;
                F[6] = 1e20;
            }
            DI1 = DI2 + ( DI2-DI1 ) /15.;
            AI = AI + DI1;
            AIH = AIH + DI2;
            AIABS = AIABS + DI3;
        } else {
            H = H/2.0;
            F[6] = F[4];
            F[5] = F[3];
            F[4] = F[2];
            F[2] = F[1];
            F[1] = 1e20;
            F[3] = 1e20;
            X = X0;
            C = 0;
        }
    }
    return AI;
}


//double sim(double, double, double, double, double, double(*)(double));

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
double sim ( double A, double B, double H, double REPS, double AEPS, double ( *fu ) ( double ) ) {
    int K;
    double F[8],P[6],S,C,X,X0,AI,AIH,AIABS,DI1,DI2,DI3,EPS,DELTA;
    H=sign ( H,B-A );
    S=sign ( 1.,H );
    AI=AIH=AIABS=0.;
    P[2]=P[4]=4.;
    P[3]=2.;
    P[5]=1.;
    if ( B-A==0. ) return ( AI );
    REPS=abs ( REPS );
    AEPS=abs ( AEPS );
    for ( K=1;K<8;F[K++]=1.e20 );
    X=A;
    C=0.;
    F[1]=fu ( X ) /3.;
L4:
    X0=X;
    if ( ( X0+4.*H-B ) *S>0. ) {
        H= ( B-X0 ) /4.;
        if ( H==0. ) return ( AI );
        for ( K=2;K<8;F[K++]=1.e20 );
        C=1.;
    }
L5:
    DI2=F[1];
    DI3=abs ( F[1] );
    for ( K=2;K<6;K++ ) {
        X+=H;
        if ( ( X-B ) *S>=0. ) X=B;
        if ( F[K]-1.e20==0. ) F[K]=fu ( X ) /3.;
        DI2+=P[K]*F[K];
        DI3+=P[K]*abs ( F[K] );
    }
    DI1= ( F[1]+4.*F[3]+F[5] ) *2.*H;
    DI2*=H;
    DI3*=H;
    if ( REPS==0.&& AEPS==0. ) goto L14;
    EPS=abs ( ( AIABS+DI3 ) *REPS );
    if ( EPS-AEPS<0 ) EPS=AEPS;
    DELTA=abs ( DI2-DI1 );
    if ( DELTA-EPS<0. ) { if ( DELTA-EPS/8.>=0. ) goto L14; }
    else goto L21;
    H*=2.;
    F[1]=F[5];
    F[2]=F[6];
    F[3]=F[7];
    for ( K=4;K<8;F[K++]=1.e20 );
    goto L18;
L14:
    F[1]=F[5];
    F[3]=F[6];
    F[5]=F[7];
    F[2]=F[4]=F[6]=F[7]=1.e20;
L18:
    DI1=DI2+ ( DI2-DI1 ) /15.;
    AI+=DI1;
    AIH+=DI2;
    AIABS+=DI3;
    goto L22;
L21:
    H/=2.;
    F[7]=F[5];
    F[6]=F[4];
    F[5]=F[3];
    F[3]=F[2];
    F[2]=F[4]=1.e20;
    X=X0;
    C=0.;
    goto L5;
L22:
    if ( C==0 ) goto L4;
    return ( AI );
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// test area
/*
double sim(double, double, double, double, double, double(*)(double));
double fu(double); //test function -appended
double fu(double x)
{
//  return(sin(x));
  return(pow(x,-10));
}
int main()
{
  double a=1., b=10* 3.14159265358979/2.;
//  cout<<sim(a,b,0.0001,1.e-5,1.e-20,&fu)<<" = "<<-(cos(b)-cos(a))<<endl;
  cout<<sim(a,b,0.0001,1.e-5,1.e-20,&fu)<<" = "<<(-pow(b,-9)+pow(a,-9))/9.<<endl;
}
*/


/*
c Fortran code
      SUBROUTINE SIM(A1,B1,H1,REPS1,AEPS1,FU,AI)
c***********************************************************************
c calculation the definite integral by Simpson's method with the automatic
c choice of the step of integration
C INPUT: A1,B1 - the limits of integration; H1 - the initial step;
C REPS1,AEPS1 - the relative and absolute precision; FU - the name of the
C user function f(x); OUTPUT: AI - the value of the integral;
C AIH - the value of integral with one more step of integration;
C AIABS - the value of the integral for module of the integrand;
C # NOTES # the subprogram returns the value of integral as one of the
C precise conditions (AEPS1,EPS1) are reached; when AEPS1=EPS1=0,
c then it is calculated with the constant step H1.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      IMPLICIT real*8 (A-H,O-Z), integer (I-N)
      DIMENSION F(7),P(5)
      H=dSIGN(H1,B1-A1)
      S=dSIGN(1.d0,H)
      A=A1
      B=B1
      AI=0.d0
      AIH=0.d0
      AIABS=0.d0
      P(2)=4.d0
      P(4)=4.d0
      P(3)=2.d0
      P(5)=1.d0
      IF(B-A) 1,2,1
    1 REPS=dABS(REPS1)
      AEPS=dABS(AEPS1)
      DO 3 K=1,7
    3    F(K)=1.d20
      X=A
      C=0.d0
      F(1)=FU(X)/3.d0
    4 X0=X
      IF((X0+4.d0*H-B)*S) 5,5,6
    6 H=(B-X0)/4.d0
      IF(H) 7,2,7
    7 DO 8 K=2,7
    8    F(K)=1.d20
      C=1.d0
    5 DI2=F(1)
      DI3=dABS(F(1))
      DO 9 K=2,5
         X=X+H
         IF((X-B)*S) 23,24,24
   24    X=B
   23    IF(F(K)-1.d20) 10,11,10
   11    F(K)=FU(X)/3.d0
   10    DI2=DI2+P(K)*F(K)
    9    DI3=DI3+P(K)*dABS(F(K))
      DI1=(F(1)+4.d0*F(3)+F(5))*2.d0*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=dABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=dABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.d0) 17,14,14
   17 H=H*2.d0
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
   19    F(K)=1.d20
      GOTO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=1.d20
      F(4)=1.d20
      F(6)=1.d20
      F(7)=1.d20
   18 DI1=DI2+(DI2-DI1)/15.d0
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GOTO 22
   21 H=H/2.d0
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=1.d20
      F(4)=1.d20
      X=X0
      C=0.d0
      GOTO 5
   22 IF(C) 2,4,2
    2 RETURN
      END
*/

