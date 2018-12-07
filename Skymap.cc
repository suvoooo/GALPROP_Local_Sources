#include "Skymap.h"
#ifdef HAVE_ASTRO
extern "C"{
#include <astro.h>
}
//Do to stupid convention in astro, we must create wrapper functions
void eq_gal2(double mj, double DEC, double RA, double *b, double *l){
   eq_gal(mj,RA,DEC,b,l);
}
void eq_ecl2(double mj, double DEC, double RA, double *lt, double *lg){
   eq_ecl(mj,RA,DEC,lt,lg);
}
void gal_eq2(double mj, double b, double l, double *DEC, double *RA){
   gal_eq(mj, b, l, RA, DEC);
}
void ecl_eq2(double mj, double lt, double lg, double *DEC, double *RA){
   ecl_eq(mj, lt, lg, RA, DEC);
}
//Functions to convert from gal to ecl and wise versa
void gal_ecl2(double mj, double b, double l, double *lt, double *lg){
   double RA, DEC;
   gal_eq(mj, b, l, &RA, &DEC);
   eq_ecl(mj, RA, DEC, lt, lg);
}
void ecl_gal2(double mj, double lt, double lg, double *b, double *l){
   double RA, DEC;
   ecl_eq(mj, lt, lg, &RA, &DEC);
   eq_gal(mj, RA, DEC, b, l);
}
void empty(double mj, double x, double y, double *z, double *w){
}
#endif

