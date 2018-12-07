#ifndef fort_interface2_h
#define fort_interface2_h

#include <string>

double wsigma_cc(int,int,int,int,double);                 // IMOS20020502
double yieldx_cc(int,int,int,int,float);                  // IMOS20020502
void sigtap_cc(int, const std::string& path);                                      // IMOS20010511
double sighad_cc(int,double,double,double,double,double); // IMOS20020502
double cfactor_cc(int,double,double,double,double,double,double,double,double);
double antiproton_cc(int,double,double,int,int,int,int);  // IMOS20010511
double synchrotron_cc(double,double,double); 
int test_antiproton_cc();
int test_synchrotron_cc();
int test_cfactor();

#endif

