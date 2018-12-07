#ifndef fort_interface1_h
#define fort_interface1_h

#include <string>

void set_sigma_cc(const std::string& path);
double e_loss_compton_cc(double,double);
double bremss_spec_cc(double,double,int,int);
double pp_meson_cc(double,double,int,int,int);
int test_e_loss_compton_cc();
int test_bremss_spec_cc();
int test_pp_meson_cc();

#endif
