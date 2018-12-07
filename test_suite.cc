
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * test_suite.cc *                               galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"

#include <fort_interface.h>

int Galprop::test_suite()//AWS20050816
{
   cout<<">>>>test_suite"<<endl;
   cout<<"     Running a series of documentary tests instead of normal run"<<endl;

   test_isotope_cs();
   test_sigma_boron_dec_heinbach_simon();
   test_kinematic();
   test_He_to_H_CS();
   test_nH();
//   test_Distribution();
//   test_Particle(); //AWS20010731
   create_gcr();
   test_e_loss_compton_cc();
   test_bremss_spec_cc();
//   test_pp_meson_cc();
   test_antiproton_cc();
   test_synchrotron_cc();
   create_galaxy();
   test_source_SNR_event();
   test_cfactor();
   test_float_accuracy();

   cout<<"     End of documentary tests"<<endl;
   cout<<"<<<<test_suite"<<endl;
   return 0;
}
