
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * fort_interface1.cc *                          galprop package * 2001/05/11 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// C++ wrapper for fortran routines 2001/05/11
// These routines are to be renamed: "_(" should be replaced with "__(" 
// and vise versa when changing platform (e.g. SUN UNIX, Linux).

using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"

#include <config.h>

#define BREMSS_SPEC_F77 F77_FUNC_(bremss_spec,BREMSS_SPEC)
#ifdef __cplusplus
extern "C" 
#endif
void BREMSS_SPEC_F77(double*, double*, int*, int*, double*);
//bremss_spec_(double*,double*,int*,int*,double*);

#define SET_SIGMA_F77 F77_FUNC_(set_sigma,SET_SIGMA)
#ifdef __cplusplus
extern "C" 
#endif
void SET_SIGMA_F77(int*, char*, int*);//set_sigma_(int*);                 // IMOS20020502

#define E_LOSS_COMPTON_F77 F77_FUNC_(e_loss_compton,E_LOSS_COMPTON)
#ifdef __cplusplus
extern "C" 
#endif
double E_LOSS_COMPTON_F77(double*, double*);//e_loss_compton_(double*,double*);

#define PP_MESON_F77 F77_FUNC_(pp_meson, PP_MESON)
#ifdef __cplusplus
extern "C" 
#endif
double PP_MESON_F77(double*, double*, int*, int*, int*);//pp_meson_(double*,double*,int*,int*,int*);

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

// initialization of Webber's code
void set_sigma_cc(const string& path)
{ 
   int  cdr=99;

   const string fullFilename = path + "/" + "WNEWTR_082693.CDR.dat";

   char* fn = const_cast<char*>(fullFilename.c_str());

   int length = fullFilename.size();

   SET_SIGMA_F77(&cdr, fn, &length);                             // IMOS20020502

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double bremss_spec_cc(double Egam, double E0, int IZ1, int Ne1)
{ 
   double dSdK;
   BREMSS_SPEC_F77(&Egam, &E0, &IZ1, &Ne1, &dSdK);
//   cout<< "bremss_spec_cc:Egam= "<<Egam<<" E0="<<E0<<" dSdK="<<dSdK<<endl;
   return dSdK;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int test_bremss_spec_cc()
{
   cout<<">>>>test_bremss_spec_cc"<<endl;
 
   int IZ1=1, Ne1=0;
   for(double Egam=10; Egam<1.e4; Egam*=10) // TOTAL electron energy in MeV
   {
      cout<<endl;
      for(double E0=10; E0<1.e7; E0*=2)
         cout<<"Egam="<<Egam<<" E0="<<E0<<" IZ1="<<IZ1<<" Ne1="<<Ne1
            <<" bremss_spec_cc="<<bremss_spec_cc(Egam,E0,IZ1,Ne1)<<endl;
   }
   Ne1=1;
   for(double Egam=10; Egam<1.e4; Egam*=10)// TOTAL electron energy in MeV
   {
      cout<<endl;
      for(double E0=10; E0<1.e7; E0*=2)
         cout<<"Egam="<<Egam<<" E0="<<E0<<" IZ1="<<IZ1<<" Ne1="<<Ne1
            <<" bremss_spec_cc="<<bremss_spec_cc(Egam,E0,IZ1,Ne1)<<endl;
   }
   cout<<"<<<<test_bremss_spec_cc"<<endl;
   return 0;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

// routine is called by FORTRAN routine antiproton (file antiproton.f)
// hence underscore is appended and extern "C" supplied           // IMOS20010511

#define NUCLEON_CS_F77 F77_FUNC_(nucleon_cs,NUCLEON_CS)
#ifdef __cplusplus
extern "C" 
#endif
void NUCLEON_CS_F77(int *option, double *Ek1, int *Zp1, int *Zt1, int *At1,
		    double *PP_inel, double *PA_inel, double *aPP_non, 
		    double *aPA_non, double *aPP_ann, double *aPA_ann) {

  double Ek = *Ek1;
  int Zp = *Zp1, Zt = *Zt1, At = *At1,opt = *option;
  nucleon_cs(opt,Ek,Zp,Zt,At,PP_inel,PA_inel,aPP_non,aPA_non,aPP_ann,aPA_ann);

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
/*  //IMOS20060411 transferred to galprop.cc file
// routine is called by FORTRAN routine emiss(r) (in file cfactor.f)
// hence underscore is appended and extern "C" supplied

extern "C" float isrf_energy_density_(float *rr, float *zz)
{
//   cout<<"isrf_energy_density_  "<<*rr<<" "<<*zz<<endl;
   float r=*rr;
   float z=*zz;
   //return ( isrf_energy_density(r,z) );AWS20050817 problem to use in Galprop class
   cout<< "float isrf_energy_density_ not implemented, exiting"<<endl;
   exit(0);
}
*/
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double e_loss_compton_cc(double w,double gam)
{ 
   return E_LOSS_COMPTON_F77(&w,&gam);
}
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int test_e_loss_compton_cc()
{
   cout<<">>>>test_e_loss_compton_cc"<<endl;

   double w,gam,eph;
   eph=1.0; // 1 eV 
   gam=10;

   for(double emev=10;emev<1.e7;emev*=2)
   {
      w=eph /.511e6;   // units of mc^2
      gam=emev/.511;   // gamma of electron
      cout<<"target photon energy in eV="<<eph
         <<"       electron energy in MeV="<<emev
         <<" e_loss_compton_cc="<<e_loss_compton_cc(w,gam)<<endl;
   }
   cout<<"<<<<test_e_loss_compton_cc"<<endl;
   return 0;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double pp_meson_cc(double Esec, double Pp1, int NA1, int NA2, int key1)
{ 
   return PP_MESON_F77(&Esec, &Pp1, &NA1, &NA2, &key1);
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int test_pp_meson_cc()
{
   cout<<">>>>test_pp_meson_cc"<<endl;
   int NA1=1, NA2=1, key1=1;
   for(double Esec=10; Esec<1.e4; Esec*=10)// TOTAL electron energy in MeV
   {
      cout<<endl;
      for(double Pp1=10; Pp1<1.e7; Pp1*=2)
         cout<<"Esec="<<Esec<<" Pp1="<<Pp1<<" NA1="<<NA1<<" NA2="<<NA2
            <<" key1"<<key1<<" pp_meson_cc="
            <<pp_meson_cc( Esec, Pp1, NA1, NA2, key1 )<<endl;
   }
   cout<<"<<<<test_pp_meson_cc"<<endl;
   return 0;
}

