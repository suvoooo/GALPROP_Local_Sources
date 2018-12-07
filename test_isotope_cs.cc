
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * test_isotope_cs.cc *                          galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

//float isotope_cs(double EMeV,int IZ,int IA, int IZF, int IAF, int kopt, int *info);
//float nucdata(int KSP,int IZ,int IA,int IZF,int IAF,int *IZL,int *IAL,double *T_HALF);
//void decayed_cross_sections(int iz,int ia,int jz,int ja, double *p,int np,double *sigma);

using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"


int Galprop::test_isotope_cs()
{ 
   int IZ = 2,IA = 4;
   int IZF= 2,IAF= 3;
   int kopt=11; kopt=2;
   int info;
   double CSmb;
   double EMeV=100;

   cout<<"test_isotope_cs_cc"<<endl;
   const string globalDataPath = DATA_PATH;
   read_nucdata(globalDataPath);

   IZ=6;IA=12;
   IZF=5;IAF=11;
   CSmb=isotope_cs(EMeV,IZ,IA,IZF,IAF,kopt,&info);
   cout<<"Z,A "<<IZ<<","<<IA<<"->"<<IZF<<","<<IAF<<" kopt="<<kopt<<" EMeV="
      <<EMeV<<" CSmb="<<CSmb<<" info="<<info<<endl;

   IZ=6 ;IA=12;
   IZF=6;IAF=11;
   CSmb=isotope_cs(EMeV,IZ,IA,IZF,IAF,kopt,&info);
   cout<<" EMeV="<<EMeV<<" CSmb="<<CSmb<<" info="<<info<<endl;

   IZ=6 ;IA=12;
   IZF=6;IAF=12;
   CSmb=isotope_cs(EMeV,IZ,IA,IZF,IAF,kopt,&info);
   cout<<" EMeV="<<EMeV<<" CSmb="<<CSmb<<" info="<<info<<endl;

   cout<<"====== systematic test of nuc_package cross sections ==="<<endl;

   IZ=6 ;IA=12;
   IZF=5;IAF=10;
   for(EMeV=10;EMeV<100000.;EMeV*=2.0)
   {
      kopt=10;
      CSmb=isotope_cs(EMeV,IZ,IA,IZF,IAF,kopt,&info);
      cout<<"Z,A "<<IZ<<","<<IA<<"->"<<IZF<<","<<IAF<<" kopt="<<kopt
         <<" EMeV="<<EMeV<<" CSmb="<<CSmb<<" info="<<info<<endl;

      kopt=11;
      CSmb=isotope_cs(EMeV,IZ,IA,IZF,IAF,kopt,&info);
      cout<<"Z,A "<<IZ<<","<<IA<<"->"<<IZF<<","<<IAF<<" kopt="<<kopt
         <<" EMeV="<<EMeV<<" CSmb="<<CSmb<<" info="<<info<<endl;
   }
   cout<<" ---------------- "<<endl;

   IZ=6 ;IA=12;
   IZF=5;IAF=11;
   for(EMeV=10;EMeV<100000.;EMeV*=2.0)
   {
      kopt=10;
      CSmb=isotope_cs(EMeV,IZ,IA,IZF,IAF,kopt,&info);
      cout<<"Z,A "<<IZ<<","<<IA<<"->"<<IZF<<","<<IAF<<" kopt="<<kopt
         <<" EMeV="<<EMeV<<" CSmb="<<CSmb<<" info="<<info<<endl;

      kopt=11;
      CSmb=isotope_cs(EMeV,IZ,IA,IZF,IAF,kopt,&info);
      cout<<"Z,A "<<IZ<<","<<IA<<"->"<<IZF<<","<<IAF<<" kopt="<<kopt
         <<" EMeV="<<EMeV<<" CSmb="<<CSmb<<" info="<<info<<endl;
   }
   cout<<" ---------------- "<<endl;
 
   IZ=6 ;IA =12;
   IZF=4;IAF= 9;
   for(EMeV=10;EMeV<100000.;EMeV*=2.0)
   {
      kopt=10;
      CSmb=isotope_cs(EMeV,IZ,IA,IZF,IAF,kopt,&info);
      cout<<"Z,A "<<IZ<<","<<IA<<"->"<<IZF<<","<<IAF<<" kopt="<<kopt
         <<" EMeV="<<EMeV<<" CSmb="<<CSmb<<" info="<<info<<endl;

      kopt=11;
      CSmb=isotope_cs(EMeV,IZ,IA,IZF,IAF,kopt,&info);
      cout<<"Z,A "<<IZ<<","<<IA<<"->"<<IZF<<","<<IAF<<" kopt="<<kopt
         <<" EMeV="<<EMeV<<" CSmb="<<CSmb<<" info="<<info<<endl;
   }
   cout<<" ---------------- "<<endl;

   IZ=6 ;IA =12;
   IZF=4;IAF=10;
   for(EMeV=10;EMeV<100000.;EMeV*=2.0)
   {
      kopt=10;
      CSmb=isotope_cs(EMeV,IZ,IA,IZF,IAF,kopt,&info);
      cout<<"Z,A "<<IZ<<","<<IA<<"->"<<IZF<<","<<IAF<<" kopt="<<kopt
         <<" EMeV="<<EMeV<<" CSmb="<<CSmb<<" info="<<info<<endl;

      kopt=11;
      CSmb=isotope_cs(EMeV,IZ,IA,IZF,IAF,kopt,&info);
      cout<<"Z,A "<<IZ<<","<<IA<<"->"<<IZF<<","<<IAF<<" kopt="<<kopt
         <<" EMeV="<<EMeV<<" CSmb="<<CSmb<<" info="<<info<<endl;
   }
   cout<<" ---------------- "<<endl;

   cout<<"testing nucdata_cc\n";

   int KSP,IZL,IAL;
   double T_HALF,B;
   int K_electron=0; //AWS20010731

   KSP=1;

   IZ=4;IA=10;
   IZF=5;IAF=10;
   B=nucdata(KSP,IZ,IA,K_electron,IZF,IAF,&IZL,&IAL,&T_HALF); // IMOS20010816
   cout<<KSP<<" "<<IZ<<" "<<IA<<" "<<IZF<<" "<<IAF<<" "<<IZL<<" "<<IAL<<" "
      <<T_HALF<<" "<<B<<endl;

   IZ=5;  IA=9 ;
   IZF=4;IAF=9 ;
   B=nucdata(KSP,IZ,IA,K_electron,IZF,IAF,&IZL,&IAL,&T_HALF);// IMOS20010816
   cout<<KSP<<" "<<IZ<<" "<<IA<<" "<<IZF<<" "<<IAF<<" "<<IZL<<" "<<IAL<<" "
      <<T_HALF<<" "<<B<<endl;

   KSP=0;

   IZ=5;  IA=9 ;
   IZF=4;IAF=9 ;
   B=nucdata(KSP,IZ,IA,K_electron,IZF,IAF,&IZL,&IAL,&T_HALF);// IMOS20010816
   cout<<KSP<<" "<<IZ<<" "<<IA<<" "<<IZF<<" "<<IAF<<" "<<IZL<<" "<<IAL<<" "
   <<T_HALF<<" "<<B<<endl;

   cout<<endl<<" full test of nucdata"<<endl;
   cout<<"IZ  IA IZF IAF IAL IZL t_half branching_ratio"<<endl;
   for(IZ =1;IZ <=28;       IZ ++)
   for(IA =1;IA <=2*IZ+10;  IA ++)
   for(IZF=1;IZF<=IZ+2;     IZF++)
   for(IAF=1;IAF<=IA  ;     IAF++)
     {
   B=nucdata(KSP,IZ,IA,K_electron,IZF,IAF,&IZL,&IAL,&T_HALF);// IMOS20010816
   if(B>0.)
     {
   cout<<IZ<<" "<<IA<<" "<<IZF<<" "<<IAF<<"         "<<IZL<<" "<<IAL<<"           "
   <<T_HALF<<" "<<B;
   if(IAL>0)cout<<"  longlived intermediate state!";
   cout<<endl;
     }
     }
   
   KSP=0;

   cout<<endl<<"testing  decayed_cross_sections\n";

   int iz=6;int ia=12;
   int jz=5;int ja=11;
   double p[]={50.,200.,300.,400.,500.,1000.,2000.,3000.,5000.,10000.}; int np=10;
   double sigma[10];

   decayed_cross_sections(iz,ia,jz,ja,p, np,sigma);//AWS20010731

   cout<<iz <<" "<<ia<<" "<<jz<<" "<<ja <<" "<<sigma<<endl;
   cout<<"    p=";for(int ip=0;ip<np;ip++)cout<<" "<<    p[ip];cout<<endl;
   cout<<"sigma=";for(int ip=0;ip<np;ip++)cout<<" "<<sigma[ip];cout<<endl;
   cout<<"test_isotope_cs_cc"<<endl;
   return 0;
}





/*
      IZ = 2
      IA = 4
cc      IZF= 2
      IAF= 3
      do i =10,2500,20
         Emev = i
         IZF= 2
         kopt = 11
         call ISOTOPE_CS(Emev,IZ,IA,IZF,IAF,kopt,info,CS10)
         IZF= 1
         kopt = 11
         call ISOTOPE_CS(Emev,IZ,IA,IZF,IAF,kopt,info,CS11)
         write(11,*)Emev,CS10,CS11, CS10+CS11
         write( *,*)Emev,CS10,CS11, CS10+CS11
      enddo
      stop


 
      subroutine ISOTOPE_CS(Emev,IZ,IA,IZF,IAF,kopt,info,CSmb)
c***********************************************************************
c                                *** I.Moskalenko, ver. 1 June, 1999 ***
c For any given primary nucleus (IZ,IA) and any final nucleus (IZF,IAF) 
c calculates the cross section of the reaction p+(Zi,Ai) -> (Zf,Af) +X.
c    INPUT:
c Emev  - energy of the primary nucleus in MeV/nucleon;
c IZ,IA - primary charge and atomic number;
c IZF,IAF - final charge and atomic number;
c kopt  =0 - uses the alghorithm described below;
c       =1 - forces to use Webber's code (no renormalization etc.);
c       =2 - forces to use ST98 code (no renormalization etc.);
c       =3 - forces to use a const cross section fitted to the data.
c       =10- forces to use Webber's code (renormalized if data exists);
c       =11- forces to use cross section fit if exists (otherwise equiv. 10);
c       =20- forces to use ST98 code (renormalized if data exists).
c       =21- forces to use cross section fit if exists (otherwise equiv. 20).
c    OUTPUT:
c info  =0 if no data exists, the Webber's or ST98 evaluation is used;
c       =1 if data exist and the renormalized Webber's formulae were used;
c       =2 if data exist and the renormalized ST98 formulae were used;
c       =3 if data exist and a const cross section is fitted to the data;
c       =10 Webber's code used in all cases (renormalized if data exists);
c       =20 ST98 code used in all cases (renormalized if data exists);
c       =11,21 if used the cross section fit;
c CSmb  - cross section of the reaction in mb.
c    COMMENTS:
c - When no data exists for the required channel, the Webber's evaluation is
c used. When Webber's code give 0, the result of ST98 code is used. 
c When the data exist, it calculates the renormalization coefficient for
c the two approximations (Webber's and ST98) using the least-square method.
c Then Xi2 value is used to choose the best approximation.
c - When both (Webber and ST98) approximations give 0, but the experimental
c value is non zero, uses a const cross section fitted to the data.
c - The CS_DATA array contents the following data:
c Zi.Ai, Zf.Af, Emev/n, CSmb, CSerr.
c if CSerr>0, it is the absolute error of the measurement;
c if CSerr<0, it is the relative error and the abs. err. will be calculated 
c automatically in the subroutine.
c    REFERENCES:
c "c$" - excluded data;
*/
