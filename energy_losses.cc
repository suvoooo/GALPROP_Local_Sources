
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * energy_losses.cc *                            galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|


//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// This file contains 2 c++ routines which evaluate nucleon and electron/positron
// energy losses; a sample driver program is in the end of the file.
//                                     ### Igor Moskalenko  ###  8/06/1999 ###
//
// Description of routines
// ^^^^^^^^^^^^^^^^^^^^^^^
//1 nucleon_loss(int,int,double,double,double,double,  double*,double*)
//1 calculates the NUCLEON energy losses in medium                                
//1 INPUT parameters:                                                           
//1   (int)    z,a  - the nucleon charge and the atomic number, correspondingly,
//1   (double) emev - the TOTAL NUCLEUS energy (MeV),                           
//1   (double) nhcm3,nhicm3 - the number densities of heutral H, and ionized H, 
//1   (double) he_to_h - the He/H ratio in the ISM (applyed for neutral gas only)
//1 OUTPUT:                                                                     
//1   returns total NUCLEUS energy loss   in eV/s,                              
//1   also partial energy losses (double* -pointers):                           
//1   aion [eV/s] - the ionization losses,                                      
//1   coul [eV/s] - the Coulomb energy losses.                                  
//
//2 electron_loss(double,double,double,double,double,double,
//2    double*,double*,double*,double*,double*,double*);
//2 calculates ELECTRON/POSITRON energy losses in medium (eV/s)               
//2 INPUT parameters:                                                           
//2   (double) emev - the total electron energy,                                
//2   (double) nhcm3,nhicm3 - the number densities of heutral H, and ionized H, 
//2   (double) he_to_h - the He/H ratio in the ISM (applyed for neutral and     
//2      ionized gas),                                                          
//2   (double) uevcm3,bevcm3 {=H^2/(8 pi)} - the energy density of photons and  
//2      magnetic field (H is the magnetic field streght).                      
//2 OUTPUT:                                                                     
//2   returns total electron energy losses in eV/s;                             
//2   also partial energy losses (double* -pointers):                           
//2   aion [eV/s] - the ionization losses (+ excitation + Cherenkov radiation), 
//2   coul [eV/s] - the Coulomb energy losses,                                  
//2   brem1, brem2 [eV/s] - the brems losses in neutral gas, and in plasma, corr.,
//2   cmptn [eV/s] - the Compton energy losses,                                 
//2   sync [eV/s] - the losses due to synchrotron emission.                     
// REFERENCES:                                                                 
// Strong A.W., Moskalenko I.V. 1998, ApJ 509, 212  and references therein    
//=="====!===="====!===="====!===="====!===="====!===="====!===="====!===="====!

using namespace std;//AWS20050624
#include <cmath>
#include <iostream>
#include "constants.h"

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double nucleon_loss(int z, int a, double emev, double nhcm3, 
   double nhicm3, double he_to_h,
   double* aion, double* coul) {

   double MA,BK,PIR0H2C2,PIR02MC2C,Te,nhe_cm3,gam,bet,qmax,bh,bhe,
      coullog,bet_e,xm,we;

// some constants
   MA = 1.e3*a*amu;                           // MeV, amu = C12/12 IMOS20061030
   BK = 1.38066e-23/1.60218e-13;              // MeV/K, Bolzmann constant
   PIR0H2C2  = Pi*Rele*H2PiC*H2PiC;           // MeV^2*cm^3
   PIR02MC2C = Pi*Rele*Rele*C*(Mele*1.e6);    // eV*cm^3/s = Pi*e^4/mc
   Te = 1.0e4;                                // K, warm ionized medium

   nhe_cm3 = he_to_h*nhcm3;                   // 1/cm^3, He number density
   gam = emev/MA;                             // nucleus Lorentz factor
   bet = sqrt(1.-1./(gam*gam));               // (=v/c) electron velosity

//## IONIZATION LOSSES in the neutral H and He
   qmax = 2.*Mele*(gam*gam-1.)/(1. +2.*gam*Mele/MA);
   bh   = log(2.*Mele*(gam*gam-1.)*qmax /(EH*EH))  -2.*bet*bet;
   bhe  = log(2.*Mele*(gam*gam-1.)*qmax /(EHe*EHe))-2.*bet*bet;
   *aion= 2.*PIR02MC2C*z*z/bet*(nhcm3*bh +nhe_cm3*bhe);

//## Coulomb energy losses in the cold H plasma limit
   coullog = -log(4.*PIR0H2C2*nhicm3 *(MA+2.*gam*Mele)
      /(4.*gam*gam*bet*bet*bet*bet*Mele*Mele*MA)) /2.;
   bet_e = sqrt(2.*BK*Te/Mele);
   xm = pow(3.*sqrt(Pi)/4.,1./3.) *bet_e;
   we = bet*bet*bet/(xm*xm*xm+bet*bet*bet);
   *coul = 4.*PIR02MC2C*z*z*nhicm3/bet *coullog *we;

// cout<<z<<" "<<a<<" "<<emev<<" "<<nhcm3<<" "<<nhicm3<<" "<<he_to_h<<" "<<*aion<<" "<<*coul<<endl;
//## TOTAL losses
   return (*aion +*coul);
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double electron_loss(double emev, double nhcm3, double nhicm3, double he_to_h, 
   double uevcm3, double bevcm3, 
   double* aion, double* coul, double* brem1, double* brem2, double* sync, 
   double* cmptn) {

   double PIR02C,PIR02MC2C,AFR02MC2C,ZHe,nhe_cm3,gam_e,bet_e,t_e,
      coullog2,MH,MHe,gam1,gam2;

// some constants
   PIR02C    = Pi*Rele*Rele*C;                // cm^3/s
   PIR02MC2C = Pi*Rele*Rele*C*(Mele*1.e6);    // eV*cm^3/s = Pi*e^4/mc
   AFR02MC2C = ALPHAf*Rele*Rele*C*(Mele*1.e6);// eV*cm^3/s = e^6/mc/hc
   ZHe = 2.;                                  // He charge

   nhe_cm3 = he_to_h*nhcm3;                   // 1/cm^3, He number density
   gam_e = emev/Mele;                         // electron Lorentz factor
   bet_e = sqrt(1.-1./(gam_e*gam_e));         // (=v/c) electron velosity
   t_e = gam_e-1;                             // electron kinetic energy, in mc^2 units IMOS20090617

//## IONIZATION LOSSES in the neutral H and He (Ginzburg 1979, p.360)
//   *aion=2.*PIR02MC2C/bet_e*(  (nhcm3+ZHe*nhe_cm3)*( log(gam_e-1.)-log(2.)+1./8.
//      +2.*log(gam_e*bet_e*Mele) )-nhcm3*2.*log(EiH)-ZHe*nhe_cm3*2.*log(EiHe)  );

// L.Pages, E.Bertel, H.Joffre, L.Sklavenitis 1972, Atomic Data 4, 1 IMOS20061221 
   double delta = 0.; // density effect =0 for now
                                                                // formula corrected by IMOS20090617
   *aion=2.*PIR02MC2C/bet_e*(  
	 (nhcm3+ZHe*nhe_cm3)*( log(pow(t_e,2)*(t_e+2.)/2.)
			       +( pow(t_e,2)/8.-(2.*t_e+1.)*log(2.) )/pow(t_e+1.,2)
			       +pow(gam_e,-2) -delta )
	 -nhcm3*2.*log(EiH/Mele)-ZHe*nhe_cm3*2.*log(EiHe/Mele) );

//## TEST ## OLD IONIZATION LOSSES in the neutral H and He (Ginzburg 1979,p.360)
/***************************************************************************
   *aion = 2.*PIR02MC2C*(  (nhcm3+ZHe*nhe_cm3)*( 3.*log(gam_e)+1./8.-log(2.)
      +2.*log(bet_e*Mele) )-nhcm3*2.*log(EiH)-ZHe*nhe_cm3*2.*log(EiHe)  );
***************************************************************************/

//## Coulomb energy losses in the cold H plasma limit (Ginzburg 1979, p.361)
   coullog2 = 0.;
   if(nhicm3 > 0.) 
      coullog2=log(emev/nhicm3*Mele/(4.*Pi*Rele*H2PiC*H2PiC))-3./4.;
   *coul = 2.*PIR02MC2C *nhicm3/bet_e *coullog2;

//## bremsstrahlung energy losses in neutral gas (Ginzburg 1979, p.386,409)
   MH = 1.67e-24;                            // grams, the mass of H atom 
   MHe= 1.66e-24 *4.;                        // grams, the mass of He atom
   gam1 = 100.;
   gam2 = 800.;
   if(gam_e < gam1) *brem1 = emev *4.*AFR02MC2C/Mele
      *(2.*nhcm3+ZHe*(ZHe+1.)*nhe_cm3)*(log(2.*gam_e)-1./3.);
   else 
      if(gam_e > gam2) *brem1 = 1.e6 *emev*C*(nhcm3*MH/TH +nhe_cm3*MHe/THe);
      else    // linear interpolation provides max 10% error to an exact value
         *brem1 = gam1*4.*AFR02MC2C
            *(2.*nhcm3+ZHe*(ZHe+1.)*nhe_cm3)*(log(2.*gam1)-1./3.)
            *(gam2-gam_e)/(gam2-gam1) +(gam_e-gam1)/(gam2-gam1) 
            *1.e6 *gam2*Mele *C*(nhcm3*MH/TH +nhe_cm3*MHe/THe);

//## bremsstrahlung energy losses in hydrogen plasma (Ginzburg 1979, p.408)
      *brem2 = emev*4.*2.*AFR02MC2C/Mele*nhicm3*(log(2.*gam_e)-1./3.);

//## TEST ## ACCURATE BREMSSTRAHLUNG CALCULATIONS
/***************************************************************************
//***  ep-bremsstrahlung losses (von J.Stickforth 1961)
   double brem_ep,ecm,pcm,qcmee;
   if(gam_e <= 2.) brem_ep = 8.*gam_e*bet_e*(1.-0.25*(gam_e-1.)
      +0.44935*pow(gam_e-1.,2.)-0.16577*pow(gam_e-1.,3.));
   else brem_ep = 1./bet_e*(6.*gam_e*log(2.*gam_e)-2.*gam_e-0.2900);
   brem_ep = brem_ep*2./3.;
//*** ee-bremsstrahlung losses (Haug 1975, Moskalenko & Jourdain 1997)
   ecm = sqrt((gam_e+1.)/2.);
   pcm = sqrt((gam_e-1.)/2.);
   qcmee = 8.*pcm*pcm/ecm*( 1.-4./3.*pcm/ecm
      +2./3.*(2.+pow(pcm/ecm,2.))*log(ecm+pcm) );
//*** TOTAL bremsstrahlung losses (ep- + ee-bremsstrahlung):
   *brem2 = AFR02MC2C *nhicm3*( qcmee*ecm+brem_ep );
***************************************************************************/

//## synchrotron energy losses (Ginzburg 1979,p.386)
   *sync  = 32./9.*PIR02C*bevcm3*(gam_e*gam_e-1.);

//## Compton energy losses in the Thomson limit
   *cmptn = 32./9.*PIR02C*uevcm3*(gam_e*gam_e-1.);

//## TOTAL
   return (*aion +*coul +*brem1 +*brem2 +*cmptn +*sync);
//     return (*brem1 +*brem2 +*cmptn +*sync);	// no loss due to sync
//   return *sync ;     	
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

/*
main()
{
   double emev,Ek,nhcm3=1.,nhicm3=1.,he_to_h=0,aion,coul;
   double brem1,brem2,sync,cmptn;
   int  z=1,a=1, i;
   for(Ek=1., i=0; i<40; i++, Ek/=1.2)
   {
      emev=a*(1.e3*Mp+Ek);
      nucleon_loss(z,a,emev,nhcm3,nhicm3,he_to_h,&aion,&coul);
      double gam=emev*1.e-3/Mp;
      double bet=sqrt(1.-1./(gam*gam));
      if(gam-1.<1.e-6) bet = sqrt(2.*(gam-1.));
      int H=0;
      if(bet>0.01) H=1;
      double eloss=1.82e-7*nhcm3*(1.+0.00185*log(bet)*H)*2.*pow(bet,2)/(1.e-6+2.*pow(bet,3));

      double eloss1=2.54e-9*a*a*nhcm3*sqrt(2./(gam-1.))*(log(gam-1.)+11.8);
      cout<<" p>> "<<Ek<<" "<<gam<<" "<<bet<<" > "<<aion<<" "<<eloss<<" "<<eloss1<<" "<<coul<<endl;
   }

   for(Ek=1., i=0; i<40; i++, Ek/=1.2)
   {
//Electron energy loss test for eprop test
      emev=(Mele+Ek);
      electron_loss(emev, nhcm3,nhicm3,0.,1.,0.,&aion,&coul,&brem1,&brem2,&sync,&cmptn);
      cout<<" e>> "<<Ek<<" "<<aion<<" "<<coul<<" "<<brem1<<" "<<brem2<<" "<<sync<<" "<<cmptn<<endl; // MeV/s
   }
}
*/
