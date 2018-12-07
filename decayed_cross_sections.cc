
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * decayed_cross_sections.cc *                   galprop package * 8/16/2001 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// decayed_cross_sections          ### A.Strong/I.Moskalenko ###  8/16/2001 ###
// Calculates the production cross section (iz,ia)->(jz,ja) summed over 
// intermediate unstable states excluding channels which lead to production 
// long-lived isotopes (others than specified)
// Input: 
// (iz,ia) and (jz,ja) initial and final isotopes, stable or long-lived
// Ekin --an array of energies, MeV/nucleon (np -size)
// Output:
// sigma --an array of cross sections, mb (np -size)
//=="====!===="====!===="====!===="====!===="====!===="====!===="====!===="====!

using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"

void Galprop::decayed_cross_sections(int iz,int ia,int jz,int ja, double *Ekin,int np,double *sigma)
{
  int ip, IZ1,IA1,IZ3,IA3,kopt,info, K_electron =0;
   int diz=3;                // maximum delta Z considered
   int galdef_network_par=0; // temporary solution; value 1 caused problems in nuc_package
   double branching_ratio,t_half,CSmb,CSratio,CStot_ratio;

   for(ip=0; ip<np; sigma[ip++]=0.);
   if(ia < ja) return;

//========================sigma_boron_dec_heinbach_simon
   if(galdef.cross_section_option>=100 
      && sigma_boron_dec_heinbach_simon( iz, ia, jz, ja, 1000.)> 0.0)
   {
      if(galdef.verbose>=1)cout<<"     using Heinbach-Simon cross sections for ("<<iz<<","<<ia
         <<")  -> ("<<jz<<","<<ja<<")"<<endl;
      for(ip=0; ip<np; ip++)
      {
         sigma[ip]=sigma_boron_dec_heinbach_simon( iz, ia, jz, ja, Ekin[ip]);
// here neglect difference for decaying channels
         He_to_H_CS(Ekin[ip]/1000.,iz,ia,jz,ja,&CSratio,&CStot_ratio);
         sigma[ip]*=(1.0+galdef.He_H_ratio*CSratio) ;
         if(galdef.verbose>=1)cout<<"     Ekin="<<Ekin[ip]<<" CSratio= "<<CSratio
            <<" Heinbach-Simon sigma="<<sigma[ip]<<endl;
      }
   return;
   }

//========================================= using NUCDATA
//  if(galdef_network_par==0) cout<<"network: NDS-I.V. Moskalenko"<<endl;
//  if(galdef_network_par==1) cout<<"network: Garcia_Munoz et al."<<endl;

   kopt= galdef.cross_section_option%100;
   if(galdef.verbose>=1)
   {
      cout<<" Cross-section algorithm: kopt="<<kopt;
      if( kopt==  0 )cout<<" Optimized choice";                 
      if( kopt==  1 )cout<<" Silberberg/Tsao cross sections";
      if( kopt==  2 )cout<<" Webber Kish and Schreier cross sections";
      if( kopt==  3 )cout<<" constant fitted to data";                    
      if( kopt== 10 )cout<<" Webber formula normalized to data";
      if( kopt== 20 )cout<<" ST     formula normalized to data";
      if( kopt== 11 )cout<<" Webber formula normalized to data or a fit if available";
      if( kopt== 21 )cout<<" ST     formula normalized to data or a fit if available";
      cout<<endl;
   }
  
// a loop over an intermediate state; final state must be as requested
   for (IZ1=(jz-diz>1) ? jz-diz: 1; IZ1<=iz && IZ1<=jz+diz; IZ1++)
      for (IA1=(2*IZ1-4>ja) ? 2*IZ1-4: ja; IA1<ia && IA1<=2.5*IZ1+4.2; IA1++)
      {
// channel selection procedure
         if(IA1<IZ1 || ia-IA1<iz-IZ1) continue;
	                                                           // IMOS20010816 line below
         branching_ratio=nucdata(galdef_network_par,IZ1,IA1,K_electron,jz,ja, &IZ3,&IA3,&t_half);
//cout<<iz<<" "<<ia<<" "<<IZ1<<" "<<IA1<<" "<<jz<<" "<<ja<<"  b="<<branching_ratio<<endl;

         if(branching_ratio==0.) continue;
         t_half=t_half/year2sec;

// skip if long-lived intermediate state
         if(t_half>=galdef.t_half_limit                            // IMOS20010816
            && 100*IZ1+IA1!=100*jz+ja && 100*IZ3+IA3!=100*jz+ja) continue;

         if(galdef.verbose>=1)
         {
	    cout<<endl<< "(" <<iz<<","<<ia<<")->>("<<jz<<","<<ja<<")";
            if(100*IZ1+IA1==100*jz+ja)     //# direct channel
               cout<<" direct channel found ("
	          <<IZ1<<","<<IA1<<")->("<<jz<<","<<ja<<")";
            else                           //# short-lived intermediate state
               if(t_half>0. && t_half<galdef.t_half_limit)         // IMOS20010816
                  cout<<" unstable channel found with intermediate ("
                     <<IZ1<<","<<IA1<<")->("<<IZ3<<","<<IA3<<")->("<<jz<<","<<ja<<")"
                     <<" t_half="<<t_half<< " branching="<<branching_ratio;
               else                        //# decays without intermediate state
                  cout<<" unstable channel found without intermediate ("
                     <<IZ1<<","<<IA1<<")->("<<jz<<","<<ja<<")"
                     <<" t_half="<<t_half<< " branching="<<branching_ratio;
            cout<<endl;
         }
         if(galdef.verbose>0) cout<<"decayed_cross_sections: CSmb=";
         for(ip=0; ip<np; ip++)
         {
            CSmb=isotope_cs(Ekin[ip],iz,ia,IZ1,IA1,kopt,&info);
            if(galdef.verbose>0) cout<<CSmb<<" ";
            He_to_H_CS(Ekin[ip]/1000.,iz,ia,IZ1,IA1,&CSratio,&CStot_ratio);
            sigma[ip]+=CSmb*branching_ratio*(1.0+galdef.He_H_ratio*CSratio);
         }
	 if(galdef.verbose>=1){
	    cout<<endl<<"decayed_cross_sections:      sigma=";
	    for(ip=0; ip<np; cout<<" "<<sigma[ip++]); cout<<endl;
	 }
      } //ja1 //jz1
   if(galdef.verbose>=1){
      cout<<"TOTAL DECAYED SIGMA "<< "(" <<iz<<","<<ia<<")->>("<<jz<<","<<ja<<"):";
      for(ip=0; ip<np; cout<<" "<<sigma[ip++]); cout<<endl;
   }
}





