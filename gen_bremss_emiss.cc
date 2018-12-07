
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_bremss_emiss.cc *                         galprop package * 2/17/2003 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// generate bremsstrahlung emissivity
/*
CR density  gcr.cr_density is in c/4pi * n(E) [ cm s^-1 sr^-1 cm^-3 MeV^-1]
emissivity (cm^-3 s^-1 sr^-1 MeV^-1)=
(c/4pi)*integral[sigma{Egamma,Eelectron}  ) n(E)E dlog(E)]
cross section from bremss_spec in barns MeV^-1
factor= 1.0e-24* log(Ekin_factor)
*/

#include <cassert>
#include <string>
#include <cstring>

using namespace std;//AWS20050624

#include "galprop_classes.h"
#include "galprop_internal.h"

#include <fort_interface.h>

#include <ErrorLogger.h>

int Galprop::gen_bremss_emiss()  {

  INFO("Entry");

  //cout<<" >>>> gen_bremss_emiss"<<endl;
  //cout<<"generating bremsstrahlung emissivity for n_spatial_dimensions="<<gcr[0].n_spatial_dimensions<<endl;

   double factor,cs_HII, cs_HI, cs_He, cs[2][3], y, Ekin=0., ye,fe;
   int stat=0, i, IZ1, Ne1, ielectrons, key,                      ippp =-1;
   Distribution electrons;

// identify the electrons/positrons 
   if(galdef.n_spatial_dimensions==2) electrons.init(gcr[0].n_rgrid,                 gcr[0].n_zgrid, gcr[0].n_pgrid);
   if(galdef.n_spatial_dimensions==3) electrons.init(gcr[0].n_xgrid, gcr[0].n_ygrid, gcr[0].n_zgrid, gcr[0].n_pgrid);
   electrons = 0.;
   for(i=0, ielectrons=-1; i<n_species; i++)  
      if(100==100*abs(gcr[i].Z)+gcr[i].A)
      {
         ielectrons=i;
 	 electrons+=gcr[ielectrons].cr_density;
         cout<<"  CR "<<gcr[ielectrons].name<<" found as species #"<<ielectrons<<endl;
      }
   if(ielectrons==-1) { cout<<"CR electrons/positrons not found!"<<endl; electrons.delete_array(); return 1; }

	 //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(cs_HII,cs_HI,cs_He,cs,y,ye,fe,Ne1,IZ1,key,i) firstprivate(Ekin)
   for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
   {
    
     //     cout << iEgamma << endl;

     for(i=0,key=0;i<3;i++)  cs[0][i] = cs[1][i] = 0.;
      for(int ip=0;ip<gcr[ielectrons].n_pgrid;ip++)
      {  // energy conservation check
	 //if(gcr[ielectrons].Etot[ip+1]<=1.02*Mele || gcr[ielectrons].Ekin[ip+1]-galaxy.E_gamma[iEgamma]<=0.) continue; //Gulli20070810 Moved to later

         if(galdef.integration_mode==1)
	 { // ### old integration ###
	 if(gcr[ielectrons].Etot[ip]<=1.02*Mele || gcr[ielectrons].Ekin[ip]-galaxy.E_gamma[iEgamma]<=0.) continue; //Gulli20070810
            IZ1=1; Ne1=0; //              ionized hydrogen H II
            cs_HII=bremss_spec_cc(galaxy.E_gamma[iEgamma],gcr[ielectrons].Etot[ip],IZ1,Ne1);
	    IZ1=1; Ne1=1; //              1-electron atom  H I
            cs_HI =bremss_spec_cc(galaxy.E_gamma[iEgamma],gcr[ielectrons].Etot[ip],IZ1,Ne1);
            IZ1=2; Ne1=2; //              2-electron atom  He
            cs_He =bremss_spec_cc(galaxy.E_gamma[iEgamma],gcr[ielectrons].Etot[ip],IZ1,Ne1);
            if(                cs_HI==0.) continue;
         } else
         { // ### new integration - analytical ###
	    if(ip>=gcr[ielectrons].n_pgrid-1) continue; // don't consider the last point since already included  //Gulli20070821
	 if(gcr[ielectrons].Etot[ip+1]<=1.02*Mele || gcr[ielectrons].Ekin[ip+1]-galaxy.E_gamma[iEgamma]<=0.) continue; //Gulli20070810
	    for(i=0;i<3;i++)  cs[0][i] = cs[1][i];   // reassign old points
// calculate new points in Etot
            IZ1=1; Ne1=0; //              ionized hydrogen H II
            cs[1][0]=bremss_spec_cc(galaxy.E_gamma[iEgamma],gcr[ielectrons].Etot[ip+1],IZ1,Ne1);
            IZ1=1; Ne1=1; //              1-electron atom  H I
            cs[1][1]=bremss_spec_cc(galaxy.E_gamma[iEgamma],gcr[ielectrons].Etot[ip+1],IZ1,Ne1);
            IZ1=2; Ne1=2; //              2-electron atom  He
            cs[1][2]=bremss_spec_cc(galaxy.E_gamma[iEgamma],gcr[ielectrons].Etot[ip+1],IZ1,Ne1);
// fool proof
            if(                cs[1][1]==0.) continue;
// lower integration limit falls between the grid points 
            if(cs[0][1]==0. && cs[1][1]!=0.) // test HI emiss. since HI has broadest distribution 
	    {
	       key=1;
	       Ekin =0.02*(gcr[ielectrons].Ekin[ip+1] +49.*galaxy.E_gamma[iEgamma]); //calc.lower int.limit
               IZ1=1; Ne1=0; //              ionized hydrogen H II
               cs[0][0]=bremss_spec_cc(galaxy.E_gamma[iEgamma],Ekin+Mele,IZ1,Ne1);
               IZ1=1; Ne1=1; //              1-electron atom  H I
               cs[0][1]=bremss_spec_cc(galaxy.E_gamma[iEgamma],Ekin+Mele,IZ1,Ne1);
               IZ1=2; Ne1=2; //              2-electron atom  He
               cs[0][2]=bremss_spec_cc(galaxy.E_gamma[iEgamma],Ekin+Mele,IZ1,Ne1);
//cout<<" lower integration limit falls between the grid points "<<cs[0][1]<<" "<<cs[1][1]<<" e= "<<gcr[ielectrons].Ekin[ip]<<" g= "<<galaxy.E_gamma[iEgamma]<<endl;
            }
         }
         if(gcr[ielectrons].n_spatial_dimensions==2)
            for(int ir=0; ir<gcr[ielectrons].n_rgrid-1; ir++)
               for(int iz=1; iz<gcr[ielectrons].n_zgrid-1; iz++)
               {
//************************************** TEST
if(galdef.verbose==-217)
{
if(galdef.integration_mode==1 && ip>0) Ekin=0.5*(gcr[ielectrons].Ekin[ip-1]+gcr[ielectrons].Ekin[ip  ]);
cs_HII =cs_HI =cs_He =1.;
for(i=0;i<3;i++)  cs[0][i] = cs[1][i] = 1.;
electrons.d2[ir][iz].s[ip  ] = pow(gcr[ielectrons].Ekin[ip  ],-3);
electrons.d2[ir][iz].s[ip+1] = pow(gcr[ielectrons].Ekin[ip+1],-3);
galdef.He_H_ratio=0.;
}
//***************************************/
                  if(galdef.integration_mode==1) // ### old integration ###
	          {
                     galaxy.bremss_emiss.        d2[ir][iz].s[iEgamma] +=
                       (cs_HI +cs_He*galdef.He_H_ratio) *electrons.d2[ir][iz].s[ip]*gcr[ielectrons].Ekin[ip]; //IMOS20020429
                     galaxy.bremss_ionized_emiss.d2[ir][iz].s[iEgamma] +=
                        cs_HII                          *electrons.d2[ir][iz].s[ip]*gcr[ielectrons].Ekin[ip]; //IMOS20020429
                     continue;
                  }
// ### new integration - analytical ###
		 if (electrons.d2[ir][iz].s[ip  ]>0.0&&  electrons.d2[ir][iz].s[ip+1]>0.0 ) //AWS20100503
		  {

		  if(key==1)  // case: lower integration limit falls between the grid points
		  { // interpolate electron spectrum (power-law)
                     ye=log(electrons.d2[ir][iz].s[ip  ]/  electrons.d2[ir][iz].s[ip+1])     // ye= power-law index
                       /log(  gcr[ielectrons].Ekin[ip  ]/    gcr[ielectrons].Ekin[ip+1]);
                     fe=    electrons.d2[ir][iz].s[ip  ]*pow(gcr[ielectrons].Ekin[ip  ],-ye);// normalization
                     fe*= pow(Ekin,ye);                                                      // electron flux @ Ekin
		  } else
		  {  // case: integration between the grid points
		     Ekin = gcr[ielectrons].Ekin[ip  ];
                     fe = electrons.d2[ir][iz].s[ip  ];
                  }

//************************************** TEST
//if(electrons.d2[ir][iz].s[ip]<0. && ippp!=ip) 
//  { for(i=0, ippp=ip;i<gcr[ielectrons].n_pgrid;i++) cout<<" "<<electrons.d2[ir][iz].s[i]; cout<<endl; }
//if(ir==0 && iz==3 && iEgamma==0) 
//cout<<" electrons @ ir,iz,iEgamma= 0  3  0 >>> "<<ip<<" "<<electrons.d2[ir][iz].s[ip]<<" "<<electrons.d2[ir][iz].s[ip+1]<<" "<<galaxy.bremss_emiss.d2[ir][iz].s[iEgamma]<<endl;
//***************************************/

// derive y (= power-law index of the total expression)
		 if (fe > 1e-40 && electrons.d2[ir][iz].s[ip+1] > 1e-40) { 
		  y =log((cs[0][1] +cs[0][2]*galdef.He_H_ratio) *fe
	               /((cs[1][1] +cs[1][2]*galdef.He_H_ratio) *electrons.d2[ir][iz].s[ip+1]))
                    /log(Ekin                              /gcr[ielectrons].       Ekin[ip+1]);
// integrate analytically
                  galaxy.bremss_emiss.        d2[ir][iz].s[iEgamma] +=(cs[1][1] +cs[1][2]*galdef.He_H_ratio) 
                              *electrons.     d2[ir][iz].s[ip+1]
                         *gcr[ielectrons].            Ekin[ip+1]/(y+1.)
                 *(1.-pow(                            Ekin/    
                          gcr[ielectrons].            Ekin[ip+1],y+1.));
								 double val = galaxy.bremss_emiss.d2[ir][iz].s[iEgamma];
								 if (isnan(val) || isinf(val)) {
									 cout<<y<<", "<<cs[0][1]<<", "<<cs[0][2]<<", "<<cs[1][1]<<", "<<cs[1][2]<<", "<<fe<<", "<<electrons.d2[ir][iz].s[ip+1]<<", "<<Ekin<<", "<<gcr[ielectrons].Ekin[ip+1]<<endl;
								 }

// derive y (= power-law index of the total expression)
	          y =log( cs[0][0]                              *fe
	               /( cs[1][0]                              *electrons.d2[ir][iz].s[ip+1]))
                    /log(Ekin                              /gcr[ielectrons].       Ekin[ip+1]);
// integrate analytically
                  galaxy.bremss_ionized_emiss.d2[ir][iz].s[iEgamma] += cs[1][0]                                   
                              *electrons.     d2[ir][iz].s[ip+1]
                         *gcr[ielectrons].            Ekin[ip+1]/(y+1.)
                 *(1.-pow(                            Ekin/    
                          gcr[ielectrons].            Ekin[ip+1],y+1.));
		 } else if ( fe > 1e-50) {
			 //Do linear interpolation, should be fine for this bin
			 galaxy.bremss_emiss.d2[ir][iz].s[iEgamma] += 0.5*(cs[0][1] + cs[0][2]*galdef.He_H_ratio)*fe/(gcr[ielectrons].Ekin[ip+1]-Ekin);
			 galaxy.bremss_ionized_emiss.d2[ir][iz].s[iEgamma] += 0.5*cs[0][0]*fe/(gcr[ielectrons].Ekin[ip+1]-Ekin);
		 } else if ( electrons.d2[ir][iz].s[ip+1] > 1e-50) {
			 //Do linear interpolation, should be fine for this bin
			 galaxy.bremss_emiss.d2[ir][iz].s[iEgamma] += 0.5*(cs[0][1] + cs[0][2]*galdef.He_H_ratio)*electrons.d2[ir][iz].s[ip+1]/(gcr[ielectrons].Ekin[ip+1]-Ekin);
			 galaxy.bremss_ionized_emiss.d2[ir][iz].s[iEgamma] += 0.5*cs[0][0]*electrons.d2[ir][iz].s[ip+1]/(gcr[ielectrons].Ekin[ip+1]-Ekin);
		 }

	      }//if electrons>0                 //AWS20100503
//cout<<" ir,iz,iEgamma= "<<ir<<"  "<<iz<<"  "<<iEgamma<<"  "<<" bremss_emiss= "<<galaxy.bremss_emiss.d2[ir][iz].s[iEgamma]<<endl;
               }//iz //ir //particle.n_spatial_dimensions==2
//************************************** TEST
if(galdef.verbose==-217)
if(iEgamma==galaxy.n_E_gammagrid/2)
  cout<<" Ekin= "<<Ekin<<"  "<<gcr[ielectrons].Ekin[ip]<<" bremss_emiss= "<<galaxy.bremss_emiss.d2[gcr[ielectrons].n_rgrid/2][gcr[ielectrons].n_zgrid/2].s[iEgamma]<<"  >>>"<<pow(Ekin,-2)/2.<<endl;
//***************************************/

         if(gcr[ielectrons].n_spatial_dimensions==3)
            for(int ix=1; ix<gcr[ielectrons].n_xgrid-1; ix++)
               for(int iy=1; iy<gcr[ielectrons].n_ygrid-1; iy++)
                  for(int iz=1; iz<gcr[ielectrons].n_zgrid-1; iz++)
                  {
                     if(galdef.integration_mode==1) // ### old integration ###
	             {
                        galaxy.bremss_emiss.        d3[ix][iy][iz].s[iEgamma] +=
                          (cs_HI +cs_He*galdef.He_H_ratio) *electrons.d3[ix][iy][iz].s[ip]*gcr[ielectrons].Ekin[ip]; //IMOS20020429
                        galaxy.bremss_ionized_emiss.d3[ix][iy][iz].s[iEgamma] +=
                           cs_HII                          *electrons.d3[ix][iy][iz].s[ip]*gcr[ielectrons].Ekin[ip]; //IMOS20020429
//cout<<"ix iy  iz  E_gamma bremss_emiss "<<ix<<" "<<iy<<" "<<iz<<" "<<" "
//<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.bremss_emiss.d3[ix][iy][iz].s[iEgamma]<<endl; 
                        continue;
                     }
// ### new integration - analytical ###
		     if (electrons.d3[ix][iy][iz].s[ip  ]>0.0&&  electrons.d3[ix][iy][iz].s[ip+1]>0.0 ) //AWS20100503
		       {
		     if(key==1)  // case: lower integration limit falls between the grid points
		     { // interpolate electron spectrum (power-law)
                        ye=log(electrons.d3[ix][iy][iz].s[ip  ]/  electrons.d3[ix][iy][iz].s[ip+1])     // ye= power-law index
                          /log(      gcr[ielectrons].Ekin[ip  ]/        gcr[ielectrons].Ekin[ip+1]);
                        fe=    electrons.d3[ix][iy][iz].s[ip  ]    *pow(gcr[ielectrons].Ekin[ip  ],-ye);// normalization
                        fe*= pow(Ekin,ye);                                                              // electron flux @ Ekin
		     } else
		     {  // case: integration between the grid points
		        Ekin =     gcr[ielectrons].Ekin[ip  ];
                        fe = electrons.d3[ix][iy][iz].s[ip  ];
                     }
// derive y (= power-law index of the total expression)
		     y =log((cs[0][1] +cs[0][2]*galdef.He_H_ratio) *fe
		          /((cs[1][1] +cs[1][2]*galdef.He_H_ratio) *electrons.d3[ix][iy][iz].s[ip+1]))
                       /log(Ekin                              /gcr[ielectrons].           Ekin[ip+1]);
// integrate analytically
                     galaxy.bremss_emiss.        d3[ix][iy][iz].s[iEgamma] +=(cs[1][1] +cs[1][2]*galdef.He_H_ratio) 
                                 *electrons.     d3[ix][iy][iz].s[ip+1]
                            *gcr[ielectrons].                Ekin[ip+1]/(y+1.)
                        *(1.-pow(                            Ekin/    
                             gcr[ielectrons].                Ekin[ip+1],y+1.));

// derive y (= power-law index of the total expression)
	             y =log( cs[0][0]                              *fe
	                  /( cs[1][0]                              *electrons.d3[ix][iy][iz].s[ip+1]))
                       /log(Ekin                              /gcr[ielectrons].           Ekin[ip+1]);
// integrate analytically
                     galaxy.bremss_ionized_emiss.d3[ix][iy][iz].s[iEgamma] += cs[1][0]                                   
                                 *electrons.     d3[ix][iy][iz].s[ip+1]
                           *gcr[ielectrons].                 Ekin[ip+1]/(y+1.)
                        *(1.-pow(                            Ekin/    
                            gcr[ielectrons].                 Ekin[ip+1],y+1.));

		     /*
		     if(isnan( galaxy.bremss_emiss.d3[ix][iy][iz].s[iEgamma] ))   //AWS20100503
		       cout<<"gen_bremss_emiss:  galaxy.bremss_emiss.d3[ix][iy][iz].s[iEgamma]="<< galaxy.bremss_emiss.d3[ix][iy][iz].s[iEgamma]
		       <<	 " ix="<<ix<<" iy="<<iy<<" iz="<<iz<<" ip="<<ip
                       <<" ye="<<ye<<" y="<<y
		       <<" electrons.d3[ix][iy][iz].s[ip  ]="<< electrons.d3[ix][iy][iz].s[ip  ]
		       <<" electrons.d3[ix][iy][iz].s[ip+1]="<< electrons.d3[ix][iy][iz].s[ip+1]
                       <<endl;
		     */

		       }//if electrons>0                 //AWS20100503


                  }//iz //iy //ix //particle.n_spatial_dimensions==3

	 if(key==1) key=2;
      }//ip
//************************************** TEST
if(galdef.verbose==-217)
if(iEgamma==galaxy.n_E_gammagrid/2) 
{
factor= log(galdef.Ekin_factor);
galaxy.bremss_emiss        *= factor;
if(galdef.integration_mode==1) cout<<"* integral= "<<galaxy.bremss_emiss.d2[gcr[ielectrons].n_rgrid/2][gcr[ielectrons].n_zgrid/2].s[iEgamma]<<endl;
//exit(1);
}
//***************************************/
   }//iEgamma
   
   factor= (galdef.integration_mode==1) ? 1.0e-24* log(galdef.Ekin_factor) : 1.0e-24;
   galaxy.bremss_emiss        *= factor;
   galaxy.bremss_ionized_emiss*= factor; //IMOS20020429
   electrons.delete_array();            // IMOS20020429

      if(galdef.verbose>=2)
   {
      cout<<"   bremsstrahlung emissivity "<<endl;
      galaxy.bremss_emiss.print();
   }//galdef.verbose

   INFO("Exit");

   return stat;
}
