
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_synch_emiss.cc *                          galprop package * 5/02/2002 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// generate synchrotron emissivity
/*
CR density  gcr.cr_density is in c/4pi * n(E) [ cm s^-1 sr^-1 cm^-3 MeV^-1]

emissivity (erg cm^-3 sr^-1 Hz^-1 s^-1)=
 (1/4pi)integral[emissivity{Eelectron,nu}  ) n(E)E dlog(E)]  AWS20050727
=  (1/c)integral[emissivity{Eelectron,nu}  ) (c/4pi)n(E)E dlog(E)] 

since emmissivity from synchrotron_cc in erg s^-1 Hz^-1
factor=  (1/c)*log(Ekin_factor)
*/
using namespace std;//AWS20050624


#include"galprop_classes.h"
#include"galprop_internal.h"

#include <fort_interface.h>

#include"synchrotron_emissivity_aws.h" // for new B field model    AWS20080227

int Galprop::gen_synch_emiss()
{
   cout<<"gen_synch_emiss"<<endl;
   cout<<"generating synchrotron    emissivity for n_spatial_dimensions="<<gcr[0].n_spatial_dimensions<<endl;

   int stat=0, i, ielectrons, ipositrons;
   Distribution electrons;

// identify the electrons/positrons IMOS20030217
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


   for(int inusynch=0; inusynch<galaxy.n_nu_synchgrid; inusynch++)
   {
      for(int ip=0; ip<gcr[ielectrons].n_pgrid; ip++)
      {

	 ////////////////////////////////////////// 2D

         if(gcr[ielectrons].n_spatial_dimensions==2)
         {
            for(int ir=0; ir<gcr[ielectrons].n_rgrid; ir++)
            {
               for(int iz=0; iz<gcr[ielectrons].n_zgrid; iz++)
               {     
		 double synch_per_electron;                                                                                                       // AWS20080227 
                                                                                                    
		 if(galdef.synchrotron==1)  
		   {                                                                                                    // AWS20080227            
                     synch_per_electron =synchrotron_cc(gcr[ielectrons].gamma[ip], // Tesla ->Gauss                                               // AWS20080227
                     galaxy.nu_synch[inusynch], galaxy.B_field.d2[ir][iz].s[0]*10000.);// IMOS20020429

                    if(galdef.verbose==-1102)// selectable debug  
		    {
		     int debug=0;
                   
		     /*
                     double  synch_per_electron_new =synchrotron_emissivity_aws(gcr[ielectrons].gamma[ip],                     
			 galaxy.nu_synch[inusynch], galaxy.r[ir],galaxy.z[iz],galdef.B_field_model,debug);  
		     */

                     double  synch_per_electron_new =synchrotron_emissivity_aws(gcr[ielectrons].gamma[ip], 
			 galaxy.nu_synch[inusynch], galaxy.r[ir],galaxy.z[iz],galdef.B_field_name, galdef.B_field_parameters,debug);  //AWS20080313 


                                                                       // AWS20080227
		     cout<<"gamma nu r z "<<gcr[ielectrons].gamma[ip] <<" "<< galaxy.nu_synch[inusynch]<<" "<< galaxy.r[ir]<<" "<<galaxy.z[iz]
                      <<" galaxy.B_field (old)="<<galaxy.B_field.d2[ir][iz].s[0]
		      <<"  sync_per_electron(old)= "<<synch_per_electron<<" (new)= "<<synch_per_electron_new
                      <<" (new-old)/old ="<<(synch_per_electron_new- synch_per_electron)/synch_per_electron <<  endl;
                    }// if verbose

		   }//if galdef.synchrotron


		  if(galdef.synchrotron==2)                                                                           // new synchrotron and B field AWS20080227
		  {
		     int debug=0;
                      if(galdef.verbose==-1100)debug=1;            // selectable debug
		      /*                                                                        
                         synch_per_electron =synchrotron_emissivity_aws(gcr[ielectrons].gamma[ip],                     
			 galaxy.nu_synch[inusynch], galaxy.r[ir],galaxy.z[iz],galdef.B_field_model,debug);        
		      */
			 synch_per_electron =synchrotron_emissivity_aws(gcr[ielectrons].gamma[ip],                                    //AWS20080313                    
			 galaxy.nu_synch[inusynch], galaxy.r[ir],galaxy.z[iz],galdef.B_field_name, galdef.B_field_parameters,debug);  //AWS20080313  
   

	          }

 		  if(galdef.synchrotron==3)                                                                           // new synchrotron with Stokes AWS20100707
		  {
		     int debug=0;
                      if(galdef.verbose==-1100)debug=1;            // selectable debug
		     
		      double I,Q,U; // Stokes parameters
		      synch_per_electron =synchrotron_emissivity_aws(gcr[ielectrons].gamma[ip],                                                                     
									galaxy.nu_synch[inusynch], galaxy.r[ir],galaxy.z[iz],galdef.B_field_name, galdef.B_field_parameters,
                                                                        I,Q,U,debug);                 

                      galaxy.synchrotron_Q_emiss.d2[ir][iz].s[inusynch] += Q 
	                              *electrons.d2[ir][iz].s[ip]*gcr[ielectrons].Ekin[ip]; 


                      galaxy.synchrotron_U_emiss.d2[ir][iz].s[inusynch] += U 
	                              *electrons.d2[ir][iz].s[ip]*gcr[ielectrons].Ekin[ip]; 

	          }

                  galaxy.synchrotron_emiss.d2[ir][iz].s[inusynch] +=synch_per_electron 
	             *electrons.d2[ir][iz].s[ip]*gcr[ielectrons].Ekin[ip];             // IMOS20020429


                if(galdef.verbose==-1100 || galdef.verbose==-1101)// selectable debug                                                                        // AWS20080227
		  cout<<"gamma nu r z sync_per_electron  "<<gcr[ielectrons].gamma[ip] <<" "<< galaxy.nu_synch[inusynch]<<" "<< galaxy.r[ir]<<" "<<galaxy.z[iz]
                      <<" galaxy.B_field (old)="<<galaxy.B_field.d2[ir][iz].s[0]
		      <<"   "<<synch_per_electron<<"cumulative sync_emiss="<< galaxy.synchrotron_emiss.d2[ir][iz].s[inusynch] <<endl;

//cout<<"ir iz  E_gamma bremss_emiss "<<ir<<" "<<iz<<" "<<" "<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.bremss_emiss.d2[ir][iz].s[iEgamma]<<endl;  
               }  //  iz
            }  //  ir
         }  //  particle.n_spatial_dimensions==2
 


	 ////////////////////////////////////////// 3D 



	 double synch_per_electron;

         if(gcr[ielectrons].n_spatial_dimensions==3)
         {
            for(int ix=0; ix<gcr[ielectrons].n_xgrid; ix++)
            {
               for(int iy=0; iy<gcr[ielectrons].n_ygrid; iy++)
               {
                  for(int iz=0; iz<gcr[ielectrons].n_zgrid; iz++)
                  {


		 if(galdef.synchrotron==1)  
		 { 
                        synch_per_electron =synchrotron_cc(gcr[ielectrons].gamma[ip],     // Tesla ->Gauss
                        galaxy.nu_synch[inusynch], galaxy.B_field.d3[ix][iy][iz].s[0]*10000.);// IMOS20020429

                     galaxy.synchrotron_emiss.d3[ix][iy][iz].s[inusynch] += synch_per_electron
		  	                *electrons.d3[ix][iy][iz].s[ip]*gcr[ielectrons].Ekin[ip];             // IMOS20020429
		 }

		 if(galdef.synchrotron==2)  
		 { 

		     int debug=0;
                      if(galdef.verbose==-1100)debug=1;            // selectable debug

		      /*
                         synch_per_electron =synchrotron_emissivity_aws(gcr[ielectrons].gamma[ip],                     
					        galaxy.nu_synch[inusynch], galaxy.x[ix ],galaxy.y[iy],galaxy.z[iz],galdef.B_field_model,debug);
		      */

		      synch_per_electron =synchrotron_emissivity_aws(gcr[ielectrons].gamma[ip],                                                                   //AWS20080313
					 galaxy.nu_synch[inusynch], galaxy.x[ix ],galaxy.y[iy],galaxy.z[iz],galdef.B_field_name,galdef.B_field_parameters,debug); //AWS20080313

                     galaxy.synchrotron_emiss.d3[ix][iy][iz].s[inusynch] += synch_per_electron
		  	                *electrons.d3[ix][iy][iz].s[ip]*gcr[ielectrons].Ekin[ip];   
		 } 


		  if(galdef.synchrotron==3)                                                                           // new synchrotron with Stokes AWS20100707
		  {
		     int debug=0;
                      if(galdef.verbose==-1100)debug=1;            // selectable debug
		     
		      double I,Q,U; // Stokes parameters
		      synch_per_electron =synchrotron_emissivity_aws(gcr[ielectrons].gamma[ip],                                                                     
					            	galaxy.nu_synch[inusynch],  galaxy.x[ix ],galaxy.y[iy],galaxy.z[iz]  ,galdef.B_field_name, galdef.B_field_parameters,
                                                                        I,Q,U,debug);  
               
                      galaxy.synchrotron_emiss.d3[ix][iy][iz].s[inusynch] += synch_per_electron
		  	            *electrons.d3[ix][iy][iz].s[ip]*gcr[ielectrons].Ekin[ip]; 

                      galaxy.synchrotron_Q_emiss.d3[ix][iy][iz].s[inusynch] += Q
		  	              *electrons.d3[ix][iy][iz].s[ip]*gcr[ielectrons].Ekin[ip]; 

                      galaxy.synchrotron_U_emiss.d3[ix][iy][iz].s[inusynch] += U
		  	              *electrons.d3[ix][iy][iz].s[ip]*gcr[ielectrons].Ekin[ip]; 
	          }


              if(galdef.verbose==-1100 || galdef.verbose==-1101)// selectable debug                                                                        // AWS20080227
		  cout<<"gamma nu x y z  "<<gcr[ielectrons].gamma[ip] <<" "<< galaxy.nu_synch[inusynch]
                      <<" "<< galaxy.x[ix]<<" "<< galaxy.y[iy]<<" " <<galaxy.z[iz]
                      <<" galaxy.B_field (old)="<<galaxy.B_field.d3[ix][iy][iz].s[0]
		      <<" sync per electron=  "<<synch_per_electron
                      <<" gcr electrons ="<<electrons.d3[ix][iy][iz].s[ip]
                      <<" gcr electrons Ekin="<<gcr[ielectrons].Ekin[ip]
                      <<" cumulative sync_emiss="<< galaxy.synchrotron_emiss.d3[ix][iy][iz].s[inusynch] <<endl;

                  }  //  iz
               }  //  iy
            }  //  ix


 

         }  //  particle.n_spatial_dimensions==3
      }  //  ip
   }  //  inusynch

   double factor= log(galdef.Ekin_factor)/C; //AWS20050727

   galaxy.synchrotron_emiss   *= factor;
   galaxy.synchrotron_Q_emiss *= factor; //AWS20100709
   galaxy.synchrotron_U_emiss *= factor; //AWS20100709

   if(galdef.verbose==-801)// AWS20050727
   {
      cout<<"   synchrotron    emissivity "<<endl;
      galaxy.synchrotron_emiss.print();
   }  //  galdef.verbose
   electrons.delete_array();  // IMOS20020429
   cout<<" <<<< gen_synch_emiss"<<endl;
   return stat;
}
