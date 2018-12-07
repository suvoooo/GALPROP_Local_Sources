
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * protri.cc *                                  galprop package *  3/25/2002 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"
/*
extern "C" int tridag_ext    ( double*, double*, double*, double*, double*,int,int);
extern "C" int tridag_sym_ext( double*, double*, double*, double*, double*,int,int);
*/
//Gulli20070810 Changed from float to double
 
int tridag_ext    ( double*, double*, double*, double*, double*,int,int);
int tridag_sym_ext( double*, double*, double*, double*, double*,int,int);




#define NMAX 4000000


/*
double Nx1[NMAX],Nx2[NMAX],Nx3[NMAX],Nx0[NMAX],Rx[NMAX];
double Ny1[NMAX],Ny2[NMAX],Ny3[NMAX],Ny0[NMAX],Ry[NMAX];
double Nz1[NMAX],Nz2[NMAX],Nz3[NMAX],Nz0[NMAX],Rz[NMAX];
double Np1[NMAX],Np2[NMAX],Np3[NMAX],Np0[NMAX],Rp[NMAX];

double alpha1_xx[NMAX],alpha1_yy[NMAX],alpha1_zz[NMAX],alpha1_pp[NMAX];
double alpha2_xx[NMAX],alpha2_yy[NMAX],alpha2_zz[NMAX],alpha2_pp[NMAX];
double alpha3_xx[NMAX],alpha3_yy[NMAX],alpha3_zz[NMAX],alpha3_pp[NMAX];

double alpha1_x0[NMAX],alpha1_y0[NMAX],alpha1_z0[NMAX]; 

double total_source_function_x[NMAX];
double total_source_function_y[NMAX];
double total_source_function_z[NMAX];
double total_source_function_p[NMAX];
*/


double *Nx1      ,*Nx2      ,*Nx3      ,*Nx0      ,*Rx      ;
double *Ny1      ,*Ny2      ,*Ny3      ,*Ny0      ,*Ry      ;
double *Nz1      ,*Nz2      ,*Nz3      ,*Nz0      ,*Rz      ;
double *Np1      ,*Np2      ,*Np3      ,*Np0      ,*Rp      ;

double *alpha1_xx      ,*alpha1_yy      ,*alpha1_zz      ,*alpha1_pp      ;
double *alpha2_xx      ,*alpha2_yy      ,*alpha2_zz      ,*alpha2_pp      ;
double *alpha3_xx      ,*alpha3_yy      ,*alpha3_zz      ,*alpha3_pp      ;

double *alpha1_x0      ,*alpha1_y0      ,*alpha1_z0       ; 

float *total_source_function_x      ;
float *total_source_function_y      ;
float *total_source_function_z      ;
float *total_source_function_p      ;

int protri_init=0;


/////////////////////////////////////////////////////////////////
void Galprop::protri(Particle &particle,

Distribution  &alpha1_x,Distribution  &alpha1_y,Distribution  &alpha1_z,Distribution  &alpha1_p,
Distribution  &alpha2_x,Distribution  &alpha2_y,Distribution  &alpha2_z,Distribution  &alpha2_p,
Distribution  &alpha3_x,Distribution  &alpha3_y,Distribution  &alpha3_z,Distribution  &alpha3_p,

Distribution &Nx1_, Distribution &Ny1_, Distribution &Nz1_, Distribution &Np1_,
Distribution &Nx2_, Distribution &Ny2_, Distribution &Nz2_, Distribution &Np2_,
Distribution &Nx3_, Distribution &Ny3_, Distribution &Nz3_, Distribution &Np3_,

Distribution &total_source_function,




double dt, int nrept_outer ,double f_use 
    
             )
{

   double t;

   int ix,iy,iz,ip;
   int nk; // number of vectors to solve in tridag_ext
   int kk; // basic array index as in tridag_ext
   int kk1;// kk+1
   int nkt; //total number of elements in array
   int i,k,ns;//counters, stride

   int irept_outer ; // for outer loop
 
 
   int irept,nrept; // for inner loops within protri
   nrept=1;

cout<<">>protri"<<endl; 

 

nkt=particle.n_xgrid*particle.n_ygrid*particle.n_zgrid*particle.n_pgrid; 

 cout<<" nkt="<<nkt<<endl;

// if (protri_init==0){ //Gulli20070810
 cout<<" Assigning arrays with dimension="<<nkt<<endl;
     protri_init=1;


       Nx1=new double[nkt]; Nx2=new double[nkt]; Nx3=new double[nkt]; Nx0=new double[nkt]; Rx=new double[nkt];       
       Ny1=new double[nkt]; Ny2=new double[nkt]; Ny3=new double[nkt]; Ny0=new double[nkt]; Ry=new double[nkt];      
       Nz1=new double[nkt]; Nz2=new double[nkt]; Nz3=new double[nkt]; Nz0=new double[nkt]; Rz=new double[nkt];       
       Np1=new double[nkt]; Np2=new double[nkt]; Np3=new double[nkt]; Np0=new double[nkt]; Rp=new double[nkt];       

       alpha1_xx=new double[nkt];  alpha1_yy=new double[nkt];  alpha1_zz=new double[nkt]; alpha1_pp=new double[nkt];       
       alpha2_xx=new double[nkt];  alpha2_yy=new double[nkt];  alpha2_zz=new double[nkt]; alpha2_pp=new double[nkt];       
       alpha3_xx=new double[nkt];  alpha3_yy=new double[nkt];  alpha3_zz=new double[nkt]; alpha3_pp=new double[nkt];         

       alpha1_x0=new double[nkt];  alpha1_y0=new double[nkt];  alpha1_z0=new double[nkt];         

       total_source_function_x=new float[nkt];       
       total_source_function_y=new float[nkt];       
       total_source_function_z=new float[nkt];       
       total_source_function_p=new float[nkt];       

// }//protri_init==0


// X propagation

              if(galdef.prop_x==1){

               kk=0;
               
                  for(ip=0; ip<particle.n_pgrid; ip++)
                  {
                     for(iz=0; iz<particle.n_zgrid; iz++)
                     {
                        for(iy=0; iy<particle.n_ygrid; iy++)
                        {
                    //       if(galdef.use_symmetry<=1||(iy<=n_ygrid_sym && iz<=n_zgrid_sym))
                 //          {
                              for(ix=0; ix<particle.n_xgrid; ix++,kk++)
                              {
                                 Nx1[kk] = Nx1_.d3[ix][iy][iz].s[ip];
                                 Nx2[kk] = Nx2_.d3[ix][iy][iz].s[ip];
                                 Nx3[kk] = Nx3_.d3[ix][iy][iz].s[ip];
                                 alpha1_xx[kk]=alpha1_x.d3[ix][iy][iz].s[ip]; 
                                 alpha2_xx[kk]=alpha2_x.d3[ix][iy][iz].s[ip];
                                 alpha3_xx[kk]=alpha3_x.d3[ix][iy][iz].s[ip];

                                 total_source_function_x[kk]=  total_source_function.d3[ix][iy][iz].s[ip];
                                 Rx[kk]=particle.cr_density.d3[ix][iy][iz].s[ip];
                              }

                     //      }   //symmetry
                        }   //iy
                     }   //iz
                  }   //ip


                nk=particle.n_pgrid*particle.n_zgrid*particle.n_ygrid;

                // if symmetry, alpha1_xx at ix=0 is needed so copy it first
                if (galdef.use_symmetry==1)
                for(kk=0                 ;kk<nkt;kk+=particle.n_xgrid)alpha1_x0[kk]=alpha1_xx[kk];

                // zero coefficients to enable loops to be vectorized
                for(kk=0                 ;kk<nkt;kk+=particle.n_xgrid)alpha1_xx[kk]=0.;
                for(kk=particle.n_xgrid-1;kk<nkt;kk+=particle.n_xgrid)alpha3_xx[kk]=0.;
                



	


	      }//galdef.prop_x==1




// Y propagation

               if(galdef.prop_y==1)
               {
                  kk=0;

                  for(ix=0; ix<particle.n_xgrid; ix++)
                  {
                     for(ip=0; ip<particle.n_pgrid; ip++)
                     {
                        for(iz=0; iz<particle.n_zgrid; iz++)
                        {
			  //     if(galdef.use_symmetry<=1 || (iz<=n_zgrid_sym && ix<=n_xgrid_sym))
                          // {
                              for(iy=0; iy<particle.n_ygrid; iy++,kk++) 
                              {
                                 Ny1[kk] = Ny1_.d3[ix][iy][iz].s[ip];
                                 Ny2[kk] = Ny2_.d3[ix][iy][iz].s[ip];
                                 Ny3[kk] = Ny3_.d3[ix][iy][iz].s[ip];

                                 alpha1_yy[kk]=alpha1_y.d3[ix][iy][iz].s[ip]; 
                                 alpha2_yy[kk]=alpha2_y.d3[ix][iy][iz].s[ip];
                                 alpha3_yy[kk]=alpha3_y.d3[ix][iy][iz].s[ip];

                                 total_source_function_y[kk]=  total_source_function.d3[ix][iy][iz].s[ip];
                                 Ry[kk]=particle.cr_density.d3[ix][iy][iz].s[ip];

                              }
  
			      //      }  //  symmetry
 
                        }  //  iz
                     }  //  ip
                  }  //  ix	 
                                 nk=particle.n_pgrid*particle.n_zgrid*particle.n_xgrid;


                // if symmetry, alpha1_yy at iy=0 is needed so copy it first
                if (galdef.use_symmetry==1)
                for(kk=0                 ;kk<nkt;kk+=particle.n_ygrid)alpha1_y0[kk]=alpha1_yy[kk];

                // zero coefficients to enable loops to be vectorized
                for(kk=0                 ;kk<nkt;kk+=particle.n_ygrid)alpha1_yy[kk]=0.;
                for(kk=particle.n_ygrid-1;kk<nkt;kk+=particle.n_ygrid)alpha3_yy[kk]=0.;


	       }//galdef.prop_y==1


// Z propagation

               if(galdef.prop_z==1)
               {
                  kk=0;

                  for(iy=0; iy<particle.n_ygrid; iy++)
                  {
                     for(ix=0; ix<particle.n_xgrid; ix++)
                     {
                        for(ip=0; ip<particle.n_pgrid; ip++)
                        {
			  //         if(galdef.use_symmetry<=1 || (ix<=n_xgrid_sym && iy<=n_ygrid_sym))
                          // {
                              for(iz=0; iz<particle.n_zgrid; iz++,kk++) 
                              {
                                 Nz1[kk]=Nz1_.d3[ix][iy][iz].s[ip];
                                 Nz2[kk]=Nz2_.d3[ix][iy][iz].s[ip];
                                 Nz3[kk]=Nz3_.d3[ix][iy][iz].s[ip];

                                 alpha1_zz[kk]=alpha1_z.d3[ix][iy][iz].s[ip]; 
                                 alpha2_zz[kk]=alpha2_z.d3[ix][iy][iz].s[ip];
                                 alpha3_zz[kk]=alpha3_z.d3[ix][iy][iz].s[ip];

                                 total_source_function_z[kk]=  total_source_function.d3[ix][iy][iz].s[ip];
                                 Rz[kk]=particle.cr_density.d3[ix][iy][iz].s[ip];
                              } 
			      //     }  //  symmetry
 
 
                        }  //  ip
                     }  //  ix
                  }  //  iy
                                 nk=particle.n_pgrid*particle.n_xgrid*particle.n_ygrid;

                // if symmetry, alpha1_zz at iz=0 is needed so copy it first
                if (galdef.use_symmetry==1)
                for(kk=0                 ;kk<nkt;kk+=particle.n_zgrid)alpha1_z0[kk]=alpha1_zz[kk];

                // zero coefficients to enable loops to be vectorized
                for(kk=0                 ;kk<nkt;kk+=particle.n_zgrid)alpha1_zz[kk]=0.;
                for(kk=particle.n_zgrid-1;kk<nkt;kk+=particle.n_zgrid)alpha3_zz[kk]=0.;


	       }//galdef.prop_z==1


// P propagation

               if(galdef.prop_p==1)
               { 

                  kk=0;

                  for(iz=0; iz<particle.n_zgrid; iz++)
                  {
                     for(iy=0; iy<particle.n_ygrid; iy++)
                     {
                        for(ix=0; ix<particle.n_xgrid; ix++)
                        {
			  //    if(galdef.use_symmetry<=1 || (ix<=n_xgrid_sym && iy<=n_ygrid_sym && iz<=n_zgrid_sym))
			  // {
                              for(ip=0; ip<particle.n_pgrid; ip++,kk++) 
                              {
                                 Np1[kk]=Np1_.d3[ix][iy][iz].s[ip];
                                 Np2[kk]=Np2_.d3[ix][iy][iz].s[ip];
                                 Np3[kk]=Np3_.d3[ix][iy][iz].s[ip];

                                 alpha1_pp[kk]=alpha1_p.d3[ix][iy][iz].s[ip]; 
                                 alpha2_pp[kk]=alpha2_p.d3[ix][iy][iz].s[ip];
                                 alpha3_pp[kk]=alpha3_p.d3[ix][iy][iz].s[ip];

                                 total_source_function_p[kk]=  total_source_function.d3[ix][iy][iz].s[ip];
			         Rp[kk]=particle.cr_density.d3[ix][iy][iz].s[ip];

                              }


			      //   }  //  symmetry
 
                        }  //  ix
                     }  //  iy
                  }  //  iz


                                 nk=particle.n_xgrid*particle.n_ygrid*particle.n_zgrid;

                // zero coefficients to enable loops to be vectorized
                for(kk=0                 ;kk<nkt;kk+=particle.n_pgrid)alpha1_pp[kk]=0.;
                for(kk=particle.n_pgrid-1;kk<nkt;kk+=particle.n_pgrid)alpha3_pp[kk]=0.;

	       }// galdef.prop_p==1

	//---------------------------------------------------------------------------------------------------
        // end of setup phase
        // -------------------------------------------------------------------------------------------


        t = 0.0;


for (irept_outer=0;irept_outer<nrept_outer;irept_outer++){

    t+=dt;

    cout<<"   irept_outer="<<irept_outer<<" t="<<t<<endl;

    //#pragma vdir nodep
    for(kk=0;kk<nkt;kk++)total_source_function_x[kk]=0.0; // zero here to vectorize
    source_SNR_event_vec(particle, t/year2sec, total_source_function_x);

// X propagation

              if(galdef.prop_x==1){

		/*
              kk=0;

                  for(ip=0; ip<particle.n_pgrid; ip++)
                  {
                     for(iz=0; iz<particle.n_zgrid; iz++)
                     {
                        for(iy=0; iy<particle.n_ygrid; iy++)
                        {
                  //         if(galdef.use_symmetry<=1||(iy<=n_ygrid_sym && iz<=n_zgrid_sym))
                    //       {

			  for(ix=0; ix<particle.n_xgrid; ix++,           kk++){
			 
			      total_source_function_[kk]=  total_source_function.d3[ix][iy][iz].s[ip]; }
 
                     //       }   //symmetry
                        }   //iy
                     }   //iz
                  }   //ip
		*/
                nk=particle.n_pgrid*particle.n_zgrid*particle.n_ygrid;



 

                

              for(irept=0;irept<nrept;irept++){//cout<<irept<<endl;
               
             
              //#pragma vdir nodep
              for(kk=1;kk<nkt-1;kk++)   Nx0[kk] =     total_source_function_x[kk]*dt/f_use;
              //#pragma vdir nodep
              for(kk=1;kk<nkt-1;kk++)   Nx0[kk]+= (1.0-0.5*alpha2_xx[kk])*Rx[kk];
              //#pragma vdir nodep
              for(kk=1;kk<nkt-1;kk++)   Nx0[kk]+=      0.5*alpha1_xx[kk] *Rx[kk-1];
              //#pragma vdir nodep
              for(kk=1,kk1=2;kk<nkt-1;kk++,kk1++)   Nx0[kk]+=      0.5*alpha3_xx[kk] *Rx[kk1]; 
        
                        
                                    
               
                 kk=0;
                 Nx0[kk] = total_source_function_x[kk]*dt/f_use;
                 Nx0[kk]+= (1.0-0.5*alpha2_xx[kk])*Rx[kk];
                 Nx0[kk]+=      0.5*alpha3_xx[kk] *Rx[kk+1];
                
                 kk=nkt-1;
                 Nx0[kk] = total_source_function_x[kk]*dt/f_use;
                 Nx0[kk]+= (1.0-0.5*alpha2_xx[kk])*Rx[kk];
                 Nx0[kk]+=      0.5*alpha1_xx[kk] *Rx[kk-1];

           
                if(galdef.use_symmetry==0)
                 tridag_ext     (Nx1,Nx2,Nx3,Nx0,Rx, particle.n_xgrid,nk );
                  
                if(galdef.use_symmetry==1){
                 //#pragma vdir nodep
                 for(kk=0,kk1=1;kk<nkt-1;kk+=particle.n_xgrid,kk1+=particle.n_xgrid)Nx0[kk]+=0.5*alpha1_x0[kk] *Rx[kk1];

                 tridag_sym_ext (Nx1,Nx2,Nx3,Nx0,Rx, particle.n_xgrid,nk );
	                         	}
                         
                 } // irept
 
} // galdef.prop_x==1


//cout<<"   x complete"<<endl;

	      // permute to y-first
                ns=particle.n_xgrid;
                kk=0;
                //#pragma vdir nodep
                for(i=0;i<ns;i++)
                 //#pragma vdir nodep
                 for(k=i;k<nkt;k+=ns,kk++)Ry[kk]=Rx[k];
		kk=0;
                //#pragma vdir nodep
                for(i=0;i<ns;i++)
                 //#pragma vdir nodep
                 for(k=i;k<nkt;k+=ns,kk++)total_source_function_y[kk]=total_source_function_x[k]; 

// Y propagation

	      if(galdef.prop_y==1){
		/*           
                    kk=0;

                  for(ix=0; ix<particle.n_xgrid; ix++)
                  {
                     for(ip=0; ip<particle.n_pgrid; ip++)
                     {
                        for(iz=0; iz<particle.n_zgrid; iz++)
			  {
			    for(iy=0; iy<particle.n_ygrid; iy++,          kk++) {
			             
				   total_source_function_[kk]=  total_source_function.d3[ix][iy][iz].s[ip]; }
			    
 
                        }  //  iz
                     }  //  ip
                  }  //  ix
		*/

                nk=particle.n_pgrid*particle.n_zgrid*particle.n_xgrid;

              for(irept=0;irept<nrept;irept++){//cout<<irept<<endl;
               
                
                 //#pragma vdir nodep
                 for(kk=1;kk<nkt-1;kk++)Ny0[kk] = total_source_function_y[kk]*dt/f_use;
                 //#pragma vdir nodep
                 for(kk=1;kk<nkt-1;kk++)Ny0[kk]+= (1.0-0.5*alpha2_yy[kk])*Ry[kk];
                 //#pragma vdir nodep
                 for(kk=1;kk<nkt-1;kk++)Ny0[kk]+=      0.5*alpha1_yy[kk] *Ry[kk-1];
                 //#pragma vdir nodep
                 for(kk=1,kk1=2;kk<nkt-1;kk++,kk1++)Ny0[kk]+=      0.5*alpha3_yy[kk] *Ry[kk1]; 
        
                            
                                    
               
                 kk=0;
                 Ny0[kk] =     total_source_function_y[kk]*dt/f_use;
                 Ny0[kk]+= (1.0-0.5*alpha2_yy[kk])*Ry[kk];
                 Ny0[kk]+=      0.5*alpha3_yy[kk] *Ry[kk+1];
                
                 kk=nkt-1;
                 Ny0[kk] =     total_source_function_y[kk]*dt/f_use;
                 Ny0[kk]+= (1.0-0.5*alpha2_yy[kk])*Ry[kk];
                 Ny0[kk]+=      0.5*alpha1_yy[kk] *Ry[kk-1];

                if(galdef.use_symmetry==0)
                 tridag_ext     (Ny1,Ny2,Ny3,Ny0,Ry, particle.n_ygrid,nk );

                if(galdef.use_symmetry==1){
                 //#pragma vdir nodep
                 for(kk=0,kk1=1;kk<nkt-1;kk+=particle.n_ygrid,kk1+=particle.n_ygrid)Ny0[kk]+=0.5*alpha1_y0[kk] *Ry[kk1];

                 tridag_sym_ext (Ny1,Ny2,Ny3,Ny0,Ry, particle.n_ygrid,nk );
	                         	}

	      }//irept



                              
 
		  /*
                  if(galdef.use_symmetry==2)
                  {
                     for(ip=0; ip<particle.n_pgrid; ip++)
                     {
                        for(iz=0; iz<=n_zgrid_sym; iz++)
                        {
                           for(ix=0; ix<=n_xgrid_sym; ix++)
                           {
 
                              for(iy=0; iy<particle.n_ygrid; iy++)
                              {
                                 double value=
                                    particle.cr_density.d3[                   ix][iy][                   iz].s[ip];
                                    particle.cr_density.d3[particle.n_xgrid-1-ix][iy][                   iz].s[ip]=value;
                                    particle.cr_density.d3[                   ix][iy][particle.n_zgrid-1-iz].s[ip]=value;
                                    particle.cr_density.d3[particle.n_xgrid-1-ix][iy][particle.n_zgrid-1-iz].s[ip]=value;
                              }  //  iy    
                           }  //  ix     
                        }  //  iz
                     }  //  ip
		         }  //  symmetry
		  */
               }  //  prop_y

	      
//cout<<"  y complete"<<endl;
 

	      // permute to z-first
                ns=particle.n_ygrid;
                kk=0;
                //#pragma vdir nodep
                for(i=0;i<ns;i++)
                 //#pragma vdir nodep
                 for(k=i;k<nkt;k+=ns,kk++)Rz[kk]=Ry[k];
                kk=0;
                //#pragma vdir nodep
                for(i=0;i<ns;i++)
                 //#pragma vdir nodep
                 for(k=i;k<nkt;k+=ns,kk++)
                 //#pragma vdir nodep
                 total_source_function_z[kk]=total_source_function_y[k]; 



// Z propagation

	      if(galdef.prop_z==1){
  

		/*
                              kk=0;
                 for(iy=0; iy<particle.n_ygrid; iy++)
                  {
                     for(ix=0; ix<particle.n_xgrid; ix++)
                     {
                        for(ip=0; ip<particle.n_pgrid; ip++)
                        {

			  for(iz=0;iz<particle.n_zgrid;iz++,             kk++){
			     
                                 total_source_function_[kk]=  total_source_function.d3[ix][iy][iz].s[ip]; }
			    
 
 
                        }  //  ip
                     }  //  ix
                  }  //  iy
		*/
                nk=particle.n_pgrid*particle.n_xgrid*particle.n_ygrid;            

              for(irept=0;irept<nrept;irept++){//cout<<irept<<endl;
               
                
                 //#pragma vdir nodep
                 for(kk=1;kk<nkt-1;kk++)Nz0[kk] =     total_source_function_z[kk]*dt/f_use;
                 //#pragma vdir nodep
                 for(kk=1;kk<nkt-1;kk++)Nz0[kk]+= (1.0-0.5*alpha2_zz[kk])*Rz[kk];
                 //#pragma vdir nodep
                 for(kk=1;kk<nkt-1;kk++)Nz0[kk]+=      0.5*alpha1_zz[kk] *Rz[kk-1];
                 //#pragma vdir nodep
                 for(kk=1,kk1=2;kk<nkt-1;kk++,kk1++)Nz0[kk]+=      0.5*alpha3_zz[kk] *Rz[kk1]; 
        
                              
                                    
               
                 kk=0;
                 Nz0[kk] =     total_source_function_z[kk]*dt/f_use;
                 Nz0[kk]+= (1.0-0.5*alpha2_zz[kk])*Rz[kk];
                 Nz0[kk]+=      0.5*alpha3_zz[kk] *Rz[kk+1];
                
                 kk=nkt-1;
                 Nz0[kk] =     total_source_function_z[kk]*dt/f_use;
                 Nz0[kk]+= (1.0-0.5*alpha2_zz[kk])*Rz[kk];
                 Nz0[kk]+=      0.5*alpha1_zz[kk] *Rz[kk-1];


                if(galdef.use_symmetry==0)
                 tridag_ext     (Nz1,Nz2,Nz3,Nz0,Rz, particle.n_zgrid,nk );

                if(galdef.use_symmetry==1){
                 //#pragma vdir nodep
                 for(kk=0,kk1=1;kk<nkt-1;kk+=particle.n_zgrid,kk1+=particle.n_zgrid)Nz0[kk]+=0.5*alpha1_z0[kk] *Rz[kk1];

                 tridag_sym_ext (Nz1,Nz2,Nz3,Nz0,Rz, particle.n_zgrid,nk );
	                               	}


	      }//irept

 
		  /*
                  if(galdef.use_symmetry==2)
                  {
                     for(ip=0; ip<particle.n_pgrid; ip++)
                     {
                        for(ix=0; ix<=n_xgrid_sym; ix++)
                        {
                           for(iy=0; iy<=n_ygrid_sym; iy++)
                           {
 
                              for(iz=0; iz<particle.n_zgrid; iz++)
                              {
                                 double value=
                                    particle.cr_density.d3[                   ix][                   iy][iz].s[ip];
                                    particle.cr_density.d3[particle.n_xgrid-1-ix][                   iy][iz].s[ip]=value;
                                    particle.cr_density.d3[                   ix][particle.n_ygrid-1-iy][iz].s[ip]=value;
                                    particle.cr_density.d3[particle.n_xgrid-1-ix][particle.n_ygrid-1-iy][iz].s[ip]=value;
                              }  //  iz
                           }  //  iy 
                        }  //  ix
                     }  //  ip
		          }  //  symmetry
		  */

               }  //  prop_z

//cout<<"  z complete"<<endl;


	      // permute to p-first
                ns=particle.n_zgrid;
                kk=0;
                //#pragma vdir nodep
                for(i=0;i<ns;i++)
                //#pragma vdir nodep
                 for(k=i;k<nkt;k+=ns,kk++)Rp[kk]=Rz[k];
                kk=0;
                //#pragma vdir nodep
                for(i=0;i<ns;i++)
                 //#pragma vdir nodep
                 for(k=i;k<nkt;k+=ns,kk++)total_source_function_p[kk]=total_source_function_z[k]; 

// P propagation

	      if(galdef.prop_p==1){
		/*
                              kk=0;   
                              for(iz=0; iz<particle.n_zgrid; iz++)
                              {
                                  for(iy=0; iy<particle.n_ygrid; iy++)
                                  {
                                     for(ix=0; ix<particle.n_xgrid; ix++)
                                     {
                                          
				       for(ip=0; ip<particle.n_pgrid; ip++,           kk++){

                                 total_source_function_[kk]=  total_source_function.d3[ix][iy][iz].s[ip]; }
			  
 
                        }  //  ix
                     }  //  iy
                  }  //  iz
		*/

               nk=particle.n_xgrid*particle.n_ygrid*particle.n_zgrid;
 

              for(irept=0;irept<nrept;irept++){//cout<<irept<<endl;
               
                
                 //#pragma vdir nodep
                 for(kk=1;kk<nkt-1;kk++)Np0[kk] =     total_source_function_p[kk]*dt/f_use;
                 //#pragma vdir nodep
                 for(kk=1;kk<nkt-1;kk++)Np0[kk]+= (1.0-0.5*alpha2_pp[kk])*Rp[kk];
                 //#pragma vdir nodep
                 for(kk=1;kk<nkt-1;kk++)Np0[kk]+=      0.5*alpha1_pp[kk] *Rp[kk-1];
                 //#pragma vdir nodep
                 for(kk=1,kk1=2;kk<nkt-1;kk++,kk1++)Np0[kk]+=      0.5*alpha3_pp[kk] *Rp[kk1]; 
        
                             
                                    
               
                 kk=0;
                 Np0[kk] =     total_source_function_p[kk]*dt/f_use;
                 Np0[kk]+= (1.0-0.5*alpha2_pp[kk])*Rp[kk];
                 Np0[kk]+=      0.5*alpha3_pp[kk] *Rp[kk+1];
                
                 kk=nkt-1;
                 Np0[kk] =     total_source_function_p[kk]*dt/f_use;
                 Np0[kk]+= (1.0-0.5*alpha2_pp[kk])*Rp[kk];
                 Np0[kk]+=      0.5*alpha1_pp[kk] *Rp[kk-1];


                 tridag_ext (Np1,Np2,Np3,Np0,Rp, particle.n_pgrid,nk );

	       }// irept

 
		  /*
                  if(galdef.use_symmetry==2)
                  {
                     for(ip=0; ip<particle.n_pgrid; ip++)
                     {
                        for(ix=0; ix<=n_xgrid_sym; ix++)
                        {
                           for(iy=0; iy<=n_ygrid_sym; iy++)
                           {
                              for(iz=0; iz<=n_zgrid_sym; iz++)
                              {
 
                                 double value=
                                    particle.cr_density.d3[                   ix][                   iy][                   iz].s[ip];
                                    particle.cr_density.d3[                   ix][particle.n_ygrid-1-iy][                   iz].s[ip]=value;
                                    particle.cr_density.d3[                   ix][                   iy][particle.n_zgrid-1-iz].s[ip]=value;
                                    particle.cr_density.d3[                   ix][particle.n_ygrid-1-iy][particle.n_zgrid-1-iz].s[ip]=value;
                                    particle.cr_density.d3[particle.n_xgrid-1-ix][                   iy][                   iz].s[ip]=value;
                                    particle.cr_density.d3[particle.n_xgrid-1-ix][                   iy][particle.n_zgrid-1-iz].s[ip]=value;
                                    particle.cr_density.d3[particle.n_xgrid-1-ix][particle.n_ygrid-1-iy][                   iz].s[ip]=value;
                                    particle.cr_density.d3[particle.n_xgrid-1-ix][particle.n_ygrid-1-iy][particle.n_zgrid-1-iz].s[ip]=value; 
                              }  //  iz   
                           }  //  iy      
                        }  //  ix
                     }  //  ip
                  }  //  symmetry

*/
               }  //  prop_p	



	      // permute to x-first
                ns=particle.n_pgrid;
                kk=0;
                //#pragma vdir nodep
                for(i=0;i<ns;i++)
                 //#pragma vdir nodep
                 for(k=i;k<nkt;k+=ns,kk++)Rx[kk]=Rp[k];	  

}//irept_outer


                              kk=0;   
                              for(iz=0; iz<particle.n_zgrid; iz++)
                              {
                                  for(iy=0; iy<particle.n_ygrid; iy++)
                                  {
                                     for(ix=0; ix<particle.n_xgrid; ix++)
                                     {
                                          
                              for(ip=0; ip<particle.n_pgrid; ip++,           kk++)
                                 particle.cr_density.d3[ix][iy][iz].s[ip]=Rp[kk];
			  
 
                        }  //  ix
                     }  //  iy
                  }  //  iz
			//Gulli20070810
       delete[] Nx1; delete[] Nx2; delete[] Nx3; delete[] Nx0; delete[] Rx;       
       delete[] Ny1; delete[] Ny2; delete[] Ny3; delete[] Ny0; delete[] Ry;      
       delete[] Nz1; delete[] Nz2; delete[] Nz3; delete[] Nz0; delete[] Rz;       
       delete[] Np1; delete[] Np2; delete[] Np3; delete[] Np0; delete[] Rp;       

       delete[] alpha1_xx;  delete[] alpha1_yy;  delete[] alpha1_zz; delete[] alpha1_pp;       
       delete[] alpha2_xx;  delete[] alpha2_yy;  delete[] alpha2_zz; delete[] alpha2_pp;       
       delete[] alpha3_xx;  delete[] alpha3_yy;  delete[] alpha3_zz; delete[] alpha3_pp;         

       delete[] alpha1_x0;  delete[] alpha1_y0;  delete[] alpha1_z0;         

       delete[] total_source_function_x;       
       delete[] total_source_function_y;       
       delete[] total_source_function_z;       
       delete[] total_source_function_p;       

cout<<"<<protri"<<endl;
}
