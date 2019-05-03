
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * D_xx.cc *                                     galprop package * 02/13/2003 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// Wave damping formalism is described in:
//
// Ptuskin, V.S., et al. 2006, ApJ 642, 902
// Ptuskin, V.S., et al. 2005, Adv. Space Res. 35, 162
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

using namespace std;//AWS20050624
#include<cstdio>
#include<cstdlib>
#include"galprop_classes.h"
#include"galprop_internal.h"

int iprotons,ir,ix,iy,iz,ip;
int damping_min_ip; // IMOS20060330

//this is to avoid problems of using Galprop class members in static function "fu" IMOS20060322
Particle *protons;
double damping_p0;
int n_spatial_dimensions, diff_reacc;
#pragma omp threadprivate(iprotons,ir,ix,iy,iz,ip,damping_min_ip,protons,damping_p0,n_spatial_dimensions,diff_reacc)

int Galprop::D_xx(Particle &particle,int iprotons_,int ir_,int ix_,int iy_,int iz_,int ip_)
{
   iprotons=iprotons_; ir=ir_; ix=ix_; iy=iy_; iz=iz_; ip=ip_;
   double L_cm, Lp_cm, tmp;
// integration parameters
   double a=particle.rigidity[ip],ai;

//this is to avoid problems of using Galprop class members in static function "fu" IMOS20060322
   protons=&gcr[iprotons];
   damping_p0=galdef.damping_p0;
   n_spatial_dimensions=galdef.n_spatial_dimensions;
   diff_reacc=galdef.diff_reacc;


// STANDARD DIFFUSION COEFFICIENT (galdef.diff_reacc =0, 1, 2, -n==beta^n Dxx)
   if(galdef.diff_reacc < 3)
   {
// test of electron propagation vs analytical calculations IMOS20061030
     if(abs(galdef.DM_int0)==99) 
       particle.Dxx.d2[ir][iz].s[ip]=galdef.D0_xx *pow(particle.Ekin[ip]/galdef.D_rigid_br, galdef.D_g_1);
       
// end of the test area
     else //IMOS20070110
       {
	 if(n_spatial_dimensions==2)
	   {
	     particle.Dxx.d2[ir][iz].s[ip] = particle.beta[ip] *galdef.D0_xx;
	     if(galdef.diff_reacc<0) particle.Dxx.d2[ir][iz].s[ip] = pow(particle.beta[ip],galdef.diff_reacc) *galdef.D0_xx;
	     if(particle.rigidity[ip]< galdef.D_rigid_br)
	       particle.Dxx.d2[ir][iz].s[ip]*= pow(particle.rigidity[ip]/galdef.D_rigid_br, galdef.D_g_1);
	     if(particle.rigidity[ip]>=galdef.D_rigid_br)
	       particle.Dxx.d2[ir][iz].s[ip]*= pow(particle.rigidity[ip]/galdef.D_rigid_br, galdef.D_g_2);
	   }
	 if(n_spatial_dimensions==3)
	   {
	     //particle.Dxx.d3[ix][iy][iz].s[ip] = particle.beta[ip] *galdef.D0_xx;
             if (abs(galaxy.z[iz]) > galdef.z_step2) { 
              particle.Dxx.d3[ix][iy][iz].s[ip] =  particle.beta[ip]*galdef.D0_xx*galdef.D_fact*2.;//SB20170323	
	      if(galdef.diff_reacc<0) particle.Dxx.d3[ix][iy][iz].s[ip] = pow(particle.beta[ip],galdef.diff_reacc) *galdef.D0_xx*galdef.D_fact*2.;
	      cout<<"Now Z2 and D1"<<galaxy.z[iz]<<""<<particle.Dxx.d3[ix][iy][iz].s[ip]<<endl; 
	      //cout<<"Now Diff:"<<galdef.D0_xx<<endl;
             }else {	
             
	       particle.Dxx.d3[ix][iy][iz].s[ip] =  particle.beta[ip]*galdef.D0_xx;
 		
	       if(galdef.diff_reacc<0) particle.Dxx.d3[ix][iy][iz].s[ip] = pow(particle.beta[ip],galdef.diff_reacc) *galdef.D0_xx;
               //cout<<"Zafter and Diff"<<galaxy.z[iz]<<galdef.D0_xx<<endl;
             }	
	     if (abs(galaxy.z[iz]) > galdef.z_step2) {
	      particle.Dxx.d3[ix][iy][iz].s[ip] =  particle.beta[ip]*galdef.D0_xx*galdef.D_fact;
	      if(galdef.diff_reacc<0) particle.Dxx.d3[ix][iy][iz].s[ip] = pow(particle.beta[ip],galdef.diff_reacc) *galdef.D0_xx*galdef.D_fact;
	      cout<<"Now Z2 and D2"<<galaxy.z[iz]<<""<<particle.Dxx.d3[ix][iy][iz].s[ip]<<endl; 		
	     }else {
	       particle.Dxx.d3[ix][iy][iz].s[ip] =  particle.beta[ip]*galdef.D0_xx;
	        if(galdef.diff_reacc<0) particle.Dxx.d3[ix][iy][iz].s[ip] = pow(particle.beta[ip],galdef.diff_reacc) *galdef.D0_xx;
               //cout<<"Zafter and Diff"<<galaxy.z[iz]<<galdef.D0_xx<<endl;
             }	 	
	       //if(galdef.diff_reacc<0) particle.Dxx.d3[ix][iy][iz].s[ip] = pow(particle.beta[ip],galdef.diff_reacc) *galdef.D0_xx;
	       if(particle.rigidity[ip]< galdef.D_rigid_br)
	         particle.Dxx.d3[ix][iy][iz].s[ip]*= pow(particle.rigidity[ip]/galdef.D_rigid_br, galdef.D_g_1);
	       if(particle.rigidity[ip]>=galdef.D_rigid_br)
	         particle.Dxx.d3[ix][iy][iz].s[ip]*= pow(particle.rigidity[ip]/galdef.D_rigid_br, galdef.D_g_2);
	     	
	   }
       //cout<<"check print D:" <<<<endl;
       }
     return 0;
   }

// WAVE DAMPING (see Ptuskin et al. astro-ph/0301420)

   if(ip==particle.n_pgrid-1) damping_min_ip=0; // IMOS20060330
   L_cm  = galdef.damping_max_path_L;                               // max free path
   if(ip<damping_min_ip) // IMOS20060330
     {
       if(n_spatial_dimensions==2) particle.Dxx.d2[ir]    [iz].s[ip] = particle.beta[ip] *C*L_cm/3.;
       if(n_spatial_dimensions==3) particle.Dxx.d3[ix][iy][iz].s[ip] = particle.beta[ip] *C*L_cm/3.;
// printout at the solar system position
       if(iz==particle.n_zgrid/2+1 && ir==9) cout<<" D_xx>>>> "<<particle.rigidity[ip]<<" "<<particle.Ekin[ip]<<" "<<particle.beta[ip] *C*L_cm/3.<<" "<<-ai<<endl;
       return 0;
     }
   Lp_cm = 3./C *galdef.D0_xx*pow(particle.rigidity[ip]/galdef.D_rigid_br,galdef.D_g_1);

   for(int i=1;i<gcr[iprotons].n_pgrid;i++) 
      if(gcr[iprotons].rigidity[i]> galdef.damping_p0)
      {
	 damping_p0 = gcr[iprotons].rigidity[i];  // re-definition of galdef.damping_p0
         break;
      }
/*******************TEST
ir=0; iz=0;
for(int i=0; i<particle.n_pgrid; i++)
{
gcr[iprotons].cr_density.d2[0][0].s[i] = pow(gcr[iprotons].p[i],-2.);
cout<<" p= "<<gcr[iprotons].p[i]<<"  y= "<<gcr[iprotons].cr_density.d2[0][0].s[i]<<endl;
}
double xxx=0.001;
fu(xxx);
*******************/

//   static double (*fuPtr)(double);
//   fuPtr = &Galprop::fu;

   ai = 0.;
   if(n_spatial_dimensions==2)
     if(iz>0 && iz<particle.n_zgrid-1 && ir<particle.n_rgrid-1)
       if(a<damping_p0) ai=sim(damping_p0,a,a/100.,0.01,1.e-30,&Galprop::fu);  //IMOS20060330
   if(n_spatial_dimensions==3)
     if(iz>0 && iz<particle.n_zgrid-1 && ix>0 && ix<particle.n_xgrid-1 && iy>0 && iy<particle.n_ygrid-1)
       if(a<damping_p0) ai=sim(damping_p0,a,a/100.,0.01,1.e-10,&Galprop::fu);  //IMOS20060330
//   if(a<damping_p0) ai=sim(damping_p0,a,h,reps,aeps,&fu);  //IMOS20060330
//   cout<<" Dxx------>"<<damping_p0<<" "<<a<<" "<<ai<<endl;

// Kolmogorov diffusion with wave damping ## 
   if(galdef.diff_reacc==11)
   {
      tmp = 1. -galdef.damping_const_G*(-ai); //    /pow(particle.Z, 5./3.)
      if(tmp>0.) L_cm = Lp_cm/pow(tmp, 2);
   }

// Kraichnan diffusion with wave damping ## 
   if(galdef.diff_reacc==12)
   {
      tmp = 1. -galdef.damping_const_G*(-ai); //    /pow(particle.Z, 3./2.)
      if(tmp>0.) L_cm = Lp_cm/tmp;
   }

   if (L_cm>galdef.damping_max_path_L && gcr[iprotons].rigidity[ip]<1.e4) // IMOS20050907
     { 
       L_cm = galdef.damping_max_path_L;
       damping_min_ip=ip;
     }
   if(n_spatial_dimensions==2) particle.Dxx.d2[ir]    [iz].s[ip] = particle.beta[ip] *C*L_cm/3.;
   if(n_spatial_dimensions==3) particle.Dxx.d3[ix][iy][iz].s[ip] = particle.beta[ip] *C*L_cm/3.;
   if(iz==particle.n_zgrid/2+1 && ir==9) 
     cout<<" D_xx>>>> "<<particle.rigidity[ip]<<" "<<particle.Ekin[ip]
	 <<" "<<particle.beta[ip] *C*L_cm/3.<<" "<<-ai<<endl;
   return 0;
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

double Galprop::fu(double x)
{
   int i,n,m;
   double int_psi=0., y;

//     cout<<" ir,iz,ip,iprotons,x= "<<ir<<" "<<iz<<" "<<ip<<" "<<iprotons<<" "<<x<<endl;
//     for(int i=0; i<protons->n_pgrid; i++) cout<<" "<<protons->cr_density.d2[ir][iz-1].s[i];
//     cout<<endl;

// search in the grid
   for(n=1;n<protons->n_pgrid;n++) if(protons->p[n]> x) break;
   for(m=1;m<protons->n_pgrid;m++) if(protons->p[m]> damping_p0) break;
//   cout<<">>>> x,n= "<<x<<" "<<n<<" gcr[iprotons].p[n]= "<<gcr[iprotons].p[n]<<endl;

// integration over rigidity (= momentum for protons)

   if(n_spatial_dimensions==2)
   { 
// fit each interval with power-law and integrate analytically
      for(int_psi=0., i=n-1; i<m; i++)
      {
	 if ( protons->cr_density.d2[ir]    [iz].s[i] > 0 && protons->cr_density.d2[ir]    [iz].s[i+1] > 0 )
	 {
// derive y (= power-law index) 
         y =log(protons->cr_density.d2[ir]    [iz].s[i]/protons->cr_density.d2[ir]    [iz].s[i+1])
           /log(protons->                          p[i]/protons->                          p[i+1]);
// integrate (psi/p) dp
         if(i>n-1) int_psi +=protons->cr_density.d2[ir][iz].s[i  ]       // fit norm.
		        *pow(protons-> p[i  ],-y)/y
                       *(pow(protons-> p[i+1], y)
                        -pow(protons-> p[i  ], y));
         else      int_psi +=protons->cr_density.d2[ir][iz].s[i  ]       // fit norm.
                        *pow(protons-> p[i  ],-y)/y
                       *(pow(protons-> p[i+1], y)
                        -pow(x,               y));
	 }
      }
   }

   if(n_spatial_dimensions==3)
   { 
// fit each interval with power-law and integrate analytically
      for(int_psi=0., i=n-1; i<m; i++)
      {
// derive y (= power-law index) 
	 if ( protons->cr_density.d3[ix][iy][iz].s[i] > 0 && protons->cr_density.d3[ix][iy][iz].s[i+1] > 0 )
	 {
         y =log(protons->cr_density.d3[ix][iy][iz].s[i]/protons->cr_density.d3[ix][iy][iz].s[i+1])
           /log(protons->                          p[i]/protons->                          p[i+1]);
// integrate (psi/p) dp
         if(i>n-1) int_psi +=protons->cr_density.d3[ix][iy][iz].s[i  ]       // fit norm.
                        *pow(protons-> p[i  ],-y)/y
                       *(pow(protons-> p[i+1], y)
                        -pow(protons-> p[i  ], y));
         else      int_psi +=protons->cr_density.d3[ix][iy][iz].s[i  ]       // fit norm.
                        *pow(protons-> p[i  ],-y)/y
                       *(pow(protons-> p[i+1], y)
                        -pow(x,               y));
	 }
      }
   }

//if(ir==0 && iz==80 && ip==40) cout<<" integral = "<<n<<" "<<y<<" "<<int_psi<<" "<<protons->p[n]<<" "<<x<<" "<<damping_p0 <<endl; //exit(1);
   if(diff_reacc==11) return ( pow(x,2./3.)*int_psi );   // Kolmogorov
   if(diff_reacc==12) return (sqrt(x)      *int_psi );   // Kraichnan
   return 0; // in case of an error return 0.
}
