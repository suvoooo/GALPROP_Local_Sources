
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_synch_skymap.cc *                         galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// generate  synchrotron skymap 

using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"

int Galprop::gen_synch_skymap()
{
   int stat=0;
   cout<<" >>>>gen_synch_skymap"<<endl;
   
   if (galdef.skymap_format == 3 || 4 == galdef.skymap_format){
#pragma omp parallel for default(shared)
     for (int ip = 0; ip < galaxy.synchrotron_hp_skymap.Npix(); ++ip){
       SM::Coordinate co(galaxy.synchrotron_hp_skymap.pix2ang(ip));
       double l=co.l();
       double b=co.b();
       vector<double> synch, Q, U; //AWS20100709
     
       //     gen_synch_skymap_pixel(l, b,synch);
       gen_synch_skymap_pixel(l, b,synch,Q,U); //AWS20100709

       for(int inusynch=0; inusynch<galaxy.n_nu_synchgrid; inusynch++)
	 {
	 galaxy.synchrotron_hp_skymap  [co][inusynch] = synch[inusynch];
	 galaxy.synchrotron_Q_hp_skymap[co][inusynch] = Q    [inusynch]; //AWS20100709
	 galaxy.synchrotron_U_hp_skymap[co][inusynch] = U    [inusynch]; //AWS20100709
	 }
     }

     galaxy.synchrotron_hp_skymap  .setSpectra(galaxy.nu_synch,galaxy.n_nu_synchgrid);
     galaxy.synchrotron_Q_hp_skymap.setSpectra(galaxy.nu_synch,galaxy.n_nu_synchgrid); //AWS20100709
     galaxy.synchrotron_U_hp_skymap.setSpectra(galaxy.nu_synch,galaxy.n_nu_synchgrid); //AWS20100709


     if(galdef.verbose>=2)
       {
	 cout<<" synchrotron skymap " <<endl;
	 galaxy.synchrotron_hp_skymap.print(cout);
       }//galdef.verbose>=2
   }else{
#pragma omp parallel for schedule(dynamic) default(shared)
     for(int i_long=0; i_long<galaxy.n_long; i_long++)
       {
	 for(int i_lat =0; i_lat <galaxy.n_lat; i_lat++)
	   {
	     double l=galaxy.long_min + i_long*galaxy.d_long;
	     double b=galaxy. lat_min + i_lat *galaxy.d_lat ;
	     vector<double> synch,Q,U;

	     //     gen_synch_skymap_pixel(l, b,synch);
             gen_synch_skymap_pixel(l, b,synch,Q,U); //AWS20100709

	     for(int inusynch=0; inusynch<galaxy.n_nu_synchgrid; inusynch++)
	     {
	       galaxy.synchrotron_skymap  .d2[i_long][i_lat].s[inusynch] = synch[inusynch];
	       galaxy.synchrotron_Q_skymap.d2[i_long][i_lat].s[inusynch] = Q    [inusynch]; //AWS20100709
	       galaxy.synchrotron_U_skymap.d2[i_long][i_lat].s[inusynch] = U    [inusynch]; //AWS20100709
	     }
	   }
       }
     if(galdef.verbose>=2)
       {
	 cout<<" synchrotron skymap " <<endl;
	 galaxy.synchrotron_skymap.print();
       }//galdef.verbose>=2
   }
   cout<<" <<<< gen_synch_skymap"<<endl;
   return stat;
}
///////////////////////////////////////////////////////////////////////////////////////////
int Galprop::gen_synch_skymap_pixel(const double l, const double b, vector<double> &synch){
  double Ro= 8.5; // Galactocentric distance of Sun, kpc
  Ro= 8.3; // to avoid discontinuity in maps
  double dd0    =0.1; // max integration step in kpc at b=90 deg.
  double ddmax = 0.5; // max integration step allowed
  double dtr=acos(-1.)/180.;
  int ir,ix,iy,iz;
  //cout<<"l b ="<<l<<" "<<b<<endl;
  double sinb=sin(b*dtr);
  double cosb=cos(b*dtr);
  double sinl=sin(l*dtr);
  double cosl=cos(l*dtr);
  double d=0;
  double dd=dd0/(fabs(sinb)+1.e-6);
  if(dd>ddmax)dd=ddmax;
  int complete=0;
  synch.resize(galaxy.n_nu_synchgrid, 0);
  while(complete==0)
    {
      d += dd;
      double RR=sqrt(Ro*Ro+pow(d*cosb,2)-2.0*Ro*d*cosb*cosl);
      double zz=d*sinb;
      double costheta=(Ro-d*cosb*cosl)/RR;
      if(costheta> 1.0)costheta= 1.0;
      if(costheta<-1.0)costheta=-1.0;
      
      if(gcr[0].n_spatial_dimensions==2)
	{
	  ir=(int)((RR-galaxy.r_min)/galaxy.dr + 0.5);//IMOS20060420
	  if(ir>galaxy.n_rgrid-1) { complete=1; ir=galaxy.n_rgrid-1; }
	  iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.5);//IMOS20060420
	  if(iz<0               ) { complete=1; iz=0;                } 
	  if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }
	  // cout<<"d RR zz ir iz "<<d<<" "<<RR<<" "<<zz<<" "<<ir<<" "<<iz<<endl;
	} // particle.n_spatial_dimensions==2
      
      if(gcr[0].n_spatial_dimensions==3)
	{
	  int LH=1; // 1=LH, 2=RH system                                                     AWS20080311
	  
	  double xx=Ro-d*cosb*cosl; // Sun on x axis at x=+Ro
	  double yy;
	  
	  if(LH==1)                                                                       // AWS20080311
	    yy=  +d*cosb*sinl; // Sun at y=0; +ve long=+ve y for Z=-X^Y (LH) system  // AWS20080311
               if(LH==2)                                                                       // AWS20080311
		 yy=  -d*cosb*sinl; // Sun at y=0; +ve long=-ve y for Z= X^Y (RH) system  // AWS20080311
	       
	       
               if(galdef.use_symmetry==1)
		 {
                  xx=fabs(xx);
                  yy=fabs(yy);
                  zz=fabs(zz);
               }
               ix=(int)((xx-galaxy.x_min)/galaxy.dx + 0.5);//IMOS20060420
               iy=(int)((yy-galaxy.y_min)/galaxy.dy + 0.5);//IMOS20060420
               iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.5);//IMOS20060420

               if(ix<0               ) { complete=1; ix=0;                }
               if(iy<0               ) { complete=1; iy=0;                }  
               if(iz<0               ) { complete=1; iz=0;                } 
               if(ix>galaxy.n_xgrid-1) { complete=1; ix=galaxy.n_xgrid-1; }
               if(iy>galaxy.n_ygrid-1) { complete=1; iy=galaxy.n_ygrid-1; }
               if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }
//  cout<<"d RR xx yy zz ix iy iz "<<d<<" "<<RR<<" "<<xx<<" "<<yy<<" "<<zz<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
            } // particle.n_spatial_dimensions==3

            for(int inusynch=0; inusynch<galaxy.n_nu_synchgrid; inusynch++)
            {
              // float delta;   AWS20100810
	        double delta; //AWS20100810
               if(gcr[0].n_spatial_dimensions==2)
                  delta =dd*kpc2cm *galaxy.synchrotron_emiss.d2[ir][iz].    s[inusynch];

               if(gcr[0].n_spatial_dimensions==3)
		  delta =dd*kpc2cm *galaxy.synchrotron_emiss.d3[ix][iy][iz].s[inusynch];

               synch[inusynch] += delta;

//cout<<"l b ir iz  E_gamma bremss_ionized_emiss "<<l<<" "<<b<<" "<<ir<<" "<<iz<<" "<<" "
//<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.bremss_ionized_emiss.d2[ir][iz].s[iEgamma]<<endl;     
            }//inusynch
         }//complete==0
	return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////// Stokes version, AWS20100709
int Galprop::gen_synch_skymap_pixel(const double l, const double b, vector<double> &synch, vector<double> &Q,  vector<double> &U )
{
  double Ro= 8.5; // Galactocentric distance of Sun, kpc
  Ro= 8.3; // to avoid discontinuity in maps
  double dd0    =0.1; // max integration step in kpc at b=90 deg.
  double ddmax = 0.5; // max integration step allowed
  double dtr=acos(-1.)/180.;
  int ir,ix,iy,iz;
  //cout<<"l b ="<<l<<" "<<b<<endl;
  double sinb=sin(b*dtr);
  double cosb=cos(b*dtr);
  double sinl=sin(l*dtr);
  double cosl=cos(l*dtr);
  double d=0;
  double dd=dd0/(fabs(sinb)+1.e-6);
  if(dd>ddmax)dd=ddmax;
  int complete=0;
  synch.resize(galaxy.n_nu_synchgrid, 0);
  Q    .resize(galaxy.n_nu_synchgrid, 0);
  U    .resize(galaxy.n_nu_synchgrid, 0);

  while(complete==0)
    {
      d += dd;
      double RR=sqrt(Ro*Ro+pow(d*cosb,2)-2.0*Ro*d*cosb*cosl);
      double zz=d*sinb;
      double costheta=(Ro-d*cosb*cosl)/RR;
      if(costheta> 1.0)costheta= 1.0;
      if(costheta<-1.0)costheta=-1.0;
      
      if(gcr[0].n_spatial_dimensions==2)
	{
	  ir=(int)((RR-galaxy.r_min)/galaxy.dr + 0.5);//IMOS20060420
	  if(ir>galaxy.n_rgrid-1) { complete=1; ir=galaxy.n_rgrid-1; }
	  iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.5);//IMOS20060420
	  if(iz<0               ) { complete=1; iz=0;                } 
	  if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }
	  // cout<<"d RR zz ir iz "<<d<<" "<<RR<<" "<<zz<<" "<<ir<<" "<<iz<<endl;
	} // particle.n_spatial_dimensions==2
      
      if(gcr[0].n_spatial_dimensions==3)
	{
	  int LH=1; // 1=LH, 2=RH system                                                     AWS20080311
	  
	  double xx=Ro-d*cosb*cosl; // Sun on x axis at x=+Ro
	  double yy;
	  
	  if(LH==1)                                                                       // AWS20080311
	    yy=  +d*cosb*sinl; // Sun at y=0; +ve long=+ve y for Z=-X^Y (LH) system  // AWS20080311
               if(LH==2)                                                                       // AWS20080311
		 yy=  -d*cosb*sinl; // Sun at y=0; +ve long=-ve y for Z= X^Y (RH) system  // AWS20080311
	       
	       
               if(galdef.use_symmetry==1)
		 {
                  xx=fabs(xx);
                  yy=fabs(yy);
                  zz=fabs(zz);
               }
               ix=(int)((xx-galaxy.x_min)/galaxy.dx + 0.5);//IMOS20060420
               iy=(int)((yy-galaxy.y_min)/galaxy.dy + 0.5);//IMOS20060420
               iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.5);//IMOS20060420

               if(ix<0               ) { complete=1; ix=0;                }
               if(iy<0               ) { complete=1; iy=0;                }  
               if(iz<0               ) { complete=1; iz=0;                } 
               if(ix>galaxy.n_xgrid-1) { complete=1; ix=galaxy.n_xgrid-1; }
               if(iy>galaxy.n_ygrid-1) { complete=1; iy=galaxy.n_ygrid-1; }
               if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }
//  cout<<"d RR xx yy zz ix iy iz "<<d<<" "<<RR<<" "<<xx<<" "<<yy<<" "<<zz<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
            } // particle.n_spatial_dimensions==3

            for(int inusynch=0; inusynch<galaxy.n_nu_synchgrid; inusynch++)
            {
	    //float  delta,deltaQ,deltaU;    AWS20100810
	      double delta,deltaQ,deltaU; // AWS10100810
               if(gcr[0].n_spatial_dimensions==2)
	       {
                  delta =dd*kpc2cm *galaxy.synchrotron_emiss  .d2[ir][iz].    s[inusynch];
                  deltaQ=dd*kpc2cm *galaxy.synchrotron_Q_emiss.d2[ir][iz].    s[inusynch];
                  deltaU=dd*kpc2cm *galaxy.synchrotron_U_emiss.d2[ir][iz].    s[inusynch];
	       }
               if(gcr[0].n_spatial_dimensions==3)
	       {
		  delta =dd*kpc2cm *galaxy.synchrotron_emiss  .d3[ix][iy][iz].s[inusynch];
		  deltaQ=dd*kpc2cm *galaxy.synchrotron_Q_emiss.d3[ix][iy][iz].s[inusynch];
		  deltaU=dd*kpc2cm *galaxy.synchrotron_U_emiss.d3[ix][iy][iz].s[inusynch];
	       }

               synch[inusynch] += delta;
               Q    [inusynch] += deltaQ;
               U    [inusynch] += deltaU;

//cout<<"l b ir iz  E_gamma bremss_ionized_emiss "<<l<<" "<<b<<" "<<ir<<" "<<iz<<" "<<" "
//<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.bremss_ionized_emiss.d2[ir][iz].s[iEgamma]<<endl; 
    
            }//inusynch
         }//complete==0

	return 0;
}
