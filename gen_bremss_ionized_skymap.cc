
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_bremss_ionized_skymap.cc *                galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// generate skymap for bremsstrahlung on ionized hydrogen 

#include <cassert>
#include <string>
#include <cstring>
#include <vector>

using namespace std;//AWS20050624

#include "galprop_classes.h"
#include "galprop_internal.h"

#include <ErrorLogger.h>

int Galprop::gen_bremss_ionized_skymap() {

  INFO("Entry");

  int stat=0;
  
  if (3 == galdef.skymap_format || 4 == galdef.skymap_format){
#pragma omp parallel for schedule(dynamic) default(shared)
    for (int ip = 0; ip < galaxy.bremss_ionized_hp_skymap.Npix(); ++ip){
      SM::Coordinate co(galaxy.bremss_ionized_hp_skymap.pix2ang(ip));
      double l=co.l();
      double b=co.b();
      vector<double> bremss_ionized;
      gen_bremss_ionized_skymap_pixel(l,b,bremss_ionized);
      for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++){
	galaxy.bremss_ionized_hp_skymap[co][iEgamma]=bremss_ionized[iEgamma];
      }
    }
    galaxy.bremss_ionized_hp_skymap.setSpectra(galaxy.E_gamma,galaxy.n_E_gammagrid);
    if(galdef.verbose>=2)
      {
        cout<<" bremsstrahlung on ionized hydrogen skymap " <<endl;
        galaxy.bremss_ionized_hp_skymap.print(cout);
      } // galdef.verbose>=2
  }else{
#pragma omp parallel for schedule(dynamic) default(shared)
    for(int i_long=0; i_long<galaxy.n_long; i_long++)
      {
        for(int i_lat =0; i_lat<galaxy.n_lat; i_lat++)
	  {
	    double l =galaxy.long_min + i_long*galaxy.d_long;
	    double b =galaxy. lat_min + i_lat *galaxy.d_lat;
	    vector<double> bremss_ionized;
	    gen_bremss_ionized_skymap_pixel(l,b,bremss_ionized);
	    for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++){
	      galaxy.bremss_ionized_skymap.d2[i_long][i_lat].s[iEgamma]=bremss_ionized[iEgamma];
	    }
	  }
      }
    
    if(galdef.verbose>=2)
      {
        cout<<" bremsstrahlung on ionized hydrogen skymap " <<endl;
        galaxy.bremss_ionized_skymap.print();
      } // galdef.verbose>=2
  }
  INFO("Exit");
  return stat;
}

int Galprop::gen_bremss_ionized_skymap_pixel(const double l, const double b, vector<double> &bremss_ionized){
  double dd0    =0.01; // max integration step in kpc at b=90 deg.
  double ddmax = 0.01; // max integration step allowed
  double dtr=acos(-1.)/180.;
  int ir,ix,iy,iz;
  if(galdef.verbose>=1) cout<<" gen_bremss_ionized_skymap l b ="<<l<<" "<<b<<endl;
  double sinb=sin(b*dtr);
  double cosb=cos(b*dtr);
  double sinl=sin(l*dtr);
  double cosl=cos(l*dtr);
  double d=0;
  double dd=dd0/(fabs(sinb)+1.e-6);
  if(dd>ddmax)dd=ddmax;
  int complete=0;
  bremss_ionized.resize(galaxy.n_E_gammagrid,0);
  while(complete==0)
    {
      d += dd;
      double RR=sqrt(Rsun*Rsun+pow(d*cosb,2)-2.0*Rsun*d*cosb*cosl); //IMOS20080114
      double zz=d*sinb;
      double costheta=(Rsun-d*cosb*cosl)/RR; //IMOS20080114
      if(costheta> 1.0) costheta= 1.0;
      if(costheta<-1.0) costheta=-1.0;
      
      if(gcr[0].n_spatial_dimensions==2)
	{
	  ir=(int)((RR-galaxy.r_min)/galaxy.dr + 0.5);//IMOS20060420
	  if(ir>galaxy.n_rgrid-1) { complete=1; ir=galaxy.n_rgrid-1; }
	  iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.5);//IMOS20060420
	  if(iz<0               ) { complete=1; iz=0;                } 
	  if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }
	  // cout<<"d RR zz ir iz "<<d<<" "<<RR<<" "<<zz<<" "<<ir<<" "<<iz<<endl;
	}
      if(gcr[0].n_spatial_dimensions==3)
	{
	  double xx=Rsun-d*cosb*cosl; // Sun on x axis at x=+Ro
	  double yy=    -d*cosb*sinl; // Sun at y=0; +ve long=-ve y since Z=X^Y system
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
	}
      
      for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
	{
	  float delta;
	  if(gcr[0].n_spatial_dimensions==2) delta =dd*kpc2cm
					       *galaxy.bremss_ionized_emiss.d2[ir][iz].    s[iEgamma]
					       *galaxy.               n_HII.d2[ir][iz].    s[0];
	  
	  if(gcr[0].n_spatial_dimensions==3) delta =dd*kpc2cm
					       *galaxy.bremss_ionized_emiss.d3[ix][iy][iz].s[iEgamma]
					       *galaxy.               n_HII.d3[ix][iy][iz].s[0];
	  
	  bremss_ionized[iEgamma]+=delta;
	  
	  //cout<<"l b ir iz  E_gamma bremss_ionized_emiss "<<l<<" "<<b<<" "<<ir<<" "<<iz<<" "<<" "
	  //<<galaxy.E_gamma[iEgamma]<<" "<<galaxy.bremss_ionized_emiss.d2[ir][iz].s[iEgamma]<<endl;     
	}
    } // complete==0
  return 0;
}
