//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_DM_skymap.cc *                            galprop package * 9/14/2005
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// generate DM photons skymaps  IMOS20050912
// Calculates for galdef.z_min<=Z<=galdef.z_max only

using namespace std;
#include"galprop_classes.h"
#include"galprop_internal.h"

#include <ErrorLogger.h>

int Galprop::gen_DM_skymap() {

  INFO("Entry");
  
  int stat=0;
  
  if (galdef.skymap_format == 3 || 4 == galdef.skymap_format){
#pragma omp parallel for default(shared)
    for (int ip = 0; ip < galaxy.DM_hp_skymap.Npix(); ++ip){
      SM::Coordinate co(galaxy.DM_hp_skymap.pix2ang(ip));
      double l=co.l();
      double b=co.b();
      vector<double> DM;
      gen_DM_skymap_pixel(l, b,DM);
      for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
	galaxy.DM_hp_skymap[co][iEgamma] = DM[iEgamma];
    }
    galaxy.DM_hp_skymap.setSpectra(galaxy.E_gamma,galaxy.n_E_gammagrid);
    if(galdef.verbose>=2)
      {
        cout<<" DM skymap "<<endl;
        galaxy.DM_hp_skymap.print(cout);
      } // galdef.verbose>=2
  }else{
#pragma omp parallel for schedule(dynamic) default(shared)
    for(int i_long=0; i_long<galaxy.n_long; i_long++)
      {
        for(int i_lat =0; i_lat<galaxy.n_lat; i_lat++)
	  {
	    double l=galaxy.long_min + i_long*galaxy.d_long;
	    double b=galaxy. lat_min + i_lat *galaxy.d_lat ;
	    vector<double> DM;
	    gen_DM_skymap_pixel(l, b,DM);
	    for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
	      galaxy.DM_skymap.d2[i_long][i_lat].s[iEgamma] = DM[iEgamma];
	  }
      }
    if(galdef.verbose>=2)
      {
        cout<<" DM skymap "<<endl;
        galaxy.DM_skymap.print();
      } // galdef.verbose>=2
  }
  INFO("Exit");
  return stat;
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// calc. of a 2D array of g-ray emission for the given E_gammagrid for a
// particular pixel (l,b)   IMOS20080114
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
int Galprop::gen_DM_skymap_pixel(const double l, const double b, vector<double> &DM){
  double dtr=acos(-1.)/180.;
  int ir,ix,iy,iz;
  
  if(galdef.verbose>=1) cout<<"  gen_DM_skymap l b ="<<l<<" "<<b<<endl;

  double sinb=sin(b*dtr);
  double cosb=cos(b*dtr);
  double sinl=sin(l*dtr);
  double cosl=cos(l*dtr);
  double d=0;
  int complete=0;

  double dd=galdef.LoS_step /galdef.LoS_substep_number /(fabs(sinb)+1.e-6);// variable step depends on b
  if(dd>galdef.LoS_step) dd=galdef.LoS_step;                               // max integration step in kpc 
  
  DM.resize(galaxy.n_E_gammagrid, 0);

  while(complete==0)
    {
      d += dd;
      double zz=d*sinb;
      double RR=sqrt(Rsun*Rsun+pow(d*cosb,2)-2.0*Rsun*d*cosb*cosl); //IMOS20080114
      double xx,yy;                                                 //IMOS20080114
      
     // checks if we got to the Galactic boundary in 2D, if so stop integration
      if(gcr[0].n_spatial_dimensions==2)
	{
	  // find the nearest grid points on the LEFT side of the current point
	  ir=(int)((RR-galaxy.r_min)/galaxy.dr);
	  iz=(int)((zz-galaxy.z_min)/galaxy.dz);

	  // checks if we got to the Galactic boundary in 2D, if so stop integration
	  if(RR>galaxy.r_max)                    complete=1;
	  if(zz<galaxy.z_min || zz>galaxy.z_max) complete=1;

	  if(ir>galaxy.n_rgrid-1) {              complete=1; ir=galaxy.n_rgrid-1; }
	  
	  if(iz<0               ) {              complete=1; iz=0; } 
	  if(iz>galaxy.n_zgrid-1) {              complete=1; iz=galaxy.n_zgrid-1; }
	  // cout<<"d RR zz ir iz "<<d<<" "<<RR<<" "<<zz<<" "<<ir<<" "<<iz<<endl;

	} // particle.n_spatial_dimensions==2
      
       // checks if we got to the Galactic boundary in 3D, if so stop integration
      if(gcr[0].n_spatial_dimensions==3) 
	{
	  xx=Rsun-d*cosb*cosl;                                   // 3D: Sun on x axis at x=+Rsun
	  yy=    -d*cosb*sinl;                                   // 3D: Sun at y=0; +ve long=-ve y since Z=X^Y system

	  if(galdef.use_symmetry==1) 
	    {	    
	      xx=fabs(xx);
	      yy=fabs(yy);
	      zz=fabs(zz);
	    }
	  // find the nearest grid points on the LEFT side of the current point
	  ix=(int)((xx-galaxy.x_min)/galaxy.dx);
	  iy=(int)((yy-galaxy.y_min)/galaxy.dy);
	  iz=(int)((zz-galaxy.z_min)/galaxy.dz);

	  // checks if we got to the Galactic boundary in 3D, if so stop integration
	  if(ix<0               ) { complete=1; ix=0;                }
	  if(iy<0               ) { complete=1; iy=0;                }  
	  if(iz<0               ) { complete=1; iz=0;                } 
	  if(ix>galaxy.n_xgrid-1) { complete=1; ix=galaxy.n_xgrid-1; }
	  if(iy>galaxy.n_ygrid-1) { complete=1; iy=galaxy.n_ygrid-1; }
	  if(iz>galaxy.n_zgrid-1) { complete=1; iz=galaxy.n_zgrid-1; }
	  
	  if(zz<galaxy.z_min || zz>galaxy.z_max) complete=1;
//	  if(fabs(zz) > zzmax                  ) complete=1;
	  //  cout<<"d RR xx yy zz ix iy iz "<<d<<" "<<RR<<" "<<xx<<" "<<yy<<" "<<zz<<" "<<ix<<" "<<iy<<" "<<iz<<endl;

	} //particle.n_spatial_dimensions==3
      
      for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
	{
	  float delta, x[8][3],f[8],y[7];

	      if(gcr[0].n_spatial_dimensions==2) 
		{
		  if (ir==galaxy.n_rgrid-1 || iz==galaxy.n_zgrid-1)
		    delta = dd*kpc2cm *galaxy.DM_emiss.d2[ir][iz]    .s[iEgamma];

		  else  // linear interpolation
		    {
		      //  x[0]=(R0,z1), x[1]=(R1,z1);  y[0] -location of the grid points
		      //  x[2]=(R0,z0), x[3]=(R1,z0);  y[1]
		      x[0][0]=galaxy.r[ir  ]; x[0][1]=galaxy.z[iz+1];  f[0]=galaxy.DM_emiss.d2[ir  ][iz+1].s[iEgamma];
		      x[1][0]=galaxy.r[ir+1]; x[1][1]=galaxy.z[iz+1];  f[1]=galaxy.DM_emiss.d2[ir+1][iz+1].s[iEgamma];
		      x[2][0]=galaxy.r[ir  ]; x[2][1]=galaxy.z[iz  ];  f[2]=galaxy.DM_emiss.d2[ir  ][iz  ].s[iEgamma];
		      x[3][0]=galaxy.r[ir+1]; x[3][1]=galaxy.z[iz  ];  f[3]=galaxy.DM_emiss.d2[ir+1][iz  ].s[iEgamma];
		      
		      y[0]=  (f[0]-f[1])/(x[0][0]-x[1][0])*(RR-x[0][0])+f[0]; // interpolation in R
		      y[1]=  (f[2]-f[3])/(x[2][0]-x[3][0])*(RR-x[2][0])+f[2];
		      
		      y[2]=(y[0]-y[1])/(x[0][1]-x[2][1])*(zz-x[0][1])+y[0];   // interpolation in z

		      delta = dd*kpc2cm *y[2];
		    }
		}
	      if(gcr[0].n_spatial_dimensions==3) 
		{ 	      
		  if (ix==galaxy.n_xgrid-1 || iy==galaxy.n_ygrid-1 || iz==galaxy.n_zgrid-1)
		    delta = dd*kpc2cm *galaxy.DM_emiss.d3[ix][iy][iz].s[iEgamma];
		  
		  else  // linear interpolation
		    {
		      //  x[0]=(x0,z1,y0), x[1]=(x1,z1,y0);   x[4]=(x0,z1,y1), x[5]=(x1,z1,y1);
		      //  x[2]=(x0,z0,y0), x[3]=(x1,z0,y0);   x[6]=(x0,z0,y1), x[7]=(x1,z0,y1);
		      x[0][0]=galaxy.x[ix  ]; x[0][1]=galaxy.z[iz+1]; x[0][2]=galaxy.y[iy  ];  f[0]=galaxy.DM_emiss.d3[ix  ][iy  ][iz+1].s[iEgamma];
		      x[1][0]=galaxy.x[ix+1]; x[1][1]=galaxy.z[iz+1]; x[1][2]=galaxy.y[iy  ];  f[1]=galaxy.DM_emiss.d3[ix+1][iy  ][iz+1].s[iEgamma];
		      x[2][0]=galaxy.x[ix  ]; x[2][1]=galaxy.z[iz  ]; x[2][2]=galaxy.y[iy  ];  f[2]=galaxy.DM_emiss.d3[ix  ][iy  ][iz  ].s[iEgamma];
		      x[3][0]=galaxy.x[ix+1]; x[3][1]=galaxy.z[iz  ]; x[3][2]=galaxy.y[iy  ];  f[3]=galaxy.DM_emiss.d3[ix+1][iy  ][iz  ].s[iEgamma];
		      x[4][0]=galaxy.x[ix  ]; x[4][1]=galaxy.z[iz+1]; x[4][2]=galaxy.y[iy+1];  f[4]=galaxy.DM_emiss.d3[ix  ][iy+1][iz+1].s[iEgamma];
		      x[5][0]=galaxy.x[ix+1]; x[5][1]=galaxy.z[iz+1]; x[5][2]=galaxy.y[iy+1];  f[5]=galaxy.DM_emiss.d3[ix+1][iy+1][iz+1].s[iEgamma];
		      x[6][0]=galaxy.x[ix  ]; x[6][1]=galaxy.z[iz  ]; x[6][2]=galaxy.y[iy+1];  f[6]=galaxy.DM_emiss.d3[ix  ][iy+1][iz  ].s[iEgamma];
		      x[7][0]=galaxy.x[ix+1]; x[7][1]=galaxy.z[iz  ]; x[7][2]=galaxy.y[iy+1];  f[7]=galaxy.DM_emiss.d3[ix+1][iy+1][iz  ].s[iEgamma];
		      
		      y[0]=  (f[0]-f[1])/(x[0][0]-x[1][0])*(xx-x[0][0])+f[0]; // interpolation in x
		      y[1]=  (f[2]-f[3])/(x[2][0]-x[3][0])*(xx-x[2][0])+f[2];
		      y[2]=  (f[4]-f[5])/(x[4][0]-x[5][0])*(xx-x[4][0])+f[4];
		      y[3]=  (f[6]-f[7])/(x[6][0]-x[7][0])*(xx-x[6][0])+f[6];
		      
		      y[4]=(y[0]-y[1])/(x[0][1]-x[2][1])*(zz-x[0][1])+y[0];   // interpolation in z
		      y[5]=(y[2]-y[3])/(x[4][1]-x[6][1])*(zz-x[4][1])+y[2];
		      
		      y[6]=(y[4]-y[5])/(x[0][2]-x[4][2])*(yy-x[0][2])+y[4];   // interpolation in y
		      
//		  cout<<" y= "; for(int j=0;j<7;j++) cout<<y[j]<<" "; cout<<" >> "<<xx+zz+yy<<endl; 
//		  for(int j=0;j<3;j++) {cout<<" x,f= "; for(int i=0;i<8;i++) cout<<x[i][j]<<" "<<f[i]<<"   "; cout<<endl;} exit(0);
		      
		      delta = dd*kpc2cm *y[6];
		    }
		}

	  DM[iEgamma] +=delta;
	}            
      
    }//complete==0
  return 0;
}
