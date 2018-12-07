
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * isrf_energy_density.cc *                      galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include<iostream>
#include"galprop_classes.h"

float Galprop::isrf_energy_density(float rr, float zz)
{
//  cout<<"isrf_energy_density  "<<rr<<" "<<zz<<endl;

  int i_comp=isrf_energy_density_i_comp; // in global
  float energy_density=0.0;
  
  int iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.5); //IMOS20060420 corrected "0.5"

  if(iz<0 || iz>galaxy.n_zgrid-1) 
    { 
      cout<<"isrf_energy_density iz out of range"<<endl;
      exit(1);
    }
  
  if(galaxy.n_spatial_dimensions==2)
    {
      int ir=(int)((rr-galaxy.r_min)/galaxy.dr + 0.5); //IMOS20060420 corrected "0.5"
//cout<<"isrf_energy_density i_comp rr zz ir iz energy_density "<<i_comp<<" "<<rr<<" "<<zz<<" "<<ir<<" "<<iz<<" "<<energy_density<<endl;
      if(ir<0 || ir>galaxy.n_rgrid-1)
	{
	  cout<<"isrf_energy_density rr ir = "<<rr<<" "<<ir<<" out of range"<<endl; 
	  exit(1);
	}
      energy_density = galaxy.ISRF_energy_density[i_comp].d2[ir][iz].s[0];
//cout<<"isrf_energy_density i_comp rr zz ir iz energy_density "<<i_comp<<" "<<rr<<" "<<zz<<" "<<ir<<" "<<iz<<" "<<energy_density<<endl;
    }
  
  if(galaxy.n_spatial_dimensions==3)
    {
      int ix=(int)((rr-galaxy.x_min)/galaxy.dx + 0.5); //IMOS20060412 corrected "0.5"
      if(ix<0 || ix>galaxy.n_xgrid-1)
	{
	  cout<<"isrf_energy_density ix out of range"<<endl;
	  exit(1);
	}
      int iy=0;
      energy_density = galaxy.ISRF_energy_density[i_comp].d3[ix][iy][iz].s[0];
      
cout<<"isrf_energy_density rr zz ix iz energy_density "<<rr<<" "<<zz<<" "<<ix<<" "<<iz<<" "<<energy_density<<endl;
    }
  return energy_density;
}
