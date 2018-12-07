
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_isrf_energy_density.cc *                  galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std; //AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"

// generate ISRF energy density for all components
/* 
ISRF is in Hz eV cm-3 Hz-1
integral energy density Hz-1  d(nu) = integral (nu* energy density Hz-1) d(log nu)
d(log nu) is constant in this ISRF
factor= LOG(nu(2)/nu(1)) 
*/

#include <ErrorLogger.h>
#include <sstream>
#include <fstream>

int Galprop::gen_isrf_energy_density() {

  INFO("Entry");

  //cout<<" >>>>gen_isrf_energy_density"<<endl;
  
  int ir,ix,iy,iz,i_comp,stat=0;
  float sum;
 
  delete[] galaxy.ISRF_energy_density;
  galaxy.ISRF_energy_density=new Distribution[galaxy.n_ISRF_components]; 
  
  double factor= log(galaxy.nu_ISRF[1]/ galaxy.nu_ISRF[0] );
  
  if(galaxy.n_spatial_dimensions==2)
    {
       ostringstream buf;
      buf<<"generating ISRF energy density for n_spatial_dimensions="<<galaxy.n_spatial_dimensions;
      INFO(buf.str());
      
      for(i_comp=0;i_comp<galaxy.n_ISRF_components;i_comp++)
	galaxy.ISRF_energy_density[i_comp].init(galaxy.n_rgrid,galaxy.n_zgrid,1);
    
      for(ir=0;ir<galaxy.n_rgrid;ir++)
	{
	  for(iz=0;iz<galaxy.n_zgrid;iz++)
	    {
	      for(i_comp=0;i_comp<galaxy.n_ISRF_components;i_comp++)
		{
//		  double sum1=0., sum2=0.;
		  sum=0.0;
		  for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++) sum+=galaxy.ISRF[i_comp].d2[ir][iz].s[inu];
		  galaxy.ISRF_energy_density[i_comp].d2[ir][iz].s[0]=sum*factor;
		  /*   //calculation of the energy density and the number density of photons for each component
		  for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++) sum1+=galaxy.ISRF[i_comp].d2[ir][iz].s[inu]*factor;
		  for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++) sum2+=galaxy.ISRF[i_comp].d2[ir][iz].s[inu]/(h_planck/eV_to_erg)/galaxy.nu_ISRF[inu]*factor;
		  if(iz==41) cout<<ir<<" isrf average energy = "<<sum1<<" "<<sum2<<" "<<sum1/sum2<<" "<<i_comp<<endl;
		  */
		}//ISRF_components
	    }//iz
	}//ir
    }//particle.n_spatial_dimensions==2

  if(galaxy.n_spatial_dimensions==3)
    { 
       ostringstream buf;
      buf<<"generating ISRF energy density for n_spatial_dimensions="<<galaxy.n_spatial_dimensions;
      INFO(buf.str());
      for(i_comp=0;i_comp<galaxy.n_ISRF_components;i_comp++)
        galaxy.ISRF_energy_density[i_comp].init(galaxy.n_xgrid,galaxy.n_ygrid,galaxy.n_zgrid,1);

      for(ix=0;ix<galaxy.n_xgrid;ix++)
	{
	  for(iy=0;iy<galaxy.n_ygrid;iy++)
	    {
	      for(iz=0;iz<galaxy.n_zgrid;iz++)
		{
 // cout<<" "<<ix<<" "<<iy<<" "<<iz<<endl;
		  for(i_comp=0;i_comp<galaxy.n_ISRF_components;i_comp++)
		    {  
		      sum=0.0;
		      for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++) sum+=galaxy.ISRF[i_comp].d3[ix][iy][iz].s[inu];
		      galaxy.ISRF_energy_density[i_comp].d3[ix][iy][iz].s[0]=sum*factor;
		    }//ISRF_components
		}//iz
	    }//iy
	}//ix
   }//particle.n_spatial_dimensions==3

  /*ofstream fout("/home/yuko/GALPROP/check.tab"); //check
  for(ix=0;ix<galaxy.n_xgrid;ix++)
    for(iy=0;iy<galaxy.n_ygrid;iy++)
      for(iz=0;iz<galaxy.n_zgrid;iz++)
      fout << galaxy.ISRF_energy_density[2].d3[ix][iy][iz].s[0] << endl;*/
  
  if(galdef.verbose>=2)  
    for(i_comp=0;i_comp<galaxy.n_ISRF_components;i_comp++)
      {
        ostringstream buf;
	buf<<"ISRF energy density for component #"<<i_comp;
	INFO(buf.str());
	galaxy.ISRF_energy_density[i_comp].print();
      }

  //check 
  /*  for(i_comp=0;i_comp<galaxy.n_ISRF_components;i_comp++)
    for(ix=0;ix<galaxy.n_xgrid;ix++)
      for(iy=0;iy<galaxy.n_ygrid;iy++)
	for(iz=0;iz<galaxy.n_zgrid;iz++)
	  cerr << "isrf_energy_density= " << galaxy.ISRF_energy_density[i_comp].d3[ix][iy][iz].s[0] << endl;
  */

  //cout<<" <<<< gen_isrf_energy_density"<<endl;
  INFO("Exit");
  return stat;
}
