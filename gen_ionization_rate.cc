
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_ionization_rate.cc *                      galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624


#include"galprop_classes.h"
#include"galprop_internal.h"

// generate ionization rate
/*
CR density  gcr.cr_density is in c/4pi * n(E) [ cm s^-1 sr^-1 cm^-3 MeV^-1]

emissivity (cm^-3 s^-1 sr^-1 MeV^-1)=
(c/4pi)*integral[sigma{Egamma,p_beam }  ) n(E)E dlog(E)]

pp_meson has Egamma in GeV, beam momentum in GeV
The particle spectra are assumed to be on equal kinetic energy per nucleon grids
which is the standard for galprop.
BUT UNITS OF density/momentum = flux/(KE/nucleon)..... CHECK CAREFULLY, ALSO factor A
cross section from pp_meson in barns GeV^-1
factor= 1.0e-24* *1.0e-3 log(Ekin_factor)
*/

int Galprop::gen_ionization_rate()
{
  cout<<" >>>>gen_ionization_rate"<<endl;
  cout<<"generating ionization rate for n_spatial_dimensions="<<gcr[0].n_spatial_dimensions<<endl;

  int stat=0;
  double ion_cs; 
  galaxy.ionization_rate = 0.0;
  
  for(int i=0;i<n_species;i++)
  {
    for(int ip=0;ip<gcr[i].n_pgrid;ip++)
      {
        if(gcr[i].A!=0) ion_cs=ionization_bethe(gcr[i].Z,gcr[i].beta[ip]);
        if(gcr[i].A==0) ion_cs=ionization_bethe(gcr[i].Z,gcr[i].beta[ip]);//AWS20000628

	//ion_cs=0.0; // for future inclusion of electrons//AWS20000628

        cout<<gcr[i].name<<" "<<"  Z Ekin beta ionization cross section: "
	    <<gcr[i].Z<<" "<<gcr[i].Ekin[ip]<<" "<<gcr[i].beta[ip]<<" "<<ion_cs<<endl;	
 
	if(galaxy.n_spatial_dimensions==2)
	  for(int ir=0;ir<gcr[i].n_rgrid;ir++)
	    for(int iz=0;iz<gcr[i].n_zgrid;iz++)
	      galaxy.ionization_rate.d2[ir][iz].s[0]+= ion_cs *gcr[i].cr_density.d2[ir][iz].s[ip] *gcr[i].Ekin[ip]; 
	
	if(galaxy.n_spatial_dimensions==3)
	  for(int ix=0;ix<gcr[i].n_xgrid;ix++)
	    for(int iy=0;iy<gcr[i].n_ygrid;iy++)
	      for(int iz=0;iz<gcr[i].n_zgrid;iz++)
		galaxy.ionization_rate.d3[ix][iy][iz].s[0]+= ion_cs *gcr[i].cr_density.d3[ix][iy][iz].s[ip] *gcr[i].Ekin[ip];  
      }//ip
  }// species
  
  double factor=4.0*Pi* log(galdef.Ekin_factor);
  galaxy.ionization_rate *= factor;
 
  if(galdef.verbose>=1)
  {
    cout<<"   ionization_rate  "<<endl;
    galaxy.ionization_rate.print();
  }//galdef.verbose

  cout<<" <<<< gen_ionization_rate"<<endl;
  return stat;
}
