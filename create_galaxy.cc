
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * create_galaxy.cc *                            galprop package * 5/22/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galprop_internal.h"

#include <ErrorLogger.h>
#include <sstream>

#include"B_field_3D_model.h" //AWS20080326

int Galprop::create_galaxy() {//AWS20050816

  INFO("Entry");

  //cout<<" >>>> create_galaxy"<<endl;
  
  int stat=0;
  double dzz=0.01; // average with this resolution in kpc
  
  if(galdef.n_spatial_dimensions==2)
    galaxy.init(galdef.r_min,  galdef.r_max, galdef.dr,  
		galdef.z_min,  galdef.z_max, galdef.dz);
  
  if(galdef.n_spatial_dimensions==3) {

    if(galdef.use_symmetry==1) 
      galdef.x_min=galdef.y_min=galdef.z_min=0.; // IMOS20020419
    
    galaxy.init(galdef.x_min,  galdef.x_max, galdef.dx,
		galdef.y_min,  galdef.y_max, galdef.dy,
		galdef.z_min,  galdef.z_max, galdef.dz);
  }                  
  // GAS and B-FIELD DISTRIBUTION
  
  nH_set_model(galdef.nHI_model, galdef.nH2_model, galdef.nHII_model, 0); //AWS20090814


  galaxy.n_HI =1.e-6; //  non-zero values to avoid problem in energy loss logarithm
  galaxy.n_H2 =1.e-6;
  galaxy.n_HII=1.e-6;
  
  if (2 == galdef.n_spatial_dimensions) { // use average at y=0
      
    for (int ir = 0; ir < galaxy.n_rgrid; ++ir) {
      
      for (int iz = 0; iz < galaxy.n_zgrid; ++iz) {

	galaxy.n_HI.d2[ir][iz].s[0]=
	  nH_av(galaxy.r[ir], 0, galaxy.z[iz], galaxy.dz, dzz, nHI); //IMOS20080114
	//		nHI_av (galaxy.r[ir],0.0,galaxy.z[iz],galaxy.dz,dzz);
	      
	galaxy.n_H2.d2[ir][iz].s[0]=
	  nH_av(galaxy.r[ir], 0, galaxy.z[iz], galaxy.dz, dzz, fX_CO(galaxy.r[ir]),nH2); //IMOS20080114
	//		nH2_av (galaxy.r[ir],0.0,galaxy.z[iz],galaxy.dz,dzz);
	      
	galaxy.n_HII.d2[ir][iz].s[0]=
	  nH_av(galaxy.r[ir], 0, galaxy.z[iz], galaxy.dz, dzz, nHII); //IMOS20080114
	//		nHII_av(galaxy.r[ir],0.0,galaxy.z[iz],galaxy.dz,dzz); 
	      
	galaxy.B_field.d2[ir][iz].s[0]=
	  B_field_model(galaxy.r[ir], 0, galaxy.z[iz], galdef.B_field_model);

	// 3D model total field  
	if (2 == galdef.synchrotron) { //AWS20080326
	  
	  galaxy.B_field.d2[ir][iz].s[0] = B_field_3D_model_tot(galdef.B_field_name, galdef.B_field_parameters, galaxy.r[ir], galaxy.z[iz])/1.0e4;     // Gauss->Tesla  

	}

      }

    }
    
  }
  
  if (3 == galdef.n_spatial_dimensions) {

    for (int ix = 0; ix < galaxy.n_xgrid; ++ix) {
      
      for (int iy = 0; iy < galaxy.n_ygrid; ++iy) {

	const double r = sqrt(galaxy.x[ix]*galaxy.x[ix] + galaxy.y[iy]*galaxy.y[iy]);

	for (int iz = 0; iz < galaxy.n_zgrid; ++iz) {

	  
	  galaxy.n_HI.d3[ix][iy][iz].s[0] =
	    nH_av(galaxy.x[ix], galaxy.y[iy], galaxy.z[iz], galaxy.dz, dzz, nHI); //IMOS20080114
	  //		    nHI_av (galaxy.x[ix],galaxy.y[iy],galaxy.z[iz],galaxy.dz,dzz);
		  
	  galaxy.n_H2.d3[ix][iy][iz].s[0] =
	    nH_av(galaxy.x[ix], galaxy.y[iy], galaxy.z[iz], galaxy.dz, dzz, fX_CO(r), nH2); //IMOS20080114
	  //		    nH2_av (galaxy.x[ix],galaxy.y[iy],galaxy.z[iz],galaxy.dz,dzz);
		  
	  galaxy.n_HII.d3[ix][iy][iz].s[0] =
	    nH_av(galaxy.x[ix], galaxy.y[iy], galaxy.z[iz], galaxy.dz, dzz, nHII);  //IMOS20080114
	  //		    nHII_av(galaxy.x[ix],galaxy.y[iy],galaxy.z[iz],galaxy.dz,dzz);  
		  
	  galaxy.B_field.d3[ix][iy][iz].s[0] =
	    B_field_model(galaxy.x[ix], galaxy.y[iy], galaxy.z[iz], galdef.B_field_model);

	  // 3D model total field       
	  if (2 == galdef.synchrotron) { //AWS20080326
	    
	    galaxy.B_field.d3[ix][iy][iz].s[0] = B_field_3D_model_tot(galdef.B_field_name, galdef.B_field_parameters, galaxy.x[ix], galaxy.y[iy], galaxy.z[iz])/1.0e4;   

	  }
	   
	}
	
      }
    
    }
  
  }

  if (galdef.verbose>=1) {

      INFO("galaxy.n_HI:   ");galaxy.n_HI.print();
      INFO("galaxy.n_H2:   ");galaxy.n_H2.print();
      INFO("galaxy.n_HII:  ");galaxy.n_HII.print();
      INFO("galaxy.B_field:");galaxy.B_field.print();

  }
  
  // SKYMAP PARAMETERS
  // moved to before gas surveys since required there                                    AWS20050913 
  
  galaxy.d_long = galdef.d_long;
  galaxy.long_min = galdef.long_min;
  galaxy.long_max = galdef.long_max;
  
  galaxy.d_lat = galdef.d_lat;
  galaxy.lat_min = galdef.lat_min;
  galaxy.lat_max = galdef.lat_max;
  
  galaxy.n_long = int((galaxy.long_max-galaxy.long_min)/galaxy.d_long + 1.001);
  galaxy.n_lat = int((galaxy. lat_max-galaxy. lat_min)/galaxy.d_lat  + 1.001);
  
  ostringstream buf;
  buf<<"    gamma-ray, synchrotron skymap arrays: n_long,nlat="<<galaxy.n_long<<" "<<galaxy.n_lat;
  DEBUGLOG(buf.str());
  
  /*Moved to Galprop before generating skymaps  //Gulli20070810
  // READING GAS SURVEYS
  read_HIR();            //AWS20041213
  read_COR();            //AWS20041213
  */
  
  
  // READING THE INTERSTELLAR RADIATION FIELD
  
  if (galdef.ISRF_factors) {

    galaxy.fISRFFactors.resize(3);
   
    for (unsigned long i = 0; i < 3; ++i)
      galaxy.fISRFFactors[i] = galdef.ISRF_factors[i];

  }

  stat = read_isrf(galdef.ISRF_filetype);
  
  if (gen_isrf_energy_density() != 0) 
    return 1;
  
  // STOCHASTIC SNR
  
  if (3 == galdef.n_spatial_dimensions && 1 == galdef.SNR_events) {

    //Moved from Galaxy::init    Gulli20070810
    galaxy.SNR_cell_time.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, 1);
    galaxy.SNR_cell_phase.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, 1);
    
    galaxy.SNR_electron_dg.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, 1); //AWS20010410
    galaxy.SNR_nuc_dg.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, 1); //AWS20010410
    create_SNR();
  
  }

  // GAMMA RAYS
  
	//Moved to correct places in Galprop   Gulli20070810
  /*if (galdef.gamma_rays >= 1) {                    //AWS20050302 IMOS20060420
    
    galaxy.E_gamma_min = galdef.E_gamma_min;
    galaxy.E_gamma_max = galdef.E_gamma_max;
    galaxy.E_gamma_factor = galdef.E_gamma_factor;
    
    galaxy.n_E_gammagrid = int(log(galaxy.E_gamma_max/galaxy.E_gamma_min)/
			       log(galaxy.E_gamma_factor) + 1.001);
    
    galaxy.E_gamma = new double[galaxy.n_E_gammagrid];

    for (int iEgamma = 0; iEgamma < galaxy.n_E_gammagrid; ++iEgamma)
      galaxy.E_gamma[iEgamma] = exp(log(galaxy.E_gamma_min) + 
				    iEgamma*log(galaxy.E_gamma_factor));

    cout<<"   creating gamma-ray emissivity arrays"<<endl;
    galaxy.IC_iso_emiss = new Distribution[galaxy.n_ISRF_components];
    galaxy.IC_aniso_emiss = new Distribution;

    if (galdef.n_spatial_dimensions == 2) {

      if (galdef.bremss) {
	
	galaxy.bremss_emiss.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid); // IMOS20060420
	galaxy.bremss_ionized_emiss.init(galaxy.n_rgrid,galaxy.n_zgrid,galaxy.n_E_gammagrid); // IMOS20060420

      }

      if (galdef.IC_isotropic + galdef.IC_anisotropic) { // IMOS20060420
	
	for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {
	  
	  galaxy.IC_iso_emiss[icomp].init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);

	}

      }

      if (galdef.pi0_decay) { 
	
	galaxy.pi0_decay_emiss.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);   // IMOS20060420

      }
	 
      if (galdef.DM_gammas) {
	
	galaxy.DM_emiss.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid); // DM: IMOS20050912 IMOS20060420
       
      }
     
    }

    if (galdef.n_spatial_dimensions == 3) {

      if (galdef.bremss) {

	galaxy.bremss_emiss.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid); // IMOS20060420
        
	galaxy.bremss_ionized_emiss.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid); // IMOS20060420

      }

      if (galdef.IC_isotropic + galdef.IC_anisotropic) {                                                                       // IMOS20060420
	for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {

	  galaxy.IC_iso_emiss[icomp].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);

	}

      }

      if (galdef.pi0_decay) { 

	galaxy.pi0_decay_emiss.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);   // IMOS20060420	 
      }

      if (galdef.DM_gammas) {

	galaxy.DM_emiss.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid); // DM: IMOS20050912 IMOS20060420
       
      }
     
    }

    // GAMMA-RAY SKYMAPS
    cout<<"   creating gamma-ray skymap arrays"<<endl;

    if (galdef.bremss) {  // IMOS20060420
	
      galaxy.bremss_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid);
      galaxy.bremss_HIR_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.HIR.n_zgrid, galaxy.n_E_gammagrid);   //AWS20041214

      galaxy.bremss_H2R_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.COR.n_zgrid, galaxy.n_E_gammagrid);   //AWS20041214
	  
      galaxy.bremss_ionized_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid);
	
    }

    if (galdef.IC_isotropic + galdef.IC_anisotropic) { // IMOS20060420
	
      galaxy.IC_iso_skymap = new Distribution[galaxy.n_ISRF_components];

      for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {

	galaxy.IC_iso_skymap[icomp].init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid);
	  
      }

      if (galdef.IC_anisotropic) { // IMOS20060420
	    
	galaxy.IC_aniso_skymap = new Distribution[galaxy.n_ISRF_components];
	
	for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {
	  
	  galaxy.IC_aniso_skymap[icomp].init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid);
	
	}
	
      }

    }

    if (galdef.pi0_decay) {// IMOS20060420
 
      galaxy.pi0_decay_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid);
      galaxy.pi0_decay_HIR_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.HIR.n_zgrid, galaxy.n_E_gammagrid);   //AWS20041214

      galaxy.pi0_decay_H2R_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.COR.n_zgrid, galaxy.n_E_gammagrid);   //AWS20041214
      
    }
    
    if(galdef.DM_gammas) {

      galaxy.DM_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid); // DM: IMOS20050912 IMOS20060420
   
    }

  }*/

// SYNCHROTRON EMISSION

  if (galdef.synchrotron >= 1) {                          //AWS20080229
     
    galaxy.nu_synch_min =   galdef.nu_synch_min;
    galaxy.nu_synch_max =   galdef.nu_synch_max;
    galaxy.nu_synch_factor = galdef.nu_synch_factor;
    galaxy.n_nu_synchgrid = int(log(galaxy.nu_synch_max/galaxy.nu_synch_min)/
				log(galaxy.nu_synch_factor) + 1.001);
    
    galaxy.nu_synch = new double[galaxy.n_nu_synchgrid];
    //    if(galdef.verbose >= 1){ bug introduced in r330, fixed by AWS20091015

    //cout<<"synchrotron frequency grid: ";

    ostringstream sBuf;

    sBuf << "Synchrotron frequency grid: ";

    for (int inusynch = 0; inusynch < galaxy.n_nu_synchgrid; ++inusynch) {
      
      galaxy.nu_synch[inusynch] = exp(log(galaxy.nu_synch_min) + 
				      inusynch*log(galaxy.nu_synch_factor));
      //cout<<galaxy.nu_synch[inusynch]<<" ";
      
      //    }

      sBuf << galaxy.nu_synch[inusynch] << " ";
      
    
      //cout<<endl;
      //cout<<"   creating synchrotron emissivity array"<<endl;
    }

    INFO(sBuf.str());

		/* Moved to Galprop.cc  Gulli20071003
    if (galdef.n_spatial_dimensions == 2) {

      galaxy.synchrotron_emiss.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_nu_synchgrid);

    }

    if (galdef.n_spatial_dimensions == 3) {
       
      galaxy.synchrotron_emiss.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_nu_synchgrid);

    }

    cout<<"   creating synchrotron skymap     array"<<endl;

    galaxy.synchrotron_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_nu_synchgrid);
		*/

  }

	/* Moved to Galprop.cc Gulli20071003
  cout<<"   creating ionization rate array"<<endl;

  if (galdef.n_spatial_dimensions == 2) {
    
    galaxy.ionization_rate.init(galaxy.n_rgrid, galaxy.n_zgrid,1);

  }

  if (galdef.n_spatial_dimensions == 3) {

    galaxy.ionization_rate.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, 1);

  }
	*/

  if(galdef.verbose >= 1){
     INFO("      galaxy.print");
     galaxy.print();
  }

  INFO("Exit");

  //cout<<" <<<< create_galaxy"<<endl;
  return stat;

}
