#include"galprop_classes.h"
#include"galprop_internal.h"

#include <fort_interface.h>

#include <ErrorLogger.h>
#include <Timer.h>

#include <config.h>

#include <iostream>
#include <sstream>

using namespace std;

Galprop* gGalprop = 0;

Galprop::Galprop() {

  gGalprop = this;
  gcr = 0;

}

Galprop::~Galprop() {

}

int Galprop::Run(const string& galdefPath,
		 const string& fitsPath,
		 const string& outputPath,
		 const string& outputPrefix,
		 const string& runNumber) {

  /*#ifdef VERSION
  const string fullVersion = VERSION;
#else
  const string fullVersion = "54.0.0";
#endif

  ostringstream bufV;
  bufV << "Entry: this is Galprop version " << fullVersion;
  INFO(bufV.str());

  //  pre-r151
  //  char version[]="53"; //AWS 20070820
  //  for changes from r151, new version defined at r158
  
  //char version[]="54"; //AWS 20080618

  const string majorVersion = fullVersion.substr(0, fullVersion.find_first_of("."));

  const string tmpVer = fullVersion.substr(fullVersion.find_first_of(".") + 1, fullVersion.size() - 1);

  const string minorVersion = tmpVer.substr(0, tmpVer.find_first_of("."));

  const string revision = tmpVer.substr(tmpVer.find_first_of(".") + 1, tmpVer.size() - 1);

  //cout << majorVersion << " " << minorVersion << " " << revision << endl;

  // Check for overrides. Otherwise grab whatever comes from the 
  // autoconfiguration, and if that fails use the stupid defaults we've
  // always had

  string galdefDir = galdefPath, fitsDir = fitsPath, outputDir = outputPath;

  if (galdefDir.empty()) {

#ifdef HAVE_GALDEF
    galdefDir = GALDEF_PATH;
#else
    galdefDir = "../GALDEF/";
#endif

  } 

  if (fitsDir.empty()) {

#ifdef HAVE_FITSDATA
    fitsDir = FITSDATA_PATH;
#else
    fitsDir = "../FITS/";
#endif

  }

  if (outputDir.empty()) {

    outputDir = "../FITS"; // *sigh*

  }
  */

  if (configure.init(galdefPath, fitsPath, outputPath, outputPrefix)) {

    FATAL("Internal error. Fix data paths!");
    return 1;

  }

  if (galdef.read(configure.fVersion, runNumber, configure.fGaldefDirectory)) {

    FATAL("Internal error. Problem reading from galdef file!");
    return 1;

  } 

  //cout << "Parameters to Run: " << galdefDir << " " << dataDir << " " << outputDir << " " << outputPrefix << " " << runNumber << endl;

  //exit(0);

  //if(galdef.read(version,argv[1],configure.galdef_directory) !=0) return 1;
  
  // electron propagation test vs analytical formula; assign galdef parameters IMOS20061030

  if (99 == abs(galdef.DM_int0)) {

    cout<<endl<<" >>>> Running electron propagation test"<<endl<<endl;
    galdef.n_spatial_dimensions =2;
    galdef.D_rigid_br           =1.e3;
    galdef.diff_reacc           =0;
    galdef.convection           =0;
    galdef.momentum_losses      =1;
    galdef.max_Z                =1;
    galdef.primary_electrons    =0;
    galdef.secondary_positrons  =0;
    galdef.secondary_electrons  =0;
    galdef.knock_on_electrons   =0;
    galdef.secondary_antiprotons=0;
    galdef.tertiary_antiprotons =0;
    galdef.gamma_rays           =0;
    galdef.DM_positrons         =1; //will contain analytical solution
    galdef.DM_electrons         =1; //will contain numerical  solution
    galdef.DM_antiprotons       =0;
    galdef.DM_gammas            =0;
    cout<<"** Test values are re-assigned in Galprop.cc: **"      <<endl;
    cout<<"  DM_int0               "<<galdef.DM_int0              <<endl;
    cout<<"  DM_double6            "<<galdef.DM_double6           <<endl;
    cout<<"  DM_double7            "<<galdef.DM_double7           <<endl;
    cout<<"  DM_double8            "<<galdef.DM_double8           <<endl;
    cout<<"  DM_double9            "<<galdef.DM_double9           <<endl;
    cout<<"  DM_positrons          "<<galdef.DM_positrons         <<endl;
    cout<<"  DM_electrons          "<<galdef.DM_electrons         <<endl;
    cout<<"  DM_antiprotons        "<<galdef.DM_antiprotons       <<endl;
    cout<<"  DM_gammas             "<<galdef.DM_gammas            <<endl;
    cout<<"  D_rigid_br            "<<galdef.D_rigid_br           <<endl;
    cout<<"  diff_reacc            "<<galdef.diff_reacc           <<endl;
    cout<<"  convection            "<<galdef.convection           <<endl;
    cout<<"  momentum_losses       "<<galdef.momentum_losses      <<endl;
    cout<<"  max_Z                 "<<galdef.max_Z                <<endl;
    cout<<"  primary_electrons     "<<galdef.primary_electrons    <<endl;
    cout<<"  secondary_positrons   "<<galdef.secondary_positrons  <<endl;
    cout<<"  secondary_electrons   "<<galdef.secondary_electrons  <<endl;
    cout<<"  knock_on_electrons    "<<galdef.knock_on_electrons   <<endl;
    cout<<"  secondary_antiproton  "<<galdef.secondary_antiprotons<<endl;
    cout<<"  tertiary_antiproton   "<<galdef.tertiary_antiprotons <<endl;
    cout<<"  gamma_rays            "<<galdef.gamma_rays           <<endl<<endl;
    cout<<"  n_spatial_dimensions  "<<galdef.n_spatial_dimensions <<endl;
  
  }
    
  // global test
  if (0 < galdef.test_suite) { 

    test_suite(); 
    return 0; 

  }
  
  //reading all isotopic cross sections, nuclear reaction network, nuclear fits
  read_nucdata(configure.fGlobalDataPath);
  
  //initialization of the Barashenkov & Polanski cross section code 
  sigtap_cc(-1, configure.fGlobalDataPath);             // IMOS20010511 AWS20010620
  
  //initialization of the Webber's routine
  set_sigma_cc(configure.fGlobalDataPath);

  Timer::initialize();

  int stat;
  TIME_FUNCTION(stat,create_galaxy);
  if (0 != stat) 
    return 1;

  TIME_FUNCTION(stat,create_gcr);
  if (0 != stat) 
    return 1;

  //major routine
  TIME_FUNCTION(stat,propagate_particles);
  if (0 != stat)
    return 1;

  cr_luminosity();
  
  TIME_FUNCTION(stat,store_gcr);
  if (0 != stat) 
     return 1;                                      //GJ20100104
  //deleting all cross section arrays etc.
  cleanup_nucdata();
  
  //Cleaning up SNR      Gulli20070810
  if (galdef.SNR_events) {

    galaxy.SNR_cell_time.delete_array();
    galaxy.SNR_cell_phase.delete_array();
    galaxy.SNR_electron_dg.delete_array();
    galaxy.SNR_nuc_dg.delete_array();
  
  }
  
  //   if(print_BC() !=0) return 1;
  //   if(store_gcr() !=0) return 1; //IMOS20030129
  
  if (galdef.output_gcr_full) {
    TIME_FUNCTION(stat,store_gcr_full);
    if (0!= stat) 
      return 1;
  }

  //exit(0);

  //Test chisq function  //Gulli20090421
  if ( galdef.verbose == -701 ) {
     bool sets[9] = {0,0,1,1,0,0,0,0,0};
     bool expmnts[16] = {0,0,0,0,1,0,0,1,0,0,0,0,0,1,1,1};
     double modulation[16] = {500.,500.,500.,500.,500.,500.,500.,500.,500.,500.,500.,500.,500.,500.,500.,500.};
     double tau[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
     double energy[16] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
     double chisq[16];
     cout<<"CHISQ: "<<calculate_GCR_chisq(sets,expmnts,modulation,tau,energy,chisq)<<endl;
     return 0;
  }
   
  
  if (galdef.synchrotron) { //AWS20080229
    
    if (galdef.primary_electrons || 
	galdef.secondary_electrons || 
	galdef.knock_on_electrons || // IMOS20060504
	galdef.secondary_positrons || 
	galdef.DM_positrons || 
	galdef.DM_electrons) {        // IMOS20050912
	
      //Moved from create_galaxy  Gulli20071003
      
      if (2 == galdef.n_spatial_dimensions) {

	galaxy.synchrotron_emiss  .init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_nu_synchgrid);
	galaxy.synchrotron_Q_emiss.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_nu_synchgrid);                  //AWS20100708
	galaxy.synchrotron_U_emiss.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_nu_synchgrid);                  //AWS20100708
      }

      if (3 == galdef.n_spatial_dimensions) {

	galaxy.synchrotron_emiss  .init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_nu_synchgrid);
	galaxy.synchrotron_Q_emiss.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_nu_synchgrid);  //AWS20100708
	galaxy.synchrotron_U_emiss.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_nu_synchgrid);  //AWS20100708
	
      }

      INFO("Creating synchrotron skymap array");
	  
      if (3 == galdef.skymap_format) {

	galaxy.synchrotron_hp_skymap  .Resize(galdef.healpix_order, galaxy.n_nu_synchgrid);
	galaxy.synchrotron_Q_hp_skymap.Resize(galdef.healpix_order, galaxy.n_nu_synchgrid);                      //AWS20100708
	galaxy.synchrotron_U_hp_skymap.Resize(galdef.healpix_order, galaxy.n_nu_synchgrid);                      //AWS20100708
	galaxy.synchrotron_P_hp_skymap.Resize(galdef.healpix_order, galaxy.n_nu_synchgrid);                      //AWS20100708
	
      } else {
	    
	galaxy.synchrotron_skymap  .init(galaxy.n_long, galaxy.n_lat, galaxy.n_nu_synchgrid);
	galaxy.synchrotron_Q_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_nu_synchgrid);                    //AWS20100708
	galaxy.synchrotron_U_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_nu_synchgrid);                    //AWS20100708
	galaxy.synchrotron_P_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_nu_synchgrid);                    //AWS20100708
	
      }
	  
      TIME_FUNCTION(stat,gen_synch_emiss);
      if (stat != 0) return 1;
      TIME_FUNCTION(stat,gen_synch_skymap);
      if (stat != 0) return 1;
      TIME_FUNCTION(stat,store_synch_skymap);
      if (stat != 0) return 1;
      TIME_FUNCTION(stat,store_synch_emiss);
      if (stat != 0) return 1;       //AWS20080314
	  
      //Delete Gulli20071003
      galaxy.synchrotron_emiss  .delete_array();
      galaxy.synchrotron_Q_emiss.delete_array();                                                                 //AWS20100708
      galaxy.synchrotron_U_emiss.delete_array();                                                                 //AWS20100708

      if (3 == galdef.skymap_format) 
      {
	galaxy.synchrotron_hp_skymap  .Resize(0,0);
	galaxy.synchrotron_Q_hp_skymap.Resize(0,0);                                                              //AWS20100708
	galaxy.synchrotron_U_hp_skymap.Resize(0,0);                                                              //AWS20100708
	galaxy.synchrotron_P_hp_skymap.Resize(0,0);                                                              //AWS20100708
      }
      else                       
      {  
	galaxy.synchrotron_skymap  .delete_array();
	galaxy.synchrotron_Q_skymap.delete_array();                                                              //AWS20100708
	galaxy.synchrotron_U_skymap.delete_array();                                                              //AWS20100708
	galaxy.synchrotron_P_skymap.delete_array();                                                              //AWS20100708
      }  
    } //electrons
    
  }//galdef.synchrotron
  
  if (galdef.gamma_rays) { //AWS20050302
    
    //Moved from create_galaxy   Gulli20070810
    galaxy.E_gamma_min = galdef.E_gamma_min;
    galaxy.E_gamma_max = galdef.E_gamma_max;
    galaxy.E_gamma_factor = galdef.E_gamma_factor;
      
    galaxy.n_E_gammagrid = int(log(galaxy.E_gamma_max/galaxy.E_gamma_min)/log(galaxy.E_gamma_factor) + 1.001);
      
    galaxy.E_gamma = new double[galaxy.n_E_gammagrid];
      
    for (int iEgamma = 0; iEgamma < galaxy.n_E_gammagrid; ++iEgamma)
      galaxy.E_gamma[iEgamma] = 
	exp(log(galaxy.E_gamma_min) +iEgamma*log(galaxy.E_gamma_factor));
      
    char cor[]={"COR"}, hir[]={"HIR"}; // IMOS20080114

    if (galdef.pi0_decay) { // IMOS20060420
      	  
      TIME_SUBROUTINE(read_gas_maps, hir);                // IMOS20080114
      TIME_SUBROUTINE(read_gas_maps, cor);                // IMOS20080114
      
      //Moved from create_galaxy    Gulli20070810
      if (galdef.n_spatial_dimensions == 2)
	galaxy.pi0_decay_emiss.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);  

      if (galdef.n_spatial_dimensions == 3)
	galaxy.pi0_decay_emiss.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);

      //New healpix format
      if (3 == galdef.skymap_format) {

	galaxy.pi0_decay_hp_skymap.Resize(galdef.healpix_order, galaxy.n_E_gammagrid);

	if (galdef.gamma_rays == 2) {

	  int n_Ring = galaxy.hpHIR.nSpectra();
	  cout<<n_Ring<<endl;
	  galaxy.pi0_decay_H2R_hp_skymap = new Skymap<double>[galaxy.hpHIR.nSpectra()];
	  galaxy.pi0_decay_HIR_hp_skymap = new Skymap<double>[galaxy.hpHIR.nSpectra()];
	  galaxy.pi0_decay_HII_hp_skymap = new Skymap<double>[galaxy.hpHIR.nSpectra()];// IMOS20080114

	  for (int i_Ring=0; i_Ring<galaxy.hpHIR.nSpectra(); ++i_Ring) { 
	    
	    galaxy.pi0_decay_H2R_hp_skymap[i_Ring].Resize(galdef.healpix_order, galaxy.n_E_gammagrid);
	    galaxy.pi0_decay_HIR_hp_skymap[i_Ring].Resize(galdef.healpix_order, galaxy.n_E_gammagrid);
	    galaxy.pi0_decay_HII_hp_skymap[i_Ring].Resize(galdef.healpix_order, galaxy.n_E_gammagrid);// IMOS20080114
	  
	  }
	  
	}
	
      } else {//l and b format
	
	galaxy.pi0_decay_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid);
	
	if (2 == galdef.gamma_rays) {
	
	  galaxy.pi0_decay_H2R_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.COR.n_zgrid, galaxy.n_E_gammagrid);
		
	  galaxy.pi0_decay_HIR_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.HIR.n_zgrid, galaxy.n_E_gammagrid);
		
	  galaxy.pi0_decay_HII_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.HIR.n_zgrid, galaxy.n_E_gammagrid);// IMOS20080114
	  
	}
	
      }
	  
      TIME_FUNCTION(stat,gen_pi0_decay_emiss);
      if (stat != 0) return 1;
      TIME_FUNCTION(stat,store_pi0_decay_emiss);
      if (stat != 0) return 1;
      TIME_FUNCTION(stat,gen_pi0_decay_skymap);
      if (stat != 0) return 1;
      TIME_FUNCTION(stat,store_pi0_decay_skymap);
      if (stat != 0) return 1;
	  
      if (2 == galdef.gamma_rays) {                           //AWS20050302
	    
        TIME_FUNCTION(stat,store_pi0_decay_H2R_skymap);
	if (stat != 0)  return 1;//AWS20041215
        TIME_FUNCTION(stat,store_pi0_decay_HIR_skymap);
	if (stat != 0)  return 1;//AWS20041215
        TIME_FUNCTION(stat,store_pi0_decay_HII_skymap);
	if (stat != 0)  return 1;//IMOS20080114
	
      }

      galaxy.pi0_decay_emiss.delete_array();  //Gulli20070810

      if (3 == galdef.skymap_format) {

	galaxy.pi0_decay_hp_skymap.Resize(0,0);

	if (2 == galdef.gamma_rays) {

	  delete[] galaxy.pi0_decay_H2R_hp_skymap;
	  galaxy.pi0_decay_H2R_hp_skymap = 0;
	  delete[] galaxy.pi0_decay_HIR_hp_skymap;
	  galaxy.pi0_decay_HIR_hp_skymap = 0;
	  delete[] galaxy.pi0_decay_HII_hp_skymap;//IMOS20080114
	  galaxy.pi0_decay_HII_hp_skymap = 0;
	  
	}
	
      } else {
	
	if (2 == galdef.gamma_rays) {                      //AWS20050302
	      
	  galaxy.pi0_decay_H2R_skymap.delete_array();  //Gulli20070810
	  galaxy.pi0_decay_HIR_skymap.delete_array();  //Gulli20070810
	  galaxy.pi0_decay_HII_skymap.delete_array();  //IMOS20080114
	  
	}
	
	galaxy.pi0_decay_skymap.delete_array();  //Gulli20070810
	
      }
	
    }
    
    if (galdef.primary_electrons || 
	galdef.secondary_electrons || 
	galdef.knock_on_electrons || // IMOS20060504
	galdef.secondary_positrons || 
	galdef.DM_positrons || 
	galdef.DM_electrons) {        // IMOS20050912 
	
      if (galdef.bremss) {                               // IMOS20060420
	    
	//No need to read if already read by pi0_decay   Gulli20070810
	
	if (!galdef.pi0_decay) {

	  //read_HIR();            //AWS20041213
	  //read_COR();            //AWS20041213
	  TIME_SUBROUTINE(read_gas_maps, hir);                // IMOS20080114
	  TIME_SUBROUTINE(read_gas_maps, cor);                // IMOS20080114

	}
	
	//Moved from create_galaxy    Gulli20070810
	if (2 == galdef.n_spatial_dimensions) {
	
	  galaxy.bremss_emiss.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);
	  
	  galaxy.bremss_ionized_emiss.init(galaxy.n_rgrid,galaxy.n_zgrid,galaxy.n_E_gammagrid);
	  
	}
	 
	if (3 == galdef.n_spatial_dimensions) {
	
	  galaxy.bremss_emiss.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);
	  
	  galaxy.bremss_ionized_emiss.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);
	  
	}
	
	if (3 == galdef.skymap_format) {

	  galaxy.bremss_hp_skymap.Resize(galdef.healpix_order, galaxy.n_E_gammagrid);
	  
	  galaxy.bremss_ionized_hp_skymap.Resize(galdef.healpix_order, galaxy.n_E_gammagrid);

		
	  if (2 == galdef.gamma_rays) {

	    galaxy.bremss_H2R_hp_skymap = new Skymap<double>[galaxy.hpHIR.nSpectra()];
		
	    galaxy.bremss_HIR_hp_skymap = new Skymap<double>[galaxy.hpHIR.nSpectra()];
		
	    galaxy.bremss_HII_hp_skymap = new Skymap<double>[galaxy.hpHIR.nSpectra()];// IMOS20080114*

		
	    for (int i_Ring=0; i_Ring<galaxy.hpHIR.nSpectra(); ++i_Ring) { 
		
	      galaxy.bremss_H2R_hp_skymap[i_Ring].Resize(galdef.healpix_order, galaxy.n_E_gammagrid);
		
	      galaxy.bremss_HIR_hp_skymap[i_Ring].Resize(galdef.healpix_order, galaxy.n_E_gammagrid);
		
	      galaxy.bremss_HII_hp_skymap[i_Ring].Resize(galdef.healpix_order, galaxy.n_E_gammagrid);// IMOS20080114*
		
	    }
		
	  }
	  
	} else {
		
	  galaxy.bremss_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid);
		
	  galaxy.bremss_ionized_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid);

		
	  if (2 == galdef.gamma_rays) {

	    galaxy.bremss_H2R_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.COR.n_zgrid, galaxy.n_E_gammagrid);
		
	    galaxy.bremss_HIR_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.HIR.n_zgrid, galaxy.n_E_gammagrid);
		
	    galaxy.bremss_HII_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.HIR.n_zgrid, galaxy.n_E_gammagrid);// IMOS20080114*
		
	  }
	  
	}
	      
        TIME_FUNCTION(stat,gen_bremss_emiss);
	if (stat != 0) return 1;
        TIME_FUNCTION(stat,store_bremss_emiss);
	if (stat != 0) return 1;
        TIME_FUNCTION(stat,gen_bremss_skymap);
	if (stat != 0) return 1;
        TIME_FUNCTION(stat,store_bremss_skymap);
	if (stat != 0)  return 1;
	      
	if (2 == galdef.gamma_rays) {                    //AWS20050302
		
	  TIME_FUNCTION(stat,store_bremss_H2R_skymap);
	  if (stat != 0) return 1;//AWS20041215
	  TIME_FUNCTION(stat,store_bremss_HIR_skymap);
	  if (stat != 0) return 1;//AWS20041215
	  TIME_FUNCTION(stat,store_bremss_HII_skymap);
	  if (stat != 0) return 1;// IMOS20080114*
		
	}
	
	TIME_FUNCTION(stat,gen_bremss_ionized_skymap);
	if (stat != 0)  return 1;
	TIME_FUNCTION(stat,store_bremss_ionized_skymap);
	if (stat != 0)  return 1;

	if (3 == galdef.skymap_format) {

	  galaxy.bremss_hp_skymap.Resize(0,0);
	  galaxy.bremss_ionized_hp_skymap.Resize(0,0);
	  
	  if (2 == galdef.gamma_rays) {

	    delete[] galaxy.bremss_H2R_hp_skymap;
	    galaxy.bremss_H2R_hp_skymap = 0;
	    delete[] galaxy.bremss_HIR_hp_skymap;
	    galaxy.bremss_HIR_hp_skymap = 0;
	    delete[] galaxy.bremss_HII_hp_skymap;// IMOS20080114*
	    galaxy.bremss_HII_hp_skymap = 0;
		
	  }
	
	} else {
		
	  galaxy.bremss_emiss.delete_array();           //Gulli20070810
	  galaxy.bremss_ionized_emiss.delete_array();   //Gulli20070810
	  galaxy.bremss_skymap.delete_array();          //Gulli20070810
	  galaxy.bremss_ionized_skymap.delete_array();  //Gulli20070810
		
	  if (2 == galdef.gamma_rays) {

	    galaxy.bremss_H2R_skymap.delete_array();  //Gulli20070810
	    galaxy.bremss_HIR_skymap.delete_array();  //Gulli20070810
	    galaxy.bremss_HII_skymap.delete_array();  // IMOS20080114*
		
	  }
	  
	}
	
      }
      
      if (galdef.IC_isotropic+galdef.IC_anisotropic) { // IMOS20060420
	    
	//Moved from create_galaxy  //Gulli20070810
	galaxy.IC_iso_emiss = new Distribution[galaxy.n_ISRF_components];
	galaxy.IC_aniso_emiss = new Distribution;

	if (2 == galdef.n_spatial_dimensions) {
	
	  for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {
		
	    galaxy.IC_iso_emiss[icomp].init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);
	  
	  }

	  if (galdef.IC_anisotropic) { 

	     galaxy.IC_aniso_emiss->init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);

	  }
		
	}
	
	if (3 == galdef.n_spatial_dimensions) {
	
	  for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {
		
	    galaxy.IC_iso_emiss[icomp].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);
		
	  }
	  
	  if (galdef.IC_anisotropic) { 

	     galaxy.IC_aniso_emiss->init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);

	  }

	}

	if (3 == galdef.skymap_format) {

	  galaxy.IC_iso_hp_skymap = new Skymap<double>[galaxy.n_ISRF_components];
	  for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {
		
	    galaxy.IC_iso_hp_skymap[icomp].Resize(galdef.healpix_order, galaxy.n_E_gammagrid);
		
	  }

	  if (galdef.IC_anisotropic) { // IMOS20060420
		
	    galaxy.IC_aniso_hp_skymap = new Skymap<double>[galaxy.n_ISRF_components];
		
	    for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {
		
	      galaxy.IC_aniso_hp_skymap[icomp].Resize(galdef.healpix_order, galaxy.n_E_gammagrid);
		
	    }
		
	  }
	  
	} else {
		
	  galaxy.IC_iso_skymap = new Distribution[galaxy.n_ISRF_components];
		
	  for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {
		
	    galaxy.IC_iso_skymap[icomp].init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid);
		
	  }

	  if (galdef.IC_anisotropic) { // IMOS20060420
		
	    galaxy.IC_aniso_skymap = new Distribution[galaxy.n_ISRF_components];
		
	    for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {
		
	      galaxy.IC_aniso_skymap[icomp].init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid);
		
	    }
		
	  } //End move Gulli20070810
	  
	}

	bool isGood;
	TIME_FUNCTION(isGood, gen_IC_emiss); // The previous check was a bit strange because it returned true if there were no electrons, then returned from the calculation which causes memory leaks -- TAP 20090311

	if (isGood) {
	  
	  string icType = "isotropic";
	  	  
	  TIME_SUBROUTINE(store_IC_emiss, icType);
	  
	  TIME_SUBROUTINE(gen_IC_skymap);
	
	  TIME_SUBROUTINE(store_IC_skymap_comp, icType);
	  TIME_SUBROUTINE(store_IC_skymap, icType); //AWS20090415
  
	  if (galdef.IC_anisotropic) {

	    icType = "anisotropic";

	    TIME_SUBROUTINE(store_IC_skymap_comp, icType);
	    TIME_SUBROUTINE(store_IC_skymap, icType); 
 
	  }

	}

	//if (               gen_IC_emiss() !=0)  return 1;
	//if (               gen_IC_skymap()!=0)  return 1;
	//if (store_IC_skymap("isotropic")!=0)  return 1;
	//if(galdef.IC_anisotropic) if(store_IC_skymap("anisotropic")!=0) return 1; // AWS20010206
	
	//Delete arrays Gulli20070810
	
	//if (store_IC_skymap_comp("isotropic")!=0)  return 1;
	//    if(galdef.IC_anisotropic) if(store_IC_skymap_comp("anisotropic")!=0) return 1;

	  
	for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {
	
	  galaxy.IC_iso_emiss[icomp].delete_array();
	  
	}

	if (3 == galdef.skymap_format) {
	
	  delete[] galaxy.IC_iso_hp_skymap;
	  galaxy.IC_iso_hp_skymap = 0;
	  delete[] galaxy.IC_aniso_hp_skymap;
	  galaxy.IC_aniso_hp_skymap = 0;
	  
	} else {
	  
	  for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {
	    
	    galaxy.IC_iso_skymap[icomp].delete_array();
		
	  }
		
	  if (galdef.IC_anisotropic) {
		
	    for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp) {
		
	      galaxy.IC_aniso_skymap[icomp].delete_array();
		
	    }
		
	  }
	  
	}
	
      }
	
    } //electrons
      //Cleaning up COR and HIR gas surveys  Gulli20070810
    if (galdef.pi0_decay || galdef.bremss) {
	
      galaxy.COR.delete_array();
      galaxy.HIR.delete_array();
      
    }
    
  } //galdef.gamma_rays
  
  if (galdef.DM_gammas) {                           // IMOS20050912
    
    //Moved from create_galaxy   Gulli20070810
    if (2 == galdef.n_spatial_dimensions) {
      
      galaxy.DM_emiss.init(galaxy.n_rgrid, galaxy.n_zgrid, galaxy.n_E_gammagrid); 
    }
    
    if (3 == galdef.n_spatial_dimensions) {
	
      galaxy.DM_emiss.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, galaxy.n_E_gammagrid);
      
    }
    
    if (3 == galdef.skymap_format) {
	
      galaxy.DM_hp_skymap.Resize(galdef.healpix_order, galaxy.n_E_gammagrid);
      
    } else {

      galaxy.DM_skymap.init(galaxy.n_long, galaxy.n_lat, galaxy.n_E_gammagrid); // DM: IMOS20050912 IMOS20060420
    }
    
    TIME_FUNCTION(stat,gen_DM_emiss);
    if (stat != 0) return 1;
    TIME_FUNCTION(stat,store_DM_emiss);
    if (stat != 0) return 1;
    TIME_FUNCTION(stat,gen_DM_skymap);
    if (stat != 0) return 1;
    TIME_FUNCTION(stat,store_DM_skymap);
    if (stat != 0) return 1;
    
    if (3 == galdef.skymap_format) {
      
      galaxy.DM_hp_skymap.Resize(0,0);
      
    } else {
	
      galaxy.DM_skymap.delete_array();//Gulli20070810
      
    }
    
    galaxy.DM_emiss.delete_array();  //Gulli20070810
    
  }
  
  if (galdef.ionization_rate) { // IMOS20060420
    
    //Moved from create_galaxy  Gulli20071003
    cout<<"   creating ionization rate array"<<endl;
    
    if (2 == galdef.n_spatial_dimensions) {
      
      galaxy.ionization_rate.init(galaxy.n_rgrid, galaxy.n_zgrid,1);
      
    }
    
    if (3 == galdef.n_spatial_dimensions) {
	
      galaxy.ionization_rate.init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, 1);
      
    }
      
    if (gen_ionization_rate() != 0) return 1;
    if (store_ionization_rate() != 0) return 1;
      
    //Delete Arrays  Gulli20071003
    
    galaxy.ionization_rate.delete_array();
    
  }

  ostringstream bufProc;
  bufProc << "Completed processing galdef_" << configure.fVersion << "_" << runNumber;
  INFO(bufProc.str());

  INFO("Exit: Galprop main procedure");
  
  //cout<<" completed processing of galdef_"<<version<<"_"<<argv[1]<<endl;
  //cout<<" <<<< galprop"<<endl;
  
  return 0;

}
