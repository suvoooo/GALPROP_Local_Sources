#include"Galprop.h"
#include"galprop_classes.h"
#include"galprop_internal.h"
#include"constants.h"
#include"fitsio.h"

#include <iostream>
#include <sstream>
#include <vector>

#include <ErrorLogger.h>

// interfaces to nH.cc
double nHI (double Rkpc, double Zkpc);
double nH2 (double Rkpc, double Zkpc);
double nH2 (double Rkpc, double Zkpc,double H2toCO);
double nHII(double Rkpc, double Zkpc);

// compute luminosity of Galaxy and other integral properties like gas mass

int Galprop::gen_luminosity() {

  INFO("Entry");

  valarray<double> gamma_bremss_luminosity_spectrum   (galaxy.n_E_gammagrid); 
  valarray<double> gamma_pi0_decay_luminosity_spectrum(galaxy.n_E_gammagrid); 
  valarray<double> gamma_IC_luminosity_spectrum       (galaxy.n_E_gammagrid); 
  valarray<double> gamma_total_luminosity_spectrum    (galaxy.n_E_gammagrid); 

  valarray<double> gamma_total_integral_luminosity    (galaxy.n_E_gammagrid); 
  valarray<double> gamma_bremss_integral_luminosity   (galaxy.n_E_gammagrid); 
  valarray<double> gamma_pi0_decay_integral_luminosity(galaxy.n_E_gammagrid); 
  valarray<double> gamma_IC_integral_luminosity       (galaxy.n_E_gammagrid); 

  valarray<double> gamma_total_integral_photons       (galaxy.n_E_gammagrid); 
  valarray<double> gamma_bremss_integral_photons      (galaxy.n_E_gammagrid); 
  valarray<double> gamma_pi0_decay_integral_photons   (galaxy.n_E_gammagrid); 
  valarray<double> gamma_IC_integral_photons          (galaxy.n_E_gammagrid); 

  valarray<double> synchrotron_luminosity_spectrum    (galaxy.n_nu_synchgrid); 
  valarray<double> synchrotron_integral_luminosity    (galaxy.n_nu_synchgrid); 

  valarray<double> nHI_column_density                 (galaxy.n_rgrid);
  valarray<double> nH2_column_density                 (galaxy.n_rgrid);
  valarray<double> nHII_column_density                (galaxy.n_rgrid);

  gamma_bremss_luminosity_spectrum = 0.;
  gamma_pi0_decay_luminosity_spectrum = 0.;
  gamma_IC_luminosity_spectrum = 0.;

  synchrotron_luminosity_spectrum = 0.;

  nHI_column_density  = 0.;
  nH2_column_density  = 0.;
  nHII_column_density = 0.;

  double nHI_total =0.;
  double nH2_total =0.;
  double nHII_total=0.;
  double nH_total  =0.;
  double volume_total = 0.;

  double pi=acos(-1.);

  const int irmax = galaxy.n_rgrid-1; // NB last point
  //  irmax=6; // for testing

  ostringstream maxRbuf;

  maxRbuf << "Maximum radius for luminosity = " << galaxy.r[irmax];

  INFO(maxRbuf.str());

  if (2 == galaxy.n_spatial_dimensions) {

    for (int ir = 0; ir < irmax; ++ir) {

      for (int iz = 0; iz < galaxy.n_zgrid; ++iz) {
	
	const double volume = pi*(pow(galaxy.r[ir+1], 2.) - pow(galaxy.r[ir], 2.))*galaxy.dz*pow(kpc2cm, 3.) ; // not exact of course since emissivity is at r[ir]. volume in kpc^3 -> cm^3
	//volume *= pow(kpc2cm, 3.); // from constants.h   volume now in cm^3

	const double nHI_ = nHI(galaxy.r[ir], galaxy.z[iz]);
	const double nH2_ = nH2(galaxy.r[ir], galaxy.z[iz], fX_CO(galaxy.r[ir]))*2.0; // Gulli20100218 Using fX_CO for X ratio. NB two atoms per molecule!
	const double nHII_ = nHII(galaxy.r[ir], galaxy.z[iz]);
	const double nH = nHI_ + nH2_ + nHII_;
     
	nHI_total += nHI_*volume; 
	nH2_total += nH2_*volume;
	nHII_total += nHII_*volume;
	nH_total += nH*volume;

	nHI_column_density[ir] += nHI_*galaxy.dz*kpc2cm;
	nH2_column_density[ir] += nH2_*galaxy.dz*kpc2cm;
	nHII_column_density[ir] += nHII_*galaxy.dz*kpc2cm;

	volume_total += volume;	

	const double volumeFactor = volume*4.*pi, gasFactor = volumeFactor*nH;

	for (int ip = 0; ip < galaxy.n_E_gammagrid; ++ip) {
	  
	  // emissivity is per sr so need 4pi factor
	  gamma_bremss_luminosity_spectrum[ip] += galaxy.bremss_emiss.d2[ir][iz].s[ip]*gasFactor;
	  gamma_pi0_decay_luminosity_spectrum[ip] += galaxy.pi0_decay_emiss.d2[ir][iz].s[ip]*gasFactor;
	  
	  for (int icomp = 0; icomp < galaxy.n_ISRF_components; ++icomp)
	    gamma_IC_luminosity_spectrum[ip] += galaxy.IC_iso_emiss[icomp].d2[ir][iz].s[ip]*volumeFactor;

	}

	for (int inu = 0; inu < galaxy.n_nu_synchgrid; ++inu) {
	  
	  synchrotron_luminosity_spectrum[inu] += galaxy.synchrotron_emiss.d2[ir][iz].s[inu]*volumeFactor;
	
	}
	
      }

    }

    gamma_total_luminosity_spectrum  = gamma_bremss_luminosity_spectrum + gamma_pi0_decay_luminosity_spectrum + gamma_IC_luminosity_spectrum;
    
    // since emissivity is MeV^2 cm-3 sr-1 MeV-1 s-1 
    // the units are  MeV s-1  as in Strong etal 2000 ApJ 537, 763 Fig 20.
    
    gamma_total_integral_luminosity = gamma_bremss_integral_luminosity = 
      gamma_pi0_decay_integral_luminosity = gamma_IC_integral_luminosity = 
      gamma_total_integral_photons = gamma_bremss_integral_photons = 
      gamma_pi0_decay_integral_photons = gamma_IC_integral_photons = 0.;
    
    for (int ip = 0; ip < galaxy.n_E_gammagrid; ++ip)
      for (int ipp = ip; ipp < galaxy.n_E_gammagrid; ++ipp) {
	
	gamma_total_integral_luminosity[ip] += gamma_total_luminosity_spectrum[ipp];
	gamma_bremss_integral_luminosity[ip] += gamma_bremss_luminosity_spectrum[ipp];
	gamma_pi0_decay_integral_luminosity [ip] += gamma_pi0_decay_luminosity_spectrum[ipp];
	gamma_IC_integral_luminosity[ip] += gamma_IC_luminosity_spectrum[ipp];
	
	gamma_total_integral_photons[ip] += gamma_total_luminosity_spectrum[ipp]/galaxy.E_gamma[ipp];  // *E instead of *E^2
	gamma_bremss_integral_photons[ip] += gamma_bremss_luminosity_spectrum[ipp]/galaxy.E_gamma[ipp];
	gamma_pi0_decay_integral_photons[ip] += gamma_pi0_decay_luminosity_spectrum[ipp]/galaxy.E_gamma[ipp];
	gamma_IC_integral_photons[ip] += gamma_IC_luminosity_spectrum[ipp]/galaxy.E_gamma[ipp];
	
      }
    
    // integral Ef(E)dE = integral E^2f(E) dlog(E) = sum E^2 f(E) delta logE.  delta logE = log(E factor)
    gamma_total_integral_luminosity *= log(galaxy.E_gamma_factor);
    gamma_total_integral_luminosity *= MEV2ERG;// from constants.h
    
    gamma_bremss_integral_luminosity *= log(galaxy.E_gamma_factor);
    gamma_bremss_integral_luminosity *= MEV2ERG;// from constants.h
    
    gamma_pi0_decay_integral_luminosity *= log(galaxy.E_gamma_factor);
    gamma_pi0_decay_integral_luminosity *= MEV2ERG;// from constants.h
    
    gamma_IC_integral_luminosity *= log(galaxy.E_gamma_factor);
    gamma_IC_integral_luminosity *= MEV2ERG;// from constants.h
    
    // integral  f(E)dE = integral Ef(E) dlog(E) = sum E f(E) delta logE.  delta logE = log(E factor)
    gamma_total_integral_photons *= log(galaxy.E_gamma_factor);
    gamma_bremss_integral_photons *= log(galaxy.E_gamma_factor);
    gamma_pi0_decay_integral_photons *= log(galaxy.E_gamma_factor);
    gamma_IC_integral_photons *= log(galaxy.E_gamma_factor);
    
    synchrotron_integral_luminosity = 0.;
    
    for (int inu = 0; inu < galaxy.n_nu_synchgrid; ++inu)
      for (int iinu = inu; iinu < galaxy.n_nu_synchgrid; ++iinu)
	synchrotron_integral_luminosity[inu] += synchrotron_luminosity_spectrum[iinu]*galaxy.nu_synch[iinu];
    
    synchrotron_integral_luminosity *= log(galaxy.nu_synch_factor);
    
    // ======================== print out the results =====================================
    
    // integral properties of Galaxy
    
    ostringstream gasBuf1, gasBuf2, volBuf;
    
    volBuf << "Total galactic volume: " << volume_total << " cm^3";
    
    INFO(volBuf.str());
    
    gasBuf1 << "Total gas atoms: " 
	    << nHI_total << " (HI) " 
	    << nH2_total << " (H2) " 
	    << nHII_total << " (HII) " 
	    << nH_total << " (total)";
    gasBuf2 << "Total gas mass (Msun): " 
	    << nHI_total*1.67e-24/2.0e33 << " (HI) " 
	    << nH2_total*1.67e-24/2.0e33 << " (H2) " 
	    << nHII_total*1.67e-24/2.0e33 << " (HII) " 
	    << nH_total*1.67e-24/2.0e33 << " (total) ";
    
    INFO(gasBuf1.str());
    INFO(gasBuf2.str());
    
    const double sigmaFactor = 1.67e-24/2.0e33*pow(kpc2cm, 2.)*1.0e-6; // atoms cm-2 to Msun/pc^2
    
    for (int ir = 0; ir < galaxy.n_rgrid-1; ++ir) {
      
      ostringstream numBuf, massBuf;

      numBuf << "z-column density: R = " << galaxy.r[ir] << " kpc, " 
	     << nHI_column_density[ir] << " (HI atoms cm^-2) "
	     << nH2_column_density[ir] << " (H2 atmos cm^-2) "
	     << nHII_column_density[ir] << " (HII atoms cm^-2)";
      
      massBuf << "z-column density: R = " << galaxy.r[ir] << " kpc, " 
	      << nHI_column_density[ir]*sigmaFactor << " (HI Msun pc^-2) "
	      << nH2_column_density[ir]*sigmaFactor << " (H2 Msun pc^-2) "
	      << nHII_column_density[ir]*sigmaFactor << " (HII Msun pc^-2)";
      
      INFO(numBuf.str());
      INFO(massBuf.str());
      
    }
    
    //cout<<endl<<" gamma-ray luminosity of Galaxy "<<endl<<endl;
    
    //for(int ip=0;ip<galaxy.n_E_gammagrid;ip++)
    //{
    //int iz=galaxy.n_zgrid/2;
    //cout<<"pi0-decay emissivity spectrum at z=0: "
    //    <<  "E="               <<galaxy.E_gamma[ip]<<" "; 
    //for(int ir=0;ir<galaxy.n_rgrid;ir++)        cout   <<galaxy.pi0_decay_emiss    .d2[ir][iz].s[ip]<<" ";
    // cout   << endl;
    //}
    
    ostringstream gammaBuf;
    gammaBuf << "Gamma-ray luminosity (erg s^-1) and photons (ph s^-1): " << endl 
	     << setw(5) << "Energy"
	     << setw(5) << "Bremss"
	     << setw(5) << "pi0"
	     << setw(5) << "IC"
	     << setw(5) << "Total"
	     << setw(5) << "Integral bremss"
	     << setw(5) << "Integral pi0"
	     << setw(5) << "Integral IC"
	     << setw(5) << "Integral Total"
	     << setw(5) << "Photons bremss"
	     << setw(5) << "Photons pi0"
	     << setw(5) << "Photons IC"
	     << setw(5) << "Photons Total" << endl;
    
    for (int ip = 0; ip < galaxy.n_E_gammagrid; ++ip)
      gammaBuf << setw(5) << galaxy.E_gamma[ip] << " "
	       << setw(5) << gamma_bremss_luminosity_spectrum[ip] << " "
	       << setw(5) << gamma_pi0_decay_luminosity_spectrum[ip] << " "
	       << setw(5) << gamma_IC_luminosity_spectrum[ip] << " "
	       << setw(5) << gamma_total_luminosity_spectrum[ip] << " "
	       << setw(5) << gamma_bremss_integral_luminosity[ip] << " "
	       << setw(5) << gamma_pi0_decay_integral_luminosity[ip] << " "
	       << setw(5) << gamma_IC_integral_luminosity[ip] << " "
	       << setw(5) << gamma_total_integral_luminosity[ip] << " "
	       << setw(5) << gamma_bremss_integral_photons[ip] << " "
	       << setw(5) << gamma_pi0_decay_integral_photons[ip] << " "
	       << setw(5) << gamma_IC_integral_photons[ip] << " "
	       << setw(5) << gamma_total_integral_photons[ip] << endl;
    
    
    //gammaBuf << "Energy = " << galaxy.E_gamma[ip] << " "  
    //	 << "bremss: " << gamma_bremss_luminosity_spectrum[ip] << " "
    //	 << "pi0: " << gamma_pi0_decay_luminosity_spectrum[ip] << " "
    //	 << "IC: " << gamma_IC_luminosity_spectrum[ip] << " "
    //	 << "total: " << gamma_total_luminosity_spectrum[ip] << " "
    //	 << "integral luminosity erg s-1 :"
    //	 << " bremss: "<< gamma_bremss_integral_luminosity   [ip]<<" "
    //	 <<  " pi0: "  <<gamma_pi0_decay_integral_luminosity [ip]<<" "
    //	 <<  " IC: "   <<gamma_IC_integral_luminosity        [ip]<<" "
    //	 <<  " total: "<<gamma_total_integral_luminosity     [ip]<<" "
    //	 << " integral photons s-1:"            
    //	 << " bremss: "<< gamma_bremss_integral_photons      [ip]<<" "
    //	 <<  " pi0: "  <<gamma_pi0_decay_integral_photons    [ip]<<" "
    //	 <<  " IC: "   <<gamma_IC_integral_photons           [ip]<<" "
    //	 <<  " total: "<<gamma_total_integral_photons        [ip]<<" "
    //<<endl;
    
    INFO(gammaBuf.str());
    
    ostringstream syncBuf;
    
    syncBuf << "Synchrotron luminosity (erg Hz^-1 s^-1): " << endl;
    syncBuf << setw(5) << "Frequency" 
	    << setw(5) << "Differential" 
	    << setw(5) << "Integral" << endl;
    
    for (int iNu = 0; iNu < galaxy.n_nu_synchgrid; ++iNu)
      syncBuf << setw(5) << galaxy.nu_synch[iNu] << " "
	      << setw(5) << synchrotron_luminosity_spectrum[iNu] << " " 
	      << setw(5) << synchrotron_integral_luminosity[iNu] << endl;
    
    INFO(syncBuf.str());
    
    //for (int inu = 0; inu < galaxy.n_nu_synchgrid; ++inu) {
    //cout<<"nu ="<<galaxy.nu_synch[inu]<<" synchrotron luminosity erg Hz-1 s-1: "<< synchrotron_luminosity_spectrum[inu]
    // <<" integral luminosity erg s-1: "<< synchrotron_integral_luminosity[inu]
    //<<endl;
    //}
    
  }
  
  INFO("Exit");

  //cout<<" <<<< gen_luminosity   "<<endl;
  return 0;
 
}
