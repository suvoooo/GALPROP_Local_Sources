
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * gen_bremss_skymap.cc *                     galprop package * 4/29/2002 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <cassert>
#include <string>
#include <cstring>
#include <vector>

using namespace std;//AWS20050624

#include "galprop_classes.h"
#include "galprop_internal.h"

#include <ErrorLogger.h>

// list of selectable debugs:  IMOS20080114
// galdef.verbose==-453 tests the line of sight integration; output skymaps HI and H2 should 
//                      match the distribution of HI and CO gas
// galdef.verbose==-454 also tests the line of sight integration in another place; 
//                      output skymaps HI and H2 should match the distribution of HI and CO gas
// galdef.verbose==-457 is the return to the old method of integration to compare with older versions (does not work for HII)


int Galprop::gen_bremss_skymap() {
  
  INFO("Entry");

  int stat=0;
  const double dtr=acos(-1.)/180.; // conversion degrees to radians
  
  //This part is Gulli's output in HEALPIX format, re-written by IMOS20080114
  if (3 == galdef.skymap_format) {

#pragma omp parallel for schedule(dynamic) default(shared) 
    for (int ii = 0; ii < galaxy.bremss_hp_skymap.Npix(); ++ii)
      {
	SM::Coordinate co(galaxy.bremss_hp_skymap.pix2ang(ii));
	double l=co.l();
	double b=co.b();
	vector< vector<double> > emiss_av; // emissivity * gas density averaged for each ring
	vector< vector<double> > emiss_HII;// emissivity * HII density averaged for each ring
	vector< vector<double> > emiss_HI ;// emissivity * HI  density averaged for each ring
	vector< vector<double> > emiss_H2 ;// emissivity * H2  density averaged for each ring
	
	vector<double> n_H_av;   //              gas density averaged for each ring
	vector<double> n_HI_av;  //         NHI  gas density averaged for each ring   
	vector<double> n_H2_av;  //         NH2  gas density averaged for each ring   
	vector<double> n_CO_av;  //         CO   density averaged for each ring   
	
	// calc. of a 2D array of g-ray emission for the given E_gammagrid and the for all rings for a particular pixel (l,b);
	// also calculates a vector of the gas column density for all rings for a particular pixel (l,b) 
	gen_bremss_skymap_pixel(l,b, n_H_av, n_HI_av, n_H2_av, n_CO_av, emiss_av, emiss_HII, emiss_HI, emiss_H2, galaxy.hpHIR.nSpectra());
	
	for (int i_Ring=0; i_Ring<galaxy.hpHIR.nSpectra(); i_Ring++) 
	  {
	    double w_HI, w_H2;
	    
	    //radius corresponding to the middle of the current ring  
	    double R_ring = (galaxy.R_bins[i_Ring] +galaxy.R_bins[galaxy.n_Ring+i_Ring])/2.;
	    double av_X_CO;
  	    if ( n_CO_av[i_Ring] != 0) {
	       av_X_CO = n_H2_av[i_Ring]/n_CO_av[i_Ring];
	    }else{
	       av_X_CO = fX_CO(R_ring);
	    }
	    
	    if( !(n_H_av[i_Ring] +emiss_HII[i_Ring][galaxy.n_E_gammagrid/2]) ) continue; // HII emissivity is added in case if there is HII gas
	    
	    w_HI =                 galaxy.hpHIR[co][i_Ring] /n_HI_av[i_Ring];
	    w_H2 = 2*av_X_CO*galaxy.hpCOR[co][i_Ring] /n_H2_av[i_Ring];
	    
	    // calc. of the g-ray spectrum for a pixel on the skymap
	    for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++) 
	      {
		if(galdef.verbose==-454) // selectable debug to test vs gas maps IMOS20080114
		  {
		    emiss_av [i_Ring][iEgamma] =n_H_av [i_Ring];
		    emiss_HI [i_Ring][iEgamma] =n_HI_av[i_Ring] *1.0e-20;
		    emiss_H2 [i_Ring][iEgamma] =n_H2_av[i_Ring] /av_X_CO/2.;
		    //			  emiss_HII[i_Ring][iEgamma] =0.; //IMOS20080114*
		  }
		
		galaxy.bremss_hp_skymap[co][iEgamma] += emiss_H2[i_Ring][iEgamma]*w_H2 +emiss_HI[i_Ring][iEgamma]*w_HI +emiss_HII[i_Ring][iEgamma];
		
		if (galdef.gamma_rays==2) // separate skymaps for HI, H2 Gulli20070810
		  {
		    galaxy.bremss_H2R_hp_skymap[i_Ring][co][iEgamma] += emiss_H2 [i_Ring][iEgamma] *w_H2;
		    //			  galaxy.bremss_HIR_hp_skymap[i_Ring][co][iEgamma] += emiss_HI [i_Ring][iEgamma] *w_HI +emiss_HII[i_Ring][iEgamma]; //IMOS20080114*
		    galaxy.bremss_HIR_hp_skymap[i_Ring][co][iEgamma] += emiss_HI [i_Ring][iEgamma] *w_HI; //IMOS20080114*
		    galaxy.bremss_HII_hp_skymap[i_Ring][co][iEgamma] += emiss_HII[i_Ring][iEgamma];       //IMOS20080114*
		  }  
		
		if(galdef.verbose==-100) // selectable debug
		  {
		    cout<<"l b i_Ring w HIR n_H_av skymap "<<l<<" "<<b<<" "<<i_Ring<<" "<<" "
			<<galaxy.hpHIR [co][i_Ring]<<" "<<n_H_av[i_Ring] <<" "<<galaxy.bremss_hp_skymap[co][iEgamma]<<endl;
		    
		    if (galdef.gamma_rays==2) //Gulli20070810
		      {
			cout<<"l b i_Ring w_HI HIR n_HI_av   skymap "<<l<<" "<<b<<" "<<i_Ring<<" "<<w_HI<<" "
			    <<galaxy.hpHIR[co][i_Ring]<<" "<<n_HI_av[i_Ring]<<" "<<galaxy.bremss_HIR_hp_skymap[i_Ring][co][iEgamma]<<endl;
			
			cout<<"l b i_Ring w_H2 COR n_H2_av   skymap "<<l<<" "<<b<<" "<<i_Ring<<" "<<w_H2<<" "
			    <<galaxy.hpCOR[co][i_Ring]<<" "<<n_H2_av[i_Ring]<<" "<<galaxy.bremss_H2R_hp_skymap[i_Ring][co][iEgamma]<<endl;
		      } 
		  }
	      }//iEgamma
	  }//i_Ring
      }//parallel for
    /*	
    //Print out points where n_HI_av is 0 in a ring the pixel spans
    for (map<double,vector<double> >::iterator it = zeroMap.begin(); it != zeroMap.end(); ++it)
    {
    for (int i = 0; i < it->second.size(); ++i)
    {
    cout<<it->first<<"  "<<it->second[i]<<endl;
    }  
    }
    */
    
    //cout<<"Set spectra"<<endl;
    
    galaxy.bremss_hp_skymap.setSpectra(galaxy.E_gamma, galaxy.n_E_gammagrid);
    
    if (2 == galdef.gamma_rays)  // separate skymaps for HI, H2 Gulli20070810
      {
	for (int i_Ring=0; i_Ring<galaxy.hpHIR.nSpectra(); ++i_Ring) 
	  { 
	    galaxy.bremss_H2R_hp_skymap[i_Ring].setSpectra(galaxy.E_gamma, galaxy.n_E_gammagrid);
	    galaxy.bremss_HIR_hp_skymap[i_Ring].setSpectra(galaxy.E_gamma, galaxy.n_E_gammagrid);
	    galaxy.bremss_HII_hp_skymap[i_Ring].setSpectra(galaxy.E_gamma, galaxy.n_E_gammagrid);//IMOS20080114*
	    }
      }
    
    if(galdef.verbose==-101)        // selectable debug
      {    
	cout<<" bremss skymap "<<endl;
	galaxy.bremss_skymap.print();
	
	if (galdef.gamma_rays==2) // separate skymaps for HI, H2 Gulli20070810
	  {
	    for (int i_Ring=0; i_Ring<galaxy.hpHIR.nSpectra(); ++i_Ring) 
	      { 
		cout<<" bremss_H2R skymap, ring: "<<i_Ring<<endl;
		galaxy.bremss_H2R_hp_skymap[i_Ring].print(cout);
		
		cout<<" bremss_HIR skymap, ring: "<<i_Ring<<endl;
		galaxy.bremss_HIR_hp_skymap[i_Ring].print(cout);
		
		cout<<" bremss_HII skymap, ring: "<<i_Ring<<endl; //IMOS20080114*
		galaxy.bremss_HII_hp_skymap[i_Ring].print(cout);  //IMOS20080114*
	      }
	  }
	
      }//galdef.verbose>=2
  }
  
  else { // standard (non-HEALPIX) format IMOS20080114
#pragma omp parallel for schedule(dynamic) default(shared) 
    
    for(int i_long=0; i_long<galaxy.n_long; i_long++) {  
      for(int i_lat =0; i_lat <galaxy.n_lat; i_lat++) {	
	double l=galaxy.long_min + i_long*galaxy.d_long;
	double b=galaxy. lat_min + i_lat *galaxy.d_lat ;
	vector< vector<double> > emiss_av; // emissivity * gas density averaged for each ring
	vector< vector<double> > emiss_HII;// emissivity * HII density averaged for each ring
	vector< vector<double> > emiss_HI ;// emissivity * HI  density averaged for each ring 
	vector< vector<double> > emiss_H2 ;// emissivity * H2  density averaged for each ring 
	
	vector<double> n_H_av;   //              gas density averaged for each ring
	vector<double> n_HI_av;  //         NHI  gas density averaged for each ring 
	vector<double> n_H2_av;  //         NH2  gas density averaged for each ring
	vector<double> n_CO_av;  //         CO   density averaged for each ring   
	
	// calc. of a 2D array of g-ray emission for the given E_gammagrid and the for all rings for a particular pixel (l,b);
	// also calculates a vector of the gas column density for all rings for a particular pixel (l,b) 
	gen_bremss_skymap_pixel(l,b, n_H_av, n_HI_av, n_H2_av, n_CO_av, emiss_av, emiss_HII, emiss_HI, emiss_H2, galaxy.n_Ring);
	
	for(int i_Ring=0; i_Ring<galaxy.n_Ring; i_Ring++) {
	  double w, w_HI, w_H2;
	  
	  //radius corresponding to the middle of the current ring  
	  double R_ring = (galaxy.R_bins[i_Ring] +galaxy.R_bins[galaxy.n_Ring+i_Ring])/2.;
	  double av_X_CO;
	  if ( n_CO_av[i_Ring] != 0) {
	     av_X_CO = n_H2_av[i_Ring]/n_CO_av[i_Ring];
	  }else{
	     av_X_CO = fX_CO(R_ring);
	  }
	  
	  if(galdef.verbose==-457 && n_H_av[i_Ring]<=0.) continue; // selectable debug (-457 -old method)
	  
	  //	  if(galdef.verbose==-455)                                 // selectable debug to test vs gas maps
	  //	    {
	  //	      if( !(galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0] +2.*fX_CO(R_ring)*galaxy.COR.d3[i_long][i_lat][i_Ring].s[0] +emiss_HII[i_Ring][galaxy.n_E_gammagrid/2]) ) continue;
	  //	    }
	  //	  else
	  if( !(n_H_av[i_Ring] +emiss_HII[i_Ring][galaxy.n_E_gammagrid/2]) ) continue; // HII emissivity is added in case if there is HII gas
	  
	  // calculation of the scaling factors w= column density from gas maps[iRing]/analytical column density[iRing] 
	  w_HI =                  galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]/n_HI_av[i_Ring];
	  
	  if(galdef.verbose==-457)  // selectable debug                                                         
	    {	  
	      if(n_H2_av[i_Ring]>0.) w_H2 = 2.*av_X_CO*galaxy.COR.d3[i_long][i_lat][i_Ring].s[0]/n_H2_av[i_Ring];
	    }
	  
	  else
	    w_H2 = 2.*av_X_CO*galaxy.COR.d3[i_long][i_lat][i_Ring].s[0]/n_H2_av[i_Ring];
	  
	  if(galdef.verbose==-457)  // selectable debug
	    w = (galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0] +2.*av_X_CO*galaxy.COR.d3[i_long][i_lat][i_Ring].s[0])/n_H_av[i_Ring];
	  
	  // calc. of the g-ray spectrum for a pixel on the skymap
	  for(int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++) 
	    {	      
	      if(galdef.verbose==-454) // selectable debug to test vs gas maps IMOS20080114
		{
		  emiss_av[i_Ring][iEgamma]  =n_H_av [i_Ring]                  /pow(galaxy.E_gamma[iEgamma],2);
		  emiss_H2[i_Ring][iEgamma]  =n_H2_av[i_Ring] /av_X_CO/2.      /pow(galaxy.E_gamma[iEgamma],2);
		  emiss_HI[i_Ring][iEgamma]  =n_HI_av[i_Ring] *1.0e-20         /pow(galaxy.E_gamma[iEgamma],2);
		  emiss_HII[i_Ring][iEgamma]/=                                  pow(galaxy.E_gamma[iEgamma],2);    //IMOS20080114*
		}
	      
	      if(galdef.verbose==-457)  // selectable debug
		galaxy.bremss_skymap.d2[i_long][i_lat].s[iEgamma] += emiss_av[i_Ring][iEgamma]*w +emiss_HII[i_Ring][iEgamma];
	      
	      else
		galaxy.bremss_skymap.d2[i_long][i_lat].s[iEgamma] += emiss_H2[i_Ring][iEgamma]*w_H2 +emiss_HI[i_Ring][iEgamma]*w_HI +emiss_HII[i_Ring][iEgamma]; 
	      
	      if (galdef.gamma_rays==2) // separate skymaps for HI, H2 Gulli20070810
		{
		  galaxy.bremss_H2R_skymap.d3[i_long][i_lat][i_Ring].s[iEgamma] =emiss_H2[i_Ring][iEgamma]*w_H2;                           
//		  galaxy.bremss_HIR_skymap.d3[i_long][i_lat][i_Ring].s[iEgamma] =emiss_HI[i_Ring][iEgamma]*w_HI+emiss_HII[i_Ring][iEgamma]; //IMOS20080114*
		  galaxy.bremss_HIR_skymap.d3[i_long][i_lat][i_Ring].s[iEgamma] =emiss_HI[i_Ring][iEgamma]*w_HI;                            //IMOS20080114*
		  galaxy.bremss_HII_skymap.d3[i_long][i_lat][i_Ring].s[iEgamma] =emiss_HII[i_Ring][iEgamma];                                //IMOS20080114*
		}  
	      	      
	      if(galdef.verbose==-100) // selectable debug
		{
		  cout<<"l b i_Ring w HIR n_H_av skymap "<<l<<" "<<b<<" "<<i_Ring<<" "<<w<<" "
		      <<galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]<<" "<<n_H_av[i_Ring]
		      <<" "<<galaxy.bremss_skymap.d2[i_long][i_lat].s[iEgamma]<<endl;
		  
		  if (galdef.gamma_rays==2) //Gulli20070810
		    {
		      
		      cout<<"l b i_Ring w_H2 COR n_H2_av   skymap "<<l<<" "<<b<<" "<<i_Ring<<" "<<w_H2<<" "
			  <<galaxy.COR.d3[i_long][i_lat][i_Ring].s[0]<<" "<<n_H2_av[i_Ring]
			  <<" "<<galaxy.bremss_H2R_skymap.d3[i_long][i_lat][i_Ring].s[iEgamma]<<endl;

		      cout<<"l b i_Ring w_HI HIR n_HI_av   skymap "<<l<<" "<<b<<" "<<i_Ring<<" "<<w_HI<<" "
			  <<galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]<<" "<<n_HI_av[i_Ring]
			  <<" "<<galaxy.bremss_HIR_skymap.d3[i_long][i_lat][i_Ring].s[iEgamma]<<endl;

		      cout<<"l b i_Ring n_HII_av skymap "<<l<<" "<<b<<" "<<i_Ring<<" " //IMOS20080114*
			  <<n_HI_av[i_Ring]
			  <<" "<<galaxy.bremss_HII_skymap.d3[i_long][i_lat][i_Ring].s[iEgamma]<<endl;
		    }	    
		}
	    }//iEgamma
	  //}if	  
	}//i_Ring
	//cout<<"i_long i_lat "<<i_long<<" "<<i_lat<<endl;
      }//i_lat
      cout<<"i_long  "<<i_long<<endl;
    }//i_long
    
    if(galdef.verbose==-101)       // selectable debug
      {
	cout<<" bremss skymap "<<endl;
	galaxy.bremss_skymap.print();
	
	if (galdef.gamma_rays==2) //Gulli20070810
	  {
	    cout<<" bremss_H2R skymap "<<endl;
	    galaxy.bremss_H2R_skymap.print();
	    cout<<" bremss_HIR skymap "<<endl;
	    galaxy.bremss_HIR_skymap.print();
	    cout<<" bremss_HII skymap "<<endl; //IMOS20080114*
	    galaxy.bremss_HII_skymap.print();  //IMOS20080114*
	  }    
      }//galdef.verbose>=2
    
  }//Skymap_format
  
  INFO("Exit");

  return stat;
    
}

//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// calc. of a 2D array of g-ray emission for the given E_gammagrid and the
// for all rings for a particular pixel (l,b); also calculates a vector of the 
// gas column density for all rings for a particular pixel (l,b)  IMOS20080114
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
int Galprop::gen_bremss_skymap_pixel(const double l, const double b, vector<double> &n_H_av, vector<double> &n_HI_av, vector<double> &n_H2_av, vector<double> &n_CO_av, vector< vector<double> > &emiss_av, vector< vector<double> > &emiss_HII, vector< vector<double> > &emiss_HI, vector< vector<double> > &emiss_H2, const int n_Ring){

  double dd    =galdef.LoS_step;                           // integration step in kpc 
  double ddd   =galdef.LoS_step /galdef.LoS_substep_number;// finer step to average the gas density (<0 -no finer steps)
  double db    =galdef.d_lat;                              // latitude step
  double ddb   =galdef.d_lat    /galdef.lat_substep_number;// latitude step to average the gas density
  double dtr=acos(-1.)/180.;                               // conversion degrees to radians
  int ir,ix,iy,iz;
  int i_Ring = gas_iRing(Rsun);                            // We always start with the same ring
  double zero=1.e-30;                    //A small number used for the analytical model if it is 0 but the data is non-0

  double w,w_HI,w_H2, nH2_av_lb, nCO_av_lb, nHI_av_lb, nHII_av_lb; 
    
  vector< vector<double> > emiss_raw(n_Ring);  //Need to store numbers for raw emissivity to correct for rings with 0 model column densities
  vector<double> distance(n_Ring,0);  //Need to store numbers for distance to correct for rings with 0 model column densities
    
  //defining vectors vs i_Ring
  n_H_av.resize(n_Ring, 0);
  emiss_av.resize(n_Ring);
  emiss_HII.resize(n_Ring);

  //separate skymaps for HI, H2
  n_HI_av.resize(n_Ring, 0);
  n_H2_av.resize(n_Ring, 0);
  n_CO_av.resize(n_Ring, 0);
  emiss_HI.resize(n_Ring);
  emiss_H2.resize(n_Ring);
  
  //defining 2D arrays vs i_Ring and iEgamma and zeroing gas density and emissivity vectors
  for (i_Ring=0; i_Ring<n_Ring; ++i_Ring)
    { 
      emiss_raw[i_Ring].resize(galaxy.n_E_gammagrid, 0);
      emiss_av [i_Ring].resize(galaxy.n_E_gammagrid, 0);
      emiss_H2 [i_Ring].resize(galaxy.n_E_gammagrid, 0);
      emiss_HI [i_Ring].resize(galaxy.n_E_gammagrid, 0);
      emiss_HII[i_Ring].resize(galaxy.n_E_gammagrid, 0);
      
      for (int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++) 
	{ 
	  emiss_raw[i_Ring][iEgamma] = 0;
	  emiss_av [i_Ring][iEgamma] = 0;
	  emiss_H2 [i_Ring][iEgamma] = 0;
	  emiss_HI [i_Ring][iEgamma] = 0;
	  emiss_HII[i_Ring][iEgamma] = 0;
	}
    }
  if(galdef.verbose>=1) cout<<"  gen_pi0_skymap l b ="<<l<<" "<<b<<endl;

  double sinb=sin(b*dtr);
  double cosb=cos(b*dtr);
  double sinl=sin(l*dtr);
  double cosl=cos(l*dtr);
  double d=0.;
  int complete=0;

  int i_long =(int) ( (l-galaxy.long_min)/galaxy.d_long );
  int i_lat  =(int) ( (b-galaxy.lat_min )/galaxy.d_lat  );
  int co;

  if (galdef.skymap_format == 3) // HEALPIX format
    co = galaxy.hpHIR.ang2pix(SM::Coordinate(l,b).healpixAng());//Get the pixel number for accessing the skymap
   
  // integration along the line of sight

  while(complete==0) 
    {
      d += dd; 
      double zz=d*sinb;                                             // altitude of the current point 
      double RR=sqrt(Rsun*Rsun+pow(d*cosb,2)-2.0*Rsun*d*cosb*cosl); // Galactocentric distance of the current point 
      double xx,yy;

      
      if(galdef.verbose==-457) // selectable debug (-457 -old method)
	{
	  double Ro=8.5;
	  RR=sqrt(Ro*Ro+pow(d*cosb,2)-2.0*Ro*d*cosb*cosl);
	}

      i_Ring=gas_iRing(RR);// ring number
      
      // check if the column density =0 for a given ring, if so do not integrate
      // HII is added in case if there is HII gas
      if (galdef.skymap_format == 3) // HEALPIX format
	{
	  if( !(galaxy.hpHIR[co][i_Ring] + 2*fX_CO(RR)*galaxy.hpCOR[co][i_Ring] +nHII(RR,zz)) )
	    {
	       /*
	      n_H_av [i_Ring] = n_HI_av[i_Ring] = n_H2_av[i_Ring] = zero;
	      
	      for (int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
		emiss_av [i_Ring][iEgamma] = emiss_HII[i_Ring][iEgamma] = emiss_HI [i_Ring][iEgamma] = emiss_H2 [i_Ring][iEgamma] = 0.; 
	      */
	      continue;
	    }
	}

      else // (non-HEALPIX) standard (l,b) format
	{
	  if( !(galaxy.HIR.d3[i_long][i_lat][i_Ring].s[0]+2.*fX_CO(RR)*galaxy.COR.d3[i_long][i_lat][i_Ring].s[0] +nHII(RR,zz)) ) 
	    {
	       /*
	      n_H_av [i_Ring] = n_HI_av[i_Ring] = n_H2_av[i_Ring] = zero;
	      
	      for (int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
		emiss_av [i_Ring][iEgamma] = emiss_HII[i_Ring][iEgamma] = emiss_HI [i_Ring][iEgamma] = emiss_H2 [i_Ring][iEgamma] = 0.; 
	      */
	      continue;
	    }
	}

			//Calculate the analytical column densities
      nHI_av_lb =    nH_av_lb(l,b,db,ddb,d,dd,ddd,nHI);
      nCO_av_lb = 2.*nH_av_lb(l,b,db,ddb,d,dd,ddd,1.0,nH2);
      nH2_av_lb =    nCO_av_lb*fX_CO(RR);
      nHII_av_lb=    nH_av_lb(l,b,db,ddb,d,dd,ddd,nHII);
    
     // checks if we got to the Galactic boundary in 2D, if so stop integration
      if(gcr[0].n_spatial_dimensions==2) 
	{
	  // find the nearest grid points on the LEFT side of the current point
	  ir=(int)((RR-galaxy.r_min)/galaxy.dr);
	  iz=(int)((zz-galaxy.z_min)/galaxy.dz);

	  if(galdef.verbose==-457)// selectable debug (-457 -old method, the nearest grid point)
	    {
	      ir=(int)((RR-galaxy.r_min)/galaxy.dr + 0.5);
	      iz=(int)((zz-galaxy.z_min)/galaxy.dz + 0.5);
	      nHI_av_lb =    galaxy.n_HI.d2[ir][iz].s[0];
	      nH2_av_lb = 2.*galaxy.n_H2.d2[ir][iz].s[0];
	    }

	  // checks if we got to the Galactic boundary in 2D, if so stop integration
	  if(RR>galaxy.r_max)                    complete=1;
	  if(zz<galaxy.z_min || zz>galaxy.z_max) complete=1;

	  if(ir>galaxy.n_rgrid-1) {              complete=1; ir=galaxy.n_rgrid-1; }
	  
	  if(iz<0               ) {              complete=1; iz=0; } 
	  if(iz>galaxy.n_zgrid-1) {              complete=1; iz=galaxy.n_zgrid-1; }
	  
	  if(galdef.verbose==-457 && abs(zz)>1.) complete=1;  // selectable debug (-457 -old method)
	  // cout<<"d RR zz ir iz "<<d<<" "<<RR<<" "<<zz<<" "<<ir<<" "<<iz<<endl;

	} // particle.n_spatial_dimensions==2
      
      // checks if we got to the Galactic boundary in 2D, if so stop integration
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

      // gas column along the line of sight for a particular ring
      n_HI_av[i_Ring] += dd*kpc2cm* nHI_av_lb;            // HI
      n_H2_av[i_Ring] += dd*kpc2cm*            nH2_av_lb; // 2H2
      n_CO_av[i_Ring] += dd*kpc2cm*            nCO_av_lb; // 2H2
      n_H_av [i_Ring] += dd*kpc2cm*(nHI_av_lb +nH2_av_lb);// HI+2H2
      distance[i_Ring] += dd*kpc2cm;

      //handling the case if the ring edge "splits" the l-pixel in two parts
      double nH2_av_lb1,nCO_av_lb1,nHI_av_lb1,nHII_av_lb1;
      int key = 0, i_Ring1;                                                                                          //IMOS20080114** 
      //Only needed when 270 < l < 90  Gulli20080626
      if (l < 90 || l > 270) {
	 double dl;
	 if (galdef.skymap_format == 3) // HEALPIX format 
	    dl = 90./galaxy.hpHIR.Nside(); //There are 4Nside pixels in each ring
	 else
	    dl = galaxy.d_long;

	 double lmin=l-dl*0.49; // 0.49 is to get close to the pixel boundary but not to the boundary
	 double lmax=l+dl*0.49;      
	 double Rmin=sqrt(Rsun*Rsun+pow(d*cosb,2)-2.0*Rsun*d*cosb*cos(lmin*dtr)); // an edge of the pixel 
	 double Rmax=sqrt(Rsun*Rsun+pow(d*cosb,2)-2.0*Rsun*d*cosb*cos(lmax*dtr)); // an edge of the pixel 

	 i_Ring1=gas_iRing(Rmin);
	 int i_Ring2=gas_iRing(Rmax); 
	 double ledge;
	 if(l > 270.) 
	    ledge=360 - asin(galaxy.R_bins[galaxy.n_Ring+i_Ring2]/Rsun)/dtr;// l corresponding to the ring edge  //IMOS20080114**
	 else 
	    ledge=asin(galaxy.R_bins[galaxy.n_Ring+i_Ring1]/Rsun)/dtr;// l corresponding to the ring edge  //IMOS20080114**

	 //The l test is needed to make sure we only perform this when the pixel 
	 //spans a tangent to the ring boundaries
	 //cout<<lmin<<", "<<lmax<<", "<<ledge<<", "<<i_Ring1<<", "<<i_Ring2<<", "<<galaxy.R_bins[galaxy.n_Ring+min(i_Ring1,i_Ring2)]<<", "<<Rsun<<endl; 
	 if(lmin < ledge && ledge < lmax && i_Ring1 != i_Ring2 && galaxy.R_bins[galaxy.n_Ring+min(i_Ring1,i_Ring2)] < Rsun) //IMOS20080114**
	 {
	    key = 1;                                                                                           //IMOS20080114**
	    if(i_Ring1 != i_Ring)
	    {
	       nHI_av_lb1 =   nH_av_lb(lmin,b,db,ddb,d,dd,ddd,nHI);
	       nCO_av_lb1 =2.*nH_av_lb(lmin,b,db,ddb,d,dd,ddd,1.0,nH2);
	       nH2_av_lb1 =   nCO_av_lb*fX_CO(Rmin);
//	      nHII_av_lb1=   nH_av_lb(lmin,b,db,ddb,d,dd,ddd,nHII); //IMOS20080114*
	    }
	    else
	    {
	       i_Ring1=i_Ring2;
	      
	       nHI_av_lb1 =   nH_av_lb(lmax,b,db,ddb,d,dd,ddd,nHI);
	       nCO_av_lb1 =2.*nH_av_lb(lmin,b,db,ddb,d,dd,ddd,1.0,nH2);
	       nH2_av_lb1 =   nCO_av_lb*fX_CO(Rmin);
//	      nHII_av_lb1=   nH_av_lb(lmax,b,db,ddb,d,dd,ddd,nHII); //IMOS20080114*
	    }
	    n_HI_av[i_Ring1] += dd*kpc2cm* nHI_av_lb1;             // HI
	    n_H2_av[i_Ring1] += dd*kpc2cm*             nH2_av_lb1; // 2H2
	    n_CO_av[i_Ring1] += dd*kpc2cm*             nCO_av_lb1; // 2H2
	    n_H_av [i_Ring1] += dd*kpc2cm*(nHI_av_lb1 +nH2_av_lb1);// HI+2H2
	    distance[i_Ring1] += dd*kpc2cm;
	 } 
      }//End l test

      //cout<<"dd d RR theta zz i_comp pi0_aniso_factor "<<dd<<" "<<d<<" "<<RR<<" "<<theta<<" "<<zz<<" "<<i_comp<<" "<<pi0_aniso_factor<<endl;
      
      for (int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++)
	{      
	  double delta, delta_HII, x[8][3],f[8],y[7]; // x[*][0,1]= (R,z), x[*][0,1,2]= (x,y,z) -array description
 
	  // g-ray emissivity per distance dd 
	  if (gcr[0].n_spatial_dimensions==2)
	    { 
	      if (ir==galaxy.n_rgrid-1 || iz==galaxy.n_zgrid-1 || galdef.verbose==-457) //whole bracket (-457 -old method)
		{
			delta = dd*kpc2cm * galaxy.bremss_emiss.d2[ir][iz].s[iEgamma];
		}
	      else  // linear interpolation
		{
		  //  x[0]=(R0,z1), x[1]=(R1,z1);  y[0] -location of the grid points
		  //  x[2]=(R0,z0), x[3]=(R1,z0);  y[1]
		  x[0][0]=galaxy.r[ir  ]; x[0][1]=galaxy.z[iz+1];  f[0]=galaxy.bremss_emiss.d2[ir  ][iz+1].s[iEgamma];
		  x[1][0]=galaxy.r[ir+1]; x[1][1]=galaxy.z[iz+1];  f[1]=galaxy.bremss_emiss.d2[ir+1][iz+1].s[iEgamma];
		  x[2][0]=galaxy.r[ir  ]; x[2][1]=galaxy.z[iz  ];  f[2]=galaxy.bremss_emiss.d2[ir  ][iz  ].s[iEgamma];
		  x[3][0]=galaxy.r[ir+1]; x[3][1]=galaxy.z[iz  ];  f[3]=galaxy.bremss_emiss.d2[ir+1][iz  ].s[iEgamma];

		  y[0]=  (f[0]-f[1])/(x[0][0]-x[1][0])*(RR-x[0][0])+f[0]; // interpolation in R
		  y[1]=  (f[2]-f[3])/(x[2][0]-x[3][0])*(RR-x[2][0])+f[2];

		  y[2]=(y[0]-y[1])/(x[0][1]-x[2][1])*(zz-x[0][1])+y[0];   // interpolation in z

//		  cout<<" y= "; for(int j=0;j<2;j++) cout<<y[j]<<" "; cout<<" >> "<<RR+zz<<endl; 
//		  for(int j=0;j<3;j++) {cout<<" x,f= "; for(int i=0;i<4;i++) cout<<x[i][j]<<" "<<f[i]<<"   "; cout<<endl;} exit(0);

		  delta = dd*kpc2cm *y[2];

		  //ionized bremsstrahlung is handled separately

		  //  x[0]=(R0,z1), x[1]=(R1,z1);  y[0] -location of the grid points
		  //  x[2]=(R0,z0), x[3]=(R1,z0);  y[1]
		  f[0]=galaxy.bremss_ionized_emiss.d2[ir  ][iz+1].s[iEgamma];
		  f[1]=galaxy.bremss_ionized_emiss.d2[ir+1][iz+1].s[iEgamma];
		  f[2]=galaxy.bremss_ionized_emiss.d2[ir  ][iz  ].s[iEgamma];
		  f[3]=galaxy.bremss_ionized_emiss.d2[ir+1][iz  ].s[iEgamma];

		  y[0]=  (f[0]-f[1])/(x[0][0]-x[1][0])*(RR-x[0][0])+f[0]; // interpolation in R
		  y[1]=  (f[2]-f[3])/(x[2][0]-x[3][0])*(RR-x[2][0])+f[2];

		  y[2]=(y[0]-y[1])/(x[0][1]-x[2][1])*(zz-x[0][1])+y[0];   // interpolation in z

		  delta_HII= dd*kpc2cm *y[2] * nHII_av_lb;//nH_av_lb(l,b,db,ddb,d,dd,ddd,nHII);
		}
	    }
	  
	  if (gcr[0].n_spatial_dimensions==3) 
	    { 	      
	      if (ix==galaxy.n_xgrid-1 || iy==galaxy.n_ygrid-1 || iz==galaxy.n_zgrid-1) //whole bracket (-457 -old method)
		{
		  delta = dd*kpc2cm*galaxy.bremss_emiss.d3[ix][iy][iz].s[iEgamma];
		  delta_HII= dd*kpc2cm*galaxy.bremss_ionized_emiss.d3[ix][iy][iz].s[iEgamma] * nHII_av_lb;//nH_av_lb(l,b,db,ddb,d,dd,ddd,nHII);
		}
	      else  // linear interpolation
		{
		  //  x[0]=(x0,z1,y0), x[1]=(x1,z1,y0);   x[4]=(x0,z1,y1), x[5]=(x1,z1,y1);
		  //  x[2]=(x0,z0,y0), x[3]=(x1,z0,y0);   x[6]=(x0,z0,y1), x[7]=(x1,z0,y1);
		  x[0][0]=galaxy.x[ix  ]; x[0][1]=galaxy.z[iz+1]; x[0][2]=galaxy.y[iy  ];  f[0]=galaxy.bremss_emiss.d3[ix  ][iy  ][iz+1].s[iEgamma];
		  x[1][0]=galaxy.x[ix+1]; x[1][1]=galaxy.z[iz+1]; x[1][2]=galaxy.y[iy  ];  f[1]=galaxy.bremss_emiss.d3[ix+1][iy  ][iz+1].s[iEgamma];
		  x[2][0]=galaxy.x[ix  ]; x[2][1]=galaxy.z[iz  ]; x[2][2]=galaxy.y[iy  ];  f[2]=galaxy.bremss_emiss.d3[ix  ][iy  ][iz  ].s[iEgamma];
		  x[3][0]=galaxy.x[ix+1]; x[3][1]=galaxy.z[iz  ]; x[3][2]=galaxy.y[iy  ];  f[3]=galaxy.bremss_emiss.d3[ix+1][iy  ][iz  ].s[iEgamma];
		  x[4][0]=galaxy.x[ix  ]; x[4][1]=galaxy.z[iz+1]; x[4][2]=galaxy.y[iy+1];  f[4]=galaxy.bremss_emiss.d3[ix  ][iy+1][iz+1].s[iEgamma];
		  x[5][0]=galaxy.x[ix+1]; x[5][1]=galaxy.z[iz+1]; x[5][2]=galaxy.y[iy+1];  f[5]=galaxy.bremss_emiss.d3[ix+1][iy+1][iz+1].s[iEgamma];
		  x[6][0]=galaxy.x[ix  ]; x[6][1]=galaxy.z[iz  ]; x[6][2]=galaxy.y[iy+1];  f[6]=galaxy.bremss_emiss.d3[ix  ][iy+1][iz  ].s[iEgamma];
		  x[7][0]=galaxy.x[ix+1]; x[7][1]=galaxy.z[iz  ]; x[7][2]=galaxy.y[iy+1];  f[7]=galaxy.bremss_emiss.d3[ix+1][iy+1][iz  ].s[iEgamma];

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

		  //ionized bremsstrahlung is handled separately

		  //  x[0]=(x0,z1,y0), x[1]=(x1,z1,y0);   x[4]=(x0,z1,y1), x[5]=(x1,z1,y1);
		  //  x[2]=(x0,z0,y0), x[3]=(x1,z0,y0);   x[6]=(x0,z0,y1), x[7]=(x1,z0,y1);
		  f[0]=galaxy.bremss_ionized_emiss.d3[ix  ][iy  ][iz+1].s[iEgamma];
		  f[1]=galaxy.bremss_ionized_emiss.d3[ix+1][iy  ][iz+1].s[iEgamma];
		  f[2]=galaxy.bremss_ionized_emiss.d3[ix  ][iy  ][iz  ].s[iEgamma];
		  f[3]=galaxy.bremss_ionized_emiss.d3[ix+1][iy  ][iz  ].s[iEgamma];
		  f[4]=galaxy.bremss_ionized_emiss.d3[ix  ][iy+1][iz+1].s[iEgamma];
		  f[5]=galaxy.bremss_ionized_emiss.d3[ix+1][iy+1][iz+1].s[iEgamma];
		  f[6]=galaxy.bremss_ionized_emiss.d3[ix  ][iy+1][iz  ].s[iEgamma];
		  f[7]=galaxy.bremss_ionized_emiss.d3[ix+1][iy+1][iz  ].s[iEgamma];

		  y[0]=  (f[0]-f[1])/(x[0][0]-x[1][0])*(xx-x[0][0])+f[0]; // interpolation in x
		  y[1]=  (f[2]-f[3])/(x[2][0]-x[3][0])*(xx-x[2][0])+f[2];
		  y[2]=  (f[4]-f[5])/(x[4][0]-x[5][0])*(xx-x[4][0])+f[4];
		  y[3]=  (f[6]-f[7])/(x[6][0]-x[7][0])*(xx-x[6][0])+f[6];

		  y[4]=(y[0]-y[1])/(x[0][1]-x[2][1])*(zz-x[0][1])+y[0];   // interpolation in z
		  y[5]=(y[2]-y[3])/(x[4][1]-x[6][1])*(zz-x[4][1])+y[2];

		  y[6]=(y[4]-y[5])/(x[0][2]-x[4][2])*(yy-x[0][2])+y[4];   // interpolation in y

		  delta_HII= dd*kpc2cm *y[6] * nHII_av_lb;//nH_av_lb(l,b,db,ddb,d,dd,ddd,nHII);
		}
	    }

	 if(galdef.verbose==-453) // selectable debug -emiss =1, no interpolation check
	 {
	    delta = dd*kpc2cm;
	    delta_HII= dd*kpc2cm * nHII_av_lb;                                    //IMOS20080114*
	 }

	 //Store raw emissivity for possible later correction
	 emiss_raw[i_Ring][iEgamma] += delta;

	  // integral g-ray emission for a particular ring 	  
	  emiss_H2 [i_Ring][iEgamma] += delta * nH2_av_lb;
	  emiss_HI [i_Ring][iEgamma] += delta * nHI_av_lb;
	  emiss_HII[i_Ring][iEgamma] += delta_HII;
	  emiss_av [i_Ring][iEgamma] += delta * (nH2_av_lb + nHI_av_lb);
	  
	  //handling the case if the ring edge "splits" the l-pixel in two parts
	  if(key == 1)                                                                //IMOS20080114**
	    {
	      emiss_H2 [i_Ring1][iEgamma] += delta * nH2_av_lb1;
	      emiss_HI [i_Ring1][iEgamma] += delta * nHI_av_lb1;
	      emiss_av [i_Ring1][iEgamma] += delta * (nH2_av_lb1 + nHI_av_lb1);
	    }

	  if(galdef.verbose==-100) // selectable debug
	    if (galdef.gamma_rays==2) //Gulli20070810
	      {
		cout<<"l b d RR zz i_Ring Egamma n_H_av n_HI_av  n_H2_av  emiss_av complete "<<l<<" "<<b<<" "<<d<<" "<<RR<<" "
		    <<zz<<" "<<i_Ring<<" "<<galaxy.E_gamma[iEgamma]<<" "
		    << n_H_av [i_Ring]<<" "
		    << n_HI_av[i_Ring]<<" "
		    << n_H2_av[i_Ring]<<" "
		    <<" "<<emiss_av[i_Ring][iEgamma]<<" "<<complete<<endl;     
	      }
	    else
	      {
		cout<<"l b d RR zz i_Ring Egamma n_H_av emiss_av complete "<<l<<" "<<b<<" "<<d<<" "<<RR<<" "
		    <<zz<<" "<<i_Ring<<" "<<galaxy.E_gamma[iEgamma]<<" "
		    << n_H_av[i_Ring]<<" "
		    <<" "<<emiss_av[i_Ring][iEgamma]<<" "<<complete<<endl;     
	      }
	  
	} // iEgamma     
    } // complete==0

	//Loop over all the rings and check for 0 values
  for (i_Ring=0; i_Ring<n_Ring; ++i_Ring)
  {
     for (int iEgamma=0; iEgamma<galaxy.n_E_gammagrid; iEgamma++) 
     {
	if (emiss_av [i_Ring][iEgamma] == 0) emiss_av [i_Ring][iEgamma] = emiss_raw[i_Ring][iEgamma]*zero;
	if (emiss_HI [i_Ring][iEgamma] == 0) emiss_HI [i_Ring][iEgamma] = emiss_raw[i_Ring][iEgamma]*zero;
	if (emiss_H2 [i_Ring][iEgamma] == 0) emiss_H2 [i_Ring][iEgamma] = emiss_raw[i_Ring][iEgamma]*zero;
	
	if (galdef.verbose == -453) {// selectable debug -emiss =1, no interpolation check
	   double delta = 1;
       	   if ( galdef.skymap_format==0 || galdef.skymap_format==2 )
	   {
	      delta = pow(galaxy.E_gamma[iEgamma],2); 
	   }
	   // integral g-ray emission for a particular ring 	  
	   emiss_H2 [i_Ring][iEgamma] /= n_H2_av[i_Ring]/n_CO_av[i_Ring]*2*delta;
	   emiss_HI [i_Ring][iEgamma] *= 1e-20/delta;
	}
     }
     if (n_HI_av[i_Ring] == 0) n_HI_av[i_Ring] = distance[i_Ring]*zero;
     if (n_H2_av[i_Ring] == 0) n_H2_av[i_Ring] = distance[i_Ring]*zero;
     if (n_H_av [i_Ring] == 0) n_H_av [i_Ring] = distance[i_Ring]*zero;
  }
  return 9;

}
