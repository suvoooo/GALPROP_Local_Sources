
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * read_isrf.cc *                                galprop package * 5/22/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include"galprop_classes.h"
#include"galprop_internal.h"
#include"fitsio.h"

#include <Units.h>
#include <PhysicalConstants.h>
#include <GalacticRadiationField.h>
#include <RadiationField.h>

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream> //check

using namespace std;
using namespace utl;
using namespace rf;

#include <ErrorLogger.h>

static double BlackBodyNumberDensity(const double energy, const double kT) {   

  const double energyOnPi = energy/kPi,
    constant = 1./(kPlanckReduced/s/eV*kSpeedOfLight_SI/cm);

  return constant*constant*constant*
    energyOnPi*energyOnPi/(exp(energy/kT) - 1.); // eV^-1 cm^-3

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//indexing of input ISRF array as read into 1D array

int isrf_index(int ir,
	       int iz,
	       int inu,
	       int icomp,
	       int nr,
	       int nz,
	       int nnu,
	       int ncomp) {

  return icomp*nnu*nz*nr + inu*nz*nr + iz*nr + ir;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
/* 
ISRF is in Hz eV cm-3 Hz-1 (or micron eV cm^-3 micron^-1)
integral  energy density Hz-1  d(nu) = integral (nu* energy density Hz-1) d(log nu)
d(log nu) is constant in this ISRF
factor=  LOG(nu(2)/nu(1)) 
*/

using namespace rf;

//void ReadISRFFormatV0(Galprop& g); // V0 format -- CMB only
//void ReadISRFFormatV1(Galprop& g); // V1 format -- FITS-based SM2000
//void ReadISRFFormatV2(Galprop& g); // V2 format -- FITS-based e.g. Porter2008
//void ReadISRFFormatV3(Galprop& g); // V3 format -- Healpix-based

static void ReadISRFFormatV0(Galprop& g) {

  INFO("Entry");

  Galaxy& galaxy = g.galaxy;
  Galdef& galdef = g.galdef;

  ostringstream buf1;
  buf1 << "Reading ISRF from no file: CMB only";
  INFO(buf1.str());

  double targetEMin = 13.6e-5, targetEMax = 13.6; // eV

  const int targetBins = 51;

  valarray<double> targetE(0., targetBins), eps1(0., targetBins), targetFreq(0., targetBins);

  for (int i = 0; i < targetBins; ++i) {

    eps1[i] = targetEMin/kElectronMass*pow(10., i*log10(targetEMax/targetEMin)/(targetBins-1));

  }

  targetE = eps1*kElectronMass; // eV

  targetFreq = targetE*1./(kPlanck_SI/e_SI);
  
  galaxy.nu_ISRF.resize(targetFreq.size());
  galaxy.nu_ISRF = targetFreq;
  
  // Only three components: stellar + scattered, infrared, CMB
  
  const unsigned long components = 3;
  
  galaxy.n_ISRF_components = components; 
  delete[] galaxy.ISRF;
  galaxy.ISRF = new Distribution[components];
 
  ostringstream fBuf;
  fBuf << "ISRF frequency grid (Hz): ";
  
  for (int inu = 0; inu < targetBins; ++inu) 
    fBuf << galaxy.nu_ISRF[inu] << " ";

  fBuf << ends;

  INFO(fBuf.str());
    
  valarray<double> cmbFlux(0., targetBins);

  for (unsigned int i = 0; i < cmbFlux.size(); ++i)
    cmbFlux[i] = BlackBodyNumberDensity(targetE[i], 2.735*kBoltzmann_SI/e_SI);

  if (galdef.n_spatial_dimensions == 2) { // 2D
    
    galaxy.ISRF[0].init(galaxy.n_rgrid, galaxy.n_zgrid, targetBins);
    galaxy.ISRF[1].init(galaxy.n_rgrid, galaxy.n_zgrid, targetBins);
    galaxy.ISRF[2].init(galaxy.n_rgrid, galaxy.n_zgrid, targetBins);
    
    ostringstream iBuf;
    iBuf << "galaxy.ISRF initialised with frequency axis dimension = " << targetBins;
    INFO(iBuf.str());
    
    for (int i = 0; i < galaxy.n_rgrid; ++i) {
      
      for (int j = 0; j < galaxy.n_zgrid; ++j) {
	
	for (int k = 0; k < targetBins; ++k) {

	  const double factor = targetE[k]*targetE[k];

	  galaxy.ISRF[0].d2[i][j].s[k] = 0;
	  
	  galaxy.ISRF[1].d2[i][j].s[k] = 0;
	  
	  galaxy.ISRF[2].d2[i][j].s[k] = cmbFlux[k]*factor;
	
	  //cout << i << " " << j << " " << r << " " << z << " " << k << " " 
	  //   << setprecision(5) << targetFreq[k] << " " 
	    //<< total*targetE[k]*targetE[k] << " " 
	    //   << stellar*targetE[k]*targetE[k] << " " 
	  // << scattered*targetE[k]*targetE[k] << " " 
	  //   << transient*targetE[k]*targetE[k] << " " 
	  //   << thermal*targetE[k]*targetE[k] << " " 
	  //   << (stellar + scattered + transient + thermal)*targetE[k]*targetE[k] << " " 
	  //   << cmbFlux[k]*targetE[k]*targetE[k] << endl;

	} // target energy bins

      } // z
      
    } // r
    
  } // 2D
    
  if (galdef.n_spatial_dimensions == 3) { // 3D
    
    galaxy.ISRF[0].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, targetBins);
    galaxy.ISRF[1].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, targetBins);
    galaxy.ISRF[2].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, targetBins);
    
    ostringstream iBuf;
    iBuf << "galaxy.ISRF initialised with frequency axis dimension = " << targetBins;
    INFO(iBuf.str());
     
    for (int i = 0; i < galaxy.n_xgrid; ++i) {
      
      for (int j = 0; j < galaxy.n_ygrid; ++j) {
		
	for (int k = 0; k < galaxy.n_zgrid; ++k) {
	  
	  // Recall when retrieving the data from the skymaps it is returned
	  // in units eV cm^-2 s^-1 sr^-1
	  
	  for (int l = 0; l < targetBins; ++l) {

	    const double factor = targetE[l]*targetE[l];
	    
	    galaxy.ISRF[0].d3[i][j][k].s[l] = 0;
	    
	    galaxy.ISRF[1].d3[i][j][k].s[l] = 0;
	    
	    galaxy.ISRF[2].d3[i][j][k].s[l] = cmbFlux[l]*factor;
	    
	  } // target energy bins

	} // z
	
      } // y
      
    } // x
    
  } // 3D
 
  INFO("Exit");

  //exit(0);

}

static void ReadISRFFormatV1(Galprop& g) {

  INFO("Entry");

  Galaxy& galaxy = g.galaxy;
  Galdef& galdef = g.galdef;

  int status = 0;

  fitsfile* fptr;
  //char ISRF_filename[200];
  int NAXIS,NAXIS1,NAXIS2,NAXIS3,NAXIS4;
  float CRVAL1,CRVAL2,CRVAL3;
  float CDELT1,CDELT2,CDELT3;
  char comment[100];
  
  //strcpy(ISRF_filename,configure.fits_directory);
  //strcat(ISRF_filename,"isrf_interp_04_000015");
  //strcat(ISRF_filename,"porter_ISRF.fits");           //AWS20050225
  //strcat(ISRF_filename,"porter_RFScattering10kpc.fits");//AWS20050301
  //strcat(ISRF_filename,galdef.ISRF_file);//AWS20050301
  
  const std::string fitsDirectory = g.configure.fFITSDataDirectory;
  const std::string isrfFilename = g.galdef.ISRF_file;
  const std::string filename = fitsDirectory + isrfFilename;
  
  //if (g.galdef.verbose>=1)cout<<" reading ISRF from "<< filename <<endl;

  ostringstream buf;
  buf << "Reading ISRF from " << filename;
  INFO(buf.str());

  if( fits_open_file(&fptr,filename.c_str(),READONLY,&status) ) {
    buf.str("");
    buf<<"read isrf open status= "<<status;
    WARNING(buf.str());
  }
  if (0 == fptr){
    buf.str("");
    buf<<"Cannot open file "<<filename;
    FATAL(buf.str());
    exit(-2);
  }
  
  if( fits_read_key(fptr,TINT,"NAXIS" ,&NAXIS ,comment,&status) ) {
    buf.str("");
    buf<<"read NAXIS open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TINT,"NAXIS1",&NAXIS1,comment,&status) ) {
    buf.str("");
    buf<<"read NAXIS1 open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TINT,"NAXIS2",&NAXIS2,comment,&status) ) {
    buf.str("");
    buf<<"read NAXIS2 open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TINT,"NAXIS3",&NAXIS3,comment,&status) ) {
    buf.str("");
    buf<<"read NAXIS3 open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TINT,"NAXIS4",&NAXIS4,comment,&status) ) {
    buf.str("");
    buf<<"read NAXIS4 open status= "<<status;
    WARNING(buf.str());
  }
  
  if( fits_read_key(fptr,TFLOAT,"CRVAL1",&CRVAL1,comment,&status) ) {
    buf.str("");
    buf<<"read CRVAL1 open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TFLOAT,"CRVAL2",&CRVAL2,comment,&status) ) {
    buf.str("");
    buf<<"read CRVAL2 open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TFLOAT,"CRVAL3",&CRVAL3,comment,&status) ) {
    buf.str("");
    buf<<"read CRVAL3 open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TFLOAT,"CDELT1",&CDELT1,comment,&status) ) {
    buf.str("");
    buf<<"read CDELT1 open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TFLOAT,"CDELT2",&CDELT2,comment,&status) ) {
    buf.str("");
    buf<<"read CDELT2 open status= "<<status;
    WARNING(buf.str());
  }
  if( fits_read_key(fptr,TFLOAT,"CDELT3",&CDELT3,comment,&status) ) {
    buf.str("");
    buf<<"read CDELT3 open status= "<<status;
    WARNING(buf.str());
  }
  
  buf.str("");
  buf<<" NAXIS = "<<NAXIS <<endl;
  DEBUGLOG(buf.str());
  buf.str("");
  buf<<" NAXIS1,2,3,4 = "<<NAXIS1<<" "<<NAXIS2<<" "<<NAXIS3<<" "<<NAXIS4<<endl;
  DEBUGLOG(buf.str());
  buf.str("");
  buf<<" CRVAL1,2,3 = "<<CRVAL1<<" "<<CRVAL2<<" "<<CRVAL3<<endl;
  DEBUGLOG(buf.str());
  buf.str("");
  buf<<" CDELT1,2,3 = "<<CDELT1<<" "<<CDELT2<<" "<<CDELT3<<endl;
  DEBUGLOG(buf.str());
  
  long nelements=NAXIS1*NAXIS2*NAXIS3*NAXIS4, felement=1;
  float *isrf_in=new float[nelements];
  float nulval=0;
  int anynul;
    
  if (fits_read_img(fptr, TFLOAT, felement, nelements, &nulval, isrf_in, 
		    &anynul, &status)){
    buf.str("");
    buf<<"read isrf open status= "<<status;
    WARNING(buf.str());
  }
  
  // for(int i=0; i<nelements; i++) cout<<isrf_in[i]<<" ";
  fits_close_file(fptr,&status);
  
  INFO("generating galaxy.ISRF:");
  
  galaxy.n_ISRF_components = NAXIS4;
  delete[] galaxy.ISRF;
  galaxy.ISRF = new Distribution[NAXIS4];
  
  for(int i=0; i<NAXIS4; i++) {
    
    if(galdef.n_spatial_dimensions==2) {             // ==== 2D ====
      
      galaxy.ISRF[i].init(galaxy.n_rgrid,galaxy.n_zgrid,NAXIS3);
      buf.str("");
      buf<<" galaxy.ISRF initialized with frequency axis dimension="
	  <<galaxy.ISRF[i].n_pgrid;
      INFO(buf.str());
      
      for(int ir=0; ir<galaxy.n_rgrid; ir++) { 
	
	for(int iz=0; iz<galaxy.n_zgrid; iz++) {
	  
	  int irr=(int)((     galaxy.r[ir] -CRVAL1) /CDELT1+0.5);//IMOS20060420
	  int izz=(int)((fabs(galaxy.z[iz])-CRVAL2) /CDELT2+0.5);//IMOS20060420
	  if(irr>NAXIS1-2) irr=NAXIS1-2;
	  if(izz>NAXIS2-2) izz=NAXIS2-2;
	  float rr=CRVAL1+irr*CDELT1;
	  float zz=CRVAL2+izz*CDELT2;
	  
	  // cout<<"r z irr izz rr zz "<< galaxy.r[ir]<<" "<<galaxy.z[iz]<<"     "<<irr<<" "<<izz<<"     "<<rr<<" "<<zz<<endl;
	  
	  for(int inu=0; inu<galaxy.ISRF[i].n_pgrid; inu++) {
	    float v1=isrf_in[isrf_index(irr  ,izz  ,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	    float v2=isrf_in[isrf_index(irr+1,izz  ,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	    float v3=isrf_in[isrf_index(irr  ,izz+1,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	    float v4=isrf_in[isrf_index(irr+1,izz+1,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	    float v5=v1+(v2-v1)*(galaxy.r[ir]-rr)/CDELT1;
	    float v6=v3+(v4-v3)*(galaxy.r[ir]-rr)/CDELT1;
	    float value=v5+(v6-v5)*(fabs(galaxy.z[iz])-zz)/CDELT2;
	    if(value<0.0) value=0.0;
	    // reverse scale from wavelength to frequency
	    galaxy.ISRF[i].d2[ir][iz].s[galaxy.ISRF[i].n_pgrid-1-inu] = value;
	    
	     // cout<<"r z irr izz rr zz "<< galaxy.r[ir]<<" "<<galaxy.z[iz]<<"   "<<irr<<" "<<izz<<"   "<<rr<<" "<<zz<<" "<< v5+(v6-v5)*(galaxy.z[iz]-zz)/CDELT2<<endl;
	  }  //  inu
	}  //  iz
      }  //  ir
    }  //  2D
    
    if(galdef.n_spatial_dimensions==3) {              // ==== 3D ====
      
      galaxy.ISRF[i].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid,NAXIS3);
      for(int ix=0; ix<galaxy.n_xgrid; ix++) {  
	
	for(int iy=0; iy<galaxy.n_ygrid; iy++) {
	  
	  for(int iz=0; iz<galaxy.n_zgrid; iz++) {
	    
	    float r=sqrt(pow(galaxy.x[ix],2)+pow(galaxy.y[iy],2));
	    int irr=(int)((            r     -CRVAL1) /CDELT1+0.5);//IMOS20060420
	    int izz=(int)((fabs(galaxy.z[iz])-CRVAL2) /CDELT2+0.5);//IMOS20060420
	    if(irr>NAXIS1-2) irr=NAXIS1-2;
	    if(izz>NAXIS2-2) izz=NAXIS2-2;
	    
	    float rr=CRVAL1+irr*CDELT1;
	    float zz=CRVAL2+izz*CDELT2;
	    
	    // cout<<"r z irr izz rr zz "<< galaxy.r[ir]<<" "<<galaxy.z[iz]<<"     "<<irr<<" "<<izz<<"     "<<rr<<" "<<zz<<endl;
	    
	    for(int inu=0; inu<galaxy.ISRF[i].n_pgrid; inu++) {
	      
	      float v1=isrf_in[isrf_index(irr  ,izz  ,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	      float v2=isrf_in[isrf_index(irr+1,izz  ,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	      float v3=isrf_in[isrf_index(irr  ,izz+1,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	      float v4=isrf_in[isrf_index(irr+1,izz+1,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	      float v5=v1+(v2-v1)*(       r    -rr)/CDELT1;
	      float v6=v3+(v4-v3)*(       r    -rr)/CDELT1;
	      float value=v5+(v6-v5)*(fabs(galaxy.z[iz])-zz)/CDELT2;
	      if(value<0.0) value=0.0;
	      // reverse scale from wavelength to frequency
	      galaxy.ISRF[i].d3[ix][iy][iz].s[galaxy.ISRF[i].n_pgrid-1-inu]= value;
	    }  //  inu   
	  }  //  ix
	}  //  iy
      }  //  iz
    }  //  3D
    
    if(galdef.verbose>=10) {
      
      buf.str("");
      buf<<"ISRF component "<<i+1;
      INFO(buf.str());
      galaxy.ISRF[i].print();
      
    }
    
  }  // isrf component i
  
  // Create the array of ISRF frequencies
  // using wavelength in microns for axis 3 of input ISRF on log10 scale.
  // Reverse scale so that frequency increases.
  
  galaxy.nu_ISRF.resize(galaxy.ISRF[0].n_pgrid);
  
  // microns -> cm; nu=c/lambda
  for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++)
    galaxy.nu_ISRF[galaxy.ISRF[0].n_pgrid-1-inu]= 
      c/(pow(10.,1.*CRVAL3+inu*CDELT3)*1.0e-4);
  
  DEBUGLOG(" ISRF frequency grid (Hz):");
  for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++) {
     buf.str("");
    buf<<"inu galaxy.nu_ISRF[inu] "<<inu<<" "<< galaxy.nu_ISRF[inu];
    DEBUGLOG(buf.str());
  }
  
  delete[] isrf_in;//AWS20010216
  
  INFO("Exit");

}

static void ReadISRFFormatV2(Galprop& g) {

  INFO("Entry");

  Galaxy& galaxy = g.galaxy;
  Galdef& galdef = g.galdef;
  gp::Configure& configure = g.configure;

  // Read in using new format
  
  const std::string fitsDirectory = configure.fFITSDataDirectory;
  const std::string isrfFilename = galdef.ISRF_file;
  const std::string filename = fitsDirectory + isrfFilename;

  ostringstream buf1;
  buf1 << "Reading ISRF from " << filename;
  INFO(buf1.str());

  //galaxy.fISRFModel = new rf::GalacticRadiationField(filename, true); // For now, no angular information is read in -- TAP 09082007

  rf::GalacticRadiationField rf(filename, true); //(galaxy.fISRFModel);

  INFO("Generating galaxy.ISRF");
  ostringstream buf;

  // Only three components: stellar + scattered, infrared, CMB
  
  const unsigned long components = 3;
  
  galaxy.n_ISRF_components = components; 
  delete[] galaxy.ISRF;
  galaxy.ISRF = new Distribution[components];
  
  const std::valarray<double>& wl = rf.GetWavelengthData();
  
  // Cheat a bit. Get what the wavelength log(delta) is and extend the 
  // range up to 10000 microns (if it already doesn't go that far). This 
  // ensures we get the full CMB as well.
  
  const unsigned long rawWlBins = wl.size();

  std::valarray<double> wavelengthData;
  
  if (wl[rawWlBins-1] < 1e4) {
    
    const double wlDelta = log10(wl[1]/wl[0]);
    
    const unsigned long bins = long(log10(1e4/wl[0])/wlDelta) + 1;
    
    wavelengthData.resize(bins);
    
    for (unsigned long i = 0; i < bins; ++i)
      wavelengthData[i] = 
	(i < rawWlBins ? wl[i] : wl[0]*pow(10.0, i*wlDelta));
    
  } else {
    
    wavelengthData.resize(rawWlBins);
    wavelengthData = wl;
    
  }
  
  const int wlBins = wavelengthData.size();
  
  // Create the array of ISRF frequencies using wavelength in microns. 
  // Reverse scale so that frequency increases.
  
  galaxy.nu_ISRF.resize(wlBins);
  
  for (int i = 0; i < wlBins; ++i)
    galaxy.nu_ISRF[wlBins - 1 - i] = 
      utl::kSpeedOfLight_SI/(wavelengthData[i]*utl::micron/utl::m);
  
  DEBUGLOG(" ISRF frequency grid (Hz):");
  for(int inu=0; inu<wlBins; inu++) {
     buf.str("");
    buf<<"inu galaxy.nu_ISRF[inu] "<<inu<<" "<< galaxy.nu_ISRF[inu];
    DEBUGLOG(buf.str());
  }
  
  // Some of this is truly, horribly, awful. However, I can't fix it until
  // Galprop undergoes a full re-write. I can't see that happening in the
  // near future -- TAP20072301
  
  if (galdef.n_spatial_dimensions == 2) { // 2D
    
    galaxy.ISRF[0].init(galaxy.n_rgrid, galaxy.n_zgrid, wlBins);
    galaxy.ISRF[1].init(galaxy.n_rgrid, galaxy.n_zgrid, wlBins);
    galaxy.ISRF[2].init(galaxy.n_rgrid, galaxy.n_zgrid, wlBins);
    
    buf.str("");
    buf << " galaxy.ISRF initialized with frequency axis dimension = "
	 << wlBins << endl;
    INFO(buf.str());
    
    for (int i = 0; i < galaxy.n_rgrid; ++i) {
      
      for (int j = 0; j < galaxy.n_zgrid; ++j) {
	
	const double r = galaxy.r[i], z = fabs(galaxy.z[j]);
	
	for (int k = 0; k < wlBins; ++k) {
	  
	  const double wl = wavelengthData[k];
	  
	  // Scale is reversed from wavelength to frequency
	  
	  const double stellar = 
	    rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::STELLAR);
	  
	  const double scattered = 
	    rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::SCATTERED);
	  
	  const double infrared = 
	    rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::INFRARED);
	  
	  const double cmb = 
	    rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::CMB);
	  
	  galaxy.ISRF[0].d2[i][j].s[galaxy.ISRF[0].n_pgrid - 1 - k] = stellar + scattered;
	  
	  galaxy.ISRF[1].d2[i][j].s[galaxy.ISRF[1].n_pgrid - 1 - k] = infrared;
	  
	  galaxy.ISRF[2].d2[i][j].s[galaxy.ISRF[2].n_pgrid - 1 - k] = cmb;

	  //cout << i << " " << j << " " << r << " " << z << " " << galaxy.ISRF[0].n_pgrid - 1 - k << " " 
	  //   << setprecision(5) << galaxy.nu_ISRF[galaxy.ISRF[0].n_pgrid - 1 - k] << " "
	  //   << stellar << " "
	  //   << scattered << " "
	  //   << infrared << " "
	  //   << (stellar + scattered + infrared) << " " 
	  //   << cmb << endl;

	} // wl
	
      } // z
      
    } // r
    
  } // 2D
    
  if (galdef.n_spatial_dimensions == 3) { // 3D
    
    galaxy.ISRF[0].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, wlBins);
    galaxy.ISRF[1].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, wlBins);
    galaxy.ISRF[2].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, wlBins);
    
    ostringstream buf;
    buf<< " galaxy.ISRF initialized with frequency axis dimension = " << wlBins;
    INFO(buf.str());
    
    for (int i = 0; i < galaxy.n_xgrid; ++i) {
      
      for (int j = 0; j < galaxy.n_ygrid; ++j) {
	
	const double r = sqrt(galaxy.x[i]*galaxy.x[i] + galaxy.y[j]*galaxy.y[j]);
	
	for (int k = 0; k < galaxy.n_zgrid; ++k) {
	  
	  const double z = fabs(galaxy.z[k]);
	  
	  for (int l = 0; l < wlBins; ++l) {
	    
	    const double wl = wavelengthData[l];
	    
	    // Scale is reversed from wavelength to frequency
	    
	    const double stellar = 
	      rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::STELLAR);
	    
	    const double scattered = 
	      rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::SCATTERED);
	    
	    const double infrared = 
	      rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::INFRARED);
	    
	    const double cmb = 
	      rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::CMB);
	    
	    galaxy.ISRF[0].d3[i][j][k].s[galaxy.ISRF[0].n_pgrid - 1 - l] = stellar + scattered;
	    
	    galaxy.ISRF[1].d3[i][j][k].s[galaxy.ISRF[1].n_pgrid - 1 - l] = infrared;
	    
	    galaxy.ISRF[2].d3[i][j][k].s[galaxy.ISRF[2].n_pgrid - 1 - l] = cmb;
	    
	  } // wl
	  
	} // z
	
      } // y
      
    } // x
    
  } // 3D
 
  INFO("Exit");

}

static void ReadISRFFormatV3(Galprop& g) {

  INFO("Entry");

  Galaxy& galaxy = g.galaxy;
  Galdef& galdef = g.galdef;
  gp::Configure& configure = g.configure;

  const std::string fitsDirectory = configure.fFITSDataDirectory;
  const std::string isrfFilename = galdef.ISRF_file;
  const std::string filename = fitsDirectory + isrfFilename;

  ostringstream buf1;
  buf1 << "Reading ISRF from " << filename;
  INFO(buf1.str());

  if (g.galdef.verbose >= 1) {

    ostringstream addInfoBuf;
    addInfoBuf << "FITS directory: " << fitsDirectory << endl 
	       << "ISRF filename: " << isrfFilename << endl 
	       << "ISRF healpix order: " << galdef.ISRF_healpixOrder;
    INFO(addInfoBuf.str());

  }

  const double targetEMin = 13.6e-5, targetEMax = 13.6; // eV

  const int targetBins = 101;

  valarray<double> targetE(0., targetBins), eps1(0., targetBins), targetFreq(0., targetBins);

  for (int i = 0; i < targetBins; ++i) {

    eps1[i] = targetEMin/kElectronMass*pow(10., i*log10(targetEMax/targetEMin)/(targetBins-1));

  }

  targetE = eps1*kElectronMass; // eV

  targetFreq = targetE*1./(kPlanck_SI/e_SI);

  RadiationField rf(filename, targetFreq, galdef.ISRF_healpixOrder);
  
  galaxy.nu_ISRF.resize(targetFreq.size());
  galaxy.nu_ISRF = targetFreq;
  
  // Only three components: stellar + scattered, infrared, CMB
  
  const unsigned int components = 3;
  
  galaxy.n_ISRF_components = components; 
  delete[] galaxy.ISRF;
  galaxy.ISRF = new Distribution[components];
 
  ostringstream fBuf;
  fBuf << "ISRF frequency grid (Hz): ";
  
  for (int inu = 0; inu < targetBins; ++inu) 
    fBuf << galaxy.nu_ISRF[inu] << " ";

  fBuf << ends;

  INFO(fBuf.str());
    
  valarray<double> cmbFlux(0., targetBins);

  for (unsigned long i = 0; i < cmbFlux.size(); ++i)
    cmbFlux[i] = BlackBodyNumberDensity(targetE[i], 2.735*kBoltzmann_SI/e_SI);

  if (2 == galdef.n_spatial_dimensions) { // 2D
    
    galaxy.ISRF[0].init(galaxy.n_rgrid, galaxy.n_zgrid, targetBins);
    galaxy.ISRF[1].init(galaxy.n_rgrid, galaxy.n_zgrid, targetBins);
    galaxy.ISRF[2].init(galaxy.n_rgrid, galaxy.n_zgrid, targetBins);
    
    ostringstream iBuf;
    iBuf << "galaxy.ISRF initialised with frequency axis dimension = " << targetBins;
    INFO(iBuf.str());
    
    for (int i = 0; i < galaxy.n_rgrid; ++i) {
      
      for (int j = 0; j < galaxy.n_zgrid; ++j) {
	
	const double r = galaxy.r[i], z = galaxy.z[j];

	rf::RadiationField::ThreeVector pos(r, 0, z);

	// Recall when retrieving the data from the skymaps it is returned
	// in units eV^-1 cm^-3 sr^-1 -- simply sum over all pixels and 
	// multiply by solid angle to get number density

	//const Skymap<double> skymapTotal = rf.GetSkymap(pos, RadiationField::TOTAL, galdef.ISRF_healpixOrder);
	
	//const Skymap<double> skymapDirect = rf.GetSkymap(pos, RadiationField::DIRECT, galdef.ISRF_healpixOrder);

	//const Skymap<double> skymapScattered = rf.GetSkymap(pos, RadiationField::SCATTERED, galdef.ISRF_healpixOrder);

	//const Skymap<double> skymapTransient = rf.GetSkymap(pos, RadiationField::TRANSIENT, galdef.ISRF_healpixOrder);

	//const Skymap<double> skymapThermal = rf.GetSkymap(pos, RadiationField::THERMAL, galdef.ISRF_healpixOrder);

	valarray<double> numberDensityDirect(0., targetBins);
	valarray<double> numberDensityScattered(0., targetBins);
	valarray<double> numberDensityTransient(0., targetBins);
	valarray<double> numberDensityThermal(0., targetBins);

	numberDensityDirect = rf.GetNumberDensity(pos, RadiationField::DIRECT);
	numberDensityScattered = rf.GetNumberDensity(pos, RadiationField::SCATTERED);
	numberDensityTransient = rf.GetNumberDensity(pos, RadiationField::TRANSIENT);
	numberDensityThermal = rf.GetNumberDensity(pos, RadiationField::THERMAL);
	  
	for (int k = 0; k < targetBins; ++k) {

	  const double factor = targetE[k]*targetE[k];

	  //const double total = skymapTotal.sum(k)*skymapTotal.solidAngle();

	  const double stellar = numberDensityDirect[k];//skymapDirect.sum(k)*skymapDirect.solidAngle();
	  const double scattered = numberDensityScattered[k];//skymapScattered.sum(k)*skymapScattered.solidAngle();
	  const double transient = numberDensityTransient[k];//skymapTransient.sum(k)*skymapTransient.solidAngle();
	  const double thermal = numberDensityThermal[k];//skymapThermal.sum(k)*skymapThermal.solidAngle();

	  galaxy.ISRF[0].d2[i][j].s[k] = (stellar + scattered)*factor;
	  
	  galaxy.ISRF[1].d2[i][j].s[k] = (transient + thermal)*factor;
	  
	  galaxy.ISRF[2].d2[i][j].s[k] = cmbFlux[k]*factor;

	  //cout << i << " " << j << " " << r << " " << z << " " << k << " " << galaxy.ISRF[0].d2[i][j].s[k] << " " << galaxy.ISRF[1].d2[i][j].s[k] << " " << galaxy.ISRF[2].d2[i][j].s[k] << endl;


	} // target energy bins

	//exit(0);

      } // z
      
    } // r
    
  } // 2D

  //exit(0);
    
  if (3 == galdef.n_spatial_dimensions) { // 3D

    //ofstream fout("/home/yuko/GALPROP/check.tab"); //check
    
    galaxy.ISRF[0].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, targetBins);
    galaxy.ISRF[1].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, targetBins);
    galaxy.ISRF[2].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, targetBins);
    
    ostringstream iBuf;
    iBuf << "galaxy.ISRF initialised with frequency axis dimension = " << targetBins;
    INFO(iBuf.str());
     
    for (int i = 0; i < galaxy.n_xgrid; ++i) {
      
      for (int j = 0; j < galaxy.n_ygrid; ++j) {
	
	const double r = sqrt(galaxy.x[i]*galaxy.x[i] + galaxy.y[j]*galaxy.y[j]);
	
	for (int k = 0; k < galaxy.n_zgrid; ++k) {
	  
	  const double z = galaxy.z[k];
	
	  rf::RadiationField::ThreeVector pos(r, 0, z);

	  // Recall when retrieving the data from the skymaps it is returned
	  // in units eV cm^-2 s^-1 sr^-1

	  //const Skymap<double> skymapTotal = rf.GetSkymap(pos, RadiationField::TOTAL, galdef.ISRF_healpixOrder);
	  
	  //const Skymap<double> skymapDirect = rf.GetSkymap(pos, RadiationField::DIRECT, galdef.ISRF_healpixOrder);

	  //const Skymap<double> skymapScattered = rf.GetSkymap(pos, RadiationField::SCATTERED, galdef.ISRF_healpixOrder);
	  
	  //const Skymap<double> skymapTransient = rf.GetSkymap(pos, RadiationField::TRANSIENT, galdef.ISRF_healpixOrder);
	
	  //const Skymap<double> skymapThermal = rf.GetSkymap(pos, RadiationField::THERMAL, galdef.ISRF_healpixOrder);
	  
	  valarray<double> numberDensityDirect(0., targetBins);
	  valarray<double> numberDensityScattered(0., targetBins);
	  valarray<double> numberDensityTransient(0., targetBins);
	  valarray<double> numberDensityThermal(0., targetBins);
	  
	  numberDensityDirect = rf.GetNumberDensity(pos, RadiationField::DIRECT);
	  numberDensityScattered = rf.GetNumberDensity(pos, RadiationField::SCATTERED);
	  numberDensityTransient = rf.GetNumberDensity(pos, RadiationField::TRANSIENT);
	  numberDensityThermal = rf.GetNumberDensity(pos, RadiationField::THERMAL);

	  for (int l = 0; l < targetBins; ++l) {

	    const double factor = targetE[l]*targetE[l];
	    
	    //const double total = skymapTotal.sum(l)*skymapTotal.solidAngle();
	    
	    const double stellar = numberDensityDirect[l];//skymapDirect.sum(l)*skymapDirect.solidAngle();
	    const double scattered = numberDensityScattered[l];//skymapScattered.sum(l)*skymapScattered.solidAngle();
	    const double transient = numberDensityTransient[l];//skymapTransient.sum(l)*skymapTransient.solidAngle();
	    const double thermal = numberDensityThermal[l];//skymapThermal.sum(l)*skymapThermal.solidAngle();
	    	    
	    galaxy.ISRF[0].d3[i][j][k].s[l] = (stellar + scattered)*factor;

	    galaxy.ISRF[1].d3[i][j][k].s[l] = (transient + thermal)*factor;

	    galaxy.ISRF[2].d3[i][j][k].s[l] = cmbFlux[l]*factor;
	    //fout << galaxy.ISRF[2].d3[i][j][k].s[l] << endl; //check	    	    
	  } // target energy bins

	} // z
	
      } // y
      
    } // x
    
  } // 3D
 
  INFO("Exit");

  //exit(0);

}

int Galprop::read_isrf(const int version) {

  INFO("Entry");

  assert (version >= 0 && version <= 3); // Update later for other versions -- only up to version 2 is supported for current release

  int status = 0;

  if (0 == version)
    ReadISRFFormatV0(*this);

  if (1 == version)
    ReadISRFFormatV1(*this);
 
  if (2 == version)
    ReadISRFFormatV2(*this);

  if (3 == version)
    ReadISRFFormatV3(*this);

  //exit(0);

  /*if (!useNewMethod) {

    fitsfile* fptr;
    //char ISRF_filename[200];
    int NAXIS,NAXIS1,NAXIS2,NAXIS3,NAXIS4;
    float CRVAL1,CRVAL2,CRVAL3;
    float CDELT1,CDELT2,CDELT3;
    char comment[100];
    
    //strcpy(ISRF_filename,configure.fits_directory);
    //strcat(ISRF_filename,"isrf_interp_04_000015");
    //strcat(ISRF_filename,"porter_ISRF.fits");           //AWS20050225
    //strcat(ISRF_filename,"porter_RFScattering10kpc.fits");//AWS20050301
    //strcat(ISRF_filename,galdef.ISRF_file);//AWS20050301
    
    const std::string fitsDirectory = configure.fits_directory;
    const std::string isrfFilename = galdef.ISRF_file;
    const std::string filename = fitsDirectory + isrfFilename;

    cout<<" reading ISRF from "<< filename <<endl;

    if( fits_open_file(&fptr,filename.c_str(),READONLY,&status) ) 
      cout<<"read isrf open status= "<<status<<endl;
    if(fptr==NULL)exit(-2);
    
    if( fits_read_key(fptr,TINT,"NAXIS" ,&NAXIS ,comment,&status) ) 
      cout<<"0read isrf status= "<<status<<endl;
    if( fits_read_key(fptr,TINT,"NAXIS1",&NAXIS1,comment,&status) ) 
      cout<<"1read isrf status= "<<status<<endl;
    if( fits_read_key(fptr,TINT,"NAXIS2",&NAXIS2,comment,&status) ) 
      cout<<"2read isrf status= "<<status<<endl;
    if( fits_read_key(fptr,TINT,"NAXIS3",&NAXIS3,comment,&status) ) 
      cout<<"3read isrf status= "<<status<<endl;
    if( fits_read_key(fptr,TINT,"NAXIS4",&NAXIS4,comment,&status) ) 
      cout<<"4read isrf status= "<<status<<endl;

    if( fits_read_key(fptr,TFLOAT,"CRVAL1",&CRVAL1,comment,&status) ) 
      cout<<"5read isrf status= "<<status<<endl;
    if( fits_read_key(fptr,TFLOAT,"CRVAL2",&CRVAL2,comment,&status) ) 
      cout<<"6read isrf status= "<<status<<endl;
    if( fits_read_key(fptr,TFLOAT,"CRVAL3",&CRVAL3,comment,&status) ) 
      cout<<"7read isrf status= "<<status<<endl;
    if( fits_read_key(fptr,TFLOAT,"CDELT1",&CDELT1,comment,&status) ) 
      cout<<"8read isrf status= "<<status<<endl;
    if( fits_read_key(fptr,TFLOAT,"CDELT2",&CDELT2,comment,&status) ) 
      cout<<"9read isrf status= "<<status<<endl;
    if( fits_read_key(fptr,TFLOAT,"CDELT3",&CDELT3,comment,&status) ) 
      cout<<"/read isrf status= "<<status<<endl;
    
    cout<<" NAXIS = "<<NAXIS <<endl;
    cout<<" NAXIS1,2,3,4 = "<<NAXIS1<<" "<<NAXIS2<<" "<<NAXIS3<<" "<<NAXIS4<<endl;
    cout<<" CRVAL1,2,3 = "<<CRVAL1<<" "<<CRVAL2<<" "<<CRVAL3<<endl;
    cout<<" CDELT1,2,3 = "<<CDELT1<<" "<<CDELT2<<" "<<CDELT3<<endl;

    long nelements=NAXIS1*NAXIS2*NAXIS3*NAXIS4, felement=1;
    float *isrf_in=new float[nelements];
    float nulval=0;
    int anynul;
    
    if (fits_read_img(fptr, TFLOAT, felement, nelements, &nulval, isrf_in, 
		      &anynul, &status))
      cout<<"#read isrf status= "<<status<<endl;

// for(int i=0; i<nelements; i++) cout<<isrf_in[i]<<" ";
    
   cout<<"generating galaxy.ISRF:"<<endl;

   galaxy.n_ISRF_components = NAXIS4;
   galaxy.ISRF = new Distribution[NAXIS4];

   for(int i=0; i<NAXIS4; i++) {

     if(galdef.n_spatial_dimensions==2) {             // ==== 2D ====
      
       galaxy.ISRF[i].init(galaxy.n_rgrid,galaxy.n_zgrid,NAXIS3);
       cout<<" galaxy.ISRF initialized with frequency axis dimension="
	   <<galaxy.ISRF[i].n_pgrid<<endl;
       
       for(int ir=0; ir<galaxy.n_rgrid; ir++) { 
	 
	 for(int iz=0; iz<galaxy.n_zgrid; iz++) {

	   int irr=(int)((     galaxy.r[ir] -CRVAL1) /CDELT1+0.5);//IMOS20060420
	   int izz=(int)((fabs(galaxy.z[iz])-CRVAL2) /CDELT2+0.5);//IMOS20060420
	   if(irr>NAXIS1-2) irr=NAXIS1-2;
	   if(izz>NAXIS2-2) izz=NAXIS2-2;
	   float rr=CRVAL1+irr*CDELT1;
	   float zz=CRVAL2+izz*CDELT2;
	   
	   // cout<<"r z irr izz rr zz "<< galaxy.r[ir]<<" "<<galaxy.z[iz]<<"     "<<irr<<" "<<izz<<"     "<<rr<<" "<<zz<<endl;
	   
	   for(int inu=0; inu<galaxy.ISRF[i].n_pgrid; inu++) {
	     float v1=isrf_in[isrf_index(irr  ,izz  ,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	     float v2=isrf_in[isrf_index(irr+1,izz  ,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	     float v3=isrf_in[isrf_index(irr  ,izz+1,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	     float v4=isrf_in[isrf_index(irr+1,izz+1,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
	     float v5=v1+(v2-v1)*(galaxy.r[ir]-rr)/CDELT1;
	     float v6=v3+(v4-v3)*(galaxy.r[ir]-rr)/CDELT1;
	     float value=v5+(v6-v5)*(fabs(galaxy.z[iz])-zz)/CDELT2;
	     if(value<0.0) value=0.0;
	     // reverse scale from wavelength to frequency
	     galaxy.ISRF[i].d2[ir][iz].s[galaxy.ISRF[i].n_pgrid-1-inu] = value;
	     
	     // cout<<"r z irr izz rr zz "<< galaxy.r[ir]<<" "<<galaxy.z[iz]<<"   "<<irr<<" "<<izz<<"   "<<rr<<" "<<zz<<" "<< v5+(v6-v5)*(galaxy.z[iz]-zz)/CDELT2<<endl;
	   }  //  inu
	 }  //  iz
       }  //  ir
     }  //  2D
     
     if(galdef.n_spatial_dimensions==3) {              // ==== 3D ====
      
       galaxy.ISRF[i].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid,NAXIS3);
         for(int ix=0; ix<galaxy.n_xgrid; ix++) {  
	   
	   for(int iy=0; iy<galaxy.n_ygrid; iy++) {

	     for(int iz=0; iz<galaxy.n_zgrid; iz++) {

	       float r=sqrt(pow(galaxy.x[ix],2)+pow(galaxy.y[iy],2));
	       int irr=(int)((            r     -CRVAL1) /CDELT1+0.5);//IMOS20060420
	       int izz=(int)((fabs(galaxy.z[iz])-CRVAL2) /CDELT2+0.5);//IMOS20060420
	       if(irr>NAXIS1-2) irr=NAXIS1-2;
	       if(izz>NAXIS2-2) izz=NAXIS2-2;
	       
	       float rr=CRVAL1+irr*CDELT1;
	       float zz=CRVAL2+izz*CDELT2;
	       
	       // cout<<"r z irr izz rr zz "<< galaxy.r[ir]<<" "<<galaxy.z[iz]<<"     "<<irr<<" "<<izz<<"     "<<rr<<" "<<zz<<endl;
	       
	       for(int inu=0; inu<galaxy.ISRF[i].n_pgrid; inu++) {
		 
		 float v1=isrf_in[isrf_index(irr  ,izz  ,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
		 float v2=isrf_in[isrf_index(irr+1,izz  ,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
		 float v3=isrf_in[isrf_index(irr  ,izz+1,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
		 float v4=isrf_in[isrf_index(irr+1,izz+1,inu,i,NAXIS1,NAXIS2,NAXIS3,NAXIS4)];
		 float v5=v1+(v2-v1)*(       r    -rr)/CDELT1;
		 float v6=v3+(v4-v3)*(       r    -rr)/CDELT1;
		 float value=v5+(v6-v5)*(fabs(galaxy.z[iz])-zz)/CDELT2;
		 if(value<0.0) value=0.0;
		 // reverse scale from wavelength to frequency
		 galaxy.ISRF[i].d3[ix][iy][iz].s[galaxy.ISRF[i].n_pgrid-1-inu]= value;
	       }  //  inu   
	     }  //  ix
	   }  //  iy
         }  //  iz
     }  //  3D

     if(galdef.verbose>=10) {

       cout<<"ISRF component "<<i+1<<endl;
       galaxy.ISRF[i].print();
     
     }
   
   }  // isrf component i

   // Create the array of ISRF frequencies
   // using wavelength in microns for axis 3 of input ISRF on log10 scale.
   // Reverse scale so that frequency increases.
   
   galaxy.nu_ISRF = new double[galaxy.ISRF[0].n_pgrid];
   
   // microns -> cm; nu=c/lambda
   for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++)
     galaxy.nu_ISRF[galaxy.ISRF[0].n_pgrid-1-inu]= 
       c/(pow(10.,1.*CRVAL3+inu*CDELT3)*1.0e-4);
   
   cout<<" ISRF frequency grid (Hz):"<<endl;
   for(int inu=0; inu<galaxy.ISRF[0].n_pgrid; inu++) 
     cout<<"inu galaxy.nu_ISRF[inu] "<<inu<<" "<< galaxy.nu_ISRF[inu]<<endl;
   
   delete[] isrf_in;//AWS20010216
      
  } else { 

    // Read in using new format

    const std::string fitsDirectory = configure.fits_directory;
    const std::string isrfFilename = galdef.ISRF_file;
    const std::string filename = fitsDirectory + isrfFilename;

    cout << "Reading ISRF from " << filename << endl;

    //galaxy.fISRFModel = new rf::GalacticRadiationField(filename, true); // For now, no angular information is read in -- TAP 09082007

    rf::GalacticRadiationField rf(filename, true); //(galaxy.fISRFModel);

    cout << "Generating galaxy.ISRF:" << endl;

    // Only three components: stellar + scattered, infrared, CMB

    const unsigned long components = 3;

    galaxy.n_ISRF_components = components; 
    galaxy.ISRF = new Distribution[components];

    // Set up the angular bins

    //const std::valarray<double> &azimuth = rf.GetAzimuthData();
    //const std::valarray<double> &cosZenith = rf.GetCosZenithData();

    //const unsigned long azBins = azimuth.size();
    
    //galaxy.phi_ISRF.resize(azBins);
    //galaxy.phi_ISRF = azimuth;

    //const unsigned long cZBins = cosZenith.size();

    //galaxy.cosTheta_ISRF.resize(cZBins);
    //galaxy.cosTheta_ISRF = cosZenith;

    cout << " galaxy.ISRF initialized with angular bins = "
	 << rf.GetAzimuthData().size() << " " << rf.GetCosZenithData().size() << endl;

    const std::valarray<double>& wl = rf.GetWavelengthData();
    
    // Cheat a bit. Get what the wavelength log(delta) is and extend the 
    // range up to 10000 microns (if it already doesn't go that far). This 
    // ensures we get the full CMB as well.

    const unsigned long rawWlBins = wl.size();

    std::valarray<double> wavelengthData;

    if (wl[rawWlBins-1] < 1e4) {

      const double wlDelta = log10(wl[1]/wl[0]);

      const unsigned long bins = long(log10(1e4/wl[0])/wlDelta) + 1;

      wavelengthData.resize(bins);

      for (unsigned long i = 0; i < bins; ++i)
	wavelengthData[i] = 
	  (i < rawWlBins ? wl[i] : wl[0]*pow(10.0, i*wlDelta));

    } else {

      wavelengthData.resize(rawWlBins);
      wavelengthData = wl;

    }

    const unsigned long wlBins = wavelengthData.size();

    // Create the array of ISRF frequencies using wavelength in microns. 
    // Reverse scale so that frequency increases.
   
    galaxy.nu_ISRF = new double[wlBins];
   
    for (unsigned long i = 0; i < wlBins; ++i)
      galaxy.nu_ISRF[wlBins - 1 - i] = 
	utl::kSpeedOfLight_SI/(wavelengthData[i]*utl::micron/utl::m);
   
    cout << " ISRF frequency grid (Hz): " << endl;

    for (int inu = 0; inu < wlBins; ++inu) 
      cout << "inu galaxy.nu_ISRF[inu] " << inu << " " << galaxy.nu_ISRF[inu] << endl;

    // Some of this is truly, horribly, awful. However, I can't fix it until
    // Galprop undergoes a full re-write. I can't see that happening in the
    // near future -- TAP20072301

    if (galdef.n_spatial_dimensions == 2) { // 2D
      
      galaxy.ISRF[0].init(galaxy.n_rgrid, galaxy.n_zgrid, wlBins);
      galaxy.ISRF[1].init(galaxy.n_rgrid, galaxy.n_zgrid, wlBins);
      galaxy.ISRF[2].init(galaxy.n_rgrid, galaxy.n_zgrid, wlBins);

      cout << " galaxy.ISRF initialized with frequency axis dimension = "
	   << wlBins << endl;

      for (unsigned long i = 0; i < galaxy.n_rgrid; ++i) {

	for (unsigned long j = 0; j < galaxy.n_zgrid; ++j) {

	  const double r = galaxy.r[i], z = fabs(galaxy.z[j]);
	  
	  for (unsigned long k = 0; k < wlBins; ++k) {
	    
	    const double wl = wavelengthData[k];

	    // Scale is reversed from wavelength to frequency

	    const double stellar = 
	      rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::STELLAR);

	    const double scattered = 
	      rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::SCATTERED);

	    const double infrared = 
	      rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::INFRARED);

	    const double cmb = 
	      rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::CMB);
	    
	    galaxy.ISRF[0].d2[i][j].s[galaxy.ISRF[0].n_pgrid - 1 - k] = stellar + scattered;

	    galaxy.ISRF[1].d2[i][j].s[galaxy.ISRF[1].n_pgrid - 1 - k] = infrared;

	    galaxy.ISRF[2].d2[i][j].s[galaxy.ISRF[2].n_pgrid - 1 - k] = cmb;
	     
	    //cout << "Read ISRF: " << i << " " << j << " " << k << " " 
	    // << galaxy.ISRF[0].d2[i][j].s[galaxy.ISRF[0].n_pgrid - 1 - k] <<" " 
	    // << rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::STELLAR) + rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::SCATTERED)  << " " 
	    //<< galaxy.ISRF[1].d2[i][j].s[galaxy.ISRF[0].n_pgrid - 1 - k] << " " 
		// << rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::INFRARED) << " " 
	    //<< galaxy.ISRF[2].d2[i][j].s[galaxy.ISRF[2].n_pgrid - 1 - k] << " " << rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::CMB) << endl;


	  } // wl

	} // z

      } // r

    } // 2D

    //exit(0);

    if (galdef.n_spatial_dimensions == 3) { // 3D

      galaxy.ISRF[0].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, wlBins);
      galaxy.ISRF[1].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, wlBins);
      galaxy.ISRF[2].init(galaxy.n_xgrid, galaxy.n_ygrid, galaxy.n_zgrid, wlBins);

      cout << " galaxy.ISRF initialized with frequency axis dimension = "
	   << wlBins << endl;

      for (unsigned long i = 0; i < galaxy.n_xgrid; ++i) {

	for (unsigned long j = 0; j < galaxy.n_ygrid; ++j) {

	  const double r = sqrt(galaxy.x[i]*galaxy.x[i] + galaxy.y[j]*galaxy.y[j]);

	  for (unsigned long k = 0; k < galaxy.n_zgrid; ++k) {

	    const double z = fabs(galaxy.z[k]);
	  
	    for (unsigned long l = 0; l < wlBins; ++l) {
	    
	      const double wl = wavelengthData[l];

	      // Scale is reversed from wavelength to frequency

	      const double stellar = 
		rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::STELLAR);
	      
	      const double scattered = 
		rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::SCATTERED);
	      
	      const double infrared = 
		rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::INFRARED);
	      
	      const double cmb = 
		rf.GetEnergyDensity(wl, r, z, rf::GalacticRadiationField::CMB);
	    
	      galaxy.ISRF[0].d3[i][j][k].s[galaxy.ISRF[0].n_pgrid - 1 - l] = stellar + scattered;

	      galaxy.ISRF[1].d3[i][j][k].s[galaxy.ISRF[1].n_pgrid - 1 - l] = infrared;
	      
	      galaxy.ISRF[2].d3[i][j][k].s[galaxy.ISRF[2].n_pgrid - 1 - l] = cmb;
	     
	    } // wl

	  } // z

	} // y

      } // x

    } // 3D

  }*/

  for (int i = 0; i < galaxy.n_ISRF_components; ++i) {

    if (galdef.n_spatial_dimensions == 2) {

      for (int iR = 0; iR < galaxy.ISRF->n_rgrid; ++iR)
	for (int iZ = 0; iZ < galaxy.ISRF->n_zgrid; ++iZ) 	  
	  for (int iP = 0; iP < galaxy.ISRF->n_pgrid; ++iP)
	    galaxy.ISRF[i].d2[iR][iZ].s[iP] *= galaxy.fISRFFactors[i]; 

    }

    if (galdef.n_spatial_dimensions == 3) {

      for (int iX = 0; iX < galaxy.ISRF->n_xgrid; ++iX)
	for (int iY = 0; iY < galaxy.ISRF->n_ygrid; ++iY)
	  for (int iZ = 0; iZ < galaxy.ISRF->n_zgrid; ++iZ)
	    for (int iP = 0; iP < galaxy.ISRF->n_pgrid; ++iP)
	      galaxy.ISRF[i].d3[iX][iY][iZ].s[iP] *= galaxy.fISRFFactors[i]; 

    }

  }
      
  INFO("Exit");

  return status;

}
