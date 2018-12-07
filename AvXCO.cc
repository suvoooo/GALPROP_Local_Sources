#include "galprop_classes.h"
#include "galprop_internal.h"
#include "SkymapFitsio.h"

#include <ErrorLogger.h>

#include <iostream>
#include <fstream>


int Galprop::AvXCO(const string& galdefPath,
		   const string& fitsPath,
		   const string& outputPath,
		   const string& outputPrefix,
		   const string& runNumber) {

  if (configure.init(galdefPath, fitsPath, outputPath, outputPrefix)) {

    FATAL("Internal error. Fix data paths!");
    return 1;

  }

  if (galdef.read(configure.fVersion, runNumber, configure.fGaldefDirectory)) {

    FATAL("Internal error. Problem reading from galdef file!");
    return 1;

  } 

	/*
  if ( 0 != create_galaxy() ) {

     FATAL("Internal error. Problem allocating memory.");
     return 1;

  }
	*/

  //Read in the CO gas maps
  read_gas_maps("COR");

  //Set up the los integrators
  std::vector< SM::LOSfunction<double>* > funcs(3);
  funcs[0] = new GasFunction("CO", 0, *this); //Plain CO
  funcs[1] = new GasFunction("H2", 0, *this); //H2, X corrected CO
  funcs[2] = new GasFunction("H2", 1.0, *this); //H2, scaled with radius to get effective radius
  std::vector<double> Rbins(galaxy.R_bins, galaxy.R_bins+galaxy.n_Ring);
  Rbins.push_back(galaxy.R_bins[2*galaxy.n_Ring-1]);
  SM::LOSintegrator<double> losInt(galdef.r_max, galdef.z_min, galdef.z_max, Rbins, Rsun, 0.0, 0.1);

  if (3 == galdef.skymap_format) {
     Skymap<double> avXCOhp(galaxy.hpCOR.Order(), galaxy.hpCOR.getSpectra(), galaxy.hpCOR.Scheme());
     Skymap<double> avRhp  (galaxy.hpCOR.Order(), galaxy.hpCOR.getSpectra(), galaxy.hpCOR.Scheme());

#pragma omp parallel for schedule(dynamic) default(shared)
     for ( int ii = 0; ii < galaxy.hpCOR.Npix(); ++ii ) {
	SM::Coordinate co ( galaxy.hpCOR.pix2ang ( ii ) );
	double l=co.l();
	double b=co.b();

	double dl = 90./galaxy.hpCOR.Nside();
	vector< vector<double> > gas = losInt.integrate(l, b, dl, funcs);

	for ( size_t j = 0; j < gas[0].size(); ++j ) {
	   if ( gas[0][j] != 0 ) {
	      avXCOhp[ii][j] = gas[1][j]/gas[0][j];
	      avRhp  [ii][j] = gas[2][j]/gas[1][j];
	   }
	}
     }

     const std::string filestart = configure.fOutputDirectory + configure.fOutputPrefix;
     const std::string fileend = "healpix_" + galdef.galdef_ID + ".gz";

     SkymapToFits(avXCOhp, filestart + "averageXCO_" + fileend,        "Radius", "kpc");
     SkymapToFits(avRhp,   filestart + "averageXCORadius_" + fileend,  "Radius", "kpc");

     // Calculate the average XCO and the corresponding radius, weighted by the
     // real CO map.
     std::vector<double> avXCO(galaxy.n_Ring, 0.0), avR(galaxy.n_Ring, 0.0);
     for (size_t i_ring = 0; i_ring < galaxy.n_Ring; ++i_ring) {
	double weight(0);
	for ( int ii = 0; ii < galaxy.hpCOR.Npix(); ++ii ) {
	   if ( galaxy.hpCOR[ii][i_ring] != 0 && avXCOhp[ii][i_ring] != 0 ) {
	      avXCO[i_ring] += avXCOhp[ii][i_ring]*galaxy.hpCOR[ii][i_ring];
	      avR  [i_ring] += avRhp  [ii][i_ring]*galaxy.hpCOR[ii][i_ring];
	      weight += galaxy.hpCOR[ii][i_ring];
	   }
	}
	avXCO[i_ring] /= weight;
	avR  [i_ring] /= weight;
     }

     // Possibly store it as a fits file in the future, but a simple text file
     // will due for now
     const std::string filename = filestart + "averageXCO_" + galdef.galdef_ID + ".txt";
     std::fstream ofs(filename.c_str(), std::fstream::trunc | std::fstream::out);
     for (size_t i_ring = 0; i_ring < galaxy.n_Ring; ++i_ring) {
	ofs << avR[i_ring] << "\t" << avXCO[i_ring] << std::endl;
     }
     ofs.close();
  }

  return 0;
} 
