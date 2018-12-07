//Program to change gas rings used in previously calculated GALPROP output
//Needs old and new gas rings to work
//Assumes to be run in the FITS directory
#include <string>
#include <iostream>
#include <sstream>
#include "Skymap.h"

int main(int argc, char * argv[]) {
   //We need new rings, old rings, galdef_id, and gas ring type on the command line
   if (argc < 6) {
      std::cout<<"Usage:  changeGasMaps  oldRings newRings oldgaldef_id newgaldef_id type(H2R or HIR)"<<std::endl;
      return 1;
   }
   std::string ogid(argv[3]);
   std::string ngid(argv[4]);
   std::string type(argv[5]);
   if (type != "H2R" && type != "HIR") {
      std::cout<<"Only types allowed are H2R and HIR (case sensitive)"<<std::endl;
      return 1;
   }

   std::string file("bremss_");
   file += type;
   file += "_ring_1_healpix_";
   file += ogid;

   //Open the first file to get the correct orders of things
   std::cout<<file<<std::endl;
   Skymap<double> bData(file);

   int order = bData.Order();

   //Load the new and old gas rings and calculate their ratio to multiply the
   //skymap
   Skymap<double> oldr(argv[1], order);
   Skymap<double> newr(argv[2], order);

   //Must not divide by 0, so we create a loop
   for (int i = 0; i < oldr.Npix(); ++i) {
      for (int j = 0; j < oldr.nSpectra(); ++j) {
	 if (oldr[i][j] != 0) {
	    newr[i][j] /= oldr[i][j];
	 }
      }
   }

   //Now we will loop over the rings, both bremss and pi0 and multiply them
   //with the ratio
   Skymap<double> pData;
   std::string bprefix("bremss_");
   std::string pprefix("pi0_decay_");
   for (int i = 1; i <= 17; ++i) {
      //Assume the number of rings are 17
      std::ostringstream os("");
      os << type <<"_ring_"<<i<<"_healpix_";
      try {
	 bData.load(bprefix+os.str()+ogid);
	 for (int j = 0; j < bData.Npix(); ++j) {
	    bData[j] *= newr[j][i-1];
	 }
	 bData.write(bprefix+os.str()+ngid);
      } catch (...) {}
      try {
	 pData.load(pprefix+os.str()+ogid);
	 for (int j = 0; j < bData.Npix(); ++j) {
	    pData[j] *= newr[j][i-1];
	 }
	 pData.write(pprefix+os.str()+ngid);
      } catch (...) {}
   }

   return 0;
}
