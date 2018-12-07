#include "SkySelection.h"
#include "Region.h"
#include "HealpixBaseExtended.h"
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <cctype>

SkySelection::SkySelection(const std::string & filename) {
   addFromFile(filename);
}

void SkySelection::addFromFile(const std::string & filename) {
   //Open the file
   std::ifstream is(filename.c_str());

   SkyRegion region;
   std::string incl, buffer;
   //Read the file and add to the vector
   while ( is.good() ) {
      //Read one line at a time
      std::getline(is, buffer);
      //Continue if the buffer is empty
      if (buffer == "")
	 continue;
      //Create a input stream from the buffer
      std::istringstream iss(buffer);
      //Read the ADD or REMOVE keywords
      iss >> std::skipws >> incl;
      //Read the region
      iss >> region;
      //Skip the line if something went wrong
      if (iss.fail() || iss.bad() ){
	 std::cerr<<"Failed parsing line \""<<buffer<<"\""<<std::endl;
	 continue;
      }
      //Check the inclusion keyword
      std::transform(incl.begin(), incl.end(), incl.begin(), toupper);
      if (incl == "ADD" ) {
	 addRegion(region);
      }else if (incl == "REMOVE" ) {
	 removeRegion(region);
      } else {
	 std::cerr<<"First word in the line has to be either ADD or REMOVE"<<std::endl;
	 std::cerr<<"Skipping line \""<<buffer<<"\""<<std::endl;
	 continue;
      }
   }

   if (fselection.size() == 0) {
      std::cerr<<"No region added after reading file "<<filename<<std::endl;
   }
}

void SkySelection::addRegion(const SkyRegion & region){
   fselection.push_back(std::pair<Inclusion,SkyRegion>(ADD,region));
}

void SkySelection::removeRegion(const SkyRegion & region){
   fselection.push_back(std::pair<Inclusion,SkyRegion>(REMOVE,region));
}

std::set<int> SkySelection::pixelsForHealpix(const HealpixBaseExtended &hp) const{
   std::set<int> pixels;
   //Loop through the selection and either add or remove the selection
   for (size_t i = 0; i < fselection.size(); ++i){
      if (fselection[i].first == ADD) {
	 std::set<int> pix = hp.regionToPixels(fselection[i].second);
	 std::set<int>::iterator it = pix.begin();
	 for ( ; it != pix.end(); ++it) { 
	    pixels.insert(*it);
	 }
      } else {
	 std::set<int> pix = hp.regionToPixels(fselection[i].second);
	 std::set<int>::iterator it = pix.begin();
	 for ( ; it != pix.end(); ++it) { 
	    pixels.erase(*it);
	 }
      }
   }
   return pixels;
}
