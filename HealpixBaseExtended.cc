#include "HealpixBaseExtended.h"

#include "config.h"
#ifdef HAVE_LSCONSTANTS_H
   #include "lsconstants.h"
#endif

#include <vector>
#include <iterator>
#include <algorithm>
#include <set>
#include <cmath>

int Healpix_Base::ring_above (double z) const
{
   double az=fabs(z);
   if (az>twothird) // polar caps
   {
      int iring = int(nside_*sqrt(3*(1-az)));
      return (z>0) ? iring : 4*nside_-iring-1;
   }
   else // ----- equatorial region ---------
      return int(nside_*(2-1.5*z));
}

//Select pixels from a region
std::set<int> HealpixBaseExtended::regionToPixels(const SkyRegion & region) const {
   //The vector to store the pixel values
   std::set<int> pixSet;
   std::vector<int> listpix;
   //Check the region type
   switch (region.type()) {
      case DISC:
	 //Use the query disk
	 if (region.radius() > 0){
   	    query_disc(region.center().healpixAng(), region.radius() * M_PI/180, listpix);
	    std::copy(listpix.begin(), listpix.end(), std::inserter(pixSet, pixSet.begin()));
	 }
	 break;
      case RECTANGLE:
	 if (region.deltal() > 0 && region.deltab() > 0) {
   	    //We need to split the region up if necessary and query them one by
   	    //one
	    double dphi = region.deltal()/2 * M_PI/180;
	    double dtheta = region.deltab()/2 * M_PI/180;
	    
	    //
	    //Now we need to check for overshooting in the theta direction
	    //North pole
	    if (region.center().healpixAng().theta - dtheta < 0){
	       //We go over the pole, add the section on the other side
	       double lowTheta = 0 - (region.center().healpixAng().theta - dtheta); 
	       double highTheta = 0;
	       double lowPhi = region.center().healpixAng().phi - dphi + M_PI;
	       double highPhi = region.center().healpixAng().phi + dphi + M_PI;
	       if (highPhi > 2*M_PI) highPhi -= 2*M_PI;
	       std::cout<<"Overshooding north pole "<<lowTheta<<", "<<highTheta<<", "<<lowPhi<<", "<<highPhi<<std::endl;
	       query_rectangle(pointing(lowTheta, lowPhi), pointing(highTheta, highPhi), listpix);
	       std::copy(listpix.begin(), listpix.end(), std::inserter(pixSet, pixSet.begin()));
	    }
	    //South pole
	    if (region.center().healpixAng().theta + dtheta > M_PI){
	       //We go over the pole, add the section on the other side
	       double lowTheta = M_PI;
	       double highTheta = 2*M_PI - (region.center().healpixAng().theta + dtheta);
	       double lowPhi = region.center().healpixAng().phi - dphi + M_PI;
	       double highPhi = region.center().healpixAng().phi + dphi + M_PI;
	       if (highPhi > 2*M_PI) highPhi -= 2*M_PI;
	       std::cout<<"Overshooding south pole "<<lowTheta<<", "<<highTheta<<", "<<lowPhi<<", "<<highPhi<<std::endl;
	       query_rectangle(pointing(lowTheta, lowPhi), pointing(highTheta, highPhi), listpix);
	       std::copy(listpix.begin(), listpix.end(), std::inserter(pixSet, pixSet.begin()));
	    }
	    //The rest, things that don't loop
	    double lowTheta = std::min(M_PI, region.center().healpixAng().theta+dtheta);
	    double highTheta = std::max(0., region.center().healpixAng().theta-dtheta);
	    double lowPhi = region.center().healpixAng().phi - dphi;
	    double highPhi = region.center().healpixAng().phi + dphi;
	    //Fix overshoot
	    if (highPhi > 2*M_PI) highPhi -= 2*M_PI;
	    if (lowPhi < 0) lowPhi += 2*M_PI;
	    query_rectangle(pointing(lowTheta, lowPhi), pointing(highTheta, highPhi), listpix);
	    std::copy(listpix.begin(), listpix.end(), std::inserter(pixSet, pixSet.begin()));
	 }
   }
   return pixSet;
}

void HealpixBaseExtended::query_pixel(const Healpix_Base &hp, int pixel, std::vector<int> &listpix) const {
   listpix.clear();
   //Find the difference in order between the two bases
   int dOrder = Order() - hp.Order();
   //If dorder is <= 0, then the input is finer binned than the output and we
   //only have a single pixel
   if (dOrder <= 0) {
      listpix.push_back(ang2pix(hp.pix2ang(pixel)));
   }else{
      //Now we must loop over all the pixels covered by the bigger pixel
      //First find its index in nested
      int nestpix = hp.Scheme() == NEST ? pixel : hp.ring2nest(pixel);
      //Find the number of pixels under each pixel
      int npix = (1<<dOrder)*(1<<dOrder);
      listpix.resize(npix);
      //loop over all the pixels, in nest scheme
      for (int p = 0; p < npix; ++p) {
	 listpix[p] = Scheme() == NEST ? p+nestpix*npix : nest2ring(p+nestpix*npix);
      }
   }

}

void HealpixBaseExtended::query_rectangle(const pointing &pointingll, const pointing &pointingur, std::vector<int> &listpix) const {
   //Clear all values from listpix
   listpix.clear();
   //Begin by finding the rings we have to iterate over
   std::vector<int> rings;
   //Check if the lower limit is above the upper limit
   if (pointingll.theta > pointingur.theta) {  //The normal way
      int ur = ring_above(cos(pointingur.theta))+1;
      int lr = ring_above(cos(pointingll.theta));
      //Use the generator algorithm to insert the rings.
      rings.resize(lr-ur+1);
      std::generate(rings.begin(), rings.end(), counter(ur-1));
   } else {  //Looping from north to south
      int ur = ring_above(cos(pointingll.theta));
      int lr = ring_above(cos(pointingur.theta))+1;
      rings.resize(ur+(4*nside_-1)-lr+1); //4*nside_-1 is the number of rings of the skymap
      std::generate(rings.begin(), rings.begin()+ur, counter(0));
      std::generate(rings.begin()+ur, rings.end(), counter(lr-1));
   }

   //Now we have to deal with the rings, selecting the correct pixels
   //The ring_info method only works for RING healpix bases, so we create one
   Healpix_Base hp(order_, RING);
   //First check if left side is really on the left side
   int startpix, ringpix;
   double theta;
   bool shifted;
   if (pointingll.phi < pointingur.phi) { //The normal way, no looping
      for (int i = 0; i < int(rings.size()); ++i) {
	 hp.get_ring_info2(rings[i], startpix, ringpix, theta, shifted);
	 //Find the location of the first pixel within the boundary.
	 //The width of each pixel is 2\pi/ringpix
	 double dl = 2.*pi/ringpix;
	 int lp = int(ceil((pointingll.phi - shifted*dl/2)/dl));
	 int rp = std::min(int((pointingur.phi - shifted*dl/2)/dl), ringpix);
	 listpix.reserve(listpix.size()+rp-lp+1);
	 std::generate_n(std::back_inserter(listpix), rp-lp+1, counter(startpix + lp-1));
	 //In the case the right boundary is 2*pi, the ring is
	 //not shifted and the left boundary is not 0, we must
	 //add the first pixel of the ring.
	 if (pointingll.phi != 0 && pointingur.phi == 2*pi && !shifted)
	    listpix.push_back(startpix);
      }
   } else { //Looping over 2 \pi
      for (int i = 0; i < int(rings.size()); ++i) {
	 hp.get_ring_info2(rings[i], startpix, ringpix, theta, shifted);
	 //Find the location of the first pixel within the boundary.
	 //The width of each pixel is 2\pi/ringpix
	 double dl = 2.*pi/ringpix;
	 int lp = int(floor((pointingur.phi - shifted*dl/2)/dl));
	 int rp = int(ceil((pointingll.phi - shifted*dl/2)/dl));
	 listpix.reserve(listpix.size()+lp+1+ringpix-rp);
	 std::generate_n(std::back_inserter(listpix), lp+1, counter(startpix-1));
	 std::generate_n(std::back_inserter(listpix), ringpix-rp, counter(startpix+rp-1));
      }
   }
   //Now we must check if the actual map is in a NESTED scheme and fix the list
   if (scheme_ == NEST) {
      std::transform(listpix.begin(), listpix.end(), listpix.begin(), std::bind1st(std::mem_fun(&Healpix_Base::ring2nest), this));
   }
}
