
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Distribution.cc *                             galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include"Distribution.h"

#include <limits>

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
/*	
Distribution::~Distribution()
{
   cout<<"~Distribution\n";
	 
   if(n_spatial_dimensions==3)
   {
      for(int ix=0; ix<n_xgrid; ix++)
      {
         for(int iy=0; iy<n_ygrid; iy++)
         {
            for(int iz=0; iz<n_zgrid; iz++) delete[]d3[ix][iy][iz].s;  
	    delete[]d3[ix][iy]; 
         }
	 delete[]d3[ix];
      }
      delete[]d3;  
   } // 3D
};
*/
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution::Distribution(const Distribution &old){
  arrays_assigned=false;
  n_spatial_dimensions=0;
  
  if (2 == old.n_spatial_dimensions) {

    init(old.n_rgrid, old.n_zgrid, old.n_pgrid);
    
    for (int ir = 0; ir < old.n_rgrid; ++ir)
      for (int iz = 0; iz < old.n_zgrid; ++iz)
	for (int ip = 0; ip < old.n_pgrid; ++ip)
	  d2[ir][iz].s[ip] = old.d2[ir][iz].s[ip];
  }
  
  if (3 == old.n_spatial_dimensions) {

    init(old.n_xgrid, old.n_ygrid, old.n_zgrid, old.n_pgrid);
    
    for (int iy = 0; iy < n_ygrid; ++iy)
      for (int ix = 0; ix < n_xgrid; ++ix)
	for (int iz = 0; iz < n_zgrid; ++iz)
	  for (int ip = 0; ip < n_pgrid; ++ip)
	    d3[ix][iy][iz].s[ip] = old.d3[ix][iy][iz].s[ip];
   
  }

}

Distribution::Distribution() {

  n_spatial_dimensions = 0;
  arrays_assigned = false;

}

Distribution::~Distribution() {

  delete_array();

}

Distribution::Distribution(int n_rgrid_, int n_zgrid_, int n_pgrid_) {

  n_spatial_dimensions=2;
  n_rgrid=n_rgrid_;
  n_zgrid=n_zgrid_;
  n_pgrid=n_pgrid_;
  
  d2=new Spectrum *[n_rgrid];
  for(int ir=0; ir<n_rgrid; ++ir) {
    d2[ir]=new Spectrum[n_zgrid];
    for(int iz=0; iz<n_zgrid; ++iz) {

      d2[ir][iz].s=new double [n_pgrid];//AWS20050624
      for(int ip=0; ip<n_pgrid; ip++) d2[ir][iz].s[ip]=0.;
      
    }
   
  }
  //cout<<"Distribution:generated 2D array of spectra "
  //<<n_rgrid<<" "<<n_zgrid<<" "<<n_pgrid<<endl; 
  arrays_assigned = true;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution::
Distribution(int n_xgrid_, int n_ygrid_, int n_zgrid_, int n_pgrid_) {

  n_spatial_dimensions=3;
  n_xgrid=n_xgrid_;
  n_ygrid=n_ygrid_;
  n_zgrid=n_zgrid_;
  n_pgrid=n_pgrid_;
  
  d3=new Spectrum **[n_xgrid];

  for (int ix = 0; ix < n_xgrid; ++ix) {

    d3[ix]=new Spectrum *[n_ygrid];
    for(int iy=0; iy<n_ygrid; ++iy) {

      d3[ix][iy]=new Spectrum [n_zgrid];
      for(int iz=0; iz<n_zgrid; ++iz) {
	d3[ix][iy][iz].s=new double [n_pgrid];//AWS20050624
	for(int ip=0; ip<n_pgrid; ++ip) 
	  d3[ix][iy][iz].s[ip]=0.;
         
      }
      
    }
   
  }
  //cout<<"Distribution:generated 3D array of spectra "
  // <<n_xgrid<<" "<<n_ygrid<<" "<<n_zgrid<<" "<<n_pgrid<<endl; 
  
  arrays_assigned = true;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void Distribution::delete_array() {

  if (arrays_assigned) {
    //cout<<"Distribution delete_array\n";
    
    if (2 == n_spatial_dimensions) {

      for(int ir=0; ir<n_rgrid; ir++) {

	for(int iz=0;iz<n_zgrid;iz++) delete[]d2[ir][iz].s;                       
	delete[]d2[ir];  
	
      }	    

      delete[]d2;  
      
    }

    if (3 == n_spatial_dimensions) { 
      
      for(int ix=0; ix<n_xgrid; ix++) {

	for(int iy=0; iy<n_ygrid; iy++) {

	  for(int iz=0; iz<n_zgrid; iz++)
	    delete[]d3[ix][iy][iz].s;  	        
	  
	  delete[]d3[ix][iy];
	  
	}
	
	delete[]d3[ix];
	
      }	    
      
      delete[]d3;  
      
    }
    
    arrays_assigned=false;
    n_spatial_dimensions=0;
   
  }

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void Distribution::init(int n_rgrid_, int n_zgrid_, int n_pgrid_) {
   
  delete_array();
  n_spatial_dimensions=2;
  n_rgrid=n_rgrid_;
  n_zgrid=n_zgrid_;
  n_pgrid=n_pgrid_;
  
  d2=new Spectrum *[n_rgrid];
  for(int ir=0; ir<n_rgrid; ir++) {

    d2[ir]=new Spectrum[n_zgrid];
    for(int iz=0; iz<n_zgrid; iz++) {

      d2[ir][iz].s=new double [n_pgrid];//AWS20050624
      for(int ip=0; ip<n_pgrid; ip++) d2[ir][iz].s[ip]=0.;
      
    }
   
  }
  //cout<<"Distribution.init:generated 2D array of spectra "
  // <<n_rgrid<<" "<<n_zgrid<<" "<<n_pgrid<<endl; 
  arrays_assigned=true;
  //return;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void 
Distribution::init(int n_xgrid_, int n_ygrid_, int n_zgrid_, int n_pgrid_) {

  delete_array();
  n_spatial_dimensions=3;
  n_xgrid=n_xgrid_;
  n_ygrid=n_ygrid_;
  n_zgrid=n_zgrid_;
  n_pgrid=n_pgrid_;
  
  d3=new Spectrum **[n_xgrid];
  for(int ix=0; ix<n_xgrid; ix++) {

    d3[ix]=new Spectrum *[n_ygrid];
    for(int iy=0; iy<n_ygrid; iy++) {

      d3[ix][iy]=new Spectrum [n_zgrid];
      for(int iz=0; iz<n_zgrid; iz++) {

	d3[ix][iy][iz].s=new double [n_pgrid];//AWS20050624
	for(int ip=0; ip<n_pgrid; ip++) d3[ix][iy][iz].s[ip]=0.;
   
      }
    
    }
  
  }
  //cout<<"Distribution.init:generated 3D array of spectra "
  //<<n_xgrid<<" "<<n_ygrid<<" "<<n_zgrid<<" "<<n_pgrid<<endl; 
  arrays_assigned=true;
  //return;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void Distribution::print() {

  //cout<<"Distribution.print: n_spatial_dimensions="<<n_spatial_dimensions<<endl;
  
  if (n_spatial_dimensions==2) {

    cout<<"Distribution with "<<n_spatial_dimensions<<" spatial dimensions: "<<endl;
    for(int ip=0; ip<n_pgrid; ip++) {

      cout<<"ip="<<ip<<endl;
      for(int iz=0; iz<n_zgrid; iz++) {

	cout<<iz<<" ";
	for(int ir=0; ir<n_rgrid; ir++) cout<<d2[ir][iz].s[ip]<<" ";
	cout<<endl;
      }
      
      cout<<endl;
      
    }
   
  }

  if (n_spatial_dimensions==3) {

    cout<<"Distribution with "<<n_spatial_dimensions<<" spatial dimensions: "<<endl;
    for(int ip=0; ip<n_pgrid; ip++) {

      for(int iz=0; iz<n_zgrid; iz++) {

	cout<<"ip="<<ip<<" iz="<<iz<<endl;
	for(int iy=0; iy<n_ygrid; iy++) {

	  cout<<iy<<" ";
	  for(int ix=0; ix<n_xgrid; ix++) cout<<d3[ix][iy][iz].s[ip]<<" ";
	  cout<<endl;
	
	}
        
	cout<<endl;
        
      }
      
      cout<<endl;
      
    }
   
  }

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

void Distribution::print(double units) {

  cout<<"Distribution.print: n_spatial_dimensions="<<n_spatial_dimensions<<endl;
  
  if (n_spatial_dimensions==2)
    {
      cout<<"Distribution with "<<n_spatial_dimensions<<" spatial dimensions: "<<endl;
      cout<<"units="<<units<<endl;
      for(int ip=0; ip<n_pgrid; ip++)
	{
	  cout<<"ip="<<ip<<endl;
	  for(int iz=0; iz<n_zgrid; iz++)
	    {
	      cout<<iz<<" ";
	      for(int ir=0; ir<n_rgrid; ir++) cout<<int(d2[ir][iz].s[ip]/units)<<" ";
	      cout<<endl;
	    }
	  cout<<endl;
	}
    }
  
  if (n_spatial_dimensions==3)
    {
      cout<<"Distribution with "<<n_spatial_dimensions<<" spatial dimensions: "<<endl;
      cout<<"units="<<units<<endl;
      for(int ip=0; ip<n_pgrid; ip++)
	{
	  for(int iz=0;iz<n_zgrid;iz++)
	    {
	      cout<<"ip="<<ip<<" iz="<<iz<<endl;
	      for(int iy=0; iy<n_ygrid; iy++)
		{
		  cout<<iy<<" ";
		  for(int ix=0; ix<n_xgrid; ix++) cout<<int(d3[ix][iy][iz].s[ip]/units)<<" ";
		  cout<<endl;
		}
	      cout<<endl;
	    }
	  cout<<endl;
	}
    }
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double Distribution::max() { //AWS20050624

  double max_= -1.*(numeric_limits<double>::max() - 1.);//-1.e30;    //AWS20050624
  
  if(n_spatial_dimensions==2)
      for(int ir=0; ir<n_rgrid; ir++)
         for(int iz=0; iz<n_zgrid; iz++)
	    for(int ip=0; ip<n_pgrid; ip++)
	       if(d2[ir][iz].s[ip]>max_) max_=d2[ir][iz].s[ip];

   if(n_spatial_dimensions==3)
      for(int ix=0; ix<n_xgrid; ix++)
         for(int iy=0; iy<n_ygrid; iy++)
            for(int iz=0; iz<n_zgrid; iz++)
	        for(int ip=0; ip<n_pgrid; ip++)
	           if(d3[ix][iy][iz].s[ip]>max_) max_=d3[ix][iy][iz].s[ip];
   return max_;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution & Distribution::operator=(double const2) 
{
   if(n_spatial_dimensions==2)
      for(int ir=0; ir<n_rgrid; ir++)
         for(int iz=0; iz<n_zgrid; iz++)
	    for(int ip=0; ip<n_pgrid; ip++) d2[ir][iz].s[ip]=const2;

   if(n_spatial_dimensions==3)
      for(int ix=0; ix<n_xgrid; ix++)
         for(int iy=0; iy<n_ygrid; iy++)
            for(int iz=0; iz<n_zgrid; iz++)
	        for(int ip=0; ip<n_pgrid; ip++) d3[ix][iy][iz].s[ip]=const2;
   return *this;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution & Distribution::operator=(const Distribution &dist)
{
  if (2 == dist.n_spatial_dimensions) {

    init(dist.n_rgrid, dist.n_zgrid, dist.n_pgrid);
    
    for (int ir = 0; ir < dist.n_rgrid; ++ir)
      for (int iz = 0; iz < dist.n_zgrid; ++iz)
	for (int ip = 0; ip < dist.n_pgrid; ++ip)
	  d2[ir][iz].s[ip] = dist.d2[ir][iz].s[ip];
  }
  
  if (3 == dist.n_spatial_dimensions) {

    init(dist.n_xgrid, dist.n_ygrid, dist.n_zgrid, dist.n_pgrid);
    
    for (int iy = 0; iy < n_ygrid; ++iy)
      for (int ix = 0; ix < n_xgrid; ++ix)
	for (int iz = 0; iz < n_zgrid; ++iz)
	  for (int ip = 0; ip < n_pgrid; ++ip)
	    d3[ix][iy][iz].s[ip] = dist.d3[ix][iy][iz].s[ip];
   
  }
   return *this;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution Distribution::operator+(double const2)
{
   if(n_spatial_dimensions==2)
   {
      Distribution tmp(n_rgrid,n_zgrid,n_pgrid);
      for(int ir=0; ir<n_rgrid; ir++)
         for(int iz=0; iz<n_zgrid; iz++)
            for(int ip=0; ip<n_pgrid; ip++)
	       tmp.d2[ir][iz].s[ip]=d2[ir][iz].s[ip]+const2;
      return tmp;
   }

   if(n_spatial_dimensions==3)
   {
      Distribution tmp( n_xgrid,n_ygrid,n_zgrid,n_pgrid); 
      for(int ix=0; ix<n_xgrid; ix++)
         for(int iy=0; iy<n_ygrid; iy++)
             for(int iz=0; iz<n_zgrid; iz++)
	        for(int ip=0; ip<n_pgrid; ip++)
	           tmp.d3[ix][iy][iz].s[ip]=d3[ix][iy][iz].s[ip]+const2;
      return tmp;
   }

   // to avoid compiler warning 'no return statement' //AWS20050627
   Distribution tmp;                                  //AWS20050627
   return tmp;                                        //AWS20050627

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution & Distribution::operator+=(double const2)
{
   if(n_spatial_dimensions==2)       	 
      for(int ir=0; ir<n_rgrid; ir++)
         for(int iz=0; iz<n_zgrid; iz++)
	    for(int ip=0; ip<n_pgrid; ip++) d2[ir][iz].s[ip]+=const2;

   if(n_spatial_dimensions==3)       	 
      for(int ix=0; ix<n_xgrid; ix++)
         for(int iy=0; iy<n_ygrid; iy++)
            for(int iz=0; iz<n_zgrid; iz++)
	       for(int ip=0; ip<n_pgrid; ip++) d3[ix][iy][iz].s[ip]+=const2;
   return *this;
}                    //     Schildt p.392

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution Distribution::operator*(double const2)
{
   if(n_spatial_dimensions==2)
   {
      Distribution tmp(n_rgrid,n_zgrid,n_pgrid);
      for(int ir=0; ir<n_rgrid; ir++)
         for(int iz=0; iz<n_zgrid; iz++)
            for(int ip=0; ip<n_pgrid; ip++)
               tmp.d2[ir][iz].s[ip]=d2[ir][iz].s[ip]*const2;
      return tmp;
   }

   if(n_spatial_dimensions==3)
   {
      Distribution tmp(n_xgrid,n_ygrid,n_zgrid,n_pgrid);
      for(int ix=0; ix<n_xgrid; ix++)
         for(int iy=0; iy<n_ygrid; iy++)
            for(int iz=0; iz<n_zgrid; iz++)
	       for(int ip=0; ip<n_pgrid; ip++)
	          tmp.d3[ix][iy][iz].s[ip]=d3[ix][iy][iz].s[ip]*const2;
      return tmp;
   }

   // to avoid compiler warning 'no return statement' //AWS20050627
   Distribution tmp;                                  //AWS20050627
   return tmp;                                        //AWS20050627
}  	 

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution &Distribution::operator*=(double const2)
{
   if(n_spatial_dimensions==2)      	 
      for(int ir=0; ir<n_rgrid; ir++)
         for(int iz=0; iz<n_zgrid; iz++)
            for(int ip=0; ip<n_pgrid; ip++) d2[ir][iz].s[ip]*=const2;

   if(n_spatial_dimensions==3)       	 
      for(int ix=0; ix<n_xgrid; ix++)
         for(int iy=0; iy<n_ygrid; iy++)
            for(int iz=0; iz<n_zgrid; iz++)
	        for(int ip=0; ip<n_pgrid; ip++) d3[ix][iy][iz].s[ip]*=const2;
   return *this;
}  

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution &Distribution::operator/=(double const2)
{
   if(n_spatial_dimensions==2)      	 
      for(int ir=0; ir<n_rgrid; ir++)
         for(int iz=0; iz<n_zgrid; iz++)
            for(int ip=0; ip<n_pgrid; ip++) d2[ir][iz].s[ip]/=const2;

   if(n_spatial_dimensions==3)       	 
      for(int ix=0; ix<n_xgrid; ix++)
         for(int iy=0; iy<n_ygrid; iy++)
            for(int iz=0; iz<n_zgrid; iz++)
	       for(int ip=0; ip<n_pgrid; ip++) d3[ix][iy][iz].s[ip]/=const2;
   return *this;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution &Distribution::operator+=(const Distribution &dist)
{
   if(n_spatial_dimensions==2)       	 
      for(int ir=0; ir<n_rgrid; ir++)
         for(int iz=0; iz<n_zgrid; iz++)
	    for(int ip=0; ip<n_pgrid; ip++) d2[ir][iz].s[ip]+=dist.d2[ir][iz].s[ip];

   if(n_spatial_dimensions==3)       	 
      for(int ix=0; ix<n_xgrid; ix++)
         for(int iy=0; iy<n_ygrid; iy++)
            for(int iz=0; iz<n_zgrid; iz++)
	       for(int ip=0; ip<n_pgrid; ip++)
	          d3[ix][iy][iz].s[ip]+=dist.d3[ix][iy][iz].s[ip];
   return *this;
}                    //     Schildt p.392

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution &Distribution::operator*=(const Distribution &dist)
{
   if(n_spatial_dimensions==2)
      for(int ir=0; ir<n_rgrid; ir++)
         for(int iz=0; iz<n_zgrid; iz++)
            for(int ip=0; ip<n_pgrid; ip++) d2[ir][iz].s[ip]*=dist.d2[ir][iz].s[ip];

   if(n_spatial_dimensions==3)
      for(int ix=0; ix<n_xgrid; ix++)
         for(int iy=0; iy<n_ygrid; iy++)
            for(int iz=0; iz<n_zgrid; iz++)
	       for(int ip=0; ip<n_pgrid; ip++)
	          d3[ix][iy][iz].s[ip]*=dist.d3[ix][iy][iz].s[ip];
   return *this;
}                    //     Schildt p.392

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution &Distribution::operator/=(const Distribution &dist) 
{
   if(n_spatial_dimensions==2)
      for(int ir=0; ir<n_rgrid; ir++)
         for(int iz=0; iz<n_zgrid; iz++)
            for(int ip=0; ip<n_pgrid; ip++) d2[ir][iz].s[ip]/=dist.d2[ir][iz].s[ip];

   if(n_spatial_dimensions==3)
      for(int ix=0; ix<n_xgrid; ix++)
         for(int iy=0; iy<n_ygrid; iy++)
            for(int iz=0; iz<n_zgrid; iz++)
	       for(int ip=0; ip<n_pgrid; ip++)
	          d3[ix][iy][iz].s[ip]/=dist.d3[ix][iy][iz].s[ip];
   return *this;
}                    //     Schildt p.392

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
//  following cases have to be non-member functions since they have 2 arguments
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution operator+(double const2,const Distribution &dist)
{
   if(dist.n_spatial_dimensions==2)
   {
      Distribution tmp(dist.n_rgrid,dist.n_zgrid,dist.n_pgrid);
      for(int ir=0; ir<dist.n_rgrid; ir++)
         for(int iz=0; iz<dist.n_zgrid; iz++)
	    for(int ip=0; ip<dist.n_pgrid; ip++)
	       tmp.d2[ir][iz].s[ip]=const2+dist.d2[ir][iz].s[ip];
      return tmp;
   }

   if(dist.n_spatial_dimensions==3)
   {
      Distribution tmp(dist.n_xgrid,dist.n_ygrid,dist.n_zgrid,dist.n_pgrid);
      for(int ix=0; ix<dist.n_xgrid; ix++)
         for(int iy=0; iy<dist.n_ygrid; iy++)
            for(int iz=0; iz<dist.n_zgrid; iz++)
	       for(int ip=0; ip<dist.n_pgrid; ip++)
	          tmp.d3[ix][iy][iz].s[ip]=const2+dist.d3[ix][iy][iz].s[ip];
      return tmp;
   }

   // to avoid compiler warning 'no return statement' //AWS20050627
   Distribution tmp;                                  //AWS20050627
   return tmp;                                        //AWS20050627

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution operator*(double const2,const Distribution &dist)
{
   if(dist.n_spatial_dimensions==2)
   {
      Distribution tmp(dist.n_rgrid,dist.n_zgrid,dist.n_pgrid);
      for(int ir=0; ir<dist.n_rgrid; ir++)
         for(int iz=0; iz<dist.n_zgrid; iz++)
	    for(int ip=0; ip<dist.n_pgrid; ip++)
	       tmp.d2[ir][iz].s[ip]=const2 * dist.d2[ir][iz].s[ip];
      return tmp;
   }

   if(dist.n_spatial_dimensions==3)
   {
      Distribution tmp(dist. n_xgrid,dist.n_ygrid,dist.n_zgrid,dist.n_pgrid); 
      for(int ix=0; ix<dist.n_xgrid; ix++)
         for(int iy=0; iy<dist.n_ygrid; iy++)
            for(int iz=0; iz<dist.n_zgrid; iz++)
	       for(int ip=0; ip<dist.n_pgrid; ip++)
	          tmp.d3[ix][iy][iz].s[ip]=const2 * dist.d3[ix][iy][iz].s[ip];
      return tmp;
   }

   // to avoid compiler warning 'no return statement' //AWS20050627
   Distribution tmp;                                  //AWS20050627
   return tmp;                                        //AWS20050627

}  

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution operator+(const Distribution &dist,const Distribution &dist2) 
{
   if(dist.n_spatial_dimensions==2)
   {
      Distribution tmp(dist.n_rgrid,dist.n_zgrid,dist.n_pgrid); 
      for(int ir=0; ir<dist.n_rgrid; ir++)
         for(int iz=0; iz<dist.n_zgrid; iz++)
            for(int ip=0; ip<dist.n_pgrid; ip++)
               tmp.d2[ir][iz].s[ip]=dist.d2[ir][iz].s[ip]+dist2.d2[ir][iz].s[ip];
      return tmp;
   }

   if(dist.n_spatial_dimensions==3)
   {
      Distribution tmp(dist. n_xgrid,dist.n_ygrid,dist.n_zgrid,dist.n_pgrid);
      for(int ix=0; ix<dist.n_xgrid; ix++)
         for(int iy=0; iy<dist.n_ygrid; iy++)
            for(int iz=0; iz<dist.n_zgrid; iz++)
	       for(int ip=0; ip<dist.n_pgrid; ip++)
	          tmp.d3[ix][iy][iz].s[ip]=dist.d3[ix][iy][iz].s[ip]
                     + dist2.d3[ix][iy][iz].s[ip];
      return tmp;
   }

   // to avoid compiler warning 'no return statement' //AWS20050627
   Distribution tmp;                                  //AWS20050627
   return tmp;                                        //AWS20050627

}  
	 
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution operator*(const Distribution &dist,const Distribution &dist2)
{
   if(dist.n_spatial_dimensions==2)
   {
      Distribution tmp(dist.n_rgrid,dist.n_zgrid,dist.n_pgrid);
      for(int ir=0; ir<dist.n_rgrid; ir++)
         for(int iz=0; iz<dist.n_zgrid; iz++)
            for(int ip=0; ip<dist.n_pgrid; ip++)
               tmp.d2[ir][iz].s[ip]=dist.d2[ir][iz].s[ip] * dist2.d2[ir][iz].s[ip];
      return tmp;
   }

   if(dist.n_spatial_dimensions==3)
   {
      Distribution tmp(dist. n_xgrid,dist.n_ygrid,dist.n_zgrid,dist.n_pgrid);
      for(int ix=0; ix<dist.n_xgrid; ix++)
         for(int iy=0; iy<dist.n_ygrid; iy++)
            for(int iz=0; iz<dist.n_zgrid; iz++)
	       for(int ip=0; ip<dist.n_pgrid; ip++)
	          tmp.d3[ix][iy][iz].s[ip]=dist.d3[ix][iy][iz].s[ip] 
                     * dist2.d3[ix][iy][iz].s[ip];
      return tmp;
   }

   // to avoid compiler warning 'no return statement' //AWS20050627
   Distribution tmp;                                  //AWS20050627
   return tmp;                                        //AWS20050627

}  

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Distribution operator/(const Distribution & dist,const Distribution &dist2) 
{
   if(dist.n_spatial_dimensions==2)
   {
      Distribution tmp(dist.n_rgrid,dist.n_zgrid,dist.n_pgrid);
      for(int ir=0; ir<dist.n_rgrid; ir++)
         for(int iz=0; iz<dist.n_zgrid; iz++)
            for(int ip=0; ip<dist.n_pgrid; ip++)
               tmp.d2[ir][iz].s[ip]=dist.d2[ir][iz].s[ip] / dist2.d2[ir][iz].s[ip];
      return tmp;
   }

   if(dist.n_spatial_dimensions==3)
   {
      Distribution tmp(dist. n_xgrid,dist.n_ygrid,dist.n_zgrid,dist.n_pgrid); 
      for(int ix=0; ix<dist.n_xgrid; ix++)
         for(int iy=0; iy<dist.n_ygrid; iy++)
            for(int iz=0; iz<dist.n_zgrid; iz++)
	       for(int ip=0; ip<dist.n_pgrid; ip++)
	          tmp.d3[ix][iy][iz].s[ip]=dist.d3[ix][iy][iz].s[ip]
                     / dist2.d3[ix][iy][iz].s[ip];
      return tmp;
   }

   // to avoid compiler warning 'no return statement' //AWS20050627
   Distribution tmp;                                  //AWS20050627
   return tmp;                                        //AWS20050627


}  
