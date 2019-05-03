//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * propel.cc *                                  galprop package * 10/18/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
// for method see notebook of 970300

#include <iostream>
#include <sstream>

#include"galprop_classes.h"
#include"galprop_internal.h"

using namespace std;

#include <ErrorLogger.h>

int propel_arrays_initialize=1;

Distribution alpha1_r, alpha1_x, alpha1_y, alpha1_z, alpha1_p;
Distribution alpha2_r, alpha2_x, alpha2_y, alpha2_z, alpha2_p;
Distribution alpha3_r, alpha3_x, alpha3_y, alpha3_z, alpha3_p;
Distribution Nr1_, Nx1_, Ny1_, Nz1_, Np1_;
Distribution Nr2_, Nx2_, Ny2_, Nz2_, Np2_;
Distribution Nr3_, Nx3_, Ny3_, Nz3_, Np3_;
Distribution total_source_function;
Distribution Work;

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Galprop::propel(Particle& particle) {

  INFO("Entry");

  int start=1;
  double end_timestep_sec = galdef.end_timestep*year2sec;
  double dt = galdef.start_timestep*year2sec;
  double factor = dt*pow(kpc2cm, -2.);  
  double t = 0., tcum =0.;  // cumulative time IMOS20061030
  //particle.cr_density=0.0; AWS20010121
  

  propel_diagnostics(); // AWS20110113 initialization for new species

  ostringstream buf;

  if (0 == particle.primary_source_function.max() && 
      0 == particle.secondary_source_function.max()) {

    buf.str("");
    buf << "Zero primary and secondary source functions. No propagation will be done";
    WARNING(buf.str());
    INFO("Exit");
    return 0;
 
  }

  if (2 == particle.n_spatial_dimensions) { // ==== 2D ====
   
    int ir,iz,ip;

    //Gulli20070810
    //      if(propel_arrays_initialize==1)                // initialization
    //      {
    propel_arrays_initialize=0;

    INFO("Generating alpha for 2D");

    alpha1_r.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    alpha1_z.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    alpha1_p.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    alpha2_r.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    alpha2_z.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    alpha2_p.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    alpha3_r.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    alpha3_z.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    alpha3_p.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    
    Nr1_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Nz1_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Np1_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Nr2_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Nz2_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Np2_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Nr3_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Nz3_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Np3_.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    
    total_source_function.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    Work.init(particle.n_rgrid, particle.n_zgrid, particle.n_pgrid);
    
    //      }     //propel_arrays_initialize==1
    
    // replaced "float" with "double"  IMOS20061214
    /*  Moved to beginning of for loop for parallization   Gulli20070821
	double *Nr0=new double[particle.n_rgrid];
	double *Nz0=new double[particle.n_zgrid];
	double *Np0=new double[particle.n_pgrid];
	double *Nr1=new double[particle.n_rgrid];
	double *Nz1=new double[particle.n_zgrid];
	double *Np1=new double[particle.n_pgrid];
	double *Nr2=new double[particle.n_rgrid];
	double *Nz2=new double[particle.n_zgrid];
	double *Np2=new double[particle.n_pgrid];
	double *Nr3=new double[particle.n_rgrid];
	double *Nz3=new double[particle.n_zgrid];
	double *Np3=new double[particle.n_pgrid];
	
	double *Rr =new double[particle.n_rgrid];
	double *Rz =new double[particle.n_zgrid];
	double *Rp =new double[particle.n_pgrid];*/
    
    alpha1_r = alpha1_z = alpha1_p = 0.;
    alpha2_r = alpha2_z = alpha2_p = 0.;
    alpha3_r = alpha3_z = alpha3_p = 0.;
    
    // DIFFUSION term: method according to ApJ 509, 212
    
    //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(ir,iz,ip)
    for (ir = 0; ir < particle.n_rgrid; ++ir) {
	  
      for (iz = 0; iz < particle.n_zgrid; ++iz) {

	for (ip = 0; ip < particle.n_pgrid; ++ip) {
		  
// test of electron propagation vs analytical calculations IMOS20061030	       
	  alpha1_z.d2[ir][iz].s[ip] = particle.Dxx.d2[ir][iz].s[ip]*pow(particle.dz,-2.);
	  alpha2_r.d2[ir][iz].s[ip] = particle.Dxx.d2[ir][iz].s[ip]*pow(particle.dr,-2.);
		  
	  //*** to compare with v50...
	  //		  if(ir==0){alpha3_r.d2[ir][iz].s[ip] =alpha2_r.d2[ir][iz].s[ip] *4.; alpha2_r.d2[ir][iz].s[ip] *=2.;}
	  //**********************/
	  if (0 == ir) // use: Dxx/(R[i-1/2] dR)*{Ri(U[i+1]-Ui)-R[i-1](Ui-U[i-1])}; R[i-1]=0
	    alpha3_r.d2[ir][iz].s[ip] = alpha2_r.d2[ir][iz].s[ip] *2.;

	  else { // use: Dxx/(Ri dR)*{R[i+1/2](U[i+1]-Ui)-R[i-1/2](Ui-U[i-1])}
	    
	    alpha1_r.d2[ir][iz].s[ip] = alpha2_r.d2[ir][iz].s[ip]*(1. - particle.dr/2./particle.r[ir]);
	    alpha3_r.d2[ir][iz].s[ip] = alpha2_r.d2[ir][iz].s[ip]*(1. + particle.dr/2./particle.r[ir]);
	
	  }

	  /* **************** v50 and earlier
	     if(ir==0)
	     {
	     alpha1_r.d2[ir][iz].s[ip] = 0.;
	     alpha2_r.d2[ir][iz].s[ip] 
	     = particle.Dxx.d2[ir][iz].s[ip] *4.0*pow(particle.dr,-2);
	     alpha3_r.d2[ir][iz].s[ip] = alpha2_r.d2[ir][iz].s[ip];
	     }
	     else
	     {
	     alpha1_r.d2[ir][iz].s[ip] = particle.Dxx.d2[ir][iz].s[ip]
	     *(1.0-particle.dr/2.0/particle.r[ir]) *pow(particle.dr,-2); 
	     alpha2_r.d2[ir][iz].s[ip] = particle.Dxx.d2[ir][iz].s[ip]
	     *2.0                                  *pow(particle.dr,-2);
	     alpha3_r.d2[ir][iz].s[ip] = particle.Dxx.d2[ir][iz].s[ip]
	     *(1.0+particle.dr/2.0/particle.r[ir]) *pow(particle.dr,-2); 
	     }
	  */
	  // method from propel.f90
	  /*	    
		    double r_1=particle.r[ir]-particle.dr*0.5;
		    if(r_1<0.0) r_1=0.0;
		    double r_2=particle.r[ir]+particle.dr*0.5;
		    double Const=particle.Dxx.d2[ir][iz].s[ip]*2.0/particle.dr
		    /(pow(r_2,2)-pow(r_1,2));
		    alpha1_r.d2[ir][iz].s[ip]= Const * r_1;
		    alpha2_r.d2[ir][iz].s[ip]= Const *(r_1+r_2);
		    alpha3_r.d2[ir][iz].s[ip]= Const * r_2;
	  */
	}   //ip
    
      }   //iz
    
    }   //ir
      
    alpha2_z = alpha1_z;   
    alpha2_z *= 2.;
    alpha3_z = alpha1_z;                                
      
    alpha2_r *= 2.; // IMOS20061030

    alpha1_r *= factor;   
    alpha1_z *= factor;
    alpha2_r *= factor;   
    alpha2_z *= factor;   
    alpha3_r *= factor;   
    alpha3_z *= factor;

    // CONVECTION term: method according to ApJ 509, 212 with corrections  IMOS20010307
    // convection introduced in v38 using imos code with modifications AWS20010325
    
    if (galdef.convection) {
			//Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ip,ir)
      for (iz = 0; iz < particle.n_zgrid; ++iz) { // numerator = abs convection velocity in cm s^-1
	
	double az1 = (particle.z[iz] > 0.) ? (galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz-1]))*1.e5/particle.dz*dt/kpc2cm: 0.;

	double az2 = (galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz]))*1.e5/particle.dz*dt/kpc2cm;
	
	double az3 = (particle.z[iz] < 0.) ? (galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz+1]))*1.e5/particle.dz*dt/kpc2cm: 0.;
	
	if (fabs(particle.z[iz]) < 0.5*particle.dz) 
	  az1 = az2 = az3 = 0.;
	
	for (ip = 1; ip < particle.n_pgrid-1; ++ip) {

	  // f90 scheme
	  double ap2 = particle.p[ip]/3.*galdef.dvdz_conv/(particle.p[ip+1]-particle.p[ip])*1.e5/kpc2cm;
	  
	  double ap3 = particle.p[ip+1]/3.*galdef.dvdz_conv/(particle.p[ip+1]-particle.p[ip])*1.e5/kpc2cm;
	  // new scheme
	  //               double ap2 = galdef.dvdz_conv/3.*(particle.p[ip]/(particle.p[ip+1]-particle.p[ip]) -1.)*1.e5/kpc2cm;
	  //               double ap3 = galdef.dvdz_conv/3.* particle.p[ip]/(particle.p[ip+1]-particle.p[ip])     *1.e5/kpc2cm;
	  // symmetrical scheme
	  //               double ap1 =particle.p[ip-1]/3. *galdef.dvdz_conv/(particle.p[ip+1]-particle.p[ip-1]) *1.e5/kpc2cm;
	  //               double ap3 =particle.p[ip+1]/3. *galdef.dvdz_conv/(particle.p[ip+1]-particle.p[ip-1]) *1.e5/kpc2cm;
	  
	  for (ir = 0; ir < particle.n_rgrid; ++ir) {

	    alpha1_z.d2[ir][iz].s[ip] += az1;
	    
	    alpha2_z.d2[ir][iz].s[ip] += az2;   
	    alpha2_p.d2[ir][iz].s[ip] += ap2;
            
	    alpha3_z.d2[ir][iz].s[ip] += az3;   
	    alpha3_p.d2[ir][iz].s[ip] += ap3;
	  }   //ir
            
	}   //ip
	
	for (ir = 0; ir < particle.n_rgrid; ++ir) {

	  alpha1_z.d2[ir][iz].s[particle.n_pgrid-1] += az1;
	  alpha2_z.d2[ir][iz].s[particle.n_pgrid-1] += az2;
	  alpha3_z.d2[ir][iz].s[particle.n_pgrid-1] += az3;
	  
	}   //ir
        
      }   //iz
      
    }

    // DIFFUSIVE REACCELERATION

    if (galdef.diff_reacc > 0) { //IMOS20060330
      //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ip,ir)
      for (ir = 0; ir < particle.n_rgrid; ++ir) {

	for (iz = 0; iz < particle.n_zgrid; ++iz) {

	  for (ip = 1; ip < particle.n_pgrid-1; ++ip) {

	    // alternative scheme #2, most detailed
	    alpha1_p.d2[ir][iz].s[ip] += 
	      (-(particle.Dpp.d2[ir][iz].s[ip] - particle.Dpp.d2[ir][iz].s[ip-1])
	       /(particle.p[ip] - particle.p[ip-1])
	       + 2.*particle.Dpp.d2[ir][iz].s[ip]
	       /(particle.p[ip+1] - particle.p[ip-1])
	       + 2.*particle.Dpp.d2[ir][iz].s[ip-1]/particle.p[ip-1]) 
	      /(particle.p[ip]-particle.p[ip-1]);

	    alpha2_p.d2[ir][iz].s[ip] +=
	      -(particle.Dpp.d2[ir][iz].s[ip]-particle.Dpp.d2[ir][iz].s[ip-1])
	      /pow(particle.p[ip] - particle.p[ip-1], 2.)
	      + 2.*particle.Dpp.d2[ir][iz].s[ip]
	      /(particle.p[ip+1] - particle.p[ip-1])
	      *(1./(particle.p[ip+1] - particle.p[ip])
		+ 1./(particle.p[ip] - particle.p[ip-1]))
	      + 2.*particle.Dpp.d2[ir][iz].s[ip]
	      /(particle.p[ip] - particle.p[ip-1])
	      /particle.p[ip];
 
	    alpha3_p.d2[ir][iz].s[ip] +=
	      2.*particle.Dpp.d2[ir][iz].s[ip]
	      /(particle.p[ip+1] - particle.p[ip-1])
	      /(particle.p[ip+1] - particle.p[ip]);

	    // old Andy's method
	    //                  double p1=(particle.p[ip]+particle.p[ip-1])/2.;
	    //                  double p2=(particle.p[ip]+particle.p[ip+1])/2.;
	    //                  double Dpp1 = particle.Dpp.d2 [ir] [iz].s[ip-1] *pow(p1/particle.p[ip-1],2.-galdef.D_g_1);
	    //                  double Dpp2 = particle.Dpp.d2 [ir] [iz].s[ip  ] *pow(p2/particle.p[ip  ],2.-galdef.D_g_1);
	    //                  alpha1_p.d2[ir] [iz].s[ip] +=    
	    //                     1./(p2-p1)* p1*p1*Dpp1/pow(particle.p[ip-1],2) /(particle.p[ip  ]-particle.p[ip-1]);
	    //                  alpha2_p.d2[ir] [iz].s[ip] +=    
	    //                     1./(p2-p1)*(p1*p1*Dpp1/pow(particle.p[ip  ],2) /(particle.p[ip  ]-particle.p[ip-1])
	    //                                +p2*p2*Dpp2/pow(particle.p[ip  ],2) /(particle.p[ip+1]-particle.p[ip  ]));
	    //                  alpha3_p.d2[ir] [iz].s[ip] +=    
	    //                     1./(p2-p1)* p2*p2*Dpp2/pow(particle.p[ip+1],2) /(particle.p[ip+1]-particle.p[ip  ]);
	  }   //ip
	
	}   //iz
      
      }   //ir
    
    }   // diffusive reacceleration

    // MOMENTUM LOSSES  IMOS20010307 minor change to make as in ApJ 509, 212
    // AWS20010622 correction to alpha2, consistent with ApJ 509, 212 Table 3 
    // AWS20010622 but note signs of alpha2, alpha3 incorrect in ApJ 509, 212 Table 3
    // AWS20010622 code is correct since dpdt = momentum loss rate is code
    
    if (galdef.momentum_losses) {
      //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ip,ir)
      for (ir = 0; ir < particle.n_rgrid; ++ir) {

	for (iz = 0; iz < particle.n_zgrid; ++iz) {

	  for (ip = 1; ip < particle.n_pgrid-1; ++ip) {

	    alpha2_p.d2[ir][iz].s[ip] += particle.dpdt.d2[ir][iz].s[ip]/(particle.p[ip+1] - particle.p[ip]);
	    
	    alpha3_p.d2[ir][iz].s[ip] += particle.dpdt.d2[ir][iz].s[ip+1]/(particle.p[ip+1] - particle.p[ip]);
              
	  }   //ip
          
	}   //iz
        
      }   //ir 
      
    }   //momentum losses

    // SET VALUES AT EXTREMES OF P GRID
    
    //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ir)
    for (ir = 0; ir < particle.n_rgrid; ++ir) {

      for (iz = 0; iz < particle.n_zgrid; ++iz) {

	alpha1_p.d2[ir][iz].s[0] = alpha1_p.d2[ir][iz].s[1];
	alpha2_p.d2[ir][iz].s[0] = alpha2_p.d2[ir][iz].s[1];
	alpha3_p.d2[ir][iz].s[0] = alpha3_p.d2[ir][iz].s[1];
	
	alpha1_p.d2[ir][iz].s[particle.n_pgrid-1] = alpha1_p.d2[ir][iz].s[particle.n_pgrid-2];
	alpha2_p.d2[ir][iz].s[particle.n_pgrid-1] = alpha2_p.d2[ir][iz].s[particle.n_pgrid-2];
	alpha3_p.d2[ir][iz].s[particle.n_pgrid-1] = alpha3_p.d2[ir][iz].s[particle.n_pgrid-2];
    
      }//iz
    
    }//ir 

    alpha1_p *= dt;
    alpha2_p *= dt;
    alpha3_p *= dt;
    
    double f_use = galdef.prop_r + galdef.prop_z + galdef.prop_p;

    // FRAGMENTATION
    
    if (galdef.fragmentation) {
       
      Work = particle.fragment;         // use Work to avoid memory leak
      Work*= (dt/f_use);
      alpha2_r += Work;                  // fragment used f_use times
      alpha2_z += Work;                         
      alpha2_p += Work;                          
      // alpha2_r+=particle.fragment*(dt/f_use);//  fragment used f_use times
      // alpha2_z+=particle.fragment*(dt/f_use);
      // alpha2_p+=particle.fragment*(dt/f_use);
    }

    // DECAY
    
    if (particle.t_half != 0 && galdef.radioactive_decay) {

      ostringstream rBuf;
      rBuf << "Radioactive decay of " << particle.name << " with half life " << particle.t_half;
      INFO(rBuf.str());
 
      Work = particle.decay;            // use Work to avoid memory leak
      Work *= (dt/f_use);
      alpha2_r += Work;                  // decay used f_use times
      alpha2_z += Work;                         
      alpha2_p += Work;
      //alpha2_r+=particle.decay   *(dt/f_use);//  decay used f_use times
      //alpha2_z+=particle.decay   *(dt/f_use);
      //alpha2_p+=particle.decay   *(dt/f_use);
    }

    buf.str("");
    buf<<"alpha1_r:"<< alpha1_r.d2[0][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha1_r.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha2_r.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha3_r.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha1_z.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha2_z.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha3_z.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha1_p.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha2_p.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_r:"<< alpha3_p.d2[1][1].s[1];
    DEBUGLOG(buf.str());
    
    if (galdef.verbose==10) {
       
      cout<<"alpha1_r:"<<endl;  alpha1_r.print();
      cout<<"alpha2_r:"<<endl;  alpha2_r.print();
      cout<<"alpha3_r:"<<endl;  alpha3_r.print();
      cout<<"alpha1_z:"<<endl;  alpha1_z.print();
      cout<<"alpha2_z:"<<endl;  alpha2_z.print();
      cout<<"alpha3_z:"<<endl;  alpha3_z.print();
      cout<<"alpha1_p:"<<endl;  alpha1_p.print();
      cout<<"alpha2_p:"<<endl;  alpha2_p.print();
      cout<<"alpha3_p:"<<endl;  alpha3_p.print();
    
    }

    ostringstream tBuf;
    tBuf << "galdef_ID " << galdef.galdef_ID << " start timestep dt = " << dt << " " << particle.name << endl;
    DEBUGLOG(tBuf.str());
    
    while (1 == start || dt > end_timestep_sec/galdef.timestep_factor*0.9999) {
      
      if (0 == start) {

	dt *= galdef.timestep_factor;
	tBuf.str("");
	tBuf<<" galdef_ID "<<galdef.galdef_ID<<" new timestep dt="<<dt<<"  "
	   <<particle.name;
	DEBUGLOG(tBuf.str());

	alpha1_r *= galdef.timestep_factor;
	alpha2_r *= galdef.timestep_factor;
	alpha3_r *= galdef.timestep_factor;
	alpha1_z *= galdef.timestep_factor;
	alpha2_z *= galdef.timestep_factor;
	alpha3_z *= galdef.timestep_factor;
	alpha1_p *= galdef.timestep_factor;
	alpha2_p *= galdef.timestep_factor;
	alpha3_p *= galdef.timestep_factor;
      
      }
      
      start = 0;

      double ADI_factor = 1;
      //  IF (ADI ==1) ADI_factor = 1.0/3.0
      
      /* leads to memory leak
	 Nr1_ =             -  0.5* alpha1_r * ADI_factor;
	 Nr2_ =   1.0 +        0.5* alpha2_r * ADI_factor;              
	 Nr3_ =             -  0.5* alpha3_r * ADI_factor;
	 
	 
	 Nz1_ =               -0.5* alpha1_z  * ADI_factor;
	 Nz2_ =   1.0 +        0.5* alpha2_z  * ADI_factor; 
	 Nz3_ =               -0.5* alpha3_z  * ADI_factor;
	 
	 Np1_ =              - 0.5* alpha1_p * ADI_factor;
	 Np2_ =   1.0 +        0.5* alpha2_p * ADI_factor;  
	 Np3_ =              - 0.5* alpha3_p * ADI_factor;
      */
      
      Nr1_ = alpha1_r;   Nr1_*= -0.5;
      Nr2_ = alpha2_r;   Nr2_*=  0.5;   Nr2_ += 1.0;   
      Nr3_ = alpha3_r;   Nr3_*= -0.5;
      Nz1_ = alpha1_z;   Nz1_*= -0.5;
      Nz2_ = alpha2_z;   Nz2_*=  0.5;   Nz2_ += 1.0;   
      Nz3_ = alpha3_z;   Nz3_*= -0.5;
      Np1_ = alpha1_p;   Np1_*= -0.5;
      Np2_ = alpha2_p;   Np2_*=  0.5;   Np2_ += 1.0;   
      Np3_ = alpha3_p;   Np3_*= -0.5;
      
      buf.str("");
      buf<<"Nr1_    ;"<< Nr1_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      buf.str("");
      buf<<"Nr2_    :"<< Nr2_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      buf.str("");
      buf<<"Nr3_    :"<< Nr3_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      buf.str("");
      buf<<"Nz1_    :"<< Nz1_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      buf.str("");
      buf<<"Nz2_    :"<< Nz2_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      buf.str("");
      buf<<"Nz3_    :"<< Nz3_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      buf.str("");
      buf<<"Np1_    :"<< Np1_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      buf.str("");
      buf<<"Np2_    :"<< Np2_.d2[1][1].s[1];
      DEBUGLOG(buf.str());
      buf.str("");
      buf<<"Np3_    :"<< Np3_.d2[1][1].s[1];
      DEBUGLOG(buf.str());

      // total_source_function = particle.primary_source_function + particle.secondary_source_function;
      
      total_source_function = particle.primary_source_function;
      total_source_function += particle.secondary_source_function;

      if (0 == total_source_function.max()) 
	return 0;

      // PROPAGATION iterations
      
      for (int irept = 1; irept <= galdef.timestep_repeat; ++irept) {
	
	// Z propagation
	//Gulli20070821
#pragma omp parallel default(shared) private(iz,ir,ip)
	if (1 == galdef.prop_z) {

	  //double* Nz0 = new double[particle.n_zgrid];    //Gulli20070821
	  //double* Nz1 = new double[particle.n_zgrid];    //Gulli20070821
	  //double* Nz2 = new double[particle.n_zgrid];    //Gulli20070821
	  //double* Nz3 = new double[particle.n_zgrid];    //Gulli20070821
	  //double* Rz = new double[particle.n_zgrid];    //Gulli20070821

	  valarray<double> Nz0(0., particle.n_zgrid);
	  valarray<double> Nz1(0., particle.n_zgrid);
	  valarray<double> Nz2(0., particle.n_zgrid);
	  valarray<double> Nz3(0., particle.n_zgrid);
	  valarray<double> Rz(0., particle.n_zgrid);

	  //Gulli20070821
#pragma omp for schedule(dynamic) 
	  for (ir = 0; ir < particle.n_rgrid; ++ir) {
	    
	    for (ip = 0; ip < particle.n_pgrid; ++ip) {
	      
	      for (iz = 0; iz < particle.n_zgrid; ++iz) {
                  
		Nz1[iz] = Nz1_.d2[ir][iz].s[ip];
		Nz2[iz] = Nz2_.d2[ir][iz].s[ip];
		Nz3[iz] = Nz3_.d2[ir][iz].s[ip];
		
		Nz0[iz] = total_source_function.d2[ir][iz].s[ip]*dt/f_use;
		Nz0[iz] += (1. - 0.5*alpha2_z.d2[ir][iz].s[ip])*particle.cr_density.d2[ir][iz].s[ip];
             
		if (iz > 0) 
		  Nz0[iz] += 0.5*alpha1_z.d2[ir][iz].s[ip] 
		    *particle.cr_density.d2[ir][iz-1].s[ip];
                    
		if (iz < particle.n_zgrid-1) 
		  Nz0[iz] += 0.5*alpha3_z.d2[ir][iz].s[ip] 
		    *particle.cr_density.d2[ir][iz+1].s[ip];
                    
	      } 

	      tridag(&Nz1[0], &Nz2[0], &Nz3[0], &Nz0[0], &Rz[0], particle.n_zgrid);  

	      for (iz = 0; iz < particle.n_zgrid; ++iz)
		particle.cr_density.d2[ir][iz].s[ip] = (Rz[iz] < 0.) ? 0. : Rz[iz]; //IMOS20030217
	      iz = 0;                                     // boundary condition
	      particle.cr_density.d2[ir][iz].s[ip] = 0;

	      iz = particle.n_zgrid-1;                    // boundary condition
	      particle.cr_density.d2[ir][iz].s[ip] = 0;
                  
	    }   //ip
            
	  }   //ir
      
	  //delete[] Nz0;    //Gulli20070821
	  //delete[] Nz1;    //Gulli20070821
	  //delete[] Nz2;    //Gulli20070821
	  //delete[] Nz3;    //Gulli20070821
	  //delete[] Rz ;    //Gulli20070821
	
	}   //prop_z

	// P propagation

	//Gulli20070821
#pragma omp parallel default(shared) private(iz,ir,ip)
	if (1 == galdef.prop_p) { 
	  
	  //double* Np0 = new double[particle.n_pgrid];    //Gulli20070821
	  //double* Np1 = new double[particle.n_pgrid];    //Gulli20070821
	  //double* Np2 = new double[particle.n_pgrid];    //Gulli20070821
	  //double* Np3 = new double[particle.n_pgrid];    //Gulli20070821
	  //double* Rp = new double[particle.n_pgrid];    //Gulli20070821
	  
	  valarray<double> Np0(0., particle.n_pgrid);
	  valarray<double> Np1(0., particle.n_pgrid);
	  valarray<double> Np2(0., particle.n_pgrid);
	  valarray<double> Np3(0., particle.n_pgrid);
	  valarray<double> Rp(0., particle.n_pgrid);

	  //Gulli20070821
#pragma omp for schedule(dynamic) 
	  for (ir = 0; ir < particle.n_rgrid; ++ir) {
	    
	    for (iz = 0; iz < particle.n_zgrid; ++iz) {
	      
	      for (ip = 0; ip < particle.n_pgrid; ++ip) {

		Np1[ip] = Np1_.d2[ir][iz].s[ip];
		Np2[ip] = Np2_.d2[ir][iz].s[ip];
		Np3[ip] = Np3_.d2[ir][iz].s[ip];
		
		Np0[ip] = total_source_function.d2[ir][iz].s[ip]*dt/f_use;
                        
		Np0[ip] += (1. - 0.5*alpha2_p.d2[ir][iz].s[ip])
		  *particle.cr_density.d2[ir] [iz].s[ip];
                        
		if (ip > 0) 
		  Np0[ip] += 0.5*alpha1_p.d2[ir][iz].s[ip] 
		    *particle.cr_density.d2[ir][iz].s[ip-1];
                    
		if (ip < particle.n_pgrid-1) 
		  Np0[ip] += 0.5*alpha3_p.d2[ir][iz].s[ip] 
		    *particle.cr_density.d2[ir][iz].s[ip+1];
                     
	      } 

	      tridag(&Np1[0], &Np2[0], &Np3[0], &Np0[0], &Rp[0], particle.n_pgrid);  

	      for (ip = 0; ip < particle.n_pgrid; ++ip)
		particle.cr_density.d2[ir][iz].s[ip] = (Rp[ip] < 0.) ? 0.: Rp[ip]; //IMOS20030217
	      
	      // ir=particle.n_rgrid-1;  particle.cr_density.d2[ir][iz].s[ip]=0.0; // boundary condition
	      // if(irept%  1==0) cout<<"\nparticle.cr_density after iteration "<<irept<<endl;
	      // particle.cr_density.print();
	    }   //iz
               
	  }   //ir
	  
	  //delete[] Np0;    //Gulli20070821
	  //delete[] Np1;    //Gulli20070821
	  //delete[] Np2;    //Gulli20070821
	  //delete[] Np3;    //Gulli20070821
	  //delete[] Rp ;    //Gulli20070821
	
	}   //prop_p

	// R propagation
				
	//Gulli20070821
#pragma omp parallel default(shared) private(iz,ir,ip)
	if (1 == galdef.prop_r) {

	  //double* Nr0 = new double[particle.n_rgrid];    //Gulli20070821
	  //double* Nr1 = new double[particle.n_rgrid];    //Gulli20070821
	  //double* Nr2 = new double[particle.n_rgrid];    //Gulli20070821
	  //double* Nr3 = new double[particle.n_rgrid];    //Gulli20070821
	  //double* Rr = new double[particle.n_rgrid];    //Gulli20070821

	  valarray<double> Nr0(0., particle.n_rgrid);
	  valarray<double> Nr1(0., particle.n_rgrid);
	  valarray<double> Nr2(0., particle.n_rgrid);
	  valarray<double> Nr3(0., particle.n_rgrid);
	  valarray<double> Rr(0.,  particle.n_rgrid);

 //Gulli20070821
#pragma omp for schedule(dynamic) 
	  for (iz = 0; iz < particle.n_zgrid; ++iz) {
	    
	    for (ip = 0; ip < particle.n_pgrid; ++ip) {

	      for (ir = 0; ir < particle.n_rgrid; ++ir) {

		Nr1[ir] = Nr1_.d2[ir][iz].s[ip];
		Nr2[ir] = Nr2_.d2[ir][iz].s[ip];
		Nr3[ir] = Nr3_.d2[ir][iz].s[ip];
		
		Nr0[ir] = total_source_function.d2[ir][iz].s[ip]*dt/f_use;
		
		Nr0[ir] += (1. - 0.5*alpha2_r.d2[ir][iz].s[ip])
		  *particle.cr_density.d2[ir][iz].s[ip];
                        
		if (ir > 0) 
		  Nr0[ir] += 0.5*alpha1_r.d2[ir][iz].s[ip] 
		    *particle.cr_density.d2[ir-1][iz].s[ip];
                    
		if (ir < particle.n_rgrid-1) 
		  Nr0[ir] += 0.5*alpha3_r.d2[ir][iz].s[ip] 
		    *particle.cr_density.d2[ir+1][iz].s[ip];
                    
	      }

	      tridag(&Nr1[0], &Nr2[0], &Nr3[0], &Nr0[0], &Rr[0], particle.n_rgrid);  

	      for (ir = 0; ir < particle.n_rgrid; ++ir)
		particle.cr_density.d2[ir][iz].s[ip] = (Rr[ir] < 0.) ? 0.: Rr[ir]; //IMOS20030217
  
	      ir = particle.n_rgrid-1;                  // boundary condition
	      particle.cr_density.d2[ir][iz].s[ip] = 0;
//if(irept%500==0)                                                                              
//cout<<"\nparticle.cr_density:"<<endl;particle.cr_density.print();
	    }   //ip
	  
	  }   //iz
      
	  //delete[] Nr0;    //Gulli20070821
	  //delete[] Nr1;    //Gulli20070821
	  //delete[] Nr2;    //Gulli20070821
	  //delete[] Nr3;    //Gulli20070821
	  //delete[] Rr ;    //Gulli20070821
            
	}   //prop_r

	if (0 == irept%galdef.timestep_print) {

	  ostringstream repBuf;
	  repBuf << " particle.cr_density after iteration "<<irept<<" with time step " <<dt/year2sec<<" yr";
	  INFO(repBuf.str());
	  particle.cr_density.print();   cout<<endl;
            
	}
        
//	cout<<"\b\b\b\b\b"<<irept;   // \b=backspace
	buf.str("");
	buf<<"maximum="<<particle.cr_density.max();
	buf<<"    propel iteration "<<irept<<endl;     

	if(irept%galdef.timestep_diagnostics==0)
	  propel_diagnostics (particle, alpha1_r, alpha1_z, alpha1_p, 
			      alpha2_r, alpha2_z, alpha2_p,
			      alpha3_r, alpha3_z, alpha3_p,
			      total_source_function, dt);
      }   //irept
    
    }   //dt>end_timestep
			//Gulli20070810
    
    alpha1_r.delete_array();
    alpha1_z.delete_array();
    alpha1_p.delete_array();
    alpha2_r.delete_array();
    alpha2_z.delete_array();
    alpha2_p.delete_array();
    alpha3_r.delete_array();
    alpha3_z.delete_array();
    alpha3_p.delete_array();
    
    Nr1_.delete_array();
    Nz1_.delete_array();
    Np1_.delete_array();
    Nr2_.delete_array();
    Nz2_.delete_array();
    Np2_.delete_array();
    Nr3_.delete_array();
    Nz3_.delete_array();
    Np3_.delete_array();
    
    total_source_function.delete_array();
    Work.delete_array();
    /* Moved to for loops to get the right scope   Gulli20070821
       delete[] Nr0;
       delete[] Nz0;
       delete[] Np0;
       delete[] Nr1;
       delete[] Nz1;
       delete[] Np1;
       delete[] Nr2;
       delete[] Nz2;
       delete[] Np2;
       delete[] Nr3;
       delete[] Nz3;
       delete[] Np3;
       
       delete[] Rr ;
       delete[] Rz ;
       delete[] Rp ;*/
  }   //2D case
  
  //////////////////////////////////////////////////////////////////////////

  if (3 == particle.n_spatial_dimensions) {             // ==== 3D ====
   
    int ix,iy,iz,ip;

    //Gulli20070810
    //if(propel_arrays_initialize==1)
    //{
    propel_arrays_initialize = 0;
    
    ostringstream aBuf;
    aBuf << "Generating alpha for 3D";
    INFO(aBuf.str());


    alpha1_x.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha1_y.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha1_z.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha1_p.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha2_x.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha2_y.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha2_z.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha2_p.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha3_x.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha3_y.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha3_z.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    alpha3_p.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    
    Nx1_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Ny1_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Nz1_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Np1_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Nx2_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Ny2_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Nz2_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Np2_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Nx3_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Ny3_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Nz3_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Np3_.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    
    total_source_function.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    Work.init(particle.n_xgrid,particle.n_ygrid,particle.n_zgrid,particle.n_pgrid);
    //}   //propel_arrays_initialize==1

    // replaced "float" with "double"  IMOS20061214
    /* Moved to for loops for parallization     Gulli20070821
       double *Nx0 = new double[particle.n_xgrid];
       double *Ny0 = new double[particle.n_ygrid];
       double *Nz0 = new double[particle.n_zgrid];
       double *Np0 = new double[particle.n_pgrid];
       double *Nx1 = new double[particle.n_xgrid];
       double *Ny1 = new double[particle.n_ygrid];
       double *Nz1 = new double[particle.n_zgrid];
       double *Np1 = new double[particle.n_pgrid];
       double *Nx2 = new double[particle.n_xgrid];
       double *Ny2 = new double[particle.n_ygrid];
       double *Nz2 = new double[particle.n_zgrid];
       double *Np2 = new double[particle.n_pgrid];
       double *Nx3 = new double[particle.n_xgrid];
       double *Ny3 = new double[particle.n_ygrid];
       double *Nz3 = new double[particle.n_zgrid];
       double *Np3 = new double[particle.n_pgrid];
       
       double *Rx = new double[particle.n_xgrid];
       double *Ry = new double[particle.n_ygrid];
       double *Rz = new double[particle.n_zgrid];
       double *Rp = new double[particle.n_pgrid];*/
    
    alpha1_x = alpha1_y = alpha1_z = alpha1_p = 0;
    alpha2_x = alpha2_y = alpha2_z = alpha2_p = 0;
    alpha3_x = alpha3_y = alpha3_z = alpha3_p = 0;
    
    // DIFFUSION term / Dxx = const in space 

    //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ix,iy,ip)
    for (ix = 0; ix < particle.n_xgrid; ++ix) {
      
      for (iy = 0; iy < particle.n_ygrid; ++iy) {
	
	for (iz = 0; iz < particle.n_zgrid; ++iz) {

	  for (ip = 0; ip < particle.n_pgrid; ++ip) {
	   int lowindex_x,lowindex_y,lowindex_z ;
	   int highindex_x,highindex_y,highindex_z ;
	   lowindex_x = ix-1;
	   lowindex_y = iy-1;
	   lowindex_z = iz-1;	
	   if (lowindex_x < 0)  lowindex_x = 0;
	   if (lowindex_y < 0)  lowindex_y = 0;
	   if (lowindex_z < 0)  lowindex_z = 0; 	  
	   highindex_x = ix+1; 
	   highindex_y = iy+1;
	   highindex_z = iz+1; 
	   if (highindex_x > particle.n_xgrid -1)  highindex_x = particle.n_xgrid -1;
	   if (highindex_y > particle.n_ygrid-1)  highindex_y = particle.n_ygrid-1;
	   if (highindex_z > particle.n_zgrid-1)  highindex_z = particle.n_zgrid-1;  	 
	    alpha1_x.d3[ix][iy][iz].s[ip]
	      = particle.Dxx.d3[ix][iy][iz].s[ip]*pow(particle.dx, -2.) - (particle.Dxx.d3[highindex_x][iy][iz].s[ip]-particle.Dxx.d3[lowindex_x][iy][iz].s[ip])/4.*pow(particle.dx, -2.);
	    
	    alpha1_y.d3[ix][iy][iz].s[ip]
	      = particle.Dxx.d3[ix][iy][iz].s[ip]*pow(particle.dy, -2.) - (particle.Dxx.d3[ix][highindex_y][iz].s[ip]-particle.Dxx.d3[ix][lowindex_y][iz].s[ip])/4.*pow(particle.dx, -2.);
              
	    alpha1_z.d3[ix][iy][iz].s[ip]
	      = particle.Dxx.d3[ix][iy][iz].s[ip]*pow(particle.dz, -2.)-(particle.Dxx.d3[ix][iy][highindex_z].s[ip]-particle.Dxx.d3[ix][iy][lowindex_z].s[ip])/4.*pow(particle.dx, -2.);
//   alpha 2  
	    alpha2_x.d3[ix][iy][iz].s[ip]
	      = particle.Dxx.d3[ix][iy][iz].s[ip]*pow(particle.dx, -2.)*2.; 
	    alpha2_y.d3[ix][iy][iz].s[ip]
	      = particle.Dxx.d3[ix][iy][iz].s[ip]*pow(particle.dy, -2.)*2.;
	    alpha2_z.d3[ix][iy][iz].s[ip]
	      = particle.Dxx.d3[ix][iy][iz].s[ip]*pow(particle.dz, -2.)*2.;
//   alpha 3
	    alpha3_x.d3[ix][iy][iz].s[ip] 
	      = particle.Dxx.d3[ix][iy][iz].s[ip]*pow(particle.dx, -2.) + (particle.Dxx.d3[highindex_x][iy][iz].s[ip]-particle.Dxx.d3[lowindex_x][iy][iz].s[ip])/4.*pow(particle.dx, -2.);
	    alpha3_y.d3[ix][iy][iz].s[ip]
	      = particle.Dxx.d3[ix][iy][iz].s[ip]*pow(particle.dy, -2.) + (particle.Dxx.d3[ix][highindex_y][iz].s[ip]-particle.Dxx.d3[ix][lowindex_y][iz].s[ip])/4.*pow(particle.dx, -2.); 
	    alpha3_z.d3[ix][iy][iz].s[ip]
	      = particle.Dxx.d3[ix][iy][iz].s[ip]*pow(particle.dz, -2.)+(particle.Dxx.d3[ix][iy][highindex_z].s[ip]-particle.Dxx.d3[ix][iy][lowindex_z].s[ip])/4.*pow(particle.dx, -2.);
              
	  }   //ip
          
	}   //iz
        
      }   //iy
      
    }   //ix

              
	 
//    alpha2_x = alpha1_x;  
//    alpha2_x *= 2.;
    
//    alpha3_x = alpha1_x;
    
//    alpha2_y = alpha1_y;  
//    alpha2_y *= 2.;
    
//    alpha3_y = alpha1_y;
    
//    alpha2_z = alpha1_z;  
//    alpha2_z *= 2.;
    
//    alpha3_z = alpha1_z;

    alpha1_x *= factor;   
    alpha1_y *= factor;   
    alpha1_z *= factor;
    
    alpha2_x *= factor;   
    alpha2_y *= factor;   
    alpha2_z *= factor; 
    
    alpha3_x *= factor;   
    alpha3_y *= factor;   
    alpha3_z *= factor;

// CONVECTION                                                      AWS20010330

    if (galdef.convection) {
      //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ix,iy,ip)
      for (iz = 0; iz < particle.n_zgrid; ++iz) { // numerator = abs convection velocity in cm s^-1
		
	double az1 = (particle.z[iz] > 0.) ? (galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz-1]))*1.e5/particle.dz*dt/kpc2cm: 0;  
	
	double az2 = (galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz]))*1.e5/particle.dz*dt/kpc2cm;    
        
	double az3 = (particle.z[iz] < 0.) ? (galdef.v0_conv + galdef.dvdz_conv*fabs(particle.z[iz+1])) *1.e5/particle.dz*dt/kpc2cm: 0;
	
	if (fabs(particle.z[iz]) < 0.5*particle.dz) 
	  az1 = az2 = az3 = 0.;

	for (ip = 0; ip < particle.n_pgrid-1; ++ip) {
	  
	  double ap2 = particle.p[ip]/3.*galdef.dvdz_conv/(particle.p[ip+1] - particle.p[ip])*1.e5/kpc2cm;
	  
	  double ap3 = particle.p[ip+1]/3.*galdef.dvdz_conv/(particle.p[ip+1] - particle.p[ip])*1.e5/kpc2cm;

	  for (ix = 0; ix < particle.n_xgrid; ++ix)
	    for (iy = 0; iy < particle.n_ygrid; ++iy) {

	      alpha1_z.d3[ix][iy][iz].s[ip] += az1;                                      
	      alpha2_z.d3[ix][iy][iz].s[ip] += az2;  
	      alpha2_p.d3[ix][iy][iz].s[ip] += ap2;
	
	      alpha3_z.d3[ix][iy][iz].s[ip] += az3;  
	      alpha3_p.d3[ix][iy][iz].s[ip] += ap3;
	      
	    }   //iy ix
            
	}   //ip
        
      }   //iz
      
    }

// DIFFUSIVE REACCELERATION   IMOS20020329

    if (galdef.diff_reacc>0) { //IMOS20060330
     
      //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ix,iy,ip)
      for (ix = 0; ix < particle.n_xgrid; ++ix) {
	
	for (iy = 0; iy < particle.n_ygrid; ++iy) {
	  
	  for (iz = 0; iz < particle.n_zgrid; ++iz) {

	    for (ip = 1; ip < particle.n_pgrid-1; ++ip) {
	      // alternative scheme #2, most detailed
	      alpha1_p.d3[ix][iy][iz].s[ip] += 
		(-(particle.Dpp.d3[ix][iy][iz].s[ip] - 
		   particle.Dpp.d3[ix][iy][iz].s[ip-1])
		 /(particle.p[ip] - particle.p[ip-1])
		 + 2.*particle.Dpp.d3[ix][iy][iz].s[ip]
		 /(particle.p[ip+1]-particle.p[ip-1])
		 + 2.*particle.Dpp.d3[ix][iy][iz].s[ip-1]
		 /particle.p[ip-1]) 
		/(particle.p[ip] - particle.p[ip-1]);

	      alpha2_p.d3[ix][iy][iz].s[ip] +=
		-(particle.Dpp.d3[ix][iy][iz].s[ip] - 
		  particle.Dpp.d3[ix][iy][iz].s[ip-1])
		/pow(particle.p[ip] - particle.p[ip-1], 2.)
		+ 2.*particle.Dpp.d3[ix][iy][iz].s[ip]
		/(particle.p[ip+1] - particle.p[ip-1])
		*(1./(particle.p[ip+1] - particle.p[ip])
		  + 1./(particle.p[ip] - particle.p[ip-1]))
		+ 2.*particle.Dpp.d3[ix][iy][iz].s[ip]
		/(particle.p[ip] - particle.p[ip-1])
		/particle.p[ip]; 
                  
	      alpha3_p.d3[ix][iy][iz].s[ip] +=
		2.*particle.Dpp.d3[ix][iy][iz].s[ip]
		/(particle.p[ip+1] - particle.p[ip-1])
		/(particle.p[ip+1] - particle.p[ip]);
	      // old Andy's method
	      /*
		double p1=(particle.p[ip]+particle.p[ip-1])/2.;
		double p2=(particle.p[ip]+particle.p[ip+1])/2.;
		double Dpp1 = particle.Dpp.d3[ix][iy][iz].s[ip-1] *pow(p1/particle.p[ip-1],2.-galdef.D_g_1);
		double Dpp2 = particle.Dpp.d3[ix][iy][iz].s[ip  ] *pow(p2/particle.p[ip  ],2.-galdef.D_g_1);
		alpha1_p.d3[ix][iy][iz].s[ip]+= 
		1./(p2-p1)* p1*p1*Dpp1/pow(particle.p[ip-1],2)/(particle.p[ip  ]-particle.p[ip-1]);
		alpha2_p.d3[ix][iy][iz].s[ip]+= 
		1./(p2-p1)*( p1*p1*Dpp1/pow(particle.p[ip  ],2)/(particle.p[ip  ]-particle.p[ip-1])
		+p2*p2*Dpp2/pow(particle.p[ip  ],2)/(particle.p[ip+1]-particle.p[ip  ]) );
		alpha3_p.d3[ix][iy][iz].s[ip]+= 
		1./(p2-p1)* p2*p2*Dpp2/pow(particle.p[ip+1],2)/(particle.p[ip+1]-particle.p[ip  ]);
	      */
	    
	    }   //ip
	  
	  }   //iz
          
	}   //iy
        
      }   //ix
      
    }   // diffusive reacceleration

    // MOMENTUM LOSSES
    
    if(galdef.momentum_losses) {
      //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ix,iy,ip)
      for (ix = 0; ix < particle.n_xgrid; ++ix) {

	for (iy = 0; iy < particle.n_ygrid; ++iy) {

	  for (iz = 0; iz < particle.n_zgrid; ++iz) {

	    for (ip = 1; ip < particle.n_pgrid-1; ++ip) {// bug fixed IMOS20020417
	      
	      alpha2_p.d3[ix][iy][iz].s[ip] += 
		particle.dpdt.d3[ix][iy][iz].s[ip]
		/(particle.p[ip+1] - particle.p[ip]);
                     
	      alpha3_p.d3[ix][iy][iz].s[ip] += 
		particle.dpdt.d3[ix][iy][iz].s[ip+1]
		/(particle.p[ip+1]-particle.p[ip]);
                
	    }   //ip
            
	  }   //iz
          
	}   //iy
        
      }   //ix 
      
    }   //momentum losses

    // fill in values at extremes of p grid

    //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ix,iy)
    for (ix = 0; ix < particle.n_xgrid; ++ix) {
      
      for (iy = 0; iy < particle.n_ygrid; ++iy) {

	for (iz = 0; iz < particle.n_zgrid; ++iz) {

	  alpha1_p.d3[ix][iy][iz].s[0] = alpha1_p.d3[ix][iy][iz].s[1];
	  alpha2_p.d3[ix][iy][iz].s[0] = alpha2_p.d3[ix][iy][iz].s[1];
	  alpha3_p.d3[ix][iy][iz].s[0] = alpha3_p.d3[ix][iy][iz].s[1];
	  
	  alpha1_p.d3[ix][iy][iz].s[particle.n_pgrid-1] = alpha1_p.d3[ix][iy][iz].s[particle.n_pgrid-2];
	  alpha2_p.d3[ix][iy][iz].s[particle.n_pgrid-1] = alpha2_p.d3[ix][iy][iz].s[particle.n_pgrid-2];
	  alpha3_p.d3[ix][iy][iz].s[particle.n_pgrid-1] = alpha3_p.d3[ix][iy][iz].s[particle.n_pgrid-2];
          
	}   //iz
      
      }   //iy
      
    }   //ix 

    alpha1_p *= dt;
    alpha2_p *= dt;
    alpha3_p *= dt;

    double f_use = galdef.prop_x + galdef.prop_y + galdef.prop_z + galdef.prop_p;

    // FRAGMENTATION

    if (galdef.fragmentation) {
       
      Work = particle.fragment;        // use Work to avoid memory leak
      Work *= (dt/f_use);
      alpha2_x += Work;                 //  fragment used f_use times
      alpha2_y += Work;  
      alpha2_z += Work;                         
      alpha2_p += Work;
      // alpha2_x+=particle.fragment*(dt/f_use);//  fragment used  f_use times
      // alpha2_y+=particle.fragment*(dt/f_use);
      // alpha2_z+=particle.fragment*(dt/f_use);
      // alpha2_p+=particle.fragment*(dt/f_use);
    
    }

    // RADIOACTIVE DECAY

    if (particle.t_half != 0.0 && galdef.radioactive_decay) {
         
      ostringstream rBuf;
      rBuf << "Radioactive decay of "<<particle.name<<" with half_life "
	   <<particle.t_half;
      INFO(rBuf.str());

      Work = particle.decay;           // use Work to avoid memory leak
      Work *= (dt/f_use);
      alpha2_x += Work;                 //  fragment used f_use times
      alpha2_y += Work; 
      alpha2_z += Work;                         
      alpha2_p += Work; 

      // alpha2_x+=particle.decay   *(dt/f_use);//  decay used  f_use times
      // alpha2_y+=particle.decay   *(dt/f_use);
      // alpha2_z+=particle.decay   *(dt/f_use);
      // alpha2_p+=particle.decay   *(dt/f_use);
    }
    
    buf.str("");
    buf<<"alpha1_x:"<< alpha1_x.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha2_x:"<< alpha2_x.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha3_x:"<< alpha3_x.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_y:"<< alpha1_y.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha2_y:"<< alpha2_y.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha3_y:"<< alpha3_y.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_z:"<< alpha1_z.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha2_z:"<< alpha2_z.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha3_z:"<< alpha3_z.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha1_p:"<< alpha1_p.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha2_p:"<< alpha2_p.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());
    buf.str("");
    buf<<"alpha3_p:"<< alpha3_p.d3[1][1][1].s[1];
    DEBUGLOG(buf.str());

    if (galdef.verbose == 10) {
       
      cout<<"alpha1_x:"<<endl;   alpha1_x.print();
      cout<<"alpha2_x:"<<endl;   alpha2_x.print();
      cout<<"alpha3_x:"<<endl;   alpha3_x.print();
      cout<<"alpha1_y:"<<endl;   alpha1_y.print();
      cout<<"alpha2_y:"<<endl;   alpha2_y.print();
      cout<<"alpha3_y:"<<endl;   alpha3_y.print();
      cout<<"alpha1_z:"<<endl;   alpha1_z.print();
      cout<<"alpha2_z:"<<endl;   alpha2_z.print();
      cout<<"alpha3_z:"<<endl;   alpha3_z.print();
      cout<<"alpha1_p:"<<endl;   alpha1_p.print();
      cout<<"alpha2_p:"<<endl;   alpha2_p.print();
      cout<<"alpha3_p:"<<endl;   alpha3_p.print();

    }
    
    ostringstream tBuf;

    tBuf <<"galdef_ID "<<galdef.galdef_ID<<" start timestep dt="<<dt<<"  "
         <<particle.name;
    DEBUGLOG(tBuf.str());

    int timestep_mode;
    int timestep_mode2done = 0;
    
    for (timestep_mode = 1; timestep_mode <= 2; ++timestep_mode) {
      
       buf.str("");
      buf<<"    timestep_mode="<<timestep_mode;
      INFO(buf.str());
      t = 0;
         
      while ((1 == start) || 
	     (dt > end_timestep_sec/galdef.timestep_factor*0.9999) || 
	     (2 == timestep_mode && 0 == timestep_mode2done)) {
         
	if (0 == start && 1 == timestep_mode) {

	  dt *= galdef.timestep_factor;
          
	  ostringstream nTBuf;
   
	  nTBuf<<" galdef_ID "<<galdef.galdef_ID<<" new timestep dt="<<dt<<"  "
	      <<particle.name;
	  DEBUGLOG(nTBuf.str());
	  
	  alpha1_x *= galdef.timestep_factor;
	  alpha2_x *= galdef.timestep_factor;
	  alpha3_x *= galdef.timestep_factor;
	  alpha1_y *= galdef.timestep_factor;
	  alpha2_y *= galdef.timestep_factor;
	  alpha3_y *= galdef.timestep_factor;
	  alpha1_z *= galdef.timestep_factor;
	  alpha2_z *= galdef.timestep_factor;
	  alpha3_z *= galdef.timestep_factor;
	  alpha1_p *= galdef.timestep_factor;
	  alpha2_p *= galdef.timestep_factor;
	  alpha3_p *= galdef.timestep_factor;
	
	}
	
	start = 0;
	double ADI_factor = 1;
	
	//  IF (ADI ==1) ADI_factor = 1.0/3.0
	/* causes memory leak
	   Nx1_ =             -  0.5* alpha1_x * ADI_factor;
	   Nx2_ =   1.0 +        0.5* alpha2_x * ADI_factor;              
	   Nx3_ =             -  0.5* alpha3_x * ADI_factor;
	   

	   Ny1_ =             -  0.5* alpha1_y * ADI_factor;
	   Ny2_ =   1.0 +        0.5* alpha2_y * ADI_factor;              
	   Ny3_ =             -  0.5* alpha3_y * ADI_factor;
	   
	   
	   Nz1_ =               -0.5* alpha1_z  * ADI_factor;
	   Nz2_ =   1.0 +        0.5* alpha2_z  * ADI_factor; 
	   Nz3_ =               -0.5* alpha3_z  * ADI_factor;
	   
	   Np1_ =              - 0.5* alpha1_p * ADI_factor;
	   Np2_ =   1.0 +        0.5* alpha2_p * ADI_factor;  
	   Np3_ =              - 0.5* alpha3_p * ADI_factor;
	*/
	
	Nx1_ = alpha1_x;   Nx1_*= -0.5;
	Nx2_ = alpha2_x;   Nx2_*=  0.5;   Nx2_+= 1.0;   
	Nx3_ = alpha3_x;   Nx3_*= -0.5;
	Ny1_ = alpha1_y;   Ny1_*= -0.5;
	Ny2_ = alpha2_y;   Ny2_*=  0.5;   Ny2_+= 1.0;   
	Ny3_ = alpha3_y;   Ny3_*= -0.5;
	Nz1_ = alpha1_z;   Nz1_*= -0.5;
	Nz2_ = alpha2_z;   Nz2_*=  0.5;   Nz2_+= 1.0;   
	Nz3_ = alpha3_z;   Nz3_*= -0.5;
	Np1_ = alpha1_p;   Np1_*= -0.5;
	Np2_ = alpha2_p;   Np2_*=  0.5;   Np2_+= 1.0;   
	Np3_ = alpha3_p;   Np3_*= -0.5;
	
	buf.str("");
	buf<<" Nx,Ny,Nz,Np recomputed ";
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Nx1_    :"<< Nx1_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Nx2_    :"<< Nx2_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Nx3_    :"<< Nx3_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Ny1_    :"<< Ny1_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
       	buf<<"Ny2_    :"<< Ny2_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Ny3_    :"<< Ny3_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Nz1_    :"<< Nz1_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Nz2_    :"<< Nz2_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Nz3_    :"<< Nz3_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Np1_    :"<< Np1_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Np2_    :"<< Np2_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	buf.str("");
	buf<<"Np3_    :"<< Np3_    .d3[1][1][1].s[1];
	DEBUGLOG(buf.str());
	
	// total_source_function= 1.0*particle.primary_source_function+ particle.secondary_source_function;
	total_source_function = particle.primary_source_function;
	total_source_function += particle.secondary_source_function;
            
	if (0 == total_source_function.max()) 
	  return 0;
	
	int nrept;
	if (1 == timestep_mode) 
	  nrept = galdef.timestep_repeat;

	if (2 == timestep_mode) 
	  nrept = galdef.timestep_repeat2;
	
	if (2 == timestep_mode && 1 == galdef.vectorized)
	  nrept = 1; // since loop is inside protri AWS20001120 20011022

	for (int irept = 1; irept <= nrept; ++irept) {
	  
	  if (1 == galdef.SNR_events && 2 == timestep_mode) {

	    t += dt;   
	    buf.str();
	    buf<<"timestep_mode 2: t="<<t;
	    INFO(buf.str());
	    source_SNR_event(particle, t/year2sec); // AWS20000807
	    total_source_function = particle.primary_source_function;
	    total_source_function += particle.secondary_source_function;

	    if(galdef.verbose>=1) {

	       buf.str("");
	      buf<<"propel:source_SNR_event: primary source function for particle " <<particle.name;
	      DEBUGLOG(buf.str());
	      particle.primary_source_function .print();
	    
	    }
            
	  }   // SNR_events

	  int n_xgrid_sym = particle.n_xgrid;
	  int n_ygrid_sym = particle.n_ygrid;
	  int n_zgrid_sym = particle.n_zgrid;

	  if (galdef.use_symmetry == 2) {
	  
	    n_xgrid_sym = particle.n_xgrid/2;
	    n_ygrid_sym = particle.n_ygrid/2;
	    n_zgrid_sym = particle.n_zgrid/2;
            
	  }
	  
	  // NON-VECTORIZED mode
	  
	  if (galdef.vectorized==0){                     //AWS20001120

	    // X propagation
	    
 //Gulli20070821
#pragma omp parallel default(shared) private(iz,ix,iy,ip)
	    if (1 == galdef.prop_x) {

	      //double *Nx0 = new double[particle.n_xgrid];    //Gulli20070821
	      //double *Nx1 = new double[particle.n_xgrid];    //Gulli20070821
	      //double *Nx2 = new double[particle.n_xgrid];    //Gulli20070821
	      //double *Nx3 = new double[particle.n_xgrid];    //Gulli20070821
	      //double *Rx = new double[particle.n_xgrid];    //Gulli20070821

	      valarray<double> Nx0(0., particle.n_xgrid);
	      valarray<double> Nx1(0., particle.n_xgrid);
	      valarray<double> Nx2(0., particle.n_xgrid);
	      valarray<double> Nx3(0., particle.n_xgrid);
	      valarray<double> Rx(0., particle.n_xgrid);
	      
 //Gulli20070821
#pragma omp for schedule(dynamic)
	      for (ip = 0; ip < particle.n_pgrid; ++ip) {
		
		for (iz = 0; iz < particle.n_zgrid; ++iz) {
		  
		  for (iy = 0; iy < particle.n_ygrid; ++iy) {

		    if (galdef.use_symmetry <= 1 || 
			(iy <= n_ygrid_sym && 
			 iz <= n_zgrid_sym)) {
                         
		      for (ix = 0; ix < particle.n_xgrid; ++ix) {

			Nx1[ix] = Nx1_.d3[ix][iy][iz].s[ip];
			Nx2[ix] = Nx2_.d3[ix][iy][iz].s[ip];
			Nx3[ix] = Nx3_.d3[ix][iy][iz].s[ip];
			Nx0[ix] = total_source_function.d3[ix][iy][iz].s[ip]*dt/f_use;

			Nx0[ix] += (1. - 0.5*alpha2_x.d3[ix][iy][iz].s[ip])*particle.cr_density.d3[ix][iy][iz].s[ip];
                     
			if (ix > 0) 
			  Nx0[ix] += 0.5*alpha1_x.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix-1][iy][iz].s[ip];
                    
			if (1 == galdef.use_symmetry && 0 == ix) 
			  Nx0[ix] += 0.5*alpha1_x.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix+1][iy][iz].s[ip];
                                 
			if (ix < particle.n_xgrid-1)
                                 
			  Nx0[ix] += 0.5*alpha3_x.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix+1][iy][iz].s[ip];
                              
		      }
                      
		      if (1 != galdef.use_symmetry) 
			tridag(&Nx1[0], &Nx2[0], &Nx3[0], &Nx0[0], &Rx[0], particle.n_xgrid);  
		      
		      if (1 == galdef.use_symmetry)
			tridag_sym(&Nx1[0], &Nx2[0], &Nx3[0], &Nx0[0], &Rx[0], particle.n_xgrid); 

		      for (ix = 0; ix < particle.n_xgrid; ++ix)
			particle.cr_density.d3[ix][iy][iz].s[ip] = (Rx[ix]< 0. ? 0.: Rx[ix]); //IMOS20030217
            
		    }   //symmetry
                    
		  }   //iy
                  
		}   //iz
                
	      }   //ip

	      if (2 == galdef.use_symmetry) {
 //Gulli20070821
#pragma omp for schedule(dynamic)
		for (ip = 0; ip < particle.n_pgrid; ++ip) {

		  for (iz = 0; iz <= n_zgrid_sym; ++iz) {

		    for (iy = 0; iy <= n_ygrid_sym; ++iy) {

		      for (ix = 0; ix < particle.n_xgrid; ++ix) {

			double value = particle.cr_density.d3[ix][iy][iz].s[ip];  //Gulli20070810
			particle.cr_density.d3[ix][particle.n_ygrid-1-iy][iz].s[ip] = value;
			
			particle.cr_density.d3[ix][iy][particle.n_zgrid-1-iz].s[ip] = value;
			
			particle.cr_density.d3[ix][particle.n_ygrid-1-iy][particle.n_zgrid-1-iz].s[ip] = value;
                         
		      }   //ix
		    
		    }   //iy
                    
		  }   //iz
                  
		}   //ip
                
	      }
	      //delete[] Nx0;    //Gulli20070821
	      //delete[] Nx1;    //Gulli20070821
	      //delete[] Nx2;    //Gulli20070821
	      //delete[] Nx3;    //Gulli20070821
	      //delete[] Rx;    //Gulli20070821
	    }   // prop_x

	    // Y propagation

	    //Gulli20070821
#pragma omp parallel default(shared) private(iz,ix,iy,ip)
	    if (1 == galdef.prop_y) {

	      //double *Ny0 = new double[particle.n_ygrid];    //Gulli20070821
	      //double *Ny1 = new double[particle.n_ygrid];    //Gulli20070821
	      //double *Ny2 = new double[particle.n_ygrid];    //Gulli20070821
	      //double *Ny3 = new double[particle.n_ygrid];    //Gulli20070821
	      //double *Ry = new double[particle.n_ygrid];    //Gulli20070821

	      valarray<double> Ny0(0., particle.n_ygrid);
	      valarray<double> Ny1(0., particle.n_ygrid);
	      valarray<double> Ny2(0., particle.n_ygrid);
	      valarray<double> Ny3(0., particle.n_ygrid);
	      valarray<double> Ry(0., particle.n_ygrid);

	      //Gulli20070821
#pragma omp for schedule(dynamic)
	      for (ip = 0; ip < particle.n_pgrid; ++ip) {
		
		for (iz = 0; iz < particle.n_zgrid; ++iz) {
		  
		  for (ix = 0; ix < particle.n_xgrid; ++ix) {

		    if (galdef.use_symmetry <= 1 || 
			(iz <= n_zgrid_sym && 
			 ix <= n_xgrid_sym)) {

		      for (iy = 0; iy < particle.n_ygrid; ++iy) {

			Ny1[iy] = Ny1_.d3[ix][iy][iz].s[ip];
			Ny2[iy] = Ny2_.d3[ix][iy][iz].s[ip];
			Ny3[iy] = Ny3_.d3[ix][iy][iz].s[ip];
			Ny0[iy] = total_source_function.d3[ix][iy][iz].s[ip]*dt/f_use;
                        
			Ny0[iy] += (1. - 0.5*alpha2_y.d3[ix][iy][iz].s[ip])*particle.cr_density.d3[ix][iy][iz].s[ip];
                         
			if (iy > 0) 
			  Ny0[iy] += 0.5*alpha1_y.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix][iy-1][iz].s[ip];

			if (1 == galdef.use_symmetry && 0 == iy) 
			  Ny0[iy] += 0.5*alpha1_y.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix][iy+1][iz].s[ip];

			if (iy < particle.n_ygrid-1) 
			  Ny0[iy] += 0.5*alpha3_y.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix][iy+1][iz].s[ip];
                              
		      }

		      if (1 != galdef.use_symmetry) 
			tridag(&Ny1[0], &Ny2[0], &Ny3[0], &Ny0[0], &Ry[0], particle.n_ygrid);  
		      
		      if (1 == galdef.use_symmetry) 
			tridag_sym(&Ny1[0], &Ny2[0], &Ny3[0], &Ny0[0], &Ry[0], particle.n_ygrid);  
		      
		      for (iy = 0; iy < particle.n_ygrid; ++iy) 
			particle.cr_density.d3[ix][iy][iz].s[ip] = (Ry[iy] < 0. ? 0.: Ry[iy]); //IMOS20030217
             
		    }  //  symmetry

		  }  //  ix
                  
		}  //  iz
                
	      }  //  ip
	      
	      if (2 == galdef.use_symmetry) {
		//Gulli20070821
#pragma omp for schedule(dynamic)
		for (ip = 0; ip < particle.n_pgrid; ++ip) {

		  for (iz = 0; iz <= n_zgrid_sym; ++iz) {
		    
		    for (ix = 0; ix <= n_xgrid_sym; ++ix) {

		      for (iy = 0; iy < particle.n_ygrid; ++iy) {

			double value =  particle.cr_density.d3[ix][iy][iz].s[ip];
			particle.cr_density.d3[particle.n_xgrid-1-ix][iy][iz].s[ip] = value;
			particle.cr_density.d3[ix][iy][particle.n_zgrid-1-iz].s[ip] = value;
			particle.cr_density.d3[particle.n_xgrid-1-ix][iy][particle.n_zgrid-1-iz].s[ip] = value;
                         
		      }  //  iy    
                      
		    }  //  ix     
                    
		  }  //  iz
                  
		}  //  ip
                
	      }  //  symmetry
	      //delete[] Ny0;    //Gulli20070821
	      //delete[] Ny1;    //Gulli20070821
	      //delete[] Ny2;    //Gulli20070821
	      //delete[] Ny3;    //Gulli20070821
	      //delete[] Ry;    //Gulli20070821
	    }  //  prop_y

	    // Z propagation
	    
	    //Gulli20070821
#pragma omp parallel default(shared) private(iz,ix,iy,ip)
	    if (1 == galdef.prop_z) {
	      
	      //double *Nz0 = new double[particle.n_zgrid];    //Gulli20070821
	      //double *Nz1 = new double[particle.n_zgrid];    //Gulli20070821
	      //double *Nz2 = new double[particle.n_zgrid];    //Gulli20070821
	      //double *Nz3 = new double[particle.n_zgrid];    //Gulli20070821
	      //double *Rz = new double[particle.n_zgrid];    //Gulli20070821
	      
	      valarray<double> Nz0(0., particle.n_zgrid);
	      valarray<double> Nz1(0., particle.n_zgrid);
	      valarray<double> Nz2(0., particle.n_zgrid);
	      valarray<double> Nz3(0., particle.n_zgrid);
	      valarray<double> Rz(0., particle.n_zgrid);
	      
	      //Gulli20070821
#pragma omp for schedule(dynamic)
	      for (ip = 0; ip < particle.n_pgrid; ++ip) {

		for (ix = 0; ix < particle.n_xgrid; ++ix) {

		  for (iy = 0; iy < particle.n_ygrid; ++iy) {

		    if (galdef.use_symmetry <= 1 || 
			(ix <= n_xgrid_sym && 
			 iy <= n_ygrid_sym)) {

		      for (iz = 0; iz < particle.n_zgrid; ++iz) {

			Nz1[iz] = Nz1_.d3[ix][iy][iz].s[ip];
			Nz2[iz] = Nz2_.d3[ix][iy][iz].s[ip];
			Nz3[iz] = Nz3_.d3[ix][iy][iz].s[ip];
			Nz0[iz] = total_source_function.d3[ix][iy][iz].s[ip]*dt/f_use;

			Nz0[iz] += (1. - 0.5*alpha2_z.d3[ix][iy][iz].s[ip])*particle.cr_density.d3[ix][iy][iz].s[ip];
                            
			if (iz > 0)
			  Nz0[iz] += 0.5*alpha1_z.d3[ix][iy][iz].s[ip] *particle.cr_density.d3[ix][iy][iz-1].s[ip];
			
			if (galdef.use_symmetry && iz==0)
			  Nz0[iz] += 0.5*alpha1_z.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix][iy][iz+1].s[ip];

			if (iz <particle.n_zgrid-1)
			  Nz0[iz] += 0.5*alpha3_z.d3[ix][iy][iz].s[ip] *particle.cr_density.d3[ix][iy][iz+1].s[ip];

		      } 
// cout<<"ix iy ip "<<ix<<" "<<iy<<" "<<ip<<endl;
                      
		      if (1 != galdef.use_symmetry) 
			tridag(&Nz1[0], &Nz2[0], &Nz3[0], &Nz0[0], &Rz[0], particle.n_zgrid);  
		      
		      if (1 == galdef.use_symmetry) 
			tridag_sym(&Nz1[0], &Nz2[0], &Nz3[0], &Nz0[0], &Rz[0], particle.n_zgrid);  

		      for (iz = 0; iz < particle.n_zgrid; ++iz) 
			particle.cr_density.d3[ix][iy][iz].s[ip] = (Rz[iz] < 0. ? 0.: Rz[iz]); //IMOS20030217
             
		    }  //  symmetry

		  }  //  iy
                  
		}  //  ix
                
	      }  //  ip
	      
	      if (2 == galdef.use_symmetry) {
		//Gulli20070821
#pragma omp for schedule(dynamic)
		for (ip = 0; ip < particle.n_pgrid; ++ip) {

		  for (ix = 0; ix <= n_xgrid_sym; ++ix) {

		    for (iy = 0; iy <= n_ygrid_sym; ++iy) {

		      for (iz = 0; iz < particle.n_zgrid; ++iz) {

			double value = particle.cr_density.d3[ix][iy][iz].s[ip];
			
			particle.cr_density.d3[particle.n_xgrid-1-ix][iy][iz].s[ip] = value;
			particle.cr_density.d3[ix][particle.n_ygrid-1-iy][iz].s[ip] = value;
			particle.cr_density.d3[particle.n_xgrid-1-ix][particle.n_ygrid-1-iy][iz].s[ip] = value;
                          
		      }  //  iz
                      
		    }  //  iy 
                    
		  }  //  ix
                  
		}  //  ip
                
	      }  //  symmetry
	      //delete[] Nz0;    //Gulli20070821
	      //delete[] Nz1;    //Gulli20070821
	      //delete[] Nz2;    //Gulli20070821
	      //delete[] Nz3;    //Gulli20070821
	      //delete[] Rz;    //Gulli20070821
	    }  //  prop_z

	    // P propagation

	    //Gulli20070821
#pragma omp parallel default(shared) private(iz,ix,iy,ip)
	    if (1 == galdef.prop_p) { 
	      //double *Np0 = new double[particle.n_pgrid];    //Gulli20070821
	      //double *Np1 = new double[particle.n_pgrid];    //Gulli20070821
	      //double *Np2 = new double[particle.n_pgrid];    //Gulli20070821
	      //double *Np3 = new double[particle.n_pgrid];    //Gulli20070821
	      //double *Rp = new double[particle.n_pgrid];    //Gulli20070821

	      valarray<double> Np0(0., particle.n_pgrid);
	      valarray<double> Np1(0., particle.n_pgrid);
	      valarray<double> Np2(0., particle.n_pgrid);
	      valarray<double> Np3(0., particle.n_pgrid);
	      valarray<double> Rp(0., particle.n_pgrid);

	      //Gulli20070821
#pragma omp for schedule(dynamic)
	      for (ix = 0; ix < particle.n_xgrid; ++ix) {

		for (iy = 0; iy < particle.n_ygrid; ++iy) {

		  for (iz = 0; iz < particle.n_zgrid; ++iz) {
		    
		    if (galdef.use_symmetry <= 1 || 
			(ix <= n_xgrid_sym && 
			 iy <= n_ygrid_sym && 
			 iz <= n_zgrid_sym)) {

		      for (ip = 0; ip < particle.n_pgrid; ++ip) {

			Np1[ip] = Np1_.d3[ix][iy][iz].s[ip];
			Np2[ip] = Np2_.d3[ix][iy][iz].s[ip];
			Np3[ip] = Np3_.d3[ix][iy][iz].s[ip];
			
			Np0[ip] = total_source_function.d3[ix][iy][iz].s[ip]*dt/f_use;
                                 
			Np0[ip] += (1. - 0.5*alpha2_p.d3[ix][iy][iz].s[ip])*particle.cr_density.d3[ix][iy][iz].s[ip];
			
			if (ip > 0)
			  Np0[ip] += 0.5*alpha1_p.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix][iy][iz].s[ip-1];
			
			if (ip < particle.n_pgrid-1)
			  Np0[ip] += 0.5*alpha3_p.d3[ix][iy][iz].s[ip]*particle.cr_density.d3[ix][iy][iz].s[ip+1];
                              
		      }
		      
		      tridag(&Np1[0], &Np2[0], &Np3[0], &Np0[0], &Rp[0], particle.n_pgrid);  
		      
		      for (ip = 0; ip < particle.n_pgrid; ++ip) 
			particle.cr_density.d3[ix][iy][iz].s[ip] = (Rp[ip] < 0. ? 0.: Rp[ip]); //IMOS20030217
		    
		    }  //  symmetry

		  }  //  iz
		
		}  //  iy
                
	      }  //  ix
	       
	      if (2 == galdef.use_symmetry) {

		//Gulli20070821
#pragma omp for schedule(dynamic)
		for (ip = 0; ip < particle.n_pgrid; ++ip) {

		  for (ix = 0; ix <= n_xgrid_sym; ++ix) {

		    for (iy = 0; iy <= n_ygrid_sym; ++iy) {

		      for (iz = 0; iz <= n_zgrid_sym; ++iz) {

			double value = particle.cr_density.d3[ix][iy][iz].s[ip];
			
			particle.cr_density.d3[ix][particle.n_ygrid-1-iy][iz].s[ip] = value;
			particle.cr_density.d3[ix][iy][particle.n_zgrid-1-iz].s[ip] = value;
			particle.cr_density.d3[ix][particle.n_ygrid-1-iy][particle.n_zgrid-1-iz].s[ip] = value;
			particle.cr_density.d3[particle.n_xgrid-1-ix][iy][iz].s[ip] = value;
			particle.cr_density.d3[particle.n_xgrid-1-ix][iy][particle.n_zgrid-1-iz].s[ip] = value;
			particle.cr_density.d3[particle.n_xgrid-1-ix][particle.n_ygrid-1-iy][iz].s[ip] = value;
			particle.cr_density.d3[particle.n_xgrid-1-ix][particle.n_ygrid-1-iy][particle.n_zgrid-1-iz].s[ip] = value; 
		      
		      }  //  iz   
		    
		    }  //  iy      
                    
		  }  //  ix
                  
		}  //  ip
                
	      }  //  symmetry
	      //delete[] Np0;    //Gulli20070821
	      //delete[] Np1;    //Gulli20070821
	      //delete[] Np2;    //Gulli20070821
	      //delete[] Np3;    //Gulli20070821
	      //delete[] Rp;    //Gulli20070821
               
	    }  //  prop_p
	    
	  }// vectorized==0     AWS20001120

	  
	  // VECTORIZED mode

	  /*if(galdef.vectorized==1){                     //AWS20001120
	    
	    INFO(" using vectorized routines");
	    
	    protri(particle,
		   
		   alpha1_x,  alpha1_y,  alpha1_z,  alpha1_p,
		   alpha2_x,  alpha2_y,  alpha2_z,  alpha2_p,
		   alpha3_x,  alpha3_y,  alpha3_z,  alpha3_p,
		   
		   Nx1_,  Ny1_,  Nz1_, Np1_,
		   Nx2_,  Ny2_,  Nz2_,  Np2_,
		   Nx3_,  Ny3_,  Nz3_, Np3_,
		   
		   total_source_function,
		   
		   dt, galdef.timestep_repeat2,  f_use );
	    
	  }// vectorized==1     AWS20001120
	  */

	  // set BOUNDARY CONDITIONS: zero on all boundaries
	  
	  //Gulli20070821
#pragma omp parallel for schedule(dynamic) default(shared) private(iz,ix,iy,ip)
	  for (ip = 0; ip < particle.n_pgrid; ++ip) {
	    
	    for (ix = 0; ix < particle.n_xgrid; ++ix) {

	      for (iy = 0; iy < particle.n_ygrid; ++iy) {

		if (1 != galdef.use_symmetry) { 
	
		  iz = 0;  
		  particle.cr_density.d3[ix][iy][iz].s[ip] = 0; 
		
		}
		  
		iz = particle.n_zgrid-1;
		particle.cr_density.d3[ix][iy][iz].s[ip] = 0;

	      }
	     
	    }
	     
	    for (ix = 0; ix < particle.n_xgrid; ++ix) {
	      
	      for (iz = 0; iz < particle.n_zgrid; ++iz) {
		  
		if (1 != galdef.use_symmetry) { 
		  
		  iy = 0;  
		  particle.cr_density.d3[ix][iy][iz].s[ip] = 0; 

		}
		
		iy = particle.n_ygrid-1;
		particle.cr_density.d3[ix][iy][iz].s[ip] = 0;

	      }	      

	    }
	     
	    for (iz = 0; iz < particle.n_zgrid; ++iz) {
	      
	      for (iy = 0; iy < particle.n_ygrid; ++iy) {
		
		if (1 != galdef.use_symmetry) { 

		  ix = 0;  
		  particle.cr_density.d3[ix][iy][iz].s[ip] = 0; 

		}
		
		ix = particle.n_xgrid-1;
		particle.cr_density.d3[ix][iy][iz].s[ip] = 0;
		
	      }
	    
	    }

	  }  //  boundary conditions
	  
	  if (0 == irept%galdef.timestep_print) {
	    
	    buf.str("");
	    buf<<"particle.cr_density after iteration "<<irept<<" with time step "<<dt/year2sec<<" yr";
	    INFO(buf.str());
	    particle.cr_density.print();
	    cout<<endl;
	  
	  }
	  
	  if (0 == irept%galdef.timestep_diagnostics) 
	    propel_diagnostics(particle,
			       alpha1_x, alpha1_y, alpha1_z, alpha1_p,
			       alpha2_x, alpha2_y, alpha2_z, alpha2_p,
			       alpha3_x, alpha3_y, alpha3_z, alpha3_p,
			       total_source_function, 
			       dt);
	}  //  irept

	if (2 == timestep_mode) 
	  timestep_mode2done = 1;
      
      }  //  dt>=end_timestep_sec
    
    }  //  timestep_mode
    
    //Delete temporary arrays         //Gulli20070810
    alpha1_x.delete_array();
    alpha1_y.delete_array();
    alpha1_z.delete_array();
    alpha1_p.delete_array();
    alpha2_x.delete_array();
    alpha2_y.delete_array();
    alpha2_z.delete_array();
    alpha2_p.delete_array();
    alpha3_x.delete_array();
    alpha3_y.delete_array();
    alpha3_z.delete_array();
    alpha3_p.delete_array();
    
    Nx1_.delete_array();
    Ny1_.delete_array();
    Nz1_.delete_array();
    Np1_.delete_array();
    Nx2_.delete_array();
    Ny2_.delete_array();
    Nz2_.delete_array();
    Np2_.delete_array();
    Nx3_.delete_array();
    Ny3_.delete_array();
    Nz3_.delete_array();
    Np3_.delete_array();
    
    total_source_function.delete_array();
    Work.delete_array();

  } //3D case
  // for method see notebook of 970300
  
  INFO("Exit");
  
  //cout<<"\n<<<<propel"<<endl;
  return 0;
}
