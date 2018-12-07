#include "galprop_classes.h"
#include "galprop_internal.h"
#include <vector>
#include <map>
#include <cstring>
#include <cstdlib>

static bool outputErrorMessage = true;

//Convenience functions
void retrieveGCRData(vector<double> &eMean, vector<double> &flux, vector<double> &errMinus, vector<double> &errPlus, const GCR_data &gcr_data, vector<int> Z_denom, vector<int> A_denom, vector<int> Z_numer, vector<int> A_numer, char *type, const vector<string> &references = vector<string>(0) ) {
   for (int i = 0; i < gcr_data.n; ++i){
      //Apply all the requirements
      bool include = strcmp(gcr_data.Y_name[i], type) == NULL;
      //The size matters
      size_t ds=0;
      while (ds < 3 && gcr_data.Z_denominator[i][ds] != 0) ++ds;
      size_t ns=0;
      while (ns < 3 && gcr_data.Z_numerator[i][ns] != 0) ++ns;
      include &= Z_denom.size() == ds;
      include &= Z_numer.size() == ns;
      //The order matters
      for (int j = 0; j < Z_denom.size(); ++j) include &= Z_denom[j] == gcr_data.Z_denominator[i][j]; 
      for (int j = 0; j < A_denom.size(); ++j) include &= A_denom[j] == gcr_data.A_denominator[i][j]; 
      for (int j = 0; j < Z_numer.size(); ++j) include &= Z_numer[j] == gcr_data.Z_numerator[i][j]; 
      for (int j = 0; j < A_numer.size(); ++j) include &= A_numer[j] == gcr_data.A_numerator[i][j]; 
      //All the references should be included
      bool instr = references.size() == 0;
      for (int j = 0; j < references.size(); ++j) instr |= references[j] == gcr_data.reference[i];
      include &= instr;
      if (include && gcr_data.value[i] != 0){
	 //Calculate the average E_mean if it does not exist
	 if (gcr_data.E_mean[i] == 0) {
	    eMean.push_back(sqrt(gcr_data.E_low[i]*gcr_data.E_high[i]));
      	 } else {
	    eMean.push_back(gcr_data.E_mean[i]);
	 }
	 flux.push_back(gcr_data.value[i]);
	 //Handle unsigned errors in either minus or plus or both
	 if (gcr_data.err_minus[i] == 0) {
	    if (gcr_data.err_plus[i] == 0) {
	       errMinus.push_back(0.1*gcr_data.value[i]);
	       errPlus.push_back(0.1*gcr_data.value[i]);
	    }else{
	       errMinus.push_back(gcr_data.err_plus[i]);
	       errPlus.push_back(gcr_data.err_plus[i]);
	    }
	 }else{
	    errMinus.push_back(gcr_data.err_minus[i]);
	    if (gcr_data.err_plus[i] == 0){
	       errPlus.push_back(gcr_data.err_minus[i]);
	    }else{
	       errPlus.push_back(gcr_data.err_plus[i]);
	    }
	 }
      }
   }
}

std::vector<double> interpolateAndModulate(vector<double> spectra, const double *energy, const vector<double> &dataEmean, int Z, int A, double phi){
   //Constants needed for modulation
   const double eMass = 0.5109990615; //Electron mass in MeV;
   const double pMass = 939.; //Proton mass in MeV

   vector<double> output(dataEmean.size(),0.0);

   for (int i = 0; i < dataEmean.size(); ++i){
      double en;
      if ( A == 0) {
	 en = dataEmean[i] + fabs(Z)*phi;
      } else {
	 en = dataEmean[i] + fabs(Z)*phi/float(A);
      }
      //Do interpolation in energy, but no extrapolation
      if (en >= energy[0] && en <= energy[spectra.size()-1]) {
	 int j(0);
	 while (en > energy[j]) ++j;
	 if (j == 0) ++j; //If en is excactly the lower limit
	 //Power law interpolation if both values exist, otherwise linear
	 if (spectra[j-1] > 0 && spectra[j] > 0) {
	    double sl = log(spectra[j]/spectra[j-1])/log(energy[j]/energy[j-1]);
	    output[i] = spectra[j]*pow(en/energy[j], sl);
	 }else if (spectra[j-1] > 0) {
	    output[i] = (energy[j] - en)*spectra[j-1]/(energy[j]-energy[j-1]);
	 }else if (spectra[j] > 0) {
	    output[i] = (en - energy[j-1])*spectra[j]/(energy[j]-energy[j-1]);
	 }
	 if ( A == 0) {
	    output[i] *= dataEmean[i]*(dataEmean[i]+2*eMass)/(en*(en+2*eMass));
	 } else {
	    output[i] *= dataEmean[i]*(dataEmean[i]+2*pMass)/(en*(en+2*pMass));
	 }
      } else {
	 output[i] = -1;
	 if (outputErrorMessage){
	    cerr<<"Energy range of model is too small for data"<<endl;
	    cerr<<en<<" is not within range "<<energy[0]<<" - "<<energy[spectra.size()-1]<<endl;
	    cerr<<"This message will not be repeated"<<endl;
	    outputErrorMessage = false;
	 }
      }
   }
   return output;
}

std::vector<double> readNuclei(const Particle *gcr, int n_species, int Z, int A, int iz1, int iz2, int ir1, int ir2, int ix1, int ix2, int iy1, int iy2, double fz1, double fz2, double fr1, double fr2, double fx1, double fx2, double fy1, double fy2){
	vector<double> output(0);
	for (int i = 0; i < n_species; ++i){
		if (gcr[i].Z == Z && gcr[i].A == A){
		 	if (output.size() == 0) output.resize(gcr[0].n_pgrid, 0.0);
			for (int j = 0; j < gcr[0].n_pgrid;++j) {
				//Do interpolation in z,r,x,y 
				if (gcr[i].n_spatial_dimensions==2){
					output[j]+=fr1*fz1*gcr[i].cr_density.d2[ir1][iz1].s[j]+fr2*fz1*gcr[i].cr_density.d2[ir2][iz1].s[j]+
                     fr1*fz2*gcr[i].cr_density.d2[ir1][iz2].s[j]+fr2*fz2*gcr[i].cr_density.d2[ir2][iz2].s[j];
				}else if (gcr[i].n_spatial_dimensions==3){
					output[j]+=fx1*fy1*fz1*gcr[i].cr_density.d3[ix1][iy1][iz1].s[j]+fx2*fy1*fz1*gcr[i].cr_density.d3[ix2][iy1][iz1].s[j]+
						         fx1*fy2*fz1*gcr[i].cr_density.d3[ix1][iy2][iz1].s[j]+fx2*fy2*fz1*gcr[i].cr_density.d3[ix2][iy2][iz1].s[j]+
						         fx1*fy1*fz2*gcr[i].cr_density.d3[ix1][iy1][iz2].s[j]+fx2*fy1*fz2*gcr[i].cr_density.d3[ix2][iy1][iz2].s[j]+
						         fx1*fy2*fz2*gcr[i].cr_density.d3[ix1][iy2][iz2].s[j]+fx2*fy2*fz2*gcr[i].cr_density.d3[ix2][iy2][iz2].s[j];
				}
			}
		}
	}
	return output;
}

//Finds the value of il and iu such that arr[il] < value < arr[iu].  Size of
//array is size, the step size in arr is da, and it should be monotomically
//increasing.
//If value lies outside of the boundaries, an error is printed and the closest
//value is chosen. il is in this case put to the closest pixel, iu to -1, fl to
//1 and fu to 0
//If the value falls within 0.1% of a arr element, it is chosen as il, fl set
//to 1, fu to 0 and iu to 0
void findInterpolationBounds(const double *arr, double da, int size, int &il, int &iu, double &fl, double &fu, double value, const std::string & type) {
	il=(int)((value-arr[0])/da);
	iu = 0;
	fl = 1; 
	fu = 0;
	if (il >= 0 && il < size) { //Should always be the case
		//Do proper interpolation if needed
		if ( fabs(value-arr[il])/da > 1e-3 ) {
			if (value-arr[il] < 0) { //Need to have a lower value for ir
				if (il == 0) {
					std::cerr<<type<<" grid not wide enough to include the sun, using closest grid point"<<std::endl;
				} else {
					iu = il;
					--il;
					fl = (value-arr[iu])/da;
					fu = 1-fl;
				}
			} else {
				if (il == size-1) {
					std::cerr<<type<<" grid not wide enough to include the sun, using closest grid point"<<std::endl;
				} else {
					iu = il+1;
					fu = (value-arr[il])/da;
					fl = 1-fu;
				}
			}
		}
	} else {
		//Just post an error and use the last grid point
		std::cerr<<type<<" grid not wide enough to include the sun, using closest grid point"<<std::endl;
		if (il < 0) {
			il = 0;
		} else {
			il = size-1;
		}
	}
}

// The bool sets determine which datasets to include.  The available ones are in
// order:
//   0: e- (HEAT, SANRIKU98, AMS )
//   1: e+ 
//   2: H (BESS[Wa02, Ha04, Sa00], ATIC2, AMS, JACEE)
//   3: pbar (BESS[Or00], PAMELA[Ad10]
//   4: He (BESS[Wa02, Ha04, Sa00], ATIC2, AMS, JACEE)
//   5: O (ACE, HEAO3)
//   6: C (ACE, HEAO3)
//   7: B/C ratio (ACE (Da00, No01), HEAO3, ATIC2, CREAM)
//   8: Be10/Be9 ratio (ACE)
// The names in parentheses indicate the experiments included.  Modulation for
// those must be set.
//
// The bool expmnts can be used to exclude/include experiments (references) in the
// calculation, in order
//   0: HEAO3 (En90)
//   1: ACE   (ACEW, Da00)
//   2: ATIC2 (Pa06, Pa07)
//   3: CREAM (Ahn8)
//   4: AMS   (Ag02)
//   5: BESS  (Ha04)
//   6: BESS  (Wa02)
//   7: JACEE (As98)
//   8: HEAT  (Du01)
//   9: SANRIKU98 (Ko99)
//  10: Fermi-LAT (Ab09, Abnn)
//  11: HESS (Ah08)
//  12: HESS09 (Ah09)
//  13: BESS (Sa00)
//  14: BESS (Or00)
//  15: PAMELA (Ad10)
// Note that no check is made for datasets with 0 data points included.  What
// happens then is currently unknown
//
// The double array modulation gives the modulation potential for each
// reference/experiment (same order as expmnts).  Suggested values are in
// parentheses
//   0: HEAO3 (600 +- 50)
//   1: ACE (325 +- 50)
//   2: ATIC2  (Not important since data at high energy, but needed as a parameter, 0
//   should be OK)
//   3: CREAM (Not very important since data at high energy, but needed as a
//   parameter, 600 +- 50 should be a safe estimate)
//   4: AMS  (650 +- 40 according to ag02)
//   5: BESS (1100 +- 100)
//   6: BESS (1300 +- 100)
//   7: JACEE (Not important since data at high energy, but needed as a
//   parameter, 0 should be OK)
//   8: HEAT (700 +- 50)
//   9: SANRIKU98 (Not important since data at high energy, but needed as a
//   parameter, 0 should be OK)
//  10: Fermi-LAT (300 should be OK)
//  11: HESS (Not important since data at high energy, but needed as a
//   parameter, 0 should be OK)
//  12: HESS09 (Not important since data at high energy, but needed as a
//   parameter, 0 should be OK)
//  13: BESS (400 +- 100 [solar minimum])
//  14: BESS (450 +- 100)
//  15: PAMELA (400 +- 100 [solar minimum])
//
// The double array tau sets a scaling factor for the error of each experiment (precision)
// such that new sigma_i^2 = sigma_i^2/tau.  The array should be set to the log10
// value of tau.  One parameter for each reference/experiment, same order as
// expmnts
// RT: modified chisquare to be -2loglike (including prefactor of log(tau)
//   0: HEAO3
//   1: ACE
//   2: ATIC2 
//   3: CREAM 
//   4: AMS 
//   5: BESS 
//   6: BESS
//   7: JACEE
//   8: HEAT
//   9: SANRIKU98
//  10: Fermi-LAT
//  11: HESS
//  12: HESS09
//  13: BESS
//  14: BESS
//  15: PAMELS
//
// Double array energyScale sets the energy scale for the experiments.  Same
// order as arrays above.
double Galprop::calculate_GCR_chisq(bool sets[], bool expmnts[], double modulation[], double tau[], double energyScale[], double chisqa[]) {

   double prt=0;

   //Read in the data if necessary
   if (gcr_data.n == 0) {                                                                                                                                                

      const string fullFilename = configure.fGlobalDataPath + "/" + galdef.GCR_data_filename;//"GCR_data_1.dat";                                                                                                
			std::cout<<"Filename for GCR database: "<<fullFilename<<std::endl;

      gcr_data.read(fullFilename.c_str(), "cm2", "MeV");                                                                                                                  
   }            


   //Find the iz value for z = 0;
   int iz1, iz2, ir1, ir2, ix1, ix2, iy1, iy2;
	 double fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2;
	 findInterpolationBounds(galaxy.z, galdef.dz, galaxy.n_zgrid, iz1, iz2, fz1, fz2, 0.0, "z");
   //3D and 2D are different
   if (gcr[0].n_spatial_dimensions==2) {
		 findInterpolationBounds(galaxy.r, galdef.dr, galaxy.n_rgrid, ir1, ir2, fr1, fr2, Rsun, "r");
   } else if (gcr[0].n_spatial_dimensions==3){
		 findInterpolationBounds(galaxy.x, galdef.dx, galaxy.n_xgrid, ix1, ix2, fx1, fx2, Rsun, "x");
		 findInterpolationBounds(galaxy.y, galdef.dy, galaxy.n_ygrid, iy1, iy2, fy1, fy2, 0.0, "y");
   }
	 //std::cout<<"Interpolation values \n z: ("<<galaxy.z[iz1]<<"["<<iz1<<"], "<<fz1<<"), ("<<galaxy.z[iz2]<<"["<<iz2<<"], "<<fz2<<")\n r: ("<<galaxy.r[ir1]<<"["<<ir1<<"], "<<fr1<<"), ("<<galaxy.r[ir2]<<"["<<ir2<<"], "<<fr2<<")\n";

   double chisq=0;
   //Check the sets
   if (sets[0]) { //electrons
      //Data storage for the electron spectrum
      vector<double> electrons = readNuclei(gcr,n_species,-1,0,iz1,iz2,ir1,ir2,ix1,ix2,iy1,iy2,fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2);
      if (electrons.size() == 0){
	 std::cerr<<"No electrons found in nuclei data, not using those in chisq"<<std::endl;
      } else {
	 //Select the appropriate data
	 std::vector<double> dataEmean, dataFlux, dataErrMinus, dataErrPlus;
	 vector<int> Z_d, A_d, Z_n, A_n;
	 vector<string> references(1); //Storage to select references
	 map<string,int> refs; //The actual references and mapping to modulation
	 if (expmnts[4]) refs["Ag02"] = 4;
	 if (expmnts[8]) refs["Du01"] = 8;
	 if (expmnts[9]) refs["Ko99"] = 9;
	 if (expmnts[10]) {
		 refs["Ab09"] = 10;
		 refs["Abnn"] = 10;
	 }
	 if (expmnts[11]) refs["Ah08"] = 11;
	 if (expmnts[12]) refs["Ah09"] = 12;
	 Z_d.push_back(-1);
	 A_d.push_back(0);
	 map<string,int>::iterator it = refs.begin();
	 for ( ; it != refs.end(); ++it) {
	    vector<double> dataEmean, dataFlux, dataErrMinus, dataErrPlus;
	    references[0] = it->first;
	    retrieveGCRData(dataEmean,dataFlux,dataErrMinus,dataErrPlus,gcr_data,Z_d,A_d,Z_n,A_n,"flux",references);
	    //Correct for energy scaling
	    for (int i = 0; i < dataEmean.size(); ++i) 
	       dataEmean[i] *= energyScale[it->second]; 
	    //Interpolate and modulate the model
	    vector<double> flux = interpolateAndModulate(electrons, gcr[0].Ekin, dataEmean, -1, 0, modulation[it->second]);
	    //Now we can calculate the chisq, including modulation
	    if (galdef.verbose == -701) cout<<"#Electron data "<<it->first<<":"<<endl;
	    //cout<<"Modulation:"<<modulation[it->second]<<endl;
	    for (int i = 0; i < dataEmean.size(); ++i){
	       if (flux[i] >= 0) {
		  if (flux[i] < dataFlux[i]) {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrMinus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrMinus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  } else {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrPlus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrPlus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  }
	       }
	    }
	    if (galdef.verbose == -701) cout<<endl<<endl;
	 }
         chisqa[0] = chisq;
         prt = chisqa[0];
      }
   }
   if (sets[1]) { //positrons
      //Data storage for the electron spectrum
      vector<double> positrons = readNuclei(gcr,n_species,1,0,iz1,iz2,ir1,ir2,ix1,ix2,iy1,iy2,fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2);
      if (positrons.size() == 0){
	 std::cerr<<"No positrons found in nuclei data, not using those in chisq"<<std::endl;
      } else {
	 //Select the appropriate data
	 std::vector<double> dataEmean, dataFlux, dataErrMinus, dataErrPlus;
	 vector<int> Z_d, A_d, Z_n, A_n;
	 vector<string> references(1); //Storage to select references
	 map<string,int> refs; //The actual references and mapping to modulation
	 if (expmnts[4]) refs["Ag02"] = 4;
	 if (expmnts[8]) refs["Du01"] = 8;
	 Z_d.push_back(1);
	 A_d.push_back(0);
	 map<string,int>::iterator it = refs.begin();
	 for ( ; it != refs.end(); ++it) {
	    vector<double> dataEmean, dataFlux, dataErrMinus, dataErrPlus;
	    references[0] = it->first;
	    retrieveGCRData(dataEmean,dataFlux,dataErrMinus,dataErrPlus,gcr_data,Z_d,A_d,Z_n,A_n,"flux",references);
	    //Correct for energy scaling
	    for (int i = 0; i < dataEmean.size(); ++i) 
	       dataEmean[i] *= energyScale[it->second]; 
	    //Interpolate and modulate the model
	    vector<double> flux = interpolateAndModulate(positrons, gcr[0].Ekin, dataEmean, -1, 0, modulation[it->second]);
	    //Now we can calculate the chisq, including modulation
	    if (galdef.verbose == -701) cout<<"#Positron data "<<it->first<<":"<<endl;
	    //cout<<"Modulation:"<<modulation[it->second]<<endl;
	    for (int i = 0; i < dataEmean.size(); ++i){
	       if (flux[i] >= 0) {
		  if (flux[i] < dataFlux[i]) {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrMinus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrMinus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  } else {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrPlus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrPlus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  }
	       }
	    }
	    if (galdef.verbose == -701) cout<<endl<<endl;
	 }
         chisqa[1] = chisq;
         prt = chisqa[1];
      }
   }
   if (sets[2]) { //protons
      //Data storage for the proton spectrum
      vector<double> protons = readNuclei(gcr,n_species,1,1,iz1,iz2,ir1,ir2,ix1,ix2,iy1,iy2,fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2);
      if (protons.size() == 0){
	 std::cerr<<"No protons found in nuclei data, not using those in chisq"<<std::endl;
      } else {
	 //Select the appropriate data
	 std::vector<double> dataEmean, dataFlux, dataErrMinus, dataErrPlus;
	 vector<int> Z_d, A_d, Z_n, A_n;
	 vector<string> references(1); //Storage to select references
	 map<string,int> refs; //The actual references and mapping to modulation
	 if (expmnts[2]) refs["Pa06"] = 2;
	 if (expmnts[4]) refs["Ag02"] = 4;
	 if (expmnts[5]) refs["Ha04"] = 5;
	 if (expmnts[6]) refs["Wa02"] = 6;
	 if (expmnts[7]) refs["As98"] = 7;
	 if (expmnts[13]) refs["Sa00"] = 13;
	 Z_d.push_back(1);
	 A_d.push_back(1);
	 map<string,int>::iterator it = refs.begin();
	 for ( ; it != refs.end(); ++it) {
	    vector<double> dataEmean, dataFlux, dataErrMinus, dataErrPlus;
	    references[0] = it->first;
	    retrieveGCRData(dataEmean,dataFlux,dataErrMinus,dataErrPlus,gcr_data,Z_d,A_d,Z_n,A_n,"flux",references);
	    //Correct for energy scaling
	    for (int i = 0; i < dataEmean.size(); ++i) 
	       dataEmean[i] *= energyScale[it->second]; 
	    //Interpolate and modulate the model
	    vector<double> flux = interpolateAndModulate(protons, gcr[0].Ekin, dataEmean, 1, 1, modulation[it->second]);
	    //Now we can calculate the chisq, including modulation
	    if (galdef.verbose == -701) cout<<"#Proton data "<<it->first<<":"<<endl;
	    //cout<<"Modulation:"<<modulation[it->second]<<endl;
	    for (int i = 0; i < dataEmean.size(); ++i){
	       if (flux[i] >= 0) {
		  if (flux[i] < dataFlux[i]) {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrMinus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrMinus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  } else {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrPlus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrPlus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  }
	       }
	    }
	    if (galdef.verbose == -701) cout<<endl<<endl;
	 }
         chisqa[2] = chisq - prt;
         prt = prt +  chisqa[2];         
      }
   }
   if (sets[3]) { //antiprotons
      //Data storage for the proton spectrum
      vector<double> antiprotons = readNuclei(gcr,n_species,-1,1,iz1,iz2,ir1,ir2,ix1,ix2,iy1,iy2,fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2);
      if (antiprotons.size() == 0){
	 std::cerr<<"No anti-protons found in nuclei data, not using those in chisq"<<std::endl;
      } else {
	 //Select the appropriate data
	 std::vector<double> dataEmean, dataFlux, dataErrMinus, dataErrPlus;
	 vector<int> Z_d, A_d, Z_n, A_n;
	 vector<string> references(1); //Storage to select references
	 map<string,int> refs; //The actual references and mapping to modulation
	 if (expmnts[14]) refs["Or00"] = 14;
	 if (expmnts[15]) refs["Ad10"] = 15;
	 Z_d.push_back(-1);
	 A_d.push_back(1);
	 map<string,int>::iterator it = refs.begin();
	 for ( ; it != refs.end(); ++it) {
	    vector<double> dataEmean, dataFlux, dataErrMinus, dataErrPlus;
	    references[0] = it->first;
	    retrieveGCRData(dataEmean,dataFlux,dataErrMinus,dataErrPlus,gcr_data,Z_d,A_d,Z_n,A_n,"flux",references);
	    //Correct for energy scaling
	    for (int i = 0; i < dataEmean.size(); ++i) 
	       dataEmean[i] *= energyScale[it->second]; 
	    //Interpolate and modulate the model
	    vector<double> flux = interpolateAndModulate(antiprotons, gcr[0].Ekin, dataEmean, -1, 1, modulation[it->second]);
	    //Now we can calculate the chisq, including modulation
	    if (galdef.verbose == -701) cout<<"#anti-proton data "<<it->first<<":"<<endl;
	    //cout<<"Modulation:"<<modulation[it->second]<<endl;
	    for (int i = 0; i < dataEmean.size(); ++i){
	       if (flux[i] >= 0) {
		  if (flux[i] < dataFlux[i]) {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrMinus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrMinus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  } else {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrPlus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrPlus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  }
	       }
	    }
	    if (galdef.verbose == -701) cout<<endl<<endl;
	 }
         chisqa[3] = chisq - prt;
         prt = prt +  chisqa[3];         
      }
   }
   if (sets[4]) { //Helium
      //Data storage for the Helium spectrum of the model
      vector<double> He4 = readNuclei(gcr,n_species,2,4,iz1,iz2,ir1,ir2,ix1,ix2,iy1,iy2,fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2);
      if (He4.size()) {
	 //Select the appropriate data
	 std::vector<double> dataEmean, dataFlux, dataErrMinus, dataErrPlus;
	 vector<int> Z_d, A_d, Z_n, A_n;
	 vector<string> references(1); //Storage to select references
	 map<string,int> refs; //The actual references and mapping to modulation
	 if (expmnts[2]) refs["Pa06"] = 2;
	 if (expmnts[4]) refs["Ag02"] = 4;
	 if (expmnts[7]) refs["As98"] = 7;
	 if (expmnts[13]) refs["Sa00"] = 13;
	 Z_d.push_back(2);
	 A_d.push_back(4);
	 map<string,int>::iterator it = refs.begin();
	 for ( ; it != refs.end(); ++it) {
	    vector<double> dataEmean, dataFlux, dataErrMinus, dataErrPlus;
	    references[0] = it->first;
	    retrieveGCRData(dataEmean,dataFlux,dataErrMinus,dataErrPlus,gcr_data,Z_d,A_d,Z_n,A_n,"flux",references);
	    //Correct for energy scaling
	    for (int i = 0; i < dataEmean.size(); ++i) 
	       dataEmean[i] *= energyScale[it->second]; 
	    //Interpolate and modulate the model
	    vector<double> flux = interpolateAndModulate(He4, gcr[0].Ekin, dataEmean, 2, 4, modulation[it->second]);
			for (int i = 0; i < flux.size(); ++i) flux[i] *= 4;
	    //Now we can calculate the chisq, including modulation
	    if (galdef.verbose == -701) cout<<"#He4 data "<<it->first<<":"<<endl;
	    //cout<<"Modulation:"<<modulation[it->second]<<endl;
	    for (int i = 0; i < dataEmean.size(); ++i){
	       if (flux[i] >= 0) {
		  if (flux[i] < dataFlux[i]) {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrMinus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrMinus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  } else {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrPlus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrPlus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  }
	       }
	    }
	    if (galdef.verbose == -701) cout<<endl<<endl;
	 }
         chisqa[4] = chisq - prt;
         prt = prt +  chisqa[4];         
      } else {
	 cerr<<"No He4 found in the model, not using that in chisq"<<endl;
      }
   }
   if (sets[5]) { //Oxygen
      //Data storage for the Oxygen spectrum of the model
      vector<double> O16 = readNuclei(gcr,n_species,8,16,iz1,iz2,ir1,ir2,ix1,ix2,iy1,iy2,fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2);
      vector<double> O17 = readNuclei(gcr,n_species,8,17,iz1,iz2,ir1,ir2,ix1,ix2,iy1,iy2,fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2);
      vector<double> O18 = readNuclei(gcr,n_species,8,18,iz1,iz2,ir1,ir2,ix1,ix2,iy1,iy2,fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2);
      //Need at least one of them
      if (O16.size() || O17.size() || O18.size() ) {
	 //Get the data
	 vector<int> Z_d, A_d, Z_n, A_n;
	 vector<string> references(1); //Storage to select references
	 map<string,int> refs; //The actual references and mapping to modulation
	 if (expmnts[0]) refs["En90"] = 0;
	 if (expmnts[1]) refs["ACEW"] = 1;
	 Z_d.push_back(8);
	 A_d.push_back(-1);
	 map<string,int>::iterator it = refs.begin();
	 for ( ; it != refs.end(); ++it) {
	    vector<double> dataEmean, dataFlux, dataErrMinus, dataErrPlus;
	    references[0] = it->first;
	    retrieveGCRData(dataEmean,dataFlux,dataErrMinus,dataErrPlus,gcr_data,Z_d,A_d,Z_n,A_n,"flux",references);
	    //Correct for energy scaling
	    for (int i = 0; i < dataEmean.size(); ++i) 
	       dataEmean[i] *= energyScale[it->second]; 
	    //Interpolate and modulate the available data
	    vector<double> flux(dataEmean.size(),0.0);
	    if (O16.size()) {
	       vector<double> flux16 = interpolateAndModulate(O16, gcr[0].Ekin, dataEmean, 8, 16, modulation[it->second]);
	       for (int i = 0; i < flux16.size(); ++i) flux16[i] *= 16; 
	       for (int i = 0; i < flux16.size(); ++i) flux[i] += flux16[i];
	    } else {
	       cerr<<"O16 not found in the model"<<endl;
	    }
	    if (O17.size()) {
	       vector<double> flux17 = interpolateAndModulate(O17, gcr[0].Ekin, dataEmean, 8, 17, modulation[it->second]);
	       for (int i = 0; i < flux17.size(); ++i) flux17[i] *= 17; 
	       for (int i = 0; i < flux17.size(); ++i) flux[i] += flux17[i];
	    } else {
	       cerr<<"O17 not found in the model"<<endl;
	    }
	    if (O18.size()) {
	       vector<double> flux18 = interpolateAndModulate(O18, gcr[0].Ekin, dataEmean, 8, 18, modulation[it->second]);
	       for (int i = 0; i < flux18.size(); ++i) flux18[i] *= 18; 
	       for (int i = 0; i < flux18.size(); ++i) flux[i] += flux18[i];
	    } else {
	       cerr<<"O18 not found in the model"<<endl;
	    }
	    //Now we can calculate the chisq, including modulation
	    if (galdef.verbose == -701) cout<<"#Oxygen data "<<it->first<<":"<<endl;
	    //cout<<"Modulation:"<<modulation[it->second]<<endl;
	    for (int i = 0; i < dataEmean.size(); ++i){
	       if (flux[i] >= 0) {
		  if (flux[i] < dataFlux[i]) {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrMinus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrMinus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  } else {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrPlus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrPlus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  }
	       }
	    }
	    if (galdef.verbose == -701) cout<<endl<<endl;
	 }
      	 chisqa[5] = chisq - prt;
         prt = prt + chisqa[5];
      } else {
	 cerr<<"No O found in the model, not using in chisq"<<endl;
      }
   }
   if (sets[6]) { //Carbon
      //Data storage for the Carbon spectrum of the model
      vector<double> C12 = readNuclei(gcr,n_species,6,12,iz1,iz2,ir1,ir2,ix1,ix2,iy1,iy2,fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2);
      vector<double> C13 = readNuclei(gcr,n_species,6,13,iz1,iz2,ir1,ir2,ix1,ix2,iy1,iy2,fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2);
      if (C12.size() || C13.size()) {
	 vector<int> Z_d, A_d, Z_n, A_n;
	 vector<string> references(1); //Storage to select references
	 map<string,int> refs; //The actual references and mapping to modulation
	 if (expmnts[0]) refs["En90"] = 0;
	 if (expmnts[1]) refs["ACEW"] = 1;
	 Z_d.push_back(6);
	 A_d.push_back(-1);
	 map<string,int>::iterator it = refs.begin();
	 for ( ; it != refs.end(); ++it) {
	    vector<double> dataEmean, dataFlux, dataErrMinus, dataErrPlus;
	    references[0] = it->first;
	    retrieveGCRData(dataEmean,dataFlux,dataErrMinus,dataErrPlus,gcr_data,Z_d,A_d,Z_n,A_n,"flux",references);
	    //Correct for energy scaling
	    for (int i = 0; i < dataEmean.size(); ++i) 
	       dataEmean[i] *= energyScale[it->second]; 
	    //Interpolate and modulate the available data
	    vector<double> flux(dataEmean.size(),0.0);
	    if (C12.size()) {
	       vector<double> flux12 = interpolateAndModulate(C12, gcr[0].Ekin, dataEmean, 6, 12, modulation[it->second]);
	       for (int i = 0; i < flux12.size(); ++i) flux12[i] *= 12; 
	       for (int i = 0; i < flux12.size(); ++i) flux[i] += flux12[i];
	    } else {
	       cerr<<"C12 not found in the model"<<endl;
	    }
	    if (C13.size()) {
	       vector<double> flux13 = interpolateAndModulate(C13, gcr[0].Ekin, dataEmean, 6, 13, modulation[it->second]);
	       for (int i = 0; i < flux13.size(); ++i) flux13[i] *= 13; 
	       for (int i = 0; i < flux13.size(); ++i) flux[i] += flux13[i];
	    } else {
	       cerr<<"C13 not found in the model"<<endl;
	    }
	    //Now we can calculate the chisq, including modulation
	    if (galdef.verbose == -701) cout<<"#Carbon data "<<it->first<<":"<<endl;
	    //cout<<"Modulation:"<<modulation[it->second]<<endl;
	    //cout<<C12[0]<<", "<<C13[0];
	    for (int i = 0; i < dataEmean.size(); ++i){
	       if (flux[i] >= 0) {
		  if (flux[i] < dataFlux[i]) {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrMinus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrMinus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  } else {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrPlus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrPlus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  }
	       }
	    }
	    if (galdef.verbose == -701) cout<<endl<<endl;
	 }
      	 chisqa[6] = chisq - prt;
	 prt = prt + chisqa[6];
      } else {
	 cerr<<"No C found in the model, not using in chisq"<<endl;
      }
   }
   if (sets[7]) { //B/C ratio
      //Data storage for the Carbon spectrum of the model
      vector<double> C12 = readNuclei(gcr,n_species,6,12,iz1,iz2,ir1,ir2,ix1,ix2,iy1,iy2,fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2);
      vector<double> C13 = readNuclei(gcr,n_species,6,13,iz1,iz2,ir1,ir2,ix1,ix2,iy1,iy2,fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2);
      //Data storage for the Boron spectrum of the model
      vector<double> B10 = readNuclei(gcr,n_species,5,10,iz1,iz2,ir1,ir2,ix1,ix2,iy1,iy2,fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2);
      vector<double> B11 = readNuclei(gcr,n_species,5,11,iz1,iz2,ir1,ir2,ix1,ix2,iy1,iy2,fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2);
      //We must have one of each at leas
      if ( (C12.size() || C13.size() ) && (B10.size() || B11.size())) {
	 //Get the data
	 vector<int> Z_d, A_d, Z_n, A_n;
	 vector<string> references(1); //Storage to select references
	 map<string,int> refs; //The actual references and mapping to modulation
	 if (expmnts[0]) refs["En90"] = 0;
	 if (expmnts[1]) {
	    refs["Da00"] = 1;
	    //refs["No01"] = 1;  /These are really bad datapoints, completely
			//incompatible with the rest
	 }
	 if (expmnts[2]) refs["Pa07"] = 2;
	 if (expmnts[3]) refs["Ahn8"] = 3;
	 Z_d.push_back(6);
	 A_d.push_back(12);
	 Z_d.push_back(6);
	 A_d.push_back(13);
	 Z_n.push_back(5);
	 A_n.push_back(10);
	 Z_n.push_back(5);
	 A_n.push_back(11);
	 map<string,int>::iterator it = refs.begin();
	 for ( ; it != refs.end(); ++it) {
	    vector<double> dataEmean, dataFlux, dataErrMinus, dataErrPlus;
	    references[0] = it->first;
	    retrieveGCRData(dataEmean,dataFlux,dataErrMinus,dataErrPlus,gcr_data,Z_d,A_d,Z_n,A_n,"ratio",references);
	    //Correct for energy scaling
	    for (int i = 0; i < dataEmean.size(); ++i) 
	       dataEmean[i] *= energyScale[it->second]; 
	    //Interpolate and modulate the data
	    vector<double> flux(dataEmean.size(),0.0), fluxC(dataEmean.size(), 0.0);


	    if (C12.size()){
	       vector<double> flux12 = interpolateAndModulate(C12, gcr[0].Ekin, dataEmean, 6, 12, modulation[it->second]);
	       for (int i = 0; i < flux12.size(); ++i) flux12[i] *= 12; 
	       for (int i = 0; i < flux12.size(); ++i) fluxC[i] += flux12[i];
	    } else {
	       cerr<<"C12 not found, not included in B/C chisq ratio"<<endl;
	    }
	    if (C13.size()){
	       vector<double> flux13 = interpolateAndModulate(C13, gcr[0].Ekin, dataEmean, 6, 13, modulation[it->second]);
	       for (int i = 0; i < flux13.size(); ++i) flux13[i] *= 13; 
	       for (int i = 0; i < flux13.size(); ++i) fluxC[i] += flux13[i];
	    } else {
	       cerr<<"C13 not found, not included in B/C chisq ratio"<<endl;
	    }
	    if (B10.size()){
	       vector<double> flux10 = interpolateAndModulate(B10, gcr[0].Ekin, dataEmean, 5, 10, modulation[it->second]);
	       for (int i = 0; i < flux10.size(); ++i) flux10[i] *= 10; 
	       for (int i = 0; i < flux10.size(); ++i) flux[i] += flux10[i];
	    } else {
	       cerr<<"B10 not found, not included in B/C chisq ratio"<<endl;
	    }
	    if (B11.size()){
	       vector<double> flux11 = interpolateAndModulate(B11, gcr[0].Ekin, dataEmean, 5, 11, modulation[it->second]);
	       for (int i = 0; i < flux11.size(); ++i) flux11[i] *= 11; 
	       for (int i = 0; i < flux11.size(); ++i) flux[i] += flux11[i];
	    } else {
	       cerr<<"B11 not found, not included in B/C chisq ratio"<<endl;
	    }
	    //Calculate the ratio and chisq
	    for (int i = 0; i < flux.size(); ++i) {
	       if (fluxC[i] > 0 && flux[i] >= 0) {
		  flux[i] /= fluxC[i];
	       }else {
		  flux[i] = -1;
	       }
	    }
	    if (galdef.verbose == -701) cout<<"#B/C ratio "<<it->first<<":"<<endl;
	    for (int i = 0; i < dataEmean.size(); ++i){
	       if (flux[i] >= 0) {
		 if (flux[i] < dataFlux[i]) {
		   chisq += pow((flux[i]-dataFlux[i])/dataErrMinus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrMinus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  } else {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrPlus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrPlus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  }
	       }
	    }
	    if (galdef.verbose == -701) cout<<endl<<endl;
	 }
         chisqa[7] = chisq - prt;
         prt = prt + chisqa[7];
      } else {
	 cerr<<"We must have both B and C to calculate B/C ratio, not included in chisq"<<endl;
      }
   }
   if (sets[8]) { //Be10/Be9 ratio
      //Data storage for the Berillium spectrum of the model
      vector<double> Be10 = readNuclei(gcr,n_species,4,10,iz1,iz2,ir1,ir2,ix1,ix2,iy1,iy2,fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2);
      vector<double> Be9 = readNuclei(gcr,n_species,4,9,iz1,iz2,ir1,ir2,ix1,ix2,iy1,iy2,fz1,fz2,fr1,fr2,fx1,fx2,fy1,fy2);
      //We must have both
      if ( Be10.size() && Be9.size() ) {
	 //Get the data
	 vector<int> Z_d, A_d, Z_n, A_n;
	 vector<string> references(1); //Storage to select references
	 map<string,int> refs; //The actual references and mapping to modulation
	 if (expmnts[1]) refs["Ya01"] = 1;
	 Z_d.push_back(4);
	 A_d.push_back(9);
	 Z_n.push_back(4);
	 A_n.push_back(10);
	 map<string,int>::iterator it = refs.begin();
	 for ( ; it != refs.end(); ++it) {
	    vector<double> dataEmean, dataFlux, dataErrMinus, dataErrPlus;
	    references[0] = it->first;
	    retrieveGCRData(dataEmean,dataFlux,dataErrMinus,dataErrPlus,gcr_data,Z_d,A_d,Z_n,A_n,"ratio",references);
	    //std::cout<<"Reading first "<<dataEmean.size()<<" points"<<std::endl;
	    //Correct for energy scaling
	    for (int i = 0; i < dataEmean.size(); ++i) 
	       dataEmean[i] *= energyScale[it->second]; 
	    //Interpolate and modulate the data
	    vector<double> flux = interpolateAndModulate(Be10,gcr[0].Ekin, dataEmean, 4, 10, modulation[it->second]);
	    vector<double> fluxD = interpolateAndModulate(Be9,gcr[0].Ekin, dataEmean, 4, 9, modulation[it->second]);
	    for (int i = 0; i < dataEmean.size(); ++i){
	       flux[i] *= 10;
	       fluxD[i] *= 9;
	    }
	    //Calculate the ratio and chisq
	    for (int i = 0; i < flux.size(); ++i) {
	       if(fluxD[i] > 0 && flux[i] >= 0) {
		  flux[i] /= fluxD[i];
	       } else {
		  flux[i] == -1;
	       }
	    }
	    if (galdef.verbose == -701) cout<<"#BE 10/9 ratio "<<it->first<<":"<<endl;
	    for (int i = 0; i < dataEmean.size(); ++i){
	       if (flux[i] >= 0) {

		  if (flux[i] < dataFlux[i]) {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrMinus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrMinus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  } else {
		    chisq += pow((flux[i]-dataFlux[i])/dataErrPlus[i], 2)*pow(10,tau[it->second]) - log(pow(10,tau[it->second]));
		     if (galdef.verbose == -701) cout<<dataEmean[i]<<" "<<flux[i]<<" "<<dataFlux[i]<<" "<<(dataErrPlus[i]/pow(10,tau[it->second]))<<" "<<chisq-prt<<endl;
		  }
	       }
	    }
	    if (galdef.verbose == -701) cout<<endl<<endl;
	 }
         chisqa[8] = chisq - prt;
      } else {
	 cerr<<"We must have both Be10 and Be9 to calculate Be10/Be9 ratio, not included in chisq"<<endl;
      }
   }



   return chisq;

}
