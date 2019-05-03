
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Galdef.cc *                                   galprop package * 10/12/2003
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>

using namespace std;

#include <ErrorLogger.h>

#include"Galdef.h"

int Bufferlength=1000;//IMOS20080114

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

Galdef::Galdef() :
  B_field_parameters ( 0 ),
  cr_source_x ( 0 ), cr_source_y ( 0 ), cr_source_z ( 0 ),
  cr_source_w ( 0 ), cr_source_L ( 0 ),
  isotopic_abundance ( 0 ),
  use_Z ( 0 ),
  ISRF_healpixOrder ( 0 ),
  deleteMemory ( false ),
  n_X_CO_values ( 0 ) {
  
}

Galdef::~Galdef() {

    //delete[] ISRF_factors;
  if ( deleteMemory ) {
    
    delete[] B_field_parameters;
    delete[] cr_source_x;
    delete[] cr_source_y;
    delete[] cr_source_z;
    delete[] cr_source_w;
    delete[] cr_source_L;
    
    if ( isotopic_abundance != 0 ) {
      
      for ( int i = 0; i < max_Z+1; ++i ) {
	
	delete [] isotopic_abundance[i];
	
      }
      
      delete [] isotopic_abundance;
      
    }
    
    delete[] use_Z;
    
  }

}

void AssignDefaultParameters ( Galdef& g );

int Galdef::read ( const string& theVersion,
                   const string& runNumber,
                   const string& galdefDirectory ) {

    INFO ( "Entry" );

    int stat=0;
    char workstring[Bufferlength];//IMOS20080114
    char parstring [Bufferlength];//IMOS20080114

    //cout<<" >>>>galdef read"<<endl;

    //strcpy(version,version_);
    //strcpy(run_no,run_no_);

    version = theVersion;
    run_no = runNumber;

    AssignDefaultParameters ( *this );

    ostringstream bufID;
    bufID << version << "_" << run_no;
    //string tmpStringBuf = bufID.str();
    //tmpStringBuf.resize(99);
    //strcpy(galdef_ID, tmpStringBuf.c_str());
    galdef_ID = bufID.str();

    ostringstream galdef_file;
    galdef_file << galdefDirectory << "galdef_" << galdef_ID;

    //char galdef_file[100];
    //strcpy(galdef_file,galdef_directory);
    //strcat(galdef_file,"galdef_");
    //strcat(galdef_file,galdef_ID);

    ostringstream oBuf;
    oBuf << "Reading from " << galdef_file.str();
    INFO ( oBuf.str() );

    //cout<<galdef_file<<endl;

    //char  filename[100];
    //strcpy(filename,galdef_file);

    const string filename = galdef_file.str();

    FILE* ft = fopen ( filename.c_str(),"r" );

    if ( 0 == ft ) {//ft==NULL) {

        ostringstream buf;
        buf << "Galdef file " << filename << " not found";
        INFO ( buf.str() );
        return -1;
    }

    fclose ( ft );

    //Start by reading and setting the verbosity
    //-2 means Errors and Fatal only
    //-1 adds Warnings
    //0 adds info
    //everything else adds debug (requires compilation with debug enabled).
    //Special values as before 
    strcpy ( parstring,                              "verbose" );
    stat= read_galdef_parameter ( filename,parstring,&verbose );
    if ( verbose == -2 ) {
       INFO("Setting logger to errors only");
       utl::ErrorLogger::GetInstance().SetSeverity(utl::ErrorLogger::eError);
    } else if ( verbose == -1 ){
       INFO("Setting logger to errors and warnings only");
       utl::ErrorLogger::GetInstance().SetSeverity(utl::ErrorLogger::eWarning);
    } else if ( verbose == 0 ) {
       INFO("Setting logger to normal (errors, warnings, and info)");
    } else {
       INFO("Setting logger to debug.  Requires compilation with --enable-debug.");
       utl::ErrorLogger::GetInstance().SetSeverity(utl::ErrorLogger::eDebug);
       utl::ErrorLogger::GetInstance().SetVerbosity(utl::ErrorLogger::eVerbose);
    }


    strcpy ( parstring,                              "n_spatial_dimensions" );
    stat= read_galdef_parameter ( filename,parstring,&n_spatial_dimensions );

    if ( stat!=0 ) {
       oBuf.str("");
       oBuf<<"no spatial dimensions specified in "<<filename; 
       FATAL(oBuf.str());
       return -1;
    }

    electron_source_normalization =1.; 
    strcpy ( parstring,                              "electron_source_norm" );
    stat= read_galdef_parameter ( filename,parstring,&electron_source_normalization );
    source_normalization =1.;          
    strcpy ( parstring,                              "source_norm" );
    stat= read_galdef_parameter ( filename,parstring,&source_normalization );

    //radial grid
    strcpy ( parstring,                              "r_min" );
    stat= read_galdef_parameter ( filename,parstring,&r_min );
    strcpy ( parstring,                              "r_max" );
    stat= read_galdef_parameter ( filename,parstring,&r_max );
    strcpy ( parstring,                              "dr" );
    stat= read_galdef_parameter ( filename,parstring,&dr );

    //z-grid
    strcpy ( parstring,                              "z_min" );
    stat= read_galdef_parameter ( filename,parstring,&z_min );
    strcpy ( parstring,                              "z_max" );
    stat= read_galdef_parameter ( filename,parstring,&z_max );
    strcpy ( parstring,                              "dz" );
    stat= read_galdef_parameter ( filename,parstring,&dz );

    //x-grid
    strcpy ( parstring,                              "x_min" );
    stat= read_galdef_parameter ( filename,parstring,&x_min );
    strcpy ( parstring,                              "x_max" );
    stat= read_galdef_parameter ( filename,parstring,&x_max );
    strcpy ( parstring,                              "dx" );
    stat= read_galdef_parameter ( filename,parstring,&dx );

    //y-grid
    strcpy ( parstring,                              "y_min" );
    stat= read_galdef_parameter ( filename,parstring,&y_min );
    strcpy ( parstring,                              "y_max" );
    stat= read_galdef_parameter ( filename,parstring,&y_max );
    strcpy ( parstring,                              "dy" );
    stat= read_galdef_parameter ( filename,parstring,&dy );

    //momentum grid
    strcpy ( parstring,                              "p_min" );
    stat= read_galdef_parameter ( filename,parstring,&p_min );
    strcpy ( parstring,                              "p_max" );
    stat= read_galdef_parameter ( filename,parstring,&p_max );
    strcpy ( parstring,                              "p_factor" );
    stat= read_galdef_parameter ( filename,parstring,&p_factor );

    //kinetic energy grid
    strcpy ( parstring,                              "Ekin_min" );
    stat= read_galdef_parameter ( filename,parstring,&Ekin_min );
    strcpy ( parstring,                              "Ekin_max" );
    stat= read_galdef_parameter ( filename,parstring,&Ekin_max );
    strcpy ( parstring,                              "Ekin_factor" );
    stat= read_galdef_parameter ( filename,parstring,&Ekin_factor );

    //p||Ekin option
    strcpy ( parstring,                              "p_Ekin_grid" );
    stat= read_galdef_parameter ( filename,parstring, p_Ekin_grid );

    //gamma-ray energy grid
    strcpy ( parstring,                              "E_gamma_min" );
    stat= read_galdef_parameter ( filename,parstring,&E_gamma_min );
    strcpy ( parstring,                              "E_gamma_max" );
    stat= read_galdef_parameter ( filename,parstring,&E_gamma_max );
    strcpy ( parstring,                              "E_gamma_factor" );
    stat= read_galdef_parameter ( filename,parstring,&E_gamma_factor );

    //mode of integration over particle spectrum
    strcpy ( parstring,                              "integration_mode" );
    stat= read_galdef_parameter ( filename,parstring,&integration_mode );

    //synchrotron grid
    strcpy ( parstring,                              "nu_synch_min" );
    stat= read_galdef_parameter ( filename,parstring,&nu_synch_min );
    strcpy ( parstring,                              "nu_synch_max" );
    stat= read_galdef_parameter ( filename,parstring,&nu_synch_max );
    strcpy ( parstring,                              "nu_synch_factor" );
    stat= read_galdef_parameter ( filename,parstring,&nu_synch_factor );

    //longitude-latitude grid
    strcpy ( parstring,                              "long_min" );
    stat= read_galdef_parameter ( filename,parstring,&long_min );
    strcpy ( parstring,                              "long_max" );
    stat= read_galdef_parameter ( filename,parstring,&long_max );
    strcpy ( parstring,                              "lat_min" );
    stat= read_galdef_parameter ( filename,parstring,&lat_min );
    strcpy ( parstring,                              "lat_max" );
    stat= read_galdef_parameter ( filename,parstring,&lat_max );

    strcpy ( parstring,                              "d_long" );
    stat= read_galdef_parameter ( filename,parstring,&d_long );
    strcpy ( parstring,                              "d_lat" );
    stat= read_galdef_parameter ( filename,parstring,&d_lat );

    //Healpix order
    strcpy ( parstring,                              "healpix_order" );
    stat= read_galdef_parameter ( filename,parstring,&healpix_order );

    //line of sight integration IMOS20080114
    strcpy ( parstring,                              "lat_substep_number" );
    stat= read_galdef_parameter ( filename,parstring,&lat_substep_number );
    strcpy ( parstring,                              "LoS_step" );
    stat= read_galdef_parameter ( filename,parstring,&LoS_step );
    strcpy ( parstring,                              "LoS_substep_number" );
    stat= read_galdef_parameter ( filename,parstring,&LoS_substep_number );

    //space diffusion
    strcpy ( parstring,                              "D0_xx" );
    stat= read_galdef_parameter ( filename,parstring,&D0_xx );

    // diffusion coefficient parametrization
    strcpy ( parstring,                              "D_rigid_br" );
    stat= read_galdef_parameter ( filename,parstring,&D_rigid_br );
    strcpy ( parstring,                              "D_g_1" );
    stat= read_galdef_parameter ( filename,parstring,&D_g_1 );
    strcpy ( parstring,                              "D_g_2" );
    stat= read_galdef_parameter ( filename,parstring,&D_g_2 );

    //diffusive reacceleration: Alfven speed
    strcpy ( parstring,                              "v_Alfven" );
    stat= read_galdef_parameter ( filename,parstring,&v_Alfven );
    strcpy ( parstring,                              "diff_reacc" );
    stat= read_galdef_parameter ( filename,parstring,&diff_reacc );
    //damping
    strcpy ( parstring,                              "damping_p0" );       //IMOS20030129
    stat= read_galdef_parameter ( filename,parstring,&damping_p0 );        //IMOS20030129
    strcpy ( parstring,                              "damping_max_path_L" );  //IMOS20030129
    stat= read_galdef_parameter ( filename,parstring,&damping_max_path_L );  //IMOS20030129
    strcpy ( parstring,                              "damping_const_G" );  //IMOS20030129
    stat= read_galdef_parameter ( filename,parstring,&damping_const_G );   //IMOS20030129

    // convection
    strcpy ( parstring,                              "convection" );     //AWS20010323
    stat= read_galdef_parameter ( filename,parstring,&convection );      //AWS20010323
    strcpy ( parstring,                              "v0_conv" );        //AWS20010323
    stat= read_galdef_parameter ( filename,parstring,&v0_conv );         //AWS20010323
    strcpy ( parstring,                              "dvdz_conv" );      //AWS20010323
    stat= read_galdef_parameter ( filename,parstring,&dvdz_conv );       //AWS20010323

    //injection spectra limitations
    rigid_min = 0;
    strcpy ( parstring,                              "rigid_min" );
    stat= read_galdef_parameter ( filename,parstring,&rigid_min );
    rigid_max = 1e100; //Some very large number that should not get in the way if rigid_max is not defined
    strcpy ( parstring,                              "rigid_max" );
    stat= read_galdef_parameter ( filename,parstring,&rigid_max );

    //injection spectra parametrization (nucleons)
    strcpy ( parstring,                              "nuc_rigid_br" );
    stat= read_galdef_parameter ( filename,parstring,&nuc_rigid_br );
    strcpy ( parstring,                              "nuc_g_1" );
    stat= read_galdef_parameter ( filename,parstring,&nuc_g_1 );
    strcpy ( parstring,                              "nuc_g_2" );
    stat= read_galdef_parameter ( filename,parstring,&nuc_g_2 );

    //rigidity||beta_rig||Etot  option
    strcpy ( parstring,                              "inj_spectrum_type" ); // IMOS20000613.2
    stat= read_galdef_parameter ( filename,parstring, inj_spectrum_type ); // IMOS20000613.3

    //injection spectra parametrization (electrons)
    strcpy ( parstring,                              "electron_g_0" );   // IMOS20031012
    stat= read_galdef_parameter ( filename,parstring,&electron_g_0 );
    strcpy ( parstring,                              "electron_rigid_br0" );
    stat= read_galdef_parameter ( filename,parstring,&electron_rigid_br0 );
    strcpy ( parstring,                              "electron_g_1" );
    stat= read_galdef_parameter ( filename,parstring,&electron_g_1 );
    strcpy ( parstring,                              "electron_rigid_br" );
    stat= read_galdef_parameter ( filename,parstring,&electron_rigid_br );
    strcpy ( parstring,                              "electron_g_2" );
    stat= read_galdef_parameter ( filename,parstring,&electron_g_2 );

    //other parameters
    strcpy ( parstring,                              "He_H_ratio" );
    stat= read_galdef_parameter ( filename,parstring,&He_H_ratio );


    //X_CO
    strcpy ( parstring,                              "n_X_CO" ); //IMOS20080114
    stat= read_galdef_parameter ( filename,parstring,&n_X_CO ); //IMOS20080114
    strcpy ( parstring,                              "X_CO" );
    stat= read_galdef_parameter ( filename,parstring,&X_CO );        //IMOS20080114
    strcpy ( parstring,                              "X_CO_parameters_0" );
    stat= read_galdef_parameter ( filename,parstring,&X_CO_parameters[0] );
    strcpy ( parstring,                              "X_CO_parameters_1" );
    stat= read_galdef_parameter ( filename,parstring,&X_CO_parameters[1] );
    strcpy ( parstring,                              "X_CO_parameters_2" );
    stat= read_galdef_parameter ( filename,parstring,&X_CO_parameters[2] );
    strcpy ( parstring,                              "X_CO_parameters_3" );
    stat= read_galdef_parameter ( filename,parstring,&X_CO_parameters[3] );

    //X_CO parameters for n_X_CO 1, using tabulated values
    strcpy ( parstring,                              "n_X_CO_values" );
    stat= read_galdef_parameter ( filename,parstring,&n_X_CO_values );
    if ( stat!=0 )
        n_X_CO_values = 0;
    X_CO_values = new double[n_X_CO_values];
    strcpy ( parstring,                              "X_CO_values" );
    stat= read_galdef_parameter ( filename,parstring,X_CO_values,n_X_CO_values );
    X_CO_radius = new double[n_X_CO_values];
    strcpy ( parstring,                              "X_CO_radius" );
    stat= read_galdef_parameter ( filename,parstring,X_CO_radius,n_X_CO_values );

    propagation_X_CO = 0;                    // default                   //AWS20090623
    strcpy ( parstring,                              "propagation_X_CO" ); //AWS20090623
    stat= read_galdef_parameter ( filename,parstring,&propagation_X_CO ); //AWS20090623

    nHI_model        = 1;                    // default                   //AWS20090814
    strcpy ( parstring,                              "nHI_model" );       //AWS20090814
    stat= read_galdef_parameter ( filename,parstring,&nHI_model );        //AWS20090814

    nH2_model        = 1;                    // default                   //AWS20090814
    strcpy ( parstring,                              "nH2_model" );       //AWS20090814
    stat= read_galdef_parameter ( filename,parstring,&nH2_model );        //AWS20090814

    nHII_model        = 1;                   // default                   //AWS20090814
    strcpy ( parstring,                              "nHII_model" );      //AWS20090814
    stat= read_galdef_parameter ( filename,parstring,&nHII_model );       //AWS20090814


    //CO gas file
    strcpy ( parstring,                              "COR_filename" );               //IMOS20080114
    stat= read_galdef_parameter ( filename,parstring, COR_filename );                //IMOS20080114
    if ( stat!=0 ) exit ( -1 );                                                      //IMOS20080114

    //HI gas file
    strcpy ( parstring,                              "HIR_filename" );               //IMOS20080114
    stat= read_galdef_parameter ( filename,parstring, HIR_filename );                //IMOS20080114
    if ( stat!=0 ) exit ( -1 );                                                      //IMOS20080114

    //CR-database, defaults to GCR_data_1.dat for backwards compatibility
    strcpy ( GCR_data_filename,                      "GCR_data_1.dat" );             //GJ20100103
    strcpy ( parstring,                              "GCR_data_filename" );          //GJ20100103
    stat= read_galdef_parameter ( filename,parstring, GCR_data_filename );           //GJ20100103

    strcpy ( parstring,                              "fragmentation" );
    stat= read_galdef_parameter ( filename,parstring,&fragmentation );
    strcpy ( parstring,                              "momentum_losses" );
    stat= read_galdef_parameter ( filename,parstring,&momentum_losses );
    strcpy ( parstring,                              "radioactive_decay" );
    stat= read_galdef_parameter ( filename,parstring,&radioactive_decay );
    strcpy ( parstring,                              "K_capture" );       //AWS20010731
    stat= read_galdef_parameter ( filename,parstring,&K_capture );        //AWS20010731
    strcpy ( parstring,                              "ionization_rate" ); // IMOS20060420
    stat= read_galdef_parameter ( filename,parstring,&ionization_rate );


    //time grid
    strcpy ( parstring,                              "start_timestep" );
    stat= read_galdef_parameter ( filename,parstring,&start_timestep );
    strcpy ( parstring,                              "end_timestep" );
    stat= read_galdef_parameter ( filename,parstring,&end_timestep );
    strcpy ( parstring,                              "timestep_factor" );
    stat= read_galdef_parameter ( filename,parstring,&timestep_factor );
    strcpy ( parstring,                              "timestep_repeat" );
    stat= read_galdef_parameter ( filename,parstring,&timestep_repeat );
    strcpy ( parstring,                              "timestep_repeat2" );
    stat= read_galdef_parameter ( filename,parstring,&timestep_repeat2 );

    strcpy ( parstring,                              "timestep_print" );
    stat= read_galdef_parameter ( filename,parstring,&timestep_print );
    strcpy ( parstring,                              "timestep_diagnostics" );
    stat= read_galdef_parameter ( filename,parstring,&timestep_diagnostics );
    strcpy ( parstring,                              "control_diagnostics" );
    stat= read_galdef_parameter ( filename,parstring,&control_diagnostics );

    //other parameters
    strcpy ( parstring,                              "network_iterations" );
    stat= read_galdef_parameter ( filename,parstring,&network_iterations );
    strcpy ( parstring,                              "network_iter_compl" );
    stat= read_galdef_parameter ( filename,parstring,&network_iter_compl );
    if ( stat!=0 )
        network_iter_compl = min ( network_iterations,2 );
    strcpy ( parstring,                              "network_iter_sec" );
    stat= read_galdef_parameter ( filename,parstring,&network_iter_sec );
    if ( stat!=0 )
        network_iter_sec = min ( network_iterations,1 );

    strcpy ( parstring,                              "prop_r" );
    stat= read_galdef_parameter ( filename,parstring,&prop_r );
    strcpy ( parstring,                              "prop_x" );
    stat= read_galdef_parameter ( filename,parstring,&prop_x );
    strcpy ( parstring,                              "prop_y" );
    stat= read_galdef_parameter ( filename,parstring,&prop_y );
    strcpy ( parstring,                              "prop_z" );
    stat= read_galdef_parameter ( filename,parstring,&prop_z );

    strcpy ( parstring,                              "prop_p" );
    stat= read_galdef_parameter ( filename,parstring,&prop_p );

    strcpy ( parstring,                              "use_symmetry" );
    stat= read_galdef_parameter ( filename,parstring,&use_symmetry );

    strcpy ( parstring,                              "vectorized" );
    stat= read_galdef_parameter ( filename,parstring,&vectorized );


    strcpy ( parstring,                              "source_specification" );
    stat= read_galdef_parameter ( filename,parstring,&source_specification );

    //nuclei to include & abundances
    strcpy ( parstring,                              "max_Z" );
    stat= read_galdef_parameter ( filename,parstring,&max_Z );

    use_Z = new int[max_Z+1];
    for ( int iZ=1; iZ<=max_Z; iZ++ ) {
        sprintf ( parstring,"use_Z_%d",iZ );
        stat= read_galdef_parameter ( filename,parstring,&use_Z[iZ] );
    }

    isotopic_abundance = new double*[max_Z+1];
    for ( int iZ=0; iZ<=max_Z; iZ++ ) {
        isotopic_abundance[iZ]=new double[max_Z*3];
        for ( int iA=0; iA<max_Z*3; isotopic_abundance[iZ][iA++]=0.0 );
    }
    for ( int iZ=1; iZ<=max_Z; iZ++ )
        if ( use_Z[iZ]==1 )
	  for ( int iA=2*iZ-2; iA<2.5*iZ+4.2; iA++) {
                sprintf ( parstring,"iso_abundance_%02d_%03d",iZ,iA ); // eg _03_007
                DEBUGLOG(std::string(parstring));

                stat= read_galdef_parameter ( filename,parstring,&isotopic_abundance[iZ][iA], true );
            }

    strcpy ( parstring,                              "total_cross_section" );//AWS20010620
    stat= read_galdef_parameter ( filename,parstring,&total_cross_section );//AWS20010620


    strcpy ( parstring,                              "cross_section_option" );
    stat= read_galdef_parameter ( filename,parstring,&cross_section_option );

    strcpy ( parstring,                              "t_half_limit" );    //AWS20010731
    stat= read_galdef_parameter ( filename,parstring,&t_half_limit );     //AWS20010731
    
    // DM positron and antiproton
    strcpy(parstring,                              "primary_DM_positron"    ); //BiXJ 2007/2/1
    stat= read_galdef_parameter(filename,parstring,&primary_DM_positron     ); //BiXJ 2007/2/1
    strcpy(parstring,                              "primary_DM_electron"    ); //BiXJ 2007/2/1
    stat= read_galdef_parameter(filename,parstring,&primary_DM_electron     ); //BiXJ 2007/2/1
    strcpy(parstring,                              "primary_DM_antip"       ); //BiXJ 2007/2/1
    stat= read_galdef_parameter(filename,parstring,&primary_DM_antip        ); //BiXJ 2007/2/1
    strcpy ( parstring,                              "DM_positron_filename" ); //HM 2014/4/24   
    stat= read_galdef_parameter ( filename,parstring, DM_positron_filename );  //HM 2014/4/24           
    strcpy ( parstring,                              "DM_electron_filename" ); //HM 2014/4/24   
    stat= read_galdef_parameter ( filename,parstring, DM_electron_filename );  //HM 2014/4/24           
    strcpy ( parstring,                              "DM_antip_filename" ); //HM 2014/4/24   
    stat= read_galdef_parameter ( filename,parstring, DM_antip_filename );  //HM 2014/4/24           

    //CR species to include
    strcpy ( parstring,                              "primary_electrons" );
    stat= read_galdef_parameter ( filename,parstring,&primary_electrons );
    strcpy ( parstring,                              "secondary_positrons" );
    stat= read_galdef_parameter ( filename,parstring,&secondary_positrons );
    strcpy ( parstring,                              "secondary_electrons" );
    stat= read_galdef_parameter ( filename,parstring,&secondary_electrons );
    strcpy ( parstring,                              "knock_on_electrons" ); //IMOS20060504
    stat= read_galdef_parameter ( filename,parstring,&knock_on_electrons ); //IMOS20060504
    strcpy ( parstring,                              "secondary_antiproton" );
    stat= read_galdef_parameter ( filename,parstring,&secondary_antiprotons );
    strcpy ( parstring,                              "tertiary_antiproton" );  // IMOS20000605
    stat= read_galdef_parameter ( filename,parstring,&tertiary_antiprotons );
    strcpy ( parstring,                              "secondary_protons" );  // IMOS20000605
    stat= read_galdef_parameter ( filename,parstring,&secondary_protons );

    strcpy ( parstring,                              "gamma_rays" );
    stat= read_galdef_parameter ( filename,parstring,&gamma_rays );

    strcpy ( parstring,                              "pi0_decay" );         // IMOS20050216
    stat= read_galdef_parameter ( filename,parstring,&pi0_decay );
    if ( stat!=0 ) exit ( -1 );                                               //AWS20050301

    strcpy ( parstring,                              "IC_isotropic" );     // IMOS20060420
    stat= read_galdef_parameter ( filename,parstring,&IC_isotropic );

    strcpy ( parstring,                              "IC_anisotropic" );
    stat= read_galdef_parameter ( filename,parstring,&IC_anisotropic );

    strcpy ( parstring,                              "bremss" );       // IMOS20060420
    stat= read_galdef_parameter ( filename,parstring,&bremss );

    strcpy ( parstring,                              "synchrotron" );
    stat= read_galdef_parameter ( filename,parstring,&synchrotron );

    //CR source parameters

    strcpy ( parstring,                              "source_model" );
    stat = read_galdef_parameter ( filename,parstring,&source_model );

    //If not defined, electron source model is the same as nuclei
    source_model_electron = source_model;
    strcpy ( parstring,                              "source_model_elec" );
    stat = read_galdef_parameter ( filename,parstring,&source_model_electron );

    //Resize the source_prameters vectors
    source_parameters.resize(10);
    source_parameters_electron.resize(10);
    for (int iS = 0; iS < 10; ++iS) {

      ostringstream buf1, buf2;
      buf1 << "source_parameters_" << iS;
      buf2 << "source_pars_elec_" << iS;

      source_parameters[iS] = source_parameters_electron[iS] = (4 == iS ? 100 : 0);

      strcpy ( parstring, buf1.str().c_str() );
      stat = read_galdef_parameter ( filename,parstring,&source_parameters[iS] );

      source_parameters_electron[iS] = source_parameters[iS]; 

      strcpy ( parstring, buf2.str().c_str() );
      stat = read_galdef_parameter ( filename,parstring,&source_parameters_electron[iS] );

    }
    
    /*source_parameters[0] = 0;
    strcpy ( parstring,                              "source_parameters_0" );
    stat = read_galdef_parameter ( filename,parstring,&source_parameters[0] );

    source_parameters[1] = 0;
    strcpy ( parstring,                              "source_parameters_1" );
    stat = read_galdef_parameter ( filename,parstring,&source_parameters[1] );
    
    source_parameters[2] = 0;
    strcpy ( parstring,                              "source_parameters_2" );
    stat = read_galdef_parameter ( filename,parstring,&source_parameters[2] );
    
    source_parameters[3] = 0;
    strcpy ( parstring,                              "source_parameters_3" );
    stat = read_galdef_parameter ( filename,parstring,&source_parameters[3] );

    source_parameters[4] = 100.;
    strcpy ( parstring,                              "source_parameters_4" );
    stat = read_galdef_parameter ( filename,parstring,&source_parameters[4] );

    source_parameters[5] = 0;
    strcpy ( parstring,                              "source_parameters_5" );
    stat = read_galdef_parameter ( filename,parstring,&source_parameters[5] );

    source_parameters[6] = 0;
    strcpy ( parstring,                              "source_parameters_6" );
    stat = read_galdef_parameter ( filename,parstring,&source_parameters[6] );

    source_parameters[7] = 0;
    strcpy ( parstring,                              "source_parameters_7" );
    stat = read_galdef_parameter ( filename,parstring,&source_parameters[7] );

    source_parameters[8] = 0;
    strcpy ( parstring,                              "source_parameters_8" );
    stat = read_galdef_parameter ( filename,parstring,&source_parameters[8] );

    source_parameters[9] = 0;
    strcpy ( parstring,                              "source_parameters_9" ); 
    stat = read_galdef_parameter ( filename,parstring,&source_parameters[9] ); 
    */

    ostringstream srcBuf;

    srcBuf << "Source model nuclei: " << source_model << ", electrons: " << source_model_electron;
    INFO(srcBuf.str());
    
    for (int iS = 0; iS < 10; ++iS) {

       srcBuf.str("");
      srcBuf << iS << " " << source_parameters[iS] << " " << source_parameters_electron[iS];
    INFO(srcBuf.str());

    }

    //INFO(srcBuf.str());

    //CR source parameters for source model 8, using tabulated values
    strcpy ( parstring,                              "n_source_values" );
    stat = read_galdef_parameter ( filename,parstring,&n_source_values );
    if ( stat!=0 )
        n_source_values = 0;
    source_values = new double[n_source_values];
    strcpy ( parstring,                              "source_values" );
    stat= read_galdef_parameter ( filename,parstring,source_values,n_source_values );
    source_radius = new double[n_source_values];
    strcpy ( parstring,                              "source_radius" );
    stat= read_galdef_parameter ( filename,parstring,source_radius,n_source_values );


    strcpy ( parstring,                              "n_cr_sources" );
    stat= read_galdef_parameter ( filename,parstring,&n_cr_sources );
    cr_source_x         =new double[n_cr_sources];
    cr_source_y         =new double[n_cr_sources];
    cr_source_z         =new double[n_cr_sources];
    cr_source_w         =new double[n_cr_sources];
    cr_source_L         =new double[n_cr_sources];

    for ( int i_cr_source=0; i_cr_source<n_cr_sources; i_cr_source++ ) {
        sprintf ( parstring,"cr_source_x_%02d", i_cr_source+1 ); // eg _01
        stat= read_galdef_parameter ( filename,parstring,&cr_source_x[i_cr_source] );

        sprintf ( parstring,"cr_source_y_%02d", i_cr_source+1 );
        stat= read_galdef_parameter ( filename,parstring,&cr_source_y[i_cr_source] );

        sprintf ( parstring,"cr_source_z_%02d", i_cr_source+1 );
        stat= read_galdef_parameter ( filename,parstring,&cr_source_z[i_cr_source] );

        sprintf ( parstring,"cr_source_w_%02d", i_cr_source+1 );
        stat= read_galdef_parameter ( filename,parstring,&cr_source_w[i_cr_source] );

        sprintf ( parstring,"cr_source_L_%02d", i_cr_source+1 );
        stat= read_galdef_parameter ( filename,parstring,&cr_source_L[i_cr_source] );
    }

//SNR
    strcpy ( parstring,                              "SNR_events" );
    stat= read_galdef_parameter ( filename,parstring,&SNR_events );
    strcpy ( parstring,                              "SNR_interval" );
    stat= read_galdef_parameter ( filename,parstring,&SNR_interval );
    strcpy ( parstring,                              "SNR_livetime" );
    stat= read_galdef_parameter ( filename,parstring,&SNR_livetime );

    strcpy ( parstring,                              "SNR_electron_sdg" );       //AWS20010410
    stat= read_galdef_parameter ( filename,parstring,&SNR_electron_sdg );        //AWS20010410
    strcpy ( parstring,                              "SNR_nuc_sdg" );            //AWS20010410
    stat= read_galdef_parameter ( filename,parstring,&SNR_nuc_sdg );             //AWS20010410

    strcpy ( parstring,                              "SNR_electron_dgpivot" );   //AWS20010410
    stat= read_galdef_parameter ( filename,parstring,&SNR_electron_dgpivot );    //AWS20010410
    strcpy ( parstring,                              "SNR_nuc_dgpivot" );        //AWS20010410
    stat= read_galdef_parameter ( filename,parstring,&SNR_nuc_dgpivot );         //AWS20010410

    //HI, CO  surveys                                                               //AWS20051309
    //HI_survey=8;                                                     //default
    //strcpy(parstring,                              "HI_survey"           );
    //stat= read_galdef_parameter(filename,parstring,&HI_survey            );
    //CO_survey=8;                                                     //default
    //strcpy(parstring,                              "CO_survey"           );
    //stat= read_galdef_parameter(filename,parstring,&CO_survey            );

    //magnetic field
    strcpy ( parstring,                              "B_field_model" );     // original definition of exponential B(R,z)
    stat= read_galdef_parameter ( filename,parstring,&B_field_model );


    strcpy ( B_field_name,"galprop_original" );                                  // default if absent in galdef file             //AWS20080312
    strcpy ( parstring,                              "B_field_name" );                                                           //AWS20080312
    stat= read_galdef_parameter ( filename,parstring, B_field_name );                                                            //AWS20080312

    n_B_field_parameters=10; // default if absent in galdef file                                                                 //AWS20080313
    strcpy ( parstring,                              "n_B_field_parameters" );                                                   //AWS20080313
    stat= read_galdef_parameter ( filename,parstring,&n_B_field_parameters );                                                    //AWS20080313

    if ( n_B_field_parameters!=10 ) {FATAL("need exactly 10 B-field parameters in this version!"); exit ( 1 );}             //AWS20080313

    B_field_parameters=new double[n_B_field_parameters];                                                                         //AWS20080312
    strcpy ( parstring,                              "B_field_parameters" );                                                     //AWS20080312
    stat = read_galdef_parameter ( filename,parstring,B_field_parameters,n_B_field_parameters );
    if ( stat!=0 )  //Default parameters if not found
        for ( int i = 0; i < n_B_field_parameters; ++i ) {
            B_field_parameters[i] = -i-1;
        }
    //strcpy(workstring,"-1,-2,-3,-4,-5,-6,-7,-8,-9,-10");                          //    default if absent in galdef file         //AWS20080312
    //stat= read_galdef_parameter(filename,parstring, workstring           );                                                      //AWS20080312
    //// this is not flexible, need to find a better way to read in arrays                                                         //AWS20080313
    //sscanf(workstring,"%le,%le,%le,%le,%le,%le,%le,%le,%le,%le",                                                                 //AWS20080312
//	 &B_field_parameters[0],&B_field_parameters[1],&B_field_parameters[2],&B_field_parameters[3],&B_field_parameters[4],   //AWS20080312
//	 &B_field_parameters[5],&B_field_parameters[6],&B_field_parameters[7],&B_field_parameters[8],&B_field_parameters[9])  ;//AWS20080312


    //ISRF file
    strcpy ( parstring, "ISRF_file" ); //AWS20050301
    stat = read_galdef_parameter ( filename, parstring, ISRF_file ); //AWS20050301

    if ( stat!=0 ) exit ( -1 ); //AWS20050301

    //ISRF file type
    strcpy ( parstring, "ISRF_filetype" ); //TAP20072401
    stat = read_galdef_parameter ( filename, parstring, &ISRF_filetype ); //TAP20072401
    if ( stat != 0 )
        exit ( -1 ); //TAP20072401

    // ISRF healpix order for skymaps
    strcpy ( parstring, "ISRF_healpixOrder" );
    stat = read_galdef_parameter ( filename, parstring, &ISRF_healpixOrder );

    if ( stat != 0 )
        exit ( -1 );

    //ISRF factors for IC calculation
    strcpy ( parstring, "ISRF_factors" );             //AWS20050301
    stat = read_galdef_parameter ( filename, parstring, workstring ); //AWS20050301

    if ( stat != 0 )
        exit ( -1 ); //AWS20050301

    //ISRF_factors = new double[3];                                                            //AWS20050301
    sscanf ( workstring,"%le,%le,%le",&ISRF_factors[0],&ISRF_factors[1],&ISRF_factors[2] );   //AWS20050301

    //spectra normalization
    strcpy ( parstring,                              "proton_norm_Ekin" );
    stat= read_galdef_parameter ( filename,parstring,&proton_norm_Ekin );
    strcpy ( parstring,                              "proton_norm_flux" );
    stat= read_galdef_parameter ( filename,parstring,&proton_norm_flux );

    strcpy ( parstring,                              "electron_norm_Ekin" );
    stat= read_galdef_parameter ( filename,parstring,&electron_norm_Ekin );
    strcpy ( parstring,                              "electron_norm_flux" );
    stat= read_galdef_parameter ( filename,parstring,&electron_norm_flux );

    //dark matter (DM) parameters  IMOS20050912
    strcpy ( parstring,                              "DM_positrons" );
    stat= read_galdef_parameter ( filename,parstring,&DM_positrons );
    strcpy ( parstring,                              "DM_electrons" );
    stat= read_galdef_parameter ( filename,parstring,&DM_electrons );
    strcpy ( parstring,                              "DM_antiprotons" );
    stat= read_galdef_parameter ( filename,parstring,&DM_antiprotons );
    strcpy ( parstring,                              "DM_gammas" );
    stat= read_galdef_parameter ( filename,parstring,&DM_gammas );

    strcpy ( parstring,                              "DM_double0" );
    stat= read_galdef_parameter ( filename,parstring,&DM_double0 );
    strcpy ( parstring,                              "DM_double1" );
    stat= read_galdef_parameter ( filename,parstring,&DM_double1 );
    strcpy ( parstring,                              "DM_double2" );
    stat= read_galdef_parameter ( filename,parstring,&DM_double2 );
    strcpy ( parstring,                              "DM_double3" );
    stat= read_galdef_parameter ( filename,parstring,&DM_double3 );
    strcpy ( parstring,                              "DM_double4" );
    stat= read_galdef_parameter ( filename,parstring,&DM_double4 );
    strcpy ( parstring,                              "DM_double5" );
    stat= read_galdef_parameter ( filename,parstring,&DM_double5 );
    strcpy ( parstring,                              "DM_double6" );
    stat= read_galdef_parameter ( filename,parstring,&DM_double6 );
    strcpy ( parstring,                              "DM_double7" );
    stat= read_galdef_parameter ( filename,parstring,&DM_double7 );
    strcpy ( parstring,                              "DM_double8" );
    stat= read_galdef_parameter ( filename,parstring,&DM_double8 );
    strcpy ( parstring,                              "DM_double9" );
    stat= read_galdef_parameter ( filename,parstring,&DM_double9 );

    strcpy ( parstring,                              "DM_int0" );
    stat= read_galdef_parameter ( filename,parstring,&DM_int0 );
    strcpy ( parstring,                              "DM_int1" );
    stat= read_galdef_parameter ( filename,parstring,&DM_int1 );
    strcpy ( parstring,                              "DM_int2" );
    stat= read_galdef_parameter ( filename,parstring,&DM_int2 );
    strcpy ( parstring,                              "DM_int3" );
    stat= read_galdef_parameter ( filename,parstring,&DM_int3 );
    strcpy ( parstring,                              "DM_int4" );
    stat= read_galdef_parameter ( filename,parstring,&DM_int4 );
    strcpy ( parstring,                              "DM_int5" );
    stat= read_galdef_parameter ( filename,parstring,&DM_int5 );
    strcpy ( parstring,                              "DM_int6" );
    stat= read_galdef_parameter ( filename,parstring,&DM_int6 );
    strcpy ( parstring,                              "DM_int7" );
    stat= read_galdef_parameter ( filename,parstring,&DM_int7 );
    strcpy ( parstring,                              "DM_int8" );
    stat= read_galdef_parameter ( filename,parstring,&DM_int8 );
    strcpy ( parstring,                              "DM_int9" );
    stat= read_galdef_parameter ( filename,parstring,&DM_int9 );

    //output controls
    strcpy ( parstring,                              "skymap_format" );
    stat= read_galdef_parameter ( filename,parstring,&skymap_format );
    strcpy ( parstring,                              "output_gcr_full" );
    stat= read_galdef_parameter ( filename,parstring,&output_gcr_full );


    strcpy ( parstring,                              "warm_start" );    //AWS20010121
    stat= read_galdef_parameter ( filename,parstring,&warm_start );     //AWS20010121

    strcpy ( parstring,                              "test_suite" );
    stat= read_galdef_parameter ( filename,parstring,&test_suite );

    print();

    deleteMemory=true;
    INFO ( "Exit" );
    //cout<<" <<<<galdef read"<<endl;
    return stat;

}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Galdef::
read_galdef_parameter ( const std::string& filename,
                        char *parstring,
                        double *value,
                        int size ) {
   if (size == 0) return 0;
    FILE *ft = fopen ( filename.c_str(),"r" );
    if ( !ft )
        return -1; //==NULL) return -1;
    char *flag;
    int   found;

    char dumstring[Bufferlength];//IMOS20080114
    char input    [Bufferlength];//IMOS20080114
    int stat=0;

    found=-1;
    while ( found!=0 ) {

        flag=fgets ( input,Bufferlength,ft ); // read string until newline (Schildt p.222)
        if ( flag==0 ) {

	   ostringstream buf;
            buf<<parstring<<" not found in galdef file!";
	    WARNING(buf);
            fclose ( ft );
            stat=1;
            return stat;
        }
        //  printf("%s",input);       // string is \0 terminated
        sscanf ( input,"%s"  ,dumstring );
        found=strcmp ( dumstring,parstring ); // search for parstring  in input

        if ( found==0 ) {
            int start = 22;
            int inc;
            for ( int i = 1; i < size; ++i ) {
                sscanf ( &input[start], "%le,%n", &value[i-1], &inc );
                start += inc;
            }
            sscanf ( &input[start], "%le", &value[size-1] );
            //  cout<<parstring<<"    found!"<<endl;
            //  cout<<endl<<"difread: value of "<< parstring <<"="<<*value <<endl;
        }
    }

    fclose ( ft );

    //  cout<<" <<<< read_galdef_parameter"<<endl;
    return stat;
}

int Galdef::
read_galdef_parameter ( const std::string& filename,
                        char *parstring,
                        double *value ) {

    FILE *ft = fopen ( filename.c_str(),"r" );
    if ( !ft )
        return -1; //==NULL) return -1;
    char *flag;
    int   found;

    char dumstring[Bufferlength];//IMOS20080114
    char input    [Bufferlength];//IMOS20080114
    int stat=0;

    found=-1;
    while ( found!=0 ) {

        flag=fgets ( input,Bufferlength,ft ); // read string until newline (Schildt p.222)
        if ( flag==0 ) {

	   ostringstream buf;
            buf<<parstring<<" not found in galdef file!";
	    WARNING(buf);
            fclose ( ft );
            stat=1;
            return stat;
        }
        //  printf("%s",input);       // string is \0 terminated
        sscanf ( input,"%s"  ,dumstring );
        found=strcmp ( dumstring,parstring ); // search for parstring  in input

        if ( found==0 ) {
            sscanf ( input,"%22c%le"  ,dumstring, value );  // le for double
            //  cout<<parstring<<"    found!"<<endl;
            //  cout<<endl<<"difread: value of "<< parstring <<"="<<*value <<endl;
        }

	//cout << parstring << " " << found << " " << dumstring << " " << *value << endl;

    }

    fclose ( ft );

    //  cout<<" <<<< read_galdef_parameter"<<endl;
    return stat;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Galdef:: // AV 2010
read_galdef_parameter ( const std::string& filename,
                        char *parstring,
                        double *value, bool mute ) {

    FILE *ft = fopen ( filename.c_str(),"r" );
    if ( !ft )
        return -1; //==NULL) return -1;
    char *flag;
    int   found;

    char dumstring[Bufferlength];//IMOS20080114
    char input    [Bufferlength];//IMOS20080114
    int stat=0;

    found=-1;
    while ( found!=0 ) {

        flag=fgets ( input,Bufferlength,ft ); // read string until newline (Schildt p.222)
        if ( flag==0 ) {

	  if (!mute)
	    {
	   ostringstream buf;
            buf<<parstring<<" not found in galdef file!";
	    WARNING(buf);
            fclose ( ft );
            stat=1;
	    }
            return stat;
        }
        //  printf("%s",input);       // string is \0 terminated
        sscanf ( input,"%s"  ,dumstring );
        found=strcmp ( dumstring,parstring ); // search for parstring  in input

        if ( found==0 ) {
            sscanf ( input,"%22c%le"  ,dumstring, value );  // le for double
            //  cout<<parstring<<"    found!"<<endl;
            //  cout<<endl<<"difread: value of "<< parstring <<"="<<*value <<endl;
        }

	//cout << parstring << " " << found << " " << dumstring << " " << *value << endl;

    }

    fclose ( ft );

    //  cout<<" <<<< read_galdef_parameter"<<endl;
    return stat;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Galdef::
read_galdef_parameter ( const std::string& filename,
                        char *parstring,
                        int *value ) {

    FILE *ft =fopen ( filename.c_str(), "r" );
    if ( !ft )
        return -1;
    char *flag;
    int   found;
    char dumstring[Bufferlength];//IMOS20080114
    char input    [Bufferlength];//IMOS20080114
    int stat=0;

    found=-1;
    while ( found!=0 ) {

        flag=fgets ( input,Bufferlength,ft );   // read string until newline (Schildt p.222)
        if ( flag==0 ) {
	   ostringstream buf;
            buf<<parstring<<" not found in galdef file!";
	    WARNING(buf);
            fclose ( ft );
            stat=1;
            return stat;
        }
        //  printf("%s",input);       // string is \0 terminated
        sscanf ( input,"%s"  ,dumstring );
        found=strcmp ( dumstring,parstring ); // search for parstring  in input

        if ( found==0 ) {
            sscanf ( input,"%22c%d"   ,dumstring, value );  // d for int
            //  cout<<parstring<<"    found!"<<endl;
            //  cout<<endl<<"difread: value of "<< parstring <<"="<<*value <<endl;
        }
    }
    fclose ( ft );
    //  cout<<" <<<< read_galdef_parameter"<<endl;
    return stat;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Galdef::
read_galdef_parameter ( const std::string& filename,
                        char *parstring,
                        char *value ) {

    FILE *ft =fopen ( filename.c_str(),"r" );
    if ( !ft )
        return -1;
    char *flag;
    int   found;
    char dumstring[Bufferlength];//IMOS20080114
    char input    [Bufferlength];//IMOS20080114
    int stat=0;

    found=-1;
    while ( found!=0 ) {

        flag=fgets ( input,Bufferlength,ft );   // read string until newline (Schildt p.222)
        if ( flag==0 ) {

	   ostringstream buf;
            buf<<parstring<<" not found in galdef file!";
	    WARNING(buf);
            fclose ( ft );
            stat=1;
            return stat;
        }
        //  printf("%s",input);       // string is \0 terminated
        sscanf ( input,"%s"  ,dumstring );
        found=strcmp ( dumstring,parstring ); // search for parstring  in input

        if ( found==0 ) {
            sscanf ( input,"%22c%s"   ,dumstring, value );    // s for string
            // cout<<parstring<<"    found!"<<endl;
            // cout<<endl<<"difread: value of "<< parstring <<"="<<*value <<endl;
        }
    }
    fclose ( ft );
    //  cout<<" <<<< read_galdef_parameter"<<endl;
    return stat;
}

void Galdef::print() {

    cout<<"  ======= galdef: "<<galdef_ID  <<endl;
    cout<<"  version  "<<version <<endl;
    cout<<"  run_no   "<<run_no      <<endl;

    cout<<"  n_spatial_dimensions "<<n_spatial_dimensions   <<endl;

    cout<<"  r_min    "<<r_min   <<endl;
    cout<<"  r_max    "<<r_max   <<endl;
    cout<<"  dr       "<<dr      <<endl;
    cout<<"  x_min    "<<x_min   <<endl;
    cout<<"  x_max    "<<x_max   <<endl;
    cout<<"  dx       "<<dx      <<endl;
    cout<<"  y_min    "<<y_min   <<endl;
    cout<<"  y_max    "<<y_max   <<endl;
    cout<<"  dy       "<<dy      <<endl;
    cout<<"  z_min    "<<z_min   <<endl;
    cout<<"  z_max    "<<z_max   <<endl;
    cout<<"  dz       "<<dz      <<endl;
    cout<<"  p_min    "<<p_min   <<endl;
    cout<<"  p_max    "<<p_max   <<endl;
    cout<<"  p_factor "<<p_factor<<endl;
    cout<<"  Ekin_min "<<   Ekin_min      <<endl;
    cout<<"  Ekin_max "<<   Ekin_max      <<endl;
    cout<<"  Ekin_factor "<<Ekin_factor   <<endl;

    cout<<"  p_Ekin_grid "<<p_Ekin_grid   <<endl;

    cout<<"  E_gamma_min   "<<   E_gamma_min      <<endl;
    cout<<"  E_gamma_max   "<<   E_gamma_max      <<endl;
    cout<<"  E_gamma_factor "<<   E_gamma_factor   <<endl;

    cout<<"  nu_synch_min    "<<  nu_synch_min      <<endl;
    cout<<"  nu_synch_max    "<<  nu_synch_max      <<endl;
    cout<<"  nu_synch_factor "<<  nu_synch_factor   <<endl;

    cout<<"  long_min      "<<   long_min         <<endl;
    cout<<"  long_max      "<<   long_max         <<endl;
    cout<<"  lat_min       "<<    lat_min         <<endl;
    cout<<"  lat_max       "<<    lat_max         <<endl;
    cout<<"  d_long        "<<   d_long           <<endl;
    cout<<"  d_lat         "<<   d_lat            <<endl;
    cout<<"  healpix_order "<<   healpix_order    <<endl;
    cout<<"  lat_substep_number "<< lat_substep_number <<endl;
    cout<<"  LoS_step      "<<      LoS_step           <<endl;
    cout<<"  LoS_substep_number "<< LoS_substep_number <<endl;

    cout<<"  D0_xx       "<<D0_xx         <<endl;
    cout<<"  D_rigid_br  "<<D_rigid_br    <<endl;
    cout<<"  D_g_1       "<<D_g_1         <<endl;
    cout<<"  D_g_2       "<<D_g_2         <<endl;
    cout<<"  diff_reacc  "<<diff_reacc    <<endl;
    cout<<"  v_Alfven    "<<v_Alfven      <<endl;
    cout<<"  convection  "<<convection    <<endl; //AWS20010323
    cout<<"  v0_conv     "<<v0_conv       <<endl; //AWS20010323
    cout<<"  dvdz_conv   "<<dvdz_conv     <<endl; //AWS20010323

    cout<<"  damping_p0         "<<damping_p0        <<endl;//IMOS20060330
    cout<<"  damping_const_G    "<<damping_const_G   <<endl;//IMOS20060330
    cout<<"  damping_max_path_L "<<damping_max_path_L<<endl;//IMOS20060330

    cout<<"  nuc_rigid_br  "<<nuc_rigid_br    <<endl;
    cout<<"  nuc_g_1       "<<nuc_g_1         <<endl;
    cout<<"  nuc_g_2       "<<nuc_g_2         <<endl;

    cout<<"  rigid_min     "<<rigid_min       <<endl;
    cout<<"  rigid_max     "<<rigid_max       <<endl;

    cout<<"  inj_spectrum_type "<<inj_spectrum_type<<endl; // IMOS20000613.4

    cout<<"  electron_g_0       "<<electron_g_0         <<endl;// IMOS20031012
    cout<<"  electron_rigid_br0 "<<electron_rigid_br0   <<endl;
    cout<<"  electron_g_1       "<<electron_g_1         <<endl;
    cout<<"  electron_rigid_br  "<<electron_rigid_br    <<endl;
    cout<<"  electron_g_2       "<<electron_g_2         <<endl;

    cout<<"  He_H_ratio        "<<He_H_ratio          <<endl;
    cout<<"  n_X_CO            "<<n_X_CO              <<endl;//IMOS20080114
    cout<<"  X_CO              "<<X_CO                <<endl;//IMOS20080114
    cout<<"  X_CO_parameters_0 "<<X_CO_parameters[0]  <<endl;
    cout<<"  X_CO_parameters_1 "<<X_CO_parameters[1]  <<endl;
    cout<<"  X_CO_parameters_2 "<<X_CO_parameters[2]  <<endl;
    cout<<"  X_CO_parameters_3 "<<X_CO_parameters[3]  <<endl;
    cout<<"  n_X_CO_values     "<<n_X_CO_values       <<endl;
    if ( n_X_CO_values > 0 ) {
        cout<<"   X_CO_values      ";
        for ( int i = 0; i < n_X_CO_values; ++i )
            cout<<X_CO_values[i]<<", ";
        cout<<endl;
        cout<<"   X_CO_radius      ";
        for ( int i = 0; i < n_X_CO_values; ++i )
            cout<<X_CO_radius[i]<<", ";
        cout<<endl;
    }
    cout<<"  propagation_X_CO  "<<propagation_X_CO    <<endl;//AWS20090623
    cout<<"  nHI_model         "<<nHI_model           <<endl;//AWS20090814
    cout<<"  nH2_model         "<<nH2_model           <<endl;//AWS20090814
    cout<<"  nHII_model        "<<nHII_model          <<endl;//AWS20090814

    cout<<"  COR_filename      "<<COR_filename        <<endl;//IMOS20080114
    cout<<"  HIR_filename      "<<HIR_filename        <<endl;//IMOS20080114

    cout<<"  GCR_data_filename "<<GCR_data_filename   <<endl;

    cout<<"  fragmentation     "<<fragmentation       <<endl;
    cout<<"  momentum_losses   "<<momentum_losses     <<endl;
    cout<<"  radioactive_decay "<<radioactive_decay   <<endl;
    cout<<"  K_capture         "<<K_capture           <<endl;
    cout<<"  ionization_rate   "<<ionization_rate     <<endl;  // IMOS20060420

    cout<<"  start_timestep   "<< start_timestep<<endl;
    cout<<"  end_timestep     "<<   end_timestep<<endl;
    cout<<"  timestep_factor  "<<timestep_factor<<endl;
    cout<<"  timestep_repeat  "<<timestep_repeat<<endl;
    cout<<"  timestep_repeat2 "<<timestep_repeat2<<endl;
    cout<<"  timestep_print   "<<timestep_print <<endl;
    cout<<"  timestep_diagnostics  "<<timestep_diagnostics <<endl;
    cout<<"  control_diagnostics  " << control_diagnostics <<endl;

    cout<<"  network_iterations "<<network_iterations<<endl;
    cout<<"  network_iter_compl "<<network_iter_compl<<endl;
    cout<<"  network_iter_sec   "<<network_iter_sec<<endl;

    cout<<"  prop_r   "<<prop_r  <<endl;
    cout<<"  prop_x   "<<prop_x  <<endl;
    cout<<"  prop_y   "<<prop_y  <<endl;
    cout<<"  prop_z   "<<prop_z  <<endl;
    cout<<"  prop_p   "<<prop_p  <<endl;

    cout<<"  use_symmetry  "<<use_symmetry<<endl;
    cout<<"  vectorized    "<<vectorized  <<endl;

    cout<<"  source_specification  "<<source_specification <<endl;
    cout<<"  source_normalization  "<<source_normalization <<endl;
    cout<<"  source_model          "<<source_model         <<endl;
    cout<<"  source_parameters_1   "<<source_parameters[1] <<endl;
    cout<<"  source_parameters_2   "<<source_parameters[2] <<endl;
    cout<<"  source_parameters_3   "<<source_parameters[3] <<endl;
    cout<<"  source_parameters_4   "<<source_parameters[4] <<endl;
    cout<<"  source_parameters_5   "<<source_parameters[5] <<endl;
    cout<<"  source_parameters_6   "<<source_parameters[6] <<endl;
    cout<<"  source_parameters_7   "<<source_parameters[7] <<endl;
    cout<<"  source_parameters_8   "<<source_parameters[8] <<endl;
    cout<<"  source_parameters_9   "<<source_parameters[9] <<endl;
    cout<<"  source_parameters_6   ="<<source_parameters[6] <<endl;
    cout<<"  source_parameters_7   ="<<source_parameters[7] <<endl;
    cout<<"  source_parameters_8   ="<<source_parameters[8] <<endl;
    cout<<"  source_parameters_9   ="<<source_parameters[9] <<endl;
    cout<<"  source_parameters_10  ="<<source_parameters[10]<<endl;

    cout<<"  source_model_electron           "<<source_model_electron          <<endl; //AWS20100507
    cout<<"  electron_source_normalization   "<<electron_source_normalization  <<endl;
    cout<<"  source_parameters_electron_1    "<<source_parameters_electron[1]  <<endl; //AWS20100507
    cout<<"  source_parameters_electron_2    "<<source_parameters_electron[2]  <<endl; //AWS20100507
    cout<<"  source_parameters_electron_3    "<<source_parameters_electron[3]  <<endl; //AWS20100507
    cout<<"  source_parameters_electron_4    "<<source_parameters_electron[4]  <<endl; //AWS20100507
    cout<<"  source_parameters_electron_5    "<<source_parameters_electron[5]  <<endl; //AWS20100507
    cout<<"  source_parameters_electron_6    "<<source_parameters_electron[6]  <<endl; //AWS20100507
    cout<<"  source_parameters_electron_7    "<<source_parameters_electron[7]  <<endl; //AWS20100507
    cout<<"  source_parameters_electron_8    "<<source_parameters_electron[8]  <<endl; //AWS20100507
    cout<<"  source_parameters_electron_9    "<<source_parameters_electron[9]  <<endl; //AWS20100507
    cout<<"  source_parameters_electron_6    ="<<source_parameters_electron[6]  <<endl; //AWS20100507
    cout<<"  source_parameters_electron_7    ="<<source_parameters_electron[7]  <<endl; //AWS20100507
    cout<<"  source_parameters_electron_8    ="<<source_parameters_electron[8]  <<endl; //AWS20100507
    cout<<"  source_parameters_electron_9    ="<<source_parameters_electron[9]  <<endl; //AWS20100507
    cout<<"  source_parameters_electron_10   ="<<source_parameters_electron[10] <<endl; //AWS20100507

    cout<<"  n_source_values       "<<n_source_values      <<endl;
    if ( n_source_values > 0 ) {
        cout<<"   source_values        ";
        for ( int i = 0; i < n_source_values; ++i )
            cout<<source_values[i]<<", ";
        cout<<endl;
        cout<<"   source_radius        ";
        for ( int i = 0; i < n_source_values; ++i )
            cout<<source_radius[i]<<", ";
        cout<<endl;
    }


    cout<<"  n_cr_sources          "<<n_cr_sources         <<endl;
    for ( int i_cr_source=0; i_cr_source<n_cr_sources; i_cr_source++ ) {
        cout<<"  cr_source_x["<<i_cr_source<<"]  "<<   cr_source_x[i_cr_source]<<endl;
        cout<<"  cr_source_y["<<i_cr_source<<"]  "<<   cr_source_y[i_cr_source]<<endl;
        cout<<"  cr_source_z["<<i_cr_source<<"]  "<<   cr_source_z[i_cr_source]<<endl;
        cout<<"  cr_source_w["<<i_cr_source<<"]  "<<   cr_source_w[i_cr_source]<<endl;
        cout<<"  cr_source_L["<<i_cr_source<<"]  "<<   cr_source_L[i_cr_source]<<endl;
    }

    cout<<"  SNR_events            "<<SNR_events           <<endl;
    cout<<"  SNR_interval          "<<SNR_interval         <<endl;
    cout<<"  SNR_livetime          "<<SNR_livetime         <<endl;

    cout<<"  SNR_electron_sdg      "<<SNR_electron_sdg     <<endl;
    cout<<"  SNR_nuc_sdg           "<<SNR_nuc_sdg          <<endl;
    cout<<"  SNR_electron_dgpivot  "<<SNR_electron_dgpivot <<endl;
    cout<<"  SNR_nuc_dgpivot       "<<SNR_nuc_dgpivot      <<endl;

    //   cout<<"  HI_survey             "<<HI_survey            <<endl;
    //   cout<<"  CO_survey             "<<CO_survey            <<endl;

    cout<<"  B_field_model         "<<B_field_model        <<endl;                                                              //AWS20080313
    cout<<"  B_field_name          "<<B_field_name         <<endl;                                                              //AWS20080313
    cout<<"n_B_field_parameters    "<<n_B_field_parameters <<endl;                                                              //AWS20080313
    cout<<"  B_field_parameters    ";for ( int i=0;i<n_B_field_parameters;i++ ) {cout<< B_field_parameters[i]<<"  ";}
    cout<<endl;   //AWS20080313

    cout<<"  ISRF_file             "<<ISRF_file            <<endl;
    cout<<"  ISRF_filetype         "<<ISRF_filetype        <<endl;                                                              //AWS20091013
    cout<<"  ISRF_healpixOrder     "<<ISRF_healpixOrder    <<endl;                                                              //AWS20091013
    cout<<"  ISRF_factors          "<<ISRF_factors[0]<<","<<ISRF_factors[1]<<","<<ISRF_factors[2]   <<endl;

    if ( ISRF_filetype == 1 ) {FATAL(" ERROR: invalid ISRF_filetype, see README for details"); exit ( 0 );}                        //AWS20091013

    cout<<"    proton_norm_Ekin    "<<  proton_norm_Ekin   <<endl;
    cout<<"    proton_norm_flux    "<<  proton_norm_flux   <<endl;
    cout<<"  electron_norm_Ekin    "<<electron_norm_Ekin   <<endl;
    cout<<"  electron_norm_flux    "<<electron_norm_flux   <<endl;

    cout<<"  max_Z                 "<<max_Z                <<endl;
    for ( int iZ=1; iZ<=max_Z; iZ++ ) cout<<"use_Z_"<<iZ<<" "<<use_Z[iZ]<<endl;

    for ( int iZ=1; iZ<=max_Z; iZ++ )
        for ( int iA=iZ; iA<max_Z*3; iA++ ) {
            if ( isotopic_abundance[iZ][iA]!=0.0 )
                cout<<"isotopic abundance for (Z,A) ("
                    <<iZ<<","<<iA<<") ="<<isotopic_abundance[iZ][iA]<<endl ;
        }

    cout<<"  total_cross_section  " <<total_cross_section  <<endl;  //AWS20010620

    cout<<"  cross_section_option  "<<cross_section_option <<endl;
    cout<<"  t_half_limit          "<<t_half_limit         <<endl;  //AWS20010731


    cout<<"  primary_electrons     "<<primary_electrons    <<endl;
    cout<<"  secondary_positrons   "<<secondary_positrons  <<endl;
    cout<<"  secondary_electrons   "<<secondary_electrons  <<endl;
    cout<<"  knock_on_electrons    "<<knock_on_electrons   <<endl;  //IMOS20060504
    cout<<"  secondary_antiprotons "<<secondary_antiprotons<<endl;
    cout<<"  tertiary_antiprotons  "<<tertiary_antiprotons <<endl;  // IMOS20000605.5
    cout<<"  secondary_protons     "<<secondary_protons    <<endl;  // IMOS20000605.6

    cout<<"  gamma_rays            "<<gamma_rays           <<endl;
    cout<<"  pi0_decay             "<<pi0_decay            <<endl;  // AWS20050218
    cout<<"  IC_isotropic          "<<IC_isotropic         <<endl;  // IMOS20060420
    cout<<"  IC_anisotropic        "<<IC_anisotropic       <<endl;
    cout<<"  bremss                "<<bremss               <<endl;  // IMOS20060420
    cout<<"  synchrotron           "<<synchrotron          <<endl;

// DM: IMOS20050912
    cout<<"  DM_positrons          "<<DM_positrons         <<endl;
    cout<<"  DM_electrons          "<<DM_electrons         <<endl;
    cout<<"  DM_antiprotons        "<<DM_antiprotons       <<endl;
    cout<<"  DM_gammas             "<<DM_gammas            <<endl;

    cout<<"  DM_double0            "<<DM_double0           <<endl;
    cout<<"  DM_double1            "<<DM_double1           <<endl;
    cout<<"  DM_double2            "<<DM_double2           <<endl;
    cout<<"  DM_double3            "<<DM_double3           <<endl;
    cout<<"  DM_double4            "<<DM_double4           <<endl;
    cout<<"  DM_double5            "<<DM_double5           <<endl;
    cout<<"  DM_double6            "<<DM_double6           <<endl;
    cout<<"  DM_double7            "<<DM_double7           <<endl;
    cout<<"  DM_double8            "<<DM_double8           <<endl;
    cout<<"  DM_double9            "<<DM_double9           <<endl;

    cout<<"  DM_int0               "<<DM_int0              <<endl;
    cout<<"  DM_int1               "<<DM_int1              <<endl;
    cout<<"  DM_int2               "<<DM_int2              <<endl;
    cout<<"  DM_int3               "<<DM_int3              <<endl;
    cout<<"  DM_int4               "<<DM_int4              <<endl;
    cout<<"  DM_int5               "<<DM_int5              <<endl;
    cout<<"  DM_int6               "<<DM_int6              <<endl;
    cout<<"  DM_int7               "<<DM_int7              <<endl;
    cout<<"  DM_int8               "<<DM_int8              <<endl;
    cout<<"  DM_int9               "<<DM_int9              <<endl;

    cout<<"  skymap_format        "<<skymap_format         <<endl;
    cout<<"  output_gcr_full      "<<output_gcr_full       <<endl;
    cout<<"  warm_start           "<<warm_start            <<endl;

    cout<<"  verbose    "<<verbose    <<endl;
    cout<<"  test_suite "<<test_suite <<endl;

    // DM annihilaton BiXJ 2007/2/2
    cout<<"  primary_DM_positron   "<<primary_DM_positron  <<endl;  //BiXJ 2007/2/2
    cout<<"  primary_DM_electron   "<<primary_DM_electron  <<endl;  //BiXJ 2007/2/2
    cout<<"  primary_DM_antiproton "<<primary_DM_antip     <<endl;  //BiXJ 2007/2/2
    cout<<"  DM_positron_filename  "<<DM_positron_filename <<endl;  //HM 2014/4/24    
    cout<<"  DM_electron_filename  "<<DM_electron_filename <<endl;  //HM 2014/4/24
    cout<<"  DM_antip_filename     "<<DM_antip_filename <<endl;  //HM 2014/4/24
    //   exit(0);
}

void AssignDefaultParameters ( Galdef& g ) {

    g.n_spatial_dimensions = 2;

    g.r_min =   0.0;
    g.r_max =  20.0;
    g.dr    =   1.0;

    g.x_min = -15.0;
    g.x_max = +15.0;
    g.dx    =   1.0;

    g.y_min = -15.0;
    g.y_max = +15.0;
    g.dy    =   1.0;

    g.z_min =  -4.0;
    g.z_max =  +4.0;
    g.dz    =   0.2;

    g.p_min =   10.0;
    g.p_max =   1.e5;
    g.p_factor= 1.3 ;

    g.Ekin_min             =1.0e1;
    g.Ekin_max             =1.0e8;
    g.Ekin_factor          =1.4;

    strcpy ( g.p_Ekin_grid,"Ekin" );

    g.E_gamma_min          = 1.e1;
    g.E_gamma_max          = 1.e6;
    g.E_gamma_factor       = 1.e1;
    g.integration_mode     = 1;

    g.nu_synch_min         = 1.0e6;
    g.nu_synch_max         = 1.0e10;
    g.nu_synch_factor      = 2.0;

    g.long_min             =   .0;
    g.long_max             =360.0;
    g.lat_min              =-90.0;
    g.lat_max              =+90.0;
    g.d_long               =  0.5 ;
    g.d_lat                =  0.5 ;
    g.lat_substep_number   =  1;
    g.LoS_step             =  0.01;
    g.LoS_substep_number   =  1;

    g.D0_xx                =5.75e28;
    g.D_rigid_br           =4.0e3;
    g.D_g_1                = 0.34;
    g.D_g_2                = 0.34;
    g.diff_reacc           =1 ;
    g.v_Alfven             =36.;

    g.damping_p0           = 1.e6;
    g.damping_const_G      = 0.02;
    g.damping_max_path_L   = 3.e21;

    g.convection           =0 ;
    g.v0_conv              =0.;
    g.dvdz_conv            =10.;

    g.nuc_rigid_br         =9.0e3;
    g.nuc_g_1              =1.82;
    g.nuc_g_2              =2.36;

    g.rigid_min            =0;
    g.rigid_max            =1e100;
    strcpy ( g.inj_spectrum_type,"rigidity" );

    g.electron_g_0         =1.60;
    g.electron_rigid_br0   =4.0e3;
    g.electron_g_1         =2.50;
    g.electron_rigid_br    =1.0e9;
    g.electron_g_2         =5.0;

    g.He_H_ratio           = 0.11;
    g.n_X_CO               = 0;
    g.X_CO                 = 1.9E20;
    g.n_X_CO_values        = 0;
    g.X_CO_parameters[0]   = 1.9E20;
    g.X_CO_parameters[1]   = 0;
    g.X_CO_parameters[2]   = 0;
    g.X_CO_parameters[3]   = 0;

    strcpy ( g.COR_filename," " );
    strcpy ( g.HIR_filename," " );

    g.B_field_model        = 050100020;

    strcpy ( g.ISRF_file," " );

    g.ISRF_filetype        = 0;
//  g.ISRF_factors[0]      = 1.0;
//  g.ISRF_factors[1]      = 1.0;
//  g.ISRF_factors[2]      = 1.0;

    g.fragmentation        =1;
    g.momentum_losses      =1;
    g.radioactive_decay    =1;
    g.K_capture            =1;
    g.ionization_rate      =0;

    g.start_timestep       =  1.e9;
    g.end_timestep         =  1.e2;
    g.timestep_factor      =  0.25;
    g.timestep_repeat      = 20;
    g.timestep_repeat2     = 0;
    g.timestep_print       =10000;
    g.timestep_diagnostics =10000;
    g.control_diagnostics  =0;

    g.network_iterations   = 1;
    g.network_iter_compl   = 1;
    g.network_iter_sec     = 1;

    g.prop_r               = 1;
    g.prop_x               = 1;
    g.prop_y               = 1;
    g.prop_z               = 1;
    g.prop_p               = 1;

    g.use_symmetry         = 0;

    g.vectorized           = 0;

    g.source_specification = 0;

    g.source_model         = 0;
    g.source_model_electron = 0;

    g.source_parameters.resize(10);
    g.source_parameters_electron.resize(10);
    for (int iS = 0; iS < 10; ++iS) {
    
      g.source_parameters[iS] = g.source_parameters_electron[iS] = 0;

    }

    //g.source_parameters[1] = 0.5;
    //g.source_parameters[2] = 1.0;
    //g.source_parameters[3] = 20.0;

    g.n_cr_sources         = 0;

    g.SNR_events           = 0;

    g.proton_norm_Ekin     = 1.00e+5;
    g.proton_norm_flux     = 4.90e-9;

    g.electron_norm_Ekin   = 34.5e3;
    g.electron_norm_flux   = .40e-9;

    g.max_Z                = 2;
//  g.use_Z[0]             = 1;
//  g.use_Z[1]             = 1;

//  g.isotopic_abundance[1][1]= 1.06e+06;
//  g.isotopic_abundance[1][2]=     34.8;
//  g.isotopic_abundance[2][3]=    9.033;
//  g.isotopic_abundance[2][4]= 7.199e+04;

    g.total_cross_section  = 1;
    g.cross_section_option = 011;

    g.t_half_limit         = 1.0e4;

    g.primary_electrons    = 1;
    g.secondary_positrons  = 0;
    g.secondary_electrons  = 0;
    g.knock_on_electrons   = 0;
    g.secondary_antiprotons= 0;
    g.tertiary_antiprotons = 0;
    g.secondary_protons    = 0;

    g.gamma_rays           = 0;
    g.pi0_decay            = 1;
    g.IC_isotropic         = 1;
    g.IC_anisotropic       = 0;
    g.synchrotron          = 1;
    g.bremss               = 1;

    g.DM_positrons         = 0;
    g.DM_electrons         = 0;
    g.DM_antiprotons       = 0;
    g.DM_gammas            = 0;

    g.skymap_format        = 0;
    g.output_gcr_full      = 0;
    g.warm_start           = 0;

    g.verbose              = 0;
    g.test_suite           = 0;


    g.primary_DM_positron    =0;			//HM 2014/4/24
    g.primary_DM_electron    =0;
    g.primary_DM_antip       =0;
    strcpy ( g.DM_positron_filename,"positron_source_2D.dat");
    strcpy ( g.DM_electron_filename,"electron_source_2D.dat");
    strcpy ( g.DM_antip_filename,"antiproton_source_2D.dat");


}










