//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Configure.cc *                                galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#include <Configure.h>

#include <ErrorLogger.h>

#include <config.h>

#include <iostream>
#include <string>

using namespace std;
using namespace gp;

int Configure::init(const string& galdefDirectory,
		    const string& fitsDataDirectory,
		    const string& outputDirectory,
		    const string& outputPrefix) {

  INFO("Entry");

#ifdef VERSION
  const string fullVersion = VERSION;
#else
  const string fullVersion = "54.0.0";
#endif

  ostringstream bufV;
  bufV << "This is Galprop version " << fullVersion;
  INFO(bufV.str());

  //  pre-r151
  //  char version[]="53"; //AWS 20070820
  //  for changes from r151, new version defined at r158
  
  //char version[]="54"; //AWS 20080618

  const string majorVersion = fullVersion.substr(0, fullVersion.find_first_of("."));

  const string tmpVer = fullVersion.substr(fullVersion.find_first_of(".")+1, fullVersion.size()-1);

  const string minorVersion = tmpVer.substr(0, tmpVer.find_first_of("."));

  const string revision = tmpVer.substr(tmpVer.find_first_of(".")+1, tmpVer.size()-1);

  //cout << majorVersion << " " << minorVersion << " " << revision << endl;

  // Check for overrides. Otherwise grab whatever comes from the 
  // autoconfiguration, and if that fails use the stupid defaults we've
  // always had

  string galdefDir = galdefDirectory, fitsDir = fitsDataDirectory, outputDir = outputDirectory;

  if (galdefDir.empty()) {

#ifdef HAVE_GALDEF
    galdefDir = GALDEF_PATH;
#else
    galdefDir = "../GALDEF";
#endif

  } 

  if (fitsDir.empty()) {

#ifdef HAVE_FITSDATA
    fitsDir = FITSDATA_PATH;
#else
    fitsDir = "../FITS";
#endif

  }

  if (outputDir.empty()) {

#ifdef HAVE_FITSDATA
    outputDir = FITSDATA_PATH;
#else
    outputDir = "../FITS"; // *sigh*
#endif

  }

  //directory_length=100;
  //galdef_directory=new char[directory_length];
  //fits_directory=new char[directory_length];
  //adjunct_directory=new char[directory_length];
  //strcpy(galdef_directory, "../GALDEF/");
  //strcpy(  fits_directory, "../FITS/"  );
  //strcpy(adjunct_directory,"../adjunct/"  );

  fGaldefDirectory = galdefDir + "/";
  fFITSDataDirectory = fitsDir + "/";
  fOutputDirectory = outputDir + "/";
  fOutputPrefix = outputPrefix;
  fVersion = majorVersion;

  fGlobalDataPath = DATA_PATH;

  //cout << fFITSDataDirectory << " " << fOutputDirectory << " " << fOutputPrefix << " " << fGaldefDirectory << " " << fGlobalDataPath << endl;

  //galdef_directory = galdefDirectory;//"../GALDEF/";
  //fits_directory = fitsDataDirectory;//"../FITS/";
  //adjunct_directory = "../adjunct/";
  
  //cout<<"Configure: galdef_directory:  "<< galdef_directory<<endl;
  //cout<<"Configure:   fits_directory:  "<<   fits_directory<<endl;
  //cout<<"Configure:adjunct_directory:  "<<adjunct_directory<<endl;
  //cout<<"<<<<Configure"<<endl;

  INFO("Exit");

  return 0;

}
