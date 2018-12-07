
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Configure.h *                                 galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#ifndef Configure_h
#define Configure_h

#include<string> //IMOS20020112

namespace gp {

  class Configure {
    
  public:
    
    std::string fGaldefDirectory;
    std::string fFITSDataDirectory;
    //std::string adjunct_directory;
    std::string fOutputDirectory;
    std::string fOutputPrefix;
    std::string fVersion;
    std::string fGlobalDataPath;

    //char *galdef_directory;
    //char *fits_directory;
    //char *adjunct_directory;
    
    //int directory_length;
    
    //interface function prototype
    int init(const std::string& galdefDirectory,
	     const std::string& fitsDataDirectory,
	     const std::string& outputDirectory,
	     const std::string& outputPrefix);

  };

}

#endif
