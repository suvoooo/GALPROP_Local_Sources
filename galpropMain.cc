#include <cstring>

#include <ErrorLogger.h>

#include <galprop_classes.h>
#include <galprop_internal.h>

#include <config.h>

void ShowHelp();

int main(int argc, char** argv) {

  //INFO("Entry");

  // Parse the command line. This enables override of built-in paths for
  // the FITS data, galdef files, output directory, and output prefix for
  // files from the current run

  string galdefPath, fitsPath, outputPath, outputPrefix, runNumber;

  bool showHelp = false;

  for (int i = 1; i < argc && *argv[i] == '-'; i++) {

    for (char *p = &argv[i][1]; *p != '\0'; p++) {
      
      switch (*p) {
	
      case 'f':
      case 'F':

	fitsPath = argv[++i];
      
      break;

      case 'g':
      case 'G':

        galdefPath = argv[++i];

      break;

      case 'h':
      case 'H':

	ShowHelp();
      exit(0);
      break;

      case 'o':
      case 'O':

	outputPath = argv[++i];
      
	break;

      case 'p':
      case 'P':

	outputPrefix = argv[++i];
      
	break;

      case 'R':
      case 'r':

	runNumber = argv[++i];
      
	break;

      default:

	ostringstream buf;
	buf << "Unspecified option " << argv[i];
	INFO(buf.str());
	exit(1);
	break;

      }

    }

  }

  // If run number is not set we can't proceed forward so report error
  // and exit

  if (runNumber.empty()) {

    INFO("Run number not specified. Exiting ...");
    exit(1);

  }

  INFO("Entry");

  Galprop galprop;
  // Reuse the input parameter detection for a program calculating the average
  // X factors

  // Compare the last proportion of the string only
  size_t iendProgramName = strlen(argv[0]);
  if ( strncmp("gpavXCO", argv[0]+iendProgramName-7, 7) == 0 ) {

     galprop.AvXCO(galdefPath, fitsPath, outputPath, outputPrefix, runNumber);

  } else {

     galprop.Run(galdefPath, fitsPath, outputPath, outputPrefix, runNumber);

  }

  INFO("Exit");

  return 0;

}

void ShowHelp() {

  ostringstream hBuf;

  hBuf << "Usage: galprop -r <run> -g <dir> -f <dir> -o <dir> -p <prefix>" << endl << "Where: " << endl << "-r <run> (the run number of the galdef file -- required)" << endl << "-g <dir> (location of the galdef directory -- default ../GALDEF)" << endl << "-f <dir> (location of the fits directory -- default ../FITS)" << endl << "-o <dir> (output directory)" << endl << "-p <prefix> (optional prefix for the output files)" << endl;

  cout << hBuf.str();

}
