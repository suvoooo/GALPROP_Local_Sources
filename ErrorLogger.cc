#include <string>
#include <iostream>
#include <sstream>

#include "ErrorLogger.h"

using std::string;
using std::ostream;
using std::ostringstream;

using namespace utl;


/** 
    Write and format an error message. For verbose message logging,
    the file name and line number are added to the beginning of the
    file. The format is such that (X)Emacs can jump directly to the
    corresponding line in the source file.
 */
//#warning lukas, why is severity parameter not used?
void
ErrorLogger::Log(const SeverityLevel severity,
                 const std::string functionName,
                 const std::string fileName,
                 const int lineNumber,
                 const std::string message,
                 VerbosityLevel verbosity) const {

  if (fVerbosity == eSilent)
    return;

  if (severity < fSeverity)
     return;

  ostream* myOStream;

  if (fOStream)
    myOStream = fOStream;
  else
    myOStream = &std::cerr;

  if (verbosity == eDefaultVerbosity) 
    verbosity = fVerbosity;

  ostringstream messageHeader;
  if (verbosity == eVerbose)
    messageHeader << fileName << ':' << lineNumber << ": ";

  messageHeader << functionName << ": ";
  string messageHeaderSpace;
  messageHeaderSpace.resize(messageHeader.str().size(), ' ');

  *myOStream << messageHeader.str() << message << std::endl;
}


// Configure (x)emacs for this file ...
// Local Variables:
// mode: c++
// compile-command: "make -k -C .. ErrorLogger/ErrorLogger.o"
// End:
