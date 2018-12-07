#ifndef SkymapFitsIO_h
#define SkymapFitsIO_h

#include "Skymap.h"
#include <string>

/** \brief Write skymap to fits file.
 *
 * \param skymap is the skymap to output
 * \param filename is the filename to write to.  This routine will overwrite
 * files without a warning.
 * \param sUnit the units of the spectra to put in the fits file
 * \param sType the type of the spectra to put in the fits file
 * \return integer value telling the cfitsio return status.  A value of 0 means
 * all went well
 *
 * The output format is slightly modified from the standard healpix output.  We
 * still use a table extension to represent the skymap, but each pixel has its
 * own row in the table and there is a single column containing a vector
 * element which has the spectra for each pixel.  The extension for the healpix
 * table is called SKYMAP.  In addition, there is a ENERGIES extension,
 * containing a table of energies (frequencies) which the spectra is evaluated
 * at.
 */
int SkymapToFits (const Skymap<double> & skymap, const std::string & filename, const std::string & sUnit, const std::string & sType);
/** \overload */
int SkymapToFits (const Skymap<long> & skymap, const std::string & filename, const std::string & sUnit, const std::string & sType);
/** \brief Read a skymap from a fits file.
 *
 * \param skymap will be replaced with the one from the fits file
 * \param filename is the name of the file to read
 */
int FitsToSkymap (Skymap<double> & skymap, const std::string & filename);
/** \overload */
int FitsToSkymap (Skymap<long> & skymap, const std::string & filename);

#endif
