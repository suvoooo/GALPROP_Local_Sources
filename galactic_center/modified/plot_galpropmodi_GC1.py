#!/usr/bin/python3

# This is version 1.03.

# New in version 1.01: the ability to multiply all spectra by a user-specified power of energy (frequency).
# New in version 1.02: this version features a new parser, "ratios", to plot isotopic ratios.
# New in version 1.03: a bug is fixed which did not allow to calculate the flux the sum of all isotopes with modulation

import astropy.io.fits as pyfits
# import pyfits
import os.path
import numpy
import math
galpath = "/media/abbl/d17db44d-1f0c-4659-891d-b7ff26b02041/suvo/galprop-54/tarfiles/galprop_check/GALDEF/compareDRAGON/boron/new/"


class galpropDATA:
    def __init__(self, fitsdir, galdefid, parser, Rsun=0.25): # change Rsun to extract flux at that point # this goes with galactic center.  
        # Storage for the cosmic ray density
        # We store it as a dictionary containing a dictionary, containing an array, where the first key is Z, the second key is A and then we have primaries, secondaries in an array

        self.Rsun = Rsun
        if ((parser == "spectra") or (parser == "abundances") or (parser == "ratios")):
            self.density = {}
            # Open the file
            hdu = pyfits.open(os.path.join(fitsdir, "nuclei_"+galdefid))
            data = hdu[0].data
            # Find out which indices to interpolate over for Rsun
            R = (numpy.arange(int(hdu[0].header["NAXIS1"]))) * \
                hdu[0].header["CDELT1"] + hdu[0].header["CRVAL1"]
            inds = []
            weights = []
            if (R[0] > self.Rsun):
                inds.append(0)
                weights.append(1)
            elif (R[-1] <= self.Rsun):
                inds.append(-1)
                weights.append(1)
            else:
                for i in range(len(R)-1):
                    if (R[i] <= self.Rsun and self.Rsun < R[i+1]):
                        inds.append(i)
                        inds.append(i+1)
                        weights.append((R[i+1]-self.Rsun)/(R[i+1]-R[i]))
                        weights.append((self.Rsun-R[i])/(R[i+1]-R[i]))
                        break
            # Calculate the energy for the spectral points
            self.energy = 10**(float(hdu[0].header["CRVAL3"]) + numpy.arange(
                int(hdu[0].header["NAXIS3"]))*float(hdu[0].header["CDELT3"]))
            print(self.energy)
            # Parse the header, looking for Nuclei definitions
            Nnuclei = hdu[0].header["NAXIS4"]
            for i in range(1, Nnuclei+1):
                id = "%03d" % i
                Z = int(hdu[0].header["NUCZ"+id])
                A = int(hdu[0].header["NUCA"+id])
                K = int(hdu[0].header["NUCK"+id])
                # Add the data to the density dictionary
                if Z not in self.density:
                    self.density[Z] = {}
                if A not in self.density[Z]:
                    self.density[Z][A] = {}
                if K not in self.density[Z][A]:
                    self.density[Z][A][K] = []
                d = ((data[i-1, :, 0, inds].swapaxes(0, 1))
                     * numpy.array(weights)).sum(axis=1)
                self.density[Z][A][K].append(d/self.energy**2)

        if (parser == "synchrotron"):
            # Search for synchrotron spectrum, start with standard output, then mapcubes, and finally healpix
            # Also search for Q and U
            self.synchrotron = {}
            if (os.path.exists(os.path.join(fitsdir, "synchrotron_skymap_"+galdefid))):
                sf = pyfits.open(os.path.join(
                    fitsdir, "synchrotron_skymap_"+galdefid))
                self.synchrotron["total"] = sf[0].data.sum(axis=3).sum(
                    axis=2)[0]/sf[0].data.shape[2]/sf[0].data.shape[3]
                self.synchrotron["nu"] = 10**(float(sf[0].header["CRVAL3"]) + numpy.arange(
                    int(sf[0].header["NAXIS3"]))*float(sf[0].header["CDELT3"]))
                if (os.path.exists(os.path.join(fitsdir, "synchrotron_Q_skymap_"+galdefid))):
                    sf = pyfits.open(os.path.join(
                        fitsdir, "synchrotron_Q_skymap_"+galdefid))
                    self.synchrotron["Q"] = sf[0].data.sum(axis=3).sum(
                        axis=2)[0]/sf[0].data.shape[2]/sf[0].data.shape[3]
                    sf = pyfits.open(os.path.join(
                        fitsdir, "synchrotron_U_skymap_"+galdefid))
                    self.synchrotron["U"] = sf[0].data.sum(axis=3).sum(
                        axis=2)[0]/sf[0].data.shape[2]/sf[0].data.shape[3]
            elif (os.path.exists(os.path.join(fitsdir, "synchrotron_mapcube_"+galdefid))):
                sf = pyfits.open(os.path.join(
                    fitsdir, "synchrotron_mapcube_"+galdefid))
                self.synchrotron["total"] = sf[0].data.sum(axis=2).sum(
                    axis=1)/sf[0].data.shape[1]/sf[0].data.shape[2]
                self.synchrotron["nu"] = sf[1].data.field(0)
                if (os.path.exists(os.path.join(fitsdir, "synchrotron_Q_mapcube_"+galdefid))):
                    sf = pyfits.open(os.path.join(
                        fitsdir, "synchrotron_Q_mapcube_"+galdefid))
                    self.synchrotron["Q"] = sf[0].data.sum(axis=2).sum(
                        axis=1)/sf[0].data.shape[1]/sf[0].data.shape[2]
                    sf = pyfits.open(os.path.join(
                        fitsdir, "synchrotron_U_mapcube_"+galdefid))
                    self.synchrotron["U"] = sf[0].data.sum(axis=2).sum(
                        axis=1)/sf[0].data.shape[1]/sf[0].data.shape[2]
            elif (os.path.exists(os.path.join(fitsdir, "synchrotron_healpix_"+galdefid))):
                sf = pyfits.open(os.path.join(
                    fitsdir, "synchrotron_healpix_"+galdefid))
                self.synchrotron["total"] = sf[1].data.field(
                    0).sum(axis=0)/float(sf[1].header["NAXIS2"])
                self.synchrotron["nu"] = sf[2].data.field(0)
                if (os.path.exists(os.path.join(fitsdir, "synchrotron_Q_healpix_"+galdefid))):
                    sf = pyfits.open(os.path.join(
                        fitsdir, "synchrotron_Q_healpix_"+galdefid))
                    self.synchrotron["Q"] = sf[1].data.field(
                        0).sum(axis=0)/float(sf[1].header["NAXIS2"])
                    sf = pyfits.open(os.path.join(
                        fitsdir, "synchrotron_U_healpix_"+galdefid))
                    self.synchrotron["U"] = sf[1].data.field(
                        0).sum(axis=0)/float(sf[1].header["NAXIS2"])

        if (parser == "gamma"):
            # Now search for gamma_rays, again starting with standard output, then mapcubes, and finally healpix
            self.gamma_rays = {}
            if (os.path.exists(os.path.join(fitsdir, "ics_isotropic_skymap_"+galdefid))
                or os.path.exists(os.path.join(fitsdir, "bremss_skymap_"+galdefid))
                    or os.path.exists(os.path.join(fitsdir, "pion_decay_skymap_"+galdefid))):
                if (os.path.exists(os.path.join(fitsdir, "ics_isotropic_skymap_"+galdefid))):
                    sf = pyfits.open(os.path.join(
                        fitsdir, "ics_isotropic_skymap_"+galdefid))
                    self.gamma_rays["energy"] = 10**(float(sf[0].header["CRVAL3"]) + numpy.arange(
                        int(sf[0].header["NAXIS3"]))*float(sf[0].header["CDELT3"]))
                    b = float(sf[0].header["CRVAL2"]) + numpy.arange(
                        int(sf[0].header["NAXIS2"]))*float(sf[0].header["CDELT2"])
                    self.gamma_rays["IC"] = (sf[0].data.sum(axis=3)*numpy.sin(numpy.pi/2.-numpy.radians(b))).sum(axis=2)[
                        0]/sf[0].data.shape[3]/numpy.sin(numpy.pi/2.-numpy.radians(b)).sum()/self.gamma_rays["energy"]**2
                if (os.path.exists(os.path.join(fitsdir, "bremss_skymap_"+galdefid))):
                    sf = pyfits.open(os.path.join(
                        fitsdir, "bremss_skymap_"+galdefid))
                    self.gamma_rays["energy"] = 10**(float(sf[0].header["CRVAL3"]) + numpy.arange(
                        int(sf[0].header["NAXIS3"]))*float(sf[0].header["CDELT3"]))
                    b = float(sf[0].header["CRVAL2"]) + numpy.arange(
                        int(sf[0].header["NAXIS2"]))*float(sf[0].header["CDELT2"])
                    self.gamma_rays["bremss"] = (sf[0].data.sum(axis=3)*numpy.sin(numpy.pi/2.-numpy.radians(b))).sum(
                        axis=2)[0]/sf[0].data.shape[3]/numpy.sin(numpy.pi/2.-numpy.radians(b)).sum()/self.gamma_rays["energy"]**2
                if (os.path.exists(os.path.join(fitsdir, "pion_decay_skymap_"+galdefid))):
                    sf = pyfits.open(os.path.join(
                        fitsdir, "pion_decay_skymap_"+galdefid))
                    self.gamma_rays["energy"] = 10**(float(sf[0].header["CRVAL3"]) + numpy.arange(
                        int(sf[0].header["NAXIS3"]))*float(sf[0].header["CDELT3"]))
                    b = float(sf[0].header["CRVAL2"]) + numpy.arange(
                        int(sf[0].header["NAXIS2"]))*float(sf[0].header["CDELT2"])
                    self.gamma_rays["pion_decay"] = (sf[0].data.sum(axis=3)*numpy.sin(numpy.pi/2.-numpy.radians(b))).sum(
                        axis=2)[0]/sf[0].data.shape[3]/numpy.sin(numpy.pi/2.-numpy.radians(b)).sum()/self.gamma_rays["energy"]**2
            elif (os.path.exists(os.path.join(fitsdir, "ics_isotropic_mapcube_"+galdefid))
                  or os.path.exists(os.path.join(fitsdir, "bremss_mapcube_"+galdefid))
                  or os.path.exists(os.path.join(fitsdir, "pion_decay_mapcube_"+galdefid))):
                if (os.path.exists(os.path.join(fitsdir, "ics_isotropic_mapcube_"+galdefid))):
                    sf = pyfits.open(os.path.join(
                        fitsdir, "ics_isotropic_mapcube_"+galdefid))
                    b = float(sf[0].header["CRVAL2"]) + (numpy.arange(int(sf[0].header["NAXIS2"]))-int(
                        sf[0].header["CRPIX2"]))*float(sf[0].header["CDELT2"])
                    self.gamma_rays["IC"] = (sf[0].data.sum(axis=2)*numpy.sin(numpy.pi/2.-numpy.radians(
                        b))).sum(axis=1)/sf[0].data.shape[2]/numpy.sin(numpy.pi/2.-numpy.radians(b)).sum()
                    self.gamma_rays["energy"] = sf[1].data.field(0)
                if (os.path.exists(os.path.join(fitsdir, "bremss_mapcube_"+galdefid))):
                    sf = pyfits.open(os.path.join(
                        fitsdir, "bremss_mapcube_"+galdefid))
                    b = float(sf[0].header["CRVAL2"]) + (numpy.arange(int(sf[0].header["NAXIS2"]))-int(
                        sf[0].header["CRPIX2"]))*float(sf[0].header["CDELT2"])
                    self.gamma_rays["bremss"] = (sf[0].data.sum(axis=2)*numpy.sin(numpy.pi/2.-numpy.radians(
                        b))).sum(axis=1)/sf[0].data.shape[2]/numpy.sin(numpy.pi/2.-numpy.radians(b)).sum()
                    self.gamma_rays["energy"] = sf[1].data.field(0)
                if (os.path.exists(os.path.join(fitsdir, "pion_decay_mapcube_"+galdefid))):
                    sf = pyfits.open(os.path.join(
                        fitsdir, "pion_decay_mapcube_"+galdefid))
                    b = float(sf[0].header["CRVAL2"]) + (numpy.arange(int(sf[0].header["NAXIS2"]))-int(
                        sf[0].header["CRPIX2"]))*float(sf[0].header["CDELT2"])
                    self.gamma_rays["pion_decay"] = (sf[0].data.sum(axis=2)*numpy.sin(numpy.pi/2.-numpy.radians(
                        b))).sum(axis=1)/sf[0].data.shape[2]/numpy.sin(numpy.pi/2.-numpy.radians(b)).sum()
                    self.gamma_rays["energy"] = sf[1].data.field(0)
            if (os.path.exists(os.path.join(fitsdir, "ics_isotropic_healpix_"+galdefid))
                or os.path.exists(os.path.join(fitsdir, "bremss_healpix_"+galdefid))
                    or os.path.exists(os.path.join(fitsdir, "pion_decay_healpix_"+galdefid))):
                if (os.path.exists(os.path.join(fitsdir, "ics_isotropic_healpix_"+galdefid))):
                    sf = pyfits.open(os.path.join(
                        fitsdir, "ics_isotropic_healpix_"+galdefid))
                    self.gamma_rays["IC"] = sf[1].data.field(
                        0).sum(axis=0)/float(sf[1].header["NAXIS2"])
                    self.gamma_rays["energy"] = sf[2].data.field(0)
                if (os.path.exists(os.path.join(fitsdir, "bremss_healpix_"+galdefid))):
                    sf = pyfits.open(os.path.join(
                        fitsdir, "bremss_healpix_"+galdefid))
                    self.gamma_rays["bremss"] = sf[1].data.field(
                        0).sum(axis=0)/float(sf[1].header["NAXIS2"])
                    self.gamma_rays["energy"] = sf[2].data.field(0)
                if (os.path.exists(os.path.join(fitsdir, "pi0_decay_healpix_"+galdefid))):
                    sf = pyfits.open(os.path.join(
                        fitsdir, "pi0_decay_healpix_"+galdefid))
                    self.gamma_rays["pion_decay"] = sf[1].data.field(
                        0).sum(axis=0)/float(sf[1].header["NAXIS2"])
                    self.gamma_rays["energy"] = sf[2].data.field(0)

            if ("energy" in self.gamma_rays):
                self.gamma_rays["total"] = numpy.zeros(
                    self.gamma_rays["energy"].shape)
                if ("IC" in self.gamma_rays):
                    self.gamma_rays["total"] += self.gamma_rays["IC"]
                if ("bremss" in self.gamma_rays):
                    self.gamma_rays["total"] += self.gamma_rays["bremss"]
                if ("pion_decay" in self.gamma_rays):
                    self.gamma_rays["total"] += self.gamma_rays["pion_decay"]


# Return the spectra for species, modulated if needed
# If joined is true (the default), then we return all secondaries, primaries summed together, otherwise we have an array of spectra
# Phi is the modulation potential, defaulting to 0
# out_energy is the output energy at which to evaluate the spectra.  Defaults to self.energy if empty.
# if A = -1, we sum up all isotopes
# if K = -1, we sum up all K-electron captures
# returns a tuple with (energies, spectra, type), spectra is a 2D numpy array
# with first index separating components and second index showing the spectra.
# type is a 2D array containing a tuple (Z,A,K) for all components (first
# dimension the same size as spectra).  A=-1, K=-1 indicate sum over all
# isotopes,K-capture electrons.  First component is always the sum.
# if joined is true, there is only one component
# if there are many components with same type, the first is secondary, the
# second primary.
# Unfortunately there is no way to tell which is which from the nuclei file.


    def CRspectra(self, Z, A, K=-1, joined=True, phi=0, out_energy=[]):
        pmass = 939.
        emass = 0.5109990615

        # Assign the energy array
        if (out_energy == []):
            out_energy = self.energy
        else:
            out_energy = numpy.array(out_energy)

        # Return -1 if Z does not exist
        if (Z not in self.density):
            return (out_energy, -numpy.ones((1, len(out_energy))), [(Z, -1, -1)])

        if (joined):
            spectra = numpy.zeros((1, self.energy.shape[0]), dtype='float')
            type = [(Z, A, K)]
        else:
            spectra = []
            spectra.append(numpy.zeros(self.energy.shape, dtype='float'))
            type = []
            type.append((Z, A, K))

        # These are the same as spectra and type, except they are only used internally
        # and not returned. They are needed to do the modulation properly when A=-1 is set
        spectra_brdn = []
        spectra_brdn.append(numpy.zeros(self.energy.shape, dtype='float'))
        type_brdn = []
        type_brdn.append((Z, A, K))

        if (A == -1):
            tmpA = self.density[Z]
        else:
            if (A not in self.density[Z]):
                return (out_energy, -numpy.ones((1, len(out_energy))), [(Z, A, -1)])
            tmpA = {A: self.density[Z][A]}

        for aA, vA in tmpA.items():
            if (K == -1):
                tmpK = vA
            else:
                if (K not in vA):
                    return (out_energy, -numpy.ones((1, len(out_energy))), [(Z, A, K)])
                tmpK = {K: vA[K]}
            for kK, vK in tmpK.items():
                for sp in vK:
                    spectra[0] += sp
                    if (not joined):
                        spectra.append(sp)
                        type.append((Z, aA, kK))
                    spectra_brdn.append(sp)
                    type_brdn.append((Z, aA, kK))

        # Do interpolation in energy for all isotopes and Ks:
        out_brdn = numpy.zeros(
            (len(spectra_brdn), len(out_energy)), dtype='float')
        for k in range(1, len(spectra_brdn)):
            aA = type_brdn[k][1]
            if (aA > 0):
                energy = out_energy + math.fabs(Z)*phi/float(aA)
            elif (aA == 0):
                energy = out_energy + math.fabs(Z)*phi
            else:
                return (out_energy, -numpy.ones((1, len(out_energy))), [(Z, A, K)])

            for (i, en) in enumerate(energy):
                if (en >= self.energy[0] and en <= self.energy[-1]):
                    j = 0
                    while (en > self.energy[j]):
                        j += 1
                    if (j == 0):
                        j = 1

                    if (spectra_brdn[k][j-1] > 0 and spectra_brdn[k][j] > 0):
                        sl = math.log(
                            spectra_brdn[k][j-1]/spectra_brdn[k][j])/math.log(self.energy[j-1]/self.energy[j])
                        out_brdn[k, i] = spectra_brdn[k][j-1] * \
                            (en/self.energy[j-1])**sl
                    elif (spectra_brdn[k][j-1] > 0):
                        out_brdn[k, i] = (
                            self.energy[j]-en)*spectra_brdn[k][j-1]/(self.energy[j]-self.energy[j-1])
                    elif (spectra_brdn[k][j] > 0):
                        out_brdn[k, i] = (
                            en - self.energy[j-1])*spectra_brdn[k][j]/(self.energy[j]-self.energy[j-1])

                    if (aA > 0):
                        m = pmass
                    else:
                        m = emass

                    out_brdn[k, i] *= out_energy[i] * \
                        (out_energy[i]+2.*m)/(en*(en+2.*m))
                    out_brdn[0, i] += out_brdn[k, i]

                else:
                    out_brdn[k, i] = -1
                    out_brdn[0, i] = -1

        # Return only the requested values
        out = numpy.zeros((len(spectra), len(out_energy)), dtype='float')
        for k in range(len(spectra)):
            for i in range(len(energy)):
                out[k, i] = out_brdn[k, i]

        return (out_energy, out, type)

    def CRIsotopes(self, energy, phi=0):
        # Returns the isotope fraction for all Z,A pairs.  Sum up all K electrons and secondaries
        # protons should always be in the table
        # Format of return array is (Z,A,fraction)

        en = self.CRspectra(1, 1, joined=True, out_energy=[
                            energy], phi=phi)[0][0]
        pr = self.CRspectra(1, 1, joined=True, out_energy=[
                            energy], phi=phi)[1][0, 0]

        out = []
        for Z, vZ in self.density.items():
            for A, vA in vZ.items():
                fr = self.CRspectra(Z, A, joined=True, out_energy=[
                                    energy], phi=phi)[1][0, 0]
                fr /= pr
                out.append((Z, A, fr))
        return out


def require(condition):
    if (condition != True):
        print("Plot Galprop version 1.03 (October 18, 2010)")
        print("\nThis is the program used by GALPROP WebRun function Quick Plots at http://galprop.stanford.edu/webrun/")
        print("It reads the FITS files produced by GALPROP and outputs the read information as data tables.")
        print("You are free to use this program as a reference, or modify it to serve your needs.")
        print("Make sure you acknowledge GALPROP by following instructions at http://galprop.stanford.edu/code.php?option=terms")
        print("if you use the GALPROP code, GALROP WebRun service, or this program in your publications.\n")
        print("Usage: ", argv[0],
              " fitsdir galdefid parser [arg1 [arg2 [ ... ]]]")
        print("   fitsdir  : directory where the data files are located")
        print("   galdefid : the common suffix of all the fits files (including .gz if archived)")
        print("   parser   : one of the following: spectra, abundances, ratios, gamma, synchrotron")
        print("   arg1, arg2, etc: depends on the parser\n")
        print("For parser 'spectra',    arg1 is the modulation potential in MV, arg2 is Z, arg3 is A and arg4 is alpha")
        print("For parser 'abundances', arg1 is the modulation potential in MV, arg2 is the kinetic energy per nucleon in MeV")
        print("For parser 'ratios', arg1 is the modulation potential, arg2 and arg3 are Z1 and A1 (numerator), arg4 and arg5 are Z2 and A2 (denominator)\n")
        print("For parsers 'gamma' and 'synchrotron', arg1 is alpha\n")
        print("Example: ", argv[0],
              " ./my_results 54_00010001.gz spectra 250 2 4 2.7")
        print(
            "         will read ./my_results/nuclei_54_00010001.gz and output the spectrum")
        print("         of Helium-4 (Z=2, A=4) at the location of the Sun, assuming Solar modulation")
        print("         with potential V=250 MV and multiplying the spectrum by E^2.7.")
        exit()


if (__name__ == "__main__"):
    from sys import argv

    require(len(argv) > 3)

    fitsdir = argv[1]
    galdefid = argv[2]
    parser = argv[3]

    gdsp = galpropDATA(fitsdir, galdefid, parser)

    if (parser == "spectra"):
        require(len(argv) == 7+1)

        phi = float(argv[4])
        Z = int(argv[5])
        A = int(argv[6])
        alpha = float(argv[7])

        (energy, spectra, type) = gdsp.CRspectra(Z, A, joined=True, phi=phi)
        spectra *= energy**alpha
        print("# Energy spectrum of a cosmic ray isotope calculated in GALPROP run # %s" % galdefid)
        print("# for Z=%3d, A=%3d, V=%10.3e MV, R=%6.2f kpc. " %
              (Z, A, phi, gdsp.Rsun))
        print("# The first column is E, the kinetic energy per nucleon after modulation, in MeV,")
        print("# and the second column contains E^%4.2f*dF/dE, where dF/dE is the spectral density of the flux of the isotope." % (alpha))
        print("# The units of dF/dE are cm^-2*s^-1*sr^-1*MeV^-1, so E^%4.2f*dF/dE is in MeV^%4.2f*cm^-2*s^-1*sr^-1*Mev^-1." % (alpha, alpha))
        print("# %10s %12s" % ("Energy", "Flux"))
        with open(galpath+galdefid+".d", "w") as o:
            for k in range(len(energy)):
                str = "%12.4e %12.4e " % (energy[k], spectra[0][k])
                for j in range(1, spectra.shape[0]):
                    str += "%12.4 " % spectra[j][k]
                print(str)
                o.write(str + "\n")
    elif (parser == "ratios"):
        require(len(argv) == 8+1)

        phi = float(argv[4])
        Z1 = int(argv[5])
        A1str = argv[6]
        Z2 = int(argv[7])
        A2str = argv[8]

        A1list = A1str.split(",")
        A2list = A2str.split(",")

        n_found_1 = 0
        n_found_2 = 0

        for A1 in A1list:
            (energy, spectra_breakdown, type) = gdsp.CRspectra(
                Z1, int(A1), joined=True, phi=phi)
            if (n_found_1 == 0):
                # This checks for non-existent isotopes
                if (spectra_breakdown[0][0] > 0):
                    spectra1 = spectra_breakdown
                    n_found_1 += 1
            else:
                if (spectra_breakdown[0][0] > 0):
                    for k in range(len(energy)):
                        for j in range(spectra1.shape[0]):
                            spectra1[j][k] += spectra_breakdown[j][k]
                    n_found_1 += 1

        for A2 in A2list:
            (energy, spectra_breakdown, type) = gdsp.CRspectra(
                Z2, int(A2), joined=True, phi=phi)
            if (n_found_2 == 0):
                if (spectra_breakdown[0][0] > 0):
                    spectra2 = spectra_breakdown
                    n_found_2 += 1
            else:
                if (spectra_breakdown[0][0] > 0):
                    for k in range(len(energy)):
                        for j in range(spectra2.shape[0]):
                            spectra2[j][k] += spectra_breakdown[j][k]
                    n_found_2 += 1

        if (n_found_1*n_found_2 == 0):
            exit()

        print("# Flux ratios cosmic ray isotopes calculated in GALPROP run # %s" % galdefid)
        print("# for Z1=%3d, A1=%s, Z2=%3d, A2=%s, V=%10.3e MV, R=%6.2f kpc. " % (
            Z1, A1str, Z2, A2str, phi, gdsp.Rsun))
        print("# The first column is E, the kinetic energy per nucleon after modulation, in MeV,")
        print("# and the second column contains the dimensionless ratio of (Z1; A1) to (Z2; A2).")
        print("# %10s %12s" % ("Energy", "Ratio"))
        for k in range(len(energy)):
            if ((spectra1[0][k] > 0) and (spectra2[0][k] > 0)):
                str = "%12.4e %12.4e " % (
                    energy[k], spectra1[0][k]/spectra2[0][k])
                for j in range(1, spectra1.shape[0]):
                    str += "%12.4 " % (spectra1[j][k]/spectra2[j][k])
                print(str)
    elif (parser == "abundances"):
        require(len(argv) == 5+1)

        phi = float(argv[4])
        E = float(argv[5])
        isotopes = gdsp.CRIsotopes(E, phi)
        pr = gdsp.CRspectra(1, 1, joined=True, out_energy=[
                            E], phi=phi)[1][0, 0]
        print("# Number fluxes of isotopes calculated in GALPROP run # %s" % galdefid)
        print("# for E=%10.3e MeV/nucleon, V=%10.3e MV, R=%6.2f kpc" %
              (E, phi, gdsp.Rsun))
        print("# All fluxes are normalized to the proton flux equal to F=%12.4e cm^-2 s^-1 sr^-1 MeV^-1" % pr)
        print("# %4s %6s %12s" % ("Z", "A", "Flux"))
        for iso in isotopes:
            print("%6d %6d %12.4e" % iso)

    elif (parser == "gamma"):
        require(len(argv) == 4+1)

        alpha = float(argv[4])

        if "energy" in gdsp.gamma_rays:
            cmp = gdsp.gamma_rays.keys()
            cmpkeys = list(cmp)
            print(cmp)
            cmpkeys.remove("energy")
            cmpkeys.remove("total")
            print(
                "# Gamma-ray emission spectra calculated in GALPROP run # %s" % galdefid)
            print(
                "# The spectra represent the averages over the region of the sky specified in the GALDEF file,")
            print("# or over the whole sky, if the HEALPIX format was used for skymaps.")
            print(
                "# The first column if the photon energy E in MeV, and columns 2, 3, 4 and 5 represent")
            print(
                "# the total, pion decay, inverse Compton and Bremsstrahlung spectra, respectively.")
            print("# These spectra are multiplied by E^%4.2f and have the units of MeV^%4.2f*cm^-2*s^-1*sr^-1*MeV^-1." % (alpha, alpha))
            str = "#Energy [MeV]  total flux*E^%4.2f [MeV^%4.2f cm^-2 s^-1 sr^-1 MeV^-1]" % (
                alpha, alpha)
            for c in cmp:
                str += " "+c
            print(str)
            print("# %10s %12s %12s %12s %12s" %
                  ("E", "Total", "Pion Decay", "Inv. Compt.", "Bremsstr."))
            for j in range(len(gdsp.gamma_rays["energy"])):
                str = "%12.4e %12.4e " % (gdsp.gamma_rays["energy"][j], (
                    gdsp.gamma_rays["total"][j]*(gdsp.gamma_rays["energy"][j])**alpha))
                str += "%12.4e" % (gdsp.gamma_rays["pion_decay"]
                                   [j]*(gdsp.gamma_rays["energy"][j])**alpha)
                # str += "%12.4e" % (gdsp.gamma_rays["IC"][j]
                                #    * (gdsp.gamma_rays["energy"][j])**alpha)
                # str += "%12.4e" % (gdsp.gamma_rays["bremss"]
                                #    [j]*(gdsp.gamma_rays["energy"][j])**alpha)
                print(str)

    elif (parser == "synchrotron"):
        require(len(argv) == 4+1)

        alpha = float(argv[4])

        if "total" in gdsp.synchrotron:
            print(
                "# Synchrotron radiation spectrum calculated in GALPROP run # %s" % galdefid)
            print(
                "# The spectrum represents the average over the region of the sky specified in the GALDEF file,")
            print("# or over the whole sky, if the HEALPIX format was used for skymaps.")
            print(
                "# The first column in the wave frequency f in Hz, and the second column is the")
            print("# product of f^%4.2f and flux density F in Hz^%4.2f*erg*cm^-2*s^-1*sr^-1*Hz^-1." % (alpha, alpha))
            str = "f^%4.2f*Flux" % (alpha)
            print("# %10s %12s" % ("Frequency", str))
            for j in range(len(gdsp.synchrotron["total"])):
                print("%12.4e %12.4e" % (
                    gdsp.synchrotron["nu"][j], (gdsp.synchrotron["total"][j])*(gdsp.synchrotron["nu"][j])**alpha))
