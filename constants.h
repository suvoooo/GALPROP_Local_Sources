
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * constants.h *                                 galprop package * 08/16/2001 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#ifndef constants_h
#define constants_h

// aws routines

const double Rsun = 8.5; // kpc - Galactocentric radius of the Sun IMOS20080114

const double year2sec=365.25*24.*3600.; 
const double kpc2cm  =3.08568e21;
const double c       =3.0e10;

const double m_electron = 0.511;      // MeV
const double m_proton   = 938. ;      // MeV 

const double h_planck=6.6260755e-27 ; // erg sec
const double eV_to_erg=1.60217733e-12;
const double erg_to_eV=1./eV_to_erg;

// imos routines

const double // PHYSICAL CONSTANTS:
   Pi  = 3.141592653589793,    // Pi number
   C   = 2.99792458e10,        // cm/s, =c speed of light
   H   = 6.626075540e-27,      // erg*s, =h Planck constant
   H2Pi= 6.582122020e-22,      // MeV*s, =h/(2Pi) Planck constant
   ALPHAf=1./137.035989561,    // =e^2/(hc/2Pi) fine-structure const.
   Rele= 2.8179409238e-13,     // cm, =e^2/mc^2 class. electron radius
   MEV2ERG = 1.6021773349e-6,  // MeV/erg, conversion const.
   H2PiC   = 0.19732705359e-10,// MeV*cm, =hc/(2Pi) conversion const.
             // IONIZATION POTENTIALS:
//L.Pages et al. 1972, Atomic Data 4, 1 IMOS20061221
   EiH  = 18.9e-6,             // MeV, eff.Hydrogen (in electron energy losses)
   EiHe = 42.0e-6,             // MeV, eff.Helium (in electron energy losses)
   EH = 19.e-6,                // MeV, H  eff. ioniz. potential (in nucl. loss)
   EHe= 44.e-6,                // MeV, He eff. ioniz. potential (in nucl. loss)
             // RADIATIVE LENGHTS:
   TH = 62.8,                  // g/cm^2, Hydrogen
   THe= 93.1,                  // g/cm^2, Helium
             // MASSES: Particle Data Group (1990, Phys.Lett. 239)
   amu  = 0.931494,            // GeV, mass 12C /12= atomic mass unit  // IMOS20010816
   Mele = 0.5109990615,        // MeV/c^2, electron rest mass
   MD   = 1.232,               // GeV, Delta(1232)-isobar mean mass 
   GD   = 0.0575,              // GeV, Delta-isobar wight (0.115/2.)
   Mp   = 0.938,               // GeV, Proton rest mass
   Mn   = 0.9396,              // GeV, Neutron rest mass
   Md   = 1.8756,              // GeV, Deutron rest mass
   Mpi0 = 0.135,               // GeV, Neutral pion rest mass
   Mpi1 = 0.1396,              // GeV, Charged pion rest mass
   Mmu  = 0.10566,             // GeV, Muon rest mass
   MK   = 0.49365,             // GeV, Kaon rest mass
             // REACTION THRESHOLDS: (from kinematics)
   Pth0 = 0.78,                // GeV/c, min momentum for pp->pi^o (+2p)
   Pth1 = 1.65,                //(1.22)GeV/c, -"- for pp->pi^-(+2p+pi^+)
   Pth2 = 0.80,                // GeV/c, -"- for pp->pi^+ (+p+n)
   Pth3 = 0.791,               // GeV/c, -"- for pp->pi^+ (+d)
   Pth4 = 3.302,               // GeV/c, -"- for pp->K^- (+2p+K^+)
   Pth5 = 1.8332,              // GeV/c, -"- for pp->K^+ (+p+n)
             // BRANCHING RATIOS:
   BR1  = 0.635,               // br. ratio for channel K -> mu +nu
   BR2  = 0.212,               // br. ratio for channel K -> pi^(+/-) +pi^0
             // Helium FORM FACTORS: [Gould 1969, Phys.Rev.185,72 (p.77)]
   DD[11] = {0., 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10.},
   PH1[11]= {134.60, 133.85, 133.11, 130.86, 127.17, 
             120.35, 104.60,  89.94,  74.19,  54.26,  40.94},
   PH2[11]= {131.40, 130.51, 130.33, 129.26, 126.76,
             120.80, 105.21,  89.46,  73.03,  51.84,  37.24};


#endif
