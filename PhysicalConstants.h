#ifndef _utl_PhysicalConstants_h_
#define _utl_PhysicalConstants_h_

#include <Units.h>

#include <cmath>

namespace utl {

  const double kPi = M_PI;
  const double kExp = exp(1.0); 
  const double kSqrtTwo = sqrt(2.0);
  const double kLogTwo = log(2.0);
  const double kSqrtThree = sqrt(3.0);
  const double kOneOnThree = 1.0/3.0;
  const double kFourPiOnThree = 4.0*kPi/3.0;
  const double kPiOnTwo = kPi/2.0;
  const double kTwoPi = 2.0*kPi;
  const double kFourPi = 2.0*kTwoPi;
  const double kOneOnPi = 1.0/kPi;
  const double kOneOnTwoPi = 1.0/kTwoPi;
  const double kOneOnFourPi = 1.0/kFourPi;
  const double kConvertDegreesToRadians = kPi/180.0;

  // All taken from PDG data tables (2002)

  // Physical constants 

  const double kSpeedOfLight_SI = 299792458.0;
  const double kSpeedOfLight = kSpeedOfLight_SI*m/s;
  const double kPlanck_SI = 6.62606876e-34;
  const double kPlanckReduced_SI = kPlanck_SI*kOneOnTwoPi;
  const double kPlanck = kPlanck_SI*joule*s;
  const double kPlanckReduced = kPlanckReduced_SI*joule*s;
  const double kMuZero_SI = 4.0*kPi*1.0e-7;
  const double kMuZero = kMuZero_SI*newton/(ampere*ampere); 
  const double kBoltzmann_SI = 1.3806503e-23;
  const double kBoltzmann = kBoltzmann_SI*joule/kelvin;

  const double kPlanckTimesSpeedOfLight_SI = kPlanck_SI*kSpeedOfLight_SI;
  const double kPlanckTimesSpeedOfLightSquared_SI = 
    kPlanckTimesSpeedOfLight_SI*kPlanckTimesSpeedOfLight_SI;
  
  // Particle and other masses

  const double kMassConversion_SI = e_SI/(kSpeedOfLight_SI*kSpeedOfLight_SI);

  const double kHydrogenMass_SI = 1.6735e-27;
  const double kHydrogenMass = kHydrogenMass_SI*kg;
  const double kSolarMass_SI = 1.98892e30;
  const double kSolarMass = kSolarMass_SI*kg;

  const double kElectronMass = 0.510998902*MeV; 
  const double kElectronMass_SI = kElectronMass*kMassConversion_SI;
  const double kMuonMass = 105.658357*MeV; 
  const double kMuonMass_SI = kMuonMass*kMassConversion_SI;
  const double kTauMass = 1776.99*MeV;
  const double kTauMass_SI = kTauMass*kMassConversion_SI;

  const double kProtonMass = 938.271998*MeV; 
  const double kProtonMass_SI = kProtonMass*kMassConversion_SI;
  const double kNeutronMass = 939.56533*MeV; 
  const double kNeutronMass_SI = kNeutronMass*kMassConversion_SI;
  const double kDeuteronMass = 1875.612762*MeV; 
  const double kDeuteronMass_SI = kDeuteronMass*kMassConversion_SI;

  const double kLambdaMass = 1115.683*MeV;
  const double kLambdaMass_SI = kLambdaMass*kMassConversion_SI;
  const double kSigmaZeroMass = 1192.642*MeV;
  const double kSigmaZeroMass_SI = kSigmaZeroMass*kMassConversion_SI;
  const double kSigmaPlusMass = 1189.37*MeV;
  const double kSigmaPlusMass_SI = kSigmaPlusMass*kMassConversion_SI;
  const double kSigmaMinusMass = 1197.449*MeV;
  const double kSigmaMinusMass_SI = kSigmaMinusMass*kMassConversion_SI;
  const double kXiZeroMass = 1314.83*MeV;
  const double kXiZeroMass_SI = kXiZeroMass*kMassConversion_SI;
  const double kXiMinusMass = 1321.31*MeV;
  const double kXiMinusMass_SI = kXiMinusMass*kMassConversion_SI;
  const double kOmegaMinusMass = 1672.45*MeV;
  const double kOmegaMinusMass_SI = kOmegaMinusMass*kMassConversion_SI;

  const double kPiZeroMass = 134.9766*MeV; 
  const double kPiZeroMass_SI = kPiZeroMass*kMassConversion_SI;
  const double kPiChargedMass = 139.57018*MeV; 
  const double kPiChargedMass_SI = kPiChargedMass*kMassConversion_SI;
  const double kKaonChargedMass = 493.677*MeV; 
  const double kKaonChargedMass_SI = kKaonChargedMass*kMassConversion_SI;

  const double kAtomicMassUnit_SI = 1.660538e-27;

  const double kCarbonMass_SI = 12.0107*kAtomicMassUnit_SI;
  const double kOxygenMass_SI = 15.9994*kAtomicMassUnit_SI;
  const double kMagnesiumMass_SI = 24.3050*kAtomicMassUnit_SI;
  const double kSiliconMass_SI = 28.0855*kAtomicMassUnit_SI;
  const double kIronMass_SI = 55.845*kAtomicMassUnit_SI;

  const double kSilicateMass_SI = kMagnesiumMass_SI + kIronMass_SI + 
    kSiliconMass_SI + 4.0*kOxygenMass_SI; // Silicate is MgFeSi0_4

  const double kGraphiteDensity_SI = 2.24*(gram/kilogram)/pow(cm/m, 3.0); 
  const double kSilicateDensity_SI = 3.50*(gram/kilogram)/pow(cm/m, 3.0); 
  
  // Particle lifetimes

  const double kMuonLifetime = 2.19703e-6*s;

  const double kNeutronLifetime = 885.7*s;

  const double kLambdaLifetime = 2.632e-10*s;
  const double kSigmaZeroLifetime = 7.4e-20*s;
  const double kSigmaPlusLifetime = 0.8018e-10*s;
  const double kSigmaMinusLifetime = 1.479e-10*s;
  const double kXiZeroLifetime = 2.9e-10*s;
  const double kXiMinusLifetime = 1.639e-10*s;
  const double kOmegaMinusLifetime = 0.821-10*s;

  const double kPiZeroLifetime = 8.4e-17*s;
  const double kPiChargedLifetime = 2.6033e-8*s;
  const double kKaonChargedLifetime = 1.2384e-8*s;

  // Derived constants

  const double kEpsilonZero_SI = 
    1.0/(kMuZero_SI*kSpeedOfLight_SI*kSpeedOfLight_SI);
  const double kAlpha = (e_SI*e_SI)/
    (4.0*kPi*kEpsilonZero_SI*kPlanckReduced_SI*kSpeedOfLight_SI); 
  const double kElectronRadius_SI = (e_SI*e_SI)/
    (4.0*kPi*kEpsilonZero_SI*kElectronMass_SI*
     kSpeedOfLight_SI*kSpeedOfLight_SI);
  const double kThomsonCrossSection_SI = 
    8.0*kPi*kElectronRadius_SI*kElectronRadius_SI/3.0;
      
  // Distance conversions

  const double kParsec = 3.0856775807e+16*m; 
  const double kKiloParsec = kParsec*1.0e+3;
  const double kMegaParsec = kKiloParsec*1.0e+3;

  const double pc = kParsec;
  const double kpc = kKiloParsec;
  const double Mpc = kMegaParsec;

  // Some other conversions and constants

  const double kGalactocentricRadiusSun = 8.5*kpc;

  const double kYearToSec = 365.25*24.0*60.0*60.0;
  const double kSolarLuminosity_SI = 3.846e26;
  const double kSolarLuminosity = kSolarLuminosity_SI*watt;

  // Conversion from LSun/kpc^3 to eV cm^-3 sr^-1 kpc^-1
  // -> LSun (J s^-1) / (kpc/cm)^3 * kpc/cm / (4*Pi*c)
  const double kWainscoatConversion = 
    (kSolarLuminosity_SI/e_SI)/pow(kpc/cm, 3.0)*
    (kpc/kSpeedOfLight_SI)*(1.0/(4.0*kPi));
    
  // Conversion from LSun/pc^3 to eV cm^-3 sr^-1 kpc^-1
  // -> LSun (J s^-1) / (pc/cm)^3 * kpc/cm / (4*Pi*c)
  const double kMathisConversion = 
    (kSolarLuminosity_SI/e_SI)/pow(pc/cm, 3.0)*
    (kpc/kSpeedOfLight_SI)*(1.0/(4.0*kPi));
    
  // Conversion from MJy sr^-1 kpc^-1 -> eV cm^-3 sr^-1 kpc^-1.
  // Take M == 10^6, Jy == 10^-26 W m^-2 Hz^-1 -> 10^-26/e_SI eV m^-2 Hz^-1 
  // 1/c converts cm^-2 s^-1 to cm^-3 if in kpc s^-1
  const double kFreudenreichConversion = 
    1.0e6*(1.0e-26/e_SI)*(cm/kSpeedOfLight_SI)*(cm*cm);
  
  // Conversion of W sr^-1 H-atom^-1 H-atom cm^-3 -> eV cm^-3 sr^-1 kpc^-1
  // W -> J s^-1/e_SI -> eV s^-1
  // 1/c converts cm^-2 s^-1 -> cm^-2 kpc^-1 if c in kpc s^-1
  // So get W sr^-1 cm^-3 -> eV cm^-3 sr^-1 kpc^-1 
  // times factor 10^-32 associated with emissivity from table 

  const double kSodroskiConversion = (1.0e-32/e_SI)*(kpc/kSpeedOfLight_SI);

  // Temperature, absolute magnitude at V, bolometric correction and magnitude 
  // for Sun 

  const double kTemperatureSun = 5770;
  const double kMVSun = 4.82;
  const double kBCSun = -0.05;
  const double kMBolometricSun = kMVSun + kBCSun;

  const double kFrequencyConstant = 
    2.0*kPlanck_SI/pow(kSpeedOfLight_SI, 2.0);
  const double kFrequencyExponentConstant = kPlanck_SI/kBoltzmann_SI;
  const double kWavelengthConstant = 
    2.0*kPlanck_SI*pow(kSpeedOfLight_SI, 2.0);
  const double kWavelengthExponentConstant = 
    kPlanck_SI*kSpeedOfLight_SI/kBoltzmann_SI;
  const double kPlanckIntegral = 2.0*pow(kBoltzmann_SI*kPi, 4.0)
    /(15.0*pow(kPlanck_SI*kSpeedOfLight_SI, 2.0)*kPlanck_SI);
  const double kStefanBoltzmann = kPlanckIntegral*kPi;

}				  

#endif 

  
