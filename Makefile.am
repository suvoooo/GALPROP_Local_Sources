###########################################################################
# Global part
###########################################################################

include $(top_srcdir)/Makefile.am.common

AM_CPPFLAGS += -I$(top_builddir)/include @CFITSIO_CPPFLAGS@ \
@CCFITS_CPPFLAGS@ @HEALPIX_CPPFLAGS@ @CLHEP_CPPFLAGS@ @GSL_CPPFLAGS@\
-DGALDEF_PATH=\"$(GALDEF_PATH)\" -DFITSDATA_PATH=\"$(FITSDATA_PATH)\" \
-DDATA_PATH=\"@datadir@\"

LIBS += @CFITSIO_LIBS@ @CCFITS_LIBS@ @HEALPIX_LIBS@ @FLIBS@ \
@GSL_LIBS@ @CLHEP_LDFLAGS@

bin_PROGRAMS = galprop changeGasMaps
lib_LIBRARIES = libgalprop.a libskymap.a

libskymap_a_SOURCES = Coordinate.cc HealpixBaseExtended.cc Region.cc SkySelection.cc Skymap.cc

libskymap_a_Headers = Skymap.h Coordinate.h HealpixBaseExtended.h Region.h SkySelection.h ArraySlice.h

libskymap_a_includedir=$(includedir)
libskymap_a_include_HEADERS = Skymap.h Coordinate.h HealpixBaseExtended.h Region.h SkySelection.h ArraySlice.h

libgalprop_a_SOURCES = B_field_model.cc Configure.cc D_pp.cc D_xx.cc \
Distribution.cc GalacticRadiationField.cc RadiationField.cc Galaxy.cc \
Galdef.cc Galprop.cc He_to_H_CS.cc IC_anisotropy_factor.cc Particle.cc \
WNEWTR_FUNC_aws.f YIELDX_011000_imos.f antiproton.f blattnig.cc bremss_spec.f \
cfactor.f cr_luminosity.cc create_SNR.cc create_galaxy.cc create_gcr.cc \
create_transport_arrays.cc crn6.f decayed_cross_sections.cc e_KN_loss.cc \
e_loss_compton.f electrons_normalize.cc energy_losses.cc eprop.cc \
fort_interface1.cc fort_interface2.cc gauss.cc gen_DM_skymap.cc \
gen_DM_source.cc gen_IC_emiss.cc gen_IC_skymap.cc gen_bremss_emiss.cc \
gen_bremss_ionized_skymap.cc gen_bremss_skymap.cc gen_ionization_rate.cc \
gen_isrf_energy_density.cc gen_knock_on_electron_source.cc \
gen_pi0_decay_emiss.cc gen_pi0_decay_skymap.cc \
gen_secondary_antiproton_source.cc gen_secondary_positron_source.cc \
gen_secondary_proton_source.cc gen_secondary_source.cc gen_synch_emiss.cc \
gen_synch_skymap.cc gen_tertiary_antiproton_source.cc integ.cc inter.f \
ionization_bethe.cc isrf_energy_density.cc kamae.cc Kcapture_cs.cc \
kinematic.cc nH.cc nuc_package.cc nuclei_normalize.cc nucleon_cs.cc \
pp_meson.f print_BC.cc propagate_particles.cc propel.cc propel_diagnostics.cc \
protri.cc read_gcr.cc read_isrf.cc read_gas_maps.cc SkymapFitsio.cc \
sigma_boron_dec_heinbach_simon.cc source_SNR_event.cc source_SNR_event_vec.cc \
source_distribution.cc store_DM_emiss.cc store_DM_skymap.cc \
store_IC_skymap.cc store_IC_skymap_comp.cc store_bremss_H2R_skymap.cc \
store_bremss_HIR_skymap.cc store_bremss_emiss.cc \
store_bremss_ionized_skymap.cc store_bremss_skymap.cc \
store_bremss_HII_skymap.cc store_pi0_decay_HII_skymap.cc \
store_gcr.cc store_gcr_full.cc store_ionization_rate.cc \
store_mapcube_skymap.cc store_pi0_decay_H2R_skymap.cc \
store_pi0_decay_HIR_skymap.cc store_IC_emiss.cc \
store_pi0_decay_emiss.cc store_pi0_decay_skymap.cc \
store_skymap.cc store_synch_skymap.cc synchrotron.f test_Distribution.cc \
test_Particle.cc test_float_accuracy.cc test_isotope_cs.cc \
test_nH.cc test_source_SNR_event.cc test_suite.cc tridag.cc tridag_ext.cc \
tridag_sym.cc tridag_sym_ext.cc ErrorLogger.cc GCR_data.cc \
calculate_GCR_chisq.cc B_field_3D_model.cc synchrotron_emissivity.cc \
synchrotron_emissivity_B_field.cc synchrotron_emissivity_aws.cc \
test_sync_package.cc store_synch_emiss.cc galprop_internal.cc \
gen_luminosity.cc Timer.cc store_gcr_source_functions.cc \
DistributionFunction.cc GasFunction.cc GasEmissFunction.cc los_integration.cc \
AvXCO.cc

libgalprop_a_Headers = Configure.h Distribution.h Galaxy.h Galdef.h Galprop.h \
Particle.h Spectrum.h constants.h \
fort_interface.h fort_interface1.h fort_interface2.h \
galprop_internal.h galprop_classes.h \
Units.h PhysicalConstants.h GalacticRadiationField.h RadiationField.h \
GCR_data.h synchrotron_emissivity_aws.h \
synchrotron_emissivity_B_field.h synchrotron_emissivity.h B_field_3D_model.h \
Singleton.h Timer.h DistributionFunction.h los_integration.h integ.h \
AmrData.h StellarGeometryDistributions.h SkymapFitsio.h ErrorLogger.h

galprop_SOURCES = galpropMain.cc 
galprop_LDADD = libgalprop.a libskymap.a
changeGasMaps_SOURCES = changeGasMaps.cc
changeGasMaps_LDADD = libgalprop.a libskymap.a

libgalprop_a_includedir=$(includedir)
libgalprop_a_include_HEADERS = Configure.h Distribution.h Galaxy.h Galdef.h \
Galprop.h Particle.h Spectrum.h constants.h \
fort_interface.h fort_interface1.h fort_interface2.h galprop_internal.h \
galprop_classes.h Units.h PhysicalConstants.h GalacticRadiationField.h \
RadiationField.h SkymapFitsio.h \
ArraySlice.h GCR_data.h \
synchrotron_emissivity_aws.h synchrotron_emissivity_B_field.h \
synchrotron_emissivity.h B_field_3D_model.h Singleton.h Timer.h \
DistributionFunction.h los_integration.h integ.h \
ArmData.h StellarGeometryDistributions.h ErrorLogger.h

EXTRA_DIST += barpol.dat eval_iso_cs.dat isotope_cs.dat j0zero500.dat \
nucdata.dat p_cs_fits.dat WNEWTR_082693.CDR.dat GCR_data_1.dat
