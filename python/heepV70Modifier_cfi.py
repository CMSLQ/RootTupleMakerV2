import FWCore.ParameterSet.Config as cms
heep_modifications = cms.VPSet(
    cms.PSet( modifierName    = cms.string('EGExtraInfoModifierFromFloatValueMaps'),
              electron_config = cms.PSet( trkPtIso = cms.InputTag("heepIDVarValueMaps","eleTrkPtIso"),
                                          electronSrc=cms.InputTag("slimmedElectrons",processName=cms.InputTag.skipCurrentProcess()),
                                          ),
              photon_config   = cms.PSet( ),
              ),
    cms.PSet( modifierName    = cms.string('EGExtraInfoModifierFromIntValueMaps'),
              electron_config = cms.PSet( nrSatCrys = cms.InputTag("heepIDVarValueMaps","eleNrSaturateIn5x5"),
                                          electronSrc=cms.InputTag("slimmedElectrons",processName=cms.InputTag.skipCurrentProcess()),
              photon_config   = cms.PSet( )
              ),
              ),
    cms.PSet( modifierName    = cms.string('EGExtraInfoModifierFromBoolToIntValueMaps'),
              electron_config = cms.PSet( heepElectronID_HEEPV70 = cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70"),
                                          electronSrc=cms.InputTag("slimmedElectrons",processName=cms.InputTag.skipCurrentProcess()),
              photon_config   = cms.PSet( )
              ),
              ),
    cms.PSet( modifierName    = cms.string('EGExtraInfoModifierFromUIntToIntValueMaps'),
              electron_config = cms.PSet( heepElectronID_HEEPV70Bitmap = cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70Bitmap"),
                                          electronSrc=cms.InputTag("slimmedElectrons",processName=cms.InputTag.skipCurrentProcess()),
              photon_config   = cms.PSet( )
              ),
              ),

    cms.PSet( modifierName    = cms.string('EGExtraInfoModifierFromVIDCutFlowResultValueMaps'),
              electron_config = cms.PSet( heepElectronID_HEEPV70 = cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV70"),
                                          electronSrc=cms.InputTag("slimmedElectrons",processName=cms.InputTag.skipCurrentProcess()),
                                          ),
              photon_config   = cms.PSet( )
              )
)
              
