import FWCore.ParameterSet.Config as cms

##----------------------------------------------------------------------------------------------------
## CaloMET collections
##----------------------------------------------------------------------------------------------------
#
## Raw CaloMET, no jet smearing
#
#rootTupleCaloMET = cms.EDProducer("RootTupleMakerV2_MET",
#    InputTag = cms.InputTag('patMETsRawCalo'),
#    Prefix = cms.string('Calo'),
#    Suffix = cms.string(''),
#    StoreUncorrectedMET = cms.bool(False),
#    StoreMETSignificance = cms.bool(False),
#    Uncertainty = cms.string('NoShift'),
#    CorrectionLevel = cms.string('Calo'))
#
## Type1 CaloMET, no jet smearing
#
#rootTupleCaloMETType1Cor = cms.EDProducer("RootTupleMakerV2_MET",
#    InputTag = cms.InputTag('patMETs'),
#    Prefix = cms.string('Calo'),
#    Suffix = cms.string('Type1Cor'),
#    StoreUncorrectedMET = cms.bool(True),
#    StoreMETSignificance = cms.bool(False),
#    Uncertainty = cms.string('NoShift'),
#    CorrectionLevel = cms.string('Calo'))
#
#----------------------------------------------------------------------------------------------------
# PFMET collections
#----------------------------------------------------------------------------------------------------

# Raw PFMET, no jet smearing

rootTuplePFMET = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string(''),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('NoShift'),
    CorrectionLevel = cms.string('Raw'))

# Type1 PFMET, with jet smearing

rootTuplePFMETType1Cor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1Cor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('NoShift'),
    CorrectionLevel = cms.string('Type1'))

# Type0+1 PFMET, with jet smearing

rootTuplePFMETType01Cor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01Cor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('NoShift'),
    CorrectionLevel = cms.string('Type01'))

# Type0+1+XY PFMET, with jet smearing

rootTuplePFMETType01XYCor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('NoShift'),
    CorrectionLevel = cms.string('Type01XY')                                       
)

# Raw PFMET, no jet smearing
rootTuplePFMETPuppi = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETsPuppi'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Puppi'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('NoShift'),
    CorrectionLevel = cms.string('Raw'))

# Type1 PFMET, with jet smearing
rootTuplePFMETPuppiType1Cor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETsPuppi'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('PuppiType1Cor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('NoShift'),
    CorrectionLevel = cms.string('Type1'))

#----------------------------------------------------------------------------------------------------
# PFMET systematics collections                                                                
#----------------------------------------------------------------------------------------------------

# Shift unclustered energy up

rootTuplePFMETType01XYCorUnclusteredUp = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorUnclusteredUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('UnclusteredEnUp'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered energy down

rootTuplePFMETType01XYCorUnclusteredDown = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorUnclusteredDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('UnclusteredEnDown'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered electron energy up

rootTuplePFMETType01XYCorElectronEnUp   = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorElectronEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('ElectronEnUp'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered electron energy down

rootTuplePFMETType01XYCorElectronEnDown = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorElectronEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('ElectronEnDown'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered muon energy up

rootTuplePFMETType01XYCorMuonEnUp       = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorMuonEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('MuonEnUp'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered muon energy down

rootTuplePFMETType01XYCorMuonEnDown     = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorMuonEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('MuonEnDown'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered tau energy up

rootTuplePFMETType01XYCorTauEnUp        = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorTauEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('TauEnUp'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered tau energy down

rootTuplePFMETType01XYCorTauEnDown      = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorTauEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('TauEnDown'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered jet energy resolution shifted up

rootTuplePFMETType01XYCorJetResUp       = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorJetResUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('JetResUp'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered jet energy resolution shifted down

rootTuplePFMETType01XYCorJetResDown     = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorJetResDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('JetResDown'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered jet energy shifted up

rootTuplePFMETType01XYCorJetEnUp        = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorJetEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('JetEnUp'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered jet energy shifted down

rootTuplePFMETType01XYCorJetEnDown      = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorJetEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('JetEnDown'),
    CorrectionLevel = cms.string('Type1'))

