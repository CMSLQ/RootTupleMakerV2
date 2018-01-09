import FWCore.ParameterSet.Config as cms

#----------------------------------------------------------------------------------------------------
# PFMET collections
#----------------------------------------------------------------------------------------------------

# Raw PFMET, no jet smearing

rootTuplePFMETRaw = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Raw'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(True),
    Uncertainty = cms.string('NoShift'),
    CorrectionLevel = cms.string('Raw'))

# Type1 PFMET, with jet smearing

rootTuplePFMETType1Cor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1Cor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(True),
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

# Type1+XY PFMET, with jet smearing

rootTuplePFMETType1XYCor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1XYCor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(True),
    Uncertainty = cms.string('NoShift'),
    CorrectionLevel = cms.string('Type1XY')) 

# Type0+1+XY PFMET, with jet smearing

rootTuplePFMETType01XYCor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('NoShift'),
    CorrectionLevel = cms.string('Type01XY'))                                       

# Raw PuppiMET, no jet smearing
rootTuplePFMETPuppiRaw = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETsPuppi'),
    Prefix = cms.string('Puppi'),
    Suffix = cms.string('Raw'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('NoShift'),
    CorrectionLevel = cms.string('Raw'))

# Type1 PuppiMET, with jet smearing
rootTuplePFMETPuppiType1Cor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETsPuppi'),
    Prefix = cms.string('Puppi'),
    Suffix = cms.string('Type1Cor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('NoShift'),
    CorrectionLevel = cms.string('Type1'))

#----------------------------------------------------------------------------------------------------
# PFMET systematics collections                                                                
#----------------------------------------------------------------------------------------------------

# Shift unclustered energy up

rootTuplePFMETType1CorUnclusteredUp = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorUnclusteredUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('UnclusteredEnUp'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered energy down

rootTuplePFMETType1CorUnclusteredDown = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorUnclusteredDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('UnclusteredEnDown'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered electron energy up

rootTuplePFMETType1CorElectronEnUp   = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorElectronEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('ElectronEnUp'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered electron energy down

rootTuplePFMETType1CorElectronEnDown = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorElectronEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('ElectronEnDown'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered muon energy up

rootTuplePFMETType1CorMuonEnUp       = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorMuonEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('MuonEnUp'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered muon energy down

rootTuplePFMETType1CorMuonEnDown     = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorMuonEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('MuonEnDown'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered tau energy up

rootTuplePFMETType1CorTauEnUp        = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorTauEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('TauEnUp'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered tau energy down

rootTuplePFMETType1CorTauEnDown      = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorTauEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('TauEnDown'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered jet energy resolution shifted up

rootTuplePFMETType1CorJetResUp       = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorJetResUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('JetResUp'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered jet energy resolution shifted down

rootTuplePFMETType1CorJetResDown     = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorJetResDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('JetResDown'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered jet energy shifted up

rootTuplePFMETType1CorJetEnUp        = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorJetEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('JetEnUp'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered jet energy shifted down

rootTuplePFMETType1CorJetEnDown      = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1CorJetEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('JetEnDown'),
    CorrectionLevel = cms.string('Type1'))

# Shift unclustered energy up

rootTuplePFMETType1XYCorUnclusteredUp = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1XYCorUnclusteredUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('UnclusteredEnUp'),
    CorrectionLevel = cms.string('Type1XY'))

# Shift unclustered energy down

rootTuplePFMETType1XYCorUnclusteredDown = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1XYCorUnclusteredDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('UnclusteredEnDown'),
    CorrectionLevel = cms.string('Type1XY'))

# Shift unclustered electron energy up

rootTuplePFMETType1XYCorElectronEnUp   = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1XYCorElectronEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('ElectronEnUp'),
    CorrectionLevel = cms.string('Type1XY'))

# Shift unclustered electron energy down

rootTuplePFMETType1XYCorElectronEnDown = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1XYCorElectronEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('ElectronEnDown'),
    CorrectionLevel = cms.string('Type1XY'))

# Shift unclustered muon energy up

rootTuplePFMETType1XYCorMuonEnUp       = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1XYCorMuonEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('MuonEnUp'),
    CorrectionLevel = cms.string('Type1XY'))

# Shift unclustered muon energy down

rootTuplePFMETType1XYCorMuonEnDown     = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1XYCorMuonEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('MuonEnDown'),
    CorrectionLevel = cms.string('Type1XY'))

# Shift unclustered tau energy up

rootTuplePFMETType1XYCorTauEnUp        = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1XYCorTauEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('TauEnUp'),
    CorrectionLevel = cms.string('Type1XY'))

# Shift unclustered tau energy down

rootTuplePFMETType1XYCorTauEnDown      = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1XYCorTauEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('TauEnDown'),
    CorrectionLevel = cms.string('Type1XY'))

# Shift unclustered jet energy resolution shifted up

rootTuplePFMETType1XYCorJetResUp       = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1XYCorJetResUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('JetResUp'),
    CorrectionLevel = cms.string('Type1XY'))

# Shift unclustered jet energy resolution shifted down

rootTuplePFMETType1XYCorJetResDown     = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1XYCorJetResDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('JetResDown'),
    CorrectionLevel = cms.string('Type1XY'))

# Shift unclustered jet energy shifted up

rootTuplePFMETType1XYCorJetEnUp        = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1XYCorJetEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('JetEnUp'),
    CorrectionLevel = cms.string('Type1XY'))

# Shift unclustered jet energy shifted down

rootTuplePFMETType1XYCorJetEnDown      = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1XYCorJetEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('JetEnDown'),
    CorrectionLevel = cms.string('Type1XY'))

rootTuplePFMETSequence = cms.Sequence(
    rootTuplePFMETRaw*
    rootTuplePFMETType1Cor*
    #rootTuplePFMETType01Cor*
    rootTuplePFMETType1XYCor*
    #rootTuplePFMETType01XYCor*
    #rootTuplePFMETPuppiRaw*
    #rootTuplePFMETPuppiType1Cor*
    #rootTuplePFMETType1CorUnclusteredUp*
    #rootTuplePFMETType1CorUnclusteredDown*
    #rootTuplePFMETType1CorElectronEnUp*
    #rootTuplePFMETType1CorElectronEnDown*
    #rootTuplePFMETType1CorMuonEnUp*
    #rootTuplePFMETType1CorMuonEnDown*
    #rootTuplePFMETType1CorTauEnUp*
    #rootTuplePFMETType1CorTauEnDown*
    #rootTuplePFMETType1CorJetResUp*
    #rootTuplePFMETType1CorJetResDown*
    #rootTuplePFMETType1CorJetEnUp*
    #rootTuplePFMETType1CorJetEnDown*
    rootTuplePFMETType1XYCorUnclusteredUp*
    rootTuplePFMETType1XYCorUnclusteredDown*
    rootTuplePFMETType1XYCorElectronEnUp*
    rootTuplePFMETType1XYCorElectronEnDown*
    rootTuplePFMETType1XYCorMuonEnUp*
    rootTuplePFMETType1XYCorMuonEnDown*
    rootTuplePFMETType1XYCorTauEnUp*
    rootTuplePFMETType1XYCorTauEnDown*
    rootTuplePFMETType1XYCorJetResUp*
    rootTuplePFMETType1XYCorJetResDown*
    rootTuplePFMETType1XYCorJetEnUp*
    rootTuplePFMETType1XYCorJetEnDown
)

    
