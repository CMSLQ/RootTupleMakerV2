import FWCore.ParameterSet.Config as cms

#----------------------------------------------------------------------------------------------------
# TCMET collections
#----------------------------------------------------------------------------------------------------

# Raw TCMET, no jet smearing

rootTupleTCMET = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETsTC'),
    Prefix = cms.string('TC'),
    Suffix = cms.string(''),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string(''),
    Level = cms.string(''))

#----------------------------------------------------------------------------------------------------
# CaloMET collections
#----------------------------------------------------------------------------------------------------

# Raw CaloMET, no jet smearing

rootTupleCaloMET = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETsRawCalo'),
    Prefix = cms.string('Calo'),
    Suffix = cms.string(''),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('NoShift'),
    Level = cms.string('Calo'))

# Type1 CaloMET, no jet smearing

rootTupleCaloMETType1Cor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETs'),
    Prefix = cms.string('Calo'),
    Suffix = cms.string('Type1Cor'),
    StoreUncorrectedMET = cms.bool(True),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('NoShift'),
    Level = cms.string('Calo'))

#----------------------------------------------------------------------------------------------------
# PFMET collections
#----------------------------------------------------------------------------------------------------

# Raw PFMET, no jet smearing

rootTuplePFMET = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    #InputTag = cms.InputTag('patMETsRawPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string(''),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('NoShift'),
    Level = cms.string('Raw'))

# Type1 PFMET, with jet smearing

rootTuplePFMETType1Cor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    #InputTag = cms.InputTag('patType1CorrectedPFMetType1Only'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1Cor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('NoShift'),
    Level = cms.string('Type1'))

# Type0+1 PFMET, with jet smearing

rootTuplePFMETType01Cor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    #InputTag = cms.InputTag('patType1CorrectedPFMetType01Only'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01Cor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('NoShift'),
    Level = cms.string('Type01'))

# Type0+1+XY PFMET, with jet smearing

rootTuplePFMETType01XYCor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    #InputTag = cms.InputTag('patType1CorrectedPFMet'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('NoShift'),
    Level = cms.string('Type01XY')                                       
)

#----------------------------------------------------------------------------------------------------
# PFMET systematics collections                                                                
#----------------------------------------------------------------------------------------------------

# Shift unclustered energy up

rootTuplePFMETType01XYCorUnclusteredUp = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    #InputTag = cms.InputTag('patType1CorrectedPFMetUnclusteredEnUp'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorUnclusteredUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('UnclusteredEnUp'),
    Level = cms.string('Type01XY'))

# Shift unclustered energy down

rootTuplePFMETType01XYCorUnclusteredDown = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    #InputTag = cms.InputTag('patType1CorrectedPFMetUnclusteredEnDown'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorUnclusteredDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('UnclusteredEnDown'),
    Level = cms.string('Type01XY'))

# Shift unclustered electron energy up

rootTuplePFMETType01XYCorElectronEnUp   = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    #InputTag = cms.InputTag('patType1CorrectedPFMetElectronEnUp'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorElectronEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('ElectronEnUp'),
    Level = cms.string('Type01XY'))

# Shift unclustered electron energy down

rootTuplePFMETType01XYCorElectronEnDown = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    #InputTag = cms.InputTag('patType1CorrectedPFMetElectronEnDown'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorElectronEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('ElectronEnDown'),
    Level = cms.string('Type01XY'))

# Shift unclustered muon energy up

rootTuplePFMETType01XYCorMuonEnUp       = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    #InputTag = cms.InputTag('patType1CorrectedPFMetMuonEnUp'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorMuonEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('MuonEnUp'),
    Level = cms.string('Type01XY'))

# Shift unclustered muon energy down

rootTuplePFMETType01XYCorMuonEnDown     = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    #InputTag = cms.InputTag('patType1CorrectedPFMetMuonEnDown'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorMuonEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('MuonEnDown'),
    Level = cms.string('Type01XY'))

# Shift unclustered tau energy up

rootTuplePFMETType01XYCorTauEnUp        = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    #InputTag = cms.InputTag('patType1CorrectedPFMetTauEnUp'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorTauEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('TauEnUp'),
    Level = cms.string('Type01XY'))

# Shift unclustered tau energy down

rootTuplePFMETType01XYCorTauEnDown      = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    #InputTag = cms.InputTag('patType1CorrectedPFMetTauEnDown'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorTauEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('TauEnDown'),
    Level = cms.string('Type01XY'))

# Shift unclustered jet energy resolution shifted up

rootTuplePFMETType01XYCorJetResUp       = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    #InputTag = cms.InputTag('patType1CorrectedPFMetJetResUp'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorJetResUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('JetResUp'),
    Level = cms.string('Type01XY'))

# Shift unclustered jet energy resolution shifted down

rootTuplePFMETType01XYCorJetResDown     = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    #InputTag = cms.InputTag('patType1CorrectedPFMetJetResDown'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorJetResDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('JetResDown'),
    Level = cms.string('Type01XY'))

# Shift unclustered jet energy shifted up

rootTuplePFMETType01XYCorJetEnUp        = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    #InputTag = cms.InputTag('patType1CorrectedPFMetJetEnUp'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorJetEnUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('JetEnUp'),
    Level = cms.string('Type01XY'))

# Shift unclustered jet energy shifted down

rootTuplePFMETType01XYCorJetEnDown      = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('slimmedMETs'),
    #InputTag = cms.InputTag('patType1CorrectedPFMetJetEnDown'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorJetEnDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False),
    Uncertainty = cms.string('JetEnDown'),
    Level = cms.string('Type01XY'))

