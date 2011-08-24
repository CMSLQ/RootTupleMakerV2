import FWCore.ParameterSet.Config as cms

rootTupleCaloMET = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETs'),
    Prefix = cms.string('Calo'),
    Suffix = cms.string(''),
    StoreUncorrectedMET = cms.bool(True),
    StoreMETSignificance = cms.bool(False)
)

rootTupleTCMET = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETsTC'),
    Prefix = cms.string('TC'),
    Suffix = cms.string(''),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

rootTuplePFMET = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETsPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string(''),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(True)
)

rootTuplePFMETType1Cor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETsPFType1Cor'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1Cor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

rootTuplePFChargedMET = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETsPFCharged'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Charged'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)
