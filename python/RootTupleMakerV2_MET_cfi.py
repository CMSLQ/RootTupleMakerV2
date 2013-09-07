import FWCore.ParameterSet.Config as cms

rootTupleCaloMET = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETsRawCalo'),
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
    InputTag = cms.InputTag('patMETsRawPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string(''),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(True)
)

rootTupleCaloMETType1Cor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETs'),
    Prefix = cms.string('Calo'),
    Suffix = cms.string('Type1Cor'),
    StoreUncorrectedMET = cms.bool(True),
    StoreMETSignificance = cms.bool(False)
)

rootTuplePFMETType1Cor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETsPF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type1Cor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(True)
)

rootTuplePFMETType01Cor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETsAK5PF'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01Cor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(True)
)

rootTuplePFMETType01XYCor = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patMETsAK5PFXYShift'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCor'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(True)
)

rootTuplePFMETType01XYCorJetSmeared = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patType1CorrectedPFMet'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorJetSmeared'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(True)
)

rootTuplePFMETType01XYCorJetSmearedUnclusteredUp = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patType1CorrectedPFMetUnclusteredEnUp'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorJetSmearedUnclusteredUp'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

rootTuplePFMETType01XYCorJetSmearedUnclusteredDown = cms.EDProducer("RootTupleMakerV2_MET",
    InputTag = cms.InputTag('patType1CorrectedPFMetUnclusteredEnDown'),
    Prefix = cms.string('PF'),
    Suffix = cms.string('Type01XYCorJetSmearedUnclusteredDown'),
    StoreUncorrectedMET = cms.bool(False),
    StoreMETSignificance = cms.bool(False)
)

