import FWCore.ParameterSet.Config as cms

rootTupleGenJetsAK5 = cms.EDProducer("RootTupleMakerV2_GenJets",
    InputTag = cms.InputTag('ak5GenJetsNoNu'),
    Prefix = cms.string('GenJet'),
    Suffix = cms.string('AK5'),
    MaxSize = cms.uint32(10)
)

rootTupleGenJetsAK4 = cms.EDProducer("RootTupleMakerV2_GenJets",
    InputTag = cms.InputTag('slimmedGenJets'),
    Prefix = cms.string('GenJet'),
    Suffix = cms.string('AK4'),
    MaxSize = cms.uint32(10)
)

rootTupleGenJetsSequence = cms.Sequence(rootTupleGenJetsAK5+rootTupleGenJetsAK4)
