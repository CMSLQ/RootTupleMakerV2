import FWCore.ParameterSet.Config as cms

rootTupleGenParticles = cms.EDProducer("RootTupleMakerV2_GenParticles",
    InputTag = cms.InputTag('genParticles'),
    Prefix = cms.string('GenParticle'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(25)
)

rootTupleGenTausFromLQs = cms.EDProducer("RootTupleMakerV2_GenParticles",
    InputTag = cms.InputTag('genTausFromLQs'),
    Prefix = cms.string('GenLQTau'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(25)
)

rootTupleGenTausFromLQTops = cms.EDProducer("RootTupleMakerV2_GenParticles",
    InputTag = cms.InputTag('genTausFromLQTops'),
    Prefix = cms.string('GenLQTopTau'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(25)
)

