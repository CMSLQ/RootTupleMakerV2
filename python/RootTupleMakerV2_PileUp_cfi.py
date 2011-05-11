import FWCore.ParameterSet.Config as cms

rootTuplePileUp = cms.EDProducer("RootTupleMakerV2_PileUp",
    pileupInfo = cms.InputTag('addPileupInfo')
)

