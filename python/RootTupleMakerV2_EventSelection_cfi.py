import FWCore.ParameterSet.Config as cms

rootTupleEventSelection = cms.EDProducer("RootTupleMakerV2_EventSelection",
    L1InputTag  = cms.InputTag('gtDigis'),
    VertexInputTag = cms.InputTag('offlinePrimaryVertices'),
    VertexMinimumNDOF = cms.uint32(4),
    VertexMaxAbsZ = cms.double(24.),
    VertexMaxd0 = cms.double(2.),
    TracksInputTag = cms.InputTag('generalTracks'),
    NumTracks = cms.uint32(10),
    HPTrackThreshold = cms.double(0.25),
    HcalNoiseInputTag = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult')
)
