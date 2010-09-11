import FWCore.ParameterSet.Config as cms

rootTupleVertex = cms.EDProducer("RootTupleMakerV2_Vertex",
    InputTag = cms.InputTag('offlinePrimaryVertices'),
    Prefix = cms.string('Vertex'),
    Suffix = cms.string('')
)

