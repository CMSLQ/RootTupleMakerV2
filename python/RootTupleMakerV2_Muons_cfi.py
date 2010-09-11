import FWCore.ParameterSet.Config as cms

rootTupleMuons = cms.EDProducer("RootTupleMakerV2_Muons",
    InputTag = cms.InputTag('cleanPatMuons'),
    Prefix = cms.string('Muon'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10),
    MuonIso = cms.double(0.05),
    MuonID = cms.string('GlobalMuonPromptTight'),
    BeamSpotCorr = cms.bool(False),
    VertexInputTag = cms.InputTag('offlinePrimaryVertices')
)
