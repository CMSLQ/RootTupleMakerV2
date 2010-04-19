import FWCore.ParameterSet.Config as cms

rootTupleCaloJets = cms.EDProducer("RootTupleMakerV2_CaloJets",
    InputTag = cms.InputTag('cleanPatJets'),
    Prefix = cms.string('CaloJet'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10),
    ElectronPt = cms.double(30.),
    ElectronIso = cms.double(0.1),
    MuonPt = cms.double(10.),
    MuonIso = cms.double(0.05)
)

