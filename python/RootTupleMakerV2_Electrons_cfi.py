import FWCore.ParameterSet.Config as cms

rootTupleElectrons = cms.EDProducer("RootTupleMakerV2_Electrons",
    InputTag = cms.InputTag('cleanPatElectrons'),
    Prefix = cms.string('Electron'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10),
    ElectronIso = cms.double(0.1),
    MuonPt = cms.double(10.),
    MuonIso = cms.double(0.05),
    MuonID = cms.string('GlobalMuonPromptTight')
)
