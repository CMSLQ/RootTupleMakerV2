import FWCore.ParameterSet.Config as cms

rootTuplePFJets = cms.EDProducer("RootTupleMakerV2_PFJets",
    InputTag = cms.InputTag('selectedPatJetsAK5PF'),
    Prefix = cms.string('PFJet'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10),
    ApplyResidualJEC = cms.bool(False),
    ResidualJEC = cms.string('CondFormats/JetMETObjects/data/Spring10DataV1_L2L3Residual_AK5PF.txt')
)

