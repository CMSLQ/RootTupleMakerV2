import FWCore.ParameterSet.Config as cms

rootTupleCaloJets = cms.EDProducer("RootTupleMakerV2_CaloJets",
    InputTag = cms.InputTag('cleanPatJets'),
    Prefix = cms.string('CaloJet'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10),
    ElectronPt = cms.double(30.),
    ElectronIso = cms.double(0.1),
    MuonPt = cms.double(10.),
    MuonIso = cms.double(0.05),
    JECUncertainty = cms.string('CondFormats/JetMETObjects/data/Spring10DataV2_Uncertainty_AK5Calo.txt'),
    ApplyResidualJEC = cms.bool(False),
    ResidualJEC = cms.string('CondFormats/JetMETObjects/data/Spring10DataV2_L2L3Residual_AK5Calo.txt')
)

