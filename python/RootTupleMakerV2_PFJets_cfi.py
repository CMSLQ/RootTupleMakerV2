import FWCore.ParameterSet.Config as cms

rootTuplePFJets = cms.EDProducer("RootTupleMakerV2_PFJets",
    InputTag = cms.InputTag('selectedPatJetsAK5PF'),
    InputTagL1Offset = cms.InputTag('selectedPatJetsAK5PFL1Offset'),
    Prefix = cms.string('PFJet'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(15),
    JECUncertainty = cms.string('AK5PF'),
    ReadJECuncertainty = cms.bool(True),
    VertexInputTag = cms.InputTag('offlinePrimaryVertices')
 
    #OLD                                
    #JECUncertainty = cms.string('CondFormats/JetMETObjects/data/Spring10_Uncertainty_AK5PF.txt'),
    #ApplyResidualJEC = cms.bool(False),
    #ResidualJEC = cms.string('CondFormats/JetMETObjects/data/Spring10_L2L3Residual_AK5PF.txt')
)

