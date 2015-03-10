import FWCore.ParameterSet.Config as cms

rootTuplePFCandidates = cms.EDProducer("RootTupleMakerV2_PFCandidates",
                                       JetInputTag = cms.InputTag('slimmedJets'),
                                       ElectronInputTag = cms.InputTag('slimmedElectrons'),
                                       MuonInputTag = cms.InputTag('slimmedMuons'),  
                                       PFCandInputTag = cms.InputTag('packedPFCandidates'),
                                       Prefix = cms.string('PFCand'),
                                       Suffix = cms.string(''),
                                       MaxSize = cms.uint32(10),
                                       DRmatch = cms.double(0.1)
                                       )
