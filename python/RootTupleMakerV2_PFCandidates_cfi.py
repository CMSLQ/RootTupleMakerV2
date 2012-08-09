import FWCore.ParameterSet.Config as cms

rootTuplePFCandidates = cms.EDProducer("RootTupleMakerV2_PFCandidates",
                                       ReducedPFCandidateInputTag = cms.InputTag('pfCandsNotInJet'),
                                       ElectronInputTag = cms.InputTag('cleanPatElectrons'),
                                       MuonInputTag = cms.InputTag('cleanPatMuons'),                                    
                                       Prefix = cms.string('PFCand'),
                                       Suffix = cms.string(''),
                                       MaxSize = cms.uint32(10),
                                       DRmatch = cms.double(0.1)
                                       )
