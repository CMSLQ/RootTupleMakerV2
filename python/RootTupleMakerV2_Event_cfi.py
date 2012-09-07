import FWCore.ParameterSet.Config as cms

rootTupleEvent = cms.EDProducer("RootTupleMakerV2_Event",
                                #FastJetForIsolationInputTag = cms.InputTag('kt6PFJetsForIsolation','rho'),#not present in 2012
                                FastJetForJECInputTag       = cms.InputTag('kt6PFJets','rho'),
                                FastJetForHEEPInputTag      = cms.InputTag('kt6PFJetsForHEEPIsolation','rho'),
                                FastJetForJECCCPUInputTag   = cms.InputTag('kt6PFJetsCentralChargedPileUp','rho'),
                                FastJetForJECCNInputTag     = cms.InputTag('kt6PFJetsCentralNeutral','rho'),
                                FastJetForJECCNTInputTag    = cms.InputTag('kt6PFJetsCentralNeutralTight','rho')
                                )

