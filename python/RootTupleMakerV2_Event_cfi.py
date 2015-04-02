import FWCore.ParameterSet.Config as cms

rootTupleEvent = cms.EDProducer("RootTupleMakerV2_Event",
                                #FastJetForIsolationInputTag = cms.InputTag('kt6PFJetsForIsolation','rho'),#not present in 2012
                                #FastJetForJECInputTag       = cms.InputTag('ak4PFJets','rho','RECO'),
                                FastJetForJECInputTag       = cms.InputTag('fixedGridRhoFastjetAll'),
                                FastJetForHEEPInputTag      = cms.InputTag('kt6PFJetsForHEEPIsolation','rho'),
                                FastJetForJECCCPUInputTag   = cms.InputTag('fixedGridRhoFastjetCentralChargedPileUp'),
                                FastJetForJECCNInputTag     = cms.InputTag('fixedGridRhoFastjetCentralNeutral'),
                                FastJetForJECCNTInputTag    = cms.InputTag('fixedGridRhoFastjetCentralNeutral')
                                )


