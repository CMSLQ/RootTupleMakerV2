import FWCore.ParameterSet.Config as cms

rootTupleEvent = cms.EDProducer("RootTupleMakerV2_Event",
                                FastJetForIsolationInputTag = cms.InputTag('kt6PFJetsForIsolation','rho')
                                )

