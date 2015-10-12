import FWCore.ParameterSet.Config as cms

rootTupleEvent = cms.EDProducer("RootTupleMakerV2_Event",
   globalTag                                       = cms.string(''),
   FixedGridRhoAllInputTag                         = cms.InputTag('fixedGridRhoAll'),
   FixedGridRhoFastjetAllInputTag                  = cms.InputTag('fixedGridRhoFastjetAll'),
   FixedGridRhoFastjetAllCaloInputTag              = cms.InputTag('fixedGridRhoFastjetAllCalo'),
   FixedGridRhoFastjetCentralCaloInputTag          = cms.InputTag('fixedGridRhoFastjetCentralCalo'),
   FixedGridRhoFastjetCentralChargedPileUpInputTag = cms.InputTag('fixedGridRhoFastjetCentralChargedPileUp'),
   FixedGridRhoFastjetCentralNeutralInputTag       = cms.InputTag('fixedGridRhoFastjetCentralNeutral'),
)
