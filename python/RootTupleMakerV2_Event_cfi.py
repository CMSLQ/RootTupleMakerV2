import FWCore.ParameterSet.Config as cms

rootTupleEvent = cms.EDProducer("RootTupleMakerV2_Event",
   globalTag                                       = cms.string(''),
   FixedGridRhoAllInputTag                         = cms.InputTag('fixedGridRhoAll'),
   FixedGridRhoFastjetAllInputTag                  = cms.InputTag('fixedGridRhoFastjetAll'),
   FixedGridRhoFastjetAllCaloInputTag              = cms.InputTag('fixedGridRhoFastjetAllCalo'),
   FixedGridRhoFastjetCentralInputTag              = cms.InputTag('fixedGridRhoFastjetCentral'),
   FixedGridRhoFastjetCentralCaloInputTag          = cms.InputTag('fixedGridRhoFastjetCentralCalo'),
   FixedGridRhoFastjetCentralNeutralInputTag       = cms.InputTag('fixedGridRhoFastjetCentralNeutral'),
   FixedGridRhoFastjetCentralChargedPileUpInputTag = cms.InputTag('fixedGridRhoFastjetCentralChargedPileUp'),
)
