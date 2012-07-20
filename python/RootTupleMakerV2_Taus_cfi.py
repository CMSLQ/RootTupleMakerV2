import FWCore.ParameterSet.Config as cms

#rootTupleTaus = cms.EDProducer("RootTupleMakerV2_Taus",
#                               InputTag = cms.InputTag('cleanPatTaus'),
#                               Prefix = cms.string('Tau'),
#                               Suffix = cms.string(''),
#                               MaxSize = cms.uint32(50),
#                               isSCTau  = cms.bool(False),
#                               isHPSTau = cms.bool(True)
#)

rootTupleSCTaus = cms.EDProducer("RootTupleMakerV2_Taus",
                                 InputTag = cms.InputTag('cleanPatTausShrinkingConePFTau'),
                                 Prefix   = cms.string('SCTau'),
                                 Suffix   = cms.string(''),
                                 isSCTau  = cms.bool(True),
                                 isHPSTau = cms.bool(False),
                                 MaxSize  = cms.uint32(50)
                                 )

rootTupleHPSTaus = cms.EDProducer("RootTupleMakerV2_Taus",
                                  #InputTag = cms.InputTag('cleanPatTausHpsPFTau'),
                                  InputTag = cms.InputTag('cleanPatTaus'),
                                  Prefix   = cms.string('HPSTau'),
                                  Suffix   = cms.string(''),
                                  isSCTau  = cms.bool(False),
                                  isHPSTau = cms.bool(True),
                                  MaxSize  = cms.uint32(50)
                                  )
