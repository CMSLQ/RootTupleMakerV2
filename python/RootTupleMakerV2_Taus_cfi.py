import FWCore.ParameterSet.Config as cms

rootTupleHPSTaus = cms.EDProducer("RootTupleMakerV2_Taus",
                                  InputTag = cms.InputTag('cleanPatTaus'),
                                  Prefix   = cms.string('HPSTau'),
                                  Suffix   = cms.string(''),
                                  VertexInputTag = cms.InputTag('offlinePrimaryVertices'),
                                  isSCTau  = cms.bool(False),
                                  isHPSTau = cms.bool(True),
                                  MaxSize  = cms.uint32(50)
                                  )
