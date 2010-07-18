import FWCore.ParameterSet.Config as cms

rootTupleTrigger = cms.EDProducer("RootTupleMakerV2_Trigger",
    L1InputTag  = cms.InputTag('gtDigis'),
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    HLTPathsOfInterest = cms.vstring('HLT_Mu9','HLT_Photon15_L1R','HLT_Photon20_L1R','HLT_Photon30_L1R_8E29','HLT_Ele15_LW_L1R','HLT_Jet30U','HLT_MET45','HLT_MinBiasBSC','HLT_Photon15_Cleaned_L1R','HLT_Photon20_Cleaned_L1R')
)
