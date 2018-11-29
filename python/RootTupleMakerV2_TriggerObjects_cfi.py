import FWCore.ParameterSet.Config as cms
from Leptoquarks.RootTupleMakerV2.RootTupleMakerV2_Trigger_cfi import *

rootTupleTriggerObjects= cms.EDProducer("RootTupleMakerV2_TriggerObjects",
    TriggerBitsTag = cms.InputTag("TriggerResults","","HLT"),
    TriggerObjectsTag = cms.InputTag("unpackedPatTrigger"),
    Prefix = cms.string("HL"),
    Suffix = cms.string(''),
    HLTPathsOfInterest = rootTupleTrigger.HLTPathsOfInterest
)

