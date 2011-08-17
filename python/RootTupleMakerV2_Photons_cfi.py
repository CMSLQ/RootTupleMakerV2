import FWCore.ParameterSet.Config as cms

rootTuplePhotons = cms.EDProducer("RootTupleMakerV2_Photons",
    InputTag = cms.InputTag('cleanPatPhotons'),
    Prefix = cms.string('Photon'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10),
    BeamSpotInputTag = cms.InputTag('offlineBeamSpot'),
    ConversionsInputTag = cms.InputTag('allConversions'),
    ElectronsInputTag = cms.InputTag('gsfElectrons'),
    EcalRecHitsEBInputTag = cms.InputTag('reducedEcalRecHitsEB'),
    EcalRecHitsEEInputTag = cms.InputTag('reducedEcalRecHitsEE')
)
