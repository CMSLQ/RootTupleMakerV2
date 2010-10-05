import FWCore.ParameterSet.Config as cms

rootTupleSuperClusters = cms.EDProducer("RootTupleMakerV2_SuperClusters",
#    EBInputTag = cms.InputTag('hybridSuperClusters'),
#    EEInputTag = cms.InputTag('multi5x5SuperClustersWithPreshower'),
    EBInputTag = cms.InputTag('correctedHybridSuperClusters'),
    EEInputTag = cms.InputTag('correctedMulti5x5SuperClustersWithPreshower'),
    Prefix = cms.string('SuperCluster'),
    Suffix = cms.string(''),
    EBMaxSize = cms.uint32(30),
    EEMaxSize = cms.uint32(30),
    EcalEBInputTag = cms.InputTag('reducedEcalRecHitsEB'),
    EcalEEInputTag = cms.InputTag('reducedEcalRecHitsEE'),
    TracksInputTag = cms.InputTag('generalTracks'),
    ElectronsInputTag = cms.InputTag('cleanPatElectrons'),
    ElectronsPrefix = cms.string('Electron'),
    ElectronsMaxSize = cms.uint32(10)
)
