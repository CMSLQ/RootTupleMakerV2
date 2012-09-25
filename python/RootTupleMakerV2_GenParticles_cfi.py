import FWCore.ParameterSet.Config as cms

rootTupleGenParticles = cms.EDProducer("RootTupleMakerV2_GenParticles",
    InputTag = cms.InputTag('genParticles'),
    Prefix = cms.string('GenParticle'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(25)
)

rootTupleGenTausFromLQs = cms.EDProducer("RootTupleMakerV2_GenParticles",
    InputTag = cms.InputTag('genTausFromLQs'),
    Prefix = cms.string('GenLQTau'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(25)
)

rootTupleGenTausFromLQTops = cms.EDProducer("RootTupleMakerV2_GenParticles",
    InputTag = cms.InputTag('genTausFromLQTops'),
    Prefix = cms.string('GenLQTopTau'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(25)
)

# -------------------------------------------------------------------------------------------------------------- #
# Optional modules                                                                                               #
# See: Leptoquarks.LeptonJetGenTools.genTauMuElFromLQs_cfi , Leptoquarks.LeptonJetGenTools.genTauMuElFromZs_cfi  #
# -------------------------------------------------------------------------------------------------------------- #

#rootTupleGenMuonsFromLQTaus = cms.EDProducer("RootTupleMakerV2_GenParticles",
#                                             InputTag = cms.InputTag('genMuonsFromLQTaus'),
#                                             Prefix = cms.string('GenLQTauMuon'),
#                                             Suffix = cms.string(''),
#                                             MaxSize = cms.uint32(25)
#                                             )

#rootTupleGenElectronsFromLQTaus = cms.EDProducer("RootTupleMakerV2_GenParticles",
#                                                 InputTag = cms.InputTag('genElectronsFromLQTaus'),
#                                                 Prefix = cms.string('GenLQTauElectron'),
#                                                 Suffix = cms.string(''),
#                                                 MaxSize = cms.uint32(25)
#                                                 )

#rootTupleGenMuonsFromLQTops = cms.EDProducer("RootTupleMakerV2_GenParticles",
#                                             InputTag = cms.InputTag('genMuonsFromLQTops'),
#                                             Prefix = cms.string('GenLQTopXMuon'),
#                                             Suffix = cms.string(''),
#                                             MaxSize = cms.uint32(25)
#                                             )

#rootTupleGenElectronsFromLQTops = cms.EDProducer("RootTupleMakerV2_GenParticles",
#                                                 InputTag = cms.InputTag('genElectronsFromLQTops'),
#                                                 Prefix = cms.string('GenLQTopXElectron'),
#                                                 Suffix = cms.string(''),
#                                                 MaxSize = cms.uint32(25)
#                                                 )

#rootTupleGenTausFromZs = cms.EDProducer("RootTupleMakerV2_GenParticles",
#                                        InputTag = cms.InputTag('genTausFromZs'),
#                                        Prefix = cms.string('GenZTau'),
#                                        Suffix = cms.string(''),
#                                        MaxSize = cms.uint32(25)
#                                        )

#rootTupleGenMuonsFromZs = cms.EDProducer("RootTupleMakerV2_GenParticles",
#                                         InputTag = cms.InputTag('genMuonsFromZs'),
#                                         Prefix = cms.string('GenZMu'),
#                                         Suffix = cms.string(''),
#                                         MaxSize = cms.uint32(25)
#                                         )

#rootTupleGenElectronsFromZs = cms.EDProducer("RootTupleMakerV2_GenParticles",
#                                             InputTag = cms.InputTag('genElectronsFromZs'),
#                                             Prefix = cms.string('GenZElectron'),
#                                             Suffix = cms.string(''),
#                                             MaxSize = cms.uint32(25)
#                                             )
