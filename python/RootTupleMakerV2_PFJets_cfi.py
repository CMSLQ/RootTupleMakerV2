import FWCore.ParameterSet.Config as cms

rootTuplePFJetsAK5 = cms.EDProducer("RootTupleMakerV2_PFJets",
    InputTag = cms.InputTag('patJetsAK5PF'),
    #InputTag = cms.InputTag('pfjetTriggerMatchEmbedderHLTEleJetJet'),
    ## InputTagL1Offset    = cms.InputTag('selectedPatJetsAK5PFL1Offset'),
    #InputTagSmearedUp   = cms.InputTag('smearedPatJetsAK5PFresUp'),                                 
    #InputTagSmearedDown = cms.InputTag('smearedPatJetsAK5PFresDown'),                                 
    #InputTagScaledUp    = cms.InputTag('shiftedPatJetsAK5PFenUpForCorrMEt'),                                 
    #InputTagScaledDown  = cms.InputTag('shiftedPatJetsAK5PFenDownForCorrMEt'),                                 
    MVAPileupIDName = cms.string('pileupJetId:fullDiscriminant'),
    Prefix = cms.string('PFJet'),
    Suffix = cms.string('AK5'),
    MaxSize = cms.uint32(30),
    JECUncertainty = cms.string('AK5PF'),
    ReadJECuncertainty = cms.bool(True),
    VertexInputTag = cms.InputTag('offlineSlimmedPrimaryVertices')
)
rootTuplePFJetsAK5CHS = cms.EDProducer("RootTupleMakerV2_PFJets",
    InputTag = cms.InputTag('patJetsAK5PFCHS'),
    #InputTag = cms.InputTag('pfjetTriggerMatchEmbedderHLTEleJetJet'),
    ## InputTagL1Offset    = cms.InputTag('selectedPatJetsAK5PFL1Offset'),
    #InputTagSmearedUp   = cms.InputTag('smearedPatJetsAK5PFresUp'),                                 
    #InputTagSmearedDown = cms.InputTag('smearedPatJetsAK5PFresDown'),                                 
    #InputTagScaledUp    = cms.InputTag('shiftedPatJetsAK5PFenUpForCorrMEt'),                                 
    #InputTagScaledDown  = cms.InputTag('shiftedPatJetsAK5PFenDownForCorrMEt'),                                 
    MVAPileupIDName = cms.string('pileupJetId:fullDiscriminant'),
    Prefix = cms.string('PFJet'),
    Suffix = cms.string('AK5CHS'),
    MaxSize = cms.uint32(30),
    JECUncertainty = cms.string('AK5PF'),
    ReadJECuncertainty = cms.bool(True),
    VertexInputTag = cms.InputTag('offlineSlimmedPrimaryVertices')
)
rootTuplePFJetsAK4CHS = cms.EDProducer("RootTupleMakerV2_PFJets",
    InputTag = cms.InputTag('slimmedJets'),
    #InputTag = cms.InputTag('pfjetTriggerMatchEmbedderHLTEleJetJet'),
    ## InputTagL1Offset    = cms.InputTag('selectedPatJetsAK5PFL1Offset'),
    #InputTagSmearedUp   = cms.InputTag('smearedPatJetsAK5PFresUp'),                                 
    #InputTagSmearedDown = cms.InputTag('smearedPatJetsAK5PFresDown'),                                 
    #InputTagScaledUp    = cms.InputTag('shiftedPatJetsAK5PFenUpForCorrMEt'),                                 
    #InputTagScaledDown  = cms.InputTag('shiftedPatJetsAK5PFenDownForCorrMEt'),                                 
    MVAPileupIDName = cms.string('pileupJetId:fullDiscriminant'),
    Prefix = cms.string('PFJet'),
    Suffix = cms.string('AK4CHS'),
    MaxSize = cms.uint32(30),
    JECUncertainty = cms.string(''),
    ReadJECuncertainty = cms.bool(False),
    VertexInputTag = cms.InputTag('offlineSlimmedPrimaryVertices')
)
rootTuplePFJetsAK4Puppi = cms.EDProducer("RootTupleMakerV2_PFJets",
    InputTag = cms.InputTag('slimmedJetsPuppi'),
    #InputTag = cms.InputTag('pfjetTriggerMatchEmbedderHLTEleJetJet'),
    ## InputTagL1Offset    = cms.InputTag('selectedPatJetsAK5PFL1Offset'),
    #InputTagSmearedUp   = cms.InputTag('smearedPatJetsAK5PFresUp'),                                 
    #InputTagSmearedDown = cms.InputTag('smearedPatJetsAK5PFresDown'),                                 
    #InputTagScaledUp    = cms.InputTag('shiftedPatJetsAK5PFenUpForCorrMEt'),                                 
    #InputTagScaledDown  = cms.InputTag('shiftedPatJetsAK5PFenDownForCorrMEt'),                                 
    MVAPileupIDName = cms.string('pileupJetId:fullDiscriminant'),
    Prefix = cms.string('PFJet'),
    Suffix = cms.string('AK4Puppi'),
    MaxSize = cms.uint32(30),
    JECUncertainty = cms.string(''),
    ReadJECuncertainty = cms.bool(False),
    VertexInputTag = cms.InputTag('offlineSlimmedPrimaryVertices')
)

# Trigger matching
#FIXME? For now we just match to the path of interest
# Maybe later we can do more general matching
# This path seems strange in terms of saved jets, though.
# cross trigger strangeness?
# Maybe this HN post could help: https://hypernews.cern.ch/HyperNews/CMS/get/physTools/2372/1/1/2.html
pfjetTriggerMatchHLTEleJetJet = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'slimmedJets' )
, matched = cms.InputTag( 'unpackedPatTrigger' )          
, matchedCuts = cms.string( 'path ( "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*" )' )
#, matchedCuts = cms.string( 'type("TriggerPhoton") || type("TriggerElectron")' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

pfjetTriggerMatchEmbedderHLTEleJetJet = cms.EDProducer(
  "PATTriggerMatchElectronEmbedder"
, src = cms.InputTag( 'slimmedJets' )
, matches = cms.VInputTag(
  cms.InputTag( 'pfjetTriggerMatchHLTEleJetJet' )
  )
)
