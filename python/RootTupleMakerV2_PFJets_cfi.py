import FWCore.ParameterSet.Config as cms

#rootTuplePFJetsAK5 = cms.EDProducer("RootTupleMakerV2_PFJets",
#    InputTag = cms.InputTag('selectedPatJetsAK5PF'),
#    #InputTag = cms.InputTag('pfjetTriggerMatchEmbedderHLTEleJetJet'),
#    ## InputTagL1Offset    = cms.InputTag('selectedPatJetsAK5PFL1Offset'),
#    #InputTagSmearedUp   = cms.InputTag('smearedPatJetsAK5PFresUp'),                                 
#    #InputTagSmearedDown = cms.InputTag('smearedPatJetsAK5PFresDown'),                                 
#    #InputTagScaledUp    = cms.InputTag('shiftedPatJetsAK5PFenUpForCorrMEt'),                                 
#    #InputTagScaledDown  = cms.InputTag('shiftedPatJetsAK5PFenDownForCorrMEt'),                                 
#    MVAPileupIDName = cms.string('pileupJetId:fullDiscriminant'),
#    Prefix = cms.string('PFJet'),
#    Suffix = cms.string('AK5'),
#    MaxSize = cms.uint32(30),
#    JECUncertainty = cms.string('AK5PF'),
#    ReadJECuncertainty = cms.bool(False),
#    ReadJERuncertainty = cms.bool(False),
#    JERUncertainty = cms.string('AK5PF'),
#    ReadJERFromGT = cms.bool(True),
#    RhoCollection = cms.InputTag('fixedGridRhoFastjetAll'),
#    JERResolutionsFile = cms.FileInPath(''),#FIXME why empty?
#    JERScaleFactorsFile = cms.FileInPath('Summer15_25nsV6_MC_JER_SF_AK4PFchs.txt'),
#    VertexInputTag = cms.InputTag('offlineSlimmedPrimaryVertices')
#)
#rootTuplePFJetsAK5CHS = cms.EDProducer("RootTupleMakerV2_PFJets",
#    InputTag = cms.InputTag('selectedPatJetsAK5PFCHS'),
#    #InputTag = cms.InputTag('pfjetTriggerMatchEmbedderHLTEleJetJet'),
#    ## InputTagL1Offset    = cms.InputTag('selectedPatJetsAK5PFL1Offset'),
#    #InputTagSmearedUp   = cms.InputTag('smearedPatJetsAK5PFresUp'),                                 
#    #InputTagSmearedDown = cms.InputTag('smearedPatJetsAK5PFresDown'),                                 
#    #InputTagScaledUp    = cms.InputTag('shiftedPatJetsAK5PFenUpForCorrMEt'),                                 
#    #InputTagScaledDown  = cms.InputTag('shiftedPatJetsAK5PFenDownForCorrMEt'),                                 
#    MVAPileupIDName = cms.string('pileupJetId:fullDiscriminant'),
#    Prefix = cms.string('PFJet'),
#    Suffix = cms.string('AK5CHS'),
#    MaxSize = cms.uint32(30),
#    JECUncertainty = cms.string('AK5PF'),
#    ReadJECuncertainty = cms.bool(True),
#    JERUncertainty = cms.string('AK5PF'),
#    ReadJERFromGT = cms.bool(True),
#    RhoCollection = cms.InputTag('fixedGridRhoFastjetAll'),
#    JERResolutionsFile = cms.FileInPath('Summer15_25nsV6_MC_PtResolution_AK4PFchs.txt'),
#    JERScaleFactorsFile = cms.FileInPath('Summer15_25nsV6_MC_JER_SF_AK4PFchs.txt'),
#    VertexInputTag = cms.InputTag('offlineSlimmedPrimaryVertices')
#)

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
    JECUncertainty = cms.string('AK4PFchs'),
    ReadJECuncertainty = cms.bool(True),
    ReadJERuncertainty = cms.bool(True),
    JERUncertainty = cms.string('AK4PFchs'),
    ReadJERFromGT = cms.bool(True),
    RhoCollection = cms.InputTag('fixedGridRhoFastjetAll'),
    JERResolutionsFile = cms.FileInPath('Summer15_25nsV6_MC_PtResolution_AK4PFchs.txt'),
    JERScaleFactorsFile = cms.FileInPath('Summer15_25nsV6_MC_JER_SF_AK4PFchs.txt'),
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
    JECUncertainty = cms.string('AK4PFPuppi'),
    ReadJECuncertainty = cms.bool(False),
    ReadJERuncertainty = cms.bool(False),
    JERUncertainty = cms.string('AK4PFPuppi'),
    ReadJERFromGT = cms.bool(True),
    RhoCollection = cms.InputTag('fixedGridRhoFastjetAll'),
    JERResolutionsFile = cms.FileInPath('Summer15_25nsV6_MC_PtResolution_AK4PFPuppi.txt'),
    JERScaleFactorsFile = cms.FileInPath('Summer15_25nsV6_MC_JER_SF_AK4PFPuppi.txt'),
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
, src     = cms.InputTag( 'slimmedElectrons' )
, matched = cms.InputTag( 'unpackedPatTrigger' )          
, matchedCuts = cms.string( 'path ( "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*" )' )
#, matchedCuts = cms.string( 'type("TriggerPhoton") || type("TriggerElectron")' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

pfjetTriggerMatchEmbedderHLTEleJetJet = cms.EDProducer(
  "PATTriggerMatchElectronEmbedder"
, src = cms.InputTag( 'slimmedElectrons' )
, matches = cms.VInputTag(
  cms.InputTag( 'pfjetTriggerMatchHLTEleJetJet' )
  )
)

rootTuplePFJetsSequence = cms.Sequence(
  pfjetTriggerMatchHLTEleJetJet*
  pfjetTriggerMatchEmbedderHLTEleJetJet*
  #rootTuplePFJetsAK5*
  #rootTuplePFJetsAK5CHS*
  rootTuplePFJetsAK4CHS*
  rootTuplePFJetsAK4Puppi
)

