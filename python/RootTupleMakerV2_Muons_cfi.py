import FWCore.ParameterSet.Config as cms

rootTupleMuons = cms.EDProducer("RootTupleMakerV2_Muons",
    #InputTag = cms.InputTag('slimmedMuons'),
    InputTag = cms.InputTag('muonsTriggerMatchAll'),
    #TriggerEventInputTag = cms.InputTag ('patTriggerEvent'),                                
    Prefix = cms.string('Muon'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10),
    MuonIso = cms.double(0.05),
    MuonID = cms.string('GlobalMuonPromptTight'),
    BeamSpotCorr = cms.bool(True),
    UseCocktailRefits = cms.bool(True),
    VertexInputTag = cms.InputTag('offlineSlimmedPrimaryVertices'),
    #SingleMuonTriggerMatch = cms.string ("cleanMuonTriggerMatchHLTSingleMuon"),
    #SingleIsoMuonTriggerMatch = cms.string ("cleanMuonTriggerMatchHLTSingleIsoMuon")
)


cleanMuonTriggerMatchHLTMuon = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'slimmedMuons' )
, matched = cms.InputTag( 'unpackedPatTrigger' )          
, matchedCuts = cms.string( 'type( "TriggerMuon" )' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)


cleanMuonTriggerMatchHLTSingleMuon = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'slimmedMuons' )
, matched = cms.InputTag( 'patTrigger' )          
, matchedCuts = cms.string( 'type( "TriggerMuon" ) && path ( "HLT_Mu40_v*" )' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

cleanMuonTriggerMatchHLTSingleIsoMuon = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'slimmedMuons' )
, matched = cms.InputTag( 'patTrigger' )          
, matchedCuts = cms.string( 'type( "TriggerMuon" ) && path ( "HLT_IsoMu24_IterTrk02_v*" )' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)


# embed the trigger objects into the muons
muonsTriggerMatchHLTSingleMuon = cms.EDProducer(
  "PATTriggerMatchMuonEmbedder"
, src = cms.InputTag( 'slimmedMuons' )
, matches = cms.VInputTag(
  cms.InputTag( 'cleanMuonTriggerMatchHLTSingleMuon' )
  )
)

muonsTriggerMatchHLTSingleIsoMuon = cms.EDProducer(
  "PATTriggerMatchMuonEmbedder"
, src = cms.InputTag( 'slimmedMuons' )
, matches = cms.VInputTag(
  cms.InputTag( 'cleanMuonTriggerMatchHLTSingleIsoMuon' )
  )
)


muonsTriggerMatchAll = cms.EDProducer(
  "PATTriggerMatchMuonEmbedder"
, src = cms.InputTag( 'slimmedMuons' )
, matches = cms.VInputTag(
  cms.InputTag( 'cleanMuonTriggerMatchHLTMuon' )
  )
)



# these will skim out (from the collections above with embedded matches) only the objects which have a trigger object match
#from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import selectedPatMuons
#muonsTriggeredHLTSingleMuon = selectedPatMuons.clone(
#   src = cms.InputTag( 'muonsTriggerMatchHLTSingleMuon' )
#, cut = 'triggerObjectMatches.size > 0'
#)#
#
#muonsTriggeredHLTSingleIsoMuon = selectedPatMuons.clone(
#   src = cms.InputTag( 'muonsTriggerMatchHLTSingleIsoMuon' )
#, cut = 'triggerObjectMatches.size > 0'
#)
