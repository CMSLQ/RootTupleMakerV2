import FWCore.ParameterSet.Config as cms

rootTupleMuons = cms.EDProducer("RootTupleMakerV2_Muons",
    InputTag = cms.InputTag('cleanPatMuons'),
    TriggerEventInputTag = cms.InputTag ('patTriggerEvent'),                                
    Prefix = cms.string('Muon'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10),
    MuonIso = cms.double(0.05),
    MuonID = cms.string('GlobalMuonPromptTight'),
    BeamSpotCorr = cms.bool(True),
    UseCocktailRefits = cms.bool(True),
    VertexInputTag = cms.InputTag('offlinePrimaryVertices'),
    SingleMuonTriggerMatch = cms.string ("cleanMuonTriggerMatchHLTSingleMuon"),
    SingleIsoMuonTriggerMatch = cms.string ("cleanMuonTriggerMatchHLTSingleIsoMuon")
)

cleanMuonTriggerMatchHLTSingleMuon = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'cleanPatMuons' )
, matched = cms.InputTag( 'patTrigger' )          
, matchedCuts = cms.string( 'type( "TriggerMuon" ) && path ( "HLT_Mu40_eta2p1_v*" )' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

cleanMuonTriggerMatchHLTSingleIsoMuon = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'cleanPatMuons' )
, matched = cms.InputTag( 'patTrigger' )          
, matchedCuts = cms.string( 'type( "TriggerMuon" ) && path ( "HLT_IsoMu24_eta2p1_v*" )' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)
