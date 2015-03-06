import FWCore.ParameterSet.Config as cms

rootTupleElectrons = cms.EDProducer("RootTupleMakerV2_Electrons",
    DCSInputTag = cms.InputTag('scalersRawToDigi'),
    InputTag = cms.InputTag('slimmedElectrons'),
    Prefix = cms.string('Electron'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10),
    ElectronIso = cms.double(0.1),
    MuonPt = cms.double(10.),
    MuonIso = cms.double(0.05),
    MuonID = cms.string('GlobalMuonPromptTight'),
    VertexInputTag = cms.InputTag('offlineSlimmedPrimaryVertices'),
    BeamSpotInputTag = cms.InputTag('offlineBeamSpot'),
    ConversionsInputTag = cms.InputTag('allConversions'),
    TriggerEventInputTag = cms.InputTag('patTriggerEvent'),                                    
    LikelihoodInputTag = cms.InputTag('egammaIDLikelihood') ,
    # rho for HEEP. See e.g. https://github.com/cms-sw/cmssw/blob/CMSSW_7_5_X/RecoEgamma/ElectronIdentification/python/Identification/heepElectronID_HEEPV50_CSA14_startup_cff.py
    RhoInputTag = cms.InputTag('fixedGridRhoFastjetAll'),
    # SIC add
    ElectronVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V1-miniAOD-standalone-veto"),
    ElectronTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V1-miniAOD-standalone-tight"),
    ElectronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V1-miniAOD-standalone-medium"),
    ElectronLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V1-miniAOD-standalone-loose"),
    ElectronHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51-miniAOD"),
)

cleanElectronTriggerMatchHLTSingleElectron = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'slimmedElectrons' )
, matched = cms.InputTag( 'patTrigger' )          
, matchedCuts = cms.string( 'type( "TriggerElectron" ) && path ( "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*" )' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

cleanElectronTriggerMatchHLTSingleElectronWP85 = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'slimmedElectrons' )
, matched = cms.InputTag( 'patTrigger' )          
, matchedCuts = cms.string( 'type( "TriggerElectron" ) && path( "HLT_Ele32_eta2p1_WP85_Gsf_v*" )' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

cleanElectronTriggerMatchHLTDoubleElectron = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'slimmedElectrons' )
, matched = cms.InputTag( 'patTrigger' ) 
, matchedCuts = cms.string( 'type( "TriggerCluster" ) && path( "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*" )' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

cleanElectronTriggerMatchHLTElectronJetJet = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'slimmedElectrons' )
, matched = cms.InputTag( 'patTrigger' ) 
, matchedCuts = cms.string( 'type( "TriggerCluster" ) && path( "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*" )' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

# embed the trigger objects into the electrons
electronsTriggerMatchHLTSingleElectron = cms.EDProducer(
  "PATTriggerMatchElectronEmbedder"
, src = cms.InputTag( 'slimmedElectrons' )
, matches = cms.VInputTag(
  cms.InputTag( 'cleanElectronTriggerMatchHLTSingleElectron' )
  )
)

electronsTriggerMatchHLTSingleElectronWP85 = cms.EDProducer(
  "PATTriggerMatchElectronEmbedder"
, src = cms.InputTag( 'slimmedElectrons' )
, matches = cms.VInputTag(
  cms.InputTag( 'cleanElectronTriggerMatchHLTSingleElectronWP85' )
  )
)

electronsTriggerMatchHLTDoubleElectron = cms.EDProducer(
  "PATTriggerMatchElectronEmbedder"
, src = cms.InputTag( 'slimmedElectrons' )
, matches = cms.VInputTag(
  cms.InputTag( 'cleanElectronTriggerMatchHLTDoubleElectron' )
  )
)

electronsTriggerMatchHLTElectronJetJet = cms.EDProducer(
  "PATTriggerMatchElectronEmbedder"
, src = cms.InputTag( 'slimmedElectrons' )
, matches = cms.VInputTag(
  cms.InputTag( 'cleanElectronTriggerMatchHLTElectronJetJet' )
  )
)

# these will skim out (from the collections above with embedded matches) only the objects which have a trigger object match
from PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi import selectedPatElectrons
electronsTriggeredHLTSingleElectron = selectedPatElectrons.clone(
   src = cms.InputTag( 'electronsTriggerMatchHLTSingleElectron' )
, cut = 'triggerObjectMatches.size > 0'
)

electronsTriggeredHLTSingleElectronWP80 = selectedPatElectrons.clone(
   src = cms.InputTag( 'electronsTriggerMatchHLTSingleElectronWP80' )
, cut = 'triggerObjectMatches.size > 0'
)

electronsTriggeredHLTDoubleElectron = selectedPatElectrons.clone(
   src = cms.InputTag( 'electronsTriggerMatchHLTDoubleElectron' )
, cut = 'triggerObjectMatches.size > 0'
)

# Extra trigger matching (for QCD estimate).  Leave commented for now.
# 
# cleanElectronTriggerMatchHLTPhotonCaloIdVL = cms.EDProducer(
#   "PATTriggerMatcherDRLessByR"
# , src     = cms.InputTag( 'cleanPatElectrons' )
# , matched = cms.InputTag( 'patTrigger' )          
# , matchedCuts = cms.string( 'type( "TriggerPhoton" ) && path( "HLT_Photon*_CaloIdVL_v*" )' )
# , maxDeltaR = cms.double( 0.5 )
# , resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
# , resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
# )
# 
# cleanElectronTriggerMatchHLTPhoton135 = cms.EDProducer(
#   "PATTriggerMatcherDRLessByR"
# , src     = cms.InputTag( 'cleanPatElectrons' )
# , matched = cms.InputTag( 'patTrigger' )          
# , matchedCuts = cms.string( 'type( "TriggerPhoton" ) && path( "HLT_Photon135_v*" )' )
# , maxDeltaR = cms.double( 0.5 )
# , resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
# , resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
# )
# 
# cleanElectronTriggerMatchHLTPhoton150 = cms.EDProducer(
#   "PATTriggerMatcherDRLessByR"
# , src     = cms.InputTag( 'cleanPatElectrons' )
# , matched = cms.InputTag( 'patTrigger' )          
# , matchedCuts = cms.string( 'type( "TriggerPhoton" ) && path( "HLT_Photon150_v*" )' )
# , maxDeltaR = cms.double( 0.5 )
# , resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
# , resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
# )
