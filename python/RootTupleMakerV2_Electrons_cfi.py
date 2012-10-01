import FWCore.ParameterSet.Config as cms

rootTupleElectrons = cms.EDProducer("RootTupleMakerV2_Electrons",
    TracksInputTag = cms.InputTag('generalTracks'),
    DCSInputTag = cms.InputTag('scalersRawToDigi'),
    InputTag = cms.InputTag('cleanPatElectrons'),
    Prefix = cms.string('Electron'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10),
    ElectronIso = cms.double(0.1),
    MuonPt = cms.double(10.),
    MuonIso = cms.double(0.05),
    MuonID = cms.string('GlobalMuonPromptTight'),
    VertexInputTag = cms.InputTag('offlinePrimaryVertices'),
    BeamSpotInputTag = cms.InputTag('offlineBeamSpot'),
    ConversionsInputTag = cms.InputTag('allConversions'),
    TriggerEventInputTag = cms.InputTag('patTriggerEvent'),                                    
    LikelihoodInputTag = cms.InputTag('egammaIDLikelihood') ,
    RhoInputTag = cms.InputTag('kt6PFJets','rho'),
    PFIsolationValues03 = cms.VInputTag ( cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                                          cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                                          cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),
                                    # In principle, these should be added too, but I get errors from PAT when I try to load them
                                    # They are not needed for EGamma PF isolation
                                    # An exception of category 'Configuration' occurred while
                                    # [0] Constructing the EventProcessor
                                    # [1] Constructing module: class=RootTupleMakerV2_Electrons label='rootTupleElectrons'
                                    # Exception Message:
                                    #     Duplicate Process The process name PAT was previously used on these products.
                                    # Please modify the configuration file to use a distinct process name.
                                    # cms.InputTag("elPFIsoValueChargedAll03PFIdPFIso"),  
                                    # cms.InputTag("elPFIsoValuePU03PFIdPFIso")),
    PFIsolationValues04 = cms.VInputTag ( cms.InputTag('elPFIsoValueCharged04PFIdPFIso'),
                                          cms.InputTag('elPFIsoValueGamma04PFIdPFIso'),
                                          cms.InputTag('elPFIsoValueNeutral04PFIdPFIso')),
                                    # In principle, these should be added too, but I get errors from PAT when I try to load them
                                    # They are not needed for EGamma PF isolation ( see above ) 
                                    # cms.InputTag("elPFIsoValueChargedAll04PFIdPFIso"), 
                                    # cms.InputTag("elPFIsoValuePU04PFIdPFIso")),
    SingleEleTriggerMatch     = cms.string ("cleanElectronTriggerMatchHLTSingleElectron"),
    SingleEleTriggerMatchWP80 = cms.string ("cleanElectronTriggerMatchHLTSingleElectronWP80"),
    DoubleEleTriggerMatch     = cms.string ("cleanElectronTriggerMatchHLTDoubleElectron")
)

cleanElectronTriggerMatchHLTSingleElectron = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'cleanPatElectrons' )
, matched = cms.InputTag( 'patTrigger' )          
, matchedCuts = cms.string( 'type( "TriggerElectron" ) && path ( "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*" )' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

cleanElectronTriggerMatchHLTSingleElectronWP80 = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'cleanPatElectrons' )
, matched = cms.InputTag( 'patTrigger' )          
, matchedCuts = cms.string( 'type( "TriggerElectron" ) && path( "HLT_Ele27_WP80_v*" )' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

cleanElectronTriggerMatchHLTDoubleElectron = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'cleanPatElectrons' )
, matched = cms.InputTag( 'patTrigger' )          
, matchedCuts = cms.string( 'type( "TriggerCluster" ) && path( "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*" )' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
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
