import FWCore.ParameterSet.Config as cms

rootTupleElectrons = cms.EDProducer("RootTupleMakerV2_Electrons",
    #DCSInputTag = cms.InputTag('scalersRawToDigi'),
    InputTag = cms.InputTag('slimmedElectrons'),
    #InputTag = cms.InputTag('electronsTriggerMatchAll'), # use all TriggerElectron/TriggerPhoton matches by default
    Prefix = cms.string('Electron'),
    Suffix = cms.string(''),
    MaxSize = cms.uint32(10),
    ElectronIso = cms.double(0.1),
    MuonPt = cms.double(10.),
    MuonIso = cms.double(0.05),
    MuonID = cms.string('GlobalMuonPromptTight'),
    VertexInputTag = cms.InputTag('offlineSlimmedPrimaryVertices'),
    BeamSpotInputTag = cms.InputTag('offlineBeamSpot'),
    LikelihoodInputTag = cms.InputTag('egammaIDLikelihood') ,
    # rho for HEEP. See e.g. https://github.com/cms-sw/cmssw/blob/CMSSW_7_5_X/RecoEgamma/ElectronIdentification/python/Identification/heepElectronID_HEEPV50_CSA14_startup_cff.py
    RhoInputTag = cms.InputTag('fixedGridRhoFastjetAll'),
    # SIC add
    ElectronVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
    ElectronLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
    ElectronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
    ElectronTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
    ElectronHLTPreselectionMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1"),
    ElectronHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"),
    ElectronMVAIdWP80Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
    ElectronMVAIdWP90Map = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
    ElectronMVAIdHZZMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-HZZ-V1-wpLoose"),
    # The map name for the full info is the same as the map name of the
    # corresponding simple pass/fail map above, they are distinguished by
    # the type of the content.
    eleVetoIdCutFlowResultMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
    eleLooseIdCutFlowResultMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
    eleMediumIdCutFlowResultMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
    eleTightIdCutFlowResultMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
    eleHLTPreselectionCutFlowResultMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1"),
    eleHEEPIdCutFlowResultMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"),
    ElectronMVAIdWP80CutFlowResultMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
    ElectronMVAIdWP90CutFlowResultMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
    ElectronMVAIdHZZCutFlowResultMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-HZZ-V1-wpLoose"),
    # new HEEP 7.0 track isolation
    heep70trkIsolMap = cms.InputTag("heepIDVarValueMaps","eleTrkPtIso"),
    EBReducedRecHitsInputTag = cms.InputTag("reducedEgamma:reducedEBRecHits"),
    EEReducedRecHitsInputTag = cms.InputTag("reducedEgamma:reducedEERecHits"),
)

# Trigger matching
# Note: typically, electrons don't come as TriggerElectron but as TriggerPhoton
# See: https://hypernews.cern.ch/HyperNews/CMS/get/egamma-hlt/201/1.html
cleanElectronTriggerMatchHLTSingleElectron = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
  #"PATTriggerMatcherDRDPtLessByR"
, src     = cms.InputTag( 'slimmedElectrons' )
, matched = cms.InputTag( 'unpackedPatTrigger' )          
#, matchedCuts = cms.string( 'path ( "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*" )' )
, matchedCuts = cms.string( 'type("TriggerPhoton") || type("TriggerElectron")' )
#, maxDeltaR = cms.double( 0.5 )
, maxDeltaR = cms.double( 1.0 )
#, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
#, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
, resolveAmbiguities    = cms.bool( False  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( False  )        # take best match found per reco object: by DeltaR here (s. above)
)

cleanElectronTriggerMatchHLTSingleElectronWP85 = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'slimmedElectrons' )
, matched = cms.InputTag( 'unpackedPatTrigger' )          
, matchedCuts = cms.string( 'type( "TriggerPhoton" ) && path( "HLT_Ele32_eta2p1_WP85_Gsf_v*" )' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

cleanElectronTriggerMatchHLTDoubleElectron = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'slimmedElectrons' )
, matched = cms.InputTag( 'unpackedPatTrigger' ) 
, matchedCuts = cms.string( 'type( "TriggerCluster" ) && path( "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*" )' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

cleanElectronTriggerMatchHLTElectronJetJet = cms.EDProducer(
  "PATTriggerMatcherDRLessByR"
, src     = cms.InputTag( 'slimmedElectrons' )
, matched = cms.InputTag( 'unpackedPatTrigger' ) 
, matchedCuts = cms.string( 'path( "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*" )' )
, maxDeltaR = cms.double( 0.5 )
, resolveAmbiguities    = cms.bool( True  )        # only one match per trigger object
, resolveByMatchQuality = cms.bool( True  )        # take best match found per reco object: by DeltaR here (s. above)
)

# embed the trigger objects into the electrons
#electronsTriggerMatchHLTSingleElectron = cms.EDProducer(
#  "PATTriggerMatchElectronEmbedder"
#, src = cms.InputTag( 'slimmedElectrons' )
#, matches = cms.VInputTag(
#  cms.InputTag( 'cleanElectronTriggerMatchHLTSingleElectron' )
#  )
#)

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

electronsTriggerMatchAll = cms.EDProducer(
  "PATTriggerMatchElectronEmbedder"
, src = cms.InputTag( 'slimmedElectrons' )
, matches = cms.VInputTag(
  cms.InputTag( 'cleanElectronTriggerMatchHLTSingleElectron' )
  )
)

# add all to the sequence
rootTupleElectronsSequence = cms.Sequence(
  cleanElectronTriggerMatchHLTSingleElectron*
  cleanElectronTriggerMatchHLTSingleElectronWP85*
  cleanElectronTriggerMatchHLTDoubleElectron*
  cleanElectronTriggerMatchHLTElectronJetJet*
  electronsTriggerMatchHLTSingleElectronWP85*
  electronsTriggerMatchHLTDoubleElectron*
  electronsTriggerMatchHLTElectronJetJet*
  electronsTriggerMatchAll*
  rootTupleElectrons
)

# these will skim out (from the collections above with embedded matches) only the objects which have a trigger object match
#from PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi import selectedPatElectrons
#electronsTriggeredHLTSingleElectron = selectedPatElectrons.clone(
#   src = cms.InputTag( 'electronsTriggerMatchHLTSingleElectron' )
#, cut = 'triggerObjectMatches.size > 0'
#)
#
#electronsTriggeredHLTSingleElectronWP80 = selectedPatElectrons.clone(
#   src = cms.InputTag( 'electronsTriggerMatchHLTSingleElectronWP80' )
#, cut = 'triggerObjectMatches.size > 0'
#)
#
#electronsTriggeredHLTDoubleElectron = selectedPatElectrons.clone(
#   src = cms.InputTag( 'electronsTriggerMatchHLTDoubleElectron' )
#, cut = 'triggerObjectMatches.size > 0'
#)

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
