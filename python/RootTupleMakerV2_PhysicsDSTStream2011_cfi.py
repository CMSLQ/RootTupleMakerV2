import FWCore.ParameterSet.Config as cms

rootTuplePhysicsDSTStream2011 = cms.EDProducer("RootTupleMakerV2_PhysicsDSTStream2011",
    Suffix = cms.string(''),
    StoreEventInfo = cms.bool(True),
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    InputTagHLTPFJets = cms.InputTag('hltAntiKT5PFJets','','HLT'),
    PrefixHLTPFJets = cms.string('HLTPFJet'),                                               
    MinPtHLTPFJets = cms.double(0),
    MaxEtaHLTPFJets = cms.double(3),
    InputTagHLTCaloJetsRaw = cms.InputTag('hltAntiKT5CaloJetsSelected','','HLT'),
    PrefixHLTCaloJetsRaw = cms.string('HLTCaloJetRaw'),                                               
    MinPtHLTCaloJetsRaw = cms.double(10),
    MaxEtaHLTCaloJetsRaw = cms.double(3),
    InputTagHLTCaloJetsCorr = cms.InputTag('hltCaloJetCorrectedSelected','','HLT'),
    PrefixHLTCaloJetsCorr = cms.string('HLTCaloJetCorr'),                                               
    MinPtHLTCaloJetsCorr = cms.double(10),
    MaxEtaHLTCaloJetsCorr = cms.double(3),
    InputTagHLTPixelVertices = cms.InputTag('hltPixelVertices','','HLT'),
    PrefixHLTPixelVertices = cms.string('HLTPixelVertex') 
)

