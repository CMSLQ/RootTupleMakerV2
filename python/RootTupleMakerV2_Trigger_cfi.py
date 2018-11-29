import FWCore.ParameterSet.Config as cms

rootTupleTrigger = cms.EDProducer("RootTupleMakerV2_Trigger",
    L1uGTInputTag = cms.InputTag('gtStage2Digis'),
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    PackedPrescalesInputTag = cms.InputTag('patTrigger'),

    # HLT config browser : http://j2eeps.cern.ch/cms-project-confdb-hltdev/browser/
    # HLT Tools https://twiki.cern.ch/twiki/bin/view/CMS/HLTriggerTools
    # ID/Iso egamma definitions: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaWorkingPointsv3

    HLTPathsOfInterest = cms.vstring(
        'HLT_Ele',
        'HLT_Mu',
        'HLT_Photon',
    )
                                     
)
