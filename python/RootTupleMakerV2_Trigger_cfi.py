import FWCore.ParameterSet.Config as cms

rootTupleTrigger = cms.EDProducer("RootTupleMakerV2_Trigger",
    L1uGTInputTag = cms.InputTag('gtStage2Digis'),
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),
    PackedPrescalesInputTag = cms.InputTag('patTrigger'),

    # HLT config browser : http://j2eeps.cern.ch/cms-project-confdb-hltdev/browser/
    # HLT Tools https://twiki.cern.ch/twiki/bin/view/CMS/HLTriggerTools
    # ID/Iso egamma definitions: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaWorkingPointsv3

    # keep if the full HLT path name contains any of these substrings
    HLTPathsOfInterest = cms.vstring(
        # Electron
        'HLT_Ele27',
        'HLT_Ele32',
        'HLT_Ele35',
        'HLT_Ele40',
        'HLT_Ele45',
        'HLT_Ele105',
        'HLT_Ele115',
        # Photon
        'HLT_Photon20',
        'HLT_Photon22',
        'HLT_Photon25',
        'HLT_Photon30',
        'HLT_Photon33',
        'HLT_Photon36',
        'HLT_Photon50',
        'HLT_Photon75',
        'HLT_Photon90',
        'HLT_Photon120',
        'HLT_Photon150',
        'HLT_Photon175',
        'HLT_Photon200',
        # Mu
        'HLT_Mu',
    )
                                     
)
