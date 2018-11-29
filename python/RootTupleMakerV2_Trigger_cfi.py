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
        'HLT_Mu'

                                     #---------------------------------------------
                                     ## SingleElectron ##
                                     #---------------------------------------------
                                     #'HLT_Ele27_WPTight_Gsf_v*',
                                     #'HLT_Ele115_CaloIdVT_GsfTrkIdT_v*',

                                     ##---------------------------------------------
                                     ### Photon ##
                                     ##---------------------------------------------
                                     ## from /cdaq/physics/Run2011/5e32/v4.2/HLT/V1
                                     ## 160404 <= run < 161216
                                     #'HLT_Photon30_CaloIdVL_IsoL_v1', #21
                                     #'HLT_Photon30_CaloIdVL_v1', #22
                                     #'HLT_Photon75_CaloIdVL_IsoL_v1', #23
                                     #'HLT_Photon75_CaloIdVL_v1', #24
                                     ## from /cdaq/physics/Run2011/5e32/v6.1/HLT/V1
                                     ## 161216 <= run < 163269
                                     #'HLT_Photon20_CaloIdVL_IsoL_v1', #25
                                     #'HLT_Photon30_CaloIdVL_IsoL_v2', #26
                                     #'HLT_Photon30_CaloIdVL_v2', #27
                                     #'HLT_Photon50_CaloIdVL_IsoL_v1', #28
                                     #'HLT_Photon75_CaloIdVL_IsoL_v2', #29
                                     #'HLT_Photon75_CaloIdVL_v2', #30
                                     ## from /cdaq/physics/Run2011/5e32/v8.1/HLT/V5
                                     ## 163269 <= run < 163876
                                     #'HLT_Photon20_CaloIdVL_IsoL_v2', #31
                                     #'HLT_Photon30_CaloIdVL_IsoL_v3', #32
                                     #'HLT_Photon30_CaloIdVL_v3', #33
                                     #'HLT_Photon50_CaloIdVL_IsoL_v2', #34
                                     #'HLT_Photon75_CaloIdVL_IsoL_v3', #35
                                     #'HLT_Photon75_CaloIdVL_v3', #36
                                     ##
                                     ## ==> nothing between 163876 and 165088 <==
                                     ##
                                     ## from /cdaq/physics/Run2011/1e33/v1.3/HLT/V2
                                     ## 165088 <= run < 165970
                                     #'HLT_Photon20_CaloIdVL_IsoL_v3', #37
                                     #'HLT_Photon30_CaloIdVL_IsoL_v4', #38
                                     #'HLT_Photon30_CaloIdVL_v4', #39
                                     #'HLT_Photon50_CaloIdVL_IsoL_v3', #40
                                     #'HLT_Photon50_CaloIdVL_v1', #41
                                     #'HLT_Photon75_CaloIdVL_IsoL_v4', #42
                                     #'HLT_Photon75_CaloIdVL_v4', #43
                                     #'HLT_Photon90_CaloIdVL_IsoL_v1', #44
                                     #'HLT_Photon90_CaloIdVL_v1', #45
                                     ## from /cdaq/physics/Run2011/1e33/v2.3/HLT/V1
                                     ## 165970 <= run < 165980
                                     #'HLT_Photon20_CaloIdVL_IsoL_v4', #46
                                     #'HLT_Photon30_CaloIdVL_IsoL_v5', #47
                                     #'HLT_Photon30_CaloIdVL_v5', #48
                                     #'HLT_Photon50_CaloIdVL_IsoL_v4', #49
                                     #'HLT_Photon50_CaloIdVL_v2', #50
                                     #'HLT_Photon75_CaloIdVL_IsoL_v5', #51
                                     #'HLT_Photon75_CaloIdVL_v5', #52
                                     #'HLT_Photon90_CaloIdVL_IsoL_v2', #53
                                     #'HLT_Photon90_CaloIdVL_v2', #54

                                     ##---------------------------------------------
                                     ### SingleMu ##
                                     ##---------------------------------------------
                                     ## from /cdaq/physics/Run2011/5e32/v4.2/HLT/V1
                                     ## 160404 <= run < 163269
                                     #'HLT_Mu20_v1', #55
                                     #'HLT_Mu24_v1', #56
                                     #'HLT_Mu30_v1', #57
                                     ## from /cdaq/physics/Run2011/5e32/v8.1/HLT/V5
                                     ## 163269 <= run < 163876
                                     #'HLT_Mu24_v2', #58
                                     #'HLT_Mu30_v2', #59
                                     ##
                                     ## ==> nothing between 163876 and 165088 <==
                                     ##
                                     ## from /cdaq/physics/Run2011/1e33/v1.3/HLT/V2
                                     ## 165088 <= run < 165980
                                     #'HLT_Mu30_v3', #60 
                                     #'HLT_Mu40_v1', #61
    )
                                     
)
