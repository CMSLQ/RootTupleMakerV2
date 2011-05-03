import FWCore.ParameterSet.Config as cms

rootTupleTrigger = cms.EDProducer("RootTupleMakerV2_Trigger",
    L1InputTag  = cms.InputTag('gtDigis'),
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),

    # HLT config browser : http://j2eeps.cern.ch/cms-project-confdb-hltdev/browser/

    HLTPathsOfInterest = cms.vstring(
                                     ##############
                                     ## 2011: run range 160329-163876
                                     ##############

                                     ## Electrons ## 
                                     # from /cdaq/physics/Run2011/5e32/v4.2/HLT/V1
                                     'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1', #0
                                     'HLT_Ele45_CaloIdVT_TrkIdT_v1', #1
                                     'HLT_Ele90_NoSpikeFilter_v1', #2
                                     # from /cdaq/physics/Run2011/5e32/v6.1/HLT/V1
                                     'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2', #3
                                     'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1', #4
                                     'HLT_Ele45_CaloIdVT_TrkIdT_v2', #5
                                     'HLT_Ele90_NoSpikeFilter_v2', #6
                                     # from /cdaq/physics/Run2011/5e32/v8.1/HLT/V5
                                     'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3', #7
                                     'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2', #8
                                     'HLT_Ele45_CaloIdVT_TrkIdT_v3', #9
                                     'HLT_Ele90_NoSpikeFilter_v3', #10
                                     ## Photons ##
                                     
                                     ## Muons ##                                     

                                     )

)
