import FWCore.ParameterSet.Config as cms

rootTupleTrigger = cms.EDProducer("RootTupleMakerV2_Trigger",
    L1InputTag  = cms.InputTag('gtDigis'),
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),

    # HLT config browser : http://j2eeps.cern.ch/cms-project-confdb-hltdev/browser/

    HLTPathsOfInterest = cms.vstring(
                                     ##############
                                     ## 2011: run range 160329-163876
                                     ##############

                                     ## Single Electrons ## 
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

                                     ## Single Photons ##
                                     # from /cdaq/physics/Run2011/5e32/v4.2/HLT/V1
                                     'HLT_Photon30_CaloIdVL_IsoL_v1', #11
                                     'HLT_Photon30_CaloIdVL_v1', #12
                                     'HLT_Photon75_CaloIdVL_IsoL_v1', #13
                                     'HLT_Photon75_CaloIdVL_v1', #14
                                     # from /cdaq/physics/Run2011/5e32/v6.1/HLT/V1
                                     'HLT_Photon20_CaloIdVL_IsoL_v1', #15
                                     'HLT_Photon30_CaloIdVL_IsoL_v2', #16
                                     'HLT_Photon30_CaloIdVL_v2', #17
                                     'HLT_Photon50_CaloIdVL_IsoL_v1', #18
                                     'HLT_Photon75_CaloIdVL_IsoL_v2', #19
                                     'HLT_Photon75_CaloIdVL_v2', #20
                                     # from /cdaq/physics/Run2011/5e32/v8.1/HLT/V5
                                     'HLT_Photon20_CaloIdVL_IsoL_v2', #21
                                     'HLT_Photon30_CaloIdVL_IsoL_v3', #22
                                     'HLT_Photon30_CaloIdVL_v3', #23
                                     'HLT_Photon50_CaloIdVL_IsoL_v2', #24
                                     'HLT_Photon75_CaloIdVL_IsoL_v3', #25
                                     'HLT_Photon75_CaloIdVL_v3', #26
                                     
                                     ## Double Photons ##
                                     # from /cdaq/physics/Run2011/5e32/v4.2/HLT/V1
                                     'HLT_DoublePhoton33_v1', #
                                     'HLT_Photon26_CaloIdL_IsoVL_Photon18_v1', #
                                     'HLT_Photon26_IsoVL_Photon18_v1', #
                                     'HLT_Photon32_CaloIdL_Photon26_CaloIdL_v1', #
                                     # from /cdaq/physics/Run2011/5e32/v6.1/HLT/V1
                                     'HLT_DoublePhoton33_v2', #
                                     'HLT_Photon26_CaloIdL_IsoVL_Photon18_v2', #
                                     'HLT_Photon26_IsoVL_Photon18_v2', #
                                     'HLT_Photon32_CaloIdL_Photon26_CaloIdL_v2', #
                                     # from /cdaq/physics/Run2011/5e32/v8.1/HLT/V5
                                     'HLT_DoublePhoton33_v3', #
                                     'HLT_Photon26_CaloIdL_IsoVL_Photon18_v3', #
                                     'HLT_Photon26_IsoVL_Photon18_IsoVL_v3', #
                                     'HLT_Photon32_CaloIdL_Photon26_CaloIdL_v3', #

                                     ## Muons ##                                     

                                     )

)
