import FWCore.ParameterSet.Config as cms

rootTupleTrigger = cms.EDProducer("RootTupleMakerV2_Trigger",
    L1InputTag  = cms.InputTag('gtDigis'),
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),

    #
    #inputs from                                 
    #1) 2E31 menu /cdaq/physics/Run2010/V7.0/HLT/V2
    #see https://twiki.cern.ch/twiki/pub/CMS/ExoticaLeptoquarkETC/Triggers.txt
    #2) 6E31 menu /cdaq/physics/Run2010/v8.0/HLT/V1
    #                              
    # see also:
    # http://indico.cern.ch/getFile.py/access?contribId=0&resId=1&materialId=slides&confId=99224
    # plans for 2E32:  https://twiki.cern.ch/twiki/bin/view/CMS/TMDLumi2E32v0
    #                  https://hypernews.cern.ch/HyperNews/CMS/get/exotica/868.html
    #3) 2E32 menu /cdaq/physics/Run2010/v9.0/HLT/V1
    #see https://twiki.cern.ch/twiki/bin/view/CMS/TMDLumi2E32v6
    #see also https://hypernews.cern.ch/HyperNews/CMS/get/commissioning/2315.html                              
                                  
    HLTPathsOfInterest = cms.vstring(
                                     ##############
                                     ## 2E31 + 6E31
                                     ##############
                                     #Electron
                                     'HLT_DoubleEle10_SW_L1R', #0
                                     'HLT_Ele12_SW_TightEleIdIsol_L1R', #1
                                     'HLT_Ele12_SW_TightEleIdIsol_NoDEtaInEE_L1R', #2
                                     'HLT_Ele17_SW_CaloEleId_L1R', #3
                                     'HLT_Ele17_SW_EleId_L1R', #4
                                     'HLT_Ele17_SW_LooseEleId_L1R', #5
                                     'HLT_Ele22_SW_CaloEleId_L1R', #6
                                     'HLT_Ele40_SW_L1R', #7
                                     'HLT_DoubleEle15_SW_L1R_v1', #8
                                     'HLT_Ele10_MET45_v1', #9
                                     'HLT_Ele12_SW_TighterEleIdIsol_L1R_v1', #10 
                                     'HLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1', #11
                                     'HLT_Ele17_SW_TightEleIdIsol_L1R_v1', #12
                                     'HLT_Ele17_SW_TightEleId_L1R', #13
                                     'HLT_Ele17_SW_TighterEleIdIsol_L1R_v1', #14
                                     'HLT_Ele17_SW_TighterEleId_L1R_v1', #15
                                     'HLT_Ele27_SW_TightCaloEleIdTrack_L1R_v1', #16
                                     'HLT_Ele32_SW_TightCaloEleIdTrack_L1R_v1', #17
                                     #Photon
                                     'HLT_DoublePhoton17_L1R', #18
                                     'HLT_DoublePhoton5_CEP_L1R', #19
                                     'HLT_Photon20_Cleaned_L1R', #20
                                     'HLT_Photon30_Cleaned_L1R', #21
                                     'HLT_Photon50_NoHE_Cleaned_L1R', #22
                                     'HLT_Photon100_NoHE_Cleaned_L1R_v1', #23
                                     'HLT_Photon17_SC17HE_L1R_v1', #24
                                     'HLT_Photon30_Isol_EBOnly_Cleaned_L1R_v1', #25
                                     'HLT_Photon35_Isol_Cleaned_L1R_v1', #26
                                     'HLT_Photon50_Cleaned_L1R_v1', #27
                                     'HLT_Photon70_NoHE_Cleaned_L1R_v1', #28
                                     #Jet
                                     'HLT_DiJetAve30U', #29
                                     'HLT_Jet30U', #30
                                     'HLT_Jet50U', #31
                                     'HLT_Jet70U', #32
                                     'HLT_Jet100U', #33
                                     'HLT_Jet70U_v2', #34
                                     'HLT_Jet100U_v2', #35
                                     'HLT_Jet140U_v1', #36
                                     #METFwd
                                     'HLT_MET45', #37
                                     'HLT_MET65', #38
                                     'HLT_MET100', #39
                                     'HLT_MET80_v1', #40 
                                     'HLT_MET100_v2', #41
                                     'HLT_MET45_HT100U_v1', #42
                                     'HLT_MET45_HT120U_v1', #43
                                     'HLT_HT100U', #44                                     
                                     'HLT_DoubleJet15U_ForwardBackward', #45
                                     'HLT_DoubleJet25U_ForwardBackward', #46
                                     #BTau
                                     'HLT_BTagMu_Jet20U', #47
                                     'HLT_DoubleIsoTau15_OneLeg_Trk5', #48
                                     'HLT_DoubleIsoTau15_Trk5', #49
                                     'HLT_SingleIsoTau20_Trk15_MET20', #50
                                     'HLT_SingleIsoTau20_Trk5_MET20', #51
                                     'HLT_SingleIsoTau30_Trk5_L120or30', #52
                                     'HLT_SingleIsoTau30_Trk5_MET20', #53
                                     'HLT_BTagMu_DiJet20U_v1', #54
                                     'HLT_SingleIsoTau30_Trk5_v2', #55
                                     #Mu
                                     'HLT_DoubleMu3', #56
                                     'HLT_Mu3', #57
                                     'HLT_Mu5', #58
                                     'HLT_Mu7', #59
                                     'HLT_Mu9', #60
                                     'HLT_Mu11', #61
                                     'HLT_IsoMu9', #62
                                     'HLT_Mu20_NoVertex', #63
                                     'HLT_DoubleMu5_v1', #64
                                     'HLT_Mu13_v1', #65
                                     'HLT_Mu15_v1', #66
                                     'HLT_Mu5_Ele9_v1', #67
                                     'HLT_Mu5_HT50U_v1', #68
                                     'HLT_Mu5_HT70U_v1', #69
                                     'HLT_Mu5_Jet35U_v1', #70
                                     'HLT_Mu5_Jet50U_v1', #71
                                     'HLT_Mu5_MET45_v1', #72
                                     'HLT_Mu5_Photon11_Cleaned_L1R_v1', #73
                                     #MinimumBias
                                     'HLT_ZeroBias', #74
                                     'HLT_L1Tech_BSC_minBias', #75
                                     'HLT_L1Tech_BSC_halo_forPhysicsBackground', #76
                                     'HLT_L1Tech_HCAL_HF', #77
                                     'HLT_L1_BPTX', #78
                                     'HLT_L1_BPTX_MinusOnly', #79
                                     'HLT_L1_BPTX_PlusOnly', #80
                                     'HLT_StoppedHSCP', #81
                                     'HLT_StoppedHSCP_v2', #82
                                     #Cosmics
                                     'HLT_L1Tech_BSC_halo', #83
                                     'HLT_TrackerCosmics', #84
                                     ##############
                                     ## 2E32
                                     ##############
                                     #Electron
                                     'HLT_DoubleEle17_SW_L1R_v1', #85
                                     'HLT_DoubleEle8_SW_HT70U_L1R_v1', #86
                                     'HLT_Ele10_SW_EleId_HT70U_L1R_v1', #87
                                     'HLT_Ele10_SW_HT100U_L1R_v1', #88
                                     'HLT_Ele10_SW_HT70U_L1R_v1', #89
                                     'HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1', #90
                                     'HLT_Ele17_SW_TighterEleIdIsol_L1R_v2', #91
                                     'HLT_Ele22_SW_TighterCaloIdIsol_L1R_v1', #92
                                     'HLT_Ele22_SW_TighterEleId_L1R_v2', #93
                                     'HLT_Ele32_SW_TighterEleId_L1R_v2', #94
                                     'HLT_IsoEle12_PFTau15_v1', #95
                                     #Photon
                                     'HLT_DoublePhoton17_SingleIsol_L1R_v1', #96
                                     'HLT_DoublePhoton22_L1R_v1', #97
                                     'HLT_DoublePhoton5_CEP_L1R_v3', #98
                                     'HLT_Photon110_NoHE_Cleaned_L1R_v1', #99
                                     'HLT_Photon17_Isol_SC17HE_L1R_v1', #100
                                     'HLT_Photon22_SC22HE_L1R_v1', #101
                                     'HLT_Photon40_CaloId_Cleaned_L1R_v1', #102
                                     'HLT_Photon40_Isol_Cleaned_L1R_v1', #103
                                     'HLT_Photon50_Cleaned_L1R_v1', #104
                                     'HLT_Photon70_Cleaned_L1R_v1', #105
                                     #Mu
                                     'HLT_IsoMu13_v3', #106
                                     'HLT_Mu17_v1', #107
                                     'HLT_Mu21_v1', #108
                                     'HLT_Mu8_Ele8_v1', #109
                                     'HLT_Mu5_HT100U_v3', #110
                                     'HLT_Mu5_HT70U_v3', #111
                                     #Jet
                                     'HLT_Jet50U_v3', #112
                                     'HLT_Jet140U_v3', #113
                                     #METFwd
                                     'HLT_MET100_v3', #114
                                     'HLT_MET120_v3', #115
                                     'HLT_MET45_DiJet30U_v3', #116
                                     ##############
                                     ## Others (from Sam, W' analysis)
                                     ##############
                                     'HLT_Ele10_LW_L1R', #117
                                     'HLT_Ele15_SW_L1R', #118
                                     'HLT_Ele15_SW_EleId_L1R', #119
                                     'HLT_Ele20_SW_L1R', #120
                                     'HLT_Ele15_SW_CaloEleId_L1R' #121                                     
                                     'HLT_Ele22_SW_TighterEleId_L1R_v3', #122                                     
                                     )

    #used in Spring10 - Summer10 analysis                                  
    #HLTPathsOfInterest = cms.vstring('HLT_Mu9','HLT_Photon15_L1R','HLT_Photon20_L1R','HLT_Photon30_L1R_8E29','HLT_Ele15_LW_L1R','HLT_Jet30U','HLT_MET45','HLT_MinBiasBSC','HLT_Photon15_Cleaned_L1R','HLT_Photon20_Cleaned_L1R')
)
