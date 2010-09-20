import FWCore.ParameterSet.Config as cms

rootTupleTrigger = cms.EDProducer("RootTupleMakerV2_Trigger",
    L1InputTag  = cms.InputTag('gtDigis'),
    HLTInputTag = cms.InputTag('TriggerResults','','HLT'),

    #from 2E31 menu /online/collisions/2010/week37/HLT/V8 (+ predictions for 6E31)
    #see https://twiki.cern.ch/twiki/pub/CMS/ExoticaLeptoquarkETC/Triggers.txt                             
    HLTPathsOfInterest = cms.vstring('HLT_DoubleEle10_SW_L1R', #0
                                     'HLT_Ele12_SW_TightEleIdIsol_NoDEtaInEE_L1R',#1
                                     'HLT_Ele12_SW_TightEleIdIsol_L1R',#2
                                     'HLT_Ele17_SW_CaloEleId_L1R',#3
                                     'HLT_Ele17_SW_LooseEleId_L1R',#4
                                     'HLT_Ele22_SW_CaloEleId_L1R',#5
                                     'HLT_Ele17_SW_EleId_L1R',#6
                                     'HLT_Ele40_SW_L1R',#7
                                     'HLT_Photon20_Cleaned_L1R',#8
                                     'HLT_DoublePhoton5_CEP_L1R',#9
                                     'HLT_Photon30_Cleaned_L1R',#10
                                     'HLT_DoublePhoton17_L1R',#11
                                     'HLT_Photon50_NoHE_Cleaned_L1R',#12
                                     'HLT_Jet30U',#13
                                     'HLT_Jet50U',#14
                                     'HLT_Jet70U',#15
                                     'HLT_Jet100U',#16
                                     'HLT_DiJetAve15U',#17
                                     'HLT_HT100U',#18
                                     'HLT_MET45',#19
                                     'HLT_MET65',#20
                                     'HLT_MET80',#21 (not available in 2E31 menu, might be in 6E31 menu)
                                     'HLT_MET100',#22
                                     'HLT_DoubleJet15U_ForwardBackward',#23
                                     'HLT_DoubleJet25U_ForwardBackward',#24
                                     'HLT_BTagMu_Jet20U',#25
                                     'HLT_SingleIsoTau20_Trk5_MET20',#26
                                     'HLT_SingleIsoTau20_Trk15_MET20',#27
                                     'HLT_DoubleIsoTau15_OneLeg_Trk5',#28
                                     'HLT_DoubleIsoTau15_Trk5',#29
                                     'HLT_SingleIsoTau30_Trk5_MET20',#30
                                     'HLT_SingleIsoTau30_Trk5_L120or30',#31
                                     'HLT_DoubleMu3',#32
                                     'HLT_Mu3',#33
                                     'HLT_Mu5',#34
                                     'HLT_Mu7',#35
                                     'HLT_Mu9',#36
                                     'HLT_Mu11',#37
                                     'HLT_Mu13',#38 (not available in 2E31 menu, might be in 6E31 menu)
                                     'HLT_Mu15',#39 (not available in 2E31 menu, might be in 6E31 menu)
                                     'HLT_IsoMu9',#40
                                     'HLT_Mu20_NoVertex',#41
                                     'HLT_ZeroBias',#42
                                     'HLT_L1Tech_BSC_minBias',#43
                                     'HLT_L1Tech_BSC_halo_forPhysicsBackground',#44
                                     'HLT_L1Tech_HCAL_HF',#45
                                     'HLT_L1_BPTX',#46
                                     'HLT_L1_BPTX_MinusOnly',#47
                                     'HLT_L1_BPTX_PlusOnly',#48
                                     'HLT_StoppedHSCP',#49
                                     'HLT_L1Tech_BSC_halo',#50
                                     'HLT_TrackerCosmics'#51
                                     )

    #used in Spring10 - Summer10 analysis                                  
    #HLTPathsOfInterest = cms.vstring('HLT_Mu9','HLT_Photon15_L1R','HLT_Photon20_L1R','HLT_Photon30_L1R_8E29','HLT_Ele15_LW_L1R','HLT_Jet30U','HLT_MET45','HLT_MinBiasBSC','HLT_Photon15_Cleaned_L1R','HLT_Photon20_Cleaned_L1R')
)
