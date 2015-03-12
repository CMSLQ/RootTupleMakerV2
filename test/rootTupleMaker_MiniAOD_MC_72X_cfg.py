import os
#----------------------------------------------------------------------------------------------------
# Load PAT template + customize
#----------------------------------------------------------------------------------------------------
# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

#process.load('PhysicsTools.PatAlgos.patSequences_cff')
#process.load('Configuration.StandardSequences.Services_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff')

# Change process name
process._Process__name="ROOTTUPLEMAKERV2"
# Options and Output Report
process.options.wantSummary = False
process.options.allowUnscheduled = cms.untracked.bool(True)

############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.default.limit = 10
#################################################################

# We should be using PFIso by default in newer CMSSW
# see: PhysicsTools/PatAlgos/python/producersLayer1/electronProducer_cff.py

#----------------------------------------------------------------------------------------------------
# Load our RootTupleMakerV2 modules
#----------------------------------------------------------------------------------------------------

process.load('Leptoquarks.RootTupleMakerV2.Ntuple_cff')

# Output ROOT file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string( "file.root" )
)

#----------------------------------------------------------------------------------------------------
# Set global settings (number of events, global tag, input files, etc)
#----------------------------------------------------------------------------------------------------
# Make sure a correct global tag is used:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Valid_Global_Tags_by_Release
# XXX SIC: below possibly needed depending on CMSSW version
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PHYS14_25_V1', '')
# override the GlobalTag, connection string and pfnPrefix
if 'GlobalTag' in process.__dict__:
    from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag as customiseGlobalTag
    process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'MCRUN2_72_V1A::All')
    #process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_CONDITIONS'
    process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
    process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
    for pset in process.GlobalTag.toGet.value():
        pset.connect = pset.connect.value().replace('frontier://FrontierProd/', 'frontier://FrontierProd/')
    # fix for multi-run processing
    process.GlobalTag.RefreshEachRun = cms.untracked.bool( False )
    process.GlobalTag.ReconnectEachRun = cms.untracked.bool( False )

# Otherwise, just plain GlobalTag
#process.GlobalTag.globaltag = 'MCRUN2_72_V3'
#process.GlobalTag.globaltag = 'PHYS14_25_V1'

# Events to process
process.maxEvents.input = 1

# Input files
process.source.fileNames = [
    # specified by InputList.txt
    # Here is a test MiniAOD RelVal in CMSSW_7_2_2_patch1
    #'/store/relval/CMSSW_7_2_2_patch1/RelValZEE_13/MINIAODSIM/PU25ns_MCRUN2_72_V3_71XGENSIM-v2/00000/4E19A622-3474-E411-A100-0026189438E1.root'
    # Here is a file from dataset=/DYJetsToLL_M-50_13TeV-madgraph-pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM
    #'/store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0432E62A-7A6C-E411-87BB-002590DB92A8.root'
    # setup: https://cms-pdmv.cern.ch/mcm/public/restapi/requests/get_setup/TSG-Phys14DR-00003
    ## Here is a file from my Phys14-like MiniAOD LQ1 production
    '/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_10_1_su7.root',
    '/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_11_1_xCf.root',
    '/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_12_1_D39.root',
    '/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_13_1_wyC.root',
    '/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_14_1_ZeD.root',
    '/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_15_1_k6e.root',
    '/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_16_1_EGW.root',
    '/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_17_1_hdA.root',
    '/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_18_1_JrI.root',
    '/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_19_1_1N4.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_1_1_Snm.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_20_1_Myo.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_21_1_n4c.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_22_1_NIS.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_23_1_aNs.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_24_1_no8.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_25_1_EFd.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_26_1_zHI.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_27_1_Mtd.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_28_1_COV.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_29_1_TLD.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_2_1_Xnu.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_30_1_EFt.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_31_1_oKN.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_32_1_gim.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_33_1_MMj.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_34_1_Dcx.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_35_1_X4p.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_36_1_D7h.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_37_1_xkq.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_38_1_c6Q.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_39_1_dcb.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_3_1_gCV.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_40_1_5wJ.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_41_1_I8M.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_42_1_DVK.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_43_1_152.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_44_1_sT1.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_45_1_uEL.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_46_1_m0b.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_47_1_6EJ.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_48_1_GJi.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_49_1_X78.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_4_1_FDs.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_50_1_gts.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_5_1_i9b.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_6_1_HFo.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_7_1_aUd.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_8_1_6kI.root',
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/721/miniAOD-noPU/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8/LQToUE_M-300_Tune4C_13TeV-pythia8-721rawToDigiL1RecoRecoEI-noPU/66cb3c103bb16b7a3fee248046c0ac3b/step3_PAT_9_1_cc9.root',
]

#----------------------------------------------------------------------------------------------------
# HEEP 4.0 (electron ID) still uses the 2011 definitions of rho for isolation corrections.
# 
# Recipe taken from here:
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaEARhoCorrection#Rho_for_2011_Effective_Areas
#----------------------------------------------------------------------------------------------------
# SIC Replace with HEEP 5.1
# Also load Egamma cut-based ID in new VID framework while we're at it
# See: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
#   and for HEEP: https://hypernews.cern.ch/HyperNews/CMS/get/egamma/1519/2/1/1/1.html
# Set up everything that is needed to compute electron IDs and
# add the ValueMaps with ID decisions into the event data stream
#
# Load tools and function definitions
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
#process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')
# Note, must be the same as input collection used for electron ntuplizer
process.egmGsfElectronIDs.physicsObjectSrc = process.rootTupleElectrons.InputTag
from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)
# Define which IDs we want to produce
# Each of these two example IDs contains all four standard
# cut-based ID working points
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V0_miniAOD_cff']
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V1_miniAOD_cff')
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff')
#Add them to the VID producer
for idmod in my_id_modules:
  setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


#----------------------------------------------------------------------------------------------------
# Turn on trigger matching
#----------------------------------------------------------------------------------------------------
# stored in pat::TriggerObjectStandAlone with default MiniAOD config
# See: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
# To access, see: https://hypernews.cern.ch/HyperNews/CMS/get/physTools/3286/1/2.html
# This is based off of the hypernews post
# FIXME check that we match to the paths we want while we're updating the cpp code
from PhysicsTools.PatAlgos.slimming.unpackedPatTrigger_cfi import unpackedPatTrigger
# need to load it into the process, or else it won't run
process.unpackedPatTrigger = unpackedPatTrigger.clone()

process.cleanElectronTriggerMatchHLTSingleElectronWP85.matched = 'unpackedPatTrigger'
process.cleanElectronTriggerMatchHLTDoubleElectron.matched = 'unpackedPatTrigger'
process.cleanElectronTriggerMatchHLTElectronJetJet.matched = 'unpackedPatTrigger'
process.pfjetTriggerMatchHLTEleJetJet.matched = 'unpackedPatTrigger'

#process.out.outputCommands += ['keep *_electronsTriggerMatchHLTSingleElectron_*_*',
#                               'keep *_electronsTriggerMatchHLTSingleElectronWP80 _*_*',
#                               'keep *_electronsTriggerMatchHLTDoubleElectron_*_*'      ]
## select the matching electrons and keep them
#process.out.outputCommands += ['keep *_ electronsTriggeredHLTSingleElectron_*_*',
#                               'keep *_ electronsTriggeredHLTSingleElectronWP80_*_*',
#                               'keep *_ electronsTriggeredHLTDoubleElectron_*_*' ]

#----------------------------------------------------------------------------------------------------
# Add PFMET and TCMET
# See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#MET_Tools
#----------------------------------------------------------------------------------------------------
# SIC: I think we can get away with using the default MiniAOD MET
# Reproduce "raw" MET from packedPFCandidates
from RecoMET.METProducers.PFMET_cfi import pfMet
process.pfMet = pfMet.clone(src = "packedPFCandidates")
process.pfMet.calculateSignificance = False # this can't be easily implemented on packed PF candidates at the moment (as of Feb 17 2015)

#----------------------------------------------------------------------------------------------------
# MET filters
#----------------------------------------------------------------------------------------------------
# SIC: A number of filters are run by default:
#Flag_HBHENoiseFilter = cms.Path(HBHENoiseFilter)
#Flag_CSCTightHaloFilter = cms.Path(CSCTightHaloFilter)
#Flag_hcalLaserEventFilter = cms.Path(hcalLaserEventFilter)
#Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(EcalDeadCellTriggerPrimitiveFilter)
#Flag_goodVertices = cms.Path(goodVertices)
#Flag_trackingFailureFilter = cms.Path(goodVertices + trackingFailureFilter)
#Flag_eeBadScFilter = cms.Path(eeBadScFilter)
#Flag_ecalLaserCorrFilter = cms.Path(ecalLaserCorrFilter)
#Flag_trkPOGFilters = cms.Path(trkPOGFilters)
#Flag_trkPOG_manystripclus53X = cms.Path(~manystripclus53X)
#Flag_trkPOG_toomanystripclus53X = cms.Path(~toomanystripclus53X)
#Flag_trkPOG_logErrorTooManyClusters = cms.Path(~logErrorTooManyClusters)
# See: https://github.com/cms-sw/cmssw/blob/CMSSW_7_2_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
# They are stored in edm::TriggerResults of the PAT process, and can be checked the same way as HLT paths.


#----------------------------------------------------------------------------------------------------
# Use default MiniAOD Taus
# Use default MiniAOD Tau IDs
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Make analysisPatTaus and add them to the cleanPatCandidates sequence
#----------------------------------------------------------------------------------------------------
#process.analysisPatTaus = process.cleanPatTaus.clone()
#process.analysisPatTaus.preselection = cms.string(
#    'tauID("decayModeFinding") > 0.5 &'
#    ' tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 &'
#    ' tauID("againstMuonLoose3") > 0.5 &'
#    ' tauID("againstElectronLooseMVA3") > 0.5'
#)
#process.analysisPatTaus.finalCut = cms.string('pt > 20. & abs(eta) < 2.3')
#process.cleanPatCandidates.replace ( process.cleanPatTaus, process.cleanPatTaus + process.analysisPatTaus )
# FIXME

#----------------------------------------------------------------------------------------------------
# Make analysisPatMuons and add them to the cleanPatCandidates sequence
#----------------------------------------------------------------------------------------------------
#process.analysisPatMuons = process.cleanPatMuons.clone()
#process.analysisPatMuons.finalCut = cms.string("isGlobalMuon & muonID('GlobalMuonPromptTight') & pt > 20")
#process.cleanPatCandidates.replace ( process.cleanPatMuons, process.cleanPatMuons + process.analysisPatMuons )
#FIXME

#----------------------------------------------------------------------------------------------------
# Make analysisPatElectrons and add them to the cleanPatCandidates sequence
#----------------------------------------------------------------------------------------------------
#process.analysisPatElectrons = process.cleanPatElectrons.clone()
#process.analysisPatElectrons.finalCut = cms.string('userInt("HEEPId") < 0.5')
#process.cleanPatCandidates.replace ( process.cleanPatElectrons, process.cleanPatElectrons + process.analysisPatElectrons )
# FIXME

#----------------------------------------------------------------------------------------------------
# Add MVA electron ID
#
# MVA electron ID details on this twiki:
# https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentification#MVA_based_Id_in_PAT
#
# Taken from the example:
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/EgammaAnalysis/ElectronTools/test/patTuple_electronId_cfg.py?revision=1.2&view=markup&pathrev=V00-00-09
#----------------------------------------------------------------------------------------------------
#process.mvaID = cms.Sequence(  process.mvaTrigV0 + process.mvaTrigNoIPV0 + process.mvaNonTrigV0 )
#process.patElectrons.electronIDSources.mvaTrigV0     = cms.InputTag("mvaTrigV0")  
#process.patElectrons.electronIDSources.mvaNonTrigV0  = cms.InputTag("mvaNonTrigV0") 
#process.patElectrons.electronIDSources.mvaTrigNoIPV0 = cms.InputTag("mvaTrigNoIPV0")
# Updated for 72X Run II
# must checkout private code from Hugues for the moment
## See: https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
#process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_CSA14_cfi')
#process.mvaID = cms.Sequence(  process.mvaTrigV025nsCSA14 + process.mvaNonTrigV025nsPHYS14 )
#FIXME: Add the stuff to calculate this in the analyzer as per Hugues' example
# See https://github.com/HuguesBrun/ExampleElectronMVAid/blob/master/plugins/ExampleElectronMVAid.cc

#process.patConversions = cms.EDProducer("PATConversionProducer",
#    electronSource = cms.InputTag("cleanPatElectrons")  
#)
# XXX FIXME still needed?

#----------------------------------------------------------------------------------------------------
# Add the PFJets
# See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#Jet_Tools
# 
# Recommended JEC for PFJets:
# - L1FastJet : https://twiki.cern.ch/twiki/bin/view/CMS/JECAnalysesRecommendations#Corrections
# - L2Relative: https://twiki.cern.ch/twiki/bin/view/CMS/IntroToJEC#Mandatory_Jet_Energy_Corrections
# - L3Absolute: https://twiki.cern.ch/twiki/bin/view/CMS/IntroToJEC#Mandatory_Jet_Energy_Corrections
#----------------------------------------------------------------------------------------------------
#process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")
# MiniAOD default: ak4PFJetsCHS with Pt > 10 GeV
#                  L1+L2+L3+residual corrections are applied;
#                  b-tagging and pileup jet id information are embedded.
#                  Links are provided to the constituent PF candidates. 
# FIXME: Do we need to rebuild AK5 jets? Or can we stick with AK4?
# AK4 will possibly have more support (for JECs, etc.)
process.load('Leptoquarks.RootTupleMakerV2.ak5pfjets_cfi')
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
addJetCollection(process, labelName = 'AK5',
                 jetSource = cms.InputTag('ak5PFJets'),
                 algo = 'AK', rParam = 0.5,
                 jetCorrections = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
                 #btagInfos = ['caTopTagInfos']
                 )
process.patJetsAK5.userData.userFloats.src = [] # start with empty list of user floats
process.selectedPatJetsAK5.cut = cms.string("pt > 10")
process.patJetGenJetMatchAK5.matched =  'slimmedGenJets'
process.patJetPartonMatchAK5.matched = 'prunedGenParticles'
process.patJetPartons.particles = 'prunedGenParticles'
process.patJetPartonsLegacy.src = 'prunedGenParticles'
process.patJetCorrFactorsAK5.primaryVertices = 'offlineSlimmedPrimaryVertices'
# needed when using btagInfos?
#process.jetTracksAssociatorAtVertexAK5.tracks = 'unpackedTracksAndVertices'
#process.jetTracksAssociatorAtVertexAK5.pvSrc = 'offlineSlimmedPrimaryVertices'

#process.out.outputCommands += ['keep *_ak5GenJets_*_EX',
#                               'keep *_ak5PFJets_*_EX',
#                               'keep *_ak5PFJetsCHS_*_EX', ]

#----------------------------------------------------------------------------------------------------
# Make analysisPatJets and add them to the patDefaultSequence
#----------------------------------------------------------------------------------------------------
#process.cleanPatJetsAK5PF = process.cleanPatJets.clone()
#process.cleanPatJetsAK5PF.src = cms.InputTag('patJetsAK5PF')
#process.patDefaultSequence.replace (process.cleanPatJets, process.cleanPatJets + process.cleanPatJetsAK5PF)
## FIXME
#process.analysisPatJetsAK5PF = process.cleanPatJetsAK5PF.clone()
#process.analysisPatJetsAK5PF.finalCut = cms.string("abs(eta)<2.5 & pt > 20")
#process.patDefaultSequence.replace ( process.cleanPatJetsAK5PF, process.cleanPatJetsAK5PF + process.analysisPatJetsAK5PF )
## FIXME

#----------------------------------------------------------------------------------------------------
# Add the pileup MVA to the PFJets
#----------------------------------------------------------------------------------------------------
#process.load("Leptoquarks.RootTupleMakerV2.pujetidsequence_cff")
# SIC FIXME
# This needs an update for 72X. Not sure what the official recipe is now.

#----------------------------------------------------------------------------------------------------
# No CaloJets in MiniAOD

#----------------------------------------------------------------------------------------------------
# Define the systematic shift correction
#----------------------------------------------------------------------------------------------------
process.load("JetMETCorrections.Type1MET.correctedMet_cff")
process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")
process.load("JetMETCorrections.Type1MET.correctionTermsPfMetShiftXY_cff")
# FIXME needs 72X update when available
process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_mc
# Type1 PFMET provided with MiniAOD default
# FIXME Check implementation of this

#----------------------------------------------------------------------------------------------------
# Use the runMetUncertainties tool for AK5PFJets
#----------------------------------------------------------------------------------------------------
from PhysicsTools.PatUtils.tools.runType1PFMEtUncertainties import runType1PFMEtUncertainties
addJetCollection(process, postfix   = "ForMetUnc", labelName = 'AK5PF', jetSource = cms.InputTag('ak5PFJets'), jetCorrections = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], ''))
process.patJetsAK5PFForMetUnc.getJetMCFlavour = False
runType1PFMEtUncertainties(process,
                           addToPatDefaultSequence=False,
                           jetCollection="selectedPatJetsAK5PFForMetUnc",
                           electronCollection="slimmedElectrons",
                           muonCollection="slimmedMuons",
                           tauCollection="slimmedTaus",
                           makeType1p2corrPFMEt=True,
                           outputModule=None)
process.patJetPartonMatchAK5PFForMetUnc.matched = 'prunedGenParticles'
process.patJetPartonsForMetUnc.particles = 'prunedGenParticles'
process.patJetPartonsLegacyForMetUnc.src = 'prunedGenParticles'
process.patJetCorrFactorsAK5PFForMetUnc.primaryVertices = 'offlineSlimmedPrimaryVertices'

#----------------------------------------------------------------------------------------------------
# Available pat::MET collections for analysis
# - process.patMETsTC                               : raw        TCMET   (NO  jet smearing)
# 
# - process.patMETsRawCalo                          : raw        CaloMET (NO  jet smearing)
# - process.patMETs                                 : Type1      CaloMET (NO  jet smearing)
# 
# - process.patMETsRawPF                            : raw        PFMET   (NO  jet smearing)
# - process.patType1CorrectedPFMet_Type1Only        : Type1      PFMET   (YES jet smearing)
# - process.patType1CorrectedPFMet_Type01Only       : Type0+1    PFMET   (YES jet smearing)
# - process.patType1CorrectedPFMet                  : Type0+1+XY PFMET   (YES jet smearing) <-- Recommended for analysis
# 
# Available pat::MET collections for systematic studies
# - process.patType1CorrectedPFMetElectronEnUp      : Type0+1+XY PFMET   (YES jet smearing), Electron energy shifted up 
# - process.patType1CorrectedPFMetElectronEnDown    : Type0+1+XY PFMET   (YES jet smearing), Electron energy shifted down
# - process.patType1CorrectedPFMetMuonEnUp          : Type0+1+XY PFMET   (YES jet smearing), Muon energy shifted up
# - process.patType1CorrectedPFMetMuonEnDown        : Type0+1+XY PFMET   (YES jet smearing), Muon energy shifted down 
# - process.patType1CorrectedPFMetTauEnUp           : Type0+1+XY PFMET   (YES jet smearing), Tau energy shifted up   
# - process.patType1CorrectedPFMetTauEnDown         : Type0+1+XY PFMET   (YES jet smearing), Tau energy shifted down 
# - process.patType1CorrectedPFMetJetResUp          : Type0+1+XY PFMET   (YES jet smearing), Jet resolution smeared up
# - process.patType1CorrectedPFMetJetResDown        : Type0+1+XY PFMET   (YES jet smearing), Jet resolution smeared down
# - process.patType1CorrectedPFMetJetEnUp           : Type0+1+XY PFMET   (YES jet smearing), Jet energy shifted up   
# - process.patType1CorrectedPFMetJetEnDown         : Type0+1+XY PFMET   (YES jet smearing), Jet energy shifted down  
# - process.patType1CorrectedPFMetUnclusteredEnUp   : Type0+1+XY PFMET   (YES jet smearing), Unclustered energy shifted up  
# - process.patType1CorrectedPFMetUnclusteredEnDown : Type0+1+XY PFMET   (YES jet smearing), Unclustered energy shifted down
# 
# Available shifted object collections:
# - process.shiftedPatElectronsEnUp                 : pat electrons, energy scale shifted up
# - process.shiftedPatElectronsEnDown               : pat electrons, energy scale shifted down
# - process.shiftedPatMuonsEnUp                     : pat muons    , energy scale shifted up
# - process.shiftedPatMuonsEnDown                   : pat muons    , energy scale shifted down
# - process.shiftedPatTausEnUp                      : pat taus     , energy scale shifted up
# - process.shiftedPatTausEnDown                    : pat taus     , energy scale shifted down
# - process.smearedPatJetsAK5PF                     : pat jets     , energy resolution smeared to match data  <-- Recommended for analysis
# - process.smearedPatJetsAK5PFresUp                : pat jets     , energy resolution smeared worse data
# - process.smearedPatJetsAK5PFresDown              : pat jets     , energy resolution sharpened better than data
# - process.shiftedPatJetsAK5PFenUpForCorrMEt       : pat jets     , energy scale shifted up  
# - process.shiftedPatJetsAK5PFenDownForCorrMEt     : pat jets     , energy scale shifted down 
# 
# Notes:
# - No Type0 corrections are available for CaloMET
# - No Type0 or Type1 corrections are available for TCMET
# - No Type2 corrections recommended for any MET, since they degrade the MET resolution
#
# Thanks to Jared Sturdy, Christian Veelken, and Guler Karapinar
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Add the raw CaloMET and raw PFMET
#----------------------------------------------------------------------------------------------------
# SIC XXX no raw calo met

# PFMET: raw
# PFMET: Type1, with jet smearing
# XXX SIC now called: pfMetT1
# PFMET: Type0+1, with jet smearing
# FIXME I think all of these can be extracted from the pat::MET object
process.load('JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi')
process.pfMetT1Txy = cms.EDProducer(
  "AddCorrectionsToPFMET",
  src = cms.InputTag('pfMet'),
  srcCorrections = cms.VInputTag(
  cms.InputTag('pfMETcorrType0'),
  cms.InputTag('corrPfMetType1', 'type1'),
  ),
) 


#----------------------------------------------------------------------------------------------------
# This is MC, so analyze the smeared PFJets by default
# But smearing should be done by default in MiniAOD ?
#----------------------------------------------------------------------------------------------------
# SIC FIXME TODO
# analyze the default jet collection in MiniAOD (AK4PFJetsCHS)
#process.rootTuplePFJets.InputTag = cms.InputTag('slimmedJets')
# analyze the remade AK5 PFJets (we also have the AK5PFJetsCHS available)
#process.rootTuplePFJets.InputTag = 'patJetsAK5'
#process.rootTuplePFJets.InputTag = 'pfjetTriggerMatchEmbedderHLTEleJetJet'
#process.pfjetTriggerMatchHLTEleJetJet.src = 'patJetsAK5' # use for trigger matching
#process.pfjetTriggerMatchEmbedderHLTEleJetJet.src = 'patJetsAK5' # also use for trigger match embedding

#process.rootTuplePFJets.InputTag = 'pfjetTriggerMatchEmbedderHLTEleJetJet'
#process.rootTuplePFJets.InputTag = cms.InputTag('smearedAnalysisPatJetsAK5PF')
#process.rootTuplePFJets.InputTagSmearedUp   = cms.InputTag('smearedAnalysisPatJetsAK5PFresUp')                                 
#process.rootTuplePFJets.InputTagSmearedDown = cms.InputTag('smearedAnalysisPatJetsAK5PFresDown')                                 
#process.rootTuplePFJets.InputTagScaledUp    = cms.InputTag('shiftedAnalysisPatJetsAK5PFenUpForCorrMEt')                                 
#process.rootTuplePFJets.InputTagScaledDown  = cms.InputTag('shiftedAnalysisPatJetsAK5PFenDownForCorrMEt')     

#----------------------------------------------------------------------------------------------------
# Set Lepton-Gen Matching Parameters
#----------------------------------------------------------------------------------------------------
# XXX FIXME SIC: Needed? At least some matching is done by default in PAT/MiniAOD
process.load("Leptoquarks.RootTupleMakerV2.leptonGenMatching_cfi")
#process.patDefaultSequence.replace( process.electronMatch, process.elMatch )
#process.patElectrons.genParticleMatch = cms.VInputTag( cms.InputTag("elMatch") )
#process.patDefaultSequence.replace( process.muonMatch, process.muMatch )
#process.patMuons.genParticleMatch = cms.VInputTag( cms.InputTag("muMatch") )
#process.patDefaultSequence.replace( process.tauMatch, process.tauLepMatch )
#process.patTaus.genParticleMatch = cms.VInputTag( cms.InputTag("tauLepMatch") )
#process.patDefaultSequence.replace( process.tauGenJetMatch, process.tauJetMatch )
#process.patTaus.genJetMatch = cms.InputTag("tauJetMatch")

#----------------------------------------------------------------------------------------------------
# Lepton + Jets filter
#----------------------------------------------------------------------------------------------------

process.load("Leptoquarks.LeptonJetFilter.leptonjetfilter_cfi")

#### Shared Muon/Electron/Tau Skim
process.LJFilter.tauLabel  = cms.InputTag("slimmedTaus")                        
process.LJFilter.muLabel   = cms.InputTag("slimmedMuons")
process.LJFilter.elecLabel = cms.InputTag("slimmedElectrons")
process.LJFilter.jetLabel  = cms.InputTag("slimmedJets")
process.LJFilter.muonsMin = 0
process.LJFilter.muPT     = 10.0
process.LJFilter.electronsMin = 0
process.LJFilter.elecPT       = 15.0
process.LJFilter.tausMin = 0
process.LJFilter.tauPT   = 15.0
process.LJFilter.jetsMin = 0
process.LJFilter.jetPT   = 15.0
process.LJFilter.counteitherleptontype = True
process.LJFilter.customfilterEMuTauJet2012 = True
# -- WARNING :
# "customfilterEMuTauJet2012" configuration is hard-coded.
# (see: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Leptoquarks/LeptonJetFilter/src/LeptonJetFilter.cc )
# "customfilterEMuTauJet2012" is the desired mode of operation for the Lepton+Jets Filter in 2012.

#----------------------------------------------------------------------------------------------------
# PDF weights
#----------------------------------------------------------------------------------------------------

process.pdfWeights = cms.EDProducer("PdfWeightProducer",
	# Fix POWHEG if buggy (this PDF set will also appear on output,
	# so only two more PDF sets can be added in PdfSetNames if not "")
	#FixPOWHEG = cms.untracked.string("CT10.LHgrid"),
        # GenTag = cms.untracked.InputTag("genParticles"),
	PdfInfoTag = cms.untracked.InputTag("generator"),
	PdfSetNames = cms.untracked.vstring(
            "CT10.LHgrid" , 
            "MSTW2008nlo68cl.LHgrid",
            "NNPDF20_100.LHgrid"
	)
)

#----------------------------------------------------------------------------------------------------
# Define the output tree for RootTupleMakerV2
#----------------------------------------------------------------------------------------------------

# RootTupleMakerV2 tree
process.rootTupleTree = cms.EDAnalyzer("RootTupleMakerV2_Tree",
    outputCommands = cms.untracked.vstring(
        'drop *',
        # Event information
        'keep *_rootTupleEvent_*_*',
        'keep *_rootTupleEventSelection_*_*',
        # Single objects
        'keep *_rootTuplePFCandidates_*_*',
        'keep *_rootTuplePFJets_*_*',
        'keep *_rootTupleElectrons_*_*',
        # FIXME ignore for now
        'keep *_rootTupleMuons_*_*',
        #'keep *_rootTupleHPSTaus_*_*',
        # FIXME ignore for now
        #'keep *_rootTuplePhotons_*_*',
        'keep *_rootTupleVertex_*_*',
        # FIXME ignore for now
        ## MET objects for analysis
        #'keep *_rootTupleTCMET_*_*',
        #'keep *_rootTupleCaloMET_*_*',
        #'keep *_rootTupleCaloMETType1Cor_*_*',
        #'keep *_rootTuplePFMET_*_*',
        #'keep *_rootTuplePFMETType1Cor_*_*',
        #'keep *_rootTuplePFMETType01Cor_*_*',
        #'keep *_rootTuplePFMETType01XYCor_*_*',
        #'keep *_rootTuplePFMETType01XYCor_*_*',
        # pdf weights
        #'keep *_rootTuplePFMETType01XYCorUnclusteredUp_*_*',
        #'keep *_rootTuplePFMETType01XYCorUnclusteredDown_*_*',
        #'keep *_rootTuplePFMETType01XYCorElectronEnUp_*_*',
        #'keep *_rootTuplePFMETType01XYCorElectronEnDown_*_*',
        #'keep *_rootTuplePFMETType01XYCorMuonEnUp_*_*',
        #'keep *_rootTuplePFMETType01XYCorMuonEnDown_*_*',
        #'keep *_rootTuplePFMETType01XYCorTauEnUp_*_*',
        #'keep *_rootTuplePFMETType01XYCorTauEnDown_*_*',
        #'keep *_rootTuplePFMETType01XYCorJetResUp_*_*',
        #'keep *_rootTuplePFMETType01XYCorJetResDown_*_*',
        #'keep *_rootTuplePFMETType01XYCorJetEnUp_*_*',
        #'keep *_rootTuplePFMETType01XYCorJetEnDown_*_*',
        # Trigger objects
        'keep *_rootTupleTrigger_*_*',
        # FIXME ignore for now
        #'keep *_rootTupleTriggerObjects_*_*',
        # GEN objects
        'keep *_rootTupleGenEventInfo_*_*',
        'keep *_rootTupleGenParticles_*_*',
        'keep *_rootTupleGenJets_*_*',
        'keep *_rootTupleGenElectronsFromWs_*_*',
        'keep *_rootTupleGenElectronsFromZs_*_*',
        'keep *_rootTupleGenMuonsFromWs_*_*',
        'keep *_rootTupleGenMuonsFromZs_*_*',
        'keep *_rootTupleGenTausFromWs_*_*',
        'keep *_rootTupleGenTausFromZs_*_*',
        # FIXME ignore for now
        #'keep *_rootTupleGenMETTrue_*_*',
        #'keep *_rootTupleGenMETCalo_*_*'       
    )
)

#----------------------------------------------------------------------------------------------------
# Define GEN particle skimmer modules
#----------------------------------------------------------------------------------------------------

process.load ('Leptoquarks.LeptonJetGenTools.genTauMuElFromZs_cfi')
process.load ('Leptoquarks.LeptonJetGenTools.genTauMuElFromWs_cfi') 

#----------------------------------------------------------------------------------------------------
# Define the path 
#----------------------------------------------------------------------------------------------------
# to see EventSetup content
#process.esContent = cms.EDAnalyzer("PrintEventSetupContent")
# to see Event content
process.load('FWCore.Modules.printContent_cfi')

process.p = cms.Path(
    # extra Phys14 VIDs (inc. HEEP)
    process.egmGsfElectronIDSequence*
    # ak5 jets
    process.ak5PFJetsCHS*
    process.ak5GenJets*
    process.ak5PFJets*
    process.patJetsAK5*
    process.patJetsAK5PFForMetUnc*
    process.patJetsAK5PFForMetUncNotOverlappingWithLeptonsForJetMEtUncertainty*
    # MVA electron ID (can't be run on MiniAOD?)
    #process.mvaID*
    # gen particle skimmer modules
    process.genTausFromWs*
    process.genMuonsFromWs*
    process.genElectronsFromWs*
    process.genTausFromZs*
    process.genMuonsFromZs*
    process.genElectronsFromZs*
    # FIXME ADD PDF WEIGHTS BACK LATER?
    # pdf weights
    process.pdfWeights*
    # Now the regular PAT default sequence
    #process.patDefaultSequence*
    # Add the pileup MVA to the jets
    #FIXME process.puJetIdSequence*
    # MET producers (included by default)
    # process.patMETsRawCalo*
    #process.patMETsRawPF*
    #process.pfMetT1*
    #process.patType1CorrectedPFMetType01Only*
    # L+J Filter
    process.LJFilter*  
    # Run PAT conversions for electrons (no longer needed, I assume)
    #process.patConversions*
    # Re-run full HPS sequence to fully profit from the fix of high pT taus (no longer needed, I assume)
    #process.recoTauClassicHPSSequence*
    #process.esContent*
    # RootTupleMakerV2
    (
    # Event information
    process.rootTupleEvent+
    process.rootTupleEventSelection+
    # Single objects
    process.rootTuplePFCandidates+
    #process.printContent*
    process.rootTuplePFJets+
    process.rootTupleElectrons+
    process.rootTupleMuons+
    # FIXME ignore for now
    #process.rootTupleHPSTaus+
    #process.rootTuplePhotons+
    process.rootTupleVertex+
    ## MET objects for analysis
    ##process.rootTupleTCMET+ # no longer made?
    # FIXME ignore for now
    #process.rootTupleCaloMET+
    #process.rootTupleCaloMETType1Cor+
    #process.rootTuplePFMET+
    #process.rootTuplePFMETType1Cor+
    #process.rootTuplePFMETType01Cor+
    #process.rootTuplePFMETType01XYCor+
    ## MET objects for systematics
    #process.rootTuplePFMETType01XYCorUnclusteredUp+
    #process.rootTuplePFMETType01XYCorUnclusteredDown+
    #process.rootTuplePFMETType01XYCorElectronEnUp+
    #process.rootTuplePFMETType01XYCorElectronEnDown+
    #process.rootTuplePFMETType01XYCorMuonEnUp+
    #process.rootTuplePFMETType01XYCorMuonEnDown+
    #process.rootTuplePFMETType01XYCorTauEnUp+
    #process.rootTuplePFMETType01XYCorTauEnDown+
    #process.rootTuplePFMETType01XYCorJetResUp+
    #process.rootTuplePFMETType01XYCorJetResDown+
    #process.rootTuplePFMETType01XYCorJetEnUp+
    #process.rootTuplePFMETType01XYCorJetEnDown+
    # Trigger objects
    process.rootTupleTrigger+
    #FIXME LATER?
    #process.rootTupleTriggerObjects+
    # GEN objects
    process.rootTupleGenEventInfo+
    process.rootTupleGenParticles+
    process.rootTupleGenJets+
    process.rootTupleGenElectronsFromWs+
    process.rootTupleGenElectronsFromZs+
    process.rootTupleGenMuonsFromWs+
    process.rootTupleGenMuonsFromZs+
    process.rootTupleGenTausFromWs+
    process.rootTupleGenTausFromZs
    #+
    # FIXME ignore for now
    #process.rootTupleGenMETTrue+
    #process.rootTupleGenMETCalo
    )*
    # Put everything into the tree
    process.rootTupleTree
)


#----------------------------------------------------------------------------------------------------
# Dump if necessary
#----------------------------------------------------------------------------------------------------

#process.dump = cms.OutputModule("PoolOutputModule",
#                                outputCommands = cms.untracked.vstring(
#                                'keep *',
#                                ),
#                                fileName = cms.untracked.string('dump.root')
#                                )
#process.DUMP    = cms.EndPath (process.dump)

# Delete predefined Endpath (needed for running with CRAB)
del process.out
del process.outpath

#----------------------------------------------------------------------------------------------------
# Run the path
#----------------------------------------------------------------------------------------------------

process.schedule = cms.Schedule(process.p)#,process.DUMP)

#print process.dumpPython()
# try to convert it to unscheduled explicitly?
#from FWCore.ParameterSet.Utilities import convertToUnscheduled
#convertToUnscheduled(process)
