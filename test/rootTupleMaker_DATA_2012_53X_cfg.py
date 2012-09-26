#----------------------------------------------------------------------------------------------------
# Load PAT template
#----------------------------------------------------------------------------------------------------

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

# Load PF isolation for muons and electrons
from PhysicsTools.PatAlgos.tools.pfTools import *
usePFIso ( process )

# Remove MC matching *and* use the right JECs for data analysis 
from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData( process )

# Options and Output Report
process.options.wantSummary = True

import os

############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.default.limit = 10
#################################################################

#----------------------------------------------------------------------------------------------------
# The ECAL laser correction filter (ecalLaserCorrFilter) occasionally (~1/20 events) 
# interpolates laser correction values of < 1 and issues a LogError message
# 
# This does not affect the laser correction that is applied: only the interpolated estimate 
# that the filter uses.  The filter runs in "Tagging Mode", so no events can be removed.
# 
# Error message comes from line 173 of:
# CalibCalorimetry/EcalLaserCorrection/src/EcalLaserDbService.cc
# function = EcalLaserDbService::getLaserCorrection
# Message = "The interpolated laser correction is <= zero!"
#
# We suppress these messages.  Suppression can be removed by commenting the following line.
#----------------------------------------------------------------------------------------------------

process.MessageLogger.suppressError = cms.untracked.vstring ('ecalLaserCorrFilter') 

#----------------------------------------------------------------------------------------------------
# Load our RootTupleMakerV2 modules
#----------------------------------------------------------------------------------------------------

process.load('Leptoquarks.RootTupleMakerV2.Ntuple_cff')

# Output ROOT file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('RootTupleMakerV2_output_DATA.root')
)

#----------------------------------------------------------------------------------------------------
# Set global settings (number of events, global tag, input files, etc)
#----------------------------------------------------------------------------------------------------
# Make sure a correct global tag is used.
#
# Correct global tags for Summer12_DR53X MC are here: 
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#2012_MC_production
#
# Correct global tags for 2012 PromptReco data are here:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Global_tags_for_prompt_reco_data
# 
# Correct global tags for 2012 ReReco data are here:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Summary_of_Global_Tags_used_in_o
#
# Recommendations as of Sept. 25, 2012:
# Summer12_DR53X  MC               (*/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM) : START53_V7F
# 2012 A          Data ReReco      (*/Run2012A-recover-06Aug2012-v1/AOD)           : FT_53_V6C_AN2
# 2012 A          Data ReReco      (*/Run2012A-13Jul2012-v1/AOD)                   : FT_53_V6_AN2
# 2012 B          Data ReReco      (*/Run2012B-13Jul2012-v1/AOD)                   : FT_53_V6_AN2
# 2012 C          Data ReReco      (*/Run2012C-24Aug2012-v1/AOD)                   : FT_53_V10_AN2
# 2012 C          Prompt Reco      (*/Run2012C-PromptReco-v2/AOD)                  : GR_P_V41_AN2
# 2012 D          Prompt Reco      Coming soon!                                    : GR_P_V42_AN2
# 
# Make sure a correct global tag is used:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Valid_Global_Tags_by_Release

process.GlobalTag.globaltag = 'GR_P_V41_AN2::All'

# Events to process
process.maxEvents.input = 10

# Input files
process.source.fileNames = [
    #'file:///afs/cern.ch/user/e/eberry/work/ZprimePSIToEE_M-2000_TuneZ2star_8TeV-pythia6_TEST.root'
    # 'file:///afs/cern.ch/user/e/eberry/work/Run2012B_ElectronHad_AOD_PromptReco-v1_TEST.root'
    #'file:///afs/cern.ch/user/e/eberry/9EE168D0-71FB-E111-B552-BCAEC5329717.root'
    'root://eoscms//eos/cms/store/user/hsaka/2012prep/Run2012C_SingleElectron_AOD_PromptReco-v2_TEST.root'
    #'rfio:///castor/cern.ch/user/h/hsaka/2012prep/Run2012B_ElectronHad_AOD_PromptReco-v1_TEST.root'
]

#----------------------------------------------------------------------------------------------------
# For 53x Data and MC, the default Jet Probability Calibration from the GlobalTag is not optimal 
# and needs to be replaced in the following way, when using CRAB:
# 
# Hypernews: https://hypernews.cern.ch/HyperNews/CMS/get/exotica/1524.html
# Twiki: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagJetProbabilityCalibration#Calibration_in_53x_Data_and_MC
#----------------------------------------------------------------------------------------------------

process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record  = cms.string          ("BTagTrackProbability2DRcd"                ),
             tag     = cms.string          ("TrackProbabilityCalibration_2D_Data53X_v2"),
             connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"    )
    ),
    cms.PSet(record  = cms.string          ("BTagTrackProbability3DRcd"                ),
             tag     = cms.string          ("TrackProbabilityCalibration_3D_Data53X_v2"),
             connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"    )
    )
)

#----------------------------------------------------------------------------------------------------
# Ensure that PF isolation for EGamma cut-based electron ID:
# - uses recommended calculation (not the one included in gsf by default)
# - uses a cone size of 0.3 (recommended)
# 
# Recommendation comes from here:
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification#Particle_Flow_Isolation
# 
# Recipe comes from here:
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolation#for_PAT_electron_users_using_sta
#
# Note: This block has to go early in the python cfg, or cmsRun won't work
#----------------------------------------------------------------------------------------------------

process.patElectrons.isolationValues = cms.PSet(
    pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03PFIdPFIso"),
    pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03PFIdPFIso"),
    pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFIdPFIso"),
    pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03PFIdPFIso"),
    pfPhotons = cms.InputTag("elPFIsoValueGamma03PFIdPFIso")
)

process.patElectrons.isolationValuesNoPFId = cms.PSet(
    pfChargedHadrons = cms.InputTag("elPFIsoValueCharged03NoPFIdPFIso"),
    pfChargedAll = cms.InputTag("elPFIsoValueChargedAll03NoPFIdPFIso"),
    pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03NoPFIdPFIso"),
    pfNeutralHadrons = cms.InputTag("elPFIsoValueNeutral03NoPFIdPFIso"),
    pfPhotons = cms.InputTag("elPFIsoValueGamma03NoPFIdPFIso")
)

#----------------------------------------------------------------------------------------------------
# HEEP 4.0 (electron ID) still uses the 2011 definitions of rho for isolation corrections.
# 
# Recipe taken from here:
# https://twiki.cern.ch/twiki/bin/view/CMS/EgammaEARhoCorrection#Rho_for_2011_Effective_Areas
#----------------------------------------------------------------------------------------------------

process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJetsForHEEPIsolation = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForHEEPIsolation.Rho_EtaMax = cms.double(2.5)

#----------------------------------------------------------------------------------------------------
# Turn on trigger matching
# 
# Example taken from PAT SWGuide twiki
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTriggerMatchExercise#TrigProduce
#
# Don't include this line (recommended by instructions): removeCleaningFromTriggerMatching( process )
# It will break the matching.
#----------------------------------------------------------------------------------------------------

# load the PAT trigger Python tools
from PhysicsTools.PatAlgos.tools.trigTools import *

# switch on the trigger matching
switchOnTriggerMatching( process, triggerMatchers = [
        # electrons 
        'cleanElectronTriggerMatchHLTSingleElectron',
        'cleanElectronTriggerMatchHLTSingleElectronWP80',
        'cleanElectronTriggerMatchHLTDoubleElectron',
        # muons
        'cleanMuonTriggerMatchHLTSingleMuon',
        'cleanMuonTriggerMatchHLTSingleIsoMuon'

] )

#----------------------------------------------------------------------------------------------------
# Add PFMET and TCMET
#----------------------------------------------------------------------------------------------------

from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')
addTcMET(process, 'TC')

#----------------------------------------------------------------------------------------------------
# Add MET filters
#----------------------------------------------------------------------------------------------------

process.load("Leptoquarks.RootTupleMakerV2.metFilters_cfi")

#----------------------------------------------------------------------------------------------------
# Add ShrinkingCone Taus
#----------------------------------------------------------------------------------------------------

from PhysicsTools.PatAlgos.tools.tauTools import *
addTauCollection(process, tauCollection = cms.InputTag('shrinkingConePFTauProducer'), algoLabel = "shrinkingCone", typeLabel = "PFTau")

#----------------------------------------------------------------------------------------------------
# Modify cleanPatTaus (HPS Taus) - loosen up a bit
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatAlgos/python/cleaningLayer1/tauCleaner_cfi.py?revision=1.11&view=markup
#----------------------------------------------------------------------------------------------------
process.cleanPatTaus.preselection = cms.string(' tauID("decayModeFinding") > 0.5 ')
process.cleanPatTaus.finalCut     = cms.string(' pt > 15.0 & abs(eta) < 2.5      ')

#----------------------------------------------------------------------------------------------------
# Add tau id sources (HPS Taus)
#----------------------------------------------------------------------------------------------------

process.load("Leptoquarks.RootTupleMakerV2.tauIDsources_cfi")

#----------------------------------------------------------------------------------------------------
# Add MVA electron ID
#
# MVA electron ID details on this twiki:
# https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentification#MVA_based_Id_in_PAT
#
# Taken from the example:
# http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/test/patTuple_electronId_cfg.py?revision=1.2&view=markup&pathrev=V00-00-21
#----------------------------------------------------------------------------------------------------

process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
process.mvaID = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
process.patElectrons.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0"   )
process.patElectrons.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")

#----------------------------------------------------------------------------------------------------
# Add PF jets --> See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#Jet_Tools
#----------------------------------------------------------------------------------------------------

# With FastJet corrections, do Type 1 MET corrections

from PhysicsTools.PatAlgos.tools.jetTools import *

addJetCollection(process,cms.InputTag('ak5PFJets'),
    'AK5', 'PF',
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'])), 
    doType1MET   = True,
    genJetCollection = cms.InputTag("ak5GenJets"),
    doJetID      = True,
)

# With Offset corrections

addJetCollection(process,cms.InputTag('ak5PFJets'),
    'AK5', 'PFL1Offset',
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5PF', cms.vstring(['L1Offset','L2Relative', 'L3Absolute','L2L3Residual'])), 
    doType1MET   = False,
    genJetCollection = cms.InputTag("ak5GenJets"),
    doJetID      = True,
)

#----------------------------------------------------------------------------------------------------
# Available pat::MET collections
# - process.patMETsTC          : raw TCMET                    (already included, no extra code needed)
# - process.patMETsRawCalo     : raw CaloMET
# - process.patMETs            : Type1-corrected CaloMET      (already included, no extra code needed)
# - process.patMETsRawPF       : raw PFMET
# - process.patMETsPF          : Type1-corrected PFMET        (already included, no extra code needed)
# - process.patMETsAK5PF       : Type0+1-corrected PFMET
# - process.patMETsAK5PFXYShift: Type0+1-corrected + XY shift-corrected PFMET
#
# No Type0 corrections are available for CaloMET
# No Type0 or Type1 corrections are available for TCMET
# Type2 corrections are not recommended, since they degrade the MET resolution
#
# Thanks to Jared Sturdy
#----------------------------------------------------------------------------------------------------

# CaloMET: raw

process.patMETsRawCalo = process.patMETsTC.clone()
process.patMETsRawCalo.metSource = cms.InputTag("met")

# PFMET: raw

process.patMETsRawPF = process.patMETsTC.clone()
process.patMETsRawPF.metSource = cms.InputTag("pfMet")

# PFMET: Type 0+1 corrections

process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")

process.AK5PFType1CorMet.applyType0Corrections = cms.bool(False)
process.AK5PFType1CorMet.srcType1Corrections   = cms.VInputTag(
    cms.InputTag('pfMETcorrType0'),
    cms.InputTag('pfJetMETcorr', 'type1')
)

# PFMET: Type 0+1 corrections + XY systematic shift correction

process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")

process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_data
process.AK5PFType1CorMetXYShift = process.AK5PFType1CorMet.clone()
process.AK5PFType1CorMetXYShift.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfMETcorrType0'),      
    cms.InputTag('pfJetMETcorr', 'type1'),
    cms.InputTag('pfMEtSysShiftCorr')    
)
process.patMETsAK5PFXYShift = process.patMETsAK5PF.clone()
process.patMETsAK5PFXYShift.metSource = cms.InputTag ("AK5PFType1CorMetXYShift")


#----------------------------------------------------------------------------------------------------
# Set Lepton-Gen Matching Parameters
#----------------------------------------------------------------------------------------------------

process.load("Leptoquarks.RootTupleMakerV2.leptonGenMatching_cfi")
process.patDefaultSequence.replace( process.electronMatch, process.elMatch )
process.patElectrons.genParticleMatch = cms.VInputTag( cms.InputTag("elMatch") )
process.patDefaultSequence.replace( process.muonMatch, process.muMatch )
process.patMuons.genParticleMatch = cms.VInputTag( cms.InputTag("muMatch") )
process.patDefaultSequence.replace( process.tauMatch, process.tauLepMatch )
process.patTaus.genParticleMatch = cms.VInputTag( cms.InputTag("tauLepMatch") )
process.patDefaultSequence.replace( process.tauGenJetMatch, process.tauJetMatch )
process.patTaus.genJetMatch = cms.InputTag("tauJetMatch")

#----------------------------------------------------------------------------------------------------
# Lepton + Jets filter
#----------------------------------------------------------------------------------------------------

process.load("Leptoquarks.LeptonJetFilter.leptonjetfilter_cfi")

#### Shared Muon/Electron/Tau Skim
process.LJFilter.tauLabel  = cms.InputTag("cleanPatTaus")                        
process.LJFilter.muLabel   = cms.InputTag("cleanPatMuons")
process.LJFilter.elecLabel = cms.InputTag("cleanPatElectrons")
process.LJFilter.jetLabel  = cms.InputTag("cleanPatJetsAK5PF")
process.LJFilter.muonsMin = 1
process.LJFilter.muPT     = 25.0
process.LJFilter.electronsMin = 1
process.LJFilter.elecPT       = 25.0
process.LJFilter.tausMin = 1
process.LJFilter.tauPT   = 25.0
process.LJFilter.counteitherleptontype = True
process.LJFilter.customfilterEMuTauJet2012 = True
# -- WARNING :
# "customfilterEMuTauJet2012" configuration is hard-coded. If enabled, other configuration parameters will NOT have any effect.
# (see: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Leptoquarks/LeptonJetFilter/src/LeptonJetFilter.cc?revision=1.12&view=markup )
# "customfilterEMuTauJet2012" is the desired mode of operation for the Lepton+Jets Filter in 2012.

#----------------------------------------------------------------------------------------------------
# Define the output tree for RootTupleMakerV2
#----------------------------------------------------------------------------------------------------

# RootTupleMakerV2 tree
process.rootTupleTree = cms.EDAnalyzer("RootTupleMakerV2_Tree",
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_rootTupleEvent_*_*',
        'keep *_rootTupleEventSelection_*_*',
        'keep *_rootTupleCaloJets_*_*',
        'keep *_rootTuplePFJets_*_*',
        'keep *_rootTupleElectrons_*_*',
        # ---
        #'keep *_rootTupleTaus_*_*',
        'keep *_rootTupleSCTaus_*_*',
        'keep *_rootTupleHPSTaus_*_*',
        # ---
        'keep *_rootTupleCaloMET_*_*',
        'keep *_rootTupleTCMET_*_*',
        'keep *_rootTuplePFMET_*_*',
        'keep *_rootTupleCaloMETType1Cor_*_*',
        'keep *_rootTuplePFMETType1Cor_*_*',
        'keep *_rootTuplePFMETType01Cor_*_*',
        'keep *_rootTuplePFMETType01XYCor_*_*',
        # 'keep *_rootTuplePFChargedMET_*_*',
        'keep *_rootTupleMuons_*_*',
        'keep *_rootTupleTrigger_*_*',
        'keep *_rootTupleTriggerObjects_*_*',
        'keep *_rootTupleVertex_*_*',
        'keep *_rootTupleGenEventInfo_*_*',
        'keep *_rootTupleGenParticles_*_*',
        'keep *_rootTupleGenTausFromLQs_*_*',
        'keep *_rootTupleGenTausFromLQTops_*_*',
        'keep *_rootTupleGenJets_*_*',
        'keep *_rootTupleGenMETTrue_*_*',
        'keep *_rootTupleGenMETCalo_*_*',       
        'keep *_rootTuplePhotons_*_*',
        'keep *_rootTuplePFCandidates_*_*'
    )
)

#----------------------------------------------------------------------------------------------------
# Define GEN particle skimmer (not used for data -- just a place-holder)
#----------------------------------------------------------------------------------------------------

process.load ('Leptoquarks.LeptonJetGenTools.genTausFromLQs_cfi')

#----------------------------------------------------------------------------------------------------
# Define the path 
#----------------------------------------------------------------------------------------------------

process.p = cms.Path(
    # gen particles
    process.genTausFromLQs*
    process.genTausFromLQTops*
    # MVA electron ID
    process.mvaID*
    # Good vertices
    process.goodVertices*
    # PFMET corrections
    process.type0PFMEtCorrection*
    process.pfMEtSysShiftCorrSequence*
    process.producePFMETCorrections*
    # PFMET producers
    process.AK5PFType1CorMet*
    process.AK5PFType1CorMetXYShift*
    # MET filters (required):
    process.EcalDeadCellTriggerPrimitiveFilter*
    process.EcalDeadCellBoundaryEnergyFilter*
    process.HBHENoiseFilterResultProducer*
    process.hcalLaserEventFilter*
    process.trackingFailureFilter*
    process.eeBadScFilter*
    process.ecalLaserCorrFilter*
    # HEEP rho for isolation correction
    process.kt6PFJetsForHEEPIsolation*
    # PAT sequence
    process.patDefaultSequence*
    # L+J filter
    process.LJFilter*    
    # PAT MET producers
    process.patMETsRawCalo*       # CaloMET: RAW
    process.patMETsRawPF*         # PFMET  : Raw
    process.patMETsAK5PF*         # PFMET  : Type 0+1 corrections
    process.patMETsAK5PFXYShift*  # PFMET  : Type 0+1 corrections, X/Y shift
    # RootTupleMakerV2
    (
    process.rootTupleEvent+
    process.rootTupleEventSelection+
    process.rootTuplePFJets+
    process.rootTupleElectrons+
    process.rootTupleSCTaus+
    process.rootTupleHPSTaus+
    process.rootTupleCaloMET+
    process.rootTupleTCMET+
    process.rootTuplePFMET+
    process.rootTupleCaloMETType1Cor+
    process.rootTuplePFMETType1Cor+
    process.rootTuplePFMETType01Cor+
    process.rootTuplePFMETType01XYCor+
    # process.rootTuplePFChargedMET+
    process.rootTupleMuons+
    process.rootTupleTrigger+
    process.rootTupleTriggerObjects+
    process.rootTupleVertex+
    process.rootTupleGenEventInfo+
    process.rootTupleGenParticles+
    process.rootTupleGenTausFromLQs+
    process.rootTupleGenTausFromLQTops+
    process.rootTupleGenJets+
    process.rootTupleGenMETTrue+
    process.rootTupleGenMETCalo+    
    process.rootTuplePhotons+
    process.rootTuplePFCandidates
    )
    *process.rootTupleTree
)

#----------------------------------------------------------------------------------------------------
# Dump if necessary
#----------------------------------------------------------------------------------------------------

# process.dump = cms.OutputModule("PoolOutputModule",
#                                 outputCommands = cms.untracked.vstring(
#                                 'keep *',
#                                 ),
#                                 fileName = cms.untracked.string('dump.root')
#                                 )
# process.DUMP    = cms.EndPath (process.dump)

# Delete predefined Endpath (needed for running with CRAB)
del process.out
del process.outpath

# Schedule definition
process.schedule = cms.Schedule(process.p)
