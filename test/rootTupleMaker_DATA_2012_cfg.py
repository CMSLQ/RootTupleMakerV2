#----------------------------------------------------------------------------------------------------
# Load PAT template
#----------------------------------------------------------------------------------------------------


# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

# load the coreTools of PAT
from PhysicsTools.PatAlgos.tools.coreTools import *

# Remove MC matching for data analysis 
removeMCMatching(process, ['All'])

import os

############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.default.limit = 10
#################################################################

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

# Make sure a correct global tag is used (please refer to https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Valid_Global_Tags_by_Release)
process.GlobalTag.globaltag = 'GR_R_52_V7::All'

# Events to process
process.maxEvents.input = 10

# Options and Output Report
process.options.wantSummary = False

# Input files
process.source.fileNames = [
    'file:///afs/cern.ch/user/e/eberry/work/Run2012B_ElectronHad_AOD_PromptReco-v1_TEST.root'
    #rfio:///castor/cern.ch/user/h/hsaka/2012prep/Run2012B_ElectronHad_AOD_PromptReco-v1_TEST.root'
]

#----------------------------------------------------------------------------------------------------
# Add PFMET and TCMET
#----------------------------------------------------------------------------------------------------

from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, 'PF')
addTcMET(process, 'TC')

#----------------------------------------------------------------------------------------------------
# Add PFMET corrections:
# - Type 0 : PU corrections
# - Type 1 : JES corrections
# - MET x/y shift correction: phi modulation
# - See: 
#----------------------------------------------------------------------------------------------------

# type-0 MET corrections
process.pfType1CorrectedMet.applyType0Corrections = cms.bool(True)

# type-1 MET corrections (see jets)

# MET x/y shift correction: phi modulation
process.load("PhysicsTools.PatUtils.patPFMETCorrections_cff")
process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")

process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_data

process.pfType1CorrectedMet.applyType0Corrections = cms.bool(True)

process.patType1CorrectedPFMet.srcType1Corrections = cms.VInputTag(
   cms.InputTag('patPFJetMETtype1p2Corr', 'type1'),
   cms.InputTag('patPFMETtype0Corr'),
   cms.InputTag('pfMEtSysShiftCorr')
)

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
process.patTaus.tauIDSources.byVLooseIsolation = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolation")
process.patTaus.tauIDSources.byLooseIsolation  = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation")
process.patTaus.tauIDSources.byMediumIsolation = cms.InputTag("hpsPFTauDiscriminationByMediumIsolation")
process.patTaus.tauIDSources.byTightIsolation  = cms.InputTag("hpsPFTauDiscriminationByTightIsolation")


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

# Skim definition
# process.load("Leptoquarks.LeptonJetFilter.leptonjetfilter_cfi")
##################################################################
#### Electron based skim
# process.LJFilter.muLabel = 'muons'
# process.LJFilter.elecLabel = 'gsfElectrons'
# process.LJFilter.jetLabel = 'ak5CaloJets'
# process.LJFilter.tauLabel = 'shrinkingConePFTauProducer'
# process.LJFilter.muonsMin = -1
# process.LJFilter.electronsMin = 1
# process.LJFilter.elecPT = 20.
# process.LJFilter.counteitherleptontype = False
##################################################################
#### Photon based skim
# process.LJFilter.muLabel = 'muons'
# process.LJFilter.elecLabel = 'gsfElectrons'
# process.LJFilter.photLabel = 'photons'
# process.LJFilter.jetLabel = 'ak5CaloJets'
# process.LJFilter.tauLabel = 'shrinkingConePFTauProducer'
# process.LJFilter.muonsMin = -1
# process.LJFilter.electronsMin = -1
# process.LJFilter.photMin = 1
# process.LJFilter.photET = 20.
# process.LJFilter.photHoE = 0.05
##################################################################


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
        'keep *_rootTuplePFMETType1Cor_*_*',
        'keep *_rootTuplePFChargedMET_*_*',
        'keep *_rootTupleMuons_*_*',
        'keep *_rootTupleTrigger_*_*',
        'keep *_rootTupleVertex_*_*',
        'keep *_rootTupleGenEventInfo_*_*',
        'keep *_rootTupleGenParticles_*_*',
        'keep *_rootTupleGenJets_*_*',
        'keep *_rootTupleGenMETTrue_*_*',
        'keep *_rootTupleGenMETCalo_*_*',       
        'keep *_rootTuplePhotons_*_*',
        'keep *_rootTuplePFCandidates_*_*'
    )
)

# Path definition
process.p = cms.Path(
    process.goodVertices*
    process.type0PFMEtCorrection*
    # MET filters (required):
    process.CSCTightHaloFilter*
    process.EcalDeadCellTriggerPrimitiveFilter*
    process.HBHENoiseFilter*
    process.HBHENoiseFilterResultProducer*
    process.hcalLaserEventFilter*
    process.trackingFailureFilter*
    # MVA electron ID
    process.mvaID + 
    # PAT sequence
    process.patDefaultSequence*
    # RootTupleMakerV2
    (
    process.rootTupleEvent+
    process.rootTupleEventSelection+
    process.rootTuplePFJets+
    process.rootTupleElectrons+
    # ---
    #process.rootTupleTaus+
    process.rootTupleSCTaus+
    process.rootTupleHPSTaus+
    # ---
    process.rootTupleCaloMET+
    process.rootTupleTCMET+
    process.rootTuplePFMET+
    # process.rootTuplePFMETType1Cor+
    # process.rootTuplePFChargedMET+
    process.rootTupleMuons+
    process.rootTupleTrigger+
    process.rootTupleVertex+
    process.rootTupleGenEventInfo+
    process.rootTupleGenParticles+
    process.rootTupleGenJets+
    process.rootTupleGenMETTrue+
    process.rootTupleGenMETCalo+    
    process.rootTuplePhotons
    # process.rootTuplePFCandidates
    )
    *process.rootTupleTree
)

##################
process.dump = cms.OutputModule("PoolOutputModule",
                                outputCommands = cms.untracked.vstring(
                                'keep *',
                                ),
                                fileName = cms.untracked.string('dump.root')
                                )
process.DUMP    = cms.EndPath (process.dump)
##################

# Delete predefined Endpath (needed for running with CRAB)
del process.out
del process.outpath

# Schedule definition
process.schedule = cms.Schedule(process.p)

##############
process.schedule.append( process.DUMP )
##############
