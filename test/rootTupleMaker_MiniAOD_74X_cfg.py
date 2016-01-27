import os
#----------------------------------------------------------------------------------------------------
# Load PAT template + customize
#----------------------------------------------------------------------------------------------------
# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from FWCore.ParameterSet.VarParsing import VarParsing
import sys

# make some options and parse
options = dict()
varOptions = VarParsing('analysis')
varOptions.register(
    "isMC",
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Customize config for MC"
)
varOptions.parseArguments()
if varOptions.isMC:
  print 'We are running on MC!'

#process.load('PhysicsTools.PatAlgos.patSequences_cff')
#process.load('Configuration.StandardSequences.Services_cff')
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
        calibratedPatElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
        ),
        calibratedElectrons = cms.PSet(
        initialSeed = cms.untracked.uint32(1),
        engineName = cms.untracked.string('TRandom3')
        ),
)
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
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.default.limit = 10
#################################################################

# We should be using PFIso by default in newer CMSSW
# see: PhysicsTools/PatAlgos/python/producersLayer1/electronProducer_cff.py

#----------------------------------------------------------------------------------------------------
# Load our RootTupleMakerV2 modules
#----------------------------------------------------------------------------------------------------
process.load('Leptoquarks.RootTupleMakerV2.Ntuple_cff')

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
# Output ROOT file
#----------------------------------------------------------------------------------------------------
process.TFileService = cms.Service("TFileService",
    fileName = cms.string( "file_m300.root" )
    #fileName = cms.string( "file_data.root" )
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
        'keep *_rootTuplePFJets*_*_*',
        'keep *_rootTupleElectrons_*_*',
        'keep *_rootTupleMuons_*_*',
        # FIXME ignore for now
        #'keep *_rootTuplePhotons_*_*',
        'keep *_rootTupleVertex_*_*',
        ## MET objects for analysis
        'keep *_rootTuplePFMET*_*_*',
        # Trigger objects
        'keep *_rootTupleTrigger_*_*',
        'keep *_rootTupleTriggerObjects_*_*',
        # GEN objects
        'keep *_rootTupleGenEventInfo_*_*',
        'keep *_rootTupleGenParticles_*_*',
        'keep *_rootTupleGenJets*_*_*',
        'keep *_rootTupleGenElectrons*_*_*',
        'keep *_rootTupleGenMuons*_*_*',
        'keep *_rootTupleGenTaus*_*_*',
        'keep *_rootTupleGenMETTrue_*_*',
    )
)

#----------------------------------------------------------------------------------------------------
# Define the path 
#----------------------------------------------------------------------------------------------------
# to see EventSetup content
#process.esContent = cms.EDAnalyzer("PrintEventSetupContent")
# to see Event content
#process.load('FWCore.Modules.printContent_cfi')
# tracer
#process.Tracer = cms.Service('Tracer')

process.p = cms.Path(
    # L+J Filter
    process.LJFilter*  
    #process.HBHENoiseFilterResultProducer* #produces HBHE baseline bools
    # Put everything into the tree
    # In unscheduled mode, anything 'kept' in the output commands above
    #  will have its producer module called automatically
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

#----------------------------------------------------------------------------------------------------
# Run the path
#----------------------------------------------------------------------------------------------------
# Delete predefined Endpath (needed for running with CRAB)
del process.out
del process.outpath

process.schedule = cms.Schedule(process.p)#,process.DUMP)

#----------------------------------------------------------------------------------------------------
# Set global settings (number of events, global tag, input files, etc)
#----------------------------------------------------------------------------------------------------
# Make sure a correct global tag is used:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Valid_Global_Tags_by_Release
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag = GlobalTag(process.GlobalTag, 'GR_P_V56', '')
# just plain GlobalTag
# MC
process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v4'
# Data
#process.GlobalTag.globaltag = '74X_dataRun2_v2'
#process.GlobalTag.globaltag = '74X_dataRun2_reMiniAOD_v0'
# feed it into the ntuple
process.rootTupleEvent.globalTag = process.GlobalTag.globaltag

# Events to process
process.maxEvents.input = 500

# Input files
process.source.fileNames = [
    # specified by InputList.txt
    # Here is a DYJetsToLL file
    #'/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/D461E72B-306D-E511-9457-90B11C04FE38.root'
    # Data files
    #'/store/data/Run2015C/SingleMuon/MINIAOD/PromptReco-v1/000/254/833/00000/220E01C3-104B-E511-837F-02163E015541.root'
    #'root://cms-xrd-global.cern.ch//store/data/Run2015D/SingleMuon/MINIAOD/PromptReco-v3/000/256/630/00000/BCD78EF7-2B5F-E511-A3A3-02163E0170B5.root'
    #'root://cms-xrd-global.cern.ch//store/data/Run2015D/SinglePhoton/MINIAOD/05Oct2015-v1/10000/006B6A67-B26F-E511-8341-002590593902.root'
    # reHLT
    #'/store/group/phys_exotica/leptonsPlusJets/leptoquarks/13TeVSamples/ReHLT/7415/LQToUE_M-300_BetaOne_TuneCUETP8M1_13TeV-pythia8/LQ300_BetaOne_7415ReHLT/151124_100054/0000/outputPhysicsEGammaCommissioning_1.root'
    #'/store/mc/RunIISpring15MiniAODv2/TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/0A69ADC9-CA6D-E511-9A59-34E6D7E05F28.root'
    # LQM300
    '/store/mc/RunIISpring15MiniAODv2/LQToUE_M-300_BetaOne_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/5E4FDA56-0E72-E511-88C0-000F530E4784.root'
    ]


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
process.cleanMuonTriggerMatchHLTSingleMuon.matched = 'unpackedPatTrigger'
process.cleanMuonTriggerMatchHLTSingleIsoMuon.matched = 'unpackedPatTrigger'
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
# in any case, below does not work in CMSSW_7_4_14 with 2015D PromptReco-v3
## Reproduce "raw" MET from packedPFCandidates
#from RecoMET.METProducers.PFMET_cfi import pfMet
#process.pfMet = pfMet.clone(src = "packedPFCandidates")
#process.pfMet.calculateSignificance = False # this can't be easily implemented on packed PF candidates at the moment (as of Feb 17 2015)

#----------------------------------------------------------------------------------------------------
# MET filters
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
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
#Flag_trkPOG_manystripclus53X = cms.Path(manystripclus53X)
#Flag_trkPOG_toomanystripclus53X = cms.Path(toomanystripclus53X)
#Flag_trkPOG_logErrorTooManyClusters = cms.Path(logErrorTooManyClusters)
# See: https://github.com/cms-sw/cmssw/blob/CMSSW_7_2_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
# They are stored in edm::TriggerResults of the PAT process, and can be checked the same way as HLT paths.
# As per https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#early_data_issues_for_the_septem
#  For jamboree, use:
#  - Flag_HBHENoiseFilter: needs to be re-run manually as below
#  - Flag_HBHENoiseIsoFilter: needs to be re-run manually as below
#  - Flag_CSCTightHaloFilter: will need to be replaced with run2 version (use txt files)
#  - Flag_goodVertices 
#  - Flag_eeBadScFilter (use txt files for one bad SC; 3 others are tagged by this flag)
# HCAL_Noise_Filter
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")
#process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
#   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
#   reverseDecision = cms.bool(False)
#)
#process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
#   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
#   reverseDecision = cms.bool(False)
#)


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

# need the egamma smearing for 74X: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMSmearer
process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
correctionType = "Prompt2015"
process.calibratedPatElectrons.isMC = varOptions.isMC
# ntuplize this corrected electron collection
process.rootTupleElectrons.InputTag = cms.InputTag('calibratedPatElectrons','')

# HEEP ID 5.1/6.0
# Also load Egamma cut-based ID in new VID framework while we're at it
# See: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
# Set up everything that is needed to compute electron IDs and
# add the ValueMaps with ID decisions into the event data stream
#
# Load tools and function definitions
#from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
#process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
#process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')
#from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
#process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
# Define which IDs we want to produce
# Each of these two example IDs contains all four standard
# cut-based ID working points
my_id_modules = []
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff')
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff')
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff') # for 50 ns, 13 TeV data
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff')
#Add them to the VID producer
for idmod in my_id_modules:
  setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
# XXX NB, must be the same as input collection used for electron ntuplizer
process.egmGsfElectronIDs.physicsObjectSrc = process.rootTupleElectrons.InputTag
process.electronMVAValueMapProducer.srcMiniAOD = process.rootTupleElectrons.InputTag


#----------------------------------------------------------------------------------------------------
# Add MVA electron ID
#
# MVA electron ID details on this twiki:
# https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
#
#----------------------------------------------------------------------------------------------------
#TODO

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
# Rebuild AK5 jets
process.load('Leptoquarks.RootTupleMakerV2.ak5pfjets_cfi')
## b-tag discriminators
#bTagDiscriminators = [
#  'pfTrackCountingHighEffBJetTags',
#  'pfTrackCountingHighPurBJetTags',
#  'pfJetProbabilityBJetTags',
#  'pfJetBProbabilityBJetTags',
#  'pfSimpleSecondaryVertexHighEffBJetTags',
#  'pfSimpleSecondaryVertexHighPurBJetTags',
#  'pfCombinedSecondaryVertexV2BJetTags',
#  'pfCombinedInclusiveSecondaryVertexV2BJetTags',
#  'pfCombinedMVABJetTags'
#  ]
bTagDiscriminatorsAK5 = [
    #'pfCombinedInclusiveSecondaryVertexV2BJetTags'
    'None' # SIC: couldn't get AK5 jets to work otherwise
]
from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,
                 labelName = 'AK5PF',
                 jetSource = cms.InputTag('ak5PFJets'),
                 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
                 pfCandidates = cms.InputTag('packedPFCandidates'),
                 svSource = cms.InputTag('slimmedSecondaryVertices'),
                 btagDiscriminators = bTagDiscriminatorsAK5,
                 jetCorrections = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
                 genJetCollection = cms.InputTag('ak5GenJetsNoNu'),
                 genParticles = cms.InputTag('prunedGenParticles'),
                 algo = 'AK',
                 rParam = 0.5,
)
#addJetCollection(process,
#                 labelName = 'AK5PFCHS',
#                 jetSource = cms.InputTag('ak5PFJetsCHS'),
#                 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
#                 pfCandidates = cms.InputTag('packedPFCandidates'),
#                 svSource = cms.InputTag('slimmedSecondaryVertices'),
#                 btagDiscriminators = bTagDiscriminatorsAK5,
#                 jetCorrections = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
#                 genJetCollection = cms.InputTag('ak5GenJetsNoNu'),
#                 genParticles = cms.InputTag('prunedGenParticles'),
#                 algo = 'AK',
#                 rParam = 0.5,
#)
#process.patJetsAK5.userData.userFloats.src = [] # start with empty list of user floats
#process.selectedPatJetsAK5PFCHS.cut = cms.string("pt > 10") # should be 15, but we can keep RunI def.
#process.patJetGenJetMatchAK5.matched =  'slimmedGenJets'
#process.patJetPartonMatchAK5.matched = 'prunedGenParticles'
#process.patJetPartons.particles = 'prunedGenParticles'
#process.patJetPartonsLegacy.src = 'prunedGenParticles'
#process.patJetCorrFactorsAK5.primaryVertices = 'offlineSlimmedPrimaryVertices'
## needed when using btagInfos?
##process.jetTracksAssociatorAtVertexAK5.tracks = 'unpackedTracksAndVertices'
##process.jetTracksAssociatorAtVertexAK5.pvSrc = 'offlineSlimmedPrimaryVertices'
from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))

#process.out.outputCommands += ['keep *_ak5GenJets_*_EX',
#                               'keep *_ak5PFJets_*_EX',
#                               'keep *_ak5PFJetsCHS_*_EX', ]


#----------------------------------------------------------------------------------------------------
# JER and JEC
#----------------------------------------------------------------------------------------------------
jerResFile = 'Summer15_25nsV6_MC_PtResolution_AK4PFchs.txt'
jerScaleFactorsFile = 'Summer15_25nsV6_DATAMCSF_AK4PFchs.txt'
# JEC from text files
# need these text files in the same directory as this cfg.py
# expand the CMSSW_SEARCH_PATH
jecUncFileData='Leptoquarks/RootTupleMakerV2/test/'+'Summer15_25nsV6_DATA_UncertaintySources_AK4PFchs.txt'
jecUncFileMC='Leptoquarks/RootTupleMakerV2/test/'+'Summer15_25nsV6_MC_UncertaintySources_AK4PFchs.txt'

dbJetMCDBFile = 'Summer15_25nsV6_MC.db'
dbJetDataDBFile = 'Summer15_25nsV6_DATA.db'

# Get JER from text files
process.rootTuplePFJetsAK5.ReadJERFromGT = False
#process.rootTuplePFJetsAK5CHS.ReadJERFromGT = False
process.rootTuplePFJetsAK4CHS.ReadJERFromGT = False
process.rootTuplePFJetsAK4Puppi.ReadJERFromGT = False
process.rootTuplePFJetsAK5.JERResolutionsFile = jerResFile
#process.rootTuplePFJetsAK5CHS.JERResolutionsFile = jerResFile
process.rootTuplePFJetsAK4CHS.JERResolutionsFile = jerResFile
process.rootTuplePFJetsAK4Puppi.JERResolutionsFile = jerResFile
process.rootTuplePFJetsAK5.JERScaleFactorsFile = jerScaleFactorsFile
#process.rootTuplePFJetsAK5CHS.JERScaleFactorsFile = jerScaleFactorsFile
process.rootTuplePFJetsAK4CHS.JERScaleFactorsFile = jerScaleFactorsFile
process.rootTuplePFJetsAK4Puppi.JERScaleFactorsFile = jerScaleFactorsFile
# XXX NB: These files are wrong for AK4Puppi and AK5, so they won't have proper scale factors.

## load JES/etc.from local db file, and into global tag
#jetDBFile = 'sqlite:'+dbJetMCDBFile if varOptions.isMC else 'sqlite:'+dbJetDataDBFile
#process.load("CondCore.DBCommon.CondDBCommon_cfi")
#from CondCore.DBCommon.CondDBSetup_cfi import *
## JEC
#jecPSetDataAK4chs = cms.PSet(
#    record = cms.string('JetCorrectionsRecord'),
#    tag    = cms.string('JetCorrectorParametersCollection_Summer15_25nsV6_DATA_AK4PFchs'),
#    label  = cms.untracked.string('AK4PFchs')
#)
#jecPSetMCAK4chs = cms.PSet(
#    record = cms.string('JetCorrectionsRecord'),
#    tag    = cms.string('JetCorrectorParametersCollection_Summer15_25nsV6_MC_AK4PFchs'),
#    label  = cms.untracked.string('AK4PFchs')
#)
#process.jec = cms.ESSource("PoolDBESSource",
#      DBParameters = cms.PSet(
#        messageLevel = cms.untracked.int32(0)
#        ),
#      timetype = cms.string('runnumber'),
#      toGet = cms.VPSet(
#        jecPSetMCAK4chs if varOptions.isMC else jecPSetDataAK4chs
#      ), 
#       connect = cms.string(jetDBFile)
#     # connect = cms.string('sqlite:Summer12_V7_DATA.db')
#     # uncomment above tag lines and this comment to use MC JEC
#     # connect = cms.string('sqlite:Summer12_V7_MC.db')
#)
### add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
#process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')
## JER from local db
#process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
#from CondCore.DBCommon.CondDBSetup_cfi import *
#process.jer = cms.ESSource("PoolDBESSource",
#        CondDBSetup,
#        toGet = cms.VPSet(
#            # Resolution
#            cms.PSet(
#                record = cms.string('JetResolutionRcd'),
#                tag    = cms.string('JER_MC_PtResolution_Summer15_25nsV6_AK4PFchs'),
#                label  = cms.untracked.string('AK4PFchs')
#                ),
#
#            # Scale factors
#            cms.PSet(
#                record = cms.string('JetResolutionScaleFactorRcd'),
#                tag    = cms.string('JER_DATAMCSF_Summer15_25nsV6_AK4PFchs'),
#                label  = cms.untracked.string('AK4PFchs')
#                ),
#        ),
#        connect = cms.string(jetDBFile)
#)
#process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')

# Apply jet energy corrections:
#   https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrPatJets
# This should load them from the global tag
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
  src = cms.InputTag("slimmedJets"),
  levels = ['L1FastJet', 
        'L2Relative', 
        'L3Absolute'],
  payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
process.patJetsReapplyJEC = patJetsUpdated.clone(
  jetSource = cms.InputTag("slimmedJets"),
  jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
)

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
# MET Uncertainties
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#MET_Systematics_Tools
#----------------------------------------------------------------------------------------------------
postfix = 'RecorrectedMet'
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
myJecUncFile = jecUncFileMC if varOptions.isMC else jecUncFileData
#default configuration for miniAOD reprocessing
#for a full met computation, remove the pfCandColl input
runMetCorAndUncFromMiniAOD(process,
                           jetCollUnskimmed='slimmedJets',
                           isData=not varOptions.isMC,
                           electronColl=cms.InputTag('calibratedPatElectrons'),
                           repro74X=True,
                           jecUncFile=myJecUncFile,
                           postfix=postfix,
)
# this should be it, but fails
#process.applyCorrections = cms.Path(getattr(process,'fullPatMetSequence{0}'.format(postfix)))
# fix things
getattr(process,'patPFMetT1T2Corr{0}'.format(postfix)).src = cms.InputTag('patJets')
getattr(process,'patPFMetT2Corr{0}'.format(postfix)).src = cms.InputTag('patJets')
getattr(process,'patPFMetTxyCorr{0}'.format(postfix)).vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices')
process.applyCorrections = cms.Path()
if varOptions.isMC: process.applyCorrections += process.genMetExtractor
process.applyCorrections += getattr(process,'patPFMet{0}'.format(postfix))
process.applyCorrections += process.patJetCorrFactorsReapplyJEC
process.applyCorrections += process.patJets
process.applyCorrections += getattr(process,'patPFMetT1Txy{0}'.format(postfix))
process.applyCorrections += getattr(process,'patPFMetT1{0}'.format(postfix))
process.applyCorrections += getattr(process,'patPFMetTxy{0}'.format(postfix))
#process.schedule.append(process.applyCorrections)

# ntuplize the newly-corrected MET
process.rootNTupleNewMETs = cms.Path()
allowedShifts = ['jres','jes','mes','ees','tes','ues']
allowedSigns = ['+','-']
signMap = {
    '+' : 'Up',
    '-' : 'Down',
}
collMap = {
    #'jres' : {'Jets'     : 'shiftedPatJetRes{sign}{postfix}'},
    'jres' : {},
    'jes'  : {'Jets'     : 'shiftedPatJetEn{sign}{postfix}'},
    'mes'  : {'Muons'    : 'shiftedPatMuonEn{sign}{postfix}'},
    'ees'  : {'Electrons': 'shiftedPatElectronEn{sign}{postfix}'},
    'tes'  : {'Taus'     : 'shiftedPatTauEn{sign}{postfix}'},
    'ues'  : {},
    'pes'  : {},
    }
metMap = {
  'jres' : 'patPFMetT1JetRes{sign}{postfix}',
  'jes'  : 'patPFMetT1JetEn{sign}{postfix}',
  'mes'  : 'patPFMetT1MuonEn{sign}{postfix}',
  'ees'  : 'patPFMetT1ElectronEn{sign}{postfix}',
  'tes'  : 'patPFMetT1TauEn{sign}{postfix}',
  'ues'  : 'patPFMetT1UnclusteredEn{sign}{postfix}',
  'pes'  : '',
}
#mettypes = ['CaloMET','CaloMETType1Cor','PFMET','PFMETType1Cor','PFMETType01Cor','PFMETType01XYCor','PFMETPuppi','PFMETPuppiType1Cor']
mettypes = ['PFMETType1Cor']
for shift in allowedShifts:
    for sign in allowedSigns:
      # FIXME TODO
        ## embed shifted objects
        #for coll in collMap[shift]:
        #    modName = '{shift}{sign}{coll}Embedding'.format(shift=shift,sign=signMap[sign],coll=coll)
        #    pluginName = 'MiniAODShifted{coll}Embedder'.format(coll=coll[:-1])
        #    dName = coll.lower()
        #    srcName = fs_daughter_inputs[dName]
        #    shiftSrcName = collMap[shift][coll].format(sign=signMap[sign],postfix=postfix)
        #    label = '{shift}{sign}{coll}'.format(shift=shift,sign=signMap[sign],coll=coll)
        #    module = cms.EDProducer(
        #        pluginName,
        #        src = cms.InputTag(srcName),
        #        shiftSrc = cms.InputTag(shiftSrcName),
        #        label = cms.string(label),
        #    )
        #    setattr(process,modName,module)
        #    fs_daughter_inputs[dName] = modName
        #    #process.embedShifts *= getattr(process,shiftSrcName)
        #    process.embedShifts *= getattr(process,modName)
        # embed shifted met
        for mettype in mettypes:
          modName = 'rootTuple{mettype}{shift}{sign}'.format(mettype=mettype,shift=shift,sign=signMap[sign])
          metName = metMap[shift].format(sign=signMap[sign],postfix=postfix)
          #label = '{shift}{sign}MET'.format(shift=shift,sign=signMap[sign])
          prefix = '{mettype}{shift}{sign}'.format(mettype=mettype,shift=shift,sign=signMap[sign])
          #prefix = modName
          #suffix = mettype[7:] if 'Calo' in mettype else mettype[5:]
          #suffix = postfix
          suffix = ''
          module = cms.EDProducer(
              'RootTupleMakerV2_MET',
              InputTag = cms.InputTag(metName),
              Prefix = cms.string(prefix),
              Suffix = cms.string(suffix),
              StoreUncorrectedMET = cms.bool(False), # this won't work either
              StoreMETSignificance = cms.bool(False),
              Uncertainty = cms.string('NoShift'),
              CorrectionLevel = cms.string('NoCorrection'), # call pt(), etc., instead of shiftedPt()
          )
          setattr(process,modName,module)
          #fs_daughter_inputs['pfmet'] = modName
          process.rootNTupleNewMETs *= getattr(process,modName)
#process.schedule.append(process.rootNTupleNewMETs)


#----------------------------------------------------------------------------------------------------
# Available pat::MET collections for analysis
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#ETmiss
#----------------------------------------------------------------------------------------------------

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
# PDF weights
#----------------------------------------------------------------------------------------------------

#process.pdfWeights = cms.EDProducer("PdfWeightProducer",
#	# Fix POWHEG if buggy (this PDF set will also appear on output,
#	# so only two more PDF sets can be added in PdfSetNames if not "")
#	#FixPOWHEG = cms.untracked.string("CT10.LHgrid"),
#       # GenTag = cms.untracked.InputTag("genParticles"),
#	PdfInfoTag = cms.untracked.InputTag("generator"),
#	PdfSetNames = cms.untracked.vstring(
#            "CT10.LHgrid" , 
#            "MSTW2008nlo68cl.LHgrid",
#            #"NNPDF20_100.LHgrid"#DMM FIXME couldn't find this file in /cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_4/src/PhysicsTools/PatUtils/data/Summer13_V1_DATA_UncertaintySources_AK5PF.txt
#	)
#)

#----------------------------------------------------------------------------------------------------
# Define GEN particle skimmer modules
#----------------------------------------------------------------------------------------------------

process.load ('Leptoquarks.LeptonJetGenTools.genTauMuElFromZs_cfi')
process.load ('Leptoquarks.LeptonJetGenTools.genTauMuElFromWs_cfi') 




#print process.dumpPython()

