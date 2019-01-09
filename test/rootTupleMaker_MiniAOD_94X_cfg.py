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
varOptions.register(
    "era",
    "2016",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Customize config for era"
)
varOptions.parseArguments()
if varOptions.isMC:
    print 'We are running on MC.'
else:
    print 'We are running on Data.'
if varOptions.era=='2017' or varOptions.era=='2016':
    print 'We are running in the',varOptions.era,'era.'
else:
    print 'ERROR: did not understand era option given as "'+varOptions.era+'"; valid values are 2016 or 2017.'
    exit(-1)

process.load('Configuration.StandardSequences.Services_cff')
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    calibratedPatPhotons = cms.PSet(
      initialSeed = cms.untracked.uint32(81),
      engineName = cms.untracked.string('TRandom3')
      ),
    calibratedPatElectrons = cms.PSet(
      initialSeed = cms.untracked.uint32(81),
      engineName = cms.untracked.string('TRandom3')
      ),
)
process.load('JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff')

# Change process name
process._Process__name="ROOTTUPLEMAKERV2"
# Options and Output Report
process.options.wantSummary = False
#Print the information about the execution of the modules, and the transitions (beginRun, etc.) they go through. For unscheduled mode debugging.
#process.Tracer = cms.Service("Tracer")

############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.default.limit = 10
process.MessageLogger.threshold = cms.untracked.string('DEBUG')
process.MessageLogger.categories.extend(["GetManyWithoutRegistration","GetByLabelWithoutRegistration"])
_messageSettings = cms.untracked.PSet(
        reportEvery = cms.untracked.int32(1),
        optionalPSet = cms.untracked.bool(True),
        limit = cms.untracked.int32(10000000)
        )

process.MessageLogger.cerr.GetManyWithoutRegistration = _messageSettings
process.MessageLogger.cerr.GetByLabelWithoutRegistration = _messageSettings
#################################################################

#----------------------------------------------------------------------------------------------------
# Load our RootTupleMakerV2 modules
#----------------------------------------------------------------------------------------------------

process.load('Leptoquarks.RootTupleMakerV2.Ntuple_cff')

#----------------------------------------------------------------------------------------------------
# Output ROOT file
#----------------------------------------------------------------------------------------------------

process.TFileService = cms.Service("TFileService",
    fileName = cms.string( "test_file_m1550.root" )
    #fileName = cms.string( "file_ntuple.root" )
)

#----------------------------------------------------------------------------------------------------
# Set global settings (number of events, global tag, input files, etc)
#----------------------------------------------------------------------------------------------------
# Make sure a correct global tag is used:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Valid_Global_Tags_by_Release
# https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag = GlobalTag(process.GlobalTag, 'GR_P_V56', '')
# just plain GlobalTag
# MC
if varOptions.isMC:
    if varOptions.era=='2017':
        process.GlobalTag.globaltag = '94X_mc2017_realistic_v17' # Fall17 MC
    elif varOptions.era=='2016':
        process.GlobalTag.globaltag = '94X_mcRun2_asymptotic_v3' # Summer16 MC
else:
    if varOptions.era=='2017':
        process.GlobalTag.globaltag = '94X_dataRun2_v11' # 17Nov2017
    elif varOptions.era=='2016':
        process.GlobalTag.globaltag = '94X_dataRun2_v10' # 07Aug2017 Rereco
# feed it into the ntuple
process.rootTupleEvent.globalTag = process.GlobalTag.globaltag


# Events to process
process.maxEvents.input = -1

# Input files
process.source.fileNames = [
    # specified by InputList.txt
    '/store/mc/RunIISummer16MiniAODv2/LQToUE_M-1550_BetaOne_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/A061D453-57C8-E611-B57E-008CFA11136C.root',
    '/store/mc/RunIISummer16MiniAODv2/LQToUE_M-1550_BetaOne_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/562CA91F-4FC8-E611-B5ED-008CFA110B08.root'
]

#----------------------------------------------------------------------------------------------------
# Turn on trigger matching
#----------------------------------------------------------------------------------------------------
# stored in pat::TriggerObjectStandAlone with default MiniAOD config
# See: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
# To access, see: https://hypernews.cern.ch/HyperNews/CMS/get/physTools/3286/1/2.html
# This is based off of the hypernews post
process.load('PhysicsTools.PatAlgos.slimming.unpackedPatTrigger_cfi')
process.unpackedPatTrigger.unpackFilterLabels = cms.bool(True)
process.unpackedPatTrigger.patTriggerObjectsStandAlone = cms.InputTag('selectedPatTrigger')
process.rootTupleEventSelection.L1InputTag  = 'gtStage2Digis'
####
#####fixme removing because no trigger info?
####process.cleanElectronTriggerMatchHLTSingleElectron.matched = 'unpackedPatTrigger'
####process.cleanElectronTriggerMatchHLTSingleElectronWP85.matched = 'unpackedPatTrigger'
####process.cleanElectronTriggerMatchHLTDoubleElectron.matched = 'unpackedPatTrigger'
####process.cleanElectronTriggerMatchHLTElectronJetJet.matched = 'unpackedPatTrigger'
process.cleanMuonTriggerMatchHLTMuon.matched = 'unpackedPatTrigger'
process.cleanMuonTriggerMatchHLTSingleMuon.matched = 'unpackedPatTrigger'
process.cleanMuonTriggerMatchHLTSingleIsoMuon.matched = 'unpackedPatTrigger'
####process.pfjetTriggerMatchHLTEleJetJet.matched = 'unpackedPatTrigger'
####
####process.out.outputCommands += ['keep *_electronsTriggerMatchHLTSingleElectron_*_*',
####                               'keep *_electronsTriggerMatchHLTSingleElectronWP85 _*_*',
####                               'keep *_electronsTriggerMatchHLTDoubleElectron_*_*',
####                               'keep *_electronsTriggerMatchHLTElectroJetJetn_*_*'      ]
##### select the matching electrons and keep them
####process.out.outputCommands += ['keep *_ electronsTriggeredHLTSingleElectron_*_*',
####                               'keep *_ electronsTriggeredHLTSingleElectronWP80_*_*',
####                               'keep *_ electronsTriggeredHLTDoubleElectron_*_*' ]

# Run2 processing for egamma
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaPostRecoRecipes
process.schedule = cms.Schedule()
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
if varOptions.era=='2017': # 2017 MiniAOD V2
    era = '2017-Nov17ReReco'
elif varOptions.era=='2016':
    era = '2016-Legacy'
setupEgammaPostRecoSeq(process,
                       runVID=True,
                       era=era)
process.egPath = cms.Path(process.egammaPostRecoSeq)
process.schedule.append(process.egPath)
process.load('Leptoquarks.RootTupleMakerV2.heepV70Modifier_cfi')
process.slimmedElectrons.modifierConfig.modifications.extend(process.heep_modifications)
# For 2016 MiniAOD v2, no additional corrections needed

process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")

#----------------------------------------------------------------------------------------------------
# JER and JEC
#----------------------------------------------------------------------------------------------------
## JER from text files
#jerResFile               = 'Leptoquarks/RootTupleMakerV2/data/Spring16_25nsV6_MC_PtResolution_AK4PFchs.txt'
#jerResFilePuppi          = 'Leptoquarks/RootTupleMakerV2/data/Spring16_25nsV6_MC_PtResolution_AK4PFPuppi.txt'
#jerScaleFactorsFile      = 'Leptoquarks/RootTupleMakerV2/data/Spring16_25nsV6_MC_SF_AK4PFchs.txt'
#jerScaleFactorsFilePuppi = 'Leptoquarks/RootTupleMakerV2/data/Spring16_25nsV6_MC_SF_AK4PFPuppi.txt'
## JEC from text files
## stored in 'data' directory, they should be available in the CMSSW_SEARCH_PATH
#jecUncFileData           ='Leptoquarks/RootTupleMakerV2/data/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt'
#jecUncFileDataPuppi      ='Leptoquarks/RootTupleMakerV2/data/Spring16_25nsV6_DATA_Uncertainty_AK4PFPuppi.txt'
#jecUncFileMC             ='Leptoquarks/RootTupleMakerV2/data/Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt'
#jecUncFileMCPuppi        ='Leptoquarks/RootTupleMakerV2/data/Spring16_25nsV6_MC_Uncertainty_AK4PFPuppi.txt'
# Load JER from text files into ntuplizer
#process.rootTuplePFJetsAK4CHS.JERResolutionsFile    = jerResFile
#process.rootTuplePFJetsAK4Puppi.JERResolutionsFile  = jerResFilePuppi
#process.rootTuplePFJetsAK4CHS.JERScaleFactorsFile   = jerScaleFactorsFile
#process.rootTuplePFJetsAK4Puppi.JERScaleFactorsFile = jerScaleFactorsFilePuppi

## Loading  JEC/JER from local DB file
#dbJetMCDBFile = 'Summer16_23Sep2016V4_MC.db'
#dbJetDataDBFile = 'Summer16_23Sep2016AllV4_DATA.db'
#
## Use JER from GR
#process.rootTuplePFJetsAK4CHS.ReadJERFromGT = True
#process.rootTuplePFJetsAK4Puppi.ReadJERFromGT = True
#
## load JES/etc.from local db file, and into global tag
#jecDBPrefix='sqlite_file:'
#jetDBFile = jecDBPrefix+dbJetMCDBFile if varOptions.isMC else jecDBPrefix+dbJetDataDBFile
#process.load("CondCore.DBCommon.CondDBCommon_cfi")
#from CondCore.DBCommon.CondDBSetup_cfi import *
## JEC
#jecPSetDataAK4chs = cms.PSet(
#    record = cms.string('JetCorrectionsRecord'),
#    tag    = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016AllV4_DATA_AK4PFchs'),
#    label  = cms.untracked.string('AK4PFchs')
#)
#jecPSetMCAK4chs = cms.PSet(
#    record = cms.string('JetCorrectionsRecord'),
#    tag    = cms.string('JetCorrectorParametersCollection_Summer16_23Sep2016V4_MC_AK4PFchs'),
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
#)
### add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
#process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')
#
## JER from local db
#process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
#from CondCore.DBCommon.CondDBSetup_cfi import *
#process.jer = cms.ESSource("PoolDBESSource",
#        CondDBSetup,
#        toGet = cms.VPSet(
#            # Resolution
#            cms.PSet(
#                record = cms.string('JetResolutionRcd'),
#                tag    = cms.string('JER_MC_PtResolution_Spring16_25nsV6_AK4PFchs'),
#                label  = cms.untracked.string('AK4PFchs')
#                ),
#
#            # Scale factors
#            cms.PSet(
#                record = cms.string('JetResolutionScaleFactorRcd'),
#                tag    = cms.string('JER_DATAMCSF_Spring16_25nsV6_AK4PFchs'),
#                label  = cms.untracked.string('AK4PFchs')
#                ),
#        ),
#        connect = cms.string(jetDBFile)
#)
#process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')
 
# Apply jet energy corrections:
#   https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrPatJets
#   https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/PhysicsTools/PatAlgos/test/patTuple_updateJets_fromMiniAOD_cfg.py
# This should load them from the global tag
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

if varOptions.isMC : 
  updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
  )
else : 
  updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
  )
process.rootTuplePFJetsAK4CHS.InputTag = cms.InputTag('updatedPatJets')
process.applyCorrections = cms.Path()
process.applyCorrections += process.patJetCorrFactors
process.applyCorrections += process.updatedPatJets
process.schedule.append(process.applyCorrections)
if varOptions.isMC : 
    process.rootTuplePFJetsAK4CHS.JECUncertaintySourcesFile = cms.FileInPath('Leptoquarks/RootTupleMakerV2/data/Summer15_25nsV7_MC_UncertaintySources_AK4PFchs.txt')
    process.rootTuplePFJetsAK4Puppi.JECUncertaintySourcesFile = cms.FileInPath('Leptoquarks/RootTupleMakerV2/data/Summer15_25nsV7_MC_UncertaintySources_AK4PFchs.txt')
else:
    process.rootTuplePFJetsAK4CHS.JECUncertaintySourcesFile = cms.FileInPath('Leptoquarks/RootTupleMakerV2/data/Summer15_25nsV7_DATA_UncertaintySources_AK4PFchs.txt')
    process.rootTuplePFJetsAK4Puppi.JECUncertaintySourcesFile = cms.FileInPath('Leptoquarks/RootTupleMakerV2/data/Summer15_25nsV7_DATA_UncertaintySources_AK4PFchs.txt')
process.rootTuplePFJetsAK4CHS.ReadJECuncertainty = cms.bool(False)
process.rootTuplePFJetsAK4Puppi.ReadJECuncertainty = cms.bool(False)

##----------------------------------------------------------------------------------------------------
## MET Re-Corrections and Uncertainties
## https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
## https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#METSysTools
##----------------------------------------------------------------------------------------------------
# XXX FIXME 94X TODO SIC

#----------------------------------------------------------------------------------------------------
# This is MC, so analyze the smeared PFJets by default
# But smearing should be done by default in MiniAOD ?
#----------------------------------------------------------------------------------------------------
# Set Lepton-Gen Matching Parameters
#----------------------------------------------------------------------------------------------------
# XXX FIXME SIC: Needed? At least some matching is done by default in PAT/MiniAOD
process.load("Leptoquarks.RootTupleMakerV2.leptonGenMatching_cfi")

#----------------------------------------------------------------------------------------------------
# Lepton + Jets filter
#----------------------------------------------------------------------------------------------------
process.load("Leptoquarks.LeptonJetFilter.leptonjetfilter_cfi")

#### Shared Muon/Electron/Tau Skim
process.LJFilter.tauLabel  = cms.InputTag("slimmedTaus")                        
process.LJFilter.muLabel   = cms.InputTag("slimmedMuons")
process.LJFilter.elecLabel = cms.InputTag("slimmedElectrons")
process.LJFilter.jetLabel  = cms.InputTag("slimmedJets")
process.LJFilter
process.LJFilter.muPT     = 8.0
process.LJFilter.electronsMin = 0
process.LJFilter.elecPT       = 11.0
process.LJFilter.tausMin = 0
process.LJFilter.tauPT   = 12.0
process.LJFilter.jetsMin = 0
process.LJFilter.jetPT   = 17.0
process.LJFilter.counteitherleptontype = True
process.LJFilter.customfilterEMuTauJet2012 = False
process.LJFilter.customfilterEMuTauJet2016 = True
process.LJFilter.debug = False
# -- WARNING :
# "customfilterEMuTauJet2012" and "customfilterEMuTauJet2016" configurations are hard-coded.
# (see: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Leptoquarks/LeptonJetFilter/src/LeptonJetFilter.cc )
# "customfilterEMuTauJet2012" is the desired mode of operation for the Lepton+Jets Filter in 2012.
# "customfilterEMuTauJet2016" is the desired mode of operation for the Lepton+Jets Filter in 2016.

#----------------------------------------------------------------------------------------------------
# PDF weights
#----------------------------------------------------------------------------------------------------

process.pdfWeights = cms.EDProducer("PdfWeightProducer",
    # Fix POWHEG if buggy (this PDF set will also appear on output,
    # so only two more PDF sets can be added in PdfSetNames if not "")
    #FixPOWHEG = cms.untracked.string("CT10.LHgrid"),
  GenTag = cms.untracked.InputTag("prunedGenParticles"),
    PdfInfoTag = cms.untracked.InputTag("generator"),
    #FIXME
    PdfSetNames = cms.untracked.vstring(
            #"PDF4LHC15_nlo_100.LHgrid" ,
            #"CT10nlo.LHgrid" , 
            #"MMHT2014nlo68cl.LHgrid",
            "NNPDF30_nlo_as_0118.LHgrid"
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
        'keep *_rootTuplePFJets*_*_*',
        'keep *_rootTupleElectrons_*_*',
        'keep *_rootTupleMuons_*_*',
        # FIXME ignore for now
        #'keep *_rootTuplePhotons_*_*',
        'keep *_rootTupleVertex_*_*',
        ## MET objects for analysis
        'keep *_rootTuplePFMET*_*_*',
        'keep *_slimmedMETs_*_Recorrected',
        # Trigger stuff
        'keep *_rootTupleTrigger_*_*',
        'keep *_rootTupleTriggerObjects_*_*',
        # GEN objects
        'keep *_rootTupleGenEventInfo_*_*',
        'keep *_rootTupleGenParticles_*_*',
        'keep *_rootTupleGenJets*_*_*',
        'keep *_rootTupleGenElectrons*_*_*',
        'keep *_rootTupleGenMuons*_*_*',
        #'keep *_rootTupleGenTaus*_*_*',
        'keep *_rootTupleGenMETTrue_*_*',
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
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)
process.p = cms.Path(
    #process.printContent*
    # supporting producers
    process.unpackedPatTrigger*
    ## L+J Filter
    process.LJFilter*  
    #FIXME process.pdfWeights*
    ## Put everything into the tree
    process.rootTuplePFCandidates*
    process.rootTupleEvent*
    process.rootTupleElectronsSequence*
    process.rootTupleMuonsSequence*
    process.rootTupleVertex*
    process.rootTuplePFJetsSequence*
    process.rootTuplePFMETType1Cor*
    #process.pfMEtMultShiftCorr2*
    process.rootTuplePFMETSequence*
    #process.rootTuplePFMETType1CorNotRecorrected*
    process.rootTupleTriggerObjects*
    process.rootTupleTrigger*
    process.rootTupleGenEventInfo*
    process.rootTupleGenParticles*
    process.rootTupleGenJetsSequence*
    process.rootTupleGenMETTrue*
    process.rootTupleEventSelection*
    process.rootTupleTree
)

process.schedule.append(process.p)

#----------------------------------------------------------------------------------------------------
# Dump if necessary
#----------------------------------------------------------------------------------------------------
#
#process.dump = cms.OutputModule("PoolOutputModule",
#                                outputCommands = cms.untracked.vstring(
#        'drop *',
#        # Event information
#        'keep *_rootTupleEvent_*_*',
#        'keep *_rootTupleEventSelection_*_*',
#        # Single objects
#        'keep *_rootTuplePFCandidates_*_*',
#        'keep *_rootTuplePFJets*_*_*',
#        'keep *_rootTupleElectrons_*_*',
#        'keep *_rootTupleMuons_*_*',
#        # FIXME ignore for now
#        #'keep *_rootTuplePhotons_*_*',
#        'keep *_rootTupleVertex_*_*',
#        ## MET objects for analysis
#        'keep *_rootTuplePFMET*_*_*',
#        # Trigger objects
#        'keep *_rootTupleTrigger_*_*',
#        'keep *_rootTupleTriggerObjects_*_*',
#        # GEN objects
#        'keep *_rootTupleGenEventInfo_*_*',
#        'keep *_rootTupleGenParticles_*_*',
#        'keep *_rootTupleGenJets*_*_*',
#        'keep *_rootTupleGenElectrons*_*_*',
#        'keep *_rootTupleGenMuons*_*_*',
#        'keep *_rootTupleGenTaus*_*_*',
#        'keep *_rootTupleGenMETTrue_*_*',
#        ),
#        )
#process.DUMP    = cms.EndPath (process.dump)
#process.schedule.append(process.DUMP)
#
#----------------------------------------------------------------------------------------------------
# Run the path
#----------------------------------------------------------------------------------------------------
# Delete predefined Endpath (needed for running with CRAB)
del process.out
del process.outpath

#print process.dumpPython()
