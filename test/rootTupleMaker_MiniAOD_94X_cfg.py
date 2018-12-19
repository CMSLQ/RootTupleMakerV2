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
else:
  print 'We are running on Data!'

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
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Valid_Global_Tags_by_Release
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag = GlobalTag(process.GlobalTag, 'GR_P_V56', '')
# just plain GlobalTag
# MC
if varOptions.isMC:
  process.GlobalTag.globaltag = '94X_mc2017_realistic_v17'
else:
  process.GlobalTag.globaltag = '94X_dataRun2_v11'
# feed it into the ntuple
process.rootTupleEvent.globalTag = process.GlobalTag.globaltag
#process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')


# Events to process
process.maxEvents.input = -1

# Input files
process.source.fileNames = [
    # specified by InputList.txt
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/120000/80DBA5F3-16BE-E811-854E-A0369FC5E530.root'
    'root://cms-xrd-global//store/mc/RunIISummer16MiniAODv3/LQToUE_M-1550_BetaOne_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/270000/5E8B6A54-47CB-E811-91C6-0025905C4264.root'
    #'file:/tmp/scooper/5E8B6A54-47CB-E811-91C6-0025905C4264.root'
]

#----------------------------------------------------------------------------------------------------
# Turn on trigger matching
#----------------------------------------------------------------------------------------------------
# stored in pat::TriggerObjectStandAlone with default MiniAOD config
# See: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
# To access, see: https://hypernews.cern.ch/HyperNews/CMS/get/physTools/3286/1/2.html
# This is based off of the hypernews post
# FIXME check that we match to the paths we want while we're updating the cpp code
#from PhysicsTools.PatAlgos.slimming.unpackedPatTrigger_cfi import unpackedPatTrigger
process.load('PhysicsTools.PatAlgos.slimming.unpackedPatTrigger_cfi')
###### need to load it into the process, or else it won't run
#process.unpackedPatTrigger = unpackedPatTrigger.clone()
process.unpackedPatTrigger.unpackFilterLabels = cms.bool(True)
process.unpackedPatTrigger.patTriggerObjectsStandAlone = cms.InputTag('slimmedPatTrigger')
#FIXME 94X SIC
#process.rootTupleEventSelection.L1InputTag  = 'gtStage2Digis'
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

process.schedule = cms.Schedule()
## EGMRegression https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMRegression
#from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
#process = regressionWeights(process)
#process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')
## Need the egamma smearing for 80X: [https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMRegression#Consistent_EGMSmearer
#process.load('EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi')
#process.calibratedPatElectrons.isMC = varOptions.isMC
#process.selectedElectrons = cms.EDFilter(
#  "PATElectronSelector",
#  src = cms.InputTag("slimmedElectrons"),
#  cut = cms.string("pt > 5 && abs(eta)<2.5")
#)
#process.calibratedPatElectrons.electrons = cms.InputTag('selectedElectrons')
## Set the lines below to True or False depending if you are correcting the scale (data) or smearing the resolution (MC) 
#process.calibratedPatElectrons.isMC = varOptions.isMC
## ntuplize this corrected electron collection
#process.rootTupleElectrons.InputTag = cms.InputTag('calibratedPatElectrons')
# HEEP 7
# Also load Egamma cut-based ID in new VID framework while we're at it
# See: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
#      https://twiki.cern.ch/twiki/bin/viewauth/CMS/HEEPElectronIdentificationRun2
#
# Load tools and function definitions
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
# Define which IDs we want to produce
# Each of these two example IDs contains all four standard cut-based ID working points
my_id_modules = []
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff')
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff')
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff')
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff')
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff')
#Add them to the VID producer
for idmod in my_id_modules:
  setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
# XXX NB, must be the same as input collection used for electron ntuplizer
process.egmGsfElectronIDs.physicsObjectSrc = process.rootTupleElectrons.InputTag
process.heepIDVarValueMaps.elesMiniAOD  = process.rootTupleElectrons.InputTag
process.electronRegressionValueMapProducer.srcMiniAOD = process.rootTupleElectrons.InputTag
process.electronMVAValueMapProducer.srcMiniAOD = process.rootTupleElectrons.InputTag#FIXME do we need this?

process.electronSupportPath = cms.Path()
#process.electronSupportPath += process.regressionApplication
#process.electronSupportPath += process.selectedElectrons
#process.electronSupportPath += process.calibratedPatElectrons
process.electronSupportPath += process.egmGsfElectronIDSequence
process.schedule.append(process.electronSupportPath)

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
#FIXME SIC 94X
#from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
#
#if varOptions.isMC : 
#  updateJetCollection(
#    process,
#    jetSource = cms.InputTag('slimmedJets'),
#    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
#  )
#else : 
#  updateJetCollection(
#    process,
#    jetSource = cms.InputTag('slimmedJets'),
#    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
#  )
#process.rootTuplePFJetsAK4CHS.InputTag = cms.InputTag('updatedPatJets')
#process.applyCorrections = cms.Path()
#process.applyCorrections += process.patJetCorrFactors
#process.applyCorrections += process.updatedPatJets
#process.schedule.append(process.applyCorrections)
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
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#myJecUncFile = jecUncFileMC if varOptions.isMC else jecUncFileData
myJecUncFile = '' # take from global tag
#default configuration for miniAOD reprocessing

#FIXME SIC 94X
#if varOptions.isMC == True:
#  from MetTools.MetPhiCorrections.tools.multPhiCorr_Summer16_MC_DY_80X_sumPt_cfi import multPhiCorr_MC_DY_sumPT_80X as multPhiCorrParams_Txy_25ns
#else :
#  from MetTools.MetPhiCorrections.tools.multPhiCorr_ReMiniAOD_Data_BCDEF_80X_sumPt_cfi import multPhiCorr_Data_BCDEF_80X as multPhiCorrParams_Txy_25ns #reMiniAOD B-F
#  #from MetTools.MetPhiCorrections.tools.multPhiCorr_ReMiniAOD_Data_GH_80X_sumPt_cfi import multPhiCorr_Data_GH_80X as multPhiCorrParams_Txy_25ns #reMiniAOD G-H
#
##Do NOT remove the following line, it is used to insert the correct XY correction from the job creator
##MetPhiCorrectionsInsertHere
#
#multPhiCorrParams_T0rtTxy_25ns     = cms.VPSet( pset for pset in multPhiCorrParams_Txy_25ns)
#multPhiCorrParams_T0rtT1Txy_25ns   = cms.VPSet( pset for pset in multPhiCorrParams_Txy_25ns)
#multPhiCorrParams_T0rtT1T2Txy_25ns = cms.VPSet( pset for pset in multPhiCorrParams_Txy_25ns)
#multPhiCorrParams_T0pcTxy_25ns     = cms.VPSet( pset for pset in multPhiCorrParams_Txy_25ns)
#multPhiCorrParams_T0pcT1Txy_25ns   = cms.VPSet( pset for pset in multPhiCorrParams_Txy_25ns)
#multPhiCorrParams_T0pcT1T2Txy_25ns = cms.VPSet( pset for pset in multPhiCorrParams_Txy_25ns)
#multPhiCorrParams_T1Txy_25ns       = cms.VPSet( pset for pset in multPhiCorrParams_Txy_25ns)
#multPhiCorrParams_T1T2Txy_25ns     = cms.VPSet( pset for pset in multPhiCorrParams_Txy_25ns)
#
#process.pfMEtMultShiftCorr2 = cms.EDProducer("MultShiftMETcorrInputProducer",
#    srcPFlow = cms.InputTag('packedPFCandidates', ''),
#    vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices'),
#    parameters = multPhiCorrParams_Txy_25ns
#)

#pfMEtSysShiftCorrSequence = cms.Sequence( pfMEtMultShiftCorr )
#process.pfMEtMultShiftCorrPath = cms.Path()
#process.pfMEtMultShiftCorrPath +=process.pfMEtMultShiftCorr
#process.schedule.append(process.pfMEtMultShiftCorrPath)
runMetCorAndUncFromMiniAOD(process,
                           #jetCollUnskimmed='slimmedJets',
                           #jetCollUnskimmed='updatedPatJets',
                           #jetColl='slimmedJets',
                           isData=not varOptions.isMC,
                           #electronColl=cms.InputTag('slimmedElectrons'),
                           #jecUncFile=myJecUncFile,
)
#FIXME SIC 94X?
#postfixMuEGClean = "MuEGClean"
#postfixMuClean = 'MuClean'
## Now we make the e/g corrected MET on top of the bad muon corrected MET (on re-miniaod)
## https://twiki.cern.ch/twiki/bin/view/CMSPublic/ReMiniAOD03Feb2017Notes#MET_Recipes
#if not varOptions.isMC:
#  from PhysicsTools.PatUtils.tools.corMETFromMuonAndEG import corMETFromMuonAndEG
#  corMETFromMuonAndEG(process,
#                   pfCandCollection="", #not needed                                                                                                                                \
#                   electronCollection="slimmedElectronsBeforeGSFix",
#                   photonCollection="slimmedPhotonsBeforeGSFix",
#                   corElectronCollection="calibratedPatElectrons",
#                   corPhotonCollection="slimmedPhotons",
#                   allMETEGCorrected=True,
#                   muCorrection=False,
#                   eGCorrection=True,
#                   runOnMiniAOD=True,
#                   postfix=postfixMuEGClean
#  )
#  #process.load('PhysicsTools.PatAlgos.slimming.slimmedMETs_cfi')
#  process.slimmedMETsMuEGClean = process.slimmedMETs.clone()
#  process.slimmedMETsMuEGClean.src = cms.InputTag("patPFMetT1MuEGClean")
#  process.slimmedMETsMuEGClean.rawVariation =  cms.InputTag("patPFMetRawMuEGClean")
#  process.slimmedMETsMuEGClean.t1Uncertainties = cms.InputTag("patPFMetT1MuEGClean")
#  del process.slimmedMETsMuEGClean.caloMET
#  process.egcorrMET = cms.Sequence(
#     process.cleanedPhotonsMuEGClean+process.cleanedCorPhotonsMuEGClean+
#     process.matchedPhotonsMuEGClean + process.matchedElectronsMuEGClean +
#     process.corMETPhotonMuEGClean+process.corMETElectronMuEGClean+
#     process.patPFMetT1MuEGClean+process.patPFMetRawMuEGClean+
#     process.patPFMetT1SmearMuEGClean+process.patPFMetT1TxyMuEGClean+
#     process.patPFMetTxyMuEGClean+process.patPFMetT1JetEnUpMuEGClean+
#     process.patPFMetT1JetResUpMuEGClean+process.patPFMetT1SmearJetResUpMuEGClean+
#     process.patPFMetT1ElectronEnUpMuEGClean+process.patPFMetT1PhotonEnUpMuEGClean+
#     process.patPFMetT1MuonEnUpMuEGClean+process.patPFMetT1TauEnUpMuEGClean+
#     process.patPFMetT1UnclusteredEnUpMuEGClean+process.patPFMetT1JetEnDownMuEGClean+
#     process.patPFMetT1JetResDownMuEGClean+process.patPFMetT1SmearJetResDownMuEGClean+
#     process.patPFMetT1ElectronEnDownMuEGClean+process.patPFMetT1PhotonEnDownMuEGClean+
#     process.patPFMetT1MuonEnDownMuEGClean+process.patPFMetT1TauEnDownMuEGClean+
#     process.patPFMetT1UnclusteredEnDownMuEGClean+process.slimmedMETsMuEGClean)
#  process.egCorrMETPath = cms.Path(process.egcorrMET)
#  # schedule the MET corrections first
#  process.metCorAndUncPath = cms.Path()
#  process.metCorAndUncPath += getattr(process,'fullPatMetSequence{0}'.format(''))
#  process.schedule.append(process.metCorAndUncPath)
#  process.schedule.append(process.egCorrMETPath)
#else: # is MC
#  # Now you are creating the bad muon corrected MET
#  process.load('RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff')
#  process.badGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)
#  process.cloneGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)
#  from PhysicsTools.PatUtils.tools.muonRecoMitigation import muonRecoMitigation
#
#  muonRecoMitigation(
#                     process = process,
#                     pfCandCollection = "packedPFCandidates", #input PF Candidate Collection
#                     runOnMiniAOD = True, #To determine if you are running on AOD or MiniAOD
#                     selection="", #You can use a custom selection for your bad muons. Leave empty if you would like to use the bad muon recipe definition.
#                     muonCollection="", #The muon collection name where your custom selection will be applied to. Leave empty if you would like to use the bad muon recipe definition.
#                     cleanCollName="cleanMuonsPFCandidates", #output pf candidate collection ame
#                     cleaningScheme="computeAllApplyClone", #Options are: "all", "computeAllApplyBad","computeAllApplyClone". Decides which (or both) bad muon collections to be used for MET cleaning coming from the bad muon recipe.
#                     #postfix=postfixMuClean #Use if you would like to add a post fix to your muon / pf collections
#                     postfix='' #Use if you would like to add a post fix to your muon / pf collections
#                     )
#
#  runMetCorAndUncFromMiniAOD(process,
#                         isData=not varOptions.isMC,
#                         pfCandColl="cleanMuonsPFCandidates",
#                         recoMetFromPFCs=True,
#                         postfix=postfixMuClean
#                         )
#
#  process.mucorMET = cms.Sequence(                     
#      process.badGlobalMuonTaggerMAOD *
#      process.cloneGlobalMuonTaggerMAOD *
#      #getattr(process,'cloneGlobalMuonTaggerMAOD{0}'.format(postfixMuClean)) *
#      #process.badMuons * # If you are using cleaning mode "all", uncomment this line
#      process.cleanMuonsPFCandidates *
#      #getattr(process,'cleanMuonsPFCandidates{0}'.format(postfixMuClean)) *
#      process.fullPatMetSequenceMuClean
#      )
#  process.mucorMETPath = cms.Path(process.mucorMET)
#  process.schedule.append(process.mucorMETPath)
#  # now schedule the MET corrections later
#  process.metCorAndUncPath = cms.Path()
#  process.metCorAndUncPath += getattr(process,'fullPatMetSequence{0}'.format(''))
#  process.schedule.append(process.metCorAndUncPath)
#
## ntuplize the newly-corrected MET
#if not varOptions.isMC:
#  postfix=postfixMuEGClean
#  #getattr(process,'patPFMetTxyMuEGClean').vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices')
#  #getattr(process,'patPFMetTxyMuEGClean').srcPFlow = cms.InputTag('packedPFCandidates')
#  getattr(process,'patPFMetTxyMuEGClean').parameters = multPhiCorrParams_Txy_25ns
#else:
#  postfix=postfixMuClean
#  #getattr(process,'patPFMetTxyCorr{0}'.format(postfix)).vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices')
#  #getattr(process,'patPFMetTxyCorr{0}'.format(postfix)).srcPFlow = cms.InputTag('packedPFCandidates')
#  getattr(process,'patPFMetTxyCorr{0}'.format(postfix)).parameters = multPhiCorrParams_Txy_25ns
#process.rootTuplePFMETType1CorNotRecorrected = process.rootTuplePFMETType1Cor.clone()
#process.rootTuplePFMETType1CorNotRecorrected.Suffix = 'Type1CorNotRecorrected'
#process.rootTuplePFMETType1XYCorNotRecorrected = process.rootTuplePFMETType1XYCor.clone()
#process.rootTuplePFMETType1XYCorNotRecorrected.Suffix = 'Type1XYCorNotRecorrected'
#process.rootTuplePFMETType1Cor.InputTag = cms.InputTag('slimmedMETs'+postfix)
#process.rootTuplePFMETType01Cor.InputTag = cms.InputTag('slimmedMETs'+postfix)
#process.rootTuplePFMETType1XYCor.InputTag = cms.InputTag('slimmedMETs'+postfix)
#process.rootTuplePFMETType01XYCor.InputTag = cms.InputTag('slimmedMETs'+postfix)
#
## ntuplize the new MET shifts
#process.rootNTupleNewMETs = cms.Path()
#allowedShifts = ['jres','jes','mes','ees','tes','ues']
#allowedSigns = ['+','-']
#signMap = {
#    '+' : 'Up',
#    '-' : 'Down',
#}
#collMap = {
#    # TODO
#    #'jres' : {'Jets'     : 'shiftedPatJetRes{sign}{postfix}'},
#    'jres' : {},
#    'jes'  : {'Jets'     : 'shiftedPatJetEn{sign}{postfix}'},
#    'mes'  : {'Muons'    : 'shiftedPatMuonEn{sign}{postfix}'},
#    'ees'  : {'Electrons': 'shiftedPatElectronEn{sign}{postfix}'},
#    'tes'  : {'Taus'     : 'shiftedPatTauEn{sign}{postfix}'},
#    'ues'  : {},
#    'pes'  : {},
#    }
#metMap = {
#  'jres' : 'patPFMetT1JetRes{sign}{postfix}',
#  'jes'  : 'patPFMetT1JetEn{sign}{postfix}',
#  'mes'  : 'patPFMetT1MuonEn{sign}{postfix}',
#  'ees'  : 'patPFMetT1ElectronEn{sign}{postfix}',
#  'tes'  : 'patPFMetT1TauEn{sign}{postfix}',
#  'ues'  : 'patPFMetT1UnclusteredEn{sign}{postfix}',
#  'pes'  : '',
#}
##mettypes = ['CaloMET','CaloMETType1Cor','PFMET','PFMETType1Cor','PFMETType01Cor','PFMETType01XYCor','PFMETPuppi','PFMETPuppiType1Cor']
##mettypes = ['PFMETType1XYCor','PFMETPuppiType1Cor']
#mettypes = ['PFMETType1XYCor']
#for shift in allowedShifts:
#    for sign in allowedSigns:
#      # FIXME TODO
#        ## embed shifted objects
#        #for coll in collMap[shift]:
#        #    modName = '{shift}{sign}{coll}Embedding'.format(shift=shift,sign=signMap[sign],coll=coll)
#        #    pluginName = 'MiniAODShifted{coll}Embedder'.format(coll=coll[:-1])
#        #    dName = coll.lower()
#        #    srcName = fs_daughter_inputs[dName]
#        #    shiftSrcName = collMap[shift][coll].format(sign=signMap[sign],postfix=postfix)
#        #    label = '{shift}{sign}{coll}'.format(shift=shift,sign=signMap[sign],coll=coll)
#        #    module = cms.EDProducer(
#        #        pluginName,
#        #        src = cms.InputTag(srcName),
#        #        shiftSrc = cms.InputTag(shiftSrcName),
#        #        label = cms.string(label),
#        #    )
#        #    setattr(process,modName,module)
#        #    fs_daughter_inputs[dName] = modName
#        #    #process.embedShifts *= getattr(process,shiftSrcName)
#        #    process.embedShifts *= getattr(process,modName)
#        # embed shifted met
#        for mettype in mettypes:
#          modName = 'rootTuple{mettype}{shift}{sign}'.format(mettype=mettype,shift=shift,sign=signMap[sign])
#          metName = metMap[shift].format(sign=signMap[sign],postfix=postfix)
#          prefix = '{mettype}{shift}{sign}'.format(mettype=mettype,shift=shift,sign=signMap[sign])
#          suffix = ''
#          module = cms.EDProducer(
#              'RootTupleMakerV2_MET',
#              InputTag = cms.InputTag(metName),
#              Prefix = cms.string(prefix),
#              Suffix = cms.string(suffix),
#              StoreUncorrectedMET = cms.bool(False), # this won't work either
#              StoreMETSignificance = cms.bool(False),
#              Uncertainty = cms.string('NoShift'),
#              CorrectionLevel = cms.string('NoCorrection'), # call pt(), etc., instead of shiftedPt()
#          )
#          setattr(process,modName,module)
#          # add shifted MET module
#          process.rootNTupleNewMETs *= getattr(process,metName)
#          process.rootNTupleNewMETs *= getattr(process,modName)
##process.schedule.append(process.rootNTupleNewMETs)#FIXME do we need this?

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
#process.load('FWCore.Modules.printContent_cfi')
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)
process.p = cms.Path(
    # supporting producers
    process.unpackedPatTrigger*
    ## L+J Filter
    process.LJFilter*  
    process.pdfWeights*
    ## Put everything into the tree
    process.rootTuplePFCandidates*
    process.rootTupleEvent*
    process.rootTupleElectronsSequence*
    process.rootTupleMuonsSequence*
    process.rootTupleVertex*
    process.rootTuplePFJetsSequence*
    #process.rootNTuplePFMET*
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

#process.endPath = cms.EndPath(
# )
##----------------------------------------------------------------------------------------------------
## Dump if necessary
##----------------------------------------------------------------------------------------------------
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
#        fileName = cms.untracked.string('dump.root')
#        )
#process.DUMP    = cms.EndPath (process.dump)

#----------------------------------------------------------------------------------------------------
# Run the path
#----------------------------------------------------------------------------------------------------
# Delete predefined Endpath (needed for running with CRAB)
del process.out
del process.outpath

#print process.dumpPython()
