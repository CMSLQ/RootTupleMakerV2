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
process.options.allowUnscheduled = cms.untracked.bool(False)
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
    #fileName = cms.string( "file_m650.root" )
    fileName = cms.string( "file_MG.root" )
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
  process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2_v1'
else:
  #process.GlobalTag.globaltag = '80X_dataRun2_Prompt_ICHEP16JEC_v0'#for 2016 H Prompt
  process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v4'# for 23Sep ReReco (2016B-G)
# feed it into the ntuple
process.rootTupleEvent.globalTag = process.GlobalTag.globaltag
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')


# Events to process
process.maxEvents.input = 500

# Input files
process.source.fileNames = [
    # specified by InputList.txt
    # LQToUE-M550
    #'/store/mc/RunIIFall15MiniAODv2/LQToUE_ENuJJFilter_M-550_BetaHalf_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/20000/9CBC6301-61B8-E511-BC05-001E67398C5A.root'
    # LQToCMu-M1000
    #'/store/mc/RunIIFall15MiniAODv2/LQToCMu_M-1000_BetaOne_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/AA32F756-7BB8-E511-96C5-001E672488A4.root'
    # Run2015D data
    #'/store/data/Run2015D/SingleMuon/MINIAOD/16Dec2015-v1/10000/00006301-CAA8-E511-AD39-549F35AD8BC9.root'

    #'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/30000/00AEB2F6-541B-E611-AF2A-0025905C42F2.root'
    #'/store/mc/RunIISpring16MiniAODv2/ZZ_TuneCUETP8M1_13TeV-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/DA6AE5C6-151B-E611-9EEF-782BCB536A50.root'
    #'/store/mc/RunIISpring16MiniAODv2/LQLQToTopMu_M-1000_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/10000/46EDE553-8624-E611-AFB2-00259073E380.root'
    #'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/80000/82E9F815-F958-E611-968C-00266CFF0234.root'
    #'file:/afs/cern.ch/user/s/scooper/work/private/cmssw/8011/TestRootNTuplizerRecipe/src/LQToUE_M-1000_BetaOne-RunIISpring16MiniAODv2.root'
    #'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/10000/0A23534F-F638-E611-8ED1-782BCB78621D.root'
    #'/store/data/Run2016G/SingleMuon/MINIAOD/PromptReco-v1/000/278/820/00000/1C7D6E1E-6564-E611-B75B-FA163EDC3FC8.root'
  #'/store/mc/RunIISpring16MiniAODv2/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/70000/0468BB7B-D53F-E611-8647-0CC47A0093EC.root'',
    #'/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/premix_withHLT_80X_mcRun2_asymptotic_v14-v1/00000/042D62D1-C597-E611-8FA4-549F3525B9A0.root'                                                                

    # sep23 rereco
    #'/store/data/Run2016B/SingleElectron/MINIAOD/23Sep2016-v3/00000/00099863-E799-E611-A876-141877343E6D.root'
    'file:0077B9D8-1E9A-E611-963A-0CC47A57CCE8.root'
    # amcatnlo DYJ inclusive
    #'/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/000FF6AC-9F2A-E611-A063-0CC47A4C8EB0.root',
    # amcatnlo WJets inclusive
    #'/store/mc/RunIISpring16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/00000/00BC765D-BD3D-E611-9343-0025905C2CD2.root'
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
process.load('PhysicsTools.PatAlgos.slimming.unpackedPatTrigger_cfi')
###### need to load it into the process, or else it won't run
process.unpackedPatTrigger = unpackedPatTrigger.clone()
#process.unpackedPatTrigger.triggerResults = 'TriggerResults::HLT2'
#process.rootTupleTriggerObjects.TriggerBitsTag = 'TriggerResults::HLT2'
#process.rootTupleTrigger.HLTInputTag = 'TriggerResults::HLT2'
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

process.schedule = cms.Schedule()
# Need the egamma smearing for 80X: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMSmearer
process.load('EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi')
correctionType = "Full2016_v1"
process.calibratedPatElectrons.isMC = varOptions.isMC
process.selectedElectrons = cms.EDFilter(
  "PATElectronSelector",
  src = cms.InputTag("slimmedElectrons"),
  cut = cms.string("pt > 5 && abs(eta)<2.5")
)
process.calibratedPatElectrons.electrons = cms.InputTag('selectedElectrons')
# ntuplize this corrected electron collection
process.rootTupleElectrons.InputTag = cms.InputTag('calibratedPatElectrons')
# HEEP 5.1/6.0
# Also load Egamma cut-based ID in new VID framework while we're at it
# See: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
#   and for HEEP: https://hypernews.cern.ch/HyperNews/CMS/get/egamma/1519/2/1/1/1.html
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
# Each of these two example IDs contains all four standard cut-based ID working points
my_id_modules = []
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff')
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff') # for 50 ns, 13 TeV data
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff')
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff')
#Add them to the VID producer
for idmod in my_id_modules:
  setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
# XXX NB, must be the same as input collection used for electron ntuplizer
process.egmGsfElectronIDs.physicsObjectSrc = process.rootTupleElectrons.InputTag
process.heepIDVarValueMaps.elesMiniAOD  = process.rootTupleElectrons.InputTag
process.electronRegressionValueMapProducer.srcMiniAOD = process.rootTupleElectrons.InputTag
process.electronMVAValueMapProducer.srcMiniAOD = process.rootTupleElectrons.InputTag#FIXME do we need this?

process.electronSupportPath = cms.Path()
process.electronSupportPath += process.selectedElectrons
process.electronSupportPath += process.calibratedPatElectrons
#process.electronSupportPath += process.egmGsfElectronIDs
process.electronSupportPath += process.egmGsfElectronIDSequence
process.schedule.append(process.electronSupportPath)


#----------------------------------------------------------------------------------------------------
# Add MVA electron ID
#
# MVA electron ID details on this twiki:
# https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
#
#----------------------------------------------------------------------------------------------------
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

# Loading  JEC/JER from local DB file
cmsswbase=os.environ['CMSSW_BASE']
dbJetMCDBFile = '{0}/src/Leptoquarks/RootTupleMakerV2/data/Spring16_23Sep2016V2_MC.db'.format(cmsswbase)
dbJetDataDBFile = '{0}/src/Leptoquarks/RootTupleMakerV2/data/Spring16_23Sep2016AllV2_DATA.db'.format(cmsswbase)

# Use JER from GR
process.rootTuplePFJetsAK4CHS.ReadJERFromGT = True
process.rootTuplePFJetsAK4Puppi.ReadJERFromGT = True

# load JES/etc.from local db file, and into global tag
jecDBPrefix='sqlite_file:'
jetDBFile = jecDBPrefix+dbJetMCDBFile if varOptions.isMC else jecDBPrefix+dbJetDataDBFile
process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
# JEC
jecPSetDataAK4chs = cms.PSet(
    record = cms.string('JetCorrectionsRecord'),
    tag    = cms.string('JetCorrectorParametersCollection_Spring16_23Sep2016AllV2_DATA_AK4PFchs'),
    label  = cms.untracked.string('AK4PFchs')
)
jecPSetMCAK4chs = cms.PSet(
    record = cms.string('JetCorrectionsRecord'),
    tag    = cms.string('JetCorrectorParametersCollection_Spring16_23Sep2016V2_MC_AK4PFchs'),
    label  = cms.untracked.string('AK4PFchs')
)
process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
        jecPSetMCAK4chs if varOptions.isMC else jecPSetDataAK4chs
      ), 
       connect = cms.string(jetDBFile)
)
## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')
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

##----------------------------------------------------------------------------------------------------
## MET Re-Corrections and Uncertainties
## https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
## https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#METSysTools
##----------------------------------------------------------------------------------------------------
postfix = "Recorrected"
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#myJecUncFile = jecUncFileMC if varOptions.isMC else jecUncFileData
myJecUncFile = '' # take from global tag
#default configuration for miniAOD reprocessing
runMetCorAndUncFromMiniAOD(process,
                           #jetCollUnskimmed='slimmedJets',
                           #jetCollUnskimmed='updatedPatJets',
                           #jetColl='slimmedJets',
                           isData=not varOptions.isMC,
                           #electronColl=cms.InputTag('slimmedElectrons'),
                           #repro74X=True,
                           #jecUncFile=myJecUncFile,
                           postfix=postfix                           
)
if varOptions.isMC:
  process.applyCorrections += getattr(process,'genMetExtractor{0}'.format(postfix))
process.applyCorrections += getattr(process,'patJetCorrFactorsReapplyJEC{0}'.format(postfix))
#process.applyCorrections += process.patJets
#process.applyCorrections += getattr(process,'selectedPatJets{0}'.format(postfix))
process.applyCorrections += getattr(process,'patJetsReapplyJEC{0}'.format(postfix))
process.applyCorrections += getattr(process,'basicJetsForMet{0}'.format(postfix))
process.applyCorrections += getattr(process,'jetSelectorForMet{0}'.format(postfix))
process.applyCorrections += getattr(process,'cleanedPatJets{0}'.format(postfix))
process.applyCorrections += getattr(process,'pfMet{0}'.format(postfix))
process.applyCorrections += getattr(process,'patPFMet{0}'.format(postfix))
process.applyCorrections += getattr(process,'patPFMetT1T2Corr{0}'.format(postfix))
process.applyCorrections += getattr(process,'patPFMetTxyCorr{0}'.format(postfix))
process.applyCorrections += getattr(process,'patPFMetT1{0}'.format(postfix))
process.applyCorrections += getattr(process,'patPFMetT1Txy{0}'.format(postfix))
process.applyCorrections += getattr(process,'patPFMetTxy{0}'.format(postfix))
process.schedule.append(process.applyCorrections)
# now schedule the MET corrections
process.metCorAndUncPath = cms.Path()
process.metCorAndUncPath += getattr(process,'fullPatMetSequence{0}'.format(postfix))
process.schedule.append(process.metCorAndUncPath)
# hacks to fix patPFMETCorrections_cff (used by runMETCorrectionsAndUncertainties)
process.patSmearedJets.src = cms.InputTag('updatedPatJets')
process.selectedPatJetsForMetT1T2Corr.src = cms.InputTag('updatedPatJets')
process.selectedPatJetsForMetT1T2CorrRecorrected.src = cms.InputTag('updatedPatJets')
process.selectedPatJetsForMetT2Corr.src = cms.InputTag('updatedPatJets')
process.selectedPatJetsForMetT2CorrRecorrected.src = cms.InputTag('updatedPatJets')
getattr(process,'patPFMetTxyCorr{0}'.format(postfix)).vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices')
getattr(process,'patPFMetTxyCorr{0}'.format(postfix)).srcPFlow = cms.InputTag('packedPFCandidates')

# ntuplize the newly-corrected MET
process.rootTuplePFMETType1CorNotRecorrected = process.rootTuplePFMETType1Cor.clone()
process.rootTuplePFMETType1CorNotRecorrected.Suffix = 'Type1CorNotRecorrected'
process.rootTuplePFMETType1XYCorNotRecorrected = process.rootTuplePFMETType1XYCor.clone()
process.rootTuplePFMETType1XYCorNotRecorrected.Suffix = 'Type1XYCorNotRecorrected'
process.rootTuplePFMETType1Cor.InputTag = cms.InputTag('slimmedMETs'+postfix)
process.rootTuplePFMETType01Cor.InputTag = cms.InputTag('slimmedMETs'+postfix)
process.rootTuplePFMETType1XYCor.InputTag = cms.InputTag('slimmedMETs'+postfix)
process.rootTuplePFMETType01XYCor.InputTag = cms.InputTag('slimmedMETs'+postfix)

# ntuplize the new MET shifts
process.rootNTupleNewMETs = cms.Path()
allowedShifts = ['jres','jes','mes','ees','tes','ues']
allowedSigns = ['+','-']
signMap = {
    '+' : 'Up',
    '-' : 'Down',
}
collMap = {
    # TODO
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
mettypes = ['PFMETType1XYCor','PFMETPuppiType1Cor']
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
          prefix = '{mettype}{shift}{sign}'.format(mettype=mettype,shift=shift,sign=signMap[sign])
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
          # add shifted MET module
          process.rootNTupleNewMETs *= getattr(process,metName)
          process.rootNTupleNewMETs *= getattr(process,modName)
#process.schedule.append(process.rootNTupleNewMETs)#FIXME do we need this?

#----------------------------------------------------------------------------------------------------
# This is MC, so analyze the smeared PFJets by default
# But smearing should be done by default in MiniAOD ?
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
process.LJFilter
process.LJFilter.muPT     = 8.0
process.LJFilter.electronsMin = 0
process.LJFilter.elecPT       = 12.0
process.LJFilter.tausMin = 0
process.LJFilter.tauPT   = 12.0
process.LJFilter.jetsMin = 0
process.LJFilter.jetPT   = 12.0
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
        #GenTag = cms.untracked.InputTag("genParticles"),
	PdfInfoTag = cms.untracked.InputTag("generator"),
	PdfSetNames = cms.untracked.vstring(
            "PDF4LHC15_nlo_100.LHgrid" ,
            "CT10nlo.LHgrid" , 
            "MMHT2014nlo68cl.LHgrid",
            #"NNPDF30_nlo_as_0118.LHgrid"
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
        'keep *_rootTupleGenTaus*_*_*',
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

process.p = cms.Path(
    ## L+J Filter
    process.LJFilter*  
    process.pdfWeights*
    # supporting producers
    process.unpackedPatTrigger*
    ## Put everything into the tree
    process.rootTupleEvent*
    process.rootTupleEventSelectionSequence*
    process.rootTuplePFCandidates*
    process.rootTupleElectronsSequence*
    process.rootTupleMuonsSequence*
    process.rootTupleVertex*
    process.rootTuplePFJetsSequence*
    #process.rootNTuplePFMET*
    process.rootTuplePFMETSequence*
    process.rootTuplePFMETType1CorNotRecorrected*
    process.rootTupleTriggerObjects*
    process.rootTupleTrigger*
    process.rootTupleGenEventInfo*
    process.rootTupleGenParticles*
    process.rootTupleGenJetsSequence*
    process.rootTupleGenMETTrue*
    process.rootTupleTree
)

process.schedule.append(process.p)

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
