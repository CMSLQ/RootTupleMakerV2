import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt

process = cms.Process('TEST')
process.load('Leptoquarks.RootTupleMakerV2.ntuple_cff')

process.load("Configuration/StandardSequences/Geometry_cff")
process.load("Configuration/StandardSequences/MagneticField_cff")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")

### global tag
#process.GlobalTag.globaltag = 'GR_R_35X_V6::All'  ##make sure to check the GT from DBS
process.GlobalTag.globaltag = 'GR10_P_V4::All'     ##or use "edmProvDump inputfile.root > f.txt"
                                                   ##and search for GlobalTag in f.txt


# process.load("Configuration/StandardSequences/RawToDigi_Data_cff")
# process.load("L1Trigger/Configuration/L1RawToDigi_cff")
process.load("Configuration/StandardSequences/ReconstructionCosmics_cff")
process.load("RecoMuon/Configuration/RecoMuon_cff")


### output
process.add_( cms.Service( "TFileService",
fileName = cms.string("DATA_RECO.root"),  ##give a name to the output file
                      closeFileFast = cms.untracked.bool(True)  ) )

### N events        
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) ) ##how many events you process?

### input 
process.source = cms.Source (
    "PoolSource",
    fileNames = cms.untracked.vstring( ## input files
    '/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/132/605/F4C295F3-B041-DF11-94EF-003048D46060.root'
    #'/store/data/Commissioning10/MinimumBias/RECO/Apr1ReReco-v2/0140/FEAE8844-6E40-DF11-92E9-0026189438E8.root'
    #'/store/mc/Summer09/MinBias/GEN-SIM-RECO/STARTUP3X_V8D_900GeV-v1/0005/E4590360-4CD7-DE11-8CB4-002618943896.root'
    ),
    
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    
    secondaryFileNames = cms.untracked.vstring() )

### message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.default.limit = 100

### jet corrections
#process.load("JetMETCorrections.Configuration.L2L3Corrections_Summer09_7TeV_ReReco332_cff")
#process.load("JetMETCorrections.Configuration.L2L3Corrections_900GeV_cff")
#process.load("JetMETCorrections.Configuration.L2L3Corrections_2360GeV_cff")

### summary report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

### trigger
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed
process.bit40OR41 = hltLevel1GTSeed.clone(L1TechTriggerSeeding = cms.bool(True), L1SeedsLogicalExpression = cms.string('40 OR 41'))

### physics declared bit - FIX THIS --> see https://hypernews.cern.ch/HyperNews/CMS/get/hlt/1421.html
#from HLTrigger.HLTfilters.hltHighLevelDev_cfi import hltHighLevelDev
#process.physDecl = hltHighLevelDev.clone(HLTPaths = ['HLT_PhysicsDeclared'], HLTPathsPrescales = [1])


##########
##### PATH
##########


process.rootTupleTree = cms.EDAnalyzer("RootTupleMakerV2_Tree",
    outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_rootTupleEvent_*_*',
    ))

process.theBigNtuple = cms.Path(
    #########process.physDecl * FIX THIS
    process.bit40OR41 *
    (
    #process.rootTupleEvent +
    process.rootTupleEvent
    )
    * process.rootTupleTree )

