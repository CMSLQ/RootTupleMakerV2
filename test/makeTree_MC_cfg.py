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
                                                   ##make sure to check the GT from DBS
process.GlobalTag.globaltag = 'START3X_V26A::All'  ##or use "edmProvDump inputfile.root > f.txt"
                                                   ##and search for GlobalTag in f.txt

# process.load("Configuration/StandardSequences/RawToDigi_Data_cff")
# process.load("L1Trigger/Configuration/L1RawToDigi_cff")
process.load("Configuration/StandardSequences/ReconstructionCosmics_cff")
process.load("RecoMuon/Configuration/RecoMuon_cff")


### output
process.add_( cms.Service( "TFileService",
fileName = cms.string("MC_RECO.root"),  ##give a name to the output file
                      closeFileFast = cms.untracked.bool(True)  ) )

### N events        
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) ) ##how many events you process?

### input 
process.source = cms.Source (
    "PoolSource",
    fileNames = cms.untracked.vstring( ## input files
    "/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V26A_356ReReco-v1/0009/FEBF7874-EF3D-DF11-910D-002354EF3BDF.root"
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

#### trigger
#process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
#from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed
#process.bit40OR41 = hltLevel1GTSeed.clone(L1TechTriggerSeeding = cms.bool(True), L1SeedsLogicalExpression = cms.string('40 OR 41'))


##########
##### PATH
##########


process.rootTupleTree = cms.EDAnalyzer("RootTupleMakerV2_Tree",
    outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_rootTupleEvent_*_*',
    ))

process.theBigNtuple = cms.Path(
    #process.bit40OR41 *
    (
    #process.rootTupleEvent +
    process.rootTupleEvent
    )
    * process.rootTupleTree )

