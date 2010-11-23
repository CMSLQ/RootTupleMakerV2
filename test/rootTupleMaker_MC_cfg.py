# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.default.limit = 10
#################################################################

# Load RootTupleMakerV2 modules
process.load('Leptoquarks.RootTupleMakerV2.Ntuple_cff')

# Output ROOT file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('RootTupleMakerV2_output_MC.root')
)

# Global tag (make sure it always matches with the global tag used to reconstruct input files)
process.GlobalTag.globaltag = 'START38_V12::All'

# Events to process
process.maxEvents.input = 100

# Options and Output Report
process.options.wantSummary = True

# Input files
process.source.fileNames = [
    '/store/mc/Fall10/TTJets_TuneZ2_7TeV-madgraph-tauola/GEN-SIM-RECO/START38_V12-v2/0006/125D724F-5EE4-DF11-816A-00145E550F3E.root'
]

# Turn off MC matching for the process
#removeMCMatching(process, ['All'])

# Add tcMET and pfMET
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

# Add pfMET type 1 corrections
from Leptoquarks.RootTupleMakerV2.tools import *
addPfMETType1Cor(process, 'PFType1Cor')

# Get the 7 TeV GeV jet corrections
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJECSet( process, "Spring10" )
process.ak5CaloL2Relative.useCondDB = False
process.ak5CaloL3Absolute.useCondDB = False
process.ak5CaloResidual.useCondDB = False
process.ak5PFL2Relative.useCondDB = False
process.ak5PFL3Absolute.useCondDB = False
process.ak5PFResidual.useCondDB = False

# Add PF jets
addJetCollection(process,cms.InputTag('ak5PFJets'),
    'AK5', 'PF',
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5','PF'),
    doType1MET   = False,
    doL1Cleaning = False,
    doL1Counters = False,
    genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = False
)

# Switch on PAT trigger
#from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
#switchOnTrigger( process )

# Restrict input to AOD
#restrictInputToAOD(process, ['All'])

# HEEPify PAT electrons
from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
process.heepPatElectrons = cms.EDProducer("HEEPAttStatusToPAT",
    eleLabel = cms.InputTag("patElectrons"),
    barrelCuts = cms.PSet(heepBarrelCuts),
    endcapCuts = cms.PSet(heepEndcapCuts)
)

# Add 'heepPatElectrons' in the right place and point 'selectedLayer1Electrons' to them
process.patDefaultSequence.replace( process.patElectrons, process.patElectrons*process.heepPatElectrons )
process.selectedPatElectrons.src = cms.InputTag("heepPatElectrons")

# Electron and jet cleaning deltaR parameters
process.cleanPatElectrons.checkOverlaps.muons.deltaR = 0.3
process.cleanPatJets.checkOverlaps.muons.deltaR = 0.5
process.cleanPatJets.checkOverlaps.electrons.deltaR = 0.5

# Skim definition
process.load("Leptoquarks.LeptonJetFilter.leptonjetfilter_cfi")
##################################################################
#### Electron based skim
process.LJFilter.muLabel = 'muons'
process.LJFilter.elecLabel = 'gsfElectrons'
process.LJFilter.jetLabel = 'ak5CaloJets'
process.LJFilter.muonsMin = -1
process.LJFilter.electronsMin = 1
process.LJFilter.elecPT = 20.
process.LJFilter.counteitherleptontype = False
##################################################################
#### SuperCluster based skim
#process.LJFilter.muLabel = 'muons'
#process.LJFilter.elecLabel = 'gsfElectrons'
#process.LJFilter.jetLabel = 'ak5CaloJets'
#process.LJFilter.muonsMin = -1
#process.LJFilter.electronsMin = -1
#process.LJFilter.scMin = 1
#process.LJFilter.scET = 20.
#process.LJFilter.scHoE = 0.05
##################################################################

# Load HBHENoiseFilterResultProducer
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')

# RootTupleMakerV2 tree
process.rootTupleTree = cms.EDAnalyzer("RootTupleMakerV2_Tree",
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_rootTupleEvent_*_*',
        'keep *_rootTupleEventSelection_*_*',
        'keep *_rootTupleCaloJets_*_*',
        'keep *_rootTuplePFJets_*_*',
        'keep *_rootTupleElectrons_*_*',
        'keep *_rootTupleCaloMET_*_*',
        'keep *_rootTupleTCMET_*_*',
        'keep *_rootTuplePFMET_*_*',
        'keep *_rootTuplePFMETType1Cor_*_*',
        'keep *_rootTupleMuons_*_*',
        'keep *_rootTupleSuperClusters_*_*',
        'keep *_rootTupleTrigger_*_*',
        'keep *_rootTupleVertex_*_*',
        'keep *_rootTupleGenEventInfo_*_*',
        'keep *_rootTupleGenParticles_*_*',
        'keep *_rootTupleGenJets_*_*',
        'keep *_rootTupleGenMETTrue_*_*'
    )
)

# Produce PDF weights (maximum is 3)
process.pdfWeights = cms.EDProducer("PdfWeightProducer",
    FixPOWHEG = cms.untracked.bool(False), # fix POWHEG (it requires cteq6m* PDFs in the list)
    PdfInfoTag = cms.untracked.InputTag("generator"),
    PdfSetNames = cms.untracked.vstring(
        "cteq65.LHgrid"
      #, "MRST2006nnlo.LHgrid"
      #, "MRST2007lomod.LHgrid"
    )
)

# In order to disable the PDF weights calculation, uncomment the line below and
# comment out the pdfWeights module in the Path 'p' below
#process.rootTupleGenEventInfo.StorePDFWeights = False

# Path definition
process.p = cms.Path(
    process.LJFilter*
    process.HBHENoiseFilterResultProducer*
    process.pdfWeights*
    process.metJESCorAK5PFJet*
    process.patDefaultSequence*
    (
    process.rootTupleEvent+
    process.rootTupleEventSelection+
    process.rootTupleCaloJets+
    process.rootTuplePFJets+
    process.rootTupleElectrons+
    process.rootTupleCaloMET+
    process.rootTupleTCMET+
    process.rootTuplePFMET+
    process.rootTuplePFMETType1Cor+
    process.rootTupleMuons+
    process.rootTupleSuperClusters+
    process.rootTupleTrigger+
    process.rootTupleVertex+
    process.rootTupleGenEventInfo+
    process.rootTupleGenParticles+
    process.rootTupleGenJets+
    process.rootTupleGenMETTrue
    )
    *process.rootTupleTree
)

# Delete predefined Endpath (needed for running with CRAB)
del process.out
del process.outpath

# Schedule definition
process.schedule = cms.Schedule(process.p)
