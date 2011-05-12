# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

import os

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

# Make sure a correct global tag is used (please refer to https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Valid_Global_Tags_by_Release)
process.GlobalTag.globaltag = 'START311_V2::All'
                               
# Events to process
process.maxEvents.input = 10

# Options and Output Report
process.options.wantSummary = True

# Input files
process.source.fileNames = [
#    '/store/relval/CMSSW_4_1_4/RelValTTbar/GEN-SIM-RECO/START311_V2-v1/0019/62AC26CD-4161-E011-BDE8-002618943857.root' #RECO
    'rfio:///castor/cern.ch/user/d/darinb/RECO_AOD_Examples/Spring11_TTbar_Example_F69DC08F-4A4E-E011-B44E-E0CB4E1A117B.root'
    #'/store/relval/CMSSW_4_1_5/RelValTTbar_Tauola/GEN-SIM-RECO/START311_V2_PU_E7TeV_AVE_2_BX156-v1/0049/EC2A2471-0472-E011-9235-0018F3D09682.root' #RECO (pile-up)
]

# Add tcMET and pfMET
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

# Add type-I corrected pfMET
from Leptoquarks.RootTupleMakerV2.tools import *
addPfMETType1Cor(process, 'PFType1Cor')



# Add JEC Fall10
# See https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCor2010
# -->remember to do from $CMSSW_BASE/src: cvs co -d JEC UserCode/KKousour/data/Jec10V3.db
# Provide external jecfile:
jecfile = os.getenv('CMSSW_BASE')+'/src/JEC/Jec10V3.db' #using env variables
#jecfile = '/afs/cern.ch/user/s/santanas/scratch0/Releases/CMSSW_4_1_5_LQ/src/JEC/Jec10V3.db' #or giving full path

process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.jec = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
        ),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
               cms.PSet(
                       record = cms.string('JetCorrectionsRecord'),
                       tag    = cms.string('JetCorrectorParametersCollection_Jec10V3_AK5PF'),
                       label  = cms.untracked.string('AK5PF')
                       ),
               cms.PSet(
                       record = cms.string('JetCorrectionsRecord'),
                       tag    = cms.string('JetCorrectorParametersCollection_Jec10V3_AK5Calo'),
                       label  = cms.untracked.string('AK5Calo')
                       )
      ),
      ## here you add as many jet types as you need (AK5Calo, AK5JPT, AK7PF, AK7Calo, KT4PF, KT4Calo, KT6PF, KT6Calo)
      connect = cms.string('sqlite_file:'+jecfile)
)
process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

from PhysicsTools.PatAlgos.tools.jetTools import *
# Add PF jets
addJetCollection(process,cms.InputTag('ak5PFJets'),
    'AK5', 'PF',
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5PF', cms.vstring(['L1Offset','L2Relative', 'L3Absolute'])),
    doType1MET   = False,
    doL1Cleaning = False,
    doL1Counters = False,
    genJetCollection=cms.InputTag("ak5GenJets"),
    doJetID      = False
)

## Modify JEC for CaloJets (default)
process.patJetCorrFactors.levels = cms.vstring('L1Offset', 
                                               'L2Relative', 
                                               'L3Absolute')
                                               
# HEEPify PAT electrons
from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
process.heepPatElectrons = cms.EDProducer("HEEPAttStatusToPAT",
    eleLabel = cms.InputTag("patElectrons"),
    barrelCuts = cms.PSet(heepBarrelCuts),
    endcapCuts = cms.PSet(heepEndcapCuts)
)

# Switch HLT Input tag for Spring11 Reprocessing
process.rootTupleTrigger.HLTInputTag = cms.InputTag('TriggerResults','','REDIGI311X')

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
#### Shared Muon/Electron Skim
process.LJFilter.muLabel = 'muons'
process.LJFilter.elecLabel = 'gsfElectrons'
process.LJFilter.jetLabel = 'ak5CaloJets'

process.LJFilter.elecPT = 20.
process.LJFilter.muPT = 20.
process.LJFilter.counteitherleptontype = True


# Load HBHENoiseFilterResultProducer
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')

# Load EcalSeverityLevelESProducer (needed only if the SuperCluster module is run)
#process.load('RecoLocalCalo/EcalRecAlgos/EcalSeverityLevelESProducer_cfi')

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
        'keep *_rootTuplePileUp_*_*',
        #'keep *_rootTupleSuperClusters_*_*', #RECO only
        'keep *_rootTupleTrigger_*_*',
        'keep *_rootTupleVertex_*_*',
        'keep *_rootTupleGenEventInfo_*_*',
        'keep *_rootTupleGenParticles_*_*',
        'keep *_rootTupleGenJets_*_*',
        'keep *_rootTupleGenMETTrue_*_*',
        'keep *_rootTuplePhotons_*_*'
    )
)

# Produce PDF weights (maximum is 3)
process.pdfWeights = cms.EDProducer("PdfWeightProducer",
    # Fix POWHEG if buggy (this PDF set will also appear on output,
    # so only two more PDF sets can be added in PdfSetNames if not "")
    #FixPOWHEG = cms.untracked.string("cteq66.LHgrid"),
    #GenTag = cms.untracked.InputTag("genParticles"),
    PdfInfoTag = cms.untracked.InputTag("generator"),
    PdfSetNames = cms.untracked.vstring(
            "cteq66.LHgrid"
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
    process.rootTuplePileUp+
    #process.rootTupleSuperClusters+ #RECO only
    process.rootTupleTrigger+
    process.rootTupleVertex+
    process.rootTupleGenEventInfo+
    process.rootTupleGenParticles+
    process.rootTupleGenJets+
    process.rootTupleGenMETTrue+
    process.rootTuplePhotons
    )
    *process.rootTupleTree
)

# Delete predefined Endpath (needed for running with CRAB)
del process.out
del process.outpath

# Schedule definition
process.schedule = cms.Schedule(process.p)