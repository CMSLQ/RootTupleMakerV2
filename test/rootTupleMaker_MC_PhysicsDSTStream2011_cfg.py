# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

from PhysicsTools.PatAlgos.tools.coreTools import *

import os

############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.default.limit = 10
#################################################################

# Load RootTupleMakerV2 modules
process.load('Leptoquarks.RootTupleMakerV2.Ntuple_cff')

# Output ROOT file
process.TFileService = cms.Service("TFileService",
                                   #fileName = cms.string('RootTupleMakerV2_output_DATA_PhysicsDSTStream2011.root')
                                   fileName = cms.string('QstarToJJ_M-500_TuneD6T_7TeV_pythia6__Fall11-PU_S6_START42_V14B-v1__GEN-RAW_1_1_bbb.root')
                                   #fileName = cms.string('QstarToJJ_M-700_TuneD6T_7TeV_pythia6__Fall11-PU_S6_START42_V14B-v1__GEN-RAW_1_1_bbb.root')
                                   #fileName = cms.string('QstarToJJ_M-1200_TuneD6T_7TeV_pythia6__Fall11-PU_S6_START42_V14B-v1__GEN-RAW_1_1_bbb.root')                                   
)

# Make sure a correct global tag is used (please refer to https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Valid_Global_Tags_by_Release)
process.GlobalTag.globaltag = 'START42_V13::All' # ===> First complete JEC set for 42x 2011 data (https://indico.cern.ch/getFile.py/access?contribId=8&resId=0&materialId=slides&confId=143981)
#process.GlobalTag.globaltag = 'START42_V12::All' # ===> for Summer11 MC analyzed in 42X (contains Jec11_V1, does not contain "residual" JEC and uncertainties yet...)
#process.GlobalTag.globaltag = 'START41_V0::All' # ===> for 41X MC analyzed in 41X (contains Jec10V3)

# Events to process
process.maxEvents.input = -1

# Options and Output Report
process.options.wantSummary = True

# Input files
#########
######### ==> IMPORTANT : MAKE SURE THAT HLT INPUT TAG IN python/RootTupleMakerV2_PhysicsDSTStream2011_cfi.py IS CORRECT ("HLT" or "TEST" or ...)
#########

process.source.fileNames = [
'file:/data/santanas/Releases/CMSSW_4_2_9_HLT3_hltpatch3/src/myoutput2.root'
# they do not include gen info (used for preliminary studies)
#'file:/data/santanas/Data/QstarToJJ_PhysicsDST_ForTest_16_11_2011/outputPhysicsDST_QstarToJJ_M-500.root',
#'file:/data/santanas/Data/QstarToJJ_PhysicsDST_ForTest_16_11_2011/outputPhysicsDST_QstarToJJ_M-700.root',
#'file:/data/santanas/Data/QstarToJJ_PhysicsDST_ForTest_16_11_2011/outputPhysicsDST_QstarToJJ_M-1200.root'
]

# RootTupleMakerV2 tree
process.rootTupleTree = cms.EDAnalyzer("RootTupleMakerV2_Tree",
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_rootTuplePhysicsDSTStream2011_*_*',
        'keep *_rootTupleGenEventInfo_*_*',
        'keep *_rootTupleGenParticles_*_*'
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

# Path definition
process.p = cms.Path(
    process.pdfWeights
    *process.rootTuplePhysicsDSTStream2011
    *process.rootTupleGenEventInfo
    *process.rootTupleGenParticles
    *process.rootTupleTree
)

###################
#process.dump = cms.OutputModule("PoolOutputModule",
#                                outputCommands = cms.untracked.vstring(
#                                'keep *',
#                                ),
#                                fileName = cms.untracked.string('dump.root')
#                                )
#process.DUMP    = cms.EndPath (process.dump)
###################

# Delete predefined Endpath (needed for running with CRAB)
del process.out
del process.outpath

# Schedule definition
process.schedule = cms.Schedule(process.p)

##############
#process.schedule.append( process.DUMP )
##############
