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
                                   fileName = cms.string('RootTupleMakerV2_output_DATA_HT150_OR_HT200_OR_HT250__PhysicsDSTContent__179959-180252.root')
                                   
)

# Make sure a correct global tag is used (please refer to https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Valid_Global_Tags_by_Release)
process.GlobalTag.globaltag = 'GR_R_42_V19::All' # ===> First complete JEC set for 42x 2011 data (https://indico.cern.ch/getFile.py/access?contribId=8&resId=0&materialId=slides&confId=143981)
#process.GlobalTag.globaltag = 'GR_R_42_V12::All' # ===> for 42X prompt reco or re-reco (contains Jec11_V1, does not contain "residual" JEC and uncertainties yet...)
#process.GlobalTag.globaltag = 'GR_R_41_V0::All' # ===> for 41X prompt reco (contains Jec10V3)

# Events to process
process.maxEvents.input = -1

# Options and Output Report
process.options.wantSummary = True

# Input files
#########
######### ==> IMPORTANT : MAKE SURE THAT HLT INPUT TAG IN python/RootTupleMakerV2_PhysicsDSTStream2011_cfi.py IS CORRECT ("HLT" or "TEST" or ...)
#########

#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("180072:80")

process.source.fileNames = [
#    'file:/tmp/santanas/CC1F559B-1800-E111-97EC-003048F01E88.root'
#    '/store/data/Run2011B/PhysicsDST/RAW/v1/000/179/959/CC1F559B-1800-E111-97EC-003048F01E88.root'
#    'file:/data/santanas/Releases/CMSSW_4_2_9_HLT3_hltpatch3/src/outputPhysicsDST.root'
#        'file:/data/santanas/Releases/CMSSW_4_2_9_HLT3_hltpatch3/src/outputPhysicsDST_lowHTskim.root',
'file:/data/santanas/Data/HT150_OR_HT200_OR_HT250__PhysicsDSTContent__179959-180252/outputPhysicsDST_179959.root',
'file:/data/santanas/Data/HT150_OR_HT200_OR_HT250__PhysicsDSTContent__179959-180252/outputPhysicsDST_179977.root',
'file:/data/santanas/Data/HT150_OR_HT200_OR_HT250__PhysicsDSTContent__179959-180252/outputPhysicsDST_180072.root',
'file:/data/santanas/Data/HT150_OR_HT200_OR_HT250__PhysicsDSTContent__179959-180252/outputPhysicsDST_180076.root',
'file:/data/santanas/Data/HT150_OR_HT200_OR_HT250__PhysicsDSTContent__179959-180252/outputPhysicsDST_180093.root',
'file:/data/santanas/Data/HT150_OR_HT200_OR_HT250__PhysicsDSTContent__179959-180252/outputPhysicsDST_180241.root',
'file:/data/santanas/Data/HT150_OR_HT200_OR_HT250__PhysicsDSTContent__179959-180252/outputPhysicsDST_180250.root',
'file:/data/santanas/Data/HT150_OR_HT200_OR_HT250__PhysicsDSTContent__179959-180252/outputPhysicsDST_180252.root'
]

# RootTupleMakerV2 tree
process.rootTupleTree = cms.EDAnalyzer("RootTupleMakerV2_Tree",
    outputCommands = cms.untracked.vstring(
        'drop *',
        #'keep *_rootTupleEvent_*_*',
        'keep *_rootTuplePhysicsDSTStream2011_*_*',
    )
)

# Path definition
process.p = cms.Path(
    process.rootTuplePhysicsDSTStream2011
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
