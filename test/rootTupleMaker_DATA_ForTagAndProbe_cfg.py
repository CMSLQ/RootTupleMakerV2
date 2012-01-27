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
    fileName = cms.string('RootTupleMakerV2_output_DATA.root')
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
process.source.fileNames = [
    #'file:/data/santanas/SingleElectron-AOD-May10ReReco-v1_42X.root'
    'file:/data/santanas/DoubleElectron_Run2011B_AOD_PromptReco-v1_180250.root' #AOD (42X)
    #'/store/data/Run2011A/SingleElectron/AOD/PromptReco-v4/000/165/121/D0162694-1382-E011-876E-003048F1BF68.root' #AOD (42X)
    #'/store/data/Run2011A/SingleElectron/AOD/PromptReco-v1/000/161/312/90646AF9-F957-E011-B0DB-003048F118C4.root' #AOD (41X)
    #'/store/data/Run2011A/SingleElectron/RECO/PromptReco-v1/000/160/405/0E58AE5B-D64F-E011-88F1-003048F024DC.root' #RECO
]

# Turn off MC matching for the process
removeMCMatching(process, ['All'], '', False)

# Add tcMET and pfMET
from PhysicsTools.PatAlgos.tools.metTools import *
addTcMET(process, 'TC')
addPfMET(process, 'PF')

# Add type-I corrected pfMET and charged pfMET
from Leptoquarks.RootTupleMakerV2.tools import *
addPfMETType1Cor(process, 'PFType1Cor')
addPfChargedMET(process, 'PFCharged')

## Add JEC Fall10 ==> Not needed if the corrections are already in the global tag <==
## See https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetEnCor2010
## -->remember to do from $CMSSW_BASE/src: cvs co -d JEC UserCode/KKousour/data/Jec10V3.db
## Provide external jecfile:
#jecfile = os.getenv('CMSSW_BASE')+'/src/JEC/Jec10V3.db' #using env variables
##jecfile = '/afs/cern.ch/user/s/santanas/scratch0/Releases/CMSSW_4_1_5_LQ/src/JEC/Jec10V3.db' #or giving full path
#
#process.load("CondCore.DBCommon.CondDBCommon_cfi")
#process.jec = cms.ESSource("PoolDBESSource",
#      DBParameters = cms.PSet(
#        messageLevel = cms.untracked.int32(0)
#        ),
#      timetype = cms.string('runnumber'),
#      toGet = cms.VPSet(
#               cms.PSet(
#                       record = cms.string('JetCorrectionsRecord'),
#                       tag    = cms.string('JetCorrectorParametersCollection_Jec10V3_AK5PF'),
#                       label  = cms.untracked.string('AK5PF')
#                       ),
#               cms.PSet(
#                       record = cms.string('JetCorrectionsRecord'),
#                       tag    = cms.string('JetCorrectorParametersCollection_Jec10V3_AK5Calo'),
#                       label  = cms.untracked.string('AK5Calo')
#                       )
#      ),
#      ## here you add as many jet types as you need (AK5Calo, AK5JPT, AK7PF, AK7Calo, KT4PF, KT4Calo, KT6PF, KT6Calo)
#      connect = cms.string('sqlite_file:'+jecfile)
#)
#process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

# For L1FastJet Corrections ( https://hypernews.cern.ch/HyperNews/CMS/get/physTools/2504/1.html 
# and https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections?redirectedfrom=CMS.WorkBookJetEnergyCorrections )

##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
##-------------------- Import the Jet RECO modules -----------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load('RecoJets.Configuration.RecoJets_cff')
##-------------------- Turn-on the FastJet density calculation -----------------------
process.kt6PFJets.doRhoFastjet = True
process.kt6CaloJets.doRhoFastjet = True
##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm -----------------------
process.ak5PFJets.doAreaFastjet = True
process.ak5CaloJets.doAreaFastjet = True
##-------------------- set rho to True to prevent Exceptions at run time, when L1FastJet is used -----------------------
process.patJetCorrFactors.useRho = True

# Add PF jets          --> See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#Jet_Tools
from PhysicsTools.PatAlgos.tools.jetTools import *
addJetCollection(process,cms.InputTag('ak5PFJets'),
    'AK5', 'PF',
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'])), 
    #jetCorrLabel = ('AK5PF', cms.vstring(['L1Offset','L2Relative', 'L3Absolute','L2L3Residual'])), 
    #jetCorrLabel = ('AK5PF', cms.vstring(['L1Offset','L2Relative', 'L3Absolute'])),
                 # check L1 corrections
                 # see https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPATDataFormats#JetCorrFactors
    doType1MET   = False,
    doL1Cleaning = False,
    doL1Counters = False,
    genJetCollection = cms.InputTag("ak5GenJets"),
    doJetID      = False
)
addJetCollection(process,cms.InputTag('ak5PFJets'),
    'AK5', 'PFL1Offset',
    doJTA        = True,
    doBTagging   = True,
    #jetCorrLabel = ('AK5PF', cms.vstring(['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'])), 
    jetCorrLabel = ('AK5PF', cms.vstring(['L1Offset','L2Relative', 'L3Absolute','L2L3Residual'])), 
    #jetCorrLabel = ('AK5PF', cms.vstring(['L1Offset','L2Relative', 'L3Absolute'])),
                 # check L1 corrections
                 # see https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPATDataFormats#JetCorrFactors
    doType1MET   = False,
    doL1Cleaning = False,
    doL1Counters = False,
    genJetCollection = cms.InputTag("ak5GenJets"),
    doJetID      = False
)
## Modify type-I corrected PFMET --> should be consistent with PF Jets (defined above)
#process.metJESCorAK5PFJet.corrector = cms.string('ak5PFL2L3') #default value
process.metJESCorAK5PFJet.corrector = cms.string('ak5PFL1FastL2L3Residual') 
#process.metJESCorAK5PFJet.corrector = cms.string('ak5PFL1L2L3Residual') 
#process.metJESCorAK5PFJet.corrector = cms.string('ak5PFL1L2L3')
process.metJESCorAK5PFJet.jetPTthreshold = cms.double(10.0) 
## Remove muons from raw pfjet collection
process.ak5PFJetsNoMuon =  cms.EDFilter("PFJetSelector",    
                                        src = cms.InputTag('ak5PFJets'),
                                        cut = cms.string("chargedMuEnergyFraction < 0.8"),
                                        filter = cms.bool(False)
                                        )
process.metJESCorAK5PFJet.inputUncorJetsLabel = cms.string('ak5PFJetsNoMuon')

## Modify JEC for CaloJets (default)
process.patJetCorrFactors.levels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual') 
#process.patJetCorrFactors.levels = cms.vstring('L1Offset', 'L2Relative', 'L3Absolute', 'L2L3Residual') 
#process.patJetCorrFactors.levels = cms.vstring('L1Offset', 'L2Relative', 'L3Absolute') 

addJetCollection(process,cms.InputTag('ak5CaloJets'),
    'AK5', 'CaloL1Offset',
    doJTA        = True,
    doBTagging   = True,
    #jetCorrLabel = ('AK5Calo', cms.vstring(['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'])), 
    jetCorrLabel = ('AK5Calo', cms.vstring(['L1Offset','L2Relative', 'L3Absolute','L2L3Residual'])), 
    #jetCorrLabel = ('AK5Calo', cms.vstring(['L1Offset','L2Relative', 'L3Absolute'])),
                 # check L1 corrections
                 # see https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPATDataFormats#JetCorrFactors
    doType1MET   = False,
    doL1Cleaning = False,
    doL1Counters = False,
    genJetCollection = cms.InputTag("ak5GenJets"),
    doJetID      = False
)
                                               
## Modify type-I corrected caloMET (default) --> should be consistent with caloJets (defined above)
#process.metJESCorAK5CaloJet.corrector = cms.string('ak5CaloL2L3') #default value
process.metJESCorAK5CaloJet.corrector = cms.string('ak5CaloL1FastL2L3Residual') 
#process.metJESCorAK5CaloJet.corrector = cms.string('ak5CaloL1L2L3Residual')
#process.metJESCorAK5CaloJet.corrector = cms.string('ak5CaloL1L2L3')
process.metJESCorAK5CaloJet.jetPTthreshold = cms.double(20.0) #default value

## Read JEC uncertainties (might not be available in some global tag)
process.rootTupleCaloJets.ReadJECuncertainty = True 
process.rootTuplePFJets.ReadJECuncertainty = True 

#OLD
# Residual jet energy corrections (only applied to real data)
#process.rootTupleCaloJets.ApplyResidualJEC = True
#process.rootTuplePFJets.ApplyResidualJEC = True

## Calculating rho to correct the isolation
process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.kt6PFJetsForIsolation = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

# HEEPify PAT electrons
#from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
#process.heepPatElectrons = cms.EDProducer("HEEPAttStatusToPAT",
#    eleLabel = cms.InputTag("patElectrons"),
#    barrelCuts = cms.PSet(heepBarrelCuts),
#    endcapCuts = cms.PSet(heepEndcapCuts)
#)

# Add 'heepPatElectrons' in the right place and point 'selectedLayer1Electrons' to them
#process.patDefaultSequence.replace( process.patElectrons, process.patElectrons*process.heepPatElectrons )
#process.selectedPatElectrons.src = cms.InputTag("heepPatElectrons")

# LikelihoodEle
process.load('RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi')
process.egammaIDLikelihood = process.eidLikelihoodExt.clone()

# Electron and jet cleaning deltaR parameters
process.cleanPatElectrons.checkOverlaps.muons.deltaR = 0.3
process.cleanPatJets.checkOverlaps.muons.deltaR = 0.5
process.cleanPatJets.checkOverlaps.electrons.deltaR = 0.5

# Add tau id sources:
process.patTaus.tauIDSources.leadingTrackPtCut              = cms.InputTag("shrinkingConePFTauDiscriminationByLeadingTrackPtCut")
process.patTaus.tauIDSources.trackIsolation                 = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolation")
process.patTaus.tauIDSources.trackIsolationUsingLeadingPion = cms.InputTag("shrinkingConePFTauDiscriminationByTrackIsolationUsingLeadingPion")
process.patTaus.tauIDSources.ecalIsolation                  = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolation")
process.patTaus.tauIDSources.ecalIsolationUsingLeadingPion  = cms.InputTag("shrinkingConePFTauDiscriminationByECALIsolationUsingLeadingPion")
process.patTaus.tauIDSources.byIsolation                    = cms.InputTag("shrinkingConePFTauDiscriminationByIsolation")

# Skim definition
process.load("Leptoquarks.LeptonJetFilter.leptonjetfilter_cfi")
##################################################################
#### Just pass all the events, but still store the histogram (use for MC signal events, for example)
#process.LJFilter.allEventsPassFilter = True 
#### Electron based skim
process.LJFilter.muLabel = 'muons'
process.LJFilter.elecLabel = 'gsfElectrons'
process.LJFilter.jetLabel = 'ak5CaloJets'
process.LJFilter.tauLabel = 'shrinkingConePFTauProducer'
process.LJFilter.muonsMin = -1
process.LJFilter.electronsMin = 1
process.LJFilter.elecPT = 20.
process.LJFilter.counteitherleptontype = False
##################################################################
#### Photon based skim
# process.LJFilter.muLabel = 'muons'
# process.LJFilter.elecLabel = 'gsfElectrons'
# process.LJFilter.photLabel = 'photons'
# process.LJFilter.jetLabel = 'ak5CaloJets'
# process.LJFilter.tauLabel = 'shrinkingConePFTauProducer'
# process.LJFilter.muonsMin = -1
# process.LJFilter.electronsMin = -1
# process.LJFilter.photMin = 1
# process.LJFilter.photET = 20.
# process.LJFilter.photHoE = 0.05
##################################################################

# HLT filter
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.skimHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.skimHLTFilter.HLTPaths = cms.vstring("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v*")
process.skimHLTFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

# Load HBHENoiseFilterResultProducer
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
# Check the latest recommendation from https://twiki.cern.ch/twiki/bin/view/CMS/HBHEAnomalousSignals2011
process.HBHENoiseFilterResultProducer.minRatio = cms.double(-999)
process.HBHENoiseFilterResultProducer.maxRatio = cms.double(999)
process.HBHENoiseFilterResultProducer.minHPDHits = cms.int32(17)
process.HBHENoiseFilterResultProducer.minRBXHits = cms.int32(999)
process.HBHENoiseFilterResultProducer.minHPDNoOtherHits = cms.int32(10)
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(10)
process.HBHENoiseFilterResultProducer.minHighEHitTime = cms.double(-9999.0)
process.HBHENoiseFilterResultProducer.maxHighEHitTime = cms.double(9999.0)
process.HBHENoiseFilterResultProducer.maxRBXEMF = cms.double(-999.0)
process.HBHENoiseFilterResultProducer.minNumIsolatedNoiseChannels = cms.int32(9999)
process.HBHENoiseFilterResultProducer.minIsolatedNoiseSumE = cms.double(9999)
process.HBHENoiseFilterResultProducer.minIsolatedNoiseSumEt = cms.double(9999)
process.HBHENoiseFilterResultProducer.useTS4TS5 = cms.bool(True)

#Load CosmicID producer
process.load('Leptoquarks.CosmicID.cosmicid_cfi')
process.cosmicCompatibility = cms.EDProducer("CosmicID",
                                             src=cms.InputTag("cosmicsVeto"),
                                             result = cms.string("cosmicCompatibility")
                                             )
process.timeCompatibility = process.cosmicCompatibility.clone(result = 'timeCompatibility')
process.backToBackCompatibility = process.cosmicCompatibility.clone(result = 'backToBackCompatibility')
process.overlapCompatibility = process.cosmicCompatibility.clone(result = 'overlapCompatibility')
process.patMuons.userData.userFloats.src = ['cosmicCompatibility',
                                            'timeCompatibility',
                                            'backToBackCompatibility',
                                            'overlapCompatibility']

# EcalDeadCell filters
#simpleDeltaR filter
process.load('JetMETAnalysis.simpleDRfilter.simpleDRfilter_cfi')
process.simpleDRfilter.debug = cms.untracked.bool(False)
process.simpleDRfilter.jetInputTag = cms.InputTag("ak5PFJets")
process.simpleDRfilter.metInputTag = cms.InputTag("metJESCorAK5PFJet")
#process.simpleDRfilter.metInputTag = cms.InputTag("pfMet")
process.simpleDRfilter.doFilter = cms.untracked.bool(False) # to enable filter or not
process.simpleDRfilter.makeProfileRoot = False

# Load EcalSeverityLevelESProducer (needed only if the SuperCluster module is run)
#process.load('RecoLocalCalo/EcalRecAlgos/EcalSeverityLevelESProducer_cfi')

# ChargedPFMET and related tools (from HWW analysis 2011)
#reduced PF Candidates (only charged + compatible with primary vertex)
process.load('Leptoquarks.MetTools.reducedCandidatesProducer_cfi')
process.reducedPFCands = cms.EDProducer("ReducedCandidatesProducer",
                                        srcCands = cms.InputTag('particleFlow'),
                                        srcVertices = cms.InputTag('offlinePrimaryVertices'),
                                        dz = cms.double(0.1)
                                        )
#chargedMET (using all reduced PF Candidates)
process.chargedPFMet = cms.EDProducer('ChargedPFMetProducer',
                                      collectionTag = cms.InputTag("reducedPFCands")
                                      )

# RootTupleMakerV2 tree
process.rootTupleTree = cms.EDAnalyzer("RootTupleMakerV2_Tree",
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_rootTupleEvent_*_*',
        'keep *_rootTupleEventSelection_*_*',
        'keep *_rootTupleCaloJets_*_*',
        'keep *_rootTuplePFJets_*_*',
        'keep *_rootTupleElectrons_*_*',
        'keep *_rootTupleTaus_*_*',
        'keep *_rootTupleCaloMET_*_*',
        'keep *_rootTupleTCMET_*_*',
        'keep *_rootTuplePFMET_*_*',
        'keep *_rootTuplePFMETType1Cor_*_*',
        'keep *_rootTuplePFChargedMET_*_*',
        'keep *_rootTupleMuons_*_*',
        'keep *_rootTupleTrigger_*_*',
        'keep *_rootTupleTriggerObjectsHLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilter_*_*',
        'keep *_rootTupleTriggerObjectsHLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter_*_*',
        'keep *_rootTupleTriggerObjectsHLTEle27WP80TrackIsoFilter_*_*',
        'keep *_rootTupleVertex_*_*',
        'keep *_rootTupleGenEventInfo_*_*',
        'keep *_rootTupleGenParticles_*_*',
        'keep *_rootTupleGenJets_*_*',
        'keep *_rootTupleGenMETTrue_*_*',
        'keep *_rootTupleGenMETCalo_*_*',       
        'keep *_rootTuplePhotons_*_*',
        'keep *_rootTuplePFCandidates_*_*'
    )
)

# Path definition
process.p = cms.Path(
    #process.LJFilter*
    process.skimHLTFilter*
    process.HBHENoiseFilterResultProducer*
    process.kt6PFJetsForIsolation*
    process.kt6PFJets*
    process.ak5PFJets*
    process.ak5PFJetsNoMuon*
    process.metJESCorAK5PFJet*
    process.egammaIDLikelihood*
    process.simpleDRfilter*
    (
    process.cosmicCompatibility +
    process.timeCompatibility +
    process.backToBackCompatibility +
    process.overlapCompatibility
    )*
    process.reducedPFCands*
    process.chargedPFMet*
    process.patDefaultSequence*
    (
    process.rootTupleEvent+
    process.rootTupleEventSelection+
    process.rootTupleCaloJets+
    process.rootTuplePFJets+
    process.rootTupleElectrons+
    process.rootTupleTaus+
    process.rootTupleCaloMET+
    process.rootTupleTCMET+
    process.rootTuplePFMET+
    process.rootTuplePFMETType1Cor+
    process.rootTuplePFChargedMET+
    process.rootTupleMuons+
    process.rootTupleTrigger+
    process.rootTupleTriggerObjectsHLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilter+
    process.rootTupleTriggerObjectsHLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter+
    process.rootTupleTriggerObjectsHLTEle27WP80TrackIsoFilter+
    process.rootTupleVertex+
    process.rootTupleGenEventInfo+
    process.rootTupleGenParticles+
    process.rootTupleGenJets+
    process.rootTupleGenMETTrue+
    process.rootTupleGenMETCalo+    
    process.rootTuplePhotons+
    process.rootTuplePFCandidates
    )
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
