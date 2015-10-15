import os
#----------------------------------------------------------------------------------------------------
# Load PAT template + customize
#----------------------------------------------------------------------------------------------------
# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

#process.load('PhysicsTools.PatAlgos.patSequences_cff')
#process.load('Configuration.StandardSequences.Services_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff')

# Change process name
process._Process__name="ROOTTUPLEMAKERV2"
# Options and Output Report
process.options.wantSummary = False
process.options.allowUnscheduled = cms.untracked.bool(True)

############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.default.limit = 10
#################################################################

# We should be using PFIso by default in newer CMSSW
# see: PhysicsTools/PatAlgos/python/producersLayer1/electronProducer_cff.py

#----------------------------------------------------------------------------------------------------
# Load our RootTupleMakerV2 modules
#----------------------------------------------------------------------------------------------------

process.load('Leptoquarks.RootTupleMakerV2.Ntuple_cff')

# Output ROOT file
process.TFileService = cms.Service("TFileService",
    fileName = cms.string( "file_m650.root" )
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
process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v2'
# feed it into the ntuple
process.rootTupleEvent.globalTag = process.GlobalTag.globaltag

# Events to process
process.maxEvents.input = 10

# Input files
process.source.fileNames = [
    # specified by InputList.txt
    # Here is a DYJetsToLL file
    '/store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/D461E72B-306D-E511-9457-90B11C04FE38.root'
    ]

# SIC Replace with HEEP 5.1/6.0
# Also load Egamma cut-based ID in new VID framework while we're at it
# See: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2
#   and for HEEP: https://hypernews.cern.ch/HyperNews/CMS/get/egamma/1519/2/1/1/1.html
# Set up everything that is needed to compute electron IDs and
# add the ValueMaps with ID decisions into the event data stream
#
# Load tools and function definitions
#from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
process.load("RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cfi")
#process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('slimmedElectrons')
#from PhysicsTools.SelectorUtils.centralIDRegistry import central_id_registry
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
#process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)
# Define which IDs we want to produce
# Each of these two example IDs contains all four standard
# cut-based ID working points
my_id_modules = []
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff')
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff')
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff') # for 50 ns, 13 TeV data
my_id_modules.append('RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff')
#Add them to the VID producer
for idmod in my_id_modules:
  setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
# XXX NB, must be the same as input collection used for electron ntuplizer
#process.egmGsfElectronIDs.physicsObjectSrc = process.rootTupleElectrons.InputTag

#----------------------------------------------------------------------------------------------------
# Turn on trigger matching
#----------------------------------------------------------------------------------------------------
# stored in pat::TriggerObjectStandAlone with default MiniAOD config
# See: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
# To access, see: https://hypernews.cern.ch/HyperNews/CMS/get/physTools/3286/1/2.html
# This is based off of the hypernews post
# FIXME check that we match to the paths we want while we're updating the cpp code
from PhysicsTools.PatAlgos.slimming.unpackedPatTrigger_cfi import unpackedPatTrigger
# need to load it into the process, or else it won't run
process.unpackedPatTrigger = unpackedPatTrigger.clone()

process.cleanElectronTriggerMatchHLTSingleElectronWP85.matched = 'unpackedPatTrigger'
process.cleanElectronTriggerMatchHLTDoubleElectron.matched = 'unpackedPatTrigger'
process.cleanElectronTriggerMatchHLTElectronJetJet.matched = 'unpackedPatTrigger'
process.cleanMuonTriggerMatchHLTSingleMuon.matched = 'unpackedPatTrigger'
process.cleanMuonTriggerMatchHLTSingleIsoMuon.matched = 'unpackedPatTrigger'
process.pfjetTriggerMatchHLTEleJetJet.matched = 'unpackedPatTrigger'

#process.out.outputCommands += ['keep *_electronsTriggerMatchHLTSingleElectron_*_*',
#                               'keep *_electronsTriggerMatchHLTSingleElectronWP80 _*_*',
#                               'keep *_electronsTriggerMatchHLTDoubleElectron_*_*'      ]
## select the matching electrons and keep them
#process.out.outputCommands += ['keep *_ electronsTriggeredHLTSingleElectron_*_*',
#                               'keep *_ electronsTriggeredHLTSingleElectronWP80_*_*',
#                               'keep *_ electronsTriggeredHLTDoubleElectron_*_*' ]

#----------------------------------------------------------------------------------------------------
# Add PFMET and TCMET
# See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#MET_Tools
#----------------------------------------------------------------------------------------------------
# SIC: I think we can get away with using the default MiniAOD MET
# in any case, below does not work in CMSSW_7_4_14 with 2015D PromptReco-v3
## Reproduce "raw" MET from packedPFCandidates
#from RecoMET.METProducers.PFMET_cfi import pfMet
#process.pfMet = pfMet.clone(src = "packedPFCandidates")
#process.pfMet.calculateSignificance = False # this can't be easily implemented on packed PF candidates at the moment (as of Feb 17 2015)

#----------------------------------------------------------------------------------------------------
# MET filters
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
#----------------------------------------------------------------------------------------------------
# SIC: A number of filters are run by default:
# FIXME NEEDS UPDATE
#Flag_HBHENoiseFilter = cms.Path(HBHENoiseFilter)
#Flag_CSCTightHaloFilter = cms.Path(CSCTightHaloFilter)
#Flag_hcalLaserEventFilter = cms.Path(hcalLaserEventFilter)
#Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(EcalDeadCellTriggerPrimitiveFilter)
#Flag_goodVertices = cms.Path(goodVertices)
#Flag_trackingFailureFilter = cms.Path(goodVertices + trackingFailureFilter)
#Flag_eeBadScFilter = cms.Path(eeBadScFilter)
#Flag_ecalLaserCorrFilter = cms.Path(ecalLaserCorrFilter)
#Flag_trkPOGFilters = cms.Path(trkPOGFilters)
#Flag_trkPOG_manystripclus53X = cms.Path(~manystripclus53X)
#Flag_trkPOG_toomanystripclus53X = cms.Path(~toomanystripclus53X)
#Flag_trkPOG_logErrorTooManyClusters = cms.Path(~logErrorTooManyClusters)
# See: https://github.com/cms-sw/cmssw/blob/CMSSW_7_2_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
# They are stored in edm::TriggerResults of the PAT process, and can be checked the same way as HLT paths.


#----------------------------------------------------------------------------------------------------
# Use default MiniAOD Taus
# Use default MiniAOD Tau IDs
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Make analysisPatTaus and add them to the cleanPatCandidates sequence
#----------------------------------------------------------------------------------------------------
#process.analysisPatTaus = process.cleanPatTaus.clone()
#process.analysisPatTaus.preselection = cms.string(
#    'tauID("decayModeFinding") > 0.5 &'
#    ' tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 &'
#    ' tauID("againstMuonLoose3") > 0.5 &'
#    ' tauID("againstElectronLooseMVA3") > 0.5'
#)
#process.analysisPatTaus.finalCut = cms.string('pt > 20. & abs(eta) < 2.3')
#process.cleanPatCandidates.replace ( process.cleanPatTaus, process.cleanPatTaus + process.analysisPatTaus )
# FIXME

#----------------------------------------------------------------------------------------------------
# Make analysisPatMuons and add them to the cleanPatCandidates sequence
#----------------------------------------------------------------------------------------------------
#process.analysisPatMuons = process.cleanPatMuons.clone()
#process.analysisPatMuons.finalCut = cms.string("isGlobalMuon & muonID('GlobalMuonPromptTight') & pt > 20")
#process.cleanPatCandidates.replace ( process.cleanPatMuons, process.cleanPatMuons + process.analysisPatMuons )
#FIXME

#----------------------------------------------------------------------------------------------------
# Make analysisPatElectrons and add them to the cleanPatCandidates sequence
#----------------------------------------------------------------------------------------------------
#process.analysisPatElectrons = process.cleanPatElectrons.clone()
#process.analysisPatElectrons.finalCut = cms.string('userInt("HEEPId") < 0.5')
#process.cleanPatCandidates.replace ( process.cleanPatElectrons, process.cleanPatElectrons + process.analysisPatElectrons )
# FIXME

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
# FIXME: Do we need to rebuild AK5 jets? Or can we stick with AK4?
# AK4 will possibly have more support (for JECs, etc.)
process.load('Leptoquarks.RootTupleMakerV2.ak5pfjets_cfi')
## b-tag discriminators
bTagDiscriminators = [
    'pfCombinedInclusiveSecondaryVertexV2BJetTags'
]
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
addJetCollection(process,
                 labelName = 'AK5PF',
                 jetSource = cms.InputTag('ak5PFJets'),
                 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
                 pfCandidates = cms.InputTag('packedPFCandidates'),
                 svSource = cms.InputTag('slimmedSecondaryVertices'),
                 btagDiscriminators = bTagDiscriminators,
                 jetCorrections = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
                 genJetCollection = cms.InputTag('ak5GenJetsNoNu'),
                 genParticles = cms.InputTag('prunedGenParticles'),
                 algo = 'AK',
                 rParam = 0.5,
)
addJetCollection(process,
                 labelName = 'AK5PFCHS',
                 jetSource = cms.InputTag('ak5PFJetsCHS'),
                 pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
                 pfCandidates = cms.InputTag('packedPFCandidates'),
                 svSource = cms.InputTag('slimmedSecondaryVertices'),
                 btagDiscriminators = bTagDiscriminators,
                 jetCorrections = ('AK5PF', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
                 genJetCollection = cms.InputTag('ak5GenJetsNoNu'),
                 genParticles = cms.InputTag('prunedGenParticles'),
                 algo = 'AK',
                 rParam = 0.5,
)
#process.patJetsAK5.userData.userFloats.src = [] # start with empty list of user floats
process.selectedPatJetsAK5PFCHS.cut = cms.string("pt > 10")
#process.patJetGenJetMatchAK5.matched =  'slimmedGenJets'
#process.patJetPartonMatchAK5.matched = 'prunedGenParticles'
#process.patJetPartons.particles = 'prunedGenParticles'
#process.patJetPartonsLegacy.src = 'prunedGenParticles'
#process.patJetCorrFactorsAK5.primaryVertices = 'offlineSlimmedPrimaryVertices'
## needed when using btagInfos?
##process.jetTracksAssociatorAtVertexAK5.tracks = 'unpackedTracksAndVertices'
##process.jetTracksAssociatorAtVertexAK5.pvSrc = 'offlineSlimmedPrimaryVertices'
from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))

#process.out.outputCommands += ['keep *_ak5GenJets_*_EX',
#                               'keep *_ak5PFJets_*_EX',
#                               'keep *_ak5PFJetsCHS_*_EX', ]

#----------------------------------------------------------------------------------------------------
# Make analysisPatJets and add them to the patDefaultSequence
#----------------------------------------------------------------------------------------------------
#process.cleanPatJetsAK5PF = process.cleanPatJets.clone()
#process.cleanPatJetsAK5PF.src = cms.InputTag('patJetsAK5PF')
#process.patDefaultSequence.replace (process.cleanPatJets, process.cleanPatJets + process.cleanPatJetsAK5PF)
## FIXME
#process.analysisPatJetsAK5PF = process.cleanPatJetsAK5PF.clone()
#process.analysisPatJetsAK5PF.finalCut = cms.string("abs(eta)<2.5 & pt > 20")
#process.patDefaultSequence.replace ( process.cleanPatJetsAK5PF, process.cleanPatJetsAK5PF + process.analysisPatJetsAK5PF )
## FIXME

#----------------------------------------------------------------------------------------------------
# Add the pileup MVA to the PFJets
#----------------------------------------------------------------------------------------------------
# already included in non-puppi jets in MiniAOD V2

#----------------------------------------------------------------------------------------------------
# No CaloJets in MiniAOD
#----------------------------------------------------------------------------------------------------

##----------------------------------------------------------------------------------------------------
## Define the systematic shift correction
##----------------------------------------------------------------------------------------------------
#process.load("JetMETCorrections.Type1MET.correctedMet_cff")
#process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")
#process.load("JetMETCorrections.Type1MET.correctionTermsPfMetShiftXY_cff")
## FIXME needs 72X update when available
#process.pfMEtSysShiftCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCDvsNvtx_data
## Type1 PFMET provided with MiniAOD default
## FIXME Check implementation of this -- no longer needed?

##----------------------------------------------------------------------------------------------------
## Use the runMetUncertainties tool for AK5PFJets
## https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription
## https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATTools#METSysTools
##----------------------------------------------------------------------------------------------------
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMETCorrectionsAndUncertainties
addJetCollection(process, postfix   = "ForMetUnc", labelName = 'AK5PF', jetSource = cms.InputTag('ak5PFJets'), jetCorrections = ('AK5PF', ['L1FastJet', 'L2Relative', 'L3Absolute'], ''))

# FIXME: this seems to rerun the PFMET producer, which requires RECO
#process.patJetsAK5PFForMetUnc.getJetMCFlavour = False
#runMETCorrectionsAndUncertainties(process,
#                                  metType='PF',
#                                  correctionLevel=["T0","T1","T2","Txy","Smear",""],
#                                  electronCollection="slimmedElectrons",
#                                  photonCollection="slimmedPhotons",
#                                  muonCollection="slimmedMuons",
#                                  tauCollection="slimmedTaus",
#                                  jetCollection="selectedPatJetsAK5PFForMetUnc",
#                                  jetCollectionUnskimmed="ak5PFJets",
#                                  addToPatDefaultSequence=False,
#                                  onMiniAOD=True,
#                                  runOnData=True,
#                                  postfix='')
#process.patJetPartonMatchAK5PFForMetUnc.matched = 'prunedGenParticles'
#process.patJetPartonsForMetUnc.particles = 'prunedGenParticles'
#process.patJetPartonsLegacyForMetUnc.src = 'prunedGenParticles'
#process.patJetCorrFactorsAK5PFForMetUnc.primaryVertices = 'offlineSlimmedPrimaryVertices'

#----------------------------------------------------------------------------------------------------
# Available pat::MET collections for analysis
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#ETmiss
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# This is MC, so analyze the smeared PFJets by default
# But smearing should be done by default in MiniAOD ?
#----------------------------------------------------------------------------------------------------
# SIC FIXME TODO
# analyze the default jet collection in MiniAOD (AK4PFJetsCHS)
#process.rootTuplePFJets.InputTag = cms.InputTag('slimmedJets')
# analyze the remade AK5 PFJets (we also have the AK5PFJetsCHS available)
#process.rootTuplePFJets.InputTag = 'patJetsAK5'
#process.rootTuplePFJets.InputTag = 'pfjetTriggerMatchEmbedderHLTEleJetJet'
#process.pfjetTriggerMatchHLTEleJetJet.src = 'patJetsAK5' # use for trigger matching
#process.pfjetTriggerMatchEmbedderHLTEleJetJet.src = 'patJetsAK5' # also use for trigger match embedding

#process.rootTuplePFJets.InputTag = 'pfjetTriggerMatchEmbedderHLTEleJetJet'
#process.rootTuplePFJets.InputTag = cms.InputTag('smearedAnalysisPatJetsAK5PF')
#process.rootTuplePFJets.InputTagSmearedUp   = cms.InputTag('smearedAnalysisPatJetsAK5PFresUp')                                 
#process.rootTuplePFJets.InputTagSmearedDown = cms.InputTag('smearedAnalysisPatJetsAK5PFresDown')                                 
#process.rootTuplePFJets.InputTagScaledUp    = cms.InputTag('shiftedAnalysisPatJetsAK5PFenUpForCorrMEt')                                 
#process.rootTuplePFJets.InputTagScaledDown  = cms.InputTag('shiftedAnalysisPatJetsAK5PFenDownForCorrMEt')     

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
process.LJFilter.muonsMin = 0
process.LJFilter.muPT     = 10.0
process.LJFilter.electronsMin = 0
process.LJFilter.elecPT       = 15.0
process.LJFilter.tausMin = 0
process.LJFilter.tauPT   = 15.0
process.LJFilter.jetsMin = 0
process.LJFilter.jetPT   = 15.0
process.LJFilter.counteitherleptontype = True
process.LJFilter.customfilterEMuTauJet2012 = True
# -- WARNING :
# "customfilterEMuTauJet2012" configuration is hard-coded.
# (see: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Leptoquarks/LeptonJetFilter/src/LeptonJetFilter.cc )
# "customfilterEMuTauJet2012" is the desired mode of operation for the Lepton+Jets Filter in 2012.

#----------------------------------------------------------------------------------------------------
# PDF weights
#----------------------------------------------------------------------------------------------------

#process.pdfWeights = cms.EDProducer("PdfWeightProducer",
#	# Fix POWHEG if buggy (this PDF set will also appear on output,
#	# so only two more PDF sets can be added in PdfSetNames if not "")
#	#FixPOWHEG = cms.untracked.string("CT10.LHgrid"),
#       # GenTag = cms.untracked.InputTag("genParticles"),
#	PdfInfoTag = cms.untracked.InputTag("generator"),
#	PdfSetNames = cms.untracked.vstring(
#            "CT10.LHgrid" , 
#            "MSTW2008nlo68cl.LHgrid",
#            #"NNPDF20_100.LHgrid"#DMM FIXME couldn't find this file in /cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_4/src/PhysicsTools/PatUtils/data/Summer13_V1_DATA_UncertaintySources_AK5PF.txt
#	)
#)

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
        # Trigger objects
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
    # L+J Filter
    process.LJFilter*  
    # Put everything into the tree
    # In unscheduled mode, anything 'kept' in the output commands above
    #  will have its producer module called automatically
    process.rootTupleTree
)


#----------------------------------------------------------------------------------------------------
# Dump if necessary
#----------------------------------------------------------------------------------------------------

#process.dump = cms.OutputModule("PoolOutputModule",
#                                outputCommands = cms.untracked.vstring(
#                                'keep *',
#                                ),
#                                fileName = cms.untracked.string('dump.root')
#                                )
#process.DUMP    = cms.EndPath (process.dump)

#----------------------------------------------------------------------------------------------------
# Run the path
#----------------------------------------------------------------------------------------------------
# Delete predefined Endpath (needed for running with CRAB)
del process.out
del process.outpath

process.schedule = cms.Schedule(process.p)#,process.DUMP)

#print process.dumpPython()

