import FWCore.ParameterSet.Config as cms

# ------------------------------------------------------------------------------------
# This python cfi file contains the HPS Tau ID sources for 2012 analysis.
#
# Out-of-the-box available IDs are:
#  'againstElectronLoose' 'againstElectronMVA' 'againstElectronMedium' 'againstElectronTight'
#  'againstMuonLoose' 'againstMuonMedium' 'againstMuonTight' 'byLooseCombinedIsolationDeltaBetaCorr'
#  'byMediumCombinedIsolationDeltaBetaCorr' 'byTightCombinedIsolationDeltaBetaCorr'
#  'byVLooseCombinedIsolationDeltaBetaCorr' 'decayModeFinding' .
# ------------------------------------------------------------------------------------

from PhysicsTools.PatAlgos.producersLayer1.tauProducer_cfi import *
# ------------------------------------------------------------------------------------------------------------------------------------------------ #
# Sep 25, 2012                                                                                                                                     #
# See: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatAlgos/python/producersLayer1/tauProducer_cfi.py?revision=1.32&view=markup  #
# Version comment on CVS:                                                                                                                          #
# " added new MVA based anti-e and tau isolation discriminators                                                                                    #
#   (NOTE: requires V01-04-13 RecoTauTag/RecoTau + V01-04-00 RecoTauTag/Configuration) "                                                           #
# ------------------------------------------------------------------------------------------------------------------------------------------------ #

patTaus.tauIDSources.byVLooseIsolation              = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolation")
patTaus.tauIDSources.byLooseIsolation               = cms.InputTag("hpsPFTauDiscriminationByLooseIsolation")
patTaus.tauIDSources.byMediumIsolation              = cms.InputTag("hpsPFTauDiscriminationByMediumIsolation")
patTaus.tauIDSources.byTightIsolation               = cms.InputTag("hpsPFTauDiscriminationByTightIsolation")
patTaus.tauIDSources.byVLooseIsolationDeltaBetaCorr = cms.InputTag("hpsPFTauDiscriminationByVLooseIsolationDBSumPtCorr")
patTaus.tauIDSources.byLooseIsolationDeltaBetaCorr  = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationDBSumPtCorr")
patTaus.tauIDSources.byMediumIsolationDeltaBetaCorr = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationDBSumPtCorr")
patTaus.tauIDSources.byTightIsolationDeltaBetaCorr  = cms.InputTag("hpsPFTauDiscriminationByTightIsolationDBSumPtCorr")
patTaus.tauIDSources.byIsolationMVAraw              = cms.InputTag("hpsPFTauDiscriminationByIsolationMVAraw")
patTaus.tauIDSources.byLooseIsolationMVA            = cms.InputTag("hpsPFTauDiscriminationByLooseIsolationMVA")
patTaus.tauIDSources.byMediumIsolationMVA           = cms.InputTag("hpsPFTauDiscriminationByMediumIsolationMVA")
patTaus.tauIDSources.byTightIsolationMVA            = cms.InputTag("hpsPFTauDiscriminationByTightIsolationMVA")
patTaus.tauIDSources.againstElectronMVA2raw         = cms.InputTag("hpsPFTauDiscriminationByMVA2rawElectronRejection")
patTaus.tauIDSources.againstElectronMVA2category    = cms.InputTag("hpsPFTauDiscriminationByMVA2rawElectronRejection:category")
patTaus.tauIDSources.againstElectronVLooseMVA2      = cms.InputTag("hpsPFTauDiscriminationByMVA2VLooseElectronRejection")
patTaus.tauIDSources.againstElectronLooseMVA2       = cms.InputTag("hpsPFTauDiscriminationByMVA2LooseElectronRejection")
patTaus.tauIDSources.againstElectronMediumMVA2      = cms.InputTag("hpsPFTauDiscriminationByMVA2MediumElectronRejection")
patTaus.tauIDSources.againstElectronTightMVA2       = cms.InputTag("hpsPFTauDiscriminationByMVA2TightElectronRejection")
