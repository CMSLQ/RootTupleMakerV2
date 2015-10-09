import FWCore.ParameterSet.Config as cms

# Remake AK5 jets from packedPFcandidates
# See: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#Advanced_topics_re_clustering_ev

## Filter out neutrinos from packed GenParticles
packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))
## Define GenJets
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
ak5GenJetsNoNu = ak5GenJets.clone(src = 'packedGenParticlesForJetsNoNu')

## Select charged hadron subtracted packed PF candidates
pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
## Define AK5PFJetsCHS
ak5PFJetsCHS = ak5PFJets.clone(src = 'pfCHS', doAreaFastjet = True)
## Define AK5PFJets
ak5PFJets = ak5PFJets.clone(src = 'packedPFCandidates')

