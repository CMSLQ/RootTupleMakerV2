import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.mcMatchLayer0.electronMatch_cfi import *
from PhysicsTools.PatAlgos.mcMatchLayer0.muonMatch_cfi import *
from PhysicsTools.PatAlgos.mcMatchLayer0.tauMatch_cfi import *

#----------------------------------------------------------------------------------------------------
# Lepton-Gen Matching 
#
# Example is provided here:  https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATMCMatchingExercise
#
# Access (Electrons - Muons - Taus):
# for(uint i = 0 ; i < it->genParticleRefs().size() ; ++i ){
#  it->genParticle(i)->status();
#  it->genParticle(i)->pt();
# }
#
# Access (Tau Jets):
# if( it->genJet() ) it->genJet()->pt();
#
#----------------------------------------------------------------------------------------------------

# Electron - Gen Particle Matching
elMatch = electronMatch.clone(mcStatus = [3])
elMatch.maxDeltaR = cms.double(0.5)
elMatch.maxDPtRel = cms.double(999.9)
elMatch.resolveAmbiguities = cms.bool(True)
elMatch.resolveByMatchQuality = cms.bool(True)

# Muon - Gen Particle Matching
muMatch = muonMatch.clone(mcStatus = [3])
muMatch.maxDeltaR = cms.double(0.5)
muMatch.maxDPtRel = cms.double(999.9)
muMatch.resolveAmbiguities = cms.bool(True)
muMatch.resolveByMatchQuality = cms.bool(True)

# Tau - Gen Particle Matching
tauLepMatch = tauMatch.clone(mcStatus = [2])
tauLepMatch.maxDeltaR = cms.double(0.7)
tauLepMatch.maxDPtRel = cms.double(999.9)
tauLepMatch.resolveAmbiguities = cms.bool(True)
tauLepMatch.resolveByMatchQuality = cms.bool(True)
# Tau - Gen Jet Matching
tauJetMatch = tauGenJetMatch.clone()
tauJetMatch.maxDeltaR = cms.double(0.7)
tauJetMatch.maxDPtRel = cms.double(999.9)
tauJetMatch.resolveAmbiguities = cms.bool(True)
tauJetMatch.resolveByMatchQuality = cms.bool(True)
