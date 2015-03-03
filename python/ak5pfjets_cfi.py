import FWCore.ParameterSet.Config as cms

# Remake AK5 jets from packedPFcandidates
# See: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#Advanced_topics_re_clustering_ev

from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
from RecoMET.METProducers.PFMET_cfi import pfMet

# Select candidates that would pass CHS requirements
chs = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))

#makes chs ak5 jets   (instead of ak4 that are default in miniAOD 70X)
ak5PFJetsCHS = ak5PFJets.clone(src = 'chs')
#let's also make non-chs ak5 jets
ak5PFJets = ak5PFJets.clone(src = 'packedPFCandidates') 

# if we really want to study the ak5 vs ak4 we may need to remake gen Jets in ak5 too
# in order to do so we use the "all status 1" gen particles (packedGenParticles)
ak5GenJets = ak5GenJets.clone(src = 'packedGenParticles')

# the following part is needed if you want to run b-tagging on the freshly made jets
# CAVEAT: it is not 100% the same b-tagging as in RECO, but performance plots are almost identical

# As tracks are not stored in miniAOD, and b-tag fwk for CMSSW < 72X does not accept candidates
# we need to recreate tracks and pv for btagging in standard reco format:
from RecoBTag.Configuration.RecoBTag_cff import *
from RecoJets.Configuration.RecoJetAssociations_cff import *
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
from PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi import *
from Configuration.StandardSequences.Geometry_cff import *
from Configuration.StandardSequences.MagneticField_38T_cff import *

ak5JetTracksAssociatorAtVertexPF.jets = cms.InputTag("ak5PFJetsCHS")
ak5JetTracksAssociatorAtVertexPF.tracks = cms.InputTag("unpackedTracksAndVertices")
impactParameterTagInfos.primaryVertex = cms.InputTag("unpackedTracksAndVertices")
inclusiveSecondaryVertexFinderTagInfos.extSVCollection = cms.InputTag("unpackedTracksAndVertices","secondary","")
combinedSecondaryVertex.trackMultiplicityMin = 1
