import FWCore.ParameterSet.Config as cms

# ------------------------------------------------------------------------------------
# This python cfi file contains the recommended MET filters for 2012 analysis.
# See: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters
# ------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------
# CSC beam halo filter: 
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters#CSC_Beam_Halo_Filter
# ------------------------------------------------------------------------------------
from RecoMET.METAnalyzers.CSCHaloFilter_cfi import *

# ------------------------------------------------------------------------------------
# HBHE noise filter:
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters#HBHE_Noise_Filter
# ------------------------------------------------------------------------------------

from CommonTools.RecoAlgos.HBHENoiseFilter_cfi import *
HBHENoiseFilter.taggingMode = cms.bool ( True ) 

# ------------------------------------------------------------------------------------
# HBHE noise filter results producer
# https://twiki.cern.ch/twiki/bin/view/CMS/HBHEAnomalousSignals2012
# ------------------------------------------------------------------------------------

from CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi import *

# HBHENoiseFilterResultProducer = cms.EDProducer( 'HBHENoiseFilterResultProducer',
#                                                 noiselabel = cms.InputTag('hcalnoise','','RECO'),
#                                                 minRatio = cms.double(-999),
#                                                 maxRatio = cms.double(999),
#                                                 minHPDHits = cms.int32(17),
#                                                 minRBXHits = cms.int32(999),
#                                                 minHPDNoOtherHits = cms.int32(10),
#                                                 minZeros = cms.int32(10),
#                                                 minHighEHitTime = cms.double(-9999.0),
#                                                 maxHighEHitTime = cms.double(9999.0),
#                                                 maxRBXEMF = cms.double(-999.0),
#                                                 minNumIsolatedNoiseChannels = cms.int32(9999),
#                                                 minIsolatedNoiseSumE = cms.double(9999),
#                                                 minIsolatedNoiseSumEt = cms.double(9999),
#                                                 useTS4TS5 = cms.bool(True)
#                                                 )


# ------------------------------------------------------------------------------------
# HCAL laser filter: 
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters#HCAL_laser_events_updated
# ------------------------------------------------------------------------------------

from Leptoquarks.MetTools.hcallasereventfilter2012_lq_cfi import * 
hcallasereventfilter2012_lq.eventFileName = cms.FileInPath('EventFilter/HcalRawToDigi/data/HCALLaser2012AllDatasets.txt.gz')

# ------------------------------------------------------------------------------------
# ECAL dead cell filter:
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters#ECAL_dead_cell_filter
# ------------------------------------------------------------------------------------
from RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi import *
EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)
# The section below is for the filter on Boundary Energy. Available in AOD in CMSSW>44x
# For releases earlier than 44x, one should make the following changes
# process.EcalDeadCellBoundaryEnergyFilter.recHitsEB = cms.InputTag("ecalRecHit","EcalRecHitsEB")
# process.EcalDeadCellBoundaryEnergyFilter.recHitsEE = cms.InputTag("ecalRecHit","EcalRecHitsEE")
from RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi import * 
EcalDeadCellBoundaryEnergyFilter.taggingMode = cms.bool(True)
EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyDeadCellsEB=cms.untracked.double(10)
EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyDeadCellsEE=cms.untracked.double(10)
EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyGapEB=cms.untracked.double(100)
EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyGapEE=cms.untracked.double(100)
EcalDeadCellBoundaryEnergyFilter.enableGap=cms.untracked.bool(False)
EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEB = cms.vint32(12,14)
EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEE = cms.vint32(12,14)
# End of Boundary Energy filter configuration 

# ------------------------------------------------------------------------------------
# Tracking failure filter:
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters#Tracking_failure_filter_updated
# ------------------------------------------------------------------------------------
from RecoMET.METFilters.trackingFailureFilter_cfi import *
goodVertices = cms.EDFilter(
  "VertexSelector",
  filter = cms.bool(False),
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)
trackingFailureFilter.taggingMode = cms.bool (True) 

# ------------------------------------------------------------------------------------
# Bad EE Supercrystal filter: 
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters#Bad_EE_Supercrystal_filter_added
# ------------------------------------------------------------------------------------
from RecoMET.METFilters.eeBadScFilter_cfi import *
eeBadScFilter.taggingMode = cms.bool (True)

# ------------------------------------------------------------------------------------
# Seema Sharma's bad (large laser) ECAL correction filter:
# https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters#EB_or_EE_Xtals_with_large_laser
# ------------------------------------------------------------------------------------

from RecoMET.METFilters.ecalLaserCorrFilter_cfi import *
ecalLaserCorrFilter.taggingMode = cms.bool (True)

# ------------------------------------------------------------------------------------
# Tracking POG filters in tagging mode
# ------------------------------------------------------------------------------------

from RecoMET.METFilters.trackingPOGFilters_cff import *

manystripclus53X.taggedMode         = cms.untracked.bool(True )
manystripclus53X.forcedValue        = cms.untracked.bool(False)

toomanystripclus53X.taggedMode      = cms.untracked.bool(True )
toomanystripclus53X.forcedValue     = cms.untracked.bool(False)

logErrorTooManyClusters.taggedMode  = cms.untracked.bool(True )
logErrorTooManyClusters.forcedValue = cms.untracked.bool(False)
