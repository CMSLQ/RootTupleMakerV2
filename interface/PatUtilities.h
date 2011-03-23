#ifndef LQ_PatUtilities_h
#define LQ_PatUtilities_h

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

reco::TrackRef pmcTrack(const pat::Muon& mu);

reco::TrackRef tevOptimized(const reco::TrackRef& trackerTrack,
			    const reco::TrackRef& gmrTrack,
			    const reco::TrackRef& fmsTrack,
			    const reco::TrackRef& pmrTrack);

#endif 
