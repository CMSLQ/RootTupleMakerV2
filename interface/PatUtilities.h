#ifndef LQ_PatUtilities_h
#define LQ_PatUtilities_h

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

reco::TrackRef pmcTrack(const pat::Muon& mu, int & refit_id );

reco::TrackRef tevOptimized(const reco::TrackRef& trackerTrack,
			    const reco::TrackRef& gmrTrack,
			    const reco::TrackRef& fmsTrack,
			    const reco::TrackRef& pmrTrack,
			    int & refit_id );

//------------------------------------------------------------------
// The code below does the work of TrackBase::validFraction()
// TrackBase::validFraction() is implemented in CMSSW_4_2_0 and later
// This function is intended for use in release prior to CMSSW_4_2_0
// When working with CMSSW_4_2_0 and later, 
//     use TrackBase::validFraction(), not this function
//------------------------------------------------------------------

template < class T > 
double validFraction ( const T & track ) { 
  
  if ( track.isNull() ) return -2 ;

  int valid   = track->hitPattern().numberOfValidTrackerHits();
  int lost    = track->hitPattern().numberOfLostTrackerHits ();
  int lostIn  = track->trackerExpectedHitsInner().numberOfLostTrackerHits();
  int lostOut = track->trackerExpectedHitsOuter().numberOfLostTrackerHits();
 
  if ((valid+lost+lostIn+lostOut)==0) return -1;
  return valid/(1.0*(valid+lost+lostIn+lostOut));
  
}

#endif 
