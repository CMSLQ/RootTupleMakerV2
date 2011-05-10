#include "Leptoquarks/RootTupleMakerV2/interface/PatUtilities.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"

reco::TrackRef pmcTrack(const pat::Muon& mu, int & refit_id ) {
  return tevOptimized(mu.innerTrack(), mu.globalTrack(), mu.tpfmsMuon(), mu.pickyMuon(), refit_id);
}

reco::TrackRef tevOptimized(const reco::TrackRef& trackerTrack,
			    const reco::TrackRef& gmrTrack,
			    const reco::TrackRef& fmsTrack,
			    const reco::TrackRef& pmrTrack,
			    int & refit_id ) {
  const reco::TrackRef refit[4] = { 
    trackerTrack, 
    gmrTrack, 
    fmsTrack, 
    pmrTrack 
  }; 

  double prob[4] = {0.}; 
  int chosen = 3; 
 
  for (unsigned int i = 0; i < 4; ++i) 
    if (refit[i].isNonnull() && refit[i]->numberOfValidHits()) 
      prob[i] = muon::trackProbability(refit[i]); 
 
  if (prob[3] == 0.) { 
    if (prob[2] > 0.) chosen=2; 
    else if (prob[1] > 0.) chosen=1; 
    else if (prob[0] > 0.) chosen=0; 
  } 
 
  if (prob[0] > 0. && prob[3] > 0. && prob[3] - prob[0] > 30.) chosen = 0; 
  if (prob[2] > 0. && prob[chosen] - prob[2] > 0.) chosen = 2;
 
  refit_id = chosen;
  return refit[chosen]; 
} 
