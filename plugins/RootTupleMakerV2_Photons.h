#ifndef RootTupleMakerV2Photons
#define RootTupleMakerV2Photons

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

class RootTupleMakerV2_Photons : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Photons(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  inline float recHitE( const  DetId id,  const EcalRecHitCollection &recHits );
  inline float recHitE( const DetId id, const EcalRecHitCollection & recHits,int di, int dj );
  inline float recHitApproxEt(  const DetId id,  const EcalRecHitCollection &recHits );
  float GetE2OverE9( const DetId id, const EcalRecHitCollection & recHits);
  const edm::EDGetTokenT<std::vector<pat::Photon> >   photonInputToken_;
  const std::string     prefix, suffix;
  const unsigned int    maxSize;
  const edm::EDGetTokenT<reco::BeamSpot>   beamSpotInputToken_;
  const edm::EDGetTokenT<reco::ConversionCollection>   conversionsInputToken_;
  const edm::EDGetTokenT<reco::GsfElectronCollection>   electronsInputToken_;
  const edm::EDGetTokenT<EBRecHitCollection>   ecalRecHitsEBInputToken_;
  const edm::EDGetTokenT<EERecHitCollection>   ecalRecHitsEEInputToken_;
};

#endif
