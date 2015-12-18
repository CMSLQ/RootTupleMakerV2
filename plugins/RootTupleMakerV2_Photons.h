#ifndef RootTupleMakerV2Photons
#define RootTupleMakerV2Photons

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

class RootTupleMakerV2_Photons : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Photons(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  inline float recHitE( const  DetId id,  const EcalRecHitCollection &recHits );
  inline float recHitE( const DetId id, const EcalRecHitCollection & recHits,int di, int dj );
  inline float recHitApproxEt(  const DetId id,  const EcalRecHitCollection &recHits );
  float GetE2OverE9( const DetId id, const EcalRecHitCollection & recHits);
  const edm::InputTag   inputTag;
  const std::string     prefix, suffix;
  const unsigned int    maxSize;
  const edm::InputTag   beamSpotInputTag, conversionsInputTag, electronsInputTag, ecalRecHitsEBInputTag, ecalRecHitsEEInputTag;
};

#endif
