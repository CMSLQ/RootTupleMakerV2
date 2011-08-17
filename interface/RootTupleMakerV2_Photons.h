#ifndef RootTupleMakerV2Photons
#define RootTupleMakerV2Photons

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class RootTupleMakerV2_Photons : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Photons(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag;
  const std::string     prefix, suffix;
  const unsigned int    maxSize;
  const edm::InputTag   beamSpotInputTag, conversionsInputTag, electronsInputTag, ecalRecHitsEBInputTag, ecalRecHitsEEInputTag;
};

#endif
