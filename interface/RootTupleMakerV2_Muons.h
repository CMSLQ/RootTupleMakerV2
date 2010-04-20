#ifndef RootTupleMakerV2Muons
#define RootTupleMakerV2Muons

#include "FWCore/Framework/interface/EDProducer.h"

class RootTupleMakerV2_Muons : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Muons(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag;
  const std::string     prefix,suffix;
  const unsigned int    maxSize;
  const double          muonIso;
  const std::string     muonID;
  const bool            beamSpotCorr;
};

#endif
