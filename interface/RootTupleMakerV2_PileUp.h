#ifndef RootTupleMakerV2PileUp
#define RootTupleMakerV2PileUp

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class RootTupleMakerV2_PileUp : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_PileUp(const edm::ParameterSet& pset);
 private:
   void produce( edm::Event &, const edm::EventSetup & );
   edm::InputTag pileupInfoSrc;
};
#endif

