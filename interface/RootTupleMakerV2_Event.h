#ifndef RootTupleMakerV2Event
#define RootTupleMakerV2Event

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class RootTupleMakerV2_Event : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Event(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
};

#endif
