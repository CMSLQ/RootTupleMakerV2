#ifndef RootTupleMakerV2MET
#define RootTupleMakerV2MET

#include "FWCore/Framework/interface/EDProducer.h"

class RootTupleMakerV2_MET : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_MET(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag;
  const std::string     prefix,suffix;
};

#endif
