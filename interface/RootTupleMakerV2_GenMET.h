#ifndef RootTupleMakerV2GenMET
#define RootTupleMakerV2GenMET

#include "FWCore/Framework/interface/EDProducer.h"

class RootTupleMakerV2_GenMET : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_GenMET(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag;
  const std::string     prefix,suffix;
};

#endif
