#ifndef RootTupleMakerV2GenJets
#define RootTupleMakerV2GenJets

#include "FWCore/Framework/interface/EDProducer.h"

class RootTupleMakerV2_GenJets : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_GenJets(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag;
  const std::string     prefix,suffix;
  const unsigned int    maxSize;
};

#endif
