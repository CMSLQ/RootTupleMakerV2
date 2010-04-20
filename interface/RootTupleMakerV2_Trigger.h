#ifndef RootTupleMakerV2Trigger
#define RootTupleMakerV2Trigger

#include "FWCore/Framework/interface/EDProducer.h"

class RootTupleMakerV2_Trigger : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Trigger(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTagL1, inputTagHLT;
  const std::vector<std::string> hltPathsOfInterest;
};

#endif
