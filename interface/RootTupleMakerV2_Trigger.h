#ifndef RootTupleMakerV2Trigger
#define RootTupleMakerV2Trigger

#include "FWCore/Framework/interface/EDProducer.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

class RootTupleMakerV2_Trigger : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Trigger(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  void beginRun( edm::Run &, const edm::EventSetup & );
  const edm::InputTag   l1InputTag, hltInputTag;
  const std::vector<std::string> hltPathsOfInterest;
  HLTConfigProvider hltConfig;
};

#endif
