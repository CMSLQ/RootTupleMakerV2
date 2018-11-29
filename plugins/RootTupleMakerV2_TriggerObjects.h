#ifndef RootTupleMakerV2TriggerObjects
#define RootTupleMakerV2TriggerObjects

#include <string> 

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

class RootTupleMakerV2_TriggerObjects : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_TriggerObjects(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::EDGetTokenT<edm::TriggerResults>   triggerBitsToken_;
  const edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection>   triggerObjectsToken_;
  const std::string     prefix, suffix;
  const std::vector<std::string> hltPathsOfInterest;

};

#endif
