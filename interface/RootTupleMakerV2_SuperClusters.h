#ifndef RootTupleMakerV2SuperClusters
#define RootTupleMakerV2SuperClusters

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class RootTupleMakerV2_SuperClusters : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_SuperClusters(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   ebInputTag,eeInputTag;
  const std::string     prefix,suffix;
  const unsigned int    ebMaxSize, eeMaxSize;
  const edm::InputTag   ecalEBInputTag,ecalEEInputTag,trkInputTag,eleInputTag;
  const std::string     elePrefix;
  const unsigned int    eleMaxSize;
};

#endif
