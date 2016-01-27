#ifndef RootTupleMakerV2EventSelection
#define RootTupleMakerV2EventSelection

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class RootTupleMakerV2_EventSelection : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_EventSelection(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   l1InputTag,vtxInputTag;
  const edm::InputTag   hbheNoiseFilterResultInputTag;
  const edm::InputTag   hbheIsoNoiseFilterResultInputTag;
  const edm::InputTag   filterResultsInputTag;

};

#endif
