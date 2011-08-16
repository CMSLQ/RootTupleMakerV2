#ifndef RootTupleMakerV2MET
#define RootTupleMakerV2MET

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class RootTupleMakerV2_MET : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_MET(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag;
  const std::string     prefix,suffix;
  const bool            store_uncorrected_MET;
  const bool            store_MET_significance;
};

#endif
