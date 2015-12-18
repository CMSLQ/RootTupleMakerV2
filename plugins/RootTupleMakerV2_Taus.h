#ifndef RootTupleMakerV2Taus
#define RootTupleMakerV2Taus

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class RootTupleMakerV2_Taus : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Taus(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag, vtxInputTag;
  const std::string     prefix, suffix;
  const unsigned int    maxSize;
  const bool            isSCTau;
  const bool            isHPSTau;
};

#endif
