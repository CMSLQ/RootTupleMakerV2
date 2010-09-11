#ifndef RootTupleMakerV2Electrons
#define RootTupleMakerV2Electrons

#include "FWCore/Framework/interface/EDProducer.h"

class RootTupleMakerV2_Electrons : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Electrons(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   trkInputTag, dcsInputTag, inputTag;
  const std::string     prefix, suffix;
  const unsigned int    maxSize;
  const double          electronIso, muonPt, muonIso;
  const std::string     muonID;
  const edm::InputTag   vtxInputTag;
};

#endif
