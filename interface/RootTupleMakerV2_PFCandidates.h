#ifndef RootTupleMakerV2PFCandidates
#define RootTupleMakerV2PFCandidates

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class RootTupleMakerV2_PFCandidates : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_PFCandidates(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   reducedPFCandidateInputTag, electronInputTag, muonInputTag;
  const std::string     prefix, suffix;
  const unsigned int    maxSize;
  const double          DRmatch; 
};

#endif
