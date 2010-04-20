#ifndef RootTupleMakerV2CaloJets
#define RootTupleMakerV2CaloJets

#include "FWCore/Framework/interface/EDProducer.h"

class RootTupleMakerV2_CaloJets : public edm::EDProducer {
 public: 
  explicit RootTupleMakerV2_CaloJets(const edm::ParameterSet&);

 private: 
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag;
  const std::string     prefix,suffix;
  const unsigned int    maxSize;
  const double          electronPt, electronIso, muonPt, muonIso;
};

#endif
