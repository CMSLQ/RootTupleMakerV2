#ifndef RootTupleMakerV2Vertex
#define RootTupleMakerV2Vertex

#include "FWCore/Framework/interface/EDProducer.h"

class RootTupleMakerV2_Vertex : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Vertex(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag;
  const std::string     prefix,suffix;
};

#endif
