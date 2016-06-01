#ifndef RootTupleMakerV2Vertex
#define RootTupleMakerV2Vertex

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class RootTupleMakerV2_Vertex : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Vertex(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::EDGetTokenT<reco::VertexCollection>   vtxInputToken_;
  const std::string     prefix,suffix;
};

#endif
