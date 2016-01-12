#ifndef RootTupleMakerV2PhysicsDSTStream2011
#define RootTupleMakerV2PhysicsDSTStream2011

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

class RootTupleMakerV2_PhysicsDSTStream2011 : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_PhysicsDSTStream2011(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  void beginRun( edm::Run &, const edm::EventSetup & );
  void printNames(const std::vector<std::string>& names);

  HLTConfigProvider hltConfig;

  const edm::InputTag       hltInputTag;
  const edm::EDGetTokenT<edm::TriggerResults>     hltInputToken_;
  const edm::EDGetTokenT<reco::PFJetCollection>   inputTokenHLTPFJets_;
  const edm::EDGetTokenT<reco::CaloJetCollection> inputTokenHLTCaloJetsRaw_;
  const edm::EDGetTokenT<reco::CaloJetCollection> inputTokenHLTCaloJetsCorr_;
  const edm::EDGetTokenT<reco::VertexCollection>  inputTokenHLTPixelVertices_;
  const std::string           suffix;
  const std::string           prefixHLTPFJets, prefixHLTCaloJetsRaw, prefixHLTCaloJetsCorr, prefixHLTPixelVertices;
  const double                minPtHLTPFJets, maxEtaHLTPFJets;
  const double                minPtHLTCaloJetsRaw, maxEtaHLTCaloJetsRaw;
  const double                minPtHLTCaloJetsCorr, maxEtaHLTCaloJetsCorr;  
  const bool                  storeEventInfo;

};

#endif
