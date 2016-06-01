#ifndef RootTupleMakerV2GenJets
#define RootTupleMakerV2GenJets

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

class RootTupleMakerV2_GenJets : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_GenJets(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::EDGetTokenT<reco::GenJetCollection> genJetsInputToken_;
  const std::string     prefix,suffix;
  const unsigned int    maxSize;
};

#endif
