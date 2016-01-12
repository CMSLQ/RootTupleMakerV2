#ifndef RootTupleMakerV2Taus
#define RootTupleMakerV2Taus

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

class RootTupleMakerV2_Taus : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Taus(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::EDGetTokenT<std::vector<pat::Tau> >  tauInputToken_;
  const edm::EDGetTokenT<reco::VertexCollection>  vtxInputToken_;
  const std::string     prefix, suffix;
  const unsigned int    maxSize;
  const bool            isSCTau;
  const bool            isHPSTau;
};

#endif
