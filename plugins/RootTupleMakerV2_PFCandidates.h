#ifndef RootTupleMakerV2PFCandidates
#define RootTupleMakerV2PFCandidates

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

class RootTupleMakerV2_PFCandidates : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_PFCandidates(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::EDGetTokenT<pat::JetCollection>   jetInputToken_;
  const edm::EDGetTokenT<pat::ElectronCollection>   electronInputToken_;
  const edm::EDGetTokenT<pat::MuonCollection>   muonInputToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection>   pfcandInputToken_;
  const std::string     prefix, suffix;
  const unsigned int    maxSize;
  const double          DRmatch; 
};

#endif
