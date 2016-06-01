#ifndef RootTupleMakerV2CaloJets
#define RootTupleMakerV2CaloJets

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

class RootTupleMakerV2_CaloJets : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_CaloJets(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::EDGetTokenT<std::vector<pat::Jet> > jetInputToken_;
  //const edm::EDGetTokenT<> inputTokenL1Offset_;
  const std::string     prefix,suffix;
  const unsigned int    maxSize;
  const double          electronPt, electronIso, muonPt, muonIso;
  const std::string     jecUncPath;
  const bool            readJECuncertainty;
  //OLD
  /*   const bool            applyResJEC; */
  /*   const std::string     resJEC; */
};

#endif
