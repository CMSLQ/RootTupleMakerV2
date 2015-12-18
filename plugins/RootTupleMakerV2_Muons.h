#ifndef RootTupleMakerV2Muons
#define RootTupleMakerV2Muons

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class RootTupleMakerV2_Muons : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Muons(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag, triggerEventInputTag;
  const std::string     prefix,suffix;
  const unsigned int    maxSize;
  const double          muonIso;
  const std::string     muonID, singleMuonTriggerMatch, singleIsoMuonTriggerMatch;
  const bool            beamSpotCorr;
  const bool            useCocktailRefits;
  const edm::InputTag   vtxInputTag;
};

#endif
