#ifndef RootTupleMakerV2Electrons
#define RootTupleMakerV2Electrons

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/ValueMap.h"

class RootTupleMakerV2_Electrons : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Electrons(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag trkInputTag, inputTag;
  const edm::InputTag vtxInputTag, beamSpotInputTag, triggerEventInputTag, rhoInputTag;
  const edm::InputTag electronVetoIdMapInputTag_, electronTightIdMapInputTag_;
  const edm::InputTag electronMediumIdMapInputTag_, electronLooseIdMapInputTag_;
  const edm::InputTag electronHEEPIdMapInputTag_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electronHEEPIdMapToken_;
  const double          electronIso, muonPt, muonIso;
  const std::string     muonID;
  const edm::InputTag   singleEleTriggerMatchTag, singleEleTriggerMatchWP80Tag, doubleEleTriggerMatchTag;
  const std::string     prefix, suffix;
  const unsigned int    maxSize;

  typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;
};

#endif
