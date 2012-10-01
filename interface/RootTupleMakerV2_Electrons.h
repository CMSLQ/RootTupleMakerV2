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
  const edm::InputTag   trkInputTag, dcsInputTag, inputTag;
  const edm::InputTag   vtxInputTag, beamSpotInputTag, conversionsInputTag, triggerEventInputTag, rhoInputTag;
  const std::vector <edm::InputTag> pfIsolation03InputTags, pfIsolation04InputTags;
  const double          electronIso, muonPt, muonIso;
  const std::string     muonID;
  const std::string     singleEleTriggerMatch, singleEleTriggerMatchWP80, doubleEleTriggerMatch;
  const std::string     prefix, suffix;
  const unsigned int    maxSize;

  typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;
};

#endif
