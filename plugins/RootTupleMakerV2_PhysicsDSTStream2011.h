#ifndef RootTupleMakerV2PhysicsDSTStream2011
#define RootTupleMakerV2PhysicsDSTStream2011

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

class RootTupleMakerV2_PhysicsDSTStream2011 : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_PhysicsDSTStream2011(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  void beginRun( edm::Run &, const edm::EventSetup & );
  void printNames(const std::vector<std::string>& names);

  HLTConfigProvider hltConfig;

  const edm::InputTag         hltInputTag;
  const edm::InputTag         inputTagHLTPFJets, inputTagHLTCaloJetsRaw, inputTagHLTCaloJetsCorr, inputTagHLTPixelVertices;
  const std::string           suffix;
  const std::string           prefixHLTPFJets, prefixHLTCaloJetsRaw, prefixHLTCaloJetsCorr, prefixHLTPixelVertices;
  const double                minPtHLTPFJets, maxEtaHLTPFJets;
  const double                minPtHLTCaloJetsRaw, maxEtaHLTCaloJetsRaw;
  const double                minPtHLTCaloJetsCorr, maxEtaHLTCaloJetsCorr;  
  const bool                  storeEventInfo;

};

#endif
