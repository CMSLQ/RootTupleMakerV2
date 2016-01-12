#ifndef RootTupleMakerV2Event
#define RootTupleMakerV2Event

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class RootTupleMakerV2_Event : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Event(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );

  //const std::string globalTag;
  const std::string globalTagLabel_;
  const edm::EDGetTokenT<std::string> globalTagToken_;
  const edm::EDGetTokenT<double> fixedGridRhoAllToken_;
  const edm::EDGetTokenT<double> fixedGridRhoFastjetAllToken_;
  const edm::EDGetTokenT<double> fixedGridRhoFastjetAllCaloToken_;
  const edm::EDGetTokenT<double> fixedGridRhoFastjetCentralCaloToken_;
  const edm::EDGetTokenT<double> fixedGridRhoFastjetCentralChargedPileUpToken_;
  const edm::EDGetTokenT<double> fixedGridRhoFastjetCentralNeutralToken_;

};

#endif
