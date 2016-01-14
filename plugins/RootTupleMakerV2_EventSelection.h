#ifndef RootTupleMakerV2EventSelection
#define RootTupleMakerV2EventSelection

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"

class RootTupleMakerV2_EventSelection : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_EventSelection(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  edm::EDGetTokenT<L1GlobalTriggerReadoutRecord>  l1InputToken_;
  //edm::EDGetTokenT<>   vtxInputToken_;
  //edm::EDGetTokenT<>  hcalNoiseInputToken_;
  edm::EDGetTokenT<edm::TriggerResults>  filterResultsInputToken_;

};

#endif
