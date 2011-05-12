#ifndef RootTupleMakerV2GenEventInfo
#define RootTupleMakerV2GenEventInfo

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class RootTupleMakerV2_GenEventInfo : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_GenEventInfo(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   genEvtInfoInputTag;
  const bool            storePDFWeights;
  const edm::InputTag   pdfWeightsInputTag;
  const edm::InputTag   pileupInfoSrc;
};

#endif
