#ifndef RootTupleMakerV2GenEventInfo
#define RootTupleMakerV2GenEventInfo

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

class RootTupleMakerV2_GenEventInfo : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_GenEventInfo(const edm::ParameterSet&);


 private:
  void produce( edm::Event &, const edm::EventSetup & );
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  const edm::EDGetTokenT<GenEventInfoProduct>  genEvtInfoInputToken_;
  const bool  storePDFWeights;
  //const edm::EDGetTokenT<std::vector<double> > pdfCTEQWeightsInputToken_;
  //const edm::EDGetTokenT<std::vector<double> > pdfMMTHWeightsInputToken_;
  const edm::EDGetTokenT<std::vector<double> > pdfNNPDFWeightsInputToken_;
  //const edm::EDGetTokenT<std::vector<double> > pdfPDF4LHCWeightsInputToken_;
  const edm::EDGetTokenT<std::vector<PileupSummaryInfo> >   pileupInfoSrcToken_;
  const edm::EDGetTokenT<LHERunInfoProduct>    LHERunInfoToken_;
  const edm::EDGetTokenT<LHEEventProduct>      LHEEventProductToken_;
};

#endif
