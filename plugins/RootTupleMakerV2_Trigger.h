#ifndef RootTupleMakerV2Trigger
#define RootTupleMakerV2Trigger

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

class RootTupleMakerV2_Trigger : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Trigger(const edm::ParameterSet&);
  ~RootTupleMakerV2_Trigger() {};

 private:
  virtual void beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup) override;
  virtual void beginJob() {};
  virtual void produce( edm::Event &, const edm::EventSetup & ) override;
  virtual void endJob() {};

  void getDataSource() ;
  void printNames(const std::vector<std::string>& names);

  enum DataSource { NOT_APPLICABLE, STREAM, DATASET };
  const edm::InputTag   l1InputTag_;
  const edm::InputTag   hltInputTag_;
  const edm::EDGetTokenT<L1GlobalTriggerReadoutRecord>   l1InputToken_;
  const edm::EDGetTokenT<edm::TriggerResults>   hltInputToken_;
  const std::vector<std::string> hltPathsOfInterest;
  HLTPrescaleProvider hltPrescaleProvider_;

  std::string                 sourceName;
  DataSource                  sourceType;
  std::vector<std::string>    dataSource;

};

#endif
