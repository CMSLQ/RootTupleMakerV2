#ifndef RootTupleMakerV2MET
#define RootTupleMakerV2MET

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "TString.h"
#include <fstream>
#include <iostream>

class RootTupleMakerV2_MET : public edm::EDProducer {
 public: 
  explicit RootTupleMakerV2_MET(const edm::ParameterSet&);

 private: 
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag;
  const std::string     prefix,suffix;
};

#endif
