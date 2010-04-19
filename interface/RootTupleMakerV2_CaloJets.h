#ifndef RootTupleMakerV2CaloJets
#define RootTupleMakerV2CaloJets

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "TString.h"
#include <fstream>
#include <iostream>

class RootTupleMakerV2_CaloJets : public edm::EDProducer {
 public: 
  explicit RootTupleMakerV2_CaloJets(const edm::ParameterSet&);

 private: 
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag;
  const std::string     prefix,suffix;
  const unsigned int    maxSize;
  const double          electronPt_, electronIso_, muonPt_, muonIso_;
};

#endif
