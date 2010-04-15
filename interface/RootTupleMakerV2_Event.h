#ifndef RootTupleMakerV2Event
#define RootTupleMakerV2Event
// DONT'USE _ IN THE NAMES ABOVE

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "TString.h"
#include <fstream>
#include <iostream>

class RootTupleMakerV2_Event : public edm::EDProducer {
 public: 
  explicit RootTupleMakerV2_Event(const edm::ParameterSet&);

 private: 
  void produce( edm::Event &, const edm::EventSetup & );
};

#endif
