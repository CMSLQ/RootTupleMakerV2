#ifndef RootTupleMakerV2PFJets
#define RootTupleMakerV2PFJets

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include <string>

class RootTupleMakerV2_PFJets : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_PFJets(const edm::ParameterSet&);

 private:
  std::string upperCase(std::string input);
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag , inputTagL1Offset;
  const edm::InputTag   inputTagSmearedUp, inputTagSmearedDown;
  const edm::InputTag   inputTagScaledUp, inputTagScaledDown;	 
  const std::string     prefix,suffix,mvaPileupIDname;
  const unsigned int    maxSize;
  const std::string     jecUncPath; 
  const bool            readJECuncertainty;
  const edm::InputTag   vtxInputTag;
  bool            isPuppiJetColl;

  //OLD
  /*   const bool            applyResJEC; */
  /*   const std::string     resJEC; */
};

#endif
