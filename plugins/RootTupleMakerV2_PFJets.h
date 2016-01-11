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
  const edm::EDGetTokenT<edm::InputTag>   inputToken_;
  const edm::InputTag   inputTag;
  const edm::EDGetTokenT<edm::InputTag>   inputTokenL1Offset_;
  const edm::EDGetTokenT<edm::InputTag>   inputTokenSmearedUp_, inputTokenSmearedDown_;
    const edm::EDGetTokenT<edm::InputTag> inputTokenScaledUp_, inputTokenScaledDown_;	
  const std::string     prefix,suffix,mvaPileupIDname;
  const unsigned int    maxSize;
  const std::string     jecUncPath; 
  const bool            readJECuncertainty;
  const edm::EDGetTokenT<edm::InputTag>   vtxInputToken_;
  bool            isPuppiJetColl;

  //OLD
  /*   const bool            applyResJEC; */
  /*   const std::string     resJEC; */
};

#endif
