#ifndef RootTupleMakerV2CaloJets
#define RootTupleMakerV2CaloJets

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class RootTupleMakerV2_CaloJets : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_CaloJets(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   inputTag  , inputTagL1Offset;
  const std::string     prefix,suffix;
  const unsigned int    maxSize;
  const double          electronPt, electronIso, muonPt, muonIso;
  const std::string     jecUncPath;
  const bool            readJECuncertainty;
  //OLD
  /*   const bool            applyResJEC; */
  /*   const std::string     resJEC; */
};

#endif
