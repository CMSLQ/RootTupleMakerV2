#ifndef RootTupleMakerV2PFJets
#define RootTupleMakerV2PFJets

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
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
  const bool            readJERuncertainty;
  const edm::InputTag   vtxInputTag;
  bool            isPuppiJetColl;
  const std::string jerUncPath;
  const bool jer_from_gt;
  const edm::FileInPath jer_resolutions_file;
  const edm::FileInPath jer_scale_factors_file;
  const edm::EDGetTokenT<double> rhoToken;

};

#endif
