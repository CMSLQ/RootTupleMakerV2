#ifndef RootTupleMakerV2PFJets
#define RootTupleMakerV2PFJets

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include <string>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "RecoJets/JetProducers/interface/PileupJetIdAlgo.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

class RootTupleMakerV2_PFJets : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_PFJets(const edm::ParameterSet&);

 private:
  std::string upperCase(std::string input);
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::EDGetTokenT<std::vector<pat::Jet> >   jetInputToken_;
  const edm::InputTag   inputTag;
  const std::string     prefix,suffix,mvaPileupIDname;
  const unsigned int    maxSize;
  const std::string     jecUncPath; 
  const bool            readJECuncertainty;
  const bool            readJERuncertainty;
  const edm::EDGetTokenT<reco::VertexCollection>   vtxInputToken_;
  bool            isPuppiJetColl;
  const std::string jerUncPath;
  const bool jer_from_gt;
  const edm::FileInPath jer_resolutions_file;
  const edm::FileInPath jer_scale_factors_file;
  const edm::FileInPath jec_sources_file;
  const edm::EDGetTokenT<double> rhoToken;
};

#endif
