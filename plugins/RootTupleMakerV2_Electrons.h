#ifndef RootTupleMakerV2Electrons
#define RootTupleMakerV2Electrons

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/EDCollection.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

class RootTupleMakerV2_Electrons : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_Electrons(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  edm::EDGetTokenT<std::vector<pat::Electron> >electronInputToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxInputToken_;
  edm::EDGetTokenT<double> rhoInputToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotInputToken_;
  edm::EDGetTokenT<EcalRecHitCollection> ebReducedRecHitsInputToken_;
  edm::EDGetTokenT<EcalRecHitCollection> eeReducedRecHitsInputToken_;

  const double          electronIso, muonPt, muonIso;
  const std::string     muonID;
  const edm::InputTag   singleEleTriggerMatchTag, singleEleTriggerMatchWP80Tag, doubleEleTriggerMatchTag;
  const std::string     prefix, suffix;
  const unsigned int    maxSize;
  
  const std::string electronVetoId;
  const std::string electronLooseId;
  const std::string electronMediumId;
  const std::string electronTightId;
  const std::string electronHEEPId;
  const std::string electronMVAIdWP80;
  const std::string electronMVAIdWP90;
  const std::string electronMVAIdHZZ;

  typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;
};

#endif
