#ifndef RootTupleMakerV2EventSelection
#define RootTupleMakerV2EventSelection

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class RootTupleMakerV2_EventSelection : public edm::EDProducer {
 public:
  explicit RootTupleMakerV2_EventSelection(const edm::ParameterSet&);

 private:
  void produce( edm::Event &, const edm::EventSetup & );
  const edm::InputTag   l1InputTag,vtxInputTag;
  const unsigned int    vtxMinNDOF;
  const double          vtxMaxAbsZ, vtxMaxd0;
  const edm::InputTag   trkInputTag;
  const unsigned int    numTracks;
  const double          hpTrackThreshold;
  const edm::InputTag   hcalNoiseInputTag;
  const edm::InputTag   beamHaloInputTag;

  const edm::InputTag trackingFilterJetInputTag   ;
  const double trackingFilterDzTrVtxMax    ;
  const double trackingFilterDxyTrVtxMax   ;
  const double trackingFilterMinSumPtOverHT;

  const edm::InputTag   ecalMaskedCellDRFilterInputTag , caloBoundaryDRFilterInputTag;
  //
  const edm::InputTag   hcalLaserEventFilterInputTag , ecalDeadCellTriggerPrimitiveFilterInputTag , ecalDeadCellBoundaryEnergyFilterInputTag;
  const edm::InputTag   trackingFailureFilterInputTag , badEESupercrystalFilterInputTag, ecalLaserCorrFilterInputTag;
  const edm::InputTag   logErrorTooManyClustersInputTag, manyStripClus53XInputTag, tooManyStripClus53XInputTag;


};

#endif
