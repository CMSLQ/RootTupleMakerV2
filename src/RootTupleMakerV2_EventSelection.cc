#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_EventSelection.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"


RootTupleMakerV2_EventSelection::RootTupleMakerV2_EventSelection(const edm::ParameterSet& iConfig) :
    l1InputTag(iConfig.getParameter<edm::InputTag>("L1InputTag")),
    vtxInputTag(iConfig.getParameter<edm::InputTag>("VertexInputTag")),
    vtxMinNDOF(iConfig.getParameter<unsigned int>("VertexMinimumNDOF")),
    vtxMaxAbsZ(iConfig.getParameter<double>("VertexMaxAbsZ")),
    vtxMaxd0(iConfig.getParameter<double>("VertexMaxd0")),
    trkInputTag(iConfig.getParameter<edm::InputTag>("TracksInputTag")),
    numTracks(iConfig.getParameter<unsigned int>("NumTracks")),
    hpTrackThreshold(iConfig.getParameter<double>("HPTrackThreshold")),
    hcalNoiseInputTag(iConfig.getParameter<edm::InputTag>("HcalNoiseInputTag"))
{
  produces <bool> ("isPhysDeclared");
  produces <bool> ("isBPTX0");
  produces <bool> ("isBSCMinBias");
  produces <bool> ("isBSCBeamHalo");
  produces <bool> ("isPrimaryVertex");
  produces <bool> ("isBeamScraping");
  produces <bool> ("passLooseNoiseFilter");
  produces <bool> ("passTightNoiseFilter");
  produces <bool> ("passHighLevelNoiseFilter");
}

void RootTupleMakerV2_EventSelection::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<bool> isphysdeclared( new bool() );
  std::auto_ptr<bool> isbptx0( new bool() );
  std::auto_ptr<bool> isbscminbias( new bool() );
  std::auto_ptr<bool> isbscbeamhalo( new bool() );
  std::auto_ptr<bool> isprimaryvertex( new bool() );
  std::auto_ptr<bool> isbeamscraping( new bool() );
  std::auto_ptr<bool> passloosenoisefilter( new bool() );
  std::auto_ptr<bool> passtightnoisefilter( new bool() );
  std::auto_ptr<bool> passhiglevelnoisefilter( new bool() );

  *isphysdeclared.get() = false;
  *isbptx0.get() = false;
  *isbscminbias.get() = false;
  *isbscbeamhalo.get() = false;
  *isprimaryvertex.get() = false;
  *isbeamscraping.get() = false;
  *passloosenoisefilter.get() = false;
  *passtightnoisefilter.get() = false;
  *passhiglevelnoisefilter.get() = false;

  //-----------------------------------------------------------------
  edm::Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  iEvent.getByLabel(l1InputTag, l1GtReadoutRecord);

  // Technical Trigger Part
  if(l1GtReadoutRecord.isValid()) {
    edm::LogInfo("RootTupleMakerV2_EventSelectionInfo") << "Successfully obtained " << l1InputTag;

    L1GtFdlWord fdlWord = l1GtReadoutRecord->gtFdlWord();
    if (fdlWord.physicsDeclared() == 1)
      *isphysdeclared.get() = true;

    // BPTX0
    if ( l1GtReadoutRecord->technicalTriggerWord()[0] )
      *isbptx0.get() = true;

    // MinBias
    if ( l1GtReadoutRecord->technicalTriggerWord()[40] || l1GtReadoutRecord->technicalTriggerWord()[41] )
      *isbscminbias.get() = true;

    // BeamHalo
    if ( l1GtReadoutRecord->technicalTriggerWord()[36] || l1GtReadoutRecord->technicalTriggerWord()[37] || 
         l1GtReadoutRecord->technicalTriggerWord()[38] || l1GtReadoutRecord->technicalTriggerWord()[39] )
      *isbscbeamhalo.get() = true;

  } else {
    edm::LogError("RootTupleMakerV2_EventSelectionError") << "Error! Can't get the product " << l1InputTag;
  }

  // Good Primary Vertex Part
  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(vtxInputTag,primaryVertices);

  if(primaryVertices.isValid()) {
    edm::LogInfo("RootTupleMakerV2_EventSelectionInfo") << "Total # Primary Vertices: " << primaryVertices->size();

    for( reco::VertexCollection::const_iterator it=primaryVertices->begin() ; it!=primaryVertices->end() ; ++it ) {
      if( !(it->isFake()) && it->ndof() > vtxMinNDOF &&
          fabs(it->z()) <= vtxMaxAbsZ && fabs(it->position().rho()) <= vtxMaxd0
        ) *isprimaryvertex.get() = true;
    }
  } else {
    edm::LogError("RootTupleMakerV2_EventSelectionError") << "Error! Can't get the product " << vtxInputTag;
  }

  // Scraping Events Part
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(trkInputTag,tracks);

  if(tracks.isValid()) {
    edm::LogInfo("RootTupleMakerV2_EventSelectionInfo") << "Total # Tracks: " << tracks->size();

    int numhighpurity = 0;
    double fraction = 1.;
    reco::TrackBase::TrackQuality trackQuality = reco::TrackBase::qualityByName("highPurity");

    if( tracks->size() > numTracks ){
      for( reco::TrackCollection::const_iterator it=tracks->begin(); it!=tracks->end(); ++it ) {
        if( it->quality(trackQuality) ) numhighpurity++;
      }
      fraction = (double)numhighpurity/(double)tracks->size();
      if( fraction < hpTrackThreshold ) *isbeamscraping.get() = true;
    }
  } else {
    edm::LogError("RootTupleMakerV2_EventSelectionError") << "Error! Can't get the product " << trkInputTag;
  }

  // Hcal Noise Part
  edm::Handle<HcalNoiseSummary> hcalNoise;
  iEvent.getByLabel(hcalNoiseInputTag,hcalNoise);

  if(hcalNoise.isValid()) {
    edm::LogInfo("RootTupleMakerV2_EventSelectionInfo") << "Successfully obtained " << hcalNoiseInputTag;;

    *passloosenoisefilter.get() = hcalNoise->passLooseNoiseFilter();
    *passtightnoisefilter.get() = hcalNoise->passTightNoiseFilter();
    *passhiglevelnoisefilter.get() = hcalNoise->passHighLevelNoiseFilter();
  } else {
    edm::LogError("RootTupleMakerV2_EventSelectionError") << "Error! Can't get the product " << hcalNoiseInputTag;
  }

  //-----------------------------------------------------------------
  iEvent.put(isphysdeclared,"isPhysDeclared");
  iEvent.put(isbptx0,"isBPTX0");
  iEvent.put(isbscminbias,"isBSCMinBias");
  iEvent.put(isbscbeamhalo,"isBSCBeamHalo");
  iEvent.put(isprimaryvertex,"isPrimaryVertex");
  iEvent.put(isbeamscraping,"isBeamScraping");
  iEvent.put(passloosenoisefilter,"passLooseNoiseFilter");
  iEvent.put(passtightnoisefilter,"passTightNoiseFilter");
  iEvent.put(passhiglevelnoisefilter,"passHighLevelNoiseFilter");
}
