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
#include "DataFormats/METReco/interface/BeamHaloSummary.h"

RootTupleMakerV2_EventSelection::RootTupleMakerV2_EventSelection(const edm::ParameterSet& iConfig) :
    l1InputTag(iConfig.getParameter<edm::InputTag>("L1InputTag")),
    vtxInputTag(iConfig.getParameter<edm::InputTag>("VertexInputTag")),
    vtxMinNDOF(iConfig.getParameter<unsigned int>("VertexMinimumNDOF")),
    vtxMaxAbsZ(iConfig.getParameter<double>("VertexMaxAbsZ")),
    vtxMaxd0(iConfig.getParameter<double>("VertexMaxd0")),
    trkInputTag(iConfig.getParameter<edm::InputTag>("TracksInputTag")),
    numTracks(iConfig.getParameter<unsigned int>("NumTracks")),
    hpTrackThreshold(iConfig.getParameter<double>("HPTrackThreshold")),
    hcalNoiseInputTag(iConfig.getParameter<edm::InputTag>("HcalNoiseInputTag")),
    beamHaloInputTag(iConfig.getParameter<edm::InputTag>("BeamHaloInputTag")),
    trackingFilterJetInputTag   (iConfig.getParameter<edm::InputTag>("TrackingFailureJets")),	      
    trackingFilterDzTrVtxMax    (iConfig.getParameter<double>       ("TrackingFailureDzTrVtzMax")),   
    trackingFilterDxyTrVtxMax   (iConfig.getParameter<double>       ("TrackingFailureDxyTrVtxMax")) ,
    trackingFilterMinSumPtOverHT(iConfig.getParameter<double>       ("TrackingFailureMinSumPtOverHT")),
    ecalMaskedCellDRFilterInputTag(iConfig.getParameter<edm::InputTag>("EcalMaskedCellDRFilterInputTag")),
    caloBoundaryDRFilterInputTag(iConfig.getParameter<edm::InputTag>("CaloBoundaryDRFilterInputTag"))
{
  produces <bool> ("isPhysDeclared");
  produces <bool> ("isBPTX0");
  produces <bool> ("isBSCMinBias");
  produces <bool> ("isBSCBeamHalo");
  produces <bool> ("isPrimaryVertex");
  produces <bool> ("isBeamScraping");
  produces <bool> ("passHBHENoiseFilter");
  produces <bool> ("passBeamHaloFilterLoose");
  produces <bool> ("passBeamHaloFilterTight");
  produces <bool> ("isTrackingFailure");
  produces <bool> ("passEcalMaskedCellDRFilter");
  produces <bool> ("passCaloBoundaryDRFilter");
}

void RootTupleMakerV2_EventSelection::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<bool> isphysdeclared( new bool() );
  std::auto_ptr<bool> isbptx0( new bool() );
  std::auto_ptr<bool> isbscminbias( new bool() );
  std::auto_ptr<bool> isbscbeamhalo( new bool() );
  std::auto_ptr<bool> isprimaryvertex( new bool() );
  std::auto_ptr<bool> isbeamscraping( new bool() );
  std::auto_ptr<bool> passhbhenoisefilter( new bool() );
  std::auto_ptr<bool> passbeamhalofilterloose( new bool() );
  std::auto_ptr<bool> passbeamhalofiltertight( new bool() );
  std::auto_ptr<bool> istrackingfailure ( new bool() ) ;
  std::auto_ptr<bool> passEcalMaskedCellDRFilter ( new bool() ) ;
  std::auto_ptr<bool> passCaloBoundaryDRFilter ( new bool() ) ;

  *isphysdeclared.get() = false;
  *isbptx0.get() = false;
  *isbscminbias.get() = false;
  *isbscbeamhalo.get() = false;
  *isprimaryvertex.get() = false;
  *isbeamscraping.get() = false;
  *passhbhenoisefilter.get() = true;
  *passbeamhalofilterloose.get() = true;
  *passbeamhalofiltertight.get() = true;
  *passEcalMaskedCellDRFilter.get() = true;
  *passCaloBoundaryDRFilter.get() = true;

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
    if ( (l1GtReadoutRecord->technicalTriggerWord()[36] || l1GtReadoutRecord->technicalTriggerWord()[37] ||
          l1GtReadoutRecord->technicalTriggerWord()[38] || l1GtReadoutRecord->technicalTriggerWord()[39]) ||
         ((l1GtReadoutRecord->technicalTriggerWord()[42] && !l1GtReadoutRecord->technicalTriggerWord()[43]) ||
          (l1GtReadoutRecord->technicalTriggerWord()[43] && !l1GtReadoutRecord->technicalTriggerWord()[42])) )
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
  edm::Handle<bool> hbheFilterResult;
  iEvent.getByLabel(hcalNoiseInputTag, hbheFilterResult);

  if(hbheFilterResult.isValid()) {
    edm::LogInfo("RootTupleMakerV2_EventSelectionInfo") << "Successfully obtained " << hcalNoiseInputTag;

      *passhbhenoisefilter.get()=*hbheFilterResult;
  } else {
    edm::LogError("RootTupleMakerV2_EventSelectionError") << "Error! Can't get the product " << hcalNoiseInputTag;
  }

  // Beam Halo part
  edm::Handle<reco::BeamHaloSummary> TheBeamHaloSummary;
  iEvent.getByLabel(beamHaloInputTag,TheBeamHaloSummary); 

  if(TheBeamHaloSummary.isValid()) {
    edm::LogInfo("RootTupleMakerV2_EventSelectionInfo") << "Successfully obtained " << beamHaloInputTag;
    const reco::BeamHaloSummary TheSummary = (*TheBeamHaloSummary.product() );
    *passbeamhalofilterloose.get() = !TheSummary.CSCLooseHaloId();
    *passbeamhalofiltertight.get() = !TheSummary.CSCTightHaloId();    
  } else {
    edm::LogError("RootTupleMakerV2_EventSelectionError") << "Error! Can't get the product " << beamHaloInputTag;
  }

  //Tracking failure filter:
  edm::Handle<edm::View<reco::Jet> > jets;
  iEvent.getByLabel(trackingFilterJetInputTag, jets);
  
  double ht = 0;
  for (edm::View<reco::Jet>::const_iterator j = jets->begin(); j != jets->end(); ++j) {
    ht += j->pt();
  }

  double sumpt = 0;
  if (primaryVertices->size() > 0) {
    const reco::Vertex * vtx = &((*primaryVertices)[0]);
    for (std::vector<reco::Track>::const_iterator tr = tracks->begin(); tr != tracks->end(); ++tr) {
      if (fabs(tr->dz(vtx->position()))  > trackingFilterDzTrVtxMax    ) continue;
      if (fabs(tr->dxy(vtx->position())) > trackingFilterDxyTrVtxMax   ) continue;
      sumpt += tr->pt();
    }
  }
  
  *istrackingfailure.get() = ((sumpt/ht) < trackingFilterMinSumPtOverHT );

  
  //Ecal maksed cell and calo boundary filter ( https://twiki.cern.ch/twiki/bin/view/CMS/SusyEcalMaskedCellSummary ):
  edm::Handle<int> EcalMaskedCellDRFilterResult;
  iEvent.getByLabel(ecalMaskedCellDRFilterInputTag, EcalMaskedCellDRFilterResult);
  
  edm::Handle<int> CaloBoundaryDRFilterResult;
  iEvent.getByLabel(caloBoundaryDRFilterInputTag, CaloBoundaryDRFilterResult);
  
  if(EcalMaskedCellDRFilterResult.isValid()) {
    *passEcalMaskedCellDRFilter.get()=!(*EcalMaskedCellDRFilterResult);
  }
  //    } else {
  //      edm::LogError("RootTupleMakerV2_EventSelectionError") << "Error! Can't get the product " << ecalMaskedCellDRFilterInputTag;
  //    }
  
  if(CaloBoundaryDRFilterResult.isValid()) {
    *passCaloBoundaryDRFilter.get()=!(*CaloBoundaryDRFilterResult);
  }
  //    } else {
  //      edm::LogError("RootTupleMakerV2_EventSelectionError") << "Error! Can't get the product " << caloBoundaryDRFilterInputTag;
  //    }
  

  //-----------------------------------------------------------------
  iEvent.put(isphysdeclared,"isPhysDeclared");
  iEvent.put(isbptx0,"isBPTX0");
  iEvent.put(isbscminbias,"isBSCMinBias");
  iEvent.put(isbscbeamhalo,"isBSCBeamHalo");
  iEvent.put(isprimaryvertex,"isPrimaryVertex");
  iEvent.put(isbeamscraping,"isBeamScraping");
  iEvent.put(passhbhenoisefilter,"passHBHENoiseFilter");
  iEvent.put(passbeamhalofilterloose,"passBeamHaloFilterLoose");
  iEvent.put(passbeamhalofiltertight,"passBeamHaloFilterTight");
  iEvent.put(istrackingfailure, "isTrackingFailure");
  iEvent.put(passEcalMaskedCellDRFilter, "passEcalMaskedCellDRFilter");
  iEvent.put(passCaloBoundaryDRFilter, "passCaloBoundaryDRFilter");

}
