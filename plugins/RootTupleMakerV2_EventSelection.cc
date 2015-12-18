#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_EventSelection.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

RootTupleMakerV2_EventSelection::RootTupleMakerV2_EventSelection(const edm::ParameterSet& iConfig) :
  l1InputTag(iConfig.getParameter<edm::InputTag>("L1InputTag")),
  hcalNoiseInputTag(iConfig.getParameter<edm::InputTag>("HcalNoiseInputTag")),
  filterResultsInputTag(iConfig.getParameter<edm::InputTag>("FilterResultsInputTag"))
{
  produces <bool> ("isPhysDeclared");
  produces <bool> ("isBPTX0");
  produces <bool> ("isBSCMinBias");
  produces <bool> ("isBSCBeamHalo");
  produces <bool> ("isPrimaryVertex");
  produces <bool> ("passHBHENoiseFilter");
  produces <bool> ("passHBHENoiseIsoFilter");
  produces <bool> ("passHcalLaserEventFilter");
  produces <bool> ("passBeamHaloFilterTight");
  produces <bool> ("isTrackingFailure");
  //
  produces <bool> ("passEcalDeadCellTriggerPrimitiveFilter");
  produces <bool> ("passTrackingFailureFilter");
  produces <bool> ("passBadEESupercrystalFilter");
  produces <bool> ("passEcalLaserCorrFilter");
  // 
  produces <bool> ("passLogErrorTooManyClusters");
  produces <bool> ("passManyStripClus53X"       );
  produces <bool> ("passTooManyStripClus53X"    );
}

void RootTupleMakerV2_EventSelection::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<bool> isphysdeclared( new bool() );
  std::auto_ptr<bool> isbptx0( new bool() );
  std::auto_ptr<bool> isbscminbias( new bool() );
  std::auto_ptr<bool> isbscbeamhalo( new bool() );
  std::auto_ptr<bool> isprimaryvertex( new bool() );
  std::auto_ptr<bool> passhbhenoisefilter( new bool() );
  std::auto_ptr<bool> passhbhenoiseisofilter( new bool() );
  std::auto_ptr<bool> passHcalLaserEventFilter( new bool() );
  std::auto_ptr<bool> passbeamhalofiltertight( new bool() );
  std::auto_ptr<bool> istrackingfailure ( new bool() ) ;
  //
  std::auto_ptr<bool> passEcalDeadCellTriggerPrimitiveFilter ( new bool() ) ;
  std::auto_ptr<bool> passTrackingFailureFilter ( new bool() ) ;
  std::auto_ptr<bool> passBadEESupercrystalFilter ( new bool() ) ;
  std::auto_ptr<bool> passEcalLaserCorrFilter (new bool() );
  //
  std::auto_ptr<bool> passLogErrorTooManyClusters(new bool());
  std::auto_ptr<bool> passManyStripClus53X       (new bool());
  std::auto_ptr<bool> passTooManyStripClus53X    (new bool());
  
  *isphysdeclared.get() = false;
  *isbptx0.get() = false;
  *isbscminbias.get() = false;
  *isbscbeamhalo.get() = false;
  *isprimaryvertex.get() = true;
  *passhbhenoisefilter.get() = true;
  *passhbhenoiseisofilter.get() = true;
  *passHcalLaserEventFilter.get() = true;
  *passbeamhalofiltertight.get() = true;
  //
  *passEcalDeadCellTriggerPrimitiveFilter.get() = true;
  *passTrackingFailureFilter.get()              = true;
  *passBadEESupercrystalFilter.get()            = true;
  *passEcalLaserCorrFilter.get()                = true;
  // 
  *passLogErrorTooManyClusters.get() = true;
  *passManyStripClus53X       .get() = true;
  *passTooManyStripClus53X    .get() = true;
  
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


  //XXX FIXME There are no tracks in MiniAOD
  //// Scraping Events Part
  //edm::Handle<reco::TrackCollection> tracks;
  //iEvent.getByLabel(trkInputTag,tracks);
  //if(tracks.isValid()) {
  //  edm::LogInfo("RootTupleMakerV2_EventSelectionInfo") << "Total # Tracks: " << tracks->size();
  //
  //  int numhighpurity = 0;
  //  double fraction = 1.;
  //  reco::TrackBase::TrackQuality trackQuality = reco::TrackBase::qualityByName("highPurity");
  //
  //  if( tracks->size() > numTracks ){
  //    for( reco::TrackCollection::const_iterator it=tracks->begin(); it!=tracks->end(); ++it ) {
  //      if( it->quality(trackQuality) ) numhighpurity++;
  //    }
  //    fraction = (double)numhighpurity/(double)tracks->size();
  //    if( fraction < hpTrackThreshold ) *isbeamscraping.get() = true;
  //  }
  //} else {
  //  edm::LogError("RootTupleMakerV2_EventSelectionError") << "Error! Can't get the product " << trkInputTag;
  //}

  // These filters are now run in MiniAOD and saved as an edm::TriggerResults of the PAT process
  // See: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
  // List of filters: https://github.com/cms-sw/cmssw/blob/CMSSW_7_0_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
  // combination of all filters:
  //  Flag_METFilters
  // individual filters:
  //  Flag_HBHENoiseFilter = cms.Path(HBHENoiseFilter)
  //  Flag_HBHENoiseIsoFilter = cms.Path(HBHENoiseIsoFilter)
  //  Flag_CSCTightHaloFilter = cms.Path(CSCTightHaloFilter)
  //  Flag_hcalLaserEventFilter = cms.Path(hcalLaserEventFilter)
  //  Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(EcalDeadCellTriggerPrimitiveFilter)
  //  Flag_goodVertices = cms.Path(goodVertices)
  //  Flag_trackingFailureFilter = cms.Path(goodVertices + trackingFailureFilter)
  //  Flag_eeBadScFilter = cms.Path(eeBadScFilter)
  //  Flag_ecalLaserCorrFilter = cms.Path(ecalLaserCorrFilter)
  //  Flag_trkPOGFilters = cms.Path(trkPOGFilters) --> This is made up of the three trkPOG filters below
  // and the individual trkPOG filters
  //  Flag_trkPOG_manystripclus53X = cms.Path(~manystripclus53X)
  //  Flag_trkPOG_toomanystripclus53X = cms.Path(~toomanystripclus53X)
  //  Flag_trkPOG_logErrorTooManyClusters = cms.Path(~logErrorTooManyClusters)

  edm::Handle<edm::TriggerResults> filterResults;
  iEvent.getByLabel(filterResultsInputTag, filterResults);
  if(filterResults.isValid())
  {
    //edm::LogInfo("RootTupleMakerV2_EventSelection") << "Successfully obtained " << filterResultsInputTag;
    const edm::TriggerNames &filterNames = iEvent.triggerNames(*filterResults);
    // print filter names
    //std::cout << "\n === FILTER NAMES === " << std::endl;
    //for (unsigned int i = 0, n = filterResults->size(); i < n; ++i)
    //{
    //  std::cout << "Filter " << filterNames.triggerName(i) << 
    //    ": " << (filterResults->accept(i) ? "PASS" : "fail (or not run)") 
    //    << std::endl;
    //}

    // CSC Beam Halo Tight
    unsigned int index = filterNames.triggerIndex("Flag_CSCTightHaloFilter");
    if(index < filterNames.size())
      *passbeamhalofiltertight.get() = filterResults->accept(index);

    // Tracking failure filter
    index = filterNames.triggerIndex("Flag_trackingFailureFilter");
    if(index < filterNames.size())
      *passTrackingFailureFilter.get() = filterResults->accept(index);

    // ECAL dead cell filter
    index = filterNames.triggerIndex("Flag_EcalDeadCellTriggerPrimitiveFilter");
    if(index < filterNames.size())
      *passEcalDeadCellTriggerPrimitiveFilter.get() = filterResults->accept(index);


    //Ecal masked cell and calo boundary filter ( https://twiki.cern.ch/twiki/bin/view/CMS/SusyEcalMaskedCellSummary )
    //XXX: FIXME? Not run in MiniAOD

    // Tracking POG filters: ManyStripClus53X
    index = filterNames.triggerIndex("Flag_trkPOG_manystripclus53X");
    if(index < filterNames.size())
      *passManyStripClus53X.get() = filterResults->accept(index);

    // Tracking POG filters: TooManyStripClus53X
    index = filterNames.triggerIndex("Flag_trkPOG_toomanystripclus53X");
    if(index < filterNames.size())
      *passTooManyStripClus53X.get() = filterResults->accept(index);
    
    // Tracking POG filters: LogErrorTooManyClusters
    index = filterNames.triggerIndex("Flag_trkPOG_logErrorTooManyClusters");
    if(index < filterNames.size())
      *passLogErrorTooManyClusters.get() = filterResults->accept(index);

    // Bad EE Supercrystal Filter  
    index = filterNames.triggerIndex("Flag_eeBadScFilter");
    if(index < filterNames.size())
      *passBadEESupercrystalFilter.get() = filterResults->accept(index);

    // large ECAL laser correction filter
    index = filterNames.triggerIndex("Flag_ecalLaserCorrFilter");
    if(index < filterNames.size())
      *passEcalLaserCorrFilter.get() = filterResults->accept(index);

    // Good Primary Vertex Part
    index = filterNames.triggerIndex("Flag_goodVertices");
    if(index < filterNames.size())
      *isprimaryvertex.get() = filterResults->accept(index);

    // Hcal Noise Part (HBHE)
    index = filterNames.triggerIndex("Flag_HBHENoiseFilter");
    if(index < filterNames.size())
      *passhbhenoisefilter.get() = filterResults->accept(index);

   // Hcal Noise Iso Part (HBHE)
    index = filterNames.triggerIndex("Flag_HBHENoiseIsoFilter");
    if(index < filterNames.size())
      *passhbhenoiseisofilter.get() = filterResults->accept(index);

    // Hcal Laser Event Filter
    index = filterNames.triggerIndex("Flag_hcalLaserEventFilter");
    if(index < filterNames.size())
      *passHcalLaserEventFilter.get() = filterResults->accept(index);

  } else {
    edm::LogError("RootTupleMakerV2_EventSelectionError") << "Error! Can't get the product " << filterResultsInputTag;
  }


  //-----------------------------------------------------------------
  iEvent.put(isphysdeclared,"isPhysDeclared");
  iEvent.put(isbptx0,"isBPTX0");
  iEvent.put(isbscminbias,"isBSCMinBias");
  iEvent.put(isbscbeamhalo,"isBSCBeamHalo");
  iEvent.put(isprimaryvertex,"isPrimaryVertex");
  iEvent.put(passhbhenoisefilter,"passHBHENoiseFilter");
  iEvent.put(passhbhenoiseisofilter,"passHBHENoiseIsoFilter");
  iEvent.put(passHcalLaserEventFilter,"passHcalLaserEventFilter");
  iEvent.put(passbeamhalofiltertight,"passBeamHaloFilterTight");
  iEvent.put(istrackingfailure, "isTrackingFailure");
  //
  iEvent.put(passEcalDeadCellTriggerPrimitiveFilter,"passEcalDeadCellTriggerPrimitiveFilter");
  iEvent.put(passTrackingFailureFilter,"passTrackingFailureFilter");
  iEvent.put(passBadEESupercrystalFilter, "passBadEESupercrystalFilter");
  iEvent.put(passEcalLaserCorrFilter, "passEcalLaserCorrFilter");
  // 
  iEvent.put(passLogErrorTooManyClusters,"passLogErrorTooManyClusters");
  iEvent.put(passManyStripClus53X       ,"passManyStripClus53X");
  iEvent.put(passTooManyStripClus53X    ,"passTooManyStripClus53X");
  
}
