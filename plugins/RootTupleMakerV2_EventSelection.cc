#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_EventSelection.h"

RootTupleMakerV2_EventSelection::RootTupleMakerV2_EventSelection(const edm::ParameterSet& iConfig) :
  l1InputToken_(consumes<L1GlobalTriggerReadoutRecord>(iConfig.getParameter<edm::InputTag>("L1InputTag"))),
  filterResultsInputToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("FilterResultsInputTag")))
{
  produces <bool> ("isPhysDeclared");
  produces <bool> ("isBPTX0");
  produces <bool> ("isBSCMinBias");
  produces <bool> ("isBSCBeamHalo");
  produces <bool> ("passGoodVertices");
  produces <bool> ("passHBHENoiseFilter");
  produces <bool> ("passHBHENoiseIsoFilter");
  produces <bool> ("passHcalLaserEventFilter");
  produces <bool> ("passCSCTightHaloFilter");
  produces <bool> ("passCSCTightHaloTrkMuUnvetoFilter");
  produces <bool> ("passCSCTightHalo2015Filter");
  produces <bool> ("passGlobalTightHalo2016Filter");
  produces <bool> ("passGlobalSuperTightHalo2016Filter");
  produces <bool> ("passHcalStripHaloFilter");
  //
  produces <bool> ("passEcalDeadCellTriggerPrimitiveFilter");
  produces <bool> ("passEcalDeadCellBoundaryEnergyFilter");
  produces <bool> ("passTrackingFailureFilter");
  produces <bool> ("passEEBadScFilter");
  produces <bool> ("passEcalLaserCorrFilter");
  // 
  produces <bool> ("passTrkPOGFilters");
  produces <bool> ("passChargedHadronTrackResolutionFilter");
  produces <bool> ("passMuonBadTrackFilter");
  produces <bool> ("passManyStripClus53X"       );
  produces <bool> ("passTooManyStripClus53X"    );
  produces <bool> ("passLogErrorTooManyClusters");
}

void RootTupleMakerV2_EventSelection::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<bool> isphysdeclared( new bool() );
  std::auto_ptr<bool> isbptx0( new bool() );
  std::auto_ptr<bool> isbscminbias( new bool() );

  std::auto_ptr<bool> passhbhenoisefilter( new bool() );
  std::auto_ptr<bool> passhbhenoiseisofilter( new bool() );
  std::auto_ptr<bool> passbeamhalofiltertight( new bool() );
  std::auto_ptr<bool> passtightHaloTrkMuUnvetoFilter( new bool() );
  std::auto_ptr<bool> passbeamhalo2015filtertight( new bool() );
  std::auto_ptr<bool> passbeamhalo2016globalfiltertight( new bool() );
  std::auto_ptr<bool> passbeamhalo2016globalfiltersupertight( new bool() );
  std::auto_ptr<bool> passhcalStripHaloFilter( new bool() );
  std::auto_ptr<bool> passhcalLaserEventFilter( new bool() );
  std::auto_ptr<bool> passecalDeadCellTriggerPrimitiveFilter ( new bool() ) ;
  std::auto_ptr<bool> passecalDeadCellBoundaryEnergyFilter( new bool() );
  std::auto_ptr<bool> isprimaryvertex( new bool() );
  std::auto_ptr<bool> passtrackingFailureFilter ( new bool() ) ;
  std::auto_ptr<bool> passbadEESupercrystalFilter ( new bool() ) ;
  std::auto_ptr<bool> passecalLaserCorrFilter (new bool() );
  std::auto_ptr<bool> passtrkPOGFilters( new bool() );
  std::auto_ptr<bool> passchargedHadronTrackResolutionFilter( new bool() );
  std::auto_ptr<bool> passmuonBadTrackFilter( new bool() );
  std::auto_ptr<bool> isbscbeamhalo( new bool() );
  std::auto_ptr<bool> passmanyStripClus53X       (new bool());
  std::auto_ptr<bool> passtooManyStripClus53X    (new bool());
  std::auto_ptr<bool> passlogErrorTooManyClusters(new bool());
  
  *isphysdeclared.get() = false;
  *isbptx0.get() = false;
  *isbscminbias.get() = false;
  *isbscbeamhalo.get() = false;
  *passhbhenoisefilter.get() = true;
  *passhbhenoiseisofilter.get() = true;
  *passbeamhalofiltertight.get() = true;
  *passtightHaloTrkMuUnvetoFilter.get() = true;
  *passbeamhalo2015filtertight.get() = true;
  *passbeamhalo2016globalfiltertight.get() = true;
  *passbeamhalo2016globalfiltersupertight.get() = true;
  *passhcalStripHaloFilter.get() = true;  
  *passhcalLaserEventFilter.get() = true;
  *passecalDeadCellTriggerPrimitiveFilter.get() = true;
  *passecalDeadCellBoundaryEnergyFilter.get() = true;
  *isprimaryvertex.get() = true;
  *passtrackingFailureFilter.get()              = true;
  *passbadEESupercrystalFilter.get()            = true;
  *passecalLaserCorrFilter.get()                = true;
  *passtrkPOGFilters.get() = true;
  *passchargedHadronTrackResolutionFilter.get() = true;
  *passmuonBadTrackFilter.get() = true;
  *passmanyStripClus53X       .get() = true;
  *passtooManyStripClus53X    .get() = true;
  *passlogErrorTooManyClusters.get() = true;
  //-----------------------------------------------------------------
  edm::Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  iEvent.getByToken(l1InputToken_, l1GtReadoutRecord);

  // Technical Trigger Part
  if(l1GtReadoutRecord.isValid()) {
    edm::LogInfo("RootTupleMakerV2_EventSelectionInfo") << "Successfully obtained l1InputTag";

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
    edm::LogError("RootTupleMakerV2_EventSelectionError") << "Error! Can't get the l1InputTag";
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
  // See: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015
  // List of filters: https://github.com/cms-sw/cmssw/blob/CMSSW_7_6_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
  // combination of all filters:
  //  Flag_METFilters
  //# individual filters
  //Flag_HBHENoiseFilter = cms.Path(HBHENoiseFilterResultProducer * HBHENoiseFilter)
  //Flag_HBHENoiseIsoFilter = cms.Path(HBHENoiseFilterResultProducer * HBHENoiseIsoFilter)
  //Flag_CSCTightHaloFilter = cms.Path(CSCTightHaloFilter)
  //Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(CSCTightHaloTrkMuUnvetoFilter)
  //Flag_CSCTightHalo2015Filter = cms.Path(CSCTightHalo2015Filter)
  //Flag_HcalStripHaloFilter = cms.Path(HcalStripHaloFilter)
  //Flag_hcalLaserEventFilter = cms.Path(hcalLaserEventFilter)
  //Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(EcalDeadCellTriggerPrimitiveFilter)
  //Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(EcalDeadCellBoundaryEnergyFilter)
  //Flag_goodVertices = cms.Path(primaryVertexFilter)
  //Flag_trackingFailureFilter = cms.Path(goodVertices + trackingFailureFilter)
  //Flag_eeBadScFilter = cms.Path(eeBadScFilter)
  //Flag_ecalLaserCorrFilter = cms.Path(ecalLaserCorrFilter)
  //Flag_trkPOGFilters = cms.Path(trkPOGFilters)
  //Flag_chargedHadronTrackResolutionFilter = cms.Path(chargedHadronTrackResolutionFilter)
  //Flag_muonBadTrackFilter = cms.Path(muonBadTrackFilter)
  //# and the sub-filters
  //Flag_trkPOG_manystripclus53X = cms.Path(~manystripclus53X)
  //Flag_trkPOG_toomanystripclus53X = cms.Path(~toomanystripclus53X)
  //Flag_trkPOG_logErrorTooManyClusters = cms.Path(~logErrorTooManyClusters)
  
  edm::Handle<edm::TriggerResults> filterResults;
  iEvent.getByToken(filterResultsInputToken_, filterResults);
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

    // Hcal Noise Part (HBHE)
    unsigned int index = filterNames.triggerIndex("Flag_HBHENoiseFilter");
    if(index < filterNames.size())
      *passhbhenoisefilter.get() = filterResults->accept(index);

    // Hcal Noise Iso Part (HBHE)
    index = filterNames.triggerIndex("Flag_HBHENoiseIsoFilter");
    if(index < filterNames.size())
      *passhbhenoiseisofilter.get() = filterResults->accept(index);

    // CSC Beam Halo Tight
    index = filterNames.triggerIndex("Flag_CSCTightHaloFilter");
    if(index < filterNames.size())
      *passbeamhalofiltertight.get() = filterResults->accept(index);

    // CSC Tight Halo Trk Mu Unveto filter
    index = filterNames.triggerIndex("Flag_CSCTightHaloTrkMuUnvetoFilter");
    if(index < filterNames.size())
      *passtightHaloTrkMuUnvetoFilter.get() = filterResults->accept(index);

    //  CSC Beam Halo 2015 Tight filter
    index = filterNames.triggerIndex("Flag_CSCTightHalo2015Filter");
    if(index < filterNames.size())
      *passbeamhalo2015filtertight.get() = filterResults->accept(index);

    //  Global Tight Halo 2016 filter
    index = filterNames.triggerIndex("Flag_globalTightHalo2016Filter");
    if(index < filterNames.size())
      *passbeamhalo2016globalfiltertight.get() = filterResults->accept(index);

    //  Global Super Tight Halo 2016 filter
    index = filterNames.triggerIndex("Flag_globalSuperTightHalo2016Filter");
    if(index < filterNames.size())
      *passbeamhalo2016globalfiltersupertight.get() = filterResults->accept(index);

    //HCAL Strip Halo filter
    index = filterNames.triggerIndex("Flag_HcalStripHaloFilter");
    if(index < filterNames.size())
      *passhcalStripHaloFilter.get() = filterResults->accept(index);

    // Hcal Laser Event Filter
    index = filterNames.triggerIndex("Flag_hcalLaserEventFilter");
    if(index < filterNames.size())
      *passhcalLaserEventFilter.get() = filterResults->accept(index);

    //Ecal dead cell and calo boundary filter
    index = filterNames.triggerIndex("Flag_EcalDeadCellTriggerPrimitiveFilter");
    if(index < filterNames.size())
      *passecalDeadCellTriggerPrimitiveFilter.get() = filterResults->accept(index);

    // ECAL dead cell filter
    index = filterNames.triggerIndex("Flag_EcalDeadCellBoundaryEnergyFilter");
    if(index < filterNames.size())
      *passecalDeadCellBoundaryEnergyFilter.get() = filterResults->accept(index);

    // Good Primary Vertex filter
    index = filterNames.triggerIndex("Flag_goodVertices");
    if(index < filterNames.size())
      *isprimaryvertex.get() = filterResults->accept(index);

    // Tracking failure filter
    index = filterNames.triggerIndex("Flag_trackingFailureFilter");
    if(index < filterNames.size())
      *passtrackingFailureFilter.get() = filterResults->accept(index);

    // Bad EE Supercrystal Filter  
    index = filterNames.triggerIndex("Flag_eeBadScFilter");
    if(index < filterNames.size())
      *passbadEESupercrystalFilter.get() = filterResults->accept(index);

    // large ECAL laser correction filter
    index = filterNames.triggerIndex("Flag_ecalLaserCorrFilter");
    if(index < filterNames.size())
      *passecalLaserCorrFilter.get() = filterResults->accept(index);

    // Tracking POG filters
    index = filterNames.triggerIndex("Flag_trkPOGFilters");
    if(index < filterNames.size())
      *passtrkPOGFilters.get() = filterResults->accept(index);

    // Charged Hadron Track Resolution filter
    index = filterNames.triggerIndex("Flag_chargedHadronTrackResolutionFilter");
    if(index < filterNames.size())
      *passchargedHadronTrackResolutionFilter.get() = filterResults->accept(index);

    // Muon Bad Track filter
    index = filterNames.triggerIndex("Flag_muonBadTrackFilter");
    if(index < filterNames.size())
      *passmuonBadTrackFilter.get() = filterResults->accept(index);

    // Tracking POG filters: ManyStripClus53X
    index = filterNames.triggerIndex("Flag_trkPOG_manystripclus53X");
    if(index < filterNames.size())
      *passmanyStripClus53X.get() = filterResults->accept(index);

    // Tracking POG filters: TooManyStripClus53X
    index = filterNames.triggerIndex("Flag_trkPOG_toomanystripclus53X");
    if(index < filterNames.size())
      *passtooManyStripClus53X.get() = filterResults->accept(index);
    
    // Tracking POG filters: LogErrorTooManyClusters
    index = filterNames.triggerIndex("Flag_trkPOG_logErrorTooManyClusters");
    if(index < filterNames.size())
      *passlogErrorTooManyClusters.get() = filterResults->accept(index);
  
  } else {
    edm::LogError("RootTupleMakerV2_EventSelectionError") << "Error! Can't get the filterResultsInputTag";
  }


  //-----------------------------------------------------------------
  iEvent.put(isphysdeclared,"isPhysDeclared");
  iEvent.put(isbptx0,"isBPTX0");
  iEvent.put(isbscminbias,"isBSCMinBias");
  iEvent.put(isbscbeamhalo,"isBSCBeamHalo");
  iEvent.put(isprimaryvertex,"passGoodVertices");
  iEvent.put(passhbhenoisefilter,"passHBHENoiseFilter");
  iEvent.put(passhbhenoiseisofilter,"passHBHENoiseIsoFilter");
  iEvent.put(passbeamhalofiltertight,"passCSCTightHaloFilter");
  iEvent.put(passtightHaloTrkMuUnvetoFilter,"passCSCTightHaloTrkMuUnvetoFilter");
  iEvent.put(passbeamhalo2015filtertight,"passCSCTightHalo2015Filter");
  iEvent.put(passbeamhalo2016globalfiltertight,"passGlobalTightHalo2016Filter");
  iEvent.put(passbeamhalo2016globalfiltersupertight,"passGlobalSuperTightHalo2016Filter");
  iEvent.put(passhcalStripHaloFilter,"passHcalStripHaloFilter");
  iEvent.put(passhcalLaserEventFilter,"passHcalLaserEventFilter");
  iEvent.put(passecalDeadCellTriggerPrimitiveFilter,"passEcalDeadCellTriggerPrimitiveFilter");
  iEvent.put(passecalDeadCellBoundaryEnergyFilter,"passEcalDeadCellBoundaryEnergyFilter");
  //
  iEvent.put(passtrackingFailureFilter,"passTrackingFailureFilter");
  iEvent.put(passbadEESupercrystalFilter, "passEEBadScFilter");
  iEvent.put(passecalLaserCorrFilter, "passEcalLaserCorrFilter");
  // 
  iEvent.put(passtrkPOGFilters,"passTrkPOGFilters");
  iEvent.put(passchargedHadronTrackResolutionFilter,"passChargedHadronTrackResolutionFilter");
  iEvent.put(passmuonBadTrackFilter,"passMuonBadTrackFilter");
  iEvent.put(passmanyStripClus53X       ,"passManyStripClus53X");
  iEvent.put(passtooManyStripClus53X    ,"passTooManyStripClus53X");
  iEvent.put(passlogErrorTooManyClusters,"passLogErrorTooManyClusters");
}
