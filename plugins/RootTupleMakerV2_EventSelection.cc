#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_EventSelection.h"

RootTupleMakerV2_EventSelection::RootTupleMakerV2_EventSelection(const edm::ParameterSet& iConfig) :
  l1InputToken_(consumes<L1GlobalTriggerReadoutRecord>(iConfig.getParameter<edm::InputTag>("L1InputTag"))),
  filterResultsInputToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("FilterResultsInputTag")))
{
  produces <bool> ("passGoodVertices");
  produces <bool> ("passHBHENoiseFilter");
  produces <bool> ("passHBHENoiseIsoFilter");
  produces <bool> ("passCSCTightHaloFilter");
  produces <bool> ("passCSCTightHaloTrkMuUnvetoFilter");
  produces <bool> ("passGlobalTightHalo2016Filter");
  produces <bool> ("passGlobalSuperTightHalo2016Filter");
  //
  produces <bool> ("passBadPFMuonFilter");
  produces <bool> ("passBadChargedCandidateFilter");
  //
  produces <bool> ("passEcalDeadCellTriggerPrimitiveFilter");
  produces <bool> ("passEEBadScFilter");
  // 
  produces <bool> ("passChargedHadronTrackResolutionFilter");
  produces <bool> ("passMuonBadTrackFilter");
  //
  produces <bool> ("badMuonsFlag");
  produces <bool> ("duplicateMuonsFlag");
  produces <bool> ("noBadMuonsFlag");
}

void RootTupleMakerV2_EventSelection::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {


  std::unique_ptr<bool> passhbhenoisefilter( new bool() );
  std::unique_ptr<bool> passhbhenoiseisofilter( new bool() );
  std::unique_ptr<bool> passbeamhalofiltertight( new bool() );
  std::unique_ptr<bool> passtightHaloTrkMuUnvetoFilter( new bool() );
  std::unique_ptr<bool> passbeamhalo2016globalfiltertight( new bool() );
  std::unique_ptr<bool> passbeamhalo2016globalfiltersupertight( new bool() );
  std::unique_ptr<bool> passbadpfmuonFilter( new bool() );
  std::unique_ptr<bool> passbadchargedcandidateFilter( new bool() );
  std::unique_ptr<bool> passecalDeadCellTriggerPrimitiveFilter ( new bool() ) ;
  std::unique_ptr<bool> isprimaryvertex( new bool() );
  std::unique_ptr<bool> passbadEESupercrystalFilter ( new bool() ) ;
  std::unique_ptr<bool> passchargedHadronTrackResolutionFilter( new bool() );
  std::unique_ptr<bool> passmuonBadTrackFilter( new bool() );
  std::unique_ptr<bool> flagBadMuons( new bool() );
  std::unique_ptr<bool> flagDuplicateMuons( new bool() );
  std::unique_ptr<bool> flagNoBadMuons( new bool() );
  
  *passhbhenoisefilter.get() = true;
  *passhbhenoiseisofilter.get() = true;
  *passbeamhalofiltertight.get() = true;
  *passtightHaloTrkMuUnvetoFilter.get() = true;
  *passbeamhalo2016globalfiltertight.get() = true;
  *passbeamhalo2016globalfiltersupertight.get() = true;
  *passbadpfmuonFilter.get() = true;
  *passbadchargedcandidateFilter.get() = true;
  *passecalDeadCellTriggerPrimitiveFilter.get() = true;
  *isprimaryvertex.get() = true;
  *passbadEESupercrystalFilter.get()            = true;
  *passchargedHadronTrackResolutionFilter.get() = true;
  *passmuonBadTrackFilter.get() = true;
  *flagBadMuons.get()       = false;
  *flagDuplicateMuons.get() = false;
  *flagNoBadMuons.get()     = true;
  //-----------------------------------------------------------------
  edm::Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  iEvent.getByToken(l1InputToken_, l1GtReadoutRecord);

  /*
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
  
  */

  edm::Handle<edm::TriggerResults> filterResults;
  iEvent.getByToken(filterResultsInputToken_, filterResults);
  if(filterResults.isValid())
  {
    //edm::LogInfo("RootTupleMakerV2_EventSelection") << "Successfully obtained " << filterResultsInputTag;
    const edm::TriggerNames &filterNames = iEvent.triggerNames(*filterResults);
    //// print filter names
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

    //  Global Tight Halo 2016 filter
    index = filterNames.triggerIndex("Flag_globalTightHalo2016Filter");
    if(index < filterNames.size())
      *passbeamhalo2016globalfiltertight.get() = filterResults->accept(index);

    //  Global Super Tight Halo 2016 filter
    index = filterNames.triggerIndex("Flag_globalSuperTightHalo2016Filter");
    if(index < filterNames.size())
      *passbeamhalo2016globalfiltersupertight.get() = filterResults->accept(index);

    ////HCAL Strip Halo filter
    //index = filterNames.triggerIndex("Flag_HcalStripHaloFilter");
    //if(index < filterNames.size())
    //  *passhcalStripHaloFilter.get() = filterResults->accept(index);

    //// Hcal Laser Event Filter
    //index = filterNames.triggerIndex("Flag_hcalLaserEventFilter");
    //if(index < filterNames.size())
    //  *passhcalLaserEventFilter.get() = filterResults->accept(index);

    //Ecal dead cell and calo boundary filter
    index = filterNames.triggerIndex("Flag_EcalDeadCellTriggerPrimitiveFilter");
    if(index < filterNames.size())
      *passecalDeadCellTriggerPrimitiveFilter.get() = filterResults->accept(index);

    //// ECAL dead cell filter
    //index = filterNames.triggerIndex("Flag_EcalDeadCellBoundaryEnergyFilter");
    //if(index < filterNames.size())
    //  *passecalDeadCellBoundaryEnergyFilter.get() = filterResults->accept(index);

    // Good Primary Vertex filter
    index = filterNames.triggerIndex("Flag_goodVertices");
    if(index < filterNames.size())
      *isprimaryvertex.get() = filterResults->accept(index);

    //// Tracking failure filter
    //index = filterNames.triggerIndex("Flag_trackingFailureFilter");
    //if(index < filterNames.size())
    //  *passtrackingFailureFilter.get() = filterResults->accept(index);

    // Bad EE Supercrystal Filter  
    index = filterNames.triggerIndex("Flag_eeBadScFilter");
    if(index < filterNames.size())
      *passbadEESupercrystalFilter.get() = filterResults->accept(index);

    //// large ECAL laser correction filter
    //index = filterNames.triggerIndex("Flag_ecalLaserCorrFilter");
    //if(index < filterNames.size())
    //  *passecalLaserCorrFilter.get() = filterResults->accept(index);

    //// Tracking POG filters
    //index = filterNames.triggerIndex("Flag_trkPOGFilters");
    //if(index < filterNames.size())
    //  *passtrkPOGFilters.get() = filterResults->accept(index);

    // Charged Hadron Track Resolution filter
    index = filterNames.triggerIndex("Flag_chargedHadronTrackResolutionFilter");
    if(index < filterNames.size())
      *passchargedHadronTrackResolutionFilter.get() = filterResults->accept(index);

    // Muon Bad Track filter
    index = filterNames.triggerIndex("Flag_muonBadTrackFilter");
    if(index < filterNames.size())
      *passmuonBadTrackFilter.get() = filterResults->accept(index);

    // Muon Bad or Duplicate Flags
    index = filterNames.triggerIndex("Flag_badMuons");
    if(index < filterNames.size())
      *flagBadMuons.get() = filterResults->accept(index);

    index = filterNames.triggerIndex("Flag_duplicateMuons");
    if(index < filterNames.size())
      *flagDuplicateMuons.get() = filterResults->accept(index);

    index = filterNames.triggerIndex("Flag_noBadMuons");
    if(index < filterNames.size())
      *flagNoBadMuons.get() = filterResults->accept(index);

    index = filterNames.triggerIndex("Flag_BadChargedCandidateFilter");
    if(index < filterNames.size())
      *passbadchargedcandidateFilter.get() = filterResults->accept(index);

    index = filterNames.triggerIndex("Flag_BadPFMuonFilter");
    if(index < filterNames.size())
      *passbadpfmuonFilter.get() = filterResults->accept(index);

  } else {
    edm::LogError("RootTupleMakerV2_EventSelectionError") << "Error! Can't get the filterResultsInputTag";
  }

  //-----------------------------------------------------------------
  iEvent.put(std::move(isprimaryvertex),"passGoodVertices");
  iEvent.put(std::move(passhbhenoisefilter),"passHBHENoiseFilter");
  iEvent.put(std::move(passhbhenoiseisofilter),"passHBHENoiseIsoFilter");
  iEvent.put(std::move(passbeamhalofiltertight),"passCSCTightHaloFilter");
  iEvent.put(std::move(passtightHaloTrkMuUnvetoFilter),"passCSCTightHaloTrkMuUnvetoFilter");
  iEvent.put(std::move(passbeamhalo2016globalfiltertight),"passGlobalTightHalo2016Filter");
  iEvent.put(std::move(passbeamhalo2016globalfiltersupertight),"passGlobalSuperTightHalo2016Filter");
  iEvent.put(std::move(passecalDeadCellTriggerPrimitiveFilter),"passEcalDeadCellTriggerPrimitiveFilter");
  //
  iEvent.put(std::move(passbadEESupercrystalFilter), "passEEBadScFilter");
  //
  iEvent.put(std::move(passbadpfmuonFilter), "passBadPFMuonFilter");
  iEvent.put(std::move(passbadchargedcandidateFilter), "passBadChargedCandidateFilter");
  //
  iEvent.put(std::move(passchargedHadronTrackResolutionFilter),"passChargedHadronTrackResolutionFilter");
  iEvent.put(std::move(passmuonBadTrackFilter),"passMuonBadTrackFilter");
  iEvent.put(std::move(flagBadMuons),"badMuonsFlag");
  iEvent.put(std::move(flagDuplicateMuons),"duplicateMuonsFlag");
  iEvent.put(std::move(flagNoBadMuons),"noBadMuonsFlag");
}
