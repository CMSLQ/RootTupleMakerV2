#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Muons.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <iostream>
#include <DataFormats/TrackReco/interface/Track.h>
#include "TLorentzVector.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

RootTupleMakerV2_Muons::RootTupleMakerV2_Muons(const edm::ParameterSet& iConfig) :
inputTag (iConfig.getParameter<edm::InputTag>("InputTag")),
triggerEventInputTag     (iConfig.getParameter<edm::InputTag>("TriggerEventInputTag"     )),
prefix   (iConfig.getParameter<std::string>  ("Prefix")),
suffix   (iConfig.getParameter<std::string>  ("Suffix")),
maxSize  (iConfig.getParameter<unsigned int> ("MaxSize")),
muonIso  (iConfig.getParameter<double>       ("MuonIso")),                                 // threshold for "passIso" and "passCTIso"
muonID   (iConfig.getParameter<std::string>  ("MuonID")),                                  // threshold for "passID"
singleMuonTriggerMatch(iConfig.getParameter<std::string>("SingleMuonTriggerMatch")),       // trigger matching string
singleIsoMuonTriggerMatch(iConfig.getParameter<std::string>("SingleIsoMuonTriggerMatch")), // trigger matching string
beamSpotCorr      (iConfig.getParameter<bool>("BeamSpotCorr")),
useCocktailRefits (iConfig.getParameter<bool>("UseCocktailRefits")),
vtxInputTag       (iConfig.getParameter<edm::InputTag>("VertexInputTag"))                  // collection of primary vertices to be used.
{
  produces <bool>                 ( "hasVeryForwardPFMuon" );
  produces <std::vector<double> > ( prefix + "Eta"                     + suffix );
  produces <std::vector<double> > ( prefix + "Phi"                     + suffix );
  produces <std::vector<double> > ( prefix + "Pt"                      + suffix );
  produces <std::vector<double> > ( prefix + "EtaError"                + suffix );
  produces <std::vector<double> > ( prefix + "PhiError"                + suffix );
  produces <std::vector<double> > ( prefix + "PtError"                 + suffix );
  produces <std::vector<double> > ( prefix + "TrkEta"                  + suffix );
  produces <std::vector<double> > ( prefix + "TrkPhi"                  + suffix );
  produces <std::vector<double> > ( prefix + "TrkPt"                   + suffix );
  produces <std::vector<double> > ( prefix + "TrkEtaError"             + suffix );
  produces <std::vector<double> > ( prefix + "TrkPhiError"             + suffix );
  produces <std::vector<double> > ( prefix + "TrkPtError"              + suffix );
  produces <std::vector<double> > ( prefix + "QOverPError"             + suffix );
  produces <std::vector<double> > ( prefix + "P"                       + suffix );
  produces <std::vector<double> > ( prefix + "Energy"                  + suffix );
  produces <std::vector<int> >    ( prefix + "Charge"                  + suffix );
  produces <std::vector<int> >    ( prefix + "TrkHits"                 + suffix );
  produces <std::vector<int> >    ( prefix + "TrkHitsTrackerOnly"      + suffix );
  produces <std::vector<int> >    ( prefix + "GlobalTrkValidHits"      + suffix );
  produces <std::vector<int> >    ( prefix + "PixelHits"               + suffix );
  produces <std::vector<int> >    ( prefix + "TrkPixelHits"            + suffix );
  produces <std::vector<int> >    ( prefix + "SegmentMatches"          + suffix );
  produces <std::vector<int> >    ( prefix + "StationMatches"          + suffix );
  produces <std::vector<double> > ( prefix + "TrkValidFractionOfHits"  + suffix );
  produces <std::vector<double> > ( prefix + "TrkD0"                   + suffix );
  produces <std::vector<double> > ( prefix + "TrkD0Error"              + suffix );
  produces <std::vector<double> > ( prefix + "TrkDz"                   + suffix );
  produces <std::vector<double> > ( prefix + "TrkDzError"              + suffix );
  produces <std::vector<double> > ( prefix + "TrkVx"                   + suffix );
  produces <std::vector<double> > ( prefix + "TrkVy"                   + suffix );
  produces <std::vector<double> > ( prefix + "TrkVz"                   + suffix );
  produces <std::vector<double> > ( prefix + "TrackChi2"               + suffix );
  produces <std::vector<double> > ( prefix + "GlobalChi2"              + suffix );
  produces <std::vector<double> > ( prefix + "TrkIso"                  + suffix );
  produces <std::vector<double> > ( prefix + "TrackerIsoSumPT"         + suffix );
  produces <std::vector<double> > ( prefix + "EcalIso"                 + suffix );
  produces <std::vector<double> > ( prefix + "HcalIso"                 + suffix );
  produces <std::vector<double> > ( prefix + "HOIso"                   + suffix );
  produces <std::vector<double> > ( prefix + "EcalVetoIso"             + suffix );
  produces <std::vector<double> > ( prefix + "HcalVetoIso"             + suffix );
  produces <std::vector<double> > ( prefix + "PFIsoR03ChargedHadron"   + suffix );
  produces <std::vector<double> > ( prefix + "PFIsoR03ChargedParticle" + suffix );
  produces <std::vector<double> > ( prefix + "PFIsoR03NeutralHadron"   + suffix );
  produces <std::vector<double> > ( prefix + "PFIsoR03Photon"          + suffix );
  produces <std::vector<double> > ( prefix + "PFIsoR03NeutralHadronHT" + suffix );
  produces <std::vector<double> > ( prefix + "PFIsoR03PhotonHT"        + suffix );
  produces <std::vector<double> > ( prefix + "PFIsoR03PU"              + suffix );
  produces <std::vector<double> > ( prefix + "PFIsoR04ChargedHadron"   + suffix );
  produces <std::vector<double> > ( prefix + "PFIsoR04ChargedParticle" + suffix );
  produces <std::vector<double> > ( prefix + "PFIsoR04NeutralHadron"   + suffix );
  produces <std::vector<double> > ( prefix + "PFIsoR04Photon"          + suffix );
  produces <std::vector<double> > ( prefix + "PFIsoR04NeutralHadronHT" + suffix );
  produces <std::vector<double> > ( prefix + "PFIsoR04PhotonHT"        + suffix );
  produces <std::vector<double> > ( prefix + "PFIsoR04PU"              + suffix );
  produces <std::vector<int> >    ( prefix + "PassID"                  + suffix );
  produces <std::vector<int> >    ( prefix + "VtxIndex"                + suffix );
  produces <std::vector<double> > ( prefix + "VtxDistXY"               + suffix );
  produces <std::vector<double> > ( prefix + "VtxDistZ"                + suffix );
  produces <std::vector<int> >    ( prefix + "BestTrackVtxIndex"      + suffix );
  produces <std::vector<double> > ( prefix + "BestTrackVtxDistXY"     + suffix );
  produces <std::vector<double> > ( prefix + "BestTrackVtxDistZ"      + suffix );
  produces <std::vector<double> > ( prefix + "PrimaryVertexDXY"        + suffix );
  produces <std::vector<double> > ( prefix + "PrimaryVertexDXYError"   + suffix );
  produces <std::vector<double> > ( prefix + "BeamSpotDXY"             + suffix );
  produces <std::vector<double> > ( prefix + "BeamSpotDXYError"        + suffix );
  produces <std::vector<int> >    ( prefix + "IsGlobal"                + suffix );
  produces <std::vector<int> >    ( prefix + "IsTracker"               + suffix );
  produces <std::vector<double> > ( prefix + "MatchedGenParticlePt"    + suffix );
  produces <std::vector<double> > ( prefix + "MatchedGenParticleEta"   + suffix );
  produces <std::vector<double> > ( prefix + "MatchedGenParticlePhi"   + suffix );
  //
  // New variables added based on CMSSW 52X recommendations for LooseMuon and TightMuon Definitions
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Basline_muon_selections_for_2012
  produces <std::vector<int> >     ( prefix + "IsPF"                       + suffix );
  produces <std::vector<int> >     ( prefix + "TrackLayersWithMeasurement" + suffix );   

  // Variables for trigger matching  
  //  HLT Single Muon
  produces <std::vector<bool  > > ( prefix + "HLTSingleMuonMatched"      + suffix );
  produces <std::vector<double> > ( prefix + "HLTSingleMuonMatchPt"      + suffix );
  produces <std::vector<double> > ( prefix + "HLTSingleMuonMatchEta"     + suffix );
  produces <std::vector<double> > ( prefix + "HLTSingleMuonMatchPhi"     + suffix );
  //  HLT Single Iso Muon
  produces <std::vector<bool  > > ( prefix + "HLTSingleIsoMuonMatched"   + suffix );
  produces <std::vector<double> > ( prefix + "HLTSingleIsoMuonMatchPt"   + suffix );
  produces <std::vector<double> > ( prefix + "HLTSingleIsoMuonMatchEta"  + suffix );
  produces <std::vector<double> > ( prefix + "HLTSingleIsoMuonMatchPhi"  + suffix );
  
  //
  if ( useCocktailRefits )
    {
      produces <std::vector<int>    > ( prefix + "CocktailRefitID"                + suffix ) ;
      produces <std::vector<double> > ( prefix + "CocktailEta"                    + suffix ) ;
      produces <std::vector<double> > ( prefix + "CocktailPhi"                    + suffix ) ;
      produces <std::vector<double> > ( prefix + "CocktailPt"                     + suffix ) ;
      produces <std::vector<double> > ( prefix + "CocktailEtaError"               + suffix ) ;
      produces <std::vector<double> > ( prefix + "CocktailPhiError"               + suffix ) ;
      produces <std::vector<double> > ( prefix + "CocktailPtError"                + suffix ) ;
      produces <std::vector<double> > ( prefix + "CocktailQOverPError"            + suffix ) ;
      produces <std::vector<double> > ( prefix + "CocktailP"                      + suffix ) ;
      produces <std::vector<int   > > ( prefix + "CocktailCharge"                 + suffix ) ;
      produces <std::vector<int   > > ( prefix + "CocktailTrkHits"                + suffix ) ;
      produces <std::vector<double> > ( prefix + "CocktailTrkValidFractionOfHits" + suffix ) ;
      produces <std::vector<double> > ( prefix + "CocktailTrkD0"                  + suffix ) ;
      produces <std::vector<double> > ( prefix + "CocktailTrkD0Error"             + suffix ) ;
      produces <std::vector<double> > ( prefix + "CocktailTrkDz"                  + suffix ) ;
      produces <std::vector<double> > ( prefix + "CocktailTrkDzError"             + suffix ) ;
      produces <std::vector<double> > ( prefix + "CocktailGlobalChi2"             + suffix ) ;
    }
  
  produces <std::vector<double> > ( prefix + "CosmicCompatibility"     + suffix );
  produces <std::vector<double> > ( prefix + "TimeCompatibility"       + suffix );
  produces <std::vector<double> > ( prefix + "BackToBackCompatibility" + suffix );
  produces <std::vector<double> > ( prefix + "OverlapCompatibility"    + suffix );
}


void RootTupleMakerV2_Muons::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::auto_ptr<bool>                  hasVeryForwardPFMuon    ( new bool() );
  std::auto_ptr<std::vector<double> >  eta                     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi                     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt                      ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  etaError                ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phiError                ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ptError                 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkEta                  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkPhi                  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkPt                   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkEtaError             ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkPhiError             ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkPtError              ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  qoverpError             ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  p                       ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy                  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     charge                  ( new std::vector<int>()     );
  std::auto_ptr<std::vector<int> >     trkHits                 ( new std::vector<int>()     );
  std::auto_ptr<std::vector<int> >     trkHitsTrackerOnly      ( new std::vector<int>()     );
  std::auto_ptr<std::vector<int> >     GlobaltrkValidHits      ( new std::vector<int>()     );
  std::auto_ptr<std::vector<int> >     pixelHits               ( new std::vector<int>()     );
  std::auto_ptr<std::vector<int> >     trkPixelHits            ( new std::vector<int>()     );
  std::auto_ptr<std::vector<int> >     segmentMatches          ( new std::vector<int>()     );
  std::auto_ptr<std::vector<int> >     stationMatches          ( new std::vector<int>()     );
  std::auto_ptr<std::vector<int> >     Valid                   ( new std::vector<int>()     );
  std::auto_ptr<std::vector<double> >  trkValidFractionOfHits  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkD0                   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkD0Error              ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkDz                   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkDzError              ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkVx                   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkVy                   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkVz                   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trackChi2               ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  globalChi2              ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkIso                  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trackerIsoSumPT         ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ecalIso                 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIso                 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hoIso                   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ecalVetoIso             ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalVetoIso             ( new std::vector<double>()  );
  //
  std::auto_ptr<std::vector<double> >  pfisor03chargedhadron   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfisor03chargedparticle ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfisor03neutralhadron   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfisor03photon          ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfisor03neutralhadronht ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfisor03photonht        ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfisor03pu              ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfisor04chargedhadron   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfisor04chargedparticle ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfisor04neutralhadron   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfisor04photon          ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfisor04neutralhadronht ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfisor04photonht        ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfisor04pu              ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     passID                  ( new std::vector<int>()     );
  std::auto_ptr<std::vector<int> >     vtxIndex                ( new std::vector<int>()     );
  std::auto_ptr<std::vector<double> >  vtxDistXY               ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  vtxDistZ                ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     bestTrackVtxIndex       ( new std::vector<int>()     );
  std::auto_ptr<std::vector<double> >  bestTrackVtxDistXY      ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  bestTrackVtxDistZ       ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  primaryVertexDXY        ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  primaryVertexDXYError   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  beamspotDXY             ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  beamspotDXYError        ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     IsGlobal                ( new std::vector<int>()     );
  std::auto_ptr<std::vector<int> >     IsTracker               ( new std::vector<int>()     );
  std::auto_ptr<std::vector<double> >  matchedgenparticlept    ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  matchedgenparticleeta   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  matchedgenparticlephi   ( new std::vector<double>()  );
  //
  // New variables added based on CMSSW 52X recommendations for LooseMuon and TightMuon Definitions
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Basline_muon_selections_for_2012 
  std::auto_ptr<std::vector<int> >     isPF                       ( new std::vector<int>()    );
  std::auto_ptr<std::vector<int> >     trackLayersWithMeasurement ( new std::vector<int>()    );   
  //
  // Trigger matching variables  
  std::auto_ptr<std::vector<bool  > >  HLTSingleMuonMatched     ( new std::vector<bool  >()  );
  std::auto_ptr<std::vector<double> >  HLTSingleMuonMatchPt 	( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  HLTSingleMuonMatchEta	( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  HLTSingleMuonMatchPhi    ( new std::vector<double>()  );
  std::auto_ptr<std::vector<bool  > >  HLTSingleIsoMuonMatched  ( new std::vector<bool  >()  );
  std::auto_ptr<std::vector<double> >  HLTSingleIsoMuonMatchPt 	( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  HLTSingleIsoMuonMatchEta	( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  HLTSingleIsoMuonMatchPhi ( new std::vector<double>()  );
  //
  std::auto_ptr<std::vector<int   > >  ctRefitID    ( new std::vector<int   > () );
  std::auto_ptr<std::vector<double> >  ctEta        ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  ctPhi        ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  ctPt         ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  ctEtaError   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ctPhiError   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ctPtError    ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ctQoverpError( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ctP          ( new std::vector<double> () );
  std::auto_ptr<std::vector<int   > >  ctCharge     ( new std::vector<int   > () );
  std::auto_ptr<std::vector<int   > >  ctTrkHits    ( new std::vector<int   > () );
  std::auto_ptr<std::vector<double> >  ctTrkValidFractionOfHits ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  ctTrkD0      ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  ctTrkD0Error ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  ctTrkDz      ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  ctTrkDzError ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  ctGlobalChi2 ( new std::vector<double> () );
  //
  std::auto_ptr<std::vector<double> >  cosmicCompatibility     ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  timeCompatibility       ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  backToBackCompatibility ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  overlapCompatibility    ( new std::vector<double> () );
  
  //-----------------------------------------------------------------
  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByLabel(inputTag, muons);
  
  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(vtxInputTag,primaryVertices);
  
  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByLabel("offlineBeamSpot", beamSpot );

  edm::Handle< pat::TriggerEvent > triggerEvent;
  iEvent.getByLabel( triggerEventInputTag, triggerEvent );
  
  // PAT trigger helper for trigger matching information
  const pat::helper::TriggerMatchHelper matchHelper;
  
  *hasVeryForwardPFMuon.get() = false;

  if(muons.isValid())
    {
      edm::LogInfo("RootTupleMakerV2_MuonsInfo") << "Total # Muons: " << muons->size();
      
      size_t iMuon = 0;

      for( std::vector<pat::Muon>::const_iterator it = muons->begin(); it != muons->end(); ++it )
	{
	  // exit from loop when you reach the required number of muons
	  if(eta->size() >= maxSize)
	    break;

	  if( it->isPFMuon() && fabs(it->eta())>2.2 ) *hasVeryForwardPFMuon.get() = true;
	  
	  //if muon is neither global nor tracker muon, continue.
	  if( !it->isGlobalMuon() && !it->isTrackerMuon() ) continue;
	  // if muonPt is less than 5GeV, continue.
	  if( it->pt()<5 ) continue;

	  /// Gen Matching
	  double genparPt = -999.;
	  double genparEta= -999.;
	  double genparPhi= -999.;
	  if ( !iEvent.isRealData() ) {
	    for(uint igen = 0 ; igen < it->genParticleRefs().size() ; ++igen ){//it->genParticleRefs().size() should be 0 or 1
	      genparPt =it->genParticle(igen)->pt();
	      genparEta=it->genParticle(igen)->eta();
	      genparPhi=it->genParticle(igen)->phi();
	    }
	  }
	  matchedgenparticlept     -> push_back ( (double)(genparPt) );
	  matchedgenparticleeta    -> push_back ( (double)(genparEta) );
	  matchedgenparticlephi    -> push_back ( (double)(genparPhi) );
	  /// Gen Matching
	  
	  double trkd0   = it->track()->d0();
	  
	  if( beamSpotCorr && beamSpot.isValid() )
	    {
	      trkd0   = -(it->track()   ->dxy( beamSpot->position()));
	    }
	  
	  else if( beamSpotCorr && !beamSpot.isValid() ) edm::LogError("RootTupleMakerV2_MuonsError") << "Error! Can't get the offlineBeamSpot";
	  
	  // Vertex association 
	  //-- using pat::muon::track()
	  double minVtxDist3D = 9999.;
	  int    vtxIndex_    = -1;
	  double vtxDistXY_   = -9999.;
	  double vtxDistZ_    = -9999.;
	  //-- using reco::muon::muonBestTrack() - 2012 recommendation  
	  double bt_minVtxDist3D = 9999.;
	  int    bt_vtxIndex_    = -1;
	  double bt_vtxDistXY_   = -9999.;
	  double bt_vtxDistZ_    = -9999.;

	  if(primaryVertices.isValid())
	    {
	      edm::LogInfo("RootTupleMakerV2_MuonsInfo") << "Total # Primary Vertices: " << primaryVertices->size();
	      
	      for( reco::VertexCollection::const_iterator v_it=primaryVertices->begin() ; v_it!=primaryVertices->end() ; ++v_it )
		{
		  //-- using pat::muon::track()
		  double distXY = it->track()->dxy(v_it->position());
		  double distZ  = it->track()->dz(v_it->position());
		  double dist3D = sqrt(pow(distXY,2) + pow(distZ,2));
		  if( dist3D < minVtxDist3D ){
		    minVtxDist3D = dist3D;
		    vtxIndex_    = int(std::distance(primaryVertices->begin(),v_it));
		    vtxDistXY_   = distXY;
		    vtxDistZ_    = distZ;
		  }
		  //-- using reco::muon::muonBestTrack() - 2012 recommendation    
		  if( (it->muonBestTrack()).isNonnull() ){
		    double bt_distXY = it->muonBestTrack()->dxy(v_it->position());
		    double bt_distZ  = it->muonBestTrack()->dz(v_it->position());
		    double bt_dist3D = sqrt( pow(bt_distXY,2) + pow(bt_distZ,2) );
		    if( bt_dist3D < bt_minVtxDist3D ){
		      bt_minVtxDist3D = bt_dist3D;
		      bt_vtxIndex_    = int(std::distance(primaryVertices->begin(),v_it));
		      bt_vtxDistXY_   = bt_distXY;
		      bt_vtxDistZ_    = bt_distZ;
		    }
		  }
		  
		}//loop over primaryVertices
	    }
	  else
	    {
	      edm::LogError("RootTupleMakerV2_MuonsError") << "Error! Can't get the product " << vtxInputTag;
	    }

	  //
	  //------------------------------------------------------------------------
	  // Trigger matching
	  // 
	  // Example taken from PatTriggerAnalyzer:
	  // http://cmslxr.fnal.gov/lxr/source/PhysicsTools/PatExamples/plugins/PatTriggerAnalyzer.cc
	  //------------------------------------------------------------------------
	  
	  // HLT Single Muon trigger matching	  
	  const pat::TriggerObjectRef singleMuonTrigRef( matchHelper.triggerMatchObject( muons, iMuon,  singleMuonTriggerMatch, iEvent, *triggerEvent ) );
	  if ( singleMuonTrigRef.isAvailable() && singleMuonTrigRef.isNonnull() ) { 
	    HLTSingleMuonMatched  -> push_back ( true ) ;
	    HLTSingleMuonMatchPt  -> push_back ( singleMuonTrigRef -> pt() );
	    HLTSingleMuonMatchEta -> push_back ( singleMuonTrigRef -> eta());
	    HLTSingleMuonMatchPhi -> push_back ( singleMuonTrigRef -> phi());
	  } else { 
	    HLTSingleMuonMatched  -> push_back ( false ) ;
	    HLTSingleMuonMatchPt  -> push_back ( -999. );
	    HLTSingleMuonMatchEta -> push_back ( -999. );
	    HLTSingleMuonMatchPhi -> push_back ( -999. );
	  }

	  // HLT Single Iso Muon trigger matching	  
	  const pat::TriggerObjectRef singleIsoMuonTrigRef( matchHelper.triggerMatchObject( muons, iMuon,  singleIsoMuonTriggerMatch, iEvent, *triggerEvent ) );
	  if ( singleIsoMuonTrigRef.isAvailable() && singleIsoMuonTrigRef.isNonnull() ) { 
	    HLTSingleIsoMuonMatched  -> push_back ( true ) ;
	    HLTSingleIsoMuonMatchPt  -> push_back ( singleIsoMuonTrigRef -> pt() );
	    HLTSingleIsoMuonMatchEta -> push_back ( singleIsoMuonTrigRef -> eta());
	    HLTSingleIsoMuonMatchPhi -> push_back ( singleIsoMuonTrigRef -> phi());
	  } else { 
	    HLTSingleIsoMuonMatched  -> push_back ( false ) ;
	    HLTSingleIsoMuonMatchPt  -> push_back ( -999. );
	    HLTSingleIsoMuonMatchEta -> push_back ( -999. );
	    HLTSingleIsoMuonMatchPhi -> push_back ( -999. );
	  }

	  eta->push_back( it->eta() );
	  phi->push_back( it->phi() );
	  pt ->push_back( it->pt()  );
	  p  ->push_back( it->p()   );
	  
	  if( it->isGlobalMuon() )
	    {
	      etaError    -> push_back ( it->globalTrack()->etaError()    );
	      phiError    -> push_back ( it->globalTrack()->phiError()    );
	      ptError     -> push_back ( it->globalTrack()->ptError ()    );
	      qoverpError -> push_back ( it->globalTrack()->qoverpError() );
	    }
	  else
	    {
	      etaError    -> push_back ( it->track()->etaError()    ); //track() returns innerTrack();
	      phiError    -> push_back ( it->track()->phiError()    ); //track() returns innerTrack();
	      ptError     -> push_back ( it->track()->ptError ()    ); //track() returns innerTrack();
	      qoverpError -> push_back ( it->track()->qoverpError() ); //track() returns innerTrack();
	    }
	  
	  trkPt  -> push_back ( it->track()->pt()  );
	  trkEta -> push_back ( it->track()->eta() );
	  trkPhi -> push_back ( it->track()->phi() );
	  
	  trkPtError  -> push_back ( it->track()->ptError()  );
	  trkEtaError -> push_back ( it->track()->etaError() );
	  trkPhiError -> push_back ( it->track()->phiError() );
	  
	  charge            ->push_back( it->charge() );
	  trkHits           ->push_back( it->track()->numberOfValidHits() );
	  trkHitsTrackerOnly->push_back( it->track()->hitPattern().numberOfValidTrackerHits() );

	  if( it->isGlobalMuon() )
	    {
	      GlobaltrkValidHits->push_back( it->globalTrack()->hitPattern().numberOfValidMuonHits()  );
	      pixelHits         ->push_back( it->globalTrack()->hitPattern().numberOfValidPixelHits() );	    
	      globalChi2        ->push_back( it->globalTrack()->normalizedChi2()                      );
	    }
	  else
	    {
	      GlobaltrkValidHits->push_back( -1 );//inner track does not have muon hits by definition.
	      pixelHits         ->push_back( -1 );//inner track pixel hit count is stored at "innerTrackPixelHits".
	      globalChi2        ->push_back( -1 );//inner track Chi squared is already stored at  "trackChi2".
	    }

	  trkPixelHits->push_back(it->track()->hitPattern().numberOfValidPixelHits());
       
	  segmentMatches  ->push_back(it->numberOfMatches());
	  stationMatches  ->push_back(it->numberOfMatchedStations());
	  
	  trkValidFractionOfHits->push_back( it->track()->validFraction() );

	  trkD0     ->push_back( trkd0                         );
	  trkD0Error->push_back( it->track()->d0Error()        );
	  trkDz     ->push_back( it->track()->dz()             );
	  trkDzError->push_back( it->track()->dzError()        );
	  trkVx     ->push_back( it->track()->vx()             );
	  trkVy     ->push_back( it->track()->vy()             );
	  trkVz     ->push_back( it->track()->vz()             );
	  trackChi2 ->push_back( it->track()->normalizedChi2() );

	  // Global High Pt Muons, aka Cocktail Muons
	  if ( useCocktailRefits && it->isGlobalMuon() )
	    {	      
	      int refit_id = -999;

	      // -----------------------------------------------------------------------------------------------------------------------------------------//
	      // HighPT Muon - New Version (recommended) 
	      // Following new instructions on https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#New_Version_recommended
	      // To use the new version of TuneP you have to follow these steps:
	      // Either run CMSSW_5_3_6_patch1 (or newer) or if you need to keep using an older CMSSW verison check out V09-04-03-02 DataFormats/MuonReco 
	      // (if you have problems getting CMSSW to compile please do also addpkg RecoMuon/MuonIdentification).
	      // #include "DataFormats/MuonReco/interface/MuonCocktails.h" add this to your analysis code
	      // reco::TrackRef cktTrack = (muon::tevOptimized(*recoMu, 200, 17., 40., 0.25)).first; call to get the optimal muon track
	      // cktTrack->pt() - to get the pT of the muon
	      //
	      //Then you can apply the new HighPT ID, which still differs from Tight Muon selection in the following points:
	      //The Particle-Flow muon id is not required
	      //The cut of dpT/pT<0.3 for the track used for momentum determination is applied, i.e. cktTrack->ptError()/cktTrack->pt()<0.3
	      //The cuts applied on recoMu.muonBestTrack (impact parameter cuts) need to be applied to cktTrack since this is the new best track now. 
	      // -----------------------------------------------------------------------------------------------------------------------------------------//

	      reco::TrackRef cocktail_track = (muon::tevOptimized(*it, 200, 17., 40., 0.25)).first;

	      double cttrkd0  = cocktail_track -> d0() ;
	      if( beamSpotCorr && beamSpot.isValid() )
		cttrkd0 = -(cocktail_track->dxy( beamSpot->position()));
	      
	      ctRefitID     -> push_back ( refit_id ) ;
	      
	      ctEtaError    -> push_back ( cocktail_track->etaError()    );
	      ctPhiError    -> push_back ( cocktail_track->phiError()    );
	      ctPtError     -> push_back ( cocktail_track->ptError ()    );
	      ctQoverpError -> push_back ( cocktail_track->qoverpError() );
	      
	      ctEta                    ->push_back( cocktail_track->eta()    );
	      ctPhi                    ->push_back( cocktail_track->phi()    );
	      ctPt                     ->push_back( cocktail_track->pt()     );
	      ctP                      ->push_back( cocktail_track->p()      );
	      ctCharge                 ->push_back( cocktail_track->charge() );
	      ctTrkHits                ->push_back( cocktail_track->hitPattern().numberOfValidTrackerHits() );
	      ctTrkValidFractionOfHits ->push_back( cocktail_track->validFraction()   );
	      ctTrkD0                  ->push_back( cttrkd0                           );
	      ctTrkD0Error             ->push_back( cocktail_track->d0Error()         );
	      ctTrkDz                  ->push_back( cocktail_track->dz()              );
	      ctTrkDzError             ->push_back( cocktail_track -> dzError()       );
	      ctGlobalChi2             ->push_back( cocktail_track ->normalizedChi2() );
	    }//cocktail fits

	  // New variables added based on CMSSW 52X recommendations for LooseMuon and TightMuon Definitions
	  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Basline_muon_selections_for_2012 	  
	  isPF                       ->push_back( it->isPFMuon()       );
	  trackLayersWithMeasurement ->push_back( it->track()->hitPattern().trackerLayersWithMeasurement() );
	  
	  energy->push_back( it->energy() );
	  //
	  // 2011 isolation parameters..
	  // See: http://cmslxr.fnal.gov/lxr/source/DataFormats/MuonReco/interface/MuonIsolation.h?v=CMSSW_5_3_3_cand1#019
	  trkIso         ->push_back( it->trackIso()           );
	  trackerIsoSumPT->push_back( it->isolationR03().sumPt );
	  ecalIso        ->push_back( it->isolationR03().emEt  );
	  hcalIso        ->push_back( it->isolationR03().hadEt );
	  hoIso          ->push_back( it->isolationR03().hoEt  );
	  //
	  ecalVetoIso    ->push_back( it->isolationR03().emVetoEt  );
	  hcalVetoIso    ->push_back( it->isolationR03().hadVetoEt );
	  //
	  // Adding PF isolation for 2012..
	  //https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
	  //The PF based isolation can be accessed by the reco::Muon using the following methods
	  //const MuonPFIsolation& pfIsolationR03() /// Cone of 0.3
	  //const MuonPFIsolation& pfIsolationR04() /// Cone of 0.4-Suggested  
	  //http://cmssdt.cern.ch/SDT/doxygen/CMSSW_5_2_5_patch3/doc/html/df/d33/structreco_1_1MuonPFIsolation.html
	  // "Inner Veto Cone" and "Iso Particle Pt Threshold" descriptions are here: 
	  // http://cmslxr.fnal.gov/lxr/source/RecoMuon/MuonIsolation/python/muonPFIsolationValues_cff.py?v=CMSSW_5_2_4
	  pfisor03chargedhadron  ->push_back( it->pfIsolationR03().sumChargedHadronPt              );
	  pfisor03chargedparticle->push_back( it->pfIsolationR03().sumChargedParticlePt            );
          pfisor03neutralhadron  ->push_back( it->pfIsolationR03().sumNeutralHadronEt              );
          pfisor03photon         ->push_back( it->pfIsolationR03().sumPhotonEt                     );
          pfisor03neutralhadronht->push_back( it->pfIsolationR03().sumNeutralHadronEtHighThreshold );
          pfisor03photonht       ->push_back( it->pfIsolationR03().sumPhotonEtHighThreshold        );
          pfisor03pu             ->push_back( it->pfIsolationR03().sumPUPt                         );
	  pfisor04chargedhadron  ->push_back( it->pfIsolationR04().sumChargedHadronPt              );
	  pfisor04chargedparticle->push_back( it->pfIsolationR04().sumChargedParticlePt            );
          pfisor04neutralhadron  ->push_back( it->pfIsolationR04().sumNeutralHadronEt              );
          pfisor04photon         ->push_back( it->pfIsolationR04().sumPhotonEt                     );
          pfisor04neutralhadronht->push_back( it->pfIsolationR04().sumNeutralHadronEtHighThreshold );
          pfisor04photonht       ->push_back( it->pfIsolationR04().sumPhotonEtHighThreshold        );
          pfisor04pu             ->push_back( it->pfIsolationR04().sumPUPt                         );

	  passID                 ->push_back( (it->muonID(muonID)) ?  1 : 0 );
	  IsGlobal               ->push_back( (it->isGlobalMuon()) ?  1 : 0 );
	  IsTracker              ->push_back( (it->isTrackerMuon()) ? 1 : 0 );
	  vtxIndex               ->push_back( vtxIndex_                     );
	  vtxDistXY              ->push_back( vtxDistXY_                    );
	  vtxDistZ               ->push_back( vtxDistZ_                     );
	  bestTrackVtxIndex      ->push_back( bt_vtxIndex_                  );
	  bestTrackVtxDistXY     ->push_back( bt_vtxDistXY_                 );
	  bestTrackVtxDistZ      ->push_back( bt_vtxDistZ_                  );
	  primaryVertexDXY       ->push_back( it->dB()                      );
	  primaryVertexDXYError  ->push_back( it->edB()                     );
	  beamspotDXY            ->push_back( it->dB(pat::Muon::BS2D)       );
	  beamspotDXYError       ->push_back( it->edB(pat::Muon::BS2D)      );
	  
	  // See https://indico.cern.ch/getFile.py/access?contribId=7&resId=0&materialId=slides&confId=102306
	  // and https://indico.cern.ch/getFile.py/access?contribId=5&resId=0&materialId=slides&confId=128840
	  cosmicCompatibility    ->push_back( it->userFloat("cosmicCompatibility")     );
	  timeCompatibility      ->push_back( it->userFloat("timeCompatibility")       );
	  backToBackCompatibility->push_back( it->userFloat("backToBackCompatibility") );
	  overlapCompatibility   ->push_back( it->userFloat("overlapCompatibility")    );
	  
	  iMuon++;

	}//end of loop over muons
    }
  else
    {
      edm::LogError("RootTupleMakerV2_MuonsError") << "Error! Can't get the product " << inputTag;
    }
  
  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( hasVeryForwardPFMuon, "hasVeryForwardPFMuon" );
  iEvent.put( eta,                        prefix + "Eta"                         + suffix );
  iEvent.put( phi,                        prefix + "Phi"                         + suffix );
  iEvent.put( pt,                         prefix + "Pt"                          + suffix );
  iEvent.put( etaError,                   prefix + "EtaError"                    + suffix );
  iEvent.put( phiError,                   prefix + "PhiError"                    + suffix );
  iEvent.put( ptError,                    prefix + "PtError"                     + suffix );
  iEvent.put( trkEta,                     prefix + "TrkEta"                      + suffix );
  iEvent.put( trkPhi,                     prefix + "TrkPhi"                      + suffix );
  iEvent.put( trkPt,                      prefix + "TrkPt"                       + suffix );
  iEvent.put( trkEtaError,                prefix + "TrkEtaError"                 + suffix );
  iEvent.put( trkPhiError,                prefix + "TrkPhiError"                 + suffix );
  iEvent.put( trkPtError,                 prefix + "TrkPtError"                  + suffix );
  iEvent.put( qoverpError,                prefix + "QOverPError"                 + suffix );
  iEvent.put( p,                          prefix + "P"                           + suffix );
  iEvent.put( energy,                     prefix + "Energy"                      + suffix );
  iEvent.put( charge,                     prefix + "Charge"                      + suffix );
  iEvent.put( trkHits,                    prefix + "TrkHits"                     + suffix );
  iEvent.put( trkHitsTrackerOnly,         prefix + "TrkHitsTrackerOnly"          + suffix );
  iEvent.put( GlobaltrkValidHits,         prefix + "GlobalTrkValidHits"          + suffix );
  iEvent.put( pixelHits,                  prefix + "PixelHits"                   + suffix );
  iEvent.put( trkPixelHits,               prefix + "TrkPixelHits"                + suffix );
  iEvent.put( segmentMatches,             prefix + "SegmentMatches"              + suffix );
  iEvent.put( stationMatches,             prefix + "StationMatches"              + suffix );
  iEvent.put( trkValidFractionOfHits,     prefix + "TrkValidFractionOfHits"      + suffix );
  iEvent.put( trkD0,                      prefix + "TrkD0"                       + suffix );
  iEvent.put( trkD0Error,                 prefix + "TrkD0Error"                  + suffix );
  iEvent.put( trkDz,                      prefix + "TrkDz"                       + suffix );
  iEvent.put( trkDzError,                 prefix + "TrkDzError"                  + suffix );
  iEvent.put( trkVx,                      prefix + "TrkVx"                       + suffix );
  iEvent.put( trkVy,                      prefix + "TrkVy"                       + suffix );
  iEvent.put( trkVz,                      prefix + "TrkVz"                       + suffix );
  iEvent.put( trackChi2,                  prefix + "TrackChi2"                   + suffix );
  iEvent.put( globalChi2,                 prefix + "GlobalChi2"                  + suffix );
  iEvent.put( trkIso,                     prefix + "TrkIso"                      + suffix );
  iEvent.put( trackerIsoSumPT,            prefix + "TrackerIsoSumPT"             + suffix );
  iEvent.put( ecalIso,                    prefix + "EcalIso"                     + suffix );
  iEvent.put( hcalIso,                    prefix + "HcalIso"                     + suffix );
  iEvent.put( hoIso,                      prefix + "HOIso"                       + suffix );
  iEvent.put( ecalVetoIso,                prefix + "EcalVetoIso"                 + suffix );
  iEvent.put( hcalVetoIso,                prefix + "HcalVetoIso"                 + suffix );
  iEvent.put( pfisor03chargedhadron,      prefix + "PFIsoR03ChargedHadron"       + suffix );
  iEvent.put( pfisor03chargedparticle,    prefix + "PFIsoR03ChargedParticle"     + suffix );
  iEvent.put( pfisor03neutralhadron,      prefix + "PFIsoR03NeutralHadron"       + suffix );
  iEvent.put( pfisor03photon,             prefix + "PFIsoR03Photon"              + suffix );
  iEvent.put( pfisor03neutralhadronht,    prefix + "PFIsoR03NeutralHadronHT"     + suffix );
  iEvent.put( pfisor03photonht,           prefix + "PFIsoR03PhotonHT"            + suffix );
  iEvent.put( pfisor03pu,                 prefix + "PFIsoR03PU"                  + suffix );
  iEvent.put( pfisor04chargedhadron,      prefix + "PFIsoR04ChargedHadron"       + suffix );
  iEvent.put( pfisor04chargedparticle,    prefix + "PFIsoR04ChargedParticle"     + suffix );
  iEvent.put( pfisor04neutralhadron,      prefix + "PFIsoR04NeutralHadron"       + suffix );
  iEvent.put( pfisor04photon,             prefix + "PFIsoR04Photon"              + suffix );
  iEvent.put( pfisor04neutralhadronht,    prefix + "PFIsoR04NeutralHadronHT"     + suffix );
  iEvent.put( pfisor04photonht,           prefix + "PFIsoR04PhotonHT"            + suffix );
  iEvent.put( pfisor04pu,                 prefix + "PFIsoR04PU"                  + suffix );
  iEvent.put( passID,                     prefix + "PassID"                      + suffix );
  iEvent.put( IsGlobal,                   prefix + "IsGlobal"                    + suffix );
  iEvent.put( IsTracker,                  prefix + "IsTracker"                   + suffix );
  iEvent.put( vtxIndex,                   prefix + "VtxIndex"                    + suffix );
  iEvent.put( vtxDistXY,                  prefix + "VtxDistXY"                   + suffix );
  iEvent.put( vtxDistZ,                   prefix + "VtxDistZ"                    + suffix );
  iEvent.put( bestTrackVtxIndex,          prefix + "BestTrackVtxIndex"           + suffix );
  iEvent.put( bestTrackVtxDistXY,         prefix + "BestTrackVtxDistXY"          + suffix );
  iEvent.put( bestTrackVtxDistZ,          prefix + "BestTrackVtxDistZ"           + suffix );
  iEvent.put( primaryVertexDXY,           prefix + "PrimaryVertexDXY"            + suffix );
  iEvent.put( primaryVertexDXYError,      prefix + "PrimaryVertexDXYError"       + suffix );
  iEvent.put( beamspotDXY,                prefix + "BeamSpotDXY"                 + suffix );
  iEvent.put( beamspotDXYError,           prefix + "BeamSpotDXYError"            + suffix );
  iEvent.put( matchedgenparticlept,       prefix + "MatchedGenParticlePt"        + suffix );
  iEvent.put( matchedgenparticleeta,      prefix + "MatchedGenParticleEta"       + suffix );
  iEvent.put( matchedgenparticlephi,      prefix + "MatchedGenParticlePhi"       + suffix );
  iEvent.put( isPF,                       prefix + "IsPF"                        + suffix );
  iEvent.put( trackLayersWithMeasurement, prefix + "TrackLayersWithMeasurement"  + suffix );
  //
  
  if ( useCocktailRefits )
    {
      iEvent.put( ctRefitID                , prefix + "CocktailRefitID"                 + suffix ) ;
      iEvent.put( ctEta                    , prefix + "CocktailEta"                     + suffix ) ;
      iEvent.put( ctPhi                    , prefix + "CocktailPhi"                     + suffix ) ;
      iEvent.put( ctPt                     , prefix + "CocktailPt"                      + suffix ) ;
      iEvent.put( ctEtaError               , prefix + "CocktailEtaError"                + suffix ) ;
      iEvent.put( ctPhiError               , prefix + "CocktailPhiError"                + suffix ) ;
      iEvent.put( ctPtError                , prefix + "CocktailPtError"                 + suffix ) ;
      iEvent.put( ctQoverpError            , prefix + "CocktailQOverPError"             + suffix ) ;
      iEvent.put( ctP                      , prefix + "CocktailP"                       + suffix ) ;
      iEvent.put( ctCharge                 , prefix + "CocktailCharge"                  + suffix ) ;
      iEvent.put( ctTrkHits                , prefix + "CocktailTrkHits"                 + suffix ) ;
      iEvent.put( ctTrkValidFractionOfHits , prefix + "CocktailTrkValidFractionOfHits"  + suffix ) ;
      iEvent.put( ctTrkD0                  , prefix + "CocktailTrkD0"                   + suffix ) ;
      iEvent.put( ctTrkD0Error             , prefix + "CocktailTrkD0Error"              + suffix ) ;
      iEvent.put( ctTrkDz                  , prefix + "CocktailTrkDz"                   + suffix ) ;
      iEvent.put( ctTrkDzError             , prefix + "CocktailTrkDzError"              + suffix ) ;
      iEvent.put( ctGlobalChi2             , prefix + "CocktailGlobalChi2"              + suffix ) ;
    }
  
  // trigger matching
  iEvent.put( HLTSingleMuonMatched     , prefix + "HLTSingleMuonMatched"      + suffix );
  iEvent.put( HLTSingleMuonMatchPt     , prefix + "HLTSingleMuonMatchPt"      + suffix );
  iEvent.put( HLTSingleMuonMatchEta    , prefix + "HLTSingleMuonMatchEta"     + suffix );
  iEvent.put( HLTSingleMuonMatchPhi    , prefix + "HLTSingleMuonMatchPhi"     + suffix );
  iEvent.put( HLTSingleIsoMuonMatched  , prefix + "HLTSingleIsoMuonMatched"   + suffix );
  iEvent.put( HLTSingleIsoMuonMatchPt  , prefix + "HLTSingleIsoMuonMatchPt"   + suffix );
  iEvent.put( HLTSingleIsoMuonMatchEta , prefix + "HLTSingleIsoMuonMatchEta"  + suffix );
  iEvent.put( HLTSingleIsoMuonMatchPhi , prefix + "HLTSingleIsoMuonMatchPhi"  + suffix );

  iEvent.put( cosmicCompatibility,     prefix + "CosmicCompatibility"     + suffix );
  iEvent.put( timeCompatibility,       prefix + "TimeCompatibility"       + suffix );
  iEvent.put( backToBackCompatibility, prefix + "BackToBackCompatibility" + suffix );
  iEvent.put( overlapCompatibility,    prefix + "OverlapCompatibility"    + suffix );
}
