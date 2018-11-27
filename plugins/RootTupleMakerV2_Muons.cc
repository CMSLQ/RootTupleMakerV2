#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_Muons.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include <iostream>
#include <DataFormats/TrackReco/interface/Track.h>
#include "TLorentzVector.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

RootTupleMakerV2_Muons::RootTupleMakerV2_Muons(const edm::ParameterSet& iConfig) :
  //inputTag (iConfig.getParameter<edm::InputTag>("InputTag")),
  muonInputToken_  (consumes<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("InputTag"))),
  //triggerEventInputTag     (iConfig.getParameter<edm::InputTag>("TriggerEventInputTag"     )),
  prefix   (iConfig.getParameter<std::string>  ("Prefix")),
  suffix   (iConfig.getParameter<std::string>  ("Suffix")),
  maxSize  (iConfig.getParameter<unsigned int> ("MaxSize")),
  // threshold for "passIso" and "passCTIso"
  muonIso  (iConfig.getParameter<double>       ("MuonIso")),
  // threshold for "passID"
  muonID   (iConfig.getParameter<std::string>  ("MuonID")),
  beamSpotCorr      (iConfig.getParameter<bool>("BeamSpotCorr")),
  useCocktailRefits (iConfig.getParameter<bool>("UseCocktailRefits")),
  // collection of primary vertices to be used.
  //vtxInputTag       (iConfig.getParameter<edm::InputTag>("VertexInputTag"))
  vtxInputToken_ (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("VertexInputTag"))),
  beamSpotToken_ (consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpotInputTag")))
{
  produces <bool>                 ( "hasVeryForwardPFMuon" );
  produces <std::vector<bool> >   ( prefix + "IsLooseMuon"             + suffix );
  produces <std::vector<bool> >   ( prefix + "IsMediumMuon"            + suffix );
  produces <std::vector<bool> >   ( prefix + "IsTightMuon"             + suffix );
  produces <std::vector<bool> >   ( prefix + "IsHighPtMuon"            + suffix );
  produces <std::vector<float> >  ( prefix + "Eta"                     + suffix );
  produces <std::vector<float> >  ( prefix + "Phi"                     + suffix );
  produces <std::vector<float> >  ( prefix + "Pt"                      + suffix );
  produces <std::vector<float> >  ( prefix + "EtaError"                + suffix );
  produces <std::vector<float> >  ( prefix + "PhiError"                + suffix );
  produces <std::vector<float> >  ( prefix + "PtError"                 + suffix );
  produces <std::vector<float> >  ( prefix + "TrkEta"                  + suffix );
  produces <std::vector<float> >  ( prefix + "TrkPhi"                  + suffix );
  produces <std::vector<float> >  ( prefix + "TrkPt"                   + suffix );
  produces <std::vector<float> >  ( prefix + "TrkEtaError"             + suffix );
  produces <std::vector<float> >  ( prefix + "TrkPhiError"             + suffix );
  produces <std::vector<float> >  ( prefix + "TrkPtError"              + suffix );
  produces <std::vector<float> >  ( prefix + "QOverPError"             + suffix );
  produces <std::vector<float> >  ( prefix + "P"                       + suffix );
  produces <std::vector<float> >  ( prefix + "Energy"                  + suffix );
  produces <std::vector<int> >    ( prefix + "Charge"                  + suffix );
  produces <std::vector<int> >    ( prefix + "TrkHits"                 + suffix );
  produces <std::vector<int> >    ( prefix + "TrkHitsTrackerOnly"      + suffix );
  produces <std::vector<int> >    ( prefix + "GlobalTrkValidHits"      + suffix );
  produces <std::vector<int> >    ( prefix + "PixelHits"               + suffix );
  produces <std::vector<int> >    ( prefix + "TrkPixelHits"            + suffix );
  produces <std::vector<int> >    ( prefix + "SegmentMatches"          + suffix );
  produces <std::vector<int> >    ( prefix + "StationMatches"          + suffix );
  produces <std::vector<float> >  ( prefix + "TrkValidFractionOfHits"  + suffix );
  produces <std::vector<float> >  ( prefix + "TrkD0"                   + suffix );
  produces <std::vector<float> >  ( prefix + "TrkD0Error"              + suffix );
  produces <std::vector<float> >  ( prefix + "TrkDz"                   + suffix );
  produces <std::vector<float> >  ( prefix + "TrkDzError"              + suffix );
  produces <std::vector<float> >  ( prefix + "TrkVx"                   + suffix );
  produces <std::vector<float> >  ( prefix + "TrkVy"                   + suffix );
  produces <std::vector<float> >  ( prefix + "TrkVz"                   + suffix );
  produces <std::vector<float> >  ( prefix + "TrackChi2"               + suffix );
  produces <std::vector<float> >  ( prefix + "GlobalChi2"              + suffix );
  produces <std::vector<float> >  ( prefix + "CombinedQualityChi2LocalPosition" + suffix );
  produces <std::vector<float> >  ( prefix + "CombinedQualityTrkKink"  + suffix );
  produces <std::vector<float> >  ( prefix + "TrkIso"                  + suffix );
  produces <std::vector<float> >  ( prefix + "TrackerIsoSumPT"         + suffix );
  produces <std::vector<float> >  ( prefix + "EcalIso"                 + suffix );
  produces <std::vector<float> >  ( prefix + "HcalIso"                 + suffix );
  produces <std::vector<float> >  ( prefix + "HOIso"                   + suffix );
  produces <std::vector<float> >  ( prefix + "EcalVetoIso"             + suffix );
  produces <std::vector<float> >  ( prefix + "HcalVetoIso"             + suffix );
  produces <std::vector<float> >  ( prefix + "SegmentCompatibility"    + suffix );
  produces <std::vector<float> >  ( prefix + "PFIsoR03ChargedHadron"   + suffix );
  produces <std::vector<float> >  ( prefix + "PFIsoR03ChargedParticle" + suffix );
  produces <std::vector<float> >  ( prefix + "PFIsoR03NeutralHadron"   + suffix );
  produces <std::vector<float> >  ( prefix + "PFIsoR03Photon"          + suffix );
  produces <std::vector<float> >  ( prefix + "PFIsoR03NeutralHadronHT" + suffix );
  produces <std::vector<float> >  ( prefix + "PFIsoR03PhotonHT"        + suffix );
  produces <std::vector<float> >  ( prefix + "PFIsoR03PU"              + suffix );
  produces <std::vector<float> >  ( prefix + "PFIsoR04ChargedHadron"   + suffix );
  produces <std::vector<float> >  ( prefix + "PFIsoR04ChargedParticle" + suffix );
  produces <std::vector<float> >  ( prefix + "PFIsoR04NeutralHadron"   + suffix );
  produces <std::vector<float> >  ( prefix + "PFIsoR04Photon"          + suffix );
  produces <std::vector<float> >  ( prefix + "PFIsoR04NeutralHadronHT" + suffix );
  produces <std::vector<float> >  ( prefix + "PFIsoR04PhotonHT"        + suffix );
  produces <std::vector<float> >  ( prefix + "PFIsoR04PU"              + suffix );
  produces <std::vector<int> >    ( prefix + "PassID"                  + suffix );
  produces <std::vector<int> >    ( prefix + "VtxIndex"                + suffix );
  produces <std::vector<float> >  ( prefix + "VtxDistXY"               + suffix );
  produces <std::vector<float> >  ( prefix + "VtxDistZ"                + suffix );
  produces <std::vector<int> >    ( prefix + "BestTrackVtxIndex"      + suffix );
  produces <std::vector<float> >  ( prefix + "BestTrackVtxDistXY"     + suffix );
  produces <std::vector<float> >  ( prefix + "BestTrackVtxDistZ"      + suffix );
  produces <std::vector<float> >  ( prefix + "PrimaryVertexDXY"        + suffix );
  produces <std::vector<float> >  ( prefix + "PrimaryVertexDXYError"   + suffix );
  produces <std::vector<float> >  ( prefix + "BeamSpotDXY"             + suffix );
  produces <std::vector<float> >  ( prefix + "BeamSpotDXYError"        + suffix );
  produces <std::vector<int> >    ( prefix + "IsGlobal"                + suffix );
  produces <std::vector<int> >    ( prefix + "IsTracker"               + suffix );
  produces <std::vector<float> >  ( prefix + "MatchedGenParticlePt"    + suffix );
  produces <std::vector<float> >  ( prefix + "MatchedGenParticleEta"   + suffix );
  produces <std::vector<float> >  ( prefix + "MatchedGenParticlePhi"   + suffix );
  //
  // New variables added based on CMSSW 52X recommendations for LooseMuon and TightMuon Definitions
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Basline_muon_selections_for_2012
  produces <std::vector<int> >     ( prefix + "IsPF"                       + suffix );
  produces <std::vector<int> >     ( prefix + "TrackLayersWithMeasurement" + suffix );

  // Variables for trigger matching
  //  HLT Single Muon
  produces <std::vector<bool  > > ( prefix + "HLTSingleMuonMatched"      + suffix );
  produces <std::vector<float> >  ( prefix + "HLTSingleMuonMatchPt"      + suffix );
  produces <std::vector<float> >  ( prefix + "HLTSingleMuonMatchEta"     + suffix );
  produces <std::vector<float> >  ( prefix + "HLTSingleMuonMatchPhi"     + suffix );
  //  HLT Single Iso Muon
  produces <std::vector<bool  > > ( prefix + "HLTSingleIsoMuonMatched"   + suffix );
  produces <std::vector<float> >  ( prefix + "HLTSingleIsoMuonMatchPt"   + suffix );
  produces <std::vector<float> >  ( prefix + "HLTSingleIsoMuonMatchEta"  + suffix );
  produces <std::vector<float> >  ( prefix + "HLTSingleIsoMuonMatchPhi"  + suffix );

  //
  if ( useCocktailRefits )
    {
      produces <std::vector<int>    > ( prefix + "CocktailRefitID"                + suffix ) ;
      produces <std::vector<float> >  ( prefix + "CocktailEta"                    + suffix ) ;
      produces <std::vector<float> >  ( prefix + "CocktailPhi"                    + suffix ) ;
      produces <std::vector<float> >  ( prefix + "CocktailPt"                     + suffix ) ;
      produces <std::vector<float> >  ( prefix + "CocktailEtaError"               + suffix ) ;
      produces <std::vector<float> >  ( prefix + "CocktailPhiError"               + suffix ) ;
      produces <std::vector<float> >  ( prefix + "CocktailPtError"                + suffix ) ;
      produces <std::vector<float> >  ( prefix + "CocktailQOverPError"            + suffix ) ;
      produces <std::vector<float> >  ( prefix + "CocktailP"                      + suffix ) ;
      produces <std::vector<int   > > ( prefix + "CocktailCharge"                 + suffix ) ;
      produces <std::vector<int   > > ( prefix + "CocktailTrkHits"                + suffix ) ;
      produces <std::vector<float> >  ( prefix + "CocktailTrkValidFractionOfHits" + suffix ) ;
      produces <std::vector<float> >  ( prefix + "CocktailTrkD0"                  + suffix ) ;
      produces <std::vector<float> >  ( prefix + "CocktailTrkD0Error"             + suffix ) ;
      produces <std::vector<float> >  ( prefix + "CocktailTrkDz"                  + suffix ) ;
      produces <std::vector<float> >  ( prefix + "CocktailTrkDzError"             + suffix ) ;
      produces <std::vector<float> >  ( prefix + "CocktailGlobalChi2"             + suffix ) ;

      produces <std::vector<float> >  ( prefix + "CocktailTrkVtxDXY"                 + suffix ) ;
      produces <std::vector<float> >  ( prefix + "CocktailTrkVtxDZ"                  + suffix ) ;
      produces <std::vector<int> >    ( prefix + "CocktailTrkVtxIndex"               + suffix ) ;

    }
  /*
  produces <std::vector<float> > ( prefix + "CosmicCompatibility"     + suffix );
  produces <std::vector<float> > ( prefix + "TimeCompatibility"       + suffix );
  produces <std::vector<float> > ( prefix + "BackToBackCompatibility" + suffix );
  produces <std::vector<float> > ( prefix + "OverlapCompatibility"    + suffix );
  */
}


void RootTupleMakerV2_Muons::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::unique_ptr<bool>                  hasVeryForwardPFMuon    ( new bool() );
  std::unique_ptr<std::vector<bool> >    isLooseMuon             ( new std::vector<bool>()    );
  std::unique_ptr<std::vector<bool> >    isMediumMuon            ( new std::vector<bool>()    );
  std::unique_ptr<std::vector<bool> >    isTightMuon             ( new std::vector<bool>()    );
  std::unique_ptr<std::vector<bool> >    isHighPtMuon            ( new std::vector<bool>()    );
  std::unique_ptr<std::vector<float> >   eta                     ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   phi                     ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   pt                      ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   etaError                ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   phiError                ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   ptError                 ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   trkEta                  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   trkPhi                  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   trkPt                   ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   trkEtaError             ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   trkPhiError             ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   trkPtError              ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   qoverpError             ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   p                       ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   energy                  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<int> >     charge                  ( new std::vector<int>()     );
  std::unique_ptr<std::vector<int> >     trkHits                 ( new std::vector<int>()     );
  std::unique_ptr<std::vector<int> >     trkHitsTrackerOnly      ( new std::vector<int>()     );
  std::unique_ptr<std::vector<int> >     GlobaltrkValidHits      ( new std::vector<int>()     );
  std::unique_ptr<std::vector<int> >     pixelHits               ( new std::vector<int>()     );
  std::unique_ptr<std::vector<int> >     trkPixelHits            ( new std::vector<int>()     );
  std::unique_ptr<std::vector<int> >     segmentMatches          ( new std::vector<int>()     );
  std::unique_ptr<std::vector<int> >     stationMatches          ( new std::vector<int>()     );
  std::unique_ptr<std::vector<int> >     Valid                   ( new std::vector<int>()     );
  std::unique_ptr<std::vector<float> >   trkValidFractionOfHits  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   trkD0                   ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   trkD0Error              ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   trkDz                   ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   trkDzError              ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   trkVx                   ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   trkVy                   ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   trkVz                   ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   trackChi2               ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   globalChi2              ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   combinedQualityChi2LocalPosition ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   combinedQualityTrkKink  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   trkIso                  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   trackerIsoSumPT         ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   ecalIso                 ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   hcalIso                 ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   hoIso                   ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   ecalVetoIso             ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   hcalVetoIso             ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   segmentCompatibility    ( new std::vector<float>()  );
  //
  std::unique_ptr<std::vector<float> >   pfisor03chargedhadron   ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   pfisor03chargedparticle ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   pfisor03neutralhadron   ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   pfisor03photon          ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   pfisor03neutralhadronht ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   pfisor03photonht        ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   pfisor03pu              ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   pfisor04chargedhadron   ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   pfisor04chargedparticle ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   pfisor04neutralhadron   ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   pfisor04photon          ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   pfisor04neutralhadronht ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   pfisor04photonht        ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   pfisor04pu              ( new std::vector<float>()  );
  std::unique_ptr<std::vector<int> >     passID                  ( new std::vector<int>()     );
  std::unique_ptr<std::vector<int> >     vtxIndex                ( new std::vector<int>()     );
  std::unique_ptr<std::vector<float> >   vtxDistXY               ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   vtxDistZ                ( new std::vector<float>()  );
  std::unique_ptr<std::vector<int> >     bestTrackVtxIndex       ( new std::vector<int>()     );
  std::unique_ptr<std::vector<float> >   bestTrackVtxDistXY      ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   bestTrackVtxDistZ       ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   primaryVertexDXY        ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   primaryVertexDXYError   ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   beamspotDXY             ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   beamspotDXYError        ( new std::vector<float>()  );
  std::unique_ptr<std::vector<int> >     IsGlobal                ( new std::vector<int>()     );
  std::unique_ptr<std::vector<int> >     IsTracker               ( new std::vector<int>()     );
  std::unique_ptr<std::vector<float> >   matchedgenparticlept    ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   matchedgenparticleeta   ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   matchedgenparticlephi   ( new std::vector<float>()  );
  //
  // New variables added based on CMSSW 52X recommendations for LooseMuon and TightMuon Definitions
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Basline_muon_selections_for_2012
  std::unique_ptr<std::vector<int> >     isPF                       ( new std::vector<int>()    );
  std::unique_ptr<std::vector<int> >     trackLayersWithMeasurement ( new std::vector<int>()    );
  //
  // Trigger matching variables
  std::unique_ptr<std::vector<bool  > >  HLTSingleMuonMatched     ( new std::vector<bool  >()  );
  std::unique_ptr<std::vector<float> >   HLTSingleMuonMatchPt     ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   HLTSingleMuonMatchEta    ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   HLTSingleMuonMatchPhi    ( new std::vector<float>()  );
  std::unique_ptr<std::vector<bool  > >  HLTSingleIsoMuonMatched  ( new std::vector<bool  >()  );
  std::unique_ptr<std::vector<float> >   HLTSingleIsoMuonMatchPt  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   HLTSingleIsoMuonMatchEta ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   HLTSingleIsoMuonMatchPhi ( new std::vector<float>()  );
  //
  std::unique_ptr<std::vector<int   > >  ctRefitID    ( new std::vector<int   > () );
  std::unique_ptr<std::vector<float> >   ctEta        ( new std::vector<float> () );
  std::unique_ptr<std::vector<float> >   ctPhi        ( new std::vector<float> () );
  std::unique_ptr<std::vector<float> >   ctPt         ( new std::vector<float> () );
  std::unique_ptr<std::vector<float> >   ctEtaError   ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   ctPhiError   ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   ctPtError    ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   ctQoverpError( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >   ctP          ( new std::vector<float> () );
  std::unique_ptr<std::vector<int   > >  ctCharge     ( new std::vector<int   > () );
  std::unique_ptr<std::vector<int   > >  ctTrkHits    ( new std::vector<int   > () );
  std::unique_ptr<std::vector<float> >   ctTrkValidFractionOfHits ( new std::vector<float> () );
  std::unique_ptr<std::vector<float> >   ctTrkD0      ( new std::vector<float> () );
  std::unique_ptr<std::vector<float> >   ctTrkD0Error ( new std::vector<float> () );
  std::unique_ptr<std::vector<float> >   ctTrkDz      ( new std::vector<float> () );
  std::unique_ptr<std::vector<float> >   ctTrkDzError ( new std::vector<float> () );
  std::unique_ptr<std::vector<float> >   ctGlobalChi2 ( new std::vector<float> () );

  std::unique_ptr<std::vector<float> >   ctTrkvtxDistXY     ( new std::vector<float> () );
  std::unique_ptr<std::vector<float> >   ctTrkvtxDistZ      ( new std::vector<float> () );
  std::unique_ptr<std::vector<int> >     ctTrkvtxIndex      ( new std::vector<int> () );


  //
  /*
  std::unique_ptr<std::vector<float> >  cosmicCompatibility     ( new std::vector<float> () );
  std::unique_ptr<std::vector<float> >  timeCompatibility       ( new std::vector<float> () );
  std::unique_ptr<std::vector<float> >  backToBackCompatibility ( new std::vector<float> () );
  std::unique_ptr<std::vector<float> >  overlapCompatibility    ( new std::vector<float> () );
  */

  //-----------------------------------------------------------------
  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByToken(muonInputToken_, muons);

  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByToken(vtxInputToken_,primaryVertices);

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByToken(beamSpotToken_,beamSpot);

  //edm::Handle< pat::TriggerEvent > triggerEvent;
  //iEvent.getByLabel( triggerEventInputTag, triggerEvent );

  // PAT trigger matches by HLT path
  // we embed these in the regular pat::Muon collection 
  // but we could make separate collections containing only the matches if we wanted
  //edm::Handle<std::vector<pat::Muon> > muonsSingleMuonHLTMatched;
  //iEvent.getByLabel(inputTag, muonsSingleMuonHLTMatched);
  //edm::Handle<std::vector<pat::Muon> > muonsSingleIsoMuonHLTMatched;
  //iEvent.getByLabel(inputTag, muonsSingleIsoMuonHLTMatched);

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
	  //FIXME DMM: Update for pythia8. should use status 23 I think for outgoing muons.  Currently miniAOD stores only status 1
			
	  float genparPt = -999.;
	  float genparEta= -999.;
	  float genparPhi= -999.;
			
	  if ( !iEvent.isRealData() ){
	    //it->genParticleRefs().size() should be 0 or 1
	    for(uint igen = 0 ; igen < it->genParticleRefs().size() ; ++igen )
	      {
		if(it->genParticleRef(igen).isNonnull()) {
		  //std::cout<<"muon igen: "<<igen<<" pdgId: "<<it->genParticle(igen)->pdgId()<<" status: "<<it->genParticle(igen)->status()<<" pt: "<<it->genParticle(igen)->pt()<<std::endl;
		  if( it->genParticle(igen)->status()==1 || it->genParticle(igen)->status()==3 || it->genParticle(igen)->status()==23){
		    genparPt =it->genParticle(igen)->pt();
		    genparEta=it->genParticle(igen)->eta();
		    genparPhi=it->genParticle(igen)->phi();
		  }
		}
		else edm::LogError("RootTupleMakerV2_MuonsError") << "genParticleRef " << igen+1 << "/" << it->genParticleRefs().size() << " is null!";
	      }
	  }
	  matchedgenparticlept     -> push_back ( (float)(genparPt) );
	  matchedgenparticleeta    -> push_back ( (float)(genparEta) );
	  matchedgenparticlephi    -> push_back ( (float)(genparPhi) );
			
	  /// Gen Matching

	  float trkd0   = it->track()->d0();

	  if( beamSpotCorr && beamSpot.isValid() )
	    {
	      trkd0   = -(it->track()   ->dxy( beamSpot->position()));
	    }

	  else if( beamSpotCorr && !beamSpot.isValid() ) edm::LogError("RootTupleMakerV2_MuonsError") << "Error! Can't get the offlineBeamSpot";

	  // Vertex association
	  //-- using pat::muon::track()
	  float minVtxDist3D = 9999.;
	  int   vtxIndex_    = -1;
	  float vtxDistXY_   = -9999.;
	  float vtxDistZ_    = -9999.;
	  //-- using reco::muon::muonBestTrack() - 2012 recommendation
	  float bt_minVtxDist3D = 9999.;
	  int   bt_vtxIndex_    = -1;
	  float bt_vtxDistXY_   = -9999.;
	  float bt_vtxDistZ_    = -9999.;

	  if(primaryVertices.isValid())
	    {
	      edm::LogInfo("RootTupleMakerV2_MuonsInfo") << "Total # Primary Vertices: " << primaryVertices->size();

	      for( reco::VertexCollection::const_iterator v_it=primaryVertices->begin() ; v_it!=primaryVertices->end() ; ++v_it )
		{
		  //-- using pat::muon::track()
		  float distXY = it->track()->dxy(v_it->position());
		  float distZ  = it->track()->dz(v_it->position());
		  float dist3D = sqrt(pow(distXY,2) + pow(distZ,2));
		  if( dist3D < minVtxDist3D )
		    {
		      minVtxDist3D = dist3D;
		      vtxIndex_    = int(std::distance(primaryVertices->begin(),v_it));
		      vtxDistXY_   = distXY;
		      vtxDistZ_    = distZ;
		    }
		  //-- using reco::muon::muonBestTrack() - 2012 recommendation
		  if( (it->muonBestTrack()).isNonnull() )
		    {
		      float bt_distXY = it->muonBestTrack()->dxy(v_it->position());
		      float bt_distZ  = it->muonBestTrack()->dz(v_it->position());
		      float bt_dist3D = sqrt( pow(bt_distXY,2) + pow(bt_distZ,2) );
		      if( bt_dist3D < bt_minVtxDist3D )
			{
			  bt_minVtxDist3D = bt_dist3D;
			  bt_vtxIndex_    = int(std::distance(primaryVertices->begin(),v_it));
			  bt_vtxDistXY_   = bt_distXY;
			  bt_vtxDistZ_    = bt_distZ;
			}
		    }

		}				 //loop over primaryVertices
	    }
	  else
	    {
	      edm::LogError("RootTupleMakerV2_MuonsError") << "Error! Can't get the primary vertices";
	    }

	  //
	  //------------------------------------------------------------------------
	  // Trigger matching
	  //
	  // Example taken from PatTriggerAnalyzer:
	  // http://cmslxr.fnal.gov/lxr/source/PhysicsTools/PatExamples/plugins/PatTriggerAnalyzer.cc
	  //------------------------------------------------------------------------

	  //We've embedded matches to selected HLT paths in the python with the PATTriggerMatchEmbedder.
	  //now just ask if we have a match to whichever HLT path in the object
			
	  // HLT Single Muon trigger matching
	  const pat::TriggerObjectStandAloneCollection matchesSingleMu   = it->triggerObjectMatchesByPath("HLT_Mu50_v*");
	  const pat::TriggerObjectStandAloneCollection matchesSingleTkMu = it->triggerObjectMatchesByPath("HLT_TkMu50_v*");
	  if (matchesSingleMu.size() > 0 || matchesSingleTkMu.size() > 0)
	    {
	      if (matchesSingleMu.size() > 0){
		HLTSingleMuonMatched  -> push_back ( true ) ;
		HLTSingleMuonMatchPt  -> push_back ( matchesSingleMu[0].pt() );
		HLTSingleMuonMatchEta -> push_back ( matchesSingleMu[0].eta());
		HLTSingleMuonMatchPhi -> push_back ( matchesSingleMu[0].phi());
	      }
	      else if (matchesSingleTkMu.size() > 0){
		HLTSingleMuonMatched  -> push_back ( true ) ;
		HLTSingleMuonMatchPt  -> push_back ( matchesSingleTkMu[0].pt() );
		HLTSingleMuonMatchEta -> push_back ( matchesSingleTkMu[0].eta());
		HLTSingleMuonMatchPhi -> push_back ( matchesSingleTkMu[0].phi());
	      }
	    }
	  else
	    {
	      HLTSingleMuonMatched  -> push_back ( false ) ;
	      HLTSingleMuonMatchPt  -> push_back ( -999. );
	      HLTSingleMuonMatchEta -> push_back ( -999. );
	      HLTSingleMuonMatchPhi -> push_back ( -999. );
	    }

	  // HLT Single Iso Muon trigger matching
	  const pat::TriggerObjectStandAloneCollection matchesSingleIsoMu = it->triggerObjectMatchesByPath("HLT_IsoMu27_v*");
	  if (matchesSingleIsoMu.size() > 0)
	    {
	      HLTSingleIsoMuonMatched  -> push_back ( true ) ;
	      HLTSingleIsoMuonMatchPt  -> push_back ( matchesSingleIsoMu[0].pt() );
	      HLTSingleIsoMuonMatchEta -> push_back ( matchesSingleIsoMu[0].eta());
	      HLTSingleIsoMuonMatchPhi -> push_back ( matchesSingleIsoMu[0].phi());
	    }
	  else
	    {
	      HLTSingleIsoMuonMatched  -> push_back ( false ) ;
	      HLTSingleIsoMuonMatchPt  -> push_back ( -999. );
	      HLTSingleIsoMuonMatchEta -> push_back ( -999. );
	      HLTSingleIsoMuonMatchPhi -> push_back ( -999. );
	    }

	  reco::VertexCollection::const_iterator v_it_muId=primaryVertices->begin();
	  isLooseMuon->push_back(  it->isLooseMuon()  );
	  isMediumMuon->push_back( it->isMediumMuon() );
	  isTightMuon->push_back(  it->isTightMuon(*v_it_muId)  );
	  isHighPtMuon->push_back( it->isHighPtMuon(*v_it_muId) );

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
	      //track() returns innerTrack();
	      etaError    -> push_back ( it->track()->etaError()    );
	      //track() returns innerTrack();
	      phiError    -> push_back ( it->track()->phiError()    );
	      //track() returns innerTrack();
	      ptError     -> push_back ( it->track()->ptError ()    );
	      //track() returns innerTrack();
	      qoverpError -> push_back ( it->track()->qoverpError() );
	    }

	  trkPt  -> push_back ( it->track()->pt()  );
	  trkEta -> push_back ( it->track()->eta() );
	  trkPhi -> push_back ( it->track()->phi() );

	  trkPtError  -> push_back ( it->track()->ptError()  );
	  trkEtaError -> push_back ( it->track()->etaError() );
	  trkPhiError -> push_back ( it->track()->phiError() );

	  charge            ->push_back( it->charge() );
	  //	trkHits           ->push_back( it->track()->numberOfValidHits() );
	  trkHits           ->push_back( it->track()->hitPattern().numberOfValidHits() );
	  trkHitsTrackerOnly->push_back( it->track()->hitPattern().numberOfValidTrackerHits() );

	  if( it->isGlobalMuon() )
	    {
	      GlobaltrkValidHits->push_back( it->globalTrack()->hitPattern().numberOfValidMuonHits()  );
	      pixelHits         ->push_back( it->globalTrack()->hitPattern().numberOfValidPixelHits() );
	      globalChi2        ->push_back( it->globalTrack()->normalizedChi2()                      );
	    }
	  else
	    {
	      //inner track does not have muon hits by definition.
	      GlobaltrkValidHits->push_back( -1 );
	      //inner track pixel hit count is stored at "innerTrackPixelHits".
	      pixelHits         ->push_back( -1 );
	      //inner track Chi squared is already stored at  "trackChi2".
	      globalChi2        ->push_back( -1 );
	    }

	  combinedQualityChi2LocalPosition->push_back(it->combinedQuality().chi2LocalPosition);
	  combinedQualityTrkKink->push_back(it->combinedQuality().trkKink);

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
	  if ( useCocktailRefits )
	    {
	      if ( it->isGlobalMuon())	
		{
		  int refit_id = -999;
				
		  // -----------------------------------------------------------------------------------------------------------------------------------------//
		  // HighPT Muon - 2015 Version
		  // Following instructions on https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonId2015#HighPT_Muon
		  //When using the High-pT selector, you should obtain the muon momentum from the muon track determined by the TuneP algorithm. To have access to the best muon track determined by the TuneP algorithm and its type, you can use the following functions, respectively:
		  //reco::TrackRef tunePBestTrack = recoMu.tunePMuonBestTrack();
		  //reco::Muon::MuonTrackType tunePBestTrackType = recoMu.tunePMuonBestTrackType();
		  //The pT assignment as obtained from the TuneP algorithm can be used anyway, independently of the adopted selection. Anyhow it is recommended not to use it directly for analyses that make use of Particle-Flow AND rely on quantities computed on the basis of the PF muon pT (e.g. PF MET). In fact, when using Particle-Flow the muon pT assignment can be further adjusted on the basis of the global event description. In this case, it is therefore advisable to stick with the estimation of the muon pT made by PF.

		  // cktTrack->pt() - to get the pT of the muon
		  //
		  //Then you can apply the new HighPT ID, which still differs from Tight Muon selection in the following points:
		  //The Particle-Flow muon id is not required
		  //The cut of dpT/pT<0.3 for the track used for momentum determination is applied, i.e. cktTrack->ptError()/cktTrack->pt()<0.3
		  //The cuts applied on recoMu.muonBestTrack (impact parameter cuts) need to be applied to cktTrack since this is the new best track now.
		  // -----------------------------------------------------------------------------------------------------------------------------------------//
			
		  reco::TrackRef cocktail_track = it->tunePMuonBestTrack();

		  float cttrkd0  = cocktail_track -> d0() ;
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
		  ctGlobalChi2            ->push_back( cocktail_track ->normalizedChi2() );
				
				
		  // Values will hold vertex distance information
		  int   bct_vtxIndex_    = -1;
		  float bct_vtxDistXY_   = -9999.;
		  float bct_vtxDistZ_    = -9999.;
				
				
		  // Loop over primary vertices
		  if(primaryVertices.isValid())
		    {
		      float bct_bestdist3D = 999999.;
				    
		      for( reco::VertexCollection::const_iterator v_it=primaryVertices->begin() ; v_it!=primaryVertices->end() ; ++v_it )
			{
			  float bct_distXY = cocktail_track->dxy(v_it->position());
			  float bct_distZ  = cocktail_track->dz(v_it->position());
			  float bct_dist3D = sqrt( pow(bct_distXY,2) + pow(bct_distZ,2) );
					
			  if( bct_dist3D < bct_bestdist3D )
			    {
			      bct_bestdist3D = bct_dist3D;
			      bct_vtxIndex_    = int(std::distance(primaryVertices->begin(),v_it));
			      bct_vtxDistXY_   = bct_distXY;
			      bct_vtxDistZ_    = bct_distZ;
			    }
					
			}				 //loop over primaryVertices
		    }
				
				
		  ctTrkvtxDistXY           ->push_back( bct_vtxDistXY_ );
		  ctTrkvtxDistZ            ->push_back( bct_vtxDistZ_ );
		  ctTrkvtxIndex            ->push_back( bct_vtxIndex_ );
				
		  // std::cout<<" Cocktail flag --> "<< bct_vtxDistXY_ <<"  "<<bct_vtxDistZ_ << std::endl;
				
				
				
		} //cocktail fits
			    
	      else	
		{
		  ctRefitID     -> push_back ( -99 ) ;
		  ctEtaError    -> push_back ( -99 );
		  ctPhiError    -> push_back ( -99 );
		  ctPtError     -> push_back ( -99 );
		  ctQoverpError -> push_back ( -99 );
		  ctEta                    ->push_back( -99 );
		  ctPhi                    ->push_back( -99 );
		  ctPt                     ->push_back( -99 );
		  ctP                      ->push_back( -99 );
		  ctCharge                 ->push_back( -99 );
		  ctTrkHits                ->push_back( -99 );
		  ctTrkValidFractionOfHits ->push_back( -99 );
		  ctTrkD0                  ->push_back( -99 );
		  ctTrkD0Error             ->push_back( -99 );
		  ctTrkDz                  ->push_back( -99 );
		  ctTrkDzError             ->push_back( -99 );
		  ctGlobalChi2             ->push_back( -99 );
				
		  ctTrkvtxDistXY           ->push_back( -99 );
		  ctTrkvtxDistZ            ->push_back( -99 );					
		  ctTrkvtxIndex            ->push_back( -99 );
		}
			    
			    
	    }					 
			
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
	  segmentCompatibility -> push_back(it->segmentCompatibility());
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

	  // std::cout<<" Regular flag --> "<<  vtxDistXY_  <<"  "<< vtxDistZ_  << std::endl;


	  // See https://indico.cern.ch/getFile.py/access?contribId=7&resId=0&materialId=slides&confId=102306
	  // and https://indico.cern.ch/getFile.py/access?contribId=5&resId=0&materialId=slides&confId=128840
	  /*
	  cosmicCompatibility    ->push_back( it->userFloat("cosmicCompatibility")     );
	  timeCompatibility      ->push_back( it->userFloat("timeCompatibility")       );
	  backToBackCompatibility->push_back( it->userFloat("backToBackCompatibility") );
	  overlapCompatibility   ->push_back( it->userFloat("overlapCompatibility")    );
	  */
	  iMuon++;

	}						 //end of loop over muons
    }
  else
    {
      edm::LogError("RootTupleMakerV2_MuonsError") << "Error! Can't get the muons";
    }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put(std::move(hasVeryForwardPFMuon), "hasVeryForwardPFMuon" );
  iEvent.put(std::move(isLooseMuon),                prefix + "IsLooseMuon"                 + suffix );
  iEvent.put(std::move(isMediumMuon),               prefix + "IsMediumMuon"                + suffix );
  iEvent.put(std::move(isTightMuon),                prefix + "IsTightMuon"                 + suffix );
  iEvent.put(std::move(isHighPtMuon),               prefix + "IsHighPtMuon"                + suffix );
  iEvent.put(std::move(eta),                        prefix + "Eta"                         + suffix );
  iEvent.put(std::move(phi),                        prefix + "Phi"                         + suffix );
  iEvent.put(std::move(pt),                         prefix + "Pt"                          + suffix );
  iEvent.put(std::move(etaError),                   prefix + "EtaError"                    + suffix );
  iEvent.put(std::move(phiError),                   prefix + "PhiError"                    + suffix );
  iEvent.put(std::move(ptError),                    prefix + "PtError"                     + suffix );
  iEvent.put(std::move(trkEta),                     prefix + "TrkEta"                      + suffix );
  iEvent.put(std::move(trkPhi),                     prefix + "TrkPhi"                      + suffix );
  iEvent.put(std::move(trkPt),                      prefix + "TrkPt"                       + suffix );
  iEvent.put(std::move(trkEtaError),                prefix + "TrkEtaError"                 + suffix );
  iEvent.put(std::move(trkPhiError),                prefix + "TrkPhiError"                 + suffix );
  iEvent.put(std::move(trkPtError),                 prefix + "TrkPtError"                  + suffix );
  iEvent.put(std::move(qoverpError),                prefix + "QOverPError"                 + suffix );
  iEvent.put(std::move(p),                          prefix + "P"                           + suffix );
  iEvent.put(std::move(energy),                     prefix + "Energy"                      + suffix );
  iEvent.put(std::move(charge),                     prefix + "Charge"                      + suffix );
  iEvent.put(std::move(trkHits),                    prefix + "TrkHits"                     + suffix );
  iEvent.put(std::move(trkHitsTrackerOnly),         prefix + "TrkHitsTrackerOnly"          + suffix );
  iEvent.put(std::move(GlobaltrkValidHits),         prefix + "GlobalTrkValidHits"          + suffix );
  iEvent.put(std::move(pixelHits),                  prefix + "PixelHits"                   + suffix );
  iEvent.put(std::move(trkPixelHits),               prefix + "TrkPixelHits"                + suffix );
  iEvent.put(std::move(segmentMatches),             prefix + "SegmentMatches"              + suffix );
  iEvent.put(std::move(stationMatches),             prefix + "StationMatches"              + suffix );
  iEvent.put(std::move(trkValidFractionOfHits),     prefix + "TrkValidFractionOfHits"      + suffix );
  iEvent.put(std::move(trkD0),                      prefix + "TrkD0"                       + suffix );
  iEvent.put(std::move(trkD0Error),                 prefix + "TrkD0Error"                  + suffix );
  iEvent.put(std::move(trkDz),                      prefix + "TrkDz"                       + suffix );
  iEvent.put(std::move(trkDzError),                 prefix + "TrkDzError"                  + suffix );
  iEvent.put(std::move(trkVx),                      prefix + "TrkVx"                       + suffix );
  iEvent.put(std::move(trkVy),                      prefix + "TrkVy"                       + suffix );
  iEvent.put(std::move(trkVz),                      prefix + "TrkVz"                       + suffix );
  iEvent.put(std::move(trackChi2),                  prefix + "TrackChi2"                   + suffix );
  iEvent.put(std::move(globalChi2),                 prefix + "GlobalChi2"                  + suffix );
  iEvent.put(std::move(combinedQualityChi2LocalPosition), prefix + "CombinedQualityChi2LocalPosition" + suffix );
  iEvent.put(std::move(combinedQualityTrkKink),     prefix + "CombinedQualityTrkKink"      + suffix );
  iEvent.put(std::move(trkIso),                     prefix + "TrkIso"                      + suffix );
  iEvent.put(std::move(trackerIsoSumPT),            prefix + "TrackerIsoSumPT"             + suffix );
  iEvent.put(std::move(ecalIso),                    prefix + "EcalIso"                     + suffix );
  iEvent.put(std::move(hcalIso),                    prefix + "HcalIso"                     + suffix );
  iEvent.put(std::move(hoIso),                      prefix + "HOIso"                       + suffix );
  iEvent.put(std::move(ecalVetoIso),                prefix + "EcalVetoIso"                 + suffix );
  iEvent.put(std::move(hcalVetoIso),                prefix + "HcalVetoIso"                 + suffix );
  iEvent.put(std::move(segmentCompatibility),       prefix + "SegmentCompatibility"        + suffix );
  iEvent.put(std::move(pfisor03chargedhadron),      prefix + "PFIsoR03ChargedHadron"       + suffix );
  iEvent.put(std::move(pfisor03chargedparticle),    prefix + "PFIsoR03ChargedParticle"     + suffix );
  iEvent.put(std::move(pfisor03neutralhadron),      prefix + "PFIsoR03NeutralHadron"       + suffix );
  iEvent.put(std::move(pfisor03photon),             prefix + "PFIsoR03Photon"              + suffix );
  iEvent.put(std::move(pfisor03neutralhadronht),    prefix + "PFIsoR03NeutralHadronHT"     + suffix );
  iEvent.put(std::move(pfisor03photonht),           prefix + "PFIsoR03PhotonHT"            + suffix );
  iEvent.put(std::move(pfisor03pu),                 prefix + "PFIsoR03PU"                  + suffix );
  iEvent.put(std::move(pfisor04chargedhadron),      prefix + "PFIsoR04ChargedHadron"       + suffix );
  iEvent.put(std::move(pfisor04chargedparticle),    prefix + "PFIsoR04ChargedParticle"     + suffix );
  iEvent.put(std::move(pfisor04neutralhadron),      prefix + "PFIsoR04NeutralHadron"       + suffix );
  iEvent.put(std::move(pfisor04photon),             prefix + "PFIsoR04Photon"              + suffix );
  iEvent.put(std::move(pfisor04neutralhadronht),    prefix + "PFIsoR04NeutralHadronHT"     + suffix );
  iEvent.put(std::move(pfisor04photonht),           prefix + "PFIsoR04PhotonHT"            + suffix );
  iEvent.put(std::move(pfisor04pu),                 prefix + "PFIsoR04PU"                  + suffix );
  iEvent.put(std::move(passID),                     prefix + "PassID"                      + suffix );
  iEvent.put(std::move(IsGlobal),                   prefix + "IsGlobal"                    + suffix );
  iEvent.put(std::move(IsTracker),                  prefix + "IsTracker"                   + suffix );
  iEvent.put(std::move(vtxIndex),                   prefix + "VtxIndex"                    + suffix );
  iEvent.put(std::move(vtxDistXY),                  prefix + "VtxDistXY"                   + suffix );
  iEvent.put(std::move(vtxDistZ),                   prefix + "VtxDistZ"                    + suffix );
  iEvent.put(std::move(bestTrackVtxIndex),          prefix + "BestTrackVtxIndex"           + suffix );
  iEvent.put(std::move(bestTrackVtxDistXY),         prefix + "BestTrackVtxDistXY"          + suffix );
  iEvent.put(std::move(bestTrackVtxDistZ),          prefix + "BestTrackVtxDistZ"           + suffix );
  iEvent.put(std::move(primaryVertexDXY),           prefix + "PrimaryVertexDXY"            + suffix );
  iEvent.put(std::move(primaryVertexDXYError),      prefix + "PrimaryVertexDXYError"       + suffix );
  iEvent.put(std::move(beamspotDXY),                prefix + "BeamSpotDXY"                 + suffix );
  iEvent.put(std::move(beamspotDXYError),           prefix + "BeamSpotDXYError"            + suffix );
  iEvent.put(std::move(matchedgenparticlept),       prefix + "MatchedGenParticlePt"        + suffix );
  iEvent.put(std::move(matchedgenparticleeta),      prefix + "MatchedGenParticleEta"       + suffix );
  iEvent.put(std::move(matchedgenparticlephi),      prefix + "MatchedGenParticlePhi"       + suffix );
  iEvent.put(std::move(isPF),                       prefix + "IsPF"                        + suffix );
  iEvent.put(std::move(trackLayersWithMeasurement), prefix + "TrackLayersWithMeasurement"  + suffix );
  //

  if ( useCocktailRefits )
    {
      iEvent.put(std::move(ctRefitID                ), prefix + "CocktailRefitID"                 + suffix ) ;
      iEvent.put(std::move(ctEta                    ), prefix + "CocktailEta"                     + suffix ) ;
      iEvent.put(std::move(ctPhi                    ), prefix + "CocktailPhi"                     + suffix ) ;
      iEvent.put(std::move(ctPt                     ), prefix + "CocktailPt"                      + suffix ) ;
      iEvent.put(std::move(ctEtaError               ), prefix + "CocktailEtaError"                + suffix ) ;
      iEvent.put(std::move(ctPhiError               ), prefix + "CocktailPhiError"                + suffix ) ;
      iEvent.put(std::move(ctPtError                ), prefix + "CocktailPtError"                 + suffix ) ;
      iEvent.put(std::move(ctQoverpError            ), prefix + "CocktailQOverPError"             + suffix ) ;
      iEvent.put(std::move(ctP                      ), prefix + "CocktailP"                       + suffix ) ;
      iEvent.put(std::move(ctCharge                 ), prefix + "CocktailCharge"                  + suffix ) ;
      iEvent.put(std::move(ctTrkHits                ), prefix + "CocktailTrkHits"                 + suffix ) ;
      iEvent.put(std::move(ctTrkValidFractionOfHits ), prefix + "CocktailTrkValidFractionOfHits"  + suffix ) ;
      iEvent.put(std::move(ctTrkD0                  ), prefix + "CocktailTrkD0"                   + suffix ) ;
      iEvent.put(std::move(ctTrkD0Error             ), prefix + "CocktailTrkD0Error"              + suffix ) ;
      iEvent.put(std::move(ctTrkDz                  ), prefix + "CocktailTrkDz"                   + suffix ) ;
      iEvent.put(std::move(ctTrkDzError             ), prefix + "CocktailTrkDzError"              + suffix ) ;
      iEvent.put(std::move(ctGlobalChi2             ), prefix + "CocktailGlobalChi2"              + suffix ) ;

      iEvent.put(std::move(ctTrkvtxDistXY           ), prefix + "CocktailTrkVtxDXY"               + suffix ) ;
      iEvent.put(std::move(ctTrkvtxDistZ            ), prefix + "CocktailTrkVtxDZ"                + suffix ) ;
      iEvent.put(std::move(ctTrkvtxIndex            ), prefix + "CocktailTrkVtxIndex"                + suffix ) ;


    }

  // trigger matching
  iEvent.put(std::move(HLTSingleMuonMatched     ), prefix + "HLTSingleMuonMatched"      + suffix );
  iEvent.put(std::move(HLTSingleMuonMatchPt     ), prefix + "HLTSingleMuonMatchPt"      + suffix );
  iEvent.put(std::move(HLTSingleMuonMatchEta    ), prefix + "HLTSingleMuonMatchEta"     + suffix );
  iEvent.put(std::move(HLTSingleMuonMatchPhi    ), prefix + "HLTSingleMuonMatchPhi"     + suffix );
  iEvent.put(std::move(HLTSingleIsoMuonMatched  ), prefix + "HLTSingleIsoMuonMatched"   + suffix );
  iEvent.put(std::move(HLTSingleIsoMuonMatchPt  ), prefix + "HLTSingleIsoMuonMatchPt"   + suffix );
  iEvent.put(std::move(HLTSingleIsoMuonMatchEta ), prefix + "HLTSingleIsoMuonMatchEta"  + suffix );
  iEvent.put(std::move(HLTSingleIsoMuonMatchPhi ), prefix + "HLTSingleIsoMuonMatchPhi"  + suffix );
  /*
  iEvent.put(std::move(cosmicCompatibility),     prefix + "CosmicCompatibility"     + suffix );
  iEvent.put(std::move(timeCompatibility),       prefix + "TimeCompatibility"       + suffix );
  iEvent.put(std::move(backToBackCompatibility), prefix + "BackToBackCompatibility" + suffix );
  iEvent.put(std::move(overlapCompatibility),    prefix + "OverlapCompatibility"    + suffix );
  */
}
