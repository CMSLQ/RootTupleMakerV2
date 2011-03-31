#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Muons.h"
#include "Leptoquarks/RootTupleMakerV2/interface/PatUtilities.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

RootTupleMakerV2_Muons::RootTupleMakerV2_Muons(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize")),
    muonIso (iConfig.getParameter<double>       ("MuonIso")),
    muonID  (iConfig.getParameter<std::string>  ("MuonID")),
    beamSpotCorr (iConfig.getParameter<bool>    ("BeamSpotCorr")),
    useCocktailRefits ( iConfig.getParameter<bool>("UseCocktailRefits")),
    vtxInputTag(iConfig.getParameter<edm::InputTag>("VertexInputTag"))
{
  produces <std::vector<double> > ( prefix + "Eta" + suffix );
  produces <std::vector<double> > ( prefix + "Phi" + suffix );
  produces <std::vector<double> > ( prefix + "Pt" + suffix );
  produces <std::vector<double> > ( prefix + "EtaError" + suffix );
  produces <std::vector<double> > ( prefix + "PhiError" + suffix );
  produces <std::vector<double> > ( prefix + "PtError" + suffix );
  produces <std::vector<double> > ( prefix + "QOverPError" + suffix );
  produces <std::vector<double> > ( prefix + "P" + suffix );
  produces <std::vector<double> > ( prefix + "Energy" + suffix );
  produces <std::vector<int> >    ( prefix + "Charge" + suffix );
  produces <std::vector<int> >    ( prefix + "TrkHits" + suffix );
  produces <std::vector<double> > ( prefix + "TrkD0" + suffix );
  produces <std::vector<double> > ( prefix + "TrkD0Error" + suffix );
  produces <std::vector<double> > ( prefix + "TrkDz" + suffix );
  produces <std::vector<double> > ( prefix + "TrkDzError" + suffix );
  produces <std::vector<double> > ( prefix + "GlobalChi2" + suffix );
  produces <std::vector<double> > ( prefix + "TrkIso" + suffix );
  produces <std::vector<double> > ( prefix + "EcalIso" + suffix );
  produces <std::vector<double> > ( prefix + "HcalIso" + suffix );
  produces <std::vector<double> > ( prefix + "HOIso" + suffix );
  produces <std::vector<double> > ( prefix + "RelIso" + suffix );
  produces <std::vector<int> >    ( prefix + "PassIso" + suffix );
  produces <std::vector<int> >    ( prefix + "PassID" + suffix );
  produces <std::vector<int> >    ( prefix + "VtxIndex" + suffix );
  produces <std::vector<double> > ( prefix + "VtxDistXY" + suffix );
  produces <std::vector<double> > ( prefix + "VtxDistZ" + suffix );

  if ( useCocktailRefits ) {
    produces <std::vector<double> > ( prefix + "CocktailEta"        + suffix ) ;
    produces <std::vector<double> > ( prefix + "CocktailPhi"        + suffix ) ;
    produces <std::vector<double> > ( prefix + "CocktailPt"         + suffix ) ;
    produces <std::vector<double> > ( prefix + "CocktailEtaError"   + suffix ) ;
    produces <std::vector<double> > ( prefix + "CocktailPhiError"   + suffix ) ;
    produces <std::vector<double> > ( prefix + "CocktailPtError"    + suffix ) ;
    produces <std::vector<double> > ( prefix + "CocktailQOverPError"+ suffix ) ;
    produces <std::vector<double> > ( prefix + "CocktailP"          + suffix ) ;
    produces <std::vector<int   > > ( prefix + "CocktailCharge"     + suffix ) ;
    produces <std::vector<int   > > ( prefix + "CocktailTrkHits"    + suffix ) ;
    produces <std::vector<double> > ( prefix + "CocktailTrkD0"      + suffix ) ;
    produces <std::vector<double> > ( prefix + "CocktailTrkD0Error" + suffix ) ;
    produces <std::vector<double> > ( prefix + "CocktailTrkDz"      + suffix ) ;
    produces <std::vector<double> > ( prefix + "CocktailTrkDzError" + suffix ) ;
    produces <std::vector<double> > ( prefix + "CocktailGlobalChi2" + suffix ) ;
    produces <std::vector<double> > ( prefix + "CocktailRelIso"     + suffix ) ;
    produces <std::vector<int   > > ( prefix + "CocktailPassIso"    + suffix ) ;
  }

}

void RootTupleMakerV2_Muons::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  etaError  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phiError  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ptError  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  qoverpError  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  p  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     charge  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >     trkHits ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  trkD0   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkD0Error ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkDz   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkDzError ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  globalChi2 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkIso   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ecalIso  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIso  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hoIso    ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  relIso   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     passIso  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >     passID   ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >     vtxIndex  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  vtxDistXY  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  vtxDistZ  ( new std::vector<double>()  );

  std::auto_ptr<std::vector<double> >  ctEta        ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  ctPhi        ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  ctPt         ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  ctEtaError  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ctPhiError  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ctPtError  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ctQoverpError  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ctP          ( new std::vector<double> () );
  std::auto_ptr<std::vector<int   > >  ctCharge     ( new std::vector<int   > () );
  std::auto_ptr<std::vector<int   > >  ctTrkHits    ( new std::vector<int   > () );
  std::auto_ptr<std::vector<double> >  ctTrkD0      ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  ctTrkD0Error ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  ctTrkDz      ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  ctTrkDzError ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  ctGlobalChi2 ( new std::vector<double> () );
  std::auto_ptr<std::vector<double> >  ctRelIso     ( new std::vector<double> () );
  std::auto_ptr<std::vector<int   > >  passCTIso    ( new std::vector<int   > () );
  
  //-----------------------------------------------------------------
  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByLabel(inputTag, muons);

  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(vtxInputTag,primaryVertices);

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByLabel("offlineBeamSpot", beamSpot );

  if(muons.isValid()) {
    edm::LogInfo("RootTupleMakerV2_MuonsInfo") << "Total # Muons: " << muons->size();

    for( std::vector<pat::Muon>::const_iterator it = muons->begin(); it != muons->end(); ++it ) {
      // exit from loop when you reach the required number of muons
      if(eta->size() >= maxSize)
        break;

      // if muon is not global muon, continue
      if(!it->isGlobalMuon())
        continue;

      double trkd0   = it->track()->d0();

      if( beamSpotCorr && beamSpot.isValid() ) {
	trkd0   = -(it->track()   ->dxy( beamSpot->position()));
      }
      
      else if( beamSpotCorr && !beamSpot.isValid() ) edm::LogError("RootTupleMakerV2_MuonsError") << "Error! Can't get the offlineBeamSpot";

      double reliso   = (it->trackIso() + it->ecalIso() + it->hcalIso())/it->pt();

      // Vertex association
      double minVtxDist3D = 9999.;
      int vtxIndex_ = -1;
      double vtxDistXY_ = -9999.;
      double vtxDistZ_ = -9999.;

      if(primaryVertices.isValid()) {
        edm::LogInfo("RootTupleMakerV2_MuonsInfo") << "Total # Primary Vertices: " << primaryVertices->size();

        for( reco::VertexCollection::const_iterator v_it=primaryVertices->begin() ; v_it!=primaryVertices->end() ; ++v_it ) {

          double distXY = it->track()->dxy(v_it->position());
          double distZ = it->track()->dz(v_it->position());
          double dist3D = sqrt(pow(distXY,2) + pow(distZ,2));

          if( dist3D<minVtxDist3D ) {
            minVtxDist3D = dist3D;
            vtxIndex_ = int(std::distance(primaryVertices->begin(),v_it));
            vtxDistXY_ = distXY;
            vtxDistZ_ = distZ;
          }
        }
      } else {
        edm::LogError("RootTupleMakerV2_MuonsError") << "Error! Can't get the product " << vtxInputTag;
      }

      eta->push_back( it->eta() );
      phi->push_back( it->phi() );
      pt->push_back( it->pt() );
      p->push_back( it->p() );

      etaError    -> push_back ( it->globalTrack()->etaError() );
      phiError    -> push_back ( it->globalTrack()->phiError() );
      ptError     -> push_back ( it->globalTrack()->ptError () );
      qoverpError -> push_back ( it->globalTrack()->qoverpError() ) ;

      charge->push_back( it->charge() );
      trkHits->push_back( it->track()->numberOfValidHits() );
      trkD0->push_back( trkd0 );
      trkD0Error->push_back( it->track()->d0Error() );
      trkDz->push_back( it->track()->dz() );
      trkDzError->push_back( it->track()->dzError() );
      globalChi2->push_back( it->track()->normalizedChi2() );
      relIso->push_back( reliso );
      passIso->push_back( (reliso<muonIso) ? 1 : 0 );

      if ( useCocktailRefits ) {

	const reco::TrackRef& cocktail_track = pmcTrack(*it);
	
	double ctreliso = (it->trackIso() + it->ecalIso() + it->hcalIso())/cocktail_track->pt();
	double cttrkd0  = cocktail_track -> d0() ;
	if( beamSpotCorr && beamSpot.isValid() ) 
	  cttrkd0 = -(cocktail_track->dxy( beamSpot->position()));
	
	ctEtaError    -> push_back ( cocktail_track->etaError() );
	ctPhiError    -> push_back ( cocktail_track->phiError() );
	ctPtError     -> push_back ( cocktail_track->ptError () );
	ctQoverpError -> push_back ( cocktail_track->qoverpError() ) ;

	ctEta       ->push_back( cocktail_track->eta() );
	ctPhi       ->push_back( cocktail_track->phi() );
	ctPt        ->push_back( cocktail_track->pt() );
	ctP         ->push_back( cocktail_track->p() );
	ctCharge    ->push_back( cocktail_track->charge() );
	ctTrkHits   ->push_back( cocktail_track->numberOfValidHits() );
	ctTrkD0     ->push_back( cttrkd0 ) ;
	ctTrkD0Error->push_back( cocktail_track->d0Error() );
	ctTrkDz     ->push_back( cocktail_track->dz() );
	ctTrkDzError->push_back( cocktail_track -> dzError() );
	ctGlobalChi2->push_back( cocktail_track ->normalizedChi2() );
	ctRelIso    ->push_back( ctreliso ) ;
	passCTIso   ->push_back( (ctreliso<muonIso) ? 1 : 0 );
      }
      
      energy->push_back( it->energy() );
      trkIso->push_back( it->trackIso() );
      ecalIso->push_back( it->ecalIso() );
      hcalIso->push_back( it->hcalIso() );
      hoIso->push_back( it->isolationR03().hoEt );

      passID->push_back( (it->muonID(muonID)) ? 1 : 0 );
      vtxIndex->push_back( vtxIndex_ );
      vtxDistXY->push_back( vtxDistXY_ );
      vtxDistZ->push_back( vtxDistZ_ );
    }
  } else {
    edm::LogError("RootTupleMakerV2_MuonsError") << "Error! Can't get the product " << inputTag;
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( pt, prefix + "Pt" + suffix );
  iEvent.put( etaError, prefix + "EtaError" + suffix );
  iEvent.put( phiError, prefix + "PhiError" + suffix );
  iEvent.put( ptError, prefix + "PtError" + suffix );
  iEvent.put( qoverpError, prefix + "QOverPError" + suffix );
  iEvent.put( p, prefix + "P" + suffix );
  iEvent.put( energy, prefix + "Energy" + suffix );
  iEvent.put( charge, prefix + "Charge" + suffix );
  iEvent.put( trkHits, prefix + "TrkHits" + suffix );
  iEvent.put( trkD0, prefix + "TrkD0" + suffix );
  iEvent.put( trkD0Error, prefix + "TrkD0Error" + suffix );
  iEvent.put( trkDz, prefix + "TrkDz" + suffix );
  iEvent.put( trkDzError, prefix + "TrkDzError" + suffix );
  iEvent.put( globalChi2, prefix + "GlobalChi2" + suffix ); 
  iEvent.put( trkIso, prefix + "TrkIso" + suffix );
  iEvent.put( ecalIso, prefix + "EcalIso" + suffix );
  iEvent.put( hcalIso, prefix + "HcalIso" + suffix );
  iEvent.put( hoIso, prefix + "HOIso" + suffix );
  iEvent.put( relIso, prefix + "RelIso" + suffix );
  iEvent.put( passIso, prefix + "PassIso" + suffix );
  iEvent.put( passID, prefix + "PassID" + suffix );
  iEvent.put( vtxIndex, prefix + "VtxIndex" + suffix );
  iEvent.put( vtxDistXY, prefix + "VtxDistXY" + suffix );
  iEvent.put( vtxDistZ, prefix + "VtxDistZ" + suffix );

  if ( useCocktailRefits ) {
    iEvent.put( ctEta        , prefix + "CocktailEta"         + suffix ) ;
    iEvent.put( ctPhi        , prefix + "CocktailPhi"         + suffix ) ;
    iEvent.put( ctPt         , prefix + "CocktailPt"          + suffix ) ;
    iEvent.put( ctEtaError   , prefix + "CocktailEtaError"    + suffix ) ;
    iEvent.put( ctPhiError   , prefix + "CocktailPhiError"    + suffix ) ;
    iEvent.put( ctPtError    , prefix + "CocktailPtError"     + suffix ) ;
    iEvent.put( ctQoverpError, prefix + "CocktailQOverPError" + suffix );
    iEvent.put( ctP          , prefix + "CocktailP"           + suffix ) ;
    iEvent.put( ctCharge     , prefix + "CocktailCharge"      + suffix ) ;
    iEvent.put( ctTrkHits    , prefix + "CocktailTrkHits"     + suffix ) ;
    iEvent.put( ctTrkD0      , prefix + "CocktailTrkD0"       + suffix ) ;
    iEvent.put( ctTrkD0Error , prefix + "CocktailTrkD0Error"  + suffix ) ;
    iEvent.put( ctTrkDz      , prefix + "CocktailTrkDz"       + suffix ) ;
    iEvent.put( ctTrkDzError , prefix + "CocktailTrkDzError"  + suffix ) ;
    iEvent.put( ctGlobalChi2 , prefix + "CocktailGlobalChi2"  + suffix ) ;
    iEvent.put( ctRelIso     , prefix + "CocktailRelIso"      + suffix ) ;
    iEvent.put( passCTIso    , prefix + "CocktailPassIso"     + suffix ) ;
  }
}
