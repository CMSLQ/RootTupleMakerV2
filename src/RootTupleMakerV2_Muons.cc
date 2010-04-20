#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Muons.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"


RootTupleMakerV2_Muons::RootTupleMakerV2_Muons(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize")),
    muonIso (iConfig.getParameter<double>       ("MuonIso")),
    muonID  (iConfig.getParameter<std::string>  ("MuonID")),
    beamSpotCorr (iConfig.getParameter<bool>    ("BeamSpotCorr"))
{
  produces <std::vector<double> > ( prefix + "Eta" + suffix );
  produces <std::vector<double> > ( prefix + "Phi" + suffix );
  produces <std::vector<double> > ( prefix + "Pt" + suffix );
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
}

void RootTupleMakerV2_Muons::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
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

  //-----------------------------------------------------------------
  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByLabel(inputTag, muons);

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByLabel("offlineBeamSpot", beamSpot );

  if(muons.isValid()) {
    edm::LogInfo("RootTupleMakerV2_MuonsInfo") << "total # Muons: " << muons->size();

    for( std::vector<pat::Muon>::const_iterator it = muons->begin(); it != muons->end();++it ) {
      // exit from loop when you reach the required number of muons
      if(eta->size() > maxSize)
        break;

      // if muon is not global muon, continue
      if(!it->isGlobalMuon())
        continue;

      double trkd0 = it->track()->d0();
      if( beamSpotCorr && beamSpot.isValid() ) trkd0 = -(it->track()->dxy( beamSpot->position()));
      else if( beamSpotCorr && !beamSpot.isValid() ) edm::LogError("RootTupleMakerV2_MuonsError") << "Error! Can't get the offlineBeamSpot";
      double reliso = (it->trackIso() + it->ecalIso() + it->hcalIso())/it->pt();

      eta->push_back( it->eta() );
      phi->push_back( it->phi() );
      pt->push_back( it->pt() );
      energy->push_back( it->energy() );
      charge->push_back( it->charge() );
      trkHits->push_back( it->track()->numberOfValidHits() );
      trkD0->push_back( trkd0 );
      trkD0Error->push_back( it->track()->d0Error() );
      trkDz->push_back( it->track()->dz() );
      trkDzError->push_back( it->track()->dzError() );
      globalChi2->push_back( it->track()->normalizedChi2() );
      trkIso->push_back( it->trackIso() );
      ecalIso->push_back( it->ecalIso() );
      hcalIso->push_back( it->hcalIso() );
      hoIso->push_back( it->isolationR03().hoEt );
      relIso->push_back( reliso );
      passIso->push_back( ((reliso<muonIso) ? 1 : 0) );
      passID->push_back( ((it->muonID(muonID)) ? 1 : 0) );
    }
  } else {
    edm::LogError("RootTupleMakerV2_MuonsError") << "Error! Can't get the product " << inputTag;
  }

  // put vectors in the event
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( pt, prefix + "Pt" + suffix );
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
}
