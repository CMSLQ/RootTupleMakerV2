#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Electrons.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"


RootTupleMakerV2_Electrons::RootTupleMakerV2_Electrons(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize")),
    electronIso (iConfig.getParameter<double>   ("ElectronIso")),
    muonPt (iConfig.getParameter<double>        ("MuonPt")),
    muonIso (iConfig.getParameter<double>       ("MuonIso")),
    muonID  (iConfig.getParameter<std::string>  ("MuonID"))
{
  produces <std::vector<double> > ( prefix + "Eta" + suffix );
  produces <std::vector<double> > ( prefix + "Phi" + suffix );
  produces <std::vector<double> > ( prefix + "Pt" + suffix );
  //produces <std::vector<double> > ( prefix + "EtHeep" + suffix );
  produces <std::vector<double> > ( prefix + "TrackPt" + suffix );
  produces <std::vector<double> > ( prefix + "Energy" + suffix );
  produces <std::vector<double> > ( prefix + "CaloEnergy" + suffix );
  produces <std::vector<int> >    ( prefix + "Charge" + suffix );
  produces <std::vector<int> >    ( prefix + "Overlaps" + suffix );
  produces <std::vector<double> > ( prefix + "HoE" + suffix );
  produces <std::vector<double> > ( prefix + "SigmaEtaEta" + suffix );
  produces <std::vector<double> > ( prefix + "SigmaIEtaIEta" + suffix );
  produces <std::vector<double> > ( prefix + "DeltaPhiTrkSC" + suffix );
  produces <std::vector<double> > ( prefix + "DeltaEtaTrkSC" + suffix );
  produces <std::vector<int> >    ( prefix + "Classif" + suffix );
  produces <std::vector<double> > ( prefix + "E1x5OverE5x5" + suffix );
  produces <std::vector<double> > ( prefix + "E2x5OverE5x5" + suffix );
  produces <std::vector<int> >    ( prefix + "HeepID" + suffix );
  produces <std::vector<int> >    ( prefix + "PassID" + suffix );
  produces <std::vector<double> > ( prefix + "TrkIso" + suffix );
  produces <std::vector<double> > ( prefix + "EcalIso" + suffix );
  produces <std::vector<double> > ( prefix + "HcalIso" + suffix );
  produces <std::vector<double> > ( prefix + "RelIso" + suffix );
  produces <std::vector<int> >    ( prefix + "PassIso" + suffix );
  produces <std::vector<double> > ( prefix + "EcalIsoHeep" + suffix );
  produces <std::vector<double> > ( prefix + "HcalIsoD1Heep" + suffix );
  produces <std::vector<double> > ( prefix + "HcalIsoD2Heep" + suffix );
  produces <std::vector<double> > ( prefix + "TrkIsoHeep" + suffix );
  produces <std::vector<double> > ( prefix + "SCEta" + suffix );
  produces <std::vector<double> > ( prefix + "SCPhi" + suffix );
  produces <std::vector<double> > ( prefix + "SCPt" + suffix );
  produces <std::vector<double> > ( prefix + "SCRawEnergy" + suffix );
}

void RootTupleMakerV2_Electrons::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
  //std::auto_ptr<std::vector<double> >  etHeep  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trackPt  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  caloEnergy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     charge  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >     overlaps  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  hoe   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  sigmaEtaEta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  sigmaIEtaIEta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  deltaPhiTrkSC  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  deltaEtaTrkSC  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     classif  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  e1x5overe5x5  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  e2x5overe5x5  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     heepID  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >     passID  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  trkIso   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ecalIso  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIso  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  relIso   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     passIso  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  ecalIsoHeep  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIsoD1Heep  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIsoD2Heep  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkIsoHeep  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scEta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scPhi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scPt  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scRawEnergy  ( new std::vector<double>()  );

  //-----------------------------------------------------------------
  edm::Handle<std::vector<pat::Electron> > electrons;
  iEvent.getByLabel(inputTag, electrons);

  if(electrons.isValid()) {
    edm::LogInfo("RootTupleMakerV2_ElectronsInfo") << "Total # Electrons: " << electrons->size();

    for( std::vector<pat::Electron>::const_iterator it = electrons->begin(); it != electrons->end(); ++it ) {
      // exit from loop when you reach the required number of electrons
      if(eta->size() >= maxSize)
        break;

      // if electron is not ECAL driven, continue
      if(!it->ecalDrivenSeed())
        continue;

      int ovrlps = 0;
      const reco::CandidatePtrVector & muons = it->overlaps("muons");
      for (size_t i = 0; i < muons.size(); ++i) {
        // try to cast into pat::Muon
        const pat::Muon *muon = dynamic_cast<const pat::Muon *>(&*muons[i]);
        if (muon) {
           if( muon->muonID(muonID) 
               && ((muon->trackIso()+muon->ecalIso()+muon->hcalIso())/muon->pt())<muonIso
               && muon->pt()>muonPt ) ovrlps = 1;
        }
      }
      double reliso = (it->trackIso() + it->ecalIso() + it->hcalIso())/it->pt();
      int passId = 0;
      /* passID for different electron IDs is assigned bitwise
         bit 0: eidRobustLoose
         bit 1: eidRobustTight
         bit 2: eidLoose
         bit 3: eidTight
         bit 4: eidRobustHighEnergy
         bit 5: HEEPId
      */
      if (it->electronID("eidRobustLoose")>0) passId = passId | 1<<0;
      if (it->electronID("eidRobustTight")>0) passId = passId | 1<<1;
      if (it->electronID("eidLoose")>0) passId = passId | 1<<2;
      if (it->electronID("eidTight")>0) passId = passId | 1<<3;
      if (it->electronID("eidRobustHighEnergy")>0) passId = passId | 1<<4;
      if (it->userInt("HEEPId")==0) passId = passId | 1<<5;

      eta->push_back( it->eta() );
      phi->push_back( it->phi() );
      pt->push_back( it->pt() );
      //etHeep->push_back( it->caloEnergy()*sin( it->p4().theta() ) );
      trackPt->push_back( it->gsfTrack()->pt() );
      energy->push_back( it->energy() );
      caloEnergy->push_back( it->caloEnergy() );
      charge->push_back( it->charge() );
      overlaps->push_back( ovrlps );
      // ID variables
      hoe->push_back( it->hadronicOverEm() );
      sigmaEtaEta->push_back( it->scSigmaEtaEta() );
      sigmaIEtaIEta->push_back( it->scSigmaIEtaIEta() );
      deltaPhiTrkSC->push_back( it->deltaPhiSuperClusterTrackAtVtx() );
      deltaEtaTrkSC->push_back( it->deltaEtaSuperClusterTrackAtVtx() );
      classif->push_back( it->classification() );
      if( it->e5x5()>0 )
	{
	  e1x5overe5x5->push_back( double ( it->e1x5() / it->e5x5() ) );
	  e2x5overe5x5->push_back( double ( it->e2x5Max() / it->e5x5() ) );
	} 
      heepID->push_back( it->userInt("HEEPId") );
      passID->push_back( passId );
      // Iso variables
      trkIso->push_back( it->trackIso() );
      ecalIso->push_back( it->ecalIso() );
      hcalIso->push_back( it->hcalIso() );
      relIso->push_back( reliso );
      passIso->push_back( (reliso<electronIso) ? 1 : 0 );
      // Iso variables (Heep)
      ecalIsoHeep->push_back( it->dr03EcalRecHitSumEt() );
      hcalIsoD1Heep->push_back( it->dr03HcalDepth1TowerSumEt() );
      hcalIsoD2Heep->push_back( it->dr03HcalDepth2TowerSumEt() );
      trkIsoHeep->push_back( it->dr03TkSumPt() );
      // SC associated with electron
      scEta->push_back( it->superCluster()->eta() );
      scPhi->push_back( it->superCluster()->phi() );
      scPt->push_back( it->superCluster()->energy()/cosh(it->superCluster()->eta()) );
      scRawEnergy->push_back( it->superCluster()->rawEnergy() );
    }
  } else {
    edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get the product " << inputTag;
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( pt, prefix + "Pt" + suffix );
  //iEvent.put( etHeep, prefix + "EtHeep" + suffix );
  iEvent.put( trackPt, prefix + "TrackPt" + suffix );
  iEvent.put( energy, prefix + "Energy" + suffix );
  iEvent.put( caloEnergy, prefix + "CaloEnergy" + suffix );
  iEvent.put( charge, prefix + "Charge" + suffix );
  iEvent.put( overlaps, prefix + "Overlaps" + suffix );
  iEvent.put( hoe, prefix + "HoE" + suffix );
  iEvent.put( sigmaEtaEta, prefix + "SigmaEtaEta" + suffix );
  iEvent.put( sigmaIEtaIEta, prefix + "SigmaIEtaIEta" + suffix );
  iEvent.put( deltaPhiTrkSC, prefix + "DeltaPhiTrkSC" + suffix );
  iEvent.put( deltaEtaTrkSC, prefix + "DeltaEtaTrkSC" + suffix );
  iEvent.put( classif, prefix + "Classif" + suffix );
  iEvent.put( e1x5overe5x5, prefix + "E1x5OverE5x5" + suffix );
  iEvent.put( e2x5overe5x5, prefix + "E2x5OverE5x5" + suffix );
  iEvent.put( heepID, prefix + "HeepID" + suffix );
  iEvent.put( passID, prefix + "PassID" + suffix );
  iEvent.put( trkIso, prefix + "TrkIso" + suffix );
  iEvent.put( ecalIso, prefix + "EcalIso" + suffix );
  iEvent.put( hcalIso, prefix + "HcalIso" + suffix );
  iEvent.put( relIso, prefix + "RelIso" + suffix );
  iEvent.put( passIso, prefix + "PassIso" + suffix );
  iEvent.put( ecalIsoHeep, prefix + "EcalIsoHeep" + suffix );
  iEvent.put( hcalIsoD1Heep, prefix + "HcalIsoD1Heep" + suffix );
  iEvent.put( hcalIsoD2Heep, prefix + "HcalIsoD2Heep" + suffix );
  iEvent.put( trkIsoHeep, prefix + "TrkIsoHeep" + suffix );
  iEvent.put( scEta, prefix + "SCEta" + suffix );
  iEvent.put( scPhi, prefix + "SCPhi" + suffix );
  iEvent.put( scPt, prefix + "SCPt" + suffix );
  iEvent.put( scRawEnergy, prefix + "SCRawEnergy" + suffix );
}
