#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Electrons.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


RootTupleMakerV2_Electrons::RootTupleMakerV2_Electrons(const edm::ParameterSet& iConfig) :
    trkInputTag(iConfig.getParameter<edm::InputTag>("TracksInputTag")),
    dcsInputTag(iConfig.getParameter<edm::InputTag>("DCSInputTag")),
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize")),
    electronIso (iConfig.getParameter<double>   ("ElectronIso")),
    muonPt (iConfig.getParameter<double>        ("MuonPt")),
    muonIso (iConfig.getParameter<double>       ("MuonIso")),
    muonID  (iConfig.getParameter<std::string>  ("MuonID")),
    vtxInputTag(iConfig.getParameter<edm::InputTag>("VertexInputTag"))
{
  produces <std::vector<double> > ( prefix + "Eta" + suffix );
  produces <std::vector<double> > ( prefix + "Phi" + suffix );
  produces <std::vector<double> > ( prefix + "Pt" + suffix );
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
  produces <std::vector<int> >    ( prefix + "MissingHits" + suffix );
  produces <std::vector<double> > ( prefix + "Dist" + suffix );
  produces <std::vector<double> > ( prefix + "DCotTheta" + suffix );
  produces <std::vector<double> > ( prefix + "SCEta" + suffix );
  produces <std::vector<double> > ( prefix + "SCPhi" + suffix );
  produces <std::vector<double> > ( prefix + "SCPt" + suffix );
  produces <std::vector<double> > ( prefix + "SCRawEnergy" + suffix );
  produces <std::vector<double> > ( prefix + "VtxDist3D" + suffix );
  produces <std::vector<int> >    ( prefix + "VtxIndex" + suffix );
}

void RootTupleMakerV2_Electrons::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
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
  std::auto_ptr<std::vector<int> >     missingHits  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  dist_vec  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  dCotTheta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scEta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scPhi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scPt  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scRawEnergy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  vtxDist3D  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     vtxIndex  ( new std::vector<int>()  );

  //-----------------------------------------------------------------
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(trkInputTag, tracks);

  edm::Handle<DcsStatusCollection> dcsHandle;
  iEvent.getByLabel(dcsInputTag, dcsHandle);

  double evt_bField = 3.8;
  // need the magnetic field
  //
  // if isRealData then derive bfield using the
  // magnet current from DcsStatus
  // otherwise take it from the IdealMagneticFieldRecord
  if(iEvent.isRealData()) {
    if(dcsHandle.isValid()) {
      edm::LogInfo("RootTupleMakerV2_ElectronsInfo") << "Successfully obtained " << dcsInputTag;
      // scale factor = 3.801/18166.0 which are
      // average values taken over a stable two-week period
      double currentToBFieldScaleFactor = 2.09237036221512717e-04;
      double current = (*dcsHandle)[0].magnetCurrent();
      evt_bField = current*currentToBFieldScaleFactor;
    } else {
      edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get the product " << dcsInputTag;
    }
  } else {
    edm::ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

    if(magneticField.isValid()) {
      edm::LogInfo("RootTupleMakerV2_ElectronsInfo") << "Successfully obtained IdealMagneticFieldRecord";

      evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
    } else {
      edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get IdealMagneticFieldRecord";
    }
  }

  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(vtxInputTag,primaryVertices);

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

      // Overlaps
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

      // Conversion
      ConversionFinder convFinder;
      double dist = -9999.;
      double dcot = -9999.;
      if(tracks.isValid()) {
        edm::LogInfo("RootTupleMakerV2_ElectronsInfo") << "Successfully obtained " << trkInputTag;

        ConversionInfo convInfo = convFinder.getConversionInfo(*it, tracks, evt_bField);
        dist = convInfo.dist();
        dcot = convInfo.dcot();
      } else {
        edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get the product " << trkInputTag;
      }

      // Vertex association
      double minVtxDist3D = 9999.;
      int indexVtx = -1;

      if(primaryVertices.isValid()) {
        edm::LogInfo("RootTupleMakerV2_ElectronsInfo") << "Total # Primary Vertices: " << primaryVertices->size();

        for( reco::VertexCollection::const_iterator v_it=primaryVertices->begin() ; v_it!=primaryVertices->end() ; ++v_it ) {

          double dist3D = sqrt(pow(it->gsfTrack()->dxy(v_it->position()),2) + pow(it->gsfTrack()->dz(v_it->position()),2));

          if( dist3D<minVtxDist3D ) {
            minVtxDist3D = dist3D;
            indexVtx = int(std::distance(primaryVertices->begin(),v_it));
          }
        }
      } else {
        edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get the product " << vtxInputTag;
      }

      eta->push_back( it->eta() );
      phi->push_back( it->phi() );
      pt->push_back( it->pt() );
      trackPt->push_back( it->gsfTrack()->pt() );
      energy->push_back( it->energy() );
      caloEnergy->push_back( it->caloEnergy() );
      charge->push_back( it->charge() );
      overlaps->push_back( ovrlps );
      // ID variables
      hoe->push_back( it->hadronicOverEm() );
      sigmaEtaEta->push_back( it->sigmaEtaEta() );
      sigmaIEtaIEta->push_back( it->sigmaIetaIeta() );
      deltaPhiTrkSC->push_back( it->deltaPhiSuperClusterTrackAtVtx() );
      deltaEtaTrkSC->push_back( it->deltaEtaSuperClusterTrackAtVtx() );
      classif->push_back( it->classification() );
      e1x5overe5x5->push_back( (it->e5x5()>0) ? (it->e1x5()/it->e5x5()) : 0 );
      e2x5overe5x5->push_back( (it->e5x5()>0) ? (it->e2x5Max()/it->e5x5()) : 0 );
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
      // Conversion variables
      missingHits->push_back( it->gsfTrack()->trackerExpectedHitsInner().numberOfHits() );
      dist_vec->push_back( dist );
      dCotTheta->push_back( dcot );
      // SC associated with electron
      scEta->push_back( it->superCluster()->eta() );
      scPhi->push_back( it->superCluster()->phi() );
      scPt->push_back( it->superCluster()->energy()/cosh(it->superCluster()->eta()) );
      scRawEnergy->push_back( it->superCluster()->rawEnergy() );
      // Vertex association variables
      vtxDist3D->push_back( minVtxDist3D );
      vtxIndex->push_back( indexVtx );
    }
  } else {
    edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get the product " << inputTag;
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( pt, prefix + "Pt" + suffix );
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
  iEvent.put( missingHits, prefix + "MissingHits" + suffix );
  iEvent.put( dist_vec, prefix + "Dist" + suffix );
  iEvent.put( dCotTheta, prefix + "DCotTheta" + suffix );
  iEvent.put( scEta, prefix + "SCEta" + suffix );
  iEvent.put( scPhi, prefix + "SCPhi" + suffix );
  iEvent.put( scPt, prefix + "SCPt" + suffix );
  iEvent.put( scRawEnergy, prefix + "SCRawEnergy" + suffix );
  iEvent.put( vtxDist3D, prefix + "VtxDist3D" + suffix );
  iEvent.put( vtxIndex, prefix + "VtxIndex" + suffix );
}
