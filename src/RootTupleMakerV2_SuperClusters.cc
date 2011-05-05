#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_SuperClusters.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/PhotonTkIsolation.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaRecHitIsolation.h"
#include "RecoEgamma/EgammaTools/interface/HoECalculator.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "TVector3.h"


RootTupleMakerV2_SuperClusters::RootTupleMakerV2_SuperClusters(const edm::ParameterSet& iConfig) :
    ebInputTag(iConfig.getParameter<edm::InputTag>("EBInputTag")),
    eeInputTag(iConfig.getParameter<edm::InputTag>("EEInputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    ebMaxSize (iConfig.getParameter<unsigned int> ("EBMaxSize")),
    eeMaxSize (iConfig.getParameter<unsigned int> ("EEMaxSize")),
    ecalEBInputTag(iConfig.getParameter<edm::InputTag>("EcalEBInputTag")),
    ecalEEInputTag(iConfig.getParameter<edm::InputTag>("EcalEEInputTag")),
    trkInputTag(iConfig.getParameter<edm::InputTag>("TracksInputTag")),
    eleInputTag(iConfig.getParameter<edm::InputTag>("ElectronsInputTag")),
    elePrefix  (iConfig.getParameter<std::string>  ("ElectronsPrefix")),
    eleMaxSize (iConfig.getParameter<unsigned int> ("ElectronsMaxSize"))
{
  produces <std::vector<double> > ( prefix + "Eta" + suffix );
  produces <std::vector<double> > ( prefix + "Phi" + suffix );
  produces <std::vector<double> > ( prefix + "Pt" + suffix );
  produces <std::vector<double> > ( prefix + "RawEnergy" + suffix );
  produces <std::vector<int> >    ( prefix + "ClustersSize" + suffix );
  produces <std::vector<double> > ( prefix + "HoE" + suffix );
  produces <std::vector<double> > ( prefix + "SigmaIEtaIEta" + suffix );
  produces <std::vector<double> > ( prefix + "EcalIso" + suffix );
  produces <std::vector<double> > ( prefix + "E1E9" + suffix );
  produces <std::vector<double> > ( prefix + "S4S1" + suffix );
  produces <std::vector<int> >    ( prefix + "kOutOfTime" + suffix );
  produces <std::vector<double> > ( prefix + "HEEPEcalIso" + suffix );
  produces <std::vector<double> > ( prefix + "HEEPTrkIso" + suffix );
  produces <std::vector<int> >    ( prefix + "TrackMatch" + suffix );
  produces <std::vector<double> > ( prefix + "DrTrack1" + suffix );
  produces <std::vector<double> > ( prefix + "Track1Eta" + suffix );
  produces <std::vector<double> > ( prefix + "Track1Phi" + suffix );
  produces <std::vector<double> > ( prefix + "Track1Pt" + suffix );
  produces <std::vector<double> > ( prefix + "DrTrack2" + suffix );
  produces <std::vector<double> > ( prefix + "Track2Eta" + suffix );
  produces <std::vector<double> > ( prefix + "Track2Phi" + suffix );
  produces <std::vector<double> > ( prefix + "Track2Pt" + suffix );
  produces <std::vector<double> > ( elePrefix + "SCE1E9" + suffix );
  produces <std::vector<double> > ( elePrefix + "SCS4S1" + suffix );
  produces <std::vector<int> >    ( elePrefix + "SCkOutOfTime" + suffix );
  produces <std::vector<double> > ( elePrefix + "SCEcalIso" + suffix );
  produces <std::vector<double> > ( elePrefix + "SCHEEPEcalIso" + suffix );
  produces <std::vector<double> > ( elePrefix + "SCHEEPTrkIso" + suffix );
}

void RootTupleMakerV2_SuperClusters::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  rawEnergy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     clustersSize  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  hoe   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  sigmaIEtaIEta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ecalIso  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  e1e9  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  s4s1  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     koutoftime  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  heepEcalIso  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  heepTrkIso  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     trackMatch  ( new std::vector<int>()  );;
  std::auto_ptr<std::vector<double> >  dRTrack1   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  track1Eta   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  track1Phi   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  track1Pt   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  dRTrack2   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  track2Eta   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  track2Phi   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  track2Pt   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scE1E9  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scS4S1  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     sckOutOfTime  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  scEcalIso   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scHEEPEcalIso   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scHEEPTrkIso   ( new std::vector<double>()  );

  //-----------------------------------------------------------------
  edm::Handle<reco::SuperClusterCollection> superClustersEBHandle;
  iEvent.getByLabel(ebInputTag, superClustersEBHandle);

  edm::Handle<reco::SuperClusterCollection> superClustersEEHandle;
  iEvent.getByLabel(eeInputTag, superClustersEEHandle);

  edm::ESHandle<CaloGeometry> caloGeometry;
  iSetup.get<CaloGeometryRecord>().get(caloGeometry);

  edm::ESHandle<CaloTopology> caloTopology;
  iSetup.get<CaloTopologyRecord>().get(caloTopology);

   edm::ESHandle<EcalSeverityLevelAlgo> sevlv;
   iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);
   const EcalSeverityLevelAlgo* sevLevel = sevlv.product();

  edm::Handle<EcalRecHitCollection> ecalBarrelRecHitHandle;
  iEvent.getByLabel(ecalEBInputTag, ecalBarrelRecHitHandle);

  edm::Handle<EcalRecHitCollection> ecalEndcapRecHitHandle;
  iEvent.getByLabel(ecalEEInputTag, ecalEndcapRecHitHandle);

  EcalClusterLazyTools EcalTool(iEvent,iSetup,ecalEBInputTag,ecalEEInputTag);

  edm::Handle<reco::TrackCollection> trackHandle;
  iEvent.getByLabel(trkInputTag, trackHandle);
  const reco::TrackCollection* tracks = trackHandle.product();
  edm::Handle<reco::BeamSpot> BeamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot",BeamSpotHandle);
  const reco::BeamSpot* spot = BeamSpotHandle.product();
  PhotonTkIsolation TrackTool(0.3,0.04,0.7,0.2,9999,tracks,math::XYZPoint(spot->x0(),spot->y0(),spot->z0()));

  double ConeOutRadius = 0.4;  // these are all needed to make an instance of EgammaRecHitIsolation
  double ConeInRadius = 0.045;
  double EtaWidth = 0.02;
  double PtMin = 0.;
  double EMin = 0.08;
  double EndcapConeInRadius = 0.07; // note that this is in number-of-crystal units
  double EndcapPtMin = 0.;
  double EndcapEMin = 0.1;
  double HeepConeOutRadius = 0.3;  // HEEP uses different iso values than Egamma/PAT default
  double HeepConeInRadius = 3.; // note that this is in num of crystals
  double HeepEtaWidth = 1.5;  // note that this is in num of crystals
  EcalRecHitMetaCollection ecalBarrelHits(*ecalBarrelRecHitHandle);// these are all needed to make an instance of EgammaRecHitIsolation
  EgammaRecHitIsolation ecalBarrelIsol(ConeOutRadius,ConeInRadius,EtaWidth,PtMin,EMin,caloGeometry,&ecalBarrelHits,sevLevel,DetId::Ecal);
  EgammaRecHitIsolation HeepEcalBarrelIsol(HeepConeOutRadius,HeepConeInRadius,HeepEtaWidth,PtMin,EMin,caloGeometry,&ecalBarrelHits,sevLevel,DetId::Ecal);
  HeepEcalBarrelIsol.setUseNumCrystals(true);
  EcalRecHitMetaCollection ecalEndcapHits(*ecalEndcapRecHitHandle);
  EgammaRecHitIsolation ecalEndcapIsol(ConeOutRadius,EndcapConeInRadius,EtaWidth,EndcapPtMin,EndcapEMin,caloGeometry,&ecalEndcapHits,sevLevel,DetId::Ecal);
  EgammaRecHitIsolation HeepEcalEndcapIsol(HeepConeOutRadius,HeepConeInRadius,HeepEtaWidth,EndcapPtMin,EndcapEMin,caloGeometry,&ecalEndcapHits,sevLevel,DetId::Ecal);
  HeepEcalEndcapIsol.setUseNumCrystals(true);

  // SuperClusters for barrel and endcap are in different collections
  // "hybrid" = barrel, "multi5x5" = endcap
  // Loop over both collections

  for( reco::SuperClusterCollection::const_iterator it = superClustersEBHandle->begin(); it != superClustersEBHandle->end(); ++it ) {
    // exit from loop when you reach the required number of barrel superclusters
    if( eta->size() >= ebMaxSize )
      break;

    eta->push_back( it->eta() );
    phi->push_back( it->phi() );
    TVector3 sc_vec;
    sc_vec.SetXYZ(it->x(),it->y(),it->z());
    pt->push_back( it->energy()*(sc_vec.Perp()/sc_vec.Mag()) );
    rawEnergy->push_back( it->rawEnergy() );
    clustersSize->push_back( it->clustersSize() );

    const reco::SuperCluster* pnt_sc = &(*it);
    HoECalculator calc_HoE; // requires HCAL RecHits that are only available in RECO files
    double schoe = calc_HoE(pnt_sc,iEvent,iSetup);
    hoe->push_back( schoe );

    std::vector<float> scLocalCov = EcalTool.scLocalCovariances(*it);
    double scSigmaiEiE= sqrt(scLocalCov[0]); // same method used in GsfElectronAlgo.cc
    sigmaIEtaIEta->push_back( scSigmaiEiE );

    reco::SuperClusterRef tempSCRef(superClustersEBHandle,std::distance(superClustersEBHandle->begin(),it));  // get SCRef to use to make ele candidate
    reco::RecoEcalCandidate ecalCand;  // make ele candidate to use Iso algorithm
    ecalCand.setSuperCluster(tempSCRef);
    ecalIso->push_back( ecalBarrelIsol.getEtSum(&ecalCand) );
    heepTrkIso->push_back( TrackTool.getPtTracks(&ecalCand) );
    heepEcalIso->push_back( HeepEcalBarrelIsol.getEtSum(&ecalCand) );

    double closestTrkDr = 99.;
    reco::Track closestTrk;
    double nextTrkDr = 99.;
    reco::Track nextTrk;

    for ( reco::TrackCollection::const_iterator trk = trackHandle->begin(); trk != trackHandle->end(); ++trk ) {
      TVector3 trk_vec;
      trk_vec.SetPtEtaPhi(trk->pt(),trk->eta(),trk->phi());
      double dR = trk_vec.DeltaR(sc_vec);
      if (dR<closestTrkDr) {
        nextTrkDr = closestTrkDr;
        nextTrk = closestTrk;
        closestTrkDr = dR;
        closestTrk = *trk;
      } else if (dR<nextTrkDr) {
        nextTrkDr = dR;
        nextTrk = *trk;
      }
    }
    if (closestTrkDr<0.04) trackMatch->push_back( 1 );  // 1=true, 0=false
    else trackMatch->push_back( 0 );
    dRTrack1->push_back( closestTrkDr );
    dRTrack2->push_back( nextTrkDr );
    track1Eta->push_back( closestTrk.eta() );
    track2Eta->push_back( nextTrk.eta() );
    track1Phi->push_back( closestTrk.phi() );
    track2Phi->push_back( nextTrk.phi() );
    track1Pt->push_back( closestTrk.pt() );
    track2Pt->push_back( nextTrk.pt() );

    double emax    = EcalClusterTools::eMax( *it, &(*ecalBarrelRecHitHandle) );
    double e9      = EcalClusterTools::e3x3( *it, &(*ecalBarrelRecHitHandle), &(*caloTopology) );
    double eright  = EcalClusterTools::eRight( *it, &(*ecalBarrelRecHitHandle), &(*caloTopology) );
    double eleft   = EcalClusterTools::eLeft( *it, &(*ecalBarrelRecHitHandle), &(*caloTopology) ) ;
    double etop    = EcalClusterTools::eTop( *it, &(*ecalBarrelRecHitHandle), &(*caloTopology) ) ;
    double ebottom = EcalClusterTools::eBottom( *it, &(*ecalBarrelRecHitHandle), &(*caloTopology) );

    e1e9->push_back( emax/e9 );
    s4s1->push_back( (eright+eleft+etop+ebottom)/emax );

    int flaggedRecHitCounter = 0;
    const std::vector<std::pair<DetId, float> > & hitsAndFractions = it->seed()->hitsAndFractions();
    std::vector<std::pair<DetId, float> >::const_iterator hitIter;

    for( hitIter=hitsAndFractions.begin(); hitIter!=hitsAndFractions.end(); ++hitIter ) {
      EcalRecHitCollection::const_iterator recHit = ecalBarrelRecHitHandle->find(hitIter->first);
      if(recHit!=ecalEndcapRecHitHandle->end()) {
        if( (hitIter->second*recHit->energy())/it->rawEnergy()>0.05 && recHit->checkFlag(EcalRecHit::kOutOfTime) ) flaggedRecHitCounter++;
      }
    }

    koutoftime->push_back( (flaggedRecHitCounter>0) ? 1 : 0 );
  }

  for( reco::SuperClusterCollection::const_iterator it = superClustersEEHandle->begin(); it != superClustersEEHandle->end(); ++it ) {
    // exit from loop when you reach the required number of endcap superclusters
    if( eta->size() >= (ebMaxSize + eeMaxSize) )
      break;

    eta->push_back( it->eta() );
    phi->push_back( it->phi() );
    TVector3 sc_vec;
    sc_vec.SetXYZ(it->x(),it->y(),it->z());
    pt->push_back( it->energy()*(sc_vec.Perp()/sc_vec.Mag()) );
    rawEnergy->push_back( it->rawEnergy() );
    clustersSize->push_back( it->clustersSize() );

    const reco::SuperCluster* pnt_sc = &(*it);
    HoECalculator calc_HoE; // requires HCAL RecHits that are only available in RECO files
    double schoe = calc_HoE(pnt_sc,iEvent,iSetup);
    hoe->push_back( schoe );

    std::vector<float> scLocalCov = EcalTool.scLocalCovariances(*it);
    double scSigmaiEiE= sqrt(scLocalCov[0]); // same method used in GsfElectronAlgo.cc
    sigmaIEtaIEta->push_back( scSigmaiEiE );

    reco::SuperClusterRef tempSCRef(superClustersEEHandle,std::distance(superClustersEEHandle->begin(),it));  // get SCRef to use to make ele candidate
    reco::RecoEcalCandidate ecalCand;  // make ele candidate to use Iso algorithm
    ecalCand.setSuperCluster(tempSCRef);
    ecalIso->push_back( ecalEndcapIsol.getEtSum(&ecalCand) );
    heepEcalIso->push_back( HeepEcalEndcapIsol.getEtSum(&ecalCand) );
    heepTrkIso->push_back( TrackTool.getPtTracks(&ecalCand) );

    double closestTrkDr = 99.;
    reco::Track closestTrk;
    double nextTrkDr = 99.;
    reco::Track nextTrk;

    for ( reco::TrackCollection::const_iterator trk = trackHandle->begin(); trk != trackHandle->end(); ++trk ) {
      TVector3 trk_vec;
      trk_vec.SetPtEtaPhi(trk->pt(),trk->eta(),trk->phi());
      double dR = trk_vec.DeltaR(sc_vec);
      if (dR<closestTrkDr) {
        nextTrkDr = closestTrkDr;
        nextTrk = closestTrk;
        closestTrkDr = dR;
        closestTrk = *trk;
      } else if (dR<nextTrkDr) {
        nextTrkDr = dR;
        nextTrk = *trk;
      }
    }
    if (closestTrkDr<0.04) trackMatch->push_back( 1 );  // 1=true, 0=false
    else trackMatch->push_back( 0 );
    dRTrack1->push_back( closestTrkDr );
    dRTrack2->push_back( nextTrkDr );
    track1Eta->push_back( closestTrk.eta() );
    track2Eta->push_back( nextTrk.eta() );
    track1Phi->push_back( closestTrk.phi() );
    track2Phi->push_back( nextTrk.phi() );
    track1Pt->push_back( closestTrk.pt() );
    track2Pt->push_back( nextTrk.pt() );

    double emax    = EcalClusterTools::eMax( *it, &(*ecalEndcapRecHitHandle) );
    double e9      = EcalClusterTools::e3x3( *it, &(*ecalEndcapRecHitHandle), &(*caloTopology) );
    double eright  = EcalClusterTools::eRight( *it, &(*ecalEndcapRecHitHandle), &(*caloTopology) );
    double eleft   = EcalClusterTools::eLeft( *it, &(*ecalEndcapRecHitHandle), &(*caloTopology) ) ;
    double etop    = EcalClusterTools::eTop( *it, &(*ecalEndcapRecHitHandle), &(*caloTopology) ) ;
    double ebottom = EcalClusterTools::eBottom( *it, &(*ecalEndcapRecHitHandle), &(*caloTopology) );

    e1e9->push_back( emax/e9 );
    s4s1->push_back( (eright+eleft+etop+ebottom)/emax );

    int flaggedRecHitCounter = 0;
    const std::vector<std::pair<DetId, float> > & hitsAndFractions = it->seed()->hitsAndFractions();
    std::vector<std::pair<DetId, float> >::const_iterator hitIter;

    for( hitIter=hitsAndFractions.begin(); hitIter!=hitsAndFractions.end(); ++hitIter ) {
      EcalRecHitCollection::const_iterator recHit = ecalEndcapRecHitHandle->find(hitIter->first);
      if(recHit!=ecalEndcapRecHitHandle->end()) {
        if( (hitIter->second*recHit->energy())/it->rawEnergy()>0.05 && recHit->checkFlag(EcalRecHit::kOutOfTime) ) flaggedRecHitCounter++;
      }
    }

    koutoftime->push_back( (flaggedRecHitCounter>0) ? 1 : 0 );
  }

  // Electrons
  edm::Handle<std::vector<pat::Electron> > electrons;
  iEvent.getByLabel(eleInputTag, electrons);

  if(electrons.isValid()) {
    edm::LogInfo("RootTupleMakerV2_SuperClustersInfo") << "Total # Electrons: " << electrons->size();

    for( std::vector<pat::Electron>::const_iterator it = electrons->begin(); it != electrons->end(); ++it ) {
      // exit from loop when you reach the required number of electrons
      if(scE1E9->size() >= eleMaxSize)
        break;

      // if electron is not ECAL driven, continue
      if(!it->ecalDrivenSeed())
        continue;

      // Ecal Spike Cleaning
      double emax    = -1.;
      double e9      = 1/999.;
      double eright  = 999.;
      double eleft   = 0.;
      double etop    = 0.;
      double ebottom = 0.;
      int ekoutoftime = -1;
      if( it->superCluster()->seed()->seed().subdetId() == EcalBarrel ) {
        emax    = EcalClusterTools::eMax( *(it->superCluster()), &(*ecalBarrelRecHitHandle) );
        e9      = EcalClusterTools::e3x3( *(it->superCluster()), &(*ecalBarrelRecHitHandle), &(*caloTopology) );
        eright  = EcalClusterTools::eRight( *(it->superCluster()), &(*ecalBarrelRecHitHandle), &(*caloTopology) );
        eleft   = EcalClusterTools::eLeft( *(it->superCluster()), &(*ecalBarrelRecHitHandle), &(*caloTopology) ) ;
        etop    = EcalClusterTools::eTop( *(it->superCluster()), &(*ecalBarrelRecHitHandle), &(*caloTopology) ) ;
        ebottom = EcalClusterTools::eBottom( *(it->superCluster()), &(*ecalBarrelRecHitHandle), &(*caloTopology) );

        int flaggedRecHitCounter = 0;
        const std::vector<std::pair<DetId, float> > & hitsAndFractions = it->superCluster()->seed()->hitsAndFractions();
        std::vector<std::pair<DetId, float> >::const_iterator hitIter;

        for( hitIter=hitsAndFractions.begin(); hitIter!=hitsAndFractions.end(); ++hitIter ) {
          EcalRecHitCollection::const_iterator recHit = ecalBarrelRecHitHandle->find(hitIter->first);
          if(recHit!=ecalEndcapRecHitHandle->end()) {
            if( (hitIter->second*recHit->energy())/it->superCluster()->rawEnergy()>0.05 && recHit->checkFlag(EcalRecHit::kOutOfTime) ) flaggedRecHitCounter++;
          }
        }

        ekoutoftime = (flaggedRecHitCounter>0) ? 1 : 0;
      } else {
        emax    = EcalClusterTools::eMax( *(it->superCluster()), &(*ecalEndcapRecHitHandle) );
        e9      = EcalClusterTools::e3x3( *(it->superCluster()), &(*ecalEndcapRecHitHandle), &(*caloTopology) );
        eright  = EcalClusterTools::eRight( *(it->superCluster()), &(*ecalEndcapRecHitHandle), &(*caloTopology) );
        eleft   = EcalClusterTools::eLeft( *(it->superCluster()), &(*ecalEndcapRecHitHandle), &(*caloTopology) ) ;
        etop    = EcalClusterTools::eTop( *(it->superCluster()), &(*ecalEndcapRecHitHandle), &(*caloTopology) ) ;
        ebottom = EcalClusterTools::eBottom( *(it->superCluster()), &(*ecalEndcapRecHitHandle), &(*caloTopology) );

        int flaggedRecHitCounter = 0;
        const std::vector<std::pair<DetId, float> > & hitsAndFractions = it->superCluster()->seed()->hitsAndFractions();
        std::vector<std::pair<DetId, float> >::const_iterator hitIter;

        for( hitIter=hitsAndFractions.begin(); hitIter!=hitsAndFractions.end(); ++hitIter ) {
          EcalRecHitCollection::const_iterator recHit = ecalEndcapRecHitHandle->find(hitIter->first);
          if(recHit!=ecalEndcapRecHitHandle->end()) {
            if( (hitIter->second*recHit->energy())/it->superCluster()->rawEnergy()>0.05 && recHit->checkFlag(EcalRecHit::kOutOfTime) ) flaggedRecHitCounter++;
          }
        }

        ekoutoftime = (flaggedRecHitCounter>0) ? 1 : 0;
      }

      scE1E9->push_back( emax/e9 );
      scS4S1->push_back( (eright+eleft+etop+ebottom)/emax );
      sckOutOfTime->push_back( ekoutoftime );

      reco::SuperClusterRef eleSCRef = it->superCluster();  // get SCRef to use to make ele candidate
      TVector3 sc_vec;
      sc_vec.SetXYZ(eleSCRef->x(),eleSCRef->y(),eleSCRef->z());
      double eleSCRefpt = eleSCRef->energy()*(sc_vec.Perp()/sc_vec.Mag());
      reco::RecoEcalCandidate ecalCand;  // make ele candidate to use Iso algorithm
      ecalCand.setSuperCluster(eleSCRef);
      const reco::Candidate::PolarLorentzVector photon_vec(eleSCRefpt, eleSCRef->eta(), eleSCRef->phi(), 0.0);
      ecalCand.setP4( photon_vec );
      if (fabs(eleSCRef->eta())<1.48) scEcalIso->push_back( ecalBarrelIsol.getEtSum(&ecalCand) );
      else scEcalIso->push_back( ecalEndcapIsol.getEtSum(&ecalCand) );
      if (fabs(eleSCRef->eta())<1.48) scHEEPEcalIso->push_back( HeepEcalBarrelIsol.getEtSum(&ecalCand) );
      else scHEEPEcalIso->push_back( HeepEcalEndcapIsol.getEtSum(&ecalCand) );
      scHEEPTrkIso->push_back( TrackTool.getPtTracks(&ecalCand) );
    }
  } else {
    edm::LogError("RootTupleMakerV2_SuperClustersError") << "Error! Can't get the product " << eleInputTag;
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( pt, prefix + "Pt" + suffix );
  iEvent.put( rawEnergy, prefix + "RawEnergy" + suffix );
  iEvent.put( clustersSize, prefix + "ClustersSize" + suffix );
  iEvent.put( hoe, prefix + "HoE" + suffix );
  iEvent.put( sigmaIEtaIEta, prefix + "SigmaIEtaIEta" + suffix );
  iEvent.put( ecalIso, prefix + "EcalIso" + suffix );
  iEvent.put( e1e9, prefix + "E1E9" + suffix );
  iEvent.put( s4s1, prefix + "S4S1" + suffix );
  iEvent.put( koutoftime, prefix + "kOutOfTime" + suffix );
  iEvent.put( heepEcalIso, prefix + "HEEPEcalIso" + suffix );
  iEvent.put( heepTrkIso, prefix + "HEEPTrkIso" + suffix );
  iEvent.put( trackMatch, prefix + "TrackMatch" + suffix );
  iEvent.put( dRTrack1, prefix + "DrTrack1" + suffix );
  iEvent.put( track1Eta, prefix + "Track1Eta" + suffix );
  iEvent.put( track1Phi, prefix + "Track1Phi" + suffix );
  iEvent.put( track1Pt, prefix + "Track1Pt" + suffix );
  iEvent.put( dRTrack2, prefix + "DrTrack2" + suffix );
  iEvent.put( track2Eta, prefix + "Track2Eta" + suffix );
  iEvent.put( track2Phi, prefix + "Track2Phi" + suffix );
  iEvent.put( track2Pt, prefix + "Track2Pt" + suffix );
  iEvent.put( scE1E9, elePrefix + "SCE1E9" + suffix );
  iEvent.put( scS4S1, elePrefix + "SCS4S1" + suffix );
  iEvent.put( sckOutOfTime, elePrefix + "SCkOutOfTime" + suffix );
  iEvent.put( scEcalIso, elePrefix + "SCEcalIso" + suffix );
  iEvent.put( scHEEPEcalIso, elePrefix + "SCHEEPEcalIso" + suffix );
  iEvent.put( scHEEPTrkIso, elePrefix + "SCHEEPTrkIso" + suffix );
}
