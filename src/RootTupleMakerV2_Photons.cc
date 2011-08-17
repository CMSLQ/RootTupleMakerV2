#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Photons.h"
#include "Leptoquarks/RootTupleMakerV2/interface/PatUtilities.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h" 
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

RootTupleMakerV2_Photons::RootTupleMakerV2_Photons(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize")),
    beamSpotInputTag(iConfig.getParameter<edm::InputTag>    ("BeamSpotInputTag")),
    conversionsInputTag(iConfig.getParameter<edm::InputTag> ("ConversionsInputTag")),
    electronsInputTag(iConfig.getParameter<edm::InputTag>   ("ElectronsInputTag")),
    ecalRecHitsEBInputTag(iConfig.getParameter<edm::InputTag>   ("EcalRecHitsEBInputTag")),
    ecalRecHitsEEInputTag(iConfig.getParameter<edm::InputTag>   ("EcalRecHitsEEInputTag"))
{
  produces <std::vector<double> > ( prefix + "Eta" + suffix );
  produces <std::vector<double> > ( prefix + "Phi" + suffix );
  produces <std::vector<double> > ( prefix + "Pt" + suffix );
  produces <std::vector<double> > ( prefix + "Energy" + suffix );
  produces <std::vector<double> > ( prefix + "EcalIsoDR04" + suffix );
  produces <std::vector<double> > ( prefix + "HcalIsoDR04" + suffix );
  produces <std::vector<double> > ( prefix + "TrkIsoSolidDR04" + suffix );
  produces <std::vector<double> > ( prefix + "TrkIsoHollowDR04" + suffix );
  produces <std::vector<double> > ( prefix + "EcalIsoDR03" + suffix );
  produces <std::vector<double> > ( prefix + "HcalIsoDR03" + suffix );
  produces <std::vector<double> > ( prefix + "TrkIsoSolidDR03" + suffix );
  produces <std::vector<double> > ( prefix + "TrkIsoHollowDR03" + suffix );
  produces <std::vector<double> > ( prefix + "HoE" + suffix );
  produces <std::vector<double> > ( prefix + "SigmaIEtaIEta" + suffix );
  produces <std::vector<bool> >   ( prefix + "HasPixelSeed" + suffix );
  produces <std::vector<double> > ( prefix + "SCseedEnergy" + suffix );
  produces <std::vector<double> > ( prefix + "SCenergy" + suffix );
  produces <std::vector<double> > ( prefix + "SCeta" + suffix );
  produces <std::vector<double> > ( prefix + "SCphi" + suffix );
  produces <std::vector<double> > ( prefix + "E3x3" + suffix );
  produces <std::vector<double> > ( prefix + "E5x5" + suffix );
  produces <std::vector<bool> >   ( prefix + "HasMatchedPromptEle" + suffix );
  produces <std::vector<bool> >   ( prefix + "IsEBGap" + suffix );
  produces <std::vector<bool> >   ( prefix + "IsEEGap" + suffix );
  produces <std::vector<bool> >   ( prefix + "IsEBEEGap" + suffix );
  produces <std::vector<bool> >   ( prefix + "HasMatchedConvPhot" + suffix );
  produces <std::vector<int> >    ( prefix + "NTracksConvPhot" + suffix );
  produces <std::vector<double> > ( prefix + "PairInvariantMassConvPhot" + suffix );
  produces <std::vector<double> > ( prefix + "PairCotThetaSeparationConvPhot" + suffix );
  produces <std::vector<double> > ( prefix + "PairMomentumxConvPhot" + suffix );
  produces <std::vector<double> > ( prefix + "PairMomentumyConvPhot" + suffix );
  produces <std::vector<double> > ( prefix + "PairMomentumzConvPhot" + suffix );
  produces <std::vector<double> > ( prefix + "Chi2ConvPhot" + suffix );
  produces <std::vector<double> > ( prefix + "NDofConvPhot" + suffix );
  produces <std::vector<double> > ( prefix + "XVtxConvPhot" + suffix );
  produces <std::vector<double> > ( prefix + "YVtxConvPhot" + suffix );
  produces <std::vector<double> > ( prefix + "ZVtxConvPhot" + suffix );
  produces <std::vector<double> > ( prefix + "EOverPConvPhot" + suffix );
  produces <std::vector<double> > ( prefix + "DistOfMinApproachConvPhot" + suffix );
  produces <std::vector<double> > ( prefix + "DPhiTracksAtVtxConvPhot" + suffix );
  produces <std::vector<double> > ( prefix + "TimeSeed" + suffix );
  produces <std::vector<double> > ( prefix + "E4SwissCross" + suffix );
}

void RootTupleMakerV2_Photons::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ecalIsoDR04  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIsoDR04  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkIsoSolidDR04  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkIsoHollowDR04  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ecalIsoDR03  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIsoDR03  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkIsoSolidDR03  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkIsoHollowDR03  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hoe  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  sigmaIetaIeta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<bool> >    hasPixelSeed  ( new std::vector<bool>()  );
  std::auto_ptr<std::vector<double> >  SCseedEnergy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  SCenergy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  SCeta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  SCphi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  E3x3  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  E5x5  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<bool> >    hasMatchedPromptEle  ( new std::vector<bool>()  );
  std::auto_ptr<std::vector<bool> >    isEBGap  ( new std::vector<bool>()  );
  std::auto_ptr<std::vector<bool> >    isEEGap  ( new std::vector<bool>()  );
  std::auto_ptr<std::vector<bool> >    isEBEEGap  ( new std::vector<bool>()  );
  std::auto_ptr<std::vector<bool> >    hasMatchedConvPhot  ( new std::vector<bool>()  );
  std::auto_ptr<std::vector<int> >     nTracksConvPhot  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  pairInvariantMassConvPhot  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pairCotThetaSeparationConvPhot  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pairMomentumxConvPhot  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pairMomentumyConvPhot  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pairMomentumzConvPhot  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  chi2ConvPhot  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  nDofConvPhot  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  xVtxConvPhot  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  yVtxConvPhot  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  zVtxConvPhot  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  eOverPConvPhot  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  distOfMinApproachConvPhot  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  dPhiTracksAtVtxConvPhot  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  timeSeed  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  e4SwissCross  ( new std::vector<double>()  );

  //-----------------------------------------------------------------

  edm::Handle<std::vector<pat::Photon> > photons;
  iEvent.getByLabel(inputTag, photons);

  edm::Handle<reco::BeamSpot> bsHandle; 
  iEvent.getByLabel(beamSpotInputTag, bsHandle); 

  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel(conversionsInputTag, hConversions); 

  edm::Handle<reco::GsfElectronCollection> electrons;
  iEvent.getByLabel(electronsInputTag, electrons);

  edm::Handle<EBRecHitCollection> ecalhitseb;
  const EBRecHitCollection* rhitseb=0;
  iEvent.getByLabel(ecalRecHitsEBInputTag, ecalhitseb);
  rhitseb = ecalhitseb.product(); // get a ptr to the product

  edm::Handle<EERecHitCollection> ecalhitsee;
  const EERecHitCollection* rhitsee=0;
  iEvent.getByLabel(ecalRecHitsEEInputTag, ecalhitsee);
  rhitsee = ecalhitsee.product(); // get a ptr to the product

  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);
  const CaloTopology *topology = pTopology.product();

  if(photons.isValid()) {
    edm::LogInfo("RootTupleMakerV2_PhotonsInfo") << "Total # Photons: " << photons->size();

    for( std::vector<pat::Photon>::const_iterator it = photons->begin(); it != photons->end(); ++it ) {
      // exit from loop when you reach the required number of photons
      if(eta->size() >= maxSize)
        break;

      bool matchesElectron = false; 
      if( hConversions.isValid() && bsHandle.isValid() && electrons.isValid()) 
	{
	  matchesElectron = ConversionTools::hasMatchedPromptElectron(it->superCluster(),electrons,hConversions,bsHandle->position());
	  //See https://indico.cern.ch/getFile.py/access?contribId=6&resId=0&materialId=slides&confId=129730
          //and https://hypernews.cern.ch/HyperNews/CMS/get/egamma/999.html ( N.3 )
	}
      else 
	{
	  if( !bsHandle.isValid() )
	    edm::LogError("RootTupleMakerV2_PhotonsError") << "Error! Can't get the product " << beamSpotInputTag;
	  if( !hConversions.isValid() )
	    edm::LogError("RootTupleMakerV2_PhotonsError") << "Error! Can't get the product " << conversionsInputTag;
	  if( !electrons.isValid() )
	    edm::LogError("RootTupleMakerV2_PhotonsError") << "Error! Can't get the product " << electronsInputTag;	  
	}

      //photon conversions
      bool   hasMatchedConvPhot_                 =   false; 
      int    nTracksConvPhot_                    =   0; 
      double pairInvariantMassConvPhot_          =   -9999.; 
      double pairCotThetaSeparationConvPhot_     =   -9999.; 
      double pairMomentumxConvPhot_              =   -9999.; 
      double pairMomentumyConvPhot_              =   -9999.; 
      double pairMomentumzConvPhot_              =   -9999.; 
      double chi2ConvPhot_                       =   -9999.; 
      double nDofConvPhot_                       =   -9999.; 
      double xVtxConvPhot_                       =   -9999.; 
      double yVtxConvPhot_                       =   -9999.; 
      double zVtxConvPhot_                       =   -9999.; 
      double eOverPConvPhot_                     =   -9999.; 
      double distOfMinApproachConvPhot_          =   -9999.; 
      double dPhiTracksAtVtxConvPhot_            =   -9999.; 

      if( hConversions.isValid() && bsHandle.isValid() ) 
	{
	  reco::ConversionRef conv = ConversionTools::matchedConversion(*it->superCluster(),hConversions,bsHandle->position());
	  hasMatchedConvPhot_ = conv.isNonnull();
          //and https://hypernews.cern.ch/HyperNews/CMS/get/egamma/999.html ( N.2 )
	  
	  if( hasMatchedConvPhot_ )
	    {
	      nTracksConvPhot_                    = conv->nTracks();
	      pairInvariantMassConvPhot_          = conv->pairInvariantMass();
	      pairCotThetaSeparationConvPhot_     = conv->pairCotThetaSeparation();
	      pairMomentumxConvPhot_              = conv->pairMomentum().x();
	      pairMomentumyConvPhot_              = conv->pairMomentum().y();
	      pairMomentumzConvPhot_              = conv->pairMomentum().z();
	      if( conv->conversionVertex().isValid() )
		{
		  chi2ConvPhot_                       = conv->conversionVertex().chi2();
		  nDofConvPhot_                       = conv->conversionVertex().ndof();
		  xVtxConvPhot_                       = conv->conversionVertex().x();
		  yVtxConvPhot_                       = conv->conversionVertex().y();
		  zVtxConvPhot_                       = conv->conversionVertex().z();
		}
	      eOverPConvPhot_                     = conv->EoverP();
	      distOfMinApproachConvPhot_          = conv->distOfMinimumApproach();
	      dPhiTracksAtVtxConvPhot_            = conv->dPhiTracksAtVtx();	      
	    }
	}
      else 
	{
	  if( !bsHandle.isValid() )
	    edm::LogError("RootTupleMakerV2_PhotonsError") << "Error! Can't get the product " << beamSpotInputTag;
	  if( !hConversions.isValid() )
	    edm::LogError("RootTupleMakerV2_PhotonsError") << "Error! Can't get the product " << conversionsInputTag;
	  if( !electrons.isValid() )
	    edm::LogError("RootTupleMakerV2_PhotonsError") << "Error! Can't get the product " << electronsInputTag;	  
	}

      //variables from ecal rechits collection
      double timeSeed_     = 0.;     
      double e4SwissCross_ = -1;     
//       double sMajMaj_      = -99.;
//       double sMinMin_      = -99.;
//       double alpha_        = -99.;
//       double sEtaEta_      = it->sigmaEtaEta();
//       double sEtaPhi_      = -99.;
//       double sPhiPhi_      = -99.;
      const reco::CaloClusterPtr theSeed = it->superCluster()->seed(); 
      const EBRecHitCollection* rechits = ( it->isEB()) ? rhitseb : rhitsee;
      if( ecalhitseb.isValid() && ecalhitsee.isValid() )
	{
	  std::pair<DetId, float> maxRH = EcalClusterTools::getMaximum( *theSeed, &(*rechits) );
	  
	  // timing:
	  DetId seedCrystalId = maxRH.first;
	  if(maxRH.second) 
	    {
	      EcalRecHitCollection::const_iterator seedRH = rechits->find(seedCrystalId);
	      timeSeed_ = seedRH->time();
	    }
	  // swiss cross:
	  e4SwissCross_ = (!it->isEB()) ? -1 :
	    ( EcalClusterTools::eLeft( *theSeed, &(*rechits), topology ) +
	      EcalClusterTools::eRight( *theSeed, &(*rechits), topology ) +
	      EcalClusterTools::eTop( *theSeed, &(*rechits), topology ) +
	      EcalClusterTools::eBottom( *theSeed, &(*rechits), topology ) );	  
// 	  //cluster shape: 	  //FIXME
// 	  if(maxRH.second) {
// 	    Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*theSeed, *rechits);
// 	    std::vector<float> etaphimoments = EcalClusterTools::localCovariances(*theSeed, &(*rechits), &(*topology));
// 	    sMajMaj_ = moments.sMaj;
// 	    sMinMin_ = moments.sMin;
// 	    alpha_   = moments.alpha;
// 	    sEtaEta_ = etaphimoments[0];
// 	    sEtaPhi_ = etaphimoments[1];
// 	    sPhiPhi_ = etaphimoments[2];
// 	  }else{
// 	    sMajMaj_ = -100.;
// 	    sMinMin_ = -100.;
// 	    alpha_   = -100.;
// 	    sEtaEta_ = it->sigmaEtaEta();//-100.;
// 	    sEtaPhi_ = -100.;
// 	    sPhiPhi_ = -100.;
// 	  }
// 	  E2OverE9Phot_=-99.;
// 	  if ( it->isEB() )
// 	    E2OverE9Phot[nPhot] = GetE2OverE9(seedCrystalId,*rechits);
	}
      else
	{
	  if( !ecalhitseb.isValid() )
	    edm::LogError("RootTupleMakerV2_PhotonsError") << "Error! Can't get the product " << beamSpotInputTag; //FIXME
	  if( !ecalhitsee.isValid() )
	    edm::LogError("RootTupleMakerV2_PhotonsError") << "Error! Can't get the product " << beamSpotInputTag; //FIXME
	}

      //-----------------------------------------------------------------      
      
      eta->push_back( it->eta() );
      phi->push_back( it->phi() );
      pt->push_back( it->pt() );
      energy->push_back( it->energy() );
      ecalIsoDR04->push_back( it->ecalRecHitSumEtConeDR04() );
      hcalIsoDR04->push_back( it->hcalTowerSumEtConeDR04() );
      trkIsoHollowDR04->push_back( it->trkSumPtHollowConeDR04() );
      trkIsoSolidDR04->push_back( it->trkSumPtSolidConeDR04() );
      ecalIsoDR03->push_back( it->ecalRecHitSumEtConeDR03() );
      hcalIsoDR03->push_back( it->hcalTowerSumEtConeDR03() );
      trkIsoHollowDR03->push_back( it->trkSumPtHollowConeDR03() );
      trkIsoSolidDR03->push_back( it->trkSumPtSolidConeDR03() );
      hoe->push_back( it->hadronicOverEm() );
      sigmaIetaIeta->push_back( it->sigmaIetaIeta() );
      hasPixelSeed->push_back( it->hasPixelSeed() );
      SCseedEnergy->push_back( it->superCluster()->seed()->energy() );
      SCenergy->push_back( it->superCluster()->energy() );
      SCeta->push_back( it->superCluster()->eta() );
      SCphi->push_back( it->superCluster()->phi() );
      E3x3->push_back( it->e3x3() );
      E5x5->push_back( it->e5x5() );
      hasMatchedPromptEle->push_back( matchesElectron );
      isEBGap->push_back( it->isEBGap() );
      isEEGap->push_back( it->isEEGap() );
      isEBEEGap->push_back( it->isEBEEGap() );
      //photon conversions
      hasMatchedConvPhot->push_back( hasMatchedConvPhot_ ); 
      nTracksConvPhot->push_back( nTracksConvPhot_ ); 
      pairInvariantMassConvPhot->push_back( pairInvariantMassConvPhot_ ); 
      pairCotThetaSeparationConvPhot->push_back( pairCotThetaSeparationConvPhot_ ); 
      pairMomentumxConvPhot->push_back( pairMomentumxConvPhot_ ); 
      pairMomentumyConvPhot->push_back( pairMomentumyConvPhot_ ); 
      pairMomentumzConvPhot->push_back( pairMomentumzConvPhot_ ); 
      chi2ConvPhot->push_back( chi2ConvPhot_ ); 
      nDofConvPhot->push_back( nDofConvPhot_ ); 
      xVtxConvPhot->push_back( xVtxConvPhot_ ); 
      yVtxConvPhot->push_back( yVtxConvPhot_ ); 
      zVtxConvPhot->push_back( zVtxConvPhot_ ); 
      eOverPConvPhot->push_back( eOverPConvPhot_ ); 
      distOfMinApproachConvPhot->push_back( distOfMinApproachConvPhot_ ); 
      dPhiTracksAtVtxConvPhot->push_back( dPhiTracksAtVtxConvPhot_ ); 
      //variables from ecal rechits
      timeSeed->push_back( timeSeed_ );
      e4SwissCross->push_back( e4SwissCross_ );
    }
  } else {
    edm::LogError("RootTupleMakerV2_PhotonsError") << "Error! Can't get the product " << inputTag;
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( pt, prefix + "Pt" + suffix );
  iEvent.put( energy, prefix + "Energy" + suffix );
  iEvent.put( ecalIsoDR04, prefix + "EcalIsoDR04" + suffix );
  iEvent.put( hcalIsoDR04, prefix + "HcalIsoDR04" + suffix );
  iEvent.put( trkIsoHollowDR04, prefix + "TrkIsoHollowDR04" + suffix );
  iEvent.put( trkIsoSolidDR04, prefix + "TrkIsoSolidDR04" + suffix );
  iEvent.put( ecalIsoDR03, prefix + "EcalIsoDR03" + suffix );
  iEvent.put( hcalIsoDR03, prefix + "HcalIsoDR03" + suffix );
  iEvent.put( trkIsoHollowDR03, prefix + "TrkIsoHollowDR03" + suffix );
  iEvent.put( trkIsoSolidDR03, prefix + "TrkIsoSolidDR03" + suffix );
  iEvent.put( hoe, prefix + "HoE" + suffix );
  iEvent.put( sigmaIetaIeta, prefix + "SigmaIEtaIEta" + suffix );
  iEvent.put( hasPixelSeed, prefix + "HasPixelSeed" + suffix );
  iEvent.put( SCseedEnergy, prefix + "SCseedEnergy" + suffix );
  iEvent.put( SCenergy, prefix + "SCenergy" + suffix );
  iEvent.put( SCeta, prefix + "SCeta" + suffix );
  iEvent.put( SCphi, prefix + "SCphi" + suffix );
  iEvent.put( E3x3, prefix + "E3x3" + suffix );
  iEvent.put( E5x5, prefix + "E5x5" + suffix );
  iEvent.put( hasMatchedPromptEle, prefix + "HasMatchedPromptEle" + suffix );
  iEvent.put( isEBGap, prefix + "IsEBGap" + suffix );
  iEvent.put( isEEGap, prefix + "IsEEGap" + suffix );
  iEvent.put( isEBEEGap, prefix + "IsEBEEGap" + suffix );
  //photon conversions
  iEvent.put( hasMatchedConvPhot, prefix + "HasMatchedConvPhot" + suffix );
  iEvent.put( nTracksConvPhot, prefix + "NTracksConvPhot" + suffix );
  iEvent.put( pairInvariantMassConvPhot, prefix + "PairInvariantMassConvPhot" + suffix );
  iEvent.put( pairCotThetaSeparationConvPhot, prefix + "PairCotThetaSeparationConvPhot" + suffix );
  iEvent.put( pairMomentumxConvPhot, prefix + "PairMomentumxConvPhot" + suffix );
  iEvent.put( pairMomentumyConvPhot, prefix + "PairMomentumyConvPhot" + suffix );
  iEvent.put( pairMomentumzConvPhot, prefix + "PairMomentumzConvPhot" + suffix );
  iEvent.put( chi2ConvPhot, prefix + "Chi2ConvPhot" + suffix );
  iEvent.put( nDofConvPhot, prefix + "NDofConvPhot" + suffix );
  iEvent.put( xVtxConvPhot, prefix + "XVtxConvPhot" + suffix );
  iEvent.put( yVtxConvPhot, prefix + "YVtxConvPhot" + suffix );
  iEvent.put( zVtxConvPhot, prefix + "ZVtxConvPhot" + suffix );
  iEvent.put( eOverPConvPhot, prefix + "EOverPConvPhot" + suffix );
  iEvent.put( distOfMinApproachConvPhot, prefix + "DistOfMinApproachConvPhot" + suffix );
  iEvent.put( dPhiTracksAtVtxConvPhot, prefix + "DPhiTracksAtVtxConvPhot" + suffix );
  iEvent.put( timeSeed, prefix + "TimeSeed" + suffix );
  iEvent.put( e4SwissCross, prefix + "E4SwissCross" + suffix );
}
