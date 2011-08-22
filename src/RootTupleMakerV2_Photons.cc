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
  produces <std::vector<double> > ( prefix + "HcalIsoDR04FullCone" + suffix );
  produces <std::vector<double> > ( prefix + "TrkIsoSolidDR04" + suffix );
  produces <std::vector<double> > ( prefix + "TrkIsoHollowDR04" + suffix );
  produces <std::vector<double> > ( prefix + "EcalIsoDR03" + suffix );
  produces <std::vector<double> > ( prefix + "HcalIsoDR03" + suffix );
  produces <std::vector<double> > ( prefix + "HcalIsoDR03FullCone" + suffix );
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
  produces <std::vector<double> > ( prefix + "SMajMaj" + suffix );
  produces <std::vector<double> > ( prefix + "SMinMin" + suffix );
  produces <std::vector<double> > ( prefix + "Alpha" + suffix );
  produces <std::vector<double> > ( prefix + "SEtaEta" + suffix );
  produces <std::vector<double> > ( prefix + "SEtaPhi" + suffix );
  produces <std::vector<double> > ( prefix + "SPhiPhi" + suffix );
  produces <std::vector<double> > ( prefix + "E2OverE9" + suffix );
}


//------ user defined functions ------

inline float RootTupleMakerV2_Photons::recHitE( const  DetId id,  const EcalRecHitCollection &recHits )
{
  if ( id == DetId(0) ) {
    return 0;
  } else {
    EcalRecHitCollection::const_iterator it = recHits.find( id );
    if ( it != recHits.end() ) return (*it).energy();
  }
  return 0;
}

inline float RootTupleMakerV2_Photons::recHitE( const DetId id, const EcalRecHitCollection & recHits,int di, int dj )
{
  // in the barrel:   di = dEta   dj = dPhi
  // in the endcap:   di = dX     dj = dY

  DetId nid;
  if( id.subdetId() == EcalBarrel) nid = EBDetId::offsetBy( id, di, dj );
  else if( id.subdetId() == EcalEndcap) nid = EEDetId::offsetBy( id, di, dj );

  return ( nid == DetId(0) ? 0 : recHitE( nid, recHits ) );
}

inline float RootTupleMakerV2_Photons::recHitApproxEt(  const DetId id,  const EcalRecHitCollection &recHits )
{
  // for the time being works only for the barrel
  if ( id.subdetId() == EcalBarrel ) {
    return recHitE( id, recHits ) / cosh( EBDetId::approxEta( id ) );
  }
  return 0;
}

float RootTupleMakerV2_Photons::GetE2OverE9( const DetId id, const EcalRecHitCollection & recHits)
{ ///////////start calculating e2/e9
  ////http://cmslxr.fnal.gov/lxr/source/RecoLocalCalo/EcalRecAlgos/src/EcalSeverityLevelAlgo.cc#240
  // compute e2overe9
  //   | | | |
  //   +-+-+-+
  //   | |1|2|
  //   +-+-+-+
  //   | | | |
  //   1 - input hit,  2 - highest energy hit in a 3x3 around 1
  //   rechit 1 must have E_t > recHitEtThreshold
  //   rechit 2 must have E_t > recHitEtThreshold2
  //   function returns value of E2/E9 centered around 1 (E2=energy of hits 1+2) if energy of 1>2
  //   if energy of 2>1 and KillSecondHit is set to true, function returns value of E2/E9 centered around 2
  //   *provided* that 1 is the highest energy hit in a 3x3 centered around 2, otherwise, function returns 0
  
float recHitEtThreshold = 10.0;
float recHitEtThreshold2 = 1.0;
bool avoidIeta85=false;
bool KillSecondHit=true;

if ( id.subdetId() == EcalBarrel ) {
  EBDetId ebId( id );
  // avoid recHits at |eta|=85 where one side of the neighbours is missing
  if ( abs(ebId.ieta())==85 && avoidIeta85) return 0;
  // select recHits with Et above recHitEtThreshold
  float e1 = recHitE( id, recHits );
  float ete1=recHitApproxEt( id, recHits );
  // check that rechit E_t is above threshold
  if (ete1 < std::min(recHitEtThreshold,recHitEtThreshold2) ) return 0;
  if (ete1 < recHitEtThreshold && !KillSecondHit ) return 0;

  float e2=-1;
  float ete2=0;
  float s9 = 0;
  // coordinates of 2nd hit relative to central hit
  int e2eta=0;
  int e2phi=0;

  // LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 1
  for ( int deta = -1; deta <= +1; ++deta ) {
    for ( int dphi = -1; dphi <= +1; ++dphi ) {
      // compute 3x3 energy
      float etmp=recHitE( id, recHits, deta, dphi );
      s9 += etmp;
      EBDetId idtmp=EBDetId::offsetBy(id,deta,dphi);
      float eapproxet=recHitApproxEt( idtmp, recHits );
      // remember 2nd highest energy deposit (above threshold) in 3x3 array
      if (etmp>e2 && eapproxet>recHitEtThreshold2 && !(deta==0 && dphi==0)) {
	e2=etmp;
	ete2=eapproxet;
	e2eta=deta;
	e2phi=dphi;
      }
    }
  }

  if ( e1 == 0 )  return 0;
  // return 0 if 2nd hit is below threshold
  if ( e2 == -1 ) return 0;
  // compute e2/e9 centered around 1st hit
  float e2nd=e1+e2;
  float e2e9=0;

  if (s9!=0) e2e9=e2nd/s9;
  // if central hit has higher energy than 2nd hit
  //  return e2/e9 if 1st hit is above E_t threshold
  if (e1 > e2 && ete1>recHitEtThreshold) return e2e9;
  // if second hit has higher energy than 1st hit
  if ( e2 > e1 ) {
    // return 0 if user does not want to flag 2nd hit, or
    // hits are below E_t thresholds - note here we
    // now assume the 2nd hit to be the leading hit.

    if (!KillSecondHit || ete2<recHitEtThreshold || ete1<recHitEtThreshold2) {
      return 0;
    }
    else {
      // LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 2
      float s92nd=0;
      float e2nd_prime=0;
      int e2prime_eta=0;
      int e2prime_phi=0;

      EBDetId secondid=EBDetId::offsetBy(id,e2eta,e2phi);

      for ( int deta = -1; deta <= +1; ++deta ) {
	for ( int dphi = -1; dphi <= +1; ++dphi ) {

	  // compute 3x3 energy
	  float etmp=recHitE( secondid, recHits, deta, dphi );
	  s92nd += etmp;

	  if (etmp>e2nd_prime && !(deta==0 && dphi==0)) {
	    e2nd_prime=etmp;
	    e2prime_eta=deta;
	    e2prime_phi=dphi;
	  }
	}
      }
      // if highest energy hit around E2 is not the same as the input hit, return 0;
      if (!(e2prime_eta==-e2eta && e2prime_phi==-e2phi))
	{
	  return 0;
	}
      // compute E2/E9 around second hit
      float e2e9_2=0;
      if (s92nd!=0) e2e9_2=e2nd/s92nd;
      //   return the value of E2/E9 calculated around 2nd hit
      return e2e9_2;
    }
  }
 } else if ( id.subdetId() == EcalEndcap ) {
  // only used for EB at the moment
  return 0;
 }
return 0;
}


//------

void RootTupleMakerV2_Photons::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ecalIsoDR04  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIsoDR04  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIsoDR04FullCone  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkIsoSolidDR04  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkIsoHollowDR04  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ecalIsoDR03  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIsoDR03  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIsoDR03FullCone  ( new std::vector<double>()  );
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
  std::auto_ptr<std::vector<double> >  sMajMaj  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  sMinMin  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  alpha  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  sEtaEta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  sEtaPhi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  sPhiPhi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  e2OverE9  ( new std::vector<double>()  );

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
	  //get the original photon ( https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATFAQs#How_can_I_retrieve_the_reference )
	  edm::Ptr<reco::Candidate> originalRef = it->originalObjectRef();
	  const reco::Photon *originalPhoton = dynamic_cast<const reco::Photon *>(originalRef.get());

	  matchesElectron = ConversionTools::hasMatchedPromptElectron(originalPhoton->superCluster(),electrons,hConversions,bsHandle->position());
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
      double sMajMaj_      = -99.;
      double sMinMin_      = -99.;
      double alpha_        = -99.;
      double sEtaEta_      = it->sigmaEtaEta();
      double sEtaPhi_      = -99.;
      double sPhiPhi_      = -99.;
      double e2OverE9_     = -99.;

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
	  // cluster shape:
	  if(maxRH.second) {
	    Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*theSeed, *rechits);
	    std::vector<float> etaphimoments = EcalClusterTools::localCovariances(*theSeed, &(*rechits), &(*topology));
	    sMajMaj_ = moments.sMaj;
	    sMinMin_ = moments.sMin;
	    alpha_   = moments.alpha;
	    sEtaEta_ = etaphimoments[0];
	    sEtaPhi_ = etaphimoments[1];
	    sPhiPhi_ = etaphimoments[2];
	  }
	  if ( it->isEB() )
	    e2OverE9_ = GetE2OverE9(seedCrystalId,*rechits);
	}
      else
	{
	  if( !ecalhitseb.isValid() )
	    edm::LogError("RootTupleMakerV2_PhotonsError") << "Error! Can't get the product " << ecalRecHitsEBInputTag; //FIXME
	  if( !ecalhitsee.isValid() )
	    edm::LogError("RootTupleMakerV2_PhotonsError") << "Error! Can't get the product " << ecalRecHitsEEInputTag; //FIXME
	}

      //-----------------------------------------------------------------      
      
      eta->push_back( it->eta() );
      phi->push_back( it->phi() );
      pt->push_back( it->pt() );
      energy->push_back( it->energy() );
      ecalIsoDR04->push_back( it->ecalRecHitSumEtConeDR04() );
      hcalIsoDR04->push_back( it->hcalTowerSumEtConeDR04() );
      hcalIsoDR04FullCone->push_back( it->hcalTowerSumEtConeDR04()
				      + ( it->hadronicOverEm() 
					  * it->superCluster()->energy() 
					  / cosh(it->superCluster()->eta())
					  )				     
				      );
      trkIsoHollowDR04->push_back( it->trkSumPtHollowConeDR04() );
      trkIsoSolidDR04->push_back( it->trkSumPtSolidConeDR04() );
      ecalIsoDR03->push_back( it->ecalRecHitSumEtConeDR03() );
      hcalIsoDR03->push_back( it->hcalTowerSumEtConeDR03() );
      hcalIsoDR03FullCone->push_back( it->hcalTowerSumEtConeDR03()
				      + ( it->hadronicOverEm() 
					  * it->superCluster()->energy() 
					  / cosh(it->superCluster()->eta())
					  )				     
				      );
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
      sMajMaj->push_back( sMajMaj_ );
      sMinMin->push_back( sMinMin_ );
      alpha->push_back( alpha_ );
      sEtaEta->push_back( sEtaEta_ );
      sEtaPhi->push_back( sEtaPhi_ );
      sPhiPhi->push_back( sPhiPhi_ );
      e2OverE9->push_back( e2OverE9_ );
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
  iEvent.put( hcalIsoDR04FullCone, prefix + "HcalIsoDR04FullCone" + suffix );
  iEvent.put( trkIsoHollowDR04, prefix + "TrkIsoHollowDR04" + suffix );
  iEvent.put( trkIsoSolidDR04, prefix + "TrkIsoSolidDR04" + suffix );
  iEvent.put( ecalIsoDR03, prefix + "EcalIsoDR03" + suffix );
  iEvent.put( hcalIsoDR03, prefix + "HcalIsoDR03" + suffix );
  iEvent.put( hcalIsoDR03FullCone, prefix + "HcalIsoDR03FullCone" + suffix );
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
  //variables from ecal rechits
  iEvent.put( timeSeed, prefix + "TimeSeed" + suffix );
  iEvent.put( e4SwissCross, prefix + "E4SwissCross" + suffix );
  iEvent.put( sMajMaj, prefix + "SMajMaj" + suffix );
  iEvent.put( sMinMin, prefix + "SMinMin" + suffix );
  iEvent.put( alpha, prefix + "Alpha" + suffix );
  iEvent.put( sEtaEta, prefix + "SEtaEta" + suffix );
  iEvent.put( sEtaPhi, prefix + "SEtaPhi" + suffix );
  iEvent.put( sPhiPhi, prefix + "SPhiPhi" + suffix );
  iEvent.put( e2OverE9, prefix + "E2OverE9" + suffix );
}

