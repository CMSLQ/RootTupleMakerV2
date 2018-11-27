//------------------------------------------------------------------------
// To do:
// - Muon isolation (for overlaps) is done with relative isolation.  Is this correct?
//------------------------------------------------------------------------

//------------------------------------------------------------------------
// Include files
//------------------------------------------------------------------------

#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_Electrons.h"
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
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/GeometryVector/interface/Vector3DBase.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

//------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------

RootTupleMakerV2_Electrons::RootTupleMakerV2_Electrons(const edm::ParameterSet& iConfig) :
  electronInputToken_          (consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("InputTag"))),
  vtxInputToken_               (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("VertexInputTag"))), 
  rhoInputToken_               (consumes<double>(iConfig.getParameter<edm::InputTag>("RhoInputTag"))),
  electronVetoIdMapToken_      (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ElectronVetoIdMap"))),
  electronLooseIdMapToken_     (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ElectronLooseIdMap"))),
  electronMediumIdMapToken_    (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ElectronMediumIdMap"))),
  electronTightIdMapToken_     (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ElectronTightIdMap"))),
  electronHLTPreselectionMapToken_ (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ElectronHLTPreselectionMap"))),
  electronHEEPIdMapToken_      (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ElectronHEEPIdMap"))),
  electronMVAIdWP80MapToken_   (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ElectronMVAIdWP80Map"))),
  electronMVAIdWP90MapToken_   (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ElectronMVAIdWP90Map"))),
  electronMVAIdHZZMapToken_    (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ElectronMVAIdHZZMap"))),
  eleVetoIdCutFlowResultMapToken_   (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleVetoIdCutFlowResultMap"))),
  eleLooseIdCutFlowResultMapToken_  (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleLooseIdCutFlowResultMap"))),
  eleMediumIdCutFlowResultMapToken_ (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleMediumIdCutFlowResultMap"))),
  eleTightIdCutFlowResultMapToken_  (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleTightIdCutFlowResultMap"))),
  eleHLTPreselectionCutFlowResultMapToken_  (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleHLTPreselectionCutFlowResultMap"))),
  eleHEEPIdCutFlowResultMapToken_   (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdCutFlowResultMap"))),
  eleMVAIdWP80CutFlowResultMapToken_   (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("ElectronMVAIdWP80CutFlowResultMap"))),
  eleMVAIdWP90CutFlowResultMapToken_   (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("ElectronMVAIdWP90CutFlowResultMap"))),
  eleMVAIdHZZCutFlowResultMapToken_    (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("ElectronMVAIdHZZCutFlowResultMap"))),
  heep70trkIsolMapToken_            (consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("heep70trkIsolMap"))),
  beamSpotInputToken_               (consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpotInputTag"))), 
  ebReducedRecHitsInputToken_       (consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBReducedRecHitsInputTag"))),
  eeReducedRecHitsInputToken_       (consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EEReducedRecHitsInputTag"))),
  electronIso                  (iConfig.getParameter<double>       ("ElectronIso"              )),
  muonPt                       (iConfig.getParameter<double>       ("MuonPt"                   )),
  muonIso                      (iConfig.getParameter<double>       ("MuonIso"                  )),
  muonID                       (iConfig.getParameter<std::string>  ("MuonID"                   )),
  prefix                       (iConfig.getParameter<std::string>  ("Prefix"                   )),
  suffix                       (iConfig.getParameter<std::string>  ("Suffix"                   )),
  maxSize                      (iConfig.getParameter<unsigned int> ("MaxSize"                  ))
 {
  
  //------------------------------------------------------------------------
  // What variables will this producer push into the event?
  //------------------------------------------------------------------------
  
  // Kinematic variables
  produces <std::vector<float> > ( prefix + "Eta"                      + suffix );
  produces <std::vector<float> > ( prefix + "Phi"                      + suffix );
  produces <std::vector<float> > ( prefix + "Pt"                       + suffix );
  produces <std::vector<float> > ( prefix + "PtHeep"                   + suffix );
  produces <std::vector<float> > ( prefix + "Energy"                   + suffix );
  produces <std::vector<float> > ( prefix + "CaloEnergy"               + suffix );
  produces <std::vector<float> > ( prefix + "EcalEnergy"               + suffix );
  produces <std::vector<int> >   ( prefix + "Charge"                   + suffix );
  produces <std::vector<float> > ( prefix + "HoE"                      + suffix );
								        
  // Supercluster kinematic variables				        
  produces <std::vector<float> > ( prefix + "ESuperClusterOverP"       + suffix );
  produces <std::vector<float> > ( prefix + "SCEta"                    + suffix );
  produces <std::vector<float> > ( prefix + "SCSeedEta"                + suffix );
  produces <std::vector<float> > ( prefix + "SCPhi"                    + suffix );
  produces <std::vector<float> > ( prefix + "SCPt"                     + suffix );
  produces <std::vector<float> > ( prefix + "SCRawEnergy"              + suffix );
  produces <std::vector<float> > ( prefix + "SCEnergy"                 + suffix );
  produces <std::vector<float> > ( prefix + "SCSeedCryEnergy"          + suffix );
  produces <std::vector<bool> >  ( prefix + "SCSeedCryIsBarrel"        + suffix );

  // ID information
  produces <std::vector<int> >     ( prefix + "PassId"                     + suffix );
  produces <std::vector<bool> >    ( prefix + "PassEGammaIDVeto"           + suffix );
  produces <std::vector<bool> >    ( prefix + "PassEGammaIDLoose"          + suffix );
  produces <std::vector<bool> >    ( prefix + "PassEGammaIDMedium"         + suffix );
  produces <std::vector<bool> >    ( prefix + "PassEGammaIDTight"          + suffix );
  produces <std::vector<bool> >    ( prefix + "PassEGammaHLTPreselection"          + suffix );
  produces <std::vector<bool> >    ( prefix + "PassHEEPID"                 + suffix );
  produces <std::vector<bool> >    ( prefix + "PassMVAIDWP80"              + suffix );
  produces <std::vector<bool> >    ( prefix + "PassMVAIDWP90"              + suffix );
  produces <std::vector<bool> >    ( prefix + "PassMVAIDHZZ"               + suffix );
  //produces <std::vector<std::string> >    ( prefix + "CutFlowNamesEGammaIDVeto"   + suffix );
  //produces <std::vector<std::string> >    ( prefix + "CutFlowNamesEGammaIDLoose"  + suffix );
  //produces <std::vector<std::string> >    ( prefix + "CutFlowNamesEGammaIDMedium" + suffix );
  //produces <std::vector<std::string> >    ( prefix + "CutFlowNamesEGammaIDTight"  + suffix );
  produces <std::vector<std::string> >    ( prefix + "CutFlowNamesEGammaIDHEEP"   + suffix );
  //produces <std::vector<std::string> >    ( prefix + "CutFlowHashesEGammaIDVeto"  + suffix );
  //produces <std::vector<std::string> >    ( prefix + "CutFlowHashesEGammaIDLoose" + suffix );
  //produces <std::vector<std::string> >    ( prefix + "CutFlowHashesEGammaIDMedium"+ suffix );
  //produces <std::vector<std::string> >    ( prefix + "CutFlowHashesEGammaIDTight" + suffix );
  produces <std::vector<std::string> >    ( prefix + "CutFlowHashesEGammaIDHEEP"  + suffix );
  produces <std::vector<float> >   ( prefix + "RhoIsoHEEP"                 + suffix );

  // Does this electron overlap with a muon?			        
  produces <std::vector<int> >    ( prefix + "Overlaps"                 + suffix );
								        
  // Number of Brems = number of basic clusters minus one	        
  produces <std::vector<int> >    ( prefix + "NumberOfBrems"            + suffix );

  // Is this ECAL driven? Or PFlow?  				        
  produces <std::vector<bool> >   ( prefix + "HasEcalDrivenSeed"        + suffix );
  produces <std::vector<bool> >   ( prefix + "IsEcalDriven"             + suffix );
  produces <std::vector<bool> >   ( prefix + "HasTrackerDrivenSeed"     + suffix );

  // Charge consistency variables - Ferdinando Giordano
  produces <std::vector<bool> >   ( prefix + "GsfCtfScPixCharge"        + suffix );
  produces <std::vector<bool> >   ( prefix + "GsfScPixCharge"           + suffix );
  produces <std::vector<bool> >   ( prefix + "GsfCtfCharge"             + suffix );

  // EE or EB - Ferdinando Giordano
  produces <std::vector<bool> >   ( prefix + "IsEB"                     + suffix );
  produces <std::vector<bool> >   ( prefix + "IsEE"                     + suffix );
  
  // ECAL eta/phi vs tracker eta/phi				        
  produces <std::vector<float> > ( prefix + "DeltaPhiTrkSC"            + suffix );
  produces <std::vector<float> > ( prefix + "DeltaEtaTrkSC"            + suffix );
  produces <std::vector<float> > ( prefix + "DeltaEtaTrkSeedSC"        + suffix );
								        
  // Shower shape						        
  produces <std::vector<float> > ( prefix + "SigmaEtaEta"              + suffix );
  //produces <std::vector<float> > ( prefix + "SigmaIEtaIEta"            + suffix );
  produces <std::vector<float> > ( prefix + "Full5x5SigmaIEtaIEta"     + suffix );
  produces <std::vector<int> >    ( prefix + "Classif"                  + suffix );
  produces <std::vector<float> > ( prefix + "R9"                       + suffix );
  //produces <std::vector<float> > ( prefix + "E1x5OverE5x5"             + suffix );
  //produces <std::vector<float> > ( prefix + "E2x5OverE5x5"             + suffix );
  produces <std::vector<float> > ( prefix + "Full5x5E1x5OverE5x5"      + suffix );
  produces <std::vector<float> > ( prefix + "Full5x5E2x5OverE5x5"      + suffix );
  								        
  //// Isolation variables: PAT					        
	//							        
  //produces <std::vector<float> > ( prefix + "TrkIsoPAT"                + suffix );
  //produces <std::vector<float> > ( prefix + "EcalIsoPAT"               + suffix );
  //produces <std::vector<float> > ( prefix + "HcalIsoPAT"               + suffix );
  //produces <std::vector<float> > ( prefix + "RelIsoPAT"                + suffix );
  //produces <std::vector<int> >    ( prefix + "PassIsoPAT"               + suffix );

  // Isolation variables: particle flow
  produces <std::vector<float> > ( prefix + "PFChargedHadronIso03"     + suffix );
  produces <std::vector<float> > ( prefix + "PFNeutralHadronIso03"     + suffix );
  produces <std::vector<float> > ( prefix + "PFPhotonIso03"            + suffix );
  produces <std::vector<float> > ( prefix + "PFPUIso03"                + suffix );

  //produces <std::vector<float> > ( prefix + "PFChargedHadronIso04"     + suffix );
  //produces <std::vector<float> > ( prefix + "PFNeutralHadronIso04"     + suffix );
  //produces <std::vector<float> > ( prefix + "PFPhotonIso04"            + suffix );
  
  // Isolation variables: DR 0.3				        
  produces <std::vector<float> > ( prefix + "EcalIsoDR03"              + suffix );
  produces <std::vector<float> > ( prefix + "HcalIsoDR03"              + suffix );
  produces <std::vector<float> > ( prefix + "HcalIsoDR03FullCone"      + suffix );
  produces <std::vector<float> > ( prefix + "HcalIsoD1DR03"            + suffix );
  produces <std::vector<float> > ( prefix + "HcalIsoD2DR03"            + suffix );
  produces <std::vector<float> > ( prefix + "TrkIsoDR03"               + suffix );
  produces <std::vector<float> > ( prefix + "EcalPFClusterIso"               + suffix );
  produces <std::vector<float> > ( prefix + "HcalPFClusterIso"               + suffix );

  // Isolation for HEEP v7.0
  produces <std::vector<float> > ( prefix + "Heep70TrkIso"             + suffix );
				        
  // Conversion variables					        
  produces <std::vector<int> >    ( prefix + "MissingHits"              + suffix );
  produces <std::vector<int> >    ( prefix + "MissingHitsEG"            + suffix );
  produces <std::vector<float> > ( prefix + "Dist"                     + suffix );
  produces <std::vector<float> > ( prefix + "DCotTheta"                + suffix );
  produces <std::vector<float> > ( prefix + "Fbrem"                    + suffix );
  produces <std::vector<bool> >   ( prefix + "HasMatchedConvPhot"       + suffix );

  // Vertex and beamspot information
  produces <std::vector<int> >   ( prefix + "VtxIndex"                 + suffix );
  produces <std::vector<float> > ( prefix + "VtxDistXY"                + suffix );
  produces <std::vector<float> > ( prefix + "VtxDistZ"                 + suffix );
  produces <std::vector<float> > ( prefix + "LeadVtxDistXY"            + suffix );
  produces <std::vector<float> > ( prefix + "LeadVtxDistZ"             + suffix );
  produces <std::vector<float> > ( prefix + "PrimaryVertexDXY"         + suffix );
  produces <std::vector<float> > ( prefix + "PrimaryVertexDXYError"    + suffix );
  produces <std::vector<float> > ( prefix + "BeamSpotDXY"              + suffix );
  produces <std::vector<float> > ( prefix + "BeamSpotDXYError"         + suffix );

  // Track information
  produces <std::vector<float> > ( prefix + "TrackPx"                  + suffix );
  produces <std::vector<float> > ( prefix + "TrackPy"                  + suffix );
  produces <std::vector<float> > ( prefix + "TrackPz"                  + suffix );
  produces <std::vector<float> > ( prefix + "TrackPt"                  + suffix );
  produces <std::vector<float> > ( prefix + "TrackValidFractionOfHits" + suffix );
  produces <std::vector<float> > ( prefix + "NormalizedChi2" + suffix );

  //// Trigger matching: float electron

  //produces <std::vector<bool  > > ( prefix + "HLTfloatEleMatched"      + suffix );
  //produces <std::vector<float> > ( prefix + "HLTfloatEleMatchPt"      + suffix );
  //produces <std::vector<float> > ( prefix + "HLTfloatEleMatchEta"     + suffix );
  //produces <std::vector<float> > ( prefix + "HLTfloatEleMatchPhi"     + suffix );

  //// Trigger matching: Single electron

  //produces <std::vector<bool  > > ( prefix + "HLTSingleEleMatched"      + suffix );
  //produces <std::vector<float> > ( prefix + "HLTSingleEleMatchPt"      + suffix );
  //produces <std::vector<float> > ( prefix + "HLTSingleEleMatchEta"     + suffix );
  //produces <std::vector<float> > ( prefix + "HLTSingleEleMatchPhi"     + suffix );

  //// Trigger matching: Single electron (WP85)

  //produces <std::vector<bool  > > ( prefix + "HLTSingleEleWP85Matched"  + suffix );
  //produces <std::vector<float> > ( prefix + "HLTSingleEleWP85MatchPt"  + suffix );
  //produces <std::vector<float> > ( prefix + "HLTSingleEleWP85MatchEta" + suffix );
  //produces <std::vector<float> > ( prefix + "HLTSingleEleWP85MatchPhi" + suffix );

  //// Trigger matching: Ele+Jet+Jet

  //produces <std::vector<bool  > > ( prefix + "HLTEleJetJetMatched"  + suffix );
  //produces <std::vector<float> > ( prefix + "HLTEleJetJetMatchPt"  + suffix );
  //produces <std::vector<float> > ( prefix + "HLTEleJetJetMatchEta" + suffix );
  //produces <std::vector<float> > ( prefix + "HLTEleJetJetMatchPhi" + suffix );

  // Gen matching: status 3 only 
  produces <std::vector<float> > ( prefix + "MatchedGenParticlePt"   + suffix );
  produces <std::vector<float> > ( prefix + "MatchedGenParticleEta"  + suffix );
  produces <std::vector<float> > ( prefix + "MatchedGenParticlePhi"  + suffix );

  // for 03Feb2017 re-miniaod
  produces <bool> ( prefix + "EGammaGSFixedDupECALClusters" + suffix );
  produces <std::vector<uint32_t> > ( prefix + "EcalMultiAndGSGlobalRecHitEBHitsNotReplaced" + suffix );
  edm::InputTag dupClustersTag("particleFlowEGammaGSFixed:dupECALClusters");
  dupEcalClustersToken_ = consumes<bool>(dupClustersTag);
  edm::InputTag ecalMultiTag("ecalMultiAndGSGlobalRecHitEB:hitsNotReplaced");
  ecalMultiAndGSGlobalRecHitEBHitsNotReplacedToken_ = consumes<edm::EDCollection<DetId> >(ecalMultiTag);
}

//------------------------------------------------------------------------
// Event-per-event processing
//------------------------------------------------------------------------

void RootTupleMakerV2_Electrons::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //------------------------------------------------------------------------
  // Declare items to push into the event
  //------------------------------------------------------------------------

  // Kinematic variables

  std::auto_ptr<std::vector<float> >  eta                       ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  phi                       ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  pt                        ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  ptHeep                    ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  energy                    ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  caloEnergy                ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  ecalEnergy                ( new std::vector<float>()  );
  std::auto_ptr<std::vector<int> >     charge                    ( new std::vector<int>   ()  );
  std::auto_ptr<std::vector<float> >  hoe                       ( new std::vector<float>()  );

  // Supercluster kinematic variables

  std::auto_ptr<std::vector<float> >  eSuperClusterOverP        ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  scEta                     ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  scSeedEta                 ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  scPhi                     ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  scPt                      ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  scRawEnergy               ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  scEnergy                  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  scSeedCryEnergy           ( new std::vector<float>()  );
  std::auto_ptr<std::vector<bool> >    scSeedCryIsBarrel         ( new std::vector<bool>()  );

  // ID information
  std::auto_ptr<std::vector<int> >     passIds                   ( new std::vector<int>   ()  );
  std::auto_ptr<std::vector<bool> >    passEGammaIDVeto          ( new std::vector<bool>   ()  );
  std::auto_ptr<std::vector<bool> >    passEGammaIDLoose         ( new std::vector<bool>   ()  );
  std::auto_ptr<std::vector<bool> >    passEGammaIDMedium        ( new std::vector<bool>   ()  );
  std::auto_ptr<std::vector<bool> >    passEGammaIDTight         ( new std::vector<bool>   ()  );
  std::auto_ptr<std::vector<bool> >    passEGammaHLTPreselection ( new std::vector<bool>   ()  );
  std::auto_ptr<std::vector<bool> >    passHEEPID                ( new std::vector<bool>   ()  ); 
  std::auto_ptr<std::vector<bool> >    passMVAIDWP80             ( new std::vector<bool>   ()  ); 
  std::auto_ptr<std::vector<bool> >    passMVAIDWP90             ( new std::vector<bool>   ()  ); 
  std::auto_ptr<std::vector<bool> >    passMVAIDHZZ              ( new std::vector<bool>   ()  ); 
  std::auto_ptr<std::vector<std::string> > cutFlowNamesEGammaIDVeto   ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowNamesEGammaIDLoose  ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowNamesEGammaIDMedium ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowNamesEGammaIDTight  ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowNamesEGammaHLTPreselection  ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowNamesEGammaIDHEEP   ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowNamesEGammaIDMVAWP80   ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowNamesEGammaIDMVAWP90   ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowNamesEGammaIDMVAHZZ    ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowHashesEGammaIDVeto   ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowHashesEGammaIDLoose  ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowHashesEGammaIDMedium ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowHashesEGammaIDTight  ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowHashesEGammaHLTPreselection  ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowHashesEGammaIDHEEP   ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowHashesEGammaIDMVAWP80   ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowHashesEGammaIDMVAWP90   ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowHashesEGammaIDMVAHZZ    ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<float> >   rhoIsoHEEP                ( new std::vector<float>   ()  ); 
  
  // Does this electron overlap with a muon?
  std::auto_ptr<std::vector<int> >     overlaps                  ( new std::vector<int>   ()  );

  // Number of Brems = number of basic clusters minus one
  std::auto_ptr<std::vector<int> >     numberOfBrems             ( new std::vector<int>   ()  );

  // Is this ECAL driven? Or PFlow?  
  std::auto_ptr<std::vector<bool> >    hasEcalDrivenSeed         ( new std::vector<bool>  ()  );
  std::auto_ptr<std::vector<bool> >    isEcalDriven              ( new std::vector<bool>  ()  );
  std::auto_ptr<std::vector<bool> >    hasTrackerDrivenSeed      ( new std::vector<bool>  ()  );
  
  // Charge consistency variables: Ferdinando Giordano
  std::auto_ptr<std::vector<bool> >    gsfCtfScPixCharge         ( new std::vector<bool>  ()  );
  std::auto_ptr<std::vector<bool> >    gsfScPixCharge            ( new std::vector<bool>  ()  );
  std::auto_ptr<std::vector<bool> >    gsfCtfCharge              ( new std::vector<bool>  ()  );
  
  // EB or EE: Ferdinando Giordano
  std::auto_ptr<std::vector<bool> >    isEB                      ( new std::vector<bool>  ()  );
  std::auto_ptr<std::vector<bool> >    isEE                      ( new std::vector<bool>  ()  );
  
  // ECAL eta/phi vs tracker eta/phi

  std::auto_ptr<std::vector<float> >  deltaPhiTrkSC             ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  deltaEtaTrkSC             ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  deltaEtaTrkSeedSC         ( new std::vector<float>()  );

  // Shower shape

  std::auto_ptr<std::vector<float> >  sigmaEtaEta               ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  sigmaIEtaIEta             ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  full5x5SigmaIEtaIEta      ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  r9                        ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  e1x5overe5x5              ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  e2x5overe5x5              ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  full5x5e1x5overe5x5       ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  full5x5e2x5overe5x5       ( new std::vector<float>()  );
  std::auto_ptr<std::vector<int> >     classif                   ( new std::vector<int>   ()  );

  // Isolation variables: PAT
  
  std::auto_ptr<std::vector<float> >  trkIsoPAT                 ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  ecalIsoPAT                ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  hcalIsoPAT                ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  relIsoPAT                 ( new std::vector<float>()  );
  std::auto_ptr<std::vector<int> >     passIsoPAT                ( new std::vector<int>   ()  );

  // Isolation variables: particle flow 
  
  std::auto_ptr<std::vector<float> >  pfChargedHadronIso03      ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  pfNeutralHadronIso03      ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  pfPhotonIso03             ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  pfPUIso03                 ( new std::vector<float>()  );

  std::auto_ptr<std::vector<float> >  pfChargedHadronIso04      ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  pfNeutralHadronIso04      ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  pfPhotonIso04             ( new std::vector<float>()  );
  
  // Isolation variables: DR 0.3

  std::auto_ptr<std::vector<float> >  ecalIsoDR03               ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  hcalIsoDR03               ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  hcalIsoDR03FullCone       ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  hcalIsoD1DR03             ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  hcalIsoD2DR03             ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  trkIsoDR03                ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  ecalPFClusterIso          ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  hcalPFClusterIso          ( new std::vector<float>()  );

  // Isolation HEEP v7.0
  std::auto_ptr<std::vector<float> >  heep70TrkIso              ( new std::vector<float>()  );

  // Conversion variables

  std::auto_ptr<std::vector<int> >     missingHits               ( new std::vector<int>   ()  );
  std::auto_ptr<std::vector<int> >     missingHitsEG             ( new std::vector<int>   ()  );
  std::auto_ptr<std::vector<float> >  dist_vec                  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  dCotTheta                 ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  fbrem                     ( new std::vector<float>()  );
  std::auto_ptr<std::vector<bool> >    hasMatchedConvPhot        ( new std::vector<bool>  ()  );
  
  // Vertex and beamspot information
  
  std::auto_ptr<std::vector<int> >     vtxIndex                  ( new std::vector<int>   ()  );
  std::auto_ptr<std::vector<float> >  vtxDistXY                 ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  vtxDistZ                  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  vtx0DistXY                ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  vtx0DistZ                 ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  primaryVertexDXY          ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  primaryVertexDXYError     ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  beamspotDXY               ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  beamspotDXYError          ( new std::vector<float>()  );

  // Track information 

  std::auto_ptr<std::vector<float> >  trackPx                   ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  trackPy                   ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  trackPz                   ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  trackPt                   ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  trackValidFractionOfHits  ( new std::vector<float>()  );  
  std::auto_ptr<std::vector<float> >  normalizedChi2            ( new std::vector<float>()  );  

  // Trigger matching: float electron

  std::auto_ptr<std::vector<bool  > >  HLTfloatEleMatched       ( new std::vector<bool  >()  );
  std::auto_ptr<std::vector<float> >  HLTfloatEleMatchPt 	 ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  HLTfloatEleMatchEta	 ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  HLTfloatEleMatchPhi      ( new std::vector<float>()  );

  // Trigger matching: Single electron

  std::auto_ptr<std::vector<bool  > >  HLTSingleEleMatched       ( new std::vector<bool  >()  );
  std::auto_ptr<std::vector<float> >  HLTSingleEleMatchPt 	 ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  HLTSingleEleMatchEta	 ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  HLTSingleEleMatchPhi      ( new std::vector<float>()  );

  // Trigger matching: Single electron (WP85)

  std::auto_ptr<std::vector<bool  > >  HLTSingleEleWP85Matched   ( new std::vector<bool  >()  );
  std::auto_ptr<std::vector<float> >  HLTSingleEleWP85MatchPt 	 ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  HLTSingleEleWP85MatchEta	 ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  HLTSingleEleWP85MatchPhi  ( new std::vector<float>()  );

  // Trigger matching: Ele+Jet+Jet

  std::auto_ptr<std::vector<bool  > >  HLTEleJetJetMatched   ( new std::vector<bool  >()  );
  std::auto_ptr<std::vector<float> >  HLTEleJetJetMatchPt 	 ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  HLTEleJetJetMatchEta	 ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  HLTEleJetJetMatchPhi  ( new std::vector<float>()  );

  // Gen matching: Status 3 only
  
  std::auto_ptr<std::vector<float> >  matchedGenParticlePt  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  matchedGenParticleEta ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  matchedGenParticlePhi ( new std::vector<float>()   );
  
  // for 03Feb2017 re-miniaod
  std::auto_ptr<bool> eGammaGSFixedDupECALClusters ( new bool() );
  std::auto_ptr<std::vector<uint32_t> > ecalMultiAndGSGlobalRecHitEBHitsNotReplaced ( new std::vector<uint32_t>() );

  //------------------------------------------------------------------------
  // Get handles for the event
  //------------------------------------------------------------------------
  
  // Tracks -- SIC note: no tracks in MiniAOD

  // track builder
  edm::ESHandle<TransientTrackBuilder> trackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);

  // DCS information -- no longer needed

  // Primary vertices
  
  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByToken(vtxInputToken_,primaryVertices);

  // Beamspot information
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamSpotInputToken_, beamSpotHandle);
  bool beamSpotValid = false;

  // Photon conversion information -- no longer needed; take what's embedded in pat::Electron

  // PAT electrons

  edm::Handle<std::vector<pat::Electron> > electrons;
  iEvent.getByToken(electronInputToken_, electrons);

  // PAT trigger matches by HLT path
  // we embed these in the regular pat::Electron collection
  // but we could make separate collections containing only the matches if we wanted
  //edm::Handle<std::vector<pat::Electron> > electronsSingleElectronHLTMatched;
  //iEvent.getByLabel(inputTag, electronsSingleElectronHLTMatched);
  //edm::Handle<std::vector<pat::Electron> > electronsSingleElectronWP85HLTMatched;
  //iEvent.getByLabel(inputTag, electronsSingleElectronWP85HLTMatched);
  //edm::Handle<std::vector<pat::Electron> > electronsDoubleElectronHLTMatched;
  //iEvent.getByLabel(inputTag, electronsDoubleElectronHLTMatched);


  // rho for EGamma isolation calculation
  edm::Handle<double> rho;
  iEvent.getByToken(rhoInputToken_, rho);
  float rhoIso = *(rho.product());
  
  // SIC add for new egamma VID framework
  // example: https://github.com/ikrav/ElectronWork/blob/master/ElectronNtupler/plugins/ElectronNtuplerIdDemoPrePHYS14miniAOD.cc
  // get electron ID maps
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  edm::Handle<edm::ValueMap<bool> > hlt_preselection_decisions;
  edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
  edm::Handle<edm::ValueMap<bool> > mva_id_wp80_decisions;
  edm::Handle<edm::ValueMap<bool> > mva_id_wp90_decisions;
  edm::Handle<edm::ValueMap<bool> > mva_id_hzz_decisions;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > veto_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > loose_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > medium_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > tight_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > hlt_preselection_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > heep_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > mva_id_wp80_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > mva_id_wp90_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > mva_id_hzz_cutflow_data;
  edm::Handle<edm::ValueMap<float> > heep70trkIsolMapHandle;
  iEvent.getByToken(electronVetoIdMapToken_,veto_id_decisions);
  iEvent.getByToken(electronLooseIdMapToken_,loose_id_decisions);
  iEvent.getByToken(electronMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(electronTightIdMapToken_,tight_id_decisions);
  iEvent.getByToken(electronHLTPreselectionMapToken_,hlt_preselection_decisions);
  iEvent.getByToken(electronHEEPIdMapToken_,heep_id_decisions);
  iEvent.getByToken(electronMVAIdWP80MapToken_,mva_id_wp80_decisions);
  iEvent.getByToken(electronMVAIdWP90MapToken_,mva_id_wp90_decisions);
  iEvent.getByToken(electronMVAIdHZZMapToken_,mva_id_hzz_decisions);
  iEvent.getByToken(eleVetoIdCutFlowResultMapToken_,veto_id_cutflow_data);
  iEvent.getByToken(eleLooseIdCutFlowResultMapToken_,loose_id_cutflow_data);
  iEvent.getByToken(eleMediumIdCutFlowResultMapToken_,medium_id_cutflow_data);
  iEvent.getByToken(eleTightIdCutFlowResultMapToken_,tight_id_cutflow_data);
  iEvent.getByToken(eleHLTPreselectionCutFlowResultMapToken_,hlt_preselection_cutflow_data);
  iEvent.getByToken(eleHEEPIdCutFlowResultMapToken_,heep_id_cutflow_data);
  iEvent.getByToken(eleMVAIdWP80CutFlowResultMapToken_,mva_id_wp80_cutflow_data);
  iEvent.getByToken(eleMVAIdWP90CutFlowResultMapToken_,mva_id_wp90_cutflow_data);
  iEvent.getByToken(eleMVAIdHZZCutFlowResultMapToken_,mva_id_hzz_cutflow_data);
  iEvent.getByToken(heep70trkIsolMapToken_,heep70trkIsolMapHandle);

  // for 03Feb2017 re-miniaod
  edm::Handle<bool> dupECALClustersHandle;
  iEvent.getByToken(dupEcalClustersToken_, dupECALClustersHandle);
  if(dupECALClustersHandle.isValid())
    *eGammaGSFixedDupECALClusters.get() = *(dupECALClustersHandle.product());
  //
  edm::Handle<edm::EDCollection<DetId> > ecalMultiAndGSGlobalRecHitEBHitsNotReplacedHandle;
  iEvent.getByToken(ecalMultiAndGSGlobalRecHitEBHitsNotReplacedToken_, ecalMultiAndGSGlobalRecHitEBHitsNotReplacedHandle);
  if(ecalMultiAndGSGlobalRecHitEBHitsNotReplacedHandle.isValid()) {
    for( auto it = ecalMultiAndGSGlobalRecHitEBHitsNotReplacedHandle->begin();
        it != ecalMultiAndGSGlobalRecHitEBHitsNotReplacedHandle->end(); ++it)
      ecalMultiAndGSGlobalRecHitEBHitsNotReplaced->push_back(*it);
  }


  bool ebRecHitsValid = false;
  edm::Handle<EcalRecHitCollection> ebRecHitsHandle;
  iEvent.getByToken(ebReducedRecHitsInputToken_, ebRecHitsHandle);
  if(ebRecHitsHandle.isValid()) {
    ebRecHitsValid = true;
  } else {
    edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get EBReducedRecHits";
  }
  bool eeRecHitsValid = false;
  edm::Handle<EcalRecHitCollection> eeRecHitsHandle;
  iEvent.getByToken(eeReducedRecHitsInputToken_, eeRecHitsHandle);
  if(eeRecHitsHandle.isValid()) {
    eeRecHitsValid = true;
  } else {
    edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get EEReducedRecHits";
  }


  //------------------------------------------------------------------------
  // Loop over electrons (finally!)
  //------------------------------------------------------------------------

  if(electrons.isValid()) {
    edm::LogInfo("RootTupleMakerV2_ElectronsInfo") << "Total # Electrons: " << electrons->size();

    size_t iElectron = 0;

    for( std::vector<pat::Electron>::const_iterator it = electrons->begin(); it != electrons->end(); ++it ) {

      //------------------------------------------------------------------------
      // Break from the loop once we have enough electrons
      //------------------------------------------------------------------------

      if(eta->size() >= maxSize) break;

      //------------------------------------------------------------------------
      // Do any of these electrons overlap with muons?
      //------------------------------------------------------------------------
      
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
      
      //------------------------------------------------------------------------
      // Look at bits for electron ID (whatever is stored in PAT)
      // passID for different electron IDs is assigned bitwise
      //  - bit 0: eidRobustLoose
      //  - bit 1: eidRobustTight
      //  - bit 2: eidLoose
      //  - bit 3: eidTight
      //  - bit 4: eidRobustHighEnergy
      //  - bit 5: MVA "trig" 
      //  - bit 6: MVA "non-trig"
      // SIC XXX FIXME: probably need to update this using the below, taken from
      //    https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/DataFormats/PatCandidates/interface/Electron.h
       /// Returns a specific electron ID associated to the pat::Electron given its name
      // For cut-based IDs, the value map has the following meaning:
      // 0: fails,
      // 1: passes electron ID only,
      // 2: passes electron Isolation only,
      // 3: passes electron ID and Isolation only,
      // 4: passes conversion rejection,
      // 5: passes conversion rejection and ID,
      // 6: passes conversion rejection and Isolation,
      // 7: passes the whole selection.
      // For more details have a look at:
      // https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID
      // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCategoryBasedElectronID
      // Note: an exception is thrown if the specified ID is not available
      //------------------------------------------------------------------------
      // In Phys14 samples, these are CSA14-tuned ID's such as: heepElectronID-HEEPV50-CSA14-25ns
      //  and the below
      
      int passId = 0;
      if (it->electronID("eidRobustLoose"     )>0) passId = passId | 1<<0;
      if (it->electronID("eidRobustTight"     )>0) passId = passId | 1<<1;
      if (it->electronID("eidLoose"           )>0) passId = passId | 1<<2;
      if (it->electronID("eidTight"           )>0) passId = passId | 1<<3;
      if (it->electronID("eidRobustHighEnergy")>0) passId = passId | 1<<4;
      //if (it->electronID("mvaTrigV0"          )>0) passId = passId | 1<<5;
      //if (it->electronID("mvaNonTrigV0"       )>0) passId = passId | 1<<6;
      
      //------------------------------------------------------------------------
      // Trigger matching
      //------------------------------------------------------------------------
      // SIC: we've embedded matches to selected HLT paths in the python with the PATTriggerMatchEmbedder.
      //      now just ask if we have a match to whichever HLT path in the object

      //TEST
      //std::cout << "size of trigger matches: " << it->triggerObjectMatches().size() << std::endl;
      //TEST

      //// Double electron
      //const pat::TriggerObjectStandAloneCollection matchesDoubleEle = it->triggerObjectMatchesByPath("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*");
      //if(matchesDoubleEle.size() > 0)
      //{
      //  HLTDoubleEleMatched  -> push_back ( true ) ;
      //  HLTDoubleEleMatchPt  -> push_back ( matchesDoubleEle[0].pt() );
      //  HLTDoubleEleMatchEta -> push_back ( matchesDoubleEle[0].eta());
      //  HLTDoubleEleMatchPhi -> push_back ( matchesDoubleEle[0].phi());
      //}
      //else
      //{
      //  HLTDoubleEleMatched  -> push_back ( false ) ;
      //  HLTDoubleEleMatchPt  -> push_back ( -999. );
      //  HLTDoubleEleMatchEta -> push_back ( -999. );
      //  HLTDoubleEleMatchPhi -> push_back ( -999. );
      //}

      // Single electron
      //FIXME: TRY WP75 path in recent menu
      //const pat::TriggerObjectStandAloneCollection matchesDoubleEle = it->triggerObjectMatchesByPath("*");
      //if(matchesSingleEle.size() > 0)
      //{
      //  HLTSingleEleMatched  -> push_back ( true ) ;
      //  HLTSingleEleMatchPt  -> push_back ( matchesSingleEle[0] -> pt() );
      //  HLTSingleEleMatchEta -> push_back ( matchesSingleEle[0] -> eta());
      //  HLTSingleEleMatchPhi -> push_back ( matchesSingleEle[0] -> phi());
      //}
      //else
      //{ 
      //  HLTSingleEleMatched  -> push_back ( false ) ;
      //  HLTSingleEleMatchPt  -> push_back ( -999. );
      //  HLTSingleEleMatchEta -> push_back ( -999. );
      //  HLTSingleEleMatchPhi -> push_back ( -999. );
      //}

      //// Single electron (WP85)
      //const pat::TriggerObjectStandAloneCollection matchesSingleEleWP85 = it->triggerObjectMatchesByPath("HLT_Ele32_eta2p1_WP85_Gsf_v*");
      //if(matchesSingleEleWP85.size() > 0)
      //{
      //  HLTSingleEleWP85Matched  -> push_back ( true ) ;
      //  HLTSingleEleWP85MatchPt  -> push_back ( matchesSingleEleWP85[0].pt() );
      //  HLTSingleEleWP85MatchEta -> push_back ( matchesSingleEleWP85[0].eta());
      //  HLTSingleEleWP85MatchPhi -> push_back ( matchesSingleEleWP85[0].phi());
      //}
      //else 
      //{ 
      //  HLTSingleEleWP85Matched  -> push_back ( false ) ;
      //  HLTSingleEleWP85MatchPt  -> push_back ( -999. );
      //  HLTSingleEleWP85MatchEta -> push_back ( -999. );
      //  HLTSingleEleWP85MatchPhi -> push_back ( -999. );
      //}

      //// E+J+J cross trigger
      //const pat::TriggerObjectStandAloneCollection matchesEleJetJet = it->triggerObjectMatchesByPath("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*");
      //if(matchesEleJetJet.size() > 0)
      //{
      //  HLTEleJetJetMatched  -> push_back ( true ) ;
      //  HLTEleJetJetMatchPt  -> push_back ( matchesEleJetJet[0].pt() );
      //  HLTEleJetJetMatchEta -> push_back ( matchesEleJetJet[0].eta());
      //  HLTEleJetJetMatchPhi -> push_back ( matchesEleJetJet[0].phi());
      //}
      //else 
      //{ 
      //  HLTEleJetJetMatched  -> push_back ( false ) ;
      //  HLTEleJetJetMatchPt  -> push_back ( -999. );
      //  HLTEleJetJetMatchEta -> push_back ( -999. );
      //  HLTEleJetJetMatchPhi -> push_back ( -999. );
      //}

      //------------------------------------------------------------------------
      // Gen matching
      // This should work for pythia8 in MiniAOD. should use status 23 I think for outgoing electrons.
      // See: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#MC_Truth
      //------------------------------------------------------------------------

      float genPartPt = -999.;
      float genPartEta= -999.;
      float genPartPhi= -999.;
      
      if ( !iEvent.isRealData() ) {
        for(uint igen = 0 ; igen < it->genParticleRefs().size() ; ++igen ){ //it->genParticleRefs().size() should be 0, 1 or 2                
          if(it->genParticleRef(igen).isNonnull()) {
            if(  it->genParticle(igen)->status()==1 || it->genParticle(igen)->status()==3 || it->genParticle(igen)->status()==23){
              genPartPt =it->genParticle(igen)->pt();
              genPartEta=it->genParticle(igen)->eta();
              genPartPhi=it->genParticle(igen)->phi();
            }
          }
          else
            edm::LogError("RootTupleMakerV2_ElectronsError") << "genParticleRef " << igen+1 << "/" << it->genParticleRefs().size() << " is null!";
        }
        //XXX FIXME TODO: can't we just use this?
        //genPartPt =it->genParticle()->pt();
        //genPartEta=it->genParticle()->eta();
        //genPartPhi=it->genParticle()->phi();
      }
      
      matchedGenParticlePt  -> push_back ( (float)(genPartPt ) );
      matchedGenParticleEta -> push_back ( (float)(genPartEta) );
      matchedGenParticlePhi -> push_back ( (float)(genPartPhi) );
      
      //------------------------------------------------------------------------
      // Relative isolation (not currently used in any analysis... remove?) 
      //------------------------------------------------------------------------
      
      float reliso = (it->trackIso() + it->ecalIso() + it->hcalIso())/it->pt();

      //------------------------------------------------------------------------
      // Conversion information
      //------------------------------------------------------------------------
      // Now take what's embedded in the pat candidate

      //------------------------------------------------------------------------
      // Vertex association
      //------------------------------------------------------------------------

      float minVtxDist3D = 9999.;
      int vtxIndex_ = -1;
      float vtxDistXY_ = -9999.;
      float vtxDistZ_ = -9999.;
      
      float vtx0DistXY_;
      float vtx0DistZ_;

      if(primaryVertices.isValid()) {
        edm::LogInfo("RootTupleMakerV2_ElectronsInfo") << "Total # Primary Vertices: " << primaryVertices->size();

        int i_vertex = 0;
        for( reco::VertexCollection::const_iterator v_it=primaryVertices->begin() ; v_it!=primaryVertices->end() ; ++v_it ) {

          float distXY = it->gsfTrack()->dxy(v_it->position());
          float distZ = it->gsfTrack()->dz(v_it->position());
          float dist3D = sqrt(pow(distXY,2) + pow(distZ,2));

          if ( i_vertex == 0 ) { 
            vtx0DistXY_ = distXY;
            vtx0DistZ_  = distZ ;
          }

          if( dist3D<minVtxDist3D ) {
            minVtxDist3D = dist3D;
            vtxIndex_ = int(std::distance(primaryVertices->begin(),v_it));
            vtxDistXY_ = distXY;
            vtxDistZ_ = distZ;
          }

          i_vertex++;
        }
      } else {
        edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get the vertex collection";
      }

      if(beamSpotHandle.isValid()) {
        beamSpotValid = true;
      } else {
        edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get the beamspot";
      }

      //------------------------------------------------------------------------ 
      // Assign values to variables
      //------------------------------------------------------------------------
      
      // Kinematic variables
      
      eta                      -> push_back ( it->eta() );
      phi                      -> push_back ( it->phi() );
      pt                       -> push_back ( it->pt() );
      ptHeep                   -> push_back ( it->caloEnergy()*sin(it->p4().theta()) );
      energy                   -> push_back ( it->energy() );
      caloEnergy               -> push_back ( it->caloEnergy() );
      ecalEnergy               -> push_back ( it->ecalEnergy() );
      charge                   -> push_back ( it->charge() );
      hoe                      -> push_back ( it->hadronicOverEm() );
      
      // Supercluster kinematic variables

      scEta                    -> push_back( it->superCluster()->eta() );
      scSeedEta                -> push_back( it->superCluster()->seed()->eta() );
      scPhi                    -> push_back( it->superCluster()->phi() );
      scPt                     -> push_back( it->superCluster()->energy()/cosh(it->superCluster()->eta()) );
      scRawEnergy              -> push_back( it->superCluster()->rawEnergy() );
      scEnergy                 -> push_back( it->superCluster()->energy() );
      eSuperClusterOverP       -> push_back( it->eSuperClusterOverP() );
      // find seed cry of supercluster and get its energy out of the reduced rechit collection
      DetId seedDetId = it->superCluster()->seed()->seed();
      if(seedDetId.subdetId() == EcalBarrel) {
        scSeedCryIsBarrel      -> push_back(true);
        if(!ebRecHitsValid)
          scSeedCryEnergy      -> push_back( -999 );
        else {
          auto rhItr =  ebRecHitsHandle->find(seedDetId);
          if(rhItr != ebRecHitsHandle->end())
            scSeedCryEnergy    -> push_back( rhItr->energy());
          else
            scSeedCryEnergy    -> push_back( -999 );
        }
      }
      else if(seedDetId.subdetId() == EcalEndcap) {
        scSeedCryIsBarrel      -> push_back(false);
        if(!eeRecHitsValid)
          scSeedCryEnergy      -> push_back( -999 );
        else {
          auto rhItr =  eeRecHitsHandle->find(seedDetId);
          if(rhItr != eeRecHitsHandle->end())
            scSeedCryEnergy    -> push_back( rhItr->energy());
          else
            scSeedCryEnergy    -> push_back( -999 );
        }
      }
      ////
      //std::cout << "Seed cry detId: " << seedDetId.rawId() << std::endl;
      //for(auto rhItr2 = it->recHits()->begin(); rhItr2 != it->recHits()->end(); ++rhItr2 ) {
      //  std::cout << "RecHit in Electron -- DetId: " << rhItr2->detid().rawId() << std::endl;
      //}
      ////

      // ID information
      passIds                  -> push_back( passId );
      // SIC: update to Run2 cut-based and HEEP ID's
      const edm::Ptr<pat::Electron> elPtr(electrons, it - electrons->begin() );
      passEGammaIDVeto           -> push_back ((*veto_id_decisions)[ elPtr ]);
      passEGammaIDLoose          -> push_back ((*loose_id_decisions)[ elPtr ]);
      passEGammaIDMedium         -> push_back ((*medium_id_decisions)[ elPtr ]);
      passEGammaIDTight          -> push_back ((*tight_id_decisions)[ elPtr ]);
      passEGammaHLTPreselection  -> push_back ((*hlt_preselection_decisions)[ elPtr ]);
      passHEEPID                 -> push_back ((*heep_id_decisions)[ elPtr ]);
      passMVAIDWP80              -> push_back ((*mva_id_wp80_decisions)[ elPtr ]);
      passMVAIDWP90              -> push_back ((*mva_id_wp90_decisions)[ elPtr ]);
      passMVAIDHZZ               -> push_back ((*mva_id_hzz_decisions)[ elPtr ]);
      cutFlowNamesEGammaIDVeto   -> push_back (((*veto_id_cutflow_data)[ elPtr ]).cutFlowName());
      cutFlowNamesEGammaIDLoose  -> push_back (((*loose_id_cutflow_data)[ elPtr ]).cutFlowName());
      cutFlowNamesEGammaIDMedium -> push_back (((*medium_id_cutflow_data)[ elPtr ]).cutFlowName());
      cutFlowNamesEGammaIDTight  -> push_back (((*tight_id_cutflow_data)[ elPtr ]).cutFlowName());
      cutFlowNamesEGammaHLTPreselection  -> push_back (((*hlt_preselection_cutflow_data)[ elPtr ]).cutFlowName());
      cutFlowNamesEGammaIDHEEP   -> push_back (((*heep_id_cutflow_data)[ elPtr ]).cutFlowName());
      cutFlowNamesEGammaIDMVAWP80-> push_back (((*mva_id_wp80_cutflow_data)[ elPtr ]).cutFlowName());
      cutFlowNamesEGammaIDMVAWP90-> push_back (((*mva_id_wp90_cutflow_data)[ elPtr ]).cutFlowName());
      cutFlowNamesEGammaIDMVAHZZ-> push_back (((*mva_id_hzz_cutflow_data)[ elPtr ]).cutFlowName());
      cutFlowHashesEGammaIDVeto  -> push_back (((*veto_id_cutflow_data)[ elPtr ]).cutFlowHash());
      cutFlowHashesEGammaIDLoose -> push_back (((*loose_id_cutflow_data)[ elPtr ]).cutFlowHash());
      cutFlowHashesEGammaIDMedium-> push_back (((*medium_id_cutflow_data)[ elPtr ]).cutFlowHash());
      cutFlowHashesEGammaIDTight -> push_back (((*tight_id_cutflow_data)[ elPtr ]).cutFlowHash());
      cutFlowHashesEGammaHLTPreselection -> push_back (((*hlt_preselection_cutflow_data)[ elPtr ]).cutFlowHash());
      cutFlowHashesEGammaIDHEEP  -> push_back (((*heep_id_cutflow_data)[ elPtr ]).cutFlowHash());
      cutFlowHashesEGammaIDMVAWP80  -> push_back (((*mva_id_wp80_cutflow_data)[ elPtr ]).cutFlowHash());
      cutFlowHashesEGammaIDMVAWP90  -> push_back (((*mva_id_wp90_cutflow_data)[ elPtr ]).cutFlowHash());
      cutFlowHashesEGammaIDMVAHZZ   -> push_back (((*mva_id_hzz_cutflow_data)[ elPtr ]).cutFlowHash());
      //
      rhoIsoHEEP               -> push_back (rhoIso);

      // Does this electron overlap with a muon?
      overlaps                 -> push_back( ovrlps );
      
      // Number of Brems = number of basic clusters minus one
      numberOfBrems            -> push_back( it->numberOfBrems() );
      
      // Is this ECAL driven? Or PFlow?  
      hasEcalDrivenSeed        -> push_back( it->ecalDrivenSeed() );
      isEcalDriven             -> push_back( it->ecalDriven() );
      hasTrackerDrivenSeed     -> push_back( it->trackerDrivenSeed() );

      // Charge consistency variables - Ferdinando Giordano
      gsfCtfScPixCharge        -> push_back( it->isGsfCtfScPixChargeConsistent() );
      gsfScPixCharge           -> push_back( it->isGsfScPixChargeConsistent() );
      gsfCtfCharge             -> push_back( it->isGsfCtfChargeConsistent() );

      // EB or EE - Ferdinando Giordano
      isEB                     -> push_back( it->isEB() );
      isEE                     -> push_back( it->isEE() );
      
      // ECAL eta/phi vs tracker eta/phi
      deltaPhiTrkSC            -> push_back ( it->deltaPhiSuperClusterTrackAtVtx() );
      deltaEtaTrkSC            -> push_back ( it->deltaEtaSuperClusterTrackAtVtx() );
      deltaEtaTrkSeedSC        -> push_back ( it->deltaEtaSeedClusterTrackAtVtx() );

      // Shower shape
      sigmaEtaEta              -> push_back ( it->sigmaEtaEta() );
      sigmaIEtaIEta            -> push_back ( it->sigmaIetaIeta() );
      full5x5SigmaIEtaIEta     -> push_back ( it->full5x5_sigmaIetaIeta() );
      classif                  -> push_back ( it->classification() );
      r9                       -> push_back ( it->r9() );
      e1x5overe5x5             -> push_back ( (it->e5x5()>0) ? (it->e1x5()/it->e5x5()) : 0 );
      e2x5overe5x5             -> push_back ( (it->e5x5()>0) ? (it->e2x5Max()/it->e5x5()) : 0 );
      full5x5e1x5overe5x5      -> push_back ( (it->full5x5_e5x5()>0) ? (it->full5x5_e1x5()/it->full5x5_e5x5()) : 0 );
      full5x5e2x5overe5x5      -> push_back ( (it->full5x5_e5x5()>0) ? (it->full5x5_e2x5Max()/it->full5x5_e5x5()) : 0 );

      // Isolation variables: PAT
      // dR 0.4, detector isolation
      trkIsoPAT                -> push_back( it->trackIso() ); // note: same as userIsolation(pat::TrackIso)
      ecalIsoPAT               -> push_back( it->ecalIso() );
      hcalIsoPAT               -> push_back( it->hcalIso() );
      relIsoPAT                -> push_back( reliso );
      passIsoPAT               -> push_back( (reliso<electronIso) ? 1 : 0 ); 

      // Isolation variables: dR 0.3, detector isolation
      // from reco::GsfElectron methods
      ecalIsoDR03              -> push_back ( it->dr03EcalRecHitSumEt() );
      hcalIsoDR03              -> push_back ( it->dr03HcalTowerSumEt() );
      hcalIsoDR03FullCone      -> push_back ( it->dr03HcalTowerSumEt() +
					      ( it->hadronicOverEm() 
						* it->superCluster()->energy() 
						/ cosh(it->superCluster()->eta())));
      hcalIsoD1DR03            -> push_back( it->dr03HcalDepth1TowerSumEt() );
      hcalIsoD2DR03            -> push_back( it->dr03HcalDepth2TowerSumEt() );
      trkIsoDR03               -> push_back( it->dr03TkSumPt() );

      // pf clusters iso
      ecalPFClusterIso         -> push_back( it->ecalPFClusterIso() );
      hcalPFClusterIso         -> push_back( it->hcalPFClusterIso() );

      // Isolation HEEP v7.0
      if(!heep70trkIsolMapHandle.isValid())
        edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get the HEEP7.0 track isolation value map";
      else
        heep70TrkIso        -> push_back((*heep70trkIsolMapHandle)[ elPtr ]);


      // Isolation variables: particle flow
      // methods from reco::GsfElectron; should be filled with dR 0.3 values from isolation value maps, e.g., gedElPFIsoValueCharged03
      // See http://cmslxr.fnal.gov/source/RecoEgamma/EgammaElectronProducers/python/gedGsfElectronFinalizer_cfi.py?v=CMSSW_7_2_3
      pfChargedHadronIso03     -> push_back ( it->pfIsolationVariables().sumChargedHadronPt );
      pfPhotonIso03            -> push_back ( it->pfIsolationVariables().sumPhotonEt );
      pfNeutralHadronIso03     -> push_back ( it->pfIsolationVariables().sumNeutralHadronEt );
      pfPUIso03                -> push_back ( it->pfIsolationVariables().sumPUPt );
      // chargedHadronIso() is same as userIsolation(pat::PfChargedHadronIso); // dR = 0.4, filled in PAT electron producer
      // see: http://cmslxr.fnal.gov/source/PhysicsTools/PatAlgos/python/producersLayer1/electronProducer_cff.py?v=CMSSW_7_2_3
      pfChargedHadronIso04     -> push_back ( it->chargedHadronIso() );
      pfPhotonIso04            -> push_back ( it->photonIso() );
      pfNeutralHadronIso04     -> push_back ( it->neutralHadronIso() );
      // Above dR 0.3/0.4 storage is confirmed by http://cmslxr.fnal.gov/source/PhysicsTools/Heppy/python/physicsobjects/Electron.py?v=CMSSW_7_3_1
      // See line 129 for chargedHadronIso

      
      // Conversion variables
      constexpr reco::HitPattern::HitCategory missingHitType = reco::HitPattern::MISSING_INNER_HITS;
      //missingHits              -> push_back ( it->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() );
      //missingHits              -> push_back ( it->gsfTrack()->numberOfLostHits() );
      //missingHitsEG            -> push_back ( it->gsfTrack()->trackerExpectedHitsInner().numberOfHits()     );
      missingHits            -> push_back ( it->gsfTrack()->hitPattern().numberOfHits(missingHitType) );
      
      dist_vec                 -> push_back ( it->convDist() ); // from reco::GsfElectron
      dCotTheta                -> push_back ( it->convDcot() ); // from reco::GsfElectron
      hasMatchedConvPhot       -> push_back ( !(it->passConversionVeto()) );
      fbrem                    -> push_back ( it->fbrem() );
      //REMOVED convFitProb_vec          -> push_back ( convFitProb );

      // Vertex association variables
      
      vtxIndex                 -> push_back( vtxIndex_  );
      vtxDistXY                -> push_back( vtxDistXY_ );
      vtxDistZ                 -> push_back( vtxDistZ_  );
      vtx0DistXY               -> push_back( vtx0DistXY_);
      vtx0DistZ                -> push_back( vtx0DistZ_ );
      primaryVertexDXY         -> push_back( fabs( it->dB() ) );      
      primaryVertexDXYError    -> push_back( fabs( it->edB() ) );
      if(beamSpotValid)
      {
        beamspotDXY              -> push_back( fabs( it->dB (pat::Electron::BS2D) ) );
        // recalculate the beamspotDXYError, since beamspot is always marked as invalid in PatElectronProducer
        // code taken from PatElectronProducer
        reco::TransientTrack tt = trackBuilder->build(it->gsfTrack());
        reco::Vertex vBeamspot(beamSpotHandle->position(),beamSpotHandle->covariance3D());
        std::pair<bool,Measurement1D> result =
          IPTools::signedTransverseImpactParameter(tt,GlobalVector(it->gsfTrack()->px(),it->gsfTrack()->py(),it->gsfTrack()->pz()),
              vBeamspot);
        beamspotDXYError         -> push_back( fabs(result.second.error()) );
      }
      else
      {
        beamspotDXY              -> push_back( -999 );
        beamspotDXYError         -> push_back( -999 );
      }
			       
      // Track information     
      trackPx                  -> push_back( it->gsfTrack()->px() );
      trackPy                  -> push_back( it->gsfTrack()->py() );
      trackPz                  -> push_back( it->gsfTrack()->pz() );
      trackPt                  -> push_back( it->gsfTrack()->pt() );
      trackValidFractionOfHits -> push_back( it->gsfTrack()->validFraction() );
      // track chi2
      normalizedChi2           ->push_back( it->gsfTrack()->normalizedChi2() );

      ++iElectron;
    }
  } else {
    edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get the electrons";
  }

  //------------------------------------------------------------------------
  // Put vectors in the event
  //------------------------------------------------------------------------

  // Kinematic variables

  iEvent.put( eta                     , prefix + "Eta"                      + suffix );
  iEvent.put( phi                     , prefix + "Phi"                      + suffix );
  iEvent.put( pt                      , prefix + "Pt"                       + suffix );
  iEvent.put( ptHeep                  , prefix + "PtHeep"                   + suffix );
  iEvent.put( energy                  , prefix + "Energy"                   + suffix );
  iEvent.put( caloEnergy              , prefix + "CaloEnergy"               + suffix );
  iEvent.put( ecalEnergy              , prefix + "EcalEnergy"               + suffix );
  iEvent.put( charge                  , prefix + "Charge"                   + suffix );
  iEvent.put( hoe                     , prefix + "HoE"                      + suffix );

  // Supercluster kinematic variables				        

  iEvent.put( eSuperClusterOverP      , prefix + "ESuperClusterOverP"       + suffix );
  iEvent.put( scEta                   , prefix + "SCEta"                    + suffix );
  iEvent.put( scSeedEta               , prefix + "SCSeedEta"                + suffix );
  iEvent.put( scPhi                   , prefix + "SCPhi"                    + suffix );
  iEvent.put( scPt                    , prefix + "SCPt"                     + suffix );
  iEvent.put( scRawEnergy             , prefix + "SCRawEnergy"              + suffix );
  iEvent.put( scEnergy                , prefix + "SCEnergy"                 + suffix );
  iEvent.put( scSeedCryEnergy         , prefix + "SCSeedCryEnergy"          + suffix );
  iEvent.put( scSeedCryIsBarrel       , prefix + "SCSeedCryIsBarrel"        + suffix );

  // ID information 
  iEvent.put( passIds                 , prefix + "PassId"                   + suffix );
  iEvent.put( passEGammaIDVeto        , prefix + "PassEGammaIDVeto"         + suffix );
  iEvent.put( passEGammaIDLoose       , prefix + "PassEGammaIDLoose"        + suffix );
  iEvent.put( passEGammaIDMedium      , prefix + "PassEGammaIDMedium"       + suffix );
  iEvent.put( passEGammaIDTight       , prefix + "PassEGammaIDTight"        + suffix );
  iEvent.put( passEGammaHLTPreselection       , prefix + "PassEGammaHLTPreselection"        + suffix );
  iEvent.put( passHEEPID              , prefix + "PassHEEPID"               + suffix );
  iEvent.put( passMVAIDWP80           , prefix + "PassMVAIDWP80"               + suffix );
  iEvent.put( passMVAIDWP90           , prefix + "PassMVAIDWP90"               + suffix );
  iEvent.put( passMVAIDHZZ            , prefix + "PassMVAIDHZZ"               + suffix );
  //iEvent.put( cutFlowNamesEGammaIDVeto   , prefix + "CutFlowNamesEGammaIDVeto"   + suffix );
  //iEvent.put( cutFlowNamesEGammaIDLoose  , prefix + "CutFlowNamesEGammaIDLoose"  + suffix );
  //iEvent.put( cutFlowNamesEGammaIDMedium , prefix + "CutFlowNamesEGammaIDMedium" + suffix );
  //iEvent.put( cutFlowNamesEGammaIDTight  , prefix + "CutFlowNamesEGammaIDTight"  + suffix );
  iEvent.put( cutFlowNamesEGammaIDHEEP   , prefix + "CutFlowNamesEGammaIDHEEP"   + suffix );
  //iEvent.put( cutFlowHashesEGammaIDVeto  , prefix + "CutFlowHashesEGammaIDVeto"  + suffix );
  //iEvent.put( cutFlowHashesEGammaIDLoose , prefix + "CutFlowHashesEGammaIDLoose" + suffix );
  //iEvent.put( cutFlowHashesEGammaIDMedium, prefix + "CutFlowHashesEGammaIDMedium"+ suffix );
  //iEvent.put( cutFlowHashesEGammaIDTight , prefix + "CutFlowHashesEGammaIDTight" + suffix );
  iEvent.put( cutFlowHashesEGammaIDHEEP  , prefix + "CutFlowHashesEGammaIDHEEP"  + suffix );
  iEvent.put( rhoIsoHEEP              , prefix + "RhoIsoHEEP"               + suffix );
  
  // Does this electron overlap with a muon?			        
  iEvent.put( overlaps                , prefix + "Overlaps"                 + suffix );

  // Number of Brems = number of basic clusters minus one	        
  iEvent.put( numberOfBrems           , prefix + "NumberOfBrems"            + suffix );
  
  // Is this ECAL driven? Or PFlow?  				        
  iEvent.put( hasEcalDrivenSeed       , prefix + "HasEcalDrivenSeed"        + suffix );
  iEvent.put( isEcalDriven            , prefix + "IsEcalDriven"             + suffix );
  iEvent.put( hasTrackerDrivenSeed    , prefix + "HasTrackerDrivenSeed"     + suffix );

  // Charge consistency variables - Ferdinando Giordano
  iEvent.put( gsfCtfScPixCharge       , prefix + "GsfCtfScPixCharge"        + suffix );
  iEvent.put( gsfScPixCharge          , prefix + "GsfScPixCharge"           + suffix );
  iEvent.put( gsfCtfCharge            , prefix + "GsfCtfCharge"             + suffix );

  // EB or EE - Ferdinando Giordano

  iEvent.put( isEB                    , prefix + "IsEB"                     + suffix );
  iEvent.put( isEE                    , prefix + "IsEE"                     + suffix );
  
  // ECAL eta/phi vs tracker eta/phi				        

  iEvent.put( deltaPhiTrkSC           , prefix + "DeltaPhiTrkSC"            + suffix );
  iEvent.put( deltaEtaTrkSC           , prefix + "DeltaEtaTrkSC"            + suffix );
  iEvent.put( deltaEtaTrkSeedSC       , prefix + "DeltaEtaTrkSeedSC"        + suffix );

  // Shower shape						        

  iEvent.put( sigmaEtaEta             , prefix + "SigmaEtaEta"              + suffix );
  //iEvent.put( sigmaIEtaIEta           , prefix + "SigmaIEtaIEta"            + suffix );
  iEvent.put( full5x5SigmaIEtaIEta    , prefix + "Full5x5SigmaIEtaIEta"     + suffix );
  iEvent.put( r9                      , prefix + "R9"                       + suffix );
  //iEvent.put( e1x5overe5x5            , prefix + "E1x5OverE5x5"             + suffix );
  //iEvent.put( e2x5overe5x5            , prefix + "E2x5OverE5x5"             + suffix );
  iEvent.put( full5x5e1x5overe5x5        , prefix + "Full5x5E1x5OverE5x5"         + suffix );
  iEvent.put( full5x5e2x5overe5x5        , prefix + "Full5x5E2x5OverE5x5"         + suffix );
  iEvent.put( classif                 , prefix + "Classif"                  + suffix );
  
  //// Isolation variables: PAT					        

  //iEvent.put( trkIsoPAT               , prefix + "TrkIsoPAT"                + suffix );
  //iEvent.put( ecalIsoPAT              , prefix + "EcalIsoPAT"               + suffix );
  //iEvent.put( hcalIsoPAT              , prefix + "HcalIsoPAT"               + suffix );
  //iEvent.put( relIsoPAT               , prefix + "RelIsoPAT"                + suffix );
  //iEvent.put( passIsoPAT              , prefix + "PassIsoPAT"               + suffix );

  // Isolation variables: DR 0.3				        

  iEvent.put( ecalIsoDR03             , prefix + "EcalIsoDR03"              + suffix );
  iEvent.put( hcalIsoDR03             , prefix + "HcalIsoDR03"              + suffix );
  iEvent.put( hcalIsoDR03FullCone     , prefix + "HcalIsoDR03FullCone"      + suffix );
  iEvent.put( hcalIsoD1DR03           , prefix + "HcalIsoD1DR03"            + suffix );
  iEvent.put( hcalIsoD2DR03           , prefix + "HcalIsoD2DR03"            + suffix );
  iEvent.put( trkIsoDR03              , prefix + "TrkIsoDR03"               + suffix );
  iEvent.put( ecalPFClusterIso        , prefix + "EcalPFClusterIso"         + suffix );
  iEvent.put( hcalPFClusterIso        , prefix + "HcalPFClusterIso"         + suffix );
  // Isolation HEEP 7.0
  iEvent.put ( heep70TrkIso           , prefix + "Heep70TrkIso"             + suffix );
  // Isolation variables: particle flow

  iEvent.put( pfChargedHadronIso03    , prefix + "PFChargedHadronIso03"     + suffix );
  iEvent.put( pfNeutralHadronIso03    , prefix + "PFNeutralHadronIso03"     + suffix );
  iEvent.put( pfPhotonIso03           , prefix + "PFPhotonIso03"            + suffix );
  iEvent.put( pfPUIso03               , prefix + "PFPUIso03"                + suffix );

  //iEvent.put( pfChargedHadronIso04    , prefix + "PFChargedHadronIso04"     + suffix );
  //iEvent.put( pfNeutralHadronIso04    , prefix + "PFNeutralHadronIso04"     + suffix );
  //iEvent.put( pfPhotonIso04           , prefix + "PFPhotonIso04"            + suffix );

  // Conversion variables					        
  
  iEvent.put( missingHits             , prefix + "MissingHits"              + suffix );
  iEvent.put( missingHitsEG           , prefix + "MissingHitsEG"            + suffix );
  iEvent.put( dist_vec                , prefix + "Dist"                     + suffix );
  iEvent.put( dCotTheta               , prefix + "DCotTheta"                + suffix );
  iEvent.put( fbrem                   , prefix + "Fbrem"                    + suffix );
  iEvent.put( hasMatchedConvPhot      , prefix + "HasMatchedConvPhot"       + suffix );

  // Vertex and beamspot information

  iEvent.put( vtxIndex                , prefix + "VtxIndex"                 + suffix );
  iEvent.put( vtxDistXY               , prefix + "VtxDistXY"                + suffix );
  iEvent.put( vtxDistZ                , prefix + "VtxDistZ"                 + suffix );
  iEvent.put( vtx0DistXY              , prefix + "LeadVtxDistXY"            + suffix );
  iEvent.put( vtx0DistZ               , prefix + "LeadVtxDistZ"             + suffix );
  iEvent.put( primaryVertexDXY        , prefix + "PrimaryVertexDXY"         + suffix );
  iEvent.put( primaryVertexDXYError   , prefix + "PrimaryVertexDXYError"    + suffix );
  iEvent.put( beamspotDXY             , prefix + "BeamSpotDXY"              + suffix );
  iEvent.put( beamspotDXYError        , prefix + "BeamSpotDXYError"         + suffix );

  // Track information

  iEvent.put( trackPx                 , prefix + "TrackPx"                  + suffix );
  iEvent.put( trackPy                 , prefix + "TrackPy"                  + suffix );
  iEvent.put( trackPz                 , prefix + "TrackPz"                  + suffix );
  iEvent.put( trackPt                 , prefix + "TrackPt"                  + suffix );
  iEvent.put( trackValidFractionOfHits, prefix + "TrackValidFractionOfHits" + suffix );
  iEvent.put( normalizedChi2          , prefix + "NormalizedChi2" + suffix );

  //// Trigger matching: Double electron

  //iEvent.put( HLTDoubleEleMatched     , prefix + "HLTDoubleEleMatched"      + suffix );
  //iEvent.put( HLTDoubleEleMatchPt     , prefix + "HLTDoubleEleMatchPt"      + suffix );
  //iEvent.put( HLTDoubleEleMatchEta    , prefix + "HLTDoubleEleMatchEta"     + suffix );
  //iEvent.put( HLTDoubleEleMatchPhi    , prefix + "HLTDoubleEleMatchPhi"     + suffix );

  //// Trigger matching: Single electron

  //iEvent.put( HLTSingleEleMatched     , prefix + "HLTSingleEleMatched"      + suffix );
  //iEvent.put( HLTSingleEleMatchPt     , prefix + "HLTSingleEleMatchPt"      + suffix );
  //iEvent.put( HLTSingleEleMatchEta    , prefix + "HLTSingleEleMatchEta"     + suffix );
  //iEvent.put( HLTSingleEleMatchPhi    , prefix + "HLTSingleEleMatchPhi"     + suffix );

  //// Trigger matching: Single electron (WP85)

  //iEvent.put( HLTSingleEleWP85Matched , prefix + "HLTSingleEleWP85Matched"  + suffix );
  //iEvent.put( HLTSingleEleWP85MatchPt , prefix + "HLTSingleEleWP85MatchPt"  + suffix );
  //iEvent.put( HLTSingleEleWP85MatchEta, prefix + "HLTSingleEleWP85MatchEta" + suffix );
  //iEvent.put( HLTSingleEleWP85MatchPhi, prefix + "HLTSingleEleWP85MatchPhi" + suffix );

  //// Trigger matching: Ele+Jet+Jet

  //iEvent.put( HLTEleJetJetMatched , prefix + "HLTEleJetJetMatched"  + suffix );
  //iEvent.put( HLTEleJetJetMatchPt , prefix + "HLTEleJetJetMatchPt"  + suffix );
  //iEvent.put( HLTEleJetJetMatchEta, prefix + "HLTEleJetJetMatchEta" + suffix );
  //iEvent.put( HLTEleJetJetMatchPhi, prefix + "HLTEleJetJetMatchPhi" + suffix );

  // Gen matching: Status 3 only

  iEvent.put( matchedGenParticlePt ,   prefix + "MatchedGenParticlePt"      + suffix );
  iEvent.put( matchedGenParticleEta,   prefix + "MatchedGenParticleEta"     + suffix );
  iEvent.put( matchedGenParticlePhi,   prefix + "MatchedGenParticlePhi"     + suffix );

  // for 03Feb2017 re-miniaod
  iEvent.put ( eGammaGSFixedDupECALClusters, prefix + "EGammaGSFixedDupECALClusters"     + suffix );
  iEvent.put ( ecalMultiAndGSGlobalRecHitEBHitsNotReplaced, prefix + "EcalMultiAndGSGlobalRecHitEBHitsNotReplaced"     + suffix );

}
