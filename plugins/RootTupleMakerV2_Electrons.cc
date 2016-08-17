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

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

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
  electronHEEPIdMapToken_      (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("ElectronHEEPIdMap"))),
  eleVetoIdCutFlowResultMapToken_   (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleVetoIdCutFlowResultMap"))),
  eleLooseIdCutFlowResultMapToken_  (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleLooseIdCutFlowResultMap"))),
  eleMediumIdCutFlowResultMapToken_ (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleMediumIdCutFlowResultMap"))),
  eleTightIdCutFlowResultMapToken_  (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleTightIdCutFlowResultMap"))),
  eleHEEPIdCutFlowResultMapToken_   (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdCutFlowResultMap"))),
  beamSpotInputToken_          (consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BeamSpotInputTag"))), 
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

  produces <std::vector<double> > ( prefix + "Eta"                      + suffix );
  produces <std::vector<double> > ( prefix + "Phi"                      + suffix );
  produces <std::vector<double> > ( prefix + "Pt"                       + suffix );
  produces <std::vector<double> > ( prefix + "PtHeep"                   + suffix );
  produces <std::vector<double> > ( prefix + "Energy"                   + suffix );
  produces <std::vector<double> > ( prefix + "CaloEnergy"               + suffix );
  produces <std::vector<double> > ( prefix + "EcalEnergy"               + suffix );
  produces <std::vector<int> >    ( prefix + "Charge"                   + suffix );
  produces <std::vector<double> > ( prefix + "HoE"                      + suffix );
								        
  // Supercluster kinematic variables				        
								        
  produces <std::vector<double> > ( prefix + "ESuperClusterOverP"       + suffix );
  produces <std::vector<double> > ( prefix + "SCEta"                    + suffix );
  produces <std::vector<double> > ( prefix + "SCPhi"                    + suffix );
  produces <std::vector<double> > ( prefix + "SCPt"                     + suffix );
  produces <std::vector<double> > ( prefix + "SCRawEnergy"              + suffix );
  produces <std::vector<double> > ( prefix + "SCEnergy"                 + suffix );

  // ID information
  
  produces <std::vector<int> >     ( prefix + "PassId"                     + suffix );
  produces <std::vector<bool> >    ( prefix + "PassEGammaIDVeto"           + suffix );
  produces <std::vector<bool> >    ( prefix + "PassEGammaIDLoose"          + suffix );
  produces <std::vector<bool> >    ( prefix + "PassEGammaIDMedium"         + suffix );
  produces <std::vector<bool> >    ( prefix + "PassEGammaIDTight"          + suffix );
  produces <std::vector<bool> >    ( prefix + "PassHEEPID"                 + suffix );
  produces <std::vector<std::string> >    ( prefix + "CutFlowNamesEGammaIDVeto"   + suffix );
  produces <std::vector<std::string> >    ( prefix + "CutFlowNamesEGammaIDLoose"  + suffix );
  produces <std::vector<std::string> >    ( prefix + "CutFlowNamesEGammaIDMedium" + suffix );
  produces <std::vector<std::string> >    ( prefix + "CutFlowNamesEGammaIDTight"  + suffix );
  produces <std::vector<std::string> >    ( prefix + "CutFlowNamesEGammaIDHEEP"   + suffix );
  produces <std::vector<std::string> >    ( prefix + "CutFlowHashesEGammaIDVeto"  + suffix );
  produces <std::vector<std::string> >    ( prefix + "CutFlowHashesEGammaIDLoose" + suffix );
  produces <std::vector<std::string> >    ( prefix + "CutFlowHashesEGammaIDMedium"+ suffix );
  produces <std::vector<std::string> >    ( prefix + "CutFlowHashesEGammaIDTight" + suffix );
  produces <std::vector<std::string> >    ( prefix + "CutFlowHashesEGammaIDHEEP"  + suffix );
  produces <std::vector<bool> >    ( prefix + "PassEGammaIDTrigTight"      + suffix );
  produces <std::vector<bool> >    ( prefix + "PassEGammaIDTrigWP70"       + suffix );
  produces <std::vector<bool> >    ( prefix + "PassEGammaIDEoP"            + suffix );
  produces <std::vector<float> >   ( prefix + "RhoIsoHEEP"                 + suffix );

  // Does this electron overlap with a muon?			        
  produces <std::vector<int> >    ( prefix + "Overlaps"                 + suffix );
								        
  // Number of Brems = number of basic clusters minus one	        
  produces <std::vector<int> >    ( prefix + "NumberOfBrems"            + suffix );

  // Is this ECAL driven? Or PFlow?  				        
  produces <std::vector<bool> >   ( prefix + "HasEcalDrivenSeed"        + suffix );
  produces <std::vector<bool> >   ( prefix + "HasTrackerDrivenSeed"     + suffix );

  // Charge consistency variables - Ferdinando Giordano

  produces <std::vector<bool> >   ( prefix + "GsfCtfScPixCharge"        + suffix );
  produces <std::vector<bool> >   ( prefix + "GsfScPixCharge"           + suffix );
  produces <std::vector<bool> >   ( prefix + "GsfCtfCharge"             + suffix );

  // EE or EB - Ferdinando Giordano

  produces <std::vector<bool> >   ( prefix + "IsEB"                     + suffix );
  produces <std::vector<bool> >   ( prefix + "IsEE"                     + suffix );
  
  // ECAL eta/phi vs tracker eta/phi				        
								        
  produces <std::vector<double> > ( prefix + "DeltaPhiTrkSC"            + suffix );
  produces <std::vector<double> > ( prefix + "DeltaEtaTrkSC"            + suffix );
  produces <std::vector<double> > ( prefix + "DeltaEtaTrkSeedSC"        + suffix );
								        
  // Shower shape						        
								        
  produces <std::vector<double> > ( prefix + "SigmaEtaEta"              + suffix );
  produces <std::vector<double> > ( prefix + "SigmaIEtaIEta"            + suffix );
  produces <std::vector<double> > ( prefix + "Full5x5SigmaIEtaIEta"     + suffix );
  produces <std::vector<int> >    ( prefix + "Classif"                  + suffix );
  produces <std::vector<double> > ( prefix + "R9"                       + suffix );
  produces <std::vector<double> > ( prefix + "E1x5OverE5x5"             + suffix );
  produces <std::vector<double> > ( prefix + "E2x5OverE5x5"             + suffix );
  produces <std::vector<double> > ( prefix + "Full5x5E1x5OverE5x5"      + suffix );
  produces <std::vector<double> > ( prefix + "Full5x5E2x5OverE5x5"      + suffix );
  								        
  // Isolation variables: PAT					        
								        
  produces <std::vector<double> > ( prefix + "TrkIsoPAT"                + suffix );
  produces <std::vector<double> > ( prefix + "EcalIsoPAT"               + suffix );
  produces <std::vector<double> > ( prefix + "HcalIsoPAT"               + suffix );
  produces <std::vector<double> > ( prefix + "RelIsoPAT"                + suffix );
  produces <std::vector<int> >    ( prefix + "PassIsoPAT"               + suffix );

  // Isolation variables: particle flow
  
  produces <std::vector<double> > ( prefix + "PFChargedHadronIso03"     + suffix );
  produces <std::vector<double> > ( prefix + "PFNeutralHadronIso03"     + suffix );
  produces <std::vector<double> > ( prefix + "PFPhotonIso03"            + suffix );
  produces <std::vector<double> > ( prefix + "PFPUIso03"                + suffix );

  produces <std::vector<double> > ( prefix + "PFChargedHadronIso04"     + suffix );
  produces <std::vector<double> > ( prefix + "PFNeutralHadronIso04"     + suffix );
  produces <std::vector<double> > ( prefix + "PFPhotonIso04"            + suffix );
  
  // Isolation variables: DR 0.3				        
								        
  produces <std::vector<double> > ( prefix + "EcalIsoDR03"              + suffix );
  produces <std::vector<double> > ( prefix + "HcalIsoDR03"              + suffix );
  produces <std::vector<double> > ( prefix + "HcalIsoDR03FullCone"      + suffix );
  produces <std::vector<double> > ( prefix + "HcalIsoD1DR03"            + suffix );
  produces <std::vector<double> > ( prefix + "HcalIsoD2DR03"            + suffix );
  produces <std::vector<double> > ( prefix + "TrkIsoDR03"               + suffix );
				        
  // Conversion variables					        
								        
  produces <std::vector<int> >    ( prefix + "MissingHits"              + suffix );
  produces <std::vector<int> >    ( prefix + "MissingHitsEG"            + suffix );
  produces <std::vector<double> > ( prefix + "Dist"                     + suffix );
  produces <std::vector<double> > ( prefix + "DCotTheta"                + suffix );
  produces <std::vector<double> > ( prefix + "Fbrem"                    + suffix );
  produces <std::vector<bool> >   ( prefix + "HasMatchedConvPhot"       + suffix );

  // Vertex and beamspot information

  produces <std::vector<int> >    ( prefix + "VtxIndex"                 + suffix );
  produces <std::vector<double> > ( prefix + "VtxDistXY"                + suffix );
  produces <std::vector<double> > ( prefix + "VtxDistZ"                 + suffix );
  produces <std::vector<double> > ( prefix + "LeadVtxDistXY"            + suffix );
  produces <std::vector<double> > ( prefix + "LeadVtxDistZ"             + suffix );
  produces <std::vector<double> > ( prefix + "PrimaryVertexDXY"         + suffix );
  produces <std::vector<double> > ( prefix + "PrimaryVertexDXYError"    + suffix );
  produces <std::vector<double> > ( prefix + "BeamSpotDXY"              + suffix );
  produces <std::vector<double> > ( prefix + "BeamSpotDXYError"         + suffix );

  // Track information

  produces <std::vector<double> > ( prefix + "TrackVx"                  + suffix );
  produces <std::vector<double> > ( prefix + "TrackVy"                  + suffix );
  produces <std::vector<double> > ( prefix + "TrackVz"                  + suffix );
  produces <std::vector<double> > ( prefix + "TrackPt"                  + suffix );
  produces <std::vector<double> > ( prefix + "TrackValidFractionOfHits" + suffix );

  // Trigger matching: Double electron

  produces <std::vector<bool  > > ( prefix + "HLTDoubleEleMatched"      + suffix );
  produces <std::vector<double> > ( prefix + "HLTDoubleEleMatchPt"      + suffix );
  produces <std::vector<double> > ( prefix + "HLTDoubleEleMatchEta"     + suffix );
  produces <std::vector<double> > ( prefix + "HLTDoubleEleMatchPhi"     + suffix );

  // Trigger matching: Single electron

  produces <std::vector<bool  > > ( prefix + "HLTSingleEleMatched"      + suffix );
  produces <std::vector<double> > ( prefix + "HLTSingleEleMatchPt"      + suffix );
  produces <std::vector<double> > ( prefix + "HLTSingleEleMatchEta"     + suffix );
  produces <std::vector<double> > ( prefix + "HLTSingleEleMatchPhi"     + suffix );

  // Trigger matching: Single electron (WP85)

  produces <std::vector<bool  > > ( prefix + "HLTSingleEleWP85Matched"  + suffix );
  produces <std::vector<double> > ( prefix + "HLTSingleEleWP85MatchPt"  + suffix );
  produces <std::vector<double> > ( prefix + "HLTSingleEleWP85MatchEta" + suffix );
  produces <std::vector<double> > ( prefix + "HLTSingleEleWP85MatchPhi" + suffix );

  // Trigger matching: Ele+Jet+Jet

  produces <std::vector<bool  > > ( prefix + "HLTEleJetJetMatched"  + suffix );
  produces <std::vector<double> > ( prefix + "HLTEleJetJetMatchPt"  + suffix );
  produces <std::vector<double> > ( prefix + "HLTEleJetJetMatchEta" + suffix );
  produces <std::vector<double> > ( prefix + "HLTEleJetJetMatchPhi" + suffix );

  // Gen matching: status 3 only 

  produces <std::vector<double> > ( prefix + "MatchedGenParticlePt"   + suffix );
  produces <std::vector<double> > ( prefix + "MatchedGenParticleEta"  + suffix );
  produces <std::vector<double> > ( prefix + "MatchedGenParticlePhi"  + suffix );
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

  std::auto_ptr<std::vector<double> >  eta                       ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi                       ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt                        ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ptHeep                    ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy                    ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  caloEnergy                ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ecalEnergy                ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     charge                    ( new std::vector<int>   ()  );
  std::auto_ptr<std::vector<double> >  hoe                       ( new std::vector<double>()  );

  // Supercluster kinematic variables

  std::auto_ptr<std::vector<double> >  eSuperClusterOverP        ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scEta                     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scPhi                     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scPt                      ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scRawEnergy               ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scEnergy                  ( new std::vector<double>()  );

  // ID information
  std::auto_ptr<std::vector<int> >     passIds                   ( new std::vector<int>   ()  );
  std::auto_ptr<std::vector<bool> >    passEGammaIDVeto          ( new std::vector<bool>   ()  );
  std::auto_ptr<std::vector<bool> >    passEGammaIDLoose         ( new std::vector<bool>   ()  );
  std::auto_ptr<std::vector<bool> >    passEGammaIDMedium        ( new std::vector<bool>   ()  );
  std::auto_ptr<std::vector<bool> >    passEGammaIDTight         ( new std::vector<bool>   ()  );
  std::auto_ptr<std::vector<bool> >    passHEEPID                ( new std::vector<bool>   ()  ); 
  std::auto_ptr<std::vector<std::string> > cutFlowNamesEGammaIDVeto   ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowNamesEGammaIDLoose  ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowNamesEGammaIDMedium ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowNamesEGammaIDTight  ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowNamesEGammaIDHEEP   ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowHashesEGammaIDVeto   ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowHashesEGammaIDLoose  ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowHashesEGammaIDMedium ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowHashesEGammaIDTight  ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<std::string> > cutFlowHashesEGammaIDHEEP   ( new std::vector<std::string> () );
  std::auto_ptr<std::vector<bool> >    passEGammaIDTrigTight     ( new std::vector<bool>   ()  );
  std::auto_ptr<std::vector<bool> >    passEGammaIDTrigWP70      ( new std::vector<bool>   ()  );
  std::auto_ptr<std::vector<bool> >    passEGammaIDEoP           ( new std::vector<bool>   ()  ); 
  std::auto_ptr<std::vector<float> >   rhoIsoHEEP                ( new std::vector<float>   ()  ); 
  
  // Does this electron overlap with a muon?
  std::auto_ptr<std::vector<int> >     overlaps                  ( new std::vector<int>   ()  );

  // Number of Brems = number of basic clusters minus one
  std::auto_ptr<std::vector<int> >     numberOfBrems             ( new std::vector<int>   ()  );

  // Is this ECAL driven? Or PFlow?  
  std::auto_ptr<std::vector<bool> >    hasEcalDrivenSeed         ( new std::vector<bool>  ()  );
  std::auto_ptr<std::vector<bool> >    hasTrackerDrivenSeed      ( new std::vector<bool>  ()  );
  
  // Charge consistency variables: Ferdinando Giordano
  std::auto_ptr<std::vector<bool> >    gsfCtfScPixCharge         ( new std::vector<bool>  ()  );
  std::auto_ptr<std::vector<bool> >    gsfScPixCharge            ( new std::vector<bool>  ()  );
  std::auto_ptr<std::vector<bool> >    gsfCtfCharge              ( new std::vector<bool>  ()  );
  
  // EB or EE: Ferdinando Giordano
  std::auto_ptr<std::vector<bool> >    isEB                      ( new std::vector<bool>  ()  );
  std::auto_ptr<std::vector<bool> >    isEE                      ( new std::vector<bool>  ()  );
  
  // ECAL eta/phi vs tracker eta/phi

  std::auto_ptr<std::vector<double> >  deltaPhiTrkSC             ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  deltaEtaTrkSC             ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  deltaEtaTrkSeedSC         ( new std::vector<double>()  );

  // Shower shape

  std::auto_ptr<std::vector<double> >  sigmaEtaEta               ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  sigmaIEtaIEta             ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  full5x5SigmaIEtaIEta      ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  r9                        ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  e1x5overe5x5              ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  e2x5overe5x5              ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  full5x5e1x5overe5x5       ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  full5x5e2x5overe5x5       ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     classif                   ( new std::vector<int>   ()  );

  // Isolation variables: PAT
  
  std::auto_ptr<std::vector<double> >  trkIsoPAT                 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ecalIsoPAT                ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIsoPAT                ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  relIsoPAT                 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     passIsoPAT                ( new std::vector<int>   ()  );

  // Isolation variables: particle flow 
  
  std::auto_ptr<std::vector<double> >  pfChargedHadronIso03      ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfNeutralHadronIso03      ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfPhotonIso03             ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfPUIso03                 ( new std::vector<double>()  );

  std::auto_ptr<std::vector<double> >  pfChargedHadronIso04      ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfNeutralHadronIso04      ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfPhotonIso04             ( new std::vector<double>()  );
  
  // Isolation variables: DR 0.3

  std::auto_ptr<std::vector<double> >  ecalIsoDR03               ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIsoDR03               ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIsoDR03FullCone       ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIsoD1DR03             ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIsoD2DR03             ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkIsoDR03                ( new std::vector<double>()  );

  // Conversion variables

  std::auto_ptr<std::vector<int> >     missingHits               ( new std::vector<int>   ()  );
  std::auto_ptr<std::vector<int> >     missingHitsEG             ( new std::vector<int>   ()  );
  std::auto_ptr<std::vector<double> >  dist_vec                  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  dCotTheta                 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  fbrem                     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<bool> >    hasMatchedConvPhot        ( new std::vector<bool>  ()  );
  
  // Vertex and beamspot information
  
  std::auto_ptr<std::vector<int> >     vtxIndex                  ( new std::vector<int>   ()  );
  std::auto_ptr<std::vector<double> >  vtxDistXY                 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  vtxDistZ                  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  vtx0DistXY                ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  vtx0DistZ                 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  primaryVertexDXY          ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  primaryVertexDXYError     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  beamspotDXY               ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  beamspotDXYError          ( new std::vector<double>()  );

  // Track information 

  std::auto_ptr<std::vector<double> >  trackVx                   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trackVy                   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trackVz                   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trackPt                   ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trackValidFractionOfHits  ( new std::vector<double>()  );  

  // Trigger matching: Double electron

  std::auto_ptr<std::vector<bool  > >  HLTDoubleEleMatched       ( new std::vector<bool  >()  );
  std::auto_ptr<std::vector<double> >  HLTDoubleEleMatchPt 	 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  HLTDoubleEleMatchEta	 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  HLTDoubleEleMatchPhi      ( new std::vector<double>()  );

  // Trigger matching: Single electron

  std::auto_ptr<std::vector<bool  > >  HLTSingleEleMatched       ( new std::vector<bool  >()  );
  std::auto_ptr<std::vector<double> >  HLTSingleEleMatchPt 	 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  HLTSingleEleMatchEta	 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  HLTSingleEleMatchPhi      ( new std::vector<double>()  );

  // Trigger matching: Single electron (WP85)

  std::auto_ptr<std::vector<bool  > >  HLTSingleEleWP85Matched   ( new std::vector<bool  >()  );
  std::auto_ptr<std::vector<double> >  HLTSingleEleWP85MatchPt 	 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  HLTSingleEleWP85MatchEta	 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  HLTSingleEleWP85MatchPhi  ( new std::vector<double>()  );

  // Trigger matching: Ele+Jet+Jet

  std::auto_ptr<std::vector<bool  > >  HLTEleJetJetMatched   ( new std::vector<bool  >()  );
  std::auto_ptr<std::vector<double> >  HLTEleJetJetMatchPt 	 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  HLTEleJetJetMatchEta	 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  HLTEleJetJetMatchPhi  ( new std::vector<double>()  );

  // Gen matching: Status 3 only
  
  std::auto_ptr<std::vector<double> >  matchedGenParticlePt  ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  matchedGenParticleEta ( new std::vector<double>()   );
  std::auto_ptr<std::vector<double> >  matchedGenParticlePhi ( new std::vector<double>()   );
  
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
  double rhoIso = *(rho.product());
  
  // SIC add for new egamma VID framework
  // example: https://github.com/ikrav/ElectronWork/blob/master/ElectronNtupler/plugins/ElectronNtuplerIdDemoPrePHYS14miniAOD.cc
  // get electron ID maps
  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > veto_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > loose_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > medium_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > tight_id_cutflow_data;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > heep_id_cutflow_data;
  iEvent.getByToken(electronVetoIdMapToken_,veto_id_decisions);
  iEvent.getByToken(electronLooseIdMapToken_,loose_id_decisions);
  iEvent.getByToken(electronMediumIdMapToken_,medium_id_decisions);
  iEvent.getByToken(electronTightIdMapToken_,tight_id_decisions);
  iEvent.getByToken(electronHEEPIdMapToken_,heep_id_decisions);
  iEvent.getByToken(eleVetoIdCutFlowResultMapToken_,veto_id_cutflow_data);
  iEvent.getByToken(eleLooseIdCutFlowResultMapToken_,loose_id_cutflow_data);
  iEvent.getByToken(eleMediumIdCutFlowResultMapToken_,medium_id_cutflow_data);
  iEvent.getByToken(eleTightIdCutFlowResultMapToken_,tight_id_cutflow_data);
  iEvent.getByToken(eleHEEPIdCutFlowResultMapToken_,heep_id_cutflow_data);

  //------------------------------------------------------------------------
  // Get magnetic field (need this for photon conversion information)
  //  - Instead, use conversion information built into PAT/MiniAOD, so no need for this
  //------------------------------------------------------------------------

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

      // Double electron
      const pat::TriggerObjectStandAloneCollection matchesDoubleEle = it->triggerObjectMatchesByPath("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*");
      if(matchesDoubleEle.size() > 0)
      {
        HLTDoubleEleMatched  -> push_back ( true ) ;
        HLTDoubleEleMatchPt  -> push_back ( matchesDoubleEle[0].pt() );
        HLTDoubleEleMatchEta -> push_back ( matchesDoubleEle[0].eta());
        HLTDoubleEleMatchPhi -> push_back ( matchesDoubleEle[0].phi());
      }
      else
      {
        HLTDoubleEleMatched  -> push_back ( false ) ;
        HLTDoubleEleMatchPt  -> push_back ( -999. );
        HLTDoubleEleMatchEta -> push_back ( -999. );
        HLTDoubleEleMatchPhi -> push_back ( -999. );
      }

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

      // Single electron (WP85)
      const pat::TriggerObjectStandAloneCollection matchesSingleEleWP85 = it->triggerObjectMatchesByPath("HLT_Ele32_eta2p1_WP85_Gsf_v*");
      if(matchesSingleEleWP85.size() > 0)
      {
        HLTSingleEleWP85Matched  -> push_back ( true ) ;
        HLTSingleEleWP85MatchPt  -> push_back ( matchesSingleEleWP85[0].pt() );
        HLTSingleEleWP85MatchEta -> push_back ( matchesSingleEleWP85[0].eta());
        HLTSingleEleWP85MatchPhi -> push_back ( matchesSingleEleWP85[0].phi());
      }
      else 
      { 
        HLTSingleEleWP85Matched  -> push_back ( false ) ;
        HLTSingleEleWP85MatchPt  -> push_back ( -999. );
        HLTSingleEleWP85MatchEta -> push_back ( -999. );
        HLTSingleEleWP85MatchPhi -> push_back ( -999. );
      }

      // E+J+J cross trigger
      const pat::TriggerObjectStandAloneCollection matchesEleJetJet = it->triggerObjectMatchesByPath("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*");
      if(matchesEleJetJet.size() > 0)
      {
        HLTEleJetJetMatched  -> push_back ( true ) ;
        HLTEleJetJetMatchPt  -> push_back ( matchesEleJetJet[0].pt() );
        HLTEleJetJetMatchEta -> push_back ( matchesEleJetJet[0].eta());
        HLTEleJetJetMatchPhi -> push_back ( matchesEleJetJet[0].phi());
      }
      else 
      { 
        HLTEleJetJetMatched  -> push_back ( false ) ;
        HLTEleJetJetMatchPt  -> push_back ( -999. );
        HLTEleJetJetMatchEta -> push_back ( -999. );
        HLTEleJetJetMatchPhi -> push_back ( -999. );
      }

      //------------------------------------------------------------------------
      // Gen matching
      // This should work for pythia8 in MiniAOD. should use status 23 I think for outgoing electrons.
      // See: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#MC_Truth
      //------------------------------------------------------------------------

      double genPartPt = -999.;
      double genPartEta= -999.;
      double genPartPhi= -999.;
      
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
      
      matchedGenParticlePt  -> push_back ( (double)(genPartPt ) );
      matchedGenParticleEta -> push_back ( (double)(genPartEta) );
      matchedGenParticlePhi -> push_back ( (double)(genPartPhi) );
      
      //------------------------------------------------------------------------
      // Relative isolation (not currently used in any analysis... remove?) 
      //------------------------------------------------------------------------
      
      double reliso = (it->trackIso() + it->ecalIso() + it->hcalIso())/it->pt();

      //------------------------------------------------------------------------
      // Conversion information
      //------------------------------------------------------------------------
      // Now take what's embedded in the pat candidate

      //------------------------------------------------------------------------
      // Vertex association
      //------------------------------------------------------------------------

      double minVtxDist3D = 9999.;
      int vtxIndex_ = -1;
      double vtxDistXY_ = -9999.;
      double vtxDistZ_ = -9999.;
      
      double vtx0DistXY_;
      double vtx0DistZ_;

      if(primaryVertices.isValid()) {
        edm::LogInfo("RootTupleMakerV2_ElectronsInfo") << "Total # Primary Vertices: " << primaryVertices->size();

        int i_vertex = 0;
        for( reco::VertexCollection::const_iterator v_it=primaryVertices->begin() ; v_it!=primaryVertices->end() ; ++v_it ) {

          double distXY = it->gsfTrack()->dxy(v_it->position());
          double distZ = it->gsfTrack()->dz(v_it->position());
          double dist3D = sqrt(pow(distXY,2) + pow(distZ,2));

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
      scPhi                    -> push_back( it->superCluster()->phi() );
      scPt                     -> push_back( it->superCluster()->energy()/cosh(it->superCluster()->eta()) );
      scRawEnergy              -> push_back( it->superCluster()->rawEnergy() );
      scEnergy                 -> push_back( it->superCluster()->energy() );
      eSuperClusterOverP       -> push_back( it->eSuperClusterOverP() );

      // ID information
      passIds                  -> push_back( passId );
      // SIC: update to Run2 cut-based and HEEP ID's
      const edm::Ptr<pat::Electron> elPtr(electrons, it - electrons->begin() );
      passEGammaIDVeto           -> push_back ((*veto_id_decisions)[ elPtr ]);
      passEGammaIDLoose          -> push_back ((*loose_id_decisions)[ elPtr ]);
      passEGammaIDMedium         -> push_back ((*medium_id_decisions)[ elPtr ]);
      passEGammaIDTight          -> push_back ((*tight_id_decisions)[ elPtr ]);
      passHEEPID                 -> push_back ((*heep_id_decisions)[ elPtr ]);
      cutFlowNamesEGammaIDVeto   -> push_back (((*veto_id_cutflow_data)[ elPtr ]).cutFlowName());
      cutFlowNamesEGammaIDLoose  -> push_back (((*loose_id_cutflow_data)[ elPtr ]).cutFlowName());
      cutFlowNamesEGammaIDMedium -> push_back (((*medium_id_cutflow_data)[ elPtr ]).cutFlowName());
      cutFlowNamesEGammaIDTight  -> push_back (((*tight_id_cutflow_data)[ elPtr ]).cutFlowName());
      cutFlowNamesEGammaIDHEEP   -> push_back (((*heep_id_cutflow_data)[ elPtr ]).cutFlowName());
      cutFlowHashesEGammaIDVeto  -> push_back (((*veto_id_cutflow_data)[ elPtr ]).cutFlowHash());
      cutFlowHashesEGammaIDLoose -> push_back (((*loose_id_cutflow_data)[ elPtr ]).cutFlowHash());
      cutFlowHashesEGammaIDMedium-> push_back (((*medium_id_cutflow_data)[ elPtr ]).cutFlowHash());
      cutFlowHashesEGammaIDTight -> push_back (((*tight_id_cutflow_data)[ elPtr ]).cutFlowHash());
      cutFlowHashesEGammaIDHEEP  -> push_back (((*heep_id_cutflow_data)[ elPtr ]).cutFlowHash());
      //
      rhoIsoHEEP               -> push_back (rhoIso);
      // XXX FIXME SIC: update with trigger updates?
      //passEGammaIDTrigTight    -> push_back (EgammaCutBasedEleId::PassTriggerCuts(EgammaCutBasedEleId::TRIGGERTIGHT, *it));
      //passEGammaIDTrigWP70     -> push_back (EgammaCutBasedEleId::PassTriggerCuts(EgammaCutBasedEleId::TRIGGERWP70 , *it));
      //passEGammaIDEoP          -> push_back (EgammaCutBasedEleId::PassEoverPCuts(*it));

      // Does this electron overlap with a muon?
      overlaps                 -> push_back( ovrlps );
      
      // Number of Brems = number of basic clusters minus one
      numberOfBrems            -> push_back( it->numberOfBrems() );
      
      // Is this ECAL driven? Or PFlow?  
      hasEcalDrivenSeed        -> push_back( it->ecalDrivenSeed() );
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
      trackVx                  -> push_back( it->gsfTrack()->vx() );
      trackVy                  -> push_back( it->gsfTrack()->vy() );
      trackVz                  -> push_back( it->gsfTrack()->vz() );
      trackPt                  -> push_back( it->gsfTrack()->pt() );
      trackValidFractionOfHits -> push_back( it->gsfTrack()->validFraction() );


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
  iEvent.put( scPhi                   , prefix + "SCPhi"                    + suffix );
  iEvent.put( scPt                    , prefix + "SCPt"                     + suffix );
  iEvent.put( scRawEnergy             , prefix + "SCRawEnergy"              + suffix );
  iEvent.put( scEnergy                , prefix + "SCEnergy"                 + suffix );

  // ID information 
  iEvent.put( passIds                 , prefix + "PassId"                   + suffix );
  iEvent.put( passEGammaIDVeto        , prefix + "PassEGammaIDVeto"         + suffix );
  iEvent.put( passEGammaIDLoose       , prefix + "PassEGammaIDLoose"        + suffix );
  iEvent.put( passEGammaIDMedium      , prefix + "PassEGammaIDMedium"       + suffix );
  iEvent.put( passEGammaIDTight       , prefix + "PassEGammaIDTight"        + suffix );
  iEvent.put( passEGammaIDTrigTight   , prefix + "PassEGammaIDTrigTight"    + suffix );
  iEvent.put( passEGammaIDTrigWP70    , prefix + "PassEGammaIDTrigWP70"     + suffix );
  iEvent.put( passEGammaIDEoP         , prefix + "PassEGammaIDEoP"          + suffix );
  iEvent.put( passHEEPID              , prefix + "PassHEEPID"               + suffix );
  iEvent.put( cutFlowNamesEGammaIDVeto   , prefix + "CutFlowNamesEGammaIDVeto"   + suffix );
  iEvent.put( cutFlowNamesEGammaIDLoose  , prefix + "CutFlowNamesEGammaIDLoose"  + suffix );
  iEvent.put( cutFlowNamesEGammaIDMedium , prefix + "CutFlowNamesEGammaIDMedium" + suffix );
  iEvent.put( cutFlowNamesEGammaIDTight  , prefix + "CutFlowNamesEGammaIDTight"  + suffix );
  iEvent.put( cutFlowNamesEGammaIDHEEP   , prefix + "CutFlowNamesEGammaIDHEEP"   + suffix );
  iEvent.put( cutFlowHashesEGammaIDVeto  , prefix + "CutFlowHashesEGammaIDVeto"  + suffix );
  iEvent.put( cutFlowHashesEGammaIDLoose , prefix + "CutFlowHashesEGammaIDLoose" + suffix );
  iEvent.put( cutFlowHashesEGammaIDMedium, prefix + "CutFlowHashesEGammaIDMedium"+ suffix );
  iEvent.put( cutFlowHashesEGammaIDTight , prefix + "CutFlowHashesEGammaIDTight" + suffix );
  iEvent.put( cutFlowHashesEGammaIDHEEP  , prefix + "CutFlowHashesEGammaIDHEEP"  + suffix );
  iEvent.put( rhoIsoHEEP              , prefix + "RhoIsoHEEP"               + suffix );
  
  // Does this electron overlap with a muon?			        
  iEvent.put( overlaps                , prefix + "Overlaps"                 + suffix );

  // Number of Brems = number of basic clusters minus one	        
  iEvent.put( numberOfBrems           , prefix + "NumberOfBrems"            + suffix );
  
  // Is this ECAL driven? Or PFlow?  				        
  iEvent.put( hasEcalDrivenSeed       , prefix + "HasEcalDrivenSeed"        + suffix );
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
  iEvent.put( sigmaIEtaIEta           , prefix + "SigmaIEtaIEta"            + suffix );
  iEvent.put( full5x5SigmaIEtaIEta    , prefix + "Full5x5SigmaIEtaIEta"     + suffix );
  iEvent.put( r9                      , prefix + "R9"                       + suffix );
  iEvent.put( e1x5overe5x5            , prefix + "E1x5OverE5x5"             + suffix );
  iEvent.put( e2x5overe5x5            , prefix + "E2x5OverE5x5"             + suffix );
  iEvent.put( full5x5e1x5overe5x5        , prefix + "Full5x5E1x5OverE5x5"         + suffix );
  iEvent.put( full5x5e2x5overe5x5        , prefix + "Full5x5E2x5OverE5x5"         + suffix );
  iEvent.put( classif                 , prefix + "Classif"                  + suffix );
  
  // Isolation variables: PAT					        

  iEvent.put( trkIsoPAT               , prefix + "TrkIsoPAT"                + suffix );
  iEvent.put( ecalIsoPAT              , prefix + "EcalIsoPAT"               + suffix );
  iEvent.put( hcalIsoPAT              , prefix + "HcalIsoPAT"               + suffix );
  iEvent.put( relIsoPAT               , prefix + "RelIsoPAT"                + suffix );
  iEvent.put( passIsoPAT              , prefix + "PassIsoPAT"               + suffix );

  // Isolation variables: DR 0.3				        

  iEvent.put( ecalIsoDR03             , prefix + "EcalIsoDR03"              + suffix );
  iEvent.put( hcalIsoDR03             , prefix + "HcalIsoDR03"              + suffix );
  iEvent.put( hcalIsoDR03FullCone     , prefix + "HcalIsoDR03FullCone"      + suffix );
  iEvent.put( hcalIsoD1DR03           , prefix + "HcalIsoD1DR03"            + suffix );
  iEvent.put( hcalIsoD2DR03           , prefix + "HcalIsoD2DR03"            + suffix );
  iEvent.put( trkIsoDR03              , prefix + "TrkIsoDR03"               + suffix );

  // Isolation variables: particle flow

  iEvent.put( pfChargedHadronIso03    , prefix + "PFChargedHadronIso03"     + suffix );
  iEvent.put( pfNeutralHadronIso03    , prefix + "PFNeutralHadronIso03"     + suffix );
  iEvent.put( pfPhotonIso03           , prefix + "PFPhotonIso03"            + suffix );
  iEvent.put( pfPUIso03               , prefix + "PFPUIso03"                + suffix );

  iEvent.put( pfChargedHadronIso04    , prefix + "PFChargedHadronIso04"     + suffix );
  iEvent.put( pfNeutralHadronIso04    , prefix + "PFNeutralHadronIso04"     + suffix );
  iEvent.put( pfPhotonIso04           , prefix + "PFPhotonIso04"            + suffix );

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

  iEvent.put( trackVx                 , prefix + "TrackVx"                  + suffix );
  iEvent.put( trackVy                 , prefix + "TrackVy"                  + suffix );
  iEvent.put( trackVz                 , prefix + "TrackVz"                  + suffix );
  iEvent.put( trackPt                 , prefix + "TrackPt"                  + suffix );
  iEvent.put( trackValidFractionOfHits, prefix + "TrackValidFractionOfHits" + suffix );

  // Trigger matching: Double electron

  iEvent.put( HLTDoubleEleMatched     , prefix + "HLTDoubleEleMatched"      + suffix );
  iEvent.put( HLTDoubleEleMatchPt     , prefix + "HLTDoubleEleMatchPt"      + suffix );
  iEvent.put( HLTDoubleEleMatchEta    , prefix + "HLTDoubleEleMatchEta"     + suffix );
  iEvent.put( HLTDoubleEleMatchPhi    , prefix + "HLTDoubleEleMatchPhi"     + suffix );

  // Trigger matching: Single electron

  iEvent.put( HLTSingleEleMatched     , prefix + "HLTSingleEleMatched"      + suffix );
  iEvent.put( HLTSingleEleMatchPt     , prefix + "HLTSingleEleMatchPt"      + suffix );
  iEvent.put( HLTSingleEleMatchEta    , prefix + "HLTSingleEleMatchEta"     + suffix );
  iEvent.put( HLTSingleEleMatchPhi    , prefix + "HLTSingleEleMatchPhi"     + suffix );

  // Trigger matching: Single electron (WP85)

  iEvent.put( HLTSingleEleWP85Matched , prefix + "HLTSingleEleWP85Matched"  + suffix );
  iEvent.put( HLTSingleEleWP85MatchPt , prefix + "HLTSingleEleWP85MatchPt"  + suffix );
  iEvent.put( HLTSingleEleWP85MatchEta, prefix + "HLTSingleEleWP85MatchEta" + suffix );
  iEvent.put( HLTSingleEleWP85MatchPhi, prefix + "HLTSingleEleWP85MatchPhi" + suffix );

  // Trigger matching: Ele+Jet+Jet

  iEvent.put( HLTEleJetJetMatched , prefix + "HLTEleJetJetMatched"  + suffix );
  iEvent.put( HLTEleJetJetMatchPt , prefix + "HLTEleJetJetMatchPt"  + suffix );
  iEvent.put( HLTEleJetJetMatchEta, prefix + "HLTEleJetJetMatchEta" + suffix );
  iEvent.put( HLTEleJetJetMatchPhi, prefix + "HLTEleJetJetMatchPhi" + suffix );

  // Gen matching: Status 3 only

  iEvent.put( matchedGenParticlePt ,   prefix + "MatchedGenParticlePt"      + suffix );
  iEvent.put( matchedGenParticleEta,   prefix + "MatchedGenParticleEta"     + suffix );
  iEvent.put( matchedGenParticlePhi,   prefix + "MatchedGenParticlePhi"     + suffix );

}
