//------------------------------------------------------------------------
// To do:
// - Muon isolation (for overlaps) is done with relative isolation.  Is this correct?
//------------------------------------------------------------------------

//------------------------------------------------------------------------
// Include files
//------------------------------------------------------------------------

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
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

//------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------s

RootTupleMakerV2_Electrons::RootTupleMakerV2_Electrons(const edm::ParameterSet& iConfig) :
  trkInputTag              (iConfig.getParameter<edm::InputTag>("TracksInputTag"           )),
  dcsInputTag              (iConfig.getParameter<edm::InputTag>("DCSInputTag"              )),
  inputTag                 (iConfig.getParameter<edm::InputTag>("InputTag"                 )),
  vtxInputTag              (iConfig.getParameter<edm::InputTag>("VertexInputTag"           )), 
  beamSpotInputTag         (iConfig.getParameter<edm::InputTag>("BeamSpotInputTag"         )),
  conversionsInputTag      (iConfig.getParameter<edm::InputTag>("ConversionsInputTag"      )),
  triggerEventInputTag     (iConfig.getParameter<edm::InputTag>("TriggerEventInputTag"     )),
  electronIso              (iConfig.getParameter<double>       ("ElectronIso"              )),
  muonPt                   (iConfig.getParameter<double>       ("MuonPt"                   )),
  muonIso                  (iConfig.getParameter<double>       ("MuonIso"                  )),
  muonID                   (iConfig.getParameter<std::string>  ("MuonID"                   )),
  singleEleTriggerMatch    (iConfig.getParameter<std::string>  ("SingleEleTriggerMatch"    )),
  singleEleTriggerMatchWP80(iConfig.getParameter<std::string>  ("SingleEleTriggerMatchWP80")),
  doubleEleTriggerMatch    (iConfig.getParameter<std::string>  ("DoubleEleTriggerMatch"    )),
  prefix                   (iConfig.getParameter<std::string>  ("Prefix"                   )),
  suffix                   (iConfig.getParameter<std::string>  ("Suffix"                   )),
  maxSize                  (iConfig.getParameter<unsigned int> ("MaxSize"                  ))
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
  produces <std::vector<int> >    ( prefix + "Charge"                   + suffix );
  produces <std::vector<double> > ( prefix + "HoE"                      + suffix );
								        
  // Supercluster kinematic variables				        
								        
  produces <std::vector<double> > ( prefix + "ESuperClusterOverP"       + suffix );
  produces <std::vector<double> > ( prefix + "SCEta"                    + suffix );
  produces <std::vector<double> > ( prefix + "SCPhi"                    + suffix );
  produces <std::vector<double> > ( prefix + "SCPt"                     + suffix );
  produces <std::vector<double> > ( prefix + "SCRawEnergy"              + suffix );

  // ID information
  produces <std::vector<int> >    ( prefix + "PassId"                   + suffix );
  								        
  // Does this electron overlap with a muon?			        
  produces <std::vector<int> >    ( prefix + "Overlaps"                 + suffix );
								        
  // Number of Brems = number of basic clusters minus one	        
  produces <std::vector<int> >    ( prefix + "NumberOfBrems"            + suffix );
								        
  // Is this ECAL driven? Or PFlow?  				        
  produces <std::vector<bool> >   ( prefix + "HasEcalDrivenSeed"        + suffix );
								        
  // ECAL eta/phi vs tracker eta/phi				        
								        
  produces <std::vector<double> > ( prefix + "DeltaPhiTrkSC"            + suffix );
  produces <std::vector<double> > ( prefix + "DeltaEtaTrkSC"            + suffix );
								        
  // Shower shape						        
								        
  produces <std::vector<double> > ( prefix + "SigmaEtaEta"              + suffix );
  produces <std::vector<double> > ( prefix + "SigmaIEtaIEta"            + suffix );
  produces <std::vector<int> >    ( prefix + "Classif"                  + suffix );
  produces <std::vector<double> > ( prefix + "E1x5OverE5x5"             + suffix );
  produces <std::vector<double> > ( prefix + "E2x5OverE5x5"             + suffix );
  								        
  // Isolation variables: PAT					        
								        
  produces <std::vector<double> > ( prefix + "TrkIsoPAT"                + suffix );
  produces <std::vector<double> > ( prefix + "EcalIsoPAT"               + suffix );
  produces <std::vector<double> > ( prefix + "HcalIsoPAT"               + suffix );
  produces <std::vector<double> > ( prefix + "RelIsoPAT"                + suffix );
  produces <std::vector<int> >    ( prefix + "PassIsoPAT"               + suffix );

  // Isolation variables: particle flow
  
  produces <std::vector<double> > ( prefix + "PFParticleIso"            + suffix );
  produces <std::vector<double> > ( prefix + "PFChargedHadronIso"       + suffix );
  produces <std::vector<double> > ( prefix + "PFNeutralHadronIso"       + suffix );
  produces <std::vector<double> > ( prefix + "PFPhotonIso"              + suffix );
  produces <std::vector<double> > ( prefix + "PFPUChargedHadronIso"     + suffix );
  
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

  // Trigger matching: Single electron (WP80)

  produces <std::vector<bool  > > ( prefix + "HLTSingleEleWP80Matched"  + suffix );
  produces <std::vector<double> > ( prefix + "HLTSingleEleWP80MatchPt"  + suffix );
  produces <std::vector<double> > ( prefix + "HLTSingleEleWP80MatchEta" + suffix );
  produces <std::vector<double> > ( prefix + "HLTSingleEleWP80MatchPhi" + suffix );

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
  std::auto_ptr<std::vector<int> >     charge                    ( new std::vector<int>   ()  );
  std::auto_ptr<std::vector<double> >  hoe                       ( new std::vector<double>()  );

  // Supercluster kinematic variables

  std::auto_ptr<std::vector<double> >  eSuperClusterOverP        ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scEta                     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scPhi                     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scPt                      ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  scRawEnergy               ( new std::vector<double>()  );

  // ID information
  std::auto_ptr<std::vector<int > >    passIds                   ( new std::vector<int>   ()  );
  
  // Does this electron overlap with a muon?
  std::auto_ptr<std::vector<int> >     overlaps                  ( new std::vector<int>   ()  );

  // Number of Brems = number of basic clusters minus one
  std::auto_ptr<std::vector<int> >     numberOfBrems             ( new std::vector<int>   ()  );

  // Is this ECAL driven? Or PFlow?  
  std::auto_ptr<std::vector<bool> >    hasEcalDrivenSeed         ( new std::vector<bool>  ()  );
  
  // ECAL eta/phi vs tracker eta/phi

  std::auto_ptr<std::vector<double> >  deltaPhiTrkSC             ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  deltaEtaTrkSC             ( new std::vector<double>()  );

  // Shower shape

  std::auto_ptr<std::vector<double> >  sigmaEtaEta               ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  sigmaIEtaIEta             ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  e1x5overe5x5              ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  e2x5overe5x5              ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     classif                   ( new std::vector<int>   ()  );

  // Isolation variables: PAT
  
  std::auto_ptr<std::vector<double> >  trkIsoPAT                 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ecalIsoPAT                ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIsoPAT                ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  relIsoPAT                 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     passIsoPAT                ( new std::vector<int>   ()  );

  // Isolation variables: particle flow 
  
  std::auto_ptr<std::vector<double> >  pfParticleIso             ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfChargedHadronIso        ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfNeutralHadronIso        ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfPhotonIso               ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pfPUChargedHadronIso      ( new std::vector<double>()  );
  
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

  // Trigger matching: Single electron (WP80)

  std::auto_ptr<std::vector<bool  > >  HLTSingleEleWP80Matched   ( new std::vector<bool  >()  );
  std::auto_ptr<std::vector<double> >  HLTSingleEleWP80MatchPt 	 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  HLTSingleEleWP80MatchEta	 ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  HLTSingleEleWP80MatchPhi  ( new std::vector<double>()  );
  
  //------------------------------------------------------------------------
  // Get handles for the event
  //------------------------------------------------------------------------
  
  // Tracks

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(trkInputTag, tracks);

  // DCS information

  edm::Handle<DcsStatusCollection> dcsHandle;
  iEvent.getByLabel(dcsInputTag, dcsHandle);
  
  // Primary vertices
  
  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(vtxInputTag,primaryVertices);

  // Beamspot information
  
  edm::Handle<reco::BeamSpot> bsHandle; 
  iEvent.getByLabel(beamSpotInputTag, bsHandle); 

  // Photon conversion information

  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel(conversionsInputTag, hConversions); 

  // PAT electrons

  edm::Handle<std::vector<pat::Electron> > electrons;
  iEvent.getByLabel(inputTag, electrons);

  // PAT trigger event

  edm::Handle< pat::TriggerEvent > triggerEvent;
  iEvent.getByLabel( triggerEventInputTag, triggerEvent );

  // PAT trigger helper for trigger matching information

  const pat::helper::TriggerMatchHelper matchHelper;
  
  //------------------------------------------------------------------------
  // Get magnetic field (need this for photon conversion information):
  //   - if isRealData then derive bfield using the magnet current from DcsStatus
  //     scale factor = 3.801 Tesla /18166.0 Amps, which are
  //     average values taken over a stable two-week period
  //   - otherwise take it from the IdealMagneticFieldRecord
  //------------------------------------------------------------------------

  double evt_bField = 3.8;
  double currentToBFieldScaleFactor = 2.09237036221512717e-04;

  if(iEvent.isRealData()) {
    if(dcsHandle.isValid()) {
      edm::LogInfo("RootTupleMakerV2_ElectronsInfo") << "Successfully obtained " << dcsInputTag;
      if( (*dcsHandle).size()>0 ) {
        double current = (*dcsHandle)[0].magnetCurrent();
        evt_bField = current*currentToBFieldScaleFactor;
      }
    } 
    else edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get the product " << dcsInputTag;
  } 

  else {
    edm::ESHandle<MagneticField> magneticField;
    iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

    if(magneticField.isValid()) {
      edm::LogInfo("RootTupleMakerV2_ElectronsInfo") << "Successfully obtained IdealMagneticFieldRecord";
      evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
    } 

    else edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get IdealMagneticFieldRecord";
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
      //------------------------------------------------------------------------
      
      int passId = 0;
      if (it->electronID("eidRobustLoose"     )>0) passId = passId | 1<<0;
      if (it->electronID("eidRobustTight"     )>0) passId = passId | 1<<1;
      if (it->electronID("eidLoose"           )>0) passId = passId | 1<<2;
      if (it->electronID("eidTight"           )>0) passId = passId | 1<<3;
      if (it->electronID("eidRobustHighEnergy")>0) passId = passId | 1<<4;
      if (it->electronID("mvaTrigV0"          )>0) passId = passId | 1<<5;
      if (it->electronID("mvaNonTrigV0"       )>0) passId = passId | 1<<6;
      
      //------------------------------------------------------------------------
      // Trigger matching
      // 
      // Example taken from PatTriggerAnalyzer:
      // http://cmslxr.fnal.gov/lxr/source/PhysicsTools/PatExamples/plugins/PatTriggerAnalyzer.cc
      //------------------------------------------------------------------------

      // Double electron

      const pat::TriggerObjectRef doubleElectronTrigRef( matchHelper.triggerMatchObject( electrons, iElectron, doubleEleTriggerMatch, iEvent, *triggerEvent ) );
      if ( doubleElectronTrigRef.isAvailable() && doubleElectronTrigRef.isNonnull() ) { 
	HLTDoubleEleMatched  -> push_back ( true ) ;
	HLTDoubleEleMatchPt  -> push_back ( doubleElectronTrigRef -> pt() );
	HLTDoubleEleMatchEta -> push_back ( doubleElectronTrigRef -> eta());
	HLTDoubleEleMatchPhi -> push_back ( doubleElectronTrigRef -> phi());
      } else { 
	HLTDoubleEleMatched  -> push_back ( false ) ;
	HLTDoubleEleMatchPt  -> push_back ( -999. );
	HLTDoubleEleMatchEta -> push_back ( -999. );
	HLTDoubleEleMatchPhi -> push_back ( -999. );
      }

      // Single electron
      
      const pat::TriggerObjectRef singleElectronTrigRef( matchHelper.triggerMatchObject( electrons, iElectron,  singleEleTriggerMatch, iEvent, *triggerEvent ) );
      if ( singleElectronTrigRef.isAvailable() && singleElectronTrigRef.isNonnull() ) { 
	HLTSingleEleMatched  -> push_back ( true ) ;
	HLTSingleEleMatchPt  -> push_back ( singleElectronTrigRef -> pt() );
	HLTSingleEleMatchEta -> push_back ( singleElectronTrigRef -> eta());
	HLTSingleEleMatchPhi -> push_back ( singleElectronTrigRef -> phi());
      } else { 
	HLTSingleEleMatched  -> push_back ( false ) ;
	HLTSingleEleMatchPt  -> push_back ( -999. );
	HLTSingleEleMatchEta -> push_back ( -999. );
	HLTSingleEleMatchPhi -> push_back ( -999. );
      }

      // Single electron (WP80)
      
      const pat::TriggerObjectRef singleElectronWP80TrigRef( matchHelper.triggerMatchObject( electrons, iElectron,  singleEleTriggerMatchWP80, iEvent, *triggerEvent ) );
      if ( singleElectronWP80TrigRef.isAvailable() && singleElectronWP80TrigRef.isNonnull() ) { 
	HLTSingleEleWP80Matched  -> push_back ( true ) ;
	HLTSingleEleWP80MatchPt  -> push_back ( singleElectronWP80TrigRef -> pt() );
	HLTSingleEleWP80MatchEta -> push_back ( singleElectronWP80TrigRef -> eta());
	HLTSingleEleWP80MatchPhi -> push_back ( singleElectronWP80TrigRef -> phi());
      } else { 
	HLTSingleEleWP80Matched  -> push_back ( false ) ;
	HLTSingleEleWP80MatchPt  -> push_back ( -999. );
	HLTSingleEleWP80MatchEta -> push_back ( -999. );
	HLTSingleEleWP80MatchPhi -> push_back ( -999. );
      }
      
      //------------------------------------------------------------------------
      // Relative isolation (not currently used in any analysis... remove?) 
      //------------------------------------------------------------------------
      
      double reliso = (it->trackIso() + it->ecalIso() + it->hcalIso())/it->pt();

      //------------------------------------------------------------------------
      // Conversion information
      //------------------------------------------------------------------------

      // Conversion (variables)
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

      // Conversion (fit) 
      bool matchesConv = false;
      if( hConversions.isValid() && bsHandle.isValid() ) 
	{
	  // See: https://twiki.cern.ch/twiki/bin/view/CMS/ConversionTools#Conversion_veto_for_electron_ID
	  matchesConv = ConversionTools::hasMatchedConversion(*it,hConversions,bsHandle->position());
	  
	}
      else 
	{
	  if( !bsHandle.isValid() )
	    edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get the product " << beamSpotInputTag;
	  if( !hConversions.isValid() )
	    edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get the product " << conversionsInputTag;
	}
      
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
        edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get the product " << vtxInputTag;
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
      charge                   -> push_back ( it->charge() );
      hoe                      -> push_back ( it->hadronicOverEm() );
      
      passIds                  -> push_back ( passId );

      // Supercluster kinematic variables

      scEta                    -> push_back( it->superCluster()->eta() );
      scPhi                    -> push_back( it->superCluster()->phi() );
      scPt                     -> push_back( it->superCluster()->energy()/cosh(it->superCluster()->eta()) );
      scRawEnergy              -> push_back( it->superCluster()->rawEnergy() );
      eSuperClusterOverP       -> push_back( it->eSuperClusterOverP() );
      
      // ID information
      passIds                  -> push_back( passId );

      // Does this electron overlap with a muon?
      overlaps                 -> push_back( ovrlps );
      
      // Number of Brems = number of basic clusters minus one
      numberOfBrems            -> push_back( it->numberOfBrems() );
      
      // Is this ECAL driven? Or PFlow?  
      hasEcalDrivenSeed        -> push_back( it->ecalDrivenSeed() );
      
      // ECAL eta/phi vs tracker eta/phi

      deltaPhiTrkSC            -> push_back ( it->deltaPhiSuperClusterTrackAtVtx() );
      deltaEtaTrkSC            -> push_back ( it->deltaEtaSuperClusterTrackAtVtx() );

      // Shower shape

      sigmaEtaEta              -> push_back ( it->sigmaEtaEta() );
      sigmaIEtaIEta            -> push_back ( it->sigmaIetaIeta() );
      classif                  -> push_back ( it->classification() );
      e1x5overe5x5             -> push_back ( (it->e5x5()>0) ? (it->e1x5()/it->e5x5()) : 0 );
      e2x5overe5x5             -> push_back ( (it->e5x5()>0) ? (it->e2x5Max()/it->e5x5()) : 0 );

      // Isolation variables: PAT
      
      trkIsoPAT                -> push_back( it->trackIso() );
      ecalIsoPAT               -> push_back( it->ecalIso() );
      hcalIsoPAT               -> push_back( it->hcalIso() );
      relIsoPAT                -> push_back( reliso );
      passIsoPAT               -> push_back( (reliso<electronIso) ? 1 : 0 ); 

      // Isolation variables: particle flow

      pfParticleIso            -> push_back ( it->particleIso       () );
      pfChargedHadronIso       -> push_back ( it->chargedHadronIso  () );
      pfNeutralHadronIso       -> push_back ( it->neutralHadronIso  () );
      pfPhotonIso              -> push_back ( it->photonIso         () );
      pfPUChargedHadronIso     -> push_back ( it->puChargedHadronIso() );
      
      // Isolation variables: DR 0.3

      ecalIsoDR03              -> push_back ( it->dr03EcalRecHitSumEt() );
      hcalIsoDR03              -> push_back ( it->dr03HcalTowerSumEt() );
      hcalIsoDR03FullCone      -> push_back ( it->dr03HcalTowerSumEt() +
					      ( it->hadronicOverEm() 
						* it->superCluster()->energy() 
						/ cosh(it->superCluster()->eta())));
      hcalIsoD1DR03            -> push_back( it->dr03HcalDepth1TowerSumEt() );
      hcalIsoD2DR03            -> push_back( it->dr03HcalDepth2TowerSumEt() );
      trkIsoDR03               -> push_back( it->dr03TkSumPt() );
      
      // Conversion variables
      
      missingHits              -> push_back ( it->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() );
      missingHitsEG            -> push_back ( it->gsfTrack()->trackerExpectedHitsInner().numberOfHits()     );
      
      dist_vec                 -> push_back ( dist );
      dCotTheta                -> push_back ( dcot );
      hasMatchedConvPhot       -> push_back ( matchesConv );
      fbrem                    -> push_back ( it->fbrem() );
      
      // Vertex association variables
      
      vtxIndex                 -> push_back( vtxIndex_  );
      vtxDistXY                -> push_back( vtxDistXY_ );
      vtxDistZ                 -> push_back( vtxDistZ_  );
      vtx0DistXY               -> push_back( vtx0DistXY_);
      vtx0DistZ                -> push_back( vtx0DistZ_ );
      primaryVertexDXY         -> push_back( fabs( it->dB() ) );      
      primaryVertexDXYError    -> push_back( fabs( it->edB() ) );
      beamspotDXY              -> push_back( fabs( it->dB (pat::Electron::BS2D) ) );
      beamspotDXYError         -> push_back( fabs( it->edB(pat::Electron::BS2D) ) );
			       
      // Track information     
      			       
      trackVx                  -> push_back( it->gsfTrack()->vx() );
      trackVy                  -> push_back( it->gsfTrack()->vy() );
      trackVz                  -> push_back( it->gsfTrack()->vz() );
      trackPt                  -> push_back( it->gsfTrack()->pt() );
      trackValidFractionOfHits -> push_back( it->gsfTrack()->validFraction() );


      ++iElectron;
    }
  } else {
    edm::LogError("RootTupleMakerV2_ElectronsError") << "Error! Can't get the product " << inputTag;
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
  iEvent.put( charge                  , prefix + "Charge"                   + suffix );
  iEvent.put( hoe                     , prefix + "HoE"                      + suffix );

  // Supercluster kinematic variables				        

  iEvent.put( eSuperClusterOverP      , prefix + "ESuperClusterOverP"       + suffix );
  iEvent.put( scEta                   , prefix + "SCEta"                    + suffix );
  iEvent.put( scPhi                   , prefix + "SCPhi"                    + suffix );
  iEvent.put( scPt                    , prefix + "SCPt"                     + suffix );
  iEvent.put( scRawEnergy             , prefix + "SCRawEnergy"              + suffix );

  // ID information 
  iEvent.put( passIds                 , prefix + "PassId"                   + suffix );
  
  // Does this electron overlap with a muon?			        
  iEvent.put( overlaps                , prefix + "Overlaps"                 + suffix );

  // Number of Brems = number of basic clusters minus one	        
  iEvent.put( numberOfBrems           , prefix + "NumberOfBrems"            + suffix );
  
  // Is this ECAL driven? Or PFlow?  				        
  iEvent.put( hasEcalDrivenSeed       , prefix + "HasEcalDrivenSeed"        + suffix );

  // ECAL eta/phi vs tracker eta/phi				        

  iEvent.put( deltaPhiTrkSC           , prefix + "DeltaPhiTrkSC"            + suffix );
  iEvent.put( deltaEtaTrkSC           , prefix + "DeltaEtaTrkSC"            + suffix );

  // Shower shape						        

  iEvent.put( sigmaEtaEta             , prefix + "SigmaEtaEta"              + suffix );
  iEvent.put( sigmaIEtaIEta           , prefix + "SigmaIEtaIEta"            + suffix );
  iEvent.put( e1x5overe5x5            , prefix + "E1x5OverE5x5"             + suffix );
  iEvent.put( e2x5overe5x5            , prefix + "E2x5OverE5x5"             + suffix );
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

  iEvent.put( pfParticleIso           , prefix + "PFParticleIso"            + suffix );
  iEvent.put( pfChargedHadronIso      , prefix + "PFChargedHadronIso"       + suffix );
  iEvent.put( pfNeutralHadronIso      , prefix + "PFNeutralHadronIso"       + suffix );
  iEvent.put( pfPhotonIso             , prefix + "PFPhotonIso"              + suffix );
  iEvent.put( pfPUChargedHadronIso    , prefix + "PFPUChargedHadronIso"     + suffix );

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

  // Trigger matching: Single electron (WP80)

  iEvent.put( HLTSingleEleWP80Matched , prefix + "HLTSingleEleWP80Matched"  + suffix );
  iEvent.put( HLTSingleEleWP80MatchPt , prefix + "HLTSingleEleWP80MatchPt"  + suffix );
  iEvent.put( HLTSingleEleWP80MatchEta, prefix + "HLTSingleEleWP80MatchEta" + suffix );
  iEvent.put( HLTSingleEleWP80MatchPhi, prefix + "HLTSingleEleWP80MatchPhi" + suffix );
}
