#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_PhysicsDSTStream2011.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

RootTupleMakerV2_PhysicsDSTStream2011::RootTupleMakerV2_PhysicsDSTStream2011(const edm::ParameterSet& iConfig) :
  hltInputTag      (iConfig.getParameter<edm::InputTag>("HLTInputTag")),
  hltInputToken_      (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLTInputTag"))),
  inputTokenHLTPFJets_ (consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("InputTagHLTPFJets"))),
  inputTokenHLTCaloJetsRaw_ (consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("InputTagHLTCaloJetsRaw"))),
  inputTokenHLTCaloJetsCorr_ (consumes<reco::CaloJetCollection>(iConfig.getParameter<edm::InputTag>("InputTagHLTCaloJetsCorr"))),
  inputTokenHLTPixelVertices_ (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("InputTagHLTPixelVertices"))),
  suffix  (iConfig.getParameter<std::string>  ("Suffix")),
  prefixHLTPFJets  (iConfig.getParameter<std::string>  ("PrefixHLTPFJets")),
  prefixHLTCaloJetsRaw  (iConfig.getParameter<std::string>  ("PrefixHLTCaloJetsRaw")),
  prefixHLTCaloJetsCorr  (iConfig.getParameter<std::string>  ("PrefixHLTCaloJetsCorr")),
  prefixHLTPixelVertices   (iConfig.getParameter<std::string>  ("PrefixHLTPixelVertices")),
  minPtHLTPFJets   (iConfig.getParameter<double> ("MinPtHLTPFJets")),
  maxEtaHLTPFJets  (iConfig.getParameter<double> ("MaxEtaHLTPFJets")),
  minPtHLTCaloJetsRaw   (iConfig.getParameter<double> ("MinPtHLTCaloJetsRaw")),
  maxEtaHLTCaloJetsRaw  (iConfig.getParameter<double> ("MaxEtaHLTCaloJetsRaw")),
  minPtHLTCaloJetsCorr   (iConfig.getParameter<double> ("MinPtHLTCaloJetsCorr")),
  maxEtaHLTCaloJetsCorr  (iConfig.getParameter<double> ("MaxEtaHLTCaloJetsCorr")),
  storeEventInfo   (iConfig.getParameter<bool> ("StoreEventInfo"))
{
  produces <unsigned int> ( "run"   );
  produces <unsigned int> ( "event" );
  produces <unsigned int> ( "bunch" );
  produces <unsigned int> ( "ls"    );
  produces <unsigned int> ( "orbit" );
  produces <double>       ( "time" );
  produces <bool>         ( "isData" );

  produces <std::string>               ("HLTKey");
  produces <std::vector<std::string> > ("HLTTriggerNames"      );
  produces <std::vector<bool > >       ("HLTTriggerDecisions"  );
  //   produces <std::vector<int> >         ("HLTTriggerPrescales"  );

  produces <std::vector<double> > ( prefixHLTPFJets + "Eta" + suffix );
  produces <std::vector<double> > ( prefixHLTPFJets + "Phi" + suffix );
  produces <std::vector<double> > ( prefixHLTPFJets + "Pt" + suffix );
  produces <std::vector<double> > ( prefixHLTPFJets + "Energy" + suffix );
  produces <std::vector<double> > ( prefixHLTPFJets + "ChargedEmEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefixHLTPFJets + "ChargedHadronEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefixHLTPFJets + "ChargedMuEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefixHLTPFJets + "ElectronEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefixHLTPFJets + "MuonEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefixHLTPFJets + "NeutralEmEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefixHLTPFJets + "NeutralHadronEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefixHLTPFJets + "PhotonEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefixHLTPFJets + "HFHadronEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefixHLTPFJets + "HFEMEnergyFraction"  + suffix );
  produces <std::vector<int> >    ( prefixHLTPFJets + "ChargedHadronMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefixHLTPFJets + "ChargedMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefixHLTPFJets + "ElectronMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefixHLTPFJets + "MuonMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefixHLTPFJets + "NeutralHadronMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefixHLTPFJets + "NeutralMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefixHLTPFJets + "PhotonMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefixHLTPFJets + "HFHadronMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefixHLTPFJets + "HFEMMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefixHLTPFJets + "NConstituents"  + suffix );
  produces <std::vector<bool> >   ( prefixHLTPFJets + "PassLooseID" + suffix);
  produces <std::vector<bool> >   ( prefixHLTPFJets + "PassTightID" + suffix);

  produces <std::vector<double> > ( prefixHLTCaloJetsRaw + "Eta" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsRaw + "Phi" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsRaw + "Pt" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsRaw + "Energy" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsRaw + "MaxEInEmTowers" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsRaw + "MaxEInHadTowers" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsRaw + "HadEnergyInHO" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsRaw + "HadEnergyInHB" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsRaw + "HadEnergyInHF" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsRaw + "HadEnergyInHE" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsRaw + "EmEnergyInEB" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsRaw + "EmEnergyInEE" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsRaw + "EmEnergyInHF" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsRaw + "EnergyFractionHadronic" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsRaw + "EnergyFractionEm" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsRaw + "TowersArea" + suffix );

  produces <std::vector<double> > ( prefixHLTCaloJetsCorr + "Eta" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsCorr + "Phi" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsCorr + "Pt" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsCorr + "Energy" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsCorr + "MaxEInEmTowers" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsCorr + "MaxEInHadTowers" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsCorr + "HadEnergyInHO" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsCorr + "HadEnergyInHB" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsCorr + "HadEnergyInHF" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsCorr + "HadEnergyInHE" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsCorr + "EmEnergyInEB" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsCorr + "EmEnergyInEE" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsCorr + "EmEnergyInHF" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsCorr + "EnergyFractionHadronic" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsCorr + "EnergyFractionEm" + suffix );
  produces <std::vector<double> > ( prefixHLTCaloJetsCorr + "TowersArea" + suffix );

  produces <std::vector<double> > ( prefixHLTPixelVertices + "Chi2" + suffix );
  produces <std::vector<double> > ( prefixHLTPixelVertices + "NormChi2" + suffix );
  produces <std::vector<double> > ( prefixHLTPixelVertices + "Ndof" + suffix );
  produces <std::vector<double> > ( prefixHLTPixelVertices + "XCoord" + suffix );
  produces <std::vector<double> > ( prefixHLTPixelVertices + "XError" + suffix );
  produces <std::vector<double> > ( prefixHLTPixelVertices + "YCoord" + suffix );
  produces <std::vector<double> > ( prefixHLTPixelVertices + "YError" + suffix );
  produces <std::vector<double> > ( prefixHLTPixelVertices + "ZCoord" + suffix );
  produces <std::vector<double> > ( prefixHLTPixelVertices + "ZError" + suffix );
  produces <std::vector<bool> > ( prefixHLTPixelVertices + "IsValid" + suffix );

}

void RootTupleMakerV2_PhysicsDSTStream2011::
printNames(const std::vector<std::string>& names) {
  for (unsigned int i = 0; i < names.size(); ++i)
    edm::LogProblem( "RootTupleMakerV2_PhysicsDSTStream2011_TriggerProblem" ) << "  " << names[i] << std::endl;
}

void RootTupleMakerV2_PhysicsDSTStream2011::
beginRun(edm::Run& iRun, const edm::EventSetup& iSetup) {

  bool changed = true;
  if (hltConfig.init(iRun, iSetup, hltInputTag.process(), changed)) {
    // if init returns TRUE, initialisation has succeeded!
    edm::LogInfo("RootTupleMakerV2_PhysicsDSTStream2011_TriggerInfo") << "HLT config with process name " << hltInputTag.process() << " successfully extracted";
  } else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    edm::LogError("RootTupleMakerV2_PhysicsDSTStream2011_TriggerError") << "Error! HLT config extraction with process name " << hltInputTag.process() << " failed";
    // In this case, all access methods will return empty values!
  }

}

void RootTupleMakerV2_PhysicsDSTStream2011::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::unique_ptr<unsigned int >  run   ( new unsigned int(iEvent.id().run()        ) );
  std::unique_ptr<unsigned int >  event ( new unsigned int(iEvent.id().event()      ) );
  std::unique_ptr<unsigned int >  ls    ( new unsigned int(iEvent.luminosityBlock() ) );
  double sec  = iEvent.time().value() >> 32 ;
  double usec = 0xFFFFFFFF & iEvent.time().value();
  double conv = 1e6;
  std::unique_ptr<unsigned int >  bunch ( new unsigned int(iEvent.bunchCrossing()   ) );
  std::unique_ptr<unsigned int >  orbit ( new unsigned int(iEvent.orbitNumber()     ) );
  std::unique_ptr<double >        time  ( new double(sec+usec/conv));
  std::unique_ptr<bool >          isdata  ( new bool(iEvent.isRealData()));

  std::unique_ptr<std::vector < std::string > > v_hlt_names             (new std::vector<std::string>  ());
  std::unique_ptr<std::vector < bool > >        v_hlt_decisions         (new std::vector<bool>         ());
  //   std::unique_ptr<std::vector < int > >         v_hlt_prescales         (new std::vector<int>          ());

  std::unique_ptr<std::vector<double> >  etaHLTPFJet     ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  phiHLTPFJet     ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  ptHLTPFJet      ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  energyHLTPFJet  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  chargedEmEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::unique_ptr<std::vector<double> >  chargedHadronEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::unique_ptr<std::vector<double> >  chargedMuEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::unique_ptr<std::vector<double> >  electronEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::unique_ptr<std::vector<double> >  muonEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::unique_ptr<std::vector<double> >  neutralEmEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::unique_ptr<std::vector<double> >  neutralHadronEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::unique_ptr<std::vector<double> >  photonEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::unique_ptr<std::vector<double> >  hfHadronEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::unique_ptr<std::vector<double> >  hfEMEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::unique_ptr<std::vector<int> >     chargedHadronMultiplicityHLTPFJet  ( new std::vector<int>()  ) ;
  std::unique_ptr<std::vector<int> >     chargedMultiplicityHLTPFJet  ( new std::vector<int>()  ) ;
  std::unique_ptr<std::vector<int> >     electronMultiplicityHLTPFJet  ( new std::vector<int>()  ) ;
  std::unique_ptr<std::vector<int> >     muonMultiplicityHLTPFJet  ( new std::vector<int>()  ) ;
  std::unique_ptr<std::vector<int> >     neutralHadronMultiplicityHLTPFJet  ( new std::vector<int>()  ) ;
  std::unique_ptr<std::vector<int> >     neutralMultiplicityHLTPFJet  ( new std::vector<int>()  ) ;
  std::unique_ptr<std::vector<int> >     photonMultiplicityHLTPFJet  ( new std::vector<int>()  ) ;
  std::unique_ptr<std::vector<int> >     hfHadronMultiplicityHLTPFJet ( new std::vector<int>()  ) ;
  std::unique_ptr<std::vector<int> >     hfEMMultiplicityHLTPFJet ( new std::vector<int>()  ) ;
  std::unique_ptr<std::vector<int> >     nConstituents  ( new std::vector<int>()  ) ;
  std::unique_ptr<std::vector<bool> >    passLooseID  ( new std::vector<bool>()  );
  std::unique_ptr<std::vector<bool> >    passTightID  ( new std::vector<bool>()  );

  std::unique_ptr<std::vector<double> >  etaHLTCaloJetRaw     ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  phiHLTCaloJetRaw     ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  ptHLTCaloJetRaw      ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  energyHLTCaloJetRaw  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  maxEInEmTowersHLTCaloJetRaw  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  maxEInHadTowersHLTCaloJetRaw  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  hadEnergyInHOHLTCaloJetRaw  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  hadEnergyInHBHLTCaloJetRaw  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  hadEnergyInHFHLTCaloJetRaw  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  hadEnergyInHEHLTCaloJetRaw  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  emEnergyInEBHLTCaloJetRaw  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  emEnergyInEEHLTCaloJetRaw  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  emEnergyInHFHLTCaloJetRaw  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  energyFractionHadronicHLTCaloJetRaw  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  energyFractionEmHLTCaloJetRaw  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  towersAreaHLTCaloJetRaw  ( new std::vector<double>()  );

  std::unique_ptr<std::vector<double> >  etaHLTCaloJetCorr     ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  phiHLTCaloJetCorr     ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  ptHLTCaloJetCorr      ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  energyHLTCaloJetCorr  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  maxEInEmTowersHLTCaloJetCorr  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  maxEInHadTowersHLTCaloJetCorr  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  hadEnergyInHOHLTCaloJetCorr  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  hadEnergyInHBHLTCaloJetCorr  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  hadEnergyInHFHLTCaloJetCorr  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  hadEnergyInHEHLTCaloJetCorr  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  emEnergyInEBHLTCaloJetCorr  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  emEnergyInEEHLTCaloJetCorr  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  emEnergyInHFHLTCaloJetCorr  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  energyFractionHadronicHLTCaloJetCorr  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  energyFractionEmHLTCaloJetCorr  ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  towersAreaHLTCaloJetCorr  ( new std::vector<double>()  );

  std::unique_ptr<std::vector<double> >  chi2HLTPixelVertices     ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  normChi2HLTPixelVertices     ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  ndofHLTPixelVertices     ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  XCoordHLTPixelVertices     ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  XErrorHLTPixelVertices     ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  YCoordHLTPixelVertices     ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  YErrorHLTPixelVertices     ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  ZCoordHLTPixelVertices     ( new std::vector<double>()  );
  std::unique_ptr<std::vector<double> >  ZErrorHLTPixelVertices     ( new std::vector<double>()  );
  std::unique_ptr<std::vector<bool> >  isValidHLTPixelVertices     ( new std::vector<bool>()  );

  //-----------------------------------------------------------------

  //---------------
  // HLT Info
  //---------------

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(hltInputToken_, triggerResults);

  if(triggerResults.isValid()) {
    edm::LogInfo("RootTupleMakerV2_PhysicsDSTStream2011_TriggerInfo") << "Successfully obtained " << triggerResults;//fixme what to put here?
    
    const edm::TriggerNames& names = iEvent.triggerNames(*triggerResults);
    
    for (int i = 0; i < (int) triggerResults->size() ; ++i) { 
      v_hlt_names->push_back ( names.triggerName(i) );
      v_hlt_decisions->push_back ( triggerResults->accept(i) );
      //v_hlt_prescales->push_back ( hltConfig.prescaleValue(iEvent,iSetup,names.triggerName(i)));	
    } 
  } else {
    edm::LogError("RootTupleMakerV2_PhysicsDSTStream2011_TriggerError") << "Error! Can't get the triggerResults";
  }

  //---------------
  // HLT PFJets
  //---------------

  edm::Handle<reco::PFJetCollection> pfjets;
  iEvent.getByToken(inputTokenHLTPFJets_, pfjets);
	
  if(pfjets.isValid())
    {
      edm::LogInfo("RootTupleMakerV2_PhysicsDSTStream2011Info") << "Total # HLTPFJets: " << pfjets->size();
      //std::cout << "Total # HLTPFJets: " << pfjets->size() << std::endl;

      for( reco::PFJetCollection::const_iterator it = pfjets->begin(); it != pfjets->end(); ++it )
	{
	  // exit from loop when you reach the required number of jets
	  if(it->pt() < minPtHLTPFJets || fabs(it->eta()) > maxEtaHLTPFJets )
	    continue;

	  etaHLTPFJet->push_back( it->eta() );
	  phiHLTPFJet->push_back( it->phi() );
	  ptHLTPFJet->push_back( it->pt() );
	  energyHLTPFJet->push_back( it->energy() );

	  int Nconstituents_var = it->chargedMultiplicity() + it->neutralMultiplicity();

 	  chargedEmEnergyFractionHLTPFJet->push_back( it->chargedEmEnergyFraction() );
 	  chargedHadronEnergyFractionHLTPFJet->push_back( it->chargedHadronEnergyFraction() );
	  chargedMuEnergyFractionHLTPFJet->push_back( it->chargedMuEnergyFraction() );
	  electronEnergyFractionHLTPFJet->push_back( it->electronEnergyFraction() );
	  muonEnergyFractionHLTPFJet->push_back( it->muonEnergyFraction() );
	  neutralEmEnergyFractionHLTPFJet->push_back( it->neutralEmEnergyFraction() );
	  neutralHadronEnergyFractionHLTPFJet->push_back( it->neutralHadronEnergyFraction() );
	  photonEnergyFractionHLTPFJet->push_back( it->photonEnergyFraction() );
	  hfHadronEnergyFractionHLTPFJet->push_back( it->HFHadronEnergyFraction() );
	  hfEMEnergyFractionHLTPFJet->push_back( it->HFEMEnergyFraction() );
	  chargedHadronMultiplicityHLTPFJet->push_back( it->chargedHadronMultiplicity() );
	  chargedMultiplicityHLTPFJet->push_back( it->chargedMultiplicity() );
	  electronMultiplicityHLTPFJet->push_back( it->electronMultiplicity() );
	  muonMultiplicityHLTPFJet->push_back( it->muonMultiplicity() );
	  neutralHadronMultiplicityHLTPFJet->push_back( it->neutralHadronMultiplicity() );
	  neutralMultiplicityHLTPFJet->push_back( it->neutralMultiplicity() );
	  photonMultiplicityHLTPFJet->push_back( it->photonMultiplicity() );
	  hfHadronMultiplicityHLTPFJet->push_back( it->HFHadronMultiplicity() );
	  hfEMMultiplicityHLTPFJet->push_back( it->HFEMMultiplicity() );
	  nConstituents->push_back( Nconstituents_var ); // check?

	  // PFJet ID ( https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Documentation ) 
	  // --> only defined for jet |eta| < 3

	  bool passLooseID_var = true;
	  bool passTightID_var = true;
	  
	  if( fabs( it->eta() )<3 
	      && 
	      (
	       it->neutralHadronEnergyFraction() >= 0.99
	       || it->neutralEmEnergyFraction() >= 0.99
	       || Nconstituents_var <= 1
	       )
	      )
	    {
	      passLooseID_var = false;
	    }

	  if( fabs( it->eta() )<3 
	      && 
	      (
	       it->neutralHadronEnergyFraction() >= 0.90
	       || it->neutralEmEnergyFraction() >= 0.90
	       || Nconstituents_var <= 1
	       )
	      )
	    {
	      passTightID_var = false;
	    }

	  if( fabs( it->eta() )<2.4
	      &&
	      (
	       it->chargedHadronEnergyFraction() <=0
	       || it->chargedMultiplicity() <=0
	       || it->chargedEmEnergyFraction() >=0.99
	       )
	      )
	    {
	      passLooseID_var = false;
	      passTightID_var = false;
	    }
	  
	  passLooseID->push_back( passLooseID_var );
	  passTightID->push_back( passTightID_var );	      
	  
	}
    }
  else
    {
      edm::LogInfo("RootTupleMakerV2_PhysicsDSTStream2011Info") << "Can't get the product " << pfjets;//fixme what to put here?
    }
  
  //---------------
  // HLT CaloJetsRaw
  //---------------

  edm::Handle<reco::CaloJetCollection> calojetsraw;
  iEvent.getByToken(inputTokenHLTCaloJetsRaw_, calojetsraw);
	
  if(calojetsraw.isValid())
    {
      edm::LogInfo("RootTupleMakerV2_PhysicsDSTStream2011Info") << "Total # HLTCaloJetsRaw: " << calojetsraw->size();
      //std::cout << "Total # HLTCaloJetsRaw: " << calojetsraw->size() << std::endl;

      for( reco::CaloJetCollection::const_iterator it = calojetsraw->begin(); it != calojetsraw->end(); ++it )
	{
	  // exit from loop when you reach the required number of jets
	  if(it->pt() < minPtHLTCaloJetsRaw || fabs(it->eta()) > maxEtaHLTCaloJetsRaw )
	    continue;

	  etaHLTCaloJetRaw->push_back( it->eta() );
	  phiHLTCaloJetRaw->push_back( it->phi() );
	  ptHLTCaloJetRaw->push_back( it->pt() );
	  energyHLTCaloJetRaw->push_back( it->energy() );

	  maxEInEmTowersHLTCaloJetRaw->push_back( it->maxEInEmTowers() );
	  maxEInHadTowersHLTCaloJetRaw->push_back( it->maxEInHadTowers() );
	  hadEnergyInHOHLTCaloJetRaw->push_back( it->hadEnergyInHO() );
	  hadEnergyInHBHLTCaloJetRaw->push_back( it->hadEnergyInHB() );
	  hadEnergyInHFHLTCaloJetRaw->push_back( it->hadEnergyInHF() );
	  hadEnergyInHEHLTCaloJetRaw->push_back( it->hadEnergyInHE() );
	  emEnergyInEBHLTCaloJetRaw->push_back( it->emEnergyInEB() );
	  emEnergyInEEHLTCaloJetRaw->push_back( it->emEnergyInEE() );
	  emEnergyInHFHLTCaloJetRaw->push_back( it->emEnergyInHF() );
	  energyFractionHadronicHLTCaloJetRaw->push_back( it->energyFractionHadronic() );
	  energyFractionEmHLTCaloJetRaw->push_back( it->emEnergyFraction() );
	  towersAreaHLTCaloJetRaw->push_back( it->towersArea() );
	}
    }
  else
    {
      edm::LogInfo("RootTupleMakerV2_PhysicsDSTStream2011Info") << "Can't get the product " << calojetsraw; //fixme what to put here?
    }


  //---------------
  // HLT CaloJetsCorr
  //---------------

  edm::Handle<reco::CaloJetCollection> calojetscorr;
  iEvent.getByToken(inputTokenHLTCaloJetsCorr_, calojetscorr);
	
  if(calojetscorr.isValid())
    {
      edm::LogInfo("RootTupleMakerV2_PhysicsDSTStream2011Info") << "Total # HLTCaloJetsCorr: " << calojetscorr->size();
      //std::cout << "Total # HLTCaloJetsCorr: " << calojetscorr->size() << std::endl;

      for( reco::CaloJetCollection::const_iterator it = calojetscorr->begin(); it != calojetscorr->end(); ++it )
	{
	  // exit from loop when you reach the required number of jets
	  if(it->pt() < minPtHLTCaloJetsCorr || fabs(it->eta()) > maxEtaHLTCaloJetsCorr )
	    continue;

	  etaHLTCaloJetCorr->push_back( it->eta() );
	  phiHLTCaloJetCorr->push_back( it->phi() );
	  ptHLTCaloJetCorr->push_back( it->pt() );
	  energyHLTCaloJetCorr->push_back( it->energy() );

	  maxEInEmTowersHLTCaloJetCorr->push_back( it->maxEInEmTowers() );
	  maxEInHadTowersHLTCaloJetCorr->push_back( it->maxEInHadTowers() );
	  hadEnergyInHOHLTCaloJetCorr->push_back( it->hadEnergyInHO() );
	  hadEnergyInHBHLTCaloJetCorr->push_back( it->hadEnergyInHB() );
	  hadEnergyInHFHLTCaloJetCorr->push_back( it->hadEnergyInHF() );
	  hadEnergyInHEHLTCaloJetCorr->push_back( it->hadEnergyInHE() );
	  emEnergyInEBHLTCaloJetCorr->push_back( it->emEnergyInEB() );
	  emEnergyInEEHLTCaloJetCorr->push_back( it->emEnergyInEE() );
	  emEnergyInHFHLTCaloJetCorr->push_back( it->emEnergyInHF() );
	  energyFractionHadronicHLTCaloJetCorr->push_back( it->energyFractionHadronic() );
	  energyFractionEmHLTCaloJetCorr->push_back( it->emEnergyFraction() );
	  towersAreaHLTCaloJetCorr->push_back( it->towersArea() );
	}
    }
  else
    {
      edm::LogInfo("RootTupleMakerV2_PhysicsDSTStream2011Info") << "Can't get the product " << calojetscorr; //fixme what to put here?
    }



  //---------------
  // HLT PixelVertices
  //---------------

  edm::Handle<reco::VertexCollection> pixelvertices;
  iEvent.getByToken(inputTokenHLTPixelVertices_, pixelvertices);
	
  if(pixelvertices.isValid())
    {
      edm::LogInfo("RootTupleMakerV2_PhysicsDSTStream2011Info") << "Total # HLTPixelVertices: " << pixelvertices->size();
      //std::cout << "Total # HLTPixelVertices: " << pixelvertices->size() << std::endl;

      for( reco::VertexCollection::const_iterator it = pixelvertices->begin(); it != pixelvertices->end(); ++it )
	{
	  chi2HLTPixelVertices->push_back( it->chi2() );
	  normChi2HLTPixelVertices->push_back( it->normalizedChi2() );
	  ndofHLTPixelVertices->push_back( it->ndof() );
	  XCoordHLTPixelVertices->push_back( it->x() );
	  XErrorHLTPixelVertices->push_back( it->xError() );
	  YCoordHLTPixelVertices->push_back( it->y() );
	  YErrorHLTPixelVertices->push_back( it->yError() );
	  ZCoordHLTPixelVertices->push_back( it->z() );
	  ZErrorHLTPixelVertices->push_back( it->zError() );
	  isValidHLTPixelVertices->push_back( it->isValid() );
	}
    }
  else
    {
      edm::LogInfo("RootTupleMakerV2_PhysicsDSTStream2011Info") << "Can't get the product " << pixelvertices; //fixme what to put here?
    }

       
  //-----------------------------------------------------------------
  // put vectors in the event

  // Event Info
  if(storeEventInfo)
    {
      iEvent.put(std::move( run),   "run"   );
      iEvent.put(std::move( event), "event" );
      iEvent.put(std::move( ls   ), "ls"    );
      iEvent.put(std::move( bunch), "bunch" );
      iEvent.put(std::move( orbit), "orbit" );
      iEvent.put(std::move( time),  "time"  );
      iEvent.put(std::move( isdata),"isData");
    }

  // HLT Info
  iEvent.put(std::move(std::unique_ptr<std::string>(new std::string(hltConfig.tableName())) ), "HLTKey");
  iEvent.put(std::move(v_hlt_names    ), "HLTTriggerNames"      ) ;
  iEvent.put(std::move(v_hlt_decisions), "HLTTriggerDecisions"  ) ;
  //   iEvent.put ( v_hlt_prescales  ), "HLTTriggerPrescales"  ) ;

  // HLT PFJets
  iEvent.put(std::move(etaHLTPFJet), prefixHLTPFJets + "Eta" + suffix );
  iEvent.put(std::move(phiHLTPFJet), prefixHLTPFJets + "Phi" + suffix );
  iEvent.put(std::move(ptHLTPFJet), prefixHLTPFJets + "Pt" + suffix );
  iEvent.put(std::move(energyHLTPFJet), prefixHLTPFJets + "Energy" + suffix );

  iEvent.put(std::move(chargedEmEnergyFractionHLTPFJet),  prefixHLTPFJets + "ChargedEmEnergyFraction"  + suffix );
  iEvent.put(std::move(chargedHadronEnergyFractionHLTPFJet),  prefixHLTPFJets + "ChargedHadronEnergyFraction"  + suffix );
  iEvent.put(std::move(chargedMuEnergyFractionHLTPFJet),  prefixHLTPFJets + "ChargedMuEnergyFraction"  + suffix );
  iEvent.put(std::move(electronEnergyFractionHLTPFJet),  prefixHLTPFJets + "ElectronEnergyFraction"  + suffix );
  iEvent.put(std::move(muonEnergyFractionHLTPFJet),  prefixHLTPFJets + "MuonEnergyFraction"  + suffix );
  iEvent.put(std::move(neutralEmEnergyFractionHLTPFJet),  prefixHLTPFJets + "NeutralEmEnergyFraction"  + suffix );
  iEvent.put(std::move(neutralHadronEnergyFractionHLTPFJet),  prefixHLTPFJets + "NeutralHadronEnergyFraction"  + suffix );
  iEvent.put(std::move(photonEnergyFractionHLTPFJet),  prefixHLTPFJets + "PhotonEnergyFraction"  + suffix );
  iEvent.put(std::move(hfHadronEnergyFractionHLTPFJet),  prefixHLTPFJets + "HFHadronEnergyFraction"  + suffix );
  iEvent.put(std::move(hfEMEnergyFractionHLTPFJet),  prefixHLTPFJets + "HFEMEnergyFraction"  + suffix );
  iEvent.put(std::move(chargedHadronMultiplicityHLTPFJet),  prefixHLTPFJets + "ChargedHadronMultiplicity"  + suffix );
  iEvent.put(std::move(chargedMultiplicityHLTPFJet),  prefixHLTPFJets + "ChargedMultiplicity"  + suffix );
  iEvent.put(std::move(electronMultiplicityHLTPFJet),  prefixHLTPFJets + "ElectronMultiplicity"  + suffix );
  iEvent.put(std::move(muonMultiplicityHLTPFJet),  prefixHLTPFJets + "MuonMultiplicity"  + suffix );
  iEvent.put(std::move(neutralHadronMultiplicityHLTPFJet),  prefixHLTPFJets + "NeutralHadronMultiplicity"  + suffix );
  iEvent.put(std::move(neutralMultiplicityHLTPFJet),  prefixHLTPFJets + "NeutralMultiplicity"  + suffix );
  iEvent.put(std::move(photonMultiplicityHLTPFJet),  prefixHLTPFJets + "PhotonMultiplicity"  + suffix );
  iEvent.put(std::move(hfHadronMultiplicityHLTPFJet),  prefixHLTPFJets + "HFHadronMultiplicity"  + suffix );
  iEvent.put(std::move(hfEMMultiplicityHLTPFJet),  prefixHLTPFJets + "HFEMMultiplicity"  + suffix );
  iEvent.put(std::move(nConstituents),  prefixHLTPFJets + "NConstituents"  + suffix );
  iEvent.put(std::move(passLooseID), prefixHLTPFJets + "PassLooseID" + suffix);
  iEvent.put(std::move(passTightID), prefixHLTPFJets + "PassTightID" + suffix);

  // HLT CaloJetsRaw
  iEvent.put(std::move(etaHLTCaloJetRaw), prefixHLTCaloJetsRaw + "Eta" + suffix );
  iEvent.put(std::move(phiHLTCaloJetRaw), prefixHLTCaloJetsRaw + "Phi" + suffix );
  iEvent.put(std::move(ptHLTCaloJetRaw), prefixHLTCaloJetsRaw + "Pt" + suffix );
  iEvent.put(std::move(energyHLTCaloJetRaw), prefixHLTCaloJetsRaw + "Energy" + suffix );

  iEvent.put(std::move(maxEInEmTowersHLTCaloJetRaw), prefixHLTCaloJetsRaw + "MaxEInEmTowers" + suffix );
  iEvent.put(std::move(maxEInHadTowersHLTCaloJetRaw), prefixHLTCaloJetsRaw + "MaxEInHadTowers" + suffix );
  iEvent.put(std::move(hadEnergyInHOHLTCaloJetRaw), prefixHLTCaloJetsRaw + "HadEnergyInHO" + suffix );
  iEvent.put(std::move(hadEnergyInHBHLTCaloJetRaw), prefixHLTCaloJetsRaw + "HadEnergyInHB" + suffix );
  iEvent.put(std::move(hadEnergyInHFHLTCaloJetRaw), prefixHLTCaloJetsRaw + "HadEnergyInHF" + suffix );
  iEvent.put(std::move(hadEnergyInHEHLTCaloJetRaw), prefixHLTCaloJetsRaw + "HadEnergyInHE" + suffix );
  iEvent.put(std::move(emEnergyInEBHLTCaloJetRaw), prefixHLTCaloJetsRaw + "EmEnergyInEB" + suffix );
  iEvent.put(std::move(emEnergyInEEHLTCaloJetRaw), prefixHLTCaloJetsRaw + "EmEnergyInEE" + suffix );
  iEvent.put(std::move(emEnergyInHFHLTCaloJetRaw), prefixHLTCaloJetsRaw + "EmEnergyInHF" + suffix );
  iEvent.put(std::move(energyFractionHadronicHLTCaloJetRaw), prefixHLTCaloJetsRaw + "EnergyFractionHadronic" + suffix );
  iEvent.put(std::move(energyFractionEmHLTCaloJetRaw), prefixHLTCaloJetsRaw + "EnergyFractionEm" + suffix );
  iEvent.put(std::move(towersAreaHLTCaloJetRaw), prefixHLTCaloJetsRaw + "TowersArea" + suffix );

  // HLT CaloJetsCorr
  iEvent.put(std::move(etaHLTCaloJetCorr), prefixHLTCaloJetsCorr + "Eta" + suffix );
  iEvent.put(std::move(phiHLTCaloJetCorr), prefixHLTCaloJetsCorr + "Phi" + suffix );
  iEvent.put(std::move(ptHLTCaloJetCorr), prefixHLTCaloJetsCorr + "Pt" + suffix );
  iEvent.put(std::move(energyHLTCaloJetCorr), prefixHLTCaloJetsCorr + "Energy" + suffix );

  iEvent.put(std::move(maxEInEmTowersHLTCaloJetCorr), prefixHLTCaloJetsCorr + "MaxEInEmTowers" + suffix );
  iEvent.put(std::move(maxEInHadTowersHLTCaloJetCorr), prefixHLTCaloJetsCorr + "MaxEInHadTowers" + suffix );
  iEvent.put(std::move(hadEnergyInHOHLTCaloJetCorr), prefixHLTCaloJetsCorr + "HadEnergyInHO" + suffix );
  iEvent.put(std::move(hadEnergyInHBHLTCaloJetCorr), prefixHLTCaloJetsCorr + "HadEnergyInHB" + suffix );
  iEvent.put(std::move(hadEnergyInHFHLTCaloJetCorr), prefixHLTCaloJetsCorr + "HadEnergyInHF" + suffix );
  iEvent.put(std::move(hadEnergyInHEHLTCaloJetCorr), prefixHLTCaloJetsCorr + "HadEnergyInHE" + suffix );
  iEvent.put(std::move(emEnergyInEBHLTCaloJetCorr), prefixHLTCaloJetsCorr + "EmEnergyInEB" + suffix );
  iEvent.put(std::move(emEnergyInEEHLTCaloJetCorr), prefixHLTCaloJetsCorr + "EmEnergyInEE" + suffix );
  iEvent.put(std::move(emEnergyInHFHLTCaloJetCorr), prefixHLTCaloJetsCorr + "EmEnergyInHF" + suffix );
  iEvent.put(std::move(energyFractionHadronicHLTCaloJetCorr), prefixHLTCaloJetsCorr + "EnergyFractionHadronic" + suffix );
  iEvent.put(std::move(energyFractionEmHLTCaloJetCorr), prefixHLTCaloJetsCorr + "EnergyFractionEm" + suffix );
  iEvent.put(std::move(towersAreaHLTCaloJetCorr), prefixHLTCaloJetsCorr + "TowersArea" + suffix );

  // HLT PixelVertices
  iEvent.put(std::move(chi2HLTPixelVertices ), prefixHLTPixelVertices + "Chi2" + suffix );
  iEvent.put(std::move(normChi2HLTPixelVertices ), prefixHLTPixelVertices + "NormChi2" + suffix );
  iEvent.put(std::move(ndofHLTPixelVertices ), prefixHLTPixelVertices + "Ndof" + suffix );
  iEvent.put(std::move(XCoordHLTPixelVertices ), prefixHLTPixelVertices + "XCoord" + suffix );
  iEvent.put(std::move(XErrorHLTPixelVertices ), prefixHLTPixelVertices + "XError" + suffix );
  iEvent.put(std::move(YCoordHLTPixelVertices ), prefixHLTPixelVertices + "YCoord" + suffix );
  iEvent.put(std::move(YErrorHLTPixelVertices ), prefixHLTPixelVertices + "YError" + suffix );
  iEvent.put(std::move(ZCoordHLTPixelVertices ), prefixHLTPixelVertices + "ZCoord" + suffix );
  iEvent.put(std::move(ZErrorHLTPixelVertices ), prefixHLTPixelVertices + "ZError" + suffix );
  iEvent.put(std::move(isValidHLTPixelVertices ), prefixHLTPixelVertices + "IsValid" + suffix );

}
