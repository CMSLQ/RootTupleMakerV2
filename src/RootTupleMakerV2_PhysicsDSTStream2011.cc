#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_PhysicsDSTStream2011.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/TriggerResults.h"

RootTupleMakerV2_PhysicsDSTStream2011::RootTupleMakerV2_PhysicsDSTStream2011(const edm::ParameterSet& iConfig) :
suffix  (iConfig.getParameter<std::string>  ("Suffix")),
storeEventInfo   (iConfig.getParameter<bool> ("StoreEventInfo")), 
hltInputTag      (iConfig.getParameter<edm::InputTag>("HLTInputTag")),
inputTagHLTPFJets(iConfig.getParameter<edm::InputTag>("InputTagHLTPFJets")),
prefixHLTPFJets  (iConfig.getParameter<std::string>  ("PrefixHLTPFJets")),
minPtHLTPFJets   (iConfig.getParameter<double> ("MinPtHLTPFJets")),
maxEtaHLTPFJets  (iConfig.getParameter<double> ("MaxEtaHLTPFJets")),
inputTagHLTCaloJetsRaw(iConfig.getParameter<edm::InputTag>("InputTagHLTCaloJetsRaw")),
prefixHLTCaloJetsRaw  (iConfig.getParameter<std::string>  ("PrefixHLTCaloJetsRaw")),
minPtHLTCaloJetsRaw   (iConfig.getParameter<double> ("MinPtHLTCaloJetsRaw")),
maxEtaHLTCaloJetsRaw  (iConfig.getParameter<double> ("MaxEtaHLTCaloJetsRaw")),
inputTagHLTCaloJetsCorr(iConfig.getParameter<edm::InputTag>("InputTagHLTCaloJetsCorr")),
prefixHLTCaloJetsCorr  (iConfig.getParameter<std::string>  ("PrefixHLTCaloJetsCorr")),
minPtHLTCaloJetsCorr   (iConfig.getParameter<double> ("MinPtHLTCaloJetsCorr")),
maxEtaHLTCaloJetsCorr  (iConfig.getParameter<double> ("MaxEtaHLTCaloJetsCorr")),
inputTagHLTPixelVertices (iConfig.getParameter<edm::InputTag>("InputTagHLTPixelVertices")),
prefixHLTPixelVertices   (iConfig.getParameter<std::string>  ("PrefixHLTPixelVertices"))
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
  std::auto_ptr<unsigned int >  run   ( new unsigned int(iEvent.id().run()        ) );
  std::auto_ptr<unsigned int >  event ( new unsigned int(iEvent.id().event()      ) );
  std::auto_ptr<unsigned int >  ls    ( new unsigned int(iEvent.luminosityBlock() ) );
  double sec  = iEvent.time().value() >> 32 ;
  double usec = 0xFFFFFFFF & iEvent.time().value();
  double conv = 1e6;
  std::auto_ptr<unsigned int >  bunch ( new unsigned int(iEvent.bunchCrossing()   ) );
  std::auto_ptr<unsigned int >  orbit ( new unsigned int(iEvent.orbitNumber()     ) );
  std::auto_ptr<double >        time  ( new double(sec+usec/conv));
  std::auto_ptr<bool >          isdata  ( new bool(iEvent.isRealData()));

  std::auto_ptr<std::vector < std::string > > v_hlt_names             (new std::vector<std::string>  ());
  std::auto_ptr<std::vector < bool > >        v_hlt_decisions         (new std::vector<bool>         ());
  //   std::auto_ptr<std::vector < int > >         v_hlt_prescales         (new std::vector<int>          ());

  std::auto_ptr<std::vector<double> >  etaHLTPFJet     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phiHLTPFJet     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ptHLTPFJet      ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energyHLTPFJet  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  chargedEmEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  chargedHadronEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  chargedMuEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  electronEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  muonEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  neutralEmEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  neutralHadronEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  photonEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  hfHadronEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  hfEMEnergyFractionHLTPFJet  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<int> >     chargedHadronMultiplicityHLTPFJet  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     chargedMultiplicityHLTPFJet  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     electronMultiplicityHLTPFJet  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     muonMultiplicityHLTPFJet  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     neutralHadronMultiplicityHLTPFJet  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     neutralMultiplicityHLTPFJet  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     photonMultiplicityHLTPFJet  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     hfHadronMultiplicityHLTPFJet ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     hfEMMultiplicityHLTPFJet ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     nConstituents  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<bool> >    passLooseID  ( new std::vector<bool>()  );
  std::auto_ptr<std::vector<bool> >    passTightID  ( new std::vector<bool>()  );

  std::auto_ptr<std::vector<double> >  etaHLTCaloJetRaw     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phiHLTCaloJetRaw     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ptHLTCaloJetRaw      ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energyHLTCaloJetRaw  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  maxEInEmTowersHLTCaloJetRaw  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  maxEInHadTowersHLTCaloJetRaw  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hadEnergyInHOHLTCaloJetRaw  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hadEnergyInHBHLTCaloJetRaw  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hadEnergyInHFHLTCaloJetRaw  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hadEnergyInHEHLTCaloJetRaw  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  emEnergyInEBHLTCaloJetRaw  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  emEnergyInEEHLTCaloJetRaw  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  emEnergyInHFHLTCaloJetRaw  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energyFractionHadronicHLTCaloJetRaw  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energyFractionEmHLTCaloJetRaw  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  towersAreaHLTCaloJetRaw  ( new std::vector<double>()  );

  std::auto_ptr<std::vector<double> >  etaHLTCaloJetCorr     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phiHLTCaloJetCorr     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ptHLTCaloJetCorr      ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energyHLTCaloJetCorr  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  maxEInEmTowersHLTCaloJetCorr  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  maxEInHadTowersHLTCaloJetCorr  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hadEnergyInHOHLTCaloJetCorr  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hadEnergyInHBHLTCaloJetCorr  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hadEnergyInHFHLTCaloJetCorr  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hadEnergyInHEHLTCaloJetCorr  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  emEnergyInEBHLTCaloJetCorr  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  emEnergyInEEHLTCaloJetCorr  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  emEnergyInHFHLTCaloJetCorr  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energyFractionHadronicHLTCaloJetCorr  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energyFractionEmHLTCaloJetCorr  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  towersAreaHLTCaloJetCorr  ( new std::vector<double>()  );

  std::auto_ptr<std::vector<double> >  chi2HLTPixelVertices     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  normChi2HLTPixelVertices     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ndofHLTPixelVertices     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  XCoordHLTPixelVertices     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  XErrorHLTPixelVertices     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  YCoordHLTPixelVertices     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  YErrorHLTPixelVertices     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ZCoordHLTPixelVertices     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ZErrorHLTPixelVertices     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<bool> >  isValidHLTPixelVertices     ( new std::vector<bool>()  );

  //-----------------------------------------------------------------

  //---------------
  // HLT Info
  //---------------

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByLabel(hltInputTag, triggerResults);

  if(triggerResults.isValid()) {
    edm::LogInfo("RootTupleMakerV2_PhysicsDSTStream2011_TriggerInfo") << "Successfully obtained " << hltInputTag;
    
    const edm::TriggerNames& names = iEvent.triggerNames(*triggerResults);
    
    for (int i = 0; i < (int) triggerResults->size() ; ++i) { 
      v_hlt_names->push_back ( names.triggerName(i) );
      v_hlt_decisions->push_back ( triggerResults->accept(i) );
      //v_hlt_prescales->push_back ( hltConfig.prescaleValue(iEvent,iSetup,names.triggerName(i)));	
    } 
  } else {
    edm::LogError("RootTupleMakerV2_PhysicsDSTStream2011_TriggerError") << "Error! Can't get the product " << hltInputTag;
  }

  //---------------
  // HLT PFJets
  //---------------

  edm::Handle<reco::PFJetCollection> pfjets;
  iEvent.getByLabel(inputTagHLTPFJets, pfjets);
	
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
      edm::LogInfo("RootTupleMakerV2_PhysicsDSTStream2011Info") << "Can't get the product " << inputTagHLTPFJets;
    }
  
  //---------------
  // HLT CaloJetsRaw
  //---------------

  edm::Handle<reco::CaloJetCollection> calojetsraw;
  iEvent.getByLabel(inputTagHLTCaloJetsRaw, calojetsraw);
	
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
      edm::LogInfo("RootTupleMakerV2_PhysicsDSTStream2011Info") << "Can't get the product " << inputTagHLTCaloJetsRaw;
    }


  //---------------
  // HLT CaloJetsCorr
  //---------------

  edm::Handle<reco::CaloJetCollection> calojetscorr;
  iEvent.getByLabel(inputTagHLTCaloJetsCorr, calojetscorr);
	
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
      edm::LogInfo("RootTupleMakerV2_PhysicsDSTStream2011Info") << "Can't get the product " << inputTagHLTCaloJetsCorr;
    }



  //---------------
  // HLT PixelVertices
  //---------------

  edm::Handle<reco::VertexCollection> pixelvertices;
  iEvent.getByLabel(inputTagHLTPixelVertices, pixelvertices);
	
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
      edm::LogInfo("RootTupleMakerV2_PhysicsDSTStream2011Info") << "Can't get the product " << inputTagHLTPixelVertices;
    }

       
  //-----------------------------------------------------------------
  // put vectors in the event

  // Event Info
  if(storeEventInfo)
    {
      iEvent.put( run,   "run"   );
      iEvent.put( event, "event" );
      iEvent.put( ls   , "ls"    );
      iEvent.put( bunch, "bunch" );
      iEvent.put( orbit, "orbit" );
      iEvent.put( time,  "time"  );
      iEvent.put( isdata,"isData");
    }

  // HLT Info
  iEvent.put ( std::auto_ptr<std::string>(new std::string(hltConfig.tableName())) , "HLTKey");
  iEvent.put ( v_hlt_names      , "HLTTriggerNames"      ) ;
  iEvent.put ( v_hlt_decisions  , "HLTTriggerDecisions"  ) ;
  //   iEvent.put ( v_hlt_prescales  , "HLTTriggerPrescales"  ) ;

  // HLT PFJets
  iEvent.put( etaHLTPFJet, prefixHLTPFJets + "Eta" + suffix );
  iEvent.put( phiHLTPFJet, prefixHLTPFJets + "Phi" + suffix );
  iEvent.put( ptHLTPFJet, prefixHLTPFJets + "Pt" + suffix );
  iEvent.put( energyHLTPFJet, prefixHLTPFJets + "Energy" + suffix );

  iEvent.put( chargedEmEnergyFractionHLTPFJet,  prefixHLTPFJets + "ChargedEmEnergyFraction"  + suffix );
  iEvent.put( chargedHadronEnergyFractionHLTPFJet,  prefixHLTPFJets + "ChargedHadronEnergyFraction"  + suffix );
  iEvent.put( chargedMuEnergyFractionHLTPFJet,  prefixHLTPFJets + "ChargedMuEnergyFraction"  + suffix );
  iEvent.put( electronEnergyFractionHLTPFJet,  prefixHLTPFJets + "ElectronEnergyFraction"  + suffix );
  iEvent.put( muonEnergyFractionHLTPFJet,  prefixHLTPFJets + "MuonEnergyFraction"  + suffix );
  iEvent.put( neutralEmEnergyFractionHLTPFJet,  prefixHLTPFJets + "NeutralEmEnergyFraction"  + suffix );
  iEvent.put( neutralHadronEnergyFractionHLTPFJet,  prefixHLTPFJets + "NeutralHadronEnergyFraction"  + suffix );
  iEvent.put( photonEnergyFractionHLTPFJet,  prefixHLTPFJets + "PhotonEnergyFraction"  + suffix );
  iEvent.put( hfHadronEnergyFractionHLTPFJet,  prefixHLTPFJets + "HFHadronEnergyFraction"  + suffix );
  iEvent.put( hfEMEnergyFractionHLTPFJet,  prefixHLTPFJets + "HFEMEnergyFraction"  + suffix );
  iEvent.put( chargedHadronMultiplicityHLTPFJet,  prefixHLTPFJets + "ChargedHadronMultiplicity"  + suffix );
  iEvent.put( chargedMultiplicityHLTPFJet,  prefixHLTPFJets + "ChargedMultiplicity"  + suffix );
  iEvent.put( electronMultiplicityHLTPFJet,  prefixHLTPFJets + "ElectronMultiplicity"  + suffix );
  iEvent.put( muonMultiplicityHLTPFJet,  prefixHLTPFJets + "MuonMultiplicity"  + suffix );
  iEvent.put( neutralHadronMultiplicityHLTPFJet,  prefixHLTPFJets + "NeutralHadronMultiplicity"  + suffix );
  iEvent.put( neutralMultiplicityHLTPFJet,  prefixHLTPFJets + "NeutralMultiplicity"  + suffix );
  iEvent.put( photonMultiplicityHLTPFJet,  prefixHLTPFJets + "PhotonMultiplicity"  + suffix );
  iEvent.put( hfHadronMultiplicityHLTPFJet,  prefixHLTPFJets + "HFHadronMultiplicity"  + suffix );
  iEvent.put( hfEMMultiplicityHLTPFJet,  prefixHLTPFJets + "HFEMMultiplicity"  + suffix );
  iEvent.put( nConstituents,  prefixHLTPFJets + "NConstituents"  + suffix );
  iEvent.put( passLooseID, prefixHLTPFJets + "PassLooseID" + suffix);
  iEvent.put( passTightID, prefixHLTPFJets + "PassTightID" + suffix);

  // HLT CaloJetsRaw
  iEvent.put( etaHLTCaloJetRaw, prefixHLTCaloJetsRaw + "Eta" + suffix );
  iEvent.put( phiHLTCaloJetRaw, prefixHLTCaloJetsRaw + "Phi" + suffix );
  iEvent.put( ptHLTCaloJetRaw, prefixHLTCaloJetsRaw + "Pt" + suffix );
  iEvent.put( energyHLTCaloJetRaw, prefixHLTCaloJetsRaw + "Energy" + suffix );

  iEvent.put( maxEInEmTowersHLTCaloJetRaw, prefixHLTCaloJetsRaw + "MaxEInEmTowers" + suffix );
  iEvent.put( maxEInHadTowersHLTCaloJetRaw, prefixHLTCaloJetsRaw + "MaxEInHadTowers" + suffix );
  iEvent.put( hadEnergyInHOHLTCaloJetRaw, prefixHLTCaloJetsRaw + "HadEnergyInHO" + suffix );
  iEvent.put( hadEnergyInHBHLTCaloJetRaw, prefixHLTCaloJetsRaw + "HadEnergyInHB" + suffix );
  iEvent.put( hadEnergyInHFHLTCaloJetRaw, prefixHLTCaloJetsRaw + "HadEnergyInHF" + suffix );
  iEvent.put( hadEnergyInHEHLTCaloJetRaw, prefixHLTCaloJetsRaw + "HadEnergyInHE" + suffix );
  iEvent.put( emEnergyInEBHLTCaloJetRaw, prefixHLTCaloJetsRaw + "EmEnergyInEB" + suffix );
  iEvent.put( emEnergyInEEHLTCaloJetRaw, prefixHLTCaloJetsRaw + "EmEnergyInEE" + suffix );
  iEvent.put( emEnergyInHFHLTCaloJetRaw, prefixHLTCaloJetsRaw + "EmEnergyInHF" + suffix );
  iEvent.put( energyFractionHadronicHLTCaloJetRaw, prefixHLTCaloJetsRaw + "EnergyFractionHadronic" + suffix );
  iEvent.put( energyFractionEmHLTCaloJetRaw, prefixHLTCaloJetsRaw + "EnergyFractionEm" + suffix );
  iEvent.put( towersAreaHLTCaloJetRaw, prefixHLTCaloJetsRaw + "TowersArea" + suffix );

  // HLT CaloJetsCorr
  iEvent.put( etaHLTCaloJetCorr, prefixHLTCaloJetsCorr + "Eta" + suffix );
  iEvent.put( phiHLTCaloJetCorr, prefixHLTCaloJetsCorr + "Phi" + suffix );
  iEvent.put( ptHLTCaloJetCorr, prefixHLTCaloJetsCorr + "Pt" + suffix );
  iEvent.put( energyHLTCaloJetCorr, prefixHLTCaloJetsCorr + "Energy" + suffix );

  iEvent.put( maxEInEmTowersHLTCaloJetCorr, prefixHLTCaloJetsCorr + "MaxEInEmTowers" + suffix );
  iEvent.put( maxEInHadTowersHLTCaloJetCorr, prefixHLTCaloJetsCorr + "MaxEInHadTowers" + suffix );
  iEvent.put( hadEnergyInHOHLTCaloJetCorr, prefixHLTCaloJetsCorr + "HadEnergyInHO" + suffix );
  iEvent.put( hadEnergyInHBHLTCaloJetCorr, prefixHLTCaloJetsCorr + "HadEnergyInHB" + suffix );
  iEvent.put( hadEnergyInHFHLTCaloJetCorr, prefixHLTCaloJetsCorr + "HadEnergyInHF" + suffix );
  iEvent.put( hadEnergyInHEHLTCaloJetCorr, prefixHLTCaloJetsCorr + "HadEnergyInHE" + suffix );
  iEvent.put( emEnergyInEBHLTCaloJetCorr, prefixHLTCaloJetsCorr + "EmEnergyInEB" + suffix );
  iEvent.put( emEnergyInEEHLTCaloJetCorr, prefixHLTCaloJetsCorr + "EmEnergyInEE" + suffix );
  iEvent.put( emEnergyInHFHLTCaloJetCorr, prefixHLTCaloJetsCorr + "EmEnergyInHF" + suffix );
  iEvent.put( energyFractionHadronicHLTCaloJetCorr, prefixHLTCaloJetsCorr + "EnergyFractionHadronic" + suffix );
  iEvent.put( energyFractionEmHLTCaloJetCorr, prefixHLTCaloJetsCorr + "EnergyFractionEm" + suffix );
  iEvent.put( towersAreaHLTCaloJetCorr, prefixHLTCaloJetsCorr + "TowersArea" + suffix );

  // HLT PixelVertices
  iEvent.put( chi2HLTPixelVertices , prefixHLTPixelVertices + "Chi2" + suffix );
  iEvent.put( normChi2HLTPixelVertices , prefixHLTPixelVertices + "NormChi2" + suffix );
  iEvent.put( ndofHLTPixelVertices , prefixHLTPixelVertices + "Ndof" + suffix );
  iEvent.put( XCoordHLTPixelVertices , prefixHLTPixelVertices + "XCoord" + suffix );
  iEvent.put( XErrorHLTPixelVertices , prefixHLTPixelVertices + "XError" + suffix );
  iEvent.put( YCoordHLTPixelVertices , prefixHLTPixelVertices + "YCoord" + suffix );
  iEvent.put( YErrorHLTPixelVertices , prefixHLTPixelVertices + "YError" + suffix );
  iEvent.put( ZCoordHLTPixelVertices , prefixHLTPixelVertices + "ZCoord" + suffix );
  iEvent.put( ZErrorHLTPixelVertices , prefixHLTPixelVertices + "ZError" + suffix );
  iEvent.put( isValidHLTPixelVertices , prefixHLTPixelVertices + "IsValid" + suffix );

}
