#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_PFJets.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"


RootTupleMakerV2_PFJets::RootTupleMakerV2_PFJets(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize")),
    applyResJEC (iConfig.getParameter<bool>     ("ApplyResidualJEC")),
    resJEC (iConfig.getParameter<std::string>   ("ResidualJEC"))
{
  produces <std::vector<double> > ( prefix + "Eta" + suffix );
  produces <std::vector<double> > ( prefix + "Phi" + suffix );
  produces <std::vector<double> > ( prefix + "Pt" + suffix );
  produces <std::vector<double> > ( prefix + "PtRaw" + suffix );
  produces <std::vector<double> > ( prefix + "Energy" + suffix );
  produces <std::vector<double> > ( prefix + "EnergyRaw" + suffix );
  produces <std::vector<int> >    ( prefix + "PartonFlavour" + suffix );
  produces <std::vector<double> > ( prefix + "ChargedHadronEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefix + "NeutralHadronEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefix + "ChargedEmEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefix + "NeutralEmEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefix + "ChargedMuEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefix + "ChargedMultiplicity"  + suffix );
  produces <std::vector<double> > ( prefix + "NeutralMultiplicity"  + suffix );
  produces <std::vector<double> > ( prefix + "MuonMultiplicity"  + suffix );
}

void RootTupleMakerV2_PFJets::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt_raw  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy_raw ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     partonFlavour  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  chargedHadronEnergyFraction  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  neutralHadronEnergyFraction  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  chargedEmEnergyFraction  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  neutralEmEnergyFraction  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  chargedMuEnergyFraction  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  chargedMultiplicity  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  neutralMultiplicity  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  muonMultiplicity  ( new std::vector<double>()  ) ;

  //-----------------------------------------------------------------
  edm::Handle<std::vector<pat::Jet> > jets;
  iEvent.getByLabel(inputTag, jets);

  if(jets.isValid()) {
    edm::LogInfo("RootTupleMakerV2_PFJetsInfo") << "Total # PFJets: " << jets->size();

    JetCorrectorParameters *ResJetCorPar = 0;
    FactorizedJetCorrector *JEC = 0;
    if(applyResJEC) {
      edm::FileInPath fipRes(resJEC);
      ResJetCorPar = new JetCorrectorParameters(fipRes.fullPath());
      std::vector<JetCorrectorParameters> vParam;
      vParam.push_back(*ResJetCorPar);
      JEC = new FactorizedJetCorrector(vParam);
    }

    for( std::vector<pat::Jet>::const_iterator it = jets->begin(); it != jets->end(); ++it ) {
      // exit from loop when you reach the required number of jets
      if(eta->size() >= maxSize)
        break;

      double corr = 1.;
      if( applyResJEC && iEvent.isRealData() ) {
        JEC->setJetEta( it->eta() );
        JEC->setJetPt( it->pt() ); // here you put the L2L3 Corrected jet pt
        corr = JEC->getCorrection();
      }

      // fill in all the vectors
      eta->push_back( it->eta() );
      phi->push_back( it->phi() );
      pt->push_back( it->pt()*corr );
      pt_raw->push_back( it->correctedJet("raw").pt() );
      energy->push_back( it->energy()*corr );
      energy_raw->push_back( it->correctedJet("raw").energy() );
      partonFlavour->push_back( it->partonFlavour() );
      chargedHadronEnergyFraction->push_back( it->chargedHadronEnergyFraction() );
      neutralHadronEnergyFraction->push_back( it->neutralHadronEnergyFraction() );
      chargedEmEnergyFraction->push_back( it->chargedEmEnergyFraction() );
      neutralEmEnergyFraction->push_back( it->neutralEmEnergyFraction() );
      chargedMuEnergyFraction->push_back( it->chargedMuEnergyFraction() );
      chargedMultiplicity->push_back( it->chargedMultiplicity() );
      neutralMultiplicity->push_back( it->neutralMultiplicity() );
      muonMultiplicity->push_back( it->muonMultiplicity() );
    }
  } else {
    edm::LogError("RootTupleMakerV2_PFJetsError") << "Error! Can't get the product " << inputTag;
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( pt, prefix + "Pt" + suffix );
  iEvent.put( pt_raw, prefix + "PtRaw" + suffix );
  iEvent.put( energy, prefix + "Energy" + suffix );
  iEvent.put( energy_raw, prefix + "EnergyRaw" + suffix );
  iEvent.put( partonFlavour, prefix + "PartonFlavour" + suffix );
  iEvent.put( chargedHadronEnergyFraction,  prefix + "ChargedHadronEnergyFraction"  + suffix );
  iEvent.put( neutralHadronEnergyFraction,  prefix + "NeutralHadronEnergyFraction"  + suffix );
  iEvent.put( chargedEmEnergyFraction,  prefix + "ChargedEmEnergyFraction"  + suffix );
  iEvent.put( neutralEmEnergyFraction,  prefix + "NeutralEmEnergyFraction"  + suffix );
  iEvent.put( chargedMuEnergyFraction,  prefix + "ChargedMuEnergyFraction"  + suffix );
  iEvent.put( chargedMultiplicity,  prefix + "ChargedMultiplicity"  + suffix );
  iEvent.put( neutralMultiplicity,  prefix + "NeutralMultiplicity"  + suffix );
  iEvent.put( muonMultiplicity,  prefix + "MuonMultiplicity"  + suffix );
}
