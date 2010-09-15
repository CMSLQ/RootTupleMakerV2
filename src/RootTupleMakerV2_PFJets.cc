#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_PFJets.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"


RootTupleMakerV2_PFJets::RootTupleMakerV2_PFJets(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize")),
    jecUncPath(iConfig.getParameter<std::string>("JECUncertainty")),
    applyResJEC (iConfig.getParameter<bool>     ("ApplyResidualJEC")),
    resJEC (iConfig.getParameter<std::string>   ("ResidualJEC"))
{
  produces <std::vector<double> > ( prefix + "Eta" + suffix );
  produces <std::vector<double> > ( prefix + "Phi" + suffix );
  produces <std::vector<double> > ( prefix + "Pt" + suffix );
  produces <std::vector<double> > ( prefix + "PtRaw" + suffix );
  produces <std::vector<double> > ( prefix + "Energy" + suffix );
  produces <std::vector<double> > ( prefix + "EnergyRaw" + suffix );
  produces <std::vector<double> > ( prefix + "JECUnc" + suffix );
  produces <std::vector<double> > ( prefix + "ResJEC" + suffix );
  produces <std::vector<int> >    ( prefix + "PartonFlavour" + suffix );
  produces <std::vector<double> > ( prefix + "ChargedEmEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefix + "ChargedHadronEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefix + "ChargedMuEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefix + "ElectronEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefix + "MuonEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefix + "NeutralEmEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefix + "NeutralHadronEnergyFraction"  + suffix );
  produces <std::vector<double> > ( prefix + "PhotonEnergyFraction"  + suffix );
  produces <std::vector<int> >    ( prefix + "ChargedHadronMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefix + "ChargedMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefix + "ElectronMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefix + "MuonMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefix + "NeutralHadronMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefix + "NeutralMultiplicity"  + suffix );
  produces <std::vector<int> >    ( prefix + "PhotonMultiplicity"  + suffix );
}

void RootTupleMakerV2_PFJets::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt_raw  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy_raw ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  jecUnc_vec ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  resJEC_vec ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     partonFlavour  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  chargedEmEnergyFraction  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  chargedHadronEnergyFraction  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  chargedMuEnergyFraction  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  electronEnergyFraction  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  muonEnergyFraction  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  neutralEmEnergyFraction  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  neutralHadronEnergyFraction  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<double> >  photonEnergyFraction  ( new std::vector<double>()  ) ;
  std::auto_ptr<std::vector<int> >     chargedHadronMultiplicity  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     chargedMultiplicity  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     electronMultiplicity  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     muonMultiplicity  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     neutralHadronMultiplicity  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     neutralMultiplicity  ( new std::vector<int>()  ) ;
  std::auto_ptr<std::vector<int> >     photonMultiplicity  ( new std::vector<int>()  ) ;

  //-----------------------------------------------------------------
  edm::FileInPath fipUnc(jecUncPath);;
  JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(fipUnc.fullPath());

  JetCorrectorParameters *ResJetCorPar = 0;
  FactorizedJetCorrector *JEC = 0;
  if(applyResJEC) {
    edm::FileInPath fipRes(resJEC);
    ResJetCorPar = new JetCorrectorParameters(fipRes.fullPath());
    std::vector<JetCorrectorParameters> vParam;
    vParam.push_back(*ResJetCorPar);
    JEC = new FactorizedJetCorrector(vParam);
  }

  edm::Handle<std::vector<pat::Jet> > jets;
  iEvent.getByLabel(inputTag, jets);

  if(jets.isValid()) {
    edm::LogInfo("RootTupleMakerV2_PFJetsInfo") << "Total # PFJets: " << jets->size();

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

      jecUnc->setJetEta( it->eta() );
      jecUnc->setJetPt( it->pt()*corr ); // the uncertainty is a function of the corrected pt

      // fill in all the vectors
      eta->push_back( it->eta() );
      phi->push_back( it->phi() );
      pt->push_back( it->pt()*corr );
      pt_raw->push_back( it->correctedJet("raw").pt() );
      energy->push_back( it->energy()*corr );
      energy_raw->push_back( it->correctedJet("raw").energy() );
      jecUnc_vec->push_back( jecUnc->getUncertainty(true) );
      resJEC_vec->push_back( corr );
      partonFlavour->push_back( it->partonFlavour() );
      chargedEmEnergyFraction->push_back( it->chargedEmEnergyFraction() );
      chargedHadronEnergyFraction->push_back( it->chargedHadronEnergyFraction() );
      chargedMuEnergyFraction->push_back( it->chargedMuEnergyFraction() );
      electronEnergyFraction->push_back( it->electronEnergy()/it->energy() );
      muonEnergyFraction->push_back( it->muonEnergyFraction() );
      neutralEmEnergyFraction->push_back( it->neutralEmEnergyFraction() );
      neutralHadronEnergyFraction->push_back( it->neutralHadronEnergyFraction() );
      photonEnergyFraction->push_back( it->photonEnergyFraction() );
      chargedHadronMultiplicity->push_back( it->chargedHadronMultiplicity() );
      chargedMultiplicity->push_back( it->chargedMultiplicity() );
      electronMultiplicity->push_back( it->electronMultiplicity() );
      muonMultiplicity->push_back( it->muonMultiplicity() );
      neutralHadronMultiplicity->push_back( it->neutralHadronMultiplicity() );
      neutralMultiplicity->push_back( it->neutralMultiplicity() );
      photonMultiplicity->push_back( it->photonMultiplicity() );
    }
  } else {
    edm::LogError("RootTupleMakerV2_PFJetsError") << "Error! Can't get the product " << inputTag;
  }

  delete jecUnc;
  delete ResJetCorPar;
  delete JEC;

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( pt, prefix + "Pt" + suffix );
  iEvent.put( pt_raw, prefix + "PtRaw" + suffix );
  iEvent.put( energy, prefix + "Energy" + suffix );
  iEvent.put( energy_raw, prefix + "EnergyRaw" + suffix );
  iEvent.put( jecUnc_vec, prefix + "JECUnc" + suffix );
  iEvent.put( resJEC_vec, prefix + "ResJEC" + suffix );
  iEvent.put( partonFlavour, prefix + "PartonFlavour" + suffix );
  iEvent.put( chargedEmEnergyFraction,  prefix + "ChargedEmEnergyFraction"  + suffix );
  iEvent.put( chargedHadronEnergyFraction,  prefix + "ChargedHadronEnergyFraction"  + suffix );
  iEvent.put( chargedMuEnergyFraction,  prefix + "ChargedMuEnergyFraction"  + suffix );
  iEvent.put( electronEnergyFraction,  prefix + "ElectronEnergyFraction"  + suffix );
  iEvent.put( muonEnergyFraction,  prefix + "MuonEnergyFraction"  + suffix );
  iEvent.put( neutralEmEnergyFraction,  prefix + "NeutralEmEnergyFraction"  + suffix );
  iEvent.put( neutralHadronEnergyFraction,  prefix + "NeutralHadronEnergyFraction"  + suffix );
  iEvent.put( photonEnergyFraction,  prefix + "PhotonEnergyFraction"  + suffix );
  iEvent.put( chargedHadronMultiplicity,  prefix + "ChargedHadronMultiplicity"  + suffix );
  iEvent.put( chargedMultiplicity,  prefix + "ChargedMultiplicity"  + suffix );
  iEvent.put( electronMultiplicity,  prefix + "ElectronMultiplicity"  + suffix );
  iEvent.put( muonMultiplicity,  prefix + "MuonMultiplicity"  + suffix );
  iEvent.put( neutralHadronMultiplicity,  prefix + "NeutralHadronMultiplicity"  + suffix );
  iEvent.put( neutralMultiplicity,  prefix + "NeutralMultiplicity"  + suffix );
  iEvent.put( photonMultiplicity,  prefix + "PhotonMultiplicity"  + suffix );
}
