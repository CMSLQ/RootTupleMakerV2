#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_GenJets.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"


RootTupleMakerV2_GenJets::RootTupleMakerV2_GenJets(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize"))
{
  produces <std::vector<double> > ( prefix + "Eta" + suffix );
  produces <std::vector<double> > ( prefix + "Phi" + suffix );
  produces <std::vector<double> > ( prefix + "P" + suffix );
  produces <std::vector<double> > ( prefix + "Pt" + suffix );
  produces <std::vector<double> > ( prefix + "Energy" + suffix );
  produces <std::vector<double> > ( prefix + "EMF" + suffix );
  produces <std::vector<double> > ( prefix + "HADF" + suffix );
}

void RootTupleMakerV2_GenJets::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  p  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  emf  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hadf  ( new std::vector<double>()  );

  //-----------------------------------------------------------------
  if( !iEvent.isRealData() ) {
    edm::Handle<reco::GenJetCollection> genJets;
    iEvent.getByLabel(inputTag, genJets);

    if( genJets.isValid() ) {
      edm::LogInfo("RootTupleMakerV2_GenJetsInfo") << "Total # GenJets: " << genJets->size();

      for( reco::GenJetCollection::const_iterator it = genJets->begin(); it != genJets->end(); ++it ) {
        // exit from loop when you reach the required number of GenJets
        if(eta->size() >= maxSize)
          break;

        // fill in all the vectors
        eta->push_back( it->eta() );
        phi->push_back( it->phi() );
        p->push_back( it->p() );
        pt->push_back( it->pt() );
        energy->push_back( it->energy() );
        emf->push_back( it->emEnergy()/it->energy() );
        hadf->push_back( it->hadEnergy()/it->energy() );
      }
    } else {
      edm::LogError("RootTupleMakerV2_GenJetsError") << "Error! Can't get the product " << inputTag;
    }
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( p, prefix + "P" + suffix );
  iEvent.put( pt, prefix + "Pt" + suffix );
  iEvent.put( energy, prefix + "Energy" + suffix );
  iEvent.put( emf, prefix + "EMF" + suffix );
  iEvent.put( hadf, prefix + "HADF" + suffix );
}
