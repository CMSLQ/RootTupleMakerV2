#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_GenJets.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


RootTupleMakerV2_GenJets::RootTupleMakerV2_GenJets(const edm::ParameterSet& iConfig) :
  genJetsInputToken_ (consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("InputTag"))),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize"))
{
  produces <std::vector<float> > ( prefix + "Eta" + suffix );
  produces <std::vector<float> > ( prefix + "Phi" + suffix );
  produces <std::vector<float> > ( prefix + "P" + suffix );
  produces <std::vector<float> > ( prefix + "Pt" + suffix );
  produces <std::vector<float> > ( prefix + "Energy" + suffix );
  produces <std::vector<float> > ( prefix + "EMF" + suffix );
  produces <std::vector<float> > ( prefix + "HADF" + suffix );
}

void RootTupleMakerV2_GenJets::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::unique_ptr<std::vector<float> >  eta  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >  phi  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >  p  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >  pt  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >  energy  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >  emf  ( new std::vector<float>()  );
  std::unique_ptr<std::vector<float> >  hadf  ( new std::vector<float>()  );

  //-----------------------------------------------------------------
  if( !iEvent.isRealData() ) {
    edm::Handle<reco::GenJetCollection> genJets;
    iEvent.getByToken(genJetsInputToken_, genJets);

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
      edm::LogError("RootTupleMakerV2_GenJetsError") << "Error! Can't get the genJetsInputToken_";
    }
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put(std::move(eta), prefix + "Eta" + suffix );
  iEvent.put(std::move(phi), prefix + "Phi" + suffix );
  iEvent.put(std::move(p), prefix + "P" + suffix );
  iEvent.put(std::move(pt), prefix + "Pt" + suffix );
  iEvent.put(std::move(energy), prefix + "Energy" + suffix );
  iEvent.put(std::move(emf), prefix + "EMF" + suffix );
  iEvent.put(std::move(hadf), prefix + "HADF" + suffix );
}
