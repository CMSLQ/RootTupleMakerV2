#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_GenParticles.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


RootTupleMakerV2_GenParticles::RootTupleMakerV2_GenParticles(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize"))
{
  produces <std::vector<double> > ( prefix + "Eta" + suffix );
  produces <std::vector<double> > ( prefix + "Phi" + suffix );
  produces <std::vector<double> > ( prefix + "P" + suffix );
  produces <std::vector<double> > ( prefix + "Px" + suffix );
  produces <std::vector<double> > ( prefix + "Py" + suffix );
  produces <std::vector<double> > ( prefix + "Pz" + suffix );
  produces <std::vector<double> > ( prefix + "Pt" + suffix );
  produces <std::vector<double> > ( prefix + "Energy" + suffix );
  produces <std::vector<int> >    ( prefix + "PdgId" + suffix );
  produces <std::vector<double> > ( prefix + "VX" + suffix );
  produces <std::vector<double> > ( prefix + "VY" + suffix );
  produces <std::vector<double> > ( prefix + "VZ" + suffix );
  produces <std::vector<int> >    ( prefix + "NumDaught" + suffix );
  produces <std::vector<int> >    ( prefix + "Status" + suffix );
  produces <std::vector<int> >    ( prefix + "MotherIndex" + suffix );
}

void RootTupleMakerV2_GenParticles::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  p  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  px  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  py  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pz  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     pdgId ( new std::vector<int>()  );
  std::auto_ptr<std::vector<double> >  vx  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  vy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  vz  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<int> >     numDaught  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >     status  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >     motherIndex  ( new std::vector<int>()  );
  
  //-----------------------------------------------------------------
  if( !iEvent.isRealData() ) {
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel(inputTag, genParticles);

    if( genParticles.isValid() ) {
      edm::LogInfo("RootTupleMakerV2_GenParticlesInfo") << "Total # GenParticles: " << genParticles->size();

      for( reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it ) {
        // exit from loop when you reach the required number of GenParticles
        if(eta->size() > maxSize)
          break;

        // fill in all the vectors
        eta->push_back( it->eta() );
        phi->push_back( it->phi() );
        p->push_back( it->p() );
        px->push_back( it->px() );
        py->push_back( it->py() );
        pz->push_back( it->pz() );
        pt->push_back( it->pt() );
        energy->push_back( it->energy() );
        pdgId->push_back( it->pdgId() );
        vx->push_back( it->vx() );
        vy->push_back( it->vy() );
        vz->push_back( it->vz() );
        numDaught->push_back( it->numberOfDaughters() );
        status->push_back( it->status() );

        int idx = -1;
        for( reco::GenParticleCollection::const_iterator mit = genParticles->begin(); mit != genParticles->end(); ++mit ) {
          if( it->mother()==&(*mit) ) {
             idx = std::distance(genParticles->begin(),mit);
             break;
          }
        }
        motherIndex->push_back( idx );
      }
    } else {
      edm::LogError("RootTupleMakerV2_GenParticlesError") << "Error! Can't get the product " << inputTag;
    }
  }
  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( p, prefix + "P" + suffix );
  iEvent.put( px, prefix + "Px" + suffix );
  iEvent.put( py, prefix + "Py" + suffix );
  iEvent.put( pz, prefix + "Pz" + suffix );
  iEvent.put( pt, prefix + "Pt" + suffix );
  iEvent.put( energy, prefix + "Energy" + suffix );
  iEvent.put( pdgId, prefix + "PdgId" + suffix );
  iEvent.put( vx, prefix + "VX" + suffix );
  iEvent.put( vy, prefix + "VY" + suffix );
  iEvent.put( vz, prefix + "VZ" + suffix );
  iEvent.put( numDaught, prefix + "NumDaught" + suffix );
  iEvent.put( status, prefix + "Status" + suffix );
  iEvent.put( motherIndex, prefix + "MotherIndex" + suffix );
}
