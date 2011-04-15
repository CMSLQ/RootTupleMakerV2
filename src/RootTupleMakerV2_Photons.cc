#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Photons.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
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


RootTupleMakerV2_Photons::RootTupleMakerV2_Photons(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize"))
{
  produces <std::vector<double> > ( prefix + "Eta" + suffix );
  produces <std::vector<double> > ( prefix + "Phi" + suffix );
  produces <std::vector<double> > ( prefix + "Pt" + suffix );
  produces <std::vector<double> > ( prefix + "Energy" + suffix );
}

void RootTupleMakerV2_Photons::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );

  //-----------------------------------------------------------------

  edm::Handle<std::vector<pat::Photon> > photons;
  iEvent.getByLabel(inputTag, photons);

  if(photons.isValid()) {
    edm::LogInfo("RootTupleMakerV2_PhotonsInfo") << "Total # Photons: " << photons->size();

    for( std::vector<pat::Photon>::const_iterator it = photons->begin(); it != photons->end(); ++it ) {
      // exit from loop when you reach the required number of photons
      if(eta->size() >= maxSize)
        break;

      eta->push_back( it->eta() );
      phi->push_back( it->phi() );
      pt->push_back( it->pt() );
      energy->push_back( it->energy() );

    }
  } else {
    edm::LogError("RootTupleMakerV2_PhotonsError") << "Error! Can't get the product " << inputTag;
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( pt, prefix + "Pt" + suffix );
  iEvent.put( energy, prefix + "Energy" + suffix );
}
