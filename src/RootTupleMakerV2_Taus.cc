#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Taus.h"
#include "Leptoquarks/RootTupleMakerV2/interface/PatUtilities.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

RootTupleMakerV2_Taus::RootTupleMakerV2_Taus(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize"))
{
  produces <std::vector<double> > ( prefix + "Eta" + suffix );
  produces <std::vector<double> > ( prefix + "Phi" + suffix );
  produces <std::vector<double> > ( prefix + "Pt"  + suffix );
}

void RootTupleMakerV2_Taus::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()   );

  edm::Handle<std::vector<pat::Tau> > taus;
  iEvent.getByLabel(inputTag, taus);
  
  if(taus.isValid()) {
    edm::LogInfo("RootTupleMakerV2_TausInfo") << "Total # Taus: " << taus->size();

    std::vector<pat::Tau>::const_iterator it     = taus -> begin();
    std::vector<pat::Tau>::const_iterator it_end = taus -> end();

    for (; it != it_end; ++it ) { 
      if ( eta->size() > maxSize ) break;

      eta -> push_back ( it -> eta() ) ;
      phi -> push_back ( it -> phi() ) ;
      pt  -> push_back ( it -> pt () ) ;
    }
  } else {
    edm::LogError("RootTupleMakerV2_TausError") << "Error! Can't get the product " << inputTag;
  }
  
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( pt , prefix + "Pt"  + suffix );
}
