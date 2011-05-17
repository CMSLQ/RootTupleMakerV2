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
  produces <std::vector<double> > ( prefix + "EcalIso" + suffix );
  produces <std::vector<double> > ( prefix + "HcalIso" + suffix );
  produces <std::vector<double> > ( prefix + "HoE" + suffix );
  produces <std::vector<double> > ( prefix + "TrkIso" + suffix );
  produces <std::vector<double> > ( prefix + "SigmaIEtaIEta" + suffix );
  produces <std::vector<bool> >   ( prefix + "TrkVeto" + suffix );
  produces <std::vector<double> > ( prefix + "SCseedEnergy" + suffix );
  produces <std::vector<double> > ( prefix + "SCenergy" + suffix );
  produces <std::vector<double> > ( prefix + "SCeta" + suffix );
  produces <std::vector<double> > ( prefix + "SCphi" + suffix );
  produces <std::vector<double> > ( prefix + "E3x3" + suffix );
  produces <std::vector<double> > ( prefix + "E5x5" + suffix );
}

void RootTupleMakerV2_Photons::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  eta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  phi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  pt  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  energy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  ecalIso  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hcalIso  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  hoe  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  trkIso  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  sigmaIetaIeta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<bool> >    trkVeto  ( new std::vector<bool>()  );
  std::auto_ptr<std::vector<double> >  SCseedEnergy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  SCenergy  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  SCeta  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  SCphi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  E3x3  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  E5x5  ( new std::vector<double>()  );

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
      ecalIso->push_back( it->ecalIso() );
      hcalIso->push_back( it->hcalIso() );
      hoe->push_back( it->hadronicOverEm() );
      trkIso->push_back( it->trkSumPtHollowConeDR04() );
      sigmaIetaIeta->push_back( it->sigmaIetaIeta() );
      trkVeto->push_back( it->hasPixelSeed() );
      SCseedEnergy->push_back( it->superCluster()->seed()->energy() );
      SCenergy->push_back( it->superCluster()->energy() );
      SCeta->push_back( it->superCluster()->eta() );
      SCphi->push_back( it->superCluster()->phi() );
      E3x3->push_back( it->e3x3() );
      E5x5->push_back( it->e5x5() );

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
  iEvent.put( ecalIso, prefix + "EcalIso" + suffix );
  iEvent.put( hcalIso, prefix + "HcalIso" + suffix );
  iEvent.put( hoe, prefix + "HoE" + suffix );
  iEvent.put( trkIso, prefix + "TrkIso" + suffix );
  iEvent.put( sigmaIetaIeta, prefix + "SigmaIEtaIEta" + suffix );
  iEvent.put( trkVeto, prefix + "TrkVeto" + suffix );
  iEvent.put( SCseedEnergy, prefix + "SCseedEnergy" + suffix );
  iEvent.put( SCenergy, prefix + "SCenergy" + suffix );
  iEvent.put( SCeta, prefix + "SCeta" + suffix );
  iEvent.put( SCphi, prefix + "SCphi" + suffix );
  iEvent.put( E3x3, prefix + "E3x3" + suffix );
  iEvent.put( E5x5, prefix + "E5x5" + suffix );
}
