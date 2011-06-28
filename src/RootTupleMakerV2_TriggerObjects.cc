#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_TriggerObjects.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"

RootTupleMakerV2_TriggerObjects::RootTupleMakerV2_TriggerObjects(const edm::ParameterSet& iConfig) :
    inputTag   (iConfig.getParameter<edm::InputTag>("InputTag")),
    filterID   (iConfig.getParameter<edm::InputTag>("FilterID")),
    prefix     (iConfig.getParameter<std::string>  ("Prefix")),
    suffix     (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize    (iConfig.getParameter<unsigned int> ("MaxSize"))
{
  produces<std::vector<double> > ( prefix + "Pt"  + suffix );
  produces<std::vector<double> > ( prefix + "Eta" + suffix );
  produces<std::vector<double> > ( prefix + "Phi" + suffix );
}

void RootTupleMakerV2_TriggerObjects::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> > pt  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> > eta ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> > phi ( new std::vector<double>()  );
  
  edm::Handle<trigger::TriggerEvent> triggerEvent;
  iEvent.getByLabel( inputTag ,triggerEvent);

  if ( triggerEvent.isValid() ){

    const trigger::TriggerObjectCollection & triggerObjects = triggerEvent -> getObjects();
    trigger::size_type filter_index = triggerEvent -> filterIndex ( filterID ) ;
    trigger::size_type n_filters    = triggerEvent -> sizeFilters();
    trigger::size_type n_colls      = triggerEvent -> sizeCollections();

    for (trigger::size_type i = 0; i < n_colls ; ++i) {
      std::cout << triggerEvent -> collectionTag(i) << std::endl;
    }

    if ( filter_index < n_filters ) {
      const trigger::Keys & triggerKeys ( triggerEvent -> filterKeys ( filter_index ) );
      const int nkeys = triggerKeys.size();

      for (int ikey = 0; ikey < nkeys; ++ikey ) {
	const trigger::TriggerObject& triggerObject = triggerObjects[ triggerKeys [ ikey ] ];
	pt  -> push_back ( triggerObject.pt () );
	eta -> push_back ( triggerObject.eta() );
	phi -> push_back ( triggerObject.phi() );
      }
    }
  } else { 
    edm::LogError("RootTupleMakerV2_TriggerObjectsError") << "Error! Can't get the product " << inputTag;
  }

  
  
  iEvent.put( pt , prefix + "Pt"  + suffix );
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
}
