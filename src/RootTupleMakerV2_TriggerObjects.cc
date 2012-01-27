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
  produces<std::vector<int   > > ( prefix + "ID"     + suffix );
  produces<std::vector<double> > ( prefix + "Pt"     + suffix );
  produces<std::vector<double> > ( prefix + "Et"     + suffix );
  produces<std::vector<double> > ( prefix + "Eta"    + suffix );
  produces<std::vector<double> > ( prefix + "Phi"    + suffix );
  produces<std::vector<double> > ( prefix + "Energy" + suffix );
  produces<std::vector<double> > ( prefix + "Mass"   + suffix );

}

void RootTupleMakerV2_TriggerObjects::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  std::auto_ptr<std::vector<int>    > id     ( new std::vector<int   >()  );
  std::auto_ptr<std::vector<double> > pt     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> > et     ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> > eta    ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> > phi    ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> > energy ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> > mass   ( new std::vector<double>()  );
  
  edm::Handle<trigger::TriggerEvent> triggerEvent;
  iEvent.getByLabel( inputTag ,triggerEvent);

  if ( triggerEvent.isValid() ){

    const trigger::TriggerObjectCollection & triggerObjects = triggerEvent -> getObjects();
    trigger::size_type filter_index = triggerEvent -> filterIndex ( filterID ) ;
    trigger::size_type n_filters    = triggerEvent -> sizeFilters();
    // trigger::size_type n_colls      = triggerEvent -> sizeCollections();
    
    if ( filter_index < n_filters ) {
      const trigger::Keys & triggerKeys ( triggerEvent -> filterKeys ( filter_index ) );
      const int nkeys = triggerKeys.size();

      for (int ikey = 0; ikey < nkeys; ++ikey ) {
	const trigger::TriggerObject& triggerObject = triggerObjects[ triggerKeys [ ikey ] ];
	id     -> push_back ( triggerObject.id    ());
	pt     -> push_back ( triggerObject.pt    ());
	et     -> push_back ( triggerObject.et    ());
	eta    -> push_back ( triggerObject.eta   ());
	phi    -> push_back ( triggerObject.phi   ());
	energy -> push_back ( triggerObject.energy());
	mass   -> push_back ( triggerObject.mass  ());
      }
    }
  } else { 
    edm::LogError("RootTupleMakerV2_TriggerObjectsError") << "Error! Can't get the product " << inputTag;
  }

  iEvent.put( id    , prefix + "ID"     + suffix );
  iEvent.put( pt    , prefix + "Pt"     + suffix );
  iEvent.put( et    , prefix + "Et"     + suffix );
  iEvent.put( eta   , prefix + "Eta"    + suffix );
  iEvent.put( phi   , prefix + "Phi"    + suffix );
  iEvent.put( energy, prefix + "Energy" + suffix );
  iEvent.put( mass  , prefix + "Mass"   + suffix );

}
