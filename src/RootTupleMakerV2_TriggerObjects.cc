#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_TriggerObjects.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "TLorentzVector.h"

RootTupleMakerV2_TriggerObjects::RootTupleMakerV2_TriggerObjects(const edm::ParameterSet& iConfig) :
    inputTag   (iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix     (iConfig.getParameter<std::string>  ("Prefix")),
    suffix     (iConfig.getParameter<std::string>  ("Suffix"))
{
  //------------------------------------------------------------------------
  // What variables will this producer push into the event?
  //------------------------------------------------------------------------
  
  produces <std::vector<std::string> >           ( prefix + "FilterName"   + suffix );
  produces <std::vector<std::vector< int   > > > ( prefix + "FilterObjId"  + suffix );
  produces <std::vector<std::vector< float > > > ( prefix + "FilterObjPt"  + suffix );
  produces <std::vector<std::vector< float > > > ( prefix + "FilterObjEta" + suffix );
  produces <std::vector<std::vector< float > > > ( prefix + "FilterObjPhi" + suffix );

}

void RootTupleMakerV2_TriggerObjects::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  //------------------------------------------------------------------------
  // Declare items to push into the event
  //------------------------------------------------------------------------
  
  std::auto_ptr<std::vector < std::string          > > v_filter_names( new std::vector< std::string          > ());
  std::auto_ptr<std::vector < std::vector< int   > > > v_filter_ids  ( new std::vector< std::vector< int   > > ());
  std::auto_ptr<std::vector < std::vector< float > > > v_filter_pts  ( new std::vector< std::vector< float > > ());
  std::auto_ptr<std::vector < std::vector< float > > > v_filter_etas ( new std::vector< std::vector< float > > ());
  std::auto_ptr<std::vector < std::vector< float > > > v_filter_phis ( new std::vector< std::vector< float > > ());
  
  //------------------------------------------------------------------------
  // Get the trigger event and make sure it's valid
  //------------------------------------------------------------------------

  edm::Handle<trigger::TriggerEvent> triggerEvent;
  iEvent.getByLabel( inputTag , triggerEvent);
  
  if ( triggerEvent.isValid() ){

    //------------------------------------------------------------------------
    // Get the trigger objects in the event
    //------------------------------------------------------------------------

    const trigger::TriggerObjectCollection & triggerObjects = triggerEvent -> getObjects();
    
    //------------------------------------------------------------------------
    // Loop over the filters in the trigger event
    // e.g.: hltEle27WP80TrackIsoFilter
    //------------------------------------------------------------------------
    
    size_t nFilters       = triggerEvent -> sizeFilters();
    size_t nFiltersPassed = 0;
    size_t iFilter        = 0;

    for (; iFilter < nFilters; ++iFilter) {

      //------------------------------------------------------------------------
      // Find information for each filter:
      // -  Name  : name of the filter, e.g. hltEle27WP80TrackIsoFilter
      // - "Keys" : std::vector<uint16_t> storing indices of trigger objects that pass the filter
      //------------------------------------------------------------------------
      
      std::string          name = triggerEvent -> filterTag ( iFilter ).label();
      const trigger::Keys& keys = triggerEvent -> filterKeys( iFilter );
      
      //------------------------------------------------------------------------
      // Loop over the keys to get to the trigger objects that pass the filter
      //------------------------------------------------------------------------
      
      trigger::Keys::const_iterator iKey   = keys.begin();
      trigger::Keys::const_iterator endKey = keys.end();

      // useful variables
      int nFilterObjects = 0;
      std::vector<TLorentzVector> triggerObjectP4s;
      std::vector<int>            triggerObjectIds;
      
      for (; iKey != endKey ; ++iKey ) {

	// Keep track of how many filter objects there are
	
	++nFilterObjects;

	// Get the trigger object associated with this key

	const trigger::TriggerObject & triggerObject = triggerObjects [*iKey];

	// Store the trigger object as a TLorentzVector

	TLorentzVector p4;
	p4.SetPtEtaPhiM ( triggerObject.pt  (),
			  triggerObject.eta (),
			  triggerObject.phi (),
			  triggerObject.mass() );

	triggerObjectP4s.push_back ( p4 ) ;
	triggerObjectIds.push_back ( int(triggerObject.id()));
	
      } // end loop over keys/trigger objects passing filters

      //------------------------------------------------------------------------
      // If the filter passed, store its information
      //------------------------------------------------------------------------

      if ( nFilterObjects > 0 ) { 
	
	v_filter_names -> push_back ( name                 );
	v_filter_ids   -> push_back ( std::vector<int  >() );
	v_filter_pts   -> push_back ( std::vector<float>() );
	v_filter_etas  -> push_back ( std::vector<float>() );
	v_filter_phis  -> push_back ( std::vector<float>() ); 
	
	for (int iFilterObject = 0; iFilterObject < nFilterObjects; ++iFilterObject) {

	  int id = int (triggerObjectIds[iFilterObject]);

	  (*v_filter_ids )[nFiltersPassed].push_back ( id ) ;
	  (*v_filter_pts )[nFiltersPassed].push_back ( (float) triggerObjectP4s[iFilterObject].Pt () );
	  (*v_filter_etas)[nFiltersPassed].push_back ( (float) triggerObjectP4s[iFilterObject].Eta() );
	  (*v_filter_phis)[nFiltersPassed].push_back ( (float) triggerObjectP4s[iFilterObject].Phi() );
	}
	
	nFiltersPassed++;
      }
    } // end loop over filters 
  } else { 
    edm::LogError("RootTupleMakerV2_TriggerObjectsError") << "Error! Can't get the product " << inputTag;
  }

  //------------------------------------------------------------------------
  // Push the information into the event
  //------------------------------------------------------------------------

  iEvent.put ( v_filter_names, prefix + "FilterName"    + suffix ) ;
  iEvent.put ( v_filter_ids  , prefix + "FilterObjId"   + suffix ) ;
  iEvent.put ( v_filter_pts  , prefix + "FilterObjPt"   + suffix ) ;
  iEvent.put ( v_filter_etas , prefix + "FilterObjEta"  + suffix ) ;
  iEvent.put ( v_filter_phis , prefix + "FilterObjPhi"  + suffix ) ;
	      

}
