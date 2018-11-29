#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_TriggerObjects.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "TLorentzVector.h"

RootTupleMakerV2_TriggerObjects::RootTupleMakerV2_TriggerObjects(const edm::ParameterSet& iConfig) :
  triggerBitsToken_    (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerBitsTag"))),
  triggerObjectsToken_ (consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("TriggerObjectsTag"))),
  prefix     (iConfig.getParameter<std::string>  ("Prefix")),
  suffix     (iConfig.getParameter<std::string>  ("Suffix")),
  hltPathsOfInterest(iConfig.getParameter<std::vector<std::string> > ("HLTPathsOfInterest"))
{
  //------------------------------------------------------------------------
  // What variables will this producer push into the event?
  //------------------------------------------------------------------------
  
  produces <std::vector<std::vector< std::string > > > ( prefix + "TriggerObjFilterNames"   + suffix );
  produces <std::vector<std::vector< std::string > > > ( prefix + "TriggerObjPathNames"   + suffix );
  produces <std::vector<std::vector< bool        > > > ( prefix + "TriggerObjPassedPathL3Filter"  + suffix );
  produces <std::vector<std::vector< bool        > > > ( prefix + "TriggerObjPassedPathLastFilter"  + suffix );
  produces <std::vector<std::vector< int         > > > ( prefix + "TriggerObjTypeIds"  + suffix );
  produces <std::vector< float                     > > ( prefix + "TriggerObjPt"  + suffix );
  produces <std::vector< float                     > > ( prefix + "TriggerObjEta" + suffix );
  produces <std::vector< float                     > > ( prefix + "TriggerObjPhi" + suffix );
  produces <std::vector< std::string               > > ( prefix + "TriggerObjCollectionName" + suffix );

}

void RootTupleMakerV2_TriggerObjects::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  //------------------------------------------------------------------------
  // Declare items to push into the event
  //------------------------------------------------------------------------

  std::unique_ptr<std::vector< std::vector< std::string > > > v_filter_names    ( new std::vector< std::vector< std::string > > ());
  std::unique_ptr<std::vector< std::vector< std::string > > > v_path_names      ( new std::vector< std::vector< std::string > > ());
  std::unique_ptr<std::vector< std::vector< bool        > > > v_passed_path_l3  ( new std::vector< std::vector< bool        > > ());
  std::unique_ptr<std::vector< std::vector< bool        > > > v_passed_path_last( new std::vector< std::vector< bool        > > ());
  std::unique_ptr<std::vector< std::vector< int         > > > v_type_ids        ( new std::vector< std::vector< int         > > ());
  std::unique_ptr<std::vector< float                      > > v_pt              ( new std::vector< float                      > ());
  std::unique_ptr<std::vector< float                      > > v_eta             ( new std::vector< float                      > ());
  std::unique_ptr<std::vector< float                      > > v_phi             ( new std::vector< float                      > ());
  std::unique_ptr<std::vector< std::string                > > v_collection      ( new std::vector< std::string                > ());

  //------------------------------------------------------------------------
  // Get the trigger bits and make sure they are valid
  //------------------------------------------------------------------------

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken( triggerBitsToken_ , triggerBits);
  if ( ! triggerBits.isValid() )
    edm::LogError("RootTupleMakerV2_TriggerObjectsError") << "Error! Can't get the trigger bits";

  //------------------------------------------------------------------------
  // Get the trigger objects and make sure they are valid
  //------------------------------------------------------------------------
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken( triggerObjectsToken_ , triggerObjects);
  if ( ! triggerObjects.isValid() )
    edm::LogError("RootTupleMakerV2_TriggerObjectsError") << "Error! Can't get the trigger objects";

  //edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

  //------------------------------------------------------------------------
  // Get the trigger names
  //------------------------------------------------------------------------
  //const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);

  int nTrigObjects = 0;
  //------------------------------------------------------------------------
  // Loop over trigger objects
  //------------------------------------------------------------------------
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    obj.unpackNamesAndLabels(iEvent,*triggerBits);
    std::vector<std::string> pathNamesAll  = obj.pathNames(false);
    // match with paths of interest
    std::vector<std::string> matchedPathNames;
    for (auto itr = hltPathsOfInterest.begin(); itr != hltPathsOfInterest.end(); ++itr) {
      for (auto objPathItr = pathNamesAll.begin(); objPathItr != pathNamesAll.end(); ++objPathItr) {
        if(objPathItr->find(*itr)!=std::string::npos)
          matchedPathNames.push_back(*objPathItr);
      }
    }
    if(matchedPathNames.empty())
      continue;

    v_path_names -> push_back(matchedPathNames);
    // Record if the object is associated to a 'l3' filter (always true for the definition used in the PAT trigger producer)
    //   and if it's associated to the last filter of a successful path,
    //   which means that this object did cause this trigger to succeed. But it doesn't work on some multi-object triggers.
    v_passed_path_l3 -> push_back(std::vector<bool>(matchedPathNames.size(),true));
    v_passed_path_last -> push_back(std::vector<bool>());
    for (unsigned h = 0, n = matchedPathNames.size(); h < n; ++h) {
      bool isLF   = obj.hasPathName( matchedPathNames[h], true, false ); 
      (*v_passed_path_last )[nTrigObjects].push_back ( isLF );
    }
    //std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
    v_pt->push_back(obj.pt());
    v_eta->push_back(obj.eta());
    v_phi->push_back(obj.phi());
    //// Print trigger object collection and type
    //std::cout << "\t   Collection: " << obj.collection() << std::endl;
    v_collection->push_back(obj.collection());
    //std::cout << "\t   Type IDs:   ";
    //for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
    v_type_ids   -> push_back ( std::vector<int>() );
    for (unsigned h = 0; h < obj.triggerObjectTypes().size(); ++h)
      (*v_type_ids )[nTrigObjects].push_back ( obj.triggerObjectTypes()[h] );
    //std::cout << std::endl;
    //// Print associated trigger filters
    //std::cout << "\t   Filters:    ";
    //for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
    v_filter_names -> push_back (std::vector<std::string>() );
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h)
      (*v_filter_names )[nTrigObjects].push_back ( obj.filterLabels()[h] );

    nTrigObjects++;

    //std::cout << std::endl;
    //std::vector<std::string> pathNamesAll  = obj.pathNames(false);
    //std::vector<std::string> pathNamesLast = obj.pathNames(true);
    //// Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
    //// definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
    //// means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
    //std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
    //for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
    //  bool isBoth = obj.hasPathName( pathNamesAll[h], true, true ); 
    //  bool isL3   = obj.hasPathName( pathNamesAll[h], false, true ); 
    //  bool isLF   = obj.hasPathName( pathNamesAll[h], true, false ); 
    //  bool isNone = obj.hasPathName( pathNamesAll[h], false, false ); 
    //  std::cout << "   " << pathNamesAll[h];
    //  if (isBoth) std::cout << "(L,3)";
    //  if (isL3 && !isBoth) std::cout << "(*,3)";
    //  if (isLF && !isBoth) std::cout << "(L,*)";
    //  if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
    //}
    //std::cout << std::endl;
  }
  //std::cout << std::endl;



  //const trigger::TriggerObjectCollection & triggerObjects = triggerEvent -> getObjects();

  ////------------------------------------------------------------------------
  //// Loop over the filters in the trigger event
  //// e.g.: hltEle27WP80TrackIsoFilter
  ////------------------------------------------------------------------------

  //size_t nFilters       = triggerEvent -> sizeFilters();
  //size_t nFiltersPassed = 0;
  //size_t iFilter        = 0;

  //for (; iFilter < nFilters; ++iFilter) {

  //  //------------------------------------------------------------------------
  //  // Find information for each filter:
  //  // -  Name  : name of the filter, e.g. hltEle27WP80TrackIsoFilter
  //  // - "Keys" : std::vector<uint16_t> storing indices of trigger objects that pass the filter
  //  //------------------------------------------------------------------------

  //  std::string          name = triggerEvent -> filterTag ( iFilter ).label();
  //  const trigger::Keys& keys = triggerEvent -> filterKeys( iFilter );
  //  const trigger::Vids& vids = triggerEvent -> filterIds ( iFilter );

  //  //------------------------------------------------------------------------
  //  // Loop over the keys to get to the trigger objects that pass the filter
  //  //------------------------------------------------------------------------

  //  int nKeys = (int) keys.size();
  //  int nVids = (int) vids.size();
  //  assert ( nKeys == nVids ) ;

  //  // useful variables
  //  std::vector<TLorentzVector> triggerObjectP4s;
  //  std::vector<int>            triggerObjectIds;

  //  for (int iTriggerObject = 0; iTriggerObject < nKeys; ++iTriggerObject ) { 

  //    // Get the object ID and key

  //    int                id  = vids[iTriggerObject];
  //    trigger::size_type key = keys[iTriggerObject];

  //    // Get the trigger object from the key

  //    const trigger::TriggerObject & triggerObject = triggerObjects [key];

  //    // Store the trigger object as a TLorentzVector (borrowed from S. Harper)

  //    TLorentzVector p4;
  //    p4.SetPtEtaPhiM ( triggerObject.pt  (),
  //        triggerObject.eta (),
  //        triggerObject.phi (),
  //        triggerObject.mass() );

  //    triggerObjectP4s.push_back ( p4 ) ;
  //    triggerObjectIds.push_back ( id ) ;

  //  } // end loop over keys/trigger objects passing filters

  //  //------------------------------------------------------------------------
  //  // If the filter passed, store its information
  //  //------------------------------------------------------------------------

  //  if ( nKeys > 0 ) { 

  //    v_filter_names -> push_back ( name                 );
  //    v_filter_ids   -> push_back ( std::vector<int  >() );
  //    v_filter_pts   -> push_back ( std::vector<float>() );
  //    v_filter_etas  -> push_back ( std::vector<float>() );
  //    v_filter_phis  -> push_back ( std::vector<float>() ); 

  //    for (int iFilterObject = 0; iFilterObject < nKeys; ++iFilterObject) {

  //      int id = int (triggerObjectIds[iFilterObject]);

  //      (*v_filter_ids )[nFiltersPassed].push_back ( id ) ;
  //      // Some trigger objects fail with the following message:
  //      // ------------------------------------------------//
  //      // Fatal Root Error: @SUB=TVector3::PseudoRapidity //
  //      // transvers momentum = 0! return +/- 10e10        //
  //      // ------------------------------------------------//
  //      // Ex/   hltElectron40CaloIdTTrkIdTCleanedPFHT300
  //      // Ex/   hltElectron60CaloIdTTrkIdTCleanedPFHT300
  //      // Hence, set Pt,Eta,Phi to -99 if Pt=0   
  //      if( (float)(triggerObjectP4s[iFilterObject].Pt())>0 )
  //      {
  //        (*v_filter_pts )[nFiltersPassed].push_back ( (float) triggerObjectP4s[iFilterObject].Pt () );
  //        (*v_filter_etas)[nFiltersPassed].push_back ( (float) triggerObjectP4s[iFilterObject].Eta() );
  //        (*v_filter_phis)[nFiltersPassed].push_back ( (float) triggerObjectP4s[iFilterObject].Phi() );
  //      }
  //      else
  //      {
  //        (*v_filter_pts )[nFiltersPassed].push_back ( (float)(-99) );
  //        (*v_filter_etas)[nFiltersPassed].push_back ( (float)(-99) );
  //        (*v_filter_phis)[nFiltersPassed].push_back ( (float)(-99) );
  //      }
  //      //
  //    }

  //    nFiltersPassed++;
  //  }
  //} // end loop over filters 

  //------------------------------------------------------------------------
  // Push the information into the event
  //------------------------------------------------------------------------

  iEvent.put(std::move( v_filter_names), prefix + "TriggerObjFilterNames"    + suffix ) ;
  iEvent.put(std::move( v_path_names), prefix + "TriggerObjPathNames" + suffix);
  iEvent.put(std::move( v_passed_path_l3), prefix + "TriggerObjPassedPathL3Filter" + suffix);
  iEvent.put(std::move( v_passed_path_last), prefix + "TriggerObjPassedPathLastFilter" + suffix);
  iEvent.put(std::move( v_type_ids ), prefix + "TriggerObjTypeIds"   + suffix ) ;
  iEvent.put(std::move( v_pt ), prefix + "TriggerObjPt"   + suffix ) ;
  iEvent.put(std::move( v_eta ), prefix + "TriggerObjEta"  + suffix ) ;
  iEvent.put(std::move( v_phi ), prefix + "TriggerObjPhi"  + suffix ) ;
  iEvent.put(std::move( v_collection), prefix + "TriggerObjCollectionName" + suffix) ;

}
