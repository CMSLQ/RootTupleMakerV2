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
  
  produces <std::vector<std::vector< std::string > > > ( prefix + "TriggerObjPathNames"   + suffix );
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

  std::unique_ptr<std::vector< std::vector< std::string > > > v_path_names      ( new std::vector< std::vector< std::string > > ());
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

  bool printStuff = false;
  int nTrigObjects = 0;
  //------------------------------------------------------------------------
  // Loop over trigger objects
  //------------------------------------------------------------------------
  if(printStuff)
    std::cout << std::endl;
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
    // Don't record if the object is associated to a 'l3' filter (always true for the definition used in the PAT trigger producer)
    // Record if it's associated to the last filter of a successful path,
    //   which means that this object did cause this trigger to succeed. But it doesn't work on some multi-object triggers.
    v_passed_path_last -> push_back(std::vector<bool>());
    for (unsigned h = 0, n = matchedPathNames.size(); h < n; ++h) {
      bool isLF   = obj.hasPathName( matchedPathNames[h], true, false ); 
      (*v_passed_path_last )[nTrigObjects].push_back ( isLF );
    }
    v_pt->push_back(obj.pt());
    v_eta->push_back(obj.eta());
    v_phi->push_back(obj.phi());
    v_collection->push_back(obj.collection());
    v_type_ids   -> push_back ( std::vector<int>() );
    for (unsigned h = 0; h < obj.triggerObjectTypes().size(); ++h)
      (*v_type_ids )[nTrigObjects].push_back ( obj.triggerObjectTypes()[h] );

    nTrigObjects++;

    if(printStuff) {
      std::cout << std::endl;
      std::cout << "---->Trigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      // Print trigger object collection and type
      std::cout << "\t   Collection: " << obj.collection() << std::endl;
      std::cout << "\t   Type IDs:   ";
      for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
      // Print associated trigger filters
      std::cout << std::endl;
      std::cout << "\t   Filters:    ";
      for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
      // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
      // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
      // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
      std::vector<std::string> pathNamesAll  = obj.pathNames(false);
      std::vector<std::string> pathNamesLast = obj.pathNames(true);
      std::cout << "\t   Paths (all: " << pathNamesAll.size()<<" last: "<<pathNamesLast.size()<<"):    ";
      for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
        bool isBoth = obj.hasPathName( pathNamesAll[h], true, true ); 
        bool isL3   = obj.hasPathName( pathNamesAll[h], false, true ); 
        bool isLF   = obj.hasPathName( pathNamesAll[h], true, false ); 
        bool isNone = obj.hasPathName( pathNamesAll[h], false, false ); 
        std::cout << "   " << pathNamesAll[h];
        if (isBoth) std::cout << "(L,3)";
        if (isL3 && !isBoth) std::cout << "(X,3)";
        if (isLF && !isBoth) std::cout << "(L,X)";
        if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(X,X)";
      }
      std::cout << std::endl;
    }
  }
  if(printStuff)
    std::cout << std::endl;

  //------------------------------------------------------------------------
  // Push the information into the event
  //------------------------------------------------------------------------

  iEvent.put(std::move( v_path_names), prefix + "TriggerObjPathNames" + suffix);
  iEvent.put(std::move( v_passed_path_last), prefix + "TriggerObjPassedPathLastFilter" + suffix);
  iEvent.put(std::move( v_type_ids ), prefix + "TriggerObjTypeIds"   + suffix ) ;
  iEvent.put(std::move( v_pt ), prefix + "TriggerObjPt"   + suffix ) ;
  iEvent.put(std::move( v_eta ), prefix + "TriggerObjEta"  + suffix ) ;
  iEvent.put(std::move( v_phi ), prefix + "TriggerObjPhi"  + suffix ) ;
  iEvent.put(std::move( v_collection), prefix + "TriggerObjCollectionName" + suffix) ;

}
