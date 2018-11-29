#include <algorithm>

#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_Trigger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"

unsigned int NmaxL1AlgoBit = 128;
unsigned int NmaxL1TechBit = 64;

RootTupleMakerV2_Trigger::RootTupleMakerV2_Trigger(const edm::ParameterSet& iConfig) :
  hltInputTag_ (iConfig.getParameter<edm::InputTag>("HLTInputTag")),
  l1uGTInputToken_ (consumes<GlobalAlgBlkBxCollection>(iConfig.getParameter<edm::InputTag>("L1uGTInputTag"))),
  hltInputToken_ (consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("HLTInputTag"))),
  packedTrigPrescalesToken_ (consumes<pat::PackedTriggerPrescales> (iConfig.getParameter<edm::InputTag>("PackedPrescalesInputTag"))),
  hltPathsOfInterest(iConfig.getParameter<std::vector<std::string> > ("HLTPathsOfInterest")),
  hltPrescaleProvider_(iConfig, consumesCollector(), *this)
{
  produces <std::string> ("HLTKey");

  produces <std::vector<std::string> > ("HLTInsideDatasetTriggerNames"      );
  produces <std::vector<bool > >       ("HLTInsideDatasetTriggerDecisions"  );
  produces <std::vector<int> >         ("HLTInsideDatasetTriggerPrescales"  );
  produces <std::vector<int> >         ("HLTInsideDatasetTriggerPackedPrescales"  );

  produces <std::vector<int> > ( "L1Bits" );
  produces <int> ( "L1PrescaleColumn" );
}


void RootTupleMakerV2_Trigger::
printNames(const std::vector<std::string>& names) {
  for (unsigned int i = 0; i < names.size(); ++i)
    edm::LogProblem( "RootTupleMakerV2_TriggerProblem" ) << "  " << names[i] << std::endl;
}

void RootTupleMakerV2_Trigger::
beginRun(const edm::Run& iRun, edm::EventSetup const& iSetup) {

  bool changed = true;
  if (hltPrescaleProvider_.init(iRun, iSetup, hltInputTag_.process(), changed)) {
    // if init returns TRUE, initialisation has succeeded!
    edm::LogInfo("RootTupleMakerV2_TriggerInfo") << "HLT config with process name " << hltInputTag_.process() << " successfully extracted";
  } else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    edm::LogError("RootTupleMakerV2_TriggerError") << "Error! HLT config extraction with process name " << hltInputTag_.process() << " failed";
    // In this case, all access methods will return empty values!
  }
}

void RootTupleMakerV2_Trigger::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::unique_ptr<std::vector<int> >  l1bits   ( new std::vector<int>() );
  std::unique_ptr<int> l1prescaleColumn  (new int);
  *l1prescaleColumn.get() = -999;

  std::unique_ptr<std::vector < std::string > > v_hlt_insideDataset_names             (new std::vector<std::string>  ());
  std::unique_ptr<std::vector < std::string > > v_hlt_outsideDataset_names            (new std::vector<std::string>  ());
  std::unique_ptr<std::vector < bool > >        v_hlt_insideDataset_decisions         (new std::vector<bool>         ());
  std::unique_ptr<std::vector < bool > >        v_hlt_outsideDataset_decisions        (new std::vector<bool>         ());
  std::unique_ptr<std::vector < int > >         v_hlt_insideDataset_prescales         (new std::vector<int>          ());
  std::unique_ptr<std::vector < int > >         v_hlt_outsideDataset_prescales        (new std::vector<int>          ());
  std::unique_ptr<std::vector < int > >         v_hlt_insideDataset_packedPrescales   (new std::vector<int>          ());
  std::unique_ptr<std::vector < int > >         v_hlt_outsideDataset_packedPrescales  (new std::vector<int>          ());

  //-----------------------------------------------------------------
  edm::Handle<GlobalAlgBlkBxCollection> l1uGtGlobalAlgBlockBXColl;
  iEvent.getByToken(l1uGTInputToken_, l1uGtGlobalAlgBlockBXColl);

  if(l1uGtGlobalAlgBlockBXColl.isValid()) {
    edm::LogInfo("RootTupleMakerV2_TriggerInfo") << "Successfully obtained " << l1uGtGlobalAlgBlockBXColl;
    GlobalAlgBlk alg = l1uGtGlobalAlgBlockBXColl->at(0,0); // look at BX==0
    int column = l1uGtGlobalAlgBlockBXColl->at(0, 0).getPreScColumn();
    *l1prescaleColumn.get() = column;
    // below will work in later CMSSWs
    //const std::vector<bool> finalDecisions = alg.getAlgoDecisionFinal();
    //for(std::vector<bool>::const_iterator bitItr = finalDecisions.begin(); bitItr != finalDecisions.end(); ++bitItr) {
    //  l1bits->push_back( *bitItr ? 1 : 0 );
    //}
    for(int algBit=0; algBit<300; algBit++) {
      l1bits->push_back( alg.getAlgoDecisionFinal(algBit) ? 1 : 0 );
    }
  } else {
    edm::LogError("RootTupleMakerV2_TriggerError") << "Error! Can't get the GlobalAlgBlkBxCollection!";
  }

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(hltInputToken_, triggerResults);

  edm::Handle<pat::PackedTriggerPrescales> trigPrescales;
  iEvent.getByToken(packedTrigPrescalesToken_,trigPrescales);

  HLTConfigProvider const& hltConfig = hltPrescaleProvider_.hltConfigProvider();

  if(triggerResults.isValid()) {
    edm::LogInfo("RootTupleMakerV2_TriggerInfo") << "Successfully obtained " << triggerResults;

    const edm::TriggerNames& names = iEvent.triggerNames(*triggerResults);
    if(!hltConfig.inited()) {
      throw edm::Exception(edm::errors::LogicError)
          << " HLTConfigProvider was not initialized! Prescales will be wrong!"
          << std::endl;
    }

    //edm::LogError("RootTupleMakerV2_Trigger") << "TRYING to do trigger name matching;" << "hltPathsOfInterest has size: " << hltPathsOfInterest.size();
    // only check paths of interest
    for (auto itr = hltPathsOfInterest.begin(); itr != hltPathsOfInterest.end(); ++itr) {
      for (unsigned int i = 0; i < names.size() ; ++i) {
        const std::string& trigName = names.triggerName(i);
        //edm::LogError("RootTupleMakerV2_Trigger") << "TRY to match " << *itr << " with " << trigName;
        if(trigName.find(*itr)==std::string::npos)
          continue;
        //edm::LogError("RootTupleMakerV2_Trigger") << "Successfully matched " << *itr << " with " << trigName;
        v_hlt_insideDataset_names->push_back ( names.triggerName(i) );
        v_hlt_insideDataset_prescales->push_back ( hltPrescaleProvider_.prescaleValue(iEvent,iSetup,names.triggerName(i)));
        v_hlt_insideDataset_decisions->push_back ( triggerResults->accept(i) );
        v_hlt_insideDataset_packedPrescales->push_back( trigPrescales->getPrescaleForIndex(i) );
      }
    }
    
  } else {
    edm::LogError("RootTupleMakerV2_TriggerError") << "Error! Can't get the product triggerResults";
  }
  
  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put(std::move( l1bits), "L1Bits" );
  iEvent.put(std::move( l1prescaleColumn), "L1PrescaleColumn" );

  iEvent.put(std::move( std::unique_ptr<std::string>(new std::string(hltConfig.tableName()))), "HLTKey");
  
  iEvent.put(std::move( v_hlt_insideDataset_names      ), "HLTInsideDatasetTriggerNames"      ) ;
  iEvent.put(std::move( v_hlt_insideDataset_decisions  ), "HLTInsideDatasetTriggerDecisions"  ) ;
  iEvent.put(std::move( v_hlt_insideDataset_prescales  ), "HLTInsideDatasetTriggerPrescales"  ) ;
  iEvent.put(std::move( v_hlt_insideDataset_packedPrescales  ), "HLTInsideDatasetTriggerPackedPrescales"  ) ;

  
}
 
