#include <algorithm>

#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_Trigger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/Common/interface/TriggerResults.h"


unsigned int NmaxL1AlgoBit = 128;
unsigned int NmaxL1TechBit = 64;

RootTupleMakerV2_Trigger::RootTupleMakerV2_Trigger(const edm::ParameterSet& iConfig) :
  hltInputTag_ (iConfig.getParameter<edm::InputTag>("HLTInputTag")),
  l1uGTInputToken_ (consumes<GlobalAlgBlkBxCollection>(iConfig.getParameter<edm::InputTag>("L1uGTInputTag"))),
  hltInputToken_ (consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("HLTInputTag"))),
  packedTrigPrescalesToken_ (consumes<pat::PackedTriggerPrescales> (iConfig.getParameter<edm::InputTag>("PackedPrescalesInputTag"))),
  hltPathsOfInterest(iConfig.getParameter<std::vector<std::string> > ("HLTPathsOfInterest")),
  hltPrescaleProvider_(iConfig, consumesCollector(), *this),
  sourceName(iConfig.getParameter<std::string>  ("SourceName")),
  sourceType(NOT_APPLICABLE)
{
  // Source is either a stream or a dataset (mutually exclusive)
  if (sourceName.length() > 0) {
    if (sourceName.length() >= 2 && sourceName[0]=='S' && sourceName[1]==':') {
      sourceType      = STREAM;
      sourceName      = sourceName.substr(2);
    }
    else if (sourceName.length() >= 3 && sourceName[0]=='D' && sourceName[1]=='S' && sourceName[2]==':') {
      sourceType      = DATASET;
      sourceName      = sourceName.substr(3);
    }
    else throw edm::Exception(edm::errors::Configuration)
	   << "Invalid SourceName = '" << sourceName 
	   << "' -- must start with either 'S:' for streams or 'DS:' for datasets."
	   << std::endl;
  }
  
  produces <std::string> ("HLTKey");

  produces <std::vector<std::string> > ("HLTInsideDatasetTriggerNames"      );
  produces <std::vector<std::string> > ("HLTOutsideDatasetTriggerNames"     );
  produces <std::vector<bool > >       ("HLTInsideDatasetTriggerDecisions"  );
  produces <std::vector<bool > >       ("HLTOutsideDatasetTriggerDecisions" );
  produces <std::vector<int> >         ("HLTInsideDatasetTriggerPrescales"  );
  produces <std::vector<int> >         ("HLTOutsideDatasetTriggerPrescales" );
  produces <std::vector<int> >         ("HLTInsideDatasetTriggerPackedPrescales"  );
  produces <std::vector<int> >         ("HLTOutsideDatasetTriggerPackedPrescales" );

  produces <std::vector<int> > ( "L1Bits" );
  produces <std::vector<int> > ( "L1PrescaleColumn" );
  //produces <std::vector<int> > ( "L1PhysBits" );
  //produces <std::vector<int> > ( "L1TechBits" );
}


void RootTupleMakerV2_Trigger::
printNames(const std::vector<std::string>& names) {
  for (unsigned int i = 0; i < names.size(); ++i)
    edm::LogProblem( "RootTupleMakerV2_TriggerProblem" ) << "  " << names[i] << std::endl;
}


void RootTupleMakerV2_Trigger::
getDataSource() {
  dataSource.clear();
  if (sourceType == NOT_APPLICABLE) return;

  HLTConfigProvider const& hltConfig = hltPrescaleProvider_.hltConfigProvider();
  
  if (sourceType == STREAM) {
    unsigned int  index   = hltConfig.streamIndex(sourceName);
    if (index >= hltConfig.streamNames().size()) {
      edm::LogError( "RootTupleMakerV2_TriggerError" ) << "Streams in '" << hltInputTag_.process() << "' HLT menu:";
      printNames(hltConfig.streamNames());
      throw edm::Exception(edm::errors::Configuration) << "Stream with name '" << sourceName << "' does not exist." << std::endl;
    }
    dataSource    = hltConfig.streamContent(sourceName);
  }
  else {
    unsigned int  index   = hltConfig.datasetIndex(sourceName);
    if (index >= hltConfig.datasetNames().size()) {
      edm::LogError( "RootTupleMakerV2_TriggerError" ) << "Datasets in '" << hltInputTag_.process() << "' HLT menu:";
      printNames(hltConfig.datasetNames());
      throw edm::Exception(edm::errors::Configuration) << "Dataset with name '" << sourceName << "' does not exist." << std::endl;
    }
    dataSource    = hltConfig.datasetContent(sourceName);
  }
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

  getDataSource();
}

void RootTupleMakerV2_Trigger::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<int> >  l1bits   ( new std::vector<int>() );
  std::auto_ptr<std::vector<int> >  l1prescaleColumn   ( new std::vector<int>() );
  //std::auto_ptr<std::vector<int> >  l1physbits   ( new std::vector<int>() );
  //std::auto_ptr<std::vector<int> >  l1techbits   ( new std::vector<int>() );
  // std::auto_ptr<std::vector<int> >  hltbits      ( new std::vector<int>() );
  // std::auto_ptr<std::vector<int> >  hltresults   ( new std::vector<int>() );
  // std::auto_ptr<std::vector<int> >  hltprescales ( new std::vector<int>() );

  std::auto_ptr<std::vector < std::string > > v_hlt_insideDataset_names             (new std::vector<std::string>  ());
  std::auto_ptr<std::vector < std::string > > v_hlt_outsideDataset_names            (new std::vector<std::string>  ());
  std::auto_ptr<std::vector < bool > >        v_hlt_insideDataset_decisions         (new std::vector<bool>         ());
  std::auto_ptr<std::vector < bool > >        v_hlt_outsideDataset_decisions        (new std::vector<bool>         ());
  std::auto_ptr<std::vector < int > >         v_hlt_insideDataset_prescales         (new std::vector<int>          ());
  std::auto_ptr<std::vector < int > >         v_hlt_outsideDataset_prescales        (new std::vector<int>          ());
  std::auto_ptr<std::vector < int > >         v_hlt_insideDataset_packedPrescales   (new std::vector<int>          ());
  std::auto_ptr<std::vector < int > >         v_hlt_outsideDataset_packedPrescales  (new std::vector<int>          ());

  /*
  std::auto_ptr<std::map<std::string,bool > > m_hlt_insideDataset_namesToDecisions  (new std::map<std::string,bool>());
  std::auto_ptr<std::map<std::string,bool > > m_hlt_outsideDataset_namesToDecisions (new std::map<std::string,bool>());
  std::auto_ptr<std::map<std::string,int  > > m_hlt_insideDataset_namesToPrescales  (new std::map<std::string,int >());
  std::auto_ptr<std::map<std::string,int  > > m_hlt_outsideDataset_namesToPrescales (new std::map<std::string,int >());
  */

  //-----------------------------------------------------------------
  edm::Handle<GlobalAlgBlkBxCollection> l1uGtGlobalAlgBlockBXColl;
  iEvent.getByToken(l1uGTInputToken_, l1uGtGlobalAlgBlockBXColl);

  if(l1uGtGlobalAlgBlockBXColl.isValid()) {
    edm::LogInfo("RootTupleMakerV2_TriggerInfo") << "Successfully obtained " << l1uGtGlobalAlgBlockBXColl;
    GlobalAlgBlk alg = l1uGtGlobalAlgBlockBXColl->at(0,0); // look at BX==0
    int column = l1uGtGlobalAlgBlockBXColl->at(0, 0).getPreScColumn();
    l1prescaleColumn->push_back(column);
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
    edm::LogInfo("RootTupleMakerV2_TriggerInfo") << "Successfully obtained " << triggerResults;//fixme what to put here?

    const edm::TriggerNames& names = iEvent.triggerNames(*triggerResults);
    if(!hltConfig.inited()) {
      throw edm::Exception(edm::errors::LogicError)
          << " HLTConfigProvider was not initialized! Prescales will be wrong!"
          << std::endl;
    }

    for (unsigned int i = 0; i < names.size() ; ++i) { 
      if (dataSource.empty() || std::find(dataSource.begin(), dataSource.end(), names.triggerName(i)) != dataSource.end()) {
        v_hlt_insideDataset_names->push_back ( names.triggerName(i) );
        v_hlt_insideDataset_prescales->push_back ( hltPrescaleProvider_.prescaleValue(iEvent,iSetup,names.triggerName(i)));
        v_hlt_insideDataset_decisions->push_back ( triggerResults->accept(i) );
        v_hlt_insideDataset_packedPrescales->push_back( trigPrescales->getPrescaleForIndex(i) );
      } else {
        v_hlt_outsideDataset_names->push_back ( names.triggerName(i) );
        v_hlt_outsideDataset_prescales->push_back ( hltPrescaleProvider_.prescaleValue(iEvent,iSetup,names.triggerName(i)));
        v_hlt_outsideDataset_decisions->push_back ( triggerResults->accept(i) );
        v_hlt_outsideDataset_packedPrescales->push_back( trigPrescales->getPrescaleForIndex(i) );
      }      
    }
    

    /*
    for( unsigned int i = 0; i < triggerResults->size(); i++ ){
      hltbits->push_back( triggerResults->at(i).accept() ? 1 : 0 );
    }

    for( std::vector<std::string>::const_iterator it = hltPathsOfInterest.begin();
         it != hltPathsOfInterest.end(); ++it ) {
      int fired = 0;
      unsigned int index = hltConfig.triggerIndex(*it);
      if( index < triggerResults->size() ) {
        if( triggerResults->accept( index ) ) fired = 1;
      } else {
        edm::LogInfo("RootTupleMakerV2_TriggerInfo") << "Requested HLT path \"" << (*it) << "\" does not exist";
      }
      hltresults->push_back( fired );

      int prescale = -1;
      if(hltConfig.prescaleSet(iEvent, iSetup)<0) {
        edm::LogError("RootTupleMakerV2_TriggerError") << "Error! The prescale set index number could not be obtained";
      } else {
        prescale = hltConfig.prescaleValue(iEvent, iSetup, *it);
      }
      hltprescales->push_back( prescale );

    }
    */
  } else {
    edm::LogError("RootTupleMakerV2_TriggerError") << "Error! Can't get the product triggerResults";
  }
  
  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( l1bits, "L1Bits" );
  iEvent.put( l1prescaleColumn, "L1PrescaleColumn" );
  //iEvent.put( l1physbits, "L1PhysBits" );
  //iEvent.put( l1techbits, "L1TechBits" );
  // iEvent.put( hltbits,    "HLTBits" );
  // iEvent.put( hltresults, "HLTResults" );
  // iEvent.put( hltprescales, "HLTPrescales" );

  iEvent.put( std::auto_ptr<std::string>(new std::string(hltConfig.tableName())), "HLTKey");
  
  /*
  iEvent.put( m_hlt_insideDataset_namesToDecisions ,"HLTInsideDatasetDecisionMap" );
  iEvent.put( m_hlt_outsideDataset_namesToDecisions,"HLTOutsideDatasetDecisionMap");
  iEvent.put( m_hlt_insideDataset_namesToPrescales ,"HLTInsideDatasetPrescaleMap" );
  iEvent.put( m_hlt_outsideDataset_namesToPrescales,"HLTOutsideDatasetPrescaleMap");
  */

  iEvent.put ( v_hlt_insideDataset_names      , "HLTInsideDatasetTriggerNames"      ) ;
  iEvent.put ( v_hlt_outsideDataset_names     , "HLTOutsideDatasetTriggerNames"     ) ;
  iEvent.put ( v_hlt_insideDataset_decisions  , "HLTInsideDatasetTriggerDecisions"  ) ;
  iEvent.put ( v_hlt_outsideDataset_decisions , "HLTOutsideDatasetTriggerDecisions" ) ;
  iEvent.put ( v_hlt_insideDataset_prescales  , "HLTInsideDatasetTriggerPrescales"  ) ;
  iEvent.put ( v_hlt_outsideDataset_prescales , "HLTOutsideDatasetTriggerPrescales" ) ;
  iEvent.put ( v_hlt_insideDataset_packedPrescales  , "HLTInsideDatasetTriggerPackedPrescales"  ) ;
  iEvent.put ( v_hlt_outsideDataset_packedPrescales , "HLTOutsideDatasetTriggerPackedPrescales" ) ;

  
}
 
