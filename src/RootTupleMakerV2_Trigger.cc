#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Trigger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

unsigned int NmaxL1AlgoBit = 128;
unsigned int NmaxL1TechBit = 64;

RootTupleMakerV2_Trigger::RootTupleMakerV2_Trigger(const edm::ParameterSet& iConfig) :
    l1InputTag(iConfig.getParameter<edm::InputTag>("L1InputTag")),
    hltInputTag(iConfig.getParameter<edm::InputTag>("HLTInputTag")),
    hltPathsOfInterest(iConfig.getParameter<std::vector<std::string> > ("HLTPathsOfInterest"))
{
  produces <std::vector<int> > ( "L1PhysBits" );
  produces <std::vector<int> > ( "L1TechBits" );
  produces <std::vector<int> > ( "HLTBits" );
  produces <std::vector<int> > ( "HLTResults" );
}

void RootTupleMakerV2_Trigger::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<int> >  l1physbits   ( new std::vector<int>() );
  std::auto_ptr<std::vector<int> >  l1techbits   ( new std::vector<int>() );
  std::auto_ptr<std::vector<int> >  hltbits      ( new std::vector<int>() );
  std::auto_ptr<std::vector<int> >  hltresults   ( new std::vector<int>() );

  //-----------------------------------------------------------------
  edm::Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  iEvent.getByLabel(l1InputTag, l1GtReadoutRecord);

  if(l1GtReadoutRecord.isValid()) {
    edm::LogInfo("RootTupleMakerV2_TriggerInfo") << "Successfully obtained " << l1InputTag;

    for (unsigned int i = 0; i < NmaxL1AlgoBit; ++i) {
      l1physbits->push_back( l1GtReadoutRecord->decisionWord()[i] ? 1 : 0 );
    }
    for (unsigned int i = 0; i < NmaxL1TechBit; ++i) {
      l1techbits->push_back( l1GtReadoutRecord->technicalTriggerWord()[i] ? 1 : 0 );
    }
  } else {
    edm::LogError("RootTupleMakerV2_TriggerError") << "Error! Can't get the product " << l1InputTag;
  }

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByLabel(hltInputTag, triggerResults);

  if(triggerResults.isValid()) {
    edm::LogInfo("RootTupleMakerV2_TriggerInfo") << "Successfully obtained " << hltInputTag;

    for( unsigned int i = 0; i < triggerResults->size(); i++ ){
      hltbits->push_back( triggerResults->at(i).accept() ? 1 : 0 );
    }

    const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);

    for( std::vector<std::string>::const_iterator it = hltPathsOfInterest.begin();
         it != hltPathsOfInterest.end(); ++it ) {
      int fired = 0;
      unsigned int index = triggerNames.triggerIndex(*it);
      if( index < triggerResults->size() ) {
        if( triggerResults->accept( index ) ) fired = 1;
      } else {
        edm::LogInfo("RootTupleMakerV2_TriggerInfo") << "Requested HLT path \"" << (*it) << "\" does not exist";
      }
      hltresults->push_back( fired );
    }
  } else {
    edm::LogError("RootTupleMakerV2_TriggerError") << "Error! Can't get the product " << hltInputTag;
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( l1physbits, "L1PhysBits" );
  iEvent.put( l1techbits, "L1TechBits" );
  iEvent.put( hltbits,    "HLTBits" );
  iEvent.put( hltresults, "HLTResults" );
}
