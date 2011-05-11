#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_PileUp.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

using namespace std;

RootTupleMakerV2_PileUp::RootTupleMakerV2_PileUp(const edm::ParameterSet& iConfig):

  pileupInfoSrc(iConfig.getParameter<edm::InputTag>("pileupInfo"))
  {
    produces <std::vector<int> > ( "PileUpInteractions"   ); 
  }

void RootTupleMakerV2_PileUp::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<int> >  Number_interactions ( new std::vector<int>()  );
  edm::Handle<std::vector<PileupSummaryInfo> >  puInfo;
  iEvent.getByLabel(pileupInfoSrc, puInfo);
  
  if(puInfo.isValid()) {
    for( std::vector<PileupSummaryInfo>::const_iterator it = puInfo->begin(); it != puInfo->end(); ++it ) {
      Number_interactions ->push_back( it->getPU_NumInteractions() );
      std::cout<<it->getPU_NumInteractions()<<std::endl;
    }
  }
  else {
    edm::LogError("RootTupleMakerV2_PileUpError") << "Error! Can't get the product " << pileupInfoSrc;
  }
  iEvent.put( Number_interactions,   "PileUpInteractions"   );
}
