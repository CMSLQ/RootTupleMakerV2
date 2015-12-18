#include "Leptoquarks/RootTupleMakerV2/plugins/RootTupleMakerV2_GenMET.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/MET.h"


RootTupleMakerV2_GenMET::RootTupleMakerV2_GenMET(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix"))
{
  produces <std::vector<double> > ( prefix + "MET" + suffix );
  produces <std::vector<double> > ( prefix + "METPhi" + suffix );
  produces <std::vector<double> > ( prefix + "SumET" + suffix );
}

void RootTupleMakerV2_GenMET::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  met  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  metphi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  sumet  ( new std::vector<double>()  );

  //-----------------------------------------------------------------
  if( !iEvent.isRealData() ) {
    edm::Handle<std::vector<pat::MET> > mets;
    iEvent.getByLabel(inputTag, mets);

    if(mets.isValid()) {
      edm::LogInfo("RootTupleMakerV2_GenMETInfo") << "Total # GenMETs: " << mets->size();

    for( std::vector<pat::MET>::const_iterator it = mets->begin(); it != mets->end(); ++it ) {

        // fill in all the vectors
        met->push_back( it->genMET()->pt() );
        metphi->push_back( it->genMET()->phi() );
        sumet->push_back( it->genMET()->sumEt() );
      }
    } else {
      edm::LogError("RootTupleMakerV2_GenMETError") << "Error! Can't get the product " << inputTag;
    }
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( met, prefix + "MET" + suffix );
  iEvent.put( metphi, prefix + "METPhi" + suffix );
  iEvent.put( sumet, prefix + "SumET" + suffix );
}
