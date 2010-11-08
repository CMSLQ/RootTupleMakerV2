#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_MET.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/MET.h"

RootTupleMakerV2_MET::RootTupleMakerV2_MET(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    store_uncorrected_MET (iConfig.getParameter<bool>  ("StoreUncorrectedMET"))
{
  produces <std::vector<double> > ( prefix + "MET" + suffix );
  produces <std::vector<double> > ( prefix + "METPhi" + suffix );
  produces <std::vector<double> > ( prefix + "SumET" + suffix );
  if ( store_uncorrected_MET ) {
    produces <std::vector<double> > ( prefix + "METUncorr" + suffix );
    produces <std::vector<double> > ( prefix + "METPhiUncorr" + suffix );
    produces <std::vector<double> > ( prefix + "SumETUncorr" + suffix );
  }
}

void RootTupleMakerV2_MET::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<double> >  met  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  metphi  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  sumet  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  metuncorr  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  metphiuncorr  ( new std::vector<double>()  );
  std::auto_ptr<std::vector<double> >  sumetuncorr  ( new std::vector<double>()  );
  
  //-----------------------------------------------------------------
  edm::Handle<std::vector<pat::MET> > mets;
  iEvent.getByLabel(inputTag, mets);

  if(mets.isValid()) {
    edm::LogInfo("RootTupleMakerV2_METInfo") << "Total # METs: " << mets->size();

    for( std::vector<pat::MET>::const_iterator it = mets->begin(); it != mets->end(); ++it ) {

      // fill in all the vectors
      met->push_back( it->pt() );
      metphi->push_back( it->phi() );
      sumet->push_back( it->sumEt() );
      
      if ( store_uncorrected_MET ) {
	metuncorr->push_back( it->uncorrectedPt(pat::MET::uncorrALL) );
	metphiuncorr->push_back( it->uncorrectedPhi(pat::MET::uncorrALL) );
	sumetuncorr->push_back( it->sumEt() - it->corSumEt(pat::MET::uncorrALL) );
      }

    }
  } else {
    edm::LogError("RootTupleMakerV2_METError") << "Error! Can't get the product " << inputTag;
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( met, prefix + "MET" + suffix );
  iEvent.put( metphi, prefix + "METPhi" + suffix );
  iEvent.put( sumet, prefix + "SumET" + suffix );
  if ( store_uncorrected_MET ) {
    iEvent.put( metuncorr, prefix + "METUncorr" + suffix );
    iEvent.put( metphiuncorr, prefix + "METPhiUncorr" + suffix );
    iEvent.put( sumetuncorr, prefix + "SumETUncorr" + suffix );
  }
}
