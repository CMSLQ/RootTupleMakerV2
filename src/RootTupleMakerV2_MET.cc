#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_MET.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/MET.h"

RootTupleMakerV2_MET::RootTupleMakerV2_MET(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    store_uncorrected_MET (iConfig.getParameter<bool>  ("StoreUncorrectedMET")),
    store_MET_significance (iConfig.getParameter<bool>  ("StoreMETSignificance")),
    uncertainty (iConfig.getParameter<std::string>  ("Uncertainty")),
    level       (iConfig.getParameter<std::string>  ("Level"))
{
  produces <std::vector<double> > ( prefix + "MET" + suffix );
  produces <std::vector<double> > ( prefix + "METPhi" + suffix );
  produces <std::vector<double> > ( prefix + "SumET" + suffix );
  if ( store_uncorrected_MET ) {
    produces <std::vector<double> > ( prefix + "METUncorr" + suffix );
    produces <std::vector<double> > ( prefix + "METPhiUncorr" + suffix );
    produces <std::vector<double> > ( prefix + "SumETUncorr" + suffix );
  }
  if ( store_MET_significance ) {
    produces <std::vector<double> > ( prefix + "METSig" + suffix );
    produces <std::vector<double> > ( prefix + "METSigMatrixDXX" + suffix );
    produces <std::vector<double> > ( prefix + "METSigMatrixDXY" + suffix );
    produces <std::vector<double> > ( prefix + "METSigMatrixDYX" + suffix );
    produces <std::vector<double> > ( prefix + "METSigMatrixDYY" + suffix );
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
  std::auto_ptr<std::vector<double> >  metsig  ( new std::vector<double>()  );  
  std::auto_ptr<std::vector<double> >  metsigmatrixdxx  ( new std::vector<double>()  );  
  std::auto_ptr<std::vector<double> >  metsigmatrixdxy  ( new std::vector<double>()  );  
  std::auto_ptr<std::vector<double> >  metsigmatrixdyx  ( new std::vector<double>()  );  
  std::auto_ptr<std::vector<double> >  metsigmatrixdyy  ( new std::vector<double>()  );  

  //-----------------------------------------------------------------
  edm::Handle<std::vector<pat::MET> > mets;
  iEvent.getByLabel(inputTag, mets);

  if(mets.isValid()) {
    edm::LogInfo("RootTupleMakerV2_METInfo") << "Total # METs: " << mets->size();

    for( std::vector<pat::MET>::const_iterator it = mets->begin(); it != mets->end(); ++it ) {

      pat::MET::METUncertainty shift = pat::MET::NoShift;
      pat::MET::METUncertaintyLevel lev = pat::MET::Type1;

      if(uncertainty=="NoShift")               shift=pat::MET::NoShift;
      else if(uncertainty=="JetEnUp")          shift=pat::MET::JetEnUp;
      else if(uncertainty=="JetEnDown")        shift=pat::MET::JetEnDown;
      else if(uncertainty=="JetResUp")         shift=pat::MET::JetResUp;
      else if(uncertainty=="JetResDown")       shift=pat::MET::JetResDown;
      else if(uncertainty=="MuonEnUp")         shift=pat::MET::MuonEnUp;
      else if(uncertainty=="MuonEnDown")       shift=pat::MET::MuonEnDown;
      else if(uncertainty=="ElectronEnUp")     shift=pat::MET::ElectronEnUp;
      else if(uncertainty=="ElectronEnDown")   shift=pat::MET::ElectronEnDown;
      else if(uncertainty=="TauEnUp")          shift=pat::MET::TauEnUp;
      else if(uncertainty=="TauEnDown")        shift=pat::MET::TauEnDown;
      else if(uncertainty=="UnclusteredEnUp")  shift=pat::MET::UnclusteredEnUp;
      else if(uncertainty=="UnclusteredEnDown")shift=pat::MET::UnclusteredEnDown;
      else edm::LogError("RootTupleMakerV2_METError") << "Error! Can't find MET uncertainty label: " << uncertainty;

      if(level=="Raw")          lev = pat::MET::Raw;
      else if(level=="Type1")   lev = pat::MET::Type1;
      else if(level=="Type1p2") lev = pat::MET::Type1p2;
      //else if(level=="Calo")    lev = pat::MET::Calo;//not until 7_4_X
      else edm::LogError("RootTupleMakerV2_METError") << "Error! Can't find MET uncertainty level label: " << level;

      met->push_back( it->shiftedPt(shift,lev) );
      metphi->push_back( it->shiftedPhi(shift,lev) );
      sumet->push_back( it->shiftedSumEt(shift,lev) );
      
      if ( store_uncorrected_MET ) {
	metuncorr->push_back( it->uncorrectedPt(pat::MET::uncorrALL) );
	metphiuncorr->push_back( it->uncorrectedPhi(pat::MET::uncorrALL) );
	sumetuncorr->push_back( it->sumEt() - it->corSumEt(pat::MET::uncorrALL) );
      }

      if ( store_MET_significance ) {
        //-- See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMETSignificance#Known_Issues
	double sigmaX2= it->getSignificanceMatrix()(0,0);
	double sigmaY2= it->getSignificanceMatrix()(1,1);
	double significance = -1;
	if(sigmaX2<1.e10 && sigmaY2<1.e10) significance = it->significance();
	//--

	metsig->push_back( significance );
	metsigmatrixdxx->push_back( it->getSignificanceMatrix()(0,0) );
	metsigmatrixdxy->push_back( it->getSignificanceMatrix()(0,1) );
	metsigmatrixdyx->push_back( it->getSignificanceMatrix()(1,0) );
	metsigmatrixdyy->push_back( it->getSignificanceMatrix()(1,1) );
	//See DataFormats/METReco/src/MET.cc
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
  if ( store_MET_significance ) {
    iEvent.put( metsig, prefix + "METSig" + suffix );
    iEvent.put( metsigmatrixdxx, prefix + "METSigMatrixDXX" + suffix );
    iEvent.put( metsigmatrixdxy, prefix + "METSigMatrixDXY" + suffix );
    iEvent.put( metsigmatrixdyx, prefix + "METSigMatrixDYX" + suffix );
    iEvent.put( metsigmatrixdyy, prefix + "METSigMatrixDYY" + suffix );
  }
}
