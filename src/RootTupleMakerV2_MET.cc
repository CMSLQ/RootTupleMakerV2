#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_MET.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

RootTupleMakerV2_MET::RootTupleMakerV2_MET(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    metToken(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("InputTag"))),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    store_uncorrected_MET (iConfig.getParameter<bool>  ("StoreUncorrectedMET")),
    store_MET_significance (iConfig.getParameter<bool>  ("StoreMETSignificance")),
    uncertainty (iConfig.getParameter<std::string>  ("Uncertainty")),
    corLevel       (iConfig.getParameter<std::string>  ("CorrectionLevel"))
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
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken, mets);

  if(mets.isValid()) {
    edm::LogInfo("RootTupleMakerV2_METInfo") << "Total # METs: " << mets->size();

    for( pat::METCollection::const_iterator it = mets->begin(); it != mets->end(); ++it ) {

      pat::MET::METUncertainty shift = pat::MET::NoShift;
      pat::MET::METCorrectionLevel level = pat::MET::Type1;

      if(     uncertainty=="NoShift")          shift=pat::MET::NoShift;
      else if(uncertainty=="JetResUp")         shift=pat::MET::JetResUp;
      else if(uncertainty=="JetResDown")       shift=pat::MET::JetResDown;
      else if(uncertainty=="JetEnUp")          shift=pat::MET::JetEnUp;
      else if(uncertainty=="JetEnDown")        shift=pat::MET::JetEnDown;
      else if(uncertainty=="MuonEnUp")         shift=pat::MET::MuonEnUp;
      else if(uncertainty=="MuonEnDown")       shift=pat::MET::MuonEnDown;
      else if(uncertainty=="ElectronEnUp")     shift=pat::MET::ElectronEnUp;
      else if(uncertainty=="ElectronEnDown")   shift=pat::MET::ElectronEnDown;
      else if(uncertainty=="TauEnUp")          shift=pat::MET::TauEnUp;
      else if(uncertainty=="TauEnDown")        shift=pat::MET::TauEnDown;
      else if(uncertainty=="UnclusteredEnUp")  shift=pat::MET::UnclusteredEnUp;
      else if(uncertainty=="UnclusteredEnDown")shift=pat::MET::UnclusteredEnDown;
      else if(uncertainty=="PhotonEnUp")       shift=pat::MET::PhotonEnUp;
      else if(uncertainty=="PhotonEnDown")     shift=pat::MET::PhotonEnDown;
      else if(uncertainty=="JetResUpSmear")    shift=pat::MET::JetResUpSmear;
      else if(uncertainty=="JetResDownSmear")  shift=pat::MET::JetResDownSmear;
      else edm::LogError("RootTupleMakerV2_METError") << "Error! Can't find MET uncertainty label: " << uncertainty;

      bool applyCorrection = true;
      if(corLevel=="Raw")                level = pat::MET::Raw;
      else if(corLevel=="Type1")         level = pat::MET::Type1;
      else if(corLevel=="Type01")        level = pat::MET::Type1;
      else if(corLevel=="TypeXY")        level = pat::MET::Type1;
      else if(corLevel=="Type1XY")       level = pat::MET::Type1;
      else if(corLevel=="Type01XY")      level = pat::MET::Type1;
      else if(corLevel=="Type1Smear")    level = pat::MET::Type1;
      else if(corLevel=="Type01Smear")   level = pat::MET::Type1;
      else if(corLevel=="Type1SmearXY")  level = pat::MET::Type1;
      else if(corLevel=="Type01SmearXY") level = pat::MET::Type1;
      else if(corLevel=="RawCalo")       level = pat::MET::RawCalo;
      // this is our adhoc way of not running the shiftedPt, shiftedPhi, shiftedSumEt on the recorrected MET collections
      //  (which we don't want to do because they aren't filled, so accessing them causes a seg fault)
      else if(corLevel=="NoCorrection")  applyCorrection = false;
      else edm::LogError("RootTupleMakerV2_METError") << "Error! Can't find MET correction level label: " << corLevel;

      if(applyCorrection)
      {
        met->push_back( it->shiftedPt(shift,level) );
        metphi->push_back( it->shiftedPhi(shift,level) );
        sumet->push_back( it->shiftedSumEt(shift,level) );
      }
      else
      {
        met->push_back( it->pt() );
        metphi->push_back( it->phi() );
        sumet->push_back( it->sumEt() );
      }
      
      if ( store_uncorrected_MET ) {
        // this will not work running in CMSSW_7_4_12+ on miniAOD v1
        // see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#ETmiss
        metuncorr->push_back( it->uncorPt() );
        metphiuncorr->push_back( it->uncorPhi() );
        sumetuncorr->push_back( it->uncorSumEt() );
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
